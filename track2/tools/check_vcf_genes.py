#!/usr/bin/env python3
"""
Local check: does the patient's VCF carry rare variants in a list of genes?
Used to test whether features not explained by MVA could have a SECOND genetic
cause (default list: antenatal Bartter / salt-wasting tubulopathy genes).

RUN THIS ONLY ON YOUR OWN MACHINE. It reads the VCF locally; only gene
coordinates are fetched from Ensembl. Share only the gene-level SUMMARY lines
(e.g. "SLC12A1: none"), never positions or genotypes, per the hackathon rules.

Usage:
    python check_vcf_genes.py --vcf path/to/patient.vcf.gz
    python check_vcf_genes.py --vcf patient.vcf --genes SLC12A1,KCNJ1 --build GRCh37 --sample CHILD_ID

Notes
  - Build is auto-detected from the VCF header (contig lengths / reference line);
    override with --build GRCh37 or GRCh38.
  - Same quality filter as Track 1 (QUAL>=30, DP>=10, GQ>=20).
  - Population frequency and consequence come from gnomAD INFO keys if present,
    otherwise from Ensembl VEP REST (only chromosome, position and alleles of
    variants inside the listed genes are sent; no sample IDs or genotypes).
    Plain INFO/AF is ignored: in GATK VCFs it is the sample's own allele fraction.
  - Only rare (max population AF <= 1%) protein- or splice-altering variants count.
  - Recessive genes need two hits (homozygous, or two different heterozygous
    variants); MAGED2 is X-linked (one hit is enough in a male).
"""
import argparse, gzip, json, re, sys, time, urllib.request

DEFAULT_GENES = ["SLC12A1", "KCNJ1", "CLCNKB", "BSND", "MAGED2"]
X_LINKED = {"MAGED2"}
AF_KEYS = ["gnomAD_AF", "gnomADe_AF", "gnomADg_AF", "gnomad_af", "AF_popmax", "gnomAD_AF_popmax", "ExAC_AF", "MAX_AF"]
# NOTE: plain INFO/AF is NOT used: in raw GATK VCFs it is the sample's own allele fraction (0.5/1.0), not population frequency.
DAMAGING = {"transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
            "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant",
            "protein_altering_variant", "splice_region_variant", "splice_donor_5th_base_variant"}
CHR1_LEN = {"248956422": "GRCh38", "249250621": "GRCh37"}
PAD = 50  # include splice-region flanks


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else \
        open(path, encoding="utf-8", errors="replace")


def detect_build(path):
    with opener(path) as f:
        for line in f:
            if not line.startswith("##"):
                break
            m = re.search(r"ID=(?:chr)?1,.*length=(\d+)", line)
            if m and m.group(1) in CHR1_LEN:
                return CHR1_LEN[m.group(1)]
            low = line.lower()
            if low.startswith("##reference") or low.startswith("##assembly"):
                if "38" in low:
                    return "GRCh38"
                if "37" in low or "hg19" in low or "b37" in low:
                    return "GRCh37"
    return None


def gene_region(symbol, build):
    host = "https://rest.ensembl.org" if build == "GRCh38" else "https://grch37.rest.ensembl.org"
    req = urllib.request.Request(f"{host}/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json",
                                 headers={"Content-Type": "application/json"})
    d = json.load(urllib.request.urlopen(req, timeout=60))
    return str(d["seq_region_name"]), d["start"] - PAD, d["end"] + PAD


def vep_annotate(variants, build):
    """Population frequency (max gnomAD/1000G) and most severe consequence via Ensembl VEP REST."""
    host = "https://rest.ensembl.org" if build == "GRCh38" else "https://grch37.rest.ensembl.org"
    out = {}
    for i in range(0, len(variants), 150):
        chunk = variants[i:i + 150]
        body = json.dumps({"variants": [f"{c} {p} . {r} {a} . . ." for c, p, r, a in chunk]}).encode()
        for attempt in range(4):
            try:
                req = urllib.request.Request(host + "/vep/homo_sapiens/region", data=body,
                                             headers={"Content-Type": "application/json", "Accept": "application/json"})
                res = json.load(urllib.request.urlopen(req, timeout=180))
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(3 * (attempt + 1))
        for r in res:
            f = r["input"].split()
            key = (f[0], int(f[1]), f[3], f[4])
            afs = []
            for cv in r.get("colocated_variants", []) or []:
                for allele_map in (cv.get("frequencies") or {}).values():
                    afs += [v for k, v in allele_map.items() if isinstance(v, (int, float)) and k in
                            ("af", "gnomade", "gnomadg", "gnomad", "gnomade_popmax", "gnomadg_popmax")]
            out[key] = {"af": max(afs) if afs else 0.0, "csq": r.get("most_severe_consequence", ""),
                        "known": bool(r.get("colocated_variants"))}
        time.sleep(0.3)
    return out


def parse_info(info):
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def allele_freq(info):
    for k in AF_KEYS:
        if k in info:
            vals = [float(x) for x in re.split(r"[,|]", info[k]) if re.fullmatch(r"[0-9.eE+-]+", x or "")]
            if vals:
                return max(vals)
    return None


def consequence(info):
    for k in ("ANN", "CSQ"):
        if k in info:
            terms = set()
            for ann in info[k].split(","):
                parts = ann.split("|")
                if len(parts) > 1:
                    terms.update(t for t in parts[1].split("&") if t)
            return ",".join(sorted(terms)[:4])
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--genes", default=",".join(DEFAULT_GENES))
    ap.add_argument("--build", choices=["GRCh37", "GRCh38"])
    ap.add_argument("--sample", help="sample column to use (default: first sample)")
    ap.add_argument("--max-af", type=float, default=0.01, help="rare-variant cut-off (default 1%%)")
    ap.add_argument("--details", action="store_true", help="also print variant lines (KEEP LOCAL, do not share)")
    ap.add_argument("--no-annotate", action="store_true", help="skip the Ensembl VEP lookup (offline)")
    ap.add_argument("--no-qc", action="store_true", help="skip the QUAL>=30, DP>=10, GQ>=20 filter")
    a = ap.parse_args()

    build = a.build or detect_build(a.vcf)
    if not build:
        sys.exit("Could not detect genome build from the VCF header; rerun with --build GRCh37 or --build GRCh38")
    genes = [g.strip().upper() for g in a.genes.split(",") if g.strip()]
    regions = {g: gene_region(g, build) for g in genes}
    print(f"Build: {build}; genes: {', '.join(genes)}")

    hits = {g: [] for g in genes}
    with opener(a.vcf) as f:
        col = 9
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if a.sample:
                    if a.sample not in header:
                        sys.exit(f"Sample {a.sample} not in VCF; samples: {', '.join(header[9:])}")
                    col = header.index(a.sample)
                continue
            fld = line.rstrip("\n").split("\t")
            chrom, pos = fld[0].replace("chr", ""), int(fld[1])
            for g, (c, s, e) in regions.items():
                if chrom == c and s <= pos <= e:
                    fmt = fld[8].split(":") if len(fld) > 8 else []
                    smp = fld[col].split(":") if len(fld) > col else []
                    gt = dict(zip(fmt, smp)).get("GT", "./.")
                    alleles = [x for x in re.split(r"[/|]", gt) if x not in (".", "")]
                    if not alleles or all(x == "0" for x in alleles):
                        continue  # no alternate allele in this sample
                    fmtd = dict(zip(fmt, smp))
                    try:
                        qual = float(fld[5]) if fld[5] != "." else 0.0
                        dp = int(fmtd.get("DP", "0") or 0)
                        gq = int(fmtd.get("GQ", "0") or 0)
                    except ValueError:
                        qual, dp, gq = 0.0, 0, 0
                    if not a.no_qc and (qual < 30 or dp < 10 or gq < 20):
                        continue  # same quality filter as Track 1
                    info = parse_info(fld[7]) if len(fld) > 7 else {}
                    af = allele_freq(info)
                    zyg = "hom" if len(set(alleles)) == 1 and len(alleles) == 2 else ("hemi" if len(alleles) == 1 else "het")
                    for alt in fld[4].split(","):
                        if alt in ("*", "<NON_REF>"):
                            continue
                        hits[g].append({"chrom": chrom, "pos": pos, "ref": fld[3], "alt": alt, "zyg": zyg, "af": af,
                                        "filter": fld[6], "csq": consequence(info)})

    todo = [(v["chrom"], v["pos"], v["ref"], v["alt"]) for g in genes for v in hits[g] if v["af"] is None]
    if todo and not a.no_annotate:
        print(f"Looking up population frequency and effect for {len(todo)} variants (Ensembl VEP) …")
        ann = vep_annotate(todo, build)
        for g in genes:
            for v in hits[g]:
                hit = ann.get((v["chrom"], v["pos"], v["ref"], v["alt"]))
                if hit and v["af"] is None:
                    v["af"], v["csq"] = hit["af"], v["csq"] or hit["csq"]

    print("\nGene-level summary (safe to share):")
    for g in genes:
        rare = [v for v in hits[g] if v["af"] is not None and v["af"] <= a.max_af]
        unknown = [v for v in hits[g] if v["af"] is None]
        damaging = [v for v in rare + unknown if any(c in DAMAGING for c in (v["csq"] or "").split(","))]
        noncoding = len(rare + unknown) - len(damaging)
        cand = damaging
        if not cand:
            verdict = "no rare protein- or splice-altering variants" + (f" ({noncoding} rare non-coding, low prior)" if noncoding else "")
        elif g in X_LINKED:
            verdict = f"{len(cand)} rare protein/splice-altering variant(s): REVIEW (X-linked, one hit can suffice)"
        elif any(v["zyg"] == "hom" for v in cand) or len(cand) >= 2:
            verdict = f"{len(cand)} rare protein/splice-altering variant(s), two-hit pattern: REVIEW (possible recessive cause)"
        else:
            verdict = "1 rare protein/splice-altering heterozygous variant: carrier pattern, unlikely to cause disease alone"
        print(f"  {g}: {verdict}")
        if a.details:
            for v in hits[g]:
                print(f"      [LOCAL ONLY] {v['pos']} {v['ref']}>{v['alt']} {v['zyg']} AF={v['af']} "
                      f"FILTER={v['filter']} {v['csq']}")
    print("\nReminder: share only the summary lines above, never positions or genotypes.")


if __name__ == "__main__":
    main()
