#!/usr/bin/env python3
"""
Local, genome-wide check: could a SECOND genetic cause explain features that the
disease gene does not? (Answers "is it genetic or not?" for those features.)

For each HPO feature it takes EVERY gene known to cause that feature (HPO),
reads each gene's inheritance pattern (HPO phenotype.hpoa), scans the coding
exons (+/-10 bp, MANE Select) of those genes in the patient's VCF, looks up
population frequency, consequence and AlphaMissense (Ensembl VEP), and applies
inheritance-aware rules:
  recessive  : rare (AF <= 1%) damaging variant, homozygous or >= 2 variants
  dominant   : ultra-rare (AF <= 0.01%) damaging heterozygous variant
  X-linked   : rare damaging variant (male patient: one copy suffices)
  damaging   = on the CANONICAL transcript: loss-of-function, AlphaMissense
               "likely pathogenic" missense, or an in-frame indel absent from gnomAD;
               variants ClinVar calls benign/likely benign are excluded

RUN ONLY ON YOUR OWN MACHINE. Share only the summary lines (gene names and
categories), never positions or genotypes. Only chromosome, position and
alleles of variants inside the listed genes are sent to Ensembl.

Usage:
    python check_feature_genes.py --vcf patient.vcf.gz --hpo HP:XXXXXXX,HP:YYYYYYY \
        --skip-genes BUB1B,CEP57,TRIP13
First run downloads ~110 MB of public reference files into tools/data/.
"""
import argparse, bisect, collections, gzip, json, os, re, sys, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from check_vcf_genes import opener, detect_build

DATA = os.path.join(HERE, "data")
P2G = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt"
HPOA = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa"
MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
INHERIT = {"HP:0000007": "AR", "HP:0000006": "AD", "HP:0001417": "XL", "HP:0001419": "XL", "HP:0001423": "XL"}
LOF = {"transcript_ablation", "stop_gained", "frameshift_variant", "splice_acceptor_variant",
       "splice_donor_variant", "start_lost"}
INFRAME = {"inframe_insertion", "inframe_deletion"}


def fetch(url, name):
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, name)
    if not os.path.exists(path):
        print(f"      downloading {name} …")
        urllib.request.urlretrieve(url, path)
    return path


def mane_path():
    path = os.path.join(DATA, "mane.gtf.gz")
    if os.path.exists(path):
        return path
    listing = urllib.request.urlopen(MANE_DIR, timeout=60).read().decode()
    name = re.findall(r"MANE\.GRCh38\.v[0-9.]+\.ensembl_genomic\.gtf\.gz", listing)[0]
    return fetch(MANE_DIR + name, "mane.gtf.gz")


def vep(variants):
    out = {}
    for i in range(0, len(variants), 150):
        chunk = variants[i:i + 150]
        body = json.dumps({"variants": [f"{c} {p} . {r} {a} . . ." for c, p, r, a in chunk],
                           "AlphaMissense": 1, "SpliceAI": 1, "canonical": 1}).encode()
        for attempt in range(5):
            try:
                req = urllib.request.Request("https://rest.ensembl.org/vep/homo_sapiens/region", data=body,
                                             headers={"Content-Type": "application/json", "Accept": "application/json"})
                res = json.load(urllib.request.urlopen(req, timeout=240))
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(5 * (attempt + 1))
        for r in res:
            f = r["input"].split()
            afs = []
            for cv in r.get("colocated_variants", []) or []:
                for amap in (cv.get("frequencies") or {}).values():
                    afs += [v for k, v in amap.items() if isinstance(v, (int, float)) and
                            k in ("af", "gnomade", "gnomadg", "gnomade_popmax", "gnomadg_popmax")]
            ams = [(t.get("alphamissense") or {}).get("am_pathogenicity") for t in r.get("transcript_consequences", []) or []]
            ams = [x for x in ams if isinstance(x, (int, float))]
            am = any((t.get("alphamissense") or {}).get("am_class") == "likely_pathogenic"
                     for t in r.get("transcript_consequences", []) or [])
            clin = sorted({c for cv in r.get("colocated_variants", []) or [] for c in (cv.get("clin_sig") or [])})
            # Judge the effect on the canonical (MANE-equivalent) transcript only: effects on minor
            # transcripts inflate false "damaging" calls.
            canon = [t for t in r.get("transcript_consequences", []) or [] if t.get("canonical")]
            terms = set()
            by_gene = collections.defaultdict(set)
            for t in canon:
                terms.update(t.get("consequence_terms", []))
                by_gene[t.get("gene_symbol", "")].update(t.get("consequence_terms", []))
            sai = []
            for t in r.get("transcript_consequences", []) or []:
                sp = t.get("spliceai") or {}
                sai += [v for k, v in sp.items() if k.startswith("DS_") and isinstance(v, (int, float))]
            out[(f[0], int(f[1]), f[3], f[4])] = {"af": max(afs) if afs else 0.0, "terms": terms, "am_lp": am,
                                                  "am": max(ams) if ams else None, "clin": clin,
                                                  "spliceai": max(sai) if sai else None, "by_gene": dict(by_gene),
                                                  "worst": r.get("most_severe_consequence", "")}
        print(f"      annotated {min(i + 150, len(variants))}/{len(variants)}")
        time.sleep(0.4)
    return out


def damaging(a):
    benign = {"benign", "likely_benign"}
    if a["clin"] and set(a["clin"]) <= benign | {"uncertain_significance"} and set(a["clin"]) & benign:
        return False  # ClinVar calls it benign/likely benign
    t = a["terms"]
    return bool(t & LOF) or ("missense_variant" in t and a["am_lp"]) or (bool(t & INFRAME) and a["af"] == 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--hpo", required=True, help="comma-separated HPO IDs of the unexplained features")
    ap.add_argument("--skip-genes", default="", help="genes already explaining the disease")
    ap.add_argument("--no-qc", action="store_true")
    ap.add_argument("--only-genes", default="", help="restrict the scan to these genes (fast re-check)")
    ap.add_argument("--details", action="store_true",
                    help="print variant-level detail for candidate genes WITHOUT positions (consequence, rarity, AlphaMissense, ClinVar)")
    ap.add_argument("--keep-hla", action="store_true", help="keep HLA genes (excluded by default: hypervariable, unreliable calls)")
    a = ap.parse_args()
    terms = [t.strip().upper().replace("HP", "HP:").replace("HP::", "HP:") for t in a.hpo.split(",") if t.strip()]
    skip = {g.strip().upper() for g in a.skip_genes.split(",") if g.strip()}
    only = {g.strip().upper() for g in a.only_genes.split(",") if g.strip()}
    if (detect_build(a.vcf) or "GRCh38") != "GRCh38":
        sys.exit("This scan uses GRCh38 annotations; the VCF looks like another build.")

    print("[1/5] Reference data (HPO genes, inheritance, MANE exons)")
    p2g = fetch(P2G, "phenotype_to_genes.txt")
    hpoa = fetch(HPOA, "phenotype.hpoa")
    gtf = mane_path()
    dis_inh = collections.defaultdict(set)
    for line in open(hpoa, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 10 and f[3] in INHERIT:
            dis_inh[f[0]].add(INHERIT[f[3]])
    term_genes, names, gene_inh = collections.defaultdict(set), {}, collections.defaultdict(set)
    for line in open(p2g, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if f[0] in terms and f[3] not in skip and (not only or f[3] in only) and (a.keep_hla or not f[3].startswith("HLA-")):
            term_genes[f[0]].add(f[3]); names[f[0]] = f[1]
            gene_inh[f[3]] |= dis_inh.get(f[4], set())
    genes = set().union(*term_genes.values()) if term_genes else set()
    print(f"      {len(genes)} genes across {len(term_genes)} features")

    exons = collections.defaultdict(list)
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if "\texon\t" not in line:
                continue
            c = line.split("\t", 8)
            m = re.search(r'gene_name "([^"]+)"', c[8])
            if m and m.group(1) in genes:
                exons[c[0].replace("chr", "")].append((int(c[3]) - 10, int(c[4]) + 10, m.group(1)))
    starts = {}
    for ch in exons:
        exons[ch].sort()
        starts[ch] = [e[0] for e in exons[ch]]
    covered = {e[2] for ch in exons for e in exons[ch]}
    print(f"      coding exons found for {len(covered)} genes")

    print("[2/5] Scanning the VCF (local)")
    hits, col = [], 9
    with opener(a.vcf) as f:
        for line in f:
            if line[0] == "#":
                continue
            fld = line.split("\t", 10)
            ch = fld[0].replace("chr", "")
            if ch not in exons:
                continue
            pos = int(fld[1])
            i = bisect.bisect_right(starts[ch], pos) - 1
            gs = set()
            j = i
            while j >= 0 and j > i - 40:
                s, e, g = exons[ch][j]
                if s <= pos <= e:
                    gs.add(g)
                j -= 1
            if not gs:
                continue
            r = line.rstrip("\n").split("\t")
            fmt = dict(zip(r[8].split(":"), r[col].split(":"))) if len(r) > col else {}
            al = [x for x in re.split(r"[/|]", fmt.get("GT", "./.")) if x not in (".", "")]
            if not al or all(x == "0" for x in al):
                continue
            try:
                ok = a.no_qc or (float(r[5]) >= 30 and int(fmt.get("DP", 0) or 0) >= 10 and int(fmt.get("GQ", 0) or 0) >= 20)
            except ValueError:
                ok = False
            if not ok:
                continue
            zyg = "hom" if len(set(al)) == 1 and len(al) == 2 else ("hemi" if len(al) == 1 else "het")
            for k, alt in enumerate(r[4].split(","), 1):
                if alt in ("*", "<NON_REF>") or str(k) not in al:
                    continue
                for g in gs:
                    hits.append((g, ch, pos, r[3], alt, zyg))
    print(f"      {len(hits)} coding/splice-region variant calls in these genes")

    print("[3/5] Population frequency, consequence and AlphaMissense (Ensembl VEP)")
    ann = vep(sorted({(h[1], h[2], h[3], h[4]) for h in hits}))

    print("[4/5] Applying inheritance-aware rules")
    per_gene = collections.defaultdict(list)
    detail = collections.defaultdict(list)
    for g, ch, pos, ref, alt, zyg in hits:
        x = ann.get((ch, pos, ref, alt))
        if x:
            x = {**x, "terms": x["by_gene"].get(g, set())}  # this gene's canonical transcript only
        if x and damaging(x):
            per_gene[g].append((x["af"], zyg))
            detail[g].append((x, zyg))
    cands = {}
    for g, vs in per_gene.items():
        inh = gene_inh.get(g) or {"AR", "AD"}
        rare = [v for v in vs if v[0] <= 0.01]
        ultra = [v for v in vs if v[0] <= 0.0001]
        if "AR" in inh and (any(z == "hom" for _, z in rare) or len(rare) >= 2):
            cands[g] = "recessive pattern (two hits)"
        elif "XL" in inh and rare:
            cands[g] = "X-linked pattern"
        elif "AD" in inh and ultra:
            cands[g] = "dominant pattern (ultra-rare damaging, inheritance unknown)"
        elif rare:
            cands.setdefault(g, None)
    print("[5/5] Summary")
    print("\nFeature-level summary (safe to share):")
    for t in terms:
        gs = term_genes.get(t, set())
        found = sorted(g for g in gs if cands.get(g))
        carriers = sorted(g for g in gs if g in cands and not cands[g])
        line = f"  {names.get(t, t)} ({t}): {len(gs)} genes scanned; "
        line += ("candidates: " + "; ".join(f"{g} [{cands[g]}]" for g in found)) if found else "no gene with a disease-consistent pattern"
        if carriers:
            line += f" (single rare damaging variants consistent with carrier status in {len(carriers)} recessive genes)"
        print(line)
    if a.details:
        print("\nVariant-level detail for candidate genes (no positions; share with care):")
        for g in sorted(c for c in cands if cands[c]):
            print(f"  {g} ({', '.join(sorted(gene_inh.get(g) or {'inheritance not annotated'}))}):")
            for x, zyg in detail[g]:
                am = f"AlphaMissense {x['am']:.2f}" if x["am"] is not None else "AlphaMissense n/a"
                am += f"; SpliceAI {x['spliceai']:.2f}" if x.get("spliceai") is not None else "; SpliceAI n/a"
                print(f"    - {', '.join(sorted(x['terms']))} (canonical transcript); {zyg}; max population AF {x['af']:.2g}; {am}; "
                      f"ClinVar: {', '.join(x['clin']) or 'no record'}")
        print("  Note: two heterozygous variants may sit on the same copy (cis); only parental testing can show they are on opposite copies.")
    print("\nCandidates need expert review (segregation, phenotype fit) before any conclusion.")
    print("Share only the lines above, never positions or genotypes.")


if __name__ == "__main__":
    main()
