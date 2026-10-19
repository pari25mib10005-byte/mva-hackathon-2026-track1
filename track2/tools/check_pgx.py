#!/usr/bin/env python3
"""
Local pharmacogenomic (PGx) safety check for drugs this patient may receive.

RUN THIS ONLY ON YOUR OWN MACHINE. It reads the VCF locally and needs no
network access. Share only the gene-level SUMMARY lines, never genotypes or
positions. Results are for clinician review (CPIC/DPWG guidance), not advice.

Why these genes: MVA carries a childhood-cancer risk, so the most likely future
drugs are chemotherapy and supportive-care agents. CEP72 and CYP3A5 relate to
vincristine nerve damage (relevant to the patient's muscle feature); TPMT/NUDT15,
DPYD and UGT1A1 to thiopurine, fluoropyrimidine and irinotecan toxicity; G6PD to
rasburicase and oxidant drugs; CYP2C19/CYP2C9/SLCO1B1 to common supportive drugs.

Usage:
    python check_pgx.py --vcf ~/mva/data/patient.vcf.gz
Notes
  - Positions and alleles are GRCh38, forward strand, verified against Ensembl
    (October 2026). For a GRCh37 VCF, rerun after lifting positions over.
  - A site missing from a single-sample VCF is presumed homozygous reference;
    a true no-call cannot be distinguished without a gVCF or the reads.
"""
import argparse, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_vcf_genes import opener, detect_build

# gene, rsID, chrom, GRCh38 pos, ref, variant allele, star/name, function
SITES = [
    ("CEP72", "rs924607", "5", 609978, "C", "T", "risk allele", "vincristine neuropathy risk"),
    ("CYP3A5", "rs776746", "7", 99672916, "T", "C", "*3", "no function"),
    ("TPMT", "rs1800462", "6", 18143724, "C", "G", "*2", "no function"),
    ("TPMT", "rs1800460", "6", 18138997, "C", "T", "*3B", "no function"),
    ("TPMT", "rs1142345", "6", 18130687, "T", "C", "*3C", "no function"),
    ("NUDT15", "rs116855232", "13", 48045719, "C", "T", "*3", "no function"),
    ("DPYD", "rs3918290", "1", 97450058, "C", "T", "*2A", "no function"),
    ("DPYD", "rs55886062", "1", 97515787, "A", "C", "*13", "no function"),
    ("DPYD", "rs67376798", "1", 97082391, "T", "A", "c.2846A>T", "decreased function"),
    ("DPYD", "rs75017182", "1", 97579893, "G", "C", "HapB3", "decreased function"),
    ("UGT1A1", "rs887829", "2", 233759924, "C", "T", "*80 (tags *28)", "decreased function"),
    ("G6PD", "rs1050828", "X", 154536002, "C", "T", "A- (202A)", "deficiency"),
    ("G6PD", "rs5030868", "X", 154534419, "G", "A", "Mediterranean", "deficiency"),
    ("CYP2C19", "rs4244285", "10", 94781859, "G", "A", "*2", "no function"),
    ("CYP2C19", "rs4986893", "10", 94780653, "G", "A", "*3", "no function"),
    ("CYP2C19", "rs12248560", "10", 94761900, "C", "T", "*17", "increased function"),
    ("CYP2C9", "rs1799853", "10", 94942290, "C", "T", "*2", "decreased function"),
    ("CYP2C9", "rs1057910", "10", 94981296, "A", "C", "*3", "no function"),
    ("SLCO1B1", "rs4149056", "12", 21178615, "T", "C", "c.521T>C", "decreased function"),
]

# Family-history sites (not drug response): inherited thrombophilia, a known contributor to
# recurrent pregnancy loss. A variant in the child means one parent carries it.
FAMILY = [
    ("F5", "rs6025", "1", 169549811, "C", "T", "Factor V Leiden"),
    ("F2", "rs1799963", "11", 46739505, "G", "A", "Prothrombin G20210A"),
]

DRUGS = {
    "CEP72": "vincristine (nerve damage)", "CYP3A5": "vincristine exposure", "TPMT": "thiopurines",
    "NUDT15": "thiopurines", "DPYD": "5-FU / capecitabine", "UGT1A1": "irinotecan",
    "G6PD": "rasburicase, oxidant drugs", "CYP2C19": "PPIs, voriconazole, clopidogrel",
    "CYP2C9": "NSAIDs, phenytoin", "SLCO1B1": "statins (myopathy)",
}


def verdict(gene, calls):
    """calls: list of (site, copies or None). Returns a one-line, genotype-free verdict."""
    copies = {s[1]: c for s, c in calls}
    total = sum(c or 0 for c in copies.values())
    if gene == "CEP72":
        return "ACTION: high-risk genotype (2 risk alleles)" if total >= 2 else "no high-risk genotype"
    if gene == "CYP3A5":
        return "non-expresser (*3/*3; linked to more vincristine neuropathy)" if total >= 2 else "expresser (at least one functional copy)"
    if gene == "CYP2C19":
        lof = (copies.get("rs4244285") or 0) + (copies.get("rs4986893") or 0)
        inc = copies.get("rs12248560") or 0
        if lof >= 2:
            return "ACTION: likely poor metaboliser"
        if lof == 1:
            return "intermediate metaboliser"
        return "rapid/ultrarapid metaboliser (*17)" if inc else "no tested reduced-function alleles"
    if gene == "UGT1A1":
        return "ACTION: likely poor metaboliser (*28-type, 2 copies)" if total >= 2 else (
            "intermediate (1 copy)" if total == 1 else "no tested reduced-function alleles")
    if gene == "G6PD":
        return "ACTION: deficiency allele present" if total >= 1 else "no tested deficiency alleles"
    if total >= 2:
        return "ACTION: two reduced/no-function alleles"
    if total == 1:
        return "ACTION: one reduced/no-function allele (intermediate)"
    return "no tested reduced-function alleles"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample", help="sample column (default: first)")
    ap.add_argument("--details", action="store_true", help="also print per-site calls (KEEP LOCAL)")
    a = ap.parse_args()

    build = detect_build(a.vcf)
    if build and build != "GRCh38":
        sys.exit(f"VCF is {build}; these positions are GRCh38. Lift over first.")
    want = {(s[2], s[3]): s for s in SITES + [f[:6] + ("", "") for f in FAMILY]}
    found, col = {}, 9
    with opener(a.vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                if a.sample:
                    col = hdr.index(a.sample)
                continue
            fld = line.split("\t", 10)
            key = (fld[0].replace("chr", ""), int(fld[1]))
            if key not in want:
                continue
            site = want[key]
            rest = line.rstrip("\n").split("\t")
            alts = rest[4].split(",")
            fmt = rest[8].split(":") if len(rest) > 8 else []
            gt = dict(zip(fmt, rest[col].split(":"))).get("GT", "./.") if len(rest) > col else "./."
            idx = [int(x) for x in re.split(r"[/|]", gt) if x.isdigit()]
            if rest[3].upper() != site[4]:
                found[key] = None  # reference mismatch: do not guess
                continue
            target = alts.index(site[5]) + 1 if site[5] in alts else None
            found[key] = sum(1 for i in idx if target is not None and i == target)

    print(f"Pharmacogenomic summary (safe to share; {len(SITES)} sites in {len(DRUGS)} genes):")
    for gene in DRUGS:
        calls = [(s, found.get((s[2], s[3]), 0)) for s in SITES if s[0] == gene]
        print(f"  {gene} [{DRUGS[gene]}]: {verdict(gene, calls)}")
        if a.details:
            for s, c in calls:
                state = "not in VCF (presumed reference)" if (s[2], s[3]) not in found else f"{c} variant copies"
                print(f"      [LOCAL ONLY] {s[1]} {s[6]}: {state}")
    print("\nFamily-history check (safe to share):")
    for gene, rs, c, pos, ref, alt, name in FAMILY:
        n = found.get((c, pos), 0) or 0
        print(f"  {name} ({gene}): " + ("PRESENT in the child: one parent carries it (relevant to recurrent pregnancy loss)"
                                       if n else "absent in the child (a parent may still carry it)"))
    print("\nFor clinician review only. Share the summary lines, never genotypes or positions.")


if __name__ == "__main__":
    main()
