import subprocess, sys, numpy as np, pandas as pd

VCF = sys.argv[1] if len(sys.argv) > 1 else "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
AUTO = [str(i) for i in range(1, 23)]

# High-quality biallelic heterozygous SNPs only
cmd = ["bcftools","query",
       "-i", 'TYPE="snp" && GT="het" && FMT/DP>=20 && FMT/GQ>=30 && QUAL>=50',
       "-f", "%CHROM\t%POS\t[%AD]\t[%DP]\n", VCF]
rows = []
for line in subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True).stdout:
    c, p, ad, dp = line.rstrip().split("\t")
    if c not in AUTO: continue
    a = ad.split(",")
    if len(a) != 2: continue
    ref, alt = int(a[0]), int(a[1])
    rows.append((c, int(p), alt/(ref+alt), int(dp)))
df = pd.DataFrame(rows, columns=["chrom","pos","baf","dp"])
df["bin"] = df.pos // 5_000_000
bad = df.groupby(["chrom","bin"]).dp.median()
bad = bad[bad > 1.15 * df.dp.median()].index
df = df[~df.set_index(["chrom","bin"]).index.isin(bad)]
print(f"excluded {len(bad)} high-depth 5 Mb bins (centromeric/repeat artefacts)")
print(f"{len(df):,} het SNPs used\n")

def summarise(g):
    dev = np.abs(g.baf - 0.5)
    return pd.Series({
        "n": len(g),
        "median_dp": g.dp.median(),
        "mean_|BAF-0.5|": dev.mean(),
        "baf_sd": g.baf.std(),
        # fraction of hets pushed outside the diploid band — sensitive to a mosaic shift
        "frac_outside_0.4-0.6": ((g.baf < 0.40) | (g.baf > 0.60)).mean(),
    })
res = df.groupby("chrom").apply(summarise)
res = res.loc[AUTO]

# Normalise depth to genome median; z-score BAF deviation against other autosomes
res["dp_ratio"] = res.median_dp / res.median_dp.median()
med, mad = res["mean_|BAF-0.5|"].median(), (res["mean_|BAF-0.5|"] - res["mean_|BAF-0.5|"].median()).abs().median()
res["baf_dev_z"] = (res["mean_|BAF-0.5|"] - med) / (1.4826*mad if mad > 0 else 1)
pd.set_option("display.width", 160)
print(res.round(4).to_string())
print("\nFlagged (|z| > 3 or dp_ratio outside 0.97-1.03):")
print(res[(res.baf_dev_z.abs() > 3) | (res.dp_ratio < 0.97) | (res.dp_ratio > 1.03)].round(4).to_string())
res.round(5).to_csv("mosaic_per_chrom.csv")
