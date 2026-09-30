"""Figure 1: 5 Mb binned depth ratio and mean het BAF on chr1/20/21/22.
Input: mosaic_bins.csv (aggregate, no genotype-level data)."""
import pandas as pd, matplotlib.pyplot as plt
d = pd.read_csv("mosaic_bins.csv")
d["artefact"] = (d.dp_ratio > 1.15) | (d.mean_baf < 0.45)
chroms = [1, 20, 21, 22]
w = {1: 4, 20: 1.6, 21: 1.2, 22: 1.3}
fig, ax = plt.subplots(2, 4, figsize=(12, 5.2), sharey="row",
                       gridspec_kw={"width_ratios": [w[c] for c in chroms]})
for j, c in enumerate(chroms):
    s = d[d.chrom == c]; ok, bad = s[~s.artefact], s[s.artefact]
    for i, (col, ref) in enumerate([("dp_ratio", 1.0), ("mean_baf", 0.5)]):
        a = ax[i, j]
        a.axhline(ref, color="grey", lw=0.8, ls="--")
        a.plot(s.bin_Mb + 2.5, s[col], color="#bbb", lw=0.8, zorder=1)
        a.scatter(ok.bin_Mb + 2.5, ok[col], s=18, color="#2b6cb0", zorder=2, label="euchromatic bin")
        a.scatter(bad.bin_Mb + 2.5, bad[col], s=26, color="#d9480f", marker="D", zorder=3,
                  label="repeat / centromeric bin")
        a.spines[["top", "right"]].set_visible(False)
        if i == 1: a.set_xlabel(f"chr{c} position (Mb)")
    ax[0, j].set_title(f"chr{c}")
ax[0, 0].set_ylabel("Depth / genome median"); ax[1, 0].set_ylabel("Mean het BAF")
ax[1, 0].set_ylim(0.28, 0.53)
ax[0, 0].legend(frameon=False, fontsize=8, loc="upper left")
fig.suptitle("Apparent chr20/21/22 gains are confined to repeat regions; euchromatin matches chr1",
             fontsize=11)
fig.tight_layout()
fig.savefig("figures/fig1_mosaic_bins.png", dpi=200)
