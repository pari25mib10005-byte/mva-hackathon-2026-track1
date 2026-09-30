"""Figure 3: variant filtering funnel (counts from report §3–4)."""
import math, matplotlib.pyplot as plt
steps = [
    ("~5,000,000", "All calls in the WGS VCF"),
    ("191", "In BUB1B, CEP57, TRIP13 (±10 kb)"),
    ("184", "Pass QC (QUAL≥30, DP≥10, GQ≥20), multi-allelics split"),
    ("2", "Protein-altering (VEP HIGH/MODERATE), rare"),
    ("1 gene", "Two rare coding hets in one AR gene: BUB1B"),
]
vals = [5_000_000, 191, 184, 2, 1]
w = [0.25 + 0.75 * math.log10(v) / math.log10(vals[0]) for v in vals]
cols = ["#cfe0f3", "#a9c9ec", "#7eaee0", "#e8835a", "#c92a2a"]
fig, ax = plt.subplots(figsize=(9, 4.2))
for i, (n, lab) in enumerate(steps):
    y = len(steps) - 1 - i
    ax.barh(y, w[i], left=-w[i] / 2, height=0.72, color=cols[i], edgecolor="white")
    ax.text(0, y, n, ha="center", va="center", fontsize=11, fontweight="bold",
            color="white" if i >= 3 else "#1a1a1a")
    ax.text(0.56, y, lab, ha="left", va="center", fontsize=9.5)
ax.set_xlim(-0.55, 1.75); ax.set_ylim(-0.6, len(steps) - 0.4); ax.axis("off")
ax.set_title("Variant filtering: from whole genome to the BUB1B compound-heterozygous pair", fontsize=11)
ax.text(-0.55, -0.75, "Bar width on log scale. Scripts: run.sh, annotate.py", fontsize=8, color="#666")
fig.tight_layout(); fig.savefig("figures/fig3_filter_funnel.png", dpi=200)
