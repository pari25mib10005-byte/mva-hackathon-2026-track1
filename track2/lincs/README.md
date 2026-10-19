# Track 2, evidence line 3: LINCS signature reversal

Asks, with **no pathway assumption**: *which compounds push human cells in the opposite direction to losing the disease gene?* It also tests our NAD⁺ hypothesis with a pre-specified statistic.

```bash
python3 signature_reversal.py --genes BUB1B,CEP57,TRIP13 --config ../pipeline/config/bub1b.json --out results/
```
Windows: `py -X utf8` instead of `python3`. The first run downloads two public LINCS L1000 consensus libraries (~15 MB each, via Enrichr) and the WHO ATC L01 list (via ChEMBL) into `data/`.

## Method
- **Signatures:** LINCS L1000 CRISPR-knockout consensus signature of the disease gene (genes up / down) vs. consensus signatures of 5,425 compounds.
- **Score:** (reversing overlaps − mimicking overlaps) / signature size. **Direction test:** binomial sign test on reversing vs mimicking overlaps; Benjamini–Hochberg FDR.
- **Pre-specified hypothesis (fixed before results):** NAD⁺/sirtuin-axis compounds (nicotinamide, NAD⁺/Nadide, SRT-1720, resveratrol, and NR/NMN/SRT-2104 if present) reverse the signature more than other compounds; one-sided Mann–Whitney U.
- **Two halves of the signature:** the knockout "up" half contains the cell's **alarm response** to chromosome errors (p53 / apoptosis / interferon genes); the "down" half contains **cell-division genes** (incl. BUB1B, BUB1, TRIP13, CDC20). Reversing the alarm half would let aneuploid cells survive, so compounds that switch off ≥ 3 alarm genes are flagged and excluded.
- **Safety review:** same target/keyword rules as `../pipeline/repurpose.py`, plus wrong-direction patterns (sirtuin / NAMPT / JAK / p53 inhibitors, from `config/bub1b.json → lincs`) and **every WHO ATC L01 antineoplastic** (catches chemotherapies whose database mechanism is missing, e.g. procarbazine, estramustine).
- **Exploratory list:** compounds that restore the division-gene half (hypergeometric, FDR < 0.05) without switching off the alarm.

Outputs: `lincs_summary.md`, `reversal_scores.csv` (all compounds), `top_reversers_reviewed.csv`, `restore_candidates.csv`, `lincs_provenance.json`.

**Limits:** consensus signatures come from cancer cell lines; transcript changes are not protein stability (our lead acts on BUBR1 protein turnover, which L1000 cannot see).
