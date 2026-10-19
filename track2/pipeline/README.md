# Track 2 — gene-agnostic repurposing pipeline

Takes a disease gene + a mechanism config, pulls drug–target evidence from **Open Targets** and **DGIdb**, removes anything predicted to worsen the disease, merges a curated literature layer, and ranks what is left. Standard-library Python 3.9+, no installs.

> Output is a list of **hypotheses for investigation, not treatment advice.**

## Run

```bash
python repurpose.py --config config/bub1b.json --out results/   # our case
python repurpose.py --config config/trip13.json --out results/  # other MVA genes
python repurpose.py --gene ANYGENE --out results/               # any gene, generic defaults
```

Outputs per gene: `*_candidates.csv` (ranked), `*_excluded.csv` (every rejected drug with its reason), `*_summary.md`, `*_provenance.json` (run time, sources, full config, gene set).

## How it decides

1. **Gene set** = index gene + curated modifier genes (with the direction the mechanism needs) + high-confidence STRING neighbours (score ≥ 0.9).
2. **Hard exclusion**: drug acts on a hazardous target or mechanism class (`config/cin_exclusions.json`: antimitotics, Aurora/PLK1/MPS1/CDK1, PP2A, proteasome, HSP90, topoisomerase, DNA-synthesis, PARP).
3. **Direction filter**: a drug must push each modifier the way the mechanism needs; for a loss-of-function gene, inhibitors of the gene or its network partners are rejected.
4. **Whole-drug safety pass** (plus a WHO ATC L01 chemotherapy-by-name pass and curated wrong-direction mechanisms for unannotated research compounds; investigational compounds of unknown direction are not proposed): every surviving drug is re-checked against its *full* Open Targets mechanism profile, so a hazardous drug can't slip in through a side annotation (this catches colchicine, which DGIdb lists against CREBBP).
5. **Not-proposable genes**: modifiers whose inhibition is itself unsafe (CREBBP/EP300) are used only to exclude activators.
6. **Curated layer**: agents the databases cannot represent (e.g. NAD⁺ precursors) are added from the literature and labelled as such.
7. **Score** = 0.5·clinical stage + 0.3·gene weight (halved if direction unknown) + 0.15 if direction confirmed + small publication bonus.

## Key result for BUB1B

Live run, 19 Oct 2026: 552 drug–gene records; 474 excluded as harmful or wrong direction, 28 CBP/p300 inhibitors and 49 investigational unknown-direction compounds set aside; **0** approved drugs act in the right direction on the BUBR1 stability axis. BUB1B's only record is a chemotherapy prodrug (excluded). The lead (nicotinamide riboside) therefore comes from the curated literature layer, after systematic search and filtering.

## Adapting to a new disease

Copy `config/bub1b.json`, change `gene`, `modifiers` (`want: up|down`, `weight`, `why`), `curated_candidates`, and point `extends` at a hazard filter suited to that disease.
