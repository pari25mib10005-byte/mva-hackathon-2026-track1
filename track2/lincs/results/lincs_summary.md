# LINCS signature reversal: BUB1B

> Hypotheses for investigation, not treatment advice.

- Knockout signature (BUB1B): 247 genes up, 244 genes down
- Compounds scored: 5425; gene universe: 9759
- Significant reversers (score > 0, FDR < 0.05): 1134
- Knockout signatures available in LINCS: BUB1B (not available: CEP57, TRIP13)

## Pre-specified test: NAD⁺/sirtuin axis

| Compound | Rank (of 5425) | Score | FDR |
|---|---|---|---|
| SRT-1720 | 283 | 0.077 | 0.00069 |
| Resveratrol | 1295 | 0.039 | 0.099 |
| Nadide | 1683 | 0.029 | 0.15 |
| Nicotinamide | 3256 | -0.002 | 0.97 |

One-sided Mann–Whitney (axis > others): U = 15142, z = 1.37, **p = 0.0848** (n = 4 axis compounds found of 7 pre-specified).

## Top reversers after safety review

- Reviewed: 100; excluded by safety rules: 94; remaining: 6; of which approved: 2

| Rank | Compound | Score | FDR | Stage | Mechanism |
|---|---|---|---|---|---|
| 63 | Prima-1-Met | 0.102 | 1.7e-08 | n/a | n/a |
| 73 | Exalamide | 0.100 | 6.2e-07 | UNKNOWN | n/a |
| 81 | Penicillin | 0.098 | 1.3e-06 | APPROVAL | Bacterial penicillin-binding protein inhibitor |
| 82 | Pempidine | 0.098 | 1.3e-06 | UNKNOWN | n/a |
| 97 | Ro-10-5824 | 0.096 | 1.9e-06 | n/a | n/a |
| 100 | Amineptine | 0.096 | 1.9e-06 | APPROVAL | Dopamine D1 and D2 receptor antagonist |

**Alarm check:** 93 of the top 100 reversers switch off ≥3 p53/interferon alarm genes (18 alarm genes are in the knockout 'up' half).

## Exploratory: restore division genes without silencing the alarm

Compounds that raise the knockout's 'down' genes (FDR < 0.05) and switch off < 3 alarm genes; after safety review: 93 of 100; approved: 27.

| Compound | Division genes restored | FDR | Stage | Mechanism |
|---|---|---|---|---|
| Procaine | 32 | 6.6e-13 | APPROVAL | Sodium channel alpha subunit blocker |
| Amisulpride | 32 | 9.6e-13 | APPROVAL | Serotonin (5-HT) receptor antagonist; Dopamine receptors; D2 |
| Bentazepam | 31 | 2.5e-12 | APPROVAL | n/a |
| Leucovorin | 31 | 3.6e-12 | APPROVAL | n/a |
| Choline | 29 | 8.7e-11 | APPROVAL | n/a |
| Levonorgestrel | 29 | 9.5e-11 | APPROVAL | Progesterone receptor agonist |
| Hydrastine | 29 | 1e-10 | APPROVAL | n/a |
| Haloperidol | 29 | 1e-10 | APPROVAL | Serotonin 2a (5-HT2a) receptor antagonist; D2-like dopamine  |
| Meloxicam | 29 | 1e-10 | APPROVAL | Cyclooxygenase-2 inhibitor |
| Resorcinol | 28 | 4.3e-10 | APPROVAL | n/a |
| Trientine | 28 | 4.6e-10 | APPROVAL | n/a |
| Medroxyprogesterone-Acetate | 28 | 4.6e-10 | APPROVAL | Progesterone receptor agonist |
| Methimazole | 27 | 1.7e-09 | APPROVAL | Thyroid peroxidase inhibitor |
| Atrasentan | 27 | 2e-09 | APPROVAL | Endothelin receptor ET-A antagonist |
| Ribostamycin | 27 | 2e-09 | APPROVAL | n/a |

Excluded among top reversers: 87 × switches off the p53/interferon alarm; 3 × wrong direction or unsafe class; 3 × antineoplastic; 1 × acts on excluded target(s) AURKA
