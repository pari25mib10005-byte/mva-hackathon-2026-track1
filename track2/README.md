# Track 2 — Drug repurposing for a *BUB1B* mosaic variegated aneuploidy case (Team PHA)

Report: [`pari25mib10005_track2_report.md`](pari25mib10005_track2_report.md). Everything here is a **hypothesis for investigation, not treatment advice.**

**No case data are included in this repository.** Variant identifiers and the patient's clinical features are withheld; patient inputs (VCF, HPO list) are read from local, git-ignored files only.

## Three independent routes
| Route | Folder | Starts from | Question |
|---|---|---|---|
| 1. Gene-first | `pipeline/` | Disease gene + mechanism config | Which drugs act on the gene, its modifiers and partners, in the right direction, safely? |
| 2. Phenotype-first | `phenotype/` | Patient HPO terms (local) | Which genes and pathways explain the features; which features are genetic; which drugs act on that biology? |
| 3. Transcriptome-first | `lincs/` | LINCS L1000 knockout signature | Which compounds reverse the effect of losing the gene, without silencing the p53/interferon alarm? Is the NAD⁺ hypothesis supported? |
| Local tool | `tools/` | Patient VCF (local) | Second-diagnosis check for features the disease network does not explain |

## Reproduce (Python 3.9+, standard library only; Windows: `py -X utf8`)
```bash
cd pipeline  && python3 repurpose.py --config config/bub1b.json --out results/
cd ../phenotype && python3 phenotype.py --hpo example_public_mva1.txt --gene BUB1B --related CEP57,TRIP13 \
                   --config ../pipeline/config/bub1b.json --out results/          # use your local patient_hpo.txt
cd ../lincs  && python3 signature_reversal.py --genes BUB1B,CEP57,TRIP13 --config ../pipeline/config/bub1b.json --out results/
python3 ../tools/check_vcf_genes.py --vcf /path/to/patient.vcf.gz          # LOCAL ONLY; share gene-level lines only
```
Databases are queried live (DGIdb, Open Targets, STRING, HPO, LINCS via Enrichr, WHO ATC via ChEMBL); counts drift slightly between runs, conclusions have not. Each run writes a provenance file.

## Safety rules shared by all routes
Antimitotics and mitotic-kinase inhibitors, PP2A/proteasome/HSP90 inhibitors, genotoxic agents, every WHO ATC L01 antineoplastic, wrong-direction agents (sirtuin, NAMPT, JAK, p53 inhibitors; CBP/p300 activators), CBP/p300 inhibitors (unsafe), investigational compounds of unknown direction, and, for signature reversal, compounds that switch off the p53/interferon alarm.

Licence: CC BY 4.0 (per hackathon rules).
