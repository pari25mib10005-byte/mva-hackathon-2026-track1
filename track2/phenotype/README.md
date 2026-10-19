# Track 2, evidence line 2: phenotype-driven analysis

Asks: *which genes and shared pathways explain this patient's symptoms, how does the Track 1 gene fit, and which drugs act on that shared biology?*

```bash
# patient terms stay in a local, git-ignored file: one HP:xxxxxxx per line; family history marked "# family"
python3 phenotype.py --hpo patient_hpo.txt --gene BUB1B --related CEP57,TRIP13 \
                     --config ../pipeline/config/bub1b.json --out results/
# public demo:
python3 phenotype.py --hpo example_public_mva1.txt --gene BUB1B --related CEP57,TRIP13 \
                     --config ../pipeline/config/bub1b.json --out results/
```

Windows: use `py -X utf8` instead of `python3`. First run downloads the HPO gene–phenotype table (~65 MB) into `data/`.

## Method
1. **Gene scoring**: each gene scores the summed information content (IC = ln(N genes / genes annotated to the term)) of the patient terms it explains, using HPO `phenotype_to_genes.txt`. Rare, specific terms outweigh common ones.
2. **Validation**: rank of the Track 1 gene using phenotype alone.
3. **Coverage**: which patient terms the disease genes explain; unexplained terms are flagged for clinicians.
4. **Module**: STRING (combined score ≥ 0.7) network among the top 200 phenotype genes; the causal gene's direct partners are the shared-biology module.
5. **Enrichment**: Reactome / KEGG / GO Biological Process via the STRING API (FDR < 1e-3).
6. **Drugs**: Open Targets + DGIdb on module genes, with the same exclusion, direction and whole-drug safety filters as `../pipeline/repurpose.py`. Links with no stated action (e.g. a gene used only as a trial biomarker) are set aside in `module_unclear_direction.csv`, not counted as evidence.

**Genetic-or-not test:** for each patient term, a hypergeometric test asks whether the causal gene's 100 strongest STRING partners (score ≥ 0.7, chosen independently of the patient) carry that term more than chance; terms are classified as genetic (annotated and/or network-enriched) or not explained by the disease network.

Outputs: `phenotype_summary.md`, `gene_scores.csv`, `module_drugs.csv`, `module_excluded.csv`, `module_unclear_direction.csv`.
