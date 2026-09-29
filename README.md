# MVA Hackathon 2026 — Track 1 (Variant Prediction)

Individual entry for *Rare Disease, Real Kid: The MVA Hackathon 2026*
(Sage Bionetworks, MVA Society, Hugging Face, BEACON).

**Result:** two rare heterozygous coding variants in *BUB1B* (p.Leu737Ter, ClinVar pathogenic;
p.Asn1002Lys, kinase-domain missense absent from gnomAD), submitted as an unphased
compound-heterozygous candidate for MVA syndrome 1. Full reasoning and limitations: [`report_track1.md`](report_track1.md).

## Reproduce

```bash
conda create -n mva -c conda-forge -c bioconda --override-channels bcftools python=3.11 -y
conda activate mva && pip install -r requirements.txt
bash run.sh /path/to/WGS_EX2312012_HGWCNDSX7.vcf.gz
```

Prints the top-ranked variants and writes `annotated.csv` (git-ignored).

## Files

| File | Purpose |
|---|---|
| `mva_genes.bed` | GRCh38 intervals for BUB1B, CEP57, TRIP13 (±10 kb), unprefixed to match the VCF |
| `run.sh` | Extraction, QC filter, multi-allelic split |
| `annotate.py` | Ensembl VEP REST annotation and ranking |
| `report_track1.md` | Methods, evidence, differential, limitations, AI disclosure |

## Data governance

The dataset is gated (WCG IRB #20252010). **No subject data is committed** — the VCF, all
intermediate VCFs, `annotated.csv` and the submission CSV are excluded via `.gitignore`.

## License

CC BY 4.0, per hackathon rules.
