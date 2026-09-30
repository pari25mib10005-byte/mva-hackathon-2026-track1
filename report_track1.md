# Track 1 — Variant Prediction: Methods and Evidence Report

**Hackathon:** Rare Disease, Real Kid: The MVA Hackathon 2026 (Sage Bionetworks, MVA Society, Hugging Face, BEACON)
**Track:** 1 — Variant prediction
**Participant:** **Team:** PHA — Pari, Harshita , Aradhya
**Proband:** PROBAND01 (sample `WGS_EX2312012`, single-sample WGS, GRCh38)
**Date:** 30 September 2026

## 1. Summary

The proband carries two rare, heterozygous, protein-altering variants in *BUB1B*, the gene responsible for Mosaic Variegated Aneuploidy syndrome 1 (MVA1; OMIM 257300):

| Rank | Variant (GRCh38) | Transcript effect | Protein | Genotype | gnomAD genomes AF | ClinVar |
|---|---|---|---|---|---|---|
| 1 | chr15:40209701 T>G | c.2210T>G, stop-gained | p.Leu737Ter | 0/1 | 3.3 × 10⁻⁵ | Pathogenic / Likely pathogenic |
| 2 | chr15:40220612 T>G | c.3006T>G, missense | p.Asn1002Lys | 0/1 | not observed | no record |

Submitted as one paired prediction (compound-heterozygous candidate). MVA1 is autosomal recessive and the clinical phenotype matches; a truncating allele plus a rare kinase-domain missense is the configuration reported in previously published MVA1 cases. This is a strongly supported genetic hypothesis, not a confirmed molecular diagnosis: with a single sample the two variants cannot be shown to be in *trans* (see §6).

## 2. Data

- `WGS_EX2312012_HGWCNDSX7.vcf.gz` + `.tbi` from the gated dataset `SageBio/mva-hackathon-2026-data`. GATK-style single-sample calls, ~5 million records, contigs unprefixed (`15`, not `chr15`).
- Build verified as GRCh38 from the header (`contig=<ID=15,length=101991189>`), not assumed from the filename.
- Raw FASTQ was not used; variants were already called.
- The clinical phenotype document supplied with the dataset.

## 3. Pipeline

All steps are in `run.sh` and `annotate.py` in the accompanying repository and can be re-run from the VCF in a few minutes.

1. **Candidate-gene extraction.** The three established MVA genes (*BUB1B* MVA1, *CEP57* MVA2, *TRIP13* MVA3) ± 10 kb were extracted with `bcftools view -R mva_genes.bed` (191 records).
2. **Quality filter.** `QUAL ≥ 30`, `FMT/DP ≥ 10`, `FMT/GQ ≥ 20`. Multi-allelic records were split with `bcftools norm -m -any` (184 records after filtering).
3. **Annotation.** Ensembl VEP REST API (`/vep/human/region`, VCF-style input, canonical transcript, `pick=1`) with gnomAD genome frequencies, SIFT, PolyPhen, ClinVar clinical significance and HGVS notation.
4. **Ranking.** Sort by VEP IMPACT (HIGH → MODERATE → LOW → MODIFIER), then ascending gnomAD allele frequency (absent-from-gnomAD ranked above observed).
5. **Inheritance reasoning.** For an autosomal-recessive disorder, a single heterozygous hit is insufficient; candidates were required to be homozygous or to occur as two rare coding hits in the same gene.

## 4. Results

Of 184 quality-filtered variants across the three genes, exactly **two** are protein-altering, both in *BUB1B*, both heterozygous, both rare. Every other variant is intronic or upstream with no predicted functional consequence and no ClinVar assertion.

- **p.Leu737Ter** — a premature stop in exon 18 of 23; predicted to trigger nonsense-mediated decay (null allele). Already classified Pathogenic/Likely pathogenic in ClinVar for MVA1 by multiple submitters.
- **p.Asn1002Lys** — lies in the C-terminal kinase domain of BUBR1. Absent from gnomAD genomes; a polar-to-basic substitution at a residue conserved across vertebrates. Not in ClinVar. (Note: a *different* nucleotide change, c.3006T>A, also produces p.Asn1002Lys and is in ClinVar as VUS; no evidence transfers between the two.)

Allele balance at both sites is close to 0.5 (from FORMAT/AD), consistent with germline heterozygosity rather than mosaicism or a copy-number artefact at these positions.

## 5. Differential diagnosis

| Gene | Disorder | Finding |
|---|---|---|
| *CEP57* | MVA2 | Intronic variants only; no coding, splice-site or ClinVar-flagged variant |
| *TRIP13* | MVA3 | Intronic variants only |
| *BUB1B* | MVA1 | Two rare coding hets (above) |

Within *BUB1B*, no homozygous coding variant and no splice-site variant was found, so a homozygous model is not supported.

## 6. Limitations

1. **Phase is unresolved.** The two variants are ~10.9 kb apart; short-read fragments cannot span them, and no parental samples are available. Compound heterozygosity is therefore inferred from the recessive model and the phenotype, not observed. Trio genotyping or long-read sequencing would resolve this.
2. **p.Asn1002Lys is a variant of uncertain significance** on formal criteria: rare, conserved and domain-located, but with no functional or ClinVar evidence. If it is not the second allele, the true second hit could be a structural variant or deep-intronic variant not called in this VCF.
3. **Copy-number and structural variants** were not assessed; the VCF contains SNV/indel calls only.
4. The primary analysis was candidate-gene driven; a gene-agnostic genome-wide pass (e.g. Exomiser with the proband's HPO terms) was not run.
5. No experimental validation was performed.

## 7. Mosaic aneuploidy scan from the proband's own WGS

MVA is defined cytogenetically by mosaic aneuploidy. I tested whether any aneuploidy is visible in the bulk blood WGS itself, using two independent signals per autosome: (a) sequencing depth relative to the genome median, and (b) the B-allele fraction (BAF) of high-quality heterozygous SNPs (`GT=het`, `DP ≥ 20`, `GQ ≥ 30`, `QUAL ≥ 50`; ~2.2 million sites). In a mosaic gain, depth rises and het BAF disperses away from 0.5. Scripts: `mosaic_check.py`, `baf_hist.py`, `baf_bins.py`.

**First pass (whole-chromosome averages) produced a false positive.** Chromosomes 20, 21 and 22 showed elevated depth (1.07–1.11×), increased BAF dispersion (robust z = 4.1, 5.6, 6.0) and a downward mean BAF shift (0.466, 0.444, 0.436 vs 0.487 on chr1). Taken alone this resembles mosaic trisomy at roughly 10–20 %.

**Binning at 5 Mb resolution showed the signal is not chromosomal.** It is confined to chr20 25–35 Mb, chr21 5–15 Mb, chr22 10–20 Mb and chr1 120–145 Mb (depth 1.4–2.1×, mean BAF 0.30–0.40), i.e. the centromeres, the acrocentric short arms of 21 and 22, the 1q21 segmental duplications and the 22q11 low-copy repeats. These are canonical short-read mapping artefacts: reads from paralogous copies pile up and skew allele counts. Every euchromatic bin on all four chromosomes sits at depth 1.00–1.04× and mean BAF 0.49–0.50, exactly matching chr1. A genuine mosaic gain would lift every bin uniformly. Chromosomes 21 and 22 are small, so their repeat regions dominated the whole-chromosome average; the identical artefact on chr1 was diluted by 240 Mb of clean sequence.

**After excluding 5 Mb bins with depth > 1.15× the genome median**, chr20 and chr21 fall to the lowest dispersion of all autosomes (z = −1.9, −2.1). A residual on chr22 (z = 6.6) traces to the 22q11 bins at 15–20 Mb, which pass the depth cut but not the BAF check; chr22 from 25 Mb onward is indistinguishable from chr1. Chromosomes 16 and 19 show depth of 1.045× with no BAF change, consistent with GC-content coverage bias rather than aneuploidy.

**Conclusion:** no mosaic aneuploidy is detectable in this bulk blood WGS. This is the expected result rather than a negative finding about the diagnosis: "variegated" means different cells carry different aneuploidies, so no single chromosome accumulates a consistent fraction, and bulk sequencing averages the signal to zero. Detection requires per-cell karyotyping or single-cell sequencing. I note that whole-chromosome summary statistics on chr20/chr22 can be mistaken for low-level mosaic gains; the per-bin check above is the safeguard. Per-chromosome and per-bin aggregates (no genotype-level data) are in `mosaic_per_chrom.csv` and `mosaic_bins.csv`.

## 8. Reproducibility

- Tools: bcftools (via bioconda), Ensembl VEP REST API (GRCh38, default release at run date), Python 3.11, `requests`.
- Repository contains the BED file, `run.sh`, `annotate.py`, the three mosaic-scan scripts, their aggregate outputs and this report. It contains **no subject data**: the VCF, intermediate VCFs, the annotated variant table and the submission CSV are all git-ignored, as the latter three contain proband genotypes.
- Submission coordinates were converted from the VCF's unprefixed contigs to the `chr`-prefixed form required by the scorer; `proband_id` is `PROBAND01`.

## 9. AI-assistance disclosure

Anthropic Claude (claude.ai, consumer plan) was used interactively to design the pipeline, debug environment setup, write the annotation and mosaic-scan scripts, interpret the mosaic-scan output (including recognising the first-pass false positive as a repeat-region artefact) and draft this report. No proband data was uploaded to the AI service; only tool output that had already been reduced to the two candidate variants and aggregate counts was discussed. All commands were run locally by the participant and all claims above were checked against the local outputs.

## 10. Acknowledgement

This work was made possible through the Hackathon organised by Sage Bionetworks in partnership with the MVA Society, Hugging Face and BEACON, with prize sponsorship from AWS and Anthropic. Thanks to the child and family who contributed their data to advance research into this rare disease.
