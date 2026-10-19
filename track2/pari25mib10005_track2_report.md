# Track 2 — Mechanism and Drug Repurposing for a *BUB1B* Mosaic Variegated Aneuploidy Case

**Team PHA** (Pari, Aradhya, Harshita) · Rare Disease, Real Kid: MVA Hackathon 2026
Code and data provenance: `track2/pipeline/` in [github.com/pari25mib10005-byte/mva-hackathon-2026-track1](https://github.com/pari25mib10005-byte/mva-hackathon-2026-track1)

> **Everything in this report is a hypothesis for investigation, not treatment advice.** Any use in the child must be decided by his clinicians under appropriate ethical oversight.

---

## In plain language

This child has a rare condition, mosaic variegated aneuploidy (MVA), in which some of his cells end up with too many or too few chromosomes when they divide. We found that it is caused by two faulty copies of a gene called *BUB1B*. That gene makes BUBR1, a protein that acts like a safety check, telling a cell to wait until its chromosomes are sorted correctly before dividing. His cells still make some BUBR1, but too little, because it is broken down too quickly.

Working backwards from his symptoms alone, our analysis independently pointed to the same gene and the same faulty cell-division safety check. We then looked for a medicine that already exists and might help his cells keep more of the BUBR1 they make. Our best idea is **nicotinamide riboside**, an inexpensive form of vitamin B3 that has been given to children with another rare condition. In mice, this approach increased BUBR1. It has **not** been tested for MVA. It might reduce new chromosome errors over time, but it cannot undo features that formed before birth. We set out how it should be tested in the lab before anyone considers it for him.

---

## Summary

Our Track 1 analysis identified compound heterozygous *BUB1B* variants (one truncating, one missense; exact variant identifiers are withheld from public materials). We argue that this genotype produces a **BUBR1 hypomorph**: the child's cells make some BUBR1, but too little, because one allele is null and the other encodes an unstable protein. The therapeutic goal is therefore to **stabilise the BUBR1 that remains**, while avoiding anything that weakens the mitotic checkpoint further or adds cancer risk.

A systematic, gene-agnostic screen of DGIdb and Open Targets across 29 related genes found **no usable drug that targets BUB1B** and **no approved drug known to push the BUBR1-stability pathway in the right direction**; all 552 drug–gene records were excluded as harmful, wrong-direction, unsafe or investigational-with-unknown-direction, except one approved compound of unknown direction. An independent phenotype-driven analysis, starting from the patient's clinical features rather than his DNA, ranked *BUB1B* 16th of 2,479 genes and converged on the same checkpoint pathway, with no directly druggable target in the helpful direction. A third, pathway-agnostic route (LINCS L1000 signature reversal) tested the NAD⁺ hypothesis with a pre-specified statistic: support was weak and not significant (p = 0.085), and the analysis showed that naive reversal of a chromosome-instability signature favours drugs that silence the cell's protective p53/interferon alarm (93 of the top 100), which we exclude. The strongest candidate is therefore taken from the literature: **nicotinamide riboside (NR)**, an NAD⁺ precursor. It acts on the NAD⁺ → SIRT2 → BUBR1 axis that raised BubR1 levels in mice, and it has open-label safety data in children with another DNA-repair/ageing disorder. We expect at most a partial benefit and propose a staged validation plan before any clinical use.

---

## 1. Mechanism

**Mechanism classification:** **loss-of-function (hypomorphic)**. Pathway disrupted: **spindle assembly checkpoint and kinetochore–microtubule attachment**. Downstream consequence: **chromosome missegregation → mosaic random aneuploidy → growth/developmental features and cancer predisposition**.

### 1.1 Starting point

Track 1 (scored 100/100; F-max 1.000) identified two *BUB1B* variants in trans:

| Allele | Variant | Type | Predicted effect |
|---|---|---|---|
| A | Truncating (premature stop codon) | Nonsense | Null (nonsense-mediated decay) |
| B | Single amino-acid substitution | Missense, kinase domain | Full-length but unstable protein |

Biallelic *BUB1B* mutations, typically one truncating plus one missense allele, are an established cause of MVA [2]. The mechanism below is therefore anchored in a confirmed genotype, not a guessed one.

### 1.2 What BUBR1 does

BUBR1 is a core component of the spindle assembly checkpoint (SAC). It has two jobs:

1. **Checkpoint signalling.** Its N-terminal KEN boxes and TPR region bind CDC20 and, within the mitotic checkpoint complex, keep APC/C–CDC20 inactive until every kinetochore is attached [3].
2. **Attachment stability.** A central region recruits PP2A-B56 to kinetochores, stabilising correct kinetochore–microtubule attachments [3].

The C-terminal (pseudo)kinase domain is largely dispensable for catalysis in vertebrates, but its structural integrity prevents degradation of the whole protein [4].

### 1.3 Predicted effect of each allele

**Allele A (truncating).** Truncating *BUB1B* alleles in MVA patients give no detectable transcript, consistent with nonsense-mediated decay [1]. Any truncated protein that escaped NMD would lack the stabilising C-terminal domain [4]. We model this allele as null.

**Allele B (kinase-domain missense).** This introduces a charge-altering substitution inside the kinase domain. In patient cells, MVA missense mutations in or near this domain consistently caused rapid BUBR1 turnover rather than loss of a specific function [1]. We predict a full-length but unstable protein.

### 1.4 Combined model: a BUBR1 hypomorph

```
Allele A (null, NMD) + Allele B (unstable protein)
        ↓
Low total BUBR1 protein (reduced, not absent)
        ↓
Weakened SAC + less stable kinetochore–microtubule attachment
        ↓
Chromosome missegregation in dividing cells
        ↓
Mosaic random aneuploidy → growth/developmental features + cancer predisposition
```

This places the case in the "low protein abundance" class of MVA mutations [1]. Reintroducing BUBR1 into MVA patient cells restored checkpoint activity [1], so the cellular defect is in principle reversible by raising BUBR1 levels. The BubR1 hypomorphic mouse (BubR1^H/H) is the closest animal model: low BubR1, aneuploidy, progeroid features and shortened lifespan [7].

### 1.5 From mechanism to drug strategy

No drug targets BUBR1 directly, so we looked for ways to **slow its degradation**. The key published observation [5]:

- With age, BubR1 is lost because NAD⁺ falls and SIRT2 can no longer keep **lysine-668** deacetylated; the acetyltransferase **CBP** acetylates this residue, priming BubR1 for degradation.
- Raising NAD⁺ with **nicotinamide mononucleotide (NMN)**, or overexpressing SIRT2, **increased BubR1 abundance in vivo**. SIRT2 overexpression also extended median lifespan of BubR1^H/H mice [5, 6].
- SIRT2 did not restore normal BubR1 levels or remove all progeroid features [5, 6], so any rescue should be expected to be partial.

This defines the **NAD⁺ → SIRT2 → BUBR1 stabilisation axis** as the primary search space, and gives four design rules used throughout:

1. Raise NAD⁺ / SIRT2 activity, or otherwise slow BUBR1 turnover.
2. Never weaken the checkpoint or increase missegregation.
3. Avoid genotoxic agents (MVA predisposes to cancer [2]).
4. Prefer agents with paediatric safety data.

---

## 2. Systematic candidate search (reproducible pipeline)

To show the candidate list is not hand-picked, we built a **gene-agnostic pipeline** (`track2/pipeline/repurpose.py`, Python standard library only). It works for any gene through a config file; we tested it unchanged on the other MVA genes, *TRIP13* and *CEP57*.

### 2.1 Method

1. **Gene set (29 genes):** *BUB1B*; curated modifiers with the direction the mechanism needs (SIRT2 ↑, NAMPT ↑, CREBBP ↓, EP300 ↓); and high-confidence STRING neighbours (score ≥ 0.9) from Open Targets.
2. **Evidence:** drug–target records from **DGIdb** and **Open Targets** (552 raw records; live query run on 19 October 2026).
3. **Hard exclusion** by target or mechanism class (shared config `cin_exclusions.json`): antimitotics/tubulin, Aurora, PLK1, MPS1/TTK, CDK1, KIF11/CENPE, PP2A, proteasome, HSP90, topoisomerases, DNA-synthesis inhibitors, PARP.
4. **Direction filter:** a drug must move each modifier in the needed direction; for a loss-of-function gene, inhibitors of the gene or its network partners are rejected.
5. **Whole-drug safety pass:** every surviving drug is re-checked against its *full* Open Targets mechanism profile, so a hazardous drug cannot enter through a side annotation.
6. **Chemotherapy-by-name pass:** every drug in WHO ATC group L01 (antineoplastic agents, via ChEMBL) is excluded by name, so a chemotherapy cannot enter when a database record lacks its mechanism or uses a variant name (e.g. a prodrug).
7. **Curated wrong-direction mechanisms:** research compounds that the databases do not annotate (e.g. sirtuin and NAMPT inhibitors) are excluded from a documented list in the config; investigational compounds with unknown direction are not proposed for a child.
8. **Not-proposable modifiers:** CREBBP/EP300 inhibition is itself unsafe (CREBBP loss causes Rubinstein–Taybi syndrome; both are tumour suppressors), so these genes are used only to exclude activators.
9. **Curated literature layer:** agents the databases cannot represent (metabolic precursors) are added and labelled as such.
10. **Outputs:** ranked candidates, a full exclusion log with reasons, and a provenance file (run time, sources, config, gene set).

### 2.2 Results

| Measure | Result |
|---|---|
| Drug records on *BUB1B* itself | **1** (a gemcitabine prodrug from a cell-sensitivity screen; excluded as chemotherapy) |
| Raw drug–gene records across the gene set | 552 |
| Excluded as harmful or wrong direction | **474** |
| Set aside: CBP/p300 inhibitors (right direction, unsafe) | 28 |
| Set aside: investigational compounds with unknown direction | 49 |
| Database hits with a confirmed, correct direction on the BUBR1 axis | **0** |
| Approved compounds with unknown direction (manual review) | 1 |

Largest exclusion groups: PLK1 inhibitors (199), CDK1 (109), Aurora B (103), NAMPT inhibitors (11, would lower NAD⁺), PP2A inhibitors (7), sirtuin inhibitors (4). Successive safety passes caught agents that would otherwise have slipped through: **colchicine** (a tubulin poison listed in DGIdb against CREBBP), Aurora inhibitors listed against NEK2, and **four cancer drugs** (a gemcitabine prodrug, gefitinib, palbociclib, pazopanib) whose records lacked a usable mechanism. Counts come from live databases and drift between runs; earlier runs gave 513 records and the same conclusion.

**Interpretation.** Drug–target databases contain no approved drug that stabilises BUBR1, and the obvious "related" drugs are almost all harmful in this setting. The NAD⁺ precursors are absent because they are metabolic substrates, not target binders. The lead candidate therefore comes from the literature layer, after systematic search and filtering.

---

## 3. Phenotype-driven analysis (independent evidence line)

Section 2 starts from the gene. This analysis starts from the **patient's clinical features** (HPO terms, used locally and not reproduced here) and asks which genes and shared pathways explain them, without using any sequencing data (`track2/phenotype/phenotype.py`).

**Method.** Every gene in the HPO gene–phenotype table (5,276 genes) was scored by the summed information content of the patient terms it explains, so specific features outweigh common ones. A STRING network (combined score ≥ 0.7) was built over the top 200 phenotype genes plus the disease genes; the causal gene's direct partners define the shared-biology module, tested for pathway enrichment (Reactome, KEGG, GO). Drugs acting on the module were retrieved and passed through the same exclusion, direction and whole-drug filters as Section 2. One family-history term was excluded from the patient's phenotype.

| Result | Value |
|---|---|
| Rank of *BUB1B* from phenotype alone | **16 of 2,479** matching genes (top 0.6%); *CEP57* 17th, *TRIP13* 104th |
| Direct network partners of *BUB1B* among phenotype genes | BUB1, BUB3, TRIP13, PPP2R5C (PP2A-B56), TP53 |
| Top enrichment of the module | KEGG cell cycle (FDR 5.9×10⁻¹⁰); GO mitotic spindle assembly checkpoint signalling (FDR 2.4×10⁻⁷); Reactome amplification of signal from unattached kinetochores (FDR 7.0×10⁻⁶) |
| Patient features explained by MVA genes | 3 of 7; **4 of 7 are not annotated to any MVA gene** |
| Drugs acting on the module with known action, correct direction and passing safety filters | 4, all via TP53, none proposable (below) |
| Association-only links set aside (e.g. TP53 used as a trial biomarker) | 466 |

**Interpretation.**
- The phenotype alone recovers the Track 1 gene near the top of 2,479 candidates, an independent validation of the genetic diagnosis.
- The biology shared between the patient's features and *BUB1B* is the **spindle assembly checkpoint**, including PP2A-B56, which matches the mechanism in Section 1.
- The four remaining module drugs do not fit this patient: eprenetapopt restores *mutant* p53 (no TP53 defect is implicated here), two are cancer vaccines, and one is a database artefact. The checkpoint core therefore cannot be drugged directly in the helpful direction, which independently supports the indirect, stabilisation strategy.
- Four of the patient's features are not explained by MVA gene annotations. Section 3.1 tests whether they are genetic.

### 3.1 Genetic or not: which features does the MVA network explain?

**Method.** *BUB1B*'s 100 strongest STRING partners (score ≥ 0.7), selected without reference to the patient's features (30 carry HPO annotations), were tested for each feature with a hypergeometric test (universe 5,276 HPO genes). Features are reported here by category only.

| Feature category (n) | Annotated to an MVA gene | Enrichment among *BUB1B* partners | Classification |
|---|---|---|---|
| Tumour (1) | Yes | 4 vs 0.1 expected, **p = 5×10⁻⁶** | Genetic (MVA) |
| Growth before and after birth (2) | Yes | 17 vs 8.4, **p = 8×10⁻⁴**; 3 vs 1.7, p = 0.25 | Genetic (MVA) |
| Perinatal (1) | No | 0 vs 1.2, p = 1.0 | Likely secondary: growth-restricted fetuses are often delivered early |
| Nutritional growth after birth (1) | No | 4 vs 6.7, p = 0.93 | Likely secondary: growth disorder plus illness |
| Muscle (1) | No | 0 vs 3.4, p = 1.0 | Likely acquired: nutrition, the cancer itself, or treatment (vincristine, standard for this tumour type, causes neuropathy and muscle wasting) |
| Renal (1) | No | 0 vs 0.6, p = 1.0 | Acquired (e.g. loop diuretics in premature infants, some chemotherapies) **or a second diagnosis** |

**The MVA network explains the tumour and growth features strongly and the other four not at all**, so those four are secondary, acquired, or due to a second cause. "Secondary" and "acquired" are hypotheses from known causes; treatment history was not available to us.

**Possible second diagnosis.** The renal, perinatal and nutritional features together form the classic pattern of antenatal Bartter syndrome (a salt-wasting tubulopathy), and the phenotype scoring independently ranked two of its genes, *SLC12A1* and *KCNJ1*, first and second of all genes. Dual molecular diagnoses occur in a minority of solved cases [22]. `track2/tools/check_vcf_genes.py` checks the patient VCF **locally** for rare variants in *SLC12A1, KCNJ1, CLCNKB, BSND* and *MAGED2* and prints only gene-level verdicts. **Result (local run, 19 October 2026): no rare protein- or splice-altering variants in any of the five genes** (three rare non-coding variants, low prior). A second genetic cause of this tubulopathy type is therefore not supported by the small-variant calls, which favours an acquired explanation for the renal and perinatal features. Caveat: the VCF holds small variants only, so whole-gene deletions (a common cause of *CLCNKB*-related Bartter syndrome) and deep-intronic effects are not assessed; a copy-number review or a clinical tubular work-up remains reasonable if the history does not explain the findings. Clinical history (diuretics, chemotherapy, urinary calcium) should be reviewed first.

**Family history.** The excluded family-history term has no overlap with MVA genes or *BUB1B* partners (its HPO genes are mainly coagulation genes); possible explanations include chance, parental carrier status, or maternal factors. It is a matter for genetic counselling, not a drug target for the child.

---

## 4. Transcriptomic signature reversal (third, pathway-agnostic route)

Routes 1 and 2 both use prior knowledge of genes and pathways. This route asks, with **no pathway assumption**, which compounds push human cells in the opposite direction to losing *BUB1B* (`track2/lincs/signature_reversal.py`).

**Method.** The LINCS L1000 CRISPR-knockout consensus signature of *BUB1B* (247 genes up, 244 down; no signatures exist for *CEP57* or *TRIP13*) was compared with consensus signatures of 5,425 compounds [18, 19]. Score = (reversing − mimicking overlaps) / signature size; direction tested with a binomial sign test, Benjamini–Hochberg FDR. **Pre-specified hypothesis, fixed before results:** NAD⁺/sirtuin-axis compounds reverse the signature more than other compounds (one-sided Mann–Whitney U). Top reversers then passed the Section 2 safety rules, wrong-direction patterns (sirtuin, NAMPT, JAK and p53 inhibitors), and exclusion of every WHO ATC L01 antineoplastic [21], which catches chemotherapies whose database mechanism is missing (e.g. procarbazine, estramustine).

**What the knockout signature contains.** The "up" half includes the cell's **protective alarm** to chromosome errors (p53/apoptosis genes such as *CDKN1A*, *GDF15*, *BAX*; interferon genes such as *ISG15*, *MX1*, a known response to chromosome instability [20]); the "down" half includes **cell-division genes**, among them *BUB1B*, *BUB1*, *TRIP13*, *CDC20*, *CCNB1* and *PLK1*.

| Result | Value |
|---|---|
| Significant reversers (FDR < 0.05) | 1,134 of 5,425 |
| NAD⁺/sirtuin axis (pre-specified) | SRT-1720 rank 283 (FDR 7×10⁻⁴); resveratrol 1,295 (FDR 0.10); NAD⁺ (Nadide) 1,683 (FDR 0.15); nicotinamide 3,256 (FDR 0.97). **Group test p = 0.085 (not significant)**; NR and NMN are not in LINCS |
| Top 100 reversers that switch off ≥ 3 alarm genes | **93** (e.g. tofacitinib, a JAK inhibitor; EX-527, a sirtuin inhibitor; APO-866, a NAMPT inhibitor) |
| Top 100 remaining after all safety rules | 6 (2 approved), with no coherent mechanism |
| Exploratory: restore division genes without silencing the alarm | 100 at FDR < 0.05; 93 pass safety (27 approved), mechanistically heterogeneous (antipsychotics, local anaesthetics, hormones, NSAIDs) |

**Interpretation.**
- **The NAD⁺ hypothesis gains no significant transcriptomic support.** Only the synthetic SIRT1 activator SRT-1720 reverses significantly; nicotinamide does not, consistent with its sirtuin-inhibiting effect at high concentration (Section 7). Because L1000 measures transcripts, it cannot see NR's proposed mechanism (BUBR1 acetylation and protein turnover), so this is an absence of support, not evidence against.
- **Naive reversal is the wrong objective for chromosome-instability disorders.** Most strong reversers act by silencing the p53/interferon alarm, which would let aneuploid cells survive in a cancer-prone child. Reversal must be split into its two halves and alarm-suppressors excluded; we report this as a methodological finding.
- **No new candidate is promoted.** The safe, division-restoring compounds share no pathway, so we treat them as exploratory comparators for the cell screen (Section 10), not as leads.

---

## 5. Convergence of evidence

| Question | Route 1: gene-first (Section 2) | Route 2: phenotype-first (Section 3) | Route 3: transcriptome-first (Section 4) | Verdict |
|---|---|---|---|---|
| Is *BUB1B* the right gene? | Track 1 genotype | Ranked 16th of 2,479 from features alone | Knockout lowers checkpoint genes (*BUB1*, *TRIP13*, *CDC20*, *PLK1*) | **Consistent** |
| Which biology is disrupted? | Spindle checkpoint; PP2A-B56 attachment | Module enriched for checkpoint signalling | Loss of division genes plus p53/interferon alarm | **Consistent** |
| Can it be drugged directly in the right direction? | 0 direction-confirmed approved drugs | 0 proposable module drugs | No safe reverser with a coherent mechanism | **Consistent (no)** |
| Is the NAD⁺/sirtuin strategy supported? | Mechanism-based (literature) | Indirect strategy required | Weak, not significant (p = 0.085) | **Literature-based; untested in transcripts** |
| Which features can a drug address? | n/a | Tumour and growth features (genetic); four others secondary, acquired or a second cause | n/a | **Scope defined** |

Three independent starting points agree on the gene, the pathway and the conclusion that the checkpoint cannot be drugged directly. The choice of an indirect BUBR1-stabilising strategy rests on that agreement; the choice of NR within it rests on mechanism and paediatric safety, and remains a hypothesis for the validation plan (Section 10).

---

## 6. Exclusion filter

Agents predicted to worsen missegregation or cancer risk are excluded regardless of other merits.

| Excluded class | Examples | Reason |
|---|---|---|
| Antimitotics | Taxanes, vinca alkaloids, colchicine | Disrupt spindles on top of a weak checkpoint |
| Mitotic kinase inhibitors | Aurora, MPS1/TTK, PLK1, CDK1 inhibitors | Weaken the same checkpoint and error-correction machinery that is already deficient |
| KIF11/CENPE inhibitors | Ispinesib, filanesib | Cause mitotic arrest and errors |
| PP2A inhibitors | Cantharidin | BUBR1 recruits PP2A-B56 to stabilise attachments |
| Proteasome inhibitors | Bortezomib, carfilzomib | Would slow BUBR1 degradation, but are cytotoxic with broad mitotic effects; unacceptable in a child |
| HSP90 inhibitors | Investigational agents | Destabilise many client proteins; not selective |
| Sirtuin / NAMPT inhibitors; CBP/p300 activators | Daporinad, AGK2 | Predicted to speed BUBR1 acetylation and loss |
| CBP/p300 inhibitors | Inobrodib, GNE-049 | Directionally consistent but unsafe (tumour suppressors; Rubinstein–Taybi) |
| Genotoxic agents | Alkylators, topoisomerase poisons, PARP inhibitors | MVA predisposes to cancer |

---

## 7. Candidates

Evidence tiers: **A** = mechanism shown in a BubR1 model and paediatric human safety data; **B** = mechanism in a BubR1 model, limited paediatric data; **C** = indirect rationale only.

| # | Agent | Mechanistic link | Key evidence | Status and paediatric safety | Main concerns | Tier |
|---|---|---|---|---|---|---|
| 1 | **Nicotinamide riboside (NR)** | NAD⁺ precursor → SIRT2 deacetylation of BUBR1 K668 → less degradation | NAD⁺ precursor (NMN) raised BubR1 in mice [5] | Supplement. Open-label studies in ataxia-telangiectasia: n=24, mean age 17.5 [8]; 2-year study, n=10 incl. 7 children aged 3+ [9]; case starting at 3.5 years, no adverse effects [10] | Partial rescue at best; theoretical concern that raising NAD⁺ could support tumour growth in a cancer-prone child | **A** |
| 2 | Nicotinamide mononucleotide (NMN) | Same axis | Directly tested: increased BubR1 in vivo [5] | Little paediatric data; regulatory status varies by country | As NR, with less human safety data | **B** |
| 3 | Nicotinamide (vitamin B3) | NAD⁺ precursor | By analogy with NR/NMN | Approved, inexpensive, long clinical use | At high concentrations it *inhibits* sirtuins, which could work against the goal | **B/C** |
| 4 | N-acetylcysteine (NAC) | Antioxidant | Mice with BUBR1 at ~20% of normal had higher mortality under oxidative stress [11]; no antioxidant was tested | Approved, long paediatric use | Weak, indirect link; exploratory only | **C** |
| 5 | Senolytic (dasatinib + quercetin) | Removes senescent cells that accumulate after missegregation | *Genetic* clearance of p16⁺ senescent cells delayed ageing disorders in BubR1 hypomorphs [12]; D+Q itself not tested there | Dasatinib is a cytotoxic kinase inhibitor | Toxicity likely outweighs benefit. **Reserve only** | **C** |

**Lead hypothesis: nicotinamide riboside.** It is the only candidate combining a mechanism that acts directly on the defect (BUBR1 stability), supporting evidence in a BubR1 mouse model, and human safety data that include young children. NR (as nicotinamide riboside chloride) also received FDA orphan-drug and rare paediatric disease designations for ataxia-telangiectasia in June 2024, as announced by the manufacturer, showing an existing regulatory path for paediatric use. It is sold widely as a dietary supplement (the ingredient itself is patent-protected), which keeps access simple and supports scalability.

**What NR cannot do.** Even SIRT2 overexpression did not restore normal BubR1 levels in mice [5, 6]. The realistic aim is to **reduce the rate of missegregation**, not to cure the condition.

---

## 8. What this could and could not help

MVA type 1 typically involves slow growth before birth, small head size, developmental differences, and a raised risk of childhood cancer [15]. The proposed mechanism acts only on **cells that are still dividing**, so its potential reach is limited and should be stated plainly.

| Feature (general MVA picture) | Driven by | Could stabilising BUBR1 help? |
|---|---|---|
| Ongoing chromosome errors in dividing tissues (e.g. blood) | Current checkpoint weakness | **Possibly.** This is the direct target and the main measurable outcome. |
| Childhood cancer risk | Accumulating aneuploid cells over time | **Possibly, in the long term**; unproven. Tumour surveillance must continue unchanged. |
| Growth after birth | Partly ongoing cell division | **Uncertain**; at most a modest effect. |
| Growth restriction before birth, small head size, structural birth differences | Development that has already happened | **No.** A treatment started now cannot reverse these. |
| Developmental and neurological features | Mostly set during brain development | **Unlikely**; mature neurons do not divide. |

*The patient-specific phenotype analysis is summarised in Section 3; individual clinical features are not reproduced in this public report.*

---

## 9. Notes for clinicians

These are considerations to support a clinical team's own judgement, not a protocol. Dosing, if any, must be decided by the treating clinicians.

**Before any trial**
- Establish a baseline: percentage of aneuploid cells in blood (karyotype or interphase FISH, or single-cell sequencing where available), blood NAD⁺ metabolites, routine bloods.
- Confirm cell-level response first (Section 10, steps 1–2), ideally in the child's own cells.

**During any trial**
- Repeat the aneuploidy measurement at pre-set intervals; this is the primary readout.
- Continue the child's existing cancer surveillance exactly as before. Because raising NAD⁺ could in theory support tumour cell metabolism, any new mass or unexplained change should prompt stopping and review.
- Use a pre-specified stop rule (adverse effects, rising aneuploidy, or no change after a defined period).
- Prefer NR over high-dose nicotinamide: at high concentrations nicotinamide inhibits sirtuins, the opposite of the intended effect.

**Drug classes to use with caution in MVA generally**
Agents that disrupt cell division or damage DNA (Section 6) may worsen chromosome errors in a child whose checkpoint is already weak. Where alternatives exist for other indications, clinicians may wish to weigh this.

---

## 10. Validation plan

1. **Cell model.** In MVA patient-derived cells (or *BUB1B*-knockdown lines), measure BUBR1 protein level, K668 acetylation and missegregation rate (live imaging, micronucleus counts) with and without NR across a dose range.
2. **Dose–response and specificity.** Confirm any BUBR1 increase is dose-dependent and SIRT2-dependent (abolished by SIRT2 knockdown), and check that high-dose nicotinamide does not show the opposite effect.
3. **Animal model.** Test NR in BubR1^H/H mice for BubR1 level, aneuploidy rate and tumour incidence over time.
4. **Only then, n-of-1.** Under the child's clinicians and an ethics framework: baseline and on-treatment measurement of the percentage of aneuploid cells in blood (karyotype or single-cell sequencing), plus safety monitoring, with a pre-specified stop rule.

---

## 11. Impact, innovation and scalability

- **Impact.** The lead is cheap, oral, widely available and has paediatric safety data, so it could be tested quickly if the cell work is positive. The exclusion table is useful on its own: it tells clinicians which commonly used drug classes to be cautious with in MVA.
- **Innovation.** The strategy targets protein *stability* (the actual defect in this genotype class), not the gene. Three independent routes (gene, phenotype, transcriptome) are combined, each feature of the patient is classified as genetic or not with a network test, and the transcriptomic route shows that reversal must exclude drugs that silence the p53/interferon alarm. The pipeline combines direction-aware filtering, a whole-drug safety pass and explicit labelling of literature-only evidence.
- **Scalability.** Both tools are disease-agnostic: the gene-first screen runs on any gene via a JSON config, and the phenotype-first analysis runs on any patient's HPO list; the chromosome-instability hazard filter is shared across MVA genes (*BUB1B*, *TRIP13*, *CEP57*). Every run writes a full exclusion log and provenance file, so results are auditable and reproducible.

Reproduce:

```bash
cd track2/pipeline
python3 repurpose.py --config config/bub1b.json --out results/
# 19 Oct 2026 run: 6 retained (5 literature + 1 for review), 551 excluded; live databases drift

cd ../phenotype
python3 phenotype.py --hpo patient_hpo.txt --gene BUB1B --related CEP57,TRIP13 \
                     --config ../pipeline/config/bub1b.json --out results/
# public demo: --hpo example_public_mva1.txt

cd ../lincs
python3 signature_reversal.py --genes BUB1B,CEP57,TRIP13 --config ../pipeline/config/bub1b.json --out results/
# expected: axis test p = 0.085; 93 of top 100 reversers switch off the alarm

# local only (never commit the output): second-diagnosis check
python3 ../tools/check_vcf_genes.py --vcf /path/to/patient.vcf.gz
```

---

## 12. Methods and AI use

- **Data used:** the hackathon case data (Track 1 genotype; patient HPO terms, processed locally) and public resources (DGIdb, Open Targets, HPO, STRING, LINCS L1000 via Enrichr, WHO ATC via ChEMBL, published literature). No subject data are reproduced in this report or the repository.
- **Code:** `track2/pipeline/repurpose.py` (gene-first screen) and `track2/phenotype/phenotype.py` (phenotype-first analysis; patient HPO terms are read from a local, git-ignored file), `track2/lincs/signature_reversal.py` (transcriptome-first route) and `track2/tools/check_vcf_genes.py` (local second-diagnosis check); configs in `track2/pipeline/config/`; every run writes a provenance file recording sources, configuration and gene set.
- **AI assistance:** Anthropic, Claude Opus 5.5 via claude.ai, Max plan; "Help improve Claude" setting **off**. Used for literature triage, drafting text and writing pipeline code; all references and outputs were checked by the team. The patient's HPO term IDs were entered into the assistant; no VCF, clinical document or other case files were uploaded, and the second-diagnosis VCF check was run locally.

## 13. Limitations

- The variant effects are predictions; protein levels in this child's cells have not been measured.
- Mouse and human BubR1 may differ in regulation; the human K668 equivalence should be confirmed.
- NR safety data come from a different disease and small, open-label studies.
- Database results depend on DGIdb and Open Targets at run time; counts can change as these update.
- LINCS signatures come from cancer cell lines and measure transcripts, not protein stability; NR and NMN are not in LINCS.
- The genetic-or-not classification uses HPO annotations, which are incomplete; "secondary" and "acquired" are hypotheses pending the clinical history.

**Future work.** Review copy number across the tubulopathy genes (especially *CLCNKB* deletions) if the clinical history does not explain the renal features; add variant-effect scores (AlphaMissense, SpliceAI) for the missense allele; test NR and the exploratory division-restoring compounds side by side in patient cells (Section 10).

---

## References

1. Suijkerbuijk SJE, et al. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res.* 2010;70(12):4891–4900.
2. Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in *BUB1B*. *Nat Genet.* 2004;36(11):1159–1161.
3. Elowe S. Bub1 and BubR1: at the interface between chromosome attachment and the spindle checkpoint. *Mol Cell Biol.* 2011;31(15):3085–3093.
4. Suijkerbuijk SJE, et al. The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. *Dev Cell.* 2012;22(6):1321–1329.
5. North BJ, et al. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO J.* 2014;33(13):1438–1453. PMID 24825348.
6. Cosentino C, Mostoslavsky R. Sirtuin to the rescue: SIRT2 extends life span of BubR1 mice. *EMBO J.* 2014;33(13):1417–1419.
7. Baker DJ, et al. BubR1 insufficiency causes early onset of aging-associated phenotypes and infertility in mice. *Nat Genet.* 2004;36(7):744–749.
8. Veenhuis SJG, et al. Nicotinamide riboside improves ataxia scores and immunoglobulin levels in ataxia telangiectasia. *Mov Disord.* 2021. doi:10.1002/mds.28788.
9. Presterud R, et al. Long-term nicotinamide riboside use improves coordination and eye movements in ataxia telangiectasia. *Mov Disord.* 2024;39(2):360–369. doi:10.1002/mds.29645. PMID 37899683.
10. Nicotinamide riboside for ataxia telangiectasia: a report of an early treated individual. *Neuropediatrics.* 2023;54(1):78–81. doi:10.1055/a-1959-9404.
11. Matsuda D, et al. BUBR1 insufficiency in mice increases their sensitivity to oxidative stress. *In Vivo.* 2016;30(6):769–776. PMID 27815460.
12. Baker DJ, et al. Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders. *Nature.* 2011;479(7372):232–236.
13. Freshour SL, et al. Integration of the Drug–Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. *Nucleic Acids Res.* 2021;49(D1):D1144–D1151.
14. Ochoa D, et al. The next-generation Open Targets Platform. *Nucleic Acids Res.* 2023;51(D1):D1353–D1359.

15. MedlinePlus Genetics. Mosaic variegated aneuploidy syndrome. U.S. National Library of Medicine. https://medlineplus.gov/genetics/condition/mosaic-variegated-aneuploidy-syndrome/
16. Gargano MA, et al. The Human Phenotype Ontology in 2024: phenotypes around the world. *Nucleic Acids Res.* 2024;52(D1):D1333–D1346.
17. Szklarczyk D, et al. The STRING database in 2023. *Nucleic Acids Res.* 2023;51(D1):D638–D646.
18. Subramanian A, et al. A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. *Cell.* 2017;171(6):1437–1452.
19. Xie Z, et al. Gene set knowledge discovery with Enrichr. *Curr Protoc.* 2021;1(3):e90.
20. Bakhoum SF, et al. Chromosomal instability drives metastasis through a cytosolic DNA response. *Nature.* 2018;553(7689):467–472.
21. WHO Collaborating Centre for Drug Statistics Methodology. ATC classification, group L01 (antineoplastic agents), accessed via ChEMBL.
22. Posey JE, et al. Resolution of disease phenotypes resulting from multilocus genomic variation. *N Engl J Med.* 2017;376(1):21–31.

*Data from the hackathon dataset are not reproduced here, per hackathon rules. This report is released under CC BY 4.0.*
