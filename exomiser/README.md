# Exomiser genome-wide confirmation (report §4.1)

Exomiser CLI 15.1.0, data release 2512 (GRCh38), `exome` preset, six HPO terms, no gene list.

| Rank | Gene | MOI | p-value |
|---|---|---|---|
| 1 | BUB1B | AR | 0.0006 |
| 2 | BUB1B | AD | 0.0048 |
| 3 | LZTR1 | AD | 0.031 |

BUB1B (AR) is driven by chr15:40209701 T>G (p.Leu737Ter, variant score 1.000) and chr15:40220612 T>G (p.Asn1002Lys, 0.923), the same pair as the primary analysis.

Filter funnel: ~4.96 M calls → 23,930 pass variant effect → 1,010 pass frequency → 588 pass inheritance, in 363 genes.

Files: `run_exomiser.sh` (command), `application.properties.snippet` (config changes), `proband_phenopacket.template.yml` (format only). The real phenotype file and all Exomiser outputs contain subject data and are not included.
