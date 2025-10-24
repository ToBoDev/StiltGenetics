Publicly available code for Price et al 2026 "Breeding behaviors and kinship in the Endangered Hawaiian Stilt (Ae‘o; Himantopus mexicanus knudseni)". 

Raw reads were processed using Trim Galore! v0.6.4 (Krueger, 2012) for trimming and adapter removal, BWA v0.7.17-r1188 (bwa-mem algorithm, Li, 2013) for read mapping, and FreeBayes v1.3.8 (Garrison and Marth, 2012) for variant calling. Reads were mapped to a <i>de novo</i> genome assembled with SPAdes v3.13.2 (Bankevich et al., 2012) (91.64% properly paired, 96% >= MAPQ 30).

A single tree file summarized from BEAST output using treeannotator, was visualized as a circular tree with the {ggtree} package (Yu et al., 2017) in R (R Core Team, 2021).
