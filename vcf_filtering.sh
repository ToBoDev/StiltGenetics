### mapped to spades denovo genome

### vcf calling with freebayes using GNU parallel, including some basic freebayes filters 

### Light filters - Maddie

# sample 32 removed (~99% missingness), max alelles 2, max depth 2000, snps only, at least 1 non-ref genome
bgzip < TotalRawSNPs.vcf > TotalRawSNPs.vcf.gz && tabix TotalRawSNPs.vcf.gz
bcftools view -s "^manu_32" --max-alleles 2 -i 'MEAN(FORMAT/DP)<=2000' -v snps -c 1 TotalRawSNPs.vcf.gz -Oz -o TotalRawSNPs_91ind_snps_maxmeanDP2000.vcf.gz

### More filters with assesspool - Evan (added to supp table as below)

mv TotalRawSNPs_91ind_snps_maxmeanDP2000.vcf.gz ./spades_denovo_light_filters.vcf.gz

vcftools --max-missing 0.5
--mac 2
--minmeanDP 3
--maxmeanDP 500

vcflib MQM > 30 MQMR > 30 
MQM / MQMR > 0.75 & MQM / MQMR < 1.25
QUAL / DP > 0.25
PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25
NS > 44.5
DP > 30
LEN < 2   #will delete this one - didnt include it in supp table, as LEN was not in vcf header so did not work
QUAL > 20
TYPE = snp
AO > 2

### Finetuning filters - Evan

vcftools --thin 1000
--max-missing 0.8
--minmeanDP 5

### Additional filters - Maddie

bcftools view -g hom
strlenREF ALT

