##### mapped to spades denovo genome - Evan #####

##### vcf calling with freebayes using GNU parallel, including some basic freebayes filters - Evan #####

##### Light filters - Maddie #####

# sample 32 removed (~99% missingness), max alelles 2, max depth 2000, snps only, at least 1 non-ref allele
bgzip < TotalRawSNPs.vcf > TotalRawSNPs.vcf.gz && tabix TotalRawSNPs.vcf.gz
bcftools view -s "^manu_32" --max-alleles 2 -i 'MEAN(FORMAT/DP)<=2000' -v snps -c 1 TotalRawSNPs.vcf.gz -Oz -o TotalRawSNPs_91ind_snps_maxmeanDP2000.vcf.gz

##### More filters with assesspool - Evan (added to supp table as below) #####

# Files renamed to spades_denovo_light_filters.vcf.gz

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

##### Finetuning filters - Evan #####

# thin to 1000, max missing 0.8, minmeanDP 5

vcftools --gzvcf filtered_spades_denovo_light_filters.vcf.gz --thin 1000 --recode --stdout final_filter_spades_denovo_light_filters.vcf

vcfin="final_filter_spades_denovo_light_filters"
vcfout="spades_20perc_miss"

vcfin="spades_20perc_miss"
vcfout="spades_20perc_5mmd"

missingness = 0.8 #0-1 with 0 allowing for all missing data at a site and 1 allowing for no missing data

vcftools --vcf ${vcfin}.vcf --max-missing $missingness --recode --stdout > ${vcfout}.vcf
vcftools --vcf ${vcfin}.vcf --min-meanDP 5 --recode --stdout >  ${vcfout}.vcf

##### Additional filters - Maddie #####

# at least 1 homozygous site (remove invariant fully heterozygous sites), remove alleles with length not equal to 1
bgzip < spades_20perc_5mmd.vcf > spades_20perc_5mmd.vcf.gz && tabix spades_20perc_5mmd.vcf.gz
bcftools view -g hom spades_20perc_5mmd.vcf.gz -Oz -o spades_20perc_5mmd_noinv.vcf.gz
bcftools view spades_20perc_5mmd_noinv.vcf.gz | bcftools filter -e 'strlen(REF)!=1 || strlen(ALT)!=1' -Oz -o spades_20perc_5mmd_noinv_LEN1.vcf.gz

### Filter for max heterozygosity based on minor allele frequency

# calculate allele frequencies
vcftools --gzvcf spades_20perc_5mmd_noinv_LEN1.vcf.gz --freq --out MAF
# determine the minor allele frequency and print to new file
awk 'BEGIN {OFS="\t";
    print "CHROM", "POS", "MAF"}
NR>1 { split($5,a,":")
split($6,b,":") 
if (a[2] < b[2]) {
maf = a[2]
} else {
maf = b[2] }
print $1, $2, maf }' MAF.frq > MAF.txt

# print ref homozygotes, alt homozygotes and heterozygotes at each loci
vcftools --gzvcf spades_20perc_5mmd_noinv_LEN1.vcf.gz --hardy --out heterozygosity
# calculate frequency of observed heterozygosity per locus, from total ref homozygotes, heterozygotes, and alt homozygotes
awk 'BEGIN {OFS="\t";
    print "CHROM", "POS", "FHET_OBS"}
NR>1 { split($3,o,"/")
fhet_o = o[2]/(o[1]+o[2]+o[3])
print $1, $2, fhet_o }' heterozygosity.hwe > heterozygosity_fhet_calc.txt

# plot fin plot in R
MAF<-read.table("MAF.txt", header=T)
het<-read.table("heterozygosity_fhet_calc.txt", header=T)
plot(MAF$MAF, het$FHET_OBS)

# Get list of loci (by CHROM, POS) with observed heterozygosity < 0.45 in R
MAF_het<-cbind(MAF,het$FHET_OBS)
MAF_het_cut<-MAF_het[MAF_het$`het$FHET_OBS`<=0.4,]
kept_SNPs<-MAF_het_cut[,c(1,2)]
write.table(kept_SNPs,"list_SNPs_het0.4.txt",quote=F, row.names=F, sep="\t") 

# use vcftools to retain only these loci with observed heterozygosity < 0.45
vcftools --gzvcf spades_20perc_5mmd_noinv_LEN1.vcf.gz --positions list_SNPs_het0.4.txt --recode --recode-INFO-all --out spades_20perc_5mmd_noinv_LEN1_fhet0.4

###

# Filter for individual missingness
vcftools --vcf spades_20perc_5mmd_noinv_LEN1_fhet0.4.recode.vcf --missing-indv --out per_ind_missingness
cat per_ind_missingness.imiss | awk '$5 > 0.3' | awk '{print $1}' > samples_fmiss0.3_rm.txt
vcftools --vcf spades_20perc_5mmd_noinv_LEN1_fhet0.4.recode.vcf --remove samples_fmiss0.3_rm.txt --recode --recode-INFO-all --out spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3

# Repeat bcftools and vcftools loci filters that would have been impacted by removing individuals
bcftools view -v snps -c 1 -g hom spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3.recode.vcf -Oz -o spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3_bcftools.vcf.gz
vcftools --gzvcf spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3_bcftools.vcf.gz --max-missing 0.8 --mac 2 --min-meanDP 5 --max-meanDP 500 --out spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3_bcftools_vcftools

