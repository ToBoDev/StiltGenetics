#!/bin/bash

NUMPROC=50
MEMORY=475
PROJ=${PWD##*/}

##rename seq files
#rename libraries to match dDocent
rename -f -e 's/STILT/STILT_/' -- ./*.f*q.gz
rename -f -e 's/\_001//' -- ./*.f*q.gz
rename -f -e 's/\_L008//' -- ./*.f*q.gz


#Trim galore: 
#trim adapter seqs, overrep seqs, qual < 20 (also dedup w/ erika's suggestion?)

if [ ! -d "trimming" ]; then mkdir trimming; fi
if [ ! -d "trimming/fqc_out" ]; then mkdir trimming/fqc_out; fi
if [ ! -d "assembly" ]; then mkdir assembly; fi
if [ ! -d "quast" ]; then mkdir quast; fi
if [ ! -d "mapping" ]; then mkdir mapping; fi
if [ ! -d "mapping/logs" ]; then mkdir mapping/logs; fi
if [ ! -d "variantcalling" ]; then mkdir variantcalling; fi
if [ ! -d "variantcalling/logs" ]; then mkdir variantcalling/logs; fi


###############################################################################
#####TRIMMING#####-------------------------------------------------------------
###############################################################################

if test -n "$(find ./raw -maxdepth 1 -name '*R2*' -print -quit)"; then
    parallel -j  $NUMPROC trim_galore --paired --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" | sort ` :::+ ` find  ./raw  -name "*_R2*.f*.gz" | sort `
    find ./trimming -name "*val_1*.gz" | sort | uniq > fwds
    find ./trimming -name "*val_2*.gz" | sort | uniq > revs
else
   parallel  -j $NUMPROC trim_galore  --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" `
fi 

###############################################################################
#####SPADES ASSEMBLY#####------------------------------------------------------
###############################################################################

if test -n "$(find ./raw -maxdepth 1 -name '*R2*' -print -quit)"; then
    echo reverse reads found. proceeding with paired end assembly
    #setspades.R
    echo "if (\"here\" %in% row.names(installed.packages())){" > spades_yaml.R
    echo "library(here)" >> spades_yaml.R
    echo "}  else {" >> spades_yaml.R
    echo "install.packages(\"here\")" >> spades_yaml.R
    echo "library(here)}" >> spades_yaml.R
    echo "workingdir <- paste0(here(),\"/\")" >> spades_yaml.R
    echo "fwds <- read.table(paste0(workingdir,\"/fwds\"))" >> spades_yaml.R
    echo "revs <- read.table(paste0(workingdir,\"/revs\"))" >> spades_yaml.R
    echo "setwd(workingdir)" >> spades_yaml.R
    echo "write(paste0('[ \n',  '   { \n', '     orientation: \"fr\", \n', '     type: \"paired-end\", \n', '     right reads: ['), file = \"libraries.yaml\", append = F)" >> spades_yaml.R
    echo "for(i in 1:nrow(fwds)){write(paste0('       \"',fwds[i,1],'\",'), file = \"libraries.yaml\", append = T)}" >> spades_yaml.R
    echo "write(paste0('       ], \n', '       left reads: ['), file = \"libraries.yaml\", append = T)" >> spades_yaml.R
    echo "for(i in 1:nrow(revs)){write(paste0('       \"',revs[i,1],'\",'), file = \"libraries.yaml\", append = T)}" >> spades_yaml.R                   
    echo "write(paste0('       ] \n', '   } \n', ']'), file = \"libraries.yaml\", append = T)" >> spades_yaml.R

    chmod 755 spades_yaml.R
    Rscript spades_yaml.R

    #paired-end assembly
    #spades.py --dataset libraries.yaml -k 71,81,91,99,121,127 -t $NUMPROC -m $MEMORY -o ./assembly/${PROJ}_usr_kmer/
    spades.py --dataset libraries.yaml -t $NUMPROC -m $MEMORY -o ./assembly/${PROJ}_default_kmer/
else
    echo no reverse reads found. proceeding with single end assebmly
    #single end assembly
    pools=(` find ./trimming -name "*merged*.gz" | sort | uniq | cat `) 
    libnum=( `seq 1 "${#pools[@]}"` )
    libnum=("${libnum[@]/#/--s}")
    unset reads
    for (( i=0; i<${#libnum[*]}; ++i)); do reads+=( ${libnum[$i]} ${pools[$i]} ); done

    #spades.py "${reads[@]}" -k 71,81,91,121,127 -t $NUMPROC -m $MEMORY -o ./assembly/${PROJ}_usr_kmer
    spades.py "${reads[@]}" -t $NUMPROC -m $MEMORY -o ./assembly/${PROJ}_default_kmer
fi

#QUAST
parallel   -j $NUMPROC cp ./assembly/{}/contigs.fasta ./quast/{}_contigs.fasta ::: \
` find ./assembly -maxdepth 1 -type d | sort | uniq | cut -d "/" -f3 `

all_default=(` find ./quast/*contigs.fasta  -type f `)

quast.py --large -e -k -t $NUMPROC -m $MEMORY "${all_default[@]}" -o ./quast/default_spades/  
#########STOP######-----------------------------------------------------------
#before proceeding: evaluate your quast results. Choose desired de novo, copy to ./assembly with a name to include chosen.fasta
cp ./assembly/${PROJ}_default_kmer/contigs.fasta ./assembly/${PROJ}_chosen.fasta

###############################################################################
#####MAPPING#####--------------------------------------------------------------
###############################################################################

#THECHOSEN=./assembly/<nameofchosenreference>
THECHOSEN=(` find ./assembly/ -maxdepth 1 -name "*chosen*" -type f `)
LIBS=(` find  ./trimming  -name "*R*.fq.gz" | cut -d "/" -f3 `)
POPS=(` find  ./trimming  -name "*.fq.gz" | cut -d "_" -f1-2 | cut -d "/" -f3 | cut -d "-" -f2 | sort | uniq `)
optA=1
optB=4
optO=6
SPLITS=2

#filter and largest squences for blast
#from:http://itrylinux.com/use-awk-to-filter-fasta-file-by-minimum-sequence-length/

for i in `find ./quast/ -maxdepth 1 -name "*fasta" -type f` ; do
    sed ':a;N;/^>/M!s/\n//;ta;P;D' ${i} | awk '/^>/ { getline seq } length(seq) >50000 { print $0 "\n" seq }' - > ${i/_kmer_contigs.fasta/}_large.fa
done

cp $THECHOSEN ./mapping/${PROJ}_reference.fasta

parallel -j $NUMPROC bwa index {} ::: \
` find ./mapping ./trimming -name "*.fasta" -o -name "*.fq.gz" `
#` find ./trimming ./mapping -name "*.f*.gz" -o -name "*.fasta" `

if test -n "$(find ./raw -maxdepth 1 -name '*R2*' -print -quit)"; then
    #PAIRED
    parallel -j "$(( ${#LIBS[@]} / ($SPLITS * 2) ))" \
        "bwa mem ./mapping/${PROJ}_reference.fasta {1} {2} "\
        "-L 20,5 -t $(( $NUMPROC / $((${#LIBS[@]} / ($SPLITS * 2))) )) "\
        "-M -T 35 -A $optA -B $optB -O $optO -R "\
        " \"@RG\\tID:{3}\\tSM:{3}\\tPL:Illumina\" 2> ./mapping/logs/bwa.{3}.log "\
        " | mawk '\$6 !~/[2-9].[SH]/ && \$6 !~ /[1-9][0-9].[SH]/' "\
        " | samtools view -@ $(( $NUMPROC / $SPLITS )) -q 1 -SbT ./mapping/${PROJ}_reference.fasta - > ./mapping/{3}.bam 2> ./mapping/logs/{3}.bwa.log" ::: \
         ` find  ./trimming  -name "*val_1*.f*.gz" | sort ` :::+ \
         ` find  ./trimming  -name "*val_2*.f*.gz" | sort ` :::+ \
         ` find  ./trimming  -name "*val_1*.f*.gz" | sort | cut -d "/" -f3 | cut -d "_" -f1-2 `

    parallel -j "$(( ${#LIBS[@]} / ($SPLITS * 4) ))" samtools sort -@ "$(( $NUMPROC / $((${#LIBS[@]} / ($SPLITS * 4))) ))" -O BAM -o {.}.sort.bam {.}.bam ::: \
    ` find ./mapping -name "*.bam" -not -name "*sort*" | sort `

else
    #Single end
    parallel -j "$(( ${#LIBS[@]} / ($SPLITS * 2) ))" \
        "bwa mem ./mapping/${PROJ}_reference.fasta {1} "\
        "-L 20,5 -t $(( $NUMPROC / $((${#LIBS[@]} / ($SPLITS * 2))) )) "\
        "-M -T 35 -A $optA -B $optB -O $optO -R "\
        " \"@RG\\tID:{2}\\tSM:{2}\\tPL:Illumina\" 2> ./mapping/logs/bwa.{2}.log "\
        " | mawk '\$6 !~/[2-9].[SH]/ && \$6 !~ /[1-9][0-9].[SH]/' "\
        " | samtools view -@ $(( $NUMPROC / $SPLITS )) -q 1 -SbT ./mapping/${PROJ}_reference.fasta - > ./mapping/{2}.bam 2> ./mapping/logs/{2}.bwa.log" ::: \
         ` find  ./trimming  -name "*R1*.f*.gz" | sort ` :::+ \
         ` find  ./trimming  -name "*R1*.f*.gz" | sort | cut -d "/" -f3 | cut -d "_" -f1-2 `    

    parallel -j "$(( ${#LIBS[@]} / ($SPLITS * 2) ))" samtools sort -@ "$(( $NUMPROC / $((${#LIBS[@]} / ($SPLITS * 4))) ))" -O BAM -o {.}.sort.bam {.}.bam ::: \
    ` find ./mapping -name "*.bam" -not -name "*sort*" | sort `
fi


###############################################################################
#####VARIANT CALLING#####------------------------------------------------------
###############################################################################
parallel -j $NUMPROC samtools index {} ::: ` find ./mapping -name "*.sort.bam" `
find ./mapping -name "*.sort.bam" > ./variantcalling/bams.list
samtools merge -@ $NUMPROC -b ./variantcalling/bams.list -f ./variantcalling/merged_pools.bam
samtools index ./variantcalling/merged_pools.bam
samtools faidx ./mapping/${PROJ}_reference.fasta

POPS=(` find  ./trimming  -name "*val_1*.f*.gz" | sort | cut -d "/" -f3 | cut -d "_" -f1-1 | sort | uniq `)

bamToBed -i ./variantcalling/merged_pools.bam | bedtools merge -i - > ./variantcalling/mapped.bed

split -n l/$NUMPROC -d --additional-suffix=.bed ./variantcalling/mapped.bed ./variantcalling/mapped.

wait

printf "%s\n" "${POPS[@]}" > poplist
paste poplist poplist > popmap

if [ ! -f "./variantcalling/reference.fasta" ]
    then cp ./mapping/${PROJ}_reference.fasta ./variantcalling/reference.fasta
    fi

if [ ! -d "./variantcalling/raw" ]; then
     mkdir ./variantcalling/raw
fi

if [ ! -d "./variantcalling/filtered" ]; then
     mkdir ./variantcalling/filtered
fi

   call_genos(){
    if [ ! -f "./variantcalling/split.$1.bam" ]
    then samtools view -b -1 -L ./variantcalling/mapped.$1.bed -o ./variantcalling/split.$1.bam ./variantcalling/merged_pools.bam
    fi
    samtools index ./variantcalling/split.$1.bam; 

    freebayes -b ./variantcalling/split.$1.bam -t ./variantcalling/mapped.$1.bed \
    -f ./variantcalling/reference.fasta \
        --min-alternate-count 2 \
        --min-mapping-quality 20 \
        --min-base-quality 20 \
        --min-repeat-entropy 0 -V \
        --min-coverage 30 \
        --min-alternate-fraction 0.05 \
        --genotype-qualities \
        --populations popmap \
         2> ./variantcalling/logs/fb.$1.error.log | grep -v "#contig" > ./variantcalling/raw/raw.$1.vcf
    }

    export -f call_genos

ulimit -s unlimited
ls ./variantcalling/mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | cut -d "/" -f3 | shuf | parallel --env call_genos -j $NUMPROC --no-notice call_genos {} 

rename -f -e 's/\d+/sprintf("%02d",$&)/e' -- ./variantcalling/raw/*.vcf
vcflib vcfcombine ./variantcalling/raw/raw.*.vcf > ./variantcalling/TotalRawSNPs.vcf

################################################ Maybe dont need this section as dont think we used this file? #####################################################

#if output is polyploid for non GT calls (ie ././.)
grep -v "^#" ./variantcalling/TotalRawSNPs.vcf | cut -f 10- | cut -d ':' -f 1 | sort | uniq
sed 's/\.\/\.\/\./\.\/\./' ./variantcalling/TotalRawSNPs.vcf > ./variantcalling/TotalRawSNPs_fix.vcf
rm ./variantcalling/TotalRawSNPs.vcf; mv ./variantcalling/TotalRawSNPs_fix.vcf ./variantcalling/${PROJ}_TotalRawSNPs.vcf

###############################################################################
#VCF_FILTERING#####-----------------------------------------------------------
###############################################################################

bgzip < TotalRawSNPs.vcf > TotalRawSNPs.vcf.gz && tabix TotalRawSNPs.vcf.gz
# Initial look at individual/sample missingness
vcftools --gzvcf TotalRawSNPs.vcf.gz --missing-indv --out initial_ind_missingness
# sample 32 removed (~99% missingness), max alelles 2, max depth 2000, snps only, at least 1 non-ref allele
bcftools view -s "^manu_32" --max-alleles 2 -i 'MEAN(FORMAT/DP)<=2000' -v snps -c 1 TotalRawSNPs.vcf.gz -Oz -o spades_denovo_light_filters.vcf.gz

vcftools --gzvcf spades_denovo_light_filters.vcf.gz --max-missing 0.5 --mac 2 --minmeanDP 3 --maxmeanDP 500 --recode --stdout filter_spades_denovo_light_filters.vcf
vcffilter -s -f "MQM > 30 & MQMR > 30" -f "MQM / MQMR > 0.75 & MQM / MQMR < 1.25" -f "QUAL / DP > 0.25" /
-f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25" -f "NS > 44.5" /
-f "DP > 30" -f "QUAL > 20" -f "TYPE = snp" -f "AO > 2" filter_spades_denovo_light_filters.vcf > final_filter_spades_denovo_light_filters.vcf

# thinning
vcftools --vcf final_filter_spades_denovo_light_filters.vcf --thin 1000 --recode --stdout final_filter_spades_denovo_light_filters.vcf

# finetuning --max-missing and --min-meanDP filters
vcfin="final_filter_spades_denovo_light_filters"
vcfout="spades_20perc_miss"
missingness=0.8 #0-1 with 0 allowing for all missing data at a site and 1 allowing for no missing data
vcftools --vcf ${vcfin}.vcf --max-missing $missingness --recode --recode-INFO-all --stdout > ${vcfout}.vcf

vcfin="spades_20perc_miss"
vcfout="spades_20perc_5mmd"
vcftools --vcf ${vcfin}.vcf --min-meanDP 5 --recode --recode-INFO-all --stdout >  ${vcfout}.vcf

# Remove invariant sites and alleles with len > 1
bgzip < spades_20perc_5mmd.vcf > spades_20perc_5mmd.vcf.gz && tabix spades_20perc_5mmd.vcf.gz
# require at least 1 homozygous site (remove invariant fully heterozygous sites)
bcftools view -g hom spades_20perc_5mmd.vcf.gz -Oz -o spades_20perc_5mmd_noinv.vcf.gz
# remove alleles with length not equal to 1 (no LEN in vcf header)
bcftools filter -e 'strlen(REF)!=1 || strlen(ALT)!=1' spades_20perc_5mmd_noinv.vcf.gz -Oz -o spades_20perc_5mmd_noinv_LEN1.vcf.gz


############### Filter for max heterozygosity based on minor allele frequency################################
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

# checked plot in R, set threshold at 0.4
Rscript << 'EOF'
MAF<-read.table("MAF.txt", header=T)
het<-read.table("heterozygosity_fhet_calc.txt", header=T)
MAF_het<-cbind(MAF,het$FHET_OBS)
# Get list of loci (by CHROM, POS) with observed heterozygosity < 0.4
MAF_het_cut<-MAF_het[MAF_het$`het$FHET_OBS`<=0.4,]
write.table(MAF_het_cut[,c(1,2)],"list_SNPs_het0.4.txt",quote=F, row.names=F, sep="\t") 
EOF

# retain only the loci with observed heterozygosity < 0.4
vcftools --gzvcf spades_20perc_5mmd_noinv_LEN1.vcf.gz --positions list_SNPs_het0.4.txt --recode --recode-INFO-all --out spades_20perc_5mmd_noinv_LEN1_fhet0.4

###########################################################################################################


# Check for individual missingness
vcftools --vcf spades_20perc_5mmd_noinv_LEN1_fhet0.4.recode.vcf --missing-indv --out per_ind_missingness
# Filter at 30%
cat per_ind_missingness.imiss | awk '$5 > 0.3' | awk '{print $1}' > samples_fmiss0.3_rm.txt
vcftools --vcf spades_20perc_5mmd_noinv_LEN1_fhet0.4.recode.vcf --remove samples_fmiss0.3_rm.txt --recode --recode-INFO-all --out spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3

# Repeat loci filters impacted by removing individuals
bcftools view -v snps -c 1 -g hom spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3.recode.vcf -Oz -o spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3_bcftools.vcf.gz
vcftools --gzvcf spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3_bcftools.vcf.gz --max-missing 0.8 --mac 2 --min-meanDP 5 --max-meanDP 500 --out spades_20perc_5mmd_noinv_LEN1_fhet0.4_imiss0.3_bcftools_vcftools




################ #plink reformatting ########### Maybe remove this section as dont think we used these files? ###########################################

vcf_in=spades_20perc_5mmd
prefix=spades
plink2 --vcf ${vcf_in}.vcf --allow-extra-chr --out $prefix
plink --vcf ${vcf_in}.vcf --allow-extra-chr --out $prefix

vcftools --vcf spades_20perc_5mmd.vcf --out ${prefix}_for_plink --plink

plink --file spades_for_plink --recode12 --out ${prefix}_plink_recode
plink --file spades --recodeA --out ${prefix}_plink_recode


vcf_in=ref2_20perc_5mmd
prefix=ref2
vcftools --vcf ${vcf_in}.vcf --out ${prefix}_for_plink --plink
plink --vcf ${vcf_in}.vcf --allow-extra-chr --out $prefix
plink --file ${prefix}_for_plink --recode 12 --out ${prefix}_plink_recode
