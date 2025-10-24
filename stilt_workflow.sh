#!/bin/bash

NUMPROC=50
MEMORY=475
PROJ=${PWD##*/}

##rename seq files
#for f in * ; do mv -- "$f" "${PROJ}_$f" ; done
#echo "for f in * ; do mv -- "$f" "${PROJ}_$f" ; done"
#rename libraries to match dDocent
rename -f -e 's/STILT/STILT_/' -- ./*.f*q.gz
#rename -f -e 's/_R1_/.F_/' -- ./*.f*q.gz
#rename -f -e 's/_R2_/.R_/' -- ./*.f*q.gz
#rename -f -e 's/fastq/fq/' -- ./*.f*q.gz
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


      
#####TRIMMING#####-----------------------------------------------------------
if test -n "$(find ./raw -maxdepth 1 -name '*R2*' -print -quit)"; then
    parallel -j  $NUMPROC trim_galore --paired --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" | sort ` :::+ ` find  ./raw  -name "*_R2*.f*.gz" | sort `
    find ./trimming -name "*val_1*.gz" | sort | uniq > fwds
    find ./trimming -name "*val_2*.gz" | sort | uniq > revs
else
   parallel  -j $NUMPROC trim_galore  --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" `
fi 


#####SPADES ASSEMBLY#####
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

#####QUAST#####-----------------------------------------------------------
# #per KMER
# for genome in `find ./assembly -maxdepth 1 -name *$PROJ* -type d | cut -d "/" -f3`; do
#     KMERS=(` find ./assembly/$genome/ -maxdepth 1 -name "K*" | cut -d "/" -f4 | sort | uniq `)
#     cp ./assembly/${genome}/contigs.fasta ./quast/"${genome}_contigs.fasta"
#          for KMER in `find ./assembly/$genome/ -maxdepth 1 -name "K*" | cut -d "/" -f4 | sort | uniq`; do
#              cp ./assembly/${genome}/${KMER}/final_contigs.fasta ./quast/${genome}_${KMER}.fasta
#         done
#     by_kmer=("${KMERS[@]/#/./quast/${genome}_}")
#     by_kmer_files=("${by_kmer[@]/%/.fasta}")
#     quast.py -e -t $NUMPROC -m $MEMORY ./quast/${genome}_contigs.fasta "${by_kmer_files[@]}" -o ./quast/${genome}_results/
# done

# #all default KMER selection
# all_user=(` find ./quast/*usr*contigs.fasta  -type f `)
# quast.py -e -t $NUMPROC -m $MEMORY  "${all_user[@]}" -o ./quast/usr_spades/

# #spades default settings
# #all default KMER selection

# all_default=(` find ./quast/*default*contigs.fasta  -type f `)
# quast.py -e -t $NUMPROC -m $MEMORY "${all_default[@]}" -o ./quast/default_spades/ 


parallel   -j $NUMPROC cp ./assembly/{}/contigs.fasta ./quast/{}_contigs.fasta ::: \
` find ./assembly -maxdepth 1 -type d | sort | uniq | cut -d "/" -f3 `

all_default=(` find ./quast/*contigs.fasta  -type f `)

quast.py --large -e -k -t $NUMPROC -m $MEMORY "${all_default[@]}" -o ./quast/default_spades/  
#########STOP######-----------------------------------------------------------
#before proceeding: evaluate your quast results. Choose desired de novo, copy to ./assembly with a name to include chosen.fasta
cp ./assembly/${PROJ}_default_kmer/contigs.fasta ./assembly/${PROJ}_chosen.fasta

#####MAPPING#####-----------------------------------------------------------
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

#sed ':a;N;/^>/M!s/\n//;ta;P;D'  {THECHOSEN} > {THECHOSEN}_sed_corr.fasta
#awk '/^>/ { getline seq } length(seq) >20000 { print $0 "\n" seq }' {THECHOSEN}_sed_corr.fasta > {THECHOSEN}_large_contigs.fasta

#grep -A 1 NODE_2_le sed_corr.fasta > NODE_2.fasta
#grep -A 1 -f desired_nodes.txt sed_corr.fasta > desired_nodes.fasta

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
#

######BamQC###################-----------------------------------------
if [ ! -d "mapping/bamqc_out" ]; then mkdir mapping/bamqc_out; fi

/10tb_leviathan/evan/bin/BamQC/bin/bamqc \
    -g ./mapping/ \
    -o ./mapping/bamqc_out \
    -t 50 \
    ./mapping/

if [ ! -d "ref_mapping/bamqc_out" ]; then mkdir ref_mapping/bamqc_out; fi
    
/10tb_leviathan/evan/bin/BamQC/bin/bamqc \
    -a ./ref_mapping/manu2_ncbi_eference.fasta \
    -o ./ref_mapping/bamqc_out \
    -t 50







NUMPROC=40
MEMORY=500
PROJ=${PWD##*/}
#####VARIANT CALLING#####-----------------------------------------------------------
parallel -j $NUMPROC samtools index {} ::: ` find ./mapping -name "*.sort.bam" `
find ./mapping -name "*.sort.bam" > ./variantcalling/bams.list
samtools merge -@ $NUMPROC -b ./variantcalling/bams.list -f ./variantcalling/merged_pools.bam
samtools index ./variantcalling/merged_pools.bam
samtools faidx ./mapping/${PROJ}_reference.fasta

###############
POPS=(` find  ./trimming  -name "*val_1*.f*.gz" | sort | cut -d "/" -f3 | cut -d "_" -f1-1 | sort | uniq `)

bamToBed -i ./variantcalling/merged_pools.bam | bedtools merge -i - > ./variantcalling/mapped.bed

split -n l/$NUMPROC -d --additional-suffix=.bed ./variantcalling/mapped.bed ./variantcalling/mapped.


# bedtools coverage -b ./variantcalling/merged_pools.bam -a ./variantcalling/mapped.bed -counts -sorted > ./variantcalling/cov.stats

# DP=$(mawk '{print $4}' ./variantcalling/cov.stats | sort -rn | perl -e '$d=.001;@l=<>;print $l[int($d*@l)]')
# CC=$( mawk -v x=$DP '$4 < x' ./variantcalling/cov.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$(( $NUMPROC * 2 ))'"}')

# find ./trimming/ -name "*.fq.gz" | cat | parallel -j $NUMPROC "gunzip -c {} | head -2 | tail -1 >> ./variantcalling/lengths.txt"
# MaxLen=$(mawk '{ print length() | "sort -rn" }' ./variantcalling/lengths.txt | head -1)
# MaxLen2=$(( $MaxLen / 2 ))
# ML1=$(( $MaxLen2 + 1 ))
# MaxCov=$(mawk '{print $4}' ./variantcalling/cov.stats | sort -rn | head -1)
# MaxCov2=$(( $MaxCov * 8 / 10 ))

# mawk '$4 > $MaxCov2' ./variantcalling/cov.stats > ./variantcalling/cov.high.stats
# mawk '$4 <= $MaxCov2' ./variantcalling/cov.stats > ./variantcalling/cov.low.stats

# #split high coverage intervals into smaller, 1/2 read-length sized intervals

# bedtools makewindows -b ./variantcalling/cov.high.stats -w $MaxLen2 -s $ML1 > ./variantcalling/temp1.bed
# bedtools intersect -a ./variantcalling/cov.stats -b ./variantcalling/temp1.bed > ./variantcalling/temp.cov.stats

# cat ./variantcalling/temp.cov.stats ./variantcalling/cov.low.stats > ./variantcalling/cov.split.stats
# rm ./variantcalling/cov.high.stats ./variantcalling/cov.low.stats ./variantcalling/temp1.bed ./variantcalling/temp.cov.stats

# TT=$(( $MaxLen2 * 1000000 ))

# mawk -v x=$DP '$4 < x' ./variantcalling/cov.split.stats |sort -V -k1,1 -k2,2 | mawk -v cutoff=$CC -v tt=$TT 'BEGIN{i=1}
# {len=$3-$2;lc=len*$4;cov = cov + lc
# if (NR == 1 && lc > tt) {x="./variantcalling/mapped."i".bed";print $1"\t"$2"\t"$3 > x; i=i+1; e=1}
# else if ( cov < cutoff && lc < tt) {x="./variantcalling/mapped."i".bed";print $1"\t"$2"\t"$3 > x; e=0}
# else if (lc > tt && e > 0 ) {x="./variantcalling/mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=0;i=i+1; e=1}
# else if (lc > tt && e < 1 ) {i=i+1; x="./variantcalling/mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=0;i=i+1;e=1}
# else if (cov > cutoff && lc < tt ) {i=i+1; x="./variantcalling/mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=lc;e=0}
# }'


wait
####################
POPS=(` find  ./trimming  -name "*val_1*.f*.gz" | sort | cut -d "/" -f3 | cut -d "_" -f1-2 | sort | uniq `)

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

# VT1="--max-missing 0.9"
# VT2="--mac 2"
# VT3="--minQ 20"
# VT4="--minDP 10"
# #VT5="--min-meanDP 5"
# #VT6="--max-meanDP 500" 

# VL1=("MQM > 39 & MQMR > 39")
# VL2=("MQM / MQMR > 0.75 & MQM / MQMR < 1.25")
# VL3=("RPR > 0 & RPL > 0")
# VL4=("QUAL / DP > 0.25")
# VL5=("PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05")
# VL6=("NS > 1")
# VL7=("LEN < 11")

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

    #vcflib vcffilter -f "QUAL > 10" raw.$1.vcf > qualfil.$1.vcf
    #add in all filtering here while the vcf is already split up?!
    /10tb_leviathan/evan/assesspool_dev/scripts/vcflib/bin/vcffilter -f \
       "QUAL > 20 & \
        DP > 30 & \
        TYPE = snp & \
        AO > 2 & \
        QUAL / DP > 0.25 & \
        PAIRED > 0.5 & \
        PAIREDR > 0.5 & \
        PAIREDR / PAIRED < 1.75 & \
        PAIREDR / PAIRED > 0.25 & \
        MQM > 39 & MQMR > 39 & \
        MQM / MQMR > 0.75 & MQM / MQMR < 1.25 & \
        NS > 2 & \
        LEN < 11" ./variantcalling/raw/raw.$1.vcf | \
    vcftools --vcf - \
        --max-missing 0.5 \
        --max-meanDP 500 \
        --min-meanDP 5 \
        --hwe 0.01 \
        --recode --recode-INFO-all \
        --stdout > ./variantcalling/filtered/filtered.$1.vcf
    #grep -v "#contig" ./variantcalling/raw.$1.vcf > ./variantcalling/raw/raw_fixed.$1.vcf
    #grep -v "#contig" ./variantcalling/filtered.$1.vcf > ./variantcalling/filtered/filtered_fixed.$1.vcf
    }

    export -f call_genos

ulimit -s unlimited
ls ./variantcalling/mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | cut -d "/" -f3 | shuf | parallel --env call_genos -j $NUMPROC --no-notice call_genos {} 
   


#mv ./mapping/reference.fasta ./mapping/${PROJ}_reference.fasta
#mv ./variantcalling/raw*.vcf ./variantcalling/raw/
#mv ./variantcalling/filtered*.vcf ./variantcalling/filtered/

rename -f -e 's/\d+/sprintf("%02d",$&)/e' -- ./variantcalling/raw/*.vcf
vcflib vcfcombine ./variantcalling/raw/raw.*.vcf > ./variantcalling/TotalRawSNPs.vcf

rename -f -e 's/\d+/sprintf("%02d",$&)/e' -- ./variantcalling/filtered/*.vcf
vcflib vcfcombine ./variantcalling/filtered/filtered.*.vcf > ./variantcalling/filtered_SNPs.vcf

######

#if output is polyploid for non GT calls (ie ././.)
grep -v "^#" ./variantcalling/TotalRawSNPs.vcf | cut -f 10- | cut -d ':' -f 1 | sort | uniq
sed 's/\.\/\.\/\./\.\/\./' ./variantcalling/TotalRawSNPs.vcf > ./variantcalling/TotalRawSNPs_fix.vcf
rm ./variantcalling/TotalRawSNPs.vcf; mv ./variantcalling/TotalRawSNPs_fix.vcf ./variantcalling/${PROJ}_TotalRawSNPs.vcf

#if output is polyploid for non GT calls (ie ././.)
grep -v "^#" ./variantcalling/filtered_SNPs.vcf | cut -f 10- | cut -d ':' -f 1 | sort | uniq
sed 's/\.\/\.\/\./\.\/\./' ./variantcalling/filtered_SNPs.vcf > ./variantcalling/filtered_SNPs_fix.vcf
rm ./variantcalling/filtered_SNPs.vcf; mv ./variantcalling/filtered_SNPs_fix.vcf ./variantcalling/${PROJ}_filtered_SNPs.vcf

rm ./variantcalling/split.*.bam*; rm ./variantcalling/mapped.*.bed
#vcftools --vcf variantcalling/${PROJ}_TotalRawSNPs.vcf --minDP 5 --thin 100 --recode --recode-INFO-all --out ${PROJ}_filteredSNPs_thin
#vcftools --vcf variantcalling/${PROJ}_TotalRawSNPs.vcf --minDP 5 --recode --recode-INFO-all --out ${PROJ}_filteredSNPs
#mv ${PROJ}_filteredSNPs.recode.vcf variantcalling/${PROJ}_filteredSNPs.vcf




#look up {skater} package
bgzip manu2_filtered_SNPs.vcf
tabix -p vcf manu2_filtered_SNPs.vcf.gz
bcftools sort -Ov -m 50G -o manu_sorted.vcf manu2_filtered_SNPs.vcf.gz

bcftools view -H manu_sorted.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > manu.chrom-map.txt
vcftools --vcf manu_sorted.vcf --plink --chrom-map manu.chrom-map.txt --out manu2ped


NUMPROC=40
MEMORY=250
PROJ=${PWD##*/}

java -Xmx${MEMORY}g -jar /10tb_leviathan/evan/bin/hap-ibd.jar \
gt=/10tb_leviathan/evan/manu2/skater/manu_sorted.vcf \
map=/10tb_leviathan/evan/manu2/skater/manu2ped.map \
out=manuibd nthreads=$NUMPROC

plink --vcf manu_sorted.vcf --recode A --out sequoia_input

###############_NCBI_REFERENCE_STACKS_######################
THECHOSEN=(` find ./assembly/ -maxdepth 2 -name "*.fa*" -type f `)
if [ ! -d "spades_mapping" ]; then mkdir spades_mapping; fi
if [ ! -d "spades_mapping/logs" ]; then mkdir spades_mapping/logs; fi

LIBS=(` find  ./trimming  -name "*R*.fq.gz" | cut -d "/" -f3 `)
POPS=(` find  ./trimming  -name "*.fq.gz" | cut -d "_" -f1-2 | cut -d "/" -f3 | cut -d "-" -f2 | sort | uniq `)
optA=1
optB=4
optO=6
SPLITS=2


cp $THECHOSEN ./spades_mapping/${PROJ}_ncbi_reference.fasta

parallel -j $NUMPROC bwa index {} ::: \
` find ./spades_mapping ./trimming -name "*.fasta" -o -name "*.fq.gz" `
#` find ./trimming ./spades_mapping -name "*.f*.gz" -o -name "*.fasta" `

    #PAIRED
    parallel -j "$(( ${#LIBS[@]} / ($SPLITS * 2) ))" \
        "bwa mem ./spades_mapping/${PROJ}_ncbi_reference.fasta {1} {2} "\
        "-L 20,5 -t $(( $NUMPROC / $((${#LIBS[@]} / ($SPLITS * 2))) )) "\
        "-M -T 35 -A $optA -B $optB -O $optO -R "\
        " \"@RG\\tID:{3}\\tSM:{3}\\tPL:Illumina\" 2> ./spades_mapping/logs/bwa.{3}.log "\
        " | mawk '\$6 !~/[2-9].[SH]/ && \$6 !~ /[1-9][0-9].[SH]/' "\
        " | samtools view -@ $(( $NUMPROC / $SPLITS )) -q 1 -SbT ./spades_mapping/${PROJ}_reference.fasta - > ./spades_mapping/{3}.bam 2> ./spades_mapping/logs/{3}.bwa.log" ::: \
         ` find  ./trimming  -name "*val_1*.f*.gz" | sort ` :::+ \
         ` find  ./trimming  -name "*val_2*.f*.gz" | sort ` :::+ \
         ` find  ./trimming  -name "*val_1*.f*.gz" | sort | cut -d "/" -f3 | cut -d "_" -f1-2 `

    parallel -j "$(( ${#LIBS[@]} / ($SPLITS * 4) ))" samtools sort -@ "$(( $NUMPROC / $((${#LIBS[@]} / ($SPLITS * 4))) ))" -O BAM -o {.}.sort.bam {.}.bam ::: \
    ` find ./spades_mapping -name "*.bam" -not -name "*sort*" | sort `


POPS=(` find  ./trimming  -name "*val_1*.f*.gz" | sort | cut -d "/" -f3 | cut -d "_" -f1-2 | sort | uniq `)

printf "%s\n" "${POPS[@]}" > ./spades_mapping/poplist
paste ./spades_mapping/poplist ./spades_mapping/poplist > ./spades_mapping/popmap
ref_map.pl -T $NUMPROC --popmap ./spades_mapping/popmap -o ./stacks/ --samples ./spades_mapping



###################################
#moving forward with best mapped genome
###################################

#
NUMPROC=40
MEMORY=500
PROJ=${PWD##*/}
ref=mapping/spades_ref.fasta

parallel -j 20 samtools sort -@ 2 -O BAM -o {.}.sort.bam {.}.bam ::: \
    ` find ./mapping -name "*.bam" -not -name "*sort*" | sort ` 

parallel -j $NUMPROC samtools index {} ::: ` find ./mapping -name "*.sort.bam" `
find ./mapping -name "*.sort.bam" > ./variantcalling/bams.list
samtools merge -@ $NUMPROC -b ./variantcalling/bams.list -f ./variantcalling/merged_pools.bam
samtools index ./variantcalling/merged_pools.bam
samtools faidx $ref

###############
POPS=(` find  ./mapping  -name "*.bam" | sort | cut -d "/" -f3 | cut -d "." -f1 | sort | uniq `)

bamToBed -i ./variantcalling/merged_pools.bam | bedtools merge -i - > ./variantcalling/mapped.bed

split -n l/$NUMPROC -d --additional-suffix=.bed ./variantcalling/mapped.bed ./variantcalling/mapped.


# bedtools coverage -b ./variantcalling/merged_pools.bam -a ./variantcalling/mapped.bed -counts -sorted > ./variantcalling/cov.stats

# DP=$(mawk '{print $4}' ./variantcalling/cov.stats | sort -rn | perl -e '$d=.001;@l=<>;print $l[int($d*@l)]')
# CC=$( mawk -v x=$DP '$4 < x' ./variantcalling/cov.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$(( $NUMPROC * 2 ))'"}')

# find ./trimming/ -name "*.fq.gz" | cat | parallel -j $NUMPROC "gunzip -c {} | head -2 | tail -1 >> ./variantcalling/lengths.txt"
# MaxLen=$(mawk '{ print length() | "sort -rn" }' ./variantcalling/lengths.txt | head -1)
# MaxLen2=$(( $MaxLen / 2 ))
# ML1=$(( $MaxLen2 + 1 ))
# MaxCov=$(mawk '{print $4}' ./variantcalling/cov.stats | sort -rn | head -1)
# MaxCov2=$(( $MaxCov * 8 / 10 ))

# mawk '$4 > $MaxCov2' ./variantcalling/cov.stats > ./variantcalling/cov.high.stats
# mawk '$4 <= $MaxCov2' ./variantcalling/cov.stats > ./variantcalling/cov.low.stats

# #split high coverage intervals into smaller, 1/2 read-length sized intervals

# bedtools makewindows -b ./variantcalling/cov.high.stats -w $MaxLen2 -s $ML1 > ./variantcalling/temp1.bed
# bedtools intersect -a ./variantcalling/cov.stats -b ./variantcalling/temp1.bed > ./variantcalling/temp.cov.stats

# cat ./variantcalling/temp.cov.stats ./variantcalling/cov.low.stats > ./variantcalling/cov.split.stats
# rm ./variantcalling/cov.high.stats ./variantcalling/cov.low.stats ./variantcalling/temp1.bed ./variantcalling/temp.cov.stats

# TT=$(( $MaxLen2 * 1000000 ))

# mawk -v x=$DP '$4 < x' ./variantcalling/cov.split.stats |sort -V -k1,1 -k2,2 | mawk -v cutoff=$CC -v tt=$TT 'BEGIN{i=1}
# {len=$3-$2;lc=len*$4;cov = cov + lc
# if (NR == 1 && lc > tt) {x="./variantcalling/mapped."i".bed";print $1"\t"$2"\t"$3 > x; i=i+1; e=1}
# else if ( cov < cutoff && lc < tt) {x="./variantcalling/mapped."i".bed";print $1"\t"$2"\t"$3 > x; e=0}
# else if (lc > tt && e > 0 ) {x="./variantcalling/mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=0;i=i+1; e=1}
# else if (lc > tt && e < 1 ) {i=i+1; x="./variantcalling/mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=0;i=i+1;e=1}
# else if (cov > cutoff && lc < tt ) {i=i+1; x="./variantcalling/mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=lc;e=0}
# }'


wait
####################
POPS=(` find  ./mapping  -name "*.bam" | sort | cut -d "/" -f3 | cut -d "." -f1 | sort | uniq `)

printf "%s\n" "${POPS[@]}" > poplist
paste poplist poplist > popmap

if [ ! -f "./variantcalling/reference.fasta" ]
    then cp $ref ./variantcalling/reference.fasta
    fi

if [ ! -d "./variantcalling/raw" ]; then
     mkdir ./variantcalling/raw
fi

if [ ! -d "./variantcalling/filtered" ]; then
     mkdir ./variantcalling/filtered
fi

# VT1="--max-missing 0.9"
# VT2="--mac 2"
# VT3="--minQ 20"
# VT4="--minDP 10"
# #VT5="--min-meanDP 5"
# #VT6="--max-meanDP 500" 

# VL1=("MQM > 39 & MQMR > 39")
# VL2=("MQM / MQMR > 0.75 & MQM / MQMR < 1.25")
# VL3=("RPR > 0 & RPL > 0")
# VL4=("QUAL / DP > 0.25")
# VL5=("PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05")
# VL6=("NS > 1")
# VL7=("LEN < 11")

   call_genos(){
    if [ ! -f "./variantcalling/split.$1.bam" ]
    then samtools view -b -1 -L ./variantcalling/mapped.$1.bed -o ./variantcalling/split.$1.bam ./variantcalling/merged_pools.bam
    fi
    
    if [ ! -f "./variantcalling/split.$1.bam.bai" ]
    then samtools index ./variantcalling/split.$1.bam
    fi

    if [ ! -f "./variantcalling/reference.fasta.fai" ]
    then samtools faidx ./variantcalling/reference.fasta
    fi

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

    #vcflib vcffilter -f "QUAL > 10" raw.$1.vcf > qualfil.$1.vcf
    #add in all filtering here while the vcf is already split up?!
    vcflib vcffilter -f \
       "QUAL > 10 & \
        DP > 10 & \
        TYPE = snp & \
        AO > 2 & \
        QUAL / DP > 0.25 & \
        PAIRED > 0.5 & \
        PAIREDR > 0.5 & \
        PAIREDR / PAIRED < 1.75 & \
        PAIREDR / PAIRED > 0.25 & \
        MQM > 30 & MQMR > 30 & \
        MQM / MQMR > 0.75 & MQM / MQMR < 1.25 & \
        NS > 2 & \
        LEN < 11" ./variantcalling/raw/raw.$1.vcf | \
    vcftools --vcf - \
        --max-missing 0.2 \
        --max-meanDP 1000 \
        --min-meanDP 3 \
        --recode --recode-INFO-all \
        --stdout > ./variantcalling/filtered/filtered.$1.vcf
    #grep -v "#contig" ./variantcalling/raw.$1.vcf > ./variantcalling/raw/raw_fixed.$1.vcf
    #grep -v "#contig" ./variantcalling/filtered.$1.vcf > ./variantcalling/filtered/filtered_fixed.$1.vcf
    }

    export -f call_genos

ulimit -s unlimited
ls ./variantcalling/mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | cut -d "/" -f3 | shuf | parallel --env call_genos -j $NUMPROC --no-notice call_genos {} 
   


#mv ./mapping/reference.fasta ./mapping/${PROJ}_reference.fasta
#mv ./variantcalling/raw*.vcf ./variantcalling/raw/
#mv ./variantcalling/filtered*.vcf ./variantcalling/filtered/

rename -f -e 's/\d+/sprintf("%02d",$&)/e' -- ./variantcalling/raw/*.vcf
vcflib vcfcombine ./variantcalling/raw/raw.*.vcf > ./variantcalling/TotalRawSNPs.vcf

rename -f -e 's/\d+/sprintf("%02d",$&)/e' -- ./variantcalling/filtered/*.vcf
vcflib vcfcombine ./variantcalling/filtered/filtered.*.vcf > ./variantcalling/filtered_SNPs.vcf

######

#if output is polyploid for non GT calls (ie ././.)
grep -v "^#" ./variantcalling/TotalRawSNPs.vcf | cut -f 10- | cut -d ':' -f 1 | sort | uniq
sed 's/\.\/\.\/\./\.\/\./' ./variantcalling/TotalRawSNPs.vcf > ./variantcalling/TotalRawSNPs_fix.vcf
rm ./variantcalling/TotalRawSNPs.vcf; mv ./variantcalling/TotalRawSNPs_fix.vcf ./variantcalling/${PROJ}_TotalRawSNPs.vcf

#if output is polyploid for non GT calls (ie ././.)
grep -v "^#" ./variantcalling/filtered_SNPs.vcf | cut -f 10- | cut -d ':' -f 1 | sort | uniq
sed 's/\.\/\.\/\./\.\/\./' ./variantcalling/filtered_SNPs.vcf > ./variantcalling/filtered_SNPs_fix.vcf
rm ./variantcalling/filtered_SNPs.vcf; mv ./variantcalling/filtered_SNPs_fix.vcf ./variantcalling/${PROJ}_filtered_SNPs.vcf


#vcffiltering
vcfin="final_filter_spades_denovo_light_filters"
vcfout="spades_20perc_miss"
missingness=0.8 #0-1 with 0 allowing for all missing data at a site and 1 allowing for no missing data
vcftools --vcf ${vcfin}.vcf --max-missing $missingness --recode --recode-INFO-all --stdout > ${vcfout}.vcf
vcfin="spades_20perc_miss"
vcfout="spades_20perc_5mmd"
vcftools --vcf ${vcfin}.vcf --min-meanDP 5 --recode --recode-INFO-all --stdout >  ${vcfout}.vcf

vcfin="final_filter_ref2_hleu_light_filters"
vcfout="ref2_20perc_miss"
missingness=0.8 #0-1 with 0 allowing for all missing data at a site and 1 allowing for no missing data
vcftools --vcf ${vcfin}.vcf --max-missing $missingness --recode --recode-INFO-all --stdout > ${vcfout}.vcf
vcfin="ref2_20perc_miss"
vcfout="ref2_20perc_5mmd"
vcftools --vcf ${vcfin}.vcf --min-meanDP 5 --recode --recode-INFO-all --stdout >  ${vcfout}.vcf


#prep sequoia input

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
