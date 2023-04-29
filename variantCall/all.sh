#!/bin/bash 
#++++++++++++++++++++++++++++++++++
FASTQ_DIR=$(pwd)
ref_fasta="MpTak_v6.1.genome.fasta"
fastq_file=($(ls ${FASTQ_DIR}/*_R1.fastq.gz))
base_name=${fastq_file%_R1.fastq.gz}
#++++++++++++++++++++++++++++++++++

## Trimming
fastp \
  -i ${base_name}_R1.fastq.gz \
  -I ${base_name}_R2.fastq.gz \
  -o ${base_name}_R1_trim.fastq.gz \
  -O ${base_name}_R2_trim.fastq.gz \
  -p 8 \
  -j ${base_name}_fastp.json \
  -h ${base_name}_fastp.html

## Mapping, sorting by name with `-n`
bwa-mem2 \
  mem \
  -t 8 \
  -R "@RG\tID:${base_name}\tLB:${base_name}\tPL:ILLUMINA\tPM:NEXTSEQ\tSM:${base_name}" \
  ${ref_fasta} \
  ${base_name}_1.fastq.gz \
  ${base_name}_2.fastq.gz \
  | samtools sort -n -m 8G -@ 8 -O bam -T ${base_name}_tmp1 - -o ${base_name}.bam

## Fixmate run with `-m`
samtools fixmate -m -@ 8 -O bam ${base_name}.bam ${base_name}_fixmate.bam

rm ${base_name}.bam

# sort witOUT `-n`
samtools sort -m 8G -@ 8 -O bam -T ${base_name}_tmp2 ${base_name}_fixmate.bam -o ${base_name}_positionsorted.bam

rm ${base_name}_fixmate.bam

# mark duplicates
samtools markdup -@ 8 -O bam ${base_name}_positionsorted.bam ${base_name}_markdup.bam

rm ${base_name}_positionsorted.bam

## samtools indexing
samtools index -@ 8 ${base_name}_markdup.bam

##=========Comment Start===========
###need -R tag for gatk
###before using fixmate, bam file have to be name-sorted "-n option"
###markdup requires fixmate run with -m option
###before using markdup, bam file have to be position-order "withOUT -n option"
###=========Comment End ============

## Variant call
gatk \
  --java-options "-Xmx16g" \
  HaplotypeCaller \
    -pairHMM LOGLESS_CACHING \
    --native-pair-hmm-threads 8 \
    -R ${ref_fasta} \
    -I ${base_name}_markdup.bam \
    -O ${base_name}_pre_bqsr_hapcall.g.vcf.gz \
    -ERC GVCF \
    -ploidy 1 

## GenotypeGVCFs
gatk \
  --java-options "-Xmx16g" \
    GenotypeGVCFs \
      -R ${ref_fasta} \
      -V ${base_name}_pre_bqsr_hapcall.g.vcf.gz \
      -O ${base_name}_pre_bqsr_genotype.vcf.gz \
      -ploidy 1 

#rm ${base_name}_pre_bqsr_hapcall.g.vcf.gz
#rm ${base_name}_pre_bqsr_hapcall.g.vcf.gz.tbi

## SelctVariants
gatk \
  --java-options "-Xmx16g" \
    SelectVariants \
      -R ${ref_fasta} \
      -V ${base_name}_pre_bqsr_genotype.vcf.gz \
      -select-type SNP \
      -O ${base_name}_pre_bqsr_SNP.vcf.gz 

gatk \
  --java-options "-Xmx16g" \
  SelectVariants \
    -R ${ref_fasta} \
    -V ${base_name}_pre_bqsr_genotype.vcf.gz \
    -select-type INDEL \
    -O ${base_name}_pre_bqsr_indel.vcf.gz

#rm ${base_name}_pre_bqsr_genotype.vcf.gz
#rm ${base_name}_pre_bqsr_genotype.vcf.gz.tbi

## Filter variants
gatk \
  --java-options "-Xmx16g" \
  VariantFiltration \
  -R ${ref_fasta} \
  -V ${base_name}_pre_bqsr_SNP.vcf.gz \
  -O ${base_name}_pre_bqsr_filtered_SNP.vcf.gz \
    -filter "QD<2.0" --filter-name "QD2" \
    -filter "QUAL<30.0" --filter-name "QUAL30" \
    -filter "SOR>3.0" --filter-name "SOR3" \
    -filter "FS>60.0" --filter-name "FS60" \
    -filter "MQ<40.0" --filter-name "MQ40" 

#rm ${base_name}_pre_bqsr_SNP.vcf.gz
#rm ${base_name}_pre_bqsr_SNP.vcf.gz.tbi

gatk \
  --java-options "-Xmx16g" \
  VariantFiltration \
  -R ${ref_fasta} \
  -V ${base_name}_pre_bqsr_indel.vcf.gz \
  -O ${base_name}_pre_bqsr_filtered_indel.vcf.gz \
    -filter "QD<2.0" --filter-name "QD2" \
    -filter "QUAL<30.0" --filter-name "QUAL30" \
    -filter "FS>200.0" --filter-name "FS200"

#rm ${base_name}_pre_bqsr_indel.vcf.gz
#rm ${base_name}_pre_bqsr_indel.vcf.gz.tbi

## VQSR is a method based on machine learning, so when you tackle on SNP detection you are recommended to recruit that. However in this case, marchantia does not possess well-known variants data, hard-filtering is preferable.
    
gatk \
  --java-options "-Xmx16g" \
  SelectVariants \
  --exclude-filtered \
  -V ${base_name}_pre_bqsr_filtered_SNP.vcf.gz \
  -O ${base_name}_for_bqsr_SNP.vcf.gz 

#rm ${base_name}_pre_bqsr_filtered_SNP.vcf.gz
#rm ${base_name}_pre_bqsr_filtered_SNP.vcf.gz.tbi

gatk \
  --java-options "-Xmx16g" \
  SelectVariants \
  --exclude-filtered \
  -V ${base_name}_pre_bqsr_filtered_indel.vcf.gz \
  -O ${base_name}_for_bqsr_indel.vcf.gz

#rm ${base_name}_pre_bqsr_filtered_indel.vcf.gz
#rm ${base_name}_pre_bqsr_filtered_indel.vcf.gz.tbi

## Base Recalibrate                                                             
#this precess requires known SNP/Indel data, whreas Marchantia does not have\      
#If you handle human data, you can merely download it                               
#When you tackle with less-known species, MAKE IT !!                                
#in this case pre-prepared alternative 'hard-filtered' variants datasets'           

gatk \
  --java-options "-Xmx16g" \
  BaseRecalibrator \
    -R ${ref_fasta} \
    -I ${base_name}_markdup.bam \
    --known-sites ${base_name}_for_bqsr_SNP.vcf.gz \
    --known-sites ${base_name}_for_bqsr_indel.vcf.gz \
    -O ${base_name}_recal_data.table 

#rm ${base_name}_for_bqsr_SNP.vcf.gz
#rm ${base_name}_for_bqsr_SNP.vcf.gz.tbi
#rm ${base_name}_for_bqsr_indel.vcf.gz
#rm ${base_name}_for_bqsr_indel.vcf.gz.tbi


## BQSR
gatk \
  --java-options "-Xmx16g" \
  ApplyBQSR \
    -R ${ref_fasta} \
    -I ${base_name}_markdup.bam \
    --bqsr-recal-file ${base_name}_recal_data.table \
    -O ${base_name}_recal.bam

#rm ${base_name}_markdup.bam
#rm ${base_name}_markdup.bam.bai

## option
## 2nd Base recalibration for analysis
gatk \
  --java-options "-Xmx16g" \
  BaseRecalibrator \
    -R ${ref_fasta} \
    -I ${base_name}_recal.bam \
    --known-sites ${base_name}_for_bqsr_SNP.vcf.gz \
    --known-sites ${base_name}_for_bqsr_indel.vcf.gz \
    -O ${base_name}_post_recal_data.table 

## analyze covariates
gatk \
  --java-options "-Xmx16g" \
  AnalyzeCovariates \
    -before ${base_name}_recal_data.table \
    -after ${base_name}_post_recal_data.table \
    -plots ${base_name}_recalibration_plots.pdf

#===haplotypecall===
gatk \
  --java-options "-Xmx16g" \
    HaplotypeCaller \
      -pairHMM LOGLESS_CACHING \
      --native-pair-hmm-threads 8 \
      -R ${ref_fasta} \
      -I ${base_name}_recal.bam \
      -O ${base_name}_raw_variants_recal.g.vcf.gz \
      -ERC GVCF \
      -ploidy 1

#===genotypeGVCF===
gatk \
  --java-options "-Xmx16g" \
  GenotypeGVCFs \
    -R ${ref_fasta} \
    -V ${base_name}_raw_variants_recal.g.vcf.gz \
    -O ${base_name}_raw_variants_recal.vcf.gz \
    -ploidy 1

#rm ${base_name}_raw_variants_recal.g.vcf.gz
#rm ${base_name}_raw_variants_recal.g.vcf.gz.tbi

#=== Select Variants 2nd ===
gatk \
  --java-options "-Xmx16g" \
  SelectVariants \
    -R ${ref_fasta} \
    -V ${base_name}_raw_variants_recal.vcf.gz \
    -select-type SNP \
    -O ${base_name}_raw_SNP_recal.vcf.gz

gatk \
  --java-options "-Xmx16g" \
  SelectVariants \
    -R ${ref_fasta} \
    -V ${base_name}_raw_variants_recal.vcf.gz \
    -select-type INDEL \
    -O ${base_name}_raw_indel_recal.vcf.gz

#rm ${base_name}_raw_variants_recal.vcf.gz
#rm ${base_name}_raw_variants_recal.vcf.gz.tbi

#=== Variants filter 2nd ===
gatk \
  --java-options "-Xmx16g" \
  VariantFiltration \
    -R ${ref_fasta} \
    -V ${base_name}_raw_SNP_recal.vcf.gz \
    -O ${base_name}_filtered_SNP_recal.vcf.gz \
    -filter "QD<2.0" --filter-name "QD2" \
    -filter "QUAL<30.0" --filter-name "QUAL30" \
    -filter "SOR>3.0" --filter-name "SOR3" \
    -filter "FS>60.0" --filter-name "FS60" \
    -filter "MQ<40.0" --filter-name "MQ40" 

gatk \
  --java-options "-Xmx16g" \
  VariantFiltration \
    -R ${ref_fasta} \
    -V ${base_name}_raw_indel_recal.vcf.gz \
    -O ${base_name}_filtered_indel_recal.vcf.gz \
    -filter "QD<2.0" --filter-name "QD2" \
    -filter "QUAL<30.0" --filter-name "QUAL30" \
    -filter "FS>200.0" --filter-name "FS200" \

#rm ${base_name}_raw_indel_recal.vcf.gz
#rm ${base_name}_raw_indel_recal.vcf.gz.tbi
#rm ${base_name}_raw_SNP_recal.vcf.gz
#rm ${base_name}_raw_SNP_recal.vcf.gz.tbi

# Remove under threshold variants
gatk \
  --java-options "-Xmx16g" \
  SelectVariants \
    --exclude-filtered \
    -V ${base_name}_filtered_SNP_recal.vcf.gz \
    -O ${base_name}_filtered_PASS_SNP_recal.vcf.gz
gatk \
  --java-options "-Xmx16g" \
  SelectVariants \
    --exclude-filtered \
    -V ${base_name}_filtered_indel_recal.vcf.gz \
    -O ${base_name}_filtered_PASS_indel_recal.vcf.gz

#rm ${base_name}_filtered_SNP_recal.vcf.gz
#rm ${base_name}_filtered_SNP_recal.vcf.gz.tbi
#rm ${base_name}_filtered_indel_recal.vcf.gz
#rm ${base_name}_filtered_indel_recal.vcf.gz.tbi

# Make .csi index for bcftools to cocatenate SNP snd indel vcf
bcftools index ${base_name}_filtered_PASS_SNP_recal.vcf.gz
bcftools index ${base_name}_filtered_PASS_indel_recal.vcf.gz

# Concatenate
bcftools concat \
  ${base_name}_filtered_PASS_SNP_recal.vcf.gz \
  ${base_name}_filtered_PASS_indel_recal.vcf.gz \
  --allow-overlaps \
  -O z --output ${base_name}_filtered_PASS_variants_recal.vcf.gz

# Version Info
fastp --version
bwa-mem2 version
samtools --version
gatk --version
bcftools --version