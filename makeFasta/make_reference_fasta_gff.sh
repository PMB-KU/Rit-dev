#!/bin/bash 
#++++++++++++++++++++++++++++++++++
ref_fasta="MpTak_v6.1.genome.fasta"
feature_vcf_file=($(ls ./*.vcf.gz))
base_name=${feature_vcf_file%%_*}
new_fasta=${base_name}.fa
# This vcf file should be bgzip compressed. not gzip.
# Please use samtools bgzip. You dont need install other bgzip.
#++++++++++++++++++++++++++++++++++

# Commands --------------------------------

gunzip ${feature_vcf_file}

bgzip ${feature_vcf_file%.gz}

picard CreateSequenceDictionary \
    R=${ref_fasta} \
    O=${ref_fasta%.fasta}.dict

gatk --java-options "-Xmx16g" \
    IndexFeatureFile \
    -F ${feature_vcf_file}

gatk --java-options "-Xmx16g" \
    FastaAlternateReferenceMaker \
    -R ${ref_fasta} \
    -V ${feature_vcf_file} \
    -O ${new_fasta}

liftoff \
    -g ${ref_gff} \
    -o ${new_fasta%.fa}.gff \
    -u ${new_fasta%.fa}_unmapped_features.txt \
    ${new_fasta} \
    ${ref_fasta}
