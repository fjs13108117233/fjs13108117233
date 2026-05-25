#!/bin/bash
#===============================================================================
# GATK Call SNP Pipeline for KN5585 (Paired-end)
# Sample: KN5585
# Input: R1.fq, R2.fq (paired-end reads)
# Reference: maize V4
#===============================================================================

set -e

# ===== Parameters =====
SAMPLE="KN5585"
WORKDIR="/data/lhz/Fang/GS/KN5585"
REF="/lthpcfs/share/db/maizev4/maizev4.fa"
BWA_INDEX="/lthpcfs/share/db/maizev4/maizev4"
THREADS=8
JAVA_MEM="-Xmx55g"

cd ${WORKDIR}

# ===== Step 1: Quality Control (fastp) =====
echo "[$(date)] Step 1: Running fastp..."
fastp \
    -i R1.fq \
    -I R2.fq \
    -o ${SAMPLE}_clean_R1.fq.gz \
    -O ${SAMPLE}_clean_R2.fq.gz \
    -w 10 \
    -l 35 \
    -h ${SAMPLE}_fastp.html \
    -j ${SAMPLE}_fastp.json

# ===== Step 2: BWA Alignment =====
echo "[$(date)] Step 2: BWA alignment..."
bwa mem -t ${THREADS} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:illumina\tPU:unit1" \
    ${BWA_INDEX} \
    ${SAMPLE}_clean_R1.fq.gz ${SAMPLE}_clean_R2.fq.gz \
    | samtools sort -@ ${THREADS} -o ${SAMPLE}.sorted.bam -

# ===== Step 3: Index BAM =====
echo "[$(date)] Step 3: Indexing BAM..."
samtools index ${SAMPLE}.sorted.bam

# ===== Step 4: Mark Duplicates =====
echo "[$(date)] Step 4: Marking duplicates..."
gatk ${JAVA_MEM} MarkDuplicates \
    -I ${SAMPLE}.sorted.bam \
    -O ${SAMPLE}.markdup.bam \
    -M ${SAMPLE}_dup_metrics.txt \
    --REMOVE_DUPLICATES false

samtools index ${SAMPLE}.markdup.bam

# ===== Step 5: HaplotypeCaller (GVCF mode) =====
echo "[$(date)] Step 5: HaplotypeCaller..."
gatk ${JAVA_MEM} HaplotypeCaller \
    -R ${REF} \
    --emit-ref-confidence GVCF \
    -I ${SAMPLE}.markdup.bam \
    -O ${SAMPLE}.g.vcf.gz \
    --native-pair-hmm-threads ${THREADS}

# ===== Step 6: GenotypeGVCFs =====
echo "[$(date)] Step 6: GenotypeGVCFs..."
gatk ${JAVA_MEM} GenotypeGVCFs \
    -R ${REF} \
    -V ${SAMPLE}.g.vcf.gz \
    -O ${SAMPLE}.vcf.gz

# ===== Step 7: Hard Filtering =====
echo "[$(date)] Step 7: Variant filtration..."
# SNP filtering
gatk ${JAVA_MEM} SelectVariants \
    -R ${REF} \
    -V ${SAMPLE}.vcf.gz \
    --select-type-to-include SNP \
    -O ${SAMPLE}.snp.vcf.gz

gatk ${JAVA_MEM} VariantFiltration \
    -R ${REF} \
    -V ${SAMPLE}.snp.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_filter" \
    -O ${SAMPLE}.snp.filtered.vcf.gz

# Extract PASS variants only
gatk ${JAVA_MEM} SelectVariants \
    -R ${REF} \
    -V ${SAMPLE}.snp.filtered.vcf.gz \
    --exclude-filtered \
    -O ${SAMPLE}.snp.PASS.vcf.gz

echo "[$(date)] Done! Final SNP file: ${SAMPLE}.snp.PASS.vcf.gz"
echo "[$(date)] Total SNPs:"
bcftools stats ${SAMPLE}.snp.PASS.vcf.gz | grep "number of SNPs:"
