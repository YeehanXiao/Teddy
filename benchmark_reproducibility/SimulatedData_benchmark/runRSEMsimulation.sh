#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: run_rsem_simulation.sh
# Purpose: Generate simulated paired-end RNA-seq reads using an empirical model
#          trained from real early-embryo RNA-seq data.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. User-defined paths
# ------------------------------------------------------------------------------
RSEM_BIN="../../RSEM"

TRAIN_FASTQ_R1="input/training_fastq/early_embryo_R1.fastq"
TRAIN_FASTQ_R2="input/training_fastq/early_embryo_R2.fastq"

SIMULATED_TRUTH_GTF="input/simulated_truth_transcripts.gtf"
REFERENCE_GENOME_FA="input/reference/mm10_reference_genome.fa"

TRAIN_PREFIX="output/rsem_training/early_embryo_training"
SIM_REFERENCE_PREFIX="output/rsem_reference/mm10_simulated_reference"
SIM_EXPRESSION_PREFIX="output/rsem_expression/simulated_truth_expression"

THREADS=40
READ_LENGTH=150
NOISE_FRACTION=0.1
DEPTHS=(5 10 25 50 100)

mkdir -p output/rsem_training output/rsem_reference output/rsem_expression output/simulated_reads

# ------------------------------------------------------------------------------
# 2. Train empirical model on real RNA-seq data
# ------------------------------------------------------------------------------
"${RSEM_BIN}/rsem-calculate-expression" \
    --paired-end \
    -p "${THREADS}" \
    --bowtie2 \
    "${TRAIN_FASTQ_R1}" \
    "${TRAIN_FASTQ_R2}" \
    "${TRAIN_PREFIX}"

# ------------------------------------------------------------------------------
# 3. Prepare RSEM reference from simulated truth GTF
# ------------------------------------------------------------------------------
"${RSEM_BIN}/rsem-prepare-reference" \
    --bowtie2 \
    --gtf "${SIMULATED_TRUTH_GTF}" \
    "${REFERENCE_GENOME_FA}" \
    "${SIM_REFERENCE_PREFIX}"

# ------------------------------------------------------------------------------
# 4. Estimate expression on the simulated truth reference
# ------------------------------------------------------------------------------
"${RSEM_BIN}/rsem-calculate-expression" \
    --paired-end \
    -p "${THREADS}" \
    --bowtie2 \
    "${TRAIN_FASTQ_R1}" \
    "${TRAIN_FASTQ_R2}" \
    "${SIM_REFERENCE_PREFIX}" \
    "${SIM_EXPRESSION_PREFIX}"

# ------------------------------------------------------------------------------
# 5. Simulate reads across sequencing depths
# ------------------------------------------------------------------------------
ISOFORM_RESULTS="${SIM_EXPRESSION_PREFIX}.isoforms.results"
MODEL_FILE="${TRAIN_PREFIX}.stat/$(basename "${TRAIN_PREFIX}").model"

total_length=$(awk 'NR>1 {sum += $3} END {printf "%.0f\n", sum}' "${ISOFORM_RESULTS}")

for depth in "${DEPTHS[@]}"; do
    total_reads=$(echo "${depth} * ${total_length} / ${READ_LENGTH}" | bc)
    out_prefix="output/simulated_reads/simulated_${depth}x_noise${NOISE_FRACTION}"

    echo "Simulating ${depth}x depth (${total_reads} reads)..."

    "${RSEM_BIN}/rsem-simulate-reads" \
        "${TRAIN_PREFIX}" \
        "${MODEL_FILE}" \
        "${ISOFORM_RESULTS}" \
        "${NOISE_FRACTION}" \
        "${total_reads}" \
        "${out_prefix}"
done