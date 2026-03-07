#!/bin/bash
# ==============================================================================
# Script: runRSEMsimulation.sh
# Purpose: Generate realistic simulated RNA-seq reads using empirical models.
# Note: The background sequencing model is trained on real mouse 2-cell stage RNA-seq
#       data to accurately mimic authentic error profiles and read distributions.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Train Background Model on Real RNA-seq Data
# ------------------------------------------------------------------------------
# Using real 2-cell stage fastq files to learn the empirical expression distribution
nohup ../../RSEM/rsem-calculate-expression --paired-end -p 40 --bowtie2 \
    /mnt/datadisk/xiaoyihan/1129/RNA/chenfei/fastq/2cellrep1_1.fastq \
    /mnt/datadisk/xiaoyihan/1129/RNA/chenfei/fastq/2cellrep1_2.fastq \
    /mnt/datadisk/xiaoyihan/manuals/compare/simulate/forStat/forStat &

# Wait for background model generation to complete before proceeding...

# ------------------------------------------------------------------------------
# 2. Prepare RSEM Reference 
# ------------------------------------------------------------------------------
# Build the reference indices based on the sampled 1000 transcripts GTF
../../RSEM/rsem-prepare-reference --bowtie2 \
    --gtf 1000_simulated_transcripts.gtf \
    /mnt/datadisk/xiaoyihan/index/mm10_no_alt_analysis_set_ENCODE.fasta \
    mm10_simulate_1000_reference

# ------------------------------------------------------------------------------
# 3. Simulate Reads Across Sequencing Depths (5x to 100x)
# ------------------------------------------------------------------------------
# Calculate base length from the training results
total_length=$(awk '{if(NR>1) sum += $3} END {print sum}' 1000simulate.isoforms.results)
read_length=150

# Loop through predefined sequencing depths
for depth in 5 10 25 50 100; do
    # Calculate the exact number of reads required for the target depth
    total_reads=$(echo "$depth * $total_length / $read_length" | bc)
    echo "Processing Depth ${depth}x: Generating ${total_reads} reads..."

    # Run RSEM simulation with a noise fraction of 0.1
    ../../RSEM/rsem-simulate-reads \
        /mnt/datadisk/xiaoyihan/manuals/compare/simulate/forStat/forStat \
        /mnt/datadisk/xiaoyihan/manuals/compare/simulate/forStat/forStat.stat/forStat.model \
        1000simulate.isoforms.results \
        0.1 \
        ${total_reads} \
        simulated_${depth}x_noise0.1_output
done

