#!/bin/bash
# ==============================================================================
# Script: 02_runLIONS_simulation.sh
# Purpose: Execute the LIONS pipeline on simulated data and extract junctions 
#          for fair benchmarking against TEDDY.
# ==============================================================================
set -euo pipefail

WD="$(pwd)"
echo "[LIONS Benchmark] Working directory: $WD"

# ------------------------------------------------------------------------------
# 1. Reference Formatting (Converting standard GTF to LIONS custom format)
# ------------------------------------------------------------------------------
echo "[LIONS Benchmark] Converting Ground Truth GTF to LIONS exon format..."
mkdir -p LIONS

awk '
  $0 !~ /^#/ && $3 == "exon" {
    split($9, a, ";");
    tid = "NA"; exon_num = "0"; gid = "NA"; gname = "NA";
    for (i in a) {
      if (a[i] ~ /transcript_id/) {
        match(a[i], /"[^"]+"/);
        if (RSTART > 0) tid = substr(a[i], RSTART+1, RLENGTH-2);
      }
      if (a[i] ~ /exon_number/) {
        match(a[i], /[0-9]+/);
        if (RSTART > 0) exon_num = substr(a[i], RSTART, RLENGTH);
      }
      if (a[i] ~ /gene_id/) {
        match(a[i], /"[^"]+"/);
        if (RSTART > 0) gid = substr(a[i], RSTART+1, RLENGTH-2);
      }
      if (a[i] ~ /gene_name/) {
        match(a[i], /"[^"]+"/);
        if (RSTART > 0) gname = substr(a[i], RSTART+1, RLENGTH-2);
      }
    }
    print $1, $4, $5, tid, ".", $7, exon_num, gid, gname;
  }
' OFS="\t" 1000_simulated_transcripts.gtf > LIONS/1000_simulated_transcripts.exons_exons.clean

# ------------------------------------------------------------------------------
# 2. LIONS Core Execution (Processing simulated BAMs)
# ------------------------------------------------------------------------------
echo "[LIONS Benchmark] Running chimericReadSearch.py across depths..."
LIONS_SCRIPT="/mnt/datadisk/xiaoyihan/software/LIONS-master/scripts/ChimericReadTool/chimericReadSearch.py"
REPEAT_BED="/mnt/datadisk/xiaoyihan/index/repeats_mm10_3.bed"
BAM_DIR="sortbam"
OUT_DIR="sortbam/LIONS"
mkdir -p $OUT_DIR

# Using a simplified array for standard depths and replicates
for depth in 5 10 25 50 100; do
  for rep in "" "_rep2"; do
    SAMPLE="simulated_${depth}x${rep}_noise0.1_output"
    echo "  -> Running LIONS for: $SAMPLE"
    
    python $LIONS_SCRIPT \
      LIONS/1000_simulated_transcripts.exons_exons.clean \
      $REPEAT_BED \
      $BAM_DIR/${SAMPLE}.bam \
      $OUT_DIR/${SAMPLE}.chimeras
  done
done

# ------------------------------------------------------------------------------
# 3. Post-Processing: Extracting Junctions for Fair Comparison
# ------------------------------------------------------------------------------
echo "[LIONS Benchmark] Extracting splice junctions and merging replicates..."
cd $OUT_DIR

declare -A PAIRS=(
  [5x]="simulated_5x_noise0.1_output simulated_5x_rep2_noise0.1_output"
  [10x]="simulated_10x_noise0.1_output simulated_10x_rep2_noise0.1_output"
  [25x]="simulated_25x_noise0.1_output simulated_25x_rep2_noise0.1_output"
  [50x]="simulated_50x_noise0.1_output simulated_50x_rep2_noise0.1_output"
  [100x]="simulated_100x_noise0.1_output simulated_100x_rep2_noise0.1_output"
)

extract_junctions() {
  local chim="$1"
  local out="$2"
  awk -F'\t' 'BEGIN{OFS="\t"}
    NF>=12 {
      chrom=$1; chromStart=$2;
      split($11,bs,","); split($12,st,",");
      seg1_end = chromStart + st[1] + bs[1];
      seg2_start = chromStart + st[2];
      start = seg1_end - 1; if(start < 0) start = 0;
      print chrom, start, seg2_start, $4
    }
  ' "$chim" > "$out"
}

for depth in "${!PAIRS[@]}"; do
  samples=${PAIRS[$depth]}
  out_merged="${depth}.junction.bed"
  tmp_list=()
  
  for s in $samples; do
    if [ -f "${s}.chimeras" ]; then
      extract_junctions "${s}.chimeras" "${s}.junction.bed"
      tmp_list+=("${s}.junction.bed")
    fi
  done

  if [ ${#tmp_list[@]} -gt 0 ]; then
    cat "${tmp_list[@]}" | sort -k1,1 -k2,2n -k3,3n | uniq > "$out_merged"
    echo "  -> Created $out_merged with $(wc -l < "$out_merged") unique junctions."
  fi
done

echo "[LIONS Benchmark] Pipeline Complete."