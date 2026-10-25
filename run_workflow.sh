#!/usr/bin/env bash
# =============================================================================
# run_workflow.sh: run the workflow without the SLURM executor plugin
#   step 1 (login node): GENCODE download + SRA prefetch (needs internet)
#   step 2 (one SLURM allocation): everything else, scheduled by snakemake
#
# usage, from the project directory on a login node (ideally inside tmux):
#   module load python star sratoolkit samtools pigz multiqc
#   bash run_workflow.sh            # both steps
#   bash run_workflow.sh submit     # step 2 only (downloads already done)
#
# version: 1
# date:    2026-10-25
# =============================================================================
set -euo pipefail

CORES=32
MEM_GB=200
TIME="3-00:00:00"      # check the general partition's max walltime
PARTITION="general"

# ---- check tools are on PATH (modules loaded, salmon installed) -------------
missing=0
for t in snakemake tools/bin/salmon STAR prefetch fasterq-dump samtools pigz multiqc; do
  if ! command -v "$t" > /dev/null 2>&1; then
    echo "missing: $t"; missing=1
  fi
done
if [[ $missing -eq 1 ]]; then
  echo "load modules: module load python star sratoolkit samtools pigz multiqc"
  echo "install salmon: bash setup_tools.sh"
  exit 1
fi

mkdir -p logs

# ---- step 1: network steps on the login node --------------------------------
if [[ "${1:-}" != "submit" ]]; then
  snakemake --cores 4 --rerun-incomplete --printshellcmds \
    --until download_gencode prefetch
fi

# ---- step 2: the rest inside one allocation ---------------------------------
# --resources keeps the sum of running rules' mem_mb under the allocation
sbatch -p "$PARTITION" -c "$CORES" --mem "${MEM_GB}G" -t "$TIME" \
  -J rnaseq-smk -o logs/snakemake-%j.out \
  --wrap "snakemake --cores $CORES --resources mem_mb=$(( (MEM_GB - 10) * 1000 )) \
          --rerun-incomplete --keep-going --printshellcmds"
