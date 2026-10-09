#!/bin/bash
# Snakemake head job: submits every rule as its own Slurm job and waits.
# From the RUN project root (logs/ is created there):
#   export BULLPEN_REPO=/path/to/smORF-Bullpen
#   sbatch --account=<acct> --partition=<part> \
#       $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/bullpen.yaml   # add -n to dry-run
# Uses the snakemake on PATH (workflow/envs/snakemake.yaml). Set PROFILE to use
# a site copy of this profile, e.g. PROFILE=$BULLPEN_REPO/profiles/mysite.
#SBATCH --job-name=bullpen_head
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
CONFIG=${1:?usage: sbatch run_head.sbatch <config.yaml> [extra snakemake args]}
shift
# sbatch runs a spooled copy of this script, so the repo cannot be found from $0.
REPO=${BULLPEN_REPO:?export BULLPEN_REPO=/path/to/smORF-Bullpen before sbatch}
PROFILE=${PROFILE:-$REPO/profiles/slurm}

# Site setup goes here if your cluster needs it, e.g.
#   module load apptainer
#   conda activate bullpen-snakemake
command -v snakemake >/dev/null || { echo "snakemake not on PATH" >&2; exit 1; }
command -v apptainer >/dev/null || command -v singularity >/dev/null \
    || { echo "apptainer/singularity not on PATH" >&2; exit 1; }

export TMPDIR=${TMPDIR:-/tmp}
mkdir -p logs

snakemake -s "$REPO/workflow/Snakefile" --configfile "$CONFIG" \
    --profile "$PROFILE" "$@"
