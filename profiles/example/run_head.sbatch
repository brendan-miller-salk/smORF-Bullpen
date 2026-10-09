#!/bin/bash
# EXAMPLE head job for the made-up cluster in profiles/example/config.yaml.
# Copy alongside it into profiles/<site>/ and replace every value marked EDIT.
# From the RUN project root:
#   sbatch <repo>/profiles/<site>/run_head.sbatch config/bullpen.yaml   # add -n to dry-run
#SBATCH --job-name=bullpen_head
#SBATCH --account=mylab_123           # EDIT
#SBATCH --partition=compute           # EDIT
#SBATCH --time=2-00:00:00             # EDIT: at most the partition's limit
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
CONFIG=${1:?usage: sbatch run_head.sbatch <config.yaml> [extra snakemake args]}
shift
# sbatch runs a spooled copy of this script, so the repo cannot be found from $0.
REPO=${BULLPEN_REPO:-/project/mylab_123/smORF-Bullpen}                     # EDIT
SNAKEMAKE=${SNAKEMAKE:-/project/mylab_123/envs/bullpen-snakemake/bin/snakemake}  # EDIT
SITE=mysite                                                                     # EDIT: the profiles/<site> name

# Whatever your cluster needs for apptainer to be on PATH. EDIT
source /etc/profile.d/modules.sh
module purge
module load apptainer

# Keep temp files off small node-local /tmp. EDIT
export TMPDIR=/scratch/$USER/tmp
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=$TMPDIR
mkdir -p "$TMPDIR" "$APPTAINER_CACHEDIR" logs

"$SNAKEMAKE" -s "$REPO/workflow/Snakefile" --configfile "$CONFIG" \
    --profile "$REPO/profiles/$SITE" "$@"
