#!/bin/bash
#SBATCH --job-name=distance_test
#SBATCH --output=debugging/distance_test_%j.out
#SBATCH --error=debugging/distance_test_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH --account=r00582

# Distance Computation Verification Test for Issue 59
# Tests whether riemtan's compute_dists() accurately computes distances

echo "=========================================="
echo "Distance Computation Verification Test"
echo "Job ID: $SLURM_JOB_ID"
echo "Started at: $(date)"
echo "=========================================="
echo ""

# Load required modules (adjust based on your cluster configuration)
# Uncomment and modify as needed:
# module load R/4.4.1
# module load gcc/11.2.0

# Set working directory to package root
cd $SLURM_SUBMIT_DIR

# Ensure debugging output directory exists
mkdir -p debugging

# Print R version and library paths
echo "R Configuration:"
Rscript --version
echo ""
echo "Library paths:"
Rscript -e ".libPaths()"
echo ""

# Print loaded package versions
echo "Package versions:"
Rscript -e "packageVersion('riemtan')"
Rscript -e "packageVersion('riemstats')"
Rscript -e "packageVersion('Matrix')"
echo ""

# Run the distance computation test
echo "Running distance computation verification..."
echo "=========================================="
Rscript debugging/test_distance_computation.R

# Capture exit code
EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job completed at: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Status: SUCCESS"
else
    echo "Status: FAILED"
fi

echo "=========================================="
echo ""
echo "Output files:"
echo "  - Log: debugging/distance_test_${SLURM_JOB_ID}.out"
echo "  - Errors: debugging/distance_test_${SLURM_JOB_ID}.err"
echo "  - Plot: debugging/distance_verification_qq.png"

exit $EXIT_CODE
