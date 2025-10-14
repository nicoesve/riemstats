#!/bin/bash
#SBATCH --job-name=decomp_test
#SBATCH --output=debugging/decomp_test_%j.out
#SBATCH --error=debugging/decomp_test_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=general
#SBATCH --account=r00582

# Decomposition Analysis for Fréchet ANOVA (Issue 59)
# Tests asymptotic behavior: Term1 → χ²(k-1), Term2 → 0

echo "=========================================="
echo "Fréchet ANOVA Asymptotic Decomposition Test"
echo "Job ID: $SLURM_JOB_ID"
echo "Started at: $(date)"
echo "=========================================="
echo ""

# Set working directory to package root
cd $SLURM_SUBMIT_DIR

# Ensure debugging output directory exists
mkdir -p debugging

# Print R version and library paths
echo "R Configuration:"
module load r/4.5.1
Rscript --version
echo ""
echo "Library paths:"
Rscript -e ".libPaths()"
echo ""

# Print loaded package versions
echo "Package versions:"
Rscript -e "cat('riemtan: '); packageVersion('riemtan')"
Rscript -e "cat('riemstats: '); packageVersion('riemstats')"
Rscript -e "cat('Matrix: '); packageVersion('Matrix')"
echo ""

# Display system info
echo "System information:"
echo "  Hostname: $(hostname)"
echo "  CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "  Memory allocated: $SLURM_MEM_PER_NODE MB"
echo ""

# Run the decomposition test
echo "Running decomposition test..."
echo "=========================================="
Rscript debugging/test_frechet_decomposition.R

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
echo "  - Log: debugging/decomp_test_${SLURM_JOB_ID}.out"
echo "  - Errors: debugging/decomp_test_${SLURM_JOB_ID}.err"
echo "  - Plot: debugging/decomposition_diagnostics.png"
echo ""
echo "To view results:"
echo "  cat debugging/decomp_test_${SLURM_JOB_ID}.out"
echo "  # Look for 'DIAGNOSTIC VERDICT' section"

exit $EXIT_CODE
