#!/bin/bash
#SBATCH --job-name=oracle_test
#SBATCH --output=debugging/oracle_test_%j.out
#SBATCH --error=debugging/oracle_test_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH --account=r00582

# Oracle Test with True Means for Issue 59
# Tests frechet_anova formula using TRUE means to isolate riemstats vs riemtan issues

echo "=========================================="
echo "Fréchet ANOVA Oracle Test"
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

# Run the oracle test
echo "Running oracle test (this may take a while with 1000 replicates)..."
echo "=========================================="
module load r/4.5.1
Rscript debugging/test_frechet_oracle.R

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
echo "  - Log: debugging/oracle_test_${SLURM_JOB_ID}.out"
echo "  - Errors: debugging/oracle_test_${SLURM_JOB_ID}.err"
echo "  - Plot: debugging/oracle_test_qq.png"
echo ""
echo "To view results:"
echo "  cat debugging/oracle_test_${SLURM_JOB_ID}.out"
echo "  # Look for 'FINAL DIAGNOSTIC VERDICT' section"

exit $EXIT_CODE
