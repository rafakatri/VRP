#!/bin/bash
#SBATCH --job-name=mpi_openmp_10
#SBATCH --nodes=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=espec
#SBATCH --mem=1G

# Capture the start time in milliseconds
start_time=$(($(date +%s%N)/1000000))

echo "Output do job $SLURM_JOB_ID"
echo "mpi_openmp n=10"

# Execute the./global file
mpirun ./mpi_openmp ../grafo_10.txt

# Capture the end time in milliseconds
end_time=$(($(date +%s%N)/1000000))

# Calculate the duration in milliseconds
duration=$((end_time - start_time))

echo "Execution time: $duration milliseconds"
