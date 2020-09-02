#!/bin/bash
#SBATCH --job-name=example
#SBATCH --account=XXX
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64

. /appl/soft/eng/OpenFOAM/OpenFOAM-8/etc/bashrc

# Run all clean
./Allclean

# Get the processor count
((N_PROCS=SLURM_NTASKS ))

# Set the processor count on blockMeshDict and decomposeParDict
foamDictionary system/decomposeParDict -entry nProc -set $N_PROCS
foamDictionary system/blockMeshDict -entry nProc -set $N_PROCS

# Set cells per core for benchmark
CELLS_PER_CORE=20
foamDictionary system/blockMeshDict -entry cellsPerCore -set $CELLS_PER_CORE

# Blockmesh, decomposePar
blockMesh > log.blockMesh
decomposePar > log.decompose

# Run the benchmark
srun chemistryBenchmark -parallel > log.chemistryBenchmark




