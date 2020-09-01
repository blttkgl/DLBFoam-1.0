#!/bin/bash
#SBATCH -p batch
#SBATCH -t 00:30:00
#SBATCH -n 4

#set -euo pipefail

. /share/apps/spack/envs/fgci-centos7-generic/software/openfoam-org/8/fadnomj/etc/bashrc
module load gcc

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


