#!/bin/bash
#
#SBATCH -J 87.2.7
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#SBATCH -D ./
#SBATCH --mail-user ryanjsfx@mpa-garching.mpg.de
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBTACH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn49ye
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:29:50
#SBATCH --partition=micro
#SBATCH --ear=off

module load spack/22.2.1
module load intel-oneapi-toolkit/2022.3.0
module load hdf5/1.8.22-intel21-impi
module list

echo "PWD: $PWD"

date
srun ./athena -i ../tst/megKH/athinputmeg.kh
date
