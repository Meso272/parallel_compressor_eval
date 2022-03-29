#!/bin/bash
#SBATCH --job-name=p1024
#SBATCH -p bdwall
#SBATCH -A ECP-EZ
#SBATCH --nodes 32
#SBATCH --ntasks-per-node=32
#SBATCH --time=5:00:00
#SBATCH -o p1024.%j.%N.out
#SBATCH -e p1024.%j.%N.error

echo wait .... 30 seconds
sleep 30
echo date=`date`


# Hurricane
srun ./parallel_sz 13 496 496 96
srun ./parallel_zfp 13 496 496 96
srun ./parallel_sz3 sz.config 13 496 496 96
srun ./parallel_qoz sz.config 13 496 496 96

# miranda
srun ./parallel_sz  7 384 384 256
srun ./parallel_zfp 7 384 384 256
srun ./parallel_sz3 sz.config 7 384 384 256
srun ./parallel_qoz sz.config 7 384 384 256

