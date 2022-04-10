#!/bin/bash
#SBATCH --job-name=p32
#SBATCH -p bdwall
#SBATCH -A ECP-EZ
#SBATCH --nodes 4
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH -o p32salt.%j.%N.out
#SBATCH -e p32salt.%j.%N.error

echo wait .... 30 seconds
sleep 30
echo date=`date`


# Salt

srun ./parallel_qoz salt.config 51 368 1010 1010
#printf "QOZ END"


