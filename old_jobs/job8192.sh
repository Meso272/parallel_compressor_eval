#!/bin/bash
#SBATCH --job-name=p8192
#SBATCH -p bdwall
#SBATCH -A ECP-EZ
#SBATCH --nodes 256
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:20:00
#SBATCH -o p8192.%j.%N.out
#SBATCH -e p8192.%j.%N.error

echo wait .... 30 seconds
sleep 30
echo date=`date`


# Hurricane
srun ./parallel_sz 13 496 496 96
srun ./parallel_zfp 13 496 496 96
srun ./parallel_mgard 13 496 496 96
srun ./parallel_sz3 hurricane.config 13 496 496 96
printf "SZ3 END\n"
srun ./parallel_qoz hurricane.config 13 496 496 96
printf "QOZ END\n"

# miranda
#srun ./parallel_sz  7 384 384 256
#srun ./parallel_zfp 7 384 384 256
#srun ./parallel_mgard 7 384 384 256
#srun ./parallel_sz3 miranda.config 7 384 384 256
#printf "SZ3 END"
#srun ./parallel_qoz miranda.config 7 384 384 256
#printf "QOZ END"

