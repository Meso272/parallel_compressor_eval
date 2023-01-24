#!/bin/bash
#SBATCH --job-name=p256ara
#SBATCH -p bdwall
#SBATCH -A ECP-EZ
#SBATCH --nodes 8
#SBATCH --ntasks-per-node=32
#SBATCH --time=8:00:00
#SBATCH -o p256ara.%j.%N.out
#SBATCH -e p256ara.%j.%N.error

echo wait .... 30 seconds
sleep 30
echo date=`date`




printf "ARAMCO\n"
srun ./parallel_sz3 sz3.config 4e-5 50 235 449 449
printf "SZ3 END\n"
srun ./parallel_qoz qoz.config 5e-5 50 235 449 449
printf "QOZ END\n"
srun ./parallel_qoz  sperr.config 8e-5 50 235 449 449
printf "SPERR END\n"
srun ./parallel_qoz  fz.config 2.2e-4 50 235 449 449
printf "FZ END\n"


