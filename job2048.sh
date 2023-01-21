#!/bin/bash
#SBATCH --job-name=p2048
#SBATCH -p bdwall
#SBATCH -A ECP-EZ
#SBATCH --nodes 64
#SBATCH --ntasks-per-node=32
#SBATCH --time=4:00:00
#SBATCH -o p2048.%j.%N.out
#SBATCH -e p2048.%j.%N.error

echo wait .... 30 seconds
sleep 30
echo date=`date`


# Hurricane
#srun ./parallel_sz 13 496 496 96
#srun ./parallel_zfp 13 496 496 96
#srun ./parallel_mgard 13 496 496 96
printf "HurricaneG1\n"
srun ./parallel_sz3 sz3.config 7e-4 13 500 500 100
printf "SZ3 END\n"
srun ./parallel_qoz  qoz.config 8e-4 13 500 500 100
printf "QOZ END\n"
srun ./parallel_qoz  sperr.config 1e-3 13 500 500 100
printf "SPERR END\n"
srun ./parallel_qoz  fz.config 1e-3 13 500 500 100
printf "FZ END\n"

printf "HurricaneG2\n"
srun ./parallel_sz3 sz3.config 3e-4 13 500 500 100
printf "SZ3 END\n"
srun ./parallel_qoz  qoz.config 3e-4 13 500 500 100
printf "QOZ END\n"
srun ./parallel_qoz  sperr.config 4e-4 13 500 500 100
printf "SPERR END\n"
srun ./parallel_qoz  fz.config 4e-4 13 500 500 100
printf "FZ END\n"

# miranda
#srun ./parallel_sz  7 384 384 256
#srun ./parallel_zfp 7 384 384 256
#srun ./parallel_mgard 7 384 384 256
#printf "Miranda\n"
#srun ./parallel_sz3 sz3.config 7 384 384 256
#printf "SZ3 END\n"
#srun ./parallel_qoz qoz.config 1e-3 7 384 384 256
#printf "QOZ END\n"
#run ./parallel_qoz  sperr.config 1e-3  7 384 384 256
#printf "SPERR END\n"
#srun ./parallel_qoz  fz.config 1e-3 7 384 384 256
#printf "FZ END\n"

# scale
#printf "SCALE\n"
#srun ./parallel_sz3 sz3.config 12 1200 1200 98
#printf "SZ3 END\n"
#srun ./parallel_qoz qoz.config 1e-3 12 1200 1200 98
#printf "QOZ END\n"
#srun ./parallel_qoz  sperr.config 1e-3  12 1200 1200 98
#printf "SPERR END\n"
#srun ./parallel_qoz  fz.config 1e-3 12 1200 1200 98
#printf "FZ END\n"

#printf "ARAMCO\n"
#srun ./parallel_sz3 sz3.config 60 235 449 449
#printf "SZ3 END\n"
#srun ./parallel_qoz qoz_b16.config 1e-3 60 235 449 449
#printf "QOZ END\n"
#srun ./parallel_qoz  sperr.config 1e-3  60 235 449 449
#printf "SPERR END\n"
#srun ./parallel_qoz  fz.config 1e-3 60 235 449 449
#printf "FZ END\n"


