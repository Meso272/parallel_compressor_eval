#! /bin/bash
echo wait .... 30 seconds
sleep 30
echo date=`date`
cd $HOME/lossycompression/parallel_compressor_eval

# Hurricane
aprun -n 16384 ./parallel_sz2 13 496 496 96
aprun -n 16384 ./parallel_zfp 13 496 496 96
#aprun -n 16384 ./parallel_mgard 13 496 496 96
aprun -n 16384 ./parallel_sz3 hurricane.config 13 496 496 96
printf "SZ3 END"
srun -n 16384 ./parallel_qoz hurricane.config 13 496 496 96
printf "QOZ END"

# miranda
aprun -n 16384 ./parallel_sz2  7 384 384 256
aprun -n 16384 ./parallel_zfp 7 384 384 256
#srun ./parallel_mgard 7 384 384 256
aprun -n 16384 ./parallel_sz3 miranda.config 7 384 384 256
printf "SZ3 END"
aprun -n 16384 ./parallel_qoz miranda.config 7 384 384 256
printf "QOZ END"

