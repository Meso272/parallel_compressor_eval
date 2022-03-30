#! /bin/bash

#module purge
#module load gcc/8.2.0-g7hppkz
#module load mpich
#module load zstd/1.4.5-uxapdkl
rm parallel_zfp parallel_sz2 parallel_sz3 parallel_mgard parallel_qoz

szsrc=$HOME/packages/SZ2/include
szlib=$HOME/packages/SZ2/lib/
sz3src=$HOME/packages/SZ3/include
sz3lib=$HOME/packages/SZ3/lib64/
zfpsrc=$HOME/lossycompression/zfp/include
zfplib=$HOME/lossycompression/zfp/lib


#mgardlib=$HOME/code/mgardp/build/lib/
mgardsrc=$HOME/lossycompression/MGARDx/include
zstdsrc=$HOME/packages/zstd/include
zstdlib=$HOME/packages/zstd/lib
metasrc=$HOME/packages/meta_compressor/include/sz_cpp   
metalib=$HOME/packages/meta_compressor/lib64
qozsrc=$HOME/packages/QOZ/include
qozlib=$HOME/packages/QOZ/lib64/
gcc -c rw.c
g++ -c rw.c -o rwx.o

# mpicc -std=c99 -O3 rw.o parallel_sz2.c -o parallel_sz2 -I$szsrc $szlib/libSZ.a $szlib/libzstd.a $szlib/libzlib.a -lm

mpicc -std=c99 -O3 rw.o parallel_sz2_theta.c -o parallel_sz2 -I $szsrc -L $szlib -l SZ -l zstd -lm

mpicxx  -O3 rwx.o parallel_sz3_theta.c -o parallel_sz3 -I $sz3src -L $sz3lib  -l zstd -lm

mpicxx  -O3 rwx.o parallel_sz3_theta.c -o parallel_qoz -I $qozsrc -L $qozlib -l zstd -lm 




mpicc -std=c99 -O3 rw.o parallel_zfp_theta.c -o parallel_zfp -I $zfpsrc -L $zfplib -l zfp -lm

mpicxx  -O3 rwx.o parallel_mgard_theta.c -o parallel_mgard -I $mgardsrc  -I $metasrc -L $metalib -l sz_cpp -l zstd -lm