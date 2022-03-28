#! /bin/bash

module purge
module load gcc/8.2.0-g7hppkz
module load mpich
rm parallel_zfp parallel_sz2

szsrc=$HOME/packages/SZ2/include
szlib=$HOME/packages/SZ2/lib/
sz3src=$HOME/packages/SZ3/include/SZ3/api
sz3lib=$HOME/packages/SZ3/lib64/
zfpsrc=$HOME/zfp/include
zfplib=$HOME/zfp/lib
mgardsrc=$HOME/MGARDx/include
#mgardlib=$HOME/code/mgardp/build/lib/


gcc -c rw.c

# mpicc -std=c99 -O3 rw.o parallel_sz2.c -o parallel_sz2 -I$szsrc $szlib/libSZ.a $szlib/libzstd.a $szlib/libzlib.a -lm

mpicc -std=c99 -O3 rw.o parallel_sz2.c -o parallel_sz2 -I$szsrc -L$szlib -lSZ -lzstd -lm

mpicc -std=c99 -O3 rw.o parallel_sz3.c -o parallel_sz3 -I$sz3src -L$sz3lib -lSZ -lzstd -lm

mpicc -std=c99 -O3 rw.o parallel_zfp.c -o parallel_zfp -I$zfpsrc -L$zfplib -lzfp -lm

mpicc -std=c99 -O3 rw.o parallel_mgard.c -o parallel_mgard -I$mgardsrc   -lm