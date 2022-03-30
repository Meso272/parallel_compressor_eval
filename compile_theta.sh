#! /bin/bash

#module purge
#module load gcc/8.2.0-g7hppkz
#module load mpich
#module load zstd/1.4.5-uxapdkl
#rm parallel_zfp parallel_sz2 parallel_sz3 parallel_mgard parallel_qoz

szsrc=$HOME/packages/SZ/include
szlib=$HOME/packages/SZ/lib/
sz3src=$HOME/packages/SZ3/include
sz3lib=$HOME/packages/SZ3/lib64/
zfpsrc=$HOME/lossycompression/zfp/include
zfplib=$HOME/lossycompression/zfp/lib


#mgardlib=$HOME/code/mgardp/build/lib/
mgardsrc=$HOME/lossycompression/MGARDx/include
zstdsrc=$HOME/packages/zstd/include
zstdlib=$HOME/packages/zstd/lib
metasrc=$HOME/packages/SZcpp/include/sz_cpp   
metalib=$HOME/packages/SZcpp/lib64
qozsrc=$HOME/packages/QOZ/include
qozlib=$HOME/packages/QOZ/lib64/
gcc -c rw.c
g++ -c rw.c -o rwx.o

# mpicc -std=c99 -O3 rw.o parallel_sz2.c -o parallel_sz2 -I$szsrc $szlib/libSZ.a $szlib/libzstd.a $szlib/libzlib.a -lm

#cc -std=c99 -O3 parallel_sz2_theta.c -o parallel_sz2 -I $szsrc -L $szlib  -l SZ -I $zstdsrc -L $zstdlib -l zlib -l zstd -lm

#CC  -O3 rwx.o parallel_sz3_theta.c -o parallel_sz3 -I $sz3src -L $sz3lib  -I $zstdsrc -L $zstdlib -l zstd -lm

#CC  -O3 rwx.o parallel_sz3_theta.c -o parallel_qoz -I $qozsrc -L $qozlib -I $zstdsrc -L $zstdlib -l zstd -lm 




#cc -std=c99 -O3 rw.o parallel_zfp_theta.c -o parallel_zfp -I $zfpsrc -L $zfplib -l zfp -lm

CC  -O3 rwx.o parallel_mgard_theta.c -o parallel_mgard -I $mgardsrc  -I $metasrc -L /home/jliu447/packages/SZcpp/lib64 -l sz_cpp -I $zstdsrc -L $zstdlib  -l zstd -lm