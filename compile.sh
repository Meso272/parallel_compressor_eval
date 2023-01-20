#! /bin/bash

module purge
module load gcc/8.2.0-g7hppkz
module load mpich
module load zstd/1.4.5-uxapdkl
#module load python/intel-parallel-studio-cluster.2019.5-zqvneip/3.6.9
rm parallel_zfp parallel_sz2 parallel_sz3 parallel_mgard parallel_qoz

szsrc=$HOME/packages/SZ2/include
szlib=$HOME/packages/SZ2/lib/
sz3src=$HOME/packages/SZ3/include
sz3lib=$HOME/packages/SZ3/lib64/
zfpsrc=$HOME/zfp/include
zfplib=$HOME/zfp/lib


#mgardlib=$HOME/code/mgardp/build/lib/
mgardsrc=$HOME/MGARDx/include
zstdsrc=$HOME/packages/zstd/include
zstdlib=$HOME/packages/zstd/lib
metasrc=$HOME/packages/meta_compressor/include/sz_cpp   
metalib=$HOME/packages/meta_compressor/lib64
qozsrc=$HOME/packages/QoZ/include
qozlib=$HOME/packages/QoZ/lib64/
pythonsrc=$HOME/miniconda3/pkgs/python-3.8.5-h7579374_1/include/python3.8
pythonlib=$HOME/miniconda3/pkgs/python-3.8.5-h7579374_1/lib
pybind11src=$HOME/miniconda3/pkgs/pybind11-2.6.1-py38h82cb98a_0/include
pybind11lib=$HOME/miniconda3/pkgs/pybind11-2.6.1-py38h82cb98a_0/lib
gcc -c rw.c
g++ -c rw.c -o rwx.o

#mpicc -std=c99 -O3 rw.o parallel_sz2.c -o parallel_sz2 -I$szsrc $szlib/libSZ.a $szlib/libzstd.a $szlib/libzlib.a -lm

mpicc -std=c99 -O3 rw.o parallel_sz2.c -o parallel_sz2 -I $szsrc -L $szlib -l SZ -l zstd -lm

mpicxx  -O3 rwx.o parallel_sz3.c -o parallel_sz3 -I $sz3src -L $sz3lib  -l zstd -lm

mpicxx  -std=c++17 -O3 rwx.o parallel_qoz.c -o parallel_qoz -I $qozsrc -L $qozlib -I $pythonsrc -L $pythonlib  -l python3.8  -I $pybind11src -L $pybind11lib  -l zstd -lm 




mpicc -std=c99 -O3 rw.o parallel_zfp.c -o parallel_zfp -I $zfpsrc -L $zfplib -l zfp -lm

mpicxx  -O3 rwx.o parallel_mgard.c -o parallel_mgard -I $mgardsrc  -I $metasrc -L $metalib -l sz_cpp -l zstd -lm