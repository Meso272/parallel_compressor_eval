/**
 *  @file test_compress_parallel.c
 *  @author Dingwen Tao
 *  @date January, 2017
 *  @brief This is an example of using compression interface in parallel
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"
#include "rw.h"
#include "sz.h"

int main(int argc, char *argv[]) {
    srand(time(0));

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    

    if (argc < 2) {
        printf("Test case: parallel_sz num_vars r1 r2 r3\n");
        printf("Example: parallel_sz 7 384 384 256 \n");
        exit(0);
    }
    char *file_folder = "/lcrc/project/ECP-EZ/public/compression/datasets/0";
    int num_vars = atoi(argv[1]);

    //Begin modify this part
    int hurricane_num_vars = 13;
    char hurricane_file[13][50] = {"Uf48_truncated.bin.dat", "Vf48_truncated.bin.dat", "Wf48_truncated.bin.dat",
                                   "TCf48_truncated.bin.dat", "Pf48_truncated.bin.dat", "QVAPORf48_truncated.bin.dat",
                                   "CLOUDf48_log10_truncated.bin.dat", "QCLOUDf48_log10_truncated.bin.dat", "QICEf48_log10_truncated.bin.dat",
                                   "QRAINf48_log10_truncated.bin.dat", "QSNOWf48_log10_truncated.bin.dat", "QGRAUPf48_log10_truncated.bin.dat",
                                   "PRECIPf48_log10_truncated.bin.dat"};
    double hurricane_rel_bound[13] = {9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4, 9e-4};
    // miranda
    int miranda_num_vars = 7;
    char miranda_file[7][50] = {"velocityy_truncated.bin.dat", "velocityx_truncated.bin.dat", "density_truncated.bin.dat",
                                "pressure_truncated.bin.dat", "velocityz_truncated.bin.dat", "viscocity_truncated.bin.dat",
                                "diffusivity_truncated.bin.dat"};
    double miranda_rel_bound[7] = {1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3};

    int aramco_num_vars = 50;
    char aramco_file[50][50];
    double aramco_rel_bound[50];
    for (int i=0;i<50;i++){
        //char name[50];
        sprintf(aramco_file[i],"aramco-snapshot-%d.f32",1000+10*i);
        //salt_file[i]=name;
        aramco_rel_bound[i]=1e-3;

    }



    char file[100][50];
    double *rel_bound;
    
    if (num_vars == hurricane_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], hurricane_file[i]);
        rel_bound = hurricane_rel_bound;
    } else if (num_vars == miranda_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], miranda_file[i]);
        rel_bound = miranda_rel_bound;
    } 
    else if (num_vars == aramco_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], aramco_file[i]);
        rel_bound = aramco_rel_bound;
    }
    else {
        printf("No such variable, exit\n");
//        SZ_Finalize();
        MPI_Finalize();
        return 0;
    }


    size_t  r1,r2,r3,r4,r5; //384 384 256, 500 500 100
    r1 = atof(argv[2]);
    r2 = atof(argv[3]);
    r3= atof(argv[4]);
    r4=0;
    r5=0;
   
    //End modify this part

    SZ_Init(NULL);

    if (world_rank == 0) printf("Start parallel compressing ... \n");
    if (world_rank == 0) printf("world_size: %d\n", world_size);
    double start, end;
    double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);

    size_t compressed_size[100];
    char zip_filename[100];
    // char out_filename[100];
    size_t inSize, outSize;
    size_t nbEle;
    int status;
    float *dataIn;

    size_t est_compressed_size = r1 * r2 * r3 * sizeof(float) * num_vars / 5;
    unsigned char *compressed_output = (unsigned char *)malloc(est_compressed_size);
    unsigned char *compressed_output_pos = compressed_output;
    for (int i = 0; i < num_vars; i++) {
        char filename[100];
        sprintf(filename, "%s/%s", file_folder, file[i]);
        // Read Input Data
        if (world_rank == 0) {
            start = MPI_Wtime();
            dataIn = readFloatData(filename, &nbEle, &status);
            end = MPI_Wtime();
            // printf("file %s read time: %.2f\n", filename, end - start);
            start = MPI_Wtime();
            MPI_Bcast(&nbEle, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
            end = MPI_Wtime();
            // printf("broadcast time: %.2f\n", end - start);
        } else {
            MPI_Bcast(&nbEle, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            dataIn = (float *)malloc(nbEle * sizeof(float));
            MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) {
            end = MPI_Wtime();
            costReadOri += end - start;
        }

        // Compress Input Data
        size_t out_size;
        // if (world_rank == 0) printf("Compressing %s\n", filename);
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) start = MPI_Wtime();
        unsigned char *bytesOut = SZ_compress_args(SZ_FLOAT, dataIn, &compressed_size[i], REL, 0, rel_bound[i], 0, r5, r4, r3, r2, r1);
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) {
            end = MPI_Wtime();
            costComp += end - start;
        }
        free(dataIn);
        memcpy(compressed_output_pos, bytesOut, compressed_size[i]);
        compressed_output_pos += compressed_size[i];
        free(bytesOut);
    }

    struct stat st = {0};
    if (stat("/lcrc/globalscratch/jinyang", &st) == -1) {
        mkdir("/lcrc/globalscratch/jinyang", 0777);
    }
    sprintf(zip_filename, "%s/sz_%d_%d_%d.out", "/lcrc/globalscratch/jinyang", world_rank, world_size, rand());  // Write Compressed Data
    size_t total_size = compressed_output_pos - compressed_output;
    // Write Compressed Data
    MPI_Barrier(MPI_COMM_WORLD);
    // if (world_rank == 0) printf("write compressed file to disk %s \n", zip_filename);
    if (world_rank == 0) start = MPI_Wtime();
    writeByteData(compressed_output, total_size, zip_filename, &status);
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {
        end = MPI_Wtime();
        costWriteZip += end - start;
    }
    free(compressed_output);

    // Read Compressed Data
    MPI_Barrier(MPI_COMM_WORLD);
    // if (world_rank == 0) printf("read compressed file from disk %s \n", zip_filename);
    if (world_rank == 0) start = MPI_Wtime();
    compressed_output = readByteData(zip_filename, &inSize, &status);
    if (inSize != total_size) {
        printf("ERROR! Broken file : %s", zip_filename);
    } else {
        remove(zip_filename);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        end = MPI_Wtime();
        costReadZip += end - start;
    }
    compressed_output_pos = compressed_output;

    for (int i = 0; i < num_vars; i++) {
        // Decompress Compressed Data
        MPI_Barrier(MPI_COMM_WORLD);
        // if (world_rank == 0) printf("decompress %d-th field\n", i);
        if (world_rank == 0) start = MPI_Wtime();
        float *dataOut = SZ_decompress(SZ_FLOAT, compressed_output_pos, compressed_size[i], r5, r4, r3, r2, r1);
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) {
            end = MPI_Wtime();
            costDecomp += end - start;
        }
        compressed_output_pos += compressed_size[i];
        free(dataOut);
    }
    free(compressed_output);

    if (world_rank == 0) {
        printf("SZ2 Finish parallel compressing, total compression ratio = %.4g\n", 1.0 * r1 * r2 * r3 * sizeof(float) * num_vars / total_size);
        printf("Separate ratios\n");
        for (int i = 0; i < num_vars; i++) {
            printf("%s = %.4g \n", file[i], 1.0 * r1 * r2 * r3 * sizeof(float) / compressed_size[i]);
        }
        // printf("\n");
        // printf("Timecost of reading original files = %.2f seconds\n", costReadOri);
        printf("Timecost of compressing using %d processes = %.4f seconds\n", world_size, costComp);
        printf("Timecost of writing compressed files = %.4f seconds\n", costWriteZip);
        printf("Timecost of reading compressed files = %.4f seconds\n", costReadZip);
        printf("Timecost of decompressing using %d processes = %.4f seconds\n", world_size, costDecomp);
        // printf("Timecost of writing decompressed files = %.2f seconds\n", costWriteOut);
        printf("Throughput of reading compressed files = %.4f GB/s\n", total_size * world_size / 1000.0 / 1000 / 1000 / costReadZip);
        printf("Throughput of writing compressed files = %.4f GB/s\n", total_size * world_size / 1000.0 / 1000 / 1000 / costWriteZip);
    }

    SZ_Finalize();

    MPI_Finalize();

    return 0;
}
