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
#include "zfp.h"

unsigned char* zfp_compress_3D(float* array, double tolerance, size_t r1, size_t r2, size_t r3, size_t* out_size) {
    int status = 0;    /* return value: 0 = success */
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_float;
    field = zfp_field_3d(array, type, r3, r2, r1);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    buffer = malloc(bufsize);

    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    zfpsize = zfp_compress(zfp, field);
    if (!zfpsize) {
        fprintf(stderr, "compression failed\n");
        status = 1;
    }

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    *out_size = zfpsize;
    return (unsigned char*)buffer;
}

float* zfp_decompress_3D(unsigned char* comp_data, double tolerance, size_t buffer_size, size_t r1, size_t r2, size_t r3) {
    int status = 0;    /* return value: 0 = success */
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    float* array = (float*)malloc(r1 * r2 * r3 * sizeof(float));
    type = zfp_type_float;
    field = zfp_field_3d(array, type, r3, r2, r1);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    // buffer = malloc(bufsize);
    buffer = (void*)comp_data;
    bufsize = buffer_size;

    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    if (!zfp_decompress(zfp, field)) {
        fprintf(stderr, "decompression failed\n");
        status = 1;
    }
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    return array;
}


int main(int argc, char* argv[]) {
    srand(time(0));
    char* cfgFile;

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   
    double * absbound;

    //Begin modify this part

    char* file_folder = "/lcrc/project/ECP-EZ/public/compression/datasets/0";
    int num_vars = atoi(argv[1]);

    //Begin modify this part
    int hurricane_num_vars = 13;
    char hurricane_file[13][50] = {"Uf48_truncated.bin.dat", "Vf48_truncated.bin.dat", "Wf48_truncated.bin.dat",
                                   "TCf48_truncated.bin.dat", "Pf48_truncated.bin.dat", "QVAPORf48_truncated.bin.dat",
                                   "CLOUDf48_log10_truncated.bin.dat", "QCLOUDf48_log10_truncated.bin.dat", "QICEf48_log10_truncated.bin.dat",
                                   "QRAINf48_log10_truncated.bin.dat", "QSNOWf48_log10_truncated.bin.dat", "QGRAUPf48_log10_truncated.bin.dat",
                                   "PRECIPf48_log10_truncated.bin.dat"};
    double hurricane_abs_bound[13] = {1.851616, 1.874008, 0.331488, 2.12402, 132.72278, 0.0004084144, 0.20304, 0.7021114, 0.373662, 0.714628, 0.695594, 0.715312, 0.716028};
    // miranda
    int miranda_num_vars = 7;
    char miranda_file[7][50] = {"velocityy_truncated.bin.dat", "velocityx_truncated.bin.dat", "density_truncated.bin.dat",
                                "pressure_truncated.bin.dat", "velocityz_truncated.bin.dat", "viscocity_truncated.bin.dat",
                                "diffusivity_truncated.bin.dat"};
    double miranda_abs_bound[7] = {0.1317327, 0.135345, 0.06, 0.1324284, 0.4246152, 0.0786399, 0.07845450000000001};

    int scale_num_vars = 12;
    char scale_file[12][50] = {"PRES-98x1200x1200.dat", "QC-98x1200x1200.log10.dat", "QG-98x1200x1200.log10.dat",
                                   "QI-98x1200x1200.log10.dat", "QR-98x1200x1200.log10.dat", "QS-98x1200x1200.log10.dat",
                                   "QV-98x1200x1200.log10.dat", "RH-98x1200x1200.dat", "T-98x1200x1200.dat",
                                   "U-98x1200x1200.dat", "V-98x1200x1200.dat", "W-98x1200x1200.dat",
                                  };
    double scale_abs_bound[12] ={1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3};

    int aramco_num_vars = 60;
    char aramco_file[60][50];
    double aramco_abs_bound[60];
    for (int i=0;i<60;i++){
        //char name[50];
        sprintf(aramco_file[i],"aramco-snapshot-%d.f32",1000+10*i);
        //salt_file[i]=name;
        aramco_abs_bound[i]=1e-3;

    }





    char file[100][50];
    
    
    if (num_vars == hurricane_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], hurricane_file[i]);
        absbound = hurricane_abs_bound;
    } else if (num_vars == miranda_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], miranda_file[i]);
        absbound = miranda_abs_bound;
    } 
    else if (num_vars == scale_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], scale_file[i]);
        absbound = scale_abs_bound;
    }
    else if (num_vars == aramco_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], aramco_file[i]);
        absbound = aramco_abs_bound;
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
    //End modify this part

    if (world_rank == 0) printf("Start parallel compressing ... \n");
    if (world_rank == 0) printf("size: %d\n", world_size);
    double start, end;
    double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);

    size_t compressed_size[100];
    char zip_filename[100];
    // char out_filename[100];
    size_t inSize, outSize;
    size_t nbEle;
    int status;
    float* dataIn;

    size_t est_compressed_size = r1 * r2 * r3 * sizeof(float) * num_vars / 3;
    unsigned char* compressed_output = (unsigned char*)malloc(est_compressed_size);
    unsigned char* compressed_output_pos = compressed_output;
    int folder_index = world_rank;
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
            dataIn = (float*)malloc(nbEle * sizeof(float));
            MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) {
            end = MPI_Wtime();
            costReadOri += end - start;
        }

        // float max = dataIn[0];
        // float min = dataIn[0];
        // for (int j = 1; j < r1 * r2 * r3; j++) {
        //     if (max < dataIn[j]) max = dataIn[j];
        //     if (min > dataIn[j]) min = dataIn[j];
        // }
        // absbound[i] = rel_bound[i] * (max - min);

        // Compress Input Data
        size_t out_size;
        // if (world_rank == 0) printf ("Compressing %s\n", filename);
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) start = MPI_Wtime();
        unsigned char* bytesOut = zfp_compress_3D(dataIn, absbound[i], r3, r2, r1, &compressed_size[i]);
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
    sprintf(zip_filename, "%s/zfp_%d_%d.out", "/lcrc/globalscratch/jinyang", folder_index, rand());  // Write Compressed Data

    size_t total_size = compressed_output_pos - compressed_output;
    // Write Compressed Data
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) start = MPI_Wtime();
    // if (world_rank == 0) printf("write compressed file to disk %s \n", zip_filename);
    writeByteData(compressed_output, total_size, zip_filename, &status);
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {
        end = MPI_Wtime();
        costWriteZip += end - start;
    }
    free(compressed_output);
    // Read Compressed Data
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) start = MPI_Wtime();
    // if (world_rank == 0) printf("read compressed file from disk %s \n", zip_filename);
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
        float* dataOut = zfp_decompress_3D(compressed_output_pos, absbound[i], compressed_size[i], r3, r2, r1);
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
        printf("ZFP Finish parallel compressing, total compression ratio = %.4g\n", 1.0 * r1 * r2 * r3 * sizeof(float) * num_vars / total_size);
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

    MPI_Finalize();

    return 0;
}
