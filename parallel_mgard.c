/**
 *  @file parallel_kai.c
 *  @author Kai Zhao
 *  @date January, 2020
 *  @brief This is an example of using compression interface in parallel
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "decompose.hpp"
#include "recompose.hpp"
#include "rw.h"
#include "mpi.h"

// USAGE
// mpirun -np 16 parallel sz.config folder_num r3 r2 r1
int main(int argc, char * argv[])
{
    srand(time(0));
	size_t r5=0,r4=0,r3=0,r2=0,r1=0;
	char *cfgFile;

	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);	
	
	if(argc < 3)
	{
		printf("Test case: parallel_mgard num_vars [dimension sizes...]\n");
        printf("Example: parallel_mgard 7 384 384 256\n");
		exit(0);
	}

	//cfgFile=argv[1];
	
	if(argc>=3)
	  r1 = atoi(argv[2]); //8
	if(argc>=4)
	  r2 = atoi(argv[3]); //8
	if(argc>=5)
	  r3 = atoi(argv[4]); //128
	if(argc>=6)
	  r4 = atoi(argv[5]);
	if(argc>=7)
	  r5 = atoi(argv[6]);
	
//	SZ_Init(NULL);

	if (world_rank == 0) printf ("Start parallel compressing ... \n");
	if (world_rank == 0) printf("size: %d\n", world_size);
	double start, end;
	double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);
    int num_vars = atoi(argv[1]);

    int qmcpack8h_num_vars = 2;
    char qmcpack8h_file[2][50] = {"spin_0_truncated.bin.dat", "spin_1_truncated.bin.dat"};
    double qmcpack8h_rel_bound[2] = {1e-6, 1e-6};

    // qmacpack6k
    int qmcpack6k_num_vars = 20;
    char qmacpack6k_file[20][50] = {"s2700l300_truncated.bin.dat", "s4500l300_truncated.bin.dat", "s1200l300_truncated.bin.dat",
                                    "s300l300_truncated.bin.dat", "s4200l300_truncated.bin.dat", "s5400l300_truncated.bin.dat",
                                    "s1800l300_truncated.bin.dat", "s5700l300_truncated.bin.dat", "s4800l300_truncated.bin.dat",
                                    "s3300l300_truncated.bin.dat", "s5100l300_truncated.bin.dat", "s1500l300_truncated.bin.dat",
                                    "s600l300_truncated.bin.dat", "s0l300_truncated.bin.dat", "s3600l300_truncated.bin.dat",
                                    "s900l300_truncated.bin.dat", "s3900l300_truncated.bin.dat", "s3000l300_truncated.bin.dat",
                                    "s2100l300_truncated.bin.dat", "s2400l300_truncated.bin.dat"};
    double qmacpack6k_rel_bound[20] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
                                       1e-6, 1e-6, 1e-6, 1e-6, 1e-6};

    // Hurricane
    int hurricane_num_vars = 13;
    char hurricane_file[13][50] = {"Uf48_truncated.bin.dat", "Vf48_truncated.bin.dat", "Wf48_truncated.bin.dat",
                                   "TCf48_truncated.bin.dat", "Pf48_truncated.bin.dat", "QVAPORf48_truncated.bin.dat",
                                   "CLOUDf48_log10_truncated.bin.dat", "QCLOUDf48_log10_truncated.bin.dat", "QICEf48_log10_truncated.bin.dat",
                                   "QRAINf48_log10_truncated.bin.dat", "QSNOWf48_log10_truncated.bin.dat", "QGRAUPf48_log10_truncated.bin.dat",
                                   "PRECIPf48_log10_truncated.bin.dat"};
    double hurricane_rel_bound[13] = {4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3};
    // miranda
    int miranda_num_vars = 7;
    char miranda_file[7][50] = {"velocityy_truncated.bin.dat", "velocityx_truncated.bin.dat", "density_truncated.bin.dat",
                                "pressure_truncated.bin.dat", "velocityz_truncated.bin.dat", "viscocity_truncated.bin.dat",
                                "diffusivity_truncated.bin.dat"};
    double miranda_rel_bound[7] = {0.02,0.02,0.02,0.02,0.02,0.02,0.02};

    // assignment
    char file[20][50];
    double *rel_bound;
    if (num_vars == qmcpack6k_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], qmacpack6k_file[i]);
        rel_bound = qmacpack6k_rel_bound;
    } else if (num_vars == qmcpack8h_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], qmcpack8h_file[i]);
        rel_bound = qmcpack8h_rel_bound;
    } else if (num_vars == hurricane_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], hurricane_file[i]);
        rel_bound = hurricane_rel_bound;
    } else if (num_vars == miranda_num_vars) {
        for (int i = 0; i < num_vars; i++) strcpy(file[i], miranda_file[i]);
        rel_bound = miranda_rel_bound;
    } else {
        printf("No such variable, exit\n");
//        SZ_Finalize();
        MPI_Finalize();
        return 0;
    }
	size_t compressed_size[20];

	char folder[50] = "/lcrc/project/ECP-EZ/public/compression/datasets";
	char filename[100];
	char zip_filename[100];
	// char out_filename[100];
	size_t inSize, outSize; 
	size_t nbEle;
	int status;
	float * dataIn;

	size_t est_compressed_size = r1 * r2 * r3 * sizeof(float) * num_vars / 3;
	unsigned char * compressed_output = (unsigned char *) malloc(est_compressed_size);
	unsigned char * compressed_output_pos = compressed_output;
	int folder_index = world_rank;
    
    std::vector<size_t> dims={r3,r2,r1};

	for(int i=0; i<num_vars; i++){
		sprintf(filename, "%s/%d/%s", folder, folder_index, file[i]);
        
		// Read Input Data
		if(world_rank == 0){
			start = MPI_Wtime();
			dataIn = readFloatData(filename, &nbEle, &status);
			end = MPI_Wtime();
			//printf("file %s read time: %.2f\n", filename, end - start);
			start = MPI_Wtime();
			MPI_Bcast(&nbEle, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
			end = MPI_Wtime();
			//printf("broadcast time: %.2f\n", end - start);
		}
		else{
			MPI_Bcast(&nbEle, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
			dataIn = (float *) malloc(nbEle * sizeof(float));
			MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0){
			end = MPI_Wtime();
			costReadOri += end - start;
		}
		
		// Compress Input Data
		size_t out_size;
		//if (world_rank == 0) printf ("Compressing %s\n", filename);
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0) start = MPI_Wtime();
        MGARD::Decomposer<float> decomposer(1);

        unsigned char *bytesOut = decomposer.compress(dataIn, dims,3,rel_bound[i],compressed_size[i]);
//		unsigned char *bytesOut = SZ_compress_args(SZ_FLOAT, dataIn, &compressed_size[i], REL, 0, rel_bound[i], 0, r5, r4, r3, r2, r1);
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0){
			end = MPI_Wtime();
			costComp += end - start;
		}
		free (dataIn);
		memcpy(compressed_output_pos, bytesOut, compressed_size[i]);
		compressed_output_pos += compressed_size[i];
		free(bytesOut);
	}
    struct stat st = {0};
    if (stat("/lcrc/globalscratch/jinyang", &st) == -1) {
        mkdir("/lcrc/globalscratch/jinyang", 0777);
    }
    sprintf(zip_filename, "%s/mgx_%d_%d.out", "/lcrc/globalscratch/jinyang", folder_index, rand());	// Write Compressed Data
    size_t total_size = compressed_output_pos - compressed_output;
	// Write Compressed Data
	MPI_Barrier(MPI_COMM_WORLD);
    //if (world_rank == 0) printf("write compressed file to disk %s \n", zip_filename);
    if(world_rank == 0) start = MPI_Wtime();
	writeByteData(compressed_output, total_size, zip_filename, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank == 0){
		end = MPI_Wtime();
		costWriteZip += end - start;
	}
	free(compressed_output);
	// Read Compressed Data
    MPI_Barrier(MPI_COMM_WORLD);
    //if (world_rank == 0) printf("read compressed file from disk %s \n", zip_filename);
    if(world_rank == 0) start = MPI_Wtime();
	compressed_output = readByteData(zip_filename, &inSize, &status);
    if (inSize != total_size) {
        printf("ERROR! Broken file : %s", zip_filename);
    } else {
        remove(zip_filename);
    }
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank == 0){
		end = MPI_Wtime();
		costReadZip += end - start;
	}
	compressed_output_pos = compressed_output;



    for(int i=0; i<num_vars; i++){
		// Decompress Compressed Data
        MPI_Barrier(MPI_COMM_WORLD);
        //if (world_rank == 0) printf("decompress %d-th field\n", i);
        if(world_rank == 0) start = MPI_Wtime();
        MGARD::Recomposer<float> recomposer;
        float *dataOut = recomposer.decompress(compressed_output_pos, compressed_size[i],dims);
//        float *dataOut = SZ_decomprescs(SZ_FLOAT, compressed_output_pos, compressed_size[i], r5, r4, r3, r2, r1);
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank == 0){
			end = MPI_Wtime();
			costDecomp += end - start; 
		}
		compressed_output_pos += compressed_size[i];
		free(dataOut);
	}
	free(compressed_output);

	if (world_rank == 0)
	{
		printf ("MGARD+ Finish parallel compressing, total compression ratio %.4g.\n", 1.0*r1*r2*r3*sizeof(float)*num_vars / total_size);
		printf("Separate ratios: ");
		for(int i=0; i<num_vars; i++){
			printf("%.4g ", 1.0*r1*r2*r3*sizeof(float) / compressed_size[i]);
		}
		printf("\n");
		printf ("Timecost of reading original files = %.4f seconds\n", costReadOri);
        printf ("Timecost of compressing using %d processes = %.4f seconds\n", world_size, costComp);
        printf ("Timecost of writing compressed files = %.4f seconds\n", costWriteZip);
        printf ("Timecost of reading compressed files = %.4f seconds\n", costReadZip);
        printf ("Timecost of decompressing using %d processes = %.4f seconds\n\n", world_size, costDecomp);
        printf ("Timecost of writing decompressed files = %.4f seconds\n", costWriteOut);
	}


//	SZ_Finalize();

	MPI_Finalize();

	return 0;
}
