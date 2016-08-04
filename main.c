#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "blackhole_lab.h"

// Include MPI Header
#include "mpi.h"

// Information for view
#define C_VISUAL_ANGLE              120.

// Information for time evolution
#define C_TIME_STEP                 20000
#define C_TIME_INTERVAL             0.05

// Information for images
#define C_IMPORT_IMAGE_FILENAME     "input.png"
#define C_EXPORT_IMAGE_FILENAME     "output.png"
#define C_EXPORT_IMAGE_WIDTH        100
#define C_EXPORT_IMAGE_HEIGHT       100          // it should be a multiple of num_procs.
#define C_EXPORT_IMAGE_TOTAL_PIXELS (C_EXPORT_IMAGE_WIDTH * C_EXPORT_IMAGE_HEIGHT)

// Information for job division
#define C_NUMBER_OF_DIVISION        4

void gather_results(RESULT *result, int number_of_results, RESULT *collection);

int main(int argc, char **argv)
{
    BLACKHOLE_LAB *blackhole_lab = NULL;
    RESULT *result = NULL;
    int number_of_results;
    RESULT *collection = NULL;
    int x_start, x_end;
    int y_start, y_end;
    int i, j;

    collection = (RESULT*)calloc(C_EXPORT_IMAGE_TOTAL_PIXELS, sizeof(RESULT));
    
    // MPI Initialize
    MPI_Init(&argc, &argv);
    int coreid, ncores;
    MPI_Comm_rank(MPI_COMM_WORLD, &coreid);
    MPI_Comm_size(MPI_COMM_WORLD, &ncores);
    
    // Work division
    int workWidth = C_EXPORT_IMAGE_WIDTH / 2;
    int workHeight = C_EXPORT_IMAGE_HEIGHT / 2;
    int workAreas[4][2] = {{0, 0}, {workWidth, 0}, {0, workHeight}, {workWidth, workHeight}};

    int* currentWork = workAreas[coreid];
    x_start = currentWork[0];
    x_end   = x_start + workWidth - 1;
    y_start = currentWork[1];
    y_end = y_start + workHeight - 1;

    blackhole_lab = BLACKHOLE_LAB_create(
        C_EXPORT_IMAGE_WIDTH, x_start, x_end,
        C_EXPORT_IMAGE_HEIGHT, y_start, y_end,
        C_VISUAL_ANGLE, C_TIME_STEP, C_TIME_INTERVAL
    );

    if (NULL == blackhole_lab) {
        BLACKHOLE_LAB_destroy(blackhole_lab);
        fprintf(stderr, "Creating BLACKHOLE_LAB is Failed.\n");
        return -1;
    }

    BLACKHOLE_LAB_start_experiment(blackhole_lab);

    number_of_results = workWidth * workHeight;

    result = (RESULT*)calloc(number_of_results, sizeof(RESULT));

    BLACKHOLE_LAB_dump_results(blackhole_lab, result);

    BLACKHOLE_LAB_destroy(blackhole_lab);

    gather_results(result, number_of_results, collection);
    
    MPI_Barrier(MPI_COMM_WORLD);

    free(result);

    if (coreid == 0) {
        BLACKHOLE_LAB_make_vision(collection, C_IMPORT_IMAGE_FILENAME, C_EXPORT_IMAGE_FILENAME, C_EXPORT_IMAGE_WIDTH, C_EXPORT_IMAGE_HEIGHT);
    }

    free(collection);
    
    MPI_Finalize();
    return 0;
}

 void gather_results(RESULT *result, int number_of_results, RESULT *collection)
 {
     int i;
     MPI_Datatype mpi_type;
     int member_block_length[4] = {1, 1, 1, 1};
     MPI_Datatype member_type[4] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
     MPI_Aint member_displacement[4] = {offsetof(RESULT, status), offsetof(RESULT, x), offsetof(RESULT, y), offsetof(RESULT, pixel_index)};
 
     MPI_Type_struct(4, member_block_length, member_displacement, member_type, &mpi_type);
     MPI_Type_commit(&mpi_type);
     
     MPI_Gather(result, number_of_results, mpi_type, collection, number_of_results, mpi_type, 0, MPI_COMM_WORLD);
     
     MPI_Type_free(&mpi_type);
 }