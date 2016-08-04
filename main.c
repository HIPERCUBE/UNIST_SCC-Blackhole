#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "blackhole_lab.h"

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

//void gather_results(RESULT *result, int number_of_results, RESULT *collection);

int main(int argc, char **argv)
{
    BLACKHOLE_LAB *blackhole_lab = NULL;
    RESULT *result = NULL;
    int number_of_results;
    RESULT *collection = NULL;
    int x_start, x_end;
    int y_start, y_end;
    int my_rank, num_procs;
    int i, j;

    collection = (RESULT*)calloc(C_EXPORT_IMAGE_TOTAL_PIXELS, sizeof(RESULT));

    num_procs = C_NUMBER_OF_DIVISION;

    x_start = 0;
    x_end   = C_EXPORT_IMAGE_WIDTH - 1;

    for (my_rank = 0, i = 0; my_rank < num_procs; my_rank++) {
        y_start = (C_EXPORT_IMAGE_HEIGHT / num_procs) * my_rank;
        y_end = (C_EXPORT_IMAGE_HEIGHT / num_procs) * (my_rank + 1) - 1;

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

        number_of_results = C_EXPORT_IMAGE_WIDTH * C_EXPORT_IMAGE_HEIGHT / num_procs;

        result = (RESULT*)calloc(number_of_results, sizeof(RESULT));

        BLACKHOLE_LAB_dump_results(blackhole_lab, result);

        BLACKHOLE_LAB_destroy(blackhole_lab);

        for (j = 0; j < number_of_results; j++, i++) {
            collection[i] = result[j];
        }

        free(result);
    }

    BLACKHOLE_LAB_make_vision(collection, C_IMPORT_IMAGE_FILENAME, C_EXPORT_IMAGE_FILENAME, C_EXPORT_IMAGE_WIDTH, C_EXPORT_IMAGE_HEIGHT);

    free(collection);

    return 0;
}

/*
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
 */
