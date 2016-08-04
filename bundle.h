//
//  bundle.h
//  geodesics
//
//  Created by Chan Park on 2016. 7. 20..
//  Copyright © 2016년 Korea Center for Numerical Relativity. All rights reserved.
//

#ifndef bundle_h
#define bundle_h

#include <stdbool.h>
//#include "worldline.h"

typedef int VARIABLE_INDEX;

typedef struct _VARIABLE_SYSTEM {
    double *variable;
    double *time;
} VARIABLE_SYSTEM;

typedef double (*VARIABLE_FUNCTION)(VARIABLE_SYSTEM*);
typedef bool (*VARIABLE_CRITERION)(VARIABLE_SYSTEM*);

typedef int COORDINATE_INDEX;

typedef struct _COORDINATE_SYSTEM {
    double *coordinate;
} COORDINATE_SYSTEM;

typedef double (*COORDINATE_FUNCTION)(COORDINATE_SYSTEM*);

typedef struct _UNIFORM_GRID {
    int number_of_cells;
    double start;
    double end;
    bool include_end_point;
} UNIFORM_GRID;

int UNIFORM_GRID_get_number_of_points(UNIFORM_GRID *grid);

typedef int BUNDLE_INDEX;

typedef struct _BUNDLE {
    int number_of_total_points;
    int dimension;
    double *coordinate;
    COORDINATE_SYSTEM *coordinate_system;
    int number_of_variables;
    double *variable;
    double time;
    VARIABLE_SYSTEM *variable_system;
    double *derivative;
    double *variable_temp;
    double time_temp;
    VARIABLE_SYSTEM variable_system_temp;
    double *derivative_temp;
    VARIABLE_FUNCTION *evaluate_derivative;
    int time_step;
    int current_time_step;
    double time_interval;
    double factor_value[4];
    double factor_delta[4];
    VARIABLE_CRITERION stop_condition;
    bool *is_stopped;
} BUNDLE;

BUNDLE *BUNDLE_create(int dimension, UNIFORM_GRID *coordinate, int number_of_variables, COORDINATE_FUNCTION *evaluate_initial_value, VARIABLE_FUNCTION *evaluate_derivative, int time_step, double time_interval, VARIABLE_CRITERION stop_condition);
void BUNDLE_destroy(BUNDLE *bundle);
void BUNDLE_evolve(BUNDLE *bundle);
bool BUNDLE_one_step(BUNDLE *bundle);

#define ALLOCATE(struct_name, instance, pointer, size, type) \
    if (0 != (size) && NULL == ((pointer) = calloc((size), sizeof(type)))) { \
        struct_name##_destroy(instance); \
        return NULL; \
    }

#define COPY_ARRAY(struct_name, instance, member_name, size, type) \
    ALLOCATE(struct_name, instance, (instance)->member_name, size, type) \
    for (i = 0; i < (size); i++) { \
        (instance)->member_name[i] = member_name[i]; \
    }

#define DEALLOCATE(pointer) \
    if (NULL != (pointer)) { \
        free(pointer); \
        (pointer) = NULL; \
    }

#endif /* bundle_h */
