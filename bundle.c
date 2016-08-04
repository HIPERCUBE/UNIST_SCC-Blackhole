//
//  bundle.c
//  geodesics
//
//  Created by Chan Park on 2016. 7. 20..
//  Copyright © 2016년 Korea Center for Numerical Relativity. All rights reserved.
//

#include "bundle.h"
#include <stdlib.h>
#include <stdio.h>

static void VARIABLE_SYSTEM_set_pointer(VARIABLE_SYSTEM *variable_system, double *variable_pointer, double *parameter_pointer);
static void COORDINATE_SYSTEM_set_pointer(COORDINATE_SYSTEM *coordinate_system, double *coordinate_pointer);
static void BUNDLE_set_coordinate(BUNDLE *bundle, UNIFORM_GRID *grid, COORDINATE_INDEX coordinate_index, BUNDLE_INDEX bundle_index, double *coordinate_value);
static double UNIFORM_GRID_get_coordinate(UNIFORM_GRID *grid, int grid_index);




void VARIABLE_SYSTEM_set_pointer(VARIABLE_SYSTEM *variable_system, double *variable_pointer, double *parameter_pointer)
{
    variable_system->variable = variable_pointer;
    variable_system->time = parameter_pointer;
}

void COORDINATE_SYSTEM_set_pointer(COORDINATE_SYSTEM *coordinate_system, double *coordinate_pointer)
{
    coordinate_system->coordinate = coordinate_pointer;
}

int UNIFORM_GRID_get_number_of_points(UNIFORM_GRID * grid)
{
    return (grid->include_end_point) ? grid->number_of_cells + 1 : grid->number_of_cells;
}

double UNIFORM_GRID_get_coordinate(UNIFORM_GRID *grid, int grid_index)
{
    return grid->start + ((grid->number_of_cells == 0) ? 0. : ((grid->end - grid->start) / grid->number_of_cells) * grid_index);
}

BUNDLE *BUNDLE_create(
    int dimension, UNIFORM_GRID *grid,
    int number_of_variables, COORDINATE_FUNCTION *evaluate_initial_value, VARIABLE_FUNCTION *evaluate_derivative,
    int time_step, double time_interval, VARIABLE_CRITERION stop_condition
)
{
    BUNDLE_INDEX i;
    VARIABLE_INDEX j;
    COORDINATE_INDEX k;
    BUNDLE *bundle = NULL;
    double coordinate_value[dimension];

    ALLOCATE(BUNDLE, bundle, bundle, 1, BUNDLE)

    bundle->number_of_total_points = 1;
    bundle->dimension = dimension;
    bundle->coordinate = NULL;
    bundle->coordinate_system = NULL;
    bundle->number_of_variables = number_of_variables;
    bundle->variable = NULL;
    bundle->time = 0.;
    bundle->variable_system = NULL;
    bundle->derivative = NULL;
    bundle->time_temp = 0.;
    bundle->variable_system_temp.variable = NULL;
    bundle->variable_system_temp.time = NULL;
    bundle->derivative_temp = NULL;
    bundle->evaluate_derivative = NULL;
    bundle->time_step = time_step;
    bundle->current_time_step = 0;
    bundle->time_interval = time_interval;
    bundle->factor_value[0] = 0.;
    bundle->factor_value[1] = 1. / 2.;
    bundle->factor_value[2] = 1. / 2.;
    bundle->factor_value[3] = 1.;
    bundle->factor_delta[0] = 1. / 6.;
    bundle->factor_delta[1] = 1. / 3.;
    bundle->factor_delta[2] = 1. / 3.;
    bundle->factor_delta[3] = 1. / 6.;
    bundle->stop_condition = stop_condition;
    
    for (k = 0; k < dimension; k++) {
        bundle->number_of_total_points *= UNIFORM_GRID_get_number_of_points(grid + k);

    }

    ALLOCATE(BUNDLE, bundle, bundle->coordinate, bundle->number_of_total_points * dimension, double)
    BUNDLE_set_coordinate(bundle, grid, 0, 0, coordinate_value);
    ALLOCATE(BUNDLE, bundle, bundle->coordinate_system, bundle->number_of_total_points, COORDINATE_SYSTEM)
    for (i = 0; i < bundle->number_of_total_points; i++) {
        COORDINATE_SYSTEM_set_pointer(bundle->coordinate_system + i, bundle->coordinate + i * dimension);
    }
    ALLOCATE(BUNDLE, bundle, bundle->variable, bundle->number_of_total_points * number_of_variables, double)
    for (i = 0; i < bundle->number_of_total_points; i++) {
        for (j = 0; j < number_of_variables; j++) {
            bundle->variable[i * number_of_variables + j] = evaluate_initial_value[j](bundle->coordinate_system + i);
        }
    }
    ALLOCATE(BUNDLE, bundle, bundle->variable_system, bundle->number_of_total_points, VARIABLE_SYSTEM)
    for (i = 0; i < bundle->number_of_total_points; i++) {
        VARIABLE_SYSTEM_set_pointer(bundle->variable_system + i, bundle->variable + i * number_of_variables, &(bundle->time));
    }
    ALLOCATE(BUNDLE, bundle, bundle->derivative, number_of_variables, double)
    ALLOCATE(BUNDLE, bundle, bundle->variable_temp, number_of_variables, double)
    VARIABLE_SYSTEM_set_pointer(&bundle->variable_system_temp, bundle->variable_temp, &(bundle->time_temp));
    ALLOCATE(BUNDLE, bundle, bundle->derivative_temp, number_of_variables, double)
    COPY_ARRAY(BUNDLE, bundle, evaluate_derivative, number_of_variables, VARIABLE_FUNCTION)
    ALLOCATE(BUNDLE, bundle, bundle->is_stopped, bundle->number_of_total_points, bool)
    for (i = 0; i < bundle->number_of_total_points; i++) {
        bundle->is_stopped[i] = false;
    }
    
    return bundle;
}

void BUNDLE_destroy(BUNDLE *bundle)
{
    if (NULL != bundle) {
        DEALLOCATE(bundle->coordinate)
        DEALLOCATE(bundle->coordinate_system)
        DEALLOCATE(bundle->variable)
        DEALLOCATE(bundle->variable_system)
        DEALLOCATE(bundle->derivative)
        DEALLOCATE(bundle->variable_temp)
        DEALLOCATE(bundle->derivative_temp)
        DEALLOCATE(bundle->evaluate_derivative)
        DEALLOCATE(bundle->is_stopped)
        DEALLOCATE(bundle)
    }
}

void BUNDLE_set_coordinate(BUNDLE *bundle, UNIFORM_GRID *grid, COORDINATE_INDEX coordinate_index, BUNDLE_INDEX bundle_index, double *coordinate_value)
{
    int i;
    int number_of_points_on_grid;

    if (coordinate_index == bundle->dimension) {
        for (i = 0; i < bundle->dimension; i++) {
            *(bundle->coordinate + bundle_index * bundle->dimension + i) = coordinate_value[i];
        }
    } else {
        number_of_points_on_grid = UNIFORM_GRID_get_number_of_points(grid + coordinate_index);
        for (i = 0; i < number_of_points_on_grid; i++) {
            coordinate_value[coordinate_index] = UNIFORM_GRID_get_coordinate(grid + coordinate_index, i);
            BUNDLE_set_coordinate(bundle, grid, coordinate_index + 1, bundle_index * number_of_points_on_grid + i, coordinate_value);
        }
    }
}

void BUNDLE_evolve(BUNDLE *bundle)
{
    while (BUNDLE_one_step(bundle))
        ;
}

bool BUNDLE_one_step(BUNDLE *bundle)
{
    bool ret = false;
    BUNDLE_INDEX i;
    VARIABLE_INDEX j;
    int k;

    if (bundle->current_time_step < bundle->time_step) {
        for (i = 0; i < bundle->number_of_total_points; i++) {
            if (NULL != bundle->stop_condition && !bundle->is_stopped[i]) {
                bundle->is_stopped[i] = bundle->stop_condition(&(bundle->variable_system[i]));
            }
            if (NULL == bundle->stop_condition || !bundle->is_stopped[i]) {
                bundle->time = bundle->current_time_step * bundle->time_interval;
                
                for (j = 0; j < bundle->number_of_variables; j++) {
                    bundle->derivative[j] = 0.;
                }
                
                for (k = 0; k < 4; k++) {
                    *(bundle->variable_system_temp.time) = bundle->time + bundle->factor_value[k] * bundle->time_interval;
                    for (j = 0; j < bundle->number_of_variables; j++) {
                        bundle->variable_system_temp.variable[j] = (bundle->variable_system + i)->variable[j] + bundle->factor_value[k] * bundle->time_interval * bundle->derivative_temp[j];
                    }
                    for (j = 0; j < bundle->number_of_variables; j++) {
                        bundle->derivative_temp[j] = bundle->evaluate_derivative[j](&(bundle->variable_system_temp));
                    }
                    for (j = 0; j < bundle->number_of_variables; j++) {
                        bundle->derivative[j] += bundle->factor_delta[k]  * bundle->derivative_temp[j];
                    }
                }
                
                for (j = 0; j < bundle->number_of_variables; j++) {
                    (bundle->variable_system + i)->variable[j] += bundle->derivative[j] * bundle->time_interval;
                }
                
                ret = true;
            }
        }
    }
    
    if (true == ret) bundle->current_time_step++;
    
    return ret;
}

