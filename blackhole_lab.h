//
//  blackhole_lab.h
//  geodesics
//
//  Created by Chan Park on 2016. 7. 20..
//  Copyright © 2016년 Korea Center for Numerical Relativity. All rights reserved.
//

#ifndef blackhole_lab_h
#define blackhole_lab_h

#include "bundle.h"

typedef enum {HIT, FALL, OUTSIDE, YET} BUNDLE_STATUS;

typedef struct {
    BUNDLE *bundle;
    int *pixel_index;
} BLACKHOLE_LAB;

typedef struct {
    BUNDLE_STATUS status;
    double x;
    double y;
    int pixel_index;
} RESULT;

BLACKHOLE_LAB *BLACKHOLE_LAB_create(int x_number_of_points, int x_start, int x_end, int y_number_of_points, int y_start, int y_end, double visual_angle, double time_step, double time_interval);
void BLACKHOLE_LAB_start_experiment(BLACKHOLE_LAB *blackhole_lab);
void BLACKHOLE_LAB_destroy(BLACKHOLE_LAB *blackhole_lab);
void BLACKHOLE_LAB_dump_results(BLACKHOLE_LAB *blackhole_lab, RESULT *result);
void BLACKHOLE_LAB_make_vision(RESULT *result, char *import_filename, char *export_filename, int export_width, int export_height);

#endif /* blackhole_lab_h */
