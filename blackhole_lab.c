//
//  blackhole_lab.c
//  geodesics
//
//  Created by Chan Park on 2016. 7. 20..
//  Copyright © 2016년 Korea Center for Numerical Relativity. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <png.h>
#include "blackhole_lab.h"
#include "KBLspacetime.h"

#define max(a, b) (((a) > (b)) ? (a) : (b) )

// Parameters for lab deployment
#define C_L1    10.                             // distance between observer and black hole
#define C_L2    10.                             // distance between black hole and picture
#define C_l     10.                             // length of side of picture
#define C_L     max(C_L1, C_L2)                 // criteria of outside boundary

// Settings for light bundle
enum _KIND_OF_COORDINATES {x, y, C_DIMENSION};

static double N(VARIABLE_SYSTEM *variable_system);
static double partial_N(int KBL_i, VARIABLE_SYSTEM *variable_system);
static double beta(int KBL_i, VARIABLE_SYSTEM *variable_system);
static double partial_beta(int KBL_i, int KBL_j, VARIABLE_SYSTEM *variable_system);
static double gamma_inverse(int KBL_i, int KBL_j, VARIABLE_SYSTEM *variable_system);
static double K(int KBL_i, int KBL_j, VARIABLE_SYSTEM *variable_system);
static double Gamma(int KBL_i, int KBL_j, int KBL_k, VARIABLE_SYSTEM *variable_system);

enum _KIND_OF_VARIABLES {X_r, X_theta, X_phi, V_r, V_theta, V_phi, C_NUMBER_OF_VARIABLES};
static double X_r_zero(COORDINATE_SYSTEM *coordinate_system);
static double X_theta_zero(COORDINATE_SYSTEM *coordinate_system);
static double X_phi_zero(COORDINATE_SYSTEM *coordinate_system);
static double V_r_zero(COORDINATE_SYSTEM *coordinate_system);
static double V_theta_zero(COORDINATE_SYSTEM *coordinate_system);
static double V_phi_zero(COORDINATE_SYSTEM *coordinate_system);
static double X_dot(VARIABLE_SYSTEM *variable_system, int KBL_i);
static double X_r_dot(VARIABLE_SYSTEM *variable_system);
static double X_theta_dot(VARIABLE_SYSTEM *variable_system);
static double X_phi_dot(VARIABLE_SYSTEM *variable_system);
static double V_dot(VARIABLE_SYSTEM *variable_system, int KBL_coordinate_index);
static double V_r_dot(VARIABLE_SYSTEM *variable_system);
static double V_theta_dot(VARIABLE_SYSTEM *variable_system);
static double V_phi_dot(VARIABLE_SYSTEM *variable_system);
COORDINATE_FUNCTION evaluate_initial_value[C_NUMBER_OF_VARIABLES] = {
    X_r_zero, X_theta_zero, X_phi_zero,
    V_r_zero, V_theta_zero, V_phi_zero,
};
VARIABLE_FUNCTION evaluate_derivative[C_NUMBER_OF_VARIABLES] = {
    X_r_dot, X_theta_dot, X_phi_dot,
    V_r_dot, V_theta_dot, V_phi_dot,
};
static bool stop_condition(VARIABLE_SYSTEM *variable_system);
static bool outside_boundary(VARIABLE_SYSTEM *variable_system);
static bool cross_picture(VARIABLE_SYSTEM *variable_system);
static bool inside_photon_sphere(VARIABLE_SYSTEM *variable_system);
static double get_X_x(VARIABLE_SYSTEM *variable_system);
static double get_X_y(VARIABLE_SYSTEM *variable_system);
static double get_X_z(VARIABLE_SYSTEM *variable_system);

static double d;    // distance between observer and view window

BLACKHOLE_LAB *BLACKHOLE_LAB_create(
    int x_number_of_points, int x_start, int x_end,
    int y_number_of_points, int y_start, int y_end,
    double visual_angle,
    double time_step, double time_interval
)
{
    BLACKHOLE_LAB *blackhole_lab = NULL;
    UNIFORM_GRID grid[C_DIMENSION] = {
        {x_end - x_start, -0.5 + x_start  / ((double)(x_number_of_points - 1)), -0.5 + x_end  / ((double)(x_number_of_points - 1)), true},
        {y_end - y_start, +0.5 - y_start  / ((double)(y_number_of_points - 1)), +0.5 - y_end  / ((double)(y_number_of_points - 1)), true}
    };
    BUNDLE_INDEX i;
    
    d = 1. / (sqrt(2.) * tan(visual_angle / 2. / 180. * M_PI));

    ALLOCATE(BLACKHOLE_LAB, blackhole_lab, blackhole_lab, 1, BLACKHOLE_LAB)

    blackhole_lab->bundle = BUNDLE_create(
        C_DIMENSION, grid,
        C_NUMBER_OF_VARIABLES, evaluate_initial_value, evaluate_derivative,
        time_step, time_interval, stop_condition
    );
    if (NULL == blackhole_lab->bundle) {
        BLACKHOLE_LAB_destroy(blackhole_lab);
        return NULL;
    }

    ALLOCATE(BLACKHOLE_LAB, blackhole_lab, blackhole_lab->pixel_index, blackhole_lab->bundle->number_of_total_points, int)

    for (i = 0; i < blackhole_lab->bundle->number_of_total_points; i++) {
        blackhole_lab->pixel_index[i] =
            + (x_start + i / UNIFORM_GRID_get_number_of_points(grid + 1)) * y_number_of_points
            + y_start + i % UNIFORM_GRID_get_number_of_points(grid + 1);
    }

    return blackhole_lab;
}

void BLACKHOLE_LAB_destroy(BLACKHOLE_LAB *blackhole_lab)
{
    if (NULL != blackhole_lab) {
        BUNDLE_destroy(blackhole_lab->bundle);
        DEALLOCATE(blackhole_lab->pixel_index)
        DEALLOCATE(blackhole_lab)
    }
}

void BLACKHOLE_LAB_start_experiment(BLACKHOLE_LAB *blackhole_lab)
{

//    FILE *fp;
//    BUNDLE_INDEX i;
//
//    fp = fopen("bundle.txt", "w");
//    do {
//        for (i = 0; i < blackhole_lab->bundle->number_of_total_points; i++) {
//            fprintf(fp, "%e %e %e "
//                ,get_X_x(blackhole_lab->bundle->worldline[i])
//                ,get_X_y(blackhole_lab->bundle->worldline[i])
//                ,get_X_z(blackhole_lab->bundle->worldline[i])
//            );
//        }
//        fprintf(fp, "\n");
//    } while (one_step_bundle(blackhole_lab->bundle));
//    fclose(fp);

    BUNDLE_evolve(blackhole_lab->bundle);
}

void BLACKHOLE_LAB_dump_results(BLACKHOLE_LAB *blackhole_lab, RESULT *result)
{
    BUNDLE_INDEX i;

    for (i = 0; i < blackhole_lab->bundle->number_of_total_points; i++) {
        if (cross_picture(blackhole_lab->bundle->variable_system + i)) {
            result[i].status = HIT;
        } else if (inside_photon_sphere(blackhole_lab->bundle->variable_system + i)) {
            result[i].status = FALL;
        } else if (outside_boundary(blackhole_lab->bundle->variable_system + i)) {
            result[i].status = OUTSIDE;
        } else {
            result[i].status = YET;
        }
        result[i].x = get_X_y(blackhole_lab->bundle->variable_system + i);
        result[i].y = get_X_z(blackhole_lab->bundle->variable_system + i);
        result[i].pixel_index = blackhole_lab->pixel_index[i];
    }
}

#define round(x) ((int)((x) + 0.5))

void BLACKHOLE_LAB_make_vision(RESULT *result, char *import_filename, char *export_filename, int export_width, int export_height)
{
    int i, j;

    FILE *import_fp;
    png_structp import_png;
    png_infop import_info;
    png_byte import_color_type;
    png_byte import_bit_depth;
    png_bytep *import_row_pointers;
    int import_width, import_height, import_size;
    int import_x, import_y;

    FILE *export_fp;
    png_structp export_png;
    png_infop export_info;
    png_bytep *export_row_pointers;
    int export_x, export_y;

    // Import process
    import_fp = fopen(import_filename, "rb");

    import_png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!import_png) abort();

    import_info = png_create_info_struct(import_png);
    if (!import_info) abort();

    if(setjmp(png_jmpbuf(import_png))) abort();

    png_init_io(import_png, import_fp);

    png_read_info(import_png, import_info);

    import_width        = png_get_image_width   (import_png, import_info);
    import_height       = png_get_image_height  (import_png, import_info);
    import_color_type   = png_get_color_type    (import_png, import_info);
    import_bit_depth    = png_get_bit_depth     (import_png, import_info);

    // Read any color_type into 8bit depth, RGBA format.
    // See http://www.libpng.org/pub/png/libpng-manual.txt

    if (import_bit_depth == 16) {
        png_set_strip_16(import_png);
    }

    if (import_color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(import_png);
    }

    // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
    if (import_color_type == PNG_COLOR_TYPE_GRAY && import_bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(import_png);
    }

    if(png_get_valid(import_png, import_info, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(import_png);
    }

    // These color_type don't have an alpha channel then fill it with 0xff.
    if (import_color_type == PNG_COLOR_TYPE_RGB || import_color_type == PNG_COLOR_TYPE_GRAY || import_color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_filler(import_png, 0xFF, PNG_FILLER_AFTER);
    }

    if (import_color_type == PNG_COLOR_TYPE_GRAY || import_color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(import_png);
    }

    png_read_update_info(import_png, import_info);

    import_row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * import_height);
    for (import_y = 0; import_y < import_height; import_y++) {
        import_row_pointers[import_y] = (png_byte*)malloc(png_get_rowbytes(import_png, import_info));
    }

    png_read_image(import_png, import_row_pointers);

    fclose(import_fp);


    // Export preparation
    export_fp = fopen(export_filename, "wb");
    if(!export_fp) abort();

    export_png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!export_png) abort();

    export_info = png_create_info_struct(export_png);
    if (!export_info) abort();

    if (setjmp(png_jmpbuf(export_png))) abort();

    png_init_io(export_png, export_fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(
        export_png,
        export_info,
        export_width, export_height,
        8,
        PNG_COLOR_TYPE_RGB,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(export_png, export_info);

    // allocation of export array
    export_row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * export_height);
    for (export_y = 0; export_y < export_height; export_y++) {
        export_row_pointers[export_y] = (png_byte*)malloc(png_get_rowbytes(export_png, export_info));
    }

    // initialize
    for (export_y = 0; export_y < export_height; export_y++) {
        for (export_x = 0; export_x < export_width; export_x++) {
            for (i = 0; i < 3; i++) {
                export_row_pointers[export_y][3 * export_x + i] = 0;
            }
        }
    }

    import_size = (import_height < import_width) ? import_height : import_width;
    for (i = 0; i < export_width * export_height; i++) {
        export_x = result[i].pixel_index / export_height;
        export_y = result[i].pixel_index % export_height;
        if (HIT == result[i].status && fabs(result[i].x) < C_l / 2. && fabs(result[i].y) < C_l / 2.) {
            import_x = round(+result[i].x / C_l * (import_size - 1) + (import_width  - 1) / 2.);
            import_y = round(-result[i].y / C_l * (import_size - 1) + (import_height - 1) / 2.);
            for (j = 0; j < 3; j++) {
                export_row_pointers[export_y][3 * export_x + j] = import_row_pointers[import_y][4 * import_x + j];
            }
        } else if (OUTSIDE == result[i].status) {
            for (j = 0; j < 3; j++) {
                export_row_pointers[export_y][3 * export_x + j] = 0xff;
            }
        }
    }

    // copy
//    for (export_y = 0; export_y < export_height; export_y++) {
//        for (export_x = 0; export_x < export_width; export_x++) {
//            for (i = 0; i < 3; i++) {
//                export_row_pointers[export_y][3 * export_x + i] = import_row_pointers[export_y][4 * export_x + i];
//            }
//        }
//    }

    // Write process
    png_write_image(export_png, export_row_pointers);
    png_write_end(export_png, NULL);

    fclose(export_fp);

    for(import_y = 0; import_y < import_height; import_y++) {
        free(import_row_pointers[import_y]);
    }
    free(import_row_pointers);

    for(export_y = 0; export_y < export_height; export_y++) {
        free(export_row_pointers[export_y]);
    }
    free(export_row_pointers);
}

void export_file(RESULT *result, char *export_filename, int export_width, int export_height)
{
    int i;
    FILE *fp;

    fp = fopen(export_filename, "w");
    for (i = 0; i < export_width * export_height; i++) {
        fprintf(fp, "%d %d %d\n", result[i].pixel_index / export_height, result[i].pixel_index % export_height, result[i].status);
    }
    fclose(fp);
}

// Functions for initial values and time derivatives

double X_r_zero(COORDINATE_SYSTEM *coordinate_system)
{
    return C_L1;
}

double X_theta_zero(COORDINATE_SYSTEM *coordinate_system)
{
    return M_PI / 2.;
}

double X_phi_zero(COORDINATE_SYSTEM *coordinate_system)
{
    return 0.;
}

double get_V_over_r(COORDINATE_SYSTEM *coordinate_system)
{
    return sqrt((1. - 2. * C_M / C_L1) / (
        + pow(coordinate_system->coordinate[x], 2.)
        + pow(coordinate_system->coordinate[y], 2.)
        + pow(d, 2.) / (1. - 2. * C_M / C_L1)
    ));
}

#define V sqrt(1. - 2. * C_M / C_L1)

double V_r_zero(COORDINATE_SYSTEM *coordinate_system)
{
    return - get_V_over_r(coordinate_system) * d;
}

double V_theta_zero(COORDINATE_SYSTEM *coordinate_system)
{
    return - get_V_over_r(coordinate_system) * coordinate_system->coordinate[y] / C_L1;
}

double V_phi_zero(COORDINATE_SYSTEM *coordinate_system)
{
    return + get_V_over_r(coordinate_system) * coordinate_system->coordinate[x] / C_L1;
}

double N(VARIABLE_SYSTEM *variable_system)
{
    return KBL_N(variable_system->variable[X_r], variable_system->variable[X_theta], variable_system->variable[X_phi]);
}

double partial_N(int KBL_i, VARIABLE_SYSTEM *variable_system)
{
    return KBL_partial_N(KBL_i, variable_system->variable[X_r], variable_system->variable[X_theta], variable_system->variable[X_phi]);
}

double beta(int KBL_i, VARIABLE_SYSTEM *variable_system)
{
    return KBL_beta(KBL_i, variable_system->variable[X_r], variable_system->variable[X_theta], variable_system->variable[X_phi]);
}

double partial_beta(int KBL_i, int KBL_j, VARIABLE_SYSTEM *variable_system)
{
    return KBL_partial_beta(KBL_i, KBL_j, variable_system->variable[X_r], variable_system->variable[X_theta], variable_system->variable[X_phi]);
}

double gamma_inverse(int KBL_i, int KBL_j, VARIABLE_SYSTEM *wvariable_system)
{
    return KBL_gamma_inverse(KBL_i, KBL_j, wvariable_system->variable[X_r], wvariable_system->variable[X_theta], wvariable_system->variable[X_phi]);
}

double Gamma(int KBL_i, int KBL_j, int KBL_k, VARIABLE_SYSTEM *variable_system)
{
    return KBL_Gamma(KBL_i, KBL_j, KBL_k, variable_system->variable[X_r], variable_system->variable[X_theta], variable_system->variable[X_phi]);
}

double K(int KBL_i, int KBL_j, VARIABLE_SYSTEM *variable_system)
{
    return KBL_K(KBL_i, KBL_j, variable_system->variable[X_r], variable_system->variable[X_theta], variable_system->variable[X_phi]);
}

double X_dot(VARIABLE_SYSTEM *variable_system, int KBL_i)
{
    return N(variable_system) * variable_system->variable[V_r + KBL_i] - beta(KBL_i, variable_system);
}

double X_r_dot(VARIABLE_SYSTEM *variable_system)
{
    return X_dot(variable_system, KBL_r);
}

double X_theta_dot(VARIABLE_SYSTEM *variable_system)
{
    return X_dot(variable_system, KBL_theta);
}

double X_phi_dot(VARIABLE_SYSTEM *variable_system)
{
    return X_dot(variable_system, KBL_phi);
}

double V_dot(VARIABLE_SYSTEM *variable_system, int KBL_i)
{
    int KBL_j, KBL_k;
    double ret = 0.;
    double temp;

    for (KBL_k = 0; KBL_k < C_NUMBER_OF_KBL_COORDINATE; KBL_k++) {
        temp = 0.;
        for (KBL_j = 0; KBL_j < C_NUMBER_OF_KBL_COORDINATE; KBL_j++) {
            temp += (
                +(2. * gamma_inverse(KBL_i, KBL_j, variable_system) - variable_system->variable[V_r + KBL_i] * variable_system->variable[V_r + KBL_j]) * K(KBL_j, KBL_k, variable_system)
                -Gamma(KBL_i, KBL_j, KBL_k, variable_system) * variable_system->variable[V_r + KBL_j]
            );
        }
        ret += temp * variable_system->variable[V_r + KBL_k];
    }

    ret *= N(variable_system);

    for (KBL_j = 0; KBL_j < C_NUMBER_OF_KBL_COORDINATE; KBL_j++) {
        ret += (variable_system->variable[V_r + KBL_i] * variable_system->variable[V_r + KBL_j] - gamma_inverse(KBL_i, KBL_j, variable_system)) * partial_N(KBL_j, variable_system);
        ret -= variable_system->variable[V_r + KBL_j] * partial_beta(KBL_j, KBL_i, variable_system);
    }

    return ret;
}

double V_r_dot(VARIABLE_SYSTEM *variable_system)
{
    return V_dot(variable_system, KBL_r);
}

double V_theta_dot(VARIABLE_SYSTEM *variable_system)
{
    return V_dot(variable_system, KBL_theta);
}

double V_phi_dot(VARIABLE_SYSTEM *variable_system)
{
    return V_dot(variable_system, KBL_phi);
}

double get_X_x(VARIABLE_SYSTEM *variable_system)
{
    return variable_system->variable[X_r] * sin(variable_system->variable[X_theta]) * cos(variable_system->variable[X_phi]);
}

double get_X_y(VARIABLE_SYSTEM *variable_system)
{
    return variable_system->variable[X_r] * sin(variable_system->variable[X_theta]) * sin(variable_system->variable[X_phi]);
}

double get_X_z(VARIABLE_SYSTEM *variable_system)
{
    return variable_system->variable[X_r] * cos(variable_system->variable[X_theta]);
}

bool outside_boundary(VARIABLE_SYSTEM *variable_system)
{
    return fabs(get_X_x(variable_system)) > C_L || fabs(get_X_y(variable_system)) > C_L || fabs(get_X_z(variable_system)) > C_L;
}

bool cross_picture(VARIABLE_SYSTEM *variable_system)
{
    return get_X_x(variable_system) < -C_L2 && fabs(get_X_y(variable_system)) < C_l / 2. && fabs(get_X_z(variable_system)) < C_l / 2.;
}

bool inside_photon_sphere(VARIABLE_SYSTEM *variable_system)
{
    return variable_system->variable[X_r] <= 3. * C_M;
}

bool stop_condition(VARIABLE_SYSTEM *variable_system)
{
    return outside_boundary(variable_system) || cross_picture(variable_system) || inside_photon_sphere(variable_system);
}
