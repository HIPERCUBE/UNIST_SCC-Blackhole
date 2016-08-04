//
//  KBLspacetime.c
//  geodesics
//
//  Created by Chan Park on 2016. 7. 27..
//  Copyright © 2016년 Korea Center for Numerical Relativity. All rights reserved.
//

#include "KBLspacetime.h"
#include <math.h>


/*
static double Sigma(double M, double a, double r, double theta, double phi)
{
    return pow(r, 2.) + pow(a * cos(theta), 2.);
}

static double Delta(double M, double a, double r, double theta, double phi)
{
    return pow(r, 2.) - 2. * M * r + pow(a, 2.);
}

static double Omega(double M, double a, double r, double theta, double phi)
{
    return Sigma(M, a, r, theta, phi) * Delta(M, a, r, theta, phi) + 2. * M * r * (pow(r, 2.) + pow(a, 2.));
}
*/

// Schwarzschild case

// lapse function
double KBL_N(double r, double theta, double phi)
{
    return sqrt(1. - 2. * C_M / r);
}

double KBL_partial_N(int i, double r, double theta, double phi)
{
    if (i == KBL_r)
        return C_M / (r * r) / KBL_N(r, theta, phi);

    return 0.;
}

// shift vector
double KBL_beta(int i, double r, double theta, double phi)
{
    return 0.;
}

double KBL_partial_beta(int i, int j, double r, double theta, double phi)
{
    return 0.;
}

// extrinsic curvature
double KBL_K(int i, int j, double r, double theta, double phi)
{
    return 0.;
}

#define IF(i, j, c1, c2, value) \
    if (((i) == (c1) && (j) == (c2)) || ((j) == (c1) && (i) == (c2))) return value;

// Induced metric
double KBL_gamma(int i, int j, double r, double theta, double phi)
{
    IF(i, j, KBL_r    , KBL_r    , 1. / (1. - 2. * C_M / r) )
    IF(i, j, KBL_theta, KBL_theta, pow(r, 2.)             )
    IF(i, j, KBL_phi  , KBL_phi  , pow(r * sin(theta), 2.))

    return 0.;
}

double KBL_gamma_inverse(int i, int j, double r, double theta, double phi)
{
    IF(i, j, KBL_r    , KBL_r    , 1. - 2. * C_M / r         )
    IF(i, j, KBL_theta, KBL_theta, pow(r, -2.)             )
    IF(i, j, KBL_phi  , KBL_phi  , pow(r * sin(theta), -2.))

    return 0.;
}

#define cot(radian) \
    (((radian) == M_PI / 2.) ? 0. : 1. / tan(radian))

// Christoffel symbol of induced metric
double KBL_Gamma(int i, int j, int k, double r, double theta, double phi)
{
    switch (i) {
    case KBL_r:
        IF(j, k, KBL_r    , KBL_r    , -C_M / r / (r - 2. * C_M)              )
        IF(j, k, KBL_theta, KBL_theta, -(r - 2. * C_M)                      )
        IF(j, k, KBL_phi  , KBL_phi  , -(r - 2. * C_M) * pow(sin(theta), 2.))
        break;
    case KBL_theta:
        IF(j, k, KBL_r    , KBL_theta, 1. / r                             )
        IF(j, k, KBL_phi  , KBL_phi  , -sin(theta) * cos(theta)           )
        break;
    case KBL_phi:
        IF(j, k, KBL_r    , KBL_phi  , 1. / r                             )
        IF(j, k, KBL_theta, KBL_phi  , cot(theta)                         )
        break;
    }

    return 0.;
}

/* Flat case

// lapse function
double KBL_N(double M, double a, double r, double theta, double phi)
{
    return 1.;
}

double KBL_partial_N(double M, double a, int i, double r, double theta, double phi)
{
    return 0.;
}

// shift vector
double KBL_beta(double M, double a, int i, double r, double theta, double phi)
{
    return 0.;
}

double KBL_partial_beta(double M, double a, int i, int j, double r, double theta, double phi)
{
    return 0.;
}

// extrinsic curvature
double KBL_K(double M, double a, int i, int j, double r, double theta, double phi)
{
    return 0.;
}

#define IF(i, j, c1, c2, value) \
    if (((i) == (c1) && (j) == (c2)) || ((j) == (c1) && (i) == (c2))) return value;

// Induced metric
double KBL_gamma(double M, double a, int i, int j, double r, double theta, double phi)
{
    IF(i, j, KBL_r    , KBL_r    , 1.                     )
    IF(i, j, KBL_theta, KBL_theta, pow(r, 2.)             )
    IF(i, j, KBL_phi  , KBL_phi  , pow(r * sin(theta), 2.))

    return 0.;
}

double KBL_gamma_inverse(double M, double a, int i, int j, double r, double theta, double phi)
{
    IF(i, j, KBL_r    , KBL_r    , 1.                     )
    IF(i, j, KBL_theta, KBL_theta, pow(r, -2.)             )
    IF(i, j, KBL_phi  , KBL_phi  , pow(r * sin(theta), -2.))

    return 0.;
}

#define cot(radian) \
    (((radian) == M_PI / 2.) ? 0. : 1. / tan(radian))

// Christoffel symbol of induced metric
double KBL_Gamma(double M, double a, int i, int j, int k, double r, double theta, double phi)
{
    switch (i) {
    case KBL_r:
        IF(j, k, KBL_theta, KBL_theta, -r                      )
        IF(j, k, KBL_phi  , KBL_phi  , -r * pow(sin(theta), 2.))
        break;
    case KBL_theta:
        IF(j, k, KBL_r    , KBL_theta, 1. / r                  )
        IF(j, k, KBL_phi  , KBL_phi  , -sin(theta) * cos(theta))
        break;
    case KBL_phi:
        IF(j, k, KBL_r    , KBL_phi  , 1. / r                  )
        IF(j, k, KBL_theta, KBL_phi  , cot(theta)              )
        break;
    }

    return 0.;
}

*/
