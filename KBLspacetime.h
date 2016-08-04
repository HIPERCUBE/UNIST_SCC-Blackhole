//
//  KBLspacetime.h
//  geodesics
//
//  Created by Chan Park on 2016. 7. 27..
//  Copyright © 2016년 Korea Center for Numerical Relativity. All rights reserved.
//

#ifndef KBLspacetime_h
#define KBLspacetime_h

#define C_M 1.
#define C_a 0.

enum KBL_COORDINATE {KBL_r, KBL_theta, KBL_phi, C_NUMBER_OF_KBL_COORDINATE};

double KBL_N(double r, double theta, double phi);
double KBL_partial_N(int i, double r, double theta, double phi);
double KBL_beta(int i, double r, double theta, double phi);
double KBL_partial_beta(int i, int j, double r, double theta, double phi);
double KBL_K(int i, int j, double r, double theta, double phi);
double KBL_gamma(int i, int j, double r, double theta, double phi);
double KBL_gamma_inverse(int i, int j, double r, double theta, double phi);
double KBL_Gamma(int i, int j, int k, double r, double theta, double phi);

#endif /* KBLspacetime_h */
