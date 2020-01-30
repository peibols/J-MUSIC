// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "data_struct.h"

#ifndef PI
#define PI (3.14159265358979324)
#endif

#ifndef hbarc
#define hbarc (0.19733)
#endif

#ifndef default_tol
#define default_tol (1.0e-8)
#endif

#define absol(a) ((a) > 0 ? (a) : (-(a)))
#define maxi(a, b) ((a) > (b) ? (a) : (b))
#define mini(a, b) ((a) < (b) ? (a) : (b))
#define sgn(x) ((x) < 0.0 ? (-1.0) : (1.0))
#define mtheta(x) ((x) < 0.0 ? (0.0) : (1.0))
#define SQR(a) ((a)*(a))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define DMAX(a,b) ((a) > (b) ? (a) : (b))
#define DMIN(a,b) ((a) < (b) ? (a) : (b))
#define SIGN(a,b) ((b) >= 0.0 ? absol(a): -absol(a)) 
/* added gmn = Minkowski metric to be used in sums */
#define gmn(a) ((a) == 0 ? (-1.0) : (1.0))

#define BT_BUF_SIZE 500

//! This is a utility class which contains a collection of helper functions.

namespace Util {
    double **mtx_malloc(int , int );

    void mtx_free(double **, int, int);
    
    int IsFile(std::string);
    
    std::string StringFind4(std::string file_name, std::string str_in);
    
    double lin_int(double x1,double x2,double f1,double f2,double x);

    double four_dimension_linear_interpolation(
            double* lattice_spacing, double fraction[2][4], double**** cube);
    double three_dimension_linear_interpolation(
            double* lattice_spacing, double fraction[2][3], double*** cube);
    int binary_search(double* array, int length, double x);
    void print_backtrace_errors();
    
    int map_2d_idx_to_1d(int a, int b);
    void map_1d_idx_to_2d(int idx_1d, int &a, int &b);

    Mat4x4 UnpackVecToMatrix(const Arr10 &in_vector);
    Mat4x4 UnpackVecToMatrix(const ViscousVec &in_vector);
    void LorentzBoost4Vector(double Tin[4], double *u, double Tout[4]);
}

#endif
