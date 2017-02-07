/** ============================================================================

    Description: Dynamic Slicing for Hypothesis Testing
         Author: Chao Ye, Bioinformatics Division, TNLIST, Tsinghua University
                 Bo Jiang, Department of Statistics, Harvard University
          Email: yec09@mails.tsinghua.edu.cn, bjiang@fas.harvard.edu
    Last Update: Feb 14, 2014

============================================================================ **/

#ifndef DYNSLICING_H_INCLUDED
#define DYNSLICING_H_INCLUDED

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

//#include <unistd.h>

/** ========================= constant definition ========================== **/

const double EPSILON = 1e-6;  //  Double variable less than it will be treated as 0
const int GETIN = 65536;      //  Maximum of buffer
const int NAMELEN = 256;      //  Length of file path, rownames or colnames
const int SIMNUM = 1000;      //  Default simulation times
const int NFILE = 10;         //  Number of files

/** ========================= function declaration ========================= **/

double dynslicing_1eqp(double *y, int len, double lambda);
double dynslicing_one(double *y, int len, double lambda, double alpha);
double dynslicing_keqp(int *x, int len, int dim, double lambda);
double dynslicing_keqp(int *x, int len, int dim, double lambda, std::vector<int> &slices);
double dynslicing_k(int *x, int len, int dim, double lambda);
double dynslicing_k(int *x, int len, int dim, double lambda, std::vector<int> &slices);

#endif // DYNSLICING_H_INCLUDED
