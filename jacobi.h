#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)<(y)?(x):(y))


int*** randomize(int*** m, int d, int h, int w, int seed);
int*** jacobi(int*** original, int d, int h, int w, int step);
int*** create(int d, int h, int w);
int*** copy(int*** dst, int*** src, int d, int h, int w);
int destroy(int*** m, int d, int h);

//#define ORDER_XYZ

#ifndef ORDER_XYZ
#ifndef ORDER_ZYX
#define ORDER_ZYX
#endif
#endif