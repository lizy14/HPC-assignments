#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifndef max
#define max(x, y) ((x)>(y)?(x):(y))
#endif
#ifndef min
#define min(x, y) ((x)<(y)?(x):(y))
#endif

#define type double

type*** randomize(type*** m, size_t d, size_t h, size_t w, int seed);
type*** jacobi(type*** original, size_t d, size_t h, size_t w, int step);
type*** create(size_t d, size_t h, size_t w);
type*** copy(type*** dst, type*** src, size_t d, size_t h, size_t w);
void destroy(type*** m, size_t d, size_t h);

//#define ORDER_XYZ

#ifndef ORDER_XYZ
#ifndef ORDER_ZYX
#define ORDER_ZYX
#endif
#endif