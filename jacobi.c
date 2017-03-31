#include "jacobi.h"

type*** create(size_t d, size_t h, size_t w) {
	type*** m = (type***)malloc(d * sizeof(type**));
	for (size_t z = 0; z < d; z++) {
		m[z] = (type**)malloc(h * sizeof(type*));
		for (size_t y = 0; y < w; y++) {
			m[z][y] = (type*)malloc(w * sizeof(type));
		}
	}
	return m;
};

type*** copy(type*** dst, type*** src, size_t d, size_t h, size_t w) {
	for (size_t z = 0; z < d; z++) {
		for (size_t y = 0; y < h; y++) {
			for (size_t x = 0; x < w; x++) {
				dst[z][y][x] = src[z][y][x];
			}
		}
	}
	return dst;
}

void destroy(type*** m, size_t d, size_t h) {
	for (size_t z = 0; z < d; z++) {
		for (size_t y = 0; y < h; y++) {
			free(m[z][y]);
		}
		free(m[z]);
	}
	return;
}

type*** randomize(type*** m, size_t d, size_t h, size_t w, int seed) {
	srand(seed);
	for (size_t z = 0; z < d; z++) {
		for (size_t y = 0; y < h; y++) {
			for (size_t x = 0; x < w; x++) {
				m[z][y][x] = rand();
			}
		}
	}
	return m;
}

type*** jacobi(type*** original, size_t d, size_t h, size_t w, int step) {

	type*** extra = create(d, h, w);
	type*** grid = extra;
	type*** other_grid = original;

	for (int t = 0; t < step; t++) {		
		for (size_t z = 1; z < d - 1; z++) {
			for (size_t y = 1; y < h - 1; y++) {
				for (size_t x = 1; x < w - 1; x++) {
					// grid[x][y][z] = avg of neighbors in other_grid*/
					type sum_neibours =
						other_grid[z][y][x + 1] +
						other_grid[z][y][x - 1] +
						other_grid[z][y + 1][x] +
						other_grid[z][y - 1][x] +
						other_grid[z + 1][y][x] +
						other_grid[z - 1][y][x];

					grid[z][y][x] = sum_neibours / 6;
				}
			}
		}
		// Exchange grid and other_grid;
		type*** tmp = grid;
		grid = other_grid;
		other_grid = tmp;
	}
	destroy(grid, d, h);
	return other_grid;
}