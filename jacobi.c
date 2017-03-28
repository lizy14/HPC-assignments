#include "jacobi.h"

int*** create(int d, int h, int w) {
	int*** m = (int***)malloc(d * sizeof(int**));
	for (int z = 0; z < d; z++) {
		m[z] = (int**)malloc(h * sizeof(int*));
		for (int y = 0; y < w; y++) {
			m[z][y] = (int*)malloc(w * sizeof(int));
		}
	}
	return m;
};

int*** copy(int*** dst, int*** src, int d, int h, int w) {
	for (int z = 0; z < d; z++) {
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				dst[z][y][x] = src[z][y][x];
			}
		}
	}
	return dst;
}

int destroy(int*** m, int d, int h) {
	for (int z = 0; z < d; z++) {
		for (int y = 0; y < h; y++) {
			free(m[z][y]);
		}
		free(m[z]);
	}
	return 0;
}

int*** randomize(int*** m, int d, int h, int w, int seed) {
	srand(seed);
	for (int z = 0; z < d; z++) {
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				m[z][y][x] = rand();
			}
		}
	}
	return m;
}

int*** jacobi(int*** original, int d, int h, int w, int step) {

	int*** extra = create(d, h, w);
	int*** grid = extra;
	int*** other_grid = original;

	for (int t = 0; t < step; t++) {		
		for (int z = 1; z < d - 1; z++) {
			for (int y = 1; y < h - 1; y++) {
				for (int x = 1; x < w - 1; x++) {
					// grid[x][y][z] = avg of neighbors in other_grid*/
					int sum_neibours =
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
		int*** tmp = grid;
		grid = other_grid;
		other_grid = tmp;
	}
	destroy(grid, d, h);
	return other_grid;
}