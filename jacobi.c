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
	for (int x = 0; x < w; x++) {
		for (int y = 0; y < h; y++) {
			for (int z = 0; z < d; z++) {
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
		for (int x = 0; x < w; x++) {
			for (int y = 0; y < h; y++) {
				for (int z = 0; z < d; z++) {
					// grid[x][y][z] = avg of neighbors in other_grid*/
					int sum_neibours =
						other_grid[z][y][min(x + 1, w-1)] +
						other_grid[z][y][max(x - 1, 0)] +
						other_grid[z][min(y + 1, h-1)][x] +
						other_grid[z][max(y - 1, 0)][x] +
						other_grid[min(z + 1, d-1)][y][x] +
						other_grid[max(z - 1, 0)][y][x];

					grid[z][y][x] = sum_neibours / 6;
				}
			}
		}
		// Exchange grid and other_grid;
		int*** tmp = create(d, h, w);
		copy(tmp, grid, d, h, w);
		copy(grid, other_grid, d, h, w);
		copy(other_grid, tmp, d, h, w);
	}
	destroy(grid, d, h);
	return other_grid;
}