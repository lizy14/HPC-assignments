

#include "jacobi.h"

#define D 128
#define D_X D
#define D_Y D
#define D_Z D
#define N_STEPS 100
#define N_SAMPLES 6

int evaluate() {
	type*** m = create(D_Z, D_Y, D_X);
	randomize(m, D_Z, D_Y, D_X, 0);

	clock_t start = clock();	
	jacobi(m, D_Z, D_Y, D_X, N_STEPS);
	clock_t diff = clock() - start;

	int msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time eclapsed: %d milliseconds\n", msec);

	return msec;
}


int main() {
	int samples[N_SAMPLES];
	for (int i = 0; i < N_SAMPLES; i++) {
		samples[i] = evaluate();
	}
	int total = 0;
	for (int i = 0; i < N_SAMPLES; i++) {
		total += samples[i];
	}
	double average = (total + 0.) / N_SAMPLES;
	double total_square = 0;
	for (int i = 0; i < N_SAMPLES; i++) {
		double delta = samples[i] - average;
		total_square += delta * delta;
	}
	double variance = (total_square / N_SAMPLES);
	printf("average %lf, variance %lf\n", average, variance);

	return 0;
}