#pragma once
static long num_steps = 100000000;
double step;

double getPi() {
	int i;
	double x, pi, bar, sum = 0.0;
	step = 1.0 / (double)num_steps;

	#pragma omp parallel
	{
		#pragma omp for private(i, x, bar)
		for (i = 0; i < num_steps; i++) {
			x = (i + 0.5) * step;
			#pragma omp critical
			{
				sum = sum + 4.0 / (1.0 + x * x);
			}
		}
	}
	pi = step * sum;

	return pi;
}