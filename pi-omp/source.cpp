#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Timing;
int main();


/**************************************************/

#include "header.h"

/**************************************************/



class Timing
{
private:
	double clockStart, clockEnd;
	vector<double> timings;
	inline double clock() {
		return omp_get_wtime();
	}
public:
	Timing() {
		clockStart = clock();
		clockEnd = clockStart;
	}
	~Timing() {}

public:
	void start() {
		clockStart = clock();
	}
	void end() {
		clockEnd = clock();
		double timeDelta = clockEnd - clockStart;
		timings.push_back(timeDelta);
	}
	void reportCPUtime(bool report_number_of_samples = false) {
		auto n = timings.size();
		double sum = 0;
		for (auto t : timings) {
			sum += t;
		}
		double average = sum / n;

		sum = 0; // sum of squared residue
		for (auto t : timings) {
			auto u = t - average;
			sum += u * u;
		}
		double stddev = sqrt(sum / n);

#ifndef NON_VERBOSE
		cout << "Elapsed wall clock time " << endl;
		if (report_number_of_samples) {
			cout << "number of samples: " << n << endl;
		};
		cout << "average: " << average << "s" << endl;
		if (n != 1) {
			cout << "sample standard deviation: " << stddev << "s" << endl;
		}
#else
        cout << average << '\t' << stddev << endl;
#endif
	}
};




int main() {
	double pi;
	for (auto number_of_threads : vector<int>({ 1,2,4,8,16,32 })) {
		omp_set_num_threads(number_of_threads);
#ifndef NON_VERBOSE
		cout << endl << "number of threads: " << number_of_threads << endl;
#else
        cout<< number_of_threads << '\t';
#endif
		auto t = Timing();

		for (int i = 0; i < 6; i++) {
			t.start();
			pi = getPi();
			t.end();
		}

		t.reportCPUtime();
	}




	cout.precision(20);
	cout << pi << endl;


	return std::system("pause");
}
