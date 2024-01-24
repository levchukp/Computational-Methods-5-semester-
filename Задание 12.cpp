#include <iostream>
#include <iomanip>

using namespace std;

bool directStep(double** extendedMatrix, int n, int step, double eps = 0.00001) {
	double coeff = extendedMatrix[step][step];

	if (-eps < coeff && coeff < eps) {
		int maxInd = -1;
		double maxCoeff = coeff;
		bool inLine = false;

		for (int i = step + 1; i < n; ++i) {
			if (abs(extendedMatrix[i][step]) > abs(maxCoeff)) {
				maxInd = i;
				maxCoeff = extendedMatrix[i][step];
			}
		}

		for (int i = step + 1; i < n; ++i) {
			if (abs(extendedMatrix[step][i]) > abs(maxCoeff)) {
				maxInd = i;
				if (!inLine) inLine = true;
				maxCoeff = extendedMatrix[step][i];
			}
		}

		if (maxInd != -1) {
			if (!inLine) {
				swap(extendedMatrix[step], extendedMatrix[maxInd]);
			}
			else {
				for (int i = 0; i < n; ++i) {
					swap(extendedMatrix[i][step], extendedMatrix[i][maxInd]);
				}
			}
			coeff = maxCoeff;
		}
		else {
			return false;
		}
	}

	for (int i = step; i < n + 1; ++i) {
		extendedMatrix[step][i] /= coeff;
	}

	for (int i = step + 1; i < n; ++i) {
		coeff = extendedMatrix[i][step];

		for (int j = step; j < n + 1; ++j) {
			extendedMatrix[i][j] -= extendedMatrix[step][j] * coeff;
		}
	}

	return true;
}

void printTable(double* nodes, double* valuesApproximate, double* valuesAccurate, int n) {
	const int PRECISION = 3;
	const int WINDOW_WIDTH = 120;
	const int LABEL_WIDTH = 20;
	const int NUM_WIDTH = 13;

	int maxNumsInLine = (WINDOW_WIDTH - LABEL_WIDTH) / NUM_WIDTH;
	int linesNumber = (n + maxNumsInLine - 1) / maxNumsInLine;
	int numsNumber;

	cout << fixed << setprecision(PRECISION) << left;

	for (int i = 0; i < linesNumber; ++i) {
		numsNumber = (n - i * maxNumsInLine >= maxNumsInLine) ? maxNumsInLine : n % maxNumsInLine;

		cout << setw(LABEL_WIDTH) << "Узлы";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << nodes[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "Приближ. значения";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << valuesApproximate[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "Точные значения";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << valuesAccurate[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "Погрешность";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << valuesApproximate[j + maxNumsInLine * i] - valuesAccurate[j + maxNumsInLine * i];
		}

		cout << "\n\n";
	}
}

double p(double x) {
	return x * x;
}

double q(double x) {
	return x;
}

const double variant = 16;

double f(double x, double t) {
	return 4 * variant * pow(x, 4) - 3 * variant * t * pow(x, 3) + 6 * variant * x - 2 * variant * t;
}

double phiOriginal(double x, double t, int k) {
	k++;
	return pow(x, k + 1) * (x - t);
}

double phiFirstDerivative(double x, double t, int k) {
	k++;
	return pow(x, k) * ((x - t) * (k + 1) + x);
}

double phiSecondDerivative(double x, double t, int k) {
	k++;
	return (k + 1) * pow(x, k - 1) * (2 * x + k * (x - t));
}

double yAccurate(double x, double t) {
	return variant * x * x * (x - t);
}

int main() {
	setlocale(LC_ALL, "RUS");


	const int n = 10;

	double t = variant;
	double step = (t - 0) / (n + 1);

	double* nodes = new double[n];
	for (int i = 0; i < n; ++i) {
		nodes[i] = step + i * step;
	}


	double** phiMatrix = new double* [n];
	double* fVector = new double[n];

	for (int i = 0; i < n; ++i) {
		phiMatrix[i] = new double[n];

		for (int j = 0; j < n; ++j) {
			phiMatrix[i][j] = phiSecondDerivative(nodes[i], t, j) + p(nodes[i]) * phiFirstDerivative(nodes[i], t, j) + q(nodes[i]) * phiOriginal(nodes[i], t, j);
		}

		fVector[i] = f(nodes[i], t);
	}


	double** extendedMatrix = new double* [n];
	for (int i = 0; i < n; ++i) {
		extendedMatrix[i] = new double[n + 1];

		for (int j = 0; j < n; ++j) {
			extendedMatrix[i][j] = phiMatrix[i][j];
		}
		extendedMatrix[i][n] = fVector[i];
	}

	bool notReduced;
	for (int i = 0; i < n; ++i) {
		notReduced = directStep(extendedMatrix, n, i);
		if (!notReduced) break;
	}


	if (notReduced) {
		double* a = new double[n];
		for (int i = 0; i < n; ++i) {
			a[n - 1 - i] = extendedMatrix[n - 1 - i][n];

			for (int j = 0; j < i; ++j) {
				a[n - 1 - i] -= extendedMatrix[n - 1 - i][n - 1 - j] * a[n - 1 - j];
			}
		}


		double* valuesApproximate = new double[n];
		for (int i = 0; i < n; ++i) {
			valuesApproximate[i] = 0;

			for (int j = 0; j < n; ++j) {
				valuesApproximate[i] += a[j] * phiOriginal(nodes[i], t, j);
			}
		}


		double* valuesAccurate = new double[n];
		for (int i = 0; i < n; ++i) {
			valuesAccurate[i] = yAccurate(nodes[i], t);
		}

		printTable(nodes, valuesApproximate, valuesAccurate, n);
	}
	else {
		cout << "Матрица вырожденная, поэтому решение уравнения найти невозможно\n";
	}


	return 0;
}