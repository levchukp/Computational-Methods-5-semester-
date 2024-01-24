#include <iostream>
#include <iomanip>

using namespace std;

void printTable(double* nodes, double* valuesApproximate, double* valuesAccurate, int n) {
	const int PRECISION = 3;
	const int WINDOW_WIDTH = 120;
	const int LABEL_WIDTH = 20;
	const int NUM_WIDTH = 10;

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

double yAccurate(double x, double t) {
	return variant * x * x * (x - t);
}

int main() {
	setlocale(LC_ALL, "RUS");

	
	const int n = 100;

	double t = variant;
	double step = (t - 0) / (n - 1);

	double* nodes = new double[n];
	for (int i = 0; i < n; ++i) {
		nodes[i] = 0 + i * step;
	}


	double stepInvertSquare = 1 / step / step;

	double** A = new double* [n];
	for (int i = 1; i < n - 1; ++i) {
		A[i] = new double[n];

		for (int j = 0; j < i - 1; ++j) {
			A[i][j] = 0;
		}

		A[i][i - 1] = stepInvertSquare - p(nodes[i]) / 2 / step;
		A[i][i] = -2 * stepInvertSquare + q(nodes[i]);
		A[i][i + 1] = stepInvertSquare + p(nodes[i]) / 2 / step;

		for (int j = i + 2; j < n; ++j) {
			A[i][j] = 0;
		}
	}

	A[0] = new double[n];
	A[n - 1] = new double[n];
	for (int i = 0; i < n; ++i) {
		A[0][i] = 0;
		A[n - 1][i] = 0;
	}
	A[0][0] = 1;
	A[n - 1][n - 1] = 1;
	
	double* b = new double[n];
	for (int i = 1; i < n - 1; ++i) {
		b[i] = f(nodes[i], t);
	}
	b[0] = 0;
	b[n - 1] = 0;


	double* P = new double[n];
	double* Q = new double[n];

	P[0] = A[0][1] / A[0][0] * (-1);
	Q[0] = -b[0] / A[0][0] * (-1);

	for (int i = 1; i < n; ++i) {
		P[i] = A[i][i + 1] / (A[i][i] * (-1) - A[i][i - 1] * P[i - 1]);
		Q[i] = (A[i][i - 1] * Q[i - 1] - b[i]) / (A[i][i] * (-1) - A[i][i - 1] * P[i - 1]);
	}


	double* y = new double[n];

	y[n - 1] = Q[n - 1];

	for (int i = n - 2; i >= 0; --i) {
		y[i] = P[i] * y[i + 1] + Q[i];
	}


	double* valuesAccurate = new double[n];
	for (int i = 0; i < n; ++i) {
		valuesAccurate[i] = yAccurate(nodes[i], t);
	}

	printTable(nodes, y, valuesAccurate, n);


	return 0;
}