#include <iostream>
#include <iomanip>

using namespace std;

double* multiplyMatrByVec(double** m, int a, int b, double* v, int n) {
	if (b != n) {
		return NULL;
	}

	double* result = new double[a];

	for (int i = 0; i < a; ++i) {
		result[i] = 0;

		for (int j = 0; j < b; ++j) {
			result[i] += m[i][j] * v[j];
		}
	}

	return result;
}

void printExtendedMatrix(double** extendedMatrix, int a, int b, string delim = "  ") {
	const int precision = 2;

	cout << fixed << setprecision(precision);

	for (int i = 0; i < a; ++i) {
		for (int j = 0; j < b; ++j) {
			cout << extendedMatrix[i][j] << delim;
		}
		cout << "|    " << extendedMatrix[i][b];

		cout << '\n';
	}

	cout << '\n';
}

void printTable(int* i_s, double* P, double* Q, int n) {
	const int PRECISION = 4;
	const int WINDOW_WIDTH = 120;
	const int LABEL_WIDTH = 8;
	const int NUM_WIDTH = 10;

	int maxNumsInLine = (WINDOW_WIDTH - LABEL_WIDTH) / NUM_WIDTH;
	int linesNumber = (n + maxNumsInLine - 1) / maxNumsInLine;
	int numsNumber;

	cout << fixed << setprecision(PRECISION) << left;

	for (int i = 0; i < linesNumber; ++i) {
		numsNumber = (n - i * maxNumsInLine >= maxNumsInLine) ? maxNumsInLine : n % maxNumsInLine;

		cout << setw(LABEL_WIDTH) << "i";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << i_s[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "P_i";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << P[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "Q_i";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << Q[j + maxNumsInLine * i];
		}

		cout << "\n\n";
	}
}

int main() {
	setlocale(LC_ALL, "RUS");


	const int n = 5;

	double** A = new double* [n];
	A[0] = new double[n] {-2, -3, 0, 0, 0};
	A[1] = new double[n] {7, -8, 4, 0, 0};
	A[2] = new double[n] {0, 2, 3, -5, 0};
	A[3] = new double[n] {0, 0, 4, -12, 15};
	A[4] = new double[n] {0, 0, 0, -1, 5};

	const int variant = 16;

	double* d = new double[n] {variant, variant + 1, variant + 2, variant + 3, variant + 4};
	d = multiplyMatrByVec(A, n, n, d, n);

	
	double* P = new double[n];
	double* Q = new double[n];

	P[0] = A[0][1] / A[0][0] * -1;
	Q[0] = -d[0] / A[0][0] * -1;

	for (int i = 1; i < n; ++i) {
		P[i] = A[i][i + 1] / (A[i][i] * -1 - A[i][i - 1] * P[i - 1]);
		Q[i] = (A[i][i - 1] * Q[i - 1] - d[i]) / (A[i][i] * -1 - A[i][i - 1] * P[i - 1]);
	}


	double* x = new double[n];

	x[n - 1] = Q[n - 1];

	for (int i = n - 2; i >= 0; --i) {
		x[i] = P[i] * x[i + 1] + Q[i];
	}

	
	double** extendedMatrix = new double* [n];
	for (int i = 0; i < n; ++i) {
		extendedMatrix[i] = new double[n + 1];

		for (int j = 0; j < n; ++j) {
			extendedMatrix[i][j] = A[i][j];
		}
		extendedMatrix[i][n] = d[i];
	}

	cout << "Расширенная матрица коэффициентов СЛАУ, появляющейся при решении данного матричного уравнения:\n";
	printExtendedMatrix(extendedMatrix, n, n, "\t");

	int* i_s = new int[n];
	for (int i = 0; i < n; ++i) {
		i_s[i] = i + 2;
	}

	cout << "Коэффициенты P и Q, полученные при прямой прогонке:\n\n";
	printTable(i_s, P, Q, n);

	cout << "Решение матричного уравнения (порядок элементов естественный):\n";
	for (int i = 0; i < n; ++i) {
		cout << x[i] << "  ";
	}
	cout << '\n';


	return 0;
}