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

void printMatrix(double** matrix, int a, int b, string delim = " ") {
	const int precision = 2;

	cout << fixed << setprecision(precision);

	for (int i = 0; i < a; ++i) {
		for (int j = 0; j < b; ++j) {
			cout << matrix[i][j] << delim;
		}

		cout << '\n';
	}

	cout << '\n';
}

int main() {
	setlocale(LC_ALL, "RUS");

	const double variant = 16;
	const int n = 5;

	double** powCoeffs = new double*[n];
	double* constCoeffs = new double[n];

	double* diagElems = new double[n];
	double nondiagElem;

	for (int i = 0; i < n; ++i) {
		powCoeffs[i] = new double[n];
		
		diagElems[i] = variant + i;
		nondiagElem = diagElems[i] * 0.01;

		for (int j = 0; j < n; ++j) {
			powCoeffs[i][j] = (i == j) ? diagElems[i] : nondiagElem;
		}
	}
	
	for (int i = 0; i < n; ++i) {
		constCoeffs[i] = 0.;

		for (int j = 0; j < n; ++j) {
			constCoeffs[i] += powCoeffs[i][j] * diagElems[j];
		}
	}

	double** extendedMatrix = new double*[n];
	for (int i = 0; i < n; ++i) {
		extendedMatrix[i] = new double[n + 1];

		for (int j = 0; j < n; ++j) {
			extendedMatrix[i][j] = powCoeffs[i][j];
		}
		extendedMatrix[i][n] = constCoeffs[i];
	}
	
	cout << "Расширенная матрица:\n";
	printMatrix(extendedMatrix, n, n + 1, "  ");

	
	bool notReduced;
	for (int i = 0; i < n; ++i) {
		notReduced = directStep(extendedMatrix, n, i);
		if (!notReduced) break;
	}

	if (notReduced) {
		cout << "Расширенная матрица после окончания прямого хода:\n";
		printMatrix(extendedMatrix, n, n + 1, "  ");

		double* x = new double[n];
		for (int i = 0; i < n; ++i) {
			x[i] = extendedMatrix[n - 1 - i][n];

			for (int j = 0; j < i; ++j) {
				x[i] -= extendedMatrix[n - 1 - i][n - 1 - j] * x[j];
			}
		}

		cout << "Решение матричного уравнения:\n";
		for (int i = n - 1; i >= 0; --i) {
			cout << x[i] << "  ";
		}
		cout << '\n';
	}
	else {
		cout << "Матрица вырожденная, поэтому решение уравнения найти невозможно\n";
	}

	return 0;
}