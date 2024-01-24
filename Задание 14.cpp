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

void printExtendedMatrix(double** extendedMatrix, int a, int b, string delim = " ") {
	const int precision = 3;

	cout << fixed << setprecision(precision);

	for (int i = 0; i < a; ++i) {
		for (int j = 0; j < b; ++j) {
			cout << extendedMatrix[i][j] << delim;
		}
		cout << '|' << delim << extendedMatrix[i][b];

		cout << '\n';
	}

	cout << '\n';
}

void printTable(double* nodes, double* valuesApproximate, double* valuesAccurate, int n) {
	const int PRECISION = 5;
	const int WINDOW_WIDTH = 120;
	const int LABEL_WIDTH = 20;
	const int NUM_WIDTH = 10;

	int maxNumsInLine = (WINDOW_WIDTH - LABEL_WIDTH) / NUM_WIDTH;
	int linesNumber = (n + maxNumsInLine - 1) / maxNumsInLine;
	int numsNumber;

	cout << fixed << setprecision(PRECISION) << left;

	for (int i = 0; i < linesNumber; ++i) {
		numsNumber = (n - i * maxNumsInLine >= maxNumsInLine) ? maxNumsInLine : n % maxNumsInLine;

		cout << setw(LABEL_WIDTH) << "����";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << nodes[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "�������. ��������";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << valuesApproximate[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "������ ��������";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << valuesAccurate[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "�����������";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << valuesApproximate[j + maxNumsInLine * i] - valuesAccurate[j + maxNumsInLine * i];
		}

		cout << "\n\n";
	}
}

const double variant = 16;

double A(double x, double t) {
	double tmp = x * t;
	return tmp + tmp * tmp + pow(tmp, 3);
}

double f(double x) {
	return variant * (4. / 3. * x + x * x / 4. + pow(x, 3) / 5.);
}

double y(double x) {
	return variant * x;
}

int main() {
	setlocale(LC_ALL, "RUS");
	

	const int n = 30;

	double lambda = 1, a = 0, b = 1;
	
	double step = (b - a) / (n - 1);
	double* nodes = new double[n];
	for (int i = 0; i < n; ++i) {
		nodes[i] = a + step * i;
	}


	double** extendedMatrix = new double* [n];
	for (int i = 0; i < n; ++i) {
		extendedMatrix[i] = new double[n + 1];

		for (int j = 0; j < n; ++j) {
			extendedMatrix[i][j] = lambda * step * A(nodes[i], nodes[j]);
		}
		extendedMatrix[i][i] += 1;

		extendedMatrix[i][n] = f(nodes[i]);
	}

	cout << "����������� ������� ������������� ����:\n";
	printExtendedMatrix(extendedMatrix, n, n, "  ");

	bool notReduced;
	for (int i = 0; i < n; ++i) {
		notReduced = directStep(extendedMatrix, n, i);
		if (!notReduced) break;
	}


	if (notReduced) {
		double* valuesApproximate = new double[n];
		double* valuesAccurate = new double[n];
		for (int i = 0; i < n; ++i) {
			valuesApproximate[n - 1 - i] = extendedMatrix[n - 1 - i][n];

			for (int j = 0; j < i; ++j) {
				valuesApproximate[n - 1 - i] -= extendedMatrix[n - 1 - i][n - 1 - j] * valuesApproximate[n - 1 - j];
			}

			valuesAccurate[i] = y(nodes[i]);
		}

		printTable(nodes, valuesApproximate, valuesAccurate, n);
	}
	else {
		cout << "������� �����������, ������� ������� ��������� ����� ����������\n";
	}


	return 0;
}