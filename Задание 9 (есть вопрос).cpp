#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

double normVector(double* v, int n) {
	double result = v[0];

	for (int i = 1; i < n; ++i) {
		if (v[i] > result) {
			result = v[i];
		}
	}

	return result;
}

double* addVectors(double* v1, double* v2, int n) {
	double* result = new double[n];

	for (int i = 0; i < n; ++i) {
		result[i] = v1[i] + v2[i];
	}

	return result;
}

double* subtractVectors(double* v1, double* v2, int n) {
	double* result = new double[n];

	for (int i = 0; i < n; ++i) {
		result[i] = v1[i] - v2[i];
	}

	return result;
}

double* multiplyVecByNum(double* v, int n, double num) {
	double* result = new double[n];

	for (int i = 0; i < n; ++i) {
		result[i] = v[i] * num;
	}

	return result;
}

double normMatrix(double** m, int a, int b, bool clearMemory = false) {
	double result = 0, lineSum;

	for (int i = 0; i < a; ++i) {
		lineSum = 0;

		for (int j = 0; j < b; ++j) {
			lineSum += abs(m[i][j]);
		}

		if (lineSum > result) {
			result = lineSum;
		}
	}

	if (clearMemory) {
		for (int i = 0; i < a; ++i) {
			delete[] m[i];
		}

		delete[] m;
	}

	return result;
}

double** addMatrices(double** m1, double** m2, int a, int b, bool clearM2Memory = false) {
	double** result = new double* [a];

	for (int i = 0; i < a; ++i) {
		result[i] = new double[b];

		for (int j = 0; j < b; ++j) {
			result[i][j] = m1[i][j] + m2[i][j];
		}
	}

	if (clearM2Memory) {
		for (int i = 0; i < a; ++i) {
			delete[] m2[i];
		}

		delete[] m2;
	}

	return result;
}

double** multiplyMatrByNum(double** m, int a, int b, double num) {
	double** result = new double* [a];

	for (int i = 0; i < a; ++i) {
		result[i] = new double[b];

		for (int j = 0; j < b; ++j) {
			result[i][j] = m[i][j] * num;
		}
	}

	return result;
}

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

void printExtendedMatrix(double** matrix, double* vector, int a, int b, string delim = "  ") {
	const int precision = 2;

	cout << fixed << setprecision(precision);

	for (int i = 0; i < a; ++i) {
		for (int j = 0; j < b; ++j) {
			cout << matrix[i][j] << delim;
		}

		cout << "|  " << vector[i] << '\n';
	}

	cout << '\n';
}

void printTable(int* nodes, double* valuesApproximate, double* valuesAccurate, int n) {
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

		cout << setw(LABEL_WIDTH) << "Номер элемента";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << nodes[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "Приближ. решение";
		for (int j = 0; j < numsNumber; ++j) {
			cout << setw(NUM_WIDTH) << valuesApproximate[j + maxNumsInLine * i];
		}

		cout << '\n' << setw(LABEL_WIDTH) << "Точное решение";
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

int main() {
	setlocale(LC_ALL, "RUS");


	int precision = 3;
	double eps = pow(10, -precision);

	const int variant = 16;

	//const int n = 5;
	// 
	//double** A = new double* [n];
	//A[0] = new double[n] {101, -3, -9, 8, 5};
	//A[1] = new double[n] {7, -102, 4, 0, 6};
	//A[2] = new double[n] {6, 2, 103, -5, -10};
	//A[3] = new double[n] {0, -4, 4, 104, 15};
	//A[4] = new double[n] {1, -3, 0, -5, 105};
	//16.38614        14.01961        16.39806        21.92308        18.76190

	random_device rd;
	int minN = 4;
	int maxN = 15;
	int n = minN + rd() % (maxN - minN + 1);

	int maxElementAbs = 20;
	double** A = new double* [n];
	double linePartialSum;
	for (int i = 0; i < n; ++i) {
		A[i] = new double[n];

		linePartialSum = 0;
		for (int j = 0; j < n; j += (j != i - 1) ? 1 : 2) {
			A[i][j] = (rd() % (maxElementAbs + 1)) * pow(-1, rd() % 2);
			linePartialSum += abs(A[i][j]);
		}
		//A[i][i] = rd() % maxElementAbs + 1;
		A[i][i] = linePartialSum + 5 * maxElementAbs + rd() % (20 * maxElementAbs + 1);
	}
	
	double* accurateSolution = new double[n];
	for (int i = 0; i < n; ++i) {
		accurateSolution[i] = variant + i;
	}

	double* b = multiplyMatrByVec(A, n, n, accurateSolution, n);

	cout << "Расширенная матрица коэффициентов СЛАУ:\n";
	printExtendedMatrix(A, b, n, n, "\t");
	cout << '\n';


	double** alpha = new double* [n];
	double* beta = new double[n];
	double* x = new double[n];
	for (int i = 0; i < n; ++i) {
		alpha[i] = new double[n];
		for (int j = 0; j < n; ++j) {
			alpha[i][j] = (-1) * A[i][j] / A[i][i];
		}
		alpha[i][i] = 0;

		beta[i] = b[i] / A[i][i];

		x[i] = 0;
	}

	cout << "Матрица alpha и столбец beta, разделённые вертикальной чертой:\n";
	printExtendedMatrix(alpha, beta, n, n, "\t");


	double alphaMatrNorm = abs(normMatrix(alpha, n, n));

	if (alphaMatrNorm >= 1 - eps) {
		double** identityMatrix = new double* [n];
		for (int i = 0; i < n; ++i) {
			identityMatrix[i] = new double[n] {0};

			for (int j = 0; j < n; ++j) {
				identityMatrix[i][j] = 0;
			}
			
			identityMatrix[i][i] = 1;
		}

		double tau = 1;

		double prevAlphaMatrNorm = alphaMatrNorm;
		alphaMatrNorm = abs(normMatrix(addMatrices(identityMatrix, multiplyMatrByNum(A, n, n, -tau), n, n, true), n, n, true));

		while (alphaMatrNorm >= 1 - eps) {
			tau -= eps;
			alphaMatrNorm = abs(normMatrix(addMatrices(identityMatrix, multiplyMatrByNum(A, n, n, -tau), n, n, true), n, n, true));
		}

		alpha = addMatrices(identityMatrix, multiplyMatrByNum(A, n, n, -tau), n, n);
		beta = multiplyVecByNum(b, n, tau);

		cout << "Была проведена параметризация МПИ, так как норма исходной матрицы alpha больше или равна 1 и составляет:\n" << prevAlphaMatrNorm << '\n';
		cout << "Коэффициент tau параметризации МПИ равен:\n" << tau << "\n\n";
		cout << "Полученные в результате параметризации МПИ матрица alpha и столбец beta, разделённые вертикальной чертой:\n";
		printExtendedMatrix(alpha, beta, n, n, "\t");
		cout << "Норма полученной матрицы alpha:\n" << alphaMatrNorm << "\n\n";
	}
	else {
		cout << "Параметризация МПИ не проводилась, т.к. норма исходной матрицы alpha меньше 1 и составляет:\n" << alphaMatrNorm << "\n\n";
	}


	double* xNext = addVectors(multiplyMatrByVec(alpha, n, n, x, n), beta, n);
	double* xTmp;
	while (normVector(subtractVectors(xNext, x, n), n) > eps) {
		xTmp = x;
		x = xNext;
		xNext = addVectors(multiplyMatrByVec(alpha, n, n, xTmp, n), beta, n);
	}

	int* indices = new int[n];
	for (int i = 0; i < n; ++i) {
		indices[i] = i;
	}

	cout << '\n';
	printTable(indices, xNext, accurateSolution, n);


	return 0;
}