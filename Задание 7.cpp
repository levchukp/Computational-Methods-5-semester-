#include <iostream>
#include <iomanip>

using namespace std;

bool directStep(double** matrix, int n, int step, double eps = 0.00001, bool accumulateChanges = false, double* detAbs = NULL, int* signChanges = NULL) {
	double coeff = matrix[step][step];

	if (-eps < coeff && coeff < eps) {
		int maxInd = -1;
		double maxCoeff = coeff;
		bool inLine = false;

		for (int i = step + 1; i < n; ++i) {
			if (abs(matrix[i][step]) > abs(maxCoeff)) {
				maxInd = i;
				maxCoeff = matrix[i][step];
			}
		}

		for (int i = step + 1; i < n; ++i) {
			if (abs(matrix[step][i]) > abs(maxCoeff)) {
				maxInd = i;
				if (!inLine) inLine = true;
				maxCoeff = matrix[step][i];
			}
		}

		if (maxInd != -1) {
			if (!inLine) {
				swap(matrix[step], matrix[maxInd]);
			}
			else {
				for (int i = 0; i < n; ++i) {
					swap(matrix[i][step], matrix[i][maxInd]);
				}
			}
			coeff = maxCoeff;
			if (accumulateChanges) *signChanges += 1;
		}
		else {
			return false;
		}
	}

	for (int i = step; i < n + (!accumulateChanges) ? 1 : 0; ++i) {
		matrix[step][i] /= coeff;
	}

	if (accumulateChanges) *detAbs *= coeff;

	for (int i = step + 1; i < n; ++i) {
		coeff = matrix[i][step];

		for (int j = step; j < n + (!accumulateChanges) ? 1 : 0; ++j) {
			matrix[i][j] -= matrix[step][j] * coeff;
		}
	}

	return true;
}

double detMinorsMethod(double** numbers, int n) {
	double result = 0.;

	if (n == 2) {
		result = numbers[0][0] * numbers[1][1] - numbers[0][1] * numbers[1][0];
	}

	else {
		double** minor = new double* [n - 1];

		for (int i = 0; i < n; ++i) {
			for (int j = 1; j < n; ++j) {
				minor[j - 1] = new double[n - 1];

				for (int k = 0; k < i; ++k) {
					minor[j - 1][k] = numbers[j][k];
				}

				for (int k = i + 1; k < n; ++k) {
					minor[j - 1][k - 1] = numbers[j][k];
				}
			}

			result += numbers[0][i] * pow(-1, i) * detMinorsMethod(minor, n - 1);
		}
	}

	return result;
}

double** multiplyMatrices(double** m1, int a1, int b1, double** m2, int a2, int b2) {
	if (b1 != a2) {
		return NULL;
	}

	double** result = new double* [a1];

	for (int i = 0; i < a1; ++i) {
		result[i] = new double[b2];

		for (int j = 0; j < b2; ++j) {
			result[i][j] = 0;

			for (int k = 0; k < b1; ++k) {
				result[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}

	return result;
}

bool matrixIsEqual(double** m1, double** m2, int a, int b, double eps = 0.00001) {
	for (int i = 0; i < a; ++i) {
		for (int j = 0; j < b; ++j) {
			if (abs(m2[i][j] - m1[i][j]) >= eps) {
				return false;
			}
		}
	}

	return true;
}

void printMatrix(double** matrix, int a, int b, string delim = "  ") {
	const int precision = 3;

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


	const int n = 4;

	double** A = new double* [n];
	A[0] = new double[n] {-2, 1, 3, 4};
	A[1] = new double[n] {3, 0, -1, 4};
	A[2] = new double[n] {-5, 2, 3, 0};
	A[3] = new double[n] {4, -1, 2, -6};

	cout << "ƒанна€ матрица:\n";
	printMatrix(A, n, n, "\t");
	cout << "\n";


	double** matrix = new double* [n];
	for (int i = 0; i < n; ++i) {
		matrix[i] = new double[n];

		for (int j = 0; j < n; ++j) {
			matrix[i][j] = A[i][j];
		}
	}

	double detAbs = 1;
	int signChanges = 0;

	bool notReduced;
	for (int i = 0; i < n; ++i) {
		notReduced = directStep(matrix, n, i, 0.00001, true, &detAbs, &signChanges);
		if (!notReduced) break;
	}


	if (notReduced) {
		double det = detAbs * pow(-1, signChanges % 2);

		cout << "ќпределитель данной матрицы, найденный с помощью пр€мого хода метода решени€ матричных уравнений √аусса, равен:\n" << det << "\n";
		double detFromMinors = detMinorsMethod(A, n);
		cout << "ƒл€ проверки полученного результата найдЄм определитель матрицы с использованием миноров. ќн равен:\n" << detFromMinors << "\n";
		if (det == detFromMinors) {
			cout << "ѕолученные значени€ равны, из чего следует, что, скорее всего, определитель с помощью пр€мого хода найден правильно.\n\n";
		}
		else {
			cout << "ѕолученные значени€ не равны, из чего следует, что, возможно, определитель с помощью пр€мого хода найден неправильно.\n\n";
		}
		cout << '\n';


		double** invertMatrix = new double* [n];
		double** extendedMatrix = new double* [n];
		for (int i = 0; i < n; ++i) {
			invertMatrix[i] = new double[n];
			extendedMatrix[i] = new double[n + 1];
		}

		for (int k = 0; k < n; ++k) {
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					extendedMatrix[i][j] = A[i][j];
				}
				extendedMatrix[i][n] = 0;
			}
			extendedMatrix[k][n] = 1;

			for (int i = 0; i < n; ++i) {
				directStep(extendedMatrix, n, i);
			}
			
			for (int i = 0; i < n; ++i) {
				invertMatrix[n - 1 - i][k] = extendedMatrix[n - 1 - i][n];

				for (int j = 0; j < i; ++j) {
					invertMatrix[n - 1 - i][k] -= extendedMatrix[n - 1 - i][n - 1 - j] * invertMatrix[n - 1 - j][k];
				}
			}
		}

		cout << "ћатрица, обратна€ к данной:\n";
		printMatrix(invertMatrix, n, n, "\t");
		

		cout << "ѕроверим правильность нахождени€ обратной матрицы, умножив исходную матрицу на полученную обратную матрицу слева и справа.\n\n";

		double** identityMatrix = new double* [n];
		for (int i = 0; i < n; ++i) {
			identityMatrix[i] = new double[n] {0};
			identityMatrix[i][i] = 1;
		}
		
		bool leftMultIsIdentityMatrix, rightMultIsIdentityMatrix;

		cout << "–езультат умножени€ исходной матрицы на найденную обратную:\n";
		double** checkMatrix = multiplyMatrices(A, n, n, invertMatrix, n, n);
		printMatrix(checkMatrix, n, n, "\t");
		leftMultIsIdentityMatrix = matrixIsEqual(checkMatrix, identityMatrix, n, n);

		cout << "–езультат умножени€ найденной обратной матрицы на исходную:\n";
		checkMatrix = multiplyMatrices(invertMatrix, n, n, A, n, n);
		printMatrix(checkMatrix, n, n, "\t");
		rightMultIsIdentityMatrix = matrixIsEqual(checkMatrix, identityMatrix, n, n);

		if (leftMultIsIdentityMatrix && rightMultIsIdentityMatrix) {
			cout << "¬ результате обоих умножений получаетс€ единичина€ матрица, значит, обратна€ матрица найдена правильно.\n";
		}
		else {
			cout << "¬ результате одного или обоих умножений не получаетс€ единична€ матрица, значит, обратна€ матрица найдена неправильно.\n";
		}
	}
	else {
		cout << "ќпределитель матрицы равен " << 0 << "\n\n";
		cout << "ћатрица вырожденна€, поэтому обратную к ней найти невозможно\n";
	}


	return 0;
}