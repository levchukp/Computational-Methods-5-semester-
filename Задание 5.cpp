#include <iostream>
#include <iomanip>

using namespace std;

double det(double** numbers, int n) {
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

            if (numbers[0][i] != 0) {
                result += numbers[0][i] * pow(-1, i) * det(minor, n - 1);
            }
        }
    }

    return result;
}

double splineValue(double x, int splineNumber, double* coefficients, double* splineBorders) {
    double xDistance = x - splineBorders[splineNumber];
    return coefficients[4 * splineNumber]
        + coefficients[4 * splineNumber + 1] * xDistance
        + coefficients[4 * splineNumber + 2] * xDistance * xDistance
        + coefficients[4 * splineNumber + 3] * pow(xDistance, 3);
}

void printTable(double* nodes, double* values, int n) {
    const int PRECISION = 3;
    const int WINDOW_WIDTH = 120;
    const int LABEL_WIDTH = 14;
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

        cout << '\n' << setw(LABEL_WIDTH) << "Значения";
        for (int j = 0; j < numsNumber; ++j) {
            cout << setw(NUM_WIDTH) << values[j + maxNumsInLine * i];
        }

        cout << "\n\n";
    }
}

int main()
{
    setlocale(LC_ALL, "RUS");


    const int NODES = 4;

    double nodes[NODES] = { -2., -1., 0., 2. };
    double values[NODES] = { -8., -1., 0., 8. };


    double steps[NODES - 1], step_squares[NODES - 1];
    for (int i = 1; i < NODES; ++i) {
        steps[i - 1] = nodes[i] - nodes[i - 1];
        step_squares[i - 1] = steps[i - 1] * steps[i - 1];
    }

    const int SPLINES = NODES - 1;

    double** matrix = new double* [4 * SPLINES];
    for (int i = 0; i < 4 * SPLINES; ++i) {
        matrix[i] = new double[4 * SPLINES];

        for (int j = 0; j < 4 * SPLINES; ++j) {
            matrix[i][j] = 0;
        }
    }

    double* fColumn = new double[4 * SPLINES];
    for (int i = 0; i < 4 * SPLINES; ++i) {
        fColumn[i] = 0;
    }


    matrix[0][4 * 0 + 2] = 2;

    matrix[1][4 * (SPLINES - 1) + 2] = 2;
    matrix[1][4 * (SPLINES - 1) + 3] = 6 * steps[SPLINES - 1];

    for (int i = 0; i < SPLINES - 1; ++i) {
        matrix[2 + i * 2][4 * i + 1] = 1;
        matrix[2 + i * 2][4 * i + 2] = 2 * steps[i];
        matrix[2 + i * 2][4 * i + 3] = 3 * step_squares[i];
        matrix[2 + i * 2][4 * (i + 1) + 1] = -1;

        matrix[3 + i * 2][4 * i + 2] = 2;
        matrix[3 + i * 2][4 * i + 3] = 6 * steps[i];
        matrix[3 + i * 2][4 * (i + 1) + 2] = -2;
    }

    for (int i = 0; i < SPLINES; ++i) {
        matrix[2 * SPLINES + i * 2][4 * i] = 1;
        fColumn[2 * SPLINES + i * 2] = values[i];

        matrix[2 * SPLINES + i * 2 + 1][4 * i] = 1;
        matrix[2 * SPLINES + i * 2 + 1][4 * i + 1] = steps[i];
        matrix[2 * SPLINES + i * 2 + 1][4 * i + 2] = step_squares[i];
        matrix[2 * SPLINES + i * 2 + 1][4 * i + 3] = step_squares[i] * steps[i];
        fColumn[2 * SPLINES + i * 2 + 1] = values[i + 1];
    }

    double mainDet = det(matrix, 4 * SPLINES);


    double coefficients[4 * SPLINES];
    double** fMatrix = new double* [4 * SPLINES];
    for (int i = 0; i < 4 * SPLINES; ++i) {
        for (int j = 0; j < 4 * SPLINES; ++j) {
            fMatrix[j] = new double[4 * SPLINES];

            for (int k = 0; k < 4 * SPLINES; ++k) {
                fMatrix[j][k] = matrix[j][k];
            }

            fMatrix[j][i] = fColumn[j];
        }

        coefficients[i] = det(fMatrix, 4 * SPLINES) / mainDet;
    }


    double newNodes[2 * NODES - 1];
    for (int i = 0; i < NODES - 1; ++i) {
        newNodes[i * 2] = nodes[i];
        newNodes[1 + i * 2] = (nodes[i] + nodes[i + 1]) / 2;
    }
    newNodes[2 * NODES - 2] = nodes[NODES - 1];

    double newValues[2 * NODES - 1];
    for (int i = 0; i < 2 * NODES - 2; ++i) {
        newValues[i] = splineValue(newNodes[i], i / 2, coefficients, nodes);
    }
    newValues[2 * NODES - 2] = splineValue(newNodes[2 * NODES - 2], NODES - 2, coefficients, nodes);


    cout << "Таблица значений интерполируемой функции в исходных узлах:\n";
    printTable(nodes, values, NODES);

    cout << "Расширенная матрица коэффициентов для СЛАУ, уравнения которой определяют коэффициенты сплайнов:\n";
    for (int i = 0; i < 4 * SPLINES; ++i) {
        for (int j = 0; j < 4 * SPLINES; ++j) {
            cout << matrix[i][j] << '\t';
        }
        cout << "|  " << fColumn[i] << '\n';
    }
    cout << '\n';

    cout << "Полученные кубические сплайны:\n";
    for (int i = 0; i < SPLINES; ++i) {
        cout << "S_" << i + 1 << " = ";
        cout << "(" << coefficients[4 * i] << ") + ";
        cout << "(" << coefficients[4 * i + 1] << ") * [x - (" << nodes[i] << ")] + ";
        cout << "(" << coefficients[4 * i + 2] << ") * [x - (" << nodes[i] << ")]^2 + ";
        cout << "(" << coefficients[4 * i + 3] << ") * [x - (" << nodes[i] << ")]^3\n";
    }
    cout << '\n';

    cout << "Таблица значений кусочно-непрерывной склейки кубических сплайнов в исходных узлах и узлах интерполяции:\n";
    printTable(newNodes, newValues, 2 * NODES - 1);


    system("pause");
    return 0;
}