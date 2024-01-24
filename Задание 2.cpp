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

            result += numbers[0][i] * pow(-1, i) * det(minor, n - 1);
        }
    }

    return result;
}

double polynomial(double x, double* coefficients, int n) {
    double result = 0.;

    for (int i = 0; i < n; ++i) {
        result += coefficients[n - 1 - i] * pow(x, i);
    }

    return result;
}

void printArray(double* array, int n) {
    for (int i = 0; i < n; ++i) {
        cout << array[i] << ' ';
    }

    cout << "\n\n";
}

void printPolynomial(double* coefficients, int n) {
    const int PRECISION = 3;

    cout << fixed << setprecision(PRECISION) << left;

    for (int i = n - 1; i > 0; --i) {
        cout << coefficients[n - 1 - i] << " * x^" << i << " + ";
    }

    cout << coefficients[n - 1] << " = P(x)\n\n";
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


    double** matrix = new double* [NODES];
    for (int i = 0; i < NODES; ++i) {
        matrix[i] = new double[NODES];

        for (int j = 0; j < NODES; ++j) {
            matrix[i][j] = pow(nodes[i], NODES - 1 - j);
        }
    }
    double mainDet = det(matrix, NODES);

    double* fColumn = new double[NODES];
    for (int i = 0; i < NODES; ++i) {
        fColumn[i] = values[i];
    }

    double coefficients[NODES];
    double** fMatrix = new double* [NODES];
    for (int i = 0; i < NODES; ++i) {
        for (int j = 0; j < NODES; ++j) {
            fMatrix[j] = new double[NODES];

            for (int k = 0; k < NODES; ++k) { // возможно, стоит проверять, верно ли k == i
                fMatrix[j][k] = matrix[j][k];
            }

            fMatrix[j][i] = fColumn[j];
        }

        coefficients[i] = det(fMatrix, NODES) / mainDet;
    }


    double newNodes[2 * NODES - 1] = { nodes[0], (nodes[0] + nodes[1]) / 2., // возможно, расчёт промежуточных узлов лучше автоматизировать
                                      nodes[1], (nodes[1] + nodes[2]) / 2.,
                                      nodes[2], (nodes[2] + nodes[3]) / 2.,
                                      nodes[3] };

    double newValues[2 * NODES - 1];

    for (int i = 0; i < 2 * NODES - 1; ++i) {
        newValues[i] = polynomial(newNodes[i], coefficients, NODES);
    }


    cout << "Таблица значений интерполируемой функции в исходных узлах:\n";
    printTable(nodes, values, NODES);

    cout << "СЛАУ, устанавливающая удовлетворение интерполяционным многочленом ГУИ:\n";
    for (int i = 0; i < NODES; ++i) {
        for (int j = 0; j < NODES; ++j) {
            cout << matrix[i][j] << "  ";
        }
        cout << "|  " << fColumn[i] << '\n';
    }
    cout << '\n';

    cout << "Интерполяционный многочлен:\n";
    printPolynomial(coefficients, NODES);
    cout << "Таблица значений интерполирующей функции в исходных узлах и узлах интерполяции:\n";
    printTable(newNodes, newValues, 2 * NODES - 1);


    system("pause");
    return 0;
}