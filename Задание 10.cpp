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

double derivative(double x, double yPrev, double parameter) {
    return 2 * parameter * x + parameter * x * x - yPrev;
}

int main() {
	setlocale(LC_ALL, "RUS");

    const int n = 30;
    double variant = 16;
    double x_0 = 1;
    double y_0 = variant;
    double step = 0.1;

    double nodes[n];
    for (int i = 0; i < n; ++i) {
        nodes[i] = x_0 + step * i;
    }

    double valuesAccurate[n];
    for (int i = 0; i < n; ++i) {
        valuesAccurate[i] = variant * nodes[i] * nodes[i];
    }

    double valuesApproximate[n];
    valuesApproximate[0] = y_0;


    for (int i = 1; i < n; ++i) {
        valuesApproximate[i] = valuesApproximate[i - 1] + step * derivative(nodes[i - 1], valuesApproximate[i - 1], variant);
    }

    cout << "Вычисления, проведённые с помощью метода Эйлера:\n\n";
    printTable(nodes, valuesApproximate, valuesAccurate, n);
    cout << '\n';


    double intermediate;
    double stepHalf = step / 2.;

    for (int i = 1; i < n; ++i) {
        intermediate = valuesApproximate[i - 1] + stepHalf * derivative(nodes[i - 1], valuesApproximate[i - 1], variant);
        valuesApproximate[i] = valuesApproximate[i - 1] + step * derivative(nodes[i - 1] + stepHalf, intermediate, variant);
    }

    cout << "Вычисления, проведённые с помощью усовершенствованного метода Эйлера:\n\n";
    printTable(nodes, valuesApproximate, valuesAccurate, n);
    cout << '\n';


    for (int i = 1; i < n; ++i) {
        intermediate = valuesApproximate[i - 1] + step * derivative(nodes[i - 1], valuesApproximate[i - 1], variant);
        valuesApproximate[i] = valuesApproximate[i - 1] + stepHalf *  \
            (derivative(nodes[i - 1], valuesApproximate[i - 1], variant) + \
             derivative(nodes[i], valuesApproximate[i], variant));
    }

    cout << "Вычисления, проведённые с помощью метода предварительного и корректирующего счёта:\n\n";
    printTable(nodes, valuesApproximate, valuesAccurate, n);


	return 0;
}