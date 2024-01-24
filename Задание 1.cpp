#include <iostream>
#include <iomanip>

using namespace std;

template <typename number>
number f(number x, double eps, int* n) {
    number result = 0;
    number xCur = x;
    *n = 1;
    
    while ((xCur <= -eps || xCur >= eps)) {
        result += xCur;
        xCur *= x * x / (*n + 1) / (*n + 2) * (-1);
        *n += 2;
    }

    return result;
}

int main()
{
    setlocale(LC_ALL, "RUS");


    const int NODES = 23;
    const double START = -5.;
    const double DELTA = 2.;
    const double EPS = .0001;

    double nodes[NODES];
    double partialSums[NODES];
    int iterations[NODES];

    for (int i = 0; i < NODES; ++i) { // начальные условия
        nodes[i] = START + i * DELTA;
    }

    for (int i = 0; i < NODES; ++i) { // вычисления
        partialSums[i] = f(nodes[i], EPS, iterations + i);
    }
    

    const int PRECISION = 3;
    const int WINDOW_WIDTH = 120;
    const int LABEL_WIDTH = 20;
    const int NUM_WIDTH = 10;

    int maxNumsInLine = (WINDOW_WIDTH - LABEL_WIDTH) / NUM_WIDTH;
    int linesNumber = (NODES + maxNumsInLine - 1) / maxNumsInLine;
    int numsNumber;

    cout << fixed << setprecision(PRECISION) << left;

    for (int i = 0; i < linesNumber; ++i) {
        numsNumber = (NODES - i * maxNumsInLine >= maxNumsInLine) ? maxNumsInLine : NODES % maxNumsInLine;

        cout << setw(LABEL_WIDTH) << "Узлы";
        for (int j = 0; j < numsNumber; ++j) {
            cout << setw(NUM_WIDTH) << nodes[j + maxNumsInLine * i];
        }

        cout << '\n' << setw(LABEL_WIDTH) << "Частич. суммы";
        for (int j = 0; j < numsNumber; ++j) {
            cout << setw(NUM_WIDTH) << partialSums[j + maxNumsInLine * i];
        }

        cout << '\n' << setw(LABEL_WIDTH) << "Кол-во слагаемых";
        for (int j = 0; j < numsNumber; ++j) {
            cout << setw(NUM_WIDTH) << iterations[j + maxNumsInLine * i];
        }

        cout << "\n\n";
    }


    system("pause");
    return 0;
}