#include <iostream>
#include <iomanip>

using namespace std;

double xSubtract(double x, double* subtrahends, int n, int k) {
	double result = 1;

	for (int i = 0; i < n; ++i) {
		if (i != k) {
			result *= x - subtrahends[i];
		}
	}

	return result;
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

int main() {
	setlocale(LC_ALL, "RUS");


	const int NODES = 4;

	double nodes[NODES] = { -2., -1., 0., 2. };
	double values[NODES] = { -8., -1., 0., 8. };


	double newNodes[2 * NODES - 1] = { nodes[0], (nodes[0] + nodes[1]) / 2., // возможно, расчёт промежуточных узлов лучше автоматизировать
								       nodes[1], (nodes[1] + nodes[2]) / 2.,
								       nodes[2], (nodes[2] + nodes[3]) / 2.,
								       nodes[3] };

	double newValues[2 * NODES - 1];

	for (int i = 0; i < 2 * NODES - 1; ++i) {
		newValues[i] = 0;

		for (int j = 0; j < NODES; ++j) {
			newValues[i] += values[j] * xSubtract(newNodes[i], nodes, NODES, j) / xSubtract(nodes[j], nodes, NODES, j);
		}
	}


	printTable(nodes, values, NODES);
	printTable(newNodes, newValues, 2 * NODES - 1);


	system("pause");
	return 0;
}