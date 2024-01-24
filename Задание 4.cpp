#include <iostream>
#include <iomanip>

using namespace std;

double different(int firstIndex, int secondIndex, double** diffs, double* nodes) {
	return (diffs[firstIndex + 1][secondIndex - 1] - diffs[firstIndex][secondIndex - 1]) / (nodes[firstIndex + secondIndex] - nodes[firstIndex]);
}

void printTable(double* nodes, double* values, int n) {
	const int PRECISION = 3;
	const int WINDOW_WIDTH = 120;
	const int LABEL_WIDTH = 15;
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


	double** diffs = new double* [NODES];

	for (int i = 0; i < NODES; ++i) {
		diffs[i] = new double[NODES - i];
		diffs[i][0] = values[i];
	}

	for (int j = 1; j < NODES; ++j) {
		for (int i = 0; i < NODES - j; ++i) {
			diffs[i][j] = different(i, j, diffs, nodes);
		}
	}


	double newNodes[2 * NODES - 1] = { nodes[0], (nodes[0] + nodes[1]) / 2., // возможно, расчёт промежуточных узлов лучше автоматизировать
								       nodes[1], (nodes[1] + nodes[2]) / 2.,
								       nodes[2], (nodes[2] + nodes[3]) / 2.,
								       nodes[3] };
	
	double newValues[2 * NODES - 1];

	double xSubtract;
	for (int i = 0; i < 2 * NODES - 1; ++i) {
		newValues[i] = 0.;
		xSubtract = 1.;
		
		for (int j = 0; j < NODES; ++j) {
			newValues[i] += diffs[0][j] * xSubtract;
			xSubtract *= newNodes[i] - nodes[j];
		}
	}


	cout << "Разделённые разности:\n";
	cout << fixed << setprecision(3) << right;
	for (int i = 0; i < NODES; ++i) {
		for (int j = 0; j < NODES - i; ++j) {
			cout << setw(7) << diffs[i][j] << ' ';
		}

		cout << '\n';
	}
	cout << '\n';

	cout << "Таблица значений интерполируемой функции в исходных узлах:\n";
	printTable(nodes, values, NODES);
	cout << "Таблица значений интерполирующей функции в исходных узлах и узлах интерполяции:\n";
	printTable(newNodes, newValues, 2 * NODES - 1);


	system("pause");
	return 0;
}