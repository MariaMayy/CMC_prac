#include <iostream>

using namespace std;

// ----- MethodGauss -----
double* MethodGauss(double** matrix, double*y, int n) {
    double* x, dMax;
    int index, k = 0;

    x = new double[n];
    while (k < n) {
        // search string with max element a[i][k]
        dMax = abs(matrix[k][k]); 
        index = k;
        for (int i = k; i < n; i++) {
            if (abs(matrix[i][k]) > dMax) {
                index = i;
                dMax = matrix[i][k];
            }
        }

        // swap 2 strings
        for (int j = 0; j < n; j++) {
            double tmp = matrix[k][j];
            matrix[k][j] = matrix[index][j];
            matrix[index][j] = tmp;
        }
        double tmpY = y[k];
        y[k] = y[index];
        y[index] = tmpY;

        // normalized equations
        for (int i = k; i < n; i++) {
            double normFactor = matrix[i][k];
            if (normFactor == 0) continue; 
            for (int j = 0; j < n; j++) {
                matrix[i][j] = matrix[i][j] / normFactor;
            }
            y[i] = y[i] / normFactor;

            if (i == k) continue; // don't subtract the equation from itself

            for (int j = 0; j < n; j++) {
                matrix[i][j] = matrix[i][j] - matrix[k][j];
            }
            y[i] = y[i] - y[k];
        }
        k++;
    }

    // reverse motion
    for (k = n - 1; k >= 0; k--) {
        x[k] = y[k]; // last equation
        for (int i = 0; i < n; i++) {
            y[i] = y[i] - matrix[i][k] * x[k]; // matrix[i][k] - normFactor on this step
        }
    }

    return x;
}

// ----- PrintSystem -----
void PrintSystem(double** matrix, double* y, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << "*x" << j;
            if (j < n - 1) cout << "+";
        }
        cout << " = " << y[i] << '\n';
    }
}


int main()
{
    int n; // count of equation
    double** matrixA, * vectorX, * vectorY;

    cout << "Write count of equation: ";
    cin >> n;
    matrixA = new double*[n];
    vectorY = new double[n];
    // enter elements of matrix A
    for (int i = 0; i < n; i++) {
        matrixA[i] = new double[n]; 
        for (int j = 0; j < n; j++) {
            cout << "a[" << i << "][" << j << "] = ";
            cin >> matrixA[i][j];
        }
    }
    // enter elements of right side (vector Y)
    for (int k = 0; k < n; k++) {
        cout << "y[" << k << "] = ";
        cin >> vectorY[k];
    }

    PrintSystem(matrixA, vectorY, n);
    vectorX = MethodGauss(matrixA, vectorY, n);

    // print result - vector X
    for (int k = 0; k < n; k++) {
        cout << "x[" << k << "] = " << vectorX[k] << endl;
    }

    return 0;
}
