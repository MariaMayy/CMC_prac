#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
const double eps = 0.01; // epsilon
const int N = 3; // number of equations

// ----- Function of calculation right part -----
// Y     - input :  array of variables
// M    - input : diffraction parameter
// A     - input : parameter
// FRight - output : array of right part
double* F(double* Y, double M, double A) {
    double* FRight = new double[N];
    FRight[0] = (-M * Y[0] + Y[1] * Y[2]);
    FRight[1] = (-M * Y[1] + Y[0] * (Y[2] - A));
    FRight[2] = (1 - Y[0] * Y[1]);

    return FRight;
}

// ----- Function Method Runge-Kytta -----
// Y     - input :  array of variables
// M    - input : diffraction parameter
// A     - input : parameter
// h     - input : step
// NextY - output : array of variables on next iteration
double* MethodRungeKytta(double M, double A, double* Y, double h) {
    double* ArrHalf = new double[N]; // additional array for intermediate calculations 
    double* NextY = new double[N];
    double* k1, * k2, * k3, * k4; // arrays of coefficients

    k1 = F(Y, M, A);

    ArrHalf[0] = Y[0] + (h / 2) * k1[0];
    ArrHalf[1] = Y[1] + (h / 2) * k1[1];
    ArrHalf[2] = Y[2] + (h / 2) * k1[2];
    k2 = F(ArrHalf, M, A);

    ArrHalf[0] = Y[0] + (h / 2) * k2[0];
    ArrHalf[1] = Y[1] + (h / 2) * k2[1];
    ArrHalf[2] = Y[2] + (h / 2) * k2[2];
    k3 = F(ArrHalf, M, A);

    ArrHalf[0] = Y[0] + h * k3[0];
    ArrHalf[1] = Y[1] + h * k3[1];
    ArrHalf[2] = Y[2] + h * k3[2];
    k4 = F(ArrHalf, M, A);

    NextY[0] = Y[0] + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    NextY[1] = Y[1] + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    NextY[2] = Y[2] + (h / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);

    return NextY;
}



int main()
{
    double* Y = new double[N]; // array of variables

    double A; // parameter
    double M; // difraction parameter
    double hStep = 0.01; // step
    double tMax = 400;

    cout << "A: ";
    cin >> A;
    cout << "M: ";
    cin >> M;

    ofstream result("Result.txt"); // open file for recording
    // check file opening and print error
    if (!result.is_open())
    {
        cout << "Oops! Couldn't open the file, bro... :с" << endl;
    }

    // setting the initial conditions
    Y[0] = 1 / sqrt((A + sqrt(4 * (A * A + M * M))) / (2 * M)) + eps; 
    cout << Y[0] << endl;
    Y[1] = sqrt((A + sqrt(4 * (A * A + M * M))) / (2 * M)) + eps;
    cout << Y[1] << endl;
    Y[2] = -M * ((A + sqrt(4 * (A * A + M * M))) / (2 * M)) + eps;
    cout << Y[2];
    result << setfill(' ');
    result << "# A = " << A << " M = " << M << " Eps = " << eps << endl;
    result << setw(10) << left << "# Time" << setw(20) << left << "Y_1" << setw(20) << left << "Y_2" << setw(20) << left << "Y_3" << endl;

    double t = 0;
    while (t < tMax) {
        for (int j = 0; j < 1; j++) {
            Y = MethodRungeKytta(M, A, Y, hStep);
            t += hStep;
        }
        
        result << setfill(' ');
        result << setw(10) << left << t << setw(20) << left << Y[0] << setw(20) << left << Y[1] << setw(20) << left << Y[2] << endl;
    }

    return 0;
}
