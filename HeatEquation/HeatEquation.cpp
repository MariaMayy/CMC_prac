#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>

using namespace std;

const double Eps = 0.01; // epsilon
const int N = 100; // num points on x
const int T = 10; // num points on t

const double Pi = 3.14159265;
const double xStep = 1.0 / N;
const double tStep = 1.0 / T;

const double tMax = 100;

double timeMax = 1.0;
double timeStep = timeMax / T;



void AnalyticalSolution() {
	double u = 0;
	double t = 0; // start time
	int tMax = 100; // max time
	double j = 0; // starting point in space
	double hStep = 1.0 / N;
	double tStep = 1.0/T;

	ofstream result("AnaliticalResult.txt"); // open file for recording
	// check file opening and print error
	if (!result.is_open()) {
		cout << "Oops! Couldn't open the file, bro... :c" << endl;
	}
	while (t < tMax) {
		j = 0;
		while (j < 1) {
			u = sin(Pi * j) / exp(t);
			result << setfill(' ');
			result << j << " " << u << endl;
			j += hStep;
		}
		result << endl << endl;
		t += tStep;
	}

}

// explicit scheme for linear case
void GetExplicitLinear() {
	int i, j;

	double coef = tStep/ (xStep * xStep);

	// allocate memory
	double** u = new double* [N+1];

	for ( i = 0; i < N + 1; i++) {
		u[i] = new double[T+1];
	}

	// initialization
	for ( i = 0; i < N + 1; i++) {
		for ( j = 0; j < T + 1; j++) {
			u[i][j] = 0.0;
		}
	}

	// boundary condition
	for ( j = 0; j < T + 1; j++) {
		u[0][j] = 0.0;
		u[N][j] = 0.0;
	}

	// fill initial state
	for ( i = 0; i < N + 1; i++) {
		u[i][0] = sin(Pi * i * xStep);
	}

	// calculate result
	
	for (j = 0; j < T; j++) {
		for ( i = 1; i < N; i++) {
			u[i][j + 1] = u[i][j] + coef * (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]);
		}
	}

	ofstream resultChisl("ExplicitLinear.txt"); // open file for recording
	// check file opening and print error
	if (!resultChisl.is_open()) {
		cout << "Oops! Couldn't open the file, bro... :c" << endl;
	}

	for (int j = 0; j < T + 1; j++) {
		for (int i = 0; i < N + 1; i++) {
			resultChisl << i * xStep << " " << u[i][j] << endl;
		}
		resultChisl << endl << endl;
	}

	resultChisl.close();

}


// progonka
void  RMethod(double* A, double* B, double* C, double* f, double* res, int N) {
	vector <double> VectorY;
	double* delta;
	double* alfa;
	double* beta;
	double* y;
	int i;

	y = new double[N + 2];
	delta = new double[N + 2];
	alfa = new double[N + 2];
	beta = new double[N + 2];

	delta[0] = C[0];
	beta[1] = f[0] / C[0];

	alfa[1] = B[0] / C[0];

	for (i = 1; i < N; i++) {
		delta[i] = C[i] - A[i] * alfa[i];
		alfa[i + 1] = B[i] / delta[i];

		beta[i + 1] = (f[i] + A[i] * beta[i]) / delta[i];
	}
	beta[N + 1] = (f[N] + A[N] * beta[N]) / (C[N] - A[N] * alfa[N]);

	y[N] = beta[N + 1];
	for (i = N - 1; i >= 0; i--) {
		y[i] = alfa[i + 1] * y[i + 1] + beta[i + 1];
	}

	for (i = 0; i <= N; i++) {
		res[i] = y[i];
	}

	/*
	for (i = 0; i <= N; i++) {
		VectorY.push_back(y[i]);
	}

	return VectorY;
	*/

	delete[]alfa;
	delete[]beta;
	delete[]delta;
}

// implicit scheme for linear case
void GetImplicitLinear() {

	double coef = tStep/ (xStep * xStep);
	double* a = new double[N + 1]; // under diagonal coefficients in system
	double* b = new double[N + 1]; // diagonal coefficients in system
	double* c = new double[N + 1]; // up diagonal coefficients in system

	double* arrAlfa = new double[N + 1];
	double* arrBeta = new double[N + 1];
	double delta;
	double* f = new double[N+1];

	double* res = new double[N + 1];

	for (int i = 0; i < N + 1; i++) {
		f[i] = sin(Pi * i * xStep);
		res[i] = 0.0;
	}
	// allocate memory
	double** u = new double* [N + 1];
	for (int i = 0; i < N + 1; i++) {
		u[i] = new double[T + 1];
	}

	// initialization
	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < T + 1; j++) {
			u[i][j] = 0.0;
		}
	}

	for (int i = 0; i < N + 1; i++) {
		if (i > 0) a[i] = coef; // under diagonal
		if (i < N) b[i] = coef; // up diagonal
		c[i] = 1 + 2 * coef; // diagonal
	}

	for (int j = 0; j < T; j++) {
		u[0][j + 1] = 0.0;
		u[N][j + 1] = 0.0;
	}

	ofstream resultImplicitChisl("ImplicitLinear.txt"); // open file for recording
// check file opening and print error
	if (!resultImplicitChisl.is_open()) {
		cout << "Oops! Couldn't open the file, bro... :c" << endl;
	}

	for (int i = 0; i < N + 1; i++) {
		u[i][0] = sin(Pi * i * xStep);
		resultImplicitChisl << i * xStep << " " << u[i][0] << endl;
	}

	for (int j = 0; j < T + 1; j++) {
			RMethod(a, b, c, f, res, N);
			for (int i = 0; i < N + 1; i++) {
				resultImplicitChisl << i * xStep << " " << res[i] << endl;
				f[i] = res[i];
			}
			resultImplicitChisl << endl << endl;
	}


	resultImplicitChisl.close();

}

// calculate K coefficient
void CalculateK(double* K, double* res, double k1, double k2, int N) {
	// K = k1*u^k2
	for (int i = 0; i < N + 1; i++) {
		K[i] = k1 * powl(res[i], k2);
	}
}

// calculate F coefficient
void CalculateF(double* F, double* res, double f1, double f2, int N) {
	// F = f1*u^f2
	for (int i = 0; i < N + 1; i++) {
		F[i] = f1 * powl(res[i],f2);
	}
}

// copy res to res_new
void copy_res(double* res, double* res_new, int N) {
	for (int i = 0; i < N + 1; i++) {
		res_new[i] = res[i];
	}
}


// calculate start distrubution u(x, 0)
void CalculateStartU(double*StartU, int N) {
	double x;
	
	for (int i = 0; i < N + 1; i++) {
		x = i * xStep;
		
		//StartU[i] = 5*sin(Pi *x); // case 01

		//if ((x < 0.2) || (x > 0.8)) StartU[i] = 0.0; // case 02
		//else StartU[i] = -60*(x-0.2)*(x-0.8);

		//StartU[i] = abs(5.65 * sin(2 * Pi * x)); // case 03

		StartU[i] = 7*x*(1-cos(2 * Pi * x)* cos(2 * Pi * x)); // case 04
	}
}

// get width of halfmaximum 
double GetWidth(double* res, int N) {
	double Max, HalfMax;
	double width, xLeft, xRight;
	int i;

	Max = res[0];
	for (i = 1; i < N + 1; i++) {
		if (res[i] > Max) Max = res[i];
	}
	HalfMax = Max / 2.0;

	i = 0;
	while (res[i] < HalfMax) {
		i++;
	}
	xLeft = (i - 0.5) * xStep;

	i = N + 1;
	while (res[i] < HalfMax) {
		i--;
	}
	xRight = (i + 0.5) * xStep;

	width = xRight - xLeft;

	return width;
}


double GetError(double* res_1, double* res_2, int N) {
	double Error = 0.0;
	for (int i = 0; i < N + 1; i++) {
		Error += (res_1[i] - res_2[i])* (res_1[i] - res_2[i]);
	}
	Error = sqrt(Error);
	return Error;
}


// implicit scheme for NONlinear case
void GetImplicitNonLinear() {

	double coef;
	double* a = new double[N + 1]; // under diagonal coefficients in system
	double* b = new double[N + 1]; // diagonal coefficients in system
	double* c = new double[N + 1]; // up diagonal coefficients in system

	double* arrAlfa = new double[N + 1];
	double* arrBeta = new double[N + 1];
	double delta;
	double* f = new double[N + 1];

	double* res = new double[N + 1];
	double* res_t1 = new double[N + 1];
	double* res_t2 = new double[N + 1];
	double* K = new double[N + 1];
	double* F = new double[N + 1];

	double width;
	double time = 0.0;
	double timeStep2;

	// initialization
	
	CalculateStartU(res, N); // calculate start distribution u(x, 0)

	ofstream resultImplicitNonLinear("ImplicitNonLinear.txt"); // open file for recording
// check file opening and print error
	if (!resultImplicitNonLinear.is_open()) {
		cout << "Oops! Couldn't open the file, bro... :c" << endl;
	}
	
	ofstream resultWidth("Width.txt"); // open file for recording
// check file opening and print error
	if (!resultWidth.is_open()) {
		cout << "Oops! Couldn't open the file, bro... :c" << endl;
	}
		
	ofstream resultTime("Time.txt"); // open file for recording
// check file opening and print error
	if (!resultTime.is_open()) {
		cout << "Oops! Couldn't open the file, bro... :c" << endl;
	}

	width = GetWidth(res, N);
	for (int i = 0; i < N + 1; i++) {
		if (i == 0) resultImplicitNonLinear << i * xStep << " " << res[i] << " t=" << 0 << " " << width << endl;
		else resultImplicitNonLinear << i * xStep << " " << res[i] << endl;
	}
	resultImplicitNonLinear << endl << endl;
	resultWidth << 0 << " " << width << endl;
	
	resultTime << time << " " << timeStep << endl;

	while (time < timeMax) {

		copy_res(res, res_t1, N); 
		copy_res(res, res_t2, N);

		// calculate with timeStep
		CalculateK(K, res_t1, 0.001, 1.0, N); // calculate thermal conductivity K
		CalculateF(F, res_t1, 0.003, 4.0, N); // calculate source heat F

		coef = timeStep / (xStep * xStep);
		for (int m = 0; m < N + 1; m++) {
			if (m > 0) a[m] = coef * (K[m - 1] + K[m]) / 2; // under diagonal
			if (m < N) b[m] = coef * (K[m] + K[m + 1]) / 2;// up diagonal
			if ((m > 0) && (m < N)) c[m] = 1 + coef * (K[m - 1] + 2 * K[m] + K[m + 1]) / 2; // diagonal
			f[m] = res_t1[m] + timeStep * F[m];
		}

		RMethod(a, b, c, f, res_t1, N);
		width = GetWidth(res_t1, N);

		// calculate with timeStep/2

		CalculateK(K, res_t2, 0.001, 1.0, N); // calculate thermal conductivity K
		CalculateF(F, res_t2, 0.003, 4.0, N); // calculate source heat F

		timeStep2 = timeStep / 2;

		for (int j = 0; j < 2; j++) {
			coef = timeStep2 / (xStep * xStep);
			for (int m = 0; m < N + 1; m++) {
				if (m > 0) a[m] = coef * (K[m - 1] + K[m]) / 2; // under diagonal
				if (m < N) b[m] = coef * (K[m] + K[m + 1]) / 2;// up diagonal
				if ((m > 0) && (m < N)) c[m] = 1 + coef * (K[m - 1] + 2 * K[m] + K[m + 1]) / 2; // diagonal
				f[m] = res_t2[m] + timeStep2 * F[m];
			}

			RMethod(a, b, c, f, res_t2, N);
			width = GetWidth(res_t2, N);
		}

		if (GetError(res_t1, res_t2, N) < Eps) {
			copy_res(res_t1, res, N);
			time += timeStep;
			
			// output
			for (int i = 0; i < N + 1; i++) {
				if (i == 0) resultImplicitNonLinear << i * xStep << " " << res[i] << " t=" <<  time << " " << width << endl;
				else resultImplicitNonLinear << i * xStep << " " << res[i] << endl;
			}
			resultImplicitNonLinear << endl << endl;
			resultWidth << time << " " << width << endl;
			resultTime << time << " " << timeStep << endl;

		}
		else {
			timeStep = timeStep2;		
		}
		
	}

	resultImplicitNonLinear.close();
	resultWidth.close();
	resultTime.close();
}


int main()
{

	// AnalyticalSolution();

	// GetExplicitLinear();

	// GetImplicitLinear();

	GetImplicitNonLinear();

	return 0;
}
