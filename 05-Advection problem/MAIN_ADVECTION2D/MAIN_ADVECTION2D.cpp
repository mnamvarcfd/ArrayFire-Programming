//////////////////////////////////////////////////////////////
// Code written by Morteza Namvar
//////////////////////////////////////////////////////////////
#pragma once
#include "arrayfire.h"

#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip> // std::setprecision

#include "src/Field.h"
#include "src/LaxWendroff.h"
#include "src/DonerCellUpWind.h"
#include "src/DonerCellUpWind.h"
#include "src/LaxWendroffDimSplit.h"

double bump(double x, double y) {

	double val = 0.0;

	if ((x * x + y * y) < 0.25) val = exp(1.0 - 0.25 / (0.25 - x * x - y * y));

	return val;
}

static double SquarWave(double x, double y) {

	double val = 0.0;

	if (abs(x) <= 0.125 && abs(y) <= 0.125) val = 1.0;

	return val;
}


int main(int argc, char* argv[])
{
	//AF_BACKEND_CPU;  AF_BACKEND_OPENCL;  AF_BACKEND_CUDA
	af::setBackend(AF_BACKEND_OPENCL);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	//Set the parameters to define the computaional domain
	int nx = 21;
	int ny = nx;
	double dx = 0.5 * (1.0 / (nx));
	double xMin = -0.5 + dx;
	double xMax = 0.5 - dx;
	double yMin = -0.5 + dx;
	double yMax = 0.5 - dx;

	//Advection equation coefficient
	double a = 1.0;
	double b = 0.0;

	//Create an object of Field class to discret the domain and define the variable of the euation
	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	
	
	//Initializing the field's variable
	U.init(&bump);

	//Writing the domain and the initial value in PLT format (uncomment bellow to write the file)
	U.write("init");

	//Creating an object of a class corresponding to solving Poison's equation (Just uncomment one of the bellow)
	LaxWendroff solver = LaxWendroff(U, a, b);
	//LaxWendroffDimX solver = LaxWendroffDimX(U, a, b);
	//LaxWendroffDimSplit solver = LaxWendroffDimSplit(U, a, b);
	//DonerCellUpWind solver = DonerCellUpWind(U, a, b);

	//Solve advection equation in serial mode
	af::timer t1 = af::timer::start();
	solver.solve(0.9, 10e6, 2.0);
	double timeSer = af::timer::stop(t1);
	printf("Time serial %g \n", timeSer);


	//If desired write the numerical solution in PLT format
	U.write("numericalSolution");


	//The analytical solution could be obtained by creating and object of Field calss in 
	//initilizing it by providing the analytical solution as the value of initilizing
	Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ua.init(&bump, a, b, 2.0);

	//If desired write the analytical solution in PLT format
	Ua.write("analyticalSolution");


	//Calculating the L1 error
	double error = 0.0;
	for (int j = 0; j < U.get_nNode(); j++) {
		error += abs(U.var[j] - Ua.var[j]);
	}
	error = error / U.get_nNode();

	printf("L1 error is = %.8f  \n", error);



	//Solve advection equation in parallel mode
	U.init(&bump);
	t1 = af::timer::start();
	solver.solveParallel(0.9, 2.0);
	double timePar = af::timer::stop(t1);
	printf("Time parallel %g \n", timePar);


	return 0;
}













//////////////////////////////////////////////////////////////////
////// Code written by Morteza Namvar
//////////////////////////////////////////////////////////////////
////
////#include "arrayfire.h"
////
////#include <stdio.h>
////#include <string>
////#include <iostream>
////#include <iomanip> // std::setprecision
////
////#include "src/Field.h"
////#include "src/ReimanSolver.h"
////
////
////class ClassA
////{
////public:
////	double* ary;
////
////public:
////	ClassA() {
////
////		ary = new double[10];
////
////		for (int i = 0; i < 10; i++)
////			ary[i] = 1.0;
////	}
////};
////
////
////class ClassB
////{
////public:
////	double* ptr;
////
////public:
////	ClassB(ClassA classA) {
////
////		ptr = classA.ary;
////
////
////		af::array ptr_d(10, ptr);
////
////		double* ary2 = new double[10];
////		ary2 = ptr_d.host<double>();
////
////		for (int i = 0; i < 10; i++)
////			ptr[i] = ary2[i];
////
////
////		ptr[0] = 12;
////	}
////
////};
////
////int main(int argc, char* argv[])
////{
////	//AF_BACKEND_CPU;  AF_BACKEND_OPENCL;  AF_BACKEND_CUDA
////	af::setBackend(AF_BACKEND_CPU);
////	af::setDevice(0);
////	af::info(); std::cout << std::endl;
////
////	ClassA classA;
////
////	ClassB classB(classA);
////
////	std::cout << "ary[0] in classA: " << classA.ary[0] << std::endl;
////
////	return 0;
////}




	//int N = 10;

	//double* A_h = new double[N];
	//for (int i = 0; i < N; i++) A_h[i] = i;

	//int idx_h[] = {0,1,0,1,1,1,1,0,1,0};

	//int indx_h[] = { 0,2,7,9 };



	//af::array A_d(N, A_h);
	//af_print(A_d);

	//gfor(af::seq i, 0, N-1) {

	//	A_d(i) = 0.0;

	//}
	//af_print(A_d);

	////A_d = af::moddims(A_d, 2, 5);
	////af_print(A_d);


	////for (int i = 0; i < N; i++) {
	////	A_h[i] = 0.0;
	////}
	////A_h = A_d.host<double>();
	////for (int i = 0; i < N; i++) {
	////	printf("%f  \n", A_h[i]);
	////}


	////af::array idx_d(N, idx_h);
	////af_print(idx_d);

	////af::array indx_d(4, indx_h);
	////af_print(indx_d);

	////af::array indx = af::where(idx_d==0);
	////af_print(indx);

	////af::array fValue = af::lookup(A_d, indx)+10;
	////af_print(fValue);

	//////fValue = 10;
	////af_print(fValue);

	////A_d(indx) = fValue;
	////af_print(A_d);