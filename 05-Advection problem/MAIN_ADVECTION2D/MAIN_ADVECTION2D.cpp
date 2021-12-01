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

//#include "src/LaxWendroffDimX.h"
//#include "src/LaxWendroffDimY.h"

double bump(double x_, double y_, double t) {

	double val = 0.0;

	double x = x_ - 0.5 * t;
	double y = y_ + 0.3 * t;

	if ((x * x + y * y - 0.25) < 0.) val = exp(1.0 - 0.25 / (0.25 - x * x - y * y));

	return val;
}

double lexico(double x, double y) {

	double val = y*6+x;

	return val;
}

double constant(double x, double y) {

	double val = 1.0;

	return val;
}

int main(int argc, char* argv[])
{
	//AF_BACKEND_CPU;  AF_BACKEND_OPENCL;  AF_BACKEND_CUDA
	af::setBackend(AF_BACKEND_OPENCL);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

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








	////////Set the parameters to define the computaional domain
	//////int nx = 10000;
	//////int ny = nx;
	//////double xMin = -0.5;
	//////double xMax = 0.5;
	//////double yMin = -0.5;
	//////double yMax = 0.5;


	////////Create an object of Field class to discret the domain and define the variable of the euation
	//////Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	//////
	//////
	////////Initializing the field's variable
	//////U.init(&bump, 0.0);

	////////Writing the domain and the initial value in PLT format (uncomment bellow to write the file)
	////////U.write("init");

	////////Creating an object of a class corresponding to solving Poison's equation (Just uncomment one of the bellow)
	////////LaxWendroff solver = LaxWendroff(U, 1.0, 0.0);

	////////LaxWendroffDimX solver = LaxWendroffDimX(U, 0.5, -0.3);

	//////DonerCellUpWind solver = DonerCellUpWind(U, 0.5, -0.3);

	////////LaxWendroffDimSplit solver = LaxWendroffDimSplit(U, 0.5, -0.3);

	////////////Solve the system of equation by invoking the function
	//////////af::timer t1 = af::timer::start();
	//////////solver.solveParallel(0.9, 2.0);
	////////////solver.solve(0.9, 2, 2.0);
	//////////printf(" %g \n", af::timer::stop(t1));


	//////af::timer t1 = af::timer::start();
	//////solver.solve(0.9, 10000, 2.0);
	//////double timeSer = af::timer::stop(t1);
	//////
	//////printf("Time serial %g \n", timeSer);


	//////U.init(&bump, 0.0);

	//////t1 = af::timer::start();
	//////solver.solveParallel(0.9, 2.0);
	//////double timePar = af::timer::stop(t1);

	//////printf("Time parallel %g \n", timePar);


	//////printf("speed up %g = \n", timeSer / timePar);



	////////////The analytical solution could be obtained by creating and object of Field calss in 
	////////////initilizing it by providing the analytical solution as the value of initilizing
	//////////Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	//////////Ua.init(&bump, 2.0);


	////////////If desired write the analytical and numerical solution in PLT format
	////////////U.write("numericalSolution");
	////////////Ua.write("analyticalSolution");


	////////////Calculating the L1 error
	//////////double error = 0.0;
	//////////for (int j = 0; j < U.get_nNode(); j++) {
	//////////	error += abs(U.var[j] - Ua.var[j]);
	//////////}
	//////////printf("L1 error is = %.8f  \n", sqrt(error));

	return 0;
}




//size_t alloc_bytes, alloc_buffers;
//size_t lock_bytes, lock_buffers;
//af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
//std::cout << "Used Memory in Mb: " << alloc_bytes / (double)1024 / (double)1024 << std::endl;
//



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