//////////////////////////////////////////////////////////////
// Code written by Morteza Namvar
//////////////////////////////////////////////////////////////

#include "arrayfire.h"
#include "IO/IO.h"

#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip> // std::setprecision

#include "Field.h"
#include "DenseSolver.h"
#include "JacobiSolver.h"

double function(double x, double y) {
	return sin(x) + cos(y);
}
int main(int argc, char *argv[])
{
	//AF_BACKEND_CPU=1;  AF_BACKEND_OPENCL
	af::setBackend(AF_BACKEND_OPENCL);
	af::setDevice(0);
	af::info(); std::cout << std::endl;


	int nx = 100;
	int ny = 100;
	double xMin = 0.0;
	double xMax = 1.0;
	double yMin = 0.0;
	double yMax = 1.0;


	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();
	//U.write("init");
	

	//JacobiSolver solver = JacobiSolver(U);
	DenseSolver solver = DenseSolver(U);

	af::timer t1 = af::timer::start();
	solver.solve();
	//U.write("numericalSolution");
	printf(" %g \n",  af::timer::stop(t1));
	//26.7462

	//Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	//Ua.init(&function);
	//Ua.write("analyticalSolution");


	//double error = 0.0;
	//for (int j = 0; j < U.get_nNode(); j++) {
	//	error += (U.var[j] - Ua.var[j])* (U.var[j] - Ua.var[j]);
	//}
	//printf("error is = %.8f  \n", sqrt(error));
	
	return 0;
}

//size_t alloc_bytes, alloc_buffers;
//size_t lock_bytes, lock_buffers;
//af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
//std::cout << "Used Memory in Mb: " << alloc_bytes / (double)1024 / (double)1024 << std::endl;
