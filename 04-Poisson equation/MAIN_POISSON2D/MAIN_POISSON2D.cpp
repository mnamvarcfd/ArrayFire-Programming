//////////////////////////////////////////////////////////////
// Code written by Morteza Namvar
//////////////////////////////////////////////////////////////

#include "arrayfire.h"

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
	//AF_BACKEND_CPU;  AF_BACKEND_OPENCL;  AF_BACKEND_CUDA
	af::setBackend(AF_BACKEND_CUDA);
	af::setDevice(0);
	af::info(); std::cout << std::endl;


	//Set the parameters to define the computaional domain
	int nx = 100;
	int ny = 100;
	double xMin = 0.0;
	double xMax = 1.0;
	double yMin = 0.0;
	double yMax = 1.0;


	//Create an object of Field class to discret the domain and define the variable of the euation
	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);

	//Initializing the field's variable
	U.init();

	//Applying BC
	U.setBC();

	//Writing the domain and the initial value in PLT format (uncomment bellow to write the file)
	//U.write("init");
	

	//Creating an object of a class corresponding to solving Poison's equation (Just uncomment one of the bellow)
	JacobiSolver solver = JacobiSolver(U);
	//DenseSolver solver = DenseSolver(U);

	//Solve the system of equation by invoking the function
	af::timer t1 = af::timer::start();
	solver.solve();
	printf(" %g \n",  af::timer::stop(t1));


	//The analytical solution could be obtained by creating and object of Field calss in 
	//initilizing it by providing the analytical solution as the value of initilizing
	Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ua.init(&function);


	//If desired write the analytical and numerical solution in PLT format
	//U.write("numericalSolution");
	//Ua.write("analyticalSolution");


	//Calculating the L2 error
	double error = 0.0;
	for (int j = 0; j < U.get_nNode(); j++) {
		error += (U.var[j] - Ua.var[j])* (U.var[j] - Ua.var[j]);
	}
	printf("L2 error is = %.8f  \n", sqrt(error));
	
	return 0;
}

//size_t alloc_bytes, alloc_buffers;
//size_t lock_bytes, lock_buffers;
//af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
//std::cout << "Used Memory in Mb: " << alloc_bytes / (double)1024 / (double)1024 << std::endl;
