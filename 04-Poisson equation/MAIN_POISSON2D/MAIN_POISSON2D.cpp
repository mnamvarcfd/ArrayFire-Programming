//////////////////////////////////////////////////////////////
// DATE: 2020 - 09 - 22
// Code written by ************
// Empty project
//////////////////////////////////////////////////////////////
#include "Grid.h"
#include "LinearSys.h"
#include "Field.h"
#include "AnalyticalSolution.h"



#include "arrayfire.h"
#include "IO/IO.h"

#undef min
#undef max

#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip> // std::setprecision

# define M_PI = atan(1.0)  /* pi */

typedef double T;

template <typename T>
af::array poisson2D_numerical_naive_c(int nx, int ny, int nt, T dx, T dy, T dt, T alpha)
{

	// Un reserve for the solution at time "n" and Unp1 reserve for time "n+1"
	T* Un = new T[nx * ny];
	T* Unp1 = new T[nx * ny];

	////Applying initial values 
	//for (int iNode = 0; iNode < nx * ny; ++iNode) {
	//	Unp1[iNode] = (T)1;
	//}

	////Applying boundary conditions 
	//int index;
	//for (int i = 0; i < nx; i++)
	//{
	//	index = getIndex(nx, ny, i, 0);
	//	Unp1[index] = (T)0;

	//	index = getIndex(nx, ny, i, ny - 1);
	//	Unp1[index] = (T)0;
	//}

	//for (int j = 0; j < ny; j++)
	//{
	//	index = getIndex(nx, ny, 0, j);
	//	Unp1[index] = (T)0;

	//	index = getIndex(nx, ny, nx - 1, j);
	//	Unp1[index] = (T)0;
	//}

	////Calculating a coeeficient to be used in the "for loop" 
	//T coeff = alpha * dt / (dx * dx);

	////Main loop which marching the time
	//for (int it = 0; it < nt; it++)
	//{

	//	//substitiud new value to the old values
	//	for (int i = 0; i < nx * ny; i++)
	//	{
	//		Un[i] = Unp1[i];
	//	}

	//	//marching in the x and y direction to solve equation for non-boundary nodes
	//	for (int i = 1; i < nx - 1; i++)
	//	{
	//		for (int j = 1; j < ny - 1; j++)
	//		{
	//			int ij = getIndex(nx, ny, i, j);
	//			int ip1j = getIndex(nx, ny, i + 1, j);
	//			int im1j = getIndex(nx, ny, i - 1, j);
	//			int ijp1 = getIndex(nx, ny, i, j + 1);
	//			int ijm1 = getIndex(nx, ny, i, j - 1);

	//			Unp1[ij] = Un[ij] + coeff * (Un[ip1j] + Un[im1j] + Un[ijp1] + Un[ijm1] - 4 * Un[ij]);
	//		}
	//	}

	//	//Applying boundary conditions 
	//	int index;
	//	for (int i = 0; i < nx; i++)
	//	{
	//		index = getIndex(nx, ny, i, 0);
	//		Unp1[index] = (T)0;

	//		index = getIndex(nx, ny, i, ny - 1);
	//		Unp1[index] = (T)0;
	//	}

	//	for (int j = 0; j < ny; j++)
	//	{
	//		index = getIndex(nx, ny, 0, j);
	//		Unp1[index] = (T)0;

	//		index = getIndex(nx, ny, nx - 1, j);
	//		Unp1[index] = (T)0;
	//	}

	//}

	//Copy solution to an ArrayFire array as the output value of the function
	af::array uSol = af::array(nx, ny, Unp1);

	//Uncomment bellow to write the solution into VTK format
	//writeScalarInVTK(nx, ny, Unp1, "numeric");

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//


	// The array "uSol" return by this function should be the numerical solution after "nt" time step.
	return uSol;
}

double function(double x, double y) {
	return -sin(x) - cos(y);
}


int main(int argc, char *argv[])
{
	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	//float v[] = { 5, 8, 3, 6 };
	//int r[] = { 0, 0, 2, 3, 4 };
	//int c[] = { 0, 1, 2, 1 };
	//const int M = 4, N = 4, nnz = 4;
	//af::array vals = af::array(af::dim4(nnz), v);
	//af::array row_ptr = af::array(af::dim4(M + 1), r);
	//af::array col_idx = af::array(af::dim4(nnz), c);
	//// Create sparse array (CSR) from af::arrays containing values,
	//// row pointers, and column indices.
	//af::array sparse = af::sparse(M, N, vals, row_ptr, col_idx, AF_STORAGE_CSR);
	//// sparse
	////     values:  [ 5.0, 8.0, 3.0, 6.0 ]
	////     row_ptr: [ 0, 0, 2, 3, 4 ]
	////     col_idx: [ 0, 1, 2, 1 ]

	//af_print(sparse);

	//af::array dens = af::dense(sparse);

	//af_print(dens);


	//////int Nx = 80;
	//////int Ny = Nx;
	//////double xMin = 0.0;
	//////double xMax = 1.0;
	//////double yMin = 0.0;
	//////double yMax = 1.0;


	//////Field u = Field(Nx, Ny, xMin, xMax, yMin, yMax);
	//////u.init();
	//////u.setBC();
	//////u.write("init");
	//////


	//////LinearSys linearSys = LinearSys(u);
	//////linearSys.creatCoeffMatrix();
	//////linearSys.writeCoeffMatrix("Matrix");

	//////linearSys.creatRigtHandSid(&function);
	////////printf("RHS[0] = %.4f ---->  %.4f\n", linearSys.RHS[0], -u.var[1] - u.var[4] - (sin(1.0) + cos(1.0)));
	////////printf("RHS[1] = %.4f ---->  %.4f\n", linearSys.RHS[1], -u.var[7] - u.var[2] - (sin(2.0) + cos(1.0)));
	////////printf("RHS[2] = %.4f ---->  %.4f\n", linearSys.RHS[2], -u.var[8] - u.var[13]- (sin(1.0) + cos(2.0)));
	////////printf("RHS[3] = %.4f ---->  %.4f\n", linearSys.RHS[3], -u.var[11]- u.var[14]- (sin(2.0) + cos(2.0)));

	//////linearSys.solve();


	//////u.write("numericalSolution");



	//////AnalyticalSolution analyticalSolution = AnalyticalSolution(Nx, Ny, xMin, xMax, yMin, yMax);
	//////analyticalSolution.write("analyticalSolution");

	//////for (int j = 0; j < Nx*Ny; j++) {
	//////	printf("error[%d] = %.8f  Numeric:Anaytic: %.8f   %.8f \n", j, u.var[j]- analyticalSolution.var[j], u.var[j] , analyticalSolution.var[j]);
	//////}
	//////double sum = 0.0;
	//////for (int j = 0; j < Nx * Ny; j++) {
	//////	sum += (u.var[j] - analyticalSolution.var[j])* (u.var[j] - analyticalSolution.var[j]);
	//////}
	//////printf("error is = %.8f  \n", sqrt(sum));
	//////



	////////double x1 = 342839 / 240000.;
	////////double x2 = 179723 / 120000.0;
	////////double x3 = 10867 / 24000.;
	////////double x4 = 125273 / 240000.;

	////////printf("u[0] = %.4f ---->  %.4f ---->  %.4f\n", u.var[5], x1, function(1.0, 1.0));
	////////printf("u[1] = %.4f ---->  %.4f ---->  %.4f\n", u.var[6], x2, function(2.0, 1.0));
	////////printf("u[2] = %.4f ---->  %.4f ---->  %.4f\n", u.var[9], x3, function(1.0, 2.0));
	////////printf("u[3] = %.4f ---->  %.4f ---->  %.4f\n", u.var[10], x4, function(2.0, 2.0));



	printf("============== end ==============\n");
	return 0;
}
