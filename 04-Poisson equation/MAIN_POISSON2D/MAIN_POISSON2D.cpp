//////////////////////////////////////////////////////////////
// DATE: 2020 - 09 - 22
// Code written by ************
// Empty project
//////////////////////////////////////////////////////////////


#include "arrayfire.h"
#include "IO/IO.h"

#undef min
#undef max

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

int main(int argc, char *argv[])
{
	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	return 0;
}
