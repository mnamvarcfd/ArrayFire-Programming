//////////////////////////////////////////////////////////////
// DATE: 2020 - 07 - 21
// Code written by Sébastien Leclaire(sebastien.leclaire@polymtl.ca)
// This code solves the 2D heat equation on a simple
// square domain using different implementations.
//////////////////////////////////////////////////////////////


#include "arrayfire.h"
#include "IO/IO.h"

#undef min
#undef max

#include <iostream>
#include <iomanip> // std::setprecision
# define M_PI 3.14159265358979323846  /* pi */

typedef double T;


int getIndex(int nx, int ny, int i, int j)
{
	int index = i * nx + j;

	return index;
}


void writeScalarInVTK(int nx, int ny, T scalar[], std::string fileName)
{
	fileName += ".vtk";

	// Declare the output flux and open the file
	std::ofstream file(fileName.c_str(), std::ios::out | std::ios::trunc);
	std::cout << "Saving file: " << fileName << std::endl;

	if (file)  // File has been opened correctly
	{
		file << "# vtk DataFile Version 2.0" << std::endl;
		file << "Test" << std::endl;
		file << "ASCII" << std::endl;
		file << "DATASET STRUCTURED_POINTS" << std::endl;
		file << "DIMENSIONS" << " " << nx << " " << ny << " " << 1 << std::endl;
		file << "ORIGIN" << " " << 0.0 << " " << 0.0 << " " << 0. << std::endl;
		file << "SPACING" << " " << 1. << " " << 1. << " " << 1. << std::endl;
		file << "POINT_DATA" << " " << nx * ny << std::endl;
		file << "SCALARS" << " " << "ScalarValue" << " " << "double" << " " << 1 << std::endl;
		file << "LOOKUP_TABLE default" << std::endl;
		for (int i = 0; i < nx * ny; ++i)
		{
			file << scalar[i] << std::endl;

		}
		file.close();  // close the file
	}

}



template <typename T>
af::array heat2D_analytical_naive_c(int nx, int ny, T u0, T L, T H, T time, T alpha)
{
	T* uAnalytical = new T[nx * ny];
	for (int iNode = 0; iNode < nx * ny; ++iNode)
		uAnalytical[iNode] = (T)0;

	for (T mm = 1; mm <= 100; mm++)
	{
		for (T nn = 1; nn <= 100; nn++)
		{
			T numerator = (T)4 * H * L * u0 * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * nn / (T)2) * std::sin(M_PI * nn / (T)2) / (mm * nn * M_PI * M_PI);
			T denominator = (H * L * (std::sin((T)2 * M_PI * mm) - (T)2 * M_PI * mm) * (std::sin((T)2 * M_PI * nn) - (T)2 * M_PI * nn)) / ((T)16 * mm * nn * M_PI * M_PI);
			T A_mmnn = numerator / denominator;
			T lambda_mmnn = (nn * M_PI / L) * (nn * M_PI / L) + (mm * M_PI / H) * (mm * M_PI / H);

			T expPart = std::exp(-alpha * lambda_mmnn * time);

			for (int iRow = 0; iRow < nx; iRow++)
				for (int iCol = 0; iCol < ny; iCol++)
				{
					int index = iRow * ny + iCol;

					T xPos = T(iRow) / (T(nx) - (T)1);
					T yPos = T(iCol) / (T(ny) - (T)1);

					T xPart = std::sin(nn * M_PI * xPos / L);
					T yPart = std::sin(mm * M_PI * yPos / H);

					uAnalytical[index] = uAnalytical[index] + A_mmnn * xPart * yPart * expPart;
				}
		}
	}

	af::array uAnalyticalArray = af::array(nx, ny, uAnalytical);

	//Uncomment bellow to write the solution into VTK format
	//writeScalarInVTK(nx, ny, uAnalytical, "Analytic");

	delete[] uAnalytical;

	return uAnalyticalArray;
}

template <typename T>
af::array heat2D_numerical_naive_c(int nx, int ny, int nt, T dx, T dy, T dt, T alpha)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// Un reserve for the solution at time "n" and Unp1 reserve for time "n+1"
	T* Un = new T[nx * ny];
	T* Unp1 = new T[nx * ny];

	//Applying initial values 
	for (int iNode = 0; iNode < nx * ny; ++iNode) {
		Unp1[iNode] = (T)1;
	}

	//Applying boundary conditions 
	int index;
	for (int i = 0; i < nx; i++)
	{
		index = getIndex(nx, ny, i, 0);
		Unp1[index] = (T)0;

		index = getIndex(nx, ny, i, ny - 1);
		Unp1[index] = (T)0;
	}

	for (int j = 0; j < ny; j++)
	{
		index = getIndex(nx, ny, 0, j);
		Unp1[index] = (T)0;

		index = getIndex(nx, ny, nx - 1, j);
		Unp1[index] = (T)0;
	}

	//Calculating a coeeficient to be used in the "for loop" 
	T coeff = alpha * dt / (dx * dx);

	//Main loop which marching the time
	for (int it = 0; it < nt; it++)
	{

		//substitiud new value to the old values
		for (int i = 0; i < nx * ny; i++)
		{
			Un[i] = Unp1[i];
		}

		//marching in the x and y direction to solve equation for non-boundary nodes
		for (int i = 1; i < nx - 1; i++)
		{
			for (int j = 1; j < ny - 1; j++)
			{
				int ij   = getIndex(nx, ny, i    , j    );
				int ip1j = getIndex(nx, ny, i + 1, j    );
				int im1j = getIndex(nx, ny, i - 1, j    );
				int ijp1 = getIndex(nx, ny, i    , j + 1);
				int ijm1 = getIndex(nx, ny, i    , j - 1);

				Unp1[ij] = Un[ij] + coeff * (Un[ip1j] + Un[im1j] + Un[ijp1] + Un[ijm1] - 4 * Un[ij]);
			}
		}

		//Applying boundary conditions 
		int index;
		for (int i = 0; i < nx; i++)
		{
			index = getIndex(nx, ny, i, 0);
			Unp1[index] = (T)0;

			index = getIndex(nx, ny, i, ny - 1);
			Unp1[index] = (T)0;
		}

		for (int j = 0; j < ny; j++)
		{
			index = getIndex(nx, ny, 0, j);
			Unp1[index] = (T)0;

			index = getIndex(nx, ny, nx - 1, j);
			Unp1[index] = (T)0;
		}

		if (it == 1) {
			size_t alloc_bytes, alloc_buffers;
			size_t lock_bytes, lock_buffers;
			af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
			std::cout << "Used Memory in Mb: " << alloc_bytes / (T)1024 / (T)1024 << std::endl;
		}
	}

	//Copy solution to an ArrayFire array as the output value of the function
	af::array uSol = af::array(nx, ny, Unp1);

	//Uncomment bellow to write the solution into VTK format
	//writeScalarInVTK(nx, ny, Unp1, "numeric");

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//


	// The array "uSol" return by this function should be the numerical solution after "nt" time step.
	return uSol;
}

template <typename T>
af::array heat2D_numerical_arrayfire(int nx, int ny, int nt, T dx, T dy, T dt, T alpha)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// Un and Unp1 reserve for the solution at time "n" for time "n+1". 
	//At the same time initialization is done. 
	af::array Un = af::constant(1.0, nx, ny, type::TYPE_AF<T>());
	af::array Unp1 = af::constant(1.0, nx, ny, type::TYPE_AF<T>());

	//Applying boundary conditions 
	Unp1(af::span, 0) = (T)0;  //Bottom nodes 
	Unp1(af::span, ny - 1) = (T)0;  //Top nodes 
	Unp1(0, af::span) = (T)0;  //left nodes 
	Unp1(nx - 1, af::span) = (T)0;  //Right nodes 

	//Calculating a coeeficient to be used in the "for loop" 
	T coeff = alpha * dt / (dx * dx);

	//Main loop which marching the time
	for (int it = 0; it < nt; it++)
	{

		//substitiud new value to the old values
		Un = Unp1;

		////marching in the x and y direction to solve equation for non-boundary nodes
		Unp1(af::seq(1, nx - 2,1), af::span) = coeff * (Un(af::seq(2, nx-1, 1), af::span) + Un(af::seq(0, nx - 3, 1), af::span));

		Unp1(af::span, af::seq(1, ny - 2, 1)) += coeff * (Un(af::span, af::seq(2, nx - 1, 1)) + Un(af::span, af::seq(0, nx - 3, 1)));
	
		Unp1 += (1- 4*coeff)*Un;

		//Applying boundary conditions 
		Unp1(af::span, 0) = (T)0;  //Bottom nodes 
		Unp1(af::span, ny - 1) = (T)0;  //Top nodes 
		Unp1(0, af::span) = (T)0;  //left nodes 
		Unp1(nx - 1, af::span) = (T)0;  //Right nodes 

	}

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	// The array "uSol" return by this function should be the numerical solution after "nt" time step.
	return Unp1;
}

template <typename T>
af::array heat2D_analytical_arrayfire(int nx, int ny, T u0, T L, T H, T time, T alpha)
{

	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	af::array uAnalytical = af::constant((T)0, nx, ny, type::TYPE_AF<T>());

	T* xPos_h = new T[nx * ny];
	T* yPos_h = new T[nx * ny];

	for (int i = 0; i <nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			int idx = getIndex(nx, ny, i, j);
			xPos_h[idx] = T(i);
			yPos_h[idx] = T(j);
		}
	}
	af::array xPos(nx, ny, xPos_h);
	af::array yPos(nx, ny, yPos_h);

	xPos /= (T(nx) - (T)1);
	yPos /= (T(ny) - (T)1);

	delete[] xPos_h;
	delete[] yPos_h;

	xPos.eval();
	yPos.eval();



	af::array xPart = af::constant((T)0, nx, ny, type::TYPE_AF<T>());
	af::array yPart = af::constant((T)0, nx, ny, type::TYPE_AF<T>());


	for (T mm = 1; mm <= 100; mm++)
	{
		for (T nn = 1; nn <= 100; nn++)
		{
			T numerator = (T)4 * H * L * u0 * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * nn / (T)2) * std::sin(M_PI * nn / (T)2) / (mm * nn * M_PI * M_PI);
			T denominator = (H * L * (std::sin((T)2 * M_PI * mm) - (T)2 * M_PI * mm) * (std::sin((T)2 * M_PI * nn) - (T)2 * M_PI * nn)) / ((T)16 * mm * nn * M_PI * M_PI);
			T A_mmnn = numerator / denominator;

			T lambda_mmnn = (nn * M_PI / L) * (nn * M_PI / L) + (mm * M_PI / H) * (mm * M_PI / H);

			T expPart = std::exp(-alpha * lambda_mmnn * time);

			xPart = af::sin(nn * M_PI * xPos / L);
			yPart = af::sin(mm * M_PI * yPos / H);

			xPart.eval();
			yPart.eval();

			uAnalytical = uAnalytical + A_mmnn * xPart * yPart * expPart;
			uAnalytical.eval();

		}
	}
	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	// The array "uSol" return by this function should be the analytical solution at time "time".
	// The variable "time" in this function also correspond to the same time as the 
	// numerical solution after "nt" time step.
	return uAnalytical;
}

template <typename T>
af::array heat2D_numerical_arrayfire_convolution(int nx, int ny, int nt, T dx, T dy, T dt, T alpha)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	//Initialazation 
	af::array Un = af::constant(1.0, nx, ny, type::TYPE_AF<T>());
	af::array Unp1 = af::constant(1.0, nx, ny, type::TYPE_AF<T>());

	Unp1(af::span, 0) = (T)0;
	Unp1(af::span, ny - 1) = (T)0;

	Unp1(0, af::span) = (T)0;
	Unp1(nx - 1, af::span) = (T)0;

	T coeff = alpha * dt / (dx * dx);

	T h_filter[] = { 0.0, coeff, 0.0,   coeff, 1.0 - 4 * coeff, coeff,   0.0, coeff, 0.0 };
	af::array filter(3, 3, h_filter);

	for (int it = 0; it < nt; it++)
	{
		Un = Unp1;

		Unp1 = convolve2(Un, filter);

		Unp1(af::span, 0) = (T)0;
		Unp1(af::span, ny - 1) = (T)0;

		Unp1(0, af::span) = (T)0;
		Unp1(nx - 1, af::span) = (T)0;

	}

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	// The array "uSol" return by this function should be the numerical solution after "nt" time step.
	return Unp1;
}

int main(int argc, char* argv[])
{
	af::setBackend(AF_BACKEND_CPU);  //AF_BACKEND_OPENCL   
	af::setDevice(0);
	af::info();

	// Physical domain size
	T L = 1;
	T H = L;
	T tFinal = (T)0.01;

	// Number of grid point in x-direction, y-direction
	int nx = 101;
	int ny = nx;

	T dx = L / (T(nx) - (T)1); //  Spatial step in x
	T dy = dx; // Spatial step in y
	T alpha = (T)1; // Thermal diffusivity

	// Stability condition on the time step
	T r = (T)1 / (T)4;
	T dt = r * dx * dx / alpha;

	// Number of time steps
	int nt = std::ceil(tFinal / dt);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Note that the content below does not require any modification.
	// Only the definition and content of the functions above the main program may require some modifications.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Computation with analytical_naive_c implementation
	af::timer timer_analytical_naive_c = af::timer::start();
	af::array uSol_analytical_naive_c = heat2D_analytical_naive_c<T>(nx, ny, (T)1, L, H, tFinal, alpha);
	uSol_analytical_naive_c.eval();
	af::sync();
	printf("Elapsed seconds for computing uSol_analytical_naive_c: %g\n", af::timer::stop(timer_analytical_naive_c));
	IO::arrayToVTK<T>("uSol_analytical_naive_c", "uSol_analytical_naive_c", nx, ny, 1, uSol_analytical_naive_c);
	IO::arrayToFile<T>("uSol_analytical_naive_c", uSol_analytical_naive_c);

	// Computation with numerical_naive_c implementation
	af::timer timer_numerical_naive_c = af::timer::start();
	af::array uSol_numerical_naive_c = heat2D_numerical_naive_c<T>(nx, ny, nt, dx, dy, dt, alpha);
	uSol_numerical_naive_c.eval();
	af::sync();
	printf("Elapsed seconds for computing uSol_numerical_naive_c: %g\n", af::timer::stop(timer_numerical_naive_c));
	IO::arrayToVTK<T>("uSol_numerical_naive_c", "uSol_numerical_naive_c", nx, ny, 1, uSol_numerical_naive_c);
	IO::arrayToFile<T>("uSol_numerical_naive_c", uSol_numerical_naive_c);

	// Computation with numerical_arrayfire implementation
	af::timer timer_numerical_arrayfire = af::timer::start();
	af::array uSol_numerical_arrayfire = heat2D_numerical_arrayfire<T>(nx, ny, nt, dx, dy, dt, alpha);
	uSol_numerical_arrayfire.eval();
	af::sync();
	printf("Elapsed seconds for computing uSol_numerical_arrayfire: %g\n", af::timer::stop(timer_numerical_arrayfire));
	IO::arrayToVTK<T>("uSol_numerical_arrayfire", "uSol_numerical_arrayfire", nx, ny, 1, uSol_numerical_arrayfire);
	IO::arrayToFile<T>("uSol_numerical_arrayfire", uSol_numerical_arrayfire);

	// Computation with numerical_arrayfire_convolution implementation
	af::timer timer_numerical_arrayfire_convolution = af::timer::start();
	af::array uSol_numerical_arrayfire_convolution = heat2D_numerical_arrayfire_convolution(nx, ny, nt, dx, dy, dt, alpha);
	uSol_numerical_arrayfire_convolution.eval();
	af::sync();
	printf("Elapsed seconds for computing uSol_numerical_arrayfire_convolution: %g\n", af::timer::stop(timer_numerical_arrayfire_convolution));
	IO::arrayToVTK<T>("uSol_numerical_arrayfire_convolution", "uSol_numerical_arrayfire_convolution", nx, ny, 1, uSol_numerical_arrayfire_convolution);
	IO::arrayToFile<T>("uSol_numerical_arrayfire_convolution", uSol_numerical_arrayfire_convolution);

	// Computation with analytical_arrayfire implementation
	af::timer timer_analytical_arrayfire = af::timer::start();
	af::array uSol_analytical_arrayfire = heat2D_analytical_arrayfire<T>(nx, ny, (T)1, L, H, tFinal, alpha);
	uSol_analytical_arrayfire.eval();
	af::sync();
	printf("Elapsed seconds for computing uSol_analytical_arrayfire: %g\n", af::timer::stop(timer_analytical_arrayfire));
	IO::arrayToVTK<T>("uSol_analytical_arrayfire", "uSol_analytical_arrayfire", nx, ny, 1, uSol_analytical_arrayfire);
	IO::arrayToFile<T>("uSol_analytical_arrayfire", uSol_analytical_arrayfire);

	// L2 error computation between the analytical solution and the numerical solution for all three implementations
	T L2error_naive_c = std::sqrt(af::sum<T>((uSol_numerical_naive_c - uSol_analytical_naive_c) * (uSol_numerical_naive_c - uSol_analytical_naive_c)) * dx * dy);
	std::cout << std::setw(35) << "L2error_naive_c: " << std::setw(20) << std::setprecision(16) << L2error_naive_c << std::endl;
	T L2error_arrayfire = std::sqrt(af::sum<T>((uSol_numerical_arrayfire - uSol_analytical_naive_c) * (uSol_numerical_arrayfire - uSol_analytical_naive_c)) * dx * dy);
	std::cout << std::setw(35) << "L2error_arrayfire: " << std::setw(20) << std::setprecision(16) << L2error_arrayfire << std::endl;
	T L2error_arrayfire_convolution = std::sqrt(af::sum<T>((uSol_numerical_arrayfire_convolution - uSol_analytical_naive_c) * (uSol_numerical_arrayfire_convolution - uSol_analytical_naive_c)) * dx * dy);
	std::cout << std::setw(35) << "L2error_arrayfire_convolution: " << std::setw(20) << std::setprecision(16) << L2error_arrayfire_convolution << std::endl;

	// L1 norm difference between naive_c implementation and your implementation
	std::cout << std::setw(35) << "L1_diff_analytical: " << std::setprecision(16) << af::sum<T>(af::abs(uSol_analytical_naive_c - uSol_analytical_arrayfire)) << std::endl;
	std::cout << std::setw(35) << "L1_diff_arrayfire: " << std::setprecision(16) << af::sum<T>(af::abs(uSol_numerical_naive_c - uSol_numerical_arrayfire)) << std::endl;
	std::cout << std::setw(35) << "L1_diff_arrayfire_convolution: " << std::setprecision(16) << af::sum<T>(af::abs(uSol_numerical_naive_c - uSol_numerical_arrayfire_convolution)) << std::endl;

	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code you want for debugging or for any other reason vvvvvvvv//

	return 0;
}
