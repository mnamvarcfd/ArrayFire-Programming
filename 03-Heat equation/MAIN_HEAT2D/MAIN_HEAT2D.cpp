/*******************************************************
 * Copyright (c) 2014, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * https://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#include <arrayfire.h>
#include <stdio.h>
#include <cstdlib>
using namespace af;
// use static variables at file scope so timeit() wrapper functions
// can reference image/kernels
// 
// image to convolve
static array img;

// 5x5 derivative with separable kernels
static float h_dx[] = { 1.f / 12,   -8.f / 12, 0,       8.f / 12, -1.f / 12 };  // five point stencil
static float h_spread[] = { 1.f / 5, 1.f / 5,  1.f / 5, 1.f / 5,  1.f / 5 };
static array dx, spread, kernel;  // device kernels
static array full_out, dsep_out, hsep_out;  // save output for value checks

// wrapper functions for timeit() below
static void full() { 
    full_out = convolve2(img, kernel); 
}

static void dsep() { 
    dsep_out = convolve(dx, spread, img); 
}

static bool fail(array& left, array& right) {
    return (max<float>(abs(left - right)) > 1e-6);
}

int main(int argc, char* argv[]) {
    try {
        //int device = argc > 1 ? atoi(argv[1]) : 0;
        //af::setDevice(device);
        //af::info();
	    af::setBackend(AF_BACKEND_CPU);
	    af::setDevice(0);
	    af::info();
 
 
        // setup image and device copies of kernels
        img = randu(640, 480);
        dx = array(5, 1, h_dx);      // 5x1 kernel
        spread = array(1, 5, h_spread);  // 1x5 kernel
        kernel = matmul(dx, spread);     // 5x5 kernel

        printf("full 2D convolution:         %.5f seconds\n", timeit(full));
        printf("separable, device pointers:  %.5f seconds\n", timeit(dsep));
        // ensure values are all the same across versions
        if (fail(full_out, dsep_out)) { throw af::exception("full != dsep"); }
    }
    catch (af::exception& e) { fprintf(stderr, "%s\n", e.what()); }
    return 0;
}



//////////////////////////////////////////////////////////////////////
////////// DATE: 2020 - 07 - 21
////////// Code written by Sébastien Leclaire(sebastien.leclaire@polymtl.ca)
////////// This code solves the 2D heat equation on a simple
////////// square domain using different implementations.
//////////////////////////////////////////////////////////////////////
////////
////////
////////#include "arrayfire.h"
////////#include "IO/IO.h"
////////
////////#undef min
////////#undef max
////////
////////#include <iostream>
////////#include <iomanip> // std::setprecision
////////# define M_PI 3.14159265358979323846  /* pi */
////////
////////typedef double T;
////////
////////
////////int getIndex(int nx, int ny, int i, int j)
////////{
////////	int index = i * nx + j;
////////
////////	return index;
////////}
////////
////////template <typename T>
////////af::array heat2D_analytical_naive_c(int nx, int ny, T u0, T L, T H, T time, T alpha)
////////{
////////	T* uAnalytical = new T[nx * ny];
////////	for (int iNode = 0; iNode < nx * ny; ++iNode)
////////		uAnalytical[iNode] = (T)0;
////////
////////	for (T mm = 1; mm <= 100; mm++)
////////	{
////////		for (T nn = 1; nn <= 100; nn++)
////////		{
////////			T numerator = (T)4 * H * L * u0 * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * nn / (T)2) * std::sin(M_PI * nn / (T)2) / (mm * nn * M_PI * M_PI);
////////			T denominator = (H * L * (std::sin((T)2 * M_PI * mm) - (T)2 * M_PI * mm) * (std::sin((T)2 * M_PI * nn) - (T)2 * M_PI * nn)) / ((T)16 * mm * nn * M_PI * M_PI);
////////			T A_mmnn = numerator / denominator;
////////			T lambda_mmnn = (nn * M_PI / L) * (nn * M_PI / L) + (mm * M_PI / H) * (mm * M_PI / H);
////////
////////			T expPart = std::exp(-alpha * lambda_mmnn * time);
////////
////////			for (int iRow = 0; iRow < nx; iRow++)
////////				for (int iCol = 0; iCol < ny; iCol++)
////////				{
////////					int index = iRow * ny + iCol;
////////
////////					T xPos = T(iRow) / (T(nx) - (T)1);
////////					T yPos = T(iCol) / (T(ny) - (T)1);
////////
////////					T xPart = std::sin(nn * M_PI * xPos / L);
////////					T yPart = std::sin(mm * M_PI * yPos / H);
////////
////////					uAnalytical[index] = uAnalytical[index] + A_mmnn * xPart * yPart * expPart;
////////				}
////////		}
////////	}
////////
////////	af::array uAnalyticalArray = af::array(nx, ny, uAnalytical);
////////	delete[] uAnalytical;
////////
////////	return uAnalyticalArray;
////////}
////////
////////template <typename T>
////////af::array heat2D_numerical_naive_c(int nx, int ny, int nt, T dx, T dy, T dt, T alpha)
////////{
////////	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
////////	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//
////////
////////	//Initialazation 
////////	T* Un = new T[nx * ny];
////////	T* Unp1 = new T[nx * ny];
////////	for (int iNode = 0; iNode < nx * ny; ++iNode)
////////		Unp1[iNode] = (T)1;
////////
////////	int index;
////////	for (int i = 0; i < nx; i++)
////////	{
////////		index = getIndex(nx, ny, i, 0);
////////		Unp1[index] = (T)0;
////////
////////		index = getIndex(nx, ny, i, ny - 1);
////////		Unp1[index] = (T)0;
////////	}
////////
////////	for (int j = 0; j < ny; j++)
////////	{
////////		index = getIndex(nx, ny, 0, j);
////////		Unp1[index] = (T)0;
////////
////////		index = getIndex(nx, ny, nx - 1, j);
////////		Unp1[index] = (T)0;
////////	}
////////
////////	T coeff = alpha * dt / (dx * dx);
////////
////////	for (int it = 0; it < nt; it++)
////////	{
////////
////////		for (int i = 0; i < nx * ny; i++)
////////		{
////////			Un[i] = Unp1[i];
////////		}
////////
////////
////////		for (int i = 1; i < nx - 1; i++)
////////		{
////////			for (int j = 1; j < ny - 1; j++)
////////			{
////////				int ij   = getIndex(nx, ny, i    , j    );
////////				int ip1j = getIndex(nx, ny, i + 1, j    );
////////				int im1j = getIndex(nx, ny, i - 1, j    );
////////				int ijp1 = getIndex(nx, ny, i    , j + 1);
////////				int ijm1 = getIndex(nx, ny, i    , j - 1);
////////
////////				Unp1[ij] = Un[ij] + coeff * (Un[ip1j] + Un[im1j] + Un[ijp1] + Un[ijm1] - 4 * Un[ij]);
////////
////////				//std::cout << ij << "   " << Unp1[ij] << std::endl;
////////				//std::cout << ij << "   " << ip1j << "   " << Un[ip1j] << std::endl;
////////				//std::cout << ij << "   " << im1j << "   " << Un[im1j] << std::endl;
////////				//std::cout << ij << "   " << ijp1 << "   " << Un[ijp1] << std::endl;
////////				//std::cout << ij << "   " << ijm1 << "   " << Un[ijm1] << std::endl;
////////				//std::cin >> ij;
////////			}
////////		}
////////
////////
////////		int index;
////////		for (int i = 0; i < nx; i++)
////////		{
////////			index = getIndex(nx, ny, i, 0);
////////			Unp1[index] = (T)0;
////////
////////			index = getIndex(nx, ny, i, ny - 1);
////////			Unp1[index] = (T)0;
////////		}
////////
////////		for (int j = 0; j < ny; j++)
////////		{
////////			index = getIndex(nx, ny, 0, j);
////////			Unp1[index] = (T)0;
////////
////////			index = getIndex(nx, ny, nx - 1, j);
////////			Unp1[index] = (T)0;
////////		}
////////
////////
////////
////////		//T sum = 0.0;
////////		//for (int i = 0; i < nx * ny; i++)
////////		//{
////////		//	sum += Un[i];
////////			//std::cout << i << "   " << Un[i] << std::endl;
////////		//}
////////		//std::cout << it << "   " << sum << std::endl;
////////		//int a;
////////		//std::cin >> a;
////////
////////	}
////////
////////	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
////////	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
////////
////////	// Dummy initialization for this homework statement to compile
////////	//////af::array uSol = af::constant(0, nx, ny, type::TYPE_AF<T>());
////////	// Instead cast your C "array" pointer to an ArrayFire "array" 
////////	// and return the ArrayFire "array" for post-processing as was done simiarly in 
////////	// the function heat2D_analytical_naive_c
////////	af::array uSol = af::array(nx, ny, Unp1/*uSol_c_pointer*/);
////////
////////
////////	// The array "uSol" return by this function should be the numerical solution after "nt" time step.
////////	return uSol;
////////}
////////
////////template <typename T>
////////af::array heat2D_numerical_arrayfire(int nx, int ny, int nt, T dx, T dy, T dt, T alpha)
////////{
////////	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
////////	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//
////////	
////////	//Initialazation 
////////	af::array Un = af::constant(1.0, nx, ny, type::TYPE_AF<T>());
////////	af::array Unp1 = af::constant(1.0, nx, ny, type::TYPE_AF<T>());
////////
////////	Unp1(af::span, 0) = (T)0;
////////	Unp1(af::span, ny - 1) = (T)0;
////////
////////	Unp1(0, af::span) = (T)0;
////////	Unp1(nx - 1, af::span) = (T)0;
////////
////////
////////	T coeff = alpha * dt / (dx * dx);
////////
////////	for (int it = 0; it < 1/*nt*/; it++)
////////	{
////////		Un = Unp1;
////////
////////		for (int i = 1; i < nx - 1; i++)
////////		{
////////			Unp1(i, af::span) = coeff * (Un(i+1, af::span) + Un(i-1, af::span) - 4 * Un(i, af::span));
////////		}
////////
////////		for (int j = 1; j < ny - 1; j++)
////////		{
////////			Unp1(af::span, j) += coeff * (Un(af::span, j + 1) + Un(af::span, j - 1));
////////		}
////////
////////		Unp1 += Un;
////////
////////		Unp1(af::span, 0) = (T)0;
////////		Unp1(af::span, ny-1) = (T)0;
////////
////////		Unp1(0, af::span) = (T)0;
////////		Unp1(nx - 1, af::span) = (T)0;
////////
////////	}
////////
////////	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
////////	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
////////
////////	// The array "uSol" return by this function should be the numerical solution after "nt" time step.
////////	return Unp1;
////////}
////////
////////template <typename T>
////////af::array heat2D_analytical_arrayfire(int nx, int ny, T u0, T L, T H, T time, T alpha)
////////{
////////	// Dummy initialization for this homework statement to compile
////////	//nx = 11;
////////	//ny = 11;
////////
////////	//af::timer timer_analytical_naive_c = af::timer::start();
////////	//af::sync();
////////	//printf("--------------------------------: %g\n", af::timer::stop(timer_analytical_naive_c));
////////
////////
////////
////////
////////	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
////////	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//
////////
////////	af::array uAnalytical = af::constant((T)0, nx, ny, type::TYPE_AF<T>());
////////
////////	////////T* xPos_h = new T[nx * ny];
////////	////////T* yPos_h = new T[nx * ny];
////////
////////	////////for (int i = 0; i <nx; i++)
////////	////////{
////////	////////	for (int j = 0; j < ny; j++)
////////	////////	{
////////	////////		int idx = getIndex(nx, ny, i, j);
////////	////////		xPos_h[idx] = T(i);
////////	////////		yPos_h[idx] = T(j);
////////	////////	}
////////	////////}
////////	////////af::array xPos(nx, ny, xPos_h);
////////	////////af::array yPos(nx, ny, yPos_h);
////////
////////	////////xPos /= (T(nx) - (T)1);
////////	////////yPos /= (T(ny) - (T)1);
////////
////////	////////delete[] xPos_h;
////////	////////delete[] yPos_h;
////////
////////	////////xPos.eval();
////////	////////yPos.eval();
////////
////////
////////
////////	////////af::array xPart = af::constant((T)0, nx, ny, type::TYPE_AF<T>());
////////	////////af::array yPart = af::constant((T)0, nx, ny, type::TYPE_AF<T>());
////////
////////
////////	////////for (T mm = 1; mm <= 100; mm++)
////////	////////{
////////	////////	for (T nn = 1; nn <= 100; nn++)
////////	////////	{
////////	////////		T numerator = (T)4 * H * L * u0 * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * mm / (T)2) * std::sin(M_PI * nn / (T)2) * std::sin(M_PI * nn / (T)2) / (mm * nn * M_PI * M_PI);
////////	////////		T denominator = (H * L * (std::sin((T)2 * M_PI * mm) - (T)2 * M_PI * mm) * (std::sin((T)2 * M_PI * nn) - (T)2 * M_PI * nn)) / ((T)16 * mm * nn * M_PI * M_PI);
////////	////////		T A_mmnn = numerator / denominator;
////////
////////	////////		T lambda_mmnn = (nn * M_PI / L) * (nn * M_PI / L) + (mm * M_PI / H) * (mm * M_PI / H);
////////
////////	////////		T expPart = std::exp(-alpha * lambda_mmnn * time);
////////
////////	////////		xPart = af::sin(nn * M_PI * xPos / L);
////////	////////		yPart = af::sin(mm * M_PI * yPos / H);
////////
////////	////////		xPart.eval();
////////	////////		yPart.eval();
////////
////////	////////		uAnalytical = uAnalytical + A_mmnn * xPart * yPart * expPart;
////////
////////	////////	}
////////	////////}
////////	////////uAnalytical.eval();
////////	//////////^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
////////	//////////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
////////	//////////af_print(xPart);
////////	//////////////af_print(yPart);
////////	//////////af_print(uAnalytical);
////////	//////////////std::cout << mm << std::endl;
////////	//////////int a;
////////	//////////std::cin >> a;
////////
////////	////////// The array "uSol" return by this function should be the analytical solution at time "time".
////////	////////// The variable "time" in this function also correspond to the same time as the 
////////	////////// numerical solution after "nt" time step.
////////	return uAnalytical;
////////}
////////
////////template <typename T>
////////af::array heat2D_numerical_arrayfire_convolution(int nx, int ny, int nt, T dx, T dy, T dt, T alpha)
////////{
////////	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
////////	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//
////////
////////	// Dummy initialization for this homework statement to compile
////////	af::array uSol = af::constant(0, nx, ny, type::TYPE_AF<T>());
////////
////////	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
////////	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
////////
////////	// The array "uSol" return by this function should be the numerical solution after "nt" time step.
////////	return uSol;
////////}
////////
////////int main(int argc, char* argv[])
////////{
////////	af::setBackend(AF_BACKEND_CPU);
////////	af::setDevice(0);
////////	af::info();
////////
////////	// Physical domain size
////////	T L = 1;
////////	T H = L;
////////	T tFinal = (T)0.01;
////////
////////	// Number of grid point in x-direction, y-direction
////////	int nx = 101;
////////	int ny = 101;
////////
////////	T dx = L / (T(nx) - (T)1); //  Spatial step in x
////////	T dy = dx; // Spatial step in y
////////	T alpha = (T)1; // Thermal diffusivity
////////
////////	// Stability condition on the time step
////////	T r = (T)1 / (T)4;
////////	T dt = r * dx * dx / alpha;
////////
////////	// Number of time steps
////////	int nt = std::ceil(tFinal / dt);
////////
////////	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////	// Note that the content below does not require any modification.
////////	// Only the definition and content of the functions above the main program may require some modifications.
////////	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////
////////	// Computation with analytical_naive_c implementation
////////	af::timer timer_analytical_naive_c = af::timer::start();
////////	af::array uSol_analytical_naive_c = heat2D_analytical_naive_c<T>(nx, ny, (T)1, L, H, tFinal, alpha);
////////	uSol_analytical_naive_c.eval();
////////	af::sync();
////////	printf("Elapsed seconds for computing uSol_analytical_naive_c: %g\n", af::timer::stop(timer_analytical_naive_c));
////////	IO::arrayToVTK<T>("uSol_analytical_naive_c", "uSol_analytical_naive_c", nx, ny, 1, uSol_analytical_naive_c);
////////	IO::arrayToFile<T>("uSol_analytical_naive_c", uSol_analytical_naive_c);
////////
////////	// Computation with numerical_naive_c implementation
////////	af::timer timer_numerical_naive_c = af::timer::start();
////////	af::array uSol_numerical_naive_c = heat2D_numerical_naive_c<T>(nx, ny, nt, dx, dy, dt, alpha);
////////	uSol_numerical_naive_c.eval();
////////	af::sync();
////////	printf("Elapsed seconds for computing uSol_numerical_naive_c: %g\n", af::timer::stop(timer_numerical_naive_c));
////////	IO::arrayToVTK<T>("uSol_numerical_naive_c", "uSol_numerical_naive_c", nx, ny, 1, uSol_numerical_naive_c);
////////	IO::arrayToFile<T>("uSol_numerical_naive_c", uSol_numerical_naive_c);
////////
////////	// Computation with numerical_arrayfire implementation
////////	af::timer timer_numerical_arrayfire = af::timer::start();
////////	af::array uSol_numerical_arrayfire = heat2D_numerical_arrayfire<T>(nx, ny, nt, dx, dy, dt, alpha);
////////	uSol_numerical_arrayfire.eval();
////////	af::sync();
////////	printf("Elapsed seconds for computing uSol_numerical_arrayfire: %g\n", af::timer::stop(timer_numerical_arrayfire));
////////	IO::arrayToVTK<T>("uSol_numerical_arrayfire", "uSol_numerical_arrayfire", nx, ny, 1, uSol_numerical_arrayfire);
////////	IO::arrayToFile<T>("uSol_numerical_arrayfire", uSol_numerical_arrayfire);
////////
////////	// Computation with numerical_arrayfire_convolution implementation
////////	af::timer timer_numerical_arrayfire_convolution = af::timer::start();
////////	af::array uSol_numerical_arrayfire_convolution = heat2D_numerical_arrayfire_convolution(nx, ny, nt, dx, dy, dt, alpha);
////////	uSol_numerical_arrayfire_convolution.eval();
////////	af::sync();
////////	printf("Elapsed seconds for computing uSol_numerical_arrayfire_convolution: %g\n", af::timer::stop(timer_numerical_arrayfire_convolution));
////////	IO::arrayToVTK<T>("uSol_numerical_arrayfire_convolution", "uSol_numerical_arrayfire_convolution", nx, ny, 1, uSol_numerical_arrayfire_convolution);
////////	IO::arrayToFile<T>("uSol_numerical_arrayfire_convolution", uSol_numerical_arrayfire_convolution);
////////
////////	// Computation with analytical_arrayfire implementation
////////	af::timer timer_analytical_arrayfire = af::timer::start();
////////	af::array uSol_analytical_arrayfire = heat2D_analytical_arrayfire<T>(nx, ny, (T)1, L, H, tFinal, alpha);
////////	uSol_analytical_arrayfire.eval();
////////	af::sync();
////////	printf("Elapsed seconds for computing uSol_analytical_arrayfire: %g\n", af::timer::stop(timer_analytical_arrayfire));
////////	IO::arrayToVTK<T>("uSol_analytical_arrayfire", "uSol_analytical_arrayfire", nx, ny, 1, uSol_analytical_arrayfire);
////////	IO::arrayToFile<T>("uSol_analytical_arrayfire", uSol_analytical_arrayfire);
////////
////////	// L2 error computation between the analytical solution and the numerical solution for all three implementations
////////	T L2error_naive_c = std::sqrt(af::sum<T>((uSol_numerical_naive_c - uSol_analytical_naive_c) * (uSol_numerical_naive_c - uSol_analytical_naive_c)) * dx * dy);
////////	std::cout << std::setw(35) << "L2error_naive_c: " << std::setw(20) << std::setprecision(16) << L2error_naive_c << std::endl;
////////	T L2error_arrayfire = std::sqrt(af::sum<T>((uSol_numerical_arrayfire - uSol_analytical_naive_c) * (uSol_numerical_arrayfire - uSol_analytical_naive_c)) * dx * dy);
////////	std::cout << std::setw(35) << "L2error_arrayfire: " << std::setw(20) << std::setprecision(16) << L2error_arrayfire << std::endl;
////////	T L2error_arrayfire_convolution = std::sqrt(af::sum<T>((uSol_numerical_arrayfire_convolution - uSol_analytical_naive_c) * (uSol_numerical_arrayfire_convolution - uSol_analytical_naive_c)) * dx * dy);
////////	std::cout << std::setw(35) << "L2error_arrayfire_convolution: " << std::setw(20) << std::setprecision(16) << L2error_arrayfire_convolution << std::endl;
////////
////////	// L1 norm difference between naive_c implementation and your implementation
////////	std::cout << std::setw(35) << "L1_diff_analytical: " << std::setprecision(16) << af::sum<T>(af::abs(uSol_analytical_naive_c - uSol_analytical_arrayfire)) << std::endl;
////////	std::cout << std::setw(35) << "L1_diff_arrayfire: " << std::setprecision(16) << af::sum<T>(af::abs(uSol_numerical_naive_c - uSol_numerical_arrayfire)) << std::endl;
////////	std::cout << std::setw(35) << "L1_diff_arrayfire_convolution: " << std::setprecision(16) << af::sum<T>(af::abs(uSol_numerical_naive_c - uSol_numerical_arrayfire_convolution)) << std::endl;
////////
////////	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
////////	//vvvvvvv Write below the code you want for debugging or for any other reason vvvvvvvv//
////////
////////	return 0;
////////}
