//////////////////////////////////////////////////////////////
// DATE: 2020 - 09 - 30
// Code written by XYZ (XYZ@polymtl.ca)
// Empty project to be used with the ArrayFire library
//////////////////////////////////////////////////////////////


#include "arrayfire.h"
#include "IO/IO.h"

#undef min
#undef max

typedef double T;

#include <iostream>
#include <iomanip> // std::setprecision
# define M_PI 3.14159265358979323846  /* pi */

int main(int argc, char *argv[])
{
	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	return 0;
}
