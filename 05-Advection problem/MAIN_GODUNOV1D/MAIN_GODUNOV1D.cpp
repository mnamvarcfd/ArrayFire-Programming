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
# define M_PI 3.14159265358979323846  /* pi */

typedef double T;

int main(int argc, char *argv[])
{
	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	return 0;
}
