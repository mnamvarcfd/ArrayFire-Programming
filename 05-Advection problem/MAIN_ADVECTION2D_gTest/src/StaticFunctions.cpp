#ifndef StaticFunctions
#define StaticFunctions

#include <cmath>


static double bump(double x, double y) {
	
	double val = 0.0;

	if ((x * x + y * y) < 0.25) val = exp(1.0 - 0.25 / (0.25 - x * x - y * y));

	return val;
}

//static double SquarWave(double x, double y) {
//
//	double val = 0.0;
//
//	if (abs(x) <= 0.125 && abs(y) <= 0.125) val = 1.0;
//
//	return val;
//}

#endif