////////////////////////////////////////////////////////////////////////////////////////////////
// DATE: 2020 - 09 - 01
// Code written by Sébastien Leclaire(sebastien.leclaire@polymtl.ca)
// Once completed, this code aims to provide an N-dimensional Monte Carlo integration algorithm.
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "arrayfire.h"
#include "IO/IO.h"

#undef min
#undef max

#include <iostream>
#include <iomanip> // std::setprecision
# define M_PI 3.14159265358979323846  /* pi */

typedef double T;

af::array myFunctionEx1(const af::array& xPos, const af::array& yPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	af::array fValue = xPos * xPos - yPos * yPos;

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

af::array myFunctionEx2(const af::array& xPos, const af::array& yPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	double radius = 0.5;

	af::array R = sqrt(xPos * xPos + yPos * yPos);
	//af_print(R);

	af::array indx = af::where(R <= radius);
	//af_print(indx);

	af::array fValue = af::lookup(R, indx);
	//af_print(fValue);

	fValue = 4.0 /* fValue * fValue*M_PI*/;
	


	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

T monteCarloIntegral2(
	af::array (*myFunction)(const af::array&, const af::array&),
	T xMin, T xMax, T yMin, T yMax, T absTol, long long int nBatch, long long int itMax)
{
	//This algorithm calculates an 2D integral using the Monte Carlo method.
	//
	//myFunction : a function pointer of a function taking 2 arguments
	//xMin, xMax, yMin, yMax : scalar defining the bounding box
	//of the Monte Carlo algorithm;
	//absTol: Absolute tolerance between the computed integralApproximation and the "last" computed integralApproximation;
	//nBatch: Number of function evaluations before testing the stopping criterion;
	//itMax: Maximum of number of iterations of the algorithm, which would lead to a maximum of nBatch*itMax function evaluations.

	// Reset ArrayFire random number generator
	af::setSeed(0);

	// "First" estimate of the integralApproximation
	T integralApproximation = std::numeric_limits<T>::quiet_NaN();

	T volume = (xMax - xMin) * (yMax - yMin);
	T fSum = 0;

	// Iterative process to compute the integral approximation

	long long int it = 1;
	while(true)
	{
		//std::cout << "iteration:  " << it<< std::endl;
		T integralApproximation_last = integralApproximation;

		af::array xPos = xMin + (xMax - xMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());
		af::array yPos = yMin + (yMax - yMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());

		af::array fValue = myFunction(xPos, yPos);

		fSum = fSum + af::sum<T>(fValue);

		long long int nTotal = nBatch * it;

		T fMean = fSum / T(nTotal);

		integralApproximation = volume * fMean;

		if( std::abs((integralApproximation - integralApproximation_last)) < absTol )
			break;
		
		if( it >= itMax)
			break;
	
		it += 1;
	}	

	return integralApproximation;
}

af::array myFunctionEx3(const af::array& xPos, const af::array& yPos, const af::array& zPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	af::array fValue = xPos * xPos - yPos * yPos + zPos * zPos * zPos;

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

af::array myFunctionEx4(const af::array& xPos, const af::array& yPos, const af::array& zPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	double radius = 0.5;

	af::array R = sqrt(xPos * xPos + yPos * yPos + zPos * zPos);
	//af_print(R);

	af::array indx = af::where(R <= radius);
	//af_print(indx);

	af::array fValue = af::lookup(R, indx);
	//af_print(fValue);

	fValue = 6.0 /* fValue * fValue*M_PI*/;


	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

af::array myFunctionEx5(const af::array& xPos, const af::array& yPos, const af::array& zPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	af::array dis = 1.0 - (1.0 - sqrt(xPos * xPos + yPos * yPos) )* (1.0 - sqrt(xPos * xPos + yPos * yPos)) - zPos * zPos;
	//af_print(R);

	af::array indx = af::where(dis > 0.0);
	//af_print(indx);

	af::array fValue = af::lookup(dis, indx);
	//af_print(fValue);

	fValue = 1.0 /* fValue * fValue*M_PI*/;


	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

T monteCarloIntegral3(
	af::array(*myFunction)(const af::array&, const af::array&, const af::array&),
	T xMin, T xMax, T yMin, T yMax, T zMin, T zMax, T absTol, long long int nBatch, long long int itMax)
{
	//This algorithm calculates an 3D integral using the Monte Carlo method.
	//
	//myFunction : a function pointer of a function taking 3 arguments
	//xMin, xMax, yMin, yMax, zMin, zMax : scalar defining the bounding box
	//of the Monte Carlo algorithm;
	//absTol: Absolute tolerance between the computed integralApproximation and the "last" computed integralApproximation;
	//nBatch: Number of function evaluations before testing the stopping criterion;
	//itMax: Maximum of number of iterations of the algorithm, which would lead to a maximum of nBatch*itMax function evaluations.

	// Reset ArrayFire random number generator
	af::setSeed(0);

	// "First" estimate of the integralApproximation
	T integralApproximation = std::numeric_limits<T>::quiet_NaN();

	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the integral approximation  vvvvvvvv//

	T volume = (xMax - xMin) * (yMax - yMin) * (zMax - zMin);
	//std::cout << "volume  " << volume << std::endl;
	T fSum = 0;

	// Iterative process to compute the integral approximation

	long long int it = 1;
	while (true)
	{
		//std::cout << "iteration:  " << it<< std::endl;
		T integralApproximation_last = integralApproximation;

		af::array xPos = xMin + (xMax - xMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());
		af::array yPos = yMin + (yMax - yMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());
		af::array zPos = zMin + (zMax - zMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());

		af::array fValue = myFunction(xPos, yPos, zPos);

		//std::cout << "====1111=====" << std::endl;
		//af_print(fValue);
		fSum = fSum + af::sum<T>(fValue);

		long long int nTotal = nBatch * it;

		T fMean = fSum / T(nTotal);
		//std::cout << "fMean  " << fMean << std::endl;

		integralApproximation = volume * fMean;

		if (std::abs((integralApproximation - integralApproximation_last)) < absTol)
			break;

		if (it >= itMax)
			break;

		it += 1;
	}




	//^^^^^^^ Write above the code that calculates the integral approximation  ^^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return integralApproximation;
}

af::array myFunctionEx6(const af::array& xPos, const af::array& yPos, const af::array& zPos, const af::array& tPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	af::array fValue = xPos * xPos - yPos * yPos + sqrt(zPos) * zPos - tPos * tPos;

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

af::array myFunctionEx7(const af::array& xPos, const af::array& yPos, const af::array& zPos, const af::array& tPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	af::array fValue = xPos * xPos * yPos * yPos * zPos * zPos * (cos(tPos/3.0) + 1.0);

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

T monteCarloIntegral4(
	af::array(*myFunction)(const af::array&, const af::array&, const af::array&, const af::array&),
	T xMin, T xMax, T yMin, T yMax, T zMin, T zMax, T tMin, T tMax, T absTol, long long int nBatch, long long int itMax)
{
	//This algorithm calculates an 4D integral using the Monte Carlo method.
	//
	//myFunction : a function pointer of a function taking 4 arguments
	//xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax : scalar defining the bounding box
	//of the Monte Carlo algorithm;
	//absTol: Absolute tolerance between the computed integralApproximation and the "last" computed integralApproximation;
	//nBatch: Number of function evaluations before testing the stopping criterion;
	//itMax: Maximum of number of iterations of the algorithm, which would lead to a maximum of nBatch*itMax function evaluations.

	// Reset ArrayFire random number generator
	af::setSeed(0);

	// "First" estimate of the integralApproximation
	T integralApproximation = std::numeric_limits<T>::quiet_NaN();

	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the integral approximation  vvvvvvvv//

	T volume = (xMax - xMin) * (yMax - yMin) * (zMax - zMin) * (tMax - tMin);
	T fSum = 0;


	// Iterative process to compute the integral approximation

	long long int it = 1;
	while (true)
	{
		//std::cout << "iteration:  " << it<< std::endl;
		T integralApproximation_last = integralApproximation;

		af::array xPos = xMin + (xMax - xMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());
		af::array yPos = yMin + (yMax - yMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());
		af::array zPos = zMin + (zMax - zMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());
		af::array tPos = tMin + (tMax - tMin) * af::randu(nBatch, 1, type::TYPE_AF<T>());

		af::array fValue = myFunction(xPos, yPos, zPos, tPos);

		//af_print(fValue);

		fSum = fSum + af::sum<T>(fValue);

		long long int nTotal = nBatch * it;

		T fMean = fSum / T(nTotal);

		integralApproximation = volume * fMean;

		if (std::abs((integralApproximation - integralApproximation_last)) < absTol)
			break;

		if (it >= itMax)
			break;

		it += 1;
	}



	//^^^^^^^ Write above the code that calculates the integral approximation  ^^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return integralApproximation;
}

af::array myFunctionExA(const af::array& xPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.
	//xPos


	af::array fValue = xPos(af::span, 0) * xPos(af::span, 0) - xPos(af::span, 1) * xPos(af::span, 1) + xPos(af::span, 2) * xPos(af::span, 2) * xPos(af::span, 2);
	//af_print(xPos(af::span, 1));
	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

af::array myFunctionExB(const af::array& xPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	af::array fValue = xPos(af::span, 0) * xPos(af::span, 0) - xPos(af::span, 1) * xPos(af::span, 1) + sqrt(xPos(af::span, 2)) * xPos(af::span, 2) - xPos(af::span, 3) * xPos(af::span, 3);

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

af::array myFunctionExC(const af::array& xPos)
{
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code that calculates the correct value of the function vvvvvvvv//

	// It is necessary to fix the calculated and returned value of this function.

	af::array fValue = xPos(af::span, 0) * xPos(af::span, 0) * xPos(af::span, 1) * xPos(af::span, 1) * xPos(af::span, 2) * xPos(af::span, 2) * (cos(xPos(af::span, 3) / 3.0) + 1.0);

	//^^^^^^^ Write above the code that calculates the correct value of the function  ^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return fValue;
}

T monteCarloIntegralN(
	af::array(*myFunction)(const af::array&),
	const af::array& xMin, const af::array& xMax, T absTol, long long int nBatch, long long int itMax)
{
	//This algorithm calculates an ND integral using the Monte Carlo method.
	//
	// myFunction: a function pointer of a function taking 1 input argument. This
	//	"input" argument need to be an array of size
	//  nBatch by nDim 
	//  so that the function is able to evaluate many points at each call.
	//	The value nDim is the number of dimension and can be obtain from the
	//  number of elements in the row arrays of the bounding box.
	// xMin, xMax : two row arrays defining the bounding box of the Monte Carlo algorithm;
	// absTol: Absolute tolerance between the computed integralApproximationand the "last" computed integralApproximation;
	// nBatch: Number of function evaluations before testing the stopping criterion;
	// itMax: Maximum of number of iterations of the algorithm, which would lead to a maximum of nBatch* itMax function evaluations.

	// Reset ArrayFire random number generator
	af::setSeed(0);

	// "First" estimate of the integralApproximation
	T integralApproximation = std::numeric_limits<T>::quiet_NaN();

	////////vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	////////vvvvvvv Write below the code that calculates the integral approximation  vvvvvvvv//

	T volume = af::product(xMax - xMin).scalar<T>();
	//std::cout << "volume  " << volume << std::endl;

	T fSum = 0;

	//af_print(af::randu(nBatch, 3, type::TYPE_AF<T>()));

	//af_print(af::tile(xMin, nBatch,1,  1, 1));


	// Iterative process to compute the integral approximation
	long long int it = 1;
	long long int nDim = xMin.dims(1);

	//std::cout << nDim << std::endl;
	//std::cout << xMin.dims(0) << std::endl;
	while (true)
	{
		//std::cout << "iteration:  " << it<< std::endl;
		T integralApproximation_last = integralApproximation;

		af::array xMinTil = af::tile(xMin, nBatch, 1, 1, 1);
		af::array xMaxTil = af::tile(xMax, nBatch, 1, 1, 1);
	

		af::array xPos = xMinTil + (xMaxTil - xMinTil) * af::randu(nBatch, nDim, type::TYPE_AF<T>());

		//af_print(xPos(af::span, 1));
		

		af::array fValue = myFunction(xPos);

		//std::cout << "=======fValue=======" << std::endl;
		//af_print(fValue);

		fSum = fSum + af::sum<T>(fValue);

		//af_print(fSum);

		long long int nTotal = nBatch * it;

		T fMean = fSum / T(nTotal);

		//std::cout << "fMean  " << fMean << std::endl;
		integralApproximation = volume * fMean;

		if (std::abs((integralApproximation - integralApproximation_last)) < absTol)
			break;

		if (it >= itMax)
			break;

		it += 1;
	}



	//^^^^^^^ Write above the code that calculates the integral approximation  ^^^^^^^^//
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	return integralApproximation;
}

int main(int argc, char *argv[])
{
	af::setBackend(AF_BACKEND_CPU);

	T absTol = (T)0.1;
	long long int nBatch = 10000;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Note that the content below does not require any modification.
	// Only the definition and content of the functions above the main program may require some modifications.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// 1. Calculation of the integral of x.^2-y.^2 with an 2D Monte Carlo algorithm over a given bounding box
	T myResultEx1 = monteCarloIntegral2(&myFunctionEx1, -1, 3, -4, 5, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultEx1: " << std::setw(20) << std::setprecision(16) << myResultEx1 << " : vs : " << -168 << std::endl;

	// 2. Calculation of PI with an 2D Monte Carlo algorithm (using the surface of a disk S=pi*R^2)
	T myResultEx2 = monteCarloIntegral2(&myFunctionEx2, -0.5, 0.5, -0.5, 0.5, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultEx2: " << std::setw(20) << std::setprecision(16) << myResultEx2 << " : vs : " << M_PI << std::endl;

	// 3. Calculation of the integral of x.^2-y.^2+z.^3 with an 3D Monte Carlo algorithm over a given bounding box
	T myResultEx3 = monteCarloIntegral3(&myFunctionEx3, -1, 3, -4, 5, -7, 8, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultEx3: " << std::setw(20) << std::setprecision(16) << myResultEx3 << " : vs : " << 12735 << std::endl;

	// 4. Calculation of PI with an 3D Monte Carlo algorithm (using the volume of a sphere V=4/3*pi*R^3)
	T myResultEx4 = monteCarloIntegral3(&myFunctionEx4, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultEx4: " << std::setw(20) << std::setprecision(16) << myResultEx4 << " : vs : " << M_PI << std::endl;

	// 5. Calculation of the volume of a torus a^2 - (c-sqrt(x^2+y^2))^2 > z^2 with a=c=1 and an 3D Monte Carlo algorithm over a given bounding box
	T myResultEx5 = monteCarloIntegral3(&myFunctionEx5, -3, 3, -3, 3, -1, 1, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultEx5: " << std::setw(20) << std::setprecision(16) << myResultEx5 << " : vs : " << (T)2 * M_PI * M_PI * (T)1 * (T)1 * (T)1 << std::endl;

	// 6. Calculation of the integral of x.^2-y.^2+z.^3-t^2 with an 4D Monte Carlo algorithm over a given bounding box
	T myResultEx6 = monteCarloIntegral4(&myFunctionEx6, -1, 3, -1, 2, 1, 2, -1, 1, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultEx6: " << std::setw(20) << std::setprecision(16) << myResultEx6 << " : vs : " << (T)68.70580079512686 << std::endl;

	// 7. Calculation of the time-averaged temperature inside a torus with an 4D Monte Carlo algorithm
	T tEnd = (T)6 * M_PI;
	T myResultEx7 = monteCarloIntegral4(&myFunctionEx7, -3, 3, -3, 3, -1, 1, 0, tEnd, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultEx7: " << std::setw(20) << std::setprecision(16) << myResultEx7 / tEnd << std::endl;

	// A. Calculation of the integral of x.^2-y.^2+z.^3 with an ND Monte-Carlo algoritm over a given bounding box
	af::array xMinA = af::constant(0, 1, 3, type::TYPE_AF<T>());
	xMinA(0) = -1; xMinA(1) = -4; xMinA(2) = -7;
	af::array xMaxA = af::constant(0, 1, 3, type::TYPE_AF<T>());
	xMaxA(0) = 3; xMaxA(1) = 5; xMaxA(2) = 8;
	T myResultExA = monteCarloIntegralN(&myFunctionExA, xMinA, xMaxA, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultExA: " << std::setw(20) << std::setprecision(16) << myResultExA << " : vs : " << 12735 << std::endl;

	// B. Calculation of the integral of x.^2-y.^2+z.^3-t^2 with an ND Monte-Carlo algoritm over a given bounding box
	af::array xMinB = af::constant(0, 1, 4, type::TYPE_AF<T>());
	xMinB(0) = -1; xMinB(1) = -1; xMinB(2) = 1; xMinB(3) = -1;
	af::array xMaxB = af::constant(0, 1, 4, type::TYPE_AF<T>());
	xMaxB(0) = 3; xMaxB(1) = 2; xMaxB(2) = 2; xMaxB(3) = 1;
	T myResultExB = monteCarloIntegralN(&myFunctionExB, xMinB, xMaxB, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultExB: " << std::setw(20) << std::setprecision(16) << myResultExB << " : vs : " << (T)68.70580079512686 << std::endl;

	// C. Calculation of the time-averaged temperature inside a torus with an ND Monte-Carlo algoritm
	af::array xMinC = af::constant(0, 1, 4, type::TYPE_AF<T>());
	xMinC(0) = -3; xMinC(1) = -3; xMinC(2) = -1; xMinC(3) = 0;
	af::array xMaxC = af::constant(0, 1, 4, type::TYPE_AF<T>());
	xMaxC(0) = 3; xMaxC(1) = 3; xMaxC(2) = 1; xMaxC(3) = tEnd;
	T myResultExC = monteCarloIntegralN(&myFunctionExC, xMinC, xMaxC, absTol, nBatch, std::numeric_limits<long long int>::max());
	std::cout << "myResultExC: " << std::setw(20) << std::setprecision(16) << myResultExC / tEnd << std::endl;

	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	//vvvvvvv Write below the code you want for debugging or for any other reason vvvvvvvv//
	int a;
	std::cin >> a;

	return 0;
}
