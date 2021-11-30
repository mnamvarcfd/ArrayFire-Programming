#include "Field.h"
#include <iostream>

Field::Field()
{
}

Field::~Field()
{
}


Field::Field(int nx, int ny, double xMin, double xMax, double yMin, double yMax):Grid(nx, ny, xMin, xMax, yMin, yMax)
{
	//allocate the main variable of class
	var = new double[nNode];

}

//initializing the main variable as zero
void Field::init(double value)
{

	for (int i = 0; i < nNode; i++) {
		var[i] = value;
	}

}

//initializing the main variable based on a given function
void Field::init(double Func(double x, double y))
{

	for (int i = 0; i < nNode; i++) {
		double x = get_xVal(i);
		double y = get_yVal(i);
		var[i] = Func(x, y);
	}

}

//initializing the main variable based on a given function
void Field::init(double Func(double x, double y, double t), double t, bool periodic)
{
	double val = 0.0;

	double *varTmp = new double[nNode];

	for (int i = 0; i < nNode; i++) {
		double x = get_xVal(i);
		double y = get_yVal(i);

		if ((x * x + y * y) < 0.25) {
			val = exp(1.0 - 0.25 / (0.25 - x * x - y * y));
		}

		varTmp[i] = val;
	}


	int shiftX = ceil(0.5 * t / get_dx());
	int shiftY = ceil(-0.3 * t / get_dy());
	//std::cout << shiftX<<  "  shiftX=====shiftY  " << shiftY << std::endl;


	for (int i = 0; i < nNode; i++) {

		double iIndex = getXInx(i) + shiftX;
		double jIndex = getYInx(i) + shiftY;


		if (iIndex >(nx-1)) {
			iIndex -= nx;
		}
		else if (iIndex < 0) {
			iIndex += nx;
		}

		if (jIndex > (ny - 1)) {
			jIndex -= ny;
		}
		else if (jIndex < 0) {
			jIndex += ny;
		}


		int iNode = getGlobInx(iIndex, jIndex);
		var[iNode] = varTmp[i];
	}

}


//initializing the main variable based on a given function
void Field::init(double Func(double x, double y), double a, double b, double t)
{
	double* varTmp = new double[nNode];

	for (int i = 0; i < nNode; i++) {
		double x = get_xVal(i);
		double y = get_yVal(i);

		varTmp[i] = Func(x,y);
	}


	int shiftX = /*ceil*/(a * t / get_dx());
	int shiftY = /*ceil*/(b * t / get_dy());

	//std::cout << a * t << "  ===== " << get_dx() << std::endl;


	//std::cout<< shiftX << "  shiftX=====shiftY "<< shiftY << std::endl;


	for (int i = 0; i < nNode; i++) {

		double iIndex = getXInx(i) + shiftX;
		double jIndex = getYInx(i) + shiftY;


		for (int i = 0; i < 10; i++) {
			if (iIndex > (nx - 1)) {
				iIndex -= nx;
			}
			else if (iIndex < 0) {
				iIndex += nx;
			}

			if (jIndex > (ny - 1)) {
				jIndex -= ny;
			}
			else if (jIndex < 0) {
				jIndex += ny;
			}
		}

		int iNode = getGlobInx(iIndex, jIndex);

		var[iNode] = varTmp[i];
	}

}




//initializing the main variable based on a given function
void Field::init(double Func(double x, double y, double t), double t)
{

	for (int i = 0; i < nNode; i++) {
		double x = get_xVal(i);
		double y = get_yVal(i);
		var[i] = Func(x, y, t);
	}

}

//applying the boundary condition
void Field::setBC()
{

	//Bottom side
	for (int i = 0; i < nx; i++) {
		//get the global index of the node
		int iNode = getGlobInx(i, 0);
		//get the x valuse of the node
		double x = get_xVal(iNode);
		//applying the BC on the node
		var[iNode] = sin(x) + cos(yMin);
	}

	//Top side
	for (int i = 0; i < nx; i++) {
		int iNode = getGlobInx(i, ny - 1);
		double x = get_xVal(iNode);
		var[iNode] = sin(x) + cos(yMax);
	}

	//Left side
	for (int j = 0; j < ny; j++) {
		int iNode = getGlobInx(0, j);
		double y = get_yVal(iNode);
		var[iNode] = sin(xMin) + cos(y);
	}

	//Right side
	for (int j = 0; j < ny; j++) {
		int iNode = getGlobInx(nx - 1, j);
		double y = get_yVal(iNode);
		var[iNode] = sin(xMax) + cos(y);
	}

}




//write the main variable of the class in a PLT file
void Field::write(std::string fileName) {

	writePLT(fileName, var);

}