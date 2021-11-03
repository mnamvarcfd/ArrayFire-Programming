#include "Field.h"

Field::Field()
{
}

Field::~Field()
{
}


Field::Field(int nx, int ny, double xMin, double xMax, double yMin, double yMax):Grid(nx, ny, xMin, xMax, yMin, yMax)
{

	var = new double[nNode];

}


void Field::init()
{

	for (int i = 0; i < nNode; i++) {
		var[i] = 0.0;
	}

}


void Field::setBC()
{

	//Bottom side
	for (int i = 0; i < nx; i++) {
		double x = i * dx;
		int iNode = getGlobInx(i, 0);
		var[iNode] = sin(x) + cos(yMin); 
	}

	//Top side
	for (int i = 0; i < nx; i++) {
		double x = i * dx;
		int iNode = getGlobInx(i, ny-1);
		var[iNode] = sin(x) + cos(yMax);
	}

	//Left side
	for (int j = 0; j < ny; j++) {
		double y = j * dy;
		int iNode = getGlobInx(0, j);
		var[iNode] = sin(xMin) + cos(y);
	}

	//Right side
	for (int j = 0; j < ny; j++) {
		double y = j * dy;
		int iNode = getGlobInx(nx-1, j);
		var[iNode] = sin(xMax) + cos(y);
	}

}

void Field::write(std::string fileName) {

	writePLT(fileName, var);

}