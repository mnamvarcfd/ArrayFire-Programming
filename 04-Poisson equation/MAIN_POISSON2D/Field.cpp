#include "Field.h"

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
void Field::init()
{

	for (int i = 0; i < nNode; i++) {
		var[i] = 0.0;
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
		int iNode = getGlobInx(i, ny-1);
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
		int iNode = getGlobInx(nx-1, j);
		double y = get_yVal(iNode);
		var[iNode] = sin(xMax) + cos(y);
	}

}

//write the main variable of the class in a PLT file
void Field::write(std::string fileName) {

	writePLT(fileName, var);

}