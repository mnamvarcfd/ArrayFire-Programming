#include "Grid.h"

Grid::Grid()
{
}

Grid::~Grid()
{
}

Grid::Grid(int nx_, int ny_)
{
	nx = nx_;
	ny = ny_;
	nNode = nx * ny;

	dx = 1.0 / ((double)nx - 1);
	dy = 1.0 / ((double)ny - 1);

}

Grid::Grid(int nx_, int ny_, double xMin_, double xMax_, double yMin_, double yMax_)
{
	nx = nx_;
	ny = ny_;
	nNode = nx * ny;

	xMin = xMin_;
	xMax = xMax_;
	yMin = yMin_;
	yMax = yMax_;

	dx = (xMax - xMin) / ((double)nx - 1);
    dy = (yMax - yMin) / ((double)ny - 1);

}

int Grid::getGlobInx(int i, int j)
{
    int iNode = j * nx + i;
    return iNode;
}

int Grid::getXInx(int iNode)
{
	int i = iNode % ny;
	return i;
}

int Grid::getYInx(int iNode)
{
	int j = iNode / ny;
	return j;
}

void Grid::get_iNeib(int iNode, int* iNeib)
{

	int i = getXInx(iNode);
	int j = getYInx(iNode);

	int ip1 = i + 1;
	int im1 = i - 1;
	int jp1 = j + 1;
	int jm1 = j - 1;

	iNeib[0] = getGlobInx(i, jp1);
	iNeib[1] = getGlobInx(ip1, j);
	iNeib[2] = getGlobInx(i, jm1);
	iNeib[3] = getGlobInx(im1, j);

	if (j == ny - 1) iNeib[0] = -1;
	if (i == nx - 1) iNeib[1] = -1;
	if (j == 0)      iNeib[2] = -1;
	if (i == 0)      iNeib[3] = -1;

}

int Grid::get_nx()
{
    return nx;
}

int Grid::get_ny()
{
    return ny;
}

int Grid::get_nNode()
{
    return nNode;
}

double Grid::get_dx()
{
	return dx;
}

double Grid::get_dy()
{
	return dy;
}

void Grid::writePLT(std::string nameOfFile, double *scalar) {

	std::string fileName = nameOfFile + ".plt";

	FILE* file;
	fopen_s(&file, fileName.c_str(), "w");

	fprintf(file, " TITLE = 'Scalar'\n");
	fprintf(file, " VARIABLES = 'X', 'Y' , 'Scalar'   \n");
	fprintf(file, "ZONE T = \"Scalar zone\", I = %d , J = %d , F = POINT \n", nx, ny);

	for (int iNode = 0; iNode < nNode; iNode++) {
		int i = getXInx(iNode);
		int j = getYInx(iNode);

		double x = i * dx;
		double y = j * dy;

		fprintf(file, "%e  %e  %e \n ", x, y, scalar[iNode]);
	}

	fclose(file);
}



