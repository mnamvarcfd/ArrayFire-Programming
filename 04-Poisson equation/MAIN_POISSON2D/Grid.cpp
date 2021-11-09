#include "Grid.h"

Grid::Grid()
{
}

Grid::~Grid()
{
}

Grid::Grid(int nx_, int ny_, double xMin_, double xMax_, double yMin_, double yMax_)
{
	//defne the required parameters
	nx = nx_;
	ny = ny_;
	nNode = nx * ny;

	xMin = xMin_;
	xMax = xMax_;
	yMin = yMin_;
	yMax = yMax_;

	dx = (xMax - xMin) / ((double)nx - 1);
    dy = (yMax - yMin) / ((double)ny - 1);

	//define if a node is budary or not
	isBound = new bool[nNode];
	findBoundNodes();

	//index of nodes after renumberig
	nodeInx = new int[nNode];
	renumbering();

}

//get global index of nodes
int Grid::getGlobInx(int i, int j)
{
    int iNode = j * nx + i;
    return iNode;
}

int Grid::getXInx(int iNode)
{
	int i = iNode % nx;
	return i;
}

int Grid::getYInx(int iNode)
{
	int j = iNode / nx;
	return j;
}

//get the neighboring node index of a given node
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

double Grid::get_xVal(int iNode)
{
	int i = getXInx(iNode);

	double x = xMin + i * dx;

	return x;
}

double Grid::get_yVal(int iNode)
{
	int j = getYInx(iNode);

	double y = yMin + j * dy;

	return y;
}

//labeling nodes as boundary or non boundary
void Grid::findBoundNodes()
{
	for (int i = 0; i < nNode; i++) {
		isBound[i] = false;
	}

	//Bottom side
	for (int i = 0; i < nx; i++) {
		int iNode = getGlobInx(i, 0);
		isBound[iNode] = true;
	}

	//Top side
	for (int i = 0; i < nx; i++) {
		int iNode = getGlobInx(i, ny - 1);
		isBound[iNode] = true;
	}

	//Left side
	for (int j = 0; j < ny; j++) {
		int iNode = getGlobInx(0, j);
		isBound[iNode] = true;
	}

	//Right side
	for (int j = 0; j < ny; j++) {
		int iNode = getGlobInx(nx - 1, j);
		isBound[iNode] = true;
	}

}

//assign an index to any nodes in the domain. Here the first all non-boundary nodes are indexed after that the boundary nodes.
void Grid::renumbering()
{

	int newInx = 0;
	int oldInx = 0;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {

			if (isBound[oldInx] == false) {
				nodeInx[oldInx] = newInx;
				newInx++;
			}

			oldInx++;

		}
	}

}

//writing an scalar value on the discrete domain (grid)
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

		double x = get_xVal(iNode);
		double y = get_yVal(iNode);

		fprintf(file, "%e  %e  %e \n ", x, y, scalar[iNode]);
	}

	fclose(file);
}



