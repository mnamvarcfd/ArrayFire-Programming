
#include "PoisonSolver.h"

PoisonSolver::PoisonSolver()
{
}


PoisonSolver::~PoisonSolver()
{
}

//Right hand side function of Poison's equation
double PoisonSolver::function(double x, double y) {
	return -sin(x) - cos(y);
}


PoisonSolver::PoisonSolver(Field & field_)
{
	field = field_;

	dx = field.get_dx();
	dy = field.get_dy();

	nNodeField = field.get_nNode();

	int nx = field.get_nx() - 2;
	int ny = field.get_ny() - 2;
	nUnknown = nx * ny; // i.e. it is the number of non boundary nodes

	RHS = new double [nUnknown]; //An array to store RHS of equation

	nCol = nUnknown;
	nRow = nUnknown;
	nNonZero = 2 * (nx * (ny - 1) + ny * (nx - 1)) + nUnknown; 
	colIdx = new int[nNonZero]; //Refere to SRC format for this
	values = new double[nNonZero]; //Refere to SRC format for this
	nNonZeroRow = new int[nUnknown + 1]; //Refere to SRC format for this
	for (int iRow = 0; iRow < nUnknown + 1; iRow++) {
		nNonZeroRow[iRow] = 0;
	}

}

//this method create the coefficient matrix
void PoisonSolver::creatCoeffMatrix()
{
	//An array to store the index of neighborin node based on the stencil used in the report
	int iNeib[4];

	//Initializing variables
	int iCol = 0;
	int iRow = 0;
	int iNonZero = 0;

	//Based on the CSR format the first value should be zero
	nNonZeroRow[0] = 0;

	//Traversing all the nodes in the discrete domain
	for (int iNode = 0; iNode < nNodeField; iNode++) {

		//Just non-boundary nodes is an unknown so the matrix is constructed for them
		if (field.isBound[iNode] == true) continue;

		//Based on the CRS format we need to count the number of non-zero element of a row
		int nNZRowCount = 0;

		//These 4 steps are done for any nodes and its non-boundary neighbors
		//1-Store the value of matrix based on the discritization in the report
		values[iNonZero] = -2 * (dx * dx + dy * dy);

		//2-Find the colomn number of the above value in the matrix (A renumbering of nodes done in the Grid class)
		colIdx[iNonZero] = field.nodeInx[iNode];

		//3-Count number of non-zero elements of each row
		nNZRowCount++;

		//4-Count number of non-zero elements
		iNonZero++;

		//Find the index of all the neighboring nodes
		field.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];

		//if the "top-neighboring node" is non-boundary 
		if (field.isBound[Tneib] == false) {
			//1-
			values[iNonZero] = dx * dx;

			//2-
			colIdx[iNonZero] = field.nodeInx[Tneib];

			//3-
			nNZRowCount++;

			//4-
			iNonZero++;
		}
		//if the "right-neighboring node" is non-boundary 
		if (field.isBound[Rneib] == false) {
			//1-
			values[iNonZero] = dy * dy;

			//2-
			colIdx[iNonZero] = field.nodeInx[Rneib];

			//3-
			nNZRowCount++;

			//4-
			iNonZero++;
		}
		//if the "bottom-neighboring node" is non-boundary 
		if (field.isBound[Bneib] == false) {
			//1-
			values[iNonZero] = dx * dx;

			//2-
			colIdx[iNonZero] = field.nodeInx[Bneib];

			//3-
			nNZRowCount++;

			//4-
			iNonZero++;
		}
		//if the "left-neighboring node" is non-boundary 
		if (field.isBound[Lneib] == false) {
			//1-
			values[iNonZero] = dy * dy;

			//2-
			colIdx[iNonZero] = field.nodeInx[Lneib];

			//3-
			nNZRowCount++;

			//4-
			iNonZero++;
		}

		//Store number of non-zero elements of a row
		nNonZeroRow[iRow + 1] = nNonZeroRow[iRow] + nNZRowCount;
		iRow++;
		
	}

}


void PoisonSolver::trnsfrCoefMatrx2AF()
{
	//Create sparse matrix in arrayfire
	af::array values_af = af::array(af::dim4(nNonZero), values);
	af::array nNonZeroRow_af = af::array(af::dim4(nRow + 1), nNonZeroRow);
	af::array colIdx_af = af::array(af::dim4(nNonZero), colIdx);
	CoefMatrxSprs = af::sparse(nRow, nCol, values_af, nNonZeroRow_af, colIdx_af, AF_STORAGE_CSR);

	//Free unnecessary allocate memory
	delete[] values;
	delete[] nNonZeroRow;
	delete[] colIdx;
}


void PoisonSolver::creatRigtHandSid()
{
	//An array to store the index of neighborin node based on the stencil used in the report
	int iNeib[4];

	int iUnknown = 0;

	//Traversing all the nodes in the discrete domain
	for (int iNode = 0; iNode < nNodeField; iNode++) {

		//Just boundary nodes participate to create rifht hand side 
		if (field.isBound[iNode] == true) continue;

		double rhs = 0.0;

		//Find the index of all the neighboring nodes
		field.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];

		//Just if a neighboring node is boundary contribute in the RHS
		if (field.isBound[Tneib] == true) {
			rhs += dx * dx * field.var[Tneib];
		}

		if (field.isBound[Rneib] == true) {
			rhs += dy * dy * field.var[Rneib];
		}

		if (field.isBound[Bneib] == true) {
			rhs += dx * dx * field.var[Bneib];
		}

		if (field.isBound[Lneib] == true) {
			rhs += dy * dy * field.var[Lneib];
		}

		double x = field.get_xVal(iNode);
		double y = field.get_yVal(iNode);

		//Create RHS
		RHS[iUnknown] = -rhs + (dx * dx * dy * dy) * function(x, y);
		iUnknown++;

	}


}


void PoisonSolver::back2Field(af::array A_d)
{
	//transfering the data to the host
	double* A_h = new double[nUnknown];
	A_h = A_d.host<double>();

	//Find the associate node in the grid class to the unknown matrix
	int cnt = 0;
	for (int iNode = 0; iNode < nNodeField; iNode++) {

		if (field.isBound[iNode] == false) {
			field.var[iNode] = A_h[cnt];
			cnt++;
		}

	}

}





