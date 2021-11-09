#include <iostream>
#include "PoisonSolver.h"

PoisonSolver::PoisonSolver()
{
}


PoisonSolver::~PoisonSolver()
{
}


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
	nUnknown = nx * ny;

	RHS = new double [nUnknown];

	nCol = nUnknown;
	nRow = nUnknown;
	nNonZero = 2 * (nx * (ny - 1) + ny * (nx - 1)) + nUnknown;
	colIdx = new int[nNonZero];
	values = new double[nNonZero];
	nNonZeroRow = new int[nUnknown + 1];
	for (int iRow = 0; iRow < nUnknown + 1; iRow++) {
		nNonZeroRow[iRow] = 0;
	}

}


void PoisonSolver::creatCoeffMatrix()
{
	int iNeib[4];

	int iCol = 0;
	int iRow = 0;
	int iNonZero = 0;

	nNonZeroRow[0] = 0;
	for (int iNode = 0; iNode < nNodeField; iNode++) {

		if (field.isBound[iNode] == true) continue;

		int nNZRowCount = 0;


		values[iNonZero] = -2 * (dx * dx + dy * dy);

		colIdx[iNonZero] = field.nodeInx[iNode];

		nNZRowCount++;

		iNonZero++;


		field.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];


		if (field.isBound[Tneib] == false) {
			values[iNonZero] = dx * dx;

			colIdx[iNonZero] = field.nodeInx[Tneib];

			nNZRowCount++;

			iNonZero++;
		}

		if (field.isBound[Rneib] == false) {
			values[iNonZero] = dy * dy;

			colIdx[iNonZero] = field.nodeInx[Rneib];

			nNZRowCount++;

			iNonZero++;
		}

		if (field.isBound[Bneib] == false) {
			values[iNonZero] = dx * dx;

			colIdx[iNonZero] = field.nodeInx[Bneib];

			nNZRowCount++;

			iNonZero++;
		}

		if (field.isBound[Lneib] == false) {
			values[iNonZero] = dy * dy;

			colIdx[iNonZero] = field.nodeInx[Lneib];

			nNZRowCount++;

			iNonZero++;
		}

		nNonZeroRow[iRow + 1] = nNonZeroRow[iRow] + nNZRowCount;
		iRow++;
		
	}

}


void PoisonSolver::trnsfrCoefMatrx2AF()
{

	af::array values_af = af::array(af::dim4(nNonZero), values);
	af::array nNonZeroRow_af = af::array(af::dim4(nRow + 1), nNonZeroRow);
	af::array colIdx_af = af::array(af::dim4(nNonZero), colIdx);

	CoefMatrxSprs = af::sparse(nRow, nCol, values_af, nNonZeroRow_af, colIdx_af, AF_STORAGE_CSR);

}


void PoisonSolver::creatRigtHandSid()
{
	int iNeib[4];

	int iUnknown = 0;
	for (int iNode = 0; iNode < nNodeField; iNode++) {
		if (field.isBound[iNode] == true) continue;

		double rhs = 0.0;

		field.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];

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

		RHS[iUnknown] = -rhs + (dx * dx * dy * dy) * function(x, y);
		iUnknown++;

	}

	
}


void PoisonSolver::back2Field(af::array A_d)
{

	double* A_h = new double[nUnknown];

	A_h = A_d.host<double>();

	int cnt = 0;
	for (int iNode = 0; iNode < nNodeField; iNode++) {

		if (field.isBound[iNode] == false) {
			field.var[iNode] = A_h[cnt];
			cnt++;
		}

	}

}





