#include "LinearSys.h"
#include <iostream>

LinearSys::LinearSys()
{
}

LinearSys::~LinearSys()
{
}

LinearSys::LinearSys(Field & field_)
{

	field = field_;
	dx = field.get_dx();
	dy = field.get_dy();
	nNodefield = field.get_nNode();

	nx = field.get_nx() - 2;
	ny = field.get_ny() - 2;

	nUnknown = nx * ny;
	nNode = nx * ny;
	coeff = new double[nUnknown * nUnknown];

	initCoeffMatrix();

	RHS = new double [nUnknown];


	grid = Grid(nx, ny);
	


}

double LinearSys::function(double x, double y) {
	return -sin(x) - cos(y);
}

void LinearSys::creatCoeffMatrixSpars()
{
	Grid matrix = Grid(nNode, nNode);

	int iNeib[4];
	int nodeIdx;

	for (int iNode = 0; iNode < nNode; iNode++) {

		nodeIdx = matrix.getGlobInx(iNode, iNode);
		coeff[nodeIdx] = -2 * (dx * dx + dy * dy);
		//printf("nodeIdx: ============ %d \n", nodeIdx);

		grid.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];

		//if (iNode == 0)printf("Bneib: %d ==== %d ==== %d ==== %d \n", Tneib, Rneib, Bneib, Lneib);

		if (Tneib != -1) {
			nodeIdx = matrix.getGlobInx(iNode, Tneib);
			coeff[nodeIdx] = dx * dx;
		}

		if (Rneib != -1) {
			nodeIdx = matrix.getGlobInx(iNode, Rneib);
			coeff[nodeIdx] = dy * dy;
		}

		if (Bneib != -1) {
			nodeIdx = matrix.getGlobInx(iNode, Bneib);
			coeff[nodeIdx] = dx * dx;
		}

		if (Lneib != -1) {
			nodeIdx = matrix.getGlobInx(iNode, Lneib);
			coeff[nodeIdx] = dy * dy;
		}

	}


	A = af::array(nNode, nNode, coeff);

}

void LinearSys::creatCoeffMatrix()
{
	int nx = field.get_nx() - 2;
	int ny = field.get_ny() - 2;

	nCol = nUnknown;
	nRow = nUnknown;
	nNonZero = 2 * ( nx * (ny - 1) + ny * (nx - 1) ) + nUnknown;
	colIdx = new int[nNonZero];
	values = new double[nNonZero];
	nNonZeroRow = new int[nUnknown +1];
	for (int iRow = 0; iRow < nUnknown; iRow++) {
		nNonZeroRow[iRow] = 0;
	}


	int iNeib[4];

	int iNonZero = 0;
	for (int iNode = 0; iNode < nNodefield; iNode++) {

		if (field.isBound[iNode] == true) continue;

		int iNonZeroRow = 0;

		values[iNonZero] = -2 * (dx * dx + dy * dy);
		colIdx[iNonZero] = iNode;
		iNonZero++;
		iNonZeroRow++;

		field.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];


		if (field.isBound[Tneib] == false) {
			values[iNonZero] = dx * dx;
			colIdx[iNonZero] = Tneib;
			iNonZero++;
			iNonZeroRow++;
		}

		if (field.isBound[Rneib] == false) {
			values[iNonZero] = dy * dy;
			colIdx[iNonZero] = Rneib;
			iNonZero++;
			iNonZeroRow++;
		}

		if (field.isBound[Bneib] == false) {
			values[iNonZero] = dx * dx;
			colIdx[iNonZero] = Bneib;
			iNonZero++;
			iNonZeroRow++;
		}

		if (field.isBound[Lneib] == false) {
			values[iNonZero] = dy * dy;
			colIdx[iNonZero] = Lneib;
			iNonZero++;
			iNonZeroRow++;
		}

		nNonZeroRow[iNode + 1] = nNonZeroRow[iNode] + iNonZeroRow;

	}

	////printf("nNonZero: %d \n", nNonZero);
	////for (int i = 0; i < nNode + 1; i++) {
	////	printf("nNonZeroRow[%d]:    %d\n", i, nNonZeroRow[i]);
	////}
	////for (int i = 0; i < nNonZero; i++) {
	////	printf("values[%d]: %d ---->  %.8f \n", i, colIdx[i] , values[i]);
	////}


	//af::array vals = af::array(af::dim4(nNonZero), values);
	//af::array row_ptr = af::array(af::dim4(nRow + 1), nNonZeroRow);
	//af::array col_idx = af::array(af::dim4(nNonZero), colIdx);
	//// Create sparse array (CSR) from af::arrays containing values,
	//// row pointers, and column indices.
	//af::array sparse = af::sparse(nRow, nCol, vals, row_ptr, col_idx, AF_STORAGE_CSR);
	//// sparse
	////     values:  [ 5.0, 8.0, 3.0, 6.0 ]
	////     row_ptr: [ 0, 0, 2, 3, 4 ]
	////     col_idx: [ 0, 1, 2, 1 ]

	////af_print(sparse);

	//A = af::dense(sparse);

	////af_print(A);

}


//void LinearSys::creatCoeffMatrixSpars()
//{
//	int nx = field.get_nx() - 2;
//	int ny = field.get_ny() - 2;
//
//	int nCol = nNode;
//	int nRow = nNode;
//	int nNonZero = 2 * (nx * (ny - 1) + ny * (nx - 1)) + nNode;
//	int* nNonZeroRow = new int[nNode + 1];
//	int* colIdx = new int[nNonZero];
//	double* values = new double[nNonZero];
//	for (int iRow = 0; iRow < nNode; iRow++) {
//		nNonZeroRow[iRow] = 0;
//	}
//
//
//
//	Grid matrix = Grid(nNode, nNode);
//
//	int iNeib[4];
//
//	int iNonZero = 0;
//	for (int iNode = 0; iNode < nNode; iNode++) {
//
//		int iNonZeroRow = 0;
//
//		values[iNonZero] = -2 * (dx * dx + dy * dy);
//		colIdx[iNonZero] = iNode;
//		iNonZero++;
//		iNonZeroRow++;
//
//		grid.get_iNeib(iNode, iNeib);
//		int Tneib = iNeib[0];
//		int Rneib = iNeib[1];
//		int Bneib = iNeib[2];
//		int Lneib = iNeib[3];
//
//
//		if (Tneib != -1) {
//			//nodeIdx = matrix.getGlobInx(iNode, Tneib);
//			//coeff[nodeIdx] = dx * dx;
//
//			values[iNonZero] = dx * dx;
//			colIdx[iNonZero] = Tneib;
//			iNonZero++;
//			iNonZeroRow++;
//		}
//
//		if (Rneib != -1) {
//			//nodeIdx = matrix.getGlobInx(iNode, Rneib);
//			//coeff[nodeIdx] = dy * dy;
//
//			values[iNonZero] = dy * dy;
//			colIdx[iNonZero] = Rneib;
//			iNonZero++;
//			iNonZeroRow++;
//		}
//
//		if (Bneib != -1) {
//			//nodeIdx = matrix.getGlobInx(iNode, Bneib);
//			//coeff[nodeIdx] = dx * dx;
//
//			values[iNonZero] = dx * dx;
//			colIdx[iNonZero] = Bneib;
//			iNonZero++;
//			iNonZeroRow++;
//		}
//
//		if (Lneib != -1) {
//			//nodeIdx = matrix.getGlobInx(iNode, Lneib);
//			//coeff[nodeIdx] = dy * dy;
//
//			values[iNonZero] = dy * dy;
//			colIdx[iNonZero] = Lneib;
//			iNonZero++;
//			iNonZeroRow++;
//		}
//
//		nNonZeroRow[iNode + 1] = nNonZeroRow[iNode] + iNonZeroRow;
//
//	}
//
//	//printf("nNonZero: %d \n", nNonZero);
//	//for (int i = 0; i < nNode + 1; i++) {
//	//	printf("nNonZeroRow[%d]:    %d\n", i, nNonZeroRow[i]);
//	//}
//	//for (int i = 0; i < nNonZero; i++) {
//	//	printf("values[%d]: %d ---->  %.8f \n", i, colIdx[i] , values[i]);
//	//}
//
//
//	af::array vals = af::array(af::dim4(nNonZero), values);
//	af::array row_ptr = af::array(af::dim4(nRow + 1), nNonZeroRow);
//	af::array col_idx = af::array(af::dim4(nNonZero), colIdx);
//	// Create sparse array (CSR) from af::arrays containing values,
//	// row pointers, and column indices.
//	af::array sparse = af::sparse(nRow, nCol, vals, row_ptr, col_idx, AF_STORAGE_CSR);
//	// sparse
//	//     values:  [ 5.0, 8.0, 3.0, 6.0 ]
//	//     row_ptr: [ 0, 0, 2, 3, 4 ]
//	//     col_idx: [ 0, 1, 2, 1 ]
//
//	//af_print(sparse);
//
//	A = af::dense(sparse);
//
//	//af_print(A);
//
//
//
//
//
//}



void LinearSys::trnsfrCoefMatrx2AF()
{

	af::array values_af = af::array(af::dim4(nNonZero), values);
	af::array nNonZeroRow_af = af::array(af::dim4(nRow + 1), nNonZeroRow);
	af::array colIdx_af = af::array(af::dim4(nNonZero), colIdx);
	// Create sparse array (CSR) from af::arrays containing values,
	// row pointers, and column indices.
	CoefMatrxSprs = af::sparse(nRow, nCol, values_af, nNonZeroRow_af, colIdx_af, AF_STORAGE_CSR);
	// sparse
	//     values:  [ 5.0, 8.0, 3.0, 6.0 ]
	//     row_ptr: [ 0, 0, 2, 3, 4 ]
	//     col_idx: [ 0, 1, 2, 1 ]

	//af_print(sparse);

	//////A = af::dense(sparse);

	////////af_print(A);

}


void LinearSys::creatRigtHandSid()
{
	int iNeib[4];

	int iNodeInMarix = 0;
	for (int iNode = 0; iNode < nNodefield; iNode++) {

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


		RHS[iNodeInMarix] = -rhs + (dx * dx * dy * dy) * function(x, y);
		iNodeInMarix++;

	}

}


void LinearSys::back2Field(af::array A_d)
{

	double* A_h = new double[nUnknown];

	A_h = A_d.host<double>();

	for (int iNode = 0; iNode < nNodefield; iNode++) {

		if (field.isBound[iNode] == false) {
			field.var[iNode] = A_h[iNode];
		}

	}

}

//
//void LinearSys::creatRigtHandSid()
//{
//	int iNeib[4];
//
//	for (int iNode = 0; iNode < nNode; iNode++) {
//
//		double rhs = 0.0;
//
//		int i = grid.getXInx(iNode);
//		int j = grid.getYInx(iNode);
//
//		grid.get_iNeib(iNode, iNeib);
//		int Tneib = iNeib[0];
//		int Rneib = iNeib[1];
//		int Bneib = iNeib[2];
//		int Lneib = iNeib[3];
//
//		if (Tneib == -1) {
//			int nodeIdx = field.getGlobInx(i+1, field.get_ny()-1);
//			rhs += dx * dx * field.var[nodeIdx];
//		}
//
//		if (Rneib == -1) {
//			int nodeIdx = field.getGlobInx(field.get_nx()-1, j+1);
//			rhs += dy * dy * field.var[nodeIdx];
//		}
//
//		if (Bneib == -1) {
//			int nodeIdx = field.getGlobInx(i+1, 0);
//			rhs += dx * dx * field.var[nodeIdx];
//			//if (iNode == 0)printf("Bneib: %.8f ---->  %.8f \n", field.var[nodeIdx], rhs);
//		}
//
//		if (Lneib == -1) {
//			int nodeIdx = field.getGlobInx(0, j+1);
//			rhs += dy * dy * field.var[nodeIdx];
//			//if (iNode == 0)printf("Lneib: %.8f ---->  %.8f \n", field.var[nodeIdx], rhs);
//		}
//
//
//		double x = (i + 1) * dx;
//		double y = (j + 1) * dy;
//
//
//		 RHS[iNode] = -rhs + (dx * dx * dy * dy) * function(x, y);
//		 //if (iNode == 0)printf("RHS[iNode]: %.8f ---->  %.8f \n", RHS[iNode], function(x, y));
//
//	}
//
//}

void LinearSys::initCoeffMatrix()
{
	int cnt = 0;
	for (int i = 0; i < nNode; i++) {
		for (int j = 0; j < nNode; j++) {
			coeff[cnt] = 0.0;
			cnt++;
		}
	}

}

void LinearSys::writeCoeffMatrix(std::string fileName) {

	Grid matrix = Grid(nNode, nNode);
	std::string filename = fileName + ".plt";

	FILE* file;
	fopen_s(&file, filename.c_str(), "w");

	fprintf(file, " TITLE = 'Scalar'\n");
	fprintf(file, " VARIABLES = 'X', 'Y' , 'Scalar'   \n");
	fprintf(file, "ZONE T = \"Scalar zone\", I = %d , J = %d , F = POINT \n", nNode, nNode);

	for (int i = 0; i < nNode; i++) {
		for (int j = 0; j < nNode; j++) {
			int nodeIdx = matrix.getGlobInx(i, j);
			fprintf(file, "%e  %e  %e \n ", (double)i, (double)nNode-j , coeff[nodeIdx]);
		}
	}

	fclose(file);
}


void LinearSys::solve()
{
	creatRigtHandSid();

	creatCoeffMatrix();

	trnsfrCoefMatrx2AF();

	af::array A = af::dense(CoefMatrxSprs);
	//A = af::array(nNode, nNode, coeff);
	//af_print(A);

	af::array b = af::array(nUnknown, RHS);
	//af_print(b);
	af::array x = af::solve(A, b);

	back2Field(x);

}


void LinearSys::solveSpars()
{

	af::array I = af::identity(nNode, nNode);
	//af_print(D);

	af::array U = af::upper(A, true) - I;
	af_print(U);

	af::array L = af::lower(A, true) - I;
	//af_print(L);

	af::array D = -4.0 * I;
	//af_print(D);


	af::array b = af::array(nNode, RHS);
	//af_print(b);

	af::array x0 = af::array(nNode, RHS);
	//af_print(x0);

	af::array x = af::array(nNode, RHS);

	double omga = 0.2;

	for (int it = 0; it < 500; it++) {

		x = (af::matmul( af::inverse(D + omga * L) , omga*b - af::matmul( (omga*U + (omga-1.0)*D),x0 )  )     );

		af_print(af::sum(x-x0));
		x0 = x;

	}

	af_print(x);

	back2Field(x);
}







	/*double A[9];
	A[0] = 0.1000;     A[1] = 3.1000;     A[2] = 6.1000;
	A[3] = 1.1000;     A[4] = 4.1000;     A[5] = 7.0000;
	A[6] = 2.0000;     A[7] = 5.0000;     A[8] = 8.0000;
	af::array af_A = af::array(3, 3, A);

	af_print(af_A);

	double b[3];
	b[0] = 21.9000;
	b[1] = 30.7000;
	b[2] = 39.0000;
	af::array af_b = af::array(3, b);

	af_print(af_b);

	af::array x = af::solve(af_A, af_b);

	af_print(x);

	af::array b2 = af::matmul(af_A, x);

	af_print(b2);*/