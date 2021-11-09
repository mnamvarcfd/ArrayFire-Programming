#include "JacobiSolver.h"

JacobiSolver::JacobiSolver()
{
}

JacobiSolver::~JacobiSolver()
{
}

JacobiSolver::JacobiSolver(Field & field_):PoisonSolver(field_)
{

}


af::array JacobiSolver::createDiagonalSpars(double val, int nUnknown)
{
	values = new double[nUnknown];
	for (int i = 0; i < nUnknown; i++) {
		values[i] = val;
	}

	colIdx = new int[nUnknown];
	for (int i = 0; i < nUnknown; i++) {
		colIdx[i] = i;
	}
	nNonZeroRow = new int[nUnknown + 1];
	nNonZeroRow[0] = 0;
	for (int i = 0; i < nUnknown; i++) {
		nNonZeroRow[i+1] = i + 1;
	}


	af::array values_af = af::array(af::dim4(nUnknown), values);
	af::array nNonZeroRow_af = af::array(af::dim4(nUnknown + 1), nNonZeroRow);
	af::array colIdx_af = af::array(af::dim4(nUnknown), colIdx);

	af::array diagonalMatyrix = af::sparse(nUnknown, nUnknown, values_af, nNonZeroRow_af, colIdx_af, AF_STORAGE_CSR);

	return diagonalMatyrix;
}


void JacobiSolver::solve()
{

	creatRigtHandSid();

	creatCoeffMatrix();

	trnsfrCoefMatrx2AF();

	double diagonalVal = values[0];
	

	af::array b = af::array(nUnknown, RHS);

	af::array x0 = b;

	af::array x;

	af::array I = createDiagonalSpars(diagonalVal, nUnknown);
	af::array LU = CoefMatrxSprs - I;


	printf("jacobi solver started....  \n");
	for (int it = 0; it < 10e5; it++) {

		x = 1 / diagonalVal * (b - af::matmul(LU, x0));

		double diff = af::sum<double>(af::abs(x - x0));

		x0 = x;

		//printf("%f  \n",diff);
		if (diff<10e-4)break;
	
	}

	////af_print(x);

	back2Field(x);
}

//void JacobiSolver::solve()
//{
//
//	creatRigtHandSid();
//
//	creatCoeffMatrix();
//
//	trnsfrCoefMatrx2AF();
//
//
//
//	//af::array A1 = CoefMatrxSprs;
//
//	//af::array A2 = CoefMatrxSprs;
//
//	//af::array b = af::array(nUnknown, RHS);
//
//	//af_print(af::inverse(A1) );
//
//
//
//
//	af::array CoefMatrxDens = af::dense(CoefMatrxSprs);
//
//	af::array D = -4.0 * af::identity(nUnknown, nUnknown);
//	//af_print(D);
//
//	af::array U = af::upper(CoefMatrxDens, false) - D;
//	//af_print(U);
//
//	af::array L = af::lower(CoefMatrxDens, false) - D;
//	//af_print(L);
//
//	af::array b = af::array(nUnknown, RHS);
//	//af_print(b);
//
//	af::array x0 = af::array(nUnknown, RHS);
//	//af_print(x0);
//
//	af::array x;
//
//	double omga = 0.2;
//
//	for (int it = 0; it < 500; it++) {
//
//		x = (af::matmul(af::inverse(D + omga * L), omga * b - af::matmul((omga * U + (omga - 1.0) * D), x0)));
//
//		af_print(af::sum(x - x0));
//		x0 = x;
//
//	}
//
//	af_print(x);
//
//	back2Field(x);
//}
