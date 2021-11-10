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

//This method creat a diagonal sparse matrix on arrayfire
af::array JacobiSolver::createDiagonalSpars(double val, int nUnknown)
{
	//All the valuse set to the input variable
	double *values = new double[nUnknown];
	for (int i = 0; i < nUnknown; i++) {
		values[i] = val;
	}

	//As we have valuse of the diagonal the index of collomns are equal to the node index in the matrix
	int *colIdx = new int[nUnknown];
	for (int i = 0; i < nUnknown; i++) {
		colIdx[i] = i;
	}

	//We know that at each row we have just one element but, in arrayfire this parameter set based on the arrayfire doc
	int* nNonZeroRow = new int[nUnknown + 1];
	nNonZeroRow[0] = 0;
	for (int i = 0; i < nUnknown; i++) {
		nNonZeroRow[i+1] = i + 1;
	}

	//reating the sparse matrix on arrayfire
	af::array values_af = af::array(af::dim4(nUnknown), values);
	af::array nNonZeroRow_af = af::array(af::dim4(nUnknown + 1), nNonZeroRow);
	af::array colIdx_af = af::array(af::dim4(nUnknown), colIdx);
	af::array diagonalMatyrix = af::sparse(nUnknown, nUnknown, values_af, nNonZeroRow_af, colIdx_af, AF_STORAGE_CSR);

	return diagonalMatyrix;
}


void JacobiSolver::solve()
{

	//Creating RHS of Poison's equation
	creatRigtHandSid();

	//Creating coefficient matrix in naive C++
	creatCoeffMatrix();

	//get the value of element on the diagonal 
	double diagonalVal = values[0];

	//coefficient matrix to arayfire
	trnsfrCoefMatrx2AF();

	//Trasfering RHS to arrayfire
	af::array b = af::array(nUnknown, RHS);

	//Set the initial value in jacobi method (However, this values could be initilize more accurate for better convergence)
	af::array x0 = b;

	//Trasfering RHS to arrayfire
	af::array I = createDiagonalSpars(diagonalVal, nUnknown);
	af::array LU = CoefMatrxSprs - I;


	//solving the system of equation by using arrayfire library and jacobi method
	af::array x;
	int it = 0;
	double diff = 10e6;
	while ( it < 10e5 || diff < 10e-4) {

		x = 1 / diagonalVal * (b - af::matmul(LU, x0));

		double diff = af::sum<double>(af::abs(x - x0));

		x0 = x;

		it++;
	}

	//Transfering data from the matrix to the discrete domain
	back2Field(x);
}

