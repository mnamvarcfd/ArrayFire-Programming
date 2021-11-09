#include "DenseSolver.h"

DenseSolver::DenseSolver()
{
}

DenseSolver::~DenseSolver()
{
}

DenseSolver::DenseSolver(Field & field_):PoisonSolver(field_)
{

}


void DenseSolver::solve()
{
	//Creating RHS of Poison's equation
	creatRigtHandSid();

	//Creating coefficient matrix in naive C++
	creatCoeffMatrix();

	//coefficient matrix to arayfire
	trnsfrCoefMatrx2AF();

	//Creatig a dense matrix (Note the coefficient matrix was in sparse shape)
	af::array A = af::dense(CoefMatrxSprs);

	//Trasfering RHS to arrayfire
	af::array b = af::array(nUnknown, RHS);

	//Free unnecessary allocate memory
	delete[] RHS;

	//solving the system of equation by using arrayfire library
	af::array x = af::solve(A, b);

	//Transfering data from the matrix to the discrete domain
	back2Field(x);

}

















//void DenseSolver::solve()
//{
//
//	af::timer t11 = af::timer::start();
//	creatRigtHandSid();
//	af::sync();
//	printf("Elapse time for creating right hand side of Eq.: %g\n", af::timer::stop(t11));
//
//
//	af::timer t12 = af::timer::start();
//	creatCoeffMatrix();
//	af::sync();
//	printf("Elapse time for creating coefficient matrix: %g\n", af::timer::stop(t12));
//
//
//	af::timer t13 = af::timer::start();
//	trnsfrCoefMatrx2AF();
//	af::sync();
//	printf("Elapse time for transfering coefficient matrix to AF: %g\n", af::timer::stop(t13));
//
//
//	af::timer t2 = af::timer::start();
//	af::array A = af::dense(CoefMatrxSprs);
//	af::sync();
//	printf("Elapse time for creating sparse matrix on AF: %g\n", af::timer::stop(t2));
//
//
//	af::timer t3 = af::timer::start();
//	af::array b = af::array(nUnknown, RHS);
//	af::sync();
//	printf("Elapse time for transferring right hand side to AF: %g\n", af::timer::stop(t3));
//
//	size_t alloc_bytes, alloc_buffers;
//	size_t lock_bytes, lock_buffers;
//	af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
//	std::cout << "Used Memory in Mb: " << alloc_bytes / (double)1024 / (double)1024 << std::endl;
//
//	af::timer t4 = af::timer::start();
//	af::array x = af::solve(A, b);
//	af::sync();
//	printf("Elapse time for solving Sys. of Eq. on AF: %g\n", af::timer::stop(t4));
//
//
//	af::timer t5 = af::timer::start();
//	back2Field(x);
//	af::sync();
//	printf("Elapse time for transfering results to host: %g\n", af::timer::stop(t5));
//
//}