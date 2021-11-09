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

	af::timer t11 = af::timer::start();
	creatRigtHandSid();
	af::sync();
	printf("Elapse time for creating right hand side of Eq.: %g\n", af::timer::stop(t11));


	af::timer t12 = af::timer::start();
	creatCoeffMatrix();
	af::sync();
	printf("Elapse time for creating coefficient matrix: %g\n", af::timer::stop(t12));


	af::timer t13 = af::timer::start();
	trnsfrCoefMatrx2AF();
	af::sync();
	printf("Elapse time for transfering coefficient matrix to AF: %g\n", af::timer::stop(t13));


	af::timer t2 = af::timer::start();
	af::array A = af::dense(CoefMatrxSprs);
	af::sync();
	printf("Elapse time for creating sparse matrix on AF: %g\n", af::timer::stop(t2));


	af::timer t3 = af::timer::start();
	af::array b = af::array(nUnknown, RHS);
	af::sync();
	printf("Elapse time for transferring right hand side to AF: %g\n", af::timer::stop(t3));

	size_t alloc_bytes, alloc_buffers;
	size_t lock_bytes, lock_buffers;
	af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
	std::cout << "Used Memory in Mb: " << alloc_bytes / (double)1024 / (double)1024 << std::endl;

	af::timer t4 = af::timer::start();
	af::array x = af::solve(A, b);
	af::sync();
	printf("Elapse time for solving Sys. of Eq. on AF: %g\n", af::timer::stop(t4));


	af::timer t5 = af::timer::start();
	back2Field(x);
	af::sync();
	printf("Elapse time for transfering results to host: %g\n", af::timer::stop(t5));

}
//Elapse time for creating right hand side of Eq.: 0.0011968
//Elapse time for creating coefficient matrix : 0.0007243
//Elapse time for transfering coefficient matrix to AF : 0.0002931
//Elapse time for creating sparse matrix on AF : 0.0007495
//Elapse time for transferring right hand side to AF : 0.182144
//Elapse time for solving Sys.of Eq.on AF : 0.0008474
//Elapse time for transfering results to host : 2.11079


//Elapse time for creating right hand side of Eq.: 0.0011476
//Elapse time for creating coefficient matrix : 0.0007274
//Elapse time for transfering coefficient matrix to AF : 0.0040418
//Elapse time for creating sparse matrix on AF : 0.0110339
//Elapse time for transferring right hand side to AF : 0.0005502
//Elapse time for solving Sys.of Eq.on AF : 10.9655
//Elapse time for transfering results to host : 0.0004493