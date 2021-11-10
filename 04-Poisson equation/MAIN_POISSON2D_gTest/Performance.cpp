#include "Performance.h"
#include "../MAIN_POISSON2D/DenseSolver.h"
#include "../MAIN_POISSON2D/JacobiSolver.h"

Performance::Performance()
{
	std::cout << "Performance" << std::endl;

	std::string filename1 = "PerformanceLog.plt";
	fopen_s(&file1, filename1.c_str(), "a");

}

Performance::~Performance() 
{
	std::cout << "~Performance" << std::endl;	
	
	fclose(file1);

}

void Performance::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void Performance::TearDown()
{
	std::cout << "TearDown" << std::endl;

}


INSTANTIATE_TEST_CASE_P(Performance, Performance, ::testing::Values(AF_BACKEND_CPU, AF_BACKEND_OPENCL, AF_BACKEND_CUDA) );


//TEST_P(Performance, denseSolver) {
//
//	int nx = 150;
//	int ny = 150;
//	double xMin = 0.0;
//	double xMax = 1.0;
//	double yMin = 0.0;
//	double yMax = 1.0;
//
//
//	af::Backend backend = GetParam();
//
//	af::setBackend(backend);
//	af::setDevice(0);
//	af::info(); std::cout << std::endl;
//
//	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
//	U.init();
//	U.setBC();
//
//
//	DenseSolver solver = DenseSolver(U);
//
//
//	af::timer t1 = af::timer::start();
//
//	solver.solve();
//
//	std::string back;
//	if (backend == AF_BACKEND_CPU) {
//		back = "CPU_backend";
//	}
//	else if (backend == AF_BACKEND_OPENCL) {
//		back = "openCL_backend";
//	}
//	else if (backend == AF_BACKEND_CUDA) {
//		back = "CUDA_backend";
//	}
//
//	fprintf(file1, "%s  %s  %d   %g \n", "DenseSolver", back.c_str(), nx, af::timer::stop(t1));
//
//
//
//
//
//	FAIL();
//}


TEST_P(Performance, jacobiSolver) {

	int nx = 100;
	int ny = 100;
	double xMin = 0.0;
	double xMax = 1.0;
	double yMin = 0.0;
	double yMax = 1.0;


	af::Backend backend = GetParam();

	af::setBackend(backend);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();


	JacobiSolver solver = JacobiSolver(U);


	af::timer t1 = af::timer::start();
	solver.solve();


	std::string back;
	if (backend == AF_BACKEND_CPU) {
		back = "CPU_backend";
	}
	else if (backend == AF_BACKEND_OPENCL) {
		back = "openCL_backend";
	}
	else if (backend == AF_BACKEND_CUDA) {
		back = "CUDA_backend";
	}

	fprintf(file1, "%s  %s  %d   %g \n", "JacobiSolver", back.c_str(), nx, af::timer::stop(t1));




	FAIL();
}










