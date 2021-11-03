#include "PerformanceTest.h"

PerformanceTest::PerformanceTest()
{
	std::cout << "PerformanceTest" << std::endl;

	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;
}

PerformanceTest::~PerformanceTest() 
{
	std::cout << "~PerformanceTest" << std::endl;
}

void PerformanceTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void PerformanceTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}

void mainFunc(int Nx, int Ny, double xMin, double xMax, double yMin, double yMax){

	Field U = Field(Nx, Ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();


	LinearSys linearSys = LinearSys(U);
	linearSys.creatCoeffMatrix();
	linearSys.creatRigtHandSid();
	linearSys.solve();


	AnalyticalSolution Ua = AnalyticalSolution(Nx, Ny, xMin, xMax, yMin, yMax);


	for (int j = 0; j < Nx * Ny; j++) {
		printf("error[%d] = %.8f  Numeric:Anaytic: %.8f   %.8f \n", j, U.var[j] - Ua.var[j], U.var[j], Ua.var[j]);
	}
	double error = 0.0;
	for (int j = 0; j < Nx * Ny/*U.get_nNode()*/; j++) {
		error += (U.var[j] - Ua.var[j]) * (U.var[j] - Ua.var[j]);
	}
	printf("error is = %.8f  \n", sqrt(error));

}





int nx_1[] = { 5, 10, 20, 40, 80 };
int ny_1[] = { 5, 10, 20, 40, 80 };

INSTANTIATE_TEST_CASE_P(2DProgram, PerformanceTest, ::testing::Values(0, 1, 2, 3, 4) );


TEST_P(PerformanceTest, Test_group1) {

	int testNumber = GetParam();
	int Nx = nx_1[testNumber];
	int Ny = ny_1[testNumber];
	double xMin = 0.0;
	double xMax = 1.0;
	double yMin = 0.0;
	double yMax = 1.0;


	mainFunc(Nx, Ny, xMin, xMax, yMin, yMax);


	FAIL();
}











