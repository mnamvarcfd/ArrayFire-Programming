#include "PerformanceTest.h"

PerformanceTest::PerformanceTest()
{
	std::cout << "PerformanceTest" << std::endl;
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

double function(double x, double y) 
{
	return -sin(x) - cos(y);
}

void mainFunc(int Nx, int Ny, double xMin, double xMax, double yMin, double yMax){

	Field u = Field(Nx, Ny, xMin, xMax, yMin, yMax);
	u.init();
	u.setBC();


	LinearSys linearSys = LinearSys(u);
	linearSys.creatCoeffMatrix();
	linearSys.creatRigtHandSid(&function);
	linearSys.solve();


	AnalyticalSolution analyticalSolution = AnalyticalSolution(Nx, Ny, xMin, xMax, yMin, yMax);


	double sum = 0.0;
	for (int j = 0; j < Nx * Ny; j++) {
		sum += (u.var[j] - analyticalSolution.var[j]) * (u.var[j] - analyticalSolution.var[j]);
	}
	printf("error is = %.8f  \n", sqrt(sum));

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











