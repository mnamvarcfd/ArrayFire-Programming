#include "ValidationTest.h"
#include "../MAIN_POISSON2D/DenseSolver.h"
#include "../MAIN_POISSON2D/JacobiSolver.h"

ValidationTest::ValidationTest()
{
	std::cout << "ValidationTest" << std::endl;

}

ValidationTest::~ValidationTest() 
{
	std::cout << "~ValidationTest" << std::endl;
}

void ValidationTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void ValidationTest::TearDown()
{
	std::cout << "TearDown" << std::endl;

}

double func(double x, double y) {
	return sin(x) + cos(y);
}


INSTANTIATE_TEST_CASE_P(Validation, ValidationTest, ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14) );


TEST_P(ValidationTest, DenseTest) {

	double pi = 4 * atan(1.0);

	int nxTest[] = { 5, 10, 20, 40, 80,   6, 12, 24, 48,  96,   9, 18, 36, 72, 144 };

	int nyTest[] = { 5, 10, 20, 40, 80,   6, 12, 24, 48,  96,   9, 18, 36, 72, 144 };

	double xMinTest[] = { 0, 0, 0, 0, 0,  -3,-3, -3, -3, -3,   -pi, -pi, -pi,-pi,-pi };
	double xMaxTest[] = { 1, 1, 1, 1, 1,  3 * pi, 3 * pi, 3 * pi,3 * pi,3 * pi,   2, 2, 2, 2, 2 };

	double yMinTest[] = { 0, 0, 0, 0, 0,  3, 3, 3, 3, 3,   -5 * pi, -5 * pi, -5 * pi,-5 * pi,-5 * pi };
	double yMaxTest[] = { 1, 1, 1, 1, 1,  4 * pi, 4 * pi, 4 * pi,4 * pi,4 * pi,   3 * pi, 3 * pi, 3 * pi,3 * pi,3 * pi };



	int testNumber = GetParam();

	int nx = nxTest[testNumber];
	int ny = nyTest[testNumber];
	double xMin = xMinTest[testNumber];
	double xMax = xMaxTest[testNumber];
	double yMin = yMinTest[testNumber];
	double yMax = yMaxTest[testNumber];


	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;



	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();


	DenseSolver denseSolver = DenseSolver(U);
	denseSolver.solve();


	Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ua.init(&func);


	double L2error = 0.0;
	for (int j = 0; j < U.get_nNode(); j++) {
		L2error += (U.var[j] - Ua.var[j]) * (U.var[j] - Ua.var[j]);
	}
	printf("L2error is = %.8f  \n", sqrt(L2error / (nx*ny)));


	//FAIL();
}

TEST_P(ValidationTest, JacobiTest) {

	double pi = 4 * atan(1.0);

	int nxTest[] = { 5, 10, 20, 40, 80,   6, 12, 24, 48,  96,   9, 18, 36, 72, 144 };

	int nyTest[] = { 5, 10, 20, 40, 80,   6, 12, 24, 48,  96,   9, 18, 36, 72, 144 };

	double xMinTest[] = { 0, 0, 0, 0, 0,  -3,-3, -3, -3, -3,   -pi, -pi, -pi,-pi,-pi };
	double xMaxTest[] = { 1, 1, 1, 1, 1,  3 * pi, 3 * pi, 3 * pi,3 * pi,3 * pi,   2, 2, 2, 2, 2 };

	double yMinTest[] = { 0, 0, 0, 0, 0,  3, 3, 3, 3, 3,   -5 * pi, -5 * pi, -5 * pi,-5 * pi,-5 * pi };
	double yMaxTest[] = { 1, 1, 1, 1, 1,  4 * pi, 4 * pi, 4 * pi,4 * pi,4 * pi,   3 * pi, 3 * pi, 3 * pi,3 * pi,3 * pi };



	int testNumber = GetParam();

	int nx = nxTest[testNumber];
	int ny = nyTest[testNumber];
	double xMin = xMinTest[testNumber];
	double xMax = xMaxTest[testNumber];
	double yMin = yMinTest[testNumber];
	double yMax = yMaxTest[testNumber];


	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;



	Field Ud = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ud.init();
	Ud.setBC();


	DenseSolver denseSolver = DenseSolver(Ud);
	denseSolver.solve();


	Field Uj = Field(nx, ny, xMin, xMax, yMin, yMax);
	Uj.init();
	Uj.setBC();
	JacobiSolver jacobiSolver = JacobiSolver(Uj);
	jacobiSolver.solve();

	for (int j = 0; j < Ud.get_nNode(); j++) {
		EXPECT_NEAR(Ud.var[j] , Uj.var[j], 10e-3);
	}


	Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ua.init(&func);


	double L2error = 0.0;
	for (int j = 0; j < Uj.get_nNode(); j++) {
		L2error += (Uj.var[j] - Ua.var[j]) * (Uj.var[j] - Ua.var[j]);
	}
	printf("L2error is = %.8f  \n", sqrt(L2error / (nx*ny)));


	FAIL();
}









