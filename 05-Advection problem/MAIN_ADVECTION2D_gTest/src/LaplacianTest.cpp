#include "LaplacianTest.h"


double lexico(double x, double y) {

	double val = y * 6 + x;

	return val;
}
LaplacianTest::LaplacianTest()
{
	std::cout << "LaplacianTest" << std::endl;

	int nx = 6;
	int ny = 6;
	double xMin = 0;
	double xMax = 5;
	double yMin = 0;
	double yMax = 5;


	field = Field(nx, ny, xMin, xMax, yMin, yMax);
	field.init(&lexico);
	field.write("init");

	stencil = D2Q4Stencil(field);

	laplacian = Laplacian(field, stencil);

	laplacian.solve(1);
	field.write("numericalSolution");
}

LaplacianTest::~LaplacianTest() 
{
	std::cout << "~LaplacianTest" << std::endl;
}

void LaplacianTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void LaplacianTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}



INSTANTIATE_TEST_CASE_P(Solver, LaplacianTest, ::testing::Values(0) );


TEST_P(LaplacianTest, nNeib) {

	EXPECT_EQ(3, 5);

}

