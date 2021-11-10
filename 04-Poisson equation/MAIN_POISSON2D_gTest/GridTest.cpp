#include "GridTest.h"
#include "../MAIN_POISSON2D/DenseSolver.h"

GridTest::GridTest()
{
	std::cout << "GridTest" << std::endl;

    tol = 10e-6;

	pi = 4.0 * atan(1.0);
	nx = 7;
	ny = 6;
	xMin = -3.0;
	xMax = 3.0 * pi;
	yMin = 3.0;
	yMax = 4.0 * pi;


	grid = Grid(nx, ny, xMin, xMax, yMin, yMax);

}

GridTest::~GridTest() 
{
	std::cout << "~GridTest" << std::endl;
}

void GridTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void GridTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}



INSTANTIATE_TEST_CASE_P(2Grid, GridTest, ::testing::Values(0) );


TEST_P(GridTest, dx) {

	double expected = (xMax - xMin) / ((double)nx - 1);

	EXPECT_NEAR(grid.get_dx(), expected, tol);

}

TEST_P(GridTest, dy) {

	double expected = (yMax - yMin) / ((double)ny - 1);

	EXPECT_NEAR(grid.get_dy(), expected, tol);

}

TEST_P(GridTest, getXInx) {

	EXPECT_EQ(grid.getXInx(6), 6);

}

TEST_P(GridTest, getXvalue) {

	EXPECT_NEAR(grid.get_xVal(0), xMin, tol);
	EXPECT_NEAR(grid.get_xVal(6), xMax, tol);

}











