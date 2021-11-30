#include "GridTest.h"

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
	grid.writGrid("gridTest");

	double* index = new double[42];
	for (int iNode = 0; iNode < 42; iNode++) {
		index[iNode] = iNode;
	}
	grid.writePLT("grid", index);
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

	EXPECT_NEAR(grid.get_dx(), 2.070796, tol);

}

TEST_P(GridTest, dy) {

	double expected = (yMax - yMin) / ((double)ny - 1);

	EXPECT_NEAR(grid.get_dy(), expected, tol);

}

TEST_P(GridTest, getXInx) {

	EXPECT_EQ(grid.getXInx(6), 6);

}

//TEST_P(GridTest, getXInx_double) {
//
//	EXPECT_EQ(grid.getXInx(5.28319), 4);
//
//}

//TEST_P(GridTest, getXInx_double_all) {
//
//	for (int i = 0; i < 42; i++) {
//		double x = grid.get_xVal(i);
//		double y = grid.get_yVal(i);
//
//		EXPECT_EQ(grid.getXInx(x), grid.getXInx(i));
//		EXPECT_EQ(grid.getYInx(y), grid.getYInx(i));
//	}
//
//}


TEST_P(GridTest, getXvalue) {

	EXPECT_NEAR(grid.get_xVal(0), xMin, tol);
	EXPECT_NEAR(grid.get_xVal(6), xMax, tol);

}

TEST_P(GridTest, iNeib_interior) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(16, iNeib);
	EXPECT_EQ(iNeib[0], 16);
	EXPECT_EQ(iNeib[1], 23);
	EXPECT_EQ(iNeib[2], 24);
	EXPECT_EQ(iNeib[3], 17);
	EXPECT_EQ(iNeib[4], 10);
	EXPECT_EQ(iNeib[5], 9);
	EXPECT_EQ(iNeib[6], 8);
	EXPECT_EQ(iNeib[7], 15);
	EXPECT_EQ(iNeib[8], 22);

}

TEST_P(GridTest, iNeib_bottomSide) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(1, iNeib);
	EXPECT_EQ(iNeib[0], 1);
	EXPECT_EQ(iNeib[1], 8);
	EXPECT_EQ(iNeib[2], 9);
	EXPECT_EQ(iNeib[3], 2);
	EXPECT_EQ(iNeib[4], 37);
	EXPECT_EQ(iNeib[5], 36);
	EXPECT_EQ(iNeib[6], 35);
	EXPECT_EQ(iNeib[7], 0);
	EXPECT_EQ(iNeib[8], 7);

}

TEST_P(GridTest, iNeib_TopSide) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(36, iNeib);
	EXPECT_EQ(iNeib[0], 36);
	EXPECT_EQ(iNeib[1], 1);
	EXPECT_EQ(iNeib[2], 2);
	EXPECT_EQ(iNeib[3], 37);
	EXPECT_EQ(iNeib[4], 30);
	EXPECT_EQ(iNeib[5], 29);
	EXPECT_EQ(iNeib[6], 28);
	EXPECT_EQ(iNeib[7], 35);
	EXPECT_EQ(iNeib[8], 0);

}

TEST_P(GridTest, iNeib_leftSide) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(7, iNeib);
	EXPECT_EQ(iNeib[0], 7);
	EXPECT_EQ(iNeib[1], 14);
	EXPECT_EQ(iNeib[2], 15);
	EXPECT_EQ(iNeib[3], 8);
	EXPECT_EQ(iNeib[4], 1);
	EXPECT_EQ(iNeib[5], 0);
	EXPECT_EQ(iNeib[6], 6);
	EXPECT_EQ(iNeib[7], 13);
	EXPECT_EQ(iNeib[8], 20);

}

TEST_P(GridTest, iNeib_rightSide) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(13, iNeib);
	EXPECT_EQ(iNeib[0], 13);
	EXPECT_EQ(iNeib[1], 20);
	EXPECT_EQ(iNeib[2], 14);
	EXPECT_EQ(iNeib[3], 7);
	EXPECT_EQ(iNeib[4], 0);
	EXPECT_EQ(iNeib[5], 6);
	EXPECT_EQ(iNeib[6], 5);
	EXPECT_EQ(iNeib[7], 12);
	EXPECT_EQ(iNeib[8], 19);

}

TEST_P(GridTest, iNeib_botLeftCorner) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(0, iNeib);
	EXPECT_EQ(iNeib[0], 0);
	EXPECT_EQ(iNeib[1], 7);
	EXPECT_EQ(iNeib[2], 8);
	EXPECT_EQ(iNeib[3], 1);
	EXPECT_EQ(iNeib[4], 36);
	EXPECT_EQ(iNeib[5], 35);
	EXPECT_EQ(iNeib[6], 33);
	EXPECT_EQ(iNeib[7], 6);
	EXPECT_EQ(iNeib[8], 13);

}

TEST_P(GridTest, iNeib_botRigtCorner) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(6, iNeib);
	EXPECT_EQ(iNeib[0], 6);
	EXPECT_EQ(iNeib[1], 13);
	EXPECT_EQ(iNeib[2], 7);
	EXPECT_EQ(iNeib[3], 0);
	EXPECT_EQ(iNeib[4], 29);
	EXPECT_EQ(iNeib[5], 41);
	EXPECT_EQ(iNeib[6], 40);
	EXPECT_EQ(iNeib[7], 5);
	EXPECT_EQ(iNeib[8], 12);

}

TEST_P(GridTest, iNeib_topRigtCorner) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(41, iNeib);
	EXPECT_EQ(iNeib[0], 41);
	EXPECT_EQ(iNeib[1], 6);
	EXPECT_EQ(iNeib[2], 8);
	EXPECT_EQ(iNeib[3], 35);
	EXPECT_EQ(iNeib[4], 28);
	EXPECT_EQ(iNeib[5], 34);
	EXPECT_EQ(iNeib[6], 33);
	EXPECT_EQ(iNeib[7], 40);
	EXPECT_EQ(iNeib[8], 5);

}

TEST_P(GridTest, iNeib_topLeftCorner) {

	int* iNeib = new int[9];

	grid.get_iNeibPeriodicD2Q9(35, iNeib);
	EXPECT_EQ(iNeib[0], 35);
	EXPECT_EQ(iNeib[1], 0);
	EXPECT_EQ(iNeib[2], 1);
	EXPECT_EQ(iNeib[3], 36);
	EXPECT_EQ(iNeib[4], 29);
	EXPECT_EQ(iNeib[5], 28);
	EXPECT_EQ(iNeib[6], 34);
	EXPECT_EQ(iNeib[7], 41);
	EXPECT_EQ(iNeib[8], 12);

}
















