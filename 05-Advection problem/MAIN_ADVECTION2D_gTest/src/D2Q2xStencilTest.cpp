#include "D2Q2xStencilTest.h"

D2Q2xStencilTest::D2Q2xStencilTest()
{
	std::cout << "D2Q2xStencilTest" << std::endl;

	int nx = 7;
	int ny = 6;
	double xMin = 0;
	double xMax = 5;
	double yMin = 0;
	double yMax = 5;

	field = Field(nx, ny, xMin, xMax, yMin, yMax);

	stencil = D2Q2xStencil(field);

	double **coeff2D; 
	coeff2D = new double* [3];
	for (int i = 0; i < 3; i++)
		coeff2D[i] = new double[3];


	coeff2D[1][1] = 0;
	coeff2D[2][1] = 1;
	coeff2D[0][1] = 2;

	stencil.set_coeff2D(coeff2D);
}

D2Q2xStencilTest::~D2Q2xStencilTest() 
{
	std::cout << "~D2Q2xStencilTest" << std::endl;
}

void D2Q2xStencilTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void D2Q2xStencilTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}


TEST_F(D2Q2xStencilTest, coeff1D) {

	EXPECT_EQ(stencil.coeff1D[0], 0.0);
	EXPECT_EQ(stencil.coeff1D[1], 1.0);
	EXPECT_EQ(stencil.coeff1D[2], 2.0);

}

TEST_F(D2Q2xStencilTest, nNeib) {

	EXPECT_EQ(stencil.nNeib, 3);

}

TEST_F(D2Q2xStencilTest, iNeib_interior) {

	int* iNeib = new int[9];

	stencil.get_iNeib(16, iNeib);
	EXPECT_EQ(iNeib[0], 16);
	EXPECT_EQ(iNeib[1], 17);
	EXPECT_EQ(iNeib[2], 15);

}

TEST_F(D2Q2xStencilTest, iNeib_bottomSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(1, iNeib);
	EXPECT_EQ(iNeib[0], 1);
	EXPECT_EQ(iNeib[1], 2);
	EXPECT_EQ(iNeib[2], 0);

}

TEST_F(D2Q2xStencilTest, iNeib_TopSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(36, iNeib);
	EXPECT_EQ(iNeib[0], 36);
	EXPECT_EQ(iNeib[1], 37);
	EXPECT_EQ(iNeib[2], 35);

}

TEST_F(D2Q2xStencilTest, iNeib_leftSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(7, iNeib);
	EXPECT_EQ(iNeib[0], 7);
	EXPECT_EQ(iNeib[1], 8);
	EXPECT_EQ(iNeib[2], 13);

}

TEST_F(D2Q2xStencilTest, iNeib_rightSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(13, iNeib);
	EXPECT_EQ(iNeib[0], 13);
	EXPECT_EQ(iNeib[1], 7);
	EXPECT_EQ(iNeib[2], 12);

}

TEST_F(D2Q2xStencilTest, iNeib_botLeftCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(0, iNeib);
	EXPECT_EQ(iNeib[0], 0);
	EXPECT_EQ(iNeib[1], 1);
	EXPECT_EQ(iNeib[2], 6);

}

TEST_F(D2Q2xStencilTest, iNeib_botRigtCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(6, iNeib);
	EXPECT_EQ(iNeib[0], 6);
	EXPECT_EQ(iNeib[1], 0);
	EXPECT_EQ(iNeib[2], 5);

}

TEST_F(D2Q2xStencilTest, iNeib_topRigtCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(41, iNeib);
	EXPECT_EQ(iNeib[0], 41);
	EXPECT_EQ(iNeib[1], 35);
	EXPECT_EQ(iNeib[2], 40);

}

TEST_F(D2Q2xStencilTest, iNeib_topLeftCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(35, iNeib);
	EXPECT_EQ(iNeib[0], 35);
	EXPECT_EQ(iNeib[1], 36);
	EXPECT_EQ(iNeib[2], 41);

}
