#include "D2Q9StencilTest.h"

D2Q9StencilTest::D2Q9StencilTest()
{
	std::cout << "D2Q9StencilTest" << std::endl;

	int nx = 7;
	int ny = 6;
	double xMin = 0;
	double xMax = 5;
	double yMin = 0;
	double yMax = 5;

	field = Field(nx, ny, xMin, xMax, yMin, yMax);

	stencil = D2Q9Stencil(field);

	double **coeff2D; 
	coeff2D = new double* [3];
	for (int i = 0; i < 3; i++)
		coeff2D[i] = new double[3];


	coeff2D[1][1] = 0;
	coeff2D[1][2] = 1;
	coeff2D[2][2] = 2;
	coeff2D[2][1] = 3;
	coeff2D[2][0] = 4;
	coeff2D[1][0] = 5;
	coeff2D[0][0] = 6;
	coeff2D[0][1] = 7;
	coeff2D[0][2] = 8;

	stencil.set_coeff2D(coeff2D);
}

D2Q9StencilTest::~D2Q9StencilTest() 
{
	std::cout << "~D2Q9StencilTest" << std::endl;
}

void D2Q9StencilTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void D2Q9StencilTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}


TEST_F(D2Q9StencilTest, coeff1D) {

	EXPECT_EQ(stencil.coeff1D[0], 0.0);
	EXPECT_EQ(stencil.coeff1D[1], 1.0);
	EXPECT_EQ(stencil.coeff1D[2], 2.0);
	EXPECT_EQ(stencil.coeff1D[3], 3.0);
	EXPECT_EQ(stencil.coeff1D[4], 4.0);
	EXPECT_EQ(stencil.coeff1D[5], 5.0);
	EXPECT_EQ(stencil.coeff1D[6], 6.0);
	EXPECT_EQ(stencil.coeff1D[7], 7.0);
	EXPECT_EQ(stencil.coeff1D[8], 8.0);

}

TEST_F(D2Q9StencilTest, nNeib) {

	EXPECT_EQ(stencil.nNeib, 9);

}

TEST_F(D2Q9StencilTest, iNeib_interior) {

	int* iNeib = new int[9];

	stencil.get_iNeib(16, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_bottomSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(1, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_TopSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(36, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_leftSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(7, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_rightSide) {

	int* iNeib = new int[9];

	stencil.get_iNeib(13, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_botLeftCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(0, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_botRigtCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(6, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_topRigtCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(41, iNeib);
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

TEST_F(D2Q9StencilTest, iNeib_topLeftCorner) {

	int* iNeib = new int[9];

	stencil.get_iNeib(35, iNeib);
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
