#include "D2Q4StencilTest.h"

D2Q4StencilTest::D2Q4StencilTest()
{
	std::cout << "D2Q4StencilTest" << std::endl;

	int nx = 6;
	int ny = 6;
	double xMin = 0;
	double xMax = 5;
	double yMin = 0;
	double yMax = 5;

	field = Field(nx, ny, xMin, xMax, yMin, yMax);

	D2Q4 = D2Q4Stencil(field);

	double **coeff2D; 
	coeff2D = new double* [3];
	for (int i = 0; i < 3; i++)
		coeff2D[i] = new double[3];


	coeff2D[1][1] = 0;
	coeff2D[1][2] = 1;
	coeff2D[2][1] = 2;
	coeff2D[1][0] = 3;
	coeff2D[0][1] = 4;

	D2Q4.set_coeff2D(coeff2D);
}

D2Q4StencilTest::~D2Q4StencilTest() 
{
	std::cout << "~D2Q4StencilTest" << std::endl;
}

void D2Q4StencilTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void D2Q4StencilTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}



TEST_F(D2Q4StencilTest, nNeib) {

	EXPECT_EQ(D2Q4.nNeib, 5);

}

TEST_F(D2Q4StencilTest, iNeib_botLeftCorner) {

	int* iNeib = new int[5];

	D2Q4.get_iNeib(0, iNeib);
	EXPECT_EQ(iNeib[0], 0);
	EXPECT_EQ(iNeib[1], 6);
	EXPECT_EQ(iNeib[2], 1);
	EXPECT_EQ(iNeib[3], 30);
	EXPECT_EQ(iNeib[4], 5);

}

TEST_F(D2Q4StencilTest, iNeib_botRigtCorner) {

	int* iNeib = new int[5];

	D2Q4.get_iNeib(5, iNeib);
	EXPECT_EQ(iNeib[0], 5);
	EXPECT_EQ(iNeib[1], 11);
	EXPECT_EQ(iNeib[2], 0);
	EXPECT_EQ(iNeib[3], 35);
	EXPECT_EQ(iNeib[4], 4);

}

TEST_F(D2Q4StencilTest, iNeib_topRigtCorner) {

	int* iNeib = new int[5];

	D2Q4.get_iNeib(35, iNeib);
	EXPECT_EQ(iNeib[0], 35);
	EXPECT_EQ(iNeib[1], 5);
	EXPECT_EQ(iNeib[2], 30);
	EXPECT_EQ(iNeib[3], 29);
	EXPECT_EQ(iNeib[4], 34);

}

TEST_F(D2Q4StencilTest, iNeib_topLeftCorner) {

	int* iNeib = new int[5];

	D2Q4.get_iNeib(30, iNeib);
	EXPECT_EQ(iNeib[0], 30);
	EXPECT_EQ(iNeib[1], 0);
	EXPECT_EQ(iNeib[2], 31);
	EXPECT_EQ(iNeib[3], 24);
	EXPECT_EQ(iNeib[4], 35);

}

TEST_F(D2Q4StencilTest, iNeib_interior) {

	int* iNeib = new int[5];

	D2Q4.get_iNeib(8, iNeib);
	EXPECT_EQ(iNeib[0], 8);
	EXPECT_EQ(iNeib[1], 14);
	EXPECT_EQ(iNeib[2], 9);
	EXPECT_EQ(iNeib[3], 2);
	EXPECT_EQ(iNeib[4], 7);

}

TEST_F(D2Q4StencilTest, coeff1D) {

	EXPECT_EQ(D2Q4.coeff1D[0], 0.0);
	EXPECT_EQ(D2Q4.coeff1D[1], 1.0);
	EXPECT_EQ(D2Q4.coeff1D[2], 2.0);
	EXPECT_EQ(D2Q4.coeff1D[3], 3.0);
	EXPECT_EQ(D2Q4.coeff1D[4], 4.0);

}












