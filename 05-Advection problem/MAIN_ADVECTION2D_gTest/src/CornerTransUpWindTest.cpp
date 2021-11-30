#include "CornerTransUpWindTest.h"

CornerTransUpWindTest::CornerTransUpWindTest()
{
	std::cout << "CornerTransUpWindTest" << std::endl;

	int nx = 21;
	int ny = 21;
	double xMin = -0.5;
	double xMax = 0.5;
	double yMin = -0.5;
	double yMax = 0.5;

	field = Field(nx, ny, xMin, xMax, yMin, yMax);

	solver = CornerTransUpWind(field, 0.5, -0.3);

	Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
}

CornerTransUpWindTest::~CornerTransUpWindTest() 
{
	std::cout << "~CornerTransUpWindTest" << std::endl;
}

void CornerTransUpWindTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void CornerTransUpWindTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}



INSTANTIATE_TEST_CASE_P(Solver, CornerTransUpWindTest, ::testing::Values(0) );


TEST_P(CornerTransUpWindTest, timeStep) {

	double initVal = 1.0;

	field.init(1.0);

	solver.solve(0.9, 1, 2.0);

	EXPECT_NEAR(solver.dt, 0.09, 10e-6);

}


TEST_P(CornerTransUpWindTest, coeff2D) {

	double initVal = 1.0;

	field.init(1.0);

	solver.solve(0.9, 1, 2.0);

	EXPECT_NEAR(solver.coeff2D[1][1], 0.046, 10e-6);
	EXPECT_NEAR(solver.coeff2D[1][2], 0.054, 10e-6);
	EXPECT_NEAR(solver.coeff2D[2][2], 0.000, 10e-6);
	EXPECT_NEAR(solver.coeff2D[2][1], 0.000, 10e-6);
	EXPECT_NEAR(solver.coeff2D[2][0], 0.000, 10e-6);
	EXPECT_NEAR(solver.coeff2D[1][0], 0.000, 10e-6);
	EXPECT_NEAR(solver.coeff2D[0][0], 0.000, 10e-6);
	EXPECT_NEAR(solver.coeff2D[0][1], 0.414, 10e-6);
	EXPECT_NEAR(solver.coeff2D[0][2], 0.486, 10e-6);

}


TEST_P(CornerTransUpWindTest, fixedValueInit) {

	double initVal = 1.0;

	field.init(1.0);

	solver.solve(0.9, 1, 2.0);

	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
		EXPECT_EQ(field.var[iNode], initVal);
	}

}


TEST_P(CornerTransUpWindTest, bumpInit_firstIter) {

	field.init(&bump);

	double val66 = field.var[66];
	double val87 = field.var[87];
	double val67 = field.var[67];
	double val45 = field.var[45];
	double val65 = field.var[65];

	solver.solve(0.9, 1, 2.0);

	double value = solver.stencil.coeff1D[0] * val66 + 
		           solver.stencil.coeff1D[1] * val87 + 
		           solver.stencil.coeff1D[2] * val67 + 
		           solver.stencil.coeff1D[3] * val45 +
				   solver.stencil.coeff1D[5] * val65;

	EXPECT_EQ(field.var[66], value);
}


TEST_P(CornerTransUpWindTest, bumpInit_firstTimeStep) {

	field.init(&bump);

	solver.solve(0.9, 1000, 9 / 160.0);

	int iNode = 66;
	EXPECT_NEAR(field.var[iNode], 0.000186806, 10e-6);

}
