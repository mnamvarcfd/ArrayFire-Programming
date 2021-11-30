#include "LaxWendroffDimSplitTest.h"

LaxWendroffDimSplitTest::LaxWendroffDimSplitTest() 
{
	std::cout << "LaxWendroffDimSplitTest" << std::endl;

	int nx = 21;
	int ny = 21;
	double xMin = -0.5;
	double xMax = 0.5;
	double yMin = -0.5;
	double yMax = 0.5;

	field = Field(nx, ny, xMin, xMax, yMin, yMax);


	Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
}

LaxWendroffDimSplitTest::~LaxWendroffDimSplitTest() 
{
	std::cout << "~LaxWendroffDimSplitTest" << std::endl;
}

void LaxWendroffDimSplitTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void LaxWendroffDimSplitTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}



INSTANTIATE_TEST_CASE_P(Solver, LaxWendroffDimSplitTest, ::testing::Values(0) );


TEST_P(LaxWendroffDimSplitTest, timeStep) {

	double initVal = 1.0;

	field.init(1.0);

	solver = LaxWendroffDimSplit(field, 0.5, -0.3);
	solver.solve(0.9, 1, 2.0);

	EXPECT_NEAR(solver.dt, 0.09, 10e-6);

}


TEST_P(LaxWendroffDimSplitTest, coeff2D) {

	double initVal = 1.0;

	field.init(1.0);

	solver = LaxWendroffDimSplit(field, 0.5, -0.3);
	solver.solve(0.9, 1, 2.0);

	EXPECT_NEAR(solver.laxWendroffDimX.coeff2D[1][1], 0.19, 10e-6);
	EXPECT_NEAR(solver.laxWendroffDimX.coeff2D[2][1], -0.045, 10e-6);
	EXPECT_NEAR(solver.laxWendroffDimX.coeff2D[0][1], 0.855, 10e-6);

	EXPECT_NEAR(solver.laxWendroffDimY.coeff2D[1][1], 0.7084, 10e-6);
	EXPECT_NEAR(solver.laxWendroffDimY.coeff2D[1][2], 0.4158, 10e-6);
	EXPECT_NEAR(solver.laxWendroffDimY.coeff2D[1][0], -0.1242, 10e-6);
}


TEST_P(LaxWendroffDimSplitTest, fixedValueInit) {

	double initVal = 1.0;

	field.init(1.0);

	solver = LaxWendroffDimSplit(field, 0.5, -0.3);
	solver.solve(0.9, 1, 2.0);


	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
		EXPECT_EQ(field.var[iNode], initVal);
	}

}


TEST_P(LaxWendroffDimSplitTest, bumpInit_firstIter) {

	field.init(&bump);

	double val44 = field.var[44];
	double val45 = field.var[45];
	double val46 = field.var[46];
	double val86 = field.var[86];
	double val65 = field.var[65];
	double val66 = field.var[66];
	double val87 = field.var[87];
	double val88 = field.var[88];
	double val67 = field.var[67];


	solver = LaxWendroffDimSplit(field, 0.5, -0.3);
	solver.solve(0.9, 1, 2.0);


	double value45star = solver.laxWendroffDimX.stencil.coeff1D[0] * val45 +
		solver.laxWendroffDimX.stencil.coeff1D[1] * val46 +
		solver.laxWendroffDimX.stencil.coeff1D[2] * val44;

	double value66star = solver.laxWendroffDimX.stencil.coeff1D[0] * val66 +
		solver.laxWendroffDimX.stencil.coeff1D[1] * val67 +
		solver.laxWendroffDimX.stencil.coeff1D[2] * val65;


	double value87star = solver.laxWendroffDimX.stencil.coeff1D[0] * val87 +
		solver.laxWendroffDimX.stencil.coeff1D[1] * val88 +
		solver.laxWendroffDimX.stencil.coeff1D[2] * val86;



	double value = solver.laxWendroffDimY.stencil.coeff1D[0] * value66star +
		solver.laxWendroffDimY.stencil.coeff1D[1] * value87star +
		solver.laxWendroffDimY.stencil.coeff1D[2] * value45star;

	EXPECT_EQ(field.var[66], value);
}

