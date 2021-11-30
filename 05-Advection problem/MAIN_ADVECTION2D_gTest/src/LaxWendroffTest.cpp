#include "LaxWendroffTest.h"
//#include "../../MAIN_ADVECTION2D/src/LaxWendroff.cpp"


LaxWendroffTest::LaxWendroffTest()
{
	std::cout << "LaxWendroffTest" << std::endl;

	int nx = 21;
	int ny = 21;
	double xMin = -0.5;
	double xMax = 0.5;
	double yMin = -0.5;
	double yMax = 0.5;

	field = Field(nx, ny, xMin, xMax, yMin, yMax);

	solver = LaxWendroff(field, 0.5, -0.3);

	Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
}

LaxWendroffTest::~LaxWendroffTest() 
{
	std::cout << "~LaxWendroffTest" << std::endl;
}

void LaxWendroffTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void LaxWendroffTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}



INSTANTIATE_TEST_CASE_P(Solver, LaxWendroffTest, ::testing::Values(0) );


TEST_P(LaxWendroffTest, timeStep) {

	double initVal = 1.0;

	field.init(1.0);

	solver.solve(0.9, 1, 2.0);

	EXPECT_NEAR(solver.dt, 0.0545705, 10e-6);

}


TEST_P(LaxWendroffTest, coeff2D) {

	double initVal = 1.0;

	field.init(1.0);

	solver.solve(0.9, 1, 2.0);


	EXPECT_NEAR(solver.coeff2D[1][1], 0.595000, 10e-6);

}


TEST_P(LaxWendroffTest, fixedValueInit) {

	double initVal = 1.0;

	field.init(1.0);

	solver.solve(0.9, 1, 2.0);


	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
		EXPECT_EQ(field.var[iNode], initVal);
	}

}


TEST_P(LaxWendroffTest, bumpInit_firstIter) {

	field.init(&bump);

	double val66 = field.var[66];
	double val87 = field.var[87];
	double val88 = field.var[88];
	double val67 = field.var[67];

	solver.solve(0.9, 1, 2.0);


	double value = solver.stencil.coeff1D[0] * val66 + 
		           solver.stencil.coeff1D[1] * val87 + 
		           solver.stencil.coeff1D[2] * val88 + 
		           solver.stencil.coeff1D[3] * val67;

	EXPECT_EQ(field.var[66], value);
}


TEST_P(LaxWendroffTest, bumpInit_firstTimeStep) {

	field.init(&bump);

	solver.solve(0.9, 1000, 0.5);

	int iNode = 66;
	EXPECT_NEAR(field.var[iNode], 0.030641, 10e-6);

}


