#include "ParallelValidation.h"

double bumpParallel(double x, double y) {

	double val = 0.0;

	if ((x * x + y * y) <0.25 ) val = exp(1.0 - 0.25 / (0.25 - x * x - y * y));

	return val;
}

ParallelValidation::ParallelValidation()
{
	std::cout << "ParallelValidation" << std::endl;	

	af::setBackend(AF_BACKEND_CPU);  //AF_BACKEND_OPENCL   
	af::setDevice(0);
	af::info();


	a = 1.0;
	b = 0.0;

	double n = 1.0;
	double lambda = min(abs(1.0 / cos(atan2(b, a))), abs(1.0 / cos(atan2(b, a))));
	tFinal = n * lambda / sqrt(a * a + b * b);
	std::cout << "tFinal: " << tFinal << std::endl;


	int nx = 21;
	int ny = nx;
	double dx = 0.5 * (1.0 / (nx));
	double xMin = -0.5 + dx;
	double xMax = 0.5 - dx;
	double yMin = -0.5 + dx;
	double yMax = 0.5 - dx;

	User = Field(nx, ny, xMin, xMax, yMin, yMax);
	Upar = Field(nx, ny, xMin, xMax, yMin, yMax);

	User.init(&bumpParallel);
	Upar.init(&bumpParallel);

}

ParallelValidation::~ParallelValidation() 
{
	std::cout << "~ParallelValidation" << std::endl;
}

void ParallelValidation::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void ParallelValidation::TearDown()
{
	std::cout << "TearDown" << std::endl;
	//FAIL();
}


INSTANTIATE_TEST_CASE_P(Validation, ParallelValidation, ::testing::Values(0) );


TEST_P(ParallelValidation, laxWendroff) {


	LaxWendroff solverSer = LaxWendroff(User, a, b);
	solverSer.solve(0.9, 10e6, tFinal);


	LaxWendroff solverPar = LaxWendroff(Upar, a, b);
	solverPar.solveParallel(0.9, tFinal);


	for (int j = 0; j < User.get_nNode(); j++) {
		EXPECT_NEAR(User.var[j] , Upar.var[j], eps);
	}
	
}

TEST_P(ParallelValidation, laxWendroffDimSplit) {

	LaxWendroffDimSplit solverSer = LaxWendroffDimSplit(User, a, b);
	solverSer.solve(0.9, 10e6, tFinal);


	LaxWendroffDimSplit solverPar = LaxWendroffDimSplit(Upar, a, b);
	solverPar.solveParallel(0.9, tFinal);


	for (int j = 0; j < User.get_nNode(); j++) {
		EXPECT_NEAR(User.var[j], Upar.var[j], eps);
	}

}

TEST_P(ParallelValidation, donerCellUpWind) {

	LaxWendroffDimSplit solverSer = LaxWendroffDimSplit(User, a, b);
	solverSer.solve(0.9, 10e6, tFinal);


	LaxWendroffDimSplit solverPar = LaxWendroffDimSplit(Upar, a, b);
	solverPar.solveParallel(0.9, tFinal);


	for (int j = 0; j < User.get_nNode(); j++) {
		EXPECT_NEAR(User.var[j], Upar.var[j], eps);
	}

}

TEST_P(ParallelValidation, cornerTransUpWind) {

	CornerTransUpWind solverSer = CornerTransUpWind(User, a, b);
	solverSer.solve(0.9, 10e6, tFinal);


	CornerTransUpWind solverPar = CornerTransUpWind(Upar, a, b);
	solverPar.solveParallel(0.9, tFinal);


	for (int j = 0; j < User.get_nNode(); j++) {
		EXPECT_NEAR(User.var[j], Upar.var[j], eps);
	}

}






TEST_P(ParallelValidation, Example) {

	//LaxWendroff solver = LaxWendroff(Upar, a, b);
	DonerCellUpWind solver = DonerCellUpWind(Upar, a, b);
	//CornerTransUpWind solver = CornerTransUpWind(Upar, a, b);
	//LaxWendroffDimSplit solver = LaxWendroffDimSplit(Upar, a, b);

	solver.solve(0.9, 1000000, tFinal);
	//solver.solveParallel(0.9, tFinal);
	Upar.write("numericalSolution");


	User.init(&bumpParallel, a, b, tFinal);
	User.write("analyticalSolution");

	//Calculating the L2 error
	double error = 0.0;
	for (int j = 0; j < Upar.get_nNode(); j++) {
		error += abs(Upar.var[j] - User.var[j]);
	}
	error = error / Upar.get_nNode();

	printf("%.16f  \n", error);

	FAIL();
}


TEST_P(ParallelValidation, Performance) {

	af::setBackend(AF_BACKEND_OPENCL);  //AF_BACKEND_CPU   
	af::setDevice(0);
	af::info();

	tFinal *= 1;

	int nx = 1000;
	int ny = nx;
	double dx = 0.5 * (1.0 / nx);
	double xMin = -0.5 + dx;
	double xMax = 0.5 - dx;
	double yMin = -0.5 + dx;
	double yMax = 0.5 - dx;

	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);

	U.init(&bumpParallel);
	DonerCellUpWind solver = DonerCellUpWind(U, a, b);

	af::timer t1 = af::timer::start();
	solver.solve(0.9, 10000, tFinal);
	double timeSer = af::timer::stop(t1);

	printf("Time serial %g \n", timeSer);


	U.init(&bumpParallel);

	t1 = af::timer::start();
	solver.solveParallel(0.9, tFinal);
	double timePar = af::timer::stop(t1);

	printf("Time parallel %g \n", timePar);


	printf("speed up %g = \n", timeSer/timePar);





	FAIL();
}








