#include "Examples.h"


double bumpExample(double x, double y) {

	double val = 0.0;

	if ((x * x + y * y) < 0.25) val = exp(1.0 - 0.25 / (0.25 - x * x - y * y));

	return val;
}

static double SquarWaveExample(double x, double y) {

	double val = 0.0;

	if (abs(x) <= 0.125 && abs(y) <= 0.125) val = 1.0;

	return val;
}

Examples::Examples()
{
	std::cout << "Examples" << std::endl;	

	a = 1.0;
	b = 0.0;

	double n = 1.0;
	double lambda = min(abs(1.0 / cos(atan2(b, a))), abs(1.0 / cos(atan2(b, a))));
	tFinal = n * lambda / sqrt(a * a + b * b);
	std::cout << "tFinal: " << tFinal << std::endl;


}

Examples::~Examples() 
{
	std::cout << "~Examples" << std::endl;
}

void Examples::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void Examples::TearDown()
{
	std::cout << "TearDown" << std::endl;
	FAIL();
}


INSTANTIATE_TEST_CASE_P(Example, Examples, ::testing::Values(0) );


TEST_P(Examples, forReport) {

	int nx = 82;
	int ny = nx;
	double dx = 0.5 * (1.0 / (nx));
	double xMin = -0.5 + dx;
	double xMax = 0.5 - dx;
	double yMin = -0.5 + dx;
	double yMax = 0.5 - dx;


	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init(&bumpExample);
	U.write("init");


	//LaxWendroff solver = LaxWendroff(U, a, b);
	DonerCellUpWind solver = DonerCellUpWind(U, a, b);
	//CornerTransUpWind solver = CornerTransUpWind(U, a, b);
	//LaxWendroffDimSplit solver = LaxWendroffDimSplit(U, a, b);

	solver.solve(0.9, 10e6, tFinal);
	//solver.solveParallel(0.9, tFinal);
	U.write("numericalSolution");


	Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ua.init(&bumpExample, a, b, tFinal);
	Ua.write("analyticalSolution");

	//Calculating the L1 error
	double error = 0.0;
	for (int j = 0; j < U.get_nNode(); j++) {
		error += abs(U.var[j] - Ua.var[j]);
	}
	error = error / U.get_nNode();

	printf("L1 error: %.8f  \n", error);
}


TEST_P(Examples, Parallel) {

	af::setBackend(AF_BACKEND_OPENCL);  //   AF_BACKEND_CPU
	af::setDevice(0);
	af::info();

	int GSize[] = { 10e4, 10e5, 10e6, 10e7 };

	for (int i = 0; i < 4; i++) {

		int nx = sqrt((double)GSize[i]);
		int ny = nx;
		double dx = 0.5 * (1.0 / (nx));
		double xMin = -0.5 + dx;
		double xMax = 0.5 - dx;
		double yMin = -0.5 + dx;
		double yMax = 0.5 - dx;


		Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&bumpExample);


		DonerCellUpWind solver = DonerCellUpWind(U, a, b);


		af::timer t1 = af::timer::start();
		solver.solveParallel(0.9, 2.0);
		printf("Overal time parallel: %g \n", af::timer::stop(t1));


	}

}


TEST_P(Examples, NaivCpp) {


	int GSize[] = { 10e4, 10e5, 10e6, 10e7, 10e8 };

	for (int i = 0; i < 5; i++) {

		int nx = sqrt( (double)GSize[i]) ;
		int ny = nx;
		double dx = 0.5 * (1.0 / (nx));
		double xMin = -0.5 + dx;
		double xMax = 0.5 - dx;
		double yMin = -0.5 + dx;
		double yMax = 0.5 - dx;


		Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&bumpExample);


		DonerCellUpWind solver = DonerCellUpWind(U, a, b);


		af::timer t1 = af::timer::start();
		solver.solve(0.9, 100, 2.0);
		printf("Overal time Serial: %g \n", af::timer::stop(t1));


	}

}




