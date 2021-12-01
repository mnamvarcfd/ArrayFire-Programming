#include "BumpValidation.h"

double bumpValid(double x, double y) {

	double val = 0.0;

	if ((x * x + y * y) <0.25 ) val = exp(1.0 - 0.25 / (0.25 - x * x - y * y));

	return val;
}

BumpValidation::BumpValidation()
{
	std::cout << "BumpValidation" << std::endl;	


	a = 0.5;
	b = -0.3;
	//a = 1.0;
	//b = 0.0;

	double n = 1.0;
	double lambda = min(abs(1.0 / cos(atan2(b, a))), abs(1.0 / cos(atan2(b, a))));
	tFinal = n * lambda / sqrt(a * a + b * b);
	std::cout << "tFinal: " << tFinal << std::endl;


}

BumpValidation::~BumpValidation() 
{
	std::cout << "~BumpValidation" << std::endl;
}

void BumpValidation::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void BumpValidation::TearDown()
{
	std::cout << "TearDown" << std::endl;
	FAIL();
}


INSTANTIATE_TEST_CASE_P(Verification, BumpValidation, ::testing::Values(0) );


int bumpGridSize[] = { 21, 41, 81, 161, 321, 641, 1281 };

TEST_P(BumpValidation, laxWendroff) {

	FILE* file;
	fopen_s(&file, "ErrorBump_LaxWendroff.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(bumpGridSize) / 4; i++) {

		int nx = bumpGridSize[i];
		int ny = nx;
		double dx = 0.5 * (1.0 / (nx));
		double xMin = -0.5 + dx;
		double xMax = 0.5 - dx;
		double yMin = -0.5 + dx;
		double yMax = 0.5 - dx;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&bumpValid);

		LaxWendroff solver = LaxWendroff(U, a, b);


		solver.solve(0.9, 10e6, tFinal);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&bumpValid, a, b, tFinal);


		//Calculating the L1 error
		double error = 0.0;
		for (int j = 0; j < U.get_nNode(); j++) {
			error += abs(U.var[j] - Ua.var[j]);
		}
		error = error / U.get_nNode();

		printf("%.8f  \n", error);
		fprintf(file, "%d   %.8f  \n ", nx, error);

	}

	fclose(file);
}

TEST_P(BumpValidation, laxWendroffDimSplit) {

	FILE* file;
	fopen_s(&file, "ErrorBump_LaxWendroffDimSplit.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(bumpGridSize) / 4; i++) {

		int nx = bumpGridSize[i];
		int ny = nx;
		double dx = 0.5 * (1.0 / (nx));
		double xMin = -0.5 + dx;
		double xMax = 0.5 - dx;
		double yMin = -0.5 + dx;
		double yMax = 0.5 - dx;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&bumpValid);


		LaxWendroffDimSplit solver = LaxWendroffDimSplit(U, a, b);


		solver.solve(0.9, 10e6, tFinal);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&bumpValid, a, b, tFinal);


		//Calculating the L1 error
		double error = 0.0;
		for (int j = 0; j < U.get_nNode(); j++) {
			error += abs(U.var[j] - Ua.var[j]);
		}
		error = error / U.get_nNode();

		printf("%.8f  \n", error);
		fprintf(file, "%d   %e  \n ", nx, error);

	}

	fclose(file);
}

TEST_P(BumpValidation, donerCellUpWind) {

	FILE* file;
	fopen_s(&file, "ErrorBump_DonerCellUpWind.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(bumpGridSize) / 4; i++) {

		int nx = bumpGridSize[i];
		int ny = nx;
		double dx = 0.5 * (1.0 / (nx));
		double xMin = -0.5 + dx;
		double xMax = 0.5 - dx;
		double yMin = -0.5 + dx;
		double yMax = 0.5 - dx;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&bumpValid);


		DonerCellUpWind solver = DonerCellUpWind(U, a, b);


		solver.solve(0.9, 10e6, 2.0);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&bumpValid, 0.5, -0.3, 2.0);


		//Calculating the L2 error
		double error = 0.0;
		for (int j = 0; j < U.get_nNode(); j++) {
			error += abs(U.var[j] - Ua.var[j]);
		}
		error = error / U.get_nNode();

		printf("%.8f  \n", error);
		fprintf(file, "%d   %e  \n ", nx, error);

	}

	fclose(file);
}

TEST_P(BumpValidation, cornerTransUpWind) {

	FILE* file;
	fopen_s(&file, "ErrorBump_CornerTransUpWind.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(bumpGridSize) / 4; i++) {

		int nx = bumpGridSize[i];
		int ny = nx;
		double dx = 0.5 * (1.0 / (nx));
		double xMin = -0.5 + dx;
		double xMax = 0.5 - dx;
		double yMin = -0.5 + dx;
		double yMax = 0.5 - dx;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&bumpValid);


		CornerTransUpWind solver = CornerTransUpWind(U, a, b);


		solver.solve(0.9, 10e6, tFinal);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&bumpValid, a, b, tFinal);


		//Calculating the L1 error
		double error = 0.0;
		for (int j = 0; j < U.get_nNode(); j++) {
			error += abs(U.var[j] - Ua.var[j]);
		}
		error = error / U.get_nNode();

		printf("%.8f  \n", error);
		fprintf(file, "%d   %e  \n ", nx, error);

	}

	fclose(file);
}









