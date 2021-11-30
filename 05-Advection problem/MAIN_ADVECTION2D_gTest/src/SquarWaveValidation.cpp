#include "SquarWaveValidation.h"

SquarWaveValidation::SquarWaveValidation() 
{
	std::cout << "SquarWaveValidation" << std::endl;

	a = 1.0;
	b = 0.0;

	double n = 1.0;
	double lambda = min(abs(1.0 / cos(atan2(b, a))), abs(1.0 / cos(atan2(b, a))));
	tFinal = n * lambda / sqrt(a * a + b * b);
	std::cout << "tFinal: " << tFinal << std::endl;

}

SquarWaveValidation::~SquarWaveValidation() 
{
	std::cout << "~SquarWaveValidation" << std::endl;
}

void SquarWaveValidation::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void SquarWaveValidation::TearDown()
{
	std::cout << "TearDown" << std::endl;
	FAIL();
}



static double SquarWave(double x, double y) {

	double val = 0.0;

	if (abs(x) <= 0.125 && abs(y) <= 0.125) val = 1.0;

	return val;
}


INSTANTIATE_TEST_CASE_P(Verification, SquarWaveValidation, ::testing::Values(0) );


int SquarWaveGridSize[] = { 24, 48, 96, 192 /*, 348, 768, 1536*/ };

TEST_P(SquarWaveValidation, laxWendroff) {

	FILE* file;
	fopen_s(&file, "ErrorSuare_LaxWendroff.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(SquarWaveGridSize) / 4; i++) {

		int nx = SquarWaveGridSize[i];
		int ny = SquarWaveGridSize[i];
		double xMin = -0.5;
		double xMax = 0.5;
		double yMin = -0.5;
		double yMax = 0.5;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&SquarWave);

		LaxWendroff solver = LaxWendroff(U, 0.5, -0.5);


		solver.solve(0.9, 10e6, 2.0);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&SquarWave, 0.5, -0.5, 2.0);


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

TEST_P(SquarWaveValidation, laxWendroffDimSplit) {

	FILE* file;
	fopen_s(&file, "ErrorSuare_LaxWendroffDimSplit.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(SquarWaveGridSize) / 4; i++) {

		int nx = SquarWaveGridSize[i];
		int ny = SquarWaveGridSize[i];
		double xMin = -0.5;
		double xMax = 0.5;
		double yMin = -0.5;
		double yMax = 0.5;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&SquarWave);


		LaxWendroffDimSplit solver = LaxWendroffDimSplit(U, 0.5, -0.5);


		solver.solve(0.9, 10e6, 2.0);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&SquarWave, 0.5, -0.5, 2.0);


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

TEST_P(SquarWaveValidation, donerCellUpWind) {

	FILE* file;
	fopen_s(&file, "ErrorSuare_DonerCellUpWind.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(SquarWaveGridSize) / 4; i++) {

		int nx = SquarWaveGridSize[i];
		int ny = nx;
		double dx = 0.5 * (1.0 / (nx/* - 1*/));
		double xMin = -0.5 + dx;
		double xMax = 0.5 - dx;
		double yMin = -0.5 + dx;
		double yMax = 0.5 - dx;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&SquarWave);


		DonerCellUpWind solver = DonerCellUpWind(U, a, b);


		solver.solve(0.5, 10e6, tFinal);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&SquarWave, a, b, tFinal);


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

TEST_P(SquarWaveValidation, cornerTransUpWind) {

	FILE* file;
	fopen_s(&file, "ErrorSuare_CornerTransUpWind.plt", "w");

	fprintf(file, "ZONE \n");


	for (int i = 0; i < sizeof(SquarWaveGridSize) / 4; i++) {

		int nx = SquarWaveGridSize[i];
		int ny = SquarWaveGridSize[i];
		double xMin = -0.5;
		double xMax = 0.5;
		double yMin = -0.5;
		double yMax = 0.5;


		U = Field(nx, ny, xMin, xMax, yMin, yMax);
		U.init(&SquarWave);


		CornerTransUpWind solver = CornerTransUpWind(U, 0.5, -0.5);


		solver.solve(0.9, 10e6, 2.0);


		Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
		Ua.init(&SquarWave, 0.5, -0.5, 2.0);


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



TEST_P(SquarWaveValidation, Example) {

	int nx = 48;
	int ny = nx;
	double dx = 0.5 * (1.0 / (nx/* - 1*/));
	double xMin = -0.5 + dx;
	double xMax = 0.5 - dx;
	double yMin = -0.5 + dx;
	double yMax = 0.5 - dx;


	U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init(&SquarWave);


	DonerCellUpWind solver = DonerCellUpWind(U, a, b);
	//LaxWendroff solver = LaxWendroff(U, 0.5, -0.5);
	//DonerCellUpWind solver = DonerCellUpWind(U, 0.5, -0.5);
	//CornerTransUpWind solver = CornerTransUpWind(U, 0.5, -0.5);
	//LaxWendroffDimSplit solver = LaxWendroffDimSplit(U, 0.5, -0.5);


	solver.solve(0.5, 10e6, tFinal);
	U.write("numericalSolution");


	Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ua.init(&SquarWave, a, b, tFinal);
	Ua.write("analyticalSolution");


	//Calculating the L1 error
	double error = 0.0;
	for (int j = 0; j < U.get_nNode(); j++) {
		error += abs(U.var[j] - Ua.var[j]);
	}
	error = error / U.get_nNode();

	printf("%.8f  \n", error);
}








