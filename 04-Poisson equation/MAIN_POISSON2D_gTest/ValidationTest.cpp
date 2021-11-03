#include "ValidationTest.h"

ValidationTest::ValidationTest()
{
	std::cout << "ValidationTest" << std::endl;

	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

}

ValidationTest::~ValidationTest() 
{
	std::cout << "~ValidationTest" << std::endl;
}

void ValidationTest::SetUp()
{
	std::cout << "SetUp" << std::endl;
}

void ValidationTest::TearDown()
{
	std::cout << "TearDown" << std::endl;
}



INSTANTIATE_TEST_CASE_P(2DProgram, ValidationTest, ::testing::Values(0) );


TEST_P(ValidationTest, 4by4_grid) {


	int Nx = 4;
	int Ny = 4;
	double xMin = 0.0;
	double xMax = 3.0;
	double yMin = 0.0;
	double yMax = 3.0;


	Field U = Field(Nx, Ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();
	U.write("init");



	LinearSys linearSys = LinearSys(U);
	linearSys.creatCoeffMatrix();
	linearSys.writeCoeffMatrix("Matrix");

	linearSys.creatRigtHandSid();
	//printf("RHS[0] = %.4f ---->  %.4f\n", linearSys.RHS[0], -U.var[1] - U.var[4] - (sin(1.0) + cos(1.0)));
	//printf("RHS[1] = %.4f ---->  %.4f\n", linearSys.RHS[1], -U.var[7] - U.var[2] - (sin(2.0) + cos(1.0)));
	//printf("RHS[2] = %.4f ---->  %.4f\n", linearSys.RHS[2], -U.var[8] - U.var[13]- (sin(1.0) + cos(2.0)));
	//printf("RHS[3] = %.4f ---->  %.4f\n", linearSys.RHS[3], -U.var[11]- U.var[14]- (sin(2.0) + cos(2.0)));

	linearSys.solve();


	U.write("numericalSolution");



	AnalyticalSolution analyticalSolution = AnalyticalSolution(Nx, Ny, xMin, xMax, yMin, yMax);
	analyticalSolution.write("analyticalSolution");

	for (int j = 0; j < Nx * Ny; j++) {
		printf("error[%d] = %.8f  Numeric:Anaytic: %.8f   %.8f \n", j, U.var[j] - analyticalSolution.var[j], U.var[j], analyticalSolution.var[j]);
	}
	double sum = 0.0;
	for (int j = 0; j < Nx * Ny; j++) {
		sum += (U.var[j] - analyticalSolution.var[j]) * (U.var[j] - analyticalSolution.var[j]);
	}
	printf("error is = %.8f  \n", sqrt(sum));




	double x1 = 342839 / 240000.;
	double x2 = 179723 / 120000.0;
	double x3 = 10867 / 24000.;
	double x4 = 125273 / 240000.;

	printf("U[0] = %.8f ---->  %.8f \n", U.var[5], x1);
	printf("U[1] = %.8f ---->  %.8f \n", U.var[6], x2);
	printf("U[2] = %.8f ---->  %.8f \n", U.var[9], x3);
	printf("U[3] = %.8f ---->  %.8f \n", U.var[10], x4);

	EXPECT_NEAR(U.var[5], 1.42849583, 10e-4);
	EXPECT_NEAR(U.var[6], 1.49769167, 10e-4);
	EXPECT_NEAR(U.var[9], 0.45279167, 10e-4);
	EXPECT_NEAR(U.var[10], 0.52197083, 10e-4);

	FAIL();
}



TEST_P(ValidationTest, 5by5_grid) {

	int Nx = 5;
	int Ny = 5;
	double xMin = 0.0;
	double xMax = 4.0;
	double yMin = 0.0;
	double yMax = 4.0;


	Field U = Field(Nx, Ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();
	U.write("init");


	LinearSys linearSys = LinearSys(U);
	linearSys.creatCoeffMatrix();
	linearSys.writeCoeffMatrix("Matrix");
	linearSys.creatRigtHandSid();
	linearSys.solve();


	U.write("numericalSolution");



	AnalyticalSolution Ua = AnalyticalSolution(Nx, Ny, xMin, xMax, yMin, yMax);
	Ua.write("analyticalSolution");

	for (int j = 0; j < Nx * Ny; j++) {
		printf("error[%d] = %.8f  Numeric:Anaytic: %.8f   %.8f \n", j, U.var[j] - Ua.var[j], U.var[j], Ua.var[j]);
	}
	double sum = 0.0;
	for (int j = 0; j < Nx * Ny; j++) {
		sum += (U.var[j] - Ua.var[j]) * (U.var[j] - Ua.var[j]);
	}
	printf("error is = %.8f  \n", sqrt(sum));


	FAIL();
}


TEST_P(ValidationTest, sparse) {

	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

	//float v[] = { 5, 8, 3, 6 };
	//int r[] = { 0, 0, 2, 3, 4 };
	//int c[] = { 0, 1, 2, 1 };
	//const int M = 4, N = 4, nnz = 4;
	//af::array vals = af::array(af::dim4(nnz), v);
	//af::array row_ptr = af::array(af::dim4(M + 1), r);
	//af::array col_idx = af::array(af::dim4(nnz), c);
	//// Create sparse array (CSR) from af::arrays containing values,
	//// row pointers, and column indices.
	//af::array sparse = af::sparse(M, N, vals, row_ptr, col_idx, AF_STORAGE_CSR);
	//// sparse
	////     values:  [ 5.0, 8.0, 3.0, 6.0 ]
	////     row_ptr: [ 0, 0, 2, 3, 4 ]
	////     col_idx: [ 0, 1, 2, 1 ]

	//af_print(sparse);

	//af::array dens = af::dense(sparse);

	//af_print(dens);

	FAIL();
}













