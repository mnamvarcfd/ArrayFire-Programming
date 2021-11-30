#include "Examples.h"
#include "../MAIN_POISSON2D/JacobiSolver.h"

Examples::Examples()
{
	std::cout << "Examples" << std::endl;

	tol = 10e-6;
	pi = 4.0 * atan(1.0);

}

Examples::~Examples() 
{
	std::cout << "~Examples" << std::endl;
}

void Examples::SetUp()
{
	std::cout << "SetUp" << std::endl;

	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;

}

void Examples::TearDown()
{
	std::cout << "TearDown" << std::endl;
}


double function(double x, double y) {
	return sin(x) + cos(y);
}

INSTANTIATE_TEST_CASE_P(Example, Examples, ::testing::Values(0) );


TEST_P(Examples, grid4by4) {

	int nx = 4;
	int ny = 4;
	double xMin = 0.0;
	double xMax = 3.0;
	double yMin = 0.0;
	double yMax = 3.0;


	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();


	DenseSolver denseSolver = DenseSolver(U);
	denseSolver.solve();


	EXPECT_NEAR(U.var[5], 1.42850855, tol);
	EXPECT_NEAR(U.var[6], 1.49770178, tol);
	EXPECT_NEAR(U.var[9], 0.45278586, tol);
	EXPECT_NEAR(U.var[10], 0.52197908, tol);

}


TEST_P(Examples, grid5by5) {

	int nx = 5;
	int ny = 5;
	double xMin = 0.0;
	double xMax = 4.0;
	double yMin = 0.0;
	double yMax = 4.0;


	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();


	DenseSolver denseSolver = DenseSolver(U);
	denseSolver.solve();

	EXPECT_NEAR(U.var[6], 1.432072, tol);
	EXPECT_NEAR(U.var[7], 1.507816, tol);
	EXPECT_NEAR(U.var[8], 0.711560, tol);

	EXPECT_NEAR(U.var[11], 0.456927, tol);
	EXPECT_NEAR(U.var[12], 0.528733, tol);
	EXPECT_NEAR(U.var[13],-0.267618, tol);

	EXPECT_NEAR(U.var[16], -0.142276, tol);
	EXPECT_NEAR(U.var[17], -0.075343, tol);
	EXPECT_NEAR(U.var[18], -0.862788, tol);
}


TEST_P(Examples, testCase6) {

	int nx = 7;
	int ny = 6;
	double xMin = -3.0;
	double xMax = 3.0 * pi;
	double yMin = 3.0;
	double yMax = 4.0 * pi;


	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();


	DenseSolver denseSolver = DenseSolver(U);
	denseSolver.solve();


	EXPECT_NEAR(U.var[5], 1.42850855, tol);
	EXPECT_NEAR(U.var[6], 1.49770178, tol);
	EXPECT_NEAR(U.var[9], 0.45278586, tol);
	EXPECT_NEAR(U.var[10], 0.52197908, tol);

}



TEST_P(Examples, forPresentation) {

	int nx = 50;
	int ny = 50;
	double xMin = 0.0;
	double xMax = 1.0;
	double yMin = 0.0;
	double yMax = 1.0;


	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();
	U.write("init");


	DenseSolver denseSolver = DenseSolver(U);
	denseSolver.solve();
	U.write("numericalSolution");


	Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ua.init(&function);
	Ua.write("analyticalSolution");


	double error = 0.0;
	for (int j = 0; j < U.get_nNode(); j++) {
	error += (U.var[j] - Ua.var[j])* (U.var[j] - Ua.var[j]);
	}
	printf("error is = %.8f  \n", sqrt(error));

}




TEST_P(Examples, grid5by5jacobi) {

	int nx = 5;
	int ny = 5;
	double xMin = 0.0;
	double xMax = 4.0;
	double yMin = 0.0;
	double yMax = 4.0;


	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
	U.init();
	U.setBC();


	JacobiSolver jacobiSolver = JacobiSolver(U);
	jacobiSolver.solve();

	EXPECT_NEAR(U.var[6], 1.432072, tol);
	EXPECT_NEAR(U.var[7], 1.507816, tol);
	EXPECT_NEAR(U.var[8], 0.711560, tol);

	EXPECT_NEAR(U.var[11], 0.456927, tol);
	EXPECT_NEAR(U.var[12], 0.528733, tol);
	EXPECT_NEAR(U.var[13], -0.267618, tol);

	EXPECT_NEAR(U.var[16], -0.142276, tol);
	EXPECT_NEAR(U.var[17], -0.075343, tol);
	EXPECT_NEAR(U.var[18], -0.862788, tol);
}


TEST_P(Examples, testCase1Jacobi) {

	int nx = 5;
	int ny = 5;
	double xMin = 0.0;
	double xMax = 1.0;
	double yMin = 0.0;
	double yMax = 1.0;


	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;



	Field Ud = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ud.init();
	Ud.setBC();
	DenseSolver denseSolver = DenseSolver(Ud);
	denseSolver.solve();


	Field Uj = Field(nx, ny, xMin, xMax, yMin, yMax);
	Uj.init();
	Uj.setBC();
	JacobiSolver jacobiSolver = JacobiSolver(Uj);
	jacobiSolver.solve();


	for (int j = 0; j < 25; j++) {
		EXPECT_NEAR(Ud.var[j], Uj.var[j], 10e-6);
	}

}


TEST_P(Examples, testCase2Jacobi) {

	int nx = 10;
	int ny = 10;
	double xMin = 0.0;
	double xMax = 1.0;
	double yMin = 0.0;
	double yMax = 1.0;


	af::setBackend(AF_BACKEND_CPU);
	af::setDevice(0);
	af::info(); std::cout << std::endl;



	Field Ud = Field(nx, ny, xMin, xMax, yMin, yMax);
	Ud.init();
	Ud.setBC();
	DenseSolver denseSolver = DenseSolver(Ud);
	denseSolver.solve();


	Field Uj = Field(nx, ny, xMin, xMax, yMin, yMax);
	Uj.init();
	Uj.setBC();
	JacobiSolver jacobiSolver = JacobiSolver(Uj);
	jacobiSolver.solve();


	for (int j = 0; j < 100; j++) {
		EXPECT_NEAR(Ud.var[j], Uj.var[j], 10e-6);
	}

}










//TEST_P(Examples, CSRformat) {

//
//	int nx = 5;
//	int ny = 5;
//	double xMin = 0.0;
//	double xMax = 4.0;
//	double yMin = 0.0;
//	double yMax = 4.0;
//
//
//	Field U = Field(nx, ny, xMin, xMax, yMin, yMax);
//	U.init();
//	U.setBC();
//	//U.write("init");
//
//
//	DenseSolver denseSolver = DenseSolver(U);
//	denseSolver.solve();
//
//
//	Field Ua = Field(nx, ny, xMin, xMax, yMin, yMax);
//	Ua.init(&function);
//
//	for (int j = 0; j < nx * ny; j++) {
//		printf("error[%d] = %.8f  Numeric:Anaytic: %.8f   %.8f \n", j, U.var[j] - Ua.var[j], U.var[j], Ua.var[j]);
//	}
//	double sum = 0.0;
//	for (int j = 0; j < nx * ny; j++) {
//		sum += (U.var[j] - Ua.var[j]) * (U.var[j] - Ua.var[j]);
//	}
//	printf("error is = %.8f  \n", sqrt(sum));
//
//	FAIL();
//}
//
//
//
//
//TEST_P(Examples, sparse) {
//
//	af::setBackend(AF_BACKEND_CPU);
//	af::setDevice(0);
//	af::info(); std::cout << std::endl;
//
//	float v[] = { 8, 5, 3, 6 };
//	int r[] = { 0, 0, 2, 3, 4 };
//	int c[] = { 1, 0, 2, 1 };
//	const int M = 4, N = 4, nnz = 4;
//	af::array vals = af::array(af::dim4(nnz), v);
//	af::array row_ptr = af::array(af::dim4(M + 1), r);
//	af::array col_idx = af::array(af::dim4(nnz), c);
//	// Create sparse array (CSR) from af::arrays containing values,
//	// row pointers, and column indices.
//	af::array sparse = af::sparse(M, N, vals, row_ptr, col_idx, AF_STORAGE_CSR);
//	// sparse
//	//     values:  [ 5.0, 8.0, 3.0, 6.0 ]
//	//     row_ptr: [ 0, 0, 2, 3, 4 ]
//	//     col_idx: [ 0, 1, 2, 1 ]
//
//	af_print(sparse);
//
//	af::array dens = af::dense(sparse);
//
//	af_print(dens);

	//FAIL();
//}
















