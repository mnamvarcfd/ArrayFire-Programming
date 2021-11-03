#include "LinearSys.h"

LinearSys::LinearSys()
{
}

LinearSys::~LinearSys()
{
}

LinearSys::LinearSys(Field & field_)
{

	field = field_;
	dx = field.get_dx();
	dy = field.get_dy();

	nx = field.get_nx() - 2;
	ny = field.get_ny() - 2;

	grid = Grid(nx, ny);
	nNode = grid.get_nNode();
	

	coeff = new double[nNode*nNode];

	initCoeffMatrix();

	RHS = new double [nNode];

}

double LinearSys::function(double x, double y) {
	return -sin(x) - cos(y);
}

void LinearSys::creatCoeffMatrix()
{
	int iNeib[4];
	int nodeIdx;

	for (int iNode = 0; iNode < nNode; iNode++) {

		nodeIdx = field.getGlobInx(iNode, iNode);
		coeff[nodeIdx] = -2 * (dx * dx + dy * dy);


		grid.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];


		if (Tneib != -1) {
			nodeIdx = field.getGlobInx(iNode, Tneib);
			coeff[nodeIdx] = dx * dx;
		}

		if (Rneib != -1) {
			nodeIdx = field.getGlobInx(iNode, Rneib);
			coeff[nodeIdx] = dy * dy;
		}

		if (Bneib != -1) {
			nodeIdx = field.getGlobInx(iNode, Bneib);
			coeff[nodeIdx] = dx * dx;
		}

		if (Lneib != -1) {
			nodeIdx = field.getGlobInx(iNode, Lneib);
			coeff[nodeIdx] = dy * dy;
		}

	}

}

void LinearSys::creatRigtHandSid()
{
	int iNeib[4];

	for (int iNode = 0; iNode < nNode; iNode++) {

		double rhs = 0.0;

		int i = grid.getXInx(iNode);
		int j = grid.getYInx(iNode);

		grid.get_iNeib(iNode, iNeib);
		int Tneib = iNeib[0];
		int Rneib = iNeib[1];
		int Bneib = iNeib[2];
		int Lneib = iNeib[3];

		if (Tneib == -1) {
			int nodeIdx = field.getGlobInx(i+1, field.get_ny()-1);
			rhs += dx * dx * field.var[nodeIdx];
		}

		if (Rneib == -1) {
			int nodeIdx = field.getGlobInx(field.get_nx()-1, j+1);
			rhs += dy * dy * field.var[nodeIdx];
		}

		if (Bneib == -1) {
			int nodeIdx = field.getGlobInx(i+1, 0);
			rhs += dx * dx * field.var[nodeIdx];
		}

		if (Lneib == -1) {
			int nodeIdx = field.getGlobInx(0, j+1);
			rhs += dy * dy * field.var[nodeIdx];
		}


		double x = (i + 1) * dx;
		double y = (j + 1) * dy;


		 RHS[iNode] = -rhs + (dx * dx * dy * dy) * function(x, y);

	}

}

void LinearSys::initCoeffMatrix()
{
	int cnt = 0;
	for (int i = 0; i < nNode; i++) {
		for (int j = 0; j < nNode; j++) {
			coeff[cnt] = 0.0;
			cnt++;
		}
	}

}

void LinearSys::writeCoeffMatrix(std::string fileName) {

	std::string filename = fileName + ".plt";

	FILE* file;
	fopen_s(&file, filename.c_str(), "w");

	fprintf(file, " TITLE = 'Scalar'\n");
	fprintf(file, " VARIABLES = 'X', 'Y' , 'Scalar'   \n");
	fprintf(file, "ZONE T = \"Scalar zone\", I = %d , J = %d , F = POINT \n", nNode, nNode);

	for (int i = 0; i < nNode; i++) {
		for (int j = 0; j < nNode; j++) {
			int nodeIdx = field.getGlobInx(i, j);
			fprintf(file, "%e  %e  %e \n ", (double)i, (double)nNode-j , coeff[nodeIdx]);
		}
	}

	fclose(file);
}


void LinearSys::solve()
{

	af::array A = af::array(nNode, nNode, coeff);

	af::array b = af::array(nNode, RHS);

	af::array x = af::solve(A, b);

	AF2CPP(x);

}

void LinearSys::AF2CPP(af::array A_d)
{

	double* A_h = new double[nNode];

	A_h = A_d.host<double>();


	for (int iNode = 0; iNode < nNode; iNode++) {

		int i = grid.getXInx(iNode);
		int j = grid.getYInx(iNode);


		int ii = i + 1;
		int jj = j + 1;

		int iNodeField = field.getGlobInx(ii, jj);

		field.var[iNodeField] = A_h[iNode];
	}

}












	/*double A[9];
	A[0] = 0.1000;     A[1] = 3.1000;     A[2] = 6.1000;
	A[3] = 1.1000;     A[4] = 4.1000;     A[5] = 7.0000;
	A[6] = 2.0000;     A[7] = 5.0000;     A[8] = 8.0000;
	af::array af_A = af::array(3, 3, A);

	af_print(af_A);

	double b[3];
	b[0] = 21.9000;
	b[1] = 30.7000;
	b[2] = 39.0000;
	af::array af_b = af::array(3, b);

	af_print(af_b);

	af::array x = af::solve(af_A, af_b);

	af_print(x);

	af::array b2 = af::matmul(af_A, x);

	af_print(b2);*/