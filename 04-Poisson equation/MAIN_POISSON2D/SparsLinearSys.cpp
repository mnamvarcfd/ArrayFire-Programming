//#include "SparsLinearSys.h"
//
//SparsLinearSys::SparsLinearSys()
//{
//}
//
//SparsLinearSys::~SparsLinearSys()
//{
//}
//
//SparsLinearSys::SparsLinearSys(Field & field_)
//{
//
//	field = field_;
//	dx = field.get_dx();
//	dy = field.get_dy();
//
//	nx = field.get_nx() - 2;
//	ny = field.get_ny() - 2;
//
//	grid = Grid(nx, ny);
//	nNode = grid.get_nNode();
//
//	//coeff = new double*[nNode];
//	//for (int i = 0; i < nNode; i++) {
//	//	coeff[i] = new double[nNode];
//	//}
//
//	//initCoeffMatrix();
//
//
//	RHS = new double [nNode];
//
//}
//
//void SparsLinearSys::creatCoeffMatrix()
//{
//
//	int iNeib[4];
//
//	for (int iNode = 0; iNode < nNode; iNode++) {
//
//		coeff[iNode][iNode] = -2 * (dx * dx + dy * dy);
//
//
//		grid.get_iNeib(iNode, iNeib);
//		int Tneib = iNeib[0];
//		int Rneib = iNeib[1];
//		int Bneib = iNeib[2];
//		int Lneib = iNeib[3];
//
//
//		if (Tneib != -1) {
//			coeff[iNode][Tneib] = dx * dx;
//		}
//
//		if (Rneib != -1) {
//			coeff[iNode][Rneib] = dy * dy;
//		}
//
//		if (Bneib != -1) {
//			coeff[iNode][Bneib] = dx * dx;
//		}
//
//		if (Lneib != -1) {
//			coeff[iNode][Lneib] = dy * dy;
//		}
//
//	}
//
//
//}
//
//void SparsLinearSys::creatRigtHandSid(double(*func)(double x, double y))
//{
//	int iNeib[4];
//
//
//	//for (int iNode = 0; iNode < 16; iNode++) {
//	//	printf("---------------------- = %.4f\n", field.var[iNode]);
//	//}
//
//	for (int iNode = 0; iNode < nNode; iNode++) {
//
//		//if (iNode == 0)printf("RHS = %.4f\n", rhs);
//
//		double rhs = 0.0;
//
//		int i = grid.getXInx(iNode);
//		int j = grid.getYInx(iNode);
//
//		grid.get_iNeib(iNode, iNeib);
//
//		int Tneib = iNeib[0];
//		if (Tneib == -1) {
//			int nodeIdx = field.getGlobInx(i+1, field.get_ny()-1);
//			rhs += dx * dx * field.var[nodeIdx];
//		}
//		//if (iNode == 0)printf("RHS_Tneib = %.4f\n", rhs);
//
//		int Rneib = iNeib[1];
//		if (Rneib == -1) {
//			int nodeIdx = field.getGlobInx(field.get_nx()-1, j+1);
//			rhs += dy * dy * field.var[nodeIdx];
//		}
//		//if (iNode == 0)printf("RHS_Rneib = %.4f\n", rhs);
//
//		int Bneib = iNeib[2];
//		if (Bneib == -1) {
//			int nodeIdx = field.getGlobInx(i+1, 0);
//			rhs += dx * dx * field.var[nodeIdx];
//			//if (iNode == 0) {
//			//	printf("ij_Bneib = %d   %d\n", i,j);
//			//	printf("idx_Bneib = %d\n", nodeIdx);
//			//	printf("RHS_Bneib = %.4f\n", rhs);
//			//}
//		}
//
//		int Lneib = iNeib[3];
//		if (Lneib == -1) {
//			int nodeIdx = field.getGlobInx(0, j+1);
//			rhs += dy * dy * field.var[nodeIdx];
//		}
//		//if (iNode == 0)printf("RHS_Lneib = %.4f\n", rhs);
//
//
//		double x = (i + 1) * dx;
//		double y = (j + 1) * dy;
//
//
//		 RHS[iNode] = -rhs + (dx * dx * dy * dy) * func(x, y);
//		 //if (iNode == 0) {
//			// printf("x y = %.4f  %.4f \n", x, y);
//			// printf("-sin(x) - cos(y) = %.4f  %.4f\n", -sin(x) ,- cos(y));
//			// printf("rhs func = %.4f  %.4f\n", rhs, func(x, y));
//		 //}
//
//	}
//
//}
//
//void SparsLinearSys::initCoeffMatrix()
//{
//
//	for (int i = 0; i < nNode; i++) {
//		for (int j = 0; j < nNode; j++) {
//			coeff[i][j] = 0.0;
//		}
//	}
//
//}
//
//void SparsLinearSys::writeCoeffMatrix(std::string fileName) {
//
//	std::string filename = fileName + ".plt";
//
//	FILE* file;
//	fopen_s(&file, filename.c_str(), "w");
//
//	fprintf(file, " TITLE = 'Scalar'\n");
//	fprintf(file, " VARIABLES = 'X', 'Y' , 'Scalar'   \n");
//	fprintf(file, "ZONE T = \"Scalar zone\", I = %d , J = %d , F = POINT \n", nNode, nNode);
//
//	for (int i = 0; i < nNode; i++) {
//		for (int j = 0; j < nNode; j++) {
//			fprintf(file, "%e  %e  %e \n ", (double)i, (double)nNode-j , coeff[i][j]);
//		}
//	}
//	fclose(file);
//}
//
//
//void SparsLinearSys::solve()
//{
//
//
//}
//
//void SparsLinearSys::AF2CPP(af::array A_d)
//{
//
//	double* A_h = new double[nNode];
//
//	A_h = A_d.host<double>();
//
//
//	for (int iNode = 0; iNode < nNode; iNode++) {
//
//		int i = grid.getXInx(iNode);
//		int j = grid.getYInx(iNode);
//
//
//		int ii = i + 1;
//		int jj = j + 1;
//
//		int iNodeField = field.getGlobInx(ii, jj);
//
//		field.var[iNodeField] = A_h[iNode];
//	}
//
//}
//	/*double A[9];
//	A[0] = 0.1000;     A[1] = 3.1000;     A[2] = 6.1000;
//	A[3] = 1.1000;     A[4] = 4.1000;     A[5] = 7.0000;
//	A[6] = 2.0000;     A[7] = 5.0000;     A[8] = 8.0000;
//	af::array af_A = af::array(3, 3, A);
//
//	af_print(af_A);
//
//	double b[3];
//	b[0] = 21.9000;
//	b[1] = 30.7000;
//	b[2] = 39.0000;
//	af::array af_b = af::array(3, b);
//
//	af_print(af_b);
//
//	af::array x = af::solve(af_A, af_b);
//
//	af_print(x);
//
//	af::array b2 = af::matmul(af_A, x);
//
//	af_print(b2);*/