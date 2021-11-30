//#include "Laplacian.h"
//
//Laplacian::Laplacian()
//{
//}
//
//Laplacian::~Laplacian()
//{
//}
//
//Laplacian::Laplacian(Field &field_, D2Q4Stencil stencil_):Solver(field_)
//{
//
//	stencil = stencil_;
//
//}
//
//
//void Laplacian::solve(int nIter)
//{
//
//	af::array Unp1(field.get_nx(), field.get_ny(), varNp1);
//		
//
//	for (int it = 0; it < nIter; it++)
//	{
//		printf("=== %d  \n", it);
//
//		af::array Un = Unp1;
//
//		setCoeff();
//
//		af::array filter = stencil.getFilter();
//
//		Unp1 = af::convolve2(Un, filter);
//
//		back2Field(Unp1);
//
//		applyBC();
//
//	}
//
//
//}
//
//
//void Laplacian::serialSolve(int nIter)
//{
//
//	for (int it = 0; it < nIter; it++)
//	{
//
//		varN = varNp1;
//
//		setCoeff();
//
//		convolve2NaiveCpp();
//
//		applyBC();
//
//	}
//
//}
//
//
//void Laplacian::setCoeff()
//{
//
//	coeff2D[1][1] = 0.0;
//	coeff2D[1][2] = 0.25;
//	coeff2D[2][1] = 0.25;
//	coeff2D[1][0] = 0.25;
//	coeff2D[0][1] = 0.25;
//
//	stencil.set_coeff2D(coeff2D);
//
//}
