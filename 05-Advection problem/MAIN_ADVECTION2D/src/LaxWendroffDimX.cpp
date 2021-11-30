#include "LaxWendroffDimX.h"

LaxWendroffDimX::LaxWendroffDimX()
{
}

LaxWendroffDimX::~LaxWendroffDimX()
{
}

LaxWendroffDimX::LaxWendroffDimX(Field &field_, double a_, double b_):Solver<D2Q2xStencil>(field_)
{

	stencil = D2Q2xStencil(field_);
	
	a = a_;
	b = b_;


}


void LaxWendroffDimX::solve(double CFL_, int nMaxIter, double totalTime)
{

	convolve2NaiveCpp();

	applyBC();

}


void LaxWendroffDimX::solveParallel(double CFL_, int nMaxIter, double totalTime)
{

		varNp1_d = convolve2(varN_d, filter);

		D2H(nReq, iReq, varN_d, varN);

		applyBC();

		H2D(field.nBoundNode, field.iBoundNode, varNp1, varNp1_d);

}


void LaxWendroffDimX::setCoeff()
{
	double dx = field.get_dx();
	double dy = field.get_dy();
	double adtdx = a * dt / dx;
	double adtdx2 = adtdx * adtdx;

	coeff2D[1][1] = 1.0 - adtdx2;
	coeff2D[2][1] = 0.5*(adtdx2 - adtdx);
	coeff2D[0][1] = 0.5*(adtdx2 + adtdx);

	//printf("set_coeff2D--LaxWendroffDimX---- is: %f \n", coeff2D[1][1]);

	stencil.set_coeff2D(coeff2D);

}


void LaxWendroffDimX::timeStep() {

	dt = CFL/max(abs(a/field.get_dx()) , abs(b/field.get_dy()) );

	//printf("Time step------------ is: %f \n", max(abs(a / field.get_dx()), abs(b / field.get_dy())));
}