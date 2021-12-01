#include "LaxWendroffDimY.h"

LaxWendroffDimY::LaxWendroffDimY()
{
}

LaxWendroffDimY::~LaxWendroffDimY()
{
}

LaxWendroffDimY::LaxWendroffDimY(Field& field_, double a_, double b_) :Solver<D2Q2yStencil>(field_)
{

	stencil = D2Q2yStencil(field_);

	a = a_;
	b = b_;
}



void LaxWendroffDimY::solve(double CFL_, int nMaxIter, double totalTime)
{

	CFL = CFL_;

	timeStep();

	convolve2NaiveCpp();

	applyBC();

}

void LaxWendroffDimY::solveParallel(double CFL_, int nMaxIter, double totalTime)
{

	CFL = CFL_;

	timeStep();

	varNp1_d = convolve2(varN_d, filter);

	D2H(nReq, iReq, varN_d, varN);

	applyBC();

	H2D(field.nBoundNode, field.iBoundNode, varNp1, varNp1_d);

}


void LaxWendroffDimY::setCoeff()
{

	double bdtdy = b * dt / (dy);
	double bdtdy2 = bdtdy * bdtdy;

	coeff2D[1][1] = 1.0 - bdtdy2;
	coeff2D[1][2] = 0.5 * (bdtdy2 - bdtdy);
	coeff2D[1][0] = 0.5 * (bdtdy2 + bdtdy);

	stencil.set_coeff2D(coeff2D);
}


void LaxWendroffDimY::timeStep() {

	dt = CFL / max(abs(a / field.get_dx()), abs(b / field.get_dy()));

}