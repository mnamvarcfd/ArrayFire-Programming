#include "CornerTransUpWind.h"

CornerTransUpWind::CornerTransUpWind()
{
}

CornerTransUpWind::~CornerTransUpWind()
{
}

CornerTransUpWind::CornerTransUpWind(Field &field_, double a_, double b_):Solver<D2Q9Stencil>(field_)
{

	stencil = D2Q9Stencil(field_);
	
	a = a_;
	b = b_;
}


void CornerTransUpWind::solve(double CFL_, int nMaxIter, double totalTime)
{

	CFL = CFL_;

	timeStep();

	int nIter = /*ceil*/(totalTime / dt);

	setCoeff();

	double time = 0.0;
	int it = 0;
	while (it <= nIter)
	//while (it < 10)
	{
		if (it == nIter) dt = totalTime - nIter * dt;


		for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
			varN[iNode] = varNp1[iNode];
		}

		convolve2NaiveCpp();

		applyBC();

		it++;
	}


}

void CornerTransUpWind::solveParallel(double CFL_, double totalTime)
{

	CFL = CFL_;

	timeStep();

	int nIter = (totalTime / dt);

	setCoeff();



	reqData();

	af::array filter = stencil.getFilter();

	af::array A_d(field.get_nNode(), varNp1);

	varNp1_d = af::moddims(A_d, field.get_nx(), field.get_ny());


	int it = 0;
	while (it <= nIter)
	//while (it < 10)
	{

		if (it == nIter) dt = totalTime - nIter * dt;

		varN_d = varNp1_d;

		varNp1_d = convolve2(varN_d, filter);

		D2H(nReq, iReq, varN_d, varN);

		applyBC();

		H2D(field.nBoundNode, field.iBoundNode, varNp1, varNp1_d);

		it++;

	}

	D2H(varNp1_d, varNp1);

}


void CornerTransUpWind::setCoeff()
{
	//double dx = field.get_dx();
	//double dy = field.get_dy();
	double dx = 1. / (field.get_nx());
	double dy = 1. / (field.get_nx());

	double dxdy = dx*dy;
	double dt2dxdy = dt*dt / (dx * dy);
	double ap = max(a , 0.0);
	double am = min(a , 0.0);
	double bp = max(b , 0.0);
	double bm = min(b , 0.0);

	coeff2D[1][1] = 1.0 - dt * (ap - am) / dx - dt * (bp - bm) / dy + dt2dxdy*(ap*bp-ap*bm-am*bp+am*bm);
	coeff2D[1][2] = -dt * bm / dy + dt2dxdy * (ap * bm - am * bm);
	coeff2D[2][2] = dt2dxdy * am * bm;
	coeff2D[2][1] = -dt * am / dx + dt2dxdy * (am * bp - am * bm);
	coeff2D[2][0] = -dt2dxdy * am * bp;
	coeff2D[1][0] = dt * bp / dy + dt2dxdy * (am * bp - ap * bp);
	coeff2D[0][0] = dt2dxdy * ap * bp;
	coeff2D[0][1] = dt * ap / dy + dt2dxdy * (ap * bm - ap * bp);
	coeff2D[0][2] = -dt2dxdy * ap * bm;

	stencil.set_coeff2D(coeff2D);

}


void CornerTransUpWind::timeStep() {

	//double dx = field.get_dx();
	//double dy = field.get_dy();
	double dx = 1. / (field.get_nx());
	double dy = 1. / (field.get_nx());

	dt = CFL / max(abs(a/dx) , abs(b/dy));

}