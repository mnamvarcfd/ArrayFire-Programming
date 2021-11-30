#include "LaxWendroffDimSplit.h"


LaxWendroffDimSplit::LaxWendroffDimSplit()
{
}


LaxWendroffDimSplit::~LaxWendroffDimSplit()
{
}


LaxWendroffDimSplit::LaxWendroffDimSplit(Field& field_, double a_, double b_)
{
	
	laxWendroffDimX = LaxWendroffDimX(field_, a_, b_);
	laxWendroffDimY = LaxWendroffDimY(field_, a_, b_);

	field = field_;

	a = a_;
	b = b_;
}



void LaxWendroffDimSplit::solve(double CFL_, int nMaxIter, double totalTime)
{
	CFL = CFL_;

	timeStep();

	laxWendroffDimX.set_dt(dt);
	laxWendroffDimY.set_dt(dt);
	laxWendroffDimX.setCoeff();
	laxWendroffDimY.setCoeff();



	int nIter = (totalTime / dt);

	int it = 0;
	while (it <= nIter)
		//while (it <10)
	{
		if (it == nIter) dt = totalTime - nIter * dt;


		laxWendroffDimX.solve(CFL_, 1, 1.0);

		for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
			laxWendroffDimY.varN[iNode] = laxWendroffDimX.varNp1[iNode];
		}

		laxWendroffDimY.solve(CFL_, 1, 1.0);

		for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
			laxWendroffDimX.varN[iNode] = laxWendroffDimY.varNp1[iNode];
		}

		it++;
	}

	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
		field.var[iNode] = laxWendroffDimX.varN[iNode];
	}

}

void LaxWendroffDimSplit::solveParallel(double CFL_, double totalTime)
{
	CFL = CFL_;

	timeStep();

	laxWendroffDimX.set_dt(dt);
	laxWendroffDimY.set_dt(dt);
	laxWendroffDimX.setCoeff();
	laxWendroffDimY.setCoeff();


	laxWendroffDimX.reqData();
	laxWendroffDimY.reqData();

	laxWendroffDimX.filter = laxWendroffDimX.stencil.getFilter();
	laxWendroffDimY.filter = laxWendroffDimY.stencil.getFilter();


	af::array A_d(field.get_nNode(), laxWendroffDimX.varNp1);

	laxWendroffDimX.varNp1_d = af::moddims(A_d, field.get_nx(), field.get_ny());
	laxWendroffDimX.varN_d = laxWendroffDimX.varNp1_d;

	int nIter = (totalTime / dt);

	int it = 0;
	while (it <= nIter)
	//while (it <10)
	{
		if (it == nIter) dt = totalTime - nIter * dt;


		laxWendroffDimX.solveParallel(CFL_, 1, 1.0);

		laxWendroffDimY.varN_d = laxWendroffDimX.varNp1_d;

		laxWendroffDimY.solveParallel(CFL_, 1, 1.0);

		laxWendroffDimX.varN_d = laxWendroffDimY.varNp1_d;

		it++;
	}

	D2H(laxWendroffDimX.varN_d, field.var);

}

void LaxWendroffDimSplit::D2H(af::array A_d, double* A_h)
{

	double* array_h = A_d.host<double>();

	for (int i = 0; i < field.get_nNode(); i++) {
		A_h[i] = array_h[i];
	}

}


void LaxWendroffDimSplit::timeStep() {

	dt = CFL / max(abs(a / field.get_dx()), abs(b / field.get_dy()));

	//printf("Time step------------ is: %f \n", CFL / max(abs(a / field.get_dx()), abs(b / field.get_dy())));
}
