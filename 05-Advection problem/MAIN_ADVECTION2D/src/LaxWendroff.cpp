#include "LaxWendroff.h"

LaxWendroff::LaxWendroff()
{
}

LaxWendroff::~LaxWendroff()
{
}

LaxWendroff::LaxWendroff(Field &field_, double a_, double b_):Solver<D2Q9Stencil>(field_)
{

	stencil = D2Q9Stencil(field_);
	
	a = a_;
	b = b_;
}


void LaxWendroff::solve(double CFL_, int nMaxIter, double totalTime)
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
		//printf("Iter: %d \n", it);

		for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
			varN[iNode] = varNp1[iNode];
		}

		//printf("Iter: %f \n", varN[220]);

		convolve2NaiveCpp();

		//printf("Iter: %f \n", dt);

		applyBC();

		it++;
		time += dt;
	}

	if(abs(totalTime-time)>10e-6 )printf("================================time: %f \n", time);

}

void LaxWendroff::solveParallel(double CFL_, double totalTime)
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
	//while (it <10)
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



void LaxWendroff::setCoeff()
{
	double dx = field.get_dx();
	double dy = field.get_dy();
	//double dx = 1. / (field.get_nx() - 2);
	//double dy = 1. / (field.get_nx() - 2);
	double a2 = a * a;
	double b2 = b * b;
	double dx2 = dx * dx;
	double dy2 = dy * dy;
	double dt2 = dt * dt;

	coeff2D[1][1] = 1.0 - dt2 * (a2/dx2 + b2/dy2);
	coeff2D[1][2] = 0.5 * (b2 * dt2/dy2 - b * dt/dy);
	coeff2D[2][2] = 0.25 * a * b * dt2 / (dx*dy);
	coeff2D[2][1] = 0.5 * (a2 * dt2/dx2 - a * dt / dx);
	coeff2D[2][0] = -0.25 * a * b * dt2 / (dx * dy);
	coeff2D[1][0] = 0.5 * (b2 * dt2 / dy2 + b * dt / dy);
	coeff2D[0][0] = 0.25 * a * b * dt2 / (dx * dy);
	coeff2D[0][1] = 0.5 * (a2 * dt2 / dx2 + a * dt / dx);
	coeff2D[0][2] = -0.25 * a * b * dt2 / (dx * dy);


	stencil.set_coeff2D(coeff2D);

	//printf("Time step is: %f \n", coeff2D[1][1]);

}


void LaxWendroff::timeStep() {

	double dx = field.get_dx();
	double dy = field.get_dy();
	//double dx = 1. / (field.get_nx() - 1);
	//double dy = 1. / (field.get_nx() - 1);

	dt = (CFL * min(dx , dy)) / sqrt(2 * (a * a + b * b));

	//printf("Time step------------ is: %f \n", a);

}


