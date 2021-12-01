#include "DonerCellUpWind.h"

DonerCellUpWind::DonerCellUpWind()
{
}

DonerCellUpWind::~DonerCellUpWind()
{
}

DonerCellUpWind::DonerCellUpWind(Field &field_, double a_, double b_):Solver<D2Q4Stencil>(field_)
{

	stencil = D2Q4Stencil(field_);
	
	a = a_;
	b = b_;
}


void DonerCellUpWind::solve(double CFL_, int nMaxIter, double totalTime)
{


	CFL = CFL_;

	timeStep();

	int nIter = totalTime / dt;

	setCoeff();

	double time = 0.0;
	int it = 0;
	//while (it <= nIter)
		while (it < 100)
	{
		if (it == nIter) dt = totalTime - nIter * dt;
	

		arrayInterchange(varNp1, varN);

		convolve2NaiveCpp();

		applyBC();

		it++;

		//time += dt;
	}

	//printf("time   %.8f  \n", time);

}

void DonerCellUpWind::solveParallel(double CFL_, double totalTime)
{

	CFL = CFL_;

	timeStep();

	int nIter = (totalTime / dt);

	setCoeff();


	reqData();

	af::array filter = stencil.getFilter();

	af::array A_d(field.get_nNode(), varNp1);

	varNp1_d = af::moddims(A_d, field.get_nx(), field.get_ny());

	double exetime = 0.0;

	af::timer t1 = af::timer::start();

	int it = 0;
	//while (it <= nIter)
	while (it <1000)
	{

		if (it == nIter) dt = totalTime - nIter * dt;

		varN_d = varNp1_d;

		varNp1_d = convolve2(varN_d, filter);

		//D2H(nReq, iReq, varN_d, varN);

		applyBCpar();

		//H2D(field.nBoundNode, field.iBoundNode, varNp1, varNp1_d);

		it++;

	}

	double timePar = af::timer::stop(t1);

	exetime += timePar;
	printf("Time parallel inside %g \n", timePar);

	//D2H(varNp1_d, varNp1);


}



void DonerCellUpWind::setCoeff()
{
	//double dx = field.get_dx();
	//double dy = field.get_dy();
	double dx = 1. / (field.get_nx());
	double dy = dx;

	double ap = max(a , 0.0);
	double am = min(a , 0.0);
	double bp = max(b , 0.0);
	double bm = min(b , 0.0);

	coeff2D[1][1] = 1.0 - dt * ( (ap - am) / dx + (bp - bm) / dy);
	coeff2D[1][2] = -dt * bm / dy;
	coeff2D[2][2] = 0.0;
	coeff2D[2][1] = -dt * am / dx;
	coeff2D[2][0] = 0.0;
	coeff2D[1][0] = dt * bp / dy;
	coeff2D[0][0] = 0.0;
	coeff2D[0][1] = dt * ap / dx;
	coeff2D[0][2] = 0.0;

	stencil.set_coeff2D(coeff2D);

}


void DonerCellUpWind::timeStep() {

	//double dx = field.get_dx();
	//double dy = field.get_dy();
	double dx = 1. / (field.get_nx());
	double dy = 1. / (field.get_nx());

	dt = CFL / (abs(a)/dx + abs(b)/dy);

}



void DonerCellUpWind::applyBCpar()
{
	varNp1_d = af::moddims(varNp1_d, field.get_nNode());


	af::array iBoundNode(field.nBoundNode, field.iBoundNode);

	//An array to store the index of neighborin node based on the used stencil
	int* iNeib = new int[stencil.nNeib];

	////Traversing all the nodes in the discrete domain
	//for (int i = 0; i < field.nBoundNode; i++) {
	gfor(af::seq i, 0, field.nBoundNode - 1) {

		//Just boundary nodes participate to create rifht hand side 
		//af::array iNode = iBoundNode(i);
		int iNode = iBoundNode(i).scalar<int>();
		//int iNode = field.iBoundNode[i];

		//Find the index of all the neighboring nodes
		stencil.get_iNeib(iNode, iNeib);

		varNp1_d(iNode) = 0.0;
		for (int j = 0; j < stencil.nNeib; j++) {
			int neibIdx = iNeib[j];
			af::array coef = stencil.coeff1D[j];
			varNp1_d(iNode) += coef * varN_d(neibIdx);
		}

	}


	//std::cout << field.nBoundNode << std::endl;

	//gfor(af::seq i, 0, field.nBoundNode-1) {

	//	varNp1_d(i) = 0.0;

	//}


	varNp1_d = af::moddims(varNp1_d, field.get_nx(), field.get_ny());
}