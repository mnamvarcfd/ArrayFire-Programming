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





//void DonerCellUpWind::applyBCpar()
//{
//	varNp1_d = af::moddims(varNp1_d, field.get_nNode());
//
//	//Transfer index of boundary nodes to the device
//	af::array iBoundNode(field.nBoundNode, field.iBoundNode);
//
//	//An array to store the index of neighborin node based on the used stencil
//	int* iNeib = new int[stencil.nNeib];
//
//	//Traversing all the nodes in the discrete domain
//	gfor(af::seq i, 0, field.nBoundNode - 1) {
//
//		//Get index of boundary node
//		int iNode = iBoundNode(i).scalar<int>();
//
//		//Find the indexes of all the neighboring nodes
//		stencil.get_iNeib(iNode, iNeib);
//
//		//Perform the computation for a bounday node
//		varNp1_d(iNode) = 0.0;
//		for (int j = 0; j < stencil.nNeib; j++) {
//			int neibIdx = iNeib[j];
//			af::array coef = stencil.coeff1D[j];
//			varNp1_d(iNode) += coef * varN_d(neibIdx);
//		}
//
//	}
//
//
//	varNp1_d = af::moddims(varNp1_d, field.get_nx(), field.get_ny());
//}




void DonerCellUpWind::setCoeff()
{

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

	dt = CFL / (abs(a)/dx + abs(b)/dy);

}
