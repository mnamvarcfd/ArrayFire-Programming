#include "D2Q2xStencil.h"

D2Q2xStencil::D2Q2xStencil()
{
}


D2Q2xStencil::~D2Q2xStencil()
{
}


D2Q2xStencil::D2Q2xStencil(Field& field_):Stencil(field_)
{

	nNeib = 3;

	coeff1D = new double[nNeib];

}


void D2Q2xStencil::set_coeff2D(double **coeff2D_)
{
	coeff2D[1][1] = coeff2D_[1][1];
	coeff2D[2][1] = coeff2D_[2][1];
	coeff2D[0][1] = coeff2D_[0][1];

	coeff1D[0] = coeff2D_[1][1];
	coeff1D[1] = coeff2D_[2][1];
	coeff1D[2] = coeff2D_[0][1];

}


void D2Q2xStencil::get_iNeib(int iNode, int* iNeib)
{
	int iNeibTmp[9];

	field.get_iNeibPeriodicD2Q9(iNode, iNeibTmp);

	iNeib[0] = iNode;
	iNeib[1] = iNeibTmp[3];
	iNeib[2] = iNeibTmp[7];

}


