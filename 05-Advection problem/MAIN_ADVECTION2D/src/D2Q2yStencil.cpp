#include "D2Q2yStencil.h"

D2Q2yStencil::D2Q2yStencil()
{
}


D2Q2yStencil::~D2Q2yStencil()
{
}


D2Q2yStencil::D2Q2yStencil(Field& field_):Stencil(field_)
{

	nNeib = 3;

	coeff1D = new double[nNeib];

}


void D2Q2yStencil::set_coeff2D(double **coeff2D_)
{
	coeff2D[1][1] = coeff2D_[1][1];
	coeff2D[1][2] = coeff2D_[1][2];
	coeff2D[1][0] = coeff2D_[1][0];

	coeff1D[0] = coeff2D_[1][1];
	coeff1D[1] = coeff2D_[1][2];
	coeff1D[2] = coeff2D_[1][0];

}


void D2Q2yStencil::get_iNeib(int iNode, int* iNeib)
{
	int iNeibTmp[9];

	field.get_iNeibPeriodicD2Q9(iNode, iNeibTmp);

	iNeib[0] = iNode;
	iNeib[1] = iNeibTmp[1];
	iNeib[2] = iNeibTmp[5];

}


