#include "D2Q9Stencil.h"

D2Q9Stencil::D2Q9Stencil()
{
}


D2Q9Stencil::~D2Q9Stencil()
{
}


D2Q9Stencil::D2Q9Stencil(Field& field_):Stencil(field_)
{

	nNeib = 9;

	coeff1D = new double[nNeib];

}


void D2Q9Stencil::set_coeff2D(double **coeff2D_)
{
	coeff2D[1][1] = coeff2D_[1][1];
	coeff2D[1][2] = coeff2D_[1][2];
	coeff2D[2][2] = coeff2D_[2][2];
	coeff2D[2][1] = coeff2D_[2][1];
	coeff2D[2][1] = coeff2D_[2][1];
	coeff2D[1][0] = coeff2D_[1][0];
	coeff2D[0][0] = coeff2D_[0][0];
	coeff2D[0][1] = coeff2D_[0][1];
	coeff2D[0][2] = coeff2D_[0][2];

	coeff1D[0] = coeff2D_[1][1];
	coeff1D[1] = coeff2D_[1][2];
	coeff1D[2] = coeff2D_[2][2];
	coeff1D[3] = coeff2D_[2][1];
	coeff1D[4] = coeff2D_[2][0];
	coeff1D[5] = coeff2D_[1][0];
	coeff1D[6] = coeff2D_[0][0];
	coeff1D[7] = coeff2D_[0][1];
	coeff1D[8] = coeff2D_[0][2];

}


void D2Q9Stencil::get_iNeib(int iNode, int* iNeib)
{
	int iNeibTmp[9];

	field.get_iNeibPeriodicD2Q9(iNode, iNeibTmp);

	iNeib[0] = iNode;
	iNeib[1] = iNeibTmp[1];
	iNeib[2] = iNeibTmp[2];
	iNeib[3] = iNeibTmp[3];
	iNeib[4] = iNeibTmp[4];
	iNeib[5] = iNeibTmp[5];
	iNeib[6] = iNeibTmp[6];
	iNeib[7] = iNeibTmp[7];
	iNeib[8] = iNeibTmp[8];

}


