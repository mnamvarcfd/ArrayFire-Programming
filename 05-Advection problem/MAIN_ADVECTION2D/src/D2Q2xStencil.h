#pragma once
#include "Stencil.h"

class D2Q2xStencil :public Stencil
{

public:
	D2Q2xStencil();
	~D2Q2xStencil();

	D2Q2xStencil(Field& field_);

	void set_coeff2D(double** coeff2D_);

	void get_iNeib(int iNode, int* iNeib);
};
