#pragma once
#include "Stencil.h"

class D2Q2yStencil :public Stencil
{

public:
	D2Q2yStencil();
	~D2Q2yStencil();

	D2Q2yStencil(Field& field_);

	void set_coeff2D(double** coeff2D_);

	void get_iNeib(int iNode, int* iNeib);
};
