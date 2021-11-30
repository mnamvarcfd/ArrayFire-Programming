#pragma once
#include "Stencil.h"

class D2Q9Stencil :public Stencil
{

public:
	D2Q9Stencil();
	~D2Q9Stencil();

	D2Q9Stencil(Field& field_);

	void set_coeff2D(double** coeff2D_);

	void get_iNeib(int iNode, int* iNeib);
};
