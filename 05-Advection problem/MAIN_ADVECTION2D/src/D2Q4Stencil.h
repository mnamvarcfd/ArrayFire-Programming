#pragma once
#include "Stencil.h"

class D2Q4Stencil :public Stencil
{

public:
	D2Q4Stencil();
	~D2Q4Stencil();

	D2Q4Stencil(Field& field_);

	//void set_coeff2D(double coeff2D_[][]);

	void set_coeff2D(double** coeff2D_);

	void get_iNeib(int iNode, int* iNeib);
};
