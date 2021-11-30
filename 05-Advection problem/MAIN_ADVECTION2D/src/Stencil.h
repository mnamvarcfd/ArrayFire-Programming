#pragma once
#include <vector>
#include <iostream>

#include "arrayfire.h"

#include "Field.h"

class Stencil
{
private:
	double coeff_af[9];

protected:
	Field field;
	double coeff2D[3][3];

public:
	int nNeib;
	double *coeff1D;

private:
	void creatAFval(int i, int j, double val);

public:
	Stencil();
	~Stencil();
	Stencil(Field& field_);
	af::array getFilter();

	virtual void set_coeff2D(double** coeff2D_) = 0;
	virtual void get_iNeib(int iNode, int* iNeib) = 0;

};
