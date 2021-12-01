#include "Stencil.h"


Stencil::Stencil()
{
}
Stencil::~Stencil()
{
}


Stencil::Stencil(Field& field_)
{
	field = field_;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			coeff2D[i][j] = 0.0;

}



af::array Stencil::getFilter()
{

	for (int i = 0; i < 9; i++)
		coeff_af[i] = 0.0;


	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 3; i++) {
			creatAFval(i, j, coeff2D[i][j]);
		}
	}

	af::array filter(3, 3, coeff_af);

	return filter;
}


void Stencil::creatAFval(int i, int j, double val)
{
	int ii = 0;
	int jj = 0;

	if (i == 0 && j == 0) {
		ii = 2;
		jj = 2;
	}
	else if (i == 1 && j == 0) {
		ii = 1;
		jj = 2;
	}
	else if (i == 2 && j == 0) {
		ii = 0;
		jj = 2;
	}

	else if (i == 0 && j == 1) {
		ii = 2;
		jj = 1;
	}
	else if (i == 1 && j == 1) {
		ii = 1;
		jj = 1;
	}
	else if (i == 2 && j == 1) {
		ii = 0;
		jj = 1;
	}

	else if (i == 0 && j == 2) {
		ii = 2;
		jj = 0;
	}
	else if (i == 1 && j == 2) {
		ii = 1;
		jj = 0;
	}
	else if (i == 2 && j == 2) {
		ii = 0;
		jj = 0;
	}

	int index = jj * 3 + ii;
	coeff_af[index] = val;
}







