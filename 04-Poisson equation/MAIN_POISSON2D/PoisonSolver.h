#pragma once
#include <stdio.h>
#include <iostream>
#include <iostream>

#include "arrayfire.h"

#include "Field.h"

class PoisonSolver
{
private:
	double dx;
	double dy;

protected:
	Field field;
	int nNodeField;

	int nUnknown;
	int nCol;
	int nRow;
	int nNonZero;
	int* nNonZeroRow;
	int* colIdx;
	double* values;

	double* RHS;
	af::array CoefMatrxSprs;


private:
	double function(double x, double y);

protected:
	void creatCoeffMatrix();
	void creatRigtHandSid();
	void trnsfrCoefMatrx2AF();
	void back2Field(af::array A_d);

public:
	PoisonSolver();
	PoisonSolver(Field& field_);
	~PoisonSolver();

	virtual void solve() = 0;
};

