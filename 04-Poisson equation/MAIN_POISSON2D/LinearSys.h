#pragma once
#include "Field.h"
#include <stdio.h>
#include "arrayfire.h"
class LinearSys
{
private:
	Field field;
	Grid grid;
	int nx;
	int ny;
	int nNodefield;

	double dx;
	double dy;
	af::array A;

	int nUnknown;
	int nCol;
	int nRow;
	int nNonZero;
	int* nNonZeroRow;
	int* colIdx;
	double* values;

	af::array CoefMatrxSprs;

public:
	int nNode;
	double* coeff;
	double* RHS;


private:
	void initCoeffMatrix();
	void back2Field(af::array A_d);
	double function(double x, double y);
	void CoeffMatrix_Cpp2AF();
	void trnsfrCoefMatrx2AF();
public:
	LinearSys();
	LinearSys(Field& field_);
	~LinearSys();

	void creatCoeffMatrix();
	void creatCoeffMatrixSpars();
	void solveSpars(af::array& dens);
	void solveSpars();
	void CoeffMatrix2AF();
	void creatRigtHandSid();
	void solve();
	void writeCoeffMatrix(std::string fileName);
};

