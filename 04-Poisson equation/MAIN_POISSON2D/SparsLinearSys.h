//#pragma once
//#include "Field.h"
//#include <stdio.h>
//#include "arrayfire.h"
//class SparsLinearSys
//{
//private:
//	Field field;
//	Grid grid;
//	int nx;
//	int ny;
//
//	double dx;
//	double dy;
//
//public:
//	int nNode;
//
//	double* values;
//	int* rowIdx;
//	int* colIdx;
//
//	int nRows;
//	int nCols;
//	int nNonZero;
//
//	double* RHS;
//
//
//private:
//	void initCoeffMatrix();
//
//
//public:
//	SparsLinearSys();
//	SparsLinearSys(Field& field_);
//	~SparsLinearSys();
//
//	void creatCoeffMatrix();
//	void creatRigtHandSid(double(*func)(double x, double y));
//
//	void solve();
//
//	void AF2CPP(af::array A_d);
//
//
//
//
//	void writeCoeffMatrix(std::string fileName);
//};
//
