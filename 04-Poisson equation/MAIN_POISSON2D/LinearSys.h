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

	double dx;
	double dy;


public:
	int nNode;
	double** coeff;
	double* RHS;


private:
	void initCoeffMatrix();


public:
	LinearSys();
	LinearSys(Field& field_);
	~LinearSys();

	void creatCoeffMatrix();
	void creatRigtHandSid(double(*func)(double x, double y));

	void solve();

	void AF2CPP(af::array A_d);

	void AF2CPP();

	void AF2CPP(af::array A_d, double* A_h);

	void AF2CPP(af::array A_d, double A_h);



	void writeCoeffMatrix(std::string fileName);
};

