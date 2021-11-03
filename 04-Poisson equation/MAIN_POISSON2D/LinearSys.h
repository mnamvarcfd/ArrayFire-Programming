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
	double* coeff;
	double* RHS;


private:
	void initCoeffMatrix();
	void AF2CPP(af::array A_d);
	double function(double x, double y);


public:
	LinearSys();
	LinearSys(Field& field_);
	~LinearSys();

	void creatCoeffMatrix();
	void creatRigtHandSid();
	void solve();
	void writeCoeffMatrix(std::string fileName);
};

