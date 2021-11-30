#pragma once
//#include "Solver.h"
//#include "Solver.cpp"
#include "LaxWendroffDimX.h"
#include "LaxWendroffDimY.h"

class LaxWendroffDimSplit/*: public Solver()*/
{

public:
	LaxWendroffDimX laxWendroffDimX;
	LaxWendroffDimY laxWendroffDimY;

	Field field;
	double CFL;
	double a;
	double b;
public:
	LaxWendroffDimSplit();
	~LaxWendroffDimSplit();
	LaxWendroffDimSplit(Field& field_, double a, double b);

	void solve(double CFL_, int nMaxIter, double totalTime);

	void solveParallel(double CFL_, double totalTime);

	void D2H(af::array A_d, double* A_h);

	void timeStep();

	double dt;
};
