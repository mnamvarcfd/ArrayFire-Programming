#pragma once
#include "Solver.h"
#include "Solver.cpp"
#include "D2Q2xStencil.h"

class LaxWendroffDimX: public Solver<D2Q2xStencil>
{

public:
	void setCoeff();
	void timeStep();
	af::array filter;

public:
	LaxWendroffDimX();
	~LaxWendroffDimX();
	LaxWendroffDimX(Field& field_, double a, double b);

	void solve(double CFL_, int nMaxIter, double totalTime);

	void solveParallel(double CFL_, int nMaxIter, double totalTime);

};
