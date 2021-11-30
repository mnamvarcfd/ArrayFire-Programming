#pragma once
#include "Solver.h"
#include "Solver.cpp"
#include "D2Q2yStencil.h"

class LaxWendroffDimY : public Solver<D2Q2yStencil>
{

public:
	virtual void setCoeff();
	void timeStep();
	af::array filter;

public:
	LaxWendroffDimY();
	~LaxWendroffDimY();
	LaxWendroffDimY(Field& field_, double a, double b);

	void solve(double CFL_, int nMaxIter, double totalTime);

	void solveParallel(double CFL_, int nMaxIter, double totalTime);

};
