#pragma once
#include "Solver.h"
#include "Solver.cpp"
#include "D2Q9Stencil.h"

class CornerTransUpWind: public Solver<D2Q9Stencil>
{

private:
	virtual void setCoeff();
	void timeStep();

public:
	CornerTransUpWind();
	~CornerTransUpWind();
	CornerTransUpWind(Field& field_, double a, double b);

	void solve(double CFL_, int nMaxIter, double totalTime);

	void solveParallel(double CFL_, double totalTime);

};
