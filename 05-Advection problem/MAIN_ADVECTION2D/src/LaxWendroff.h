#pragma once
#include "Solver.h"
#include "Solver.cpp"
#include "D2Q9Stencil.h"

class LaxWendroff: public Solver<D2Q9Stencil>
{

private:
	virtual void setCoeff();
	void timeStep();

	//void applyBCpar();

	//void applyBCpar();


public:
	LaxWendroff();
	~LaxWendroff();
	LaxWendroff(Field& field_, double a, double b);

	void solve(double CFL_, int nMaxIter, double totalTime);

	void solveParallel(double CFL_, double totalTime);

};
