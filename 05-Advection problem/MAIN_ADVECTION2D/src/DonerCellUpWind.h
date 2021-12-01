#pragma once
#include "Solver.h"
#include "Solver.cpp"
#include "D2Q4Stencil.h"

class DonerCellUpWind: public Solver<D2Q4Stencil>
{

private:
	virtual void setCoeff();
	void timeStep();

	//void applyBCpar();

public:
	DonerCellUpWind();
	~DonerCellUpWind();
	DonerCellUpWind(Field& field_, double a, double b);

	//void solve(double CFL_, int nMaxIter, double totalTime);

	//void solveParallel(double CFL_, double totalTime);

};
