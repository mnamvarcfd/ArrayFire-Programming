#pragma once
#include "Solver.h"
#include "Solver.cpp"
#include "D2Q9Stencil.h"

class LaxWendroff: public Solver<D2Q9Stencil>
{

private:
	virtual void setCoeff();
	void timeStep();

public:
	LaxWendroff();
	~LaxWendroff();
	LaxWendroff(Field& field_, double a, double b);

};
