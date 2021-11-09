#pragma once
#include "PoisonSolver.h"

class DenseSolver: public PoisonSolver
{

public:
	DenseSolver();
	DenseSolver(Field& field_);
	~DenseSolver();

	void solve() override;

};

