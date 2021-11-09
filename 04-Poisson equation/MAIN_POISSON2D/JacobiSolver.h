#pragma once
#include "PoisonSolver.h"

class JacobiSolver: public PoisonSolver
{

public:
	JacobiSolver();
	JacobiSolver(Field& field_);
	af::array createDiagonalSpars(double val, int nUnknown);
	af::array createSparsidentity(int nUnknown);
	~JacobiSolver();

	void solve() override;

};

