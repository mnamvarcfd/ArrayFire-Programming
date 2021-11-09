#pragma once
#include <gtest/gtest.h>
#include "../MAIN_POISSON2D/Field.h"
#include "../MAIN_POISSON2D/DenseSolver.h"

class Examples :public ::testing::TestWithParam<int> {

protected:
	double tol;
	double pi;

public:
	Examples();
	~Examples();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction
	//double function(double x, double y);
};


