#pragma once
#include <gtest/gtest.h>

#include "../../MAIN_ADVECTION2D/src/LaxWendroff.h"
#include "../../MAIN_ADVECTION2D/src/DonerCellUpWind.h"
#include "../../MAIN_ADVECTION2D/src/CornerTransUpWind.h"
#include "../../MAIN_ADVECTION2D/src/LaxWendroffDimSplit.h"

class ParallelValidation :public ::testing::TestWithParam<int>{

protected:
	Field User;
	Field Upar;

	double a;
	double b;
	double tFinal;

	double eps = 10e-14;

public:
	ParallelValidation();
	~ParallelValidation();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


