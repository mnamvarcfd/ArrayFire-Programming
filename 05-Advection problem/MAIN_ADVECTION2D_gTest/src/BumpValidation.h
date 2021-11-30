#pragma once
#include <gtest/gtest.h>

#include "../../MAIN_ADVECTION2D/src/LaxWendroff.h"
#include "../../MAIN_ADVECTION2D/src/DonerCellUpWind.h"
#include "../../MAIN_ADVECTION2D/src/CornerTransUpWind.h"
#include "../../MAIN_ADVECTION2D/src/LaxWendroffDimSplit.h"

class BumpValidation :public ::testing::TestWithParam<int>{

protected:
	Field U;
	Field Ua;

	double a;
	double b;
	double tFinal;

public:
	BumpValidation();
	~BumpValidation();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


