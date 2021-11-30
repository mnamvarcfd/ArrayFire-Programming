#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/LaxWendroffDimSplit.h"
#include "../../src/StaticFunctions.cpp"

class LaxWendroffDimSplitTest :public ::testing::TestWithParam<int>{

protected:
	Field field;
	LaxWendroffDimSplit solver;
	Field Ua;

public:
	LaxWendroffDimSplitTest();
	~LaxWendroffDimSplitTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


