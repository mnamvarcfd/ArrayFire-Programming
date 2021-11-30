#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/CornerTransUpWind.h"
#include "../../src/StaticFunctions.cpp"

class CornerTransUpWindTest :public ::testing::TestWithParam<int>{

protected:
	Field field;
	CornerTransUpWind solver;
	Field Ua;

public:
	CornerTransUpWindTest();
	~CornerTransUpWindTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


