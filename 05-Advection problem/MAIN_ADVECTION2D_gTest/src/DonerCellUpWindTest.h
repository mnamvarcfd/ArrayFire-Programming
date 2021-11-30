#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/DonerCellUpWind.h"
#include "../../src/StaticFunctions.cpp"

class DonerCellUpWindTest :public ::testing::TestWithParam<int>{

protected:
	Field field;
	DonerCellUpWind solver;
	Field Ua;

public:
	DonerCellUpWindTest();
	~DonerCellUpWindTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


