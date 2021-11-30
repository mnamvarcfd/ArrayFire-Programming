#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/D2Q9Stencil.h"
#include "../../MAIN_ADVECTION2D/src/LaxWendroff.h"
#include "../../src/StaticFunctions.cpp"

class LaxWendroffTest :public ::testing::TestWithParam<int>/* , public LaxWendroff*/{

protected:
	Field field;
	D2Q9Stencil stencil;
	LaxWendroff solver;

	Field Ua;
public:
	LaxWendroffTest();
	~LaxWendroffTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


