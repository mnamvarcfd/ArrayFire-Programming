#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/D2Q4Stencil.h"
#include "../../MAIN_ADVECTION2D/src/Laplacian.h"

class LaplacianTest :public ::testing::TestWithParam<int> {

protected:
	Field field;
	D2Q4Stencil stencil;
	Laplacian laplacian;

public:
	LaplacianTest();
	~LaplacianTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


