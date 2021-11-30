#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/D2Q2xStencil.h"


class D2Q2xStencilTest :public ::testing::Test {

protected:
	Field field;
	D2Q2xStencil stencil;

public:
	D2Q2xStencilTest();
	~D2Q2xStencilTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


