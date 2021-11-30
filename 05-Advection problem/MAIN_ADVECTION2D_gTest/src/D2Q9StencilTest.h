#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/D2Q9Stencil.h"


class D2Q9StencilTest :public ::testing::Test {

protected:
	Field field;
	D2Q9Stencil stencil;

public:
	D2Q9StencilTest();
	~D2Q9StencilTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


