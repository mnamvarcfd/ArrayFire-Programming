#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/D2Q2yStencil.h"


class D2Q2yStencilTest :public ::testing::Test {

protected:
	Field field;
	D2Q2yStencil stencil;

public:
	D2Q2yStencilTest();
	~D2Q2yStencilTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


