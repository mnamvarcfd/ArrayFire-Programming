#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/D2Q4Stencil.h"


class D2Q4StencilTest :public ::testing::Test {

protected:
	Field field;
	D2Q4Stencil D2Q4;

public:
	D2Q4StencilTest();
	~D2Q4StencilTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


