#pragma once
#include <gtest/gtest.h>
#include "../MAIN_POISSON2D/Field.h"
#include "../MAIN_POISSON2D/LinearSys.h"
#include "../MAIN_POISSON2D/AnalyticalSolution.h"

class Performance :public ::testing::TestWithParam<af::Backend> {

protected:
	FILE* file1;
public:
	Performance();
	~Performance();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


