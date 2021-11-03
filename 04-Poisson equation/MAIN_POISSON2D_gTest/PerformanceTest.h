#pragma once
#include <gtest/gtest.h>
#include "../MAIN_POISSON2D/Field.h"
#include "../MAIN_POISSON2D/LinearSys.h"
#include "../MAIN_POISSON2D/AnalyticalSolution.h"

class PerformanceTest :public ::testing::TestWithParam<int> {




public:
	PerformanceTest();
	~PerformanceTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


