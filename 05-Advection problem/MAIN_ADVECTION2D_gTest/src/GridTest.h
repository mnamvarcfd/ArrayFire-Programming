#pragma once
#include <gtest/gtest.h>
#include "../../MAIN_ADVECTION2D/src/Grid.h"


class GridTest :public ::testing::TestWithParam<int> {

protected:
	double tol;
	double pi;

	int nx;
	int ny;
	double xMin;
	double xMax;
	double yMin;
	double yMax;

	Grid grid;

public:
	GridTest();
	~GridTest();
	void SetUp() override;     // gtest construction
	void TearDown() override;  // gtest destruction

};


