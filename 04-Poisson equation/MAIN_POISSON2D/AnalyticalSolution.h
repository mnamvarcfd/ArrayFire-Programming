#pragma once
#include "Field.h"
class AnalyticalSolution: public Field
{

public:
	AnalyticalSolution();
	~AnalyticalSolution();
	AnalyticalSolution(int nx, int ny, double xMin, double xMax, double yMin, double yMax);
};

