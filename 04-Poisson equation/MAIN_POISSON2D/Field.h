#pragma once
#include "Grid.h"
class Field: public Grid
{
public:
	double *var;

public:
	Field();
	~Field();
	Field(int nx, int ny, double xMin, double xMax, double yMin, double yMax);

	void init();
	void init(double Func(double x, double y));
	void setBC();
	void write(std::string fileName);

};




