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

	void init(double value);

	void init();
	void init(double Func(double x, double y));
	void init(double Func(double x, double y, double t), double t, bool periodic);
	void init(double Func(double x, double y), double a, double b, double t);
	void init(double Func(double x, double y), double a, double b, double t, bool periodic);
	void init(double Func(double x, double y), double t, bool periodic);
	void init(double Func(double x, double y, double t), double t);
	void setBC();
	void write(std::string fileName);

};




