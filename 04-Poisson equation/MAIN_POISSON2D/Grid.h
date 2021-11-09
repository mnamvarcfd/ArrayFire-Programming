#pragma once
#include <stdio.h>
#include <string>
class Grid
{
protected:
	int nx;
	int ny;
	int nNode;

	double dx;
	double dy;

	double xMin; 
	double xMax; 
	double yMin; 
	double yMax;

public:
	bool *isBound;
	int* nodeInx;

public:
	Grid();
	~Grid();
	Grid(int nx_, int ny_);
	Grid(int nx_, int ny_, double xMin_, double xMax_, double yMin_, double yMax_);

	int get_nx();
	int get_ny();
	int get_nNode();
	double get_dx();
	double get_dy();

	double get_xVal(int iNode);

	double get_yVal(int iNode);

	void findBoundNodes();

	void renumbering();

	int getGlobInx(int i, int j);
	int getXInx(int iNode);
	int getYInx(int iNode);
	void get_iNeib(int iNode, int* iNeib);
	void writePLT(std::string fileName, double *scalar);
};

