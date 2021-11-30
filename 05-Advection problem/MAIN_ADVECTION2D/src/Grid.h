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
	int nBoundNode;
	bool *isBound;
	int* nodeInx;
	int* iBoundNode;

public:
	Grid();
	~Grid();
	Grid(int nx_, int ny_);
	Grid(int nx_, int ny_, double xMin_, double xMax_, double yMin_, double yMax_);

	void get_iNeibPeriodicD2Q9(int iNode, int* iNeib);

	int get_nx();
	int get_ny();
	int get_nNode();
	double get_dx();
	double get_dy();

	double get_xVal(int iNode);

	double get_yVal(int iNode);

	void findBoundNodes();

	void renumbering();

	void findnumberBoundNode();

	void findBoundNodeIndex();

	int getGlobInx(int i, int j);
	int getXInx(int iNode);
	int getYInx(int iNode);
	int getXInx(double x);
	int getYInx(double y);
	void get_iNeib(int iNode, int* iNeib);
	void writePLT(std::string fileName, double *scalar);
	void writGrid(std::string nameOfFile);
};

