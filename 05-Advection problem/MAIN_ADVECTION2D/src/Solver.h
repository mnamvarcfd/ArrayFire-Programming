#pragma once
#include "arrayfire.h"
#include "Field.h"

template <class T>
class Solver
{
public:
	Field field;

	double* varNp1;
	double* varN;
	double** coeff2D;

	double CFL;
	double dt;
	double a;
	double b;

	T stencil;
    



	int nReq;
	int* iReq;
	af::array varNp1_d;
	af::array varN_d;

protected:
	virtual void setCoeff() = 0;
	virtual void timeStep() = 0;
	void convolve2NaiveCpp();
	void applyBC();



	void arrayInterchange(double* ary1, double* ary2);

	void D2H(af::array A_d, double* A_h);
	void D2H(int nData, int* iData, af::array Data_d, double* Data_h);
	void H2D(int nData, int* iData, double* Data_h, af::array &Data_d);
	void applyBCpar();


public:
	Solver<T>();
	~Solver();
	Solver<T>(Field &field_);

	double get_dt();
	void set_dt(double dt_);

	void solveParallel(double CFL_, double totalTime);

	void reqData();
	virtual void solve(double CFL_, int nMaxIter, double totalTime) = 0;
};
