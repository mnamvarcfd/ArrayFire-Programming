#ifndef Solver_H
#define Solver_H


#include "Solver.h"

template<class T>
Solver<T>::Solver()
{
}


template<class T>
Solver<T>::~Solver()
{
}


template<class T>
Solver<T>::Solver(Field &field_)
{

	field = field_;
	varNp1 = field_.var;

	varN = new double[field.get_nNode()];

	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
		varN[iNode] = varNp1[iNode];
	}


	coeff2D = new double* [3];
	for (int i = 0; i < 3; i++)
		coeff2D[i] = new double[3];

	dx = field.get_dx();
	dy = field.get_dx();
}


template<class T>
void Solver<T>::convolve2NaiveCpp()
{

	int* iNeib = new int[stencil.nNeib];

	//Traversing all the nodes in the discrete domain
	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {

		//Find the index of all the neighboring nodes
		stencil.get_iNeib(iNode, iNeib);

		varNp1[iNode] = 0.0;
		for (int i = 0; i < stencil.nNeib; i++) {
			int j = iNeib[i];
			varNp1[iNode] += stencil.coeff1D[i] * varN[j];
		}
		
	}
}


template<class T>
void Solver<T>::solve(double CFL_, int nMaxIter, double totalTime)
{

	CFL = CFL_;

	timeStep();

	int nIter = totalTime / dt;

	setCoeff();

	double time = 0.0;
	int it = 0;
	while (it <= nIter)
	{
		//if (it == nIter) dt = totalTime - nIter * dt;


		arrayInterchange(varNp1, varN);

		convolve2NaiveCpp();

		applyBC();

		it++;

	}

}


template<class T>
void Solver<T>::solveParallel(double CFL_, double totalTime)
{

	CFL = CFL_;

	timeStep();

	int nIter = (totalTime / dt);

	setCoeff();


	reqData();

	//Transfer index of boundary nodes to the device
	af::array d(field.nBoundNode, field.iBoundNode);

	iBoundNode = d;

	af::array filter = stencil.getFilter();

	af::array A_d(field.get_nNode(), varNp1);

	varNp1_d = af::moddims(A_d, field.get_nx(), field.get_ny());

	double exetime = 0.0;

	af::timer t1 = af::timer::start();

	int it = 0;
	while (it <= nIter)
		//while (it <1000)
	{

		//if (it == nIter) dt = totalTime - nIter * dt;

		varN_d = varNp1_d;

		varNp1_d = convolve2(varN_d, filter);

		applyBCtransData();

		it++;

	}

	D2H(varNp1_d, varNp1);


}


template<class T>
void Solver<T>::applyBC()
{
	//An array to store the index of neighborin node based on the used stencil
	int* iNeib = new int[stencil.nNeib];

	//Traversing all the nodes in the discrete domain
	for (int i = 0; i < field.nBoundNode; i++) {

		int iNode = field.iBoundNode[i];

		//Find the index of all the neighboring nodes
		stencil.get_iNeib(iNode, iNeib);

		varNp1[iNode] = 0.0;
		for (int ii = 0; ii < stencil.nNeib; ii++) {
			int j = iNeib[ii];
			varNp1[iNode] += stencil.coeff1D[ii] * varN[j];
		}

	}


}

template<class T>
void Solver<T>::applyBCtransData()
{
	D2H(nReq, iReq, varN_d, varN);

	applyBC();

	H2D(field.nBoundNode, field.iBoundNode, varNp1, varNp1_d);

}



template<class T>
void Solver<T>::applyBCgfor()
{
	varNp1_d = af::moddims(varNp1_d, field.get_nNode());

	//An array to store the index of neighborin node based on the used stencil
	int* iNeib = new int[stencil.nNeib];

	//Traversing all the nodes in the discrete domain
	gfor(af::seq i, 0, field.nBoundNode - 1) {

		//Get index of boundary node
		int iNode = iBoundNode(i).scalar<int>();

		//Find the indexes of all the neighboring nodes
		stencil.get_iNeib(iNode, iNeib);

		//Perform the computation for a bounday node
		varNp1_d(iNode) = 0.0;
		for (int j = 0; j < stencil.nNeib; j++) {
			int neibIdx = iNeib[j];
			af::array coef = stencil.coeff1D[j];
			varNp1_d(iNode) += coef * varN_d(neibIdx);
		}

	}


	varNp1_d = af::moddims(varNp1_d, field.get_nx(), field.get_ny());
}


template<class T>
double Solver<T>::get_dt()
{

	return dt;

}

template<class T>
void Solver<T>::set_dt(double dt_)
{
	dt = dt_;
}


template<class T>
void Solver<T>::arrayInterchange(double* ary1, double* ary2)
{

	double* tmp = varNp1;
	varNp1 = varN;
	varN = tmp;

}

template<class T>
void Solver<T>::D2H(af::array A_d, double * A_h)
{

	double* array_h = A_d.host<double>();

	for (int i = 0; i < field.get_nNode(); i++) {
		A_h[i] = array_h[i];
	}

}

template<class T>
void Solver<T>::D2H(int nData, int *iData, af::array Data_d, double* Data_h)
{

	af::array iData_d(nData, iData);

	af::array vData_d = Data_d(iData_d);

	double* array_h = vData_d.host<double>();
	for (int i = 0; i < nData; i++) {
		int indx = iData[i];
		Data_h[indx] = array_h[i];
	}

}


template<class T>
void Solver<T>::H2D(int nData, int* iData, double* Data_h, af::array &Data_d)
{

	double* vData_h = new double[nData];

	for (int i = 0; i < nData; i++) {

		int idx = iData[i];
		vData_h[i] = Data_h[idx];

	}

	af::array vData_d(nData, vData_h);

	af::array iData_d(nData, iData);

	Data_d(iData_d) = vData_d;

}




template<class T>
void Solver<T>::reqData()
{
	int* iNeib = new int[stencil.nNeib];

	nReq = 0;
	bool* isReq = new bool[field.get_nNode()];

	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {

		if (field.isBound[iNode] == false) continue;

		stencil.get_iNeib(iNode, iNeib);


		for (int i = 0; i < stencil.nNeib; i++) {
			int j = iNeib[i];
			isReq[j] = true;
			nReq++;
		}

	}

	nReq = 0;
	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {

		if (isReq[iNode] == true) nReq++;;

	}

	iReq = new int[nReq];

	int cnt = 0;
	for (int i = 0; i < field.get_nNode(); i++) {

		if (isReq[i] == true) {
			iReq[cnt] = i;
			cnt++;
		}

	}

}

#endif