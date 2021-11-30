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
void Solver<T>::applyBC()
{
	//An array to store the index of neighborin node based on the used stencil
	int* iNeib = new int[stencil.nNeib];

	//Traversing all the nodes in the discrete domain
	for (int i = 0; i < field.nBoundNode; i++) {

		//Just boundary nodes participate to create rifht hand side 
		int iNode = field.iBoundNode[i];

		//Find the index of all the neighboring nodes
		stencil.get_iNeib(iNode, iNeib);

		varNp1[iNode] = 0.0;
		for (int i = 0; i < stencil.nNeib; i++) {
			int neibIdx = iNeib[i];
			varNp1[iNode] += stencil.coeff1D[i] * varN[neibIdx];
		}

	}


}

//template<class T>
//void Solver<T>::applyBC()
//{
//	//An array to store the index of neighborin node based on the used stencil
//	int* iNeib = new int[stencil.nNeib];
//
//	//Traversing all the nodes in the discrete domain
//	for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
//
//		//Just boundary nodes participate to create rifht hand side 
//		if (field.isBound[iNode] == false) continue;
//
//		//Find the index of all the neighboring nodes
//		stencil.get_iNeib(iNode, iNeib);
//
//		varNp1[iNode] = 0.0;
//		for (int i = 0; i < stencil.nNeib; i++) {
//			int j = iNeib[i];
//			varNp1[iNode] += stencil.coeff1D[i] * varN[j];
//		}
//
//	}
//
//
//}

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
//
//template<class T>
//void Solver<T>::solveParallel(double CFL_, int nMaxIter, double totalTime)
//{
//	int nx = 21;
//	int ny = 21;
//
//	af::array Unp1(nx, ny, varNp1);
//	af::array Un(nx, ny, varN);
//
//	af::array filter(3, 3, stencil.coeff1D);
//
//
//	CFL = CFL_;
//
//	timeStep();
//
//	int nIter = /*ceil*/(totalTime / dt);
//
//	setCoeff();
//
//	double time = 0.0;
//	int it = 0;
//	while (it <= nIter && it < nMaxIter)
//	{
//		if (it == nIter) dt = totalTime - nIter * dt;
//
//		Un = Unp1;
//
//		Unp1 = convolve2(Un, filter);
//
//		varNp1 = Un.host<double>();
//		for (int i = 0; i < field.get_nNode(); i++) {
//			varN[i] = varNp1[i];
//		}
//
//
//		applyBC();
//
//
//		Unp1 = varNp1;
//
//		it++;
//		time += dt;
//	}
//
//	/*if (abs(totalTime - time) > 10e-6)*/printf("================================time: %f \n", time);
//
//	back2Field(Unp1);
//
//}
//


template<class T>
void Solver<T>::arrayInterchange(double* ary1, double* ary2)
{

	//for (int iNode = 0; iNode < field.get_nNode(); iNode++) {
	//	ary2[iNode] = ary1[iNode];
	//}

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


	//printf("nReq %d \n", nReq);

	iReq = new int[nReq];

	int cnt = 0;
	for (int i = 0; i < field.get_nNode(); i++) {

		if (isReq[i] == true) {
			iReq[cnt] = i;
			//printf("nReq %d \n", i);

			cnt++;
		}

	}

}


template<class T>
void Solver<T>::applyBCpar() {

	//reqData();


	//D2H(nReq, iReq, varN_d, varN);


	//applyBC();


	//H2D(field.nBoundNode, field.iBoundNode, varNp1, varNp1_d);

}








#endif