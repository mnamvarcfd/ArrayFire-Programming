//////////////////////////////////////////////////////////////
//
//  Contributors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __FILE_H__
#define __FILE_H__

#include <arrayfire.h>

namespace IO
{
	template <typename T>
	void scalarToFile(std::string fileName, const T scalarToWrite);

	template <typename T>
	T fileToScalar(std::string fileName);

	template <typename T>
	int readAsciiHeaderFile(const std::string fileName, long long int& dim0, long long int& dim1, long long int& dim2, long long int& dim3);

	template <typename T>
	int writeAsciiHeaderFile(const std::string fileName, const long long int dim0, const long long int dim1, const long long int dim2, const long long int dim3);

	template <typename T>
	af::array readBinaryColumnMajorArrayAtPosition(const std::string fileName, const long long int dim0, const long long int dim1, const long long int dim2, const long long int dim3, const int filePositonAfterAsciiHeader, const bool isSwapEnd);

	template <typename T>
	void writeBinaryColumnMajorArrayAtPosition(const std::string fileName, af::array columnMajorArray, const int filePositonAfterAsciiHeader, const bool isSwapEnd);

	template<typename T>
	void arrayToFile(std::string fileName, const af::array& arrayToWrite);

	template<typename T>
	af::array fileToArray(std::string fileName);
}
#endif