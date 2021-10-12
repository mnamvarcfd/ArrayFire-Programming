//////////////////////////////////////////////////////////////
//
//  Contributors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __VTK_H__
#define __VTK_H__

#include <arrayfire.h>
#include <stdio.h>

namespace IO
{
	template <typename T>
	int readAsciiHeaderVTK(const std::string fileName, long long int& dim0, long long int& dim1, long long int& dim2);

	template <typename T>
	int writeAsciiHeaderVTK(const std::string fileName, const std::string tag, const long long int dim0, const long long int dim1, const long long int dim2);
	
	template <typename T>
	void arrayToVTK(std::string filePrefix, const std::string tag,
		long long int dim0, long long int dim1, long long int dim2,
		const af::array& columnMajorArray, const bool isSwapEnd = true);

	template <typename T>
	af::array VTKToArray(std::string filePrefix, const bool isSwapEnd = true);
}
#endif