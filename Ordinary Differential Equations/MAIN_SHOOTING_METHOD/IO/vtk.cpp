//////////////////////////////////////////////////////////////
//
//  Contributors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __VTK_CPP__
#define __VTK_CPP__

#include "vtk.h"
#include "file.h"
#include "../type/type.h"
#include "../stackoverflow/stackoverflow.h"

#include <iostream> // std::cout, std::fixed
#include <fstream>
#include <istream>
#include <string>
#include <stdlib.h> // std::abort()

namespace IO
{
	template <typename T>
	int readAsciiHeaderVTK(const std::string fileName, long long int& dim0, long long int& dim1, long long int& dim2)
	{
		int filePositonAfterAsciiHeader;

		std::fstream myFile(fileName.c_str(), std::ios::in | std::ios::binary);
		if (myFile.good())
		{
			std::string myLine;
			std::getline(myFile, myLine); // (eg. "vtk DataFile Version 3.0")
			std::getline(myFile, myLine); // (eg. "My VTK file!") 

			std::string BINARY_OR_ASCII;
			std::getline(myFile, BINARY_OR_ASCII); // (eg. "BINARY OR ASCII") 

			std::getline(myFile, myLine); // (eg. "DATASET STRUCTURED_POINTS") 

			// Get DIMENSIONS:
			std::string myWord;
			std::getline(myFile, myWord, ' ');

			std::getline(myFile, myWord, ' ');
			dim0 = std::atol(myWord.c_str());

			std::getline(myFile, myWord, ' ');
			dim1 = std::atol(myWord.c_str());

			std::getline(myFile, myWord, ' ');
			dim2 = std::atol(myWord.c_str());

			std::getline(myFile, myLine); // (eg. "ORIGIN") 
			std::getline(myFile, myLine); // (eg. "SPACING") 
			std::getline(myFile, myLine); // (eg. "POINT_DATA") 
			std::getline(myFile, myLine); // (eg. "SCALARS") 
			std::getline(myFile, myLine); // (eg. "LOOKUP_TABLE default") 

			filePositonAfterAsciiHeader = myFile.tellg();

			myFile.close();
		}
		else
		{
			std::cout << fileName + " does not exist." << std::endl;
			std::abort();
		}

		return filePositonAfterAsciiHeader;
	}

	template <typename T>
	int writeAsciiHeaderVTK(const std::string fileName, const std::string tag, const long long int dim0, const long long int dim1, const long long int dim2)
	{
		int filePositonAfterAsciiHeader;

		std::fstream myFile(fileName.c_str(), std::ios::out | std::ios::trunc);

		// Header ASCII part of the vtk
		myFile << "# vtk DataFile Version 3.0" << std::endl;
		myFile << "My VTK file!" << std::endl;
		myFile << "BINARY" << std::endl;
		myFile << "DATASET STRUCTURED_POINTS" << std::endl;
		myFile << "DIMENSIONS" << " " << dim0 << " " << dim1 << " " << dim2 << std::endl;
		myFile << "ORIGIN" << " " << 0. << " " << 0. << " " << 0. << std::endl;
		myFile << "SPACING" << " " << 1. << " " << 1. << " " << 1. << std::endl;
		myFile << "POINT_DATA" << " " << dim0 * dim1 * dim2 << std::endl;
		myFile << "SCALARS" << " " << tag.c_str() << " " << "float" << " " << 1 << std::endl;
		myFile << "LOOKUP_TABLE default" << std::endl;

		filePositonAfterAsciiHeader = myFile.tellp();

		myFile.close();

		return filePositonAfterAsciiHeader;
	}

	template <typename T>
	void arrayToVTK(std::string fileName, const std::string tag,
		long long int dim0, long long int dim1, long long int dim2,
		const af::array& columnMajorArray, const bool isSwapEnd)
	{
		if (fileName.substr(fileName.find_last_of(".") + 1) != "vtk")
			fileName = fileName + ".vtk";

		if (columnMajorArray.elements() != dim0 * dim1 * dim2)
		{
			std::cout << "Not writing af::array to vtk because: columnMajorArray.elements() != dim0*dim1*dim2 ." << std::endl;
			std::abort();
		}
		
		int filePositonAfterAsciiHeader = IO::writeAsciiHeaderVTK<T>(fileName, tag, dim0, dim1, dim2);

		IO::writeBinaryColumnMajorArrayAtPosition<float>(fileName, columnMajorArray, filePositonAfterAsciiHeader, isSwapEnd);
	}

	template <typename T>
	af::array VTKToArray(std::string fileName, const bool isSwapEnd)
	{
		if (fileName.substr(fileName.find_last_of(".") + 1) != "vtk")
			fileName = fileName + ".vtk";
		
		long long int dim0;
		long long int dim1;
		long long int dim2;

		int filePositonAfterAsciiHeader = readAsciiHeaderVTK<T>(fileName, dim0, dim1, dim2);

		af::array columnMajorArray = IO::readBinaryColumnMajorArrayAtPosition<float>(fileName, dim0, dim1, dim2, 1, filePositonAfterAsciiHeader, isSwapEnd);
					
		return columnMajorArray.as(type::TYPE_AF<T>());
	}
}
#endif