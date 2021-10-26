//////////////////////////////////////////////////////////////
//
//  Contributors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __FILE_CPP__
#define __FILE_CPP__

#include "console.h"
#include "file.h"
#include "../type/type.h"

#include "../stackoverflow/stackoverflow.h"
#include "../stackoverflow/stackoverflow.cpp"

#include <iomanip>
#include <iostream> // std::cout, std::fixed
#include <fstream>
#include <limits>
#include <stdlib.h>  // std::abort()

#undef min
#undef max

namespace IO
{
	template <typename T>
	void scalarToFile(std::string fileName, const T scalarToWrite) {
		std::ofstream file(fileName);
		file << std::setprecision(16) << scalarToWrite;
		file.close();
	}
	
	template <typename T>
	T fileToScalar(std::string fileName) {		
		T scalarToRead;

		std::fstream myFile(fileName.c_str(), std::ios::in | std::ios::binary);
		if (myFile.good())
		{
			std::string myWord;
			std::getline(myFile, myWord, ' ');
			scalarToRead = std::atof(myWord.c_str());
			myFile.close();
		}
		else
		{
			std::cout << fileName + " does not exist." << std::endl;
			std::abort();
		}

		return scalarToRead;
	}

	template <typename T>
	int readAsciiHeaderFile(const std::string fileName, long long int& dim0, long long int& dim1, long long int& dim2, long long int& dim3)
	{
		int filePositonAfterAsciiHeader;

		std::fstream myFile(fileName.c_str(), std::ios::in | std::ios::binary);
		if (myFile.good())
		{
			// Reading DIMENSIONS from the ascii part:
			std::string myWord;

			std::getline(myFile, myWord, ' ');
			dim0 = std::atol(myWord.c_str());

			std::getline(myFile, myWord, ' ');
			dim1 = std::atol(myWord.c_str());

			std::getline(myFile, myWord, ' ');
			dim2 = std::atol(myWord.c_str());

			std::getline(myFile, myWord);
			dim3 = std::atol(myWord.c_str());

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
	int writeAsciiHeaderFile(const std::string fileName, const long long int dim0, const long long int dim1, const long long int dim2, const long long int dim3)
	{
		int filePositonAfterAsciiHeader;

		std::fstream myFile(fileName.c_str(), std::ios::out | std::ios::trunc);

		// Writing header ASCII part of the .array file
		myFile << dim0 << " " << dim1 << " " << dim2 << " " << dim3 << std::endl;

		filePositonAfterAsciiHeader = myFile.tellp();

		myFile.close();

		return filePositonAfterAsciiHeader;
	}
	
	template <typename T>
	af::array readBinaryColumnMajorArrayAtPosition(const std::string fileName, const long long int dim0, const long long int dim1, const long long int dim2, const long long int dim3, const int filePositonAfterAsciiHeader, const bool isSwapEnd)
	{
		// start timer
		af::timer start = af::timer::start();

		af::array columnMajorArray;

		std::fstream myFile(fileName.c_str(), std::ios::in | std::ios::binary);
		if (myFile.good())
		{
			myFile.seekg(filePositonAfterAsciiHeader);

			T* columnMajorArray_host = new T[dim0*dim1*dim2*dim3];

			myFile.read((char*)&columnMajorArray_host[0], dim0 * dim1 * dim2 * dim3 * sizeof(T));

			if (isSwapEnd == true)
				for (long long int iSize = 0; iSize < dim0 * dim1 * dim2 * dim3; ++iSize)
					stackoverflow::SwapEnd<T>(columnMajorArray_host[iSize]);

			columnMajorArray = af::array(dim0, dim1, dim2, dim3, columnMajorArray_host);

			delete columnMajorArray_host;

			myFile.close();
		}
		else
		{
			std::cout << fileName + " does not exist." << std::endl;
			std::abort();
		}

		double megaBytes = dim0 * dim1 * dim2 * dim3 * sizeof(T) / (double)1024 / (double)1024;
		IO::printToConsole<double>("File loaded (" + stackoverflow::dbl2str<double>(std::floor(megaBytes / af::timer::stop(start))) + " MB/s)", fileName, 15);

		return columnMajorArray;
	}

	template <typename T>
	void writeBinaryColumnMajorArrayAtPosition(const std::string fileName, af::array columnMajorArray, const int filePositonAfterAsciiHeader, const bool isSwapEnd)
	{
		// start timer
		af::timer start = af::timer::start();

		columnMajorArray = columnMajorArray.as(type::TYPE_AF<T>());

		long long int nElement = columnMajorArray.elements();

		T* columnMajorArray_host = columnMajorArray.host<T>();

		// Swap between big-endian and little-endian
		if (isSwapEnd == true)
			for (long long int iSize = 0; iSize < nElement; ++iSize)
				stackoverflow::SwapEnd<T>(columnMajorArray_host[iSize]);

		std::fstream myFile(fileName.c_str(), std::ios::out | std::ios::binary | std::ios::app);
		myFile.seekp(filePositonAfterAsciiHeader);

		myFile.write((char*)&columnMajorArray_host[0], nElement * sizeof(T));

		myFile.close();

		af::freeHost(columnMajorArray_host);

		double megaBytes = nElement * sizeof(T) / (double)1024 / (double)1024;
		IO::printToConsole<double>("File saved (" + stackoverflow::dbl2str<double>(std::floor(megaBytes / af::timer::stop(start))) + " MB/s)", fileName, 15);
	}

	template<typename T>
	void arrayToFile(std::string fileName, const af::array& arrayToWrite)
	{
		if (fileName.substr(fileName.find_last_of(".") + 1) != "array")
			fileName = fileName + ".array";

		long long int dim0 = arrayToWrite.dims(0);
		long long int dim1 = arrayToWrite.dims(1);
		long long int dim2 = arrayToWrite.dims(2);
		long long int dim3 = arrayToWrite.dims(3);

		int filePositonAfterAsciiHeader = IO::writeAsciiHeaderFile<T>(fileName, dim0, dim1, dim2, dim3);

		IO::writeBinaryColumnMajorArrayAtPosition<T>(fileName, arrayToWrite, filePositonAfterAsciiHeader, false);
	}

	template<typename T>
	af::array fileToArray(std::string fileName)
	{
		if (fileName.substr(fileName.find_last_of(".") + 1) != "array")
			fileName = fileName + ".array";

		long long int dim0;
		long long int dim1;
		long long int dim2;
		long long int dim3;
		int filePositonAfterAsciiHeader = IO::readAsciiHeaderFile<T>(fileName, dim0, dim1, dim2, dim3);

		af::array arrayToReturn = IO::readBinaryColumnMajorArrayAtPosition<T>(fileName, dim0, dim1, dim2, dim3, filePositonAfterAsciiHeader, false);

		return arrayToReturn;
	}

}
#endif