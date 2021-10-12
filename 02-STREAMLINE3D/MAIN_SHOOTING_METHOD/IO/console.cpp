//////////////////////////////////////////////////////////////
//
//  Contributors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __CONSOLE_CPP__
#define __CONSOLE_CPP__

#include "console.h"
#include <iomanip>
#include <string>

namespace IO
{
	template<typename T>
	void printToConsole(std::string name, const T value, const int precision, const int width)
	{
		std::cout << std::left << std::setw(width) << name << ": " << std::setw(width) << std::setprecision(precision) << value << std::endl;
	}

	template<typename T>
	void printToConsole(std::string name, std::string value, const int width)
	{
		std::cout << std::left << std::setw(width) << name << ": " << std::setw(width) << value << std::endl;
	}
}
#endif