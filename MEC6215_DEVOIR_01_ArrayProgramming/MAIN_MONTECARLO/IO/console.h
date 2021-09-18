//////////////////////////////////////////////////////////////
//
//  Contributors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __CONSOLE_H__
#define __CONSOLE_H__

#include <iostream>

namespace IO
{
	template<typename T>
	void printToConsole(std::string name, const T value, const int precision = 16, const int width = 25);

	template<typename T>
	void printToConsole(std::string name, std::string value, const int width = 25);
}
#endif