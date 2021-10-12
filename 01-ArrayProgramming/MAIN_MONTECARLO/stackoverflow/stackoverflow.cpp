//////////////////////////////////////////////////////////////
//
//  Contributors: 
//  1) Contribution from http://stackoverflow.com/questions/15165502/double-to-string-without-scientific-notation-or-trailing-zeros-efficiently
//  2) Contribution from https://stackoverflow.com/questions/105252
//	3) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  4) ...
//
//////////////////////////////////////////////////////////////

#ifndef __STACKOVERFLOW_CPP__
#define __STACKOVERFLOW_CPP__

#include "stackoverflow.h"

#include <iomanip> // std::setprecision
#include <sstream>

namespace stackoverflow
{
	// Thanks to useful code from http://stackoverflow.com/questions/15165502/double-to-string-without-scientific-notation-or-trailing-zeros-efficiently
	template <typename T>
	std::string dbl2str(double d)
	{
		std::stringstream ss;
		ss << std::fixed << std::setprecision(10) << d;              //convert double to string w fixed notation, hi precision
		std::string s = ss.str();                                    //output to std::string
		s.erase(s.find_last_not_of('0') + 1, std::string::npos);     //remove trailing 000s    (123.1200 => 123.12,  123.000 => 123.)
		return (s[s.size() - 1] == '.') ? s.substr(0, s.size() - 1) : s; //remove dangling decimal (123. => 123)
	}

	// Thanks to useful code from https://stackoverflow.com/questions/105252
	template <typename T>
	void SwapEnd(T& var)
	{
		char* varArray = reinterpret_cast<char*>(&var);
		for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
			std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
	}


}
#endif