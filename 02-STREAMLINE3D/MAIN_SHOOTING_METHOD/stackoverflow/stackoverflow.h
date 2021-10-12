//////////////////////////////////////////////////////////////
//
//  Contributors: 
//  1) Contribution from http://stackoverflow.com/questions/15165502/double-to-string-without-scientific-notation-or-trailing-zeros-efficiently
//  2) Contribution from https://stackoverflow.com/questions/105252
//	3) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  4) ...
//
//////////////////////////////////////////////////////////////

#ifndef __STACKOVERFLOW_H__
#define __STACKOVERFLOW_H__

#include <string>

namespace stackoverflow
{
	// Useful code from http://stackoverflow.com/questions/15165502/double-to-string-without-scientific-notation-or-trailing-zeros-efficiently
	template <typename T>
	std::string dbl2str(double d);

	// Thanks to https://stackoverflow.com/questions/105252
	template <typename T>
	void SwapEnd(T& var);
}
#endif