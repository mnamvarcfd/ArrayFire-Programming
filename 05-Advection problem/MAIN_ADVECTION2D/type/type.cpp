//////////////////////////////////////////////////////////////
//
//  Authors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __TYPE_CPP__
#define __TYPE_CPP__

#include "type.h"

#include <iostream>
#include <stdlib.h>  // std::abort()

namespace type
{
	template<typename T>
	af::dtype TYPE_AF() { return NULL; } // Actual value NULL is never returned because of template specialization, need to return a value to remove warning on g++

	template<>
	af::dtype TYPE_AF<double>()
	{
		return af::dtype::f64;
	}

	template<>
	af::dtype TYPE_AF<float>()
	{
		return af::dtype::f32;
	}

	template<>
	af::dtype TYPE_AF<long long int>()
	{
		return af::dtype::s64;
	}

	template<>
	af::dtype TYPE_AF<int>()
	{
		return af::dtype::s32;
	}

	template<>
	af::dtype TYPE_AF<short int>()
	{
		return af::dtype::s16;
	}
}
#endif