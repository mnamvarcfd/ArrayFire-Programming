//////////////////////////////////////////////////////////////
//
//  Authors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#ifndef __TYPE_H__
#define __TYPE_H__

#include <arrayfire.h>

namespace type
{
	template<typename T>
	af::dtype TYPE_AF();

	template<>
	af::dtype TYPE_AF<double>();

	template<>
	af::dtype TYPE_AF<float>();

	template<>
	af::dtype TYPE_AF <long long int> ();

	template<>
	af::dtype TYPE_AF <int>();

	template<>
	af::dtype TYPE_AF<short int>();
}
#endif