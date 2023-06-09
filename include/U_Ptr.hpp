#pragma once

#ifndef __UPTR__
#define __UPTR__
#include <iostream>
#include <utility>
#include "Type.hpp"

class U_Ptr
{
public:
	std::vector<double4> DEM;
	std::vector<double> albedoMap;
	int use;

	U_Ptr(std::vector<double4>  _DEM, std::vector<double>  _albedoMap) :DEM(std::move(_DEM)), albedoMap(std::move(_albedoMap)), use(1) {}
	~U_Ptr() = default;
};
#endif // !__UPTR__#pragma once
