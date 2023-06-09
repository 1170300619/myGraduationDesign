//
// Created by admin on 2023/2/28.
//
#ifndef SPC_RAY_HPP
#define SPC_RAY_HPP

#include <utility>
#include "../include/Type.hpp"

namespace Render{
    static const float Epsilon = 0.0005f;
}

class Ray{
public:
    float4 origin;
    float4 direction;

    Ray() {}
    Ray(float4 _origin, float4 _direction):origin(std::move(_origin)), direction(std::move(_direction)){}
    ~Ray() = default;
};
#endif //SPC_RAY_HPP
