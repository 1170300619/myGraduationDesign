//
// Created by admin on 2023/2/28.
//

#include "../include/Type.hpp"

#ifndef SPC_HITRECORD_HPP
#define SPC_HITRECORD_HPP


class HitRecord {
public:
    float t; // intersection t
    float4 position; // intersection point position in world space, e + t * d
    float4 color; // color buffer
    float4 flatNormal; // plane normal
    float4 normal; // intersection normal
//    float4 tangent; // tangent vector
//    float2 textureCoord; // texture coordinate
//    Material *material;
//    std::vector<bool> isInShadow; // is intersection point is in shadow of each light source

    HitRecord(): t(0), position(0, 0, 0, 1.0f), color(0, 0, 0, 1.0f), normal(0, 0, 0, 0){}

    ~HitRecord() = default;
};


#endif //SPC_HITRECORD_HPP
