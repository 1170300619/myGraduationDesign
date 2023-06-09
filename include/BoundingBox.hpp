//
// Created by admin on 2023/2/28.
//

#ifndef SPC_BOUNDINGBOX_HPP
#define SPC_BOUNDINGBOX_HPP

#include "Ray.hpp"

/*
 * AABB
 */

class BoundingBox{
public:
    float3 minPoint;
    float3 maxPoint;

    BoundingBox();

    explicit BoundingBox(float bound[6]);
    void setBound(float bound[6]);
    ~BoundingBox() = default;

    void setBound(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax);
    bool hit(const Ray &ray, float tMin, float tMax);
    static BoundingBox combine(const BoundingBox& box_0, const BoundingBox& box_1);
};

#endif //SPC_BOUNDINGBOX_HPP
