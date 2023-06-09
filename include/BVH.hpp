//
// Created by admin on 2023/2/28.
//

#ifndef SPC_BVH_HPP
#define SPC_BVH_HPP

#include "Type.hpp"
#include "Primitive.hpp"
#include "BoundingBox.hpp"
#include "Ray.hpp"
#include "Triangle.hpp"
#include <iostream>

class BVHNode: public Primitive
{
public:
    Primitive *left;
    Primitive *right;

    BVHNode();

    /**
     *
     * @param primitives
     * @param axis
     */
    BVHNode(std::vector<Primitive *> primitives, int axis);

    ~BVHNode() override;

    /**
     *
     * @param ray
     * @param tMin
     * @param tMax
     * @param record
     * @return
     */
    bool hit(Ray &ray, float tMin, float tMax, HitRecord &record) override;

    /**
     * Destroy the node.
     */
    void destroy() override;
};

#endif //SPC_BVH_HPP
