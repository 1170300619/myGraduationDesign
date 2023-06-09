//
// Created by admin on 2023/2/28.
//

#ifndef SPC_PRIMITIVE_HPP
#define SPC_PRIMITIVE_HPP

#include "../include/BoundingBox.hpp"
#include "../include/HitRecord.hpp"
#include "../include/Vertex.hpp"

class Primitive
{
protected:
    Primitive() {}

public:
    BoundingBox box;
    float3 center; // primitive center coordinate in world space

    virtual ~Primitive() = default;

    /**
     *
     * @param ray
     * @param tMin
     * @param tMax
     * @param record
     * @return if hit a primitive object
     */
    virtual bool hit(Ray &ray, float tMin, float tMax, HitRecord &record)
    {
        return true;
    }

    /**
     * Destroy the node.
     * @param node
     */
    virtual void destroy()
    {
        delete this;
    }
};

#endif //SPC_PRIMITIVE_HPP
