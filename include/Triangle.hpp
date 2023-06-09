

#ifndef RENDERER_TRIANGLE_HPP
#define RENDERER_TRIANGLE_HPP

#include "Vertex.hpp"
#include "Primitive.hpp"
#include <vector>

class Mesh;

class Triangle : public Primitive
{
public:
    //    Eigen::Vector3f vertex_0;
    //    Eigen::Vector3f vertex_1;
    //    Eigen::Vector3f vertex_2;

    // Index of vertices
    int vertex_0;
    int vertex_1;
    int vertex_2;

    Mesh *mesh;

    /**
     * Parameters of plane equation Ax + By + Cz + D = 0
     */
    //    float A;
    //    float B;
    //    float C;
    //    float D;

    float4 normal; // Normal vector of triangle plane

//    bool is_clip; // if is clipped

    /**
     * A world on plane
     */
    //    Eigen::Vector3f plane_point;
    Triangle();
    Triangle(const std::vector<Vertex> &vertices, int v_0, int v_1, int v_2);
    ~Triangle() override;
    bool hit(Ray &ray, float tMin, float tMax, HitRecord &record) override;
    void destroy() override;
};

#endif
