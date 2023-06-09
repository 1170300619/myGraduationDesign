

#ifndef RENDERER_MODEL_HPP
#define RENDERER_MODEL_HPP

#include "Vertex.hpp"
#include "Triangle.hpp"
#include "BVH.hpp"

class Material;

class Mesh
{
public:
    uint32_t numVertices;
    uint32_t numTriangles;

    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;

    Material *material;
    BVHNode *BVH;

    Mesh();

    ~Mesh();

    /**
     *
     * @param vertex
     */
    void add_vertex(Vertex &vertex);

    /**
     *
     * @param triangle
     */
    void add_triangle(Triangle &triangle);
	void createBVH();
};

class Model
{
public:
    std::vector<Mesh *> meshes;
    BVHNode* BVH;
    Model();

    ~Model();

    void createBVH();
    bool hit(Ray &ray, float tMin, float tMax, HitRecord &record) const;
    /**
     *
     * @return
     */
    int num_model() const;
};

#endif //RENDERER_MODEL_HPP
