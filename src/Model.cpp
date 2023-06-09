#include "../include/Model.hpp"

using namespace std;

Mesh::Mesh(): numVertices(0), numTriangles(0), material(nullptr), BVH(nullptr){}

Mesh::~Mesh()
{
    delete BVH;
    delete material;
}

void Mesh::add_vertex(Vertex &vertex){
    vertices.push_back(vertex);
    numVertices++;
}

void Mesh::add_triangle(Triangle &triangle){
    triangles.emplace_back(triangle);
    numTriangles++;
//    triangle.mesh = this;
}

int Model::num_model() const{
    return meshes.size();
}

void Mesh::createBVH()
{
    if (!BVH)
    {
        std::vector<Primitive *> primitives;
        for (auto &tri: triangles)
        {
            primitives.push_back(&tri);
        }
        BVH = new BVHNode(primitives, 0);
    }
}

Model::Model(){
    BVH = nullptr;
}

Model::~Model()
{
    for (auto mesh: meshes)
    {
        delete mesh;
    }
}

bool Model::hit(Ray &ray, float tMin, float tMax, HitRecord &record) const {
    return BVH->hit(ray, tMin, tMax, record);
}

void Model::createBVH()
{
    if (!BVH)
    {
        std::vector<Primitive *> primitives;
        for (const auto mesh: meshes)
        {
            mesh->createBVH();
            primitives.push_back(mesh->BVH);
        }
        BVH = new BVHNode(primitives, 0);
    }
}