#pragma warning(disable:4996);
#include "../include/SphericalSurfICQ.hpp"

using namespace std;
using namespace cv;
using namespace Render;

SphericalSurfICQ::SphericalSurfICQ(){ vertexNum = 0; }

SphericalSurfICQ::~SphericalSurfICQ() = default;

void SphericalSurfICQ::initFlat(int n, int landmarkSize) {
    vertexNum = n;
    fibonacciGrid(vertexList, vertexNum);
    for(int i = 0; i < landmarkSize; ++i){
        auto *model = new Model();
        string path = kFirstOptimationPath + "lmap_" + to_string(i) + "/procModel.ply";
        errno_t err = 0;
        if ((err = _access_s(path.c_str(), 0)) == 0){
            iodata::readPly(path, model);
            modelList.emplace_back(model);
        }
    }

#if 1 // 融合最终模型时，将南北极高于70°的区域以初始化Maplet的地形融合（即没有优化的地形）
    for(int i = 0; i < 94; ++i){
        auto *model = new Model();
        string path = "C:/SPC_Files/SPC_generate_file/bigger_than_70/lmap_" + to_string(i) + "/procModel.ply";
        errno_t err = 0;
        if ((err = _access_s(path.c_str(), 0)) == 0){
            iodata::readPly(path, model);
            modelList.emplace_back(model);
        }
    }
#endif

}


void SphericalSurfICQ::updataLandmarkVp() {
    Ray ray = Ray(kZeroFloat4, kZeroFloat4);
    float4 dir;
    for(auto & i : vertexList){
        float4 vp = float4(0, 0, 0, 0);
        int num = 0;
        HitRecord record;
        for(auto & j : modelList){
            dir = i;
            ray = Ray(float4(0, 0, 0, 1.0f), dir);
            if(j->hit(ray, 0.001, 100000.0f, record)){
                record.position = ray.origin + (record.t) * ray.direction;
                vp += record.position;
                ++num;
            }
        }
        if(num > 0){
            vp.x() /= float(num), vp.y() /= float(num), vp.z() /= float(num), vp.w() = 1.0f;
            i = vp;
        }
    }
}