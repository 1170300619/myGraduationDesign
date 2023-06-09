#pragma once
#ifndef SPHERICALSURFICQ_H_
#define SPHERICALSURFICQ_H_

#include <Eigen/Dense>
#include <cmath>
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <iostream>
#include "Model.hpp"
#include "Util.hpp"
#include "iodata.hpp"

class SphericalSurfICQ {
public:
    SphericalSurfICQ();
	~SphericalSurfICQ();

	void initFlat(int n, int landmarkSize);

	void updataLandmarkVp();

public:
	vector<float4> vertexList;
    vector<Model*> modelList;
    int vertexNum; // 代表球面上多少分布点
};


#endif