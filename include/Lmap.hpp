#pragma once
#ifndef LMAP_H_
#define LMAP_H_
#include <fstream>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include "U_Ptr.hpp"
#include "Type.hpp"
#include "Util.hpp"
#include "Model.hpp"
#include <ctime>

using namespace std;
using namespace Eigen;

class Lmap {
public:
	float4 landmark;//landmark for this Lmap area
    float4 u1;//x direction of local coordinate system
    float4 u2;//y direction of local coordinate system
    float4 u3;//z direction of local coordinate system
	U_Ptr *ptr;//smart pointer


	~Lmap();

	Lmap();

	explicit Lmap(float4 &controlPoint);

	Lmap(const Lmap &lmap) :ptr(lmap.ptr), landmark(lmap.landmark), u1(lmap.u1), u2(lmap.u2), u3(lmap.u3) {
		++ptr->use;
	}

	Lmap& operator=(const Lmap &lmap);
	/**
	* 保存lmap信息
	* FullName:  Lmap::writeLmap
	* Access:    public
	* Returns:   void
	* Param: const string & LmapPath 保存路径
	*/
	void writeLmap(const string & LmapPath);
	/**
	* 读取lmap信息
	* FullName:  Lmap::readLmap
	* Access:    public
	* Returns:   void
	* Param: const string & LmapPath 保存路径
	*/
	void readLmap(const string &LmapPath);
	/**
	* 插值获取初始DEM
	* FullName:  Lmap::InterpolatePoints
	* Access:    public
	* Returns:   void
	* Param: const string & file_name 路径
	* Param: Model * model 模型信息
	* Param: Lmap & lmap lmap信息
	* Param: const Camera * camera 相机参数
	*/
	bool GenIniDEMFromBennu(Model* model, float modelGSD);
//    float4 localToBody(float4& localPoint);
	//void InterpolatePoints_v2(const string &modelPath);
	//void calcLkPixelCoordInImage(const string &renderMapPath, Mat &image,
		//float &ObjectX, float &ObjectY);

};


#endif // !LMAP_H_