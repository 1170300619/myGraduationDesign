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
	* ����lmap��Ϣ
	* FullName:  Lmap::writeLmap
	* Access:    public
	* Returns:   void
	* Param: const string & LmapPath ����·��
	*/
	void writeLmap(const string & LmapPath);
	/**
	* ��ȡlmap��Ϣ
	* FullName:  Lmap::readLmap
	* Access:    public
	* Returns:   void
	* Param: const string & LmapPath ����·��
	*/
	void readLmap(const string &LmapPath);
	/**
	* ��ֵ��ȡ��ʼDEM
	* FullName:  Lmap::InterpolatePoints
	* Access:    public
	* Returns:   void
	* Param: const string & file_name ·��
	* Param: Model * model ģ����Ϣ
	* Param: Lmap & lmap lmap��Ϣ
	* Param: const Camera * camera �������
	*/
	bool GenIniDEMFromBennu(Model* model, float modelGSD);
//    float4 localToBody(float4& localPoint);
	//void InterpolatePoints_v2(const string &modelPath);
	//void calcLkPixelCoordInImage(const string &renderMapPath, Mat &image,
		//float &ObjectX, float &ObjectY);

};


#endif // !LMAP_H_