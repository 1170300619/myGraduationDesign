#pragma once
#ifndef OPTIMATION_H_
#define OPTIMATION_H_
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <ceres/ceres.h>
#include <gflags/gflags.h>
#include <Eigen/Dense>
#include <iostream>
#include <glog/logging.h>
#include "Lmap.hpp"
#include "AnimationConfig.hpp"
#include "iodata.hpp"
#include <ceres/problem.h>

class OptimizeTool {
public:
    OptimizeTool();
	~OptimizeTool();

    static void optiSlope(Lmap& lmap, const vector<AnimationConfig*>& configVec, const string& path, int index);

	/**
	* ���������Ż�֮��DEMƽ���仯
	* FullName:  Optimate_tool::calcAveDEMChange
	* Access:    public 
	* Returns:   float
	* Param: Lmap & lmap lmap��Ϣ
	* Param: string & preDEMPath ǰһ��DEM�Ż�����·��
	*/
	static float calcAveAbsDiff(Lmap &lmap, string &preDEMPath);

	/**
	* ���¸߶ȣ�ʹ��landmark��ΪDEM�����
	* FullName:  Optimate_tool::updateHeight
	* Access:    public 
	* Returns:   void
	* Param: Lmap & lmap lmap��Ϣ
	*/
	static void updateHeight(Lmap &lmap);
};

#endif // !OPTIMATION_H_