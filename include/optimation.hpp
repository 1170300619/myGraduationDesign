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
	* 计算两次优化之间DEM平均变化
	* FullName:  Optimate_tool::calcAveDEMChange
	* Access:    public 
	* Returns:   float
	* Param: Lmap & lmap lmap信息
	* Param: string & preDEMPath 前一次DEM优化保存路径
	*/
	static float calcAveAbsDiff(Lmap &lmap, string &preDEMPath);

	/**
	* 更新高度，使得landmark作为DEM的零点
	* FullName:  Optimate_tool::updateHeight
	* Access:    public 
	* Returns:   void
	* Param: Lmap & lmap lmap信息
	*/
	static void updateHeight(Lmap &lmap);
};

#endif // !OPTIMATION_H_