#ifndef BUNDLEADJUSTMENT_H_
#define BUNDLEADJUSTMENT_H_

#include <Eigen/Dense>
#include "Util.hpp"
#include "AnimationConfig.hpp"
#include "Camera.hpp"
#include "../include/Lmap.hpp"

class BundleAdjustment {
public:
	static std::tuple<bool, float, float> calcMatching(const cv::Mat& orthoImage, const cv::Mat &image);

    static void OptSingleCamParam(const vector<cv::Mat> &coResImage, const vector<cv::Mat> &newOrthoImage,
                           AnimationConfig *animationConfig, const vector<float4> &landmarks, int index,
                           Camera *camera);

    //	void AddNoise2Cam(string OptResultCamParamPath, vector<float4> &cam_pos_list, vector<float4> &cam_hpb_list);
//	void writeCamParam(string camParamPath, vector<float4> &cam_pos_list, vector<float4> &cam_hpb_list);
};

#endif // !BUNDLEADJUSTMENT_H_