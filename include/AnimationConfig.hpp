#pragma once
#include "Type.hpp"
#include <utility>
#include <vector>
#include "Lmap.hpp"
#include "Util.hpp"


class FrameConfig {
public:
	float4 camera_pos;//相机位置
    float4 camera_focal;//相机拍摄位置（landmark位置）
    float4 camera_up;//相机Z轴
    float4 light_pos;//光线位置（为方便计算的虚拟位置）
    float4 light_focal;//光线正对位置（landmark位置）
    float4 light_up;//光线Z轴
    float4x4 P;
	FrameConfig(float4 _camera_pos, float4 _camera_focal, float4 _camera_up, float4 _light_pos, float4 _light_focal, float4 _light_up, float4x4 _P) : camera_pos(std::move(_camera_pos)), camera_focal(std::move(_camera_focal)), camera_up(std::move(_camera_up)), light_pos(std::move(_light_pos)), light_focal(std::move(_light_focal)), light_up(std::move(_light_up)), P(std::move(_P)) {}
};

class AnimationConfig {
public:
	std::vector<FrameConfig> frames_config;//记录每一张图像对应的相机和光线参数
    vector<vector<int>> coRes;
    Render::processType type;
	AnimationConfig();
    explicit AnimationConfig(Render::processType type);

	/*
	* 析构函数
	* FullName:  AnimationConfig::~AnimationConfig
	* Access:    public 
	*/
	~AnimationConfig();


	/**
	* 对landmark点进行序列化拍照，并记录对应参数
	* FullName:  AnimationConfig::build
	* Access:    public 
	* Returns:   void
	* Param: std::vector<float4> landmarks landmarks位置
	*/
//	void build(std::vector<float4> &landmarks);
//	void albedo_build(std::vector<float4> &landmarks);
//	void test(float4 &landmarks);
    std::tuple<bool, float> orthoRectifyAndCalcTk(Lmap& lmap, float* orthoMes, int* angleMes, const cv::Mat &image, int index);
    bool orthoRectify(Lmap& lmap, float* orthoMes, const cv::Mat& image, int index) const;
    float calcTk(Lmap& lmap, const float* orthoMes, int* angleMes, vector<vector<float>>& Rp, int index);
    cv::Mat simulateLmap(const Lmap& lmap, int index);
};


