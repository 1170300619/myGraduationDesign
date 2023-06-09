#pragma once
#include "Type.hpp"
#include <utility>
#include <vector>
#include "Lmap.hpp"
#include "Util.hpp"


class FrameConfig {
public:
	float4 camera_pos;//���λ��
    float4 camera_focal;//�������λ�ã�landmarkλ�ã�
    float4 camera_up;//���Z��
    float4 light_pos;//����λ�ã�Ϊ������������λ�ã�
    float4 light_focal;//��������λ�ã�landmarkλ�ã�
    float4 light_up;//����Z��
    float4x4 P;
	FrameConfig(float4 _camera_pos, float4 _camera_focal, float4 _camera_up, float4 _light_pos, float4 _light_focal, float4 _light_up, float4x4 _P) : camera_pos(std::move(_camera_pos)), camera_focal(std::move(_camera_focal)), camera_up(std::move(_camera_up)), light_pos(std::move(_light_pos)), light_focal(std::move(_light_focal)), light_up(std::move(_light_up)), P(std::move(_P)) {}
};

class AnimationConfig {
public:
	std::vector<FrameConfig> frames_config;//��¼ÿһ��ͼ���Ӧ������͹��߲���
    vector<vector<int>> coRes;
    Render::processType type;
	AnimationConfig();
    explicit AnimationConfig(Render::processType type);

	/*
	* ��������
	* FullName:  AnimationConfig::~AnimationConfig
	* Access:    public 
	*/
	~AnimationConfig();


	/**
	* ��landmark��������л����գ�����¼��Ӧ����
	* FullName:  AnimationConfig::build
	* Access:    public 
	* Returns:   void
	* Param: std::vector<float4> landmarks landmarksλ��
	*/
//	void build(std::vector<float4> &landmarks);
//	void albedo_build(std::vector<float4> &landmarks);
//	void test(float4 &landmarks);
    std::tuple<bool, float> orthoRectifyAndCalcTk(Lmap& lmap, float* orthoMes, int* angleMes, const cv::Mat &image, int index);
    bool orthoRectify(Lmap& lmap, float* orthoMes, const cv::Mat& image, int index) const;
    float calcTk(Lmap& lmap, const float* orthoMes, int* angleMes, vector<vector<float>>& Rp, int index);
    cv::Mat simulateLmap(const Lmap& lmap, int index);
};


