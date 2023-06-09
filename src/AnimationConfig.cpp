#include "../include/AnimationConfig.hpp"

using namespace std;
using namespace Render;

AnimationConfig::AnimationConfig() = default;


AnimationConfig::~AnimationConfig() = default;

bool AnimationConfig::orthoRectify(Lmap &lmap, float *orthoMes, const cv::Mat &image, int index) const {
    int Xleft, Ydown, Xright, Yup;
    float lu, ld, ru, rd;
    float f1, f2;

    for(int i = kBorderSize; i < kDEMSize - kBorderSize; ++i){
        for(int j = kBorderSize; j < kDEMSize - kBorderSize; ++j){
            float4 localPoint = double4ToFloat4(lmap.ptr->DEM[i * kDEMSize + j]);
            float4 bodyPoint = localToBody(lmap.u1, lmap.u2, lmap.u3, lmap.landmark, localPoint);
            float4 result = frames_config[index].P * bodyPoint; // 这里需要使用本体坐标系（因为P矩阵是本体坐标系下的），下面角度和Rk的计算都用局部坐标系
            result /= result.w();
            float _X = result.x(); // 图像上x轴方向坐标，代表图像上的列
            float _Y = result.y(); // 图像上y轴方向坐标，代表图像上的行
//            vector<int> a(2, 0);
//            a[0] = floor(_X), a[1] = floor(_Y);
//            axis.emplace_back(a);

            if(_X > 0 && _X < kImageSize - 1 && _Y > 0 && _Y < kImageSize - 1){
                Xleft = floor(_X);
                Ydown = floor(_Y);
                Xright = ceil(_X);
                Yup = ceil(_Y);
                //使用观测图像和阴影图像对应像素点相乘后双线性插值
                lu = (float)image.at<uchar>(Ydown, Xleft);
                ld = (float)image.at<uchar>(Yup, Xleft);
                ru = (float)image.at<uchar>(Ydown, Xright);
                rd = (float)image.at<uchar>(Yup, Xright);
                f1 = (Xright - _X) * lu + (_X - Xleft) * ru;
                f2 = (Xright - _X) * ld + (_X - Xleft) * rd;
                orthoMes[i * kDEMSize + j] = (Yup - _Y) * f1 + (_Y - Ydown) * f2;
            }
            else {
                return false;
            }
        }
    }
    return true;
}


float AnimationConfig::calcTk(Lmap &lmap, const float* orthoMes, int *angleMes, vector<vector<float>>& Rp, int index) {
    float Ik = 0, Rk = 0, Tk = 0;
    float4 t1, t2, normal;
    float4x4 axis;
    axis << lmap.u1, lmap.u2, lmap.u3, float4(0, 0, 0, 0);
    float4 localSunDir = frames_config[index].light_pos;
    float4 localCamDir = (frames_config[index].camera_pos - frames_config[index].camera_focal);
    localSunDir = (localSunDir.transpose() * axis).transpose();
    localCamDir = (localCamDir.transpose() * axis).transpose();
    localSunDir.normalize(), localCamDir.normalize();
    float cameraAngle, sunAngle, phaseAngle;

    for(int i = kBorderSize; i < kDEMSize - kBorderSize; ++i){
        for(int j = kBorderSize; j < kDEMSize - kBorderSize; ++j){
            Ik += orthoMes[i * kDEMSize + j];

            t1 = double4ToFloat4(lmap.ptr->DEM[(i + 1) * kDEMSize + j] - lmap.ptr->DEM[i * kDEMSize + j]);
            t2 = double4ToFloat4(lmap.ptr->DEM[i * kDEMSize + j + 1] - lmap.ptr->DEM[i * kDEMSize + j]);

            normal = t1.cross3(t2).normalized();
            cameraAngle = rad2deg(acosf(normal.dot(localCamDir)));
            sunAngle = rad2deg(acosf(normal.dot(localSunDir)));
            phaseAngle = rad2deg(acosf(localCamDir.dot(localSunDir)));
//            cout << "cameraAngle = " << cameraAngle << " " << "sunAngle = " << sunAngle << endl;
            if(cameraAngle > 0 && cameraAngle < 70 && sunAngle > 10 && sunAngle < 70){
                angleMes[i * kDEMSize + j] = 1;
                float fa = 1.0 - 0.019 * phaseAngle + 2.42e-4 * pow(phaseAngle,2) - 1.46e-6 * pow(phaseAngle,3);
                float r = fa * 2 * cos(deg2rad(sunAngle)) / (cos(deg2rad(sunAngle)) + cos(deg2rad(cameraAngle))) + (1 - fa)*cos(deg2rad(sunAngle));
                Rk += r, Rp[i][j] = r;
            }
        }
    }
    if(Rk == 0)return 0;
    Tk = Ik / Rk;
    cout << "Tk = " << Tk << endl;
    return Tk;
}


AnimationConfig::AnimationConfig(processType type) {
    this->type = type;
}


//void AnimationConfig::build(vector<float4> &landmarks)
//{
//	int size = landmarks.size();
//	float camera_distance = 5000.0f;//单位:m
//	float sun_distance = 1000.0f;//单位:m
//	for (int i = 2; i < 55; i+=3) {
//		for (int j = 2; j < 22; j+=3) {
//			//沿x轴拍摄
//			for (int k = 0; k < 10; k++) {
//				float deg = -22.5f + float(k) * 5.0f;//每次角度变化5度
//				float4 camera_focal = landmarks[j * 55 + i];
//				float4 camera_position(camera_distance*sin(Render::deg2rad(deg)) + camera_focal.x(), camera_focal.y(), camera_distance*cos(Render::deg2rad(deg)) + camera_focal.z(), 1.0f);
//				float4 camera_up(0, 1.0f, 0, 0);
//				float4 light_focal = landmarks[j * 55 + i];
//				float4 light_position(sun_distance*sin(Render::deg2rad(deg + 30.0f)) + light_focal.x(), light_focal.y(), sun_distance*cos(Render::deg2rad(deg + 30.0f)) + light_focal.z(), 1.0f);
//				float4 light_up(0, 1.0f, 0, 0);
//				frames_config.emplace_back(camera_position, camera_focal, camera_up, light_position, light_focal, light_up);
//			}
//			//沿y轴拍摄
//			for (int k = 0; k < 10; k++) {
//				float deg = -22.5f + float(k) * 5.0f;
//				float4 camera_focal = landmarks[j * 55 + i];
//				float4 camera_position(camera_focal.x(), camera_distance*sin(Render::deg2rad(deg)) + camera_focal.y(), camera_distance*cos(Render::deg2rad(deg)) + camera_focal.z(), 1.0f);
//				float4 camera_up(1.0f, 0, 0, 0);
//				float4 light_focal = landmarks[j * 55 + i];
//				float4 light_position(light_focal.x(), sun_distance*sin(Render::deg2rad(deg + 30.0f)) + light_focal.y(), sun_distance*cos(Render::deg2rad(deg + 30.0f)) + light_focal.z(), 1.0f);
//				float4 light_up(1.0f, 0, 0, 0);
//				frames_config.emplace_back(camera_position, camera_focal, camera_up, light_position, light_focal, light_up);
//			}
//		}
//	}
//}
//
//
////0度出射角，10度入射角
//void AnimationConfig::albedo_build(std::vector<float4> &landmarks) {
//	int size = landmarks.size();
//	float camera_distance = 5000.f;
//	float sun_distance = 1000.f;
//	for (int i = 0; i < size; i++) {
//		float4 camera_focal = landmarks[i];
//		float4 camera_pos(camera_focal.x(), camera_focal.y(), camera_distance, 1.0f);
//		float4 camera_up(0, 1.0f, 0, 0);
//		float4 sun_focal = landmarks[i];
//		float4 sun_pos(sun_focal.x() + sun_distance * sin(Render::deg2rad(10)), sun_focal.y(), sun_focal.z() + sun_distance * cos(Render::deg2rad(10)), 1.0f);
//		float4 sun_up(0, 1.f, 0, 0);
//		frames_config.emplace_back(camera_pos, camera_focal, camera_up, sun_pos, sun_focal, sun_up);
//	}
//}
//

