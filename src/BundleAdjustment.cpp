#include <eigen/core>
#include <opencv2/opencv.hpp>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <opencv2/core/eigen.hpp>
#include "../include/Bundleadjustment.hpp"

using namespace std;
using namespace Eigen;
using namespace cv;
using namespace Render;

struct reprojectionerror {
	reprojectionerror(const float4& _landmark, Camera *_camera, float _objectx, float _objecty, const float4& _cameraFocal, const float4& _cameraUp) :point(_landmark.template cast<double>()), camera(_camera), objectx((double)_objectx), objecty((double)_objecty), cameraFocal(_cameraFocal.template cast<double>()), cameraUp(_cameraUp.template cast<double>()) {}

	template<typename T>
	bool operator()(const T * const camPos, T * residual)const
	{
        Eigen::Matrix<T, 3, 1> oldCamPos((T)camPos[0], (T)camPos[1], (T)camPos[2]);
        Eigen::Matrix<T, 3, 1> Z3(T(cameraFocal[0]) - oldCamPos[0], T(cameraFocal[1]) - oldCamPos[1], T(cameraFocal[2]) - oldCamPos[2]);
        Z3.normalize();
        Eigen::Matrix<T, 3, 1> camUp((T)cameraUp[0], (T)cameraUp[1], (T)cameraUp[2]);
        Eigen::Matrix<T, 3, 1> X3 = Z3.cross(camUp).normalized();
        Eigen::Matrix<T, 3, 1> Y3 = Z3.cross(X3).normalized();
        Eigen::Matrix<T, 4, 1> X4(X3[0], X3[1], X3[2], T(0));
        Eigen::Matrix<T, 4, 1> Y4(Y3[0], Y3[1], Y3[2], T(0));
        Eigen::Matrix<T, 4, 1> Z4(Z3[0], Z3[1], Z3[2], T(0));
        Eigen::Matrix<T, 4, 4> MView;
        Eigen::Matrix<T, 4, 1> temp(T(0), T(0), T(0), T(0));
        MView << X4, Y4, Z4, temp;
        MView.transposeInPlace();
        Eigen::Matrix<T, 4, 1> translation(-X3.dot(oldCamPos), -Y3.dot(oldCamPos), -Z3.dot(oldCamPos), T(1.0f));
        MView.col(3) = translation;
        Eigen::Matrix<T, 4, 1> p = camera->M_viewport.template cast<T>() * camera->M_per.template cast<T>() * MView * point.template cast<T>();
        p /= p.w();
        residual[0] = p[0] - T(objectx);
        residual[1] = p[1] = T(objecty);
        return true;

//        // 旧版本
//		Eigen::Matrix<T, 4, 1> R;
//		Eigen::Matrix<T, 4, 1> p;
//		Eigen::Matrix<T, 4, 4> V;
//		Eigen::Matrix<T, 4, 1> pos((T)camera_pos[0], (T)camera_pos[1], (T)camera_pos[2], (T)1.0);
//		V = M_view.template cast<T>();
//		V.col(3) << T(0), T(0), T(0), T(0);
//		R = -V * pos;
//		V.col(3) << R[0], R[1], R[2], 1.0f;;
//		p = camera->M_viewport.template cast<T>() * camera->M_per.template cast<T>() * V * point.template cast<T>();
//		p = p / p.w();
//		residual[0] = p[0] - T(objectx);
//		residual[1] = p[1] - T(objecty);
//		return true;
	}

	static ceres::CostFunction* create(const float4& landmark, Camera *camera, float objectx, float objecty, const float4& cameraFocal, const float4& cameraUp) {
		return (new ceres::AutoDiffCostFunction<reprojectionerror, 2, 4>
			(new reprojectionerror(landmark, camera, objectx, objecty, cameraFocal, cameraUp)));
	}
	const Eigen::Vector4d point;
	const Camera *camera;
	double objectx;
	double objecty;
    const Eigen::Vector4d cameraFocal;
    const Eigen::Vector4d cameraUp;
};

////此项约束：相机参数不能偏离初始值太多
//struct errorBetweenTrueValues {
//	explicit errorBetweenTrueValues(const float4& _camInitialPos)
//		:camInitialPos(_camInitialPos.template cast<double>()) {}
//
//	template<typename T>
//	bool operator()(const T *const camera_OpPos, T *residual)const {
//		residual[0] = T(0.01) * (camera_OpPos[0] - T(camInitialPos[0])) * (camera_OpPos[0] - T(camInitialPos[0]));
//		residual[1] = T(0.01) * (camera_OpPos[1] - T(camInitialPos[1])) * (camera_OpPos[1] - T(camInitialPos[1]));
//		residual[2] = T(0.01) * (camera_OpPos[2] - T(camInitialPos[2])) * (camera_OpPos[2] - T(camInitialPos[2]));
//		return true;
//	}
//
//	static ceres::CostFunction* create(const float4& camInitialPos) {
//		return (new ceres::AutoDiffCostFunction<errorBetweenTrueValues, 3, 4>
//			(new errorBetweenTrueValues(camInitialPos)));
//	}
//	Eigen::Vector4d camInitialPos;//相机初始位置
//};


tuple<bool, float, float> BundleAdjustment::calcMatching(const Mat& image, const Mat &orthoImage) {
	Mat srcImage;
	image.copyTo(srcImage);

	// 计算出 renderMap 的中心在 image (原始拍摄图像)上的像素坐标，作为 ObjectX、ObjectY
	Mat resultImage;
	int resImageCols = image.cols - orthoImage.cols + 1;
	int resImageRows = image.rows - orthoImage.rows + 1;
	resultImage.create(resImageCols, resImageRows, CV_8UC3);
	matchTemplate(image, orthoImage, resultImage, TM_CCOEFF_NORMED);
	normalize(resultImage, resultImage, 0, 1, NORM_MINMAX, -1, Mat());
	double minVal, maxVal;

	Point minLocation, maxLocation, matchLocation;
	minMaxLoc(resultImage, &minVal, &maxVal, &minLocation, &maxLocation);
    if(maxVal < 0.65)return std::make_tuple(false, 0, 0);
    matchLocation = maxLocation;
    return std::make_tuple(true, matchLocation.x + orthoImage.cols / 2, matchLocation.y + orthoImage.rows / 2);
}


// index 表示相机序号
void BundleAdjustment::OptSingleCamParam(const vector<cv::Mat>& coResImage, const vector<cv::Mat>& newOrthoImage, AnimationConfig* animationConfig, const vector<float4>& landmarks, int index, Camera* camera)
{
    ceres::Problem problem;
    double camPos[] {animationConfig->frames_config[index].camera_pos[0], animationConfig->frames_config[index].camera_pos[1], animationConfig->frames_config[index].camera_pos[2], 1.0f};
    // i 表示该相机的 coRes 的第几个
    for(int i = 0; i < coResImage.size(); ++i){
        // tmp 表示具体哪个 landmark
        int tmp = 0;
        if(animationConfig->type == ALBEDO)tmp = animationConfig->coRes[index][i];
        else tmp = animationConfig->coRes[index + landmarks.size()][i];
        //todo: 比对coResImage 和 newOrthoImage 中的差异，得到 landmark 的 x和 y 坐标
        auto tuple = calcMatching(coResImage[i], newOrthoImage[i]);
        if(std::get<0>(tuple)){
            ceres::CostFunction* costFunction = reprojectionerror::create(landmarks[tmp], camera, std::get<1>(tuple), std::get<2>(tuple), animationConfig->frames_config[index].camera_focal, animationConfig->frames_config[index].camera_up);
            problem.AddResidualBlock(costFunction, nullptr, camPos);
        }
    }


//	float objectX, objectY;
//	ceres::Problem problem;
//	double *camparamnoise = {animationConfig->frames_config->};
//	camera->set_look_at(noiseConfig.camera_pos, noiseConfig.camera_focal, noiseConfig.camera_up);
//	for (int i = 0; i < 11; i++) {
//		float4 result = camera->P * landmarks[i];
//		result /= result.w();
//		float _X = result.x();
//		float _Y = result.y();
//		if ((_X >= kDEMSize / 2) && (_X < kImageSize - kDEMSize * 0.5) && (_Y >= kDEMSize / 2) && (_Y < kImageSize - 0.5 * kDEMSize)) {
//			string orthoImagePath = prepath + to_string(i) + "\\Ortho_image" + to_string(index) + ".png";//对应该相机的正射校正图像
//			if (_access(orthoImagePath.c_str(), 0) == -1) {
//				cout << "can't find" << endl;
//				continue;
//			}
//			calcMatching(orthoImagePath, image, objectX, objectY);
//			cout << "objectX = " << objectX << "," << "objectY = " << objectY << endl;
//			ceres::CostFunction* cost_func = reprojectionerror::create(landmarks[i], camera, objectX, objectY, noiseConfig, camera->M_view);
//			problem.AddResidualBlock(cost_func, nullptr, camparamnoise);
//		}
//	}
	/*ceres::CostFunction* cost_func2 = errorBetweenTrueValues::create(noiseConfig.camera_pos);
	problem.AddResidualBlock(cost_func2, NULL, camparamnoise);*/

	ceres::Solver::Options solver_options;
	solver_options.linear_solver_type = ceres::DENSE_QR;
	solver_options.minimizer_progress_to_stdout = true;
	solver_options.max_num_iterations = 100;
	ceres::Solver::Summary summary;

	ceres::Solve(solver_options, &problem, &summary);
	cout << "summary:\n" << summary.BriefReport() << endl;
    cout << animationConfig->frames_config[index].camera_pos << endl;
	animationConfig->frames_config[index].camera_pos << (float)camPos[0], (float)camPos[1], (float)camPos[2], (float)1.0f;
    cout << animationConfig->frames_config[index].camera_pos << endl;
}
