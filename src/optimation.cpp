#include "../include/optimation.hpp"

using namespace std;
using namespace Eigen;
using namespace cv;
using namespace Render;
using namespace ceres;

OptimizeTool::OptimizeTool() = default;
OptimizeTool::~OptimizeTool() = default;

struct CostFunctor1 {
    const Eigen::Vector4d DEMDown;
    const Eigen::Vector4d DEMRight;
    const double localX;
    const double localY;
	double orthoValue;
	double Tk;
	const Eigen::Vector4d sunDir;
	const Eigen::Vector4d cameraDir;
	double A;

	CostFunctor1(const double4& _DEMDown, const double4& _DEMRight, const double& _localX, const double _localY, const float& _orthoValue, const float& _Tk, const float4& _sunDir, const float4& _cameraDir, const double& _A) : DEMDown(_DEMDown), DEMRight(_DEMRight), localX(_localX), localY(_localY), orthoValue((double)_orthoValue), Tk((double)_Tk), sunDir(_sunDir.cast<double>()), cameraDir(_cameraDir.cast<double>()), A(_A) {}

	template <typename T>
	bool operator()(const T * const DEM, T * residual)const
	{
        Eigen::Matrix<T, 4, 1> normal;
        Eigen::Matrix<T, 3, 1> t1, t2, normal3f;
        t1 << T(DEMDown[0] - localX), T(DEMDown[1] -localY), T(DEMDown[2] - DEM[0]);
        t2 << T(DEMRight[0] - localX), T(DEMRight[1] - localY), T(DEMRight[2] - DEM[0]);
        normal3f = t1.cross(t2).template cast<T>();
        normal3f = normal3f.normalized();
        normal << normal3f[0], normal3f[1], normal3f[2], T(0);
        T camera_angle = acos(normal.dot(cameraDir)) / T(kPi) * T(180.0f);
        T light_angle = acos(normal.dot(sunDir.cast<T>())) / T(kPi) * T(180.0f);
        T phase_angle = acos(cameraDir.dot(sunDir)) / T(kPi) * T(180.0f);
        T fa = T(1) - T(0.019) * phase_angle + T(2.42e-4) * pow(phase_angle, T(2)) - T(1.46e-6) * pow(phase_angle, T(3));
        T Rk = (fa * T(2.0f) * cos(light_angle * T(kPi) / T(180)))/ (cos(light_angle * T(kPi) / T(180.0f)) + cos(camera_angle * T(kPi) / T(180.0f))) + (T(1.0f) - fa) * cos(light_angle * T(kPi) / T(180.0f));
        residual[0] = (orthoValue - Tk * A * Rk) * (orthoValue - Tk * A * Rk);
        return true;
	}

	static ceres::CostFunction* Create(const double4& DEMDown, const double4& DEMRight, const double& localX, const double& localY, const float& orthoValue, const float& Tk, const float4& sunDir, const float4& cameraDir, const double& A)
	{
		return (new ceres::AutoDiffCostFunction<CostFunctor1, 1, 1>
			(new CostFunctor1(DEMDown, DEMRight, localX, localY, orthoValue, Tk, sunDir, cameraDir, A)));
	}
};

// 第二项残差：局部高度二阶偏导数平方和
struct CostFunctor2 {
    const Eigen::Vector4d DEMLeft;
    const Eigen::Vector4d DEMRight;
    const Eigen::Vector4d DEMUp;
    const Eigen::Vector4d DEMDown;
	explicit CostFunctor2(const double4& _DEMLeft, const double4& _DEMRight, const double4& _DEMUp, const double4& _DEMDown) :DEMLeft(_DEMLeft), DEMRight(_DEMRight), DEMUp(_DEMUp), DEMDown(_DEMDown) {}

	template<typename T>
	bool operator()(const T* const DEM, T* residual)const
	{
		T d1 = (T(2.0) * DEM[0] - DEMRight[2] - DEMLeft[2]) / T(kModelGSD * kModelGSD);
		T d2 = (T(2.0) * DEM[0] - DEMUp[2] - DEMDown[2]) / T(kModelGSD * kModelGSD);
		residual[0] = T(0.6 * 0.6) * (d1 * d1 + d2 * d2);
//        residual[0] = (d1 * d1 + d2 * d2);
		return true;
	}

    static ceres::CostFunction *Create(const double4& DEMLeft, const double4& DEMRight, const double4& DEMUp, const double4& DEMDown){
        return (new ceres::AutoDiffCostFunction<CostFunctor2, 1, 1>(new CostFunctor2(DEMLeft, DEMRight, DEMUp, DEMDown)));
    }
};

// 第三项残差，保证优化后的地形和初始地形差别不大
struct CostFunctor3 {
	explicit CostFunctor3(const double& _init_DEM) :init_DEM(_init_DEM) {}

	template<typename T>
	bool operator()(const T* const DEM, T* residual)const
	{
		residual[0] = T(0.0008 * 0.0008) * (DEM[0] - init_DEM) * (DEM[0] - init_DEM);
		return true;
	}
	static ceres::CostFunction* Create(const double& init_DEM)
	{
		return (new ceres::AutoDiffCostFunction<CostFunctor3, 1, 1>(new CostFunctor3(init_DEM)));
	}
    double init_DEM;
};

void OptimizeTool::optiSlope(Lmap &lmap, const vector<AnimationConfig*>& configVec, const string& path, int index) {
    ceres::Problem problem;
    float4x4 axis;
    axis << lmap.u1, lmap.u2, lmap.u3, float4(0, 0, 0, 0);

    for(int ele = 1; ele < configVec.size(); ++ele){
        vector<cv::Mat> coResImage;
        cv::String pattern;
        if(!configVec[ele]->coRes[index].empty()){
            for (int cosImage : configVec[ele]->coRes[index]) {
                pattern = getPathFromEnum(configVec[ele]->type) + formatNum(cosImage) + ".png";
                iodata::readImages(pattern, coResImage);
            }
            cout << "Landmark " << to_string(index) << " read image done, coResImage size is " << coResImage.size() << endl;
        }
        else{
            cout << "Landmark " << to_string(index) << " didn't appear in any images" << endl;
            continue;
        }

        for(int k = 0; k < configVec[ele]->coRes[index].size(); ++k){
            int tmp = configVec[ele]->coRes[index][k]; // 相机位置的序号
//            if(ele == FB1 && tmp == 450)continue;
            if(ele == FB5A)continue;
            if(ele == FB6B)continue;
            if(coResImage[k].cols != kImageSize || coResImage[k].rows != kImageSize){
                cout << "error image" << endl;
                continue;
            }
            auto *orthoDEM = new float[kDEMSize * kDEMSize];
            auto *angleMes = new int[kDEMSize * kDEMSize];
            for(int i = 0; i < kDEMSize * kDEMSize; ++i){
                angleMes[i] = 0;
            }

            bool flag = configVec[ele]->orthoRectify(lmap, orthoDEM, coResImage[k], tmp);
            if(flag){
                vector<vector<float>> Rp(kDEMSize, vector<float>(kDEMSize, 0));
                iodata::writeImage(orthoDEM, path + "lmap_" + to_string(index) + "/" + getTextFromEnum(configVec[ele]->type) + formatNum(tmp) + ".png");

                float Tk = configVec[ele]->calcTk(lmap, orthoDEM, angleMes, Rp, tmp);
                if(Tk == 0)continue;

                float4 localSunDir = configVec[ele]->frames_config[tmp].light_pos;
                float4 localCamDir = configVec[ele]->frames_config[tmp].camera_pos - configVec[ele]->frames_config[tmp].camera_focal;
                localSunDir = (localSunDir.transpose() * axis).transpose();
                localCamDir = (localCamDir.transpose() * axis).transpose();
                localSunDir.normalize(), localCamDir.normalize();

                for(int i = kBorderSize; i < kDEMSize - kBorderSize; ++i){
                    for(int j = kBorderSize; j < kDEMSize - kBorderSize; ++j){
                        if (orthoDEM[i * kDEMSize + j] != 0 && angleMes[i * kDEMSize + j] == 1) {
                            ceres::CostFunction *costFunc1 = CostFunctor1::Create(lmap.ptr->DEM[(i + 1) * kDEMSize + j], lmap.ptr->DEM[i * kDEMSize + j + 1], lmap.ptr->DEM[i * kDEMSize + j][0], lmap.ptr->DEM[i * kDEMSize + j][1], orthoDEM[i * kDEMSize + j], Tk, localSunDir, localCamDir, lmap.ptr->albedoMap[i * kDEMSize + j]);
                            problem.AddResidualBlock(costFunc1, nullptr, &lmap.ptr->DEM[i * kDEMSize + j][2]);
                        }
                    }
                }
            }
            delete[] orthoDEM, delete[] angleMes;
        }
    }

    //添加第二项残差:局部高度二阶偏导数平方和
    for (int i = kBorderSize; i < kDEMSize - kBorderSize; i++)
    {
        for (int j = kBorderSize; j < kDEMSize - kBorderSize; j++)
        {
            ceres::CostFunction * costFunc2 = CostFunctor2::Create(lmap.ptr->DEM[i * kDEMSize + j - 1], lmap.ptr->DEM[i * kDEMSize + j + 1],lmap.ptr->DEM[(i - 1) * kDEMSize + j],lmap.ptr->DEM[(i + 1) * kDEMSize + j]);
            problem.AddResidualBlock(costFunc2, nullptr, &lmap.ptr->DEM[i * kDEMSize + j][2]);
        }
    }

    //添加第三项残差：使得优化后的DEM和初值不能相差过大
    for (int i = kBorderSize; i < kDEMSize - kBorderSize; i++)
    {
        for (int j = kBorderSize; j < kDEMSize - kBorderSize; j++)
        {
            ceres::CostFunction *costFunc3 = CostFunctor3::Create(lmap.ptr->DEM[i * kDEMSize + j][2]);
            problem.AddResidualBlock(costFunc3, nullptr, &lmap.ptr->DEM[i * kDEMSize + j][2]);
        }
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::CGNR;
    options.num_threads = 12;
    options.function_tolerance = 0.001;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 100;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout << summary.BriefReport() << endl;
    updateHeight(lmap);
}


void OptimizeTool::updateHeight(Lmap &lmap) {
	/*
		当每次迭代停止计算出DEM之后，使用该方法计算高度
		DEM就是得到的高度
		size代表DEM尺寸
		w1表示DEM每个点的四周高度权重取值
		w2表示该点本身高度所占权重
	*/
//    float w1 = 1.0f, w2 = 0.1;
//	for (int i = 1; i < kDEMSize - 1; i++) {
//			for (int j = 1; j < kDEMSize - 1; j++) {
//			lmap.ptr->DEM[i][j] = w1 * (lmap.ptr->DEM[i - 1][j] + lmap.ptr->DEM[i + 1][j] + lmap.ptr->DEM[i][j - 1] + lmap.ptr->DEM[i][j + 1]) + w2 * lmap.ptr->DEM[i][j];
//            lmap.ptr->DEM[i][j] = lmap.ptr->DEM[i][j] / (4 * w1 + w2);
//		}
//	}

	//应该将DEM中心的高度转为0
	double mid = lmap.ptr->DEM[kDEMSize / 2 * kDEMSize + kDEMSize / 2][2];
	int size = kDEMSize * kDEMSize;
	for (int i = 0; i < size; i++) {
        lmap.ptr->DEM[i][2] -= mid;
	}
    // 这里实际应该把这个向量转为本体坐标系下再加到landmark向量上，但是由于只改变高度，且改变的方向和landmark的向量方向一致，所以直接加上去就可以
	lmap.landmark(2) += (float)mid;
}

float OptimizeTool::calcAveAbsDiff(Lmap &lmap, string &preDEMPath) {
	float change = 0;
	Lmap preLmap;
	preLmap.readLmap(preDEMPath);
	int size = kDEMSize * kDEMSize;
	for (int i = 0; i < size; i++) {
		change += fabs((float)lmap.ptr->DEM[i][2] - (float)preLmap.ptr->DEM[i][2]);
	}
	return change / (float)size;
}
