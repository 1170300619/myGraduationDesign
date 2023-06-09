#pragma warning(disable:4996);
#include "../include/Renderer.hpp"
#include "../include/Material.hpp"
#include "../include/Mipmap.hpp"
#include "opencv2/imgproc.hpp"
#include <opencv2/opencv.hpp>
#include "../include/BundleAdjustment.hpp"
#include "../include/Optimation.hpp"
#include "../include/SphericalSurfICQ.hpp"

using namespace std;
using namespace Eigen;
using namespace Render;

int num_frag = 0;

void Renderer::read_data(const string &modelFilePath, const string &config_file)
{
    cout << "read data... " << endl;
    s = iodata::load_config(config_file);
    iodata::readPly(modelFilePath, s->model);

#if 0
	vector<float4> randomPoints;
	vector<float4> initLandmarks;
    fibonacciGrid(randomPoints, 500);
//	iodata::generateWellDistributedPoints(randomPoints, 500);
	iodata::chooseLandmarksFromSpherome(s->model, randomPoints, initLandmarks);
#endif

	cout << "Read data done!" << endl;
}


void Renderer::render()
{
    iodata::readLandmarks(kLandmarksPath + "Landmarks_0.txt", landmarks);

#if 0 // 查找landmark最近点和均值分布
//    ofstream out;
//    string nearestPath = "./data/feibonaqi_landmark_nearest_4.txt";
//    out.open(nearestPath.c_str());
//    iodata::readLandmarks(kLandmarksPath + "Landmarks_0.txt", landmarks);
//    float maxDis = INT_MIN, minDis = INT_MAX;
//    int oneGSD = 0, threeGSD = 0;
//    int aveOne = 0, aveThree = 0;
//    for(int i = 0; i < landmarks.size(); ++i){
//        vector<float> disVec;
//        for(int j = 0; j < landmarks.size(); ++j){
//            if(j != i){
//                float tmp = Dist(landmarks[i], landmarks[j]);
//                maxDis = max(maxDis, tmp);
//                minDis = min(minDis, tmp);
//                disVec.emplace_back(tmp);
//            }
//        }
//        sort(disVec.begin(), disVec.end());
//        out << disVec[0] << " " << disVec[1] << " " << disVec[2] << " " << disVec[3] << endl;
//        if((disVec[0] + disVec[1] + disVec[2] + disVec[3]) / 4 <= 2.0f)aveOne++;
//        if((disVec[0] + disVec[1] + disVec[2] + disVec[3]) / 4 <= 4.0f)aveThree++;
//        if(disVec[0] < 2.0f)oneGSD++;
//        if(disVec[0] < 4.0f)threeGSD++;
//    }
//    cout << oneGSD << " " << threeGSD << endl;
//    cout << aveOne << " " << aveThree << endl;
#endif

#if 0 // 计算因为图像中暗像素被筛选掉的图像个数
    for(int i = 0; i < 990; ++i){
        int blackNum = 0;
        cv::Mat indexImage = cv::imread(kFB6BPath + formatNum(i) + ".png", -1);
        auto row = indexImage.rows, col = indexImage.cols;
        for(int j = 0; j < row; ++j){
            for(int k = 0; k < col; ++k){
                if(indexImage.at<uchar>(j, k) == 0)blackNum++;
            }
        }
        if(((float)blackNum / (float)row / (float)col) >= 0.3)cout << i << endl;
    }
#endif

    for(processType i = ALBEDO; i <= FB6B; i = processType(i + 1)){
        configVec.emplace_back(new AnimationConfig(i));
    }

#if 0
    // 先读取真实的相机参数，然后加噪声后保存，之后使用的都是这些带误差的相机参数
    for(auto & i : configVec){
        iodata::readCamAndSunMes(i, landmarks, getPathFromEnum(i->type) + '/', s);
    }
    iodata::addNoiseToCamPos(configVec);
#endif

    for(auto & i : configVec){
        iodata::readCamAndSunMes(i, landmarks, kNoiseCamAndSunMes + getTextFromEnum(i->type) + '/', s);
    }
    cout << "Read camera and sun meassage done" << endl;

#if 0
    // 计算对应关系并保存
    for(auto & i : configVec){
        calcCoRES(i);
    }
#endif

    // i 从 1 开始，ALBEDO 类型不需要计算关联
    for(auto & i : configVec){
        iodata::readCorrespnse("./data/coRes/" + getTextFromEnum(i->type) + "Res.txt", i->coRes);
    }

#if 0 // 计算每个landmark对应的图像数量
    for(int i = 1; i < configVec.size(); ++i){
        vector<int> picNum(6, 0);
        for(int j = 0; j < landmarks.size(); ++j){
            if(configVec[i]->coRes[j].size() > 0 && configVec[i]->coRes[j].size() <= 20)picNum[0]++;
            else if(configVec[i]->coRes[j].size() > 20 && configVec[i]->coRes[j].size() <= 40)picNum[1]++;
            else if(configVec[i]->coRes[j].size() > 40 && configVec[i]->coRes[j].size() <= 60)picNum[2]++;
            else if(configVec[i]->coRes[j].size() > 60 && configVec[i]->coRes[j].size() <= 80)picNum[3]++;
            else if(configVec[i]->coRes[j].size() > 80 && configVec[i]->coRes[j].size() <= 100)picNum[4]++;
            else if(configVec[i]->coRes[j].size() > 100)picNum[5]++;
        }
        for(auto & ele : picNum)cout << ele << endl;
        cout << "----------------------------" << endl;
    }
#endif

#if 0
    vector<int> pngNum(5, 0);
    for(int i = 0; i < landmarks.size(); ++i) {
        string path = kFirstOptimationPath + "lmap_" + to_string(i) + "/*.png";
        getPngFiles(path, pngNum);
    }
    for(auto & item : pngNum)cout << item << endl;
#endif

#if 1
    for (int seq = 47; seq < 48; ++seq) {
        string fileName = kFirstOptimationPath + "lmap_" + to_string(seq) + "/v0.txt";
        auto lmap = Lmap(landmarks[seq]);
        if(!lmap.GenIniDEMFromBennu(s->model, kModelGSD)){
            cout << "The Lmap " + to_string(seq) + " genInitDEM failed" << endl;
            continue;
        }

        cv::Mat rgbObAlbedoImage = cv::imread(kAlbedoPath + formatNum(seq) + ".png");
        cv::Mat grayObAlbedoImage;
        cv::cvtColor(rgbObAlbedoImage, grayObAlbedoImage, cv::COLOR_RGB2GRAY);
        if(grayObAlbedoImage.cols != kImageSize || grayObAlbedoImage.rows != kImageSize){
            cout << "albedo image error" << endl;
            continue;
        }
        cout << "imread " + to_string(seq) + " albedo image done" << endl;
        auto *albedoMes = new float[kDEMSize * kDEMSize];
        auto *angleMes = new int[kDEMSize * kDEMSize];
        for(int i = 0; i < kDEMSize * kDEMSize; ++i){
            angleMes[i] = 0;
        }
        bool flag = configVec[0]->orthoRectify(lmap, albedoMes, grayObAlbedoImage, seq);
        cout << "albedo orthoRectify done" << endl;
        // todo: 怎么判断一张正射校正后的反照率图像能不能使用
        if(flag){
            lmap.writeLmap(fileName); // 先写入 Lmap 数据，这个函数内有路径检查的函数，所以可以创建路径，下面的反照率正射图像就可以写入了
            float maxBright = 0, minBright = 255.0f;
            int size = kDEMSize * kDEMSize;
            for(auto i = 0; i < size; ++i){
                minBright = min(minBright, albedoMes[i]);
                maxBright = max(maxBright, albedoMes[i]);
            }
            float conver = (kMaxAlbedo - kMinAlbedo) / (maxBright - minBright);
            for (int i = 0; i < size; i++) {
                lmap.ptr->albedoMap[i] = kMinAlbedo + (albedoMes[i] - minBright) * conver;
            }
            iodata::writeImage(albedoMes, kFirstOptimationPath + "lmap_" + to_string(seq) + "/albedo.png");
        }
        else{
            cout << "this albedo image doesn't meet conditions" << endl;
            delete[] albedoMes;
            delete[] angleMes;
            continue;
        }
        delete[] albedoMes;
        delete[] angleMes;
        int step = 0;
        float aveAbsDiff = 1.0f;
        string preLmapPath, lmapPath;
        while (aveAbsDiff > 0.05) {
            preLmapPath = kFirstOptimationPath + "lmap_" + to_string(seq) + "/v" + to_string(step) + ".txt";
            lmapPath = kFirstOptimationPath + "lmap_" + to_string(seq) + "/v" + to_string(step + 1) + ".txt";
            if(step > 0)lmap.readLmap(preLmapPath);
            cout << "Optimizing hight......" << endl;
            OptimizeTool::optiSlope(lmap, configVec, kFirstOptimationPath, seq);
            lmap.writeLmap(lmapPath);
            aveAbsDiff = OptimizeTool::calcAveAbsDiff(lmap, preLmapPath);
            cout << "Iteration = " << aveAbsDiff << endl;
            step++;
        }

        // 保存最新的 DEM 的值，更新 landmark 容器内的值
        lmap.writeLmap(kFirstOptimationPath + "lmap_" + to_string(seq) + "/" + "finalDEM.txt");
        landmarks[seq] = lmap.landmark;

#if 1
        // 测试：将初始DEM保存为模型，方便比对前后的差别
        auto *preModel = new Model();
        iodata::readDEMWithoutEdge(fileName, preModel); // 即v0.txt
        iodata::write_ply(preModel, kFirstOptimationPath + "lmap_" + to_string(seq) + '/' + "preModel.ply");
        cout << "Write preModel done, the model path is " + kFirstOptimationPath + "lmap_" + to_string(seq) + '/' + "preModel.ply" << endl;
#endif

        auto *procModel = new Model();
        iodata::readDEMWithoutEdge(lmapPath, procModel); // 即最后保存的一个txt
        iodata::write_ply(procModel, kFirstOptimationPath + "lmap_" + to_string(seq) + '/' + "procModel.ply");
        cout << "Write newModel done, the model path is " + kFirstOptimationPath + "lmap_" + to_string(seq) + '/' + "procModel.ply" << endl;
        delete procModel;

    }
#endif

    cout << "End of all lmap optimization" << endl;

#if 0
    // 将第一次更新的 landmark 保存一下，方式是读取每一个优化后的 finalDEM.txt文件，如果没有改路径的则不更新
    for(int i = 0; i < landmarks.size(); ++i){
        ifstream in;
        string opLandmarkPath = kFirstOptimationPath + "lmap_" + to_string(i) + "/finalDEM.txt";
        if (_access(opLandmarkPath.c_str(), 0) == -1){
            continue;
        }
        in.open(opLandmarkPath.c_str());
        int demSize;
        in >> demSize >> demSize;
        in >> landmarks[i][0] >> landmarks[i][1] >> landmarks[i][2] >> landmarks[i][3];
        in.close();
    }

    iodata::writeLandmarks(landmarks, kLandmarksPath, 1);
    cout << "Save landmarks after first optimazion done" << endl;
#endif

#if 0 // BA阶段
    iodata::readLandmarks(kLandmarksPath + "Landmarks_1.txt", landmarks);

    for(auto & i : configVec){
        int interval = 0;
        if(i->type != ALBEDO)interval = landmarks.size();
        for(int j = 0; j < i->frames_config.size(); ++j){ // j 表示第几个相机
            vector<cv::Mat> coResImage;
            vector<cv::Mat> newOrthoImage;
            if(i->coRes[j + interval].empty())continue;
            for(auto & ele : i->coRes[j + interval]){ // ele 代表第几个 landmark 点
                cout << "ele = " << ele << endl;
                string pattern = getPathFromEnum(i->type) + formatNum(ele) + ".png";
                if (_access(pattern.c_str(), 0) == -1){
                    cout << "can not find " + to_string(j) + " image path" << endl;
                    continue;
                }
                iodata::readImages(pattern, coResImage); // 将原始观测图像读入，作为匹配的原图像

                Lmap lmap;
                string finalLmapPath = kFirstOptimationPath + "lmap_" + to_string(ele) + "/finalDEM.txt";
                if(_access(finalLmapPath.c_str(), 0) == -1)continue;
                iodata::readLmap(finalLmapPath, lmap);

                auto *orthoMes = new float[kDEMSize * kDEMSize];
                auto *angleMes = new int[kDEMSize * kDEMSize];
                for(int k = 0; k < kDEMSize * kDEMSize; ++k){
                    orthoMes[k] = 0;
                    angleMes[k] = 0;
                }

                bool flag = i->orthoRectify(lmap, orthoMes, coResImage[coResImage.size() - 1], j);
                if(!flag)continue;
                string newOrthoImagePath = kBundleAdjustmentPath + getTextFromEnum(i->type) + "/camera_" + to_string(j) + "/landmark_" + to_string(ele) + ".png";
                if(!savePathCheack(newOrthoImagePath)){
                    cout << "Failed to open newOrthoImage path" << endl;
                }
                iodata::writeImage(orthoMes, newOrthoImagePath);
//                float Tk = i->calcTk(lmap, orthoMes, angleMes, Rp, j);
//
//                cv::Mat tmp = cv::Mat(kDEMSize, kDEMSize, CV_8UC1, cv::Scalar::all(0));
//                float maxBright = 0, minBright = 255.0f;
//                for(int k = 0; k < kDEMSize; ++k){
//                    for(int t = 0; t < kDEMSize; ++t){
//                        float value = Rp[k][t] * Tk; // 这里删除了(float)lmap.ptr->albedoMap[k * kDEMSize + t]
//                        tmp.at<uchar>(k, t) = value;
//                        minBright = min(minBright, value);
//                        maxBright = max(maxBright, value);
//                    }
//                }
//                float conver = 255.0f / (maxBright - minBright);
//                for(int k = 0; k < kDEMSize; ++k){
//                    for(int t = 0; t < kDEMSize; ++t){
//                        auto value = tmp.at<uchar>(k ,t);
//                        tmp.at<uchar>(k, t) = ((float)value - minBright) * conver;
//                    }
//                }
//                string simulateImagePath = kBundleAdjustmentPath + getTextFromEnum(i->type) + "/camera_" + to_string(j) + "/landmark_" +
//                        to_string(ele) + ".png";
//                if(!savePathCheack(simulateImagePath)){
//                    cout << "Failed to open newOrthoImage path" << endl;
//                }
//                auto writeFlag = cv::imwrite(simulateImagePath, tmp);
//                if(!writeFlag)continue;
                iodata::readImages(newOrthoImagePath, newOrthoImage);
                delete[] orthoMes, delete[] angleMes;
            }
            BundleAdjustment::OptSingleCamParam(coResImage, newOrthoImage, i, landmarks, j, s->camera);
        }

        string bundleAdjustmentCamPosPath = kNoiseCamAndSunMes + getTextFromEnum(i->type) + "/bundleAdjustmentCamPos.txt";
        ofstream out;
        out.open(bundleAdjustmentCamPosPath.c_str());
        for(int j = 0; j < i->frames_config.size(); ++j){
            out << i->frames_config[j].camera_pos[0] << " " << i->frames_config[j].camera_pos[1] << " " << i->frames_config[j].camera_pos[2] << " " << i->frames_config[j].camera_pos[3];
            if(j < i->frames_config.size() - 1)out << endl;
        }
    }
#endif


#if 0  // 生成最终模型
    SphericalSurfICQ surfIcq;
    surfIcq.initFlat(4000000, landmarks.size());
    surfIcq.updataLandmarkVp();
    iodata::write3DPoints(surfIcq.vertexList, "C:/SPC_Files/SPC_generate_file/all_Maplets_model/400_all_Maplets_model.ply");
#endif

}

void Renderer::calcCoRES(AnimationConfig* animationConfig) {
    auto n = landmarks.size();
    auto m = animationConfig->frames_config.size();
    vector<vector<int>> existMat(n, vector<int>(m, 0));

    if(animationConfig->type == ALBEDO){
        for(int i = 0; i < landmarks.size(); ++i){
            existMat[i][i] = 1;
        }
    }
    else{
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                float4 result = animationConfig->frames_config[j].P * landmarks[i];
                float4 sunDir = animationConfig->frames_config[j].light_pos.normalized();
                float4 camDir = (animationConfig->frames_config[j].camera_pos - animationConfig->frames_config[j].camera_focal).normalized();
                float4 normal = landmarks[i].normalized();
                normal.w() = 0; // 这里很重要，计算角度时，如果使用的是 float4，一定要把 w 分量置为 0
                float cameraAngle = rad2deg(acosf(normal.dot(camDir)));
                float sunAngle = rad2deg(acosf(normal.dot(sunDir)));
                if(cameraAngle <= 0 || cameraAngle >= 70 || sunAngle <= 10 || sunAngle >= 70)continue;
                result /= result.w();
                float X = result.x();
                float Y = result.y();
                if ((X >= kDEMSize * 0.5) && (X < kImageSize - kDEMSize * 0.5) && (Y >= kDEMSize * 0.5) && (Y < kImageSize - 0.5 * kDEMSize)) {
                    existMat[i][j] = 1;
                }
            }
        }
    }

    cout << "Calc " << getTextFromEnum(animationConfig->type) << " res done" << endl;

	//以下根据存在矩阵将图像与landmark对应关系写入文件
	ofstream out;
    string resPath = "./data/coRes/" + getTextFromEnum(animationConfig->type) + "Res.txt";
    out.open(resPath.c_str());
    if(!out.is_open()){
        cerr << "Fail to open " + resPath << endl;
        throw exception();
    }

	//每个landmark对应的图像有哪些
	for (int i = 0; i < n; i++) {
        out << i;
		for (int j = 0; j < m; j++) {
			if (existMat[i][j] == 1)
				out << " " << j;
		}
		out << "\n";
	}
	//每个图像中包含哪些landmarks
	for (int j = 0; j < m; j++) {
		out << j;
		for (int i = 0; i < n; i++) {
			if (existMat[i][j] == 1)
                out << " " << i;
		}
		if(j < m - 1)out << "\n";
	}
	out.close();
    cout << "Write correRes file done" << endl;
}

//void Renderer::view_transform()
//{
//    //    cout << "view transform..." << endl;
//
//    for (auto mesh: s->model->meshes)
//    {
//        //#pragma omp parallel for
//        for (int i = 0; i < mesh->num_triangle; i++)
//        {
//            mesh->triangles[i].normal = s->camera->Q * mesh->triangles[i].normal; // normal vector transform
//        }
//    }
//}

inline bool Renderer::clipping(const VertexP &v_0, const VertexP &v_1, const VertexP &v_2)
{
    //    if (v_0.clip.w() < 1e-5 || v_1.clip.w() < 1e-5 || v_2.clip.w() < 1e-5) // to avoid divide by 0
    //    {
    //        return true;
    //    }
    if (v_0.clip.x() < -v_0.clip.w() && v_1.clip.x() < -v_1.clip.w() && v_2.clip.x() < -v_2.clip.w())
    {
        return true;
    }
    if (v_0.clip.x() > v_0.clip.w() && v_1.clip.x() > v_1.clip.w() && v_2.clip.x() > v_2.clip.w())
    {
        return true;
    }
    if (v_0.clip.y() < -v_0.clip.w() && v_1.clip.y() < -v_1.clip.w() && v_2.clip.y() < -v_2.clip.w())
    {
        return true;
    }
    if (v_0.clip.y() > v_0.clip.w() && v_1.clip.y() > v_1.clip.w() && v_2.clip.y() > v_2.clip.w())
    {
        return true;
    }
    if (v_0.clip.z() < 0 || v_1.clip.z() < 0 || v_2.clip.z() < 0)
    {
        return true;
    }
    if (v_0.clip.z() > v_0.clip.w() && v_1.clip.z() > v_1.clip.w() && v_2.clip.z() > v_2.clip.w())
    {
        return true;
    }
    return false;
}

vector<VertexP> Renderer::clip_near(const VertexP &v_0, const VertexP &v_1, const VertexP &v_2)
{
    VertexP vertices[3] = {v_0, v_1, v_2};
    vector<VertexP> output;
    for (int i = 0; i < 3; i++)
    {
        const VertexP &start = vertices[i];
        const VertexP &end = vertices[(i + 1) % 3];

        if (end.clip.z() > 0)
        {
            if (start.clip.z() < 0)
            {
                float a = start.clip.z();
                float b = end.clip.z();
                float t = b / (b - a);
                output.push_back(Render::lerp(t, end, start));
            }
            output.push_back(end);
        } else if (start.clip.z() > 0)
        {
            float a = end.clip.z();
            float b = start.clip.z();
            float t = b / (b - a);
            output.push_back(Render::lerp(t, start, end));
        }
    }
    return output;
}

//void Renderer::draw(int shader_mode)
//{
//    for (auto mesh: s->model->meshes)
//    {
//        Uniform u = mesh->material->get_uniform();
//		if (shader_mode == 0)
//		{
//			s->shader = mesh->material->get_shader(); // imageShader
//		}
//		else if(shader_mode == 1)
//		{
//			s->shader = ShadowShader::instance();//shadowShader
//		}
//        s->shader->set_uniform(&u); // set uniform
//        s->shader->set_light(&s->light); // set light source
//        s->shader->set_camera(s->camera); // set camera
//        s->shader->set_texture_type(s->texture_type); // set texture type
//        s->shader->set_sampler(s->sampler); // set sampler
//
//        //#pragma omp parallel for
//        //        for (int i = 0; i < mesh->num_vertex; i++)
//        //        {
//        //            s->shader->vertex_shader(mesh->vertices[i]);
//        //            perspective_division(mesh->vertices[i]);
//        //            mesh->vertices[i].screen = s->camera->M_viewport * mesh->vertices[i].screen;
//        //        }
//
//#pragma omp parallel for
//        for (int i = 0; i < mesh->num_triangle; i++)
//        {
//            auto v_0 = VertexP(mesh->vertices[mesh->triangles[i].vertex_0]);
//            auto v_1 = VertexP(mesh->vertices[mesh->triangles[i].vertex_1]);
//            auto v_2 = VertexP(mesh->vertices[mesh->triangles[i].vertex_2]);
//            s->shader->vertex_shader(v_0);
//            s->shader->vertex_shader(v_1);
//            s->shader->vertex_shader(v_2);
//
//            if (!clipping(v_0, v_1, v_2))
//            {
//                // perspective division
//				Render::perspective_division(v_0);
//				Render::perspective_division(v_1);
//				Render::perspective_division(v_2);
//                // screen space coordinate
//                v_0.screen = s->camera->M_viewport * v_0.clip * v_0.z_rec;
//                v_1.screen = s->camera->M_viewport * v_1.clip * v_1.z_rec;
//                v_2.screen = s->camera->M_viewport * v_2.clip * v_2.z_rec;
//
//                // face culling
//                if (s->face_cull_mode == BACK)
//                {
//                    float AB_x = v_1.screen.x() - v_0.screen.x();
//                    float AB_y = v_1.screen.y() - v_0.screen.y();
//                    float AC_x = v_2.screen.x() - v_0.screen.x();
//                    float AC_y = v_2.screen.y() - v_0.screen.y();
//                    if (AB_x * AC_y - AB_y * AC_x > 0)
//                    {
//                        continue;
//                    }
//                } else if (s->face_cull_mode == FRONT)
//                {
//                    float AB_x = v_1.screen.x() - v_0.screen.x();
//                    float AB_y = v_1.screen.y() - v_0.screen.y();
//                    float AC_x = v_2.screen.x() - v_0.screen.x();
//                    float AC_y = v_2.screen.y() - v_0.screen.y();
//                    if (AB_x * AC_y - AB_y * AC_x < 0)
//                    {
//                        continue;
//                    }
//                }
//                num_triangle++;
//                draw_triangle(v_0, v_1, v_2, mesh->triangles[i].normal);
//            } else
//            {
//                vector<VertexP> vertices = clip_near(v_0, v_1, v_2);
//                for(auto &v: vertices)
//                {
//					Render::perspective_division(v);
//                    v.screen = s->camera->M_viewport * v.clip * v.z_rec;
//                }
//
//                int num_clip_tri = vertices.size();
//                for (int j = 0; j < num_clip_tri; j += 2)
//                {
//                    const VertexP &v0 = vertices[j % num_clip_tri];
//                    const VertexP &v1 = vertices[(j + 1) % num_clip_tri];
//                    const VertexP &v2 = vertices[(j + 2) % num_clip_tri];
//
//                    // face culling
//                    if (s->face_cull_mode == BACK)
//                    {
//                        float AB_x = v1.screen.x() - v0.screen.x();
//                        float AB_y = v1.screen.y() - v0.screen.y();
//                        float AC_x = v2.screen.x() - v0.screen.x();
//                        float AC_y = v2.screen.y() - v0.screen.y();
//                        if (AB_x * AC_y - AB_y * AC_x > 0)
//                        {
//                            continue;
//                        }
//                    } else if (s->face_cull_mode == FRONT)
//                    {
//                        float AB_x = v1.screen.x() - v0.screen.x();
//                        float AB_y = v1.screen.y() - v0.screen.y();
//                        float AC_x = v2.screen.x() - v0.screen.x();
//                        float AC_y = v2.screen.y() - v0.screen.y();
//                        if (AB_x * AC_y - AB_y * AC_x < 0)
//                        {
//                            continue;
//                        }
//                    }
//                    num_triangle++;
//                    draw_triangle(v0, v1, v2, mesh->triangles[i].normal);
//                }
//            }
//        }
//    }
//
//    //#pragma omp parallel for
//    //        for (int i = 0; i < mesh->num_vertex; i++)
//    //        {
//    //            perspective_restore(mesh->vertices[i]); // perspective restore
//    //        }
//}

void Renderer::draw_triangle(const VertexP &v_0, const VertexP &v_1, const VertexP &v_2, const float4 &normal)
{
    // real bounding box
    int min_x = max((int)min(v_0.screen.x(), min(v_1.screen.x(), v_2.screen.x())), 0);
    int min_y = max((int)min(v_0.screen.y(), min(v_1.screen.y(), v_2.screen.y())), 0);
    int max_x = min((int)max(v_0.screen.x(), max(v_1.screen.x(), v_2.screen.x())) + 1, s->camera->x - 1);
    int max_y = min((int)max(v_0.screen.y(), max(v_1.screen.y(), v_2.screen.y())) + 1, s->camera->y - 1);

    int index = min_y * s->camera->x;
    for (int i = min_y; i <= max_y; i++)
    {
        for (int j = min_x; j <= max_x; j++)
        {
            float center_x = float(j) + 0.5f;
            float center_y = float(i) + 0.5f;
            float v0x = v_0.screen.x() - center_x;
            float v0y = v_0.screen.y() - center_y;
            float v1x = v_1.screen.x() - center_x;
            float v1y = v_1.screen.y() - center_y;
            float v2x = v_2.screen.x() - center_x;
            float v2y = v_2.screen.y() - center_y;
            float AB = v1x * v0y - v1y * v0x;
            float BC = v2x * v1y - v2y * v1x;
            float CA = v0x * v2y - v0y * v2x;
            if (AB > 0 && BC > 0 && CA > 0)
            {
                float S = 1.0f / (AB + BC + CA);
                float u = BC * S;
                float v = CA * S;
                float w = AB * S;
                float z = Render::lerp(v_0.screen.z(), v_1.screen.z(), v_2.screen.z(), u, v, w);
                // early-z test
                if (z < s->z_buffer[index + j])
                {
                    num_frag++;
                    // perspective correct interpolation
                    Fragment frag = Render::lerp(v_0, v_1, v_2, u, v, w);
                    frag.x = j;
                    frag.y = i;
                    frag.z = z;
                    frag.flat_normal = normal;

                    // calculate texture space mapping
                    if (s->texture_type == MIPMAP)
                    {
                        center_x = float(j) + 1.5f;
                        center_y = float(i) + 0.5f;
                        v0x = v_0.screen.x() - center_x;
                        v0y = v_0.screen.y() - center_y;
                        v1x = v_1.screen.x() - center_x;
                        v1y = v_1.screen.y() - center_y;
                        v2x = v_2.screen.x() - center_x;
                        v2y = v_2.screen.y() - center_y;
                        AB = v1x * v0y - v1y * v0x;
                        BC = v2x * v1y - v2y * v1x;
                        CA = v0x * v2y - v0y * v2x;
                        S = 1.0f / (AB + BC + CA);
                        u = BC * S;
                        v = CA * S;
                        w = AB * S;
                        frag.texture_x = Render::lerp(v_0.texture_uv, v_1.texture_uv, v_2.texture_uv, u, v, w);

                        center_x = float(j) + 0.5f;
                        center_y = float(i) + 1.5f;
                        v0x = v_0.screen.x() - center_x;
                        v0y = v_0.screen.y() - center_y;
                        v1x = v_1.screen.x() - center_x;
                        v1y = v_1.screen.y() - center_y;
                        v2x = v_2.screen.x() - center_x;
                        v2y = v_2.screen.y() - center_y;
                        AB = v1x * v0y - v1y * v0x;
                        BC = v2x * v1y - v2y * v1x;
                        CA = v0x * v2y - v0y * v2x;
                        S = 1.0f / (AB + BC + CA);
                        u = BC * S;
                        v = CA * S;
                        w = AB * S;
                        frag.texture_y = Render::lerp(v_0.texture_uv, v_1.texture_uv, v_2.texture_uv, u, v, w);
                    }

                    // fragment perspective restore
					Render::perspective_restore(frag);
                    s->shader->fragment_shader(frag);
                    write_color(frag);
                }
            }
        }
        index += s->camera->x;
    }
}


inline void Renderer::write_color(Fragment &frag) const
{
    int index = frag.y * s->camera->x + frag.x;
    if (frag.z < s->z_buffer[index])
    {
        s->z_buffer[index] = frag.z;
        s->frame_buffer->write_color(frag);
    }
}

