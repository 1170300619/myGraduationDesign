#pragma warning(disable:4996);
#include "../include/Lmap.hpp"
#include <direct.h>

using namespace std;
using namespace Render;

Lmap::Lmap() {
	landmark << 0, 0, 0, 1.0f;
	u1 << 1.0f, 0, 0, 0;
	u2 << 0, 1.0f, 0, 0;
	u3 << 0, 0, 1.0f, 0;

	int size = kDEMSize * kDEMSize;
	std::vector<double4> DEM(size);
	std::vector<double> albedoMap(size);
	ptr = new U_Ptr(DEM, albedoMap);
}


Lmap::Lmap(float4 &controlPoint) {
	landmark = controlPoint;
	//u1 << 1.0f, 0, 0, 0;
	//u2 << 0, 1.0f, 0, 0;
	//u3 << 0, 0, 1.0f, 0;

	//ʹ�������еľֲ�����ϵ���巽��
    float3 _u3 = float3(controlPoint[0], controlPoint[1], controlPoint[2]) / sqrt(controlPoint[0] * controlPoint[0] + controlPoint[1] * controlPoint[1] + controlPoint[2] * controlPoint[2]);
	u3 = float4(_u3[0], _u3[1], _u3[2], 0);

	u1 << -controlPoint[1], controlPoint[0], 0, 0;
	u1 = u1 / sqrt(controlPoint[0] * controlPoint[0] + controlPoint[1] * controlPoint[1]);
	u2 = u3.cross3(u1);
	u2[3] = 0;

	int size = kDEMSize * kDEMSize;
	vector<double4> DEM(size);
	vector<double> albedoMap(size);
	ptr = new U_Ptr(DEM, albedoMap);
}


Lmap& Lmap::operator=(const Lmap &lmap) {
	++lmap.ptr->use;
	if (--ptr->use == 0)delete ptr;
	ptr = lmap.ptr;
	landmark = lmap.landmark;
	u1 = lmap.u1;
	u2 = lmap.u2;
	u3 = lmap.u3;
	return *this;
}


Lmap::~Lmap() {
	if (--ptr->use == 0)
		delete ptr;
}


void Lmap::writeLmap(const string& lmapPath) {
	ofstream out;
	if (!savePathCheack(lmapPath)) {
		cerr << "save path check failed" << endl;
		throw exception();
	}
	out.open(lmapPath.c_str());
	if (!out.is_open()) {
		cerr << "Fail to open lmapPath" << endl;
		throw exception();
	}

	int size = kDEMSize * kDEMSize;
	out << kDEMSize << " " << kDEMSize << endl;
	out << landmark[0] << " " << landmark[1] << " " << landmark[2] << " " << landmark[3] << endl;
	out << u1[0] << " " << u1[1] << " " << u1[2] << " " << u1[3] << endl;
	out << u2[0] << " " << u2[1] << " " << u2[2] << " " << u2[3] << endl;
	out << u3[0] << " " << u3[1] << " " << u3[2] << " " << u3[3] << endl;

	for (int i = 0; i < size; i++) {
		out << ptr->DEM[i][0] << " " << ptr->DEM[i][1] << " " << ptr->DEM[i][2] << " " << ptr->DEM[i][3] << endl;
	}
	for (int i = 0; i < size - 1; i++) {
		out << ptr->albedoMap[i] <<  endl;
	}
    out << ptr->albedoMap[size - 1]; //��ֹ������
	out.close();
}


void Lmap::readLmap(const string &LmapPath) {
	ifstream infile;
	infile.open(LmapPath.c_str());
	if (!infile)
	{
		cout << "open file error!" << endl;
	}
	vector<float> v;
	double x, y, z, w;
    double albedo;
	string s;
	getline(infile, s);
	infile >> landmark(0) >> landmark(1) >> landmark(2) >> landmark(3);
	//ע���landmark�޸�
	landmark << landmark(0), landmark(1), landmark(2), landmark(3);
	infile >> u1(0) >> u1(1) >> u1(2) >> u1(3);
	infile >> u2(0) >> u2(1) >> u2(2) >> u2(3);
	infile >> u3(0) >> u3(1) >> u3(2) >> u3(3);

	for(int i = 0; i < kDEMSize * kDEMSize; ++i){
        infile >> x >> y >> z >> w;
        ptr->DEM[i] << x, y, z, w;
    }
    for(int i = 0; i < kDEMSize * kDEMSize; ++i){
        infile >> albedo;
        ptr->albedoMap[i] = albedo;
    }
	infile.close();
}


bool Lmap::GenIniDEMFromBennu(Model* model, float modelGSD)
{
    float L = modelGSD * kDEMSize;
	vector<float4> worldPoints;
	for (int i = 0; i < kDEMSize; ++i) {
		for (int j = 0; j < kDEMSize; ++j) {
			float4 p;
			p << (j * modelGSD - L / 2), (-i * modelGSD + L / 2), 0, 1.0f;
            worldPoints.emplace_back(localToBody(u1, u2, u3, landmark,p));
		}
	}

    Ray ray = Ray(kZeroFloat4, kZeroFloat4);
    HitRecord record;
    float4 direction;

	for (int j = 0; j < kDEMSize; j++) {
		for (int k = 0; k < kDEMSize; k++) {
			direction << worldPoints[j * kDEMSize + k][0], worldPoints[j * kDEMSize + k][1], worldPoints[j * kDEMSize + k][2], 0;
            direction.normalize();
            ray = Ray(float4(0, 0, 0, 1.0f), direction);
            if(!model->hit(ray, 0.001, 100000.0f, record))return false;
            record.position = ray.origin + (record.t + 0.01) * ray.direction;
            auto x = u1.dot(record.position - landmark);
            auto y = u2.dot(record.position - landmark);
            auto z = u3.dot(record.position - landmark);
            ptr->DEM[j * kDEMSize + k] << x, y, z, 1.0f;
		}
	}
	cout << "generate initial DEM from Bennu success" << endl;
    return true;
}


//void Lmap::calcLkPixelCoordInImage(
//	const string &renderMapPath, Mat &image,
//	float &ObjectX, float &ObjectY) {
//	/*��Ҫ�ǽ��ж���֮���ƥ��
//	�����landmark��imageData�е���������ObjectX��ObjectY
//	*/
//	Mat renderMap = cv::imread(renderMapPath, 0);//����У��ͼ��
//
//	Mat srcImage;
//	image.copyTo(srcImage);
//
//	//�����renderMap��������image(ԭʼ����ͼ��)�ϵ��������꣬��ΪObjectX��ObjectY
//	Mat resultImage;
//	int match_method = 0;
//	int resImageCols = image.cols - renderMap.cols + 1;
//	int resImageRows = image.rows - renderMap.rows + 1;
//	resultImage.create(resImageCols, resImageRows, CV_32FC1);
//	matchTemplate(image, renderMap, resultImage, match_method);
//	normalize(resultImage, resultImage, 0, 1, NORM_MINMAX, -1, Mat());
//	float minVal, maxVal;
//
//	Point minLocation, maxLocation, matchLocation;
//	minMaxLoc(resultImage, &minVal, &maxVal, &minLocation, &maxLocation);
//
//	if (match_method == TM_SQDIFF || match_method == CV_TM_SQDIFF_NORMED) {
//		matchLocation = minLocation;
//	}
//	else
//		matchLocation = maxLocation;
//	ObjectX = matchLocation.x + renderMap.cols / 2;
//	ObjectY = matchLocation.y + renderMap.rows / 2;
//
//	/*rectangle(image, matchLocation, Point(matchLocation.x + renderMap.cols, matchLocation.y + renderMap.rows), Scalar(0, 0, 255), 2, 8, 0);
//	rectangle(resultImage, matchLocation, Point(matchLocation.x + renderMap.cols, matchLocation.y + renderMap.rows), Scalar(0, 0, 255), 2, 8, 0);
//
//	namedWindow("ԭʼͼ", WINDOW_AUTOSIZE);
//	imshow("ԭʼͼ", srcImage);
//	waitKey(0);
//	namedWindow("ƥ��ͼ", WINDOW_AUTOSIZE);
//	imshow("ƥ��ͼ", image);
//	waitKey(0);
//	namedWindow("Ч��ͼ", WINDOW_AUTOSIZE);
//	imshow("Ч��ͼ", resultImage);*/
//
//}