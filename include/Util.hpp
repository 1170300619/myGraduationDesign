#ifndef RENDERER_UTIL_HPP
#define RENDERER_UTIL_HPP

#include "Fragment.hpp"
#include "Vertex.hpp"
#include <string>
#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Type.hpp"
#include <io.h>
#include <direct.h>
#include <sys/timeb.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

namespace Render {
	static const float kPi = 3.1415926536;
	static const int kImageSize = 1024;
	static const int kDEMSize = 99;
	static const float kMinIncludedAngle = 0.000000001;
	static const float kOrthoRemoveRatio = 0.1;
	static const float kMinAlbedo = 0.9;
	static const float kMaxAlbedo = 1.1;
	static const int kBorderSize = 1;
    static const float kModelGSD = 0.04;
    static const float4 kZeroFloat4(0, 0, 0, 0);

    static const string kInitModelPath = "C:/SPC_Files/Bennu_Preliminary_Model/50x.ply";
	static const string kLandmarksPath = "./data/Landmark/";
	static const string kRandomPointsPath = "./data/randomPoints.txt";
	static const string kPlyPath = "./data/Model.ply";
//	static const string kTestDEMPath = "./data/testDEM.txt";
    static const string kFB1Path = "C:/SPC_Files/SPC_generate_file/FB1/";
    static const string kFB5APath = "C:/SPC_Files/SPC_generate_file/FB5-A/";
    static const string kFB6BPath = "C:/SPC_Files/SPC_generate_file/FB6-B/";
    static const string kFirstOptimationPath = "C:/SPC_Files/SPC_generate_file/optimation_0/";
//    static const string kSecondOptimationPath = "C:/SPC_Files/SPC_generate_file/optimation_1/";
    static const string kAlbedoPath = "C:/SPC_Files/SPC_generate_file/albedo_image/";
    static const string kNoiseCamAndSunMes = "./data/noiseCamAndSunMes/";
    static const string kBundleAdjustmentPath = "C:/SPC_Files/SPC_generate_file/BundleAdjustment/";

	enum {
		FIRST = 0,
		SECOND
	};

    enum processType{
        ALBEDO = 0,
        FB1,
        FB5A,
        FB6B
    };

    static const string EnumStrings[] = {"ALBEDO", "FB1", "FB5A", "FB6B"};

    inline string getTextFromEnum(processType enumVal)
    {
        return EnumStrings[enumVal];
    }

    inline string getPathFromEnum(processType enumVal){
        if(enumVal == ALBEDO)return kAlbedoPath;
        else if(enumVal == FB1)return kFB1Path;
        else if(enumVal == FB5A)return kFB5APath;
        else return kFB6BPath;
    }

	/**
	*
	* @param deg
	* @return
	*/
	inline float deg2rad(float deg)
	{
		return deg * kPi / 180.0f;
	}

    inline void copyDEM(const string& srcPath, const string& dstPath){
        if(_access(srcPath.c_str(), 0) == -1) {
            cout << "srcPath: " << srcPath << "doesn't exist" << endl;
            return;
        }
        ifstream in;
        ofstream out;
        in.open((srcPath).c_str());
        out.open(dstPath.c_str());
        out << in.rdbuf();
        in.close(), out.close();
    }

    inline float4 double4ToFloat4(const double4& dp){
        float4 fp;
        fp << (float)dp[0], (float)dp[1], (float)dp[2], (float)dp[3];
        return fp;
    }

    inline string formatNum(int num){
        stringstream ss;
        ss << setw(5) << setfill('0') << num ;
        string str;
        ss >> str;
        return str;
    }

    inline float4 localToBody(const float4& u1, const float4& u2, const float4& u3, const float4& landmark, const float4& localPoint){
        float3x3 axis;
        axis << u1[0], u1[1], u1[2], u2[0], u2[1], u2[2], u3[0], u3[1], u3[2];
        axis.transposeInPlace();
        float3 lp = float3(localPoint[0], localPoint[1], localPoint[2]);
        float3 worldPoint = axis * lp + float3(landmark[0], landmark[1], landmark[2]);
        return {worldPoint[0], worldPoint[1], worldPoint[2], 1.0f};
    }

	/**
	*
	* @param rad
	* @return
	*/
	inline float rad2deg(float rad)
	{
		return rad * 180.0f / kPi;
	}

	/**
	*
	* @tparam T
	* @param x
	* @param v_0
	* @param v_1
	* @return
	*/
	template<typename T> inline T lerp(float x, T v_0, T v_1)
	{
		return v_0 + x * (v_1 - v_0);
	}

	/**
	*
	* @tparam T
	* @param v0
	* @param v1
	* @param v2
	* @param u
	* @param v
	* @param w
	* @return
	*/
	template<typename T> inline T lerp(const T &v0, const T &v1, const T &v2, float u, float v, float w)
	{
		return u * v0 + v * v1 + w * v2;
	}

	/**
	*
	* @param v_0
	* @param v_1
	* @param v_2
	* @param u
	* @param v
	* @param w
	* @return
	*/
	inline Fragment lerp(const VertexP &v_0, const VertexP &v_1, const VertexP &v_2, float u, float v, float w)
	{
		Fragment frag;
		frag.world = lerp(v_0.position, v_1.position, v_2.position, u, v, w);
		frag.normal = lerp(v_0.normal, v_1.normal, v_2.normal, u, v, w);
		frag.color = lerp(v_0.color, v_1.color, v_2.color, u, v, w);
		frag.texture_uv = lerp(v_0.texture_uv, v_1.texture_uv, v_2.texture_uv, u, v, w);
		frag.clip_z = 1.0f / lerp(v_0.z_rec, v_1.z_rec, v_2.z_rec, u, v, w);
		return frag;
	}

	/**
	*
	* @param x
	* @param v_0
	* @param v_1
	* @return
	*/
	inline VertexP lerp(float x, const VertexP &v_0, const VertexP &v_1)
	{
		VertexP v;
		v.position = lerp(x, v_0.position, v_1.position);
		v.clip = lerp(x, v_0.clip, v_1.clip);
		v.color = lerp(x, v_0.color, v_1.color);
		v.texture_uv = lerp(x, v_0.texture_uv, v_1.texture_uv);
		v.z_rec = 1.0f / v.clip.w();
		return v;
	}

	/**
	* Judge if a point is in a triangle.
	* @param AB
	* @param BC
	* @param CA
	* @param u
	* @param v
	* @param w
	* @return
	*/
	inline bool is_in_triangle(float AB, float BC, float CA, float &u, float &v, float &w)
	{
		float S = 1.0f / (AB + BC + CA);

		u = BC * S;
		if (u < 0 || u > 1.0f)
		{
		    return false;
		}
	
		v = CA * S;
		if (v < 0 || v > 1.0f)
		{
		    return false;
		}
	
		w = AB * S;
		if (w < 0 || w > 1.0f)
		{
		    return false;
		}

		return true;
	}

	/**
	*
	* @param vertex
	*/
	inline void perspective_division(VertexP &vertex)
	{
		vertex.position *= vertex.z_rec;
		vertex.normal *= vertex.z_rec;
		vertex.color *= vertex.z_rec;
		vertex.texture_uv *= vertex.z_rec;
	}

	/**
	*
	* @param frag
	*/
	inline void perspective_restore(Fragment &frag)
	{
		frag.world *= frag.clip_z;
		frag.normal *= frag.clip_z;
		frag.color *= frag.clip_z;
		frag.texture_uv *= frag.clip_z;
		frag.texture_x *= frag.clip_z;
		frag.texture_y *= frag.clip_z;
	}

	/**
	*
	* @param vertex
	*/
	inline void perspective_restore(VertexP &vertex)
	{
		float z = 1.0f / vertex.z_rec;
		vertex.position *= z;
		vertex.normal *= z;
		vertex.color *= z;
		vertex.texture_uv *= z;
	}

	/**
	* 检测射线和三角面是否有交点
	* FullName:  IntersectTriangle
	* Access:    inline
	* Returns:   bool
	* Param:     const float4 orig 射线出发点
	* Param:	 const float4 dir 射线方向
	* Param:     Vertex v0 三角面角点
	* Param:     Vertex v1 三角面角点
	* Param:     Vertex v1 三角面角点
	* Param:     float t[] 中间参数(如果有交点，代表从射线出发点到交点的距离)
	*/
	inline bool IntersectTriangle(const float4 orig, const float4 dir, const Vertex v0, const Vertex v1, const Vertex v2, float t[]) {
		float u, v;
		float4 _v0 = v0.position;
		float4 _v1 = v1.position;
		float4 _v2 = v2.position;
		float4 E1 = _v1 - _v0;
		float4 E2 = _v2 - _v0;
		float4 P = dir.cross3(E2);
		// determinant
		float det = E1.dot(P);
		// keep det > 0, modify T accordingly
		float4 _T;
		if (det >0)
		{
			_T = orig - _v0;
		}
		else
		{
			_T = _v0 - orig;
			det = -det;
		}
		// Calculate u and make sure u <= 1
		u = _T.dot(P);
		if (u < 0.0f || u > det)
			return false;
		// Q
		float4 Q = _T.cross3(E1);
		// Calculate v and make sure u + v <= 1
		v = dir.dot(Q);
		if (v < 0.0f || u + v > det)
			return false;
		// Calculate t, scale parameters, ray intersects triangle
		t[0] = E2.dot(Q);
		float fInvDet = 1.0f / det;
		t[0] *= fInvDet;
		u *= fInvDet;
		v *= fInvDet;
		return true;
	}

	/**
	* 计算三维点间距离
	* FullName:  Dist
	* Access:    inline
	* Returns:   float
	* Param:     float4 &v1 A点
	* Param:	 float4 &v2 B点
	*/
	inline float Dist(float4 &v1, float4 &v2)
	{
		return (v1 - v2).lpNorm<2>();
	}


	// 生成均值为 avg, 标准差为 standarddeviation 的高斯分布
	inline float gaussrand(float avg, float staDev)
	{
        std::random_device rd;
        std::default_random_engine rng(rd());
        normal_distribution<> u{avg, staDev};
        return (float)u(rng);
	}

	// 创建多级目录
	inline bool create_dir(const string& path)
	{
		int m = 0, n;
		string str1, str2;

		str1 = path;
		str2 = str1.substr(0, 2);
		str1 = str1.substr(3, str1.size());
		cout << str2 << " " << str1 << endl;
		while (m >= 0)
		{
			m = str1.find("/");
			str2 += "/" + str1.substr(0, m);
			if (m == string::npos) {
				break;
			}
			n = _access(str2.c_str(), 0); // 判断该目录是否存在
			if (n == -1)
			{
				if (_mkdir(str2.c_str()) != 0)     // 创建目录
				{
					cout << "fail to mkdir" << endl;
					return false;
				}
			}
			str1 = str1.substr(m + 1, str1.size());
		}
		return true;
	}

	// 检查绝对路径是否存在，如果不存在则创建该路径
	inline bool savePathCheack(const string& path)
	{
		if (path.length() <= 0)
			return false;

		errno_t err = 0;  
		if ((err = _access_s(path.c_str(), 0)) != 0)
		{
			if (!create_dir(path)) {
				std::cout << "create filepath failed" << endl;
				return false;
			}
			else
			{
				std::cout << "create filepath success : " << path << std::endl;
				return true;
			}
		}
		else
			std::cout << path << " already exist" << std::endl;

		return true;
	}

    inline void TimeMilliSecond()
    {
        timeb now;
        ftime(&now);
        cout <<  now.time * 1000 + now.millitm << endl;
    }

    inline void fibonacciGrid(vector<float4>& randomPoint, int n){
        int index = 0;
        auto theta = float(0.5 * (sqrt(5) - 1)); // 黄金分割数
        while(index++ < n){
            float z = float(2 * index - 1) / float(n) - 1;
            float x = sqrt(1 - z * z) * cos(2 * kPi * float(index) * theta);
            float y = sqrt(1 - z * z) * sin(2 * kPi * float(index) * theta);
            randomPoint.emplace_back(x, y, z, 1.0f);
        }
    }

    inline void getPngFiles(string& filePath, vector<int>& fileNum){
        intptr_t hFile = 0;
        _finddata_t fileInfo;
        hFile = _findfirst(filePath.c_str(), &fileInfo);
        int indexNum = 0;
        if (hFile != -1){
            do{
                //如果为文件夹，可以通过递归继续遍历，此处我们不需要
                if ((fileInfo.attrib &  _A_SUBDIR)){
                    continue;
                }
                    //如果为单个文件直接push_back
                else{
                    indexNum++;
                }

            } while (_findnext(hFile, &fileInfo) ==0);
            if(indexNum > 0 && indexNum <= 4)fileNum[0]++;
            else if(indexNum > 4 && indexNum <= 8)fileNum[1]++;
            else if(indexNum > 8 && indexNum <= 12)fileNum[2]++;
            else if(indexNum > 12 && indexNum <= 16)fileNum[3]++;
            else if(indexNum > 16)fileNum[4]++;
            _findclose(hFile);
        }
    }
}


#endif //RENDERER_UTIL_HPP