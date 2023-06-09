//#pragma once
//#ifndef ICQ_H_
//#define ICQ_H_
//
//#include <vector>
//#include <Eigen/Dense>
//#include <cmath>
//#include "opencv2/imgcodecs.hpp"
//#include "opencv2/highgui.hpp"
//#include <iostream>
//#include "Lmap.hpp"
//
//
//class ICQ {
//public:
//	ICQ();
//	ICQ(float4 _center, int _q, float _side_l);
//	~ICQ();
//
//	vector<float4> getVectorList();
//	float4 getVertex(int i, int j, int f);
//	int getVertexLabel(int i, int j, int f) const;
//	int getFaceLabel(int i, int j, int f) const;
//	float4 getNormal(int i, int j, int f);
//
//	static vector<float4> read_dem(const string &path);
//	float getIntersectHeight(int i, int j, int f,int size, double *DEM, Lmap &lmap, float Dstep);
//	float getIntersectH(int i, int j, int f, string prefix);
//
//	//void updateVp(vector<Vector3d> &pointList, int size,
//		//vector<double**> DEMlist, double Dstep, int scale);
//	void updateVp(int size, double *DEM, Lmap &lmap, float Dstep);
//	void updateVp(string prefix);
//
//	void saveModel(const char * DistName);
//	void saveModel_txt(const char * DistName);
//
//	float vp_change_aver(string preVerPath);
//	void densify_model();
//	void read_model(string path);
//
//	int get_q() { return q; }
//private:
//	float4 center;
//	vector<float4> vectorList;
//	int q;
//	float sideL;
//};
//
//#endif // !MODEL_H_
//
