//#include "../include/ICQ.hpp"
//#include <map>
//#include <cmath>
//#include <algorithm>
//#include <utility>
//#include "../include/Util.hpp"
//
//using namespace std;
//
//ICQ::ICQ() {
//    q = 0;
//    sideL = 0;
//}
//
//ICQ::~ICQ() = default;
//
//ICQ::ICQ(float4 _center, int _q, float _sideL) {
//	center = std::move(_center);
//	q = _q;
//    sideL = _sideL;
//	//保存点的三维坐标
//	float step = sideL / (float)q;
//
//	//第一个面（上面）
//	for (int j = 0; j < q + 1; j++) {
//		for (int i = 0; i < q + 1; i++) {
//			float4 v(float(i - q * 0.5) * step, float(sideL * 0.5), float(j - q * 0.5) * step, 0);
//			vectorList.push_back(v);
//		}
//	}
//
//	//第二个面
//	for (int j = 0; j < q + 1; j++) {
//		for (int i = 0; i < q + 1; i++) {
//			float4 v(float(i - q * 0.5) * step, float(q * 0.5 - j) * step, float(sideL * 0.5), 0);
//			vectorList.push_back(v);
//		}
//	}
//
//	//第三个面
//	for (int j = 0; j < q + 1; j++) {
//		for (int i = 0; i < q + 1; i++) {
//			float4 v(-sideL * 0.5f, float(q * 0.5 - j) * step, float(i - q * 0.5) * step, 0);
//			vectorList.push_back(v);
//		}
//	}
//
//	//第四个面
//	for (int j = 0; j < q + 1; j++) {
//		for (int i = 0; i < q + 1; i++) {
//			float4 v(sideL * 0.5f, float(q * 0.5 - j) * step, float(q * 0.5 - i) * step, 0);
//			vectorList.push_back(v);
//		}
//	}
//
//	//第五个面
//	for (int j = 0; j < q + 1; j++) {
//		for (int i = 0; i < q + 1; i++) {
//			float4 v(float(q * 0.5 - i) * step, float(q * 0.5 - j) * step, float(-sideL * 0.5), 0);
//			vectorList.push_back(v);
//		}
//	}
//
//	//第六个面
//	for (int j = 0; j < q + 1; j++) {
//		for (int i = 0; i < q + 1; i++) {
//			float4 v(float(i - q * 0.5) * step, float(-sideL * 0.5), float(q * 0.5 - j) * step, 0);
//			vectorList.push_back(v);
//		}
//	}
//}
//
//// 获取空间三维位置
//float4 ICQ::getVertex(int i, int j, int f) {
//	int index = getVertexLabel(i, j, f);
//	if ((index < 0) || (index > vectorList.size())) {
//		return {0, 0, 0, 0};
//	}
//	else {
//		return vectorList[index];
//	}
//}
//
//// 节点标签
//int ICQ::getVertexLabel(int i, int j, int f) const {
//	return i + (q + 1) * j + (q + 1) * (q + 1) * (f - 1);
//}
//
//// 面标签
//int ICQ::getFaceLabel(int i, int j, int f) const {
//	return i + q * (j - 1) + q * q * (f - 1);
//}
//
//// 向量vector
//vector<float4> ICQ::getVectorList() {
//	return vectorList;
//}
//
////计算法线
//float4 ICQ::getNormal(int i, int j, int f) {
//	//上层4个重合顶点
//	if ((i == 0 && j == 0 && f == 2) || (i == q && j == 0 && f == 3))
//		return getNormal(0, q, 1);
//	if ((i == 0 && j == 0 && f == 5) || (i == q && j == 0 && f == 2))
//		return getNormal(q, q, 1);
//	if ((i == 0 && j == 0 && f == 4) || (i == q && j == 0 && f == 5))
//		return getNormal(q, 0, 1);
//	if ((i == 0 && j == 0 && f == 3) || (i == q && j == 0 && f == 4))
//		return getNormal(0, 0, 1);
//	//下层4个重合顶点
//	if ((i == 0 && j == q && f == 2) || (i == q && j == q && f == 3))
//		return getNormal(0, 0, 6);
//	if ((i == 0 && j == q && f == 3) || (i == q && j == q && f == 4))
//		return getNormal(0, q, 6);
//	if ((i == 0 && j == q && f == 5) || (i == q && j == q && f == 2))
//		return getNormal(q, 0, 6);
//	if ((i == 0 && j == q && f == 4) || (i == q && j == q && f == 5))
//		return getNormal(q, q, 6);
//
//	//处理侧边四条棱
//	if ((i == q) && (f == 3))return getNormal(0, j, 2);
//	if ((i == 0) && (f == 5))return getNormal(q, j, 2);
//	if ((i == 0) && (f == 3))return getNormal(q, j, 4);
//	if ((i == q) && (f == 5))return getNormal(0, j, 4);
//
//	float4 v1, v2;
//	if (f == 1) {
//		if (i == 0) {
//			if (j == 0) {
//				v1 = getVertex(1, 0, 1) - getVertex(0, 1, 1);
//				v2 = getVertex(0, 1, 3) - getVertex(0, 1, 1);
//			}
//			else if (j == q) {
//				v1 = getVertex(1, 0, 2) - getVertex(0, 1, 2);
//				v2 = getVertex(0, q - 1, 1) - getVertex(0, 1, 2);
//			}
//			else {//(i,j,f)=(0,n,1)
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(j, 1, 3);
//			}
//		}
//		else if (i == q) {
//			if (j == 0) {
//				v1 = getVertex(q, 1, 5) - getVertex(q, 1, 1);
//				v2 = getVertex(q - 1, 0, 1) - getVertex(q, 1, 1);
//			}
//			else if (j == q) {
//				v1 = getVertex(1, 0, 5) - getVertex(q, 1, 2);
//				v2 = getVertex(q - 1, 0, 2) - getVertex(q, 1, 2);
//			}
//			else {//(i,j,f)=(q,n,1)
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(q - j, 1, 5) - getVertex(i - 1, j, f);
//			}
//		}
//		else {//i!=0 && i!=q
//			if (j == 0) {
//				v1 = getVertex(i, j + 1, f) - getVertex(q - i, 1, 4);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//			else if (j == q) {
//				v1 = getVertex(i, 1, 2) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//			else {
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//		}
//	}
//	else if (f == 6) {
//		if (i == 0) {
//			if (j == 0) {
//				v1 = getVertex(0, 1, 6) - getVertex(0, q - 1, 2);
//				v2 = getVertex(1, 0, 6) - getVertex(0, q - 1, 2);
//			}
//			else if (j == q) {
//				v1 = getVertex(1, q, 6) - getVertex(0, q - 1, 3);
//				v2 = getVertex(0, q - 1, 6) - getVertex(0, q - 1, 3);
//			}
//			else {
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(q - j, q - 1, 3);
//			}
//		}
//		else if (i == q) {
//			if (j == 0) {
//				v1 = getVertex(q - 1, 0, 6) - getVertex(0, q - 1, 5);
//				v2 = getVertex(1, q, 5) - getVertex(0, q - 1, 5);
//			}
//			else if (j == q) {
//				v1 = getVertex(q - 1, q, 5) - getVertex(0, q - 1, 4);
//				v2 = getVertex(1, q, 4) - getVertex(0, q - 1, 4);
//			}
//			else {
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(j, q - 1, 5) - getVertex(i - 1, j, f);
//			}
//		}
//		else {
//			if (j == 0) {
//				v1 = getVertex(i, j + 1, f) - getVertex(i, q - 1, 2);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//			else if (j == q) {
//				v1 = getVertex(q - i, q - 1, 4) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//			else {
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//		}
//	}
//	else if (f == 2) {
//		if (i == 0) {
//			v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//			//cout << "v1=" << v1 << endl;
//			v2 = getVertex(1, j, f) - getVertex(q - 1, j, 3);
//			//cout << "v2=" << v2 << endl;
//		}
//		else if (i == q) {
//			v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//			v2 = getVertex(1, j, 5) - getVertex(q - 1, j, f);
//		}
//		else {
//			if (j == 0) {
//				return getNormal(i, q, 1);
//			}
//			else if (j == q) {
//				return getNormal(i, 0, 6);
//			}
//			else {
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//		}
//	}
//
//	else if (f == 4) {
//		if (i == 0) {
//			v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//			v2 = getVertex(1, j, f) - getVertex(q - 1, j, 5);
//		}
//		else if (i == q) {
//			v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//			v2 = getVertex(1, j, 3) - getVertex(q - 1, j, f);
//		}
//		else {
//			if (j == 0) {
//				return getNormal(q - i, 0, 1);
//			}
//			else if (j == q) {
//				return getNormal(q - i, q, 6);
//			}
//			else {
//				v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//				v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//			}
//		}
//	}
//	else if(f==3){
//		if (j == 0) {
//			return getNormal(0, i, 1);
//		}
//		else if (j == q) {
//			return getNormal(0, q - i, 6);
//		}
//		else {
//			v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//			v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//		}
//	}
//	else if (f == 5) {
//		if (j == 0) {
//			return getNormal(q, q - i, 1);
//		}
//		else if (j == q) {
//			return getNormal(q, i, 6);
//		}
//		else {
//			v1 = getVertex(i, j + 1, f) - getVertex(i, j - 1, f);
//			v2 = getVertex(i + 1, j, f) - getVertex(i - 1, j, f);
//		}
//	}
//	float4 n = (v1.cross3(v2)).normalized();
//	return n;
//}
//
//
//vector<float4> ICQ::read_dem(const string &path) {
//	ifstream infile;
//	infile.open(path);
//	if (!infile) {
//		cout << "open file error!" << endl;
//	}
//	int num_point;
//	infile >> num_point;
//	vector<float4> dem;
//	for (int i = 0; i < num_point; i++) {
//		float4 tmp;
//		for (int j = 0; j < 3; j++) {
//			infile >> tmp[j];
//		}
//		dem.push_back(tmp);
//	}
//	return dem;
//}
//
//float ICQ::getIntersectHeight(int i, int j, int f,int size, double *DEM, Lmap &lmap, float Dstep) {
//	float4 np = getNormal(i, j, f);
//	float4 vp = getVertex(i, j, f);
//	double t = (lmap.landmark.dot(lmap.landmark) - lmap.landmark.dot(vp))/(lmap.landmark.dot(np));
//	//求出与DEM交点
//	Vector3d inter_point = t*np + vp;
//	//求出交点在局部坐标系的坐标
//	double res_x = lmap.u1.dot(inter_point - lmap.landmark);
//	double res_y = lmap.u2.dot(inter_point - lmap.landmark);
//	double res_z = lmap.u3.dot(inter_point - lmap.landmark);
//	cout << "res_x=" << res_x << endl;
//	cout << "res_y=" << res_y << endl;
//	cout << "res_z=" << res_z << endl;
//
//	int i_index = size / 2 - res_y / Dstep;
//	int j_index = size / 2 + res_x / Dstep;
//	cout << "i_index=" << i_index << endl;
//	cout << "j_index=" << j_index << endl;
//
//	if (i_index < 0 || i_index >= size || j_index < 0 || j_index >= size) {
//		return 0;
//	}
//	res_z = lmap.ptr->DEM[i_index * DEMsize + j_index];
//	float4 res_w = lmap.landmark + res_x * lmap.u1 + res_y * lmap.u2 + res_z * lmap.u3;
//	return Dist(vp, res_w);
//
//	/*Vector3d n1 = getNormal(i, j, f);
//	cout << "法向量是：" << n1(0) << "," << n1(1) << "," << n1(2) << endl;
//	Vector3d n2 = Cdem.V.normalized();
//	double t = (n2.dot(Cdem.V)) / (n2.dot(n1));
//
//	cout << "t=" << t << endl;
//	Vector3d res = t*n1;
//	cout << "与DEM的交点是：" << res(0) << "," << res(1) << "," << res(2) << endl;
//	//res只是法线n1与DEM平面的交点，实际上DEM还有高度，应该求得DEM在res处的高度
//	//先将res转换到局部坐标系得到其对应的行列号，得到DEM[i][j]
//	//得到res带有第三维信息的坐标，再转换到世界坐标系下
//	double res_x = Cdem.u1.dot(res - Cdem.V);
//	double res_y = Cdem.u2.dot(res - Cdem.V);
//	double res_z = Cdem.u3.dot(res - Cdem.V);
//	cout << "res_x=" << res_x << endl;
//	cout << "res_y=" << res_y << endl;
//	cout << "res_z=" << res_z << endl;
//
//	int i_index = size / 2 - res_y / Dstep;
//	int j_index = size / 2 + res_x / Dstep;
//	cout << "i_index=" << i_index << endl;
//	cout << "j_index=" << j_index << endl;
//
//	if (i_index < 0 || i_index >= size || j_index < 0 || j_index >= size) {
//		return 0;
//	}
//	res_z = DEM[i_index][j_index];
//	Vector3d res_w = Cdem.V + res_x*Cdem.u1 + res_y*Cdem.u2 + res_z*Cdem.u3;
//
//	Vector3d vp = getVertex(i, j, f);
//	return getDistance(vp, res_w);*/
//}
//
//typedef pair<float, float4> PAIR;
//vector<PAIR> all_min_theta;
//
///*struct compare_by_value {
//bool operator()(PAIR &p1, PAIR &p2) {
//return p1.second < p2.second;
//}
//};*/
//
//
//float ICQ::getIntersectH(int i, int j, int f, string prefix) {
//	//(i,j,f)对应模型上一点，file_path：保存DEM三维点的文件路径
//	float4 np = getNormal(i, j, f);//对应(i,j,f)的法向量
//	float4 vp = getVertex(i, j, f);
//	double thres_angle = 5;
//	map<float, float4> angles;//保存dem的每一点与np的夹角:(dem索引，夹角)
//	for (int i = 1; i <= 50; i++) {//总共50个dem
//
//		string file_path = prefix + to_string(i) + ".txt";//第i个dem对应的路径
//		vector<float4> dem = read_dem(file_path);//
//
//		float min_theta = 180;
//		float4 min_point;
//		for (int j = 0; j < dem.size(); j++) {
//			float4 tmp = dem[j] - vp;
//			//np与dem[j]-vp的夹角
//			float tt = tmp.dot(np) / (sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]));
//			float theta = acos(tt) * 180 / PI; //acos返回的是弧度！！！！
//			if (theta < min_theta) {
//				min_theta = theta;
//				min_point = dem[j];
//			}
//		}
//		angles[min_theta] = min_point;
//	}
//	//将angles按照value从小到大排序
//	//vector<PAIR> vec(angles.begin(), angles.end());
//	//sort(vec.begin(), vec.end());
//	all_min_theta.push_back(*angles.begin());
//	cout << (*angles.begin()).first << endl;
//	if ((*angles.begin()).first <= thres_angle) {//若最小夹角满足阈值要求，则认为是np与该dem的交点
//		return Dist(vp, (*angles.begin()).second);
//	}
//	return 0;
//}
//
//void ICQ::updateVp(int size, double *DEM, Lmap &lmap, float Dstep) {
//	for (int i = 1; i < q; i++) {
//		for (int j = 1; j < q; j++) {
//			for (int f = 1; f < 7; f++) {
//				// cout << getVertexLabel(i, j, f) << endl;
//				int index = getVertexLabel(i, j, f);
//				cout << index << endl;
//				float4 &vp = vectorList[index];
//				float4 np = getNormal(i, j, f);
//				double h = getIntersectHeight(i, j, f, size, DEM, lmap, Dstep);
//				cout << "h=:" << h << endl;
//				vp += np*h;
//			}
//		}
//	}
//}
//
//void ICQ::updateVp(string prefix) {
//	vector<float> height_vec;
//	for (int f = 1; f < 7; f++) {
//		for (int j = 0; j <= q; j++) {
//			for (int i = 0; i <= q; i++) {
//				string path = prefix + "dem";
//				double h = getIntersectH(i, j, f, path);
//				height_vec.push_back(h);
//			}
//		}
//	}
//	for (int f = 1; f < 7; f++) {
//		for (int j = 0; j <= q; j++) {
//			for (int i = 0; i <= q; i++) {
//				int index = getVertexLabel(i, j, f);
//				float4 &vp = vectorList[index];
//				float4 np = getNormal(i, j, f);
//				float h = height_vec[index];
//				vp += np*h;
//			}
//		}
//	}
//}
//
///*double Model::getIntersectHeight(int i, int j, int f, Vector3d point,
//	int size, double **DEM, double Dstep) {
//	Vector3d n1 = getNormal(i, j, f);
//	Vector3d n2 = point.normalized();
//	double t = (n2.dot(point)) / (n2.dot(n1));
//	Vector3d result = t*n1;
//	int x = floor(size / 2 - (result(1) - point(1)) / Dstep);
//	int y = floor((result(0) - point(0)) / Dstep + size / 2);
//	if (x<0 || x >= size || y<0 || y >= size) {
//		return 0;
//	}
//	else {
//		Vector3d vp = getVertex(i, j, f);
//		Vector3d v_inter(result(0), result(1), point(2) + DEM[x][y]);
//		return getDistance(vp, v_inter);
//		//return DEM[x][y];
//	}
//}*/
//
///*void Model::updateVp(vector<Vector3d> pointList, int size,
//	vector<double**> DEMlist, double Dstep, int scale) {
//	for (int i = 1; i<q; i++) {
//		for (int j = 1; j<q; j++) {
//			for (int f = 1; f<7; f++) {
//				int index = getVertexLabel(i, j, f);
//				Vector3d &vp = vectorList[index];//这里要传递引用，否则返回的是拷贝，并没有修改原数组！！！！！！！！
//				Vector3d np = getNormal(i, j, f);
//				double h = 0;
//				int count = DEMlist.size();
//				for (int k = 0; k<count; k++) {
//					h += getIntersectHeight(i, j, f, pointList[k], size, DEMlist[k], Dstep);
//				}
//				cout << "h:" + to_string(h) << endl;
//				double hp = h / count;
//				vp += np*hp*scale;
//			}
//		}
//	}
//}*/
//
//void ICQ::saveModel(const char * DistName) {
//	FILE *fp;
//	fp = fopen(DistName, "w");
//	fprintf(fp, "%s", "ply\n");
//	fprintf(fp, "%s", "format ascii 1.0\n");
//	fprintf(fp, "%s", "comment created by MATLAB ply_write\n");
//	fprintf(fp, "%s %d%s", "element vertex", int(vectorList.size()), "\n");
//	fprintf(fp, "%s", "property float x\n");
//	fprintf(fp, "%s", "property float y\n");
//	fprintf(fp, "%s", "property float z\n");
//	fprintf(fp, "%s", "end_header\n");
//
//	for (int i = 0; i < vectorList.size(); i++) {
//		//cout << vectorList[i](0)<<" " << vectorList[i](1)<<" " << vectorList[i](2) << endl;
//		fprintf(fp, "%f %f %f", vectorList[i](0), vectorList[i](1), vectorList[i](2));
//		fprintf(fp, "\n");
//	}
//
//	fclose(fp);
//	cout << "Write ply Complete!" << endl;
//}
//
//void ICQ::saveModel_txt(const char * DistName) {
//	FILE *fp;
//	fp = fopen(DistName, "w");
//	for (int i = 0; i<vectorList.size(); i++) {
//		//cout << vectorList[i](0)<<" " << vectorList[i](1)<<" " << vectorList[i](2) << endl;
//		fprintf(fp, "%f %f %f", vectorList[i](0), vectorList[i](1), vectorList[i](2));
//		fprintf(fp, "\n");
//	}
//	fclose(fp);
//	cout << "Write txt Complete!" << endl;
//}
//
//void ICQ::densify_model() {
//	vector<float4> densifiedVertexList;
//	for (int f = 1; f <= 6; f++) {
//		for (int j = 0; j <= 2 * q; j++) {
//			for (int i = 0; i <= 2 * q; i++) {
//				if (i % 2 == 0 && j % 2 == 0){
//					densifiedVertexList.push_back(vectorList[getVertexLabel(i / 2, j / 2, f)]);
//				}
//				else if (i % 2 == 1 && j % 2 == 0) {
//					float4 left = vectorList[getVertexLabel(i / 2, j / 2, f)];
//					float4 right = vectorList[getVertexLabel(i / 2 + 1, j / 2, f)];
//					densifiedVertexList.push_back((left + right) / 2);
//				}
//				else if (i % 2 == 0 && j % 2 == 1) {
//					float4 up = vectorList[getVertexLabel(i / 2, j / 2, f)];
//					float4 down = vectorList[getVertexLabel(i / 2, j / 2 + 1, f)];
//					densifiedVertexList.push_back((up + down) / 2);
//				}
//				else {
//					float4 left_up = vectorList[getVertexLabel(i / 2, j / 2, f)];
//					float4 left_down = vectorList[getVertexLabel(i / 2, j / 2 + 1, f)];
//					float4 right_up = vectorList[getVertexLabel(i / 2 + 1, j / 2, f)];
//					float4 right_down = vectorList[getVertexLabel(i / 2 + 1, j / 2 + 1, f)];
//					densifiedVertexList.push_back((left_up+ left_down+ right_up + right_down) / 4);
//				}
//			}
//		}
//	}
//	q = q * 2;
//	vectorList = densifiedVertexList;
//}
//
//void saveMinTheta(string prefix) {
//	FILE *fp;
//	string path = prefix + "min_theta.txt";
//	fp = fopen(path.c_str(), "w");
//	fprintf(fp, "min_theta\tdem_point\n");
//	for (int i = 0; i < all_min_theta.size(); i++) {
//		fprintf(fp, "%f\t(%f,%f,%f)\n",all_min_theta[i].first, all_min_theta[i].second[0], all_min_theta[i].second[1], all_min_theta[i].second[2]);
//	}
//	fclose(fp);
//}
//
//float ICQ::vp_change_aver(string preVerPath) {
//	vector<float4> vertices =  getVectorList();
//	ifstream infile;
//	infile.open(preVerPath.c_str());
//	if (!infile)
//		cout << "open file error!" << endl;
//	float sum_dist = 0.0;
//	int count = 0;
//	while (!infile.eof()) {
//		float4 tmp;
//		for (int j = 0; j < 3; j++) {
//			infile >> tmp[j];
//		}
//		//cout << tmp << endl;
//		sum_dist += Dist(tmp, getVectorList()[count++]);
//	}
//	return sum_dist / count;
//}
//
//void ICQ::read_model(string path) {
//	ifstream infile;
//	infile.open(path.c_str());
//	if (!infile) {
//		cout << "open file error!\n";
//	}
//	int count = 0;
//	while (!infile.eof()) {
//		float4 tmp;
//		for (int j = 0; j < 3; j++) {
//			infile >> tmp[j];
//		}
//		vectorList[count++] = tmp;
//	}
//}
//
//
////int main() {
////
////	float4 center(0, 0, 0, 0);
////	int q = 32;
////	double sideL = 32;
////
////	//初始选择q=32
////	ICQ model(center, q, sideL);
////
////	//q为32时，需要多次更新
////	double thres = 1;
////	string prefix = "E:\\小行星\\ICQ方法201910\\ICQDEM\\Version1\\";
////	string ply_path = prefix + "model" + to_string(model.get_q()) + "_0.ply";
////	string txt_path = prefix + "model" + to_string(model.get_q()) + "_0.txt";
////	model.saveModel(ply_path.c_str());
////	model.saveModel_txt(txt_path.c_str());
////
////	string preVerPath;
////	int count = 1;
////	do {
////		preVerPath = txt_path;
////		model.updateVp(prefix);
////		saveMinTheta(prefix);
////		ply_path = prefix + "model" + to_string(model.get_q()) + "_" + to_string(count) + ".ply";
////		model.saveModel(ply_path.c_str());
////		txt_path = prefix + "model" + to_string(model.get_q()) + "_" + to_string(count) + ".txt";
////		model.saveModel_txt(txt_path.c_str());
////		count++;
////	} while (model.vp_change_aver(preVerPath) <= thres);
////
////	//接着执行将q加倍
////	string prefix = "E:\\小行星\\ICQ方法201910\\ICQDEM\\Version3\\";
////	cout << "读取q=32的模型" << endl;
////	model.read_model("E:\\小行星\\ICQ方法201910\\ICQDEM\\Version3\\model32_1.txt");
////	cout << "读取结束\n";
////	cout << "开始加密模型" << endl;
////	model.densify_model();
////	cout << "模型加密结束\n";
////	string ply_path = prefix + "model" + to_string(model.get_q()) + "_0.ply";
////	string txt_path = prefix + "model" + to_string(model.get_q()) + "_0.txt";
////	cout << "将加密后的模型保持成点云格式...." << endl;
////	model.saveModel(ply_path.c_str());
////	cout << "将加密的模型保存成txt格式...." << endl;
////	model.saveModel_txt(txt_path.c_str());
////
////	string dem_prefix = "E:\\小行星\\ICQ方法201910\\ICQDEM\\";
////	cout << "开始更新加密后的模型...." << endl;
////	model.updateVp(dem_prefix);
////	cout << "更新结束....\n";
////	saveMinTheta(prefix);
////	cout << "开始保存模型....\n";
////	ply_path = prefix + "model" + to_string(model.get_q()) + "_" + "1.ply";
////	model.saveModel(ply_path.c_str());
////	txt_path = prefix + "model" + to_string(model.get_q()) + "_" + "1.txt";
////	model.saveModel_txt(txt_path.c_str());
////	return 0;
////}
//
//
