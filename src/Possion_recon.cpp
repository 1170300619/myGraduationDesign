//#pragma warning(disable:4996);
//#include <pcl/point_types.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/io/ply_io.h>
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/filters/radius_outlier_removal.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/statistical_outlier_removal.h>
//#include <pcl/surface/mls.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/surface/gp3.h>
//#include <pcl/surface/poisson.h>
//#include "Possion_recon.hpp"
//using namespace std;
//
//Possion_recon::Possion_recon()
//{}
//
//Possion_recon::~Possion_recon() = default;
//
//void Possion_recon::possion_recon(std::string filename) {
//	// Load input file
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_downSampled(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_smoothed(new pcl::PointCloud<pcl::PointXYZ>);
//	if (pcl::io::loadPLYFile("E:\\SPC\\lmap_0.ply", *cloud) == -1)
//	{
//		cout << "点云数据读取失败！" << endl;
//	}
//	std::cout << "Orginal points number: " << cloud->points.size() << std::endl;
//
//
//	// 统计滤波
//	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> statisOutlierRemoval;       //创建滤波器对象
//	statisOutlierRemoval.setInputCloud(cloud);            //设置待滤波的点云
//	statisOutlierRemoval.setMeanK(50);                                //设置在进行统计时考虑查询点临近点数
//	statisOutlierRemoval.setStddevMulThresh(3.0);                     //设置判断是否为离群点的阀值:均值+1.0*标准差
//	statisOutlierRemoval.filter(*cloud_filtered);                     //滤波结果存储到cloud_filtered
//	pcl::io::savePLYFile("E:\\SPC\\filtered.ply", *cloud_filtered);
//
//	// 对点云重采样  
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr treeSampling(new pcl::search::KdTree<pcl::PointXYZ>);
//	pcl::PointCloud<pcl::PointXYZ> mls_point; 
//	pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointXYZ> mls; // 定义最小二乘实现的对象mls
//	mls.setComputeNormals(false);  //设置在最小二乘计算中是否需要存储计算的法线
//	mls.setInputCloud(cloud_filtered);         //设置待处理点云
//	mls.setPolynomialOrder(3);            // 拟合3阶多项式拟合
//	mls.setPolynomialFit(false);// 设置为false可以加速smooth
//	mls.setSearchMethod(treeSampling);         
//	mls.setSearchRadius(0.05);         
//	mls.process(mls_point);     
//	cloud_smoothed = mls_point.makeShared();
//	std::cout << "cloud_smoothed: " << cloud_smoothed->size() << std::endl;
//	//save cloud_smoothed
//	pcl::io::savePLYFileASCII("E:\\SPC\\cloud_smoothed.ply", *cloud_smoothed);
//
//	// 法线估计
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;         
//	normalEstimation.setInputCloud(cloud_smoothed);               
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
//	normalEstimation.setSearchMethod(tree);
//	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>); 																		
//	normalEstimation.setKSearch(30);
//	//normalEstimation.setRadiusSearch(0.03);//对于每一个点都用半径为3cm的近邻搜索方式
//	normalEstimation.compute(*normals);
//	std::cout << "normals: " << normals->size() << ", " << "normals fields: " << pcl::getFieldsList(*normals) << std::endl;
//	pcl::io::savePLYFileASCII("E:\\SPC\\normal.ply", *normals);
//
//
//	// 将点云位姿、颜色、法线信息连接到一起
//	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>);
//	pcl::concatenateFields(*cloud_smoothed, *normals, *cloud_with_normals);
//	pcl::io::savePLYFileASCII("E:\\SPC\\cloud_with_normals.ply", *cloud_with_normals);
//
//	// 贪心投影三角化
//	pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(new pcl::search::KdTree<pcl::PointNormal>);
//	tree2->setInputCloud(cloud_with_normals);
//
//	// 三角化
//	pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;   // 定义三角化对象
//	pcl::PolygonMesh triangles; //存储最终三角化的网络模型
//
//	// 设置三角化参数
//	gp3.setSearchRadius(1);  //设置搜索时的半径，也就是KNN的球半径
//	gp3.setMu(2.5);  //设置样本点搜索其近邻点的最远距离为2.5倍（典型值2.5-3），这样使得算法自适应点云密度的变化
//	gp3.setMaximumNearestNeighbors(200); 
//
//	gp3.setMinimumAngle(M_PI / 18); // 三角化后得内角的最小角度
//	gp3.setMaximumAngle(2 * M_PI / 3); // 三角化后内角的最大角度
//
//	gp3.setMaximumSurfaceAngle(M_PI / 4); // 设置某点法线方向偏离样本点法线的最大角度45°，如果超过，连接时不考虑该点
//	gp3.setNormalConsistency(false);  //设置该参数为true保证法线朝向一致，设置为false的话不会进行法线一致性检查
//
//	gp3.setInputCloud(cloud_with_normals);
//	gp3.setSearchMethod(tree2);
//	gp3.reconstruct(triangles);
//
//	pcl::io::savePLYFile("result.ply", triangles);
//}
//
////int main() {
////	string filename = "E:\\SPC\\lmap_0.ply";
////	Possion_recon recon;
////	recon.possion_recon(filename);
////}