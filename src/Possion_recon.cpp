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
//		cout << "�������ݶ�ȡʧ�ܣ�" << endl;
//	}
//	std::cout << "Orginal points number: " << cloud->points.size() << std::endl;
//
//
//	// ͳ���˲�
//	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> statisOutlierRemoval;       //�����˲�������
//	statisOutlierRemoval.setInputCloud(cloud);            //���ô��˲��ĵ���
//	statisOutlierRemoval.setMeanK(50);                                //�����ڽ���ͳ��ʱ���ǲ�ѯ���ٽ�����
//	statisOutlierRemoval.setStddevMulThresh(3.0);                     //�����ж��Ƿ�Ϊ��Ⱥ��ķ�ֵ:��ֵ+1.0*��׼��
//	statisOutlierRemoval.filter(*cloud_filtered);                     //�˲�����洢��cloud_filtered
//	pcl::io::savePLYFile("E:\\SPC\\filtered.ply", *cloud_filtered);
//
//	// �Ե����ز���  
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr treeSampling(new pcl::search::KdTree<pcl::PointXYZ>);
//	pcl::PointCloud<pcl::PointXYZ> mls_point; 
//	pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointXYZ> mls; // ������С����ʵ�ֵĶ���mls
//	mls.setComputeNormals(false);  //��������С���˼������Ƿ���Ҫ�洢����ķ���
//	mls.setInputCloud(cloud_filtered);         //���ô��������
//	mls.setPolynomialOrder(3);            // ���3�׶���ʽ���
//	mls.setPolynomialFit(false);// ����Ϊfalse���Լ���smooth
//	mls.setSearchMethod(treeSampling);         
//	mls.setSearchRadius(0.05);         
//	mls.process(mls_point);     
//	cloud_smoothed = mls_point.makeShared();
//	std::cout << "cloud_smoothed: " << cloud_smoothed->size() << std::endl;
//	//save cloud_smoothed
//	pcl::io::savePLYFileASCII("E:\\SPC\\cloud_smoothed.ply", *cloud_smoothed);
//
//	// ���߹���
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;         
//	normalEstimation.setInputCloud(cloud_smoothed);               
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
//	normalEstimation.setSearchMethod(tree);
//	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>); 																		
//	normalEstimation.setKSearch(30);
//	//normalEstimation.setRadiusSearch(0.03);//����ÿһ���㶼�ð뾶Ϊ3cm�Ľ���������ʽ
//	normalEstimation.compute(*normals);
//	std::cout << "normals: " << normals->size() << ", " << "normals fields: " << pcl::getFieldsList(*normals) << std::endl;
//	pcl::io::savePLYFileASCII("E:\\SPC\\normal.ply", *normals);
//
//
//	// ������λ�ˡ���ɫ��������Ϣ���ӵ�һ��
//	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>);
//	pcl::concatenateFields(*cloud_smoothed, *normals, *cloud_with_normals);
//	pcl::io::savePLYFileASCII("E:\\SPC\\cloud_with_normals.ply", *cloud_with_normals);
//
//	// ̰��ͶӰ���ǻ�
//	pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(new pcl::search::KdTree<pcl::PointNormal>);
//	tree2->setInputCloud(cloud_with_normals);
//
//	// ���ǻ�
//	pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;   // �������ǻ�����
//	pcl::PolygonMesh triangles; //�洢�������ǻ�������ģ��
//
//	// �������ǻ�����
//	gp3.setSearchRadius(1);  //��������ʱ�İ뾶��Ҳ����KNN����뾶
//	gp3.setMu(2.5);  //������������������ڵ����Զ����Ϊ2.5��������ֵ2.5-3��������ʹ���㷨����Ӧ�����ܶȵı仯
//	gp3.setMaximumNearestNeighbors(200); 
//
//	gp3.setMinimumAngle(M_PI / 18); // ���ǻ�����ڽǵ���С�Ƕ�
//	gp3.setMaximumAngle(2 * M_PI / 3); // ���ǻ����ڽǵ����Ƕ�
//
//	gp3.setMaximumSurfaceAngle(M_PI / 4); // ����ĳ�㷨�߷���ƫ�������㷨�ߵ����Ƕ�45�㣬�������������ʱ�����Ǹõ�
//	gp3.setNormalConsistency(false);  //���øò���Ϊtrue��֤���߳���һ�£�����Ϊfalse�Ļ�������з���һ���Լ��
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