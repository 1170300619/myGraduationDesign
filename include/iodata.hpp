
#ifndef RENDERER_IODATA_HPP
#define RENDERER_IODATA_HPP

#include "State.hpp"
#include "RandomPoint.hpp"
#include <opencv2/opencv.hpp>
#include "AnimationConfig.hpp"

#define GAMMA_VALUE 1/2.2
#define Z_GAMMA_VALUE 1/2.2


class iodata
{
public:
    /**
     *
     * @param file_name
     */
    static void read_model(const std::string &file_name);

    /**
     * Read ply file data.
     * @param file_name ply file name
     * @return mesh
     */
    static void readPly(const std::string &modelPath, Model *model);

	/**
	* Read obj file data.
	* @param file_name obj file name
	* @return mesh
	*/
	static void read_obj(const std::string &file_name, Model *model);

	/**
	* 插值选取landmark点坐标
	* FullName:  iodata::chooseLandmarks_V2
	* Access:    public static 
	* Returns:   void
	* Param: Model * model 模型信息
	* Param: const string & landmarksPath 路径信息
	* Param: float coverRatio 相邻lmap区域重叠率
	* Param: Camera * camera 相机参数
	*/
	static void writeLandmarks(const vector<float4>& landmarks, const string& path, int mode);

	static void chooseLandmarksFromSpherome(Model *model, vector<float4>& randomPoints, vector<float4>& Landmarks);

	static void  generateWellDistributedPoints(vector<float4>& randomPoints, const int& num);
    /**
     * Read DEM data, the vertex number must be over or equal 2x2.
     * @param file_name DEM file name
     * @param model
     * @return mesh
     */
    static void readDEMWithoutEdge(const std::string &file_name, Model *model);

    static void readLmap(const std::string& path, Lmap& lmap);

    /**
     * Load config file and initialize the state.
     * @param config
     * @return
     */
    static State *load_config(const std::string &config);

	/**
	* 读取观测图像
	* FullName:  iodata::readImages
	* Access:    public static
	* Returns:   void
	* Param: String pattern 图像格式
	* Param: vector<cv::Mat> & images 图像容器
	*/
	static void readImages(const cv::String& pattern, vector<cv::Mat> &images);
    /**
     *
     * @param mesh
     * @param file_name
     */
    static void write_ply(Model *model, const std::string &file_name);

    /**
     *
     */
    void clear()
    {
    }

    /**
     * Write result depth image. The gray level manifests the depth, 0(black) means farthest and 255(white) means nearest.
     * @param depth_buffer
     * @param c
     */
    static void write_depth_image(const float *depth_buffer, Camera *c);

    /**
     * Write result depth image. The gray level manifests the depth, 0(black) means farthest and 255(white) means nearest.
     * @param light
     */
    static void write_depth_image(SunLight *light);

    /**
     *
     * @param frame frame buffer
     */
    static void write_result_image(const FrameBuffer &frame);

	/**
	* 读取landmark信息
	* FullName:  iodata::readLandmarks
	* Access:    public static 
	* Returns:   void
	* Param: const string & path 路径信息
	* Param: vector<float4> & landmarks landmark位置容器
	*/
	static void readLandmarks(const string &path, vector<float4> &landmarks);

    /**
     *
     * @param file_name
     * @param frame
     */
    static void write_result_image(const std::string& file_name, const FrameBuffer &frame);

	//static void calcCorreponse(vector<float4> &camPos, vector<float4>_focal_center, vector<float3> &landmarks, vector<float4>_up, const string &resPath);

	/**
	* 读取相关性文件信息
	* FullName:  iodata::readCorrespnse
	* Access:    public static 
	* Returns:   void
	* Param: const string & path 路径信息
	* Param: vector<vector<int> > & res 相关性容器
	*/
	static void readCamAndSunMes(AnimationConfig* animationConfig, const vector<float4>& landmarks, const string& filePath, State* s);

	static void readCorrespnse(const string&path, vector<vector<int> > &res);

    static void writeImage(const float* orthoDEM, const string& path);

    static void writeUcharImage(const uchar* orthoDEM, const string& path);

    static void addNoiseToCamPos(vector<AnimationConfig*>& configVec);

    static void write3DPoints(const vector<float4>& points, const string& path);

private:
};

#endif
