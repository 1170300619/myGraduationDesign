#pragma warning(disable:4996);
#include "../include/iodata.hpp"
#include <string>
#include "../include/PhongMaterial.hpp"
#include "opencv2/imgproc.hpp"

using namespace std;
using namespace Render;
using namespace cv;

void iodata::readDEMWithoutEdge(const string &file_name, Model *model)
{
	int demSize;
	float4 landmark;
	float4 u1, u2, u3;

	ifstream in;
	in.open(file_name.c_str());
	in >> demSize >> demSize;
	in >> landmark[0] >> landmark[1] >> landmark[2] >> landmark[3];
	in >> u1[0] >> u1[1] >> u1[2] >> u1[3];
	in >> u2[0] >> u2[1] >> u2[2] >> u2[3];
	in >> u3[0] >> u3[1] >> u3[2] >> u3[3];

	auto *mesh = new Mesh();
	float _x, _y, _z, _w;
	for (int i = 0; i < demSize; i++) {
		for (int j = 0; j < demSize; j++) {
			in >> _x >> _y >> _z >> _w;
            if((i == 0) || (i == demSize - 1) || (j == 0) || (j == demSize - 1))continue;
			float4 c(_x, _y, _z, _w);
            float4 ver = localToBody(u1, u2, u3, landmark, c); // 将 Lmap 的局部坐标转为本体坐标
			Vertex v(ver);
			mesh->add_vertex(v);
		}
	}

    for (int i = 0; i < demSize - 3; i++)
    {
        for (int j = 0; j < demSize - 3; j++)
        {
            Triangle triangle_0(mesh->vertices, i * (demSize - 2) + j, i * (demSize - 2) + demSize - 2 + j, i * (demSize - 2) + j + 1);
            Triangle triangle_1(mesh->vertices, i * (demSize - 2) + j + 1, i * (demSize - 2) + demSize - 2 + j, i * (demSize - 2) + demSize - 2 + j + 1);
            mesh->add_triangle(triangle_0);
            mesh->add_triangle(triangle_1);
        }
    }
	in.close();

    cout << "vertex number: " << mesh->numVertices<< endl;
    cout << "triangles number: " << mesh->numTriangles << endl;

    model->meshes.push_back(mesh);
}


void iodata::readLmap(const std::string &path, Lmap &lmap) {
    int demSize;

    ifstream in;
    in.open(path.c_str());
    in >> demSize >> demSize;
    in >> lmap.landmark[0] >> lmap.landmark[1] >> lmap.landmark[2] >> lmap.landmark[3];
    in >> lmap.u1[0] >> lmap.u1[1] >> lmap.u1[2] >> lmap.u1[3];
    in >> lmap.u2[0] >> lmap.u2[1] >> lmap.u2[2] >> lmap.u2[3];
    in >> lmap.u3[0] >> lmap.u3[1] >> lmap.u3[2] >> lmap.u3[3];

    for(int i = 0; i < demSize * demSize; ++i){
        in >> lmap.ptr->DEM[i][0] >> lmap.ptr->DEM[i][1] >> lmap.ptr->DEM[i][2] >> lmap.ptr->DEM[i][3];
    }
    for(int i = 0; i < demSize * demSize; ++i){
        in >> lmap.ptr->albedoMap[i];
    }
}

State *iodata::load_config(const string &config)
{
    auto *s = new State();

    ifstream in;
    int pixel_x, pixel_y;
    float ccd_size_x, ccd_size_y, focal;
    float camera_center_x, camera_center_y, camera_center_z;
    float focal_center_x, focal_center_y, focal_center_z;
    float up_x, up_y, up_z;

    in.open(config.c_str());
    if (!in.is_open())
    {
        cerr << "FAIL TO OPEN FILE" << endl;
        throw exception();
    }

    // load camera message
    in >> pixel_x >> pixel_y >> ccd_size_x >> ccd_size_y >> focal;
    in >> camera_center_x >> camera_center_y >> camera_center_z;
    in >> focal_center_x >> focal_center_y >> focal_center_z;
    in >> up_x >> up_y >> up_z;

    auto c = new Camera();
    float4 camera_center(camera_center_x / 100.0f, camera_center_y / 100.0f, camera_center_z / 100.0f, 1.0f);
    float4 focal_center(focal_center_x / 100.0f, focal_center_y / 100.0f, focal_center_z / 100.0f, 1.0f);
    float4 up(up_x, up_y, up_z, 0);
    c->set_viewport(pixel_x, pixel_y, ccd_size_x, ccd_size_y, focal);
    c->set_look_at(camera_center, focal_center, up);
    s->camera = c;

    // load light source parameters
    string type;
    int shadow_x, shadow_y;
    float range_x, luminance;
    float center_x, center_y, center_z;

    while (in >> type)
    {
        if (type == "POINT")
        {
            cout << "light: POINT" << endl;
            in >> shadow_x; // point light case: x == y
            in >> luminance;
            in >> center_x >> center_y >> center_z;
            auto point_light = new PointLight(luminance, shadow_x, center_x, center_y, center_z);
//            s->light.push_back(point_light);
            continue;
        } else if (type == "SUN")
        {
            cout << "light: SUN" << endl;
            in >> shadow_x >> shadow_y >> range_x;
            in >> luminance;
            in >> center_x >> center_y >> center_z;
            in >> focal_center_x >> focal_center_y >> focal_center_z;
            in >> up_x >> up_y >> up_z;

            auto sun_light = new SunLight();
            float4 center(center_x / 100.0f, center_y / 100.0f, center_z / 100.0f, 1.0f);
            focal_center << focal_center_x / 100.0f, focal_center_y / 100.0f, focal_center_z / 100.0f, 1.0f;
            up << up_x, up_y, up_z, 0;

            sun_light->set_luminance(luminance);
            sun_light->set_viewport(shadow_x, shadow_y, range_x);
            sun_light->set_look_at(center, focal_center, up);
//            s->light.push_back(sun_light);
            continue;
        } else if (type == "SHADOW")
        {
            cout << "SHADOW ";
            in >> type;
            if (type == "ON")
            {
                cout << "ON" << endl;
                s->shadow = SHADOW;
                continue;
            } else if (type == "OFF")
            {
                cout << "OFF" << endl;
                s->shadow = NO_SHADOW;
                continue;
            } else
            {
                cerr << "SHADOW TYPE ERROR: " + type << endl;
                throw exception();
            }
        } else if (type == "CULL_MODE")
        {
            cout << "CULL_MODE ";
            in >> type;
            if (type == "BACK")
            {
                cout << "BACK" << endl;
                s->face_cull_mode = BACK;
                continue;
            } else if (type == "FRONT")
            {
                cout << "FRONT" << endl;
                s->face_cull_mode = FRONT;
                continue;
            } else if (type == "NONE")
            {
                cout << "CULL_MODE NONE" << endl;
                s->face_cull_mode = NONE;
                continue;
            } else
            {
                cerr << "CULL MODE TYPE ERROR: " + type << endl;
                throw exception();
            }
        } else if (type == "TEXTURE_TYPE")
        {
            cout << "TEXTURE_TYPE ";
            in >> type;
            if (type == "NORMAL")
            {
                cout << "NORMAL" << endl;
                s->texture_type = NORMAL_TEXTURE;
                continue;
            } else if (type == "MIPMAP")
            {
                cout << "MIPMAP" << endl;
                s->texture_type = MIPMAP;
                continue;
            } else
            {
                cerr << "TEXTURE TYPE ERROR: " + type << endl;
                throw exception();
            }
        } else if (type == "TEXTURE_SAMPLER")
        {
            cout << "TEXTURE_SAMPLER ";
            in >> type;
            if (type == "NORMAL")
            {
                cout << "NORMAL" << endl;
                s->sampler = NORMAL;
                continue;
            } else if (type == "BILINEAR")
            {
                cout << "BILINEAR" << endl;
                s->sampler = BILINEAR;
                continue;
            } else if (type == "TRILINEAR")
            {
                cout << "TRILINEAR" << endl;
                s->sampler = TRILINEAR;
                continue;
            } else if (type == "ANISOTROPIC")
            {
                cout << "ANISOTROPIC" << endl;
                s->sampler = ANISOTROPIC;
                continue;
            } else
            {
                cerr << "TEXTURE SAMPLER TYPE ERROR: " + type << endl;
                throw exception();
            }
        } else
        {
            cerr << "TYPE ERROR: " + type << endl;
            throw exception();
        }
    }
    in.close();
    return s;
}

void iodata::read_obj(const string &file_name, Model *model) {
	string obj_file;
	int vertex_attr_size = 0;
	int face_attr_size = 0;
	char buffer[256], str[255];
	float x, y, z;
	int index_0, index_1, index_2;
	auto *mesh = new Mesh();
	ifstream in;

	in.open(file_name.c_str());
	if (!in.is_open()) {
		cerr << "Fail to open data/model.txt file" << endl;
		throw exception();
	}
	in >> obj_file;
	in.close();

	in.open(obj_file.c_str());
	if (!in.is_open()) {
		cerr << "Fail to open obj file" << endl;
		throw exception();
	}

	float4 world, color;
	while (!in.getline(buffer, 255).eof()) {
		buffer[255] = '\0';
		sscanf_s(buffer, "%s", str, 255);

		if (buffer[0] == 'v' && (buffer[1] == ' ' || buffer[1] == 32)) {
			if (sscanf(buffer, "v %f %f %f", &x, &y, &z) == 3) {
				world << x, y, z, 1.0f;
				color << 1.0f, 1.0f, 1.0f, 1.0f;
				Vertex vertex(world, color);
				mesh->add_vertex(vertex);
			}
		}
		else if (buffer[0] == 'f' && (buffer[1] == ' ' || buffer[1] == 32)) {
			if (sscanf(buffer, "f %d %d %d", &index_0, &index_1, &index_2) == 3) {
				Triangle triangle(mesh->vertices, index_0, index_1, index_2);
				mesh->add_triangle(triangle);
			}
		}
	}
	float4 *normal_sum = new float4[mesh->numVertices];
	for (int i = 0; i < mesh->numVertices; i++)
	{
		normal_sum[i] << 0, 0, 0, 0;
	}
	for (int i = 0; i < mesh->numTriangles; ++i) {
		normal_sum[mesh->triangles[i].vertex_0] += mesh->triangles[i].normal;
		normal_sum[mesh->triangles[i].vertex_1] += mesh->triangles[i].normal;
		normal_sum[mesh->triangles[i].vertex_2] += mesh->triangles[i].normal;
	}
	for (int i = 0; i < mesh->numVertices; ++i) {
		mesh->vertices[i].normal = normal_sum[i].normalized();
	}
	in.close();
	delete[] normal_sum;

	cout << "vertex number: " << mesh->numVertices << endl;
	cout << "triangles number: " << mesh->numTriangles << endl;

	// load material source parameters
	in.open(file_name.c_str());
	if (!in.is_open())
	{
		cerr << "FAIL TO OPEN FILE" << endl;
		throw exception();
	}

	string type;
	float ka, kd, ks, spec_rank;
	in >> obj_file;
	in >> type;
	in >> ka >> kd >> ks >> spec_rank;
	cout << "ka: " << ka << " kd: " << kd << " ks: " << ks << endl;
	cout << "specular rank: " << spec_rank << endl;

	Material *m;
	if (type == "FLAT")
	{
//		m = new FlatMaterial();
//		m->ambient = ka;
//		m->diffuse = kd;
//		m->specular = ks;
//		m->spec_rank = spec_rank;
//		mesh->material = m;
	}
	else if (type == "PHONG")
	{
		m = new PhongMaterial();
		m->ambient = ka;
		m->diffuse = kd;
		m->specular = ks;
		m->spec_rank = spec_rank;
		mesh->material = m;
	}
	else
	{
		cerr << "MATERIAL TYPE ERROR" << endl;
		throw exception();
	}

	// load base_texture
	string texture_file;
	while (!in.eof())
	{
		in >> type;
		if (type == "TEXTURE")
		{
			in >> texture_file;
			auto *texture = new Texture2D(texture_file);
			m->set_texture(texture);
			cout << "base_texture: " << texture_file << endl;
		}
		else
		{
			break;
		}
	}

	in.close();
	model->meshes.push_back(mesh);
}

void iodata::readPly(const string &modelPath, Model *model)
{
    int numVertices, numTriangles;
    int vertex_attr_size = 0;
    int face_attr_size = 0;
    ifstream in;

    in.open(modelPath.c_str());
    if (!in.is_open())
    {
        cerr << "FAIL TO OPEN FILE" << endl;
        throw exception();
    }

//    in >> ply_file;
//	cout << "ply_file is " << ply_file << endl;
//    in.close();
//
//    in.open(ply_file.c_str());
//    if (!in.is_open())
//    {
//        cerr << "FAIL TO OPEN FILE" << endl;
//        throw exception();
//    }

    string temp;
    in >> temp;
    while (!in.eof())
    {
        //        cout << temp << endl;
        if (temp == "ply")
        {
            in >> temp;
            continue;
        }
        if (temp == "format")
        {
            getline(in, temp, '\n');
            in >> temp;
            continue;
        }
        if (temp == "comment")
        {
            getline(in, temp, '\n');
            in >> temp;
            continue;
        }
        if (temp == "element")
        {
            in >> temp;
            if (temp == "vertex")
            {
                in >> numVertices;
                while (!in.eof())
                {
                    in >> temp;
                    if (temp == "property")
                    {
                        getline(in, temp, '\n');
                        vertex_attr_size++;
                    } else
                    {
                        break;
                    }
                }
                continue;
            }

            if (temp == "face")
            {
                in >> numTriangles;
                while (!in.eof())
                {
                    in >> temp;
                    if (temp == "property")
                    {
                        getline(in, temp, '\n');
                        face_attr_size++;
                    } else
                    {
                        break;
                    }
                }
                continue;
            }
        }

        if (temp == "end_header")
        {
            break;
        }

        cerr << "FILE FORMAT ERROR" << endl;
        throw exception();
    }

    auto *mesh = new Mesh();
    float x, y, z, u, v;
    int r, g, b, a;
    float4 world, color;
    float2 uv;
    // to calculate vertex normal
    float4 *normal_sum = new float4[numVertices]; // sum of triangle normal vector contains vertex
    mesh->vertices.reserve(numVertices);
    mesh->triangles.reserve(numTriangles);
    for (int i = 0; i < numVertices; i++)
    {
        normal_sum[i] << 0, 0, 0, 0;
    }
	cout << numVertices << endl;
    for (int i = 0; i < numVertices; i++)
    {
        switch (vertex_attr_size)
        {
            case 3:
            {
                in >> x >> y >> z;
                world << x, y, z, 1.0f;
                color << 1.0f, 1.0f, 1.0f, 1.0f;
                Vertex vertex(world, color);
                mesh->add_vertex(vertex);
                break;
            }
            case 5:
            {
                in >> x >> y >> z >> u >> v;
                world << x, y, z, 1.0f;
                color << 1.0f, 1.0f, 1.0f, 1.0f;
                uv << u, 1 - v; // OpenCV read from left-up, to convert the coordinate center to left-down
                Vertex vertex(world, uv);
                mesh->add_vertex(vertex);
                break;
            }
            case 6:
            {
                in >> x >> y >> z >> r >> g >> b;
                world << x, y, z, 1.0f;
                color << (float)r / 255.0f, float(g) / 255.0f, float(b) / 255.0f, 1.0f;
                Vertex vertex(world, color);
                mesh->add_vertex(vertex);
                break;
            }
            case 7:
            {
                in >> x >> y >> z >> r >> g >> b >> a;
                world << x, y, z, 1.0f;
                color << (float)r / 255.0f, float(g) / 255.0f, float(b) / 255.0f, float(a) / 255.0f;
                Vertex vertex(world, color);
                mesh->add_vertex(vertex);
                break;
            }
            default:
                break;
        }
    }
    for (int i = 0; i < numTriangles; i++)
    {
        int num, index_0, index_1, index_2, index_3;
        in >> num;
        if (num == 3)
        {
            in >> index_0 >> index_1 >> index_2;
            Triangle triangle(mesh->vertices, index_0, index_1, index_2);
            mesh->add_triangle(triangle);
            normal_sum[index_0] += triangle.normal;
            normal_sum[index_1] += triangle.normal;
            normal_sum[index_2] += triangle.normal;
        } else if (num == 4)
        {
            in >> index_0 >> index_1 >> index_2 >> index_3;
            Triangle triangle_0(mesh->vertices, index_0, index_1, index_2);
            Triangle triangle_1(mesh->vertices, index_0, index_2, index_3);
            mesh->add_triangle(triangle_0);
            mesh->add_triangle(triangle_1);
            normal_sum[index_0] += triangle_0.normal;
            normal_sum[index_1] += triangle_0.normal;
            normal_sum[index_2] += triangle_0.normal;
            normal_sum[index_0] += triangle_1.normal;
            normal_sum[index_2] += triangle_1.normal;
            normal_sum[index_3] += triangle_1.normal;
        }
    }
    in.close();

    // calculate vertex normal
    for (int i = 0; i < numVertices; i++)
    {
        mesh->vertices[i].normal = normal_sum[i].normalized();
    }
    delete[] normal_sum;

    cout << "vertex number: " << mesh->numVertices << endl;
    cout << "triangles number: " << mesh->numTriangles << endl;

//    // load material source parameters
//    in.open(modelPath.c_str());
//    if (!in.is_open())
//    {
//        cerr << "FAIL TO OPEN FILE" << endl;
//        throw exception();
//    }
//
//    string type;
//    float ka, kd, ks, spec_rank;
//    in >> ply_file;
//    in >> type;
//    in >> ka >> kd >> ks >> spec_rank;
//    cout << "ka: " << ka << " kd: " << kd << " ks: " << ks << endl;
//    cout << "specular rank: " << spec_rank << endl;
//
//    Material *m;
//    if (type == "FLAT")
//    {
////        m = new FlatMaterial();
////        m->ambient = ka;
////        m->diffuse = kd;
////        m->specular = ks;
////        m->spec_rank = spec_rank;
////        mesh->material = m;
//    } else if (type == "PHONG")
//    {
//        m = new PhongMaterial();
//        m->ambient = ka;
//        m->diffuse = kd;
//        m->specular = ks;
//        m->spec_rank = spec_rank;
//        mesh->material = m;
//    } else
//    {
//        cerr << "MATERIAL TYPE ERROR" << endl;
//        throw exception();
//    }
//
//    // load base_texture
//    string texture_file;
//    while (!in.eof())
//    {
//        in >> type;
//        if (type == "TEXTURE")
//        {
//            in >> texture_file;
//            auto *texture = new Texture2D(texture_file);
//            m->set_texture(texture);
//            cout << "base_texture: " << texture_file << endl;
//        } else
//        {
//            break;
//        }
//    }
//    in.close();

    model->meshes.push_back(mesh);
    model->createBVH();
    for(auto &tri : model->meshes[0]->triangles){
        tri.mesh = model->meshes[0];
    }
}

void iodata::generateWellDistributedPoints(vector<float4>& randomPoints, const int& num) {
	int z1 = 20, z2 = 112;
	RandomPoint sita(z1, 0.0f);
	RandomPoint pesi(z2, 0.0f);
	float u, v, r2;
	for (int i = 0; i < num; ++i) {
		sita.calcRandom(pesi.seed);
		pesi.calcRandom(sita.seed);

		u = 2 * sita.random - 1.0;
		v = 2 * pesi.random - 1.0;
		r2 = pow(u, 2) + pow(v, 2);
		if (r2 < 1)
		{
			float4 index(2 * u*sqrt(1 - r2), 2 * v*sqrt(1 - r2), 1 - 2 * r2, 1.0f);
			randomPoints.emplace_back(index);
		}
	}
	ofstream out;
	out.open(kRandomPointsPath.c_str());
	if (!out.is_open()) {
		cerr << "Fail to open randomPointsFilePath" << endl;
		throw exception();
	}
	for (auto i : randomPoints) {
		out << i[0] << " " << i[1] << " " << i[2] <<  " " << i[3] << endl;
	}
	out.close();
}


// 利用射线和三角面相交原理对Landmark初始点的选取优化版本
//void iodata::chooseLandmarks_V2(Model *model, const string& landmarksPath, float coverRatio,Camera *camera) {
//	/*功能：读取obj的model，将其根据点云坐标划分成多个网格，每个网格内选取一个最接近中心点的三维点作为landmark。
//	参数：modelPath读取的模型路径
//	landmarksPath:保存的landmarks路径
//	coverRatio:相邻的landmark覆盖率
//	*/
//	Mesh *mesh = model->meshes[0];
//	int _numVertices = mesh->vertices.size();
//	int _numTriangles = mesh->triangles.size();
//	//去掉模型边缘区域10%
//	float edgeRatio = 0.2;
//	mesh->min_X = mesh->min_X*(1 - edgeRatio);
//	mesh->max_X = mesh->max_X*(1 - edgeRatio);
//	mesh->min_Y = mesh->min_Y*(1 - edgeRatio);
//	mesh->max_Y = mesh->max_Y*(1 - edgeRatio);
//	float length = mesh->max_X - mesh->min_X;
//	float width = mesh->max_Y - mesh->min_Y;
//	float modelLenOfOnePixel = camera->ccd_x * Render::kCamDist / (camera->focal * Render::kImageSize);
//	//3、假设一个像素代表长度len,相邻landmark的间距是(len-len*重叠率)
//	// 得到每个网格内的中心点
//	float intervalBetweenLandmark = Render::kDEMSize * modelLenOfOnePixel *(1 - coverRatio);
//	cout << intervalBetweenLandmark << endl;
//	int numLandmarkX = (length - modelLenOfOnePixel * Render::kDEMSize) / intervalBetweenLandmark + 1;//记录X方向上分配多少个landmark点
//	int numLandmarkY = (width - modelLenOfOnePixel * Render::kDEMSize) / intervalBetweenLandmark + 1;//记录Y方向上分配多少个landmark点
//
//	cout << (numLandmarkX * numLandmarkY) << endl;
//	int size = numLandmarkX * numLandmarkY;
//	std::vector<std::vector<float4> > centerPoints(numLandmarkY, std::vector<float4>(numLandmarkX));
//	for (int i = 0; i < numLandmarkX; i++) {
//		float x = mesh->min_X + intervalBetweenLandmark*i;
//		for (int j = 0; j < numLandmarkY; j++) {
//			float y = mesh->min_Y + intervalBetweenLandmark*j;
//			centerPoints[j][i] = float4(x, y, floor(mesh->min_Z),1.0f);//记录中心点坐标,保证初始Landmark点在点云坐标下方
//		}
//	}
//
//	//找点云文件中距离landmark点最近的点，将其作为landmark点坐标
//	float4 ray = { 0,0,1.0f,0 };//方向设置为沿z轴向上
//	float T[1] = { 0 };//记录参数
//	float4 interPoint;//记录Landmark射线和三角面的交点
//	float maxDist = FLT_MIN;//记录初始Landmark点和交点之间的距离
//	float4 tmp;
//	//储存原始Landmark点坐标（利用最近点原理)
//	/*vector<vector<float4> > landmarks(numLandmarkX, vector<float4>(numLandmarkY, float4(FLT_MIN, FLT_MIN, FLT_MIN,1.0f)));*/
//	vector<float4> landmarks;
//	landmarks.resize(size, float4(FLT_MIN, FLT_MIN, FLT_MIN, 1.0f));
//	for (int j = 0; j < numLandmarkY; j++) {
//		for (int k = 0; k < numLandmarkX; k++) {
//			for (int i = 0; i < _numVertices; i++) {
//				tmp = mesh->vertices[i].position;
//				if (Render::Dist(centerPoints[j][k], tmp) < Render::Dist(landmarks[j * numLandmarkX + k], centerPoints[j][k])) {
//					landmarks[j * numLandmarkX + k] = tmp;
//				}
//			}
//		}
//		cout << j << endl;
//	}
//	//将和三角面有交点的射线产生的Landmark进行更新
//	for (int j = 0; j < numLandmarkY; j++) {
//		for (int k = 0; k < numLandmarkX; k++) {
//			for (int i = 0; i < _numTriangles; i++) {
//				float4 index_normal = mesh->triangles[i].normal;
//				if (index_normal.dot(ray) < Render::kMinIncludedAngle)continue;
//				const Vertex &v0 = mesh->vertices[mesh->triangles[i].vertex_0];
//				const Vertex &v1 = mesh->vertices[mesh->triangles[i].vertex_1];
//				const Vertex &v2 = mesh->vertices[mesh->triangles[i].vertex_2];
//				if (Render::IntersectTriangle(centerPoints[j][k], ray, v0, v1, v2, T)) {
//					interPoint = centerPoints[j][k] + ray * T[0];
//					float tmp = T[0];
//					if (tmp > maxDist) {
//						landmarks[j * numLandmarkX + k] = interPoint;
//						maxDist = tmp;
//					}
//				}
//			}
//			maxDist = FLT_MIN;
//		}
//	}
//
//	//5、将landmarks写入文件
//	FILE *fp;
//	fp = fopen(landmarksPath.c_str(), "w");
//	if (fp == NULL) {
//		cout << "can't open!" << endl;
//	}
//	for (int i = 0; i < size; i++) {
//		fprintf(fp, "%f %f %f %f\n", landmarks[i](0), landmarks[i](1), landmarks[i](2), landmarks[i](3));
//	}
//	fclose(fp);
//}

void iodata::writeLandmarks(const vector<float4>& landmarks, const string& path, int mode){
    std::ofstream out;
    out.open((path + "landmarks_" + to_string(mode) + ".txt").c_str());
    if (!out.is_open()) {
        cerr << "Fail to open landmark filePath" << endl;
        throw exception();
    }
    for (auto & l : landmarks) {
        out << l[0] << " " << l[1] << " " << l[2] << " " << 1.0f;
        if(l != landmarks[landmarks.size() - 1])out << endl;
    }
    out.close();

    out.open((path + "landmarksMesh_" + to_string(mode) + ".pp").c_str());
    if (!out.is_open()) {
        cerr << "Fail to open landmarkMesh filepath" << endl;
        throw exception();
    }
    out << "<!DOCTYPE PickedPoints>" << endl;
    out << "<PickedPoints>" << endl;
    out << "<DocumentData>" << endl;
    out << R"(<DateTime date = "2023 - 02 - 24" time = "19:15 : 54"/>)" << endl;
    out << "<User name = \"admin\"/>" << endl;
    out << "<DataFileName name = \"5x.ply\"/>" << endl;
    out << "<templateName name = \"\"/>" << endl;
    out << "</DocumentData>" << endl;
    for (int i = 0; i < landmarks.size(); ++i) {
        out << "<point x = " << "\"" << landmarks[i][0] << "\"" << " name = \"" << i << "\"" << " y = \"" << landmarks[i][1] << "\"" << " active = \"" << "1\"" <<  " z = \"" << landmarks[i][2] << "\"" << "/>" << endl;
    }
    out << "</PickedPoints>" << endl;
    out.close();
}

void iodata::chooseLandmarksFromSpherome(Model *model, vector<float4>& randomPoints, vector<float4>& landmarks)
{
    HitRecord record;
	//将和三角面有交点的射线产生的Landmark进行更新

	for (auto & randomPoint : randomPoints) {
        float4 direction = {randomPoint[0], randomPoint[1], randomPoint[2], 0};
        direction.normalize();
        Ray ray = Ray(kZeroFloat4, direction);
        if(model->hit(ray, 0.001, 100000.0f, record)){
            record.position = ray.origin + record.t * ray.direction;
            float4 dir = record.position.normalized();
            float4 flatDir = float4(dir.x(), dir.y(), 0, 0);
            if((acosf(dir.dot(flatDir)) / kPi * 180.0f) <= 70)landmarks.emplace_back(record.position);
        }
	}
    writeLandmarks(landmarks, kLandmarksPath, 0);
}

void iodata::readImages(const cv::String& pattern, vector<cv::Mat> &images) {
    cv::Mat rgbImage = cv::imread(pattern);
    cv::Mat grayImage;
    cv::cvtColor(rgbImage, grayImage, COLOR_RGB2GRAY);
    images.emplace_back(grayImage);
}


void iodata::write_ply(Model *model, const string &file_name)
{
    cout << "Writing to " << file_name << endl;
    std::ofstream out;
	Mesh *mesh = model->meshes[0];
    out.open(file_name.c_str());
    out << "ply" << endl;
    out << "format ascii 1.0" << endl;
    out << "element vertex " << model->meshes[0]->numVertices << endl;
    out << "property float x" << endl;
    out << "property float y" << endl;
    out << "property float z" << endl;
    out << "property uchar red" << endl;
    out << "property uchar green" << endl;
    out << "property uchar blue" << endl;
    out << "property uchar alpha" << endl;
    //    out << "property ushort gray" << endl;
    out << "element face " << model->meshes[0]->numTriangles << endl;
    out << "property list uchar int vertex_indices" << endl;
    out << "end_header" << endl;

    for (auto vertex: model->meshes[0]->vertices)
    {
        out << vertex.position.x() << " " << vertex.position.y() << " " << vertex.position.z() << " " << 255 << " " << 255 << " "
            << 255 << " " << 255 << endl;
    }

    for (const auto &triangle: model->meshes[0]->triangles)
    {
        out << "3 " << triangle.vertex_0 << " " << triangle.vertex_1 << " " << triangle.vertex_2 << endl;
    }

    out.close();
}

void iodata::write_depth_image(const float *z_buffer, Camera *c)
{
    cout << "rendering..." << endl;
    string file = "Data/Output/depth_image.png";

    uint16_t z_lut[65536];
    // update the lut
    for (int i = 0; i < 65536; i++)
    {
        z_lut[i] = (uint16_t)(pow(float(i) / 65535.0, GAMMA_VALUE) * 65535.0); // 16 bit gray
    }

    cv::Mat image = cv::Mat::zeros(c->y, c->x, CV_16U);
    auto *p = (uint16_t *)image.data;
    for (int i = 0; i < c->y * c->x; i++)
    {
        *p = 65535 - z_lut[(uint16_t)(z_buffer[i] * 65535.0f)];
        p++;
    }
    cv::imwrite(file, image);
}

void iodata::write_depth_image(SunLight *light)
{
    cout << "rendering..." << endl;
    string file = "Data/Output/depth_image.png";

    uint16_t z_lut[65536]; // look up table
    // update the lut
    for (int i = 0; i < 65536; i++)
    {
        z_lut[i] = (uint16_t)(pow(float(i) / 65535.0, Z_GAMMA_VALUE) * 65535.0); // 16 bit gray
    }

    cv::Mat image = cv::Mat::zeros(light->y, light->x, CV_16U);
    auto *p = (uint16_t *)image.data;
    for (int i = 0; i < light->y * light->x; i++)
    {
        *p = 65535 - z_lut[(uint16_t)(light->shadow_map->map[i] * 65535.0f)];
        p++;
    }
    cv::imwrite(file, image);
}

//void iodata::write_result_image(const FrameBuffer &frame)
//{
//    cout << "output image rendering..." << endl;
//    string file = "data/image.png";
//
//    uint16_t lut[65536]; // look up table
//    // update the lut
//    for (int i = 0; i < 65536; i++)
//    {
//        lut[i] = (uint16_t)(pow(float(i) / 65535.0, GAMMA_VALUE) * 65535.0); // 16 bit gray
//    }
//
//    cv::Mat image = cv::Mat::zeros(frame.y, frame.x, CV_16U);
//    auto *p = (uint16_t *)image.data;
//    for (int i = 0; i < frame.x * frame.y; i++)
//    {
//        *p = lut[(uint16_t)(frame.buffer[i] * 65535.0f)];
//        p++;
//    }
//    cv::imwrite(file, image);
//}

void iodata::write_result_image(const FrameBuffer &frame)
{
    cout << "output image rendering..." << endl;
    string file = "Data/Output/image.png";

    unsigned char lut[256]; // look up table
    // update the gamma correction lut
    for (int i = 0; i < 256; i++)
    {
        lut[i] = (uint8_t)(pow(float(i) / 255.0, GAMMA_VALUE) * 255.0); // 8 bit
    }

    cv::Mat image = cv::Mat::zeros(frame.y, frame.x, CV_8UC4);
    auto *p = (uint8_t *)image.data;
    int size = 4 * frame.x * frame.y;
    for (int i = 0; i < size; i += 4)
    {
        //        *p = lut[frame.buffer[i + 2]]; // B
        //        p++;
        //        *p = lut[frame.buffer[i + 1]]; // G
        //        p++;
        //        *p = lut[frame.buffer[i]]; // R
        //        p++;
        //        *p = frame.buffer[i + 3]; // Alpha channel
        //        p++;
        *p = frame.buffer[i + 2]; // B
        p++;
        *p = frame.buffer[i + 1]; // G
        p++;
        *p = frame.buffer[i]; // R
        p++;
        *p = frame.buffer[i + 3]; // Alpha channel
        p++;
    }
    cv::imwrite(file, image);
}

void iodata::write_result_image(const string &file_name, const FrameBuffer &frame)
{
    cout << "writing " << file_name << endl;
    //    unsigned char lut[256]; // look up table
    // update the gamma correction lut
    //    for (int i = 0; i < 256; i++)
    //    {
    //        lut[i] = (uint8_t)(pow(float(i) / 255.0, GAMMA_VALUE) * 255.0); // 8 bit
    //    }

    cv::Mat image = cv::Mat::zeros(frame.y, frame.x, CV_8UC4);
    auto *p = (uint8_t *)image.data;
    int size = 4 * frame.x * frame.y;
    for (int i = 0; i < size; i += 4)
    {
        //        *p = lut[frame.buffer[i + 2]]; // B
        //        p++;
        //        *p = lut[frame.buffer[i + 1]]; // G
        //        p++;
        //        *p = lut[frame.buffer[i]]; // R
        //        p++;
        //        *p = frame.buffer[i + 3]; // Alpha channel
        //        p++;
        *p = frame.buffer[i + 2]; // B
        p++;
        *p = frame.buffer[i + 1]; // G
        p++;
        *p = frame.buffer[i]; // R
        p++;
        *p = frame.buffer[i + 3]; // Alpha channel
        p++;
    }
	/*int row = frame.y / 2;
	int col = frame.x / 2;
	cv::pyrDown(image, image, cv::Size(row, col));
	row = row / 2;
	col = col / 2;
	cv::pyrDown(image, image, cv::Size(row ,col));*/
    cv::imwrite(file_name, image);
}

void iodata::readLandmarks(const string &path, vector<float4> &landmarks) {
	/*
	path:保存landmarks.txt的文件路径
	landmarks:保存最终读取的结果
	*/
	ifstream infile;
	infile.open(path);
	if (!infile) {
		cout << "open file error!" << endl;
	}

	float x, y, z, w;
	while (!infile.eof()) {
		//读取的点坐标是cm,需要换算成m，且换算尺度
		infile >> x >> y >> z >> w;
		infile.get(); // 读取最后的回车符
		if (infile.peek() == '\n')
			break;
		float4 landmark(x, y, z, w);
		landmarks.push_back(landmark);
	}
	infile.close();
}


void iodata::readCorrespnse(const string&path, vector<vector<int> > &res) {
	ifstream infile;
	infile.open(path.data());   //将文件流对象与文件连接起来 
	assert(infile.is_open());   //若失败,则输出错误消息,并终止程序运行 
	vector<int> tmp;
	string s;
	while (getline(infile, s)) {
		istringstream is(s); //将读出的一行转成数据流进行操作
		int d;
		is >> d;
		while (!is.eof()) {
			is >> d;
			tmp.emplace_back(d);
		}
		res.emplace_back(tmp);
		tmp.clear();
		s.clear();
	}
	infile.close();
}

void iodata::writeImage(const float *orthoDEM, const string& path) {
    cv::Mat temp = Mat(kDEMSize, kDEMSize, CV_8UC1, Scalar::all(0));
    for (int i = 0; i < kDEMSize; i++) {
        for (int j = 0; j < kDEMSize; j++) {
            temp.at<uchar>(i, j) = orthoDEM[i * kDEMSize + j]; // 实际上是做了向下取整的操作
        }
    }
    try {
        imwrite(path, temp);
    }
    catch (cv::Exception& ex) {
        cerr << "Failed to write image" << endl;
        throw exception();
    }
    cout << "write ortho image: " << path << endl;
}

void iodata::writeUcharImage(const uchar *orthoDEM, const string& path) {
    cv::Mat temp = Mat(kDEMSize, kDEMSize, CV_8UC1, Scalar::all(0));
    for (int i = 0; i < kDEMSize; i++) {
        for (int j = 0; j < kDEMSize; j++) {
            temp.at<uchar>(i, j) = orthoDEM[i * kDEMSize + j]; // 实际上是做了向下取整的操作
        }
    }
    try {
        imwrite(path, temp);
    }
    catch (cv::Exception& ex) {
        cerr << "Failed to write image" << endl;
        throw exception();
    }
    cout << "write ortho image: " << path << endl;
}

void iodata::readCamAndSunMes(AnimationConfig *animationConfig, const vector<float4>& landmarks, const string& filePath, State* s) {
    if(animationConfig->frames_config.capacity() > 0)animationConfig->frames_config.clear();

    ifstream in;
//    string filePath = getPathFromEnum(animationConfig->type);
    in.open((filePath + "cameraMessage.txt").c_str());
    if(!in.is_open()){
        cerr << "Fail to open " + filePath + "cameraMessage.txt" << endl;
        throw exception();
    }
    float cameraX, cameraY, cameraZ, focalX, focalY, focalZ;
    int index = 0;
    float4 cameraCenter, cameraFocal;
    float4x4 P;
    while(!in.eof()){
        if(animationConfig->type == ALBEDO){
            in >> cameraX >> cameraY >> cameraZ;
            cameraFocal = landmarks[index++];
        }
        else{
            in >> cameraX >> cameraY >> cameraZ >> focalX >> focalY >> focalZ;
            cameraFocal << focalX, focalY, focalZ, 1.0f;
        }
        cameraCenter << cameraX, cameraY, cameraZ, 1.0f;
        animationConfig->frames_config.emplace_back(cameraCenter, cameraFocal, float4(0, 0, 1.0f, 0), float4 (0, 0, 0, 1.0f), float4 (0, 0, 0, 1.0f), float4 (0, 0, 1.0f, 0), P);
    }
    cout << animationConfig->frames_config.size() << endl;
    in.close();

    in.open((filePath + "sunMessage.txt").c_str());
    if(!in.is_open()){
        cerr << "Fail to open " + filePath + "sunMessage.txt" << endl;
        throw exception();
    }
    float sunX, sunY, sunZ;
    index = 0;
    while(!in.eof()){
        in >> sunX >> sunY >> sunZ;
        float4 sunCenter(sunX, sunY, sunZ, 1.0f);
        animationConfig->frames_config[index++].light_pos << sunCenter;
    }
    in.close();

    cout << "Read " << animationConfig->type << " camera and sun message done" << endl;

    for(auto & i : animationConfig->frames_config){
        s->camera->set_look_at(i.camera_pos, i.camera_focal, i.camera_up);
        i.P = s->camera->P;
    }
}

void iodata::addNoiseToCamPos(vector<AnimationConfig*>& configVec) {
    for(auto & i : configVec){
        for(auto & j : i->frames_config){
            j.camera_pos.x() *= (1 + gaussrand(0, 0.01));
            j.camera_pos.y() *= (1 + gaussrand(0, 0.01));
            j.camera_pos.z() *= (1 + gaussrand(0, 0.01));
        }
        ofstream out;
        out.open(kNoiseCamAndSunMes + getTextFromEnum(i->type) + "/cameraMessage.txt");
        if(!out.is_open()){
            cerr << "fail to open" + getTextFromEnum(i->type) +  " noiseCamMes" << endl;
            throw exception();
        }
        auto tmp = configVec[i->type]->frames_config.size() - 1;
        for(auto & j : i->frames_config){
            out << j.camera_pos.x() << " " << j.camera_pos.y() << " " << j.camera_pos.z();
            if(i->type != ALBEDO){
                out << " " << j.camera_focal.x() << " " << j.camera_focal.y() << " " << j.camera_focal.z();
            }
            if(tmp > 0){
                out << "\n";
                --tmp;
            }
        }
        out.close();
        out.open(kNoiseCamAndSunMes + getTextFromEnum(i->type) + "/sunMessage.txt");
        if(!out.is_open()){
            cerr << "fail to open" + getTextFromEnum(i->type) +  " noiseSunMes" << endl;
            throw exception();
        }
        tmp = configVec[i->type]->frames_config.size() - 1;
        for(auto & j : i->frames_config){
            out << j.light_pos.x() << " " << j.light_pos.y() << " " << j.light_pos.z();
            if(tmp > 0){
                out << "\n";
                --tmp;
            }
        }
        out.close();
    }
}

void iodata::write3DPoints(const vector<float4>& points, const string& path) {
    std::ofstream out;
    out.open(path.c_str());
    out << "ply" << endl;
    out << "format ascii 1.0" << endl;
    out << "element vertex " << points.size() << endl;
    out << "property float x" << endl;
    out << "property float y" << endl;
    out << "property float z" << endl;
    out << "property uchar red" << endl;
    out << "property uchar green" << endl;
    out << "property uchar blue" << endl;
    out << "property uchar alpha" << endl;
    out << "end_header" << endl;

    for (auto & vertex: points)
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << " " << 255 << " " << 255 << " "
            << 255 << " " << 255 << endl;
    }
    out.close();
    cout << "Write 3D points done" << endl;
}

