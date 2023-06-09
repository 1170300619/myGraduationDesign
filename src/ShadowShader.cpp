#include "../include/ShadowShader.hpp"
#include "../include/Light.hpp"
#include <vector>

#define MIN_LUM 0.003921569

ShadowShader *ShadowShader::shader = nullptr;

ShadowShader::~ShadowShader() = default;

ShadowShader *ShadowShader::instance()
{
	if (shader == nullptr)
	{
		shader = new ShadowShader();
	}
	return shader;
}

void ShadowShader::destroy()
{
	uniform = nullptr;
	delete shader;
	shader = nullptr;
}

void ShadowShader::vertex_shader(VertexP &vertex)
{
	vertex.clip = uniform->camera->VP * vertex.position;
	vertex.z_rec = 1.0f / vertex.clip.w();
}

void ShadowShader::fragment_shader(Fragment &frag)
{
	float s = shadow((*uniform->light_source)[0], frag);
	frag.color << s, s, s, 0;
}

float ShadowShader::shadow(Light *light, const Fragment &frag)
{
	if (!light->shadow_map)
	{
		return 1.0f;
	}
	float4 shadow_map_index;
	float bias;
	
	auto *sun = (SunLight *)light;
	shadow_map_index = sun->MO * frag.world;
	bias = shadow_map_index.z() * 0.005f;
	int index_x = (int)shadow_map_index.x();
	int index_y = (int)shadow_map_index.y();
	if (index_x < 0 || index_x >= sun->x || index_y < 0 || index_y >= sun->y)
	{
		return 1.0f;
	}

	// shadow blur by 3x3 PCF
	//        int min_x = index_x - 1 < 0 ? 0 : index_x - 1;
	//        int max_x = index_x + 1 >= sun->x ? sun->x - 1 : index_x + 1;
	//        int min_y = index_y - 1 < 0 ? 0 : index_y - 1;
	//        int max_y = index_y + 1 >= sun->y ? sun->y - 1 : index_y + 1;

	// shadow blur by 5x5 PCF
	//        int min_x = index_x - 2 < 0 ? 0 : index_x - 2;
	//        int max_x = index_x + 2 >= sun->x ? sun->x - 1 : index_x + 2;
	//        int min_y = index_y - 2 < 0 ? 0 : index_y - 2;
	//        int max_y = index_y + 2 >= sun->y ? sun->y - 1 : index_y + 2;
	//
	//        float shadow = 0;
	//        for (int i = min_y; i <= max_y; i++)
	//        {
	//            for (int j = min_x; j <= max_x; j++)
	//            {
	//                float depth = shadow_map_index.z() - bias; // bias
	//                //            float deviation = abs(depth - light->shadow_map->map[index_y * sun->x + index_x]);
	//                float deviation = depth - sun->shadow_map->map[i * sun->x + j];
	//                if (deviation > 0)
	//                {
	//                    shadow += 0.04f;
	//                } else
	//                {
	//                    shadow += 1.0f;
	//                }
	//            }
	//        }
	//        return shadow / float((max_x - min_x + 1) * (max_y - min_y + 1));
	return shadow_map_index.z() - bias - sun->shadow_map->map[index_y * sun->x + index_x] > 0 ? 0 : 1.0f;
}
