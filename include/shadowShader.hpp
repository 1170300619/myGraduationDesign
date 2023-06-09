#pragma once

#ifndef RENDERER_FLATSHADER_HPP
#define RENDERER_FLATSHADER_HPP

#include "Shader.hpp"

class ShadowShader :public Shader
{
	static ShadowShader *shader;

	ShadowShader() = default;
public:
	~ShadowShader() override;

	/**
	* Get a instance and reset.
	*/
	static ShadowShader *instance();

	/**
	* Shader. Automatically select the type.
	* @param vertex
	*/
	void vertex_shader(VertexP &vertex) override;

	/**
	* Shader. Automatically select the type.
	* @param frag
	*/
	void fragment_shader(Fragment &frag) override;

	void destroy() override;

	static float shadow(Light *light, const Fragment &frag);
};

#endif //RENDERER_FLATSHADER_HPP
