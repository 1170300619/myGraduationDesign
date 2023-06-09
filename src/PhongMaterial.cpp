#include "../include/PhongMaterial.hpp"
#include "../include/PhongShader.hpp"

PhongMaterial::PhongMaterial(){
    set_shader(PhongShader::instance());
//    shadow_shader = ShadowShader::instance();
}

void PhongMaterial::destroy()
{
    base_texture.reset();
}
