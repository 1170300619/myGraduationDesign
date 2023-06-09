

#ifndef RENDERER_RENDERER_HPP
#define RENDERER_RENDERER_HPP

#include "Type.hpp"
#include "State.hpp"
#include "AnimationConfig.hpp"

class Renderer
{
public:
    State *s;
    bool stop; // if loop should stop rendering
    vector<AnimationConfig*> configVec;
    vector<float4> landmarks;

    Renderer():stop(false), s(nullptr){}

    explicit Renderer(bool _stop):stop(_stop), s(nullptr){}

    ~Renderer()
    {
        delete s;
    }

    /**
     *
     * @param model_file
     * @param camera_file
     */
    void read_data(const std::string &model_file, const std::string &config_file);

    /**
     *
     */
    void render();
private:
    /**
     *
     */
    void view_transform();

    /**
     * Homogeneous space clipping.
     * @param v_0
     * @param v_1
     * @param v_2
     * @return
     */
    static bool clipping(const VertexP &v_0, const VertexP &v_1, const VertexP &v_2);

    /**
     * Homogeneous space clipping
     * @param v_0
     * @param v_1
     * @param v_2
     * @return
     */
    static std::vector<VertexP> clip_near(const VertexP &v_0, const VertexP &v_1, const VertexP &v_2);

    /**
     *
     */
    void draw(int shader_mode);

    void calcCoRES(AnimationConfig* animationConfig) ;

    /**
     *
     * @param v_0
     * @param v_1
     * @param v_2
     * @param normal
     */
    inline void draw_triangle(const VertexP &v_0, const VertexP &v_1, const VertexP &v_2, const float4 &normal);

    /**
     * Write fragment color to buffer and test depth.
     * @param frag
     */
    void write_color(Fragment &frag) const;

};

#endif //RENDERER_RENDERER_HPP
