#include <iostream>
#include "../include/Renderer.hpp"

using namespace std;
using namespace Render;

int main(int argc, char *argv[])
{
    string config_file = argv[1]; // config file name
    auto *render = new Renderer();
    render->read_data(kInitModelPath, config_file);
    render->render();
    delete render;
    return 0;
}
