cmake_minimum_required(VERSION 3.23)
project(SPC)
set(CMAKE_CXX_STANDARD 14)

include_directories("C:/ceres/newInstall/eigen-3.3.9/include/eigen3")

set(OpenCV_DIR C:\\opencv\\opencv\\mingw-build\\install)
find_package(OpenCV REQUIRED)
include_directories("C:/opencv/opencv/mingw-build/install/include")

#find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
include_directories("C:/SPC_Files/ceres/install/ceres/include")
include_directories("C:/SPC_Files/ceres/install/glog/include")
include_directories("C:/SPC_Files/ceres/install/gflags/Include")
link_directories("C:/SPC_Files/ceres/install/ceres/bin")
link_directories("C:/SPC_Files/ceres/install/glog/lib")
link_directories("C:/SPC_Files/ceres/install/gflags/Lib")

add_executable(SPC src/main.cpp src/AnimationConfig.cpp src/BundleAdjustment.cpp src/SphericalSurfICQ.cpp src/Camera.cpp src/FlatShader.cpp src/ICQ.cpp src/iodata.cpp src/Light.cpp src/Lmap.cpp src/Material.cpp src/Mipmap.cpp src/Model.cpp src/optimation.cpp src/PhongMaterial.cpp src/PhongShader.cpp src/Possion_recon.cpp src/Renderer.cpp  src/ShadowShader.cpp src/State.cpp src/Texture.cpp include/AnimationConfig.hpp include/BundleAdjustment.hpp include/SphericalSurfICQ.hpp include/Camera.hpp include/FlatShader.hpp include/ICQ.hpp  include/iodata.hpp include/Light.hpp include/Lmap.hpp include/Material.hpp include/Mipmap.hpp include/Model.hpp include/optimation.hpp include/PhongMaterial.hpp include/PhongShader.hpp include/Possion_recon.hpp include/Renderer.hpp  include/ShadowShader.hpp include/State.hpp include/Texture.hpp include/Fragment.hpp include/RandomPoint.hpp include/Shader.hpp  include/Triangle.hpp include/Type.hpp include/U_Ptr.hpp include/Util.hpp include/Vertex.hpp include/BoundingBox.hpp include/Ray.hpp src/BoundingBox.cpp include/HitRecord.hpp include/Primitive.hpp src/Triangle.cpp src/BVH.cpp include/BVH.hpp
        )

#add_executable(ceresCurve src/optimation.cpp)
#target_link_libraries(ceresCurve ${CERES_LIBRARIES})
target_link_libraries(SPC ${catkin_LIBRARIES} ${OpenCV_LIBS} ${CERES_LIBRARIES} ${GLOG_LIBRARIES} ${Eigen_LIBS})
target_link_libraries(SPC "C:/SPC_Files/ceres/install/ceres/lib/libceres.dll.a")
#target_link_libraries(SPC "C:/SPC_Files/ceres/install/ceres/lib/libgtest.dll.a")
target_link_libraries(SPC "C:/SPC_Files/ceres/install/glog/lib/libglog.dll.a")
target_link_libraries(SPC "C:/SPC_Files/ceres/install/gflags/Lib/libgflags.dll.a")
target_link_libraries(SPC "C:/SPC_Files/ceres/install/gflags/Lib/libgflags_nothreads.dll.a")

add_compile_definitions(CERES_USE_CXX_THREADS)
