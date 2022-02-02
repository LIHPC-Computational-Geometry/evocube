
cmake_minimum_required(VERSION 3.1)
project(evocube)

add_definitions(-DIGL_VIEWER_VIEWER_QUIET)
set(LIBIGL_USE_STATIC_LIBRARY ON)
add_subdirectory(lib/libigl)

find_package(OpenMP REQUIRED)

include_directories(include)

add_executable(igl_core app/igl_core.cpp)
target_link_libraries(igl_core igl::opengl_glfw igl::core igl::opengl_glfw_imgui)

add_executable(jacob_test app/jacob_test.cpp)
target_link_libraries(jacob_test igl::opengl_glfw igl::core igl::opengl_glfw_imgui OpenMP::OpenMP_CXX)

add_executable(openmp_test app/openmp_test.cpp)
target_link_libraries(openmp_test OpenMP::OpenMP_CXX)

add_executable(tetrahedralize_app app/tetrahedralize_app.cpp)
target_include_directories(tetrahedralize_app PUBLIC lib/libigl/external/tetgen)
target_link_libraries(tetrahedralize_app igl::core igl::tetgen)