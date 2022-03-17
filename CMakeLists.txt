
cmake_minimum_required(VERSION 3.1)
project(evocube)
set(CMAKE_CXX_STANDARD 17)
file(GLOB SRCFILES src/*.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/libHexEx/cmake)
set(OPENVOLUMEMESH_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/OpenVolumeMesh/src)

add_definitions(-DIGL_VIEWER_VIEWER_QUIET)
set(LIBIGL_USE_STATIC_LIBRARY ON)
add_subdirectory(lib/libigl)
add_subdirectory(lib/gco)
add_subdirectory(lib/OpenNL_psm)
add_subdirectory(lib/OpenVolumeMesh)
add_subdirectory(lib/libHexEx)

find_package(OpenMP REQUIRED)

include_directories(include)
include_directories(lib/json/include)
include_directories(lib/gco/include)
include_directories(lib/OpenNL_psm)
include_directories(lib/OpenVolumeMesh)
include_directories(lib/libHexEx/src)

SET(IGL_VIEWER_LIBS
    igl::core
    igl::opengl_glfw
    igl::opengl_glfw_imgui
)

SET(EVOCUBE_LIBS
    gco
    OpenMP::OpenMP_CXX
    OpenNL_psm
    OpenVolumeMesh 
    HexEx
)

add_executable(viz_tet app/viz_tet.cpp src/flagging_utils.cpp src/tet_boundary.cpp src/logging.cpp)
target_link_libraries(viz_tet ${IGL_VIEWER_LIBS} igl::tetgen igl::png)

add_executable(jacob_test app/jacob_test.cpp ${SRCFILES})
target_link_libraries(jacob_test ${IGL_VIEWER_LIBS} ${EVOCUBE_LIBS})

add_executable(flagging_viewer app/flagging_viewer.cpp ${SRCFILES})
target_link_libraries(flagging_viewer ${IGL_VIEWER_LIBS} OpenMP::OpenMP_CXX)

add_executable(edit_labeling app/edit_labeling.cpp ${SRCFILES})
target_link_libraries(edit_labeling ${EVOCUBE_LIBS} ${IGL_VIEWER_LIBS})

add_executable(test_greedy app/test_greedy.cpp ${SRCFILES})
target_link_libraries(test_greedy ${IGL_VIEWER_LIBS} ${EVOCUBE_LIBS} igl::tetgen)

add_executable(init_from_folder app/init_from_folder.cpp ${SRCFILES})
target_link_libraries(init_from_folder igl::core igl::tetgen ${EVOCUBE_LIBS})

add_executable(polycube_withHexEx app/polycube_withHexEx.cpp ${SRCFILES})
target_link_libraries(polycube_withHexEx igl::core igl::tetgen ${EVOCUBE_LIBS})

add_executable(test_scaled_jacobian app/test_scaled_jacobian.cpp ${SRCFILES})
target_link_libraries(test_scaled_jacobian igl::core igl::tetgen ${EVOCUBE_LIBS})