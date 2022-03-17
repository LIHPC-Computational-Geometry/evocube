#include <omp.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>
#include <chrono>

#include "distortion.h"

int main(int argc, char *argv[]){
    using std::chrono::steady_clock;
    using std::chrono::duration_cast;
    using std::chrono::microseconds; 

    Eigen::MatrixXd V, V_def;
    Eigen::MatrixXi F;

    igl::readOBJ("../data/DATASET2/medium_mambo/M1/boundary.obj", V, F);
    igl::readOBJ("../data/DATASET2/medium_mambo/M1/fast_polycube_surf.obj", V_def, F);

    Eigen::MatrixXd N;
    igl::per_face_normals(V, F, N);
    Eigen::MatrixXd N_def;
    igl::per_face_normals(V_def, F, N_def);
    Eigen::VectorXd A;
    igl::doublearea(V, F, A);

    auto time_disto_start = steady_clock::now();
    Eigen::VectorXd disto;
    computeDisto(V, V_def, F, N, N_def, disto);
    auto time_disto_end = steady_clock::now();
    std::cout << "Total disto time (μs): "<< duration_cast<microseconds>(time_disto_end - time_disto_start).count() << std::endl;

    std::cout << "disto.mean(): " << disto.mean() << std::endl;    
    double final_disto = integrateDistortion(A, disto);
    std::cout << "final_disto: " << final_disto << std::endl;    

    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_def, F);

    Eigen::MatrixXd colors = Eigen::MatrixXd::Constant(disto.rows(), 3, 1.0);
    colors.col(1) = 1.0 - disto.array();
    colors.col(2) = 1.0 - disto.array();
    viewer.data().set_colors(colors);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);
    menu.callback_draw_custom_window = [&](){};

    viewer.data().show_lines = true;
    viewer.core().lighting_factor = 0.0;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}