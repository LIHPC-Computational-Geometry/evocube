

#include <omp.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/per_face_normals.h>

#include <chrono>

#include "distortion.h"


int main(int argc, char *argv[]){
    using std::chrono::steady_clock;
    using std::chrono::duration_cast;
    using std::chrono::microseconds; 

    int NUM_THREADS = omp_get_num_procs();
    omp_set_num_threads(NUM_THREADS);
    omp_set_num_threads(10);

    std::cout << "# of threads: " << omp_get_num_threads() << std::endl;

    Eigen::MatrixXd V, V_def;
    Eigen::MatrixXi F;

    //igl::readOBJ("../data/cross/cross_flat.obj", V, F);
    //igl::readOBJ("../data/cross/cross_flat_bad.obj", V_def, F);


    //igl::readOBJ("../data/sphere.obj", V, F);
    //igl::readOBJ("../data/sphere_def.obj", V_def, F);
    
    
    //igl::readOBJ("../data/abc/boundary.obj", V, F);
    //igl::readOBJ("../data/abc/boundary_ce.obj", V_def, F);
    
    igl::readOBJ("../data/S1/boundary.obj", V, F);
    igl::readOBJ("../data/S1/boundary_ce.obj", V_def, F);

    /*V_def.resize(4, 3);
    V.resize(4, 3);
    F.resize(2, 3);

    double dx = 1.0; //0.8 * std::sqrt(2.0)/2.0;
    double dy = 0.0;
    V_def <<  0,   0, 0,
           1.0,   dy, 0,
             dx,   1.0, 0,
           1.0 + dx, 1.0 + dy, 0;

    V <<  0,   0, 0,
           1.0,   0, 0,
             0, 1.0, 0,
           1.0, 1.0, 0;

    F << 0, 1, 2,
         1, 3, 2;*/

    Eigen::MatrixXd N;
    igl::per_face_normals(V, F, N);
    Eigen::MatrixXd N_def;
    igl::per_face_normals(V_def, F, N_def);

    auto time_disto_start = steady_clock::now();
    Eigen::VectorXd disto;
    computeDisto(V, V_def, F, N, N_def, disto);
    auto time_disto_end = steady_clock::now();
    std::cout << "Total disto time (μs): "<< duration_cast<microseconds>(time_disto_end - time_disto_start).count() << std::endl;

    std::cout << "disto.mean(): " << disto.mean() << std::endl;    
    double final_disto = integrateDistortion(V, F, disto);
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

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();  
    };

    //helper function for menu
    auto make_checkbox = [&](const char *label, unsigned int &option) {
        return ImGui::Checkbox(
            label,
            [&]() { return viewer.core().is_set(option); },
            [&](bool value) { return viewer.core().set(option, value); });
    };

    menu.callback_draw_custom_window = [&]() {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(350, -1), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("IGL")) {
            
            
            ImGui::End();
        }
    
    };

    updateViz();

    viewer.data().show_lines = true;
    viewer.core().lighting_factor = 0.5;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}