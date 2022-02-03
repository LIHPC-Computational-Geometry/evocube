#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>

/*
#include <utility>
#include <imgui.h>
#include <chrono>
#include <random>
*/

#include "flagging_utils.h"

int main(int argc, char *argv[]){
    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    int orig_id = viewer.data_list[0].id;
    int hud_id = viewer.data_list[1].id;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    /*std::string input_name = "bunny"; //"cow";
    std::string input_off = "../data/"+input_name+".off";
    igl::readOFF(input_off, V, F);*/
    
    std::string input_obj = "../data/S1/boundary.obj";
    std::string input_labeling = "../data/S1/labeling.txt"; 
    
    igl::readOBJ(input_obj, V, F);

    viewer.data(orig_id).set_mesh(V, F);


    Eigen::VectorXi flagging = openFlagging(input_labeling, F.rows());
    Eigen::MatrixXd flagging_colors = colorsFromFlagging(flagging);
    viewer.data(orig_id).set_colors(flagging_colors);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    menu.callback_draw_custom_window = [&]() {
        bool show = true;
        ImGui::SetNextWindowPos(ImVec2(0.f * menu.menu_scaling(), 0),
                                    ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(350, 400));
        if (ImGui::Begin("Test flagging")){

            ImGui::SliderFloat("Lighting", &viewer.core().lighting_factor, 0.0f, 5.0f, "%.3f");

            ImGui::End();
        }
    };

    viewer.data(orig_id).line_width = 8;
    viewer.data(hud_id).line_width = 1;
    viewer.data(orig_id).point_size = 4;


    viewer.data(orig_id).show_lines = false;
    viewer.core().background_color.setZero();
    viewer.core().lighting_factor=0.5;
    viewer.launch();
}