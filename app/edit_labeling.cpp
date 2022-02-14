#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/unproject_onto_mesh.h>

#include "flagging_utils.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    igl::readOBJ("../data/S1/boundary.obj", V, F);



    Eigen::VectorXi labeling = Eigen::VectorXi::Constant(F.rows(), 1);
    Eigen::MatrixXd colors;
    

    // --- VISUALIZATION ---

    int insert_label = 0;

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();
        colors = colorsFromFlagging(labeling);
        viewer.data().set_colors(colors);
        
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
            
            if (ImGui::Button("Open labeling")){
                std::string file = igl::file_dialog_open();
                labeling = openFlagging(file, F.rows());
                updateViz();
            }

            if (ImGui::Button("Save labeling")){
                std::string file = igl::file_dialog_save();
                saveFlagging(file, labeling);
            }

            ImGui::SliderInt("Insert label", &insert_label, 0, 5);   
            ImGui::End();
        }
    };

    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
        if (key == ' '){
            int fid;
            Eigen::Vector3f bc;
            // Cast a ray in the view direction starting from the mouse position
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core().view,
                viewer.core().proj, viewer.core().viewport, V, F, fid, bc)){
                labeling(fid) = insert_label;
                
                updateViz();
                return true;
            }
        }
        return false;
    };

    updateViz();

    viewer.core().lighting_factor = 0.0;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}