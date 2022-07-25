#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/avg_edge_length.h>

#include "flagging_utils.h"
#include "labeling_individual.h"
#include "evocube.h"
#include "quick_label_ev.h"
#include "evaluator.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::string input_path = argv[1];
    igl::readOBJ(input_path, V, F);
    double l_avg = igl::avg_edge_length(V, F);

    Eigen::VectorXi labeling = Eigen::VectorXi::Constant(F.rows(), 1);
    Eigen::MatrixXd colors;

    std::shared_ptr<Evocube> evo = std::make_shared<Evocube>(Evocube(V, F));
    std::shared_ptr<const QuickLabelEv> qle = std::make_shared<const QuickLabelEv>(QuickLabelEv(V, F));
    const Evaluator evaluator(evo, qle); 
    std::shared_ptr<LabelingIndividual> indiv = std::make_shared<LabelingIndividual>(LabelingIndividual(evo, qle, labeling));

    double labeling_score = -1.0;
    double fast_poly_score, invalid_score;
    int n_fail_invert;
    auto updateLabeling = [&](){
        indiv = std::make_shared<LabelingIndividual>(LabelingIndividual(evo, qle, labeling));
        indiv->updateChartsAndTPs(true);
        labeling_score = evaluator.evaluate(*indiv, fast_poly_score, invalid_score, n_fail_invert);
    };
    

    // --- VISUALIZATION ---

    int insert_label = 0;
    float insert_size = 1.0;

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
            if (ImGui::Button("Quick open labeling")){
                std::string file = input_path.substr(0, input_path.find_last_of("/\\") + 1) + "labeling.txt";
                std::cout << "Attempting to open: " << file << std::endl;
                labeling = openFlagging(file, F.rows());
                updateLabeling();
                updateViz();
            }
            
            if (ImGui::Button("Open labeling")){
                std::string file = igl::file_dialog_open();
                labeling = openFlagging(file, F.rows());
                updateLabeling();
                updateViz();
            }

            if (ImGui::Button("Save labeling")){
                std::string file = igl::file_dialog_save();
                saveFlagging(file, labeling);
            }

            if (ImGui::Button("Measure labeling")){
                std::cout << "Measuring" << std::endl;
            }

            ImGui::Text("Total labeling score: %f", labeling_score);
            ImGui::Text("\tFastPoly score: %f", fast_poly_score);
            ImGui::Text("\tInvalid score: %f", invalid_score);
            ImGui::Text("\t# inverted tris: %i", n_fail_invert);
            
            ImGui::SliderInt("Insert label", &insert_label, 0, 5);   
            ImGui::SliderFloat("Insert size", &insert_size, 0.0*l_avg, 10.0*l_avg);   
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
                //labeling(fid) = insert_label;
                Eigen::VectorXi old_labeling = labeling;
                growFromTri(old_labeling, labeling, evo->TT_, evo->dists_, insert_label, fid, insert_size);
        
                updateLabeling();
                updateViz();
                return true;
            }
        }
        return false;
    };

    updateViz();

    viewer.core().lighting_factor = 0.0;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    //viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.core().background_color = Eigen::Vector4f(255.0/255.0, 255.0/255.0, 255.0/255.0, 1.0);
    viewer.launch();
}