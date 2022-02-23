#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <ctime>
#include <random>
#include <queue>

#include "graphcut_labeling.h"
#include "flagging_utils.h"
#include "logging.h"
#include "chart.h"
#include "quick_label_ev.h"
#include "evocube.h"
#include "evaluator.h"
#include "labeling_individual.h"

int main(int argc, char *argv[]){

    std::string input_tris = "../data/bunny/boundary.obj";
    if (argc > 1) input_tris = argv[1];

    //std::srand(std::time(nullptr));

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //igl::readOBJ("../data/Basic/B23/boundary.obj", V, F);
    igl::readOBJ(input_tris, V, F);

    const Evocube evo(V, F);
    const QuickLabelEv qle(V, F);
    const Evaluator evaluator(std::make_shared<Evocube>(evo),
                              std::make_shared<QuickLabelEv>(qle)); 

    
    int compact_coeff = 1;
    int fidelity_coeff = 3;
    Eigen::VectorXi labeling_init = graphcutFlagging(V, F, evo.N_, evo.TT_, compact_coeff, fidelity_coeff);
    LabelingIndividual ancestor(std::make_shared<Evocube>(evo), std::make_shared<QuickLabelEv>(qle), labeling_init);

    std::cout << "Fix init... " << std::endl;
    ancestor.updateChartsAndTPs();
    ancestor.repairHighValenceCorner();
    ancestor.updateChartsAndTPs();
    ancestor.repairOppositeLabels();
    ancestor.updateChartsAndTPs();

    std::cout << "Unspike init... " << std::endl;

    ancestor.repairUnspikeLabeling();
    ancestor.updateChartsAndTPs(true);

    int mut_count = 1;
    
    double current_score = evaluator.evaluate(ancestor);
    Eigen::MatrixXd def_V = qle.computeDeformedV(ancestor.getLabeling());

    std::shared_ptr<LabelingIndividual> best_indiv = std::make_shared<LabelingIndividual>(ancestor); 

    auto pickMutation = [](double progress_percentage){
        double r = static_cast<double>(std::rand() % 1000)/1000.0;
        double prob0 = 0.7 - 0.4 * progress_percentage;
        double prob1 = 0.0 + 0.4 * progress_percentage;
        double prob2 = 0.3 - 0.0 * progress_percentage;
        if (r < prob0) return 0;
        if (r < prob0 + prob1) return 1;
        else return 2;
    };

    int max_mut = 100;
    for (int i=0; i<max_mut; i++){
        std::cout << "Attempt: " << i << std::endl;

	    auto time_new_indiv = std::chrono::steady_clock::now();    
        LabelingIndividual new_indiv(*best_indiv);
        new_indiv.updateChartsAndTPs(true); // OPTIM: could be skipped if we kept info from parent
   
	    auto time_before_mutation = std::chrono::steady_clock::now();     

        int mutation_type = std::rand() % 3;
        mutation_type = pickMutation(static_cast<double>(i)/max_mut);
        //mutation_type = 0;
        if (mutation_type == 0) new_indiv.mutationVertexGrow();
        if (mutation_type == 1) new_indiv.mutationRemoveChart();
        if (mutation_type == 2) new_indiv.mutationGreedyPath();

	    auto time_after_mutation = std::chrono::steady_clock::now();

        new_indiv.updateChartsAndTPs();
        new_indiv.repairUnspikeLabeling();
        new_indiv.updateChartsAndTPs(true);

	    auto time_after_charts = std::chrono::steady_clock::now();
        double time_create_indiv = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_before_mutation - time_new_indiv).count()) / 1000;
        double time_mut = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_mutation - time_before_mutation).count()) / 1000;
        double time_chart = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_charts - time_after_mutation).count()) / 1000;
        double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_charts - time_before_mutation).count()) / 1000;
	    std::cout << "Create indiv + mutation + chart time: " << time_create_indiv << " + " << time_mut << " + " << time_chart << " = " << totaltime << "sec." << std::endl;

        double new_score = evaluator.evaluate(new_indiv);
        if (new_score < current_score){
            current_score = new_score; 
            best_indiv = std::make_shared<LabelingIndividual>(new_indiv);
            mut_count ++;
        }
    }

    LabelingIndividual final_indiv = *best_indiv;
    final_indiv.updateChartsAndTPs(true);
    def_V = qle.computeDeformedV(final_indiv.getLabeling());
    

    std::vector<Eigen::MatrixXd> vec_border_begs, vec_border_ends;
    std::vector<Eigen::RowVector3d> vec_border_colors;
    final_indiv.getBordersViz(vec_border_begs, vec_border_ends, vec_border_colors);

    Eigen::MatrixXd threshold_colors = qle.distoAboveThreshold(final_indiv.getLabeling(), 10e10);

    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    int orig_id = viewer.data_list[0].id;
    int hud_id = viewer.data_list[1].id;

    viewer.data(orig_id).set_mesh(V, F);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data(orig_id).clear_edges();
        viewer.data(orig_id).clear_points();

        Eigen::MatrixXd colors = colorsFromFlagging(final_indiv.getLabeling());
        viewer.data(orig_id).set_colors(colors);
        
        for (int i=0; i<vec_border_begs.size(); i++){
            viewer.data(hud_id).add_edges(vec_border_begs[i], vec_border_ends[i], vec_border_colors[i]);
        }

        Eigen::MatrixXd tp_points = final_indiv.getTurningPointsMat();
        viewer.data(hud_id).add_points(tp_points, Eigen::RowVector3d(1.0, 1.0, 0.0));
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
            if (ImGui::Button("Init labeling", ImVec2(-1, 0))){
                viewer.data(orig_id).set_mesh(V, F);
                Eigen::MatrixXd colors = colorsFromFlagging(labeling_init);
                viewer.data(orig_id).set_colors(colors);
            }
            if (ImGui::Button("Labeling", ImVec2(-1, 0))){
                viewer.data(orig_id).set_mesh(V, F);
                Eigen::MatrixXd colors = colorsFromFlagging(final_indiv.getLabeling());
                viewer.data(orig_id).set_colors(colors);
            }
            if (ImGui::Button("Fastbndpolycube", ImVec2(-1, 0))){
                viewer.data(orig_id).set_mesh(def_V, F);
            }
            if (ImGui::Button("Threshold dist", ImVec2(-1, 0))){
                viewer.data(orig_id).set_colors(threshold_colors);
            }

            if (ImGui::Button("Save labeling to B23", ImVec2(-1, 0))){
                std::string folder = "../data/Basic/B23/";
                Eigen::VectorXi save_labeling = final_indiv.getLabeling();
                saveFlagging(folder + "labeling.txt", save_labeling);
                saveFlaggingOnTets(folder + "labeling_on_tets.txt", folder + "tris_to_tets.txt", save_labeling);
            }

            make_checkbox("Show mesh", viewer.data(orig_id).show_lines);
            ImGui::End();
        }
    };

    updateViz();

    viewer.data(hud_id).line_width = 10.0;
    viewer.core().lighting_factor = 0.0;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();

    return 0;
}