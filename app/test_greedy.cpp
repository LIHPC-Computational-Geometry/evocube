#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/file_dialog_save.h>
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
#include "archive.h"

#define READ_DOT_MESH
#ifdef READ_DOT_MESH
#include "mesh_io.h" // requires tetgen
#endif


int main(int argc, char *argv[]){

    std::string input_tris = "../data/bunny/boundary.obj";
    if (argc > 1) input_tris = argv[1];

    std::srand(std::time(nullptr));

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //igl::readOBJ("../data/Basic/B23/boundary.obj", V, F);
    igl::readOBJ(input_tris, V, F);
    //readDotMeshTri(input_tris, V, F);

    std::shared_ptr<Evocube> evo = std::make_shared<Evocube>(Evocube(V, F));
    std::shared_ptr<const QuickLabelEv> qle = std::make_shared<const QuickLabelEv>(QuickLabelEv(V, F));
    const Evaluator evaluator(evo, qle); 

    int compact_coeff = 1;
    int fidelity_coeff = 3;
    Eigen::VectorXi labeling_init = graphcutFlagging(V, F, evo->N_, evo->TT_, compact_coeff, fidelity_coeff);
    std::shared_ptr<LabelingIndividual> ancestor = std::make_shared<LabelingIndividual>(LabelingIndividual(evo, qle, labeling_init));

    std::cout << "Fix init... " << std::endl;
    ancestor->updateChartsAndTPs();
    ancestor->repairHighValenceCorner();
    ancestor->updateChartsAndTPs();
    ancestor->repairOppositeLabels();
    ancestor->updateChartsAndTPs();

    std::cout << "Unspike init... " << std::endl;

    ancestor->repairUnspikeLabeling();
    ancestor->updateChartsAndTPs(true);
    ancestor->updateTimestamps();
    
    int mut_count = 1;
    
    double current_score = evaluator.evaluate(*ancestor);
    Eigen::MatrixXd def_V;// = qle->computeDeformedV(ancestor.getLabeling());

    Archive<std::shared_ptr<LabelingIndividual>> archive(10);
    archive.insert(ancestor, current_score);

    auto pickMutation = [](double progress_percentage){
        double r = static_cast<double>(std::rand() % 1000)/1000.0;
        double prob0 = 0.7 - 0.3 * progress_percentage;
        double prob1 = 0.0 + 0.3 * progress_percentage;
        double prob2 = 0.3 - 0.0 * progress_percentage;
        if (r < prob0) return 0;
        if (r < prob0 + prob1) return 1;
        else return 2;
    };

    int n_generations = 20;
    for (int generation=0; generation<n_generations; generation++){
        int max_mut = 30;
        Archive<std::shared_ptr<LabelingIndividual>> gen_archive(5);

        for (int i=0; i<max_mut; i++){
            evo->timestamp_ ++;
            std::cout << "Generation " << generation << ", Mutation: " << i << std::endl;

            auto time_new_indiv = std::chrono::steady_clock::now();    
            //LabelingIndividual new_indiv(*best_indiv);
            //std::shared_ptr<LabelingIndividual> new_indiv = std::make_shared<LabelingIndividual>(LabelingIndividual(*best_indiv));
            
            int pick = archive.probabilisticIndividual();
            std::shared_ptr<LabelingIndividual> new_indiv = std::make_shared<LabelingIndividual>(*archive.getIndiv(pick));
            new_indiv->updateChartsAndTPs(true); // OPTIM: could be skipped if we kept info from parent
    
            auto time_before_mutation = std::chrono::steady_clock::now();     

            int mutation_type = std::rand() % 3;
            mutation_type = pickMutation(static_cast<double>(generation)/n_generations);
            //mutation_type = 0;
            if (mutation_type == 0) new_indiv->mutationVertexGrow();
            if (mutation_type == 1) new_indiv->mutationRemoveChart();
            if (mutation_type == 2) new_indiv->mutationGreedyPath();

            auto time_after_mutation = std::chrono::steady_clock::now();

            new_indiv->updateChartsAndTPs();
            new_indiv->repairUnspikeLabeling();
            new_indiv->updateChartsAndTPs(true);
            new_indiv->updateTimestamps();

            auto time_after_charts = std::chrono::steady_clock::now();
            
            double new_score = evaluator.evaluate(*new_indiv);

            auto time_after_eval = std::chrono::steady_clock::now();

            gen_archive.insert(new_indiv, new_score);
            /*if (new_score < current_score){
                current_score = new_score; 
                best_indiv = new_indiv;
                mut_count ++;
            }*/

            auto time_after_archive = std::chrono::steady_clock::now();
            /*double time_create_indiv = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_before_mutation - time_new_indiv).count()) / 1000;
            double time_mut = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_mutation - time_before_mutation).count()) / 1000;
            double time_chart = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_charts - time_after_mutation).count()) / 1000;
            double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_charts - time_before_mutation).count()) / 1000;
            */
            double time_create_indiv = measureTime(time_new_indiv, time_before_mutation);
            double time_mut = measureTime(time_before_mutation, time_after_mutation);
            double time_chart = measureTime(time_after_mutation, time_after_charts);
            double time_eval = measureTime(time_after_charts, time_after_eval);
            double time_archive = measureTime(time_after_eval, time_after_archive);
            double totaltime = measureTime(time_new_indiv, time_after_archive);

            std::cout << "Individual mutation, time repartition:" << std::endl;
            std::cout << "\tCreation \t" << time_create_indiv << std::endl;
            std::cout << "\tMutation \t" << time_mut << std::endl;
            std::cout << "\tCharts   \t" << time_chart << std::endl;
            std::cout << "\tEval.    \t" << time_eval << std::endl;
            std::cout << "\tArchive  \t" << time_archive << std::endl;
            std::cout << "\tTotal:   \t" << totaltime << std::endl;
        }

        auto time_after_mutations = std::chrono::steady_clock::now();

        int n_cross = 10;
        for (int i=0; i<n_cross; i++){
            std::cout << "Generation " << generation << ", Crossing: " << i << std::endl;
            int p1 = gen_archive.probabilisticIndividual();
            int p2 = gen_archive.probabilisticIndividual();
            if (p1 == p2) continue;
            std::cout << "cross " << p1 << " + " << p2 << " /" << gen_archive.getSize() << std::endl;
            std::shared_ptr<LabelingIndividual> child = std::make_shared<LabelingIndividual>(
                                                                LabelingIndividual(*gen_archive.getIndiv(p1), 
                                                                                   *gen_archive.getIndiv(p2)));
            child->updateChartsAndTPs();
            child->repairUnspikeLabeling();
            child->updateChartsAndTPs(true); // TPs not needed?
            child->updateTimestamps(); // not needed?
            double new_score = evaluator.evaluate(*child);
            gen_archive.insert(child, new_score);
        }

        // Insert generation into general archive
        for (int id=0; id<gen_archive.getSize(); id++){
            archive.insert(gen_archive.getIndiv(id), gen_archive.getScore(id));
        }

        auto time_after_archive = std::chrono::steady_clock::now();
        double time_insert_archive = measureTime(time_after_mutations, time_after_archive);
        std::cout << "End of generation, additional time:" << std::endl;
        std::cout << "\tGeneral archive \t" << time_insert_archive << std::endl;
    }

    LabelingIndividual final_indiv = *archive.getIndiv(0);
    final_indiv.updateChartsAndTPs(true);
    def_V = qle->computeDeformedV(final_indiv.getLabeling());

    coloredPrint("Final validity: " + std::to_string(final_indiv.invalidityScore()), "cyan");

    Eigen::MatrixXd threshold_colors = qle->distoAboveThreshold(final_indiv.getLabeling(), 100.0);

    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    int orig_id = viewer.data_list[0].id;
    int hud_id = viewer.data_list[1].id;

    viewer.data(orig_id).set_mesh(V, F);

    bool flip_label_normals = false;
    int arch_elem = 0;

    std::shared_ptr<LabelingIndividual> best_indiv = archive.getIndiv(0); 
    std::shared_ptr<LabelingIndividual> displayed_indiv = best_indiv;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data(hud_id).clear_edges();
        viewer.data(orig_id).clear_edges();
        viewer.data(orig_id).clear_points();

        Eigen::MatrixXd colors = colorsFromFlagging(displayed_indiv->getLabeling());
        viewer.data(orig_id).set_colors(colors);
        
        std::vector<Eigen::MatrixXd> vec_border_begs, vec_border_ends;
        std::vector<Eigen::RowVector3d> vec_border_colors;
        displayed_indiv->getBordersViz(vec_border_begs, vec_border_ends, vec_border_colors);
        for (int i=0; i<vec_border_begs.size(); i++){
            viewer.data(hud_id).add_edges(vec_border_begs[i], vec_border_ends[i], vec_border_colors[i]);
        }

        Eigen::MatrixXd tp_points = displayed_indiv->getTurningPointsMat();
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
                displayed_indiv = ancestor;
                viewer.data(orig_id).set_mesh(V, F);
                //Eigen::MatrixXd colors = colorsFromFlagging(labeling_init);
                //viewer.data(orig_id).set_colors(colors);
                updateViz();
            }
            if (ImGui::Button("Labeling", ImVec2(-1, 0))){
                displayed_indiv = best_indiv;
                viewer.data(orig_id).set_mesh(V, F);
                updateViz();
            }
            if (ImGui::Button("Fastbndpolycube", ImVec2(-1, 0))){
                viewer.data(orig_id).set_mesh(def_V, F);
            }
            if (ImGui::Button("Threshold dist", ImVec2(-1, 0))){
                viewer.data(orig_id).set_colors(threshold_colors);
            }

            if (ImGui::SliderInt("Archive labeling", &arch_elem, 0, archive.getSize() - 1)){
                displayed_indiv = archive.getIndiv(arch_elem);
                viewer.data(orig_id).set_mesh(V, F);
                updateViz();
            }

            if (ImGui::Button("Save labeling to folder", ImVec2(-1, 0))){
                std::string folder = igl::file_dialog_save();
                Eigen::VectorXi save_labeling = final_indiv.getLabeling();
                saveFlagging(folder + "labeling.txt", save_labeling);
                saveFlaggingOnTets(folder + "labeling_on_tets.txt", folder + "tris_to_tets.txt", save_labeling);
            }

            if (make_checkbox("Show timestamps", viewer.data(orig_id).show_custom_labels)){
                viewer.data(orig_id).clear_labels();
                Eigen::VectorXi timestamps = final_indiv.getTimestamps();
                Eigen::MatrixXd N;
                igl::per_face_normals(V, F, N);
                if (flip_label_normals) N = -N;
                double l_avg = igl::avg_edge_length(V, F);
                for (int i=0; i<F.rows(); i++){
                    Eigen::RowVector3d p = (V.row(F(i,0)) + V.row(F(i,1)) + V.row(F(i,2)))/3.0;
                    p -= N.row(i) * l_avg / 3.0;
                    viewer.data(orig_id).add_label(p, std::to_string(timestamps(i)));
                }
            }

            ImGui::Checkbox("Flip label normals", &flip_label_normals);

            

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