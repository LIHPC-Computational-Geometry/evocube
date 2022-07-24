#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/file_dialog_save.h>
#include <ctime>
#include <random>
#include <queue>
#include <nlohmann/json.hpp>

#include "graphcut_labeling.h"
#include "flagging_utils.h"
#include "logging.h"
#include "chart.h"
#include "quick_label_ev.h"
#include "evocube.h"
#include "evaluator.h"
#include "labeling_individual.h"
#include "archive.h"

//#define DEBUG_EVOCUBE

//#define PRINT_EVOCUBE_TIMINGS
//#define VERBOSE_BIRTH_DEATH

#define READ_DOT_MESH
#ifdef READ_DOT_MESH
#include "mesh_io.h" // requires tetgen
#endif


int main(int argc, char *argv[]){
    std::string input_tris = "../data/bunny/boundary.obj";
    if (argc > 1) input_tris = argv[1];

    bool show_ui = (argc <= 2);

    //std::srand(std::time(nullptr));
    std::srand(1);

    auto time_before_init = std::chrono::steady_clock::now();

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //igl::readOBJ("../data/Basic/B23/boundary.obj", V, F);
    igl::readOBJ(input_tris, V, F);
    //readDotMeshTri(input_tris, V, F);

    std::shared_ptr<Evocube> evo = std::make_shared<Evocube>(Evocube(V, F));
    std::shared_ptr<const QuickLabelEv> qle = std::make_shared<const QuickLabelEv>(QuickLabelEv(V, F));
    const Evaluator evaluator(evo, qle); 

    // ---- GENERATE INITIAL SOLUTION ---- //

    int compact_coeff = 1;
    int fidelity_coeff = 3;
    Eigen::VectorXi labeling_init = graphcutFlagging(V, F, evo->N_, evo->TT_, compact_coeff, fidelity_coeff);
    std::shared_ptr<LabelingIndividual> ancestor = std::make_shared<LabelingIndividual>(LabelingIndividual(evo, qle, labeling_init));

    ancestor->updateChartsAndTPs();
    //ancestor->repairHighValenceCorner();
    ancestor->updateChartsAndTPs();
    ancestor->repairOppositeLabels();
    ancestor->updateChartsAndTPs();

    ancestor->repairUnspikeLabeling();
    ancestor->updateChartsAndTPs(true);
    ancestor->updateTimestamps();
    
    double current_score = evaluator.evaluate(*ancestor);
    Eigen::MatrixXd def_V;// = qle->computeDeformedV(ancestor.getLabeling());

    Archive<std::shared_ptr<LabelingIndividual>> archive(10);
    archive.insert(ancestor, current_score);

    auto pickMutation = [](double progress_percentage){
        double r = static_cast<double>(std::rand() % 1000)/1000.0;
        double prob0 = 0.5 - 0.2 * progress_percentage;
        double prob1 = 0.0 + 0.2 * progress_percentage;
        double prob2 = 0.3 - 0.0 * progress_percentage;
        double prob3 = 0.2 - 0.0 * progress_percentage;
        if (r < prob0) return 0;
        if (r < prob0 + prob1) return 1;
        if (r < prob0 + prob1 + prob2) return 2;
        else return 3;
    };

    double sum_time_create_indiv = 0;
    double sum_time_cross = 0;
    double sum_time_mut = 0;
    double sum_time_chart = 0;
    double sum_time_archive = 0;
    double sum_time_eval = 0;


    // ---- GENETIC OPTIMIZATION ---- //

    auto time_before_evocube = std::chrono::steady_clock::now();

    int n_generations = 40;
    int max_mut = 100;

    int convergence_stop = 3; // If score doesn't improve for convergence_stop consecutive iterations, stop
    int convergence_count = 0;
    double previous_best_score = -1; 

    for (int generation=0; generation<n_generations; generation++){

        evo->timestamp_ ++;

        Archive<std::shared_ptr<LabelingIndividual>> gen_archive(5);

        std::vector<std::shared_ptr<LabelingIndividual>> new_gen;
        std::vector<double> new_scores;

        new_gen.resize(max_mut);
        new_scores.resize(max_mut);

        if (archive.bestScore() != previous_best_score){
            convergence_count = 0;
            previous_best_score = archive.bestScore(); 
        }
        else {
            convergence_count ++;
            if (convergence_count > convergence_stop){
                coloredPrint("Reached convergence, stopping labeling optimization", "green");
                n_generations = generation;
                break;
            }
        }

        // -- Generate new individuals through MUTATIONS -- //
        #pragma omp parallel for
        for (int i=0; i<max_mut; i++){
            std::cout << "Generation " << generation << ", Mutation: " << i << std::endl;

            auto time_new_indiv = std::chrono::steady_clock::now();    
            int pick = archive.probabilisticIndividual();
            std::shared_ptr<LabelingIndividual> new_indiv = std::make_shared<LabelingIndividual>(*archive.getIndiv(pick));
            new_indiv->updateChartsAndTPs(true); // OPTIM: could be skipped if we kept info from parent

            auto time_before_mutation = std::chrono::steady_clock::now();     

            int mutation_type = pickMutation(static_cast<double>(generation)/n_generations);
            if (mutation_type == 0) new_indiv->mutationVertexGrow();
            if (mutation_type == 1) new_indiv->mutationRemoveChart();
            if (mutation_type == 2) new_indiv->mutationGreedyPath();
            if (mutation_type == 3) new_indiv->mutationBorderGrow();

            auto time_after_mutation = std::chrono::steady_clock::now();

            new_indiv->updateChartsAndTPs();
            new_indiv->repairUnspikeLabeling();
            new_indiv->updateChartsAndTPs(true);
            new_indiv->updateTimestamps();

            auto time_after_charts = std::chrono::steady_clock::now();
            
            double new_score = evaluator.evaluate(*new_indiv);
            new_gen[i] = new_indiv;
            new_scores[i] = new_score;

            auto time_after_eval = std::chrono::steady_clock::now();

            double time_create_indiv = measureTime(time_new_indiv, time_before_mutation);
            double time_mut = measureTime(time_before_mutation, time_after_mutation);
            double time_chart = measureTime(time_after_mutation, time_after_charts);
            double time_eval = measureTime(time_after_charts, time_after_eval);
            double totaltime = measureTime(time_new_indiv, time_after_eval);

            sum_time_create_indiv += time_create_indiv;
            sum_time_mut += time_mut;
            sum_time_chart += time_chart;
            sum_time_eval += time_eval;

            #ifdef PRINT_EVOCUBE_TIMINGS
            std::cout << "Individual mutation, time repartition:" << std::endl;
            std::cout << "\tCreation \t" << time_create_indiv << std::endl;
            std::cout << "\tMutation \t" << time_mut << std::endl;
            std::cout << "\tCharts   \t" << time_chart << std::endl;
            std::cout << "\tEval.    \t" << time_eval << std::endl;
            std::cout << "\tTotal:   \t" << totaltime << std::endl;
            #endif
        }

        // Insert new indivs into archive (not done in the loop because it's not //)
        for (int i=0; i<new_gen.size(); i++){
            gen_archive.insert(new_gen[i], new_scores[i]);
        }

        auto time_after_mutations = std::chrono::steady_clock::now();

        // -- Generate new individuals through CROSSING -- //
        int n_cross = 10;
        for (int i=0; i<n_cross; i++){
            std::cout << "Generation " << generation << ", Crossing: " << i << std::endl;

            auto time_before_cross = std::chrono::steady_clock::now();    
            int p1 = gen_archive.probabilisticIndividual();
            int p2 = gen_archive.probabilisticIndividual();
            if (p1 == p2) continue;
            std::shared_ptr<LabelingIndividual> child = std::make_shared<LabelingIndividual>(
                                                                LabelingIndividual(*gen_archive.getIndiv(p1), 
                                                                                   *gen_archive.getIndiv(p2)));

            auto time_after_cross = std::chrono::steady_clock::now();    

            child->updateChartsAndTPs();
            child->repairUnspikeLabeling();
            child->updateChartsAndTPs(true); // TPs not needed?
            child->updateTimestamps(); // not needed?

            auto time_after_charts = std::chrono::steady_clock::now();   

            double new_score = evaluator.evaluate(*child);
            auto time_after_eval = std::chrono::steady_clock::now();   

            gen_archive.insert(child, new_score);
            auto time_after_archive = std::chrono::steady_clock::now();   


            sum_time_cross += measureTime(time_before_cross, time_after_cross);
            sum_time_chart += measureTime(time_after_cross, time_after_charts);
            sum_time_eval += measureTime(time_after_charts, time_after_eval);
            sum_time_archive += measureTime(time_after_eval, time_after_archive);
        }

        auto time_before_archive = std::chrono::steady_clock::now();
        // Insert generation into general archive
        for (int id=0; id<gen_archive.getSize(); id++){
            archive.insert(gen_archive.getIndiv(id), gen_archive.getScore(id));
        }
        auto time_after_archive = std::chrono::steady_clock::now();
        sum_time_archive += measureTime(time_before_archive, time_after_archive);
        //std::cout << "End of generation, additional time:" << std::endl;
        //std::cout << "\tGeneral archive \t" << time_insert_archive << std::endl;
    }

    auto time_after_evocube = std::chrono::steady_clock::now();


    // ---- FINAL SOLUTION ---- //
    std::shared_ptr<LabelingIndividual> final_indiv = archive.getIndiv(0);

    final_indiv->updateChartsAndTPs();
    final_indiv->removeChartsWithTooFewNeighbors();
    final_indiv->updateChartsAndTPs();
    final_indiv->repairHighValenceCorner();
    final_indiv->updateChartsAndTPs();
    final_indiv->repairOppositeLabels();
    final_indiv->updateChartsAndTPs();
    final_indiv->updateChartsAndTPs(true);
    def_V = qle->computeDeformedV(final_indiv->getLabeling());

    coloredPrint("Final validity: " + std::to_string(final_indiv->invalidityScore()), "cyan");

    Eigen::MatrixXd threshold_colors = qle->distoAboveThreshold(final_indiv->getLabeling(), 100.0);

    auto time_after_post = std::chrono::steady_clock::now();
    double time_init_evo = measureTime(time_before_init, time_before_evocube);
    double time_evocube = measureTime(time_before_evocube, time_after_evocube);
    double time_post_evo = measureTime(time_after_evocube, time_after_post);
    coloredPrint("Evocube time: " + std::to_string(time_evocube), "cyan");

    
    if (!show_ui) {
        std::string save_path = argv[2];
        Eigen::VectorXi save_labeling = final_indiv->getLabeling();

        std::cout << "Saving to:" << save_path << std::endl;
        saveFlagging(save_path + "/labeling.txt", save_labeling);
        saveFlaggingOnTets(save_path + "/labeling_on_tets.txt", save_path + "/tris_to_tets.txt", save_labeling);

        saveFlagging(save_path + "/labeling_init.txt", labeling_init);
        igl::writeOBJ(save_path + "/fast_polycube_surf.obj", def_V, F);

        std::string logs_path = save_path + "/logs.json";
        final_indiv->fillIndivLogInfo(logs_path, "LabelingFinal");
        evaluator.fillIndivLogInfo(logs_path, *final_indiv, "LabelingFinal");

        std::shared_ptr<LabelingIndividual> graph_cut_indiv = std::make_shared<LabelingIndividual>(LabelingIndividual(evo, qle, labeling_init));
        graph_cut_indiv->updateChartsAndTPs(true);
        graph_cut_indiv->fillIndivLogInfo(logs_path, "LabelingGraphCut");
        evaluator.fillIndivLogInfo(logs_path, *graph_cut_indiv, "LabelingGraphCut");

        Eigen::VectorXi normal_labeling = normalFlagging(V, F);
        std::shared_ptr<LabelingIndividual> normal_indiv = std::make_shared<LabelingIndividual>(LabelingIndividual(evo, qle, normal_labeling));
        normal_indiv->updateChartsAndTPs(true);
        normal_indiv->fillIndivLogInfo(logs_path, "LabelingNormal");
        evaluator.fillIndivLogInfo(logs_path, *normal_indiv, "LabelingNormal");

        auto start_time = std::chrono::system_clock::now();
        std::time_t start_timet = std::chrono::system_clock::to_time_t(start_time);
        fillLogInfo("FinishTime", logs_path, std::ctime(&start_timet));

        // ALL THESE ARE IN CPU TIME! Since it's //, the sum is > to time_evocube
        fillLogInfo("Timing", "CreateIndiv", logs_path, sum_time_create_indiv);
        fillLogInfo("Timing", "Cross", logs_path, sum_time_cross);
        fillLogInfo("Timing", "Mutations", logs_path, sum_time_mut);
        fillLogInfo("Timing", "ChartsAndTps", logs_path, sum_time_chart);
        fillLogInfo("Timing", "Archive", logs_path, sum_time_archive);
        fillLogInfo("Timing", "Eval", logs_path, sum_time_eval);

        // Real-world time
        fillLogInfo("Timing", "PreGenetics", logs_path, time_init_evo);
        fillLogInfo("Timing", "Genetics", logs_path, time_evocube);
        fillLogInfo("Timing", "PostGenetics", logs_path, time_post_evo);

        fillLogInfo("#generations", logs_path, std::to_string(n_generations));
        fillLogInfo("#mutations_per_gen", logs_path, std::to_string(max_mut));

        fillLogInfo("GraphCutParams", "CompactCoeff", logs_path, compact_coeff);
        fillLogInfo("GraphCutParams", "FidelityCoeff", logs_path, fidelity_coeff);
        evo->fillMeshLogInfo(logs_path); 
    }
    else { // ---- VISUALIZATION ---- //
        igl::opengl::glfw::Viewer viewer;
        viewer.append_mesh();
        int orig_id = viewer.data_list[0].id;
        int hud_id = viewer.data_list[1].id;

        viewer.data(orig_id).set_mesh(V, F);

        bool flip_label_normals = false;
        bool debug_option = false;
        bool show_indiv = true;

        std::shared_ptr<LabelingIndividual> best_indiv = final_indiv; 
        std::shared_ptr<LabelingIndividual> displayed_indiv = best_indiv;

        igl::opengl::glfw::imgui::ImGuiMenu menu;
        menu.callback_draw_viewer_window = []() {};
        viewer.plugins.push_back(&menu);

        auto updateViz = [&](){
            viewer.data(hud_id).clear_edges();
            viewer.data(hud_id).clear_points();
            viewer.data(orig_id).clear_edges();
            viewer.data(orig_id).clear_points();

            Eigen::MatrixXd colors = colorsFromFlagging(displayed_indiv->getLabeling());
            viewer.data(orig_id).set_colors(colors);

            if (show_indiv){
                std::vector<Eigen::MatrixXd> vec_border_begs, vec_border_ends;
                std::vector<Eigen::RowVector3d> vec_border_colors;
                displayed_indiv->getBordersViz(vec_border_begs, vec_border_ends, vec_border_colors);
                for (int i=0; i<vec_border_begs.size(); i++){
                    viewer.data(hud_id).add_edges(vec_border_begs[i], vec_border_ends[i], vec_border_colors[i]);
                }

                Eigen::MatrixXd tp_points = displayed_indiv->getTurningPointsMat();
                viewer.data(hud_id).add_points(tp_points, Eigen::RowVector3d(1.0, 1.0, 0.0));
            }
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
                    show_indiv = true;
                    displayed_indiv = ancestor;
                    viewer.data(orig_id).set_mesh(V, F);
                    //Eigen::MatrixXd colors = colorsFromFlagging(labeling_init);
                    //viewer.data(orig_id).set_colors(colors);
                    updateViz();
                }
                if (ImGui::Button("Labeling", ImVec2(-1, 0))){
                    show_indiv = true;
                    displayed_indiv = best_indiv;
                    viewer.data(orig_id).set_mesh(V, F);
                    updateViz();
                }
                if (ImGui::Button("Fast boundary polycube", ImVec2(-1, 0))){
                    show_indiv = false;
                    displayed_indiv = best_indiv;
                    viewer.data(orig_id).set_mesh(def_V, F);
                    updateViz();
                }

                ImGui::Separator();
                
                make_checkbox("Show mesh", viewer.data(orig_id).show_lines);
                ImGui::Checkbox("Debug mode", &debug_option);
                
                if (debug_option){
                    if (ImGui::Button("Threshold dist", ImVec2(-1, 0))){
                        viewer.data(orig_id).set_colors(threshold_colors);
                    }

                    if (ImGui::Button("Save labeling to folder", ImVec2(-1, 0))){
                        std::string folder = igl::file_dialog_save();
                        folder = folder.substr(0, folder.find_last_of("/\\") + 1);
                        Eigen::VectorXi save_labeling = final_indiv->getLabeling();
                        saveFlagging(folder + "labeling.txt", save_labeling);
                        saveFlaggingOnTets(folder + "labeling_on_tets.txt", folder + "tris_to_tets.txt", save_labeling);
                    }

                    if (make_checkbox("Show timestamps", viewer.data(orig_id).show_custom_labels)){
                        viewer.data(orig_id).clear_labels();
                        Eigen::VectorXi timestamps = final_indiv->getTimestamps();
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


                    if (ImGui::Button("Quick save labeling", ImVec2(-1, 0))){
                        std::string folder = input_tris.substr(0, input_tris.find_last_of("/\\") + 1);
                        folder = folder.substr(0, folder.find_last_of("/\\") + 1);
                        Eigen::VectorXi save_labeling = final_indiv->getLabeling();
                        saveFlagging(folder + "/labeling.txt", save_labeling);
                        saveFlaggingOnTets(folder + "/labeling_on_tets.txt", folder + "/tris_to_tets.txt", save_labeling);
                    }
                    
                    if (ImGui::Button("Quick hexex")){
                        std::string folder = input_tris.substr(0, input_tris.find_last_of("/\\") + 1);
                        std::string tet = folder + "/tetra.mesh";
                        std::string label_tet = folder + "/labeling_on_tets.txt";
                        std::string hexes = folder + "/hexes.mesh";
                        double scale = 1.4;
                        std::string command = "./polycube_withHexEx  " + tet + " " + label_tet + " " + hexes + " 1.4";
                        std::cout << "Running command: " << command << std::endl;
                        int result_command = system(command.c_str());
                    }
                }
                ImGui::End();
            }
        };

        updateViz();
        viewer.data(hud_id).line_width = 10.0;
        viewer.core().lighting_factor = 0.0;
        viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
        viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
        viewer.launch();

    }

    return 0;
}