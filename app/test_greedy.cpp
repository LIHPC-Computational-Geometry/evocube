#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/per_face_normals.h>
#include <igl/avg_edge_length.h>
#include <queue>
#include <random>

#include "graphcut_labeling.h"
#include "logging.h"
#include "chart.h"
#include "quick_label_ev.h"

void growMutation(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                  const Eigen::MatrixXi& TT, const Eigen::MatrixXd& dists, 
                  int introduced_label, int start_tri_id,
                  double threshold_dist){
    
    int initial_label = old_labeling(start_tri_id);
    new_labeling = old_labeling;

    std::unordered_map<int, double> seen;
    std::queue<int> queue;

    queue.push(start_tri_id);
    seen.insert(std::make_pair(start_tri_id, 0.0));

    while (!queue.empty()){
        int current = queue.front(); 
        #ifdef DEBUG_MUTATIONS
        std::cout << "queue.size() " << queue.size() << std::endl;
        std::cout << "seen.size() " << seen.size() << std::endl;
        std::cout << "current " << current << std::endl;
        #endif
        queue.pop();
        new_labeling(current) = introduced_label;
        for (int neigh=0; neigh<3; neigh++){
            int neigh_id = TT(current, neigh);
            if (neigh_id == -1) coloredPrint("Error, invalid neighboring triangle", "red");
            if (old_labeling(neigh_id) != initial_label) continue;
            double dist = seen[current] + dists(current, neigh);
            if (dist > threshold_dist) continue;
            if (seen.find(neigh_id) != seen.end()){ // face seen before
                if (seen.count(neigh_id) > 1) coloredPrint("Error, several = keys", "red");
                if (seen[neigh_id] < dist) continue; // old path was better, skip
                // new path is better, update dist and reinsert
                seen.find(neigh_id)->second = dist;
            }
            else { // face not seen before

                #ifdef DEBUG_MUTATIONS
                //std::cout << "New face " << dist << std::endl;
                #endif
                seen.insert(std::pair<int, double>(neigh_id, dist));
            }
            queue.push(neigh_id);
        }
    }
}

void removeChartMutation(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                         const Eigen::VectorXi& charts, 
                         const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, 
                         const Eigen::MatrixXd& N, 
                         const Eigen::MatrixXi& TT, int chart_id){

    Eigen::VectorXi locked_flags = old_labeling;
    Eigen::VectorXi forbidden_flags = Eigen::VectorXi::Constant(F.rows(), -1);
    #pragma omp parallel for
    for (int i = 0; i < F.rows(); i ++){
        if (charts(i) == chart_id) {
            locked_flags(i) = -1;
            forbidden_flags(i) = old_labeling(i);
        }
    }

    int compact_coeff = 3;
    int fidelity_coeff = 1;
    new_labeling = graphcutFlagging(V, F, N, TT, locked_flags, forbidden_flags, compact_coeff, fidelity_coeff);
}

// might want to add normal angle to distance

// Note: such a distance will lead to "round" regions
// but we're doing polycubes
// how to generate square regions?
// Maybe by considering max(dist_x, dist_y, dist_z) instead?
void triangle_triangle_dist(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, 
                            const Eigen::MatrixXi& TT, 
                            Eigen::MatrixXd& dists){
    dists.resize(F.rows(), 3);

    dists = Eigen::MatrixXd::Constant(F.rows(), 3, -1);

    for (int f_id = 0; f_id < F.rows(); f_id++){
        Eigen::RowVector3d G1 = (V.row(F(f_id, 0)) + V.row(F(f_id, 1)) + V.row(F(f_id, 2))) / 3.0;
        for (int neigh = 0; neigh < 3; neigh ++){
            int neigh_id = TT(f_id, neigh);
            if (neigh_id == -1) {
                coloredPrint("Triangle doesn't have a neighbor", "red");
                continue;
            }
            Eigen::RowVector3d G2 = (V.row(F(neigh_id, 0)) + V.row(F(neigh_id, 1)) + V.row(F(neigh_id, 2))) / 3.0;
            // middle point on shared edge:
            Eigen::RowVector3d mp = (V.row(F(f_id, neigh)) + V.row(F(f_id, (neigh + 1) % 3))) / 2.0;
            dists(f_id, neigh) = (G1 - mp).norm() + (mp - G2).norm();
        }
    }
}

double invalidChartsScore(const std::vector<std::vector<int>>& chart_adj){
    int sum = 0;
    for (int i = 0; i < chart_adj.size(); i++){
        if (chart_adj[i].size() < 4){
            sum += 4 - chart_adj[i].size();
        }
    }
    return static_cast<double>(sum);
}

int main(int argc, char *argv[]){

    Eigen::MatrixXd V, N, def_V;
    Eigen::MatrixXi F, TT;
    igl::readOBJ("../data/bunny/boundary.obj", V, F);
    igl::triangle_triangle_adjacency(F, TT);
    igl::per_face_normals(V, F, N);
    
    int compact_coeff = 1;
    int fidelity_coeff = 3;
    Eigen::VectorXi labeling = graphcutFlagging(V, F, N, TT, compact_coeff, fidelity_coeff);

    Eigen::VectorXi charts;
    int n_charts = updateCharts(TT, labeling, charts);
    std::cout << "n_charts " << n_charts << std::endl;
    std::vector<std::vector<int>> adj = chartAdjacency(TT, charts);

    Eigen::MatrixXd dists;
    triangle_triangle_dist(V, F, TT, dists);

    double l_avg = igl::avg_edge_length(V, F);
    int mut_count = 0;

    QuickLabelEv qle(V, F);

    double current_score = qle.evaluate(labeling) + 10000.0 * invalidChartsScore(adj);
    def_V = qle.getDefV();

    for (int i=0; i<10; i++){
        std::cout << "Attempt: " << i << std::endl;
        Eigen::VectorXi new_labeling;

	    auto time_before_mutation = std::chrono::steady_clock::now();        
        
        int mutation_type = 1;
        if (mutation_type == 0){
            int start_tri_id = std::rand() % TT.rows();
            int introduced_label = std::rand() % 6;
            growMutation(labeling, new_labeling, TT, dists, introduced_label, start_tri_id, 5 * l_avg);
        }

        if (mutation_type == 1){
            int chart_id = std::rand() % (charts.maxCoeff() + 1);
            removeChartMutation(labeling, new_labeling, charts, V, F, N, TT, chart_id);
        }

	    auto time_after_mutation = std::chrono::steady_clock::now();

        Eigen::VectorXi new_charts;
        int new_n_charts = updateCharts(TT, new_labeling, new_charts);
        std::vector<std::vector<int>> new_adj = chartAdjacency(TT, new_charts);

	    auto time_after_charts = std::chrono::steady_clock::now();
        double time_mut = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_mutation - time_before_mutation).count()) / 1000;
        double time_chart = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_charts - time_after_mutation).count()) / 1000;
        double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (time_after_charts - time_before_mutation).count()) / 1000;
	    std::cout << "Mutation + chart time: " << time_mut << " + " << time_chart << " = " << totaltime << "sec." << std::endl;

        double new_score = qle.evaluate(new_labeling) + 10000.0 * invalidChartsScore(new_adj);

        if (new_score < current_score){
            current_score = new_score; 
            labeling = new_labeling;
            def_V = qle.getDefV();
            charts = new_charts;
            adj = new_adj;
            mut_count ++;
        }
    }


    std::vector<std::vector<int>> borders;
    std::vector<std::pair<int, int>> patches_per_border;
    computeBoundaries(TT, F, charts, borders, patches_per_border);
    std::vector<Eigen::MatrixXd> vec_border_begs; //(border.size(), 3);
    std::vector<Eigen::MatrixXd> vec_border_ends; //(border.size(), 3);

    Eigen::VectorXi per_chart_labels = perChartLabels(charts, labeling);
    std::vector<Eigen::RowVector3d> boundary_axis;
    for (int i=0; i<borders.size(); i++){
        Eigen::RowVector3d axis;
        axis << 0, 0, 0;
        int label1 = per_chart_labels[patches_per_border[i].first];
        int label2 = per_chart_labels[patches_per_border[i].second];
        axis(3 - label1/2 - label2/2) = 1.0;
        boundary_axis.push_back(axis);
    }

    for (int b=0; b<borders.size(); b++){
        Eigen::MatrixXd border_begs(borders[b].size() - 1, 3);
        Eigen::MatrixXd border_ends(borders[b].size() - 1, 3);
        for (int i=0; i<borders[b].size() - 1; i++){
            border_begs.row(i) = V.row(borders[b][i]);
            border_ends.row(i) = V.row(borders[b][i + 1]);
        }
        vec_border_begs.push_back(border_begs);
        vec_border_ends.push_back(border_ends);
    }

    std::vector<std::vector<int>> vec_tps;
    int n_tps = 0;
    for (int i=0; i<borders.size(); i++){
        std::vector<int> tps = graphcutTurningPoints(borders[i], V, boundary_axis[i]);
        vec_tps.push_back(tps);
        n_tps += tps.size();
    }

    std::cout << "#borders: " << borders.size() << std::endl;
    std::cout << "turning points: " << n_tps << std::endl;
    std::cout << "Mut accepted: " << mut_count << std::endl;

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

        Eigen::MatrixXd colors = colorsFromFlagging(labeling);
        viewer.data(orig_id).set_colors(colors);

        for (int i=0; i<vec_border_begs.size(); i++){
            Eigen::RowVector3d bnd_color = Eigen::RowVector3d::Random();
            bnd_color = boundary_axis[i];
            viewer.data(hud_id).add_edges(vec_border_begs[i], vec_border_ends[i], bnd_color);
        }

        for (int i=0; i<vec_tps.size(); i++){
            for (int j=0; j<vec_tps[i].size(); j++){
                viewer.data(hud_id).add_points(V.row(borders[i][vec_tps[i][j]]), Eigen::RowVector3d(1.0, 1.0, 0.0));
            }
        }
    };

    menu.callback_draw_custom_window = [&]() {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(350, -1), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("IGL")) {
            if (ImGui::Button("Labeling", ImVec2(-1, 0))){
                viewer.data(orig_id).set_mesh(V, F);
            }
            if (ImGui::Button("Fastbndpolycube", ImVec2(-1, 0))){
                viewer.data(orig_id).set_mesh(def_V, F);
            }
            
            if (ImGui::Button("Grow", ImVec2(-1, 0))){
                mut_count ++;
                std::cout << "Mutation: " << mut_count << std::endl;
                Eigen::VectorXi old_labeling = labeling;
                int start_tri_id = std::rand() % TT.rows();
                int introduced_label = (labeling(start_tri_id) + 2) % 6;
                growMutation(old_labeling, labeling, TT, dists, introduced_label, start_tri_id, 5 * l_avg);
                updateViz();
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

    return 0;
}