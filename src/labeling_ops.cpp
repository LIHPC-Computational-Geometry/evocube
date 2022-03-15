#include "labeling_ops.h"
#include <queue>
#include <vector>
#include <unordered_map>

#include "graphcut_labeling.h"
#include "flagging_utils.h"
#include "logging.h"
#include "chart.h"

//#define DEBUG_MUTATIONS

void fixHighValenceCorners(const Eigen::MatrixXd& dists, 
                           const Eigen::MatrixXi& TT, 
                           const std::vector<std::vector<int>>& VT, 
                           const Eigen::VectorXi& charts,
                           const std::vector<std::vector<int>>& borders,
                           double l_avg,
                           const QuickLabelEv& qle,
                           Eigen::VectorXi& labeling){

    std::vector<int> corner_ids;
    for (auto b: borders){
        corner_ids.push_back(b[0]);
        corner_ids.push_back(b[b.size() - 1]);
    }

    std::vector<int> corner_uni = corner_ids; // corners without duplicates
    std::sort(corner_uni.begin(), corner_uni.end());
    corner_uni.erase(std::unique(corner_uni.begin(), corner_uni.end()), corner_uni.end());

    for (int c: corner_uni){
        #ifdef DEBUG_LABELING_REPAIRS
        if (std::count(corner_ids.begin(), corner_ids.end(), c) < 3){
            coloredPrint("ERROR: Incorrect border connectivity", "red");
        }
        if (std::count(corner_ids.begin(), corner_ids.end(), c) == 3){
            coloredPrint("Correct corner: " + std::to_string(c), "cyan");
        }
        #endif

        if (std::count(corner_ids.begin(), corner_ids.end(), c) > 3){
            coloredPrint("High valence corner detected: id " + std::to_string(c), "green");

            std::vector<int> start_tri_ids = VT[c];
        
            Eigen::VectorXi best_labeling = labeling;
            double best_score = -1;

            for (int introduced_label=0; introduced_label<6; introduced_label++){
                for (int dist_coeff = 1; dist_coeff<10; dist_coeff += 2){
                    Eigen::VectorXi labeling_attempt = labeling;
                    double threshold_dist = dist_coeff * l_avg;
                    growFromTris(labeling, labeling_attempt, TT, dists, introduced_label, start_tri_ids, threshold_dist, false);
                    
                    double score = qle.evaluate(labeling_attempt);
                    if (best_score == -1 || score < best_score){
                        best_score = score;
                        best_labeling = labeling_attempt;
                    }
                }
            }

            labeling = best_labeling;
        }
    }
}

void directionalPathMutation(const std::vector<std::vector<int>>& VT, 
                             const Eigen::MatrixXi& TT, 
                             const Eigen::MatrixXd& V, 
                             const Eigen::MatrixXi& F, 
                             const Eigen::MatrixXd& dists,
                             const std::vector<std::vector<int>>& borders,
                             const std::vector<std::pair<int, int>>& patches_per_border,
                             const Eigen::VectorXi& charts,
                             const Eigen::RowVector3d& direction,
                             int starting_vertex,  
                             int chart_id,
                             double distance_thresold,
                             int introduced_label,
                             Eigen::VectorXi& labeling){

    std::cout << "directionalPathMutation..." << std::endl;
    // precompute: sorted vector of patch boundary vertices
    std::vector<int> border_vs;
    for (int i=0; i<borders.size(); i++){
        if (patches_per_border[i].first != chart_id &&
            patches_per_border[i].second != chart_id) continue;
        border_vs.insert(border_vs.end(), borders[i].begin(), borders[i].end());
    }
    std::sort(border_vs.begin(), border_vs.end());
    border_vs.erase(std::unique(border_vs.begin(), border_vs.end()), border_vs.end());

    auto isBorderVertex = [&border_vs](int v_id){
        for (int i: border_vs){
            if (i == v_id) return true; // OPTIM dichotomy
        }
        return false;
    };

    if (!isBorderVertex(starting_vertex)){
        coloredPrint("Directional path starting from a non-border vertex", "red");
    }

    // GREEDY PICK

    std::vector<int> path = {starting_vertex};
    int current = starting_vertex;

    auto alreadyPicked = [&path](int v_id){
        for (int i: path) if (v_id == i) return true;
        return false;
    };

    while (!isBorderVertex(current) || current == starting_vertex){
        std::vector<int> candidates;
        for (int tri: VT[current]){
            if (charts(tri) != chart_id) continue;
            for (int k=0; k<3; k++){
                int cand = F(tri, k);
                if (alreadyPicked(cand)) continue;
                candidates.push_back(cand);
            }
        }

        if (candidates.size() == 0){
            coloredPrint("No more candidates in greedy path selection", "red");
            return;
        }

        double best = -1;
        int best_cand = -1;
        for (int cand: candidates){
            Eigen::RowVector3d path_vec = V.row(cand) - V.row(starting_vertex);
            double score = path_vec.dot(direction) / (path_vec.norm() * direction.norm());
            if (score < best || best_cand == -1){
                best = score;
                best_cand = cand;
            }
        }

        current = best_cand;
        path.push_back(best_cand);
    }

    

    // Once we have our path of vertices, we get a path of triangles using
    // VT
    // and grow around that region
    std::vector<int> start_tri_ids = {};
    for (int v_id: path){
        for (int tri_id: VT[v_id]){
            if (charts(tri_id) == chart_id){
                start_tri_ids.push_back(tri_id);
            }
        }
    }

    std::sort(start_tri_ids.begin(), start_tri_ids.end());
    start_tri_ids.erase(std::unique(start_tri_ids.begin(), start_tri_ids.end()), start_tri_ids.end());

    Eigen::VectorXi old_labeling = labeling;
    growFromTris(old_labeling, labeling, TT, dists, introduced_label, start_tri_ids, distance_thresold);

    #ifdef DEBUG_LABELING_MUTATIONS
    std::cout << "path: ";
    for (int i: path) std::cout << i << " ";
    std::cout << std::endl;
    std::cout << "end of directionalPathMutation" << std::endl;
    #endif
}

void unspikeLabeling(const Eigen::MatrixXi& TT, const std::vector<int>& boundary_triangles, 
                     const Eigen::VectorXi& charts, Eigen::VectorXi& labeling){
    // Initialize queue with all triangles adjacent to boundaries
    // pop and check triangle for unspike operation
    // if triangle changes, add its neighbors to queue

    // possible optim: only add borders of relevant patches to initial queue

    std::queue<int> queue;
    for (const int tri: boundary_triangles) queue.push(tri);

    auto unspikeTriangle = [&labeling, &charts, &TT](int k){
        for (int j=0; j<TT.cols() - 1; j++){
            for (int l=j+1; l<TT.cols(); l++){
                if ((labeling(TT(k,j)) == labeling(TT(k,l)) && labeling(TT(k,j)) != labeling(k))
                    && charts(TT(k,j)) == charts(TT(k,l))){
                    // Found two neighboring triangles with the same label
                    // AND their label is different from my own
                    // So change my label to their's
                    labeling(k) = labeling(TT(k,l));
                    return true;
                } 
            }
        }
        return false;
    };

    int safety_net = 0;
    while (queue.size() > 0){
        int current = queue.front(); 
        queue.pop();

        if (unspikeTriangle(current)){
            queue.push(TT(current, 0));
            queue.push(TT(current, 1));
            queue.push(TT(current, 2));
        }

        if (safety_net > 1000000){
            coloredPrint("ERROR infinite unspiking", "red");
            break;
        }
    }
}

// Improvement: remove need for "old_labeling"?
void growFromTris(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                 const Eigen::MatrixXi& TT, const Eigen::MatrixXd& dists, 
                 int introduced_label, const std::vector<int>& start_tri_ids,
                 double threshold_dist, bool single_chart){
    
    int initial_label = old_labeling(start_tri_ids[0]);
    new_labeling = old_labeling;

    std::unordered_map<int, double> seen;
    std::queue<int> queue;

    for (int id: start_tri_ids){
        queue.push(id);
        seen.insert(std::make_pair(id, 0.0));
    }

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
            if (single_chart && old_labeling(neigh_id) != initial_label) continue;
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

void growFromTri(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                const Eigen::MatrixXi& TT, const Eigen::MatrixXd& dists, 
                int introduced_label, int start_tri_id,
                double threshold_dist){
    growFromTris(old_labeling, new_labeling, TT, dists, introduced_label, {start_tri_id}, threshold_dist);
}


// Given a border and a chart_id, propagate new label on the nearby region on that chart

void growAroundBorder(const Eigen::MatrixXd& dists, 
                       const Eigen::MatrixXi& TT,
                       const std::vector<std::vector<int>>& VT,
                       const Eigen::VectorXi& charts,
                       const std::vector<int>& border,
                       double threshold_dist,
                       int chart_id,
                       int introduced_label,
                       Eigen::VectorXi& labeling){

    std::vector<int> start_tri_ids = {};
    for (int v_id: border){
        for (int tri_id: VT[v_id]){
            if (charts(tri_id) == chart_id){
                start_tri_ids.push_back(tri_id);
            }
        }
    }

    std::sort(start_tri_ids.begin(), start_tri_ids.end());
    start_tri_ids.erase(std::unique(start_tri_ids.begin(), start_tri_ids.end()), start_tri_ids.end());

    Eigen::VectorXi old_labeling = labeling;
    growFromTris(old_labeling, labeling, TT, dists, introduced_label, start_tri_ids, threshold_dist);
}

void fixOppositeLabels(const Eigen::MatrixXd& dists, 
                       const Eigen::MatrixXi& TT,
                       const std::vector<std::vector<int>>& VT,
                       const Eigen::VectorXi& charts,
                       const std::vector<std::vector<int>>& borders,
                       const std::vector<std::pair<int, int>>& patches_per_border,
                       const Eigen::VectorXi& per_chart_labels,
                       double l_avg,
                       const QuickLabelEv& qle,
                       Eigen::VectorXi& labeling){

    // TODO some chart_dirty flag to check it's up to date
    
    for (int b=0; b<patches_per_border.size(); b++){
        int p1 = patches_per_border[b].first;
        int p2 = patches_per_border[b].second;
        if (per_chart_labels(p1) != oppositeLabel(per_chart_labels(p2))){
            continue;
        }
        coloredPrint("Found opposite patches, fixing them...", "cyan");
        const std::vector<int> border = borders[b];

        Eigen::VectorXi best_labeling = labeling;
        double best_score = -1;

        for (int side: {p1,p2}){
            for (int dist_coeff = 1; dist_coeff<10; dist_coeff += 2){
                for (int introduced_label = 0; introduced_label<6; introduced_label++){
                    // OPTIM : we could skip some non eligible labels here
                    if (introduced_label == per_chart_labels(p1) 
                     || introduced_label == per_chart_labels(p2)) continue;

                    Eigen::VectorXi labeling_attempt = labeling;
                    
                    double threshold_dist = dist_coeff * l_avg;
                    int chart_id = side;
                    
                    growAroundBorder(dists, TT, VT, charts, border, threshold_dist, chart_id, introduced_label, labeling_attempt);
                    double score = qle.evaluate(labeling_attempt);
                    if (best_score == -1 || score < best_score){
                        best_score = score;
                        best_labeling = labeling_attempt;
                    }
                }
            }
        }

        labeling = best_labeling;

    }
}

void vertexGrowMutation(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                        const Eigen::VectorXi& charts,
                        const Eigen::MatrixXi& TT, const std::vector<std::vector<int>>& VT,
                        const Eigen::MatrixXd& dists, double threshold_dist,
                        const std::vector<std::vector<int>>& borders,
                        const std::vector<std::pair<int, int>>& patches_per_border,
                        const Eigen::VectorXi& per_chart_labels,
                        int border_id,
                        int vertex_start){

    /*int border_id = std::rand() % borders.size(); // TODO big border big odds

    int n_tps = 0;
    for (auto v: vec_tps) n_tps += v.size();
    bool force_turning_point = true;
    if (force_turning_point && n_tps > 0){
        while (vec_tps[border_id].size() == 0){
            border_id = std::rand() % borders.size();
        }
    }

    int vertex_start = borders[border_id][std::rand() % borders[border_id].size()];

    if (vec_tps[border_id].size() > 0){
        int rand_tp = std::rand() % vec_tps[border_id].size(); 
        vertex_start = borders[border_id][vec_tps[border_id][rand_tp]];
    }*/



    int rand_side = std::rand() % 2;
    
    int chart_pick = rand_side ? patches_per_border[border_id].first : patches_per_border[border_id].second;
    int other_chart = rand_side ? patches_per_border[border_id].second : patches_per_border[border_id].first;
    int introduced_label = per_chart_labels(other_chart);

    std::vector<int> eligible;
    for (int tri: VT[vertex_start]){
        if (charts(tri) == chart_pick) eligible.push_back(tri);
    }
    int start_tri_id = eligible[std::rand() % eligible.size()];

    growFromTri(old_labeling, new_labeling, TT, dists, introduced_label, start_tri_id, threshold_dist);
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