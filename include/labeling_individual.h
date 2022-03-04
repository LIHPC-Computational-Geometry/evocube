#pragma once
#include <memory>
#include <Eigen/Core>

#include "evocube.h"
#include "chart.h"
#include "quick_label_ev.h"
#include "labeling_ops.h"
#include "flagging_utils.h"
#include "graphcut_labeling.h"

class LabelingIndividual {
public:
    LabelingIndividual(std::shared_ptr<Evocube> evo, std::shared_ptr<const QuickLabelEv> qle, Eigen::VectorXi labeling) // first generation only
        : evo_(evo), qle_(qle), labeling_(labeling){ 
        coloredPrint("An individual has come into existence...", "green");
        prev_labeling_ = labeling;
        timestamps_ = Eigen::VectorXi::Zero(labeling.rows());
    }
    LabelingIndividual(const LabelingIndividual& indiv) : evo_(indiv.evo_), qle_(indiv.qle_){
        coloredPrint("An individual rose from the shadows...", "green");
        labeling_ = indiv.getLabeling();
        prev_labeling_ = labeling_;
        timestamps_ = indiv.getTimestamps();
    }

    LabelingIndividual(const LabelingIndividual& parent1, const LabelingIndividual& parent2)
        : evo_(parent1.evo_), qle_(parent1.qle_) { // Cross operator
        coloredPrint("An individual was brought into the world...", "green");
        labeling_ = parent1.getLabeling();
        timestamps_ = parent1.getTimestamps();

        const Eigen::VectorXi labeling2 = parent2.getLabeling();
        const Eigen::VectorXi timestamps2 = parent2.getTimestamps();

        for (int i=0; i<labeling_.rows(); i++){ // OPTIM omp
            if (timestamps2(i) > timestamps_(i)){
                labeling_(i) = labeling2(i);
                timestamps_(i) = timestamps2(i);
            }
        }
        prev_labeling_ = labeling_;

    }

    virtual ~LabelingIndividual(){
        coloredPrint("An individual breathed his last...", "yellow");
    }
    
    Eigen::VectorXi getLabeling() const {return labeling_;}

    void updateChartsAndTPs(bool also_update_turning_points = true){
        // OPTIM if flags are clean, no need to recompute
        int n_charts = updateCharts(evo_->TT_, labeling_, charts_);
        adj = chartAdjacency(evo_->TT_, charts_);
        computeBoundaries(evo_->TT_, evo_->F_, charts_, borders, patches_per_border, boundary_triangles);

        per_chart_labels = perChartLabels(charts_, labeling_);
        per_border_axis = perBorderAxis(patches_per_border, per_chart_labels);
        
        charts_dirty_ = false;

        

        if (also_update_turning_points){
            vec_tps.clear();
            borders_with_tps.clear();
            n_tps = 0;
            for (int i=0; i<borders.size(); i++){
                std::vector<int> tps = graphcutTurningPoints(borders[i], evo_->V_, per_border_axis[i]);
                vec_tps.push_back(tps);
                n_tps += tps.size();
                if (tps.size() > 0) borders_with_tps.push_back(i);
            }
            turning_points_dirty_ = false;
        }
    };

    void updateTimestamps(){
        latest_change_ = evo_->timestamp_;
        for (int i=0; i<labeling_.rows(); i++){ // OPTIM omp
            if (labeling_(i) != prev_labeling_(i)){
                //timestamps_(i) ++;
                timestamps_(i) = evo_->timestamp_;
            }   
        }
        prev_labeling_ = labeling_;
        timestamps_dirty_ = false;
    }

    void repairHighValenceCorner(){
        checkClean(charts_dirty_);
        fixHighValenceCorners(evo_->dists_, evo_->TT_, evo_->VT_, charts_, borders, 
                              evo_->l_avg_, *qle_, labeling_);
        setFlagsDirty();
    }

    void repairOppositeLabels(){
        checkClean(charts_dirty_);
        fixOppositeLabels(evo_->dists_, evo_->TT_, evo_->VT_, charts_, borders, 
                          patches_per_border, per_chart_labels, evo_->l_avg_, *qle_, labeling_);
        setFlagsDirty();
    }

    void repairUnspikeLabeling(){
        unspikeLabeling(evo_->TT_, boundary_triangles, labeling_);
        setFlagsDirty();
        spikes_dirty_ = false;
    }

    void mutationVertexGrow(bool on_turning_point = false){
        checkClean(charts_dirty_);
        int mut_size = 1 + (std::rand() % 7);
        Eigen::VectorXi old_labeling = labeling_;
        double dist_thresh = mut_size * evo_->l_avg_;
        /*vertexGrowMutation(old_labeling, labeling_, charts_, evo_->TT_, evo_->VT_, 
                           evo_->dists_, dist_thresh , borders, 
                           vec_tps, patches_per_border, per_chart_labels);*/

        int b_id = std::rand() % borders.size(); // IMPROVEMENT bigger border bigger odds?
        int vertex_start = borders[b_id][std::rand() % borders[b_id].size()];

        if (on_turning_point){
            checkClean(turning_points_dirty_);
            if (n_tps > 0){
                b_id = borders_with_tps[std::rand() % borders_with_tps.size()];
                vertex_start = borders[b_id][vec_tps[b_id][std::rand() % vec_tps[b_id].size()]];
            }
        }

        vertexGrowMutation(old_labeling, labeling_, charts_, evo_->TT_, evo_->VT_, 
                        evo_->dists_, dist_thresh ,
                        borders, patches_per_border, per_chart_labels, b_id, vertex_start);
        setFlagsDirty();
    }

    void mutationRemoveChart(bool on_invalid_chart = true){
        checkClean(charts_dirty_);
        int chart_id = std::rand() % (charts_.maxCoeff() + 1);
        if (on_invalid_chart){
            std::vector<int> candidates;
            for (int i = 0; i < adj.size(); i++){
                if (adj[i].size() < 4){
                    candidates.push_back(i);
                }
            }
            if (candidates.size() > 0) chart_id = candidates[std::rand() % candidates.size()];
        }
        Eigen::VectorXi old_labeling = labeling_;
        removeChartMutation(old_labeling, labeling_, charts_, evo_->V_, evo_->F_, 
                            evo_->N_, evo_->TT_, chart_id);
        setFlagsDirty();
    }

    void mutationGreedyPath(bool on_turning_point = true){
        checkClean(charts_dirty_);
        int b = std::rand() % borders.size();
        int starting_vertex = borders[b][std::rand() % borders[b].size()];

        if (on_turning_point){
            checkClean(turning_points_dirty_);
            if (n_tps > 0){
                b = borders_with_tps[std::rand() % borders_with_tps.size()];
                starting_vertex = borders[b][vec_tps[b][std::rand() % vec_tps[b].size()]];
            }
        }

        int chart_id = patches_per_border[b].first;
        int other_chart_id = patches_per_border[b].second;
        if (std::rand() % 2) {
            chart_id = patches_per_border[b].second;
            other_chart_id = patches_per_border[b].first;
        }

        int axis1 = per_chart_labels[chart_id]/2;
        int axis2 = per_chart_labels[other_chart_id]/2;

        if (axis1 != axis2){ // only try if the border is valid

            bool follow_opp = std::rand() % 2 ? true : false; // either pick the other axis, or the same one as the opp chart (only the one on the chart itself is forbidden)
            int other_axis = 3 - axis1 - axis2;
            double distance_threshold = 1 * evo_->l_avg_;

            Eigen::RowVector3d direction = Eigen::RowVector3d(0, 0, 0);
            if (!follow_opp) direction(axis2) = 1.0;
            if (follow_opp) direction(other_axis) = 1.0;
            if (std::rand() % 2) direction *= -1;

            int introduced_label;
            if (!follow_opp){
                introduced_label = 2 * other_axis;
                if (std::rand() % 2) introduced_label += 1;
            }
            else {
                introduced_label = per_chart_labels[other_chart_id];
            }
            directionalPathMutation(evo_->VT_, evo_->TT_, evo_->V_, evo_->F_, evo_->dists_, 
                                    borders, patches_per_border, charts_, direction,
                                    starting_vertex, chart_id, distance_threshold, 
                                    introduced_label, labeling_);
            setFlagsDirty();
        }
    }
    
    int invalidChartsScore() const {
        checkClean(charts_dirty_);
        int sum = 0;
        for (int i = 0; i < adj.size(); i++){
            if (adj[i].size() < 4){
                sum += 4 - adj[i].size();
            }
        }
        return sum;
    }

    int invalidBordersScore() const {
        checkClean(charts_dirty_);
        int sum = 0;
        for (int b=0; b<patches_per_border.size(); b++){
            int c1 = patches_per_border[b].first;
            int c2 = patches_per_border[b].second;
            int f1 = per_chart_labels[c1];
            int f2 = per_chart_labels[c2];
            if (f1 == oppositeLabel(f2)) sum ++;
        }
        return sum;
    }

    int invalidCornersScore() const {
        checkClean(charts_dirty_);
        int sum = 0;
        std::vector<int> corner_ids;
        for (auto b: borders){
            corner_ids.push_back(b[0]);
            corner_ids.push_back(b[b.size() - 1]);
        }

        std::vector<int> corner_uni = corner_ids; // corners without duplicates
        std::sort(corner_uni.begin(), corner_uni.end());
        corner_uni.erase(std::unique(corner_uni.begin(), corner_uni.end()), corner_uni.end());

        for (int c: corner_uni){
            if (std::count(corner_ids.begin(), corner_ids.end(), c) > 3){
                sum ++;
            }
        }
        return sum;
    }

    int invalidityScore() const {
        return invalidChartsScore() + invalidBordersScore() + invalidCornersScore();
    }

    // -- viz code -- //

    void getBordersViz(std::vector<Eigen::MatrixXd>& vec_border_begs, 
                       std::vector<Eigen::MatrixXd>& vec_border_ends,
                       std::vector<Eigen::RowVector3d>& vec_border_colors) const {
        checkClean(charts_dirty_);
        for (int b=0; b<borders.size(); b++){
            Eigen::MatrixXd border_begs(borders[b].size() - 1, 3);
            Eigen::MatrixXd border_ends(borders[b].size() - 1, 3);
            for (int i=0; i<borders[b].size() - 1; i++){
                border_begs.row(i) = evo_->V_.row(borders[b][i]);
                border_ends.row(i) = evo_->V_.row(borders[b][i + 1]);
            }
            vec_border_begs.push_back(border_begs);
            vec_border_ends.push_back(border_ends);
            vec_border_colors.push_back(per_border_axis[b]);
        }
    }

    Eigen::MatrixXd getTurningPointsMat() const {
        checkClean(turning_points_dirty_);
        Eigen::MatrixXd out(n_tps, 3);
        int next_id = 0;
        for (int i=0; i<vec_tps.size(); i++){
            for (int j=0; j<vec_tps[i].size(); j++){
                out.row(next_id) = evo_->V_.row(borders[i][vec_tps[i][j]]); 
                next_id ++;
            }
        }
        if (next_id != n_tps) coloredPrint("Turning points viz error, TP not up to date?", "red");
        return out;
    }

    Eigen::VectorXi getTimestamps() const {
        checkClean(timestamps_dirty_);
        return timestamps_;
    }

    const std::shared_ptr<Evocube> evo_;    
    const std::shared_ptr<const QuickLabelEv> qle_;

private:
    Eigen::VectorXi labeling_;

    Eigen::VectorXi prev_labeling_;
    Eigen::VectorXi timestamps_;

    // safety flags 
    bool spikes_dirty_ = true;
    bool charts_dirty_ = true;
    bool turning_points_dirty_ = true;
    bool timestamps_dirty_ = false;
    int latest_change_ = -1;

    // Chart connectivity info
    Eigen::VectorXi charts_;
    std::vector<std::vector<int>> adj; // Chart adjacency graph
    std::vector<std::vector<int>> borders; // ordered borders, represented by consecutive vertex ids
    std::vector<std::vector<int>> vec_tps;
    std::vector<int> borders_with_tps;
    std::vector<Eigen::RowVector3d> per_border_axis; // Axis colinear to border (ideally)
    std::vector<std::pair<int, int>> patches_per_border; // patches on each side of a border
    Eigen::VectorXi per_chart_labels;
    std::vector<int> boundary_triangles; // triangles who have >= 1 edge on a boundary
    int n_tps = 0;

    void setFlagsDirty(){
        // Labeling was changed and chart information is outdated
        spikes_dirty_ = true;
        charts_dirty_ = true;
        turning_points_dirty_ = true;
        timestamps_dirty_ = true;
    }

    void checkClean(bool flag_dirty) const {
        if (flag_dirty) coloredPrint("LabelingIndividual: Flag check failed!", "red");
    }
};
