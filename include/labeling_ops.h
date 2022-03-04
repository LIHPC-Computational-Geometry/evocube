#pragma once

#include <Eigen/Core>
#include "quick_label_ev.h"

void fixHighValenceCorners(const Eigen::MatrixXd& dists, 
                           const Eigen::MatrixXi& TT, 
                           const std::vector<std::vector<int>>& VT, 
                           const Eigen::VectorXi& charts,
                           const std::vector<std::vector<int>>& borders,
                           double l_avg,
                           const QuickLabelEv& qle,
                           Eigen::VectorXi& labeling);

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
                             Eigen::VectorXi& labeling);

void unspikeLabeling(const Eigen::MatrixXi& TT, const std::vector<int>& boundary_triangles, Eigen::VectorXi& labeling);
 
void growFromTris(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                 const Eigen::MatrixXi& TT, const Eigen::MatrixXd& dists, 
                 int introduced_label, const std::vector<int>& start_tri_ids,
                 double threshold_dist, bool single_chart = true);

// Remove this one ?
/**
 * @brief Same as previously but with a single starting triangle
 * 
 * @param old_labeling 
 * @param new_labeling 
 * @param TT 
 * @param dists 
 * @param introduced_label 
 * @param start_tri_id 
 * @param threshold_dist 
 */
void growFromTri(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                const Eigen::MatrixXi& TT, const Eigen::MatrixXd& dists, 
                int introduced_label, int start_tri_id,
                double threshold_dist);

void growAroundBorder(const Eigen::MatrixXd& dists, 
                       const Eigen::MatrixXi& TT,
                       const std::vector<std::vector<int>>& VT,
                       const Eigen::VectorXi& charts,
                       const std::vector<int>& border,
                       double threshold_dist,
                       int chart_id,
                       int introduced_label,
                       Eigen::VectorXi& labeling);

void fixOppositeLabels(const Eigen::MatrixXd& dists, 
                       const Eigen::MatrixXi& TT,
                       const std::vector<std::vector<int>>& VT,
                       const Eigen::VectorXi& charts,
                       const std::vector<std::vector<int>>& borders,
                       const std::vector<std::pair<int, int>>& patches_per_border,
                       const Eigen::VectorXi& per_chart_labels,
                       double l_avg,
                       const QuickLabelEv& qle,
                       Eigen::VectorXi& labeling);

void vertexGrowMutation(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                        const Eigen::VectorXi& charts,
                        const Eigen::MatrixXi& TT, const std::vector<std::vector<int>>& VT,
                        const Eigen::MatrixXd& dists, double threshold_dist,
                        const std::vector<std::vector<int>>& borders,
                        const std::vector<std::pair<int, int>>& patches_per_border,
                        const Eigen::VectorXi& per_chart_labels,
                        int border_id,
                        int vertex_start);

void removeChartMutation(const Eigen::VectorXi& old_labeling, Eigen::VectorXi& new_labeling, 
                         const Eigen::VectorXi& charts, 
                         const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, 
                         const Eigen::MatrixXd& N, 
                         const Eigen::MatrixXi& TT, int chart_id);