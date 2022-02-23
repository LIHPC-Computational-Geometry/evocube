#pragma once
#include <Eigen/Core>


int updateCharts(const Eigen::MatrixXi& TT, const Eigen::VectorXi& labeling, Eigen::VectorXi& charts);

Eigen::VectorXi perChartLabels(const Eigen::VectorXi& charts, const Eigen::VectorXi& labeling);

std::vector<Eigen::RowVector3d> perBorderAxis(const std::vector<std::pair<int, int>>& patches_per_border,
                                              const Eigen::VectorXi& per_chart_labels);

std::vector<std::vector<int>> chartAdjacency(const Eigen::MatrixXi& TT, const Eigen::VectorXi& charts);

void computeBoundaries(const Eigen::MatrixXi& TT, const Eigen::MatrixXi& F, 
                       const Eigen::VectorXi& charts,
                       std::vector<std::vector<int>>& ordered_borders,
                       std::vector<std::pair<int,int>>& patches_per_border,
                       std::vector<int>& border_triangles);