#pragma once

#include <Eigen/Core>

Eigen::VectorXi graphcutFlagging(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& N, const Eigen::MatrixXi& TT,
        const Eigen::VectorXi& locked_flags, const Eigen::VectorXi& forbidden_flags,
        int compact_coeff, int fidelity_coeff);

Eigen::VectorXi graphcutFlagging(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                                 const Eigen::MatrixXd& N, const Eigen::MatrixXi& TT,
                                 int compact_coeff, int fidelity_coeff);

std::vector<int> graphcutTurningPoints(const std::vector<int>& bnd, const Eigen::MatrixXd& V,
                                       const Eigen::RowVector3d& desired_dir);