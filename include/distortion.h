#pragma once

#include <Eigen/Geometry>

Eigen::RowVector3d crossProd(const Eigen::RowVector3d& v1,
                             const Eigen::RowVector3d& v2);

void localTransfo(const Eigen::MatrixXd& V,
                  const Eigen::MatrixXi& F,
                  const Eigen::MatrixXd& N,
                  std::vector<Eigen::Matrix2d>& vecA);

void computeJacobians(const Eigen::MatrixXd& V1,
                      const Eigen::MatrixXd& V2,
                      const Eigen::MatrixXi& F,
                      const Eigen::MatrixXd& N,
                      const Eigen::MatrixXd& N_def,
                      std::vector<Eigen::Matrix2d>& jacobians);

void computeDisto(const Eigen::MatrixXd& V1,
                  const Eigen::MatrixXd& V2,
                  const Eigen::MatrixXi& F,
                  const Eigen::MatrixXd& N,
                  const Eigen::MatrixXd& N_def,
                  Eigen::VectorXd& disto);

double integrateDistortion(const Eigen::VectorXd& A,
                          const Eigen::VectorXd& disto);

// --- Stand-alone functions for polycube metrics --- // 

// Spherical Parametrization and Remeshing
// Praun & Hoppe
double computeStretch(const Eigen::VectorXd& A, double A_m, double A_d,
                      std::vector<std::pair<double, double>> per_tri_singular_values);

// PolyCube-Maps
// Tarini & Hormann & Cignoni & Montani
double computeAreaDisto(const Eigen::VectorXd& A, std::vector<std::pair<double, double>> per_tri_singular_values);
double computeAngleDisto(const Eigen::VectorXd& A, std::vector<std::pair<double, double>> per_tri_singular_values);
