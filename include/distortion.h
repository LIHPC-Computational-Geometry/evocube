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