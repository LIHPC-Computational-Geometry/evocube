#pragma once

#include <omp.h>
#include <igl/doublearea.h>
#include <Eigen/Geometry>


Eigen::RowVector3d crossProd(const Eigen::RowVector3d& v1,
                             const Eigen::RowVector3d& v2){
    return v1.transpose().cross(v2.transpose()).transpose();
}

void localTransfo(const Eigen::MatrixXd& V,
                  const Eigen::MatrixXi& F,
                  const Eigen::MatrixXd& N,
                  std::vector<Eigen::Matrix2d>& vecA){

    vecA.resize(F.rows());

#pragma omp parallel for
    for (int f_id = 0; f_id < F.rows(); f_id ++){
        Eigen::RowVector3d local_axis1 = V.row(F(f_id, 1)) - V.row(F(f_id, 0));
        local_axis1.normalize();
        Eigen::RowVector3d local_axis2 = crossProd(N.row(f_id), local_axis1); //N.row(f_id).cross(local_axis1);

        Eigen::RowVector2d local_coords1 = Eigen::RowVector2d((V.row(F(f_id, 1)) - V.row(F(f_id, 0))).norm(), 0);
        Eigen::RowVector2d local_coords2;
        local_coords2(0) = (V.row(F(f_id, 2)) - V.row(F(f_id, 0))).dot(local_axis1);
        local_coords2(1) = (V.row(F(f_id, 2)) - V.row(F(f_id, 0))).dot(local_axis2);

        Eigen::Matrix2d A;
        A.col(0) = local_coords1;
        A.col(1) = local_coords2;
        vecA[f_id] = A;
    }
}

void computeJacobians(const Eigen::MatrixXd& V1,
                      const Eigen::MatrixXd& V2,
                      const Eigen::MatrixXi& F,
                      const Eigen::MatrixXd& N,
                      const Eigen::MatrixXd& N_def,
                      std::vector<Eigen::Matrix2d>& jacobians){
    
    std::vector<Eigen::Matrix2d> vecA1, vecA2;
    localTransfo(V1, F, N, vecA1);
    localTransfo(V2, F, N_def, vecA2);

    jacobians.resize(F.rows());

#pragma omp parallel for
    for (int f_id = 0; f_id < F.rows(); f_id++){
        jacobians[f_id] = vecA2[f_id] * vecA1[f_id].inverse(); 
    }
}

void computeDisto(const Eigen::MatrixXd& V1,
                  const Eigen::MatrixXd& V2,
                  const Eigen::MatrixXi& F,
                  const Eigen::MatrixXd& N,
                  const Eigen::MatrixXd& N_def,
                  Eigen::VectorXd& disto){
    disto.resize(F.rows());

    std::vector<Eigen::Matrix2d> jacobians;
    computeJacobians(V1, V2, F, N, N_def, jacobians);

#pragma omp parallel for
    for (int f_id = 0; f_id < F.rows(); f_id++){
        Eigen::JacobiSVD<Eigen::Matrix2d> svd(jacobians[f_id], Eigen::ComputeThinU | Eigen::ComputeThinV);
        //disto(f_id) = svd.singularValues()(0) + 1.0 / svd.singularValues()(1) - 2.0;
        double s1 = svd.singularValues()(0);
        double s2 = svd.singularValues()(1);
        disto(f_id) =  s1 * s2 + 1.0 / (s1 * s2) - 2.0;
        //disto(f_id) =  s1 / s2 + s2 / s1 - 2.0;
    }
}

double integrateDistortion(const Eigen::MatrixXd& V,
                          const Eigen::MatrixXi& F,
                          const Eigen::VectorXd& disto){
    Eigen::VectorXd A;
    igl::doublearea(V, F, A);

    int n_inf = disto.array().isInf().count();

    std::cout << "infinite values: " << n_inf << std::endl;

    return (A * disto).sum() / A.sum() + static_cast<double>(n_inf * 10);

}