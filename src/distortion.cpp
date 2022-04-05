#include "distortion.h"
#include <omp.h>
#include <igl/doublearea.h>
#include <iostream>

#include "logging.h"

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
        double s1 = svd.singularValues()(0);
        double s2 = svd.singularValues()(1);
        disto(f_id) =  s1 * s2 + 1.0 / (s1 * s2) - 2.0; // area
        disto(f_id) +=  s1 / s2 + s2 / s1 - 2.0; // angle
    }
}

double integrateDistortion(const Eigen::VectorXd& A,
                          const Eigen::VectorXd& disto){

    if (disto.minCoeff() < 0) coloredPrint("ERROR: distortion is < 0", "red");
    
    /*std::cout << "disto range" << std::endl;
    std::cout << disto.minCoeff() << std::endl;
    std::cout << disto.maxCoeff() << std::endl;*/

    int max_disto = 10e8;

    Eigen::VectorXd distop = disto.array().min(static_cast<double>(max_disto));

    return (A.array() * distop.array()).sum() / A.sum(); //+ static_cast<double>(n_inf * max_disto);

}

// Spherical Parametrization and Remeshing
// Praun & Hoppe
double computeStretch(const Eigen::VectorXd& A, double A_m, double A_d,
                      std::vector<std::pair<double, double>> per_tri_singular_values){
    
    double sum = 0;
    for (int f_id = 0; f_id < per_tri_singular_values.size(); f_id++){
        double s1 = per_tri_singular_values[f_id].first;
        double s2 = per_tri_singular_values[f_id].second;
        sum += A(f_id) * (s1 * s1 + s2 * s2) / 2.0;
    }
    sum /= A.sum();

    return (A_m / A_d) * (1.0 / std::pow(sum, 2)); 
}

// PolyCube-Maps
// Tarini & Hormann & Cignoni & Montani
double computeAreaDisto(const Eigen::VectorXd& A, std::vector<std::pair<double, double>> per_tri_singular_values){
    double sum = 0;
    for (int f_id = 0; f_id < per_tri_singular_values.size(); f_id++){
        double s1 = per_tri_singular_values[f_id].first;
        double s2 = per_tri_singular_values[f_id].second;
        sum += A(f_id) * 0.5 * (s1 * s2 + 1.0 / (s1 * s2));
    }

    return sum / A.sum();
}

double computeAngleDisto(const Eigen::VectorXd& A, std::vector<std::pair<double, double>> per_tri_singular_values){
    double sum = 0;
    for (int f_id = 0; f_id < per_tri_singular_values.size(); f_id++){
        double s1 = per_tri_singular_values[f_id].first;
        double s2 = per_tri_singular_values[f_id].second;
        sum += A(f_id) * 0.5 * (s1 / s2 + s2 / s1);
    }

    return sum / A.sum();
}