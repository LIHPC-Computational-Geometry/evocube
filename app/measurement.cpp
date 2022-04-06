#include <igl/readOBJ.h>
#include <igl/doublearea.h>
#include <igl/per_face_normals.h>
#include <iostream>

#include "distortion.h"

int main(int argc, char *argv[]){

    // read input name

    std::string base_name = "../data/DATASET2/OM_smooth/bunny_input_tri";
    if (argc > 1) base_name = argv[1];

    Eigen::MatrixXd V1, V2; // V1 reference triangle mesh, V2 deformed
    Eigen::MatrixXi F;

    igl::readOBJ(base_name + "/boundary.obj", V1, F);
    //igl::readOBJ(base_name + "/boundary.obj", V2, F);
    igl::readOBJ(base_name + "/fast_polycube_surf.obj", V2, F);
    
    Eigen::VectorXd A1, A2;
    igl::doublearea(V1, F, A1);
    igl::doublearea(V2, F, A2);
    A1 /= 2.0;
    A2 /= 2.0;
    double A_m = A1.sum();
    double A_d = A2.sum();

    V2 = V2 * std::sqrt(A_m / A_d);

    igl::doublearea(V1, F, A1);
    igl::doublearea(V2, F, A2);
    A1 /= 2.0;
    A2 /= 2.0;
    A_m = A1.sum();
    A_d = A2.sum();

    Eigen::MatrixXd N, N_def;
    igl::per_face_normals(V1, F, N);
    igl::per_face_normals(V2, F, N_def);

    std::vector<std::pair<double, double>> per_tri_singular_values;
    per_tri_singular_values.resize(F.rows());
    std::vector<Eigen::Matrix2d> jacobians;
    computeJacobians(V1, V2, F, N, N_def, jacobians);

    for (int f_id = 0; f_id < F.rows(); f_id++){
        Eigen::JacobiSVD<Eigen::Matrix2d> svd(jacobians[f_id], Eigen::ComputeThinU | Eigen::ComputeThinV);
        per_tri_singular_values[f_id].first = svd.singularValues()(0);
        per_tri_singular_values[f_id].second = svd.singularValues()(1);
    }

    double stretch = computeStretch(A1, A_m, A_d, per_tri_singular_values);
    double area_disto = computeAreaDisto(A1, per_tri_singular_values);
    double angle_disto =  computeAngleDisto(A1, per_tri_singular_values);

    std::cout << "Area:\t reference " << A_m << " vs deformed " << A_d << std::endl;
    std::cout << "Stretch:\t" << stretch << std::endl;
    std::cout << "AreaDisto:\t" << area_disto << std::endl;
    std::cout << "AngleDisto:\t" << angle_disto << std::endl;

    Eigen::VectorXd disto;
    computeDisto(V1, V2, F, N, N_def, disto);
    std::cout << "Integrated:\t" << integrateDistortion(A1, disto) << std::endl;

    // fill logs

    // measure similarity % between initial graph cut labeling and final labeling?
}