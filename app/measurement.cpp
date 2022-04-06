#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/doublearea.h>
#include <igl/per_face_normals.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <iostream>

#include "distortion.h"
#include "mesh_io.h"
#include "logging.h"
#include "tet_boundary.h"

void tetToBnd(const Eigen::MatrixXd& V_tets, const Eigen::MatrixXi& tets, 
              Eigen::MatrixXd& NV, Eigen::MatrixXi& NF){
    Eigen::VectorXi Fb_to_TT, K;
    Eigen::MatrixXi Fb;
    igl::boundary_facets(tets, Fb, Fb_to_TT, K);
    Eigen::VectorXi I;
    igl::remove_unreferenced(V_tets, Fb, NV, NF, I);
}

int main(int argc, char *argv[]){

    std::string base_name = "../data/DATASET2/OM_smooth/bunny_input_tri";
    if (argc > 1) base_name = argv[1];

    // TODO cleanup and have two modes: one for tri meshes, and one for tet meshes

    Eigen::MatrixXi F, F1, F2;
    Eigen::MatrixXd Vb, V1, V2, V_tets1, V_tets2; // V1 reference triangle mesh, V2 deformed

    igl::readOBJ(base_name + "/boundary.obj", Vb, F);
    V1 = Vb;
    igl::readOBJ(base_name + "/fast_polycube_surf.obj", V2, F2);

    // For tet meshes comparison:
    Eigen::MatrixXi tets1, tets2;
    /*
    readDotMeshTet(base_name + "/vol.mesh", V_tets1, tets1);
    readDotMeshTet(base_name + "/poly.mesh", V_tets2, tets2);
    tetToBnd(V_tets1, tets1, V1, F1);
    tetToBnd(V_tets2, tets2, V2, F2);
    F = F1;
    //*/

    
    ///*
    // Checking input validity
    if (V1.rows() != V2.rows()) coloredPrint("ERROR: vertex sizes don't match", "red");
    if (F.rows() != F2.rows()) coloredPrint("ERROR: face sizes don't match", "red");
    else {
        for (int i=0; i<F.rows(); i+= F.rows()/15){
            if (F(i,1) != F2(i,1)){
                coloredPrint("ERROR: Face mismatch", "red");
                break;
            }
        }
    }
    //*/

    //igl::writeOBJ("../debug1.obj", V1, F);
    //igl::writeOBJ("../debug2.obj", V2, F2);

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
    double isometric_disto =  computeIsometricDisto(A1, per_tri_singular_values);

    std::cout << "Area:\t reference " << A_m << " vs deformed " << A_d << std::endl;
    std::cout << "Stretch:\t" << stretch << std::endl;
    std::cout << "AreaDisto:\t" << area_disto << std::endl;
    std::cout << "AngleDisto:\t" << angle_disto << std::endl;
    std::cout << "IsometricDisto:\t" << isometric_disto << std::endl;

    Eigen::VectorXd disto;
    computeDisto(V1, V2, F, N, N_def, disto);
    std::cout << "Integrated:\t" << integrateDistortion(A1, disto) << std::endl;

    // TODO fill logs

    // measure similarity % between initial graph cut labeling and final labeling?
}