#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/doublearea.h>
#include <igl/per_face_normals.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <map>

#include "distortion.h"
#include "mesh_io.h"
#include "logging.h"
#include "flagging_utils.h"
#include "tet_boundary.h"

void tetToBnd(const Eigen::MatrixXd& V_tets, const Eigen::MatrixXi& tets, 
              Eigen::MatrixXd& NV, Eigen::MatrixXi& NF){
    Eigen::VectorXi Fb_to_TT, K;
    Eigen::MatrixXi Fb;
    igl::boundary_facets(tets, Fb, Fb_to_TT, K);
    Eigen::VectorXi I;
    igl::remove_unreferenced(V_tets, Fb, NV, NF, I);
}

const std::map<std::string,std::string> polycube_filename_to_JSON_tag {
    {"fast_polycube_surf.obj", "FastPolycubeFloat"},
    {"polycube_surf_int.obj",  "FastPolycubeInt"}
};

#define DEFAULT_TRI_INPUT   "../data/DATASET2/OM_smooth/bunny_input_tri/boundary.obj"
#define DEFAULT_POLYCUBE    "../data/DATASET2/OM_smooth/bunny_input_tri/fast_polycube_surf.obj"
#define LABELING_GRAPHCUT_FILENAME  "/labeling_init.txt" //expected to be in the same folder as polycube_filepath
#define LABELING_FINAL_FILENAME     "/labeling.txt"      //expected to be in the same folder as polycube_filepath



int main(int argc, char *argv[]){

    // read input name

    std::string input_filepath = DEFAULT_TRI_INPUT;
    if (argc > 1) input_filepath = argv[1];

    std::string polycube_filepath = DEFAULT_POLYCUBE, polycube_filename;
    if (argc > 2) polycube_filepath = argv[2];
    polycube_filename = std::filesystem::path(polycube_filepath).filename();

    std::string JSON_tag;
    try {
        JSON_tag = polycube_filename_to_JSON_tag.at(polycube_filename);
    }
    catch (const std::out_of_range&) {
        JSON_tag = polycube_filename;
    }

    std::string save_path = std::filesystem::path(polycube_filepath).parent_path().string();//by default, write the logs inside the same folder as the polycube file
    if (argc > 3) save_path = argv[3];

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

    igl::readOBJ(input_filepath, V1, F);
    igl::readOBJ(polycube_filepath, V2, F);
    
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

    // fill logs
    std::string logs_path = save_path + "/logs.json";
    std::stringstream stretch_rounded, area_disto_rounded, angle_disto_rounded;
    fillLogInfo(JSON_tag, "ReferenceArea", logs_path, A_m);
    fillLogInfo(JSON_tag, "DeformedArea", logs_path, A_d);
    //reduce the precision : the values will be printed in the supplemental and too many digits is superfluous
    stretch_rounded << std::fixed << std::setprecision(3) << stretch;
    area_disto_rounded << std::fixed << std::setprecision(3) << area_disto;
    angle_disto_rounded << std::fixed << std::setprecision(3) << angle_disto;
    fillLogInfo(JSON_tag, "Stretch", logs_path, stretch_rounded.str());
    fillLogInfo(JSON_tag, "AreaDistortion", logs_path, area_disto_rounded.str());
    fillLogInfo(JSON_tag, "AngleDistortion", logs_path, angle_disto_rounded.str());

    //compute the similarity between the graphcut and the final labeling
    std::string labeling_graphcut_filepath  = save_path + LABELING_GRAPHCUT_FILENAME,
                labeling_final_filepath     = save_path + LABELING_FINAL_FILENAME;
    Eigen::VectorXi labeling_graphcut   = openFlagging(labeling_graphcut_filepath,F.rows()),
                    labeling_final      = openFlagging(labeling_final_filepath,F.rows());
    if(labeling_graphcut.isZero() || labeling_final.isZero()) {
        //case one of the labeling file don't have the expected size
        coloredPrint("Error : invalid number of labels in one of the labeling files","red");
    }
    else {
        double similarity = flaggingSimilarity(labeling_graphcut,labeling_final);
        if(similarity >= 0.0) {
            std::cout << "LabelingSimilarity:\t" << similarity << std::endl;
            std::stringstream similarity_rounded;
            similarity_rounded << std::fixed << std::setprecision(3) << similarity;
            fillLogInfo("LabelingSimilarity", logs_path, similarity_rounded.str());
        }
        else {
            coloredPrint("Error : flaggingSimilarity() says invalid vector lengths","red");
        }
    }

    std::cout << logs_path << " updated" << std::endl;
}
