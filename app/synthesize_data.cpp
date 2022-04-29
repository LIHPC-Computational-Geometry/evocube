
#include <Eigen/Core>
#include <igl/readOBJ.h>
#include <iostream>
#include <fstream>
#include <set>
#include <filesystem>

#include "labeling_individual.h"
#include "mesh_io.h"
#include "scaled_jacobian.h"
#include "tet_boundary.h"
#include "distortion.h"

enum DATA_FORMAT {POLYCUT_OUTPUT, OM_HEXES, OURS_HEXES, OURS_POLYCUBES, OURS_ABC};
std::vector<std::string> polycut_tags = {"red", "darkred", "green", "darkgreen", "blue", "darkblue"};

// TODO put somewhere else

double surfacePolycubeMeasure(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F,
                                    Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2){
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

    Eigen::VectorXd A1, A2;
    igl::doublearea(V1, F, A1);
    igl::doublearea(V2, F, A2);
    A1 /= 2.0;
    A2 /= 2.0;
    double A_m = A1.sum();
    double A_d = A2.sum();

    V2 = V2 * std::sqrt(A_m / A_d); // TODO do this before saving the polycube instead?

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

    std::cout << "Polycube stats: " << std::endl;
    std::cout << "\tStretch:\t" << stretch << std::endl;
    std::cout << "\tAreaDisto:\t" << area_disto << std::endl;
    std::cout << "\tAngleDisto:\t" << angle_disto << std::endl;
    std::cout << "\tIsometricDisto:\t" << isometric_disto << std::endl;

    return area_disto;
}

double volumetricPolycubeMeasure(std::string base_name){
    Eigen::MatrixXi F, F1, F2;
    Eigen::MatrixXd V1, V2, V_tets1, V_tets2; // V1 reference triangle mesh, V2 deformed

    // For tet meshes comparison:
    Eigen::MatrixXi tets1, tets2;
    readDotMeshTet(base_name + "_remesh.mesh", V_tets1, tets1);
    readDotMeshTet(base_name + "_defo.mesh", V_tets2, tets2);
    tetToBnd(V_tets1, tets1, V1, F1);
    tetToBnd(V_tets2, tets2, V2, F2);
    F = F1;

    return surfacePolycubeMeasure(V1, F, V2, F2);
}

int main(int argc, char *argv[]){
    DATA_FORMAT data_format = OURS_ABC;
    std::string folder_path = argv[1];

    std::set<std::string> entries; // set orders entries
    for (auto &entry : std::filesystem::directory_iterator(folder_path)){
        std::string str = std::string(entry.path());
        if (data_format == OURS_HEXES && str.substr(str.size() - 15, str.size()) != "_outputhex.mesh") continue;
        if (data_format == OURS_POLYCUBES && str.substr(str.size() - 12, str.size()) != "_remesh.mesh") continue;
        entries.insert(str);
    }

    Eigen::VectorXi labeling_success = Eigen::VectorXi::Zero(entries.size());
    Eigen::VectorXi hex_success = Eigen::VectorXi::Zero(entries.size());
    Eigen::VectorXd min_sj_vec = Eigen::VectorXd::Constant(entries.size(), -1.0);
    Eigen::VectorXd avg_sj_vec = Eigen::VectorXd::Constant(entries.size(), -1.0);
    Eigen::VectorXd poly_area_disto_vec = Eigen::VectorXd::Constant(entries.size(), 10000000.0);
    Eigen::VectorXi poly_success = Eigen::VectorXi::Zero(entries.size());

    auto fillHexInfo = [&hex_success, &min_sj_vec, &avg_sj_vec](std::string hex_path, int input_id){
        Eigen::MatrixXi hexes;
        Eigen::MatrixXd V_hexes;
        Eigen::VectorXd min_sj;
        double overall_min_sj;
        
        readDotMeshHex(hex_path, V_hexes, hexes);
        if (hexes.rows() == 0) return;
        overall_min_sj = compute_min_scaled_jacobian(hexes, V_hexes, min_sj);

        if (overall_min_sj >= -10e-5) hex_success(input_id) = 1;

        std::cout << "minSJ: " << overall_min_sj << std::endl;
        std::cout << "avgSJ: " << min_sj.mean() << std::endl;
        min_sj_vec(input_id) = overall_min_sj;
        avg_sj_vec(input_id) = min_sj.mean();
    };

    int input_id = -1;
    for (std::string folder: entries){
        std::cout << "Input: " << folder << std::endl;
        input_id ++;

        if (data_format == POLYCUT_OUTPUT){
            std::string obj_path = folder + "/segmentation_XYZ.obj";
            if (!std::ifstream(obj_path).good()){
                std::cout << "NO LABELING FILE: " << folder << std::endl;
                continue;
            }

            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            igl::readOBJ(obj_path, V, F);

            Eigen::VectorXi labeling = Eigen::VectorXi::Constant(F.rows(), -1);

            int next_id = 0;
            std::fstream file;
            std::string word;
            file.open(obj_path.c_str());
            while (file >> word) {
                if (word == "usemtl"){
                    file >> word;
                    for (int i=0; i<6; i++){
                        if (word == polycut_tags[i]){
                            labeling(next_id) = i;
                            next_id++;
                            break;
                        }
                    }
                }
            }
            file.close();

            std::shared_ptr<Evocube> evo = std::make_shared<Evocube>(Evocube(V, F));
            std::shared_ptr<const QuickLabelEv> qle = std::make_shared<const QuickLabelEv>(QuickLabelEv(V, F));
            LabelingIndividual indiv(evo, qle, labeling);


            std::cout << indiv.invalidityScore() << std::endl;
            if (indiv.invalidityScore() == 0) labeling_success(input_id) = 1;
            else std::cout << "INVALID LABELING: " << folder << std::endl;

            std::string hex_path = folder + "/hexa_presmooth_highest_quality.mesh";
            if (!std::ifstream(hex_path).good()){
                std::cout << "No hex output " << folder << std::endl;
            }
            else {
                fillHexInfo(hex_path, input_id);
            }
        }

        else if (data_format == OM_HEXES){
            // labeling_success is N/A
            // for the rest, use fillHexInfo(hex_path, input_id);
        }
        else if (data_format == OURS_HEXES){
            // need to find hexes_highest_quality.mesh first, like polycut

            std::vector<std::string> candidates = {"_outputhex.mesh", 
                                                   //"_outputhex_coarse01.mesh", too coarse to make sense
                                                   "_outputhex_coarse05.mesh",
                                                   "_outputhex_fine.mesh",
                                                   "_result.mesh", 
                                                   //"_result_coarse01.mesh",
                                                   "_result_coarse05.mesh",
                                                   "_result_fine.mesh"};

            if (folder.substr(folder.size() - 15, folder.size()) != "_outputhex.mesh") continue;

            std::string model_name = folder.substr(0, folder.size() - 15);
            std::string model = model_name.substr(model_name.find_last_of("/\\") + 1,
                                 model_name.size());

            double best_sj = -1;
            std::string best_candidate = "";

            for (std::string name: candidates){
                Eigen::MatrixXi hexes;
                Eigen::MatrixXd V_hexes;
                Eigen::VectorXd min_sj;
                double overall_min_sj;
                readDotMeshHex(model_name + name, V_hexes, hexes);
                overall_min_sj = compute_min_scaled_jacobian(hexes, V_hexes, min_sj);
                std::cout << overall_min_sj << std::endl;

                if (overall_min_sj > best_sj) {
                    best_sj = overall_min_sj;
                    best_candidate = name;
                }
            }

            std::cout << "Best " << best_sj << " from " << best_candidate << std::endl;
            //std::filesystem::copy(model_name + best_candidate, "../best/" + model + "_best.mesh");
            fillHexInfo(model_name + best_candidate, input_id);
        }
        else if (data_format == OURS_POLYCUBES){
            if (folder.substr(folder.size() - 12, folder.size()) != "_remesh.mesh") continue;

            std::string model_name = folder.substr(0, folder.size() - 12);
            std::string model = model_name.substr(model_name.find_last_of("/\\") + 1,
                                 model_name.size());

            poly_area_disto_vec(input_id) = volumetricPolycubeMeasure(model_name);
            if (poly_area_disto_vec(input_id) < 2.0) poly_success(input_id) = 1;
        }
        else if (data_format == OURS_ABC) {
            std::string hex_path = folder + "/hexes.mesh";
            if (!std::ifstream(hex_path).good()){
                std::cout << "No hex output " << folder << std::endl;
            }
            else {
                fillHexInfo(hex_path, input_id);
            }

            Eigen::MatrixXi F1, F2;
            Eigen::MatrixXd V1, V2; // V1 reference triangle mesh, V2 deformed
            igl::readOBJ(folder + "/boundary.obj", V1, F1);
            igl::readOBJ(folder + "/fast_polycube_surf.obj", V2, F2);
            

            poly_area_disto_vec(input_id) = surfacePolycubeMeasure(V1, F1, V2, F2);;
            if (poly_area_disto_vec(input_id) < 2.0) poly_success(input_id) = 1;
        }
    }

    std::cout << "Labeling success: " << labeling_success.sum() << " / " << labeling_success.rows() << std::endl;
    std::cout << "hex success: " << hex_success.sum() << " / " << hex_success.rows() << std::endl;
    std::cout << "avg minSJ: " << min_sj_vec.mean() << std::endl;
    std::cout << "avg avgSJ: " << avg_sj_vec.mean() << std::endl;
    //std::cout << "poly_area_disto_vec: " << poly_area_disto_vec << std::endl;
    std::cout << "poly success: " << poly_success.sum() << " / " << poly_success.rows() << std::endl;
    
}