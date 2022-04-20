
#include <Eigen/Core>
#include <igl/readOBJ.h>
#include <iostream>
#include <fstream>
#include <set>
#include <filesystem>

#include "labeling_individual.h"
#include "mesh_io.h"
#include "scaled_jacobian.h"

enum DATA_FORMAT {POLYCUT_OUTPUT, OM_HEXES, OURS_HEXES, GRAPHCUT_INIT};
std::vector<std::string> polycut_tags = {"red", "darkred", "green", "darkgreen", "blue", "darkblue"};

int main(int argc, char *argv[]){
    DATA_FORMAT data_format = POLYCUT_OUTPUT;
    std::string folder_path = argv[1];

    std::set<std::string> entries; // set orders entries
    for (auto &entry : std::filesystem::directory_iterator(folder_path)){
        entries.insert(std::string(entry.path()));
    }

    Eigen::VectorXi labeling_success = Eigen::VectorXi::Zero(entries.size());
    Eigen::VectorXi hex_success = Eigen::VectorXi::Zero(entries.size());
    Eigen::VectorXd min_sj_vec = Eigen::VectorXd::Constant(entries.size(), -1.0);
    Eigen::VectorXd avg_sj_vec = Eigen::VectorXd::Constant(entries.size(), -1.0);

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

            }
            else {
                Eigen::MatrixXi hexes;
                Eigen::MatrixXd V_hexes;
                Eigen::VectorXd min_sj;
                double overall_min_sj;
                
                readDotMeshHex(hex_path, V_hexes, hexes);
                overall_min_sj = compute_min_scaled_jacobian(hexes, V_hexes, min_sj);

                if (overall_min_sj >= -10e-5) hex_success(input_id) = 1;

                std::cout << "minSJ: " << overall_min_sj << std::endl;
                std::cout << "avgSJ: " << min_sj.mean() << std::endl;
                min_sj_vec(input_id) = overall_min_sj;
                avg_sj_vec(input_id) = min_sj.mean();
            }
        }

        else if (data_format == OM_HEXES){

        }
        else if (data_format == OURS_HEXES){

        } 
        else if (data_format == GRAPHCUT_INIT){

        }
    }

    std::cout << "Labeling success: " << labeling_success.sum() << " / " << labeling_success.rows() << std::endl;
    std::cout << "hex success: " << hex_success.sum() << " / " << hex_success.rows() << std::endl;
    std::cout << "avg minSJ: " << min_sj_vec.mean() << std::endl;
    std::cout << "avg avgSJ: " << avg_sj_vec.mean() << std::endl;
}