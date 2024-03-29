#include <string>
#include <iostream>
#include <filesystem>
#include <Eigen/Core>
#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/upsample.h>
#include <set>

#include "mesh_io.h"
#include "tet_boundary.h"
#include "logging.h"

#define PREPROCESS_APP ""
//"../../preprocess_polycube/build/preprocess"

int main(){
    enum INPUT_TYPE {TRI_OBJ, TRI_STL, CAD_STEP, TET_MESH};

    // INPUT FOLDER
    std::string input_path = "../../quadqs_models/";
    INPUT_TYPE input_type = CAD_STEP;
    std::string expected_extension = "step";

    // OUTPUT FOLDER
    std::string output_path = "../data/DATASET2/OM_smooth/";

    // COMPUTING OPTIONS
    bool tet_mesh_already_computed = true; // computes tets from surface/CAD 
                                           // + preprocess tet mesh (split some tets and edges, very slow) 
    bool labeling_already_computed = true; // Calls evolabel
    bool hexes_already_computed = true; // Simple hexex-based method, no padding, no smoothing

    // DEBUG OPTIONS
    bool break_after_first = false;
    bool skip_first = false;
    bool skip_if_folder_exists = false;
    bool skip_already_valid = false;
    bool skip_logs_exist = false;
    bool remesh_invalid = false;

    if (remesh_invalid){
        skip_already_valid = true;
        tet_mesh_already_computed = false;
    }

    //input_path = "../data/DATASET/OM_cad_meshes/";
    // "../data/2019-OctreeMeshing/input/smooth/"


    /*
    input_path = "../data/Inputs/smooth/";
    input_type = TRI_OBJ;
    expected_extension = "obj";
    output_path = "../data/DATASET2/OM_smooth/";
    //*/

    /*
    input_path = "../data/Inputs/octreemeshes_cad_christophed/";
    input_type = TET_MESH;
    expected_extension = "mesh";
    output_path = "../data/DATASET2/OM_cad/";
    //*/

    /*
    input_path = "../../quadqs_models/";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/abc/";
    //*/

    /*
    input_path = "../data/Inputs/mambo/Basic";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/basic_mambo/";
    //*/

    //*
    input_path = "../data/Inputs/mambo/Simple";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/simple_mambo/";
    //*/

    /*
    input_path = "../data/Inputs/mambo/Medium";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/medium_mambo/";
    //*/
    
    /*
    input_path = "../data/custom/";
    input_type = TRI_OBJ;
    expected_extension = "obj";
    output_path = "../data/custom_out/";
    //*/
    
    // -------------------- INIT_FROM_FOLDER --------------------

    std::string path_to_preprocess_app = PREPROCESS_APP;

    // Output names (you probably shouldn't change these)
    std::string tet_output = "/tetra.mesh";
    std::string step_input = "/input.step";
    std::string obj_input = "/input.obj";
    std::string stl_input = "/input.stl";
    std::string output_tris_to_tets = "/tris_to_tets.txt";
    std::string output_labeling_on_tets = "/labeling_on_tets.txt";
    std::string output_bnd = "/boundary.obj";
    std::string output_hex = "/hexes.mesh";
    std::string logs_file = "/logs.json";
    std::string polycube_bnd_float = "/fast_polycube_surf.obj";//fast floating-point surface of the polycube
    std::string polycube_bnd_int = "/polycube_surf_int.obj";//fast integer surface of the polycube

    std::set<std::filesystem::path> entries; // set orders entries
    for (auto &entry : std::filesystem::directory_iterator(input_path)){
        entries.insert(entry.path());
    }

    for (const auto & entry : entries){
        if (skip_first){
            skip_first = false;
            continue;
        }
        std::cout << "Dealing with: " << entry << std::endl;
        std::string base_filename = std::string(entry).substr(std::string(entry).find_last_of("/\\") + 1);
        std::string::size_type const p(base_filename.find_first_of('.'));
        std::string file_without_extension = base_filename.substr(0, p);
        std::string extension = base_filename.substr(base_filename.find_last_of('.')+1, base_filename.size());
        if (extension != expected_extension) continue;

        std::string new_folder = output_path + file_without_extension;
        std::string output_mesh = new_folder + tet_output;
        std::string boundary_obj_path = new_folder + output_bnd;
        std::string logs_path = new_folder + logs_file;

        if (skip_if_folder_exists && std::filesystem::is_directory(new_folder)){
            std::cout << "Folder already exists, skipping " << new_folder << std::endl;
            continue;
        }

        if (skip_logs_exist && readLogsValue("LabelingFinal", "InvalidTotal", logs_path) != "null"){
            std::cout << "Logs already exist, skipping " << new_folder << std::endl;
            continue;
        }

        if (skip_already_valid && readLogsValue("LabelingFinal", "InvalidTotal", logs_path) == "0"){
            std::cout << "Labeling already valid, skipping " << new_folder << std::endl;
            continue;
        }

        // ---- COMPUTE TET MESH + SPLIT DOUBLE BOUNDARY TETS ---- //
        if (!tet_mesh_already_computed){
            if (std::filesystem::is_directory(new_folder)){
                std::filesystem::remove_all(new_folder);
            }
            std::filesystem::create_directory(new_folder);

            Eigen::MatrixXd TV(0,0);
            Eigen::MatrixXi TT(0,0);

            switch (input_type){
            case TRI_OBJ: {
                std::string input_obj = new_folder + obj_input;
                std::filesystem::copy(entry, input_obj);
                triObjtoMeshTet(input_obj, output_mesh, TV, TT);
                break;
            }
            case TRI_STL:{ // This one is not recommended 
                std::string input_stl = new_folder + stl_input;
                std::filesystem::copy(entry, input_stl);

                bool refine_mesh = true;
                if (refine_mesh){
                    Eigen::MatrixXd V, N, NV;
                    Eigen::MatrixXi F, NF; 
                    FILE * fp = fopen(input_stl.c_str(), "rb");
                    igl::readSTL(fp, V, F, N);
                    std::cout << F.rows() << std::endl;
                    input_stl = new_folder + "/refined_input.stl";

                    NV = V;
                    NF = F;
                    while (NV.rows() < 3000){
                        coloredPrint("NV.rows() ", "blue");
                        std::cout << NV.rows() << std::endl;
                        V = NV;
                        F = NF;
                        igl::upsample(V, F, NV, NF);
                        break;
                    }
                    igl::writeSTL(input_stl, NV, NF);
                }

                std::string command = "/usr/bin/python3 ../scripts/stl_to_tet.py " + input_stl + " " + output_mesh;
                int result = system(command.c_str());
                if (result) {
                    coloredPrint("FAILURE on " + input_stl, "red");
                    break;
                }
                readDotMeshTet(output_mesh, TV, TT);
                break;
            }
            case CAD_STEP:{
                std::string input_step = new_folder + step_input;
                std::filesystem::copy(entry, input_step);
                std::string command = "/usr/bin/python3 ../scripts/step_to_tet.py " + input_step + " " + output_mesh;
                int result = system(command.c_str());
                if (result) {
                    coloredPrint("FAILURE on " + input_step, "red");
                    continue;
                }
                readDotMeshTet(output_mesh, TV, TT);
                break;
            }
            case TET_MESH:{
                std::filesystem::copy(entry, output_mesh);

                std::string file_without_mesh = std::string(entry) // name with potentially several "."
                                            .substr(std::string(entry).find_last_of("/\\") + 1,
                                                    std::string(entry).find_last_of('.') + 1);
                file_without_mesh = file_without_mesh.substr(0, file_without_mesh.size() - 5);
                
                std::string png_screenshot = std::string(entry)
                                            .substr(0, std::string(entry).find_last_of("/\\") + 1)
                                            + "../octreemeshes_cad_screenshots/" 
                                            + file_without_mesh
                                            + ".png";
                //std::filesystem::copy(png_screenshot, new_folder + "/screenshot.png");
                readDotMeshTet(output_mesh, TV, TT);
                break;
            }
            default:
                coloredPrint("Case not handled", "red");
                break;
            }

            if (path_to_preprocess_app != ""){
                std::string command = path_to_preprocess_app + " " + output_mesh + " " + output_mesh;
                int result = system(command.c_str());
                if (result) {
                    coloredPrint("FAILURE (preprocess) on " + output_mesh, "red");
                    std::cout << "Can you make sure that variable path_to_preprocess_app is properly set?" << std::endl;
                    continue;
                }
                readDotMeshTet(output_mesh, TV, TT);
            }
            
            computeTetMeshBoundary(TT, TV, new_folder + output_tris_to_tets, boundary_obj_path);
        } // !tet_mesh_already_computed
        
        // ---- COMPUTE LABELING ---- //
        if (!labeling_already_computed){
            resetLogFile(logs_path);
            std::string command_labeling = "./evolabel " + boundary_obj_path + " " + new_folder;
            int result_labeling = system(command_labeling.c_str());
            if (result_labeling == 2) return 2;
        }

        // ---- COMPUTE HEX MESH (SIMPLE METHOD) ---- //
        if (!hexes_already_computed){
            std::string scale = "1.3";
            std::string command_hexex = "./polycube_withHexEx " + new_folder + " " + scale;
            int result_hexex = system(command_hexex.c_str());
            if (result_hexex == 2) return 2;
        }

        // ---- POST-PROCESSING MEASUREMENTS ---- //
        //compute polycube distortion. 2 inputs : boundary.obj and polycube surface. Will update the JSON in the same folder as the polycube
        std::cout << "Measures for the fast float polycube surf" << std::endl;
        std::string command_measurement_1 = "./measurement " + boundary_obj_path + " " + new_folder + polycube_bnd_float;
        int result_measurement_1 = system(command_measurement_1.c_str());
        if (result_measurement_1 == 2) return 2;

        std::cout << "Measures for the fast int polycube surf" << std::endl;
        std::string command_measurement_2 = "./measurement " + boundary_obj_path + " " + new_folder + polycube_bnd_int;
        int result_measurement_2 = system(command_measurement_2.c_str());
        if (result_measurement_2 == 2) return 2;

        std::cout << "Measures for the final polycube:" << std::endl;
        std::string command_measurement_3 = "./measurement " + boundary_obj_path + " " + new_folder + "/polycube_final.obj";
        int result_measurement_3 = system(command_measurement_3.c_str());
        if (result_measurement_3 == 2) return 2;

        // ---- GENERATE FIGURES ---- //
        std::string command_figure = "./figure_generator " + file_without_extension + " 0 " + output_path;
        command_figure += " && ./figure_generator " + file_without_extension + " 1 " + output_path;
        command_figure += " && ./figure_generator " + file_without_extension + " 4 " + output_path;
        int result_figures = system(command_figure.c_str());
        
        if (break_after_first){
            std::cout << "BREAKING" << std::endl;
            break;
        }
    }
        
}
