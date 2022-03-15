#include <string>
#include <iostream>
#include <filesystem>
#include <Eigen/Core>
#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/upsample.h>

#include "mesh_io.h"
#include "tet_boundary.h"
#include "logging.h"

int main(){
    enum INPUT_TYPE {TRI_OBJ, TRI_STL, CAD_STEP, TET_MESH};

    // INPUT FOLDER
    std::string input_path = "../../quadqs_models/";
    INPUT_TYPE input_type = CAD_STEP;
    std::string expected_extension = "step";

    // OUTPUT FOLDER
    std::string output_path = "../data/DATASET2/OM_smooth/";

    //input_path = "../data/DATASET/OM_cad_meshes/";
    // "../data/2019-OctreeMeshing/input/smooth/"


    /*input_path = "../data/2019-OctreeMeshing/input/smooth/";
    input_type = TRI_OBJ;
    expected_extension = "obj";
    output_path = "../data/DATASET2/OM_smooth/";*/

    /*input_path = "../../quadqs_models/";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/abc/";*/

    /*
    input_path = "../../mambo/Basic";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/basic_mambo/";
    //*/

    /*
    input_path = "../../mambo/Simple";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/simple_mambo/";
    //*/

    input_path = "../../mambo/Medium";
    input_type = CAD_STEP;
    expected_extension = "step";
    output_path = "../data/DATASET2/medium_mambo/";

    

    std::string path_to_preprocess_app = "../../preprocess_polycube/build/preprocess";

    bool tet_mesh_already_computed = true;


    // Output names (you probably shouldn't change these)
    std::string tet_output = "/tetra.mesh";
    std::string step_input = "/input.step";
    std::string obj_input = "/input.obj";
    std::string stl_input = "/input.stl";
    std::string output_tris_to_tets = "/tris_to_tets.txt";
    std::string output_labeling_on_tets = "/labeling_on_tets.txt";
    std::string output_bnd = "/boundary.obj";
    std::string output_hex = "/hexes.mesh";

    bool skip_first = true;

    for (const auto & entry : std::filesystem::directory_iterator(input_path)){
        if (skip_first) {
            skip_first = false;
            continue;
        }
        std::cout << "Dealing with: " << entry << std::endl;
        std::string base_filename = std::string(entry.path()).substr(std::string(entry.path()).find_last_of("/\\") + 1);
        std::string::size_type const p(base_filename.find_first_of('.'));
        std::string file_without_extension = base_filename.substr(0, p);
        std::string extension = base_filename.substr(base_filename.find_last_of('.')+1, base_filename.size());
        if (extension != expected_extension) continue;

        std::string new_folder = output_path + file_without_extension;
        std::string output_mesh = new_folder + tet_output;
        std::string boundary_obj_path = new_folder + output_bnd;

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
                std::filesystem::copy(entry.path(), input_obj);
                triObjtoMeshTet(input_obj, output_mesh, TV, TT);
                break;
            }
            case TRI_STL:{ // This one is not recommended 
                std::string input_stl = new_folder + stl_input;
                std::filesystem::copy(entry.path(), input_stl);

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
                    continue;
                }
                readDotMeshTet(output_mesh, TV, TT);
                break;
            }
            case CAD_STEP:{
                std::string input_step = new_folder + step_input;
                std::filesystem::copy(entry.path(), input_step);
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
                std::filesystem::copy(entry.path(), output_mesh);

                std::string file_without_mesh = std::string(entry.path()) // name with potentially several "."
                                            .substr(std::string(entry.path()).find_last_of("/\\") + 1,
                                                    std::string(entry.path()).find_last_of('.') + 1);
                file_without_mesh = file_without_mesh.substr(0, file_without_mesh.size() - 5);
                std::cout << "gotcha: " << file_without_mesh << std::endl;
                
                std::string png_screenshot = std::string(entry.path())
                                            .substr(0, std::string(entry.path()).find_last_of("/\\") + 1)
                                            + "../OM_cad_screenshots/" 
                                            + file_without_mesh
                                            + ".png";
                std::filesystem::copy(png_screenshot, new_folder + "/screenshot.png");
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

        std::string command_labeling = "./test_greedy " + boundary_obj_path + " " + new_folder;
        int result_labeling = system(command_labeling.c_str());

        std::string scale = "1.3";
        std::string command_hexex = "./polycube_withHexEx " + output_mesh + " " + new_folder + output_labeling_on_tets + " " + new_folder + output_hex + " " + scale;
        int result_hexex = system(command_hexex.c_str());

        //./test_greedy ../data/DATASET2/OM_smoothscrewdriver_input_tri/boundary.obj 0
        //./polycube_withHexEx ../data/DATASET2/medium_mambo/M6/tetra.mesh ../data/DATASET2/medium_mambo/M6/labeling_on_tets.txt ../data/DATASET2/medium_mambo/M6/hexes.mesh 1.4

        //break;
    }
        
}