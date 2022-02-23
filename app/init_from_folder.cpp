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

    // INPUT FOLDER
    std::string input_type = "obj";
    std::string input_path = "../data/2019-OctreeMeshing/input/smooth/";//"../../mambo/Simple/";
    input_type = "step";
    input_path = "../../mambo/Basic/";

    // OUTPUT FOLDER
    std::string output_path = "../data/mambo/";

    // Output names (you probably shouldn't change these)
    std::string tet_output = "/tetra.mesh";
    std::string step_input = "/input.step";
    std::string obj_input = "/input.obj";
    std::string stl_input = "/input.stl";
    std::string output_tris_to_tets = "/tris_to_tets.txt";
    std::string output_bnd = "/boundary.obj";


    for (const auto & entry : std::filesystem::directory_iterator(input_path)){
        std::string base_filename = std::string(entry.path()).substr(std::string(entry.path()).find_last_of("/\\") + 1);
        std::string::size_type const p(base_filename.find_last_of('.'));
        std::string file_without_extension = base_filename.substr(0, p);
        std::string extension = base_filename.substr(p+1, base_filename.size());
        if (extension != input_type) continue;

        std::string new_folder = output_path + file_without_extension;
        std::string output_mesh = new_folder + tet_output;

        if (std::filesystem::is_directory(new_folder)){
            std::filesystem::remove_all(new_folder);
        }
        std::filesystem::create_directory(new_folder);


        Eigen::MatrixXd TV(0,0);
        Eigen::MatrixXi TT(0,0);

        if (input_type == "step"){
            std::string input_step = new_folder + step_input;
            std::filesystem::copy(entry.path(), input_step);
            std::string command = "/usr/bin/python3 ../scripts/step_to_tet.py " + input_step + " " + output_mesh;
            int result = system(command.c_str());
            if (result) {
                coloredPrint("FAILURE on " + input_step, "red");
                continue;
            }
            readDotMeshTet(output_mesh, TV, TT);
        }

        if (input_type == "obj"){
            std::string input_obj = new_folder + obj_input;
            std::filesystem::copy(entry.path(), input_obj);
            triObjtoMeshTet(input_obj, output_mesh, TV, TT);
        }

        if (input_type == "stl"){
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
        }
        
        computeTetMeshBoundary(TT, TV, new_folder + output_tris_to_tets, new_folder + output_bnd);

        //break;
    }
        
}