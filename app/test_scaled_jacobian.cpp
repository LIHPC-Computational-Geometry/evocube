#include <Eigen/Core>
#include <string>
#include <iostream>
#include <cassert>
#include <filesystem>
#include <fstream>

#include "scaled_jacobian.h"
#include "logging.h"
#include "mesh_io.h"

int main(int argc, char* argv[]) {

	if (argc < 3) {
		std::cout << "Usage is: " << argv[0] << " path/to/hexes.mesh path/to/SJ.csv" << std::endl;
		return 1;
	}

    std::string hex_mesh_file = argv[1], csv_file = argv[2];

    if(!std::filesystem::exists(hex_mesh_file)) {
        coloredPrint((hex_mesh_file) + " does not exist","red");
        return 1;
    }

    Eigen::MatrixXi hexes;
	Eigen::MatrixXd V_hexes;
    Eigen::VectorXd min_sj;
    double overall_min_sj;
	
	readDotMeshHex(hex_mesh_file, V_hexes, hexes);

    overall_min_sj = compute_min_scaled_jacobian(hexes,V_hexes,min_sj);

    std::cout << "Overall minSJ = " << overall_min_sj << std::endl;

    //export as CSV
    std::ofstream ofs(csv_file.c_str(),std::ofstream::out);
    ofs << "hex index,minSJ" << std::endl;
    for(int h = 0; h < min_sj.rows(); h++) {
        ofs << h << "," << min_sj(h) << std::endl;
    }
    ofs.close();

    std::cout << "table write into " << csv_file << std::endl;

    return 0;
}