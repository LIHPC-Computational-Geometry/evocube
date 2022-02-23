#pragma once
#include <Eigen/Core>
#include <string>

class BndToTetConverter {
public:
    BndToTetConverter(std::string input_file);
    BndToTetConverter(const Eigen::VectorXi table, int n_tets);
    void writeTo(std::string write_path);
    
    int n_tris_;
    int n_tets_;
    Eigen::VectorXi table_;
};

void computeTetMeshBoundary(const Eigen::MatrixXi& tets, const Eigen::MatrixXd& V,
                            std::string output_tris_to_tets, std::string output_bnd);