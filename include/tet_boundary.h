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

// "inverse" operation: From tet mesh, compute boundary, matching the order of (Fb, Vb) 
// and the correspondences defined in corres_path 
void tetVerticesToBoundaryVertices(const Eigen::MatrixXd& Vb, const Eigen::MatrixXi& Fb, 
                                   const Eigen::MatrixXd& V_tets, const Eigen::MatrixXi& tets, 
                                   std::string corres_path,
                                   Eigen::MatrixXd& Vf);

// TODO avoid duplicate with computeTetMeshBoundary
void tetToBnd(const Eigen::MatrixXd& V_tets, const Eigen::MatrixXi& tets, 
              Eigen::MatrixXd& NV, Eigen::MatrixXi& NF);