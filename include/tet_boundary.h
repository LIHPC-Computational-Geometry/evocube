#pragma once
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>

/**
 * @brief Utility to 
 * 1) recover tet face id from boundary face id 
 * 2) save information for this purpose
 */
class BndToTetConverter {
public:
    BndToTetConverter(std::string input_file){
        int line;
        std::ifstream file(input_file);
        
        if (file.is_open()){
            file >> n_tris_;
            file >> n_tets_;
            table_.resize(n_tris_);
            for (int i=0; i < n_tris_; i++){
                file >> table_(i);
            }
            file.close();
        }

        else std::cout << "BndToTetConverter: Unable to open file";
    }

    BndToTetConverter(const Eigen::VectorXi table, int n_tets) : table_(table) {
        n_tris_ = table_.rows();
        n_tets_ = n_tets;
        
        std::cout << "ok1 " << n_tets << std::endl;
    }

    void writeTo(std::string write_path){
        std::cout << "ok2 " << n_tris_ << std::endl;
        std::cout << "ok3 " << n_tets_ << std::endl;
        std::ofstream file(write_path);
        if (file.is_open()){
            file << n_tris_;
            file << "\n";
            file << n_tets_;
            file << "\n";
            for (int i=0; i<table_.rows(); i++){
                file << table_(i);
                file << "\n";
            }
            file.close();
        }
        else std::cout << "BndToTetConverter: Unable to write to file";
    }
    
    int n_tris_;
    int n_tets_;
    Eigen::VectorXi table_;
};

void computeTetMeshBoundary(const Eigen::MatrixXi& tets, const Eigen::MatrixXd& V,
                            std::string output_tris_to_tets, std::string output_bnd){
    // BOUNDARY_FACETS Determine boundary faces (edges) of tetrahedra (triangles)
    // stored in T (analogous to qptoolbox's `outline` and `boundary_faces`).
    //
    // Input:
    //  T  tetrahedron (triangle) index list, m by 4 (3), where m is the number of tetrahedra
    // Output:
    //  F  list of boundary faces, n by 3 (2), where n is the number of boundary faces
    //  J  list of indices into T, n by 1
    //  K  list of indices revealing across from which vertex is this facet
    //

    Eigen::VectorXi Fb_to_TT, K;
    Eigen::MatrixXi Fb;
    igl::boundary_facets(tets, Fb, Fb_to_TT, K);

    Eigen::VectorXi corres(Fb_to_TT.rows());
    for (int i=0; i<corres.rows(); i++){
        corres(i) = 4 * Fb_to_TT(i) + K(i);
    }

    BndToTetConverter conv(corres, tets.rows());
    conv.writeTo(output_tris_to_tets);

    // Remove unreferenced vertices from V, updating F accordingly
    //
    // Input:
    //   V  #V by dim list of mesh vertex positions
    //   F  #F by ss list of simplices 
    // Outputs:
    //   NV  #NV by dim list of mesh vertex positions
    //   NF  #NF by ss list of simplices
    //   I   #V by 1 list of indices such that: NF = IM(F) and NT = IM(T)
    //      and V(find(IM<=size(NV,1)),:) = NV
    //   J  #NV by 1 list, such that NV = V(J,:)
    Eigen::MatrixXd NV;
    Eigen::MatrixXi NF;
    Eigen::VectorXi I;
    igl::remove_unreferenced(V, Fb, NV, NF, I);

    igl::writeOBJ(output_bnd, NV, NF);
}



