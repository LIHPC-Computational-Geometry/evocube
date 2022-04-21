#include "tet_boundary.h"

#include <iostream>
#include <fstream>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>

#include "logging.h"

BndToTetConverter::BndToTetConverter(std::string input_file){
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

BndToTetConverter::BndToTetConverter(const Eigen::VectorXi table, int n_tets) : table_(table) {
    n_tris_ = table_.rows();
    n_tets_ = n_tets;
}

void BndToTetConverter::writeTo(std::string write_path){
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

void tetVerticesToBoundaryVertices(const Eigen::MatrixXd& Vb, const Eigen::MatrixXi& Fb, 
                                   const Eigen::MatrixXd& V_tets, const Eigen::MatrixXi& tets, 
                                   std::string corres_path,
                                   Eigen::MatrixXd& Vf){

    Vf = Eigen::MatrixXd::Zero(Vb.rows(), Vb.cols());
    BndToTetConverter conv(corres_path);

    std::vector<std::vector<int>> corres = {
        {1, 3, 2},
        {0, 2, 3},
        {0, 3, 1},
        {0, 1, 2}	
    };

    if (Fb.rows() != conv.table_.rows()) coloredPrint("Correspondence mismatch", "red");

    for (int i=0; i<Fb.rows(); i++){
        int f_tet = conv.table_(i);
        int tet_id = f_tet / 4;
        int f_in_tet = f_tet % 4;
        Vf.row(Fb(i, 0)) = V_tets.row(tets(tet_id, corres[f_in_tet][0]));
        Vf.row(Fb(i, 1)) = V_tets.row(tets(tet_id, corres[f_in_tet][1]));
        Vf.row(Fb(i, 2)) = V_tets.row(tets(tet_id, corres[f_in_tet][2]));
    }
}

void tetToBnd(const Eigen::MatrixXd& V_tets, const Eigen::MatrixXi& tets, 
              Eigen::MatrixXd& NV, Eigen::MatrixXi& NF){
    Eigen::VectorXi Fb_to_TT, K;
    Eigen::MatrixXi Fb;
    igl::boundary_facets(tets, Fb, Fb_to_TT, K);
    Eigen::VectorXi I;
    igl::remove_unreferenced(V_tets, Fb, NV, NF, I);
};