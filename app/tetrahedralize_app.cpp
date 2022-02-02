
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <iostream>

#include "mesh_io.h"

// see https://github.com/libigl/libigl/blob/main/tutorial/605_Tetgen/main.cpp


// TODO move
#include <fstream>
class BndToTetConverter {
public:
    BndToTetConverter(std::string input_file){
        int n_values;
        int line;
        std::ifstream file(input_file);
        
        if (file.is_open())
        {
            file >> n_values;

            table_.resize(n_values);
            for (int i=0; i < n_values; i++){
                file >> table_(i);
            }
            file.close();
        }

        else std::cout << "Unable to open file";

        std::cout << table_ << std::endl;
    }

    BndToTetConverter(const Eigen::VectorXi table) : table_(table) {}

    void writeTo(std::string write_path){
        std::ofstream file(write_path);
        if (file.is_open()){
            file << table_.rows();
            file << "\n";
            for (int i=0; i<table_.rows(); i++){
                file << table_(i);
                file << "\n";
            }
            file.close();
        }
        else std::cout << "Unable to open file";
    }
private:
    Eigen::VectorXi table_;
};

int main(){

    std::string input_path = "../data/S1/input_boundary.obj";
    std::string output_path = "../data/S1/tetra.mesh";
    std::string output_bnd = "../data/S1/boundary.obj";

    Eigen::MatrixXd V, TV;
    Eigen::MatrixXi F, TF, TT;
    igl::readOBJ(input_path, V, F);

    V.rowwise() -= V.colwise().minCoeff();
    V /= V.maxCoeff() / 10.0;
    std::cout << "min: " << std::endl;
    std::cout << V.colwise().minCoeff() << std::endl;
    std::cout << "max: " << std::endl;
    std::cout << V.colwise().maxCoeff() << std::endl;


    // Mesh the interior of a surface mesh (V,F) using tetgen
    //
    // Inputs:
    //   V  #V by 3 vertex position list
    //   F  #F list of polygon face indices into V (0-indexed)
    //   switches  string of tetgen options (See tetgen documentation) e.g.
    //     "pq1.414a0.01" tries to mesh the interior of a given surface with
    //       quality and area constraints
    //     "" will mesh the convex hull constrained to pass through V (ignores F)
    // Outputs:
    //   TV  #V by 3 vertex position list
    //   TT  #T by 4 list of tet face indices
    //   TF  #F by 3 list of triangle face indices
    // Returns status:
    //   0 success
    //   1 tetgen threw exception
    //   2 tetgen did not crash but could not create any tets (probably there are
    //     holes, duplicate faces etc.)
    //   -1 other error

    // tetgen manual https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual.pdf
    // -p Tetrahedralizes a piecewise linear complex (PLC)
    // -q Refines mesh, see 4.2.3
    // -a Applies a maximum tetrahedron volume constraint.
    // -Y keep boundary mesh from input

    int mesh_refinement = 100;
    double max_tet_volume = 1.0;
    std::string arguments = "pqma" + std::to_string(max_tet_volume);

    int success = igl::copyleft::tetgen::tetrahedralize(V, F, "pYq1000000a10000.0", TV, TT, TF);
    writeDotMeshTet(output_path, TV, TT);
    if (TV.rows() > 70000) std::cout << "Warning! Generated mesh is really large." << std::endl;

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
    igl::boundary_facets(TT, Fb, Fb_to_TT, K);

    BndToTetConverter conv(Fb_to_TT);
    conv.writeTo("../from_tet_to_tris_with_extra_vs.txt");

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
    igl::remove_unreferenced(TV, Fb, NV, NF, I);
    need to code your own remove unref

    igl::writeOBJ(output_bnd, NV, NF);

    return 0;
}