#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>

// .mesh documentation : https://www.ljll.math.upmc.fr/frey/logiciels/Docmedit.dir/index.html
// (spelling mistake in HexaHedra)

void writeDotMeshTet(std::string filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& tets);

void writeDotMeshTet(std::string filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& tets){
    // tets : matrix of size (n_tets, 4)

    std::ofstream out_mesh;
    out_mesh.open (filename);

    out_mesh << "MeshVersionFormatted 1\n";
    out_mesh << "Dimension\n3\n";

    out_mesh << "Vertices\n";
    out_mesh << V.rows() << "\n";
    for (int i=0; i<V.rows(); i++){
        out_mesh << V(i,0) << " " << V(i, 1) << " " << V(i,2) << " 0\n";
    }

    out_mesh << "Tetrahedra\n";
    out_mesh << tets.rows() << "\n";
    for (int i=0; i<tets.rows(); i++){
        for (int j=0; j<4; j++){
            // Warning : .mesh indices start at 1, not 0
            out_mesh << tets(i,j) + 1 << " ";
        }
        out_mesh << " 0\n";
    }

    out_mesh << "End\n";
    out_mesh.close();
}
