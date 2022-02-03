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

void readDotMeshTet(std::string input_tets, Eigen::MatrixXd &V, Eigen::MatrixXi &tets);
void readDotMeshTet(std::string input_tets, Eigen::MatrixXd &V, Eigen::MatrixXi &tets){
    
    // Note: INDICES START AT 1 !!

    std::ifstream input_file(input_tets);

    std::string header;
    for (int i=0; i<4; i++){
        input_file >> header;
        std::cout << header <<" ";
    }
    std::cout << std::endl;

    std::string section_name;
    input_file >> section_name;

    if (section_name != "Vertices"){
        std::cout << "Error : expected Vertices" << std::endl;
    }

    int n_vertices;
    input_file >> n_vertices;
    std::cout << n_vertices << " vertices" << std::endl;
    V.resize(n_vertices, 3);
    for (int i=0; i<n_vertices; i++){
        double a,b,c;
        input_file >> a;
        input_file >> b;
        input_file >> c;
        V.row(i) = Eigen::RowVector3d(a, b, c);
        input_file >> a; // discard ref element
    }

    int discarded = 0;
    input_file >> section_name;
    while (section_name != "Tetrahedra" && input_file >> section_name){
        discarded ++;
    }

    std::cout << discarded << " elements discarded" << std::endl;

    if (section_name != "Tetrahedra"){
        std::cout << "Error : expected Tetrahedra" << std::endl;
    }
    else {
        int n_tets;
        input_file >> n_tets;
        std::cout << n_tets << " tets" << std::endl;
        tets.resize(n_tets, 4);

        for (int i=0; i<n_tets; i++){
            int a,b,c,d;
            input_file >> a;
            input_file >> b;
            input_file >> c;
            input_file >> d;
            tets.row(i) = Eigen::RowVector4i(a-1, b-1, c-1, d-1);
            input_file >> a; // discard ref element
        }
    }
    input_file.close();
}
