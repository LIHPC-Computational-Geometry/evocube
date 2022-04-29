#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <igl/readOBJ.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

#include "logging.h"

// .mesh documentation : https://www.ljll.math.upmc.fr/frey/logiciels/Docmedit.dir/index.html
// (spelling mistake in HexaHedra)

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
        coloredPrint("Error : expected Vertices", "red");
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
        coloredPrint("Error : expected Tetrahedra", "red");
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

void triObjtoMeshTet(std::string input_path, std::string output_path, 
                     Eigen::MatrixXd& TV, Eigen::MatrixXi& TT){
    Eigen::MatrixXd V;
    Eigen::MatrixXi F, TF;
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

    int success = igl::copyleft::tetgen::tetrahedralize(V, F, "pqa", TV, TT, TF);
    writeDotMeshTet(output_path, TV, TT);
    if (TV.rows() > 70000) std::cout << "Warning! Generated mesh is really large." << std::endl;
}

void triObjtoMeshTet(std::string input_path, std::string output_path){
    Eigen::MatrixXd TV;
    Eigen::MatrixXi TT;
    triObjtoMeshTet(input_path, output_path, TV, TT);
}

void readDotMeshTri(std::string input_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F){
    std::ifstream input_file(input_path);

    if (input_file.fail()){
        std::cout << "\033[31mWarning\033[0m: no triangles\n" << std::endl;
        V = Eigen::MatrixXd(0, 3);
        F = Eigen::MatrixXi(0, 3);
    }

    // TODO vertices part is a copy past from readTet function, factorize ?
    std::string header;
    for (int i=0; i<4; i++){
        input_file >> header;
    }

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
    while (section_name != "Triangles" && input_file >> section_name){
        discarded ++;
    }

    std::cout << discarded << " elements discarded" << std::endl;

    if (section_name != "Triangles"){
        std::cout << "Error : expected Triangles" << std::endl;
    }
    else {
        int n_tri;
        input_file >> n_tri;
        std::cout << n_tri << " tri" << std::endl;
        F.resize(n_tri, 3);

        for (int i=0; i<n_tri; i++){
            int a, b, c;
            input_file >> a;
            input_file >> b;
            input_file >> c;
            F.row(i) = Eigen::RowVector3i(a-1, b-1, c-1);
            input_file >> a; // discard ref element
        }
    }
    input_file.close();
}

void writeDotMeshHex(std::string path, Eigen::MatrixXd V, Eigen::MatrixXi hexahedra) {
    // hexahedra : matrix of size (n_hexa, 8)

    std::cout << "Saving hex mesh: " << path << std::endl;

    std::ofstream out_mesh;
    out_mesh.open(path);

    out_mesh << "MeshVersionFormatted 1\n";
    out_mesh << "Dimension\n3\n";

    out_mesh << "Vertices\n";
    out_mesh << V.rows() << "\n";
    for (int i=0; i<V.rows(); i++){
        out_mesh << V(i,0) << " " << V(i, 1) << " " << V(i,2) << " 0\n";
    }

    out_mesh << "Hexahedra\n";
    out_mesh << hexahedra.rows() << "\n";
    
    for (int i=0; i<hexahedra.rows(); i++){
        for (int j=0; j<8; j++){
            // Warning : .mesh indices start at 1, not 0
            out_mesh << hexahedra(i,j) + 1 << " ";
        }
        out_mesh << " 0\n";
    }

    out_mesh << "End\n";
    out_mesh.close();
}

void readDotMeshHex(std::string input_hex, Eigen::MatrixXd &V, Eigen::MatrixXi &hexes) {
    // Note: INDICES START AT 1 !!

    std::ifstream input_file(input_hex);

    std::string header;
    for (int i=0; i<4; i++){
        input_file >> header;
        //std::cout << header <<" ";
    }
    //std::cout << std::endl;

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
    while (section_name != "Hexahedra" && input_file >> section_name){
        discarded ++;
    }

    std::cout << discarded << " elements discarded" << std::endl;

    if (section_name != "Hexahedra"){
        std::cout << "Error : expected Hexahedra" << std::endl;
    }
    else {
        int n_hex;
        input_file >> n_hex;
        std::cout << n_hex << " hexes" << std::endl;
        hexes.resize(n_hex, 8);

        for (int i=0; i<n_hex; i++){
            int a,b,c,d,e,f,g,h;
            input_file >> a;
            input_file >> b;
            input_file >> c;
            input_file >> d;
            input_file >> e;
            input_file >> f;
            input_file >> g;
            input_file >> h;
            Eigen::RowVectorXi row(8);
            row << a-1, b-1, c-1, d-1, e-1, f-1, g-1, h-1;
            hexes.row(i) = row;
            input_file >> a; // discard ref element
        }
    }
    input_file.close();
}

void readDotMeshEdges(std::string input_edges, Eigen::MatrixXd &V, Eigen::MatrixXi &edges) {
    std::ifstream input_file(input_edges);

    if (input_file.fail()){
        std::cout << "\033[31mWarning\033[0m: no feature edges.\n" << std::endl;
        V = Eigen::MatrixXd(0, 3);
        edges = Eigen::MatrixXi(0, 2);
    }

    // TODO vertices part is a copy past from readTet function, factorize ?
    std::string header;
    for (int i=0; i<4; i++){
        input_file >> header;
    }

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
    while (section_name != "Edges" && input_file >> section_name){
        discarded ++;
    }

    std::cout << discarded << " elements discarded" << std::endl;

    if (section_name != "Edges"){
        std::cout << "Error : expected Edges" << std::endl;
    }
    else {
        int n_edges;
        input_file >> n_edges;
        std::cout << n_edges << " edges" << std::endl;
        edges.resize(n_edges, 2);

        for (int i=0; i<n_edges; i++){
            int a,b;
            input_file >> a;
            input_file >> b;
            edges.row(i) = Eigen::RowVector2i(a-1, b-1);
            input_file >> a; // discard ref element
        }
    }
    input_file.close();
}