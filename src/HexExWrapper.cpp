// HexExWrapper.cpp
// https://github.com/fprotais/polycube_withHexEx with Eigen instead of ultimaille

#include <HexEx.hh>
// #include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh> //OpenVolumeMesh::VertexHandle

#include "HexExWrapper.h"

//preserve face orientation
//face order according to README of https://github.com/fprotais/polycube_withHexEx/blob/main/README.md
const Eigen::MatrixXi TETRAHEDRON_TO_FACES = (Eigen::MatrixXi(4,3) << 1, 2, 3,
																	  0, 3, 2,
																	  0, 1, 3,
																	  0, 2, 1).finished();

inline Eigen::RowVector3d geom_swap(const HexEx::Vec3d& v) {
	return Eigen::RowVector3d(v[0], v[1], v[2]);
}

inline HexEx::Vec3d geom_swap(const Eigen::RowVector3d& v) {
	return { v[0], v[1], v[2] };
}

inline void HexEx2Eigen(const HexEx::HexahedralMesh& in, Eigen::MatrixXi& hexes_out, Eigen::MatrixXd& V_hexes_out) {

	//fill V_hexes_out (vertices)
	V_hexes_out.resize(in.n_vertices(),3);
	int ct = 0;
	for (auto it = in.vertices_begin(); it != in.vertices_end(); it++) {
		V_hexes_out.row(ct++) = geom_swap(in.vertex(*it));
	}

	//fill hexes_out (cells)
	hexes_out.resize(in.n_cells(),8);
	// /!\ WARNING : MEDIT convention here (.mesh). Different from UM, OVM conventions
	//      5-------6
	//     /|      /|
	//    / |     / |
	//   1-------2  |
	//   |  4----|--7
	//   | /     | /
	//   |/      |/
	//   0-------3
	constexpr int transition[8] = { 0, 3, 2, 1, 4, 5, 6, 7 };
	for (auto c_it = in.cells_begin(); c_it != in.cells_end(); ++c_it) {
		ct = 0;
		for (auto cv_it = OpenVolumeMesh::HexVertexIter(*c_it, &in); cv_it.valid(); ++cv_it) {
			hexes_out(c_it->idx(), transition[ct++]) = cv_it->idx(); 
		}
	}
}

bool run_HexEx(const Eigen::MatrixXi& tets, const Eigen::MatrixXd& V_tets, const Eigen::MatrixXd& corner_param, Eigen::MatrixXi& hexes, Eigen::MatrixXd& V_hexes) {

	//fill the tetrahedral mesh in the HexEx format
    HexEx::TetrahedralMesh tetMesh;
	std::vector<OpenVolumeMesh::VertexHandle> corners;
	for (int i = 0; i < V_tets.rows(); i++) { //for each vertex of the tetrahedral mesh
		corners.push_back(tetMesh.add_vertex(geom_swap(V_tets.row(i))));
	}

	//fill the parametrization in the HexEx format
	auto parametrization = tetMesh.request_cell_property<std::map<OpenVolumeMesh::VertexHandle, HexEx::Vec3d> >("Parametrization");
	tetMesh.set_persistent(parametrization);
	for(int c = 0; c < tets.rows(); c++) { //for each tet cell
		//get the 4 vertices of c
		int T[4];
		for (int cv = 0; cv < 4; cv++) {
			T[cv] = tets(c,cv); 
		}
		auto chexex = tetMesh.add_cell(corners[T[0]], corners[T[1]], corners[T[2]], corners[T[3]]);
		for (int cv = 0; cv < 4; cv++) {
			parametrization[chexex][corners[T[cv]]] = geom_swap(corner_param.row(4 * c + cv));
		}
	}

	//fill the hexahedral mesh in the HexEx format
	HexEx::HexahedralMesh hexMesh;
	std::cerr << "Running Hexex... ";
	HexEx::extractHexMesh(tetMesh, parametrization, hexMesh);
	std::cerr << "Done!" << std::endl;
	std::cout << "The extracted hex mesh has " << hexMesh.n_cells() << " cells, "
		<< hexMesh.n_faces() << " faces, " << hexMesh.n_edges() << " edges, and "
		<< hexMesh.n_vertices() << " vertices." << std::endl;

	HexEx2Eigen(hexMesh, hexes, V_hexes);
	return true;
}