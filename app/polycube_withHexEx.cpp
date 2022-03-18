// polycube_withHexEx.cpp
// https://github.com/fprotais/polycube_withHexEx with Eigen instead of ultimaille

#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <OpenNL_psm.h>
#include <cassert>
#include <filesystem>
#include <Eigen/LU>

#include "HexExWrapper.h"
#include "mesh_io.h"
#include "tet_boundary.h"
#include "disjointset.h"
#include "logging.h"
#include "tet_boundary.h"

#define DEFAULT_SCALE			1.0

#define CHECK_IF_FILE_EXISTS(filepath) if(!std::filesystem::exists(filepath)) { coloredPrint((filepath) + " does not exist","red"); return 1; }

void float_cubecover(const Eigen::MatrixXi& tets, const Eigen::MatrixXd& V, const Eigen::VectorXi& tets_labeling, Eigen::MatrixXd& U) {
	U=V;
	std::cerr << "Integrating with float boundary...\n";
	DisjointSet ds(V.rows() * 3);
	for(int c = 0; c < tets.rows(); c++) { //for each cell
		for(int cf = 0; cf < 4; cf++) { //for each face
			if(tets_labeling(4 * c + cf) != -1) {
				int d = tets_labeling(4 * c + cf) / 2;//get direction
				for(int cfv = 0; cfv < 3; cfv++) { //for each vertex
					ds.merge(d * V.rows() + TETRAHEDRON_TO_VERTEX(tets, c, cf, cfv), d * V.rows() + TETRAHEDRON_TO_VERTEX(tets, c, cf, (cfv + 1) % 3));
				}
			}
		}
	}
	std::vector<int> idmap;
	int nb_var = ds.get_sets_id(idmap);

	auto context = nlNewContext();
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(nb_var));

	nlBegin(NL_SYSTEM);
	nlEnable(NL_VERBOSE);
	nlBegin(NL_MATRIX);
	for(int c = 0; c < tets.rows(); c++) { //for each cell (tetrahedron)
		int v[4] = { tets(c,0) , tets(c,1), tets(c,2), tets(c,3) };//get the 4 vertices of c
		Eigen::Matrix3d M;
		M.row(0) = V.row(v[1]) - V.row(v[0]);
		M.row(1) = V.row(v[2]) - V.row(v[0]);
		M.row(2) = V.row(v[3]) - V.row(v[0]);
		Eigen::Matrix3d invM = M.inverse();
		invM.transposeInPlace();
		Eigen::MatrixXd grad_coef(4,3);
		grad_coef.row(0) = -invM.row(0) - invM.row(1) - invM.row(2);
		grad_coef.row(1) = invM.row(0);
		grad_coef.row(2) = invM.row(1);
		grad_coef.row(3) = invM.row(2);
		for(int dim = 0; dim < 3; dim++) {
			for(int dim2 = 0; dim2 < 3; dim2++) {
				Eigen::Vector3i e(0,0,0);
				e(dim2) = 1;
				nlBegin(NL_ROW);
				for(int dim_e = 0; dim_e < 3; dim_e++) {
					for(int point = 0; point < 4; point++) {
						nlCoefficient(idmap[dim * V.rows() + v[point]], e(dim_e) * grad_coef(point,dim_e));
					}
				}
				if (dim == dim2) nlRightHandSide(1);
				nlEnd(NL_ROW);
			}
		}
	}

	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);
	nlSolve();

	for(int v = 0; v < V.rows(); v++) {
		for(int dim = 0; dim < 3; dim++) {
			U(v,dim) = nlGetVariable(idmap[V.rows() * dim + v]);
		}
	}

	nlDeleteContext(context);
	std::cerr << " Done.\n";

}

void integer_cubecover(const Eigen::MatrixXi& tets, const Eigen::MatrixXd& V, const Eigen::VectorXi& tets_labeling, const Eigen::MatrixXd& U, Eigen::MatrixXd& int_U) {
	
	int_U.resize(V.rows(),3);
	
	std::cerr << "Integrating with int boundary...\n";
	DisjointSet ds(V.rows() * 3);
	for(int c = 0; c < tets.rows(); c++) { //for each cell
		for(int cf = 0; cf < 4; cf++) { //for each face
			if(tets_labeling(4 * c + cf) != -1) {
				int d = tets_labeling(4 * c + cf) / 2;//get direction
				for(int cfv = 0; cfv < 3; cfv++) { //for each vertex
					ds.merge(d * V.rows() + TETRAHEDRON_TO_VERTEX(tets, c, cf, cfv), d * V.rows() + TETRAHEDRON_TO_VERTEX(tets, c, cf, (cfv + 1) % 3));
				}
			}
		}
	}
	std::vector<int> idmap;
	int nb_var = ds.get_sets_id(idmap);

	auto context = nlNewContext();
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(nb_var));

	nlBegin(NL_SYSTEM);
	for(int c = 0; c < tets.rows(); c++) { //for each cell
		for(int cf = 0; cf < 4; cf++) { //for each face
			if(tets_labeling(4 * c + cf) != -1) {
				int d = tets_labeling(4 * c + cf) / 2;//get direction
				for(int cfv = 0; cfv < 3; cfv++) { //for each vertex
					nlSetVariable(idmap[d * V.rows() + TETRAHEDRON_TO_VERTEX(tets, c, cf, cfv)], std::round(U(TETRAHEDRON_TO_VERTEX(tets, c, cf, cfv),d)));
					nlLockVariable(idmap[d * V.rows() + TETRAHEDRON_TO_VERTEX(tets, c, cf, cfv)]);
				}
			}
		}
	}
	nlEnable(NL_VERBOSE);
	nlBegin(NL_MATRIX);
	for(int c = 0; c < tets.rows(); c++) { //for each cell (tetrahedron)
		int v[4] = { tets(c,0) , tets(c,1), tets(c,2), tets(c,3) };//get the 4 vertices of c
		Eigen::Matrix3d M;
		M.row(0) = V.row(v[1]) - V.row(v[0]);
		M.row(1) = V.row(v[2]) - V.row(v[0]);
		M.row(2) = V.row(v[3]) - V.row(v[0]);
		Eigen::Matrix3d invM = M.inverse();
		invM.transposeInPlace();

		Eigen::MatrixXd grad_coef(4,3);
		grad_coef.row(0) = -invM.row(0) - invM.row(1) - invM.row(2);
		grad_coef.row(1) = invM.row(0);
		grad_coef.row(2) = invM.row(1);
		grad_coef.row(3) = invM.row(2);
		for(int dim = 0; dim < 3; dim++) {
			for(int dim2 = 0; dim2 < 3; dim2++) {
				Eigen::Vector3i e(0,0,0);
				e(dim2) = 1;
				nlBegin(NL_ROW);
				for(int dim_e = 0; dim_e < 3; dim_e++) {
					for(int point = 0; point < 4; point++) {
						nlCoefficient(idmap[dim * V.rows() + v[point]], e(dim_e) * grad_coef(point,dim_e));
					}
				}
				if (dim == dim2) nlRightHandSide(1);
				nlEnd(NL_ROW);
			}
		}
	}

	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);
	nlSolve();

	for(int v = 0; v < V.rows(); v++) {
		for(int dim = 0; dim < 3; dim++) {
			int_U(v,dim) = nlGetVariable(idmap[V.rows() * dim + v]);
		}
	}

	nlDeleteContext(context);
	std::cerr << " Done.\n";
}

double rescaling(const Eigen::MatrixXi& tets, Eigen::MatrixXd& V, double scale) {
	double size = 0.0;
	for(int c = 0; c < tets.rows(); c++) { //for each tetrahedron
		for(int cf = 0; cf < 4; cf++) { //for each face of the current tetrahedron
			for(int cfv = 0; cfv < 3; cfv++) { //for each vertex of the current face
				//distance between vertex cfv and the next one
				Eigen::RowVector3d v1 = V.row(TETRAHEDRON_TO_VERTEX(tets,c,cf,cfv));
				Eigen::RowVector3d v2 = V.row(TETRAHEDRON_TO_VERTEX(tets,c,cf,(cfv + 1) % 3));
				size += (v1 - v2).norm();
			}
		}
	}
	size /= tets.rows() * 12 * scale;
	for(int vertex = 0; vertex < V.rows(); vertex++) {
		V.row(vertex) /= size;
	}
	return size;
}

void revert_rescaling(Eigen::MatrixXd& V, double sizing) {
	for(int vertex = 0; vertex < V.rows(); vertex++) {
		V.row(vertex) *= sizing;
	}
}

int main(int argc, char* argv[]) {

	std::string folder, tetra_file, tets_labeling_file, hex_mesh_file;
	double hexscale = DEFAULT_SCALE;

	if (argc == 3){
		folder = argv[1];
		tetra_file = folder + "/tetra.mesh";
		tets_labeling_file = folder + "/labeling_on_tets.txt";
		hex_mesh_file = folder + "/hexes.mesh";
		hexscale = std::stod(argv[2]);
	}
	else if (argc == 5){
		folder = "";
		tetra_file = argv[1];
		tets_labeling_file = argv[2];
		hex_mesh_file = argv[3];
		hexscale = std::stod(argv[4]);
	}
	else {
		std::cout << "Usage is: " << argv[0] << " tetra.mesh labeling_on_tets.txt output_hexes.mesh scale" << std::endl;
		std::cout << "exemple: " << argv[0] << " ../data/S1/tetra.mesh ../data/S1/labeling_on_tets.txt ../data/S1/hexes.mesh 1." << std::endl;
		return 1;
	}

	CHECK_IF_FILE_EXISTS(tetra_file)
	CHECK_IF_FILE_EXISTS(tets_labeling_file)

	Eigen::MatrixXd V_tets, V_boundary;
    Eigen::MatrixXi tets, F_boundary;
    Eigen::VectorXi tri_labeling, tets_labeling;
	
	//read tets
	readDotMeshTet(tetra_file, V_tets, tets);

	//read labels (flags)
	tets_labeling = Eigen::VectorXi::Constant(4 * tets.rows(), -1);
	std::ifstream ifs(tets_labeling_file);
	if (!ifs.is_open()) {
		std::cerr << "Failed opening of flags at : " << tets_labeling_file << std::endl;
		abort();
	}
	for(int tet_face = 0; tet_face < 4 * tets.rows(); tet_face++) { //for each face of each tetrahedron
		if (ifs.eof()) tets_labeling(tet_face) = -1;
		else ifs >> tets_labeling(tet_face);
	}
	ifs.close();

	Eigen::MatrixXd U, int_U;
	double sizing = rescaling(tets, V_tets, hexscale);

	// writeDotMeshTet("tetra_scaled.mesh",V_tets,tets);

	float_cubecover(tets, V_tets, tets_labeling, U);

	integer_cubecover(tets, V_tets, tets_labeling, U, int_U);

	revert_rescaling(V_tets, sizing);


	Eigen::MatrixXd corner_param(tets.rows()*4,3);
	for(int c=0; c < tets.rows(); c++) { //for each cell (tetrahedron)
		for(int cc = 0; cc < 4; cc++) { //for each vertex of the cell c
			corner_param.row(4 * c + cc) = int_U.row(tets(c,cc));
		}
	}

	Eigen::MatrixXi hexes;
	Eigen::MatrixXd V_hexes;
	run_HexEx(tets, V_tets, corner_param, hexes, V_hexes);
	writeDotMeshHex(hex_mesh_file, V_hexes, hexes);

	// Save final polycube boundary
	if (folder != ""){
		writeDotMeshTet(folder + "/polycube_tets_int.mesh", int_U, tets);
		Eigen::MatrixXd Vbi, Vbf; // initial vs final boundary
		Eigen::MatrixXi Fb;
		igl::readOBJ(folder + "/boundary.obj", Vbi, Fb);
		Vbf.resize(Vbi.rows(), Vbi.cols());

		BndToTetConverter conv(folder + "/tris_to_tets.txt");

		for (int i=0; i<Fb.rows(); i++){
			int f_tet = conv.table_(i);
			int tet_id = f_tet / 4;
			int f_in_tet = f_tet % 4;

			std::vector<std::vector<int>> corres = {
				{1, 3, 2},
				{0, 2, 3},
				{0, 3, 1},
				{0, 1, 2}	
			};
			Vbf.row(Fb(i, 0)) = int_U.row(tets(tet_id, corres[f_in_tet][0]));
			Vbf.row(Fb(i, 1)) = int_U.row(tets(tet_id, corres[f_in_tet][1]));
			Vbf.row(Fb(i, 2)) = int_U.row(tets(tet_id, corres[f_in_tet][2]));
		}

		igl::writeOBJ(folder + "/polycube_surf_int.obj", Vbf, Fb);

	}

	return 0;
}

