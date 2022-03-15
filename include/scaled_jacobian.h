#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>

// /!\ WARNING : MEDIT convention here (.mesh). Different from UM, OVM conventions
//      5-------6
//     /|      /|
//    / |     / |
//   1-------2  |
//   |  4----|--7
//   | /     | /
//   |/      |/
//   0-------3
constexpr int HEX_CORNER_SPLITTING[8][4] = {
	{0,3,4,1},//bottom-front-left corner
	{3,7,0,2},//bottom-front-right corner
	{4,0,7,5},//bottom-back-left corner
	{7,4,3,6},//bottom-back-right corner
	{1,5,2,0},//top-front-left corner
	{2,1,6,3},//top-front-right corner
	{5,6,1,4},//top-back-left corner
	{6,2,5,7},//top-back-right corner
};

/**
 * @brief Compte the min Scaled Jacobian on an hexahedral mesh
 * @param hexes             #hexes by 8 matrix of vertex id (integers). MEDIT vertex ordering.
 * @param V_hexes           #vertices by 3 matrix of 3D coordinates (doubles)
 * @param per_cell_min_sj   Output. #hexes by 1 vector of min Scaled Jacobian values (doubles)
 * @return the global min Scaled Jacobian
*/
double compute_min_scaled_jacobian(const Eigen::MatrixXi& hexes, const Eigen::MatrixXd& V_hexes, Eigen::VectorXd& per_cell_min_sj) {

    assert(( hexes.cols() == 8 ));
    assert(( hexes.maxCoeff() <= V_hexes.rows() ));

    per_cell_min_sj.resize(hexes.rows());
	double glob_min = 1;
	double min_sj = 1;
	Eigen::Vector3d n1, n2, n3;

	for (int h = 0; h < hexes.rows(); h++) { //for each hexahedron
		
		min_sj = 1;
		for (int hv = 0; hv < 8; hv++) { //for each vertex of the current hexahedron

            //get the coordinates of the 4 vertices constituting the corner n°hv
            Eigen::MatrixXd v(4,3);
			for (int i = 0; i < 4; i++) {
                v.row(i) = V_hexes.row(hexes(h,HEX_CORNER_SPLITTING[hv][i]));
            }

			n1 = (v.row(1) - v.row(0)).normalized();
            n2 = (v.row(2) - v.row(0)).normalized();
            n3 = (v.row(3) - v.row(0)).normalized();
			min_sj = std::min(min_sj, n3.dot(n1.cross(n2)));
		}
		per_cell_min_sj(h) = min_sj;
		glob_min = std::min(glob_min, min_sj);
	}

	return glob_min;
}