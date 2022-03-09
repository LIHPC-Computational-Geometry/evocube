// HexExWrapper.h
// https://github.com/fprotais/polycube_withHexEx with Eigen instead of ultimaille

#pragma once
#include <Eigen/Core>

extern const Eigen::MatrixXi TETRAHEDRON_TO_FACES;

//get the vertex id of the v_th vertex of the f_th face of the t_th tetrahedron
#define TETRAHEDRON_TO_VERTEX(tets,t,f,v) (tets(t,TETRAHEDRON_TO_FACES(f,v)))

bool run_HexEx(const Eigen::MatrixXi& tets, const Eigen::MatrixXd& V_tets, const Eigen::MatrixXd& corner_param, Eigen::MatrixXi& hexes, Eigen::MatrixXd& V_hexes);
