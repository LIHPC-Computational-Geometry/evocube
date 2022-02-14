#pragma once

// see https://github.com/fprotais/fastbndpolycube

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

#include <OpenNL_psm.h>
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>

#include "disjointset.h"
#include "distortion.h"

#define DEBUG_QUICK_LABEL_EV
#ifdef DEBUG_QUICK_LABEL_EV
#include <igl/writeOBJ.h>
#endif

#define FOR(i, n) for(int i = 0; i < n; i++)

class QuickLabelEv {
public:
	QuickLabelEv(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V_(V), F_(F) {
		def_V_ = Eigen::MatrixXd(V.rows(), V.cols()); 
		igl::per_face_normals(V_, F_, N_);
		igl::doublearea(V_, F_, A_);
		A_ /= 2.0;
	}
	double evaluate(const Eigen::VectorXi& labeling);

	Eigen::MatrixXd getDefV(){return def_V_;}

private:
	// Set by constructor
	const Eigen::MatrixXd V_;
	const Eigen::MatrixXi F_;
	Eigen::MatrixXd N_;
	Eigen::VectorXd A_;

	// Each evaluation computes a new deformation
	Eigen::MatrixXd def_V_;
	
	void hardDeformation(const Eigen::VectorXi& labeling);

};

void QuickLabelEv::hardDeformation(const Eigen::VectorXi& labeling) {
	int N_vs = V_.rows();
	int N_fs = labeling.rows();

	auto begin = std::chrono::steady_clock::now();
	auto context = nlNewContext();

	DisjointSet ds(3 * N_vs);
	FOR(f, N_fs) {
		int dim = labeling(f) / 2;
		FOR(fv, 3) ds.merge(dim * N_vs + F_(f, fv), dim * N_vs + F_(f, (fv + 1) % 3));
	}

	std::vector<int> idmap;
	int nb_variables = ds.get_sets_id(idmap);

	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(nb_variables));
	//nlEnable(NL_VERBOSE);
	nlBegin(NL_SYSTEM);

	nlBegin(NL_MATRIX);
	FOR(f, N_fs) {
		int dim = labeling(f) / 2;
		FOR(d, 3) {
			if (dim == d) continue;
			FOR(fv, 3) {
				nlRowScaling(A_(f));
				nlBegin(NL_ROW);
				nlCoefficient(idmap[d * N_vs + F_(f, fv)], 1);
				nlCoefficient(idmap[d * N_vs + F_(f, (fv + 1) % 3)], -1);
				nlRightHandSide(V_(F_(f, fv), d) - V_(F_(f, (fv + 1) % 3), d));
				nlEnd(NL_ROW);
			}
		}
	}
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);

	auto buildend = std::chrono::steady_clock::now();

	nlSolve();
	FOR(v, N_vs) FOR(d, 3) def_V_(v, d) = nlGetVariable(idmap[d * N_vs + v]);
	nlDeleteContext(context);
	auto solveend = std::chrono::steady_clock::now();

	double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - begin).count()) / 1000;
	double buildtime = double(std::chrono::duration_cast<std::chrono::milliseconds> (buildend - begin).count()) / 1000;
	double solvetime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - buildend).count()) / 1000;
	std::cerr << "Fast bnd poly build + solve: " << buildtime << " + " << solvetime << " = " << totaltime << "sec." << std::endl;


	#ifdef DEBUG_QUICK_LABEL_EV
	igl::writeOBJ("../debug_quick_label_ev.obj", def_V_, F_);
	#endif

}

double QuickLabelEv::evaluate(const Eigen::VectorXi& labeling) {
	
	hardDeformation(labeling);


    Eigen::MatrixXd N_def;
    igl::per_face_normals(def_V_, F_, N_def);

    Eigen::VectorXd disto;
    computeDisto(V_, def_V_, F_, N_, N_def, disto);
    //std::cout << "disto: " << disto << std::endl;    
	    
    double final_disto = integrateDistortion(V_, F_, disto);
    std::cout << "final_disto: " << final_disto << std::endl;
    std::cout << "disto.maxCoeff(): " << disto.maxCoeff() << std::endl;

	return final_disto;

}