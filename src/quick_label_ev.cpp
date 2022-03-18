#include "quick_label_ev.h"

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

#include <Eigen/Sparse>
#include <OpenNL_psm.h>

#include "disjointset.h"
#include "distortion.h"
#include "flagging_utils.h"

//#define DEBUG_QUICK_LABEL_EV
#ifdef DEBUG_QUICK_LABEL_EV
#include <igl/writeOBJ.h>
#endif

#define FOR(i, n) for(int i = 0; i < n; i++)

Eigen::MatrixXd QuickLabelEv::LDLTDeformation(const Eigen::VectorXi& labeling) const {
	int N_vs = V_.rows();
	int N_fs = labeling.rows();
	auto begin = std::chrono::steady_clock::now();

	DisjointSet ds(3 * N_vs);
	FOR(f, N_fs) {
		int dim = labeling(f) / 2;
		FOR(fv, 3) ds.merge(dim * N_vs + F_(f, fv), dim * N_vs + F_(f, (fv + 1) % 3));
	}

	std::vector<int> idmap;
	int nb_variables = ds.get_sets_id(idmap);

	Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Eigen::VectorXd x;

	{ // makeSparseMatrix
		int expected_eqs = 3 * 2 * N_fs;
		int n_coords = nb_variables; // TODO
		A.resize(expected_eqs, n_coords);
		x.resize(n_coords);
		b.resize(expected_eqs);

		std::vector<Eigen::Triplet<double>> triplet_list; 
		triplet_list.resize(2 * expected_eqs);
		int next_equation_id = 0;

		FOR(f, N_fs) {
			int dim = labeling(f) / 2;
			FOR(d, 3) {
				if (dim == d) continue;
				FOR(fv, 3) {
					triplet_list[2 * next_equation_id] = Eigen::Triplet<double>(next_equation_id, idmap[d * N_vs + F_(f, fv)], 1.0);
					triplet_list[2 * next_equation_id + 1] = Eigen::Triplet<double>(next_equation_id, idmap[d * N_vs + F_(f, (fv + 1) % 3)], -1.0);
					b(next_equation_id) = V_(F_(f, fv), d) - V_(F_(f, (fv + 1) % 3), d);
					next_equation_id ++;
				}
			}
		}
		A.setFromTriplets(triplet_list.begin(), triplet_list.end());
	}

	auto buildend = std::chrono::steady_clock::now();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	/*Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
	solver.setMaxIterations(20);
	solver.setTolerance(0.001);*/
	solver.compute(A.transpose() * A);
	if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: decomposition failed" << std::endl;
    }

    x = solver.solve(A.transpose() * b);
	if(solver.info() != Eigen::Success) {
        std::cout << "ERROR, solving failed: ";
        if(solver.info() == Eigen::NumericalIssue) 
            std::cout << "NumericalIssue" << std::endl;
        if(solver.info() == Eigen::NoConvergence) 
            std::cout << "NoConvergence" << std::endl;
        if(solver.info() == Eigen::InvalidInput) 
            std::cout << "InvalidInput" << std::endl;
		//x = Eigen::VectorXd::Constant(x.rows(), 10e6);
    }
	
	Eigen::MatrixXd def_V(N_vs, 3);
	FOR(v, N_vs) FOR(d, 3) def_V(v, d) = x(idmap[d * N_vs + v]);

	auto solveend = std::chrono::steady_clock::now();

	double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - begin).count()) / 1000;
	double buildtime = double(std::chrono::duration_cast<std::chrono::milliseconds> (buildend - begin).count()) / 1000;
	double solvetime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - buildend).count()) / 1000;
	//std::cerr << "LDLT build + solve: " << buildtime << " + " << solvetime << " = " << totaltime << "sec." << std::endl;

	#ifdef DEBUG_QUICK_LABEL_EV
	//igl::writeOBJ("../debug_quick_label_ev.obj", def_V, F_);
	std::cout << "hardDeformation done" << std::endl;
	#endif

	return def_V;
}

Eigen::MatrixXd QuickLabelEv::hardDeformation(const Eigen::VectorXi& labeling) const {
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
	Eigen::MatrixXd def_V(N_vs, 3);
	FOR(v, N_vs) FOR(d, 3) def_V(v, d) = nlGetVariable(idmap[d * N_vs + v]);
	nlDeleteContext(context);

	auto solveend = std::chrono::steady_clock::now();

	double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - begin).count()) / 1000;
	double buildtime = double(std::chrono::duration_cast<std::chrono::milliseconds> (buildend - begin).count()) / 1000;
	double solvetime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - buildend).count()) / 1000;
	std::cerr << "Fast bnd poly build + solve: " << buildtime << " + " << solvetime << " = " << totaltime << "sec." << std::endl;

	#ifdef DEBUG_QUICK_LABEL_EV
	//igl::writeOBJ("../debug_quick_label_ev.obj", def_V, F_);
	std::cout << "hardDeformation done" << std::endl;
	#endif

	return def_V;
}

Eigen::MatrixXd QuickLabelEv::computeDeformedV(const Eigen::VectorXi& labeling) const {
	return LDLTDeformation(labeling);
}

double QuickLabelEv::evaluate(const Eigen::VectorXi& labeling) const {
	int n_fail_invert;
	return evaluate(labeling, n_fail_invert);
}

double QuickLabelEv::evaluate(const Eigen::VectorXi& labeling, int& n_fail_invert) const {
	
	Eigen::MatrixXd def_V = LDLTDeformation(labeling);

    Eigen::MatrixXd N_def;
    igl::per_face_normals(def_V, F_, N_def);

    Eigen::VectorXd disto;
    computeDisto(V_, def_V, F_, N_, N_def, disto);

	disto = disto * disto; // OPTIONAL
	    
    double final_disto = integrateDistortion(A_, disto);

	n_fail_invert = 0;
	Eigen::MatrixXd axes_matrix = axesMatrix();
	for (int i=0; i<labeling.rows(); i++){
		if (N_def.row(i).dot(axes_matrix.row(labeling(i))) < 0){
			n_fail_invert ++;
		}
	}

	//double invert_coeff = 0 * 100000.0 * 1.0 / F_.rows();
	double invert_coeff = 100000.0 * 1.0 / F_.rows();
	final_disto += invert_coeff * static_cast<double>(n_fail_invert);
	
    /*std::cout << "n_fail_invert: " << n_fail_invert << std::endl;
    std::cout << "final_disto: " << final_disto << std::endl;
    std::cout << "disto.maxCoeff(): " << disto.maxCoeff() << std::endl;
	std::cout << "A_.minCoeff() : " << A_.minCoeff() << std::endl;*/
	return final_disto;
}

Eigen::MatrixXd QuickLabelEv::distoAboveThreshold(const Eigen::VectorXi& labeling, double threshold) const {
	Eigen::MatrixXd def_V = LDLTDeformation(labeling);

    Eigen::MatrixXd N_def;
    igl::per_face_normals(def_V, F_, N_def);

    Eigen::VectorXd disto;
    computeDisto(V_, def_V, F_, N_, N_def, disto);

	Eigen::MatrixXd colors = Eigen::MatrixXd::Constant(labeling.rows(), 3, 1.0);

	int n_fail_threshold = 0;
	for (int i=0; i<labeling.rows(); i++){
		if (disto(i) > threshold){
			n_fail_threshold ++;
			colors.row(i) = Eigen::RowVector3d(1.0, 0, 0);
		}
	}

	coloredPrint("Elements above threshold: " + std::to_string(n_fail_threshold), "cyan");
	return colors;
}