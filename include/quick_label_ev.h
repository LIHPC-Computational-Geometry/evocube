#pragma once

// see https://github.com/fprotais/fastbndpolycube

#include <igl/per_face_normals.h>
#include <igl/doublearea.h>

class QuickLabelEv {
public:
	QuickLabelEv(){}
	QuickLabelEv(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V_(V), F_(F) {
		igl::per_face_normals(V_, F_, N_);
		igl::doublearea(V_, F_, A_);
		A_ /= 2.0;
	}

	double evaluate(const Eigen::VectorXi& labeling) const;
	
	Eigen::MatrixXd computeDeformedV(const Eigen::VectorXi& labeling) const;
	

	//Eigen::MatrixXd getDefV(){return def_V_;}

	Eigen::MatrixXd distoAboveThreshold(const Eigen::VectorXi& labeling, double threshold) const;

	/*void copyData(Eigen::MatrixXd& V, Eigen::MatrixXi& F,	
				  Eigen::MatrixXd& N, Eigen::VectorXd& A,
				  Eigen::MatrixXd& def_V) const {
		V = V_;
		F = F_;
		N = N_;
		A = A_; 
		def_V = def_V_;
	};*/

private:
	// Set by constructor
	const Eigen::MatrixXd V_;
	const Eigen::MatrixXi F_;
	Eigen::MatrixXd N_;
	Eigen::VectorXd A_;

	// Each evaluation computes a new deformation
	// Eigen::MatrixXd def_V_;
	
	Eigen::MatrixXd hardDeformation(const Eigen::VectorXi& labeling) const;

};