#pragma once

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <cassert>

#include "logging.h"

Eigen::MatrixXd axesMatrix();

inline int flagToAxis(const int flag){ //just use flag/2 instead
    assert(( (0<=flag) && (flag<=5) ));
    //  flag | axis
    //  +X=0 | X=0
    //  -X=1 | X=0
    //  +Y=2 | Y=1
    //  -Y=3 | Y=1
    //  +Z=4 | Z=2
    //  -Z=5 | Z=2
    return flag/2;//integer division
}

int oppositeLabel(int label);

Eigen::VectorXi openFlagging(std::string file_name, int expected_size);

Eigen::MatrixXd colorsFromFlagging(const Eigen::VectorXi& flagging);

void saveFlagging(std::string file_name, const Eigen::VectorXi& assignment);

void saveFlaggingOnTets(std::string file_name, std::string tris_to_tets_path, const Eigen::VectorXi& assignment);

Eigen::VectorXi normalFlagging(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

// % of same labels between 2 labeling vectors
// if different lengths, return -1.0
double flaggingSimilarity(const Eigen::VectorXi& labeling1, const Eigen::VectorXi& labeling2);