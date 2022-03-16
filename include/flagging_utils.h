#pragma once

#include <Eigen/Core>
#include <iostream>
#include <fstream>

#include "logging.h"

Eigen::MatrixXd axesMatrix();

int oppositeLabel(int label);

Eigen::VectorXi openFlagging(std::string file_name, int expected_size);

Eigen::MatrixXd colorsFromFlagging(const Eigen::VectorXi& flagging);

void saveFlagging(std::string file_name, const Eigen::VectorXi& assignment);

void saveFlaggingOnTets(std::string file_name, std::string tris_to_tets_path, const Eigen::VectorXi& assignment);

Eigen::VectorXi normalFlagging(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);