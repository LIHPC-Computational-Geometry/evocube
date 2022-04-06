#pragma once

#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/per_face_normals.h>
#include <igl/avg_edge_length.h>

#include "logging.h"

// might want to add normal angle to distance

// Note: such a distance will lead to "round" regions, also it's not mesh independent
// but we're doing polycubes
// how to generate square regions?
// Maybe by considering max(dist_x, dist_y, dist_z) instead?
void triangle_triangle_dist(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, 
                            const Eigen::MatrixXi& TT, 
                            Eigen::MatrixXd& dists){
    dists.resize(F.rows(), 3);

    dists = Eigen::MatrixXd::Constant(F.rows(), 3, -1);

    for (int f_id = 0; f_id < F.rows(); f_id++){
        Eigen::RowVector3d G1 = (V.row(F(f_id, 0)) + V.row(F(f_id, 1)) + V.row(F(f_id, 2))) / 3.0;
        for (int neigh = 0; neigh < 3; neigh ++){
            int neigh_id = TT(f_id, neigh);
            if (neigh_id == -1) {
                coloredPrint("Triangle doesn't have a neighbor", "red");
                continue;
            }
            Eigen::RowVector3d G2 = (V.row(F(neigh_id, 0)) + V.row(F(neigh_id, 1)) + V.row(F(neigh_id, 2))) / 3.0;
            // middle point on shared edge:
            Eigen::RowVector3d mp = (V.row(F(f_id, neigh)) + V.row(F(f_id, (neigh + 1) % 3))) / 2.0;
            dists(f_id, neigh) = (G1 - mp).norm() + (mp - G2).norm();
        }
    }
}

class Evocube {
public:
    Evocube(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
        : V_(V), F_(F){
        
        igl::triangle_triangle_adjacency(F_, TT_);
        igl::vertex_triangle_adjacency(V_.rows(), F_, VT_, VI_);
        igl::per_face_normals(V_, F_, N_);
        triangle_triangle_dist(V, F, TT_, dists_);
        l_avg_ = igl::avg_edge_length(V, F);

    }

    virtual ~Evocube(){
        coloredPrint("An Evocube bites the dust...", "yellow");
    }

    void fillMeshLogInfo(std::string logs_path){
        fillLogInfo("InputTris", "vertices", logs_path, std::to_string(V_.rows()));
        fillLogInfo("InputTris", "faces", logs_path, std::to_string(F_.rows()));
        fillLogInfo("InputTris", "AvgEdgeLength", logs_path, l_avg_);
    }

//private:
    // Set by constructor
    const Eigen::MatrixXd V_; 
    const Eigen::MatrixXi F_; 
    Eigen::MatrixXi TT_;
    Eigen::MatrixXd N_;
    Eigen::MatrixXd dists_;
    std::vector<std::vector<int>> VT_, VI_;
    double l_avg_;

    int timestamp_ = 0;
};