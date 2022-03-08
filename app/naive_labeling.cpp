
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>

#include "tet_boundary.h"
#include "flagging_utils.h"
#include "graphcut_labeling.h"

int main(int argc, char *argv[]){

    std::string folder = "../data/testf/screwdriver_input_tri/";

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(folder + "boundary.obj", V, F);
    
    Eigen::MatrixXd N;
    igl::per_face_normals(V, F, N);

    Eigen::MatrixXi TT;
    igl::triangle_triangle_adjacency(F, TT);

    /*Eigen::VectorXi labeling(F.rows());
    for (int f_id = 0; f_id < F.rows(); f_id ++){
        double best_match = -1;
        double best_label = -1;
        double match;

        for (int axis = 0; axis < 3; axis ++){
            Eigen::RowVector3d vec(0.0, 0.0, 0.0);
            vec(axis) = 1.0;
            match = N.row(f_id).dot(vec);
            if (match > best_match){
                best_match = match;
                best_label = 2 * axis;
            }

            vec(axis) = -1;
            match = N.row(f_id).dot(vec);
            if (match > best_match){
                best_match = match;
                best_label = 2 * axis + 1;
            }
        }

        labeling(f_id) = best_label;
    }*/

    Eigen::VectorXi labeling = graphcutFlagging(V, F, N, TT, 3, 1);

    //Eigen::MatrixXd colors = colorsFromFlagging(labeling);

    saveFlagging(folder + "labeling.txt", labeling);
    saveFlaggingOnTets(folder + "labeling_on_tets.txt", folder + "tris_to_tets.txt", labeling);

}