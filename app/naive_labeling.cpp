
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/per_face_normals.h>
#include "bnd_to_tet.h"
#include "flagging_utils.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/S1/boundary.obj", V, F);
    
    Eigen::MatrixXd N;
    igl::per_face_normals(V, F, N);

    Eigen::VectorXi labeling(F.rows());
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
    }

    //Eigen::MatrixXd colors = colorsFromFlagging(labeling);

    saveFlagging("../data/S1/labeling.txt", labeling);
    int n_tets = 113683;
    coloredPrint("TODO n_tets hardcoded", "red");
    saveFlaggingOnTets("../data/S1/labeling_on_tets.txt", "../data/S1/tris_to_tets.txt", n_tets, labeling);

    //BndToTetConverter conv("../from_tris_to_tets.txt");    
}