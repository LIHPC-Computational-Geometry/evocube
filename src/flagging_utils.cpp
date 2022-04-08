
#include "flagging_utils.h"
#include <Eigen/Geometry>

Eigen::MatrixXd axesMatrix(){
    Eigen::MatrixXd axes = Eigen::MatrixXd::Zero(6,3);
    for (int axis=0; axis<3; axis++){
        for (int dir=0; dir<2; dir++){
            axes(2*axis + dir, axis) = 1 - 2*dir;
        }
    }
    return axes;
}

int oppositeLabel(int label){
    if (label % 2) return label -1;
    return label +1;
}

Eigen::VectorXi openFlagging(std::string file_name, int expected_size){
    Eigen::VectorXi assignment(expected_size);
    if (expected_size == 0) return assignment;
    std::ifstream input_file(file_name);
    if (!input_file.is_open()) {
        coloredPrint("Could not open assignment file: " + file_name, "red");
        assignment = Eigen::VectorXi::Zero(expected_size);
    }
    else {
        int number;
        int n=0;
        while (input_file >> number) {
            assignment(n) = number;
            n++;
        }
        if (n != expected_size){
            std::cout << "Error reading assignment" << std::endl;
            assignment = Eigen::VectorXi::Zero(expected_size);
            return assignment;
        }
    }

    return assignment;
}

Eigen::MatrixXd colorsFromFlagging(const Eigen::VectorXi& flagging){
    if (flagging.rows() == 0){
        std::cout << "WARNING: computing colors from empty flagging" << std::endl;
    }
    Eigen::MatrixXd flagging_colors(flagging.rows(), 3);

    Eigen::MatrixXd color_map = Eigen::MatrixXd::Zero(7,3);
    for (int axis=0; axis<3; axis++){
        for (int dir=0; dir<2; dir++){
            color_map(2*axis + dir, axis) = 1 - 0.4*dir;
        }
    }

    // replace green with white
    color_map.row(2) = Eigen::RowVector3d(1.0, 1.0, 1.0);
    color_map.row(3) = Eigen::RowVector3d(0.7, 0.7, 0.7);


    double diff = 40.0;
    color_map.row(0) = Eigen::RowVector3d(196.0/255, 32.0/255, 33.0/255);
    //color_map.row(1) = Eigen::RowVector3d(166.0/255, 15.0/255, 16.0/255);
    color_map.row(2) = Eigen::RowVector3d(235.0/255, 245.0/255, 238.0/255);
    color_map.row(2) = Eigen::RowVector3d(235.0/255, 235.0/255, 235.0/255);
    //color_map.row(3) = Eigen::RowVector3d(185.0/255, 195.0/255, 198.0/255);
    color_map.row(4) = Eigen::RowVector3d(5.0/255, 74.0/255, 145.0/255);
    color_map.row(4) = Eigen::RowVector3d(5.0/255, 90.0/255, 195.0/255);
    //color_map.row(5) = Eigen::RowVector3d(0.0/255, 44.0/255, 115.0/255);

    // German theme
    //color_map.row(0) = Eigen::RowVector3d(0.0/255, 0.0/255, 0.0/255);
    //color_map.row(2) = Eigen::RowVector3d(221.0/255, 0.0/255, 0.0/255);
    //color_map.row(4) = Eigen::RowVector3d(255.0/255, 206.0/255, 0.0/255);


    // updated theme
    /*
    color_map.row(0) = Eigen::RowVector3d(255.0/255, 15.0/255, 60.0/255);
    color_map.row(2) = Eigen::RowVector3d(245.0/255, 245.0/255, 245.0/255);
    color_map.row(4) = Eigen::RowVector3d(0.0/255, 130.0/255, 255.0/255);
    //*/

    color_map.row(1) = color_map.row(0).array() - 0.9 * diff/255.0;
    color_map.row(3) = color_map.row(2).array() - 0.6 * diff/255.0;
    color_map.row(5) = color_map.row(4).array() - 0.7 * diff/255.0;

    color_map = color_map.array() * 1.1;


    // TODO REMOVE: flipping for a figure
    bool reverse_whites = false;
    if (reverse_whites){
        color_map.row(3) = Eigen::RowVector3d(235.0/255, 235.0/255, 235.0/255);
        color_map.row(2) = color_map.row(3).array() - 1.2 * diff/255.0;
    }

    for (int i=0; i<flagging.rows(); i++){
        flagging_colors.row(i) = color_map.row(flagging(i));
    }
    return flagging_colors;
}


void saveFlagging(std::string file_name, const Eigen::VectorXi& assignment){
    // TODO improve this, use binary
    std::ofstream out_assig;
    out_assig.open(file_name);
    for (int i=0; i<assignment.rows(); i++) {
        out_assig << assignment(i) << "\n";
    }
    out_assig.close();
}

#include "tet_boundary.h"

void saveFlaggingOnTets(std::string file_name, std::string tris_to_tets_path, const Eigen::VectorXi& assignment){

    BndToTetConverter conv(tris_to_tets_path);
    int n_tets = conv.n_tets_;

    Eigen::VectorXi assigment_tets = Eigen::VectorXi::Constant(4 * n_tets, -1);
    if (conv.table_.rows() != assignment.rows()){
        std::cout << "ERROR sizes don't match in saveFlaggingOnTets" << std::endl;
        std::cout << conv.table_.rows() << " vs " << assignment.rows() << std::endl;
    }

    for (int i=0; i<assignment.rows(); i++){
        assigment_tets(conv.table_(i)) = assignment(i);
    }

    std::ofstream out_assig;
    out_assig.open(file_name);
    for (int i=0; i<assigment_tets.rows(); i++) {
        out_assig << assigment_tets(i) << "\n";
    }
    out_assig.close();
}

Eigen::VectorXi normalFlagging(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F){
    Eigen::MatrixXd axes = axesMatrix();
    Eigen::VectorXi flagging = Eigen::VectorXi::Zero(F.rows());

    for (int i=0; i<flagging.rows(); i++){
        Eigen::RowVector3d seg1 = V.row(F(i,1))-V.row(F(i,0)); 
        Eigen::RowVector3d seg2 = V.row(F(i,2))-V.row(F(i,0)); 
        Eigen::RowVector3d cross_p = seg1.cross(seg2);
        cross_p /= cross_p.norm();

        double best = 0;
        double mini = cross_p.dot(axes.row(0));
        for (int ax=1; ax<axes.rows(); ax++){
            if (cross_p.dot(axes.row(ax)) > mini){
                best = ax;
                mini = cross_p.dot(axes.row(ax));
            }
        }
        flagging(i) = best;
    }
    return flagging;
}

double flaggingSimilarity(const Eigen::VectorXi& labeling1, const Eigen::VectorXi& labeling2) {
    if(labeling1.rows()!=labeling2.rows()) {
        return -1.0;
    }
    int same_label_counter = 0;
    for(int face_number = 0; face_number < labeling1.rows(); face_number++) {
        if(labeling1(face_number) == labeling2(face_number)) {
            same_label_counter++;
        }
    }
    return ((double) same_label_counter) / labeling1.rows();
}