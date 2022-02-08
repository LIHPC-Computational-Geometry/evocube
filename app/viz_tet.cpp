#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/barycenter.h>

#include "mesh_io.h"
#include "flagging_utils.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXd tet_colors;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
    using namespace std;
    using namespace Eigen;

    if (key >= '1' && key <= '9')
    {
        double t = double((key - '1')+1) / 9.0;

        VectorXd v = B.col(2).array() - B.col(2).minCoeff();
        v /= v.col(0).maxCoeff();

        vector<int> s;

        for (unsigned i=0; i<v.size();++i)
        if (v(i) < t || true)
            s.push_back(i);

        MatrixXd V_temp(s.size()*4,3);
        MatrixXi F_temp(s.size()*4,3);

        for (unsigned i=0; i<s.size();++i)
        {
            V_temp.row(i*4+0) = TV.row(TT(s[i],0));
            V_temp.row(i*4+1) = TV.row(TT(s[i],1));
            V_temp.row(i*4+2) = TV.row(TT(s[i],2));
            V_temp.row(i*4+3) = TV.row(TT(s[i],3));
            F_temp.row(i*4+0) << (i*4)+1, (i*4)+2, (i*4)+3;
            F_temp.row(i*4+1) << (i*4)+0, (i*4)+3, (i*4)+2;
            F_temp.row(i*4+2) << (i*4)+0, (i*4)+1, (i*4)+3;
            F_temp.row(i*4+3) << (i*4)+0, (i*4)+2, (i*4)+1;
        }

        viewer.data().clear();
        viewer.data().set_mesh(V_temp,F_temp);
        viewer.data().set_face_based(true);
        viewer.data().set_colors(tet_colors);

        std::cout << "F_temp.rows(): " << F_temp.rows() << std::endl;
        std::cout << "tet_colors.rows(): " << tet_colors.rows() << std::endl;
    }

    return false;
}

int main(int argc, char *argv[]){
    using namespace Eigen;
    using namespace std;

    std::string folder = "../data/mambo/base1-Part_1.hh.sat/";
    folder = "../data/S1/";
    std::string input_tets = folder + "tetra.mesh";
    readDotMeshTet(input_tets, TV, TT);

    Eigen::VectorXi tet_labeling = openFlagging(folder + "labeling_on_tets.txt", 4 * TT.rows());
    tet_colors = colorsFromFlagging(tet_labeling);

    // Compute barycenters
    igl::barycenter(TV,TT,B);

    // Plot the generated mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_down = &key_down;
    key_down(viewer,'5',0);
    viewer.launch();
}
