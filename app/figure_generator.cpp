/**
 * @brief Generates figures from a processed mesh. see FIGURE possible values
 * @date 2022-03-18
 */

//file IO
#include <igl/read_triangle_mesh.h>
#include <igl/readOBJ.h>
#include <igl/png/writePNG.h>
//math
#include <igl/PI.h>
//geometry
#include <Eigen/Geometry>
#include <igl/gaussian_curvature.h>
#include <igl/per_corner_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/barycenter.h>
#include <igl/triangle_triangle_adjacency.h>
// embree
#include "MyEmbreeRenderer.h"
#include <igl/embree/ambient_occlusion.h>
//openGL
#include <igl/opengl/glfw/background_window.h>
//SDL
#include <iostream>
#include <fstream>
#include <cassert>

#include "flagging_utils.h"
#include "mesh_io.h"

// **********************************************************
//         PARAMETERS

#define PICTURES_WIDTH  360
#define PICTURES_HEIGHT 360

#define DEFAULT_ZOOM                0.9
#define DEFAULT_TRANSLATION_VECTOR  Eigen::Vector3d(0,0,0)
#define USE_ORTHOGRAPHIC_PROJ       false
#define TILT                        0.6

// #define USE_SPECIAL_ROTATION_MATRIX__BUNNY
// #define USE_SPECIAL_ROTATION_MATRIX__BIMBA
// #define USE_SPECIAL_ROTATION_MATRIX__B52
// #define USE_SPECIAL_ROTATION_MATRIX__MAX
// #define USE_SPECIAL_ROTATION_MATRIX__MOAI
// #define USE_SPECIAL_ROTATION_MATRIX__ONI

#define DEFAULT_DATA_FODLER     "../data/"
#define DEFAULT_MODEL           "bunny"
#define BOUNDARY_PATH(data_folder,mesh_name)        ((data_folder) + (mesh_name) + "/boundary.obj")
#define POLYCUBE_PATH(data_folder,mesh_name)        ((data_folder) + (mesh_name) + "/polycube_final.obj")
#define FLAGGING_PATH(data_folder,mesh_name)        ((data_folder) + (mesh_name) + "/labeling.txt")
#define FEATURE_EDGES_PATH(data_folder,mesh_name)   ((data_folder) + (mesh_name) + "/input_edges.mesh")

#define DEFAULT_FIGURE_TYPE                             FIG_0_POLYCUBE //see enum FIGURE
#define DEFAULT_OUTPUT_PATH(fig_type,number)                ("../temp_picture" + std::to_string(fig_type) + "_" + std::to_string(number) + ".png")
#define OUTPUT_PATH(data_folder,mesh_name,fig_type,number)  ((data_folder) + (mesh_name) + "/fig" + std::to_string(fig_type) + "_" + std::to_string(number) + ".png")

// **********************************************************

int main(int argc, char *argv[]) {

    // usage: ./figure_generator [model] [figure_type] [data_folder]
    // examples :
    //
    // ./figure_generator B21 1 ../custom/data/folder/
    //      Will open the B21 folder in ../custom/data/folder/,
    //      search for the input files,
    //      generate figures of type 1 (see enum FIGURE),
    //      and write the pictures there
    //
    // ./figure_generator B21 1
    //      Same, but will open the DEFAULT_DATA_FODLER/B21 folder
    //
    // ./figure_generator B21
    //      Same, but will generate the DEFAULT_FIGURE_TYPE
    //
    // ./figure_generator
    //      Same, but will use the DEFAULT_MODEL

    std::string mesh_name = DEFAULT_MODEL, data_folder = DEFAULT_DATA_FODLER;
    double zoom_value = DEFAULT_ZOOM;
    Eigen::Vector3d trans_vector = DEFAULT_TRANSLATION_VECTOR;
    int n_figs = -1;

    if (argc >= 2){
        mesh_name = argv[1];
    }

    if (argc >= 4){
        data_folder = argv[3];
    }

    // --- options for the figure types --- //

    enum MESH_MODE {INITIAL,    //input triangular mesh
                    POLYCUBE    //polycube triangular mesh
                    };
    MESH_MODE mesh_mode = INITIAL;

    enum COLORS_MODE {WHITE,            //solid white color
                      LABELLING,        //per triangle label (white/red/blue)
                      NORMALS,          //colors according to the normals
                      NAIVE_LABELLING   //naive labeling, see flagging_utils.cpp normalFlagging()
                      };
    COLORS_MODE colors_mode = WHITE;
    bool colors_face_AO = false;

    enum EDGES_MODE {BOUNDARIES,    //show edges between two different labels
                     FEATURES,      //show feature edges
                     NO_EDGES       //don't show any edges
                     };
    EDGES_MODE edges_mode = NO_EDGES;

    enum ROTATION_MODE {FOUR_POV,       //4 points of view
                        SINGLE_PRESET   //only 1 picture of predefined rotation
                        };
    ROTATION_MODE rotation_mode = FOUR_POV;
    Eigen::Matrix3d rot_preset = Eigen::MatrixXd::Identity(3, 3);

    enum FIGURE {FIG_0_POLYCUBE,
                 FIG_1_INPUT,
                 FIG_2_LABELING_CUSTOM_VIEW,
                 FIG_3_NORMALS,
                 FIG_4_LABELING};
    FIGURE figure = DEFAULT_FIGURE_TYPE;

    if (argc >= 3){
        figure = (FIGURE) std::atoi(argv[2]);
    }
    
    switch(figure) {

        case FIG_0_POLYCUBE:
            colors_mode = LABELLING;
            edges_mode = BOUNDARIES;
            mesh_mode = POLYCUBE;
            n_figs = 4;
            colors_face_AO = true;
            break;
        
        case FIG_1_INPUT:
            colors_mode = WHITE;
            edges_mode = FEATURES;
            mesh_mode = INITIAL;
            n_figs = 4;
            colors_face_AO = true;
            break;

        case FIG_2_LABELING_CUSTOM_VIEW:
            colors_mode = NAIVE_LABELLING;
            mesh_mode = INITIAL;
            n_figs = 1;
            rotation_mode = SINGLE_PRESET;
            rot_preset = Eigen::MatrixXd(3, 3);
            rot_preset << 0.886888, -0.459069, 0.0518283,
                        0.161148,  0.412548,  0.896569,
                        -0.432969, -0.786804,  0.439861;
            break;

        case FIG_3_NORMALS:
            // This one also has slight modifications further down
            colors_mode = NORMALS;
            mesh_mode = INITIAL;
            n_figs = 3;
            rotation_mode = SINGLE_PRESET;
            rot_preset = Eigen::MatrixXd(3, 3);
            rot_preset << 0.886888, -0.459069, 0.0518283,
                        0.161148,  0.412548,  0.896569,
                        -0.432969, -0.786804,  0.439861;
            break;

        case FIG_4_LABELING:
            colors_mode = LABELLING;
            edges_mode = NO_EDGES; //BOUNDARIES;
            mesh_mode = INITIAL;
            n_figs = 4;
            colors_face_AO = true;
            break;
    }

    // ---  --- //

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXi TT;
    Eigen::MatrixXd N;
    Eigen::MatrixXi edges;
    Eigen::MatrixXd V_edges;
    std::string mesh_file, feature_edges_file;

    //open the mesh file according to the mesh_mode
    mesh_file = (mesh_mode == INITIAL) ? BOUNDARY_PATH(data_folder,mesh_name) : POLYCUBE_PATH(data_folder,mesh_name);
    std::cout << "Reading " << mesh_file << "..." << std::endl;
    igl::read_triangle_mesh(mesh_file, V, F);

    //open the feature edges file
    feature_edges_file = FEATURE_EDGES_PATH(data_folder,mesh_name);
    if (edges_mode == FEATURES) {
        if (!std::ifstream(feature_edges_file.c_str()).good()){
            coloredPrint("Feature edges can't be read.", "yellow");
        }
        else {
            readDotMeshEdges(feature_edges_file, V_edges, edges);
        }
    }

    igl::triangle_triangle_adjacency(F, TT);
    igl::per_vertex_normals(V, F, N);

    V.col(0) = V.col(0).array() - V.col(0).mean();
    V.col(1) = V.col(1).array() - V.col(1).mean();
    V.col(2) = V.col(2).array() - V.col(2).mean();

    std::cout << "Mean: " << V.colwise().mean() << std::endl;

    // --------------------------

    Eigen::VectorXd AO;
    igl::embree::ambient_occlusion(V, F, V, N, 30, AO);

    //open flagging file
    std::string flagging_file = FLAGGING_PATH(data_folder,mesh_name);
    std::cout << "Reading " << flagging_file << "..." << std::endl;
    Eigen::VectorXi flagging = openFlagging(flagging_file,F.rows());

    Eigen::VectorXd face_AO(F.rows());
    for (int i=0; i<F.rows(); i++){
        face_AO(i) = (AO(F(i,0)) + AO(F(i,1)) + AO(F(i,2)))/3.0;
        face_AO(i) /= 5.0;
    }
    Eigen::VectorXd smooth_AO = face_AO;
    for (int i=0; i<TT.rows(); i++){
        for (int j=0; j<3; j++){
            smooth_AO(i) += face_AO(TT(i,j));
        }
    }
    smooth_AO = smooth_AO.array() / 4.0;

    //compute colors
    Eigen::MatrixXd colors;
    switch (colors_mode) {
        
        case LABELLING:
            colors = colorsFromFlagging(flagging);
            break;

        case WHITE:
            colors = Eigen::MatrixXd::Constant(F.rows(), 3, 0.8);
            break;

        case NORMALS:
            colors = N.array().abs();
            break;

        case NAIVE_LABELLING: {
            Eigen::VectorXi naive_flagging = normalFlagging(V, F);
            colors = colorsFromFlagging(naive_flagging);
            break; }
            
    }

    if (colors_face_AO){
        if (colors.rows() != smooth_AO.rows()){
            coloredPrint("Error: mixing face-based and vertex-based colors ?", "red");
        }
        colors = colors.colwise() - smooth_AO;
    }

    std::cout << "COLORS set" << std::endl;

    Eigen::MatrixXf begs(0,0);
    Eigen::MatrixXf ends(0,0);
    Eigen::MatrixXf edge_colors(0,0);

    //compute edges
    if (edges_mode == BOUNDARIES){

        Eigen::MatrixXf bnd_flag_colors(3,3);
        bnd_flag_colors.row(0) = Eigen::RowVector3f(196.0/255, 32.0/255, 33.0/255);
        bnd_flag_colors.row(1) = Eigen::RowVector3f(175.0/255, 185.0/255, 178.0/255);;
        bnd_flag_colors.row(2) = Eigen::RowVector3f(5.0/255, 74.0/255, 145.0/255);
        std::vector<std::pair<int, int>> edges_vec;
        std::vector<int> edge_axes;
        for (int i=0; i<F.rows(); i++){
            for (int j=0; j<3; j++){
                if (flagging(i) != flagging(TT(i,j)) && i < TT(i,j)){
                    edges_vec.push_back(std::make_pair(F(i,(j+0)%3), F(i,(j+1)%3)));
                    int axis = 3 - flagToAxis(flagging(i)) - flagToAxis(flagging(TT(i,j)));
                    edge_axes.push_back(axis);
                }
            }
        }
        begs = Eigen::MatrixXf(edges_vec.size(), 3);
        ends = Eigen::MatrixXf(edges_vec.size(), 3);
        edge_colors = Eigen::MatrixXf::Constant(edges_vec.size(), 3, 0.0);
      

        for (int i=0; i<edges_vec.size(); i++){
            begs.row(i) = V.row(edges_vec[i].first).cast<float>();
            ends.row(i) = V.row(edges_vec[i].second).cast<float>();
            //edge_colors(i, edge_axes[i]) = 1.0;
            //edge_colors.row(i) = bnd_flag_colors.row(edge_axes[i]);
        }
    }
    else if (edges_mode == FEATURES){
        assert(( edges.rows()>0 ));
        begs = Eigen::MatrixXf(edges.rows(), 3);
        ends = Eigen::MatrixXf(edges.rows(), 3);
        edge_colors = Eigen::MatrixXf::Constant(edges.rows(), 3, 0.0);

        for (int i=0; i<edges.rows(); i++){
            begs.row(i) = V.row(edges(i,0)).cast<float>();
            ends.row(i) = V.row(edges(i,1)).cast<float>();
            edge_colors.row(i) = Eigen::Vector3f(250.0/255.0, 200.0/255.0, 85.0/255.0);
        }
    }

    // Pre-computing edges per tri - needed if you want to display any edges
    std::vector<std::vector<int>> edges_per_tri;
    if (edges_mode != NO_EDGES){
        std::cout << "Computing edges per tri..." << std::endl;
        double threshold = 20 * 0.3;
        for (int i=0; i<F.rows(); i++){
            edges_per_tri.push_back({});
            Eigen::RowVector3f p1 = ((V.row(F(i,0)) + V.row(F(i,1)) + V.row(F(i,2)))/3.0).cast<float>();
            for (int j=0; j<begs.rows(); j++){
                Eigen::RowVector3f p2 = (begs.row(j) + ends.row(j))/2.0;
                if ((p1-p2).norm() < threshold){
                    edges_per_tri[i].push_back(j);
                } 
            }
        }
        std::cout << "Done computing edges per tri" << std::endl;
    }

    std::cout << "EDGES set" << std::endl;

    // --------------------------

    // embree object
    igl::embree::EmbreeRenderer er;
    std::cout << "Sending mesh to renderer" << std::endl;
    er.set_mesh(V, F, true);
    
    if (edges_mode != NO_EDGES) {
        std::cout << "Sending edges to renderer" << std::endl;
        er.set_edges(begs, ends, edge_colors);
        er.set_edges_per_tri(edges_per_tri);
    }

    std::cout << "Sending colors to renderer" << std::endl;
    er.set_colors(colors);

    std::vector<double> rotX = {45 * TILT, -45 * TILT, -(180 - 45 * TILT), 180 - 45 * TILT};
    std::vector<double> rotY = {-45, 45, 45, -45};
    std::vector<double> rotZ = {0, 0, 180, 180};
    Eigen::Matrix3d rot_matrix;

    for (int pic=0; pic<n_figs; pic++) {
        std::cout << "Pic " << pic << std::endl;

        if (rotation_mode == SINGLE_PRESET) {
            rot_matrix = rot_preset;
        }
        else {
            rot_matrix =  Eigen::AngleAxisd((rotX[pic])*igl::PI/180.0, Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd((rotY[pic])*igl::PI/180.0, Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd((rotZ[pic])*igl::PI/180.0, Eigen::Vector3d::UnitZ());
        }

        #if defined(USE_SPECIAL_ROTATION_MATRIX__BUNNY)
        rot_matrix << 0.857921,  -0.513781, -0.00074251,
                      0.167674,  0.278619,  0.945652,
                      -0.485651, -0.811419, 0.325181;
        #elif defined(USE_SPECIAL_ROTATION_MATRIX__BIMBA)
        rot_matrix << -0.917247, -0.397913, -0.0179899,
                      0.149208,  -0.385121, 0.910725,
                      -0.369317, 0.832675,  0.412623;
        #elif defined(USE_SPECIAL_ROTATION_MATRIX__B52)
        rot_matrix << -0.907966, -0.416895, 0.0423878,
                      -0.218289, 0.384208,  -0.89707,
                      0.357698,  -0.823761, -0.439852;
        #elif defined(USE_SPECIAL_ROTATION_MATRIX__MAX)
        rot_matrix << 0.726065,  0.687486,   -0.0139064,
                      0.0716724, -0.0555499, 0.99588,
                      0.683881,  -0.72407,   -0.0896065;
        #elif defined(USE_SPECIAL_ROTATION_MATRIX__MOAI)
        rot_matrix << 0.621212,  -0.779063, 0.0845922,
                      0.292434,  0.330615,  0.897316,
                      -0.727033, -0.532686, 0.433206;
        #elif defined(USE_SPECIAL_ROTATION_MATRIX__ONI)
        rot_matrix << 0.784854,  -0.619391, -0.0189354,
                      0.120424,  0.122477,  0.985138,
                      -0.607867, -0.77547,  0.170716;
        #endif

        if (figure == FIG_3_NORMALS){
            //for this mode the colors depend on the picture number
            Eigen::MatrixXd new_colors = colors;
            new_colors.col(pic) = new_colors.col(pic).array() * new_colors.col(pic).array(); 
            new_colors.col((pic+1)%3) = Eigen::VectorXd::Constant(new_colors.rows(), 0.0);
            new_colors.col((pic+2)%3) = Eigen::VectorXd::Constant(new_colors.rows(), 0.0);
            if (pic == 1){
                new_colors.col(0) = new_colors.col(1);
                new_colors.col(2) = new_colors.col(1);
            }
            er.set_colors(new_colors);
        }

        er.set_rot(rot_matrix);
        er.set_translation(trans_vector);
        er.set_zoom(zoom_value);
        er.set_orthographic(USE_ORTHOGRAPHIC_PROJ);
        
        //variables for the Red/Green/Blue/Alpha channels
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(PICTURES_WIDTH, PICTURES_HEIGHT);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(PICTURES_WIDTH, PICTURES_HEIGHT);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(PICTURES_WIDTH, PICTURES_HEIGHT);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(PICTURES_WIDTH, PICTURES_HEIGHT);

        // render view using embree
        std::cout << "Rendering..." << std::endl;
        er.render_buffer(R,G,B,A);
        std::cout << "Rendered !" << std::endl;

        // save to PNG file
        std::string png_file = (argc >= 2) ? OUTPUT_PATH(data_folder,mesh_name,(int) figure,pic) : DEFAULT_OUTPUT_PATH((int) figure,pic);
        igl::png::writePNG(R, G, B, A, png_file);
        std::cout << "Rendered scene saved to "<< png_file << std::endl;
    }

    return 0;
}
