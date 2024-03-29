#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerCore.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/writeOBJ.h>
#include <igl/png/readPNG.h>
#include <filesystem>
#include <imgui.h>

#include "mesh_io.h"
#include "flagging_utils.h"
#include "logging.h"
#include "tet_boundary.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXd tet_colors;

bool display_labeling = true;
bool show_only_invalid = false;
std::string folder;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
    using namespace std;
    using namespace Eigen;

    if (key >= '1' && key <= '9')
    {
        double t = double((key - '1')+1) / 9.0;

        MatrixXd V_temp(TT.rows()*4, 3);
        MatrixXi F_temp(TT.rows()*4, 3);

        for (unsigned i=0; i<TT.rows();++i)
        {
            V_temp.row(i*4+0) = TV.row(TT(i,0));
            V_temp.row(i*4+1) = TV.row(TT(i,1));
            V_temp.row(i*4+2) = TV.row(TT(i,2));
            V_temp.row(i*4+3) = TV.row(TT(i,3));
            F_temp.row(i*4+0) << (i*4)+1, (i*4)+2, (i*4)+3;
            F_temp.row(i*4+1) << (i*4)+0, (i*4)+3, (i*4)+2;
            F_temp.row(i*4+2) << (i*4)+0, (i*4)+1, (i*4)+3;
            F_temp.row(i*4+3) << (i*4)+0, (i*4)+2, (i*4)+1;
        }

        V_temp.col(0) = V_temp.col(0).array() - V_temp.col(0).minCoeff();
        V_temp.col(1) = V_temp.col(1).array() - V_temp.col(1).minCoeff();
        V_temp.col(2) = V_temp.col(2).array() - V_temp.col(2).minCoeff(); 

        std::cout << "F_temp.rows(): " << F_temp.rows() << std::endl;
        viewer.data().clear();
        viewer.data().set_mesh(V_temp, F_temp);
        viewer.core().align_camera_center(V_temp);
        //viewer.core().get_scale_and_shift_to_fit_mesh(V_temp, F_temp, viewer.core().camera_zoom, viewer.core().camera_translation);
        //viewer.core().camera_zoom = 100.0 * viewer.core().camera_zoom;
        viewer.data().set_face_based(true);

        if (display_labeling){
            Eigen::VectorXi tet_labeling = openFlagging(folder + "/labeling_on_tets.txt", 4 * TT.rows());
            tet_colors = colorsFromFlagging(tet_labeling);
            viewer.data().set_colors(tet_colors);
            std::cout << "tet_colors.rows(): " << tet_colors.rows() << std::endl;
        }
    }

    return false;
}

std::vector<std::string> listDirs(const std::string& s){
    std::vector<std::string> r;
    for(auto& p : std::filesystem::directory_iterator(s))
        if (p.is_directory()){
            if (!show_only_invalid || readLogsValue("LabelingFinal", "InvalidTotal", p.path().string() + "/logs.json") != "0")
                r.push_back(p.path().string());
        }
    std::sort(r.begin(), r.end());
    return r;
}

int main(int argc, char *argv[]){
    using namespace Eigen;
    using namespace std;

    folder = "../data/DATASET2/";
    std::vector<std::string> dirs1 = listDirs(folder);
    folder = dirs1[0];

    std::vector<std::string> dirs1_names;
    for (int i=0; i<dirs1.size(); i++){
        std::string path = dirs1[i];
        dirs1_names.push_back(path.substr(path.find_last_of("/\\") + 1,
                                 path.size()));
    }
    

    std::vector<std::string> dirs2 = listDirs(folder);
    folder = dirs2[0];
    std::string input_tets = folder + "/tetra.mesh";

    std::vector<std::string> dirs2_names;
    for (int i=0; i<dirs2.size(); i++){
        std::string path = dirs2[i];
        dirs2_names.push_back(path.substr(path.find_last_of("/\\") + 1,
                                 path.size()));
    }

    int picked_dir1 = 0;
    int picked_dir2 = 0;

    if (argc > 1) input_tets = argv[1];
    readDotMeshTet(input_tets, TV, TT);

    // Plot the generated mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_down = &key_down;
    key_down(viewer,'5',0);


    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    menu.callback_draw_custom_window = [&]() {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(350, -1), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Viz tet")) {
            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;

            ImGui::Checkbox("Display only invalid meshes", &show_only_invalid);
            
            /*if (ImGui::ListBox("Folder 1", &folder_id, folder1, IM_ARRAYSIZE(folder1), 40)){
                int z = 0;
                key_down(viewer,'5',0);
            }*/

            /*const char* dir_items[4];
            dir_items[0] = "Mesh";
            dir_items[1] = "Woven texture";
            dir_items[2] = "Stretch";
            dir_items[3] = "Alignment";
            int display_mode = 0;
            if (ImGui::ListBox("Display", &display_mode, dir_items, IM_ARRAYSIZE(dir_items), 4)){
            }*/

            const char* dir_items[dirs1_names.size()];
            for (int dir=0; dir<dirs1_names.size(); dir++){
                dir_items[dir] = dirs1_names[dir].c_str();
            }
            if (ImGui::ListBox("Main dir", &picked_dir1, dir_items, IM_ARRAYSIZE(dir_items), 6)){
                /*try {
                    std::cout << "Reading "<<file_names[item_current]<<"..." << std::endl;
                    visualizeHexex(file_names[item_current]);
                }
                catch (...) {
                    std::cout << "Reading error !!" << std::endl;
                }*/
                
                dirs2 = listDirs(dirs1[picked_dir1]);
                dirs2_names.clear();
                for (int i=0; i<dirs2.size(); i++){
                    std::string path = dirs2[i];
                    dirs2_names.push_back(path.substr(path.find_last_of("/\\") + 1,
                                            path.size()));
                }
            }

            const char* dir2_items[dirs2_names.size()];
            for (int dir=0; dir<dirs2_names.size(); dir++){
                dir2_items[dir] = dirs2_names[dir].c_str();
            }
            if (ImGui::ListBox("Mesh", &picked_dir2, dir2_items, IM_ARRAYSIZE(dir2_items), 30)){
                folder = dirs2[picked_dir2];
                input_tets = folder + "/tetra.mesh";
                std::cout << input_tets << std::endl; 

                readDotMeshTet(input_tets, TV, TT);
                key_down(viewer,'5',0);
                /*try {
                    std::cout << "Reading "<<file_names[item_current]<<"..." << std::endl;
                    visualizeHexex(file_names[item_current]);
                }
                catch (...) {
                    std::cout << "Reading error !!" << std::endl;
                }*/
            }



            if (ImGui::Button("Show init", ImVec2(-1, 0))){
                viewer.data().clear();
                igl::readOBJ(folder + "/boundary.obj", V, F);
                viewer.data().set_mesh(V, F);
                Eigen::MatrixXd N;
                igl::per_face_normals(V, F, N);
                viewer.data().set_normals(-N);
                Eigen::VectorXi labeling = openFlagging(folder + "/labeling_init.txt", F.rows());
                viewer.data().set_colors(colorsFromFlagging(labeling));

                Eigen::VectorXi axes(F.rows());
                for (int i=0; i<F.rows(); i++) axes(i) = flagToAxis(labeling(i));

                for (int axis=0; axis<3; axis++){
                    int n_f_axis = (axes.array() == axis).count();
                    Eigen::MatrixXd V_ax(3 * n_f_axis, 3); // no need to dupl?
                    Eigen::MatrixXi F_ax(n_f_axis, 3);
                    int next_f = 0;
                    for (int i=0; i<F.rows(); i++){
                        if (axis != axes(i)) continue;
                        V_ax.row(3 * next_f) = V.row(F(i, 0));
                        V_ax.row(3 * next_f + 1) = V.row(F(i, 1));
                        V_ax.row(3 * next_f + 2) = V.row(F(i, 2));
                        F_ax.row(next_f) = Eigen::RowVector3i(3 * next_f, 3 * next_f + 1, 3 * next_f + 2);
                        next_f ++;
                    }
                    if (next_f != n_f_axis) std::cout << "UNEXPECTED" << std::endl;

                    std::string file = "../todel_ax" + std::to_string(axis) + ".obj";
                    igl::writeOBJ(file, V_ax, F_ax);
                }
            }


            if (ImGui::Button("Show polycube", ImVec2(-1, 0))){
                viewer.data().clear();

                Eigen::MatrixXd V, V_poly;
                Eigen::MatrixXi F;

                std::string model_name = dirs2_names[picked_dir2];
                std::string dataset_name = dirs1_names[picked_dir1];

                if (dataset_name == "abc"){
                    igl::readOBJ(folder + "/boundary.obj", V, F);
                    igl::readOBJ(folder + "/polycube_surf_int.obj", V_poly, F); // reference polycube here
                }
                else {
                    std::string results_path = "../../RESULTS_Evocube/"; // TODO CHANGE THIS
                    std::string remesh_path = results_path + "/" + dataset_name + "/" + model_name + "_remesh.mesh";
                    std::string polycube_path = results_path + "/" + dataset_name + "/" + model_name + "_final_polycube.mesh";

                    std::cout << "Looking for final results here:" << std::endl;
                    std::cout << remesh_path << std::endl;
                    std::cout << polycube_path << std::endl;

                    Eigen::MatrixXi F1, F2;
                    Eigen::MatrixXd V1, V2, V_tets1, V_tets2; // V1 reference triangle mesh, V2 deformed

                    // For tet meshes comparison:
                    Eigen::MatrixXi tets1, tets2;
                    readDotMeshTet(remesh_path, V_tets1, tets1);
                    readDotMeshTet(polycube_path, V_tets2, tets2);
                    tetToBnd(V_tets1, tets1, V1, F1);
                    tetToBnd(V_tets2, tets2, V2, F2);

                    V = V1;
                    V_poly = V2;
                    F = F1;
                }


                Eigen::MatrixXd N_poly, N;
                igl::per_face_normals(V_poly, F, N_poly);
                igl::per_face_normals(V, F, N);

                Eigen::MatrixXd V_dupl(F.rows() * 3, 3);
                Eigen::MatrixXi F_dupl(F.rows(), 3);
                for (int f = 0; f<F.rows(); f++) {
                    F_dupl.row(f) << 3 * f, 3 * f + 1, 3 * f + 2;
                    V_dupl.row(3 * f) = V.row(F(f, 0));
                    V_dupl.row(3 * f + 1) = V.row(F(f, 1));
                    V_dupl.row(3 * f + 2) = V.row(F(f, 2));
                }

                Eigen::MatrixXd U(3 * F.rows(), 2);
                Eigen::VectorXi axes(F.rows());
                for (int f = 0; f<F.rows(); f++) {
                    Eigen::RowVector3d n = N_poly.row(f);
                    if (std::abs(std::abs(n(0)) - 1) < 1E-8) axes(f) = 0;
                    else if (std::abs(std::abs(n(1)) - 1) < 1E-8) axes(f) = 1;
                    else if (std::abs(std::abs(n(2)) - 1) < 1E-8) axes(f) = 2;
                    
                    for (int fc=0; fc<3; fc++){
                        for (int d=0; d<2; d++) {
                            double u = V_poly(F(f, fc), (d + axes(f) + 1) % 3);
                            //u /= 2;
                            U(3*f+fc, d) = u;
                        }
                    }
                }
                viewer.data().set_mesh(V_dupl, F_dupl);
                viewer.data().set_uv(U);
                viewer.data().show_lines = false;
                viewer.data().set_normals(-N);
                viewer.data().show_texture = true;
                Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
                igl::png::readPNG("../data/grey_checkerboard.png",R,G,B,A);
                viewer.data().set_texture(R,G,B,A);

                Eigen::VectorXi labeling = openFlagging(folder + "/labeling.txt", F.rows());
                Eigen::MatrixXd colors = colorsFromFlagging(labeling);
                viewer.data().set_colors(colors);

                for (int axis=0; axis<3; axis++){
                    int n_f_axis = (axes.array() == axis).count();
                    Eigen::MatrixXd V_ax(3 * n_f_axis, 3); // no need to dupl?
                    Eigen::MatrixXd V_poly_ax(3 * n_f_axis, 3); // no need to dupl?
                    Eigen::MatrixXi F_ax(n_f_axis, 3);
                    Eigen::MatrixXd CN_ax(n_f_axis, 3);
                    Eigen::MatrixXi FN_ax(n_f_axis, 3);
                    Eigen::MatrixXd U_ax(3 * n_f_axis, 2); // per-vertex texture coord
                    int next_f = 0;
                    for (int i=0; i<F.rows(); i++){
                        if (axis != axes(i)) continue;
                        V_ax.row(3 * next_f) = V.row(F(i, 0));
                        V_ax.row(3 * next_f + 1) = V.row(F(i, 1));
                        V_ax.row(3 * next_f + 2) = V.row(F(i, 2));
                        V_poly_ax.row(3 * next_f) = V_poly.row(F(i, 0));
                        V_poly_ax.row(3 * next_f + 1) = V_poly.row(F(i, 1));
                        V_poly_ax.row(3 * next_f + 2) = V_poly.row(F(i, 2));
                        F_ax.row(next_f) = Eigen::RowVector3i(3 * next_f, 3 * next_f + 1, 3 * next_f + 2);
                        CN_ax.row(next_f) = N.row(i);
                        FN_ax.row(next_f) = Eigen::RowVector3i(next_f, next_f, next_f); // replace with vertex normals here?
                        for (int fc=0; fc<3; fc++){
                            for (int d=0; d<2; d++) {
                                double u = V_poly(F(i, fc), (d + axes(i) + 1) % 3);
                                u /= 1;
                                U_ax(3 * next_f + fc, d) = u;
                            }
                        }
                        next_f ++;
                    }
                    if (next_f != n_f_axis) std::cout << "UNEXPECTED" << std::endl;

                    std::string labeling_file = "../labeling_ax" + std::to_string(axis) + ".obj";
                    std::string polycube_file = "../polycube_ax" + std::to_string(axis) + ".obj";

                    // writeOBJ DOC:
                    //   str  path to outputfile
                    //   V  #V by 3 mesh vertex positions
                    //   F  #F by 3|4 mesh indices into V
                    //   CN #CN by 3 normal vectors
                    //   FN  #F by 3|4 corner normal indices into CN
                    //   TC  #TC by 2|3 texture coordinates
                    //   FTC #F by 3|4 corner texture coord indices into TC
                    igl::writeOBJ(labeling_file, V_ax, F_ax, CN_ax, FN_ax, U_ax, F_ax);
                    igl::writeOBJ(polycube_file, V_poly_ax, F_ax, CN_ax, FN_ax, U_ax, F_ax);
                }
                
            }


            ImGui::End();
        }
    
    };

    //viewer.core().lighting_factor = 0.7;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(255.0/255.0, 255.0/255.0, 255.0/255.0, 1.0);
    viewer.launch();
}
