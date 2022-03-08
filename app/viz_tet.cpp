#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerCore.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOFF.h>
#include <filesystem>
#include <imgui.h>

#include "mesh_io.h"
#include "flagging_utils.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXd tet_colors;

bool display_labeling = false;

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
            viewer.data().set_colors(tet_colors);
            std::cout << "tet_colors.rows(): " << tet_colors.rows() << std::endl;
        }
    }

    return false;
}

std::vector<std::string> listDirs(const std::string& s){
    std::vector<std::string> r;
    for(auto& p : std::filesystem::directory_iterator(s))
        if (p.is_directory())
            r.push_back(p.path().string());

    return r;
}

int main(int argc, char *argv[]){
    using namespace Eigen;
    using namespace std;

    std::string folder = "../data/DATASET/";
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


    if (display_labeling){
        Eigen::VectorXi tet_labeling = openFlagging(folder + "labeling_on_tets.txt", 4 * TT.rows());
        tet_colors = colorsFromFlagging(tet_labeling);
    }

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


            ImGui::End();
        }
    
    };

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}
