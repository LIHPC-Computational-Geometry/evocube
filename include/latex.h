#include <iostream>
#include <fstream>
#include <string>

class LatexDoc {

public:
    LatexDoc(std::string filename);
    ~LatexDoc();

    bool add_mesh(std::filesystem::path path_to_mesh_folder);

private:
    std::ofstream ofs;

    bool add_pictures(std::filesystem::path path_to_mesh_folder, int figure_id, std::string caption);
};