#include <iostream>
#include <fstream>
#include <string>
#include <utility> //std::pair
#include <vector>

class LatexDoc {

public:
    LatexDoc(std::string filename);
    ~LatexDoc();

    void add_subpage(std::filesystem::path path_to_subpage);
    bool add_mesh(std::filesystem::path path_to_mesh_folder);

private:
    std::ofstream ofs;

    bool add_pictures(std::filesystem::path path_to_mesh_folder, int figure_id, std::string caption);
    void add_table(const std::vector<std::pair<std::string,std::string>>& values);

};