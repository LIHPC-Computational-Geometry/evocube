#include <iostream>
#include <fstream>
#include <string>
#include <utility> //std::pair
#include <vector>

class LatexDoc {

public:
    /**
     * @brief Create a LatexDoc object, wrapping a LaTeX file
     * @param filename Path and name of the LaTeX file to create
     */
    LatexDoc(std::string filename);
    ~LatexDoc();

    /**
     * @brief Insert a LaTeX file with the \c \\input command
     * @param path_to_subpage Path to the LaTeX file to insert
     */
    void add_subpage(std::filesystem::path path_to_subpage);

    /**
     * @brief Insert all the figures of a mesh (input, labelling, polycube) on a new page
     * @param path_to_mesh_folder Path to the mesh folder containing the pictures
     * @return False if the page is complete, True if at least 1 figure is missing
     */
    bool add_mesh(std::filesystem::path path_to_mesh_folder);

private:
    std::ofstream ofs;

    /**
     * @brief Insert 4 subfigures of a mesh. \e figure_id select the type of pictures.
     * @param path_to_mesh_folder   Path to the mesh folder containing the pictures
     * @param figure_id             See \c PATH_TO_FIGURE in \c latex.cpp and the \c figure_generator app
     * @param caption               The caption to put bellow the figure
     * @return False if the figure (= 4 subfigures) is complete, True if at least 1 subfigure is missing
     */
    bool add_pictures(std::filesystem::path path_to_mesh_folder, int figure_id, std::string caption);

    /**
     * @brief Insert a 2 by \c values.size() table
     * @param values List of { column name , value }
     */
    void add_table(const std::vector<std::pair<std::string,std::string>>& values);

    /**
     * @brief Insert a table containing metrics of a labelling
     * @param found_valid_labelling If a valid labelling is found
     * @param invalid_patches       The number of invalid patches
     * @param invalid_corners       The number of invalid corners
     * @param invalid_boundaries    The number of invalid boundaries
     * @param nb_turning_points     The number of turning points
     * @param fidelity              The fidelity value
     */
    void add_metrics_table(bool found_valid_labelling, int invalid_patches, int invalid_corners, int invalid_boundaries, int nb_turning_points, double fidelity);

};