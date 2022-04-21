#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>

#define REMOVE_INPUT_TRI_SUFFIX

struct time_plot_entry {
    double pre_optimization;
    double individual_selection;
    double individual_mutations;
    double charts_and_turning_points;
    double fitness_evaluation;
    double crossing;
    double insertion_in_archive;
    double post_optimization;
    double genetics;
};
#define INIT_TIME_PLOT_ENTRY(struct) \
    (struct).pre_optimization = 0.0; \
    (struct).individual_selection = 0.0; \
    (struct).individual_mutations = 0.0; \
    (struct).charts_and_turning_points = 0.0; \
    (struct).fitness_evaluation = 0.0; \
    (struct).crossing = 0.0; \
    (struct).insertion_in_archive = 0.0; \
    (struct).post_optimization = 0.0; \
    (struct).genetics = 0.0;

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
     * @param path_to_mesh_folder   Path to the mesh folder containing the pictures
     * @param polycube_tagname      Name of the JSON tag in which the distortion measures will be read. Must start with '/'.
     *                              For now, should be "/FastPolycubeFloat" or "/FastPolycubeInt"
     * @return  0 if good                                   -> mesh added
     *          1 if some pictures are missing              -> mesh added anyway
     *          2 if labeling invalid (given the log file)  -> create a page saying no valid labeling was found
     *          3 if the log file is not found              -> mesh skipped
     */
    int add_mesh(std::filesystem::path path_to_mesh_folder, std::string polycube_tagname, time_plot_entry& timings);

    void add_time_plot(const time_plot_entry& cpu, const time_plot_entry& real);

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
     * @param values List of rows. A row being a list of strings
     *               All sub-vectors should have the same size
     *               If a string contains '#' or '_', they will be escaped
     */
    void add_table(const std::vector<std::vector<std::string>>& values);

};

//replace '#' with '\\#' and '_' with '\\_'
std::string escape_special_chars(const std::string input);
std::string remove_special_chars(const std::string input);

std::string double2string(double value, int precision);