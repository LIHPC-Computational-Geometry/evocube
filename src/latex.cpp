#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> //std::replace
#include <utility> //std::pair
#include <vector>

#include "latex.h"
#include "logging.h"

#define PATH_TO_FIGURE(mesh_folder,fig_type,fig_number) ((mesh_folder) + "fig" + std::to_string(fig_type) + "_" + std::to_string(fig_number) + ".png")

LatexDoc::LatexDoc(std::string filename)
    : ofs(filename)
{
    ofs << "\\documentclass{article}%" << std::endl;
    ofs << "\\usepackage[T1]{fontenc}%" << std::endl;
    ofs << "\\usepackage[utf8]{inputenc}%" << std::endl;
    ofs << "\\usepackage{lmodern}%" << std::endl;
    ofs << "\\usepackage{textcomp}%" << std::endl;
    ofs << "\\usepackage{lastpage}%" << std::endl;
    ofs << "\\usepackage{geometry}%" << std::endl;
    ofs << "\\geometry{rmargin=0.5cm,lmargin=0.5cm,tmargin=1cm,bmargin=2cm}%" << std::endl;
    ofs << "\\usepackage{hyperref}%" << std::endl;
    ofs << "\\usepackage{subcaption}%" << std::endl;
    ofs << "\\usepackage{graphicx}%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "\\begin{document}%" << std::endl;
}

LatexDoc::~LatexDoc() {
    ofs << "\\end{document}" << std::endl;
    ofs.close();
}

void LatexDoc::add_subpage(std::filesystem::path path_to_subpage) {
    ofs << "\\input{" << path_to_subpage.string() << "}%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
}

bool LatexDoc::add_mesh(std::filesystem::path path_to_mesh_folder) {

    bool incomplete_page = false;
    std::string mesh_name = path_to_mesh_folder.filename(), label = "";

    //replace '_' by "\\_" for mesh_name
    //remove '_' for label
    size_t pos = 0;
    while(pos < mesh_name.length()) {
        if(mesh_name[pos] == '_') {
            mesh_name.replace(pos, 1, "\\_");
            pos += 2;
        }
        else {
            label += mesh_name[pos];
            pos++;
        }
    }
    
    ofs << "\\section{" << mesh_name << "}%" << std::endl;
    ofs << "\\label{sec:" << label << "}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    //TODO read metrics and put them in a table
    // add_metrics_table(1,0,0,0,2,0.052177);

    incomplete_page |= add_pictures(path_to_mesh_folder,1,mesh_name + ", input");

    ofs << "\\par%" << std::endl;
    ofs << "}%" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    incomplete_page |= add_pictures(path_to_mesh_folder,4,mesh_name + ", labelling");

    ofs << "\\par%" << std::endl;
    ofs << "}%" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    incomplete_page |= add_pictures(path_to_mesh_folder,0,mesh_name + ", polycube");

    ofs << "\\par%" << std::endl;
    ofs << "}%" << std::endl;
    ofs << "\\clearpage%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;

    if(incomplete_page)
        coloredPrint(std::string("Some/all figures of ") + mesh_name + " are missing","red");
    return incomplete_page;
}

bool LatexDoc::add_pictures(std::filesystem::path path_to_mesh_folder, int figure_id, std::string caption) {

    std::string fig0 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,0),
                fig1 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,0),
                fig2 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,0),
                fig3 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,0);

    if( std::filesystem::exists(fig0) && 
        std::filesystem::exists(fig1) && 
        std::filesystem::exists(fig2) && 
        std::filesystem::exists(fig3)
    ) {
        ofs << "\\begin{figure}[ht]%" << std::endl;

        ofs << "\\begin{subfigure}{0.25\\linewidth}%" << std::endl;
        ofs << "\\includegraphics[width=\\linewidth]{" + fig0 + "}%" << std::endl;
        ofs << "\\end{subfigure}%" << std::endl;

        ofs << "\\begin{subfigure}{0.25\\linewidth}%" << std::endl;
        ofs << "\\includegraphics[width=\\linewidth]{" + fig1 + "}%" << std::endl;
        ofs << "\\end{subfigure}%" << std::endl;

        ofs << "\\begin{subfigure}{0.25\\linewidth}%" << std::endl;
        ofs << "\\includegraphics[width=\\linewidth]{" + fig2 + "}%" << std::endl;
        ofs << "\\end{subfigure}%" << std::endl;

        ofs << "\\begin{subfigure}{0.25\\linewidth}%" << std::endl;
        ofs << "\\includegraphics[width=\\linewidth]{" + fig3 + "}%" << std::endl;
        ofs << "\\end{subfigure}%" << std::endl;

        ofs << "\\caption{" << caption << "}%" << std::endl;
        ofs << "\\end{figure}%" << std::endl;
        return false;
    }
    else 
        return true;
}

void LatexDoc::add_table(const std::vector<std::pair<std::string,std::string>>& values) {

    if(values.empty())
        return;

    //declare a tabular environment with values.size() centered columns
    ofs << "\\begin{tabular}{|";
    for(int index = 0; index < values.size(); index++) { ofs << "c|"; }
    ofs << "}%" << std::endl;
    ofs << "\\hline%" << std::endl;

    //write column names
    for(int index = 0; index < values.size(); index++) {
        ofs << values[index].first;
        if(index != values.size()-1)
            ofs << " & ";
        else
            ofs << " \\\\%" << std::endl;
    }
    ofs << "\\hline%" << std::endl;

    //write column values
    for(int index = 0; index < values.size(); index++) {
        ofs << values[index].second;
        if(index != values.size()-1)
            ofs << " & ";
        else
            ofs << " \\\\%" << std::endl;
    }
    ofs << "\\hline%" << std::endl;
    ofs << "\\end{tabular}%" << std::endl;
}

void LatexDoc::add_metrics_table(bool found_valid_labelling, int invalid_patches, int invalid_corners, int invalid_boundaries, int nb_turning_points, double fidelity) {
    std::vector<std::pair<std::string,std::string>> table_content;
    table_content.push_back(std::make_pair<std::string,std::string>("Found valid labelling",std::to_string(found_valid_labelling)));
    table_content.push_back(std::make_pair<std::string,std::string>("Inv. Patches",         std::to_string(invalid_patches)));
    table_content.push_back(std::make_pair<std::string,std::string>("Inv. Corners",         std::to_string(invalid_corners)));
    table_content.push_back(std::make_pair<std::string,std::string>("Inv. Boundaries",      std::to_string(invalid_boundaries)));
    table_content.push_back(std::make_pair<std::string,std::string>("\\# turning points",     std::to_string(nb_turning_points)));
    table_content.push_back(std::make_pair<std::string,std::string>("Fidelity",             std::to_string(fidelity)));
    add_table(table_content);
}