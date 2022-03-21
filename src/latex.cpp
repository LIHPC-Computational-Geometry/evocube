#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>

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

    // cover page

    ofs << "\\normalsize%" << std::endl;
    ofs << "\\title{Evocube {-} Supplemental material}%" << std::endl;
    ofs << "\\date{}%" << std::endl;
    ofs << "\\maketitle%" << std::endl;
    ofs << "\\large%" << std::endl;
    ofs << "This file includes all the supplemental material submitted along with \\textit{Evocube: A Genetic Labelling Framework for Polycube-Maps}.%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "\\vspace{20pt}%" << std::endl;
    ofs << "\\normalsize%" << std::endl;
    ofs << "\\clearpage%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
}

LatexDoc::~LatexDoc() {
    ofs << "\\end{document}" << std::endl;
    ofs.close();
}

bool LatexDoc::add_mesh(std::filesystem::path path_to_mesh_folder) {

    bool incomplete_page = false;
    std::string mesh_name = path_to_mesh_folder.parent_path().filename();//extract the folder name from the path
    
    ofs << "\\section{" << mesh_name << "}%" << std::endl;
    ofs << "\\label{sec:" << mesh_name << "}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    //TODO read metrics and put them in a table

    incomplete_page |= add_pictures(path_to_mesh_folder,1,mesh_name + ", input");

    ofs << "\\par%" << std::endl;
    ofs << "}" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    incomplete_page |= add_pictures(path_to_mesh_folder,4,mesh_name + ", labelling");

    ofs << "\\par%" << std::endl;
    ofs << "}" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    incomplete_page |= add_pictures(path_to_mesh_folder,0,mesh_name + ", polycube");

    ofs << "\\par%" << std::endl;
    ofs << "}" << std::endl;
    ofs << "\\clearpage%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;

    if(incomplete_page)
        coloredPrint(std::string("Some/all figures of ") + mesh_name + " are missing","red");
    return incomplete_page;
}

bool LatexDoc::add_pictures(std::filesystem::path path_to_mesh_folder, int figure_id, std::string caption) {

    std::string fig0 = PATH_TO_FIGURE(path_to_mesh_folder.string(),figure_id,0),
                fig1 = PATH_TO_FIGURE(path_to_mesh_folder.string(),figure_id,0),
                fig2 = PATH_TO_FIGURE(path_to_mesh_folder.string(),figure_id,0),
                fig3 = PATH_TO_FIGURE(path_to_mesh_folder.string(),figure_id,0);

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