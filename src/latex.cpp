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
    ofs << "\\documentclass{article}" << std::endl;
    ofs << "\\usepackage{pdfpages}" << std::endl;
    ofs << "\\begin{document}" << std::endl;
}

LatexDoc::~LatexDoc() {
    ofs << "\\end{document}" << std::endl;
    ofs.close();
}

bool LatexDoc::add_mesh(std::string path_to_mesh_folder) {

    bool incomplete_page = false;
    
    ofs << "\\section{" << path_to_mesh_folder << "}" << std::endl;
    ofs << "\\centering" << std::endl;

    //TODO read metrics and put them in a table

    //add input mesh pictures : fig1_X.png
    incomplete_page |= add_pictures(path_to_mesh_folder,1);
    ofs << "\\vspace{-20pt}" << std::endl;

    //add labeling pictures : fig4_X.png
    incomplete_page |= add_pictures(path_to_mesh_folder,4);
    ofs << "\\vspace{-20pt}" << std::endl;

    //add polycube pictures : fig0_X.png
    incomplete_page |= add_pictures(path_to_mesh_folder,0);

    ofs << "\\clearpage" << std::endl;
    ofs << "\\pagebreak" << std::endl;

    if(incomplete_page)
        coloredPrint(std::string("Some/all figures of") + path_to_mesh_folder + " are missing","red");
    return incomplete_page;
}

bool LatexDoc::add_pictures(std::string path_to_mesh_folder, int figure_id) {

    std::string fig0 = PATH_TO_FIGURE(path_to_mesh_folder,figure_id,0),
                fig1 = PATH_TO_FIGURE(path_to_mesh_folder,figure_id,0),
                fig2 = PATH_TO_FIGURE(path_to_mesh_folder,figure_id,0),
                fig3 = PATH_TO_FIGURE(path_to_mesh_folder,figure_id,0);

    if( std::filesystem::exists(fig0) && 
        std::filesystem::exists(fig1) && 
        std::filesystem::exists(fig2) && 
        std::filesystem::exists(fig3)
    ) {
        ofs << "\\includegraphics[width=0.24\\linewidth]{" + fig0 + "}" << std::endl;
        ofs << "\\includegraphics[width=0.24\\linewidth]{" + fig1 + "}" << std::endl;
        ofs << "\\includegraphics[width=0.24\\linewidth]{" + fig2 + "}" << std::endl;
        ofs << "\\includegraphics[width=0.24\\linewidth]{" + fig3 + "}" << std::endl;
        return false;
    }
    else 
        return true;
}