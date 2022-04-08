#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> //std::replace
#include <vector>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <sstream>

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

int LatexDoc::add_mesh(std::filesystem::path path_to_mesh_folder, std::string polycube_tagname) {

    bool figures_are_missing = false;
    std::string mesh_name = path_to_mesh_folder.filename(), section_name = "", label = "";

    std::ifstream logs_path = path_to_mesh_folder / "logs.json";
    if(!logs_path.good())
        return 3;
    nlohmann::json j;
    logs_path >> j;

    if (j.empty()) 
        return 3;

    //replace '_' by "\\_" for section_name
    //remove '_' for label
    size_t pos = 0;
    while(pos < mesh_name.length()) {
        if(mesh_name[pos] == '_') {
            section_name += "\\_";
            pos += 2;
        }
        else {
            section_name += mesh_name[pos];
            label += mesh_name[pos];
            pos++;
        }
    }
    
    ofs << "\\section{" << section_name << "}%" << std::endl;
    ofs << "\\label{sec:" << label << "}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    if(j.value<std::string>("/LabelingFinal/InvalidTotal"_json_pointer,"")!="0") {
        ofs << "Failed to find a valid labeling" << std::endl;
        ofs << "}%" << std::endl;
        ofs << "\\clearpage%" << std::endl;
        ofs << "%" << std::endl;
        ofs << "%" << std::endl;
        ofs << "%" << std::endl;
        return 2;
    }

    // TABLE ABOUT THE INPUT MESH

    std::vector<std::vector<std::string>> table;
    table.push_back({"#vertices","#faces","avg. edge length"});
    table.push_back({
        j.value<std::string>("/InputTris/vertices"_json_pointer,""),
        j.value<std::string>("/InputTris/faces"_json_pointer,""),
        j.value<std::string>("/InputTris/AvgEdgeLength"_json_pointer,"")
    });
    add_table(table);

    // INPUT MESH FIGURE

    ofs << "\\par" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    figures_are_missing |= add_pictures(path_to_mesh_folder,1,section_name + ", input");

    // TABLES ABOUT THE LABELING

    ofs << "\\par%" << std::endl;
    ofs << "\\vspace{15pt}%" << std::endl;
    double graphcut_score, final_score;
    graphcut_score = std::stod(j.value<std::string>("/LabelingGraphCut/ScoreFinal"_json_pointer,"NAN"));
    final_score = std::stod(j.value<std::string>("/LabelingFinal/ScoreFinal"_json_pointer,"NAN"));
    table.clear();
    table.push_back({"","#charts","#corners","#TP","#invalidities","fidelity","score"});
    table.push_back({
        "Initial",
        j.value<std::string>("/LabelingGraphCut/#charts"_json_pointer,""),
        j.value<std::string>("/LabelingGraphCut/#corners"_json_pointer,""),
        j.value<std::string>("/LabelingGraphCut/#tps"_json_pointer,""),
        j.value<std::string>("/LabelingGraphCut/InvalidTotal"_json_pointer,""),
        j.value<std::string>("/LabelingGraphCut/fidelity"_json_pointer,""),
        double2string(graphcut_score,3)
    });
    table.push_back({
        "Final",
        j.value<std::string>("/LabelingFinal/#charts"_json_pointer,""),
        j.value<std::string>("/LabelingFinal/#corners"_json_pointer,""),
        j.value<std::string>("/LabelingFinal/#tps"_json_pointer,""),
        j.value<std::string>("/LabelingFinal/InvalidTotal"_json_pointer,""),
        j.value<std::string>("/LabelingFinal/fidelity"_json_pointer,""),
        double2string(final_score,3)
    });
    add_table(table);

    ofs << "\\par" << std::endl;
    ofs << "\\vspace{5pt}%" << std::endl;
    double graphcut_coeff = std::stod(j.value<std::string>("/GraphCutParams/FidelityCoeff"_json_pointer,"NAN")) / 
                            std::stod(j.value<std::string>("/GraphCutParams/CompactCoeff"_json_pointer,"NAN"));
    double total_time = 0.0;
    total_time += std::stod(j.value<std::string>("/Timing/PreGenetics"_json_pointer,"NAN"));
    total_time += std::stod(j.value<std::string>("/Timing/Genetics"_json_pointer,"NAN"));
    total_time += std::stod(j.value<std::string>("/Timing/PostGenetics"_json_pointer,"NAN"));
    table.clear();
    table.push_back({"GraphCut fidelity/compactness","#generations","duration (s)","LabelingSimilarity"});
    table.push_back({
        double2string(graphcut_coeff,3),
        j.value<std::string>("/#generations"_json_pointer,""),
        std::to_string(total_time),
        j.value<std::string>("/LabelingSimilarity"_json_pointer,"")
    });
    add_table(table);

    // LABELLING FIGURE

    ofs << "\\par" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    figures_are_missing |= add_pictures(path_to_mesh_folder,4,section_name + ", labelling");

    // TABLE ABOUT THE POLYCUBE

    ofs << "\\par%" << std::endl;
    ofs << "\\vspace{15pt}%" << std::endl;
    table.clear();
    table.push_back({"angle/area dist.","stretch"});
    table.push_back({
        j.value<std::string>(nlohmann::json::json_pointer(polycube_tagname+"/AngleDistortion"),"") + "/" +
        j.value<std::string>(nlohmann::json::json_pointer(polycube_tagname+"/AreaDistortion"),""),
        j.value<std::string>(nlohmann::json::json_pointer(polycube_tagname+"/Stretch"),"")
    });
    add_table(table);

    // POLYCUBE FIGURE

    ofs << "\\par" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    figures_are_missing |= add_pictures(path_to_mesh_folder,0,section_name + ", polycube");

    ofs << "\\par%" << std::endl;
    ofs << "}%" << std::endl;
    ofs << "\\clearpage%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;

    return figures_are_missing;// 0 = good, 1 = some figs missing
}

bool LatexDoc::add_pictures(std::filesystem::path path_to_mesh_folder, int figure_id, std::string caption) {

    std::string fig0 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,0),
                fig1 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,1),
                fig2 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,2),
                fig3 = PATH_TO_FIGURE(path_to_mesh_folder.string()+"/",figure_id,3);

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

        ofs << "\\vspace{-20pt}%" << std::endl;
        ofs << "\\caption{" << caption << "}%" << std::endl;
        ofs << "\\end{figure}%" << std::endl;
        return false;
    }
    else 
        return true;
}

void LatexDoc::add_table(const std::vector<std::vector<std::string>>& values) {

    if(values.empty())
        return;

    //declare a tabular environment with values[0].size() centered columns
    ofs << "\\begin{tabular}{|";
    for(int col = 0; col < values[0].size(); col++) { ofs << "c|"; }
    ofs << "}%" << std::endl;
    ofs << "\\hline%" << std::endl;

    //write values
    for(int row = 0; row < values.size(); row++) {
        for(int col = 0; col < values[0].size(); col++) {
            ofs << escape_number_sign(values[row].at(col));
            if(col != values[0].size()-1)
                ofs << " & ";
            else
                ofs << " \\\\%" << std::endl;
        }
        ofs << "\\hline%" << std::endl;
    }
    ofs << "\\end{tabular}%" << std::endl;
}

std::string escape_number_sign(std::string input) {
    std::string output;
    int pos=0;
    while(pos < input.length()) {
        if(input[pos] == '#')
            output += "\\#";
        else
            output += input[pos];
        pos++;
    }
    return output;
}

std::string double2string(double value, int precision) {
    std::stringstream rounded;
    rounded << std::fixed << std::setprecision(precision) << value;
    return rounded.str();
}