#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> //std::replace
#include <utility> //std::pair
#include <vector>
#include <nlohmann/json.hpp>

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
    
    //write the logs in a table. put an empty string if the key doesn't exist.
    //to get all the items in LabelingFinal
    // for (auto& el : j["/LabelingFinal"_json_pointer].items() ) {
    //     std::cout << el.key() << " " << el.value() << std::endl;
    // }
    std::vector<std::pair<std::string,std::string>> table1, table2, table3;
    table1.push_back(std::make_pair("#charts",            j.value<std::string>("/LabelingFinal/#charts"_json_pointer,         "")));
    table1.push_back(std::make_pair("#corners",           j.value<std::string>("/LabelingFinal/#corners"_json_pointer,        "")));
    table1.push_back(std::make_pair("#tps",               j.value<std::string>("/LabelingFinal/#tps"_json_pointer,            "")));
    table1.push_back(std::make_pair("fidelity",           j.value<std::string>("/LabelingFinal/fidelity"_json_pointer,        "")));
    table1.push_back(std::make_pair("score",              j.value<std::string>("/LabelingFinal/ScoreFinal"_json_pointer,      "")));
    table1.push_back(std::make_pair("angle/area dist.",   j.value<std::string>(nlohmann::json::json_pointer(polycube_tagname+"/AngleDistortion"), "") + "/" +
                                                                j.value<std::string>(nlohmann::json::json_pointer(polycube_tagname+"/AreaDistortion"),  "")));
    table1.push_back(std::make_pair("stretch",            j.value<std::string>(nlohmann::json::json_pointer(polycube_tagname+"/Stretch"),         "")));
    add_table(table1);

    ofs << "\\par" << std::endl;
    ofs << "\\vspace{10pt}%" << std::endl;

    table2.push_back(std::make_pair("GraphCutCompactness",j.value<std::string>("/GraphCutParams/CompactCoeff"_json_pointer,   "")));
    table2.push_back(std::make_pair("GraphCutFidelity",   j.value<std::string>("/GraphCutParams/FidelityCoeff"_json_pointer,  "")));
    table2.push_back(std::make_pair("GraphCutInvalidies", j.value<std::string>("/LabelingGraphCut/InvalidTotal"_json_pointer, "")));
    table2.push_back(std::make_pair("GraphCutScore",      j.value<std::string>("/LabelingGraphCut/ScoreFinal"_json_pointer,   "")));
    add_table(table2);

    ofs << "\\par" << std::endl;
    ofs << "\\vspace{10pt}%" << std::endl;

    double total_time = 0.0;
    total_time += std::stod(j.value<std::string>("/Timing/PreGenetics"_json_pointer,   "NAN"));
    total_time += std::stod(j.value<std::string>("/Timing/Genetics"_json_pointer,      "NAN"));
    total_time += std::stod(j.value<std::string>("/Timing/PostGenetics"_json_pointer,  "NAN"));
    table3.push_back(std::make_pair("TotalTime (s)",      std::to_string(total_time)));
    table3.push_back(std::make_pair("#generations",       j.value<std::string>("/#generations"_json_pointer,"")));
    add_table(table3);

    figures_are_missing |= add_pictures(path_to_mesh_folder,1,section_name + ", input");

    ofs << "\\par%" << std::endl;
    ofs << "}%" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

    figures_are_missing |= add_pictures(path_to_mesh_folder,4,section_name + ", labelling");

    ofs << "\\par%" << std::endl;
    ofs << "}%" << std::endl;
    ofs << "\\vspace{-20pt}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\centering%" << std::endl;

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

        std::string column_name = values[index].first;
        //replace '#' with '\\#'
        int pos=0;
        while(pos < column_name.length()) {
            if(column_name[pos] == '#') {
                column_name.replace(pos, 1, "\\#");
                pos += 2;
            }
            else
                pos++;
        }

        ofs << column_name;
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