#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> //std::replace, std::max
#include <vector>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <sstream>
#include <array>

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
    ofs << "% horizontal stacked bar chart" << std::endl;
    ofs << "\\usepackage{xcolor}%" << std::endl;
    ofs << "\\usepackage{tikz}%" << std::endl;
    ofs << "\\usepackage{pgfplots}%" << std::endl;
    ofs << "\\pgfplotsset{compat=1.17}%" << std::endl;
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

int LatexDoc::add_mesh(std::filesystem::path path_to_mesh_folder, std::string polycube_tagname, time_plot_entry& timings) {

    bool figures_are_missing = false;
    std::string mesh_name = path_to_mesh_folder.filename(), section_name = "", label = "";

    std::ifstream logs_path = path_to_mesh_folder / "logs.json";
    if(!logs_path.good())
        return 3;
    nlohmann::json j;
    logs_path >> j;

    if (j.empty()) 
        return 3;

#ifdef REMOVE_INPUT_TRI_SUFFIX
    if(mesh_name.size() > 10) {
        //if ends with "_input_tri", remove this suffix
        if( mesh_name.substr(mesh_name.size()-10,10) == "_input_tri" ) {
            mesh_name.resize(mesh_name.size()-10);
        }
    }
#endif

    section_name = escape_special_chars(mesh_name);
    
    ofs << "\\section{" << section_name << "}%" << std::endl;
    ofs << "{%" << std::endl;
    ofs << "\\small%" << std::endl;
    ofs << "\\centering%" << std::endl;

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
    figures_are_missing |= add_pictures(path_to_mesh_folder,1,section_name + ", input");

    // TABLES ABOUT THE LABELING

    if(j.value<std::string>("/LabelingFinal/InvalidTotal"_json_pointer,"")!="0") {
        ofs << "\\par" << std::endl;
        ofs << "\\textit{Failed to find a valid labeling}" << std::endl;
        //show the invalid labeling anyway, but not the polycube
        ofs << "\\par" << std::endl;
        figures_are_missing |= add_pictures(path_to_mesh_folder,4,section_name + ", labelling");
        ofs << "}%" << std::endl;
        ofs << "\\clearpage%" << std::endl;
        ofs << "%" << std::endl;
        ofs << "%" << std::endl;
        ofs << "%" << std::endl;
        return 2;
    }

    ofs << "\\par%" << std::endl;
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

    //update the timing counters
    timings.pre_optimization             += std::stod(j.value<std::string>("/Timing/PreGenetics"_json_pointer,"NAN"));
    timings.individual_selection         += std::stod(j.value<std::string>("/Timing/CreateIndiv"_json_pointer,"NAN"));
    timings.individual_mutations         += std::stod(j.value<std::string>("/Timing/Mutations"_json_pointer,"NAN"));
    timings.charts_and_turning_points    += std::stod(j.value<std::string>("/Timing/ChartsAndTps"_json_pointer,"NAN"));
    timings.fitness_evaluation           += std::stod(j.value<std::string>("/Timing/Eval"_json_pointer,"NAN"));
    timings.crossing                     += std::stod(j.value<std::string>("/Timing/Cross"_json_pointer,"NAN"));
    timings.insertion_in_archive         += std::stod(j.value<std::string>("/Timing/Archive"_json_pointer,"NAN"));
    timings.post_optimization            += std::stod(j.value<std::string>("/Timing/PostGenetics"_json_pointer,"NAN"));
    timings.genetics                     += std::stod(j.value<std::string>("/Timing/Genetics"_json_pointer,"NAN"));

    // LABELLING FIGURE

    ofs << "\\par" << std::endl;
    figures_are_missing |= add_pictures(path_to_mesh_folder,4,section_name + ", labelling");

    // TABLE ABOUT THE POLYCUBE

    ofs << "\\par%" << std::endl;
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
    figures_are_missing |= add_pictures(path_to_mesh_folder,0,section_name + ", polycube");

    ofs << "\\par%" << std::endl;
    ofs << "}%" << std::endl;
    ofs << "\\clearpage%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;
    ofs << "%" << std::endl;

    return figures_are_missing;// 0 = good, 1 = some figs missing
}

void LatexDoc::add_time_plot(const time_plot_entry& cpu, const time_plot_entry& real) {
    ofs << "\\definecolor{preliminary}{HTML}{dfff00}%" << std::endl;
    ofs << "\\definecolor{createindiv}{HTML}{ffbf00}%" << std::endl;
    ofs << "\\definecolor{mutations}{HTML}{ff7f50}%" << std::endl;
    ofs << "\\definecolor{chartscolor}{HTML}{de3163}%" << std::endl;
    ofs << "\\definecolor{fitnesscolor}{HTML}{9fe2bf}%" << std::endl;
    ofs << "\\definecolor{crossingcolor}{HTML}{40e0d0}%" << std::endl;
    ofs << "\\definecolor{archivecolor}{HTML}{6495ed}%" << std::endl;
    ofs << "\\definecolor{postcolor}{HTML}{ccccff}%" << std::endl;

    ofs << "\\begin{figure}%" << std::endl;
    ofs << "\\centering%" << std::endl;
    ofs << "\\begin{tikzpicture}%" << std::endl;
    ofs << "\\begin{axis}[%" << std::endl;
    ofs << "    xbar stacked,%" << std::endl;
    ofs << "    legend style={%" << std::endl;
    ofs << "        legend columns=2,%" << std::endl;
    ofs << "        at={(xticklabel cs:0.5)},%" << std::endl;
    ofs << "        anchor=north,%" << std::endl;
    ofs << "        draw=none,%" << std::endl;
    ofs << "        cells={anchor=west}, % left-align cell content" << std::endl;
    ofs << "        align=left%" << std::endl;
    ofs << "    },%" << std::endl;
    ofs << "    ytick=data,%" << std::endl;
    ofs << "    axis y line*=none, % line here for a vertical bar at 0" << std::endl;
    ofs << "    axis x line*=bottom,%" << std::endl;
    ofs << "    tick label style={font=\\footnotesize},%" << std::endl;
    ofs << "    legend style={font=\\footnotesize},%" << std::endl;
    ofs << "    label style={font=\\footnotesize},%" << std::endl;
    ofs << "    xtick={0,1000,2000,3000,4000},%" << std::endl;
    ofs << "    width=0.8\\linewidth,%" << std::endl;
    ofs << "    bar width=6mm,%" << std::endl;
    ofs << "    xlabel={Time in hours},%" << std::endl;
    ofs << "    yticklabels={CPU time, Real time},%" << std::endl;
    ofs << "    xmin=0,%" << std::endl;
    ofs << "    xmax=" << std::max(cpu.sum(),real.sum())*1.1 << ",%" << std::endl;
    ofs << "    area legend,%" << std::endl;
    ofs << "    y=8mm,%" << std::endl;
    ofs << "    enlarge y limits={abs=0.625},%" << std::endl;
    ofs << "]%" << std::endl;

    ofs << "\\addplot[preliminary,fill=preliminary] coordinates%" << std::endl;
    ofs << "{(" << cpu.pre_optimization << ",0) (" << real.pre_optimization << ",1)};%" << std::endl;
    ofs << "\\addplot[createindiv,fill=createindiv] coordinates%" << std::endl;
    ofs << "{(" << cpu.individual_selection << ",0) (" << real.individual_selection << ",1)};%" << std::endl;
    ofs << "\\addplot[mutations,fill=mutations] coordinates%" << std::endl;
    ofs << "{(" << cpu.individual_mutations << ",0) (" << real.individual_mutations << ",1)};%" << std::endl;
    ofs << "\\addplot[chartscolor,fill=chartscolor] coordinates%" << std::endl;
    ofs << "{(" << cpu.charts_and_turning_points << ",0) (" << real.charts_and_turning_points << ",1)};%" << std::endl;
    ofs << "\\addplot[fitnesscolor,fill=fitnesscolor] coordinates%" << std::endl;
    ofs << "{(" << cpu.fitness_evaluation << ",0) (" << real.fitness_evaluation << ",1)};%" << std::endl;
    ofs << "\\addplot[crossingcolor,fill=crossingcolor] coordinates%" << std::endl;
    ofs << "{(" << cpu.crossing << ",0) (" << real.crossing << ",1)};%" << std::endl;
    ofs << "\\addplot[archivecolor,fill=archivecolor] coordinates%" << std::endl;
    ofs << "{(" << cpu.insertion_in_archive << ",0) (" << real.insertion_in_archive << ",1)};%" << std::endl;
    ofs << "\\addplot[postcolor,fill=postcolor] coordinates%" << std::endl;
    ofs << "{(" << cpu.post_optimization << ",0) (" << real.post_optimization << ",1)};%" << std::endl;
    ofs << "\\legend{Pre-optimization, Individual selection,  Individual mutations, Charts and turning points, Fitness evaluation, Crossing, Insertion in archive, Post-optimization}%" << std::endl;
    ofs << "\\end{axis}  %" << std::endl;
    ofs << "\\end{tikzpicture}%" << std::endl;
    ofs << "\\caption{Time plot of our labeling optimization over all tested models}%" << std::endl;
    ofs << "\\end{figure}%" << std::endl;
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
            ofs << escape_special_chars(values[row].at(col));
            if(col != values[0].size()-1)
                ofs << " & ";
            else
                ofs << " \\\\%" << std::endl;
        }
        ofs << "\\hline%" << std::endl;
    }
    ofs << "\\end{tabular}%" << std::endl;
}

std::string escape_special_chars(const std::string input) {
    std::string output;
    for(int pos = 0; pos < input.length(); pos++) {
        if( (input[pos] == '#') ||
            (input[pos] == '_') ) { output += "\\"; }
        output += input[pos];
    }
    return output;
}

std::string remove_special_chars(const std::string input) {
    std::string output;
    for(int pos = 0; pos < input.length(); pos++) {
        if( (input[pos] != '#') && 
            (input[pos] != '_') ) { output += input[pos]; }
    }
    return output;
}

std::string double2string(double value, int precision) {
    std::stringstream rounded;
    rounded << std::fixed << std::setprecision(precision) << value;
    return rounded.str();
}