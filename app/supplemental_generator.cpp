#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "latex.h"
#include "logging.h"

#define DEFAULT_INPUT_PATH      "../data/DATASET_12_Avril/basic_mambo/"
#define INPUT_COVER_PAGE        "../supplemental/cover_page.tex"
#define INPUT_POLYCUBE_TAGNAME  "/FastPolycubeFloat" //which polycube distortions (from the the log file) to insert. "/FastPolycubeFloat" or "/FastPolycubeInt"
#define DEFAULT_OUTPUT_PATH    "../supplemental/"

//usage : ./supplemental_generator [output] [input1] [input2 input3 ...]
//  
//  output is where the files will be written. Default is DEFAULT_OUTPUT_PATH
//  input1 is the path to the fist dataset to include. default is DEFAULT_INPUT_PATH
//  other datasets can be included to assemble an unique supplemental document
//
//  ex : ./supplemental_generator ../paper/supplemental ../data/first_dataset_CAD ../data/first_dataset_Smooth ../data/second_dataset ../another_dataset

int main(int argc, char** argv) {

    // get arguments

    std::string output_path = DEFAULT_OUTPUT_PATH;
    std::vector<std::string> input_paths;

    if(argc >= 2)
        output_path = argv[1];

    for(int arg_parser = 2; arg_parser < argc; arg_parser++) {
        input_paths.push_back(argv[arg_parser]);
    }
    if(input_paths.empty())
        input_paths.push_back(DEFAULT_INPUT_PATH);

    // create output file

    if(!std::filesystem::exists(output_path))
        std::filesystem::create_directory(output_path);

    std::vector<std::string> no_logs_meshes, invalid_labeling_meshes, missing_figs_meshes;

    LatexDoc latex(std::string(output_path) + "/supplemental.tex");
    latex.add_subpage(std::filesystem::relative(INPUT_COVER_PAGE,output_path));//argument = where is the cover page relative to the output folder

    // parse input datasets

    int success_count = 0, total_meshes_number = 0;

    for(const std::string input_path: input_paths) { //for each given input dataset
        std::string dataset_short_name = std::filesystem::path(input_path).parent_path().filename();
        coloredPrint("~~~~~~ DATASET " + dataset_short_name + " ~~~~~~","cyan");

        std::set<std::filesystem::directory_entry> entries_sorted_by_name;
        for(const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(input_path))//interators have no order -> sort by alphabetical order with a set
            entries_sorted_by_name.insert(dir_entry);

        for(const std::filesystem::directory_entry& dir_entry : entries_sorted_by_name) {
            if(!dir_entry.is_directory()) continue;
            total_meshes_number++;
            std::string current_filepath = dir_entry.path().string();
            std::cout << "Working on " << current_filepath << " : ";
            int status_code = latex.add_mesh(dir_entry.path(), INPUT_POLYCUBE_TAGNAME);
            if (status_code == 1) {
                coloredPrint("Some/all figures are missing","red");
                missing_figs_meshes.push_back(current_filepath);
            }
            else if (status_code == 2) {
                coloredPrint("The final labeling is invalid","red");
                invalid_labeling_meshes.push_back(current_filepath);
            }
            else if (status_code == 3) {
                coloredPrint("No log file, skipped","red");
                no_logs_meshes.push_back(current_filepath);
            }
            else if (status_code != 0) {
                coloredPrint(std::string("Unknown status code of LatexDoc::add_mesh() : ") + std::to_string(status_code) + " with " + current_filepath, "red");
                return 1;
            }
            else {
                std::cout << "Done" << std::endl;
                success_count ++;
            }
        }
    }

    LatexDoc dummy_timeplot(std::string(output_path) + "/time_plot.tex");
    std::array<double,8> cpu = {20,1.4,72,6,20,1.4,72,6};
    std::array<double,8> real = {10,2,20,10,10,2,20,10};
    dummy_timeplot.add_time_plot(cpu,real);

    std::cout << std::endl << "-- SUMMARY ----------------" << std::endl;
    if(!missing_figs_meshes.empty()) {
        coloredPrint(std::to_string(missing_figs_meshes.size()) + " mesh(es) have some figures missing:", "red");
        for(auto& name : missing_figs_meshes)
            std::cout << "\t" << name << std::endl;
    }
    if(!invalid_labeling_meshes.empty()) {
        coloredPrint(std::to_string(invalid_labeling_meshes.size()) + " mesh(es) have an invalid final labeling:", "red");
        for(auto& name : invalid_labeling_meshes)
            std::cout << "\t" << name << std::endl;
    }
    if(!no_logs_meshes.empty()) {
        coloredPrint(std::to_string(no_logs_meshes.size()) + " mesh(es) don't have a log file:", "red");
        for(auto& name : no_logs_meshes)
            std::cout << "\t" << name << std::endl;
    }
    if(missing_figs_meshes.empty() && invalid_labeling_meshes.empty() && no_logs_meshes.empty())
        coloredPrint("No figure missing", "green");
    coloredPrint("Success: " + std::to_string(success_count) + " /"
                             + std::to_string(total_meshes_number), "green");
    std::cout << "---------------------------" << std::endl;
    return 0;
}

