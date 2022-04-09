#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "latex.h"
#include "logging.h"

#define INPUT_DATA_PATH             "../output_medium_mambo/"
#define INPUT_COVER_PAGE            "../supplemental/cover_page.tex"
#define INPUT_POLYCUBE_TAGNAME      "/FastPolycubeFloat" //which polycube distortions (from the the log file) to insert. "/FastPolycubeFloat" or "/FastPolycubeInt"
#define OUTPUT_SUPPLEMENTAL_PATH    "../supplemental/"

int main(int argc, char** argv) {

    if(!std::filesystem::exists(OUTPUT_SUPPLEMENTAL_PATH))
        std::filesystem::create_directory(OUTPUT_SUPPLEMENTAL_PATH);

    std::vector<std::string> no_logs_meshes, invalid_labeling_meshes, missing_figs_meshes;

    LatexDoc latex(std::string(OUTPUT_SUPPLEMENTAL_PATH) + "supplemental.tex");
    latex.add_subpage(std::filesystem::relative(INPUT_COVER_PAGE,OUTPUT_SUPPLEMENTAL_PATH));//argument = where is the cover page relative to the output folder

    std::set<std::filesystem::directory_entry> entries_sorted_by_name;
    for(const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(INPUT_DATA_PATH))//interators have no order -> sort by alphabetical order with a set
        entries_sorted_by_name.insert(dir_entry);

    int success_count = 0;

    for(const std::filesystem::directory_entry& dir_entry : entries_sorted_by_name) {
        if(!dir_entry.is_directory()) continue;
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
                             + std::to_string(entries_sorted_by_name.size()), "green");
    std::cout << "---------------------------" << std::endl;
    return 0;
}

