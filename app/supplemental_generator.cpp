#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "latex.h"
#include "logging.h"

#define INPUT_DATA_PATH             "../data/"
#define INPUT_COVER_PAGE            "../supplemental/cover_page.tex"
#define OUTPUT_SUPPLEMENTAL_PATH    "../supplemental/"

int main(int argc, char** argv) {

    if(!std::filesystem::exists(OUTPUT_SUPPLEMENTAL_PATH))
        std::filesystem::create_directory(OUTPUT_SUPPLEMENTAL_PATH);

    unsigned int nb_incomplete_meshes = 0;

    LatexDoc latex(std::string(OUTPUT_SUPPLEMENTAL_PATH) + "supplemental.tex");
    latex.add_subpage(std::filesystem::relative(INPUT_COVER_PAGE,OUTPUT_SUPPLEMENTAL_PATH));//argument = where is the cover page relative to the output folder

    std::set<std::filesystem::directory_entry> entries_sorted_by_name;
    for(const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(INPUT_DATA_PATH))//interators have no order -> sort by alphabetical order with a set
        entries_sorted_by_name.insert(dir_entry);

    for(const std::filesystem::directory_entry& dir_entry : entries_sorted_by_name) {
        if(!dir_entry.is_directory()) continue;
        nb_incomplete_meshes += latex.add_mesh(dir_entry.path());
    }

    if(nb_incomplete_meshes > 0)
        coloredPrint( std::to_string(nb_incomplete_meshes) + " mesh(es) have some figures missing", "red");
    else
        coloredPrint( "No figure missing", "green");

    return 0;
}

