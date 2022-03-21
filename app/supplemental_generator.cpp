#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <array>
#include <map>

#include "latex.h"
#include "logging.h"

#define OUTPUT_SUPPLEMENTAL_PATH    "../supplemental/"

int main(int argc, char** argv) {

    if(!std::filesystem::exists(OUTPUT_SUPPLEMENTAL_PATH))
        std::filesystem::create_directory(OUTPUT_SUPPLEMENTAL_PATH);

    unsigned int nb_incomplete_meshes = 0;

    LatexDoc latex(std::string(OUTPUT_SUPPLEMENTAL_PATH) + "supplemental.tex");
    nb_incomplete_meshes += latex.add_mesh("../data/B21/");
    nb_incomplete_meshes += latex.add_mesh("../data/M8/");
    nb_incomplete_meshes += latex.add_mesh("../data/M9/");

    if(nb_incomplete_meshes > 0)
        coloredPrint( std::to_string(nb_incomplete_meshes) + " mesh(es) have some figures missing", "red");
    else
        coloredPrint( "No figure missing", "green");

    return 0;
}

