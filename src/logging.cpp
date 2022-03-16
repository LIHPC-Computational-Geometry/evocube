
#include "logging.h"
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <nlohmann/json.hpp>

void coloredPrint(std::string text, std::string color){
    std::cout << "\033[1;";
    if (color == "black")
        std::cout << "30m";
    else if (color == "red")
        std::cout << "31m";
    else if (color == "green")
        std::cout << "32m";
    else if (color == "yellow")
        std::cout << "33m";
    else if (color == "blue")
        std::cout << "34m";
    else if (color == "magenta")
        std::cout << "35m";
    else if (color == "cyan")
        std::cout << "36m";
    else std::cout << "37m";

    std::cout << text << "\033[0m" << std::endl;
}

nlohmann::json readJSON(std::string filepath){
    std::ifstream i(filepath);
    nlohmann::json j;
    if (i.good())
        i >> j;
    i.close();
    return j;
}

void writeJSON(nlohmann::json j, std::string filepath){
    std::ofstream o(filepath);
    o << std::setw(4) << j << std::endl;
    o.close();
}

void resetLogFile(std::string filepath){
    writeJSON(nlohmann::json(), filepath);
}

void fillLogInfo(std::string tag, std::string filepath, std::string value){
    nlohmann::json j = readJSON(filepath);
    j[tag] = value;
    writeJSON(j, filepath);
}

void fillLogInfo(std::string tag1, std::string tag2, std::string filepath, std::string value){
    nlohmann::json j = readJSON(filepath);
    j[tag1][tag2] = value;
    writeJSON(j, filepath);
}

void removeLogInfo(std::string tag, std::string filepath){
    nlohmann::json j = readJSON(filepath);
    j.erase(tag);
    writeJSON(j, filepath);
}

std::string readLogInfo(std::string tag, std::string filepath){
    nlohmann::json j = readJSON(filepath);
    return j[tag];
}

std::map<std::string, std::string> readAllLogInfo(std::string filepath){
    std::map<std::string, std::string> result;
    nlohmann::json j = readJSON(filepath);
    for(auto it: j.items()){
        result[it.key()] = it.value();
    }
    return result;
}

std::vector<std::string> getLogTags(std::string filepath){
    std::vector<std::string> result;
    nlohmann::json j = readJSON(filepath);
    for(auto it: j.items()){
        result.push_back(it.key());
    }
    return result;
}

void printLog(std::string filepath){
    nlohmann::json j = readJSON(filepath);
    for(auto it: j.items()){
        std::cout << it.key() << "\t| " << it.value() << '\n';
    }
}


double measureTime(std::chrono::time_point<std::chrono::steady_clock> t1,
                   std::chrono::time_point<std::chrono::steady_clock> t2){
    return static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds> 
                              (t2 - t1).count()) / 1000.0;
}