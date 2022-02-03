#pragma once
#include <Eigen/Core>
#include <iostream>
#include <fstream>

class BndToTetConverter {
public:
    BndToTetConverter(std::string input_file){
        int n_values;
        int line;
        std::ifstream file(input_file);
        
        if (file.is_open())
        {
            file >> n_values;

            table_.resize(n_values);
            for (int i=0; i < n_values; i++){
                file >> table_(i);
            }
            file.close();
        }

        else std::cout << "BndToTetConverter: Unable to open file";
    }

    BndToTetConverter(const Eigen::VectorXi table) : table_(table) {}

    void writeTo(std::string write_path){
        std::ofstream file(write_path);
        if (file.is_open()){
            file << table_.rows();
            file << "\n";
            for (int i=0; i<table_.rows(); i++){
                file << table_(i);
                file << "\n";
            }
            file.close();
        }
        else std::cout << "BndToTetConverter: Unable to write to file";
    }
//private:
    Eigen::VectorXi table_;
};
