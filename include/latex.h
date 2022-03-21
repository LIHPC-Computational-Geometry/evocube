#include <iostream>
#include <fstream>
#include <string>

class MainLatexDoc {

public:
    MainLatexDoc(std::string filename);
    ~MainLatexDoc();

    void include_pdf(std::string pdf_filename);

private:
    std::ofstream ofs;

};

class SubLatexDoc {

public:
    SubLatexDoc(std::string filename);
    ~SubLatexDoc();

private:
    std::ofstream ofs;
};