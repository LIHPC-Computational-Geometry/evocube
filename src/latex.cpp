#include <iostream>
#include <fstream>
#include <string>

#include "latex.h"

MainLatexDoc::MainLatexDoc(std::string filename)
    : ofs(filename)
{
    ofs << "\\documentclass{article}" << std::endl;
    ofs << "\\usepackage{pdfpages}" << std::endl;
    ofs << "\\begin{document}" << std::endl;
}

MainLatexDoc::~MainLatexDoc() {
    ofs << "\\end{document}" << std::endl;
    ofs.close();
}

void MainLatexDoc::include_pdf(std::string pdf_filename) {
    ofs << "\\includepdf[pages=-]{"<< pdf_filename <<"}" << std::endl;
}

SubLatexDoc::SubLatexDoc(std::string filename)
    : ofs(filename)
{
    ofs << "\\let\\mypdfximage\\pdfximage" << std::endl;
    ofs << "\\def\\pdfximage{ \\immediate\\mypdfximage }" << std::endl;
    ofs << std::endl;
    ofs << "\\documentclass{article}" << std::endl;
    ofs << std::endl;
    ofs << "\\usepackage[top = 1cm, bottom = 2cm, left = 1cm, right = 1cm]{geometry}" << std::endl;
    ofs << "\\usepackage{graphicx}" << std::endl;
    ofs << "\\usepackage{array}" << std::endl;
    ofs << "\\usepackage{tabularx}" << std::endl;
    ofs << "\\newcolumntype{Y}{>{\\centering\\arraybackslash}X}" << std::endl;
    ofs << "\\usepackage{xcolor}" << std::endl;
    ofs << std::endl;
    ofs << "\\begin{document}" << std::endl;
    ofs << "\\small" << std::endl;
}

SubLatexDoc::~SubLatexDoc() {
    ofs << "\\end{document}" << std::endl;
    ofs.close();
}