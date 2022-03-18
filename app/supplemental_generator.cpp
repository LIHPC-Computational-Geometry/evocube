#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <array>
#include <map>

namespace fs = std::filesystem;


struct run_result {
    int nb_tets = 0;
    int nb_tet_flipped = 0;
    double time_input_polycube = 0;
    double time_algo = 0;
    double time_solo_milp = 0;
    int nb_solver_runs_coarse = 0;
    int nb_solver_runs_fine = 0;
    bool unsolvable_milp = false;
    bool obviously_bad_flagging = false;
    bool GARBAGE_IN_GARBAGE_OUT = false;
};

constexpr run_result TIMEDOUT = { -1, -1, -1, -1, -1, -1,-1,0,0,0 };

inline bool res_is_timeout(run_result res) {
    return res.nb_tets == -1;
}

static run_result parse_results(const std::string& filename) {

    run_result RES;

    std::ifstream res_file(filename);
    std::string line;
    std::getline(res_file, line);
    std::getline(res_file, line);
    std::getline(res_file, line);

    while (!res_file.eof()) {
        std::getline(res_file, line); std::stringstream sst(line);
        if (sst.str().empty()) break;
        std::string part_name; double data;
        sst >> data; sst >> part_name;
        if (part_name == "setConstraints" || part_name == "bendonlyFrameField" || part_name == "compute_charts" || part_name == "compute_Po")
            RES.time_input_polycube += data;
        else if (part_name == "compute_Pf") {
            RES.time_algo += data;
            RES.time_solo_milp += data;
        }
        else if (part_name == "inverse_Pf" || part_name == "compute_Pr")
            RES.time_algo += data;
    }
    if (res_file.eof()) exit(1);
    std::getline(res_file, line);
    std::getline(res_file, line);
    std::getline(res_file, line);
    while (!res_file.eof()) {
        std::getline(res_file, line); std::stringstream sst(line);
        if (sst.str().empty()) break;
        std::string part_name; double data;
        sst >> data; sst >> part_name;
        if (part_name == "nbTets") RES.nb_tets = (int)data;
        if (part_name == "nbFlipped") RES.nb_tet_flipped = (int)data;
        if (part_name == "nb_MILP_SOLVER_run(scale=0.000000,ovelaps=1)") RES.nb_solver_runs_coarse = (int)data;
        if (part_name == "nb_MILP_SOLVER_run(scale=1.000000,ovelaps=1)") RES.nb_solver_runs_fine = (int)data;
        if (part_name == "Unsolvable") RES.unsolvable_milp = (bool)data;
        if (part_name == "BadCharts") RES.obviously_bad_flagging = (bool)data;
        if (part_name == "chartToBad") RES.GARBAGE_IN_GARBAGE_OUT = (bool)data;
    }
    res_file.close();
    return RES;
}

std::string make_table(const run_result& res, const std::string& filename) {
    std::stringstream sst;

    std::string tmp = filename;
    tmp.replace(tmp.find("_"), 1, "\\_");
    
    sst << "\\begin{tabularx}{\\linewidth} { | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | } " << std::endl;
    sst << "\\hline" << std::endl;
    sst << "Model Name& \\#tets& \\#flipped tets& Input gen.time& Method time& Milp time& \\#solver call (coarse)& \\#solver call (fine)& Obv. bad flagging& Unsolv. Milp& Lost case \\\\ " << std::endl;
    sst << "\\hline" << std::endl;
    sst << "{ \\tiny " << tmp << "} &";
    if (res_is_timeout(res)) {
        sst << "{\\color{red} -1}" << "&";
        sst << "{\\color{red} -1}" << "&";
        sst << "{\\color{red}TIMEDOUT}" << "&";
        sst << "{\\color{red}TIMEDOUT}" << "&";
        sst << "{\\color{red}TIMEDOUT}" << "&";
        sst << "{\\color{red} -1}" << "&";
        sst << "{\\color{red} -1}" << "&";
        sst << "{\\color{red} 0}" << "&";
        sst << "{\\color{red} 0}" << "&";
        sst << "{\\color{red} 1}" << "\\\\" << std::endl;
    }
    else {
        sst << res.nb_tets << "&";
        sst << res.nb_tet_flipped << "&";
        sst << res.time_input_polycube << "&";
        sst << res.time_algo << "&";
        sst << res.time_solo_milp << "&";
        sst << res.nb_solver_runs_coarse << "&";
        sst << res.nb_solver_runs_fine << "&";
        sst << res.obviously_bad_flagging << "&";
        sst << res.unsolvable_milp << "&";
        sst << res.GARBAGE_IN_GARBAGE_OUT << "\\\\" << std::endl;
    }
    sst << "\\hline" << std::endl;
    sst << "\\end{tabularx}" << std::endl;

    return sst.str();
}

/*
* \usepackage{ array }
\usepackage{ tabularx }
\newcolumntype{ Y }{ > {\centering\arraybackslash}X}
\begin{ document }

\begin{ figure }
\centering
\begin{ tabularx }{\linewidth} { | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | }
\hline
Model Name& \#tets& \#flipped tets& Input gen.time& Method time& Milp time& \#solver call(coarse)& \#solver call(fine)& Obv.bad flagging& Unsol.Milp& Lost case \\
\hline
1.1 & 1.2 & 1.3 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
\hline
\end{ tabularx } \\
*/


static std::string make_array_data(const std::string& dir_name) {
    int nb_timeout = 0;
    for (const auto& entry : fs::directory_iterator(dir_name)) {
        std::string ext = entry.path().extension().string();
        if (ext != ".step") continue;
        std::string filename = entry.path().stem().string();


        std::string path2file = fs::absolute(entry.path().parent_path()).string();
        std::replace(path2file.begin(), path2file.end(), '\\', '/');
        std::string abs_path = path2file + "/" + filename;

        std::cerr << "Adding: " << abs_path << std::endl;

        run_result res = TIMEDOUT;
        if (!fs::exists(abs_path + ".timeout") && fs::exists(abs_path + "_timing.txt"))
            res = parse_results(abs_path + "_timing.txt");
        else {
            nb_timeout++;
            continue;
        }


    }
    std::stringstream sst;




    return sst.str();

}


int main(int argc, char** argv) {
    if (argc != 3){
        std::cout << "Usage: " << argv[0] << " dir doc " << std::endl;
        return 1;
    } 



    std::string dir_name = argv[1]; 
    std::string latex_file_name = argv[2];

    std::vector<fs::path> step_meshes_path;
    for (const auto& entry : fs::directory_iterator(dir_name)) {
        if (entry.path().extension().string() != ".step") continue;
        step_meshes_path.push_back(entry.path());
    }
    constexpr int nb_of_mesh_per_pdf = 100;
    int nb_of_latex = (int)step_meshes_path.size() / nb_of_mesh_per_pdf;
    for (int i = 0; i < nb_of_latex; i++) {
        
        std::string small_pdf = latex_file_name + "_" + std::to_string(i) + ".tex";
        
        std::cerr << "-> " << small_pdf << std::endl;
        
        std::ofstream latex_file(small_pdf);
        latex_file << "\\let\\mypdfximage\\pdfximage" << std::endl;
        latex_file << "\\def\\pdfximage{ \\immediate\\mypdfximage }" << std::endl;
        latex_file << std::endl;
        latex_file << "\\documentclass{article}" << std::endl;
        latex_file << std::endl;
        latex_file << "\\usepackage[top = 1cm, bottom = 2cm, left = 1cm, right = 1cm]{geometry}" << std::endl;
        latex_file << "\\usepackage{graphicx}" << std::endl;
        latex_file << "\\usepackage{array}" << std::endl;
        latex_file << "\\usepackage{tabularx}" << std::endl;
        latex_file << "\\newcolumntype{Y}{>{\\centering\\arraybackslash}X}" << std::endl;
        latex_file << "\\usepackage{xcolor}" << std::endl;
        latex_file << std::endl;
        latex_file << "\\begin{document}" << std::endl;
        latex_file << "\\small" << std::endl;
        latex_file << "\\setcounter{page}{" << i*nb_of_mesh_per_pdf << "}" << std::endl;

        for (int j = i * nb_of_mesh_per_pdf; j < std::min((int)step_meshes_path.size(), (i + 1) * nb_of_mesh_per_pdf); j++) {
            fs::path entry = step_meshes_path[j];
            std::string filename = entry.stem().string();
            std::string path2file = fs::absolute(entry.parent_path()).string();
            std::replace(path2file.begin(), path2file.end(), '\\', '/');
            std::string abs_path = path2file + "/" + filename;

            std::cerr << "Adding: " << abs_path << std::endl;

            run_result res = TIMEDOUT;
            if (!fs::exists(abs_path + ".timeout") && fs::exists(abs_path + "_timing.txt"))
                res = parse_results(abs_path + "_timing.txt");

            latex_file << "\\begin{figure}" << std::endl;
            latex_file << "\\centering" << std::endl;
            latex_file << make_table(res, filename);
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_0.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_1.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_2.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_3.png}\\\\" << std::endl;

            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_pf_0.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_pf_1.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_pf_2.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_pf_3.png}\\\\" << std::endl;

            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_polycube_0.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_polycube_1.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_polycube_2.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_polycube_3.png}\\\\" << std::endl;

            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_hexmesh_0.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_hexmesh_1.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_hexmesh_2.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_coarse_hexmesh_3.png}\\\\" << std::endl;

            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_hexmesh_0.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_hexmesh_1.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_hexmesh_2.png}" << std::endl;
            latex_file << "\\includegraphics[width=0.24\\linewidth]{" + abs_path + "_hexmesh_3.png}\\\\" << std::endl;

            latex_file << "\\end{figure}" << std::endl;
            latex_file << "\\clearpage" << std::endl;
            latex_file << "\\pagebreak" << std::endl;

        }
        latex_file << "\\end{document}" << std::endl;
        latex_file.close();
    }


    std::ofstream latex_file(latex_file_name + ".tex");
    

    latex_file << "\\documentclass{article}" << std::endl;
    latex_file << "\\usepackage{pdfpages}" << std::endl;
    latex_file << "\\begin{document}" << std::endl;
    for (int i = 0; i < nb_of_latex; i++) {
        std::string small_pdf = latex_file_name + "_" + std::to_string(i) + ".pdf";

        latex_file << "\\includepdf[pages=-]{"<< small_pdf <<"}" << std::endl;
    }
    latex_file << "\\end{document}" << std::endl;

    return 0;
}

