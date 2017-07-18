//
// Created by jmh on 6/23/17.
//

#include "common.h"
#include <fstream>
#include <dirent.h>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include "rouxfileio.h"

using boost::str;
using boost::format;

std::string strip_filename(std::string filename) {
    std::string stripped_filename;
    boost::replace_all(filename, "_cntr", "");
    boost::replace_all(filename, "_whole", "");
    stripped_filename = filename.substr(0, filename.size() - 4);
    return stripped_filename;
}

template<typename T>
std::vector<T> unique_elements(std::vector<T> vec) {
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    return vec;
}

std::tuple<std::vector<int>, std::vector<int>> unique_from_pairs(std::vector<std::pair<int, int>> pairs) {
    std::vector<int> temp_first_resids, temp_second_resids, first_resids, second_resids;
    for (int i = 0; i < pairs.size(); ++i) {
        temp_first_resids.push_back(pairs.at(i).first);
        temp_second_resids.push_back(pairs.at(i).second);
    }
    first_resids = unique_elements(temp_first_resids);
    second_resids = unique_elements(temp_second_resids);
    return std::make_tuple(first_resids, second_resids);
}

void backup_file(std::string filename, int ensemble_number) {
    auto length = filename.size();
    char buffer[BUFFER_LENGTH];
    snprintf(buffer, BUFFER_LENGTH, "mv %s %s.part%04i.%s",
             filename.c_str(),
             filename.substr(0, length - 4).c_str(),
             ensemble_number,
             filename.substr(length - 3, length).c_str());
    system(buffer);
}

int find_last_run_number(std::string ensemble_path,
                         std::string directory_prefix,
                         int replica) {

    int max_ensemble_number = 0;

    char directory_name[BUFFER_LENGTH];
    std::snprintf(directory_name, BUFFER_LENGTH, "%s/%s_%02i/",
                  ensemble_path.c_str(), directory_prefix.c_str(), replica);

    std::string pattern = "^.+part([0-9][0-9][0-9][0-9])[.]log$";
    boost::regex reg{pattern};
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(directory_name)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            boost::cmatch cm;
            if (boost::regex_search(ent->d_name, cm, reg)) {
                std::string lognum = cm.str(1);
                int intlognum = std::stoi(lognum);
                max_ensemble_number = std::max(max_ensemble_number, intlognum);
            }
        }
        closedir(dir);
    }
    return max_ensemble_number + 1;
}

gromacs_files generate_gromacs_filenames(int ensemble_number, prefixes prefs, int replica,
                                         bool start_tpr, bool protein_ndx) {
    gromacs_files names;
    auto ensemble_path = prefs.ensemble_path.c_str();
    auto directory_prefix = prefs.directory_prefix.c_str();
    auto prod_prefix = prefs.prod_prefix.c_str();
    auto start_prefix = prefs.start_prefix.c_str();

    if (ensemble_number > 1) {
        names.xtc = str(format("%s/%s_%02i/%s.part%04i.xtc") % ensemble_path % directory_prefix % replica % prod_prefix
                        % (ensemble_number - 1));
        names.tpr = str(format("%s/%s_%02i/%s.tpr") % ensemble_path % directory_prefix % replica % prod_prefix);
        names.xvg = str(format("%s/%s_%02i/%s_pullx.part%04i.xvg")
                        % ensemble_path % directory_prefix % replica % prod_prefix % (ensemble_number - 1));
        names.gro = str(format("%s/%s_%02i/%s.part%04i.gro")
                        % ensemble_path % directory_prefix % replica % prod_prefix % (ensemble_number - 1));

    } else if (ensemble_number == 1) {
        names.xtc = str(format("%s/%s_%02i/%s.gro") % ensemble_path % directory_prefix % replica % start_prefix);
        if (start_tpr)
            names.tpr = str(format("%s/%s_%02i/%s.tpr") % ensemble_path % directory_prefix % replica % start_prefix);
        else
            names.tpr = str(format("%s/%s_%02i/%s.tpr") % ensemble_path % directory_prefix % replica % prod_prefix);

        names.xvg = str(format("%s/%s_%02i/%s_dist.xvg") % ensemble_path % directory_prefix % replica % start_prefix);
        names.gro = str(format("%s/%s_%02i/%s.gro") % ensemble_path % directory_prefix % replica % start_prefix);
    }
    if (protein_ndx)
        names.ndx = str(format("%s/%s_%02i/protein.ndx") % ensemble_path % directory_prefix % replica);
    else
        names.ndx = str(format("%s/%s_%02i/roux_groups.ndx") % ensemble_path % directory_prefix % replica);

    return names;
}

std::vector<gromacs_files> generate_gromacs_filenames(int ensemble_number,
                                                      prefixes prefs,
                                                      std::vector<int> replicas,
                                                      bool start_tpr,
                                                      bool protein_ndx) {

    std::vector<gromacs_files> gro_files;
    for (auto &replica: replicas) {
        gro_files.push_back(generate_gromacs_filenames(ensemble_number,
                                                       prefs,
                                                       replica,
                                                       start_tpr,
                                                       protein_ndx));
    }
    return gro_files;
}


void make_ndx(std::string tpr_filename,
              std::string ndx_filename,
              std::string selection,
              std::string gmx_exe,
              bool aa,
              bool rewrite) {

    if (rewrite || !boost::filesystem::exists(ndx_filename)) {
        if (selection.empty()) {
            if (aa) {
                selection = "mol 1 and name CA; mol 2 and name CA; "
                        "(mol 1 and name CA) or (mol 2 and name CA); mol 1 or mol 2";
            } else {
                selection = "mol 1 and name BB; mol 2 and name BB; "
                        "(mol 1 and name BB) or (mol 2 and name BB); mol 1 or mol 2";
            }
        }

        char command[BUFFER_LENGTH];
        snprintf(command, BUFFER_LENGTH,
                 "%s select -s %s -select \"%s\" -on %s",
                 gmx_exe.c_str(),
                 tpr_filename.c_str(),
                 selection.c_str(),
                 ndx_filename.c_str());
        system(command);

    } else warn_file_exists(ndx_filename.c_str());
}

void make_xvg(std::string tpr_filename,
              std::string xtc_filename,
              std::string xvg_filename,
              int chains,
              bool aa,
              bool rewrite,
              std::vector<std::pair<int, int>> pairs,
              std::string gmx_exe) {

    std::string selection, str_1AA, str_2AA, str_1BB, str_2BB;
    std::vector<int> first, second;

    if (rewrite or !boost::filesystem::exists(xvg_filename)) {
        for (auto &pair: pairs) {
            first.push_back(pair.first);
            second.push_back(pair.second);
        }

        if (chains == 2) {
            selection = "";
            str_1AA = "mol 1 and name CA and resid ";
            str_2AA = "mol 2 and name CA and resid ";
            str_1BB = "mol 1 and name BB and resid ";
            str_2BB = "mol 2 and name BB and resid ";
        } else {
            selection = "";
            str_1AA = "name CA and resid ";
            str_2AA = str_1AA;
            str_1BB = "name BB and resid ";
            str_2BB = str_1BB;

        }

        auto num_pairs = pairs.size();
        for (uint i = 0; i < num_pairs; i++) {
            std::string pair_selection{""};
            if (aa)
                pair_selection = str(format("(%s%i) or (%s%i)")
                                     % str_1AA.c_str() % first.at(i)
                                     % str_2AA.c_str() % second.at(i));
            else
                pair_selection = str(format("(%s%i) or (%s%i)\"")
                                     % str_1BB.c_str() % first.at(i)
                                     % str_2BB.c_str() % second.at(i));

            if (num_pairs == 1) {
                selection += "\"" + pair_selection + "\"";
            } else {
                if (i == 0)
                    selection += "\"" + pair_selection + ";";
                else if (i == num_pairs - 1)
                    selection += " " + pair_selection + "\"";
                else {
                    selection += " " + pair_selection + ";";
                }
            }

        }

        std::string command = str(format("%s distance -f %s -s %s -oall %s -select %s")
                                  % gmx_exe % xtc_filename % tpr_filename % xvg_filename % selection);
        system(command.c_str());
    } else warn_file_exists(xvg_filename.c_str());
}

void pre_process(gromacs_files name,
                 int chains,
                 bool aa,
                 std::vector<std::pair<int, int>> pairs,
                 std::string gmx_exe,
                 bool rewrite) {

    std::string echo1, echo2, whole_name, cntr_name, selection;

    if (chains == 2) { //TODO: We got some problems here: where is selection defined?
        make_ndx(name.tpr, name.ndx, selection, gmx_exe, aa, rewrite);
        echo1 = "echo 3 3 | ";
        echo2 = "echo 1 3 | ";
    } else {
        make_ndx(name.tpr, name.ndx, "group Protein", gmx_exe, aa, rewrite);
        echo1 = "";
        echo2 = "";
    }

    whole_name = strip_filename(name.xtc) + "_whole.xtc";
    cntr_name = strip_filename(name.xtc) + "_cntr.xtc";

    char buffer[BUFFER_LENGTH];

    if (boost::filesystem::exists(whole_name) && !rewrite) {
        warn_file_exists(whole_name.c_str());
    } else {
        snprintf(buffer, BUFFER_LENGTH, "%s %s trjconv -f %s -s %s -pbc whole -n %s -o %s",
                 echo1.c_str(),
                 gmx_exe.c_str(),
                 name.xtc.c_str(),
                 name.tpr.c_str(),
                 name.ndx.c_str(),
                 whole_name.c_str());
        system(buffer);
    }

    if (boost::filesystem::exists(cntr_name) && !rewrite) {
        warn_file_exists(cntr_name.c_str());
    } else {
        snprintf(buffer, BUFFER_LENGTH,
                 "%s %s trjconv -f %s -s %s -pbc mol -ur compact -center -n %s -o %s",
                 echo2.c_str(),
                 gmx_exe.c_str(),
                 whole_name.c_str(),
                 name.tpr.c_str(),
                 name.ndx.c_str(),
                 cntr_name.c_str());
        system(buffer);
    }
    make_xvg(name.tpr, cntr_name, name.xvg, chains, aa, rewrite, pairs, gmx_exe);

}

void calculate_histogram(std::vector<pair_data> vec_pd,
                         const char *out_filename,
                         parameters params,
                         int ensemble_number) {


    if (boost::filesystem::exists(out_filename)) {
        backup_file(out_filename, ensemble_number - 1);
    }

    std::ofstream hist_file;
    hist_file.open(out_filename);

    hist_file << params.bin_width << "," << params.sigma << ",";
    hist_file << params.min_dist << "," << params.max_dist << "\n";

    auto norm = 1.0 / sqrt(2.0 * M_PI * pow(params.sigma, 2.0));

    for (auto &pd: vec_pd) {
        auto &exp_data = pd.exp_distribution;
        if (exp_data.empty()) {
            char error[BUFFER_LENGTH];
            snprintf(error, BUFFER_LENGTH, "Cannot calculate histogram if experimental distribution is empty");
            throw std::invalid_argument(error);
        }

        auto &sim_data = pd.sim_dist_data;
        auto num_bins = exp_data.size();
        auto num_samples = sim_data.size();
        auto sample_norm = (1.0 / num_samples) * norm;

        for (int n = 0; n < num_bins; ++n) {
            double h_ij_n{0};
            for (auto &sample_dist: sim_data) {
                h_ij_n += sample_norm * exp(-pow(n * params.bin_width - sample_dist,
                                                 2.0) / (2 * pow(params.sigma, 2)));
                if (h_ij_n < pow(10, -10))
                    h_ij_n = 0;
            }
            if (n == 0) hist_file << h_ij_n;
            else hist_file << "," << h_ij_n;
        }
        hist_file << "\n";
    }
    hist_file.close();
}

std::vector<std::string> mdrun_chararray(parameters params, prefixes prefs, int &argc) {
    std::array<std::string, 10> flags = {"-x", "-cpi", "-cpo", "-e", "-c", "-g", "-s", "-o", "-px", "-pf"};
    std::vector<std::string> argv;

    argv.push_back("mdrun");
    for (auto flag: flags) {
        argv.push_back(flag);
        if (flag.find("px") != std::string::npos) {
            argv.push_back(prefs.prod_prefix + "_pullx");
        } else if (flag.find("pf") != std::string::npos) {
            argv.push_back(prefs.prod_prefix + "_pullf");
        } else {
            argv.push_back(prefs.prod_prefix);
        }
    }
    argv.push_back("-v");
    //argv.push_back("1");
    argv.push_back("-multidir");

// MULTIDIR STRING

    for (auto i: params.replicas) {
        argv.push_back(str(format("%s/%s_%02i/") % prefs.ensemble_path.c_str() % prefs.directory_prefix.c_str() % i));
    }
    argv.push_back("-noappend");
    argv.push_back("-pin");
    argv.push_back("on");
    argv.push_back("-nsteps");
    argv.push_back(std::to_string(params.num_steps));
//    argv.push_back("-nb");
//    argv.push_back("cpu");

    if (!params.aa) {
        argv.push_back("-rdd");
        argv.push_back("2.0");
    }

    argc = (int) argv.size();

    return argv;
}