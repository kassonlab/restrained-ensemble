//
// Created by jmh on 6/23/17.
//

#include "common.h"
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

using boost::str;
using boost::format;

template <typename T>
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

gromacs_files generate_gromacs_filenames(int ensemble_number, prefixes prefs, int replica,
                                         bool start_tpr, bool protein_ndx) {
    gromacs_files names;
    auto ensemble_path = prefs.ensemble_path;
    auto directory_prefix = prefs.directory_prefix;
    auto prod_prefix = prefs.prod_prefix;
    auto start_prefix = prefs.start_prefix;

    if (ensemble_number > 1) {
        names.xtc = str(format("%s/%s_%02i/%s.part%04i.xtc")
                        % ensemble_path % directory_prefix % replica % prod_prefix % (ensemble_number - 1));
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

void calculate_histogram(std::vector<pair_data> vec_pd,
                         const char *out_filename,
                         parameters params,
                         int ensemble_number) {

    if (boost::filesystem::exists(out_filename)) {
        backup_file(out_filename, ensemble_number);
    }

    std::ofstream hist_file;
    hist_file.open(out_filename);

    hist_file << params.bin_width << "," << params.sigma << ",";
    hist_file << params.min_dist << "," << params.max_dist << "\n";

    auto norm = 1.0 / sqrt(2.0 * M_PI * pow(params.sigma, 2.0));

    for (auto &pd: vec_pd) {
        auto exp_data = &pd.exp_distribution;
        if (exp_data->empty()) {
            char error[BUFFER_LENGTH];
            snprintf(error, BUFFER_LENGTH, "Cannot calculate histogram if experimental distribution is empty");
            throw std::invalid_argument(error);
        }

        auto sim_data = &pd.sim_dist_data;
        auto num_bins = exp_data->size();
        auto num_samples = sim_data->size();
        auto sample_norm = (1.0 / num_samples) * norm;

        for (int n = 0; n < num_bins; ++n) {
            double h_ij_n{0};
            for (auto &sample_dist: *sim_data) {
                h_ij_n += sample_norm * exp(-pow(n * params.bin_width - sample_dist,
                                                 2.0) / (2 * pow(params.sigma, 2)));
            }
            if (n == 0) hist_file << h_ij_n;
            else hist_file << "," << h_ij_n;
        }
        hist_file << "\n";
    }
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
        argv.push_back(str(format("%s/%s_%02i/") % prefs.ensemble_path % prefs.directory_prefix % i));
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