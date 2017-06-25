//
// Created by jmh on 6/23/17.
//

#ifndef ROUX_COMMON_H
#define ROUX_COMMON_H

#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

constexpr int BUFFER_LENGTH = 500;
typedef std::vector<std::vector<float>> simdata;
typedef std::vector<std::string> vecofstrings;

struct pair_data {
    int xvg_col_num;
    double k;
    std::pair<int, int> residue_ids;
    std::vector<float> exp_distribution;
    std::vector<float> sim_dist_data;
    std::vector<float> sim_time_data;
};

template<class Archive>
void serialize(Archive &ar, pair_data &pd, const unsigned int version) {
    ar & pd.xvg_col_num;
    ar & pd.k;
    ar & pd.residue_ids;
    ar & pd.exp_distribution;
    ar & pd.sim_dist_data;
    ar & pd.sim_time_data;
}

struct gromacs_files {
    std::string xtc;
    std::string tpr;
    std::string ndx;
    std::string xvg;
    std::string gro;
};

struct prefixes {
    std::string ensemble_path;
    std::string directory_prefix;
    std::string prod_prefix;
    std::string start_prefix;
};

struct parameters {
    int                 ref;
    int                 chains;
    bool                aa;
    double              bin_width;
    double              sigma;
    double              min_dist;
    double              max_dist;
    unsigned long       num_bins;
    unsigned long       num_steps;
    unsigned long       boxcar_parts;
    unsigned long       num_pairs;
    unsigned long       num_replicas;
    std::vector<int>    replicas;
    std::string         lipid;
};

struct input_filenames{
    std::string exp_filename;
    std::string mdp_template;
    std::string roux_mdp;
    std::string gmx_exe;
    std::string log_dir;
    std::string differences;
};

typedef std::vector<pair_data> setofpairs;
#endif //ROUX_COMMON_H
