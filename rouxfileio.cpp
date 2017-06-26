//
// Created by jmh on 6/23/17.
//

#include <fstream>
#include "json/json.h"
#include "rouxfileio.h"
#include "rouxalgorithms.h"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>

namespace mpi = boost::mpi;


template<typename T>
std::stringstream &operator>>(std::stringstream &stream, std::pair<T, T> &p) {
    stream >> p.first >> p.second;
    return stream;
}

template<typename T>
std::vector<T> str2vec(std::string tokens) {
    std::vector<T> result;
    std::stringstream stream(tokens);
    T n;
    while (stream >> n) {
        result.push_back(n);
    }
    return result;
}

void warn_file_exists(const char *filename) {
    printf("WARNING: The file %s already exists: use rewrite=true to override\n", filename);
}

void parse_ini(const char *ini_filename,
               std::vector<pair_data> &vec_pd,
               input_filenames &in_filenames,
               prefixes &prefs, parameters &pars) {
    using namespace boost::property_tree;
    ptree pt;
    ini_parser::read_ini(ini_filename, pt);
    in_filenames.exp_filename = pt.get<std::string>("files.experimental-filename");
    in_filenames.mdp_template = pt.get<std::string>("files.mdp-template");
    in_filenames.roux_mdp = pt.get<std::string>("files.roux-mdp");
    in_filenames.gmx_exe = pt.get<std::string>("files.gmx-executable");
    in_filenames.log_dir = pt.get<std::string>("files.log-directory");
    prefs.ensemble_path = pt.get<std::string>("prefixes.ensemble-path");
    prefs.directory_prefix = pt.get<std::string>("prefixes.directory-prefix");
    prefs.start_prefix = pt.get<std::string>("prefixes.start-prefix");
    prefs.prod_prefix = pt.get<std::string>("prefixes.production-prefix");
    pars.num_bins = (unsigned long) pt.get<int>("parameters.num-bins");
    pars.bin_width = pt.get<double>("parameters.bin-width");
    pars.sigma = pt.get<double>("parameters.sigma");
    pars.min_dist = pt.get<double>("parameters.min-distance");
    pars.max_dist = pt.get<double>("parameters.max-distance");
    pars.num_steps = (unsigned long) pt.get<int>("parameters.num-steps");
    pars.boxcar_parts = (unsigned long) pt.get<int>("parameters.num-parts-boxcar");
    std::string temp_replicas = pt.get<std::string>("parameters.replica-numbers");
    pars.replicas = str2vec<int>(temp_replicas);
    pars.lipid = pt.get<std::string>("parameters.lipid");
    pars.ref = pt.get<int>("parameters.reference-residue");

    // Get the number of pairs to resize pair data
    auto temp_pairs = pt.get<std::string>("parameters.pairs");
    auto vec_pairs = str2vec<std::pair<int, int>>(temp_pairs);
    std::string temp_k = pt.get<std::string>("parameters.k");
    auto k = str2vec<double>(temp_k);
    pars.num_pairs = vec_pairs.size();
    vec_pd.resize(pars.num_pairs);

    for (int i = 0; i < pars.num_pairs; ++i) {
        vec_pd[i].k = k[i];
        vec_pd[i].residue_ids = vec_pairs[i];
    }

    pars.aa = (pt.get<std::string>("parameters.aa") == "yes");
    pars.chains = pt.get<int>("parameters.chains");
    pars.num_replicas = pars.replicas.size();
}

void read_exp_json(std::string exp_filename, std::vector<pair_data> &pd) {

    std::ifstream exp_stream(exp_filename);
    Json::Value root;
    Json::Reader reader;
    bool parsingSuccessful = reader.parse(exp_stream, root);
    if (!parsingSuccessful) {
        printf("Failed to parse config file\n%s",
               reader.getFormattedErrorMessages().c_str());
    }
    for (auto &pair: pd) {
        char pair_string[100];
        snprintf(pair_string, 100, "%03i %03i", pair.residue_ids.first, pair.residue_ids.second);
        const Json::Value &distribution = root[pair_string];
        assert(distribution.size() != 0);
        for (int i = 0; i < distribution.size(); ++i) {
            pair.exp_distribution.push_back(distribution[i].asFloat());
        }
        assert(pair.exp_distribution.size() != 0);
    }
}

simdata read_sim_xvgs(vecofstrings filenames, unsigned long num_pairs, bool skip_time) {
    simdata result;
    if (skip_time) result.resize(num_pairs); // each row of result corresponds to data for a particular pair
    else result.resize(num_pairs + 1); // each row of result corresponds to data for a particular pair OR the time data

    for (auto &filename : filenames) {
        std::ifstream infile(filename);
        std::string line;
        ulong col{0};

        while (std::getline(infile, line)) {
            col = 0;
            if (line.at(0) != '@' && line.at(0) != '#') {
                auto tokens = strtok(const_cast<char *>(line.c_str()), " \t");
                while (tokens != NULL) {
                    if (skip_time && col != 0) {
                        result[col - 1].push_back(std::stof(tokens));
                    } else if (!skip_time) {
                        result[col].push_back(std::stof(tokens));
                    }
                    ++col;
                    tokens = strtok(NULL, " \t");
                }
                auto num_pairs_xvg = col - 1;
                if (num_pairs_xvg != num_pairs) {
                    char error_message[BUFFER_LENGTH];
                    snprintf(error_message, BUFFER_LENGTH,
                             "The number of pairs (%lu) in the experimental json does not match the number of "
                                     "pairs (%lu) in the simulation xvgs", num_pairs, num_pairs_xvg);
                    throw std::invalid_argument(error_message);
                }
            }

        }
        infile.close();
    }

    return result;
}

std::vector<vecofstrings> scatter_files(vecofstrings filenames, unsigned long num_ranks) {
    std::vector<vecofstrings> result;
    auto num_files = filenames.size();

    /* If we have more files than processes, we need to ~evenly distribute files
     * among the different ranks. If we have fewer files than processes, each
     * process will get one (or zero) files. */

    unsigned long num_used_ranks{0};
    if (num_files < num_ranks) num_used_ranks = num_files;
    else num_used_ranks = num_ranks;
    result.resize(num_used_ranks);

    int rank{0};
    for (int i = 0; i < num_files; ++i) {
        if (rank >= num_used_ranks) rank = 0;
        result[rank].push_back(filenames[i]);
        ++rank;
    }
    return result;
}

void vec2pd(simdata &sim_data, std::vector<pair_data> &vec_pd, bool skip_time) {
    // Do a bit of checking first
    auto sim_len = sim_data.size();
    auto num_pairs = vec_pd.size();
    if (!skip_time) --sim_len;

    if (sim_len != num_pairs) {
        char error_message[BUFFER_LENGTH];
        snprintf(error_message, BUFFER_LENGTH,
                 "The number of pairs (%lu) in the vector of pair data does not match the number of "
                         "pairs (%lu) in the simulation xvgs", num_pairs, sim_len);
        throw std::invalid_argument(error_message);
    }

    // Now start storing data
    for (int i = 0; i < num_pairs; ++i) {
        auto pd = &vec_pd[i];
        std::vector<float> *sim_pd;

        if (skip_time) sim_pd = &sim_data[i];
        else {
            pd->sim_time_data.insert(pd->sim_time_data.end(),
                                     sim_data[0].begin(),
                                     sim_data[0].end());
            sim_pd = &sim_data[i + 1];
        }
        pd->sim_dist_data.insert(pd->sim_dist_data.end(), sim_pd->begin(), sim_pd->end());
    }

}

void mpi_read_xvgs(boost::mpi::communicator &world,
                   setofpairs &vec_pd,
                   vecofstrings filenames,
                   unsigned long num_pairs,
                   bool skip_time) {

    std::vector<simdata> vec_global_sim_data;
    simdata local_sim_data;

    int rank{world.rank()}, num_ranks{world.size()};

    std::vector<vecofstrings> scattered_files;
    scattered_files = scatter_files(filenames, num_ranks);
    if (rank < scattered_files.size()) {
        local_sim_data = read_sim_xvgs(scattered_files.at(rank), num_pairs, skip_time);
    }

    boost::mpi::all_gather(world, local_sim_data, vec_global_sim_data);
    if (rank == 0) {
        for (auto &global_sim_data: vec_global_sim_data) { // Get the data from a particular rank
            if (!global_sim_data.empty())
                vec2pd(global_sim_data, vec_pd, skip_time);
        }
    }
    mpi::broadcast(world, vec_pd, 0);
}

void generate_ndx_files(std::string gmx_exe,
                        gromacs_files name,
                        std::string dat,
                        bool rewrite) {
    if (boost::filesystem::exists(name.ndx) && !rewrite)
        warn_file_exists(name.ndx.c_str());
    else {
        char buffer[BUFFER_LENGTH];
        snprintf(buffer, BUFFER_LENGTH, "%s select -sf %s -f %s -s %s -on %s",
                 gmx_exe.c_str(),
                 dat.c_str(),
                 name.xtc.c_str(),
                 name.tpr.c_str(),
                 name.ndx.c_str());
        system(buffer);
    }

}

vecofstrings make_dat(const char *dat_filename, std::vector<pair_data> vec_pd, parameters params) {
    vecofstrings groups, ionnames;
    std::vector<int> first_resids, second_resids, temp_first_reids, temp_second_resids;
    int first_num, second_num;
    std::string atomname, water;

    std::vector<std::pair<int, int>> pairs;

    for (auto &pd: vec_pd) {
        pairs.push_back(pd.residue_ids);
    }

    std::tie(first_resids, second_resids) = unique_from_pairs(pairs);
    if (params.aa) {
        atomname = {"CB"};
        ionnames = {"NA", "CL"};
        water = "SOL";
    } else {
        atomname = "BB";
        ionnames = {"ION", "ION"};
        water = "W";
    }
    FILE *dat_file = fopen(dat_filename, "w");
    if (params.chains == 1) {
        fprintf(dat_file, "Pull_ref = name %s and resid %i;\n", atomname.c_str(), params.ref);
        groups.push_back("Pull_ref");
        fprintf(dat_file, "System = "
                        "group \"Protein\" or resname %s or resname %s or resname %s or resname %s;\n\n",
                water.c_str(), ionnames.at(0).c_str(), ionnames.at(1).c_str(), params.lipid.c_str());
        first_num = 0;
        second_num = 0;
        for (auto first_resid: first_resids) {
            fprintf(dat_file, "first_%i = name %s and resid %i;\n",
                    first_num, atomname.c_str(), first_resid);
            groups.push_back("first_" + std::to_string(first_num));
            ++first_num;
        }
        for (auto second_resid: second_resids) {
            fprintf(dat_file, "second_%i = name %s and resid %i;\n",
                    second_num, atomname.c_str(), second_resid);
            groups.push_back("second_" + std::to_string(second_num));
            ++second_num;
        }
    }

    if (params.chains == 2) {
        fprintf(dat_file, "Pull_ref = chain B and name %s and resid %i;\n",
                atomname.c_str(), params.ref);
        groups.push_back("Pull_ref");
        fprintf(dat_file, "System = "
                        "group \"Protein\" or resname %s or resname %s or resname %s or resname %s;\n\n",
                water.c_str(), ionnames.at(0).c_str(), ionnames.at(1).c_str(), params.lipid.c_str());
        first_num = 0;
        second_num = 0;
        for (auto first_resid: first_resids) {
            fprintf(dat_file, "first_%i = chain A and name %s and resid %i;\n",
                    first_num, atomname.c_str(), first_resid);
            groups.push_back("first_" + std::to_string(first_num));
            ++first_num;
        }
        for (auto second_resid: second_resids) {
            fprintf(dat_file, "second_%i = chain B and name %s and resid %i;\n",
                    second_num, atomname.c_str(), second_resid);
            groups.push_back("second_" + std::to_string(second_num));
            ++second_num;
        }
    }
    groups.push_back("System");
    for (auto group: groups) {
        fprintf(dat_file, "\n%s;\n", group.c_str());
    }
    fclose(dat_file);
    groups.pop_back();
    return groups;
}

std::vector<std::pair<long, long>> resis_to_groups(std::vector<std::pair<int, int>> pairs) {
    auto n_pairs = pairs.size();
    std::vector<int> first_unique, second_unique;
    std::vector<std::pair<long, long>> groups(n_pairs);
    std::pair<int, int> zero_pair = {0, 0};
    std::fill(groups.begin(), groups.end(), zero_pair);

    std::tie(first_unique, second_unique) = unique_from_pairs(pairs);
    auto n_first = first_unique.size();

    for (int i = 0; i < n_pairs; ++i) {
        auto first_group = std::distance(first_unique.begin(),
                                         std::find(first_unique.begin(),
                                                   first_unique.end(),
                                                   pairs[i].first));
        auto second_group = std::distance(second_unique.begin(),
                                          std::find(second_unique.begin(),
                                                    second_unique.end(),
                                                    pairs[i].second));

        groups[i].first = first_group + 2;
        groups[i].second = second_group + n_first + 2;
    }

    return groups;
}

void write_roux_pull_entry(std::pair<int, int> a_pair, FILE *roux_file, int coord_ind, double k) {
    fprintf(roux_file, "\npull-coord%i-type = roux\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-geometry = distance-reference\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-groups = %i %i 1\n", coord_ind, a_pair.first, a_pair.second);
    fprintf(roux_file, "pull-coord%i-dim = Y Y Y\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-origin = 0.0 0.0 0.0\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-vec = 0.0 0.0 0.0\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-start = no\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-init = 0.0\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-rate = 0.0\n", coord_ind);
    fprintf(roux_file, "pull-coord%i-k = %f\n", coord_ind, k);
    fprintf(roux_file, "pull-coord%i-kB = 0\n", coord_ind);
}

void make_mdp(std::vector<pair_data> vec_pd,
              input_filenames input_files,
              parameters params,
              vecofstrings pull_coord) {

    std::vector<std::pair<int, int>> pairs;
    std::vector<double> k;

    for (auto &pd: vec_pd) {
        pairs.push_back(pd.residue_ids);
        k.push_back(pd.k);
    }

    FILE *rouxfile = fopen(input_files.roux_mdp.c_str(), "w");

    std::ifstream infile(input_files.mdp_template);
    std::string line;
    auto n_groups = pull_coord.size();
    auto n_coords = pairs.size();

    while (std::getline(infile, line)) {
        fprintf(rouxfile, "%s\n", line.c_str());
        if (line.find("pull-nstfout") != std::string::npos) {
            fprintf(rouxfile, "pull_ngroups = %lu \n", n_groups);
            fprintf(rouxfile, "pull_ncoords = %lu \n", n_coords);
            auto group_pairs = resis_to_groups(pairs);
            for (int i = 0; i < n_groups; ++i) {
                fprintf(rouxfile, "pull-group%i-name = %s\n", i + 1, pull_coord.at(i).c_str());
            }
            for (int i = 0; i < n_coords; ++i) {
                write_roux_pull_entry(group_pairs[i], rouxfile, i + 1, k[i]);
            }

        }
    }
    infile.close();
    fclose(rouxfile);
}