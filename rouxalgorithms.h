//
// Created by jmh on 6/23/17.
//

#ifndef ROUX_ROUXALGORITHMS_H
#define ROUX_ROUXALGORITHMS_H

std::tuple<std::vector<int>, std::vector<int>> unique_from_pairs(std::vector<std::pair<int, int>> pairs);

void calculate_histogram(std::vector<pair_data> vec_pd,
                         const char *out_filename,
                         parameters params,
                         int ensemble_number);


int find_last_run_number(std::string ensemble_path,
                         std::string directory_prefix,
                         int replica);

void pre_process(gromacs_files name,
                 int chains,
                 bool aa,
                 std::vector<std::pair<int, int>> pairs,
                 std::string gmx_exe,
                 bool rewrite);

gromacs_files generate_gromacs_filenames(int ensemble_number,
                                prefixes prefs,
                                int replica,
                                bool start_tpr,
                                bool protein_ndx);

std::vector<gromacs_files> generate_gromacs_filenames(int ensemble_number,
                                                      prefixes prefs,
                                                      std::vector<int> replicas,
                                                      bool start_tpr,
                                                      bool protein_ndx);

std::vector<std::string> mdrun_chararray(parameters params, prefixes prefs, int &argc);


#endif //ROUX_ROUXALGORITHMS_H
