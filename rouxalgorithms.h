//
// Created by jmh on 6/23/17.
//

#ifndef ROUX_ROUXALGORITHMS_H
#define ROUX_ROUXALGORITHMS_H

void calculate_histogram(std::vector<pair_data> vec_pd,
                         const char *out_filename,
                         parameters params,
                         int ensemble_number);

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


#endif //ROUX_ROUXALGORITHMS_H
