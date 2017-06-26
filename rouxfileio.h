//
// Created by jmh on 6/23/17.
//

#ifndef ROUX_ROUXFILEIO_H
#define ROUX_ROUXFILEIO_H

#include <boost/mpi.hpp>
#include "common.h"

void warn_file_exists(const char *filename);

void parse_ini(const char *ini_filename,
               std::vector<pair_data>& vec_pd,
               input_filenames& in_filenames,
               prefixes& prefs,
               parameters& pars);

void read_exp_json(std::string exp_filename, std::vector<pair_data> &pd);

void mpi_read_xvgs(boost::mpi::communicator &world,
                   setofpairs &vec_pd,
                   vecofstrings filenames,
                   unsigned long num_pairs,
                   bool skip_time);

vecofstrings make_dat(const char *dat_filename,
                      std::vector<pair_data> vec_pd,
                      parameters params);

void generate_ndx_files(std::string gmx_exe,
                        gromacs_files name,
                        std::string dat,
                        bool rewrite);

void make_mdp(std::vector<pair_data> vec_pd,
              input_filenames input_files,
              parameters params,
              vecofstrings pull_coord);

#endif //ROUX_ROUXFILEIO_H
