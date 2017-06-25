//
// Created by jmh on 6/23/17.
//

#ifndef ROUX_ENSEMBLE_H
#define ROUX_ENSEMBLE_H

#include <boost/mpi.hpp>
#include "common.h"

class Ensemble {
public:
    input_filenames input_names;
    prefixes prefs;
    parameters params;
    std::vector<pair_data> vec_pd;
    boost::mpi::communicator world;

    Ensemble(const char * ini_filename, boost::mpi::communicator& comm);
    void do_histogram(int ensemble_number);
    void do_mdp(int ensemble_number);
    void do_grompp(int ensemble_number);
    void do_mdrun();
    int setup_restart(bool check_forces=false);
};


#endif //ROUX_ENSEMBLE_H
