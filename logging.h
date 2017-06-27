//
// Created by jmh on 6/26/17.
//

#ifndef ROUX_LOGGING_H
#define ROUX_LOGGING_H

#include <fstream>
#include <boost/mpi/communicator.hpp>
#include "common.h"

class Logging {
public:
    void read_summary_data(int part);
    void calculate_forces();
    void write_header(std::ofstream& outfile);
    void write_summary(int part, bool check_forces);

    boost::mpi::communicator world;
    parameters params;
    prefixes prefs;
    input_filenames input_names;
    std::vector<summary_data> vec_sd;
};


#endif //ROUX_LOGGING_H
