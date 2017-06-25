//
// Created by jmh on 6/23/17.
//

#include <boost/filesystem.hpp>
#include "ensemble.h"
#include "rouxalgorithms.h"
#include "rouxfileio.h"

namespace mpi = boost::mpi;

Ensemble::Ensemble(const char * ini_filename, boost::mpi::communicator &comm){
    parse_ini(ini_filename, vec_pd, input_names, prefs, params);
    world = comm;
}

void Ensemble::do_histogram(int ensemble_number) {
    int rank = world.rank();
    vecofstrings xvgs;

    if (rank == 0) {
        if (input_names.differences.empty()){
            throw std::invalid_argument("You have not set the environment variable HISTDIF");
        }
        if (ensemble_number > 1) {
            for (int part = 0; part < params.boxcar_parts; part++) {

                auto vec_names = generate_gromacs_filenames(ensemble_number-part,
                                                            prefs, params.replicas,
                                                            false, false);
                for (auto &name: vec_names) {
                    if (name.xtc.empty()) break;
                    if (boost::filesystem::exists(name.xvg)) {
                        xvgs.push_back(name.xvg);
                    }
                }
            }
        } else {
            std::vector<gromacs_files> vec_names = generate_gromacs_filenames(
                    ensemble_number, prefs, params.replicas, true, true);
            for (auto &names: vec_names) {
                xvgs.push_back(names.xvg);
            }
        }
    }
    mpi::broadcast(world, xvgs, 0);
    mpi_read_xvgs(world, vec_pd, xvgs, params.num_pairs, true);
    world.barrier();
//    std::cout << vec_pd[0].sim_dist_data.at(2) << std::endl;
    if (rank == 0) {
        calculate_histogram(vec_pd,
                            input_names.differences.c_str(),
                            params,
                            ensemble_number);
    }
}