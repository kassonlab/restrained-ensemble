//
// Created by jmh on 6/23/17.
//

#include <boost/filesystem.hpp>
#include "ensemble.h"
#include "rouxalgorithms.h"
#include "rouxfileio.h"
#include "mdrun_main.h"

namespace mpi = boost::mpi;

Ensemble::Ensemble(const char *ini_filename, boost::mpi::communicator &comm) {
    parse_ini(ini_filename, vec_pd, input_names, prefs, params);
    world = comm;
}

void Ensemble::do_histogram(int ensemble_number) {
    int rank = world.rank();
    vecofstrings xvgs;

    if (rank == 0) {
        if (input_names.differences.empty()) {
            throw std::invalid_argument("You have not set the environment variable HISTDIF");
        }
        if (ensemble_number > 1) {
            for (int part = 0; part < params.boxcar_parts; part++) {

                auto vec_names = generate_gromacs_filenames(ensemble_number - part,
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

void Ensemble::do_mdp(int ensemble_number) {
    int rank{world.rank()};
    std::vector<std::string> pull_coord;
    if (rank >= params.num_replicas) {
        return;
    }
    std::string selection_filename = prefs.ensemble_path + "/select_roux_groups.dat";
    int replica = params.replicas[rank];
    if (rank == 0) {
        pull_coord = make_dat(selection_filename.c_str(), vec_pd, params);
    }
    mpi::broadcast(world, pull_coord, 0);
    auto names = generate_gromacs_filenames(ensemble_number, prefs, replica, true, false);

    generate_ndx_files(input_names.gmx_exe, names, selection_filename, false);

    make_mdp(vec_pd, input_names, params, pull_coord);
}

void Ensemble::do_grompp(int ensemble_number) {
    int rank{world.rank()};
    if (rank >= params.num_replicas) {
        return;
    }
    int replica = params.replicas[rank];
    auto name = generate_gromacs_filenames(ensemble_number, prefs, replica, false, false);

    char buffer[BUFFER_LENGTH];
    snprintf(buffer, BUFFER_LENGTH,
             "%s grompp -f %s -c %s -p %s/tops/topol%02i.top -o %s -n %s -po %s/%s_%02i/rank.mdp",
             input_names.gmx_exe.c_str(),
             input_names.roux_mdp.c_str(),
             name.gro.c_str(),
             prefs.ensemble_path.c_str(),
             replica,
             name.tpr.c_str(),
             name.ndx.c_str(),
             prefs.ensemble_path.c_str(),
             prefs.directory_prefix.c_str(),
             replica);
    system(buffer);
}

void Ensemble::do_mdrun() {
    const char **argv;
    std::vector<std::string> argv_string;
    int argc;
    argv_string = mdrun_chararray(params, prefs, argc);

    std::vector<const char *> v(argv_string.size());
    for (int i = 0; i < v.size(); ++i) {
        v[i] = argv_string[i].c_str();
    }
    v.push_back(NULL);

    argv = &v[0];

    gmx_mdrun(argc, const_cast<char **>(argv));
}