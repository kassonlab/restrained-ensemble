#include <boost/mpi.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "ensemble.h"
#include "rouxfileio.h"

namespace mpi = boost::mpi;

int main(int argc, char **argv) {
    mpi::environment env(argc, argv);
    mpi::communicator world;

    Ensemble ensemble("/home/jmh/Research/test-roux-c/test_opa_opa.ini", world);
    ensemble.input_names.differences = "/home/jmh/Research/test-roux-c/difference-files/testdiff.csv";

    int rank = world.rank();
    if (rank == 0) {
        read_exp_json(ensemble.input_names.exp_filename, ensemble.vec_pd);
    }
    mpi::broadcast(world, ensemble.vec_pd, 0);

    ensemble.do_histogram(5);

    return EXIT_SUCCESS;
}