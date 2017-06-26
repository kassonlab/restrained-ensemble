#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include "ensemble.h"
#include "rouxfileio.h"

namespace mpi = boost::mpi;

int main(int argc, char **argv) {
    mpi::environment env(argc, argv);
    mpi::communicator world;
    int rank = world.rank();
    int root{0}, n{1};
    std::string config_filename;
    bool grompp{false}, check_forces{false};

    if (rank == root) {
        using namespace boost::program_options;
        options_description desc{
                "\n\t\t===== ROUX ENSEMBLE SIMULATOR =====\n"
                        "\tThis software package will perform a modified version\n"
                        "\tof Roux's restrained-ensemble simulations \n"
                        "\t(http://pubs.acs.org/doi/abs/10.1021/jp3110369)\n"
                        "\t\t===================================\n"
        };
        desc.add_options()("help,h", "Help screen")(
                "num,n", value<int>()->default_value(1),
                "Number of mdrun iterations to perform")
                ("config,f",
                 value<std::string>()->default_value("roux.ini"),
                 "Configuration file")
                ("grompp,g", "Only run grompp. "
                        "This will generate tprs, but will not call mdrun: "
                        "if you use this option, be aware that -n and -c become meaningless")
                ("checkf,c", "Check forces?");

        variables_map vm;
        store(parse_command_line(argc, const_cast<const char **>(argv), desc), vm);

        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n" << desc << std::endl;
            return EXIT_SUCCESS;
        }

        config_filename = vm.count("config") ? vm["config"].as<std::string>() :
                          "roux.ini";
        n = vm.count("num") ? vm["num"].as<int>() : 1;
        grompp = (bool) vm.count("grompp");
        check_forces = (bool) vm.count("checkf");
    }

    mpi::broadcast(world, config_filename, root);
    mpi::broadcast(world, n, root);
    mpi::broadcast(world, grompp, root);
    mpi::broadcast(world, check_forces, root);

    Ensemble ensemble(config_filename.c_str(), world);
    ensemble.input_names.differences = "/home/jmh/Research/test-roux-c/difference-files/testdiff.csv";

//    if (rank == 0) {
        read_exp_json(ensemble.input_names.exp_filename, ensemble.vec_pd);
//    }
//    mpi::broadcast(world, ensemble.vec_pd, 0);

    ensemble.do_histogram(5);
    if (grompp){
        ensemble.do_mdp(5);
        ensemble.do_grompp(5);
    }
    ensemble.do_mdrun();


    return EXIT_SUCCESS;
}