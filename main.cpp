
#include <boost/mpi.hpp>
#include <boost/filesystem.hpp>
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

        if (!boost::filesystem::exists(config_filename)) {
            char error[BUFFER_LENGTH];
            snprintf(error, BUFFER_LENGTH, "The configuration file %s does not exist: "
                    "please provide a valid configuration file.", config_filename.c_str());
            throw std::invalid_argument(error);
        }
        n = vm.count("num") ? vm["num"].as<int>() : 1;
        grompp = (bool) vm.count("grompp");
        check_forces = (bool) vm.count("checkf");
    }

    mpi::broadcast(world, config_filename, root);
    mpi::broadcast(world, n, root);
    mpi::broadcast(world, grompp, root);
    mpi::broadcast(world, check_forces, root);
    Ensemble ensemble(config_filename.c_str(), world);
    Logging logger;
    int ensemble_number;

    ensemble.input_names.differences = getenv("HISTDIF");
    read_exp_json(ensemble.input_names.exp_filename, ensemble.vec_pd);
    ensemble.link_to_logging(logger);
    ensemble_number = ensemble.setup_restart(logger, check_forces);
    mpi::broadcast(world, ensemble_number, root);

    ensemble.do_histogram(ensemble_number);

    if (grompp) {
        ensemble.do_mdp(ensemble_number);
        ensemble.do_grompp(ensemble_number);
    }

    int iter{0};
    while (iter < n) {
        ensemble.do_mdrun();
        ensemble.do_histogram(ensemble_number + 1);
        logger.write_summary(ensemble_number, check_forces);
        ++ensemble_number;
        ++iter;
    }

    return EXIT_SUCCESS;
}