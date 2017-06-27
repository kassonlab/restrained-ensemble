//
// Created by jmh on 6/26/17.
//

#include "logging.h"
#include "rouxfileio.h"
#include <boost/format.hpp>

using boost::format;
using boost::str;
namespace mpi=boost::mpi;

void Logging::read_summary_data(int part) {
    std::vector<std::string> dist_files, force_files;
    std::string hist_difs;
    hist_difs = str(format("%s.part%04i.csv") %
                                    input_names.differences.substr(0, input_names.differences.length() - 4) %
                                    part);

    for (auto ensemble_member: params.replicas) {
        dist_files.push_back(str(format("%s/%s_%02i/%s_pullx.part%04i.xvg")
                                 % prefs.ensemble_path
                                 % prefs.directory_prefix
                                 % ensemble_member
                                 % prefs.prod_prefix
                                 % part));
        force_files.push_back(str(format("%s/%s_%02i/%s_pullf.part%04i.xvg")
                                  % prefs.ensemble_path
                                  % prefs.directory_prefix
                                  % ensemble_member
                                  % prefs.prod_prefix
                                  % part));
    }
    mpi_read_xvgs(world, vec_sd, dist_files, force_files, params.num_pairs);
    if (world.rank() == 0) read_histograms(hist_difs, vec_sd);
    mpi::broadcast(world, vec_sd, 0);
}

void Logging::calculate_forces() {
    for (auto &sd: vec_sd) {
        auto n_samples = sd.sim_dist_data.size();
        sd.calc_forc_data.resize(n_samples);
        double force, x, exponential;
        auto norm = sd.k / (pow(params.sigma, 3) * sqrt(2 * M_PI));
        for (int sample = 0; sample < n_samples; ++sample) {
            force = 0;
            for (int bin_num = 0; bin_num < params.num_bins; ++bin_num) {
                x = bin_num * params.bin_width - sd.sim_dist_data[sample];
                exponential = exp(-pow(x, 2) / (2 * pow(params.sigma, 2)));
                force -= sd.hist_difference[bin_num] * x * exponential;
            }
            sd.calc_forc_data[sample] = norm * force;
        }
    }
}

void Logging::write_header(std::ofstream &outfile) {

    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    char timebuf[BUFFER_LENGTH], timedata[BUFFER_LENGTH], inputdata[BUFFER_LENGTH], head[BUFFER_LENGTH];

    std::strftime(timebuf, BUFFER_LENGTH, "# Log file opened on %D at %I:%M%p\n", timeinfo);

    std::string pair_string, k_string;
    for (auto &sd: vec_sd) {
        pair_string += str(format("%03i %03i ") % sd.residue_ids.first % sd.residue_ids.second);
        k_string += str(format("%f ") % sd.k);
    }
    snprintf(inputdata, BUFFER_LENGTH,
             " Input parameters:\n@ experimental-filename:\t%s\n@ mdp:\t\t\t\t%s\n@ num-bins:"
                     "\t\t\t%lu\n@ bin-width:\t\t\t%.2f\n@ sigma:\t\t\t%.2f\n@ boxcar-parts:\t\t\t%lu\n@ "
                     "pairs:\t\t\t%s\n@ k:\t\t\t\t%s\n@ out-filename:\t\t\t%s",
             input_names.exp_filename.c_str(),
             input_names.roux_mdp.c_str(),
             params.num_bins,
             params.bin_width,
             params.sigma,
             params.boxcar_parts,
             pair_string.c_str(),
             k_string.c_str(),
             input_names.differences.c_str());

    snprintf(head, BUFFER_LENGTH, "\n#\n@ Time\t(X) %s\t(F) %s\t(CF) %s\n",
             pair_string.c_str(), pair_string.c_str(), pair_string.c_str());

    outfile << timebuf << "\n#\n#\tR O U X  R E S T R A I N E D - E N S E M B L E\n#\tsimulation logfile\n"
            "#\timplemented by J.M. Hays\n#\n#" << timedata << inputdata << head;
}

void Logging::write_summary(int part, bool check_forces) {
//    simdata dists, time_and_forces, dummy, histogram_difference, forces_check;
//    printf("INFO: Writing summary to file\n");

    read_summary_data(part);
    auto num_samples = vec_sd[0].sim_time_data.size();

    if (check_forces) {
        calculate_forces();
    }

    char buffer[BUFFER_LENGTH];
    for (auto &sd: vec_sd) {
        std::snprintf(buffer, BUFFER_LENGTH, "%s/pair%03i_%03i.part%04i.txt",
                      input_names.log_dir.c_str(),
                      sd.residue_ids.first,
                      sd.residue_ids.second,
                      part);
        sd.summary_filename = buffer;
        std::ofstream summary_file(sd.summary_filename);
        write_header(summary_file);
        for (int i = 0; i < num_samples; ++i) {
            summary_file.precision(4);
            summary_file << sd.sim_time_data[i] << "\t";
            summary_file.precision(6);
            summary_file << sd.sim_dist_data[i] << "\t";
            summary_file << sd.sim_forc_data[i] << "\t"; // j+1 since first column is time

            if (check_forces) {
                summary_file << sd.calc_forc_data[i] << "\t";
            }
            summary_file << "\n";
        }
        summary_file.close();
    }
//    printf("INFO: Finished writing summary to file\n");
}
