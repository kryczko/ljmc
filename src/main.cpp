#include <iostream>
#include <fstream>
#include <cstdlib>

#include "config.h"
#include "helper.h"

int main(int argc, char** argv) {
    Config config;
    std::vector<std::vector<double>> coords = set_cell(config);
    double toten = totalEnergy(coords, config);
    int counter = 0;
    int accepts = 0;
    std::ofstream traj("trajectory.xyz");


    double next_en, delta_e, delta_e_beta;
    while (counter < config.n_trials) {
        for (int i = 0; i < config.n_particles; i ++) {
            std::vector<std::vector<double>> test_coords = coords;
            test_coords[i][0] += (2.0 * config.rng.random() - 1) * config.max_move;
            test_coords[i][1] += (2.0 * config.rng.random() - 1) * config.max_move;
            test_coords[i][2] += (2.0 * config.rng.random() - 1) * config.max_move;

            next_en = toten + (singleAtomUpdate(test_coords, i, config) - singleAtomUpdate(coords, i, config));
            delta_e = next_en - toten;
            delta_e_beta = delta_e * config.beta;
            if ( delta_e_beta < 75.0) {
                if (delta_e_beta < 0.0) {
                    toten = next_en;
                    coords = test_coords;
                    accepts ++;
                } else if (std::exp(- delta_e_beta) > config.rng.random()) {
                    toten = next_en;
                    coords = test_coords;
                    accepts ++;
                }
            }
            counter ++;
            if (counter % 1000 == 0) {
                toten = totalEnergy(coords, config);
                std::cout << "Step: " << counter << "  Accepts: " << accepts << " Energy: " << toten  / config.n_particles << std::endl;
                traj << config.n_particles << "\n\n";
                for (int i = 0; i < config.n_particles; i ++) {
                    traj << "0  " << coords[i][0] << "  " << coords[i][1] << "  " << coords[i][2] << std::endl;  
                }
            }
        }
    }
    traj.close();
    std::string restart = "tail -n " + std::to_string(config.n_particles + 2) + " trajectory.xyz > restart.xyz";
    system(restart.c_str());
    return 0;
}