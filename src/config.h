#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <cmath>
#include "rng.h"

struct Config {
    double beta = 1.0;
    double epsilon = 2.0;
    double cutoff = 2.5;
    double density = 0.2;
    double cell_len = 10.0;
    double max_move = 0.2;

    int n_particles = 1000;
    int n_trials = 1e6;

    Rng rng;


    Config() {
        this->cell_len = std::pow((this->n_particles * M_PI) / (6.0 * this->density), 1.0 / 3.0);
    }
};

#endif