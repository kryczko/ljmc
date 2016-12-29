#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "config.h"

int pbc_round(double in) {
    int i = in;
    if (std::abs(in - i) >= 0.5) {
        if (in > 0) {
            i += 1;
        }
        if (in < 0) {
            i -= 1;
        }
    }
    return i;
}

double dist(std::vector<double> r1, std::vector<double> r2, Config& c) {
    
    double x = r1[0] - r2[0];
    double y = r1[1] - r2[1];
    double z = r1[2] - r2[2];

    x -= pbc_round( x / c.cell_len ) * c.cell_len;
    y -= pbc_round( y / c.cell_len ) * c.cell_len;
    z -= pbc_round( z / c.cell_len ) * c.cell_len;

    return std::sqrt( x * x + y * y + z * z);
}

std::vector<std::vector<double>> set_cell(Config& c) {
    std::vector<std::vector<double>> particles(c.n_particles, {0.0, 0.0, 0.0});
    std::ifstream traj("restart.xyz");
    if (traj.good()) {
        std::cout << "Found restart file.\n";
        int n_parts;
        traj >> n_parts;
        if (n_parts != c.n_particles) {
            std::cout << "ERROR: Number of particles is not the same in restart xyz file.\n";
            std::exit(-1);
        }
        double atom, x, y, z;
        int counter = 0;
        while(counter < n_parts) {
            traj >> atom >> x >> y >> z;
            particles[counter][0] = x;
            particles[counter][1] = y;
            particles[counter][2] = z;
            counter ++;
        }
        return particles;
    }

    for (int i = 0; i < c.n_particles; i ++) {
        int counter = 0;
        double x,y,z;
        while (counter < c.n_particles - 1) {
            x = c.rng.random() * c.cell_len;
            y = c.rng.random() * c.cell_len;
            z = c.rng.random() * c.cell_len;
            counter = 0;
            for (int j = 0; j < c.n_particles; j ++) {
                double r = dist({x, y, z}, particles[j], c);
                if (r > 0.75 and i != j) {
                    counter ++;
                }
            }
        }
        particles[i][0] = x;
        particles[i][1] = y;
        particles[i][2] = z;
    }
    return particles;
}

double energy(double dr, Config& c) {
    if (dr < c.cutoff) {
        return 4 * c.epsilon * (std::pow((1 / dr), 12) - std::pow((1 / dr), 6));
    }
    return 0.0;
}

double singleAtomUpdate(std::vector<std::vector<double>> coords, int index, Config& c) {
    double toten = 0.0;
    for (int i = 0; i < coords.size(); i ++) {
        if (i != index) {
            double dr = dist(coords[i], coords[index], c);
            toten += energy(dr, c);
        }
    }
    return toten;
}

double totalEnergy(std::vector<std::vector<double>> coords, Config& c) {
    double toten = 0.0;
    for (int i = 0; i < coords.size(); i ++) {
        for (int j = 0; j < coords.size(); j ++) {
            if (i != j) {
                double dr = dist(coords[i], coords[j], c);
                toten += energy(dr, c);
            }
        }
    }
    return toten * 0.5;
}