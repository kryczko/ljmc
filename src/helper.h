#ifndef _HELPER_H_
#define _HELPER_H_

int pbc_round(double in);
double dist(std::vector<double> r1, std::vector<double> r2, double cell_len);
std::vector<std::vector<double>> set_cell(Config& c);
double energy(double dr, double epsilon, double cutoff, Config&);
double singleAtomUpdate(std::vector<std::vector<double>> coords, int index, Config&);
double totalEnergy(std::vector<std::vector<double>> coords, Config&);

#endif