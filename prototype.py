#!/usr/bin/env python

import numpy as np
import copy
import time

beta = 1.0
interaction_strength = 2.0
interaction_cutoff = 2.5
n_particles = 256
density = 0.8
cell_length = ((n_particles*np.pi) / (6.0 * density))**(1/3.)
n_trials = 1e6
max_move = 0.5

def pbc_round(input):
     i = int(input)
     if (abs(input - i) >= 0.5):
         if (input > 0):
             i += 1
         if (input < 0):
             i -= 1
     return i

def dist(r1, r2):
    x = r1[0] - r2[0]
    y = r1[1] - r2[1]
    z = r1[2] - r2[2]
    x -= pbc_round(x / cell_length) * cell_length
    y -= pbc_round(y / cell_length) * cell_length
    z -= pbc_round(z / cell_length) * cell_length
    return np.sqrt(x * x + y * y + z * z)

def set_cell():
    particles = np.zeros((n_particles, 3))
    for i, p in enumerate(particles):
        keep_trying = True
        counter = 0
        while counter < n_particles - 1:
            x = np.random.random() * cell_length
            y = np.random.random() * cell_length
            z = np.random.random() * cell_length
            counter = 0
            for j, p2 in enumerate(particles):
                r = dist(np.array([x, y ,z]), p2)
                if r > 0.75 and i != j:
                    counter += 1

        particles[i][0] = x
        particles[i][1] = y
        particles[i][2] = z
    return particles

def energy(dr):
    if dr < interaction_cutoff:
        return 4 * interaction_strength * ((1 / dr)**12 - (1/dr)**6)
    return 0.

def singleAtomUpdate(coords, index):
    toten = 0.
    for i, c in enumerate(coords):
        if i != index:
            dr = dist(c, coords[index])
            toten += energy(dr)
    return toten
    
def totalEnergy(coords):
    toten = 0.
    for i in range(n_particles):
        for j in range(n_particles):
            if i != j:
                dr = dist(coords[i], coords[j])
                toten += energy(dr)
    return toten * 0.5

def main():
       coords = set_cell()
       print "Initialized"
       print "Cell size:", cell_length
       toten = totalEnergy(coords)
       counter = 0
       n_accepts = 0
       f = open('trajectory.xyz', 'w')
       start_time = time.clock()
       while counter < n_trials:
           # print "try",  counter
           atom = np.random.randint(n_particles)
           test_coords = copy.copy(coords)
           test_coords[atom][0] += (2.0 * np.random.random() - 1) * max_move
           test_coords[atom][1] += (2.0 * np.random.random() - 1) * max_move
           test_coords[atom][2] += (2.0 * np.random.random() - 1) * max_move
           next_en = toten + (singleAtomUpdate(test_coords, atom) - singleAtomUpdate(coords, atom))
           delta_e = next_en - toten
           delta_e_beta = delta_e * beta
           if delta_e_beta < 75.:
               if delta_e_beta < 0.:
                   # print "accepted"
                   toten = next_en
                   coords = copy.copy(test_coords)
                   n_accepts += 1
               elif np.exp( - delta_e_beta ) > np.random.random():
                   # print "accepted"
                   toten = next_en
                   coords = copy.copy(test_coords)
                   n_accepts += 1
           if counter % 1000 == 0:
               print "Step:", counter, "Accepts: ", n_accepts,"Energy:", toten / n_particles, "Elapsed time:", time.clock() - start_time

               f.write(str(n_particles) + "\n\n")
               for c in coords:
                   f.write("0  %f   %f   %f\n" % (c[0], c[1], c[2]))

           counter += 1
       f.close()


if __name__ == '__main__':
    main()