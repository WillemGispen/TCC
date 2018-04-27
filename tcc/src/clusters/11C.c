#include <clusters/simple_cluster_methods.h>
#include "globals.h"
#include "bonds.h"
#include "tools.h"
#include "11C.h"

void Clusters_Get11C() {
    //!  An 11C cluster is the intersection of two 7A clusters with a common spindle
    /*!
   *  Find 11C clusters
   *  An 11C is two 7A clusters where:
   *      - There are two common particles between the two sp5 rings. These are a bonded pair.
   *      - There are two more bonds between two pairs of distinct particles in the sp5 ring.
   *
   *  Cluster output: SOOBBBBBBBB
   *  Storage order: common_spindle x 1, uncommon_spindles x 2, common_ring_particles x 2,
   */

    int common_ring_particles[5],uncommon_spindle[2];
    int k, l, m;
    int common_spindle[2];
    int break_out;

    int second_7A_pointer;

    for (int first_7A_id = 0; first_7A_id < nsp5c - 1; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for(int first_spindle_pointer = 5; first_spindle_pointer < 7; first_spindle_pointer++) {
            int first_spindle_id = first_7A_cluster[first_spindle_pointer];
            for (second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[first_spindle_id]; second_7A_pointer++) {
                int second_7A_id = mem_sp5c[first_spindle_id][second_7A_pointer];
                if(second_7A_id > first_7A_id) {
                    int *second_7A_cluster = hcsp5c[second_7A_id];

                    if(count_common_spindle_particles(first_7A_cluster, second_7A_cluster, 7, 7, common_spindle) == 1) {

                        uncommon_spindle[0] = get_uncommon_spindle(first_7A_cluster, 7, common_spindle[0]);
                        uncommon_spindle[1] = get_uncommon_spindle(second_7A_cluster, 7, common_spindle[0]);

                        // need two common particles from SP5 rings
                        if (count_common_ring_particles(first_7A_cluster, second_7A_cluster, 5, 5, common_ring_particles) == 2) {

                            // the common ring particles need to be bonded
                            if (Bonds_BondCheck(common_ring_particles[0], common_ring_particles[1]) == 1) {

                                // two bonds between non-common SP5 ring particles
                                if (count_bonded_ring_particles_11C(common_ring_particles, first_7A_cluster, second_7A_cluster) == 2) {

                                    int clusSize = 11;

                                    if (n11C == m11C) {
                                        hc11C = resize_2D_int(hc11C, m11C, m11C + incrStatic, clusSize, -1);
                                        m11C = m11C + incrStatic;
                                    }

                                    // hc11C key: (s_com, s_i, s_j, r_ca, r_cb, d_i, d_i, d_j, d_j, unc_i, unc_j)
                                    hc11C[n11C][0] = common_spindle[0];
                                    hc11C[n11C][1] = uncommon_spindle[0];
                                    hc11C[n11C][2] = uncommon_spindle[1];
                                    hc11C[n11C][3] = common_ring_particles[0];
                                    hc11C[n11C][4] = common_ring_particles[1];

                                    int first_bonded_particles[3];
                                    int second_bonded_particles[3];
                                    if (count_particles_bonded_to_common(first_7A_cluster, common_ring_particles, first_bonded_particles) == 2) {
                                        if (count_particles_bonded_to_common(second_7A_cluster, common_ring_particles, second_bonded_particles) == 2) {

                                            hc11C[n11C][5] = first_bonded_particles[0];
                                            hc11C[n11C][6] = first_bonded_particles[1];
                                            hc11C[n11C][7] = second_bonded_particles[0];
                                            hc11C[n11C][8] = second_bonded_particles[1];

                                            // Check that the bonded non-common particles are bonded
                                            if (Bonds_BondCheck(hc11C[n11C][5], hc11C[n11C][7]) == 0) continue;
                                            if (Bonds_BondCheck(hc11C[n11C][6], hc11C[n11C][8]) == 0) continue;

                                            // Get the ID's of the non-common particles
                                            for (k = 0; k < 5; ++k) {
                                                if (Bonds_BondCheck(first_7A_cluster[k], hc11C[n11C][5]) &&
                                                    Bonds_BondCheck(first_7A_cluster[k], hc11C[n11C][6])) {
                                                    hc11C[n11C][9] = first_7A_cluster[k];
                                                }
                                                if (Bonds_BondCheck(second_7A_cluster[k], hc11C[n11C][7]) &&
                                                    Bonds_BondCheck(second_7A_cluster[k], hc11C[n11C][8])) {
                                                    hc11C[n11C][10] = second_7A_cluster[k];
                                                }
                                            }
                                            quickSort(&hc11C[n11C][1], 2);
                                            quickSort(&hc11C[n11C][3], 2);
                                            quickSort(&hc11C[n11C][5], 4);
                                            quickSort(&hc11C[n11C][9], 2);

                                            Cluster_Write_11C();

                                            ++n11C;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

int count_particles_bonded_to_common(const int *cluster, const int *common_particles, int *bonded_particles) {

    int num_bonds = 0;
    for (int ring_pointer = 0; ring_pointer < 5; ++ring_pointer) {
        int ring_particle_id = cluster[ring_pointer];
        if (is_particle_in_cluster(common_particles, 2, ring_particle_id) == 0) {
            if (Bonds_BondCheck(ring_particle_id, common_particles[0]) == 1) {
                bonded_particles[0] = ring_particle_id;
                num_bonds++;
            } else if (Bonds_BondCheck(ring_particle_id, common_particles[1]) == 1) {
                bonded_particles[1] = ring_particle_id;
                num_bonds++;
            }
        }
    }
    return num_bonds;
}

int count_bonded_ring_particles_11C(const int *common_ring_particles, const int *first_7A_cluster, const int *second_7A_cluster) {
    int num_bonded_ring_particles = 0;
    for (int first_ring_pointer = 0; first_ring_pointer < 5; ++first_ring_pointer) {
        int first_particle_id = first_7A_cluster[first_ring_pointer];

        if (is_particle_in_cluster(common_ring_particles, 2, first_particle_id) == 0) {
            for (int second_ring_pointer = 0; second_ring_pointer < 5; ++second_ring_pointer) {
                int second_particle_id = second_7A_cluster[second_ring_pointer];

                if (is_particle_in_cluster(common_ring_particles, 2, second_particle_id) == 0) {
                    if (Bonds_BondCheck(first_particle_id, second_particle_id)) {
                        ++num_bonded_ring_particles;
                    }
                }
            }
        }
    }
    return num_bonded_ring_particles;
}

void Cluster_Write_11C() {
    int i;

    s11C[hc11C[n11C][0]] = 'S';
    if(s11C[hc11C[n11C][1]] != 'S') s11C[hc11C[n11C][1]] = 'O';
    if(s11C[hc11C[n11C][2]] != 'S') s11C[hc11C[n11C][2]] = 'O';
    for(i=3; i< 11; i++) {
        if (s11C[hc11C[n11C][i]] == 'C') s11C[hc11C[n11C][i]] = 'B';
    }
}