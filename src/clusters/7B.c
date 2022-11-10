#include "simple_cluster_methods.h"
#include "7B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

//!  An 7B cluster is an 6A cluster with an extra particle attached.
/*!
*  Find 7B clusters
*  An 7B is a 6A cluster with an extra particle bonded to two of the ring particles and a spindle of the 6A cluster
*
*  Cluster output: BBBBOOB
*  Storage order: as for 6A x 6, extra_particle)
*/
void Clusters_Get7B() {
    for (int first_6A_id = 0; first_6A_id < nsp4c; ++first_6A_id) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for (int spindle_pointer = 0; spindle_pointer < 2; ++spindle_pointer) {
            int primary_spindle = first_6A_cluster[4 + spindle_pointer];

            for (int new_particle_pointer = 0; new_particle_pointer < num_bonds[primary_spindle]; ++new_particle_pointer) {
                int new_particle_id = bond_list[primary_spindle][new_particle_pointer];

                if(is_particle_in_cluster(first_6A_cluster, 6, new_particle_id) == 0) {

                    if (count_bonds_to_6A_ring(first_6A_id, new_particle_id) == 2) {

                        Cluster_Write_7B(first_6A_cluster, new_particle_id);
                    }
                }
            }
        }
    }
}

int count_bonds_to_6A_ring(int first_6A_id, int new_particle_id) {
    int num_bonds_to_ring = 0;
    for (int ring_pointer = 0; ring_pointer < 4; ++ring_pointer) {
        if (Bonds_BondCheck(new_particle_id, hcsp4c[first_6A_id][ring_pointer])) {
            ++num_bonds_to_ring;
        }
    }
    return num_bonds_to_ring;
}

void Cluster_Write_7B(int *first_6A_cluster, int new_particle_id) {
    int clusSize = 7;

    // Now we have found the 7B Cs cluster
    if (n7B == m7B) {
        hc7B = resize_2D_int(hc7B, m7B, m7B + incrStatic, clusSize, -1);
        m7B = m7B + incrStatic;
    }

    for (int i = 0; i < 6; i++) {
        hc7B[n7B][i] = first_6A_cluster[i];
    }
    hc7B[n7B][7] = new_particle_id;

    if (s7B[hc7B[n7B][6]] == 'C') s7B[hc7B[n7B][6]] = 'B';
    if (s7B[hc7B[n7B][0]] == 'C') s7B[hc7B[n7B][0]] = 'B';
    if (s7B[hc7B[n7B][1]] == 'C') s7B[hc7B[n7B][1]] = 'B';
    if (s7B[hc7B[n7B][2]] == 'C') s7B[hc7B[n7B][2]] = 'B';
    if (s7B[hc7B[n7B][3]] == 'C') s7B[hc7B[n7B][3]] = 'B';
    s7B[hc7B[n7B][4]] = 'O';
    s7B[hc7B[n7B][5]] = 'O';
    ++n7B;
}