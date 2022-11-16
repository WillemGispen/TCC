#include "simple_cluster_methods.h"
#include "BCC15.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

//!  An BCC15 cluster is an BCC9 cluster with four extra particles attached.
/*!
*  Find BCC15 clusters
*  A BCC15 is a BCC9 cluster with four extra particles attached.
*  Each extra particle is bonded to four of the non-spindle particles of the BCC9 cluster.
*
*  Cluster output: SBBBBBBBBEEEE
*  Storage order: as for BCC9 x 9, extra_particles)
*/
void Clusters_GetBCC_15() {
    for (int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for (int spindle_pointer = 0; spindle_pointer < 2; ++spindle_pointer) {
            int primary_spindle = first_7A_cluster[5 + spindle_pointer];

            for (int new_particle_pointer = 0; new_particle_pointer < num_bonds[primary_spindle]; ++new_particle_pointer) {
                int new_particle_id = bond_list[primary_spindle][new_particle_pointer];

                if(is_particle_in_cluster(first_7A_cluster, 7, new_particle_id) == 0) {

                    if (count_bonds_to_7A_ring(first_7A_id, new_particle_id) == 2) {

                        Cluster_Write_BCC15(first_7A_cluster, new_particle_id);
                    }
                }
            }
        }
    }
}

int count_bonds_to_7A_ring(int first_7A_id, int new_particle_id) {
    int num_bonds_to_ring = 0;
    for (int ring_pointer = 0; ring_pointer < 5; ++ring_pointer) {
        if (Bonds_BondCheck(new_particle_id, hcsp5c[first_7A_id][ring_pointer])) {
            ++num_bonds_to_ring;
        }
    }
    return num_bonds_to_ring;
}

void Cluster_Write_BCC15(int *first_7A_cluster, int new_particle_id) {
    int clusSize = 8;

    // Now we have found the BCC15 cluster
    if (nBCC_15 == mBCC_15) {
        hcBCC_15 = resize_2D_int(hcBCC_15, mBCC_15, mBCC_15 + incrStatic, clusSize, -1);
        mBCC_15 = mBCC_15 + incrStatic;
    }

    for (int i = 0; i < 7; i++) {
        hcBCC_15[nBCC_15][i] = first_7A_cluster[i];
    }
    hcBCC_15[nBCC_15][7] = new_particle_id;

    if (sBCC_15[hcBCC_15[nBCC_15][7]] == 'C') sBCC_15[hcBCC_15[nBCC_15][7]] = 'B';
    if (sBCC_15[hcBCC_15[nBCC_15][0]] == 'C') sBCC_15[hcBCC_15[nBCC_15][0]] = 'B';
    if (sBCC_15[hcBCC_15[nBCC_15][1]] == 'C') sBCC_15[hcBCC_15[nBCC_15][1]] = 'B';
    if (sBCC_15[hcBCC_15[nBCC_15][2]] == 'C') sBCC_15[hcBCC_15[nBCC_15][2]] = 'B';
    if (sBCC_15[hcBCC_15[nBCC_15][3]] == 'C') sBCC_15[hcBCC_15[nBCC_15][3]] = 'B';
    if (sBCC_15[hcBCC_15[nBCC_15][4]] == 'C') sBCC_15[hcBCC_15[nBCC_15][4]] = 'B';
    sBCC_15[hcBCC_15[nBCC_15][5]] = 'O';
    sBCC_15[hcBCC_15[nBCC_15][6]] = 'O';
    ++nBCC_15;
}