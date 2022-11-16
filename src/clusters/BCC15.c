#include "simple_cluster_methods.h"
#include "BCC15.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

//!  An BCC15 cluster is an arrangement of six sp4c clusters.
/*!
*  Find BCC15 clusters
*  A BCC15 is an arrangement of six sp4c clusters.
*
*  Cluster output: SBBBBBBBBBBBB
*  Storage order: central particle, other particles
*/

void Clusters_GetBCC_15() {
    int i, j, j2, k, l, m;
    int flg;
    int s_com=-1;
    int trial[15];
    int clusSize=15;

    for (i=0; i < nsp4c - 1; i++) { // loop over sp4c
        for (j2=4; j2<6; j2++) { // loop over spindles
            int central_particle = hcsp4c[i][j2];
            int *first_sp4c_cluster = hcsp4c[i];

            // check that central spindle is part of at least six sp4c clusters
            if (nmem_sp4c[central_particle] < 6) continue;

            int sp4c_ids[6];
            int idx_sp4c = 0;
            sp4c_ids[idx_sp4c] = i;
            idx_sp4c++;

             // loop over sp4c's the central particle is part of
            flg = 0;
            for (j=0; j < nmem_sp4c[central_particle]; ++j) {
                int other_sp4c_id = mem_sp4c[central_particle][j];
                if (other_sp4c_id <= i) continue;
                int *other_sp4c_cluster = hcsp4c[other_sp4c_id];

                // make sure central particle is spindle of other sp4c
                if (central_particle != other_sp4c_cluster[4] && central_particle != other_sp4c_cluster[5]) continue;

                // count common particles
                int ncom = 0;
                for (k=0; k<6; k++) {
                    for (l=0; l<6; l++) {
                        if (first_sp4c_cluster[k] == other_sp4c_cluster[l]) {
                            ncom++;
                        }
                    }
                }
                if (ncom!=1 || ncom!=3) continue;

                // count common ring particles
                int ncom_ring = 0;
                for (k=0; k<4; k++) {
                    for (l=0; l<4; l++) {
                        if (first_sp4c_cluster[k] == other_sp4c_cluster[l]) {
                            ncom_ring++;
                        }
                    }
                }
                if (ncom-1 != ncom_ring) continue;

                // add to possible sp4c's
                sp4c_ids[idx_sp4c] = i;
                idx_sp4c++;
                if (idx_sp4c == 6) {
                    flg = 1;
                    break;
                }
            }

            if (flg==1) continue;


            for (j=0; j < nmem_sp4c[central_particle]; ++j) {
                // sort particle ids to check uniqueness
                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]= first_sp4c_cluster[k];
                    trial[k+5]= hcsp4c[mem_sp4c[central_particle][j]][k];
                }
                quickSort(&trial[1],14);

                flg=0;  // check trial cluster not already found
                for (k=0; k < nBCC_15; ++k) {
                    for (l=0; l<15; ++l) {
                        if (trial[l] != hcBCC_15[k][l]) break;
                    }
                    if (l==15) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_15 == mBCC_15) {
                    hcBCC_15= resize_2D_int(hcBCC_15, mBCC_15, mBCC_15 + incrStatic, clusSize, -1);
                    mBCC_15= mBCC_15 + incrStatic;
                }
                for (k=0; k<15; ++k) hcBCC_15[nBCC_15][k]=trial[k];

                Cluster_Write_BCC15();
            }
        }
    }
}

int count_shared_sp4cs(int *neighbouring_sp4c_ids, const int first_sp4c_id, const int center_id) {
    // Count how many sp4c's have a spindle in common with first_sp4c_id and get their ids
    int num_shared_sp5b = 0;
    for (int other_sp5b_pointer = 0; other_sp5b_pointer < nmem_sp5b[center_id]; ++other_sp5b_pointer) {
        int other_sp5_id = mem_sp5b[center_id][other_sp5b_pointer];
        if (other_sp5_id > first_sp5b_id) {
            if (num_shared_sp5b < 9) {
                neighbouring_sp5_ids[num_shared_sp5b] = other_sp5_id;
            }
            num_shared_sp5b++;
        }
    }
    return num_shared_sp5b;
}


void Cluster_Write_BCC15() {
    int i;
    sBCC_15[hcBCC_15[nBCC_15][0]] = 'S';
    for (i = 1; i< 15; i++){
        if (sBCC_15[hcBCC_15[nBCC_15][i]] == 'C') sBCC_15[hcBCC_15[nBCC_15][i]] = 'B';
    }

    ++nBCC_15;
}