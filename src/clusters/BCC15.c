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
    int count1,count2,count3,count4;
    count1=count2=count3=count4=0;

    for (int central_particle = 0; central_particle < particles_in_current_frame; central_particle++) {

        // check that central spindle is part of at least six sp4c clusters
        if (nmem_sp4c[central_particle] < 6) continue;

        // check that central particle is a spindle of exactly six sp4c clusters
        int sp4c_ids[6];
        int num_sp4cs = 0;

        for (j=0; j < nmem_sp4c[central_particle]; ++j) {
            int sp4c_id = mem_sp4c[central_particle][j];
            int *sp4c_cluster = hcsp4c[sp4c_id];

            if (central_particle == sp4c_cluster[4] || central_particle == sp4c_cluster[5]) {
                // add to possible sp4c's
                if (num_sp4cs < 6) {
                    sp4c_ids[num_sp4cs] = sp4c_id;
                }
                num_sp4cs++;
            }
        }
        if (num_sp4cs != 6) {
            count1++;
            continue;
        }

        // check common particles between sp4c clusters
        flg = 0;
        for (i=0; i<6; ++i) {
            int ncom1 = 0;
            int ncom3 = 0;
            for (j=0; j<6; ++j) {
                if (i == j) continue;
                int *sp4c_cluster = hcsp4c[sp4c_ids[i]];
                int *sp4c_cluster_ = hcsp4c[sp4c_ids[j]];

                // count common particles: should be 1 or 3
                int ncom = 0;
                for (k=0; k<6; k++) {
                    for (l=0; l<6; l++) {
                        if (sp4c_cluster[k] == sp4c_cluster_[l]) {
                            ncom++;
                        }
                    }
                }
                if (ncom!=1 && ncom!=3) {
                    flg = 1;
                    count2++;
                    break;
                }

                // count common ring particles: should be ncom-1
                int ncom_ring = 0;
                for (k=0; k<4; k++) {
                    for (l=0; l<4; l++) {
                        if (sp4c_cluster[k] == sp4c_cluster_[l]) {
                            ncom_ring++;
                        }
                    }
                }
                if (ncom-1 != ncom_ring) {
                    count3++;
                    flg = 1;
                    break;
                }

                // count number of sp4cs with one common particle: should be 1
                if (ncom == 1) {
                    ncom1++;
                }
            }
            if (ncom1 != 1) {
                count4++;
                flg = 1;
                break;
            }
            if (flg==1) break;
        }
        if (flg == 1) continue;

        // TODO: output right particle ID's
        trial[0] = central_particle;
        for (i=1; i<15; i++) {
            trial[i] = central_particle;
        }

        if (nBCC_15 == mBCC_15) {
            hcBCC_15= resize_2D_int(hcBCC_15, mBCC_15, mBCC_15 + incrStatic, clusSize, -1);
            mBCC_15= mBCC_15 + incrStatic;
        }
        for (k=0; k<15; ++k) hcBCC_15[nBCC_15][k]=trial[k];

        Cluster_Write_BCC15();
    }
    // printf("%d %d %d %d", count1, count2, count3, count4);
}


void Cluster_Write_BCC15() {
    int i;
    sBCC_15[hcBCC_15[nBCC_15][0]] = 'S';
    for (i = 1; i< 15; i++){
        if (sBCC_15[hcBCC_15[nBCC_15][i]] == 'C') sBCC_15[hcBCC_15[nBCC_15][i]] = 'B';
    }
    ++nBCC_15;
}