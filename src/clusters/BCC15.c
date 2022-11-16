#include "simple_cluster_methods.h"
#include "BCC15.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

//!  An BCC15 cluster is an BCC15 cluster with six extra particles attached.
/*!
*  Find BCC15 clusters
*  A BCC15 is a BCC15 cluster with six extra particles attached.
*  Each extra particle is bonded to four of the non-spindle particles of the BCC15 cluster.
*
*  Cluster output: SBBBBBBBBEEEE
*  Storage order: as for BCC15 x 9, extra_particles)
*/

void Clusters_GetBCC_15() {
    int i, j, j2, k, l, m;
    int flg;
    int s_com=-1;
    int trial[15];
    int clusSize=15;

    for (i=0; i < nsp4c - 1; i++) {
        for (j2=4; j2<6; j2++) {
            for (j=0; j < nmem_sp4c[hcsp4c[i][j2]]; ++j) { // loop over all sp3c_j
                if (mem_sp4c[hcsp4c[i][j2]][j] <= i) continue;
                m=0;
                for (k=4; k<6; k++) {
                    for (l=4; l<6; l++) {
                        if (hcsp4c[i][k] == hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l]) {
                            s_com= hcsp4c[i][k];
                            m++;
                        }
                    }
                }
                if (m==0 || m>1) continue;

                flg=0;
                for (k=0; k<6; k++) {
                    if (hcsp4c[i][k] == s_com) continue;
                    for (l=0; l<6; l++) {
                        if (hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l] == s_com) continue;
                        if (hcsp4c[i][k] == hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l]) {
                            flg=1;
                            break;
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[i][k], hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][k], hcsp4c[i][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]= hcsp4c[i][k];
                    trial[k+5]= hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

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

void Cluster_Write_BCC15() {
    int i;
    sBCC_15[hcBCC_15[nBCC_15][0]] = 'S';
    for (i = 1; i< 15; i++){
        if (sBCC_15[hcBCC_15[nBCC_15][i]] == 'C') sBCC_15[hcBCC_15[nBCC_15][i]] = 'B';
    }

    ++nBCC_15;
}