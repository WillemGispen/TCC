#include "stats.h"
#include "globals.h"
#include "tools.h"

void Stats_Init() {
    int i;
    char errMsg[1000];

    a5=malloc(N*sizeof(int));	if (a5==NULL) { sprintf(errMsg,"Stats_Init(): a5[] malloc out of memory\n");	Error(errMsg); }
    a6=malloc(N*sizeof(int));	if (a6==NULL) { sprintf(errMsg,"Stats_Init(): a6[] malloc out of memory\n");	Error(errMsg); }
    a7=malloc(N*sizeof(int));	if (a7==NULL) { sprintf(errMsg,"Stats_Init(): a7[] malloc out of memory\n");	Error(errMsg); }
    a8=malloc(N*sizeof(int));	if (a8==NULL) { sprintf(errMsg,"Stats_Init(): a8[] malloc out of memory\n");	Error(errMsg); }
    a9=malloc(N*sizeof(int));	if (a9==NULL) { sprintf(errMsg,"Stats_Init(): a9[] malloc out of memory\n");	Error(errMsg); }
    a10=malloc(N*sizeof(int));	if (a10==NULL) { sprintf(errMsg,"Stats_Init(): a10[] malloc out of memory\n");	Error(errMsg); }
    a11=malloc(N*sizeof(int));	if (a11==NULL) { sprintf(errMsg,"Stats_Init(): a11[] malloc out of memory\n");	Error(errMsg); }
    a12=malloc(N*sizeof(int));	if (a12==NULL) { sprintf(errMsg,"Stats_Init(): a12[] malloc out of memory\n");	Error(errMsg); }
    a13=malloc(N*sizeof(int));	if (a13==NULL) { sprintf(errMsg,"Stats_Init(): a13[] malloc out of memory\n");	Error(errMsg); }
    a15=malloc(N*sizeof(int));	if (a15==NULL) { sprintf(errMsg,"Stats_Init(): a15[] malloc out of memory\n");	Error(errMsg); }

    for(i=0;i<N;i++) {
        a5[i]=a6[i]=a7[i]=a8[i]=a9[i]=a10[i]=a11[i]=a12[i]=a13[i]=a15[i]=0;
    }

    ngsp3=ngsp3a=ngsp3b=ngsp3c=0;
    ngsp4=ngsp4a=ngsp4b=ngsp4c=0;
    ngsp5=ngsp5a=ngsp5b=ngsp5c=0;
    ng6Z=ng7K=0;
    ng8A=ng8B=ng8K=0;
    ng9A=ng9B=ng9K=0;
    ng10A=ng10B=ng10K=ng10W=0;
    ng11A=ng11B=ng11C=ng11E=ng11F=ng11W=0;	// number of particles in clusers gross
    ng12A=ng12B=ng12D=ng12E=ng12K=0;
    ng13A=ng13B=ng13K=0;
    ngFCC=ngHCP=ngBCC_9=ngBCC_15=0;

    nnsp3c=0;
    nnsp4c=0;
    nnsp5c=0;
    nn6Z=nn7K=0;
    nn8A=nn8B=nn8K=0;
    nn9A=nn9B=nn9K=0;
    nn10A=nn10B=nn10K=nn10W=0;
    nn11A=nn11B=nn11C=nn11E=nn11F=nn11W=0;	// number of particles in clusers gross
    nn12A=nn12B=nn12D=nn12E=nn12K=0;
    nn13A=nn13B=nn13K=0;
    nnFCC=nnHCP=nnBCC_9=nnBCC_15=0;

    ncsp3=ncsp3a=ncsp3b=ncsp3c=0;
    ncsp4=ncsp4a=ncsp4b=0;
    ncsp5=ncsp5a=ncsp5b=ncsp5c=0;
    nc6Z=nc7K=0;
    nc8A=nc8B=nc8K=0;
    nc9A=nc9B=nc9K=0;
    nc10A=nc10B=nc10K=nc10W=0;
    nc11A=nc11B=nc11C=nc11E=nc11F=nc11W=0;	// number of particles in clusers gross
    nc12A=nc12B=nc12D=nc12E=nc12K=0;
    nc13A=nc13B=nc13K=0;
    ncFCC=ncHCP=ncBCC_9=ncBCC_15=0;
}

void Stats_Reset() {
    int i;

    for(i=0;i<N;i++) {
        a5[i]=a6[i]=a7[i]=a8[i]=a9[i]=a10[i]=a11[i]=a12[i]=a13[i]=a15[i]=0;
    }
}

void Stats_FreeMem() {
    free(a5);
    free(a6);
    free(a7);
    free(a8);
    free(a9);
    free(a10);
    free(a11);
    free(a12);
    free(a13);
    free(a15);
}

void Stats_SetA() { // Set arrays to true if the ith particle is a member of any clusters with this or a larger number of particles
    int flg1, flg2;
    int i;

    for (i=0;i<N;i++) {
        flg1 = sFCC[i] != 'C' || sHCP[i] != 'C' || sBCC_15[i] != 'C' || s13A[i] != 'C' || s13B[i] != 'C' || s13K[i] != 'C';
        if(flg1==1) a13[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s12A[i] != 'C' || s12B[i] != 'C' || s12D[i] != 'C' || s12E[i] != 'C' || s12K[i] != 'C';
        flg2 = a13[i];
        if(flg1==1 || flg2==1) a12[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s11A[i] != 'C' || s11B[i] != 'C' || s11C[i] != 'C' || s11E[i] != 'C' || s11F[i] != 'C' || s11W[i] != 'C';
        flg2 = a12[i];
        if(flg1==1 || flg2==1) a11[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s10A[i] != 'C' || s10B[i] != 'C' || s10K[i] != 'C' || s10W[i] != 'C';
        flg2 = a11[i];
        if(flg1==1 || flg2==1) a10[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s9A[i] != 'C' || s9B[i] != 'C' || s9K[i] != 'C' || sBCC_9[i] != 'C';
        flg2 = a10[i];
        if(flg1==1 || flg2==1) a9[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s8A[i] != 'C' || s8B[i] != 'C' || s8K[i] != 'C';
        flg2 = a9[i];
        if(flg1==1 || flg2==1) a8[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = ssp5c[i] != 'C' || s7K[i] != 'C';
        flg2 = a8[i];
        if(flg1==1 || flg2==1) a7[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = ssp4c[i] != 'C' || s6Z[i] != 'C';
        flg2 = a7[i];
        if(flg1==1 || flg2==1) a6[i] = 1;
    }
}

void Stats_Analyse() {
    int i;

    for(i=0; i<N; ++i){
        if(ssp3[i] != 'C') ++ngsp3;
        if(ssp3a[i] != 'C') ++ngsp3a;
        if(ssp3b[i] != 'C') ++ngsp3b;
        if(ssp3c[i] != 'C') ++ngsp3c;
        if(ssp4[i] != 'C') ++ngsp4;
        if(ssp4a[i] != 'C') ++ngsp4a;
        if(ssp4b[i] != 'C') ++ngsp4b;
        if(ssp4c[i] != 'C') ++ngsp4c;
        if(s6Z[i] != 'C') ++ng6Z;
        if(s7K[i] != 'C') ++ng7K;
        if(ssp5[i] != 'C') ++ngsp5;
        if(ssp5a[i] != 'C') ++ngsp5a;
        if(ssp5b[i] != 'C') ++ngsp5b;
        if(ssp5c[i] != 'C') ++ngsp5c;
        if(s8A[i] != 'C') ++ng8A;
        if(s8B[i] != 'C') ++ng8B;
        if(s8K[i] != 'C') ++ng8K;
        if(s9A[i] != 'C') ++ng9A;
        if(s9B[i] != 'C') ++ng9B;
        if(s9K[i] != 'C') ++ng9K;
        if(s10A[i] != 'C') ++ng10A;
        if(s10B[i] != 'C') ++ng10B;
        if(s10K[i] != 'C') ++ng10K;
        if(s10W[i] != 'C') ++ng10W;
        if(s11A[i] != 'C') ++ng11A;
        if(s11B[i] != 'C') ++ng11B;
        if(s11C[i] != 'C') ++ng11C;
        if(s11E[i] != 'C') ++ng11E;
        if(s11F[i] != 'C') ++ng11F;
        if(s11W[i] != 'C') ++ng11W;
        if(s12A[i] != 'C') ++ng12A;
        if(s12B[i] != 'C') ++ng12B;
        if(s12D[i] != 'C') ++ng12D;
        if(s12E[i] != 'C') ++ng12E;
        if(s12K[i] != 'C') ++ng12K;
        if(s13A[i] != 'C') ++ng13A;
        if(s13B[i] != 'C') ++ng13B;
        if(s13K[i] != 'C') ++ng13K;
        if(sFCC[i] != 'C') ++ngFCC;
        if(sHCP[i] != 'C') ++ngHCP;
        if(sBCC_9[i] != 'C') ++ngBCC_9;
        if(sBCC_15[i] != 'C') ++ngBCC_15;
    }
    Stats_SetA();
    for(i=0; i<N; ++i){
        if(ssp3c[i] != 'C' && !a6[i]) ++nnsp3c;
        if(ssp4c[i] != 'C' && !a7[i]) ++nnsp4c;
        if(s6Z[i] != 'C' && !a7[i]) ++nn6Z;
        if(s7K[i] != 'C' && !a7[i]) ++nn7K;
        if(ssp5c[i] != 'C' && !a8[i]) ++nnsp5c;
        if(s8A[i] != 'C' && !a9[i]) ++nn8A;
        if(s8B[i] != 'C' && !a9[i]) ++nn8B;
        if(s8K[i] != 'C' && !a9[i]) ++nn8K;
        if(s9A[i] != 'C' && !a10[i]) ++nn9A;
        if(s9B[i] != 'C' && !a10[i]) ++nn9B;
        if(s9K[i] != 'C' && !a10[i]) ++nn9K;
        if(s10A[i] != 'C' && !a11[i]) ++nn10A;
        if(s10B[i] != 'C' && !a11[i]) ++nn10B;
        if(s10K[i] != 'C' && !a11[i]) ++nn10K;
        if(s10W[i] != 'C' && !a11[i]) ++nn10W;
        if(s11A[i] != 'C' && !a12[i]) ++nn11A;
        if(s11B[i] != 'C' && !a12[i]) ++nn11B;
        if(s11C[i] != 'C' && !a12[i]) ++nn11C;
        if(s11E[i] != 'C' && !a12[i]) ++nn11E;
        if(s11F[i] != 'C' && !a12[i]) ++nn11F;
        if(s11W[i] != 'C' && !a12[i]) ++nn11W;
        if(s12A[i] != 'C' && !a13[i]) ++nn12A;
        if(s12B[i] != 'C' && !a13[i]) ++nn12B;
        if(s12D[i] != 'C' && !a13[i]) ++nn12D;
        if(s12E[i] != 'C' && !a13[i]) ++nn12E;
        if(s12K[i] != 'C' && !a13[i]) ++nn12K;
        if(s13A[i] != 'C') ++nn13A;
        if(s13B[i] != 'C') ++nn13B;
        if(s13K[i] != 'C') ++nn13K;
        if(sFCC[i] != 'C') ++nnFCC;
        if(sHCP[i] != 'C') ++nnHCP;
        if(sBCC_9[i] != 'C') ++nnBCC_9;
        if(sBCC_15[i] != 'C') ++nnBCC_15;
    }
}

void Accuumlate_Stats() {

    ncsp3+=nsp3;
    ncsp3a+=nsp3a;
    ncsp3b+=nsp3b;
    ncsp3c+=nsp3c;
    ncsp4+=nsp4;
    ncsp4a+=nsp4a;
    ncsp4b+=nsp4b;
    ncsp4c+=nsp4c;
    nc6Z+=n6Z;
    nc7K+=n7K;
    ncsp5+=nsp5;
    ncsp5a+=nsp5a;
    ncsp5b+=nsp5b;
    ncsp5c+=nsp5c;
    nc8A+=n8A;
    nc8B+=n8B;
    nc8K+=n8K;
    nc9A+=n9A;
    nc9B+=n9B;
    nc9K+=n9K;
    nc10A+=n10A;
    nc10B+=n10B;
    nc10K+=n10K;
    nc10W+=n10W;
    nc11A+=n11A;
    nc11B+=n11B;
    nc11C+=n11C;
    nc11E+=n11E;
    nc11F+=n11F;
    nc11W+=n11W;
    nc12A+=n12A;
    nc12B+=n12B;
    nc12D+=n12D;
    nc12E+=n12E;
    nc12K+=n12K;
    nc13A+=n13A;
    nc13B+=n13B;
    nc13K+=n13K;
    ncFCC+=nFCC;
    ncHCP+=nHCP;
    ncBCC_9+=nBCC_9;
    ncBCC_15+=nBCC_15;
}

void Stats_Report(char *filename) {
    char errMsg[1000];
    FILE *writeout;

    printf("Clust	Number	Gross	Net	Mean Pop Per Frame\n");
    printf("sp3	%d	%d	%d	%.5lg\n",ncsp3,ngsp3,0,mean_pop_per_frame[0]);
    printf("sp3a	%d	%d	%d	%.5lg\n",ncsp3a,ngsp3a,0,mean_pop_per_frame[1]);
    printf("sp3b	%d	%d	%d	%.5lg\n",ncsp3b,ngsp3b,0,mean_pop_per_frame[2]);
    printf("sp3c	%d	%d	%d	%.5lg\n",ncsp3c,ngsp3c,nnsp3c,mean_pop_per_frame[3]);
    printf("sp4	%d	%d	%d	%.5lg\n",ncsp4,ngsp4,0,mean_pop_per_frame[4]);
    printf("sp4a	%d	%d	%d	%.5lg\n",ncsp4a,ngsp4a,0,mean_pop_per_frame[5]);
    printf("sp4b	%d	%d	%d	%.5lg\n",ncsp4b,ngsp4b,0,mean_pop_per_frame[6]);
    printf("sp4c	%d	%d	%d	%.5lg\n",ncsp4c,ngsp4c,0,mean_pop_per_frame[7]);
    printf("sp5	%d	%d	%d	%.5lg\n",ncsp5,ngsp5,0,mean_pop_per_frame[8]);
    printf("sp5a	%d	%d	%d	%.5lg\n",ncsp5a,ngsp5a,0,mean_pop_per_frame[9]);
    printf("sp5b	%d	%d	%d	%.5lg\n",ncsp5b,ngsp5b,0,mean_pop_per_frame[10]);
    printf("sp5c	%d	%d	%d	%.5lg\n",ncsp5c,ngsp5c,nnsp5c,mean_pop_per_frame[11]);
    printf("6Z_C2v	%d	%d	%d	%.5lg\n",nc6Z,ng6Z,nn6Z,mean_pop_per_frame[12]);
    printf("7K	%d	%d	%d	%.5lg\n",nc7K,ng7K,nn7K,mean_pop_per_frame[13]);
    printf("8A_D2d	%d	%d	%d	%.5lg\n",nc8A,ng8A,nn8A,mean_pop_per_frame[14]);
    printf("8B_Cs	%d	%d	%d	%.5lg\n",nc8B,ng8B,nn8B,mean_pop_per_frame[15]);
    printf("8K	%d	%d	%d	%.5lg\n",nc8K,ng8K,nn8K,mean_pop_per_frame[16]);
    printf("9A_D3h	%d	%d	%d	%.5lg\n",nc9A,ng9A,nn9A,mean_pop_per_frame[17]);
    printf("9B_C2v	%d	%d	%d	%.5lg\n",nc9B,ng9B,nn9B,mean_pop_per_frame[18]);
    printf("9K	%d	%d	%d	%.5lg\n",nc9K,ng9K,nn9K,mean_pop_per_frame[19]);
    printf("10A_D4d	%d	%d	%d	%.5lg\n",nc10A,ng10A,nn10A,mean_pop_per_frame[20]);
    printf("10B_C3v	%d	%d	%d	%.5lg\n",nc10B,ng10B,nn10B,mean_pop_per_frame[21]);
    printf("10K	%d	%d	%d	%.5lg\n",nc10K,ng10K,nn10K,mean_pop_per_frame[22]);
    printf("10W	%d	%d	%d	%.5lg\n",nc10W,ng10W,nn10W,mean_pop_per_frame[23]);
    printf("11A_D4d	%d	%d	%d	%.5lg\n",nc11A,ng11A,nn11A,mean_pop_per_frame[24]);
    printf("11B_C2v	%d	%d	%d	%.5lg\n",nc11B,ng11B,nn11B,mean_pop_per_frame[25]);
    printf("11CD	%d	%d	%d	%.5lg\n",nc11C,ng11C,nn11C,mean_pop_per_frame[26]);
    printf("11E_C2	%d	%d	%d	%.5lg\n",nc11E,ng11E,nn11E,mean_pop_per_frame[27]);
    printf("11F_C2v	%d	%d	%d	%.5lg\n",nc11F,ng11F,nn11F,mean_pop_per_frame[28]);
    printf("11W_Cs	%d	%d	%d	%.5lg\n",nc11W,ng11W,nn11W,mean_pop_per_frame[29]);
    printf("12A_C2v	%d	%d	%d	%.5lg\n",nc12A,ng12A,nn12A,mean_pop_per_frame[30]);
    printf("12BC	%d	%d	%d	%.5lg\n",nc12B,ng12B,nn12B,mean_pop_per_frame[31]);
    printf("12D_D2d	%d	%d	%d	%.5lg\n",nc12D,ng12D,nn12D,mean_pop_per_frame[32]);
    printf("12E_D3h	%d	%d	%d	%.5lg\n",nc12E,ng12E,nn12E,mean_pop_per_frame[33]);
    printf("12K	%d	%d	%d	%.5lg\n",nc12K,ng12K,nn12K,mean_pop_per_frame[34]);
    printf("13A_Ih	%d	%d	%d	%.5lg\n",nc13A,ng13A,nn13A,mean_pop_per_frame[35]);
    printf("13B_D5h	%d	%d	%d	%.5lg\n",nc13B,ng13B,nn13B,mean_pop_per_frame[36]);
    printf("13K	%d	%d	%d	%.5lg\n",nc13K,ng13K,nn13K,mean_pop_per_frame[37]);
    printf("FCC_m13	%d	%d	%d	%.5lg\n",ncFCC,ngFCC,nnFCC,mean_pop_per_frame[38]);
    printf("HCP_m13	%d	%d	%d	%.5lg\n",ncHCP,ngHCP,nnHCP,mean_pop_per_frame[39]);
    printf("BCC_m9	%d	%d	%d	%.5lg\n",ncBCC_9,ngBCC_9,nnBCC_9,mean_pop_per_frame[40]);
    printf("BCC_m15	%d	%d	%d	%.5lg\n",ncBCC_15,ngBCC_15,nnBCC_15,mean_pop_per_frame[41]);
    printf("maxnB	%d\n",maxnb);

    printf("correctedBonds %d per frame %.5lg per part per frame %.5lg\n",correctedBonds,(double)correctedBonds/FRAMES,(double)correctedBonds/(FRAMES*N));


    writeout=fopen(filename,"w");
    if (writeout==NULL)  {
        sprintf(errMsg,"Stats_Report(): Error opening file %s",filename);	// Always test file open
        Error(errMsg);
    }
    fprintf(writeout,"%s\n",filename);
    fprintf(writeout,"Clust	Number	Gross	Net	Mean Pop Per Frame\n");
    fprintf(writeout,"sp3	%d	%d	%d	%.5lg\n",ncsp3,ngsp3,0,mean_pop_per_frame[0]);
    fprintf(writeout,"sp3a	%d	%d	%d	%.5lg\n",ncsp3a,ngsp3a,0,mean_pop_per_frame[1]);
    fprintf(writeout,"sp3b	%d	%d	%d	%.5lg\n",ncsp3b,ngsp3b,0,mean_pop_per_frame[2]);
    fprintf(writeout,"sp3c	%d	%d	%d	%.5lg\n",ncsp3c,ngsp3c,nnsp3c,mean_pop_per_frame[3]);
    fprintf(writeout,"sp4	%d	%d	%d	%.5lg\n",ncsp4,ngsp4,0,mean_pop_per_frame[4]);
    fprintf(writeout,"sp4a	%d	%d	%d	%.5lg\n",ncsp4a,ngsp4a,0,mean_pop_per_frame[5]);
    fprintf(writeout,"sp4b	%d	%d	%d	%.5lg\n",ncsp4b,ngsp4b,0,mean_pop_per_frame[6]);
    fprintf(writeout,"sp4c	%d	%d	%d	%.5lg\n",ncsp4c,ngsp4c,0,mean_pop_per_frame[7]);
    fprintf(writeout,"sp5	%d	%d	%d	%.5lg\n",ncsp5,ngsp5,0,mean_pop_per_frame[8]);
    fprintf(writeout,"sp5a	%d	%d	%d	%.5lg\n",ncsp5a,ngsp5a,0,mean_pop_per_frame[9]);
    fprintf(writeout,"sp5b	%d	%d	%d	%.5lg\n",ncsp5b,ngsp5b,0,mean_pop_per_frame[10]);
    fprintf(writeout,"sp5c	%d	%d	%d	%.5lg\n",ncsp5c,ngsp5c,nnsp5c,mean_pop_per_frame[11]);
    fprintf(writeout,"6Z_C2v	%d	%d	%d	%.5lg\n",nc6Z,ng6Z,nn6Z,mean_pop_per_frame[12]);
    fprintf(writeout,"7K	%d	%d	%d	%.5lg\n",nc7K,ng7K,nn7K,mean_pop_per_frame[13]);
    fprintf(writeout,"8A_D2d	%d	%d	%d	%.5lg\n",nc8A,ng8A,nn8A,mean_pop_per_frame[14]);
    fprintf(writeout,"8B_Cs	%d	%d	%d	%.5lg\n",nc8B,ng8B,nn8B,mean_pop_per_frame[15]);
    fprintf(writeout,"8K	%d	%d	%d	%.5lg\n",nc8K,ng8K,nn8K,mean_pop_per_frame[16]);
    fprintf(writeout,"9A_D3h	%d	%d	%d	%.5lg\n",nc9A,ng9A,nn9A,mean_pop_per_frame[17]);
    fprintf(writeout,"9B_C2v	%d	%d	%d	%.5lg\n",nc9B,ng9B,nn9B,mean_pop_per_frame[18]);
    fprintf(writeout,"9K	%d	%d	%d	%.5lg\n",nc9K,ng9K,nn9K,mean_pop_per_frame[19]);
    fprintf(writeout,"10A_D4d	%d	%d	%d	%.5lg\n",nc10A,ng10A,nn10A,mean_pop_per_frame[20]);
    fprintf(writeout,"10B_C3v	%d	%d	%d	%.5lg\n",nc10B,ng10B,nn10B,mean_pop_per_frame[21]);
    fprintf(writeout,"10K	%d	%d	%d	%.5lg\n",nc10K,ng10K,nn10K,mean_pop_per_frame[22]);
    fprintf(writeout,"10W	%d	%d	%d	%.5lg\n",nc10W,ng10W,nn10W,mean_pop_per_frame[23]);
    fprintf(writeout,"11A_D4d	%d	%d	%d	%.5lg\n",nc11A,ng11A,nn11A,mean_pop_per_frame[24]);
    fprintf(writeout,"11B_C2v	%d	%d	%d	%.5lg\n",nc11B,ng11B,nn11B,mean_pop_per_frame[25]);
    fprintf(writeout,"11CD	%d	%d	%d	%.5lg\n",nc11C,ng11C,nn11C,mean_pop_per_frame[26]);
    fprintf(writeout,"11E_C2	%d	%d	%d	%.5lg\n",nc11E,ng11E,nn11E,mean_pop_per_frame[27]);
    fprintf(writeout,"11F_C2v	%d	%d	%d	%.5lg\n",nc11F,ng11F,nn11F,mean_pop_per_frame[28]);
    fprintf(writeout,"11W_Cs	%d	%d	%d	%.5lg\n",nc11W,ng11W,nn11W,mean_pop_per_frame[29]);
    fprintf(writeout,"12A_C2v	%d	%d	%d	%.5lg\n",nc12A,ng12A,nn12A,mean_pop_per_frame[30]);
    fprintf(writeout,"12BC	%d	%d	%d	%.5lg\n",nc12B,ng12B,nn12B,mean_pop_per_frame[31]);
    fprintf(writeout,"12D_D2d	%d	%d	%d	%.5lg\n",nc12D,ng12D,nn12D,mean_pop_per_frame[32]);
    fprintf(writeout,"12E_D3h	%d	%d	%d	%.5lg\n",nc12E,ng12E,nn12E,mean_pop_per_frame[33]);
    fprintf(writeout,"12K	%d	%d	%d	%.5lg\n",nc12K,ng12K,nn12K,mean_pop_per_frame[34]);
    fprintf(writeout,"13A_Ih	%d	%d	%d	%.5lg\n",nc13A,ng13A,nn13A,mean_pop_per_frame[35]);
    fprintf(writeout,"13B_D5h	%d	%d	%d	%.5lg\n",nc13B,ng13B,nn13B,mean_pop_per_frame[36]);
    fprintf(writeout,"13K	%d	%d	%d	%.5lg\n",nc13K,ng13K,nn13K,mean_pop_per_frame[37]);
    fprintf(writeout,"FCC_m13	%d	%d	%d	%.5lg\n",ncFCC,ngFCC,nnFCC,mean_pop_per_frame[38]);
    fprintf(writeout,"HCP_m13	%d	%d	%d	%.5lg\n",ncHCP,ngHCP,nnHCP,mean_pop_per_frame[39]);
    fprintf(writeout,"BCC_m9	%d	%d	%d	%.5lg\n",ncBCC_9,ngBCC_9,nnBCC_9,mean_pop_per_frame[40]);
    fprintf(writeout,"BCC_m15	%d	%d	%d	%.5lg\n",ncBCC_15,ngBCC_15,nnBCC_15,mean_pop_per_frame[41]);
    fprintf(writeout,"maxnB	%d\n",maxnb);

    fprintf(writeout,"correctedBonds	%d	per frame	%.15lg	per part per frame	%.15lg\n",correctedBonds,(double)correctedBonds/FRAMES,(double)correctedBonds/(FRAMES*N));

    fclose(writeout);
}

void Pop_Per_Frame(int f) {
    int i;

    for(i=0; i<N; ++i){
        if(ssp3[i] != 'C') pop_per_frame[0][f]+=1.0;
        if(ssp3a[i] != 'C') pop_per_frame[1][f]+=1.0;
        if(ssp3b[i] != 'C') pop_per_frame[2][f]+=1.0;
        if(ssp3c[i] != 'C') pop_per_frame[3][f]+=1.0;
        if(ssp4[i] != 'C') pop_per_frame[4][f]+=1.0;
        if(ssp4a[i] != 'C') pop_per_frame[5][f]+=1.0;
        if(ssp4b[i] != 'C') pop_per_frame[6][f]+=1.0;
        if(ssp4c[i] != 'C') pop_per_frame[7][f]+=1.0;
        if(ssp5[i] != 'C') pop_per_frame[8][f]+=1.0;
        if(ssp5a[i] != 'C') pop_per_frame[9][f]+=1.0;
        if(ssp5b[i] != 'C') pop_per_frame[10][f]+=1.0;
        if(ssp5c[i] != 'C') pop_per_frame[11][f]+=1.0;
        if(s6Z[i] != 'C') pop_per_frame[12][f]+=1.0;
        if(s7K[i] != 'C') pop_per_frame[13][f]+=1.0;
        if(s8A[i] != 'C') pop_per_frame[14][f]+=1.0;
        if(s8B[i] != 'C') pop_per_frame[15][f]+=1.0;
        if(s8K[i] != 'C') pop_per_frame[16][f]+=1.0;
        if(s9A[i] != 'C') pop_per_frame[17][f]+=1.0;
        if(s9B[i] != 'C') pop_per_frame[18][f]+=1.0;
        if(s9K[i] != 'C') pop_per_frame[19][f]+=1.0;
        if(s10A[i] != 'C') pop_per_frame[20][f]+=1.0;
        if(s10B[i] != 'C') pop_per_frame[21][f]+=1.0;
        if(s10K[i] != 'C') pop_per_frame[22][f]+=1.0;
        if(s10W[i] != 'C') pop_per_frame[23][f]+=1.0;
        if(s11A[i] != 'C') pop_per_frame[24][f]+=1.0;
        if(s11B[i] != 'C') pop_per_frame[25][f]+=1.0;
        if(s11C[i] != 'C') pop_per_frame[26][f]+=1.0;
        if(s11E[i] != 'C') pop_per_frame[27][f]+=1.0;
        if(s11F[i] != 'C') pop_per_frame[28][f]+=1.0;
        if(s11W[i] != 'C') pop_per_frame[29][f]+=1.0;
        if(s12A[i] != 'C') pop_per_frame[30][f]+=1.0;
        if(s12B[i] != 'C') pop_per_frame[31][f]+=1.0;
        if(s12D[i] != 'C') pop_per_frame[32][f]+=1.0;
        if(s12E[i] != 'C') pop_per_frame[33][f]+=1.0;
        if(s12K[i] != 'C') pop_per_frame[34][f]+=1.0;
        if(s13A[i] != 'C') pop_per_frame[35][f]+=1.0;
        if(s13B[i] != 'C') pop_per_frame[36][f]+=1.0;
        if(s13K[i] != 'C') pop_per_frame[37][f]+=1.0;
        if(sFCC[i] != 'C') pop_per_frame[38][f]+=1.0;
        if(sHCP[i] != 'C') pop_per_frame[39][f]+=1.0;
        if(sBCC_9[i] != 'C') pop_per_frame[40][f]+=1.0;
        if(sBCC_15[i] != 'C') pop_per_frame[41][f]+=1.0;
    }

    // Add pop per frame to running total over all frames
    for(i=0; i<num_cluster_types; i++) {
        pop_per_frame[i][f]/=N;
        mean_pop_per_frame[i] += pop_per_frame[i][f];
    }
}

void Normalise_Populations() {
    int i;
    for(i=0; i<num_cluster_types; i++) {
        mean_pop_per_frame[i]/=FRAMES;
    }
}