#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define MAX 10000

void linreg(double X[MAX],double Y[MAX],double A[1],double B[1],double R2[1],int N);

int main(){
    char dummy[2];              // Dummy character to store atoms names
    int i,j,k;                  // i, j and k counters        
    int n_particles;            // Number of particles in the system   
    int tmax;                   // Maximum time and origin time limit
    int omax;                   // Origin for multiple origins
    int omax_cut;               // Origin for multiple origins
    int t,o;                    // t and o counters  
    double dt;                  // dt between samples
    double B[1],R2[1],A[1];     // Variables for linear regression        
    double dif;                 // Self-diffusion coefficient (m²/s)
    double sig;                 // Stantard deviation
    double time[MAX];                // Time array
    double msd[MAX];                 // MSD array
    double msc[MAX];                 // msc array
    double msc_cut[MAX];             // msc_cut array
    double time_cut[MAX];            // time_cut array
    double tmpx,tmpy,tmpz;           // Auxiliar variables for position
    double **rx = malloc(5000 * sizeof(double *)); // Position arrays
    double **ry = malloc(5000 * sizeof(double *)); // Position arrays
    double **rz = malloc(5000 * sizeof(double *)); // Position arrays

    //!~~ Time related inputs ~~!
    tmax    = 5000;                     // Simulation time outputs
    dt      = 10.e0;                    // dt between outputs
    omax    = (int)((double)(tmax)/2.e0);   // Setting limit for time origins

    //  Memory allocation
    for(i=0;i<5000;i++){
        rx[i] = malloc(1000 * sizeof(double));
        ry[i] = malloc(1000 * sizeof(double));
        rz[i] = malloc(1000 * sizeof(double));
    }

    // Opening the .xyz file
    FILE *in;
    in = fopen("traj.xyz","r");

    // Reading positions of each time frame
    for(t=0;t<tmax;t++){
        printf("\rStep = %d",t);
        fflush(stdout);
        fscanf(in,"%d",&n_particles);
        fscanf(in,"%s",dummy);
        for(i=0;i<n_particles;i++){
            fscanf(in,"%s %lf %lf %lf",dummy,&tmpx,&tmpy,&tmpz);
            rx[t][i] = tmpx;
            ry[t][i] = tmpy;
            rz[t][i] = tmpz;
        }
    }

    printf("\n");
    fclose(in);

    // Opening the msd file
    FILE *out1,*out2;
    out1 = fopen("msd_c.dat","w");
    out2 = fopen("msd_cut_c.dat","w");

    // Emptying the arrays of msd and msc
    for(o=0;o<omax;o++){
        msd[o] = 0.e0;
        msc[o] = 0.e0;
    }
    
    printf("Calculation the MSD with multiple time origins...\n");
    //~~ Calculating the MSD correlation with multiple t=0
    for(o=0;o<omax;o++){
        if((o % 100) == 0){
            printf("\r %lf %% completo",((double)(o)/(double)(omax))*100.0);
            fflush(stdout);
        }
        time[o] = (double)(o)*(double)(dt);
        for(t=0;t<omax;t++){
            for(k=0;k<n_particles;k++){ 
                msd[o] = msd[o] + (rx[t+o][k]-rx[t][k])*(rx[t+o][k]-rx[t][k]) 
                                + (ry[t+o][k]-ry[t][k])*(ry[t+o][k]-ry[t][k]) 
                                + (rz[t+o][k]-rz[t][k])*(rz[t+o][k]-rz[t][k]);
            }
        }
        msc[o] = msd[o]/(double)(n_particles)/(double)(omax)/3.e0;
        fprintf(out1,"%lf %lf\n",time[o],msc[o]);
    }
    printf("\n");
    fclose(out1);

    for(t=(int)(0.2e0*omax);t<(int)(0.8e0*omax);t++){
        msc_cut[(int)(t-0.2e0*omax)]  = msc[t];
        time_cut[(int)(t-0.2e0*omax)] = time[t];
        fprintf(out2,"%lf %lf\n",time_cut[(int)(t-0.2e0*omax)],msc_cut[(int)(t-0.2e0*omax)]);
    }
    fclose(out2);

    linreg(time,msc,A,B,R2,omax);
    dif = B[0]/2.e0/1.0e5;
    sig = dif*sqrt(1.e0/R2[0]-1.e0)/sqrt(1.e0*(omax-1));
    printf("\nD = %le +/- %le m²/s\n",dif,sig);

    omax_cut = (int)((0.6e0*omax)-1);
    linreg(time_cut,msc_cut,A,B,R2,omax_cut);
    dif = B[0]/2.e0/1.0e5;
    sig = dif*sqrt(1.e0/R2[0]-1.e0)/sqrt(1.e0*(omax_cut-1));
    printf("\nD = %le +/- %le m²/s\n",dif,sig);

//    out2 = fopen("diff.dat","w");
//    fprintf(out2,"%le %le\n",dif,sig);
//    printf("\nD = %le +/- %le m²/s\n",dif,sig);
//    fclose(out2);
return 0;
}

void linreg(double X[MAX],double Y[MAX], double A[1],double B[1],double R2[1], int N){
    int K;
    double SX,SY,SXX,SXY,SYY;
    double NN,R;

    SX  = 0.0;
    SY  = 0.0;
    SXX = 0.0;
    SXY = 0.0;
    SYY = 0.0;
    for(K=0;K<N;K++){
        SX  = SX+X[K];
        SY  = SY+Y[K];
        SXX = SXX+X[K]*X[K];
        SXY = SXY+X[K]*Y[K];
        SYY = SYY+Y[K]*Y[K];
    }
    NN      = (1.e0*N);
    B[0]    = (NN*SXY-SX*SY)/(NN*SXX-SX*SX);
    A[0]    = (SY-B[0]*SX)/NN;
    R       = (NN*SXY-SX*SY)/sqrt((NN*SXX-SX*SX)*(NN*SYY-SY*SY));
    R2[0]   = R*R;
}





