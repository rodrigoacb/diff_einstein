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
