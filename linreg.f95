subroutine linreg(DX,DY,DA,DB,DR2,DN)
    implicit none
    integer*8 :: l,DN
    real*8 :: SX,SY,SXX,SXY,SYY,R,NN
    real*8,intent(out) :: DX(DN),DY(DN)
    real*8,intent(out) :: DB,DR2,DA
    

    SX  = 0.0
    SY  = 0.0
    SXX = 0.0
    SXY = 0.0
    SYY = 0.0

    do l=1,DN
        SX  = SX + DX(l);
        SY  = SY + DY(l);
        SXX = SXX + DX(l)*DX(l);
        SXY = SXY + DX(l)*DY(l);
        SYY = SYY + DY(l)*DY(l);
    end do
    NN  = (1.d0*DN)   
    DB   = (NN*SXY-SX*SY)/(NN*SXX-SX*SX);
    DA   = (SY-DB*SX)/NN;
    R   = (NN*SXY-SX*SY)/sqrt((NN*SXX-SX*SX)*(NN*SYY-SY*SY));
    DR2  = R*R;    
 
end subroutine linreg
