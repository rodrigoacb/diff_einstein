!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Student: Rodrigo Amaral Coutinho Bartolomeu RA:209688
! Supervisor: Dr. Luís Fernando Mercier Franco
! Universidade Estadual de Campinas -UNICAMP 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM diff
    !~~ Variable declaration ~~!
    implicit none
    character(4)   :: dummy             ! Dummy character to store atoms names
    integer*8 :: i,j,k                  ! i, j and k counters        
    integer*8 :: n_particles            ! Number of particles in the system   
    integer*8 :: tmax                   ! Maximum time and origin time limit
    integer*8 :: omax                   ! Origin for multiple origins
    integer*8 :: omax_cut               ! Origin for multiple origins
    integer*8 :: t,o                    ! t and o counters  
    real*8 :: dt                        ! dt between samples
    real*8 :: B,R2,A                    ! Variables for linear regression        
    real*8 :: dif                       ! Self-diffusion coefficient (m²/s)
    real*8 :: sig                       ! Stantard deviation
    real*8,dimension(:),allocatable   :: time       ! Time array
    real*8,dimension(:),allocatable   :: msd        ! MSD array
    real*8,dimension(:),allocatable   :: msc        ! msc array
    real*8,dimension(:),allocatable   :: msc_cut    ! msc_cut array
    real*8,dimension(:),allocatable   :: time_cut   ! time_cut array
    real*8,dimension(:,:),allocatable :: rx,ry,rz   ! Position arrays
    
    !~~ Time related inputs ~~!
    tmax    = 5000  ! Simulation time outputs
    dt      = 10.d0 ! dt btween outputs
    omax    = int(dble(tmax)/2.d0) ! Setting limit for tiem origins

    !~~ Memory allocation ~~!
    allocate(rx(1000,tmax),ry(1000,tmax),rz(1000,tmax))
    allocate(time(omax),msd(omax),msc(omax))
    allocate(time_cut((int(omax*0.6d0))),msc_cut((int(omax*0.6d0))))
    
    !~~ Opening the .xyz file ~~!
    open(unit=101,status="old",file="traj.xyz")
    
    !~~ Reading positions of each time frame ~~!
    do t=1,tmax
        !print*,'Step = ',t
        write(*,'(A7,I6)')'Step = ',t
        read(101,*)n_particles
        !read(101,*)dummy
        do i=1,n_particles
            read(101,*)dummy,rx(i,t),ry(i,t),rz(i,t)
        end do    
    end do
    close(101)
    write(*,*)''

    !~~ Opening the msd file ~~!
    open(unit=102,file="msd_f.dat")
    open(unit=103,file="msd_cut_f.dat")

    !~~ Emptying the arrays of msd and msc~~!
    do o=1,omax
        msd(o) = 0.d0
        msc(o) = 0.d0
    end do
    
    print*,"Calculation the MSD with multiple time origins..."
    !~~ Calculating the MSD correlation with multiple t=0 ~~!
    do o=1,omax
        if(mod(o,100) == 0)then
            write(*,'(F6.2)')dble(dble(o)/dble(omax)*100.0)
        end if
        
        time(o) = dble(o)*dt
        do t=1,omax
            do k=1,n_particles
                msd(o) = msd(o) + (rx(k,t+o)-rx(k,t))**2 + (ry(k,t+o)-ry(k,t))**2 + (rz(k,t+o)-rz(k,t))**2
            end do    
        end do 
        msc(o) = msd(o)/dble(n_particles)/dble(omax)/3.d0
        write(102,*)time(o),msc(o)  
    end do
    close(102)
   
    do t=int(1+0.2d0*omax),int(0.8d0*omax)
        msc_cut(int(t-0.2d0*omax))  = msc(t)
        time_cut(int(t-0.2d0*omax)) = time(t)
        write(103,*)time_cut(int(t-0.2d0*omax)),msc_cut(int(t-0.2d0*omax))  
    end do
    close(103)
    !~~ Using linear regression to calculate the slope of the msc and correcting the unit to m²/s ~~!
    call linreg(time,msc,A,B,R2,omax)
    dif = B/2.d0/1.d5 !Conversion to m2/s ps -> 1D6
    sig = dif*sqrt(1.d0/R2-1.d0)/sqrt(1.d0*(omax-1))
    write(*,*) ' D =',dif,' +/-',sig,' m2/s'
    
    omax_cut = int(0.6d0*omax)-1
    call linreg(time_cut,msc_cut,A,B,R2,omax_cut)
    dif = B/2.d0/1.d5 !Conversion to m2/s
    sig = dif*sqrt(1.d0/R2-1.d0)/sqrt(1.d0*(omax_cut-1))
    write(*,*) ' D =',dif,' +/-',sig,' m2/s'

    !~~ Deallocating memory ~~!
    deallocate(rx,ry,rz)
    deallocate(time,msd,msc,time_cut,msc_cut)
    
    end program diff

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
