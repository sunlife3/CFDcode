program pltdist
        implicit none
        integer :: fo=20,i=0,j=0,GridNum=18,Grid=19,Qbin=100,index=1,Nx,Ny
        character filenameraw*128,filename*128
        real(8),allocatable :: X(:,:),Y(:,:),Q(:,:,:)
        real(8) :: rho,u,v,p,t,xcenter,ycenter,hoge,beta,theta,rhoinf,pinf,Minf,Vn,Vt
        !beta ... shock wave angle theta ... deflection angle
        real(8),parameter :: GAMMA=1.4d0,PI = 3.14159265359d0

        open(GridNum,file = 'GridNum.txt')
        read(GridNum,*)Nx,Ny
        close(GridNum)
        write(*,*)Nx,Ny

        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),Q(-1:Nx+2,-1:Ny+2,4))
        
        !open(Grid,file = 'Gridtube.txt')
        open(Grid,file = 'MESH_rampGrid.txt')
        do j=-3,Ny+4
            do i=-3,Nx+4
                read(Grid,*)X(i,j),Y(i,j)
            enddo
        enddo
        close(Grid)

        rhoinf = 1.0d0
        Minf   = 2.0d0
        pinf   = 1.0d0/GAMMA
        beta = 34.4d0*PI/180d0
        theta = 5.0d0*PI/180d0
        
        do index=1,1
            !write(filenameraw,'("QbinSD-SLAU",i3.3,".dat")')index
            write(filenameraw,'("QbinAUSMDV",i3.3,".dat")')index
            !write(filenameraw,'("QbinSLAU2",i3.3,".dat")')index
            !write(filenameraw,'("QbinSD-SLAU2",i3.3,".dat")')index
            open(Qbin,file = filenameraw)
            do j=-1,Ny+2
                do i=-1,Nx+2
                    read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                enddo
            enddo
            close(Qbin)
            write(filename,'("ramp",f3.1,"deg_M",f4.1,"_AUSMDV",i3.3,".dat")')theta*180d0/PI,Minf,index
            open(fo,file = filename)
            write(fo,*)'x',' ','AUSMDV',' ','AUSMDV',' ',' AUSMDV',' ','AUSMDV'
            !write(fo,*)'x',' ','SLAU2',' ','SLAU2',' ',' SLAU2',' ','SLAU2'
            !write(fo,*)'x',' ','SD-SLAU',' ','SD-SLAU',' ',' SD-SLAU',' ','SD-SLAU'
            !write(fo,*)'x',' ','SD-SLAU2',' ','SD-SLAU2',' ',' SD-SLAU2',' ','SD-SLAU2'
            do i=1,Nx
                j = 1
                rho = Q(i,j,1)
                u   = Q(i,j,2)/Q(i,j,1)
                v   = Q(i,j,3)/Q(i,j,1)
                p   = (GAMMA - 1.0d0)*(Q(i,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                Vn  = u*sin(beta) - v*cos(beta)
                Vt  = u*cos(beta) + v*sin(beta)
                write(fo,*)hoge/Nx-0.5d0,rho,Vn,Vt,p
                hoge = hoge + 1.0d0
            enddo
            close(fo)
        enddo
        
        !!output exact solution
        !fo = fo+1
        !hoge =0.0d0
        !write(filename,'("exactRamp",f3.1,"_M",f4.1,".dat")')theta*180d0/PI,Minf
        !open(fo,file = filename)
        !write(fo,*)'x',' ','rho_exact',' ','Vn_exact',' ',' Vt_exact',' ','p_exact'
        !do i=1,Nx
        !    if(i<=20)then
        !        write(fo,*)hoge/Nx-0.5d0,1.0d0,Minf*sin(beta),Minf*cos(beta),1.0d0/1.4d0
        !    else
        !        rho = (((GAMMA + 1.0d0)*(Minf*sin(beta))**2.0d0)/((GAMMA-1.0d0)*(Minf*sin(beta))**2.0d0 + 2.0d0)) * rhoinf
        !        p   = ((2.0d0*GAMMA*((Minf*sin(beta))**2.0d0))/(GAMMA + 1.0d0) - (GAMMA - 1.0d0)/(GAMMA + 1.0d0)) * pinf
        !        Vn  = ((GAMMA - 1.0d0)/(GAMMA + 1.0d0) + 2.0d0/((GAMMA + 1.0d0)*((Minf*sin(beta))**2.0d0))) * Minf*sin(beta)
        !        Vt  = Minf*cos(beta)
        !       write(fo,*)hoge/Nx-0.5d0,rho,Vn,Vt,p
        !    endif
        !    hoge = hoge + 1.0d0
        !enddo
        !close(fo)

end program pltdist