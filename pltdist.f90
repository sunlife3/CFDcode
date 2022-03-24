program pltdist
        implicit none
        integer :: fo=20,i=0,j=0,GridNum=18,Grid=19,Qbin=100,index=1,Nx,Ny,is_cyl=0
        character filenameraw*128,filename*128,meshfile*64
        real(8),allocatable :: X(:,:),Y(:,:),Q(:,:,:)
        real(8) :: rho,u,v,p,t,c,xcenter,ycenter,beta,theta,rhoinf,uinf,pinf,Minf,Vn,Vt,dx,hoge,pp,pm,&
                    rho1,p1,t1,rho2,p2,t2,Sx,Sy,f
        !real(8) Vn2,T0,T1
        !beta ... shock wave angle theta ... deflection angle
        real(8),parameter :: GAMMA=1.4d0,PI = 3.14159265359d0

        open(GridNum,file = 'GridNum.txt')
        read(GridNum,*)Nx,Ny
        close(GridNum)
        write(*,*)Nx,Ny


        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),Q(-2:Nx+3,-2:Ny+3,4))
        
        !open(Grid,file = 'Gridtube.txt')
        !open(Grid,file = 'MESH_cylinder(100100).txt')
        !do j=-3,Ny+4
        !    do i=-3,Nx+4
        !        read(Grid,*)X(i,j),Y(i,j)
        !    enddo
        !enddo
        !close(Grid)

        dx = 1.0d0/Nx
        Minf  = 5.0d0
        rhoinf = 1.0d0
        uinf   = 1.0d0
        pinf   = 1.0d0/GAMMA
        !beta = 7.15d0*PI/180d0!M20
        beta = 7.9d0*PI/180d0!M15
        !beta = 8.1d0*PI/180d0!M8
        !beta = 9.5d0*PI/180d0!M10
        !beta = 10.9d0*PI/180d0!M8
        !beta = 11.8d0*PI/180d0!M7
        !beta = 15.1d0*PI/180d0!M5
        !beta = 18.0d0*PI/180d0!M4
        !beta = 23.1d0*PI/180d0!M3
        !beta = 34.3d0*PI/180d0!M2
        !beta = 90.0d0*PI/180.0d0! Normal Shock
        theta = 5.0d0*PI/180d0
        
        do index=1,1
            !write(filename,'("Qbin HLLAUSM Cylinder(M20,100)_",i3.3,".dat")')index
            !write(filename,'("QbinSD-AUSMDV SNS(M20,SFver2)_",i3.3,".dat")')index
            !write(filenameraw,'("Qbin_SD&2_ramp5deg(Case3)_",i1.1,".dat")')index
            write(filename,'("Qbin_SLAU2_ramp5deg(M15,100,1vs3)_",i1.1,".dat")')index
            !write(filename,'("Qbin_SWBLI(HLLAUSM,WCNS,2000,Re1000,M15)_",i2.2,".dat")')index
            !write(filename,'("Qbin_AUSM_duct5deg(M15,100)_",i1.1,".dat")')index
            !write(filename,'("Qbin_HLL_1stratio_50(HLL)&25(HLL)&0(HLL)_",i1.1,".dat")')index
            !write(filename,'("Qbin_HLLAUSM_SOD_MUSCL(200)_",i1.1,".dat")')index
            !write(filename,'("Qbin_SD_duct5deg(M15,100)_",i1.1,".dat")')index
            !write(filename,'("Qbin_HLL(WCNS)_ramp5deg(M15,100)_",i1.1,".dat")')index
            write(*,*)filename
            open(Qbin,file = filename)
            read(Qbin,*)meshfile
            write(*,*)meshfile
            open(Grid,file = meshfile)
            read(Grid,*)Nx,Ny
            write(*,*)Nx,Ny

            do j=-3,Ny+4
                do i=-3,Nx+4
                    read(Grid,*)X(i,j),Y(i,j)
                enddo
            enddo
            close(Grid)

            open(Qbin,file = filename)
            do j=-2,Ny+3
                do i=-2,Nx+3
                    read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                    !if(2.9d0 < Q(i,j,1))then
                    !    write(*,*)j
                    !endif
                enddo
            enddo
            write(*,*)'done'
            close(Qbin)
            !write(filename,'("ramp5deg_HLL(MUSCL)_(M15)(100)",i1.1,".dat")')index
            !write(filename,'("ramp5deg_HLLAUSMmix_(M15,Case3-2,triple)(100)",i1.1,".dat")')index
            write(filename,'("ramp5deg_SLAU2(M15,100,1vs3)",i1.1,".dat")')index
            !write(filename,'("HLLAUSM_Shu-Osher3rd_MUSCL(500)_",i1.1,".dat")')index
            !write(filename,'("duct5deg_SD-AUSMDV_(M15)",i1.1,".dat")')index
            !write(filename,'("SNS_SD-AUSMDV_(M20,SFver2)",i1.1,".dat")')index
            !write(filename,'("ramp5 grad Vn SD-AUSMDV_(M15,C,100)",i1.1,".dat")')index
            !write(filename,'("Cylinder_HLLAUSM(M20,100)_",i1.1,".dat")')index
            !write(filename,'("SOD_HLLAUSM_(200)",i1.1,".dat")')index
            !write(filename,'("SWBLI_HLLAUSM(Re1000,2000,WCNS,M15)_bottom",i2.2,".dat")')index
            open(fo,file = filename)
            !write(fo,*)'x,',' ','HLLAUSM,',' ','HLLAUSM,',' ','HLLAUSM,',' ','HLLAUSM'
            !write(fo,*)'x,',' ','rho,',' ','Vn,',' ','Vt,',' ','p'
            !write(fo,*)'x',' ','Vn',' ','Vn_i-Vn_(i-1)',' ','threshold'
            !write(fo,*)'x,',' ','SD-SLAU,',' ','SD-SLAU,',' ',' SD-SLAU,',' ','SD-SLAU'
            !write(fo,*)'x,',' ','SLAU2,',' ','SLAU2,',' ','SLAU2,',' ','SLAU2'
            !write(fo,*)'x,',' ','HLL,',' ','HLL,',' ','HLL,',' ','HLL'
            !write(fo,*)'x,',' ',&
            !write(fo,*)'x,',' ','3000,',' ','3000,',' ',' 3000,',' ','3000'
            !'HLL[1st:3rd=5:5](Shock)HLL[1st:3rd=1:3](Around)AUSMDV3rd(Other),',' ',&
            !'HLL[1st:3rd=5:5](Shock)HLL[1st:3rd=1:3](Around)AUSMDV3rd(Other),',' ',&
            !'HLL[1st:3rd=5:5](Shock)HLL[1st:3rd=1:3](Around)AUSMDV3rd(Other),',' ',&
            !'HLL[1st:3rd=5:5](Shock)HLL[1st:3rd=1:3](Around)AUSMDV3rd(Other)'
            !write(fo,*)'x,','M5'
            write(fo,*)'x,','rho,','Vn,','Vt,','p,'
            !write(fo,*)'x,','rho,','u,','v,','p'
            j = 1

            do i=1,Nx    
                rho = Q(i,j,1)
                u   = Q(i,j,2)/Q(i,j,1)
                v   = Q(i,j,3)/Q(i,j,1)
                p   = (GAMMA - 1.0d0)*(Q(i,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                c   = sqrt(GAMMA*p/rho)
                pm  = (GAMMA - 1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                Vn  = u*sin(beta) - v*cos(beta)
                Vt  = u*cos(beta) + v*sin(beta)
                pp  = (p - pinf)/pinf * ((GAMMA + 1.0d0)/(2.0d0*GAMMA)) * 1.0d0/((Minf*(sin(beta)))**2.0d0 -1.0d0)

                rho1 = Q(i,j,1)
                p1   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
                t1   = p1/rho1

                rho2 = Q(i+1,j,1)
                p2   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                t2   = p2/rho2
                Sx = p/(rho**GAMMA) !log((t2/t1)**(GAMMA/(GAMMA-1.0d0))*((p1**rho1)/(p2**rho2)))
                !write(*,*)i,log((p2/p1)*(rho1/rho2)**1.40d0)!
                f = min(p1/p2,p2/p1)**3.0d0

                rho2 = Q(i,j+1,1)
                p2   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
                t2   = p2/rho2
                Sy = p1/(rho1**GAMMA)!p2/(rho2**GAMMA) - p1/(rho1**GAMMA)!log((t1/t2)**(GAMMA/(GAMMA-1.0d0))*((p2**rho2)/(p1**rho1)))

                !write(fo,*)X(i,j)+0.1d0,',',pp
                write(fo,*)X(i,j),',',rho,',',Vn,',',Vt,',',p
                !write(fo,*)X(i,j)-0.5d0,',',rho,',',u,',',v,',',p
                
                !rho = Q(47,1,1)
                !u = Q(47,1,2)/Q(47,1,1)
                !v = Q(47,1,3)/Q(47,1,1)
                !p   = (GAMMA - 1.0d0)*(Q(47,1,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                !Vn  = u*sin(1.0d0*i*PI/180) - v*cos(1.0d0*i*PI/180)
                !Vt  = u*cos(1.0d0*i*PI/180) + v*sin(1.0d0*i*PI/180)
                !write(fo,*)i,',',u,',',v,',',Vn,',',Vt,',',sqrt(GAMMA*p/rho)
                
            enddo

            close(fo)

            !if(is_cyl == 1)then
            !    write(filename,'("Cyl gradT_SD-SLAU(M20,100)_",i1.1,".dat")')index
            !    open(fo,file = filename)

            !    write(fo,*)'deg',' ','T1-T0'

            !    do j=1,Ny

            !        rho = Q(Nx,j,1)
            !        u   = Q(Nx,j,2)/Q(Nx,j,1)
            !        v   = Q(Nx,j,3)/Q(Nx,j,1)
            !        T0   = (GAMMA - 1.0d0)*(Q(Nx,j,4)/Q(Nx,j,1) - 0.5d0*(u**2.0d0 + v**2.0d0))

            !        rho = Q(Nx-1,j,1)
            !        u   = Q(Nx-1,j,2)/Q(Nx-1,j,1)
            !        v   = Q(Nx-1,j,3)/Q(Nx-1,j,1)
            !        T1   = (GAMMA - 1.0d0)*(Q(Nx-1,j,4)/Q(Nx-1,j,1) - 0.5d0*(u**2.0d0 + v**2.0d0))

            !        write(fo,*)i*(180.0d0/100.0d0)-90.0d0,T1-T0

                    !u   = Q(i-1,j,2)/Q(i-1,j,1)
                    !v   = Q(i-1,j,3)/Q(i-1,j,1)
                    !Vn2 = u*sin(beta) - v*cos(beta)
                    !write(fo,*)i*dx-0.5d0,Vn,abs(Vn-Vn2),1.0d0
            !    enddo
            !    close(fo)
            !endif
        enddo
        
        !output exact solution
        !fo = fo+1
        !hoge =0.0d0
        !write(filename,'("exactRamp",f3.1,".dat")')theta*180d0/PI
        !open(fo,file = filename)
        !write(fo,*)'x,',' ','Exact,',' ','Exact,',' ',' Exact,',' ','Exact'
        !do i=1,Nx
        !    if(i<=30)then
        !        p=0
        !        !write(fo,*)hoge/Nx-0.5d0,',',1.0d0,',',Minf*sin(beta),',',Minf*cos(beta),',',1.0d0/1.4d0
        !        write(fo,*)hoge/Nx-0.5d0,',',p
        !    else
        !        rho = (((GAMMA + 1.0d0)*(Minf*sin(beta))**2.0d0)/((GAMMA-1.0d0)*(Minf*sin(beta))**2.0d0 + 2.0d0)) * rhoinf
        !        p   = ((2.0d0*GAMMA*((Minf*sin(beta))**2.0d0))/(GAMMA + 1.0d0) - (GAMMA - 1.0d0)/(GAMMA + 1.0d0)) * pinf
        !        Vn  = ((GAMMA - 1.0d0)/(GAMMA + 1.0d0) + 2.0d0/((GAMMA + 1.0d0)*((Minf*sin(beta))**2.0d0))) * Minf*sin(beta)
        !        Vt  = Minf*cos(beta)
        !        p=1
        !    write(fo,*)hoge/Nx-0.5d0,',',p
        !    !write(fo,*)hoge/Nx-0.5d0,',',rho,',',Vn,',',Vt,',',p
        !    endif
        !    hoge = hoge + 1.0d0
        !enddo
        !close(fo)


end program pltdist