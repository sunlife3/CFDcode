module plot
implicit none

contains
    subroutine print_profile(Q,Nx,Ny,cnt)
        real(8),intent(in) :: Q(-2:Nx+2,-2:Ny+2,4)
        real(8) dx,beta
        integer,intent(in) :: Nx,Ny,cnt
        integer i,j
        character(len=*),parameter :: base1 = "./TimeSeriesDat_rho.",base2 = "./TimeSeriesDat_Vt."
        character(len=2) :: SerialNum
        character(len=40) filename

        write(SerialNum,'(i2.2)') cnt
        dx = 1.0d0/Nx
        beta = 7.2d0
        j=1

        write(filename,"(a,i3.3,a)") base1//SerialNum
        open(256,file=filename)
        do i = 1, Nx
           write(256,*) i*dx-0.5d0,Q(i,j,1) 
        enddo
        close(256)

        write(filename,"(a,i3.3,a)") base2//SerialNum
        open(257,file=filename)
        do i = 1, Nx
           write(257,*) i*dx-0.5d0,Q(i,j,2)/Q(i,j,1)*cos(beta) + Q(i,j,3)/Q(i,j,1)*sin(beta) 
        enddo
        close(257)


    end subroutine print_profile

    
    subroutine AnimationScript
        integer :: cnt,FileNum=20
        character(len=*),parameter :: base = "./TimeSeriesDat_rho."
        character(len=2) :: SerialNum

        print *, "#"
        print *, "# gnuplot script generated by subroutine AnimationScript"
        print *, "#"
        print *, " "
        print *, "set xlabel 'x'"
        print *, "set ylabel 'Vt'"       !something else(rho,Vn,Vt,p,...)
        print *, "set xrange[-0.5:0.5]"
        print *, "set yrange[19.82:19.85]"            ! Be changed by what you plot
        print *, "set xlabel font 'Arial,20' "
        print *, "set ylabel font 'Arial,20' "
        print *, "set tics font 'Arial,20' "

        do cnt = 0,FileNum
            write(SerialNum,'(i3.3)')cnt
            print *,"plot '"//base//SerialNum//"' w lp"

            if(cnt == 0)then
                print *,"pause 0.2"
            else
                print *,"pause 0.1"
            end if
        enddo
   
    end subroutine AnimationScript

    subroutine output2file(Q,Nx,Ny,index,meshfile)
        implicit none
        integer,intent(in) :: Nx,Ny,index
        integer :: fo=20,i=0,j=0,Grid=19
        real(8),intent(in) :: Q(-2:Nx+2,-2:Ny+2,4)
        character filename*128
        character,intent(in) :: meshfile*64

        !write(filename,'("Qbin HLLAUSM Cylinder(M20,100)_",i3.3,".dat")')index
        !write(filename,'("QbinSLAU2 SNS(M20,Buff,3timesLT)_",i3.3,".dat")')index
        write(filename,'("Qbin_AUSM_duct5deg(M15,100)_",i1.1,".dat")')index
        !write(filename,'("Qbin_HLLAUSM_ramp5deg(M15,100,case1)_",i1.1,".dat")')index
        !write(filename,'("Qbin_HLL_1stratio_50(HLL)&25(HLL)&0(HLL)_",i1.1,".dat")')index
        !write(filename,'("Qbin_AUSM_Shu-Osher3rd_50(AUSM)&25(AUSM)&0(AUSM)_",i1.1,".dat")')index
        !write(filename,'("Qbin_SD_DoubleMach(M5,200)_",i1.1,".dat")')index
        !write(filename,'("Qbin_Exact_ramp5deg(M15,100)_",i1.1,".dat")')index
        !write(filename,'("CSP.dat")')
        !write(filename,'("Qbin_SWBLI(HLLAUSM,3000,Re1000)_",i2.2,".dat")')index
        !open(fo,file='Qbin.dat')
        open(fo,file = filename)
        write(fo,*)meshfile
        do j=-2,Ny+2
            do i=-2,Nx+2
                write(fo,*)i,j,Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            enddo
        enddo
        close(fo)
    end subroutine output2file

    subroutine pltdist(Q,Nx,Ny,meshfile)
        implicit none
        integer :: fo=20,i=0,j=0,GridNum=18,Grid=19,Qbin=100,index=1,is_cyl=0
        integer,intent(in) :: Nx,Ny
        character filenameraw*128,meshfile*64, filename*128
        real(8),allocatable :: X(:,:),Y(:,:)
        real(8),intent(in) :: Q(-2:Nx+2,-2:Ny+2,4)
        real(8) :: rho,u,v,p,t,c,xcenter,ycenter,beta,theta,rhoinf,uinf,pinf,Minf,Vn,Vt,dx,hoge,pp,pm,&
                    rho1,p1,t1,rho2,p2,t2,Sx,Sy
        !beta ... shock wave angle theta ... deflection angle
        real(8),parameter :: GAMMA=1.4d0,PI = 3.14159265359d0


        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4))
        
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
        !beta = 7.2d0*PI/180d0!M15
        !beta = 8.1d0*PI/180d0!M8
        !beta = 9.5d0*PI/180d0!M10
        !beta = 10.9d0*PI/180d0!M8
        beta = 11.8d0*PI/180d0!M7
        !beta = 15.1d0*PI/180d0!M5
        !beta = 18.0d0*PI/180d0!M4
        !beta = 23.1d0*PI/180d0!M3
        !beta = 34.3d0*PI/180d0!M2
        theta = 5.0d0*PI/180d0
        
        do index=1,1

            open(Grid,file = meshfile)

            do j=-3,Ny+4
                do i=-3,Nx+4
                    read(Grid,*)X(i,j),Y(i,j)
                enddo
            enddo
            close(Grid)

            
            close(Qbin)
            !write(filename,'("ramp5deg_SD-AUSMDV_(M15)(100)",i1.1,".dat")')index
            !write(filename,'("ramp5deg_HLLAUSMmix_(M15,Case3-2,triple)(100)",i1.1,".dat")')index
            write(filename,'("ramp5deg_HLLAUSM(M15,100,case1)",i3,".dat")')index
            !write(filename,'("AUSM_Shu-Osher3rd_50(AUSM)&25(AUSM)&0(AUSM)_",i1.1,".dat")')index
            !write(filename,'("duct5deg_SD-AUSMDV_(M15)",i1.1,".dat")')index
            !write(filename,'("SNS_SD-AUSMDV_(M20,SFver2)",i1.1,".dat")')index
            !write(filename,'("ramp5 grad Vn SD-AUSMDV_(M15,C,100)",i1.1,".dat")')index
            !write(filename,'("Cylinder_SLAU2(M20,100)_",i1.1,".dat")')index
            !write(filename,'("LAX_AUSMDV_(M15,1st)(100)",i1.1,".dat")')index
            open(fo,file = filename)
            write(fo,*)'x,',' ','HLLAUSM,',' ','HLLAUSM,',' ','HLLAUSM,',' ','HLLAUSM'
            !write(fo,*)'x,',' ','AUSMDV,',' ','AUSMDV,',' ',' AUSMDV,',' ','AUSMDV'
            !write(fo,*)'x',' ','Vn',' ','Vn_i-Vn_(i-1)',' ','threshold'
            !write(fo,*)'x,',' ','SD-SLAU,',' ','SD-SLAU,',' ',' SD-SLAU,',' ','SD-SLAU'
            !write(fo,*)'x,',' ','SLAU2,',' ','SLAU2,',' ','SLAU2,',' ','SLAU2'
            !write(fo,*)'x,',' ','HLL(1st),',' ','HLL(1st),',' ','HLL(1st),',' ','HLL(1st)'
            !write(fo,*)'x,',' ',&
            !'SLAU2[1st:3rd=2:8](Shock)SLAU2[1st:3rd=1:9](Around)AUSMDV3rd(Other),',' ',&
            !'SLAU2[1st:3rd=2:8](Shock)SLAU2[1st:3rd=1:9](Around)AUSMDV3rd(Other),',' ',&
            !'SLAU2[1st:3rd=2:8](Shock)SLAU2[1st:3rd=1:9](Around)AUSMDV3rd(Other),',' ',&
            !'SLAU2[1st:3rd=2:8](Shock)SLAU2[1st:3rd=1:9](Around)AUSMDV3rd(Other)'
            !write(fo,*)'x,','rho,','Vn,','Vt,','p,'
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

                rho1 = Q(i,j,1)
                p1   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
                t1   = p1/rho1

                rho2 = Q(i+1,j,1)
                p2   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                t2   = p2/rho2
                Sx = p/(rho**GAMMA) !log((t2/t1)**(GAMMA/(GAMMA-1.0d0))*((p1**rho1)/(p2**rho2)))
                !write(*,*)i,log((p2/p1)*(rho1/rho2)**1.40d0)!

                rho2 = Q(i,j+1,1)
                p2   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
                t2   = p2/rho2
                Sy = p1/(rho1**GAMMA)!p2/(rho2**GAMMA) - p1/(rho1**GAMMA)!log((t1/t2)**(GAMMA/(GAMMA-1.0d0))*((p2**rho2)/(p1**rho1)))

                write(fo,*)X(i,j)-0.01d0,',',rho,',',Vn,',',Vt,',',p
                !write(fo,*)X(i,j)-5.0d0,',',rho,',',u,',',v,',',p

                
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
        !write(filename,'("exactRamp",f3.1,"_M",f4.1,".dat")')theta*180d0/PI,Minf
        !open(fo,file = filename)
        !write(fo,*)'x,',' ','Exact,',' ','Exact,',' ',' Exact,',' ','Exact'
        !do i=1,Nx
        !    if(i<=20)then
        !        write(fo,*)hoge/Nx-0.5d0,',',1.0d0,',',Minf*sin(beta),',',Minf*cos(beta),',',1.0d0/1.4d0
        !    else
        !        rho = (((GAMMA + 1.0d0)*(Minf*sin(beta))**2.0d0)/((GAMMA-1.0d0)*(Minf*sin(beta))**2.0d0 + 2.0d0)) * rhoinf
        !        p   = ((2.0d0*GAMMA*((Minf*sin(beta))**2.0d0))/(GAMMA + 1.0d0) - (GAMMA - 1.0d0)/(GAMMA + 1.0d0)) * pinf
        !        Vn  = ((GAMMA - 1.0d0)/(GAMMA + 1.0d0) + 2.0d0/((GAMMA + 1.0d0)*((Minf*sin(beta))**2.0d0))) * Minf*sin(beta)
        !        Vt  = Minf*cos(beta)
        !    write(fo,*)hoge/Nx-0.5d0,',',rho,',',Vn,',',Vt,',',p
        !    endif
        !    hoge = hoge + 1.0d0
        !enddo
        !close(fo)


    end subroutine pltdist
end module plot