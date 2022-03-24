    subroutine setInitCond(Q,Nx,Ny,X,Y,argc)
        implicit none
        integer,intent(in) :: Nx,Ny,X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4)
        integer :: i=0,j=0,Qbin=100,index = 1
        character filename*64,meshfile*64
        character,intent(in) :: argc*64
        real(8),intent(out) :: Q(-2:Nx+3,-2:Ny+3,4)
        real(8),parameter :: PI=acos(-1.0d0)
        real(8) :: rho(-2:Nx+3),p(-2:Nx+3),u(-2:Nx+3),v(-2:Nx+3),e,a1,a2,a3,Minf,&
                   rhoL,rhoR,rhoM,uL,uR,uM,pL,pR,pM

        
        if(argc == "SNS")then
            !Steady Normal Shock
            write(*,*)"SNS InitCond"
            Minf=20d0
            do j=-2,Ny+3
                do i=-2,Nx+3
                    if(i <= 40)then
                        Q(i,j,1) = 1.0d0
                        Q(i,j,2) = 1.0d0
                        Q(i,j,3) = 0.0d0
                        Q(i,j,4) = 1.0d0/(1.4d0*(1.4d0 - 1.0d0)*Minf**2.0d0) + 0.5d0
                    else
                        Q(i,j,1) = (2.0d0/(2.4d0*Minf**2.0d0) + 0.4d0/2.4d0)**(-1.0d0)
                        Q(i,j,2) = 1.0d0
                        Q(i,j,3) = 0.0d0
                        Q(i,j,4) = (2.0d0*1.4d0*Minf**2.0d0/2.4d0 - 0.4d0/2.4d0)/(1.4d0*0.4d0*Minf**2.0d0)&
                                    + 1.0d0/(2.0d0*(2.0d0/(2.4d0*Minf**2.0d0) + 0.4d0/2.4d0)**(-1.0d0))
                    endif
                enddo
            enddo
            !e = 1.0d0
            !a1 = e
            !rhoL = Q(12,1,1)
            !rhoR = Q(14,1,1)
            !rhoM = (1.0d0 - a1)*rhoL + a1*rhoR

            !uL = Q(12,1,2)/Q(12,1,1)
            !uR = Q(14,1,2)/Q(14,1,1)
            !a2 = e*(1.0d0 - e)*(1.0d0 + (1.0d0 - e)*(Minf**2.0d0 - 1.0d0)/(1.0d0 + 0.5d0*(1.4d0 - 1.0d0)*Minf**2.0d0))**(-0.5d0)*&
            !                (1.0d0 + (1.0d0 - e)*(Minf**2.0d0 - 1.0d0)/(2.0d0*1.4d0*Minf**2.0d0/(1.4d0 - 1.0d0) - 1.0d0))**(-0.5d0)
            !uM = (1.0d0 - a2)*uL + a2*uR
            !    pL =  (1.4d0 - 1.0d0)*(Q(12,1,4) - 0.5d0*(Q(12,1,2)**2.0d0 + Q(12,1,3)**2.0d0)/Q(12,1,1))
            !pR =  (1.4d0 - 1.0d0)*(Q(14,1,4) - 0.5d0*(Q(14,1,2)**2.0d0 + Q(14,1,3)**2.0d0)/Q(14,1,1))
            !a3 =  1.0d0 - (1.0d0 - e)*(1.0d0 + e*(1.4d0 + 1.0d0)/(1.4d0 - 1.0d0) * (Minf**2.0d0 - 1.0d0)/(Minf**2.0d0))**(-0.5d0)
            !pM = (1.0d0 - a3)*pL + a3*pR
            !do j=-1,Ny+2
            !    Q(13,j,1) = rhoM
            !    Q(13,j,2) = rhoM*uM
            !    Q(13,j,3) = 0.0d0
            !    Q(13,j,4) = pM/(1.4d0 - 1.0d0) + 0.5d0*rhoM*(uM**2.0d0)
            !enddo


        else if(argc == "ramp")then
            !Ramp
            do j=-2,Ny+3
                do i=-2,Nx+3
                    rho(i) = 1.0d0
                    p(i)   = 1.0d0/1.4d0
                    u(i)   = 15.0d0
                    v(i)   = 0.0d0
                    Q(i,j,1) = rho(i)
                    Q(i,j,2) = rho(i)*u(i)
                    Q(i,j,3) = rho(i)*v(i)
                    Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                enddo
            enddo

        else if(argc == "Bow")then
            !Bow
            
            do j=-2,Ny+3
                do i=-2,Nx+3
                    rho(i) = 1.0d0
                    p(i)   = 1.0d0/1.4d0
                    u(i)   = 0.0d0
                    v(i)   = 0.0d0
                    Q(i,j,1) = rho(i)
                    Q(i,j,2) = rho(i)*u(i)
                    Q(i,j,3) = rho(i)*v(i)
                    Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                enddo
            enddo

        else if(argc == "PS")then
            !Propagating Shock Wave
            do j=-2,Ny+3
                do i=-2,Nx+3
                    rho(i) = 1.0d0
                    p(i)   = 1.0d0/1.4d0
                    u(i)   = 3.0d0
                    v(i)   = 0.0d0
                    Q(i,j,1) = rho(i)
                    Q(i,j,2) = rho(i)*u(i)
                    Q(i,j,3) = rho(i)*v(i)
                    Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                enddo
            enddo


        else if(argc == "SO")then
            !Shu Osher
            write(*,*)'SO'
            do j=-2,Ny+3
                do i=-2,Nx+3
                    if(i <= Nx/10)then
                        rho(i) = 3.857143d0
                        p(i)   = 31.0d0/3
                        u(i)   = 2.629369
                        v(i)   = 0.0d0
                    else
                        rho(i) = 1.0d0 + 0.2d0*sin(5.0*i*10.0/Nx)
                        p(i)   = 1.0d0
                        u(i)   = 0.0d0
                        v(i)   = 0.0d0
                    endif
                    Q(i,j,1) = rho(i)
                    Q(i,j,2) = rho(i)*u(i)
                    Q(i,j,3) = rho(i)*v(i)
                    Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                enddo
            enddo

        else if(argc == "SOD")then
            !Shu Osher
            write(*,*)'SOD'
            do j=-2,Ny+3
                do i=-2,Nx+3
                    if(i <= Nx/2)then
                        rho(i) = 1.0d0
                        p(i)   = 1.0d0
                        u(i)   = 0.0d0
                        v(i)   = 0.0d0
                    else
                        rho(i) = 0.125d0
                        p(i)   = 0.1d0
                        u(i)   = 0.0d0
                        v(i)   = 0.0d0
                    endif
                    Q(i,j,1) = rho(i)
                    Q(i,j,2) = rho(i)*u(i)
                    Q(i,j,3) = rho(i)*v(i)
                    Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                enddo
            enddo
        
        else if(argc == "LAX")then
            !LAX problem 
            do j=-2,Ny+3
                do i=-2,Nx+3
                    if(i < Nx/2)then
                        rho(i) = 0.445d0
                        u(i) = 0.698d0
                        v(i) = 0.0d0
                        p(i) = 3.528d0
                        Q(i,j,1) = rho(i)
                        Q(i,j,2) = rho(i)*u(i)
                        Q(i,j,3) = rho(i)*v(i)
                        Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                    else
                        rho(i) = 0.5d0
                        u(i) = 0.0d0
                        v(i) = 0.0d0
                        p(i) = 0.571d0
                        Q(i,j,1) = rho(i)
                        Q(i,j,2) = rho(i)*u(i)
                        Q(i,j,3) = rho(i)*v(i)
                        Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                    endif
                enddo
            enddo
        
        else if(argc == "visplate")then
            !write(*,*) "plate"
            do j=-2,Ny+3
                do i=-2,Nx+3
                    rho(i) = 1.0d0
                    u(i) = 5.0d0
                    v(i) = 0.0d0
                    p(i) = 1.0d0/1.4d0
                    Q(i,j,1) = rho(i)
                    Q(i,j,2) = rho(i)*u(i)
                    Q(i,j,3) = rho(i)*v(i)
                    Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                enddo
            enddo
            !write(filename,'("Qbin_AUSM_BL(M5)_",i1.1,".dat")')9
            !open(Qbin,file = filename)
            !read(Qbin,*)meshfile
            !write(*,*)meshfile
            !do j=-2,Ny+3
            !    do i=-2,Nx+3
            !        read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            !        !write(*,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            !    enddo
            !enddo
            !close(Qbin)

        else if(argc == "SWBLI")then
            !write(filename,'("Qbin_SWBLI_",i1.1,".dat")')9
            !open(Qbin,file = filename)
            !read(Qbin,*)meshfile
            !write(*,*)meshfile
            !do j=-1,Ny+2
            !    do i=-1,Nx+2
            !        read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            !        !write(*,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            !    enddo
            !enddo
            !close(Qbin)
            do j=-2,Ny+3
                do i=-2,Nx+3
                    if(i < Nx/2)then
                        rho(i) = 120d0
                        u(i) = 0.0d0
                        v(i) = 0.0d0
                        p(i) = rho(i)/1.4d0
                        Q(i,j,1) = rho(i)
                        Q(i,j,2) = rho(i)*u(i)
                        Q(i,j,3) = rho(i)*v(i)
                        Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                    else
                        rho(i) = 1.2d0
                        u(i) = 0.0d0
                        v(i) = 0.0d0
                        p(i) = rho(i)/1.4d0
                        Q(i,j,1) = rho(i)
                        Q(i,j,2) = rho(i)*u(i)
                        Q(i,j,3) = rho(i)*v(i)
                        Q(i,j,4) = p(i)/(1.4d0 - 1.0d0) + 0.5d0*rho(i)*(u(i)**2.0d0 + v(i)**2.0d0)
                    endif
                enddo
            enddo
        
        else if(argc == "DMR")then
            !Double Mach Reflection
            do j=-2,Ny+3
                do i=-2,Nx+3
                    if(0.01d0*i <= (1.0d0+0.9d0)/tan(60.0d0*PI/180.0d0))then
                        if(tan(60.0d0*PI/180.0d0)*0.01d0*i - 0.9d0 < 0.01d0*j)then
                            Q(i,j,1) = 8.0d0
                            Q(i,j,2) = 8.0d0*7.145d0
                            Q(i,j,3) = 8.0d0*(-4.125d0)
                            Q(i,j,4) = 116.5d0/(1.4d0 - 1.0d0) + 0.5d0*8.0d0*(7.145d0**2.0d0 + 4.125d0**2.0d0)
                        else
                            Q(i,j,1) = 1.4d0
                            Q(i,j,2) = 0.0d0
                            Q(i,j,3) = 0.0d0
                            Q(i,j,4) = 1.0d0/(1.4d0 - 1.0d0)
                        endif
                    else
                        Q(i,j,1) = 1.4d0
                        Q(i,j,2) = 0.0d0
                        Q(i,j,3) = 0.0d0
                        Q(i,j,4) = 1.0d0/(1.4d0 - 1.0d0)
                    endif
                enddo
            enddo

        endif  

    end subroutine setInitCond
