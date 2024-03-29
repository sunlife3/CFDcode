
    subroutine FixedInFlow(Q,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer j
        real(8) rho,u,v,p
        real(8),intent(out) :: Q(-2:Nx+3,-2:Ny+3,4)
        !rho = 1.0d0
        !p   = 1.0/1.4d0
        !u   = 20.0d0
        !v   = 0.0d0
        !!! Shu-osher
        !rho = 3.857143
        !p   = 31.0/3.0d0
        !u   = 2.629369
        !v   = 0.0d0
        !!! LAX
        !rho = 0.445d0
        !u = 0.698d0
        !v = 0.0d0
        !p = 3.528d0
        !!! SOD
        rho = 1.0d0
        u = 20.0d0 * cos(30.0d0*3.141592d0/180.0d0)
        v = 20.0d0 * sin(30.0d0*3.141592d0/180.0d0)
        p = 1.0d0/1.4d0
        do j=-2,Ny+3
            Q(0,j,1) = rho
            Q(0,j,2) = rho*u
            Q(0,j,3) = rho*v
            Q(0,j,4) = p/(1.4d0 - 1.0d0) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)
            !
            Q(-1,j,1) = rho
            Q(-1,j,2) = rho*u
            Q(-1,j,3) = rho*v
            Q(-1,j,4) = p/(1.4d0 - 1.0d0) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

            Q(-2,j,1) = rho
            Q(-2,j,2) = rho*u
            Q(-2,j,3) = rho*v
            Q(-2,j,4) = p/(1.4d0 - 1.0d0) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)
        enddo
    end subroutine FixedInFlow

    subroutine FreeOutFlow(Q,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer j
        real(8),intent(out) :: Q(-2:Nx+3,-2:Ny+3,4)
        do j=-2,Ny+3
            Q(Nx+1,j,1) = Q(Nx,j,1)
            Q(Nx+1,j,2) = Q(Nx,j,2)
            Q(Nx+1,j,3) = Q(Nx,j,3)
            Q(Nx+1,j,4) = Q(Nx,j,4)
            !
            Q(Nx+2,j,1) = Q(Nx,j,1)
            Q(Nx+2,j,2) = Q(Nx,j,2)
            Q(Nx+2,j,3) = Q(Nx,j,3)
            Q(Nx+2,j,4) = Q(Nx,j,4)
            !
            Q(Nx+3,j,1) = Q(Nx,j,1)
            Q(Nx+3,j,2) = Q(Nx,j,2)
            Q(Nx+3,j,3) = Q(Nx,j,3)
            Q(Nx+3,j,4) = Q(Nx,j,4)
        enddo
    end subroutine FreeOutFlow

    subroutine setbdycond(Q,m,n,Nx,Ny,is_SF_xi,argc)
        implicit none
        integer,intent(in) :: Nx,Ny,is_SF_xi(-1:Nx+2,-1:Ny+2)
        integer :: i=0,j=0
        character,intent(in) :: argc*64
        real(8),intent(in)  :: m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2)
        real(8),intent(out) :: Q(-2:Nx+3,-2:Ny+3,4)
        real(8),parameter :: PI = acos(-1.0d0),GAMMA = 1.4d0
        real(8) :: Ut=0,Un=0,u=0,v=0,rho,p,Minf,nxctr,nyctr,mxctr,myctr,&
                    rhoNx,uNx,vNx,pNx

        !do j=-1,Ny+2
        !    do i=-1,Nx+2
        !        write(*,*)i,j,Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
        !    enddo
        !enddo
        
        
        
        !!!!!!!!!!!!!!!!!!!      LEFT TERMINAL AND RIGHT TERMINAL      !!!!!!!!!!!!!!!!!!!!!!
        if(argc == "SNS")then
            Minf = 20d0
            !write(*,*)"SNS BdyCond"
            do j=-2,Ny+3!B.C.@i=-1,0,Nx+1,Nx+2
                !!!!! left terminal : Fixed Inflow(Steady Normal Shock) !!!!!!!
                Q(0,j,1) = 1.0d0
                Q(0,j,2) = 1.0d0
                Q(0,j,3) = 0.0d0
                Q(0,j,4) = 1.0d0/(1.4d0*(1.4d0 - 1.0d0)*Minf**2.0d0) + 0.5d0
                
                Q(-1,j,1) = 1.0d0
                Q(-1,j,2) = 1.0d0
                Q(-1,j,3) = 0.0d0
                Q(-1,j,4) = 1.0d0/(1.4d0*(1.4d0 - 1.0d0)*Minf**2.0d0) + 0.5d0
                !!!!! Right terminal : Fixed mass flux (Steady Normal Shock) !!!!!!!
                Q(Nx+1,j,1) = Q(Nx,j,1)
                Q(Nx+1,j,2) = 1.0d0!Q(Nx,j,2)
                Q(Nx+1,j,3) = Q(Nx,j,3)
                rhoNx = Q(Nx,j,1)
                uNx   = 1.0d0/rhoNx
                vNx   = Q(Nx,j,3)/rhoNx
                pNx   = (GAMMA - 1.0d0)*( Q(Nx,j,4) - 0.5d0* ( Q(Nx,j,2)**2.0d0 + Q(Nx,j,3)**2.0d0)/Q(Nx,j,1))
                Q(Nx+1,j,4) = pNx/(GAMMA - 1.0d0) + 0.5d0*rhoNx*(uNx**2.0d0 + vNx**2.0d0)
                
                Q(Nx+2,j,1) = Q(Nx,j,1)
                Q(Nx+2,j,2) = 1.0d0!Q(Nx,j,2)
                Q(Nx+2,j,3) = Q(Nx,j,3)
                Q(Nx+2,j,4) = pNx/(GAMMA - 1.0d0) + 0.5d0*rhoNx*(uNx**2.0d0 + vNx**2.0d0)
            enddo
            
        
        else if(argc == "PS")then
            !!!! For Propagating shock !!!!!!!!!!!
            do j=1,Ny
                Q(0,j,1) = ((1.4d0 + 1.0d0)*Minf**2.0d0)/((1.4d0 - 1.0d0)*Minf**2.0d0 + 2.0d0)
                Q(0,j,2) = 0.0d0
                Q(0,j,3) = 0.0d0
                Q(0,j,4) = ((2.0d0*1.4d0*Minf**2.0d0)-(1.4d0 - 1.0d0))/(1.4d0 +1.0d0)*(1.0d0/1.4d0)
                
                Q(-1,j,1) = ((1.4d0 + 1.0d0)*Minf**2.0d0)/((1.4d0 - 1.0d0)*Minf**2.0d0 + 2.0d0)
                Q(-1,j,2) = 0.0d0
                Q(-1,j,3) = 0.0d0
                Q(-1,j,4) = ((2.0d0*1.4d0*Minf**2.0d0)-(1.4d0 - 1.0d0))/(1.4d0 +1.0d0)*(1.0d0/1.4d0)
                !call FixedInFlow(Q,Nx,Ny)
                call FreeOutFlow(Q,Nx,Ny)
                
            enddo
        
        else if(argc == "SWBLI")then
            !!!!!!! Reflect wall (left and right)!!!!!!!!!!!
            do j=1,Ny
                Q(0,j,1) = Q(1,j,1)
                Q(0,j,2) = -Q(1,j,2)
                Q(0,j,3) = -Q(1,j,3)
                Q(0,j,4) = Q(1,j,4)
                
                Q(-1,j,1) = Q(2,j,1)
                Q(-1,j,2) = -Q(2,j,2)
                Q(-1,j,3) = -Q(2,j,3)
                Q(-1,j,4) = Q(2,j,4)
                
                Q(Nx+1,j,1) = Q(Nx,j,1)
                Q(Nx+1,j,2) = -Q(Nx,j,2)
                Q(Nx+1,j,3) = -Q(Nx,j,3)
                Q(Nx+1,j,4) = Q(Nx,j,4)
                
                Q(Nx+2,j,1) = Q(Nx-1,j,1)
                Q(Nx+2,j,2) = -Q(Nx-1,j,2)
                Q(Nx+2,j,3) = -Q(Nx-1,j,3)
                Q(Nx+2,j,4) = Q(Nx-1,j,4)
                
            enddo
        
        else if(argc == "Bow")then

            call FixedInFlow(Q,Nx,Ny)
            do j=-2,Ny+3
                mxctr = 0.5d0*(m(Nx,j,1) + m(Nx-1,j,1))
                myctr = 0.5d0*(m(Nx,j,2) + m(Nx-1,j,2))
                Un = Q(Nx,j,2)/Q(Nx,j,1)*mxctr + Q(Nx,j,3)/Q(Nx,j,1)*myctr
                Ut = -Q(Nx,j,2)/Q(Nx,j,1)*myctr + Q(Nx,j,3)/Q(Nx,j,1)*mxctr
                mxctr = 0.5d0*(m(Nx+1,j,1) + m(Nx,j,1))
                myctr = 0.5d0*(m(Nx+1,j,2) + m(Nx,j,2))
                u  = (-mxctr*Un - myctr*Ut)/(mxctr**2.0d0 + myctr**2.0d0)
                v  = (-myctr*Un + mxctr*Ut)/(mxctr**2.0d0 + myctr**2.0d0)
                Q(Nx+1,j,1) = Q(Nx,j,1)
                Q(Nx+1,j,4) = Q(Nx,j,4)
                Q(Nx+1,j,2) = Q(Nx+1,j,1)*u
                Q(Nx+1,j,3) = Q(Nx+1,j,1)*v
                !write(*,*)Nx,j,Q(Nx,j,1),Q(Nx,j,2),Q(Nx,j,3),Q(Nx,j,4)
                !write(*,*)Nx+1,j,Q(Nx+1,j,1),Q(Nx+1,j,2),Q(Nx+1,j,3),Q(Nx+1,j,4) 

                mxctr = 0.5d0*(m(Nx-1,j,1) + m(Nx-2,j,1))
                myctr = 0.5d0*(m(Nx-1,j,2) + m(Nx-2,j,2))
                Un = Q(Nx-1,j,2)/Q(Nx-1,j,1)*mxctr + Q(Nx-1,j,3)/Q(Nx-1,j,1)*myctr
                Ut = -Q(Nx-1,j,2)/Q(Nx-1,j,1)*myctr + Q(Nx-1,j,3)/Q(Nx-1,j,1)*mxctr
                mxctr = 0.5d0*(m(Nx+2,j,1) + m(Nx+1,j,1))
                myctr = 0.5d0*(m(Nx+2,j,2) + m(Nx+1,j,2))
                u  = (-mxctr*Un - myctr*Ut)/(mxctr**2.0d0 + myctr**2.0d0)
                v  = (-myctr*Un + mxctr*Ut)/(mxctr**2.0d0 + myctr**2.0d0)
                Q(Nx+2,j,1) = Q(Nx-1,j,1)
                Q(Nx+2,j,4) = Q(Nx-1,j,4)
                Q(Nx+2,j,2) = Q(Nx+2,j,1)*u
                Q(Nx+2,j,3) = Q(Nx+2,j,1)*v
            enddo
        else
            call FixedInFlow(Q,Nx,Ny)
            call FreeOutFlow(Q,Nx,Ny)
        endif

        !!!!!!!!!!!!!!!!!!!      UPSIDE TEMINAL AND DOWNSIDE TERMINAL      !!!!!!!!!!!!!!!!!!!!!!
        if(argc == "Bow")then
            do i=-2,Nx+3!B.C.@i=-1,0,Nx+1,Nx+2
                !!!!! downside terminal : free Outflow !!!!!!!!!!!!!
                Q(i,0,1) = Q(i,1,1)
                Q(i,0,2) = Q(i,1,2)
                Q(i,0,3) = Q(i,1,3)
                Q(i,0,4) = Q(i,1,4)

                Q(i,-1,1) = Q(i,1,1)
                Q(i,-1,2) = Q(i,1,2)
                Q(i,-1,3) = Q(i,1,3)
                Q(i,-1,4) = Q(i,1,4)
                !write(*,*)i,j,Q(i,0,1),Q(i,0,2),Q(i,0,3),Q(i,0,4)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !!!!! upside terminal : free OutFlow !!!!!!
                Q(i,Ny+1,1) = Q(i,Ny,1)
                Q(i,Ny+1,2) = Q(i,Ny,2)
                Q(i,Ny+1,3) = Q(i,Ny,3)
                Q(i,Ny+1,4) = Q(i,Ny,4)

                Q(i,Ny+2,1) = Q(i,Ny,1)
                Q(i,Ny+2,2) = Q(i,Ny,2)
                Q(i,Ny+2,3) = Q(i,Ny,3)
                Q(i,Ny+2,4) = Q(i,Ny,4)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo
        else
            do i=-2,Nx+3!B.C.@j=-1,0,Ny+1,Ny+2
                !!!!!!!!  SLIP CONDITION  !!!!!!!!
                nxctr = 0.5d0*(n(i,0,1) + n(i,1,1))
                nyctr = 0.5d0*(n(i,0,2) + n(i,1,2))
                Un =  Q(i,1,2)/Q(i,1,1)*nxctr + Q(i,1,3)/Q(i,1,1)*nyctr
                Ut = -Q(i,1,2)/Q(i,1,1)*nyctr + Q(i,1,3)/Q(i,1,1)*nxctr
                nxctr = 0.5d0*(n(i,0,1) + n(i,-1,1))
                nyctr = 0.5d0*(n(i,0,2) + n(i,-1,2))
                u  = (-nxctr*Un - nyctr*Ut)/(nxctr**2 + nyctr**2)
                v  = (-nyctr*Un + nxctr*Ut)/(nxctr**2 + nyctr**2)
                Q(i,0,1) = Q(i,1,1)
                Q(i,0,4) = Q(i,1,4)
                Q(i,0,2) = Q(i,0,1)*u
                Q(i,0,3) = Q(i,0,1)*v

                nxctr = 0.5d0*(n(i,1,1) + n(i,2,1))
                nyctr = 0.5d0*(n(i,1,2) + n(i,2,2))
                Un =  Q(i,2,2)/Q(i,2,1)*nxctr + Q(i,2,3)/Q(i,2,1)*nyctr
                Ut = -Q(i,2,2)/Q(i,2,1)*nyctr + Q(i,2,3)/Q(i,2,1)*nxctr
                nxctr = 0.5d0*(n(i,-1,1) + n(i,-2,1))
                nyctr = 0.5d0*(n(i,-1,2) + n(i,-2,2))
                u  = (-nxctr*Un - nyctr*Ut)/(nxctr**2 + nyctr**2)
                v  = (-nyctr*Un + nxctr*Ut)/(nxctr**2 + nyctr**2)
                Q(i,-1,1) = Q(i,2,1)
                Q(i,-1,4) = Q(i,2,4)
                Q(i,-1,2) = Q(i,-1,1)*u
                Q(i,-1,3) = Q(i,-1,1)*v

                nxctr = 0.5d0*(n(i,2,1) + n(i,3,1))
                nyctr = 0.5d0*(n(i,2,2) + n(i,3,2))
                Un =  Q(i,3,2)/Q(i,3,1)*nxctr + Q(i,3,3)/Q(i,3,1)*nyctr
                Ut = -Q(i,3,2)/Q(i,3,1)*nyctr + Q(i,3,3)/Q(i,3,1)*nxctr
                nxctr = 0.5d0*(n(i,-1,1) + n(i,-2,1))
                nyctr = 0.5d0*(n(i,-1,2) + n(i,-2,2))
                u  = (-nxctr*Un - nyctr*Ut)/(nxctr**2 + nyctr**2)
                v  = (-nyctr*Un + nxctr*Ut)/(nxctr**2 + nyctr**2)
                Q(i,-2,1) = Q(i,3,1)
                Q(i,-2,4) = Q(i,3,4)
                Q(i,-2,2) = Q(i,-2,1)*u
                Q(i,-2,3) = Q(i,-2,1)*v

                nxctr = 0.5d0*(n(i,Ny,1) + n(i,Ny-1,1))
                nyctr = 0.5d0*(n(i,Ny,2) + n(i,Ny-1,2))
                Un =  Q(i,Ny,2)/Q(i,Ny,1)*nxctr + Q(i,Ny,3)/Q(i,Ny,1)*nyctr
                Ut = -Q(i,Ny,2)/Q(i,Ny,1)*nyctr + Q(i,Ny,3)/Q(i,Ny,1)*nxctr
                nxctr = 0.5d0*(n(i,Ny,1) + n(i,Ny+1,1))
                nyctr = 0.5d0*(n(i,Ny,2) + n(i,Ny+1,2))
                u  = (-nxctr*Un - nyctr*Ut)/(nxctr**2 + nyctr**2)
                v  = (-nyctr*Un + nxctr*Ut)/(nxctr**2 + nyctr**2)
                Q(i,Ny+1,1) = Q(i,Ny,1)
                Q(i,Ny+1,4) = Q(i,Ny,4)
                Q(i,Ny+1,2) = Q(i,Ny+1,1)*u
                Q(i,Ny+1,3) = Q(i,Ny+1,1)*v

                nxctr = 0.5d0*(n(i,Ny-1,1) + n(i,Ny-2,1))
                nyctr = 0.5d0*(n(i,Ny-1,2) + n(i,Ny-2,2))
                Un =  Q(i,Ny-1,2)/Q(i,Ny-1,1)*nxctr + Q(i,Ny-1,3)/Q(i,Ny-1,1)*nyctr
                Ut = -Q(i,Ny-1,2)/Q(i,Ny-1,1)*nyctr + Q(i,Ny-1,3)/Q(i,Ny-1,1)*nxctr
                nxctr = 0.5d0*(n(i,Ny+2,1) + n(i,Ny+1,1))
                nyctr = 0.5d0*(n(i,Ny+2,2) + n(i,Ny+1,2))
                u  = (-nxctr*Un - nyctr*Ut)/(nxctr**2 + nyctr**2)
                v  = (-nyctr*Un + nxctr*Ut)/(nxctr**2 + nyctr**2)
                Q(i,Ny+2,1) = Q(i,Ny-1,1)
                Q(i,Ny+2,4) = Q(i,Ny-1,4)
                Q(i,Ny+2,2) = Q(i,Ny+2,1)*u
                Q(i,Ny+2,3) = Q(i,Ny+2,1)*v

                nxctr = 0.5d0*(n(i,Ny-2,1) + n(i,Ny-3,1))
                nyctr = 0.5d0*(n(i,Ny-2,2) + n(i,Ny-3,2))
                Un =  Q(i,Ny-2,2)/Q(i,Ny-2,1)*nxctr + Q(i,Ny-2,3)/Q(i,Ny-2,1)*nyctr
                Ut = -Q(i,Ny-2,2)/Q(i,Ny-2,1)*nyctr + Q(i,Ny-2,3)/Q(i,Ny-2,1)*nxctr
                nxctr = 0.5d0*(n(i,Ny+3,1) + n(i,Ny+2,1))
                nyctr = 0.5d0*(n(i,Ny+3,2) + n(i,Ny+2,2))
                u  = (-nxctr*Un - nyctr*Ut)/(nxctr**2 + nyctr**2)
                v  = (-nyctr*Un + nxctr*Ut)/(nxctr**2 + nyctr**2)
                Q(i,Ny+3,1) = Q(i,Ny-2,1)
                Q(i,Ny+3,4) = Q(i,Ny-2,4)
                Q(i,Ny+3,2) = Q(i,Ny+3,1)*u
                Q(i,Ny+3,3) = Q(i,Ny+3,1)*v
            enddo
        endif
    endsubroutine setbdycond
