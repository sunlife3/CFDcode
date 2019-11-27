    subroutine setbdycond(Q,m,n,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0
        real(8),intent(in)  :: m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2)
        real(8),intent(out) :: Q(-1:Nx+2,-1:Ny+2,4)
        real(8) :: Ut=0,Un=0,u=0,v=0,rho,p,Minf=2.0d0,nxctr,nyctr

        !do j=-1,Ny+2
        !    do i=-1,Nx+2
        !        write(*,*)i,j,Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
        !    enddo
        !enddo
        rho = 1.0d0
        p   = 1.0d0/1.4d0
        u   = 2.0d0
        v   = 0.0d0

        do j=1,Ny!B.C.@i=-1,0,Nx+1,Nx+2
            !!!!! left terminal : Inflow !!!!!!!!!!!!!
            Q(0,j,1) = rho
            Q(0,j,2) = rho*u
            Q(0,j,3) = rho*v
            Q(0,j,4) = p/(1.4d0-1.0d0) + 0.5d0*(rho*(u**2 + v**2))
            
            Q(-1,j,1) = rho
            Q(-1,j,2) = rho*u
            Q(-1,j,3) = rho*v
            Q(-1,j,4) = p/(1.4d0-1.0d0) + 0.5d0*(rho*(u**2 + v**2))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!! left terminal : Reflection wall !!!!
            !Q(0,j,1) = Q(1,j,1)
            !Q(0,j,2) = -Q(1,j,2)
            !Q(0,j,3) = Q(1,j,3)
            !Q(0,j,4) = Q(1,j,4)

            !Q(-1,j,1) = Q(2,j,1)
            !Q(-1,j,2) = -Q(2,j,2)
            !Q(-1,j,3) = Q(2,j,3)
            !Q(-1,j,4) = Q(2,j,4)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            !!!!! Right terminal : free OutFlow !!!!!!
            Q(Nx+1,j,1) = Q(Nx,j,1)
            Q(Nx+1,j,2) = Q(Nx,j,2)
            Q(Nx+1,j,3) = Q(Nx,j,3)
            Q(Nx+1,j,4) = Q(Nx,j,4)

            Q(Nx+2,j,1) = Q(Nx,j,1)
            Q(Nx+2,j,2) = Q(Nx,j,2)
            Q(Nx+2,j,3) = Q(Nx,j,3)
            Q(Nx+2,j,4) = Q(Nx,j,4)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!Right temrminal : Reflection wall !!!!!!!!!!
            !Q(Nx+1,j,1) = Q(Nx,j,1)
            !Q(Nx+1,j,2) = -Q(Nx,j,2)
            !Q(Nx+1,j,3) = Q(Nx,j,3)
            !Q(Nx+1,j,4) = Q(Nx,j,4)

            !Q(Nx+2,j,1) = Q(Nx-1,j,1)
            !Q(Nx+2,j,2) = -Q(Nx-1,j,2)
            !Q(Nx+2,j,3) = Q(Nx-1,j,3)
            !Q(Nx+2,j,4) = Q(Nx-1,j,4)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo

        do i=-1,Nx+2!B.C.@j=-1,0,Ny+1,Ny+2
            !SLIP CONDITION
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
            !write(*,*)i,Q(i,0,2),Q(i,1,2),Q(i,0,3),Q(i,1,3)
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
            !write(*,*)i,Q(i,-1,2),Q(i,2,2),Q(i,-1,3),Q(i,2,3) 
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
            !write(*,*)i,Q(i,Ny+1,1),Q(i,Ny+1,2),Q(i,Ny+1,3),Q(i,Ny+1,4) 
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
            !write(*,*)i,Q(i,Ny+2,1),Q(i,Ny+2,2),Q(i,Ny+2,3),Q(i,Ny+2,4) 
        enddo
    end subroutine setbdycond