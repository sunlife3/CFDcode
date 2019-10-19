    subroutine setInitCond(Q,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0
        real(8),intent(out) :: Q(-1:Nx+2,-1:Ny+2,4)
        real(8) :: Minf=2.0d0,rho(1:Nx),p(1:Nx),u(1:Nx),v(1:Nx)
        do i=1,Nx
                rho(i) = 1.0d0
                p(i)   = 1.0d0/1.4d0
                u(i)   = 2.0d0
                v(i)   = 0.0d0
        enddo

        do j=1,Ny
            do i=1,Nx
                Q(i,j,1) = rho(i)
                Q(i,j,2) = rho(i)*u(i)
                Q(i,j,3) = rho(i)*v(i)
                Q(i,j,4) = p(i)/(1.40d0-1.0d0) + 0.5d0*rho(i)*((u(i)**2 + v(i)**2))
                
            enddo
        enddo
    end subroutine setInitCond