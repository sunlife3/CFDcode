module param
    implicit none
    real(8) :: GAMMA = 1.4d0,Re = 5000.0d0
contains
    function minmod(x,y) result(fluxLimiter)
        real(8),intent(in) :: x,y
        real(8) fluxLimiter
        fluxLimiter = sign(1.0d0,x)*max(0.0d0,min(abs(x),y*sign(1.0d0,x)))
    end function minmod
end module param

module subprog
    use param
    implicit none
contains
    subroutine setMetrics(m,n,S,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0,Grid=12
        real(8),intent(out) :: m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
                               S(-2:Nx+3,-2:Ny+3)
        real(8) X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4)
        
        open(Grid,file = 'MESH_rampGrid.txt')
        !open(Grid,file = 'MESH_GridVisPlate.txt')
        !open(Grid,file = 'Gridtube.txt')
        do j=-3,Ny+4
            do i=-3,Nx+4
                read(Grid,*)X(i,j),Y(i,j)
            enddo
        enddo
        close(Grid)

        do j=-2,Ny+3
            do i=-2,Nx+3
                m(i,j,1) = Y(i,j) - Y(i,j-1)
                m(i,j,2) = -(X(i,j) - X(i,j-1))
                n(i,j,1) = -(Y(i,j) - Y(i-1,j))
                n(i,j,2) = X(i,j) - X(i-1,j)
                if(j<=2)then
                !write(*,*)m(i,j,1),m(i,j,2),n(i,j,1),n(i,j,2)           
                endif
                S(i,j) = 0.5d0*(n(i,j-1,2)+n(i,j,2))*0.5d0*(m(i-1,j,1)+m(i,j,1))-&
                         0.5d0*(n(i,j-1,1)+n(i,j,1))*0.5d0*(m(i-1,j,2)+m(i,j,2))
            enddo
        enddo
    end subroutine setMetrics


    function calcSpeed(Q,Nx,Ny,m,n,S) result(qmax)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0
        real(8),intent(in) ::Q(-1:Nx+2,-1:Ny+2,4),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
                             S(-2:Nx+3,-2:Ny+3)
        real(8) :: tmp=0,c=0,rho=0,p=0,u=0,v=0,U1=0,U2=0,Cxi=0,Ceta=0,CxiMax=0,CetaMax=0
        real(8) qmax

        do j=1,Ny
            do i=1,Nx
                rho = Q(i,j,1)
                u   = Q(i,j,2)/Q(i,j,1)
                v   = Q(i,j,3)/Q(i,j,1)
                p   = (GAMMA-1)*(Q(i,j,4) - 0.5d0*rho*(u**2+v**2))
                U1  = m(i,j,1)*u + m(i,j,2)*v
                U2  = n(i,j,1)*u + n(i,j,2)*v
                c   = sqrt(GAMMA*p/rho)
                Cxi = 1.0d0/S(i,j)*(abs(U1) + c*sqrt(m(i,j,1)**2 + m(i,j,2)**2))
                Ceta = 1.0d0/S(i,j)*(abs(U2) + c*sqrt(n(i,j,1)**2 + n(i,j,2)**2))

                if(CxiMax < Cxi)then
                    CxiMax = Cxi
                endif
                if(CetaMax < Ceta)then
                    CetaMax = Ceta
                endif
            enddo
        enddo
        qmax = max(CetaMax,CxiMax)
        
    end function calcSpeed

    subroutine output2file(Q,Nx,Ny,index)
        implicit none
        integer,intent(in) :: Nx,Ny,index
        integer :: fo=20,i=0,j=0,Grid=19
        real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4)
        character filename*128

        write(filename,'("QbinSF-SD-SLAU",i3.3,".dat")')index
        !open(fo,file='Qbin.dat')
        open(fo,file = filename)
        do j=-1,Ny+2
            do i=-1,Nx+2
                write(fo,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            enddo
        enddo
        close(fo)
    end subroutine output2file

end module subprog

program main
    use subprog
    use AUSMDVFlux
    use SLAUFlux
    implicit none
    integer :: Nx,Ny,i=0,j=0,k=0,time=0,GridNum=11,foo=12,index=1
    real(8),allocatable :: Q(:,:,:),Qn(:,:,:),Qast(:,:,:),E(:,:,:),F(:,:,:),S(:,:),Ev(:,:,:),Fv(:,:,:)
    real(8),allocatable :: m(:,:,:),n(:,:,:)
    real(8) dt,qmax
    real(8),parameter :: CFL = 0.15d0
    real(8) :: elapsedTime=0.0d0,Q1old=0.0d0,Q2old=0.0d0,res=0.0d0,ressum=0.0d0


    open(GridNum,file = 'GridNum.txt')
    read(GridNum,*)Nx,Ny
    close(GridNum)
    write(*,*)Nx,Ny

    allocate(Q(-1:Nx+2,-1:Ny+2,4),Qn(-1:Nx+2,-1:Ny+2,4),Qast(-1:Nx+2,-1:Ny+2,4),E(-1:Nx+2,-1:Ny+2,4),F(-1:Nx+2,-1:Ny+2,4),&
                S(-2:Nx+3,-2:Ny+3),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
                Ev(-1:Nx+2,-1:Ny+2,4),Fv(-1:Nx+2,-1:Ny+2,4))

    call setInitCond(Q,Nx,Ny)
    call setMetrics(m,n,S,Nx,Ny)
    
    !open(foo,file = 'residualerror.txt')

    do !time=1,1!500!time
        
        call setbdycond(Q,m,n,Nx,Ny)

        call AusmDV_XFlux(Q,E,m,Nx,Ny)
        call AusmDV_YFlux(Q,F,n,Nx,Ny)
        !call SLAU2_Xflux(Q,E,m,Nx,Ny)
        !call SLAU2_YFlux(Q,F,n,Nx,Ny)
        !call Xviscflux(Q,Ev,m,n,S,Nx,Ny)
        !call YviscFlux(Q,Fv,m,n,S,Nx,Ny)

        
        qmax = calcSpeed(Q,Nx,Ny,m,n,S)
        dt =CFL/qmax

        do j=1,Ny
            do i=1,Nx !space
                do k=1,4
                    Qn(i,j,k) = Q(i,j,k)
                    Q(i,j,k) = Q(i,j,k) - 0.50d0*(dt/S(i,j))*(E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k)&          !inertia
                                                        - (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re) !viscous
                enddo
            enddo
        enddo

        do j=-1,Ny+2
            do i=-1,Nx+2
                do k=1,4
                    Qast(i,j,k) = Q(i,j,k)
                enddo
            enddo
        enddo

        !call setbdycond(Qast,m,n,Nx,Ny)
        call AusmDV_XFlux(Qast,E,m,Nx,Ny)
        call AusmDV_YFlux(Qast,F,n,Nx,Ny)
        !call SLAU2_Xflux(Qast,E,m,Nx,Ny)
        !call SLAU2_YFlux(Qast,F,n,Nx,Ny)
        !call Xviscflux(Qast,Ev,m,n,S,Nx,Ny)
        !call YviscFlux(Qast,Fv,m,n,S,Nx,Ny)

        do j=1,Ny
            do i=1,Nx !space
                Q1old = Qn(i,j,1)
                do k=1,4
                    Q(i,j,k) = Qn(i,j,k) - (dt/S(i,j))*(E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k)& 
                                                        - (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re)
                enddo

                res = abs(Q(i,j,1) - Q1old)
                ressum = ressum + res
            enddo
        enddo
        ressum = ressum/(Nx*Ny)
        elapsedTime = elapsedTime + dt
        
        !if(modulo(time,50) == 0)then
        !    write(foo,*)elapsedTime,ressum
        !    call output2file(Q,Nx,Ny,index)
        !    index = index + 1
        !endif

        write(*,*)ressum,elapsedTime,dt

        if(1.2d0<=elapsedTime)then!(ressum <= 0.0000000001d0 .and. 500 <time)then
            exit
        endif
        ressum = 0.0d0
        time = time + 1
    enddo

    close(foo)

    call output2file(Q,Nx,Ny,index)

    deallocate(Q,E,F,S,m,n,Ev,Fv)

end program main