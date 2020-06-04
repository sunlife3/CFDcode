module param
    implicit none
    real(8) :: GAMMA = 1.4d0,Re = 5000.0d0,CFL = 0.10d0
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
        real(8) X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),yeta,yxi,xeta,xxi
        
        open(Grid,file = 'MESH_cylinder(100100).txt')
        !open(Grid,file = 'MESH_rampGridM10.txt')
        !open(Grid,file = 'MESH_GridVisPlate.txt')
        
        !open(Grid,file = 'MESH_tube(5025).txt')
        do j=-3,Ny+4
            do i=-3,Nx+4
                read(Grid,*)X(i,j),Y(i,j)
            enddo
        enddo
        close(Grid)

        do j=-2,Ny+3
            do i=-2,Nx+3
                m(i,j,1) =  0.5d0*(0.5d0*(Y(i,j+1) - Y(i,j-1)) + 0.5d0*(Y(i+1,j+1) - Y(i+1,j-1))) 
                m(i,j,2) = -0.5d0*(0.5d0*(X(i,j+1) - X(i,j-1)) + 0.5d0*(X(i+1,j+1) - X(i+1,j-1)))
                n(i,j,1) = -0.5d0*(0.5d0*(Y(i+1,j) - Y(i-1,j)) + 0.5d0*(Y(i+1,j+1) - Y(i-1,j+1)))
                n(i,j,2) =  0.5d0*(0.5d0*(X(i+1,j) - X(i-1,j)) + 0.5d0*(X(i+1,j+1) - X(i-1,j+1)))

                S(i,j) = 0.5d0*( (X(i,j) - X(i-1,j-1)) * (Y(i-1,j) - Y(i,j-1)) -&
                         (X(i-1,j) - X(i,j-1)) * (Y(i,j) - Y(i-1,j-1)) )
            enddo
        enddo
    end subroutine setMetrics

    subroutine detectSonicPoint(Q,Nx,Ny,is_SF_xi,is_SF_eta)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0,cnt
        integer,intent(out) :: is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2)
        real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4)
        real(8) :: rho,u,v,p,t,VnL,VnR,cL,cR,beta,f,f1,f2,g
        real(8),parameter :: PI=acos(-1.0d0)

        beta = 0.0d0!9.5d0*PI/180d0

        do j=1,Ny
            do i=1,Nx
                do cnt=0,90
                    beta = beta + 1.0d0

                    VnL = Q(i-2,j,2)/Q(i-2,j,1)*sin(beta) - Q(i-2,j,3)/Q(i-2,j,1)*cos(beta)
                    VnR = Q(i+2,j,2)/Q(i+2,j,1)*sin(beta) - Q(i+2,j,3)/Q(i+2,j,1)*cos(beta)
                    p   = (GAMMA-1.0d0)*(Q(i-2,j,4) - 0.5d0*(Q(i-2,j,2)**2.0d0 + Q(i-2,j,3)**2.0d0)/Q(i-2,j,1))
                    cL = sqrt(GAMMA*p/Q(i-2,j,1))
                    p   = (GAMMA-1.0d0)*(Q(i+2,j,4) - 0.5d0*(Q(i+2,j,2)**2.0d0 + Q(i+2,j,3)**2.0d0)/Q(i+2,j,1))
                    cR = sqrt(GAMMA*p/Q(i+2,j,1))
                    f1 = sign(1.0d0,VnR-cR) * max(abs(VnL/cL),1.0d0)
                    f2 = sign(1.0d0,VnL-cL) * max(abs(VnR/cR),1.0d0)
                    f = min(0.0d0,f1*f2)
                    g = (Q(i+2,j,1) - 2.0d0*Q(i+1,j,1) + Q(i,j,1))*&
                        (Q(i,j,1) - 2.0d0*Q(i-1,j,1) + Q(i-2,j,1))
                    if(abs(f*g) > 0.0d0)then
                       is_SF_xi(i,j) = 1
                       is_SF_xi(i-1,j) = 1
                       is_SF_eta(i,j) = 1
                       is_SF_eta(i,j-1) = 1
                        !write(*,*)i,j
                    endif
            !!!!! detection direction Leftside and Rightside !!!!!!
                    rho = Q(i-1,j,1)
                    u   = Q(i-1,j,2)/Q(i-1,j,1)
                    v   = Q(i-1,j,3)/Q(i-1,j,1)
                    t   = (GAMMA-1.0d0)*(Q(i-1,j,4)/rho - 0.5d0*(u**2.0d0 + v**2.0d0))
                    p   = rho*t
                    cL = sqrt(GAMMA*p/rho)
                    VnL = u*sin(beta) - v*cos(beta)
                    !if(j<2)then
                    !    write(*,*)i,j,VnL,cL
                    !endif

                    rho = Q(i+1,j,1)
                    u   = Q(i+1,j,2)/Q(i+1,j,1)
                    v   = Q(i+1,j,3)/Q(i+1,j,1)
                    t   = (GAMMA-1.0d0)*(Q(i+1,j,4)/rho - 0.5d0*(u**2 + v**2))
                    p   = rho*t
                    cR = sqrt(GAMMA*p/rho)
                    VnR = u*sin(beta) - v*cos(beta)
         
                    if((VnL-cL>0 .and. VnR-cR<0) .or. (VnL+cL>0 .and. VnR+cR<0))then
                        !case3-2
                        !write(*,*)i,j
                        !is_SF_xi(i,j+1) = 1
                        !is_SF_xi(i-1,j+1) = 1
                        !is_SF_xi(i,j-1) = 1
                        !is_SF_xi(i-1,j-1) = 1
                        is_SF_xi(i,j) = 1
                        is_SF_xi(i-1,j) = 1
                        !is_SF_xi(i+1,j) = 1
                        !is_SF_xi(i-2,j) = 1
                        is_SF_eta(i,j) = 1
                        is_SF_eta(i,j-1) = 1                        
                        !is_SF_eta(i-1,j) = 1
                        !is_SF_eta(i-1,j-1) = 1
                        !is_SF_eta(i+1,j) = 1
                        !is_SF_eta(i+1,j-1) = 1
                        !is_SF_eta(i,j-2) = 1
                        !is_SF_eta(i,j+1) = 1
                        !write(*,*)i,j
                    endif 
            
            !!!!! detection direction Upside and Downside !!!!!
                    VnL = Q(i,j-2,2)/Q(i,j-2,1)*sin(beta) - Q(i,j-2,3)/Q(i,j-2,1)*cos(beta)
                    VnR = Q(i,j+2,2)/Q(i,j+2,1)*sin(beta) - Q(i,j+2,3)/Q(i,j+2,1)*cos(beta)
                    p   = (GAMMA-1.0d0)*(Q(i,j-2,4) - 0.5d0*(Q(i,j-2,2)**2.0d0 + Q(i,j-2,3)**2.0d0)/Q(i,j-2,1))
                    cL = sqrt(GAMMA*p/Q(i,j-2,1))
                    p   = (GAMMA-1.0d0)*(Q(i,j+2,4) - 0.5d0*(Q(i,j+2,2)**2.0d0 + Q(i,j+2,3)**2.0d0)/Q(i,j+2,1))
                    cR = sqrt(GAMMA*p/Q(i,j+2,1))
                    f1 = sign(1.0d0,VnR-cR) * max(VnL/cL,1.0d0)
                    f2 = sign(1.0d0,VnL-cL) * max(VnR/cR,1.0d0)
                    f = min(0.0d0,f1*f2)
                    g = (Q(i,j+2,1) - 2.0d0*Q(i,j+2,1) + Q(i,j+2,1))*&
                        (Q(i,j,1) - 2.0d0*Q(i,j-2,1) + Q(i,j-2,1))
                    if(abs(f*g) > 0.0d0)then
                        is_SF_eta(i,j) = 1
                        is_SF_eta(i,j-1) = 1
                        is_SF_xi(i,j) = 1
                        is_SF_xi(i-1,j) = 1
                        !write(*,*)i,j
                    endif

                    rho = Q(i,j-1,1)
                    u   = Q(i,j-1,2)/Q(i,j-1,1)
                    v   = Q(i,j-1,3)/Q(i,j-1,1)
                    t   = (GAMMA-1.0d0)*(Q(i,j-1,4)/rho - 0.5d0*(u**2.0d0 + v**2.0d0))
                    p   = rho*t
                    cL = sqrt(GAMMA*p/rho)
                    VnL = u*sin(beta) - v*cos(beta)

                    rho = Q(i,j+1,1)
                    u   = Q(i,j+1,2)/Q(i,j+1,1)
                    v   = Q(i,j+1,3)/Q(i,j+1,1)
                    t   = (GAMMA-1.0d0)*(Q(i,j+1,4)/rho - 0.5d0*(u**2 + v**2))
                    p   = rho*t
                    cR = sqrt(GAMMA*p/rho)
                    VnR = u*sin(beta) - v*cos(beta)
         
                    if((VnL-cL>0 .and. VnR-cR<0) .or. (VnL+cL>0 .and. VnR+cR<0))then
                        !write(*,*)i,j
                        !case3-2
                        !is_SF_xi(i,j+1) = 1
                        !is_SF_xi(i-1,j+1) = 1
                        !is_SF_xi(i,j-1) = 1
                        !is_SF_xi(i-1,j-1) = 1
                        is_SF_xi(i,j) = 1
                        is_SF_xi(i-1,j) = 1
                        !is_SF_xi(i+1,j) = 1
                        !is_SF_xi(i-2,j) = 1
                        is_SF_eta(i,j) = 1
                        is_SF_eta(i,j-1) = 1                        
                        !is_SF_eta(i-1,j) = 1
                        !is_SF_eta(i-1,j-1) = 1
                        !is_SF_eta(i+1,j) = 1
                        !is_SF_eta(i+1,j-1) = 1
                        !is_SF_eta(i,j-2) = 1
                        !is_SF_eta(i,j+1) = 1
                    endif
                enddo                           
            enddo
        enddo
    end subroutine detectSonicPoint

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
                Cxi = 1.0d0/S(i,j)*(abs(U1) + c*sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0))
                Ceta = 1.0d0/S(i,j)*(abs(U2) + c*sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0))
                !write(*,*)sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0),m(i,j,1),m(i,j,2)

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

end module subprog

program main
    use subprog
    use AUSMDVFlux
    use SD_SLAUFlux
    !use SLAU2Flux
    use plot
    implicit none
    integer :: Nx,Ny,i=0,j=0,k=0,time=0,GridNum=11,foo=12,index=1,Qbin=100
    integer,allocatable :: is_SF_xi(:,:),is_SF_eta(:,:)
    character filenameraw*128,resfile*64
    real(8),allocatable :: Q(:,:,:),Qn(:,:,:),Qast(:,:,:),E(:,:,:),F(:,:,:),S(:,:),Ev(:,:,:),Fv(:,:,:)
    real(8),allocatable :: m(:,:,:),n(:,:,:)
    real(8) dt,qmax
    real(8) :: elapsedTime=0.0d0,Q1old=0.0d0,res=0.0d0,ressum=0.0d0,piyo=0.0d0


    open(GridNum,file = 'GridNum.txt')
    read(GridNum,*)Nx,Ny
    close(GridNum)

    allocate(Q(-1:Nx+2,-1:Ny+2,4),Qn(-1:Nx+2,-1:Ny+2,4),Qast(-1:Nx+2,-1:Ny+2,4),E(-1:Nx+2,-1:Ny+2,4),F(-1:Nx+2,-1:Ny+2,4),&
                S(-2:Nx+3,-2:Ny+3),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
                Ev(-1:Nx+2,-1:Ny+2,4),Fv(-1:Nx+2,-1:Ny+2,4),is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2))

    call setInitCond(Q,Nx,Ny)
    call setMetrics(m,n,S,Nx,Ny)
    
    write(resfile,'("residual error.dat")')
    open(foo,file = resfile)


    do !time=0,50        
        !call setbdycond(Q,m,n,Nx,Ny)
        call setbdycondCyl(Q,m,n,Nx,Ny)   
        call detectSonicPoint(Q,Nx,Ny,is_SF_xi,is_SF_eta)

        call evalFlux(Q,E,F,m,n,Nx,Ny,is_SF_xi,is_SF_eta)
        !do j=0,Ny
        !    do i=0,Nx
        !        if(is_SF_xi(i,j) == 1)then
        !            call SD_SLAU_Xflux(Q,E,m,Nx,Ny,i,j)
        !        else
        !            call AusmDV_XFlux(Q,E,m,Nx,Ny,i,j)
        !            !call SD_SLAU_Xflux(Q,E,m,Nx,Ny,i,j)
        !        endif
        !        
        !        
        !        if(is_SF_eta(i,j) == 1)then
        !            call SD_SLAU_YFlux(Q,F,n,Nx,Ny,i,j)
        !        else
        !            call AusmDV_YFlux(Q,F,n,Nx,Ny,i,j)
        !            !call SD_SLAU_YFlux(Q,F,n,Nx,Ny,i,j)
        !        endif
        !
        !        !call Xviscflux(Q,Ev,m,n,S,Nx,Ny)
        !        !call YviscFlux(Q,Fv,m,n,S,Nx,Ny)
        !    enddo
        !enddo


        do j=1,Ny
            do i=1,Nx !space
                do k=1,4
                    Qn(i,j,k) = Q(i,j,k)
                    Q(i,j,k) = Q(i,j,k) - 0.50d0*(dt/S(i,j))*(E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k))!&          !inertia
                                                        !- (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re) !viscous
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
        call setbdycondCyl(Qast,m,n,Nx,Ny)

        call evalFlux(Qast,E,F,m,n,Nx,Ny,is_SF_xi,is_SF_eta)


        !do j=0,Ny
        !    do i=0,Nx
        !        if(is_SF_xi(i,j) == 1)then
        !            !call AusmDV_XFlux_dum(Q,E,n,Nx,Ny,i,j)
        !            call SD_SLAU_Xflux(Qast,E,m,Nx,Ny,i,j)
        !        else
        !            call AusmDV_XFlux(Qast,E,m,Nx,Ny,i,j)
        !            !call SD_SLAU_XFlux(Qast,E,m,Nx,Ny,i,j)
        !        endif
        !        
        !        
        !        if(is_SF_eta(i,j) == 1)then
        !            !call AusmDV_YFlux_dum(Q,F,n,Nx,Ny,i,j)
        !            call SD_SLAU_YFlux(Qast,F,n,Nx,Ny,i,j)
        !        else
        !            call AusmDV_YFlux(Qast,F,n,Nx,Ny,i,j)
        !            !call SD_SLAU_YFlux(Qast,F,n,Nx,Ny,i,j)
        !        endif
        !        !call Xviscflux(Q,Ev,m,n,S,Nx,Ny)
        !        !call YviscFlux(Q,Fv,m,n,S,Nx,Ny)
        !    enddo
        !enddo

        qmax = calcSpeed(Q,Nx,Ny,m,n,S)
        dt =CFL/qmax
        
        do j=1,Ny
            do i=1,Nx !space
                Q1old = Qn(i,j,1)
                do k=1,4
                    Q(i,j,k) = Qn(i,j,k) - (dt/S(i,j))*(E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k))!& 
                                                        !- (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re)
                enddo

                res = abs(Q(i,j,1) - Q1old)
                ressum = ressum + res
            enddo
        enddo
        ressum = ressum/(Nx*Ny)
        elapsedTime = elapsedTime + dt
        
        if(modulo(time,50) == 0)then
            write(foo,*)elapsedTime,ressum
        !    call print_profile(Q,Nx,Ny,index)
        !    index = index + 1
            !piyo = piyo + 0.006d0
        endif
        time = time + 1

        write(*,*)ressum,elapsedTime
        
        if(isNaN(ressum))then
            write(*,*)'NaN'
            exit
        endif

        if(3.0d0 <= elapsedTime .or. ressum <= 0.0000000001d0)then!(ressum <= 0.0000000001d0 .and. 500 <time)then
            exit
        endif

        ressum = 0.0d0
        
        do j=0,Ny
            do i=0,Nx
                is_SF_xi(i,j) = 0
                is_SF_eta(i,j) = 0
            enddo
        enddo         
    enddo

    close(foo)


    call output2file(Q,Nx,Ny,index)

    deallocate(Q,E,F,S,m,n,Ev,Fv)

end program main