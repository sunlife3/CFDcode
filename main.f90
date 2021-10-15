module param
    implicit none
    real(8) :: GAMMA = 1.4d0,Re = 1000.0d0,CFL = 0.15d0
contains
    function minmod(x,y) result(fluxLimiter)
        real(8),intent(in) :: x,y
        real(8) fluxLimiter
        fluxLimiter = sign(1.0d0,x)*max(0.0d0,min(abs(x),y*sign(1.0d0,x)))
    end function minmod
end module param

program main
    use subprog
    use plot

    implicit none
    integer :: Nx,Ny,i=0,j=0,k=0,time=0,GridNum=11,index=1,Qbin=100,Grid=12,argidx=1
    integer,allocatable :: is_SF_xi(:,:),is_SF_eta(:,:)
    character filenameraw*128,resfile*64,meshfile*64,filename*128,argc*64
    real(8),allocatable :: Q(:,:,:),Qn(:,:,:),Qast(:,:,:),E(:,:,:),F(:,:,:),E1(:,:,:),F1(:,:,:),E2(:,:,:),F2(:,:,:),&
                           S(:,:),Ev(:,:,:),Fv(:,:,:)
    real(8),allocatable :: m(:,:,:),n(:,:,:),X(:,:),Y(:,:)!,is_SF_xi(:,:),is_SF_eta(:,:)
    real(8) dt,qmax
    real(8) :: elapsedTime=0.0d0,Q1old=0.0d0,res=0.0d0,ressum=0.0d0,piyo=0.0d0

    call getarg(argidx,argc)
    write(*,*)argc
    !!!!!!!!!!!!!! Grid Read !!!!!!!!!!!!!!!!!!!!!!
    if(argc == "SNS")then
        !meshfile = 'MESH_tube(5025).txt'
        meshfile = 'MESH_Tube_Buff(100025).txt'
    else if(argc == "Bow")then
        !meshfile = 'MESH_cylinderM2(100100).txt'
        meshfile = 'MESH_cylinder(100100).txt'
    else if(argc == "ramp")then
        meshfile = 'MESH_rampGridM15(60)(100100).txt'
        !meshfile = 'MESH_ductGridM15(100100).txt'
    else if(argc == "SWBLI")then
        meshfile = 'MESH_SWBLI(15000750).txt'
        !meshfile = 'MESH_ductGridM15(100100).txt'
    else if(argc == "visplate")then
        meshfile = 'MESH_GridVisPlate.txt'
        !meshfile = 'MESH_ductGridM15(100100).txt'
    else
        !open(Grid,file = 'MESH_TaperGridM15(100100).txt')
        !open(Grid,file = 'MESH_GridVisPlate.txt')
        !open(Grid,file = 'MESH_OddEven(35019).txt')
        meshfile = 'MESH_SOTube(500020).txt'
        !meshfile = 'MESH_LaxTube(200020).txt'
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open(Grid,file = meshfile)
    read(Grid,*)Nx,Ny
    write(*,*)Nx,Ny

    allocate(Q(-2:Nx+2,-2:Ny+2,4),Qn(-2:Nx+2,-2:Ny+2,4),Qast(-2:Nx+2,-2:Ny+2,4),E(-1:Nx+2,-1:Ny+2,4),F(-1:Nx+2,-1:Ny+2,4),&
                E1(-1:Nx+2,-1:Ny+2,4),F1(-1:Nx+2,-1:Ny+2,4),E2(-1:Nx+2,-1:Ny+2,4),F2(-1:Nx+2,-1:Ny+2,4),&
                S(-2:Nx+3,-2:Ny+3),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),&
                Ev(-1:Nx+2,-1:Ny+2,4),Fv(-1:Nx+2,-1:Ny+2,4),is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2))

    do j=-3,Ny+4
        do i=-3,Nx+4
            read(Grid,*)X(i,j),Y(i,j)
        enddo
    enddo
    close(Grid)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call setInitCond(Q,Nx,Ny,X,Y,argc)
    !call output2file(Q,Nx,Ny,index,meshfile)
    !index = index + 1
    call setMetrics(m,n,S,Nx,Ny,X,Y)
    

    do !time=0,0
        
        call setbdycond(Q,m,n,Nx,Ny,is_SF_xi,argc)
        !call detectSonicPoint(Q,Nx,Ny,is_SF_xi,is_SF_eta)
        !call ShockDetection(Q,Nx,Ny,is_SF_xi,is_SF_eta,X,Y)
        call evalFlux(Q,E,F,E1,E2,F1,F2,m,n,Nx,Ny,is_SF_xi,is_SF_eta)
        !call Xviscflux(Q,Ev,m,n,S,Nx,Ny)
        !call YviscFlux(Q,Fv,m,n,S,Nx,Ny)
        
        qmax = calcSpeed(Q,Nx,Ny,m,n,S)
        dt =CFL/qmax

        do j=1,Ny
            do i=1,Nx !space
                do k=1,4
                    Qn(i,j,k) = Q(i,j,k)
                    Qast(i,j,k) = Qn(i,j,k) - 0.50d0*(dt/S(i,j))*((E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k)) &      
                                                        - (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re) 
                enddo
            enddo
        enddo

        call setbdycond(Qast,m,n,Nx,Ny,is_SF_xi,argc)
        call evalFlux(Qast,E,F,E1,E2,F1,F2,m,n,Nx,Ny,is_SF_xi,is_SF_eta)
        !call Xviscflux(Q,Ev,m,n,S,Nx,Ny)
        !call YviscFlux(Q,Fv,m,n,S,Nx,Ny)
        
        do j=1,Ny
            do i=1,Nx !space
                !Q1old = Qn(i,j,1)
                do k=1,4
                    Q(i,j,k) = Qn(i,j,k) - (dt/S(i,j))*((E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k)) & 
                                                        - (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re)
                enddo
                res = abs(Q(i,j,1) - Q1old)
                ressum = ressum + res
            enddo
        enddo
        ressum = ressum/(Nx*Ny)
        elapsedTime = elapsedTime + dt
        
        if(0.1d0*index < elapsedtime)then
            !call print_profile(Q,Nx,Ny,index)
            call output2file(Q,Nx,Ny,index,meshfile)
            index = index + 1
            !piyo = piyo + 0.006d0
        endif

        write(*,*)elapsedTime,',',ressum
        
        if(isNaN(ressum))then
            write(*,*)'NaN'
            do j=1,Ny
                do i=1,Nx
                    do k=1,4
                        if(isNaN(Q(i,j,k)))then
                            !write(*,*)i,j,k,Q(i,j,k),Ev(i,j,k),Fv(i,j,k)
                        endif
                    enddo
                enddo
            enddo
            exit
        endif

        if(1.8d0 <= elapsedTime)then
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

    call output2file(Q,Nx,Ny,index,meshfile)
    !call pltdist(Q,Nx,Ny,meshfile)

    deallocate(Q,E,F,S,m,n,Ev,Fv,X,Y)

end program main
