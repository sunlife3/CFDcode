program cylinderGrid
    !generate 2Dgrid of cylinder
    implicit none
    integer :: Nx,Ny,i,j,fi=11,fo=12
    real(8),allocatable ::  X(:,:),Y(:,:)
    character filename*128
    real(8) :: theta=0.0d0,PI=acos(-1.0d0),R=1.0d0,degRange=180.0d0,dtheta,dx,dR
    real(8),parameter :: Ratio = 1.0092d0

    open(fi,file = 'GridNum.txt')
    read(fi,*)Nx,Ny
    close(fi)
    write(*,*)Nx,Ny
    
    allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4))

    dtheta = degRange/Ny
    theta  = -degRange/2.0d0
    !R = R*comRatio**(Nx - 80 + 4)!R*comRatio**(Nx+4)
    dR = 3.0d0/Nx!todo...衝撃波離脱距離を格子数で割る

    do i=-4,Nx+3
        do j=1,Ny
            X(Nx-i,j) = -R*cos(theta*PI/180.0d0)
            Y(Nx-i,j) = R*sin(theta*PI/180.0d0)
            theta = theta + dtheta
            !write(*,*)theta
        enddo
        theta  = -degRange/2.0d0
       ! if(abs(X(Nx,Ny/2) - X(Nx-i,Ny/2)) < 0.8d0)then
       !     R = R + dR
       ! else
       !     R = R*Ratio
       ! endif
       	
       !R = R + dR ! 等間隔
       R = R*Ratio ! 等比
        
    enddo

    dx = R*sin(0.5d0*dtheta*PI/180.0d0)
    do i=-3,Nx+4
        do j=Ny+1,Ny+4
            X(i,j) = X(i,Ny) + (j-Ny)*dx
            Y(i,j) = Y(i,Ny)
        enddo
    enddo

    do i=-3,Nx+4
        do j=-3,0
            X(i,j) = X(i,1) + abs(j-1)*dx
            Y(i,j) = Y(i,1)
        enddo
    enddo


    write(filename,'("MESH_cylinderM20_eqsp(100100).txt")')
    open(fo,file = filename)
    do j=-3,Ny+4
        do i=-3,Nx+4
            write(fo,*)X(i,j),',',Y(i,j)
        enddo
    enddo
    close(fo)
    !write(*,*)X(100,1)-X(0,1)

    deallocate(X,Y)
end program cylinderGrid