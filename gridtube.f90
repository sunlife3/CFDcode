program tubeGrid
    !generate 2Dgrid of flat board
    implicit none
    integer :: Nx,Ny,i,j,fi=11,fo=12
    integer,parameter :: idxLeadEdge=10
    real(8),allocatable ::  X(:,:),Y(:,:)
    real(8) :: thickBL,dx,dy,intervalCoffY=1.0d0,intervalCoffX=5.0d0
    real(8) ,parameter :: XL = 5.0d0,Re = 5000.0d0,eta_d = 8.0d0,M_inf = 2.0d0

    open(fi,file = 'GridNum.txt')
    read(fi,*)Nx,Ny
    close(fi)
    write(*,*)Nx,Ny
    thickBL = eta_d*sqrt(XL/(Re*M_inf))
    !write(*,*)thickBL
    dx = XL/100.0d0
    dy = thickBL/60.0d0
        
    allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4))
    
    do i=-3,Nx+4
        Y(i,-3) = 0.0d0
    enddo
    
    do j=-3,Ny+4
        X(-3,j) = 0.0d0
    enddo

    do i=-2,Nx+4
        j=-3
        if(i < idxLeadEdge)then
            X(i,j) = X(i-1,j) + intervalCoffX*dx
        else
            X(i,j) = X(i-1,j) + dx
        endif
        
        if(1 <= i .and. i < idxLeadEdge-1)then
            intervalCoffX = intervalCoffX - 5.0d0/10.0d0
        else if(i == idxLeadEdge-1)then
            intervalCoffX = 5.0d0
        endif
    enddo


    do j=-2,Ny+4
        do i=-2,Nx+4
            if(Y(Nx+4,j-1) < thickBL+0.01d0)then
                Y(i,j) = Y(i,j-1) + dy
                Y(-3,j)= Y(-3,j-1)+ dy
            else
                Y(i,j) = Y(i,j-1) + intervalCoffY*dy
                Y(-3,j)= Y(-3,j-1)+ intervalCoffY*dy
            endif

            if(i < idxLeadEdge)then
                X(i,j) = X(i-1,j) + intervalCoffX*dx
            else
                X(i,j) = X(i-1,j) + dx
            endif

            if(1 <= i .and. i < idxLeadEdge-1)then
                intervalCoffX = intervalCoffX - 5.0d0/10.0d0
            else if(i == idxLeadEdge-1)then
                intervalCoffX = 5.0d0
            endif
        !if(j<=3)then
        !write(*,*)j,Y(Nx+4,j-1)
        !endif
        enddo
        if(thickBL+0.01d0 <= Y(Nx+4,j-1))then
            intervalCoffY = intervalCoffY + 0.6d0
        endif

    enddo
    
    open(fo,file = 'Gridtube.txt')
    !open(fo,file = 'Gridtube.csv')
    do j=-3,Ny+4
        do i=-3,Nx+4
            write(fo,*)X(i,j),Y(i,j)
        enddo
    enddo
    close(fo)
    write(*,*)Y(110,120)

    deallocate(X,Y)
end program tubeGrid