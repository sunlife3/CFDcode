program rampGrid5
    !generate 2Dgrid of 5 deg ramp
    implicit none
    integer :: Nx,Ny,i,j,fi=11,fo=12,is_change_grad=1
    real(8),allocatable ::  X(:,:),Y(:,:)
    character filename*128
    real(8) :: dx=0,dy=0,theta=0,PI=3.141592,restheta=0.0,Ythreshold=0.0d0,Minf,beta

    open(fi,file = 'GridNum.txt')
    read(fi,*)Nx,Ny
    close(fi)
    write(*,*)Nx,Ny
    dx = 1.0d0/Nx
    theta = 5.0d0*PI/180.0d0
    beta  = 34.3d0
    Ythreshold = 0.80d0*tan(beta*PI/180.0d0)
    !M=2.0  beta = 34.3
    !M=5.0  beta = 15.1
    !M=10.0 beta = 9.5
    !M=15.0 beta = 7.9
    dy = Ythreshold/80.0d0
    
    
    allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4))
    Y(-3,-3) = - 0.50d0*3.0d0/100.0d0 

    do j=-3,Ny+4
        X(-3,j) = -0.5d0 - 0.50d0*3.0d0/100.0d0 
    enddo

    do i=-2,Nx+4
        if(i<20)then
            Y(i,-3) = Y(i-1,-3)
        else
            Y(i,-3) = Y(i-1,-3) + dx*tan(theta)
        endif
        X(i,-3) = X(i-1,-3) + dx
    enddo

    do j=-2,Ny+4
        do i=-2,Nx+4
            Y(-3,j) = Y(-3,j-1) + dy
            if(i < 20)then
                !dx = dx * 1.0d0/(4.0d0**0.05d0)
                !dx = 0.5/20 から dx = 0.5/80　へ滑らかに縮まるような公比が4^-0.05
                !0.5/20 * r^20 = 0.5/80
                !仮想セル分のマージンとって4^-0.04にしてみた
  
            else
                !dx = 0.50d0/80.0d0
            endif
            
            X(i,j) = X(i-1,j) + dx
            !write(*,*)i,j,dx,X(i,j)-X(i-1,j)            
            
            if(i < 20)then
                Y(i,j) = Y(i,j-1) + dy
            else
                Y(i,j) = Y(i-1,j) + dx*tan(theta)
            endif

        enddo

        if(1<=j)then
            if(is_change_grad == 1)then
                theta = theta - 0.1d0*PI/180
                restheta = theta - 0.0d0
                if(restheta <= 0.0001d0)then
                    is_change_grad = 0
                    theta = 0.0d0
                endif
            endif
        endif
        if(Ythreshold*1.20d0 <= Y(Nx+3,j))then
            dy = 1.1d0*dy
        endif
    enddo
    write(*,*)Y(Nx,Ny),Ythreshold

    write(filename,'("MESH_ramp",f3.1,"deg_MachWaveAngle",f4.1,".txt")')theta*180d0/PI,beta
    open(fo,file = filename)
    do j=-3,Ny+4
        do i=-3,Nx+4
            write(fo,*)X(i,j),',',Y(i,j)
        enddo
    enddo
    close(fo)
    !write(*,*)X(100,1)-X(0,1)

    deallocate(X,Y)
end program rampGrid5