    program hogehoge
        implicit none
        !integer,intent(in) :: Nx,Ny
        integer :: fo=22,i=0,j=0,GridNum=20,Grid=21,Qbin=19,index=1,is_output=0,Nx,Ny
        character filename*128
        !real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4)
        real(8),allocatable :: X(:,:),Y(:,:),Xcenter(:,:),Ycenter(:,:),Q(:,:,:)
        real(8) :: rho,u,v,p,T,eta,Minf,Tinf,pinf,einf
        real(8),parameter :: GAMMA=1.4d0,Re=5000.0d0

        open(GridNum,file = 'GridNum.txt')
        read(GridNum,*)Nx,Ny
        close(GridNum)


        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),&
                Xcenter(-2:Nx+4,-2:Ny+4),Ycenter(-2:Nx+4,-2:Ny+4),Q(-1:Nx+2,-1:Ny+2,4))

        Minf = 2.0d0
        pinf = 1.0d0/1.4d0
        einf = pinf/(1.40d0-1.0d0) + 0.5d0*(2.0d0**2.0d0) 
        Tinf = (GAMMA-1.0d0)*(einf - 0.5d0*(2.0d0**2.0d0))

        open(Qbin,file = 'Qbin001.dat')
        do j=-1,Ny+2
            do i=-1,Nx+2
                read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            enddo
        enddo
        close(Qbin)

        open(Grid,file = 'Gridtube.txt')
        do j=-3,Ny+4
            do i=-3,Nx+4
                read(Grid,*)X(i,j),Y(i,j)
            enddo
        enddo
        close(Grid)

        do j=-2,Ny+4
            do i=-2,Nx+4
                Xcenter(i,j)=0.25d0*(X(i-1,j-1)+X(i-1,j)+X(i,j-1)+X(i,j))           !cell center
                Ycenter(i,j)=0.25d0*(Y(i-1,j-1)+Y(i-1,j)+Y(i,j-1)+Y(i,j))
            enddo
        enddo


        do i=1,Nx
            if((X(i-1,0)-X(10,0)) < 1.0d0 .and. 1.0d0 <= (X(i,0)-X(10,0)))then
                is_output = 1
                
            else if((X(i-1,0)-X(10,0)) < 2.0d0 .and. 2.0d0 <= (X(i,0)-X(10,0)))then
                is_output = 1
                
            else if((X(i-1,0)-X(10,0)) < 2.5d0 .and. 2.5d0 <= (X(i,0)-X(10,0)))then
                is_output = 1
                
            else if((X(i-1,0)-X(10,0)) < 3.0d0 .and. 3.0d0 <= (X(i,0)-X(10,0)))then
                is_output = 1
                
            else if((X(i-1,0)-X(10,0)) < 4.0d0 .and. 4.0d0 <= (X(i,0)-X(10,0)))then
                is_output = 1
                
            else
                is_output = 0
            endif

            if(is_output == 1)then
                write(filename,'("BoundaryLayer",i4.4,".dat")')index
                open(fo,file = filename,status='replace')
                write(fo,*)'eta',' ','U/Minf',' ','T/Tinf'
                do j=1,Ny
                    eta = (Ycenter(i,j+1)-Y(i,1))*sqrt(Re*Minf/(Xcenter(i+1,j)-X(10,j)))
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    T   = (GAMMA - 1.0d0)*(Q(i,j,4)/Q(i,j,1) - 0.5d0*(u**2.0d0 + v**2.0d0))

                    write(fo,*)eta,u/Minf,T/Tinf
                    write(*,*)i,j,eta,T/Tinf
                    if(22.4d0 <= eta)then
                        exit
                    endif
                enddo
                close(fo)
                index = index + 1
                fo = fo +1
            endif
        enddo
        
    end program hogehoge