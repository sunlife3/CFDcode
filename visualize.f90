program visualize
        implicit none
        !integer,intent(in) :: Nx,Ny
        integer :: fo=20,i=0,j=0,GridNum=18,Grid=19,Qbin=100,index=1,Nx,Ny
        character filename*128
        !real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4)
        real(8),allocatable :: X(:,:),Y(:,:),Q(:,:,:),m(:,:,:),n(:,:,:),Flag(:,:)
        real(8) :: rho,u,v,p,t,xcenter,ycenter,ax,ay,bx,by,cL1,cL2,cR1,cR2,VnL,VnR,beta!UxiL,UetaL,UxiR,UetaR,
        real(8),parameter :: GAMMA=1.4d0,PI=3.14159265359d0

        open(GridNum,file = 'GridNum.txt')
        read(GridNum,*)Nx,Ny
        close(GridNum)
        write(*,*)Nx,Ny

        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),Q(-1:Nx+2,-1:Ny+2,4),&
                    m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),Flag(-1:Nx+2,-1:Ny+2))

        !do i = -1,Nx+2
        !    hoge = sqrt((Q(i,Ny-4,2)/Q(i,Ny-4,1))**2 + (Q(i,Ny-4,3)/Q(i,Ny-4,1))**2)
        !    write(*,*)i,Ny-4,hoge
        !enddo
        
        !open(Grid,file = 'Gridtube.txt')
        open(Grid,file = 'MESH_rampGrid.txt')
        !open(Grid,file = 'MESH_GridVisPlate.txt')
        
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
            enddo
        enddo
        
        do index=1,20
            write(filename,'("QbinAUSMDV",i3.3,".dat")')index
            open(Qbin,file = filename)
            do j=-1,Ny+2
                do i=-1,Nx+2
                    read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                enddo
            enddo
            close(Qbin)
        
            !write(filename,'("BoundaryLayer",i3.3,".vtk")')index
            !open(fo,file = 'BoundaryLayer.vtk')
            write(filename,'("Ramp5deg_AUSMDV_M2_",i3.3,".vtk")')index
            open(fo,file = filename)
            write(fo,"('# vtk DataFile Version 3.0')")
            write(fo,"('ramp')")
            write(fo,"('ASCII')")
            write(fo,"('DATASET STRUCTURED_GRID')")
            write(fo,"('DIMENSIONS',3(1x,i3))") Nx,Ny,1
            write(fo,"('POINTS ',i9,' float')") Nx*Ny

            do j=1,Ny
                do i=1,Nx
                    xcenter  = 0.25d0 *(X(i,j) + X(i-1,j) + X(i,j-1) +X(i-1,j-1))
                    ycenter  = 0.25d0 *(Y(i,j) + Y(i-1,j) + Y(i,j-1) +Y(i-1,j-1))
                    write(fo,'(3(f9.4,1x))')xcenter,ycenter,0.0d0
                enddo
            enddo

            write(fo,"('POINT_DATA',i9)")Nx*Ny
            write(fo,"('SCALARS rho float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    write(fo,"(f10.7)")Q(i,j,1)
                enddo
            enddo

            write(fo,"('SCALARS velocity float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    write(fo,"(f10.7)")u
                enddo
            enddo

            write(fo,"('SCALARS p float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    rho = Q(i,j,1)
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    t   = (GAMMA-1.0d0)*(Q(i,j,4)/rho - 0.5d0*(u**2 + v**2))
                    p   = rho*t
                    write(fo,"(f10.7)")p
                enddo
            enddo        

            write(fo,"('SCALARS Flag float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    beta = 9.5d0*PI/180d0
                    rho = Q(i-1,j,1)
                    u   = Q(i-1,j,2)/Q(i-1,j,1)
                    v   = Q(i-1,j,3)/Q(i-1,j,1)
                    t   = (GAMMA-1.0d0)*(Q(i-1,j,4)/rho - 0.5d0*(u**2.0d0 + v**2.0d0))
                    p   = rho*t
                    cL1 = sqrt(GAMMA*p/rho)
                    VnL = u*sin(beta) - v*cos(beta)
                    
                    rho = Q(i,j-1,1)
                    u   = Q(i,j-1,2)/Q(i,j-1,1)
                    v   = Q(i,j-1,3)/Q(i,j-1,1)
                    t   = (GAMMA-1.0d0)*(Q(i,j-1,4)/rho - 0.5d0*(u**2.0d0 + v**2.0d0))
                    p   = rho*t
                    cL2 = sqrt(GAMMA*p/rho)

                    !ax = m(i-1,j,1)
                    !ay = m(i-1,j,2)
                    !ax = ax/sqrt(m(i-1,j,1)**2 + m(i-1,j,2)**2)
                    !ay = ay/sqrt(m(i-1,j,1)**2 + m(i-1,j,2)**2)
                    !UxiL = (Q(i-1,j,2)*ax + Q(i-1,j,3)*ay)/Q(i-1,j,1)


                    !ax = n(i,j-1,1)
                    !ay = n(i,j-1,2)
                    !ax = ax/sqrt(n(i,j-1,1)**2 + n(i,j-1,2)**2)
                    !ay = ay/sqrt(n(i,j-1,1)**2 + n(i,j-1,2)**2)
                    !UetaL= (Q(i,j-1,2)*ax + Q(i,j-1,3)*ay)/Q(i,j-1,1)

                    !ax = m(i+1,j,1)
                    !ay = m(i+1,j,2)
                    !ax = ax/sqrt(m(i+1,j,1)**2 + m(i+1,j,2)**2)
                    !ay = ay/sqrt(m(i+1,j,1)**2 + m(i+1,j,2)**2)
                    !UxiR = (Q(i+1,j,2)*ax + Q(i+1,j,3)*ay)/Q(i+1,j,1)

                    !ax = n(i,j+1,1)
                    !ay = n(i,j+1,2)
                    !ax = ax/sqrt(n(i,j+1,1)**2 + n(i,j+1,2)**2)
                    !ay = ay/sqrt(n(i,j+1,1)**2 + n(i,j+1,2)**2)
                    !UetaR= (Q(i,j+1,2)*ax + Q(i,j+1,3)*ay)/Q(i,j+1,1)

                    rho = Q(i+2,j,1)
                    u   = Q(i+2,j,2)/Q(i+2,j,1)
                    v   = Q(i+2,j,3)/Q(i+2,j,1)
                    t   = (GAMMA-1.0d0)*(Q(i+2,j,4)/rho - 0.5d0*(u**2 + v**2))
                    p   = rho*t
                    cR1 = sqrt(GAMMA*p/rho)
                    VnR = u*sin(beta) - v*cos(beta)

                    rho = Q(i,j+2,1)
                    u   = Q(i,j+2,2)/Q(i,j+2,1)
                    v   = Q(i,j+2,3)/Q(i,j+2,1)
                    t   = (GAMMA-1.0d0)*(Q(i,j+2,4)/rho - 0.5d0*(u**2 + v**2))
                    p   = rho*t
                    cR2 = sqrt(GAMMA*p/rho)

                    !if(j==1)then
                    !    write(*,*)i,cL1,VnL,cR1,VnR
                    !endif
                    
                    if((VnL-cL1>0 .and. VnR-cR1<0) .or. (VnL+cL1>0 .and. VnR+cR1<0))then
                        Flag(i,j) = 1.0d0
                        !write(*,*)i,j
                    endif
                    !if((UetaL-cL2>0 .and. UxiR-cR2<0) .or. (UetaL+cL2>0 .and. UetaR+cR2<0))then
                    !    Flag(i,j) = 1.0d0
                    !endif

                    write(fo,"(f10.7)")Flag(i,j)
                enddo
            enddo            
            
            close(fo)

            do j=1,Ny
                do i=1,Nx
                    Flag(i,j) = 0.0d0
                enddo
            enddo            

        enddo
        
            deallocate(X,Y)
            index = index + 1
end program visualize