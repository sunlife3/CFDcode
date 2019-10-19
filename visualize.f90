program visualize
        implicit none
        !integer,intent(in) :: Nx,Ny
        integer :: fo=20,i=0,j=0,GridNum=18,Grid=19,Qbin=100,index=1,Nx,Ny
        character filename*128
        !real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4)
        real(8),allocatable :: X(:,:),Y(:,:),Q(:,:,:)
        real(8) :: rho,u,v,p,t,xcenter,ycenter
        real(8),parameter :: GAMMA=1.4d0

        open(GridNum,file = 'GridNum.txt')
        read(GridNum,*)Nx,Ny
        close(GridNum)
        write(*,*)Nx,Ny

        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),Q(-1:Nx+2,-1:Ny+2,4))

        !do i = -1,Nx+2
        !    hoge = sqrt((Q(i,Ny-4,2)/Q(i,Ny-4,1))**2 + (Q(i,Ny-4,3)/Q(i,Ny-4,1))**2)
        !    write(*,*)i,Ny-4,hoge
        !enddo
        
        open(Grid,file = 'Gridtube.txt')
        do j=-3,Ny+4
            do i=-3,Nx+4
                read(Grid,*)X(i,j),Y(i,j)
            enddo
        enddo
        close(Grid)        
        
        do index=1,49
            write(filename,'("Qbin",i3.3,".dat")')index
            open(Qbin,file = filename)
            do j=-1,Ny+2
                do i=-1,Nx+2
                    read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                enddo
            enddo
            close(Qbin)
        
            write(filename,'("BoundaryLayer",i3.3,".vtk")')index
            !open(fo,file = 'BoundaryLayer.vtk')
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

            write(fo,"('SCALARS t float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    rho = Q(i,j,1)
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    t   = (GAMMA-1.0d0)*(Q(i,j,4)/rho - 0.5d0*(u**2 + v**2))
                    write(fo,"(f10.7)")t
                enddo
            enddo        
            close(fo)
        
        enddo
        
            deallocate(X,Y)
            index = index + 1
end program visualize