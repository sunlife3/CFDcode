
program visualize 
        implicit none
        !integer,intent(in) :: Nx,Ny
        integer :: fo=20,i=0,j=0,k=0,l=0,cnt,GridNum=18,Grid=19,Qbin=100,index=1,Nx,Ny,iold=0,jold=0
        character filename*128,meshfile*64,generalcoord*32,pltfile*64
        !real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4)
        real(8),allocatable :: X(:,:),Y(:,:),Q(:,:,:),m(:,:,:),n(:,:,:),Flag(:,:),prf(:,:)
        real(8) :: rho1,u1,v1,p1,t1,rho2,u2,v2,p2,t2,xcenter,ycenter,ax,ay,bx,by,cL1,cL2,cR1,cR2,VnL,VnR,VtL,VtR,&
                    betadeg,beta,betadeg2,beta2,kaku(2),theta,theta1,theta2,M1,hoge,fuga,tantheta,uave1,vave1,uave2,vave2,&
                    rho,u,v,p,t,mx,my,gradP,Vn(4),Vabs1,Vabs2,Entropy
        real(8),parameter :: GAMMA=1.4d0,PI=3.14159265359d0

        !open(GridNum,file = 'GridNum.txt')
        !read(GridNum,*)Nx,Ny
        !close(GridNum)
        !write(*,*)Nx,Ny

        
        
        do index=10,10
        !open(Grid,file = 'MESH_tube.txt')
        !write(filename,'("QbinSLAU2 SNS(M20,Buff,3timesLT)_",i3.3,".dat")')index
        !write(filename,'("QbinSF-SD-SLAU(meshM5).dat")')
        !write(filename,'("Qbin_SD_duct5deg(M10,100)_",i1.1,".dat")')index
        !write(filename,'("Qbin_SLAU2_ramp5deg(M2,100)_",i1.1,".dat")')index
        !write(filename,'("Qbin_AUSMDV3rd_duct5deg(M10,100,2step)_",i1.1,".dat")')index
        !write(filename,'("QbinHLL SNS(M2,b=25,e=10)_",i3.3,".dat")')index
        !write(filename,'("Qbin HLL Cylinder(M10,100)_",i3.3,".dat")')index
        !write(filename,'("Qbin_Exact_ramp5deg(M15,100)_",i1.1,".dat")')index
        !write(filename,'("Qbin_SD_DoubleMach(M5,200)_",i1.1,".dat")')index!
        write(filename,'("Qbin_SWBLI(1000)_",i2.2,".dat")')index
        open(Qbin,file = filename)
        read(Qbin,*)meshfile
        write(*,*)meshfile
        open(Grid,file = meshfile)
        read(Grid,*)Nx,Ny
        write(*,*)Nx,Ny
        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),Q(-1:Nx+2,-1:Ny+2,4),Flag(-1:Nx+2,-1:Ny+2),&
                    m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),prf(-1:Nx+2,-1:Ny+2))

        do j=-3,Ny+4
            do i=-3,Nx+4
                read(Grid,*)X(i,j),Y(i,j)
            enddo
        enddo
        close(Grid)

        !do j=-2,Ny+3
        !    do i=-2,Nx+3
        !        m(i,j,1) =  0.5d0*(0.5d0*(Y(i,j+1) - Y(i,j-1)) + 0.5d0*(Y(i+1,j+1) - Y(i+1,j-1))) 
        !        m(i,j,2) = -0.5d0*(0.5d0*(X(i,j+1) - X(i,j-1)) + 0.5d0*(X(i+1,j+1) - X(i+1,j-1)))
        !        n(i,j,1) = -0.5d0*(0.5d0*(Y(i+1,j) - Y(i-1,j)) + 0.5d0*(Y(i+1,j+1) - Y(i-1,j+1)))
        !        n(i,j,2) =  0.5d0*(0.5d0*(X(i+1,j) - X(i-1,j)) + 0.5d0*(X(i+1,j+1) - X(i-1,j+1)))
        !    enddo
        !enddo
        
        
            
            do j=-1,Ny+2
                do i=-1,Nx+2
                    read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                    !write(*,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                enddo
            enddo
            close(Qbin)
            !write(filename,'("BoundaryLayer",i3.3,".vtk")')index
            !open(fo,file = 'BoundaryLayer.vtk')
            !write(filename,'("duct5deg_HLLAUSM_M10(100)_",i3.3,".vtk")')index
            !write(filename,'("duct5deg_sd_M10(200100)_",i3.3,".vtk")')index
            !write(filename,'("ramp5deg_SD-SLAU_M2_",i3.3,".vtk")')index
            !write(filename,'("Bow Shock_HLL M20_100",i3.3,".vtk")')index
            !write(filename,'("Steady Normal Shock_SLAU2(M20,Buff,3LT)_",i3.3,".vtk")')index
            !write(filename,'("DoubleMach_60deg_sd_M10(200100)_",i3.3,".vtk")')index
            write(filename,'("Shock wave boundary layer interaction(1000)",i3.3,".vtk")')index
            open(fo,file = filename)
            write(fo,"('# vtk DataFile Version 3.0')")
            write(fo,"('ramp')")
            write(fo,"('ASCII')")
            write(fo,"('DATASET STRUCTURED_GRID')")
            write(fo,"('DIMENSIONS',3(1x,i4))") Nx,Ny,1
            write(fo,"('POINTS ',i9,' float')") Nx*Ny

            do j=1,Ny
                do i=1,Nx
                    xcenter  = 0.25d0 *(X(i,j) + X(i-1,j) + X(i,j-1) +X(i-1,j-1))
                    ycenter  = 0.25d0 *(Y(i,j) + Y(i-1,j) + Y(i,j-1) +Y(i-1,j-1))
                    !xcenter = m(i,j,2) * 
                    write(fo,'(3(f9.4,1x))')xcenter,ycenter,0.0d0
                enddo
            enddo

            write(fo,"('POINT_DATA ',i9)") Nx*Ny
            write(fo,"('SCALARS rho float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    write(fo,"(f10.6)")Q(i,j,1)
                    
                enddo
            enddo

            write(fo,"('SCALARS U float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    rho = Q(i,j,1)
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    p   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                    t   = p/rho
                    write(fo,"(f10.7)")u!*m(i,j,1)+&!/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0) + &
                                       !v*m(i,j,2)!&/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)&
                                        ! - sqrt(GAMMA * p/rho)!*sqrt((m(i,j,1)/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0))**2.0d0 + &
                                                                !     (m(i,j,2)/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0))**2.0d0))
                enddo
            enddo

            write(fo,"('SCALARS V float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    rho = Q(i,j,1)
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    p   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                    t   = p/rho
                    write(fo,"(f10.7)")v!u*n(i,j,1)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0) + &
            !                           v*n(i,j,2)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)!&
                                         !- (sqrt(GAMMA * p/rho)*sqrt((n(i,j,1)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0))**2.0d0 + &
                                         !                            (n(i,j,2)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0))**2.0d0))
                enddo
            enddo

            write(fo,"('SCALARS p float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    rho = Q(i,j,1)
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    p   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                    t   = p/rho
                    write(fo,"(f10.6)")p
                    !write(*,*)'C',i,j,sqrt(GAMMA*p/rho)
                enddo
            enddo

            write(fo,"('SCALARS f float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx

                    rho = Q(i-1,j,1)
                    u   = Q(i-1,j,2)/Q(i-1,j,1)
                    v   = Q(i-1,j,3)/Q(i-1,j,1)
                    p   = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))

                    rho1 = Q(i,j,1)
                    u1   = Q(i,j,2)/Q(i,j,1)
                    v1   = Q(i,j,3)/Q(i,j,1)
                    p1   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*rho1*(u1**2.0d0 + v1**2.0d0))
                    t1   = p1/rho1

                    rho2 = Q(i+1,j,1)
                    u2   = Q(i+1,j,2)/Q(i+1,j,1)
                    v2   = Q(i+1,j,3)/Q(i+1,j,1)
                    p2   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*rho2*(u2**2.0d0 + v2**2.0d0))
                    t2   = p2/rho2
                    
                    write(fo,"(f10.7)")min(p1/p2,p2/p1)**3.0d0
                    !write(fo,"(f10.7)")min(p2/p1,p1/p2)**3.0d0
                    !write(*,*)'C',i,j,sqrt(GAMMA*p/rho)
                enddo
            enddo

            write(fo,"('SCALARS Flag float')")
            write(fo,"('LOOKUP_TABLE default')")
            betadeg=0.0d0!asin(PI/3600.0d0)
            betadeg2=1.0d0
            do j=1,Ny
                do i=1,Nx
                    Flag(i,j) = 0.0d0
                enddo
            enddo
            do j=1,Ny
                do i=1,Nx
                    mx = m(i-2,j,1)!/sqrt(m(i-1,j,1)**2.0d0 + m(i-1,j,2)**2.0d0)
                    my = m(i-2,j,2)!/sqrt(m(i-1,j,1)**2.0d0 + m(i-1,j,2)**2.0d0)
                    rho1 = Q(i-1,j,1)
                    u1   = Q(i-1,j,2)/Q(i-1,j,1)!(sqrt(Q(i,j,1))*uave1 + sqrt(Q(i-1,j,1))*uave2)/(sqrt(Q(i,j,1))+sqrt(Q(i-1,j,1)))
                    v1   = Q(i-1,j,3)/Q(i-1,j,1)!(sqrt(Q(i,j,1))*vave1 + sqrt(Q(i-1,j,1))*vave2)/(sqrt(Q(i,j,1))+sqrt(Q(i-1,j,1)))
                    theta1 = atan(v1/u1)
                    p1   = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1))
                    t1   = p1/rho1
                    cL1 = sqrt(GAMMA*p1/rho1)!*sqrt(mx**2.0d0 + my**2.0d0)
                    VnL = u1*mx + v1*my

                    mx = m(i,j,1)/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)
                    my = m(i,j,2)/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)

                    rho2 = Q(i+1,j,1)
                    u2   = Q(i+1,j,2)/Q(i+1,j,1)!(sqrt(Q(i,j,1))*uave1 + sqrt(Q(i+1,j,1))*uave2)/(sqrt(Q(i,j,1))+sqrt(Q(i+1,j,1)))!(Q(i,j,2) + Q(i+1,j,2))/(Q(i,j,1) + Q(i+1,j,1))
                    v2   = Q(i+1,j,3)/Q(i+1,j,1)!(sqrt(Q(i,j,1))*vave1 + sqrt(Q(i+1,j,1))*vave2)/(sqrt(Q(i,j,1))+sqrt(Q(i+1,j,1)))!(Q(i,j,3) + Q(i+1,j,3))/(Q(i,j,1) + Q(i+1,j,1))
                    theta2 = atan(v2/u2)
                    p2   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                    t2   = p2/rho2
                    cR1 = sqrt(GAMMA*p2/rho2)
                    
                    !if((u1 - cL1 > 0 .and. u2 - cR1 < 0) .or. (u1 + cL1 > 0 .and. u2 + cR1 <0))then
                    !    Flag(i,j) = 1.0d0
                    !endif

                    prf(i,j) = (min(p1/p2,p2/p1))**3.0d0
                    !write(*,*)i,j,prf(i,j)
                    !Flag(i,j) = prf(i,j)
                    
                    if(prf(i,j) < 0.55d0)then
                        ! Shock Detection by velocity
                        do l=1,Nx 
                            p1   = (GAMMA-1.0d0)*(Q(l-1,j-2,4) - 0.5d0*(Q(l-1,j-2,2)**2.0d0 + Q(l-1,j-2,3)**2.0d0)/Q(l-1,j-2,1))
                            p2   = (GAMMA-1.0d0)*(Q(l+1,j-2,4) - 0.5d0*(Q(l+1,j-2,2)**2.0d0 + Q(l+1,j-2,3)**2.0d0)/Q(l+1,j-2,1))
                            prf(l,j-2) = (min(p1/p2,p2/p1))**3.0d0

                            if(prf(l,j-2) < 0.55d0)then
                                
                                theta = atan((Y(i,j) - Y(l,j-2))/(X(i,j) - X(l,j-2)))
                                if(theta < 0.0d0)then
                                    VnL = -u1*sin(theta) + v1*cos(theta)
                                    VnR = -u2*sin(theta) + v2*cos(theta)
                                else
                                    VnL = u1*sin(theta) - v1*cos(theta)
                                    VnR = u2*sin(theta) - v2*cos(theta)
                                endif
                                !if(13.0d0 < theta*180/PI .and. theta*180/PI < 17.0d0)then
                                !    write(*,*)i,j,theta*180/PI,VnL-cL1,VnR-cR1,"! 1"
                                !endif
                                if((VnL - cL1 > 0 .and. VnR - cR1 < 0) .or. (VnL + cL1 > 0 .and. VnR + cR1 <0))then
                                    !write(*,*)i,j,theta*180/PI,VnL-cL1,VnR-cR1,"! 1"
                                    if(sqrt(u1**2.0d0 + v1**2.0d0) < sqrt(u2**2.0d0 + v2**2.0d0))then
                                        !Flag(i,j) = 0.0d0
                                    else
                                        Flag(i,j) = 1.0d0                 
                                        
                                        !write(*,*)'xi',i,j,theta*180.0d0/PI,u1,v1,u2,v2,&
                                        !sqrt(u1**2.0d0 + v1**2.0d0),sqrt(u2**2.0d0 + v2**2.0d0)
                                    endif
                                endif
                                !write(*,*)i,j,theta*180.0d0/PI,Flag(i,j),VnL-cL1,VnR-cR1
                            endif
                        enddo
                    else
                        !Flag(i,j) = (min(p1/p2,p2/p1))**3.0d0
                        !Flag(i,j) = 1.0d0
                    endif

                    mx = n(i,j-1,1)/sqrt(n(i,j-1,1)**2.0d0 + n(i,j-1,2)**2.0d0)
                    my = n(i,j-1,2)/sqrt(n(i,j-1,1)**2.0d0 + n(i,j-1,2)**2.0d0)
                    rho1 = Q(i,j-1,1)
                    u1   = Q(i,j-1,2)/Q(i,j-1,1)!(sqrt(Q(i,j+1,1))*uave1 + sqrt(Q(i,j,1))*uave2)/(sqrt(Q(i,j+1,1))+sqrt(Q(i,j,1)))!(Q(i,j,2) + Q(i-1,j,2))/(Q(i,j,1) + Q(i-1,j,1))
                    v1   = Q(i,j-1,3)/Q(i,j-1,1)!(sqrt(Q(i,j+1,1))*vave1 + sqrt(Q(i,j,1))*vave2)/(sqrt(Q(i,j+1,1))+sqrt(Q(i,j,1)))!(Q(i,j,3) + Q(i-1,j,3))/(Q(i,j,1) + Q(i-1,j,1))
                    theta1 = atan(v1/u1)
                    p1   = (GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1))
                    t1   = p1/rho1
                    !cL2 = sqrt(GAMMA*p1/rho1)*sqrt(n(i,j-1,1)**2.0d0 + n(i,j-1,2)**2.0d0)!sqrt(GAMMA*p1/rho1)*sqrt(n(i,j-1,1)**2.0d0 + n(i,j-1,2)**2.0d0)
                    cL2 = sqrt(GAMMA*p1/rho1)!*sqrt(mx**2.0d0 + my**2.0d0)
                    M1 = 20.0d0  
                    !VnL = u1*n(i,j-1,1) + v1*n(i,j-1,2)
                    !VnL = u1*mx + v1*my

                    mx = n(i,j,1)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)
                    my = n(i,j,2)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)
                    rho2 = Q(i,j+1,1)
                    u2   = Q(i,j+1,2)/Q(i,j+1,1)!(sqrt(Q(i,j-1,1))*uave1 + sqrt(Q(i,j,1))*uave2)/(sqrt(Q(i,j-1,1))+sqrt(Q(i,j,1)))!(Q(i,j,2) + Q(i-1,j,2))/(Q(i,j,1) + Q(i-1,j,1))
                    v2   = Q(i,j+1,3)/Q(i,j+1,1)!(sqrt(Q(i,j-1,1))*vave1 + sqrt(Q(i,j,1))*vave2)/(sqrt(Q(i,j-1,1))+sqrt(Q(i,j,1)))!(Q(i,j,3) + Q(i-1,j,3))/(Q(i,j,1) + Q(i-1,j,1))
                    theta2 = atan(v2/u2)
                    p2   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
                    t2   = p2/rho2
                    cR2 = sqrt(GAMMA*p2/rho2)!*sqrt(mx**2.0d0 + my**2.0d0)
                    !VnR = u2*mx + v2*my
                    

                    prf(i,j) = (min(p1/p2,p2/p1))**3.0d0

                    if(prf(i,j) < 0.5d0)then

                        ! Shock Detection by compare Vn
                        do l=1,Ny
                            p1   = (GAMMA-1.0d0)*(Q(i-2,l-1,4) - 0.5d0*(Q(i-2,l-1,2)**2.0d0 + Q(i-2,l-1,3)**2.0d0)/Q(i-2,l-1,1))
                            p2   = (GAMMA-1.0d0)*(Q(i-2,l+1,4) - 0.5d0*(Q(i-2,l+1,2)**2.0d0 + Q(i-2,l+1,3)**2.0d0)/Q(i-2,l+1,1))
                            prf(i-2,l) = (min(p1/p2,p2/p1))**3.0d0
                            
                            if(prf(i-2,l) < 0.5d0)then
                                theta = atan((Y(i,j) - Y(i-2,l))/(X(i,j) - X(i-2,l)))
                                
                                if(90.0d0*PI/180.0d0 < theta .or. theta < -90.0d0*PI/180.0d0)then
                                endif
                                if(theta < 0.0d0)then
                                    VnL = -u1*sin(theta) + v1*cos(theta)
                                    VnR = -u2*sin(theta) + v2*cos(theta)
                                else
                                    VnL = u1*sin(theta) - v1*cos(theta)
                                    VnR = u2*sin(theta) - v2*cos(theta)
                                endif
                                
                                if(v1 < 0.0d0)then
                                    if((VnL - cL2 > 0 .and. VnR - cR2 < 0) .or. (VnL + cL2 > 0 .and. VnR + cR2 <0))then!Compresible sonic point
                                        if(sqrt(u1**2.0d0 + v1**2.0d0) < sqrt(u2**2.0d0 + v2**2.0d0))then!NOT accerelation
                                            !Flag(i,j) = 0.0d0
                                        else
                                            if(theta < 0.0d0)then 
                                            Flag(i,j) = 1.0d0
                                            endif
                                            
                                            !write(*,*)i,j,l,theta*180/PI,VnL,cL2,sqrt(u1**2.0d0 + v1**2.0d0),VnR,cR2,sqrt(u2**2.0d0 + v2**2.0d0)
                                            !write(*,*)i,j,theta*180/PI,VnL,u1,v1,u2,v2,VnL,VnR
                                        endif
                                    endif
                                else
                                    !if(13.0d0 < theta*180/PI .and. theta*180/PI < 17.0d0)then
                                    !    write(*,*)i,j,theta*180/PI,VnL-cL2,VnR-cR2,"! 3"
                                    !endif
                
                                    if((VnR - cR2 > 0 .and. VnL - cL2 < 0) .or. (VnR + cR2 > 0 .and. VnL + cL2 <0))then!Compresible sonic point
                                        
                                        if(sqrt(u1**2.0d0 + v1**2.0d0) > sqrt(u2**2.0d0 + v2**2.0d0))then!NOT accerelation
                                            !Flag(i,j) = 0.0d0
                                        else
                                            !if(theta > 0.0d0)then
                                            Flag(i,j) = 1.0d0
                                            !write(*,*)i,j,theta*180/PI,VnL,u1,v1,u2,v2,VnL,VnR
                                            !endif
                                            
                                            !write(*,*)i,j,l,theta*180/PI,VnR,cR2,sqrt(u2**2.0d0 + v2**2.0d0),VnL,cL2,sqrt(u1**2.0d0 + v1**2.0d0)
                                        endif
                                    endif
                                endif                        
                                !write(*,*)X(i,j),',',Y(i,j),',',prf(i,j),',',theta*180.0d0/PI
                            endif
                        enddo
                    else
                        
                    endif
                enddo
            enddo

            do j=1,Ny
                do i=1,Nx
                    write(fo,"(f10.7)")Flag(i,j)
                enddo
            enddo
            close(fo)
           
            deallocate(X,Y,Q,Flag,m,n,prf)
        enddo
        
            index = index + 1
end program visualize