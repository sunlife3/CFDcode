
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

        open(GridNum,file = 'GridNum.txt')
        read(GridNum,*)Nx,Ny
        close(GridNum)
        write(*,*)Nx,Ny

        allocate(X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4),Q(-1:Nx+2,-1:Ny+2,4),Flag(-1:Nx+2,-1:Ny+2),&
                    m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),prf(-1:Nx+2,-1:Ny+2))
        
        
        !open(Grid,file = 'MESH_tube.txt')
        !write(filename,'("QbinSLAU2_SteadyNormalShock(1storder)",i3.3,".dat")')index
        !write(filename,'("QbinSF-SD-SLAU(meshM5).dat")')
        write(filename,'("Qbin_SD_duct5deg(M15,100)_",i1.1,".dat")')index
        !write(filename,'("Qbin_SD_ramp5deg(60)(M15,100)_",i1.1,".dat")')index
        !write(filename,'("QbinSD- Cylinder(M20,100)_",i3.3,".dat")')index
        open(Qbin,file = filename)
        read(Qbin,*)meshfile
        write(*,*)meshfile
        open(Grid,file = meshfile)

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
            enddo
        enddo
        
        !open(77,file = 'MESH_GeneralCoord(100100).txt')
        !do j=-3,Ny+4
        !    do i=-3,Nx+4
        !        read(77,*)X(i,j),Y(i,j)
        !    enddo
        !enddo
        !close(77)

        do index=1,1
            
            do j=-1,Ny+2
                do i=-1,Nx+2
                    read(Qbin,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                    !write(*,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
                enddo
            enddo
            close(Qbin)
        
            !write(filename,'("BoundaryLayer",i3.3,".vtk")')index
            !open(fo,file = 'BoundaryLayer.vtk')
            !write(filename,'("ramp5deg_sd_M15(5,60)(100)_",i3.3,".vtk")')index
            write(filename,'("duct5deg_sd_M15(200100)_",i3.3,".vtk")')index
            !write(filename,'("Bow Shock_SD-AUSMDV M20_100_valid",i3.3,".vtk")')index
            !write(filename,'("Steady Normal Shock_SLAU2(1storder)_",i3.3,".vtk")')index
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
                    !xcenter = m(i,j,2) * 
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

            !write(fo,"('SCALARS V float')")
            !write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    rho = Q(i,j,1)
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    p   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                    t   = p/rho
            !        write(fo,"(f10.7)")u*n(i,j,1)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0) + &
            !                           v*n(i,j,2)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)!&
                                         !- (sqrt(GAMMA * p/rho)*sqrt((n(i,j,1)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0))**2.0d0 + &
                                         !                            (n(i,j,2)/sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0))**2.0d0))
                enddo
            enddo

            write(fo,"('SCALARS C float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
                    rho = Q(i,j,1)
                    u   = Q(i,j,2)/Q(i,j,1)
                    v   = Q(i,j,3)/Q(i,j,1)
                    p   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*rho*(u**2.0d0 + v**2.0d0))
                    t   = p/rho
                    write(fo,"(f10.7)")sqrt(GAMMA*p/rho)
                    !write(*,*)'C',i,j,sqrt(GAMMA*p/rho)
                enddo
            enddo

            write(fo,"('SCALARS Entropy float')")
            write(fo,"('LOOKUP_TABLE default')")
            do j=1,Ny
                do i=1,Nx
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
                    
                    write(fo,"(f10.7)")log(t2/(rho2**(GAMMA-1.0d0))*(rho1**(GAMMA-1.0d0)/t1))
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
                    mx = m(i-1,j,1)!/sqrt(m(i-1,j,1)**2.0d0 + m(i-1,j,2)**2.0d0)
                    my = m(i-1,j,2)!/sqrt(m(i-1,j,1)**2.0d0 + m(i-1,j,2)**2.0d0)

                    rho1 = Q(i-1,j,1)
                    u1   = Q(i-1,j,2)/Q(i-1,j,1)!(sqrt(Q(i,j,1))*uave1 + sqrt(Q(i-1,j,1))*uave2)/(sqrt(Q(i,j,1))+sqrt(Q(i-1,j,1)))
                    v1   = Q(i-1,j,3)/Q(i-1,j,1)!(sqrt(Q(i,j,1))*vave1 + sqrt(Q(i-1,j,1))*vave2)/(sqrt(Q(i,j,1))+sqrt(Q(i-1,j,1)))
                    theta1 = atan(v1/u1)
                    p1   = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1))
                    t1   = p1/rho1
                    !cL1 = sqrt(GAMMA*p1/rho1)*sqrt(m(i-1,j,1)**2.0d0 + m(i-1,j,2)**2.0d0)
                    cL1 = sqrt(GAMMA*p1/rho1)!*sqrt(mx**2.0d0 + my**2.0d0)
                    M1 = 20.0d0  
                    !VnL = (u1*m(i-1,j,1) + v1*m(i-1,j,2))
                    VnL = u1*mx + v1*my
                    !Vabs1 = sqrt(u1**2.0d0 + v1**2.0d0)

                    mx = m(i,j,1)/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)
                    my = m(i,j,2)/sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)

                    rho2 = Q(i+1,j,1)
                    u2   = Q(i+1,j,2)/Q(i+1,j,1)!(sqrt(Q(i,j,1))*uave1 + sqrt(Q(i+1,j,1))*uave2)/(sqrt(Q(i,j,1))+sqrt(Q(i+1,j,1)))!(Q(i,j,2) + Q(i+1,j,2))/(Q(i,j,1) + Q(i+1,j,1))
                    v2   = Q(i+1,j,3)/Q(i+1,j,1)!(sqrt(Q(i,j,1))*vave1 + sqrt(Q(i+1,j,1))*vave2)/(sqrt(Q(i,j,1))+sqrt(Q(i+1,j,1)))!(Q(i,j,3) + Q(i+1,j,3))/(Q(i,j,1) + Q(i+1,j,1))
                    theta2 = atan(v2/u2)
                    p2   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                    t2   = p2/rho2
                    !cR1 = sqrt(GAMMA*p2/rho2)*sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)
                    cR1 = sqrt(GAMMA*p2/rho2)!*sqrt(mx**2.0d0 + my**2.0d0)
                    !VnR = (u2*m(i,j,1) + v2*m(i,j,2))
                    !VnR = u2*mx + v2*my
                    !Vabs2 = sqrt(u2**2.0d0 + v2**2.0d0) 

                    !if(i==47 .and. j==1)then
                    !    do cnt=1,90
                    !        VnL = u1*sin(1.0d0*cnt*PI/180.0d0) - v1*cos(1.0d0*cnt*PI/180.0d0)
                    !        VnR = u2*sin(1.0d0*cnt*PI/180.0d0) - v2*cos(1.0d0*cnt*PI/180.0d0)
                    !        write(*,*)cnt,VnL-cL1,u1,v1,VnR-cR1,u2,v2
                    !    enddo
                    !endif

                    prf(i,j) = (min(p1/p2,p2/p1))**3.0d0
                    !Flag(i,j) = prf(i,j)
                    
                    if(prf(i,j) < 0.5d0)then

                        ! Shock Detection by Entropy
                        !Entropy = log((t2/t1)**(GAMMA/(GAMMA-1.0d0))*((p1**rho1)/(p2**rho2)))
                        !if(Entropy > 0.01d0)then
                        !    !Flag(i,j) = 1.0d0
                        !endif
                        !Entropy = log((t1/t2)**(GAMMA/(GAMMA-1.0d0))*((p2**rho1)/(p1**rho1)))
                        !if(Entropy > 0.0d0)then
                        !    !Flag(i,j) = 1.0d0
                        !endif

                        ! Shock Detection by velocity

                        do l=1,Nx
                            if(prf(l,j-1) < 0.6d0)then
                                theta = atan((Y(i,j) - Y(l,j-1))/(X(i,j) - X(l,j-1)))
                                
                                if(theta < 0.0d0)then
                                    VnL = -u1*sin(theta) + v1*cos(theta)
                                    VnR = -u2*sin(theta) + v2*cos(theta)
                                else
                                    VnL = u1*sin(theta) - v1*cos(theta)
                                    VnR = u2*sin(theta) - v2*cos(theta)
                                endif
                                
                                if((VnL - cL1 > 0 .and. VnR - cR1 < 0) .or. (VnL + cL1 > 0 .and. VnR + cR1 <0))then
                                    if(sqrt(u1**2.0d0 + v1**2.0d0) < sqrt(u2**2.0d0 + v2**2.0d0))then
                                        Flag(i,j) = 0.0d0
                                    else
                                        Flag(i,j) = 1.0d0
                                        !write(*,*)'True',i,j,theta*180/PI,VnL,cL1,VnR,cR1
                                        write(*,*)'xi',i,j,theta*180.0d0/PI,u1,v1,u2,v2,&
                                        sqrt(u1**2.0d0 + v1**2.0d0),sqrt(u2**2.0d0 + v2**2.0d0)
                                    endif
                                endif
                                !write(*,*)i,j,theta*180.0d0/PI,Flag(i,j),VnL-cL1,VnR-cR1
                            endif
                        enddo
                    else
                        !Flag(i,j) = (min(p1/p2,p2/p1))**3.0d0
                        !Flag(i,j) = 1.0d0
                    endif

                    
                    theta = theta2 - theta1
                    tantheta = tan(theta)

                    !do cnt=-90,90
                    !    betadeg = betadeg + 1.0d0
                    !    beta = PI*betadeg/180.d0
                    !    betadeg2 = betadeg2 + 1.0d0
                    !    beta2 = PI*betadeg2/180.d0
                    !    hoge = (M1**2.0d0*sin(beta)**2.0d0 - 1.0d0)/(M1**2.0d0*(GAMMA+cos(2.0d0*beta)+2.0d0))*(2.0d0/tan(beta))
                    !    fuga = (M1**2.0d0*sin(beta2)**2.0d0 - 1.0d0)/(M1**2.0d0*(GAMMA+cos(2.0d0*beta2)+2.0d0))*(2.0d0/tan(beta2))
                    !    if((hoge-tantheta)*(fuga-tantheta) < 0.0d0)then
                    !        VnL = u1*sin(beta) - v1*cos(beta)
                    !        VnR = u2*sin(beta) - v2*cos(beta)
                    !        if((VnL-cL1>0.0d0 .and. VnR-cR1<0.0d0))then! .or. (VnL+cL1>0.0d0 .and. VnR+cR1<0.0d0))then
                    !            write(*,*)i,j,cnt,atan(beta)*180.d0/PI
                    !            Flag(i,j) = 1
                    !        endif
                    !    endif
                    !enddo
                    
                    betadeg = 0.0d0

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

                        !write(*,*)i,',',j,',',X(i,j),',',Y(i,j),',',180.0d0*atan((Y(i,j) - Y(iold,jold))/(X(i,j) - X(iold,jold)))/PI
                        !iold = i
                        !jold = j

                        ! Shock Detection wave by Entropy
                        !Entropy = log((t2/t1)**(GAMMA/(GAMMA-1.0d0))*((p1**rho1)/(p2**rho2)))
                        !if(Entropy > 0.0d0)then
                        !    !Flag(i,j) = 1.0d0
                        !endif
                        !Entropy = log((t1/t2)**(GAMMA/(GAMMA-1.0d0))*((p2**rho2)/(p1**rho1)))
                        !if(Entropy > 0.0d0)then
                        !    !Flag(i,j) = 1.0d0
                        !endif

                        ! Shock Detection Wave by compare Vn
                        do l=1,Ny
                            if(prf(i-1,l) < 0.5d0)then                                
                                theta = atan((Y(i,j) - Y(i-1,l))/(X(i,j) - X(i-1,l)))
                                if(theta < 0.0d0)then
                                    VnL = -u1*sin(theta) + v1*cos(theta)
                                    VnR = -u2*sin(theta) + v2*cos(theta)
                                else
                                    VnL = u1*sin(theta) - v1*cos(theta)
                                    VnR = u2*sin(theta) - v2*cos(theta)
                                endif
                            if((VnL - cL2 > 0 .and. VnR - cR2 < 0) .or. (VnL + cL2 > 0 .and. VnR + cR2 <0))then
                                    !if(sqrt(u1**2.0d0 + v1**2.0d0) )
                                    Flag(i,j) = 1.0d0
                                    write(*,*)'et',i,j,theta*180.0d0/PI,u1,v1,u2,v2,&
                                    sqrt(u1**2.0d0 + v1**2.0d0),sqrt(u2**2.0d0 + v2**2.0d0)

                                    !write(*,*)X(i,j),',',Y(i,j),',',prf(i,j),',',theta*180.0d0/PI
                                else
                                    !Flag(i,j) = 0.0d0
                                    !write(*,*)'0',i,j,prf(i,j),theta*180.0d0/PI
                                endif
                            endif
                        enddo
                    else
                        !Flag(i,j) = 1.0d0
                    endif

                    
                

                    !do cnt=-90,90
                    !    betadeg = betadeg + 1.0d0
                    !    beta = PI*betadeg/180.d0
                    !    betadeg2 = betadeg2 + 1.0d0
                    !    beta2 = PI*betadeg2/180.d0
                    !    hoge = (M1**2.0d0*sin(beta)**2.0d0 - 1.0d0)/(M1**2.0d0*(GAMMA+cos(2.0d0*beta)+2.0d0))*(2.0d0/tan(beta))
                    !    fuga = (M1**2.0d0*sin(beta2)**2.0d0 - 1.0d0)/(M1**2.0d0*(GAMMA+cos(2.0d0*beta2)+2.0d0))*(2.0d0/tan(beta2))
                    !    if((hoge-tantheta)*(fuga-tantheta) < 0.0d0)then
                    !        VnL = u1*sin(beta) - v1*cos(beta)
                    !        VnR = u2*sin(beta) - v2*cos(beta)
                    !        if((VnL-cL1>0.0d0 .and. VnR-cR1<0.0d0))then! .or. (VnL+cL1>0.0d0 .and. VnR+cR1<0.0d0))then
                    !            !write(*,*)i,j,cnt,atan(beta)*180.d0/PI
                    !            Flag(i,j) = 1
                    !        endif
                    !    endif
                    !enddo
                    
                    beta = 0.0d0
                enddo
            enddo
            
            do j=1,Ny
                do i=1,Nx
                    write(fo,"(f10.7)")Flag(i,j)
                enddo
            enddo  
            
            close(fo)
           

        enddo
        
            deallocate(X,Y)
            index = index + 1
end program visualize