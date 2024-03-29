module subprog
    use param
    implicit none
contains
    subroutine setMetrics(m,n,S,Nx,Ny,X,Y)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0
        real(8),intent(in) :: X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4)
        real(8),intent(out) :: m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
                               S(-2:Nx+3,-2:Ny+3)
        real(8) yeta,yxi,xeta,xxi
        character meshfile*64
        

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

    subroutine switchEnable(Nx,Ny,is_SF_xi,is_SF_eta,i,j)
        implicit none
        integer,intent(in) :: i,j,Nx,Ny
        integer,intent(out) :: is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2)

        is_SF_xi(i,j+1) = 2
        is_SF_xi(i-1,j+1) = 2
        is_SF_xi(i,j-1) = 2
        is_SF_xi(i-1,j-1) = 2
        is_SF_xi(i,j) = 1     
        is_SF_xi(i-1,j) = 1   
        is_SF_xi(i+1,j) = 2
        is_SF_xi(i-2,j) = 2
        is_SF_eta(i,j) = 1    
        is_SF_eta(i,j-1) = 1  
        is_SF_eta(i-1,j) = 2
        is_SF_eta(i-1,j-1) = 2
        is_SF_eta(i+1,j) = 2
        is_SF_eta(i+1,j-1) = 2
        is_SF_eta(i,j-2) = 2
        is_SF_eta(i,j+1) = 2

    end subroutine switchEnable

    subroutine ShockDetection(Q,Nx,Ny,is_SF_xi,is_SF_eta,X,Y)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0,l=0,cnt,idxBndUpper=0,idxBndLower=0
        integer,intent(out) :: is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2)
        !real(8),intent(out) :: is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2)
        real(8),intent(in) :: Q(-2:Nx+3,-2:Ny+3,4),X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4)
        real(8) :: rho1,u1,v1,p1,t1,rho2,u2,v2,p2,t2,VnL,VnR,cL,cR,beta,theta,f,f1,f2,g,mx,my
        real(8),allocatable :: prf(:,:)
        real(8),parameter :: PI=acos(-1.0d0)

        allocate(prf(-1:Nx+2,-1:Ny+2))

        do j=1,Ny
            do i=1,Nx
                is_SF_xi(i,j) = 0
                is_SF_eta(i,j) = 0
            enddo
        enddo

        do j=1,Ny-1
            do i=1,Nx-1
                rho1 = Q(i-1,j,1)
                u1   = Q(i-1,j,2)/Q(i-1,j,1)
                v1   = Q(i-1,j,3)/Q(i-1,j,1)
                p1   = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1))
                t1   = p1/rho1
                cL = sqrt(GAMMA*p1/rho1)
                rho2 = Q(i+1,j,1)
                u2   = Q(i+1,j,2)/Q(i+1,j,1)
                v2   = Q(i+1,j,3)/Q(i+1,j,1)
                p2   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                t2   = p2/rho2
                cR = sqrt(GAMMA*p2/rho2)
                
                
                prf(i,j) = (min(p1/p2,p2/p1))**3.0d0
                
                if(prf(i,j) < 0.5d0)then
                    ! Shock Detection by velocity
                    do l=1,Nx 
                        p1   = (GAMMA-1.0d0)*(Q(l,j-2,4) - 0.5d0*(Q(l,j-2,2)**2.0d0 + Q(l,j-2,3)**2.0d0)/Q(l,j-2,1))
                        p2   = (GAMMA-1.0d0)*(Q(l+1,j-2,4) - 0.5d0*(Q(l+1,j-2,2)**2.0d0 + Q(l+1,j-2,3)**2.0d0)/Q(l+1,j-2,1))
                        prf(l,j-1) = (min(p1/p2,p2/p1))**3.0d0
                        if(prf(l,j-2) < 0.5d0)then
                            
                            theta = atan((Y(i,j) - Y(l,j-2))/(X(i,j) - X(l,j-2)))
if(90.0d0 < theta .or. theta < -90.0d0)then
    GOTO 1000
endif                            
                            if(theta < 0.0d0)then
                                VnL = -u1*sin(theta) + v1*cos(theta)
                                VnR = -u2*sin(theta) + v2*cos(theta)
                            else
                                VnL = u1*sin(theta) - v1*cos(theta)
                                VnR = u2*sin(theta) - v2*cos(theta)
                            endif
                            
                            if((VnL - cL > 0 .and. VnR - cR < 0) .or. (VnL + cL > 0 .and. VnR + cR <0))then
                                if(sqrt(u1**2.0d0 + v1**2.0d0) < sqrt(u2**2.0d0 + v2**2.0d0))then
                                else
                                    call switchEnable(Nx,Ny,is_SF_xi,is_SF_eta,i,j)
                                endif
                            endif
                        endif
                    enddo
                else
                    !Flag(i,j) = 0.0d0
                endif
1000 CONTINUE

                rho1 = Q(i,j-1,1)
                u1   = Q(i,j-1,2)/Q(i,j-1,1)
                v1   = Q(i,j-1,3)/Q(i,j-1,1)
                p1   = (GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1))
                t1   = p1/rho1
                cL = sqrt(GAMMA*p1/rho1)
                
                rho2 = Q(i,j+1,1)
                u2   = Q(i,j+1,2)/Q(i,j+1,1)
                v2   = Q(i,j+1,3)/Q(i,j+1,1)
                p2   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
                t2   = p2/rho2
                cR = sqrt(GAMMA*p2/rho2)
                
                prf(i,j) = (min(p1/p2,p2/p1))**3.0d0
                if(prf(i,j) < 0.5d0)then
                    ! Shock Detection Wave by compare Vn
                    do l=1,Ny
                        p1   = (GAMMA-1.0d0)*(Q(i-2,l-1,4) - 0.5d0*(Q(i-2,l-1,2)**2.0d0 + Q(i-2,l-1,3)**2.0d0)/Q(i-2,l-1,1))
                        p2   = (GAMMA-1.0d0)*(Q(i-2,l+1,4) - 0.5d0*(Q(i-2,l+1,2)**2.0d0 + Q(i-2,l+1,3)**2.0d0)/Q(i-2,l+1,1))
                        prf(i-1,l) = (min(p1/p2,p2/p1))**3.0d0
                        
                        if(prf(i-2,l) < 0.5d0)then
                            theta = atan((Y(i,j) - Y(i-2,l))/(X(i,j) - X(i-2,l)))

if(90.0d0 < theta .or. theta < -90.0d0)then
    GOTO 1001
endif                            
                            if(theta < 0.0d0)then
                                VnL = -u1*sin(theta) + v1*cos(theta)
                                VnR = -u2*sin(theta) + v2*cos(theta)
                            else
                                VnL = u1*sin(theta) - v1*cos(theta)
                                VnR = u2*sin(theta) - v2*cos(theta)
                            endif
                        
                            if(v1 < 0.0d0)then
                                if((VnL - cL > 0 .and. VnR - cR < 0) .or. (VnL + cL > 0 .and. VnR + cR <0))then!Compresible sonic point
                                    if(sqrt(u1**2.0d0 + v1**2.0d0) < sqrt(u2**2.0d0 + v2**2.0d0))then!NOT accerelation
                                    else
                                        if(theta < 0.0d0)then 
                                            call switchEnable(Nx,Ny,is_SF_xi,is_SF_eta,i,j)
                                        endif
                                       
                                    endif
                                endif
                            else
                                if((VnR - cR > 0 .and. VnL - cL < 0) .or. (VnR + cR > 0 .and. VnL + cL <0))then!Compresible sonic point
                                    if(sqrt(u1**2.0d0 + v1**2.0d0) > sqrt(u2**2.0d0 + v2**2.0d0))then!NOT accerelation
                                        
                                    else
                                        if(theta > 0.0d0)then
                                            call switchEnable(Nx,Ny,is_SF_xi,is_SF_eta,i,j)
                                        endif
                                        
                                    endif
                                endif
                            endif                        
                            
                        endif
                    enddo
                else
                endif
1001 CONTINUE
            enddo
        enddo

        deallocate(prf)

    end subroutine ShockDetection

    function calcSpeed(Q,Nx,Ny,m,n,S) result(qmax)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0
        real(8),intent(in) ::Q(-2:Nx+3,-2:Ny+3,4),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
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