module param
    implicit none
    real(8) :: GAMMA = 1.4d0,Re = 5000.0d0
contains
    function minmod(x,y) result(fluxLimiter)
        real(8),intent(in) :: x,y
        real(8) fluxLimiter
        fluxLimiter = sign(1.0d0,x)*max(0.0d0,min(abs(x),y*sign(1.0d0,x)))
    end function minmod
end module param

module subprog
    use param
    implicit none
contains
    subroutine setMetrics(m,n,S,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0,Grid=12
        real(8),intent(out) :: m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
                               S(-2:Nx+3,-2:Ny+3)
                               !MetricsXx(i,j,1) = mx, MetricsXx(i,j,2) = my
                               !MetricsYy(i,j,1) = nx, MetricsYy(i,j,2) = ny
                               !
                               !         ↑Yy(nx,ny)
                               !      ------
                               !     |  i,j |
                               !     |      |→Xx(mx,my)
                               !      ------
        real(8) X(-3:Nx+4,-3:Ny+4),Y(-3:Nx+4,-3:Ny+4)
        
        !open(Grid,file = 'rampGrid.txt')
        open(Grid,file = 'Gridtube.txt')
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
                
                S(i,j) = 0.5d0*(n(i,j-1,2)+n(i,j,2))*0.5d0*(m(i-1,j,1)+m(i,j,1))-&
                         0.5d0*(n(i,j-1,1)+n(i,j,1))*0.5d0*(m(i-1,j,2)+m(i,j,2))
                
                !m(i,j,1) = 0.5d0*(Y(i,j+1) - Y(i,j-1))
                !m(i,j,2) = -0.5d0*(X(i,j+1) - X(i,j-1))
                !n(i,j,1) = -0.5d0*(Y(i+1,j) - Y(i-1,j))
                !n(i,j,2) = 0.5d0*(X(i+1,j) - X(i-1,j))
                !S(i,j) = (n(i,j,2)*m(i,j,1) - n(i,j,1)*m(i,j,2))

                !if(j<=3)then
                !    write(*,*)i,j,m(i,j,1),m(i,j,2),n(i,j,1),n(i,j,2),S(i,j)
                !endif
            enddo
        enddo
    end subroutine setMetrics

    subroutine evalXFlux(Q,E,m,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0,k=0,l=0
        real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4),m(-2:Nx+3,-2:Ny+3,2)
        real(8),intent(out) :: E(-1:Nx+2,-1:Ny+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0,largeK = 10.0d0
        real(8) :: beta=0,fluxLimMin=0,fluxLimPl=0, &
                rhoL=0,UL=0,VL=0,pL=0,rhoR=0,UR=0,VR=0,pR=0,&
                pLp=0,pRm=0,ULp=0,URm=0,alphaL=0,alphaR=0,HL=0,HR=0,cm=0,s=0,&
                massFlow=0,momentumV=0,momentumD=0,&
                hoge=0,fuga=0,tmp1=0,tmp2=0,&
                mx=0,my=0,kx=0,ky=0,mmag=0,kmag=0,&
                E1=0,E2=0,E3=0,E4=0
        
        beta = 0.5d0*(1.0d0 + (3.0d0-phi)/(1.0d0-phi))
        l = l+1
        do j=0,Ny    
            do i=0,Nx
                mx = m(i,j,1)
                my = m(i,j,2)
                mmag = sqrt(mx**2+my**2)
                kx = -my
                ky = mx
                kmag = sqrt(kx**2+ky**2)

                mx = mx/mmag!mx hat
                my = my/mmag!my hat
                kx = kx/kmag!kx hat
                ky = ky/kmag!ky hat
                
                rhoL = Q(i,j,1)
                tmp1 = Q(i+1,j,1) - Q(i,j,1)
                tmp2 = Q(i,j,1)   - Q(i-1,j,1)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                rhoL = rhoL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
    
                UL   = (Q(i,j,2)/Q(i,j,1))*mx + (Q(i,j,3)/Q(i,j,1))*my
                tmp1 = Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my&
                       - (Q(i,j,2)/Q(i,j,1)*mx + Q(i,j,3)/Q(i,j,1)*my)
                tmp2 = Q(i,j,2)/Q(i,j,1)*mx + Q(i,j,3)/Q(i,j,1)*my&
                       - (Q(i-1,j,2)/Q(i-1,j,1)*mx + Q(i-1,j,3)/Q(i-1,j,1)*my)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                UL = UL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl

                VL   = (Q(i,j,2)/Q(i,j,1))*kx + (Q(i,j,3)/Q(i,j,1))*ky
                tmp1 = Q(i+1,j,2)/Q(i+1,j,1)*kx + Q(i+1,j,3)/Q(i+1,j,1)*ky&
                       - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky)
                tmp2 = Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky&
                       - (Q(i-1,j,2)/Q(i-1,j,1)*kx + Q(i-1,j,3)/Q(i-1,j,1)*ky)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                VL = VL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl

                pL   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*(Q(i,j,2)**2 + Q(i,j,3)**2)/Q(i,j,1))
                hoge = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2 + Q(i+1,j,3)**2)/Q(i+1,j,1))
                fuga = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2 + Q(i+1,j,3)**2)/Q(i-1,j,1))
                tmp1 = hoge-pL
                tmp2 = pL-fuga
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                pL = pL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
                    
                rhoR = Q(i+1,j,1)
                tmp1 = Q(i+2,j,1) - Q(i+1,j,1)
                tmp2 = Q(i+1,j,1) - Q(i,j,1)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                rhoR = rhoR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

                UR = Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my
                tmp1 = Q(i+2,j,2)/Q(i+2,j,1)*mx + Q(i+2,j,3)/Q(i+2,j,1)*my& 
                        - (Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my)
                tmp2 = Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my&
                        - (Q(i,j,2)/Q(i,j,1)*mx + Q(i,j,3)/Q(i,j,1)*my)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                UR = UR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

                VR = (Q(i+1,j,2)/Q(i+1,j,1))*kx + (Q(i+1,j,3)/Q(i+1,j,1))*ky
                tmp1 = Q(i+2,j,2)/Q(i+2,j,1)*kx + Q(i+2,j,3)/Q(i+2,j,1)*ky&
                        - (Q(i+1,j,2)/Q(i+1,j,1)*kx + Q(i+1,j,3)/Q(i+1,j,1)*ky)
                tmp2 = Q(i+1,j,2)/Q(i+1,j,1)*kx + Q(i+1,j,3)/Q(i+1,j,1)*ky&
                        - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                VR = VR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

                pR   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2 + Q(i+1,j,3)**2)/Q(i+1,j,1))
                hoge = (GAMMA-1.0d0)*(Q(i+2,j,4) - 0.5d0*(Q(i+2,j,2)**2 + Q(i+2,j,3)**2)/Q(i+2,j,1))
                fuga = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2 + Q(i,j,3)**2)/Q(i,j,1))
                tmp1 = hoge - pR
                tmp2 = pR - fuga
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                pR = pR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

                alphaL = 2.0d0*(pL/rhoL)/((pL/rhoL) + (pR/rhoR))
                alphaR = 2.0d0*(pR/rhoR)/((pL/rhoL) + (pR/rhoR))
                cm     = max(sqrt(GAMMA*pL/rhoL),sqrt(GAMMA*pR/rhoR))         
                s      = 0.5d0*min(1.0d0,largeK*(abs(pR-pL)/min(pR,pL)))

                HL     = (GAMMA*pL/(GAMMA-1) + 0.5d0*rhoL*(UL**2+VL**2))/rhoL
                HR     = (GAMMA*pR/(GAMMA-1) + 0.5d0*rhoR*(UR**2+VR**2))/rhoR

                !calculation uLplus,uRminus,pLplus,pRminus            
                if(abs(UL)/cm <= 1)then
                    ULp = 0.25d0*alphaL*((UL + cm)**2)/cm + 0.5d0*(1.0d0 - alphaL)*(UL + abs(UL))
                    pLp = 0.25d0*pL*((UL/cm + 1.0d0)**2)*(2.0d0 - UL/cm)
                else
                    ULp = 0.5d0*(UL + abs(UL))
                    pLp = 0.5d0*pL*(UL + abs(UL))/UL
                endif

                if(abs(UR)/cm <= 1)then
                    URm = -0.25d0*alphaR*((UR - cm)**2)/cm + 0.5d0*(1.0d0 - alphaR)*(UR - abs(UR))
                    pRm = 0.25d0*pR*((UR/cm - 1.0d0)**2)*(2.0d0 + UR/cm)
                else
                    URm = 0.5d0*(UR - abs(UR))
                    pRm = 0.5d0*pR*(UR - abs(UR))/UR
                endif

                massFlow  = ULp*rhoL + URm*rhoR
                momentumV = ULp*rhoL*UL + URm*rhoR*UR
                momentumD = 0.5d0*(massFlow*(UL + UR) - abs(massFlow)*(UR - UL))

                E1 = massFlow
                E2 = (0.5d0 + s)*momentumV + (0.5d0 - s)*momentumD + pLp + pRm
                E3 = 0.5d0*(massFlow*(VL+VR) - abs(massFlow)*(VR-VL))
                E4 = 0.5d0*(massFlow*(HL+HR) - abs(massFlow)*(HR-HL))
                E(i,j,1) = sqrt(m(i,j,1)**2 + m(i,j,2)**2)*E1
                E(i,j,2) = sqrt(m(i,j,1)**2 + m(i,j,2)**2)*(mx*E2 + kx*E3)
                E(i,j,3) = sqrt(m(i,j,1)**2 + m(i,j,2)**2)*(my*E2 + ky*E3)
                E(i,j,4) = sqrt(m(i,j,1)**2 + m(i,j,2)**2)*E4
            enddo
        enddo
    end subroutine evalXFlux

    subroutine evalYFlux(Q,F,n,Nxnum,Nynum)
        implicit none
        integer,intent(in) :: Nxnum,Nynum
        integer :: i=0,j=0,k=0,l=0
        real(8),intent(in) :: Q(-1:Nxnum+2,-1:Nynum+2,4),n(-2:Nxnum+3,-2:Nynum+3,2)
        real(8),intent(out) :: F(-1:Nxnum+2,-1:Nynum+2,4)
        real(8),parameter :: phi = 1.0/3.0,largeK = 10.0
        real(8) :: beta=0,fluxLimMin=0,fluxLimPl=0, &
                rhoL=0,UL=0,VL=0,pL=0,rhoR=0,UR=0,VR=0,pR=0,&
                pLp=0,pRm=0,ULp=0,URm=0,alphaL=0,alphaR=0,HL=0,HR=0,cm=0,s=0,&
                massFlow=0,momentumV=0,momentumD=0,&
                hoge=0,fuga=0,tmp1=0,tmp2=0,&
                nx=0,ny=0,kx=0,ky=0,nmag=0,kmag=0,&
                F1=0,F2=0,F3=0,F4=0
        
        l=l+1
        beta = 0.5d0*(1.0d0 + (3.0d0-phi)/(1.0d0-phi))

        do j=0,Nynum
            do i=0,Nxnum                
                nx = n(i,j,1)
                ny = n(i,j,2)
                nmag = sqrt(nx**2+ny**2)
                kx = ny
                ky = -nx
                kmag = sqrt(kx**2+ky**2)
                
                nx = nx/nmag
                ny = ny/nmag
                kx = kx/kmag
                ky = ky/kmag

                rhoL = Q(i,j,1)
                tmp1 = Q(i,j+1,1) - Q(i,j,1)
                tmp2 = Q(i,j,1)   - Q(i,j-1,1)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                rhoL = rhoL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
    
                UL   = (Q(i,j,2)/Q(i,j,1))*kx + (Q(i,j,3)/Q(i,j,1))*ky
                tmp1 = Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j,1)*ky&
                        - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky)
                tmp2 = Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky&
                        - (Q(i,j-1,2)/Q(i,j-1,1)*kx + Q(i,j-1,3)/Q(i,j-1,1)*ky)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                UL = UL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl

                VL   = (Q(i,j,2)/Q(i,j,1))*nx + (Q(i,j,3)/Q(i,j,1))*ny
                tmp1 = Q(i,j+1,2)/Q(i,j+1,1)*nx + Q(i,j+1,3)/Q(i,j,1)*ny&
                        - (Q(i,j,2)/Q(i,j,1)*nx + Q(i,j,3)/Q(i,j,1)*ny)
                tmp2 = Q(i,j,2)/Q(i,j,1)*nx + Q(i,j,3)/Q(i,j,1)*ny&
                        - (Q(i,j-1,2)/Q(i,j-1,1)*nx + Q(i,j-1,3)/Q(i,j-1,1)*ny)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                VL = VL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl

                pL   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*(Q(i,j,2)**2 + Q(i,j,3)**2)/Q(i,j,1))
                hoge = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2 + Q(i,j+1,3)**2)/Q(i,j+1,1))
                fuga = (GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2 + Q(i,j-1,3)**2)/Q(i,j-1,1))
                tmp1 = hoge-pL
                tmp2 = pL-fuga
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                pL = pL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl

                rhoR = Q(i,j+1,1)
                tmp1 = Q(i,j+2,1) - Q(i,j+1,1)
                tmp2 = Q(i,j+1,1) - Q(i,j,1)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                rhoR = rhoR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

                UR   = Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j+1,1)*ky
                tmp1 = Q(i,j+2,2)/Q(i,j+2,1)*kx + Q(i,j+2,3)/Q(i,j+2,1)*ky&
                         - (Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j+1,1)*ky)
                tmp2 = Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j+1,1)*ky&
                         - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                UR = UR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

                VR = (Q(i,j+1,2)/Q(i,j+1,1))*nx + (Q(i,j+1,3)/Q(i,j+1,1))*ny
                tmp1 = Q(i,j+2,2)/Q(i,j+2,1)*nx + Q(i,j+2,3)/Q(i,j+2,1)*ny&
                         - (Q(i,j+1,2)/Q(i,j+1,1)*nx + Q(i,j+1,3)/Q(i,j+1,1)*ny)
                tmp2 = Q(i,j+1,2)/Q(i,j+1,1)*nx + Q(i,j+1,3)/Q(i,j+1,1)*ny&
                         - (Q(i,j,2)/Q(i,j,1)*nx + Q(i,j,3)/Q(i,j,1)*ny)
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                VR = VR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

                pR   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2 + Q(i,j+1,3)**2)/Q(i,j+1,1))
                hoge = (GAMMA-1.0d0)*(Q(i,j+2,4) - 0.5d0*(Q(i,j+2,2)**2 + Q(i,j+2,3)**2)/Q(i,j+2,1))
                fuga = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2 + Q(i,j,3)**2)/Q(i,j,1))
                tmp1 = hoge - pR
                tmp2 = pR - fuga
                fluxLimMin = minmod(tmp2,beta*tmp1)
                fluxLimPl  = minmod(tmp1,beta*tmp2)
                pR = pR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
                
                alphaL = 2.0d0*(pL/rhoL)/((pL/rhoL) + (pR/rhoR))
                alphaR = 2.0d0*(pR/rhoR)/((pL/rhoL) + (pR/rhoR))
                cm     = max(sqrt(GAMMA*pL/rhoL),sqrt(GAMMA*pR/rhoR))         
                s      = 0.5d0*min(1.0d0,largeK*(abs(pR-pL)/min(pR,pL)))

                HL     = (GAMMA*pL/(GAMMA-1) + 0.5d0*rhoL*(UL**2+VL**2))/rhoL
                HR     = (GAMMA*pR/(GAMMA-1) + 0.5d0*rhoR*(UR**2+VR**2))/rhoR

                !calculation uLplus,uRminus,pLplus,pRminus            
                if(abs(VL)/cm <= 1)then
                    ULp = 0.25d0*alphaL*((VL + cm)**2)/cm + 0.5d0*(1.0d0 - alphaL)*(VL + abs(VL))
                    pLp = 0.25d0*pL*((VL/cm + 1.0d0)**2)*(2.0d0 - VL/cm)
                else
                    ULp = 0.5d0*(VL + abs(VL))
                    pLp = 0.5d0*pL*(VL + abs(VL))/VL
                endif

                if(abs(VR)/cm <= 1)then
                    URm = -0.25d0*alphaR*((VR - cm)**2)/cm + 0.5d0*(1.0d0 - alphaR)*(VR - abs(VR))
                    pRm = 0.25d0*pR*((VR/cm - 1.0d0)**2)*(2.0d0 + VR/cm)
                else
                    URm = 0.5d0*(VR - abs(VR))
                    pRm = 0.5d0*pR*(VR - abs(VR))/VR
                endif
                
                massFlow  = ULp*rhoL + URm*rhoR
                momentumV = ULp*rhoL*VL + URm*rhoR*VR
                momentumD = 0.5d0*(massFlow*(VL + VR) - abs(massFlow)*(VR - VL))

                F1 = massFlow
                F2 = 0.5d0*(massFlow*(UL+UR) - abs(massFlow)*(UR-UL))
                F3 = (0.5d0 + s)*momentumV + (0.5d0 - s)*momentumD + pLp + pRm
                F4 = 0.5d0*(massFlow*(HL+HR) - abs(massFlow)*(HR-HL))
                F(i,j,1) = sqrt(n(i,j,1)**2 + n(i,j,2)**2)*F1
                F(i,j,2) = sqrt(n(i,j,1)**2 + n(i,j,2)**2)*(kx*F2 + nx*F3)
                F(i,j,3) = sqrt(n(i,j,1)**2 + n(i,j,2)**2)*(ky*F2 + ny*F3)
                F(i,j,4) = sqrt(n(i,j,1)**2 + n(i,j,2)**2)*F4
            enddo
        enddo
    end subroutine evalYFlux

    function calcSpeed(Q,Nx,Ny,m,n,S) result(qmax)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0
        real(8),intent(in) ::Q(-1:Nx+2,-1:Ny+2,4),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
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
                Cxi = 1.0d0/S(i,j)*(abs(U1) + c*sqrt(m(i,j,1)**2 + m(i,j,2)**2))
                Ceta = 1.0d0/S(i,j)*(abs(U2) + c*sqrt(n(i,j,1)**2 + n(i,j,2)**2))

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

    subroutine output2file(Q,Nx,Ny,index)
        implicit none
        integer,intent(in) :: Nx,Ny,index
        integer :: fo=20,i=0,j=0,Grid=19
        real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4)
        character filename*128

        write(filename,'("Qbin",i3.3,".dat")')index
        !open(fo,file='Qbin.dat')
        open(fo,file = filename)
        do j=-1,Ny+2
            do i=-1,Nx+2
                write(fo,*)Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
            enddo
        enddo
        close(fo)
    end subroutine output2file

end module subprog

program ausmdv
    use subprog
    implicit none
    integer :: Nx,Ny,i=0,j=0,k=0,time=0,GridNum=11,foo=12,index=1
    real(8),allocatable :: Q(:,:,:),Qn(:,:,:),Qast(:,:,:),E(:,:,:),F(:,:,:),S(:,:),Ev(:,:,:),Fv(:,:,:)
    real(8),allocatable :: m(:,:,:),n(:,:,:)
    real(8) dt,qmax
    real(8),parameter :: CFL = 0.15d0
    real(8) :: elapsedTime=0.0d0,Q1old=0.0d0,Q2old=0.0d0,res=0.0d0,ressum=0.0d0

    open(GridNum,file = 'GridNum.txt')
    read(GridNum,*)Nx,Ny
    close(GridNum)
    allocate(Q(-1:Nx+2,-1:Ny+2,4),Qn(-1:Nx+2,-1:Ny+2,4),Qast(-1:Nx+2,-1:Ny+2,4),E(-1:Nx+2,-1:Ny+2,4),F(-1:Nx+2,-1:Ny+2,4),&
                S(-2:Nx+3,-2:Ny+3),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),&
                Ev(-1:Nx+2,-1:Ny+2,4),Fv(-1:Nx+2,-1:Ny+2,4))
        
    call setInitCond(Q,Nx,Ny)
    call setMetrics(m,n,S,Nx,Ny)
    
    open(foo,file = 'residualerror.txt')

    do time=1,1!500!time
        
        call setbdycond(Q,m,n,Nx,Ny)
        call evalXFlux(Q,E,m,Nx,Ny)
        call Xviscflux(Q,Ev,m,n,S,Nx,Ny)
        call evalYFlux(Q,F,n,Nx,Ny)
        call YviscFlux(Q,Fv,m,n,S,Nx,Ny)
        
        qmax = calcSpeed(Q,Nx,Ny,m,n,S)
        dt =CFL/qmax

        do j=1,Ny
            do i=1,Nx !space
                do k=1,4
                    Qn(i,j,k) = Q(i,j,k)
                    Q(i,j,k) = Q(i,j,k) - 0.50d0*(dt/S(i,j))*(E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k)& 
                                                        - (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re)
                enddo
            enddo
        enddo

        do j=-1,Ny+2
            do i=-1,Nx+2
                do k=1,4
                    Qast(i,j,k) = Q(i,j,k)
                enddo
            enddo
        enddo

        !call setbdycond(Qast,m,n,Nx,Ny)

        call evalXFlux(Qast,E,m,Nx,Ny)
        call Xviscflux(Qast,Ev,m,n,S,Nx,Ny)
        call evalYFlux(Qast,F,n,Nx,Ny)
        call YviscFlux(Qast,Fv,m,n,S,Nx,Ny)

        do j=1,Ny
            do i=1,Nx !space
                Q1old = Qn(i,j,1)
                do k=1,4
                    Q(i,j,k) = Qn(i,j,k) - (dt/S(i,j))*(E(i,j,k)-E(i-1,j,k) + F(i,j,k)-F(i,j-1,k)& 
                                                        - (Ev(i,j,k)-Ev(i-1,j,k))/Re - (Fv(i,j,k)-Fv(i,j-1,k))/Re)
                enddo

                res = abs(Q(i,j,1) - Q1old)
                ressum = ressum + res
            enddo
        enddo
        ressum = ressum/(Nx*Ny)
        elapsedTime = elapsedTime + dt
        
        if(modulo(time,500) == 0)then
        !    write(foo,*)elapsedTime,ressum
            call output2file(Q,Nx,Ny,index)
            index = index + 1
        endif

        write(*,*)ressum,elapsedTime,dt

        if(ressum < 0.0000005d0)then
            exit
        endif
        ressum = 0.0d0
        !time = time + 1
    enddo

    close(foo)

    call output2file(Q,Nx,Ny,index)

    deallocate(Q,E,F,S,m,n,Ev,Fv)

end program ausmdv