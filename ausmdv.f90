module AUSMDVFlux
    use param
    implicit none
contains
    subroutine AusmDV_XFlux(Q,E,m,Nx,Ny,i,j)
        implicit none
        integer,intent(in) :: i,j,Nx,Ny
        integer :: k=0,l=0
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

        mx = m(i,j,1)
        my = m(i,j,2)
        mmag = sqrt(mx**2+my**2)
        kx = -my
        ky = mx
        kmag = sqrt(kx**2+ky**2)

        !normalized metrics
        mx = mx/mmag!mx hat
        my = my/mmag!my hat
        kx = kx/kmag!kx hat
        ky = ky/kmag!ky hat
        
        !Higher resolution by MUSCL approach
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

        pL   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        hoge = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
        fuga = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i-1,j,1))
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

        pR   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
        hoge = (GAMMA-1.0d0)*(Q(i+2,j,4) - 0.5d0*(Q(i+2,j,2)**2.0d0 + Q(i+2,j,3)**2.0d0)/Q(i+2,j,1))
        fuga = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        tmp1 = hoge - pR
        tmp2 = pR - fuga
        fluxLimMin = minmod(tmp2,beta*tmp1)
        fluxLimPl  = minmod(tmp1,beta*tmp2)
        pR = pR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin

        alphaL = 2.0d0*(pL/rhoL)/((pL/rhoL) + (pR/rhoR))
        alphaR = 2.0d0*(pR/rhoR)/((pL/rhoL) + (pR/rhoR))
        cm     = max(sqrt(GAMMA*pL/rhoL),sqrt(GAMMA*pR/rhoR))         
        s      = 0.5d0*min(1.0d0,largeK*(abs(pR-pL)/min(pR,pL)))

        HL     = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0))/rhoL
        HR     = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0))/rhoR

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

    end subroutine AusmDV_XFlux

    subroutine AusmDV_YFlux(Q,F,n,Nxnum,Nynum,i,j)
        implicit none
        integer,intent(in) :: i,j,Nxnum,Nynum
        integer :: k=0,l=0
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
      
        nx = n(i,j,1)
        ny = n(i,j,2)
        nmag = sqrt(nx**2+ny**2)
        kx = ny
        ky = -nx
        kmag = sqrt(kx**2+ky**2)
        
        !normalized metrics
        nx = nx/nmag
        ny = ny/nmag
        kx = kx/kmag
        ky = ky/kmag
        
        !Higher resolution by MUSCL approach
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

        pL   = (GAMMA-1.0d0)*(Q(i,j,4) - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        hoge = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
        fuga = (GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1))
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

        pR   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
        hoge = (GAMMA-1.0d0)*(Q(i,j+2,4) - 0.5d0*(Q(i,j+2,2)**2.0d0 + Q(i,j+2,3)**2.0d0)/Q(i,j+2,1))
        fuga = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        tmp1 = hoge - pR
        tmp2 = pR - fuga
        fluxLimMin = minmod(tmp2,beta*tmp1)
        fluxLimPl  = minmod(tmp1,beta*tmp2)
        pR = pR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        
        alphaL = 2.0d0*(pL/rhoL)/((pL/rhoL) + (pR/rhoR))
        alphaR = 2.0d0*(pR/rhoR)/((pL/rhoL) + (pR/rhoR))
        cm     = max(sqrt(GAMMA*pL/rhoL),sqrt(GAMMA*pR/rhoR))         
        s      = 0.5d0*min(1.0d0,largeK*(abs(pR-pL)/min(pR,pL)))

        HL     = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0+VL**2.0d0))/rhoL
        HR     = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0+VR**2.0d0))/rhoR

        !calculation uLplus,uRminus,pLplus,pRminus            
        if(abs(VL)/cm <= 1)then
            ULp = 0.25d0*alphaL*((VL + cm)**2.0d0)/cm + 0.5d0*(1.0d0 - alphaL)*(VL + abs(VL))
            pLp = 0.25d0*pL*((VL/cm + 1.0d0)**2.0d0)*(2.0d0 - VL/cm)
        else
            ULp = 0.5d0*(VL + abs(VL))
            pLp = 0.5d0*pL*(VL + abs(VL))/VL
        endif

        if(abs(VR)/cm <= 1)then
            URm = -0.25d0*alphaR*((VR - cm)**2.0d0)/cm + 0.5d0*(1.0d0 - alphaR)*(VR - abs(VR))
            pRm = 0.25d0*pR*((VR/cm - 1.0d0)**2.0d0)*(2.0d0 + VR/cm)
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

    end subroutine AusmDV_YFlux
end module AUSMDVFlux