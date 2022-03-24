module SD_SLAUFlux
    use param
    implicit none
contains
    subroutine SD_SLAU_XFlux(Q,E,m,Nx,Ny,i,j,b)
        implicit none
        integer,intent(in) :: i,j,Nx,Ny
        integer :: n=0,foo=0,k=1
        real(8),intent(in) :: Q(-2:Nx+3,-2:Ny+3,4),m(-2:Nx+3,-2:Ny+3,2),b
        real(8),intent(out) :: E(-1:Nx+2,-1:Ny+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0,Csd1=0.10d0,Csd2=10.0d0,&
                            CL1=1.0d0/16.0d0,CL2=10.0d0/16.0d0,CL3=5.0d0/16.0d0,epsiron= 0.000001d0
        real(8) rhoL,rhoR,rhoBar,UL,u_l,UR,u_r,VL,v_l,VR,v_r,pL,pR,HL,HR,HBar,cL,cR,cBar,cStar,cHalf,MP,MM,MHat,VnL,VnR,MBar,&
                absVnBar,absVnBar_M,absVnBar_P,chi,g,beta_P,beta_M,massFlow,Ptilde,thetaSD,chiSD,p1,p2,p3,p4,p5,p6,dPmax,PBar,&
                fluxLmtP,fluxlmtM,beta,deltaP,deltaM,hoge,fuga,f_(3),s_(3),IS(3),bWCNS(3),omega(3),qtilde(3),h,qtmp,&
                mx,my,kx,ky,mmag,kmag,E1,E2,E3,E4
        
        n=n+1
        h = 1.0d0/Nx
        beta = b!0.5d0*(3.0d0 - phi)/(1.0d0 - phi)
        
        mx = m(i,j,1)
        my = m(i,j,2)
        mmag = sqrt(mx**2.0d0 + my**2.0d0)
        kx = -my
        ky = mx
        kmag = sqrt(kx**2.0d0 + ky**2.0d0)

        mx = mx/mmag!mx hat
        my = my/mmag!my hat
        kx = kx/kmag!kx hat
        ky = ky/kmag!ky hat

        !!!!!!!!! WCNS !!!!!!!!!!!!!!
        !do k =1,4   
        !    f_(1) = 0.5d0*(Q(i-2,j,k) - 4.0d0*Q(i-1,j,k) + 3.0d0*Q(i,j,k))
        !    f_(2) = 0.5d0*(Q(i+1,j,k) - Q(i-1,j,k))
        !    f_(3) = 0.5d0*(-3.0d0*Q(i,j,k) + 4.0d0*Q(i+1,j,k) - Q(i+2,j,k))
        !    s_(1) = (Q(i-2,j,k) - 2.0d0*Q(i-1,j,k) + Q(i,j,k))
        !    s_(2) = (Q(i-1,j,k) - 2.0d0*Q(i,j,k)   + Q(i+1,j,k))
        !    s_(3) = (Q(i,j,k)   - 2.0d0*Q(i+1,j,k) + Q(i+2,j,k))
        !    IS(1) = f_(1)**2.0d0 + s_(1)**2.0d0
        !    IS(2) = f_(2)**2.0d0 + s_(2)**2.0d0
        !    IS(3) = f_(3)**2.0d0 + s_(3)**2.0d0
        !    bWCNS(1) = CL1/((epsiron + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((epsiron + IS(2))**2.0d0)
        !    bWCNS(3) = CL3/((epsiron + IS(3))**2.0d0)
        !    omega(1) = bWCNS(1)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(2) = bWCNS(2)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(3) = bWCNS(3)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    qtilde(1) = Q(i,j,k) + 0.5d0*f_(1) + 0.125d0*s_(1)
        !    qtilde(2) = Q(i,j,k) + 0.5d0*f_(2) + 0.125d0*s_(2)
        !    qtilde(3) = Q(i,j,k) + 0.5d0*f_(3) + 0.125d0*s_(3)
        !    if(k == 1)then
        !        rhoL = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !    else if(k == 2)then
        !        u_l = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoL
        !    else if(k == 3)then
        !        v_l = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoL
        !    else if(k == 4)then
        !        qtmp = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !        pL = (GAMMA - 1.0d0)*(qtmp - 0.5d0*rhoL*(u_l**2.0d0 + v_l**2.0d0))
        !    endif
        !    f_(1) = 0.5d0*(Q(i-1,j,k) - 4.0d0*Q(i,j,k) + 3.0d0*Q(i+1,j,k))
        !    f_(2) = 0.5d0*(Q(i+2,j,k) - Q(i,j,k))
        !    f_(3) = 0.5d0*(-3.0d0*Q(i+1,j,k) + 4.0d0*Q(i+2,j,k) - Q(i+3,j,k))
        !    s_(1) = (Q(i-1,j,k) - 2.0d0*Q(i,j,k)   + Q(i+1,j,k))
        !    s_(2) = (Q(i,j,k)   - 2.0d0*Q(i+1,j,k) + Q(i+2,j,k))
        !    s_(3) = (Q(i+1,j,k) - 2.0d0*Q(i+2,j,k) + Q(i+3,j,k))
        !    IS(1) = f_(1)**2.0d0 + s_(1)**2.0d0
        !    IS(2) = f_(2)**2.0d0 + s_(2)**2.0d0
        !    IS(3) = f_(3)**2.0d0 + s_(3)**2.0d0
        !    bWCNS(1) = CL3/((epsiron + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((epsiron + IS(2))**2.0d0)
        !    bWCNS(3) = CL1/((epsiron + IS(3))**2.0d0)
        !    omega(1) = bWCNS(1)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(2) = bWCNS(2)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(3) = bWCNS(3)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    qtilde(1) = Q(i+1,j,k) - 0.5d0*f_(1) + 0.125d0*s_(1)
        !    qtilde(2) = Q(i+1,j,k) - 0.5d0*f_(2) + 0.125d0*s_(2)
        !    qtilde(3) = Q(i+1,j,k) - 0.5d0*f_(3) + 0.125d0*s_(3)
        !    if(k == 1)then
        !        rhoR = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !    else if(k == 2)then
        !        u_r = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoR
        !    else if(k == 3)then
        !        v_r = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoR
        !    else if(k == 4)then
        !        qtmp = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !        pR = (GAMMA - 1.0d0)*(qtmp - 0.5d0*rhoR*(u_r**2.0d0 + v_r**2.0d0))
        !    endif
        !enddo

        !!!!!!!! MUSCL
        deltaP = Q(i+1,j,1) - Q(i,j,1)
        deltaM = Q(i,j,1) - Q(i-1,j,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)
        rhoL = Q(i,j,1) + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
        !
        deltaP = Q(i+2,j,1) - Q(i+1,j,1)
        deltaM = Q(i+1,j,1) - Q(i,j,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)
        rhoR = Q(i+1,j,1) - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
        !
        u_l = Q(i,j,2)/Q(i,j,1)
        deltaP = Q(i+1,j,2)/Q(i+1,j,1) - Q(i,j,2)/Q(i,j,1) 
        deltaM = Q(i,j,2)/Q(i,j,1) - Q(i-1,j,2)/Q(i-1,j,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)            
        u_l   = u_l + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
        !
        v_l = Q(i,j,3)/Q(i,j,1)
        deltaP = Q(i+1,j,3)/Q(i+1,j,1) - Q(i,j,3)/Q(i,j,1)
        deltaM = Q(i,j,3)/Q(i,j,1) - Q(i-1,j,3)/Q(i-1,j,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)
        v_l   = v_l + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
        !
        u_r = Q(i+1,j,2)/Q(i+1,j,1)                
        deltaP = Q(i+2,j,2)/Q(i+2,j,1) - Q(i+1,j,2)/Q(i+1,j,1)
        deltaM = Q(i+1,j,2)/Q(i+1,j,1) - Q(i,j,2)/Q(i,j,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)            
        u_r   = u_r - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
        !
        v_r = Q(i+1,j,3)/Q(i+1,j,1)
        deltaP = Q(i+2,j,3)/Q(i+2,j,1) - Q(i+1,j,3)/Q(i+1,j,1)
        deltaM = Q(i+1,j,3)/Q(i+1,j,1) - Q(i,j,3)/Q(i,j,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)
        v_r = v_r - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM    
        !
        hoge   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
        pL     = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        fuga   = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1))
        deltaP = hoge - pL
        deltaM = pL - fuga
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)            
        pL   = pL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
        !
        hoge   = (GAMMA-1.0d0)*(Q(i+2,j,4) - 0.5d0*(Q(i+2,j,2)**2.0d0 + Q(i+2,j,3)**2.0d0)/Q(i+2,j,1))
        pR     = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
        fuga   = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        deltaP = hoge - pR
        deltaM = pR - fuga
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)                
        pR   = pR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM


        rhoBar = 0.5d0*(rhoL + rhoR)
        PBar = 0.5d0*(pL + pR)
        UL = u_l*mx + v_l*my
        UR = u_r*mx + v_r*my
        VL = u_l*kx + v_l*ky
        VR = u_r*kx + v_r*ky
        VnL  = UL
        VnR  = UR

        HL   = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0))/rhoL
        HR   = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0))/rhoR
        !HBar = 0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0) + (GAMMA/(GAMMA -1.0d0))*PBar/rhoBar
        HBar = 0.5d0*(HL+HR)

        cStar = sqrt(2.0d0*(GAMMA - 1.0d0)*HBar/(GAMMA + 1.0d0))
        !cL   = (cStar**2.0d0)/max(cStar,abs(VnL)) !SLAU2
        !cR   = (cStar**2.0d0)/max(cStar,abs(VnR)) !SLAU2   
        cL   = sqrt(GAMMA*pL/rhoL)               !SLAU
        cR   = sqrt(GAMMA*pR/rhoR)               !SLAU
        cBar = 0.5d0*(cL + cR)
        cHalf = min(cL,cR)
        
        MP   = VnL/cBar
        MM   = VnR/cBar
        
        ! mass flux parameter for SLAU/SLAU2
        absVnBar  = (rhoL*abs(VnL) + rhoR*abs(VnR))/(rhoL+rhoR)
        MHat  = min(1.0d0,sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))/cBar)
        chi = (1.0d0 - MHat)**2.0d0                                                      
        g = - max(min(MP,0.0d0),-1.0d0) * min(max(MM,0.0d0),1.0d0)

        ! mass flux parameter for SD-SLAU
        p1 = abs((GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1)) - pL)
        p2 = abs((GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1)) - pL)
        p3 = abs((GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1)) - pL)
        p4 = abs((GAMMA-1.0d0)*(Q(i+1,j-1,4) - 0.5d0*(Q(i+1,j-1,2)**2.0d0 + Q(i+1,j-1,3)**2.0d0)/Q(i+1,j-1,1)) - pR)
        p5 = abs((GAMMA-1.0d0)*(Q(i+2,j,4) - 0.5d0*(Q(i+2,j,2)**2.0d0 + Q(i+2,j,3)**2.0d0)/Q(i+2,j,1)) - pR)
        p6 = abs((GAMMA-1.0d0)*(Q(i+1,j+1,4) - 0.5d0*(Q(i+1,j+1,2)**2.0d0 + Q(i+1,j+1,3)**2.0d0)/Q(i+1,j+1,1)) - pR)
        dPmax = max(p1,p2,p3,p4,p5,p6)
        MBar = ((rhoL*abs(UL) + rhoR*abs(UR))/(rhoL + rhoR))/cBar
        thetaSD = min(1.0d0,((Csd2*(pR-pL) + Csd1*PBar)/(dPmax + Csd1*PBar))**2.0d0)*max(0.0d0,1.0d0 - abs(MBar))
        chiSD   = 0.5d0*thetaSD*(abs(MBar + 1.0d0) + abs(MBar - 1.0d0) - 2.0d0*abs(MBar))

        !massFlow = 0.5d0*(rhoL*VnL + rhoR*VnR - absVnBar*(rhoR - rhoL))*(1.0d0 - g) - 0.5d0*chiSD*(pR - pL)/cBar
        massFlow = 0.5d0*(rhoL*VnL + rhoR*VnR - absVnBar*(rhoR - rhoL))*(1.0d0 - g) - 0.5d0*chi*(pR - pL)/cBar

        if(1.0d0 <= abs(MP))then
            beta_P = 0.5d0*(1.0d0 + sign(1.0d0,MP))
        else
            beta_P = 0.25d0*(2.0d0 - MP)*(MP + 1.0d0)**2.0d0
        endif

        if(1.0d0 <= abs(MM))then
            beta_M = 0.5d0*(1.0d0 + sign(1.0d0,-MM))
        else
            beta_M = 0.25d0*(2.0d0 + MM)*(MM - 1.0d0)**2.0d0
        endif

        ! pressure flux parameter for SLAU2
        !absVnBar_P = (1.0d0 - g)*absVnBar + g*abs(VnL)
        !absVnBar_M = (1.0d0 - g)*absVnBar + g*abs(VnR)             
        !cHalf = min(cL,cR)
        !write(*,*)i,j,'X',HBar,cStar,cl,cr,cHalf

        
        Ptilde = 0.5d0*(pL + pR) + 0.5d0*(beta_P - beta_M)*(pL - pR)&
                    !+ (beta_P + beta_M - 1.0d0)*sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))*rhoBar*cHalf!SLAU2
                    + 0.5d0*(1.0d0 - chi)*(beta_P + beta_M - 1.0d0)*(pL+pR)                   !SLAU

        E1 = massFlow
        E2 = 0.5d0*(massFlow*(UR + UL) - abs(massFlow)*(UR - UL)) + Ptilde
        E3 = 0.5d0*(massFlow*(VR + VL) - abs(massFlow)*(VR - VL))
        E4 = 0.5d0*(massFlow*(HR + HL) - abs(massFlow)*(HR - HL))
        !write(*,*)i,j,'E',E1,E2,E3,E4
        
        E(i,j,1) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*E1!*(1.0d0 - is_SF_xi(i,j))
        E(i,j,2) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*(mx*E2 + kx*E3)!*(1.0d0 - is_SF_xi(i,j))
        E(i,j,3) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*(my*E2 + ky*E3)!*(1.0d0 - is_SF_xi(i,j))
        E(i,j,4) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*E4!*(1.0d0 - is_SF_xi(i,j))
        !write(*,*)i,E(i,j,1),E(i,j,2),E(i,j,3),E(i,j,4)
        !if(i<50 .and. j==1)then
        !    write(*,*)i,rhoL,u_l,v_l,pL
        !endif

    end subroutine SD_SLAU_XFlux

    subroutine SD_SLAU_YFlux(Q,F,n,Nxnum,Nynum,i,j,b)
        implicit none
        integer,intent(in) :: i,j,Nxnum,Nynum
        integer :: iter = 0,k=1
        real(8),intent(in) :: Q(-2:Nxnum+3,-2:Nynum+3,4),n(-2:Nxnum+3,-2:Nynum+3,2),b
        real(8),intent(out) :: F(-1:Nxnum+2,-1:Nynum+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0,Csd1=0.1d0,Csd2=10.0d0,&
                            CL1=1.0d0/16.0d0,CL2=10.0d0/16.0d0,CL3=5.0d0/16.0d0,epsiron = 0.000001d0
        real(8) rhoL,rhoR,rhoBar,UL,u_l,UR,u_r,VL,v_l,VR,v_r,pL,pR,HL,HR,HBar,cL,cR,cBar,cStar,cHalf,MP,MM,MHat,VnL,VnR,MBar,&
                absVnBar,absVnBar_M,absVnBar_P,chi,g,beta_P,beta_M,massFlow,Ptilde,thetaSD,chiSD,p1,p2,p3,p4,p5,p6,dPmax,PBar,&
                fluxLmtP,fluxlmtM,beta,deltaP,deltaM,hoge,fuga,f_(3),s_(3),IS(3),bWCNS(3),omega(3),qtilde(3),h,qtmp,&
                nx,ny,kx,ky,nmag,kmag,F1,F2,F3,F4

        beta = b!0.5d0*(3.0d0 - phi)/(1.0d0 - phi)
        h = 1.0d0/Nxnum

        nx = n(i,j,1)
        ny = n(i,j,2)
        nmag = sqrt(nx**2.0d0 + ny**2.0d0)
        kx = ny
        ky = -nx
        kmag = sqrt(kx**2.0d0 + ky**2.0d0)

        nx = nx/nmag!mx hat
        ny = ny/nmag!my hat
        kx = kx/kmag!kx hat
        ky = ky/kmag!ky hat

        !!!!!!!!! WCNS !!!!!!!
        !do k =1,4
        !    f_(1) = 0.5d0*(Q(i,j-2,k) - 4.0d0*Q(i,j-1,k) + 3.0d0*Q(i,j,k))
        !    f_(2) = 0.5d0*(Q(i,j+1,k) - Q(i,j-1,k))
        !    f_(3) = 0.5d0*(-3.0d0*Q(i,j,k) + 4.0d0*Q(i,j+1,k) - Q(i,j+2,k))
        !    s_(1) = (Q(i,j-2,k) - 2.0d0*Q(i,j-1,k) + Q(i,j,k))
        !    s_(2) = (Q(i,j-1,k) - 2.0d0*Q(i,j,k)   + Q(i,j+1,k))
        !    s_(3) = (Q(i,j,k)   - 2.0d0*Q(i,j+1,k) + Q(i,j+2,k))
        !    IS(1) = f_(1)**2.0d0 + s_(1)**2.0d0
        !    IS(2) = f_(2)**2.0d0 + s_(2)**2.0d0
        !    IS(3) = f_(3)**2.0d0 + s_(3)**2.0d0
        !    bWCNS(1) = CL1/((epsiron + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((epsiron + IS(2))**2.0d0)
        !    bWCNS(3) = CL3/((epsiron + IS(3))**2.0d0)
        !    omega(1) = bWCNS(1)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(2) = bWCNS(2)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(3) = bWCNS(3)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    qtilde(1) = Q(i,j,k) + 0.5d0*f_(1) + 0.125d0*s_(1)
        !    qtilde(2) = Q(i,j,k) + 0.5d0*f_(2) + 0.125d0*s_(2)
        !    qtilde(3) = Q(i,j,k) + 0.5d0*f_(3) + 0.125d0*s_(3)
        !    if(k == 1)then
        !        rhoL = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !    else if(k == 2)then
        !        u_l = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoL
        !    else if(k == 3)then
        !        v_l = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoL
        !    else if(k == 4)then
        !        qtmp = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !        pL = (GAMMA - 1.0d0)*(qtmp - 0.5d0*rhoL*(u_l**2.0d0 + v_l**2.0d0))
        !    endif
        !    f_(1) = 0.5d0*(Q(i,j-1,k) - 4.0d0*Q(i,j,k) + 3.0d0*Q(i,j+1,k))
        !    f_(2) = 0.5d0*(Q(i,j+2,k) - Q(i,j,k))
        !    f_(3) = 0.5d0*(-3.0d0*Q(i,j+1,k) + 4.0d0*Q(i,j+2,k) - Q(i,j+3,k))
        !    s_(1) = (Q(i,j-1,k) - 2.0d0*Q(i,j,k)   + Q(i,j+1,k))
        !    s_(2) = (Q(i,j,k)   - 2.0d0*Q(i,j+1,k) + Q(i,j+2,k))
        !    s_(3) = (Q(i,j+1,k) - 2.0d0*Q(i,j+2,k) + Q(i,j+3,k))
        !    IS(1) = f_(1)**2.0d0 + s_(1)**2.0d0
        !    IS(2) = f_(2)**2.0d0 + s_(2)**2.0d0
        !    IS(3) = f_(3)**2.0d0 + s_(3)**2.0d0
        !    bWCNS(1) = CL3/((epsiron + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((epsiron + IS(2))**2.0d0)
        !    bWCNS(3) = CL1/((epsiron + IS(3))**2.0d0)
        !    omega(1) = bWCNS(1)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(2) = bWCNS(2)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    omega(3) = bWCNS(3)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !    qtilde(1) = Q(i,j+1,k) - 0.5d0*f_(1) + 0.125d0*s_(1)
        !    qtilde(2) = Q(i,j+1,k) - 0.5d0*f_(2) + 0.125d0*s_(2)
        !    qtilde(3) = Q(i,j+1,k) - 0.5d0*f_(3) + 0.125d0*s_(3)
        !    if(k == 1)then
        !        rhoR = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !    else if(k == 2)then
        !        u_r = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoR
        !    else if(k == 3)then
        !        v_r = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoR
        !    else if(k == 4)then
        !        qtmp = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !        pR = (GAMMA - 1.0d0)*(qtmp - 0.5d0*rhoR*(u_r**2.0d0 + v_r**2.0d0))
        !    endif
        !enddo

        !!!!!!!!!! MUSCL
        deltaP = Q(i,j+1,1) - Q(i,j,1)
        deltaM = Q(i,j,1) - Q(i,j-1,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)
        rhoL = Q(i,j,1) + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
        !
        deltaP = Q(i,j+2,1) - Q(i,j+1,1)
        deltaM = Q(i,j+1,1) - Q(i,j,1)
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)
        rhoR = Q(i,j+1,1) - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
        !
        u_l   = Q(i,j,2)/Q(i,j,1)
        deltaP = Q(i,j+1,2)/Q(i,j+1,1) - Q(i,j,2)/Q(i,j,1)
        deltaM = Q(i,j,2)/Q(i,j,1) - Q(i,j-1,2)/Q(i,j-1,1)
        fluxLmtM = minmod(deltaP,beta*deltaM)
        fluxLmtP  = minmod(deltaM,beta*deltaP)
        u_l = u_l + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP     
        
        !
        v_l   =  Q(i,j,3)/Q(i,j,1)
        deltaP = Q(i,j+1,3)/Q(i,j+1,1) - Q(i,j,3)/Q(i,j,1)
        deltaM = Q(i,j,3)/Q(i,j,1) - Q(i,j-1,3)/Q(i,j-1,1)
        fluxLmtM = minmod(deltaP,beta*deltaM)
        fluxLmtP  = minmod(deltaM,beta*deltaP)
        v_l = v_l + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
        !
        u_r   = Q(i,j+1,2)/Q(i,j+1,1)
        deltaP = Q(i,j+2,2)/Q(i,j+2,1) - Q(i,j+1,2)/Q(i,j+1,1)
        deltaM = Q(i,j+1,2)/Q(i,j+1,1) - Q(i,j,2)/Q(i,j,1)
        fluxLmtM = minmod(deltaP,beta*deltaM)
        fluxLmtP  = minmod(deltaM,beta*deltaP)
        u_r = u_r - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
        !     
        v_r   = Q(i,j+1,3)/Q(i,j+1,1)
        deltaP = Q(i,j+2,3)/Q(i,j+2,1) - Q(i,j+1,3)/Q(i,j+1,1)
        deltaM = Q(i,j+1,3)/Q(i,j+1,1) - Q(i,j,3)/Q(i,j,1)
        fluxLmtM = minmod(deltaP,beta*deltaM)
        fluxLmtP  = minmod(deltaM,beta*deltaP)
        v_r = v_r - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
        !
        hoge   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
        pL     = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        fuga   = (GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1))
        deltaP = hoge - pL
        deltaM = pL - fuga
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)
        pL   = pL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP 
        !
        hoge   = (GAMMA-1.0d0)*(Q(i,j+2,4) - 0.5d0*(Q(i,j+2,2)**2.0d0 + Q(i,j+2,3)**2.0d0)/Q(i,j+2,1))
        pR     = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
        fuga   = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        deltaP = hoge - pR
        deltaM = pR - fuga
        fluxLmtP = minmod(deltaP,beta*deltaM)
        fluxLmtM = minmod(deltaM,beta*deltaP)          
        pR = pR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
                
        rhoBar = 0.5d0*(rhoL + rhoR)
        UL = u_l*kx + v_l*ky
        UR = u_r*kx + v_r*ky
        VL = u_l*nx + v_l*ny
        VR = u_r*nx + v_r*ny
        VnL  = VL
        VnR  = VR
        PBar = 0.5d0*(pL + pR)

        HL   = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0))/rhoL
        HR   = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0))/rhoR
        !HBar = 0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0) + (GAMMA/(GAMMA -1.0d0))*PBar/rhoBar
        HBar = 0.5d0*(HL+HR)

        cStar = sqrt(2.0d0*(GAMMA - 1.0d0)*HBar/(GAMMA + 1.0d0))  
        cL   = sqrt(GAMMA*pL/rhoL)             !SLAU
        cR   = sqrt(GAMMA*pR/rhoR)             !SLAU
        !cL   = (cStar**2.0d0)/max(cStar,abs(VnL)) !SLAU2
        !cR   = (cStar**2.0d0)/max(cStar,abs(VnR)) !SLAU2   
        cBar = 0.5d0*(cL + cR)
        cHalf= min(cL,cR)
        
        MP   = VnL/cBar
        MM   = VnR/cBar

        ! mass flux parameter for SLAU / SLAU2
        absVnBar  = (rhoL*abs(VnL) + rhoR*abs(VnR))/(rhoL + rhoR)
        MHat  = min(1.0d0,sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))/cBar)
        chi = (1.0d0 - MHat)**2.0d0                                                      !SLAU / SLAU2
        g = - max(min(MP,0.0d0),-1.0d0)*min(max(MM,0.0d0),1.0d0)

        ! mass flux parameter for SD-SLAU
        p1 = abs((GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1)) - pL)
        p2 = abs((GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1)) - pL)
        p3 = abs((GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1)) - pL)
        p4 = abs((GAMMA-1.0d0)*(Q(i+1,j+1,4) - 0.5d0*(Q(i+1,j+1,2)**2.0d0 + Q(i+1,j+1,3)**2.0d0)/Q(i+1,j+1,1)) - pR)
        p5 = abs((GAMMA-1.0d0)*(Q(i,j+2,4) - 0.5d0*(Q(i,j+2,2)**2.0d0 + Q(i,j+2,3)**2.0d0)/Q(i,j+2,1)) - pR)
        p6 = abs((GAMMA-1.0d0)*(Q(i-1,j+1,4) - 0.5d0*(Q(i-1,j+1,2)**2.0d0 + Q(i-1,j+1,3)**2.0d0)/Q(i-1,j+1,1)) - pR)
        dPmax = max(p1,p2,p3,p4,p5,p6)
        MBar = ((rhoL*abs(VL) + rhoR*abs(VR))/(rhoL + rhoR))/cBar
        thetaSD = min(1.0d0,((Csd2*(pR-pL) + Csd1*PBar)/(dPmax + Csd1*PBar))**2.0d0)*max(0.0d0,1.0d0 - abs(MBar))
        chiSD   = 0.5d0*thetaSD*(abs(MBar + 1.0d0) + abs(MBar - 1.0d0) - 2.0d0*abs(MBar))!SD-SLAU


        massFlow = 0.5d0*(rhoL*VnL + rhoR*VnR - absVnBar*(rhoR - rhoL))*(1.0d0 - g) - 0.5d0*chiSD*(pR - pL)/cBar
        !massFlow = 0.5d0*(rhoL*VnL + rhoR*VnR - absVnBar*(rhoR - rhoL))*(1.0d0 - g) - 0.5d0*chi*(pR - pL)/cBar

        if(1.0d0 <= abs(MP))then
            beta_P = 0.5d0*(1.0d0 + sign(1.0d0,MP))
        else
            beta_P = 0.25d0*(2.0d0 - MP)*(MP + 1.0d0)**2.0d0
        endif

        if(1.0d0 <= abs(MM))then
            beta_M = 0.5d0*(1.0d0 + sign(1.0d0,-MM))
        else
            beta_M = 0.25d0*(2.0d0 + MM)*(MM - 1.0d0)**2.0d0
        endif

        !pressure flux parameter for SLAU2

        !absVnBar_P = (1.0d0 - g)*absVnBar + g*abs(VnL)
        !absVnBar_M = (1.0d0 - g)*absVnBar + g*abs(VnR)
        !cHalf = min(cL,cR)

        Ptilde = 0.5d0*(pL + pR) + 0.5d0*(beta_P - beta_M)*(pL - pR)&
                    !+ (beta_P + beta_M - 1.0d0)*sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))*rhoBar*cHalf!SLAU2
                    + 0.5d0*(1.0d0 - chi)*(beta_P + beta_M - 1.0d0)*(pL + pR)                   !SLAU

        F1 = massFlow
        F2 = 0.5d0*(massFlow*(UR + UL) - abs(massFlow)*(UR - UL))
        F3 = 0.5d0*(massFlow*(VR + VL) - abs(massFlow)*(VR - VL)) + Ptilde
        F4 = 0.5d0*(massFlow*(HR + HL) - abs(massFlow)*(HR - HL))
        !write(*,*)i,j,'F',F1,F2,F3,F4
        
        F(i,j,1) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*F1!*(1.0d0 - is_SF_eta(i,j))
        F(i,j,2) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*(kx*F2 + nx*F3)!*(1.0d0 - is_SF_eta(i,j))
        F(i,j,3) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*(ky*F2 + ny*F3)!*(1.0d0 - is_SF_eta(i,j))
        F(i,j,4) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*F4!*(1.0d0 - is_SF_eta(i,j))
        !if(i<50 .and. j==1)then
        !    write(*,*)i,rhoL,u_l,v_l,pL
        !endif
        
    end subroutine SD_SLAU_YFlux

end module SD_SLAUFlux