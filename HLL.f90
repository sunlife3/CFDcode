module HLLFlux
    use param
    implicit none
contains
    subroutine HLL_XFlux(Q,E,m,Nx,Ny,i,j,b)
        implicit none
        integer,intent(in) :: i,j,Nx,Ny
        integer :: n=0,foo=0,k=1
        real(8),intent(in) :: Q(-2:Nx+3,-2:Ny+3,4),m(-2:Nx+3,-2:Ny+3,2),b
        real(8),intent(out) :: E(-1:Nx+2,-1:Ny+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0, eps = 0.000001d0, CL1=1.0d0/16.0d0,CL2=10.0d0/16.0d0,CL3=5.0d0/16.0d0
        real(8) rhoL,rhoR,rhoBar,UL,u_l,UR,u_r,VL,v_l,VR,v_r,pL,pR,HL,HR,HBar,cL,cR,cBar,ubar,MHat,VnL,VnR,MBar,&
                massFlow,PBar,fluxLmtP,fluxlmtM,beta,deltaP,deltaM,hoge,fuga,eL,eR,&
                mx,my,kx,ky,mmag,kmag,E1,E2,E3,E4,SL,SR,Sast,hast,past,Uast(4),phiL,phiR,epsiron(4),&
                f_(3),s_(3),IS(3),bWCNS(3),omega(3),qtilde(3),h,qtmp
        n=n+1
        
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

        !!!!!!!! WCNS !!!!!!!
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
        !    bWCNS(1) = CL1/((eps + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((eps + IS(2))**2.0d0)
        !    bWCNS(3) = CL3/((eps + IS(3))**2.0d0)
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
        !    bWCNS(1) = CL3/((eps + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((eps + IS(2))**2.0d0)
        !    bWCNS(3) = CL1/((eps + IS(3))**2.0d0)
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

        !!!!! MUSCL !!!
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
        UL = u_l*mx + v_l*my
        UR = u_r*mx + v_r*my
        VL = u_l*kx + v_l*ky
        VR = u_r*kx + v_r*ky

        VnL  = UL
        VnR  = UR
        HL   = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0))/rhoL
        HR   = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0))/rhoR
        cL   = sqrt(GAMMA*pL/rhoL)
        cR   = sqrt(GAMMA*pR/rhoR)
        cBar = sqrt((sqrt(rhoL)*cL**2.0d0 + sqrt(rhoR)*cR**2.0d0)/(sqrt(rhoL)+sqrt(rhoR))) !&
                !+ 0.5d0*(sqrt(rhoL)*sqrt(rhoR))*(UR - UL)**2.0d0/((sqrt(rhoL) + sqrt(rhoR))**2.0d0))
        ubar = (sqrt(rhoL)*UL + sqrt(rhoR)*UR)/(sqrt(rhoL) + sqrt(rhoR))
        
        eL = pL/(1.4d0 - 1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0)
        eR = pR/(1.4d0 - 1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0)
        
        massFlow = 0.5d0*(rhoL*UL + rhoR*UR)
        
        !! HLL
        SL = min(UL-cL*sqrt(mx**2.0d0 + my**2.0d0),UR-cR*sqrt(mx**2.0d0 + my**2.0d0))
        SR = max(UL+cL*sqrt(mx**2.0d0 + my**2.0d0),UR+cR*sqrt(mx**2.0d0 + my**2.0d0))
        
        ! HLLC
        !SL = min(UL - cL*sqrt(nx**2.0d0 + ny**2.0d0),ubar - cbar*sqrt(nx**2.0d0 + ny**2.0d0))
        !SR = max(UR + cR*sqrt(nx**2.0d0 + ny**2.0d0),ubar + cbar*sqrt(nx**2.0d0 + ny**2.0d0))
        Sast = (pR - pL + rhoL*UL*(SL - UL) - rhoR*UR*(SR - UR))/(rhoL*(SL - UL) - rhoR*(SR - UR))

        !!!!!!!!!!!!!!!!!!!! HLL Flux !!!!!!!!!!!!!!!!!!!!!!!!!!! 
        if(0.0d0 < SL)then
            E1 = rhoL*UL
            E2 = rhoL*UL**2.0d0 + pL
            E3 = rhoL*UL*VL
            E4 = (eL + pL)*UL
        else if(SL <= 0.0d0 .and. 0.0d0 <= SR)then
            E1 = (SR*rhoL*UL - SL*rhoR*UR + SL*SR*(rhoR - rhoL))/(SR-SL)
            E2 = (SR*(rhoL*UL**2.0d0 + pL) - SL*(rhoR*UR**2.0d0 + pR) + SL*SR*(rhoR*UR - rhoL*UL))/(SR-SL)
            E3 = (SR*(rhoL*UL*VL) - SL*(rhoR*UR*VR) + SL*SR*(rhoR*VR - rhoL*VL))/(SR-SL)
            E4 = (SR*(eL+pL)*UL - SL*(eR+pR)*UR + SL*SR*(eR - eL))/(SR-SL)
        else
            E1 = rhoR*UR
            E2 = rhoR*UR**2.0d0 + pR
            E3 = rhoR*UR*VR
            E4 = (eR + pR)*UR
        endif

        !!!!!!!!!!!!!!!!!!!!!! HLLC Flux !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !if(0.0d0 < SL)then
        !    E1 = rhoL*UL
        !    E2 = rhoL*UL**2.0d0 + pL
        !    E3 = rhoL*UL*VL
        !    E4 = (eL + pL)*UL
        !else if(SL <= 0.0d0 .and. 0.0d0 < Sast)then
        !    past = pL + rhoL*(SL - UL)*(Sast - UL)
        !    Uast(1) = (SL - UL)*rhoL/(SL - Sast)
        !    Uast(2) = (SL - UL)*rhoL*Sast/(SL - Sast)
        !    Uast(3) = (SL - UL)*rhoL*VL/(SL - Sast)
        !    !Uast(4) = (SL - UL)*(eL + (Sast - UL)*(rhoL*Sast + pL/(SL - UL)))/(SL - Sast)
        !    Uast(4) = (SL - UL)*eL/(SL - Sast)
        !    E1 = rhoL*UL + SL*(Uast(1) - rhoL)
        !    E2 = rhoL*UL**2.0d0 + pL + SL*(Uast(2) - rhoL*UL)
        !    E3 = rhoL*UL*VL + SL*(Uast(3) - rhoL*VL)
        !    E4 = (eL + pL)*UL + SL*(Uast(4) - eL)
        !else if(Sast <= 0.0d0 .and. 0.0d0 < SR)then
        !    past = pR + rhoR*(SR - UR)*(Sast - UR)
        !    Uast(1) = (SR - UR)*rhoR/(SR - Sast)
        !    Uast(2) = (SR - UR)*rhoR*Sast/(SR - Sast)
        !    Uast(3) = (SR - UR)*rhoR*VR/(SR - Sast)
        !    !Uast(4) = (SR - UR)*(eR + (Sast - UR)*(rhoR*Sast + pR/(SR - UR)))/(SR - Sast)
        !    Uast(4) = (SR - UR)*eR/(SR - Sast)
        !    E1 = rhoR*UR + SR*(Uast(1) - rhoR)
        !    E2 = rhoR*UR**2.0d0 + pR + SR*(Uast(2) - rhoR*UR)
        !    E3 = rhoR*UR*VR + SR*(Uast(3) - rhoR*VR)
        !    E4 = (eR + pR)*UR + SR*(Uast(4) - eR)    
        !else if(SR <= 0.0d0)then
        !    E1 = rhoR*UR
        !    E2 = rhoR*UR**2.0d0 + pR
        !    E3 = rhoR*UR*VR
        !    E4 = (eR + pR)*UR
        !endif
        !!       
        !epsiron(1) = 0.5d0*(SR + SL)/(SR - SL)*(rhoL*UL - rhoR*UR) - &
        !                 SL*SR/(SR - SL)*(rhoL - rhoR)
        !epsiron(2) = 0.5d0*(SR + SL)/(SR - SL)*((rhoL*UL**2.0d0 + pL) - (rhoR*UR**2.0d0 + pR)) - &
        !                 SL*SR/(SR - SL)*(rhoL*UL - rhoR*UR)
        !epsiron(3) = 0.5d0*(SR + SL)/(SR - SL)*(rhoL*UL*VL - rhoR*UR*VR) - &
        !                 SL*SR/(SR - SL)*(rhoL*VL - rhoR*VR)
        !epsiron(4) = 0.5d0*(SR + SL)/(SR - SL)*((UL*(eL + pL)) - (UR*(eR + pR))) - &
        !                 SL*SR/(SR - SL)*(eL - eR) ! HLL Riemann Solver
        
        !E1 = massFlow + epsiron(1)
        !E2 = 0.5d0*((rhoL*UL**2.0d0 + pL) + (rhoR*UR**2.0d0 + pR)) + epsiron(2)
        !E3 = 0.5d0*(rhoL*UL*VL + rhoR*UR*VR) + epsiron(3)
        !E4 = 0.5d0*(rhoL*UL*HL + rhoR*UR*HR) + epsiron(4)

        E(i,j,1) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*E1
        E(i,j,2) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*(mx*E2 + kx*E3)
        E(i,j,3) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*(my*E2 + ky*E3)
        E(i,j,4) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*E4

        !write(*,*)i,j,E(i,j,1),E(i,j,2),E(i,j,3),E(i,j,4)

    end subroutine HLL_XFlux

    subroutine HLL_YFlux(Q,F,n,Nxnum,Nynum,i,j,b)
        implicit none
        integer,intent(in) :: i,j,Nxnum,Nynum
        integer :: iter = 0, k=0
        real(8),intent(in) :: Q(-2:Nxnum+3,-2:Nynum+3,4),n(-2:Nxnum+3,-2:Nynum+3,2),b
        real(8),intent(out) :: F(-1:Nxnum+2,-1:Nynum+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0,CL1=1.0d0/16.0d0,CL2=10.0d0/16.0d0,CL3=5.0d0/16.0d0,eps=0.000001d0
        real(8) rhoL,rhoR,rhoBar,UL,u_l,UR,u_r,VL,v_l,VR,v_r,pL,pR,HL,HR,HBar,cL,cR,cBar,ubar,VnL,VnR,&
                massFlow,fluxLmtP,fluxlmtM,beta,deltaP,deltaM,hoge,fuga,eL,eR,&
                nx,ny,kx,ky,nmag,kmag,F1,F2,F3,F4,SL,SR,Sast,past,Uast(4),epsiron(4),&
                f_(3),s_(3),IS(3),bWCNS(3),omega(3),qtilde(3),h,qtmp

        beta = b!0.5d0*(3.0d0 - phi)/(1.0d0 - phi)

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

        !!!!!!!!! WCNS
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
        !    bWCNS(1) = CL1/((eps + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((eps + IS(2))**2.0d0)
        !    bWCNS(3) = CL3/((eps + IS(3))**2.0d0)
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
        !    bWCNS(1) = CL3/((eps + IS(1))**2.0d0)
        !    bWCNS(2) = CL2/((eps + IS(2))**2.0d0)
        !    bWCNS(3) = CL1/((eps + IS(3))**2.0d0)
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

        !!!!! MUSCL
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
        HL   = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0))/rhoL
        HR   = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0))/rhoR
        HBar = 0.5d0*(HL+HR)
        eL = pL/(1.4d0 - 1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0)
        eR = pR/(1.4d0 - 1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0)

        cL   = sqrt(GAMMA*pL/rhoL)
        cR   = sqrt(GAMMA*pR/rhoR)
        cBar = sqrt((sqrt(rhoL)*cL**2.0d0 + sqrt(rhoR)*cR**2.0d0)/(sqrt(rhoL) + sqrt(rhoR))) !&
                !+ 0.5d0*(sqrt(rhoL)*sqrt(rhoR))*(VR - VL)**2.0d0/((sqrt(rhoL) + sqrt(rhoR))**2.0d0))
        ubar = (sqrt(rhoL)*VL + sqrt(rhoR)*VR)/(sqrt(rhoL) + sqrt(rhoR))
    
        massFlow = 0.5d0*(rhoL*VL + rhoR*VR)
        
        !! HLL
        SL = min(VL - cL*sqrt(nx**2.0d0 + ny**2.0d0),VR - cR*sqrt(nx**2.0d0 + ny**2.0d0))
        SR = max(VL + cL*sqrt(nx**2.0d0 + ny**2.0d0),VR + cR*sqrt(nx**2.0d0 + ny**2.0d0))
        
        ! HLLC
        !SL = min(VL - cL*sqrt(nx**2.0d0 + ny**2.0d0), ubar - cbar*sqrt(nx**2.0d0 + ny**2.0d0))
        !SR = max(VR + cR*sqrt(nx**2.0d0 + ny**2.0d0), ubar + cbar*sqrt(nx**2.0d0 + ny**2.0d0))
        Sast = (pR - pL + rhoL*VL*(SL - VL) - rhoR*VR*(SR - VR))/(rhoL*(SL - VL) - rhoR*(SR - VR))
        !write(*,*)i,j,SL,SR,Sast

        !!!!!!!!!!!! HLL !!!!!!!!!!!!!!!!!!!
        if(0.0d0 < SL)then
            F1 = rhoL*VL
            F2 = rhoL*VL*UL
            F3 = rhoL*VL**2.0d0 + pL
            F4 = (eL + pL)*VL
        else if(SL <= 0.0d0 .and. 0.0d0 <= SR)then
            F1 = (SR*rhoL*VL - SL*rhoR*VR + SL*SR*(rhoR - rhoL))/(SR-SL)
            F2 = (SR*(rhoL*VL*UL) - SL*(rhoR*VR*UR) + SL*SR*(rhoR*UR - rhoL*UL))/(SR-SL)
            F3 = (SR*(rhoL*VL**2.0d0 + pL) - SL*(rhoR*VR**2.0d0 + pR) + SL*SR*(rhoR*VR - rhoL*VL))/(SR-SL)
            F4 = (SR*(eL+pL)*VL - SL*(eR+pR)*VR + SL*SR*(eR - eL))/(SR-SL)
        else
            F1 = rhoR*VR
            F2 = rhoR*VR*UR
            F3 = rhoR*VR**2.0d0 + pR
            F4 = (eR + pR)*VR
        endif

        !!!!!!!!!!!!!! HLLC !!!!!!!!!!!!!!!!!!!
        !if(0.0d0 <= SL)then
        !    F1 = rhoL*VL
        !    F2 = rhoL*VL*UL
        !    F3 = rhoL*VL**2.0d0 + pL
        !    F4 = (eL + pL)*VL
        !else if(SL < 0.0d0 .and. 0.0d0 <= Sast)then
        !    past = pL + rhoL*(SL - VL)*(Sast - VL)
        !    Uast(1) = (SL - VL)*rhoL/(SL - Sast)
        !    Uast(2) = (SL - VL)*rhoL*UL/(SL - Sast)
        !    Uast(3) = (SL - VL)*rhoL*Sast/(SL - Sast)
        !    !Uast(4) = (SL - VL)*(eL + (Sast - VL)*(rhoL*Sast + pL/(SL - VL)))/(SL - Sast)
        !    Uast(4) = (SL - VL)*eL/(SL - Sast)
        !    F1 = rhoL*VL + SL*(Uast(1) - rhoL)
        !    F2 = rhoL*VL*UL + SL*(Uast(2) - rhoL*UL)
        !    F3 = rhoL*VL**2.0d0 + pL + SL*(Uast(3) - rhoL*VL)
        !    F4 = (eL + pL)*VL + SL*(Uast(4) - eL)
        !else if(Sast < 0.0d0 .and. 0.0d0 <= SR)then
        !    past = pR + rhoR*(SR - VR)*(Sast - VR)
        !    Uast(1) = (SR - VR)*rhoR/(SR - Sast)
        !    Uast(2) = (SR - VR)*rhoR*UR/(SR - Sast)
        !    Uast(3) = (SR - VR)*rhoR*Sast/(SR - Sast)
        !    !Uast(4) = (SR - VR)*(eR + (Sast - VR)*(rhoR*Sast + pR/(SR - VR)))/(SR - Sast)
        !    Uast(4) = (SR - VR)*eR/(SR - Sast)
        !    F1 = rhoR*VR + SR*(Uast(1) - rhoR)
        !    F2 = rhoR*VR*UR + SR*(Uast(2) - rhoR*UR)
        !    F3 = rhoR*VR**2.0d0 + pR + SR*(Uast(3) - rhoR*VR)
        !    F4 = (eR + pR)*VR + SR*(Uast(4) - eR)
        !else if(SR < 0.0d0)then
        !    F1 = rhoR*VR
        !    F2 = rhoR*VR*UR
        !    F3 = rhoR*VR**2.0d0 + pR
        !    F4 = (eR + pR)*VR
        !endif
        !!
        !epsiron(1) = 0.5d0*(SR + SL)/(SR - SL)*(rhoL*VL - rhoR*VR) - &
        !                 SL*SR/(SR - SL)*(rhoL - rhoR)
        !epsiron(2) = 0.5d0*(SR + SL)/(SR - SL)*((rhoL*VL*UL) - (rhoR*VR*UR)) -&
        !                 SL*SR/(SR - SL)*(rhoL*UL - rhoR*UR)
        !epsiron(3) = 0.5d0*(SR + SL)/(SR - SL)*((rhoL*VL**2.0d0 + pL) - (rhoR*VR**2.0d0 + pR)) -&
        !                 SL*SR/(SR - SL)*(rhoL*VL - rhoR*VR)
        !epsiron(4) = 0.5d0*(SR + SL)/(SR - SL)*((VL*(eL + pL)) - (VR*(eR + pR))) - &
        !                 SL*SR/(SR - SL)*(eL - eR) ! HLL Riemann Solver
        !F1 = massFlow + epsiron(1)
        !F2 = 0.5d0*(rhoL*VL*UL + rhoR*VR*UR) + epsiron(2)
        !F3 = 0.5d0*((rhoL*VL**2.0d0 + pL) + (rhoR*VR**2.0d0 + pR)) + epsiron(3)
        !F4 = 0.5d0*(rhoL*VL*HL + rhoR*VR*HR) + epsiron(4)
        !write(*,*)i,j,'F',F1,F2,F3,F4
        
        F(i,j,1) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*F1 
        F(i,j,2) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*(kx*F2 + nx*F3) 
        F(i,j,3) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*(ky*F2 + ny*F3) 
        F(i,j,4) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*F4 
        
        
    end subroutine HLL_YFlux

end module HLLFlux