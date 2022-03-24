module AUSMDVFlux
    use param
    implicit none
contains
    subroutine AusmDV_XFlux(Q,E,m,Nx,Ny,i,j,b)
        implicit none
        integer,intent(in) :: i,j,Nx,Ny
        integer :: k=1,l=0
        real(8),intent(in) :: Q(-2:Nx+3,-2:Ny+3,4),m(-2:Nx+3,-2:Ny+3,2),b
        real(8),intent(out) :: E(-1:Nx+2,-1:Ny+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0,largeK = 10.0d0,&
                            CL1=1.0d0/16.0d0,CL2=10.0d0/16.0d0,CL3=5.0d0/16.0d0,epsiron= 0.000001d0
        real(8) :: beta=0,fluxLimMin=0,fluxLimPl=0, &
                rhoL=0,u_l=0,v_l=0,UL=0,VL=0,pL=0,rhoR=0,u_r=0,v_r=0,UR=0,VR=0,pR=0,&
                pLp=0,pRm=0,ULp=0,URm=0,alphaL=0,alphaR=0,HL=0,HR=0,cm=0,s=0,&
                massFlow=0,momentumV=0,momentumD=0,&
                hoge=0,fuga=0,tmp1=0,tmp2=0,f_(3),s_(3),IS(3),bWCNS(3),omega(3),qtilde(3),h,qtmp,&
                mx=0,my=0,kx=0,ky=0,mmag=0,kmag=0,&
                E1=0,E2=0,E3=0,E4=0
        
        !beta = 3.0d0!(1.0d0 + (3.0d0-phi)/(1.0d0-phi))
        beta = b
        l = l+1
        !h = X(i+1,j)-X(i,j)!1.0d0/Nx

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
        
        !!!!!! WCNS
        !do k =1,4   
        !   f_(1) = 0.5d0*(Q(i-2,j,k) - 4.0d0*Q(i-1,j,k) + 3.0d0*Q(i,j,k))
        !   f_(2) = 0.5d0*(Q(i+1,j,k) - Q(i-1,j,k))
        !   f_(3) = 0.5d0*(-3.0d0*Q(i,j,k) + 4.0d0*Q(i+1,j,k) - Q(i+2,j,k))
        !   s_(1) = (Q(i-2,j,k) - 2.0d0*Q(i-1,j,k) + Q(i,j,k))
        !   s_(2) = (Q(i-1,j,k) - 2.0d0*Q(i,j,k)   + Q(i+1,j,k))
        !   s_(3) = (Q(i,j,k)   - 2.0d0*Q(i+1,j,k) + Q(i+2,j,k))
        !   IS(1) = f_(1)**2.0d0 + s_(1)**2.0d0
        !   IS(2) = f_(2)**2.0d0 + s_(2)**2.0d0
        !   IS(3) = f_(3)**2.0d0 + s_(3)**2.0d0
        !   bWCNS(1) = CL1/((epsiron + IS(1))**2.0d0)
        !   bWCNS(2) = CL2/((epsiron + IS(2))**2.0d0)
        !   bWCNS(3) = CL3/((epsiron + IS(3))**2.0d0)
        !   omega(1) = bWCNS(1)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !   omega(2) = bWCNS(2)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !   omega(3) = bWCNS(3)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !   qtilde(1) = Q(i,j,k) + 0.5d0*f_(1) + 0.125d0*s_(1)
        !   qtilde(2) = Q(i,j,k) + 0.5d0*f_(2) + 0.125d0*s_(2)
        !   qtilde(3) = Q(i,j,k) + 0.5d0*f_(3) + 0.125d0*s_(3)
        !   if(k == 1)then
        !       rhoL = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !   else if(k == 2)then
        !       u_l = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoL
        !   else if(k == 3)then
        !       v_l = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoL
        !   else if(k == 4)then
        !       qtmp = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !       pL = (GAMMA - 1.0d0)*(qtmp - 0.5d0*rhoL*(u_l**2.0d0 + v_l**2.0d0))
        !   endif
        !   f_(1) = 0.5d0*(Q(i-1,j,k) - 4.0d0*Q(i,j,k) + 3.0d0*Q(i+1,j,k))
        !   f_(2) = 0.5d0*(Q(i+2,j,k) - Q(i,j,k))
        !   f_(3) = 0.5d0*(-3.0d0*Q(i+1,j,k) + 4.0d0*Q(i+2,j,k) - Q(i+3,j,k))
        !   s_(1) = (Q(i-1,j,k) - 2.0d0*Q(i,j,k)   + Q(i+1,j,k))
        !   s_(2) = (Q(i,j,k)   - 2.0d0*Q(i+1,j,k) + Q(i+2,j,k))
        !   s_(3) = (Q(i+1,j,k) - 2.0d0*Q(i+2,j,k) + Q(i+3,j,k))
        !   IS(1) = f_(1)**2.0d0 + s_(1)**2.0d0
        !   IS(2) = f_(2)**2.0d0 + s_(2)**2.0d0
        !   IS(3) = f_(3)**2.0d0 + s_(3)**2.0d0
        !   bWCNS(1) = CL3/((epsiron + IS(1))**2.0d0)
        !   bWCNS(2) = CL2/((epsiron + IS(2))**2.0d0)
        !   bWCNS(3) = CL1/((epsiron + IS(3))**2.0d0)
        !   omega(1) = bWCNS(1)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !   omega(2) = bWCNS(2)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !   omega(3) = bWCNS(3)/(bWCNS(1) + bWCNS(2) + bWCNS(3))
        !   qtilde(1) = Q(i+1,j,k) - 0.5d0*f_(1) + 0.125d0*s_(1)
        !   qtilde(2) = Q(i+1,j,k) - 0.5d0*f_(2) + 0.125d0*s_(2)
        !   qtilde(3) = Q(i+1,j,k) - 0.5d0*f_(3) + 0.125d0*s_(3)
        !   if(k == 1)then
        !       rhoR = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !   else if(k == 2)then
        !       u_r = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoR
        !   else if(k == 3)then
        !       v_r = (qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3))/rhoR
        !   else if(k == 4)then
        !       qtmp = qtilde(1)*omega(1) + qtilde(2)*omega(2) + qtilde(3)*omega(3)
        !       pR = (GAMMA - 1.0d0)*(qtmp - 0.5d0*rhoR*(u_r**2.0d0 + v_r**2.0d0))
        !   endif
        !enddo

        !!!!!!Higher presision by MUSCL approach
        tmp1 = Q(i+1,j,1) - Q(i,j,1)
        tmp2 = Q(i,j,1) - Q(i-1,j,1)
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        rhoL = Q(i,j,1) + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
        !
        tmp1 = Q(i+2,j,1) - Q(i+1,j,1)
        tmp2 = Q(i+1,j,1) - Q(i,j,1)
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        rhoR = Q(i+1,j,1) - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        !
        u_l = Q(i,j,2)/Q(i,j,1)
        tmp1 = Q(i+1,j,2)/Q(i+1,j,1) - Q(i,j,2)/Q(i,j,1) 
        tmp2 = Q(i,j,2)/Q(i,j,1) - Q(i-1,j,2)/Q(i-1,j,1)
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)            
        u_l   = u_l + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
        !
        v_l = Q(i,j,3)/Q(i,j,1)
        tmp1 = Q(i+1,j,3)/Q(i+1,j,1) - Q(i,j,3)/Q(i,j,1)
        tmp2 = Q(i,j,3)/Q(i,j,1) - Q(i-1,j,3)/Q(i-1,j,1)
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        v_l   = v_l + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
        !
        u_r = Q(i+1,j,2)/Q(i+1,j,1)                
        tmp1 = Q(i+2,j,2)/Q(i+2,j,1) - Q(i+1,j,2)/Q(i+1,j,1)
        tmp2 = Q(i+1,j,2)/Q(i+1,j,1) - Q(i,j,2)/Q(i,j,1)
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)            
        u_r   = u_r - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        !
        v_r = Q(i+1,j,3)/Q(i+1,j,1)
        tmp1 = Q(i+2,j,3)/Q(i+2,j,1) - Q(i+1,j,3)/Q(i+1,j,1)
        tmp2 = Q(i+1,j,3)/Q(i+1,j,1) - Q(i,j,3)/Q(i,j,1)
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        v_r = v_r - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin    
        !
        hoge   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
        pL     = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        fuga   = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1))
        tmp1 = hoge - pL
        tmp2 = pL - fuga
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)            
        pL   = pL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
        !
        hoge   = (GAMMA-1.0d0)*(Q(i+2,j,4) - 0.5d0*(Q(i+2,j,2)**2.0d0 + Q(i+2,j,3)**2.0d0)/Q(i+2,j,1))
        pR     = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
        fuga   = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
        tmp1 = hoge - pR
        tmp2 = pR - fuga
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)                
        pR   = pR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        
        UL = u_l*mx + v_l*my
        UR = u_r*mx + v_r*my
        VL = u_l*kx + v_l*ky
        VR = u_r*kx + v_r*ky

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

    subroutine AusmDV_YFlux(Q,F,n,Nxnum,Nynum,i,j,b)
        implicit none
        integer,intent(in) :: i,j,Nxnum,Nynum
        integer :: k=1,l=0
        real(8),intent(in) :: Q(-2:Nxnum+3,-2:Nynum+3,4),n(-2:Nxnum+3,-2:Nynum+3,2),b
        real(8),intent(out) :: F(-1:Nxnum+2,-1:Nynum+2,4)
        real(8),parameter :: phi = 1.0/3.0,largeK = 10.0,&
                            CL1=1.0d0/16.0d0,CL2=10.0d0/16.0d0,CL3=5.0d0/16.0d0,epsiron= 0.000001d0
        real(8) :: beta=0,fluxLimMin=0,fluxLimPl=0, &
                rhoL=0,u_l=0,v_l=0,UL=0,VL=0,pL=0,rhoR=0,u_r=0,v_r=0,UR=0,VR=0,pR=0,&
                pLp=0,pRm=0,ULp=0,URm=0,alphaL=0,alphaR=0,HL=0,HR=0,cm=0,s=0,&
                massFlow=0,momentumV=0,momentumD=0,&
                hoge=0,fuga=0,tmp1=0,tmp2=0,f_(3),s_(3),IS(3),bWCNS(3),omega(3),qtilde(3),h,qtmp,&
                nx=0,ny=0,kx=0,ky=0,nmag=0,kmag=0,&
                F1=0,F2=0,F3=0,F4=0
        
        l=l+1
        beta = b!3.0d0!(1.0d0 + (3.0d0-phi)/(1.0d0-phi))
        !h = Y(i,j+1)-Y(i,j)!0.178d0/Nynum
      
        nx = n(i,j,1)
        ny = n(i,j,2)
        nmag = sqrt(nx**2.0d0+ny**2.0d0)
        kx = ny
        ky = -nx
        kmag = sqrt(kx**2.0d0+ky**2.0d0)
        
        !normalized metrics
        nx = nx/nmag
        ny = ny/nmag
        kx = kx/kmag
        ky = ky/kmag
        
        !!!!!!!!!! WCNS !!!!!!!
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
        
        !!!! Higher resolution by MUSCL approach
        tmp1 = Q(i,j+1,1) - Q(i,j,1)
        tmp2 = Q(i,j,1) - Q(i,j-1,1)
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        rhoL = Q(i,j,1) + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
        !
        tmp1 = Q(i,j+2,1) - Q(i,j+1,1)
        tmp2 = Q(i,j+1,1) - Q(i,j,1)
        fluxLimPl = minmod(tmp1,beta*tmp1)
        fluxLimMin = minmod(tmp2,beta*tmp2)
        rhoR = Q(i,j+1,1) - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        !
        u_l   = Q(i,j,2)/Q(i,j,1)
        tmp1 = Q(i,j+1,2)/Q(i,j+1,1) - Q(i,j,2)/Q(i,j,1)
        tmp2 = Q(i,j,2)/Q(i,j,1) - Q(i,j-1,2)/Q(i,j-1,1)
        fluxLimPl  = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        u_l = u_l + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl     
        !
        v_l   =  Q(i,j,3)/Q(i,j,1)
        tmp1 = Q(i,j+1,3)/Q(i,j+1,1) - Q(i,j,3)/Q(i,j,1)
        tmp2 = Q(i,j,3)/Q(i,j,1) - Q(i,j-1,3)/Q(i,j-1,1)
        fluxLimPl  = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        v_l = v_l + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl
        !
        u_r   = Q(i,j+1,2)/Q(i,j+1,1)
        tmp1 = Q(i,j+2,2)/Q(i,j+2,1) - Q(i,j+1,2)/Q(i,j+1,1)
        tmp2 = Q(i,j+1,2)/Q(i,j+1,1) - Q(i,j,2)/Q(i,j,1)
        fluxLimPl  = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        u_r = u_r - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        !     
        v_r   = Q(i,j+1,3)/Q(i,j+1,1)
        tmp1 = Q(i,j+2,3)/Q(i,j+2,1) - Q(i,j+1,3)/Q(i,j+1,1)
        tmp2 = Q(i,j+1,3)/Q(i,j+1,1) - Q(i,j,3)/Q(i,j,1)
        fluxLimPl  = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        v_r = v_r - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        !
        hoge   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
        pL     = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0   + Q(i,j,3)**2.0d0)/Q(i,j,1))
        fuga   = (GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1))
        tmp1 = hoge - pL
        tmp2 = pL - fuga
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)
        pL   = pL + 0.25d0*(1.0d0 - phi)*fluxLimMin + 0.25d0*(1.0d0 + phi)*fluxLimPl 
        !
        hoge   = (GAMMA-1.0d0)*(Q(i,j+2,4) - 0.5d0*(Q(i,j+2,2)**2.0d0 + Q(i,j+2,3)**2.0d0)/Q(i,j+2,1))
        pR     = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
        fuga   = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0   + Q(i,j,3)**2.0d0)/Q(i,j,1))
        tmp1 = hoge - pR
        tmp2 = pR - fuga
        fluxLimPl = minmod(tmp1,beta*tmp2)
        fluxLimMin = minmod(tmp2,beta*tmp1)          
        pR = pR - 0.25d0*(1.0d0 - phi)*fluxLimPl - 0.25d0*(1.0d0 + phi)*fluxLimMin
        
        UL = u_l*kx + v_l*ky
        UR = u_r*kx + v_r*ky
        VL = u_l*nx + v_l*ny
        VR = u_r*nx + v_r*ny
        
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
        !if(j==0 .or. j==1)then
        !    write(*,*)j,F(i,j,1),F(i,j,2),F(i,j,3),F(i,j,4)
        !endif

    end subroutine AusmDV_YFlux
end module AUSMDVFlux