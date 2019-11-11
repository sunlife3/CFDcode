module SLAUFlux
    use param
    implicit none
contains
    subroutine SLAU2_XFlux(Q,E,m,Nx,Ny)
        implicit none
        integer,intent(in) :: Nx,Ny
        integer :: i=0,j=0,n=0
        real(8),intent(in) :: Q(-1:Nx+2,-1:Ny+2,4),m(-2:Nx+3,-2:Ny+3,2)
        real(8),intent(out) :: E(-1:Nx+2,-1:Ny+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0
        real(8) rhoL,rhoR,rhoBar,UL,UR,VL,VR,pL,pR,HL,HR,cL,cR,cBar,MP,MM,MHat,VnL,VnR,&
                absVnBar,chi,g,beta_P,beta_M,massFlow,Ptilde,&
                fluxLmtP,fluxlmtM,beta,deltaP,deltaM,hoge,fuga,&
                mx,my,kx,ky,mmag,kmag,E1,E2,E3,E4
        
        n=n+1
        beta = 0.5d0*(1.0d0 + (3.0d0-phi)/(1.0d0-phi))
        !do j=0,Ny
        !    do i=0,Nx
        !        write(*,*)i,j,Q(i,j,1),Q(i,j,2),Q(i,j,3),Q(i,j,4)
        !    enddo
        !enddo

        do j=0,Ny
            do i=0,Nx
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
                !write(*,*)i,j,mx,my,kx,ky

                deltaP = Q(i+1,j,1) - Q(i,j,1)
                deltaM = Q(i,j,1) - Q(i-1,j,1)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)
                rhoL = Q(i,j,1) + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP

                deltaP = Q(i+2,j,1) - Q(i+1,j,1)
                deltaM = Q(i+1,j,1) - Q(i,j,1)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)
                rhoR = Q(i+1,j,1) - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
                rhoBar = 0.5d0*(rhoL + rhoR)

                UL   = (Q(i,j,2)/Q(i,j,1))*mx + (Q(i,j,3)/Q(i,j,1))*my                
                deltaP = Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my&
                        - (Q(i,j,2)/Q(i,j,1)*mx + Q(i,j,3)/Q(i,j,1)*my) 
                deltaM = Q(i,j,2)/Q(i,j,1)*mx + Q(i,j,3)/Q(i,j,1)*my&
                        - (Q(i-1,j,2)/Q(i-1,j,1)*mx + Q(i-1,j,3)/Q(i-1,j,1)*my)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)            
                UL   = UL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP

                UR = Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my
                deltaP = Q(i+2,j,2)/Q(i+2,j,1)*mx + Q(i+2,j,3)/Q(i+2,j,1)*my& 
                        - (Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my)
                deltaM = Q(i+1,j,2)/Q(i+1,j,1)*mx + Q(i+1,j,3)/Q(i+1,j,1)*my&
                        - (Q(i,j,2)/Q(i,j,1)*mx + Q(i,j,3)/Q(i,j,1)*my)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)   
                UR = UR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM

                VL   = (Q(i,j,2)/Q(i,j,1))*kx + (Q(i,j,3)/Q(i,j,1))*ky
                deltaP = Q(i+1,j,2)/Q(i+1,j,1)*kx + Q(i+1,j,3)/Q(i+1,j,1)*ky&
                       - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky)
                deltaM = Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky&
                       - (Q(i-1,j,2)/Q(i-1,j,1)*kx + Q(i-1,j,3)/Q(i-1,j,1)*ky)
                fluxLmtM = minmod(deltaP,beta*deltaM)
                fluxLmtP  = minmod(deltaM,beta*deltaP)
                VL = VL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP

                VR = (Q(i+1,j,2)/Q(i+1,j,1))*kx + (Q(i+1,j,3)/Q(i+1,j,1))*ky
                deltaP = Q(i+2,j,2)/Q(i+2,j,1)*kx + Q(i+2,j,3)/Q(i+2,j,1)*ky&
                        - (Q(i+1,j,2)/Q(i+1,j,1)*kx + Q(i+1,j,3)/Q(i+1,j,1)*ky)
                deltaM = Q(i+1,j,2)/Q(i+1,j,1)*kx + Q(i+1,j,3)/Q(i+1,j,1)*ky&
                        - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)   
                VR = VR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
                
                VnL  = UL
                VnR  = UR

                hoge   = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                pL     = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
                fuga   = (GAMMA-1.0d0)*(Q(i-1,j,4) - 0.5d0*(Q(i-1,j,2)**2.0d0 + Q(i-1,j,3)**2.0d0)/Q(i-1,j,1))
                deltaP = hoge - pL
                deltaM = pL - fuga
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)            
                pL   = pL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
                
                hoge   = (GAMMA-1.0d0)*(Q(i+2,j,4) - 0.5d0*(Q(i+2,j,2)**2.0d0 + Q(i+2,j,3)**2.0d0)/Q(i+2,j,1))
                pR     = (GAMMA-1.0d0)*(Q(i+1,j,4) - 0.5d0*(Q(i+1,j,2)**2.0d0 + Q(i+1,j,3)**2.0d0)/Q(i+1,j,1))
                fuga   = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
                deltaP = hoge - pR
                deltaM = pR - fuga
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)                
                pR   = pR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
                
                HL   = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0))/rhoL
                HR   = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0))/rhoR

                cL   = sqrt(GAMMA*pL/rhoL)
                cR   = sqrt(GAMMA*pR/rhoR)
                cBar = 0.5d0*(cL + cR)
                MP   = VnL/cBar
                MM   = VnR/cBar

                g = - max(min(MP,0.0d0),-1.0d0) * min(max(MM,0.0d0),1.0d0)
                absVnBar  = (rhoL*abs(VnL) + rhoR*abs(VnR))/(rhoL+rhoR)
                !absVnBar_P = (1.0d0 - g)*absVnBar + g*abs(VnL)
                !absVnBar_M = (1.0d0 - g)*absVnBar + g*abs(VnR)
                MHat  = min(1.0d0,sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))/cBar)
                chi = (1.0d0 - MHat)**2.0d0
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

                Ptilde = 0.5d0*(pL + pR) + 0.5d0*(beta_P - beta_M)*(pL - pR)&
                            + (beta_P + beta_M - 1.0d0)*sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))*rhoBar*cBar!SLAU2
                            !+ 0.5d0*(1.0d0 - chi)*(beta_P + beta_M - 1.0d0)*(pL+pR)                   !SLAU


                E1 = massFlow
                E2 = 0.5d0*(massFlow*(UR + UL) - abs(massFlow)*(UR - UL)) + Ptilde
                E3 = 0.5d0*(massFlow*(VR + VL) - abs(massFlow)*(VR - VL))
                E4 = 0.5d0*(massFlow*(HR + HL) - abs(massFlow)*(HR - HL))
                !write(*,*)i,j,'E',E1,E2,E3,E4
                
                E(i,j,1) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*E1
                E(i,j,2) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*(mx*E2 + kx*E3)
                E(i,j,3) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*(my*E2 + ky*E3)
                E(i,j,4) = sqrt(m(i,j,1)**2.0d0 + m(i,j,2)**2.0d0)*E4
                !write(*,*)i,E(i,j,1),E(i,j,2),E(i,j,3),E(i,j,4)
            enddo
        enddo
    end subroutine SLAU2_XFlux

    subroutine SLAU2_YFlux(Q,F,n,Nxnum,Nynum)
        implicit none
        integer,intent(in) :: Nxnum,Nynum
        integer :: i=0,j=0
        real(8),intent(in) :: Q(-1:Nxnum+2,-1:Nynum+2,4),n(-2:Nxnum+3,-2:Nynum+3,2)
        real(8),intent(out) :: F(-1:Nxnum+2,-1:Nynum+2,4)
        real(8),parameter :: phi = 1.0d0/3.0d0
        real(8) rhoL,rhoR,rhoBar,UL,UR,VL,VR,pL,pR,HL,HR,cL,cR,cBar,MP,MM,MHat,VnL,VnR,&
                absVnBar,chi,g,beta_P,beta_M,massFlow,Ptilde,&
                fluxLmtP,fluxlmtM,beta,deltaP,deltaM,hoge,fuga,&
                nx,ny,kx,ky,nmag,kmag,F1,F2,F3,F4
        
        beta = 0.5d0*(1.0d0 + (3.0d0-phi)/(1.0d0-phi))

        do j=0,Nynum
            do i=0,Nxnum
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
                !write(*,*)i,j,mx,my,kx,ky

                deltaP = Q(i,j+1,1) - Q(i,j,1)
                deltaM = Q(i,j,1) - Q(i,j-1,1)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)
                rhoL = Q(i,j,1) + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP

                deltaP = Q(i,j+2,1) - Q(i,j+1,1)
                deltaM = Q(i,j+1,1) - Q(i,j,1)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)
                rhoR = Q(i,j+1,1) - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
                rhoBar = 0.5d0*(rhoL + rhoR)

                UL   = (Q(i,j,2)/Q(i,j,1))*kx + (Q(i,j,3)/Q(i,j,1))*ky                
                deltaP = Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j+1,1)*ky&
                        - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky) 
                deltaM = Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky&
                        - (Q(i,j-1,2)/Q(i,j-1,1)*kx + Q(i,j-1,3)/Q(i,j-1,1)*ky)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)            
                UL   = UL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP

                UR = Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j+1,1)*ky
                deltaP = Q(i,j+2,2)/Q(i,j+2,1)*kx + Q(i,j+2,3)/Q(i,j+2,1)*ky& 
                        - (Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j+1,1)*ky)
                deltaM = Q(i,j+1,2)/Q(i,j+1,1)*kx + Q(i,j+1,3)/Q(i,j+1,1)*ky&
                        - (Q(i,j,2)/Q(i,j,1)*kx + Q(i,j,3)/Q(i,j,1)*ky)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)   
                UR = UR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM

                VL   = (Q(i,j,2)/Q(i,j,1))*nx + (Q(i,j,3)/Q(i,j,1))*ny
                deltaP = Q(i,j+1,2)/Q(i,j+1,1)*nx + Q(i,j+1,3)/Q(i,j+1,1)*ny&
                       - (Q(i,j,2)/Q(i,j,1)*nx + Q(i,j,3)/Q(i,j,1)*ny)
                deltaM = Q(i,j,2)/Q(i,j,1)*nx + Q(i,j,3)/Q(i,j,1)*ny&
                       - (Q(i,j-1,2)/Q(i,j-1,1)*nx + Q(i,j-1,3)/Q(i,j-1,1)*ny)
                fluxLmtM = minmod(deltaP,beta*deltaM)
                fluxLmtP  = minmod(deltaM,beta*deltaP)
                VL = VL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP

                VR = (Q(i,j+1,2)/Q(i,j+1,1))*nx + (Q(i,j+1,3)/Q(i,j+1,1))*ny
                deltaP = Q(i,j+2,2)/Q(i,j+2,1)*nx + Q(i,j+2,3)/Q(i,j+2,1)*ny&
                        - (Q(i,j+1,2)/Q(i,j+1,1)*nx + Q(i,j+1,3)/Q(i,j+1,1)*ny)
                deltaM = Q(i,j+1,2)/Q(i,j+1,1)*nx + Q(i,j+1,3)/Q(i,j+1,1)*ny&
                        - (Q(i,j,2)/Q(i,j,1)*nx + Q(i,j,3)/Q(i,j,1)*ny)
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)
                VR = VR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
                
                VnL  = VL
                VnR  = VR

                hoge   = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
                pL     = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
                fuga   = (GAMMA-1.0d0)*(Q(i,j-1,4) - 0.5d0*(Q(i,j-1,2)**2.0d0 + Q(i,j-1,3)**2.0d0)/Q(i,j-1,1))
                deltaP = hoge - pL
                deltaM = pL - fuga
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)            
                pL   = pL + 0.25d0*(1.0d0 - phi)*fluxLmtM + 0.25d0*(1.0d0 + phi)*fluxLmtP
                
                hoge   = (GAMMA-1.0d0)*(Q(i,j+2,4) - 0.5d0*(Q(i,j+2,2)**2.0d0 + Q(i,j+2,3)**2.0d0)/Q(i,j+2,1))
                pR     = (GAMMA-1.0d0)*(Q(i,j+1,4) - 0.5d0*(Q(i,j+1,2)**2.0d0 + Q(i,j+1,3)**2.0d0)/Q(i,j+1,1))
                fuga   = (GAMMA-1.0d0)*(Q(i,j,4)   - 0.5d0*(Q(i,j,2)**2.0d0 + Q(i,j,3)**2.0d0)/Q(i,j,1))
                deltaP = hoge - pR
                deltaM = pR - fuga
                fluxLmtP = minmod(deltaP,beta*deltaM)
                fluxLmtM = minmod(deltaM,beta*deltaP)                
                pR   = pR - 0.25d0*(1.0d0 - phi)*fluxLmtP - 0.25d0*(1.0d0 + phi)*fluxLmtM
                
                !write(*,*)i,j,'L',rhoL,UL,VL,pL
                !write(*,*)i,j,'R',rhoR,UR,VR,pR

                HL   = (GAMMA*pL/(GAMMA-1.0d0) + 0.5d0*rhoL*(UL**2.0d0 + VL**2.0d0))/rhoL
                HR   = (GAMMA*pR/(GAMMA-1.0d0) + 0.5d0*rhoR*(UR**2.0d0 + VR**2.0d0))/rhoR

                cL   = sqrt(GAMMA*pL/rhoL)
                cR   = sqrt(GAMMA*pR/rhoR)
                cBar = 0.5d0*(cL + cR)
                MP   = VnL/cBar
                MM   = VnR/cBar

                g = - max(min(MP,0.0d0),-1.0d0) * min(max(MM,0.0d0),1.0d0)
                absVnBar  = (rhoL*abs(VnL) + rhoR*abs(VnR))/(rhoL+rhoR)
                !absVnBar_P = (1.0d0 - g)*absVnBar + g*abs(VnL)
                !absVnBar_M = (1.0d0 - g)*absVnBar + g*abs(VnR)
                MHat  = min(1.0d0,sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))/cBar)
                chi = (1.0d0 - MHat)**2.0d0
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

                Ptilde = 0.5d0*(pL + pR) + 0.5d0*(beta_P - beta_M)*(pL - pR)&
                            + (beta_P + beta_M - 1.0d0)*sqrt(0.5d0*(UL**2.0d0 + VL**2.0d0 + UR**2.0d0 + VR**2.0d0))*rhoBar*cBar!SLAU2
                            !+ 0.5d0*(1.0d0 - chi)*(beta_P + beta_M - 1.0d0)*(pL+pR)                   !SLAU


                F1 = massFlow
                F2 = 0.5d0*(massFlow*(UR + UL) - abs(massFlow)*(UR - UL))
                F3 = 0.5d0*(massFlow*(VR + VL) - abs(massFlow)*(VR - VL)) + Ptilde
                F4 = 0.5d0*(massFlow*(HR + HL) - abs(massFlow)*(HR - HL))
                !write(*,*)i,j,'F',F1,F2,F3,F4
                
                F(i,j,1) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*F1
                F(i,j,2) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*(kx*F2 + nx*F3)
                F(i,j,3) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*(ky*F2 + ny*F3)
                F(i,j,4) = sqrt(n(i,j,1)**2.0d0 + n(i,j,2)**2.0d0)*F4
                !write(*,*)i,F(i,j,1)-F(i,j-1,1),F(i,j,2)-F(i,j-1,2),F(i,j,3)-F(i,j-1,3),F(i,j,4)-F(i,j-1,4)
            enddo
        enddo
    end subroutine SLAU2_YFlux

end module SLAUFlux