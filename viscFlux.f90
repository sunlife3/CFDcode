subroutine XviscFlux(Q,Ev,m,n,S,Nxnum,Nynum)
    implicit none
    integer,intent(in) :: Nxnum,Nynum
    integer,parameter :: idxLeadEdge = 10,idxTrailEdge=500
    integer :: i=0,j=0
    real(8),intent(in) :: Q(-2:Nxnum+3,-2:Nynum+3,4),S(-2:Nxnum+3,-2:Nynum+3),&
                            m(-2:Nxnum+3,-2:Nynum+3,2),n(-2:Nxnum+3,-2:Nynum+3,2)
    real(8),intent(out) :: Ev(-1:Nxnum+2,-1:Nynum+2,4)
    real(8),parameter :: GAMMA = 1.4d0,Pr = 1.0d0
    real(8) :: mu,uxi,ueta,vxi,veta,Txi,Teta,ux,uy,vx,vy,Tx,Ty,&
                uifp,uifm,vifp,vifm,Tifp,Tifm,nx,ny,mxHat,myHat,kxHat,kyHat,nxHat,nyHat,Jac,&
                tauxx,tauxy,tauyy,uif,vif,Tif
    real(8) rho(-1:Nxnum+2,-1:Nynum+2),u(-1:Nxnum+2,-1:Nynum+2),v(-1:Nxnum+2,-1:Nynum+2),T(-1:Nxnum+2,-1:Nynum+2)

    
    do j=-1,Nynum+2
        do i=-1,Nxnum+2
            rho(i,j) = Q(i,j,1)
            u(i,j)   = Q(i,j,2)/Q(i,j,1)
            v(i,j)   = Q(i,j,3)/Q(i,j,1)
            T(i,j)   = (GAMMA-1.0d0)*(Q(i,j,4)/Q(i,j,1) - 0.5d0*(u(i,j)**2.0d0 + v(i,j)**2.0d0))
        enddo
    enddo


    do j=0,Nynum
        do i=0,Nxnum
            if(j == 1 .and. idxLeadEdge-1 <= i)then !if(j == 1)then!
            !viscous Flux at wall (j = 1)
                !calculation ...xi,...eta
                uxi = u(i+1,j) - u(i,j)
                vxi = v(i+1,j) - v(i,j)
                Txi = T(i+1,j) - T(i,j)

                !ここのuをUにする    
                uifp = (u(i+1,j+1)*sqrt(rho(i+1,j+1)) + u(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i+1,j+1)) + sqrt(rho(i,j+1))) !uC  |    |C __j=2
                uifm = (u(i+1,j)*sqrt(rho(i+1,j)) + u(i,j)*sqrt(rho(i,j)))/(sqrt(rho(i+1,j)) + sqrt(rho(i,j)))             !uB  |____|
                ueta = (uifp + 3.0d0*uifm)/3.0d0                                                                           !    |    |B __j=1
                vifp = (v(i+1,j+1)*sqrt(rho(i+1,j+1)) + v(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i+1,j+1)) + sqrt(rho(i,j+1))) !vC  |____!A
                vifm = (v(i+1,j)*sqrt(rho(i+1,j)) + v(i,j)*sqrt(rho(i,j)))/(sqrt(rho(i+1,j)) + sqrt(rho(i,j)))             !vB  |////|
                veta = (vifp + 3.0d0*vifm)/3.0d0                                                                           !    |wall|
                Tifp = (T(i+1,j+1)*sqrt(rho(i+1,j+1)) + T(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i+1,j+1)) + sqrt(rho(i,j+1))) !TC 
                Tifm = (T(i+1,j)*sqrt(rho(i+1,j)) + T(i,j)*sqrt(rho(i,j)))/(sqrt(rho(i+1,j)) + sqrt(rho(i,j)))             !TB
                Teta = 0.5d0*(Tifp - Tifm)

            else
                !calculation ...xi,...eta
                !ここのuをUにする
                uxi = u(i+1,j) - u(i,j)
                vxi = v(i+1,j) - v(i,j)
                Txi = T(i+1,j) - T(i,j)

                !ここのuをUにする                
                uifp = (u(i+1,j+1)*sqrt(rho(i+1,j+1)) + u(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i+1,j+1)) + sqrt(rho(i,j+1)))
                uifm = (u(i+1,j-1)*sqrt(rho(i+1,j-1)) + u(i,j-1)*sqrt(rho(i,j-1)))/(sqrt(rho(i+1,j-1)) + sqrt(rho(i,j-1)))
                ueta = 0.5d0*(uifp - uifm)
                vifp = (v(i+1,j+1)*sqrt(rho(i+1,j+1)) + v(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i+1,j+1)) + sqrt(rho(i,j+1)))
                vifm = (v(i+1,j-1)*sqrt(rho(i+1,j-1)) + v(i,j-1)*sqrt(rho(i,j-1)))/(sqrt(rho(i+1,j-1)) + sqrt(rho(i,j-1)))
                veta = 0.5d0*(vifp - vifm)
                Tifp = (T(i+1,j+1)*sqrt(rho(i+1,j+1)) + T(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i+1,j+1)) + sqrt(rho(i,j+1)))
                Tifm = (T(i+1,j-1)*sqrt(rho(i+1,j-1)) + T(i,j-1)*sqrt(rho(i,j-1)))/(sqrt(rho(i+1,j-1)) + sqrt(rho(i,j-1)))
                Teta = 0.5d0*(Tifp - Tifm)                
            endif

            !calculation ...x,...y
            Jac = 0.5d0*(1.0d0/S(i+1,j) + 1.0d0/S(i,j))
            nx  = 0.25d0*(n(i,j,1) + n(i,j-1,1) + n(i+1,j,1) + n(i+1,j-1,1))
            ny  = 0.25d0*(n(i,j,2) + n(i,j-1,2) + n(i+1,j,2) + n(i+1,j-1,2)) 
            ux  = (m(i,j,1)*uxi + nx*ueta)*Jac
            vx  = (m(i,j,1)*vxi + nx*veta)*Jac
            Tx  = (m(i,j,1)*Txi + nx*Teta)*Jac
            uy  = (m(i,j,2)*uxi + ny*ueta)*Jac
            vy  = (m(i,j,2)*vxi + ny*veta)*Jac
            Ty  = (m(i,j,2)*Txi + ny*Teta)*Jac
        
            !calculation shear stress
            uif  = (u(i+1,j)*sqrt(rho(i+1,j)) + u(i,j)*sqrt(rho(i,j)))/(sqrt(rho(i+1,j)) + sqrt(rho(i,j)))
            vif  = (v(i+1,j)*sqrt(rho(i+1,j)) + v(i,j)*sqrt(rho(i,j)))/(sqrt(rho(i+1,j)) + sqrt(rho(i,j)))
            Tif  = (T(i+1,j)*sqrt(rho(i+1,j)) + T(i,j)*sqrt(rho(i,j)))/(sqrt(rho(i+1,j)) + sqrt(rho(i,j)))
            mu   = GAMMA*Tif
            !mu = (GAMMA*Tif)**(3.0d0/2.0d0)*((1/GAMMA + S)/(T + S))
            tauxx = mu*(2.0d0/3.0d0)*(2.0d0*ux-vy)
            tauxy = mu*(uy+vx)
            tauyy = mu*(2.0d0/3.0d0)*(2.0d0*vy-ux)
            !write(*,*)tauxx,tauxy,tauyy
        
            !calculation viscous flux
            Ev(i,j,1) = 0.0d0
            Ev(i,j,2) = m(i,j,1)*tauxx + m(i,j,2)*tauxy
            Ev(i,j,3) = m(i,j,1)*tauxy + m(i,j,2)*tauyy
            Ev(i,j,4) = m(i,j,1)*(uif*tauxx + vif*tauxy + (GAMMA/(GAMMA-1.0d0)*mu/Pr*Tx))+&
                        m(i,j,2)*(uif*tauxy + vif*tauyy + (GAMMA/(GAMMA-1.0d0)*mu/Pr*Ty))
            !if(0<=j .and. j<=1)then
            !    write(*,*)i,j,Ev(i,j,1)/5000.0d0,Ev(i,j,2)/5000.0d0,Ev(i,j,3)/5000.0d0,Ev(i,j,4)/5000.0d0
            !endif
        enddo
    enddo

end subroutine XviscFlux

subroutine YviscFlux(Q,Fv,m,n,S,Nxnum,Nynum)
    implicit none
    integer,intent(in) :: Nxnum,Nynum
    integer :: i=0,j=0,HOGE=0
    integer,parameter :: idxLeadEdge = 10,idxTrailEdge=500
    real(8),intent(in) :: Q(-2:Nxnum+3,-2:Nynum+3,4),S(-2:Nxnum+3,-2:Nynum+3),&
                            m(-2:Nxnum+3,-2:Nynum+3,2),n(-2:Nxnum+3,-2:Nynum+3,2)
    real(8),intent(out) :: Fv(-1:Nxnum+2,-1:Nynum+2,4)
    real(8),parameter :: GAMMA = 1.4d0,Pr = 1.0d0
    real(8) :: mu,uxi,ueta,vxi,veta,Txi,Teta,ux,uy,vx,vy,Tx,Ty,&
                uifp,uifm,vifp,vifm,Tifp,Tifm,Ta,mx,my,Jac,&
                tauxx,tauxy,tauyy,Tif,uif,vif
    real(8) rho(-1:Nxnum+2,-1:Nynum+2),u(-1:Nxnum+2,-1:Nynum+2),v(-1:Nxnum+2,-1:Nynum+2),T(-1:Nxnum+2,-1:Nynum+2)
    
    do j=-1,Nynum+2
        do i=-1,Nxnum+2
            rho(i,j) = Q(i,j,1)
            u(i,j)   = Q(i,j,2)/Q(i,j,1)
            v(i,j)   = Q(i,j,3)/Q(i,j,1)
            T(i,j)   = (GAMMA-1.0d0)*(Q(i,j,4)/Q(i,j,1) - 0.5d0*(u(i,j)**2.0d0 + v(i,j)**2.0d0))
            !if(j <= 1)then
            !    write(*,*)i,j,rho(i,j),u(i,j),v(i,j),T(i,j)
            !endif
        enddo
    enddo

    do j=0,Nynum
        do i=0,Nxnum
            if(j == 0)then
            !viscous Flux at wall(j = 0)
                if(1<=i .and. i<=idxLeadEdge-1)then
                    !calculation ...xi,...eta
                    ueta = u(i,j+1) - u(i,j)
                    veta = v(i,j+1) - v(i,j)
                    Teta = T(i,j+1) - T(i,j)

                    uifp = (u(i+1,j)*sqrt(rho(i+1,j)) + u(i+1,j+1)*sqrt(rho(i+1,j+1)))/(sqrt(rho(i+1,j)) + sqrt(rho(i+1,j+1)))
                    uifm = (u(i-1,j)*sqrt(rho(i-1,j)) + u(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                    uxi  = 0.5d0*(uifp - uifm)
                    vifp = (v(i+1,j)*sqrt(rho(i+1,j)) + v(i+1,j+1)*sqrt(rho(i+1,j+1)))/(sqrt(rho(i+1,j)) + sqrt(rho(i+1,j+1)))
                    vifm = (v(i-1,j)*sqrt(rho(i-1,j)) + v(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                    vxi  = 0.5d0*(vifp - vifm)
                    Tifp = (T(i+1,j)*sqrt(rho(i+1,j)) + T(i+1,j+1)*sqrt(rho(i+1,j+1)))/(sqrt(rho(i+1,j)) + sqrt(rho(i+1,j+1)))
                    Tifm = (T(i-1,j)*sqrt(rho(i-1,j)) + T(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                    Txi  = 0.5d0*(Tifp - Tifm)

                    Tif  = (T(i,j)*sqrt(rho(i,j)) + T(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))
                    uif  = (u(i,j)*sqrt(rho(i,j)) + u(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))
                    vif  = (v(i,j)*sqrt(rho(i,j)) + v(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))

                else
                    !calculation ...xi
                    if(i == idxLeadEdge-1)then ! at the leading edge
                        uifp = 0.0d0
                        uifm = (u(i-1,j)*sqrt(rho(i-1,j)) + u(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                        uxi  = 0.5d0*(uifp - uifm)
                        vifp = 0.0d0
                        vifm = (v(i-1,j)*sqrt(rho(i-1,j)) + v(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                        vxi  = 0.5d0*(vifp - vifm)
                        Tifp = 0.125d0*(9.0d0*T(i+1,j+1) - T(i+1,j+2))
                        Tifm = (T(i-1,j)*sqrt(rho(i-1,j)) + T(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                        Txi  = 0.5d0*(Tifp - Tifm)
                    else                                            ! on the wall
                        uxi = 0.0d0
                        vxi = 0.0d0
                        Tifp = 0.125d0*(9.0d0*T(i+1,j+1) - T(i+1,j+2))
                        Tifm = 0.125d0*(9.0d0*T(i-1,j+1) - T(i-1,j+2))
                        Txi  = 0.5d0*(Tifp - Tifm)
                    endif
                    
                    !calculation ...eta
                    if(i == idxLeadEdge-1)then
                        ueta = u(i,j+1) - u(i,j)
                        veta = v(i,j+1) - v(i,j)
                        Teta = T(i,j+1) - T(i,j)

                        Tif  = (T(i,j)*sqrt(rho(i,j)) + T(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))
                        uif  = (u(i,j)*sqrt(rho(i,j)) + u(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))
                        vif  = (v(i,j)*sqrt(rho(i,j)) + v(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))

                    else ! on the wall
                        ueta = (9.0d0*u(i,j+1) - u(i,j+2))/3.0d0
                        veta = (9.0d0*v(i,j+1) - v(i,j+2))/3.0d0
                        Teta = 0.0d0
                        
                        Tif  = 0.125d0*(9.0d0*T(i,j+1) - T(i,j+2))
                        uif  = 0.0d0
                        vif  = 0.0d0
                    endif   
                endif

            else
                !calculation ...xi,...eta
                ueta = u(i,j+1) - u(i,j)
                veta = v(i,j+1) - v(i,j)
                Teta = T(i,j+1) - T(i,j)

                uifp = (u(i+1,j)*sqrt(rho(i+1,j)) + u(i+1,j+1)*sqrt(rho(i+1,j+1)))/(sqrt(rho(i+1,j)) + sqrt(rho(i+1,j+1)))
                uifm = (u(i-1,j)*sqrt(rho(i-1,j)) + u(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                uxi  = 0.5d0*(uifp - uifm)
                vifp = (v(i+1,j)*sqrt(rho(i+1,j)) + v(i+1,j+1)*sqrt(rho(i+1,j+1)))/(sqrt(rho(i+1,j)) + sqrt(rho(i+1,j+1)))
                vifm = (v(i-1,j)*sqrt(rho(i-1,j)) + v(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                vxi  = 0.5d0*(vifp - vifm)
                Tifp = (T(i+1,j)*sqrt(rho(i+1,j)) + T(i+1,j+1)*sqrt(rho(i+1,j+1)))/(sqrt(rho(i+1,j)) + sqrt(rho(i+1,j+1)))
                Tifm = (T(i-1,j)*sqrt(rho(i-1,j)) + T(i-1,j+1)*sqrt(rho(i-1,j+1)))/(sqrt(rho(i-1,j)) + sqrt(rho(i-1,j+1)))
                Txi  = 0.5d0*(Tifp - Tifm)

                Tif  = (T(i,j)*sqrt(rho(i,j)) + T(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))
                uif  = (u(i,j)*sqrt(rho(i,j)) + u(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))
                vif  = (v(i,j)*sqrt(rho(i,j)) + v(i,j+1)*sqrt(rho(i,j+1)))/(sqrt(rho(i,j)) + sqrt(rho(i,j+1)))
            endif

            !calculation ...x,...y
            Jac = 0.5d0*(1.0d0/S(i,j) + 1.0d0/S(i,j+1))
            mx  = 0.25d0*(m(i-1,j+1,1) + m(i-1,j,1) + m(i,j+1,1) + m(i,j,1))
            my  = 0.25d0*(m(i-1,j+1,2) + m(i-1,j,2) + m(i,j+1,2) + m(i,j,2))
            ux  = Jac*(mx*uxi + n(i,j,1)*ueta)
            vx  = Jac*(mx*vxi + n(i,j,1)*veta)
            Tx  = Jac*(mx*Txi + n(i,j,1)*Teta)
            uy  = Jac*(my*uxi + n(i,j,2)*ueta)
            vy  = Jac*(my*vxi + n(i,j,2)*veta)
            Ty  = Jac*(my*Txi + n(i,j,2)*Teta)

            !calculation shear stress
            mu = GAMMA*Tif
            !mu = (GAMMA*Tif)**(3.0d0/2.0d0)*((1/GAMMA + S)/(T + S))
            tauxx = mu*(2.0d0/3.0d0)*(2.0d0*ux-vy)
            tauxy = mu*(uy+vx)
            tauyy = mu*(2.0d0/3.0d0)*(2.0d0*vy-ux)
            
            !if(j <= 3)then
            !    write(*,*)i,j,uif,vif,Tif
            !endif
            !calculation viscous flux
            Fv(i,j,1) = 0.0d0 
            Fv(i,j,2) = n(i,j,1)*tauxx + n(i,j,2)*tauxy
            Fv(i,j,3) = n(i,j,1)*tauxy + n(i,j,2)*tauyy
            Fv(i,j,4) = n(i,j,1)*(uif*tauxx + vif*tauxy + (GAMMA/(GAMMA-1.0d0))*mu*Tx/Pr)+&
                        n(i,j,2)*(uif*tauxy + vif*tauyy + (GAMMA/(GAMMA-1.0d0))*mu*Ty/Pr)
            
            !if(0<=j .and. j<=1)then
            !    write(*,*)i,j,Fv(i,j,1)/5000.0d0,Fv(i,j,2)/5000.0d0,Fv(i,j,3)/5000.0d0,Fv(i,j,4)/5000.0d0
            !endif
        enddo
    enddo
end subroutine YviscFlux