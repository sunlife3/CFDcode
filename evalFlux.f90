subroutine evalFlux(Q,E,F,E1,E2,F1,F2,m,n,Nx,Ny,is_SF_xi,is_SF_eta)
    use AUSMDVFlux
    use SD_SLAUFlux
    use HLLFlux
    implicit none
    integer,intent(in)  :: Nx,Ny
    integer i,j
    real(8),intent(in)  :: Q(-2:Nx+3,-2:Ny+3,4),m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2)!,&
                        !   is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2)
    integer,intent(in)  :: is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2)
    real(8),intent(inout) :: E(-1:Nx+2,-1:Ny+2,4),F(-1:Nx+2,-1:Ny+2,4),&
                             E1(-1:Nx+2,-1:Ny+2,4),F1(-1:Nx+2,-1:Ny+2,4),E2(-1:Nx+2,-1:Ny+2,4),F2(-1:Nx+2,-1:Ny+2,4)

    do j=0,Ny
        do i=0,Nx

            if(is_SF_xi(i,j) == 1)then
                !call HLL_XFlux(Q,E,m,Nx,Ny,i,j)
                !call AusmDV_XFlux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                !call AusmDV_XFlux(Q,E2,m,Nx,Ny,i,j,0.0d0)
                call HLL_XFlux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                call HLL_XFlux(Q,E2,m,Nx,Ny,i,j,2.5d0)
                !call SD_SLAU_Xflux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                !call SD_SLAU_Xflux(Q,E2,m,Nx,Ny,i,j,4.0d0)
                E(i,j,1) = 0.5d0*E1(i,j,1) + 0.5d0*E2(i,j,1)
                E(i,j,2) = 0.5d0*E1(i,j,2) + 0.5d0*E2(i,j,2)
                E(i,j,3) = 0.5d0*E1(i,j,3) + 0.5d0*E2(i,j,3)
                E(i,j,4) = 0.5d0*E1(i,j,4) + 0.5d0*E2(i,j,4)
            else if(is_SF_xi(i,j) == 2)then
                !call AusmDV_XFlux(Q,E,m,Nx,Ny,i,j,4.0d0)
                !call HLL_XFlux(Q,E,m,Nx,Ny,i,j)
                !call AusmDV_XFlux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                !call AusmDV_XFlux(Q,E2,m,Nx,Ny,i,j,2.5d0)
                call HLL_XFlux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                call HLL_XFlux(Q,E2,m,Nx,Ny,i,j,2.5d0)
                !call SD_SLAU_Xflux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                !call SD_SLAU_Xflux(Q,E2,m,Nx,Ny,i,j,4.0d0)
                E(i,j,1) = 0.25d0*E1(i,j,1) + 0.75d0*E2(i,j,1)
                E(i,j,2) = 0.25d0*E1(i,j,2) + 0.75d0*E2(i,j,2)
                E(i,j,3) = 0.25d0*E1(i,j,3) + 0.75d0*E2(i,j,3)
                E(i,j,4) = 0.25d0*E1(i,j,4) + 0.75d0*E2(i,j,4)
            else
                !call AusmDV_XFlux(Q,E,m,Nx,Ny,i,j,4.0d0)
                !call HLL_XFlux(Q,E,m,Nx,Ny,i,j,2.5d0)
                !call AusmDV_XFlux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                !call AusmDV_XFlux(Q,E,m,Nx,Ny,i,j,2.0d0)
                !call HLL_XFlux(Q,E1,m,Nx,Ny,i,j,0.0d0)
                call HLL_XFlux(Q,E,m,Nx,Ny,i,j,2.5d0)
                !call SD_SLAU_Xflux(Q,E1,m,Nx,Ny,i,j,2.5d0)
                !call SD_SLAU_Xflux(Q,E,m,Nx,Ny,i,j,2.5d0)
                !E(i,j,1) = 0.25d0*E1(i,j,1) + 0.75d0*E2(i,j,1)
                !E(i,j,2) = 0.25d0*E1(i,j,2) + 0.75d0*E2(i,j,2)
                !E(i,j,3) = 0.25d0*E1(i,j,3) + 0.75d0*E2(i,j,3)
                !E(i,j,4) = 0.25d0*E1(i,j,4) + 0.75d0*E2(i,j,4)
            endif

            if(is_SF_eta(i,j) == 1)then
               !call AusmDV_YFlux(Q,F1,n,Nx,Ny,i,j,4.0d0)
               !call AusmDV_YFlux(Q,F2,n,Nx,Ny,i,j,4.0d0)
               call HLL_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
               call HLL_YFlux(Q,F2,n,Nx,Ny,i,j,2.5d0)
               !call SD_SLAU_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
               !call SD_SLAU_YFlux(Q,F2,n,Nx,Ny,i,j,4.0d0)    
               F(i,j,1) = 0.5d0*F1(i,j,1) + 0.5d0*F2(i,j,1)
               F(i,j,2) = 0.5d0*F1(i,j,2) + 0.5d0*F2(i,j,2)
               F(i,j,3) = 0.5d0*F1(i,j,3) + 0.5d0*F2(i,j,3)
               F(i,j,4) = 0.5d0*F1(i,j,4) + 0.5d0*F2(i,j,4)
               !call HLL_YFlux(Q,F,n,Nx,Ny,i,j)
            
            else if(is_SF_eta(i,j) == 2)then
                !call AusmDV_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
                !call AusmDV_YFlux(Q,F2,n,Nx,Ny,i,j,2.5d0)
                call HLL_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
                call HLL_YFlux(Q,F2,n,Nx,Ny,i,j,2.5d0)
                !call SD_SLAU_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
                !call SD_SLAU_YFlux(Q,F2,n,Nx,Ny,i,j,4.0d0)    
                F(i,j,1) = 0.25d0*F1(i,j,1) + 0.75d0*F2(i,j,1)
                F(i,j,2) = 0.25d0*F1(i,j,2) + 0.75d0*F2(i,j,2)
                F(i,j,3) = 0.25d0*F1(i,j,3) + 0.75d0*F2(i,j,3)
                F(i,j,4) = 0.25d0*F1(i,j,4) + 0.75d0*F2(i,j,4)
            else
                !call AusmDV_YFlux(Q,F,n,Nx,Ny,i,j,4.0d0)
                !call HLL_YFlux(Q,F,n,Nx,Ny,i,j,2.5d0)
                !call AusmDV_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
                !call AusmDV_YFlux(Q,F,n,Nx,Ny,i,j,2.0d0)
                !call HLL_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
                call HLL_YFlux(Q,F,n,Nx,Ny,i,j,2.5d0)
                !call SD_SLAU_YFlux(Q,F1,n,Nx,Ny,i,j,0.0d0)
                !call SD_SLAU_YFlux(Q,F,n,Nx,Ny,i,j,2.5d0)
                !F(i,j,1) = 0.25d0*F1(i,j,1) + 0.75d0*F2(i,j,1)
                !F(i,j,2) = 0.25d0*F1(i,j,2) + 0.75d0*F2(i,j,2)
                !F(i,j,3) = 0.25d0*F1(i,j,3) + 0.75d0*F2(i,j,3)
                !F(i,j,4) = 0.25d0*F1(i,j,4) + 0.75d0*F2(i,j,4)
            endif

        enddo
    enddo

end subroutine evalFlux
