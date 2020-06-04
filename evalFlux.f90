subroutine evalFlux(Q,E,F,m,n,Nx,Ny)
    implicit none
    integer,intent(in)  :: Nx,Ny
    real(8),intent(in)  :: m(-2:Nx+3,-2:Ny+3,2),n(-2:Nx+3,-2:Ny+3,2),is_SF_xi(-1:Nx+2,-1:Ny+2),is_SF_eta(-1:Nx+2,-1:Ny+2)
    real(8),intent(out) :: Q(-1:Nx+2,-1:Ny+2,4)
    real(8),intent(inout) :: E(-1:Nx+2,-1:Ny+2,4),F(-1:Nxnum+2,-1:Nynum+2,4)
    
    do j=0,Ny
        do i=0,Nx
            if(is_SF_xi(i,j) == 1)then
                call SD_SLAU_Xflux(Q,E,m,Nx,Ny,i,j)
            else
                call AusmDV_XFlux(Q,E,m,Nx,Ny,i,j)
                !call SD_SLAU_Xflux(Q,E,m,Nx,Ny,i,j)
            endif
            
            
            if(is_SF_eta(i,j) == 1)then
                call SD_SLAU_YFlux(Q,F,n,Nx,Ny,i,j)
            else
                call AusmDV_YFlux(Q,F,n,Nx,Ny,i,j)
                !call SD_SLAU_YFlux(Q,F,n,Nx,Ny,i,j)
            endif
    
            !call Xviscflux(Q,Ev,m,n,S,Nx,Ny)
            !call YviscFlux(Q,Fv,m,n,S,Nx,Ny)
        enddo
    enddo
end subroutine evalFlux