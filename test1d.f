!>
!>    @brief Determines if the model is 1D
!>
!>    @param[out] is3d    True if model is 3D.
!>                        False is model is 1D.
!>
      SUBROUTINE TEST1D(is3d) !tests if a model is 1d
      USE MODEL_MODULE, ONLY : da, mu, rho, nx, ny, nz
      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: is3d
      COMPLEX da1, mu1, rh1 
      INTEGER ix, iz, iy

      is3d = .TRUE.
      DO iz = 1, nz !loop on depth
         da1 = da(iz,1,1)
         mu1 = mu(iz,1,1)
         rh1 = rho(iz,1,1)
         DO iy = 2, ny !scan horizontally  
            DO ix = 2, nx
               IF (da(iz,ix,iy).ne.da1)  RETURN
               IF (mu(iz,ix,iy).ne.mu1)  RETURN
               IF (rho(iz,ix,iy).ne.rh1) RETURN
            ENDDO
         ENDDO
      ENDDO
      is3d = .FALSE. !model is 1d
      return
      end

