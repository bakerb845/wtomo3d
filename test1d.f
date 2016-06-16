      SUBROUTINE test1d(is3d) !tests if a model is 1d
c     include 'dimension.inc'
c     include 'init.inc'
c     include 'common.inc'
      USE MODEL_MODULE, ONLY : da, mu, rho, nx, ny, nz
      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: is3d
      COMPLEX da1, mu1, rh1 
      INTEGER ix, iz, iy

      do iz = 1, nz !loop on depth
         da1 = da(iz,1,1)
         mu1 = mu(iz,1,1)
         rh1 = rho(iz,1,1)
         do iy = 2, ny !scan horizontally  
            do ix = 2, nx
               if (da(iz,ix,iy).ne.da1) go to 5
               if (mu(iz,ix,iy).ne.mu1) go to 5
               if (rho(iz,ix,iy).ne.rh1) go to 5
            enddo
         enddo
      enddo

      is3d = .false. !model is 1d
      return

    5 is3d = .true. !model is 3d

      return
      end

