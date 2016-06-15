      SUBROUTINE makeuh(isrc, sfld, ierr)
c
c     the driver program has called mumps and solved su_h=-(s-s_0)u_0 
c     for u_h.  we finish by adding the background field to u_h 
c     that is, u_t = u_h + u_0.  here sfld holds the heterogeneous 
c     displacement, u_h.   
c 
c 
      USE SOURCE_MODULE, ONLY : srctyp
      USE WAVEFIELD_MODULE, ONLY : utbkr, vtbkr, wtbkr, ufac, vfac, wfac
      USE RECPRM_MODULE, ONLY : xg, yg, zg, xr, yr, zr,
     ;                          sg, sr, receiver,
     ;                          uest, vest, west,
     ;                          utest,
     ;                          igreg, irreg, ng, nr, recin,
     ;                          gspread, rspread 
      USE MODEL_MODULE, ONLY : da, mu, iom, omega
      USE INIT_MODULE, ONLY : dx, dy, dz, nx, ny, nz, totfld
      USE INTERFACE_MODULE, ONLY : srcsetup 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: isrc
      COMPLEX, DIMENSION(:), INTENT(IN) :: sfld
      INTEGER, INTENT(OUT) :: ierr
      COMPLEX utmp, vtmp, wtmp, cfact, um1, up1,
     ;        dudx, dvdy, dwdz
      INTEGER ix, iy, iz, indx, igrdw, igrdv, igrdu,
     ;        ir, ig
c 
c----------------------------------------------------------------------c
c 
c.... copy solution to ?fac 
      ierr = 0
      !sfac = srcscl(isrc)*source(isrc,iom)
      !open(unit=31,file='mysln.txt')
      do iy = 1,ny
         do ix = 1,nx
            do iz = 1,nz
               indx = (iy-1)*nx*nz + (ix-1)*nz + iz
               igrdu = 3*indx - 2
               igrdv = igrdu + 1
               igrdw = igrdu + 2
               ufac(iz,ix,iy) = sfld(igrdu)
               vfac(iz,ix,iy) = sfld(igrdv)
               wfac(iz,ix,iy) = sfld(igrdw)
c
c.......... add on the background field if needed
               if (srctyp(isrc)(1:1).eq.'p'. or.
     ;             srctyp(isrc)(1:1).eq.'s') then
c
c............. note that we convolve with background wavefield with the stf in hetfldsrc so this 
c............. multiplication by sfac is no longer needed
c               ufac(iz,ix,iy) = ufac(iz,ix,iy) + sfac*utbkr(iz,ix,iy)
c               vfac(iz,ix,iy) = vfac(iz,ix,iy) + sfac*vtbkr(iz,ix,iy)
c               wfac(iz,ix,iy) = wfac(iz,ix,iy) + sfac*wtbkr(iz,ix,iy)
c............. set totfld to false if you want to see only
c............. the waves from the heterogeneous sources
                   if (totfld) then
                      ufac(iz,ix,iy) = ufac(iz,ix,iy) + utbkr(iz,ix,iy)
                      vfac(iz,ix,iy) = vfac(iz,ix,iy) + vtbkr(iz,ix,iy)
                      wfac(iz,ix,iy) = wfac(iz,ix,iy) + wtbkr(iz,ix,iy)
c--for a test, have look at the background field only
                   else
                      ufac(iz,ix,iy) = utbkr(iz,ix,iy)
                      vfac(iz,ix,iy) = vtbkr(iz,ix,iy)
                      wfac(iz,ix,iy) = wtbkr(iz,ix,iy)
                   endif
               endif
            enddo
         enddo
      enddo
c 
c.... turn off displacement output 
      do ig = 1,ng
c
c....... set up the receiver region 
         if (sg(isrc,ig).ne.0) then
            CALL srcsetup(.false.,
     ;                    xg(ig), yg(ig), zg(ig), gspread,
     ;                    igreg, utmp, ufac, ierr)
            CALL srcsetup(.false.,
     ;                    xg(ig), yg(ig), zg(ig), gspread,
     ;                    igreg, vtmp, vfac, ierr)
            CALL srcsetup(.false.,
     ;                    xg(ig), yg(ig), zg(ig), gspread,
     ;                    igreg, wtmp, wfac, ierr)
            uest(iom,ig,isrc) = utmp
            vest(iom,ig,isrc) = vtmp
            west(iom,ig,isrc) = wtmp

c
c.......... modify the receiver response if available
            if (recin.ne.0) then
               cfact = receiver(ig, iom)
               uest(iom,ig,isrc) = uest(iom,ig,isrc)*cfact
               vest(iom,ig,isrc) = vest(iom,ig,isrc)*cfact
               west(iom,ig,isrc) = west(iom,ig,isrc)*cfact
            endif
         else !set to zero if non-active source-receiver pair
            uest(iom,ig,isrc) = CMPLX(0.,0.)
            vest(iom,ig,isrc) = CMPLX(0.,0.)
            west(iom,ig,isrc) = CMPLX(0.,0.)
         endif
      enddo

c***********************************************************************
c
c pressure output -> projnm.utest
c
c       pressure is the mean normal stress = (s11 + s22 + s33)/3
c
c       s11 = lambda*(du/dx + dv/dy + dw/dz) + 2mu*du/dx
c       s22 = lambda*(du/dx + dv/dy + dw/dz) + 2mu*dv/dy
c       s33 = lambda*(du/dx + dv/dy + dw/dz) + 2mu*dw/dz
c
c       (s11 + s22 + s33)/3 = (lambda+2mu/3)*(du/dx + dv/dy + dw/dz)
c
c***********************************************************************

c.... use nr here b/c it applies to the number of hydrophones, but 
c.... fdatatio may not be sensitive to differences between nr and ng
      cfact = CMPLX(0.0,-1.0)*omega
      do ir = 1, nr
c 
c....... compute pressure at receivers from mean compressive stress
         if (sr(isrc,ir).ne.0.) then
            ix = INT(1.5 + xr(ir)/dx)
            iy = INT(1.5 + yr(ir)/dy)
            iz = INT(1.5 + zr(ir)/dz)

c.......... make sure ix and iz are in bounds
            if (ix.lt.1) ix = 1
            if (iy.lt.1) iy = 1
            if (iz.lt.1) iz = 1
            if (ix.gt.nx) ix = nx
            if (iy.gt.ny) iy = ny
            if (iz.gt.nz) iz = nz
            CALL srcsetup(.false.,
     ;                    xr(ir)-dx, yr(ir), zr(ir), rspread,
     ;                    irreg, um1, ufac, ierr)
            CALl srcsetup(.false.,
     ;                    xr(ir)+dx, yr(ir), zr(ir), rspread,
     ;                    irreg, up1, ufac, ierr)
            dudx = -(up1-um1)/(2.0*dx*cfact)

            CALL srcsetup(.false.,
     ;                    xr(ir), yr(ir)-dy, zr(ir), rspread,
     ;                    irreg, um1, vfac, ierr)
            CALL srcsetup(.false.,
     ;                    xr(ir), yr(ir)+dy, zr(ir), rspread,
     ;                    irreg, up1, vfac, ierr)
            dvdy = -(up1-um1)/(2.0*dy*cfact)

            CALL srcsetup(.false.,
     ;                    xr(ir), yr(ir), zr(ir)-dz, rspread,
     ;                    irreg, um1, wfac, ierr)
            CALL srcsetup(.false.,
     ;                    xr(ir),yr(ir),zr(ir)+dz,rspread,
     ;                    irreg, up1, wfac, ierr)
            dwdz = -(up1-um1)/(2.0*dz*cfact)

            utmp = (da(iz,ix,iy)+2.*mu(iz,ix,iy)/3.)
     ;            *(dudx + dvdy + dwdz)
            utest(iom,ir,isrc) = utmp
         else !set to zero if non-active source-hydrophone pair
            utest(iom,ir,isrc) = CMPLX(0.,0.)
         endif
      enddo

      return
      end

