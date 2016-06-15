      SUBROUTINE srcsetup(lsetup,
     ;                    xloc, yloc, zloc, sprd, ireg,
     ;                    sterm, swave, ierr)

c**********************************************************************
c Setup the source/receiver region and either enter a source term into 
c the modelling grid or extract a pressure value at a receiver from a 
c wavefield.
c
c If lsetup=T : setup source term vector -> swave(nz,nx,ny)
c If lsetup=F : extract wavefiled value at receiver -> sterm
c
c
c New option added if ireg=0.
c Does the same as for ireg>0 but all the sources/receivers are treated
c as POINT sources/receivers LOCATED AT THE NEAREST MODELLING NODE.
c If ireg>0 this subroutine treats sources/receivers located between
c nodes correctly (and enables the use of sources/receivers with a 
c finite spatial extent) but is quite slow when ns and nr are large.
c NB This subroutine is called ns*nr times during the computation of the 
c backpropagated wavefields.      G.Hicks, April 2001.
c**********************************************************************
      USE INIT_MODULE, ONLY : freesurf, dx, dy, dz, nx, ny, nz
      USE INTERFACE_MODULE, ONLY : srcreg
      IMPLICIT NONE

c       include 'dimension.inc'
c       include 'init.inc'

      LOGICAL, INTENT(IN) ::  lsetup
      REAL, INTENT(IN) :: xloc, yloc, zloc, sprd
      INTEGER, INTENT(IN) :: ireg
      COMPLEX, DIMENSION(:,:,:), INTENT(INOUT) :: swave
      COMPLEX, INTENT(INOUT) :: sterm
      INTEGER, INTENT(OUT) :: ierr

c local variables:
      REAL fs(9, 9, 9)
      INTEGER  ixs,iys,izs,ixmaxs,ixmins,iymaxs,iymins,izmaxs,izmins
      INTEGER  ixmax,ixmin,iymax,iymin,izmax,izmin,izsrg,iysrg,ixsrg
      INTEGER  ix,iy,iz,ixfs,iyfs,izfs
      REAL     xoff, yoff, zoff
      INTEGER, PARAMETER :: nrgg = 9

c       print*,'DEBUG lsetup,xloc,zloc,sprd,ireg',
c     +   lsetup,xloc,zloc,sprd,ireg
      ierr = 0
      fs(:,:,:) = 0.0

      if (ireg.eq.0) then

C Source location is now specified as a real number - 
C convert this to THE NEAREST grid location:
C
         ixs = INT( 0.5 + (xloc/dx) ) + 1
         iys = INT( 0.5 + (yloc/dy) ) + 1
         izs = INT( 0.5 + (zloc/dz) ) + 1

C
C Setup source term or find wavefield value at current reciever:
C
         if (lsetup) then
            do iy=1,ny
               do ix=1,nx
                  do iz=1,nz
                     swave(iz,ix,iy)=CMPLX(0.,0.)
                  enddo
               enddo
            enddo
            swave(izs,ixs,iys) = sterm
         else
            sterm = swave(izs,ixs,iys)
         endif

c       write(*,*) ' ixs, iys, izs = ', ixs, iys, izs
c       write(*,*) ' sterm = ', sterm

      else
c
C Source location is now specified as a real number - 
C convert this to a grid location (top left):
C
         ixs = INT(xloc/dx) + 1
         iys = INT(yloc/dy) + 1
         izs = INT(zloc/dz) + 1
         xoff = xloc - (FLOAT(ixs-1)*dx)
         yoff = yloc - (FLOAT(iys-1)*dy)
         zoff = zloc - (FLOAT(izs-1)*dz)

c       write(*,*) ' ixs, iys, izs = ', ixs, iys, izs
c       write(*,*) ' sterm = ', sterm
C
C Setup source region diemsions:
C
         if (sprd.le.(0.5)) then
            ixmaxs = ixs + ireg
            ixmins = ixs - ireg
            iymaxs = iys + ireg
            iymins = iys - ireg
            izmaxs = izs + ireg
            izmins = izs - ireg
         else
            ixmaxs = ixs + INT(2*sprd)
            ixmins = ixs - INT(2*sprd)
            iymaxs = iys + INT(2*sprd)
            iymins = iys - INT(2*sprd)
            izmaxs = izs + INT(2*sprd)
            izmins = izs - INT(2*sprd)
         endif
C
C Compute the spatial source distribution, fs:
C
         CALL srcreg (fs, ireg, sprd, dx, dy, dz,
     ;                xoff, yoff, zoff, ierr)
         IF (ierr /= 0) RETURN
C
C Clear initial term(s) to zero:
C
         if (lsetup) then
            do iy = 1, ny
               do ix = 1, nx
                  do iz = 1, nz
                     swave(iz,ix,iy) = CMPLX(0.,0.)
                  enddo
              enddo
            enddo
         else
            sterm = CMPLX(0.,0.)
         endif
C
C If part of the source region goes out of the model, then:
C if an absorbing boundary - ignore,
C if above a free surface - reflect source terms back inside the model
C with a change of sign.
C
         izmin=izmins
         izmax=izmaxs
         iymin=iymins
         iymax=iymaxs
         ixmin=ixmins
         ixmax=ixmaxs
         if (izmin.lt.1.and.(.not.freesurf(1))) izmin=1
         if (izmax.gt.nz.and.(.not.freesurf(3))) izmax=nz
         if (ixmin.lt.1.and.(.not.freesurf(4))) ixmin=1
         if (ixmax.gt.nx.and.(.not.freesurf(2))) ixmax=nx
         if (iymin.lt.1.and.(.not.freesurf(6))) iymin=1
         if (iymax.gt.ny.and.(.not.freesurf(5))) iymax=ny

         do iz = izmin, izmax
            izsrg = iz - izmins + 1

            if (iz.lt.1) then
               izfs = 2-iz
            else if (iz.gt.nz) then
               izfs = (2*nz)-iz
            else
               izfs = iz
            endif

            do ix = ixmin, ixmax
               ixsrg = ix - ixmins + 1

               if (ix.lt.1) then
                  ixfs = 2-ix
               else if (ix.gt.nx) then
                  ixfs = (2*nx)-ix
               else
                  ixfs = ix
               endif

               do iy = iymin, iymax
                  iysrg = iy - iymins + 1

                  if (iy.lt.1) then
                     iyfs = 2-iy
                  else if (iy.gt.ny) then
                     iyfs = (2*ny)-iy
                  else
                     iyfs = iy
                  endif

                  if (lsetup) then

                     if (iz.eq.izfs.and.ix.eq.ixfs.and.iy.eq.iyfs) then
                         swave(izfs,ixfs,iyfs) = swave(izfs,ixfs,iyfs) +
     &                          fs(izsrg,ixsrg,iysrg)*sterm
                     else
                       swave(izfs,ixfs,iyfs) = swave(izfs,ixfs,iyfs) - 
     &                            fs(izsrg,ixsrg,iysrg)*sterm
                     endif

                  else

                     if (iz.eq.izfs.and.ix.eq.ixfs.and.iy.eq.iyfs) then
                        sterm = sterm
     ;                    + fs(izsrg,ixsrg,iysrg)*swave(izfs,ixfs,iyfs)
                     else
                        sterm = sterm
     ;                    - fs(izsrg,ixsrg,iysrg)*swave(izfs,ixfs,iyfs)
                     endif

c       if (ix.eq.ixmin.and.iy.eq.iymin.and.iz.eq.izmin) then
c       write(*,*) " in srcsetup, sterm, fs, swave = ", sterm, fs(izsrg, ixsrg, iysrg), swave(izfs, ixfs, iysrg)
c       endif

                  endif

               enddo
            enddo
         enddo
      endif

      return
      end
