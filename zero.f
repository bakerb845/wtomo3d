c
      SUBROUTINE zero()

c***********************************************************************
c certain large arrays are zeroed.
c***********************************************************************
c 27 / 11 / 1991 updated by z.song, for 3d case.
c 27/4/1993 updated by z.song.
c 26/5/98 separated from initsub. rgp
c 1/20/10 moved precomp,precons,gp1,gs1 into first loop - baker
c 6/20/16 check if array is allocated and if so null it out - baker

      USE GRADIENTS_MODULE, ONLY : gp1, gs1
      USE PRECON_MODULE, ONLY : preconp, precons
      USE DEMUX_MODULE, ONLY : prof
      USE MAT_MODULE
      USE INTERFACE_MODULE, ONLY : zeroas
      IMPLICIT NONE
!     include 'dimension.inc'
!     include 'init.inc'
!     include 'common.inc'
!     include 'oxford.inc'
!     include 'demux.inc'
      COMPLEX, PARAMETER :: zzero = CMPLX(0.0, 0.0)
      REAL, PARAMETER :: rzero = 0.0

c_______________________________________________________________________
        
c executable code:

c-----------------------------------------------------------------------
c initialize certain arrays to zero:
c-----------------------------------------------------------------------
      IF (ALLOCATED(gp1)) gp1(:,:,:) = zzero
      IF (ALLOCATED(gs1)) gs1(:,:,:) = zzero
      IF (ALLOCATED(preconp)) preconp(:,:,:) = rzero
      IF (ALLOCATED(precons)) precons(:,:,:) = rzero
      IF (ALLOCATED(prof)) prof(:,:) = rzero

!       if (inverse) then 
!         do iy = 1,ny
!         do ix = 1,nx
!           do iz = 1,nz
!             gp1(iz,ix,iy)=CMPLX(0.,0.)
!             gs1(iz,ix,iy)=CMPLX(0.,0.)
!           enddo
!         enddo
!         enddo
!         if (preconspatial.ne.0) then
!           do iy = 1,ny
!           do ix = 1,nx
!             do iz = 1,nz
!               preconp(iz,ix,iy)=0.
!               precons(iz,ix,iy)=0.
!             enddo
!           enddo
!           enddo
!         endif
!       endif

c-----nb: prof is a generic i/o file for profio
!       do j = 1,ng
!         do i = 1,nsam
!           prof(i,j)=0.
!         enddo
!       enddo
      CALL zeroas()

      return
      end
