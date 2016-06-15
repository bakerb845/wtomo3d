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
!TODO: call zeroas()

      a5muv = CMPLX(0.0,0.0)
      a5muu = CMPLX(0.0,0.0)
      a5muw = CMPLX(0.0,0.0)
      a5mvu = CMPLX(0.0,0.0)
      a5mvv = CMPLX(0.0,0.0)
      a5mvw = CMPLX(0.0,0.0)
      a5mwu = CMPLX(0.0,0.0)
      a5mwv = CMPLX(0.0,0.0)
      a5mww = CMPLX(0.0,0.0)

      a5nuv = CMPLX(0.0,0.0)
      a5nuu = CMPLX(0.0,0.0)
      a5nuw = CMPLX(0.0,0.0)
      a5nvu = CMPLX(0.0,0.0)
      a5nvv = CMPLX(0.0,0.0)
      a5nvw = CMPLX(0.0,0.0)
      a5nwu = CMPLX(0.0,0.0)
      a5nwv = CMPLX(0.0,0.0)
      a5nww = CMPLX(0.0,0.0)

      a5puv = CMPLX(0.0,0.0)
      a5puu = CMPLX(0.0,0.0)
      a5puw = CMPLX(0.0,0.0)
      a5pvu = CMPLX(0.0,0.0)
      a5pvv = CMPLX(0.0,0.0)
      a5pvw = CMPLX(0.0,0.0)
      a5pwu = CMPLX(0.0,0.0)
      a5pwv = CMPLX(0.0,0.0)
      a5pww = CMPLX(0.0,0.0)

      a2muv = CMPLX(0.0,0.0)
      a2muu = CMPLX(0.0,0.0)
      a2muw = CMPLX(0.0,0.0)
      a2mvu = CMPLX(0.0,0.0)
      a2mvv = CMPLX(0.0,0.0)
      a2mvw = CMPLX(0.0,0.0)
      a2mwu = CMPLX(0.0,0.0)
      a2mwv = CMPLX(0.0,0.0)
      a2mww = CMPLX(0.0,0.0)

      a2nuv = CMPLX(0.0,0.0)
      a2nuu = CMPLX(0.0,0.0)
      a2nuw = CMPLX(0.0,0.0)
      a2nvu = CMPLX(0.0,0.0)
      a2nvv = CMPLX(0.0,0.0)
      a2nvw = CMPLX(0.0,0.0)
      a2nwu = CMPLX(0.0,0.0)
      a2nwv = CMPLX(0.0,0.0)
      a2nww = CMPLX(0.0,0.0)

      a2puv = CMPLX(0.0,0.0)
      a2puu = CMPLX(0.0,0.0)
      a2puw = CMPLX(0.0,0.0)
      a2pvu = CMPLX(0.0,0.0)
      a2pvv = CMPLX(0.0,0.0)
      a2pvw = CMPLX(0.0,0.0)
      a2pwu = CMPLX(0.0,0.0)
      a2pwv = CMPLX(0.0,0.0)
      a2pww = CMPLX(0.0,0.0)

      a8muv = CMPLX(0.0,0.0)
      a8muu = CMPLX(0.0,0.0)
      a8muw = CMPLX(0.0,0.0)
      a8mvu = CMPLX(0.0,0.0)
      a8mvv = CMPLX(0.0,0.0)
      a8mvw = CMPLX(0.0,0.0)
      a8mwu = CMPLX(0.0,0.0)
      a8mwv = CMPLX(0.0,0.0)
      a8mww = CMPLX(0.0,0.0)

      a8nuv = CMPLX(0.0,0.0)
      a8nuu = CMPLX(0.0,0.0)
      a8nuw = CMPLX(0.0,0.0)
      a8nvu = CMPLX(0.0,0.0)
      a8nvv = CMPLX(0.0,0.0)
      a8nvw = CMPLX(0.0,0.0)
      a8nwu = CMPLX(0.0,0.0)
      a8nwv = CMPLX(0.0,0.0)
      a8nww = CMPLX(0.0,0.0)

      a8puv = CMPLX(0.0,0.0)
      a8puu = CMPLX(0.0,0.0)
      a8puw = CMPLX(0.0,0.0)
      a8pvu = CMPLX(0.0,0.0)
      a8pvv = CMPLX(0.0,0.0)
      a8pvw = CMPLX(0.0,0.0)
      a8pwu = CMPLX(0.0,0.0)
      a8pwv = CMPLX(0.0,0.0)
      a8pww = CMPLX(0.0,0.0)

      a4muv = CMPLX(0.0,0.0)
      a4muu = CMPLX(0.0,0.0)
      a4muw = CMPLX(0.0,0.0)
      a4mvu = CMPLX(0.0,0.0)
      a4mvv = CMPLX(0.0,0.0)
      a4mvw = CMPLX(0.0,0.0)
      a4mwu = CMPLX(0.0,0.0)
      a4mwv = CMPLX(0.0,0.0)
      a4mww = CMPLX(0.0,0.0)

      a4nuv = CMPLX(0.0,0.0)
      a4nuu = CMPLX(0.0,0.0)
      a4nuw = CMPLX(0.0,0.0)
      a4nvu = CMPLX(0.0,0.0)
      a4nvv = CMPLX(0.0,0.0)
      a4nvw = CMPLX(0.0,0.0)
      a4nwu = CMPLX(0.0,0.0)
      a4nwv = CMPLX(0.0,0.0)
      a4nww = CMPLX(0.0,0.0)

      a4puv = CMPLX(0.0,0.0)
      a4puu = CMPLX(0.0,0.0)
      a4puw = CMPLX(0.0,0.0)
      a4pvu = CMPLX(0.0,0.0)
      a4pvv = CMPLX(0.0,0.0)
      a4pvw = CMPLX(0.0,0.0)
      a4pwu = CMPLX(0.0,0.0)
      a4pwv = CMPLX(0.0,0.0)
      a4pww = CMPLX(0.0,0.0)

      a6muv = CMPLX(0.0,0.0)
      a6muu = CMPLX(0.0,0.0)
      a6muw = CMPLX(0.0,0.0)
      a6mvu = CMPLX(0.0,0.0)
      a6mvv = CMPLX(0.0,0.0)
      a6mvw = CMPLX(0.0,0.0)
      a6mwu = CMPLX(0.0,0.0)
      a6mwv = CMPLX(0.0,0.0)
      a6mww = CMPLX(0.0,0.0)

      a6nuv = CMPLX(0.0,0.0)
      a6nuu = CMPLX(0.0,0.0)
      a6nuw = CMPLX(0.0,0.0)
      a6nvu = CMPLX(0.0,0.0)
      a6nvv = CMPLX(0.0,0.0)
      a6nvw = CMPLX(0.0,0.0)
      a6nwu = CMPLX(0.0,0.0)
      a6nwv = CMPLX(0.0,0.0)
      a6nww = CMPLX(0.0,0.0)

      a6puv = CMPLX(0.0,0.0)
      a6puu = CMPLX(0.0,0.0)
      a6puw = CMPLX(0.0,0.0)
      a6pvu = CMPLX(0.0,0.0)
      a6pvv = CMPLX(0.0,0.0)
      a6pvw = CMPLX(0.0,0.0)
      a6pwu = CMPLX(0.0,0.0)
      a6pwv = CMPLX(0.0,0.0)
      a6pww = CMPLX(0.0,0.0)

      a1muv = CMPLX(0.0,0.0)
      a1muu = CMPLX(0.0,0.0)
      a1muw = CMPLX(0.0,0.0)
      a1mvu = CMPLX(0.0,0.0)
      a1mvv = CMPLX(0.0,0.0)
      a1mvw = CMPLX(0.0,0.0)
      a1mwu = CMPLX(0.0,0.0)
      a1mwv = CMPLX(0.0,0.0)
      a1mww = CMPLX(0.0,0.0)

      a1nuv = CMPLX(0.0,0.0)
      a1nuu = CMPLX(0.0,0.0)
      a1nuw = CMPLX(0.0,0.0)
      a1nvu = CMPLX(0.0,0.0)
      a1nvv = CMPLX(0.0,0.0)
      a1nvw = CMPLX(0.0,0.0)
      a1nwu = CMPLX(0.0,0.0)
      a1nwv = CMPLX(0.0,0.0)
      a1nww = CMPLX(0.0,0.0)

      a1puv = CMPLX(0.0,0.0)
      a1puu = CMPLX(0.0,0.0)
      a1puw = CMPLX(0.0,0.0)
      a1pvu = CMPLX(0.0,0.0)
      a1pvv = CMPLX(0.0,0.0)
      a1pvw = CMPLX(0.0,0.0)
      a1pwu = CMPLX(0.0,0.0)
      a1pwv = CMPLX(0.0,0.0)
      a1pww = CMPLX(0.0,0.0)

      a9muv = CMPLX(0.0,0.0)
      a9muu = CMPLX(0.0,0.0)
      a9muw = CMPLX(0.0,0.0)
      a9mvu = CMPLX(0.0,0.0)
      a9mvv = CMPLX(0.0,0.0)
      a9mvw = CMPLX(0.0,0.0)
      a9mwu = CMPLX(0.0,0.0)
      a9mwv = CMPLX(0.0,0.0)
      a9mww = CMPLX(0.0,0.0)

      a9nuv = CMPLX(0.0,0.0)
      a9nuu = CMPLX(0.0,0.0)
      a9nuw = CMPLX(0.0,0.0)
      a9nvu = CMPLX(0.0,0.0)
      a9nvv = CMPLX(0.0,0.0)
      a9nvw = CMPLX(0.0,0.0)
      a9nwu = CMPLX(0.0,0.0)
      a9nwv = CMPLX(0.0,0.0)
      a9nww = CMPLX(0.0,0.0)

      a9puv = CMPLX(0.0,0.0)
      a9puu = CMPLX(0.0,0.0)
      a9puw = CMPLX(0.0,0.0)
      a9pvu = CMPLX(0.0,0.0)
      a9pvv = CMPLX(0.0,0.0)
      a9pvw = CMPLX(0.0,0.0)
      a9pwu = CMPLX(0.0,0.0)
      a9pwv = CMPLX(0.0,0.0)
      a9pww = CMPLX(0.0,0.0)

      a3muv = CMPLX(0.0,0.0)
      a3muu = CMPLX(0.0,0.0)
      a3muw = CMPLX(0.0,0.0)
      a3mvu = CMPLX(0.0,0.0)
      a3mvv = CMPLX(0.0,0.0)
      a3mvw = CMPLX(0.0,0.0)
      a3mwu = CMPLX(0.0,0.0)
      a3mwv = CMPLX(0.0,0.0)
      a3mww = CMPLX(0.0,0.0)

      a3nuv = CMPLX(0.0,0.0)
      a3nuu = CMPLX(0.0,0.0)
      a3nuw = CMPLX(0.0,0.0)
      a3nvu = CMPLX(0.0,0.0)
      a3nvv = CMPLX(0.0,0.0)
      a3nvw = CMPLX(0.0,0.0)
      a3nwu = CMPLX(0.0,0.0)
      a3nwv = CMPLX(0.0,0.0)
      a3nww = CMPLX(0.0,0.0)

      a3puv = CMPLX(0.0,0.0)
      a3puu = CMPLX(0.0,0.0)
      a3puw = CMPLX(0.0,0.0)
      a3pvu = CMPLX(0.0,0.0)
      a3pvv = CMPLX(0.0,0.0)
      a3pvw = CMPLX(0.0,0.0)
      a3pwu = CMPLX(0.0,0.0)
      a3pwv = CMPLX(0.0,0.0)
      a3pww = CMPLX(0.0,0.0)

      a7muv = CMPLX(0.0,0.0)
      a7muu = CMPLX(0.0,0.0)
      a7muw = CMPLX(0.0,0.0)
      a7mvu = CMPLX(0.0,0.0)
      a7mvv = CMPLX(0.0,0.0)
      a7mvw = CMPLX(0.0,0.0)
      a7mwu = CMPLX(0.0,0.0)
      a7mwv = CMPLX(0.0,0.0)
      a7mww = CMPLX(0.0,0.0)

      a7nuv = CMPLX(0.0,0.0)
      a7nuu = CMPLX(0.0,0.0)
      a7nuw = CMPLX(0.0,0.0)
      a7nvu = CMPLX(0.0,0.0)
      a7nvv = CMPLX(0.0,0.0)
      a7nvw = CMPLX(0.0,0.0)
      a7nwu = CMPLX(0.0,0.0)
      a7nwv = CMPLX(0.0,0.0)
      a7nww = CMPLX(0.0,0.0)

      a7puv = CMPLX(0.0,0.0)
      a7puu = CMPLX(0.0,0.0)
      a7puw = CMPLX(0.0,0.0)
      a7pvu = CMPLX(0.0,0.0)
      a7pvv = CMPLX(0.0,0.0)
      a7pvw = CMPLX(0.0,0.0)
      a7pwu = CMPLX(0.0,0.0)
      a7pwv = CMPLX(0.0,0.0)
      a7pww = CMPLX(0.0,0.0)

      return
      end
