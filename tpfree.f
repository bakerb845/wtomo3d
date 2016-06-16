      SUBROUTINE tpfree (lbl, ix, iy, ierr)

      USE COEFFS_MODULE, ONLY : daloc, muloc, roloc,
     ;                          da5, da6, da7, da8,
     ;                          lp2u5, lp2u6, lp2u7, lp2u8,
     ;                          u5, u6, u7, u8
      USE MODEL_MODULE, ONLY : omega, dx, dy, nx, ny
      USE PML_MODULE, ONLY : pmld, pmlf, ipml
      USE MAT_MODULE
      IMPLICIT NONE
c
c     Input
c
c     ibl        is the three letter code for position in the 27 point cube
c     iz,ix,iy   are the model node positions
c

c parameters:
c     include 'dimension.inc'
c arguments:
      INTEGER, INTENT(IN) :: ix, iy
      CHARACTER(3), INTENT(IN) :: lbl
      INTEGER, INTENT(OUT) :: ierr
c common variables:
c     include 'init.inc'
c     include 'common.inc'
c     include 'coeffs.inc'
c     include 'pml.inc'

      COMPLEX r1x, r1y, r1z

      REAL dxo8

      COMPLEX aro, vp
      COMPLEX eyeom
      COMPLEX dnx, dny

      COMPLEX alphax, alphay
      COMPLEX kappax, kappay
      COMPLEX fdom
      REAL    dpmlx, dpmly
      REAL    nodx, nody
      REAL    kmax
      REAL    vpreal


c---this is a trial for kmax; Henrik Bernth suggests using 1 for most cases, and
c---maybe 1.001 for solid-fluid interface.
      ierr = 0
      kmax = 1.0
c     fdom = omega/2.0
c---these should revert back to classic pml
c     kmax = 1.0
      fdom = CMPLX(0., 0.)
c--- this works with kmax = 1.
c       fdom = CMPLX(pi*1.0, 0.)
c--- this also works with kmax = 1, and kmax = 2.
c NB: tests of this setting with Rayleigh waves shows only minor improvement over classic
c PML.   There is still a "dribble" down the side of the model, almost as if there
c were evanescent waves trapped along the edge.

c---reconstruct vp
      aro = roloc(2,2,2)
      vp = (daloc(2,2,2) + 2.0*muloc(2,2,2))/aro
      vp = sqrt(vp)
      vpreal = real(vp)

c----x distance into PML
      dpmlx = 0.0
      if (ix.le.ipml) then
         dpmlx = pmld - (ix-1)*dx
      else if ((nx-ix+1).le.ipml) then
         dpmlx = pmld - (nx-ix)*dx
      endif
c---check to be sure we are in the region
      if (dpmlx.lt.0) then
         write(*,*) ' ***ERROR*** dpmlx = ', dpmlx
         ierr = 1
         return
      endif

      dnx = dpmlx*dpmlx*vpreal*pmlf
      nodx = dpmlx/pmld
      kappax = 1.0 + (kmax - 1)*nodx*nodx
      alphax = fdom*(1.0 - nodx)

c----y distance into PML
      dpmly = 0.0
      if (iy.le.ipml) then
         dpmly = pmld - (iy-1)*dy
      else if ((ny-iy+1).le.ipml) then
         dpmly = pmld - (ny-iy)*dy
      endif
c---check to be sure we are in the region
      if (dpmly.lt.0) then
         write(*,*) ' ***ERROR*** dpmly = ', dpmly
         ierr = 1
         return
      endif

      dny = dpmly*dpmly*vpreal*pmlf
      nody = dpmly/pmld
      kappay = 1.0 + (kmax - 1)*nody*nody
      alphay = fdom*(1.0 - nody)

      eyeom  = CMPLX(0.0, -1.0)*omega

      r1x = (alphax + eyeom)/(kappax*(alphax + eyeom) + dnx)
      r1y = (alphay + eyeom)/(kappay*(alphay + eyeom) + dny)

c--we call this routine only when iz = 1, so r1z is meaningless in this case
c  it is retained here only for asthetics.
      r1z = 1.0

      dxo8 = 1.0/(dx*8.0)

      if (lbl.eq.'adm') then
         a1mUU =  CMPLX(0.0,0.0)
         a1mUV =  CMPLX(0.0,0.0)
         a1mUW =  CMPLX(0.0,0.0)
         a1mVU =  CMPLX(0.0,0.0)
         a1mVV =  CMPLX(0.0,0.0)
         a1mVW =  CMPLX(0.0,0.0)
         a1mWU =  CMPLX(0.0,0.0)
         a1mWV =  CMPLX(0.0,0.0)
         a1mWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'adn') then
         a1nUU =  CMPLX(0.0,0.0)
         a1nUV =  CMPLX(0.0,0.0)
         a1nUW =  CMPLX(0.0,0.0)
         a1nVU =  CMPLX(0.0,0.0)
         a1nVV =  CMPLX(0.0,0.0)
         a1nVW =  CMPLX(0.0,0.0)
         a1nWU =  CMPLX(0.0,0.0)
         a1nWV =  CMPLX(0.0,0.0)
         a1nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'adp') then
         a1pUU =  CMPLX(0.0,0.0)
         a1pUV =  CMPLX(0.0,0.0)
         a1pUW =  CMPLX(0.0,0.0)
         a1pVU =  CMPLX(0.0,0.0)
         a1pVV =  CMPLX(0.0,0.0)
         a1pVW =  CMPLX(0.0,0.0)
         a1pWU =  CMPLX(0.0,0.0)
         a1pWV =  CMPLX(0.0,0.0)
         a1pWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'aam') then
         a2mUU =  CMPLX(0.0,0.0)
         a2mUV =  CMPLX(0.0,0.0)
         a2mUW =  CMPLX(0.0,0.0)
         a2mVU =  CMPLX(0.0,0.0)
         a2mVV =  CMPLX(0.0,0.0)
         a2mVW =  CMPLX(0.0,0.0)
         a2mWU =  CMPLX(0.0,0.0)
         a2mWV =  CMPLX(0.0,0.0)
         a2mWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'aan') then
         a2nUU =  CMPLX(0.0,0.0)
         a2nUV =  CMPLX(0.0,0.0)
         a2nUW =  CMPLX(0.0,0.0)
         a2nVU =  CMPLX(0.0,0.0)
         a2nVV =  CMPLX(0.0,0.0)
         a2nVW =  CMPLX(0.0,0.0)
         a2nWU =  CMPLX(0.0,0.0)
         a2nWV =  CMPLX(0.0,0.0)
         a2nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'aap') then
         a2pUU =  CMPLX(0.0,0.0)
         a2pUV =  CMPLX(0.0,0.0)
         a2pUW =  CMPLX(0.0,0.0)
         a2pVU =  CMPLX(0.0,0.0)
         a2pVV =  CMPLX(0.0,0.0)
         a2pVW =  CMPLX(0.0,0.0)
         a2pWU =  CMPLX(0.0,0.0)
         a2pWV =  CMPLX(0.0,0.0)
         a2pWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'afm') then
         a3mUU =  r1z*dxo8*u5
         a3mUV =  CMPLX(0.0,0.0)
         a3mUW = -r1x*dxo8*u5
         a3mVU =  CMPLX(0.0,0.0)
         a3mVV =  r1z*dxo8*u5
         a3mVW = -r1y*dxo8*u5
         a3mWU = -r1x*dxo8*da5
         a3mWV = -r1y*dxo8*da5
         a3mWW =  r1z*dxo8*lp2u5

      elseif (lbl.eq.'afn') then
         a3nUU =  CMPLX(0.0,0.0)
         a3nUV =  CMPLX(0.0,0.0)
         a3nUW =  CMPLX(0.0,0.0)
         a3nVU =  CMPLX(0.0,0.0)
         a3nVV =  CMPLX(0.0,0.0)
         a3nVW =  CMPLX(0.0,0.0)
         a3nWU =  CMPLX(0.0,0.0)
         a3nWV =  CMPLX(0.0,0.0)
         a3nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'afp') then
         a3pUU =  r1z*dxo8*u6
         a3pUV =  CMPLX(0.0,0.0)
         a3pUW = -r1x*dxo8*u6
         a3pVU =  CMPLX(0.0,0.0)
         a3pVV =  r1z*dxo8*u6
         a3pVW =  r1y*dxo8*u6
         a3pWU = -r1x*dxo8*da6
         a3pWV =  r1y*dxo8*da6
         a3pWW =  r1z*dxo8*lp2u6

      elseif (lbl.eq.'ddm') then
         a4mUU =  CMPLX(0.0,0.0)
         a4mUV =  CMPLX(0.0,0.0)
         a4mUW =  CMPLX(0.0,0.0)
         a4mVU =  CMPLX(0.0,0.0)
         a4mVV =  CMPLX(0.0,0.0)
         a4mVW =  CMPLX(0.0,0.0)
         a4mWU =  CMPLX(0.0,0.0)
         a4mWV =  CMPLX(0.0,0.0)
         a4mWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ddn') then
         a4nUU =  CMPLX(0.0,0.0)
         a4nUV =  CMPLX(0.0,0.0)
         a4nUW =  CMPLX(0.0,0.0)
         a4nVU =  CMPLX(0.0,0.0)
         a4nVV =  CMPLX(0.0,0.0)
         a4nVW =  CMPLX(0.0,0.0)
         a4nWU =  CMPLX(0.0,0.0)
         a4nWV =  CMPLX(0.0,0.0)
         a4nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ddp') then
         a4pUU =  CMPLX(0.0,0.0)
         a4pUV =  CMPLX(0.0,0.0)
         a4pUW =  CMPLX(0.0,0.0)
         a4pVU =  CMPLX(0.0,0.0)
         a4pVV =  CMPLX(0.0,0.0)
         a4pVW =  CMPLX(0.0,0.0)
         a4pWU =  CMPLX(0.0,0.0)
         a4pWV =  CMPLX(0.0,0.0)
         a4pWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'bem') then
         a5mUU =  CMPLX(0.0,0.0)
         a5mUV =  CMPLX(0.0,0.0)
         a5mUW =  CMPLX(0.0,0.0)
         a5mVU =  CMPLX(0.0,0.0)
         a5mVV =  CMPLX(0.0,0.0)
         a5mVW =  CMPLX(0.0,0.0)
         a5mWU =  CMPLX(0.0,0.0)
         a5mWV =  CMPLX(0.0,0.0)
         a5mWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ben') then
         a5nUU = -r1z*dxo8*(u8 + u7 + u6 + u5)
         a5nUV =  CMPLX(0.0,0.0)
         a5nUW =  r1x*dxo8*(u5 - u8 + u6 - u7)
         a5nVU =  CMPLX(0.0,0.0)
         a5nVV = -r1z*dxo8*(u8 + u7 + u6 + u5)
         a5nVW =  r1y*dxo8*(da5 - da8 - da6 + da7)
         a5nWU =  r1x*dxo8*(da5 - da8 + da6 - da7)
         a5nWV =  r1y*dxo8*(da5 - da8 - da6 + da7)
         a5nWW = -r1z*dxo8*(lp2u8 + lp2u7 + lp2u6 + lp2u5)

      elseif (lbl.eq.'bep') then
         a5pUU =  CMPLX(0.0,0.0)
         a5pUV =  CMPLX(0.0,0.0)
         a5pUW =  CMPLX(0.0,0.0)
         a5pVU =  CMPLX(0.0,0.0)
         a5pVV =  CMPLX(0.0,0.0)
         a5pVW =  CMPLX(0.0,0.0)
         a5pWU =  CMPLX(0.0,0.0)
         a5pWV =  CMPLX(0.0,0.0)
         a5pWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ffm') then
         a6mUU =  CMPLX(0.0,0.0)
         a6mUV =  CMPLX(0.0,0.0)
         a6mUW =  CMPLX(0.0,0.0)
         a6mVU =  CMPLX(0.0,0.0)
         a6mVV =  CMPLX(0.0,0.0)
         a6mVW =  CMPLX(0.0,0.0)
         a6mWU =  CMPLX(0.0,0.0)
         a6mWV =  CMPLX(0.0,0.0)
         a6mWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ffn') then
         a6nUU =  CMPLX(0.0,0.0)
         a6nUV =  CMPLX(0.0,0.0)
         a6nUW =  CMPLX(0.0,0.0)
         a6nVU =  CMPLX(0.0,0.0)
         a6nVV =  CMPLX(0.0,0.0)
         a6nVW =  CMPLX(0.0,0.0)
         a6nWU =  CMPLX(0.0,0.0)
         a6nWV =  CMPLX(0.0,0.0)
         a6nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ffp') then
         a6pUU =  CMPLX(0.0,0.0)
         a6pUV =  CMPLX(0.0,0.0)
         a6pUW =  CMPLX(0.0,0.0)
         a6pVU =  CMPLX(0.0,0.0)
         a6pVV =  CMPLX(0.0,0.0)
         a6pVW =  CMPLX(0.0,0.0)
         a6pWU =  CMPLX(0.0,0.0)
         a6pWV =  CMPLX(0.0,0.0)
         a6pWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'cdm') then
         a7mUU =  CMPLX(0.0,0.0)
         a7mUV =  CMPLX(0.0,0.0)
         a7mUW =  CMPLX(0.0,0.0)
         a7mVU =  CMPLX(0.0,0.0)
         a7mVV =  CMPLX(0.0,0.0)
         a7mVW =  CMPLX(0.0,0.0)
         a7mWU =  CMPLX(0.0,0.0)
         a7mWV =  CMPLX(0.0,0.0)
         a7mWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'cdn') then
         a7nUU =  CMPLX(0.0,0.0)
         a7nUV =  CMPLX(0.0,0.0)
         a7nUW =  CMPLX(0.0,0.0)
         a7nVU =  CMPLX(0.0,0.0)
         a7nVV =  CMPLX(0.0,0.0)
         a7nVW =  CMPLX(0.0,0.0)
         a7nWU =  CMPLX(0.0,0.0)
         a7nWV =  CMPLX(0.0,0.0)
         a7nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'cdp') then
         a7pUU =  CMPLX(0.0,0.0)
         a7pUV =  CMPLX(0.0,0.0)
         a7pUW =  CMPLX(0.0,0.0)
         a7pVU =  CMPLX(0.0,0.0)
         a7pVV =  CMPLX(0.0,0.0)
         a7pVW =  CMPLX(0.0,0.0)
         a7pWU =  CMPLX(0.0,0.0)
         a7pWV =  CMPLX(0.0,0.0)
         a7pWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ccm') then
         a8mUU =  CMPLX(0.0,0.0)
         a8mUV =  CMPLX(0.0,0.0)
         a8mUW =  CMPLX(0.0,0.0)
         a8mVU =  CMPLX(0.0,0.0)
         a8mVV =  CMPLX(0.0,0.0)
         a8mVW =  CMPLX(0.0,0.0)
         a8mWU =  CMPLX(0.0,0.0)
         a8mWV =  CMPLX(0.0,0.0)
         a8mWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ccn') then
         a8nUU =  CMPLX(0.0,0.0)
         a8nUV =  CMPLX(0.0,0.0)
         a8nUW =  CMPLX(0.0,0.0)
         a8nVU =  CMPLX(0.0,0.0)
         a8nVV =  CMPLX(0.0,0.0)
         a8nVW =  CMPLX(0.0,0.0)
         a8nWU =  CMPLX(0.0,0.0)
         a8nWV =  CMPLX(0.0,0.0)
         a8nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'ccp') then
         a8pUU =  CMPLX(0.0,0.0)
         a8pUV =  CMPLX(0.0,0.0)
         a8pUW =  CMPLX(0.0,0.0)
         a8pVU =  CMPLX(0.0,0.0)
         a8pVV =  CMPLX(0.0,0.0)
         a8pVW =  CMPLX(0.0,0.0)
         a8pWU =  CMPLX(0.0,0.0)
         a8pWV =  CMPLX(0.0,0.0)
         a8pWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'cfm') then
         a9mUU =  r1z*dxo8*u7
         a9mUV =  CMPLX(0.0,0.0)
         a9mUW =  r1x*dxo8*u7
         a9mVU =  CMPLX(0.0,0.0)
         a9mVV =  r1z*dxo8*u7
         a9mVW = -r1y*dxo8*u7
         a9mWU =  r1x*dxo8*da7
         a9mWV = -r1y*dxo8*da7
         a9mWW =  r1z*dxo8*lp2u7

      elseif (lbl.eq.'cfn') then
         a9nUU =  CMPLX(0.0,0.0)
         a9nUV =  CMPLX(0.0,0.0)
         a9nUW =  CMPLX(0.0,0.0)
         a9nVU =  CMPLX(0.0,0.0)
         a9nVV =  CMPLX(0.0,0.0)
         a9nVW =  CMPLX(0.0,0.0)
         a9nWU =  CMPLX(0.0,0.0)
         a9nWV =  CMPLX(0.0,0.0)
         a9nWW =  CMPLX(0.0,0.0)

      elseif (lbl.eq.'cfp') then
         a9pUU =  r1z*dxo8*u8
         a9pUV =  CMPLX(0.0,0.0)
         a9pUW =  r1x*dxo8*u8
         a9pVU =  CMPLX(0.0,0.0)
         a9pVV =  r1z*dxo8*u8
         a9pVW =  r1y*dxo8*u8
         a9pWU =  r1x*dxo8*da8
         a9pWV =  r1y*dxo8*da8
         a9pWW =  r1z*dxo8*lp2u8

      endif


      return
      end

