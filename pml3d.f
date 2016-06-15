      SUBROUTINE pml3d (lbl, iz, ix, iy, ierr)

c----- This routine generates the entries on the 3D elastic FDFD system of equations.
c       S. Roecker 2/2012
c
c
c       Input:
c       iz, iy, ix  are the model indices for the point in question
c                       Note that these indices are used only for finding the position
c                       in the PML, so if any of them are outide the PML then there is no associated
c                       decay in that particular direction.
c       lbl  jis a three letter code to specify which of the 27 elements is desired
c
c       Forms for the BcU system
c       Conventions
c
c       The indices on Lame parameters are
c       m    -1
c       mh -1/2
c       n     0 (null)
c       ph +1/2
c       p    +1
c
c       and are in order of i,j,k.   So, for example
c
c       lp2upphn        would be (lambda +2 mu) (i+1, j+1/2, k)
c
c       For the matrix entries in the 27 point star we adopt the taxonomy used
c       in the 2D problem and extend to 3D by adding an additional letter (m,n,p).
c       We retain the x-z meaning of the original 2D frame, and add (m,n,p) to
c       indicate the index on j.   For example, the standard 2D x-z 9 star frame would be
c
c               a1n  a4n  a7n           adn  ddn  cdn
c               a2n  a5n  a8n     or    aan  ben  ccn
c               a3n  a6n  a9n           afn  ffn  cfn
c
c       for j+1 we would use "p" and for j-1 we use "m". For example, the point
c       uu(i+1, j-1, k) = a8puu
c
c       Weights:
c
c       We use three systems and mass weighting, and weight the contributions of each
c       to the total. The system weights are w1, w2, and w3; the mass dump weights are
c       wm1, wm2, wm3, and wm4.
c
c       Mass weighting:
c
c       w^2pU = w^2*(wm1*nnn*a5n 
c                 +  wm2*(pnn*a8n + mnn*a2n + npn*a5p + nmn*a5m + nnp*a6n + nnm*a4n)/6.0
c                 +  wm3*(ppn*a8p + mpn*a2p + pmn*a8m + mmn*a2m 
c                       + pnp*a9n + mnp*a3n + pnm*a7n + mnm*a1n 
c                       + npp*a6p + nmp*a6m + npm*a4p + nmm*a4m)/12.0
c                 +  wm4*(ppp*a9p + ppm*a7p + pmp*a9m + pmm*a7m 
c                       + mpp*a3p + mpm*a1p + mmp*a3m + mmm*a1m)/8.0)
c       Bring the densities in a roloc(3,3,3) array,  1 = n, 2 = m, 3 = p.   Note the ordering
c       above is (ix, iy, iz), so, for example npm = (ix, iy+1, iz-1)
c
c
      USE INIT_MODULE, ONLY : dx, nx, ny, nz
      USE COEFFS_MODULE, ONLY : daloc, muloc, roloc,
     ;                          lmu1, lmu2, lmu3, lmu4,
     ;                          lmu5, lmu6, lmu7, lmu8,
     ;                          lpu1, lpu2, lpu3, lpu4,
     ;                          lpu5, lpu6, lpu7, lpu8,
     ;                          lp2u1, lp2u2, lp2u3, lp2u4, lp2u5,
     ;                          lp2u6, lp2u7, lp2u8,
     ;                          u1, u2, u3, u4, u5, u6, u7, u8,
     ;                          wm1, wm2, wm3, wm4
      USE MODEL_MODULE, ONLY : omega
      USE PML_MODULE, ONLY : pmld, pmlf, ipml
      USE MAT_MODULE
      IMPLICIT NONE


c parameters:
c       include 'dimension.inc'
c
c arguments:
      INTEGER, INTENT(IN) :: ix, iy, iz
      CHARACTER(3), INTENT(IN) :: lbl
      INTEGER, INTENT(OUT) :: ierr

c local variables:
      REAL    dx2, dxx, dxx2, dxx4
      REAL    dpmlx, dpmly, dpmlz
      REAL    nodx, nody, nodz
      REAL    kmax
      REAL    vpreal

      REAL    wm2o6, wm3o12, wm4o8

c----Note:  aro and vp are declared complex because they are
c    involved with complex operations, but the density (aro) is
c    always real and vp SHOULD be real for the terms the use it.
c    If we have finite Q, however, vp will be complex, and we
c    force it to be real.
      COMPLEX aro, vp   
      COMPLEX eye, eyeom
      COMPLEX dnx, dny, dnz
      COMPLEX r1x, r1xsq, r2x
      COMPLEX r1y, r1ysq, r2y
      COMPLEX r1z, r1zsq, r2z
      COMPLEX r2xo2d, r1xqodq, r1xqo2dq
      COMPLEX r2yo2d, r1yqodq, r1yqo2dq
      COMPLEX r2zo2d, r1zqodq, r1zqo2dq

      COMPLEX r1xr1yo4dq, r1xr1zo4dq, r1yr1zo4dq
      COMPLEX r2xpr2y, r2xmr2y
      COMPLEX r2xpr2z, r2xmr2z
      COMPLEX r2ypr2z, r2ymr2z

      COMPLEX r1xpr1y, r1xpr1z, r1ypr1z

      COMPLEX alphax, alphay, alphaz
      COMPLEX kappax, kappay, kappaz
      COMPLEX fdom
      COMPLEX omegasq
      COMPLEX massf
      COMPLEX lsum, usum

      INTEGER isnx, isny, isnz

c common variables:

c       include 'common.inc'
c       include 'init.inc'
c       include 'coeffs.inc'
c       include 'pml.inc'

c  Contents of pml.inc:
c     pmld    width of PML layer in meters
c     pmlr    desired reflection coefficient (nominally 0.001)
c     pmlf    3/(2*pmld^3)*log(1/pmlr)
c     pmlf is computed in initio
c      
c     Thus:   
c      d(n)  = n*n*vp*pmlf
c      d'(n) = 2*n*vp*pmlf
c_______________________________________________________________________
c
c 9-point fd star  - Values estimated from the graph published by 
c Stekl and Pratt.  It would be better to produce determine these as a 
c function of Poisson's ratio although it is more important for acoef 
c than bcoef
c
c---this is now computed in setcoefs, and passed through coeffs.inc
c     w1 = 0.1
c     w2 = 0.1
c     w3 = 1.0 - w1 - w3
c     wm1 = 0.2
c     wm2 = 0.2
c     wm3 = 0.2
c     wm4 = 1.0 - (wm1 + wm2 + wm3)

c---turn off mass lumping for a test
c       wm1 = 1.0
c       wm2 = 0.0
c       wm3 = 0.0
c       wm4 = 0.0
      ierr = 0
c---this is a trial for kmax; Henrik Bernth suggests using 1 for most cases, and
c---maybe 1.001 for solid-fluid interface. 
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
c

          
c executable code:
      eye    = CMPLX(0.0, -1.0)
      eyeom  = CMPLX(0.0, -1.0)*omega

      omegasq = omega * omega

c TODO: check this
      wm2o6  = REAL(omegasq*wm2/6.0)
      wm3o12 = REAL(omegasq*wm3/12.0)
      wm4o8  = REAL(omegasq*wm4/8.0)

      dx2 = dx + dx
      dxx = dx*dx
      dxx2 = dxx + dxx
      dxx4 = dxx2 + dxx2

c---reconstruct vp
      aro = roloc(2,2,2)
      vp = (daloc(2,2,2) + 2.0*muloc(2,2,2))/aro
      vp = CSQRT(vp)
      vpreal = REAL(vp)

c----x distance into PML
      dpmlx = 0.0
      if (ix.le.ipml) then
         dpmlx = pmld - (ix-1)*dx
         isnx = 1
      else if ((nx-ix+1).le.ipml) then
         dpmlx = pmld - (nx-ix)*dx
         isnx = -1
      endif
c---check to be sure we are in the region
      if (dpmlx.lt.0) then
         write(*,*) ' ***ERROR*** dpmlx = ', dpmlx
         ierr = 1 
         RETURN 
      endif

      dnx = dpmlx*dpmlx*vpreal*pmlf
      nodx = dpmlx/pmld
      kappax = 1.0 + (kmax - 1)*nodx*nodx
      alphax = fdom*(1.0 - nodx)

c----y distance into PML
      dpmly = 0.0
      if (iy.le.ipml) then
         dpmly = pmld - (iy-1)*dx
         isny = 1
      else if ((ny-iy+1).le.ipml) then
         dpmly = pmld - (ny-iy)*dx
         isny = -1
      endif
c---check to be sure we are in the region
      if (dpmly.lt.0) then
         write(*,*) ' ***ERROR*** dpmly = ', dpmly
         ierr = 1 
         RETURN 
      endif

      dny = dpmly*dpmly*vpreal*pmlf
      nody = dpmly/pmld
      kappay = 1.0 + (kmax - 1)*nody*nody
      alphay = fdom*(1.0 - nody)

c----z distance into PML
      dpmlz = 0.0
      if (iz.le.ipml) then
         dpmlz = pmld - (iz-1)*dx
         isnz = 1
      else if ((nz-iz+1).le.ipml) then
         dpmlz = pmld - (nz-iz)*dx
         isnz = -1
      endif
c---check to be sure we are in the region
      if (dpmlz.lt.0) then
         write(*,*) ' ***ERROR*** dpmlz = ', dpmlz
         ierr = 1 
         RETURN 
      endif

      dnz = dpmlz*dpmlz*vpreal*pmlf
      nodz = dpmlz/pmld
      kappaz = 1.0 + (kmax - 1)*nodz*nodz
      alphaz = fdom*(1.0 - nodz)

      r1x = (alphax + eyeom)/(kappax*(alphax + eyeom) + dnx)
      r1xsq = r1x*r1x
c--build r2x
      r2x = fdom*(1.-nodx) + eyeom
      r2x = (fdom*(1.-0.5*nodx) + eyeom)/(r2x*r2x)
      r2x = isnx*2.0*r1xsq*r1x*dpmlx*((kmax-1)/(pmld*pmld)
     ;    + vpreal*pmlf*r2x)

      r1y = (alphay + eyeom)/(kappay*(alphay + eyeom) + dny)
      r1ysq = r1y*r1y
c--build r2y
      r2y = fdom*(1.-nody) + eyeom
      r2y = (fdom*(1.-0.5*nody) + eyeom)/(r2y*r2y)
      r2y = isny*2.0*r1ysq*r1y*dpmly*((kmax-1)/(pmld*pmld)
     ;    + vpreal*pmlf*r2y)

      r1z = (alphaz + eyeom)/(kappaz*(alphaz + eyeom) + dnz)
      r1zsq = r1z*r1z
c--build r2z
      r2z = fdom*(1.-nodz) + eyeom
      r2z = (fdom*(1.-0.5*nodz) + eyeom)/(r2z*r2z)
      r2z = isnz*2.0*r1zsq*r1z*dpmlz*((kmax-1)/(pmld*pmld)
     ;    + vpreal*pmlf*r2z)

c---for a test of symmetry, turn off the r2's
c       r2x = CMPLX(0.0,0.0)
c       r2y = CMPLX(0.0,0.0)
c       r2z = CMPLX(0.0,0.0)
c---for a test of symmetry, turn off the r1's
c       r1x = CMPLX(0.0,0.0)
c       r1y = CMPLX(0.0,0.0)
c       r1z = CMPLX(0.0,0.0)

      r2xo2d  = r2x/dx2
      r2yo2d  = r2y/dx2
      r2zo2d  = r2z/dx2

      r1xqodq = (r1x*r1x)/dxx
      r1yqodq = (r1y*r1y)/dxx
      r1zqodq = (r1z*r1z)/dxx

      r1xr1yo4dq = r1x*r1y/dxx4
      r1xr1zo4dq = r1x*r1z/dxx4
      r1yr1zo4dq = r1y*r1z/dxx4

c---Divide the above by 4 because B1-4 is an average over 4 systems
      r2xo2d  = r2xo2d/4.0
      r2yo2d  = r2yo2d/4.0
      r2zo2d  = r2zo2d/4.0
      r1xqodq = r1xqodq/4.0
      r1yqodq = r1yqodq/4.0
      r1zqodq = r1zqodq/4.0 
      r1xr1yo4dq = r1xr1yo4dq/4.0
      r1xr1zo4dq = r1xr1zo4dq/4.0
      r1yr1zo4dq = r1yr1zo4dq/4.0

      r1xqo2dq = r1xqodq/2.0
      r1yqo2dq = r1yqodq/2.0
      r1zqo2dq = r1zqodq/2.0

      r2xpr2y = r2xo2d + r2yo2d
      r2xmr2y = r2xo2d - r2yo2d
      r2xpr2z = r2xo2d + r2zo2d
      r2xmr2z = r2xo2d - r2zo2d
      r2ypr2z = r2yo2d + r2zo2d
      r2ymr2z = r2yo2d - r2zo2d

      r1xpr1y = r1xqo2dq + r1yqo2dq
      r1xpr1z = r1xqo2dq + r1zqo2dq
      r1ypr1z = r1yqo2dq + r1zqo2dq

      if (lbl.eq.'adm') then

         massf = wm4o8*roloc(1,1,1)

         a1mUU =  r1xqo2dq*lp2u1 + r1ypr1z*u1 - r2xo2d*lp2u1
     ;         - r2ypr2z*u1 + massf
         a1mVV =  r1yqo2dq*lp2u1 + r1xpr1z*u1 - r2yo2d*lp2u1
     ;         - r2xpr2z*u1 + massf
         a1mWW =  r1zqo2dq*lp2u1 + r1xpr1y*u1 - r2zo2d*lp2u1
     ;         - r2xpr2y*u1 + massf

         a1mUV =  r1xr1yo4dq*lpu1
         a1mVU =  a1mUV
         a1mUW =  r1xr1zo4dq*lpu1
         a1mWU =  a1mUW
         a1mVW =  r1yr1zo4dq*lpu1
         a1mWV =  a1mVW

      else if (lbl.eq.'adn') then

         massf = wm3o12*roloc(1,2,1)

         a1nUU = -r1yqo2dq*(u1 + u2) + massf
         a1nVV =  a1nUU
         a1nWW =  a1nUU

         a1nUV =  r1xr1yo4dq*(lmu2 - lmu1)
         a1nVU = -a1nUV

         a1nUW =  r1xr1zo4dq*(lpu1 + lpu2)
         a1nWU =  a1nUW

         a1nVW =  r1yr1zo4dq*(lmu1 - lmu2)
         a1nWV = -a1nVW

      else if (lbl.eq.'adp') then

         massf = wm4o8*roloc(1,3,1)

         a1pUU =  r1xqo2dq*lp2u2 + r1ypr1z*u2 - r2xo2d*lp2u2
     ;         + r2ymr2z*u2 + massf
         a1pVV =  r1yqo2dq*lp2u2 + r1xpr1z*u2 + r2yo2d*lp2u2
     ;         - r2xpr2z*u2 + massf
         a1pWW =  r1zqo2dq*lp2u2 + r1xpr1y*u2 - r2zo2d*lp2u2
     ;         - r2xmr2y*u2 + massf

         a1pUV = -r1xr1yo4dq*lpu2
         a1pVU =  a1pUV
         a1pUW =  r1xr1zo4dq*lpu2
         a1pWU =  a1pUW
         a1pVW = -r1yr1zo4dq*lpu2
         a1pWV =  a1pVW

      else if (lbl.eq.'aam') then

         massf = wm3o12*roloc(1,1,2)

         a2mUU = -r1zqo2dq*(u1 + u5) + massf
         a2mVV =  a2mUU
         a2mWW = -r1zqo2dq*(lp2u1 + lp2u5) + massf

         a2mUV =  r1xr1yo4dq*(lpu1 + lpu5)
         a2mVU =  a2mUV
         a2mUW = -r1xr1zo4dq*(lmu1 - lmu5)
         a2mWU = -a2mUW
         a2mVW = -r1yr1zo4dq*(lmu1 - lmu5)
         a2mWV = -a2mVW

      else if (lbl.eq.'aan') then

         massf = wm2o6*roloc(1,2,2)

         a2nUU =  r1xqo2dq*(lp2u1 + lp2u6 + lp2u2 + lp2u5) + massf
         a2nVV =  a2nUU
         a2nWW =  r1xqo2dq*(u1 + u6 + u2 + u5)             + massf

         a2nUV =  r1xr1yo4dq*(lmu2 - lmu5 + lmu6 - lmu1)
         a2nVU = -a2nUV
         a2nUW = -r1xr1zo4dq*(lmu2 - lmu5 + lmu1 - lmu6)
         a2nWU = -a2nUW
         a2nVW =  r1yr1zo4dq*(lpu2 + lpu5 - lpu1 - lpu6)
         a2nWV =  a2nVW

      else if (lbl.eq.'aap') then

         massf = wm3o12*roloc(1,3,2)

         a2pUU = -r1zqo2dq*(u2 + u6) + massf
         a2pVV =  a2pUU
         a2pWW = -r1zqo2dq*(lp2u2 + lp2u6) + massf

         a2pUV = -r1xr1yo4dq*(lpu2 + lpu6)
         a2pVU =  a2pUV
         a2pUW = -r1xr1zo4dq*(lmu2 - lmu6)
         a2pWU = -a2pUW
         a2pVW = -r1yr1zo4dq*(lmu6 - lmu2)
         a2pWV = -a2pVW

      else if (lbl.eq.'afm') then

         massf = wm4o8*roloc(1,1,3)

         a3mUU =  r1xqo2dq*lp2u5 + r1ypr1z*u5 - r2xo2d*lp2u5
     ;         - r2ymr2z*u5 + massf
         a3mVV =  r1yqo2dq*lp2u5 + r1xpr1z*u5 - r2yo2d*lp2u5
     ;         - r2xmr2z*u5 + massf
         a3mWW =  r1zqo2dq*lp2u5 + r1xpr1y*u5 + r2zo2d*lp2u5
     ;         - r2xmr2y*u5 + massf

         a3mUV =  r1xr1yo4dq*lpu5
         a3mVU =  a3mUV
         a3mUW = -r1xr1zo4dq*lpu5
         a3mWU =  a3mUW
         a3mVW = -r1yr1zo4dq*lpu5
         a3mWV =  a3mVW

      else if (lbl.eq.'afn') then

         massf = wm3o12*roloc(1,2,3)

         a3nUU = -r1yqo2dq*(u5 + u6) + massf
         a3nVV =  a3nUU
         a3nWW =  a3nUU

         a3nUV = -r1xr1yo4dq*(lmu5 - lmu6)
         a3nVU = -a3nUV
         a3nUW = -r1xr1zo4dq*(lpu5 + lpu6)
         a3nWU =  a3nUW
         a3nVW = -r1yr1zo4dq*(lmu5 - lmu6)
         a3nWV = -a3nVW

      else if (lbl.eq.'afp') then

         massf = wm4o8*roloc(1,3,3)

         a3pUU =  r1xqo2dq*lp2u6 + r1ypr1z*u6 - r2xo2d*lp2u6
     ;         + r2ypr2z*u6 + massf
         a3pVV =  r1yqo2dq*lp2u6 + r1xpr1z*u6 + r2yo2d*lp2u6
     ;         - r2xmr2z*u6 + massf
         a3pWW =  r1zqo2dq*lp2u6 + r1xpr1y*u6 + r2zo2d*lp2u6
     ;         - r2xpr2y*u6 + massf

         a3pUV = -r1xr1yo4dq*lpu6
         a3pVU =  a3pUV
         a3pUW = -r1xr1zo4dq*lpu6
         a3pWU =  a3pUW
         a3pVW =  r1yr1zo4dq*lpu6
         a3pWV =  a3pVW

      else if (lbl.eq.'ddm') then

         massf = wm3o12*roloc(2,1,1)

         a4mUU = -r1xqo2dq*(lp2u1 + lp2u3) + massf
         a4mVV =  a4mUU
         a4mWW = -r1xqo2dq*(u1 + u3)       + massf

         a4mUV =  r1xr1yo4dq*(lmu1 - lmu3)
         a4mVU = -a4mUV
         a4mUW =  r1xr1zo4dq*(lmu1 - lmu3)
         a4mWU = -a4mUW
         a4mVW =  r1yr1zo4dq*(lpu1 + lpu3)
         a4mWV =  a4mVW

      else if (lbl.eq.'ddn') then

         massf = wm2o6*roloc(2,2,1)

         a4nUU =  r1zqo2dq*(u1 + u4 + u2 + u3) + massf
         a4nVV =  a4nUU
         a4nWW =  r1zqo2dq*(lp2u1 + lp2u4 + lp2u2 + lp2u3) + massf

         a4nUV = -r1xr1yo4dq*(lpu1 + lpu4 - lpu2 - lpu3)
         a4nVU =  a4nUV
         a4nUW =  r1xr1zo4dq*(lmu1 - lmu4 + lmu2 - lmu3)
         a4nWU = -a4nUW
         a4nVW =  r1yr1zo4dq*(lmu1 - lmu4 + lmu3 - lmu2)
         a4nWV = -a4nVW

      else if (lbl.eq.'ddp') then

         massf = wm3o12*roloc(2,3,1)

         a4pUU = -r1xqo2dq*(lp2u2 + lp2u4) + massf
         a4pVV =  a4pUU
         a4pWW = -r1xqo2dq*(u2 + u4)       + massf

         a4pUV =  r1xr1yo4dq*(lmu4 - lmu2)
         a4pVU = -a4pUV
         a4pUW =  r1xr1zo4dq*(lmu2 - lmu4)
         a4pWU = -a4pUW
         a4pVW = -r1yr1zo4dq*(lpu2 + lpu4)
         a4pWV =  a4pVW


       else if (lbl.eq.'bem') then

         massf = wm2o6*roloc(2,1,2)

         a5mUU =  r1yqo2dq*(u1 + u7 + u5 + u3) + massf
         a5mVV =  a5mUU
         a5mWW =  a5mUU

         a5mUV =  r1xr1yo4dq*(lmu1 - lmu7 + lmu5 - lmu3)
         a5mVU = -a5mUV
         a5mUW = -r1xr1zo4dq*(lpu1 + lpu7 - lpu5 - lpu3)
         a5mWU =  a5mUW
         a5mVW = -r1yr1zo4dq*(lmu1 - lmu7 - lmu5 + lmu3)
         a5mWV = -a5mVW

      else if (lbl.eq.'ben') then

         massf = wm1*omegasq*roloc(2,2,2)

         lsum = lp2u1 + lp2u2 + lp2u3 + lp2u4 + lp2u5 + lp2u6 + lp2u7
     ;        + lp2u8
         usum = u1 + u2 + u3 + u4 + u5 + u6 + u7 + u8

         a5nUU = -r1xqo2dq*lsum - r1ypr1z*usum
     +         +  r2xo2d*(lp2u1 - lp2u8 + lp2u2 - lp2u7 + lp2u5 - lp2u4
     ;                  + lp2u6 - lp2u3) 
     +         + r2ypr2z*(u1 - u8 + u3 - u6)
     ;         + r2ymr2z*(u7 - u2 + u5 - u4)
     +         + massf
         a5nVV = -r1yqo2dq*lsum - r1xpr1z*usum
     +         +  r2yo2d*(lp2u1 - lp2u2 + lp2u5 - lp2u6 + lp2u7 - lp2u8
     ;                  + lp2u3 - lp2u4) 
     +         + r2xpr2z*(u1 - u8 + u2 - u7)
     ;         + r2xmr2z*(u6 - u3 + u5 - u4)
     +         + massf
         a5nWW = -r1zqo2dq*lsum - r1xpr1y*usum
     +         +  r2zo2d*(lp2u1 - lp2u8 + lp2u2 - lp2u7 + lp2u4 - lp2u5
     ;                  + lp2u3 - lp2u6) 
     +         + r2xpr2y*(u1 - u8 + u5 - u4)
     ;         + r2xmr2y*(u2 - u7 + u6 - u3)
     +         + massf

         a5nUV = r1xr1yo4dq*(lpu2 - lpu1 + lpu7 - lpu8 + lpu3 - lpu4
     ;                     + lpu6 - lpu5)
         a5nVU = a5nUV
         a5nUW = r1xr1zo4dq*(lpu4 - lpu1 + lpu5 - lpu8 + lpu3 - lpu2
     ;                     + lpu6 - lpu7)
         a5nWU = a5nUW
         a5nVW = r1yr1zo4dq*(lpu2 - lpu1 + lpu7 - lpu8 + lpu4 - lpu3
     ;                     + lpu5 - lpu6)
         a5nWV = a5nVW

      else if (lbl.eq.'bep') then

          massf = wm2o6*roloc(2,3,2)

          a5pUU =  r1yqo2dq*(u8 + u2 + u4 + u6) + massf
          a5pVV =  a5pUU
          a5pWW =  a5pUU

          a5pUV =  r1xr1yo4dq*(lmu8 - lmu2 + lmu4 - lmu6)
          a5pVU = -a5pUV
          a5pUW = -r1xr1zo4dq*(lpu8 + lpu2 - lpu4 - lpu6)
          a5pWU =  a5pUW
          a5pVW = -r1yr1zo4dq*(lmu8 - lmu2 - lmu4 + lmu6)
          a5pWV = -a5pVW

c       write(*,*) ' in pml3d: a5pUU = ', a5pUU
c       write(*,*) ' in pml3d: ', r1yqo2dq, u2, u4, u6, u8, massf

      else if (lbl.eq.'ffm') then

         massf = wm3o12*roloc(2,1,3)

         a6mUU = -r1xqo2dq*(lp2u7 + lp2u5) + massf
         a6mVV =  a6mUU
         a6mWW = -r1xqo2dq*(u7 + u5)       + massf

         a6mUV =  r1xr1yo4dq*(lmu5 - lmu7)
         a6mVU = -a6mUV
         a6mUW =  r1xr1zo4dq*(lmu7 - lmu5)
         a6mWU = -a6mUW
         a6mVW = -r1yr1zo4dq*(lpu7 + lpu5)
         a6mWV =  a6mVW

      else if (lbl.eq.'ffn') then

         massf = wm2o6*roloc(2,2,3)

         a6nUU =  r1zqo2dq*(u8 + u5 + u7 + u6)             + massf
         a6nVV =  a6nUU
         a6nWW =  r1zqo2dq*(lp2u8 + lp2u5 + lp2u7 + lp2u6) + massf

         a6nUV = -r1xr1yo4dq*(lpu8 + lpu5 - lpu7 - lpu6)
         a6nVU =  a6nUV
         a6nUW =  r1xr1zo4dq*(lmu8 - lmu5 + lmu7 - lmu6)
         a6nWU = -a6nUW
         a6nVW =  r1yr1zo4dq*(lmu8 - lmu5 + lmu6 - lmu7)
         a6nWV = -a6nVW

      else if (lbl.eq.'ffp') then

         massf = wm3o12*roloc(2,3,3)

         a6pUU = -r1xqo2dq*(lp2u8 + lp2u6) + massf
         a6pVV =  a6pUU
         a6pWW = -r1xqo2dq*(u8 + u6)       + massf

         a6pUV =  r1xr1yo4dq*(lmu8 - lmu6)
         a6pVU = -a6pUV
         a6pUW =  r1xr1zo4dq*(lmu8 - lmu6)
         a6pWU = -a6pUW
         a6pVW =  r1yr1zo4dq*(lpu8 + lpu6)
         a6pWV =  a6pVW

      else if (lbl.eq.'cdm') then

         massf = wm4o8*roloc(3,1,1)

         a7mUU =  r1xqo2dq*lp2u3 + r1ypr1z*u3 + r2xo2d*lp2u3
     ;         - r2ypr2z*u3 + massf
         a7mVV =  r1yqo2dq*lp2u3 + r1xpr1z*u3 - r2yo2d*lp2u3
     ;         + r2xmr2z*u3 + massf
         a7mWW =  r1zqo2dq*lp2u3 + r1xpr1y*u3 - r2zo2d*lp2u3
     ;         + r2xpr2y*u3 + massf

         a7mUV = -r1xr1yo4dq*lpu3
         a7mVU =  a7mUV
         a7mUW = -r1xr1zo4dq*lpu3
         a7mWU =  a7mUW
         a7mVW =  r1yr1zo4dq*lpu3
         a7mWV =  a7mVW

      else if (lbl.eq.'cdn') then

         massf = wm3o12*roloc(3,2,1)

         a7nUU = -r1yqo2dq*(u4 + u3) + massf
         a7nVV =  a7nUU
         a7nWW =  a7nUU

         a7nUV = -r1xr1yo4dq*(lmu4 - lmu3)
         a7nVU = -a7nUV
         a7nUW = -r1xr1zo4dq*(lpu4 + lpu3)
         a7nWU =  a7nUW
         a7nVW = -r1yr1zo4dq*(lmu4 - lmu3)
         a7nWV = -a7nVW

      else if (lbl.eq.'cdp') then

         massf = wm4o8*roloc(3,3,1)

         a7pUU =  r1xqo2dq*lp2u4 + r1ypr1z*u4 + r2xo2d*lp2u4
     ;         + r2ymr2z*u4 + massf
         a7pVV =  r1yqo2dq*lp2u4 + r1xpr1z*u4 + r2yo2d*lp2u4
     ;         + r2xmr2z*u4 + massf
         a7pWW =  r1zqo2dq*lp2u4 + r1xpr1y*u4 - r2zo2d*lp2u4
     ;         + r2xmr2y*u4 + massf

         a7pUV =  r1xr1yo4dq*lpu4
         a7pVU =  a7pUV
         a7pUW = -r1xr1zo4dq*lpu4
         a7pWU =  a7pUW
         a7pVW = -r1yr1zo4dq*lpu4
         a7pWV =  a7pVW

      else if (lbl.eq.'ccm') then

         massf = wm3o12*roloc(3,1,2)

         a8mUU = -r1zqo2dq*(u7 + u3)       + massf
         a8mVV =  a8mUU
         a8mWW = -r1zqo2dq*(lp2u7 + lp2u3) + massf

         a8mUV = -r1xr1yo4dq*(lpu7 + lpu3)
         a8mVU =  a8mUV
         a8mUW = -r1xr1zo4dq*(lmu7 - lmu3)
         a8mWU = -a8mUW
         a8mVW = -r1yr1zo4dq*(lmu3 - lmu7)
         a8mWV = -a8mVW

      else if (lbl.eq.'ccn') then

         massf = wm2o6*roloc(3,2,2)

         a8nUU =  r1xqo2dq*(lp2u8 + lp2u3 + lp2u7 + lp2u4) + massf
         a8nVV =  a8nUU
         a8nWW =  r1xqo2dq*(u8 + u3 + u7 + u4)             + massf

         a8nUV =  r1xr1yo4dq*(lmu7 - lmu4 + lmu3 - lmu8)
         a8nVU = -a8nUV
         a8nUW = -r1xr1zo4dq*(lmu7 - lmu4 + lmu8 - lmu3)
         a8nWU = -a8nUW
         a8nVW =  r1yr1zo4dq*(lpu7 + lpu4 - lpu8 - lpu3)
         a8nWV =  a8nVW

      else if (lbl.eq.'ccp') then

         massf = wm3o12*roloc(3,3,2)

         a8pUU = -r1zqo2dq*(u8 + u4)       + massf
         a8pVV =  a8pUU
         a8pWW = -r1zqo2dq*(lp2u8 + lp2u4) + massf

         a8pUV =  r1xr1yo4dq*(lpu8 + lpu4)
         a8pVU =  a8pUV
         a8pUW = -r1xr1zo4dq*(lmu8 - lmu4)
         a8pWU = -a8pUW
         a8pVW = -r1yr1zo4dq*(lmu8 - lmu4)
         a8pWV = -a8pVW

      else if (lbl.eq.'cfm') then

         massf = wm4o8*roloc(3,1,3)

         a9mUU =  r1xqo2dq*lp2u7 + r1ypr1z*u7 + r2xo2d*lp2u7
     ;         - r2ymr2z*u7 + massf
         a9mVV =  r1yqo2dq*lp2u7 + r1xpr1z*u7 - r2yo2d*lp2u7
     ;         + r2xpr2z*u7 + massf
         a9mWW =  r1zqo2dq*lp2u7 + r1xpr1y*u7 + r2zo2d*lp2u7
     ;         + r2xmr2y*u7 + massf

         a9mUV = -r1xr1yo4dq*lpu7
         a9mVU =  a9mUV
         a9mUW =  r1xr1zo4dq*lpu7
         a9mWU =  a9mUW
         a9mVW = -r1yr1zo4dq*lpu7
         a9mWV =  a9mVW

      else if (lbl.eq.'cfn') then

         massf = wm3o12*roloc(3,2,3)

         a9nUU = -r1yqo2dq*(u8 + u7) + massf
         a9nVV =  a9nUU
         a9nWW =  a9nUU

         a9nUV =  r1xr1yo4dq*(lmu7 - lmu8)
         a9nVU = -a9nUV
         a9nUW =  r1xr1zo4dq*(lpu8 + lpu7)
         a9nWU =  a9nUW
         a9nVW =  r1yr1zo4dq*(lmu8 - lmu7)
         a9nWV = -a9nVW

      else if (lbl.eq.'cfp') then

         massf = wm4o8*roloc(3,3,3)

         a9pUU = r1xqo2dq*lp2u8 + r1ypr1z*u8 + r2xo2d*lp2u8
     ;         + r2ypr2z*u8 + massf
         a9pVV = r1yqo2dq*lp2u8 + r1xpr1z*u8 + r2yo2d*lp2u8
     ;         + r2xpr2z*u8 + massf
         a9pWW = r1zqo2dq*lp2u8 + r1xpr1y*u8 + r2zo2d*lp2u8
     ;         + r2xpr2y*u8 + massf

         a9pUV = r1xr1yo4dq*lpu8
         a9pVU = a9pUV
         a9pUW = r1xr1zo4dq*lpu8
         a9pWU = a9pUW
         a9pVW = r1yr1zo4dq*lpu8
         a9pWV = a9pVW

      endif

      return
      end


