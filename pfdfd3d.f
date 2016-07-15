      SUBROUTINE pfdfd3d (lbl, ierr)

c----- This routine generates the entries on the 3D elastic FDFD system of equations.
c       S. Roecker 2/2012
c
c
c       Input:
c       iz, iy, ix  are the model indices for the point in question
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
      USE COEFFS_MODULE, ONLY : daloc, muloc, roloc,
     ;                          lmu1, lmu2, lmu3, lmu4,
     ;                          lmu5, lmu6, lmu7, lmu8,
     ;                          lpu1, lpu2, lpu3, lpu4,
     ;                          lpu5, lpu6, lpu7, lpu8,
     ;                          lp2u1, lp2u2, lp2u3, lp2u4, lp2u5,
     ;                          lp2u6, lp2u7, lp2u8,
     ;                          u1, u2, u3, u4, u5, u6, u7, u8,
     ;                          wm1, wm2, wm3, wm4
      USE FD_ENUM_MODULE, ONLY : bem, ben, bep,
     ;                           adm, adn, adp,
     ;                           aam, aan, aap,
     ;                           afm, afn, afp,
     ;                           ddm, ddn, ddp,
     ;                           ffm, ffn, ffp,
     ;                           cdm, cdn, cdp,
     ;                           ccm, ccn, ccp,
     ;                           cfm, cfn, cfp
      USE MODEL_MODULE, ONLY : omega, dx
      USE MAT_MODULE
      IMPLICIT NONE


c parameters:
c       include 'dimension.inc'
c
c arguments:
      INTEGER, INTENT(IN) :: lbl
      INTEGER, INTENT(OUT) :: ierr

c local variables:
      REAL    dx2, dxx, dxx2, dxx4
      REAL    vpreal

      REAL    wm2o6, wm3o12, wm4o8

c----Note:  aro and vp are declared complex because they are
c    involved with complex operations, but the density (aro) is
c    always real and vp SHOULD be real for the terms the use it.
c    If we have finite Q, however, vp will be complex, and we
c    force it to be real.
      COMPLEX aro, vp   

      REAL d1dxsq, d2dxsq, d4dxsq

      COMPLEX omegasq
      COMPLEX massf

c common variables:

c       include 'common.inc'
c       include 'init.inc'
c       include 'coeffs.inc'

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
c

c---turn off mass lumping for a test
c       wm1 = 1.0
c       wm2 = 0.0
c       wm3 = 0.0
c       wm4 = 0.0

      ierr = 0 
          
c executable code:

      omegasq = omega * omega

!TODO: check this
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

c----Note that the B1-4 frame is the sum of 4 coordinate systems, so we
c    should divide these terms by 4.
      d1dxsq = 1.0/(dx*dx)
      d1dxsq = d1dxsq/4.0
      d2dxsq = d1dxsq/2.0
      d4dxsq = d2dxsq/2.0

      if (lbl.eq.adm) then

         massf = wm4o8*roloc(1,1,1)

         a1mUU =  d2dxsq*lp2u1 + d1dxsq*u1 + massf
         a1mVV =  a1mUU
         a1mWW =  a1mUU

         a1mUV =  d4dxsq*lpu1
         a1mVU =  a1mUV
         a1mUW =  a1mUV
         a1mWU =  a1mUV
         a1mVW =  a1mUV
         a1mWV =  a1mUV

      else if (lbl.eq.adn) then

         massf = wm3o12*roloc(1,2,1)

         a1nUU = -d2dxsq*(u1 + u2) + massf
         a1nVV =  a1nUU
         a1nWW =  a1nUU

         a1nUV =  d4dxsq*(lmu2 - lmu1)
         a1nVU =  -a1nUV

         a1nUW =  d4dxsq*(lpu1 + lpu2)
         a1nWU =  a1nUW

         a1nVW =  d4dxsq*(lmu1 - lmu2)
         a1nWV = -a1nVW


      else if (lbl.eq.adp) then

         massf = wm4o8*roloc(1,3,1)

         a1pUU =  d2dxsq*lp2u2 + d1dxsq*u2 + massf
         a1pVV =  a1pUU
         a1pWW =  a1pUU

         a1pUV = -d4dxsq*lpu2
         a1pVU =  a1pUV
         a1pUW = -a1pUV
         a1pWU = -a1pUV
         a1pVW =  a1pUV
         a1pWV =  a1pUV

      else if (lbl.eq.aam) then

         massf = wm3o12*roloc(1,1,2)

         a2mUU = -d2dxsq*(u1 + u5)       + massf
         a2mVV =  a2mUU
         a2mWW = -d2dxsq*(lp2u1 + lp2u5) + massf

         a2mUV =  d4dxsq*(lpu1 + lpu5)
         a2mVU =  a2mUV
         a2mUW = -d4dxsq*(lmu1 - lmu5)
         a2mWU = -a2mUW
         a2mVW =  a2mUW
         a2mWV = -a2mUW 

      else if (lbl.eq.aan) then

         massf = wm2o6*roloc(1,2,2)

         a2nUU =  d2dxsq*(lp2u1 + lp2u6 + lp2u2 + lp2u5) + massf
         a2nVV =  a2nUU
         a2nWW =  d2dxsq*(u1 + u6 + u2 + u5)             + massf

         a2nUV =  d4dxsq*(lmu2 - lmu5 + lmu6 - lmu1)
         a2nVU = -a2nUV
         a2nUW = -d4dxsq*(lmu2 - lmu5 + lmu1 - lmu6)
         a2nWU = -a2nUW
         a2nVW =  d4dxsq*(lpu2 + lpu5 - lpu1 - lpu6)
         a2nWV =  a2nVW

      else if (lbl.eq.aap) then

         massf = wm3o12*roloc(1,3,2)

         a2pUU = -d2dxsq*(u2 + u6)       + massf
         a2pVV =  a2pUU
         a2pWW = -d2dxsq*(lp2u2 + lp2u6) + massf

         a2pUV = -d4dxsq*(lpu2 + lpu6)
         a2pVU =  a2pUV
         a2pUW = -d4dxsq*(lmu2 - lmu6)
         a2pWU = -a2pUW
         a2pVW = -a2pUW
         a2pWV =  a2pUW

      else if (lbl.eq.afm) then

         massf = wm4o8*roloc(1,1,3)

         a3mUU =  d2dxsq*lp2u5 + d1dxsq*u5 + massf
         a3mVV =  a3mUU
         a3mWW =  a3mUU

         a3mUV =  d4dxsq*lpu5
         a3mVU =  a3mUV
         a3mUW = -a3mUV
         a3mWU = -a3mUV
         a3mVW = -a3mUV
         a3mWV = -a3mUV

      else if (lbl.eq.afn) then

         massf = wm3o12*roloc(1,2,3)

         a3nUU = -d2dxsq*(u5 + u6) + massf
         a3nVV =  a3nUU
         a3nWW =  a3nUU

         a3nUV = -d4dxsq*(lmu5 - lmu6)
         a3nVU = -a3nUV
         a3nUW = -d4dxsq*(lpu5 + lpu6)
         a3nWU =  a3nUW
         a3nVW =  a3nUV
         a3nWV = -a3nUV

      else if (lbl.eq.afp) then

         massf = wm4o8*roloc(1,3,3)

         a3pUU =  d2dxsq*lp2u6 + d1dxsq*u6 + massf
         a3pVV =  a3pUU
         a3pWW =  a3pUU

         a3pUV = -d4dxsq*lpu6
         a3pVU =  a3pUV
         a3pUW =  a3pUV
         a3pWU =  a3pUV
         a3pVW = -a3pUV
         a3pWV = -a3pUV

      else if (lbl.eq.ddm) then

         massf = wm3o12*roloc(2,1,1)

         a4mUU = -d2dxsq*(lp2u1 + lp2u3) + massf
         a4mVV =  a4mUU
         a4mWW = -d2dxsq*(u1 + u3)       + massf

         a4mUV =  d4dxsq*(lmu1 - lmu3)
         a4mVU = -a4mUV
         a4mUW =  a4mUV
         a4mWU = -a4mUV
         a4mVW =  d4dxsq*(lpu1 + lpu3)
         a4mWV =  a4mVW

      else if (lbl.eq.ddn) then

         massf = wm2o6*roloc(2,2,1)

         a4nUU =  d2dxsq*(u1 + u4 + u2 + u3)             + massf
         a4nVV =  a4nUU
         a4nWW =  d2dxsq*(lp2u1 + lp2u4 + lp2u2 + lp2u3) + massf

         a4nUV = -d4dxsq*(lpu1 + lpu4 - lpu2 - lpu3)
         a4nVU =  a4nUV
         a4nUW =  d4dxsq*(lmu1 - lmu4 + lmu2 - lmu3)
         a4nWU = -a4nUW
         a4nVW =  d4dxsq*(lmu1 - lmu4 + lmu3 - lmu2)
         a4nWV = -a4nVW

      else if (lbl.eq.ddp) then

         massf = wm3o12*roloc(2,3,1)

         a4pUU = -d2dxsq*(lp2u2 + lp2u4) + massf
         a4pVV =  a4pUU
         a4pWW = -d2dxsq*(u2 + u4)       + massf

         a4pUV =  d4dxsq*(lmu4 - lmu2)
         a4pVU = -a4pUV
         a4pUW =  d4dxsq*(lmu2 - lmu4)
         a4pWU = -a4pUW
         a4pVW = -d4dxsq*(lpu2 + lpu4)
         a4pWV =  a4pVW

      else if (lbl.eq.bem) then

         massf = wm2o6*roloc(2,1,2)

         a5mUU =  d2dxsq*(u1 + u7 + u5 + u3) + massf
         a5mVV =  a5mUU
         a5mWW =  a5mUU

         a5mUV =  d4dxsq*(lmu1 - lmu7 + lmu5 - lmu3)
         a5mVU = -a5mUV
         a5mUW = -d4dxsq*(lpu1 + lpu7 - lpu5 - lpu3)
         a5mWU =  a5mUW
         a5mVW = -d4dxsq*(lmu1 - lmu7 - lmu5 + lmu3)
         a5mWV = -a5mVW

      else if (lbl.eq.ben) then

         massf = wm1*omegasq*roloc(2,2,2)

         a5nUU = -d2dxsq*(lp2u1 + lp2u2 + lp2u3 + lp2u4 + lp2u5 + lp2u6 
     ;                  + lp2u7 + lp2u8) 
     +           - d1dxsq*(u1 + u2 + u3 + u4 + u5 + u6 + u7 + u8)
     +           + massf
         a5nVV =  a5nUU
         a5nWW =  a5nUU

         a5nUV = d4dxsq*(lpu2 - lpu1 + lpu7 - lpu8 + lpu3 - lpu4 + lpu6
     ;                 - lpu5)
         a5nVU = a5nUV
         a5nUW = d4dxsq*(lpu4 - lpu1 + lpu5 - lpu8 + lpu3 - lpu2 + lpu6
     ;                 - lpu7)
         a5nWU = a5nUW
         a5nVW = d4dxsq*(lpu2 - lpu1 + lpu7 - lpu8 + lpu4 - lpu3 + lpu5
     ;                 - lpu6)
         a5nWV = a5nVW

      else if (lbl.eq.bep) then

         massf = wm2o6*roloc(2,3,2)

         a5pUU =  d2dxsq*(u8 + u2 + u4 + u6) + massf
         a5pVV =  a5pUU
         a5pWW =  a5pUU

         a5pUV =  d4dxsq*(lmu8 - lmu2 + lmu4 - lmu6)
         a5pVU = -a5pUV
         a5pUW = -d4dxsq*(lpu8 + lpu2 - lpu4 - lpu6)
         a5pWU =  a5pUW
         a5pVW = -d4dxsq*(lmu8 - lmu2 - lmu4 + lmu6)
         a5pWV = -a5pVW

      else if (lbl.eq.ffm) then

         massf = wm3o12*roloc(2,1,3)

         a6mUU = -d2dxsq*(lp2u7 + lp2u5) + massf
         a6mVV =  a6mUU
         a6mWW = -d2dxsq*(u7 + u5)       + massf

         a6mUV =  d4dxsq*(lmu5 - lmu7)
         a6mVU = -a6mUV
         a6mUW = -a6mUV
         a6mWU =  a6mUV
         a6mVW = -d4dxsq*(lpu7 + lpu5)
         a6mWV =  a6mVW

      else if (lbl.eq.ffn) then

         massf = wm2o6*roloc(2,2,3)

         a6nUU =  d2dxsq*(u8 + u5 + u7 + u6)             + massf
         a6nVV =  a6nUU
         a6nWW =  d2dxsq*(lp2u8 + lp2u5 + lp2u7 + lp2u6) + massf

         a6nUV = -d4dxsq*(lpu8 + lpu5 - lpu7 - lpu6)
         a6nVU =  a6nUV
         a6nUW =  d4dxsq*(lmu8 - lmu5 + lmu7 - lmu6)
         a6nWU = -a6nUW
         a6nVW =  d4dxsq*(lmu8 - lmu5 + lmu6 - lmu7)
         a6nWV = -a6nVW

      else if (lbl.eq.ffp) then

         massf = wm3o12*roloc(2,3,3)

         a6pUU = -d2dxsq*(lp2u8 + lp2u6) + massf
         a6pVV =  a6pUU
         a6pWW = -d2dxsq*(u8 + u6)       + massf

         a6pUV =  d4dxsq*(lmu8 - lmu6)
         a6pVU = -a6pUV
         a6pUW =  a6pUV
         a6pWU = -a6pUV
         a6pVW =  d4dxsq*(lpu8 + lpu6)
         a6pWV =  a6pVW

      else if (lbl.eq.cdm) then

         massf = wm4o8*roloc(3,1,1)

         a7mUU =  d2dxsq*lp2u3 + d1dxsq*u3 + massf
         a7mVV =  a7mUU
         a7mWW =  a7mUU

         a7mUV = -d4dxsq*lpu3
         a7mVU =  a7mUV
         a7mUW =  a7mUV
         a7mWU =  a7mUV
         a7mVW = -a7mUV
         a7mWV = -a7mUV

      else if (lbl.eq.cdn) then

         massf = wm3o12*roloc(3,2,1)

         a7nUU = -d2dxsq*(u4 + u3) + massf
         a7nVV =  a7nUU
         a7nWW =  a7nUU

         a7nUV = -d4dxsq*(lmu4 - lmu3)
         a7nVU = -a7nUV
         a7nUW = -d4dxsq*(lpu4 + lpu3)
         a7nWU =  a7nUW
         a7nVW =  a7nUV
         a7nWV = -a7nUV

      else if (lbl.eq.cdp) then

         massf = wm4o8*roloc(3,3,1)

         a7pUU =  d2dxsq*lp2u4 + d1dxsq*u4 + massf
         a7pVV =  a7pUU
         a7pWW =  a7pUU

         a7pUV =  d4dxsq*lpu4
         a7pVU =  a7pUV
         a7pUW = -a7pUV
         a7pWU = -a7pUV
         a7pVW = -a7pUV
         a7pWV = -a7pUV

      else if (lbl.eq.ccm) then

         massf = wm3o12*roloc(3,1,2)

         a8mUU = -d2dxsq*(u7 + u3)       + massf
         a8mVV =  a8mUU
         a8mWW = -d2dxsq*(lp2u7 + lp2u3) + massf

         a8mUV = -d4dxsq*(lpu7 + lpu3)
         a8mVU =  a8mUV
         a8mUW = -d4dxsq*(lmu7 - lmu3)
         a8mWU = -a8mUW
         a8mVW = -a8mUW
         a8mWV =  a8mUW

      else if (lbl.eq.ccn) then

         massf = wm2o6*roloc(3,2,2)

         a8nUU =  d2dxsq*(lp2u8 + lp2u3 + lp2u7 + lp2u4) + massf
         a8nVV =  a8nUU
         a8nWW =  d2dxsq*(u8 + u3 + u7 + u4)             + massf

         a8nUV =  d4dxsq*(lmu7 - lmu4 + lmu3 - lmu8)
         a8nVU = -a8nUV
         a8nUW = -d4dxsq*(lmu7 - lmu4 + lmu8 - lmu3)
         a8nWU = -a8nUW
         a8nVW =  d4dxsq*(lpu7 + lpu4 - lpu8 - lpu3)
         a8nWV =  a8nVW

      else if (lbl.eq.ccp) then

         massf = wm3o12*roloc(3,3,2)

         a8pUU = -d2dxsq*(u8 + u4)       + massf
         a8pVV =  a8pUU
         a8pWW = -d2dxsq*(lp2u8 + lp2u4) + massf

         a8pUV =  d4dxsq*(lpu8 + lpu4)
         a8pVU =  a8pUV
         a8pUW = -d4dxsq*(lmu8 - lmu4)
         a8pWU = -a8pUW
         a8pVW =  a8pUW
         a8pWV = -a8pUW

      else if (lbl.eq.cfm) then

         massf = wm4o8*roloc(3,1,3)

         a9mUU =  d2dxsq*lp2u7 + d1dxsq*u7 + massf
         a9mVV =  a9mUU
         a9mWW =  a9mUU

         a9mUV = -d4dxsq*lpu7
         a9mVU =  a9mUV
         a9mUW = -a9mUV
         a9mWU = -a9mUV
         a9mVW =  a9mUV
         a9mWV =  a9mUV

      else if (lbl.eq.cfn) then

         massf = wm3o12*roloc(3,2,3)

         a9nUU = -d2dxsq*(u8 + u7) + massf
         a9nVV =  a9nUU
         a9nWW =  a9nUU

         a9nUV =  d4dxsq*(lmu7 - lmu8)
         a9nVU = -a9nUV
         a9nUW =  d4dxsq*(lpu8 + lpu7)
         a9nWU =  a9nUW
         a9nVW = -a9nUV
         a9nWV =  a9nUV

      else if (lbl.eq.cfp) then

         massf = wm4o8*roloc(3,3,3)

         a9pUU =  d2dxsq*lp2u8 + d1dxsq*u8 + massf
         a9pVV =  a9pUU
         a9pWW =  a9pUU

         a9pUV =  d4dxsq*lpu8
         a9pVU =  a9pUV
         a9pUW =  a9pUV
         a9pWU =  a9pUV
         a9pVW =  a9pUV
         a9pWV =  a9pUV

      endif

      return
      end
