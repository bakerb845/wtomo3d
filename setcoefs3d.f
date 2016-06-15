      SUBROUTINE setcoefs3d(iz, ix, iy, ibg)

c----This routine assigns all the Lame parameters and densities 
c       related to the node located at point (iz, ix, iy) to 3x3x3
c       matrices for use in pfdfd3d and pml3d.
c
c       ix, iy, iz are the indices of the central grid point
c       ibg = 0 means use the current model for assignments
c       ibg = 1 means use the background model for assignments
c       

      USE INIT_MODULE, ONLY : freesurf, nx, ny, nz
      USE MODEL_MODULE, ONLY : da, dap, mu, mup, rho, rhop
      USE COEFFS_MODULE, ONLY : daloc, muloc, roloc,
     ;                          da1, da2, da3, da4,
     ;                          da5, da6, da7, da8,
     ;                          lmu1, lmu2, lmu3, lmu4,
     ;                          lmu5, lmu6, lmu7, lmu8,
     ;                          lpu1, lpu2, lpu3, lpu4,
     ;                          lpu5, lpu6, lpu7, lpu8,
     ;                          lp2u1, lp2u2, lp2u3, lp2u4,
     ;                          lp2u5, lp2u6, lp2u7, lp2u8,
     ;                          u1, u2, u3, u4, u5, u6, u7, u8,
     ;                          w1, w2, w3, wm1, wm2, wm3, wm4
      USE INTERFACE_MODULE, ONLY : calcwts
      IMPLICIT NONE

c       include 'dimension.inc'            
c       include 'common.inc'            
c       include 'init.inc'            
c       include 'coeffs.inc'            

c-- arguments
      INTEGER, INTENT(IN) :: iz, ix, iy, ibg

c-- local variables
c   Note that in this case we requre the real part of the velocities in order
c       to estimate the a and b coefficients
      REAL vp, vs
c---system weights
      REAL w1c, w2c, w3c
c---mass lumping weights
      REAL wm1c, wm2c, wm3c, wm4c

      INTEGER ixg, iyg, izg
      INTEGER i, j, k

      COMPLEX mukm11, mukm12, mukm13, mukm14
      COMPLEX mukm21, mukm22, mukm23, mukm24
      COMPLEX mukm31, mukm32, mukm33, mukm34
      COMPLEX dakm11, dakm12, dakm13, dakm14
      COMPLEX dakm21, dakm22, dakm23, dakm24
      COMPLEX dakm31, dakm32, dakm33, dakm34

c----assigns general coefficients for model node iz, ix, iy
c     mu holds the shear moduli, da the lamda's.  
c       The naming conventions follow my notes which order
c       indexes by (ix, iz) rather than the more computationally
c       efficient  (iz, ix) that is used in this routine.  
c       SWR 2/12
c
c       Nomenclature:
c
c       u, l    first character denoting mu or lambda
c       p, m, n plus, minus, null operations
c       h       half step
c       i, j, k i (x), j (y), or k (z) index.  If not present, then no displacement.
c
c       Generally we presume an i,j,k order on indices, and these will be present
c       on each variable.  So, for example
c
c       lp2upphn        would be (lambda +2 mu) (i+1, j+1/2, k)
c

c----assign the local arrays.  Note the the ordering of the local arrays is (ix,jy,kz) while
c    that of global arrays is (kz, ix, jy).  The reason for this difference is code legacy for
c    global arryas and fidelity with my notes for local arrays.
c
      do k = 1, 3
         izg = iz + k - 2
         if (izg.lt.1)  izg = 1
         if (izg.gt.nz) izg = nz
         do i = 1, 3
            ixg = ix + i - 2
            if (ixg.lt.1)  ixg = 1
            if (ixg.gt.nx) ixg = nx
            do j = 1, 3
               iyg = iy + j - 2
               if (iyg.lt.1)  iyg = 1
               if (iyg.gt.ny) iyg = ny
               if (ibg.eq.0) then
                  daloc(i,j,k) =  da(izg, ixg, iyg)
                  muloc(i,j,k) =  mu(izg, ixg, iyg)
                  roloc(i,j,k) = rho(izg, ixg, iyg)
               else
                  daloc(i,j,k) =  dap(izg, ixg, iyg)
                  muloc(i,j,k) =  mup(izg, ixg, iyg)
                  roloc(i,j,k) = rhop(izg, ixg, iyg)
               endif
            enddo
         enddo
      enddo

c---generate values for the 8 Min cells
      mukm11 = (muloc(1,1,1) + muloc(2,1,1)
     ;        + muloc(1,2,1) + muloc(2,2,1))/4.0
      mukm12 = (muloc(1,3,1) + muloc(2,3,1)
     ;        + muloc(1,2,1) + muloc(2,2,1))/4.0
      mukm13 = (muloc(3,1,1) + muloc(2,1,1)
     ;        + muloc(3,2,1) + muloc(2,2,1))/4.0
      mukm14 = (muloc(3,3,1) + muloc(2,3,1)
     ;        + muloc(3,2,1) + muloc(2,2,1))/4.0
      mukm21 = (muloc(1,1,2) + muloc(2,1,2)
     ;        + muloc(1,2,2) + muloc(2,2,2))/4.0
      mukm22 = (muloc(1,3,2) + muloc(2,3,2)
     ;        + muloc(1,2,2) + muloc(2,2,2))/4.0
      mukm23 = (muloc(3,1,2) + muloc(2,1,2)
     ;        + muloc(3,2,2) + muloc(2,2,2))/4.0
      mukm24 = (muloc(3,3,2) + muloc(2,3,2)
     ;        + muloc(3,2,2) + muloc(2,2,2))/4.0
      mukm31 = (muloc(1,1,3) + muloc(2,1,3)
     ;        + muloc(1,2,3) + muloc(2,2,3))/4.0
      mukm32 = (muloc(1,3,3) + muloc(2,3,3)
     ;        + muloc(1,2,3) + muloc(2,2,3))/4.0
      mukm33 = (muloc(3,1,3) + muloc(2,1,3)
     ;        + muloc(3,2,3) + muloc(2,2,3))/4.0
      mukm34 = (muloc(3,3,3) + muloc(2,3,3)
     ;        + muloc(3,2,3) + muloc(2,2,3))/4.0

      u1 = (mukm11 + mukm21)/2.0
      u2 = (mukm12 + mukm22)/2.0
      u3 = (mukm13 + mukm23)/2.0
      u4 = (mukm14 + mukm24)/2.0
      u5 = (mukm31 + mukm21)/2.0
      u6 = (mukm32 + mukm22)/2.0
      u7 = (mukm33 + mukm23)/2.0
      u8 = (mukm34 + mukm24)/2.0

      dakm11 = (daloc(1,1,1) + daloc(2,1,1)
     ;        + daloc(1,2,1) + daloc(2,2,1))/4.0
      dakm12 = (daloc(1,3,1) + daloc(2,3,1)
     ;        + daloc(1,2,1) + daloc(2,2,1))/4.0
      dakm13 = (daloc(3,1,1) + daloc(2,1,1)
     ;        + daloc(3,2,1) + daloc(2,2,1))/4.0
      dakm14 = (daloc(3,3,1) + daloc(2,3,1)
     ;        + daloc(3,2,1) + daloc(2,2,1))/4.0
      dakm21 = (daloc(1,1,2) + daloc(2,1,2)
     ;        + daloc(1,2,2) + daloc(2,2,2))/4.0
      dakm22 = (daloc(1,3,2) + daloc(2,3,2)
     ;        + daloc(1,2,2) + daloc(2,2,2))/4.0
      dakm23 = (daloc(3,1,2) + daloc(2,1,2)
     ;        + daloc(3,2,2) + daloc(2,2,2))/4.0
      dakm24 = (daloc(3,3,2) + daloc(2,3,2)
     ;        + daloc(3,2,2) + daloc(2,2,2))/4.0
      dakm31 = (daloc(1,1,3) + daloc(2,1,3)
     ;        + daloc(1,2,3) + daloc(2,2,3))/4.0
      dakm32 = (daloc(1,3,3) + daloc(2,3,3)
     ;        + daloc(1,2,3) + daloc(2,2,3))/4.0
      dakm33 = (daloc(3,1,3) + daloc(2,1,3)
     ;        + daloc(3,2,3) + daloc(2,2,3))/4.0
      dakm34 = (daloc(3,3,3) + daloc(2,3,3)
     ;        + daloc(3,2,3) + daloc(2,2,3))/4.0

      da1 = (dakm11 + dakm21)/2.0
      da2 = (dakm12 + dakm22)/2.0
      da3 = (dakm13 + dakm23)/2.0
      da4 = (dakm14 + dakm24)/2.0
      da5 = (dakm31 + dakm21)/2.0
      da6 = (dakm32 + dakm22)/2.0
      da7 = (dakm33 + dakm23)/2.0
      da8 = (dakm34 + dakm24)/2.0

c---free surface:  the condition at the free surface is that values in cells 1-4 are set to zero.
      if (freesurf(1).and.(iz.eq.1)) then
         u1  = CMPLX(0., 0.)
         u2  = CMPLX(0., 0.)
         u3  = CMPLX(0., 0.)
         u4  = CMPLX(0., 0.)
         da1 = CMPLX(0., 0.)
         da2 = CMPLX(0., 0.)
         da3 = CMPLX(0., 0.)
         da4 = CMPLX(0., 0.)
c---set the upper layer densities to zero so they are not included in mass lumping
         do j = 1, 3
            do i = 1, 3
              roloc(i,j,1) = 0.0
            enddo
         enddo
      endif

c---set up some specialty values from these cells

      lmu1 = da1 - u1
      lmu2 = da2 - u2
      lmu3 = da3 - u3
      lmu4 = da4 - u4
      lmu5 = da5 - u5
      lmu6 = da6 - u6
      lmu7 = da7 - u7
      lmu8 = da8 - u8

      lpu1 = da1 + u1
      lpu2 = da2 + u2
      lpu3 = da3 + u3
      lpu4 = da4 + u4
      lpu5 = da5 + u5
      lpu6 = da6 + u6
      lpu7 = da7 + u7
      lpu8 = da8 + u8

      lp2u1 = lpu1 + u1
      lp2u2 = lpu2 + u2
      lp2u3 = lpu3 + u3
      lp2u4 = lpu4 + u4
      lp2u5 = lpu5 + u5
      lp2u6 = lpu6 + u6
      lp2u7 = lpu7 + u7
      lp2u8 = lpu8 + u8

c---finally, set the system and mass lumping weights
      vp = REAL(CSQRT((daloc(2,2,2) + 2.0*muloc(2,2,2))/roloc(2,2,2)))
      vs = REAL(CSQRT(muloc(2,2,2)/roloc(2,2,2)))
      CALL calcwts(vp, vs, w1c, w2c, w3c, wm1c, wm2c, wm3c, wm4c)
      w1 = w1c
      w2 = w2c
      w3 = w3c
      wm1 = wm1c
      wm2 = wm2c
      wm3 = wm3c
      wm4 = wm4c

      return
      end

