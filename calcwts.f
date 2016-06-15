      SUBROUTINE calcwts(vp, vs, w1, w2, w3, wm1, wm2, wm3, wm4)
c
c     calculates optimal weighting coefficients
c
c 
c  
c     input    type     meaning 
c     -----    ------   ------- 
c     vp         real*4   p wavespeed
c     vs         real*4   s wavespeed (same units as p; we only use vp/vs here)
c 
c     output    type    meaning 
c     ------    ------  ------- 
c     w1-w3     real*4  frame weights; we take w3 = 1 and the others to be zero 
c     wm1-wm4   real*4  mass lumping weights 
c 
c-----variable declarations 
      IMPLICIT NONE
      REAL, INTENT(IN) :: vp, vs
      REAL, INTENT(OUT) :: w1, w2, w3, wm1, wm2, wm3, wm4

c     include 'mpif.h'
c     integer*4 mpierr
c
c---------------------------------------------------------------------c
c 

      DOUBLE PRECISION wm1p(4), wm2p(4), wm3p(4),
     ;                 x, w1d, w2d, w3d, wmd1, wmd2, wmd3, wmd4

      DATA wm1p / 0.2909366D+00, -0.1792402D-01,
     ;            0.7329184D-02, -0.1055034D-02/
      DATA wm2p / 0.6367377D+00, -0.7982022D-02,
     ;            0.3385164D-02, -0.5065399D-03/
      DATA wm3p / 0.7244875D-01,  0.2569959D-01,
     ;           -0.1060011D-01,  0.1540534D-02/

      x = vp/vs

      wmd2 = 0.d0
      wmd3 = 0.d0
      wmd4 = 0.d0

      wmd1 = wm1p(1) + x*(wm1p(2) + x*(wm1p(3) + x*wm1p(4)))
      if (wmd1.gt.1.d0) then
         wmd1 = 1.d0
      else
         wmd2 = wm2p(1) + x*(wm2p(2) + x*(wm2p(3) + x*wm2p(4)))
      endif
      if ((wmd1 + wmd2).gt.1.d0) then
         wmd2 = 0.d0
      else
         wmd3 = wm3p(1) + x*(wm3p(2) + x*(wm3p(3) + x*wm3p(4)))
      endif
      if ((wmd1 + wmd2 + wmd3).gt.1.d0) then
         wmd3 = 1.d0 - (wmd1 + wmd2)
      endif

      wmd4 = 1.d0 - (wmd1 + wmd2 + wmd3)

      w1d = 1.d0
      w2d = 0.d0
      w3d = 0.d0

      w1 = SNGL(w1d)
      w2 = SNGL(w2d)
      w3 = SNGL(w3d)
      wm1 = SNGL(wmd1)
      wm2 = SNGL(wmd2)
      wm3 = SNGL(wmd3)
      wm4 = SNGL(wmd4)

      return 
      end 

