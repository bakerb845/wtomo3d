c----------------------------------------------------------------------
c
c       Compute the Green's functions for a 1D medium using TH (Thomson-Haskell).  
c
c       In this version we use the model contained in vp1disc and vs1disc
c       which nonmialy is the 1D model corrected for dipsersion (i.e., they
c       are the wavespeeds in the discrete model, as opposed to the continuous
c       model.
c
c       This version differs from the orignal haskgrn only in the use of
c       vp1disc and vs1disc instead of vp1d and vs1d.
c
c       This routine is the computative "engine" of the T-H method.
c       It impliments the concepts in Haskell (1953; H53) and Haskell
c       (1962; H62) to compute the Green's functions in layered media
c       do to an external plane wave source of unit amplitude
c
c       Arguments:
c       df      The frequency 
c       c       The apparent velocity
c       ro      The density array
c       h       The thickness array
c       ctype   = P for incident P wave, S for incident S wave
c
c       This version for 1 frequency = df
c
c       This version allows the addition of a fluid layer (vs = 0)
c
c       Author:  S. Roecker 4/09
c
      SUBROUTINE haskgrn(df, isrc, ierr)
      USE MODEL_MODULE, ONLY : dz, nz, zorig
      USE SOURCE_MODULE, ONLY : ain, srctyp
      USE HASKELL_MODULE, ONLY : nl, ugrn1f, wgrn1f, vp1d, vs1d,
     ;                           rh1d, h1d, vs1disc, vp1disc
      IMPLICIT NONE

c     include 'dimension.inc'
c     include 'common.inc'
c     include 'init.inc'
      REAL, INTENT(IN) :: df
      INTEGER, INTENT(IN) :: isrc
      INTEGER, INTENT(OUT) :: ierr
c     dimension r1(nzz),r2(nzz),aag(nzz),aag1(nzz),
c    *       tt(nzz),xh(nzz),xh1(nzz)
c     dimension e1(nzz),e2(nzz),s01(nzz),c01(nzz),
c    +       s02(nzz),c02(nzz)

      DOUBLE PRECISION, ALLOCATABLE :: s01(:), c01(:), c02(:), s02(:),
     ;                                 aag(:), aag1(:), e1(:), e2(:),
     ;                                 h(:), r1(:), r2(:), tt(:),
     ;                                 xh(:), xh1(:)

      CHARACTER(1) ctype

      DOUBLE PRECISION  x,y,ag,ag1,tm,wr,wi,ur,ui,pic,hc,
     &                  pd,qd,pm,qm,cp,cq
      DOUBLE PRECISION  c11,c12,c21,c22,
     &                  c31,c32,c41,c42
      DOUBLE PRECISION  b11,b12,b21,b22,
     &                  b31,b32,b41,b42
      DOUBLE PRECISION  d11,d12,d21,d22,
     &                  d31,d32,d41,d42
      DOUBLE PRECISION  a11,a12,a13,a14,a21,a22,a23,a24,
     &                  a31,a32,a33,a34,a41,a42,a43,a44

c     REAL aag, aag1 
      DOUBLE PRECISION cova, covb, dept

      DOUBLE PRECISION e11, e13, e22, e24, e33, e31, e42, e44
      DOUBLE PRECISION dloc
      DOUBLE PRECISION qiu, qiw, qru, qrw
c     REAL r1, r2
      DOUBLE PRECISION scl1, scl2
c     REAL tt, xh, xh1
      DOUBLE PRECISION z01, z02
      DOUBLE PRECISION cv

      INTEGER i, m, iz
      INTEGER ibov, ibovo, ifin

c     INTEGER itest

c       data pi /3.141592654/
!     allocate space
      ierr = 0
      ALLOCATE(s01(nl))
      ALLOCATE(c01(nl))
      ALLOCATE(c02(nl))
      ALLOCATE(s02(nl))
      ALLOCATE(aag(nl))
      ALLOCATE(aag1(nl))
      ALLOCATE(e1(nl))
      ALLOCATE(e2(nl))
      ALLOCATE(h(nl))
      ALLOCATE(r1(nl))
      ALLOCATE(r2(nl))
      ALLOCATE(tt(nl))
      ALLOCATE(xh(nl))
      ALLOCATE(xh1(nl))
 
c----incident wave type
      ierr = 0
      ctype = 'P'
      if (srctyp(isrc).eq.'ps') ctype = 'S'


c----horizontal wavespeed
c  NB:  it seems to me that if we do the dispersion correction we use the apparent
c  wavespeed from the original vp1d/vs1d and ain(isrc), so to be consistent we should
c  use vp1d and vs1d here.  Otherwise we should either redefine ain, or define cv elsewhere.
      if (ctype.eq.'P') then
c        cv = vp1disc(nl)/sin(ain(isrc))
         cv = vp1d(nl)/sin(ain(isrc))
      else
c        cv = vs1disc(nl)/sin(ain(isrc))
         cv = vs1d(nl)/sin(ain(isrc))
      endif

c----df is already computed as omega before this routine is called
c     pic = 2.*pi*df/cv
      pic = DBLE(df)/cv

c       write(*,*) ' cv, pic = ', cv, pic

c----thickness array
      do i = 1, nl-1
         h(i) = h1d(i+1) - h1d(i)
      enddo

c----Loop over the model, setting up layer specific values
      do m = 1,nl

c--- c/alpha and c/beta
         cova = cv/vp1disc(m)
c---set up a bogus value for covb in case of vs = 0
         if (vs1disc(m).ne.0.0) then
            covb = cv/vs1disc(m)
         else
            covb = 1.0
         endif 

c----r1 and r2 are r(alpha) and r(beta)
c----aag is gamma, aag1 is gamma-1
         x = cova*cova - 1.
         r1(m) = DSIGN(DSQRT(DABS(x)),x)
         x = covb*covb - 1.
         r2(m) = DSIGN(DSQRT(DABS(x)),x)
         aag(m)  = 2./(covb*covb)
         aag1(m) = aag(m) - 1.

         if(m.lt.nl) then
            tt(m) = rh1d(m)*cv*cv
            xh(m) = aag(m)*aag(m)
            xh1(m) = aag1(m)*aag1(m)
            hc = pic*h(m)

c----z01, z02 are Pm and Qm in H53
            z01 = hc*r1(m)
            z02 = hc*r2(m)

c----Case of r(alpha) being imaginary
            if(cv.gt.vp1disc(m)) then
               e1(m) = 1.
               c01(m) = DCOS(z01)
               s01(m) = DSIN(z01)
            else
c----for imaginary case, use hyperbolic sin and cosine
               x = DEXP(z01)
               e1(m) = -1.
               c01(m) = 0.5*(x + 1./x)
               s01(m) = 0.5*(x - 1./x)
            end if

c----Case of r(beta) being imaginary
            if(cv.gt.vs1disc(m)) then
              e2(m) = 1.
              c02(m) = DCOS(z02)
              s02(m) = DSIN(z02)
            else
c----for imaginary case, use hyperbolic sin and cosine
              x = DEXP(z02)
              e2(m) = -1.
              c02(m) = 0.5*(x + 1./x)
              s02(m) = 0.5*(x - 1./x)
            end if
c           write(*,*) m, e1(m), c01(m), s01(m), e2(m), c02(m), s02(m)
c            read(*,*) pause
         endif
      enddo

c---Now propagate to the surface to compute the surface displacements

c---The following should be equivalent to the inverse E matrix for
c       the half space (see Eqn 2.16 of H53) but are not exactly.  The
c       first two rows of E inverse are multiplied here by
c       -(alpha*alpha) and the last two rows by -2(beta*beta).  Because
c       we will eventually use displacement ratios, these factors will
c       cancel out, so it should be ok.

      x = -aag1(nl)*cv*cv

      e11 = 2.*vs1disc(nl)*vs1disc(nl)
      e13 = -1./rh1d(nl)
      e22 = x/r1(nl)
      e24 = e13/r1(nl)
      e31 = x/r2(nl)
      e33 = -e13/r2(nl)
      e42 = -e11
      e44 = e13


c       write(*,*) ' nl, x, e11, e33 = ', nl, x, e11, e33
c----The 4 x 2 b matrix will hold the multiplication
c    of all the a matrices.  This initialization allows
c    the first two columns of the a matrix to be saved
c    the first time around.

      b11 = 1.
      b12 = 0.
      b21 = 0.
      b22 = 1.
      b31 = 0.
      b32 = 0.
      b41 = 0.
      b42 = 0.

c----Loop over all layers
      do 7 m = 1,nl-1
         pd = s01(m)/r1(m)
         qd = s02(m)/r2(m)
         pm = e1(m)*s01(m)*r1(m)
         qm = e2(m)*s02(m)*r2(m)
c
c       c1, cp  = cos(Pm)
c       c2, cq  = cos(Qm)
c       s1      = sin(Pm)
c       s2      = sin(Qm)
c  The following lines are a clever way of
c  upgrading the sin and cosine for the next
c  frequency.  Essentally we are evaluating:
c       cos(f + df) = cos(f)cos(df) - sin(f)sin(df)
c       sin(f + df) = sin(f)cos(df) + cos(f)sin(df)
c
         cp = c01(m)
         cq = c02(m)
         tm = tt(m)
         ag = aag(m)
         ag1 = aag1(m)
         x = xh(m)
         y = xh1(m)
c
c  These "a" terms those shown on page 21 of H53
c
         if (vs1disc(m).ne.0.0) then
            a11 = ag*cp - ag1*cq
            a12 = ag1*pd + ag*qm
            a13 = (cq - cp)/tm
            a14 = (pd + qm)/tm
            a21 = -ag*pm - ag1*qd
            a22 = ag*cq - ag1*cp
            a23 = (pm + qd)/tm
            a24 = a13
            a31 = tm*ag*ag1*(cp - cq)
            a32 = tm*(y*pd + x*qm)
            a33 = a22
            a34 = a12
            a41 = tm*(x*pm + y*qd)
            a42 = a31
            a43 = a21
            a44 = a11
         else
c---special forms for fluid layer from eqn 6.3 of H53.  Note there
c   seems to be a sign error on a12 in his paper, however.
            a11 = 0.0
            a12 = -pd
            a13 = -cp/tm
            a14 = 0.0
            a21 = 0.0
            a22 = cp
            a23 = pm/tm
            a24 = 0.0
            a31 = 0.0
            a32 = tm*pd
            a33 = cp
            a34 = 0.0
            a41 = 0.0
            a42 = 0.0
            a43 = 0.0
            a44 = 0.0
         endif
c
c       D is the accumulation of the multiplication of the A
c       matrices as we descend through the layers.  Note that
c       because we start at the surface, the stresses are 0 and
c       so the vector we multiply is (u, w, 0, 0).  This means 
c       that for the entire problem we need save only the first
c       two columns of the multiplication.
c
c       The + and - signs below are a consequence of every other
c       term in A being real and imaginary.  Thus, we multiply:
c               R I R I  R I
c               I R I R  I R
c               R I R I  R I
c               I R I R  I R
c       See H53, page 23 for details.
c
         d11 =  a11*b11 - a12*b21 + a13*b31 - a14*b41
         d12 =  a11*b12 + a12*b22 + a13*b32 + a14*b42
         d21 =  a21*b11 + a22*b21 + a23*b31 + a24*b41
         d22 = -a21*b12 + a22*b22 - a23*b32 + a24*b42
         d31 =  a31*b11 - a32*b21 + a33*b31 - a34*b41
         d32 =  a31*b12 + a32*b22 + a33*b32 + a34*b42
         d41 =  a41*b11 + a42*b21 + a43*b31 + a44*b41
         d42 = -a41*b12 + a42*b22 - a43*b32 + a44*b42
         b11 = d11
         b12 = d12
         b21 = d21
         b22 = d22
         b31 = d31
         b32 = d32
         b41 = d41
         b42 = d42
7     continue
c
c       Now multiply A by inverse E of the last layer (with
c       the scaling exception discussed above.  C is then the
c       same as J in H53 (again with the different scaling).

      c11 = e11*b11 + e13*b31
      c12 = e11*b12 + e13*b32
      c21 = e22*b21 + e24*b41
      c22 = e22*b22 + e24*b42
      c31 = e31*b11 + e33*b31
      c32 = e31*b12 + e33*b32
      c41 = e42*b21 + e44*b41
      c42 = e42*b22 + e44*b42

c---reintroduce the scaling: Note that technically neither these terms
c       nor the "t" scaling below are really necessary because we are
c       only interested in the ratio of U/W in this application, but we include
c       them here anyway for completeness.
c
      scl1 = -1./(vp1disc(nl)*vp1disc(nl))
      scl2 = -1./(2.0*vs1disc(nl)*vs1disc(nl))
      c11 = c11*scl1
      c12 = c12*scl1
      c21 = c21*scl1
      c22 = c22*scl1
      c31 = c31*scl2
      c32 = c32*scl2
      c41 = c41*scl2
      c42 = c42*scl2
c
c       x and y are the real and imaginary parts of D in
c       equation 7 of H62. ag is the magnitude of D squared
c
      x = -c11*c42 + c21*c32 + c31*c22 - c41*c12
      y =  c11*c32 + c21*c42 - c31*c12 - c41*c22
      ag = x*x + y*y

c---t for P waves
      if (ctype.eq.'P') then
         tm = 2.*cv/(vp1disc(nl)*ag)
c---t for S waves
      else
         tm = cv/(vs1disc(nl)*ag)
      endif
c
c       wr, wi are the real and imaginary components of vertical
c               velocity
c       ur, ui are the real and imaginary components of horizontal
c               velocity
c       See equations 5 and 6 of H62.  Note that D = |D| e(i theta), so
c       1/D = e(-i theta)/|D| = (x - iy)/|D|^2 = D*/(DD*)
c
c
c       rs1 and rs2 are the real and imaginary parts of the
c       velocity ratio u/w for P and w/u for S.  Note that
c       u/w = uw*/(ww*) = (ur*wr + ui*wi + i(ui*wr - ur*wi))/ww*
c       w/u = wu*/(uu*) = (wr*ur + wi*ui + i(wi*ur - wr*ui))/uu*
c
c---P terms: The "W" terms are reversed in sign because a postive
c       "L" component will will be in the (-Z, +R) direction as defined by
c       H53.  Note that this does NOT change the sense of Z - it is still 
c       positive down.  
c
c       The idea is that the impulse is (sind, cosd) which would be a pulse in
c       the +X, +Z direction, but we want the response to (+x, -z) or (sind, -cosd)
c       so we revers the sign on the Z component.
c
      if (ctype.eq.'P') then
         wr = tm*( c31*x - c41*y)
         wi = tm*(-c31*y - c41*x)
         ur = tm*( c42*x - c32*y)
         ui = tm*(-c42*y - c32*x)
      else
c---S terms: these are the Haskell conventions on W and U.  A source in the +H
c   direction would be (+Z, +R) in the Haskell sense, so we retain the signs
c   on W and U.
         wr = tm*(-c11*x + c21*y)
         wi = tm*( c11*y + c21*x)
         ur = tm*(-c22*x + c12*y)
         ui = tm*( c22*y + c12*x)
      endif

c---these are the surface displacements
      qrw = wr
      qiw = wi
      qru = ur
      qiu = ui

c       write(*,*) ' wr, wi, ur, ui = ', wr, wi, ur, ui
c       read(*,*) pause

c---now propagate backwards, looping over the grid points
      ibovo  = 1
      ifin = 0

c---copy over response at the surface
      ugrn1f(1) = CMPLX(SNGL(qru), SNGL(qiu))
      wgrn1f(1) = CMPLX(SNGL(qrw), SNGL(qiw))

      do iz = 2, nz
         dept = dz*(iz-1) + zorig
          
c----Loop over the model, setting up layer specific values
c----dep(1) should the top of the model (e.g., zero) so start at dep(2)
         do i = 2,nl
            if (dept.lt.h1d(i)) go to 4
         enddo
         i = nl
4        ibov = i - 1
         dloc = dept - h1d(i-1)
         if (dloc.lt.0.) then
            WRITE(*,*) 'Error in haskgrn dloc is < 0! '
            !WRITE(LUNLOG,*) 'Error in haskgrn dloc is < 0! '
            ierr = 1
            RETURN
         endif
c         write(*,*) ' iz, dept, dloc, ibov = ', iz, dept, dloc, ibov

c----see if we have entered a new layer.  If so, finish off the previous one
         if (ibov.ne.ibovo) then
            ibovo = ibov
            if (ibov.gt.1) then
               ifin = ifin + 1
               hc = pic*h(ifin)
c----z01, z02 are Pm and Qm in H53
               z01 = hc*r1(ifin)
               z02 = hc*r2(ifin)

c----Case of r(alpha) being imaginary
               if (cv.gt.vp1disc(m)) then
                  e1(ifin) = 1.
                  c01(ifin) = DCOS(z01)
                  s01(ifin) = DSIN(z01)
               else
c----for imaginary case, use hyperbolic sin and cosine
                  x = DEXP(z01)
                  e1(ifin) = -1.
                  c01(ifin) = 0.5*(x + 1./x)
                  s01(ifin) = 0.5*(x - 1./x)
               end if
c----Case of r(beta) being imaginary
               if (cv.gt.vs1disc(m)) then
                  e2(ifin) = 1.
                  c02(ifin) = DCOS(z02)
                  s02(ifin) = DSIN(z02)
               else
c----for imainary case, use hyperbolic sin and cosine
                  x = DEXP(z02)
                  e2(ifin) = -1.
                  c02(ifin) = 0.5*(x + 1./x)
                  s02(ifin) = 0.5*(x - 1./x)
               end if
            end if
         end if

c----current layer terms
         hc = pic*dloc
c----z01, z02 are Pm and Qm in H53
         z01 = hc*r1(ibov)
         z02 = hc*r2(ibov)
c----Case of r(alpha) being imaginary
         if (cv.gt.vp1disc(ibov)) then
            e1(ibov) = 1.
            c01(ibov) = DCOS(z01)
            s01(ibov) = DSIN(z01)
         else
c----for imaginary case, use hyperbolic sin and cosine
            x = DEXP(z01)
            e1(ibov) = -1.
            c01(ibov) = 0.5*(x + 1./x)
            s01(ibov) = 0.5*(x - 1./x)
         end if

c----Case of r(beta) being imaginary
         if(cv.gt.vs1disc(ibov)) then
            e2(ibov) = 1.
            c02(ibov) = DCOS(z02)
            s02(ibov) = DSIN(z02)
         else
c----for imaginary case, use hyperbolic sin and cosine
            x = DEXP(z02)
            e2(ibov) = -1.
            c02(ibov) = 0.5*(x + 1./x)
            s02(ibov) = 0.5*(x - 1./x)
         end if

c---The following should be equivalent to the inverse E matrix for
c         the half space (see Eqn 2.16 of H53) but are not exactly.  The
c         first two rows of E inverse are multiplied here by
c         -(alpha*alpha) and the last two rows by -2(beta*beta).  Because
c         we will eventually use displacement ratios, these factors will
c         cancel out, so it should be ok.

         if (ibov.eq.nl) then
            x = -aag1(nl)*cv*cv
            e11 = 2.*vs1disc(nl)*vs1disc(nl)
            e13 = -1./rh1d(nl)
            e22 = x/r1(nl)
            e24 = e13/r1(nl)
            e31 = x/r2(nl)
            e33 = -e13/r2(nl)
            e42 = -e11
            e44 = e13
         endif


c----The 4 x 2 b matrix will hold the multiplication
c    of all the a matrices.  This initialization allows
c    the first two columns of the a matrix to be saved
c    the first time around.

c---this part you only need do when you enter a new layer.  Save
c   results up to ifin then then do one ore iteration for the current layer.

         b11 = 1.
         b12 = 0.
         b21 = 0.
         b22 = 1.
         b31 = 0.
         b32 = 0.
         b41 = 0.
         b42 = 0.

c----Loop over all layers
         do  m = 1,ibov
            pd = s01(m)/r1(m)
            qd = s02(m)/r2(m)
            pm = e1(m)*s01(m)*r1(m)
            qm = e2(m)*s02(m)*r2(m)
c
c       c1, cp  = cos(Pm)
c       c2, cq  = cos(Qm)
c       s1      = sin(Pm)
c       s2      = sin(Qm)
c  The following lines are a clever way of
c  upgrading the sin and cosine for the next
c  frequency.  Essentally we are evaluating:
c       cos(f + df) = cos(f)cos(df) - sin(f)sin(df)
c       sin(f + df) = sin(f)cos(df) + cos(f)sin(df)
c
            cp = c01(m)
            cq = c02(m)
            tm = tt(m)
            ag = aag(m)
            ag1 = aag1(m)
            x = xh(m)
            y = xh1(m)
c
c  These "a" terms those shown on page 21 of H53
c
            if (vs1disc(m).ne.0.0) then
               a11 = ag*cp - ag1*cq
               a12 = ag1*pd + ag*qm
               a13 = (cq - cp)/tm
               a14 = (pd + qm)/tm
               a21 = -ag*pm - ag1*qd
               a22 = ag*cq - ag1*cp
               a23 = (pm + qd)/tm
               a24 = a13
               a31 = tm*ag*ag1*(cp - cq)
               a32 = tm*(y*pd + x*qm)
               a33 = a22
               a34 = a12
               a41 = tm*(x*pm + y*qd)
               a42 = a31
               a43 = a21
               a44 = a11
            else
c---special forms for fluid layer from eqn 6.3 of H53.  Note there
c   seems to be a sign error on a12 in his paper, however.
               a11 = 0.0
               a12 = -pd
               a13 = -cp/tm
               a14 = 0.0
               a21 = 0.0
               a22 = cp
               a23 = pm/tm
               a24 = 0.0
               a31 = 0.0
               a32 = tm*pd
               a33 = cp
               a34 = 0.0
               a41 = 0.0
               a42 = 0.0
               a43 = 0.0
               a44 = 0.0
            endif
c
c       D is the accumulation of the multiplication of the A
c       matrices as we descend through the layers.  Note that
c       because we start at the surface, the stresses are 0 and
c       so the vector we multiply is (u, w, 0, 0).  This means 
c       that for the entire problem we need save only the first
c       two columns of the multiplication.
c
c       The + and - signs below are a consequence of every other
c       term in A being real and imaginary.  Thus, we multiply:
c               R I R I    R I
c               I R I R    I R
c               R I R I    R I
c               I R I R    I R
c       See H53, page 23 for details.
c
            d11 =  a11*b11 - a12*b21 + a13*b31 - a14*b41
            d12 =  a11*b12 + a12*b22 + a13*b32 + a14*b42
            d21 =  a21*b11 + a22*b21 + a23*b31 + a24*b41
            d22 = -a21*b12 + a22*b22 - a23*b32 + a24*b42
            d31 =  a31*b11 - a32*b21 + a33*b31 - a34*b41
            d32 =  a31*b12 + a32*b22 + a33*b32 + a34*b42
            d41 =  a41*b11 + a42*b21 + a43*b31 + a44*b41
            d42 = -a41*b12 + a42*b22 - a43*b32 + a44*b42
            b11 = d11
            b12 = d12
            b21 = d21
            b22 = d22
            b31 = d31
            b32 = d32
            b41 = d41
            b42 = d42
         enddo

         ur =  b11*qru - b12*qiw
         ui =  b11*qiu + b12*qrw
         wr = -b21*qiu + b22*qrw
         wi =  b21*qru + b22*qiw

         ugrn1f(iz) = CMPLX(SNGL(ur), SNGL(ui))
         wgrn1f(iz) = CMPLX(SNGL(wr), SNGL(wi))

      enddo
c---end loop over nodes

      DEALLOCATE(s01)
      DEALLOCATE(c01)
      DEALLOCATE(c02)
      DEALLOCATE(s02)
      DEALLOCATE(aag)
      DEALLOCATE(aag1)
      DEALLOCATE(e1)
      DEALLOCATE(e2)
      DEALLOCATE(h)
      DEALLOCATE(r1)
      DEALLOCATE(r2)
      DEALLOCATE(tt)
      DEALLOCATE(xh)
      DEALLOCATE(xh1)
      return
      end

