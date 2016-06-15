c
c---------------------------------------------------------------------------
c
      SUBROUTINE srcreg(fs, ireg, sspread, dx, dy, dz,
     ;                  xoff, yoff, zoff, ierr)

c***************************************************************************
c Computes source weights for a small region surrounding each source. 
C The region is defined using the array fs(nsgg, nsgg). The half-width
C of the region is defined by sspread. The grid spacing is defined by dx, dz. 
C The offsets xoff, zoff are the distances to the source from centre 
C of the region.
c***************************************************************************
C
C Revision history:
C
C Gerhard Pratt, June 14th, 1998
C Removed this piece of code from srcsub, made it stand alone. This allows
C the source to move around inside the region from shot to shot (ie, in 
C order that the shot be located in-between grid points).
C
C G Hicks. 6/99
C Cleaned up. Normalisation of sum of 'fs' values to 1.0 re-installed, and 
C dimensions of 'fs' checked.
C
C G Hicks. 3/00
C This subroutine now has two possible functions: If sspread>0.5 a 
C distributed ('big') source is employed, using the dirac function to spread
C the source terms over a region of the modelling grid. If sspread=0.5 a 
C point source is employed, using a tapered sinc function to allow the 
C source to be located between nodes. All distances are now in nodes.
C___________________________________________________________________________
C
      IMPLICIT NONE
      INTERFACE
         REAL FUNCTION dirac(h, x, ierr)
         IMPLICIT NONE
         REAL, INTENT(IN) ::  h, x
         INTEGER, INTENT(OUT) :: ierr
         END FUNCTION dirac

         REAL FUNCTION sinc(ireg, x)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ireg
         REAL, INTENT(IN) :: x
         END FUNCTION sinc
      END INTERFACE


c     include 'dimension.inc'
C
C Arguments:
C
      REAL, INTENT(OUT) :: fs(9, 9, 9)
      REAL, INTENT(IN) :: sspread, dx, dy, dz, xoff, yoff, zoff
      INTEGER, INTENT(IN) :: ireg
      INTEGER, INTENT(OUT) :: ierr
C
C local variables:
C
      INTEGER      ix, iy, iz, ixmax, iymax, izmax
      REAL         amp, deltax, deltay, deltaz, xsrc, ysrc, zsrc
C
C Functions:
C  
C
C trc file variables:
C
c      include 'oxford.inc'
c      character*80 filnam
c      integer ierr
c____________________________________________________________________
C      
c executable code:
c--------------------------------------------------------------------
c calculate spatial variation within source region.
      ierr = 0
      if (sspread.le.(0.5)) then

c Use point source:

         ixmax = 2*ireg + 1
         iymax = ixmax ! TODO: check this
         izmax = ixmax
C      if (ixmax.gt.nrgg) then
C         write (*,10) ixmax
C         sTop
C      endif
         xsrc = ireg + 1 + (xoff/dx)
         ysrc = ireg + 1 + (yoff/dy)
         zsrc = ireg + 1 + (zoff/dz)
         do iy = 1, iymax
            deltay = ABS(iy-ysrc)
            do ix = 1, ixmax
               deltax = ABS(ix-xsrc)
               do iz = 1, izmax
                  deltaz = ABS(iz-zsrc)
                  fs(iz,ix,iy) = sinc(ireg,deltaz)
     ;                         * sinc(ireg,deltay)
     ;                         * sinc(ireg,deltax)
               enddo
            enddo
         enddo

      else

c Use distributed source:

         ixmax = INT(4*sspread) + 1
         izmax = ixmax
         iymax = ixmax
C      if (ixmax.gt.nrgg) then
C         write (*,10) ixmax
C         sTop
C      endif
         xsrc = INT(2*sspread) + 1 + (xoff/dx)
         ysrc = INT(2*sspread) + 1 + (yoff/dy)
         zsrc = INT(2*sspread) + 1 + (zoff/dz)
         amp = 0.0
         do iy = 1, iymax
            deltay = ABS(iy-ysrc)
            do ix = 1, ixmax
               deltax = ABS(ix-xsrc)
               do iz = 1, izmax
                  deltaz = ABS(iz-zsrc)
                  fs(iz,ix,iy) = dirac(sspread,deltaz,ierr)
     ;                         * dirac(sspread,deltay,ierr)
     ;                         * dirac(sspread,deltax,ierr)
                  IF (ierr /= 0) RETURN
                  amp = amp + fs(iz,ix,iy)
               enddo
            enddo
         enddo
c Normalize amplitude to 1.0 (this is a kludge)
         do iy=1, iymax
            do ix=1, ixmax
               do iz=1, izmax
                  fs(iz,ix,iy) = fs(iz,ix,iy) / amp
               enddo
            enddo
         enddo

      end if

c Output source region as trc file (for visualization/debug)
C      filnam = 'sregion.trc'
C      call dopen (75, filnam, 2, ixmax, ierr)
C      cpx = .false.
C      call whead(75,title,izmax,ixmax,cpx,deltat,buf,ierr)
C      do iz = 1, izmax
C         write (75,rec=iz+2) (real(fs(iz,ix)),ix=1,ixmax)
C      end do
C      STOP

C      do iz = 1, izmax
C         write(*,*) (real(fs(iz,ix)),ix=1,ixmax)
C      end do

c  10 format(/,' ERROR - NRGG is too small. Must be at least,',I3,
c    &      //,' Check file dimension.inc and recompile!')

      return
      end

c**************************************************************************

      REAL FUNCTION sinc(ireg, x)

c**************************************************************************
c     include 'dimension.inc'
      IMPLICIT NONE
      INTERFACE
         REAL FUNCTION bessi0(x)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         END FUNCTION
      END INTERFACE
      INTEGER, INTENT(IN) :: ireg
      REAL, INTENT(IN) :: x
      REAL    taper,a,lx,b

c Function:
      !REAL    bessi0
      REAL, PARAMETER :: pi = 3.141592653589793

c Evaluates an optimally windowed sinc function which is used as a phase
c shift operator in order to allow sources and recievers to be positioned
c between nodes.
c 'ireg' is the half length of the tapered sinc function (in nodes).
c 'x' is the distance along the x-axis at which to evaluate the function.
c G Hicks 2/3/00

c Form sinc function:
      if (x.le.(1.0e-10)) then
         sinc=1.0
      else
         sinc=SIN(pi*x)/(pi*x)
      end if

c These coefficients define the optimal Kaiser window. See "Arbitrary
c source and receiver positioning in finte-difference schemes, using
c Kaiser windowed sinc functions", Hicks 2000.

c Use this block for optimal upto half the Nyquist wavenumber 
c (ie, 4 nodes per wavelength):
      if (ireg.eq.1) then
        a = 1.24
      else if (ireg.eq.2) then
        a = 2.94
      else if (ireg.eq.3) then
        a = 4.53
      else if (ireg.eq.4) then
        a = 6.31
      else if (ireg.eq.5) then
        a = 7.91
      else if (ireg.eq.6) then
        a = 9.42
      else if (ireg.eq.7) then
        a = 10.95
      else if (ireg.eq.8) then
        a = 12.53
      else if (ireg.eq.9) then
        a = 14.09
      else
        a = 14.18
      endif

c Use this block for optimal upto two thirds the Nyquist wavenumber 
c (ie, 3 nodes per wavelength):
C      if (ireg.eq.1) then
C        a = 0.0
C      else if (ireg.eq.2) then
C        a = 1.84
C      else if (ireg.eq.3) then
C        a = 3.04
C      else if (ireg.eq.4) then
C        a = 4.14
C      else if (ireg.eq.5) then
C        a = 5.26
C      else if (ireg.eq.6) then
C        a = 6.40
C      else if (ireg.eq.7) then
C        a = 7.51
C      else if (ireg.eq.8) then
C        a = 8.56
C      else if (ireg.eq.9) then
C        a = 9.60
C      else
C        a = 10.64
C      endif

c Kaiser windowing:
      lx = x/FLOAT(ireg)
      if (lx.lt.1) then
         b = a * (1.0-lx**2)**0.5
         taper = bessi0(b) / bessi0(a)
      else
         taper = 0.0
      endif

c Apply optimal windowing factor:
      sinc = sinc * taper

      return
      end

c**************************************************************************

      REAL FUNCTION dirac(h, x, ierr)

c**************************************************************************
c discrete approximation of the dirac delta function, where h is the spatial
c half width of the pulse. normally 0.5 < h < 1.5. (alterman & aboudi, 1970;
c m. korn, 1983). x is the absolute distance to the pulse centre.
C 
C     I'm too lazy to rewrite the necessary code, so I'll call 
C     EXIT_MPI here with a failure - BB 08/2010
C
c     COMMON/LOGFILE/LUNLOG
c     INTEGER*4 LUNLOG 
      REAL, INTENT(IN) ::  h, x
      INTEGER, INTENT(OUT) :: ierr
      ierr = 0
      dirac = 0.0
      if (h.le.(1.0e-10)) then
c        WRITE (LUNLOG,1000) H
         write (*,1000) h
 1000    format(/,' dirac - fatal error, spatial extent too small (',
     &            e10.2,').')
         ierr = 1
         RETURN
c        CALL EXIT_MPI('Error in srcreg',LUNLOG,MPIERR)  
      end if

      if (x.lt.h) then
         dirac = ( 1.0 - (x**2/(2.0*h**2)) ) / (2.0*h)
      else if (x.lt.(2.0*h) ) then
         dirac = ( 1.0 + (x**2 - 4.0*x*h)/(4.0*h**2) ) / h
      else
         dirac = 0.0
      end if

      return
      end

c******************************************************************

! TODO: bessel functions are now fortran standard
      REAL FUNCTION bessi0(x)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,
     *-0.1647633d-1,0.392377d-2/
      if (ABS(x).lt.3.75) then
         y=(x/3.75)**2
         bessi0=SNGL(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         ax=ABS(x)
         y=3.75/ax
         bessi0=(EXP(ax)/SQRT(ax))*
     ;         SNGL(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     ;         (q7+y*(q8+y*q9))))))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 'LQ2$0iy.
