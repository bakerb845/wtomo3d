
      SUBROUTINE dispersion()
c************************************************************************
c
c    This routine contains the parts of modsub that were relevant to
c    dispersion calculation of wavespeeds in the case of finite Q
c
c    Generally, this routine will be called when Qp and Qs are read in as
c    separate files.   If viscous medium without dispersion is desired,
c    read in Vpi and Vsi instead.
c
c    This is the dispersion model from Aki and Richards, P177. 
c    -RGP March 2000
c
c    This algorithm uses the relation in equation 5.88 of A&R on page 182:
c
c       c(w) = c(1) (1 + 1/(piQ) ln(w/2pi) - i/2Q)
c
c    were c(1) is the wavespeed at 1 hz (w = 2pi).  A&R start with a plane
c    wave form of exp(i(Kx - wt)), which is the same as the default in Pratt,
c    so we retain the signs on the imaginary part here.   But note that
c    Haskell uses the reverse, so we need to conjugate for the background.
c    We do that in calc1d. - SWR 10/2011
c 
c************************************************************************
      USE INIT_MODULE, ONLY : ncom, freqbase, nx, ny, nz
      USE MODEL_MODULE, ONLY : da, mu, qs, qp, vsr, vpr, rho, omega 
      IMPLICIT NONE

c---need to check which of these are really necessary
c     include 'dimension.inc'
c     include 'init.inc'
c     include 'common.inc'
c     include 'oxford.inc'

c     local variables:
        
      INTEGER        ix, iy, iz
      REAL           fact
      REAL           pdiffmax, pdiffmin 
      REAL           sdiffmax, sdiffmin 
      REAL           vprmax, vprmin, vsrmax, vsrmin
      REAL           fac1
      REAL           omegabase
      REAL           omega1
      REAL    vsrold, vsrnew, vsinew
      REAL    vprold, vprnew, vpinew
      COMPLEX        vc

      REAL, PARAMETER :: eps = EPSILON(1.0)*10.0 !.1e-5
      REAL, PARAMETER :: pi = 3.141592653589793

c executable code:
      if(ncom.ge.3)write(*,*)'=== subroutine dispersion ==='

c Remove stretch from nz if anisotropy (restored before return).   
c This has not yet been adapted to elastic case, so next 3 lines
c commented out
c     if (aniso.gt.eps) then
c        nz = nshrink(nz, aniso)
c     end if

      pdiffmax = -1e30
      pdiffmin =  1e30
      vprmax = -1e30
      vprmin =  1e30

      sdiffmax = -1e30
      sdiffmin =  1e30
      vsrmax = -1e30
      vsrmin =  1e30

      omegabase = freqbase * 2.0 * pi
      if  (omegabase.le.0.0.and.ncom.ge.1) 
     +  print*,'Warning omegabase less than zero - no dispersion.'

c---Dont forget, if there is anti-time aliasing, omega is complex
      omega1 = REAL(omega)
      if (omega1.gt.eps.and.omegabase.gt.0.0) then
         fac1 = LOG(omega1/omegabase)
         do iy = 1, ny
            do ix = 1, nx
               do iz = 1,nz
c----S section
                  vsrold = vsr(iz, ix, iy)
                  fact = 1.+((qs(iz,ix,iy)/pi)*fac1)
c---Avoid ridiculously low velocity, especially negative velocities
                  if (fact.lt.0.1) fact = 0.1
                  vsrnew = vsrold*fact
                  sdiffmax = MAX(sdiffmax, vsrold-vsrnew)
                  sdiffmin = MIN(sdiffmin, vsrold-vsrnew)
                  vsrmax = MAX(vsrmax, vsrnew)
                  vsrmin = MIN(vsrmin, vsrnew)
                  vsinew = -vsrnew*qs(iz,ix,iy)/2.0
                  vc = CMPLX(vsrnew, vsinew)
                  mu(iz,ix,iy) = vc*vc*rho(iz,ix,iy)
c----P section
                  vprold = vpr(iz, ix, iy)
                  fact = 1.+((qp(iz,ix,iy)/pi)*fac1)
c---Avoid ridiculously low velocity, especially negative velocities
                  if (fact.lt.0.1) fact = 0.1
                  vprnew = vpr(iz,ix,iy)*fact
                  pdiffmax = MAX(pdiffmax, vprold-vprnew)
                  pdiffmin = MIN(pdiffmin, vprold-vprnew)
                  vprmax = MAX(vprmax, vprnew)
                  vprmin = MIN(vprmin, vprnew)
                  vpinew = -vprnew*qp(iz,ix,iy)/2.0
                  vc = CMPLX(vprnew, vpinew)
                  da(iz,ix,iy) = vc*vc*rho(iz,ix,iy) - 2.0*mu(iz,ix,iy)
               enddo
            enddo
         enddo
      endif

      if (ncom.ge.1) then
         print*, ' Frequency omega, omegabase ', omega1, omegabase
         print*, ' fac1, fact ', fac1, fact
         print*, ' New (dispersed) maximum P velocity; ',vprmax
         print*, ' New (dispersed) minimum P velocity; ',vprmin
         print*, ' New (dispersed) maximum S velocity; ',vsrmax
         print*, ' New (dispersed) minimum S velocity; ',vsrmin
         print*, ' Old velocity - new velocity: '
         print*, ' Maximum difference in P velocity at this frequency ',
     ;           pdiffmax
         print*, ' Minimum difference in P velocity at this frequency ',
     ;           pdiffmin
         print*, ' Maximum difference in S velocity at this frequency ',
     ;           sdiffmax
         print*, ' Minimum difference in S velocity at this frequency ',
     ;           sdiffmin
      end if

      if(ncom.ge.3)write(*,*)'=== FIN dispersion ==='

C Before return, the model must be interpolated onto stretched coordinate system.
C Also restore stretch to nz if anisotropy (restored before return)
C  Not implimented for elastic case so ignored for now
c      if (aniso.gt.eps) then
c         nzsml = nz
c         nz = nstretch(nz, aniso)
c         call stretch(rho, nx, nz, nzsml, aniso)
c         call cstretch(m, nx, nz, nzsml, aniso)
c      end if

      return
      end

