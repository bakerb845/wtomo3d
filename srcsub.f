        SUBROUTINE srcsub(projnm,
     ;                    ncom, nsam, nom, nomi,
     ;                    nst, modsrcp,
     ;                    srcin, nsg, isg,
     ;                    deltatt, tau, freq,
     ;                    dsd, source, ierr) 
c 
c     Reads the source time functions 
c 
c     Input      Meaning 
c     -----      ------- 
c     deltatt    sampling rate (ms) 
c     freq       frequency list 
c     freq_ex    extrapolation frequency
c     isg        holds source group numbers
c     lunlog     log unit number
c     modsrcp    tells how to update source in inversion; likely=-1 
c     ncom       comment level
c     nom        number of frequencies 
c     nsam       number of samples
c     nst        number of sources from .ini file 
c     srcin      =1 then time domain source; .src file 
c                =2 then a frequency domain source; .gs file
c                =3 then a delta function (white spectrum)
c     syndat     =0 strongly recommend 
c                =1 assume input data are temporal delta fn convolved w/
c                 travel times 
c                =2 extrapolate response from nearby frequencies, 
c     tau        anti-aliasing number 
c 
c     output     meaning 
c     ------     ------- 
c     ierr       != 0 then an error was encountered 
c     dsd 
c     source     frequency domain source response for all sources 
c 
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: projnm 
      REAL, DIMENSION(:), INTENT(IN) :: freq !(nom)
      INTEGER, DIMENSION(:), INTENT(INOUT) :: isg
      REAL, INTENT(IN) :: deltatt, tau
      INTEGER, INTENT(IN) :: nsg, nom, nomi, nst, nsam,
     ;                       srcin, modsrcp, ncom!, lunlog
      COMPLEX, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: source
      REAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: dsd !(nsg,*) 
      INTEGER, INTENT(OUT) :: ierr

c     integer*4       isg(*), nsg, iunit, nom, nst, nsam
c     integer*4       syndat, srcin, modsrcp, ncom, lunlog
c     integer*4       ierr 

c       complex*8       source(nst,*) 

      REAL, ALLOCATABLE :: trace(:,:) 
      COMPLEX omega, cfact, source1
      REAL dt, startt, time, zheng, xiang, deltat, dtms
      INTEGER ns, nsuse, iii, i, j, kg, ifreq, nomin, nsgin, lenos
      INTEGER itest
      REAL st
      REAL, PARAMETER :: small = 1.0e-5
      REAL, PARAMETER :: pi = 3.141592653589793
 
c----------------------------------------------------------------------c
 
c.... read appropriate source file
      ierr = 0 
      IF (.NOT.ALLOCATED(source)) ALLOCATE(source(nsg, nom))
      IF (.NOT.ALLOCATED(dsd)) ALLOCATE(dsd(nsg,nom))
      source(:,:) = CMPLX(1.0, 0.0)
      dsd(:,:) = 1.0
      if (ncom.ge.3) then 
         write(*,*)      ''
         write(*,*)      'srcsub: '
c        write(lunlog,*) ''
c        write(lunlog,*) 'srcsub: '
      endif
      if (srcin.eq.1) then 
         WRITE(*,*) 'srcin == 1 not yet done'
         ierr = 1
         RETURN
c....... read # of source, # of samples, s-rate, start time (ms)
         allocate(trace(nsam,nst)) 
         trace(:,:) = 0.0
         ns = nst
         lenos = LEN_TRIM(projnm)
c        call dblank(80,projnm,lenos)
         write(*,*)      'srcsub: Opening ', projnm(1:lenos)//'src.h5'
c        write(lunlog,*) 'srcsub: Opening ', flname(1:lenos)//'src.h5'
         !call rdsrc_c(flname, nsam, ns, nsam, dt, startt, trace, ierr)
         dtms = dt/1000. !dt comes in as micro-secs; convert to milli-seconds 
         deltat = deltatt 
         if (dtms.ne.deltatt) then 
            write(*,*)      'srcsub Error: Sampling rates do not match!'
c           write(lunlog,*) 'srcsub Error: Sampling rates do not match!'
            write(*,*)      ' From file: ', dtms,' from ini: ', deltatt
c           write(lunlog,*) ' From file: ', dtms,' from ini: ', deltatt
            ierr = 1 
            return 
         endif 
         nsuse = ns  

         if (ns.ne.nst) then 
            write(*,*)  'srcsub Warning: Number of sources doesnt match'
c           write(lunlog,*) 'srcsub Warning: Number of sources doesnt match'
            if (ns.lt.nst) then 
               nsuse = 1 
               write(*,*) '       Using first record for all sources'
c              write(lunlog,*) '       Using first record for all sources' 
            else 
               nsuse = nst 
               write(*,*) '       Number of sources exceeds ini file'
c              write(lunlog,*) '       Source file has more sources than .ini file'
               write(*,*) '       Will use first ',ns,' source records'
c             write(lunlog,*) '       Will use first ',ns,' source records'
            endif 
         endif 

         if (startt.ne.0.) then 
            write(*,*) '       srcsub Warning: Source starts at ',
     ;                  startt,' (s)'
c           write(lunlog,*) '       srcsub Warning: Source starts at ',startt,' (s)'
            write(*,*)      '       Setting startt = 0'
c           write(lunlog,*) '       Setting startt = 0'
            startt = 0. 
         endif 
 
c....... loop over the number of sources we're using
         do 1 i=1,nsuse 
c.......... not implemented (?) but this would allow extrapolation from 
c.......... nearby frequencies 
!           if (syndat.eq.2) then 
!              omega_ex = 2.*pi*freq_ex
!              source_ex = CMPLX(0.0,0.0)
!              do j = 1,nsam 
!                 time = FLOAT(j - 1)*deltat/1000. + startt/1000.0
!                 cfact = CEXP(CMPLX(0.0,1.0)*(omega_ex*time))
!                 cfact = cfact/FLOAT(nsam)
!                 source_ex = source_ex + trace(j,i)*cfact
!              enddo
!           endif 
 
c.......... begin discrete fourier transform 
            do 4 ifreq=1,nom
               source(i,ifreq) = cmplx(0.,0.) 
               omega = 2.*pi*freq(ifreq)

c............ anti aliasing 
               if (tau.ne.999.999) omega = omega + CMPLX(0.0,1.0/tau)

               do j=1,nsam 
                  time = FLOAT(j - 1)*deltat/1000. + (startt/1000.0)
                  cfact = CEXP(CMPLX(0.0,1.0)*(omega*time))
                  cfact = cfact/FLOAT(nsam) 
                  source(i,ifreq) = source(i,ifreq) + trace(j,i)*cfact
               enddo
 
c............. modification to source 
!              if (syndat.eq.1) then !delta fn. 
!                 time = 0. 
!                 source(i,ifreq) = CEXP(CMPLX(0.0,1.0)*omega*time) 
!              elseif (syndat.eq.2) then !extrapolate 
!                 twidth = 1.0/freq_ex
!                 fdiff = freq(ifreq) - freq_ex
!                 omega_ex = 2.0*pi*freq_ex
!                 if (tau.ne.999.999)
!    ;            omega_ex = omega_ex + CMPLX(0.0,1.0/tau) 
!                 omega_diff = omega - omega_ex
!                 time = twidth/2.0
!                 cfact = CEXP(CMPLX(0.0,1.0)*omega_diff*time)
!                 if (CABS(omega_diff).ge.small) then
!                    rfact = 2.0/omega_diff*SIN(omega_diff*twidth/2.0)
!                    !normalize so rfact is close to a sinc
!                    rfact = rfact/twidth 
!                 else 
!                    rfact = 1. 
!                 endif 
!                 source(i,ifreq) = source_ex*cfact*rfact
!              endif 

               if (ifreq.eq.1) then
c                 write(lunlog,900) i,omega/(2.*pi),source(i,ifreq)  
                  IF (ncom.ge.3)
     ;            WRITE(*,900) i, omega/(2.*pi), source(i,ifreq)
900               format(/,'     Source ',i4,
     ;                   /,'     Frequency ',e10.4,1x,e10.4,1x,
     ;                   /,'     Hz, Freq. computed',e13.5,1x,e13.5)
               else 
c                 write(lunlog,901) i,omega/(2.*pi),source(i,ifreq)  
                  IF (ncom.ge.3) 
     ;            WRITE(*,901) i,omega/(2.*pi),source(i,ifreq)
901               format(/,'     Source ',i4,
     ;                   /,'     Frequency ',e10.4,1x,e10.4,1x,
     ;                   /,'     Hz, Freq. computed',e13.5,1x,e13.5)
               endif 
4           continue !end loop on frequencies 
c---test to recreate the source
            itest = 0
            if (itest.eq.1) then
               do j = 1, nsam
                  time = FLOAT(j-1)*deltatt/1000.
                  st = 0.0
                  do ifreq = 1, nom
                     omega = 2.0*pi*freq(ifreq)
                     cfact = CEXP(CMPLX(0.0,-1.0)*(omega*time))
                     st = st + REAL(2.0*source(i,ifreq)*cfact)
                  enddo
                  write(68,*) time, st
               enddo
               stop
            endif

            if (ncom.ge.3) write(*,*) ''  
1        continue !end loop on sources 

         deallocate(trace) 
 
c.... input a frequency source component 
      elseif (srcin.eq.2) then  
         write(*,*) 'srcsub Error: IO not yet programmed'
         ierr = 1 
         return
         !read(iunit,end=60) nomin, nsgin 
         if (nsgin.ne.nsg) then 
            write(*,*) 'srcsub Error: Number of source group mismatch'
c           write(lunlog,*) 'srcsub Error: Number of source group mismatch'
            ierr = 1 
            return 
         endif 
         if (nomin.ne.nom) then 
            write(*,*) 'srcsub error: frequency mismatch'
c           write(lunlog,*) 'srcsub error: frequency mismatch'
            ierr = 1 
            return 
         endif 
         isg(nsg + 1) = ns + 1  
         do ifreq=1,nom 
            do kg=1,nsg 
               !read(iunit,end=60) source1 
               iii = iii + 1 
               zheng = CABS(source1) 
               xiang = ATAN2(AIMAG(source1),REAL(source1))
               do j=isg(kg),isg(kg+1)-1 
                  source(j,ifreq) = source1
               enddo
            enddo
         enddo
c 
c.... Greens Function computation 
      elseif (srcin.eq.3) then 
         ns = nst
         write(*,*)      'srcsub: Assigning delta function to STF'
c        write(lunlog,*) 'srcsub: Assigning delta function to STF'
         do ifreq = 1,nom 
            do j = 1,ns 
               source(j,ifreq) = CMPLX(1.0,0.0)
            enddo
         enddo
      else 
         write(*,*)      'srcsub Error: Invalid choice for srcin'
c        write(lunlog,*) 'srcsub Error: Invalid choice for srcin'
         ierr = 1
         return 
      endif 

      return 

      isg(nsg+1) = ns + 1 
 
c.... calculate 'dsd' to compute weighted misfit fn, if the source is
c.... not going to be updated 
      if (modsrcp.eq.0) then 
         if (nom.gt.nomi) then 
            write(*,*)      'srcsub Error: nomi < nom'
c           write(lunlog,*) 'srcsub Error: nomi < nom'
            ierr = 1 
            return 
         endif 
         do ifreq=1,nom 
            do kg=1,nsg 
               dsd(kg,ifreq) = (CABS(source(1,1))/
     ;                          CABS(source(isg(kg),ifreq)))**2
            enddo
         enddo
      endif 

      return 

c  60 CONTINUE
      ierr = 1 
      write(*,*)      'srcsub Error: Premature end of source file'
c     write(lunlog,*) 'srcsub Error: Premature end of source file'

      return 
      end  
c                                                                      c
c======================================================================c
