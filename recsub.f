        subroutine recsub(projnm,
     ;                    ncom, nsam, nom, ngt,
     ;                    recin,
     ;                    deltatt,freq,
     ;                    receiver, ierr)  
c 
c     Reads the receiver time functions 
c 
c     Input      Meaning 
c     -----      ------- 
c     deltatt    sampling rate (ms) 
c     freq       frequency list 
c     iunit      unit number of file for reading 
c     lunlog     log unit number
c     modrecp    tells how to update receiver in inversion; likely=-1 
c     ncom       comment level
c     nom        number of frequencies 
c     nsam       number of samples
c     ngt        number of receivers from .ini file 
c     recin      =1 then time domain receiver; .rec file 
c                =2 then a frequency domain receiver; .gs file
c                =3 then a delta function (white)
c     tau        anti-aliasing number 
c 
c     output     meaning 
c     ------     ------- 
c     ierr       != 0 then an error was encountered 
c     receiver   frequency domain receiver response for all receivers 
c 
      IMPLICIT NONE

c       include 'dimension.inc'

      CHARACTER(*), INTENT(IN) :: projnm 
      REAL, DIMENSION(:), INTENT(IN) :: freq
      REAL, INTENT(IN) :: deltatt
      INTEGER, INTENT(IN) :: nom, ngt, nsam, recin, ncom 
      COMPLEX, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: receiver
      INTEGER, INTENT(OUT) :: ierr 

      REAL, ALLOCATABLE :: trace(:,:) 
      COMPLEX omega, cfact, receiver1
      REAL  dt, startt, time, zheng, xiang,
     ;      deltat, dtms
      INTEGER ng, nguse, i, j, ifreq, nomin
      REAL, PARAMETER :: small = 1.0e-5
      REAL, PARAMETER :: pi = 3.141592653589793
 
c----------------------------------------------------------------------c
 
      IF (.NOT.ALLOCATED(receiver)) ALLOCATE(receiver(ngt, nom))
      receiver(:,:) = CMPLX(1.0, 0.0)
c.... read appropriate receiver file
      ierr = 0 
      if (ncom.ge.3) then 
         write(*,*)      ''
         write(*,*)      'RECSUB: '
c        write(lunlog,*) ''
c        write(lunlog,*) 'RECSUB: '
      endif
      if (recin.eq.1) then 
         WRITE(*,*) 'recin == 1 not yet done'
         ierr = 1
         RETURN
c....... read # of receiver, # of samples, s-rate, start time (ms)
         allocate(trace(nsam,ngt)) 
         ng = ngt
c....... use the same procedure to read receivers as sources
         print *, len_trim(projnm)
c        call rdsrc_c(flname, nsam, ng, nsam, dt, startt ,trace, ierr)
c-------we use microseconds in the rec file
         dtms = dt/1000. !convert to milli-seconds 
         deltat = deltatt 
         if (dtms.ne.deltatt) then 
            write(*,*)      'RECSUB Error: Sampling rates do not match!'
c           write(lunlog,*) 'RECSUB Error: Sampling rates do not match!'
            write(*,*) ' From file: ', dtms,' from ini: ', deltatt
            ierr = 1 
            return 
         endif 
         nguse = ng  

         if (ng.ne.ngt) then 
            write(*,*)'RECSUB Warning: Number of receivers doesnt match'
c           write(lunlog,*) 'RECSUB Warning: Number of receivers doesnt match'
            if (ng.lt.ngt) then 
               nguse = 1 
               write(*,*) '       Using first record for all receivers'
c              write(lunlog,*) '       Using first record for all receivers' 
            else 
               nguse = ngt 
               write(*,*) '       Receiver exceeds .ini file'
c              write(lunlog,*) '       Receiver file has more receivers than .ini file'
               write(*,*) '       Will use first',ng,' receiver records'
c              write(lunlog,*) '       Will use first ',ng,' receiverrecords'
            endif 
         endif 

         if (startt.ne.0.) then 
            write(*,*) '       RECSUB Warning: Receiver starts at ',
     ;                 startt,' (s)'
c           write(lunlog,*) '       RECSUB Warning: Receiver starts at ',startt,' (s)'
            write(*,*)      '       Setting startt = 0'
c           write(lunlog,*) '       Setting startt = 0'
            startt = 0. 
         endif 
 
c....... loop over the number of receivers we're using
         do i = 1,nguse 
 
c.......... begin discrete fourier transform 
            do ifreq = 1,nom
               receiver(i,ifreq) = CMPLX(0.,0.) 
               omega = 2.*pi*freq(ifreq)

               do j = 1,nsam 
                  time = FLOAT(j - 1)*deltat/1000. + (startt/1000.0)
                  cfact = CEXP(CMPLX(0.0,1.0)*(omega*time))
                  cfact = cfact/FLOAT(nsam) 
                  receiver(i,ifreq) = receiver(i,ifreq)
     ;                              + trace(j,i)*cfact
               enddo

               if (ifreq.eq.1) then
c                 write(lunlog,900) i,omega/(2.*pi),receiver(i,ifreq)  
                  if (ncom.ge.3)
     ;            write(*,900) i,omega/(2.*pi),receiver(i,ifreq)
900               format(/,'     Receiver ',i4,
     ;                   /,'     Frequency ',e10.4,1x,e10.4,1x,
     ;                   /,'     Hz, Freq. computed',e13.5,1x,e13.5)
               else 
c                 write(lunlog,901) i,omega/(2.*pi),receiver(i,ifreq)  
                  if (ncom.ge.3)
     ;            write(*,901) i,omega/(2.*pi),receiver(i,ifreq)
901               format(/,'     Receiver ',i4,
     ;                   /,'     Frequency ',e10.4,1x,e10.4,1x,
     ;                   /,'     Hz, Freq. computed',e13.5,1x,e13.5)
               endif 
            enddo       !end loop on frequencies 
            if (ncom.ge.3) write(*,*) ''  
         enddo         !end loop on receivers 

         deallocate(trace) 
 
c.... input a frequency receiver component 
      elseif (recin.eq.2) then  
         ng = ngt
         write(*,*) 'RECSUB Error: IO not yet programmed'
         ierr = 1 
         return
         if (nomin.ne.nom) then 
            write(*,*)      'RECSUB error: frequency mismatch'
c           write(lunlog,*) 'RECSUB error: frequency mismatch'
            ierr = 1 
            return 
         endif 
         do ifreq = 1, nom 
            zheng = CABS(receiver1) 
            xiang = ATAN2(AIMAG(receiver1),REAL(receiver1))
            do j = 1, ng
              receiver(j,ifreq) = receiver1
            enddo
         enddo
 
c.... Greens Function computation 
      elseif (recin.eq.3) then 
          ng = ngt
          do ifreq = 1,nom 
             do j = 1,ngt 
                receiver(j,ifreq) = CMPLX(1.0,0.0)
             enddo
          enddo
      else 
          write(*,*)      'RECSUB Error: Invalid choice for recin'
c         write(lunlog,*) 'RECSUB Error: Invalid choice for recin'
          ierr = 1
          return 
      endif 

      return 

c  60 CONTINUE
      ierr = 1 
      write(*,*)      'RECSUB Error: Premature end of receiver file'
c     write(lunlog,*) 'RECSUB Error: Premature end of receiver file'

      return 
      end  
c                                                                      c
c======================================================================c
