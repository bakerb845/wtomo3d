      SUBROUTINE INFO(nx, ny, nz, ns, nr, ng, nom)
c
c Provide the user with basic information about the run and the compiled
c version of the software.  
c 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, ny, nz, ns, nr, ng, nom
      REAL array, ressz
    
      !ressz = float(imax1*imax2*imax3*8)/1000000.
      array = FLOAT(6*nx*ny*nz*4 + 4*nx*ny*nz*8)/1000000.0
      ressz = FLOAT(3*nom*ng*ns*8)/1000000.0
c     write (lunlog, 2) nx,nx, ny,ny, nz,nz, ns,ns, nr,nr, ng,nr,
c    ;                  nom, nom, ressz
      write (*, 2) nx, nx, ny,ny, nz,nz, ns,ns, nr,nr, ng,ng,
     ;             nom, nom, array, ressz

2     format(  /,15X,' ******************************************** ',
     &         /,15X,' *                                          * ',
     &         /,15X,' * Run time information:                    * ',
     &         /,15X,' *                                          * ',
     &         /,15X,' * Dimensions       .ini file      compiler * ',
     &         /,15X,' *                                          * ',
     &         /,15X,' * NX               ',I5,'        ',I5,'      * ',
     &         /,15X,' * NY               ',I5,'        ',I5,'      * ',
     &         /,15X,' * NZ               ',I5,'        ',I5,'      * ',
     &         /,15X,' * NS (sources)     ',I5,'        ',I5,'      * ',
     &         /,15X,' * NR (hydrophones) ',I5,'        ',I5,'      * ',
     &         /,15X,' * NG (geophones)   ',I5,'        ',I5,'      * ',
     &         /,15X,' * NOM(frequencies) ',I5,'        ',I5,'      * ',
     &         /,15X,' *                                          * ',
     &         /,15X,' * Model arrays:  ',F12.4,' Mbytes       * ',
     &         /,15X,' * Results array: ',F12.4,' Mbytes       * ',
     &         /,15X,' *                                          * ',
     &         /,15X,' ******************************************** ')

      return
      end 
