      SUBROUTINE READINI_FORWARD(projnm, ierr)
      USE ISO_C_BINDING
      USE MODEL_MODULE, ONLY : freq, freesurf, azm, dx, dy, dz, freqbase, &
                               xorig, yorig, zorig, nom, nx, ny, nz
      USE PML_MODULE, ONLY : pmld, pmlf, pmlr, ipml
      USE INTERFACE_MODULE, ONLY : iniparser_getbool, iniparser_getint, &
                                   iniparser_getfloat, iniparser_getint, &
                                   iniparser_init, iniparser_finalize
      IMPLICIT NONE
      CHARACTER(256), INTENT(IN) :: projnm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(256) var
      CHARACTER(6) ck
      INTEGER k
      LOGICAL(C_BOOL), PARAMETER :: true = .TRUE.
      LOGICAL(C_BOOL), PARAMETER :: false = .FALSE.
      ierr = 0
      ! set the ini file name (notice iniparser is C program and must be null terminated)
!     CALL SET_C_STRING(TRIM(ADJUSTL(projnm))//'.ini', inifl)
!     DO i=1,256
!        inifl(i) = cwork(i:i)
!     ENDDO
      ierr = INIPARSER_INIT(TRIM(ADJUSTL(projnm))//'.ini'//CHAR(0))
      IF (ierr /= 0) THEN
         WRITE(*,*) 'readini_forward: Error reading ini file: ', &
                    TRIM(ADJUSTL(projnm))//'.ini'
         RETURN
      ENDIF
      ! get the pertinent variables for forward modeling
      nx = INIPARSER_GETINT('model:nx'//CHAR(0), 0, ierr)
      ny = INIPARSER_GETINT('model:ny'//CHAR(0), 0, ierr)
      nz = INIPARSER_GETINT('model:nz'//CHAR(0), 0, ierr)
      IF (nx < 1 .OR. ny < 1 .OR. nz < 1) THEN
         WRITE(*,*) 'readini_forward: Invalid nx, ny, or nz:', nx, ny, nz
         ierr = 1
         RETURN
      ENDIF
      ! model origin
      xorig = INIPARSER_GETFLOAT('model:xorig'//CHAR(0), 0.0, ierr)
      yorig = INIPARSER_GETFLOAT('model:yorig'//CHAR(0), 0.0, ierr)
      zorig = INIPARSER_GETFLOAT('model:zorig'//CHAR(0), 0.0, ierr) 
      ! grid spacing
      dx = INIPARSER_GETFLOAT('model:dx'//CHAR(0), 0.0, ierr)
      dy = INIPARSER_GETFLOAT('model:dy'//CHAR(0), dx, ierr)
      dz = INIPARSER_GETFLOAT('model:dz'//CHAR(0), dx, ierr)
      IF (dx <= 0.0 .OR. dy <= 0.0 .OR. dz <= 0.0) THEN
         WRITE(*,*) 'readini_forward: Invalid dx, dy, or dz:', dx, dy, dz
         ierr = 1
         RETURN
      ENDIF
      IF (dx /= dy .OR. dx /= dz) THEN
         WRITE(*,*) 'readini_forward: Irregular grid spacing not permitted'
         ierr = 1
         RETURN
      ENDIF
      ! model azimuth
      azm = INIPARSER_GETFLOAT('model:azm'//CHAR(0), 0.0, ierr)
      ! attenuation
      freqbase = INIPARSER_GETFLOAT('model:freqbase'//CHAR(0), 20.0, ierr)
      ! free surface/boundary condition information
      freesurf(:) = .FALSE.
      freesurf(1) = INIPARSER_GETBOOL('model:fst'//CHAR(0), true,  ierr)
      freesurf(2) = INIPARSER_GETBOOL('model:fsr'//CHAR(0), false, ierr)
      freesurf(3) = INIPARSER_GETBOOL('model:fsb'//CHAR(0), false, ierr) 
      freesurf(4) = INIPARSER_GETBOOL('model:fsl'//CHAR(0), false, ierr)
      freesurf(5) = INIPARSER_GETBOOL('model:fsm'//CHAR(0), false, ierr)
      freesurf(6) = INIPARSER_GETBOOL('model:fsp'//CHAR(0), false, ierr)
      ! PML information
      pmlr = INIPARSER_GETFLOAT('pml:pmlr'//CHAR(0), 0.0010, ierr) 
      pmld = INIPARSER_GETFLOAT('pml:pmld'//CHAR(0), 10.0*dx, ierr)
      ipml = INT(pmld/dx) + 1
      pmlf = 3.0*LOG(1.0/pmlr)/(2.0*pmld*pmld*pmld)
      ! read the modeling frequencies 
      nom = INIPARSER_GETINT('model:nom'//CHAR(0), 0, ierr)
      IF (nom < 1) THEN
         WRITE(*,*) 'readini_forward: No frequencies to model1'
         ierr = 1
         RETURN
      ENDIF
      ALLOCATE(freq(nom))
      freq(:) = 0.0
      DO k=1,nom
         var(:) = ' '
         ck(:) = ' '
         WRITE(ck, '(I6)') k
         var = 'model:freq_'//TRIM(ADJUSTL(ck))//CHAR(0)
         freq(k) = INIPARSER_GETFLOAT(var, 0.0, ierr)
         IF (ierr /= 0 .OR. freq(k) <= 0.0) THEN
            WRITE(*,*) 'readini_forward: Error reading frequency', k
            ierr = 1
            RETURN
         ENDIF
      ENDDO
      ! read the receivers
      
      ! free the memory associated with iniparser
      CALL iniparser_finalize() 
      RETURN
      END 

