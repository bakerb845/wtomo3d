      SUBROUTINE READINI_FORWARD(projnm, ierr)
      USE ISO_C_BINDING
      USE INIT_MODULE, ONLY : tau, ncom, infspace, totfld, usemin
      USE MODEL_MODULE, ONLY : freq, freesurf, azm, dx, dy, dz, freqbase, &
                               vpvs, xorig, yorig, zorig, nom, nomi, nx, ny, nz
      USE RECPRM_MODULE, ONLY : gwght, rwght, xg, yg, zg, xr, yr, zr, &
                                igreg, irreg, modrecp, gspread, rspread, nr, ng, &
                                recin, usegwt, userwt
      USE SRCPRM_MODULE, ONLY : isg, isreg, modsrcp, ns, nsg, srcin, sspread
      USE PML_MODULE, ONLY : pml, pmld, pmlf, pmlr, ipml
      USE INTERFACE_MODULE, ONLY : iniparser_getbool, iniparser_getchar, &
                                   iniparser_getfloat, iniparser_getint, &
                                   iniparser_getRecvWXYZ, &
                                   iniparser_init, iniparser_finalize
      IMPLICIT NONE
      CHARACTER(256), INTENT(IN) :: projnm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(256) var, cvar
      CHARACTER(6) ck
      INTEGER ig, ir, is, k
      INTEGER, PARAMETER :: nrgg = 9
      LOGICAL(C_BOOL), PARAMETER :: true = .TRUE.
      LOGICAL(C_BOOL), PARAMETER :: false = .FALSE.
      !----------------------------------------------------------------------------------!
      !
      ! parse the ini file 
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
      ! Vp/Vs ratio (presumably used in inversion)
      vpvs = INIPARSER_GETFLOAT('model:vpvs'//CHAR(0), SQRT(3.0), ierr)
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
      pml = INIPARSER_GETBOOL('pml:pml'//CHAR(0), true, ierr)
      IF (.NOT.pml) THEN
         WRITE(*,*) 'readini_forward: PMLs only!'
         ierr = 1
         RETURN
      ENDIF
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
      ! read the geophones - these are 3c displacement
      ng = INIPARSER_GETINT('receivers:ng'//CHAR(0), 0, ierr)
      IF (ng > 0) THEN
         ALLOCATE(gwght(ng))
         ALLOCATE(xg(ng))
         ALLOCATE(yg(ng))
         ALLOCATE(zg(ng))
         DO ig=1,ng
            cvar(:) = ' '
            var(:) = ' '
            ck(:) = ' '
            WRITE(ck, '(I6)') ig 
            var = 'receivers:geo_'//TRIM(ADJUSTL(ck))//CHAR(0)
            ierr = iniparser_getRecvWXYZ(var, gwght(ig), xg(ig), yg(ig), zg(ig))
            IF (ierr /= 0) THEN
               WRITE(*,*) 'readini_forward: Error reading geophone', ig
               RETURN
            ENDIF
         ENDDO
         IF (ABS(xorig) > 0.0 .OR. ABS(yorig) > 0.0 .OR. ABS(zorig) > 0.0) THEN
            WRITE(*,*) 'readini_forward: Shifting geophone positions'
         ENDIF
         xg(:) = xg(:) + xorig
         yg(:) = yg(:) + yorig
         zg(:) = zg(:) + zorig
      ENDIF
      igreg = INIPARSER_GETINT('receivers:igreg'//CHAR(0), 0, ierr)
      gspread = INIPARSER_GETFLOAT('receivers:gspread', 0.5, ierr)
      usegwt = INIPARSER_GETBOOL('receivers:usegwt', false, ierr)
      IF (ng > 0 .AND. &
          ((gspread <= 0.5 .AND. 2*igreg + 1 > nrgg) .OR. &
           (gspread >  0.5 .AND. INT(4.0*gspread) + 1 > nrgg))) THEN
         WRITE(*,*) 'nrgg too small for geophones'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.usegwt .AND. ng > 0) gwght(:) = 1.0
      ! read the receivers - these are pressure
      nr = INIPARSER_GETINT('receivers:nr'//CHAR(0), 0, ierr) 
      IF (nr > 0) THEN
         ALLOCATE(rwght(nr))
         ALLOCATE(xr(nr))
         ALLOCATE(yr(nr))
         ALLOCATE(zr(nr))
         DO ir=1,nr
            cvar(:) = ' ' 
            var(:) = ' ' 
            ck(:) = ' ' 
            WRITE(ck, '(I6)') ir
            var = 'receivers:rec_'//TRIM(ADJUSTL(ck))//CHAR(0)
            ierr = iniparser_getRecvWXYZ(var, rwght(ir), xr(ir), yr(ir), zr(ir))
            IF (ierr /= 0) THEN
               WRITE(*,*) 'readini_forward: Error reading receiver', ir
               RETURN
            ENDIF
         ENDDO
         IF (ABS(xorig) > 0.0 .OR. ABS(yorig) > 0.0 .OR. ABS(zorig) > 0.0) THEN
            WRITE(*,*) 'readini_forward: Shifting receiver positions'
         ENDIF
         xr(:) = xr(:) + xorig
         yr(:) = yr(:) + yorig
         zr(:) = zr(:) + zorig
      ENDIF
      irreg = INIPARSER_GETINT('receivers:irreg'//CHAR(0), igreg, ierr)
      recin = INIPARSER_GETINT('receivers:recin'//CHAR(0), 3, ierr)
      rspread = INIPARSER_GETFLOAT('receivers:rspread'//CHAR(0), gspread, ierr)
      userwt = INIPARSER_GETBOOL('receivers:userwt'//CHAR(0), false, ierr)
      IF (nr > 0 .AND. &
          ((rspread <= 0.5 .AND. 2*irreg + 1 > nrgg) .OR. &
           (rspread >  0.5 .AND. INT(4.0*rspread) + 1 > nrgg))) THEN
         WRITE(*,*) 'nrgg too small for receivers'
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.userwt .AND. nr > 0) rwght(:) = 1.0
      ! read the sources
      ns = INIPARSER_GETINT('sources:ns'//CHAR(0), 0, ierr)
      IF (ns < 1) THEN
         WRITE(*,*) 'readini_forward: Error no sources'
         ierr = 1
         RETURN
      ENDIF
      DO is=1,ns
         cvar(:) = ' ' 
         var(:) = ' ' 
         ck(:) = ' ' 
         WRITE(ck, '(I6)') is
         var = 'sources:src_'//TRIM(ADJUSTL(ck))//CHAR(0)

      ENDDO
      isreg = INIPARSER_GETINT('sources:irreg'//CHAR(0), 0, ierr)
      srcin = INIPARSER_GETINT('sources:srcin'//CHAR(0), 3, ierr)
      sspread = INIPARSER_GETFLOAT('sources:sspread'//CHAR(0), 0.5, ierr)
      IF (ns > 0 .AND. &
          ((sspread <= 0.5 .AND. 2*isreg + 1 > nrgg) .OR. &
           (sspread >  0.5 .AND. INT(4.0*sspread) + 1 > nrgg))) THEN
         WRITE(*,*) 'readini_forward: nrgg too small for sources' 
         ierr = 1
         RETURN
      ENDIF 
      nsg = INIPARSER_GETINT('sources:nsg'//CHAR(0), 1, ierr)
      IF (nsg /= 1) THEN
         WRITE(*,*) 'readini_forward: Error nsg = 1 only allowed'
         ierr = 1
         RETURN
      ENDIF
      ALLOCATE(isg(nsg+1))
      ! optional parameters
      tau = INIPARSER_GETFLOAT('optional:tau'//CHAR(0), 999.999, ierr)
      ncom = INIPARSER_GETINT('optional:comment'//CHAR(0), 8, ierr)
      usemin = INIPARSER_GETBOOL('optional:usemin'//CHAR(0), true, ierr)
      totfld = INIPARSER_GETBOOL('optional:totfld'//CHAR(0), true, ierr)
      infspace = INIPARSER_GETBOOL('optional:infspace'//CHAR(0), false, ierr) 
      ! miscellaneous inversion parameters
      nomi = INIPARSER_GETINT('inversion:nomi'//CHAR(0), nom, ierr)
      modsrcp = INIPARSER_GETINT('inversion:modsrcp'//CHAR(0), -1, ierr)
      modrecp = INIPARSER_GETINT('inversion:modrecp'//CHAR(0), -1, ierr)
      ! free the memory associated with iniparser
      CALL iniparser_finalize() 
      RETURN
      END 

