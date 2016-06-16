      SUBROUTINE READINI_FORWARD(projnm, ierr)
      USE ISO_C_BINDING
      USE MODEL_MODULE, ONLY : xorig, yorig, zorig, nx, ny, nz
      USE INTERFACE_MODULE, ONLY : iniparser_getchar, iniparser_getfloat, &
                                   iniparser_getint, iniparser_init, iniparser_finalize
      INTERFACE
         SUBROUTINE SET_C_STRING(variable, cstring)
         USE ISO_C_BINDING
         CHARACTER(*), INTENT(IN) :: variable
         CHARACTER(C_CHAR), INTENT(OUT) :: cstring(256)
         END SUBROUTINE SET_C_STRING
      END INTERFACE
      CHARACTER(256), INTENT(IN) :: projnm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(C_CHAR) inifl(256), var(256)
      ierr = 0
      ! set the ini file name (notice iniparser is C program and must be null terminated)
      CALL SET_C_STRING(TRIM(ADJUSTL(projnm))//'.ini', inifl)
!     DO i=1,256
!        inifl(i) = cwork(i:i)
!     ENDDO
      ierr = INIPARSER_INIT(inifl) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'readini_forward: Error reading ini file: ', inifl
         RETURN
      ENDIF
      ! get the pertinent variables for forward modeling
      CALL SET_C_STRING('model:nx', var)
      nx = INIPARSER_GETINT(var, 0, ierr)
      CALL SET_C_STRING('model:ny', var) 
      ny = INIPARSER_GETINT(var, 0, ierr)
      CALL SET_C_STRING('model:nz', var)
      nz = INIPARSER_GETINT(var, 0, ierr)

      CALL SET_C_STRING('model:xorig', var)
      xorig = INIPARSER_GETFLOAT(var, 0.0, ierr)
      CALL SET_C_STRING('model:yorig', var)
      yorig = INIPARSER_GETFLOAT(var, 0.0, ierr)
      CALL SET_C_STRING('model:zorig', var)
      zorig = INIPARSER_GETFLOAT(var, 0.0, ierr) 

      CALL SET_C_STRING('model:dx', var)
      dx = INIPARSER_GETFLOAT(var, 0.0, ierr)
      CALL SET_C_STRING('model:dy', var)
      dy = INIPARSER_GETFLOAT(var, 0.0, ierr)
      CALL SET_C_STRING('model:dz', var)
      dz = INIPARSER_GETFLOAT(var, 0.0, ierr)
 print *, nx, ny, nz, xorig, yorig, zorig, dx, dy, dz
      ! free the memory associated with iniparser
      CALL iniparser_finalize() 
      RETURN
      END 

      SUBROUTINE SET_C_STRING(variable, cstring)
      USE ISO_C_BINDING
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: variable
      CHARACTER(C_CHAR), INTENT(OUT) :: cstring(256)
      INTEGER i
      cstring(1:256) = C_NULL_CHAR 
      DO i=1,LEN_TRIM(variable)
         cstring(i) = variable(i:i)
      ENDDO
      RETURN
      END
