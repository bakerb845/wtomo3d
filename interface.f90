      MODULE INTERFACE_MODULE
         INTERFACE

             SUBROUTINE asmble(lunlog, nzero, nx, ny, nz, nzeror, &
                               irn, jcn, a, ierr)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: nzero, lunlog 
             INTEGER, INTENT(IN) :: nzeror(3*nx*ny*nz), irn(nzero), &
                                    jcn(nzero), nx, ny, nz
             COMPLEX, INTENT(OUT) :: a(nzero)
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE asmble

             SUBROUTINE calcwts(vp, vs, w1, w2, w3, wm1, wm2, wm3, wm4)
             IMPLICIT NONE
             REAL, INTENT(IN) :: vp, vs
             REAL, INTENT(OUT) :: w1, w2, w3, wm1, wm2, wm3, wm4
             END SUBROUTINE calcwts

             SUBROUTINE dispersion()
             IMPLICIT NONE
             END SUBROUTINE dispersion

             SUBROUTINE fdfd3d(lbl, iz, ix, iy, ierr)
             IMPLICIT NONE
             INTEGER, INTENT(IN) ::      iz, ix, iy
             CHARACTER(3), INTENT(IN) :: lbl 
             INTEGER, INTENT(OUT) ::     ierr
             END SUBROUTINE fdfd3d

             SUBROUTINE fixed3d(lbl)
             IMPLICIT NONE
             CHARACTER(3), INTENT(IN) :: lbl 
             END SUBROUTINE fixed3d

             SUBROUTINE hetfldsrc (srcfn, sfld, myid, master, ierr)
             IMPLICIT NONE
             COMPLEX, INTENT(IN) :: srcfn
             COMPLEX, DIMENSION(:), INTENT(OUT) :: sfld
             INTEGER, INTENT(IN) :: myid, master
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE hetfldsrc

             REAL(C_DOUBLE) FUNCTION iniparser_getDouble(variable, def, ierr) &
                            BIND(C, NAME='iniparser_getDouble')
             USE ISO_C_BINDING
             IMPLICIT NONE
             CHARACTER(C_CHAR), INTENT(IN) :: variable(256)
             REAL(C_DOUBLE), INTENT(IN) :: def
             INTEGER(C_INT), INTENT(OUT) :: ierr             
             END FUNCTION  iniparser_getDouble

             REAL(C_FLOAT) FUNCTION iniparser_getFloat(variable, def, ierr) &
                           BIND(C, NAME='iniparser_getFloat')
             USE ISO_C_BINDING
             IMPLICIT NONE
             CHARACTER(C_CHAR), INTENT(IN) :: variable(256)
             REAL(C_FLOAT), INTENT(IN) :: def 
             INTEGER(C_INT), INTENT(OUT) :: ierr  
             END FUNCTION iniparser_getFloat

             INTEGER(C_INT) FUNCTION iniparser_getInt(variable, def, ierr) &
                            BIND(C, NAME='iniparser_getInt')
             USE ISO_C_BINDING
             IMPLICIT NONE
             CHARACTER(C_CHAR), INTENT(IN) :: variable(256)
             INTEGER(C_INT), INTENT(IN) :: def 
             INTEGER(C_INT), INTENT(OUT) :: ierr
             END FUNCTION iniparser_getInt

             INTEGER(C_INT) FUNCTION iniparser_getChar(variable, def, ret) &
                            BIND(C, NAME='iniparser_getChar')
             USE ISO_C_BINDING
             IMPLICIT NONE
             CHARACTER(C_CHAR), INTENT(IN) :: variable(256), def(256)
             CHARACTER(C_CHAR), INTENT(OUT) :: ret(256)
             END FUNCTION iniparser_getChar

             INTEGER(C_INT) FUNCTION iniparser_init(ini_file) &
                            BIND(C, NAME='iniparser_init')
             USE ISO_C_BINDING
             IMPLICIT NONE
             CHARACTER(C_CHAR), INTENT(IN) :: ini_file(256)
             END FUNCTION iniparser_init

             SUBROUTINE iniparser_finalize() &
                        BIND(C, NAME='iniparser_finalize')
             END SUBROUTINE iniparser_finalize
 
             SUBROUTINE makeuh(isrc, sfld, ierr)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: isrc
             COMPLEX, DIMENSION(:), INTENT(IN) :: sfld
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE makeuh

             SUBROUTINE pfdfd3d (lbl, ierr)
             IMPLICIT NONE
             CHARACTER(3), INTENT(IN) :: lbl
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE

             SUBROUTINE pml3d (lbl, iz, ix, iy, ierr)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: ix, iy, iz
             CHARACTER(3), INTENT(IN) :: lbl
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE pml3d

             SUBROUTINE setcoefs3d(iz, ix, iy, ibg)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: iz, ix, iy, ibg 
             END SUBROUTINE setcoefs3d

             SUBROUTINE srcreg(fs, ireg, sspread, dx, dy, dz,  &
                               xoff, yoff, zoff, ierr)
             IMPLICIT NONE
             REAL, INTENT(OUT) :: fs(9, 9, 9)
             REAL, INTENT(IN) :: sspread, dx, dy, dz, xoff, yoff, zoff
             INTEGER, INTENT(IN) :: ireg
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE srcreg

             SUBROUTINE srcsetup(lsetup,                       &
                                 xloc, yloc, zloc, sprd, ireg, &
                                 sterm, swave, ierr)
             LOGICAL, INTENT(IN) ::  lsetup
             REAL, INTENT(IN) :: xloc, yloc, zloc, sprd
             INTEGER, INTENT(IN) :: ireg
             COMPLEX, DIMENSION(:,:,:), INTENT(INOUT) :: swave
             COMPLEX, INTENT(INOUT) :: sterm
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE srcsetup

             SUBROUTINE srel3d_getsize(nx, ny, nz, numvert, nedge)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: nx, ny, nz 
             INTEGER, INTENT(OUT) :: numvert, nedge
             END SUBROUTINE srel3d_getsize

             SUBROUTINE srel3d(nx, ny, nz, numvert, nedge, xadj, adjncy)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: nedge, numvert, nx, ny, nz
             INTEGER, INTENT(OUT) :: xadj(numvert+1), adjncy(2*nedge)
             END SUBROUTINE srel3d

             SUBROUTINE READINI_FORWARD(projnm, ierr)
             USE ISO_C_BINDING
             IMPLICIT NONE
             CHARACTER(256), INTENT(IN) :: projnm
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE READINI_FORWARD

             SUBROUTINE test1d(is3d)
             IMPLICIT NONE
             LOGICAL, INTENT(OUT) :: is3d
             END SUBROUTINE test1d

             SUBROUTINE tpfree(lbl, ix, iy, ierr)
             IMPLICIT NONE
             INTEGER, INTENT(IN) :: ix, iy
             CHARACTER(3), INTENT(IN) :: lbl 
             INTEGER, INTENT(OUT) :: ierr
             END SUBROUTINE tpfree

             SUBROUTINE zero()
             IMPLICIT NONE
             END SUBROUTINE zero

             SUBROUTINE zeroas()
             IMPLICIT NONE
             END SUBROUTINE zeroas

         END INTERFACE

      END MODULE INTERFACE_MODULE
