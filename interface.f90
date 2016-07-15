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
            INTEGER, INTENT(IN) ::      lbl, iz, ix, iy
            INTEGER, INTENT(OUT) ::     ierr
            END SUBROUTINE fdfd3d

            SUBROUTINE fixed3d(lbl)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: lbl 
            END SUBROUTINE fixed3d

            SUBROUTINE hetfldsrc (srcfn, sfld, myid, master, ierr)
            IMPLICIT NONE
            COMPLEX, INTENT(IN) :: srcfn
            COMPLEX, DIMENSION(:), INTENT(OUT) :: sfld
            INTEGER, INTENT(IN) :: myid, master
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE hetfldsrc

            SUBROUTINE INFO(nx, ny, nz, ns, nr, ng, nom)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nx, ny, nz, ns, nr, ng, nom 
            END SUBROUTINE INFO

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
            CHARACTER(C_CHAR), INTENT(IN) :: variable(*)
            REAL(C_FLOAT), INTENT(IN) :: def 
            INTEGER(C_INT), INTENT(OUT) :: ierr  
            END FUNCTION iniparser_getFloat

            INTEGER(C_INT) FUNCTION iniparser_getInt(variable, def, ierr) &
                           BIND(C, NAME='iniparser_getInt')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: variable(*)
            INTEGER(C_INT), INTENT(IN) :: def 
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END FUNCTION iniparser_getInt

            INTEGER(C_INT) FUNCTION iniparser_getChar(variable, def, ret) &
                           BIND(C, NAME='iniparser_getChar')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: variable(*), def(*)
            CHARACTER(C_CHAR), INTENT(OUT) :: ret(*)
            END FUNCTION iniparser_getChar

            LOGICAL(C_BOOL) FUNCTION iniparser_getBool(variable, def, ierr) &
                            BIND(C, NAME='iniparser_getBool')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: variable(*)
            LOGICAL(C_BOOL), INTENT(IN) :: def
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END FUNCTION iniparser_getBool

            INTEGER(C_INT) FUNCTION iniparser_getRecvWXYZ(variable,  &
                                                          wt, x, y, z) &
                           BIND(C, NAME='iniparser_getRecvWXYZ')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: variable
            REAL(C_FLOAT), INTENT(OUT) :: wt, x, y, z
            END FUNCTION iniparser_getRecvWXYZ

            INTEGER(C_INT) FUNCTION iniparser_init(ini_file) &
                           BIND(C, NAME='iniparser_init')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(C_CHAR), INTENT(IN) :: ini_file(*)
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

            SUBROUTINE NULLUEST(nom, ng, ns, uest, vest, west)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nom, ng, ns
            COMPLEX, DIMENSION(:,:,:), INTENT(OUT) :: uest, vest, west
            END SUBROUTINE NULLUEST

            SUBROUTINE NULLUTEST(nom, nr, ns, utest)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nom, nr, ns
            COMPLEX, DIMENSION(:,:,:), INTENT(OUT) :: utest
            END SUBROUTINE NULLUTEST

            SUBROUTINE pfdfd3d (lbl, ierr)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: lbl
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE

            SUBROUTINE pml3d (lbl, iz, ix, iy, ierr)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: lbl, ix, iy, iz
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
            IMPLICIT NONE
            CHARACTER(256), INTENT(IN) :: projnm
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE READINI_FORWARD

            SUBROUTINE RECSUB(projnm,               &
                              ncom, nsam, nom, ngt, &
                              recin,                &
                              deltatt, freq,        &
                              receiver, ierr)

            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: projnm 
            REAL, DIMENSION(:), INTENT(IN) :: freq
            REAL, INTENT(IN) :: deltatt
            INTEGER, INTENT(IN) :: nom, ngt, nsam, recin, ncom
            COMPLEX, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: receiver
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE RECSUB

            SUBROUTINE SRCSUB(projnm, &
                              ncom, nsam, nom, nomi, &
                              nst, modsrcp, &
                              srcin, nsg, isg, &
                              deltatt, tau, freq, &
                              dsd, source, ierr)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: projnm
            REAL, DIMENSION(:), INTENT(IN) :: freq
            INTEGER, DIMENSION(:), INTENT(INOUT) :: isg
            REAL deltatt, tau
            INTEGER, INTENT(IN) :: nsg, nom, nomi, nst, nsam, &
                                   srcin, modsrcp, ncom
            COMPLEX, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: source
            REAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: dsd
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE SRCSUB


            SUBROUTINE TEST1D(is3d)
            IMPLICIT NONE
            LOGICAL, INTENT(OUT) :: is3d
            END SUBROUTINE TEST1D

            SUBROUTINE TPFREE(lbl, ix, iy, ierr)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: lbl, ix, iy
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE TPFREE

            SUBROUTINE ZERO()
            IMPLICIT NONE
            END SUBROUTINE ZERO

            SUBROUTINE ZEROAS()
            IMPLICIT NONE
            END SUBROUTINE ZEROAS

         END INTERFACE

      END MODULE INTERFACE_MODULE

      MODULE SORT_MODULE
         INTERFACE !SORT_INTER
            SUBROUTINE SORT_INT_FINTER(isort, nn) &
                       BIND(C,NAME='sort_int_finter')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(INOUT) :: isort(nn)
            INTEGER(C_INT), INTENT(IN) :: nn
            END SUBROUTINE SORT_INT_FINTER

            SUBROUTINE ARGSORT_INT_FINTER(iperm, isort, nn) &
                       BIND(C,NAME='argsort_int_finter')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: isort(nn), nn
            INTEGER(C_INT), INTENT(OUT) :: iperm(nn)
            END SUBROUTINE ARGSORT_INT_FINTER

            INTEGER(C_INT) FUNCTION BSEARCH_INT_FINTER(keyIn, values, n, iverb) &
                                        BIND(C,NAME='bsearch_int_finter')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: values(n), keyIn, n, iverb
            END FUNCTION BSEARCH_INT_FINTER
         END INTERFACE
      END MODULE SORT_MODULE

      MODULE SPARSE_MODULE
         INTERFACE !SPARSE_INTER
            SUBROUTINE SPARSE_ADJ2CRS(nzero, n, xadj, adjncy, &
                                      irptr, jcptr)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nzero, n
            INTEGER, INTENT(IN) :: xadj(n+1), adjncy(nzero-n)
            INTEGER, INTENT(OUT) :: irptr(n+1), jcptr(nzero)
            END SUBROUTINE SPARSE_ADJ2CRS
         END INTERFACE
      END MODULE SPARSE_MODULE
