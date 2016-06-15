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
