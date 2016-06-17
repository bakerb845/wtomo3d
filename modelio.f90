      SUBROUTINE MODELIO_READ_MODEL(projnm, ierr)
      USE MODEL_MODULE, ONLY : nx, ny, nz, vpr, vsr, rho, qp, qs, &
                               da, mu, qpex, qsex
      USE HDF5_MODULE, ONLY : file_exists, h5_close, h5_item_exists, &
                              h5_open_rdonly, h5_read_float
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTERFACE
         SUBROUTINE MODELIO_H52WT_MODEL(nx, ny, nz, h5mod, wmod)
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN) :: h5mod
         INTEGER, INTENT(IN) :: nx, ny, nz
         REAL, DIMENSION(:,:,:), INTENT(OUT) :: wmod
         END SUBROUTINE MODELIO_H52WT_MODEL
      END INTERFACE
      CHARACTER(*), INTENT(IN) :: projnm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(256) h5fl
      INTEGER(C_INT) file_id
      REAL(C_FLOAT), ALLOCATABLE :: x(:)
      COMPLEX vc
      INTEGER ix, iy, iz, nxyz
      ! open h5 model file
      ierr = 0
      h5fl(:) = ' '
      h5fl = TRIM(ADJUSTL(projnm))//'_model.h5'//CHAR(0)
      IF (.NOT. FILE_EXISTS(h5fl)) THEN
         WRITE(*,*) 'modelio_read_model: Error file doesnt exist ', TRIM(ADJUSTL(h5fl))
         ierr = 1
         RETURN
      ENDIF
      file_id = H5_OPEN_RDONLY(h5fl)
      ! set space
      IF (.NOT.ALLOCATED(vpr)) ALLOCATE(vpr(nz, nx, ny), STAT=ierr)
      IF (.NOT.ALLOCATED(vsr)) ALLOCATE(vsr(nz, nx, ny), STAT=ierr)
      IF (.NOT.ALLOCATED(rho)) ALLOCATE(rho(nz, nx, ny), STAT=ierr)
      IF (.NOT.ALLOCATED(qp))  ALLOCATE(qp(nz, nx, ny), STAT=ierr)
      IF (.NOT.ALLOCATED(qs))  ALLOCATE(qs(nz, nx, ny), STAT=ierr)
      IF (.NOT.ALLOCATED(da))  ALLOCATE(da(nz, nx, ny), STAT=ierr)
      IF (.NOT.ALLOCATED(mu))  ALLOCATE(mu(nz, nx, ny), STAT=ierr)
      ! get vp model 
      CALL H5_READ_FLOAT('/Model/Vp'//CHAR(0), file_id, x, ierr)
      nxyz = SIZE(x)
      IF (nxyz /= nx*ny*nz .OR. ierr /= 0) THEN
         WRITE(*,*) 'model_read_model: Error inconsistent sizes', nx, ny, nz, nxyz
         ierr = 1
         RETURN
      ENDIF
      CALL MODELIO_H52WT_MODEL(nx, ny, nz, x, vpr)
      IF (ALLOCATED(x)) DEALLOCATE(x)
      ! get the vs model
      CALL H5_READ_FLOAT('/Model/Vs'//CHAR(0), file_id, x, ierr)
      CALL MODELIO_H52WT_MODEL(nx, ny, nz, x, vsr)
      IF (ALLOCATED(x)) DEALLOCATE(x)
      ! get the density model
      CALL H5_READ_FLOAT('/Model/Density'//CHAR(0), file_id, x, ierr)
      CALL MODELIO_H52WT_MODEL(nx, ny, nz, x, rho)
      IF (ALLOCATED(x)) DEALLOCATE(x)
      ! check if there is a qp and qs model
      qpex = H5_ITEM_EXISTS('/Model/Qp'//CHAR(0), file_id)
      qsex = H5_ITEM_EXISTS('/Model/Qs'//CHAR(0), file_id)
      qp(:,:,:) = 0.0
      qs(:,:,:) = 0.0
      IF (qpex) THEN
         WRITE(*,*) 'modelio_read_model: make function to read qp'
         ierr = 1
         RETURN
      ENDIF
      IF (qsex) THEN
         WRITE(*,*) 'modelio_read_model: make function to read qs'
         ierr = 1
         RETURN
      ENDIF
      ! compute lambda and mu
      DO 101 iy=1,ny
         DO 102 ix=1,nx
            DO 103 iz=1,nz
               vc = CMPLX(vsr(iz,ix,iy), 0.0) !TODO fix here with p5 mkmod3d.f
               mu(iz,ix,iy) = vc*vc*rho(iz,ix,iy)
               vc = CMPLX(vpr(iz,ix,iy), 0.0) !TODO fix here with p2 mkmod3d.f
               da(iz,ix,iy) = vc*vc*rho(iz,ix,iy) - 2.0*mu(iz,ix,iy)
               !vc = CMPLX(vsr(iz,ix,iy), 0.0)
               !mup(iz,ix,iy) = vc*vc*rhop(iz,ix,iy)
               !vc = CMPLX(vpr(iz,ix,iy) 0.0)
               !dap(iz,ix,iy) = vc*vc*rhop(iz,ix,iy) - 20*mup(iz,ix,iy)
  103       CONTINUE
  102    CONTINUE
  101 CONTINUE 
      ! close the model file
      ierr = H5_CLOSE(file_id)
      RETURN
      END

      SUBROUTINE MODELIO_H52WT_MODEL(nx, ny, nz, h5mod, wmod)
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: h5mod
      INTEGER, INTENT(IN) :: nx, ny, nz
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: wmod
      INTEGER indx, ix, iy, iz 
      DO 1 iy=1,ny
         DO 2 ix=1,nx
            DO 3 iz=1,nz
               indx = (iz - 1)*nx*ny + (iy - 1)*nx + ix
               wmod(iz,ix,iy) = h5mod(indx)
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN
      END

      SUBROUTINE MODELIO_WT2H5_MODEL(nx, ny, nz, wmod, h5mod)
      IMPLICIT NONE
      REAL, DIMENSION(:,:,:), INTENT(IN) :: wmod
      INTEGER, INTENT(IN) :: nx, ny, nz
      REAL, DIMENSION(:), INTENT(OUT) :: h5mod
      INTEGER indx, ix, iy, iz
      DO 1 iz=1,nz
         DO 2 iy=1,ny
            DO 3 ix=1,nx
               indx = (iz - 1)*nx*ny + (iy - 1)*nx + ix
               h5mod(indx) = wmod(iz,ix,iy)
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
