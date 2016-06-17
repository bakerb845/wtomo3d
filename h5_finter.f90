      MODULE HDF5_MODULE
         INTERFACE
            SUBROUTINE H5_READ_DOUBLE(dset_name, file_id, x, ierr)
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: dset_name
            INTEGER(C_INT), INTENT(IN) :: file_id
            REAL(C_DOUBLE), ALLOCATABLE, INTENT(OUT) :: x(:)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE H5_READ_DOUBLE

            SUBROUTINE H5_READ_FLOAT(dset_name, file_id, x, ierr)
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: dset_name
            INTEGER(C_INT), INTENT(IN) :: file_id
            REAL(C_FLOAT), ALLOCATABLE, INTENT(OUT) :: x(:)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE H5_READ_FLOAT

            SUBROUTINE H5_READ_INT(dset_name, file_id, x, ierr)
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: dset_name
            INTEGER(C_INT), INTENT(IN) :: file_id
            INTEGER(C_INT), ALLOCATABLE, INTENT(OUT) :: x(:)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE H5_READ_INT

            !============================================================================!
            !                            This is the C interface                         !
            !============================================================================!
            LOGICAL(C_BOOL) FUNCTION DIRECTORY_EXISTS(dirnm) &
                                     BIND(C,name='os_path_isdir')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: dirnm(*)
            END FUNCTION DIRECTORY_EXISTS

            LOGICAL(C_BOOL) FUNCTION FILE_EXISTS(filenm) &
                                     BIND(C,name='os_path_isfile')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: filenm(*)
            END FUNCTION FILE_EXISTS

!           INTEGER(C_INT) FUNCTION MAKE_DIRECTORY(dirnm) &
!                                   BIND(C,name='utilsFilesMakeDirectory')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           CHARACTER(KIND=C_CHAR), INTENT(IN) :: dirnm(*)
!           END FUNCTION MAKE_DIRECTORY

            INTEGER(C_INT) FUNCTION H5_READ_DOUBLE_ARRAY(dset_name, file_id, n, x) &
                                    BIND(C,name='h5_read_array__double')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: dset_name(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            INTEGER(C_INT), INTENT(IN) :: n
            REAL(C_DOUBLE), INTENT(OUT) :: x(n)
            END FUNCTION H5_READ_DOUBLE_ARRAY

            LOGICAL(C_BOOL) FUNCTION H5_ITEM_EXISTS(citem, file_id) &
                                     BIND(C,name='h5_item_exists')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: citem
            INTEGER(C_INT), INTENT(IN) :: file_id
            END FUNCTION H5_ITEM_EXISTS

            INTEGER(C_INT) FUNCTION H5_READ_FLOAT_ARRAY(dset_name, file_id, n, x) &
                                    BIND(C,name='h5_read_array__float')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: dset_name(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            INTEGER(C_INT), INTENT(IN) :: n
            REAL(C_FLOAT), INTENT(OUT) :: x(n)
            END FUNCTION H5_READ_FLOAT_ARRAY

!           INTEGER(C_INT) FUNCTION H5_READ_FLOAT_ELEMENTS(dset_name, file_id, ne, &
!                                                          locs, n, x) &
!                                       BIND(C,name='h5_read_float_elements')
!           USE ISO_C_BINDING
!           IMPLICIT NONE
!           CHARACTER(KIND=C_CHAR), INTENT(IN) :: dset_name(*)
!           INTEGER(C_INT), INTENT(IN) :: file_id
!           INTEGER(C_INT), INTENT(IN) :: locs(ne), ne, n
!           REAL(C_FLOAT), INTENT(OUT) :: x(n)
!           END FUNCTION H5_READ_FLOAT_ELEMENTS

            INTEGER(C_INT) FUNCTION H5_READ_INT_ARRAY(dset_name, file_id, n, x) &
                                    BIND(C,name='h5_read_array__int')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: dset_name(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            INTEGER(C_INT), INTENT(IN) :: n
            INTEGER(C_INT), INTENT(OUT) :: x(n)
            END FUNCTION H5_READ_INT_ARRAY

            INTEGER(C_INT) FUNCTION H5_GET_ARRAY_SIZE(citem, file_id) &
                                    BIND(C,name='h5_get_array_size')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: file_id
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: citem(*)
            END FUNCTION H5_GET_ARRAY_SIZE

            INTEGER(C_INT) FUNCTION H5_N_GROUP_MEMBERS(group_name, file_id) &
                                    BIND(C,name='h5_n_group_members')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: group_name(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            END FUNCTION H5_N_GROUP_MEMBERS

            INTEGER(C_INT) FUNCTION H5_OPEN_RDONLY(filenm) &
                                    BIND(C,name='h5_open_rdonly')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: filenm(*)
            END FUNCTION H5_OPEN_RDONLY

            INTEGER(C_INT) FUNCTION H5_OPEN_RDWT(filenm) &
                                    BIND(C,name='h5_open_rdwt')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: filenm(*)
            END FUNCTION H5_OPEN_RDWT

            INTEGER(C_INT) FUNCTION H5_OPEN_NEW(filenm) &
                                    BIND(C,name='h5_open_new')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: filenm(*)
            END FUNCTION H5_OPEN_NEW

            INTEGER(C_INT) FUNCTION H5_CLOSE(file_id) &
                                    BIND(C,name='h5_close')
            USE ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: file_id
            END FUNCTION H5_CLOSE

            INTEGER(C_INT) FUNCTION H5_WRITE_DOUBLE_ARRAY(dset_name, file_id, n, x) &
                                    BIND(C,name='h5_write_array__double')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: dset_name(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            INTEGER(C_INT), INTENT(IN) :: n
            REAL(C_DOUBLE), INTENT(IN) :: x(n)
            END FUNCTION H5_WRITE_DOUBLE_ARRAY

            INTEGER(C_INT) FUNCTION H5_WRITE_FLOAT_ARRAY(dset_name, file_id, n, x) &
                                    BIND(C,name='h5_write_array__float')
            USE ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: dset_name(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            INTEGER(C_INT), INTENT(IN) :: n
            REAL(C_FLOAT), INTENT(IN) :: x(n)
            END FUNCTION H5_WRITE_FLOAT_ARRAY

            INTEGER(C_INT) FUNCTION H5_WRITE_INT_ARRAY(dset_name, file_id, n, x) &
                                    BIND(C,name='h5_write_array__int')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(IN) :: dset_name(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            INTEGER(C_INT), INTENT(IN) :: n
            INTEGER(C_INT), INTENT(IN) :: x(n)
            END FUNCTION H5_WRITE_INT_ARRAY

            INTEGER(C_INT) FUNCTION H5_CREATE_GROUP(file_id, cgroup) &
                                    BIND(C,name='h5_create_group')
            USE ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR), INTENT(In) :: cgroup(*)
            INTEGER(C_INT), INTENT(IN) :: file_id
            END FUNCTION H5_CREATE_GROUP

         END INTERFACE
      END MODULE HDF5_MODULE
!>    Reads a double array from an HDF5 file
!>
!>    @param[in] dset_name     name of dataset to read (null terminated)
!>    @param[in] file_id       HDF5 file handle
!>
!>    @param[out] x            array x from HDF5 file
!>    @param[out] ierr         0 indicates success
!> 
!>    @author Ben Baker, ISTI
!>
      SUBROUTINE H5_READ_DOUBLE(dset_name, file_id, x, ierr)
      USE ISO_C_BINDING
      USE HDF5_MODULE, ONLY : h5_get_array_size, h5_read_double_array
      CHARACTER(*), INTENT(IN) :: dset_name
      INTEGER(C_INT), INTENT(IN) :: file_id
      REAL(C_DOUBLE), ALLOCATABLE, INTENT(OUT) :: x(:)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! local variables
      INTEGER(C_INT) n
      n = h5_get_array_size(dset_name, file_id)
      IF (n < 1) THEN
         WRITE(*,*) 'h5_read_double: No dataset:',TRIM(ADJUSTL(dset_name))
         ierr = 1
         RETURN
      ENDIF
      IF (ALLOCATED(x)) DEALLOCATE(x)
      ALLOCATE(x(n))
      ierr = h5_read_double_array(dset_name, file_id, n, x)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'h5_read_double: Error reading:',TRIM(ADJUSTL(dset_name))
         ierr = 1
         RETURN
      ENDIF
      RETURN 
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    Reads a float array from an HDF5 file
!>
!>    @param[in] dset_name     name of dataset to read (null terminated)
!>    @param[in] file_id       HDF5 file handle
!>
!>    @param[out] x            array x from HDF5 file
!>    @param[out] ierr         0 indicates success
!> 
!>    @author Ben Baker, ISTI
!>
      SUBROUTINE H5_READ_FLOAT(dset_name, file_id, x, ierr)
      USE ISO_C_BINDING
      USE HDF5_MODULE, ONLY : h5_get_array_size, h5_read_float_array
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: dset_name
      INTEGER(C_INT), INTENT(IN) :: file_id
      REAL(C_FLOAT), ALLOCATABLE, INTENT(OUT) :: x(:)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! local variables
      INTEGER(C_INT) n
      n = h5_get_array_size(dset_name, file_id)
      IF (n < 1) THEN
         WRITE(*,*) 'h5_read_float: No dataset:',TRIM(ADJUSTL(dset_name))
         ierr = 1 
         RETURN
      ENDIF
      IF (ALLOCATED(x)) DEALLOCATE(x)
      ALLOCATE(x(n))
      ierr = h5_read_float_array(dset_name, file_id, n, x)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'h5_read_float: Error reading:',TRIM(ADJUSTL(dset_name))
         ierr = 1 
         RETURN
      ENDIF
      RETURN 
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    Reads an integer array from an HDF5 file
!>
!>    @param[in] dset_name     name of dataset to read (null terminated)
!>    @param[in] file_id       HDF5 file handle
!>
!>    @param[out] x            array x from HDF5 file
!>    @param[out] ierr         0 indicates success
!> 
!>    @author Ben Baker, ISTI
!>
      SUBROUTINE H5_READ_INT(dset_name, file_id, x, ierr)
      USE ISO_C_BINDING
      USE HDF5_MODULE, ONLY : h5_get_array_size, h5_read_int_array
      CHARACTER(*), INTENT(IN) :: dset_name
      INTEGER(C_INT), INTENT(IN) :: file_id
      INTEGER(C_INT), ALLOCATABLE, INTENT(OUT) :: x(:)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! local variables
      INTEGER(C_INT) n
      n = h5_get_array_size(dset_name, file_id)
      IF (n < 1) THEN
         WRITE(*,*) 'h5_read_int: No dataset:',TRIM(ADJUSTL(dset_name))
         ierr = 1 
         RETURN
      ENDIF
      IF (ALLOCATED(x)) DEALLOCATE(x)
      ALLOCATE(x(n))
      ierr = h5_read_int_array(dset_name, file_id, n, x)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'h5_read_int: Error reading:',TRIM(ADJUSTL(dset_name))
         ierr = 1
         RETURN
      ENDIF
      RETURN
      END SUBROUTINE
 
