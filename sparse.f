!>
!>    @brief Converts the George and Liu xadj/adjncy sparse matrix
!>           representation to the Compressed Row Storage format
!>
!>    @param[in] nzero    number of non-zeros in graph
!>    @param[in] n        number of rows in matrix
!>    @param[in] xadj     George and Liu row pointer [n+1]
!>    @param[in] adjncy   George and Liu column pointer [xadj(n+1)-1]
!>
!>    @param[out] irptr   CRS row pointer [n+1]
!>    @param[out] jcptr   CRS column pointer [irptr(n+1)-1]
!>
      SUBROUTINE SPARSE_ADJ2CRS(nzero, n, xadj, adjncy,
     ;                          irptr, jcptr)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nzero, n
      INTEGER, INTENT(IN) :: xadj(n+1), adjncy(nzero-n)
      INTEGER, INTENT(OUT) :: irptr(n+1), jcptr(nzero)
      ! local variables
      INTEGER i, izero, j, jbeg, jend
      LOGICAL ldiag
      !----------------------------------------------------------------!
      !
      ! loop on rows
      irptr(1) = 1
      izero = 0
      DO 1 i=1,n
         jbeg = xadj(i)
         jend = xadj(i+1) - 1
         ldiag = .FALSE.
         ! loop on colums
         DO 2 j=jbeg,jend
            IF (.NOT.ldiag) THEN
               IF (adjncy(j) > i) THEN
                  izero = izero + 1
                  jcptr(izero) = i
                  ldiag = .true.
               ENDIF
            ENDIF
            izero = izero + 1
            jcptr(izero) = adjncy(j)
    2    CONTINUE !loop on column
         irptr(i+1) = izero + 1
    1 CONTINUE
      izero = izero + 1
      jcptr(izero) = n
      irptr(n+1) = nzero + 1
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
!>    @brief Converts a compressed sparse row row to compressed ordinate
!>
!>    @param[in] nzero    number of non-zeros in matrix
!>    @param[in] n        number of rows in matrix
!>    @param[in] irptr    maps from i'th row to start index of jcptr
!>                        [n+1]
!>    @param[in] jcptr    maps from iz'th non-zero to column number
!>                        [nzero] 
!>
!>    @param[out] irn     maps from iz'th non-zero to row index [nzero]
!>    @param[out] jcn     maps from iz'th non-zero to column index
!>                        [nzero]
!>
!>    @author Ben Baker
!>
      SUBROUTINE SPARSE_CRS2COO(nzero, n, irptr, jcptr, irn, jcn)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nzero, n
      INTEGER, INTENT(IN) :: irptr(n+1), jcptr(nzero) 
      INTEGER, INTENT(OUT) :: irn(nzero), jcn(nzero)
      INTEGER i, j, jbeg, jend
      DO 1 i=1,n
         jbeg = irptr(i)
         jend = irptr(i+1) - 1
         DO 2 j=jbeg,jend
            irn(j) = i
            jcn(j) = jcptr(j)
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
