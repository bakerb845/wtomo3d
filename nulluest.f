!>    @brief Null the three-component frequency responses 
!>
!>    @param[in] nom     number of frequencies
!>    @param[in] ng      number of geophones
!>    @param[in] ns      number of sources
!>
!>    @param[out] uest   U component set to zero
!>    @param[out] vest   V component set to zero
!>    @param[out] west   W component set to zero
!>
!>    @author Ben Baker
!>
      SUBROUTINE NULLUEST(nom, ng, ns, uest, vest, west) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nom, ng, ns
      COMPLEX, DIMENSION(:,:,:), INTENT(OUT) :: uest, vest, west
      COMPLEX, PARAMETER :: zero = CMPLX(0.0, 0.0)
      INTEGER i, j, k
      DO 1 k=1,ns
         DO 2 j=1,ng
            DO 3 i=1,nom      
               uest(i,j,k) = zero
               vest(i,j,k) = zero
               west(i,j,k) = zero
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
!>    @brief Null the receiver pressure frequency responses
!>
!>    @param[in] nom     number of frequencies
!>    @param[in] ng      number of geophones
!>    @param[in] ns      number of sources
!>
!>    @param[out] utest  pressure response to zero
!>
!>    @author Ben Baker
!>
      SUBROUTINE NULLUTEST(nom, nr, ns, utest)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nom, nr, ns
      COMPLEX, DIMENSION(:,:,:), INTENT(OUT) :: utest
      COMPLEX, PARAMETER :: zero = CMPLX(0.0, 0.0)
      INTEGER i, j, k
      DO 1 k=1,ns
         DO 2 j=1,nr
            DO 3 i=1,nom
               utest(i,j,k) = zero
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN 
      END 

