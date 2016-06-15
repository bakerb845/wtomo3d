c======================================================================c
c                                                                      c
      subroutine qbuf3d(ld1, ld2, nvars, job, l, m, n, qbuf, amat)
c     puts or/extracts a 3d complex matrix to a buffer 
      INTEGER, INTENT(IN) :: l, ld1, ld2, m, n, nvars, job
      COMPLEX, INTENT(INOUT) :: amat(ld1,ld2,*), qbuf(nvars)
      INTEGER i, j, k, iter 
      iter = 0
      if (job.eq.0) then !matrix -> buffer
         do k=1,n
            do j=1,m
               do i=1,l
                  iter = iter + 1
                  qbuf(iter) = amat(i,j,k)
               enddo
            enddo
         enddo
      else !buffer -> matrix
         do k=1,n
            do j=1,m
               do i=1,l
                  iter = iter + 1
                  amat(i,j,k) = qbuf(iter)
               enddo
            enddo
         enddo
      endif
      return
      end

c======================================================================c

