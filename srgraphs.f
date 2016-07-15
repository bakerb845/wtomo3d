c 
c     program srgraph: 
c     this routine makes a graph of the finite difference matrix for 
c     input to the metis package to compute an optimal reordering.
c
c     nx, nz number of nodes in the x and z directions
c
c     ordering of nodes (or "vertices" in graph speak) is z fast and 
c     x slow, so that the vertex corresponding to the node at (iz, ix) 
c     is iz + nz*(ix-1)
c
      SUBROUTINE srac_getsize(nx, nz, numvert, nedge)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, nz
      INTEGER, INTENT(OUT) :: numvert, nedge
      numvert = nx*nz
      nedge = (nz-1)*(5 + 4*(nx-2)) + (nx - 1)
      RETURN
      END

      SUBROUTINE srac(nx,nz, numvert,nedge, xadj,adjncy)
c 
c     acoustic graph  
c 
c     input      meaning 
c     -----      ------- 
c     nedge      number of edges 
c     numvert    number of vertices  
c     nx         number of x grid points 
c     nz         number of z grid points 
c   
c     output     meaning 
c     ------     ------- 
c     adjncy     edge vector 
c     xadj       verticy vector  
c
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nedge, numvert, nx, nz
      INTEGER, INTENT(OUT) :: xadj(numvert+1), adjncy(2*nedge)
      INTEGER iedge(8), i, iter, ix, iz, j, k, nedger, nvt, nvto
c 
c----------------------------------------------------------------------c
c
c-----compute the total number of independent edges
c     numvert = nx*nz
c     nedge = (nz-1)*(5 + 4*(nx-2)) + (nx - 1)
c-----loop over all vertices and collect the edges for each one
      iter = 0 
      nvto = 1 
      nedger = 0  
      do i = 1, nx
         do j = 1, nz
            nvt = 0
c-----------left side
            if (i.gt.1) then
               ix = i - 1
               do k = 1, 3
                  iz = j + k - 2
                  if (iz.gt.0.and.iz.le.nz) then
                     nvt = nvt + 1
                     iedge(nvt) = iz + nz*(ix-1)
                  endif
               enddo
            endif
c-----------middle
            ix = i
            do k = 1, 3
               iz = j + k - 2
               if (k.ne.2) then
                  if (iz.gt.0.and.iz.le.nz) then
                     nvt = nvt + 1
                     iedge(nvt) = iz + nz*(ix-1)
                  endif
               endif
            enddo
c-----------right side
            if (i.lt.nx) then
               ix = i + 1
               do k = 1, 3
                  iz = j + k - 2
                  if (iz.gt.0.and.iz.le.nz) then
                     nvt = nvt + 1
                     iedge(nvt) = iz + nz*(ix-1)
                  endif
               enddo
            endif
            iter = iter + 1
            xadj(iter) = nvto 
            nvto = nvto + nvt
            do k=1,nvt
               nedger = nedger + 1 
               adjncy(nedger) = iedge(k)
            enddo
         enddo !end loop on z
      enddo !end loop on x
      iter = iter + 1 
      xadj(iter) = nvto 
      return  
      end
c 
c======================================================================c
c 
c     program srgraph_elastic:
c     this routine makes a graph of the finite difference matrix for 
c     input to the metis package to compute an optimal reordering.
c
c     nx, nz number of nodes in the x and z directions
c
c     ordering of nodes (or "vertices" in graph speak) is z fast and x 
c     slow, so that the vertex corresponding to the node at (iz, ix) is 
c     iz + nz*(ix-1)
c
c     the metis graph file consists of a line for each diagonal term 
c     that lists the positions of all the off diagonal elements for that
c     row.  the first line in the file has the total number of diagonals
c     and the number of independant edges, by whic is meant the number
c     of nonrepeating combinations of diagonal and offdiagonals.   
c
c     for the acoustic case, the off diagonals are listed in 
c     ad, aa, af, dd, ff, cd, cc, cf order for each element.  the number
c     of diagonals is nx*nz.   to determine the number of unique 
c     combinations of diags and offdiags, consider that each element
c     will uniquely be associated with the 4 offdiagonals to the left 
c     of the diagonal (ff, cd, cc, and cf).  in the first row of the 
c     model, there will be 3 of these (ff, cc, cf), and the last row
c     will have two (cd and cc), and the last column only one (ff).  
c     to sum:
c
c       nz - 1 first row elements with 3 unique offdiags
c       nz - 1 last row elements with 2 unique offdiag
c       nx - 1 last column elements with 1 unique offdiag
c       (nx - 2)*(nz - 1) "interior" elements with 4 unique offdiags
c
c     summing:
c
c       3*(nz-1) + 2*(nz-1) + (nx-1) + (nz-1)*(nx-2)*4 
c     or
c       (nz-1)*(5 + 4*(nx-2)) + (nx-1)
c
c     the elastic case is analogous but we replace the individual 
c     elements by 2 x 2 matrices.  each node has 2 rows associated with
c     it, and each row has 17 off diagonal elements, which differ only 
c     switching the diagonal and its offdiagonal pair.  to compute the
c     number of unique associations, do the same trick as the acoustic 
c     case, but consider the coupling between rows.  a typical pair of 
c     (u,v) rows in a looks like this:
c
c     ...a1u a1v a2u a2v a3u a3v ... a4u a4v a5u a5v a6u a6v ... a7u a7v a8u a8v a9u a9v ...
c     ...a1u a1v a2u a2v a3u a3v ... a4u a4v a5u a5v a6u a6v ... a7u a7v a8u a8v a9u a9v ...
c
c     the diagonals are a5u and a5v.  if we do the acoustic trick for 
c     the first row, then the dependencies of a5u are 
c            a5v a6u a6v ... a7u a7v a8u a8v a9u a9v ...
c
c     for each node, there will be u->u terms, v->v terms, u->v terms, 
c     and v->u terms (i.e. 4 total per node) plus the reference to 
c     itself.  thus we should have
c
c       nx - 1 first row elements with 4*3 + 1 = 13 unique offdiags
c       nx - 1 last row elements with  4*2 + 1 =  9 unique offdiags
c       nz - 1 last col elements with  4*1 + 1 =  5 unique offdiags
c       (nz - 2)*(nx - 1) "interior" elements with 4*4 + 1 = 17 unique 
c       offdiags
c     the very last element with 1 offdiag
c
c     summing:
c       13*(nx-1) + 9*(nx-1) + 5*(nz-1) + 17*(nz-2)*(nx-1) + 1
c       (nx-1)*(22 + 17*(nz-2)) + 5*(nz-1) + 1
c     or
c       17*nx*nz -12(nx+nz) + 8
c
      SUBROUTINE srel_getsize(nx, nz, numvert, nedge)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, nz
      INTEGER, INTENT(OUT) :: numvert, nedge
      numvert = 2*nx*nz
      nedge = 17*nx*nz -12*(nx+nz) + 8
      RETURN
      END

      SUBROUTINE srel(nx, nz, numvert, nedge, xadj, adjncy)
c 
c     2d elastic graph 
c 
c     input      meaning 
c     -----      ------- 
c     nedge      number of edges 
c     numvert    number of vertices  
c     nx         number of x grid points 
c     nz         number of z grid points 
c   
c     output     meaning 
c     ------     ------- 
c     adjncy     edge vector 
c     xadj       verticy vector  c  
c
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nedge, numvert, nx, nz
      INTEGER, INTENT(OUT) :: xadj(numvert+1), adjncy(2*nedge)
      INTEGER iedge(18), i, iu, iv, ix, iz, j, iter,
     ;        k, n, nac, nedger, nvt, nvto
c 
c---------------------------------------------------------------------c
c
c     numvert = 2*nx*nz
c-----compute the total number of independant edges
c     nedge = 17*nx*nz -12*(nx+nz) + 8
c
c-----loop over all vertices and collect the edges for each one
c
c     note on indexing:
c
c     in the acoustic case, a given element located at row iz and column
c     ix in the model will have the following offdiags;
c       (iz-1, ix-1)  (iz-1, ix)  (iz-1, ix+1)
c       (iz  , ix-1)  (iz  , ix)  (iz  , ix+1)
c       (iz+1, ix-1)  (iz+1, ix)  (iz+1, ix+1)
c
c     and the absolute index is the z index + nz * (x index - 1)
c
c     for the elastic case, we have u and v parts for each model node.
c
c     the u element at node (iz, ix) will be 
c       2*[iz + nz*(ix-1)] -1
c     the v position is just 1+ that:
c       2*[iz + nz*(ix-1)]
c
      iter = 0 
      nvto = 1 
      nedger = 0 
      do i = 1, nx
         do j = 1, nz
            do n = 1, 2
               nvt = 0
c--------------left side
               if (i.gt.1) then
                  ix = i - 1
                  do k = 1, 3
                     iz = j + k - 2
                     if (iz.gt.0.and.iz.le.nz) then
                        nvt = nvt + 1
                        nac = nz*(ix-1) + iz
                        iu = 2*nac - 1    
                        iedge(nvt) = iu
                        nvt = nvt + 1
                        iv = iu + 1    
                        iedge(nvt) = iv
                     endif
                  enddo
               endif
c--------------middle
               ix = i
               do k = 1, 3
                  iz = j + k - 2
                  if (iz.gt.0.and.iz.le.nz) then
c--------------------condition to ignore the diagaonal for u
                     if (n.eq.1) then
                        if (k.ne.2) then
                           nvt = nvt + 1
                           nac = nz*(ix-1) + iz
                           iu = 2*nac - 1    
                           iedge(nvt) = iu 
                        endif
                        nvt = nvt + 1
                        nac = nz*(ix-1) + iz
                        iv = 2*nac
                        iedge(nvt) = iv
                     endif
c--------------------condition to ignore the diagaonal for v
                     if (n.eq.2) then
                        nvt = nvt + 1
                        nac = nz*(ix-1) + iz
                        iu = 2*nac - 1    
                        iedge(nvt) = iu 
                        if (k.ne.2) then
                           nvt = nvt + 1
                           nac = nz*(ix-1) + iz
                           iv = 2*nac
                           iedge(nvt) = iv
                        endif
                     endif
                  endif
               enddo
c--------------right side
               if (i.lt.nx) then
                  ix = i + 1
                  do k = 1, 3
                     iz = j + k - 2
                     if (iz.gt.0.and.iz.le.nz) then
                        nvt = nvt + 1
                        nac = iz + nz*(ix-1)    
                        iu = 2*nac - 1    
                        iedge(nvt) = iu
                        nvt = nvt + 1
                        iv = iu + 1    
                        iedge(nvt) = iv
                     endif
                  enddo
               endif
               iter = iter + 1
               xadj(iter) = nvto 
               nvto = nvto + nvt 
               do k=1,nvt
                  nedger = nedger +1 
                  adjncy(nedger) = iedge(k)
               enddo
c 200          format(17i10)
            enddo
         enddo !end loop on z 
      enddo !end loop on x 
      iter = iter + 1 
      xadj(iter) = nvto 
      return  
      end
c 
c=====================================================================c
c 
c     program srgraph_elasticpl:
c     this routine makes a graph of the finite difference matrix for 
c     input to the metis package to compute an optimal reordering.
c
c     this version for the 2.5d teleseismic case.  quick note: a good
c     way to think about ordering is that all the unique offdiags are
c     in the upper (or lower if you like) triangle of the symmetric 
c     s matrix.  for a 9 point star there will be generally be 4 
c     points in the upper triangle (less if you are at the edge).  if
c     there are n components of motion, then there will be 4*n unique
c     offidags.  each point will be represented by n rows, which also
c     means that each point in the star is represented by an nxn 
c     matix.  we need to specify the upper triangle of the middle 
c     matrix, which would b n*(n-1)/2 elements.  hence the total 
c     number for the general case will be:
c        4*n*n + n*(n-1)/2
c     again, this number decreases at points along the edge of the 
c     model, but computation is straightforward once you understand 
c     this concept. 
c
c     nb:  at this point i have only written comments - code has not 
c     yet been altered.
c       
c     nx, nz number of nodes in the x and z directions
c
c     ordering of nodes (or "vertices" in graph speak) is z fast and 
c     x slow, so that the vertex corresponding to the node at (iz, ix)
c     is iz + nz*(ix-1)
c
c     the metis graph file consists of a line for each diagonal term 
c     that lists the positions of all the off diagonal elements for 
c     that row.  the first line in the file has the total number of 
c     diagonals and the number of independant edges, by which is meant
c     the number of nonrepeating combinations of diagonal and 
c     offdiagonals.   
c
c     for the acoustic case, the off diagonals are listed in 
c     ad, aa, af, dd, ff, cd, cc, cf order for each element.  the 
c     number of diagonals is nx*nz.   to determine the number of 
c     unique combinations of diags and offdiags, consider that each 
c     element will uniquely be associated with the 4 offdiagonals to 
c     the left of the diagonal (ff, cd, cc, and cf).  
c 
c     in the first row of the model, there will be 3 of these 
c     (ff, cc, cf), and the last row will have two (cd and cc), and 
c     the last column only one (ff).  to sum:
c
c       nz - 1 first row elements with 3 unique offdiags
c       nz - 1 last row elements with 2 unique offdiag
c       nx - 1 last column elements with 1 unique offdiag
c       (nx - 2)*(nz - 1) "interior" elements with 4 unique offdiags
c
c     summing:
c
c       3*(nz-1) + 2*(nz-1) + (nx-1) + (nz-1)*(nx-2)*4 
c     or
c       (nz-1)*(5 + 4*(nx-2)) + (nx-1)
c
c     the elastic case is analogous but we replace the individual 
c     elements by 2 x 2 matrices.  each node has 2 rows associated 
c     with it, and each row has 17 off diagonal elements, which differ
c     only switching the diagonal and its offdiagonal pair.  to 
c     compute the number of unique associations, do the same trick as
c     the acoustic case, but consider the coupling between rows.  a 
c     typical pair of (u,v) rows in a looks like this:
c
c       ...a1u a1v a2u a2v a3u a3v ... a4u a4v a5u a5v a6u a6v ... a7u a7v a8u a8v a9u a9v ...
c       ...a1u a1v a2u a2v a3u a3v ... a4u a4v a5u a5v a6u a6v ... a7u a7v a8u a8v a9u a9v ...
c
c       the diagonals are a5u and a5v.  if we do the acoustic trick for
c       the first row, then the dependencies of a5u are 
c            a5v a6u a6v ... a7u a7v a8u a8v a9u a9v ...
c
c     for each node, there will be u->u terms, v->v terms, u->v terms, 
c     and v->u terms (i.e. 4 total per node) plus the reference to 
c     itself.  thus we should have
c
c       nx - 1 first row elements with 4*3 + 1 = 13 unique offdiags
c       nx - 1 last row elements with  4*2 + 1 =  9 unique offdiags
c       nz - 1 last col elements with  4*1 + 1 =  5 unique offdiags
c       (nz - 2)*(nx - 1) "interior" elements with 4*4 + 1 = 17 
c       unique offdiags
c
c     the very last element with 1 offdiag
c
c     summing:
c      13*(nx-1) + 9*(nx-1) + 5*(nz-1) + 17*(nz-2)*(nx-1) + 1
c      (nx-1)*(22 + 17*(nz-2)) + 5*(nz-1) + 1
c     or
c      17*nx*nz -12(nx+nz) + 8
c
c     for the 2.5d case, there are 3 components of motion (u,v,w), so 
c     the invidiual elements are replaced by 3 x 3 matrices.   each row
c     as 26 off diagonal elements. 
c
c     for each node, there will be 9 total (per node) plus the reference
c     to itself.  thus we should have
c
c       nx - 1 first row elements with 9*3 + 3 = 30 unique offdiags
c       nx - 1 last row elements with  9*2 + 3 = 21 unique offdiags
c       nz - 1 last col elements with  9*1 + 3 = 12 unique offdiags
c       (nz - 2)*(nx - 1) "interior" elements with 4*9 + 3 = 39 unique 
c       offdiags
c
c     the very last element with 2 offdiags
c
c     summing:
c       30*(nx-1) + 21*(nx-1) + 12*(nz-1) + 39*(nz-2)*(nx-1) + 2 =
c
c       nx *  30 + 21 - 78 = -27*nx
c       nz *  12 - 39      = -27*nz
c       nx*nz * 39
c       -30 - 21 - 12 + 78 + 3 = 18
c     or
c       39*nx*nz - 27*(nx+nz) + 18
c
      SUBROUTINE srelpl_getsize(nx, nz, numvert, nedge)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, nz
      INTEGER, INTENT(OUT) :: numvert, nedge
      numvert = 3*nx*nz
      nedge = 39*nx*nz - 27*(nx+nz) + 18
      RETURN
      END

      SUBROUTINE srelpl(nx, nz, numvert, nedge, xadj, adjncy)
c 
c     2.5d elastic graph.  
c 
c     input      meaning 
c     -----      ------- 
c     nedge      number of edges 
c     numvert    number of vertices  
c     nx         number of x grid points 
c     nz         number of z grid points 
c   
c     output     meaning 
c     ------     ------- 
c     adjncy     edge vector 
c     xadj       verticy vector  
c
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nedge, numvert, nx, nz
      INTEGER, INTENT(OUT) :: xadj(numvert+1), adjncy(2*nedge)
      INTEGER iedge(27), i, iu, iv, iw, ix, iz, j, iter,
     ;        k, n, nac, nedger, nvt, nvto
c 
c---------------------------------------------------------------------c
c
c     numvert = 3*nx*nz
c     nedge = 39*nx*nz - 27*(nx+nz) + 18
c
c-----loop over all vertices and collect the edges for each one
c
c     note on indexing:
c
c     in the acoustic case, a given element located at row iz and column
c     ix in the model will have the following offdiags;
c       (iz-1, ix-1)  (iz-1, ix)  (iz-1, ix+1)
c       (iz  , ix-1)  (iz  , ix)  (iz  , ix+1)
c       (iz+1, ix-1)  (iz+1, ix)  (iz+1, ix+1)
c     and the absolute index is the z index + nz * (x index - 1)
c
c     for the elastic case, we have u and v parts for each model node.
c
c     the u element at node (iz, ix) will be 
c       3*[iz + nz*(ix-1)] -2
c     the v position is just 1+ that:
c       3*[iz + nz*(ix-1)] -1
c     the w position is just 1+ that:
c       3*[iz + nz*(ix-1)]
c
      iter = 0 
      nvto = 1 
      nedger = 0 
      do i = 1, nx
         do j = 1, nz
c-----------loop over rows;  3 rows per variable; order is u, v, w
            do n = 1, 3
               nvt = 0
c--------------left side
               if (i.gt.1) then
                  ix = i - 1
c-----------------loop over the 3 elements to the left (a1, a2, a3)
                  do k = 1, 3
                     iz = j + k - 2
                     if (iz.gt.0.and.iz.le.nz) then
                        nvt = nvt + 1
                        nac = nz*(ix-1) + iz
                        iu = 3*nac - 2
                        iedge(nvt) = iu
                        nvt = nvt + 1
                        iv = iu + 1
                        iedge(nvt) = iv
                        nvt = nvt + 1
                        iw = iv + 1
                        iedge(nvt) = iw
                     endif
                  enddo
               endif
c--------------middle
               ix = i
c--------------loop over a4, a5, a6; but need to skip over the self 
c--------------reference in a5 (when k = 2)
               do k = 1, 3
                  iz = j + k - 2
                  if (iz.gt.0.and.iz.le.nz) then
c--------------------condition to ignore the diagaonal for u
                     if (n.eq.1) then
                        if (k.ne.2) then
                           nvt = nvt + 1
                           nac = nz*(ix-1) + iz
                           iu = 3*nac - 2    
                           iedge(nvt) = iu 
                        endif
                        nvt = nvt + 1
                        nac = nz*(ix-1) + iz
                        iv = 3*nac - 1
                        iedge(nvt) = iv
                        nvt = nvt + 1
                        iw = iv + 1  
                        iedge(nvt) = iw
                     endif
c--------------------condition to ignore the diagaonal for v
                     if (n.eq.2) then
                        nvt = nvt + 1
                        nac = nz*(ix-1) + iz
                        iu = 3*nac - 2  
                        iedge(nvt) = iu 
                        if (k.ne.2) then
                           nvt = nvt + 1
                           iv = 3*nac - 1
                           iedge(nvt) = iv
                        endif
                        nvt = nvt + 1
                        iw = 3*nac
                        iedge(nvt) = iw
                     endif
c--------------------condition to ignore the diagaonal for w
                     if (n.eq.3) then
                        nvt = nvt + 1
                        nac = nz*(ix-1) + iz
                        iu = 3*nac - 2
                        iedge(nvt) = iu 
                        nvt = nvt + 1
                        iv = iu + 1
                        iedge(nvt) = iv
                        if (k.ne.2) then
                           nvt = nvt + 1
                           iw = 3*nac
                           iedge(nvt) = iw
                        endif
                     endif
                  endif
               enddo
c--------------right side
               if (i.lt.nx) then
                  ix = i + 1
                  do k = 1, 3
                     iz = j + k - 2
                     if (iz.gt.0.and.iz.le.nz) then
                        nvt = nvt + 1
                        nac = iz + nz*(ix-1)   
                        iu = 3*nac - 2    
                        iedge(nvt) = iu
                        nvt = nvt + 1
                        iv = iu + 1    
                        iedge(nvt) = iv
                        nvt = nvt + 1
                        iw = iv + 1    
                        iedge(nvt) = iw
                     endif
                  enddo
               endif
               iter = iter + 1 
               xadj(iter) = nvto 
               nvto = nvto + nvt 
               do k=1,nvt 
                  nedger = nedger + 1 
                  adjncy(nedger) = iedge(k) 
               enddo 
c 200          format(30i10)
            enddo
         enddo !end loop on z 
      enddo !end loop on x 
      iter = iter + 1 
      xadj(iter) = nvto 
      return  
      end


c----This routine makes a graph of the finite difference matrix for input
c       to the Metis package to compute an optimal reordering.
c
c       This version for the 3D elastic case.   Quick note:  a good
c       way to think about ordering is that all the unique offdiags are in the upper
c       (or lower if you like) triangle of the symmetric S matrix.  For a 9
c       point star there will be generally be 4 points in the upper triangle (less
c       if you are at the edge).  If there are N components of motion, then there
c       will be 4*N unique offidags.  Each point will be represented by N rows, which
c       also means that each point in the star is represented by an NxN matix.  
c       We need to specify the upper triangle of the middle matrix, which would b
c       N*(N-1)/2 elements. 
c       Hence the total number for the general case will be:
c               4*N*N + N*(N-1)/2
c       Again, this number decreases at points along the edge of the model, but computation
c       is straightforward once you understand this concept. 
c
      SUBROUTINE srel3d_getsize(nx, ny, nz, numvert, nedge)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, ny, nz
      INTEGER, INTENT(OUT) :: numvert, nedge
      numvert = 3*nx*ny*nz
      nedge = 120*nx*ny*nz
     ;      - 81*(nx*ny + nx*nz + ny*nz)
     ;      + 54*(nx + nz + ny)
     ;      - 36
      RETURN
      END

      SUBROUTINE srel3d(nx, ny, nz, numvert, nedge, xadj, adjncy)
c     input      meaning
c     -----      -------
c     nedge      number of edges
c     numvert    number of vertices
c     nx         number of x grid points
c     nz         number of z grid points
c
c     output     meaning
c     ------     -------
c     adjncy     edge vector
c     xadj       verticy vector
c
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nedge, numvert, nx, ny, nz
      INTEGER, INTENT(OUT) :: xadj(numvert+1), adjncy(2*nedge)

      INTEGER iedge(81), i, igrid, irow, iter, iu, ix, iy, iz,
     ;        j, k, nedger, ngrid, nvt, nvto, xm, xp, ym, yp, zm, zp


c---compute the total number of independant edges

c       numvert = 3*nx*ny*nz
c       nedge = 120*nx*ny*nz - 81*(nx*ny + nx*nz + ny*nz) + 54*(nx + nz + ny) - 36

c----Loop over all vertices and collect the edges for each one
c
c       Note on indexing:
c       
c       The U element at node (iz, ix, iy) will be 
c               3*[iz + nz*(ix-1) + nx*nz*(iy-1)] - 2
c       The V position is just 1+ that:
c               3*[iz + nz*(ix-1) + nx*nz*(iy-1)] - 1
c       The W position is just 1+ that:
c               3*[iz + nz*(ix-1) + nx*nz*(iy-1)]
c
      iter = 0
      nvto = 1
      nedger = 0

c---Loop over 3 x 3 submatrices
      do iy = 1, ny
         ym = MAX0(iy - 1,1)
         yp = MIN0(iy + 1,ny)
         do ix = 1, nx
            xm = MAX0(ix - 1,1)
            xp = MIN0(ix + 1,nx)
            do iz = 1, nz 
               zm = MAX0(iz - 1,1)
               zp = MIN0(iz + 1,nz)
               ngrid = iz + (ix-1)*nz + (iy-1)*nx*nz
c---loop over the 3 rows of the submatrix
               do irow = 1, 3
                  nvt = 0
                  do j = ym, yp
                     do i = xm, xp
                        do k = zm, zp
                           igrid = k + (i-1)*nz + (j-1)*nx*nz
                           iu = 3*igrid - 2
                           if (igrid.ne.ngrid.or.irow.ne.1) then
                              nvt = nvt + 1
                              iedge(nvt) = iu
                           endif
                           iu = iu + 1
                           if (igrid.ne.ngrid.or.irow.ne.2) then
                              nvt = nvt + 1
                              iedge(nvt) = iu
                           endif
                           iu = iu + 1
                           if (igrid.ne.ngrid.or.irow.ne.3) then
                              nvt = nvt + 1
                              iedge(nvt) = iu
                           endif
                        enddo
                     enddo
                  enddo
c----end loop over row elements for this row; save result
                  iter = iter + 1
                  xadj(iter) = nvto
                  nvto = nvto + nvt
                  do k = 1, nvt
                     nedger = nedger + 1
                     adjncy(nedger) = iedge(k)
                  enddo
               enddo
c---end loop over rows of submatrix
            enddo
         enddo
      enddo
c----end loop over grid points

      iter = iter + 1
      xadj(iter) = nvto

      return
      end
