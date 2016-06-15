      SUBROUTINE asmble(lunlog, nzero, nx, ny, nz, nzeror,
     ;                  irn, jcn, a, ierr)
c 
c     input      meaning 
c     -----      ------- 
c     irn        row pointer 
c     jcn        column pointer 
c     lunlog     log file id number
c     nx         number of x grid points 
c     ny         number of y grid points 
c     nz         number of z grid points 
c     nzero      number of non-zeros 
c     nzeror     number of non zeros in each row of matrix
c  
c     output     meaning 
c     ------     ------- 
c     a          assembled matrix 
c     ierr       error flag, != 0 error occurred 
c  
c     Handles the assembly for the finite difference matrix.  
c     What we do is loop on the grid, with the fast direction in 
c     z and next direction in x, then slow in y, then we loop over the components, 
c     u,v,and w for that grid point.  We therefore generate 
c     nx*ny*nz equations as will be dictated by our row counter 
c     irow. :then for each row we wish to generate the up to 81
c     values associated with the finite difference, 
c
c       In the 2D and 2.5D cases we do the following:
c 
c     For the u component we generate a row in the matrix: 
c     ...a1uu a1uv a1uw a2uu a2uv a2uw a3uu a3uv a3uw... 
c     ...a4uu a4uv a4uw a5uu a5uv a5uw a6uu a6uv a6uw... 
c     ...a7uu a7uv a7uw a8uu a8uv a8uw a9uu a9uv a9uw...
c 
c     For the v compnoent we generate a row in the matrix: 
c     ...a1vu a1vv a1vw a2vu a2vv a2vw a3vu a3uv a3vw... 
c     ...a4vu a4vv a4vw a5vu a5vv a5vw a6vu a6vv a6vw...
c     ...a7vu a7vv z7vw a8vu a8vv a8vw a9vu a9vv a9vw... 
c 
c     and for the w component we generate a row in the matrix: 
c     ...a1wu a1wv a1ww a2wu a2wv a2ww a3wu a3uv a3ww...
c     ...a4wu a4wv a4ww a5wu a5wv a5ww a6wu a6wv a6ww... 
c     ...a7wu a7wv z7ww a8wu a8wv a8ww a9wu a9wv a9ww...
c
c       In the 3D case, we expand the 2D by adding  3 x 3 matrices of elements
c       at j-1 and j+1 (the original 2D would be at j).  So for example a row would
c       look like this:
c
c     ...a1m a2m a3m...a4m a5m a6m...a7m a8m a9m...
c         ...a1n a2n a3n...a4n a5n a6n...a7n a8n a9n...
c               ...a1p a2p a3p...a4p a5p a6p...a7p a8p a9p...
c
c       Each of the above 27 entries is actauly a 3 x 3 matrix, so in fact there
c       are three real rows shown above.   Each also has three columns, for example
c
c             |a1mUU  a1mUV  a1mUW|
c       a1m = |a1mVU  a1mVV  a1mVW|
c             |a1mWU  a1mWV  a1mWW|
c 
c       So 27 matrices will result in 81 elements per row.

c     Implicitly assumed is that irn is ordered from 1 to 3*nx*ny*nz 
c     and jcn is ordered from low to high, realize, that since metis 
c     is called from mumps it is not necessary to use the permutation 
c     vector iperm 
c 
      USE MAT_MODULE
      USE INTERFACE_MODULE, ONLY : fdfd3d, setcoefs3d, zeroas
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nzero, lunlog 
      INTEGER, INTENT(IN) :: irn(nzero), jcn(nzero),
     ;                       nzeror(3*nx*ny*nz), nx, ny, nz
      COMPLEX, INTENT(OUT) :: a(nzero)
      INTEGER, INTENT(OUT) :: ierr
      INTEGER jstar(27), igrd,irow,igrd1,igrd2, ix,iy,iz, 
     ;        i,nnz,izero, inz,loc,
     ;        itest  
c     integer mpierr
      INTEGER, PARAMETER :: nsd = 3


      COMPLEX uline(81), vline(81), wline(81) 

c 
c---------------------------------------------------------------------c
   
c.... loop on grid 
      ierr = 0 
      izero = 0 
      igrd = 0 
      irow = 0
      igrd1 = 1 
      igrd2 = nx*ny*nz 
      do 3 iy = 1,ny 
         do 2 ix = 1,nx 
            do 1 iz = 1,nz 
 
c.......... null out the lines of the matrix 
               CALL zeroas()
               do i = 1,  81
                  uline(i) = CMPLX(0.,0.) 
                  vline(i) = CMPLX(0.,0.) 
                  wline(i) = CMPLX(0.,0.) 
               enddo
               do i = 1, 27 
                  jstar(i) = 0 
               enddo

c.......... now generate values for each of the matrix  
               igrd = igrd + 1 
               nnz = 0
               CALL setcoefs3d (iz,ix,iy,0)   !find constants for this grid point

               if (iy.gt.1) then
                  CALL fdfd3d('bem',iz, ix, iy, ierr) !finite difference 27 pt. cube
                  if (ierr.ne.0) goto 900
                  uline(13) = a5muu 
                  uline(14) = a5muv 
                  uline(15) = a5muw 
                  vline(13) = a5mvu 
                  vline(14) = a5mvv 
                  vline(15) = a5mvw 
                  wline(13) = a5mwu 
                  wline(14) = a5mwv 
                  wline(15) = a5mww 
                  jstar(5) = 1 
                  nnz = nnz + 3
               endif

               CALL fdfd3d('ben',iz, ix, iy, ierr)
               uline(40) = a5nuu 
               uline(41) = a5nuv 
               uline(42) = a5nuw 
               vline(40) = a5nvu 
               vline(41) = a5nvv 
               vline(42) = a5nvw 
               wline(40) = a5nwu 
               wline(41) = a5nwv 
               wline(42) = a5nww 
               jstar(14) = 1 
               nnz = nnz + 3

               if (iy.lt.ny) then
                  CALL fdfd3d('bep',iz, ix, iy, ierr)
                  uline(67) = a5puu 
                  uline(68) = a5puv 
                  uline(69) = a5puw 
                  vline(67) = a5pvu 
                  vline(68) = a5pvv 
                  vline(69) = a5pvw 
                  wline(67) = a5pwu 
                  wline(68) = a5pwv 
                  wline(69) = a5pww 
                  jstar(23) = 1 
                  nnz = nnz + 3
               endif
c 
c.......... not on the left side of the model   
               if (ix.gt.1) then 
                  if (iz.gt.1) then
                     if (iy.gt.1) then
                        CALL fdfd3d ('adm', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900     
                        uline(1)  = a1muu
                        uline(2)  = a1muv
                        uline(3)  = a1muw
                        vline(1)  = a1mvu
                        vline(2)  = a1mvv
                        vline(3)  = a1mvw 
                        wline(1)  = a1mwu 
                        wline(2)  = a1mwv
                        wline(3)  = a1mww  
                        jstar(1) = 1 
                        nnz = nnz + 3 
                     endif

                     CALL fdfd3d ('adn', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900     
                     uline(28) = a1nuu
                     uline(29) = a1nuv
                     uline(30) = a1nuw
                     vline(28) = a1nvu
                     vline(29) = a1nvv
                     vline(30) = a1nvw 
                     wline(28) = a1nwu 
                     wline(29) = a1nwv
                     wline(30) = a1nww  
                     jstar(10) = 1 
                     nnz = nnz + 3 

                     if (iy.lt.ny) then
                        CALL fdfd3d ('adp', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900     
                        uline(55) = a1puu
                        uline(56) = a1puv
                        uline(57) = a1puw
                        vline(55) = a1pvu
                        vline(56) = a1pvv
                        vline(57) = a1pvw 
                        wline(55) = a1pwu 
                        wline(56) = a1pwv
                        wline(57) = a1pww  
                        jstar(19) = 1 
                        nnz = nnz + 3 
                     endif
                  endif 

                  if (iy.gt.1) then 
                     CALL fdfd3d ('aam', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(4)  = a2muu
                     uline(5)  = a2muv
                     uline(6)  = a2muw
                     vline(4)  = a2mvu
                     vline(5)  = a2mvv
                     vline(6)  = a2mvw 
                     wline(4)  = a2mwu 
                     wline(5)  = a2mwv
                     wline(6)  = a2mww  
                     jstar(2) = 1 
                     nnz = nnz + 3 
                  endif

                  CALL fdfd3d ('aan', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  uline(31) = a2nuu
                  uline(32) = a2nuv
                  uline(33) = a2nuw
                  vline(31) = a2nvu
                  vline(32) = a2nvv
                  vline(33) = a2nvw 
                  wline(31) = a2nwu 
                  wline(32) = a2nwv
                  wline(33) = a2nww  
                  jstar(11) = 1 
                  nnz = nnz + 3 

                  if (iy.lt.ny) then 
                     CALL fdfd3d ('aap', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(58) = a2puu
                     uline(59) = a2puv
                     uline(60) = a2puw
                     vline(58) = a2pvu
                     vline(59) = a2pvv
                     vline(60) = a2pvw 
                     wline(58) = a2pwu 
                     wline(59) = a2pwv
                     wline(60) = a2pww  
                     jstar(20) = 1 
                     nnz = nnz + 3 
                  endif 

                  if (iz.lt.nz) then
                     if (iy.gt.1) then 
                        CALL fdfd3d ('afm', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900
                        uline(7) = a3muu
                        uline(8) = a3muv
                        uline(9) = a3muw
                        vline(7) = a3mvu
                        vline(8) = a3mvv
                        vline(9) = a3mvw 
                        wline(7) = a3mwu 
                        wline(8) = a3mwv
                        wline(9) = a3mww  
                        jstar(3) = 1 
                        nnz = nnz + 3 
                     endif

                     CALL fdfd3d ('afn', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(34) = a3nuu
                     uline(35) = a3nuv
                     uline(36) = a3nuw
                     vline(34) = a3nvu
                     vline(35) = a3nvv
                     vline(36) = a3nvw 
                     wline(34) = a3nwu 
                     wline(35) = a3nwv
                     wline(36) = a3nww  
                     jstar(12) = 1 
                     nnz = nnz + 3 

                     if (iy.lt.ny) then 
                        CALL fdfd3d ('afp', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900
                        uline(61) = a3puu
                        uline(62) = a3puv
                        uline(63) = a3puw
                        vline(61) = a3pvu
                        vline(62) = a3pvv
                        vline(63) = a3pvw 
                        wline(61) = a3pwu 
                        wline(62) = a3pwv
                        wline(63) = a3pww  
                        jstar(21) = 1 
                        nnz = nnz + 3 
                     endif
                  endif 
               endif
               if (iz.gt.1) then 
                  if (iy.gt.1) then 
                     CALL fdfd3d ('ddm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(10) = a4muu
                     uline(11) = a4muv
                     uline(12) = a4muw
                     vline(10) = a4mvu
                     vline(11) = a4mvv
                     vline(12) = a4mvw 
                     wline(10) = a4mwu 
                     wline(11) = a4mwv
                     wline(12) = a4mww 
                     jstar(4) = 1 
                     nnz = nnz + 3 
                  endif

                  CALL fdfd3d ('ddn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  uline(37) = a4nuu
                  uline(38) = a4nuv
                  uline(39) = a4nuw
                  vline(37) = a4nvu
                  vline(38) = a4nvv
                  vline(39) = a4nvw 
                  wline(37) = a4nwu 
                  wline(38) = a4nwv
                  wline(39) = a4nww 
                  jstar(13) = 1 
                  nnz = nnz + 3 

                  if (iy.lt.ny) then 
                     CALL fdfd3d ('ddp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(64) = a4puu
                     uline(65) = a4puv
                     uline(66) = a4puw
                     vline(64) = a4pvu
                     vline(65) = a4pvv
                     vline(66) = a4pvw 
                     wline(64) = a4pwu 
                     wline(65) = a4pwv
                     wline(66) = a4pww 
                     jstar(22) = 1 
                     nnz = nnz + 3 
                  endif
               endif 
c 
c.......... not on bottom
               if (iz.lt.nz) then
                  if (iy.gt.1) then 
                     CALL fdfd3d ('ffm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(16) = a6muu
                     uline(17) = a6muv
                     uline(18) = a6muw
                     vline(16) = a6mvu
                     vline(17) = a6mvv
                     vline(18) = a6mvw
                     wline(16) = a6mwu
                     wline(17) = a6mwv
                     wline(18) = a6mww
                     jstar(6) = 1 
                     nnz = nnz + 3 
                  endif

                  CALL fdfd3d ('ffn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  uline(43) = a6nuu
                  uline(44) = a6nuv
                  uline(45) = a6nuw
                  vline(43) = a6nvu
                  vline(44) = a6nvv
                  vline(45) = a6nvw
                  wline(43) = a6nwu
                  wline(44) = a6nwv
                  wline(45) = a6nww
                  jstar(15) = 1 
                  nnz = nnz + 3 

                  if (iy.lt.ny) then 
                     CALL fdfd3d ('ffp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(70) = a6puu
                     uline(71) = a6puv
                     uline(72) = a6puw
                     vline(70) = a6pvu
                     vline(71) = a6pvv
                     vline(72) = a6pvw
                     wline(70) = a6pwu
                     wline(71) = a6pwv
                     wline(72) = a6pww
                     jstar(24) = 1 
                     nnz = nnz + 3 
                  endif
               endif 
c 
c.......... not on right side 
               if (ix.lt.nx) then 
                  if (iz.gt.1) then 
                     if (iy.gt.1) then
                        CALL fdfd3d ('cdm', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900
                        uline(19) = a7muu
                        uline(20) = a7muv
                        uline(21) = a7muw
                        vline(19) = a7mvu
                        vline(20) = a7mvv
                        vline(21) = a7mvw 
                        wline(19) = a7mwu 
                        wline(20) = a7mwv
                        wline(21) = a7mww 
                        jstar(7) = 1 
                        nnz = nnz + 3 
                     endif

                     CALL fdfd3d ('cdn', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(46) = a7nuu
                     uline(47) = a7nuv
                     uline(48) = a7nuw
                     vline(46) = a7nvu
                     vline(47) = a7nvv
                     vline(48) = a7nvw 
                     wline(46) = a7nwu 
                     wline(47) = a7nwv
                     wline(48) = a7nww 
                     jstar(16) = 1 
                     nnz = nnz + 3 

                     if (iy.lt.ny) then
                        CALL fdfd3d ('cdp', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900
                        uline(73) = a7puu
                        uline(74) = a7puv
                        uline(75) = a7puw
                        vline(73) = a7pvu
                        vline(74) = a7pvv
                        vline(75) = a7pvw 
                        wline(73) = a7pwu 
                        wline(74) = a7pwv
                        wline(75) = a7pww 
                        jstar(25) = 1 
                        nnz = nnz + 3 
                     endif
                  endif   
                  if (iy.gt.1) then
                     CALL fdfd3d ('ccm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(22) = a8muu
                     uline(23) = a8muv
                     uline(24) = a8muw
                     vline(22) = a8mvu
                     vline(23) = a8mvv
                     vline(24) = a8mvw
                     wline(22) = a8mwu
                     wline(23) = a8mwv
                     wline(24) = a8mww
                     jstar(8) = 1 
                     nnz = nnz + 3 
                  endif

                  CALL fdfd3d ('ccn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  uline(49) = a8nuu
                  uline(50) = a8nuv
                  uline(51) = a8nuw
                  vline(49) = a8nvu
                  vline(50) = a8nvv
                  vline(51) = a8nvw
                  wline(49) = a8nwu
                  wline(50) = a8nwv
                  wline(51) = a8nww
                  jstar(17) = 1 
                  nnz = nnz + 3 

                  if (iy.lt.ny) then
                     CALL fdfd3d ('ccp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(76) = a8puu
                     uline(77) = a8puv
                     uline(78) = a8puw
                     vline(76) = a8pvu
                     vline(77) = a8pvv
                     vline(78) = a8pvw
                     wline(76) = a8pwu
                     wline(77) = a8pwv
                     wline(78) = a8pww
                     jstar(26) = 1 
                     nnz = nnz + 3 
                  endif
                  if (iz.lt.nz) then
                     if (iy.gt.1) then
                        CALL fdfd3d ('cfm', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900
                        uline(25) = a9muu
                        uline(26) = a9muv
                        uline(27) = a9muw
                        vline(25) = a9mvu
                        vline(26) = a9mvv
                        vline(27) = a9mvw 
                        wline(25) = a9mwu 
                        wline(26) = a9mwv
                        wline(27) = a9mww 
                        jstar(9) = 1 
                        nnz = nnz + 3 
                     endif

                     CALL fdfd3d ('cfn', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     uline(52) = a9nuu
                     uline(53) = a9nuv
                     uline(54) = a9nuw
                     vline(52) = a9nvu
                     vline(53) = a9nvv
                     vline(54) = a9nvw 
                     wline(52) = a9nwu 
                     wline(53) = a9nwv
                     wline(54) = a9nww 
                     jstar(18) = 1 
                     nnz = nnz + 3 

                     if (iy.lt.ny) then
                        CALL fdfd3d ('cfp', iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 900
                        uline(79) = a9puu
                        uline(80) = a9puv
                        uline(81) = a9puw
                        vline(79) = a9pvu
                        vline(80) = a9pvv
                        vline(81) = a9pvw 
                        wline(79) = a9pwu 
                        wline(80) = a9pwv
                        wline(81) = a9pww 
                        jstar(27) = 1 
                        nnz = nnz + 3 
                    endif
                 endif 
              endif 
c 
c.......... for each row we have have 3 components so fill the line  
              do 5 i=1,nsd 
c 
c............. intially we check that this row fits, and pointers are
c............. correctly ordered  
                 irow = irow + 1 
                 if (nzeror(irow).ne.nnz) then 
                    write(*,*) 'Error in asmble: nnz incorrect', irow
                    !write(*,*) 'Error in asmble: incorrect numbers of non-zeros in row', irow
                    ierr = 1 
                    return 
                 endif  
c 
c............. for each non-zero insert into appropriate location in 
c............. global stiffness 
                 do 7 inz = 1, 27 
                    if (jstar(inz).eq.0) go to 60 !term inactive 
                    loc = inz*3 - 2 
                    if (i.eq.1) then 
                       izero = izero + 1 
                       a(izero) = uline(loc)  
                       izero = izero + 1 
                       a(izero) = uline(loc + 1) 
                       izero = izero + 1 
                       a(izero) = uline(loc + 2) 
                    endif
                    if (i.eq.2) then 
                       izero = izero + 1 
                       a(izero) = vline(loc) 
                       izero = izero + 1 
                       a(izero) = vline(loc + 1)
                       izero = izero + 1 
                       a(izero) = vline(loc + 2)
                   endif 
                   if (i.eq.3) then 
                      izero = izero + 1 
                      a(izero) = wline(loc)  
                      izero = izero + 1 
                      a(izero) = wline(loc + 1) 
                      izero = izero + 1 
                      a(izero) = wline(loc + 2)
                   endif
60                 continue 
7                continue 
5             continue 

1           continue 
2        continue 
3     continue 

      itest = 0
      if (itest.eq.1) then
         open(unit=30,file='mymat.txt')
         print *, 'nzero =',nzero
         do i=1,nzero
            write(30,*) irn(i), jcn(i), a(i)
         enddo
         close(30)
         return
         !CALL exit_mpi(' Done with test in asmble.f so stopping ', 15, mpierr)
      endif
c 
c.... error checks
      if (nzero.ne.izero) then 
         write(*,*)      'Predicted number of non-zeros is',nzero
         write(*,*)      'Number of non-zeros found in matrix is ',izero
         write(lunlog,*) 'Predicted number of non-zeros is',nzero
         write(lunlog,*) 'Number of non-zeros found in matrix is ',izero
         ierr = 1 
         return
      endif 

  900 continue 

      if (ierr.ne.0) then 
         write(*,*)      'ERROR calling fdfd3d from ASBMLE'
         write(lunlog,*) 'ERROR calling fdfd3d from ASBMLE' 
         return
      endif

      return 
      end 

