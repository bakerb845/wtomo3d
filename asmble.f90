      SUBROUTINE ASMBLE(lunlog, nzero, nx, ny, nz, &
                        irptr, jcptr, a, ierr)
! 
!     input      meaning 
!     -----      ------- 
!     irptr      CSR row pointer 
!     jcptr      CSR column pointer 
!     lunlog     log file id number
!     nx         number of x grid points 
!     ny         number of y grid points 
!     nz         number of z grid points 
!     nzero      number of non-zeros 
!     irptr      CSR row pointer (3*nx*ny*nz+1)
!  
!     output     meaning 
!     ------     ------- 
!     a          assembled matrix 
!     ierr       error flag, != 0 error occurred 
!  
!     Handles the assembly for the finite difference matrix.  
!     What we do is loop on the grid, with the fast direction in 
!     z and next direction in x, then slow in y, then we loop over the components, 
!     u,v,and w for that grid point.  We therefore generate 
!     nx*ny*nz equations as will be dictated by our row counter 
!     irow. :then for each row we wish to generate the up to 81
!     values associated with the finite difference, 
!
!       In the 2D and 2.5D cases we do the following:
! 
!     For the u component we generate a row in the matrix: 
!     ...a1uu a1uv a1uw a2uu a2uv a2uw a3uu a3uv a3uw... 
!     ...a4uu a4uv a4uw a5uu a5uv a5uw a6uu a6uv a6uw... 
!     ...a7uu a7uv a7uw a8uu a8uv a8uw a9uu a9uv a9uw...
! 
!     For the v compnoent we generate a row in the matrix: 
!     ...a1vu a1vv a1vw a2vu a2vv a2vw a3vu a3uv a3vw... 
!     ...a4vu a4vv a4vw a5vu a5vv a5vw a6vu a6vv a6vw...
!     ...a7vu a7vv z7vw a8vu a8vv a8vw a9vu a9vv a9vw... 
! 
!     and for the w component we generate a row in the matrix: 
!     ...a1wu a1wv a1ww a2wu a2wv a2ww a3wu a3uv a3ww...
!     ...a4wu a4wv a4ww a5wu a5wv a5ww a6wu a6wv a6ww... 
!     ...a7wu a7wv z7ww a8wu a8wv a8ww a9wu a9wv a9ww...
!
!       In the 3D case, we expand the 2D by adding  3 x 3 matrices of elements
!       at j-1 and j+1 (the original 2D would be at j).  So for example a row would
!       look like this:
!
!     ...a1m a2m a3m...a4m a5m a6m...a7m a8m a9m...
!         ...a1n a2n a3n...a4n a5n a6n...a7n a8n a9n...
!               ...a1p a2p a3p...a4p a5p a6p...a7p a8p a9p...
!
!       Each of the above 27 entries is actauly a 3 x 3 matrix, so in fact there
!       are three real rows shown above.   Each also has three columns, for example
!
!             |a1mUU  a1mUV  a1mUW|
!       a1m = |a1mVU  a1mVV  a1mVW|
!             |a1mWU  a1mWV  a1mWW|
! 
!       So 27 matrices will result in 81 elements per row.
!
!     Implicitly assumed is that irn is ordered from 1 to 3*nx*ny*nz 
!     and jcn is ordered from low to high, realize, that since metis 
!     is called from mumps it is not necessary to use the permutation 
!     vector iperm 
! 
      USE MAT_MODULE
      USE FD_ENUM_MODULE, ONLY : bem, ben, bep, &
                                 adm, adn, adp, &
                                 aam, aan, aap, &
                                 afm, afn, afp, &
                                 ddm, ddn, ddp, &
                                 ffm, ffn, ffp, &
                                 cdm, cdn, cdp, &
                                 ccm, ccn, ccp, &
                                 cfm, cfn, cfp
      USE INTERFACE_MODULE, ONLY : fdfd3d, setcoefs3d, zeroas
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nzero, lunlog 
      INTEGER, INTENT(IN) :: irptr(3*nx*ny*nz+1), jcptr(nzero), &
                             nx, ny, nz
      COMPLEX, INTENT(OUT) :: a(nzero)
      INTEGER, INTENT(OUT) :: ierr
      INTEGER jstar(27), irow, ix, iy, iz,    & 
              i, ierr1, nnz, izero, inz, j, loc, itest 
      INTEGER, PARAMETER :: nsd = 3


      COMPLEX uline(81), vline(81), wline(81) 
      COMPLEX, PARAMETER :: zzero = CMPLX(0.0, 0.0) 

! 
!---------------------------------------------------------------------!
   
!.... loop on grid 
      ierr = 0 
      izero = 0 
      irow = 0
      do 1 iy = 1,ny 
         do 2 ix = 1,nx 
            do 3 iz = 1,nz 
 
!.......... null out the lines of the matrix 
               CALL ZEROAS()
               do i = 1,  81
                  uline(i) = zzero !CMPLX(0.,0.) 
                  vline(i) = zzero !CMPLX(0.,0.) 
                  wline(i) = zzero !CMPLX(0.,0.) 
               enddo
               do i = 1, 27 
                  jstar(i) = 0 
               enddo

!.......... now generate values for each of the matrix  
               !igrd = igrd + 1 
               !igrd = (iy - 1)*nx*nz + (ix - 1)*nz + iz
               nnz = 0
               CALL SETCOEFS3D(iz, ix, iy, 0)   !find constants for this grid point

               if (iy.gt.1) then
                  CALL FDFD3D(bem, iz, ix, iy, ierr1) !finite difference 27 pt. cube
                  IF (ierr1 /= 0) ierr = ierr + 1 !goto 900
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

               CALL FDFD3D(ben, iz, ix, iy, ierr1)
               IF (ierr1 /= 0) ierr = ierr + 1
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
                  CALL FDFD3D(bep, iz, ix, iy, ierr1)
                  IF (ierr1 /= 0) ierr = ierr + 1
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
! 
!.......... not on the left side of the model   
               if (ix.gt.1) then 
                  if (iz.gt.1) then
                     if (iy.gt.1) then
                        CALL FDFD3D(adm, iz, ix, iy, ierr1)
                        IF (ierr1 /= 0) ierr = ierr + 1
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

                     CALL FDFD3D(adn, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
                        CALL FDFD3D(adp, iz, ix, iy, ierr1)
                        IF (ierr1 /= 0) ierr = ierr + 1
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
                     CALL FDFD3D(aam, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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

                  CALL FDFD3D(aan, iz, ix, iy, ierr1)
                  IF (ierr1 /= 0) ierr = ierr + 1
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
                     CALL FDFD3D(aap, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
                        CALL FDFD3D(afm, iz, ix, iy, ierr1)
                        IF (ierr1 /= 0) ierr = ierr + 1
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

                     CALL FDFD3D(afn, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
                        CALL FDFD3D(afp, iz, ix, iy, ierr1)
                        IF (ierr1 /= 0) ierr = ierr + 1
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
                     CALL FDFD3D(ddm, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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

                  CALL FDFD3D(ddn, iz, ix, iy, ierr1)
                  IF (ierr1 /= 0) ierr = ierr + 1
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
                     CALL FDFD3D(ddp, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
! 
!.......... not on bottom
               if (iz.lt.nz) then
                  if (iy.gt.1) then 
                     CALL FDFD3D(ffm, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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

                  CALL FDFD3D(ffn, iz, ix, iy, ierr1)
                  IF (ierr1 /= 0) ierr = ierr + 1
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
                     CALL FDFD3D(ffp, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
! 
!.......... not on right side 
               if (ix.lt.nx) then 
                  if (iz.gt.1) then 
                     if (iy.gt.1) then
                        CALL FDFD3D(cdm, iz, ix, iy, ierr1)
                        IF (ierr1 /= 0) ierr = ierr + 1
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

                     CALL FDFD3D(cdn, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
                        CALL FDFD3D(cdp, iz, ix, iy, ierr1)
                        IF (ierr1 /= 0) ierr = ierr + 1
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
                     CALL FDFD3D(ccm, iz, ix, iy, ierr1)
                     IF (ierr1 /= 0) ierr = ierr + 1
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

                  CALL FDFD3D(ccn, iz, ix, iy, ierr1)
                  IF (ierr1 /= 0) ierr = ierr + 1
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
                     CALL FDFD3D(ccp, iz, ix, iy, ierr)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
                        CALL FDFD3D(cfm, iz, ix, iy, ierr1)
                        IF (ierr1 /= 0) ierr = ierr + 1
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

                     CALL FDFD3D(cfn, iz, ix, iy, ierr)
                     IF (ierr1 /= 0) ierr = ierr + 1
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
                        CALL FDFD3D(cfp, iz, ix, iy, ierr)
                        IF (ierr1 /= 0) ierr = ierr + 1
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
! 
!.......... for each row we have have 3 components so fill the line  
              do 5 i=1,nsd 
! 
!............. intially we check that this row fits, and pointers are
!............. correctly ordered  
                 irow = irow + 1 
                 !if (nzeror(irow).ne.nnz) then 
                 if (irptr(irow+1) - irptr(irow) /= nnz) then
                    ierr = ierr + 1
                    write(*,*) 'Error in asmble: nnz incorrect', irow
                    !write(*,*) 'Error in asmble: incorrect numbers of non-zeros in row', irow
                    !ierr = 1 
                    !return
                 endif  
! 
!............. for each non-zero insert into appropriate location in 
!............. global stiffness 
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
   60              continue 
    7            continue 
    5         continue 

    3       continue ! loop on z
    2    continue ! loop on x
    1 continue ! loop on y 

      itest = 0
      if (itest.eq.1) then
         open(unit=30,file='mymat.txt')
         print *, 'nzero =',nzero
         !do i=1,nzero
         do i=1,3*nx*ny*nz
            do j=irptr(i),irptr(i+1)-1
               write(30,*) i, jcptr(j), a(j) !irn(i), jcn(i), a(i)
            enddo
         enddo
         close(30)
         return
         !CALL exit_mpi(' Done with test in asmble.f so stopping ', 15, mpierr)
      endif
! 
!.... error checks
      if (nzero.ne.izero) then 
         write(*,*)      'Predicted number of non-zeros is',nzero
         write(*,*)      'Number of non-zeros found in matrix is ',izero
         write(lunlog,*) 'Predicted number of non-zeros is',nzero
         write(lunlog,*) 'Number of non-zeros found in matrix is ',izero
         ierr = 1 
         return
      endif 

      if (ierr.ne.0) then 
         write(*,*)      'ERROR calling fdfd3d from ASBMLE'
         write(lunlog,*) 'ERROR calling fdfd3d from ASBMLE' 
         return
      endif

      return 
      end 

