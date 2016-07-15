c
      SUBROUTINE fdfd3d(lbl, iz, ix, iy, ierr)

      USE INIT_MODULE, ONLY : usemin
      USE MODEL_MODULE, ONLY : freesurf, nx, ny, nz
      USE PML_MODULE, ONLY : ipml
      USE INTERFACE_MODULE, ONLY : fixed3d, pfdfd3d, pml3d, tpfree
      IMPLICIT NONE

c arguments:
      INTEGER, INTENT(IN) ::      lbl, iz, ix, iy
      INTEGER, INTENT(OUT) ::     ierr

c       include 'dimension.inc'            
c       include 'init.inc'
c       include 'common.inc'
c       include 'pml.inc'

c
c       input:
c       lbl          three letter designator for the 27 point cube
c       iz, ix, iy   the model node indices (these are (u,v) row pairs in the elastic case)
c
c       Locations of region numbers used in logic below. 
c       These are the parts where iy is not in the PML, correspond to the normal
c       2D and 2.5D case:
c
c        1n       2n       3n
c
c        4n       5n       6n
c
c        7n       8n       9n
c
c       We use "m" for iy < ipml (back) and "p" for iy.le. ny-ipml (front)
c
c---freesurf 1 = top    (iz = 1);   2 = right (ix = 1); 
c            3 = bottom (iz = nz);  4 = left  (ix = nx); 
c            5 = front  (iy = ny);  6 = back  (iy = ny)
c
c       NOTE: This program currently is coded only for the case of freesurf = 1 or 0 (top free or fixed).
c
c      Dirichlet condition (set displacement to zero) used at the far edge of the pml
c NB:  Currently if one side is free and the adjoining one has a pml, we use the pml
c       going into the pml only as the corner BC.  This should be ok, except that the
c       zero displacement condition set at the surface may cause some problems (because
c       we go from free to fixed in one grid interval.   Need to check this.
c_______________________________________________________________________
          
c executable code:
      ierr = 0
c       write(*,*) ' in fdfd3d:, ipml, iz, ix, iy, lbl = ', ipml, iz, ix, iy, lbl
c If this is not an absorbing region (region 5n), then call pfdfd3d
      if ( (ix.gt.ipml).and.(ix.le.(nx-ipml)) 
     +      .and.(iy.gt.ipml).and.(iy.le.(ny-ipml)) 
     +      .and.(iz.gt.ipml).and.(iz.le.(nz-ipml)) ) then
c--if you want to check fidelity of the pml routine, then call it here
c         call pml3d (lbl, iz, ix, iy, ierr)
         call pfdfd3d (lbl, ierr) !iz, ix, iy, ierr)
         if (ierr.ne.0) goto 100

c---fix all the edges but leave the top optional
      elseif (iz.eq.nz .or.
     +        ix.eq.1  .or. ix.eq.nx. or.
     +        iy.eq.1  .or. iy.eq.ny) then
         call fixed3d(lbl) !, iz, ix, iy)

      else
c--otherwise we are in the PML; but we need to ignore it at the free surface
         if (iz.le.ipml) then
            if (iz.eq.1) then
               if (freesurf(1)) then
c  In this case we use Min cells for the free surface   
                  if (usemin) then
                     if (     (ix.gt.ipml).and.(ix.le.(nx-ipml)) 
     +                   .and.(iy.gt.ipml).and.(iy.le.(ny-ipml))) then
                        call pfdfd3d (lbl, ierr) !, iz, ix, iy, ierr)
                        if (ierr.ne.0) goto 100
                     else
c  This call deactivates the iz pml but retains ix, iy.
                        call pml3d (lbl, ipml+1, ix, iy, ierr)
                        if (ierr.ne.0) goto 100
                     endif
                  else
c  Use the traction free conditions instead
                     call tpfree(lbl, ix, iy, ierr) !iz, ix, iy)
                     if (ierr /= 0) goto 100
                  endif
               else
c  Fixed surface
                  call fixed3d(lbl)!, iz, ix, iy)
               endif
            else
c   Case of  ipml >= iz > 1
               if (freesurf(1)) then
                  if (    (ix.gt.ipml).and.(ix.le.(nx-ipml)) 
     +               .and.(iy.gt.ipml).and.(iy.le.(ny-ipml))) then
                     call pfdfd3d (lbl, ierr) !iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 100
                  else
c  This call deactivates the iz pml but retains ix, iy.
                     call pml3d (lbl, ipml+1, ix, iy, ierr)
                     if (ierr.ne.0) goto 100
                  endif
               else
c-- surface not free, use normal PML decay
                  call pml3d (lbl, iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 100
               endif
            endif
         else
c  All other boundaries; use normal PML
            call pml3d (lbl, iz, ix, iy, ierr)
            if (ierr.ne.0) goto 100
         endif
      endif

      return
c 
c-----error section: 
100   continue  
      write(*,*) 'Error occurred in fdfd3d.f, stopping!' 

      return
      end


