c-------------------------------------------------------------------------------------
c       hetfldsrc.f
c
c       author: s. roecker 9/09
c
c       this routine calculates the source distribution required for 
c       a scattered field from
c
c               suh = f = -(s - so)uo = (so - s)uo
c
c       we solve the above for uh and then add on uo to retrive the total field.
c
c       all terms are in general complex, although we read the model files in as real to begin with.
c
c

      SUBROUTINE hetfldsrc (srcfn, sfld, myid, master, ierr)

      USE INIT_MODULE, ONLY : nx, ny, nz
      USE WAVEFIELD_MODULE, ONLY : utbkr, vtbkr, wtbkr
      USE MAT_MODULE
      USE MATBG_MODULE
      USE INTERFACE_MODULE, ONLY : fdfd3d, setcoefs3d
      IMPLICIT NONE

c---need to check to see which of these is really needed

      COMPLEX, INTENT(IN) :: srcfn
      COMPLEX, DIMENSION(:), INTENT(OUT) :: sfld
      INTEGER, INTENT(IN) :: myid, master
      INTEGER, INTENT(OUT) :: ierr 

c     INTEGER itest
        
c---local variables
      COMPLEX us, vs, ws
      INTEGER igrd, indx, iz, ix, iy

c---collect all the terms for the current model
      ierr = 0
      if (myid.eq.master) then
         write(*,*)'HETFLDSRC: Looping over background coefficients...'
         !write(lunlog,*) 'HETFLDSRC: Looping over background coefficients'
      endif
      sfld(:) = CMPLX(0.0, 0.0)
      do iy = 1, ny
         do ix = 1, nx
            do iz = 1, nz
               us = CMPLX(0.0, 0.0)
               vs = CMPLX(0.0, 0.0)
               ws = CMPLX(0.0, 0.0)

               CALL setcoefs3d (iz, ix, iy, 0)

               CALL fdfd3d ('ben', iz, ix, iy, ierr)
               if (ierr.ne.0) goto 900
               a5nuu1 = a5nuu
               a5nuv1 = a5nuv
               a5nuw1 = a5nuw
               a5nvu1 = a5nvu
               a5nvv1 = a5nvv
               a5nvw1 = a5nvw
               a5nwu1 = a5nwu
               a5nwv1 = a5nwv
               a5nww1 = a5nww
               if (iy.gt.1) then
                  CALL fdfd3d ('bem', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a5muu1 = a5muu
                  a5muv1 = a5muv
                  a5muw1 = a5muw
                  a5mvu1 = a5mvu
                  a5mvv1 = a5mvv
                  a5mvw1 = a5mvw
                  a5mwu1 = a5mwu
                  a5mwv1 = a5mwv
                  a5mww1 = a5mww
               endif
               if (iy.lt.ny) then
                  CALL fdfd3d ('bep', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a5puu1 = a5puu
                  a5puv1 = a5puv
                  a5puw1 = a5puw
                  a5pvu1 = a5pvu
                  a5pvv1 = a5pvv
                  a5pvw1 = a5pvw
                  a5pwu1 = a5pwu
                  a5pwv1 = a5pwv
                  a5pww1 = a5pww
               endif
               if (ix.gt.1.and.iz.gt.1) then
                  CALL fdfd3d ('adn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a1nuu1 = a1nuu
                  a1nuv1 = a1nuv
                  a1nuw1 = a1nuw
                  a1nvu1 = a1nvu
                  a1nvv1 = a1nvv
                  a1nvw1 = a1nvw
                  a1nwu1 = a1nwu
                  a1nwv1 = a1nwv
                  a1nww1 = a1nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('adm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a1muu1 = a1muu
                     a1muv1 = a1muv
                     a1muw1 = a1muw
                     a1mvu1 = a1mvu
                     a1mvv1 = a1mvv
                     a1mvw1 = a1mvw
                     a1mwu1 = a1mwu
                     a1mwv1 = a1mwv
                     a1mww1 = a1mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('adp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a1puu1 = a1puu
                     a1puv1 = a1puv
                     a1puw1 = a1puw
                     a1pvu1 = a1pvu
                     a1pvv1 = a1pvv
                     a1pvw1 = a1pvw
                     a1pwu1 = a1pwu
                     a1pwv1 = a1pwv
                     a1pww1 = a1pww
                  endif
               endif
               if (ix.gt.1) then
                  CALL fdfd3d ('aan', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a2nuu1 = a2nuu
                  a2nuv1 = a2nuv
                  a2nuw1 = a2nuw
                  a2nvu1 = a2nvu
                  a2nvv1 = a2nvv
                  a2nvw1 = a2nvw
                  a2nwu1 = a2nwu
                  a2nwv1 = a2nwv
                  a2nww1 = a2nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('aam', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a2muu1 = a2muu
                     a2muv1 = a2muv
                     a2muw1 = a2muw
                     a2mvu1 = a2mvu
                     a2mvv1 = a2mvv
                     a2mvw1 = a2mvw
                     a2mwu1 = a2mwu
                     a2mwv1 = a2mwv
                     a2mww1 = a2mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('aap', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a2puu1 = a2puu
                     a2puv1 = a2puv
                     a2puw1 = a2puw
                     a2pvu1 = a2pvu
                     a2pvv1 = a2pvv
                     a2pvw1 = a2pvw
                     a2pwu1 = a2pwu
                     a2pwv1 = a2pwv
                     a2pww1 = a2pww
                  endif
               endif
               if (ix.gt.1.and.iz.lt.nz) then
                  CALL fdfd3d ('afn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a3nuu1 = a3nuu
                  a3nuv1 = a3nuv
                  a3nuw1 = a3nuw
                  a3nvu1 = a3nvu
                  a3nvv1 = a3nvv
                  a3nvw1 = a3nvw
                  a3nwu1 = a3nwu
                  a3nwv1 = a3nwv
                  a3nww1 = a3nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('afm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a3muu1 = a3muu
                     a3muv1 = a3muv
                     a3muw1 = a3muw
                     a3mvu1 = a3mvu
                     a3mvv1 = a3mvv
                     a3mvw1 = a3mvw
                     a3mwu1 = a3mwu
                     a3mwv1 = a3mwv
                     a3mww1 = a3mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('afp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a3puu1 = a3puu
                     a3puv1 = a3puv
                     a3puw1 = a3puw
                     a3pvu1 = a3pvu
                     a3pvv1 = a3pvv
                     a3pvw1 = a3pvw
                     a3pwu1 = a3pwu
                     a3pwv1 = a3pwv
                     a3pww1 = a3pww
                  endif
               endif
               if (iz.gt.1) then
                  CALL fdfd3d ('ddn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a4nuu1 = a4nuu
                  a4nuv1 = a4nuv
                  a4nuw1 = a4nuw
                  a4nvu1 = a4nvu
                  a4nvv1 = a4nvv
                  a4nvw1 = a4nvw
                  a4nwu1 = a4nwu
                  a4nwv1 = a4nwv
                  a4nww1 = a4nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('ddm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a4muu1 = a4muu
                     a4muv1 = a4muv
                     a4muw1 = a4muw
                     a4mvu1 = a4mvu
                     a4mvv1 = a4mvv
                     a4mvw1 = a4mvw
                     a4mwu1 = a4mwu
                     a4mwv1 = a4mwv
                     a4mww1 = a4mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('ddp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a4puu1 = a4puu
                     a4puv1 = a4puv
                     a4puw1 = a4puw
                     a4pvu1 = a4pvu
                     a4pvv1 = a4pvv
                     a4pvw1 = a4pvw
                     a4pwu1 = a4pwu
                     a4pwv1 = a4pwv
                     a4pww1 = a4pww
                  endif
               endif
               if (iz.lt.nz) then
                  CALL fdfd3d ('ffn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a6nuu1 = a6nuu
                  a6nuv1 = a6nuv
                  a6nuw1 = a6nuw
                  a6nvu1 = a6nvu
                  a6nvv1 = a6nvv
                  a6nvw1 = a6nvw
                  a6nwu1 = a6nwu
                  a6nwv1 = a6nwv
                  a6nww1 = a6nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('ffm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a6muu1 = a6muu
                     a6muv1 = a6muv
                     a6muw1 = a6muw
                     a6mvu1 = a6mvu
                     a6mvv1 = a6mvv
                     a6mvw1 = a6mvw
                     a6mwu1 = a6mwu
                     a6mwv1 = a6mwv
                     a6mww1 = a6mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('ffp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a6puu1 = a6puu
                     a6puv1 = a6puv
                     a6puw1 = a6puw
                     a6pvu1 = a6pvu
                     a6pvv1 = a6pvv
                     a6pvw1 = a6pvw
                     a6pwu1 = a6pwu
                     a6pwv1 = a6pwv
                     a6pww1 = a6pww
                  endif
               endif
               if (ix.lt.nx.and.iz.gt.1) then
                  CALL fdfd3d ('cdn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a7nuu1 = a7nuu
                  a7nuv1 = a7nuv
                  a7nuw1 = a7nuw
                  a7nvu1 = a7nvu
                  a7nvv1 = a7nvv
                  a7nvw1 = a7nvw
                  a7nwu1 = a7nwu
                  a7nwv1 = a7nwv
                  a7nww1 = a7nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('cdm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a7muu1 = a7muu
                     a7muv1 = a7muv
                     a7muw1 = a7muw
                     a7mvu1 = a7mvu
                     a7mvv1 = a7mvv
                     a7mvw1 = a7mvw
                     a7mwu1 = a7mwu
                     a7mwv1 = a7mwv
                     a7mww1 = a7mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('cdp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a7puu1 = a7puu
                     a7puv1 = a7puv
                     a7puw1 = a7puw
                     a7pvu1 = a7pvu
                     a7pvv1 = a7pvv
                     a7pvw1 = a7pvw
                     a7pwu1 = a7pwu
                     a7pwv1 = a7pwv
                     a7pww1 = a7pww
                  endif
               endif
               if (ix.lt.nx) then
                  CALL fdfd3d ('ccn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a8nuu1 = a8nuu
                  a8nuv1 = a8nuv
                  a8nuw1 = a8nuw
                  a8nvu1 = a8nvu
                  a8nvv1 = a8nvv
                  a8nvw1 = a8nvw
                  a8nwu1 = a8nwu
                  a8nwv1 = a8nwv
                  a8nww1 = a8nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('ccm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a8muu1 = a8muu
                     a8muv1 = a8muv
                     a8muw1 = a8muw
                     a8mvu1 = a8mvu
                     a8mvv1 = a8mvv
                     a8mvw1 = a8mvw
                     a8mwu1 = a8mwu
                     a8mwv1 = a8mwv
                     a8mww1 = a8mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('ccp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a8puu1 = a8puu
                     a8puv1 = a8puv
                     a8puw1 = a8puw
                     a8pvu1 = a8pvu
                     a8pvv1 = a8pvv
                     a8pvw1 = a8pvw
                     a8pwu1 = a8pwu
                     a8pwv1 = a8pwv
                     a8pww1 = a8pww
                  endif
               endif
               if (ix.lt.nx.and.iz.lt.nz) then
                  CALL fdfd3d ('cfn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a9nuu1 = a9nuu
                  a9nuv1 = a9nuv
                  a9nuw1 = a9nuw
                  a9nvu1 = a9nvu
                  a9nvv1 = a9nvv
                  a9nvw1 = a9nvw
                  a9nwu1 = a9nwu
                  a9nwv1 = a9nwv
                  a9nww1 = a9nww
                  if (iy.gt.1) then
                     CALL fdfd3d ('cfm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a9muu1 = a9muu
                     a9muv1 = a9muv
                     a9muw1 = a9muw
                     a9mvu1 = a9mvu
                     a9mvv1 = a9mvv
                     a9mvw1 = a9mvw
                     a9mwu1 = a9mwu
                     a9mwv1 = a9mwv
                     a9mww1 = a9mww
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('cfp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a9puu1 = a9puu
                     a9puv1 = a9puv
                     a9puw1 = a9puw
                     a9pvu1 = a9pvu
                     a9pvv1 = a9pvv
                     a9pvw1 = a9pvw
                     a9pwu1 = a9pwu
                     a9pwv1 = a9pwv
                     a9pww1 = a9pww
                  endif
               endif

c----background  model
               CALL setcoefs3d (iz, ix, iy, 1)
          
               CALL fdfd3d ('ben', iz, ix, iy, ierr)
               if (ierr.ne.0) goto 900
               a5nuu = a5nuu - a5nuu1
               a5nuv = a5nuv - a5nuv1
               a5nuw = a5nuw - a5nuw1
               a5nvu = a5nvu - a5nvu1
               a5nvv = a5nvv - a5nvv1
               a5nvw = a5nvw - a5nvw1
               a5nwu = a5nwu - a5nwu1
               a5nwv = a5nwv - a5nwv1
               a5nww = a5nww - a5nww1
               us = us + a5nuu*utbkr(iz,ix,iy)
               us = us + a5nuv*vtbkr(iz,ix,iy)
               us = us + a5nuw*wtbkr(iz,ix,iy)
               vs = vs + a5nvu*utbkr(iz,ix,iy)
               vs = vs + a5nvv*vtbkr(iz,ix,iy)
               vs = vs + a5nvw*wtbkr(iz,ix,iy)
               ws = ws + a5nwu*utbkr(iz,ix,iy)
               ws = ws + a5nwv*vtbkr(iz,ix,iy)
               ws = ws + a5nww*wtbkr(iz,ix,iy)
               if (iy.gt.1) then
                  CALL fdfd3d ('bem', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a5muu = a5muu - a5muu1
                  a5muv = a5muv - a5muv1
                  a5muw = a5muw - a5muw1
                  a5mvu = a5mvu - a5mvu1
                  a5mvv = a5mvv - a5mvv1
                  a5mvw = a5mvw - a5mvw1
                  a5mwu = a5mwu - a5mwu1
                  a5mwv = a5mwv - a5mwv1
                  a5mww = a5mww - a5mww1
                  us = us + a5muu*utbkr(iz,ix,iy-1)
                  us = us + a5muv*vtbkr(iz,ix,iy-1)
                  us = us + a5muw*wtbkr(iz,ix,iy-1)
                  vs = vs + a5mvu*utbkr(iz,ix,iy-1)
                  vs = vs + a5mvv*vtbkr(iz,ix,iy-1)
                  vs = vs + a5mvw*wtbkr(iz,ix,iy-1)
                  ws = ws + a5mwu*utbkr(iz,ix,iy-1)
                  ws = ws + a5mwv*vtbkr(iz,ix,iy-1)
                  ws = ws + a5mww*wtbkr(iz,ix,iy-1)
               endif
               if (iy.lt.ny) then
                  CALL fdfd3d ('bep', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a5puu = a5puu - a5puu1
                  a5puv = a5puv - a5puv1
                  a5puw = a5puw - a5puw1
                  a5pvu = a5pvu - a5pvu1
                  a5pvv = a5pvv - a5pvv1
                  a5pvw = a5pvw - a5pvw1
                  a5pwu = a5pwu - a5pwu1
                  a5pwv = a5pwv - a5pwv1
                  a5pww = a5pww - a5pww1
                  us = us + a5puu*utbkr(iz,ix,iy+1)
                  us = us + a5puv*vtbkr(iz,ix,iy+1)
                  us = us + a5puw*wtbkr(iz,ix,iy+1)
                  vs = vs + a5pvu*utbkr(iz,ix,iy+1)
                  vs = vs + a5pvv*vtbkr(iz,ix,iy+1)
                  vs = vs + a5pvw*wtbkr(iz,ix,iy+1)
                  ws = ws + a5pwu*utbkr(iz,ix,iy+1)
                  ws = ws + a5pwv*vtbkr(iz,ix,iy+1)
                  ws = ws + a5pww*wtbkr(iz,ix,iy+1)
               endif
               if (ix.gt.1.and.iz.gt.1) then
                  CALL fdfd3d ('adn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a1nuu = a1nuu - a1nuu1
                  a1nuv = a1nuv - a1nuv1
                  a1nuw = a1nuw - a1nuw1
                  a1nvu = a1nvu - a1nvu1
                  a1nvv = a1nvv - a1nvv1
                  a1nvw = a1nvw - a1nvw1
                  a1nwu = a1nwu - a1nwu1
                  a1nwv = a1nwv - a1nwv1
                  a1nww = a1nww - a1nww1
                  us = us + a1nuu*utbkr(iz-1,ix-1,iy)
                  us = us + a1nuv*vtbkr(iz-1,ix-1,iy)
                  us = us + a1nuw*wtbkr(iz-1,ix-1,iy)
                  vs = vs + a1nvu*utbkr(iz-1,ix-1,iy)
                  vs = vs + a1nvv*vtbkr(iz-1,ix-1,iy)
                  vs = vs + a1nvw*wtbkr(iz-1,ix-1,iy)
                  ws = ws + a1nwu*utbkr(iz-1,ix-1,iy)
                  ws = ws + a1nwv*vtbkr(iz-1,ix-1,iy)
                  ws = ws + a1nww*wtbkr(iz-1,ix-1,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('adm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a1muu = a1muu - a1muu1
                     a1muv = a1muv - a1muv1
                     a1muw = a1muw - a1muw1
                     a1mvu = a1mvu - a1mvu1
                     a1mvv = a1mvv - a1mvv1
                     a1mvw = a1mvw - a1mvw1
                     a1mwu = a1mwu - a1mwu1
                     a1mwv = a1mwv - a1mwv1
                     a1mww = a1mww - a1mww1
                     us = us + a1muu*utbkr(iz-1,ix-1,iy-1)
                     us = us + a1muv*vtbkr(iz-1,ix-1,iy-1)
                     us = us + a1muw*wtbkr(iz-1,ix-1,iy-1)
                     vs = vs + a1mvu*utbkr(iz-1,ix-1,iy-1)
                     vs = vs + a1mvv*vtbkr(iz-1,ix-1,iy-1)
                     vs = vs + a1mvw*wtbkr(iz-1,ix-1,iy-1)
                     ws = ws + a1mwu*utbkr(iz-1,ix-1,iy-1)
                     ws = ws + a1mwv*vtbkr(iz-1,ix-1,iy-1)
                     ws = ws + a1mww*wtbkr(iz-1,ix-1,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('adp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a1puu = a1puu - a1puu1
                     a1puv = a1puv - a1puv1
                     a1puw = a1puw - a1puw1
                     a1pvu = a1pvu - a1pvu1
                     a1pvv = a1pvv - a1pvv1
                     a1pvw = a1pvw - a1pvw1
                     a1pwu = a1pwu - a1pwu1
                     a1pwv = a1pwv - a1pwv1
                     a1pww = a1pww - a1pww1
                     us = us + a1puu*utbkr(iz-1,ix-1,iy+1)
                     us = us + a1puv*vtbkr(iz-1,ix-1,iy+1)
                     us = us + a1puw*wtbkr(iz-1,ix-1,iy+1)
                     vs = vs + a1pvu*utbkr(iz-1,ix-1,iy+1)
                     vs = vs + a1pvv*vtbkr(iz-1,ix-1,iy+1)
                     vs = vs + a1pvw*wtbkr(iz-1,ix-1,iy+1)
                     ws = ws + a1pwu*utbkr(iz-1,ix-1,iy+1)
                     ws = ws + a1pwv*vtbkr(iz-1,ix-1,iy+1)
                     ws = ws + a1pww*wtbkr(iz-1,ix-1,iy+1)
                  endif
               endif
               if (ix.gt.1) then
                  CALL fdfd3d ('aan', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a2nuu = a2nuu - a2nuu1
                  a2nuv = a2nuv - a2nuv1
                  a2nuw = a2nuw - a2nuw1
                  a2nvu = a2nvu - a2nvu1
                  a2nvv = a2nvv - a2nvv1
                  a2nvw = a2nvw - a2nvw1
                  a2nwu = a2nwu - a2nwu1
                  a2nwv = a2nwv - a2nwv1
                  a2nww = a2nww - a2nww1
                  us = us + a2nuu*utbkr(iz,ix-1,iy)
                  us = us + a2nuv*vtbkr(iz,ix-1,iy)
                  us = us + a2nuw*wtbkr(iz,ix-1,iy)
                  vs = vs + a2nvu*utbkr(iz,ix-1,iy)
                  vs = vs + a2nvv*vtbkr(iz,ix-1,iy)
                  vs = vs + a2nvw*wtbkr(iz,ix-1,iy)
                  ws = ws + a2nwu*utbkr(iz,ix-1,iy)
                  ws = ws + a2nwv*vtbkr(iz,ix-1,iy)
                  ws = ws + a2nww*wtbkr(iz,ix-1,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('aam', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a2muu = a2muu - a2muu1
                     a2muv = a2muv - a2muv1
                     a2muw = a2muw - a2muw1
                     a2mvu = a2mvu - a2mvu1
                     a2mvv = a2mvv - a2mvv1
                     a2mvw = a2mvw - a2mvw1
                     a2mwu = a2mwu - a2mwu1
                     a2mwv = a2mwv - a2mwv1
                     a2mww = a2mww - a2mww1
                     us = us + a2muu*utbkr(iz,ix-1,iy-1)
                     us = us + a2muv*vtbkr(iz,ix-1,iy-1)
                     us = us + a2muw*wtbkr(iz,ix-1,iy-1)
                     vs = vs + a2mvu*utbkr(iz,ix-1,iy-1)
                     vs = vs + a2mvv*vtbkr(iz,ix-1,iy-1)
                     vs = vs + a2mvw*wtbkr(iz,ix-1,iy-1)
                     ws = ws + a2mwu*utbkr(iz,ix-1,iy-1)
                     ws = ws + a2mwv*vtbkr(iz,ix-1,iy-1)
                     ws = ws + a2mww*wtbkr(iz,ix-1,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('aap', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a2puu = a2puu - a2puu1
                     a2puv = a2puv - a2puv1
                     a2puw = a2puw - a2puw1
                     a2pvu = a2pvu - a2pvu1
                     a2pvv = a2pvv - a2pvv1
                     a2pvw = a2pvw - a2pvw1
                     a2pwu = a2pwu - a2pwu1
                     a2pwv = a2pwv - a2pwv1
                     a2pww = a2pww - a2pww1
                     us = us + a2puu*utbkr(iz,ix-1,iy+1)
                     us = us + a2puv*vtbkr(iz,ix-1,iy+1)
                     us = us + a2puw*wtbkr(iz,ix-1,iy+1)
                     vs = vs + a2pvu*utbkr(iz,ix-1,iy+1)
                     vs = vs + a2pvv*vtbkr(iz,ix-1,iy+1)
                     vs = vs + a2pvw*wtbkr(iz,ix-1,iy+1)
                     ws = ws + a2pwu*utbkr(iz,ix-1,iy+1)
                     ws = ws + a2pwv*vtbkr(iz,ix-1,iy+1)
                     ws = ws + a2pww*wtbkr(iz,ix-1,iy+1)
                  endif
               endif
               if (ix.gt.1.and.iz.lt.nz) then
                  CALL fdfd3d ('afn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a3nuu = a3nuu - a3nuu1
                  a3nuv = a3nuv - a3nuv1
                  a3nuw = a3nuw - a3nuw1
                  a3nvu = a3nvu - a3nvu1
                  a3nvv = a3nvv - a3nvv1
                  a3nvw = a3nvw - a3nvw1
                  a3nwu = a3nwu - a3nwu1
                  a3nwv = a3nwv - a3nwv1
                  a3nww = a3nww - a3nww1
                  us = us + a3nuu*utbkr(iz+1,ix-1,iy)
                  us = us + a3nuv*vtbkr(iz+1,ix-1,iy)
                  us = us + a3nuw*wtbkr(iz+1,ix-1,iy)
                  vs = vs + a3nvu*utbkr(iz+1,ix-1,iy)
                  vs = vs + a3nvv*vtbkr(iz+1,ix-1,iy)
                  vs = vs + a3nvw*wtbkr(iz+1,ix-1,iy)
                  ws = ws + a3nwu*utbkr(iz+1,ix-1,iy)
                  ws = ws + a3nwv*vtbkr(iz+1,ix-1,iy)
                  ws = ws + a3nww*wtbkr(iz+1,ix-1,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('afm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a3muu = a3muu - a3muu1
                     a3muv = a3muv - a3muv1
                     a3muw = a3muw - a3muw1
                     a3mvu = a3mvu - a3mvu1
                     a3mvv = a3mvv - a3mvv1
                     a3mvw = a3mvw - a3mvw1
                     a3mwu = a3mwu - a3mwu1
                     a3mwv = a3mwv - a3mwv1
                     a3mww = a3mww - a3mww1
                     us = us + a3muu*utbkr(iz+1,ix-1,iy-1)
                     us = us + a3muv*vtbkr(iz+1,ix-1,iy-1)
                     us = us + a3muw*wtbkr(iz+1,ix-1,iy-1)
                     vs = vs + a3mvu*utbkr(iz+1,ix-1,iy-1)
                     vs = vs + a3mvv*vtbkr(iz+1,ix-1,iy-1)
                     vs = vs + a3mvw*wtbkr(iz+1,ix-1,iy-1)
                     ws = ws + a3mwu*utbkr(iz+1,ix-1,iy-1)
                     ws = ws + a3mwv*vtbkr(iz+1,ix-1,iy-1)
                     ws = ws + a3mww*wtbkr(iz+1,ix-1,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('afp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a3puu = a3puu - a3puu1
                     a3puv = a3puv - a3puv1
                     a3puw = a3puw - a3puw1
                     a3pvu = a3pvu - a3pvu1
                     a3pvv = a3pvv - a3pvv1
                     a3pvw = a3pvw - a3pvw1
                     a3pwu = a3pwu - a3pwu1
                     a3pwv = a3pwv - a3pwv1
                     a3pww = a3pww - a3pww1
                     us = us + a3puu*utbkr(iz+1,ix-1,iy+1)
                     us = us + a3puv*vtbkr(iz+1,ix-1,iy+1)
                     us = us + a3puw*wtbkr(iz+1,ix-1,iy+1)
                     vs = vs + a3pvu*utbkr(iz+1,ix-1,iy+1)
                     vs = vs + a3pvv*vtbkr(iz+1,ix-1,iy+1)
                     vs = vs + a3pvw*wtbkr(iz+1,ix-1,iy+1)
                     ws = ws + a3pwu*utbkr(iz+1,ix-1,iy+1)
                     ws = ws + a3pwv*vtbkr(iz+1,ix-1,iy+1)
                     ws = ws + a3pww*wtbkr(iz+1,ix-1,iy+1)
                  endif
               endif
               if (iz.gt.1) then
                  CALL fdfd3d ('ddn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a4nuu = a4nuu - a4nuu1
                  a4nuv = a4nuv - a4nuv1
                  a4nuw = a4nuw - a4nuw1
                  a4nvu = a4nvu - a4nvu1
                  a4nvv = a4nvv - a4nvv1
                  a4nvw = a4nvw - a4nvw1
                  a4nwu = a4nwu - a4nwu1
                  a4nwv = a4nwv - a4nwv1
                  a4nww = a4nww - a4nww1
                  us = us + a4nuu*utbkr(iz-1,ix,iy)
                  us = us + a4nuv*vtbkr(iz-1,ix,iy)
                  us = us + a4nuw*wtbkr(iz-1,ix,iy)
                  vs = vs + a4nvu*utbkr(iz-1,ix,iy)
                  vs = vs + a4nvv*vtbkr(iz-1,ix,iy)
                  vs = vs + a4nvw*wtbkr(iz-1,ix,iy)
                  ws = ws + a4nwu*utbkr(iz-1,ix,iy)
                  ws = ws + a4nwv*vtbkr(iz-1,ix,iy)
                  ws = ws + a4nww*wtbkr(iz-1,ix,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('ddm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a4muu = a4muu - a4muu1
                     a4muv = a4muv - a4muv1
                     a4muw = a4muw - a4muw1
                     a4mvu = a4mvu - a4mvu1
                     a4mvv = a4mvv - a4mvv1
                     a4mvw = a4mvw - a4mvw1
                     a4mwu = a4mwu - a4mwu1
                     a4mwv = a4mwv - a4mwv1
                     a4mww = a4mww - a4mww1
                     us = us + a4muu*utbkr(iz-1,ix,iy-1)
                     us = us + a4muv*vtbkr(iz-1,ix,iy-1)
                     us = us + a4muw*wtbkr(iz-1,ix,iy-1)
                     vs = vs + a4mvu*utbkr(iz-1,ix,iy-1)
                     vs = vs + a4mvv*vtbkr(iz-1,ix,iy-1)
                     vs = vs + a4mvw*wtbkr(iz-1,ix,iy-1)
                     ws = ws + a4mwu*utbkr(iz-1,ix,iy-1)
                     ws = ws + a4mwv*vtbkr(iz-1,ix,iy-1)
                     ws = ws + a4mww*wtbkr(iz-1,ix,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('ddp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a4puu = a4puu - a4puu1
                     a4puv = a4puv - a4puv1
                     a4puw = a4puw - a4puw1
                     a4pvu = a4pvu - a4pvu1
                     a4pvv = a4pvv - a4pvv1
                     a4pvw = a4pvw - a4pvw1
                     a4pwu = a4pwu - a4pwu1
                     a4pwv = a4pwv - a4pwv1
                     a4pww = a4pww - a4pww1
                     us = us + a4puu*utbkr(iz-1,ix,iy+1)
                     us = us + a4puv*vtbkr(iz-1,ix,iy+1)
                     us = us + a4puw*wtbkr(iz-1,ix,iy+1)
                     vs = vs + a4pvu*utbkr(iz-1,ix,iy+1)
                     vs = vs + a4pvv*vtbkr(iz-1,ix,iy+1)
                     vs = vs + a4pvw*wtbkr(iz-1,ix,iy+1)
                     ws = ws + a4pwu*utbkr(iz-1,ix,iy+1)
                     ws = ws + a4pwv*vtbkr(iz-1,ix,iy+1)
                     ws = ws + a4pww*wtbkr(iz-1,ix,iy+1)
                  endif
               endif
               if (iz.lt.nz) then
                  CALL fdfd3d ('ffn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a6nuu = a6nuu - a6nuu1
                  a6nuv = a6nuv - a6nuv1
                  a6nuw = a6nuw - a6nuw1
                  a6nvu = a6nvu - a6nvu1
                  a6nvv = a6nvv - a6nvv1
                  a6nvw = a6nvw - a6nvw1
                  a6nwu = a6nwu - a6nwu1
                  a6nwv = a6nwv - a6nwv1
                  a6nww = a6nww - a6nww1
                  us = us + a6nuu*utbkr(iz+1,ix,iy)
                  us = us + a6nuv*vtbkr(iz+1,ix,iy)
                  us = us + a6nuw*wtbkr(iz+1,ix,iy)
                  vs = vs + a6nvu*utbkr(iz+1,ix,iy)
                  vs = vs + a6nvv*vtbkr(iz+1,ix,iy)
                  vs = vs + a6nvw*wtbkr(iz+1,ix,iy)
                  ws = ws + a6nwu*utbkr(iz+1,ix,iy)
                  ws = ws + a6nwv*vtbkr(iz+1,ix,iy)
                  ws = ws + a6nww*wtbkr(iz+1,ix,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('ffm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a6muu = a6muu - a6muu1
                     a6muv = a6muv - a6muv1
                     a6muw = a6muw - a6muw1
                     a6mvu = a6mvu - a6mvu1
                     a6mvv = a6mvv - a6mvv1
                     a6mvw = a6mvw - a6mvw1
                     a6mwu = a6mwu - a6mwu1
                     a6mwv = a6mwv - a6mwv1
                     a6mww = a6mww - a6mww1
                     us = us + a6muu*utbkr(iz+1,ix,iy-1)
                     us = us + a6muv*vtbkr(iz+1,ix,iy-1)
                     us = us + a6muw*wtbkr(iz+1,ix,iy-1)
                     vs = vs + a6mvu*utbkr(iz+1,ix,iy-1)
                     vs = vs + a6mvv*vtbkr(iz+1,ix,iy-1)
                     vs = vs + a6mvw*wtbkr(iz+1,ix,iy-1)
                     ws = ws + a6mwu*utbkr(iz+1,ix,iy-1)
                     ws = ws + a6mwv*vtbkr(iz+1,ix,iy-1)
                     ws = ws + a6mww*wtbkr(iz+1,ix,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('ffp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a6puu = a6puu - a6puu1
                     a6puv = a6puv - a6puv1
                     a6puw = a6puw - a6puw1
                     a6pvu = a6pvu - a6pvu1
                     a6pvv = a6pvv - a6pvv1
                     a6pvw = a6pvw - a6pvw1
                     a6pwu = a6pwu - a6pwu1
                     a6pwv = a6pwv - a6pwv1
                     a6pww = a6pww - a6pww1
                     us = us + a6puu*utbkr(iz+1,ix,iy+1)
                     us = us + a6puv*vtbkr(iz+1,ix,iy+1)
                     us = us + a6puw*wtbkr(iz+1,ix,iy+1)
                     vs = vs + a6pvu*utbkr(iz+1,ix,iy+1)
                     vs = vs + a6pvv*vtbkr(iz+1,ix,iy+1)
                     vs = vs + a6pvw*wtbkr(iz+1,ix,iy+1)
                     ws = ws + a6pwu*utbkr(iz+1,ix,iy+1)
                     ws = ws + a6pwv*vtbkr(iz+1,ix,iy+1)
                     ws = ws + a6pww*wtbkr(iz+1,ix,iy+1)
                  endif
               endif
               if (ix.lt.nx.and.iz.gt.1) then
                  CALL fdfd3d ('cdn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a7nuu = a7nuu - a7nuu1
                  a7nuv = a7nuv - a7nuv1
                  a7nuw = a7nuw - a7nuw1
                  a7nvu = a7nvu - a7nvu1
                  a7nvv = a7nvv - a7nvv1
                  a7nvw = a7nvw - a7nvw1
                  a7nwu = a7nwu - a7nwu1
                  a7nwv = a7nwv - a7nwv1
                  a7nww = a7nww - a7nww1
                  us = us + a7nuu*utbkr(iz-1,ix+1,iy)
                  us = us + a7nuv*vtbkr(iz-1,ix+1,iy)
                  us = us + a7nuw*wtbkr(iz-1,ix+1,iy)
                  vs = vs + a7nvu*utbkr(iz-1,ix+1,iy)
                  vs = vs + a7nvv*vtbkr(iz-1,ix+1,iy)
                  vs = vs + a7nvw*wtbkr(iz-1,ix+1,iy)
                  ws = ws + a7nwu*utbkr(iz-1,ix+1,iy)
                  ws = ws + a7nwv*vtbkr(iz-1,ix+1,iy)
                  ws = ws + a7nww*wtbkr(iz-1,ix+1,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('cdm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a7muu = a7muu - a7muu1
                     a7muv = a7muv - a7muv1
                     a7muw = a7muw - a7muw1
                     a7mvu = a7mvu - a7mvu1
                     a7mvv = a7mvv - a7mvv1
                     a7mvw = a7mvw - a7mvw1
                     a7mwu = a7mwu - a7mwu1
                     a7mwv = a7mwv - a7mwv1
                     a7mww = a7mww - a7mww1
                     us = us + a7muu*utbkr(iz-1,ix+1,iy-1)
                     us = us + a7muv*vtbkr(iz-1,ix+1,iy-1)
                     us = us + a7muw*wtbkr(iz-1,ix+1,iy-1)
                     vs = vs + a7mvu*utbkr(iz-1,ix+1,iy-1)
                     vs = vs + a7mvv*vtbkr(iz-1,ix+1,iy-1)
                     vs = vs + a7mvw*wtbkr(iz-1,ix+1,iy-1)
                     ws = ws + a7mwu*utbkr(iz-1,ix+1,iy-1)
                     ws = ws + a7mwv*vtbkr(iz-1,ix+1,iy-1)
                     ws = ws + a7mww*wtbkr(iz-1,ix+1,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('cdp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a7puu = a7puu - a7puu1
                     a7puv = a7puv - a7puv1
                     a7puw = a7puw - a7puw1
                     a7pvu = a7pvu - a7pvu1
                     a7pvv = a7pvv - a7pvv1
                     a7pvw = a7pvw - a7pvw1
                     a7pwu = a7pwu - a7pwu1
                     a7pwv = a7pwv - a7pwv1
                     a7pww = a7pww - a7pww1
                     us = us + a7puu*utbkr(iz-1,ix+1,iy+1)
                     us = us + a7puv*vtbkr(iz-1,ix+1,iy+1)
                     us = us + a7puw*wtbkr(iz-1,ix+1,iy+1)
                     vs = vs + a7pvu*utbkr(iz-1,ix+1,iy+1)
                     vs = vs + a7pvv*vtbkr(iz-1,ix+1,iy+1)
                     vs = vs + a7pvw*wtbkr(iz-1,ix+1,iy+1)
                     ws = ws + a7pwu*utbkr(iz-1,ix+1,iy+1)
                     ws = ws + a7pwv*vtbkr(iz-1,ix+1,iy+1)
                     ws = ws + a7pww*wtbkr(iz-1,ix+1,iy+1)
                  endif
               endif
               if (ix.lt.nx) then
                  CALL fdfd3d ('ccn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a8nuu = a8nuu - a8nuu1
                  a8nuv = a8nuv - a8nuv1
                  a8nuw = a8nuw - a8nuw1
                  a8nvu = a8nvu - a8nvu1
                  a8nvv = a8nvv - a8nvv1
                  a8nvw = a8nvw - a8nvw1
                  a8nwu = a8nwu - a8nwu1
                  a8nwv = a8nwv - a8nwv1
                  a8nww = a8nww - a8nww1
                  us = us + a8nuu*utbkr(iz,ix+1,iy)
                  us = us + a8nuv*vtbkr(iz,ix+1,iy)
                  us = us + a8nuw*wtbkr(iz,ix+1,iy)
                  vs = vs + a8nvu*utbkr(iz,ix+1,iy)
                  vs = vs + a8nvv*vtbkr(iz,ix+1,iy)
                  vs = vs + a8nvw*wtbkr(iz,ix+1,iy)
                  ws = ws + a8nwu*utbkr(iz,ix+1,iy)
                  ws = ws + a8nwv*vtbkr(iz,ix+1,iy)
                  ws = ws + a8nww*wtbkr(iz,ix+1,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('ccm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a8muu = a8muu - a8muu1
                     a8muv = a8muv - a8muv1
                     a8muw = a8muw - a8muw1
                     a8mvu = a8mvu - a8mvu1
                     a8mvv = a8mvv - a8mvv1
                     a8mvw = a8mvw - a8mvw1
                     a8mwu = a8mwu - a8mwu1
                     a8mwv = a8mwv - a8mwv1
                     a8mww = a8mww - a8mww1
                     us = us + a8muu*utbkr(iz,ix+1,iy-1)
                     us = us + a8muv*vtbkr(iz,ix+1,iy-1)
                     us = us + a8muw*wtbkr(iz,ix+1,iy-1)
                     vs = vs + a8mvu*utbkr(iz,ix+1,iy-1)
                     vs = vs + a8mvv*vtbkr(iz,ix+1,iy-1)
                     vs = vs + a8mvw*wtbkr(iz,ix+1,iy-1)
                     ws = ws + a8mwu*utbkr(iz,ix+1,iy-1)
                     ws = ws + a8mwv*vtbkr(iz,ix+1,iy-1)
                     ws = ws + a8mww*wtbkr(iz,ix+1,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('ccp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a8puu = a8puu - a8puu1
                     a8puv = a8puv - a8puv1
                     a8puw = a8puw - a8puw1
                     a8pvu = a8pvu - a8pvu1
                     a8pvv = a8pvv - a8pvv1
                     a8pvw = a8pvw - a8pvw1
                     a8pwu = a8pwu - a8pwu1
                     a8pwv = a8pwv - a8pwv1
                     a8pww = a8pww - a8pww1
                     us = us + a8puu*utbkr(iz,ix+1,iy+1)
                     us = us + a8puv*vtbkr(iz,ix+1,iy+1)
                     us = us + a8puw*wtbkr(iz,ix+1,iy+1)
                     vs = vs + a8pvu*utbkr(iz,ix+1,iy+1)
                     vs = vs + a8pvv*vtbkr(iz,ix+1,iy+1)
                     vs = vs + a8pvw*wtbkr(iz,ix+1,iy+1)
                     ws = ws + a8pwu*utbkr(iz,ix+1,iy+1)
                     ws = ws + a8pwv*vtbkr(iz,ix+1,iy+1)
                     ws = ws + a8pww*wtbkr(iz,ix+1,iy+1)
                  endif
               endif
               if (ix.lt.nx.and.iz.lt.nz) then
                  CALL fdfd3d ('cfn', iz, ix, iy, ierr)
                  if (ierr.ne.0) goto 900
                  a9nuu = a9nuu - a9nuu1
                  a9nuv = a9nuv - a9nuv1
                  a9nuw = a9nuw - a9nuw1
                  a9nvu = a9nvu - a9nvu1
                  a9nvv = a9nvv - a9nvv1
                  a9nvw = a9nvw - a9nvw1
                  a9nwu = a9nwu - a9nwu1
                  a9nwv = a9nwv - a9nwv1
                  a9nww = a9nww - a9nww1
                  us = us + a9nuu*utbkr(iz+1,ix+1,iy)
                  us = us + a9nuv*vtbkr(iz+1,ix+1,iy)
                  us = us + a9nuw*wtbkr(iz+1,ix+1,iy)
                  vs = vs + a9nvu*utbkr(iz+1,ix+1,iy)
                  vs = vs + a9nvv*vtbkr(iz+1,ix+1,iy)
                  vs = vs + a9nvw*wtbkr(iz+1,ix+1,iy)
                  ws = ws + a9nwu*utbkr(iz+1,ix+1,iy)
                  ws = ws + a9nwv*vtbkr(iz+1,ix+1,iy)
                  ws = ws + a9nww*wtbkr(iz+1,ix+1,iy)
                  if (iy.gt.1) then
                     CALL fdfd3d ('cfm', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a9muu = a9muu - a9muu1
                     a9muv = a9muv - a9muv1
                     a9muw = a9muw - a9muw1
                     a9mvu = a9mvu - a9mvu1
                     a9mvv = a9mvv - a9mvv1
                     a9mvw = a9mvw - a9mvw1
                     a9mwu = a9mwu - a9mwu1
                     a9mwv = a9mwv - a9mwv1
                     a9mww = a9mww - a9mww1
                     us = us + a9muu*utbkr(iz+1,ix+1,iy-1)
                     us = us + a9muv*vtbkr(iz+1,ix+1,iy-1)
                     us = us + a9muw*wtbkr(iz+1,ix+1,iy-1)
                     vs = vs + a9mvu*utbkr(iz+1,ix+1,iy-1)
                     vs = vs + a9mvv*vtbkr(iz+1,ix+1,iy-1)
                     vs = vs + a9mvw*wtbkr(iz+1,ix+1,iy-1)
                     ws = ws + a9mwu*utbkr(iz+1,ix+1,iy-1)
                     ws = ws + a9mwv*vtbkr(iz+1,ix+1,iy-1)
                     ws = ws + a9mww*wtbkr(iz+1,ix+1,iy-1)
                  endif
                  if (iy.lt.ny) then
                     CALL fdfd3d ('cfp', iz, ix, iy, ierr)
                     if (ierr.ne.0) goto 900
                     a9puu = a9puu - a9puu1
                     a9puv = a9puv - a9puv1
                     a9puw = a9puw - a9puw1
                     a9pvu = a9pvu - a9pvu1
                     a9pvv = a9pvv - a9pvv1
                     a9pvw = a9pvw - a9pvw1
                     a9pwu = a9pwu - a9pwu1
                     a9pwv = a9pwv - a9pwv1
                     a9pww = a9pww - a9pww1
                     us = us + a9puu*utbkr(iz+1,ix+1,iy+1)
                     us = us + a9puv*vtbkr(iz+1,ix+1,iy+1)
                     us = us + a9puw*wtbkr(iz+1,ix+1,iy+1)
                     vs = vs + a9pvu*utbkr(iz+1,ix+1,iy+1)
                     vs = vs + a9pvv*vtbkr(iz+1,ix+1,iy+1)
                     vs = vs + a9pvw*wtbkr(iz+1,ix+1,iy+1)
                     ws = ws + a9pwu*utbkr(iz+1,ix+1,iy+1)
                     ws = ws + a9pwv*vtbkr(iz+1,ix+1,iy+1)
                     ws = ws + a9pww*wtbkr(iz+1,ix+1,iy+1)
                  endif
               endif

c---this was a temporary test to be sure that homogeneous media did not produce any sources
c       if (myid.eq.master) then
c         if (us.ne.0.0.or.vs.ne.0.0.or.ws.ne.0.0) then
c           write(*,*) ' in hetfldsrc:  ', ix, iy, iz, us, vs, ws
c           stop
c         endif
c       endif
        
c---convolve the source time function and assign to the source array
c---note on sign: above we formed (so-s)uo, so the source term should be positive
               igrd = (iy-1)*nx*ny + (ix-1)*nz + iz
               indx = 3*igrd-2
               sfld(indx)   =  us*srcfn
               sfld(indx+1) =  vs*srcfn
               sfld(indx+2) =  ws*srcfn

            enddo
         enddo
      enddo

900   continue
      if (ierr.ne.0) then
         if (myid.eq.master) then
            write(*,*)      'Error in hetfldsrc calling fdfd3d'
            !write(lunlog,*) 'Error in hetfldsrc calling fdfd3d' 
         endif
      endif 
      return
      end
        
