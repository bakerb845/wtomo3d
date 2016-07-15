c***********************************************************************
c
      SUBROUTINE fixed3d(lbl)
      USE MAT_MODULE
      USE FD_ENUM_MODULE, ONLY : bem, ben, bep,
     ;                           adm, adn, adp,
     ;                           aam, aan, aap,
     ;                           afm, afn, afp,
     ;                           ddm, ddn, ddp,
     ;                           ffm, ffn, ffp,
     ;                           cdm, cdn, cdp,
     ;                           ccm, ccn, ccp,
     ;                           cfm, cfn, cfp
      IMPLICIT NONE
C
C     Input:
C     iz, ix  are the model indices for the point in question
C     lbl  is a two letter code to specify which of the 9 elements is 
C     desired
C
C     The condition here is simply (u,v,w) = 0, so set the diagonals 
c     to 1 and the remainder to 0.
c
C     Change 03/94 I. Stekl by adding a 9 point star approximation to
C              get 4 grid points per wavelength  with .5% accurarcy but 
C              condition is to get dx=dz
C 
C     Change 02/10 B. Baker removed type conversions 
C  
c parameters:
!       include 'dimension.inc'            
c arguments:
      INTEGER, INTENT(IN) :: lbl
c common variables:
!       include 'init.inc'
!       include 'common.inc'
!       include 'coeffs.inc'
c_______________________________________________________________________
c executable code:

c----a5
      if (lbl.eq.bem) then
         a5muv = CMPLX(0.0,0.0)
         a5muu = CMPLX(0.0,0.0)
         a5muw = CMPLX(0.0,0.0)
         a5mvu = CMPLX(0.0,0.0)
         a5mvv = CMPLX(0.0,0.0)
         a5mvw = CMPLX(0.0,0.0)
         a5mwu = CMPLX(0.0,0.0)
         a5mwv = CMPLX(0.0,0.0)
         a5mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.ben) then
         a5nuv = CMPLX(0.0,0.0)
         a5nuu = CMPLX(1.0,0.0)
         a5nuw = CMPLX(0.0,0.0)
         a5nvu = CMPLX(0.0,0.0)
         a5nvv = CMPLX(1.0,0.0)
         a5nvw = CMPLX(0.0,0.0)
         a5nwu = CMPLX(0.0,0.0)
         a5nwv = CMPLX(0.0,0.0)
         a5nww = CMPLX(1.0,0.0)

      elseif (lbl.eq.bep) then
         a5puv = CMPLX(0.0,0.0)
         a5puu = CMPLX(0.0,0.0)
         a5puw = CMPLX(0.0,0.0)
         a5pvu = CMPLX(0.0,0.0)
         a5pvv = CMPLX(0.0,0.0)
         a5pvw = CMPLX(0.0,0.0)
         a5pwu = CMPLX(0.0,0.0)
         a5pwv = CMPLX(0.0,0.0)
         a5pww = CMPLX(0.0,0.0)
c----a2
      elseif (lbl.eq.aam) then
         a2muv = CMPLX(0.0,0.0)
         a2muu = CMPLX(0.0,0.0)
         a2muw = CMPLX(0.0,0.0)
         a2mvu = CMPLX(0.0,0.0)
         a2mvv = CMPLX(0.0,0.0)
         a2mvw = CMPLX(0.0,0.0)
         a2mwu = CMPLX(0.0,0.0)
         a2mwv = CMPLX(0.0,0.0)
         a2mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.aan) then
         a2nuv = CMPLX(0.0,0.0)
         a2nuu = CMPLX(0.0,0.0)
         a2nuw = CMPLX(0.0,0.0)
         a2nvu = CMPLX(0.0,0.0)
         a2nvv = CMPLX(0.0,0.0)
         a2nvw = CMPLX(0.0,0.0)
         a2nwu = CMPLX(0.0,0.0)
         a2nwv = CMPLX(0.0,0.0)
         a2nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.aap) then
         a2puv = CMPLX(0.0,0.0)
         a2puu = CMPLX(0.0,0.0)
         a2puw = CMPLX(0.0,0.0)
         a2pvu = CMPLX(0.0,0.0)
         a2pvv = CMPLX(0.0,0.0)
         a2pvw = CMPLX(0.0,0.0)
         a2pwu = CMPLX(0.0,0.0)
         a2pwv = CMPLX(0.0,0.0)
         a2pww = CMPLX(0.0,0.0)
c----a8
      elseif (lbl.eq.ccm) then
         a8muv = CMPLX(0.0,0.0)
         a8muu = CMPLX(0.0,0.0)
         a8muw = CMPLX(0.0,0.0)
         a8mvu = CMPLX(0.0,0.0)
         a8mvv = CMPLX(0.0,0.0)
         a8mvw = CMPLX(0.0,0.0)
         a8mwu = CMPLX(0.0,0.0)
         a8mwv = CMPLX(0.0,0.0)
         a8mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.ccn) then
         a8nuv = CMPLX(0.0,0.0)
         a8nuu = CMPLX(0.0,0.0)
         a8nuw = CMPLX(0.0,0.0)
         a8nvu = CMPLX(0.0,0.0)
         a8nvv = CMPLX(0.0,0.0)
         a8nvw = CMPLX(0.0,0.0)
         a8nwu = CMPLX(0.0,0.0)
         a8nwv = CMPLX(0.0,0.0)
         a8nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.ccp) then
         a8puv = CMPLX(0.0,0.0)
         a8puu = CMPLX(0.0,0.0)
         a8puw = CMPLX(0.0,0.0)
         a8pvu = CMPLX(0.0,0.0)
         a8pvv = CMPLX(0.0,0.0)
         a8pvw = CMPLX(0.0,0.0)
         a8pwu = CMPLX(0.0,0.0)
         a8pwv = CMPLX(0.0,0.0)
         a8pww = CMPLX(0.0,0.0)
c----a4
      elseif (lbl.eq.ddm) then
         a4muv = CMPLX(0.0,0.0)
         a4muu = CMPLX(0.0,0.0)
         a4muw = CMPLX(0.0,0.0)
         a4mvu = CMPLX(0.0,0.0)
         a4mvv = CMPLX(0.0,0.0)
         a4mvw = CMPLX(0.0,0.0)
         a4mwu = CMPLX(0.0,0.0)
         a4mwv = CMPLX(0.0,0.0)
         a4mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.ddn) then
         a4nuv = CMPLX(0.0,0.0)
         a4nuu = CMPLX(0.0,0.0)
         a4nuw = CMPLX(0.0,0.0)
         a4nvu = CMPLX(0.0,0.0)
         a4nvv = CMPLX(0.0,0.0)
         a4nvw = CMPLX(0.0,0.0)
         a4nwu = CMPLX(0.0,0.0)
         a4nwv = CMPLX(0.0,0.0)
         a4nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.ddp) then
         a4puv = CMPLX(0.0,0.0)
         a4puu = CMPLX(0.0,0.0)
         a4puw = CMPLX(0.0,0.0)
         a4pvu = CMPLX(0.0,0.0)
         a4pvv = CMPLX(0.0,0.0)
         a4pvw = CMPLX(0.0,0.0)
         a4pwu = CMPLX(0.0,0.0)
         a4pwv = CMPLX(0.0,0.0)
         a4pww = CMPLX(0.0,0.0)
c----a6
      elseif (lbl.eq.ffm) then
         a6muv = CMPLX(0.0,0.0)
         a6muu = CMPLX(0.0,0.0)
         a6muw = CMPLX(0.0,0.0)
         a6mvu = CMPLX(0.0,0.0)
         a6mvv = CMPLX(0.0,0.0)
         a6mvw = CMPLX(0.0,0.0)
         a6mwu = CMPLX(0.0,0.0)
         a6mwv = CMPLX(0.0,0.0)
         a6mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.ffn) then
         a6nuv = CMPLX(0.0,0.0)
         a6nuu = CMPLX(0.0,0.0)
         a6nuw = CMPLX(0.0,0.0)
         a6nvu = CMPLX(0.0,0.0)
         a6nvv = CMPLX(0.0,0.0)
         a6nvw = CMPLX(0.0,0.0)
         a6nwu = CMPLX(0.0,0.0)
         a6nwv = CMPLX(0.0,0.0)
         a6nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.ffp) then
         a6puv = CMPLX(0.0,0.0)
         a6puu = CMPLX(0.0,0.0)
         a6puw = CMPLX(0.0,0.0)
         a6pvu = CMPLX(0.0,0.0)
         a6pvv = CMPLX(0.0,0.0)
         a6pvw = CMPLX(0.0,0.0)
         a6pwu = CMPLX(0.0,0.0)
         a6pwv = CMPLX(0.0,0.0)
         a6pww = CMPLX(0.0,0.0)
c----a1
      elseif (lbl.eq.adm) then
         a1muv = CMPLX(0.0,0.0)
         a1muu = CMPLX(0.0,0.0)
         a1muw = CMPLX(0.0,0.0)
         a1mvu = CMPLX(0.0,0.0)
         a1mvv = CMPLX(0.0,0.0)
         a1mvw = CMPLX(0.0,0.0)
         a1mwu = CMPLX(0.0,0.0)
         a1mwv = CMPLX(0.0,0.0)
         a1mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.adn) then
         a1nuv = CMPLX(0.0,0.0)
         a1nuu = CMPLX(0.0,0.0)
         a1nuw = CMPLX(0.0,0.0)
         a1nvu = CMPLX(0.0,0.0)
         a1nvv = CMPLX(0.0,0.0)
         a1nvw = CMPLX(0.0,0.0)
         a1nwu = CMPLX(0.0,0.0)
         a1nwv = CMPLX(0.0,0.0)
         a1nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.adp) then
         a1puv = CMPLX(0.0,0.0)
         a1puu = CMPLX(0.0,0.0)
         a1puw = CMPLX(0.0,0.0)
         a1pvu = CMPLX(0.0,0.0)
         a1pvv = CMPLX(0.0,0.0)
         a1pvw = CMPLX(0.0,0.0)
         a1pwu = CMPLX(0.0,0.0)
         a1pwv = CMPLX(0.0,0.0)
         a1pww = CMPLX(0.0,0.0)
c----a9
      elseif (lbl.eq.cfm) then
         a9muv = CMPLX(0.0,0.0)
         a9muu = CMPLX(0.0,0.0)
         a9muw = CMPLX(0.0,0.0)
         a9mvu = CMPLX(0.0,0.0)
         a9mvv = CMPLX(0.0,0.0)
         a9mvw = CMPLX(0.0,0.0)
         a9mwu = CMPLX(0.0,0.0)
         a9mwv = CMPLX(0.0,0.0)
         a9mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.cfn) then
         a9nuv = CMPLX(0.0,0.0)
         a9nuu = CMPLX(0.0,0.0)
         a9nuw = CMPLX(0.0,0.0)
         a9nvu = CMPLX(0.0,0.0)
         a9nvv = CMPLX(0.0,0.0)
         a9nvw = CMPLX(0.0,0.0)
         a9nwu = CMPLX(0.0,0.0)
         a9nwv = CMPLX(0.0,0.0)
         a9nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.cfp) then
         a9puv = CMPLX(0.0,0.0)
         a9puu = CMPLX(0.0,0.0)
         a9puw = CMPLX(0.0,0.0)
         a9pvu = CMPLX(0.0,0.0)
         a9pvv = CMPLX(0.0,0.0)
         a9pvw = CMPLX(0.0,0.0)
         a9pwu = CMPLX(0.0,0.0)
         a9pwv = CMPLX(0.0,0.0)
         a9pww = CMPLX(0.0,0.0)
c----a3
      elseif (lbl.eq.afm) then
         a3muv = CMPLX(0.0,0.0)
         a3muu = CMPLX(0.0,0.0)
         a3muw = CMPLX(0.0,0.0)
         a3mvu = CMPLX(0.0,0.0)
         a3mvv = CMPLX(0.0,0.0)
         a3mvw = CMPLX(0.0,0.0)
         a3mwu = CMPLX(0.0,0.0)
         a3mwv = CMPLX(0.0,0.0)
         a3mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.afn) then
         a3nuv = CMPLX(0.0,0.0)
         a3nuu = CMPLX(0.0,0.0)
         a3nuw = CMPLX(0.0,0.0)
         a3nvu = CMPLX(0.0,0.0)
         a3nvv = CMPLX(0.0,0.0)
         a3nvw = CMPLX(0.0,0.0)
         a3nwu = CMPLX(0.0,0.0)
         a3nwv = CMPLX(0.0,0.0)
         a3nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.afp) then
         a3puv = CMPLX(0.0,0.0)
         a3puu = CMPLX(0.0,0.0)
         a3puw = CMPLX(0.0,0.0)
         a3pvu = CMPLX(0.0,0.0)
         a3pvv = CMPLX(0.0,0.0)
         a3pvw = CMPLX(0.0,0.0)
         a3pwu = CMPLX(0.0,0.0)
         a3pwv = CMPLX(0.0,0.0)
         a3pww = CMPLX(0.0,0.0)
c----a7
      elseif (lbl.eq.cdm) then
         a7muv = CMPLX(0.0,0.0)
         a7muu = CMPLX(0.0,0.0)
         a7muw = CMPLX(0.0,0.0)
         a7mvu = CMPLX(0.0,0.0)
         a7mvv = CMPLX(0.0,0.0)
         a7mvw = CMPLX(0.0,0.0)
         a7mwu = CMPLX(0.0,0.0)
         a7mwv = CMPLX(0.0,0.0)
         a7mww = CMPLX(0.0,0.0)

      elseif (lbl.eq.cdn) then
         a7nuv = CMPLX(0.0,0.0)
         a7nuu = CMPLX(0.0,0.0)
         a7nuw = CMPLX(0.0,0.0)
         a7nvu = CMPLX(0.0,0.0)
         a7nvv = CMPLX(0.0,0.0)
         a7nvw = CMPLX(0.0,0.0)
         a7nwu = CMPLX(0.0,0.0)
         a7nwv = CMPLX(0.0,0.0)
         a7nww = CMPLX(0.0,0.0)

      elseif (lbl.eq.cdp) then
         a7puv = CMPLX(0.0,0.0)
         a7puu = CMPLX(0.0,0.0)
         a7puw = CMPLX(0.0,0.0)
         a7pvu = CMPLX(0.0,0.0)
         a7pvv = CMPLX(0.0,0.0)
         a7pvw = CMPLX(0.0,0.0)
         a7pwu = CMPLX(0.0,0.0)
         a7pwv = CMPLX(0.0,0.0)
         a7pww = CMPLX(0.0,0.0)

      else
         write(*,100) lbl
  100    format(' warning: fixed cannot identify coefficient ',i2,'!')

      endif

      return
      end

