      SUBROUTINE zeroas()
c 
c     Zero out each of the a's
      USE MAT_MODULE

      COMPLEX, PARAMETER :: zero = CMPLX(0.0, 0.0)

      a1muu = zero 
      a1muv = zero 
      a1muw = zero
      a1mvu = zero 
      a1mvv = zero
      a1mvw = zero
      a1mwu = zero
      a1mwv = zero
      a1mww = zero
      a2muu = zero
      a2muv = zero
      a2muw = zero
      a2mvu = zero
      a2mvv = zero
      a2mvw = zero
      a2mwu = zero
      a2mwv = zero
      a2mww = zero
      a3muu = zero
      a3muv = zero
      a3muw = zero
      a3mvu = zero
      a3mvv = zero
      a3mvw = zero
      a3mwu = zero
      a3mwv = zero
      a3mww = zero
      a4muu = zero
      a4muv = zero
      a4muw = zero
      a4mvu = zero 
      a4mvv = zero
      a4mvw = zero
      a4mwu = zero
      a4mwv = zero
      a4mww = zero
      a5muu = zero
      a5muv = zero 
      a5muw = zero
      a5mvu = zero
      a5mvv = zero
      a5mvw = zero
      a5mwu = zero
      a5mwv = zero
      a5mww = zero
      a6muu = zero
      a6muv = zero
      a6muw = zero
      a6mvu = zero
      a6mvv = zero
      a6mvw = zero
      a6mwu = zero
      a6mwv = zero
      a6mww = zero
      a7muu = zero
      a7muv = zero
      a7muw = zero
      a7mvu = zero
      a7mvv = zero
      a7mvw = zero
      a7mwu = zero
      a7mwv = zero
      a7mww = zero
      a8muu = zero
      a8muv = zero
      a8muw = zero
      a8mvu = zero
      a8mvv = zero
      a8mvw = zero
      a8mwu = zero
      a8mwv = zero
      a8mww = zero
      a9muu = zero
      a9muv = zero
      a9muw = zero
      a9mvu = zero
      a9mvv = zero
      a9mvw = zero
      a9mwu = zero
      a9mwv = zero
      a9mww = zero
      a1nuu = zero 
      a1nuv = zero 
      a1nuw = zero
      a1nvu = zero 
      a1nvv = zero
      a1nvw = zero
      a1nwu = zero
      a1nwv = zero
      a1nww = zero
      a2nuu = zero
      a2nuv = zero
      a2nuw = zero
      a2nvu = zero
      a2nvv = zero
      a2nvw = zero
      a2nwu = zero
      a2nwv = zero
      a2nww = zero
      a3nuu = zero
      a3nuv = zero
      a3nuw = zero
      a3nvu = zero
      a3nvv = zero
      a3nvw = zero
      a3nwu = zero
      a3nwv = zero
      a3nww = zero
      a4nuu = zero
      a4nuv = zero
      a4nuw = zero
      a4nvu = zero 
      a4nvv = zero
      a4nvw = zero
      a4nwu = zero
      a4nwv = zero
      a4nww = zero
      a5nuu = zero
      a5nuv = zero 
      a5nuw = zero
      a5nvu = zero
      a5nvv = zero
      a5nvw = zero
      a5nwu = zero
      a5nwv = zero
      a5nww = zero
      a6nuu = zero
      a6nuv = zero
      a6nuw = zero
      a6nvu = zero
      a6nvv = zero
      a6nvw = zero
      a6nwu = zero
      a6nwv = zero
      a6nww = zero
      a7nuu = zero
      a7nuv = zero
      a7nuw = zero
      a7nvu = zero
      a7nvv = zero
      a7nvw = zero
      a7nwu = zero
      a7nwv = zero
      a7nww = zero
      a8nuu = zero
      a8nuv = zero
      a8nuw = zero
      a8nvu = zero
      a8nvv = zero
      a8nvw = zero
      a8nwu = zero
      a8nwv = zero
      a8nww = zero
      a9nuu = zero
      a9nuv = zero
      a9nuw = zero
      a9nvu = zero
      a9nvv = zero
      a9nvw = zero
      a9nwu = zero
      a9nwv = zero
      a9nww = zero
      a1puu = zero 
      a1puv = zero 
      a1puw = zero
      a1pvu = zero 
      a1pvv = zero
      a1pvw = zero
      a1pwu = zero
      a1pwv = zero
      a1pww = zero
      a2puu = zero
      a2puv = zero
      a2puw = zero
      a2pvu = zero
      a2pvv = zero
      a2pvw = zero
      a2pwu = zero
      a2pwv = zero
      a2pww = zero
      a3puu = zero
      a3puv = zero
      a3puw = zero
      a3pvu = zero
      a3pvv = zero
      a3pvw = zero
      a3pwu = zero
      a3pwv = zero
      a3pww = zero
      a4puu = zero
      a4puv = zero
      a4puw = zero
      a4pvu = zero 
      a4pvv = zero
      a4pvw = zero
      a4pwu = zero
      a4pwv = zero
      a4pww = zero
      a5puu = zero
      a5puv = zero 
      a5puw = zero
      a5pvu = zero
      a5pvv = zero
      a5pvw = zero
      a5pwu = zero
      a5pwv = zero
      a5pww = zero
      a6puu = zero
      a6puv = zero
      a6puw = zero
      a6pvu = zero
      a6pvv = zero
      a6pvw = zero
      a6pwu = zero
      a6pwv = zero
      a6pww = zero
      a7puu = zero
      a7puv = zero
      a7puw = zero
      a7pvu = zero
      a7pvv = zero
      a7pvw = zero
      a7pwu = zero
      a7pwv = zero
      a7pww = zero
      a8puu = zero
      a8puv = zero
      a8puw = zero
      a8pvu = zero
      a8pvv = zero
      a8pvw = zero
      a8pwu = zero
      a8pwv = zero
      a8pww = zero
      a9puu = zero
      a9puv = zero
      a9puw = zero
      a9pvu = zero
      a9pvv = zero
      a9pvw = zero
      a9pwu = zero
      a9pwv = zero
      a9pww = zero

      return 
      end