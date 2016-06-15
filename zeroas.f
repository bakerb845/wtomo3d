      SUBROUTINE zeroas()
c 
c     Zero out each of the a's
      USE MAT_MODULE

      COMPLEX, PARAMETER :: zzero = CMPLX(0.0, 0.0)

      a1muu = zzero 
      a1muv = zzero 
      a1muw = zzero
      a1mvu = zzero 
      a1mvv = zzero
      a1mvw = zzero
      a1mwu = zzero
      a1mwv = zzero
      a1mww = zzero

      a2muu = zzero
      a2muv = zzero
      a2muw = zzero
      a2mvu = zzero
      a2mvv = zzero
      a2mvw = zzero
      a2mwu = zzero
      a2mwv = zzero
      a2mww = zzero

      a3muu = zzero
      a3muv = zzero
      a3muw = zzero
      a3mvu = zzero
      a3mvv = zzero
      a3mvw = zzero
      a3mwu = zzero
      a3mwv = zzero
      a3mww = zzero

      a4muu = zzero
      a4muv = zzero
      a4muw = zzero
      a4mvu = zzero 
      a4mvv = zzero
      a4mvw = zzero
      a4mwu = zzero
      a4mwv = zzero
      a4mww = zzero

      a5muu = zzero
      a5muv = zzero 
      a5muw = zzero
      a5mvu = zzero
      a5mvv = zzero
      a5mvw = zzero
      a5mwu = zzero
      a5mwv = zzero
      a5mww = zzero

      a6muu = zzero
      a6muv = zzero
      a6muw = zzero
      a6mvu = zzero
      a6mvv = zzero
      a6mvw = zzero
      a6mwu = zzero
      a6mwv = zzero
      a6mww = zzero

      a7muu = zzero
      a7muv = zzero
      a7muw = zzero
      a7mvu = zzero
      a7mvv = zzero
      a7mvw = zzero
      a7mwu = zzero
      a7mwv = zzero
      a7mww = zzero

      a8muu = zzero
      a8muv = zzero
      a8muw = zzero
      a8mvu = zzero
      a8mvv = zzero
      a8mvw = zzero
      a8mwu = zzero
      a8mwv = zzero
      a8mww = zzero

      a9muu = zzero
      a9muv = zzero
      a9muw = zzero
      a9mvu = zzero
      a9mvv = zzero
      a9mvw = zzero
      a9mwu = zzero
      a9mwv = zzero
      a9mww = zzero

      a1nuu = zzero 
      a1nuv = zzero 
      a1nuw = zzero
      a1nvu = zzero 
      a1nvv = zzero
      a1nvw = zzero
      a1nwu = zzero
      a1nwv = zzero
      a1nww = zzero

      a2nuu = zzero
      a2nuv = zzero
      a2nuw = zzero
      a2nvu = zzero
      a2nvv = zzero
      a2nvw = zzero
      a2nwu = zzero
      a2nwv = zzero
      a2nww = zzero

      a3nuu = zzero
      a3nuv = zzero
      a3nuw = zzero
      a3nvu = zzero
      a3nvv = zzero
      a3nvw = zzero
      a3nwu = zzero
      a3nwv = zzero
      a3nww = zzero

      a4nuu = zzero
      a4nuv = zzero
      a4nuw = zzero
      a4nvu = zzero 
      a4nvv = zzero
      a4nvw = zzero
      a4nwu = zzero
      a4nwv = zzero
      a4nww = zzero

      a5nuu = zzero
      a5nuv = zzero 
      a5nuw = zzero
      a5nvu = zzero
      a5nvv = zzero
      a5nvw = zzero
      a5nwu = zzero
      a5nwv = zzero
      a5nww = zzero

      a6nuu = zzero
      a6nuv = zzero
      a6nuw = zzero
      a6nvu = zzero
      a6nvv = zzero
      a6nvw = zzero
      a6nwu = zzero
      a6nwv = zzero
      a6nww = zzero

      a7nuu = zzero
      a7nuv = zzero
      a7nuw = zzero
      a7nvu = zzero
      a7nvv = zzero
      a7nvw = zzero
      a7nwu = zzero
      a7nwv = zzero
      a7nww = zzero

      a8nuu = zzero
      a8nuv = zzero
      a8nuw = zzero
      a8nvu = zzero
      a8nvv = zzero
      a8nvw = zzero
      a8nwu = zzero
      a8nwv = zzero
      a8nww = zzero

      a9nuu = zzero
      a9nuv = zzero
      a9nuw = zzero
      a9nvu = zzero
      a9nvv = zzero
      a9nvw = zzero
      a9nwu = zzero
      a9nwv = zzero
      a9nww = zzero

      a1puu = zzero 
      a1puv = zzero 
      a1puw = zzero
      a1pvu = zzero 
      a1pvv = zzero
      a1pvw = zzero
      a1pwu = zzero
      a1pwv = zzero
      a1pww = zzero

      a2puu = zzero
      a2puv = zzero
      a2puw = zzero
      a2pvu = zzero
      a2pvv = zzero
      a2pvw = zzero
      a2pwu = zzero
      a2pwv = zzero
      a2pww = zzero

      a3puu = zzero
      a3puv = zzero
      a3puw = zzero
      a3pvu = zzero
      a3pvv = zzero
      a3pvw = zzero
      a3pwu = zzero
      a3pwv = zzero
      a3pww = zzero

      a4puu = zzero
      a4puv = zzero
      a4puw = zzero
      a4pvu = zzero 
      a4pvv = zzero
      a4pvw = zzero
      a4pwu = zzero
      a4pwv = zzero
      a4pww = zzero

      a5puu = zzero
      a5puv = zzero 
      a5puw = zzero
      a5pvu = zzero
      a5pvv = zzero
      a5pvw = zzero
      a5pwu = zzero
      a5pwv = zzero
      a5pww = zzero

      a6puu = zzero
      a6puv = zzero
      a6puw = zzero
      a6pvu = zzero
      a6pvv = zzero
      a6pvw = zzero
      a6pwu = zzero
      a6pwv = zzero
      a6pww = zzero

      a7puu = zzero
      a7puv = zzero
      a7puw = zzero
      a7pvu = zzero
      a7pvv = zzero
      a7pvw = zzero
      a7pwu = zzero
      a7pwv = zzero
      a7pww = zzero

      a8puu = zzero
      a8puv = zzero
      a8puw = zzero
      a8pvu = zzero
      a8pvv = zzero
      a8pvw = zzero
      a8pwu = zzero
      a8pwv = zzero
      a8pww = zzero

      a9puu = zzero
      a9puv = zzero
      a9puw = zzero
      a9pvu = zzero
      a9pvv = zzero
      a9pvw = zzero
      a9pwu = zzero
      a9pwv = zzero
      a9pww = zzero

      return 
      end
