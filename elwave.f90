      USE INIT_MODULE
      USE MODEL_MODULE
      USE STRUMPACK_MODULE
      USE RECPRM_MODULE
      USE SRCPRM_MODULE
      USE INTERFACE_MODULE
      USE SPARSE_MODULE
      !USE MPI
      CHARACTER(256) projnm
      COMPLEX, DIMENSION(:), ALLOCATABLE :: a, rhs, sol
      INTEGER, DIMENSION(:), ALLOCATABLE :: xadj, adjncy, irptr, jcptr, part
      INTEGER ierr, mpierr, myid, nedge, numvert, nprocs
      LOGICAL is3d
      INTEGER, PARAMETER :: master = 0
      REAL, PARAMETER :: pi = 3.141592653589793 
mpierr = 0
myid =0
nprocs = 0
nparts = 1
nfreq_groups = 1
myfreq_group = 0
      ierr = 0
! 
!.... initialize mpi 
!     CALL MPI_INIT(mpierr)
!     CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierr)
!     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierr)
      IF (myid == master) THEN
         projnm(:) = ' '
projnm(:) = 'test' !TODO: fix
         WRITE(*,*) 'elwave: Enter project name'
         !READ(*,'(A)') projnm
         WRITE(*,*) 'elwave: Reading ini file...'
         CALL READINI_FORWARD(projnm, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'elwave: Error reading ini file!'
            GOTO 500
         ENDIF
         WRITE(*,*) 'elwave: Loading model file...'
         CALL MODELIO_READ_MODEL(projnm, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'elwave: Error reading model!'
            GOTO 500
         ENDIF
         ! get the source information
         WRITE(*,*) 'elwave: Setting source time functions...'
         nst = ns
         CALL SRCSUB(projnm, &
                     ncom, nsam, nom, nomi, &
                     nst, modsrcp, &
                     srcin, nsg, isg, &
                     deltatt, tau, freq, &
                     dsd, source, ierr) 
         IF (ierr /= 0) THEN
            WRITE(*,*) 'elwave: Error initializing STFs!'
            GOTO 500
         ENDIF
         WRITE(*,*) 'elwave: Initializing receiver response functions...'
         CALL RECSUB(projnm, &
                      ncom, nsam, nom, ng, &
                      recin, &
                      deltatt, freq, &
                      receiver, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'elwave: Error initializing RRFs!'
            GOTO 500
         ENDIF
         ! initialize space for estiamtes
         WRITE(*,*) 'elwave: Nulling estimates...'
         IF (ng > 0) THEN
            ALLOCATE(uest(nom, ng, ns))
            ALLOCATE(vest(nom, ng, ns))
            ALLOCATE(west(nom, ng, ns))
            CALL NULLUEST(nom, ng, ns, uest, vest, west)
         ENDIF
         IF (nr > 0) THEN
            ALLOCATE(utest(nom, nr, ns))
            CALL NULLUTEST(nom, nr, ns, utest)
         ENDIF
         CALL ZERO()
         ! display model information
         IF (ncom > 1) CALL INFO(nx, ny, nz, ns, nr, ng, nom)
         ! determine if model is 3d
         CALL TEST1D(is3d)
         IF (is3d) THEN
            WRITE(*,*) 'elwave: Model is 3D'
         ELSE
            WRITE(*,*) 'elwave: Model is 1D'
         ENDIF
         ! determine the full graph size
         WRITE(*,*) 'elwave: Computing graph...'
         CALL SREL3D_GETSIZE(nx, ny, nz, numvert, nedge)
         ! compute the graph
         ALLOCATE(xadj(numvert+1))
         ALLOCATE(adjncy(2*nedge))
         CALL SREL3D(nx, ny, nz, numvert, nedge, xadj, adjncy)
         ! compute the edge cuts
         WRITE(*,*) 'elwave: Partitioning graph...' 
         ALLOCATE(part(numvert))
         nxyz = nx*ny*nz
         nvloc = 3*MAX0(nxyz/nparts, 1) ! want process to have all 3 rows of a node
         ipart = 1
         part(:) = 0
         DO i=1,numvert
            IF (i > ipart*nvloc) ipart = ipart + 1 ! new partition
            IF (ipart > nparts) ipart = ipart - 1  ! don't exceed number of processes
            part(i) = ipart
         ENDDO
         IF (MINVAL(part) == 0 .OR. MAXVAL(part) > nparts) THEN
            WRITE(*,*) 'elwave: Error partitioning array!'
            ierr = 1
            GOTO 500
         ENDIF
         IF (ierr /= 0) THEN
            WRITE(*,*) 'elwave: Error partitioning mesh!'
            GOTO 500
         ENDIF
         ! now copy the subblocks and make a local to global map
         ALLOCATE(irptr(numvert+1))
         nzero = 2*nedge + numvert
         ALLOCATE(jcptr(nzero))
         CALL SPARSE_ADJ2CRS(nzero, numvert, xadj, adjncy, irptr, jcptr)
         IF (ALLOCATED(xadj))   DEALLOCATE(xadj)
         IF (ALLOCATED(adjncy)) DEALLOCATE(adjncy)

      ENDIF
      ! create the processor groups
      IF (myid == master) WRITE(*,*) 'elwave: Generating process groups...'

      IF (ierr /= 0) THEN
         WRITE(*,*) 'elwave: An error has occurred on process', myid
         !CALL MPI_ABORT(MPI_COMM_WORLD, 30, mpierr) 
      ENDIF
      ! initialize the matrix
      ALLOCATE(a(nzero))
      ALLOCATE(rhs(numvert))
      ALLOCATE(sol(numvert))
      a(:) = CMPLX(0.0, 0.0)
      rhs(:) = CMPLX(0.0, 0.0)
      sol(:) = CMPLX(0.0, 0.0) 
      CALL STRUMPACK_SERIAL_FINTER(1, numvert, nzero, &
                                   irptr, jcptr, a,   &
                                   rhs, sol, ierr)
      !----------------------------------------------------------------------------------!
      !                                begin the frequency loop                          !
      !----------------------------------------------------------------------------------!
      DO 1000 ifreql=1,nom
         ! get the frequency
         iom = (ifreql - 1)*nfreq_groups + myfreq_group + 1
         IF (iom > nom) CYCLE
         omega = CMPLX(2.0*pi*freq(iom), 0.0)
         ! tell user program is doing something 
         IF (myid == master) THEN
            DO if_group=1,nfreq_groups
               jom = iom + if_group - 1
               IF (jom <= nom) WRITE(*,999) if_group, freq(jom)
  999          FORMAT(/,15x,' ********************************************',     &
                         /,15x,' *    Group ',i5'                           *',  &
                         /,15x,' *    Processing ',e12.5, ' hz            *',    &   
                         /,15x,' ********************************************',/)
            ENDDO
         ENDIF
         ! anti-alias filter
         IF (tau /= 999.999) omega = omega + CMPLX(0.0, 1.0/tau)
         ! apply the dispersion model
         IF (qpex) THEN
            IF (myid == master) WRITE(*,*) 'elwave: Warning dispersion is not yet tested'
            CALL DISPERSION()
         ENDIF
         IF (myid == master) WRITE(*,*) 'elwave: Assembling model...'
         CALL ASMBLE(lunlog, nzero, nx, ny, nz, &
                     irptr, jcptr, a, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'elwave: Error assembling matrix'
            GOTO 500
         ENDIF
         IF (myid == master) WRITE(*,*) 'elwave: Factoring matrix...'
         CALL STRUMPACK_SERIAL_FINTER(2, numvert, nzero, &
                                      irptr, jcptr, a,   &
                                      rhs, sol, ierr)
         ! solve phase
 1000 CONTINUE
  500 CONTINUE
      IF (ALLOCATED(irptr)) DEALLOCATE(irptr)
      IF (ALLOCATED(jcptr)) DEALLOCATE(jcptr)
      IF (ALLOCATED(part))  DEALLOCATE(part)
      IF (ALLOCATED(a))     DEALLOCATE(a)
      CALL STRUMPACK_SERIAL_FINTER(-1, numvert, nzero, &
                                   irptr, jcptr, a,   &
                                   rhs, sol, ierr)
      CALL MEMORY_FREE_HASKELL()
      CALL MEMORY_FREE_MODEL()
      CALL MEMORY_FREE_RECEIVERS()
      CALL MEMORY_FREE_SOURCE()
      !CALL MPI_FINALIZE()
      STOP
      END
