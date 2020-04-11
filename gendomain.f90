!	gendomain.f90 generates a domain file and mass input file to be used in FDPML method

!	input file structure
!	&system
!		PD(3) - size of the primary domain
!		crystal_coordinates - .true. if the domain is to be generated in crystal coordinates.
!								if .false. the domain is generated in conventional lattice
!		domain_prefix - file where the domain specification is stored
!		mass_prefix - file where the mass specification is stored
!		mass_input - logical - indicates if the mass specification file needs to generated
!								if .false. FDPML can generated masses on an individual unit cell basis
!		flfrc1 - force constant file for the material type 1. Generated via Quantum Espresso (NOT IN XML)
!		flfrc2 - force constant file for the material type 2. Generated via Quantum Espresso (NOT IN XML)
!	/
!	&domain
!		nanoparticle - .true. if the domain defines a nanoparticle
!						if .false. the domain is of type interface
!		Radius - radius of the nanoparticle in unit cells
!		random - if .true.  the domain generated is an alloy of materials type 1 and 2
!		c - alloy concentration for domain type 2 i.e. type_1_{1-c}type_2_{c}
!	/


PROGRAM gendomain

	USE mp_module
	USE io_module
	USE kinds
	USE lapack95
	USE blas95
	USE f95_precision
	USE essentials
	use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
											  stdout=>output_unit, &
											  stderr=>error_unit
	
	IMPLICIT NONE
	
	INCLUDE 'mpif.h' ! MPI header file
	INCLUDE 'mkl.fi' ! MKL header file
	
	CHARACTER(LEN = 256) :: flfrc1, flfrc2, mass_prefix, domain_prefix
	REAL(KIND = 4) :: PD(3), center(3), zmid, r(3), rn
	LOGICAL :: nanoparticle, random, core_shell ! if nanoparticle is false interface = .true.
	REAL(KIND = RP) :: Radius, c
	REAL(KIND = RP), ALLOCATABLE :: amass_PD(:,:,:,:)
	INTEGER, ALLOCATABLE :: ityp_PD(:,:,:,:)
	INTEGER :: i, n1,n2, n3, na
	INTEGER :: nr1, nr2, nr3, nat(2), ibrav(2), ntyp(2), natsc
	REAL(KIND = RP) :: epsil(3,3,2), at1(3,3), at2(3,3), at_conv(3,3)
	REAL(KIND = RP) :: alat1, alat2, omega1, omega2, q(3), distance
	REAL(KIND = RP), ALLOCATABLE :: amass1(:), amass2(:)
	LOGICAL :: has_zstar(2), fd, na_ifc
	LOGICAL :: crystal_coordinates, mass_input
	REAL(KIND = RP), ALLOCATABLE :: frc1(:,:,:,:,:,:,:), tau1(:,:),  zeu1(:,:,:), m_loc(:,:)
	INTEGER, ALLOCATABLE  :: ityp1(:)
	REAL(KIND = RP), ALLOCATABLE :: frc2(:,:,:,:,:,:,:), tau2(:,:),  zeu2(:,:,:)
	INTEGER, ALLOCATABLE  :: ityp2(:)
	INTEGER, ALLOCATABLE :: itypsc(:)
	REAL(KIND = RP), ALLOCATABLE :: tausc(:,:)
	REAL :: yplane
	REAL :: counter
	REAL :: natoms
	
	
	CALL mp_init ( ) ! initialize MPI
	
	IF (io_node) THEN
		write (stdout, *) 'Genrating Domain for Frequency Domain Perfectly Matched Layer'
		write (stdout, *) '	'
	ENDIF
	
	
	
	NAMELIST /system/ PD, crystal_coordinates, domain_prefix, mass_prefix, mass_input, flfrc1, flfrc2, &
			 /domain/ nanoparticle, random, Radius, c

	nanoparticle = .false.
	Radius = 0.D0
	random = .false.
			 
	IF (root_node) THEN
		READ(stdin, system)
		READ(stdin, domain)
	ENDIF
	
	IF (io_node) THEN
		IF (.not. crystal_coordinates) THEN
			WRITE (stdout, *) 'Generating domain in convetional co-ordinate system (cartesian)'
		ELSE
			WRITE (stdout, *) 'Generating domain in crystal co-ordinate system'
		ENDIF
		
	ENDIF
	
	IF (nanoparticle .and. (Radius.eq.0)) THEN
		WRITE(stdout, '(a)') ' Radius == 0'
		CALL MPI_ABORT(comm, ierr)
	ENDIF
	
	CALL readfc( flfrc1, frc1, tau1, zeu1, m_loc, ityp1, nr1, nr2, nr3, epsil(:,:,1),&
				 nat(1), ibrav(1), alat1, at1, ntyp(1), amass1, omega1, has_zstar(1) )
	CALL readfc( flfrc2, frc2, tau2, zeu2, m_loc, ityp2, nr1, nr2, nr3, epsil(:,:,2),&
				 nat(2), ibrav(2), alat2, at2, ntyp(2), amass2, omega2, has_zstar(2) )

	nr1 = PD(1)
	nr2 = PD(2)
	nr3 = PD(3)
	
	IF (nanoparticle .or. core_shell) THEN
		center = PD/2 + (/1.0, 1.0, 1.0/)
	ELSE
		zmid = nr3/2
	ENDIF

	IF (.not. crystal_coordinates) THEN
		at_conv(:,:) = 0.D0
		DO i = 1, 3
			at_conv(i,i) = 1.D0
		ENDDO
		CALL Supercell(at_conv, at1, tau1, ityp1, nat(1), tausc, itypsc, natsc)
		DEALLOCATE (tau1, ityp1, ityp2)
		ALLOCATE(tau1(3, natsc), ityp1(natsc), ityp2(natsc))
		at1 = at_conv
		tau1 = tausc
		ityp1 = itypsc
		ityp2 = itypsc
	ENDIF
	
	ALLOCATE(ityp_PD(nr1,nr2,nr3,natsc), amass_PD(nr1,nr2,nr3,natsc))
	
	DO n1 = 1, nr1
		DO n2 = 1, nr2
			DO n3 = 1, nr3
				DO na = 1, natsc
					r = n1*at_conv(:,1) + n2*at_conv(:,2) + n3*at_conv(:,3) + tausc(:,na)
					IF (nanoparticle) THEN
						distance = nrm2((center-r(:)))
						IF (distance.le.Radius) THEN 
							ityp_PD(n1,n2,n3,na) = 2
							amass_PD(n1,n2,n3,na) = amass2(itypsc(na))
						ELSE
							ityp_PD(n1,n2,n3,na) = 1
							amass_PD(n1,n2,n3,na)  = amass1(itypsc(na))
						ENDIF
!					ELSE IF (core_shell) THEN
!						distance = nrm2((center-r(:)))
!						IF (distance.le.Radius(1)) THEN 
!							ityp_PD(n1,n2,n3,na) = 2
!							amass_PD(n1,n2,n3,na) = amass2(itypsc(na))
! 						ELSE IF (distance.le.Radius(2)) THEN 
! 							ityp_PD(n1,n2,n3,na) = 3
! 							amass_PD(n1,n2,n3,na) = amass2(itypsc(na))
!						ELSE
!							ityp_PD(n1,n2,n3,na) = 1
!							amass_PD(n1,n2,n3,na)  = amass1(itypsc(na))
!						ENDIF
					ELSE
						IF (r(3).gt.zmid) THEN
							ityp_PD(n1,n2,n3,na) = 2
							amass_PD(n1,n2,n3,na) = amass2(itypsc(na))
						ELSE
							ityp_PD(n1,n2,n3,na) = 1
							amass_PD(n1,n2,n3,na)  = amass1(itypsc(na))
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	
	counter = 0
	
	CALL init_random_seed()  

	
	IF (random) THEN
		DO n1 = 1, nr1
			DO n2 = 1, nr2
				DO n3 = 1, nr3
					DO na = 1, natsc
						IF (itypsc(na).eq.1) THEN
							call random_number(rn)
							IF (rn.le.c) THEN
								amass_PD(n1,n2,n3,na) = amass2(itypsc(na))
								counter = counter + 1
							ENDIF
						ENDIF
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF
	
	print *, counter
	
	natoms = nr1*nr2*nr3*natsc
	
	print *, natoms
	
	WRITE(stdout, '(a, F10.3)') 'Actual concentration = ', counter/natoms
	
!	amass_PD(nr1/2, nr2/2, nr3/2, 2) = amass2(itypsc(2))

	print *, itypsc
	
	IF(nanoparticle) THEN
		yplane = int(center(2))
	ELSE
		yplane = nr2/2
	ENDIF
	
	IF (io_node) THEN
		WRITE (stdout, *) 'Cross-sectional image of the domain'
		WRITE (stdout, *) '--------------------------------------------------------'
		DO n1 = 1, nr1 
			DO n3 = 1, nr3
				IF (n3.eq.nr3) THEN
					write(stdout, fmt = '(I2)') ityp_PD(n1,yplane,n3,1)
				ELSE
					write(stdout, fmt = '(I2)', advance = 'no') ityp_PD(n1,yplane,n3,1)
				ENDIF
			ENDDO
		ENDDO
	ENDIF
	
	IF (io_node) THEN
		WRITE (stdout, fmt = '(a)') 'Storing primary domain in ', domain_prefix
		IF (mass_input) WRITE (stdout, fmt = '(a)') 'Storing mass of every atom in domain in ', mass_prefix
		
		OPEN (unit = 1, file = domain_prefix, form = 'unformatted')
		write (unit = 1) ityp_PD
		close (unit = 1)
		
		IF (mass_input) THEN
			OPEN (unit = 1, file = mass_prefix, form = 'unformatted')
			write (unit = 1) amass_PD
			close (unit = 1)
		ENDIF
	ENDIF
	
	CALL mp_finalize( )
	
END PROGRAM gendomain

SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE
