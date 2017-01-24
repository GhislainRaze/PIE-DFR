program Main2D

	! include here the different modules...
	use structure
	use ioGMSH


	! implicit none necessary to declare all variables
	implicit none
	integer :: ierr

	ierr = 0

	! ###
	! ### Step 1: read the input file name - 1st argument of command line...
	! ###
	call getarg(1, Global%inputfile)

	! This indexfile is used to open the mesh several time with the same ID.
	! Never use ID=200 for writing anything !
	global%indexfile = 200

	! ###
	! ### Step 2: read the input file / boundary file
	! ###	
	call readParamFile

	! ###
	! ### Step 3: read the mesh
	! ###	
	call readMeshGMSH

contains

	subroutine readParamFile

		implicit none

		integer :: ios

		open (111,file=trim(adjustl(Global%inputfile)), status="old", &
         action="read", iostat=ios)

		if (ios /=0) then
			print *,"Input file "//trim(adjustl(Global%inputfile))//" NOT found!"; stop
		end if

		rewind(111)
		read(111,*) Global%gridFile
		read(111,*) Global%RiemannSolver
		read(111,*) Global%p
		read(111,*) 
		read(111,*) Global%timeIntegration
		read(111,*) Global%isTimeStepConstant
		read(111,*) Global%CFL
		read(111,*) Global%timeStep
		read(111,*) ! now, time to treat bocos

		read(111,*) InputMarkers%nMarkers
		if(InputMarkers%nMarkers .gt. 0) then
			allocate(InputMarkers%name(InputMarkers%nMarkers))
			allocate(InputMarkers%type(InputMarkers%nMarkers))

			! From nom on, depending on the type, possibility to allocate dedicated 
			! data for the boundary conditions. 
			! allocate(InputMarkers%markerData(InputMarkers%nMarkers))
			! And for any marker, associate the number of data to read (nData) using 
			! InputMarkers%markerData(iData)%nData and allocate and read the data stored in 
			! InputMarkers%markerData(iData)%dataName(1::nData) and 
			! InputMarkers%markerData(iData)%dataValue(1::nData)
		end if

		read(111,*) ! join on the block, for instance periodicity
		read(111,*) InputPerMarkers%nPerMarkers
		if (InputPerMarkers%nPerMarkers .gt. 0) then
			allocate(InputPerMarkers%joinName  (InputPerMarkers%nPerMarkers))
			allocate(InputPerMarkers%opJoinName(InputPerMarkers%nPerMarkers))
			allocate(InputPerMarkers%type      (InputPerMarkers%nPerMarkers))
			allocate(InputPerMarkers%markerData(InputPerMarkers%nPerMarkers))

			! From nom on, depending on the type, possibility to allocate dedicated 
			! data for the join condition. 
			! allocate(InputPerMarkers%markerData(InputPerMarkers%nMarkers))
			! And for any marker, associate the number of data to read (nData) using 
			! InputPerMarkers%markerData(iData)%nData and allocate and read the data stored in 
			! InputPerMarkers%markerData(iData)%dataName(1::nData) and 
			! InputPerMarkers%markerData(iData)%dataValue(1::nData)
		end if

		close(111)


	end subroutine readParamFile




end program Main2D


