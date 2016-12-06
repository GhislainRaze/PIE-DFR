module structure

  use constants

  implicit none

  !************************
  !*** Global variables ***
  !************************
  type strGlobal
    character(len = 240):: inputfile
    character(len = 240):: gridFile
    integer             :: indexfile
    integer             :: RiemannSolver
    integer             :: p
    integer             :: timeIntegration
    real(DP)            :: CFL
    real(DP)            :: maxTime
    logical             :: isTimeStepConstant
    real(DP)            :: timeStep

  end type strGlobal
  
  !**************************
  !*** Gas chemistry data ***
  !**************************
  type strGasChemistry
    real(DP) :: gamma
    real(DP) :: Cp
    real(DP) :: Cv
    real(DP) :: invCv
    real(DP) :: Rgas
  end type strGasChemistry

  
  !*******************************
  !*** Computational loop data ***
  !*******************************
  type strUserData
    integer              :: initIter
    integer              :: maxIter
    real(DP)             :: initialPhysicalTime
    real(DP)             :: time
    real(DP)             :: timeStep
    character(len = 240) :: solutionName
    character(len = 240) :: restartName
    integer              :: currentIteration
  end type strUserData


  !*****************************
  !*** Polynomial definition ***
  !*****************************
  type strPolynomialRoots
  !!!  real(DP), dimension(:),    allocatable :: solToFlux
  !!!  real(DP), dimension(:,:),  allocatable :: solToNode
  !!!  real(DP), dimension(:),    allocatable :: weightOfQuadRule
  !!!  real(DP), dimension(:,:),  allocatable :: interpolationScheme
  !!!  real(DP), dimension(:),    allocatable :: fluxToSol
  !!!  real(DP), dimension(:),    allocatable :: fluxPosition
  !!!  real(DP), dimension(:),    allocatable :: solPosition
  !!!  real(DP), dimension(:,:),  allocatable :: solToOutputPoint
  !!!  logical, dimension(:,:),   allocatable :: faceFlags
  end type strPolynomialRoots


  !***************************
  !*** Mesh data structure ***
  !***************************
  type strMesh
    ! number of local vertices
    integer                                   :: n_vert
    ! number of vertex (on the global mesh)
    integer                                   :: n_vert_tot
    ! number of local elements
    integer                                   :: n_cell
    ! number of elements (on the global mesh)
    integer                                   :: n_cell_tot
    ! Total number of solution points on the whole mesh
    integer                                   :: nTotSolutionPoints
    ! Total number of Flux points on the whole mesh
    integer                                   :: nTotFluxPoints
    ! current timestep
    real(DP)                                  :: dtMin
    ! Type of the element (GMSH)
    integer, dimension(:),        allocatable :: ctype
  end type strMesh

  !**********************
  !*** Output buffers ***
  !**********************
  type strBuf
  !!!   ! Iteration of the last computation of output
  !!!   integer :: sol_ite
  !!!   ! Solution stored on output points
  !!!   real(DP), dimension(:), allocatable :: solOP
  !!!   ! Solution stored on output points (temporary buffer for one cell)
  !!!   real(DP), dimension(:), allocatable :: elt_solOP
  !!!   ! Solution stored on solution points (temporary buffer for one cell)
  !!!   real(DP), dimension(:), allocatable :: elt_solSP
  !!!   ! Buffer to compute cell physical solutions
  !!!   real(DP), dimension(:), allocatable :: elt_phySolSP
  !!!   ! Buffer to store solution on one face
  !!!   real(DP), dimension(:), allocatable :: face_solOP
  end type strBuf

  !*********************************
  !*** Face-based data structure ***
  !*********************************
  type strFace
    ! total number of faces = nFacesInBoco+nFacesIntern
    integer                                :: nFaces
    ! total number of faces for bocos and joins
    integer                                :: nFacesInBoco
    ! Total number of intern faces
    integer                                :: nFacesIntern
    ! List of the four nodes defining the face
    integer,  dimension(:, :), allocatable :: indicNodes
    ! intern = 1, boundary and join for others...
    integer,  dimension(:),    allocatable :: status
    ! index of the left volume
    integer,  dimension(:),    allocatable :: leftCell
    ! index of the right volume
    integer,  dimension(:),    allocatable :: rightCell
    ! -1 = imin, 1 = imax,
    ! -2 = jmin, 2 = jmin,
    ! -3 = kmin, 3 = kmax in the L/R cells
    integer,  dimension(:, :), allocatable :: faceType
    integer,  dimension(:, :), allocatable :: faceLoc
    ! Left and right solution fields before Riemann solver
    real(DP), dimension(:, :), allocatable :: sol_l
    real(DP), dimension(:, :), allocatable :: sol_r
    ! flux values due to riemann solver
    real(DP), dimension(:, :), allocatable :: flux
    ! for face n, represents sum of NoTotSolP for face 1 to n-1
    integer,  dimension(:),    allocatable :: faceDataAddress
    ! matrix to perform index transformation from left to right
    integer,  dimension(:, :), allocatable :: transformation
    ! Face flux point locations on the left side 
    integer, dimension(:, :), allocatable :: Flxpts_L
    ! Face flux point locations on the right side 
    integer, dimension(:, :), allocatable :: Flxpts_R
    ! normal vector components
    real(DP), dimension(:, :), allocatable :: normalVector
    ! face area
    real(DP), dimension(:),    allocatable :: area
  end type strFace

  type strJoinBoundaryCondition
    ! G2N: en cas de p variable, la structure est incomplete...
    integer                                        :: nFaces
    integer                                        :: movement
    integer                                        :: oppositeZoneIndex
    integer                                        :: oppositeJoinIndex
    integer                                        :: convSendRequest
    integer                                        :: convRecvRequest
    integer                                        :: coorSendRequest
    integer                                        :: sizeBufferSConv
    integer                                        :: sizeBufferRConv
    integer                                        :: sizeSendBufferHexa
    integer                                        :: sizeRecBufferHexa
    integer                                        :: sizeSendBufferFace
    integer                                        :: sizeRecBufferFace
    integer                                        :: sizeSendBufferCoor
    integer                                        :: sizeRecBufferCoor
    integer                                        :: idGmsh
    integer,          dimension(:),    allocatable :: indicFaceLocal
    integer,          dimension(:),    allocatable :: indicFaceopposite
    integer,          dimension(:),    allocatable :: indicCellLocal
    integer,          dimension(:, :), allocatable :: localHexa
    integer,          dimension(:, :), allocatable :: Flxpts_L
    integer,          dimension(:, :), allocatable :: Flxpts_R
    real(JR),         dimension(:),    allocatable :: bufferRConv
    real(JR),         dimension(:),    allocatable :: bufferSConv
    integer,          dimension(:),    allocatable :: SendBufferHexa
    integer,          dimension(:),    allocatable :: RecBufferHexa
    double precision, dimension(:),    allocatable :: SendBufferCoor
    double precision, dimension(:),    allocatable :: RecBufferCoor
    integer,          dimension(:),    allocatable :: bufferSMesh
    integer,          dimension(:),    allocatable :: bufferRMesh
    ! not checked yet!
    real(JR),         dimension(:),    allocatable :: rotationCenter
    ! not checked yet!
    real(JR),         dimension(:),    allocatable :: rotationAngle
    ! not checked yet!
    real(JR),         dimension(:),    allocatable :: translation
  end type strJoinBoundaryCondition


  !**********************************
  !*** Runge-Kutta data structure ***
  !**********************************
  type strRungeKutta ! Runge Kutta stepping
    ! kind of Runge Kutta time integation
    integer                             :: type
    ! number of Runge Kutta steps
    integer                             :: nSteps
    ! current Runge Kutta step
    integer                             :: iStep
    ! alpha coefficient for integration
    real(DP), dimension(:), allocatable :: alpha
    ! beta coefficient for integration
    real(DP), dimension(:), allocatable :: beta
  end type strRungeKutta


  !*******************
  !*** Inflow data ***
  !*******************
  type strInfiniteState
    real(DP) :: ro
    real(DP) :: rou
    real(DP) :: rov
    real(DP) :: row
    real(DP) :: roE
    real(DP) :: u
    real(DP) :: v
    real(DP) :: w
    real(DP) :: gamma
    real(DP) :: Cv
    real(DP) :: uNorm
    real(DP) :: Mach
    real(DP) :: Reynolds
    real(DP) :: p
    real(DP) :: T
    real(DP) :: Rgas
    real(DP) :: mu
    real(DP) :: alpha
    real(DP) :: beta
  end type strInfiniteState


  !******************************************************************
  !*** Structures dedicated to the input file boundary conditions ***
  !******************************************************************
  type strMarkerData
    integer                                       :: nData
    character(len=512), dimension(:), allocatable :: dataName
    real(DP)          , dimension(:), allocatable :: dataValue
  end type strMarkerData

  type strMarkersFromInput
    integer                                         :: nMarkers
    character(len=512)  , dimension(:), allocatable :: name
    integer             , dimension(:), allocatable :: type
    type (strMarkerData), dimension(:), allocatable :: markerData
  end type strMarkersFromInput

  type strPeriodicMarkersFromInput
    integer                                         :: nPerMarkers
    character(len=512)  , dimension(:), allocatable :: joinName
    character(len=512)  , dimension(:), allocatable :: opJoinName
    integer             , dimension(:), allocatable :: type
    type (strMarkerData), dimension(:), allocatable :: markerData
  end type strPeriodicMarkersFromInput

  !**********************************************************************
  !*** Structures dedicated to the GMSH input for boundary conditions ***
  !**********************************************************************
  type strGMSH_Family
    character(len=512) :: name
    integer            :: flag
    integer            :: idMarker
    integer            :: nFaces
  end type strGMSH_Family


  !****************************************************
  !*** Dedicated structures for boundary conditions ***
  !****************************************************
  !==> Reference state
  type strReferenceState
    real(DP) :: ro
    real(DP) :: rou
    real(DP) :: rov
    real(DP) :: row
    real(DP) :: roE
  end type strReferenceState


  ! ==> boundary condition structure filled when reading the mesh
  type strPhysBoundaryCondition
    ! kind of boundary condition
    integer                                           :: bocoType
    integer                                           :: idGmsh

    ! number of faces / local to global face index
    integer                                           :: nFaces
    integer , dimension(:)           ,    allocatable :: indicFace

    ! list of adjacent volume (true one is always on left side)
    ! ghost cell is necessary on right side
    integer,             dimension(:),    allocatable :: indicCell
    integer,             dimension(:),    allocatable :: indicGhostCell

    ! structures dedicated to data needed by the boundary condition
    type(strReferenceState)                           :: referenceState
    type(strInjection)                                :: referenceInjection
    type(strTurbInjection)                            :: referenceTurbInjection
    type(strViscousWall)                              :: viscousWall
    real(DP),            dimension(:),    allocatable :: outflowPressure
    real(DP),            dimension(:, :), allocatable :: bndSolution

    ! vector of the derivative of the basis function on the face
    real(DP),            dimension(:, :), allocatable :: dLflux
  end type strPhysBoundaryCondition



  !**************************************************************
  !*** save all variables for a static-like memory management ***
  !**************************************************************
  type (strGlobal)                                     :: global
  type (strGasChemistry)                               :: gas
  type (strUserData)                                   :: userData
  type (strMarkersFromInput)                           :: InputMarker
  type (strPeriodicMarkersFromInput)                   :: InputPerMarker

  !size of pmax, permits to use different cube types
  type (strPolynomialRoots), dimension(:), allocatable :: polyRoots

  !size of pmax, permits to use different cube types
  type (strMesh)                                             :: Mesh
  type (strFace)                                             :: Face
  type (strGMSH_Family),           dimension(:), allocatable :: GMSH_Family

  ! number of boundary conditions recovered from Mesh file
  integer                                              :: nBocos

  ! number of joins  recovered from Mesh file
  integer                                              :: nJoins

  ! number of output surfaces
  integer                                              :: nOutSurface

  real(DP), dimension(:),   allocatable                  :: physGradientInSP
  real(DP), dimension(:),   allocatable                  :: tmpGradientInFP
  real(DP), dimension(:),   allocatable                  :: compGradientInSP
  real(DP), dimension(:),   allocatable                  :: correctedSolution
!!$  real(DP), dimension(:),   allocatable                  :: physicalSolution
!!$  real(DP), dimension(:,:),   allocatable                :: extraSolution
!!$  real(DP), dimension(:,:), allocatable                  :: coordOP
!!$  real(DP), dimension(:,:), allocatable                  :: solutionInOP
!!$  real(DP), dimension(:), allocatable                    :: outputPoints

  logical :: isCWnodeOrder = .false.

  !=============================================================================!
end module structure