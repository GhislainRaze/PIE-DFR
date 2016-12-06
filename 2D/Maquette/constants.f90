module constants

  implicit none
  
  integer, parameter :: DP = selected_real_kind(15,307)

  real(DP),   parameter :: MINREAL         = SQRT(EPSILON(1.0_DP))
  real(DP),   parameter :: DP_PI          = 4.0_DP * atan(1._DP)
  real(DP),   parameter :: DP_ZERO        = 0.0_DP
  real(DP),   parameter :: DP_FOURTH      = 0.25_DP
  real(DP),   parameter :: DP_THIRD       = 1.0_DP / 3.0_DP
  real(DP),   parameter :: DP_TWOTHIRD    = 2.0_DP / 3.0_DP
  real(DP),   parameter :: DP_HALF        = 0.5_DP
  real(DP),   parameter :: DP_ONE         = 1.0_DP
  real(DP),   parameter :: DP_TWO         = 2.0_DP
  real(DP),   parameter :: DP_THREE       = 3.0_DP
  real(DP),   parameter :: DP_FOUR        = 4.0_DP
  real(DP),   parameter :: DP_FIVE        = 5.0_DP
  real(DP),   parameter :: DP_SIX         = 6.0_DP
  real(DP),   parameter :: DP_cutoff      = 1.e-15_DP
  real(DP),   parameter :: DP_min_surface = 1.e-30_DP
  real(DP),   parameter :: DP_min_volume  = 1.e-20_DP
  real(DP),   parameter :: DP_DIST_TOL    = 1.e-6_DP
  integer,    parameter :: DP_YES         = 1
  integer,    parameter :: DP_NO          = 0


  integer, parameter :: nVars             = 5
  integer, parameter :: PERIODIC          = 10

  !***************************************************************************
  !****** GMSH element-type values                                                        *
  !****** Check http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format *  
  !***************************************************************************
  integer, parameter   :: GMSH_POINT = 15
  integer, parameter   :: GMSH_LINE  = 1
  integer, parameter   :: GMSH_QUAD  = 3
  integer, parameter   :: GMSH_HEXA  = 5
  integer, parameter   :: GMSH_LINE_3N  = 8
  integer, parameter   :: GMSH_QUAD_9N  = 10
  integer, parameter   :: GMSH_HEXA_27N  = 12

  integer, parameter   :: CONV_ROE      = 1
  integer, parameter   :: CONV_AUSMP    = 2
  integer, parameter   :: CONV_SLAU     = 3
  integer, parameter   :: CONV_RUSANOV  = 4

  !=============================================================================!
end module constants