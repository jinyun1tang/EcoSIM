module GridDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  integer, PARAMETER :: JX=4
  integer, PARAMETER :: JY=4
  integer, PARAMETER :: JZ=20
  integer, PARAMETER :: JH=JX+1
  integer, PARAMETER :: JV=JY+1
  integer, PARAMETER :: JD=JZ+1
  integer, PARAMETER :: JP=5
  integer, PARAMETER :: JC=10
  integer, PARAMETER :: JS=5
  integer, PARAMETER :: JLI=4   !# of sectors for the leaf zenith [0,pi/2]
  integer, PARAMETER :: JLA=4   !# of sectors for the leaf azimuth, [0,pi]
  integer, PARAMETER :: JSA=4   !# of sectors for the sky azimuth  [0,2*pi]

  integer :: NPX                  !number of E-W grid cells
  integer :: NPY                  !number of N-S grid cells

  integer :: NU(JY,JX)            !soil surface layer number
  integer :: NUI(JY,JX)           !initial soil surface layer number
  integer :: NJ(JY,JX)            !maximum root layer number
  integer :: NK(JY,JX)            !additional soil lower boundary layers
  integer :: NLI(JV,JH)           !initial lowest soil layer number
  integer :: NL(JV,JH)            !lowest soil layer number
  integer :: NUM(JY,JX)           !new surface layer number
  real(r8) :: CDPTH(0:JZ,JY,JX)   !depth to bottom of soil layer [m]
  real(r8) :: CDPTHI(JY,JX)       !initial depth to bottom of soil layer [m]
  real(r8) :: DLYR(3,0:JZ,JY,JX)  !thickness of soil layer [m]
  real(r8) :: DLYRI(3,0:JZ,JY,JX) !thickness of soil layer [m]
  real(r8) :: DPTH(JZ,JY,JX)      !depth to middle of soil layer [m]
  real(r8) :: XDPTH(3,JZ,JY,JX)   !cross-sectional area / distance between adjacent grid cells [m]
  real(r8) :: CDPTHZ(0:JZ,JY,JX)  !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8) :: DPTHZ(JZ,JY,JX)     !depth to middle of soil layer from  surface of grid cell [m]
  real(r8) :: AREA(3,0:JZ,JY,JX)  !cross-sectional area  [m2 d-2]
  real(r8) :: DIST(3,JD,JV,JH)    !distance between adjacent layers:1=EW,2=NS,3=vertical [m]
  REAL(R8) :: ALAT(JY,JX)         !latitude	[degrees]
  real(r8) :: DH(JY,JX)           !number of EW grid cells, [-]
  real(r8) :: DV(JY,JX)           !number of EW grid cells, [-]
  real(r8) :: TAREA               !total area of landscape	[m2]
  integer :: NCN(JY,JX)           !number of dimensions for grid cell connections
  integer :: LSG(JZ,JY,JX)        !match PFT from different scenarios
  integer :: NP(JY,JX)            !number of plant species
  integer :: NP0(JY,JX)           !intitial number of plant species
end module GridDataType
