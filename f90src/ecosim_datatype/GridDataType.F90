module GridDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
implicit none
  public
  save

  integer, PARAMETER :: JX=4
  integer, PARAMETER :: JY=4
  integer, PARAMETER :: JZ=20
  integer, PARAMETER :: JH=JX+1
  integer, PARAMETER :: JV=JY+1
  integer, PARAMETER :: JD=JZ+1
  integer, PARAMETER :: JP=5
  integer, PARAMETER :: JC=10
  integer, PARAMETER :: JS=5
  integer, PARAMETER :: JLI=4
  integer, PARAMETER :: JLA=4   !# of sectors for the leaf zimuth, [0,pi]
  integer, PARAMETER :: JSA=4   !# of sectors for the sky azimuth  [0,2*pi]

  integer :: NU(JY,JX)
  integer :: NUI(JY,JX)
  integer :: NJ(JY,JX)
  integer :: NK(JY,JX)
  integer :: NLI(JV,JH)
  integer :: NL(JV,JH)
  integer :: NUM(JY,JX)
  real(r8) :: CDPTH(0:JZ,JY,JX)   !depth to bottom of soil layer [m]
  real(r8) :: CDPTHI(JY,JX)
  real(r8) :: DLYR(3,0:JZ,JY,JX)  !thickness of soil layer [m]
  real(r8) :: DLYRI(3,0:JZ,JY,JX) !
  real(r8) :: DPTH(JZ,JY,JX)      !depth to middle of soil layer [m]
  real(r8) :: XDPTH(3,JZ,JY,JX)   !cross-sectional area / distance between adjacent grid cells [m]
  real(r8) :: CDPTHZ(0:JZ,JY,JX)  !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8) :: DPTHZ(JZ,JY,JX)     !depth to middle of soil layer from  surface of grid cell [m]
  real(r8) :: AREA(3,0:JZ,JY,JX)  !cross-sectional area  [m2 d-2]
  real(r8) :: DIST(3,JD,JV,JH)    !distance between adjacent layers:1=EW,2=NS,3=vertical [m]
  REAL(R8) :: ALAT(JY,JX)
end module GridDataType
