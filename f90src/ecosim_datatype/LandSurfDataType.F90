module LandSurfDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: ZS(JY,JX)                         !initial soil surface roughness height, [m]
  real(r8) :: ZD(JY,JX)                         !zero plane displacement height, [m]
  real(r8) :: ZR(JY,JX)                         !canopy surface roughness height, [m]
  real(r8) :: ZM(JY,JX)                         ! soil surface roughness height for calculating runoff velocity, [m]
  real(r8) :: Z0(JY,JX)                         !wind speed measurement height, [m]
  real(r8) :: ALT(JY,JX)                        !altitude of grid cell, [m]
  real(r8) :: RAB(JY,JX)                        !isothermal boundary layer resistance, [h m-1]
  real(r8) :: RIB(JY,JX)                        !Richardson number for calculating boundary layer resistance, [-]
  real(r8) :: ALTI(JY,JX)                       !altitude of landscape, [m]
  real(r8) :: GSIN(JY,JX)                       !sine of slope, [-]
  real(r8) :: GCOS(JY,JX)                       !cosine of slope, [-]
  real(r8) :: GAZI(JY,JX)                       !azimuth of slope, [-]
  real(r8) :: ALTIG                             !altitude of landscape, [m]
  real(r8) :: ALTZ(JY,JX)                       !altitude, [m]
  real(r8) :: SL(JY,JX)                         !slope, [o]
  real(r8) :: ASP(JY,JX)                        !aspect , [o]

  real(r8) :: XCODFS(JY,JX)                     !surface - atmosphere CO2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XCHDFS(JY,JX)                     !surface - atmosphere CH4 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XOXDFS(JY,JX)                     !surface - atmosphere O2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XNGDFS(JY,JX)                     !surface - atmosphere N2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XN2DFS(JY,JX)                     !surface - atmosphere N2O dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8) :: XN3DFS(JY,JX)                     !surface - atmosphere NH3 dissolution (+ve) - volatilization (-ve) non-band, [g d-2 h-1]
  real(r8) :: XNBDFS(JY,JX)                     !surface - atmosphere NH3 dissolution (+ve) - volatilization (-ve) band, [g d-2 h-1]
  real(r8) :: XHGDFS(JY,JX)                     !surface - atmosphere H2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]

end module LandSurfDataType
