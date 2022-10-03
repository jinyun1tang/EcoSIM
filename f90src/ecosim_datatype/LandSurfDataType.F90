module LandSurfDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__
  real(r8) :: ALTIG                             !altitude of landscape, [m]

  real(r8),allocatable ::  ZS(:,:)                            !initial soil surface roughness height, [m]
  real(r8),allocatable ::  ZD(:,:)                            !zero plane displacement height, [m]
  real(r8),allocatable ::  ZR(:,:)                            !canopy surface roughness height, [m]
  real(r8),allocatable ::  ZM(:,:)                            ! soil surface roughness height for calculating runoff velocity, [m]
  real(r8),allocatable ::  Z0(:,:)                            !wind speed measurement height, [m]
  real(r8),allocatable ::  ALT(:,:)                           !altitude of grid cell, [m]
  real(r8),allocatable ::  RAB(:,:)                           !isothermal boundary layer resistance, [h m-1]
  real(r8),allocatable ::  RIB(:,:)                           !Richardson number for calculating boundary layer resistance, [-]
  real(r8),allocatable ::  ALTI(:,:)                          !altitude of landscape, [m]
  real(r8),allocatable ::  GSIN(:,:)                          !sine of slope, [-]
  real(r8),allocatable ::  GCOS(:,:)                          !cosine of slope, [-]
  real(r8),allocatable ::  GAZI(:,:)                          !azimuth of slope, [-]
  real(r8),allocatable ::  ALTZ(:,:)                          !altitude, [m]
  real(r8),allocatable ::  SL(:,:)                            !slope, [o]
  real(r8),allocatable ::  ASP(:,:)                           !aspect , [o]
  real(r8),allocatable ::  XCODFS(:,:)                        !surface - atmosphere CO2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XCHDFS(:,:)                        !surface - atmosphere CH4 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XOXDFS(:,:)                        !surface - atmosphere O2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XNGDFS(:,:)                        !surface - atmosphere N2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XN2DFS(:,:)                        !surface - atmosphere N2O dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),allocatable ::  XN3DFS(:,:)                        !surface - atmosphere NH3 dissolution (+ve) - volatilization (-ve) non-band, [g d-2 h-1]
  real(r8),allocatable ::  XNBDFS(:,:)                        !surface - atmosphere NH3 dissolution (+ve) - volatilization (-ve) band, [g d-2 h-1]
  real(r8),allocatable ::  XHGDFS(:,:)                        !surface - atmosphere H2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
!----------------------------------------------------------------------

contains
  subroutine InitLandSurfData

  implicit none
  allocate(ZS(JY,JX));          ZS=0._r8
  allocate(ZD(JY,JX));          ZD=0._r8
  allocate(ZR(JY,JX));          ZR=0._r8
  allocate(ZM(JY,JX));          ZM=0._r8
  allocate(Z0(JY,JX));          Z0=0._r8
  allocate(ALT(JY,JX));         ALT=0._r8
  allocate(RAB(JY,JX));         RAB=0._r8
  allocate(RIB(JY,JX));         RIB=0._r8
  allocate(ALTI(JY,JX));        ALTI=0._r8
  allocate(GSIN(JY,JX));        GSIN=0._r8
  allocate(GCOS(JY,JX));        GCOS=0._r8
  allocate(GAZI(JY,JX));        GAZI=0._r8
  allocate(ALTZ(JY,JX));        ALTZ=0._r8
  allocate(SL(JY,JX));          SL=0._r8
  allocate(ASP(JY,JX));         ASP=0._r8
  allocate(XCODFS(JY,JX));      XCODFS=0._r8
  allocate(XCHDFS(JY,JX));      XCHDFS=0._r8
  allocate(XOXDFS(JY,JX));      XOXDFS=0._r8
  allocate(XNGDFS(JY,JX));      XNGDFS=0._r8
  allocate(XN2DFS(JY,JX));      XN2DFS=0._r8
  allocate(XN3DFS(JY,JX));      XN3DFS=0._r8
  allocate(XNBDFS(JY,JX));      XNBDFS=0._r8
  allocate(XHGDFS(JY,JX));      XHGDFS=0._r8
  end subroutine InitLandSurfData

!----------------------------------------------------------------------
  subroutine DestructLandSurfData
  use abortutils, only : destroy
  implicit none
  call destroy(ZS)
  call destroy(ZD)
  call destroy(ZR)
  call destroy(ZM)
  call destroy(Z0)
  call destroy(ALT)
  call destroy(RAB)
  call destroy(RIB)
  call destroy(ALTI)
  call destroy(GSIN)
  call destroy(GCOS)
  call destroy(GAZI)
  call destroy(ALTZ)
  call destroy(SL)
  call destroy(ASP)
  call destroy(XCODFS)
  call destroy(XCHDFS)
  call destroy(XOXDFS)
  call destroy(XNGDFS)
  call destroy(XN2DFS)
  call destroy(XN3DFS)
  call destroy(XNBDFS)
  call destroy(XHGDFS)
  end subroutine DestructLandSurfData

end module LandSurfDataType
