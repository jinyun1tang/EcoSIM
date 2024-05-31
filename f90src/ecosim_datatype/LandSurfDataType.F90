module LandSurfDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod, only : idg_beg,idg_end
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8) :: ALTIG                             !altitude of landscape, [m]

  real(r8),target,allocatable ::  SoiSurfRoughnesst0(:,:)                            !initial soil surface roughness height, [m]
  real(r8),target,allocatable ::  ZERO4Groth_pftlanDisp(:,:)                            !zero plane displacement height, [m]
  real(r8),target,allocatable ::  RoughHeight(:,:)                   !canopy surface roughness height, [m]
  real(r8),target,allocatable ::  SoiSurfRoughness(:,:)                            ! soil surface roughness height for calculating runoff velocity, [m]
  real(r8),target,allocatable ::  WindMesHeight(:,:)                            !wind speed measurement height, [m]
  real(r8),target,allocatable ::  ALT(:,:)                           !altitude of grid cell, [m]
  real(r8),target,allocatable ::  BndlResistAboveCanG(:,:)                           !isothermal boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  RIB(:,:)                           !Richardson number for calculating boundary layer resistance, [-]
  real(r8),target,allocatable ::  ALTI(:,:)                          !altitude of landscape, [m]
  real(r8),target,allocatable ::  SineGrndSlope_col(:,:)                          !sine of slope, [-]
  real(r8),target,allocatable ::  CosineGrndSlope_col(:,:)                          !cosine of slope, [-]
  real(r8),target,allocatable ::  GroundSurfAzimuth_col(:,:)                          !azimuth of slope, [-]
  real(r8),target,allocatable ::  ALTZ(:,:)                          !altitude, [m]
  real(r8),target,allocatable ::  SL(:,:)                            !slope, [o]
  real(r8),target,allocatable ::  ASP(:,:)                           !aspect , [o]
  real(r8),target,allocatable :: GasSfAtmFlx_col(:,:,:)   ! surface-atmosphere gas exchange flux , [g d-2 h-1]
  real(r8),target,allocatable ::  XCODFS(:,:)                        !surface - atmosphere CO2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XCHDFS(:,:)                        !surface - atmosphere CH4 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XOXDFS(:,:)                        !surface - atmosphere O2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XNGDFS(:,:)                        !surface - atmosphere N2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XN2DFS(:,:)                        !surface - atmosphere N2O dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
  real(r8),target,allocatable ::  XN3DFS(:,:)                        !surface - atmosphere NH3 dissolution (+ve) - volatilization (-ve) non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  XNBDFS(:,:)                        !surface - atmosphere NH3 dissolution (+ve) - volatilization (-ve) band, [g d-2 h-1]
  real(r8),target,allocatable ::  XHGDFS(:,:)                        !surface - atmosphere H2 dissolution (+ve) - volatilization (-ve) , [g d-2 h-1]
!----------------------------------------------------------------------

contains
  subroutine InitLandSurfData

  implicit none
  allocate(SoiSurfRoughnesst0(JY,JX));          SoiSurfRoughnesst0=0._r8
  allocate(ZERO4Groth_pftlanDisp(JY,JX));          ZERO4Groth_pftlanDisp=0._r8
  allocate(RoughHeight(JY,JX));          RoughHeight=0._r8
  allocate(SoiSurfRoughness(JY,JX));          SoiSurfRoughness=0._r8
  allocate(WindMesHeight(JY,JX));          WindMesHeight=0._r8
  allocate(ALT(JY,JX));         ALT=0._r8
  allocate(BndlResistAboveCanG(JY,JX));         BndlResistAboveCanG=0._r8
  allocate(RIB(JY,JX));         RIB=0._r8
  allocate(ALTI(JY,JX));        ALTI=0._r8
  allocate(SineGrndSlope_col(JY,JX));        SineGrndSlope_col=0._r8
  allocate(CosineGrndSlope_col(JY,JX));        CosineGrndSlope_col=0._r8
  allocate(GroundSurfAzimuth_col(JY,JX));        GroundSurfAzimuth_col=0._r8
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

  allocate(GasSfAtmFlx_col(idg_beg:idg_end,JY,JX)); GasSfAtmFlx_col=0._r8
  end subroutine InitLandSurfData

!----------------------------------------------------------------------
  subroutine DestructLandSurfData
  use abortutils, only : destroy
  implicit none
  call destroy(SoiSurfRoughnesst0)
  call destroy(ZERO4Groth_pftlanDisp)
  call destroy(RoughHeight)
  call destroy(SoiSurfRoughness)
  call destroy(WindMesHeight)
  call destroy(ALT)
  call destroy(BndlResistAboveCanG)
  call destroy(RIB)
  call destroy(ALTI)
  call destroy(SineGrndSlope_col)
  call destroy(CosineGrndSlope_col)
  call destroy(GroundSurfAzimuth_col)
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
