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

  real(r8),target,allocatable ::  SoilSurfRoughnesst0_col(:,:)                            !initial soil surface roughness height, [m]
  real(r8),target,allocatable ::  ZERO4PlantDisplace_col(:,:)                            !zero plane displacement height, [m]
  real(r8),target,allocatable ::  RoughHeight_col(:,:)                   !canopy surface roughness height, [m]
  real(r8),target,allocatable ::  SoiSurfRoughness(:,:)              ! soil surface roughness height for calculating runoff velocity, [m]
  real(r8),target,allocatable ::  WindMesHeight(:,:)                 !wind speed measurement height, [m]
  real(r8),target,allocatable ::  ALT(:,:)                           !altitude of grid cell, [m]
  real(r8),target,allocatable ::  AbvCanopyBndlResist_col(:,:)                           !isothermal boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  RIB(:,:)                           !Richardson number for calculating boundary layer resistance, [-]
  real(r8),target,allocatable ::  ALTI(:,:)                          !altitude of landscape, [m]
  real(r8),target,allocatable ::  SineGrndSlope_col(:,:)                          !sine of slope, [-]
  real(r8),target,allocatable ::  CosineGrndSlope_col(:,:)                          !cosine of slope, [-]
  real(r8),target,allocatable ::  GroundSurfAzimuth_col(:,:)                          !azimuth of slope, [-]
  real(r8),target,allocatable ::  ALTZ(:,:)                          !altitude, [m]
  real(r8),target,allocatable ::  SL(:,:)                            !slope, [o]
  real(r8),target,allocatable ::  ASP_col(:,:)                           !aspect , [o]
  real(r8),target,allocatable ::  Gas_Flx_atmDif2soil_col(:,:,:)              ! surface-atmosphere gas exchange flux , >0 into soil [g d-2 h-1]

!----------------------------------------------------------------------

contains
  subroutine InitLandSurfData

  implicit none
  allocate(SoilSurfRoughnesst0_col(JY,JX));          SoilSurfRoughnesst0_col=0._r8
  allocate(ZERO4PlantDisplace_col(JY,JX));          ZERO4PlantDisplace_col=0._r8
  allocate(RoughHeight_col(JY,JX));          RoughHeight_col=0._r8
  allocate(SoiSurfRoughness(JY,JX));          SoiSurfRoughness=0._r8
  allocate(WindMesHeight(JY,JX));          WindMesHeight=0._r8
  allocate(ALT(JY,JX));         ALT=0._r8
  allocate(AbvCanopyBndlResist_col(JY,JX));         AbvCanopyBndlResist_col=0._r8
  allocate(RIB(JY,JX));         RIB=0._r8
  allocate(ALTI(JY,JX));        ALTI=0._r8
  allocate(SineGrndSlope_col(JY,JX));        SineGrndSlope_col=0._r8
  allocate(CosineGrndSlope_col(JY,JX));        CosineGrndSlope_col=0._r8
  allocate(GroundSurfAzimuth_col(JY,JX));        GroundSurfAzimuth_col=0._r8
  allocate(ALTZ(JY,JX));        ALTZ=0._r8
  allocate(SL(JY,JX));          SL=0._r8
  allocate(ASP_col(JY,JX));         ASP_col=0._r8

  allocate(Gas_Flx_atmDif2soil_col(idg_beg:idg_end,JY,JX)); Gas_Flx_atmDif2soil_col=0._r8
  end subroutine InitLandSurfData

!----------------------------------------------------------------------
  subroutine DestructLandSurfData
  use abortutils, only : destroy
  implicit none
  call destroy(SoilSurfRoughnesst0_col)
  call destroy(ZERO4PlantDisplace_col)
  call destroy(RoughHeight_col)
  call destroy(SoiSurfRoughness)
  call destroy(WindMesHeight)
  call destroy(ALT)
  call destroy(AbvCanopyBndlResist_col)
  call destroy(RIB)
  call destroy(ALTI)
  call destroy(SineGrndSlope_col)
  call destroy(CosineGrndSlope_col)
  call destroy(GroundSurfAzimuth_col)
  call destroy(ALTZ)
  call destroy(SL)
  call destroy(ASP_col)

  end subroutine DestructLandSurfData

end module LandSurfDataType
