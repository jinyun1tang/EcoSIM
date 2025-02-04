module LandSurfDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod, only : idg_beg,idg_end
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8) :: ALTIG                                                      !Altitude of landscape, [m]
  real(r8),target,allocatable ::  SoilSurfRoughnesst0_col(:,:)           !Initial soil surface roughness height, [m]
  real(r8),target,allocatable ::  ZERO4PlantDisplace_col(:,:)            !Zero plane displacement height, [m]
  real(r8),target,allocatable ::  RoughHeight_col(:,:)                   !Canopy surface roughness height, [m]
  real(r8),target,allocatable ::  SoiSurfRoughness(:,:)                  !Soil surface roughness height for calculating runoff velocity, [m]
  real(r8),target,allocatable ::  WindMesureHeight_col(:,:)              !Wind speed measurement height, [m]
  real(r8),target,allocatable ::  ALT(:,:)                               !Altitude of grid cell, [m]
  real(r8),target,allocatable ::  AbvCanopyBndlResist_col(:,:)           !Isothermal boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  RIB_col(:,:)                           !Richardson number for calculating boundary layer resistance, [-]
  real(r8),target,allocatable ::  ALTI(:,:)                              !Altitude of landscape, [m]
  real(r8),target,allocatable ::  SineGrndSlope_col(:,:)                 !Sine of slope, [-]
  real(r8),target,allocatable ::  CosineGrndSlope_col(:,:)               !Cosine of slope, [-]
  real(r8),target,allocatable ::  GroundSurfAzimuth_col(:,:)             !Azimuth of slope, [-]
  real(r8),target,allocatable ::  ALTZ_col(:,:)                          !Altitude, [m]
  real(r8),target,allocatable ::  SL(:,:)                                !Slope, [o]
  real(r8),target,allocatable ::  ASP_col(:,:)                           !Aspect , [o]
  real(r8),target,allocatable ::  trcg_DisolEvap_Atm2Soil_flx(:,:,:)     !Surface-atmosphere gas exchange through diffusion, >0 into soil [g d-2 h-1]

!----------------------------------------------------------------------

contains
  subroutine InitLandSurfData

  implicit none
  allocate(SoilSurfRoughnesst0_col(JY,JX));          SoilSurfRoughnesst0_col=0._r8
  allocate(ZERO4PlantDisplace_col(JY,JX));          ZERO4PlantDisplace_col=0._r8
  allocate(RoughHeight_col(JY,JX));          RoughHeight_col=0._r8
  allocate(SoiSurfRoughness(JY,JX));          SoiSurfRoughness=0._r8
  allocate(WindMesureHeight_col(JY,JX));          WindMesureHeight_col=0._r8
  allocate(ALT(JY,JX));         ALT=0._r8
  allocate(AbvCanopyBndlResist_col(JY,JX));         AbvCanopyBndlResist_col=0._r8
  allocate(RIB_col(JY,JX));         RIB_col=0._r8
  allocate(ALTI(JY,JX));        ALTI=0._r8
  allocate(SineGrndSlope_col(JY,JX));        SineGrndSlope_col=0._r8
  allocate(CosineGrndSlope_col(JY,JX));        CosineGrndSlope_col=0._r8
  allocate(GroundSurfAzimuth_col(JY,JX));        GroundSurfAzimuth_col=0._r8
  allocate(ALTZ_col(JY,JX));        ALTZ_col=0._r8
  allocate(SL(JY,JX));          SL=0._r8
  allocate(ASP_col(JY,JX));         ASP_col=0._r8

  allocate(trcg_DisolEvap_Atm2Soil_flx(idg_beg:idg_end,JY,JX)); trcg_DisolEvap_Atm2Soil_flx=0._r8
  end subroutine InitLandSurfData

!----------------------------------------------------------------------
  subroutine DestructLandSurfData
  use abortutils, only : destroy
  implicit none
  call destroy(SoilSurfRoughnesst0_col)
  call destroy(ZERO4PlantDisplace_col)
  call destroy(RoughHeight_col)
  call destroy(SoiSurfRoughness)
  call destroy(WindMesureHeight_col)
  call destroy(ALT)
  call destroy(AbvCanopyBndlResist_col)
  call destroy(RIB_col)
  call destroy(ALTI)
  call destroy(SineGrndSlope_col)
  call destroy(CosineGrndSlope_col)
  call destroy(GroundSurfAzimuth_col)
  call destroy(ALTZ_col)
  call destroy(SL)
  call destroy(ASP_col)

  end subroutine DestructLandSurfData

end module LandSurfDataType
