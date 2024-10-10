module SoilHeatDatatype


  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  TKSZ(:,:,:)                        !
  real(r8),target,allocatable ::  TKS_vr(:,:,:)                         !
  real(r8),target,allocatable ::  TLIceThawMicP(:,:,:)               !hourly accumulated freeze-thaw flux in micropores
  real(r8),target,allocatable ::  TLPhaseChangeHeat2Soi(:,:,:)        !hourly accumulated freeze-thaw latent heat flux
  real(r8),target,allocatable ::  TLIceThawMacP(:,:,:)               !hourly accumulated freeze-thaw flux in macropores
  real(r8),target,allocatable ::  XPhaseChangeHeatL_snvr(:,:,:)           !hourly accumulated latent heat flux from freeze-thaw
  real(r8),target,allocatable ::  VHeatCapacity_vr(:,:,:)               !soil heat capacity [MJ m-3 K-1]
  real(r8),target,allocatable ::  TCS(:,:,:)                         !soil temperature [oC]
  real(r8),target,allocatable ::  HeatStore_col(:,:)                 !heat stored over the grid MJ d-2, including soil, litter and canopy
  real(r8),target,allocatable ::  NumerSolidThermCond(:,:,:)         !numerator for soil solid thermal conductivity [MJ m h-1 K-1]
  real(r8),target,allocatable ::  DenomSolidThermCond(:,:,:)         !denominator for soil solid thermal conductivity
  real(r8),target,allocatable ::  HeatFlx2Grnd_col(:,:)              !heat flux into ground, computed from surface energy balance model 
  real(r8),target,allocatable ::  THeatFlow2Soil_vr(:,:,:)              !hourly heat flux into soil layer  [MJ m-3]
  real(r8),target,allocatable ::  tHeatUptk_col(:,:)                 !Heat utake by plant through transpiration [MJ/d2/h] 
  real(r8),target,allocatable ::  HeatDrain_col(:,:)                 !heat loss through drainage [MJ/d2/h]
  real(r8),target,allocatable ::  HeatRunSurf_col(:,:)               !heat loss through surface runoff [MJ/d2/h]
!----------------------------------------------------------------------

contains
  subroutine InitSoilHeatData

  implicit none
  allocate(TKSZ(366,24,JZ));    TKSZ=0._r8
  allocate(TKS_vr(0:JZ,JY,JX));    TKS_vr=0._r8
  allocate(HeatFlx2Grnd_col(JY,JX));   HeatFlx2Grnd_col=0._r8
  allocate(HeatStore_col(JY,JX));  HeatStore_col=0._r8
  allocate(TLIceThawMicP(JZ,JY,JX));     TLIceThawMicP=0._r8
  allocate(TLPhaseChangeHeat2Soi(JZ,JY,JX));    TLPhaseChangeHeat2Soi=0._r8
  allocate(TLIceThawMacP(JZ,JY,JX));    TLIceThawMacP=0._r8
  allocate(XPhaseChangeHeatL_snvr(JS,JY,JX));   XPhaseChangeHeatL_snvr=0._r8
  allocate(VHeatCapacity_vr(0:JZ,JY,JX));   VHeatCapacity_vr=0._r8
  allocate(TCS(0:JZ,JY,JX));    TCS=0._r8
  allocate(NumerSolidThermCond(JZ,JY,JX));      NumerSolidThermCond=0._r8
  allocate(DenomSolidThermCond(JZ,JY,JX));      DenomSolidThermCond=0._r8
  allocate(THeatFlow2Soil_vr(JZ,JY,JX));    THeatFlow2Soil_vr=0._r8  
  allocate(tHeatUptk_col(JY,JX));   tHeatUptk_col=0._r8
  allocate(HeatDrain_col(JY,JX));  HeatDrain_col = 0._r8
  allocate(HeatRunSurf_col(JY,JX)); HeatRunSurf_col=0._r8
  end subroutine InitSoilHeatData

!----------------------------------------------------------------------
  subroutine DestructSoilHeatData
  use abortutils, only : destroy
  implicit none
  call destroy(TKSZ)
  call destroy(TKS_vr)
  call destroy(TLIceThawMicP)
  call destroy(TLPhaseChangeHeat2Soi)
  call destroy(TLIceThawMacP)
  call destroy(HeatStore_col)
  call destroy(XPhaseChangeHeatL_snvr)
  call destroy(VHeatCapacity_vr)
  call destroy(TCS)
  call destroy(NumerSolidThermCond)
  call destroy(DenomSolidThermCond)
  call destroy(THeatFlow2Soil_vr)  
  call destroy(tHeatUptk_col)
  call destroy(HeatDrain_col)
  call destroy(HeatRunSurf_col)
  end subroutine DestructSoilHeatData

end module SoilHeatDatatype
