module SoilHeatDatatype


  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  TKSZ(:,:,:)                        !
  real(r8),target,allocatable ::  TKS(:,:,:)                         !
  real(r8),target,allocatable ::  TLIceThawMicP(:,:,:)               !hourly accumulated freeze-thaw flux in micropores
  real(r8),target,allocatable ::  TLPhaseChangeHeat2Soi(:,:,:)        !hourly accumulated freeze-thaw latent heat flux
  real(r8),target,allocatable ::  TLIceThawMacP(:,:,:)               !hourly accumulated freeze-thaw flux in macropores
  real(r8),target,allocatable ::  XPhaseChangeHeatL(:,:,:)                      !hourly accumulated latent heat flux from freeze-thaw
  real(r8),target,allocatable ::  VHeatCapacity(:,:,:)                        !soil heat capacity [MJ m-3 K-1]
  real(r8),target,allocatable ::  TCS(:,:,:)                         !soil temperature [oC]
  real(r8),target,allocatable ::  STC(:,:,:)                         !numerator for soil solid thermal conductivity [MJ m h-1 K-1]
  real(r8),target,allocatable ::  DTC(:,:,:)                         !denominator for soil solid thermal conductivity
!----------------------------------------------------------------------

contains
  subroutine InitSoilHeatData

  implicit none
  allocate(TKSZ(366,24,JZ));    TKSZ=0._r8
  allocate(TKS(0:JZ,JY,JX));    TKS=0._r8
  allocate(TLIceThawMicP(JZ,JY,JX));     TLIceThawMicP=0._r8
  allocate(TLPhaseChangeHeat2Soi(JZ,JY,JX));    TLPhaseChangeHeat2Soi=0._r8
  allocate(TLIceThawMacP(JZ,JY,JX));    TLIceThawMacP=0._r8
  allocate(XPhaseChangeHeatL(JS,JY,JX));   XPhaseChangeHeatL=0._r8
  allocate(VHeatCapacity(0:JZ,JY,JX));   VHeatCapacity=0._r8
  allocate(TCS(0:JZ,JY,JX));    TCS=0._r8
  allocate(STC(JZ,JY,JX));      STC=0._r8
  allocate(DTC(JZ,JY,JX));      DTC=0._r8
  end subroutine InitSoilHeatData

!----------------------------------------------------------------------
  subroutine DestructSoilHeatData
  use abortutils, only : destroy
  implicit none
  call destroy(TKSZ)
  call destroy(TKS)
  call destroy(TLIceThawMicP)
  call destroy(TLPhaseChangeHeat2Soi)
  call destroy(TLIceThawMacP)
  call destroy(XPhaseChangeHeatL)
  call destroy(VHeatCapacity)
  call destroy(TCS)
  call destroy(STC)
  call destroy(DTC)
  end subroutine DestructSoilHeatData

end module SoilHeatDatatype
