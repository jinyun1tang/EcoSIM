module SoilHeatDatatype


  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target,allocatable ::  TKSZ(:,:,:)                        !
  real(r8),target,allocatable ::  TKS(:,:,:)                         !
  real(r8),target,allocatable ::  THAW(:,:,:)                        !hourly accumulated freeze-thaw flux in micropores
  real(r8),target,allocatable ::  HTHAW(:,:,:)                       !hourly accumulated freeze-thaw latent heat flux
  real(r8),target,allocatable ::  THAWH(:,:,:)                       !hourly accumulated freeze-thaw flux in macropores
  real(r8),target,allocatable ::  XTHAWW(:,:,:)                      !hourly accumulated latent heat flux from freeze-thaw
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
  allocate(THAW(JZ,JY,JX));     THAW=0._r8
  allocate(HTHAW(JZ,JY,JX));    HTHAW=0._r8
  allocate(THAWH(JZ,JY,JX));    THAWH=0._r8
  allocate(XTHAWW(JS,JY,JX));   XTHAWW=0._r8
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
  call destroy(THAW)
  call destroy(HTHAW)
  call destroy(THAWH)
  call destroy(XTHAWW)
  call destroy(VHeatCapacity)
  call destroy(TCS)
  call destroy(STC)
  call destroy(DTC)
  end subroutine DestructSoilHeatData

end module SoilHeatDatatype
