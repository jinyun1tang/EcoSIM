module SoilPhysDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  SLOPE_col(:,:,:)                             !slope	in four directions, [o]
  real(r8),target,allocatable ::  FieldCapacity_vr(:,:,:)                      !water contents at field capacity,[m3 d-2]
  real(r8),target,allocatable ::  WiltPoint_vr(:,:,:)                          !water contents at wilting point,[m3 d-2]
  real(r8),target,allocatable ::  SatHydroCondVert_vr(:,:,:)                   !soil vertical saturated hydraulic conductivity, [mm h-1]
  real(r8),target,allocatable ::  SatHydroCondHrzn_vr(:,:,:)                   !soil horizontal saturated hydraulic conductivity, [mm h-1]
  real(r8),target,allocatable ::  PSIAtFldCapacity_col(:,:)                    !water potentials at field capacity, [MPa]
  real(r8),target,allocatable ::  PSIAtWiltPoint_col(:,:)                      !water potentials at wilting point [MPa]
  real(r8),target,allocatable ::  THW_vr(:,:,:)                                !initial soil water content, [m3 m-3]
  real(r8),target,allocatable ::  THI_vr(:,:,:)                                !initial ice content, [m3 m-3]
  REAL(R8),target,allocatable ::  SurfAlbedo_col(:,:)                          !Surface albedo,[-]
  real(r8),target,allocatable ::  LOGPOROS_vr(:,:,:)                           !log soil porosity,	[-]
  real(r8),target,allocatable ::  LOGFldCapacity_vr(:,:,:)                     !log water content at field capacity,[-]
  real(r8),target,allocatable ::  LOGWiltPoint_vr(:,:,:)                       !log water content at wilting point,[-]
  real(r8),target,allocatable ::  PSD_vr(:,:,:)                                !log (soil porosity /water content at field capacity),[-]
  real(r8),target,allocatable ::  FCD_vr(:,:,:)                                !log water content at field capacity,[-]
  real(r8),target,allocatable ::  SRP_vr(:,:,:)                                !shape parameter for water desorption,[-]
  real(r8),target,allocatable ::  FSLOPE_2DH(:,:,:)                            !fraction of slope in 1 and 2,[-]
  REAL(R8),target,allocatable ::  VLMicPt0_col(:,:,:)                          !initial total soil micropore porosity,	[m3 d-2]
  REAL(R8),target,allocatable ::  LOGPSIAtSat(:,:)                             !log water potential at saturation,	[MPa]
  REAL(R8),target,allocatable ::  LOGPSIFLD_col(:,:)                               !log water potential at field capacity,	[-]
  REAL(R8),target,allocatable ::  LOGPSIMN_col(:,:)                                !log water potential at wilting point,[-]
  REAL(R8),target,allocatable ::  LOGPSIMXD_col(:,:)                               !log water potential at field capacity ,[-]
  REAL(R8),target,allocatable ::  LOGPSIMND_col(:,:)                               !log water potential at saturation - log water potential at field capacity,[-]
  real(r8),target,allocatable ::  VHeatCapacitySoilM_vr(:,:,:)                 !soil solid heat capacity, [MPa m-3 K-1]
  real(r8),target,allocatable ::  ActiveLayDepZ_col(:,:)                       !active layer depth of a permafrost soil, [m]
  real(r8),target,allocatable :: CondGasXSurf_col(:,:)                         !gas conductance for soil-atmosphere exchange, [m/h]
!----------------------------------------------------------------------

contains
  subroutine InitSoilPhysData

  implicit none
  allocate(SLOPE_col(0:3,JY,JX));   SLOPE_col=0._r8
  allocate(FieldCapacity_vr(0:JZ,JY,JX));     FieldCapacity_vr=0._r8
  allocate(WiltPoint_vr(0:JZ,JY,JX));     WiltPoint_vr=0._r8
  allocate(SatHydroCondVert_vr(0:JZ,JY,JX));   SatHydroCondVert_vr=0._r8
  allocate(SatHydroCondHrzn_vr(JZ,JY,JX));     SatHydroCondHrzn_vr=0._r8
  allocate(PSIAtFldCapacity_col(JY,JX));       PSIAtFldCapacity_col=0._r8
  allocate(PSIAtWiltPoint_col(JY,JX));       PSIAtWiltPoint_col=0._r8
  allocate(THW_vr(JZ,JY,JX));      THW_vr=0._r8
  allocate(THI_vr(JZ,JY,JX));      THI_vr=0._r8
  allocate(SurfAlbedo_col(JY,JX));        SurfAlbedo_col=0._r8
  allocate(LOGPOROS_vr(0:JZ,JY,JX));    LOGPOROS_vr=0._r8
  allocate(LOGFldCapacity_vr(0:JZ,JY,JX));    LOGFldCapacity_vr=0._r8
  allocate(LOGWiltPoint_vr(0:JZ,JY,JX));    LOGWiltPoint_vr=0._r8
  allocate(PSD_vr(0:JZ,JY,JX));    PSD_vr=0._r8
  allocate(FCD_vr(0:JZ,JY,JX));    FCD_vr=0._r8
  allocate(SRP_vr(0:JZ,JY,JX));    SRP_vr=0._r8
  allocate(FSLOPE_2DH(2,JY,JX));    FSLOPE_2DH=0._r8
  allocate(VLMicPt0_col(0:JZ,JY,JX));  VLMicPt0_col=0._r8
  allocate(LOGPSIAtSat(JY,JX));       LOGPSIAtSat=0._r8
  allocate(LOGPSIFLD_col(JY,JX));       LOGPSIFLD_col=0._r8
  allocate(LOGPSIMN_col(JY,JX));       LOGPSIMN_col=0._r8
  allocate(LOGPSIMXD_col(JY,JX));       LOGPSIMXD_col=0._r8
  allocate(LOGPSIMND_col(JY,JX));       LOGPSIMND_col=0._r8
  allocate(VHeatCapacitySoilM_vr(0:JZ,JY,JX));   VHeatCapacitySoilM_vr=0._r8
  allocate(ActiveLayDepZ_col(JY,JX));       ActiveLayDepZ_col=0._r8
  allocate(CondGasXSurf_col(JY,JX)); CondGasXSurf_col=0._r8
  end subroutine InitSoilPhysData

!----------------------------------------------------------------------
  subroutine DestructSoilPhysData
  use abortutils, only : destroy
  implicit none
  call destroy(CondGasXSurf_col)
  call destroy(SLOPE_col)
  call destroy(FieldCapacity_vr)
  call destroy(WiltPoint_vr)
  call destroy(SatHydroCondVert_vr)
  call destroy(SatHydroCondHrzn_vr)
  call destroy(PSIAtFldCapacity_col)
  call destroy(PSIAtWiltPoint_col)
  call destroy(THW_vr)
  call destroy(THI_vr)
  call destroy(SurfAlbedo_col)
  call destroy(LOGPOROS_vr)
  call destroy(LOGFldCapacity_vr)
  call destroy(LOGWiltPoint_vr)
  call destroy(PSD_vr)
  call destroy(FCD_vr)
  call destroy(SRP_vr)
  call destroy(FSLOPE_2DH)
  call destroy(VLMicPt0_col)
  call destroy(LOGPSIAtSat)
  call destroy(LOGPSIFLD_col)
  call destroy(LOGPSIMN_col)
  call destroy(LOGPSIMXD_col)
  call destroy(LOGPSIMND_col)
  call destroy(VHeatCapacitySoilM_vr)
  call destroy(ActiveLayDepZ_col)
  end subroutine DestructSoilPhysData

end module SoilPhysDataType
