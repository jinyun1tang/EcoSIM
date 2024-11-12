module SoilPhysDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  SLOPE(:,:,:)                              !slope	in four directions [o]
  real(r8),target,allocatable ::  FieldCapacity_vr(:,:,:)                      !water contents at field capacity
  real(r8),target,allocatable ::  WiltPoint_vr(:,:,:)                          !water contents at wilting point
  real(r8),target,allocatable ::  SatHydroCondVert_vr(:,:,:)                   !soil vertical saturated hydraulic conductivity [mm h-1]
  real(r8),target,allocatable ::  SatHydroCondHrzn_vr(:,:,:)                   !soil horizontal saturated hydraulic conductivity, [mm h-1]
  real(r8),target,allocatable ::  PSIAtFldCapacity(:,:)                     !water potentials at field capacity, [MPa]
  real(r8),target,allocatable ::  PSIAtWiltPoint(:,:)                       !water potentials at wilting point [MPa]
  real(r8),target,allocatable ::  THW(:,:,:)                         !initial soil water content
  real(r8),target,allocatable ::  THI(:,:,:)                         !initial ice content
  REAL(R8),target,allocatable ::  SurfAlbedo_col(:,:)                          !Surface albedo
  real(r8),target,allocatable ::  LOGPOROS_vr(:,:,:)                         !log soil porosity	-
  real(r8),target,allocatable ::  LOGFldCapacity_vr(:,:,:)                         !log water content at field capacity
  real(r8),target,allocatable ::  LOGWiltPoint_vr(:,:,:)                         !log water content at wilting point
  real(r8),target,allocatable ::  PSD(:,:,:)                         !log soil porosity - log water content at field capacity
  real(r8),target,allocatable ::  FCD(:,:,:)                         !log water content at field capacity
  real(r8),target,allocatable ::  SRP(:,:,:)                         !shape parameter for water desorption
  real(r8),target,allocatable ::  FSLOPE(:,:,:)                      !fraction of slope in 1 and 2
  REAL(R8),target,allocatable ::  VLMicPt0_col(:,:,:)                       !initial total soil micropore porosity	m3 d-2
  REAL(R8),target,allocatable ::  LOGPSIAtSat(:,:)                         !log water potential at saturation	MPa
  REAL(R8),target,allocatable ::  LOGPSIFLD(:,:)                         !log water potential at field capacity	-
  REAL(R8),target,allocatable ::  LOGPSIMN(:,:)                         !log water potential at wilting point
  REAL(R8),target,allocatable ::  LOGPSIMXD(:,:)                         !log water potential at field capacity 	-
  REAL(R8),target,allocatable ::  LOGPSIMND(:,:)                         !log water potential at saturation - log water potential at field capacity
  real(r8),target,allocatable ::  VHeatCapacitySoilM_vr(:,:,:)                        !soil solid heat capacity [MPa m-3 K-1]
  real(r8),target,allocatable ::  ActiveLayDepth(:,:)                         !active layer depth, [m]
!----------------------------------------------------------------------

contains
  subroutine InitSoilPhysData

  implicit none
  allocate(SLOPE(0:3,JY,JX));   SLOPE=0._r8
  allocate(FieldCapacity_vr(0:JZ,JY,JX));     FieldCapacity_vr=0._r8
  allocate(WiltPoint_vr(0:JZ,JY,JX));     WiltPoint_vr=0._r8
  allocate(SatHydroCondVert_vr(0:JZ,JY,JX));   SatHydroCondVert_vr=0._r8
  allocate(SatHydroCondHrzn_vr(JZ,JY,JX));     SatHydroCondHrzn_vr=0._r8
  allocate(PSIAtFldCapacity(JY,JX));       PSIAtFldCapacity=0._r8
  allocate(PSIAtWiltPoint(JY,JX));       PSIAtWiltPoint=0._r8
  allocate(THW(JZ,JY,JX));      THW=0._r8
  allocate(THI(JZ,JY,JX));      THI=0._r8
  allocate(SurfAlbedo_col(JY,JX));        SurfAlbedo_col=0._r8
  allocate(LOGPOROS_vr(0:JZ,JY,JX));    LOGPOROS_vr=0._r8
  allocate(LOGFldCapacity_vr(0:JZ,JY,JX));    LOGFldCapacity_vr=0._r8
  allocate(LOGWiltPoint_vr(0:JZ,JY,JX));    LOGWiltPoint_vr=0._r8
  allocate(PSD(0:JZ,JY,JX));    PSD=0._r8
  allocate(FCD(0:JZ,JY,JX));    FCD=0._r8
  allocate(SRP(0:JZ,JY,JX));    SRP=0._r8
  allocate(FSLOPE(2,JY,JX));    FSLOPE=0._r8
  allocate(VLMicPt0_col(0:JZ,JY,JX));  VLMicPt0_col=0._r8
  allocate(LOGPSIAtSat(JY,JX));       LOGPSIAtSat=0._r8
  allocate(LOGPSIFLD(JY,JX));       LOGPSIFLD=0._r8
  allocate(LOGPSIMN(JY,JX));       LOGPSIMN=0._r8
  allocate(LOGPSIMXD(JY,JX));       LOGPSIMXD=0._r8
  allocate(LOGPSIMND(JY,JX));       LOGPSIMND=0._r8
  allocate(VHeatCapacitySoilM_vr(0:JZ,JY,JX));   VHeatCapacitySoilM_vr=0._r8
  allocate(ActiveLayDepth(JY,JX));       ActiveLayDepth=0._r8
  end subroutine InitSoilPhysData

!----------------------------------------------------------------------
  subroutine DestructSoilPhysData
  use abortutils, only : destroy
  implicit none
  call destroy(SLOPE)
  call destroy(FieldCapacity_vr)
  call destroy(WiltPoint_vr)
  call destroy(SatHydroCondVert_vr)
  call destroy(SatHydroCondHrzn_vr)
  call destroy(PSIAtFldCapacity)
  call destroy(PSIAtWiltPoint)
  call destroy(THW)
  call destroy(THI)
  call destroy(SurfAlbedo_col)
  call destroy(LOGPOROS_vr)
  call destroy(LOGFldCapacity_vr)
  call destroy(LOGWiltPoint_vr)
  call destroy(PSD)
  call destroy(FCD)
  call destroy(SRP)
  call destroy(FSLOPE)
  call destroy(VLMicPt0_col)
  call destroy(LOGPSIAtSat)
  call destroy(LOGPSIFLD)
  call destroy(LOGPSIMN)
  call destroy(LOGPSIMXD)
  call destroy(LOGPSIMND)
  call destroy(VHeatCapacitySoilM_vr)
  call destroy(ActiveLayDepth)
  end subroutine DestructSoilPhysData

end module SoilPhysDataType
