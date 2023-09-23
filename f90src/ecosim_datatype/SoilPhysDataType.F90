module SoilPhysDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none

  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target,allocatable ::  SLOPE(:,:,:)                       !slope	in four directions [o]
  real(r8),target,allocatable ::  FieldCapacity(:,:,:)                          !water contents at field capacity
  real(r8),target,allocatable ::  WiltPoint(:,:,:)                          !water contents at wilting point
  real(r8),target,allocatable ::  SatHydroCondVert(:,:,:)                        !soil vertical saturated hydraulic conductivity [mm h-1]
  real(r8),target,allocatable ::  SatHydroCondHrzn(:,:,:)                        !soil horizontal saturated hydraulic conductivity, [mm h-1]
  real(r8),target,allocatable ::  PSIAtFldCapacity(:,:)                         !water potentials at field capacity, [MPa]
  real(r8),target,allocatable ::  PSIAtWiltPoint(:,:)                         !water potentials at wilting point [MPa]
  real(r8),target,allocatable ::  THW(:,:,:)                         !initial soil water content
  real(r8),target,allocatable ::  THI(:,:,:)                         !initial ice content
  REAL(R8),target,allocatable ::  ALBX(:,:)                          !Surface albedo
  real(r8),target,allocatable ::  LOGPOROS(:,:,:)                         !log soil porosity	-
  real(r8),target,allocatable ::  LOGFldCapacity(:,:,:)                         !log water content at field capacity
  real(r8),target,allocatable ::  LOGWiltPoint(:,:,:)                         !log water content at wilting point
  real(r8),target,allocatable ::  PSD(:,:,:)                         !log soil porosity - log water content at field capacity
  real(r8),target,allocatable ::  FCD(:,:,:)                         !log water content at field capacity
  real(r8),target,allocatable ::  SRP(:,:,:)                         !shape parameter for water desorption
  real(r8),target,allocatable ::  FSLOPE(:,:,:)                      !fraction of slope in 1 and 2
  REAL(R8),target,allocatable ::  VLMicPt0(:,:,:)                       !initial total soil micropore porosity	m3 d-2
  REAL(R8),target,allocatable ::  LOGPSIAtSat(:,:)                         !log water potential at saturation	MPa
  REAL(R8),target,allocatable ::  LOGPSIFLD(:,:)                         !log water potential at field capacity	-
  REAL(R8),target,allocatable ::  LOGPSIMN(:,:)                         !log water potential at wilting point
  REAL(R8),target,allocatable ::  LOGPSIMXD(:,:)                         !log water potential at field capacity 	-
  REAL(R8),target,allocatable ::  LOGPSIMND(:,:)                         !log water potential at saturation - log water potential at field capacity
  real(r8),target,allocatable ::  VHeatCapacitySoilM(:,:,:)                        !soil solid heat capacity [MPa m-3 K-1]
  real(r8),target,allocatable ::  ActiveLayDepth(:,:)                         !active layer depth, [m]
!----------------------------------------------------------------------

contains
  subroutine InitSoilPhysData

  implicit none
  allocate(SLOPE(0:3,JY,JX));   SLOPE=0._r8
  allocate(FieldCapacity(0:JZ,JY,JX));     FieldCapacity=0._r8
  allocate(WiltPoint(0:JZ,JY,JX));     WiltPoint=0._r8
  allocate(SatHydroCondVert(0:JZ,JY,JX));   SatHydroCondVert=0._r8
  allocate(SatHydroCondHrzn(JZ,JY,JX));     SatHydroCondHrzn=0._r8
  allocate(PSIAtFldCapacity(JY,JX));       PSIAtFldCapacity=0._r8
  allocate(PSIAtWiltPoint(JY,JX));       PSIAtWiltPoint=0._r8
  allocate(THW(JZ,JY,JX));      THW=0._r8
  allocate(THI(JZ,JY,JX));      THI=0._r8
  allocate(ALBX(JY,JX));        ALBX=0._r8
  allocate(LOGPOROS(0:JZ,JY,JX));    LOGPOROS=0._r8
  allocate(LOGFldCapacity(0:JZ,JY,JX));    LOGFldCapacity=0._r8
  allocate(LOGWiltPoint(0:JZ,JY,JX));    LOGWiltPoint=0._r8
  allocate(PSD(0:JZ,JY,JX));    PSD=0._r8
  allocate(FCD(0:JZ,JY,JX));    FCD=0._r8
  allocate(SRP(0:JZ,JY,JX));    SRP=0._r8
  allocate(FSLOPE(2,JY,JX));    FSLOPE=0._r8
  allocate(VLMicPt0(0:JZ,JY,JX));  VLMicPt0=0._r8
  allocate(LOGPSIAtSat(JY,JX));       LOGPSIAtSat=0._r8
  allocate(LOGPSIFLD(JY,JX));       LOGPSIFLD=0._r8
  allocate(LOGPSIMN(JY,JX));       LOGPSIMN=0._r8
  allocate(LOGPSIMXD(JY,JX));       LOGPSIMXD=0._r8
  allocate(LOGPSIMND(JY,JX));       LOGPSIMND=0._r8
  allocate(VHeatCapacitySoilM(0:JZ,JY,JX));   VHeatCapacitySoilM=0._r8
  allocate(ActiveLayDepth(JY,JX));       ActiveLayDepth=0._r8
  end subroutine InitSoilPhysData

!----------------------------------------------------------------------
  subroutine DestructSoilPhysData
  use abortutils, only : destroy
  implicit none
  call destroy(SLOPE)
  call destroy(FieldCapacity)
  call destroy(WiltPoint)
  call destroy(SatHydroCondVert)
  call destroy(SatHydroCondHrzn)
  call destroy(PSIAtFldCapacity)
  call destroy(PSIAtWiltPoint)
  call destroy(THW)
  call destroy(THI)
  call destroy(ALBX)
  call destroy(LOGPOROS)
  call destroy(LOGFldCapacity)
  call destroy(LOGWiltPoint)
  call destroy(PSD)
  call destroy(FCD)
  call destroy(SRP)
  call destroy(FSLOPE)
  call destroy(VLMicPt0)
  call destroy(LOGPSIAtSat)
  call destroy(LOGPSIFLD)
  call destroy(LOGPSIMN)
  call destroy(LOGPSIMXD)
  call destroy(LOGPSIMND)
  call destroy(VHeatCapacitySoilM)
  call destroy(ActiveLayDepth)
  end subroutine DestructSoilPhysData

end module SoilPhysDataType
