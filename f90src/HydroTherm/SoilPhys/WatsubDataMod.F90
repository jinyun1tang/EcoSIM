module WatsubDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  TMLiceThawMicP(:,:,:)              !micropore layer integrated ice mass loss due to thaw
  real(r8),allocatable ::  TMLiceThawMacP(:,:,:)              !macropore layer integrated ice mass loss due to thaw

  real(r8),allocatable ::  AREAU(:,:,:)                       !
  real(r8),allocatable ::  AreaUnderWaterTBL(:,:,:)                      !

  real(r8),allocatable ::  VLairMacP(:,:,:)                     !
  real(r8),allocatable ::  TTFLXL(:,:,:)                      !

  real(r8),allocatable ::  FWatExMacP2MicP(:,:,:)             !pressure-driven water flow from macpore to micpore

  real(r8),allocatable ::  TWatCharge2MicP(:,:,:)                       !
  real(r8),allocatable ::  TConvectWaterFlowMacP(:,:,:)                      !
  real(r8),allocatable ::  THeatFlowi(:,:,:)                      !
  real(r8),allocatable ::  FWatFreezeThawMicP(:,:,:)                       !
  real(r8),allocatable ::  TFLXL(:,:,:)                       !
  real(r8),allocatable ::  AVCNHL(:,:,:,:)                    !

  real(r8),allocatable ::  TWatXChange2WatTableX(:,:,:)                      !
  real(r8),allocatable ::  FWatIrrigate2MicP1(:,:,:)                        !
  real(r8),allocatable ::  HeatIrrigation1(:,:,:)                      !

  real(r8),allocatable ::  FWatFreezeThawMacP(:,:,:)                      !

  real(r8),allocatable ::  HydroCondMacP1(:,:,:)                       !
  real(r8),allocatable ::  VLMicP1(:,:,:)                       !
  real(r8),allocatable ::  VLMacP1(:,:,:)                      !

  real(r8),allocatable ::  FLWNX(:,:)                         !
  real(r8),allocatable ::  FLWXNX(:,:)                        !
  real(r8),allocatable ::  FLWHNX(:,:)                        !
  real(r8),allocatable ::  HFLWNX(:,:)                        !
  real(r8),allocatable ::  PSISoilMatricPtmp(:,:,:)                      !


  integer, allocatable ::  N6X(:,:)

!----------------------------------------------------------------------
  public :: InitWatsubData
contains
  subroutine InitWatSubData

  implicit none

  allocate(N6X(JY,JX));         N6X=0

  allocate(TMLiceThawMicP(JZ,JY,JX));   TMLiceThawMicP=0._r8

  allocate(TMLiceThawMacP(JZ,JY,JX));   TMLiceThawMacP=0._r8

  allocate(AREAU(JZ,JY,JX));    AREAU=0._r8
  allocate(AreaUnderWaterTBL(JZ,JY,JX));   AreaUnderWaterTBL=0._r8


  allocate(VLairMacP(JZ,JY,JX));  VLairMacP=0._r8
  allocate(TTFLXL(JZ,JY,JX));   TTFLXL=0._r8

  allocate(FWatExMacP2MicP(JZ,JY,JX));    FWatExMacP2MicP=0._r8

  allocate(TWatCharge2MicP(JZ,JY,JX));    TWatCharge2MicP=0._r8
  allocate(TConvectWaterFlowMacP(JZ,JY,JX));   TConvectWaterFlowMacP=0._r8
  allocate(THeatFlowi(JZ,JY,JX));   THeatFlowi=0._r8
  allocate(FWatFreezeThawMicP(JZ,JY,JX));    FWatFreezeThawMicP=0._r8
  allocate(TFLXL(JZ,JY,JX));    TFLXL=0._r8
  allocate(AVCNHL(3,JD,JV,JH)); AVCNHL=0._r8

  allocate(TWatXChange2WatTableX(JZ,JY,JX));   TWatXChange2WatTableX=0._r8
  allocate(FWatIrrigate2MicP1(JZ,JY,JX));     FWatIrrigate2MicP1=0._r8
  allocate(HeatIrrigation1(JZ,JY,JX));   HeatIrrigation1=0._r8

  allocate(FWatFreezeThawMacP(JZ,JY,JX));   FWatFreezeThawMacP=0._r8

  allocate(HydroCondMacP1(JZ,JY,JX));    HydroCondMacP1=0._r8
  allocate(VLMicP1(0:JZ,JY,JX));  VLMicP1=0._r8
  allocate(VLMacP1(JZ,JY,JX));   VLMacP1=0._r8


  allocate(FLWNX(JY,JX));       FLWNX=0._r8
  allocate(FLWXNX(JY,JX));      FLWXNX=0._r8
  allocate(FLWHNX(JY,JX));      FLWHNX=0._r8
  allocate(HFLWNX(JY,JX));      HFLWNX=0._r8
  allocate(PSISoilMatricPtmp(JZ,JY,JX));   PSISoilMatricPtmp=0._r8

  end subroutine InitWatSubData

!----------------------------------------------------------------------
  subroutine DestructWatSubData
  use abortutils, only : destroy
  implicit none

  call destroy(N6X)

  call destroy(TMLiceThawMicP)

  call destroy(TMLiceThawMacP)

  call destroy(AREAU)
  call destroy(AreaUnderWaterTBL)

  call destroy(VLairMacP)
  call destroy(TTFLXL)

  call destroy(FWatExMacP2MicP)

  call destroy(TWatCharge2MicP)
  call destroy(TConvectWaterFlowMacP)
  call destroy(THeatFlowi)
  call destroy(FWatFreezeThawMicP)
  call destroy(TFLXL)
  call destroy(AVCNHL)

  call destroy(TWatXChange2WatTableX)
  call destroy(FWatIrrigate2MicP1)
  call destroy(HeatIrrigation1)

  call destroy(FWatFreezeThawMacP)

  call destroy(HydroCondMacP1)
  call destroy(VLMicP1)

  call destroy(VLMacP1)

  call destroy(FLWNX)
  call destroy(FLWXNX)
  call destroy(FLWHNX)
  call destroy(HFLWNX)
  call destroy(PSISoilMatricPtmp)

  end subroutine DestructWatSubData

end module WatsubDataMod
