module WatsubDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  TWFLXL(:,:,:)                      !

  real(r8),allocatable ::  TWFLXH(:,:,:)                      !

  real(r8),allocatable ::  AREAU(:,:,:)                       !
  real(r8),allocatable ::  AreaUnderWaterTable(:,:,:)                      !

  real(r8),allocatable ::  VOLPH1Z(:,:,:)                     !
  real(r8),allocatable ::  TTFLXL(:,:,:)                      !

  real(r8),allocatable ::  FINHL(:,:,:)                       !

  real(r8),allocatable ::  TFLWL(:,:,:)                       !
  real(r8),allocatable ::  TConvectWaterFlowMacP(:,:,:)                      !
  real(r8),allocatable ::  THeatFlowi(:,:,:)                      !
  real(r8),allocatable ::  WFLXL(:,:,:)                       !
  real(r8),allocatable ::  TFLXL(:,:,:)                       !
  real(r8),allocatable ::  AVCNHL(:,:,:,:)                    !

  real(r8),allocatable ::  TWatXChange2WatTableX(:,:,:)                      !
  real(r8),allocatable ::  FLU1(:,:,:)                        !
  real(r8),allocatable ::  HeatIrrigation1(:,:,:)                      !

  real(r8),allocatable ::  WFLXLH(:,:,:)                      !

  real(r8),allocatable ::  CNDH1(:,:,:)                       !
  real(r8),allocatable ::  VMicP1(:,:,:)                       !
  real(r8),allocatable ::  VAirMacP1(:,:,:)                      !

  real(r8),allocatable ::  FLWNX(:,:)                         !
  real(r8),allocatable ::  FLWXNX(:,:)                        !
  real(r8),allocatable ::  FLWHNX(:,:)                        !
  real(r8),allocatable ::  HFLWNX(:,:)                        !
  real(r8),allocatable ::  PSISA1(:,:,:)                      !


  integer, allocatable ::  N6X(:,:)

!----------------------------------------------------------------------
  public :: InitWatsubData
contains
  subroutine InitWatSubData

  implicit none

  allocate(N6X(JY,JX));         N6X=0

  allocate(TWFLXL(JZ,JY,JX));   TWFLXL=0._r8

  allocate(TWFLXH(JZ,JY,JX));   TWFLXH=0._r8

  allocate(AREAU(JZ,JY,JX));    AREAU=0._r8
  allocate(AreaUnderWaterTable(JZ,JY,JX));   AreaUnderWaterTable=0._r8


  allocate(VOLPH1Z(JZ,JY,JX));  VOLPH1Z=0._r8
  allocate(TTFLXL(JZ,JY,JX));   TTFLXL=0._r8

  allocate(FINHL(JZ,JY,JX));    FINHL=0._r8

  allocate(TFLWL(JZ,JY,JX));    TFLWL=0._r8
  allocate(TConvectWaterFlowMacP(JZ,JY,JX));   TConvectWaterFlowMacP=0._r8
  allocate(THeatFlowi(JZ,JY,JX));   THeatFlowi=0._r8
  allocate(WFLXL(JZ,JY,JX));    WFLXL=0._r8
  allocate(TFLXL(JZ,JY,JX));    TFLXL=0._r8
  allocate(AVCNHL(3,JD,JV,JH)); AVCNHL=0._r8

  allocate(TWatXChange2WatTableX(JZ,JY,JX));   TWatXChange2WatTableX=0._r8
  allocate(FLU1(JZ,JY,JX));     FLU1=0._r8
  allocate(HeatIrrigation1(JZ,JY,JX));   HeatIrrigation1=0._r8

  allocate(WFLXLH(JZ,JY,JX));   WFLXLH=0._r8

  allocate(CNDH1(JZ,JY,JX));    CNDH1=0._r8
  allocate(VMicP1(0:JZ,JY,JX));  VMicP1=0._r8

  allocate(VAirMacP1(JZ,JY,JX));   VAirMacP1=0._r8


  allocate(FLWNX(JY,JX));       FLWNX=0._r8
  allocate(FLWXNX(JY,JX));      FLWXNX=0._r8
  allocate(FLWHNX(JY,JX));      FLWHNX=0._r8
  allocate(HFLWNX(JY,JX));      HFLWNX=0._r8
  allocate(PSISA1(JZ,JY,JX));   PSISA1=0._r8

  end subroutine InitWatSubData

!----------------------------------------------------------------------
  subroutine DestructWatSubData
  use abortutils, only : destroy
  implicit none

  call destroy(N6X)

  call destroy(TWFLXL)

  call destroy(TWFLXH)

  call destroy(AREAU)
  call destroy(AreaUnderWaterTable)

  call destroy(VOLPH1Z)
  call destroy(TTFLXL)

  call destroy(FINHL)

  call destroy(TFLWL)
  call destroy(TConvectWaterFlowMacP)
  call destroy(THeatFlowi)
  call destroy(WFLXL)
  call destroy(TFLXL)
  call destroy(AVCNHL)

  call destroy(TWatXChange2WatTableX)
  call destroy(FLU1)
  call destroy(HeatIrrigation1)

  call destroy(WFLXLH)

  call destroy(CNDH1)
  call destroy(VMicP1)

  call destroy(VAirMacP1)

  call destroy(FLWNX)
  call destroy(FLWXNX)
  call destroy(FLWHNX)
  call destroy(HFLWNX)
  call destroy(PSISA1)

  end subroutine DestructWatSubData

end module WatsubDataMod
