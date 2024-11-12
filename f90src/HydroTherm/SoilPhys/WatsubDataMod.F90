module WatsubDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),allocatable ::  TMLiceThawMicP(:,:,:)              !micropore layer integrated ice mass loss due to thaw
  real(r8),allocatable ::  TMLiceThawMacP(:,:,:)              !macropore layer integrated ice mass loss due to thaw

  real(r8),allocatable ::  AREAU(:,:,:)                       !
  real(r8),allocatable ::  AreaUnderWaterTBL(:,:,:)           !

  real(r8),allocatable ::  VLairMacP_vr(:,:,:)                   !
  real(r8),allocatable ::  TLPhaseChangeHeat2Soi1(:,:,:)        !total soil layer latent heat release from melting
  real(r8),allocatable ::  TLPhaseChangeHeat2Soi1s(:,:,:)        !total soil layer latent heat release from melting

  real(r8),allocatable ::  FWatExMacP2MicPi(:,:,:)             !pressure-driven water flow from macpore to micpore

  real(r8),allocatable ::  TWatCharge2MicP(:,:,:)                       !
  real(r8),allocatable ::  TConvWaterFlowMacP_3D_vr(:,:,:)                      !
  real(r8),allocatable ::  THeatFlow2Soili_vr(:,:,:)                      !
  real(r8),allocatable ::  FIceThawMicP(:,:,:)                       !
  real(r8),allocatable ::  SoiPLIceHeatFlxFrez(:,:,:)                       !
  real(r8),allocatable ::  AVCNHL(:,:,:,:)                    !

  real(r8),allocatable ::  TWatXChange2WatTableX(:,:,:)                      !
  real(r8),allocatable ::  FWatIrrigate2MicP1_vr(:,:,:)                        !
  real(r8),allocatable ::  HeatIrrigation1(:,:,:)                      !

  real(r8),allocatable ::  FIceThawMacP(:,:,:)                      !

  real(r8),allocatable ::  HydroCondMacP1_vr(:,:,:)                       !
  real(r8),allocatable ::  VLMicP1_vr(:,:,:)                       !
  real(r8),allocatable ::  VLMacP1_vr(:,:,:)                      !

  real(r8),allocatable ::  H2OFlow2TopSoiMicP_col(:,:)                         !
  real(r8),allocatable ::  H2OFlow2TopSoiMicPX_col(:,:)                        !
  real(r8),allocatable ::  H2OFlow2TopSoiMacP_col(:,:)                        !
  real(r8),allocatable ::  HeatFlow2TopSoi_col(:,:)                        !
  real(r8),allocatable ::  PSISoilMatricPtmp_vr(:,:,:)                      !


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


  allocate(VLairMacP_vr(JZ,JY,JX));  VLairMacP_vr=0._r8
  allocate(TLPhaseChangeHeat2Soi1(JZ,JY,JX));   TLPhaseChangeHeat2Soi1=0._r8
  allocate(TLPhaseChangeHeat2Soi1s(JZ,JY,JX));   TLPhaseChangeHeat2Soi1s=0._r8

  allocate(FWatExMacP2MicPi(JZ,JY,JX));    FWatExMacP2MicPi=0._r8

  allocate(TWatCharge2MicP(JZ,JY,JX));    TWatCharge2MicP=0._r8
  allocate(TConvWaterFlowMacP_3D_vr(JZ,JY,JX));   TConvWaterFlowMacP_3D_vr=0._r8
  allocate(THeatFlow2Soili_vr(JZ,JY,JX));   THeatFlow2Soili_vr=0._r8
  allocate(FIceThawMicP(JZ,JY,JX));    FIceThawMicP=0._r8
  allocate(SoiPLIceHeatFlxFrez(JZ,JY,JX));    SoiPLIceHeatFlxFrez=0._r8
  allocate(AVCNHL(3,JD,JV,JH)); AVCNHL=0._r8

  allocate(TWatXChange2WatTableX(JZ,JY,JX));   TWatXChange2WatTableX=0._r8
  allocate(FWatIrrigate2MicP1_vr(JZ,JY,JX));     FWatIrrigate2MicP1_vr=0._r8
  allocate(HeatIrrigation1(JZ,JY,JX));   HeatIrrigation1=0._r8

  allocate(FIceThawMacP(JZ,JY,JX));   FIceThawMacP=0._r8

  allocate(HydroCondMacP1_vr(JZ,JY,JX));    HydroCondMacP1_vr=0._r8
  allocate(VLMicP1_vr(0:JZ,JY,JX));  VLMicP1_vr=0._r8
  allocate(VLMacP1_vr(JZ,JY,JX));   VLMacP1_vr=0._r8


  allocate(H2OFlow2TopSoiMicP_col(JY,JX));       H2OFlow2TopSoiMicP_col=0._r8
  allocate(H2OFlow2TopSoiMicPX_col(JY,JX));      H2OFlow2TopSoiMicPX_col=0._r8
  allocate(H2OFlow2TopSoiMacP_col(JY,JX));      H2OFlow2TopSoiMacP_col=0._r8
  allocate(HeatFlow2TopSoi_col(JY,JX));      HeatFlow2TopSoi_col=0._r8
  allocate(PSISoilMatricPtmp_vr(JZ,JY,JX));   PSISoilMatricPtmp_vr=0._r8

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

  call destroy(VLairMacP_vr)
  call destroy(TLPhaseChangeHeat2Soi1)
  call destroy(TLPhaseChangeHeat2Soi1s)

  call destroy(FWatExMacP2MicPi)

  call destroy(TWatCharge2MicP)
  call destroy(TConvWaterFlowMacP_3D_vr)
  call destroy(THeatFlow2Soili_vr)
  call destroy(FIceThawMicP)
  call destroy(SoiPLIceHeatFlxFrez)
  call destroy(AVCNHL)

  call destroy(TWatXChange2WatTableX)
  call destroy(FWatIrrigate2MicP1_vr)
  call destroy(HeatIrrigation1)

  call destroy(FIceThawMacP)

  call destroy(HydroCondMacP1_vr)
  call destroy(VLMicP1_vr)

  call destroy(VLMacP1_vr)

  call destroy(H2OFlow2TopSoiMicP_col)
  call destroy(H2OFlow2TopSoiMicPX_col)
  call destroy(H2OFlow2TopSoiMacP_col)
  call destroy(HeatFlow2TopSoi_col)
  call destroy(PSISoilMatricPtmp_vr)

  end subroutine DestructWatSubData

end module WatsubDataMod
