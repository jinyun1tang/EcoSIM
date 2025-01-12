module WatsubDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),allocatable ::  TMLiceThawMicP_vr(:,:,:)              !micropore layer integrated ice mass loss due to thaw, 
  real(r8),allocatable ::  TMLiceThawMacP_vr(:,:,:)              !macropore layer integrated ice mass loss due to thaw, 

  real(r8),allocatable ::  AREAU(:,:,:)                       !
  real(r8),allocatable ::  AreaUnderWaterTBL(:,:,:)           !
  real(r8),allocatable ::  VLairMacP_vr(:,:,:)                   !
  real(r8),allocatable ::  TLPhaseChangeHeat2Soi1_vr(:,:,:)         !total soil layer latent heat release from melting

  real(r8),allocatable ::  FWatExMacP2MicPi_vr(:,:,:)             !pressure-driven water flow from macpore to micpore

  real(r8),allocatable ::  TWatFlow2MicP_3DM_vr(:,:,:)           !total water charge to micropore in iteration M
  real(r8),allocatable ::  TWaterFlow2Macpt_3DM_vr(:,:,:)          !total water charge to macropore in iteration M
  real(r8),allocatable ::  THeatFlow2Soil_3DM_vr(:,:,:)            !total heat flux add to soil layer in iteration M
  real(r8),allocatable ::  FIceThawMicP(:,:,:)                       !
  real(r8),allocatable ::  SoiPLIceHeatFlxFrez(:,:,:)                       !
  real(r8),allocatable ::  AVCNHL(:,:,:,:)                    !

  real(r8),allocatable ::  TWaterFlow2MicptX_3DM_vr(:,:,:)                      !
  real(r8),allocatable ::  FWatIrrigate2MicP1_vr(:,:,:)                        !
  real(r8),allocatable ::  HeatIrrigation1_vr(:,:,:)                  !heat flux due to irrigation, [MJ/d2]

  real(r8),allocatable ::  FIceThawMacP(:,:,:)                      !

  real(r8),allocatable ::  HydroCondMacP1_vr(:,:,:)                       !
  real(r8),allocatable ::  VLMicP1_vr(:,:,:)                       !
  real(r8),allocatable ::  VLMacP1_vr(:,:,:)                      !

  real(r8),allocatable ::  H2OFlow2TopSoiMicP_col(:,:)                         !
  real(r8),allocatable ::  H2OFlow2TopSoiMicPX_col(:,:)                        !
  real(r8),allocatable ::  H2OFlow2TopSoiMacP_col(:,:)                        !
  real(r8),allocatable ::  HeatFlow2TopSoi_col(:,:)                        !
  real(r8),allocatable ::  PSISoilMatricPtmp_vr(:,:,:)                      !
  real(r8),allocatable :: QDischarM_col(:,:)
  real(r8),allocatable :: QDrainM_col(:,:)

  integer, allocatable ::  N6X(:,:)

!----------------------------------------------------------------------
  public :: InitWatsubData
contains
  subroutine InitWatSubData

  implicit none

  allocate(N6X(JY,JX));         N6X=0
  allocate(TMLiceThawMicP_vr(JZ,JY,JX));   TMLiceThawMicP_vr=0._r8
  allocate(TMLiceThawMacP_vr(JZ,JY,JX));   TMLiceThawMacP_vr=0._r8
  allocate(AREAU(JZ,JY,JX));    AREAU=0._r8
  allocate(AreaUnderWaterTBL(JZ,JY,JX));   AreaUnderWaterTBL=0._r8
  allocate(QDischarM_col(JY,JX)); QDischarM_col(:,:)=0._r8
  allocate(QDrainM_col(JY,JX)); QDrainM_col(:,:)=0._r8
  allocate(VLairMacP_vr(JZ,JY,JX));  VLairMacP_vr=0._r8
  allocate(TLPhaseChangeHeat2Soi1_vr(JZ,JY,JX));   TLPhaseChangeHeat2Soi1_vr=0._r8
  allocate(FWatExMacP2MicPi_vr(JZ,JY,JX));    FWatExMacP2MicPi_vr=0._r8

  allocate(TWatFlow2MicP_3DM_vr(JZ,JY,JX));    TWatFlow2MicP_3DM_vr=0._r8
  allocate(TWaterFlow2Macpt_3DM_vr(JZ,JY,JX));   TWaterFlow2Macpt_3DM_vr=0._r8
  allocate(THeatFlow2Soil_3DM_vr(JZ,JY,JX));   THeatFlow2Soil_3DM_vr=0._r8
  allocate(FIceThawMicP(JZ,JY,JX));    FIceThawMicP=0._r8
  allocate(SoiPLIceHeatFlxFrez(JZ,JY,JX));    SoiPLIceHeatFlxFrez=0._r8
  allocate(AVCNHL(3,JD,JV,JH)); AVCNHL=0._r8

  allocate(TWaterFlow2MicptX_3DM_vr(JZ,JY,JX));   TWaterFlow2MicptX_3DM_vr=0._r8
  allocate(FWatIrrigate2MicP1_vr(JZ,JY,JX));     FWatIrrigate2MicP1_vr=0._r8
  allocate(HeatIrrigation1_vr(JZ,JY,JX));   HeatIrrigation1_vr=0._r8

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
  call destroy(TMLiceThawMicP_vr)
  call destroy(TMLiceThawMacP_vr)
  call destroy(AREAU)
  call destroy(AreaUnderWaterTBL)
  call destroy(VLairMacP_vr)
  call destroy(TLPhaseChangeHeat2Soi1_vr)
  call destroy(FWatExMacP2MicPi_vr)
  call destroy(TWatFlow2MicP_3DM_vr)
  call destroy(TWaterFlow2Macpt_3DM_vr)
  call destroy(THeatFlow2Soil_3DM_vr)
  call destroy(FIceThawMicP)
  call destroy(SoiPLIceHeatFlxFrez)
  call destroy(AVCNHL)
  call destroy(TWaterFlow2MicptX_3DM_vr)
  call destroy(FWatIrrigate2MicP1_vr)
  call destroy(HeatIrrigation1_vr)
  call destroy(FIceThawMacP)
  call destroy(HydroCondMacP1_vr)
  call destroy(VLMicP1_vr)
  call destroy(VLMacP1_vr)
  call destroy(H2OFlow2TopSoiMicP_col)
  call destroy(H2OFlow2TopSoiMicPX_col)
  call destroy(H2OFlow2TopSoiMacP_col)
  call destroy(HeatFlow2TopSoi_col)
  call destroy(PSISoilMatricPtmp_vr)
  call destroy(QDischarM_col)
  call destroy(QDrainM_col)

  end subroutine DestructWatSubData

end module WatsubDataMod
