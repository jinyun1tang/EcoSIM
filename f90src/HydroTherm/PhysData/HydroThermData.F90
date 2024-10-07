module HydroThermData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),allocatable ::  FracSoiPAsWat_vr(:,:,:)                      !
  real(r8),allocatable ::  PSISM1_vr(:,:,:)                      !  
  real(r8),allocatable ::  TKSoi1_vr(:,:,:)                         !  
  real(r8),allocatable ::  DLYRR_COL(:,:)                         !  
  real(r8),allocatable ::  FracSoiPAsIce_vr(:,:,:)               !  
  real(r8),allocatable ::  FracSoiPAsAir_vr(:,:,:)               !  
  real(r8),allocatable ::  VHeatCapacity1_vr(:,:,:)              !whole layer heat capacity (not snow)
  real(r8),allocatable ::  FracSoilAsAirt(:,:,:)              !fraction of soil volume as air, normalized using current pore volume
  real(r8),allocatable ::  VLWatMicP1_vr(:,:,:)                  !    
  real(r8),allocatable ::  VLiceMicP1_vr(:,:,:)                  !  
  real(r8),allocatable ::  SnowFallt(:,:)                         !  
  real(r8),allocatable ::  Rain2Snowt(:,:)                         !  
  real(r8),allocatable ::  VLairMacP1_vr(:,:,:)                  !  positively corrected macropore volume
  real(r8),allocatable ::  EVAPW(:,:)                         !  
  real(r8),allocatable ::  EVAPS(:,:)                         !  
  real(r8),allocatable ::  VapXAir2Sno(:,:)                   ! air-snow exchange of water vapor through evaporation/condensation, sublimation/deposition
  real(r8),allocatable ::  SoilFracAsMacP1_vr(:,:,:)                        !  
  real(r8),allocatable ::  HeatFall2Snowt(:,:)                        !  
  real(r8),allocatable ::  VPQ_col(:,:)                           !  
  real(r8),allocatable ::  PARSW(:,:)                         !  
  real(r8),allocatable ::  Ice2Snowt(:,:)                         !  
  real(r8),allocatable ::  TKQ(:,:)                           !  
  real(r8),allocatable ::  LWEmscefSnow_col(:,:)                         !  
  real(r8),allocatable ::  RAGW(:,:)                          !  
  real(r8),allocatable ::  LWRad2Snow(:,:)                         !  
  real(r8),allocatable ::  VLairMicP1_vr(:,:,:)                       ! corrected air-filled micropore volume 
  real(r8),allocatable ::  TKSnow0_snvr(:,:,:)                         !  
  real(r8),allocatable ::  VLWatMacP1_vr(:,:,:)                      !
  real(r8),allocatable ::  SoilFracAsMicP_vr(:,:,:)                        !
  real(r8),allocatable ::  VLHeatCapacityA_vr(:,:,:)                      !
  real(r8),allocatable ::  VLHeatCapacityB_vr(:,:,:)                      !
  real(r8),allocatable ::  RAG(:,:)                           !
  real(r8),allocatable ::  PAREW(:,:)                         ! 
  real(r8),allocatable ::  Altitude_grid(:,:)                          ! 
  real(r8),allocatable ::  VLiceMacP1_vr(:,:,:)                      !
  real(r8),allocatable ::  RadSWonSno(:,:)                         !  
  real(r8),allocatable ::  WatFlx2LitRByRunoff(:,:,:,:)                       !  
  real(r8),allocatable ::  HeatFlx2LitRByRunoff(:,:,:,:)                      !  
  real(r8),allocatable ::  DrySnoFlxBySnowRedistribut(:,:,:)                         !  
  real(r8),allocatable ::  WatFlxBySnowRedistribut(:,:,:)                         !
  real(r8),allocatable ::  IceFlxBySnowRedistribut(:,:,:)                         !  
  real(r8),allocatable ::  HeatFlxBySnowRedistribut(:,:,:)                        !  
  real(r8),allocatable ::  HeatFlow2Soili(:,:,:,:)                     !  
  real(r8),allocatable ::  WatXChange2WatTable(:,:,:,:)                      !  
  real(r8),allocatable ::  ConvWaterFlowMacP_3D(:,:,:,:)                     !  
  real(r8),allocatable ::  WatXChange2WatTableX(:,:,:,:)                     !  
  real(r8),allocatable ::  VLWatMicP2_vr(:,:,:)                       !
  real(r8),allocatable ::  VLairMicP_vr(:,:,:)                      !
  real(r8),allocatable ::  VLWatMicPX1_vr(:,:,:)               !micropore water volume behind wetting front

  public :: InitHydroThermData
  public :: DestructHydroThermData
  contains

!------------------------------------------------------------------------------------------
  subroutine InitHydroThermData
  implicit none


  allocate(FracSoiPAsWat_vr(0:JZ,JY,JX)); FracSoiPAsWat_vr=0._r8
  allocate(PSISM1_vr(0:JZ,JY,JX)); PSISM1_vr=0._r8
  allocate(TKSoi1_vr(0:JZ,JY,JX));    TKSoi1_vr=0._r8    
  allocate(DLYRR_COL(JY,JX));       DLYRR_col=0._r8  
  allocate(FracSoiPAsIce_vr(0:JZ,JY,JX)); FracSoiPAsIce_vr=0._r8 
  allocate(FracSoiPAsAir_vr(0:JZ,JY,JX)); FracSoiPAsAir_vr=0._r8   
  allocate(VHeatCapacity1_vr(0:JZ,JY,JX));  VHeatCapacity1_vr=0._r8  
  allocate(FracSoilAsAirt(0:JZ,JY,JX)); FracSoilAsAirt=0._r8  
  allocate(VLWatMicP1_vr(0:JZ,JY,JX));  VLWatMicP1_vr=0._r8  
  allocate(VLiceMicP1_vr(0:JZ,JY,JX));  VLiceMicP1_vr=0._r8
  allocate(SnowFallt(JY,JX));       SnowFallt=0._r8
  allocate(Rain2Snowt(JY,JX));       Rain2Snowt=0._r8
  allocate(VLairMacP1_vr(JZ,JY,JX));   VLairMacP1_vr=0._r8 
  allocate(VapXAir2Sno(JY,JX));      VapXAir2Sno=0._r8  
  allocate(EVAPW(JY,JX));       EVAPW=0._r8    
  allocate(EVAPS(JY,JX));       EVAPS=0._r8  
  allocate(SoilFracAsMacP1_vr(JZ,JY,JX));     SoilFracAsMacP1_vr=0._r8  
  allocate(HeatFall2Snowt(JY,JX));      HeatFall2Snowt=0._r8  
  allocate(VPQ_col(JY,JX));         VPQ_col=0._r8  
  allocate(PARSW(JY,JX));       PARSW=0._r8
  allocate(Ice2Snowt(JY,JX));       Ice2Snowt=0._r8    
  allocate(TKQ(JY,JX));         TKQ=0._r8  
  allocate(LWEmscefSnow_col(JY,JX));       LWEmscefSnow_col=0._r8  
  allocate(RAGW(JY,JX));        RAGW=0._r8  
  allocate(LWRad2Snow(JY,JX));       LWRad2Snow=0._r8  
  allocate(VLairMicP1_vr(0:JZ,JY,JX));  VLairMicP1_vr=0._r8  
  allocate(TKSnow0_snvr(JS,JY,JX));      TKSnow0_snvr=0._r8  
  allocate(VLWatMacP1_vr(JZ,JY,JX));   VLWatMacP1_vr=0._r8
  allocate(SoilFracAsMicP_vr(JZ,JY,JX));     SoilFracAsMicP_vr=0._r8
  allocate(VLHeatCapacityA_vr(JZ,JY,JX));   VLHeatCapacityA_vr=0._r8
  allocate(VLHeatCapacityB_vr(JZ,JY,JX));   VLHeatCapacityB_vr=0._r8  
  allocate(RAG(JY,JX));         RAG=0._r8
  allocate(PAREW(JY,JX));       PAREW=0._r8
  allocate(Altitude_grid(JY,JX));        Altitude_grid=0._r8  
  allocate(VLiceMacP1_vr(JZ,JY,JX));   VLiceMacP1_vr=0._r8  
  allocate(RadSWonSno(JY,JX));       RadSWonSno=0._r8  
  allocate(WatFlx2LitRByRunoff(2,2,JV,JH));     WatFlx2LitRByRunoff=0._r8  
  allocate(HeatFlx2LitRByRunoff(2,2,JV,JH));    HeatFlx2LitRByRunoff=0._r8
  allocate(DrySnoFlxBySnowRedistribut(2,JV,JH));       DrySnoFlxBySnowRedistribut=0._r8
  allocate(WatFlxBySnowRedistribut(2,JV,JH));       WatFlxBySnowRedistribut=0._r8
  allocate(IceFlxBySnowRedistribut(2,JV,JH));       IceFlxBySnowRedistribut=0._r8  
  allocate(HeatFlxBySnowRedistribut(2,JV,JH));      HeatFlxBySnowRedistribut=0._r8
  allocate(HeatFlow2Soili(3,JD,JV,JH));  HeatFlow2Soili=0._r8
  allocate(WatXChange2WatTable(3,JD,JV,JH));   WatXChange2WatTable=0._r8  
  allocate(ConvWaterFlowMacP_3D(3,JD,JV,JH));  ConvWaterFlowMacP_3D=0._r8
  allocate(WatXChange2WatTableX(3,JD,JV,JH));  WatXChange2WatTableX=0._r8
  allocate(VLWatMicP2_vr(JZ,JY,JX));    VLWatMicP2_vr=0._r8
  allocate(VLairMicP_vr(JZ,JY,JX));   VLairMicP_vr=0._r8
  allocate(VLWatMicPX1_vr(JZ,JY,JX));   VLWatMicPX1_vr=0._r8

  end subroutine InitHydroThermData

!------------------------------------------------------------------------------------------
  subroutine DestructHydroThermData
  use abortutils, only : destroy  
  implicit none

  call destroy(WatFlx2LitRByRunoff)
  call destroy(FracSoiPAsWat_vr)
  call destroy(PSISM1_vr)  
  call destroy(TKSoi1_vr)  
  call destroy(DLYRR_COL)  
  call destroy(FracSoiPAsIce_vr)  
  call destroy(FracSoiPAsAir_vr)  
  call destroy(VHeatCapacity1_vr)  
  call destroy(FracSoilAsAirt)  
  call destroy(VLWatMicP1_vr)  
  call destroy(VLiceMicP1_vr)
  call destroy(SnowFallt)
  call destroy(Rain2Snowt)
  call destroy(VLairMacP1_vr)    
  call destroy(EVAPW)  
  call destroy(EVAPS)  
  call destroy(VapXAir2Sno)  
  call destroy(SoilFracAsMacP1_vr)  
  call destroy(HeatFall2Snowt)
  call destroy(VPQ_col)    
  call destroy(PARSW)  
  call destroy(Ice2Snowt)
  call destroy(TKQ)
  call destroy(LWEmscefSnow_col)
  call destroy(RAGW) 
  call destroy(LWRad2Snow)
  call destroy(VLairMicP1_vr)
  call destroy(TKSnow0_snvr)  
  call destroy(VLWatMacP1_vr)
  call destroy(SoilFracAsMicP_vr)
  call destroy(VLHeatCapacityA_vr)
  call destroy(VLHeatCapacityB_vr)  
  call destroy(RAG)
  call destroy(PAREW)
  call destroy(Altitude_grid)
  call destroy(VLiceMacP1_vr)
  call destroy(RadSWonSno)
  call destroy(HeatFlx2LitRByRunoff)    
  call destroy(DrySnoFlxBySnowRedistribut)  
  call destroy(WatFlxBySnowRedistribut)
  call destroy(IceFlxBySnowRedistribut)
  call destroy(HeatFlxBySnowRedistribut)    
  call destroy(HeatFlow2Soili)  
  call destroy(WatXChange2WatTable)  
  call destroy(ConvWaterFlowMacP_3D)
  call destroy(WatXChange2WatTableX)
  call destroy(VLWatMicP2_vr)
  call destroy(VLairMicP_vr)
  call destroy(VLWatMicPX1_vr)

  end subroutine DestructHydroThermData
end module HydroThermData