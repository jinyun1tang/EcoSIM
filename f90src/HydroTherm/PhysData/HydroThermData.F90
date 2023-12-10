module HydroThermData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__
  real(r8),allocatable ::  FracSoiPAsWat(:,:,:)                      !
  real(r8),allocatable ::  PSISM1(:,:,:)                      !  
  real(r8),allocatable ::  TKSoi1(:,:,:)                         !  
  real(r8),allocatable ::  DLYRR(:,:)                         !  
  real(r8),allocatable ::  FracSoiPAsIce(:,:,:)               !  
  real(r8),allocatable ::  FracSoiPAsAir(:,:,:)               !  
  real(r8),allocatable ::  VLHeatCapacity(:,:,:)              !whole layer heat capacity (not snow)
  real(r8),allocatable ::  FracSoilAsAirt(:,:,:)              !fraction of soil volume as air, normalized using current pore volume
  real(r8),allocatable ::  VLWatMicP1(:,:,:)                  !    
  real(r8),allocatable ::  VLiceMicP1(:,:,:)                  !  
  real(r8),allocatable ::  SnowFallt(:,:)                         !  
  real(r8),allocatable ::  Rain2Snowt(:,:)                         !  
  real(r8),allocatable ::  VLairMacP1(:,:,:)                  !  positively corrected macropore volume
  real(r8),allocatable ::  EVAPW(:,:)                         !  
  real(r8),allocatable ::  EVAPS(:,:)                         !  
  real(r8),allocatable ::  VapXAir2Sno(:,:)                   ! air-snow exchange of water vapor through evaporation/condensation, sublimation/deposition
  real(r8),allocatable ::  SoilFracAsMacP1(:,:,:)                        !  
  real(r8),allocatable ::  HeatFall2Snowt(:,:)                        !  
  real(r8),allocatable ::  VPQ(:,:)                           !  
  real(r8),allocatable ::  PARSW(:,:)                         !  
  real(r8),allocatable ::  Ice2Snowt(:,:)                         !  
  real(r8),allocatable ::  TKQ(:,:)                           !  
  real(r8),allocatable ::  LWEmscefSnow(:,:)                         !  
  real(r8),allocatable ::  RAGW(:,:)                          !  
  real(r8),allocatable ::  LWRad2Snow(:,:)                         !  
  real(r8),allocatable ::  VLairMicP1(:,:,:)                       ! corrected air-filled micropore volume 
  real(r8),allocatable ::  TKSnow0(:,:,:)                         !  
  real(r8),allocatable ::  VLWatMacP1(:,:,:)                      !
  real(r8),allocatable ::  SoilFracAsMicP(:,:,:)                        !
  real(r8),allocatable ::  VLHeatCapacityA(:,:,:)                      !
  real(r8),allocatable ::  VLHeatCapacityB(:,:,:)                      !
  real(r8),allocatable ::  RAG(:,:)                           !
  real(r8),allocatable ::  PAREW(:,:)                         ! 
  real(r8),allocatable ::  Altitude_grid(:,:)                          ! 
  real(r8),allocatable ::  VLiceMacP1(:,:,:)                      !
  real(r8),allocatable ::  RadSWonSno(:,:)                         !  
  real(r8),allocatable ::  WatFlx2LitRByRunoff(:,:,:,:)                       !  
  real(r8),allocatable ::  HeatFlx2LitRByRunoff(:,:,:,:)                      !  
  real(r8),allocatable ::  DrySnoFlxBySnowRedistribut(:,:,:)                         !  
  real(r8),allocatable ::  WatFlxBySnowRedistribut(:,:,:)                         !
  real(r8),allocatable ::  IceFlxBySnowRedistribut(:,:,:)                         !  
  real(r8),allocatable ::  HeatFlxBySnowRedistribut(:,:,:)                        !  
  real(r8),allocatable ::  HeatFlow2Soili(:,:,:,:)                     !  
  real(r8),allocatable ::  WatXChange2WatTable(:,:,:,:)                      !  
  real(r8),allocatable ::  ConvectWaterFlowMacP(:,:,:,:)                     !  
  real(r8),allocatable ::  WatXChange2WatTableX(:,:,:,:)                     !  
  real(r8),allocatable ::  VLWatMicP2(:,:,:)                       !
  real(r8),allocatable ::  VLairMicP(:,:,:)                      !
  real(r8),allocatable ::  VLWatMicPX1(:,:,:)               !micropore water volume behind wetting front

  public :: InitHydroThermData
  public :: DestructHydroThermData
  contains

!------------------------------------------------------------------------------------------
  subroutine InitHydroThermData
  implicit none


  allocate(FracSoiPAsWat(0:JZ,JY,JX)); FracSoiPAsWat=0._r8
  allocate(PSISM1(0:JZ,JY,JX)); PSISM1=0._r8
  allocate(TKSoi1(0:JZ,JY,JX));    TKSoi1=0._r8    
  allocate(DLYRR(JY,JX));       DLYRR=0._r8  
  allocate(FracSoiPAsIce(0:JZ,JY,JX)); FracSoiPAsIce=0._r8 
  allocate(FracSoiPAsAir(0:JZ,JY,JX)); FracSoiPAsAir=0._r8   
  allocate(VLHeatCapacity(0:JZ,JY,JX));  VLHeatCapacity=0._r8  
  allocate(FracSoilAsAirt(0:JZ,JY,JX)); FracSoilAsAirt=0._r8  
  allocate(VLWatMicP1(0:JZ,JY,JX));  VLWatMicP1=0._r8  
  allocate(VLiceMicP1(0:JZ,JY,JX));  VLiceMicP1=0._r8
  allocate(SnowFallt(JY,JX));       SnowFallt=0._r8
  allocate(Rain2Snowt(JY,JX));       Rain2Snowt=0._r8
  allocate(VLairMacP1(JZ,JY,JX));   VLairMacP1=0._r8 
  allocate(VapXAir2Sno(JY,JX));      VapXAir2Sno=0._r8  
  allocate(EVAPW(JY,JX));       EVAPW=0._r8    
  allocate(EVAPS(JY,JX));       EVAPS=0._r8  
  allocate(SoilFracAsMacP1(JZ,JY,JX));     SoilFracAsMacP1=0._r8  
  allocate(HeatFall2Snowt(JY,JX));      HeatFall2Snowt=0._r8  
  allocate(VPQ(JY,JX));         VPQ=0._r8  
  allocate(PARSW(JY,JX));       PARSW=0._r8
  allocate(Ice2Snowt(JY,JX));       Ice2Snowt=0._r8    
  allocate(TKQ(JY,JX));         TKQ=0._r8  
  allocate(LWEmscefSnow(JY,JX));       LWEmscefSnow=0._r8  
  allocate(RAGW(JY,JX));        RAGW=0._r8  
  allocate(LWRad2Snow(JY,JX));       LWRad2Snow=0._r8  
  allocate(VLairMicP1(0:JZ,JY,JX));  VLairMicP1=0._r8  
  allocate(TKSnow0(JS,JY,JX));      TKSnow0=0._r8  
  allocate(VLWatMacP1(JZ,JY,JX));   VLWatMacP1=0._r8
  allocate(SoilFracAsMicP(JZ,JY,JX));     SoilFracAsMicP=0._r8
  allocate(VLHeatCapacityA(JZ,JY,JX));   VLHeatCapacityA=0._r8
  allocate(VLHeatCapacityB(JZ,JY,JX));   VLHeatCapacityB=0._r8  
  allocate(RAG(JY,JX));         RAG=0._r8
  allocate(PAREW(JY,JX));       PAREW=0._r8
  allocate(Altitude_grid(JY,JX));        Altitude_grid=0._r8  
  allocate(VLiceMacP1(JZ,JY,JX));   VLiceMacP1=0._r8  
  allocate(RadSWonSno(JY,JX));       RadSWonSno=0._r8  
  allocate(WatFlx2LitRByRunoff(2,2,JV,JH));     WatFlx2LitRByRunoff=0._r8  
  allocate(HeatFlx2LitRByRunoff(2,2,JV,JH));    HeatFlx2LitRByRunoff=0._r8
  allocate(DrySnoFlxBySnowRedistribut(2,JV,JH));       DrySnoFlxBySnowRedistribut=0._r8
  allocate(WatFlxBySnowRedistribut(2,JV,JH));       WatFlxBySnowRedistribut=0._r8
  allocate(IceFlxBySnowRedistribut(2,JV,JH));       IceFlxBySnowRedistribut=0._r8  
  allocate(HeatFlxBySnowRedistribut(2,JV,JH));      HeatFlxBySnowRedistribut=0._r8
  allocate(HeatFlow2Soili(3,JD,JV,JH));  HeatFlow2Soili=0._r8
  allocate(WatXChange2WatTable(3,JD,JV,JH));   WatXChange2WatTable=0._r8  
  allocate(ConvectWaterFlowMacP(3,JD,JV,JH));  ConvectWaterFlowMacP=0._r8
  allocate(WatXChange2WatTableX(3,JD,JV,JH));  WatXChange2WatTableX=0._r8
  allocate(VLWatMicP2(JZ,JY,JX));    VLWatMicP2=0._r8
  allocate(VLairMicP(JZ,JY,JX));   VLairMicP=0._r8
  allocate(VLWatMicPX1(JZ,JY,JX));   VLWatMicPX1=0._r8

  end subroutine InitHydroThermData

!------------------------------------------------------------------------------------------
  subroutine DestructHydroThermData
  use abortutils, only : destroy  
  implicit none

  call destroy(WatFlx2LitRByRunoff)
  call destroy(FracSoiPAsWat)
  call destroy(PSISM1)  
  call destroy(TKSoi1)  
  call destroy(DLYRR)  
  call destroy(FracSoiPAsIce)  
  call destroy(FracSoiPAsAir)  
  call destroy(VLHeatCapacity)  
  call destroy(FracSoilAsAirt)  
  call destroy(VLWatMicP1)  
  call destroy(VLiceMicP1)
  call destroy(SnowFallt)
  call destroy(Rain2Snowt)
  call destroy(VLairMacP1)    
  call destroy(EVAPW)  
  call destroy(EVAPS)  
  call destroy(VapXAir2Sno)  
  call destroy(SoilFracAsMacP1)  
  call destroy(HeatFall2Snowt)
  call destroy(VPQ)    
  call destroy(PARSW)  
  call destroy(Ice2Snowt)
  call destroy(TKQ)
  call destroy(LWEmscefSnow)
  call destroy(RAGW) 
  call destroy(LWRad2Snow)
  call destroy(VLairMicP1)
  call destroy(TKSnow0)  
  call destroy(VLWatMacP1)
  call destroy(SoilFracAsMicP)
  call destroy(VLHeatCapacityA)
  call destroy(VLHeatCapacityB)  
  call destroy(RAG)
  call destroy(PAREW)
  call destroy(Altitude_grid)
  call destroy(VLiceMacP1)
  call destroy(RadSWonSno)
  call destroy(HeatFlx2LitRByRunoff)    
  call destroy(DrySnoFlxBySnowRedistribut)  
  call destroy(WatFlxBySnowRedistribut)
  call destroy(IceFlxBySnowRedistribut)
  call destroy(HeatFlxBySnowRedistribut)    
  call destroy(HeatFlow2Soili)  
  call destroy(WatXChange2WatTable)  
  call destroy(ConvectWaterFlowMacP)
  call destroy(WatXChange2WatTableX)
  call destroy(VLWatMicP2)
  call destroy(VLairMicP)
  call destroy(VLWatMicPX1)

  end subroutine DestructHydroThermData
end module HydroThermData