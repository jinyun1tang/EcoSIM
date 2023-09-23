module SnowPhysData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod

  implicit none

  character(len=*), private, parameter :: mod_filename=__FILE__

  real(r8),allocatable ::  trcg_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcn_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcsa_TBLS(:,:,:,:)                      !
  real(r8),allocatable ::  TQS(:,:)                           !
  real(r8),allocatable ::  trcsa_TQS(:,:,:)  
  real(r8),allocatable ::  trcn_TQR(:,:,:)                        !
  real(r8),allocatable ::  trcsa_TQR(:,:,:)                         !
  real(r8),allocatable ::  trcg_QSS(:,:,:)
  real(r8),allocatable ::  trcn_QSS(:,:,:)
  real(r8),allocatable ::  trcg_TQR(:,:,:)                        !
  real(r8),allocatable ::  XSnowThawMassL(:,:,:)              !hourly convective heat flux from snow transfer
  real(r8),allocatable ::  XIceThawMassL(:,:,:)                      !hourly convective heat flux from ice transfer  
  real(r8),allocatable ::  THQS(:,:)                          !
  real(r8),allocatable ::  TQI(:,:)                           !
  real(r8),allocatable ::  TQW(:,:)                           !
  real(r8),allocatable ::  TQR1(:,:)                          !
  real(r8),allocatable ::  THQR1(:,:)                         !
  real(r8),allocatable ::  TQS1(:,:)                          !
  real(r8),allocatable ::  TQW1(:,:)                          !
  real(r8),allocatable ::  TQI1(:,:)                          !
  real(r8),allocatable ::  THQS1(:,:)                         !
  real(r8),allocatable ::  TKSnow1(:,:,:)                     ! temporary snow layer temperature in Kelvin
  real(r8),allocatable ::  VLHeatCapSnowM1(:,:,:)             ! temporary snow layer volumetric heat capacity 
  real(r8),allocatable ::  VLIceSnow0M(:,:,:)                 !
  real(r8),allocatable ::  VLWatSnow0M(:,:,:)                 !snow-held water volume during iteration
  real(r8),allocatable ::  CumSno2SnowLay(:,:,:)                       !
  real(r8),allocatable ::  CumWat2SnowLay(:,:,:)                       !
  real(r8),allocatable ::  CumIce2SnowLay(:,:,:)                       !
  real(r8),allocatable ::  CumHeat2SnowLay(:,:,:)                      !
  real(r8),allocatable ::  PhaseChangeHeatL(:,:,:)                       !
  real(r8),allocatable ::  SnowThawMassL(:,:,:)                       !
  real(r8),allocatable ::  IceThawMassL(:,:,:)                       !
  real(r8),allocatable ::  WatX2SnoLay(:,:,:)                       !
  real(r8),allocatable ::  SnoX2SnoLay(:,:,:)                       !
  real(r8),allocatable ::  IceX2SnoLay(:,:,:)                       !
  real(r8),allocatable ::  HeatX2SnoLay(:,:,:)                      !
  real(r8),allocatable ::  VLDrySnoWE0M(:,:,:)                !dry snow layer volume during iteration
  real(r8),allocatable ::  VLDrySnoWE0(:,:,:)                 !dry snow layer volume before update
  real(r8),allocatable ::  VLIceSnow0(:,:,:)                  !snow-held ice layer volume before update
  real(r8),allocatable ::  VLWatSnow0(:,:,:)                       !snow-held water layer volume before update
  real(r8),allocatable ::  VLSnoDWI1(:,:,:)                   !
  real(r8),allocatable ::  SnowLayerThick0(:,:,:)                      !
   
  public :: InitSnowPhysData
  public :: DestructSnowPhysData
  contains
!----------------------------------------------------------------------  
  subroutine  InitSnowPhysData
  implicit none

  allocate(trcg_TBLS(idg_beg:idg_end-1,JS,JY,JX));        trcg_TBLS=0._r8
  allocate(trcn_TBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_TBLS=0._r8
  allocate(trcsa_TBLS(idsa_beg:idsa_end,JS,JY,JX));       trcsa_TBLS=0._r8  
  allocate(trcsa_TQS(idsa_beg:idsa_end,JY,JX));           trcsa_TQS=0._r8  
  allocate(trcn_TQR(ids_nut_beg:ids_nuts_end,JY,JX));     trcn_TQR=0._r8  
  allocate(trcsa_TQR(idsa_beg:idsa_end,JY,JX));           trcsa_TQR=0._r8  
  allocate(trcg_QSS(idg_beg:idg_end-1,JY,JX));            trcg_QSS=0._r8
  allocate(trcn_QSS(ids_nut_beg:ids_nuts_end,JY,JX));trcn_QSS=0._r8
  allocate(trcg_TQR(idg_beg:idg_end-1,JY,JX));      trcg_TQR=0._r8

  allocate(TQI(JY,JX));         TQI=0._r8
  allocate(TQW(JY,JX));         TQW=0._r8
  allocate(TQS(JY,JX));         TQS=0._r8
  allocate(THQS(JY,JX));        THQS=0._r8
  allocate(TKSnow1(JS,JY,JX));     TKSnow1=0._r8
  allocate(THQR1(JY,JX));       THQR1=0._r8
  allocate(THQS1(JY,JX));       THQS1=0._r8
  allocate(TQR1(JY,JX));        TQR1=0._r8
  allocate(TQS1(JY,JX));        TQS1=0._r8
  allocate(TQW1(JY,JX));        TQW1=0._r8
  allocate(TQI1(JY,JX));        TQI1=0._r8
  allocate(VLHeatCapSnowM1(JS,JY,JX));  VLHeatCapSnowM1=0._r8
  allocate(VLIceSnow0M(JS,JY,JX));   VLIceSnow0M=0._r8
  allocate(VLWatSnow0M(JS,JY,JX));   VLWatSnow0M=0._r8
  allocate(CumSno2SnowLay(JS,JY,JX));    CumSno2SnowLay=0._r8
  allocate(CumWat2SnowLay(JS,JY,JX));    CumWat2SnowLay=0._r8
  allocate(CumIce2SnowLay(JS,JY,JX));    CumIce2SnowLay=0._r8
  allocate(CumHeat2SnowLay(JS,JY,JX));   CumHeat2SnowLay=0._r8
  allocate(PhaseChangeHeatL(JS,JY,JX));    PhaseChangeHeatL=0._r8
  allocate(SnowThawMassL(JS,JY,JX));    SnowThawMassL=0._r8
  allocate(IceThawMassL(JS,JY,JX));    IceThawMassL=0._r8
  allocate(WatX2SnoLay(JS,JY,JX));    WatX2SnoLay=0._r8
  allocate(SnoX2SnoLay(JS,JY,JX));    SnoX2SnoLay=0._r8
  allocate(IceX2SnoLay(JS,JY,JX));    IceX2SnoLay=0._r8
  allocate(HeatX2SnoLay(JS,JY,JX));   HeatX2SnoLay=0._r8
  allocate(VLDrySnoWE0M(JS,JY,JX));   VLDrySnoWE0M=0._r8
  allocate(VLDrySnoWE0(JS,JY,JX));    VLDrySnoWE0=0._r8
  allocate(VLIceSnow0(JS,JY,JX));    VLIceSnow0=0._r8
  allocate(VLWatSnow0(JS,JY,JX));    VLWatSnow0=0._r8
  allocate(VLSnoDWI1(JS,JY,JX));    VLSnoDWI1=0._r8
  allocate(SnowLayerThick0(JS,JY,JX));   SnowLayerThick0=0._r8
  allocate(XSnowThawMassL(JS,JY,JX));   XSnowThawMassL=0._r8
  allocate(XIceThawMassL(JS,JY,JX));   XIceThawMassL=0._r8  
  end subroutine  InitSnowPhysData
!----------------------------------------------------------------------  
  subroutine DestructSnowPhysData
  use abortutils, only : destroy
  implicit none

  call destroy(trcg_TBLS)
  call destroy(trcn_TBLS)
  call destroy(trcsa_TBLS)
  call destroy(trcn_TQR)  
  call destroy(trcsa_TQR)
  call destroy(trcg_QSS)
  call destroy(trcn_QSS)  
  call destroy(trcg_TQR)
  call destroy(TQW)  
  call destroy(TQR1)
  call destroy(TQS1)
  call destroy(TQW1)
  call destroy(TQI1)
  call destroy(TQI)  
  call destroy(THQS)  
  call destroy(THQS1)
  call destroy(TKSnow1)  
  call destroy(VLHeatCapSnowM1)
  call destroy(VLIceSnow0M)
  call destroy(VLWatSnow0M)
  call destroy(CumSno2SnowLay)
  call destroy(CumWat2SnowLay)
  call destroy(CumIce2SnowLay)
  call destroy(CumHeat2SnowLay)
  call destroy(PhaseChangeHeatL)
  call destroy(SnowThawMassL)
  call destroy(IceThawMassL)
  call destroy(WatX2SnoLay)
  call destroy(SnoX2SnoLay)
  call destroy(IceX2SnoLay)
  call destroy(HeatX2SnoLay)
  call destroy(VLDrySnoWE0M)
  call destroy(VLDrySnoWE0)
  call destroy(VLIceSnow0)
  call destroy(VLWatSnow0)
  call destroy(VLSnoDWI1)
  call destroy(SnowLayerThick0)
  call destroy(XSnowThawMassL)  
  call destroy(trcsa_TQS)
  call destroy(TQS)
  call destroy(XIceThawMassL)  
  end subroutine DestructSnowPhysData
end module SnowPhysData