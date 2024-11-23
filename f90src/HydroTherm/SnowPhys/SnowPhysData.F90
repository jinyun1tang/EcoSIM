module SnowPhysData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod

  implicit none

  character(len=*), private, parameter :: mod_filename=&
  __FILE__

  real(r8),allocatable ::  trcg_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcn_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcSalt_TBLS(:,:,:,:)                      !
  real(r8),allocatable ::  TDrysnoBySnowRedist(:,:)                           !
  real(r8),allocatable ::  trcSalt_TQS(:,:,:)  
  real(r8),allocatable ::  trcn_TFloXSurRunoff_2D(:,:,:)                        !
  real(r8),allocatable ::  trcSalt_TQR(:,:,:)                         !
  real(r8),allocatable ::  trcg_QSS(:,:,:)
  real(r8),allocatable ::  trcn_QSS(:,:,:)
  real(r8),allocatable ::  trcg_TFloXSurRunoff(:,:,:)                        !
  real(r8),allocatable ::  XSnowThawMassL_snvr(:,:,:)              !hourly convective heat flux from snow transfer
  real(r8),allocatable ::  XIceThawMassL_snvr(:,:,:)                      !hourly convective heat flux from ice transfer  
  real(r8),allocatable ::  THeatBySnowRedist_col(:,:)                          !
  real(r8),allocatable ::  TIceBySnowRedist(:,:)                           !
  real(r8),allocatable ::  TWatBySnowRedist(:,:)                           !
  real(r8),allocatable ::  cumWatFlx2LitRByRunoff_col(:,:)                          !
  real(r8),allocatable ::  cumHeatFlx2LitRByRunoff_col(:,:)                         !
  real(r8),allocatable ::  cumDrySnoFlxByRedistribut(:,:)                          !
  real(r8),allocatable ::  cumWatFlxBySnowRedistribut(:,:)                          !
  real(r8),allocatable ::  cumIceFlxBySnowRedistribut(:,:)                          !
  real(r8),allocatable ::  cumHeatFlxBySnowRedistribut(:,:)                         !
  real(r8),allocatable ::  TKSnow1_snvr(:,:,:)                     ! temporary snow layer temperature in Kelvin
  real(r8),allocatable ::  VLHeatCapSnowM1_snvr(:,:,:)             ! temporary snow layer volumetric heat capacity 
  real(r8),allocatable ::  VLIceSnow0M_snvr(:,:,:)                 !
  real(r8),allocatable ::  VLWatSnow0M_snvr(:,:,:)                 !snow-held water volume during iteration
  real(r8),allocatable ::  CumSno2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  CumWat2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  CumIce2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  CumHeat2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  XSnowThawMassLM_snvr(:,:,:)
  real(r8),allocatable ::  XIceThawMassLM_snvr(:,:,:)
  real(r8),allocatable ::  tEnGYM_snvr(:,:,:)  
  real(r8),allocatable ::  CumSno2SnowL_snvr(:,:,:)                       !
  real(r8),allocatable ::  CumWat2SnowL_snvr(:,:,:)                       !
  real(r8),allocatable ::  CumIce2SnowL_snvr(:,:,:)                       !
  real(r8),allocatable ::  CumHeat2SnowL_snvr(:,:,:)                      !
  real(r8),allocatable ::  XPhaseChangeHeatLM_snvr(:,:,:)                       !
  real(r8),allocatable ::  WatX2SnoLay_snvr(:,:,:)                       !
  real(r8),allocatable ::  SnoX2SnoLay_snvr(:,:,:)                       !
  real(r8),allocatable ::  IceX2SnoLay_snvr(:,:,:)                       !
  real(r8),allocatable ::  HeatX2SnoLay_snvr(:,:,:)                      !
  real(r8),allocatable ::  VLDrySnoWE0M_snvr(:,:,:)                !dry snow layer volume during iteration
  real(r8),allocatable ::  VLDrySnoWE0_snvr(:,:,:)                 !dry snow layer volume before update
  real(r8),allocatable ::  VLIceSnow0_snvr(:,:,:)                  !snow-held ice layer volume before update
  real(r8),allocatable ::  VLWatSnow0_snvr(:,:,:)                       !snow-held water layer volume before update
  real(r8),allocatable ::  VLSnoDWI1_snvr(:,:,:)                   !
  real(r8),allocatable ::  SnowThickL0_snvr(:,:,:)                      !
  public :: InitSnowPhysData
  public :: DestructSnowPhysData
  contains
!----------------------------------------------------------------------  
  subroutine  InitSnowPhysData
  implicit none

  allocate(trcg_TBLS(idg_beg:idg_end-1,JS,JY,JX));        trcg_TBLS=0._r8
  allocate(trcn_TBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_TBLS=0._r8
  allocate(trcSalt_TBLS(idsalt_beg:idsalt_end,JS,JY,JX));       trcSalt_TBLS=0._r8  
  allocate(trcSalt_TQS(idsalt_beg:idsalt_end,JY,JX));           trcSalt_TQS=0._r8  
  allocate(trcn_TFloXSurRunoff_2D(ids_nut_beg:ids_nuts_end,JY,JX));     trcn_TFloXSurRunoff_2D=0._r8  
  allocate(trcSalt_TQR(idsalt_beg:idsalt_end,JY,JX));           trcSalt_TQR=0._r8  
  allocate(trcg_QSS(idg_beg:idg_end-1,JY,JX));            trcg_QSS=0._r8
  allocate(trcn_QSS(ids_nut_beg:ids_nuts_end,JY,JX));trcn_QSS=0._r8
  allocate(trcg_TFloXSurRunoff(idg_beg:idg_end-1,JY,JX));      trcg_TFloXSurRunoff=0._r8

  allocate(TIceBySnowRedist(JY,JX));         TIceBySnowRedist=0._r8
  allocate(TWatBySnowRedist(JY,JX));         TWatBySnowRedist=0._r8
  allocate(TDrysnoBySnowRedist(JY,JX));         TDrysnoBySnowRedist=0._r8
  allocate(THeatBySnowRedist_col(JY,JX));        THeatBySnowRedist_col=0._r8
  allocate(TKSnow1_snvr(JS,JY,JX));     TKSnow1_snvr=0._r8
  allocate(cumHeatFlx2LitRByRunoff_col(JY,JX));       cumHeatFlx2LitRByRunoff_col=0._r8
  allocate(cumHeatFlxBySnowRedistribut(JY,JX));       cumHeatFlxBySnowRedistribut=0._r8
  allocate(cumWatFlx2LitRByRunoff_col(JY,JX));        cumWatFlx2LitRByRunoff_col=0._r8
  allocate(cumDrySnoFlxByRedistribut(JY,JX));        cumDrySnoFlxByRedistribut=0._r8
  allocate(cumWatFlxBySnowRedistribut(JY,JX));        cumWatFlxBySnowRedistribut=0._r8
  allocate(cumIceFlxBySnowRedistribut(JY,JX));        cumIceFlxBySnowRedistribut=0._r8
  allocate(VLHeatCapSnowM1_snvr(JS,JY,JX));  VLHeatCapSnowM1_snvr=0._r8
  allocate(VLIceSnow0M_snvr(JS,JY,JX));   VLIceSnow0M_snvr=0._r8
  allocate(VLWatSnow0M_snvr(JS,JY,JX));   VLWatSnow0M_snvr=0._r8
  allocate(CumSno2SnowL_snvr(JS,JY,JX));    CumSno2SnowL_snvr=0._r8
  allocate(tEnGYM_snvr(JS,JY,JX));  tEnGYM_snvr=0._r8
  allocate(CumSno2SnowLM_snvr(JS,JY,JX)); CumSno2SnowLM_snvr=0._r8
  allocate(CumWat2SnowLM_snvr(JS,JY,JX)); CumWat2SnowLM_snvr=0._r8
  allocate(CumIce2SnowLM_snvr(JS,JY,JX)); CumIce2SnowLM_snvr=0._r8
  allocate(CumHeat2SnowLM_snvr(JS,JY,JX)); CumHeat2SnowLM_snvr=0._r8
  allocate(XSnowThawMassLM_snvr(JS,JY,JX)); XSnowThawMassLM_snvr=0._r8
  allocate(XIceThawMassLM_snvr(JS,JY,JX)); XIceThawMassLM_snvr=0._r8

  allocate(CumWat2SnowL_snvr(JS,JY,JX));    CumWat2SnowL_snvr=0._r8
  allocate(CumIce2SnowL_snvr(JS,JY,JX));    CumIce2SnowL_snvr=0._r8
  allocate(CumHeat2SnowL_snvr(JS,JY,JX));   CumHeat2SnowL_snvr=0._r8
  allocate(XPhaseChangeHeatLM_snvr(JS,JY,JX));    XPhaseChangeHeatLM_snvr=0._r8
  allocate(WatX2SnoLay_snvr(JS,JY,JX));    WatX2SnoLay_snvr=0._r8
  allocate(SnoX2SnoLay_snvr(JS,JY,JX));    SnoX2SnoLay_snvr=0._r8
  allocate(IceX2SnoLay_snvr(JS,JY,JX));    IceX2SnoLay_snvr=0._r8
  allocate(HeatX2SnoLay_snvr(JS,JY,JX));   HeatX2SnoLay_snvr=0._r8
  allocate(VLDrySnoWE0M_snvr(JS,JY,JX));   VLDrySnoWE0M_snvr=0._r8
  allocate(VLDrySnoWE0_snvr(JS,JY,JX));    VLDrySnoWE0_snvr=0._r8
  allocate(VLIceSnow0_snvr(JS,JY,JX));    VLIceSnow0_snvr=0._r8
  allocate(VLWatSnow0_snvr(JS,JY,JX));    VLWatSnow0_snvr=0._r8
  allocate(VLSnoDWI1_snvr(JS,JY,JX));    VLSnoDWI1_snvr=0._r8
  allocate(SnowThickL0_snvr(JS,JY,JX));   SnowThickL0_snvr=0._r8
  allocate(XSnowThawMassL_snvr(JS,JY,JX));   XSnowThawMassL_snvr=0._r8
  allocate(XIceThawMassL_snvr(JS,JY,JX));   XIceThawMassL_snvr=0._r8  
  end subroutine  InitSnowPhysData
!----------------------------------------------------------------------  
  subroutine DestructSnowPhysData
  use abortutils, only : destroy
  implicit none

  call destroy(trcg_TBLS)
  call destroy(trcn_TBLS)
  call destroy(trcSalt_TBLS)
  call destroy(trcn_TFloXSurRunoff_2D)  
  call destroy(trcSalt_TQR)
  call destroy(tEnGYM_snvr)
  call destroy(trcg_QSS)
  call destroy(trcn_QSS)  
  call destroy(trcg_TFloXSurRunoff)
  call destroy(TWatBySnowRedist)  
  call destroy(cumWatFlx2LitRByRunoff_col)
  call destroy(cumDrySnoFlxByRedistribut)
  call destroy(cumWatFlxBySnowRedistribut)
  call destroy(cumIceFlxBySnowRedistribut)
  call destroy(TIceBySnowRedist)  
  call destroy(THeatBySnowRedist_col)  
  call destroy(cumHeatFlxBySnowRedistribut)
  call destroy(TKSnow1_snvr)  
  call destroy(VLHeatCapSnowM1_snvr)
  call destroy(VLIceSnow0M_snvr)
  call destroy(VLWatSnow0M_snvr)
  call destroy(CumSno2SnowL_snvr)
  call destroy(CumSno2SnowLM_snvr)
  call destroy(CumWat2SnowLM_snvr)
  call destroy(CumIce2SnowLM_snvr)
  call destroy(CumHeat2SnowLM_snvr)
  call destroy(XSnowThawMassLM_snvr)
  call destroy(XIceThawMassLM_snvr)
  call destroy(CumWat2SnowL_snvr)
  call destroy(CumIce2SnowL_snvr)
  call destroy(CumHeat2SnowL_snvr)
  call destroy(XPhaseChangeHeatLM_snvr)
  call destroy(WatX2SnoLay_snvr)
  call destroy(SnoX2SnoLay_snvr)
  call destroy(IceX2SnoLay_snvr)
  call destroy(HeatX2SnoLay_snvr)
  call destroy(VLDrySnoWE0M_snvr)
  call destroy(VLDrySnoWE0_snvr)
  call destroy(VLIceSnow0_snvr)
  call destroy(VLWatSnow0_snvr)
  call destroy(VLSnoDWI1_snvr)
  call destroy(SnowThickL0_snvr)
  call destroy(XSnowThawMassL_snvr)  
  call destroy(trcSalt_TQS)
  call destroy(cumHeatFlx2LitRByRunoff_col)
  call destroy(TDrysnoBySnowRedist)
  call destroy(XIceThawMassL_snvr)  
  end subroutine DestructSnowPhysData
end module SnowPhysData