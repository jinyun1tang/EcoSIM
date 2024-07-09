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
  real(r8),allocatable ::  trcn_TFloXSurRunoff(:,:,:)                        !
  real(r8),allocatable ::  trcSalt_TQR(:,:,:)                         !
  real(r8),allocatable ::  trcg_QSS(:,:,:)
  real(r8),allocatable ::  trcn_QSS(:,:,:)
  real(r8),allocatable ::  trcg_TFloXSurRunoff(:,:,:)                        !
  real(r8),allocatable ::  XSnowThawMassL(:,:,:)              !hourly convective heat flux from snow transfer
  real(r8),allocatable ::  XIceThawMassL(:,:,:)                      !hourly convective heat flux from ice transfer  
  real(r8),allocatable ::  THeatBySnowRedist_col(:,:)                          !
  real(r8),allocatable ::  TIceBySnowRedist(:,:)                           !
  real(r8),allocatable ::  TWatBySnowRedist(:,:)                           !
  real(r8),allocatable ::  cumWatFlx2LitRByRunoff(:,:)                          !
  real(r8),allocatable ::  cumHeatFlx2LitRByRunoff(:,:)                         !
  real(r8),allocatable ::  cumDrySnoFlxByRedistribut(:,:)                          !
  real(r8),allocatable ::  cumWatFlxBySnowRedistribut(:,:)                          !
  real(r8),allocatable ::  cumIceFlxBySnowRedistribut(:,:)                          !
  real(r8),allocatable ::  cumHeatFlxBySnowRedistribut(:,:)                         !
  real(r8),allocatable ::  TKSnow1(:,:,:)                     ! temporary snow layer temperature in Kelvin
  real(r8),allocatable ::  VLHeatCapSnowM1(:,:,:)             ! temporary snow layer volumetric heat capacity 
  real(r8),allocatable ::  VLIceSnow0M(:,:,:)                 !
  real(r8),allocatable ::  VLWatSnow0M(:,:,:)                 !snow-held water volume during iteration
  real(r8),allocatable ::  CumSno2SnowLayM(:,:,:)
  real(r8),allocatable ::  CumWat2SnowLayM(:,:,:)
  real(r8),allocatable ::  CumIce2SnowLayM(:,:,:)
  real(r8),allocatable ::  CumHeat2SnowLayM(:,:,:)
  real(r8),allocatable ::  XSnowThawMassLM(:,:,:)
  real(r8),allocatable ::  XIceThawMassLM(:,:,:)

  real(r8),allocatable ::  CumSno2SnowLay(:,:,:)                       !
  real(r8),allocatable ::  CumWat2SnowLay(:,:,:)                       !
  real(r8),allocatable ::  CumIce2SnowLay(:,:,:)                       !
  real(r8),allocatable ::  CumHeat2SnowLay(:,:,:)                      !
  real(r8),allocatable ::  XPhaseChangeHeatLM(:,:,:)                       !
  real(r8),allocatable ::  WatX2SnoLay(:,:,:)                       !
  real(r8),allocatable ::  SnoX2SnoLay(:,:,:)                       !
  real(r8),allocatable ::  IceX2SnoLay(:,:,:)                       !
  real(r8),allocatable ::  HeatX2SnoLay(:,:,:)                      !
  real(r8),allocatable ::  VLDrySnoWE0M(:,:,:)                !dry snow layer volume during iteration
  real(r8),allocatable ::  VLDrySnoWE0(:,:,:)                 !dry snow layer volume before update
  real(r8),allocatable ::  VLIceSnow0(:,:,:)                  !snow-held ice layer volume before update
  real(r8),allocatable ::  VLWatSnow0(:,:,:)                       !snow-held water layer volume before update
  real(r8),allocatable ::  VLSnoDWI1(:,:,:)                   !
  real(r8),allocatable ::  SnowThickL_snvr0(:,:,:)                      !
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
  allocate(trcn_TFloXSurRunoff(ids_nut_beg:ids_nuts_end,JY,JX));     trcn_TFloXSurRunoff=0._r8  
  allocate(trcSalt_TQR(idsalt_beg:idsalt_end,JY,JX));           trcSalt_TQR=0._r8  
  allocate(trcg_QSS(idg_beg:idg_end-1,JY,JX));            trcg_QSS=0._r8
  allocate(trcn_QSS(ids_nut_beg:ids_nuts_end,JY,JX));trcn_QSS=0._r8
  allocate(trcg_TFloXSurRunoff(idg_beg:idg_end-1,JY,JX));      trcg_TFloXSurRunoff=0._r8

  allocate(TIceBySnowRedist(JY,JX));         TIceBySnowRedist=0._r8
  allocate(TWatBySnowRedist(JY,JX));         TWatBySnowRedist=0._r8
  allocate(TDrysnoBySnowRedist(JY,JX));         TDrysnoBySnowRedist=0._r8
  allocate(THeatBySnowRedist_col(JY,JX));        THeatBySnowRedist_col=0._r8
  allocate(TKSnow1(JS,JY,JX));     TKSnow1=0._r8
  allocate(cumHeatFlx2LitRByRunoff(JY,JX));       cumHeatFlx2LitRByRunoff=0._r8
  allocate(cumHeatFlxBySnowRedistribut(JY,JX));       cumHeatFlxBySnowRedistribut=0._r8
  allocate(cumWatFlx2LitRByRunoff(JY,JX));        cumWatFlx2LitRByRunoff=0._r8
  allocate(cumDrySnoFlxByRedistribut(JY,JX));        cumDrySnoFlxByRedistribut=0._r8
  allocate(cumWatFlxBySnowRedistribut(JY,JX));        cumWatFlxBySnowRedistribut=0._r8
  allocate(cumIceFlxBySnowRedistribut(JY,JX));        cumIceFlxBySnowRedistribut=0._r8
  allocate(VLHeatCapSnowM1(JS,JY,JX));  VLHeatCapSnowM1=0._r8
  allocate(VLIceSnow0M(JS,JY,JX));   VLIceSnow0M=0._r8
  allocate(VLWatSnow0M(JS,JY,JX));   VLWatSnow0M=0._r8
  allocate(CumSno2SnowLay(JS,JY,JX));    CumSno2SnowLay=0._r8

  allocate(CumSno2SnowLayM(JS,JY,JX)); CumSno2SnowLayM=0._r8
  allocate(CumWat2SnowLayM(JS,JY,JX)); CumWat2SnowLayM=0._r8
  allocate(CumIce2SnowLayM(JS,JY,JX)); CumIce2SnowLayM=0._r8
  allocate(CumHeat2SnowLayM(JS,JY,JX)); CumHeat2SnowLayM=0._r8
  allocate(XSnowThawMassLM(JS,JY,JX)); XSnowThawMassLM=0._r8
  allocate(XIceThawMassLM(JS,JY,JX)); XIceThawMassLM=0._r8

  allocate(CumWat2SnowLay(JS,JY,JX));    CumWat2SnowLay=0._r8
  allocate(CumIce2SnowLay(JS,JY,JX));    CumIce2SnowLay=0._r8
  allocate(CumHeat2SnowLay(JS,JY,JX));   CumHeat2SnowLay=0._r8
  allocate(XPhaseChangeHeatLM(JS,JY,JX));    XPhaseChangeHeatLM=0._r8
  allocate(WatX2SnoLay(JS,JY,JX));    WatX2SnoLay=0._r8
  allocate(SnoX2SnoLay(JS,JY,JX));    SnoX2SnoLay=0._r8
  allocate(IceX2SnoLay(JS,JY,JX));    IceX2SnoLay=0._r8
  allocate(HeatX2SnoLay(JS,JY,JX));   HeatX2SnoLay=0._r8
  allocate(VLDrySnoWE0M(JS,JY,JX));   VLDrySnoWE0M=0._r8
  allocate(VLDrySnoWE0(JS,JY,JX));    VLDrySnoWE0=0._r8
  allocate(VLIceSnow0(JS,JY,JX));    VLIceSnow0=0._r8
  allocate(VLWatSnow0(JS,JY,JX));    VLWatSnow0=0._r8
  allocate(VLSnoDWI1(JS,JY,JX));    VLSnoDWI1=0._r8
  allocate(SnowThickL_snvr0(JS,JY,JX));   SnowThickL_snvr0=0._r8
  allocate(XSnowThawMassL(JS,JY,JX));   XSnowThawMassL=0._r8
  allocate(XIceThawMassL(JS,JY,JX));   XIceThawMassL=0._r8  
  end subroutine  InitSnowPhysData
!----------------------------------------------------------------------  
  subroutine DestructSnowPhysData
  use abortutils, only : destroy
  implicit none

  call destroy(trcg_TBLS)
  call destroy(trcn_TBLS)
  call destroy(trcSalt_TBLS)
  call destroy(trcn_TFloXSurRunoff)  
  call destroy(trcSalt_TQR)
  call destroy(trcg_QSS)
  call destroy(trcn_QSS)  
  call destroy(trcg_TFloXSurRunoff)
  call destroy(TWatBySnowRedist)  
  call destroy(cumWatFlx2LitRByRunoff)
  call destroy(cumDrySnoFlxByRedistribut)
  call destroy(cumWatFlxBySnowRedistribut)
  call destroy(cumIceFlxBySnowRedistribut)
  call destroy(TIceBySnowRedist)  
  call destroy(THeatBySnowRedist_col)  
  call destroy(cumHeatFlxBySnowRedistribut)
  call destroy(TKSnow1)  
  call destroy(VLHeatCapSnowM1)
  call destroy(VLIceSnow0M)
  call destroy(VLWatSnow0M)
  call destroy(CumSno2SnowLay)
  call destroy(CumSno2SnowLayM)
  call destroy(CumWat2SnowLayM)
  call destroy(CumIce2SnowLayM)
  call destroy(CumHeat2SnowLayM)
  call destroy(XSnowThawMassLM)
  call destroy(XIceThawMassLM)
  call destroy(CumWat2SnowLay)
  call destroy(CumIce2SnowLay)
  call destroy(CumHeat2SnowLay)
  call destroy(XPhaseChangeHeatLM)
  call destroy(WatX2SnoLay)
  call destroy(SnoX2SnoLay)
  call destroy(IceX2SnoLay)
  call destroy(HeatX2SnoLay)
  call destroy(VLDrySnoWE0M)
  call destroy(VLDrySnoWE0)
  call destroy(VLIceSnow0)
  call destroy(VLWatSnow0)
  call destroy(VLSnoDWI1)
  call destroy(SnowThickL_snvr0)
  call destroy(XSnowThawMassL)  
  call destroy(trcSalt_TQS)
  call destroy(cumHeatFlx2LitRByRunoff)
  call destroy(TDrysnoBySnowRedist)
  call destroy(XIceThawMassL)  
  end subroutine DestructSnowPhysData
end module SnowPhysData