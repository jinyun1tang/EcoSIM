module SnowDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMCtrlMod, only : salt_model
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target, allocatable ::  VLSnowHeatCapM(:,:,:,:)            !volumetric heat capacity of snowpack
  real(r8),target, allocatable ::  WatFlowInSnowM(:,:,:,:)            !snowpack water flux
  real(r8),target, allocatable ::  DrySnoFlxBySnoRedistM(:,:,:,:)                       !runoff snow flux, [m3 d-2 t-1]
  REAL(R8),target, allocatable ::  SoilAlbedo(:,:)                          !snowpack albedo
  real(r8),target, allocatable ::  NewSnowDens(:,:)                   !new snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  TCSnow(:,:,:)                         !snow temperature, [oC]
  real(r8),target, allocatable ::  TKSnow(:,:,:)                         !snow temperature, [K]
  real(r8),target, allocatable ::  VLHeatCapSnow_col(:,:,:)                       !snowpack heat capacity, [MJ m-3 K-1]
  real(r8),target, allocatable ::  VLDrySnoWE_col(:,:,:)                  !water equivalent dry snow in snowpack layer
  real(r8),target, allocatable ::  VLWatSnow_col(:,:,:)                      !snow water volume in snowpack layer
  real(r8),target, allocatable ::  VLIceSnow_col(:,:,:)                      !snow ice volume in snowpack layer
  real(r8),target, allocatable ::  VLSnoDWIprev_col(:,:,:)                       !snow volume in snowpack layer
  real(r8),target, allocatable ::  SnoDensL(:,:,:)                       !snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  SnowThickL_col(:,:,:)                 !snowpack layer thickness
  real(r8),target, allocatable ::  WatXfer2SnoLay(:,:,:)                       !hourly snow water transfer
  real(r8),target, allocatable ::  SnoXfer2SnoLay(:,:,:)                       !hourly snow transfer
  real(r8),target, allocatable ::  IceXfer2SnoLay(:,:,:)                       !hourly snow ice transfer
  real(r8),target, allocatable ::  HeatXfer2SnoLay(:,:,:)                      !hourly convective heat flux from water transfer
  integer ,target, allocatable ::  nsnol_col(:,:)                     !number of snow layers in column
  real(r8),target, allocatable ::  cumSnowDepz_col(:,:,:)                !cumulative depth to bottom of snowpack layer
  real(r8),target, allocatable ::  VLSnoDWIMax_col(:,:,:)                    !Initial snowpack volume, [m3 d-2]
  real(r8),target, allocatable ::  SnowDepth(:,:)                     !snowpack depth, [m]
  real(r8),target, allocatable ::  VcumDrySnoWE_col(:,:)                  !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8),target, allocatable ::  VcumWatSnow_col(:,:)                   !water volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VcumIceSnow_col(:,:)                   !ice volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VcumSnoDWI(:,:)                    !snowpack volume, [m3 d-2]
  real(r8),target, allocatable ::  VcumSnowWE(:,:)                    !water equivalent snowpack [m3 d-2]
  real(r8),target, allocatable ::  VLHeatCapSnowMin_col(:,:)              !minimum layer integrated snowpack heat capacity  [MJ d-2 K-1]
  real(r8),target, allocatable ::  WatConvSno2MicP(:,:,:)                        !water from snowpack to soil micropores
  real(r8),target, allocatable ::  WatConvSno2MacP(:,:,:)             !water from snowpack to soil macropores
  real(r8),target, allocatable ::  HeatConvSno2Soi(:,:,:)                       !convective heat from snowpack to soil
  real(r8),target, allocatable ::  WatConvSno2LitR(:,:,:)                       !water flux from snowpack to litter
  real(r8),target, allocatable ::  HeatConvSno2LitR(:,:,:)                      !convective heat flux from snowpack to litter
  real(r8),target, allocatable ::  DrysnoBySnowRedistrib(:,:,:)                          !snowpack runoff snow, [m3 d-2 h-1]
  real(r8),target, allocatable ::  WatBySnowRedistrib(:,:,:)                          !snowpack runoff water, [m3 d-2 h-1]
  real(r8),target, allocatable ::  IceBySnowRedistrib(:,:,:)                          !snowpack runoff ice, [m3 d-2 h-1]
  real(r8),target, allocatable ::  HeatBySnowRedistrib_2DH(:,:,:)                         !snowpack runoff heat, [MJ d-2 h-1]
  real(r8),target, allocatable ::  trcg_FloXSnow_2DH(:,:,:,:)                      !snowpack runoff CO2 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  trcn_FloXSnow(:,:,:,:)                      !snowpack runoff NH4 flux, [g d-2 h-1]

  real(r8),target, allocatable ::  trcg_solsml(:,:,:,:)               ! snowpack dual phase disolved tracers
  real(r8),target, allocatable ::  trcn_solsml(:,:,:,:)               ! snowpack nutrient dissolved tracers
  real(r8),target, allocatable ::  trcs_solsml(:,:,:,:)               ! snowpack salt dissolved tracers


  real(r8),target, allocatable ::  trcSalt_XQS(:,:,:,:)                       !total salt in snow drift, [mol d-2 h-1]
!----------------------------------------------------------------------

contains
  subroutine InitSnowData

  implicit none
  allocate(VLSnowHeatCapM(60,JS,JY,JX));VLSnowHeatCapM=0._r8
  allocate(WatFlowInSnowM(60,JS,JY,JX)); WatFlowInSnowM=0._r8
  allocate(DrySnoFlxBySnoRedistM(60,2,JV,JH));    DrySnoFlxBySnoRedistM=0._r8
  allocate(SoilAlbedo(JY,JX));        SoilAlbedo=0._r8
  allocate(NewSnowDens(JY,JX));       NewSnowDens=0._r8
  allocate(TCSnow(JS,JY,JX));      TCSnow=0._r8
  allocate(TKSnow(JS,JY,JX));      TKSnow=0._r8
  allocate(VLHeatCapSnow_col(JS,JY,JX));    VLHeatCapSnow_col=0._r8
  allocate(VLDrySnoWE_col(JS,JY,JX));   VLDrySnoWE_col=0._r8
  allocate(VLWatSnow_col(JS,JY,JX));   VLWatSnow_col=0._r8
  allocate(VLIceSnow_col(JS,JY,JX));   VLIceSnow_col=0._r8
  allocate(VLSnoDWIprev_col(JS,JY,JX));    VLSnoDWIprev_col=0._r8
  allocate(SnoDensL(JS,JY,JX));    SnoDensL=0._r8
  allocate(SnowThickL_col(JS,JY,JX));    SnowThickL_col=0._r8
  allocate(WatXfer2SnoLay(JS,JY,JX));    WatXfer2SnoLay=0._r8
  allocate(SnoXfer2SnoLay(JS,JY,JX));    SnoXfer2SnoLay=0._r8
  allocate(IceXfer2SnoLay(JS,JY,JX));    IceXfer2SnoLay=0._r8
  allocate(HeatXfer2SnoLay(JS,JY,JX));   HeatXfer2SnoLay=0._r8
  allocate(nsnol_col(JY,JX)); nsnol_col=0
  allocate(cumSnowDepz_col(0:JS,JY,JX)); cumSnowDepz_col=0._r8
  allocate(VLSnoDWIMax_col(JS,JY,JX));    VLSnoDWIMax_col=0._r8
  allocate(SnowDepth(JY,JX));       SnowDepth=0._r8
  allocate(VcumSnowWE(JY,JX));   VcumSnowWE=0._r8
  allocate(VcumDrySnoWE_col(JY,JX));       VcumDrySnoWE_col=0._r8
  allocate(VcumWatSnow_col(JY,JX));       VcumWatSnow_col=0._r8
  allocate(VcumIceSnow_col(JY,JX));       VcumIceSnow_col=0._r8
  allocate(VcumSnoDWI(JY,JX));        VcumSnoDWI=0._r8
  allocate(VLHeatCapSnowMin_col(JY,JX));      VLHeatCapSnowMin_col=0._r8
  allocate(WatConvSno2MicP(JS,JY,JX));     WatConvSno2MicP=0._r8
  allocate(WatConvSno2MacP(JS,JY,JX));    WatConvSno2MacP=0._r8
  allocate(HeatConvSno2Soi(JS,JY,JX));    HeatConvSno2Soi=0._r8
  allocate(WatConvSno2LitR(JS,JY,JX));    WatConvSno2LitR=0._r8
  allocate(HeatConvSno2LitR(JS,JY,JX));   HeatConvSno2LitR=0._r8
  allocate(DrysnoBySnowRedistrib(2,JV,JH));        DrysnoBySnowRedistrib=0._r8
  allocate(WatBySnowRedistrib(2,JV,JH));        WatBySnowRedistrib=0._r8
  allocate(IceBySnowRedistrib(2,JV,JH));        IceBySnowRedistrib=0._r8
  allocate(HeatBySnowRedistrib_2DH(2,JV,JH));       HeatBySnowRedistrib_2DH=0._r8
  allocate(trcg_FloXSnow_2DH(idg_beg:idg_NH3,2,JV,JH));    trcg_FloXSnow_2DH=0._r8

  allocate(trcn_FloXSnow(ids_nut_beg:ids_nuts_end,2,JV,JH));    trcn_FloXSnow=0._r8

! exclude NH3B
  allocate(trcg_solsml(idg_beg:idg_NH3,JS,JY,JX));trcg_solsml=0._r8
  allocate(trcn_solsml(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml=0._r8
  if(salt_model)then
    allocate(trcs_solsml(idsalt_beg:idsalt_end,JS,JY,JX)); trcs_solsml=0._r8
  endif

  allocate(trcSalt_XQS(idsalt_beg:idsalt_end,2,JV,JH));     trcSalt_XQS=0._r8
  end subroutine InitSnowData

!----------------------------------------------------------------------
  subroutine DestructSnowData
  use abortutils, only : destroy
  implicit none
  if(salt_model)then
    call destroy(trcs_solsml)
    call destroy(trcSalt_XQS)
  endif
  call destroy(VLSnowHeatCapM)
  call destroy(WatFlowInSnowM)
  call destroy(DrySnoFlxBySnoRedistM)
  call destroy(SoilAlbedo)
  call destroy(NewSnowDens)
  call destroy(TCSnow)
  call destroy(TKSnow)
  call destroy(VLHeatCapSnow_col)
  call destroy(VLDrySnoWE_col)
  call destroy(VLWatSnow_col)
  call destroy(VLIceSnow_col)
  call destroy(VLSnoDWIprev_col)
  call destroy(SnoDensL)
  call destroy(SnowThickL_col)
  call destroy(WatXfer2SnoLay)
  call destroy(SnoXfer2SnoLay)
  call destroy(IceXfer2SnoLay)
  call destroy(HeatXfer2SnoLay)
  call destroy(nsnol_col)
  call destroy(cumSnowDepz_col)
  call destroy(VLSnoDWIMax_col)
  call destroy(SnowDepth)
  call destroy(VcumSnowWE)
  call destroy(VcumDrySnoWE_col)
  call destroy(VcumWatSnow_col)
  call destroy(VcumIceSnow_col)
  call destroy(VcumSnoDWI)
  call destroy(VLHeatCapSnowMin_col)
  call destroy(WatConvSno2MicP)
  call destroy(WatConvSno2MacP)
  call destroy(HeatConvSno2Soi)
  call destroy(WatConvSno2LitR)
  call destroy(HeatConvSno2LitR)
  call destroy(DrysnoBySnowRedistrib)
  call destroy(WatBySnowRedistrib)
  call destroy(IceBySnowRedistrib)
  call destroy(HeatBySnowRedistrib_2DH)
  call destroy(trcg_FloXSnow_2DH)
  call destroy(trcn_FloXSnow)
  if(salt_model)then
    call destroy(trcSalt_XQS)
  endif
  end subroutine DestructSnowData

end module SnowDataType
