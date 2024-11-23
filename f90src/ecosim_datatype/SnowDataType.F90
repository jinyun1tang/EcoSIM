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

  real(r8),target, allocatable ::  VLSnowHeatCapM_snvr(:,:,:,:)            !volumetric heat capacity of snowpack
  real(r8),target, allocatable ::  WatFlowInSnowM_snvr(:,:,:,:)            !snowpack water flux
  real(r8),target, allocatable ::  DrySnoFlxBySnoRedistM_2DH(:,:,:,:)                       !runoff snow flux, [m3 d-2 t-1]
  REAL(R8),target, allocatable ::  SoilAlbedo_col(:,:)                          !snowpack albedo
  real(r8),target, allocatable ::  NewSnowDens_col(:,:)                   !new snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  TCSnow_snvr(:,:,:)                         !snow temperature, [oC]
  real(r8),target, allocatable ::  TKSnow_snvr(:,:,:)                         !snow temperature, [K]
  real(r8),target, allocatable ::  VLHeatCapSnow_snvr(:,:,:)                       !snowpack heat capacity, [MJ m-3 K-1]
  real(r8),target, allocatable ::  VLDrySnoWE_snvr(:,:,:)                  !water equivalent dry snow in snowpack layer
  real(r8),target, allocatable ::  VLWatSnow_snvr(:,:,:)                      !snow water volume in snowpack layer
  real(r8),target, allocatable ::  VLIceSnow_snvr(:,:,:)                      !snow ice volume in snowpack layer
  real(r8),target, allocatable ::  VLSnoDWIprev_snvr(:,:,:)                       !snow volume in snowpack layer
  real(r8),target, allocatable ::  SnoDens_snvr(:,:,:)                       !snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  SnowThickL_snvr(:,:,:)                 !snowpack layer thickness
  real(r8),target, allocatable ::  WatXfer2SnoLay_snvr(:,:,:)                       !hourly snow water transfer
  real(r8),target, allocatable ::  SnoXfer2SnoLay_snvr(:,:,:)                       !hourly snow transfer
  real(r8),target, allocatable ::  IceXfer2SnoLay_snvr(:,:,:)                       !hourly snow ice transfer
  real(r8),target, allocatable ::  HeatXfer2SnoLay_snvr(:,:,:)                      !hourly convective heat flux from water transfer
  integer ,target, allocatable ::  nsnol_col(:,:)                     !number of snow layers in column
  real(r8),target, allocatable ::  cumSnowDepz_col(:,:,:)                !cumulative depth to bottom of snowpack layer
  real(r8),target, allocatable ::  VLSnoDWIMax_col(:,:,:)                    !Initial snowpack volume, [m3 d-2]
  real(r8),target, allocatable ::  SnowDepth_col(:,:)                     !snowpack depth, [m]
  real(r8),target, allocatable ::  VcumDrySnoWE_col(:,:)                  !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8),target, allocatable ::  VcumWatSnow_col(:,:)                   !water volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VcumIceSnow_col(:,:)                   !ice volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VcumSnoDWI_col(:,:)                    !snowpack volume, [m3 d-2]
  real(r8),target, allocatable ::  VcumSnowWE_col(:,:)                    !water equivalent snowpack [m3 d-2]
  real(r8),target, allocatable ::  VLHeatCapSnowMin_col(:,:)              !minimum layer integrated snowpack heat capacity  [MJ d-2 K-1]
  real(r8),target, allocatable ::  WatConvSno2MicP_snvr(:,:,:)                        !water from snowpack to soil micropores
  real(r8),target, allocatable ::  WatConvSno2MacP_snvr(:,:,:)             !water from snowpack to soil macropores
  real(r8),target, allocatable ::  HeatConvSno2Soi_snvr(:,:,:)                       !convective heat from snowpack to soil
  real(r8),target, allocatable ::  WatConvSno2LitR_snvr(:,:,:)                       !water flux from snowpack to litter
  real(r8),target, allocatable ::  HeatConvSno2LitR_snvr(:,:,:)                      !convective heat flux from snowpack to litter
  real(r8),target, allocatable ::  DrysnoBySnowRedistrib(:,:,:)                          !snowpack runoff snow, [m3 d-2 h-1]
  real(r8),target, allocatable ::  WatBySnowRedistrib_2DH(:,:,:)                          !snowpack runoff water, [m3 d-2 h-1]
  real(r8),target, allocatable ::  IceBySnowRedistrib_2DH(:,:,:)                          !snowpack runoff ice, [m3 d-2 h-1]
  real(r8),target, allocatable ::  HeatBySnowRedistrib_2DH(:,:,:)                         !snowpack runoff heat, [MJ d-2 h-1]
  real(r8),target, allocatable ::  trcg_FloXSnow_2DH(:,:,:,:)                      !snowpack runoff CO2 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  trcn_FloXSnow_2DH(:,:,:,:)                      !snowpack runoff NH4 flux, [g d-2 h-1]

  real(r8),target, allocatable ::  trcg_solsml_snvr(:,:,:,:)               ! Disolved volatile tracers in snow [g d-2]
  real(r8),target, allocatable ::  trcn_solsml_snvr(:,:,:,:)               ! Dissolved nutrient tracers in snow [g d-2]
  real(r8),target, allocatable ::  trc_Saltml_snvr(:,:,:,:)                    ! snowpack salt dissolved tracers


  real(r8),target, allocatable ::  trcSalt_XQS(:,:,:,:)                       !total salt in snow drift, [mol d-2 h-1]
!----------------------------------------------------------------------

contains
  subroutine InitSnowData

  implicit none
  allocate(VLSnowHeatCapM_snvr(60,JS,JY,JX));VLSnowHeatCapM_snvr=0._r8
  allocate(WatFlowInSnowM_snvr(60,JS,JY,JX)); WatFlowInSnowM_snvr=0._r8
  allocate(DrySnoFlxBySnoRedistM_2DH(60,2,JV,JH));    DrySnoFlxBySnoRedistM_2DH=0._r8
  allocate(SoilAlbedo_col(JY,JX));        SoilAlbedo_col=0._r8
  allocate(NewSnowDens_col(JY,JX));       NewSnowDens_col=0._r8
  allocate(TCSnow_snvr(JS,JY,JX));      TCSnow_snvr=0._r8
  allocate(TKSnow_snvr(JS,JY,JX));      TKSnow_snvr=0._r8
  allocate(VLHeatCapSnow_snvr(JS,JY,JX));    VLHeatCapSnow_snvr=0._r8
  allocate(VLDrySnoWE_snvr(JS,JY,JX));   VLDrySnoWE_snvr=0._r8
  allocate(VLWatSnow_snvr(JS,JY,JX));   VLWatSnow_snvr=0._r8
  allocate(VLIceSnow_snvr(JS,JY,JX));   VLIceSnow_snvr=0._r8
  allocate(VLSnoDWIprev_snvr(JS,JY,JX));    VLSnoDWIprev_snvr=0._r8
  allocate(SnoDens_snvr(JS,JY,JX));    SnoDens_snvr=0._r8
  allocate(SnowThickL_snvr(JS,JY,JX));    SnowThickL_snvr=0._r8
  allocate(WatXfer2SnoLay_snvr(JS,JY,JX));    WatXfer2SnoLay_snvr=0._r8
  allocate(SnoXfer2SnoLay_snvr(JS,JY,JX));    SnoXfer2SnoLay_snvr=0._r8
  allocate(IceXfer2SnoLay_snvr(JS,JY,JX));    IceXfer2SnoLay_snvr=0._r8
  allocate(HeatXfer2SnoLay_snvr(JS,JY,JX));   HeatXfer2SnoLay_snvr=0._r8
  allocate(nsnol_col(JY,JX)); nsnol_col=0
  allocate(cumSnowDepz_col(0:JS,JY,JX)); cumSnowDepz_col=0._r8
  allocate(VLSnoDWIMax_col(JS,JY,JX));    VLSnoDWIMax_col=0._r8
  allocate(SnowDepth_col(JY,JX));       SnowDepth_col=0._r8
  allocate(VcumSnowWE_col(JY,JX));   VcumSnowWE_col=0._r8
  allocate(VcumDrySnoWE_col(JY,JX));       VcumDrySnoWE_col=0._r8
  allocate(VcumWatSnow_col(JY,JX));       VcumWatSnow_col=0._r8
  allocate(VcumIceSnow_col(JY,JX));       VcumIceSnow_col=0._r8
  allocate(VcumSnoDWI_col(JY,JX));        VcumSnoDWI_col=0._r8
  allocate(VLHeatCapSnowMin_col(JY,JX));      VLHeatCapSnowMin_col=0._r8
  allocate(WatConvSno2MicP_snvr(JS,JY,JX));     WatConvSno2MicP_snvr=0._r8
  allocate(WatConvSno2MacP_snvr(JS,JY,JX));    WatConvSno2MacP_snvr=0._r8
  allocate(HeatConvSno2Soi_snvr(JS,JY,JX));    HeatConvSno2Soi_snvr=0._r8
  allocate(WatConvSno2LitR_snvr(JS,JY,JX));    WatConvSno2LitR_snvr=0._r8
  allocate(HeatConvSno2LitR_snvr(JS,JY,JX));   HeatConvSno2LitR_snvr=0._r8
  allocate(DrysnoBySnowRedistrib(2,JV,JH));        DrysnoBySnowRedistrib=0._r8
  allocate(WatBySnowRedistrib_2DH(2,JV,JH));        WatBySnowRedistrib_2DH=0._r8
  allocate(IceBySnowRedistrib_2DH(2,JV,JH));        IceBySnowRedistrib_2DH=0._r8
  allocate(HeatBySnowRedistrib_2DH(2,JV,JH));       HeatBySnowRedistrib_2DH=0._r8
  allocate(trcg_FloXSnow_2DH(idg_beg:idg_NH3,2,JV,JH));    trcg_FloXSnow_2DH=0._r8

  allocate(trcn_FloXSnow_2DH(ids_nut_beg:ids_nuts_end,2,JV,JH));    trcn_FloXSnow_2DH=0._r8

! exclude NH3B
  allocate(trcg_solsml_snvr(idg_beg:idg_NH3,JS,JY,JX));trcg_solsml_snvr=0._r8
  allocate(trcn_solsml_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml_snvr=0._r8
  if(salt_model)then
    allocate(trc_Saltml_snvr(idsalt_beg:idsalt_end,JS,JY,JX)); trc_Saltml_snvr=0._r8
  endif

  allocate(trcSalt_XQS(idsalt_beg:idsalt_end,2,JV,JH));     trcSalt_XQS=0._r8
  end subroutine InitSnowData

!----------------------------------------------------------------------
  subroutine DestructSnowData
  use abortutils, only : destroy
  implicit none
  if(salt_model)then
    call destroy(trc_Saltml_snvr)
    call destroy(trcSalt_XQS)
  endif
  call destroy(VLSnowHeatCapM_snvr)
  call destroy(WatFlowInSnowM_snvr)
  call destroy(DrySnoFlxBySnoRedistM_2DH)
  call destroy(SoilAlbedo_col)
  call destroy(NewSnowDens_col)
  call destroy(trcg_solsml_snvr)
  call destroy(TCSnow_snvr)
  call destroy(TKSnow_snvr)
  call destroy(VLHeatCapSnow_snvr)
  call destroy(VLDrySnoWE_snvr)
  call destroy(VLWatSnow_snvr)
  call destroy(VLIceSnow_snvr)
  call destroy(VLSnoDWIprev_snvr)
  call destroy(SnoDens_snvr)
  call destroy(SnowThickL_snvr)
  call destroy(WatXfer2SnoLay_snvr)
  call destroy(SnoXfer2SnoLay_snvr)
  call destroy(IceXfer2SnoLay_snvr)
  call destroy(HeatXfer2SnoLay_snvr)
  call destroy(nsnol_col)
  call destroy(cumSnowDepz_col)
  call destroy(VLSnoDWIMax_col)
  call destroy(SnowDepth_col)
  call destroy(VcumSnowWE_col)
  call destroy(VcumDrySnoWE_col)
  call destroy(VcumWatSnow_col)
  call destroy(VcumIceSnow_col)
  call destroy(VcumSnoDWI_col)
  call destroy(VLHeatCapSnowMin_col)
  call destroy(WatConvSno2MicP_snvr)
  call destroy(WatConvSno2MacP_snvr)
  call destroy(HeatConvSno2Soi_snvr)
  call destroy(WatConvSno2LitR_snvr)
  call destroy(HeatConvSno2LitR_snvr)
  call destroy(DrysnoBySnowRedistrib)
  call destroy(WatBySnowRedistrib_2DH)
  call destroy(IceBySnowRedistrib_2DH)
  call destroy(HeatBySnowRedistrib_2DH)
  call destroy(trcg_FloXSnow_2DH)
  call destroy(trcn_FloXSnow_2DH)
  call destroy(trcn_solsml_snvr)
  if(salt_model)then
    call destroy(trcSalt_XQS)
  endif
  end subroutine DestructSnowData

end module SnowDataType
