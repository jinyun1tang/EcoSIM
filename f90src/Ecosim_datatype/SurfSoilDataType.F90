module SurfSoilDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable ::  FracSurfAsSnow_col(:,:)                !fraction of snow cover
  real(r8),target,allocatable ::  FracSurfSnoFree_col(:,:)               !fraction of snow-free cover
  real(r8),target,allocatable ::  FracSurfBareSoil_col(:,:)              !fraction of exposed soil surface, [-]
  real(r8),target,allocatable ::  LWRadBySurf_col(:,:)                   !longwave radiation emitted from ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatByRad2Surf_col(:,:)                !total net radiation at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatEvapAir2Surf_col(:,:)              !total latent heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatSensAir2Surf_col(:,:)              !total sensible heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatSensVapAir2Surf_col(:,:)           !total convective heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatNet2Surf_col(:,:)                  !total ground heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  VapXAir2GSurf_col(:,:)                 !negative of total evaporation at ground surface, [m3 d-2 t-1]
  real(r8),target,allocatable ::  VWatStoreCapSurf_col(:,:)              !surface water storage capacity, [m3 d-2]
  real(r8),target,allocatable ::  VLWatheldCapSurf_col(:,:)              !soil surface water retention capacity
  real(r8),target,allocatable ::  VHCPNX_col(:,:)                        !minimum heat capacities,[MJ k-1 d-2]
  real(r8),target,allocatable ::  CondGasXSnowM_col(:,:,:)               !area upscaled soil surface boundary layer conductance, [m d-2]
  real(r8),target,allocatable ::  Rain2SoilSurf_col(:,:)                 !precipitation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  Irrig2SoilSurf_col(:,:)                !irrifation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  LakeSurfFlowMicP_col(:,:)              !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfFlowMicPX_col(:,:)             !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfFlowMacP_col(:,:)              !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfHeatFlux_col(:,:)              !lake surface heat flux, outgoing positive
!----------------------------------------------------------------------

contains
  subroutine InitSurfSoilData

  implicit none
  allocate(FracSurfAsSnow_col(JY,JX));         FracSurfAsSnow_col              = 0._r8
  allocate(FracSurfSnoFree_col(JY,JX));        FracSurfSnoFree_col             = 0._r8
  allocate(LWRadBySurf_col(JY,JX));             LWRadBySurf_col        = 0._r8
  allocate(HeatByRad2Surf_col(JY,JX));          HeatByRad2Surf_col   = 0._r8
  allocate(HeatEvapAir2Surf_col(JY,JX));       HeatEvapAir2Surf_col    = 0._r8
  allocate(HeatSensAir2Surf_col(JY,JX));       HeatSensAir2Surf_col    = 0._r8
  allocate(HeatSensVapAir2Surf_col(JY,JX));    HeatSensVapAir2Surf_col = 0._r8
  allocate(HeatNet2Surf_col(JY,JX));           HeatNet2Surf_col        = 0._r8
  allocate(VapXAir2GSurf_col(JY,JX));          VapXAir2GSurf_col       = 0._r8
  allocate(FracSurfBareSoil_col(JY,JX));      FracSurfBareSoil_col     = 0._r8
  allocate(VWatStoreCapSurf_col(JY,JX));       VWatStoreCapSurf_col            = 0._r8
  allocate(VLWatheldCapSurf_col(JY,JX));       VLWatheldCapSurf_col        = 0._r8
  allocate(VHCPNX_col(JY,JX));      VHCPNX_col                                 = 0._r8
  allocate(CondGasXSnowM_col(60,JY,JX));     CondGasXSnowM_col                                   = 0._r8
  allocate(Rain2SoilSurf_col(JY,JX));       Rain2SoilSurf_col          = 0._r8
  allocate(Irrig2SoilSurf_col(JY,JX));       Irrig2SoilSurf_col                = 0._r8
  allocate(LakeSurfFlowMicP_col(JY,JX));       LakeSurfFlowMicP_col    = 0._r8
  allocate(LakeSurfFlowMicPX_col(JY,JX));      LakeSurfFlowMicPX_col           = 0._r8
  allocate(LakeSurfFlowMacP_col(JY,JX));      LakeSurfFlowMacP_col     = 0._r8
  allocate(LakeSurfHeatFlux_col(JY,JX));      LakeSurfHeatFlux_col     = 0._r8
  end subroutine InitSurfSoilData

!----------------------------------------------------------------------
  subroutine DestructSurfSoilData
  use abortutils, only : destroy
  implicit none
  call destroy(FracSurfAsSnow_col)
  call destroy(FracSurfSnoFree_col)
  call destroy(LWRadBySurf_col)
  call destroy(HeatByRad2Surf_col)
  call destroy(HeatEvapAir2Surf_col)
  call destroy(HeatSensAir2Surf_col)
  call destroy(HeatSensVapAir2Surf_col)
  call destroy(HeatNet2Surf_col)
  call destroy(VapXAir2GSurf_col)
  call destroy(FracSurfBareSoil_col)
  call destroy(VWatStoreCapSurf_col)
  call destroy(VLWatheldCapSurf_col)
  call destroy(VHCPNX_col)
  call destroy(CondGasXSnowM_col)
  call destroy(Rain2SoilSurf_col)
  call destroy(Irrig2SoilSurf_col)
  call destroy(LakeSurfFlowMicP_col)
  call destroy(LakeSurfFlowMicPX_col)
  call destroy(LakeSurfFlowMacP_col)
  call destroy(LakeSurfHeatFlux_col)
  end subroutine DestructSurfSoilData

end module SurfSoilDataType
