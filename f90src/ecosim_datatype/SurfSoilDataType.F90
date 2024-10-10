module SurfSoilDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable ::  FracSurfAsSnow(:,:)                !fraction of snow cover
  real(r8),target,allocatable ::  FracSurfSnoFree(:,:)               !fraction of snow-free cover
  real(r8),target,allocatable ::  FracSurfBareSoil_col(:,:)             !fraction of exposed soil surface, [-]
  real(r8),target,allocatable ::  LWRadBySurf_col(:,:)                         !longwave radiation emitted from ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatByRadiation_col(:,:)                 !total net radiation at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatEvapAir2Surf_col(:,:)              !total latent heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatSensAir2Surf_col(:,:)              !total sensible heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatSensVapAir2Surf_col(:,:)           !total convective heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatNet2Surf_col(:,:)                  !total ground heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  VapXAir2GSurf_col(:,:)                 !negative of total evaporation at ground surface, [m3 d-2 t-1]
  real(r8),target,allocatable ::  VWatStoreCapSurf(:,:)                         !surface water storage capacity, [m3 d-2]
  real(r8),target,allocatable ::  MaxVLWatByLitR_col(:,:)                         !soil surface water retention capacity
  real(r8),target,allocatable ::  VHCPNX(:,:)                          !minimum heat capacities
  real(r8),target,allocatable ::  PARG(:,:,:)                          !soil surface boundary layer conductance, [m t-1]
  real(r8),target,allocatable ::  Rain2SoilSurf_col(:,:)               !precipitation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  Irrig2SoilSurf(:,:)                  !irrifation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  LakeSurfFlowMicP_col(:,:)                   !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfFlowMicPX(:,:)                      !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfFlowMacP_col(:,:)                       !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfHeatFlux_col(:,:)              !lake surface heat flux, outgoing positive
!----------------------------------------------------------------------

contains
  subroutine InitSurfSoilData

  implicit none
  allocate(FracSurfAsSnow(JY,JX));         FracSurfAsSnow              = 0._r8
  allocate(FracSurfSnoFree(JY,JX));        FracSurfSnoFree             = 0._r8
  allocate(LWRadBySurf_col(JY,JX));             LWRadBySurf_col        = 0._r8
  allocate(HeatByRadiation_col(JY,JX));          HeatByRadiation_col   = 0._r8
  allocate(HeatEvapAir2Surf_col(JY,JX));       HeatEvapAir2Surf_col    = 0._r8
  allocate(HeatSensAir2Surf_col(JY,JX));       HeatSensAir2Surf_col    = 0._r8
  allocate(HeatSensVapAir2Surf_col(JY,JX));    HeatSensVapAir2Surf_col = 0._r8
  allocate(HeatNet2Surf_col(JY,JX));           HeatNet2Surf_col        = 0._r8
  allocate(VapXAir2GSurf_col(JY,JX));          VapXAir2GSurf_col       = 0._r8
  allocate(FracSurfBareSoil_col(JY,JX));      FracSurfBareSoil_col     = 0._r8
  allocate(VWatStoreCapSurf(JY,JX));       VWatStoreCapSurf            = 0._r8
  allocate(MaxVLWatByLitR_col(JY,JX));       MaxVLWatByLitR_col        = 0._r8
  allocate(VHCPNX(JY,JX));      VHCPNX                                 = 0._r8
  allocate(PARG(60,JY,JX));     PARG                                   = 0._r8
  allocate(Rain2SoilSurf_col(JY,JX));       Rain2SoilSurf_col          = 0._r8
  allocate(Irrig2SoilSurf(JY,JX));       Irrig2SoilSurf                = 0._r8
  allocate(LakeSurfFlowMicP_col(JY,JX));       LakeSurfFlowMicP_col    = 0._r8
  allocate(LakeSurfFlowMicPX(JY,JX));      LakeSurfFlowMicPX           = 0._r8
  allocate(LakeSurfFlowMacP_col(JY,JX));      LakeSurfFlowMacP_col     = 0._r8
  allocate(LakeSurfHeatFlux_col(JY,JX));      LakeSurfHeatFlux_col     = 0._r8
  end subroutine InitSurfSoilData

!----------------------------------------------------------------------
  subroutine DestructSurfSoilData
  use abortutils, only : destroy
  implicit none
  call destroy(FracSurfAsSnow)
  call destroy(FracSurfSnoFree)
  call destroy(LWRadBySurf_col)
  call destroy(HeatByRadiation_col)
  call destroy(HeatEvapAir2Surf_col)
  call destroy(HeatSensAir2Surf_col)
  call destroy(HeatSensVapAir2Surf_col)
  call destroy(HeatNet2Surf_col)
  call destroy(VapXAir2GSurf_col)
  call destroy(FracSurfBareSoil_col)
  call destroy(VWatStoreCapSurf)
  call destroy(MaxVLWatByLitR_col)
  call destroy(VHCPNX)
  call destroy(PARG)
  call destroy(Rain2SoilSurf_col)
  call destroy(Irrig2SoilSurf)
  call destroy(LakeSurfFlowMicP_col)
  call destroy(LakeSurfFlowMicPX)
  call destroy(LakeSurfFlowMacP_col)
  call destroy(LakeSurfHeatFlux_col)
  end subroutine DestructSurfSoilData

end module SurfSoilDataType
