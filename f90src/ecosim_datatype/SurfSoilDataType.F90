module SurfSoilDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__
  real(r8),target,allocatable ::  FracSurfAsSnow(:,:)                !fraction of snow cover
  real(r8),target,allocatable ::  FracSurfSnoFree(:,:)               !fraction of snow-free cover
  real(r8),target,allocatable ::  FracSurfAsBareSoi(:,:)             !fraction of exposed soil surface, [-]
  real(r8),target,allocatable ::  LWRadBySurf(:,:)                         !longwave radiation emitted from ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatByRadiation(:,:)                 !total net radiation at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatEvapAir2Surf(:,:)              !total latent heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatSensAir2Surf(:,:)              !total sensible heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatSensVapAir2Surf(:,:)           !total convective heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HeatNet2Surf(:,:)                  !total ground heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  VapXAir2GSurf(:,:)                 !negative of total evaporation at ground surface, [m3 d-2 t-1]
  real(r8),target,allocatable ::  VWatStoreCapSurf(:,:)                         !surface water storage capacity, [m3 d-2]
  real(r8),target,allocatable ::  MaxVLWatByLitR(:,:)                         !soil surface water retention capacity
  real(r8),target,allocatable ::  VHCPNX(:,:)                        !minimum heat capacities
  real(r8),target,allocatable ::  PARG(:,:,:)                        !soil surface boundary layer conductance, [m t-1]
  real(r8),target,allocatable ::  Rain2SoilSurf(:,:)                         !precipitation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  Irrig2SoilSurf(:,:)                         !irrifation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  LakeSurfFlowMicP(:,:)                         !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfFlowMicPX(:,:)                        !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfFlowMacP(:,:)                        !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfHeatFlux(:,:)              !lake surface heat flux, outgoing positive
!----------------------------------------------------------------------

contains
  subroutine InitSurfSoilData

  implicit none
  allocate(FracSurfAsSnow(JY,JX));         FracSurfAsSnow=0._r8
  allocate(FracSurfSnoFree(JY,JX));        FracSurfSnoFree=0._r8
  allocate(LWRadBySurf(JY,JX));             LWRadBySurf=0._r8
  allocate(HeatByRadiation(JY,JX));          HeatByRadiation=0._r8
  allocate(HeatEvapAir2Surf(JY,JX));       HeatEvapAir2Surf=0._r8
  allocate(HeatSensAir2Surf(JY,JX));       HeatSensAir2Surf=0._r8
  allocate(HeatSensVapAir2Surf(JY,JX));    HeatSensVapAir2Surf=0._r8
  allocate(HeatNet2Surf(JY,JX));           HeatNet2Surf=0._r8
  allocate(VapXAir2GSurf(JY,JX));          VapXAir2GSurf=0._r8
  allocate(FracSurfAsBareSoi(JY,JX));      FracSurfAsBareSoi=0._r8
  allocate(VWatStoreCapSurf(JY,JX));       VWatStoreCapSurf=0._r8
  allocate(MaxVLWatByLitR(JY,JX));       MaxVLWatByLitR=0._r8
  allocate(VHCPNX(JY,JX));      VHCPNX=0._r8
  allocate(PARG(60,JY,JX));     PARG=0._r8
  allocate(Rain2SoilSurf(JY,JX));       Rain2SoilSurf=0._r8
  allocate(Irrig2SoilSurf(JY,JX));       Irrig2SoilSurf=0._r8
  allocate(LakeSurfFlowMicP(JY,JX));       LakeSurfFlowMicP=0._r8
  allocate(LakeSurfFlowMicPX(JY,JX));      LakeSurfFlowMicPX=0._r8
  allocate(LakeSurfFlowMacP(JY,JX));      LakeSurfFlowMacP=0._r8
  allocate(LakeSurfHeatFlux(JY,JX));      LakeSurfHeatFlux=0._r8
  end subroutine InitSurfSoilData

!----------------------------------------------------------------------
  subroutine DestructSurfSoilData
  use abortutils, only : destroy
  implicit none
  call destroy(FracSurfAsSnow)
  call destroy(FracSurfSnoFree)
  call destroy(LWRadBySurf)
  call destroy(HeatByRadiation)
  call destroy(HeatEvapAir2Surf)
  call destroy(HeatSensAir2Surf)
  call destroy(HeatSensVapAir2Surf)
  call destroy(HeatNet2Surf)
  call destroy(VapXAir2GSurf)
  call destroy(FracSurfAsBareSoi)
  call destroy(VWatStoreCapSurf)
  call destroy(MaxVLWatByLitR)
  call destroy(VHCPNX)
  call destroy(PARG)
  call destroy(Rain2SoilSurf)
  call destroy(Irrig2SoilSurf)
  call destroy(LakeSurfFlowMicP)
  call destroy(LakeSurfFlowMicPX)
  call destroy(LakeSurfFlowMacP)
  call destroy(LakeSurfHeatFlux)
  end subroutine DestructSurfSoilData

end module SurfSoilDataType
