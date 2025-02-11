module SurfLitterDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx=> jcplxc, jcplx1=>jcplxcm1
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) ,target,allocatable ::   BulkDensLitR(:)                   !surface litter bulk density	[Mg m-3]
  real(r8) ,target,allocatable ::   PARR_col(:,:)                     !surface litter boundary layer conductance, [m t-1]
  integer  ,target,allocatable ::   iLitrType_col(:,:,:)              !surface litter type:1 = plant, 2 = manure
  real(r8) ,target,allocatable ::   XTillCorp_col(:,:)                !factor for surface litter incorporation and soil mixing
  real(r8) ,target,allocatable ::   WatFLo2LitrM(:,:,:)               !water transfer between soil surface and surface litter, [g d-2 t-1]
  real(r8) ,target,allocatable ::   WatFlowSno2LitRM_col(:,:,:)       !meltwater flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   FracSurfByLitR_col(:,:)           !fraction of soil surface covered by surface litter, [-]
  real(r8) ,target,allocatable ::   HeatFLoByWat2LitR_col(:,:)        !net heat transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,target,allocatable ::   VLitR_col(:,:)                    !surface litter volume, [m3 d-2]
  real(r8) ,target,allocatable ::   VHeatCapLitRMin_col(:,:)          !threshold surface litter heat capacity, [MJ d-2 K-1]
  real(r8) ,target,allocatable ::   VWatLitRHoldCapcity_col(:,:)      !surface litter water holding capacity, [m3 d-2]
  real(r8) ,target,allocatable ::   WatFLo2LitR_col(:,:)              !net water transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,target,allocatable ::   TLitrIceFlxThaw_col(:,:)          !water from ice thaw in surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   TLitrIceHeatFlxFrez_col(:,:)      !latent heat released from water freeze in surface litter, [MJ d-2 h-1]
  real(r8) ,target,allocatable ::   Rain2LitRSurf_col(:,:)            !precipitation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   Irrig2LitRSurf_col(:,:)               !irrigation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   POROS0_col(:,:)                   !litter porosity, [m3 pore m-3 litr]
  real(r8) ,target,allocatable ::   RC0(:,:,:)                        !surface litter in each complex,	[g d-2]
  real(r8) ,target,allocatable ::   RC0ff(:,:)
  real(r8) ,target,allocatable ::   LitWatMassBeg_col(:,:)            !total inital water mass in litter layer, [m3 H2O d-2]
  real(r8) ,target,allocatable ::   LitWatMassEnd_col(:,:)            !total inital water mass in litter layer, [m3 H2O d-2]
  real(r8) ,target,allocatable ::   Rain2LitR_col(:,:)                !total precipiation reaches the litter layer [m3 H3O d-2 h-1]
  private :: InitAllocate
  contains
!----------------------------------------------------------------------

  subroutine InitSurfLitter(NumOfLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  call InitAllocate(NumOfLitrCmplxs)
  end subroutine InitSurfLitter
!----------------------------------------------------------------------

  subroutine InitAllocate(NumOfLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  allocate(Rain2LitR_col(JY,JX)); Rain2LitR_col=0._r8
  allocate(LitWatMassBeg_col(JY,JX)); LitWatMassBeg_col = 0._r8
  allocate(LitWatMassEnd_col(JY,JX)); LitWatMassEnd_col = 0._r8
  allocate(BulkDensLitR(1:NumOfLitrCmplxs));    BulkDensLitR =0._r8
  allocate(PARR_col(JY,JX));         PARR_col=0._r8
  allocate(iLitrType_col(2,JY,JX));      iLitrType_col=0
  allocate(XTillCorp_col(JY,JX));        XTillCorp_col=0._r8
  allocate(WatFLo2LitrM(60,JY,JX));     WatFLo2LitrM=0._r8
  allocate(WatFlowSno2LitRM_col(60,JY,JX));     WatFlowSno2LitRM_col=0._r8
  allocate(FracSurfByLitR_col(JY,JX));         FracSurfByLitR_col=0._r8
  allocate(HeatFLoByWat2LitR_col(JY,JX));        HeatFLoByWat2LitR_col=0._r8
  allocate(VLitR_col(JY,JX));         VLitR_col=0._r8
  allocate(VHeatCapLitRMin_col(JY,JX));       VHeatCapLitRMin_col=0._r8
  allocate(VWatLitRHoldCapcity_col(JY,JX));       VWatLitRHoldCapcity_col=0._r8
  allocate(WatFLo2LitR_col(JY,JX));         WatFLo2LitR_col=0._r8
  allocate(TLitrIceFlxThaw_col(JY,JX));        TLitrIceFlxThaw_col=0._r8
  allocate(TLitrIceHeatFlxFrez_col(JY,JX));       TLitrIceHeatFlxFrez_col=0._r8
  allocate(Rain2LitRSurf_col(JY,JX));        Rain2LitRSurf_col=0._r8
  allocate(Irrig2LitRSurf_col(JY,JX));        Irrig2LitRSurf_col=0._r8
  allocate(POROS0_col(JY,JX));       POROS0_col=0._r8
  allocate(RC0(1:NumOfLitrCmplxs,JY,JX));      Rc0=0._r8
  allocate(RC0ff(JY,JX)); RC0ff=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSurfLitter
  use abortutils, only : destroy
  implicit none

  call destroy(Rain2LitR_col)
  call destroy(LitWatMassBeg_col)
  call destroy(LitWatMassEnd_col)
  call destroy(BulkDensLitR)
  call destroy(PARR_col)
  call destroy(iLitrType_col)
  call destroy(XTillCorp_col)
  call destroy(WatFLo2LitrM)
  call destroy(WatFlowSno2LitRM_col)
  call destroy(FracSurfByLitR_col)
  call destroy(HeatFLoByWat2LitR_col)
  call destroy(VLitR_col)
  call destroy(VHeatCapLitRMin_col)
  call destroy(VWatLitRHoldCapcity_col)
  call destroy(WatFLo2LitR_col)
  call destroy(TLitrIceFlxThaw_col)
  call destroy(TLitrIceHeatFlxFrez_col)
  call destroy(Rain2LitRSurf_col)
  call destroy(Irrig2LitRSurf_col)
  call destroy(POROS0_col)
  call destroy(RC0)
  call destroy(RC0ff)

  end subroutine DestructSurfLitter

end module SurfLitterDataType
