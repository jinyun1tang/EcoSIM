module EcosimBGCFluxType

!
! Ecosystm fluxes for C, N, and P budget
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  Eco_NetRad_col(:,:)                           !ecosystem net radiation, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_Heat_Latent_col(:,:)                           !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_Heat_Sens_col(:,:)                           !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_Heat_Grnd_col(:,:)                           !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_GPP_CumYr_col(:,:)                          !ecosystem GPP, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_AutoR_CumYr_col(:,:)                          !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_NPP_CumYr_col(:,:)                          !ecosystem NPP, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_HR_CumYr_col(:,:)                          !ecosystem heterotrophic respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  EcoHavstElmnt_CumYr_col(:,:,:)                      !ecosystem harvest , [g d-2]
  real(r8),target,allocatable ::  NetNH4Mineralize_CumYr_col(:,:)                        !total NH4 net mineraln (-ve) or immobiln (+ve)
  real(r8),target,allocatable ::  NetPO4Mineralize_CumYr_col(:,:)                        !total H2PO4 net mineraln (-ve) or immobiln (+ve)
  real(r8),target,allocatable ::  GPP(:,:)                           !gross primary productivity, [g d-2 h-1]
  real(r8),target,allocatable ::  Canopy_NEE_col(:,:)                         !total net CO2 fixation
  real(r8),target,allocatable ::  LitrFallStrutElms_col(:,:,:)                       !total LitrFall element, [g d-2 h-1]
  real(r8),target,allocatable ::  ECO_ER_col(:,:)                          !ecosystem respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_NBP_CumYr_col(:,:)                          !total NBP, [g d-2]
  real(r8),target,allocatable ::  REcoH1PO4DmndSoil_vr(:,:,:)                       !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RH1PO4EcoDmndSoilPrev_vr(:,:,:)                       !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8),target,allocatable ::  REcoH1PO4DmndBand_vr(:,:,:)                       !HPO4 demand in band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RH1PO4EcoDmndBandPrev_vr(:,:,:)                       !HPO4 demand in band by all microbial,root,myco populations
!----------------------------------------------------------------------

contains
  subroutine InitEcosimBGCFluxData

  implicit none
  allocate(Eco_NetRad_col(JY,JX));         Eco_NetRad_col=0._r8
  allocate(Eco_Heat_Latent_col(JY,JX));         Eco_Heat_Latent_col=0._r8
  allocate(Eco_Heat_Sens_col(JY,JX));         Eco_Heat_Sens_col=0._r8
  allocate(Eco_Heat_Grnd_col(JY,JX));         Eco_Heat_Grnd_col=0._r8
  allocate(Eco_GPP_CumYr_col(JY,JX));        Eco_GPP_CumYr_col=0._r8
  allocate(Eco_AutoR_CumYr_col(JY,JX));        Eco_AutoR_CumYr_col=0._r8
  allocate(Eco_NPP_CumYr_col(JY,JX));        Eco_NPP_CumYr_col=0._r8
  allocate(Eco_HR_CumYr_col(JY,JX));        Eco_HR_CumYr_col=0._r8
  allocate(EcoHavstElmnt_CumYr_col(NumPlantChemElms,JY,JX));      EcoHavstElmnt_CumYr_col=0._r8
  allocate(NetNH4Mineralize_CumYr_col(JY,JX));      NetNH4Mineralize_CumYr_col=0._r8
  allocate(NetPO4Mineralize_CumYr_col(JY,JX));      NetPO4Mineralize_CumYr_col=0._r8
  allocate(GPP(JY,JX));         GPP=0._r8
  allocate(Canopy_NEE_col(JY,JX));       Canopy_NEE_col=0._r8
  allocate(LitrFallStrutElms_col(NumPlantChemElms,JY,JX));       LitrFallStrutElms_col=0._r8
  allocate(ECO_ER_col(JY,JX));        ECO_ER_col=0._r8
  allocate(Eco_NBP_CumYr_col(JY,JX));        Eco_NBP_CumYr_col=0._r8
  allocate(REcoH1PO4DmndSoil_vr(0:JZ,JY,JX));  REcoH1PO4DmndSoil_vr=0._r8
  allocate(RH1PO4EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RH1PO4EcoDmndSoilPrev_vr=0._r8
  allocate(REcoH1PO4DmndBand_vr(0:JZ,JY,JX));  REcoH1PO4DmndBand_vr=0._r8
  allocate(RH1PO4EcoDmndBandPrev_vr(0:JZ,JY,JX));  RH1PO4EcoDmndBandPrev_vr=0._r8
  end subroutine InitEcosimBGCFluxData

!----------------------------------------------------------------------
  subroutine DestructEcosimBGCFluxData
  use abortutils, only : destroy
  implicit none
  call destroy(Eco_NetRad_col)
  call destroy(Eco_Heat_Latent_col)
  call destroy(Eco_Heat_Sens_col)
  call destroy(Eco_Heat_Grnd_col)
  call destroy(Eco_GPP_CumYr_col)
  call destroy(Eco_AutoR_CumYr_col)
  call destroy(Eco_NPP_CumYr_col)
  call destroy(Eco_HR_CumYr_col)
  call destroy(EcoHavstElmnt_CumYr_col)
  call destroy(NetNH4Mineralize_CumYr_col)
  call destroy(NetPO4Mineralize_CumYr_col)
  call destroy(GPP)
  call destroy(Canopy_NEE_col)
  call destroy(LitrFallStrutElms_col)
  call destroy(ECO_ER_col)
  call destroy(Eco_NBP_CumYr_col)
  call destroy(REcoH1PO4DmndSoil_vr)
  call destroy(RH1PO4EcoDmndSoilPrev_vr)
  call destroy(REcoH1PO4DmndBand_vr)
  call destroy(RH1PO4EcoDmndBandPrev_vr)
  end subroutine DestructEcosimBGCFluxData

end module EcosimBGCFluxType
