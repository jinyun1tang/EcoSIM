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
  real(r8),target,allocatable ::  Eco_GPP_col(:,:)                          !ecosystem GPP, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_AutoR_col(:,:)                          !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_NPP_col(:,:)                          !ecosystem NPP, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_HR_col(:,:)                          !ecosystem heterotrophic respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  EcoHavstElmnt_col(:,:,:)                      !ecosystem harvest , [g d-2]
  real(r8),target,allocatable ::  NetNH4Mineralize_col(:,:)                        !total NH4 net mineraln (-ve) or immobiln (+ve)
  real(r8),target,allocatable ::  NetPO4Mineralize_col(:,:)                        !total H2PO4 net mineraln (-ve) or immobiln (+ve)
  real(r8),target,allocatable ::  GPP(:,:)                           !gross primary productivity, [g d-2 h-1]
  real(r8),target,allocatable ::  Canopy_NEE_col(:,:)                         !total net CO2 fixation
  real(r8),target,allocatable ::  LitterFallChemElmnt_col(:,:,:)                       !total litterfall element, [g d-2 h-1]
  real(r8),target,allocatable ::  ECO_ER_col(:,:)                          !ecosystem respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  Eco_NBP_col(:,:)                          !total NBP, [g d-2]
  real(r8),target,allocatable ::  RP14X(:,:,:)                       !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RP14Y(:,:,:)                       !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RP1BX(:,:,:)                       !HPO4 demand in band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RP1BY(:,:,:)                       !HPO4 demand in band by all microbial,root,myco populations
!----------------------------------------------------------------------

contains
  subroutine InitEcosimBGCFluxData

  implicit none
  allocate(Eco_NetRad_col(JY,JX));         Eco_NetRad_col=0._r8
  allocate(Eco_Heat_Latent_col(JY,JX));         Eco_Heat_Latent_col=0._r8
  allocate(Eco_Heat_Sens_col(JY,JX));         Eco_Heat_Sens_col=0._r8
  allocate(Eco_Heat_Grnd_col(JY,JX));         Eco_Heat_Grnd_col=0._r8
  allocate(Eco_GPP_col(JY,JX));        Eco_GPP_col=0._r8
  allocate(Eco_AutoR_col(JY,JX));        Eco_AutoR_col=0._r8
  allocate(Eco_NPP_col(JY,JX));        Eco_NPP_col=0._r8
  allocate(Eco_HR_col(JY,JX));        Eco_HR_col=0._r8
  allocate(EcoHavstElmnt_col(NumPlantChemElms,JY,JX));      EcoHavstElmnt_col=0._r8
  allocate(NetNH4Mineralize_col(JY,JX));      NetNH4Mineralize_col=0._r8
  allocate(NetPO4Mineralize_col(JY,JX));      NetPO4Mineralize_col=0._r8
  allocate(GPP(JY,JX));         GPP=0._r8
  allocate(Canopy_NEE_col(JY,JX));       Canopy_NEE_col=0._r8
  allocate(LitterFallChemElmnt_col(NumPlantChemElms,JY,JX));       LitterFallChemElmnt_col=0._r8
  allocate(ECO_ER_col(JY,JX));        ECO_ER_col=0._r8
  allocate(Eco_NBP_col(JY,JX));        Eco_NBP_col=0._r8
  allocate(RP14X(0:JZ,JY,JX));  RP14X=0._r8
  allocate(RP14Y(0:JZ,JY,JX));  RP14Y=0._r8
  allocate(RP1BX(0:JZ,JY,JX));  RP1BX=0._r8
  allocate(RP1BY(0:JZ,JY,JX));  RP1BY=0._r8
  end subroutine InitEcosimBGCFluxData

!----------------------------------------------------------------------
  subroutine DestructEcosimBGCFluxData
  use abortutils, only : destroy
  implicit none
  call destroy(Eco_NetRad_col)
  call destroy(Eco_Heat_Latent_col)
  call destroy(Eco_Heat_Sens_col)
  call destroy(Eco_Heat_Grnd_col)
  call destroy(Eco_GPP_col)
  call destroy(Eco_AutoR_col)
  call destroy(Eco_NPP_col)
  call destroy(Eco_HR_col)
  call destroy(EcoHavstElmnt_col)
  call destroy(NetNH4Mineralize_col)
  call destroy(NetPO4Mineralize_col)
  call destroy(GPP)
  call destroy(Canopy_NEE_col)
  call destroy(LitterFallChemElmnt_col)
  call destroy(ECO_ER_col)
  call destroy(Eco_NBP_col)
  call destroy(RP14X)
  call destroy(RP14Y)
  call destroy(RP1BX)
  call destroy(RP1BY)
  end subroutine DestructEcosimBGCFluxData

end module EcosimBGCFluxType
