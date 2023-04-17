module EcosimBGCFluxType

!
! Ecosystm fluxes for C, N, and P budget
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target,allocatable ::  TCAN(:,:)                          !total net CO2 fixation
  real(r8),target,allocatable ::  TRN(:,:)                           !ecosystem net radiation, [MJ d-2 h-1]
  real(r8),target,allocatable ::  TLE(:,:)                           !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  TSH(:,:)                           !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  TGH(:,:)                           !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  TGPP(:,:)                          !ecosystem GPP, [g d-2 h-1]
  real(r8),target,allocatable ::  TRAU(:,:)                          !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  TNPP(:,:)                          !ecosystem NPP, [g d-2 h-1]
  real(r8),target,allocatable ::  THRE(:,:)                          !ecosystem heterotrophic respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  XHVSTE(:,:,:)                      !ecosystem harvest , [g d-2]
  real(r8),target,allocatable ::  TRINH4(:,:)                        !total NH4 net mineraln (-ve) or immobiln (+ve)
  real(r8),target,allocatable ::  TRIPO4(:,:)                        !total H2PO4 net mineraln (-ve) or immobiln (+ve)
  real(r8),target,allocatable ::  GPP(:,:)                           !gross primary productivity, [g d-2 h-1]
  real(r8),target,allocatable ::  TCCAN(:,:)                         !total net CO2 fixation
  real(r8),target,allocatable ::  ZESNC(:,:,:)                       !total litterfall element, [g d-2 h-1]
  real(r8),target,allocatable ::  RECO(:,:)                          !ecosystem respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  TNBP(:,:)                          !total NBP, [g d-2]
  real(r8),target,allocatable ::  RP14X(:,:,:)                       !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RP14Y(:,:,:)                       !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RP1BX(:,:,:)                       !HPO4 demand in band by all microbial,root,myco populations
  real(r8),target,allocatable ::  RP1BY(:,:,:)                       !HPO4 demand in band by all microbial,root,myco populations
!----------------------------------------------------------------------

contains
  subroutine InitEcosimBGCFluxData

  implicit none
  allocate(TCAN(JY,JX));        TCAN=0._r8
  allocate(TRN(JY,JX));         TRN=0._r8
  allocate(TLE(JY,JX));         TLE=0._r8
  allocate(TSH(JY,JX));         TSH=0._r8
  allocate(TGH(JY,JX));         TGH=0._r8
  allocate(TGPP(JY,JX));        TGPP=0._r8
  allocate(TRAU(JY,JX));        TRAU=0._r8
  allocate(TNPP(JY,JX));        TNPP=0._r8
  allocate(THRE(JY,JX));        THRE=0._r8
  allocate(XHVSTE(npelms,JY,JX));      XHVSTE=0._r8
  allocate(TRINH4(JY,JX));      TRINH4=0._r8
  allocate(TRIPO4(JY,JX));      TRIPO4=0._r8
  allocate(GPP(JY,JX));         GPP=0._r8
  allocate(TCCAN(JY,JX));       TCCAN=0._r8
  allocate(ZESNC(npelms,JY,JX));       ZESNC=0._r8
  allocate(RECO(JY,JX));        RECO=0._r8
  allocate(TNBP(JY,JX));        TNBP=0._r8
  allocate(RP14X(0:JZ,JY,JX));  RP14X=0._r8
  allocate(RP14Y(0:JZ,JY,JX));  RP14Y=0._r8
  allocate(RP1BX(0:JZ,JY,JX));  RP1BX=0._r8
  allocate(RP1BY(0:JZ,JY,JX));  RP1BY=0._r8
  end subroutine InitEcosimBGCFluxData

!----------------------------------------------------------------------
  subroutine DestructEcosimBGCFluxData
  use abortutils, only : destroy
  call destroy(TCAN)
  call destroy(TRN)
  call destroy(TLE)
  call destroy(TSH)
  call destroy(TGH)
  call destroy(TGPP)
  call destroy(TRAU)
  call destroy(TNPP)
  call destroy(THRE)
  call destroy(XHVSTE)
  call destroy(TRINH4)
  call destroy(TRIPO4)
  call destroy(GPP)
  call destroy(TCCAN)
  call destroy(ZESNC)
  call destroy(RECO)
  call destroy(TNBP)
  call destroy(RP14X)
  call destroy(RP14Y)
  call destroy(RP1BX)
  call destroy(RP1BY)
  end subroutine DestructEcosimBGCFluxData

end module EcosimBGCFluxType
