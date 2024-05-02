module SOMDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc, ndbiomcp=>NumDeadMicrbCompts
  use data_const_mod, only : spval => DAT_CONST_SPVAL  
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: RSC(:,:,:,:)                       !initial surface litter [C	g m-2]
  real(r8),target,allocatable :: RSN(:,:,:,:)                       !initial surface litter [N	g m-2]
  real(r8),target,allocatable :: RSP(:,:,:,:)                       !initial surface litter [P	g m-2]
  real(r8),target,allocatable :: CFOSC(:,:,:,:,:)                   !fraction of SOC in kinetic components
  real(r8),target,allocatable :: CNOSC(:,:,:,:,:)                   !N:C,ratios of SOC kinetic components
  real(r8),target,allocatable :: CPOSC(:,:,:,:,:)                   !P:C ratios of SOC kinetic components
  real(r8),target,allocatable :: SolidOM_vr(:,:,:,:,:,:)                   !humus soil OM	[g d-2]
  real(r8),target,allocatable :: SorbedOM_vr(:,:,:,:,:)                     !adsorbed soil C	[g d-2]
  real(r8),target,allocatable :: OMBioResdu_vr(:,:,:,:,:,:)                     !microbial residue [C	g d-2]
  real(r8),target,allocatable :: DOM(:,:,:,:,:)                       !dissolved organic C micropore	[g d-2]
  real(r8),target,allocatable :: DOM_Macp(:,:,:,:,:)                      !dissolved organic C macropore	[g d-2]
  real(r8),target,allocatable :: ORGC_vr(:,:,:)                        !total soil organic C [g d-2]
  real(r8),target,allocatable :: ORGN_vr(:,:,:)                        !total soil organic N [g d-2]
  real(r8),target,allocatable :: ORGP_vr(:,:,:)                        !total soil organic P [g d-2]
  real(r8),target,allocatable :: ORGCX(:,:,:)                       !SOC concentration	[g Mg-1]
  real(r8),target,allocatable :: SolidOMAct_vr(:,:,:,:,:)                     !colonized humus C in each complex [g d-2]
  real(r8),target,allocatable :: OMLitrC_vr(:,:,:)                        !total particulate organic C [g d-2]
  real(r8),target,allocatable :: CORGC(:,:,:)                       !soil organic C content [g kg-1]
  real(r8),target,allocatable :: CORGN(:,:,:)                       !soil organic N content [mg kg-1]
  real(r8),target,allocatable :: CORGP(:,:,:)                       !soil organic P content  [mg kg-1]
  real(r8),target,allocatable :: COMLitrC_vr(:,:,:)                       !soil particulate C content [g kg-1]
  real(r8),target,allocatable :: CFOMC(:,:,:,:)                     !allocation coefficient to humus fractions
  real(r8),target,allocatable ::  TOMET(:,:,:)                         !total micriobial C, [g d-2]
  real(r8),target,allocatable ::  URSDM(:,:,:)                        !total litter C, [g d-2]
  real(r8),target,allocatable ::  UORGM(:,:,:)                        !total humus C, [g d-2]
  real(r8),target,allocatable ::  EPOC(:,:,:)                       !partitioning coefficient between POC and litter, []
  real(r8),target,allocatable ::  EHUM(:,:,:)                       !partitioning coefficient between humus and microbial residue, []
  real(r8),target,allocatable ::  CDOM(:,:,:,:,:)                     !DOC concentration, [g m-3]
  real(r8),target,allocatable ::  FracBulkSOMC_vr(:,:,:,:)                    !fraction of total organic C in complex, [-]
  real(r8),target,allocatable ::  DIC_mass_col(:,:)                        !total soil DIC, [g d-2]
  real(r8),target,allocatable ::  UNH4(:,:)                         !total soil NH4 + NH3 content, [g d-2]
  real(r8),target,allocatable ::  UNO3(:,:)                         !total soil NO3 + NO2 content, [g d-2]
  real(r8),target,allocatable ::  UPO4(:,:)                         !total soil PO4 content, [g d-2]
  real(r8),target,allocatable ::  OMCL(:,:,:)
  real(r8),target,allocatable ::  OMNL(:,:,:)

  private :: InitAllocate
  contains

  subroutine InitSOMData(NumOfLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  call InitAllocate(NumOfLitrCmplxs)

  end subroutine InitSOMData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate(NumOfLitrCmplxs)
  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  allocate(RSC(1:NumOfLitrCmplxs,0:JZ,JY,JX))
  allocate(RSN(1:NumOfLitrCmplxs,0:JZ,JY,JX))
  allocate(RSP(1:NumOfLitrCmplxs,0:JZ,JY,JX))
  allocate(CFOSC(jsken,1:jcplx,0:JZ,JY,JX))
  allocate(CNOSC(jsken,1:jcplx,0:JZ,JY,JX))
  allocate(CPOSC(jsken,1:jcplx,0:JZ,JY,JX))
  allocate(SolidOM_vr(NumPlantChemElms,jsken,1:jcplx,0:JZ,JY,JX));SolidOM_vr=spval
  allocate(SorbedOM_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));SorbedOM_vr=spval
  allocate(OMBioResdu_vr(1:NumPlantChemElms,ndbiomcp,1:jcplx,0:JZ,JY,JX)); OMBioResdu_vr=spval
  allocate(DOM(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX))
  allocate(DOM_Macp(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));DOM_MacP=0._r8
  allocate(ORGC_vr(0:JZ,JY,JX));ORGC_vr=0._r8
  allocate(ORGN_vr(0:JZ,JY,JX));ORGN_vr=0._r8
  allocate(ORGP_vr(0:JZ,JY,JX));ORGP_vr=0._r8
  allocate(ORGCX(0:JZ,JY,JX))
  allocate(SolidOMAct_vr(jsken,1:jcplx,0:JZ,JY,JX));SolidOMAct_vr=0._r8
  allocate(OMLitrC_vr(0:JZ,JY,JX));OMLitrC_vr=0._r8
  allocate(CORGC(0:JZ,JY,JX))
  allocate(CORGN(JZ,JY,JX))
  allocate(CORGP(JZ,JY,JX))
  allocate(COMLitrC_vr(JZ,JY,JX));COMLitrC_vr=0._r8
  allocate(CFOMC(2,JZ,JY,JX))
  allocate(TOMET(NumPlantChemElms,JY,JX));        TOMET=0._r8
  allocate(URSDM(NumPlantChemElms,JY,JX));       URSDM=0._r8
  allocate(UORGM(NumPlantChemElms,JY,JX));       UORGM=0._r8
  allocate(EPOC(0:JZ,JY,JX));   EPOC=0._r8
  allocate(EHUM(0:JZ,JY,JX));   EHUM=0._r8
  allocate(CDOM(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));CDOM=0._r8
  allocate(FracBulkSOMC_vr(1:jcplx,0:JZ,JY,JX));FracBulkSOMC_vr=0._r8
  allocate(DIC_mass_col(JY,JX));       DIC_mass_col=0._r8
  allocate(UNH4(JY,JX));        UNH4=0._r8
  allocate(UNO3(JY,JX));        UNO3=0._r8
  allocate(UPO4(JY,JX));        UPO4=0._r8
  ALLOCATE(OMCL(0:JZ,JY,JX));   OMCL=0._r8
  ALLOCATE(OMNL(0:JZ,JY,JX));   OMNL=0._r8
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------

  subroutine DestructSOMData

  use abortutils, only : destroy
  implicit none

  call destroy(RSC)
  call destroy(RSN)
  call destroy(RSP)
  call destroy(CFOSC)
  call destroy(CNOSC)
  call destroy(CPOSC)
  call destroy(SolidOM_vr)
  call destroy(SorbedOM_vr)
  call destroy(OMBioResdu_vr)
  call destroy(DOM)
  call destroy(DOM_MacP)
  call destroy(ORGC_vr)
  call destroy(ORGN_vr)
  call destroy(ORGP_vr)
  call destroy(ORGCX)
  call destroy(SolidOMAct_vr)
  call destroy(OMLitrC_vr)
  call destroy(CORGC)
  call destroy(CORGN)
  call destroy(CORGP)
  call destroy(COMLitrC_vr)
  call destroy(CFOMC)
  call destroy(TOMET)
  call destroy(URSDM)
  call destroy(UORGM)
  call destroy(EPOC)
  call destroy(EHUM)

  call destroy(CDOM)
  call destroy(FracBulkSOMC_vr)
  call destroy(DIC_mass_col)
  call destroy(UNH4)
  call destroy(UNO3)
  call destroy(UPO4)
  call destroy(OMCL)
  call destroy(OMNL)
  end subroutine DestructSOMData
end module SOMDataType
