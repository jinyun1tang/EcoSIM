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
  real(r8),target,allocatable :: DOM_vr(:,:,:,:,:)                       !dissolved organic C micropore	[g d-2]
  real(r8),target,allocatable :: DOM_MacP_vr(:,:,:,:,:)                      !dissolved organic C macropore	[g d-2]
  real(r8),target,allocatable :: SoilOrgM_vr(:,:,:,:)                        !total soil organic C [g d-2]
  real(r8),target,allocatable :: ORGCX_vr(:,:,:)                       !SOC concentration	[g Mg-1]
  real(r8),target,allocatable :: SolidOMAct_vr(:,:,:,:,:)                     !colonized humus C in each complex [g d-2]
  real(r8),target,allocatable :: OMLitrC_vr(:,:,:)                        !total particulate organic C [g d-2]
  real(r8),target,allocatable :: CSoilOrgM_vr(:,:,:,:)                       !soil organic matter content [g kg-1]
  real(r8),target,allocatable :: COMLitrC_vr(:,:,:)                       !soil particulate C content [g kg-1]
  real(r8),target,allocatable :: ElmAllocmatMicrblitr2POM_vr(:,:,:,:)                     !allocation coefficient to humus fractions
  real(r8),target,allocatable :: tMicBiome_col(:,:,:)                         !total micriobial C, [g d-2]
  real(r8),target,allocatable :: tSoilOrgM_col(:,:,:)                  !total soil organic matter, include everything organic (exclude live roots) [g d-2]
  real(r8),target,allocatable :: tLitrOM_col(:,:,:)                        !total litter C, [g d-2]
  real(r8),target,allocatable :: litrOM_vr(:,:,:,:)
  real(r8),target,allocatable :: tHumOM_col(:,:,:)                        !total humus C, [g d-2]
  real(r8),target,allocatable :: EPOC(:,:,:)                       !partitioning coefficient between POC and litter, []
  real(r8),target,allocatable :: EHUM(:,:,:)                       !partitioning coefficient between humus and microbial residue, []
  real(r8),target,allocatable :: CDOM_vr(:,:,:,:,:)                     !DOC concentration, [g m-3]
  real(r8),target,allocatable :: FracBulkSOMC_vr(:,:,:,:)                    !fraction of total organic C in complex, [-]
  real(r8),target,allocatable :: DIC_mass_col(:,:)                        !total soil DIC, [g d-2]
  real(r8),target,allocatable :: tNH4_col(:,:)                         !total soil NH4 + NH3 content, [g d-2]
  real(r8),target,allocatable :: tNO3_col(:,:)                         !total soil NO3 + NO2 content, [g d-2]
  real(r8),target,allocatable :: tHxPO4_col(:,:)                         !total soil PO4 content, [g d-2]

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

  allocate(RSC(1:NumOfLitrCmplxs,0:JZ,JY,JX)); RSC=0._r8
  allocate(RSN(1:NumOfLitrCmplxs,0:JZ,JY,JX)); RSN=0._r8
  allocate(RSP(1:NumOfLitrCmplxs,0:JZ,JY,JX)); RSP=0._r8
  allocate(CFOSC(jsken,1:jcplx,0:JZ,JY,JX)); CFOSC=0._r8
  allocate(CNOSC(jsken,1:jcplx,0:JZ,JY,JX)); CNOSC=0._r8
  allocate(CPOSC(jsken,1:jcplx,0:JZ,JY,JX)); CPOSC=0._r8
  allocate(SolidOM_vr(NumPlantChemElms,jsken,1:jcplx,0:JZ,JY,JX));SolidOM_vr=0._r8
  allocate(SorbedOM_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));SorbedOM_vr=0._r8
  allocate(OMBioResdu_vr(1:NumPlantChemElms,ndbiomcp,1:jcplx,0:JZ,JY,JX)); OMBioResdu_vr=0._r8
  allocate(DOM_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));DOM_vr=0._r8
  allocate(DOM_MacP_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));DOM_MacP_vr=0._r8
  allocate(SoilOrgM_vr(1:NumPlantChemElms,0:JZ,JY,JX));SoilOrgM_vr=0._r8
  allocate(ORGCX_vr(0:JZ,JY,JX)); ORGCX_vr=0._r8
  allocate(SolidOMAct_vr(jsken,1:jcplx,0:JZ,JY,JX));SolidOMAct_vr=0._r8
  allocate(OMLitrC_vr(0:JZ,JY,JX));OMLitrC_vr=0._r8
  allocate(CSoilOrgM_vr(1:NumPlantChemElms,0:JZ,JY,JX));CSoilOrgM_vr=0._r8
  allocate(COMLitrC_vr(JZ,JY,JX));COMLitrC_vr=0._r8
  allocate(ElmAllocmatMicrblitr2POM_vr(2,JZ,JY,JX));; ElmAllocmatMicrblitr2POM_vr=0._r8
  allocate(tMicBiome_col(NumPlantChemElms,JY,JX));        tMicBiome_col=0._r8
  allocate(tSoilOrgM_col(NumPlantChemElms,JY,JX)); tSoilOrgM_col=0._r8
  allocate(tLitrOM_col(NumPlantChemElms,JY,JX));       tLitrOM_col=0._r8
  allocate(litrOM_vr(NumPlantChemElms,0:JZ,JY,JX)); litrOM_vr=0._r8
  allocate(tHumOM_col(NumPlantChemElms,JY,JX));       tHumOM_col=0._r8
  allocate(EPOC(0:JZ,JY,JX));   EPOC=0._r8
  allocate(EHUM(0:JZ,JY,JX));   EHUM=0._r8
  allocate(CDOM_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));CDOM_vr=0._r8
  allocate(FracBulkSOMC_vr(1:jcplx,0:JZ,JY,JX));FracBulkSOMC_vr=0._r8
  allocate(DIC_mass_col(JY,JX));       DIC_mass_col=0._r8
  allocate(tNH4_col(JY,JX));        tNH4_col=0._r8
  allocate(tNO3_col(JY,JX));        tNO3_col=0._r8
  allocate(tHxPO4_col(JY,JX));        tHxPO4_col=0._r8

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------

  subroutine DestructSOMData

  use abortutils, only : destroy
  implicit none

  call destroy(litrOM_vr)
  call destroy(RSC)
  call destroy(RSN)
  call destroy(RSP)
  call destroy(CFOSC)
  call destroy(CNOSC)
  call destroy(CPOSC)
  call destroy(SolidOM_vr)
  call destroy(SorbedOM_vr)
  call destroy(OMBioResdu_vr)
  call destroy(DOM_vr)
  call destroy(DOM_MacP_vr)
  call destroy(SoilOrgM_vr)
  call destroy(ORGCX_vr)
  call destroy(SolidOMAct_vr)
  call destroy(OMLitrC_vr)
  call destroy(CSoilOrgM_vr)
  call destroy(COMLitrC_vr)
  call destroy(ElmAllocmatMicrblitr2POM_vr)
  call destroy(tMicBiome_col)
  call destroy(tSoilOrgM_col)
  call destroy(tLitrOM_col)
  call destroy(tHumOM_col)
  call destroy(EPOC)
  call destroy(EHUM)

  call destroy(CDOM_vr)
  call destroy(FracBulkSOMC_vr)
  call destroy(DIC_mass_col)
  call destroy(tNH4_col)
  call destroy(tNO3_col)
  call destroy(tHxPO4_col)
  end subroutine DestructSOMData
end module SOMDataType
