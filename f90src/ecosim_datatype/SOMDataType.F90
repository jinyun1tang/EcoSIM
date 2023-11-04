module SOMDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc, ndbiomcp=>NumOfDeadMicrobiomComponents
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
  real(r8),target,allocatable :: OSC(:,:,:,:,:)                     !humus soil C	[g d-2]
  real(r8),target,allocatable :: OSN(:,:,:,:,:)                     !humus soil N	[g d-2]
  real(r8),target,allocatable :: OSP(:,:,:,:,:)                     !humus soil P	[g d-2]
  real(r8),target,allocatable :: OHC(:,:,:,:)                       !adsorbed soil C	[g d-2]
  real(r8),target,allocatable :: OHN(:,:,:,:)                       !adsorbed soil N	[g d-2]
  real(r8),target,allocatable :: OHP(:,:,:,:)                       !adsorbed soil P	[g d-2]
  real(r8),target,allocatable :: OHA(:,:,:,:)                       !adsorbed soil acetate	[g d-2]
  real(r8),target,allocatable :: ORC(:,:,:,:,:)                     !microbial residue [C	g d-2]
  real(r8),target,allocatable :: ORN(:,:,:,:,:)                     !microbial residue N	[g d-2]
  real(r8),target,allocatable :: ORP(:,:,:,:,:)                     !microbial residue P	[g d-2]
  real(r8),target,allocatable :: DOM(:,:,:,:,:)                       !dissolved organic C micropore	[g d-2]
  real(r8),target,allocatable :: DOM_Macp(:,:,:,:,:)                      !dissolved organic C macropore	[g d-2]
  real(r8),target,allocatable :: ORGC(:,:,:)                        !total soil organic C [g d-2]
  real(r8),target,allocatable :: ORGN(:,:,:)                        !total soil organic N [g d-2]
  real(r8),target,allocatable :: ORGCX(:,:,:)                       !SOC concentration	[g Mg-1]
  real(r8),target,allocatable :: OSA(:,:,:,:,:)                     !colonized humus C in each complex [g d-2]
  real(r8),target,allocatable :: ORGR(:,:,:)                        !total particulate organic C [g d-2]
  real(r8),target,allocatable :: CORGC(:,:,:)                       !soil organic C content [g kg-1]
  real(r8),target,allocatable :: CORGN(:,:,:)                       !soil organic N content [mg kg-1]
  real(r8),target,allocatable :: CORGP(:,:,:)                       !soil organic P content  [mg kg-1]
  real(r8),target,allocatable :: CORGR(:,:,:)                       !soil particulate C content [g kg-1]
  real(r8),target,allocatable :: CFOMC(:,:,:,:)                     !allocation coefficient to humus fractions
  real(r8),target,allocatable ::  TOMT(:,:)                         !total micriobial C, [g d-2]
  real(r8),target,allocatable ::  TONT(:,:)                         !total micriobial N, [g d-2]
  real(r8),target,allocatable ::  TOPT(:,:)                         !total micriobial P, [g d-2]
  real(r8),target,allocatable ::  URSDC(:,:)                        !total litter C, [g d-2]
  real(r8),target,allocatable ::  UORGC(:,:)                        !total humus C, [g d-2]
  real(r8),target,allocatable ::  URSDN(:,:)                        !total litter N, [g d-2]
  real(r8),target,allocatable ::  URSDP(:,:)                        !total litter P, [g d-2]
  real(r8),target,allocatable ::  UORGN(:,:)                        !total humus N, [g d-2]
  real(r8),target,allocatable ::  UORGP(:,:)                        !total humus P, [g d-2]
  real(r8),target,allocatable ::  EPOC(:,:,:)                       !partitioning coefficient between POC and litter, []
  real(r8),target,allocatable ::  EHUM(:,:,:)                       !partitioning coefficient between humus and microbial residue, []
  real(r8),target,allocatable ::  CDOM(:,:,:,:,:)                     !DOC concentration, [g m-3]
  real(r8),target,allocatable ::  FOSRH(:,:,:,:)                    !fraction of total organic C in complex, [-]
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
  allocate(OSC(jsken,1:jcplx,0:JZ,JY,JX))
  allocate(OSN(jsken,1:jcplx,0:JZ,JY,JX))
  allocate(OSP(jsken,1:jcplx,0:JZ,JY,JX))
  allocate(OHC(1:jcplx,0:JZ,JY,JX))
  allocate(OHN(1:jcplx,0:JZ,JY,JX))
  allocate(OHP(1:jcplx,0:JZ,JY,JX))
  allocate(OHA(1:jcplx,0:JZ,JY,JX))
  allocate(ORC(ndbiomcp,1:jcplx,0:JZ,JY,JX))
  allocate(ORN(ndbiomcp,1:jcplx,0:JZ,JY,JX))
  allocate(ORP(ndbiomcp,1:jcplx,0:JZ,JY,JX))
  allocate(DOM(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX))
  allocate(DOM_Macp(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));DOM_MacP=0._r8
  allocate(ORGC(0:JZ,JY,JX))
  allocate(ORGN(0:JZ,JY,JX))
  allocate(ORGCX(0:JZ,JY,JX))
  allocate(OSA(jsken,1:jcplx,0:JZ,JY,JX))
  allocate(ORGR(0:JZ,JY,JX))
  allocate(CORGC(0:JZ,JY,JX))
  allocate(CORGN(JZ,JY,JX))
  allocate(CORGP(JZ,JY,JX))
  allocate(CORGR(JZ,JY,JX))
  allocate(CFOMC(2,JZ,JY,JX))
  allocate(TOMT(JY,JX));        TOMT=0._r8
  allocate(TONT(JY,JX));        TONT=0._r8
  allocate(TOPT(JY,JX));        TOPT=0._r8
  allocate(URSDC(JY,JX));       URSDC=0._r8
  allocate(UORGC(JY,JX));       UORGC=0._r8
  allocate(URSDN(JY,JX));       URSDN=0._r8
  allocate(URSDP(JY,JX));       URSDP=0._r8
  allocate(UORGN(JY,JX));       UORGN=0._r8
  allocate(UORGP(JY,JX));       UORGP=0._r8
  allocate(EPOC(0:JZ,JY,JX));   EPOC=0._r8
  allocate(EHUM(0:JZ,JY,JX));   EHUM=0._r8
  allocate(CDOM(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));CDOM=0._r8
  allocate(FOSRH(1:jcplx,0:JZ,JY,JX));FOSRH=0._r8
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
  call destroy(OSC)
  call destroy(OSN)
  call destroy(OSP)
  call destroy(OHC)
  call destroy(OHN)
  call destroy(OHP)
  call destroy(OHA)
  call destroy(ORC)
  call destroy(ORN)
  call destroy(ORP)
  call destroy(DOM)
  call destroy(DOM_MacP)
  call destroy(ORGC)
  call destroy(ORGN)
  call destroy(ORGCX)
  call destroy(OSA)
  call destroy(ORGR)
  call destroy(CORGC)
  call destroy(CORGN)
  call destroy(CORGP)
  call destroy(CORGR)
  call destroy(CFOMC)
  call destroy(TOMT)
  call destroy(TONT)
  call destroy(TOPT)
  call destroy(URSDC)
  call destroy(UORGC)
  call destroy(URSDN)
  call destroy(URSDP)
  call destroy(UORGN)
  call destroy(UORGP)
  call destroy(EPOC)
  call destroy(EHUM)

  call destroy(CDOM)
  call destroy(FOSRH)
  call destroy(DIC_mass_col)
  call destroy(UNH4)
  call destroy(UNO3)
  call destroy(UPO4)
  call destroy(OMCL)
  call destroy(OMNL)
  end subroutine DestructSOMData
end module SOMDataType
