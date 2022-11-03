module SOMDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc, ndbiomcp=>ndbiomcpc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),pointer :: RSC(:,:,:,:)                       !initial surface litter [C	g m-2]
  real(r8),pointer :: RSN(:,:,:,:)                       !initial surface litter [N	g m-2]
  real(r8),pointer :: RSP(:,:,:,:)                       !initial surface litter [P	g m-2]
  real(r8),pointer :: CFOSC(:,:,:,:,:)                   !fraction of SOC in kinetic components
  real(r8),pointer :: CNOSC(:,:,:,:,:)                   !N:C,ratios of SOC kinetic components
  real(r8),pointer :: CPOSC(:,:,:,:,:)                   !P:C ratios of SOC kinetic components
  real(r8),pointer :: OSC(:,:,:,:,:)                     !humus soil C	[g d-2]
  real(r8),pointer :: OSN(:,:,:,:,:)                     !humus soil N	[g d-2]
  real(r8),pointer :: OSP(:,:,:,:,:)                     !humus soil P	[g d-2]
  real(r8),pointer :: OHC(:,:,:,:)                       !adsorbed soil C	[g d-2]
  real(r8),pointer :: OHN(:,:,:,:)                       !adsorbed soil N	[g d-2]
  real(r8),pointer :: OHP(:,:,:,:)                       !adsorbed soil P	[g d-2]
  real(r8),pointer :: OHA(:,:,:,:)                       !adsorbed soil acetate	[g d-2]
  real(r8),pointer :: ORC(:,:,:,:,:)                     !microbial residue [C	g d-2]
  real(r8),pointer :: ORN(:,:,:,:,:)                     !microbial residue N	[g d-2]
  real(r8),pointer :: ORP(:,:,:,:,:)                     !microbial residue P	[g d-2]
  real(r8),pointer :: OQC(:,:,:,:)                       !dissolved organic C micropore	[g d-2]
  real(r8),pointer :: OQN(:,:,:,:)                       !dissolved organic N micropore	[g d-2]
  real(r8),pointer :: OQP(:,:,:,:)                       !dissolved organic P micropore	[g d-2]
  real(r8),pointer :: OQA(:,:,:,:)                       !dissolved acetate micropore	[g d-2]
  real(r8),pointer :: OQCH(:,:,:,:)                      !dissolved organic C macropore	[g d-2]
  real(r8),pointer :: OQNH(:,:,:,:)                      !dissolved organic N macropore	[g d-2]
  real(r8),pointer :: OQPH(:,:,:,:)                      !dissolved organic P macropore	[g d-2]
  real(r8),pointer :: OQAH(:,:,:,:)                      !dissolved acetate macropore	[g d-2]
  real(r8),pointer :: ORGC(:,:,:)                        !total soil organic C [g d-2]
  real(r8),pointer :: ORGN(:,:,:)                        !total soil organic N [g d-2]
  real(r8),pointer :: ORGCX(:,:,:)                       !SOC concentration	[g Mg-1]
  real(r8),pointer :: OSA(:,:,:,:,:)                     !colonized humus C in each complex [g d-2]
  real(r8),pointer :: ORGR(:,:,:)                        !total particulate organic C [g d-2]
  real(r8),pointer :: CORGC(:,:,:)                       !soil organic C content [g kg-1]
  real(r8),pointer :: CORGN(:,:,:)                       !soil organic N content [mg kg-1]
  real(r8),pointer :: CORGP(:,:,:)                       !soil organic P content  [mg kg-1]
  real(r8),pointer :: CORGR(:,:,:)                       !soil particulate C content [g kg-1]
  real(r8),pointer :: CFOMC(:,:,:,:)                     !allocation coefficient to humus fractions
  real(r8),pointer ::  TOMT(:,:)                         !total micriobial C, [g d-2]
  real(r8),pointer ::  TONT(:,:)                         !total micriobial N, [g d-2]
  real(r8),pointer ::  TOPT(:,:)                         !total micriobial P, [g d-2]
  real(r8),pointer ::  URSDC(:,:)                        !total litter C, [g d-2]
  real(r8),pointer ::  UORGC(:,:)                        !total humus C, [g d-2]
  real(r8),pointer ::  URSDN(:,:)                        !total litter N, [g d-2]
  real(r8),pointer ::  URSDP(:,:)                        !total litter P, [g d-2]
  real(r8),pointer ::  UORGN(:,:)                        !total humus N, [g d-2]
  real(r8),pointer ::  UORGP(:,:)                        !total humus P, [g d-2]
  real(r8),pointer ::  EPOC(:,:,:)                       !partitioning coefficient between POC and litter, []
  real(r8),pointer ::  EHUM(:,:,:)                       !partitioning coefficient between humus and microbial residue, []
  real(r8),pointer ::  COQC(:,:,:,:)                     !DOC concentration, [g m-3]
  real(r8),pointer ::  COQA(:,:,:,:)                     !acetate concentration, [g m-3]
  real(r8),pointer ::  FOSRH(:,:,:,:)                    !fraction of total organic C in complex, [-]
  real(r8),pointer ::  UCO2S(:,:)                        !total soil DIC, [g d-2]
  real(r8),pointer ::  UNH4(:,:)                         !total soil NH4 + NH3 content, [g d-2]
  real(r8),pointer ::  UNO3(:,:)                         !total soil NO3 + NO2 content, [g d-2]
  real(r8),pointer ::  UPO4(:,:)                         !total soil PO4 content, [g d-2]
  real(r8),pointer ::  OMCL(:,:,:)
  real(r8),pointer ::  OMNL(:,:,:)

  private :: InitAllocate
  contains

  subroutine InitSOMData(n_litrsfk)

  implicit none
  integer, intent(in) :: n_litrsfk

  call InitAllocate(n_litrsfk)

  end subroutine InitSOMData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate(n_litrsfk)
  implicit none
  integer, intent(in) :: n_litrsfk

  allocate(RSC(1:n_litrsfk,0:JZ,JY,JX))
  allocate(RSN(1:n_litrsfk,0:JZ,JY,JX))
  allocate(RSP(1:n_litrsfk,0:JZ,JY,JX))
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
  allocate(OQC(1:jcplx,0:JZ,JY,JX))
  allocate(OQN(1:jcplx,0:JZ,JY,JX))
  allocate(OQP(1:jcplx,0:JZ,JY,JX))
  allocate(OQA(1:jcplx,0:JZ,JY,JX))
  allocate(OQCH(1:jcplx,0:JZ,JY,JX))
  allocate(OQNH(1:jcplx,0:JZ,JY,JX))
  allocate(OQPH(1:jcplx,0:JZ,JY,JX))
  allocate(OQAH(1:jcplx,0:JZ,JY,JX))
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
  allocate(COQC(1:jcplx,0:JZ,JY,JX));COQC=0._r8
  allocate(COQA(1:jcplx,0:JZ,JY,JX));COQA=0._r8
  allocate(FOSRH(1:jcplx,0:JZ,JY,JX));FOSRH=0._r8
  allocate(UCO2S(JY,JX));       UCO2S=0._r8
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
  call destroy(OQC)
  call destroy(OQN)
  call destroy(OQP)
  call destroy(OQA)
  call destroy(OQCH)
  call destroy(OQNH)
  call destroy(OQPH)
  call destroy(OQAH)
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
  call destroy(COQC)
  call destroy(COQA)
  call destroy(FOSRH)
  call destroy(UCO2S)
  call destroy(UNH4)
  call destroy(UNO3)
  call destroy(UPO4)
  call destroy(OMCL)
  call destroy(OMNL)
  end subroutine DestructSOMData
end module SOMDataType
