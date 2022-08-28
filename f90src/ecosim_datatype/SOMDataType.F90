module SOMDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable :: RSC(:,:,:,:)                       !initial surface litter [C	g m-2]
  real(r8),allocatable :: RSN(:,:,:,:)                       !initial surface litter [N	g m-2]
  real(r8),allocatable :: RSP(:,:,:,:)                       !initial surface litter [P	g m-2]
  real(r8),allocatable :: CFOSC(:,:,:,:,:)                   !fraction of SOC in kinetic components
  real(r8),allocatable :: CNOSC(:,:,:,:,:)                   !N:C,ratios of SOC kinetic components
  real(r8),allocatable :: CPOSC(:,:,:,:,:)                   !P:C ratios of SOC kinetic components
  real(r8),allocatable :: OSC(:,:,:,:,:)                     !humus soil C	[g d-2]
  real(r8),allocatable :: OSN(:,:,:,:,:)                     !humus soil N	[g d-2]
  real(r8),allocatable :: OSP(:,:,:,:,:)                     !humus soil P	[g d-2]
  real(r8),allocatable :: OHC(:,:,:,:)                       !adsorbed soil C	[g d-2]
  real(r8),allocatable :: OHN(:,:,:,:)                       !adsorbed soil N	[g d-2]
  real(r8),allocatable :: OHP(:,:,:,:)                       !adsorbed soil P	[g d-2]
  real(r8),allocatable :: OHA(:,:,:,:)                       !adsorbed soil acetate	[g d-2]
  real(r8),allocatable :: ORC(:,:,:,:,:)                     !microbial residue [C	g d-2]
  real(r8),allocatable :: ORN(:,:,:,:,:)                     !microbial residue N	[g d-2]
  real(r8),allocatable :: ORP(:,:,:,:,:)                     !microbial residue P	[g d-2]
  real(r8),allocatable :: OQC(:,:,:,:)                       !dissolved organic C micropore	[g d-2]
  real(r8),allocatable :: OQN(:,:,:,:)                       !dissolved organic N micropore	[g d-2]
  real(r8),allocatable :: OQP(:,:,:,:)                       !dissolved organic P micropore	[g d-2]
  real(r8),allocatable :: OQA(:,:,:,:)                       !dissolved acetate micropore	[g d-2]
  real(r8),allocatable :: OQCH(:,:,:,:)                      !dissolved organic C macropore	[g d-2]
  real(r8),allocatable :: OQNH(:,:,:,:)                      !dissolved organic N macropore	[g d-2]
  real(r8),allocatable :: OQPH(:,:,:,:)                      !dissolved organic P macropore	[g d-2]
  real(r8),allocatable :: OQAH(:,:,:,:)                      !dissolved acetate macropore	[g d-2]
  real(r8),allocatable :: ORGC(:,:,:)                        !total soil organic C [g d-2]
  real(r8),allocatable :: ORGN(:,:,:)                        !total soil organic N [g d-2]
  real(r8),allocatable :: ORGCX(:,:,:)                       !SOC concentration	[g Mg-1]
  real(r8),allocatable :: OSA(:,:,:,:,:)                     !colonized humus C in each complex [g d-2]
  real(r8),allocatable :: ORGR(:,:,:)                        !total particulate organic C [g d-2]
  real(r8),allocatable :: CORGC(:,:,:)                       !soil organic C content [g kg-1]
  real(r8),allocatable :: CORGN(:,:,:)                       !soil organic N content [mg kg-1]
  real(r8),allocatable :: CORGP(:,:,:)                       !soil organic P content  [mg kg-1]
  real(r8),allocatable :: CORGR(:,:,:)                       !soil particulate C content [g kg-1]
  real(r8),allocatable :: CFOMC(:,:,:,:)                     !allocation coefficient to humus fractions
  real(r8),allocatable ::  TOMT(:,:)                         !total micriobial C, [g d-2]
  real(r8),allocatable ::  TONT(:,:)                         !total micriobial N, [g d-2]
  real(r8),allocatable ::  TOPT(:,:)                         !total micriobial P, [g d-2]
  real(r8),allocatable ::  URSDC(:,:)                        !total litter C, [g d-2]
  real(r8),allocatable ::  UORGC(:,:)                        !total humus C, [g d-2]
  real(r8),allocatable ::  URSDN(:,:)                        !total litter N, [g d-2]
  real(r8),allocatable ::  URSDP(:,:)                        !total litter P, [g d-2]
  real(r8),allocatable ::  UORGN(:,:)                        !total humus N, [g d-2]
  real(r8),allocatable ::  UORGP(:,:)                        !total humus P, [g d-2]
  real(r8),allocatable ::  EPOC(:,:,:)                       !partitioning coefficient between POC and litter, []
  real(r8),allocatable ::  EHUM(:,:,:)                       !partitioning coefficient between humus and microbial residue, []
  real(r8),allocatable ::  COQC(:,:,:,:)                     !DOC concentration, [g m-3]
  real(r8),allocatable ::  COQA(:,:,:,:)                     !acetate concentration, [g m-3]
  real(r8),allocatable ::  FOSRH(:,:,:,:)                    !fraction of total organic C in complex, [-]
  real(r8),allocatable ::  UCO2S(:,:)                        !total soil DIC, [g d-2]
  real(r8),allocatable ::  UNH4(:,:)                         !total soil NH4 + NH3 content, [g d-2]
  real(r8),allocatable ::  UNO3(:,:)                         !total soil NO3 + NO2 content, [g d-2]
  real(r8),allocatable ::  UPO4(:,:)                         !total soil PO4 content, [g d-2]
  REAL(R8),ALLOCATABLE ::  OMCL(:,:,:)
  REAL(R8),ALLOCATABLE ::  OMNL(:,:,:)

  private :: InitAllocate
  contains

  subroutine InitSOMData

  implicit none

  call InitAllocate

  end subroutine InitSOMData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  allocate(RSC(0:2,0:JZ,JY,JX))
  allocate(RSN(0:2,0:JZ,JY,JX))
  allocate(RSP(0:2,0:JZ,JY,JX))
  allocate(CFOSC(jsken,0:jcplx1,0:JZ,JY,JX))
  allocate(CNOSC(jsken,0:jcplx1,0:JZ,JY,JX))
  allocate(CPOSC(jsken,0:jcplx1,0:JZ,JY,JX))
  allocate(OSC(jsken,0:jcplx1,0:JZ,JY,JX))
  allocate(OSN(jsken,0:jcplx1,0:JZ,JY,JX))
  allocate(OSP(jsken,0:jcplx1,0:JZ,JY,JX))
  allocate(OHC(0:jcplx1,0:JZ,JY,JX))
  allocate(OHN(0:jcplx1,0:JZ,JY,JX))
  allocate(OHP(0:jcplx1,0:JZ,JY,JX))
  allocate(OHA(0:jcplx1,0:JZ,JY,JX))
  allocate(ORC(2,0:jcplx1,0:JZ,JY,JX))
  allocate(ORN(2,0:jcplx1,0:JZ,JY,JX))
  allocate(ORP(2,0:jcplx1,0:JZ,JY,JX))
  allocate(OQC(0:jcplx1,0:JZ,JY,JX))
  allocate(OQN(0:jcplx1,0:JZ,JY,JX))
  allocate(OQP(0:jcplx1,0:JZ,JY,JX))
  allocate(OQA(0:jcplx1,0:JZ,JY,JX))
  allocate(OQCH(0:jcplx1,0:JZ,JY,JX))
  allocate(OQNH(0:jcplx1,0:JZ,JY,JX))
  allocate(OQPH(0:jcplx1,0:JZ,JY,JX))
  allocate(OQAH(0:jcplx1,0:JZ,JY,JX))
  allocate(ORGC(0:JZ,JY,JX))
  allocate(ORGN(0:JZ,JY,JX))
  allocate(ORGCX(0:JZ,JY,JX))
  allocate(OSA(jsken,0:jcplx1,0:JZ,JY,JX))
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
  allocate(COQC(0:jcplx1,0:JZ,JY,JX));COQC=0._r8
  allocate(COQA(0:jcplx1,0:JZ,JY,JX));COQA=0._r8
  allocate(FOSRH(0:jcplx1,0:JZ,JY,JX));FOSRH=0._r8
  allocate(UCO2S(JY,JX));       UCO2S=0._r8
  allocate(UNH4(JY,JX));        UNH4=0._r8
  allocate(UNO3(JY,JX));        UNO3=0._r8
  allocate(UPO4(JY,JX));        UPO4=0._r8
  ALLOCATE(OMCL(0:JZ,JY,JX));   OMCL=0._r8
  ALLOCATE(OMNL(0:JZ,JY,JX));   OMNL=0._r8
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------

  subroutine DestructSOMData

  implicit none

  if(allocated(RSC))deallocate(RSC)
  if(allocated(RSN))deallocate(RSN)
  if(allocated(RSP))deallocate(RSP)
  if(allocated(CFOSC))deallocate(CFOSC)
  if(allocated(CNOSC))deallocate(CNOSC)
  if(allocated(CPOSC))deallocate(CPOSC)
  if(allocated(OSC))deallocate(OSC)
  if(allocated(OSN))deallocate(OSN)
  if(allocated(OSP))deallocate(OSP)
  if(allocated(OHC))deallocate(OHC)
  if(allocated(OHN))deallocate(OHN)
  if(allocated(OHP))deallocate(OHP)
  if(allocated(OHA))deallocate(OHA)
  if(allocated(ORC))deallocate(ORC)
  if(allocated(ORN))deallocate(ORN)
  if(allocated(ORP))deallocate(ORP)
  if(allocated(OQC))deallocate(OQC)
  if(allocated(OQN))deallocate(OQN)
  if(allocated(OQP))deallocate(OQP)
  if(allocated(OQA))deallocate(OQA)
  if(allocated(OQCH))deallocate(OQCH)
  if(allocated(OQNH))deallocate(OQNH)
  if(allocated(OQPH))deallocate(OQPH)
  if(allocated(OQAH))deallocate(OQAH)
  if(allocated(ORGC))deallocate(ORGC)
  if(allocated(ORGN))deallocate(ORGN)
  if(allocated(ORGCX))deallocate(ORGCX)
  if(allocated(OSA))deallocate(OSA)
  if(allocated(ORGR))deallocate(ORGR)
  if(allocated(CORGC))deallocate(CORGC)
  if(allocated(CORGN))deallocate(CORGN)
  if(allocated(CORGP))deallocate(CORGP)
  if(allocated(CORGR))deallocate(CORGR)
  if(allocated(CFOMC))deallocate(CFOMC)
  if (allocated(TOMT))     deallocate(TOMT)
  if (allocated(TONT))     deallocate(TONT)
  if (allocated(TOPT))     deallocate(TOPT)
  if (allocated(URSDC))    deallocate(URSDC)
  if (allocated(UORGC))    deallocate(UORGC)
  if (allocated(URSDN))    deallocate(URSDN)
  if (allocated(URSDP))    deallocate(URSDP)
  if (allocated(UORGN))    deallocate(UORGN)
  if (allocated(UORGP))    deallocate(UORGP)
  if (allocated(EPOC))     deallocate(EPOC)
  if (allocated(EHUM))     deallocate(EHUM)
  if (allocated(COQC))     deallocate(COQC)
  if (allocated(COQA))     deallocate(COQA)
  if (allocated(FOSRH))    deallocate(FOSRH)
  if (allocated(UCO2S))    deallocate(UCO2S)
  if (allocated(UNH4))     deallocate(UNH4)
  if (allocated(UNO3))     deallocate(UNO3)
  if (allocated(UPO4))     deallocate(UPO4)
  if(allocated(OMCL))   deallocate(OMCL)
  if(allocated(OMNL)) deallocate(OMNL)
  end subroutine DestructSOMData
end module SOMDataType
