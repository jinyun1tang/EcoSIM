module PlantTraitDataType

!
!!
! data types of plant trait characteristics that cannot be grouped into canopy
! or roots
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__


!allocation parameter
  real(r8) :: FVRN(0:5)              !allocation parameter
  real(r8) :: FWOOD(0:1)             !woody C allocation
  real(r8) :: FWOODN(0:1)            !woody N allocation
  real(r8) :: FWOODP(0:1)            !woody P allocation
  REAL(R8) :: FWODB(0:1)             !woody C allocation
  real(r8) :: FWODLN(0:1)            !allocation
  real(r8) :: FWODLP(0:1)
  REAL(R8) :: FWODSN(0:1)            !N woody fraction in petiole
  real(r8) :: FWODSP(0:1)            !P woody fraction in petiole
  real(r8) :: FWODR(0:1)
  real(r8) :: FWODRN(0:1)
  real(r8) :: FWODRP(0:1)


  real(r8),allocatable ::  ARSTK(:,:,:,:,:)                   !stem layer area, [m2 d-2]
  real(r8),allocatable ::  ARLFP(:,:,:)                       !plant leaf area, [m2 d-2]
  real(r8),allocatable ::  ARLFS(:,:,:)                       !plant leaf area, [m2 d-2]
  real(r8),allocatable ::  ARLFX(:,:)                         !total canopy leaf area, [m2 d-2]
  real(r8),allocatable ::  ARSTP(:,:,:)                       !plant stem area, [m2 d-2]
  real(r8),allocatable ::  ZC(:,:,:)                          !canopy height, [m]
  real(r8),allocatable ::  ARSTX(:,:)                         !total canopy stem area, [m2 d-2]
  real(r8),allocatable ::  ARLFT(:,:,:)                       !total leaf area, [m2 d-2]
  real(r8),allocatable ::  ARSTT(:,:,:)                       !total stem area, [m2 d-2]
  real(r8),allocatable ::  ARLFC(:,:)                         !total canopy leaf area, [m2 d-2]
  real(r8),allocatable ::  ARSTC(:,:)                         !total canopy stem area, [m2 d-2]
  real(r8),allocatable ::  ARLSS(:,:)                         !stalk area of combined,each PFT canopy
  integer ,allocatable ::  NG(:,:,:)                           !soil layer at planting depth, [-]
  real(r8),allocatable ::  SDPTHI(:,:,:)                      !planting depth, [m]
  real(r8),allocatable ::  SDPTH(:,:,:)                       !seeding depth, [m]
  real(r8),allocatable ::  SDVL(:,:,:)                        !seed volume, [m3 ]
  real(r8),allocatable ::  SDLG(:,:,:)                        !seed length, [m]
  real(r8),allocatable ::  SDAR(:,:,:)                        !seed surface area, [m2]
  real(r8),allocatable ::  HTCTL(:,:,:)                       !cotyledon height, [m]
  real(r8),allocatable ::  ZT(:,:)                            !canopy height , [m]
  real(r8),allocatable ::  ZL(:,:,:)                          !canopy layer height , [m]
  real(r8),allocatable ::  ANGBR(:,:,:)                       !branching angle, [degree from horizontal]
  real(r8),allocatable ::  ANGSH(:,:,:)                       !sheath angle, [degree from horizontal]
  real(r8),allocatable ::  RAZ(:,:,:)                         !canopy roughness height, [m]
  real(r8),allocatable ::  HTSTZ(:,:,:)                       !canopy height, [m]
  real(r8),allocatable ::  ARLF(:,:,:,:,:)                    !leaf area, [m2 d-2]
  real(r8),allocatable ::  HTSHE(:,:,:,:,:)                   !sheath height, [m]
  real(r8),allocatable ::  HTNODE(:,:,:,:,:)                  !internode height, [m]
  real(r8),allocatable ::  ARLFB(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),allocatable ::  ARLFZ(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),allocatable ::  HTSHEX(:,:,:,:)                    !branch height, [m]
  real(r8),allocatable ::  GRNOB(:,:,:,:)                     !branch grain number, [d-2]
  real(r8),allocatable ::  GRNXB(:,:,:,:)                     !branch potential grain number, [d-2]
  real(r8),allocatable ::  GRNO(:,:,:)                        !canopy grain number, [d-2]
  real(r8),allocatable ::  PP(:,:,:)                          !plant population, [d-2]
  real(r8),allocatable ::  HTNODX(:,:,:,:,:)                  !internode height, [m]
  real(r8),allocatable ::  CNLF(:,:,:)                        !maximum leaf N:C ratio, [g g-1]
  real(r8),allocatable ::  CPLF(:,:,:)                        !maximum leaf P:C ratio, [g g-1]
  real(r8),allocatable ::  CNSHE(:,:,:)                       !sheath N:C ratio, [g g-1]
  real(r8),allocatable ::  CNSTK(:,:,:)                       !stalk N:C ratio, [g g-1]
  real(r8),allocatable ::  CNRSV(:,:,:)                       !reserve N:C ratio, [g g-1]
  real(r8),allocatable ::  CNHSK(:,:,:)                       !husk N:C ratio, [g g-1]
  real(r8),allocatable ::  CNEAR(:,:,:)                       !ear N:C ratio, [g g-1]
  real(r8),allocatable ::  CNGR(:,:,:)                        !grain N:C ratio, [g g-1]
  real(r8),allocatable ::  CNND(:,:,:)                        !nodule N:C ratio, [g g-1]
  real(r8),allocatable ::  CPSHE(:,:,:)                       !sheath P:C ratio, [g g-1]
  real(r8),allocatable ::  CPSTK(:,:,:)                       !stalk P:C ratio, [g g-1]
  real(r8),allocatable ::  CPRSV(:,:,:)                       !reserve P:C ratio, [g g-1]
  real(r8),allocatable ::  CPHSK(:,:,:)                       !husk P:C ratio, [g g-1]
  real(r8),allocatable ::  CPEAR(:,:,:)                       !ear P:C ratio, [g g-1]
  real(r8),allocatable ::  CPGR(:,:,:)                        !grain P:C ratio, [g g-1]
  real(r8),allocatable ::  CPND(:,:,:)                        !nodule P:C ratio, [g g-1]
  real(r8),allocatable ::  CNWS(:,:,:)                        !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8),allocatable ::  CPWS(:,:,:)                        !C:P ratio in remobilizable nonstructural biomass, [-]
  real(r8),allocatable ::  OSMO(:,:,:)                        !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8),allocatable ::  TCX(:,:,:)                         !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8),allocatable ::  ZTYPI(:,:,:)                       !initial plant thermal adaptation zone, [-]
  real(r8),allocatable ::  ZTYP(:,:,:)                        !plant thermal adaptation zone, [-]
  real(r8),allocatable ::  GROUP(:,:,:,:)                     !plant maturity group, [-]
  real(r8),allocatable ::  GROUPI(:,:,:)                      !acclimated plant maturity group, [-]
  real(r8),allocatable ::  GROUPX(:,:,:)                      !initial plant maturity group, [-]
  real(r8),allocatable ::  PPI(:,:,:)                         !initial plant population, [m-2]
  real(r8),allocatable ::  WTSTDI(:,:,:)                      !initial standing dead C, [g C m-2]
  real(r8),allocatable ::  PPX(:,:,:)                         !plant population, [m-2]
  integer,allocatable ::  IFLGT(:,:)                          !number of active PFT
  real(r8),allocatable ::  PPT(:,:)                           !total plant population, [d-2]
  real(r8),allocatable ::  PPZ(:,:,:)                         !plant population at seeding, [m-2]
  real(r8),allocatable ::  WSTR(:,:,:)                        !canopy plant water stress indicator, number of hours PSILT < PSILY, []
  real(r8),allocatable ::  OSTR(:,:,:)                        !plant O2 stress indicator, []
  real(r8),allocatable ::  TFN3(:,:,:)                        !canopy temperature growth function, [-]
  real(r8),allocatable ::  TCG(:,:,:)                         !canopy growth temperature, [oC]
  real(r8),allocatable ::  TKG(:,:,:)                         !canopy growth temperature, [K]
  real(r8),allocatable ::  DMSHE(:,:,:)                       !sheath growth yield, [g g-1]
  real(r8),allocatable ::  DMSTK(:,:,:)                       !stalk growth yield, [g g-1]
  real(r8),allocatable ::  DMRSV(:,:,:)                       !reserve growth yield, [g g-1]
  real(r8),allocatable ::  DMHSK(:,:,:)                       !husk growth yield, [g g-1]
  real(r8),allocatable ::  DMEAR(:,:,:)                       !ear growth yield, [g g-1]
  real(r8),allocatable ::  DMGR(:,:,:)                        !grain growth yield, [g g-1]
  real(r8),allocatable ::  DMND(:,:,:)                        !nodule growth yield, [g g-1]
  real(r8),allocatable ::  DMLF(:,:,:)                        !leaf growth yield, [g g-1]
  real(r8),allocatable ::  VRNY(:,:,:,:)                      !initial heat requirement for spring leafout/dehardening, [h]
  real(r8),allocatable ::  VRNZ(:,:,:,:)                      !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8),allocatable ::  VSTG(:,:,:,:)                      !leaf number, [-]
  real(r8),allocatable ::  VSTGX(:,:,:,:)                     !leaf number at floral initiation, [-]
  real(r8),allocatable ::  VRNS(:,:,:,:)                      !heat requirement for spring leafout/dehardening, [h]
  real(r8),allocatable ::  VRNF(:,:,:,:)                      !cold requirement for autumn leafoff/hardening, [h]
  integer,allocatable ::  KLEAF(:,:,:,:)                      !leaf number, [-]
  integer,allocatable ::  KVSTGN(:,:,:,:)                     !leaf growth stage counter, [-]
  integer,allocatable ::  KLEAFX(:,:,:,:)                     !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,allocatable ::  KVSTG(:,:,:,:)                      !leaf growth stage counter, [-]
  real(r8),allocatable ::  XRLA(:,:,:)                        !rate of leaf initiation, [h-1 at 25 oC]
  real(r8),allocatable ::  WDLF(:,:,:)                        !leaf length:width ratio, [-]
  real(r8),allocatable ::  SLA1(:,:,:)                        !leaf area:mass during growth, [m2 g-1]
  real(r8),allocatable ::  TCZ(:,:,:)                         !threshold temperature for spring leafout/dehardening, [oC]
  real(r8),allocatable ::  SSL1(:,:,:)                        !petiole length:mass during growth, [m g-1]
  real(r8),allocatable ::  VRNL(:,:,:,:)                      !hours above threshold temperature required for spring leafout/dehardening, [-]
  integer,allocatable ::  NBR(:,:,:)                          !branch number, [-]
  integer,allocatable ::  NBT(:,:,:)                          !branch number, [-]
  integer,allocatable ::  NBTB(:,:,:,:)                       !branch number, [-]
  integer,allocatable ::  NB1(:,:,:)                          !number of main branch, [-]
  integer,allocatable ::  IFLGR(:,:,:,:)                      !branch phenology flag, [-]
  integer,allocatable ::  IFLGQ(:,:,:,:)                      !branch phenology flag, [h]
  integer,allocatable ::  IFLGG(:,:,:,:)                      !branch phenology flag, [-]
  integer,allocatable ::  IFLGP(:,:,:,:)                      !branch phenology flag, [-]
  integer,allocatable ::  IFLGA(:,:,:,:)                      !branch phenology flag, [-]
  integer,allocatable ::  IFLGE(:,:,:,:)                      !branch phenology flag, [-]
  integer,allocatable ::  IFLGF(:,:,:,:)                      !branch phenology flag, [-]
  integer,allocatable ::  IDTHB(:,:,:,:)                      !flag to detect branch death , [-]
  real(r8),allocatable ::  PB(:,:,:)                          !branch nonstructural C content required for new branch, [g g-1]
  real(r8),allocatable ::  GSTGI(:,:,:,:)                     !normalized node number during vegetative growth stages , [-]
  real(r8),allocatable ::  DGSTGI(:,:,:,:)                    !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8),allocatable ::  DGSTGF(:,:,:,:)                    !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8),allocatable ::  PSTG(:,:,:,:)                      !node number, [-]
  real(r8),allocatable ::  PSTGI(:,:,:,:)                     !node number at floral initiation, [-]
  real(r8),allocatable ::  GSTGF(:,:,:,:)                     !normalized node number during reproductive growth stages, [-]
  real(r8),allocatable ::  PSTGF(:,:,:,:)                     !node number at anthesis, [-]
  real(r8),allocatable ::  TGSTGI(:,:,:,:)                    !normalized node number during vegetative growth stages , [-]
  real(r8),allocatable ::  TGSTGF(:,:,:,:)                    !normalized node number during reproductive growth stages , [-]
  real(r8),allocatable ::  XRNI(:,:,:)                        !rate of node initiation, [h-1 at 25 oC]
  real(r8),allocatable ::  SNL1(:,:,:)                        !internode length:mass during growth, [m g-1]
  real(r8),allocatable ::  FNOD(:,:,:)                        !parameter for allocation of growth to nodes, [-]
  integer,allocatable ::  NNOD(:,:,:)                         !number of concurrently growing nodes
  real(r8),allocatable ::  PSILZ(:,:,:)                       !minimum daily canopy water potential, [MPa]
  real(r8),allocatable ::  CFX(:,:,:)                         !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8),allocatable ::  CF(:,:,:)                          !clumping factor for self-shading in canopy layer, [-]
  integer,allocatable ::  IDTHP(:,:,:)                        !flag to detect canopy death , [-]
  real(r8),allocatable ::  STMX(:,:,:)                        !maximum grain node number per branch, [-]
  real(r8),allocatable ::  SDMX(:,:,:)                        !maximum grain number per node , [-]
  real(r8),allocatable ::  GRMX(:,:,:)                        !maximum grain size   , [g]
  real(r8),allocatable ::  XTLI(:,:,:)                        !number of nodes in seed, [-]
  real(r8),allocatable ::  GRDM(:,:,:)                        !grain size at seeding, [g]
  real(r8),allocatable ::  GFILL(:,:,:)                       !maximum rate of fill per grain, [g h-1]
  real(r8),allocatable ::  FLG4(:,:,:,:)                      !flag to detect physiological maturity from  grain fill , [-]
  real(r8),allocatable ::  ATRP(:,:,:,:)                      !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8),allocatable ::  FLGZ(:,:,:,:)                      !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  integer,allocatable ::  IDAY(:,:,:,:,:)                     !plant growth stage, [-]
  real(r8),allocatable ::  CTC(:,:,:)                         !temperature below which seed set is adversely affected, [oC]
  real(r8),allocatable ::  HTC(:,:,:)                         !temperature above which seed set is adversely affected, [oC]
  real(r8),allocatable ::  SSTX(:,:,:)                        !sensitivity to HTC (seeds oC-1 above HTC)
  real(r8),allocatable ::  XDL(:,:,:)                         !critical daylength for phenological progress, [h]
  real(r8),allocatable ::  XPPD(:,:,:)                        !difference between current and critical daylengths used to calculate  phenological progress, [h]
  real(r8),allocatable ::  CFI(:,:,:)                         !initial clumping factor for self-shading in canopy layer, [-]
  real(r8),allocatable ::  VRNX(:,:,:,:)                      !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8),allocatable ::  OFFST(:,:,:)                       !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
!----------------------------------------------------------------------

contains
  subroutine InitPlantTraits

  implicit none

  FVRN =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)

  allocate(ARSTK(JC,JC,JP,JY,JX));ARSTK=0._r8
  allocate(ARLFP(JP,JY,JX));    ARLFP=0._r8
  allocate(ARLFS(JP,JY,JX));    ARLFS=0._r8
  allocate(ARLFX(JY,JX));       ARLFX=0._r8
  allocate(ARSTP(JP,JY,JX));    ARSTP=0._r8
  allocate(ZC(JP,JY,JX));       ZC=0._r8
  allocate(ARSTX(JY,JX));       ARSTX=0._r8
  allocate(ARLFT(JC,JY,JX));    ARLFT=0._r8
  allocate(ARSTT(JC,JY,JX));    ARSTT=0._r8
  allocate(ARLFC(JY,JX));       ARLFC=0._r8
  allocate(ARSTC(JY,JX));       ARSTC=0._r8
  allocate(ARLSS(JY,JX));       ARLSS=0._r8
  allocate(NG(JP,JY,JX));       NG=0
  allocate(SDPTHI(JP,JY,JX));   SDPTHI=0._r8
  allocate(SDPTH(JP,JY,JX));    SDPTH=0._r8
  allocate(SDVL(JP,JY,JX));     SDVL=0._r8
  allocate(SDLG(JP,JY,JX));     SDLG=0._r8
  allocate(SDAR(JP,JY,JX));     SDAR=0._r8
  allocate(HTCTL(JP,JY,JX));    HTCTL=0._r8
  allocate(ZT(JY,JX));          ZT=0._r8
  allocate(ZL(0:JC,JY,JX));     ZL=0._r8
  allocate(ANGBR(JP,JY,JX));    ANGBR=0._r8
  allocate(ANGSH(JP,JY,JX));    ANGSH=0._r8
  allocate(RAZ(JP,JY,JX));      RAZ=0._r8
  allocate(HTSTZ(JP,JY,JX));    HTSTZ=0._r8
  allocate(ARLF(0:JNODS,JC,JP,JY,JX));ARLF=0._r8
  allocate(HTSHE(0:JNODS,JC,JP,JY,JX));HTSHE=0._r8
  allocate(HTNODE(0:JNODS,JC,JP,JY,JX));HTNODE=0._r8
  allocate(ARLFB(JC,JP,JY,JX)); ARLFB=0._r8
  allocate(ARLFZ(JC,JP,JY,JX)); ARLFZ=0._r8
  allocate(HTSHEX(JC,JP,JY,JX));HTSHEX=0._r8
  allocate(GRNOB(JC,JP,JY,JX)); GRNOB=0._r8
  allocate(GRNXB(JC,JP,JY,JX)); GRNXB=0._r8
  allocate(GRNO(JP,JY,JX));     GRNO=0._r8
  allocate(PP(JP,JY,JX));       PP=0._r8
  allocate(HTNODX(0:JNODS,JC,JP,JY,JX));HTNODX=0._r8
  allocate(CNLF(JP,JY,JX));     CNLF=0._r8
  allocate(CPLF(JP,JY,JX));     CPLF=0._r8
  allocate(CNSHE(JP,JY,JX));    CNSHE=0._r8
  allocate(CNSTK(JP,JY,JX));    CNSTK=0._r8
  allocate(CNRSV(JP,JY,JX));    CNRSV=0._r8
  allocate(CNHSK(JP,JY,JX));    CNHSK=0._r8
  allocate(CNEAR(JP,JY,JX));    CNEAR=0._r8
  allocate(CNGR(JP,JY,JX));     CNGR=0._r8
  allocate(CNND(JP,JY,JX));     CNND=0._r8
  allocate(CPSHE(JP,JY,JX));    CPSHE=0._r8
  allocate(CPSTK(JP,JY,JX));    CPSTK=0._r8
  allocate(CPRSV(JP,JY,JX));    CPRSV=0._r8
  allocate(CPHSK(JP,JY,JX));    CPHSK=0._r8
  allocate(CPEAR(JP,JY,JX));    CPEAR=0._r8
  allocate(CPGR(JP,JY,JX));     CPGR=0._r8
  allocate(CPND(JP,JY,JX));     CPND=0._r8
  allocate(CNWS(JP,JY,JX));     CNWS=0._r8
  allocate(CPWS(JP,JY,JX));     CPWS=0._r8
  allocate(OSMO(JP,JY,JX));     OSMO=0._r8
  allocate(TCX(JP,JY,JX));      TCX=0._r8
  allocate(ZTYPI(JP,JY,JX));    ZTYPI=0._r8
  allocate(ZTYP(JP,JY,JX));     ZTYP=0._r8
  allocate(GROUP(JC,JP,JY,JX)); GROUP=0._r8
  allocate(GROUPI(JP,JY,JX));   GROUPI=0._r8
  allocate(GROUPX(JP,JY,JX));   GROUPX=0._r8
  allocate(PPI(JP,JY,JX));      PPI=0._r8
  allocate(WTSTDI(JP,JY,JX));   WTSTDI=0._r8
  allocate(PPX(JP,JY,JX));      PPX=0._r8
  allocate(IFLGT(JY,JX));       IFLGT=0
  allocate(PPT(JY,JX));         PPT=0._r8
  allocate(PPZ(JP,JY,JX));      PPZ=0._r8
  allocate(WSTR(JP,JY,JX));     WSTR=0._r8
  allocate(OSTR(JP,JY,JX));     OSTR=0._r8
  allocate(TFN3(JP,JY,JX));     TFN3=0._r8
  allocate(TCG(JP,JY,JX));      TCG=0._r8
  allocate(TKG(JP,JY,JX));      TKG=0._r8
  allocate(DMSHE(JP,JY,JX));    DMSHE=0._r8
  allocate(DMSTK(JP,JY,JX));    DMSTK=0._r8
  allocate(DMRSV(JP,JY,JX));    DMRSV=0._r8
  allocate(DMHSK(JP,JY,JX));    DMHSK=0._r8
  allocate(DMEAR(JP,JY,JX));    DMEAR=0._r8
  allocate(DMGR(JP,JY,JX));     DMGR=0._r8
  allocate(DMND(JP,JY,JX));     DMND=0._r8
  allocate(DMLF(JP,JY,JX));     DMLF=0._r8
  allocate(VRNY(JC,JP,JY,JX));  VRNY=0._r8
  allocate(VRNZ(JC,JP,JY,JX));  VRNZ=0._r8
  allocate(VSTG(JC,JP,JY,JX));  VSTG=0._r8
  allocate(VSTGX(JC,JP,JY,JX)); VSTGX=0._r8
  allocate(VRNS(JC,JP,JY,JX));  VRNS=0._r8
  allocate(VRNF(JC,JP,JY,JX));  VRNF=0._r8
  allocate(KLEAF(JC,JP,JY,JX)); KLEAF=0
  allocate(KVSTGN(JC,JP,JY,JX));KVSTGN=0
  allocate(KLEAFX(JC,JP,JY,JX));KLEAFX=0
  allocate(KVSTG(JC,JP,JY,JX)); KVSTG=0
  allocate(XRLA(JP,JY,JX));     XRLA=0._r8
  allocate(WDLF(JP,JY,JX));     WDLF=0._r8
  allocate(SLA1(JP,JY,JX));     SLA1=0._r8
  allocate(TCZ(JP,JY,JX));      TCZ=0._r8
  allocate(SSL1(JP,JY,JX));     SSL1=0._r8
  allocate(VRNL(JC,JP,JY,JX));  VRNL=0._r8
  allocate(NBR(JP,JY,JX));      NBR=0
  allocate(NBT(JP,JY,JX));      NBT=0
  allocate(NBTB(JC,JP,JY,JX));  NBTB=0
  allocate(NB1(JP,JY,JX));      NB1=0
  allocate(IFLGR(JC,JP,JY,JX)); IFLGR=0
  allocate(IFLGQ(JC,JP,JY,JX)); IFLGQ=0
  allocate(IFLGG(JC,JP,JY,JX)); IFLGG=0
  allocate(IFLGP(JC,JP,JY,JX)); IFLGP=0
  allocate(IFLGA(JC,JP,JY,JX)); IFLGA=0
  allocate(IFLGE(JC,JP,JY,JX)); IFLGE=0
  allocate(IFLGF(JC,JP,JY,JX)); IFLGF=0
  allocate(IDTHB(JC,JP,JY,JX)); IDTHB=0
  allocate(PB(JP,JY,JX));       PB=0._r8
  allocate(GSTGI(JC,JP,JY,JX)); GSTGI=0._r8
  allocate(DGSTGI(JC,JP,JY,JX));DGSTGI=0._r8
  allocate(DGSTGF(JC,JP,JY,JX));DGSTGF=0._r8
  allocate(PSTG(JC,JP,JY,JX));  PSTG=0._r8
  allocate(PSTGI(JC,JP,JY,JX)); PSTGI=0._r8
  allocate(GSTGF(JC,JP,JY,JX)); GSTGF=0._r8
  allocate(PSTGF(JC,JP,JY,JX)); PSTGF=0._r8
  allocate(TGSTGI(JC,JP,JY,JX));TGSTGI=0._r8
  allocate(TGSTGF(JC,JP,JY,JX));TGSTGF=0._r8
  allocate(XRNI(JP,JY,JX));     XRNI=0._r8
  allocate(SNL1(JP,JY,JX));     SNL1=0._r8
  allocate(FNOD(JP,JY,JX));     FNOD=0._r8
  allocate(NNOD(JP,JY,JX));     NNOD=0
  allocate(PSILZ(JP,JY,JX));    PSILZ=0._r8
  allocate(CFX(JP,JY,JX));      CFX=0._r8
  allocate(CF(JP,JY,JX));       CF=0._r8
  allocate(IDTHP(JP,JY,JX));    IDTHP=0
  allocate(STMX(JP,JY,JX));     STMX=0._r8
  allocate(SDMX(JP,JY,JX));     SDMX=0._r8
  allocate(GRMX(JP,JY,JX));     GRMX=0._r8
  allocate(XTLI(JP,JY,JX));     XTLI=0._r8
  allocate(GRDM(JP,JY,JX));     GRDM=0._r8
  allocate(GFILL(JP,JY,JX));    GFILL=0._r8
  allocate(FLG4(JC,JP,JY,JX));  FLG4=0._r8
  allocate(ATRP(JC,JP,JY,JX));  ATRP=0._r8
  allocate(FLGZ(JC,JP,JY,JX));  FLGZ=0._r8
  allocate(IDAY(jpstgs,JC,JP,JY,JX));IDAY=0
  allocate(CTC(JP,JY,JX));      CTC=0._r8
  allocate(HTC(JP,JY,JX));      HTC=0._r8
  allocate(SSTX(JP,JY,JX));     SSTX=0._r8
  allocate(XDL(JP,JY,JX));      XDL=0._r8
  allocate(XPPD(JP,JY,JX));     XPPD=0._r8
  allocate(CFI(JP,JY,JX));      CFI=0._r8
  allocate(VRNX(JC,JP,JY,JX));  VRNX=0._r8
  allocate(OFFST(JP,JY,JX));    OFFST=0._r8
  end subroutine InitPlantTraits

!----------------------------------------------------------------------
  subroutine DestructPlantTraits
  use abortutils, only : destroy
  implicit none
  call destroy(ARSTK)
  call destroy(ARLFP)
  call destroy(ARLFS)
  call destroy(ARLFX)
  call destroy(ARSTP)
  call destroy(ZC)
  call destroy(ARSTX)
  call destroy(ARLFT)
  call destroy(ARSTT)
  call destroy(ARLFC)
  call destroy(ARSTC)
  call destroy(ARLSS)
  call destroy(NG)
  call destroy(SDPTHI)
  call destroy(SDPTH)
  call destroy(SDVL)
  call destroy(SDLG)
  call destroy(SDAR)
  call destroy(HTCTL)
  call destroy(ZT)
  call destroy(ZL)
  call destroy(ANGBR)
  call destroy(ANGSH)
  call destroy(RAZ)
  call destroy(HTSTZ)
  call destroy(ARLF)
  call destroy(HTSHE)
  call destroy(HTNODE)
  call destroy(ARLFB)
  call destroy(ARLFZ)
  call destroy(HTSHEX)
  call destroy(GRNOB)
  call destroy(GRNXB)
  call destroy(GRNO)
  call destroy(PP)
  call destroy(HTNODX)
  call destroy(CNLF)
  call destroy(CPLF)
  call destroy(CNSHE)
  call destroy(CNSTK)
  call destroy(CNRSV)
  call destroy(CNHSK)
  call destroy(CNEAR)
  call destroy(CNGR)
  call destroy(CNND)
  call destroy(CPSHE)
  call destroy(CPSTK)
  call destroy(CPRSV)
  call destroy(CPHSK)
  call destroy(CPEAR)
  call destroy(CPGR)
  call destroy(CPND)
  call destroy(CNWS)
  call destroy(CPWS)
  call destroy(OSMO)
  call destroy(TCX)
  call destroy(ZTYPI)
  call destroy(ZTYP)
  call destroy(GROUP)
  call destroy(GROUPI)
  call destroy(GROUPX)
  call destroy(PPI)
  call destroy(WTSTDI)
  call destroy(PPX)
  call destroy(IFLGT)
  call destroy(PPT)
  call destroy(PPZ)
  call destroy(WSTR)
  call destroy(OSTR)
  call destroy(TFN3)
  call destroy(TCG)
  call destroy(TKG)
  call destroy(DMSHE)
  call destroy(DMSTK)
  call destroy(DMRSV)
  call destroy(DMHSK)
  call destroy(DMEAR)
  call destroy(DMGR)
  call destroy(DMND)
  call destroy(DMLF)
  call destroy(VRNY)
  call destroy(VRNZ)
  call destroy(VSTG)
  call destroy(VSTGX)
  call destroy(VRNS)
  call destroy(VRNF)
  call destroy(KLEAF)
  call destroy(KVSTGN)
  call destroy(KLEAFX)
  call destroy(KVSTG)
  call destroy(XRLA)
  call destroy(WDLF)
  call destroy(SLA1)
  call destroy(TCZ)
  call destroy(SSL1)
  call destroy(VRNL)
  call destroy(NBR)
  call destroy(NBT)
  call destroy(NBTB)
  call destroy(NB1)
  call destroy(IFLGR)
  call destroy(IFLGQ)
  call destroy(IFLGG)
  call destroy(IFLGP)
  call destroy(IFLGA)
  call destroy(IFLGE)
  call destroy(IFLGF)
  call destroy(IDTHB)
  call destroy(PB)
  call destroy(GSTGI)
  call destroy(DGSTGI)
  call destroy(DGSTGF)
  call destroy(PSTG)
  call destroy(PSTGI)
  call destroy(GSTGF)
  call destroy(PSTGF)
  call destroy(TGSTGI)
  call destroy(TGSTGF)
  call destroy(XRNI)
  call destroy(SNL1)
  call destroy(FNOD)
  call destroy(NNOD)
  call destroy(PSILZ)
  call destroy(CFX)
  call destroy(CF)
  call destroy(IDTHP)
  call destroy(STMX)
  call destroy(SDMX)
  call destroy(GRMX)
  call destroy(XTLI)
  call destroy(GRDM)
  call destroy(GFILL)
  call destroy(FLG4)
  call destroy(ATRP)
  call destroy(FLGZ)
  call destroy(IDAY)
  call destroy(CTC)
  call destroy(HTC)
  call destroy(SSTX)
  call destroy(XDL)
  call destroy(XPPD)
  call destroy(CFI)
  call destroy(VRNX)
  call destroy(OFFST)
  end subroutine DestructPlantTraits

end module PlantTraitDataType
