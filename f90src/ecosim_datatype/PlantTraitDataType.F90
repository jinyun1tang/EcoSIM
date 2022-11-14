module PlantTraitDataType

!
!!
! data types of plant trait characteristics that cannot be grouped into canopy
! or roots
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use ElmIDMod
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__


!allocation parameter
  real(r8) :: FVRN(0:5)              !allocation parameter
  REAL(R8),target,allocatable :: FWODBE(:,:)             !woody element allocation
  real(r8),target,allocatable :: FWODLE(:,:)             !element allocation for leaf
  real(r8),target,allocatable :: FWODRE(:,:)
  real(r8),target,allocatable :: FWOODE(:,:)             !woody C allocation
  real(r8),target,allocatable ::  ARSTK(:,:,:,:,:)                   !stem layer area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFP(:,:,:)                       !plant leaf area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFS(:,:,:)                       !plant leaf area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFX(:,:)                         !total canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  ARSTP(:,:,:)                       !plant stem area, [m2 d-2]
  real(r8),target,allocatable ::  ZC(:,:,:)                          !canopy height, [m]
  real(r8),target,allocatable ::  ARSTX(:,:)                         !total canopy stem area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFT(:,:,:)                       !total leaf area, [m2 d-2]
  real(r8),target,allocatable ::  ARSTT(:,:,:)                       !total stem area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFC(:,:)                         !total canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  ARSTC(:,:)                         !total canopy stem area, [m2 d-2]
  real(r8),target,allocatable ::  ARLSS(:,:)                         !stalk area of combined,each PFT canopy
  integer ,target,allocatable ::  NG(:,:,:)                           !soil layer at planting depth, [-]
  real(r8),target,allocatable ::  SDPTHI(:,:,:)                      !planting depth, [m]
  real(r8),target,allocatable ::  SDPTH(:,:,:)                       !seeding depth, [m]
  real(r8),target,allocatable ::  SDVL(:,:,:)                        !seed volume, [m3 ]
  real(r8),target,allocatable ::  SDLG(:,:,:)                        !seed length, [m]
  real(r8),target,allocatable ::  SDAR(:,:,:)                        !seed surface area, [m2]
  real(r8),target,allocatable ::  HTCTL(:,:,:)                       !cotyledon height, [m]
  real(r8),target,allocatable ::  ZT(:,:)                            !canopy height , [m]
  real(r8),target,allocatable ::  ZL(:,:,:)                          !canopy layer height , [m]
  real(r8),target,allocatable ::  ANGBR(:,:,:)                       !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  ANGSH(:,:,:)                       !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  RAZ(:,:,:)                         !canopy roughness height, [m]
  real(r8),target,allocatable ::  HTSTZ(:,:,:)                       !canopy height, [m]
  real(r8),target,allocatable ::  ARLF(:,:,:,:,:)                    !leaf area, [m2 d-2]
  real(r8),target,allocatable ::  HTSHE(:,:,:,:,:)                   !sheath height, [m]
  real(r8),target,allocatable ::  HTNODE(:,:,:,:,:)                  !internode height, [m]
  real(r8),target,allocatable ::  ARLFB(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFZ(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  HTSHEX(:,:,:,:)                    !branch height, [m]
  real(r8),target,allocatable ::  GRNOB(:,:,:,:)                     !branch grain number, [d-2]
  real(r8),target,allocatable ::  GRNXB(:,:,:,:)                     !branch potential grain number, [d-2]
  real(r8),target,allocatable ::  GRNO(:,:,:)                        !canopy grain number, [d-2]
  real(r8),target,allocatable ::  PP(:,:,:)                          !plant population, [d-2]
  real(r8),target,allocatable ::  HTNODX(:,:,:,:,:)                  !internode height, [m]
  real(r8),target,allocatable ::  CNLF(:,:,:)                        !maximum leaf N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPLF(:,:,:)                        !maximum leaf P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNSHE(:,:,:)                       !sheath N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNSTK(:,:,:)                       !stalk N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNRSV(:,:,:)                       !reserve N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNHSK(:,:,:)                       !husk N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNEAR(:,:,:)                       !ear N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNGR(:,:,:)                        !grain N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNND(:,:,:)                        !nodule N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPSHE(:,:,:)                       !sheath P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPSTK(:,:,:)                       !stalk P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPRSV(:,:,:)                       !reserve P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPHSK(:,:,:)                       !husk P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPEAR(:,:,:)                       !ear P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPGR(:,:,:)                        !grain P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPND(:,:,:)                        !nodule P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNWS(:,:,:)                        !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  CPWS(:,:,:)                        !C:P ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  OSMO(:,:,:)                        !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8),target,allocatable ::  TCX(:,:,:)                         !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8),target,allocatable ::  ZTYPI(:,:,:)                       !initial plant thermal adaptation zone, [-]
  real(r8),target,allocatable ::  ZTYP(:,:,:)                        !plant thermal adaptation zone, [-]
  real(r8),target,allocatable ::  GROUP(:,:,:,:)                     !plant maturity group, [-]
  real(r8),target,allocatable ::  GROUPI(:,:,:)                      !acclimated plant maturity group, [-]
  real(r8),target,allocatable ::  GROUPX(:,:,:)                      !initial plant maturity group, [-]
  real(r8),target,allocatable ::  PPI(:,:,:)                         !initial plant population, [m-2]
  real(r8),target,allocatable ::  WTSTDI(:,:,:)                      !initial standing dead C, [g C m-2]
  real(r8),target,allocatable ::  PPX(:,:,:)                         !plant population, [m-2]
  integer,target,allocatable ::  IFLGT(:,:)                          !number of active PFT
  real(r8),target,allocatable ::  PPT(:,:)                           !total plant population, [d-2]
  real(r8),target,allocatable ::  PPZ(:,:,:)                         !plant population at seeding, [m-2]
  real(r8),target,allocatable ::  WSTR(:,:,:)                        !canopy plant water stress indicator, number of hours PSILT < PSILY, []
  real(r8),target,allocatable ::  OSTR(:,:,:)                        !plant O2 stress indicator, []
  real(r8),target,allocatable ::  TFN3(:,:,:)                        !canopy temperature growth function, [-]
  real(r8),target,allocatable ::  TCG(:,:,:)                         !canopy growth temperature, [oC]
  real(r8),target,allocatable ::  TKG(:,:,:)                         !canopy growth temperature, [K]
  real(r8),target,allocatable ::  DMSHE(:,:,:)                       !sheath growth yield, [g g-1]
  real(r8),target,allocatable ::  DMSTK(:,:,:)                       !stalk growth yield, [g g-1]
  real(r8),target,allocatable ::  DMRSV(:,:,:)                       !reserve growth yield, [g g-1]
  real(r8),target,allocatable ::  DMHSK(:,:,:)                       !husk growth yield, [g g-1]
  real(r8),target,allocatable ::  DMEAR(:,:,:)                       !ear growth yield, [g g-1]
  real(r8),target,allocatable ::  DMGR(:,:,:)                        !grain growth yield, [g g-1]
  real(r8),target,allocatable ::  DMND(:,:,:)                        !nodule growth yield, [g g-1]
  real(r8),target,allocatable ::  DMLF(:,:,:)                        !leaf growth yield, [g g-1]
  real(r8),target,allocatable ::  VRNY(:,:,:,:)                      !initial heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  VRNZ(:,:,:,:)                      !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8),target,allocatable ::  VSTG(:,:,:,:)                      !leaf number, [-]
  real(r8),target,allocatable ::  VSTGX(:,:,:,:)                     !leaf number at floral initiation, [-]
  real(r8),target,allocatable ::  VRNS(:,:,:,:)                      !heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  VRNF(:,:,:,:)                      !cold requirement for autumn leafoff/hardening, [h]
  integer,target,allocatable ::  KLEAF(:,:,:,:)                      !leaf number, [-]
  integer,target,allocatable ::  KVSTGN(:,:,:,:)                     !leaf growth stage counter, [-]
  integer,target,allocatable ::  KLEAFX(:,:,:,:)                     !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,target,allocatable ::  KVSTG(:,:,:,:)                      !leaf growth stage counter, [-]
  real(r8),target,allocatable ::  XRLA(:,:,:)                        !rate of leaf initiation, [h-1 at 25 oC]
  real(r8),target,allocatable ::  WDLF(:,:,:)                        !leaf length:width ratio, [-]
  real(r8),target,allocatable ::  SLA1(:,:,:)                        !leaf area:mass during growth, [m2 g-1]
  real(r8),target,allocatable ::  TCZ(:,:,:)                         !threshold temperature for spring leafout/dehardening, [oC]
  real(r8),target,allocatable ::  SSL1(:,:,:)                        !petiole length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  VRNL(:,:,:,:)                      !hours above threshold temperature required for spring leafout/dehardening, [-]
  integer,target,allocatable ::  NBR(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  NBT(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  NBTB(:,:,:,:)                       !branch number, [-]
  integer,target,allocatable ::  NB1(:,:,:)                          !number of main branch, [-]
  integer,target,allocatable ::  IFLGR(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IFLGQ(:,:,:,:)                      !branch phenology flag, [h]
  integer,target,allocatable ::  IFLGG(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IFLGP(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IFLGA(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IFLGE(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IFLGF(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IDTHB(:,:,:,:)                      !flag to detect branch death , [-]
  real(r8),target,allocatable ::  PB(:,:,:)                          !branch nonstructural C content required for new branch, [g g-1]
  real(r8),target,allocatable ::  GSTGI(:,:,:,:)                     !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  DGSTGI(:,:,:,:)                    !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8),target,allocatable ::  DGSTGF(:,:,:,:)                    !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8),target,allocatable ::  PSTG(:,:,:,:)                      !node number, [-]
  real(r8),target,allocatable ::  PSTGI(:,:,:,:)                     !node number at floral initiation, [-]
  real(r8),target,allocatable ::  GSTGF(:,:,:,:)                     !normalized node number during reproductive growth stages, [-]
  real(r8),target,allocatable ::  PSTGF(:,:,:,:)                     !node number at anthesis, [-]
  real(r8),target,allocatable ::  TGSTGI(:,:,:,:)                    !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  TGSTGF(:,:,:,:)                    !normalized node number during reproductive growth stages , [-]
  real(r8),target,allocatable ::  XRNI(:,:,:)                        !rate of node initiation, [h-1 at 25 oC]
  real(r8),target,allocatable ::  SNL1(:,:,:)                        !internode length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  FNOD(:,:,:)                        !parameter for allocation of growth to nodes, [-]
  integer,target,allocatable ::  NNOD(:,:,:)                         !number of concurrently growing nodes
  real(r8),target,allocatable ::  PSILZ(:,:,:)                       !minimum daily canopy water potential, [MPa]
  real(r8),target,allocatable ::  CFX(:,:,:)                         !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8),target,allocatable ::  CF(:,:,:)                          !clumping factor for self-shading in canopy layer, [-]
  integer,target,allocatable ::  IDTHP(:,:,:)                        !flag to detect canopy death , [-]
  real(r8),target,allocatable ::  STMX(:,:,:)                        !maximum grain node number per branch, [-]
  real(r8),target,allocatable ::  SDMX(:,:,:)                        !maximum grain number per node , [-]
  real(r8),target,allocatable ::  GRMX(:,:,:)                        !maximum grain size   , [g]
  real(r8),target,allocatable ::  XTLI(:,:,:)                        !number of nodes in seed, [-]
  real(r8),target,allocatable ::  GRDM(:,:,:)                        !grain size at seeding, [g]
  real(r8),target,allocatable ::  GFILL(:,:,:)                       !maximum rate of fill per grain, [g h-1]
  real(r8),target,allocatable ::  FLG4(:,:,:,:)                      !flag to detect physiological maturity from  grain fill , [-]
  real(r8),target,allocatable ::  ATRP(:,:,:,:)                      !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  FLGZ(:,:,:,:)                      !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  integer,target,allocatable ::  IDAY(:,:,:,:,:)                     !plant growth stage, [-]
  real(r8),target,allocatable ::  CTC(:,:,:)                         !temperature below which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  HTC(:,:,:)                         !temperature above which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  SSTX(:,:,:)                        !sensitivity to HTC (seeds oC-1 above HTC)
  real(r8),target,allocatable ::  XDL(:,:,:)                         !critical daylength for phenological progress, [h]
  real(r8),target,allocatable ::  XPPD(:,:,:)                        !difference between current and critical daylengths used to calculate  phenological progress, [h]
  real(r8),target,allocatable ::  CFI(:,:,:)                         !initial clumping factor for self-shading in canopy layer, [-]
  real(r8),target,allocatable ::  VRNX(:,:,:,:)                      !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8),target,allocatable ::  OFFST(:,:,:)                       !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
!----------------------------------------------------------------------

contains
  subroutine InitPlantTraits(n_pltlitrk)

  implicit none
  integer, intent(in) :: n_pltlitrk

  FVRN =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)
  allocate(FWODLE(npelms,1:n_pltlitrk));  FWODLE=0._r8
  allocate(FWODBE(npelms,1:n_pltlitrk));  FWODBE=0._r8
  allocate(FWODRE(npelms,1:n_pltlitrk));  FWODRE=0._r8         !
  allocate(FWOODE(npelms,1:n_pltlitrk));  FWOODE=0._r8         !woody element allocation
  allocate(ARSTK(JC,JBR,JP,JY,JX));ARSTK=0._r8
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
  allocate(ARLF(0:JNODS,JBR,JP,JY,JX));ARLF=0._r8
  allocate(HTSHE(0:JNODS,JBR,JP,JY,JX));HTSHE=0._r8
  allocate(HTNODE(0:JNODS,JBR,JP,JY,JX));HTNODE=0._r8
  allocate(ARLFB(JBR,JP,JY,JX)); ARLFB=0._r8
  allocate(ARLFZ(JBR,JP,JY,JX)); ARLFZ=0._r8
  allocate(HTSHEX(JBR,JP,JY,JX));HTSHEX=0._r8
  allocate(GRNOB(JBR,JP,JY,JX)); GRNOB=0._r8
  allocate(GRNXB(JBR,JP,JY,JX)); GRNXB=0._r8
  allocate(GRNO(JP,JY,JX));     GRNO=0._r8
  allocate(PP(JP,JY,JX));       PP=0._r8
  allocate(HTNODX(0:JNODS,JBR,JP,JY,JX));HTNODX=0._r8
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
  allocate(GROUP(JBR,JP,JY,JX)); GROUP=0._r8
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
  allocate(VRNY(JBR,JP,JY,JX));  VRNY=0._r8
  allocate(VRNZ(JBR,JP,JY,JX));  VRNZ=0._r8
  allocate(VSTG(JBR,JP,JY,JX));  VSTG=0._r8
  allocate(VSTGX(JBR,JP,JY,JX)); VSTGX=0._r8
  allocate(VRNS(JBR,JP,JY,JX));  VRNS=0._r8
  allocate(VRNF(JBR,JP,JY,JX));  VRNF=0._r8
  allocate(KLEAF(JBR,JP,JY,JX)); KLEAF=0
  allocate(KVSTGN(JBR,JP,JY,JX));KVSTGN=0
  allocate(KLEAFX(JBR,JP,JY,JX));KLEAFX=0
  allocate(KVSTG(JBR,JP,JY,JX)); KVSTG=0
  allocate(XRLA(JP,JY,JX));     XRLA=0._r8
  allocate(WDLF(JP,JY,JX));     WDLF=0._r8
  allocate(SLA1(JP,JY,JX));     SLA1=0._r8
  allocate(TCZ(JP,JY,JX));      TCZ=0._r8
  allocate(SSL1(JP,JY,JX));     SSL1=0._r8
  allocate(VRNL(JC,JP,JY,JX));  VRNL=0._r8
  allocate(NBR(JP,JY,JX));      NBR=0
  allocate(NBT(JP,JY,JX));      NBT=0
  allocate(NBTB(JBR,JP,JY,JX));  NBTB=0
  allocate(NB1(JP,JY,JX));      NB1=0
  allocate(IFLGR(JBR,JP,JY,JX)); IFLGR=0
  allocate(IFLGQ(JBR,JP,JY,JX)); IFLGQ=0
  allocate(IFLGG(JBR,JP,JY,JX)); IFLGG=0
  allocate(IFLGP(JBR,JP,JY,JX)); IFLGP=0
  allocate(IFLGA(JBR,JP,JY,JX)); IFLGA=0
  allocate(IFLGE(JBR,JP,JY,JX)); IFLGE=0
  allocate(IFLGF(JBR,JP,JY,JX)); IFLGF=0
  allocate(IDTHB(JBR,JP,JY,JX)); IDTHB=0
  allocate(PB(JP,JY,JX));       PB=0._r8
  allocate(GSTGI(JBR,JP,JY,JX)); GSTGI=0._r8
  allocate(DGSTGI(JBR,JP,JY,JX));DGSTGI=0._r8
  allocate(DGSTGF(JBR,JP,JY,JX));DGSTGF=0._r8
  allocate(PSTG(JBR,JP,JY,JX));  PSTG=0._r8
  allocate(PSTGI(JBR,JP,JY,JX)); PSTGI=0._r8
  allocate(GSTGF(JBR,JP,JY,JX)); GSTGF=0._r8
  allocate(PSTGF(JBR,JP,JY,JX)); PSTGF=0._r8
  allocate(TGSTGI(JBR,JP,JY,JX));TGSTGI=0._r8
  allocate(TGSTGF(JBR,JP,JY,JX));TGSTGF=0._r8
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
  allocate(FLG4(JBR,JP,JY,JX));  FLG4=0._r8
  allocate(ATRP(JBR,JP,JY,JX));  ATRP=0._r8
  allocate(FLGZ(JBR,JP,JY,JX));  FLGZ=0._r8
  allocate(IDAY(jpstgs,JBR,JP,JY,JX));IDAY=0
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

  call destroy(FWODBE)
  call destroy(FWODRE)
  call destroy(FWOODE)
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
