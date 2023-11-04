module PlantTraitDataType

!
!!
! data types of plant trait characteristics that cannot be grouped into canopy
! or roots
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__


!allocation parameter
  real(r8) :: FVRN(0:5)              !allocation parameter
  REAL(R8),target,allocatable :: FWODBE(:,:)             !woody element allocation
  real(r8),target,allocatable :: FWODLE(:,:)             !element allocation for leaf
  real(r8),target,allocatable :: FWODRE(:,:)
  real(r8),target,allocatable :: FWOODE(:,:)             !woody C allocation
  real(r8),target,allocatable ::  CanopyBranchStemApft_lyr(:,:,:,:,:)              !stem layer area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLeafA_pft(:,:,:)                       !plant leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyArea_pft(:,:,:)                       !plant canopy leaf+stem/stalk area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFX(:,:)                         !total canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyStemA_pft(:,:,:)                      !plant stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyHeight(:,:,:)                          !canopy height, [m]
  real(r8),target,allocatable ::  ARSTX(:,:)                         !total canopy stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLAgrid_lyr(:,:,:)                       !total leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyStemA_lyr(:,:,:)                       !total stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLA_grd(:,:)                         !grid level plant canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  StemAreag(:,:)                     !total canopy stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyArea_grid(:,:)                         !canopy area of combined over the grid [m2 d-2]
  integer ,target,allocatable ::  NGTopRootLayer(:,:,:)                           !soil layer at planting depth, [-]
  real(r8),target,allocatable ::  PlantinDepth(:,:,:)                      !planting depth, [m]
  real(r8),target,allocatable ::  SeedinDepth(:,:,:)                       !seeding depth, [m]
  real(r8),target,allocatable ::  SeedVolume(:,:,:)                        !seed volume, [m3 ]
  real(r8),target,allocatable ::  SeedLength(:,:,:)                        !seed length, [m]
  real(r8),target,allocatable ::  SeedArea(:,:,:)                        !seed surface area, [m2]
  real(r8),target,allocatable ::  HypoctoylHeight(:,:,:)                       !cotyledon height, [m]
  real(r8),target,allocatable ::  GridMaxCanopyHeight(:,:)                            !canopy height , [m]
  real(r8),target,allocatable ::  CanopyHeightz(:,:,:)                          !canopy layer height , [m]
  real(r8),target,allocatable ::  ANGBR(:,:,:)                       !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  ANGSH(:,:,:)                       !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  RAZ(:,:,:)                         !canopy roughness height, [m]
  real(r8),target,allocatable ::  CanPHeight4WatUptake(:,:,:)                       !canopy height, [m]
  real(r8),target,allocatable ::  ARLF(:,:,:,:,:)                    !leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanPSheathHeight(:,:,:,:,:)                   !sheath height, [m]
  real(r8),target,allocatable ::  HTNODE(:,:,:,:,:)                  !internode height, [m]
  real(r8),target,allocatable ::  CanopyBranchLeafA_pft(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFZ(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanPBranchHeight(:,:,:,:)                    !branch height, [m]
  real(r8),target,allocatable ::  GRNOB(:,:,:,:)                     !branch grain number, [d-2]
  real(r8),target,allocatable ::  GRNXB(:,:,:,:)                     !branch potential grain number, [d-2]
  real(r8),target,allocatable ::  GRNO(:,:,:)                        !canopy grain number, [d-2]
  real(r8),target,allocatable ::  pftPlantPopulation(:,:,:)             !plant population, [d-2]
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
  integer,target,allocatable ::  NumActivePlants(:,:)                          !number of active PFT
  real(r8),target,allocatable ::  PPT(:,:)                           !total plant population, [d-2]
  real(r8),target,allocatable ::  PPZ(:,:,:)                         !plant population at seeding, [m-2]
  real(r8),target,allocatable ::  WSTR(:,:,:)                        !canopy plant water stress indicator, number of hours PSILT < PSILY, []
  real(r8),target,allocatable ::  PlantO2Stress(:,:,:)                        !plant O2 stress indicator, []
  real(r8),target,allocatable ::  fTgrowCanP(:,:,:)                  !canopy temperature growth function, [-]
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
  real(r8),target,allocatable ::  Hours4LeafOff(:,:,:,:)                      !cold requirement for autumn leafoff/hardening, [h]
  integer,target,allocatable ::  KLEAF(:,:,:,:)                      !leaf number, [-]
  integer,target,allocatable ::  KVSTGN(:,:,:,:)                     !leaf growth stage counter, [-]
  integer,target,allocatable ::  KLEAFX(:,:,:,:)                     !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,target,allocatable ::  KVSTG(:,:,:,:)                      !leaf growth stage counter, [-]
  real(r8),target,allocatable ::  XRLA(:,:,:)                        !rate of leaf initiation, [h-1 at 25 oC]
  real(r8),target,allocatable ::  WDLF(:,:,:)                        !leaf length:width ratio, [-]
  real(r8),target,allocatable ::  SLA1(:,:,:)                        !leaf area:mass during growth, [m2 g-1]
  real(r8),target,allocatable ::  TCZ(:,:,:)                         !threshold temperature for spring leafout/dehardening, [oC]
  real(r8),target,allocatable ::  SSL1(:,:,:)                        !petiole length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  HourThreshold4LeafOut(:,:,:,:)                      !hours above threshold temperature required for spring leafout/dehardening, [-]
  integer,target,allocatable ::  NumOfBranches_pft(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  NBT(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  BranchNumber_brchpft(:,:,:,:)                       !branch number, [-]
  integer,target,allocatable ::  NB1(:,:,:)                          !number of main branch, [-]
  integer,target,allocatable ::  IFLGR(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IFLGQ(:,:,:,:)                      !branch phenology flag, [h]
  integer,target,allocatable ::  doPlantSenescence(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doPlantRemobilization(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doInitLeafOut(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doPlantLeafOut(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doPlantLeaveOff(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  iPlantBranchState(:,:,:,:)                      !flag to detect branch death , [-]
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
  real(r8),target,allocatable ::  PSICanPDailyMin(:,:,:)                       !minimum daily canopy water potential, [MPa]
  real(r8),target,allocatable ::  ClumpFactort(:,:,:)                         !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8),target,allocatable ::  ClumpFactor(:,:,:)                          !clumping factor for self-shading in canopy layer, [-]
  integer,target,allocatable ::  iPlantShootState(:,:,:)                        !flag to detect canopy death , [-]
  real(r8),target,allocatable ::  STMX(:,:,:)                        !maximum grain node number per branch, [-]
  real(r8),target,allocatable ::  SDMX(:,:,:)                        !maximum grain number per node , [-]
  real(r8),target,allocatable ::  MaxSeedCMass(:,:,:)                        !maximum grain size   , [g]
  real(r8),target,allocatable ::  XTLI(:,:,:)                        !number of nodes in seed, [-]
  real(r8),target,allocatable ::  SeedCMass(:,:,:)                        !grain size at seeding, [g]
  real(r8),target,allocatable ::  GFILL(:,:,:)                       !maximum rate of fill per grain, [g h-1]
  real(r8),target,allocatable ::  FLG4(:,:,:,:)                      !flag to detect physiological maturity from  grain fill , [-]
  real(r8),target,allocatable ::  HourCounter4LeafOut_brch(:,:,:,:)                      !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  FLGZ(:,:,:,:)                      !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  integer,target,allocatable ::  iPlantCalendar(:,:,:,:,:)                     !plant growth stage, [-]
  real(r8),target,allocatable ::  CTC(:,:,:)                         !temperature below which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  HTC(:,:,:)                         !temperature above which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  SSTX(:,:,:)                        !sensitivity to HTC (seeds oC-1 above HTC)
  real(r8),target,allocatable ::  XDL(:,:,:)                         !critical daylength for phenological progress, [h]
  real(r8),target,allocatable ::  XPPD(:,:,:)                        !difference between current and critical daylengths used to calculate  phenological progress, [h]
  real(r8),target,allocatable ::  ClumpFactort0(:,:,:)                         !initial clumping factor for self-shading in canopy layer, [-]
  real(r8),target,allocatable ::  HourThreshold4LeafOff(:,:,:,:)                      !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8),target,allocatable ::  OFFST(:,:,:)                       !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
!----------------------------------------------------------------------

contains
  subroutine InitPlantTraits(NumOfPlantLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs

  FVRN =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)
  allocate(FWODLE(NumOfPlantChemElements,1:NumOfPlantLitrCmplxs));  FWODLE=0._r8
  allocate(FWODBE(NumOfPlantChemElements,1:NumOfPlantLitrCmplxs));  FWODBE=0._r8
  allocate(FWODRE(NumOfPlantChemElements,1:NumOfPlantLitrCmplxs));  FWODRE=0._r8         !
  allocate(FWOODE(NumOfPlantChemElements,1:NumOfPlantLitrCmplxs));  FWOODE=0._r8         !woody element allocation
  allocate(CanopyBranchStemApft_lyr(JC,JBR,JP,JY,JX));CanopyBranchStemApft_lyr=0._r8
  allocate(CanopyLeafA_pft(JP,JY,JX));    CanopyLeafA_pft=0._r8
  allocate(CanopyArea_pft(JP,JY,JX));    CanopyArea_pft=0._r8
  allocate(ARLFX(JY,JX));       ARLFX=0._r8
  allocate(CanopyStemA_pft(JP,JY,JX));    CanopyStemA_pft=0._r8
  allocate(CanopyHeight(JP,JY,JX));       CanopyHeight=0._r8
  allocate(ARSTX(JY,JX));       ARSTX=0._r8
  allocate(CanopyLAgrid_lyr(JC,JY,JX));    CanopyLAgrid_lyr=0._r8
  allocate(CanopyStemA_lyr(JC,JY,JX));    CanopyStemA_lyr=0._r8
  allocate(CanopyLA_grd(JY,JX));       CanopyLA_grd=0._r8
  allocate(StemAreag(JY,JX));       StemAreag=0._r8
  allocate(CanopyArea_grid(JY,JX));       CanopyArea_grid=0._r8
  allocate(NGTopRootLayer(JP,JY,JX));       NGTopRootLayer=0
  allocate(PlantinDepth(JP,JY,JX));   PlantinDepth=0._r8
  allocate(SeedinDepth(JP,JY,JX));    SeedinDepth=0._r8
  allocate(SeedVolume(JP,JY,JX));     SeedVolume=0._r8
  allocate(SeedLength(JP,JY,JX));     SeedLength=0._r8
  allocate(SeedArea(JP,JY,JX));     SeedArea=0._r8
  allocate(HypoctoylHeight(JP,JY,JX));    HypoctoylHeight=0._r8
  allocate(GridMaxCanopyHeight(JY,JX));          GridMaxCanopyHeight=0._r8
  allocate(CanopyHeightz(0:JC,JY,JX));     CanopyHeightz=0._r8
  allocate(ANGBR(JP,JY,JX));    ANGBR=0._r8
  allocate(ANGSH(JP,JY,JX));    ANGSH=0._r8
  allocate(RAZ(JP,JY,JX));      RAZ=0._r8
  allocate(CanPHeight4WatUptake(JP,JY,JX));    CanPHeight4WatUptake=0._r8
  allocate(ARLF(0:MaxCanopyNodes,JBR,JP,JY,JX));ARLF=0._r8
  allocate(CanPSheathHeight(0:MaxCanopyNodes,JBR,JP,JY,JX));CanPSheathHeight=0._r8
  allocate(HTNODE(0:MaxCanopyNodes,JBR,JP,JY,JX));HTNODE=0._r8
  allocate(CanopyBranchLeafA_pft(JBR,JP,JY,JX)); CanopyBranchLeafA_pft=0._r8
  allocate(ARLFZ(JBR,JP,JY,JX)); ARLFZ=0._r8
  allocate(CanPBranchHeight(JBR,JP,JY,JX));CanPBranchHeight=0._r8
  allocate(GRNOB(JBR,JP,JY,JX)); GRNOB=0._r8
  allocate(GRNXB(JBR,JP,JY,JX)); GRNXB=0._r8
  allocate(GRNO(JP,JY,JX));     GRNO=0._r8
  allocate(pftPlantPopulation(JP,JY,JX));       pftPlantPopulation=0._r8
  allocate(HTNODX(0:MaxCanopyNodes,JBR,JP,JY,JX));HTNODX=0._r8
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
  allocate(NumActivePlants(JY,JX));       NumActivePlants=0
  allocate(PPT(JY,JX));         PPT=0._r8
  allocate(PPZ(JP,JY,JX));      PPZ=0._r8
  allocate(WSTR(JP,JY,JX));     WSTR=0._r8
  allocate(PlantO2Stress(JP,JY,JX));     PlantO2Stress=0._r8
  allocate(fTgrowCanP(JP,JY,JX));     fTgrowCanP=0._r8
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
  allocate(Hours4LeafOff(JBR,JP,JY,JX));  Hours4LeafOff=0._r8
  allocate(KLEAF(JBR,JP,JY,JX)); KLEAF=0
  allocate(KVSTGN(JBR,JP,JY,JX));KVSTGN=0
  allocate(KLEAFX(JBR,JP,JY,JX));KLEAFX=0
  allocate(KVSTG(JBR,JP,JY,JX)); KVSTG=0
  allocate(XRLA(JP,JY,JX));     XRLA=0._r8
  allocate(WDLF(JP,JY,JX));     WDLF=0._r8
  allocate(SLA1(JP,JY,JX));     SLA1=0._r8
  allocate(TCZ(JP,JY,JX));      TCZ=0._r8
  allocate(SSL1(JP,JY,JX));     SSL1=0._r8
  allocate(HourThreshold4LeafOut(JC,JP,JY,JX));  HourThreshold4LeafOut=0._r8
  allocate(NumOfBranches_pft(JP,JY,JX));      NumOfBranches_pft=0
  allocate(NBT(JP,JY,JX));      NBT=0
  allocate(BranchNumber_brchpft(JBR,JP,JY,JX));  BranchNumber_brchpft=0
  allocate(NB1(JP,JY,JX));      NB1=0
  allocate(IFLGR(JBR,JP,JY,JX)); IFLGR=0
  allocate(IFLGQ(JBR,JP,JY,JX)); IFLGQ=0
  allocate(doPlantSenescence(JBR,JP,JY,JX)); doPlantSenescence=0
  allocate(doPlantRemobilization(JBR,JP,JY,JX)); doPlantRemobilization=0
  allocate(doInitLeafOut(JBR,JP,JY,JX)); doInitLeafOut=0
  allocate(doPlantLeafOut(JBR,JP,JY,JX)); doPlantLeafOut=0
  allocate(doPlantLeaveOff(JBR,JP,JY,JX)); doPlantLeaveOff=0
  allocate(iPlantBranchState(JBR,JP,JY,JX)); iPlantBranchState=0
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
  allocate(PSICanPDailyMin(JP,JY,JX));    PSICanPDailyMin=0._r8
  allocate(ClumpFactort(JP,JY,JX));      ClumpFactort=0._r8
  allocate(ClumpFactor(JP,JY,JX));       ClumpFactor=0._r8
  allocate(iPlantShootState(JP,JY,JX));    iPlantShootState=0
  allocate(STMX(JP,JY,JX));     STMX=0._r8
  allocate(SDMX(JP,JY,JX));     SDMX=0._r8
  allocate(MaxSeedCMass(JP,JY,JX));     MaxSeedCMass=0._r8
  allocate(XTLI(JP,JY,JX));     XTLI=0._r8
  allocate(SeedCMass(JP,JY,JX));     SeedCMass=0._r8
  allocate(GFILL(JP,JY,JX));    GFILL=0._r8
  allocate(FLG4(JBR,JP,JY,JX));  FLG4=0._r8
  allocate(HourCounter4LeafOut_brch(JBR,JP,JY,JX));  HourCounter4LeafOut_brch=0._r8
  allocate(FLGZ(JBR,JP,JY,JX));  FLGZ=0._r8
  allocate(iPlantCalendar(NumGrowthStages,JBR,JP,JY,JX));iPlantCalendar=0
  allocate(CTC(JP,JY,JX));      CTC=0._r8
  allocate(HTC(JP,JY,JX));      HTC=0._r8
  allocate(SSTX(JP,JY,JX));     SSTX=0._r8
  allocate(XDL(JP,JY,JX));      XDL=0._r8
  allocate(XPPD(JP,JY,JX));     XPPD=0._r8
  allocate(ClumpFactort0(JP,JY,JX));      ClumpFactort0=0._r8
  allocate(HourThreshold4LeafOff(JC,JP,JY,JX));  HourThreshold4LeafOff=0._r8
  allocate(OFFST(JP,JY,JX));    OFFST=0._r8
  end subroutine InitPlantTraits

!----------------------------------------------------------------------
  subroutine DestructPlantTraits
  use abortutils, only : destroy
  implicit none

  call destroy(FWODBE)
  call destroy(FWODRE)
  call destroy(FWOODE)
  call destroy(CanopyBranchStemApft_lyr)
  call destroy(CanopyLeafA_pft)
  call destroy(CanopyArea_pft)
  call destroy(ARLFX)
  call destroy(CanopyStemA_pft)
  call destroy(CanopyHeight)
  call destroy(ARSTX)
  call destroy(CanopyLAgrid_lyr)
  call destroy(CanopyStemA_lyr)
  call destroy(CanopyLA_grd)
  call destroy(StemAreag)
  call destroy(CanopyArea_grid)
  call destroy(NGTopRootLayer)
  call destroy(PlantinDepth)
  call destroy(SeedinDepth)
  call destroy(SeedVolume)
  call destroy(SeedLength)
  call destroy(SeedArea)
  call destroy(HypoctoylHeight)
  call destroy(GridMaxCanopyHeight)
  call destroy(CanopyHeightz)
  call destroy(ANGBR)
  call destroy(ANGSH)
  call destroy(RAZ)
  call destroy(CanPHeight4WatUptake)
  call destroy(ARLF)
  call destroy(CanPSheathHeight)
  call destroy(HTNODE)
  call destroy(CanopyBranchLeafA_pft)
  call destroy(ARLFZ)
  call destroy(CanPBranchHeight)
  call destroy(GRNOB)
  call destroy(GRNXB)
  call destroy(GRNO)
  call destroy(pftPlantPopulation)
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
  call destroy(NumActivePlants)
  call destroy(PPT)
  call destroy(PPZ)
  call destroy(WSTR)
  call destroy(PlantO2Stress)
  call destroy(fTgrowCanP)
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
  call destroy(Hours4LeafOff)
  call destroy(KLEAF)
  call destroy(KVSTGN)
  call destroy(KLEAFX)
  call destroy(KVSTG)
  call destroy(XRLA)
  call destroy(WDLF)
  call destroy(SLA1)
  call destroy(TCZ)
  call destroy(SSL1)
  call destroy(HourThreshold4LeafOut)
  call destroy(NumOfBranches_pft)
  call destroy(NBT)
  call destroy(BranchNumber_brchpft)
  call destroy(NB1)
  call destroy(IFLGR)
  call destroy(IFLGQ)
  call destroy(doPlantSenescence)
  call destroy(doPlantRemobilization)
  call destroy(doInitLeafOut)
  call destroy(doPlantLeafOut)
  call destroy(doPlantLeaveOff)
  call destroy(iPlantBranchState)
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
  call destroy(PSICanPDailyMin)
  call destroy(ClumpFactort)
  call destroy(ClumpFactor)
  call destroy(iPlantShootState)
  call destroy(STMX)
  call destroy(SDMX)
  call destroy(MaxSeedCMass)
  call destroy(XTLI)
  call destroy(SeedCMass)
  call destroy(GFILL)
  call destroy(FLG4)
  call destroy(HourCounter4LeafOut_brch)
  call destroy(FLGZ)
  call destroy(iPlantCalendar)
  call destroy(CTC)
  call destroy(HTC)
  call destroy(SSTX)
  call destroy(XDL)
  call destroy(XPPD)
  call destroy(ClumpFactort0)
  call destroy(HourThreshold4LeafOff)
  call destroy(OFFST)
  end subroutine DestructPlantTraits

end module PlantTraitDataType
