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
  real(r8),target,allocatable ::  BranchAngle_pft(:,:,:)                       !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  PetioleAngle_pft(:,:,:)                       !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  SineBranchAngle_pft(:,:,:)                       !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  SinePetioleAngle_pft(:,:,:)                       !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  RAZ(:,:,:)                         !canopy roughness height, [m]
  real(r8),target,allocatable ::  CanPHeight4WatUptake(:,:,:)                       !canopy height, [m]
  real(r8),target,allocatable ::  LeafAreaNode_brch(:,:,:,:,:)                    !leaf area, [m2 d-2]
  real(r8),target,allocatable ::  PetioleLengthNode_brch(:,:,:,:,:)                   !sheath height, [m]
  real(r8),target,allocatable ::  InternodeHeightLive_brch(:,:,:,:,:)                  !internode height, [m]
  real(r8),target,allocatable ::  LeafAreaLive_brch(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  LeafAreaDying_brch(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanPBranchHeight(:,:,:,:)                    !branch height, [m]
  real(r8),target,allocatable ::  GRNOB(:,:,:,:)                     !branch grain number, [d-2]
  real(r8),target,allocatable ::  GRNXB(:,:,:,:)                     !branch potential grain number, [d-2]
  real(r8),target,allocatable ::  GRNO(:,:,:)                        !canopy grain number, [d-2]
  real(r8),target,allocatable ::  pftPlantPopulation(:,:,:)             !plant population, [d-2]
  real(r8),target,allocatable ::  InternodeHeightDying_brch(:,:,:,:,:)                  !internode height, [m]
  real(r8),target,allocatable ::  CNLF(:,:,:)                        !maximum leaf N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPLF(:,:,:)                        !maximum leaf P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNSHE(:,:,:)                       !sheath N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCStalk_pft(:,:,:)                       !stalk N:C ratio, [g g-1]
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
  real(r8),target,allocatable ::  iPlantInitThermoAdaptZone(:,:,:)                       !initial plant thermal adaptation zone, [-]
  real(r8),target,allocatable ::  iPlantThermoAdaptZone(:,:,:)                        !plant thermal adaptation zone, [-]
  real(r8),target,allocatable ::  MatureGroup_brch(:,:,:,:)                     !plant maturity group, [-]
  real(r8),target,allocatable ::  MatureGroup_pft(:,:,:)                      !acclimated plant maturity group, [-]
  real(r8),target,allocatable ::  GROUPX(:,:,:)                      !initial plant maturity group, [-]
  real(r8),target,allocatable ::  PPI(:,:,:)                         !initial plant population, [m-2]
  real(r8),target,allocatable ::  StandingDeadInitC_pft(:,:,:)                      !initial standing dead C, [g C m-2]
  real(r8),target,allocatable ::  PPX(:,:,:)                         !plant population, [m-2]
  integer,target,allocatable ::  NumActivePlants(:,:)                          !number of active PFT
  real(r8),target,allocatable ::  PPT(:,:)                           !total plant population, [d-2]
  real(r8),target,allocatable ::  PPZ(:,:,:)                         !plant population at seeding, [m-2]
  real(r8),target,allocatable ::  HoursCanopyPSITooLow(:,:,:)                        !canopy plant water stress indicator, number of hours PSILT < PSILY, []
  real(r8),target,allocatable ::  PlantO2Stress(:,:,:)                        !plant O2 stress indicator, []
  real(r8),target,allocatable ::  fTgrowCanP(:,:,:)                  !canopy temperature growth function, [-]
  real(r8),target,allocatable ::  TCG(:,:,:)                         !canopy growth temperature, [oC]
  real(r8),target,allocatable ::  TKG(:,:,:)                         !canopy growth temperature, [K]
  real(r8),target,allocatable ::  PetioleBiomGrowthYield(:,:,:)                       !sheath growth yield, [g g-1]
  real(r8),target,allocatable ::  StalkBiomGrowthYield(:,:,:)                       !stalk growth yield, [g g-1]
  real(r8),target,allocatable ::  ReserveBiomGrowthYield(:,:,:)                       !reserve growth yield, [g g-1]
  real(r8),target,allocatable ::  HuskBiomGrowthYield(:,:,:)                       !husk growth yield, [g g-1]
  real(r8),target,allocatable ::  EarBiomGrowthYield(:,:,:)                       !ear growth yield, [g g-1]
  real(r8),target,allocatable ::  GrainBiomGrowthYield(:,:,:)                        !grain growth yield, [g g-1]
  real(r8),target,allocatable ::  DMND(:,:,:)                        !nodule growth yield, [g g-1]
  real(r8),target,allocatable ::  LeafBiomGrowthYield(:,:,:)                        !leaf growth yield, [g g-1]
  real(r8),target,allocatable ::  Hours4LenthenPhotoPeriod(:,:,:,:)                      !initial heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  Hours4ShortenPhotoPeriod(:,:,:,:)                      !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8),target,allocatable ::  NumOfLeaves_brch(:,:,:,:)                      !leaf number, [-]
  real(r8),target,allocatable ::  VSTGX(:,:,:,:)                     !leaf number at floral initiation, [-]
  real(r8),target,allocatable ::  Hours4Leafout(:,:,:,:)                      !heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  Hours4LeafOff(:,:,:,:)                      !cold requirement for autumn leafoff/hardening, [h]
  integer,target,allocatable ::  KLEAF(:,:,:,:)                      !leaf number, [-]
  integer,target,allocatable ::  KVSTGN(:,:,:,:)                     !leaf growth stage counter, [-]
  integer,target,allocatable ::  KLEAFX(:,:,:,:)                     !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,target,allocatable ::  KLeafNodeNumber(:,:,:,:)                      !leaf growth stage counter, [-]
  real(r8),target,allocatable ::  RefLeafAppearRate(:,:,:)                        !rate of leaf initiation, [h-1 at 25 oC]
  real(r8),target,allocatable ::  WDLF(:,:,:)                        !leaf length:width ratio, [-]
  real(r8),target,allocatable ::  SLA1(:,:,:)                        !leaf area:mass during growth, [m2 g-1]
  real(r8),target,allocatable ::  TCelciusChill4Leaf(:,:,:)                         !threshold temperature for spring leafout/dehardening, [oC]
  real(r8),target,allocatable ::  SSL1(:,:,:)                        !petiole length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  HourThreshold4LeafOut(:,:,:,:)                      !hours above threshold temperature required for spring leafout/dehardening, [-]
  integer,target,allocatable ::  NumOfBranches_pft(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  BranchNumber_pft(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  BranchNumber_brch(:,:,:,:)                       !branch number, [-]
  integer,target,allocatable ::  NumOfMainBranch_pft(:,:,:)                          !number of main branch, [-]
  integer,target,allocatable ::  IFLGR(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  IFLGQ(:,:,:,:)                      !branch phenology flag, [h]
  integer,target,allocatable ::  doSenescence_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doRemobilization_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doInitLeafOut(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doPlantLeafOut(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doPlantLeaveOff(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  iPlantBranchState(:,:,:,:)                      !flag to detect branch death , [-]
  real(r8),target,allocatable ::  MinNonstructalC4InitBranch(:,:,:)                          !branch nonstructural C content required for new branch, [g g-1]
  real(r8),target,allocatable ::  GSTGI(:,:,:,:)                     !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  DGSTGI(:,:,:,:)                    !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8),target,allocatable ::  DGSTGF(:,:,:,:)                    !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8),target,allocatable ::  ShootNodeNumber(:,:,:,:)                      !node number, [-]
  real(r8),target,allocatable ::  NodeNumberToInitFloral(:,:,:,:)                     !node number at floral initiation, [-]
  real(r8),target,allocatable ::  GSTGF(:,:,:,:)                     !normalized node number during reproductive growth stages, [-]
  real(r8),target,allocatable ::  NodeNumberAtAnthesis(:,:,:,:)                     !node number at anthesis, [-]
  real(r8),target,allocatable ::  TGSTGI(:,:,:,:)                    !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  TGSTGF(:,:,:,:)                    !normalized node number during reproductive growth stages , [-]
  real(r8),target,allocatable ::  RefNodeInitRate(:,:,:)                        !rate of node initiation, [h-1 at 25 oC]
  real(r8),target,allocatable ::  SNL1(:,:,:)                        !internode length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  FNOD(:,:,:)                        !parameter for allocation of growth to nodes, [-]
  integer,target,allocatable ::  NumConCurrentGrowinNode(:,:,:)                         !number of concurrently growing nodes
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
  real(r8),target,allocatable ::  TCelciusChill4Seed(:,:,:)                         !temperature below which seed set is adversely affected, [oC]
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
  allocate(FWODLE(NumOfPlantChemElmnts,1:NumOfPlantLitrCmplxs));  FWODLE=0._r8
  allocate(FWODBE(NumOfPlantChemElmnts,1:NumOfPlantLitrCmplxs));  FWODBE=0._r8
  allocate(FWODRE(NumOfPlantChemElmnts,1:NumOfPlantLitrCmplxs));  FWODRE=0._r8         !
  allocate(FWOODE(NumOfPlantChemElmnts,1:NumOfPlantLitrCmplxs));  FWOODE=0._r8         !woody element allocation
  allocate(CanopyBranchStemApft_lyr(JC,MaxNumBranches,JP,JY,JX));CanopyBranchStemApft_lyr=0._r8
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
  allocate(BranchAngle_pft(JP,JY,JX));    BranchAngle_pft=0._r8
  allocate(PetioleAngle_pft(JP,JY,JX));    PetioleAngle_pft=0._r8
  allocate(SineBranchAngle_pft(JP,JY,JX));    SineBranchAngle_pft=0._r8
  allocate(SinePetioleAngle_pft(JP,JY,JX));    SinePetioleAngle_pft=0._r8
  allocate(RAZ(JP,JY,JX));      RAZ=0._r8
  allocate(CanPHeight4WatUptake(JP,JY,JX));    CanPHeight4WatUptake=0._r8
  allocate(LeafAreaNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafAreaNode_brch=0._r8
  allocate(PetioleLengthNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleLengthNode_brch=0._r8
  allocate(InternodeHeightLive_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));InternodeHeightLive_brch=0._r8
  allocate(LeafAreaLive_brch(MaxNumBranches,JP,JY,JX)); LeafAreaLive_brch=0._r8
  allocate(LeafAreaDying_brch(MaxNumBranches,JP,JY,JX)); LeafAreaDying_brch=0._r8
  allocate(CanPBranchHeight(MaxNumBranches,JP,JY,JX));CanPBranchHeight=0._r8
  allocate(GRNOB(MaxNumBranches,JP,JY,JX)); GRNOB=0._r8
  allocate(GRNXB(MaxNumBranches,JP,JY,JX)); GRNXB=0._r8
  allocate(GRNO(JP,JY,JX));     GRNO=0._r8
  allocate(pftPlantPopulation(JP,JY,JX));       pftPlantPopulation=0._r8
  allocate(InternodeHeightDying_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));InternodeHeightDying_brch=0._r8
  allocate(CNLF(JP,JY,JX));     CNLF=0._r8
  allocate(CPLF(JP,JY,JX));     CPLF=0._r8
  allocate(CNSHE(JP,JY,JX));    CNSHE=0._r8
  allocate(rNCStalk_pft(JP,JY,JX));    rNCStalk_pft=0._r8
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
  allocate(iPlantInitThermoAdaptZone(JP,JY,JX));    iPlantInitThermoAdaptZone=0._r8
  allocate(iPlantThermoAdaptZone(JP,JY,JX));     iPlantThermoAdaptZone=0._r8
  allocate(MatureGroup_brch(MaxNumBranches,JP,JY,JX)); MatureGroup_brch=0._r8
  allocate(MatureGroup_pft(JP,JY,JX));   MatureGroup_pft=0._r8
  allocate(GROUPX(JP,JY,JX));   GROUPX=0._r8
  allocate(PPI(JP,JY,JX));      PPI=0._r8
  allocate(StandingDeadInitC_pft(JP,JY,JX));   StandingDeadInitC_pft=0._r8
  allocate(PPX(JP,JY,JX));      PPX=0._r8
  allocate(NumActivePlants(JY,JX));       NumActivePlants=0
  allocate(PPT(JY,JX));         PPT=0._r8
  allocate(PPZ(JP,JY,JX));      PPZ=0._r8
  allocate(HoursCanopyPSITooLow(JP,JY,JX));     HoursCanopyPSITooLow=0._r8
  allocate(PlantO2Stress(JP,JY,JX));     PlantO2Stress=0._r8
  allocate(fTgrowCanP(JP,JY,JX));     fTgrowCanP=0._r8
  allocate(TCG(JP,JY,JX));      TCG=0._r8
  allocate(TKG(JP,JY,JX));      TKG=0._r8
  allocate(PetioleBiomGrowthYield(JP,JY,JX));    PetioleBiomGrowthYield=0._r8
  allocate(StalkBiomGrowthYield(JP,JY,JX));    StalkBiomGrowthYield=0._r8
  allocate(ReserveBiomGrowthYield(JP,JY,JX));    ReserveBiomGrowthYield=0._r8
  allocate(HuskBiomGrowthYield(JP,JY,JX));    HuskBiomGrowthYield=0._r8
  allocate(EarBiomGrowthYield(JP,JY,JX));    EarBiomGrowthYield=0._r8
  allocate(GrainBiomGrowthYield(JP,JY,JX));     GrainBiomGrowthYield=0._r8
  allocate(DMND(JP,JY,JX));     DMND=0._r8
  allocate(LeafBiomGrowthYield(JP,JY,JX));     LeafBiomGrowthYield=0._r8
  allocate(Hours4LenthenPhotoPeriod(MaxNumBranches,JP,JY,JX));  Hours4LenthenPhotoPeriod=0._r8
  allocate(Hours4ShortenPhotoPeriod(MaxNumBranches,JP,JY,JX));  Hours4ShortenPhotoPeriod=0._r8
  allocate(NumOfLeaves_brch(MaxNumBranches,JP,JY,JX));  NumOfLeaves_brch=0._r8
  allocate(VSTGX(MaxNumBranches,JP,JY,JX)); VSTGX=0._r8
  allocate(Hours4Leafout(MaxNumBranches,JP,JY,JX));  Hours4Leafout=0._r8
  allocate(Hours4LeafOff(MaxNumBranches,JP,JY,JX));  Hours4LeafOff=0._r8
  allocate(KLEAF(MaxNumBranches,JP,JY,JX)); KLEAF=0
  allocate(KVSTGN(MaxNumBranches,JP,JY,JX));KVSTGN=0
  allocate(KLEAFX(MaxNumBranches,JP,JY,JX));KLEAFX=0
  allocate(KLeafNodeNumber(MaxNumBranches,JP,JY,JX)); KLeafNodeNumber=0
  allocate(RefLeafAppearRate(JP,JY,JX));     RefLeafAppearRate=0._r8
  allocate(WDLF(JP,JY,JX));     WDLF=0._r8
  allocate(SLA1(JP,JY,JX));     SLA1=0._r8
  allocate(TCelciusChill4Leaf(JP,JY,JX));      TCelciusChill4Leaf=0._r8
  allocate(SSL1(JP,JY,JX));     SSL1=0._r8
  allocate(HourThreshold4LeafOut(JC,JP,JY,JX));  HourThreshold4LeafOut=0._r8
  allocate(NumOfBranches_pft(JP,JY,JX));      NumOfBranches_pft=0
  allocate(BranchNumber_pft(JP,JY,JX));      BranchNumber_pft=0
  allocate(BranchNumber_brch(MaxNumBranches,JP,JY,JX));  BranchNumber_brch=0
  allocate(NumOfMainBranch_pft(JP,JY,JX));      NumOfMainBranch_pft=0
  allocate(IFLGR(MaxNumBranches,JP,JY,JX)); IFLGR=0
  allocate(IFLGQ(MaxNumBranches,JP,JY,JX)); IFLGQ=0
  allocate(doSenescence_brch(MaxNumBranches,JP,JY,JX)); doSenescence_brch=0
  allocate(doRemobilization_brch(MaxNumBranches,JP,JY,JX)); doRemobilization_brch=0
  allocate(doInitLeafOut(MaxNumBranches,JP,JY,JX)); doInitLeafOut=0
  allocate(doPlantLeafOut(MaxNumBranches,JP,JY,JX)); doPlantLeafOut=0
  allocate(doPlantLeaveOff(MaxNumBranches,JP,JY,JX)); doPlantLeaveOff=0
  allocate(iPlantBranchState(MaxNumBranches,JP,JY,JX)); iPlantBranchState=0
  allocate(MinNonstructalC4InitBranch(JP,JY,JX));       MinNonstructalC4InitBranch=0._r8
  allocate(GSTGI(MaxNumBranches,JP,JY,JX)); GSTGI=0._r8
  allocate(DGSTGI(MaxNumBranches,JP,JY,JX));DGSTGI=0._r8
  allocate(DGSTGF(MaxNumBranches,JP,JY,JX));DGSTGF=0._r8
  allocate(ShootNodeNumber(MaxNumBranches,JP,JY,JX));  ShootNodeNumber=0._r8
  allocate(NodeNumberToInitFloral(MaxNumBranches,JP,JY,JX)); NodeNumberToInitFloral=0._r8
  allocate(GSTGF(MaxNumBranches,JP,JY,JX)); GSTGF=0._r8
  allocate(NodeNumberAtAnthesis(MaxNumBranches,JP,JY,JX)); NodeNumberAtAnthesis=0._r8
  allocate(TGSTGI(MaxNumBranches,JP,JY,JX));TGSTGI=0._r8
  allocate(TGSTGF(MaxNumBranches,JP,JY,JX));TGSTGF=0._r8
  allocate(RefNodeInitRate(JP,JY,JX));     RefNodeInitRate=0._r8
  allocate(SNL1(JP,JY,JX));     SNL1=0._r8
  allocate(FNOD(JP,JY,JX));     FNOD=0._r8
  allocate(NumConCurrentGrowinNode(JP,JY,JX));     NumConCurrentGrowinNode=0
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
  allocate(FLG4(MaxNumBranches,JP,JY,JX));  FLG4=0._r8
  allocate(HourCounter4LeafOut_brch(MaxNumBranches,JP,JY,JX));  HourCounter4LeafOut_brch=0._r8
  allocate(FLGZ(MaxNumBranches,JP,JY,JX));  FLGZ=0._r8
  allocate(iPlantCalendar(NumGrowthStages,MaxNumBranches,JP,JY,JX));iPlantCalendar=0
  allocate(TCelciusChill4Seed(JP,JY,JX));      TCelciusChill4Seed=0._r8
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
  call destroy(BranchAngle_pft)
  call destroy(PetioleAngle_pft)
  call destroy(SineBranchAngle_pft)
  call destroy(SinePetioleAngle_pft)
  call destroy(RAZ)
  call destroy(CanPHeight4WatUptake)
  call destroy(LeafAreaNode_brch)
  call destroy(PetioleLengthNode_brch)
  call destroy(InternodeHeightLive_brch)
  call destroy(LeafAreaLive_brch)
  call destroy(LeafAreaDying_brch)
  call destroy(CanPBranchHeight)
  call destroy(GRNOB)
  call destroy(GRNXB)
  call destroy(GRNO)
  call destroy(pftPlantPopulation)
  call destroy(InternodeHeightDying_brch)
  call destroy(CNLF)
  call destroy(CPLF)
  call destroy(CNSHE)
  call destroy(rNCStalk_pft)
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
  call destroy(iPlantInitThermoAdaptZone)
  call destroy(iPlantThermoAdaptZone)
  call destroy(MatureGroup_brch)
  call destroy(MatureGroup_pft)
  call destroy(GROUPX)
  call destroy(PPI)
  call destroy(StandingDeadInitC_pft)
  call destroy(PPX)
  call destroy(NumActivePlants)
  call destroy(PPT)
  call destroy(PPZ)
  call destroy(HoursCanopyPSITooLow)
  call destroy(PlantO2Stress)
  call destroy(fTgrowCanP)
  call destroy(TCG)
  call destroy(TKG)
  call destroy(PetioleBiomGrowthYield)
  call destroy(StalkBiomGrowthYield)
  call destroy(ReserveBiomGrowthYield)
  call destroy(HuskBiomGrowthYield)
  call destroy(EarBiomGrowthYield)
  call destroy(GrainBiomGrowthYield)
  call destroy(DMND)
  call destroy(LeafBiomGrowthYield)
  call destroy(Hours4LenthenPhotoPeriod)
  call destroy(Hours4ShortenPhotoPeriod)
  call destroy(NumOfLeaves_brch)
  call destroy(VSTGX)
  call destroy(Hours4Leafout)
  call destroy(Hours4LeafOff)
  call destroy(KLEAF)
  call destroy(KVSTGN)
  call destroy(KLEAFX)
  call destroy(KLeafNodeNumber)
  call destroy(RefLeafAppearRate)
  call destroy(WDLF)
  call destroy(SLA1)
  call destroy(TCelciusChill4Leaf)
  call destroy(SSL1)
  call destroy(HourThreshold4LeafOut)
  call destroy(NumOfBranches_pft)
  call destroy(BranchNumber_pft)
  call destroy(BranchNumber_brch)
  call destroy(NumOfMainBranch_pft)
  call destroy(IFLGR)
  call destroy(IFLGQ)
  call destroy(doSenescence_brch)
  call destroy(doRemobilization_brch)
  call destroy(doInitLeafOut)
  call destroy(doPlantLeafOut)
  call destroy(doPlantLeaveOff)
  call destroy(iPlantBranchState)
  call destroy(MinNonstructalC4InitBranch)
  call destroy(GSTGI)
  call destroy(DGSTGI)
  call destroy(DGSTGF)
  call destroy(ShootNodeNumber)
  call destroy(NodeNumberToInitFloral)
  call destroy(GSTGF)
  call destroy(NodeNumberAtAnthesis)
  call destroy(TGSTGI)
  call destroy(TGSTGF)
  call destroy(RefNodeInitRate)
  call destroy(SNL1)
  call destroy(FNOD)
  call destroy(NumConCurrentGrowinNode)
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
  call destroy(TCelciusChill4Seed)
  call destroy(HTC)
  call destroy(SSTX)
  call destroy(XDL)
  call destroy(XPPD)
  call destroy(ClumpFactort0)
  call destroy(HourThreshold4LeafOff)
  call destroy(OFFST)
  end subroutine DestructPlantTraits

end module PlantTraitDataType
