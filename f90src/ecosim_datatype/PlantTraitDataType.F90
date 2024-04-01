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
  real(r8) :: FracHour4LeafoffRemob(0:5)              !allocation parameter
  REAL(R8),target,allocatable :: FWODBE(:,:)             !woody element allocation
  real(r8),target,allocatable :: FWODLE(:,:)             !element allocation for leaf
  real(r8),target,allocatable :: FWODRE(:,:)
  real(r8),target,allocatable :: FWOODE(:,:)             !woody C allocation
  real(r8),target,allocatable :: PARTS_brch(:,:,:,:,:)       !C partitioning coefficient
  real(r8),target,allocatable ::  CanopyStemArea_lbrch(:,:,:,:,:)              !stem layer area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLeafArea_pft(:,:,:)                       !plant leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyArea_pft(:,:,:)                       !plant canopy leaf+stem/stalk area, [m2 d-2]
  real(r8),target,allocatable ::  ARLFX(:,:)                         !total canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyStemArea_pft(:,:,:)                      !plant stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyHeight_pft(:,:,:)                          !canopy height, [m]
  real(r8),target,allocatable ::  ARSTX(:,:)                         !total canopy stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLAgrid_lyr(:,:,:)                       !total leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyStemA_lyr(:,:,:)                       !total stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLeafArea_grd(:,:)                         !grid level plant canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  StemArea_grd(:,:)                     !total canopy stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyArea_grd(:,:)                         !canopy area of combined over the grid [m2 d-2]
  integer ,target,allocatable ::  NGTopRootLayer_pft(:,:,:)                           !soil layer at planting depth, [-]
  real(r8),target,allocatable ::  PlantinDepth(:,:,:)                      !planting depth, [m]
  real(r8),target,allocatable ::  SeedDepth_pft(:,:,:)                       !seeding depth, [m]
  real(r8),target,allocatable ::  SeedVolumeMean_pft(:,:,:)                        !seed volume, [m3 ]
  real(r8),target,allocatable ::  SeedMeanLen_pft(:,:,:)                        !seed length, [m]
  real(r8),target,allocatable ::  SeedAreaMean_pft(:,:,:)                        !seed surface area, [m2]
  real(r8),target,allocatable ::  HypoctoHeight_pft(:,:,:)                       !cotyledon height, [m]
  real(r8),target,allocatable ::  MaxCanopyHeight_grd(:,:)                            !canopy height , [m]
  real(r8),target,allocatable ::  CanopyHeightz_col(:,:,:)                          !canopy layer height , [m]
  real(r8),target,allocatable ::  BranchAngle_pft(:,:,:)                       !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  PetioleAngle_pft(:,:,:)                       !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  SineBranchAngle_pft(:,:,:)                       !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  SinePetioleAngle_pft(:,:,:)                       !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  RAZ(:,:,:)                         !canopy roughness height, [m]
  real(r8),target,allocatable ::  CanPHeight4WatUptake(:,:,:)                       !canopy height, [m]
  real(r8),target,allocatable ::  LeafAreaNode_brch(:,:,:,:,:)                    !leaf area, [m2 d-2]
  real(r8),target,allocatable ::  PetioleLengthNode_brch(:,:,:,:,:)                   !sheath height, [m]
  real(r8),target,allocatable ::  LiveInterNodeHight_brch(:,:,:,:,:)                  !internode height, [m]
  real(r8),target,allocatable ::  LeafAreaLive_brch(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  LeafAreaDying_brch(:,:,:,:)                     !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanPBranchHeight(:,:,:,:)                    !branch height, [m]
  real(r8),target,allocatable ::  SeedNumSet_brch(:,:,:,:)                     !branch grain number, [d-2]
  real(r8),target,allocatable ::  PotentialSeedSites_brch(:,:,:,:)                     !branch potential grain number, [d-2]
  real(r8),target,allocatable ::  CanopySeedNum_pft(:,:,:)                        !canopy grain number, [d-2]
  real(r8),target,allocatable ::  PlantPopulation_pft(:,:,:)             !plant population, [d-2]
  real(r8),target,allocatable ::  InternodeHeightDying_brch(:,:,:,:,:)                  !internode height, [m]
  real(r8),target,allocatable ::  CNLF(:,:,:)                        !maximum leaf N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPLF(:,:,:)                        !maximum leaf P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNSHE(:,:,:)                       !sheath N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCStalk_pft(:,:,:)                       !stalk N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCReserve_pft(:,:,:)                       !reserve N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCHusk_pft(:,:,:)                       !husk N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCEar_pft(:,:,:)                       !ear N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CNGR(:,:,:)                        !grain N:C ratio, [g g-1]
  real(r8),target,allocatable ::  NodulerNC_pft(:,:,:)                        !nodule N:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPSHE(:,:,:)                       !sheath P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCStalk_pft(:,:,:)                       !stalk P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCReserve_pft(:,:,:)                       !reserve P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCHusk_pft(:,:,:)                       !husk P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCEar_pft(:,:,:)                       !ear P:C ratio, [g g-1]
  real(r8),target,allocatable ::  CPGR(:,:,:)                        !grain P:C ratio, [g g-1]
  real(r8),target,allocatable ::  NodulerPC_pft(:,:,:)                        !nodule P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rCNNonstructRemob_pft(:,:,:)                        !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  rCPNonstructRemob_pft(:,:,:)                        !C:P ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  OSMO(:,:,:)                        !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8),target,allocatable ::  TC4LeafOff_pft(:,:,:)                         !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8),target,allocatable ::  PlantInitThermoAdaptZone(:,:,:)                       !initial plant thermal adaptation zone, [-]
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
  real(r8),target,allocatable ::  HoursCanopyPSITooLow_pft(:,:,:)                        !canopy plant water stress indicator, number of hours PSILT < PSILY, []
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
  real(r8),target,allocatable ::  NoduGrowthYield_pft(:,:,:)                        !nodule growth yield, [g g-1]
  real(r8),target,allocatable ::  LeafBiomGrowthYield(:,:,:)                        !leaf growth yield, [g g-1]
  real(r8),target,allocatable ::  Hours4LenthenPhotoPeriod_brch(:,:,:,:)                      !initial heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  Hours4ShortenPhotoPeriod_brch(:,:,:,:)                      !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8),target,allocatable ::  NumOfLeaves_brch(:,:,:,:)                      !leaf number, [-]
  real(r8),target,allocatable ::  LeafNumberAtFloralInit_brch(:,:,:,:)                     !leaf number at floral initiation, [-]
  real(r8),target,allocatable ::  Hours4Leafout_brch(:,:,:,:)                      !heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  Hours4LeafOff_brch(:,:,:,:)                      !cold requirement for autumn leafoff/hardening, [h]
  integer,target,allocatable ::  KLeafNumber_brch(:,:,:,:)                      !leaf number, [-]
  integer,target,allocatable ::  KLowestGroLeafNode_brch(:,:,:,:)                     !leaf growth stage counter, [-]
  integer,target,allocatable ::  KLEAFX(:,:,:,:)                     !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,target,allocatable ::  KHiestGroLeafNode_brch(:,:,:,:)                      !leaf growth stage counter, [-]
  real(r8),target,allocatable ::  RefLeafAppearRate_pft(:,:,:)                        !rate of leaf initiation, [h-1 at 25 oC]
  real(r8),target,allocatable ::  WDLF(:,:,:)                        !leaf length:width ratio, [-]
  real(r8),target,allocatable ::  SLA1(:,:,:)                        !leaf area:mass during growth, [m2 g-1]
  real(r8),target,allocatable ::  TC4LeafOut_pft(:,:,:)                         !threshold temperature for spring leafout/dehardening, [oC]
  real(r8),target,allocatable ::  PetoLen2Mass_pft(:,:,:)                        !petiole length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  HourReq4LeafOut_brch(:,:,:,:)                      !hours above threshold temperature required for spring leafout/dehardening, [-]
  integer,target,allocatable ::  NumOfBranches_pft(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  BranchNumber_pft(:,:,:)                          !branch number, [-]
  integer,target,allocatable ::  BranchNumber_brch(:,:,:,:)                       !branch number, [-]
  integer,target,allocatable ::  MainBranchNum_pft(:,:,:)                          !number of main branch, [-]
  integer,target,allocatable ::  Prep4Literfall_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  Hours4LiterfalAftMature_brch(:,:,:,:)                      !branch phenology flag, [h]
  integer,target,allocatable ::  doSenescence_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doRemobilization_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doInitLeafOut_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doPlantLeafOut_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  doPlantLeaveOff_brch(:,:,:,:)                      !branch phenology flag, [-]
  integer,target,allocatable ::  iPlantBranchState_brch(:,:,:,:)                      !flag to detect branch death , [-]
  real(r8),target,allocatable ::  MinNonstC2InitBranch_pft(:,:,:)                          !branch nonstructural C content required for new branch, [g g-1]
  real(r8),target,allocatable ::  NodeNumNormByMatgrp_brch(:,:,:,:)                     !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  HourlyNodeNumNormByMatgrp_brch(:,:,:,:)                    !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8),target,allocatable ::  dReproNodeNumNormByMatG_brch(:,:,:,:)                    !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8),target,allocatable ::  ShootNodeNum_brch(:,:,:,:)                      !node number, [-]
  real(r8),target,allocatable ::  NodeNum2InitFloral_brch(:,:,:,:)                     !node number at floral initiation, [-]
  real(r8),target,allocatable ::  ReprodNodeNumNormByMatrgrp_brch(:,:,:,:)                     !normalized node number during reproductive growth stages, [-]
  real(r8),target,allocatable ::  NodeNumberAtAnthesis_brch(:,:,:,:)                     !node number at anthesis, [-]
  real(r8),target,allocatable ::  TotalNodeNumNormByMatgrp_brch(:,:,:,:)                    !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  TotReproNodeNumNormByMatrgrp_brch(:,:,:,:)                    !normalized node number during reproductive growth stages , [-]
  real(r8),target,allocatable ::  RefNodeInitRate_pft(:,:,:)                        !rate of node initiation, [h-1 at 25 oC]
  real(r8),target,allocatable ::  NodeLenPergC(:,:,:)                        !internode length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  FNOD(:,:,:)                        !parameter for allocation of growth to nodes, [-]
  integer,target,allocatable ::  NumCogrowNode(:,:,:)                         !number of concurrently growing nodes
  real(r8),target,allocatable ::  PSICanPDailyMin(:,:,:)                       !minimum daily canopy water potential, [MPa]
  real(r8),target,allocatable ::  ClumpFactorCurrent_pft(:,:,:)                         !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8),target,allocatable ::  ClumpFactor(:,:,:)                          !clumping factor for self-shading in canopy layer, [-]
  integer,target,allocatable ::  iPlantShootState_pft(:,:,:)                        !flag to detect canopy death , [-]
  real(r8),target,allocatable ::  MaxPotentSeedNumber_pft(:,:,:)                        !maximum grain node number per branch, [-]
  real(r8),target,allocatable ::  MaxSeedNumPerSite_pft(:,:,:)                        !maximum grain number per node , [-]
  real(r8),target,allocatable ::  MaxSeedCMass(:,:,:)                        !maximum grain size   , [g]
  real(r8),target,allocatable ::  ShootNodeNumAtPlanting_pft(:,:,:)                        !number of nodes in seed, [-]
  real(r8),target,allocatable ::  SeedCMass(:,:,:)                        !grain size at seeding, [g]
  real(r8),target,allocatable ::  GrainFillRate25C_pft(:,:,:)                       !maximum rate of fill per grain, [g h-1]
  real(r8),target,allocatable ::  HourFailGrainFill_brch(:,:,:,:)                      !flag to detect physiological maturity from  grain fill , [-]
  real(r8),target,allocatable ::  Hours2LeafOut_brch(:,:,:,:)                      !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  HoursDoingRemob_brch(:,:,:,:)                      !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  integer,target,allocatable ::  iPlantCalendar_brch(:,:,:,:,:)                     !plant growth stage, [-]
  real(r8),target,allocatable ::  TCChill4Seed_pft(:,:,:)                         !temperature below which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  HighTCLimtSeed_pft(:,:,:)                         !temperature above which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  SSTX(:,:,:)                        !sensitivity to HTC (seeds oC-1 above HTC)
  real(r8),target,allocatable ::  CriticPhotoPeriod_pft(:,:,:)                         !critical daylength for phenological progress, [h]
  real(r8),target,allocatable ::  PhotoPeriodSens_pft(:,:,:)                        !difference between current and critical daylengths used to calculate  phenological progress, [h]
  real(r8),target,allocatable ::  ClumpFactorInit_pft(:,:,:)                         !initial clumping factor for self-shading in canopy layer, [-]
  real(r8),target,allocatable ::  HourReq4LeafOff_brch(:,:,:,:)                      !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8),target,allocatable ::  OFFST(:,:,:)                       !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
!----------------------------------------------------------------------

contains
  subroutine InitPlantTraits(NumOfPlantLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs

  FracHour4LeafoffRemob =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)
  allocate(FWODLE(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FWODLE=0._r8
  allocate(FWODBE(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FWODBE=0._r8
  allocate(FWODRE(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FWODRE=0._r8         !
  allocate(FWOODE(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FWOODE=0._r8         !woody element allocation
  allocate(CanopyStemArea_lbrch(NumOfCanopyLayers,MaxNumBranches,JP,JY,JX));CanopyStemArea_lbrch=0._r8
  allocate(CanopyLeafArea_pft(JP,JY,JX));    CanopyLeafArea_pft=0._r8
  allocate(CanopyArea_pft(JP,JY,JX));    CanopyArea_pft=0._r8
  allocate(ARLFX(JY,JX));       ARLFX=0._r8
  allocate(CanopyStemArea_pft(JP,JY,JX));    CanopyStemArea_pft=0._r8
  allocate(CanopyHeight_pft(JP,JY,JX));       CanopyHeight_pft=0._r8
  allocate(ARSTX(JY,JX));       ARSTX=0._r8
  allocate(CanopyLAgrid_lyr(NumOfCanopyLayers,JY,JX));    CanopyLAgrid_lyr=0._r8
  allocate(CanopyStemA_lyr(NumOfCanopyLayers,JY,JX));    CanopyStemA_lyr=0._r8
  allocate(CanopyLeafArea_grd(JY,JX));       CanopyLeafArea_grd=0._r8
  allocate(StemArea_grd(JY,JX));       StemArea_grd=0._r8
  allocate(CanopyArea_grd(JY,JX));       CanopyArea_grd=0._r8
  allocate(NGTopRootLayer_pft(JP,JY,JX));       NGTopRootLayer_pft=0
  allocate(PlantinDepth(JP,JY,JX));   PlantinDepth=0._r8
  allocate(SeedDepth_pft(JP,JY,JX));    SeedDepth_pft=0._r8
  allocate(SeedVolumeMean_pft(JP,JY,JX));     SeedVolumeMean_pft=0._r8
  allocate(SeedMeanLen_pft(JP,JY,JX));     SeedMeanLen_pft=0._r8
  allocate(SeedAreaMean_pft(JP,JY,JX));     SeedAreaMean_pft=0._r8
  allocate(HypoctoHeight_pft(JP,JY,JX));    HypoctoHeight_pft=0._r8
  allocate(MaxCanopyHeight_grd(JY,JX));          MaxCanopyHeight_grd=0._r8
  allocate(CanopyHeightz_col(0:NumOfCanopyLayers,JY,JX));     CanopyHeightz_col=0._r8
  allocate(BranchAngle_pft(JP,JY,JX));    BranchAngle_pft=0._r8
  allocate(PetioleAngle_pft(JP,JY,JX));    PetioleAngle_pft=0._r8
  allocate(SineBranchAngle_pft(JP,JY,JX));    SineBranchAngle_pft=0._r8
  allocate(SinePetioleAngle_pft(JP,JY,JX));    SinePetioleAngle_pft=0._r8
  allocate(RAZ(JP,JY,JX));      RAZ=0._r8
  allocate(CanPHeight4WatUptake(JP,JY,JX));    CanPHeight4WatUptake=0._r8
  allocate(PARTS_brch(NumOfPlantMorphUnits,MaxNumBranches,JP,JY,JX));PARTS_brch=0._r8
  allocate(LeafAreaNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafAreaNode_brch=0._r8
  allocate(PetioleLengthNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleLengthNode_brch=0._r8
  allocate(LiveInterNodeHight_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LiveInterNodeHight_brch=0._r8
  allocate(LeafAreaLive_brch(MaxNumBranches,JP,JY,JX)); LeafAreaLive_brch=0._r8
  allocate(LeafAreaDying_brch(MaxNumBranches,JP,JY,JX)); LeafAreaDying_brch=0._r8
  allocate(CanPBranchHeight(MaxNumBranches,JP,JY,JX));CanPBranchHeight=0._r8
  allocate(SeedNumSet_brch(MaxNumBranches,JP,JY,JX)); SeedNumSet_brch=0._r8
  allocate(PotentialSeedSites_brch(MaxNumBranches,JP,JY,JX)); PotentialSeedSites_brch=0._r8
  allocate(CanopySeedNum_pft(JP,JY,JX));     CanopySeedNum_pft=0._r8
  allocate(PlantPopulation_pft(JP,JY,JX));       PlantPopulation_pft=0._r8
  allocate(InternodeHeightDying_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));InternodeHeightDying_brch=0._r8
  allocate(CNLF(JP,JY,JX));     CNLF=0._r8
  allocate(CPLF(JP,JY,JX));     CPLF=0._r8
  allocate(CNSHE(JP,JY,JX));    CNSHE=0._r8
  allocate(rNCStalk_pft(JP,JY,JX));    rNCStalk_pft=0._r8
  allocate(rNCReserve_pft(JP,JY,JX));    rNCReserve_pft=0._r8
  allocate(rNCHusk_pft(JP,JY,JX));    rNCHusk_pft=0._r8
  allocate(rNCEar_pft(JP,JY,JX));    rNCEar_pft=0._r8
  allocate(CNGR(JP,JY,JX));     CNGR=0._r8
  allocate(NodulerNC_pft(JP,JY,JX));     NodulerNC_pft=0._r8
  allocate(CPSHE(JP,JY,JX));    CPSHE=0._r8
  allocate(rPCStalk_pft(JP,JY,JX));    rPCStalk_pft=0._r8
  allocate(rPCReserve_pft(JP,JY,JX));    rPCReserve_pft=0._r8
  allocate(rPCHusk_pft(JP,JY,JX));    rPCHusk_pft=0._r8
  allocate(rPCEar_pft(JP,JY,JX));    rPCEar_pft=0._r8
  allocate(CPGR(JP,JY,JX));     CPGR=0._r8
  allocate(NodulerPC_pft(JP,JY,JX));     NodulerPC_pft=0._r8
  allocate(rCNNonstructRemob_pft(JP,JY,JX));     rCNNonstructRemob_pft=0._r8
  allocate(rCPNonstructRemob_pft(JP,JY,JX));     rCPNonstructRemob_pft=0._r8
  allocate(OSMO(JP,JY,JX));     OSMO=0._r8
  allocate(TC4LeafOff_pft(JP,JY,JX));      TC4LeafOff_pft=0._r8
  allocate(PlantInitThermoAdaptZone(JP,JY,JX));    PlantInitThermoAdaptZone=0._r8
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
  allocate(HoursCanopyPSITooLow_pft(JP,JY,JX));     HoursCanopyPSITooLow_pft=0._r8
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
  allocate(NoduGrowthYield_pft(JP,JY,JX));     NoduGrowthYield_pft=0._r8
  allocate(LeafBiomGrowthYield(JP,JY,JX));     LeafBiomGrowthYield=0._r8
  allocate(Hours4LenthenPhotoPeriod_brch(MaxNumBranches,JP,JY,JX));  Hours4LenthenPhotoPeriod_brch=0._r8
  allocate(Hours4ShortenPhotoPeriod_brch(MaxNumBranches,JP,JY,JX));  Hours4ShortenPhotoPeriod_brch=0._r8
  allocate(NumOfLeaves_brch(MaxNumBranches,JP,JY,JX));  NumOfLeaves_brch=0._r8
  allocate(LeafNumberAtFloralInit_brch(MaxNumBranches,JP,JY,JX)); LeafNumberAtFloralInit_brch=0._r8
  allocate(Hours4Leafout_brch(MaxNumBranches,JP,JY,JX));  Hours4Leafout_brch=0._r8
  allocate(Hours4LeafOff_brch(MaxNumBranches,JP,JY,JX));  Hours4LeafOff_brch=0._r8
  allocate(KLeafNumber_brch(MaxNumBranches,JP,JY,JX)); KLeafNumber_brch=0
  allocate(KLowestGroLeafNode_brch(MaxNumBranches,JP,JY,JX));KLowestGroLeafNode_brch=0
  allocate(KLEAFX(MaxNumBranches,JP,JY,JX));KLEAFX=0
  allocate(KHiestGroLeafNode_brch(MaxNumBranches,JP,JY,JX)); KHiestGroLeafNode_brch=0
  allocate(RefLeafAppearRate_pft(JP,JY,JX));     RefLeafAppearRate_pft=0._r8
  allocate(WDLF(JP,JY,JX));     WDLF=0._r8
  allocate(SLA1(JP,JY,JX));     SLA1=0._r8
  allocate(TC4LeafOut_pft(JP,JY,JX));      TC4LeafOut_pft=0._r8
  allocate(PetoLen2Mass_pft(JP,JY,JX));     PetoLen2Mass_pft=0._r8
  allocate(HourReq4LeafOut_brch(NumOfCanopyLayers,JP,JY,JX));  HourReq4LeafOut_brch=0._r8
  allocate(NumOfBranches_pft(JP,JY,JX));      NumOfBranches_pft=0
  allocate(BranchNumber_pft(JP,JY,JX));      BranchNumber_pft=0
  allocate(BranchNumber_brch(MaxNumBranches,JP,JY,JX));  BranchNumber_brch=0
  allocate(MainBranchNum_pft(JP,JY,JX));      MainBranchNum_pft=0
  allocate(Prep4Literfall_brch(MaxNumBranches,JP,JY,JX)); Prep4Literfall_brch=ifalse
  allocate(Hours4LiterfalAftMature_brch(MaxNumBranches,JP,JY,JX)); Hours4LiterfalAftMature_brch=0
  allocate(doSenescence_brch(MaxNumBranches,JP,JY,JX)); doSenescence_brch=0
  allocate(doRemobilization_brch(MaxNumBranches,JP,JY,JX)); doRemobilization_brch=0
  allocate(doInitLeafOut_brch(MaxNumBranches,JP,JY,JX)); doInitLeafOut_brch=0
  allocate(doPlantLeafOut_brch(MaxNumBranches,JP,JY,JX)); doPlantLeafOut_brch=0
  allocate(doPlantLeaveOff_brch(MaxNumBranches,JP,JY,JX)); doPlantLeaveOff_brch=0
  allocate(iPlantBranchState_brch(MaxNumBranches,JP,JY,JX)); iPlantBranchState_brch=0
  allocate(MinNonstC2InitBranch_pft(JP,JY,JX));       MinNonstC2InitBranch_pft=0._r8
  allocate(NodeNumNormByMatgrp_brch(MaxNumBranches,JP,JY,JX)); NodeNumNormByMatgrp_brch=0._r8
  allocate(HourlyNodeNumNormByMatgrp_brch(MaxNumBranches,JP,JY,JX));HourlyNodeNumNormByMatgrp_brch=0._r8
  allocate(dReproNodeNumNormByMatG_brch(MaxNumBranches,JP,JY,JX));dReproNodeNumNormByMatG_brch=0._r8
  allocate(ShootNodeNum_brch(MaxNumBranches,JP,JY,JX));  ShootNodeNum_brch=0._r8
  allocate(NodeNum2InitFloral_brch(MaxNumBranches,JP,JY,JX)); NodeNum2InitFloral_brch=0._r8
  allocate(ReprodNodeNumNormByMatrgrp_brch(MaxNumBranches,JP,JY,JX)); ReprodNodeNumNormByMatrgrp_brch=0._r8
  allocate(NodeNumberAtAnthesis_brch(MaxNumBranches,JP,JY,JX)); NodeNumberAtAnthesis_brch=0._r8
  allocate(TotalNodeNumNormByMatgrp_brch(MaxNumBranches,JP,JY,JX));TotalNodeNumNormByMatgrp_brch=0._r8
  allocate(TotReproNodeNumNormByMatrgrp_brch(MaxNumBranches,JP,JY,JX));TotReproNodeNumNormByMatrgrp_brch=0._r8
  allocate(RefNodeInitRate_pft(JP,JY,JX));     RefNodeInitRate_pft=0._r8
  allocate(NodeLenPergC(JP,JY,JX));     NodeLenPergC=0._r8
  allocate(FNOD(JP,JY,JX));     FNOD=0._r8
  allocate(NumCogrowNode(JP,JY,JX));     NumCogrowNode=0
  allocate(PSICanPDailyMin(JP,JY,JX));    PSICanPDailyMin=0._r8
  allocate(ClumpFactorCurrent_pft(JP,JY,JX));      ClumpFactorCurrent_pft=0._r8
  allocate(ClumpFactor(JP,JY,JX));       ClumpFactor=0._r8
  allocate(iPlantShootState_pft(JP,JY,JX));    iPlantShootState_pft=0
  allocate(MaxPotentSeedNumber_pft(JP,JY,JX));     MaxPotentSeedNumber_pft=0._r8
  allocate(MaxSeedNumPerSite_pft(JP,JY,JX));     MaxSeedNumPerSite_pft=0._r8
  allocate(MaxSeedCMass(JP,JY,JX));     MaxSeedCMass=0._r8
  allocate(ShootNodeNumAtPlanting_pft(JP,JY,JX));     ShootNodeNumAtPlanting_pft=0._r8
  allocate(SeedCMass(JP,JY,JX));     SeedCMass=0._r8
  allocate(GrainFillRate25C_pft(JP,JY,JX));    GrainFillRate25C_pft=0._r8
  allocate(HourFailGrainFill_brch(MaxNumBranches,JP,JY,JX));  HourFailGrainFill_brch=0._r8
  allocate(Hours2LeafOut_brch(MaxNumBranches,JP,JY,JX));  Hours2LeafOut_brch=0._r8
  allocate(HoursDoingRemob_brch(MaxNumBranches,JP,JY,JX));  HoursDoingRemob_brch=0._r8
  allocate(iPlantCalendar_brch(NumGrowthStages,MaxNumBranches,JP,JY,JX));iPlantCalendar_brch=0
  allocate(TCChill4Seed_pft(JP,JY,JX));      TCChill4Seed_pft=0._r8
  allocate(HighTCLimtSeed_pft(JP,JY,JX));      HighTCLimtSeed_pft=0._r8
  allocate(SSTX(JP,JY,JX));     SSTX=0._r8
  allocate(CriticPhotoPeriod_pft(JP,JY,JX));      CriticPhotoPeriod_pft=0._r8
  allocate(PhotoPeriodSens_pft(JP,JY,JX));     PhotoPeriodSens_pft=0._r8
  allocate(ClumpFactorInit_pft(JP,JY,JX));      ClumpFactorInit_pft=0._r8
  allocate(HourReq4LeafOff_brch(NumOfCanopyLayers,JP,JY,JX));  HourReq4LeafOff_brch=0._r8
  allocate(OFFST(JP,JY,JX));    OFFST=0._r8
  end subroutine InitPlantTraits

!----------------------------------------------------------------------
  subroutine DestructPlantTraits
  use abortutils, only : destroy
  implicit none

  call destroy(FWODBE)
  call destroy(FWODRE)
  call destroy(FWOODE)
  call destroy(CanopyStemArea_lbrch)
  call destroy(CanopyLeafArea_pft)
  call destroy(CanopyArea_pft)
  call destroy(ARLFX)
  call destroy(CanopyStemArea_pft)
  call destroy(CanopyHeight_pft)
  call destroy(ARSTX)
  call destroy(CanopyLAgrid_lyr)
  call destroy(CanopyStemA_lyr)
  call destroy(CanopyLeafArea_grd)
  call destroy(StemArea_grd)
  call destroy(CanopyArea_grd)
  call destroy(NGTopRootLayer_pft)
  call destroy(PlantinDepth)
  call destroy(SeedDepth_pft)
  call destroy(SeedVolumeMean_pft)
  call destroy(SeedMeanLen_pft)
  call destroy(SeedAreaMean_pft)
  call destroy(HypoctoHeight_pft)
  call destroy(MaxCanopyHeight_grd)
  call destroy(CanopyHeightz_col)
  call destroy(BranchAngle_pft)
  call destroy(PetioleAngle_pft)
  call destroy(SineBranchAngle_pft)
  call destroy(SinePetioleAngle_pft)
  call destroy(RAZ)
  call destroy(CanPHeight4WatUptake)
  call destroy(LeafAreaNode_brch)
  call destroy(PetioleLengthNode_brch)
  call destroy(LiveInterNodeHight_brch)
  call destroy(LeafAreaLive_brch)
  call destroy(LeafAreaDying_brch)
  call destroy(CanPBranchHeight)
  call destroy(SeedNumSet_brch)
  call destroy(PotentialSeedSites_brch)
  call destroy(CanopySeedNum_pft)
  call destroy(PlantPopulation_pft)
  call destroy(InternodeHeightDying_brch)
  call destroy(PARTS_brch)
  call destroy(CNLF)
  call destroy(CPLF)
  call destroy(CNSHE)
  call destroy(rNCStalk_pft)
  call destroy(rNCReserve_pft)
  call destroy(rNCHusk_pft)
  call destroy(rNCEar_pft)
  call destroy(CNGR)
  call destroy(NodulerNC_pft)
  call destroy(CPSHE)
  call destroy(rPCStalk_pft)
  call destroy(rPCReserve_pft)
  call destroy(rPCHusk_pft)
  call destroy(rPCEar_pft)
  call destroy(CPGR)
  call destroy(NodulerPC_pft)
  call destroy(rCNNonstructRemob_pft)
  call destroy(rCPNonstructRemob_pft)
  call destroy(OSMO)
  call destroy(TC4LeafOff_pft)
  call destroy(PlantInitThermoAdaptZone)
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
  call destroy(HoursCanopyPSITooLow_pft)
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
  call destroy(NoduGrowthYield_pft)
  call destroy(LeafBiomGrowthYield)
  call destroy(Hours4LenthenPhotoPeriod_brch)
  call destroy(Hours4ShortenPhotoPeriod_brch)
  call destroy(NumOfLeaves_brch)
  call destroy(LeafNumberAtFloralInit_brch)
  call destroy(Hours4Leafout_brch)
  call destroy(Hours4LeafOff_brch)
  call destroy(KLeafNumber_brch)
  call destroy(KLowestGroLeafNode_brch)
  call destroy(KLEAFX)
  call destroy(KHiestGroLeafNode_brch)
  call destroy(RefLeafAppearRate_pft)
  call destroy(WDLF)
  call destroy(SLA1)
  call destroy(TC4LeafOut_pft)
  call destroy(PetoLen2Mass_pft)
  call destroy(HourReq4LeafOut_brch)
  call destroy(NumOfBranches_pft)
  call destroy(BranchNumber_pft)
  call destroy(BranchNumber_brch)
  call destroy(MainBranchNum_pft)
  call destroy(Prep4Literfall_brch)
  call destroy(Hours4LiterfalAftMature_brch)
  call destroy(doSenescence_brch)
  call destroy(doRemobilization_brch)
  call destroy(doInitLeafOut_brch)
  call destroy(doPlantLeafOut_brch)
  call destroy(doPlantLeaveOff_brch)
  call destroy(iPlantBranchState_brch)
  call destroy(MinNonstC2InitBranch_pft)
  call destroy(NodeNumNormByMatgrp_brch)
  call destroy(HourlyNodeNumNormByMatgrp_brch)
  call destroy(dReproNodeNumNormByMatG_brch)
  call destroy(ShootNodeNum_brch)
  call destroy(NodeNum2InitFloral_brch)
  call destroy(ReprodNodeNumNormByMatrgrp_brch)
  call destroy(NodeNumberAtAnthesis_brch)
  call destroy(TotalNodeNumNormByMatgrp_brch)
  call destroy(TotReproNodeNumNormByMatrgrp_brch)
  call destroy(RefNodeInitRate_pft)
  call destroy(NodeLenPergC)
  call destroy(FNOD)
  call destroy(NumCogrowNode)
  call destroy(PSICanPDailyMin)
  call destroy(ClumpFactorCurrent_pft)
  call destroy(ClumpFactor)
  call destroy(iPlantShootState_pft)
  call destroy(MaxPotentSeedNumber_pft)
  call destroy(MaxSeedNumPerSite_pft)
  call destroy(MaxSeedCMass)
  call destroy(ShootNodeNumAtPlanting_pft)
  call destroy(SeedCMass)
  call destroy(GrainFillRate25C_pft)
  call destroy(HourFailGrainFill_brch)
  call destroy(Hours2LeafOut_brch)
  call destroy(HoursDoingRemob_brch)
  call destroy(iPlantCalendar_brch)
  call destroy(TCChill4Seed_pft)
  call destroy(HighTCLimtSeed_pft)
  call destroy(SSTX)
  call destroy(CriticPhotoPeriod_pft)
  call destroy(PhotoPeriodSens_pft)
  call destroy(ClumpFactorInit_pft)
  call destroy(HourReq4LeafOff_brch)
  call destroy(OFFST)
  end subroutine DestructPlantTraits

end module PlantTraitDataType
