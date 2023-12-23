module CanopyDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use EcoSIMConfig, only : jsken => jskenc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  CanopyPARalbedo_pft(:,:,:)                        !canopy PAR albedo , [-]
  real(r8),target,allocatable ::  TAUP(:,:,:)                        !canopy PAR transmissivity , [-]
  real(r8),target,allocatable ::  CanopySWabsorpty_pft(:,:,:)                        !canopy shortwave absorptivity , [-]
  real(r8),target,allocatable ::  CanopyPARabsorpty_pft(:,:,:)                        !canopy PAR absorptivity , [-]
  real(r8),target,allocatable ::  RSMX(:,:,:)                        !maximum stomatal resistance to vapor, [s m-1]
  real(r8),target,allocatable ::  CO2CuticleResist_pft(:,:,:)                        !maximum stomatal resistance to CO2, [s h-1]
  real(r8),target,allocatable ::  MaxCanPStomaResistH2O_pft(:,:,:)      !maximum stomatal resistance to vapor, [s h-1]
  real(r8),target,allocatable ::  RCS(:,:,:)                         !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8),target,allocatable ::  CanPStomaResistH2O_pft(:,:,:)         !canopy stomatal resistance, [h m-1]
  real(r8),target,allocatable ::  MinCanPStomaResistH2O_pft(:,:,:)      !canopy minimum stomatal resistance, [s m-1]
  real(r8),target,allocatable ::  BndlResistCanG(:,:)                           !canopy boundary layer resistance, [m h-1]
  real(r8),target,allocatable ::  O2I(:,:,:)                         !leaf gaseous O2 concentration, [umol m-3]
  real(r8),target,allocatable ::  LeafIntracellularCO2_pft(:,:,:)                        !leaf gaseous CO2 concentration, [umol m-3]
  real(r8),target,allocatable ::  AirConc_pft(:,:,:)                        !total gas concentration, [mol m-3]
  real(r8),target,allocatable ::  DiffCO2Atmos2Intracel_pft(:,:,:)                        !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8),target,allocatable ::  CanopyGasCO2_pft(:,:,:)                        !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8),target,allocatable ::  aquCO2Intraleaf_pft(:,:,:)                        !leaf aqueous CO2 concentration, [uM]
  real(r8),target,allocatable ::  O2L(:,:,:)                         !leaf aqueous O2 concentration, [uM]
  real(r8),target,allocatable ::  CO2Solubility_pft(:,:,:)                        !leaf CO2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  LeafO2Solubility_pft(:,:,:)                         !leaf O2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  Km4LeafaqCO2_pft(:,:,:)                      !leaf aqueous CO2 Km no O2, [uM]
  real(r8),target,allocatable ::  Km4RubiscoCarboxy_pft(:,:,:)                      !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8),target,allocatable ::  CHILL(:,:,:)                       !chilling effect on CO2 fixation, [-]
  real(r8),target,allocatable ::  Vmax4RubiscoCarboxy_pft(:,:,:,:,:)                   !maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  CO2lmtRubiscoCarboxyRate_node(:,:,:,:,:)                    !carboxylation rate, [umol m-2 s-1]
  real(r8),target,allocatable ::  CO2CompenPoint_node(:,:,:,:,:)                   !CO2 compensation point, [uM]
  real(r8),target,allocatable ::  LigthSatCarboxyRate_node(:,:,:,:,:)                   !maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  RubiscoCarboxyEff_node(:,:,:,:,:)                    !carboxylation efficiency, [umol umol-1]
  real(r8),target,allocatable ::  CMassCO2BundleSheath_node(:,:,:,:,:)                    !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  Vmax4PEPCarboxy_pft(:,:,:,:,:)                   !maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  CO2lmtPEPCarboxyRate_node(:,:,:,:,:)                   !C4 carboxylation rate, [umol m-2 s-1]
  real(r8),target,allocatable ::  LigthSatC4CarboxyRate_node(:,:,:,:,:)                   !maximum  light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  C4CarboxyEff_node(:,:,:,:,:)                   !C4 carboxylation efficiency, [umol umol-1]
  real(r8),target,allocatable ::  CPOOL4(:,:,:,:,:)                  !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  CMassHCO3BundleSheath_node(:,:,:,:,:)                    !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  RubiscoActivity_brch(:,:,:,:)                      !branch down-regulation of CO2 fixation, [-]
  real(r8),target,allocatable ::  NutrientCtrlonC4Carboxy_node(:,:,:,:,:)                   !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  C4PhotosynDowreg_brch(:,:,:,:)                     !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  CNETX(:,:)                         !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  VmaxRubCarboxyRef_pft(:,:,:)                        !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VmaxRubOxyRef_pft(:,:,:)                        !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VmaxPEPCarboxyRef_pft(:,:,:)                       !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  XKCO2(:,:,:)                       !Km for rubisco carboxylase activity, [uM]
  real(r8),target,allocatable ::  XKO2(:,:,:)                        !Km for rubisco oxygenase activity, [uM]
  real(r8),target,allocatable ::  Km4PEPCarboxy_pft(:,:,:)                      !Km for PEP carboxylase activity, [uM]
  real(r8),target,allocatable ::  LeafRuBPConc_pft(:,:,:)                        !leaf rubisco content, [g g-1]
  real(r8),target,allocatable ::  FracLeafProtinAsPEPCarboxyl_pft(:,:,:)                        !leaf PEP carboxylase content, [g g-1]
  real(r8),target,allocatable ::  SpecChloryfilAct_pft(:,:,:)                        !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  LeafC3ChlorofilConc_pft(:,:,:)                         !leaf C3 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  LeafC4ChlorofilConc_pft(:,:,:)                        !leaf C4 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  CanPCi2CaRatio(:,:,:)                        !Ci:Ca ratio, [-]
  real(r8),target,allocatable ::  RadNet2CanP(:,:,:)                 !canopy net radiation , [MJ d-2 h-1] >0
  real(r8),target,allocatable ::  LWRadCanP(:,:,:)                   !canopy longwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadSWbyCanopy_pft(:,:,:)                 !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadPARbyCanopy_pft(:,:,:)                   !canopy absorbed PAR , [umol m-2 s-1]
  real(r8),target,allocatable ::  FracRadPARbyCanopy_pft(:,:,:)                       !fraction of incoming PAR absorbed by canopy, [-]
  real(r8),target,allocatable ::  TAU_RadThru(:,:,:)                        !fraction of radiation transmitted by canopy layer, [-]
  real(r8),target,allocatable ::  TAU_RadCapt(:,:,:)                        !fraction of radiation intercepted by canopy layer, [-]
  real(r8),target,allocatable ::  FracSWRad2Grnd(:,:)                         !fraction of radiation intercepted by ground surface, [-]
  real(r8),target,allocatable ::  SWRadOnGrnd(:,:)                          !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  LWRadCanGPrev(:,:)                        !longwave radiation emitted by canopy, [MJ m-2 h-1]
  real(r8),target,allocatable ::  LWRadGrnd(:,:)                        !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  CanH2OHeldVg(:,:)                  !grid canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  TFLWCI(:,:)                        !net ice transfer to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  PrecIntcptByCanG(:,:)              !grid net precipitation water interception to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  EvapTransHeatP(:,:,:)                       !canopy latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatXAir2PCan(:,:,:)               !air to canopy sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatStorCanP(:,:,:)                       !canopy storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  ENGYX(:,:,:)                       !canopy heat storage from previous time step, [MJ d-2]
  real(r8),target,allocatable ::  VHeatCapCanP(:,:,:)                       !canopy heat capacity, [MJ d-2 K-1]
  real(r8),target,allocatable ::  PSICanopy_pft(:,:,:)                     !plant canopy total water potential , [Mpa]
  real(r8),target,allocatable ::  PSICanopyTurg_pft(:,:,:)                       !plant canopy turgor water potential, [Mpa]
  real(r8),target,allocatable ::  PSICanopyOsmo_pft(:,:,:)                 !platn canopy osmotic water potential, [Mpa]
  real(r8),target,allocatable ::  CanopyBndlResist_pft(:,:,:)                          !canopy boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  Transpiration_pft(:,:,:)                          !canopy transpiration, [m2 d-2 h-1]
  real(r8),target,allocatable ::  VapXAir2Canopy_pft(:,:,:)                !negative of canopy evaporation, [m2 d-2 h-1]
  real(r8),target,allocatable ::  CanopyWater_pft(:,:,:)                       !canopy water content associated with dry matter, [m3 d-2]
  real(r8),target,allocatable ::  TEVAPP(:,:)                        !total canopy evaporation + transpiration, [m3 d-2]
  real(r8),target,allocatable ::  VapXAir2CanG(:,:)                        !total canopy evaporation, [m3 d-2]
  real(r8),target,allocatable ::  TENGYC(:,:)                        !total canopy heat content, [MJ  d-2]
  real(r8),target,allocatable ::  THFLXC(:,:)                        !total canopy heat flux, [MJ  d-2]
  real(r8),target,allocatable ::  CanWatg(:,:)                       !total canopy water content stored in dry matter, [m3 d-2]
  real(r8),target,allocatable ::  LWRadCanG(:,:)                         !total canopy LW emission, [MJ d-2 h-1]
  real(r8),target,allocatable ::  CanopySWAlbedo_pft(:,:,:)                        !canopy shortwave albedo , [-]
  real(r8),target,allocatable ::  TAUR(:,:,:)                        !canopy shortwave transmissivity , [-]
  real(r8),target,allocatable ::  PrecIntcptByCanopy_pft(:,:,:)                        !water flux into plant canopy, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WatByPCanopy(:,:,:)                   !canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  TKC(:,:,:)                         !canopy temperature, [K]
  real(r8),target,allocatable ::  TCelciusCanopy_pft(:,:,:)                         !canopy temperature, [oC]
  real(r8),target,allocatable ::  DTKC(:,:,:)                        !change in canopy temperature, [K]
  real(r8),target,allocatable ::  TKCanopy_pft(:,:,:)                        !canopy temperature, [K]
  real(r8),target,allocatable ::  CPOOL3(:,:,:,:,:)                  !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8),target,allocatable ::  NetCumElmntFlx2Plant_pft(:,:,:,:)                     !effect of canopy element status on seed set , []
  real(r8),target,allocatable ::  WGLFT(:,:,:)                       !total leaf mass, [g d-2]
  real(r8),target,allocatable ::  CFOPE(:,:,:,:,:,:)                 !litter kinetic fraction, [-]
  real(r8),target,allocatable ::  ShootChemElmnts_pft(:,:,:,:)                    !canopy shoot element, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmnts_pft(:,:,:,:)                     !canopy leaf element, [g d-2]
  real(r8),target,allocatable ::  PetioleChemElmnts_pft(:,:,:,:)                    !canopy sheath element , [g d-2]
  real(r8),target,allocatable ::  StalkChemElmnts_pft(:,:,:,:)                    !canopy stalk element, [g d-2]
  real(r8),target,allocatable ::  CanopyStalkC_pft(:,:,:)                       !canopy active stalk C, [g d-2]
  real(r8),target,allocatable ::  ReserveChemElmnts_pft(:,:,:,:)                    !canopy reserve element, [g d-2]
  real(r8),target,allocatable ::  HuskChemElmnts_pft(:,:,:,:)                    !canopy husk element, [g d-2]
  real(r8),target,allocatable ::  EarChemElmnts_pft(:,:,:,:)                    !canopy ear element, [g d-2]
  real(r8),target,allocatable ::  GrainChemElmnts_pft(:,:,:,:)                     !canopy grain element, [g d-2]
  real(r8),target,allocatable ::  CanopyLeafShethC_pft(:,:,:)              !plant canopy leaf + sheath C, [gC d-2]
  real(r8),target,allocatable ::  CanopyLeafApft_lyr(:,:,:,:)                     !canopy layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CO2NetFix_pft(:,:,:)                        !canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  CanopyLeafCpft_lyr(:,:,:,:)                     !canopy layer leaf C, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstructElements_pft(:,:,:,:)   !canopy nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstructElementConc_pft(:,:,:,:)                    !canopy nonstructural element concentration, [g d-2]
  real(r8),target,allocatable ::  CanopyStemApft_lyr(:,:,:,:)                   !plant canopy layer stem area, [m2 d-2]
  real(r8),target,allocatable ::  NoduleNonstructElmnt_pft(:,:,:,:)                    !canopy nodule nonstructural element, [g d-2]
  real(r8),target,allocatable ::  StalkBiomassC_brch(:,:,:,:)                    !branch active stalk C, [g d-2]
  real(r8),target,allocatable ::  NonstructElmnt_brch(:,:,:,:,:)                   !branch nonstructural element, [g d-2]
  real(r8),target,allocatable ::  LeafPetolBiomassC_brch(:,:,:,:)           !plant branch leaf + sheath C, [g d-2]
  real(r8),target,allocatable ::  ShootChemElmnt_brch(:,:,:,:,:)                 !branch shoot C, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmnts_brch(:,:,:,:,:)                  !branch leaf element, [g d-2]
  real(r8),target,allocatable ::  PetoleChemElmnt_brch(:,:,:,:,:)                 !branch sheath element , [g d-2]
  real(r8),target,allocatable ::  StalkChemElmnts_brch(:,:,:,:,:)                  !branch stalk element, [g d-2]
  real(r8),target,allocatable ::  ReserveElmnts_brch(:,:,:,:,:)                  !branch reserve element, [g d-2]
  real(r8),target,allocatable ::  HuskChemElmnts_brch(:,:,:,:,:)                  !branch husk element, [g d-2]
  real(r8),target,allocatable ::  EarChemElmnts_brch(:,:,:,:,:)                 !branch ear element, [g d-2]
  real(r8),target,allocatable ::  GrainChemElmnts_brch(:,:,:,:,:)                  !branch grain element, [g d-2]
  real(r8),target,allocatable ::  LeafPetoNonstructElmntConc_brch(:,:,:,:,:)                    !branch nonstructural C concentration, [g d-2]
  real(r8),target,allocatable ::  NoduleNonstructElmnt_brch(:,:,:,:,:)                  !branch nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  CanopyNoduleChemElmnt_brch(:,:,:,:,:)                  !branch nodule element, [g d-2]
  real(r8),target,allocatable ::  PetioleChemElmntRemob_brch(:,:,:,:,:)                  !branch sheath structural element, [g d-2]
  real(r8),target,allocatable ::  BranchStalkChemElmnts_pft_pft(:,:,:,:,:)                    !branch stalk structural C, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmntRemob_brch(:,:,:,:,:)                     !branch leaf structural element, [g d-2]
  real(r8),target,allocatable ::  LeafElmntNode_brch(:,:,:,:,:,:)                    !leaf element, [g d-2]
  real(r8),target,allocatable ::  PetioleElmntNode_brch(:,:,:,:,:,:)                 !sheath element , [g d-2]
  real(r8),target,allocatable ::  InternodeChemElmnt_brch(:,:,:,:,:,:)                  !internode element, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmntByLayer_pft(:,:,:,:,:,:,:)                 !layer leaf element, [g d-2]
  real(r8),target,allocatable ::  CanopyLeafAreaByLayer_pft(:,:,:,:,:,:)                 !layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  LeafProteinCNode_brch(:,:,:,:,:)                    !layer leaf protein C, [g d-2]
  real(r8),target,allocatable ::  PetioleProteinCNode_brch(:,:,:,:,:)                   !layer sheath protein C, [g d-2]
  real(r8),target,allocatable ::  NoduleNonstructCconc_pft(:,:,:)                      !nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  GrainSeedBiomCMean_brch(:,:,:,:)                     !maximum grain C during grain fill, [g d-2]
  real(r8),target,allocatable ::  StandingDeadKCompChemElmnts_pft(:,:,:,:,:)                  !standing dead element fraction, [g d-2]
  real(r8),target,allocatable ::  StandingDeadChemElmnts_pft(:,:,:,:)                    !standing dead element, [g d-2]
  real(r8),target,allocatable ::  NonstructalElmnts_pft(:,:,:,:)                     !plant stored nonstructural element, [g d-2]
  real(r8),target,allocatable ::  SeedCPlanted_pft(:,:,:)                       !plant stored nonstructural C at planting, [g d-2]
  REAL(R8),target,allocatable ::  AvgCanopyBiomC2Graze_pft(:,:,:)                      !landscape average canopy shoot C, [g d-2]
  contains
!----------------------------------------------------------------------

  subroutine InitCanopyData

  implicit none
  allocate(CanopyPARalbedo_pft(JP,JY,JX));     CanopyPARalbedo_pft=0._r8
  allocate(TAUP(JP,JY,JX));     TAUP=0._r8
  allocate(CanopySWabsorpty_pft(JP,JY,JX));     CanopySWabsorpty_pft=0._r8
  allocate(CanopyPARabsorpty_pft(JP,JY,JX));     CanopyPARabsorpty_pft=0._r8
  allocate(RSMX(JP,JY,JX));     RSMX=0._r8
  allocate(CO2CuticleResist_pft(JP,JY,JX));     CO2CuticleResist_pft=0._r8
  allocate(MaxCanPStomaResistH2O_pft(JP,JY,JX));     MaxCanPStomaResistH2O_pft=0._r8
  allocate(RCS(JP,JY,JX));      RCS=0._r8
  allocate(CanPStomaResistH2O_pft(JP,JY,JX));       CanPStomaResistH2O_pft=0._r8
  allocate(MinCanPStomaResistH2O_pft(JP,JY,JX));     MinCanPStomaResistH2O_pft=0._r8
  allocate(BndlResistCanG(JY,JX));         BndlResistCanG=0._r8
  allocate(O2I(JP,JY,JX));      O2I=0._r8
  allocate(LeafIntracellularCO2_pft(JP,JY,JX));     LeafIntracellularCO2_pft=0._r8
  allocate(AirConc_pft(JP,JY,JX));     AirConc_pft=0._r8
  allocate(DiffCO2Atmos2Intracel_pft(JP,JY,JX));     DiffCO2Atmos2Intracel_pft=0._r8
  allocate(CanopyGasCO2_pft(JP,JY,JX));     CanopyGasCO2_pft=0._r8
  allocate(aquCO2Intraleaf_pft(JP,JY,JX));     aquCO2Intraleaf_pft=0._r8
  allocate(O2L(JP,JY,JX));      O2L=0._r8
  allocate(CO2Solubility_pft(JP,JY,JX));     CO2Solubility_pft=0._r8
  allocate(LeafO2Solubility_pft(JP,JY,JX));      LeafO2Solubility_pft=0._r8
  allocate(Km4LeafaqCO2_pft(JP,JY,JX));   Km4LeafaqCO2_pft=0._r8
  allocate(Km4RubiscoCarboxy_pft(JP,JY,JX));   Km4RubiscoCarboxy_pft=0._r8
  allocate(CHILL(JP,JY,JX));    CHILL=0._r8
  allocate(Vmax4RubiscoCarboxy_pft(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));Vmax4RubiscoCarboxy_pft=0._r8
  allocate(CO2lmtRubiscoCarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CO2lmtRubiscoCarboxyRate_node=0._r8
  allocate(CO2CompenPoint_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CO2CompenPoint_node=0._r8
  allocate(LigthSatCarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LigthSatCarboxyRate_node=0._r8
  allocate(RubiscoCarboxyEff_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));RubiscoCarboxyEff_node=0._r8
  allocate(CMassCO2BundleSheath_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CMassCO2BundleSheath_node=0._r8
  allocate(Vmax4PEPCarboxy_pft(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));Vmax4PEPCarboxy_pft=0._r8
  allocate(CO2lmtPEPCarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CO2lmtPEPCarboxyRate_node=0._r8
  allocate(LigthSatC4CarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LigthSatC4CarboxyRate_node=0._r8
  allocate(C4CarboxyEff_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));C4CarboxyEff_node=0._r8
  allocate(CPOOL4(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CPOOL4=0._r8
  allocate(CMassHCO3BundleSheath_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CMassHCO3BundleSheath_node=0._r8
  allocate(RubiscoActivity_brch(MaxNumBranches,JP,JY,JX));  RubiscoActivity_brch=0._r8
  allocate(NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));NutrientCtrlonC4Carboxy_node=0._r8
  allocate(C4PhotosynDowreg_brch(MaxNumBranches,JP,JY,JX)); C4PhotosynDowreg_brch=0._r8
  allocate(CNETX(JY,JX));       CNETX=0._r8
  allocate(VmaxRubCarboxyRef_pft(JP,JY,JX));     VmaxRubCarboxyRef_pft=0._r8
  allocate(VmaxRubOxyRef_pft(JP,JY,JX));     VmaxRubOxyRef_pft=0._r8
  allocate(VmaxPEPCarboxyRef_pft(JP,JY,JX));    VmaxPEPCarboxyRef_pft=0._r8
  allocate(XKCO2(JP,JY,JX));    XKCO2=0._r8
  allocate(XKO2(JP,JY,JX));     XKO2=0._r8
  allocate(Km4PEPCarboxy_pft(JP,JY,JX));   Km4PEPCarboxy_pft=0._r8
  allocate(LeafRuBPConc_pft(JP,JY,JX));     LeafRuBPConc_pft=0._r8
  allocate(FracLeafProtinAsPEPCarboxyl_pft(JP,JY,JX));     FracLeafProtinAsPEPCarboxyl_pft=0._r8
  allocate(SpecChloryfilAct_pft(JP,JY,JX));     SpecChloryfilAct_pft=0._r8
  allocate(LeafC3ChlorofilConc_pft(JP,JY,JX));      LeafC3ChlorofilConc_pft=0._r8
  allocate(LeafC4ChlorofilConc_pft(JP,JY,JX));     LeafC4ChlorofilConc_pft=0._r8
  allocate(CanPCi2CaRatio(JP,JY,JX));     CanPCi2CaRatio=0._r8
  allocate(RadNet2CanP(JP,JY,JX));     RadNet2CanP=0._r8
  allocate(LWRadCanP(JP,JY,JX));    LWRadCanP=0._r8
  allocate(RadSWbyCanopy_pft(JP,JY,JX));     RadSWbyCanopy_pft=0._r8
  allocate(RadPARbyCanopy_pft(JP,JY,JX));     RadPARbyCanopy_pft=0._r8
  allocate(FracRadPARbyCanopy_pft(JP,JY,JX));    FracRadPARbyCanopy_pft=0._r8
  allocate(TAU_RadThru(NumOfCanopyLayers+1,JY,JX));   TAU_RadThru=0._r8
  allocate(TAU_RadCapt(NumOfCanopyLayers+1,JY,JX));   TAU_RadCapt=0._r8
  allocate(FracSWRad2Grnd(JY,JX));       FracSWRad2Grnd=0._r8
  allocate(SWRadOnGrnd(JY,JX));        SWRadOnGrnd=0._r8
  allocate(LWRadCanGPrev(JY,JX));      LWRadCanGPrev=0._r8
  allocate(LWRadGrnd(JY,JX));      LWRadGrnd=0._r8
  allocate(CanH2OHeldVg(JY,JX));      CanH2OHeldVg=0._r8
  allocate(TFLWCI(JY,JX));      TFLWCI=0._r8
  allocate(PrecIntcptByCanG(JY,JX));       PrecIntcptByCanG=0._r8
  allocate(EvapTransHeatP(JP,JY,JX));    EvapTransHeatP=0._r8
  allocate(HeatXAir2PCan(JP,JY,JX));    HeatXAir2PCan=0._r8
  allocate(HeatStorCanP(JP,JY,JX));    HeatStorCanP=0._r8
  allocate(ENGYX(JP,JY,JX));    ENGYX=0._r8
  allocate(VHeatCapCanP(JP,JY,JX));    VHeatCapCanP=0._r8
  allocate(PSICanopy_pft(JP,JY,JX));    PSICanopy_pft=0._r8
  allocate(PSICanopyTurg_pft(JP,JY,JX));    PSICanopyTurg_pft=0._r8
  allocate(PSICanopyOsmo_pft(JP,JY,JX));    PSICanopyOsmo_pft=0._r8
  allocate(CanopyBndlResist_pft(JP,JY,JX));       CanopyBndlResist_pft=0._r8
  allocate(Transpiration_pft(JP,JY,JX));       Transpiration_pft=0._r8
  allocate(VapXAir2Canopy_pft(JP,JY,JX));    VapXAir2Canopy_pft=0._r8
  allocate(CanopyWater_pft(JP,JY,JX));    CanopyWater_pft=0._r8
  allocate(TEVAPP(JY,JX));      TEVAPP=0._r8
  allocate(VapXAir2CanG(JY,JX));      VapXAir2CanG=0._r8
  allocate(TENGYC(JY,JX));      TENGYC=0._r8
  allocate(THFLXC(JY,JX));      THFLXC=0._r8
  allocate(CanWatg(JY,JX));      CanWatg=0._r8
  allocate(LWRadCanG(JY,JX));       LWRadCanG=0._r8
  allocate(CanopySWAlbedo_pft(JP,JY,JX));     CanopySWAlbedo_pft=0._r8
  allocate(TAUR(JP,JY,JX));     TAUR=0._r8
  allocate(PrecIntcptByCanopy_pft(JP,JY,JX));     PrecIntcptByCanopy_pft=0._r8
  allocate(WatByPCanopy(JP,JY,JX));    WatByPCanopy=0._r8
  allocate(TKC(JP,JY,JX));      TKC=0._r8
  allocate(TCelciusCanopy_pft(JP,JY,JX));      TCelciusCanopy_pft=0._r8
  allocate(DTKC(JP,JY,JX));     DTKC=0._r8
  allocate(TKCanopy_pft(JP,JY,JX));     TKCanopy_pft=0._r8
  allocate(CPOOL3(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CPOOL3=0._r8
  allocate(NetCumElmntFlx2Plant_pft(NumPlantChemElmnts,JP,JY,JX));    NetCumElmntFlx2Plant_pft=0._r8
  allocate(WGLFT(NumOfCanopyLayers,JY,JX));    WGLFT=0._r8
  allocate(CFOPE(NumPlantChemElmnts,0:NumLitterGroups,jsken,JP,JY,JX));CFOPE=0._r8
  allocate(ShootChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX)); ShootChemElmnts_pft=0._r8
  allocate(LeafChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX));  LeafChemElmnts_pft=0._r8
  allocate(PetioleChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX)); PetioleChemElmnts_pft=0._r8
  allocate(StalkChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX)); StalkChemElmnts_pft=0._r8
  allocate(CanopyStalkC_pft(JP,JY,JX));    CanopyStalkC_pft=0._r8
  allocate(ReserveChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX));    ReserveChemElmnts_pft=0._r8
  allocate(HuskChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX));    HuskChemElmnts_pft=0._r8
  allocate(EarChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX));    EarChemElmnts_pft=0._r8
  allocate(GrainChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX));     GrainChemElmnts_pft=0._r8
  allocate(CanopyLeafShethC_pft(JP,JY,JX));     CanopyLeafShethC_pft=0._r8
  allocate(CanopyLeafApft_lyr(NumOfCanopyLayers,JP,JY,JX)); CanopyLeafApft_lyr=0._r8
  allocate(CO2NetFix_pft(JP,JY,JX));     CO2NetFix_pft=0._r8
  allocate(CanopyLeafCpft_lyr(NumOfCanopyLayers,JP,JY,JX)); CanopyLeafCpft_lyr=0._r8
  allocate(CanopyNonstructElements_pft(NumPlantChemElmnts,JP,JY,JX));   CanopyNonstructElements_pft=0._r8
  allocate(CanopyNonstructElementConc_pft(NumPlantChemElmnts,JP,JY,JX));   CanopyNonstructElementConc_pft=0._r8
  allocate(CanopyStemApft_lyr(NumOfCanopyLayers,JP,JY,JX)); CanopyStemApft_lyr=0._r8
  allocate(NoduleNonstructElmnt_pft(NumPlantChemElmnts,JP,JY,JX));   NoduleNonstructElmnt_pft=0._r8
  allocate(StalkBiomassC_brch(MaxNumBranches,JP,JY,JX));StalkBiomassC_brch=0._r8
  allocate(NonstructElmnt_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX)); NonstructElmnt_brch=0._r8
  allocate(LeafPetolBiomassC_brch(MaxNumBranches,JP,JY,JX)); LeafPetolBiomassC_brch=0._r8
  allocate(ShootChemElmnt_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));ShootChemElmnt_brch=0._r8
  allocate(LeafChemElmnts_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX)); LeafChemElmnts_brch=0._r8
  allocate(PetoleChemElmnt_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));PetoleChemElmnt_brch=0._r8
  allocate(StalkChemElmnts_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));StalkChemElmnts_brch=0._r8
  allocate(ReserveElmnts_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));ReserveElmnts_brch=0._r8
  allocate(HuskChemElmnts_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));HuskChemElmnts_brch=0._r8
  allocate(EarChemElmnts_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));EarChemElmnts_brch=0._r8
  allocate(GrainChemElmnts_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX)); GrainChemElmnts_brch=0._r8
  allocate(LeafPetoNonstructElmntConc_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));LeafPetoNonstructElmntConc_brch=0._r8
  allocate(NoduleNonstructElmnt_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));NoduleNonstructElmnt_brch=0._r8
  allocate(CanopyNoduleChemElmnt_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX)); CanopyNoduleChemElmnt_brch=0._r8
  allocate(PetioleChemElmntRemob_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));PetioleChemElmntRemob_brch=0._r8
  allocate(BranchStalkChemElmnts_pft_pft(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX));BranchStalkChemElmnts_pft_pft=0._r8
  allocate(LeafChemElmntRemob_brch(NumPlantChemElmnts,MaxNumBranches,JP,JY,JX)); LeafChemElmntRemob_brch=0._r8
  allocate(LeafElmntNode_brch(NumPlantChemElmnts,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafElmntNode_brch=0._r8
  allocate(PetioleElmntNode_brch(NumPlantChemElmnts,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleElmntNode_brch=0._r8
  allocate(InternodeChemElmnt_brch(NumPlantChemElmnts,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));InternodeChemElmnt_brch=0._r8
  allocate(LeafChemElmntByLayer_pft(NumPlantChemElmnts,NumOfCanopyLayers,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafChemElmntByLayer_pft=0._r8
  allocate(CanopyLeafAreaByLayer_pft(NumOfCanopyLayers,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CanopyLeafAreaByLayer_pft=0._r8
  allocate(LeafProteinCNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafProteinCNode_brch=0._r8
  allocate(PetioleProteinCNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleProteinCNode_brch=0._r8
  allocate(NoduleNonstructCconc_pft(JP,JY,JX));   NoduleNonstructCconc_pft=0._r8
  allocate(GrainSeedBiomCMean_brch(MaxNumBranches,JP,JY,JX)); GrainSeedBiomCMean_brch=0._r8
  allocate(StandingDeadKCompChemElmnts_pft(NumPlantChemElmnts,jsken,JP,JY,JX)); StandingDeadKCompChemElmnts_pft=0._r8
  allocate(StandingDeadChemElmnts_pft(NumPlantChemElmnts,JP,JY,JX));    StandingDeadChemElmnts_pft=0._r8
  allocate(NonstructalElmnts_pft(NumPlantChemElmnts,JP,JY,JX));  NonstructalElmnts_pft=0._r8
  allocate(SeedCPlanted_pft(JP,JY,JX));    SeedCPlanted_pft=0._r8
  allocate(AvgCanopyBiomC2Graze_pft(JP,JY,JX));   AvgCanopyBiomC2Graze_pft=0._r8
  end subroutine InitCanopyData

!----------------------------------------------------------------------
  subroutine DestructCanopyData
  use abortutils, only : destroy
  implicit none

  call destroy(CanopyPARalbedo_pft)
  call destroy(TAUP)
  call destroy(CanopySWabsorpty_pft)
  call destroy(CanopyPARabsorpty_pft)
  call destroy(RSMX)
  call destroy(CO2CuticleResist_pft)
  call destroy(MaxCanPStomaResistH2O_pft)
  call destroy(RCS)
  call destroy(CanPStomaResistH2O_pft)
  call destroy(MinCanPStomaResistH2O_pft)
  call destroy(BndlResistCanG)
  call destroy(O2I)
  call destroy(LeafIntracellularCO2_pft)
  call destroy(AirConc_pft)
  call destroy(DiffCO2Atmos2Intracel_pft)
  call destroy(CanopyGasCO2_pft)
  call destroy(aquCO2Intraleaf_pft)
  call destroy(O2L)
  call destroy(CO2Solubility_pft)
  call destroy(LeafO2Solubility_pft)
  call destroy(Km4LeafaqCO2_pft)
  call destroy(Km4RubiscoCarboxy_pft)
  call destroy(CHILL)
  call destroy(Vmax4RubiscoCarboxy_pft)
  call destroy(CO2lmtRubiscoCarboxyRate_node)
  call destroy(CO2CompenPoint_node)
  call destroy(LigthSatCarboxyRate_node)
  call destroy(RubiscoCarboxyEff_node)
  call destroy(CMassCO2BundleSheath_node)
  call destroy(Vmax4PEPCarboxy_pft)
  call destroy(CO2lmtPEPCarboxyRate_node)
  call destroy(LigthSatC4CarboxyRate_node)
  call destroy(C4CarboxyEff_node)
  call destroy(CPOOL4)
  call destroy(CMassHCO3BundleSheath_node)
  call destroy(RubiscoActivity_brch)
  call destroy(NutrientCtrlonC4Carboxy_node)
  call destroy(C4PhotosynDowreg_brch)
  call destroy(CNETX)
  call destroy(VmaxRubCarboxyRef_pft)
  call destroy(VmaxRubOxyRef_pft)
  call destroy(VmaxPEPCarboxyRef_pft)
  call destroy(XKCO2)
  call destroy(XKO2)
  call destroy(Km4PEPCarboxy_pft)
  call destroy(LeafRuBPConc_pft)
  call destroy(FracLeafProtinAsPEPCarboxyl_pft)
  call destroy(SpecChloryfilAct_pft)
  call destroy(LeafC3ChlorofilConc_pft)
  call destroy(LeafC4ChlorofilConc_pft)
  call destroy(CanPCi2CaRatio)
  call destroy(RadNet2CanP)
  call destroy(LWRadCanP)
  call destroy(RadSWbyCanopy_pft)
  call destroy(RadPARbyCanopy_pft)
  call destroy(FracRadPARbyCanopy_pft)
  call destroy(TAU_RadThru)
  call destroy(TAU_RadCapt)
  call destroy(FracSWRad2Grnd)
  call destroy(SWRadOnGrnd)
  call destroy(LWRadCanGPrev)
  call destroy(LWRadGrnd)
  call destroy(CanH2OHeldVg)
  call destroy(TFLWCI)
  call destroy(PrecIntcptByCanG)
  call destroy(EvapTransHeatP)
  call destroy(HeatXAir2PCan)
  call destroy(HeatStorCanP)
  call destroy(ENGYX)
  call destroy(VHeatCapCanP)
  call destroy(PSICanopy_pft)
  call destroy(PSICanopyTurg_pft)
  call destroy(PSICanopyOsmo_pft)
  call destroy(CanopyBndlResist_pft)
  call destroy(Transpiration_pft)
  call destroy(VapXAir2Canopy_pft)
  call destroy(CanopyWater_pft)
  call destroy(TEVAPP)
  call destroy(VapXAir2CanG)
  call destroy(TENGYC)
  call destroy(THFLXC)
  call destroy(CanWatg)
  call destroy(LWRadCanG)
  call destroy(CanopySWAlbedo_pft)
  call destroy(TAUR)
  call destroy(PrecIntcptByCanopy_pft)
  call destroy(WatByPCanopy)
  call destroy(TKC)
  call destroy(TCelciusCanopy_pft)
  call destroy(DTKC)
  call destroy(TKCanopy_pft)
  call destroy(CPOOL3)
  call destroy(NetCumElmntFlx2Plant_pft)
  call destroy(WGLFT)
  call destroy(CFOPE)
  call destroy(ShootChemElmnts_pft)
  call destroy(LeafChemElmnts_pft)
  call destroy(PetioleChemElmnts_pft)
  call destroy(StalkChemElmnts_pft)
  call destroy(CanopyStalkC_pft)
  call destroy(ReserveChemElmnts_pft)
  call destroy(HuskChemElmnts_pft)
  call destroy(EarChemElmnts_pft)
  call destroy(GrainChemElmnts_pft)
  call destroy(CanopyLeafShethC_pft)
  call destroy(CanopyLeafApft_lyr)
  call destroy(CO2NetFix_pft)
  call destroy(CanopyLeafCpft_lyr)
  call destroy(CanopyNonstructElements_pft)
  call destroy(CanopyNonstructElementConc_pft)
  call destroy(CanopyStemApft_lyr)
  call destroy(NoduleNonstructElmnt_pft)
  call destroy(StalkBiomassC_brch)
  call destroy(NonstructElmnt_brch)
  call destroy(LeafPetolBiomassC_brch)
  call destroy(ShootChemElmnt_brch)
  call destroy(LeafChemElmnts_brch)
  call destroy(PetoleChemElmnt_brch)
  call destroy(StalkChemElmnts_brch)
  call destroy(ReserveElmnts_brch)
  call destroy(HuskChemElmnts_brch)
  call destroy(EarChemElmnts_brch)
  call destroy(GrainChemElmnts_brch)
  call destroy(LeafPetoNonstructElmntConc_brch)
  call destroy(NoduleNonstructElmnt_brch)
  call destroy(CanopyNoduleChemElmnt_brch)
  call destroy(PetioleChemElmntRemob_brch)
  call destroy(BranchStalkChemElmnts_pft_pft)
  call destroy(LeafChemElmntRemob_brch)
  call destroy(LeafElmntNode_brch)
  call destroy(PetioleElmntNode_brch)
  call destroy(InternodeChemElmnt_brch)
  call destroy(LeafChemElmntByLayer_pft)
  call destroy(CanopyLeafAreaByLayer_pft)
  call destroy(LeafProteinCNode_brch)
  call destroy(PetioleProteinCNode_brch)
  call destroy(NoduleNonstructCconc_pft)
  call destroy(GrainSeedBiomCMean_brch)
  call destroy(StandingDeadKCompChemElmnts_pft)
  call destroy(StandingDeadChemElmnts_pft)
  call destroy(NonstructalElmnts_pft)
  call destroy(SeedCPlanted_pft)
  call destroy(AvgCanopyBiomC2Graze_pft)
  end subroutine DestructCanopyData

end module CanopyDataType
