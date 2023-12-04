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
  real(r8),target,allocatable ::  MaxCanPStomaResistH2O(:,:,:)      !maximum stomatal resistance to vapor, [s h-1]
  real(r8),target,allocatable ::  RCS(:,:,:)                         !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8),target,allocatable ::  CanPStomaResistH2O(:,:,:)         !canopy stomatal resistance, [h m-1]
  real(r8),target,allocatable ::  MinCanPStomaResistH2O(:,:,:)      !canopy minimum stomatal resistance, [s m-1]
  real(r8),target,allocatable ::  BndlResistCanG(:,:)                           !canopy boundary layer resistance, [m h-1]
  real(r8),target,allocatable ::  O2I(:,:,:)                         !leaf gaseous O2 concentration, [umol m-3]
  real(r8),target,allocatable ::  LeafIntracellularCO2_pft(:,:,:)                        !leaf gaseous CO2 concentration, [umol m-3]
  real(r8),target,allocatable ::  AirConc_pft(:,:,:)                        !total gas concentration, [mol m-3]
  real(r8),target,allocatable ::  DiffCO2Atmos2Intracel_pft(:,:,:)                        !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8),target,allocatable ::  CanopyGasCO2_pft(:,:,:)                        !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8),target,allocatable ::  CO2L(:,:,:)                        !leaf aqueous CO2 concentration, [uM]
  real(r8),target,allocatable ::  O2L(:,:,:)                         !leaf aqueous O2 concentration, [uM]
  real(r8),target,allocatable ::  CO2Solubility_pft(:,:,:)                        !leaf CO2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  SO2(:,:,:)                         !leaf O2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  XKCO2L(:,:,:)                      !leaf aqueous CO2 Km no O2, [uM]
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
  real(r8),target,allocatable ::  RubiscoActivity_brpft(:,:,:,:)                      !branch down-regulation of CO2 fixation, [-]
  real(r8),target,allocatable ::  NutrientCtrlonC4Carboxy_node(:,:,:,:,:)                   !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  C4PhotosynDowreg_brch(:,:,:,:)                     !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  CNETX(:,:)                         !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  VCMX(:,:,:)                        !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VOMX(:,:,:)                        !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VCMX4(:,:,:)                       !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  XKCO2(:,:,:)                       !Km for rubisco carboxylase activity, [uM]
  real(r8),target,allocatable ::  XKO2(:,:,:)                        !Km for rubisco oxygenase activity, [uM]
  real(r8),target,allocatable ::  Km4PEPCarboxy_pft(:,:,:)                      !Km for PEP carboxylase activity, [uM]
  real(r8),target,allocatable ::  RUBP(:,:,:)                        !leaf rubisco content, [g g-1]
  real(r8),target,allocatable ::  PEPC(:,:,:)                        !leaf PEP carboxylase content, [g g-1]
  real(r8),target,allocatable ::  ETMX(:,:,:)                        !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  CHL(:,:,:)                         !leaf C3 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  CHL4(:,:,:)                        !leaf C4 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  CanPCi2CaRatio(:,:,:)                        !Ci:Ca ratio, [-]
  real(r8),target,allocatable ::  RadNet2CanP(:,:,:)                 !canopy net radiation , [MJ d-2 h-1] >0
  real(r8),target,allocatable ::  LWRadCanP(:,:,:)                   !canopy longwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadSWbyCanopy_pft(:,:,:)                 !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadPARbyCanopy_pft(:,:,:)                   !canopy absorbed PAR , [umol m-2 s-1]
  real(r8),target,allocatable ::  FracRadPARbyCanopy_pft(:,:,:)                       !fraction of incoming PAR absorbed by canopy, [-]
  real(r8),target,allocatable ::  TAU0(:,:,:)                        !fraction of radiation transmitted by canopy layer, [-]
  real(r8),target,allocatable ::  TAUS(:,:,:)                        !fraction of radiation intercepted by canopy layer, [-]
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
  real(r8),target,allocatable ::  PSICanPTurg(:,:,:)                       !plant canopy turgor water potential, [Mpa]
  real(r8),target,allocatable ::  PSICanPOsmo(:,:,:)                 !platn canopy osmotic water potential, [Mpa]
  real(r8),target,allocatable ::  CanPbndlResist(:,:,:)                          !canopy boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  Transpiration_pft(:,:,:)                          !canopy transpiration, [m2 d-2 h-1]
  real(r8),target,allocatable ::  VapXAir2Canopy_pft(:,:,:)                !negative of canopy evaporation, [m2 d-2 h-1]
  real(r8),target,allocatable ::  CanWatP(:,:,:)                       !canopy water content associated with dry matter, [m3 d-2]
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
  real(r8),target,allocatable ::  TKCZ(:,:,:)                        !canopy temperature, [K]
  real(r8),target,allocatable ::  CPOOL3(:,:,:,:,:)                  !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8),target,allocatable ::  RSETE(:,:,:,:)                     !effect of canopy element status on seed set , []
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
  real(r8),target,allocatable ::  LeafPetioleBiomassC_brch(:,:,:,:)           !plant branch leaf + sheath C, [g d-2]
  real(r8),target,allocatable ::  ShootChemElmnt_brch(:,:,:,:,:)                 !branch shoot C, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmnts_brch(:,:,:,:,:)                  !branch leaf element, [g d-2]
  real(r8),target,allocatable ::  PetioleChemElmnts_brch(:,:,:,:,:)                 !branch sheath element , [g d-2]
  real(r8),target,allocatable ::  StalkChemElmnts_brch(:,:,:,:,:)                  !branch stalk element, [g d-2]
  real(r8),target,allocatable ::  ReserveChemElmnts_brch(:,:,:,:,:)                  !branch reserve element, [g d-2]
  real(r8),target,allocatable ::  HuskChemElmnts_brch(:,:,:,:,:)                  !branch husk element, [g d-2]
  real(r8),target,allocatable ::  EarChemElmnts_brch(:,:,:,:,:)                 !branch ear element, [g d-2]
  real(r8),target,allocatable ::  GrainChemElmnts_brch(:,:,:,:,:)                  !branch grain element, [g d-2]
  real(r8),target,allocatable ::  LeafPetioNonstructElmntConc_brch(:,:,:,:,:)                    !branch nonstructural C concentration, [g d-2]
  real(r8),target,allocatable ::  NoduleNonstructElmnt_brch(:,:,:,:,:)                  !branch nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  CanopyNoduleChemElmnt_brch(:,:,:,:,:)                  !branch nodule element, [g d-2]
  real(r8),target,allocatable ::  PetioleChemElmntRemob_brch(:,:,:,:,:)                  !branch sheath structural element, [g d-2]
  real(r8),target,allocatable ::  BranchStalkChemElmnts_pft_pft(:,:,:,:,:)                    !branch stalk structural C, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmntRemob_brch(:,:,:,:,:)                     !branch leaf structural element, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmntNode_brch(:,:,:,:,:,:)                    !leaf element, [g d-2]
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
  real(r8),target,allocatable ::  NonstructalChemElmnts_pft(:,:,:,:)                     !plant stored nonstructural element, [g d-2]
  real(r8),target,allocatable ::  SeedCPlanted_pft(:,:,:)                       !plant stored nonstructural C at planting, [g d-2]
  REAL(R8),target,allocatable ::  WTSHTA(:,:,:)                      !landscape average canopy shoot C, [g d-2]
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
  allocate(MaxCanPStomaResistH2O(JP,JY,JX));     MaxCanPStomaResistH2O=0._r8
  allocate(RCS(JP,JY,JX));      RCS=0._r8
  allocate(CanPStomaResistH2O(JP,JY,JX));       CanPStomaResistH2O=0._r8
  allocate(MinCanPStomaResistH2O(JP,JY,JX));     MinCanPStomaResistH2O=0._r8
  allocate(BndlResistCanG(JY,JX));         BndlResistCanG=0._r8
  allocate(O2I(JP,JY,JX));      O2I=0._r8
  allocate(LeafIntracellularCO2_pft(JP,JY,JX));     LeafIntracellularCO2_pft=0._r8
  allocate(AirConc_pft(JP,JY,JX));     AirConc_pft=0._r8
  allocate(DiffCO2Atmos2Intracel_pft(JP,JY,JX));     DiffCO2Atmos2Intracel_pft=0._r8
  allocate(CanopyGasCO2_pft(JP,JY,JX));     CanopyGasCO2_pft=0._r8
  allocate(CO2L(JP,JY,JX));     CO2L=0._r8
  allocate(O2L(JP,JY,JX));      O2L=0._r8
  allocate(CO2Solubility_pft(JP,JY,JX));     CO2Solubility_pft=0._r8
  allocate(SO2(JP,JY,JX));      SO2=0._r8
  allocate(XKCO2L(JP,JY,JX));   XKCO2L=0._r8
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
  allocate(RubiscoActivity_brpft(MaxNumBranches,JP,JY,JX));  RubiscoActivity_brpft=0._r8
  allocate(NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));NutrientCtrlonC4Carboxy_node=0._r8
  allocate(C4PhotosynDowreg_brch(MaxNumBranches,JP,JY,JX)); C4PhotosynDowreg_brch=0._r8
  allocate(CNETX(JY,JX));       CNETX=0._r8
  allocate(VCMX(JP,JY,JX));     VCMX=0._r8
  allocate(VOMX(JP,JY,JX));     VOMX=0._r8
  allocate(VCMX4(JP,JY,JX));    VCMX4=0._r8
  allocate(XKCO2(JP,JY,JX));    XKCO2=0._r8
  allocate(XKO2(JP,JY,JX));     XKO2=0._r8
  allocate(Km4PEPCarboxy_pft(JP,JY,JX));   Km4PEPCarboxy_pft=0._r8
  allocate(RUBP(JP,JY,JX));     RUBP=0._r8
  allocate(PEPC(JP,JY,JX));     PEPC=0._r8
  allocate(ETMX(JP,JY,JX));     ETMX=0._r8
  allocate(CHL(JP,JY,JX));      CHL=0._r8
  allocate(CHL4(JP,JY,JX));     CHL4=0._r8
  allocate(CanPCi2CaRatio(JP,JY,JX));     CanPCi2CaRatio=0._r8
  allocate(RadNet2CanP(JP,JY,JX));     RadNet2CanP=0._r8
  allocate(LWRadCanP(JP,JY,JX));    LWRadCanP=0._r8
  allocate(RadSWbyCanopy_pft(JP,JY,JX));     RadSWbyCanopy_pft=0._r8
  allocate(RadPARbyCanopy_pft(JP,JY,JX));     RadPARbyCanopy_pft=0._r8
  allocate(FracRadPARbyCanopy_pft(JP,JY,JX));    FracRadPARbyCanopy_pft=0._r8
  allocate(TAU0(JC+1,JY,JX));   TAU0=0._r8
  allocate(TAUS(JC+1,JY,JX));   TAUS=0._r8
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
  allocate(PSICanPTurg(JP,JY,JX));    PSICanPTurg=0._r8
  allocate(PSICanPOsmo(JP,JY,JX));    PSICanPOsmo=0._r8
  allocate(CanPbndlResist(JP,JY,JX));       CanPbndlResist=0._r8
  allocate(Transpiration_pft(JP,JY,JX));       Transpiration_pft=0._r8
  allocate(VapXAir2Canopy_pft(JP,JY,JX));    VapXAir2Canopy_pft=0._r8
  allocate(CanWatP(JP,JY,JX));    CanWatP=0._r8
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
  allocate(TKCZ(JP,JY,JX));     TKCZ=0._r8
  allocate(CPOOL3(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CPOOL3=0._r8
  allocate(RSETE(NumOfPlantChemElmnts,JP,JY,JX));    RSETE=0._r8
  allocate(WGLFT(JC,JY,JX));    WGLFT=0._r8
  allocate(CFOPE(NumOfPlantChemElmnts,0:NumLitterGroups,jsken,JP,JY,JX));CFOPE=0._r8
  allocate(ShootChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX)); ShootChemElmnts_pft=0._r8
  allocate(LeafChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));  LeafChemElmnts_pft=0._r8
  allocate(PetioleChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX)); PetioleChemElmnts_pft=0._r8
  allocate(StalkChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX)); StalkChemElmnts_pft=0._r8
  allocate(CanopyStalkC_pft(JP,JY,JX));    CanopyStalkC_pft=0._r8
  allocate(ReserveChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));    ReserveChemElmnts_pft=0._r8
  allocate(HuskChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));    HuskChemElmnts_pft=0._r8
  allocate(EarChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));    EarChemElmnts_pft=0._r8
  allocate(GrainChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));     GrainChemElmnts_pft=0._r8
  allocate(CanopyLeafShethC_pft(JP,JY,JX));     CanopyLeafShethC_pft=0._r8
  allocate(CanopyLeafApft_lyr(JC,JP,JY,JX)); CanopyLeafApft_lyr=0._r8
  allocate(CO2NetFix_pft(JP,JY,JX));     CO2NetFix_pft=0._r8
  allocate(CanopyLeafCpft_lyr(JC,JP,JY,JX)); CanopyLeafCpft_lyr=0._r8
  allocate(CanopyNonstructElements_pft(NumOfPlantChemElmnts,JP,JY,JX));   CanopyNonstructElements_pft=0._r8
  allocate(CanopyNonstructElementConc_pft(NumOfPlantChemElmnts,JP,JY,JX));   CanopyNonstructElementConc_pft=0._r8
  allocate(CanopyStemApft_lyr(JC,JP,JY,JX)); CanopyStemApft_lyr=0._r8
  allocate(NoduleNonstructElmnt_pft(NumOfPlantChemElmnts,JP,JY,JX));   NoduleNonstructElmnt_pft=0._r8
  allocate(StalkBiomassC_brch(MaxNumBranches,JP,JY,JX));StalkBiomassC_brch=0._r8
  allocate(NonstructElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX)); NonstructElmnt_brch=0._r8
  allocate(LeafPetioleBiomassC_brch(MaxNumBranches,JP,JY,JX)); LeafPetioleBiomassC_brch=0._r8
  allocate(ShootChemElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));ShootChemElmnt_brch=0._r8
  allocate(LeafChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX)); LeafChemElmnts_brch=0._r8
  allocate(PetioleChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));PetioleChemElmnts_brch=0._r8
  allocate(StalkChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));StalkChemElmnts_brch=0._r8
  allocate(ReserveChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));ReserveChemElmnts_brch=0._r8
  allocate(HuskChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));HuskChemElmnts_brch=0._r8
  allocate(EarChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));EarChemElmnts_brch=0._r8
  allocate(GrainChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX)); GrainChemElmnts_brch=0._r8
  allocate(LeafPetioNonstructElmntConc_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));LeafPetioNonstructElmntConc_brch=0._r8
  allocate(NoduleNonstructElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));NoduleNonstructElmnt_brch=0._r8
  allocate(CanopyNoduleChemElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX)); CanopyNoduleChemElmnt_brch=0._r8
  allocate(PetioleChemElmntRemob_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));PetioleChemElmntRemob_brch=0._r8
  allocate(BranchStalkChemElmnts_pft_pft(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX));BranchStalkChemElmnts_pft_pft=0._r8
  allocate(LeafChemElmntRemob_brch(NumOfPlantChemElmnts,MaxNumBranches,JP,JY,JX)); LeafChemElmntRemob_brch=0._r8
  allocate(LeafChemElmntNode_brch(NumOfPlantChemElmnts,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafChemElmntNode_brch=0._r8
  allocate(PetioleElmntNode_brch(NumOfPlantChemElmnts,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleElmntNode_brch=0._r8
  allocate(InternodeChemElmnt_brch(NumOfPlantChemElmnts,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));InternodeChemElmnt_brch=0._r8
  allocate(LeafChemElmntByLayer_pft(NumOfPlantChemElmnts,JC,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafChemElmntByLayer_pft=0._r8
  allocate(CanopyLeafAreaByLayer_pft(JC,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CanopyLeafAreaByLayer_pft=0._r8
  allocate(LeafProteinCNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafProteinCNode_brch=0._r8
  allocate(PetioleProteinCNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleProteinCNode_brch=0._r8
  allocate(NoduleNonstructCconc_pft(JP,JY,JX));   NoduleNonstructCconc_pft=0._r8
  allocate(GrainSeedBiomCMean_brch(MaxNumBranches,JP,JY,JX)); GrainSeedBiomCMean_brch=0._r8
  allocate(StandingDeadKCompChemElmnts_pft(NumOfPlantChemElmnts,jsken,JP,JY,JX)); StandingDeadKCompChemElmnts_pft=0._r8
  allocate(StandingDeadChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));    StandingDeadChemElmnts_pft=0._r8
  allocate(NonstructalChemElmnts_pft(NumOfPlantChemElmnts,JP,JY,JX));  NonstructalChemElmnts_pft=0._r8
  allocate(SeedCPlanted_pft(JP,JY,JX));    SeedCPlanted_pft=0._r8
  allocate(WTSHTA(JP,JY,JX));   WTSHTA=0._r8
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
  call destroy(MaxCanPStomaResistH2O)
  call destroy(RCS)
  call destroy(CanPStomaResistH2O)
  call destroy(MinCanPStomaResistH2O)
  call destroy(BndlResistCanG)
  call destroy(O2I)
  call destroy(LeafIntracellularCO2_pft)
  call destroy(AirConc_pft)
  call destroy(DiffCO2Atmos2Intracel_pft)
  call destroy(CanopyGasCO2_pft)
  call destroy(CO2L)
  call destroy(O2L)
  call destroy(CO2Solubility_pft)
  call destroy(SO2)
  call destroy(XKCO2L)
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
  call destroy(RubiscoActivity_brpft)
  call destroy(NutrientCtrlonC4Carboxy_node)
  call destroy(C4PhotosynDowreg_brch)
  call destroy(CNETX)
  call destroy(VCMX)
  call destroy(VOMX)
  call destroy(VCMX4)
  call destroy(XKCO2)
  call destroy(XKO2)
  call destroy(Km4PEPCarboxy_pft)
  call destroy(RUBP)
  call destroy(PEPC)
  call destroy(ETMX)
  call destroy(CHL)
  call destroy(CHL4)
  call destroy(CanPCi2CaRatio)
  call destroy(RadNet2CanP)
  call destroy(LWRadCanP)
  call destroy(RadSWbyCanopy_pft)
  call destroy(RadPARbyCanopy_pft)
  call destroy(FracRadPARbyCanopy_pft)
  call destroy(TAU0)
  call destroy(TAUS)
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
  call destroy(PSICanPTurg)
  call destroy(PSICanPOsmo)
  call destroy(CanPbndlResist)
  call destroy(Transpiration_pft)
  call destroy(VapXAir2Canopy_pft)
  call destroy(CanWatP)
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
  call destroy(TKCZ)
  call destroy(CPOOL3)
  call destroy(RSETE)
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
  call destroy(LeafPetioleBiomassC_brch)
  call destroy(ShootChemElmnt_brch)
  call destroy(LeafChemElmnts_brch)
  call destroy(PetioleChemElmnts_brch)
  call destroy(StalkChemElmnts_brch)
  call destroy(ReserveChemElmnts_brch)
  call destroy(HuskChemElmnts_brch)
  call destroy(EarChemElmnts_brch)
  call destroy(GrainChemElmnts_brch)
  call destroy(LeafPetioNonstructElmntConc_brch)
  call destroy(NoduleNonstructElmnt_brch)
  call destroy(CanopyNoduleChemElmnt_brch)
  call destroy(PetioleChemElmntRemob_brch)
  call destroy(BranchStalkChemElmnts_pft_pft)
  call destroy(LeafChemElmntRemob_brch)
  call destroy(LeafChemElmntNode_brch)
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
  call destroy(NonstructalChemElmnts_pft)
  call destroy(SeedCPlanted_pft)
  call destroy(WTSHTA)
  end subroutine DestructCanopyData

end module CanopyDataType
