module CanopyDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL     
  use GridConsts
  use ElmIDMod
  use EcoSIMConfig, only : jsken => jskenc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  canopy_growth_pft(:,:,:)                !canopy structural growth rate [gC/h]
  real(r8),target,allocatable ::  StomatalStress_pft(:,:,:)               !stomatal stress from water/turgor, [0-1]
  real(r8),target,allocatable ::  CanopyPARalbedo_pft(:,:,:)                        !canopy PAR albedo , [-]
  real(r8),target,allocatable ::  RadPARLeafTransmis_pft(:,:,:)                        !canopy PAR transmissivity , [-]
  real(r8),target,allocatable ::  LeafSWabsorpty_pft(:,:,:)                        !canopy shortwave absorptivity , [-]
  real(r8),target,allocatable ::  LeafPARabsorpty_pft(:,:,:)                        !canopy PAR absorptivity , [-]
  real(r8),target,allocatable ::  CuticleResist_pft(:,:,:)                        !maximum stomatal resistance to vapor, [s m-1]
  real(r8),target,allocatable ::  CO2CuticleResist_pft(:,:,:)                        !maximum stomatal resistance to CO2, [s h-1]
  real(r8),target,allocatable ::  H2OCuticleResist_pft(:,:,:)      !maximum stomatal resistance to vapor, [s h-1]
  real(r8),target,allocatable ::  RCS(:,:,:)                         !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8),target,allocatable ::  CanPStomaResistH2O_pft(:,:,:)         !canopy stomatal resistance, [h m-1]
  real(r8),target,allocatable ::  MinCanPStomaResistH2O_pft(:,:,:)      !canopy minimum stomatal resistance, [s m-1]
  real(r8),target,allocatable ::  CanopyBndlResist_col(:,:)                           !canopy boundary layer resistance, [m h-1]
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
  real(r8),target,allocatable ::  ChillHours_pft(:,:,:)                       !chilling effect on CO2 fixation, [-]
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
  real(r8),target,allocatable ::  CPOOL4_node(:,:,:,:,:)                  !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  CMassHCO3BundleSheath_node(:,:,:,:,:)                    !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  RubiscoActivity_brch(:,:,:,:)                      !branch down-regulation of CO2 fixation, [-]
  real(r8),target,allocatable ::  NutrientCtrlonC4Carboxy_node(:,:,:,:,:)                   !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  C4PhotosynDowreg_brch(:,:,:,:)                     !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  NetCO2Flx2Canopy_col(:,:)                         !total net canopy CO2 exchange, [g d-2 h-1]
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
  real(r8),target,allocatable ::  RadNet2Canopy_pft(:,:,:)                 !canopy net radiation , [MJ d-2 h-1] >0
  real(r8),target,allocatable ::  LWRadCanopy_pft(:,:,:)                   !canopy longwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadSWbyCanopy_pft(:,:,:)                 !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadPARbyCanopy_pft(:,:,:)                   !canopy absorbed PAR , [umol m-2 s-1]
  real(r8),target,allocatable ::  FracPARads2Canopy_pft(:,:,:)                       !fraction of incoming PAR absorbed by canopy, [-]
  real(r8),target,allocatable ::  TAU_RadThru(:,:,:)                        !fraction of radiation transmitted by canopy layer, [-]
  real(r8),target,allocatable ::  TAU_DirRadTransm(:,:,:)                        !fraction of radiation intercepted by canopy layer, [-]
  real(r8),target,allocatable ::  FracSWRad2Grnd_col(:,:)                         !fraction of radiation intercepted by ground surface, [-]
  real(r8),target,allocatable ::  RadSWGrnd_col(:,:)                        !shortwave radiation incident on ground surface, [MJ h-1]
  real(r8),target,allocatable ::  LWRadCanGPrev_col(:,:)                        !longwave radiation emitted by canopy, [MJ h-1]
  real(r8),target,allocatable ::  LWRadGrnd(:,:)                            !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  CanH2OHeldVg_col(:,:)                     !canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  Prec2Canopy_col(:,:)                             !net ice transfer to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  PrecIntceptByCanopy_col(:,:)            !grid net precipitation water interception to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  EvapTransHeat_pft(:,:,:)                       !canopy latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatXAir2PCan_pft(:,:,:)               !air to canopy sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatStorCanopy_pft(:,:,:)                       !canopy storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  ENGYX_pft(:,:,:)                       !canopy heat storage from previous time step, [MJ d-2]
  real(r8),target,allocatable ::  VHeatCapCanP_pft(:,:,:)                       !canopy heat capacity, [MJ d-2 K-1]
  real(r8),target,allocatable ::  PSICanopy_pft(:,:,:)                     !plant canopy total water potential , [Mpa]
  real(r8),target,allocatable ::  PSICanopyTurg_pft(:,:,:)                       !plant canopy turgor water potential, [Mpa]
  real(r8),target,allocatable ::  PSICanopyOsmo_pft(:,:,:)                 !platn canopy osmotic water potential, [Mpa]
  real(r8),target,allocatable ::  CanopyBndlResist_pft(:,:,:)                          !canopy boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  Transpiration_pft(:,:,:)                          !canopy transpiration, [m2 d-2 h-1]
  real(r8),target,allocatable ::  VapXAir2Canopy_pft(:,:,:)                !negative of canopy evaporation, [m2 d-2 h-1]
  real(r8),target,allocatable ::  CanopyWater_pft(:,:,:)                       !canopy water content associated with dry matter, [m3 d-2]
  real(r8),target,allocatable ::  QvET_col(:,:)                        !total canopy evaporation + transpiration, [m3 d-2]
  real(r8),target,allocatable ::  VapXAir2Canopy_col(:,:)                        !total canopy evaporation, [m3 d-2]
  real(r8),target,allocatable ::  CanopyHeatStor_col(:,:)                        !total canopy heat content, [MJ  d-2]
  real(r8),target,allocatable ::  HeatFlx2Canopy_col(:,:)                        !total canopy heat flux, [MJ  d-2]
  real(r8),target,allocatable ::  CanWat_col(:,:)                       !total canopy water content stored in dry matter, [m3 d-2]
  real(r8),target,allocatable ::  LWRadCanG(:,:)                         !total canopy LW emission, [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadSWLeafAlbedo_pft(:,:,:)                        !canopy shortwave albedo , [-]
  real(r8),target,allocatable ::  RadSWLeafTransmis_pft(:,:,:)                        !canopy shortwave transmissivity , [-]
  real(r8),target,allocatable ::  PrecIntcptByCanopy_pft(:,:,:)                        !water flux into plant canopy, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WatByPCanopy_pft(:,:,:)                   !canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  TKC(:,:,:)                         !canopy temperature, [K]
  real(r8),target,allocatable ::  TCelciusCanopy_pft(:,:,:)                         !canopy temperature, [oC]
  real(r8),target,allocatable ::  DeltaTKC_pft(:,:,:)                        !change in canopy temperature, [K]
  real(r8),target,allocatable ::  TKCanopy_pft(:,:,:)                        !canopy temperature, [K]
  real(r8),target,allocatable ::  CPOOL3_node(:,:,:,:,:)                  !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8),target,allocatable ::  NetCumElmntFlx2Plant_pft(:,:,:,:)                     !effect of canopy element status on seed set , []
  real(r8),target,allocatable ::  tCanLeafC_cl(:,:,:)                       !total leaf mass, [g d-2]
  real(r8),target,allocatable ::  ElmAllocmat4Litr(:,:,:,:,:,:)                 !litter kinetic fraction, [-]
  real(r8),target,allocatable ::  ShootElms_pft(:,:,:,:)                 !shoot structural element, [g d-2]
  real(r8),target,allocatable ::  ShootC4NonstC_brch(:,:,:,:)
  real(r8),target,allocatable ::  ShootStrutElms_pft(:,:,:,:)                    !canopy shoot element, [g d-2]
  real(r8),target,allocatable ::  LeafStrutElms_pft(:,:,:,:)                     !canopy leaf element, [g d-2]
  real(r8),target,allocatable ::  PetoleStrutElms_pft(:,:,:,:)                    !canopy sheath element , [g d-2]
  real(r8),target,allocatable ::  StalkStrutElms_pft(:,:,:,:)                    !canopy stalk element, [g d-2]
  real(r8),target,allocatable ::  CanopyStalkC_pft(:,:,:)                       !canopy active stalk C, [g d-2]
  real(r8),target,allocatable ::  StalkRsrvElms_pft(:,:,:,:)                    !canopy reserve element, [g d-2]
  real(r8),target,allocatable ::  HuskStrutElms_pft(:,:,:,:)                    !canopy husk element, [g d-2]
  real(r8),target,allocatable ::  EarStrutElms_pft(:,:,:,:)                    !canopy ear element, [g d-2]
  real(r8),target,allocatable ::  GrainStrutElms_pft(:,:,:,:)                     !canopy grain element, [g d-2]
  real(r8),target,allocatable ::  CanopyLeafShethC_pft(:,:,:)              !plant canopy leaf + sheath C, [gC d-2]
  real(r8),target,allocatable ::  CanopyLeafAreaZ_pft(:,:,:,:)                     !canopy layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CO2NetFix_pft(:,:,:)                        !canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  CanopyLeafCLyr_pft(:,:,:,:)                     !canopy layer leaf C, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstElms_pft(:,:,:,:)   !canopy nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstElmConc_pft(:,:,:,:)                    !canopy nonstructural element concentration, [g d-2]
  real(r8),target,allocatable ::  CanopyStemAreaZ_pft(:,:,:,:)                   !plant canopy layer stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyNodulNonstElms_pft(:,:,:,:)                    !canopy nodule nonstructural element, [g d-2]
  real(r8),target,allocatable ::  CanopyNodulElms_pft(:,:,:,:)                  !canopy nodule elemental biomass [g d-2]
  real(r8),target,allocatable ::  StalkBiomassC_brch(:,:,:,:)                    !branch active stalk C, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstElms_brch(:,:,:,:,:)                   !branch nonstructural element, [g d-2]
  real(r8),target,allocatable ::  LeafPetolBiomassC_brch(:,:,:,:)           !plant branch leaf + sheath C, [g d-2]
  real(r8),target,allocatable ::  ShootStrutElms_brch(:,:,:,:,:)                 !branch shoot C, [g d-2]
  real(r8),target,allocatable ::  LeafStrutElms_brch(:,:,:,:,:)                  !branch leaf element, [g d-2]
  real(r8),target,allocatable ::  PetoleStrutElms_brch(:,:,:,:,:)                 !branch sheath element , [g d-2]
  real(r8),target,allocatable ::  StalkStrutElms_brch(:,:,:,:,:)                  !branch stalk element, [g d-2]
  real(r8),target,allocatable ::  StalkRsrvElms_brch(:,:,:,:,:)                  !branch reserve element, [g d-2]
  real(r8),target,allocatable ::  HuskStrutElms_brch(:,:,:,:,:)                  !branch husk element, [g d-2]
  real(r8),target,allocatable ::  EarStrutElms_brch(:,:,:,:,:)                 !branch ear element, [g d-2]
  real(r8),target,allocatable ::  GrainStrutElms_brch(:,:,:,:,:)                  !branch grain element, [g d-2]
  real(r8),target,allocatable ::  LeafPetoNonstElmConc_brch(:,:,:,:,:)                    !branch nonstructural C concentration, [g d-2]
  real(r8),target,allocatable ::  CanopyNodulNonstElms_brch(:,:,:,:,:)                  !branch nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  CanopyNodulStrutElms_brch(:,:,:,:,:)                  !branch nodule element, [g d-2]
  real(r8),target,allocatable ::  PetioleChemElmRemob_brch(:,:,:,:,:)                  !branch sheath structural element, [g d-2]
  real(r8),target,allocatable ::  SenecStalkStrutElms_brch(:,:,:,:,:)                    !branch stalk structural C, [g d-2]
  real(r8),target,allocatable ::  LeafChemElmRemob_brch(:,:,:,:,:)                     !branch leaf structural element, [g d-2]
  real(r8),target,allocatable ::  LeafElmntNode_brch(:,:,:,:,:,:)                    !leaf element, [g d-2]
  real(r8),target,allocatable ::  PetioleElmntNode_brch(:,:,:,:,:,:)                 !sheath element , [g d-2]
  real(r8),target,allocatable ::  InternodeStrutElms_brch(:,:,:,:,:,:)                  !internode element, [g d-2]
  real(r8),target,allocatable ::  LeafElmsByLayerNode_brch(:,:,:,:,:,:,:)                 !layer leaf element, [g d-2]
  real(r8),target,allocatable ::  CanopyLeafArea_lpft(:,:,:,:,:,:)                 !layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  LeafProteinCNode_brch(:,:,:,:,:)                    !layer leaf protein C, [g d-2]
  real(r8),target,allocatable ::  PetoleProteinCNode_brch(:,:,:,:,:)                   !layer sheath protein C, [g d-2]
  real(r8),target,allocatable ::  NoduleNonstructCconc_pft(:,:,:)                      !nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  GrainSeedBiomCMean_brch(:,:,:,:)                     !maximum grain C during grain fill, [g d-2]
  real(r8),target,allocatable ::  StandDeadKCompElms_pft(:,:,:,:,:)                  !standing dead element fraction, [g d-2]
  real(r8),target,allocatable ::  StandDeadStrutElms_pft(:,:,:,:)                    !standing dead element, [g d-2]
  real(r8),target,allocatable ::  SeasonalNonstElms_pft(:,:,:,:)                     !plant stored nonstructural element, [g d-2]
  real(r8),target,allocatable ::  SeedCPlanted_pft(:,:,:)                       !plant stored nonstructural C at planting, [g d-2]
  REAL(R8),target,allocatable ::  AvgCanopyBiomC2Graze_pft(:,:,:)                      !landscape average canopy shoot C, [g d-2]
  real(r8),target,allocatable :: CO2FixCL_pft(:,:,:)
  real(r8),target,allocatable :: CO2FixLL_pft(:,:,:)
  contains
!----------------------------------------------------------------------

  subroutine InitCanopyData

  implicit none
  allocate(CO2FixCL_pft(JP,JY,JX)); CO2FixCL_pft=spval
  allocate(CO2FixLL_pft(JP,JY,JX)); CO2FixLL_pft=spval
  allocate(canopy_growth_pft(JP,JY,JX)); canopy_growth_pft=spval
  allocate(StomatalStress_pft(JP,JY,JX)); StomatalStress_pft=spval     !no stress when equals to one
  allocate(CanopyPARalbedo_pft(JP,JY,JX));     CanopyPARalbedo_pft=0._r8
  allocate(RadPARLeafTransmis_pft(JP,JY,JX));     RadPARLeafTransmis_pft=0._r8
  allocate(LeafSWabsorpty_pft(JP,JY,JX));     LeafSWabsorpty_pft=0._r8
  allocate(LeafPARabsorpty_pft(JP,JY,JX));     LeafPARabsorpty_pft=0._r8
  allocate(CuticleResist_pft(JP,JY,JX));     CuticleResist_pft=0._r8
  allocate(CO2CuticleResist_pft(JP,JY,JX));     CO2CuticleResist_pft=0._r8
  allocate(H2OCuticleResist_pft(JP,JY,JX));     H2OCuticleResist_pft=0._r8
  allocate(RCS(JP,JY,JX));      RCS=0._r8
  allocate(CanPStomaResistH2O_pft(JP,JY,JX));       CanPStomaResistH2O_pft=0._r8
  allocate(MinCanPStomaResistH2O_pft(JP,JY,JX));     MinCanPStomaResistH2O_pft=0._r8
  allocate(CanopyBndlResist_col(JY,JX));         CanopyBndlResist_col=0._r8
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
  allocate(ChillHours_pft(JP,JY,JX));    ChillHours_pft=0._r8
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
  allocate(CPOOL4_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CPOOL4_node=0._r8
  allocate(CMassHCO3BundleSheath_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CMassHCO3BundleSheath_node=0._r8
  allocate(RubiscoActivity_brch(MaxNumBranches,JP,JY,JX));  RubiscoActivity_brch=0._r8
  allocate(NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));NutrientCtrlonC4Carboxy_node=0._r8
  allocate(C4PhotosynDowreg_brch(MaxNumBranches,JP,JY,JX)); C4PhotosynDowreg_brch=0._r8
  allocate(NetCO2Flx2Canopy_col(JY,JX));       NetCO2Flx2Canopy_col=0._r8
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
  allocate(RadNet2Canopy_pft(JP,JY,JX));     RadNet2Canopy_pft=0._r8
  allocate(LWRadCanopy_pft(JP,JY,JX));    LWRadCanopy_pft=0._r8
  allocate(RadSWbyCanopy_pft(JP,JY,JX));     RadSWbyCanopy_pft=0._r8
  allocate(RadPARbyCanopy_pft(JP,JY,JX));     RadPARbyCanopy_pft=0._r8
  allocate(FracPARads2Canopy_pft(JP,JY,JX));    FracPARads2Canopy_pft=0._r8
  allocate(TAU_RadThru(NumOfCanopyLayers+1,JY,JX));   TAU_RadThru=0._r8
  allocate(TAU_DirRadTransm(NumOfCanopyLayers+1,JY,JX));   TAU_DirRadTransm=0._r8
  allocate(FracSWRad2Grnd_col(JY,JX));       FracSWRad2Grnd_col=0._r8
  allocate(RadSWGrnd_col(JY,JX));        RadSWGrnd_col=0._r8
  allocate(LWRadCanGPrev_col(JY,JX));      LWRadCanGPrev_col=0._r8
  allocate(LWRadGrnd(JY,JX));      LWRadGrnd=0._r8
  allocate(CanH2OHeldVg_col(JY,JX));      CanH2OHeldVg_col=0._r8
  allocate(Prec2Canopy_col(JY,JX));      Prec2Canopy_col=0._r8
  allocate(PrecIntceptByCanopy_col(JY,JX));       PrecIntceptByCanopy_col=0._r8
  allocate(EvapTransHeat_pft(JP,JY,JX));    EvapTransHeat_pft=0._r8
  allocate(HeatXAir2PCan_pft(JP,JY,JX));    HeatXAir2PCan_pft=0._r8
  allocate(HeatStorCanopy_pft(JP,JY,JX));    HeatStorCanopy_pft=0._r8
  allocate(ENGYX_pft(JP,JY,JX));    ENGYX_pft=0._r8
  allocate(VHeatCapCanP_pft(JP,JY,JX));    VHeatCapCanP_pft=0._r8
  allocate(PSICanopy_pft(JP,JY,JX));    PSICanopy_pft=0._r8
  allocate(PSICanopyTurg_pft(JP,JY,JX));    PSICanopyTurg_pft=0._r8
  allocate(PSICanopyOsmo_pft(JP,JY,JX));    PSICanopyOsmo_pft=0._r8
  allocate(CanopyBndlResist_pft(JP,JY,JX));       CanopyBndlResist_pft=0._r8
  allocate(Transpiration_pft(JP,JY,JX));       Transpiration_pft=0._r8
  allocate(VapXAir2Canopy_pft(JP,JY,JX));    VapXAir2Canopy_pft=0._r8
  allocate(CanopyWater_pft(JP,JY,JX));    CanopyWater_pft=0._r8
  allocate(QvET_col(JY,JX));      QvET_col=0._r8
  allocate(VapXAir2Canopy_col(JY,JX));      VapXAir2Canopy_col=0._r8
  allocate(CanopyHeatStor_col(JY,JX));      CanopyHeatStor_col=0._r8
  allocate(HeatFlx2Canopy_col(JY,JX));      HeatFlx2Canopy_col=0._r8
  allocate(CanWat_col(JY,JX));      CanWat_col=0._r8
  allocate(LWRadCanG(JY,JX));       LWRadCanG=0._r8
  allocate(RadSWLeafAlbedo_pft(JP,JY,JX));     RadSWLeafAlbedo_pft=0._r8
  allocate(RadSWLeafTransmis_pft(JP,JY,JX));     RadSWLeafTransmis_pft=0._r8
  allocate(PrecIntcptByCanopy_pft(JP,JY,JX));     PrecIntcptByCanopy_pft=0._r8
  allocate(WatByPCanopy_pft(JP,JY,JX));    WatByPCanopy_pft=0._r8
  allocate(TKC(JP,JY,JX));      TKC=0._r8
  allocate(TCelciusCanopy_pft(JP,JY,JX));      TCelciusCanopy_pft=0._r8
  allocate(DeltaTKC_pft(JP,JY,JX));     DeltaTKC_pft=0._r8
  allocate(TKCanopy_pft(JP,JY,JX));     TKCanopy_pft=0._r8
  allocate(CPOOL3_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CPOOL3_node=0._r8
  allocate(NetCumElmntFlx2Plant_pft(NumPlantChemElms,JP,JY,JX));    NetCumElmntFlx2Plant_pft=0._r8
  allocate(tCanLeafC_cl(NumOfCanopyLayers,JY,JX));    tCanLeafC_cl=0._r8
  allocate(ElmAllocmat4Litr(NumPlantChemElms,0:NumLitterGroups,jsken,JP,JY,JX));ElmAllocmat4Litr=0._r8
  allocate(ShootStrutElms_pft(NumPlantChemElms,JP,JY,JX)); ShootStrutElms_pft=0._r8
  allocate(LeafStrutElms_pft(NumPlantChemElms,JP,JY,JX));  LeafStrutElms_pft=0._r8
  allocate(PetoleStrutElms_pft(NumPlantChemElms,JP,JY,JX)); PetoleStrutElms_pft=0._r8
  allocate(StalkStrutElms_pft(NumPlantChemElms,JP,JY,JX)); StalkStrutElms_pft=0._r8
  allocate(CanopyStalkC_pft(JP,JY,JX));    CanopyStalkC_pft=0._r8
  allocate(StalkRsrvElms_pft(NumPlantChemElms,JP,JY,JX));    StalkRsrvElms_pft=0._r8
  allocate(HuskStrutElms_pft(NumPlantChemElms,JP,JY,JX));    HuskStrutElms_pft=0._r8
  allocate(EarStrutElms_pft(NumPlantChemElms,JP,JY,JX));    EarStrutElms_pft=0._r8
  allocate(GrainStrutElms_pft(NumPlantChemElms,JP,JY,JX));     GrainStrutElms_pft=0._r8
  allocate(CanopyLeafShethC_pft(JP,JY,JX));     CanopyLeafShethC_pft=0._r8
  allocate(CanopyLeafAreaZ_pft(NumOfCanopyLayers,JP,JY,JX)); CanopyLeafAreaZ_pft=0._r8
  allocate(CO2NetFix_pft(JP,JY,JX));     CO2NetFix_pft=0._r8
  allocate(CanopyLeafCLyr_pft(NumOfCanopyLayers,JP,JY,JX)); CanopyLeafCLyr_pft=0._r8
  allocate(CanopyNonstElms_pft(NumPlantChemElms,JP,JY,JX));   CanopyNonstElms_pft=0._r8
  allocate(CanopyNonstElmConc_pft(NumPlantChemElms,JP,JY,JX));   CanopyNonstElmConc_pft=0._r8
  allocate(CanopyStemAreaZ_pft(NumOfCanopyLayers,JP,JY,JX)); CanopyStemAreaZ_pft=0._r8
  allocate(CanopyNodulElms_pft(NumPlantChemElms,JP,JY,JX));CanopyNodulElms_pft=0._r8
  allocate(CanopyNodulNonstElms_pft(NumPlantChemElms,JP,JY,JX));   CanopyNodulNonstElms_pft=0._r8
  allocate(StalkBiomassC_brch(MaxNumBranches,JP,JY,JX));StalkBiomassC_brch=0._r8
  allocate(ShootElms_pft(NumPlantChemElms,JP,JY,JX));ShootElms_pft=0._r8
  allocate(ShootC4NonstC_brch(MaxNumBranches,JP,JY,JX));ShootC4NonstC_brch=0._r8
  allocate(CanopyNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); CanopyNonstElms_brch=0._r8
  allocate(LeafPetolBiomassC_brch(MaxNumBranches,JP,JY,JX)); LeafPetolBiomassC_brch=0._r8
  allocate(ShootStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));ShootStrutElms_brch=0._r8
  allocate(LeafStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); LeafStrutElms_brch=0._r8
  allocate(PetoleStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));PetoleStrutElms_brch=0._r8
  allocate(StalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));StalkStrutElms_brch=0._r8
  allocate(StalkRsrvElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));StalkRsrvElms_brch=0._r8
  allocate(HuskStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));HuskStrutElms_brch=0._r8
  allocate(EarStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));EarStrutElms_brch=0._r8
  allocate(GrainStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); GrainStrutElms_brch=0._r8
  allocate(LeafPetoNonstElmConc_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));LeafPetoNonstElmConc_brch=0._r8
  allocate(CanopyNodulNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));CanopyNodulNonstElms_brch=0._r8
  allocate(CanopyNodulStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); CanopyNodulStrutElms_brch=0._r8
  allocate(PetioleChemElmRemob_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));PetioleChemElmRemob_brch=0._r8
  allocate(SenecStalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));SenecStalkStrutElms_brch=0._r8
  allocate(LeafChemElmRemob_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); LeafChemElmRemob_brch=0._r8
  allocate(LeafElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafElmntNode_brch=0._r8
  allocate(PetioleElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleElmntNode_brch=0._r8
  allocate(InternodeStrutElms_brch(NumPlantChemElms,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));InternodeStrutElms_brch=0._r8
  allocate(LeafElmsByLayerNode_brch(NumPlantChemElms,NumOfCanopyLayers,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));
  LeafElmsByLayerNode_brch=0._r8
  allocate(CanopyLeafArea_lpft(NumOfCanopyLayers,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CanopyLeafArea_lpft=0._r8
  allocate(LeafProteinCNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafProteinCNode_brch=0._r8
  allocate(PetoleProteinCNode_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetoleProteinCNode_brch=0._r8
  allocate(NoduleNonstructCconc_pft(JP,JY,JX));   NoduleNonstructCconc_pft=0._r8
  allocate(GrainSeedBiomCMean_brch(MaxNumBranches,JP,JY,JX)); GrainSeedBiomCMean_brch=0._r8
  allocate(StandDeadKCompElms_pft(NumPlantChemElms,jsken,JP,JY,JX)); StandDeadKCompElms_pft=0._r8
  allocate(StandDeadStrutElms_pft(NumPlantChemElms,JP,JY,JX));    StandDeadStrutElms_pft=0._r8
  allocate(SeasonalNonstElms_pft(NumPlantChemElms,JP,JY,JX));  SeasonalNonstElms_pft=0._r8
  allocate(SeedCPlanted_pft(JP,JY,JX));    SeedCPlanted_pft=0._r8
  allocate(AvgCanopyBiomC2Graze_pft(JP,JY,JX));   AvgCanopyBiomC2Graze_pft=0._r8
  end subroutine InitCanopyData

!----------------------------------------------------------------------
  subroutine DestructCanopyData
  use abortutils, only : destroy
  implicit none
  call destroy(canopy_growth_pft)
  call destroy(CO2FixCL_pft)
  call destroy(CO2FixLL_pft)
  call destroy(StomatalStress_pft)
  call destroy(CanopyPARalbedo_pft)
  call destroy(RadPARLeafTransmis_pft)
  call destroy(LeafSWabsorpty_pft)
  call destroy(LeafPARabsorpty_pft)
  call destroy(CuticleResist_pft)
  call destroy(CO2CuticleResist_pft)
  call destroy(H2OCuticleResist_pft)
  call destroy(RCS)
  call destroy(CanPStomaResistH2O_pft)
  call destroy(MinCanPStomaResistH2O_pft)
  call destroy(CanopyBndlResist_col)
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
  call destroy(ChillHours_pft)
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
  call destroy(CPOOL4_node)
  call destroy(CMassHCO3BundleSheath_node)
  call destroy(RubiscoActivity_brch)
  call destroy(NutrientCtrlonC4Carboxy_node)
  call destroy(C4PhotosynDowreg_brch)
  call destroy(NetCO2Flx2Canopy_col)
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
  call destroy(RadNet2Canopy_pft)
  call destroy(LWRadCanopy_pft)
  call destroy(RadSWbyCanopy_pft)
  call destroy(RadPARbyCanopy_pft)
  call destroy(FracPARads2Canopy_pft)
  call destroy(TAU_RadThru)
  call destroy(TAU_DirRadTransm)
  call destroy(FracSWRad2Grnd_col)
  call destroy(RadSWGrnd_col)
  call destroy(LWRadCanGPrev_col)
  call destroy(LWRadGrnd)
  call destroy(CanH2OHeldVg_col)
  call destroy(Prec2Canopy_col)
  call destroy(PrecIntceptByCanopy_col)
  call destroy(EvapTransHeat_pft)
  call destroy(HeatXAir2PCan_pft)
  call destroy(HeatStorCanopy_pft)
  call destroy(ENGYX_pft)
  call destroy(VHeatCapCanP_pft)
  call destroy(PSICanopy_pft)
  call destroy(PSICanopyTurg_pft)
  call destroy(PSICanopyOsmo_pft)
  call destroy(CanopyBndlResist_pft)
  call destroy(Transpiration_pft)
  call destroy(VapXAir2Canopy_pft)
  call destroy(CanopyWater_pft)
  call destroy(QvET_col)
  call destroy(VapXAir2Canopy_col)
  call destroy(CanopyHeatStor_col)
  call destroy(HeatFlx2Canopy_col)
  call destroy(CanWat_col)
  call destroy(LWRadCanG)
  call destroy(RadSWLeafAlbedo_pft)
  call destroy(RadSWLeafTransmis_pft)
  call destroy(PrecIntcptByCanopy_pft)
  call destroy(WatByPCanopy_pft)
  call destroy(TKC)
  call destroy(TCelciusCanopy_pft)
  call destroy(DeltaTKC_pft)
  call destroy(TKCanopy_pft)
  call destroy(CPOOL3_node)
  call destroy(ShootElms_pft)
  call destroy(NetCumElmntFlx2Plant_pft)
  call destroy(tCanLeafC_cl)
  call destroy(ElmAllocmat4Litr)
  call destroy(ShootStrutElms_pft)
  call destroy(LeafStrutElms_pft)
  call destroy(PetoleStrutElms_pft)
  call destroy(StalkStrutElms_pft)
  call destroy(CanopyStalkC_pft)
  call destroy(StalkRsrvElms_pft)
  call destroy(HuskStrutElms_pft)
  call destroy(EarStrutElms_pft)
  call destroy(GrainStrutElms_pft)
  call destroy(CanopyLeafShethC_pft)
  call destroy(CanopyLeafAreaZ_pft)
  call destroy(CO2NetFix_pft)
  call destroy(CanopyLeafCLyr_pft)
  call destroy(CanopyNonstElms_pft)
  call destroy(CanopyNonstElmConc_pft)
  call destroy(CanopyStemAreaZ_pft)
  call destroy(CanopyNodulNonstElms_pft)
  call destroy(CanopyNodulElms_pft)
  call destroy(StalkBiomassC_brch)
  call destroy(CanopyNonstElms_brch)
  call destroy(ShootC4NonstC_brch)
  call destroy(LeafPetolBiomassC_brch)
  call destroy(ShootStrutElms_brch)
  call destroy(LeafStrutElms_brch)
  call destroy(PetoleStrutElms_brch)
  call destroy(StalkStrutElms_brch)
  call destroy(StalkRsrvElms_brch)
  call destroy(HuskStrutElms_brch)
  call destroy(EarStrutElms_brch)
  call destroy(GrainStrutElms_brch)
  call destroy(LeafPetoNonstElmConc_brch)
  call destroy(CanopyNodulNonstElms_brch)
  call destroy(CanopyNodulStrutElms_brch)
  call destroy(PetioleChemElmRemob_brch)
  call destroy(SenecStalkStrutElms_brch)
  call destroy(LeafChemElmRemob_brch)
  call destroy(LeafElmntNode_brch)
  call destroy(PetioleElmntNode_brch)
  call destroy(InternodeStrutElms_brch)
  call destroy(LeafElmsByLayerNode_brch)
  call destroy(CanopyLeafArea_lpft)
  call destroy(LeafProteinCNode_brch)
  call destroy(PetoleProteinCNode_brch)
  call destroy(NoduleNonstructCconc_pft)
  call destroy(GrainSeedBiomCMean_brch)
  call destroy(StandDeadKCompElms_pft)
  call destroy(StandDeadStrutElms_pft)
  call destroy(SeasonalNonstElms_pft)
  call destroy(SeedCPlanted_pft)
  call destroy(AvgCanopyBiomC2Graze_pft)
  end subroutine DestructCanopyData

end module CanopyDataType
