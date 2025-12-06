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

  real(r8),target,allocatable ::  canopy_growth_pft(:,:,:)                   !canopy structural growth rate, [gC/h]
  real(r8),target,allocatable ::  StomatalStress_pft(:,:,:)                  !stomatal stress from water/turgor,(0,1), [-]
  real(r8),target,allocatable ::  CanopyPARalbedo_pft(:,:,:)                 !canopy PAR albedo , [-]
  real(r8),target,allocatable ::  RadPARLeafTransmis_pft(:,:,:)              !canopy PAR transmissivity , [-]
  real(r8),target,allocatable ::  LeafSWabsorpty_pft(:,:,:)                  !canopy shortwave absorptivity , [-]
  real(r8),target,allocatable ::  LeafPARabsorpty_pft(:,:,:)                 !canopy PAR absorptivity , [-]
  real(r8),target,allocatable ::  CuticleResist_pft(:,:,:)                   !maximum stomatal resistance to vapor, [s m-1]
  real(r8),target,allocatable ::  CO2CuticleResist_pft(:,:,:)                !maximum stomatal resistance to CO2, [s h-1]
  real(r8),target,allocatable ::  H2OCuticleResist_pft(:,:,:)                !maximum stomatal resistance to vapor, [s h-1]
  real(r8),target,allocatable ::  RCS_pft(:,:,:)                             !e-folding turgor pressure for stomatal resistance, [MPa]
  real(r8),target,allocatable ::  CanPStomaResistH2O_pft(:,:,:)              !canopy stomatal resistance, [h m-1]
  real(r8),target,allocatable ::  CanopyMinStomaResistH2O_pft(:,:,:)           !canopy minimum stomatal resistance, [s m-1]
  real(r8),target,allocatable ::  CanopyBndlResist_col(:,:)                  !canopy boundary layer resistance, [m h-1]
  real(r8),target,allocatable ::  O2I_pft(:,:,:)                             !leaf gaseous O2 concentration, [umol m-3]
  real(r8),target,allocatable ::  LeafIntracellularCO2_pft(:,:,:)            !leaf gaseous CO2 concentration, [umol m-3]
  real(r8),target,allocatable ::  AirConc_pft(:,:,:)                         !total gas concentration, [mol m-3]
  real(r8),target,allocatable ::  DiffCO2Atmos2Intracel_pft(:,:,:)           !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8),target,allocatable ::  CanopyGasCO2_pft(:,:,:)                    !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8),target,allocatable ::  aquCO2Intraleaf_pft(:,:,:)                 !leaf aqueous CO2 concentration, [uM]
  real(r8),target,allocatable ::  O2L_pft(:,:,:)                                 !leaf aqueous O2 concentration, [uM]
  real(r8),target,allocatable ::  CO2Solubility_pft(:,:,:)                   !leaf CO2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  LeafO2Solubility_pft(:,:,:)                !leaf O2 solubility, [uM /umol mol-1]
  real(r8),target,allocatable ::  Km4LeafaqCO2_pft(:,:,:)                    !leaf aqueous CO2 Km no O2, [uM]
  real(r8),target,allocatable ::  Km4RubiscoCarboxy_pft(:,:,:)               !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8),target,allocatable ::  ChillHours_pft(:,:,:)                      !chilling effect on CO2 fixation, [-]
  real(r8),target,allocatable ::  Vmax4RubiscoCarboxy_node(:,:,:,:,:)         !maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  CO2lmtRubiscoCarboxyRate_node(:,:,:,:,:)   !carboxylation rate, [umol m-2 s-1]
  real(r8),target,allocatable ::  CO2CompenPoint_node(:,:,:,:,:)             !CO2 compensation point, [uM]
  real(r8),target,allocatable ::  LigthSatCarboxyRate_node(:,:,:,:,:)        !maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  RubiscoCarboxyEff_node(:,:,:,:,:)          !carboxylation efficiency, [umol umol-1]
  real(r8),target,allocatable ::  ProteinCperm2LeafArea_node(:,:,:,:,:)      !Protein C per m2 of leaf aera,            [gC (leaf area m-2)]
  real(r8),target,allocatable ::  CMassCO2BundleSheath_node(:,:,:,:,:)       !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  Vmax4PEPCarboxy_node(:,:,:,:,:)             !maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  CO2lmtPEPCarboxyRate_node(:,:,:,:,:)       !C4 carboxylation rate, [umol m-2 s-1]
  real(r8),target,allocatable ::  LigthSatC4CarboxyRate_node(:,:,:,:,:)      !maximum  light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8),target,allocatable ::  C4CarboxyEff_node(:,:,:,:,:)               !C4 carboxylation efficiency, [umol umol-1]
  real(r8),target,allocatable ::  CPOOL4_node(:,:,:,:,:)                     !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  CMassHCO3BundleSheath_node(:,:,:,:,:)      !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8),target,allocatable ::  LeafProteinC_brch(:,:,:,:)                  !Protein C for the branches, [gC protein d-2]
  real(r8),target,allocatable ::  LeafProteinCperm2LA_pft(:,:,:)             !Protein C for the plant, [gC protein m-2 leaf area]
  real(r8),target,allocatable ::  RubiscoActivity_brch(:,:,:,:)              !branch down-regulation of CO2 fixation, [-]
  real(r8),target,allocatable ::  NutrientCtrlonC4Carboxy_node(:,:,:,:,:)    !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  C4PhotosynDowreg_brch(:,:,:,:)             !down-regulation of C4 photosynthesis, [-]
  real(r8),target,allocatable ::  NetCO2Flx2Canopy_col(:,:)                  !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  VmaxSpecRubCarboxyRef_pft(:,:,:)               !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VmaxRubOxyRef_pft(:,:,:)                   !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  VmaxPEPCarboxyRef_pft(:,:,:)               !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  XKCO2_pft(:,:,:)                           !Km for rubisco carboxylase activity, [uM]
  real(r8),target,allocatable ::  XKO2_pft(:,:,:)                            !Km for rubisco oxygenase activity, [uM]
  real(r8),target,allocatable ::  Km4PEPCarboxy_pft(:,:,:)                   !Km for PEP carboxylase activity, [uM]
  real(r8),target,allocatable ::  LeafRubisco2Protein_pft(:,:,:)                    !leaf rubisco content, [g g-1]
  real(r8),target,allocatable ::  LeafPEP2Protein_pft(:,:,:)     !leaf PEP carboxylase content, [g g-1]
  real(r8),target,allocatable ::  SpecLeafChlAct_pft(:,:,:)                !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8),target,allocatable ::  LeafC3Chl2Protein_pft(:,:,:)             !leaf C3 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  LeafC4Chl2Protein_pft(:,:,:)             !leaf C4 chlorophyll content, [g g-1]
  real(r8),target,allocatable ::  CanopyCi2CaRatio_pft(:,:,:)                      !Ci:Ca ratio, [-]
  real(r8),target,allocatable ::  RadNet2Canopy_pft(:,:,:)                   !canopy net radiation , [MJ d-2 h-1] >0
  real(r8),target,allocatable ::  LWRadCanopy_pft(:,:,:)                     !canopy longwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadSWbyCanopy_pft(:,:,:)                   !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadPARbyCanopy_pft(:,:,:)                  !canopy absorbed PAR , [umol m-2 s-1]
  real(r8),target,allocatable ::  FracPARads2Canopy_pft(:,:,:)               !fraction of incoming PAR absorbed by canopy, [-]
  real(r8),target,allocatable ::  TAU_RadThru(:,:,:)                         !fraction of radiation transmitted by canopy layer, [-]
  real(r8),target,allocatable ::  TAU_DirectRTransmit(:,:,:)                 !fraction of radiation intercepted by canopy layer, [-]
  real(r8),target,allocatable ::  FracSWRad2Grnd_col(:,:)                    !fraction of radiation intercepted by ground surface, [-]
  real(r8),target,allocatable ::  RadSWGrnd_col(:,:)                         !shortwave radiation incident on ground surface, [MJ h-1]
  real(r8),target,allocatable ::  LWRadCanGPrev_col(:,:)                     !longwave radiation emitted by canopy, [MJ h-1]
  real(r8),target,allocatable ::  LWRadGrnd_col(:,:)                         !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8),target,allocatable ::  WatHeldOnCanopy_col(:,:)                   !canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  Prec2Canopy_col(:,:)                       !precipitation to canopy over the grid, [MJ d-2 h-1]
  real(r8),target,allocatable ::  PrecIntceptByCanopy_col(:,:)               !grid net precipitation water interception to canopy, [MJ d-2 t-1]
  real(r8),target,allocatable ::  EvapTransLHeat_pft(:,:,:)                  !canopy latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatXAir2PCan_pft(:,:,:)                   !air to canopy sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatStorCanopy_pft(:,:,:)                  !canopy storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  ENGYX_pft(:,:,:)                           !canopy heat storage from previous time step, [MJ d-2]
  real(r8),target,allocatable ::  VHeatCapCanopy_pft(:,:,:)                  !canopy heat capacity, [MJ d-2 K-1]
  real(r8),target,allocatable ::  PSICanopy_pft(:,:,:)                       !plant canopy total water potential , [Mpa]
  real(r8),target,allocatable ::  PSICanopyTurg_pft(:,:,:)                   !plant canopy turgor water potential, [Mpa]
  real(r8),target,allocatable ::  PSICanopyOsmo_pft(:,:,:)                   !platn canopy osmotic water potential, [Mpa]
  real(r8),target,allocatable ::  CanopyBndlResist_pft(:,:,:)                !canopy boundary layer resistance, [h m-1]
  real(r8),target,allocatable ::  Transpiration_pft(:,:,:)                   !canopy transpiration, [m3 d-2 h-1]
  real(r8),target,allocatable ::  VapXAir2Canopy_pft(:,:,:)                  !negative of canopy evaporation, [m2 d-2 h-1]
  real(r8),target,allocatable ::  CanopyBiomWater_pft(:,:,:)                 !canopy water content associated with dry matter, [m3 d-2]
  real(r8),target,allocatable ::  CanopyWaterMassBeg_col(:,:)                !Canopy water before mass balance check [m3 d-2]
  real(r8),target,allocatable ::  CanopyWaterMassEnd_col(:,:)                !Canopy water at mass balance check [m3 d-2]
  real(r8),target,allocatable ::  HeatCanopy2Dist_col(:,:)                   !Canopy heat content loss to disturbance, [MJ d-2]
  real(r8),target,allocatable ::  QCanopyWat2Dist_col(:,:)                   !canopy water loss to disturbance, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QVegET_col(:,:)                            !total canopy evaporation + transpiration, [m3 d-2 h-1]
  real(r8),target,allocatable ::  VapXAir2Canopy_col(:,:)                    !total canopy evaporation, [m3 d-2]
  real(r8),target,allocatable ::  CanopyHeatStor_col(:,:)                    !total canopy heat content, [MJ  d-2]
  real(r8),target,allocatable ::  HeatFlx2Canopy_col(:,:)                    !total canopy heat flux, [MJ  d-2 h-1]
  real(r8),target,allocatable ::  CanopyWat_col(:,:)                         !total canopy water content stored in dry matter, [m3 d-2]
  real(r8),target,allocatable ::  LWRadCanG_col(:,:)                         !total canopy LW emission, [MJ d-2 h-1]
  real(r8),target,allocatable ::  RadSWLeafAlbedo_pft(:,:,:)                 !canopy shortwave albedo , [-]
  real(r8),target,allocatable ::  RadSWLeafTransmis_pft(:,:,:)               !canopy shortwave transmissivity , [-]
  real(r8),target,allocatable ::  PrecIntcptByCanopy_pft(:,:,:)              !water flux into plant canopy, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WatHeldOnCanopy_pft(:,:,:)                 !canopy held water content, [m3 d-2]
  real(r8),target,allocatable ::  TKC_pft(:,:,:)                             !canopy temperature after energy iteration, [K]
  real(r8),target,allocatable ::  TdegCCanopy_pft(:,:,:)                     !canopy temperature, [oC]
  real(r8),target,allocatable ::  DeltaTKC_pft(:,:,:)                        !change in canopy temperature, [K]
  real(r8),target,allocatable ::  TKCanopy_pft(:,:,:)                        !canopy temperature during canopy energy iteration, [K]
  real(r8),target,allocatable ::  CPOOL3_node(:,:,:,:,:)                     !bundle sheath C4 carbon product to support C3 photosynthesis, [gC d-2]
  real(r8),target,allocatable ::  NetCumElmntFlx2Plant_pft(:,:,:,:)          !effect of canopy chemical element status on seed setting, []
  real(r8),target,allocatable ::  tCanLeafC_clyr(:,:,:)                      !total leaf mass, [g d-2]
  real(r8),target,allocatable ::  PlantElmAllocMat4Litr(:,:,:,:,:,:)              !litter kinetic fraction, [-]
  real(r8),target,allocatable ::  C4PhotoShootNonstC_brch(:,:,:,:)           !C4 specific nonstructural shoot C in branch, [gC d-2]
  real(r8),target,allocatable ::  ShootElms_pft(:,:,:,:)                     !canopy shoot chemical element, [g d-2]
  real(r8),target,allocatable ::  LeafStrutElms_pft(:,:,:,:)                 !canopy leaf chemical element, [g d-2]
  real(r8),target,allocatable ::  PetoleStrutElms_pft(:,:,:,:)               !canopy sheath chemical element , [g d-2]
  real(r8),target,allocatable ::  StalkStrutElms_pft(:,:,:,:)                !canopy stalk chemical element, [g d-2]
  real(r8),target,allocatable ::  CanopySapwoodC_pft(:,:,:)                    !canopy active stalk C, [g d-2]
  real(r8),target,allocatable ::  StalkRsrvElms_pft(:,:,:,:)                 !canopy reserve chemical element, [g d-2]
  real(r8),target,allocatable ::  HuskStrutElms_pft(:,:,:,:)                 !canopy husk chemical element, [g d-2]
  real(r8),target,allocatable ::  EarStrutElms_pft(:,:,:,:)                  !canopy ear chemical element, [g d-2]
  real(r8),target,allocatable ::  GrainStrutElms_pft(:,:,:,:)                !canopy grain chemical element, [g d-2]
  real(r8),target,allocatable ::  CanopyLeafSheathC_pft(:,:,:)                !plant canopy leaf + sheath C, [gC d-2]
  real(r8),target,allocatable ::  CanopyLeafAreaZ_pft(:,:,:,:)               !canopy layer-distributed leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CO2NetFix_pft(:,:,:)                       !canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RCanMaintDef_CO2_pft(:,:,:)                !canopy maintenance respiraiton deficit as CO2,  [gC d-2 h-1]  
  real(r8),target,allocatable ::  CanopyLeafCLyr_pft(:,:,:,:)                !canopy layer leaf C, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstElms_pft(:,:,:,:)               !canopy nonstructural chemical element, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstElmConc_pft(:,:,:,:)            !canopy nonstructural chemical element concentration, [g d-2]
  real(r8),target,allocatable ::  CanopyStemAreaZ_pft(:,:,:,:)               !plant canopy layer stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyNodulNonstElms_pft(:,:,:,:)          !canopy nodule nonstructural chemical element, [g d-2]
  real(r8),target,allocatable ::  ShootNoduleElms_pft(:,:,:,:)               !canopy nodule chemical elemental biomass [g d-2]
  real(r8),target,allocatable ::  SapwoodBiomassC_brch(:,:,:,:)              !branch active stalk C, [g d-2]
  real(r8),target,allocatable ::  CanopyNonstElms_brch(:,:,:,:,:)            !branch nonstructural chemical element, [g d-2]
  real(r8),target,allocatable ::  CanopyLeafSheathC_brch(:,:,:,:)            !plant branch leaf + sheath C, [g d-2]
  real(r8),target,allocatable ::  ShootElms_brch(:,:,:,:,:)             !branch shoot C, [g d-2]
  real(r8),target,allocatable ::  LeafStrutElms_brch(:,:,:,:,:)              !branch leaf chemical element, [g d-2]
  real(r8),target,allocatable ::  PetoleStrutElms_brch(:,:,:,:,:)            !branch sheath chemical element , [g d-2]
  real(r8),target,allocatable ::  StalkStrutElms_brch(:,:,:,:,:)             !branch stalk chemical element, [g d-2]
  real(r8),target,allocatable ::  StalkRsrvElms_brch(:,:,:,:,:)              !branch reserve chemical element, [g d-2]
  real(r8),target,allocatable ::  HuskStrutElms_brch(:,:,:,:,:)              !branch husk chemical element, [g d-2]
  real(r8),target,allocatable ::  EarStrutElms_brch(:,:,:,:,:)               !branch ear chemical element, [g d-2]
  real(r8),target,allocatable ::  GrainStrutElms_brch(:,:,:,:,:)             !branch grain chemical element, [g d-2]
  real(r8),target,allocatable ::  LeafPetoNonstElmConc_brch(:,:,:,:,:)       !branch nonstructural C concentration, [g d-2]
  real(r8),target,allocatable ::  CanopyNodulNonstElms_brch(:,:,:,:,:)       !branch nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  CanopyNodulStrutElms_brch(:,:,:,:,:)       !branch nodule chemical element, [g d-2]
  real(r8),target,allocatable ::  SenecStalkStrutElms_brch(:,:,:,:,:)        !senescing branch stalk structural chemical elements, [g d-2]
  real(r8),target,allocatable ::  LeafElmntNode_brch(:,:,:,:,:,:)            !leaf chemical element, [g d-2]
  real(r8),target,allocatable ::  PetioleElmntNode_brch(:,:,:,:,:,:)         !sheath chemical element , [g d-2]
  real(r8),target,allocatable ::  StructInternodeElms_brch(:,:,:,:,:,:)       !internode chemical element, [g d-2]
  real(r8),target,allocatable ::  LeafLayerElms_node(:,:,:,:,:,:,:)    !layer leaf chemical element, [g d-2]
  real(r8),target,allocatable ::  CanopyLeafArea_lnode(:,:,:,:,:,:)           !layer leaf area, [m2 d-2]
  real(r8),target,allocatable ::  LeafProteinC_node(:,:,:,:,:)           !layer leaf protein C, [g d-2]
  real(r8),target,allocatable ::  PetoleProteinC_node(:,:,:,:,:)         !layer sheath protein C, [g d-2]
  real(r8),target,allocatable ::  CanopyNoduleNonstCConc_pft(:,:,:)            !nodule nonstructural C, [g d-2]
  real(r8),target,allocatable ::  GrainSeedBiomCMean_brch(:,:,:,:)           !maximum grain C during grain fill, [g d-2]
  real(r8),target,allocatable ::  CanopyNLimFactor_brch(:,:,:,:)             !Canopy N-limitation factor, [0->1] weaker limitation,[-]
  real(r8),target,allocatable ::  CanopyPLimFactor_brch(:,:,:,:)             !Canopy P-limitation factor, [0->1] weaker limitation,[-]
  real(r8),target,allocatable ::  StandDeadKCompElms_pft(:,:,:,:,:)          !standing dead chemical element fraction, [g d-2]
  real(r8),target,allocatable ::  StandDeadStrutElms_pft(:,:,:,:)            !standing dead chemical element, [g d-2]
  real(r8),target,allocatable ::  SeasonalNonstElms_pft(:,:,:,:)             !plant stored nonstructural chemical element, [g d-2]
  real(r8),target,allocatable ::  SeedPlantedElm_pft(:,:,:,:)              !plant stored nonstructural chemical elements at planting, [g d-2]
  REAL(R8),target,allocatable ::  AvgCanopyBiomC2Graze_pft(:,:,:)            !landscape average canopy shoot C, [g d-2]
  real(r8),target,allocatable :: CO2FixCL_pft(:,:,:)                         !CO2-limited carboxylation rate, [gC d2 h-1] 
  real(r8),target,allocatable :: CO2FixLL_pft(:,:,:)                         !Light-limited carboxylation rate,[gC d2 h-1]
  real(r8),target,allocatable :: CanopyMassC_pft(:,:,:)                      !Canopy biomass, [gC d-2]
  real(r8),target,allocatable :: fNCLFW_pft(:,:,:)                            !NC ratio of growing leaf, [gN/gC]
  real(r8),target,allocatable :: fPCLFW_pft(:,:,:)                           !PC ratio of growing leaf, [gP/gC]
  real(r8),target,allocatable :: CanopyVcMaxRubisco25C_pft(:,:,:)               !Canopy VcMax for rubisco carboxylation, [umol h-1 m-2]
  real(r8),target,allocatable :: CanopyVoMaxRubisco25C_pft(:,:,:)               !Canopy VoMax for rubisco oxygenation, [umol h-1 m-2]
  real(r8),target,allocatable :: CanopyVcMaxPEP25C_pft(:,:,:)                   !Canopy VcMax in PEP C4 fixation, [umol h-1 m-2]
  real(r8),target,allocatable :: ElectronTransptJmax25C_pft(:,:,:)           !Canopy Jmax at reference temperature, [umol e- s-1 m-2]
  real(r8),target,allocatable :: TFN_Carboxy_pft(:,:,:)                      ! temperature dependence of carboxylation, [-]
  real(r8),target,allocatable :: TFN_Oxygen_pft(:,:,:)                       !temperature dependence of oxygenation, [-] 
  real(r8),target,allocatable :: TFN_eTranspt_pft(:,:,:)                     !temperature dependence of electron transport, [-]
  real(r8),target,allocatable :: LeafAreaSunlit_pft(:,:,:)                   !leaf irradiated surface area, [m2 d-2]    
  real(r8),target,allocatable :: PARSunlit_pft(:,:,:)                        !PAR absorbed by sunlit leaf, [umol m-2 s-1]
  real(r8),target,allocatable :: PARSunsha_pft(:,:,:)                        !PAR absorbed by sun-shaded leaf, [umol m-2 s-1]
  real(r8),target,allocatable :: CH2OSunlit_pft(:,:,:)                       !carbon fixation by sun-lit leaf, [gC d-2 h-1]
  real(r8),target,allocatable :: CH2OSunsha_pft(:,:,:)                       !carbon fixation by sun-shaded leaf, [gC d-2 h-1]    
!  real(r8),target,allocatable :: LeafC3ChlC_brch(:,:,:,:)                    !Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]     
!  real(r8),target,allocatable :: LeafC4ChlC_brch(:,:,:,:)                    !Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]     
!  real(r8),target,allocatable :: LeafRubiscoC_brch(:,:,:,:)                  !Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
!  real(r8),target,allocatable :: LeafPEPC_brch(:,:,:,:)                      !PEP C for the branches, [gC PEP d-2]
  real(r8),target,allocatable :: LeafC3ChlCperm2LA_pft(:,:,:)                       !Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]     
  real(r8),target,allocatable :: LeafC4ChlCperm2LA_pft(:,:,:)                       !Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]     
  real(r8),target,allocatable :: LeafRubiscoCperm2LA_pft(:,:,:)                     !Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
  real(r8),target,allocatable :: LeafPEPCperm2LA_pft(:,:,:)                         !PEP C for the branches, [gC PEP d-2]
  real(r8),target,allocatable :: SpecificLeafArea_pft(:,:,:)                   !specifc leaf area per g C of leaf mass, [m2 leaf area (gC leaf C)-1]

  contains
!----------------------------------------------------------------------

  subroutine InitCanopyData

  implicit none

  allocate(PARSunlit_pft(JP,JY,JX));PARSunlit_pft=0._r8
  allocate(PARSunsha_pft(JP,JY,JX));PARSunsha_pft=0._r8
  allocate(CH2OSunlit_pft(JP,JY,JX));CH2OSunlit_pft=0._r8
  allocate(CH2OSunsha_pft(JP,JY,JX));CH2OSunsha_pft=0._r8
  ALLOCATE(SpecificLeafArea_pft(JP,JY,JX)); SpecificLeafArea_pft=0._r8
  allocate(LeafAreaSunlit_pft(JP,JY,JX)); LeafAreaSunlit_pft=0._r8
  allocate(TFN_Carboxy_pft(JP,JY,JX));TFN_Carboxy_pft=0._r8
  allocate(TFN_Oxygen_pft(JP,JY,JX)); TFN_Oxygen_pft=0._r8
  allocate(TFN_eTranspt_pft(JP,JY,JX)); TFN_eTranspt_pft=0._r8
  allocate(CanopyVcMaxRubisco25C_pft(JP,JY,JX));CanopyVcMaxRubisco25C_pft=0._r8
  allocate(CanopyVoMaxRubisco25C_pft(JP,JY,JX));CanopyVoMaxRubisco25C_pft=0._r8  
  allocate(CanopyVcMaxPEP25C_pft(JP,JY,JX)); CanopyVcMaxPEP25C_pft=0._r8
  allocate(ElectronTransptJmax25C_pft(JP,JY,JX));ElectronTransptJmax25C_pft=0._r8
  allocate(CanopyNLimFactor_brch(MaxNumBranches,JP,JY,JX));CanopyNLimFactor_brch=1._r8
  allocate(CanopyPLimFactor_brch(MaxNumBranches,JP,JY,JX));CanopyPLimFactor_brch=1._r8
  allocate(CanopyMassC_pft(JP,JY,jX)) ; CanopyMassC_pft=0._r8
  allocate(fNCLFW_pft(JP,JY,JX)); fNCLFW_pft=0._r8
  allocate(fPCLFW_pft(JP,JY,JX)); fPCLFW_pft=0._r8
  allocate(CanopyWaterMassBeg_col(JY,JX)); CanopyWaterMassBeg_col=0._r8
  allocate(CanopyWaterMassEnd_col(JY,JX)); CanopyWaterMassEnd_col=0._r8
  allocate(HeatCanopy2Dist_col(JY,JX)); HeatCanopy2Dist_col=0._r8
  allocate(QCanopyWat2Dist_col(JY,JX)); QCanopyWat2Dist_col=0._r8
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
  allocate(RCS_pft(JP,JY,JX));      RCS_pft=0._r8
  allocate(CanPStomaResistH2O_pft(JP,JY,JX));       CanPStomaResistH2O_pft=0._r8
  allocate(CanopyMinStomaResistH2O_pft(JP,JY,JX));     CanopyMinStomaResistH2O_pft=0._r8
  allocate(CanopyBndlResist_col(JY,JX));         CanopyBndlResist_col=0._r8
  allocate(O2I_pft(JP,JY,JX));      O2I_pft=0._r8
  allocate(LeafIntracellularCO2_pft(JP,JY,JX));     LeafIntracellularCO2_pft=0._r8
  allocate(AirConc_pft(JP,JY,JX));     AirConc_pft=0._r8
  allocate(DiffCO2Atmos2Intracel_pft(JP,JY,JX));     DiffCO2Atmos2Intracel_pft=0._r8
  allocate(CanopyGasCO2_pft(JP,JY,JX));     CanopyGasCO2_pft=0._r8
  allocate(aquCO2Intraleaf_pft(JP,JY,JX));     aquCO2Intraleaf_pft=0._r8
  allocate(O2L_pft(JP,JY,JX));      O2L_pft=0._r8
  allocate(CO2Solubility_pft(JP,JY,JX));     CO2Solubility_pft=0._r8
  allocate(LeafO2Solubility_pft(JP,JY,JX));      LeafO2Solubility_pft=0._r8
  allocate(Km4LeafaqCO2_pft(JP,JY,JX));   Km4LeafaqCO2_pft=0._r8
  allocate(Km4RubiscoCarboxy_pft(JP,JY,JX));   Km4RubiscoCarboxy_pft=0._r8
  allocate(ChillHours_pft(JP,JY,JX));    ChillHours_pft=0._r8
  allocate(Vmax4RubiscoCarboxy_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));Vmax4RubiscoCarboxy_node=0._r8
  allocate(CO2lmtRubiscoCarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CO2lmtRubiscoCarboxyRate_node=0._r8
  allocate(CO2CompenPoint_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CO2CompenPoint_node=0._r8
  allocate(LigthSatCarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LigthSatCarboxyRate_node=0._r8
  allocate(RubiscoCarboxyEff_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));RubiscoCarboxyEff_node=0._r8
  allocate(ProteinCperm2LeafArea_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));ProteinCperm2LeafArea_node=0._r8
  allocate(CMassCO2BundleSheath_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CMassCO2BundleSheath_node=0._r8
  allocate(Vmax4PEPCarboxy_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));Vmax4PEPCarboxy_node=0._r8
  allocate(CO2lmtPEPCarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CO2lmtPEPCarboxyRate_node=0._r8
  allocate(LigthSatC4CarboxyRate_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LigthSatC4CarboxyRate_node=0._r8
  allocate(C4CarboxyEff_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));C4CarboxyEff_node=0._r8
  allocate(CPOOL4_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CPOOL4_node=0._r8
  allocate(CMassHCO3BundleSheath_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CMassHCO3BundleSheath_node=0._r8
  allocate(RubiscoActivity_brch(MaxNumBranches,JP,JY,JX));  RubiscoActivity_brch=0._r8
  allocate(LeafProteinC_brch(MaxNumBranches,JP,JY,JX)); LeafProteinC_brch=0._r8
  allocate(LeafProteinCperm2LA_pft(JP,JY,JX)); LeafProteinCperm2LA_pft=0._r8
  allocate(NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));NutrientCtrlonC4Carboxy_node=0._r8
  allocate(C4PhotosynDowreg_brch(MaxNumBranches,JP,JY,JX)); C4PhotosynDowreg_brch=0._r8
  allocate(NetCO2Flx2Canopy_col(JY,JX));       NetCO2Flx2Canopy_col=0._r8
  allocate(VmaxSpecRubCarboxyRef_pft(JP,JY,JX));     VmaxSpecRubCarboxyRef_pft=0._r8
  allocate(VmaxRubOxyRef_pft(JP,JY,JX));     VmaxRubOxyRef_pft=0._r8
  allocate(VmaxPEPCarboxyRef_pft(JP,JY,JX));    VmaxPEPCarboxyRef_pft=0._r8
  allocate(XKCO2_pft(JP,JY,JX));    XKCO2_pft=0._r8
  allocate(XKO2_pft(JP,JY,JX));     XKO2_pft=0._r8
  allocate(Km4PEPCarboxy_pft(JP,JY,JX));   Km4PEPCarboxy_pft=0._r8
  allocate(LeafRubisco2Protein_pft(JP,JY,JX));     LeafRubisco2Protein_pft=0._r8
  allocate(LeafPEP2Protein_pft(JP,JY,JX));     LeafPEP2Protein_pft=0._r8
  allocate(SpecLeafChlAct_pft(JP,JY,JX));     SpecLeafChlAct_pft=0._r8
  allocate(LeafC3Chl2Protein_pft(JP,JY,JX));      LeafC3Chl2Protein_pft=0._r8
  allocate(LeafC4Chl2Protein_pft(JP,JY,JX));     LeafC4Chl2Protein_pft=0._r8
  allocate(CanopyCi2CaRatio_pft(JP,JY,JX));     CanopyCi2CaRatio_pft=0._r8
  allocate(RadNet2Canopy_pft(JP,JY,JX));     RadNet2Canopy_pft=0._r8
  allocate(LWRadCanopy_pft(JP,JY,JX));    LWRadCanopy_pft=0._r8
  allocate(RadSWbyCanopy_pft(JP,JY,JX));     RadSWbyCanopy_pft=0._r8
  allocate(RadPARbyCanopy_pft(JP,JY,JX));     RadPARbyCanopy_pft=0._r8
  allocate(FracPARads2Canopy_pft(JP,JY,JX));    FracPARads2Canopy_pft=0._r8
  allocate(TAU_RadThru(NumCanopyLayers+1,JY,JX));   TAU_RadThru=0._r8
  allocate(TAU_DirectRTransmit(NumCanopyLayers+1,JY,JX));   TAU_DirectRTransmit=0._r8
  allocate(FracSWRad2Grnd_col(JY,JX));       FracSWRad2Grnd_col=0._r8
  allocate(RadSWGrnd_col(JY,JX));        RadSWGrnd_col=0._r8
  allocate(LWRadCanGPrev_col(JY,JX));      LWRadCanGPrev_col=0._r8
  allocate(LWRadGrnd_col(JY,JX));      LWRadGrnd_col=0._r8
  allocate(WatHeldOnCanopy_col(JY,JX));      WatHeldOnCanopy_col=0._r8
  allocate(Prec2Canopy_col(JY,JX));      Prec2Canopy_col=0._r8
  allocate(PrecIntceptByCanopy_col(JY,JX));       PrecIntceptByCanopy_col=0._r8
  allocate(EvapTransLHeat_pft(JP,JY,JX));    EvapTransLHeat_pft=0._r8
  allocate(HeatXAir2PCan_pft(JP,JY,JX));    HeatXAir2PCan_pft=0._r8
  allocate(HeatStorCanopy_pft(JP,JY,JX));    HeatStorCanopy_pft=0._r8
  allocate(ENGYX_pft(JP,JY,JX));    ENGYX_pft=0._r8
  allocate(VHeatCapCanopy_pft(JP,JY,JX));    VHeatCapCanopy_pft=0._r8
  allocate(PSICanopy_pft(JP,JY,JX));    PSICanopy_pft=0._r8
  allocate(PSICanopyTurg_pft(JP,JY,JX));    PSICanopyTurg_pft=0._r8
  allocate(PSICanopyOsmo_pft(JP,JY,JX));    PSICanopyOsmo_pft=0._r8
  allocate(CanopyBndlResist_pft(JP,JY,JX));       CanopyBndlResist_pft=0._r8
  allocate(Transpiration_pft(JP,JY,JX));       Transpiration_pft=0._r8
  allocate(VapXAir2Canopy_pft(JP,JY,JX));    VapXAir2Canopy_pft=0._r8
  allocate(CanopyBiomWater_pft(JP,JY,JX));    CanopyBiomWater_pft=0._r8
  allocate(QVegET_col(JY,JX));      QVegET_col=0._r8
  allocate(VapXAir2Canopy_col(JY,JX));      VapXAir2Canopy_col=0._r8
  allocate(CanopyHeatStor_col(JY,JX));      CanopyHeatStor_col=0._r8
  allocate(HeatFlx2Canopy_col(JY,JX));      HeatFlx2Canopy_col=0._r8
  allocate(CanopyWat_col(JY,JX));      CanopyWat_col=0._r8
  allocate(LWRadCanG_col(JY,JX));       LWRadCanG_col=0._r8
  allocate(RadSWLeafAlbedo_pft(JP,JY,JX));     RadSWLeafAlbedo_pft=0._r8
  allocate(RadSWLeafTransmis_pft(JP,JY,JX));     RadSWLeafTransmis_pft=0._r8
  allocate(PrecIntcptByCanopy_pft(JP,JY,JX));     PrecIntcptByCanopy_pft=0._r8
  allocate(WatHeldOnCanopy_pft(JP,JY,JX));    WatHeldOnCanopy_pft=0._r8
  allocate(TKC_pft(JP,JY,JX));      TKC_pft=0._r8
  allocate(TdegCCanopy_pft(JP,JY,JX));      TdegCCanopy_pft=0._r8
  allocate(DeltaTKC_pft(JP,JY,JX));     DeltaTKC_pft=0._r8
  allocate(TKCanopy_pft(JP,JY,JX));     TKCanopy_pft=0._r8
  allocate(CPOOL3_node(MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CPOOL3_node=0._r8
  allocate(NetCumElmntFlx2Plant_pft(NumPlantChemElms,JP,JY,JX));    NetCumElmntFlx2Plant_pft=0._r8
  allocate(tCanLeafC_clyr(NumCanopyLayers,JY,JX));    tCanLeafC_clyr=0._r8
  allocate(PlantElmAllocMat4Litr(NumPlantChemElms,0:NumLitterGroups,jsken,JP,JY,JX));PlantElmAllocMat4Litr=0._r8
  allocate(LeafStrutElms_pft(NumPlantChemElms,JP,JY,JX));  LeafStrutElms_pft=0._r8
  allocate(PetoleStrutElms_pft(NumPlantChemElms,JP,JY,JX)); PetoleStrutElms_pft=0._r8
  allocate(StalkStrutElms_pft(NumPlantChemElms,JP,JY,JX)); StalkStrutElms_pft=0._r8
  allocate(CanopySapwoodC_pft(JP,JY,JX));    CanopySapwoodC_pft=0._r8
  allocate(StalkRsrvElms_pft(NumPlantChemElms,JP,JY,JX));    StalkRsrvElms_pft=0._r8
  allocate(HuskStrutElms_pft(NumPlantChemElms,JP,JY,JX));    HuskStrutElms_pft=0._r8
  allocate(EarStrutElms_pft(NumPlantChemElms,JP,JY,JX));    EarStrutElms_pft=0._r8
  allocate(GrainStrutElms_pft(NumPlantChemElms,JP,JY,JX));     GrainStrutElms_pft=0._r8
  allocate(CanopyLeafSheathC_pft(JP,JY,JX));     CanopyLeafSheathC_pft=0._r8
  allocate(CanopyLeafAreaZ_pft(NumCanopyLayers,JP,JY,JX)); CanopyLeafAreaZ_pft=0._r8
  allocate(CO2NetFix_pft(JP,JY,JX));     CO2NetFix_pft=0._r8
  allocate(RCanMaintDef_CO2_pft(JP,JY,JX)); RCanMaintDef_CO2_pft=0._r8
  allocate(CanopyLeafCLyr_pft(NumCanopyLayers,JP,JY,JX)); CanopyLeafCLyr_pft=0._r8
  allocate(CanopyNonstElms_pft(NumPlantChemElms,JP,JY,JX));   CanopyNonstElms_pft=0._r8
  allocate(CanopyNonstElmConc_pft(NumPlantChemElms,JP,JY,JX));   CanopyNonstElmConc_pft=0._r8
  allocate(CanopyStemAreaZ_pft(NumCanopyLayers,JP,JY,JX)); CanopyStemAreaZ_pft=0._r8
  allocate(ShootNoduleElms_pft(NumPlantChemElms,JP,JY,JX));ShootNoduleElms_pft=0._r8
  allocate(CanopyNodulNonstElms_pft(NumPlantChemElms,JP,JY,JX));   CanopyNodulNonstElms_pft=0._r8
  allocate(SapwoodBiomassC_brch(MaxNumBranches,JP,JY,JX));SapwoodBiomassC_brch=0._r8
  allocate(ShootElms_pft(NumPlantChemElms,JP,JY,JX));ShootElms_pft=0._r8
!  allocate(LeafC3ChlC_brch(MaxNumBranches,JP,JY,JX)); LeafC3ChlC_brch=0._r8             
!  allocate(LeafC4ChlC_brch(MaxNumBranches,JP,JY,JX)); LeafC4ChlC_brch=0._r8
!  allocate(LeafRubiscoC_brch(MaxNumBranches,JP,JY,JX)); LeafRubiscoC_brch=0._r8           
!  allocate(LeafPEPC_brch(MaxNumBranches,JP,JY,JX));LeafPEPC_brch=0._r8               
  allocate(LeafC3ChlCperm2LA_pft(JP,JY,JX));LeafC3ChlCperm2LA_pft=0._r8                
  allocate(LeafC4ChlCperm2LA_pft(JP,JY,JX));LeafC4ChlCperm2LA_pft=0._r8
  allocate(LeafRubiscoCperm2LA_pft(JP,JY,JX));LeafRubiscoCperm2LA_pft=0._r8             
  allocate(LeafPEPCperm2LA_pft(JP,JY,JX));LeafPEPCperm2LA_pft=0._r8                  
  allocate(C4PhotoShootNonstC_brch(MaxNumBranches,JP,JY,JX));C4PhotoShootNonstC_brch=0._r8
  allocate(CanopyNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); CanopyNonstElms_brch=0._r8
  allocate(CanopyLeafSheathC_brch(MaxNumBranches,JP,JY,JX)); CanopyLeafSheathC_brch=0._r8
  allocate(ShootElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));ShootElms_brch=0._r8
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
  allocate(SenecStalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX));SenecStalkStrutElms_brch=0._r8
  allocate(LeafElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafElmntNode_brch=0._r8
  allocate(PetioleElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetioleElmntNode_brch=0._r8
  allocate(StructInternodeElms_brch(NumPlantChemElms,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));StructInternodeElms_brch=0._r8
  allocate(LeafLayerElms_node(NumPlantChemElms,NumCanopyLayers,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));
  LeafLayerElms_node=0._r8
  allocate(CanopyLeafArea_lnode(NumCanopyLayers,0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));CanopyLeafArea_lnode=0._r8
  allocate(LeafProteinC_node(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafProteinC_node=0._r8
  allocate(PetoleProteinC_node(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetoleProteinC_node=0._r8
  allocate(CanopyNoduleNonstCConc_pft(JP,JY,JX));   CanopyNoduleNonstCConc_pft=0._r8
  allocate(GrainSeedBiomCMean_brch(MaxNumBranches,JP,JY,JX)); GrainSeedBiomCMean_brch=0._r8
  allocate(StandDeadKCompElms_pft(NumPlantChemElms,jsken,JP,JY,JX)); StandDeadKCompElms_pft=0._r8
  allocate(StandDeadStrutElms_pft(NumPlantChemElms,JP,JY,JX));    StandDeadStrutElms_pft=0._r8
  allocate(SeasonalNonstElms_pft(NumPlantChemElms,JP,JY,JX));  SeasonalNonstElms_pft=0._r8
  allocate(SeedPlantedElm_pft(NumPlantChemElms,JP,JY,JX));    SeedPlantedElm_pft=0._r8
  allocate(AvgCanopyBiomC2Graze_pft(JP,JY,JX));   AvgCanopyBiomC2Graze_pft=0._r8
  end subroutine InitCanopyData

!----------------------------------------------------------------------
  subroutine DestructCanopyData
  use abortutils, only : destroy
  implicit none
  call destroy(fNCLFW_pft)
  call destroy(fPCLFW_pft)
  call destroy(SpecificLeafArea_pft)
  call destroy(LeafC3ChlCperm2LA_pft)
  call destroy(LeafC4ChlCperm2LA_pft)
  call destroy(LeafRubiscoCperm2LA_pft)
  call destroy(LeafPEPCperm2LA_pft)
  call destroy(PARSunlit_pft)
  call destroy(PARSunsha_pft)
  call destroy(CH2OSunlit_pft)
  call destroy(CH2OSunsha_pft)
  call destroy(TFN_Carboxy_pft)
  call destroy(TFN_Oxygen_pft)
  call destroy(TFN_eTranspt_pft)
  call destroy(LeafAreaSunlit_pft)
  call destroy(CanopyVcMaxRubisco25C_pft)
  call destroy(CanopyVoMaxRubisco25C_pft)
  call destroy(CanopyVcMaxPEP25C_pft)
  call destroy(CanopyNLimFactor_brch)
  call destroy(CanopyPLimFactor_brch)
  call destroy(CanopyMassC_pft)
  call destroy(CanopyWaterMassBeg_col)
  call destroy(CanopyWaterMassEnd_col)
  call destroy(HeatCanopy2Dist_col)
  call destroy(QCanopyWat2Dist_col)
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
  call destroy(RCS_pft)
  call destroy(CanPStomaResistH2O_pft)
  call destroy(CanopyMinStomaResistH2O_pft)
  call destroy(CanopyBndlResist_col)
  call destroy(O2I_pft)
  call destroy(LeafIntracellularCO2_pft)
  call destroy(AirConc_pft)
  call destroy(DiffCO2Atmos2Intracel_pft)
  call destroy(CanopyGasCO2_pft)
  call destroy(aquCO2Intraleaf_pft)
  call destroy(O2L_pft)
  call destroy(CO2Solubility_pft)
  call destroy(LeafO2Solubility_pft)
  call destroy(Km4LeafaqCO2_pft)
  call destroy(Km4RubiscoCarboxy_pft)
  call destroy(ChillHours_pft)
  call destroy(Vmax4RubiscoCarboxy_node)
  call destroy(CO2lmtRubiscoCarboxyRate_node)
  call destroy(CO2CompenPoint_node)
  call destroy(LigthSatCarboxyRate_node)
  call destroy(RubiscoCarboxyEff_node)
  call destroy(ProteinCperm2LeafArea_node)
  call destroy(CMassCO2BundleSheath_node)
  call destroy(Vmax4PEPCarboxy_node)
  call destroy(CO2lmtPEPCarboxyRate_node)
  call destroy(LigthSatC4CarboxyRate_node)
  call destroy(C4CarboxyEff_node)
  call destroy(CPOOL4_node)
  call destroy(CMassHCO3BundleSheath_node)
  call destroy(RubiscoActivity_brch)
  call destroy(NutrientCtrlonC4Carboxy_node)
  call destroy(C4PhotosynDowreg_brch)
  call destroy(NetCO2Flx2Canopy_col)
  call destroy(VmaxSpecRubCarboxyRef_pft)
  call destroy(VmaxRubOxyRef_pft)
  call destroy(VmaxPEPCarboxyRef_pft)
  call destroy(XKCO2_pft)
  call destroy(XKO2_pft)
  call destroy(Km4PEPCarboxy_pft)
  call destroy(LeafRubisco2Protein_pft)
  call destroy(LeafPEP2Protein_pft)
  call destroy(SpecLeafChlAct_pft)
  call destroy(LeafC3Chl2Protein_pft)
  call destroy(LeafC4Chl2Protein_pft)
  call destroy(CanopyCi2CaRatio_pft)
  call destroy(RadNet2Canopy_pft)
  call destroy(LWRadCanopy_pft)
  call destroy(RadSWbyCanopy_pft)
  call destroy(RadPARbyCanopy_pft)
  call destroy(FracPARads2Canopy_pft)
  call destroy(TAU_RadThru)
  call destroy(TAU_DirectRTransmit)
  call destroy(FracSWRad2Grnd_col)
  call destroy(RadSWGrnd_col)
  call destroy(LWRadCanGPrev_col)
  call destroy(LWRadGrnd_col)
  call destroy(WatHeldOnCanopy_col)
  call destroy(Prec2Canopy_col)
  call destroy(PrecIntceptByCanopy_col)
  call destroy(EvapTransLHeat_pft)
  call destroy(HeatXAir2PCan_pft)
  call destroy(HeatStorCanopy_pft)
  call destroy(ENGYX_pft)
  call destroy(VHeatCapCanopy_pft)
  call destroy(PSICanopy_pft)
  call destroy(PSICanopyTurg_pft)
  call destroy(PSICanopyOsmo_pft)
  call destroy(CanopyBndlResist_pft)
  call destroy(Transpiration_pft)
  call destroy(VapXAir2Canopy_pft)
  call destroy(CanopyBiomWater_pft)
  call destroy(QVegET_col)
  call destroy(VapXAir2Canopy_col)
  call destroy(CanopyHeatStor_col)
  call destroy(HeatFlx2Canopy_col)
  call destroy(CanopyWat_col)
  call destroy(LWRadCanG_col)
  call destroy(RadSWLeafAlbedo_pft)
  call destroy(RadSWLeafTransmis_pft)
  call destroy(PrecIntcptByCanopy_pft)
  call destroy(WatHeldOnCanopy_pft)
  call destroy(TKC_pft)
  call destroy(TdegCCanopy_pft)
  call destroy(DeltaTKC_pft)
  call destroy(TKCanopy_pft)
  call destroy(CPOOL3_node)
  call destroy(ShootElms_pft)
  call destroy(NetCumElmntFlx2Plant_pft)
  call destroy(tCanLeafC_clyr)
  call destroy(PlantElmAllocMat4Litr)
  call destroy(LeafStrutElms_pft)
  call destroy(PetoleStrutElms_pft)
  call destroy(StalkStrutElms_pft)
  call destroy(CanopySapwoodC_pft)
  call destroy(StalkRsrvElms_pft)
  call destroy(HuskStrutElms_pft)
  call destroy(EarStrutElms_pft)
  call destroy(GrainStrutElms_pft)
  call destroy(CanopyLeafSheathC_pft)
  call destroy(CanopyLeafAreaZ_pft)
  call destroy(CO2NetFix_pft)
  call destroy(RCanMaintDef_CO2_pft)
  call destroy(CanopyLeafCLyr_pft)
  call destroy(CanopyNonstElms_pft)
  call destroy(CanopyNonstElmConc_pft)
  call destroy(CanopyStemAreaZ_pft)
  call destroy(CanopyNodulNonstElms_pft)
  call destroy(ShootNoduleElms_pft)
  call destroy(SapwoodBiomassC_brch)
  call destroy(CanopyNonstElms_brch)
  call destroy(C4PhotoShootNonstC_brch)
  call destroy(CanopyLeafSheathC_brch)
  call destroy(ShootElms_brch)
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
  call destroy(SenecStalkStrutElms_brch)
  call destroy(LeafElmntNode_brch)
  call destroy(PetioleElmntNode_brch)
  call destroy(StructInternodeElms_brch)
  call destroy(LeafLayerElms_node)
  call destroy(CanopyLeafArea_lnode)
  call destroy(LeafProteinC_node)
  call destroy(PetoleProteinC_node)
  call destroy(CanopyNoduleNonstCConc_pft)
  call destroy(GrainSeedBiomCMean_brch)
  call destroy(StandDeadKCompElms_pft)
  call destroy(StandDeadStrutElms_pft)
  call destroy(SeasonalNonstElms_pft)
  call destroy(SeedPlantedElm_pft)
  call destroy(AvgCanopyBiomC2Graze_pft)
  end subroutine DestructCanopyData

end module CanopyDataType
