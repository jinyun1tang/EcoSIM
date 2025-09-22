module PlantAPIData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL   
  use ElmIDMod
  use abortutils, only : destroy
  use EcoSiMParDataMod, only : pltpar
  use TracerIDMod
implicit none
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public

  integer, pointer :: NumGrowthStages              !number of growth stages
  integer, pointer :: MaxNumRootAxes               !maximum number of root layers
  integer, pointer :: MaxNumBranches               !maximum number of plant branches
  integer, pointer :: JP1                          !number of plants
  integer, pointer :: NumOfSkyAzimuthSects1        !number of sectors for the sky azimuth  [0,2*pi]
  integer, pointer :: jcplx                        !number of organo-microbial complexes
  integer, pointer :: NumOfLeafAzimuthSectors1     !number of sectors for the leaf azimuth, [0,pi]
  integer, pointer :: NumCanopyLayers1           !number of canopy layers
  integer, pointer :: JZ1                          !number of soil layers
  integer, pointer :: NumLeafZenithSectors1      !number of sectors for the leaf zenith [0,pi/2]
  integer, pointer :: MaxNodesPerBranch1           !number of canopy nodes
  integer, pointer :: jsken                        !number of kinetic components in litter, PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  integer, pointer :: NumLitterGroups              !number of litter groups nonstructural(0,*),foliar(1,*),non-foliar(jroots,*),stalk(3,*),root(4,*), coarse woody (5,*)
  integer, pointer :: NumOfPlantMorphUnits         !number of organs involved in partition
  integer, pointer :: NumOfPlantLitrCmplxs         !number of plant litter complexes
  integer, pointer :: jroots                       !number of root types, root, myco

  !begin_derived_data_type
  type, public :: plant_siteinfo_type
  real(r8) :: ALAT                                 !latitude,	[degrees north]
  real(r8) :: ATCA                                 !mean annual air temperature, [oC]
  real(r8) :: ALT                                  !altitude of the grid cell, [m]
  real(r8) :: CCO2EI_gperm3                        !initial atmospheric CO2 concentration, [g m-3]
  real(r8) :: CO2EI                                !initial atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: COXYE                                !current atmospheric O2 concentration, [g m-3]
  real(r8) :: SolarNoonHour_col                    !time of solar noon, [h]
  real(r8) :: CO2E                                 !atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: DayLenthPrev                         !daylength of previous day, [h]
  real(r8) :: DayLenthCurrent                      !current daylength of the grid, [h]
  real(r8) :: DayLenthMax_col                      !maximum daylength of the grid, [h]
  real(r8) :: SoilSurfRoughnesst0_col              !initial soil surface roughness height, [m]
  real(r8) :: OXYE                                 !atmospheric O2 concentration, [umol mol-1]
  real(r8) :: PlantPopu_col                        !total plant population, [plants d-2]
  real(r8) :: POROS1                               !top layer soil porosity, [m3 m-3]
  real(r8) :: WindSpeedAtm_col                     !wind speed, [m h-1]
  real(r8) :: WindMesureHeight_col                 !wind speed measurement height, [m]
  real(r8) :: ZEROS2                               !threshold zero for numerical stability,[-]
  real(r8) :: ZEROS                                !threshold zero for numerical stability,[-]
  real(r8) :: ZERO                                 !threshold zero for numerical stability, [-]
  real(r8) :: ZERO2                                !threshold zero for numerical stability,[-]
  integer :: KoppenClimZone                        !Koppen climate zone for the grid,[-]
  integer :: NumActivePlants                       !number of active PFT in the grid, [-]
  integer :: NL                                    !lowest soil layer number,[-]
  integer :: NP0                                   !intitial number of plant species,[-]
  integer :: MaxNumRootLays                        !maximum root layer number,[-]
  integer :: NP                                    !current number of plant species,[-]
  integer :: NU                                    !current soil surface layer number, [-]
  integer :: NK                                    !current hydrologically active layer, [-]
  integer :: DazCurrYear                           !number of days in current year,[-]
  integer :: iYearCurrent                          !current year,[-]
  real(r8) :: QH2OLoss_lnds                        !total subsurface water loss flux over the landscape,	[m3 d-2]

  character(len=16), pointer :: DATAP(:)             => null()    !parameter file name,[-]
  CHARACTER(len=16), pointer :: DATA(:)              => null()    !pft file,[-]
  logical ,     pointer :: flag_active_pft(:)        => null()    !flag for active plants,  [-]
  real (r8), pointer :: PlantElemntStoreLandscape(:) => null()    !total plant element balance,                                              [g d-2]
  real (r8), pointer :: AtmGasc(:)                   => null()    !atmospheric gas concentrations,                                           [g m-3]
  real (r8), pointer :: AREA3(:)                     => null()    !soil cross section area (vertical plane defined by its normal direction), [m2]
  real (r8), pointer :: PlantElmBalCum_pft(:,:)       => null()    !cumulative plant element balance,                                         [g d-2]
  real (r8), pointer :: CumSoilThickness_vr(:)       => null()    !depth to bottom of soil layer from  surface of grid cell,                 [m]
  real (r8), pointer :: PPI_pft(:)                   => null()    !initial plant population,                                                 [plants d-2]
  real (r8), pointer :: PPatSeeding_pft(:)           => null()    !plant population at seeding,                                              [plants d-2]
  real (r8), pointer :: PPX_pft(:)                   => null()    !plant population,                                                         [plants m-2]
  real (r8), pointer :: PlantPopulation_pft(:)       => null()    !plant population,                                         [d-2]
  real (r8), pointer :: CumSoilThickMidL_vr(:)       => null()    !depth to middle of soil layer from  surface of grid cell, [m]
  real (r8), pointer :: FracSoiAsMicP_vr(:)          => null()    !micropore fraction,                                       [-]
  real (r8), pointer :: DLYR3(:)                     => null()    !vertical thickness of soil layer,                         [m]
  real (r8), pointer :: VLWatMicPM_vr(:,:)           => null()    !soil micropore water content,                             [m3 d-2]
  real (r8), pointer :: VLsoiAirPM_vr(:,:)           => null()    !soil air content,                                         [m3 d-2]
  real (r8), pointer :: TortMicPM_vr(:,:)            => null()    !micropore soil tortuosity,                                [m3 m-3]
  real (r8), pointer :: FILMM_vr(:,:)                => null()    !soil water film thickness,                                [m]
  contains
    procedure, public :: Init =>  plt_site_Init
    procedure, public :: Destroy => plt_site_destroy
  end type plant_siteinfo_type

  type, public :: plant_photosyns_type
  integer,  pointer :: iPlantPhotosynthesisType(:)          => null()  !plant photosynthetic type (C3 or C4),[-]  
  real(r8), pointer :: SpecLeafChlAct_pft(:)              => null()  !cholorophyll activity  at 25 oC,                                           [umol g-1 h-1]
  real(r8), pointer :: LeafC3Chl2Protein_pft(:)           => null()  !leaf C3 chlorophyll content,                                               [gC gC-1]
  real(r8), pointer :: LeafPEP2Protein_pft(:)   => null()  !leaf PEP carboxylase content,                                              [gC gC-1]
  real(r8), pointer :: LeafC4Chl2Protein_pft(:)           => null()  !leaf C4 chlorophyll content,                                               [gC gC-1]
  real(r8), pointer :: LeafRubisco2Protein_pft(:)                  => null()  !leaf rubisco content,                                                      [gC gC-1]
  real(r8), pointer :: VmaxPEPCarboxyRef_pft(:)             => null()  !PEP carboxylase activity at 25 oC                                          [umol g-1 h-1]
  real(r8), pointer :: VmaxRubOxyRef_pft(:)                 => null()  !rubisco oxygenase activity  at 25 oC,                                      [umol g-1 h-1]
  real(r8), pointer :: VmaxSpecRubCarboxyRef_pft(:)             => null()  !rubisco carboxylase activity  at 25 oC,                                    [umol g-1 h-1]
  real(r8), pointer :: XKCO2_pft(:)                         => null()  !Km for rubisco carboxylase activity,                                       [uM]
  real(r8), pointer :: XKO2_pft(:)                          => null()  !Km for rubisco oxygenase activity,                                         [uM]
  real(r8), pointer :: RubiscoActivity_brch(:,:)            => null()   !branch down-regulation of CO2 fixation,                                   [-]
  real(r8), pointer :: C4PhotosynDowreg_brch(:,:)           => null()   !down-regulation of C4 photosynthesis,                                     [-]
  real(r8), pointer :: aquCO2Intraleaf_pft(:)               => null()   !leaf aqueous CO2 concentration,                                           [uM]
  real(r8), pointer :: O2I_pft(:)                           => null()   !leaf gaseous O2 concentration,                                            [umol m-3]
  real(r8), pointer :: LeafIntracellularCO2_pft(:)          => null()   !leaf gaseous CO2 concentration,                                           [umol m-3]
  real(r8), pointer :: Km4RubiscoCarboxy_pft(:)             => null()   !leaf aqueous CO2 Km ambient O2,                                           [uM]
  real(r8), pointer :: Km4LeafaqCO2_pft(:)                  => null()   !leaf aqueous CO2 Km no O2,                                                [uM]
  real(r8), pointer :: RCS_pft(:)                           => null()   !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8), pointer :: CanopyCi2CaRatio_pft(:)                    => null()   !Ci:Ca ratio,                                                              [-]
  real(r8), pointer :: H2OCuticleResist_pft(:)              => null()   !maximum stomatal resistance to vapor,                                     [s h-1]
  real(r8), pointer :: ChillHours_pft(:)                    => null()   !chilling effect on CO2 fixation,                                          [-]
  real(r8), pointer :: CO2Solubility_pft(:)                 => null()   !leaf CO2 solubility,                                                      [uM /umol mol-1]
  real(r8), pointer :: CanopyGasCO2_pft(:)                  => null()   !canopy gaesous CO2 concentration,                                         [umol mol-1]
  real(r8), pointer :: O2L_pft(:)                           => null()   !leaf aqueous O2 concentration,                                            [uM]
  real(r8), pointer :: Km4PEPCarboxy_pft(:)                 => null()   !Km for PEP carboxylase activity,                                          [uM]
  real(r8), pointer :: CanopyMinStomaResistH2O_pft(:)         => null()   !canopy minimum stomatal resistance,                          [s m-1]
  real(r8), pointer :: CuticleResist_pft(:)                 => null()   !maximum stomatal resistance to vapor,                        [s m-1]
  real(r8), pointer :: CanPStomaResistH2O_pft(:)            => null()   !canopy stomatal resistance,                                  [h m-1]
  real(r8), pointer :: CanopyBndlResist_pft(:)              => null()   !canopy boundary layer resistance,                            [h m-1]
  real(r8), pointer :: LeafO2Solubility_pft(:)              => null()   !leaf O2 solubility,                                          [uM /umol mol-1]
  real(r8), pointer :: CO2CuticleResist_pft(:)              => null()   !maximum stomatal resistance to CO2,                          [s h-1]
  real(r8), pointer :: DiffCO2Atmos2Intracel_pft(:)         => null()   !gaesous CO2 concentration difference across stomates,        [umol m-3]
  real(r8), pointer :: LeafAreaSunlit_zsec(:,:,:,:,:)        => null()  !leaf irradiated surface area in different leaf sector,       [m2 d-2]
  real(r8), pointer :: CPOOL3_node(:,:,:)                   => null()   !minimum sink strength for nonstructural C transfer,          [g d-2]
  real(r8), pointer :: CPOOL4_node(:,:,:)                   => null()   !leaf nonstructural C4 content in C4 photosynthesis,          [g d-2]
  real(r8), pointer :: CMassCO2BundleSheath_node(:,:,:)     => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: CO2CompenPoint_node(:,:,:)           => null()   !CO2 compensation point,                                      [uM]
  real(r8), pointer :: VoMaxRubiscoRef_brch(:,:)            => null()   !branch maximum rubisco oxygenation rate at reference temperature, [umol g-1 h-1]
  real(r8), pointer :: VcMaxRubiscoRef_brch(:,:)            => null()   !branch maximum rubisco carboxylation rate at reference temperature, [umol g-1 h-1]
  real(r8), pointer :: RubiscoCarboxyEff_node(:,:,:)        => null()   !carboxylation efficiency,                                    [umol umol-1]
  real(r8), pointer :: VcMaxPEPCarboxyRef_brch(:,:)         => null()   !branch reference maximum dark C4 carboxylation rate under saturating CO2, [umol s-1]
  real(r8), pointer :: ElectronTransptJmaxRef_brch(:,:)     => null()   !branch Jmax at reference temperature, [umol e- s-1]
  real(r8), pointer :: C4CarboxyEff_node(:,:,:)             => null()   !C4 carboxylation efficiency,                                 [umol umol-1]
  real(r8), pointer :: LigthSatCarboxyRate_node(:,:,:)      => null()   !maximum light carboxylation rate under saturating CO2,       [umol m-2 s-1]
  real(r8), pointer :: LigthSatC4CarboxyRate_node(:,:,:)    => null()   !maximum  light C4 carboxylation rate under saturating CO2,   [umol m-2 s-1]
  real(r8), pointer :: NutrientCtrlonC4Carboxy_node(:,:,:)  => null()   !down-regulation of C4 photosynthesis,                        [-]
  real(r8), pointer :: CMassHCO3BundleSheath_node(:,:,:)    => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: Vmax4RubiscoCarboxy_node(:,:,:)      => null()   !maximum dark carboxylation rate under saturating CO2,        [umol m-2 s-1]
  real(r8), pointer :: ProteinCperm2LeafArea_node(:,:,:)    => null()   !Protein C per m2 of leaf aera,                               [gC (leaf area m-2)]
  real(r8), pointer :: CO2lmtRubiscoCarboxyRate_node(:,:,:) => null()   !carboxylation rate,                                          [umol m-2 s-1]
  real(r8), pointer :: Vmax4PEPCarboxy_node(:,:,:)          => null()   !maximum dark C4 carboxylation rate under saturating CO2,     [umol m-2 s-1]
  real(r8), pointer :: CO2lmtPEPCarboxyRate_node(:,:,:)     => null()   !C4 carboxylation rate,                                       [umol m-2 s-1]
  real(r8), pointer :: AirConc_pft(:)                       => null()   !total gas concentration,                                     [mol m-3]
  real(r8), pointer :: CanopyVcMaxRubisco25C_pft(:)         => null()   !Canopy VcMax for rubisco carboxylation, [umol h-1 m-2]
  real(r8), pointer :: CanopyVoMaxRubisco25C_pft(:)         => null()   !Canopy VoMax for rubisco oxygenation, [umol h-1 m-2]
  real(r8), pointer :: CanopyVcMaxPEP25C_pft(:)             => null()   !Canopy VcMax in PEP C4 fixation, [umol h-1 m-2]
  real(r8), pointer :: ElectronTransptJmax25C_pft(:)        => null()   !Canopy Jmax at reference temperature, [umol e- s-1 m-2]
  real(r8), pointer :: TFN_Carboxy_pft(:)                   => null()   !temperature dependence of carboxylation, [-]
  real(r8), pointer :: TFN_Oxygen_pft(:)                    => null()   !temperature dependence of oxygenation, [-]
  real(r8), pointer :: TFN_eTranspt_pft(:)                  => null()   !temperature dependence of electron transport, [-]
  real(r8), pointer :: LeafAreaSunlit_pft(:)                => null()   !leaf irradiated surface area, [m2 d-2]    
  real(r8), pointer :: PARSunlit_pft(:)                     => null()   !PAR absorbed by sunlit leaf, [umol m-2 s-1]
  real(r8), pointer :: PARSunsha_pft(:)                     => null()   !PAR absorbed by sun-shaded leaf, [umol m-2 s-1]
  real(r8), pointer :: CH2OSunlit_pft(:)                    => null()   !carbon fixation by sun-lit leaf, [gC d-2 h-1]
  real(r8), pointer :: CH2OSunsha_pft(:)                    => null()   !carbon fixation by sun-shaded leaf, [gC d-2 h-1]    

  contains
    procedure, public :: Init    =>  plt_photo_init
    procedure, public :: Destroy => plt_photo_destroy
  end type plant_photosyns_type

  type, public ::  plant_radiation_type
  real(r8) :: TotSineSkyAngles_grd           !sine of sky angles,[-]
  real(r8) :: SoilAlbedo                     !soil albedo,[-]
  real(r8) :: SurfAlbedo_col                 !Surface albedo,[-]
  real(r8) :: RadPARSolarBeam_col            !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8) :: RadSWDiffus_col                !diffuse shortwave radiation, [W m-2]
  real(r8) :: RadSWDirect_col                !direct shortwave radiation, [W m-2]
  real(r8) :: RadPARDiffus_col               !diffuse PAR, [umol m-2 s-1]
  real(r8) :: RadDirectPAR_col               !direct PAR, [umol m-2 s-1]
  real(r8) :: Eco_NetRad_col                 !ecosystem net radiation, [MJ d-2 h-1]
  real(r8) :: RadSWSolarBeam_col             !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8) :: FracSWRad2Grnd_col             !fraction of radiation intercepted by ground surface, [-]
  real(r8) :: RadSWGrnd_col                  !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8) :: SineGrndSlope_col              !sine of slope, [-]
  real(r8) :: GroundSurfAzimuth_col          !azimuth of slope, [-]
  real(r8) :: CosineGrndSlope_col            !cosine of slope, [-]
  real(r8) :: LWRadGrnd_col                      !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8) :: LWRadSky_col                   !sky longwave radiation , [MJ d-2 h-1]
  real(r8) :: SineSunInclAnglNxtHour_col     !sine of solar angle next hour, [-]
  real(r8) :: SineSunInclAngle_col           !sine of solar angle, [-]

  integer,  pointer :: iScatteringDiffus(:,:,:)  => null() !flag for calculating backscattering of radiation in canopy,[-]
  real(r8), pointer :: RadSWLeafAlbedo_pft(:)    => null() !canopy shortwave albedo,                           [-]
  real(r8), pointer :: CanopyPARalbedo_pft(:)    => null() !canopy PAR albedo,                                 [-]
  real(r8), pointer :: TAU_DirectRTransmit(:)    => null() !fraction of radiation intercepted by canopy layer, [-]
  real(r8), pointer :: TAU_RadThru(:)            => null() !fraction of radiation transmitted by canopy layer, [-]
  real(r8), pointer :: LWRadCanopy_pft(:)        => null() !canopy longwave radiation,                         [MJ d-2 h-1]
  real(r8), pointer :: RadSWbyCanopy_pft(:)      => null() !canopy absorbed shortwave radiation,               [MJ d-2 h-1]
  real(r8), pointer :: OMEGX(:,:,:)              => null() !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation,[-]
  real(r8), pointer :: OMEGAG(:)                 => null() !sine of solar beam on leaf surface,                [-]
  real(r8), pointer :: OMEGA(:,:,:)              => null() !sine of indirect sky radiation on leaf surface,[-]
  real(r8), pointer :: SineLeafAngle(:)          => null() !sine of leaf angle,[-]
  real(r8), pointer :: CosineLeafAngle(:)        => null() !cosine of leaf angle,[-]
  real(r8), pointer :: RadNet2Canopy_pft(:)      => null() !canopy net radiation,                              [MJ d-2 h-1]
  real(r8), pointer :: LeafSWabsorpty_pft(:)     => null() !canopy shortwave absorptivity,                     [-]
  real(r8), pointer :: LeafPARabsorpty_pft(:)    => null() !canopy PAR absorptivity,[-]
  real(r8), pointer :: RadSWLeafTransmis_pft(:)  => null() !canopy shortwave transmissivity,                   [-]
  real(r8), pointer :: RadPARLeafTransmis_pft(:) => null() !canopy PAR transmissivity,                         [-]
  real(r8), pointer :: RadPARbyCanopy_pft(:)     => null() !canopy absorbed PAR,                               [umol m-2 s-1]
  real(r8), pointer :: FracPARads2Canopy_pft(:)  => null() !fraction of incoming PAR absorbed by canopy,       [-]
  real(r8), pointer :: RadTotPAR_zsec(:,:,:,:)      => null()     !direct incoming PAR,                           [umol m-2 s-1]
  real(r8), pointer :: RadDifPAR_zsec(:,:,:,:)   => null()  !diffuse incoming PAR,                             [umol m-2 s-1]
  contains
    procedure, public :: Init    => plt_rad_init
    procedure, public :: Destroy => plt_rad_destroy
  end type plant_radiation_type

  type, public :: plant_morph_type
  real(r8) :: LeafStalkArea_col                       !stalk area of combined, each PFT canopy,[m^2 d-2]
  real(r8) :: CanopyLeafArea_col                      !grid canopy leaf area, [m2 d-2]
  real(r8) :: StemArea_col                            !grid canopy stem area, [m2 d-2]
  real(r8) :: CanopyHeight_col                        !canopy height , [m]
  real(r8), pointer :: tlai_day_pft(:)                 => null() !prescribed leaf area, [m2 m-2]
  real(r8), pointer :: tsai_day_pft(:)                 => null() !prescribed stem area, [m2 m-2]
  real(r8), pointer :: PARTS_brch(:,:,:)               => null() !fraction of C allocated to each morph unit,                                 [-]
  real(r8), pointer :: RootVolPerMassC_pft(:,:)        => null() !root volume:mass ratio,                                                     [m3 g-1]
  real(r8), pointer :: RootPorosity_pft(:,:)           => null() !root porosity,                                                              [m3 m-3]
  real(r8), pointer :: Root2ndXSecArea_pft(:,:)        => null() !root  cross-sectional area  secondary axes,                                 [m2]
  real(r8), pointer :: Root1stXSecArea_pft(:,:)        => null() !root cross-sectional area primary axes,                                     [m2]
  real(r8), pointer :: Root1stMaxRadius1_pft(:,:)      => null() !root diameter primary axes,                                                 [m]
  real(r8), pointer :: Root2ndMaxRadius1_pft(:,:)      => null() !root diameter secondary axes,                                               [m]
  real(r8), pointer :: SeedCMass_pft(:)                => null() !grain size at seeding,                                                      [g]
  real(r8), pointer :: RootPoreTortu4Gas_pft(:,:)      => null() !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8), pointer :: Root1stLen_rpvr(:,:,:,:)        => null() !root layer length primary axes,                                             [m d-2]
  real(r8), pointer :: Root2ndLen_rpvr(:,:,:,:)        => null() !root layer length secondary axes,                                           [m d-2]
  real(r8), pointer :: RootTotLenPerPlant_pvr(:,:,:)   => null() !total root length per plant,                                                [m p-1]
  real(r8), pointer :: RootLenPerPlant_pvr(:,:,:)      => null() !fine root length per plant, [m p-1]
  real(r8), pointer :: Root2ndEffLen4uptk_rpvr(:,:,:)  => null() !Layer effective root length four resource uptake, [m]
  real(r8), pointer :: Root1stSpecLen_pft(:,:)         => null() !specific root length primary axes,                                          [m g-1]
  real(r8), pointer :: Root2ndSpecLen_pft(:,:)         => null() !specific root length secondary axes,                                        [m g-1]
  real(r8), pointer :: Root2ndXNum_rpvr(:,:,:,:)       => null() !root layer number secondary axes,                                           [d-2]
  real(r8), pointer :: Root1stDepz_pft(:,:,:)          => null() !root layer depth,                                                           [m]
  real(r8), pointer :: ClumpFactorInit_pft(:)          => null() !initial clumping factor for self-shading in canopy layer,                   [-]
  real(r8), pointer :: ClumpFactorNow_pft(:)           => null() !clumping factor for self-shading in canopy layer at current LAI,            [-]
  real(r8), pointer :: RootBranchFreq_pft(:)           => null() !root brancing frequency,                                                    [m-1]
  real(r8), pointer :: HypoctoHeight_pft(:)            => null() !cotyledon height,                                                           [m]
  real(r8), pointer :: CanopyHeight4WatUptake_pft(:)   => null() !canopy height,                                                              [m]
  real(r8), pointer :: CanopyHeightZ_col(:)            => null() !canopy layer height,                                                        [m]
  real(r8), pointer :: CanopyStemArea_pft(:)           => null() !plant stem area,                                                            [m2 d-2]
  real(r8), pointer :: CanopyLeafArea_pft(:)           => null() !plant canopy leaf area,                                                     [m2 d-2]
  real(r8), pointer :: ShootNodeNum_brch(:,:)          => null() !shoot node number,                                                          [-]
  real(r8), pointer :: NodeNum2InitFloral_brch(:,:)    => null() !shoot node number at floral initiation,                                     [-]
  real(r8), pointer :: NodeNumberAtAnthesis_brch(:,:)  => null() !shoot node number at anthesis,                                              [-]
  real(r8), pointer :: SineBranchAngle_pft(:)          => null() !branching angle,                                                            [degree from horizontal]
  real(r8), pointer :: LeafAreaZsec_brch(:,:,:,:,:)    => null() !leaf surface area,                                                          [m2 d-2]
  real(r8), pointer :: PotentialSeedSites_brch(:,:)    => null() !branch potential grain number,                                              [d-2]
  real(r8), pointer :: SeedSitesSet_brch(:,:)            => null() !branch grain number,                                                        [d-2]
  real(r8), pointer :: CanPBranchHeight(:,:)           => null() !branch height,                                                              [m]
  real(r8), pointer :: LeafAreaDying_brch(:,:)         => null() !branch leaf area,                                                           [m2 d-2]
  real(r8), pointer :: LeafAreaLive_brch(:,:)          => null() !branch leaf area,                                                           [m2 d-2]
  real(r8), pointer :: SinePetioleAngle_pft(:)         => null() !sheath angle,                                                               [degree from horizontal]
  real(r8), pointer :: LeafAngleClass_pft(:,:)         => null() !fractionction of leaves in different angle classes,                         [-]
  real(r8), pointer :: CanopyStemAreaZ_pft(:,:)        => null() !plant canopy layer stem area,                                               [m2 d-2]
  real(r8), pointer :: CanopyLeafAreaZ_pft(:,:)        => null() !canopy layer leaf area,                                                     [m2 d-2]
  real(r8), pointer :: LeafArea_node(:,:,:)        => null() !leaf area,                                                                  [m2 d-2]
  real(r8), pointer :: CanopySeedNum_pft(:)            => null() !canopy grain number,                                                        [d-2]
  real(r8), pointer :: SeedDepth_pft(:)                => null() !seeding depth,                                                              [m]
  real(r8), pointer :: PlantinDepz_pft(:)              => null() !planting depth,                                                             [m]
  real(r8), pointer :: SeedMeanLen_pft(:)              => null() !seed length,                                                                [m]
  real(r8), pointer :: SeedVolumeMean_pft(:)           => null() !seed volume,                                                                [m3 ]
  real(r8), pointer :: SeedAreaMean_pft(:)             => null() !seed surface area,                                                          [m2]
  real(r8), pointer :: CanopyStemAareZ_col(:)          => null() !total stem area,                                                            [m2 d-2]
  real(r8), pointer :: PetoLen2Mass_pft(:)             => null() !petiole length:mass during growth,                                          [m gC-1]
  real(r8), pointer :: NodeLenPergC_pft(:)             => null() !internode length:mass during growth,                                        [m gC-1]
  real(r8), pointer :: SLA1_pft(:)                     => null() !leaf area:mass during growth,                                               [m2 gC-1]
  real(r8), pointer :: CanopyLeafAareZ_col(:)          => null() !total leaf area,                                                            [m2 d-2]
  real(r8), pointer :: LeafStalkArea_pft(:)            => null() !plant leaf+stem/stalk area,                                                 [m2 d-2]
  real(r8), pointer :: StalkNodeVertLength_brch(:,:,:) => null() !internode height,                                                           [m]
  real(r8), pointer :: PetoleLensNode_brch(:,:,:)      => null() !sheath height,                                                              [m]
  real(r8), pointer :: StalkNodeHeight_brch(:,:,:)  => null() !internode height,                                                           [m]
  real(r8), pointer :: StemAreaZsec_brch(:,:,:,:)      => null() !stem surface area,                                                          [m2 d-2]
  real(r8), pointer :: CanopyLeafArea_lnode(:,:,:,:)   => null() !layer/node/branch leaf area,                                                [m2 d-2]
  real(r8), pointer :: CanopyStalkArea_lbrch(:,:,:)    => null() !plant canopy layer branch stem area,                                        [m2 d-2]
  real(r8), pointer :: ClumpFactor_pft(:)              => null() !clumping factor for self-shading in canopy layer,                           [-]
  real(r8), pointer :: ShootNodeNumAtPlanting_pft(:)   => null() !number of nodes in seed,                                                    [-]
  real(r8), pointer :: CanopyHeight_pft(:)             => null() !canopy height,                                                              [m]
  real(r8), pointer :: TreeRingAveRadius_pft(:)        => null() !tree ring radius,[m]
  integer , pointer :: iPlantGrainType_pft(:)        => null() !grain type (below or above-ground),[-]
  integer,  pointer :: iPlantNfixType_pft(:)          => null() !N2 fixation type,[-]
  integer,  pointer :: Myco_pft(:)                    => null() !mycorrhizal type (no or yes),[-]
  integer,  pointer :: MainBranchNum_pft(:)           => null() !number of main branch,[-]
  integer,  pointer :: MaxSoiL4Root_pft(:)            => null() !maximum soil layer number for all root axes,[-]
  integer,  pointer :: NMaxRootBotLayer_pft(:)         => null() !maximum soil layer number for all root axes, [-]
  integer,  pointer :: NumPrimeRootAxes_pft(:)             => null() !root primary axis number,[-]
  integer,  pointer :: NumCogrowthNode_pft(:)         => null() !number of concurrently growing nodes,[-]
  integer,  pointer :: BranchNumber_pft(:)            => null() !main branch numeric id,[-]
  integer,  pointer :: NumOfBranches_pft(:)           => null() !number of branches,[-]
  integer,  pointer :: NIXBotRootLayer_raxes(:,:)      => null() !maximum soil layer number for root axes,     [-]
  integer,  pointer :: KMinNumLeaf4GroAlloc_brch(:,:) => null() !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION,[-]
  integer,  pointer :: BranchNumerID_brch(:,:)         => null() !branch meric id,                             [-]
  integer,  pointer :: NGTopRootLayer_pft(:)          => null() !soil layer at planting depth,                [-]
  integer,  pointer :: KLeafNumber_brch(:,:)          => null() !leaf number,                                 [-]
  
  real(r8), pointer :: NumOfLeaves_brch(:,:)          => null() !leaf number,                          [-]
  real(r8), pointer :: GrothStalkMaxSeedSites_pft(:)     => null() !maximum grain node number per branch, [-]
  real(r8), pointer :: MaxSeedNumPerSite_pft(:)       => null() !maximum grain number per node,        [-]
  real(r8), pointer :: rLen2WidthLeaf_pft(:)          => null() !leaf length:width ratio,              [-]
  real(r8), pointer :: SeedCMassMax_pft(:)            => null() !maximum grain size,                   [g]
  real(r8), pointer :: Root1stRadius_pvr(:,:,:)       => null() !root layer diameter primary axes,     [m]
  real(r8), pointer :: Root2ndRadius_rpvr(:,:,:)       => null() !root layer diameter secondary axes,   [m]
  real(r8), pointer :: RootRaidus_rpft(:,:)           => null() !root internal radius,                 [m]
  real(r8), pointer :: Root1stMaxRadius_pft(:,:)      => null() !maximum radius of primary roots,      [m]
  real(r8), pointer :: Root2ndMaxRadius_pft(:,:)      => null() !maximum radius of secondary roots,    [m]
  real(r8), pointer :: RootRadialResist_pft(:,:)      => null() !root radial resistivity,              [MPa h m-2]
  real(r8), pointer :: RootAxialResist_pft(:,:)       => null() !root axial resistivity,               [MPa h m-4]
  real(r8), pointer :: totRootLenDens_vr(:)           => null() !total root length density,            [m m-3]
  real(r8), pointer :: Root1stXNumL_rpvr(:,:,:)        => null() !root layer number primary axes,       [d-2]
  real(r8), pointer :: Root2ndXNumL_rpvr(:,:,:)         => null() !root layer number axes,               [d-2]
  real(r8), pointer :: RootLenDensPerPlant_pvr(:,:,:) => null() !layer root length density,            [m m-3]
  real(r8), pointer :: RootPoreVol_rpvr(:,:,:)        => null() !root layer volume air,                [m2 d-2]
  real(r8), pointer :: RootVH2O_pvr(:,:,:)            => null() !root layer volume water,              [m2 d-2]
  real(r8), pointer :: RootAreaPerPlant_pvr(:,:,:)    => null() !layer root area per plant,            [m2 p-1]
  real(r8), pointer :: RootArea1stPP_pvr(:,:,:)       => null() !layer 1st root area per plant, [m2 plant-1]
  real(r8), pointer :: RootArea2ndPP_pvr(:,:,:)       => null() !layer 2nd root area per plant, [m2 plant-1]

  contains
    procedure, public :: Init    => plt_morph_init
    procedure, public :: Destroy => plt_morph_destroy
  end type plant_morph_type

  type, public :: plant_pheno_type
  real(r8), pointer :: fTgrowRootP_vr(:,:)                => null()     !root layer temperature growth functiom,                              [-]
  real(r8), pointer :: ShootRootNonstElmConduts_pft(:)   => null()     !shoot-root rate constant for nonstructural C exchange,               [h-1]
  real(r8), pointer :: GrainFillRate25C_pft(:)            => null()     !maximum rate of fill per grain,                                      [g h-1]
  real(r8), pointer :: TempOffset_pft(:)                  => null()     !adjustment of Arhhenius curves for plant thermal acclimation,        [oC]
  real(r8), pointer :: PlantO2Stress_pft(:)               => null()     !plant O2 stress indicator,                                           [-]
  real(r8), pointer :: NonstCMinConc2InitBranch_pft(:)    => null()     !branch nonstructural C content required for new branch,              [gC gC-1]
  real(r8), pointer :: MinNonstC2InitRoot_pft(:)          => null()     !threshold root nonstructural C content for initiating new root axis, [gC gC-1]
  real(r8), pointer :: LeafElmntRemobFlx_brch(:,:,:)      => null()    !element translocated from leaf during senescence,                     [g d-2 h-1]
  real(r8), pointer :: PetioleChemElmRemobFlx_brch(:,:,:) => null()    !element translocated from sheath during senescence,                   [g d-2 h-1]
  real(r8), pointer :: TC4LeafOut_pft(:)                  => null()     !threshold temperature for spring leafout/dehardening,                [oC]
  real(r8), pointer :: TCGroth_pft(:)                     => null()     !canopy growth temperature,                                           [oC]
  real(r8), pointer :: TC4LeafOff_pft(:)                  => null()     !threshold temperature for autumn leafoff/hardening,                  [oC]
  real(r8), pointer :: TKGroth_pft(:)                     => null()     !canopy growth temperature,                                           [K]
  real(r8), pointer :: fTCanopyGroth_pft(:)               => null()     !canopy temperature growth function,                                  [-]
  real(r8), pointer :: HoursTooLowPsiCan_pft(:)           => null()     !canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY), [h]
  real(r8), pointer :: TCChill4Seed_pft(:)           => null()     !temperature below which seed set is adversely affected, [oC]
  real(r8), pointer :: rPlantThermoAdaptZone_pft(:)  => null()     !plant thermal adaptation zone,                          [-]
  real(r8), pointer :: PlantInitThermoAdaptZone_pft(:)   => null()     !initial plant thermal adaptation zone,                  [-]
  real(r8), pointer :: HighTempLimitSeed_pft(:)      => null()     !temperature above which seed set is adversely affected, [oC]
  real(r8), pointer :: SeedTempSens_pft(:)           => null()     !sensitivity to HTC (seeds oC-1 above HTC),[oC-1]
  real(r8), pointer :: NetCumElmntFlx2Plant_pft(:,:) => null()     !effect of canopy element status on seed set,            [-]
  real(r8), pointer :: RefNodeInitRate_pft(:)        => null()     !rate of node initiation,                                [h-1 at 25 oC]
  real(r8), pointer :: RefLeafAppearRate_pft(:)      => null()     !rate of leaf initiation,                                [h-1 at 25 oC]
  real(r8), pointer :: CriticPhotoPeriod_pft(:)      => null()     !critical daylength for phenological progress,           [h]
  real(r8), pointer :: PhotoPeriodSens_pft(:)        => null()     !difference between current and critical daylengths used to calculate  phenological progress, [h]

  integer,  pointer :: iPlantState_pft(:)                => null()  !flag for species death, [-]
  integer,  pointer :: IsPlantActive_pft(:)              => null()  !flag for living pft, [-]
  integer,  pointer :: iPlantBranchState_brch(:,:)       => null()  !flag to detect branch death,                      [-]
  integer,  pointer :: doRemobilization_brch(:,:)        => null()  !branch phenology flag,                            [-]
  integer,  pointer :: doPlantLeaveOff_brch(:,:)         => null()  !branch phenology flag,                            [-]
  integer,  pointer :: doPlantLeafOut_brch(:,:)          => null()  !branch phenology flag,                            [-]
  integer,  pointer :: doInitLeafOut_brch(:,:)           => null()  !branch phenology flag,                            [-]
  logical,  pointer :: doReSeed_pft(:)                   => null()  !flag to do annual plant reseeding, [-]
  integer,  pointer :: doSenescence_brch(:,:)            => null()  !branch phenology flag,                            [-]
  integer,  pointer :: Prep4Literfall_brch(:,:)          => null()  !branch phenology flag,                            [-]
  integer,  pointer :: Hours4LiterfalAftMature_brch(:,:) => null()  !branch phenology flag,                            [h]
  integer,  pointer :: KHiestGroLeafNode_brch(:,:)       => null()  !leaf growth stage counter,                        [-]
  integer,  pointer :: iPlantPhenolType_pft(:)           => null()  !climate signal for phenological progress: none,   temperature, water stress,[-]
  integer,  pointer :: iPlantTurnoverPattern_pft(:)      => null()  !phenologically-driven above-ground turnover: all, foliar only, none,[-]
  integer,  pointer :: iPlantShootState_pft(:)           => null()  !flag to detect canopy death,[-]
  integer,  pointer :: iPlantPhenolPattern_pft(:)        => null()  !plant growth habit: annual or perennial,[-]
  integer,  pointer :: iPlantRootState_pft(:)            => null()  !flag to detect root system death,[-]
  integer,  pointer :: iPlantDevelopPattern_pft(:)       => null()  !plant growth habit (determinate or indeterminate),[-]
  integer,  pointer :: iPlantPhotoperiodType_pft(:)      => null()  !photoperiod type (neutral, long day, short day),[-]
  integer,  pointer :: iPlantRootProfile_pft(:)          => null()  !plant growth type (vascular, non-vascular),[-]
  integer,  pointer :: doInitPlant_pft(:)                => null()  !PFT initialization flag:0=no,1=yes,[-]
  integer,  pointer :: KLowestGroLeafNode_brch(:,:)      => null()  !leaf growth stage counter,                        [-]
  integer,  pointer :: iPlantCalendar_brch(:,:,:)        => null()  !plant growth stage,                               [-]

  real(r8), pointer :: fRootGrowPSISense_pvr(:,:,:)           => null()  !water stress to plant root growth, [-]
  real(r8), pointer :: TotalNodeNumNormByMatgrp_brch(:,:)     => null()  !normalized node number during vegetative growth stages,                      [-]
  real(r8), pointer :: TotReproNodeNumNormByMatrgrp_brch(:,:) => null()  !normalized node number during reproductive growth stages,                    [-]
  real(r8), pointer :: LeafNumberAtFloralInit_brch(:,:)       => null()  !leaf number at floral initiation,                                            [-]
  real(r8), pointer :: Hours4LenthenPhotoPeriod_brch(:,:)     => null()  !initial heat requirement for spring leafout/dehardening,                     [h]
  real(r8), pointer :: Hours4ShortenPhotoPeriod_brch(:,:)     => null()  !initial cold requirement for autumn leafoff/hardening,                       [h]
  real(r8), pointer :: Hours4Leafout_brch(:,:)                => null()  !heat requirement for spring leafout/dehardening,                             [h]
  real(r8), pointer :: HourReq4LeafOut_brch(:,:)              => null()  !hours above threshold temperature required for spring leafout/dehardening,   [-]
  real(r8), pointer :: Hours4LeafOff_brch(:,:)                => null()  !cold requirement for autumn leafoff/hardening,                               [h]
  real(r8), pointer :: HourReq4LeafOff_brch(:,:)              => null()  !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8), pointer :: Hours2LeafOut_brch(:,:)                => null()  !counter for mobilizing nonstructural C during spring leafout/dehardening,    [h]
  real(r8), pointer :: HoursDoingRemob_brch(:,:)              => null()  !counter for mobilizing nonstructural C during autumn leafoff/hardening,      [h]
  real(r8), pointer :: HourlyNodeNumNormByMatgrp_brch(:,:)    => null()  !gain in normalized node number during vegetative growth stages,              [h-1]
  real(r8), pointer :: dReproNodeNumNormByMatG_brch(:,:)      => null()  !gain in normalized node number during reproductive growth stages,            [h-1]
  real(r8), pointer :: MatureGroup_brch(:,:)                  => null()  !plant maturity group,                                                        [-]
  real(r8), pointer :: NodeNumNormByMatgrp_brch(:,:)          => null()  !normalized node number during vegetative growth stages,                      [-]
  real(r8), pointer :: ReprodNodeNumNormByMatrgrp_brch(:,:)   => null()  !normalized node number during reproductive growth stages,                    [-]
  real(r8), pointer :: HourFailGrainFill_brch(:,:)            => null()  !flag to detect physiological maturity from  grain fill,                      [-]
  real(r8), pointer :: MatureGroup_pft(:)                     => null()  !acclimated plant maturity group,                                             [-]
  contains
    procedure, public :: Init    =>  plt_pheno_init
    procedure, public :: Destroy =>  plt_pheno_destroy
  end type plant_pheno_type

  type, public :: plant_soilchem_type
  real(r8), pointer :: FracBulkSOMC_vr(:,:)           => null()  !fraction of total organic C in complex,       [-]
  real(r8), pointer :: PlantElmAllocMat4Litr(:,:,:,:)      => null() !litter kinetic fraction,                       [-]
  real(r8), pointer :: TScal4Difsvity_vr(:)           => null()  !temperature effect on diffusivity,[-]
  real(r8), pointer :: FracAirFilledSoilPoreM_vr(:,:) => null()  !soil air-filled porosity,                     [m3 m-3]
  real(r8), pointer :: DiffusivitySolutEffM_vr(:,:)   => null()  !coefficient for dissolution - volatilization, [-]
  real(r8), pointer :: SoilResit4RootPentrate_vr(:)   => null()  !soil hydraulic resistance,                    [MPa h m-2]
  real(r8), pointer :: SoilBulkDensity_vr(:)          => null()  !soil bulk density,                            [Mg m-3]
  real(r8), pointer :: trc_solcl_vr(:,:)              => null()  !aqueous tracer concentration, [g m-3]
  real(r8), pointer :: trcg_gascl_vr(:,:)             => null()  !gaseous tracer concentration, [g m-3]
  real(r8), pointer :: CSoilOrgM_vr(:,:)              => null()  !soil organic C content, [gC kg soil-1]
  real(r8), pointer :: HYCDMicP4RootUptake_vr(:) => null()  !soil micropore hydraulic conductivity for root water uptake, [m MPa-1 h-1]
  real(r8), pointer :: GasDifc_vr(:,:)                => null()  !gaseous diffusivity, [m2 h-1]
  real(r8), pointer :: SoluteDifusvty_vr(:,:)         => null()  !aqueous diffusivity, [m2 h-1]
  real(r8), pointer :: trcg_gasml_vr(:,:)             => null()  !gas layer mass, [g d-2]
  real(r8), pointer :: GasSolbility_vr(:,:)           => null()  !gas solubility,                                [m3 m-3]
  real(r8), pointer :: THETW_vr(:)                    => null()  !volumetric water content, [m3 m-3]
  real(r8), pointer :: SoilWatAirDry_vr(:)            => null()  !air-dry water content,                        [m3 m-3]
  real(r8), pointer :: VLSoilPoreMicP_vr(:)           => null()  !volume of soil layer,	[m3 d-2]
  real(r8), pointer :: trcs_VLN_vr(:,:)               => null()  !effective relative tracer volume, [-]
  real(r8), pointer :: VLSoilMicP_vr(:)               => null()  !total micropore volume in layer, [m3 d-2]
  real(r8), pointer :: VLiceMicP_vr(:)                => null()  !soil micropore ice content,   [m3 d-2]
  real(r8), pointer :: VLWatMicP_vr(:)                => null()  !soil micropore water content, [m3 d-2]
  real(r8), pointer :: VLMicP_vr(:)                   => null()  !total volume in micropores, [m3 d-2]
  real(r8), pointer :: DOM_MicP_vr(:,:,:)             => null()  !dissolved organic matter in micropore,	[g d-2]
  real(r8), pointer :: DOM_MicP_drib_vr(:,:,:)        => null()  !dribbling flux for micropore dom,[g d-2]
  real(r8), pointer :: trcs_solml_vr(:,:)             => null() !aqueous tracer, [g d-2]
  contains
    procedure, public :: Init => plt_soilchem_init
    procedure, public :: Destroy => plt_soilchem_destroy
  end type plant_soilchem_type

  type, public :: plant_allometry_type
  real(r8), pointer :: CPRTS_pft(:)                     => null()  !root P:C ratio x root growth yield,                      [-]
  real(r8), pointer :: CNRTS_pft(:)                     => null()  !root N:C ratio x root growth yield,                      [-]
  real(r8), pointer :: rNCNodule_pft(:)                 => null()  !nodule N:C ratio,                                        [gN gC-1]
  real(r8), pointer :: rPCNoduler_pft(:)                 => null()  !nodule P:C ratio,                                        [gP gC-1]
  real(r8), pointer :: rNCRoot_pft(:)                   => null()  !root N:C ratio,                                          [gN gC-1]
  real(r8), pointer :: rPCRootr_pft(:)                   => null()  !root P:C ratio,                                          [gP gC-1]
  real(r8), pointer :: rProteinC2N_pft(:)             => null()  !C:N ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: rProteinC2P_pft(:)             => null()  !C:P ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: NoduGrowthYield_pft(:)           => null()  !nodule growth yield,                                     [g g-1]
  real(r8), pointer :: RootBiomGrosYld_pft(:)           => null()  !root growth yield,                                       [g g-1]
  real(r8), pointer :: rPCEar_pft(:)                    => null()  !ear P:C ratio,                                           [gP gC-1]
  real(r8), pointer :: PetioleBiomGrowthYld_pft(:)      => null()  !sheath growth yield,                                     [g g-1]
  real(r8), pointer :: rPCHusk_pft(:)                   => null()  !husk P:C ratio,                                          [gP gC-1]
  real(r8), pointer :: StalkBiomGrowthYld_pft(:)        => null()  !stalk growth yield,                                      [gC gC-1]
  real(r8), pointer :: HuskBiomGrowthYld_pft(:)         => null()  !husk growth yield,                                       [gC gC-1]
  real(r8), pointer :: ReserveBiomGrowthYld_pft(:)      => null()  !reserve growth yield,                                    [gC gC-1]
  real(r8), pointer :: GrainBiomGrowthYld_pft(:)        => null()  !grain growth yield,                                      [gC gC-1]
  real(r8), pointer :: EarBiomGrowthYld_pft(:)          => null()  !ear growth yield,                                        [gC gC-1]
  real(r8), pointer :: rNCHusk_pft(:)                   => null()  !husk N:C ratio,                                          [gN gC-1]
  real(r8), pointer :: rNCReserve_pft(:)                => null()  !reserve N:C ratio,                                       [gN gC-1]
  real(r8), pointer :: rNCEar_pft(:)                    => null()  !ear N:C ratio,                                           [gN gC-1]
  real(r8), pointer :: rPCReserve_pft(:)                => null()  !reserve P:C ratio,                                       [gP gC-1]
  real(r8), pointer :: rPCGrain_pft(:)                      => null()  !grain P:C ratio,                                         [gP gP-1]
  real(r8), pointer :: rNCStalk_pft(:)                  => null()  !stalk N:C ratio,                                         [gN gC-1]
  real(r8), pointer :: FracShootLeafAlloc2Litr(:,:)  => null()  !woody element allocation, [-]
  real(r8), pointer :: FracShootPetolAlloc2Litr(:,:) => null()  !leaf element allocation,[-]
  real(r8), pointer :: FracRootElmAlloc2Litr(:,:)       => null()  !C woody fraction in root,[-]
  real(r8), pointer :: FracWoodStalkElmAlloc2Litr(:,:)  => null()  !woody element allocation,[-]
  real(r8), pointer :: LeafBiomGrowthYld_pft(:)         => null()  !leaf growth yield,                                       [g g-1]
  real(r8), pointer :: rNCGrain_pft(:)                      => null()  !grain N:C ratio,                                         [g g-1]
  real(r8), pointer :: rPCLeaf_pft(:)                      => null()  !maximum leaf P:C ratio,                                  [g g-1]
  real(r8), pointer :: rPCSheath_pft(:)                     => null()  !sheath P:C ratio,                                        [g g-1]
  real(r8), pointer :: rNCSheath_pft(:)                     => null()  !sheath N:C ratio,                                        [g g-1]
  real(r8), pointer :: rPCStalk_pft(:)                  => null()  !stalk P:C ratio,                                         [g g-1]
  real(r8), pointer :: rNCLeaf_pft(:)                      => null()  !maximum leaf N:C ratio,                                  [g g-1]
  real(r8), pointer :: GrainSeedBiomCMean_brch(:,:)     => null()  !maximum grain C during grain fill,                       [g d-2]
  real(r8), pointer :: FracGroth2Node_pft(:)            => null()  !parameter for allocation of growth to nodes,             [-]
  real(r8), pointer ::RootFracRemobilizableBiom_pft(:)      => null()  !fraction of remobilizable nonstructural biomass in root, [-]

  contains
    procedure, public :: Init => plt_allom_init
    procedure, public :: Destroy => plt_allom_destroy
  end type plant_allometry_type

  type, public :: plant_biom_type
  real(r8), pointer :: TotBegVegE_pft(:,:)                  => null()    !total vegetation biomass at the beginning of time step, [g d-2]
  real(r8), pointer :: TotEndVegE_pft(:,:)                  => null()    !total vegetation biomass at the end of time step, [g d-2]
  real(r8), pointer :: StomatalStress_pft(:)                => null()    !stomatal stress from root turgor (0-1),             [-]
  real(r8), pointer :: RootMycoMassElm_pvr(:,:,:,:)         => null()    !root biomass in chemical elements,                  [g d-2]
  real(r8), pointer :: StandingDeadStrutElms_col(:)         => null()    !total standing dead biomass chemical element,       [g d-2]
  real(r8), pointer :: ZERO4LeafVar_pft(:)                  => null()    !threshold zero for leaf calculation,                [-]
  real(r8), pointer :: LeafProteinC_brch(:,:)               => null()    !Protein C for the branches, [gC protein d-2]
  real(r8), pointer :: LeafProteinCperm2LA_pft(:)           => null()    !Protein C for the plant, [gC protein m-2 leaf area]
  real(r8), pointer :: ZERO4Groth_pft(:)                    => null()    !threshold zero for plang growth calculation,        [-]
  real(r8), pointer :: RootNodulNonstElms_rpvr(:,:,:)       => null()    !root  layer nonstructural element,                  [g d-2]
  real(r8), pointer :: RootNodulStrutElms_rpvr(:,:,:)       => null()    !root layer nodule element,                          [g d-2]
  real(r8), pointer :: CanopyLeafCLyr_pft(:,:)              => null()    !canopy layer leaf C,                                [g d-2]
  real(r8), pointer :: RootMycoNonstElms_pft(:,:,:)         => null()    !nonstructural root-myco chemical element,           [g d-2]
  real(r8), pointer :: RootMyco2ndStrutElms_rpvr(:,:,:,:,:) => null()    !root layer element secondary axes,                  [g d-2]
  real(r8), pointer :: RootMyco1stStrutElms_rpvr(:,:,:,:,:) => null()    !root layer element primary axes,                    [g d-2]
  real(r8), pointer :: RootMyco1stElm_raxs(:,:,:,:)         => null()    !root C primary axes,                                [g d-2]
  real(r8), pointer :: StandDeadKCompElms_pft(:,:,:)        => null()    !standing dead element fraction,                     [g d-2]
  real(r8), pointer :: CanopyNonstElmConc_pft(:,:)          => null()    !canopy nonstructural element concentration,         [g d-2]
  real(r8), pointer :: CanopyNonstElms_pft(:,:)             => null()    !canopy nonstructural element concentration,         [g d-2]
  real(r8), pointer :: CanopyNodulNonstElms_pft(:,:)        => null()    !canopy nodule nonstructural element,                [g d-2]
  real(r8), pointer :: CanopyNoduleNonstCConc_pft(:)          => null()    !nodule nonstructural C,                             [gC d-2]
  real(r8), pointer :: RootMycoActiveBiomC_pvr(:,:,:)       => null()    !root layer structural C,                            [gC d-2]
  real(r8), pointer :: PopuRootMycoC_pvr(:,:,:)             => null()    !root layer C,                                       [gC d-2]
  real(r8), pointer :: RootProteinC_pvr(:,:,:)              => null()    !root layer protein C,                               [gC d-2]
  real(r8), pointer :: RootProteinConc_rpvr(:,:,:)          => null()    !root layer protein C concentration,                 [g g-1]
  real(r8), pointer :: RootMycoNonstElms_rpvr(:,:,:,:)      => null()    !root  layer nonstructural element,                  [g d-2]
  real(r8), pointer :: RootNonstructElmConc_rpvr(:,:,:,:)   => null()    !root  layer nonstructural C concentration,          [g g-1]
  real(r8), pointer :: LeafPetoNonstElmConc_brch(:,:,:)     => null()    !branch nonstructural C concentration,               [g d-2]
  real(r8), pointer :: StructInternodeElms_brch(:,:,:,:)     => null()    !internode C,                                        [g d-2]
  real(r8), pointer :: LeafElmntNode_brch(:,:,:,:)          => null()    !leaf element,                                       [g d-2]
  real(r8), pointer :: LeafProteinC_node(:,:,:)             => null()    !layer leaf protein C,                               [g d-2]
  real(r8), pointer :: PetioleElmntNode_brch(:,:,:,:)       => null()    !sheath chemical element,                            [g d-2]
  real(r8), pointer :: PetoleProteinCNode_brch(:,:,:)       => null()    !layer sheath protein C,                             [g d-2]
  real(r8), pointer :: LeafLayerElms_node(:,:,:,:,:)  => null()    !layer leaf element,                                 [g d-2]
  real(r8), pointer :: tCanLeafC_clyr(:)                    => null()    !total leaf carbon mass in canopy layers,            [gC d-2]
  real(r8), pointer :: StandingDeadInitC_pft(:)             => null()    !initial standing dead C,                            [g C m-2]
  real(r8), pointer :: RootElms_pft(:,:)                    => null()    !plant root element mass,                            [g d-2]
  real(r8), pointer :: RootNoduleElms_pft(:,:)              => null()    !plant root nodule element mass,                     [g d-2]
  real(r8), pointer :: RootNoduleElmsBeg_pft(:,:)           => null()    !previous time step plant root nodule element mass,  [g d-2]  
  real(r8), pointer :: RootElmsBeg_pft(:,:)                 => null()    !plant root element at previous time step,           [g d-2]
  real(r8), pointer :: StandDeadStrutElmsBeg_pft(:,:)       => null()    !standing dead element at previous time step,        [g d-2]
  real(r8), pointer :: RootStrutElms_pft(:,:)               => null()    !plant root structural element mass,                 [g d-2]
  real(r8), pointer :: SeedPlantedElm_pft(:,:)              => null()    !plant stored nonstructural chemical elements at planting,           [gC d-2]
  real(r8), pointer :: SeasonalNonstElms_pft(:,:)           => null()    !plant stored nonstructural element at current step, [g d-2]
  real(r8), pointer :: SeasonalNonstElmsbeg_pft(:,:)        => null()    !plant stored nonstructural element at prev step,    [g d-2]
  real(r8), pointer :: CanopyLeafSheathC_pft(:)              => null()    !canopy leaf + sheath C,                             [g d-2]
  real(r8), pointer :: AvgCanopyBiomC2Graze_pft(:)          => null()    !landscape average canopy shoot C,                   [g d-2]
  real(r8), pointer :: StandDeadStrutElms_pft(:,:)          => null()    !standing dead element,                              [g d-2]
  real(r8), pointer :: CanopyNonstElms_brch(:,:,:)          => null()    !branch nonstructural element,                       [g d-2]
  real(r8), pointer :: C4PhotoShootNonstC_brch(:,:)         => null()    !branch shoot nonstrucal elelment,                   [g d-2]
  real(r8), pointer :: CanopyNodulNonstElms_brch(:,:,:)     => null()    !branch nodule nonstructural element,                [g d-2]
  real(r8), pointer :: CanopyLeafSheathC_brch(:,:)          => null()    !plant branch leaf + sheath C,                       [g d-2]
  real(r8), pointer :: StalkRsrvElms_brch(:,:,:)            => null()    !branch reserve element mass,                        [g d-2]
  real(r8), pointer :: LeafStrutElms_brch(:,:,:)            => null()    !branch leaf structural element mass,                [g d-2]
  real(r8), pointer :: CanopyNodulStrutElms_brch(:,:,:)     => null()   !branch nodule structural element,                    [g d-2]
  real(r8), pointer :: PetoleStrutElms_brch(:,:,:)          => null()   !branch sheath structural element,                    [g d-2]
  real(r8), pointer :: EarStrutElms_brch(:,:,:)             => null()   !branch ear structural chemical element mass,         [g d-2]
  real(r8), pointer :: HuskStrutElms_brch(:,:,:)            => null()   !branch husk structural element mass,                 [g d-2]
  real(r8), pointer :: GrainStrutElms_brch(:,:,:)           => null()   !branch grain structural element mass,                [g d-2]
  real(r8), pointer :: StalkStrutElms_brch(:,:,:)           => null()   !branch stalk structural element mass,                [g d-2]
  real(r8), pointer :: ShootElms_brch(:,:,:)                => null()   !branch shoot structural element mass,                [g d-2]
  real(r8), pointer :: LeafChemElmRemob_brch(:,:,:)         => null()   !branch leaf structural element,                      [g d-2]
  real(r8), pointer :: PetioleChemElmRemob_brch(:,:,:)      => null()   !branch sheath structural element,                    [g d-2]
  real(r8), pointer :: SenecStalkStrutElms_brch(:,:,:)      => null()   !branch stalk structural element,                     [g d-2]
  real(r8), pointer :: SapwoodBiomassC_brch(:,:)            => null()   !branch live stalk C,                                 [gC d-2]
  real(r8), pointer :: StalkStrutElms_pft(:,:)              => null()   !canopy stalk structural element mass,                [g d-2]
  real(r8), pointer :: ShootElms_pft(:,:)                   => null()   !current time whole plant shoot element mass,         [g d-2]
  real(r8), pointer :: ShootElmsBeg_pft(:,:)                => null()   !previous whole plant shoot element mass,             [g d-2]
  real(r8), pointer :: CanopySapwoodC_pft(:)                  => null()   !canopy active stalk C,                               [g d-2]
  real(r8), pointer :: LeafStrutElms_pft(:,:)               => null()   !canopy leaf structural element mass,                 [g d-2]
  real(r8), pointer :: PetoleStrutElms_pft(:,:)             => null()   !canopy sheath structural element mass,               [g d-2]
  real(r8), pointer :: StalkRsrvElms_pft(:,:)               => null()   !canopy reserve element mass,                         [g d-2]
  real(r8), pointer :: HuskStrutElms_pft(:,:)               => null()   !canopy husk structural element mass,                 [g d-2]
  real(r8), pointer :: RootBiomCPerPlant_pft(:)             => null()   !root C biomass per plant,                            [g p-1]
  real(r8), pointer :: GrainStrutElms_pft(:,:)              => null()   !canopy grain structural element,                     [g d-2]
  real(r8), pointer :: EarStrutElms_pft(:,:)                => null()   !canopy ear structural element,                       [g d-2]
  real(r8), pointer :: CanopyMassC_pft(:)                   => null()   !Canopy biomass C,                                    [g d-2]
  real(r8), pointer :: ROOTNLim_rpvr(:,:,:)                 => null()   !root N-limitation, 0->1 weaker limitation, [-]     
  real(r8), pointer :: ROOTPLim_rpvr(:,:,:)                 => null()   !root P-limitation, 0->1 weaker limitation, [-]         
  real(r8), pointer :: LeafC3ChlC_brch(:,:)                 => null()   !Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]     
  real(r8), pointer :: LeafC4ChlC_brch(:,:)                 => null()   !Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]     
  real(r8), pointer :: LeafRubiscoC_brch(:,:)               => null()   !Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
  real(r8), pointer :: LeafPEPC_brch(:,:)                   => null()   !PEP C for the branches, [gC PEP d-2]
  real(r8), pointer :: LeafC3ChlCperm2LA_pft(:)             => null()   !Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]     
  real(r8), pointer :: LeafC4ChlCperm2LA_pft(:)             => null()   !Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]     
  real(r8), pointer :: LeafRubiscoCperm2LA_pft(:)           => null()   !Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
  real(r8), pointer :: LeafPEPCperm2LA_pft(:)               => null()   !PEP C for the branches, [gC PEP d-2]
  real(r8), pointer :: SpecificLeafArea_pft(:)              => null()   !specifc leaf area per g C of leaf mass, [m2 leaf area (gC leaf C)-1]  
  real(r8), pointer :: ShootNoduleElms_pft(:,:)             => null()   !current time canopy nodule element mass, [g d-2]      
  real(r8), pointer :: ShootNoduleElmsBeg_pft(:,:)          => null()   !previous time canopy nodule element mass, [g d-2]        
  contains
    procedure, public :: Init => plt_biom_init
    procedure, public :: Destroy => plt_biom_destroy
  end type plant_biom_type

  type, public :: plant_ew_type
  real(r8) :: SnowDepth                     !snowpack depth, [m]
  real(r8) :: VcumWatSnow_col               !water volume in snowpack, [m3 d-2]
  real(r8) :: VcumDrySnoWE_col              !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8) :: Air_Heat_Sens_store_col       !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VLHeatCapSnowMin_col          !minimum snowpack heat capacity, [MJ d-2 K-1]
  real(r8) :: VLHeatCapSurfSnow_col         !snowpack heat capacity, [MJ m-3 K-1]
  real(r8) :: CanopyHeatStor_col            !total canopy heat content, [MJ  d-2]
  real(r8) :: LWRadCanG                     !grid total canopy LW emission, [MJ d-2 h-1]
  real(r8) :: QVegET_col                    !total canopy evaporation + transpiration, [m3 d-2]
  real(r8) :: VapXAir2Canopy_col            !grid canopy evaporation, [m3 d-2]
  real(r8) :: HeatFlx2Canopy_col            !total canopy heat flux, [MJ  d-2]
  real(r8) :: H2OLoss_CumYr_col             !total subsurface water flux, [m3 d-2]
  real(r8) :: VPA                           !vapor concentration, [m3 m-3]
  real(r8) :: TairK                         !air temperature, [K]
  real(r8) :: CanopyWat_col                 !total canopy water content stored with dry matter, [m3 d-2]
  real(r8) :: Eco_Heat_Latent_col           !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8) :: Air_Heat_Latent_store_col     !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VcumIceSnow_col               !ice volume in snowpack, [m3 d-2]
  real(r8) :: TKSnow                        !snow temperature, [K]
  real(r8) :: WatHeldOnCanopy_col           !canopy surface water content, [m3 d-2]
  real(r8) :: Eco_Heat_Sens_col             !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8) :: RoughHeight                   !canopy surface roughness height, [m]
  real(r8) :: ZERO4PlantDisplace_col        !zero plane displacement height, [m]
  real(r8) :: AbvCanopyBndlResist_col       !isothermal boundary layer resistance, [h m-1]
  real(r8) :: RIB                           !Richardson number for calculating boundary layer resistance, [-]
  real(r8) :: Eco_Heat_GrndSurf_col         !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8) :: HeatCanopy2Dist_col           !canopy energy +/- due to disturbance, [MJ /d2]
  real(r8) :: QCanopyWat2Dist_col           !canopy water +/- due to disturbance, [m3 H2O/d2]
  real(r8) :: TPlantRootH2OUptake_col       !total water uptake by roots, [m3 H2O/d2]

  real(r8), pointer :: Transpiration_pft(:)           => null()    !canopy transpiration,                                         [m3 d-2 h-1]
  real(r8), pointer :: PSICanopyTurg_pft(:)           => null()    !plant canopy turgor water potential,                          [MPa]
  real(r8), pointer :: PrecIntcptByCanopy_pft(:)      => null()    !water flux into canopy,                                       [m3 d-2 h-1]
  real(r8), pointer :: PSICanopy_pft(:)               => null()    !canopy total water potential,                                 [Mpa]
  real(r8), pointer :: VapXAir2Canopy_pft(:)          => null()    !canopy evaporation,                                           [m3 d-2 h-1]
  real(r8), pointer :: HeatStorCanopy_pft(:)          => null()    !canopy storage heat flux,                                     [MJ d-2 h-1]
  real(r8), pointer :: EvapTransLHeat_pft(:)          => null()    !canopy latent heat flux,                                      [MJ d-2 h-1]
  real(r8), pointer :: CanopyIsothBndlResist_pft(:)   => null()    !canopy isothermal boundary later resistance,                  [h m-1]
  real(r8), pointer :: TKS_vr(:)                      => null()    !mean annual soil temperature,                                 [K]
  real(r8), pointer :: PSICanPDailyMin_pft(:)         => null()    !minimum daily canopy water potential,                         [MPa]
  real(r8), pointer :: TdegCCanopy_pft(:)             => null()    !canopy temperature,                                           [oC]
  real(r8), pointer :: DeltaTKC_pft(:)                => null()    !change in canopy temperature,                                 [K]
  real(r8), pointer :: ENGYX_pft(:)                   => null()    !canopy heat storage from previous time step,                  [MJ d-2]
  real(r8), pointer :: TKC_pft(:)                     => null()    !canopy temperature,                                           [K]
  real(r8), pointer :: PSICanopyOsmo_pft(:)           => null()    !canopy osmotic water potential,                               [Mpa]
  real(r8), pointer :: CanOsmoPsi0pt_pft(:)           => null()    !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8), pointer :: HeatXAir2PCan_pft(:)           => null()    !canopy sensible heat flux,                                    [MJ d-2 h-1]
  real(r8), pointer :: TKCanopy_pft(:)                => null()    !canopy temperature,                                           [K]
  real(r8), pointer :: PSIRoot_pvr(:,:,:)             => null()    !root total water potential,                                   [Mpa]
  real(r8), pointer :: PSIRootOSMO_vr(:,:,:)          => null()    !root osmotic water potential,                                 [Mpa]
  real(r8), pointer :: PSIRootTurg_vr(:,:,:)          => null()    !root turgor water potential,                                  [Mpa]
  real(r8), pointer :: RPlantRootH2OUptk_pvr(:,:,:)   => null()    !root water uptake,                                            [m3 d-2 h-1]
  real(r8), pointer :: RootH2OUptkStress_pvr(:,:,:)   => null()    !root water uptake stress indicated by rate,                   [m3 d-2 h-1] 
  real(r8), pointer :: THeatLossRoot2Soil_vr(:)       => null()    !total root heat uptake,                                       [MJ d-2]
  real(r8), pointer :: TWaterPlantRoot2Soil_vr(:)     => null()    !total root water uptake,                                      [m3 d-2]
  real(r8), pointer :: WatHeldOnCanopy_pft(:)         => null()    !canopy surface water content,                                 [m3 d-2]
  real(r8), pointer :: CanopyBiomWater_pft(:)         => null()    !canopy water content,                                         [m3 d-2]
  real(r8), pointer :: VHeatCapCanopy_pft(:)          => null()    !canopy heat capacity,                                         [MJ d-2 K-1]
  real(r8), pointer :: ElvAdjstedSoilH2OPSIMPa_vr(:)  => null()    !soil micropore total water potential,                         [MPa]
  real(r8), pointer :: ETCanopy_CumYr_pft(:)          => null()    !total transpiration,                                          [m H2O d-2]
  real(r8), pointer :: QdewCanopy_pft(:)              => null()    !dew fall on to canopy,                                        [m3 H2O d-2 h-1]
  real(r8), pointer :: RootResist4H2O_pvr(:,:,:)      => null()    !total root (axial+radial) resistance for water uptake,        [MPa-1 h-1]     
  real(r8), pointer :: RootRadialKond2H2O_pvr(:,:,:)  => null()    !radial root conductance for water uptake, [m3 H2O h-1 MPa-1]
  real(r8), pointer :: RootAxialKond2H2O_pvr(:,:,:)   => null()    !axial root conductance for water uptake, [m3 H2O h-1 MPa-1]

  contains
    procedure, public :: Init => plt_ew_init
    procedure, public :: Destroy=> plt_ew_destroy
  end type plant_ew_type

  type, public :: plant_disturb_type

  real(r8) :: XCORP                     !factor for surface litter incorporation and soil mixing,[-]
  real(r8) :: DCORP                     !soil mixing fraction with tillage, [-]
  integer  :: IYTYP                     !fertilizer release type from fertilizer input file,[-]
  integer  :: iSoilDisturbType_col      !soil disturbance type, [-]

  integer,  pointer :: iHarvestDay_pft(:)       => null()  !day of harvest, [-]
  integer,  pointer :: iYearPlanting_pft(:)     => null()  !year of planting,[-]
  integer,  pointer :: iPlantingDay_pft(:)      => null()  !day of planting,[-]
  integer,  pointer :: iPlantingYear_pft(:)     => null()  !year of planting,[-]
  integer,  pointer :: iHarvestYear_pft(:)      => null()  !year of harvest,[-]
  integer,  pointer :: iDayPlanting_pft(:)      => null()  !day of planting,[-]
  integer,  pointer :: iDayPlantHarvest_pft(:)  => null()  !day of harvest,[-]
  integer,  pointer :: iYearPlantHarvest_pft(:) => null()  !year of harvest,[-]
  integer,  pointer :: iHarvstType_pft(:)       => null()  !type of harvest,[-]
  integer,  pointer :: jHarvstType_pft(:)           => null()  !flag for stand replacing disturbance,[-]

  real(r8), pointer :: EcoHavstElmnt_CumYr_col(:)   => null()  !ecosystem harvest element,                                    [gC d-2]
  real(r8), pointer :: O2ByFire_CumYr_pft(:)        => null()  !plant O2 uptake from fire,                                    [g d-2 ]
  real(r8), pointer :: FERT(:)                      => null()  !fertilizer application,                                       [g m-2]
  real(r8), pointer :: FracBiomHarvsted(:,:,:)      => null()  !harvest efficiency,                                           [-]
  real(r8), pointer :: CanopyHeightCut_pft(:)   => null()  !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  real(r8), pointer :: THIN_pft(:)                  => null()  !thinning of plant population,                                 [-]
  real(r8), pointer :: CH4ByFire_CumYr_pft(:)       => null()  !plant CH4 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: CO2ByFire_CumYr_pft(:)       => null()  !plant CO2 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: N2ObyFire_CumYr_pft(:)       => null()  !plant N2O emission from fire,                                 [g d-2 ]
  real(r8), pointer :: NH3byFire_CumYr_pft(:)       => null()  !plant NH3 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: PO4byFire_CumYr_pft(:)       => null()  !plant PO4 emission from fire,                                 [g d-2 ]
  real(r8), pointer :: EcoHavstElmntCum_pft(:,:)    => null()  !total plant element harvest,                                  [gC d-2 ]
  real(r8), pointer :: EcoHavstElmnt_CumYr_pft(:,:) => null()  !plant element harvest,                                        [g d-2 ]
  real(r8), pointer :: PlantElmDistLoss_pft(:,:)    => null()  !plant element loss to disturbance,                            [g d-2 h-1]
  contains
    procedure, public :: Init    =>  plt_disturb_init
    procedure, public :: Destroy => plt_disturb_destroy
  end type plant_disturb_type

  type, public :: plant_bgcrate_type
  real(r8) :: Eco_NBP_CumYr_col        !total NBP, [g d-2]
  real(r8) :: NetCO2Flx2Canopy_col     !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8) :: ECO_ER_col               !ecosystem respiration, [g d-2 h-1]
  real(r8) :: Eco_AutoR_CumYr_col      !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8) :: TRootH2Flx_col           !total root H2 flux, [g d-2]
  real(r8) :: Canopy_NEE_col           !total net CO2 fixation, [gC d-2]
  real(r8), pointer :: Nutruptk_fClim_rpvr(:,:,:)          => null()  !Carbon limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(r8), pointer :: Nutruptk_fNlim_rpvr(:,:,:)          => null()  !Nitrogen limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(r8), pointer :: Nutruptk_fPlim_rpvr(:,:,:)          => null()  !Phosphorus limitation for root nutrient uptake,(0->1),stronger limitation, [-]
  real(r8), pointer :: Nutruptk_fProtC_rpvr(:,:,:)         => null()  !transporter scalar indicated by protein for root nutrient uptake, greater value greater capacity, [-] 
  real(r8), pointer :: LitrFallStrutElms_col(:)            => null()  !total LitrFall structural element mass,     [g d-2 h-1]
  real(r8), pointer :: NetPrimProduct_pft(:)               => null()  !total net primary productivity,             [gC d-2]
  real(r8), pointer :: NH3Dep2Can_pft(:)                   => null()  !canopy NH3 flux,                            [g d-2 h-1]
  real(r8), pointer :: tRootMycoExud2Soil_vr(:,:,:)        => null()  !total root element exchange,                 [g d-2 h-1]
  real(r8), pointer :: RCO2Nodule_pvr(:,:)                 => null()  !layered root nodule respiration, [gC d-2 h-1]
  real(r8), pointer :: RootN2Fix_pvr(:,:)                  => null()  !root N2 fixation,                            [gN d-2 h-1]
  real(r8), pointer :: RootO2_TotSink_pvr(:,:,:)           => null()  !root O2 sink for autotrophic respiraiton,     [gC d-2 h-1]
  real(r8), pointer :: RootO2_TotSink_vr(:)                => null()  !all root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
  real(r8), pointer :: RCanMaintDef_CO2_pft(:)             => null()  !canopy maintenance respiraiton deficit as CO2, [gC d-2 h-1]  
  real(r8), pointer :: CanopyRespC_CumYr_pft(:)            => null()  !total autotrophic respiration,               [gC d-2 ]
  real(r8), pointer :: CanopyNLimFactor_brch(:,:)          => null()  !Canopy N-limitation factor, [0->1] weaker limitation,[-]  
  real(r8), pointer :: CanopyPLimFactor_brch(:,:)          => null()  !Canopy P-limitation factor, [0->1] weaker limitation,[-]    
  real(r8), pointer :: LitrfallElms_pvr(:,:,:,:,:)         => null()  !plant LitrFall element,                      [g d-2 h-1]
  real(r8), pointer :: LitrFallElms_brch(:,:,:)            => null()  !litterfall from the branch, [g d-2 h-1]  
  real(r8), pointer :: RootMaintDef_CO2_pvr(:,:,:)         => null()  !plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]  
  real(r8), pointer :: REcoO2DmndResp_vr(:)                => null()  !total root + microbial O2 uptake,            [g d-2 h-1]
  real(r8), pointer :: REcoNH4DmndBand_vr(:)               => null()   !total root + microbial NH4 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: REcoH1PO4DmndSoil_vr(:)             => null()   !HPO4 demand in non-band by all microbial,   root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: REcoH1PO4DmndBand_vr(:)             => null()   !HPO4 demand in band by all microbial,       root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: REcoNO3DmndSoil_vr(:)               => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: REcoNH4DmndSoil_vr(:)               => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: REcoNO3DmndBand_vr(:)               => null()   !total root + microbial NO3 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: REcoH2PO4DmndSoil_vr(:)             => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: REcoH2PO4DmndBand_vr(:)             => null()   !total root + microbial PO4 uptake band,     [gP d-2 h-1]
  real(r8), pointer :: RH2PO4EcoDmndSoilPrev_vr(:)         => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: RH2PO4EcoDmndBandPrev_vr(:)         => null()   !total root + microbial PO4 uptake band,     [gP d-2 h-1]
  real(r8), pointer :: RH1PO4EcoDmndSoilPrev_vr(:)         => null()   !HPO4 demand in non-band by all microbial,   root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: RH1PO4EcoDmndBandPrev_vr(:)         => null()   !HPO4 demand in band by all microbial,       root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: RNO3EcoDmndSoilPrev_vr(:)           => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4EcoDmndSoilPrev_vr(:)           => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4EcoDmndBandPrev_vr(:)           => null()   !total root + microbial NH4 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: RNO3EcoDmndBandPrev_vr(:)           => null()   !total root + microbial NO3 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: RGasTranspFlxPrev_vr(:,:)           => null()   !net gaseous flux,                           [g d-2 h-1]
  real(r8), pointer :: RO2AquaSourcePrev_vr(:)             => null()   !net aqueous O2 flux,                        [g d-2 h-1]
  real(r8), pointer :: RO2EcoDmndPrev_vr(:)                => null()   !total root + microbial O2 uptake,           [g d-2 h-1]
  real(r8), pointer :: RootCO2Emis2Root_vr(:)              => null()   !total root CO2 flux,                        [gC d-2 h-1]
  real(r8), pointer :: RUptkRootO2_vr(:)                   => null()   !total root internal O2 flux,                [g d-2 h-1]
  real(r8), pointer :: LitrfalStrutElms_vr(:,:,:,:)        => null()   !total LitrFall element,                       [g d-2 h-1]
  real(r8), pointer :: REcoDOMProd_vr(:,:,:)               => null()   !net microbial DOC flux,                      [gC d-2 h-1]
  real(r8), pointer :: CO2NetFix_pft(:)                    => null()   !canopy net CO2 exchange,                     [gC d-2 h-1]
  real(r8), pointer :: GrossCO2Fix_pft(:)                  => null()   !total gross CO2 fixation,                    [gC d-2 ]
  real(r8), pointer :: LitrfallElms_pft(:,:)           => null()   !plant element Litrfall,                      [g d-2 h-1]
  real(r8), pointer :: LitrfallAbvgElms_pft(:,:)           => null()   !aboveground litterfall, [g d-2 h-1]
  real(r8), pointer :: LitrfallBlgrElms_pft(:,:)           => null()   !belwground litterfall, [g d-2 h-1]
  real(r8), pointer :: RootGasLossDisturb_pft(:,:)         => null()   !gaseous flux fron root disturbance,           [g d-2 h-1]
  real(r8), pointer :: SurfLitrfalStrutElms_CumYr_pft(:,:) => null()   !total surface LitrFall element,              [g d-2]
  real(r8), pointer :: LitrfalStrutElms_CumYr_pft(:,:)     => null()   !total plant element LitrFall,                [g d-2 ]
  real(r8), pointer :: GrossResp_pft(:)                    => null()   !total plant respiration,                     [gC d-2 ]
  real(r8), pointer :: CanopyGrosRCO2_pft(:)               => null()   !canopy plant+nodule autotrophic respiraiton, [gC d-2]
  real(r8), pointer :: CanopyResp_brch(:,:)                => null()   !canopy respiration for a branch, [gC d-2 h-1]      
  real(r8), pointer :: RootAutoCO2_pft(:)                  => null()   !root autotrophic respiraiton, [gC d-2]
  real(r8), pointer :: NodulInfectElms_pft(:,:)            => null()   !nodule infection chemical element mass, [g d-2]
  real(r8), pointer :: NH3Emis_CumYr_pft(:)                => null()   !total canopy NH3 flux,                       [gN d-2 ]
  real(r8), pointer :: PlantN2Fix_CumYr_pft(:)             => null()   !total plant N2 fixation,                     [g d-2 ]

  contains
    procedure, public :: Init  => plt_bgcrate_init
    procedure, public :: Destroy  => plt_bgcrate_destroy
  end type plant_bgcrate_type

  type, public :: plant_rootbgc_type
  real(r8), pointer :: canopy_growth_pft(:)              => null()  !canopy structural C growth rate,                                [gC d-2 h-1]
  real(r8), pointer :: TRootGasLossDisturb_col(:)        => null()  !total root gas content,                                         [g d-2]
  real(r8), pointer :: trcs_Soil2plant_uptake_vr(:,:)    => null()  !total root-soil solute flux non-band,                           [g d-2 h-1]
  real(r8), pointer :: RootMycoExudElms_pft(:,:)         => null()  !total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
  real(r8), pointer :: CanopyN2Fix_pft(:)                => null()  !total canopy N2 fixation, [g d-2 h-1]
  real(r8), pointer :: RootN2Fix_pft(:)                  => null()  !total root N2 fixation,                                         [g d-2 h-1]
  real(r8), pointer :: RootNO3Uptake_pft(:)              => null()  !total root uptake of NO3,                                       [g d-2 h-1]
  real(r8), pointer :: RootNH4Uptake_pft(:)              => null()  !total root uptake of NH4,                                       [g d-2 h-1]
  real(r8), pointer :: RootHPO4Uptake_pft(:)             => null()  !total root uptake of HPO4,                                      [g d-2 h-1]
  real(r8), pointer :: RootH2PO4Uptake_pft(:)            => null()  !total root uptake of PO4,                                       [g d-2 h-1]
  real(r8), pointer :: RootMycoExudEUptk_pvr(:,:,:,:,:)  => null()  !root uptake (+ve) - exudation (-ve) of DOE,                     [g d-2 h-1]
  real(r8), pointer :: PlantRootSoilElmNetX_pft(:,:)     => null()  !net root element uptake (+ve) - exudation (-ve),                [gC d-2 h-1]
  real(r8), pointer :: REcoUptkSoilO2M_vr(:,:)              => null()  !total O2 sink,                                                  [g d-2 t-1]
  real(r8), pointer :: ZERO4Uptk_pft(:)                  => null()  !threshold zero for uptake calculation,                          [-]
  real(r8), pointer :: CMinPO4Root_pft(:,:)              => null()  !minimum PO4 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxPO4Root_pft(:,:)              => null()  !maximum root PO4 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmPO4Root_pft(:,:)                => null()  !Km for root PO4 uptake,                                         [g m-3]
  real(r8), pointer :: CminNO3Root_pft(:,:)              => null()  !minimum NO3 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxNO3Root_pft(:,:)              => null()  !maximum root NO3 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmNO3Root_pft(:,:)                => null()  !Km for root NO3 uptake,                                         [g m-3]
  real(r8), pointer :: CMinNH4Root_pft(:,:)              => null()  !minimum NH4 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxNH4Root_pft(:,:)              => null()  !maximum root NH4 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmNH4Root_pft(:,:)                => null()  !Km for root NH4 uptake,                                         [g m-3]
  real(r8), pointer :: RCO2Emis2Root_pvr(:,:,:)          => null()  !aqueous CO2 flux from roots to root water,                      [g d-2 h-1]
  real(r8), pointer :: RootO2Uptk_pvr(:,:,:)             => null()  !aqueous O2 flux from roots to root water,                       [g d-2 h-1]
  real(r8), pointer :: RootUptkSoiSol_pvr(:,:,:,:)       => null()  !aqueous CO2 flux from roots to soil water,                      [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx_pvr(:,:,:,:)    => null()  !gaseous tracer flux through roots,                              [g d-2 h-1]
  real(r8), pointer :: trcg_Root_gas2aqu_flx_vr(:,:,:,:) => null()  !dissolution (+ve) - volatilization (-ve) gas flux in roots,     [g d-2 h-1]
  real(r8), pointer :: RootO2Dmnd4Resp_pvr(:,:,:)        => null()  !root  O2 demand from respiration,                               [g d-2 h-1]

  real(r8), pointer :: RootNH4DmndSoilPrev_pvr(:,:,:)    => null()  !previous root uptake of NH4 non-band unconstrained by NH4,      [g d-2 h-1]
  real(r8), pointer :: RootNH4DmndBandPrev_pvr(:,:,:)    => null()  !previous root uptake of NO3 band unconstrained by NO3,          [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndSoilPrev_pvr(:,:,:)    => null()  !previous root uptake of NH4 band unconstrained by NH4,          [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndBandPrev_pvr(:,:,:)    => null()  !previous root uptake of NO3 non-band unconstrained by NO3,      [g d-2 h-1]
  real(r8), pointer :: VmaxNH4Root_pvr(:,:,:)            => null()  !maximum NH4 uptake rate,                                        [gN h-1 (gC root)-1]
  real(r8), pointer :: VmaxNO3Root_pvr(:,:,:)            => null()  !maximum NO3 uptake rate,                                        [gN h-1 (gC root)-1]
  real(r8), pointer :: RootNH4DmndSoil_pvr(:,:,:)        => null()  !root uptake of NH4 non-band unconstrained by NH4,               [g d-2 h-1]
  real(r8), pointer :: RootNH4DmndBand_pvr(:,:,:)        => null()  !root uptake of NO3 band unconstrained by NO3,                   [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndSoil_pvr(:,:,:)        => null()  !root uptake of NH4 band unconstrained by NH4,                   [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndBand_pvr(:,:,:)        => null()  !root uptake of NO3 non-band unconstrained by NO3,               [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndSoil_pvr(:,:,:)      => null()  !root uptake of H2PO4 non-band,                                  [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndBand_pvr(:,:,:)      => null()  !root uptake of H2PO4 band,                                      [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndSoil_pvr(:,:,:)      => null()  !HPO4 demand in non-band by each root population,                [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndBand_pvr(:,:,:)      => null()  !HPO4 demand in band by each root population,                    [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndSoilPrev_pvr(:,:,:)  => null()  !previous step root uptake of H2PO4 non-band,                    [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndBandPrev_pvr(:,:,:)  => null()  !previous step root uptake of H2PO4 band,                        [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndSoilPrev_pvr(:,:,:)  => null()  !previous step HPO4 demand in non-band by each root population,  [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndBandPrev_pvr(:,:,:)  => null()  !previous step HPO4 demand in band by each root population,      [g d-2 h-1]

  real(r8), pointer :: RAutoRootO2Limter_rpvr(:,:,:)     => null()  !O2 constraint to root respiration (0-1),                        [-]
  real(r8), pointer :: trcg_rootml_pvr(:,:,:,:)          => null() !root gas content,                                                [g d-2]
  real(r8), pointer :: trcs_rootml_pvr(:,:,:,:)          => null() !root aqueous content,                                            [g d-2]
  real(r8), pointer :: RootGasConductance_rpvr(:,:,:,:)   => null()  !Conductance for gas diffusion                                   [m3 d-2 h-1]
  real(r8), pointer :: NH3Dep2Can_brch(:,:)              => null()  !gaseous NH3 flux fron root disturbance band,                    [g d-2 h-1]
  real(r8), pointer :: GPP_brch(:,:)                    => null()  !dGPP (C4-C3 product) over branch, [gC d-2 h-1]
  real(r8), pointer :: RootNutUptake_pvr(:,:,:,:)        => null()  !root uptake of Nutrient band,                                   [g d-2 h-1]
  real(r8), pointer :: RootOUlmNutUptake_pvr(:,:,:,:)    => null()  !root uptake of NH4 band unconstrained by O2,                    [g d-2 h-1]
  real(r8), pointer :: RootCUlmNutUptake_pvr(:,:,:,:)    => null()  !root uptake of NH4 band unconstrained by root nonstructural C,  [g d-2 h-1]
  real(r8), pointer :: RootRespPotent_pvr(:,:,:)         => null()  !root respiration unconstrained by O2,                           [g d-2 h-1]
  real(r8), pointer :: RootCO2EmisPot_pvr(:,:,:)         => null()  !root CO2 efflux unconstrained by root nonstructural C,          [g d-2 h-1]
  real(r8), pointer :: RootCO2Autor_pvr(:,:,:)           => null()  !root respiration constrained by O2,                             [g d-2 h-1]
  real(r8), pointer :: RootCO2AutorX_pvr(:,:,:)          => null()  !root respiration from previous time step,                       [g d-2 h-1]
  real(r8), pointer :: PlantExudElm_CumYr_pft(:,:)       => null()  !total net root element uptake (+ve) - exudation (-ve),          [gC d-2 ]
  real(r8), pointer :: trcg_root_vr(:,:)                 => null()   !total root internal gas flux,                                  [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx_vr(:,:)         => null()   !total internal root gas flux,                                  [gC d-2 h-1]
  real(r8), pointer :: CO2FixCL_pft(:)                   => null()   !Rubisco-limited CO2 fixation,                                  [gC d-2 h-1]
  real(r8), pointer :: RootNutUptakeN_pft(:)             => null()   !total N uptake by plant roots,                                 [gN d-h2 h-1]
  real(r8), pointer :: RootNutUptakeP_pft(:)             => null()   !total P uptake by plant roots,                                 [gP d-h2 h-1]
  real(r8), pointer :: CO2FixLL_pft(:)                   => null()   !Light-limited CO2 fixation,                                    [gC d-h2 h-1]
  real(r8), pointer :: RootUptk_N_CumYr_pft(:)           => null()  !cumulative plant N uptake,                                      [gN d-2]
  real(r8), pointer :: RootUptk_P_CumYr_pft(:)           => null()  !cumulative plant P uptake,                                      [gP d-2]
  real(r8), pointer :: RootCO2Ar2Soil_pvr(:,:)           => null()  !root respiration released to soil,                              [gC d-2 h-1]
  real(r8), pointer :: RootCO2Ar2RootX_pvr(:,:)          => null()  !root respiration released to root,                              [gC d-2 h-1]
  real(r8), pointer :: trcs_deadroot2soil_pvr(:,:,:)     => null()  !gases released to soil upong dying roots,                       [g d-2 h-1]
  contains
    procedure, public :: Init => plt_rootbgc_init
    procedure, public :: Destroy  => plt_rootbgc_destroy
  end type plant_rootbgc_type

  type(plant_siteinfo_type) , public, target :: plt_site      !site info
  type(plant_rootbgc_type)  , public, target :: plt_rbgc      !root bgc
  type(plant_bgcrate_type)  , public, target :: plt_bgcr      !bgc reaction
  type(plant_disturb_type)  , public, target :: plt_distb     !plant disturbance type
  type(plant_ew_type)       , public, target :: plt_ew        !plant energy and water type
  type(plant_allometry_type), public, target :: plt_allom     !plant allometric parameters
  type(plant_biom_type)     , public, target :: plt_biom      !plant biomass variables
  type(plant_soilchem_type) , public, target :: plt_soilchem  !soil bgc interface with plant root
  type(plant_pheno_type)    , public, target :: plt_pheno     !plant phenology
  type(plant_morph_type)    , public, target :: plt_morph     !plant morphology
  type(plant_radiation_type), public, target :: plt_rad       !plant radiation type
  type(plant_photosyns_type), public, target :: plt_photo     !plant photosynthesis type

  contains

  subroutine plt_rootbgc_init(this)

  implicit none
  class(plant_rootbgc_type) :: this
  allocate(this%trcs_Soil2plant_uptake_vr(ids_beg:ids_end,JZ1)); this%trcs_Soil2plant_uptake_vr=0._r8
  allocate(this%trcg_rootml_pvr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcg_rootml_pvr=spval
  allocate(this%trcs_rootml_pvr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcs_rootml_pvr=spval
  allocate(this%RootGasConductance_rpvr(idg_beg:idg_NH3,jroots,JZ1,JP1)); this%RootGasConductance_rpvr=0._r8
  allocate(this%TRootGasLossDisturb_col(idg_beg:idg_NH3));this%TRootGasLossDisturb_col=spval
  allocate(this%REcoUptkSoilO2M_vr(60,0:JZ1)); this%REcoUptkSoilO2M_vr=spval
  allocate(this%RootMycoExudEUptk_pvr(NumPlantChemElms,jroots,1:jcplx,0:JZ1,JP1));this%RootMycoExudEUptk_pvr=spval
  allocate(this%PlantRootSoilElmNetX_pft(NumPlantChemElms,JP1)); this%PlantRootSoilElmNetX_pft=spval
  allocate(this%PlantExudElm_CumYr_pft(NumPlantChemElms,JP1));this%PlantExudElm_CumYr_pft=spval
  allocate(this%RootUptk_N_CumYr_pft(JP1)); this%RootUptk_N_CumYr_pft=spval
  allocate(this%RootUptk_P_CumYr_pft(JP1)); this%RootUptk_P_CumYr_pft=spval  
  allocate(this%RootMycoExudElms_pft(NumPlantChemElms,JP1));this%RootMycoExudElms_pft=spval
  allocate(this%RootN2Fix_pft(JP1)); this%RootN2Fix_pft=spval
  allocate(this%CanopyN2Fix_pft(JP1));this%CanopyN2Fix_pft=spval
  allocate(this%canopy_growth_pft(JP1)); this%canopy_growth_pft=spval
  allocate(this%RootNO3Uptake_pft(JP1)); this%RootNO3Uptake_pft=0._r8
  allocate(this%RootNH4Uptake_pft(JP1)); this%RootNH4Uptake_pft=0._r8
  allocate(this%RootHPO4Uptake_pft(JP1)); this%RootHPO4Uptake_pft=0._r8
  allocate(this%RootH2PO4Uptake_pft(JP1)); this%RootH2PO4Uptake_pft=0._r8

  allocate(this%ZERO4Uptk_pft(JP1)); this%ZERO4Uptk_pft=spval
  allocate(this%RootRespPotent_pvr(jroots,JZ1,JP1)); this%RootRespPotent_pvr=spval
  allocate(this%RootCO2EmisPot_pvr(jroots,JZ1,JP1)); this%RootCO2EmisPot_pvr=spval
  allocate(this%RootCO2Autor_pvr(jroots,JZ1,JP1)); this%RootCO2Autor_pvr=0._r8
  allocate(this%trcs_deadroot2soil_pvr(idg_beg:idg_NH3,JZ1,JP1));this%trcs_deadroot2soil_pvr=0._r8
  allocate(this%RootCO2Ar2Soil_pvr(JZ1,JP1)); this%RootCO2Ar2Soil_pvr=0._r8
  allocate(this%RootCO2Ar2RootX_pvr(JZ1,JP1)); this%RootCO2Ar2RootX_pvr=0._r8
  allocate(this%RootCO2AutorX_pvr(jroots,JZ1,JP1)); this%RootCO2AutorX_pvr=spval
  allocate(this%RootNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1)); this%RootNutUptake_pvr=0._r8
  allocate(this%RootOUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1));this%RootOUlmNutUptake_pvr=spval
  allocate(this%RootCUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1));this%RootCUlmNutUptake_pvr=spval
  allocate(this%NH3Dep2Can_brch(MaxNumBranches,JP1));this%NH3Dep2Can_brch=spval
  allocate(this%GPP_brch(MaxNumBranches,JP1)); this%GPP_brch=spval
  allocate(this%CO2FixCL_pft(JP1)); this%CO2FixCL_pft=spval
  allocate(this%CO2FixLL_pft(JP1)); this%CO2FixLL_pft=spval
  allocate(this%RootNutUptakeN_pft(JP1));this%RootNutUptakeN_pft=spval
  allocate(this%RootNutUptakeP_pft(JP1));this%RootNutUptakeP_pft=spval
  allocate(this%trcg_air2root_flx_vr(idg_beg:idg_NH3,JZ1));this%trcg_air2root_flx_vr=spval
  allocate(this%trcg_root_vr(idg_beg:idg_NH3,JZ1));this%trcg_root_vr=spval

  allocate(this%trcg_air2root_flx_pvr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcg_air2root_flx_pvr=spval
  allocate(this%trcg_Root_gas2aqu_flx_vr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcg_Root_gas2aqu_flx_vr=spval
  allocate(this%RootO2Dmnd4Resp_pvr(jroots,JZ1,JP1));this%RootO2Dmnd4Resp_pvr=spval
  allocate(this%RootNH4DmndSoil_pvr(jroots,JZ1,JP1));this%RootNH4DmndSoil_pvr=spval
  allocate(this%VmaxNH4Root_pvr(jroots,JZ1,JP1)); this%VmaxNH4Root_pvr=0._r8
  allocate(this%VmaxNO3Root_pvr(jroots,JZ1,JP1)); this%VmaxNO3Root_pvr=0._r8
  allocate(this%RootNH4DmndBand_pvr(jroots,JZ1,JP1));this%RootNH4DmndBand_pvr=spval
  allocate(this%RootNO3DmndSoil_pvr(jroots,JZ1,JP1));this%RootNO3DmndSoil_pvr=spval
  allocate(this%RootNO3DmndBand_pvr(jroots,JZ1,JP1));this%RootNO3DmndBand_pvr=spval
  allocate(this%RootH2PO4DmndSoil_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndSoil_pvr=spval
  allocate(this%RootH2PO4DmndBand_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndBand_pvr=spval
  allocate(this%RootH1PO4DmndSoil_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndSoil_pvr=spval
  allocate(this%RootH1PO4DmndBand_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndBand_pvr=spval

  allocate(this%RootNH4DmndSoilPrev_pvr(jroots,JZ1,JP1)); this%RootNH4DmndSoilPrev_pvr=0._r8
  allocate(this%RootNH4DmndBandPrev_pvr(jroots,JZ1,JP1)); this%RootNH4DmndBandPrev_pvr=0._r8
  allocate(this%RootNO3DmndSoilPrev_pvr(jroots,JZ1,JP1)); this%RootNO3DmndSoilPrev_pvr=0._r8
  allocate(this%RootNO3DmndBandPrev_pvr(jroots,JZ1,JP1)); this%RootNO3DmndBandPrev_pvr=0._r8
  allocate(this%RootH2PO4DmndSoilPrev_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndSoilPrev_pvr=0._r8
  allocate(this%RootH2PO4DmndBandPrev_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndBandPrev_pvr=0._r8
  allocate(this%RootH1PO4DmndSoilPrev_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndSoilPrev_pvr=0._r8
  allocate(this%RootH1PO4DmndBandPrev_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndBandPrev_pvr=0._r8

  allocate(this%RAutoRootO2Limter_rpvr(jroots,JZ1,JP1));this%RAutoRootO2Limter_rpvr=spval
  allocate(this%CMinPO4Root_pft(jroots,JP1));this%CMinPO4Root_pft=spval
  allocate(this%VmaxPO4Root_pft(jroots,JP1));this%VmaxPO4Root_pft=spval
  allocate(this%KmPO4Root_pft(jroots,JP1));this%KmPO4Root_pft=spval
  allocate(this%CminNO3Root_pft(jroots,JP1));this%CminNO3Root_pft=spval
  allocate(this%VmaxNO3Root_pft(jroots,JP1));this%VmaxNO3Root_pft=spval
  allocate(this%KmNO3Root_pft(jroots,JP1));this%KmNO3Root_pft=spval
  allocate(this%CMinNH4Root_pft(jroots,JP1));this%CMinNH4Root_pft=spval
  allocate(this%VmaxNH4Root_pft(jroots,JP1));this%VmaxNH4Root_pft=spval
  allocate(this%KmNH4Root_pft(jroots,JP1));this%KmNH4Root_pft=spval
  allocate(this%RCO2Emis2Root_pvr(jroots,JZ1,JP1));this%RCO2Emis2Root_pvr=spval
  allocate(this%RootO2Uptk_pvr(jroots,JZ1,JP1));this%RootO2Uptk_pvr=spval
  allocate(this%RootUptkSoiSol_pvr(idg_beg:idg_end,jroots,JZ1,JP1));this%RootUptkSoiSol_pvr=0._r8
  end subroutine plt_rootbgc_init
!----------------------------------------------------------------------

  subroutine plt_rootbgc_destroy(this)

  implicit none
  class(plant_rootbgc_type) :: this

!  call destroy(this%trcg_rootml_pvr)
!  call destroy(this%trcs_rootml_pvr)
!  if(allocated(CO2P))deallocate(CO2P)
!  if(allocated(CO2A))deallocate(CO2A)
!  if(allocated(H2GP))deallocate(H2GP)
!  if(allocated(H2GA))deallocate(H2GA)
!  if(allocated(OXYP))deallocate(OXYP)
!  if(allocated(OXYA))deallocate(OXYA)
!  if(allocated(TRootN2Fix_pft))deallocate(TRootN2Fix_pft)
!  if(allocated(REcoUptkSoilO2M_vr))deallocate(REcoUptkSoilO2M_vr)
!  if(allocated(RootMycoExudEUptk_pvr))deallocate(RootMycoExudEUptk_pvr)
!  if(allocated(PlantRootSoilElmNetX_pft))deallocate(PlantRootSoilElmNetX_pft)
!  if(allocated(PlantExudElm_CumYr_pft))deallocate(PlantExudElm_CumYr_pft)
!  if(allocated(RootMycoExudElms_pft))deallocate(RootMycoExudElms_pft)
!  if(allocated(RootN2Fix_pft))deallocate(RootN2Fix_pft)
!  if(allocated(RootNO3Uptake_pft))deallocate(RootNO3Uptake_pft)
!  if(allocated(RootNH4Uptake_pft))deallocate(RootNH4Uptake_pft)
!  if(allocated(RootHPO4Uptake_pft))deallocate(RootHPO4Uptake_pft)
!  if(allocated(RootH2PO4Uptake_pft))deallocate(RootH2PO4Uptake_pft)
!  if(allocated(NH3Dep2Can_brch))deallocate(NH3Dep2Can_brch)
!  if(allocated(ZERO4Uptk_pft))deallocate(ZERO4Uptk_pft)
!  if(allocated(RootRespPotent_pvr))deallocate(RootRespPotent_pvr)
!  if(allocated(RootCO2EmisPot_pvr))deallocate(RootCO2EmisPot_pvr)
!  if(allocated(RootCO2Autor_pvr))deallocate(RootCO2Autor_pvr)
!  if(allocated(RCO2P))deallocate(RCO2P)
!  if(allocated(RootO2Uptk_pvr))deallocate(RootO2Uptk_pvr)
!  if(allocated(RCO2S))deallocate(RCO2S)
!  if(allocated(RO2UptkHeterS))deallocate(RO2UptkHeterS)
!  if(allocated(RUPCHS))deallocate(RUPCHS)
!  if(allocated(RUPN2S))deallocate(RUPN2S)
!  if(allocated(RUPN3S))deallocate(RUPN3S)
!  if(allocated(RUPN3B))deallocate(RUPN3B)
!  if(allocated(RUPHGS))deallocate(RUPHGS)
!  if(allocated(RCOFLA))deallocate(RCOFLA)
!  if(allocated(ROXFLA))deallocate(ROXFLA)
!  if(allocated(RCHFLA))deallocate(RCHFLA)

!  if(allocated(RootNH4BUptake_pvr))deallocate(RootNH4BUptake_pvr)
!  if(allocated(RootNH4Uptake_pvr))deallocate(RootNH4Uptake_pvr)
!  if(allocated(RootH2PO4Uptake_pvr))deallocate(RootH2PO4Uptake_pvr)
!  if(allocated(RootNO3BUptake_pvr))deallocate(RootNO3BUptake_pvr)
!  if(allocated(RootNO3Uptake_pvr))deallocate(RootNO3Uptake_pvr)
!  if(allocated(RootH1PO4BUptake_pvr))deallocate(RootH1PO4BUptake_pvr)
!  if(allocated(RootHPO4Uptake_pvr))deallocate(RootHPO4Uptake_pvr)
!  if(allocated(RootH2PO4BUptake_pvr))deallocate(RootH2PO4BUptake_pvr)
!  if(allocated(RootOUlmNutUptake_pvr))deallocate(RootOUlmNutUptake_pvr)
!  if(allocated(RootCUlmNutUptake_pvr))deallocate(RootCUlmNutUptake_pvr)



!  if(allocated(TN2FLA))deallocate(TN2FLA)
!  if(allocated(TNHFLA))deallocate(TNHFLA)
!  if(allocated(TCOFLA))deallocate(TCOFLA)
!  if(allocated(TOXFLA))deallocate(TOXFLA)
!  if(allocated(TCHFLA))deallocate(TCHFLA)
!  if(allocated(TLCH4P))deallocate(TLCH4P)
!  if(allocated(TLCO2P))deallocate(TLCO2P)
!  if(allocated(TLNH3P))deallocate(TLNH3P)

!  if(allocated(TLOXYP))deallocate(TLOXYP)
!  if(allocated(TLN2OP))deallocate(TLN2OP)

!  if(allocated(RN2FLA))deallocate(RN2FLA)
!  if(allocated(RNHFLA))deallocate(RNHFLA)
!  if(allocated(RHGFLA))deallocate(RHGFLA)
!  if(allocated(RCODFA))deallocate(RCODFA)
!  if(allocated(ROXDFA))deallocate(ROXDFA)
!  if(allocated(RCHDFA))deallocate(RCHDFA)
!  if(allocated(RN2DFA))deallocate(RN2DFA)

!  if(allocated(RNHDFA))deallocate(RNHDFA)
!  if(allocated(RHGDFA))deallocate(RHGDFA)
!  if(allocated(RootO2Dmnd4Resp_pvr))deallocate(RootO2Dmnd4Resp_pvr)
!  if(allocated(RootNH4DmndSoil_pvr))deallocate(RootNH4DmndSoil_pvr)
!  if(allocated(RootNH4DmndBand_pvr))deallocate(RootNH4DmndBand_pvr)
!  if(allocated(RootNO3DmndSoil_pvr))deallocate(RootNO3DmndSoil_pvr)
!  if(allocated(RootNO3DmndBand_pvr))deallocate(RootNO3DmndBand_pvr)
!  if(allocated(RootH2PO4DmndSoil_pvr))deallocate(RootH2PO4DmndSoil_pvr)
!  if(allocated(RootH2PO4DmndBand_pvr))deallocate(RootH2PO4DmndBand_pvr)
!  if(allocated(RootH1PO4DmndSoil_pvr))deallocate(RootH1PO4DmndSoil_pvr)
!  if(allocated(RootH1PO4DmndBand_pvr))deallocate(RootH1PO4DmndBand_pvr)
!  if(allocated(RAutoRootO2Limter_rpvr))deallocate(RAutoRootO2Limter_rpvr)
!  if(allocated(CMinPO4Root_pft))deallocate(CMinPO4Root_pft)
!  if(allocated(VmaxPO4Root_pft))deallocate(VmaxPO4Root_pft)
!  if(allocated(KmPO4Root_pft))deallocate(KmPO4Root_pft)
!  if(allocated(CminNO3Root_pft))deallocate(CminNO3Root_pft)
!  if(allocated(VmaxNO3Root_pft))deallocate(VmaxNO3Root_pft)
!  if(allocated(KmNO3Root_pft))deallocate(KmNO3Root_pft)
!  if(allocated(CMinNH4Root_pft))deallocate(CMinNH4Root_pft)
!  if(allocated(VmaxNH4Root_pft))deallocate(VmaxNH4Root_pft)
!  if(allocated(KmNH4Root_pft))deallocate(KmNH4Root_pft)
  end subroutine plt_rootbgc_destroy
!----------------------------------------------------------------------
  subroutine plt_site_Init(this)
  implicit none
  class(plant_siteinfo_type) :: this

  allocate(this%PlantElemntStoreLandscape(NumPlantChemElms));this%PlantElemntStoreLandscape=spval
  allocate(this%FracSoiAsMicP_vr(0:JZ1));this%FracSoiAsMicP_vr=spval
  allocate(this%AtmGasc(idg_beg:idg_NH3));this%AtmGasc=spval
  allocate(this%DATAP(JP1)); this%DATAP=''
  allocate(this%DATA(30)); this%DATA=''
  allocate(this%AREA3(0:JZ1));this%AREA3=spval
  allocate(this%DLYR3(0:JZ1)); this%DLYR3=spval
  allocate(this%PlantElmBalCum_pft(NumPlantChemElms,JP1));this%PlantElmBalCum_pft=spval
  allocate(this%CumSoilThickness_vr(0:JZ1));this%CumSoilThickness_vr=spval
  allocate(this%CumSoilThickMidL_vr(0:JZ1));this%CumSoilThickMidL_vr=spval
  allocate(this%PPI_pft(JP1));this%PPI_pft=spval
  allocate(this%PPatSeeding_pft(JP1));this%PPatSeeding_pft=spval
  allocate(this%PPX_pft(JP1));this%PPX_pft=spval
  allocate(this%PlantPopulation_pft(JP1));this%PlantPopulation_pft=spval
  allocate(this%flag_active_pft(JP1));  this%flag_active_pft=.false.
  allocate(this%VLWatMicPM_vr(60,0:JZ1));this%VLWatMicPM_vr=spval
  allocate(this%VLsoiAirPM_vr(60,0:JZ1));this%VLsoiAirPM_vr=spval
  allocate(this%TortMicPM_vr(60,0:JZ1));this%TortMicPM_vr=spval
  allocate(this%FILMM_vr(60,0:JZ1)); this%FILMM_vr=spval

  end subroutine plt_site_Init
!----------------------------------------------------------------------
  subroutine plt_site_destroy(this)
  implicit none
  class(plant_siteinfo_type) :: this



  end subroutine plt_site_destroy
!----------------------------------------------------------------------

  subroutine plt_bgcrate_init(this)
  implicit none
  class(plant_bgcrate_type) :: this

  allocate(this%RootCO2Emis2Root_vr(JZ1)); this%RootCO2Emis2Root_vr=spval
  allocate(this%RUptkRootO2_vr(JZ1)); this%RUptkRootO2_vr=0._r8
  allocate(this%RH2PO4EcoDmndSoilPrev_vr(0:JZ1)); this%RH2PO4EcoDmndSoilPrev_vr=spval
  allocate(this%RH2PO4EcoDmndBandPrev_vr(0:JZ1)); this%RH2PO4EcoDmndBandPrev_vr=spval
  allocate(this%RH1PO4EcoDmndSoilPrev_vr(0:JZ1)); this%RH1PO4EcoDmndSoilPrev_vr=spval
  allocate(this%RH1PO4EcoDmndBandPrev_vr(0:JZ1)); this%RH1PO4EcoDmndBandPrev_vr=spval
  allocate(this%RNO3EcoDmndSoilPrev_vr(0:JZ1)); this%RNO3EcoDmndSoilPrev_vr=spval
  allocate(this%RNH4EcoDmndSoilPrev_vr(0:JZ1)); this%RNH4EcoDmndSoilPrev_vr =spval
  allocate(this%RNH4EcoDmndBandPrev_vr(0:JZ1)); this%RNH4EcoDmndBandPrev_vr=spval
  allocate(this%RNO3EcoDmndBandPrev_vr(0:JZ1)); this%RNO3EcoDmndBandPrev_vr=spval
  allocate(this%RGasTranspFlxPrev_vr(idg_beg:idg_end,0:JZ1)); this%RGasTranspFlxPrev_vr=spval
  allocate(this%RO2AquaSourcePrev_vr(0:JZ1)); this%RO2AquaSourcePrev_vr=spval
  allocate(this%RO2EcoDmndPrev_vr(0:JZ1)); this%RO2EcoDmndPrev_vr=spval
  allocate(this%LitrfalStrutElms_vr(NumPlantChemElms,jsken,NumOfPlantLitrCmplxs,0:JZ1));this%LitrfalStrutElms_vr=spval
  allocate(this%GrossCO2Fix_pft(JP1));this%GrossCO2Fix_pft=spval
  allocate(this%REcoDOMProd_vr(1:NumPlantChemElms,1:jcplx,0:JZ1));this%REcoDOMProd_vr=spval
  allocate(this%CO2NetFix_pft(JP1));this%CO2NetFix_pft=spval
  allocate(this%RCanMaintDef_CO2_pft(JP1));this%RCanMaintDef_CO2_pft=spval
  allocate(this%RootGasLossDisturb_pft(idg_beg:idg_NH3,JP1));this%RootGasLossDisturb_pft=spval
  allocate(this%GrossResp_pft(JP1));this%GrossResp_pft=spval
  allocate(this%CanopyGrosRCO2_pft(JP1));this%CanopyGrosRCO2_pft=spval
  allocate(this%CanopyResp_brch(MaxNumBranches,JP1));this%CanopyResp_brch=spval
  allocate(this%RootAutoCO2_pft(JP1));this%RootAutoCO2_pft=spval
  allocate(this%PlantN2Fix_CumYr_pft(JP1));this%PlantN2Fix_CumYr_pft=spval
  allocate(this%NH3Emis_CumYr_pft(JP1));this%NH3Emis_CumYr_pft=spval
  allocate(this%NodulInfectElms_pft(NumPlantChemElms,JP1));this%NodulInfectElms_pft=spval
  allocate(this%SurfLitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP1));this%SurfLitrfalStrutElms_CumYr_pft=0._r8
  allocate(this%LitrFallStrutElms_col(NumPlantChemElms));this%LitrFallStrutElms_col=0._r8
  allocate(this%NetPrimProduct_pft(JP1));this%NetPrimProduct_pft=spval
  allocate(this%Nutruptk_fClim_rpvr(jroots,JZ1,JP1));this%Nutruptk_fClim_rpvr=0._r8
  allocate(this%Nutruptk_fNlim_rpvr(jroots,JZ1,JP1));this%Nutruptk_fNlim_rpvr=0._r8
  allocate(this%Nutruptk_fPlim_rpvr(jroots,JZ1,JP1));this%Nutruptk_fPlim_rpvr=0._r8
  allocate(this%Nutruptk_fProtC_rpvr(jroots,JZ1,JP1));this%Nutruptk_fProtC_rpvr=0._r8
  allocate(this%NH3Dep2Can_pft(JP1));this%NH3Dep2Can_pft=spval
  allocate(this%tRootMycoExud2Soil_vr(NumPlantChemElms,1:jcplx,JZ1));this%tRootMycoExud2Soil_vr=spval
  allocate(this%RootN2Fix_pvr(JZ1,JP1));this%RootN2Fix_pvr=0._r8
  allocate(this%RCO2Nodule_pvr(JZ1,JP1)); this%RCO2Nodule_pvr=0._r8
  allocate(this%RootO2_TotSink_pvr(jroots,JZ1,JP1)); this%RootO2_TotSink_pvr=0._r8
  allocate(this%RootO2_TotSink_vr(JZ1)); this%RootO2_TotSink_vr=0._r8
  allocate(this%CanopyRespC_CumYr_pft(JP1));this%CanopyRespC_CumYr_pft=spval
  allocate(this%CanopyNLimFactor_brch(MaxNumBranches,JP1));this%CanopyNLimFactor_brch=1._r8
  allocate(this%CanopyPLimFactor_brch(MaxNumBranches,JP1));this%CanopyPLimFactor_brch=1._r8
  allocate(this%REcoH1PO4DmndBand_vr(0:JZ1));this%REcoH1PO4DmndBand_vr=spval
  allocate(this%REcoNO3DmndSoil_vr(0:JZ1));this%REcoNO3DmndSoil_vr=spval
  allocate(this%REcoH2PO4DmndSoil_vr(0:JZ1));this%REcoH2PO4DmndSoil_vr=spval
  allocate(this%REcoNH4DmndSoil_vr(0:JZ1));this%REcoNH4DmndSoil_vr=spval
  allocate(this%REcoO2DmndResp_vr(0:JZ1));this%REcoO2DmndResp_vr=spval
  allocate(this%REcoH2PO4DmndBand_vr(0:JZ1));this%REcoH2PO4DmndBand_vr=spval
  allocate(this%REcoNO3DmndBand_vr(0:JZ1));this%REcoNO3DmndBand_vr=spval
  allocate(this%REcoNH4DmndBand_vr(0:JZ1));this%REcoNH4DmndBand_vr=spval
  allocate(this%REcoH1PO4DmndSoil_vr(0:JZ1));this%REcoH1PO4DmndSoil_vr=spval
  allocate(this%LitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP1));this%LitrfalStrutElms_CumYr_pft=0._r8
  allocate(this%LitrfallElms_pft(NumPlantChemElms,JP1));this%LitrfallElms_pft=spval
  allocate(this%LitrfallElms_pvr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ1,JP1));this%LitrfallElms_pvr=spval
  allocate(this%LitrFallElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LitrFallElms_brch=spval
  allocate(this%LitrfallAbvgElms_pft(NumPlantChemElms,JP1)); this%LitrfallAbvgElms_pft=spval
  allocate(this%LitrfallBlgrElms_pft(NumPlantChemElms,JP1)); this%LitrfallBlgrElms_pft=spval
  allocate(this%RootMaintDef_CO2_pvr(jroots,JZ1,JP1));this%RootMaintDef_CO2_pvr=0._r8  

  end subroutine plt_bgcrate_init
!----------------------------------------------------------------------
  subroutine plt_bgcrate_destroy(this)

  implicit none
  class(plant_bgcrate_type) :: this


  end subroutine plt_bgcrate_destroy
!----------------------------------------------------------------------
  subroutine plt_disturb_init(this)

  implicit none
  class(plant_disturb_type) :: this

  allocate(this%THIN_pft(JP1));this%THIN_pft=spval
  allocate(this%EcoHavstElmnt_CumYr_col(NumPlantChemElms));this%EcoHavstElmnt_CumYr_col=spval
  allocate(this%EcoHavstElmntCum_pft(NumPlantChemElms,JP1));this%EcoHavstElmntCum_pft=spval
  allocate(this%EcoHavstElmnt_CumYr_pft(NumPlantChemElms,JP1));this%EcoHavstElmnt_CumYr_pft=spval
  allocate(this%CH4ByFire_CumYr_pft(JP1));this%CH4ByFire_CumYr_pft=spval
  allocate(this%CO2ByFire_CumYr_pft(JP1));this%CO2ByFire_CumYr_pft=spval
  allocate(this%PlantElmDistLoss_pft(NumPlantChemElms,JP1)); this%PlantElmDistLoss_pft=spval
  allocate(this%N2ObyFire_CumYr_pft(JP1));this%N2ObyFire_CumYr_pft=spval
  allocate(this%NH3byFire_CumYr_pft(JP1));this%NH3byFire_CumYr_pft=spval
  allocate(this%PO4byFire_CumYr_pft(JP1));this%PO4byFire_CumYr_pft=spval

  allocate(this%iHarvestDay_pft(JP1)); this%iHarvestDay_pft=0
  allocate(this%O2ByFire_CumYr_pft(JP1));this%O2ByFire_CumYr_pft=spval
  allocate(this%FracBiomHarvsted(1:2,1:4,JP1));this%FracBiomHarvsted=spval
  allocate(this%CanopyHeightCut_pft(JP1)); this%CanopyHeightCut_pft=spval
  allocate(this%iYearPlantHarvest_pft(JP1));this%iYearPlantHarvest_pft=0
  allocate(this%FERT(1:20));this%FERT=spval
  allocate(this%iYearPlanting_pft(JP1));this%iYearPlanting_pft=0
  allocate(this%iPlantingDay_pft(JP1));this%iPlantingDay_pft=0
  allocate(this%iPlantingYear_pft(JP1));this%iPlantingYear_pft=0
  allocate(this%iHarvestYear_pft(JP1));this%iHarvestYear_pft=0
  allocate(this%iDayPlanting_pft(JP1));this%iDayPlanting_pft=0
  allocate(this%iDayPlantHarvest_pft(JP1));this%iDayPlantHarvest_pft=0
  allocate(this%iHarvstType_pft(JP1));this%iHarvstType_pft=-1
  allocate(this%jHarvstType_pft(JP1));this%jHarvstType_pft=0

  end subroutine plt_disturb_init
!----------------------------------------------------------------------
  subroutine plt_disturb_destroy(this)
  implicit none
  class(plant_disturb_type) :: this


!  if(allocated(THIN_pft))deallocate(THIN_pft)
!  if(allocated(EcoHavstElmnt_CumYr_pft))deallocate(EcoHavstElmnt_CumYr_pft)
!  if(allocated(EcoHavstElmntCum_pft))deallocate(EcoHavstElmntCum_pft)
!  if(allocated(CH4ByFire_CumYr_pft))deallocate(CH4ByFire_CumYr_pft)
!  if(allocated(CO2ByFire_CumYr_pft))deallocate(CO2ByFire_CumYr_pft)
!  if(allocated(N2ObyFire_CumYr_pft))deallocate(N2ObyFire_CumYr_pft)
!  if(allocated(NH3byFire_CumYr_pft))deallocate(NH3byFire_CumYr_pft)
!  if(allocated(PO4byFire_CumYr_pft))deallocate(PO4byFire_CumYr_pft)

!  if(allocated(iHarvestDay_pft))deallocate(iHarvestDay_pft)
!  if(allocated(O2ByFire_CumYr_pft))deallocate(O2ByFire_CumYr_pft)
!  if(allocated(FracBiomHarvsted))deallocate(FracBiomHarvsted)
!  if(allocated(HVST))deallocate(HVST)
!  if(allocated(iYearPlantHarvest_pft))deallocate(iYearPlantHarvest_pft)
!  if(allocated(iDayPlanting_pft))deallocate(iDayPlanting_pft)
!  if(allocated(iDayPlantHarvest_pft))deallocate(iDayPlantHarvest_pft)
!  if(allocated(iHarvestYear_pft))deallocate(iHarvestYear_pft)
!  if(allocated(iPlantingDay_pft))deallocate(iPlantingDay_pft)
!  if(allocated(iYearPlanting_pft))deallocate(iYearPlanting_pft)
!  if(allocated(FERT))deallocate(FERT)
!  if(allocated(iPlantingYear_pft))deallocate(iPlantingYear_pft)
!  if(allocated(iHarvstType_pft))deallocate(iHarvstType_pft)
!  if(allocated(jHarvstType_pft))deallocate(jHarvstType_pft)

  end subroutine plt_disturb_destroy
!----------------------------------------------------------------------
  subroutine  plt_ew_init(this)

  implicit none
  class(plant_ew_type) :: this

  allocate(this%RootResist4H2O_pvr(jroots,JZ1,JP1)); this%RootResist4H2O_pvr=spval
  allocate(this%RootRadialKond2H2O_pvr(jroots,JZ1,JP1));this%RootRadialKond2H2O_pvr=spval
  allocate(this%RootAxialKond2H2O_pvr(jroots,JZ1,JP1));this%RootAxialKond2H2O_pvr=spval
  allocate(this%QdewCanopy_pft(JP1)); this%QdewCanopy_pft=spval
  allocate(this%ETCanopy_CumYr_pft(JP1));this%ETCanopy_CumYr_pft=spval
  allocate(this%ElvAdjstedSoilH2OPSIMPa_vr(0:JZ1));this%ElvAdjstedSoilH2OPSIMPa_vr=spval
  allocate(this%THeatLossRoot2Soil_vr(0:JZ1));this%THeatLossRoot2Soil_vr=spval
  allocate(this%TKCanopy_pft(JP1));this%TKCanopy_pft=spval
  allocate(this%HeatXAir2PCan_pft(JP1));this%HeatXAir2PCan_pft=spval
  allocate(this%PrecIntcptByCanopy_pft(JP1));this%PrecIntcptByCanopy_pft=spval
  allocate(this%PSICanopyTurg_pft(JP1));this%PSICanopyTurg_pft=spval
  allocate(this%PSICanopy_pft(JP1));this%PSICanopy_pft=spval
  allocate(this%VapXAir2Canopy_pft(JP1));this%VapXAir2Canopy_pft=spval
  allocate(this%HeatStorCanopy_pft(JP1));this%HeatStorCanopy_pft=spval
  allocate(this%EvapTransLHeat_pft(JP1));this%EvapTransLHeat_pft=spval
  allocate(this%WatHeldOnCanopy_pft(JP1));this%WatHeldOnCanopy_pft=spval
  allocate(this%VHeatCapCanopy_pft(JP1));this%VHeatCapCanopy_pft=spval
  allocate(this%CanopyBiomWater_pft(JP1));this%CanopyBiomWater_pft=spval
  allocate(this%PSIRoot_pvr(jroots,JZ1,JP1));this%PSIRoot_pvr=spval
  allocate(this%PSIRootOSMO_vr(jroots,JZ1,JP1));this%PSIRootOSMO_vr=spval
  allocate(this%PSIRootTurg_vr(jroots,JZ1,JP1));this%PSIRootTurg_vr=spval
  allocate(this%RPlantRootH2OUptk_pvr(jroots,JZ1,JP1));this%RPlantRootH2OUptk_pvr=spval
  allocate(this%RootH2OUptkStress_pvr(jroots,JZ1,JP1)); this%RootH2OUptkStress_pvr=spval
  allocate(this%TWaterPlantRoot2Soil_vr(0:JZ1));this%TWaterPlantRoot2Soil_vr=spval
  allocate(this%Transpiration_pft(JP1));this%Transpiration_pft=spval
  allocate(this%PSICanopyOsmo_pft(JP1));this%PSICanopyOsmo_pft=spval
  allocate(this%TKS_vr(0:JZ1));this%TKS_vr=spval
  allocate(this%CanOsmoPsi0pt_pft(JP1));this%CanOsmoPsi0pt_pft=spval
  allocate(this%CanopyIsothBndlResist_pft(JP1));this%CanopyIsothBndlResist_pft=spval
  allocate(this%DeltaTKC_pft(JP1));this%DeltaTKC_pft=spval
  allocate(this%TKC_pft(JP1));this%TKC_pft=spval
  allocate(this%ENGYX_pft(JP1));this%ENGYX_pft=spval
  allocate(this%TdegCCanopy_pft(JP1));this%TdegCCanopy_pft=spval
  allocate(this%PSICanPDailyMin_pft(JP1));this%PSICanPDailyMin_pft=spval

  end subroutine plt_ew_init
!----------------------------------------------------------------------

  subroutine plt_ew_destroy(this)
  implicit none
  class(plant_ew_type) :: this



  end subroutine plt_ew_destroy

!----------------------------------------------------------------------

  subroutine plt_allom_init(this)
  implicit none
  class(plant_allometry_type) :: this


  allocate(this%RootFracRemobilizableBiom_pft(JP1));this%RootFracRemobilizableBiom_pft=spval
  allocate(this%FracGroth2Node_pft(JP1));this%FracGroth2Node_pft=spval
  allocate(this%GrainSeedBiomCMean_brch(MaxNumBranches,JP1));this%GrainSeedBiomCMean_brch=spval
  allocate(this%NoduGrowthYield_pft(JP1));this%NoduGrowthYield_pft=spval
  allocate(this%RootBiomGrosYld_pft(JP1));this%RootBiomGrosYld_pft=spval
  allocate(this%rPCRootr_pft(JP1));this%rPCRootr_pft=spval
  allocate(this%rProteinC2N_pft(JP1));this%rProteinC2N_pft=spval
  allocate(this%rProteinC2P_pft(JP1));this%rProteinC2P_pft=spval
  allocate(this%CPRTS_pft(JP1));this%CPRTS_pft=spval
  allocate(this%CNRTS_pft(JP1));this%CNRTS_pft=spval
  allocate(this%rNCNodule_pft(JP1));this%rNCNodule_pft=spval
  allocate(this%rPCNoduler_pft(JP1));this%rPCNoduler_pft=spval
  allocate(this%rNCRoot_pft(JP1));this%rNCRoot_pft=spval
  allocate(this%rPCLeaf_pft(JP1));this%rPCLeaf_pft=spval
  allocate(this%rPCSheath_pft(JP1));this%rPCSheath_pft=spval
  allocate(this%rNCLeaf_pft(JP1));this%rNCLeaf_pft=spval
  allocate(this%rNCSheath_pft(JP1));this%rNCSheath_pft=spval
  allocate(this%rNCGrain_pft(JP1));this%rNCGrain_pft=spval
  allocate(this%rPCStalk_pft(JP1));this%rPCStalk_pft=spval
  allocate(this%rNCStalk_pft(JP1));this%rNCStalk_pft=spval
  allocate(this%rPCGrain_pft(JP1));this%rPCGrain_pft=spval
  allocate(this%rPCEar_pft(JP1));this%rPCEar_pft=spval
  allocate(this%rPCReserve_pft(JP1));this%rPCReserve_pft=spval
  allocate(this%rNCReserve_pft(JP1));this%rNCReserve_pft=spval
  allocate(this%rPCHusk_pft(JP1));this%rPCHusk_pft=spval
  allocate(this%FracWoodStalkElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracWoodStalkElmAlloc2Litr=spval
  allocate(this%FracShootPetolAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracShootPetolAlloc2Litr=spval
  allocate(this%FracRootElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracRootElmAlloc2Litr=spval
  allocate(this%FracShootLeafAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracShootLeafAlloc2Litr=spval

  allocate(this%PetioleBiomGrowthYld_pft(JP1));this%PetioleBiomGrowthYld_pft=spval
  allocate(this%HuskBiomGrowthYld_pft(JP1));this%HuskBiomGrowthYld_pft=spval
  allocate(this%StalkBiomGrowthYld_pft(JP1));this%StalkBiomGrowthYld_pft=spval
  allocate(this%ReserveBiomGrowthYld_pft(JP1));this%ReserveBiomGrowthYld_pft=spval
  allocate(this%EarBiomGrowthYld_pft(JP1));this%EarBiomGrowthYld_pft=spval
  allocate(this%GrainBiomGrowthYld_pft(JP1));this%GrainBiomGrowthYld_pft=spval
  allocate(this%rNCHusk_pft(JP1));this%rNCHusk_pft=spval
  allocate(this%rNCEar_pft(JP1));this%rNCEar_pft=spval
  allocate(this%LeafBiomGrowthYld_pft(JP1));this%LeafBiomGrowthYld_pft=spval

  end subroutine plt_allom_init
!----------------------------------------------------------------------

  subroutine plt_allom_destroy(this)
  implicit none

  class(plant_allometry_type) :: this

!  if(allocated(FracGroth2Node_pft))deallocate(FracGroth2Node_pft)
!  if(allocated(GrainSeedBiomCMean_brch))deallocate(GrainSeedBiomCMean_brch)
!  if(allocated(NoduGrowthYield_pft))deallocate(NoduGrowthYield_pft)
!  if(allocated(RootBiomGrosYld_pft))deallocate(RootBiomGrosYld_pft)
!  if(allocated(rPCRootr_pft))deallocate(rPCRootr_pft)
!  if(allocated(rProteinC2N_pft))deallocate(rProteinC2N_pft)
!  if(allocated(rProteinC2P_pft))deallocate(rProteinC2P_pft)
!  if(allocated(CPRTS_pft))deallocate(CPRTS_pft)
!  if(allocated(CNRTS_pft))deallocate(CNRTS_pft)
!  if(allocated(rNCNodule_pft))deallocate(rNCNodule_pft)
!  if(allocated(rPCNoduler_pft))deallocate(rPCNoduler_pft)
!  if(allocated(rNCRoot_pft))deallocate(rNCRoot_pft)
!  if(allocated(RootFracRemobilizableBiom_pft))deallocate(RootFracRemobilizableBiom_pft)
!  if(allocated(CPLF))deallocate(CPLF)
!  if(allocated(CPSHE))deallocate(CPSHE)
!  if(allocated(CNSHE))deallocate(CNSHE)
!  if(allocated(CNLF))deallocate(CNLF)
!  if(allocated(LeafBiomGrowthYld_pft))deallocate(LeafBiomGrowthYld_pft)
!  if(allocated(rPCStalk_pft))deallocate(rPCStalk_pft)
!  if(allocated(CNGR))deallocate(CNGR)
!  if(allocated(rNCStalk_pft))deallocate(rNCStalk_pft)
!  if(allocated(CPGR))deallocate(CPGR)
!  if(allocated(PetioleBiomGrowthYld_pft))deallocate(PetioleBiomGrowthYld_pft)
!  if(allocated(StalkBiomGrowthYld_pft))deallocate(StalkBiomGrowthYld_pft)
!  if(allocated(GrainBiomGrowthYld_pft))deallocate(GrainBiomGrowthYld_pft)
!  if(allocated(ReserveBiomGrowthYld_pft))deallocate(ReserveBiomGrowthYld_pft)
!  if(allocated(EarBiomGrowthYld_pft))deallocate(EarBiomGrowthYld_pft)
!  if(allocated(HuskBiomGrowthYld_pft))deallocate(HuskBiomGrowthYld_pft)
!  if(allocated(FracShootPetolAlloc2Litr))deallocate(FracShootPetolAlloc2Litr)
!  if(allocated(FracWoodStalkElmAlloc2Litr))deallocate(FracWoodStalkElmAlloc2Litr)
!  if(allocated(FracRootElmAlloc2Litr))deallocate(FracRootElmAlloc2Litr)
!  if(allocated(FracShootLeafAlloc2Litr))deallocate(FracShootLeafAlloc2Litr)
!  if(allocated(rPCEar_pft))deallocate(rPCEar_pft)
!  if(allocated(rPCHusk_pft))deallocate(rPCHusk_pft)
!  if(allocated(rNCHusk_pft))deallocate(rNCHusk_pft)
!  if(allocated(rNCEar_pft))deallocate(rNCEar_pft)
!  if(allocated(rNCReserve_pft))deallocate(rNCReserve_pft)
!  if(allocated(rPCReserve_pft))deallocate(rPCReserve_pft)

  end subroutine plt_allom_destroy
!----------------------------------------------------------------------
  subroutine plt_biom_init(this)
  implicit none
  class(plant_biom_type) :: this

  allocate(this%StomatalStress_pft(JP1));  this%StomatalStress_pft=spval
  allocate(this%ZERO4LeafVar_pft(JP1));this%ZERO4LeafVar_pft=spval
  allocate(this%ZERO4Groth_pft(JP1));this%ZERO4Groth_pft=spval
  allocate(this%StandingDeadStrutElms_col(NumPlantChemElms));this%StandingDeadStrutElms_col=spval
  allocate(this%RootNodulStrutElms_rpvr(NumPlantChemElms,JZ1,JP1));this%RootNodulStrutElms_rpvr=spval
  allocate(this%CanopyLeafCLyr_pft(NumCanopyLayers1,JP1));this%CanopyLeafCLyr_pft=spval
  allocate(this%RootNodulNonstElms_rpvr(NumPlantChemElms,JZ1,JP1));this%RootNodulNonstElms_rpvr=spval
  allocate(this%StandDeadKCompElms_pft(NumPlantChemElms,jsken,JP1));this%StandDeadKCompElms_pft=0._r8
  allocate(this%RootMyco2ndStrutElms_rpvr(NumPlantChemElms,jroots,JZ1,MaxNumRootAxes,JP1))
  this%RootMyco2ndStrutElms_rpvr=spval
  allocate(this%RootMycoNonstElms_pft(NumPlantChemElms,jroots,JP1));this%RootMycoNonstElms_pft=spval
  allocate(this%RootMyco1stStrutElms_rpvr(NumPlantChemElms,jroots,JZ1,MaxNumRootAxes,JP1))
  this%RootMyco1stStrutElms_rpvr=spval
  allocate(this%CanopyNonstElmConc_pft(NumPlantChemElms,JP1));this%CanopyNonstElmConc_pft=spval
  allocate(this%CanopyNonstElms_pft(NumPlantChemElms,JP1));this%CanopyNonstElms_pft=spval
  allocate(this%CanopyNodulNonstElms_pft(NumPlantChemElms,JP1));this%CanopyNodulNonstElms_pft=spval
  allocate(this%CanopyNoduleNonstCConc_pft(JP1));this%CanopyNoduleNonstCConc_pft=spval
  allocate(this%RootProteinConc_rpvr(jroots,JZ1,JP1));this%RootProteinConc_rpvr=spval
  allocate(this%RootProteinC_pvr(jroots,JZ1,JP1));this%RootProteinC_pvr=spval
  allocate(this%RootMycoActiveBiomC_pvr(jroots,JZ1,JP1));this%RootMycoActiveBiomC_pvr=spval
  allocate(this%RootMycoMassElm_pvr(NumPlantChemElms,jroots,JZ1,JP1)); this%RootMycoMassElm_pvr = 0._r8
  allocate(this%PopuRootMycoC_pvr(jroots,JZ1,JP1));this%PopuRootMycoC_pvr=spval
  allocate(this%RootMycoNonstElms_rpvr(NumPlantChemElms,jroots,JZ1,JP1));this%RootMycoNonstElms_rpvr=spval
  allocate(this%RootNonstructElmConc_rpvr(NumPlantChemElms,jroots,JZ1,JP1));this%RootNonstructElmConc_rpvr=spval
  allocate(this%ROOTNLim_rpvr(jroots,JZ1,JP1)); this%ROOTNLim_rpvr=0._r8
  allocate(this%ROOTPLim_rpvr(jroots,JZ1,JP1)); this%ROOTPLim_rpvr=0._r8
  allocate(this%CanopyNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNonstElms_brch=spval
  allocate(this%C4PhotoShootNonstC_brch(MaxNumBranches,JP1));this%C4PhotoShootNonstC_brch=spval
  allocate(this%CanopySapwoodC_pft(JP1));this%CanopySapwoodC_pft=spval
  allocate(this%ShootElms_pft(NumPlantChemElms,JP1));this%ShootElms_pft=spval
  allocate(this%ShootElmsBeg_pft(NumPlantChemElms,JP1));this%ShootElmsBeg_pft=spval
  allocate(this%LeafProteinC_node(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafProteinC_node=spval
  allocate(this%PetoleProteinCNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%PetoleProteinCNode_brch=spval
  allocate(this%StructInternodeElms_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%StructInternodeElms_brch=spval
  allocate(this%LeafElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafElmntNode_brch=spval
  allocate(this%PetioleElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%PetioleElmntNode_brch=spval
  allocate(this%LeafLayerElms_node(NumPlantChemElms,NumCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafLayerElms_node=spval
  allocate(this%LeafProteinC_brch(MaxNumBranches,JP1)); this%LeafProteinC_brch=spval
  allocate(this%LeafProteinCperm2LA_pft(JP1)); this%LeafProteinCperm2LA_pft=spval
  allocate(this%SapwoodBiomassC_brch(MaxNumBranches,JP1));this%SapwoodBiomassC_brch=spval
  allocate(this%CanopyNodulNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNodulNonstElms_brch=spval
  allocate(this%LeafPetoNonstElmConc_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafPetoNonstElmConc_brch=spval
  allocate(this%RootStrutElms_pft(NumPlantChemElms,JP1));this%RootStrutElms_pft=spval
  allocate(this%StandDeadStrutElmsBeg_pft(NumPlantChemElms,JP1));this%StandDeadStrutElmsBeg_pft=spval
  allocate(this%tCanLeafC_clyr(NumCanopyLayers1));this%tCanLeafC_clyr=spval
  allocate(this%RootElms_pft(NumPlantChemElms,JP1));this%RootElms_pft=spval
  allocate(this%RootNoduleElms_pft(NumPlantChemElms,JP1));this%RootNoduleElms_pft=spval
  allocate(this%RootNoduleElmsBeg_pft(NumPlantChemElms,JP1));this%RootNoduleElmsBeg_pft=spval  
  allocate(this%RootElmsBeg_pft(NumPlantChemElms,JP1));this%RootElmsBeg_pft=spval
  allocate(this%SeedPlantedElm_pft(NumPlantChemElms,JP1));this%SeedPlantedElm_pft=spval
  allocate(this%SeasonalNonstElms_pft(NumPlantChemElms,JP1));this%SeasonalNonstElms_pft=spval
  allocate(this%SeasonalNonstElmsbeg_pft(NumPlantChemElms,JP1));this%SeasonalNonstElmsbeg_pft=spval
  allocate(this%TotBegVegE_pft(NumPlantChemElms,JP1)); this%TotBegVegE_pft=spval
  allocate(this%TotEndVegE_pft(NumPlantChemElms,JP1)); this%TotEndVegE_pft=spval
  allocate(this%CanopyLeafSheathC_pft(JP1));this%CanopyLeafSheathC_pft=spval
  allocate(this%StandDeadStrutElms_pft(NumPlantChemElms,JP1));this%StandDeadStrutElms_pft=spval
  allocate(this%RootBiomCPerPlant_pft(JP1));this%RootBiomCPerPlant_pft=spval
  allocate(this%PetoleStrutElms_pft(NumPlantChemElms,JP1));this%PetoleStrutElms_pft=spval
  allocate(this%StalkStrutElms_pft(NumPlantChemElms,JP1));this%StalkStrutElms_pft=spval
  allocate(this%StalkRsrvElms_pft(NumPlantChemElms,JP1));this%StalkRsrvElms_pft=spval
  allocate(this%GrainStrutElms_pft(NumPlantChemElms,JP1));this%GrainStrutElms_pft=0._r8
  allocate(this%HuskStrutElms_pft(NumPlantChemElms,JP1));this%HuskStrutElms_pft=spval
  allocate(this%EarStrutElms_pft(NumPlantChemElms,JP1));this%EarStrutElms_pft=spval
  allocate(this%CanopyMassC_pft(JP1)); this%CanopyMassC_pft=0._r8
  allocate(this%LeafChemElmRemob_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafChemElmRemob_brch=spval
  allocate(this%PetioleChemElmRemob_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetioleChemElmRemob_brch=spval
  allocate(this%ShootElms_pft(NumPlantChemElms,JP1));this%ShootElms_pft=spval
  allocate(this%AvgCanopyBiomC2Graze_pft(JP1));this%AvgCanopyBiomC2Graze_pft=spval
  allocate(this%LeafStrutElms_pft(NumPlantChemElms,JP1));this%LeafStrutElms_pft=spval
  allocate(this%StandingDeadInitC_pft(JP1));this%StandingDeadInitC_pft=spval
  allocate(this%CanopyLeafSheathC_brch(MaxNumBranches,JP1));this%CanopyLeafSheathC_brch=spval
  allocate(this%StalkRsrvElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%StalkRsrvElms_brch=spval
  allocate(this%LeafStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafStrutElms_brch=spval
  allocate(this%CanopyNodulStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNodulStrutElms_brch=spval
  allocate(this%PetoleStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetoleStrutElms_brch=spval
  allocate(this%EarStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%EarStrutElms_brch=spval
  allocate(this%HuskStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%HuskStrutElms_brch=spval
  allocate(this%GrainStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%GrainStrutElms_brch=0._r8
  allocate(this%StalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%StalkStrutElms_brch=spval
  allocate(this%ShootElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%ShootElms_brch=spval
  allocate(this%SenecStalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%SenecStalkStrutElms_brch=spval
  allocate(this%RootMyco1stElm_raxs(NumPlantChemElms,jroots,MaxNumRootAxes,JP1));this%RootMyco1stElm_raxs=spval

  allocate(this%LeafC3ChlC_brch(MaxNumBranches,JP1));this%LeafC3ChlC_brch=0._r8             
  allocate(this%LeafC4ChlC_brch(MaxNumBranches,JP1));this%LeafC4ChlC_brch=0._r8             
  allocate(this%LeafRubiscoC_brch(MaxNumBranches,JP1));this%LeafRubiscoC_brch=0._r8           
  allocate(this%LeafPEPC_brch(MaxNumBranches,JP1));this%LeafPEPC_brch=0._r8               
  allocate(this%LeafC3ChlCperm2LA_pft(JP1))     ;this%LeafC3ChlCperm2LA_pft=0._r8           
  allocate(this%LeafC4ChlCperm2LA_pft(JP1))     ;this%LeafC4ChlCperm2LA_pft=0._r8                      
  allocate(this%LeafRubiscoCperm2LA_pft(JP1))   ;this%LeafRubiscoCperm2LA_pft=0._r8           
  allocate(this%LeafPEPCperm2LA_pft(JP1))       ;this%LeafPEPCperm2LA_pft=0._r8           
  allocate(this%SpecificLeafArea_pft(JP1)) ; this%SpecificLeafArea_pft=0._r8
  allocate(this%ShootNoduleElms_pft(NumPlantChemElms,JP1)); this%ShootNoduleElms_pft=0._r8
  allocate(this%ShootNoduleElmsBeg_pft(NumPlantChemElms,JP1));this%ShootNoduleElmsBeg_pft=0._r8
  end subroutine plt_biom_init
!----------------------------------------------------------------------

  subroutine plt_biom_destroy(this)

  implicit none
  class(plant_biom_type) :: this

!  if(allocated(ZERO4LeafVar_pft))deallocate(ZERO4LeafVar_pft)
!  if(allocated(ZERO4Groth_pft))deallocate(ZERO4Groth_pft)
!  if(allocated(RootNodulNonstElms_rpvr))deallocate(RootNodulNonstElms_rpvr)
!  if(allocated(RootNodulStrutElms_rpvr))deallocate(RootNodulStrutElms_rpvr)
!  if(allocated(CanopyLeafCLyr_pft))deallocate(CanopyLeafCLyr_pft)
!  call destroy(RootMyco1stElm_raxs)
!  if(allocated(RootMyco1stStrutElms_rpvr))deallocate(RootMyco1stStrutElms_rpvr)
!  if(allocated(StandDeadKCompElms_pft))deallocate(StandDeadKCompElms_pft)
!  if(allocated(RootMyco2ndStrutElms_rpvr))deallocate(RootMyco2ndStrutElms_rpvr)
!  if(allocated(CanopyNonstElmConc_pft))deallocate(CanopyNonstElmConc_pft)
!  if(allocated(CanopyNonstElms_pft))deallocate(CanopyNonstElms_pft)
!  if(allocated(CanopyNodulNonstElms_pft))deallocate(CanopyNodulNonstElms_pft)
!  if(allocated(CanopyNoduleNonstCConc_pft))deallocate(CanopyNoduleNonstCConc_pft)
!  if(allocated(RootProteinConc_rpvr))deallocate(RootProteinConc_rpvr)
!  if(allocated(RootProteinC_pvr))deallocate(RootProteinC_pvr)
!  if(allocated(RootMycoActiveBiomC_pvr))deallocate(RootMycoActiveBiomC_pvr)
!  if(allocated( PopuRootMycoC_pvr))deallocate( PopuRootMycoC_pvr)
!  if(allocated(RootMycoNonstElms_rpvr))deallocate(RootMycoNonstElms_rpvr)
!  if(allocated(WVSTK))deallocate(WVSTK)
!  if(allocated(LeafElmntNode_brch))deallocate(LeafElmntNode_brch)
!  if(allocated(LeafProteinC_node))deallocate(LeafProteinC_node)
!  if(allocated(PetoleProteinCNode_brch))deallocate(PetoleProteinCNode_brch)
!  if(allocated(LeafLayerElms_node))deallocate(LeafLayerElms_node)
!  if(allocated(StructInternodeElms_brch))deallocate(StructInternodeElms_brch)
!  if(allocated(PetioleElmntNode_brch))deallocate(PetioleElmntNode_brch)
!  if(allocated(SapwoodBiomassC_brch))deallocate(SapwoodBiomassC_brch)
!  if(allocated(PPOOL))deallocate(PPOOL)
!  if(allocated(CanopyNodulNonstElms_brch))deallocate(CanopyNodulNonstElms_brch)
!  if(allocated(LeafPetoNonstElmConc_brch))deallocate(LeafPetoNonstElmConc_brch)
!  if(allocated(RootStrutElms_pft))deallocate(RootStrutElms_pft)
!  if(allocated(tCanLeafC_clyr))deallocate(tCanLeafC_clyr)
!  if(allocated(CanopyNonstElms_brch))deallocate(CanopyNonstElms_brch)
!  if(allocated(RootBiomCPerPlant_pft))deallocate(RootBiomCPerPlant_pft)
!  if(allocated(WTLF))deallocate(WTLF)
!  if(allocated(WTSHE))deallocate(WTSHE)
!  if(allocated(WTRSV))deallocate(WTRSV)
!  if(allocated(WTSTK))deallocate(WTSTK)
!  if(allocated(WTLFP))deallocate(WTLFP)
!  if(allocated(WTGR))deallocate(WTGR)
!  if(allocated(WTLFN))deallocate(WTLFN)
!  if(allocated(WTEAR))deallocate(WTEAR)
!  if(allocated(HuskStrutElms_pft))deallocate(HuskStrutElms_pft)
!  if(allocated(LeafChemElmRemob_brch))deallocate(LeafChemElmRemob_brch)
!  if(allocated(PetioleChemElmRemob_brch))deallocate(PetioleChemElmRemob_brch)
!  if(allocated(CanopyLeafSheathC_brch))deallocate(CanopyLeafSheathC_brch)
!  if(allocated(StalkRsrvElms_brch))deallocate(StalkRsrvElms_brch)
!  if(allocated(LeafStrutElms_brch))deallocate(LeafStrutElms_brch)
!  if(allocated(CanopyNodulStrutElms_brch))deallocate(CanopyNodulStrutElms_brch)
!  if(allocated(PetoleStrutElms_brch))deallocate(PetoleStrutElms_brch)
!  if(allocated(EarStrutElms_brch))deallocate(EarStrutElms_brch)
!  if(allocated(HuskStrutElms_brch))deallocate(HuskStrutElms_brch)
!  if(allocated(GrainStrutElms_brch))deallocate(GrainStrutElms_brch)
!  if(allocated(StalkStrutElms_brch))deallocate(StalkStrutElms_brch)
!  if(allocated(ShootElms_brch))deallocate(ShootElms_brch)
!  if(allocated(SenecStalkStrutElms_brch))deallocate(SenecStalkStrutElms_brch)
!  if(allocated(WTSTDI))deallocate(WTSTDI)
!  if(allocated(SeedPlantedElm_pft))deallocate(SeedPlantedElm_pft)
!  if(allocated(WTLS))deallocate(WTLS)
!  if(allocated(ShootElms_pft))deallocate(ShootElms_pft)
!  if(allocated(AvgCanopyBiomC2Graze_pft))deallocate(AvgCanopyBiomC2Graze_pft)
!  if(allocated(StandDeadStrutElms_pft)deallocate(StandDeadStrutElms_pft)
!  if(allocated(NodulStrutElms_pft))deallocate(NodulStrutElms_pft)
  end subroutine plt_biom_destroy


!----------------------------------------------------------------------

  subroutine  plt_soilchem_init(this)

  implicit none

  class(plant_soilchem_type) :: this

  allocate(this%FracBulkSOMC_vr(1:jcplx,0:JZ1));this%FracBulkSOMC_vr=spval
  allocate(this%PlantElmAllocMat4Litr(NumPlantChemElms,0:NumLitterGroups,jsken,JP1));this%PlantElmAllocMat4Litr=spval
  allocate(this%TScal4Difsvity_vr(0:JZ1));this%TScal4Difsvity_vr=spval
  allocate(this%FracAirFilledSoilPoreM_vr(60,0:JZ1));this%FracAirFilledSoilPoreM_vr=spval
  allocate(this%DiffusivitySolutEffM_vr(60,0:JZ1));this%DiffusivitySolutEffM_vr=spval
  allocate(this%VLSoilMicP_vr(0:JZ1));this%VLSoilMicP_vr=spval
  allocate(this%VLiceMicP_vr(0:JZ1));this%VLiceMicP_vr=spval
  allocate(this%VLWatMicP_vr(0:JZ1));this%VLWatMicP_vr=spval
  allocate(this%VLMicP_vr(0:JZ1));this%VLMicP_vr=spval
  allocate(this%trcs_VLN_vr(ids_nuts_beg:ids_nuts_end,0:JZ1));this%trcs_VLN_vr=spval
  allocate(this%DOM_MicP_vr(idom_beg:idom_end,1:jcplx,0:JZ1));this%DOM_MicP_vr=spval
  allocate(this%DOM_MicP_drib_vr(idom_beg:idom_end,1:jcplx,0:JZ1));this%DOM_MicP_drib_vr=spval
  allocate(this%trcs_solml_vr(ids_beg:ids_end,0:JZ1));this%trcs_solml_vr=spval
  allocate(this%trcg_gasml_vr(idg_beg:idg_NH3,0:JZ1));this%trcg_gasml_vr=spval
  allocate(this%trcg_gascl_vr(idg_beg:idg_NH3,0:JZ1));this%trcg_gascl_vr=spval
  allocate(this%CSoilOrgM_vr(1:NumPlantChemElms,0:JZ1));this%CSoilOrgM_vr=spval

  allocate(this%trc_solcl_vr(ids_beg:ids_end,0:jZ1));this%trc_solcl_vr=spval

  allocate(this%VLSoilPoreMicP_vr(0:JZ1));this%VLSoilPoreMicP_vr=spval
  allocate(this%THETW_vr(0:JZ1));this%THETW_vr=spval
  allocate(this%SoilWatAirDry_vr(0:JZ1));this%SoilWatAirDry_vr=spval

  allocate(this%GasSolbility_vr(idg_beg:idg_end,0:JZ1));this%GasSolbility_vr=spval
  allocate(this%GasDifc_vr(idg_beg:idg_end,0:JZ1));this%GasDifc_vr=spval
  allocate(this%SoluteDifusvty_vr(ids_beg:ids_end,0:JZ1));this%SoluteDifusvty_vr=spval
  allocate(this%SoilResit4RootPentrate_vr(JZ1));this%SoilResit4RootPentrate_vr=spval
  allocate(this%SoilBulkDensity_vr(0:JZ1));this%SoilBulkDensity_vr=spval
  allocate(this%HYCDMicP4RootUptake_vr(JZ1));this%HYCDMicP4RootUptake_vr=spval

  end subroutine plt_soilchem_init
!----------------------------------------------------------------------

  subroutine plt_soilchem_destroy(this)
  implicit none
  class(plant_soilchem_type) :: this

!  if(allocated(FracBulkSOMC_vr))deallocate(FracBulkSOMC_vr)

!  if(allocated(PlantElmAllocMat4Litr))deallocate(PlantElmAllocMat4Litr)
!  if(allocated(TScal4Difsvity_vr))deallocate(TScal4Difsvity_vr)
!  if(allocated(FracAirFilledSoilPoreM_vr))deallocate(FracAirFilledSoilPoreM_vr)
!  if(allocated(DiffusivitySolutEff))deallocate(DiffusivitySolutEff)
!  if(allocated(ZVSGL))deallocate(ZVSGL)
!  if(allocated(O2GSolubility))deallocate(O2GSolubility)

!  if(allocated(HYCDMicP4RootUptake_vr))deallocate(HYCDMicP4RootUptake_vr)
!  if(allocated(CGSGL))deallocate(CGSGL)
!  if(allocated(CHSGL))deallocate(CHSGL)
!  if(allocated(HGSGL))deallocate(HGSGL)
!  if(allocated(OGSGL))deallocate(OGSGL)
!  if(allocated(SoilResit4RootPentrate_vr))deallocate(SoilResit4RootPentrate_vr)

!   call destroy(this%GasSolbility_vr)
!  if(allocated(THETW))deallocate(THETW)
!  if(allocated(THETY))deallocate(THETY)
!  if(allocated(VLSoilPoreMicP_vr))deallocate(VLSoilPoreMicP_vr)
!  if(allocated(CCH4G))deallocate(CCH4G)
!  if(allocated(CZ2OG))deallocate(CZ2OG)
!  if(allocated(CNH3G))deallocate(CNH3G)
!  if(allocated(CH2GG))deallocate(CH2GG)
!  if(allocated(CORGC))deallocate(CORGC)
!  if(allocated(H2PO4))deallocate(H2PO4)
!  if(allocated(HLSGL))deallocate(HLSGL)

!   call destroy(this%trcs_VLN_vr)
!  if(allocated(O2AquaDiffusvity))deallocate(O2AquaDiffusvity)
!  if(allocated(POSGL))deallocate(POSGL)
!  if(allocated(VLSoilMicP))deallocate(VLSoilMicP)
!  if(allocated(VOLI))deallocate(VOLI)
!  if(allocated(VOLW))deallocate(VOLW)
!  if(allocated(VLMicP))deallocate(VLMicP)
!  if(allocated(ZOSGL))deallocate(ZOSGL)

!  if(allocated(OQC))deallocate(OQC)
!  if(allocated(OQN))deallocate(OQN)
!  if(allocated(OQP))deallocate(OQP)

!  if(allocated(ZNSGL))deallocate(ZNSGL)

!  if(allocated(Z2SGL))deallocate(Z2SGL)
!  if(allocated(ZHSGL))deallocate(ZHSGL)
!  if(allocated(SoilBulkDensity_vr))deallocate(SoilBulkDensity_vr)
!  if(allocated(CO2G))deallocate(CO2G)
!  if(allocated(CLSGL))deallocate(CLSGL)
!  if(allocated(CQSGL))deallocate(CQSGL)

  end subroutine plt_soilchem_destroy
!----------------------------------------------------------------------


  subroutine InitPlantAPIData()

  implicit none

  JZ1                      => pltpar%JZ1
  NumCanopyLayers1         => pltpar%NumCanopyLayers1
  JP1                      => pltpar%JP1
  NumOfLeafAzimuthSectors1 => pltpar%NumOfLeafAzimuthSectors
  NumOfSkyAzimuthSects1    => pltpar%NumOfSkyAzimuthSects1
  NumLeafZenithSectors1    => pltpar%NumLeafZenithSectors1
  MaxNodesPerBranch1       => pltpar%MaxNodesPerBranch1
  !the following variable should be consistent with the soil bgc model
  jcplx                => pltpar%jcplx
  jsken                => pltpar%jsken
  NumLitterGroups      => pltpar%NumLitterGroups
  MaxNumBranches       => pltpar%MaxNumBranches
  MaxNumRootAxes       => pltpar%MaxNumRootAxes
  NumOfPlantMorphUnits => pltpar%NumOfPlantMorphUnits
  NumOfPlantLitrCmplxs => pltpar%NumOfPlantLitrCmplxs
  NumGrowthStages      => pltpar%NumGrowthStages
  jroots               => pltpar%jroots

  call plt_site%Init()

  call plt_rbgc%Init()

  call plt_bgcr%Init()
  
  call plt_pheno%Init()

  call plt_ew%Init()

  call plt_distb%Init()

  call plt_allom%Init()

  call plt_biom%Init()

  call plt_soilchem%Init()

  call plt_rad%Init()

  call plt_photo%Init()

  call plt_morph%Init()

  call InitAllocate()
  end subroutine InitPlantAPIData
!----------------------------------------------------------------------
  subroutine InitAllocate()
  implicit none


  end subroutine InitAllocate
!----------------------------------------------------------------------
  subroutine DestructPlantAPIData
  implicit none

  call plt_pheno%Destroy()

  call plt_bgcr%Destroy()

  call plt_ew%Destroy()

  call plt_distb%Destroy()

  call plt_allom%Destroy()

  call plt_biom%Destroy()

  call plt_soilchem%Destroy()

  call plt_rad%Destroy()

  call plt_photo%Destroy()

  call plt_morph%Destroy()

  call plt_site%Destroy()



  end subroutine DestructPlantAPIData

!------------------------------------------------------------------------
  subroutine plt_rad_init(this)
! DESCRIPTION
! initialize data type for plant_radiation_type
  implicit none
  class(plant_radiation_type) :: this

  allocate(this%RadTotPAR_zsec(NumLeafZenithSectors1,NumOfSkyAzimuthSects1,NumCanopyLayers1,JP1));this%RadTotPAR_zsec=0._r8
  allocate(this%RadDifPAR_zsec(NumLeafZenithSectors1,NumOfSkyAzimuthSects1,NumCanopyLayers1,JP1));this%RadDifPAR_zsec=0._r8
  allocate(this%RadSWLeafAlbedo_pft(JP1))
  allocate(this%CanopyPARalbedo_pft(JP1))
  allocate(this%TAU_DirectRTransmit(NumCanopyLayers1+1));this%TAU_DirectRTransmit=0._r8
  allocate(this%TAU_RadThru(NumCanopyLayers1+1));this%TAU_RadThru=0._r8
  allocate(this%LWRadCanopy_pft(JP1))
  allocate(this%RadSWbyCanopy_pft(JP1))
  allocate(this%OMEGX(NumOfSkyAzimuthSects1,NumLeafZenithSectors1,NumOfLeafAzimuthSectors1));this%OMEGX=0._r8
  allocate(this%OMEGAG(NumOfSkyAzimuthSects1));this%OMEGAG=0._r8
  allocate(this%OMEGA(NumOfSkyAzimuthSects1,NumLeafZenithSectors1,NumOfLeafAzimuthSectors1));this%OMEGA=0._r8
  allocate(this%SineLeafAngle(NumLeafZenithSectors1));this%SineLeafAngle=0._r8
  allocate(this%CosineLeafAngle(NumLeafZenithSectors1));this%CosineLeafAngle=0._r8
  allocate(this%iScatteringDiffus(NumOfSkyAzimuthSects1,NumLeafZenithSectors1,NumOfLeafAzimuthSectors1))
  allocate(this%RadNet2Canopy_pft(JP1))
  allocate(this%LeafSWabsorpty_pft(JP1))
  allocate(this%LeafPARabsorpty_pft(JP1))
  allocate(this%RadPARLeafTransmis_pft(JP1))
  allocate(this%RadSWLeafTransmis_pft(JP1))
  allocate(this%RadPARbyCanopy_pft(JP1))
  allocate(this%FracPARads2Canopy_pft(JP1))
  end subroutine plt_rad_init
!------------------------------------------------------------------------
  subroutine plt_rad_destroy(this)
! DESCRIPTION
! deallocate memory for plant_radiation_type
  implicit none
  class(plant_radiation_type) :: this


  end subroutine plt_rad_destroy

!------------------------------------------------------------------------

  subroutine plt_photo_init(this)
  class(plant_photosyns_type) :: this

  allocate(this%CuticleResist_pft(JP1));this%CuticleResist_pft=spval
  allocate(this%CanopyMinStomaResistH2O_pft(JP1));this%CanopyMinStomaResistH2O_pft=spval
  allocate(this%LeafO2Solubility_pft(JP1));this%LeafO2Solubility_pft=spval
  allocate(this%CanopyBndlResist_pft(JP1));this%CanopyBndlResist_pft=spval
  allocate(this%CanPStomaResistH2O_pft(JP1));this%CanPStomaResistH2O_pft=spval
  allocate(this%DiffCO2Atmos2Intracel_pft(JP1));this%DiffCO2Atmos2Intracel_pft=spval
  allocate(this%AirConc_pft(JP1));this%AirConc_pft=spval
  allocate(this%CanopyVcMaxRubisco25C_pft(JP1));this%CanopyVcMaxRubisco25C_pft=0._r8
  allocate(this%CanopyVoMaxRubisco25C_pft(JP1));this%CanopyVoMaxRubisco25C_pft=0._r8
  allocate(this%CanopyVcMaxPEP25C_pft(JP1)); this%CanopyVcMaxPEP25C_pft=0._r8
  allocate(this%ElectronTransptJmax25C_pft(JP1));this%ElectronTransptJmax25C_pft=0._r8
  allocate(this%TFN_Carboxy_pft(JP1));this%TFN_Carboxy_pft=0._r8
  allocate(this%TFN_Oxygen_pft(JP1));this%TFN_Oxygen_pft=0._r8
  allocate(this%TFN_eTranspt_pft(JP1));this%TFN_eTranspt_pft=0._r8
  allocate(this%LeafAreaSunlit_pft(JP1)); this%LeafAreaSunlit_pft=0._r8
  allocate(this%PARSunlit_pft(JP1));this%PARSunlit_pft=0._r8
  allocate(this%PARSunsha_pft(JP1));this%PARSunsha_pft=0._r8
  allocate(this%CH2OSunlit_pft(JP1));this%CH2OSunlit_pft=0._r8
  allocate(this%CH2OSunsha_pft(JP1));this%CH2OSunsha_pft=0._r8
  allocate(this%CO2CuticleResist_pft(JP1));this%CO2CuticleResist_pft=spval
  allocate(this%LeafAreaSunlit_zsec(NumLeafZenithSectors1,NumCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafAreaSunlit_zsec=0._r8
  allocate(this%CPOOL3_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CPOOL3_node=spval
  allocate(this%CPOOL4_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CPOOL4_node=spval
  allocate(this%CMassCO2BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CMassCO2BundleSheath_node=spval
  allocate(this%CO2CompenPoint_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2CompenPoint_node=spval
  allocate(this%RubiscoCarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%RubiscoCarboxyEff_node=spval
  allocate(this%VcMaxPEPCarboxyRef_brch(MaxNumBranches,JP1));this%VcMaxPEPCarboxyRef_brch=0._r8
  allocate(this%ElectronTransptJmaxRef_brch(MaxNumBranches,JP1));this%ElectronTransptJmaxRef_brch=0._r8
  allocate(this%VcMaxRubiscoRef_brch(MaxNumBranches,JP1)); this%VcMaxRubiscoRef_brch=0._r8
  allocate(this%VoMaxRubiscoRef_brch(MaxNumBranches,JP1)); this%VoMaxRubiscoRef_brch=0._r8
  allocate(this%C4CarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%C4CarboxyEff_node=spval
  allocate(this%LigthSatCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%LigthSatCarboxyRate_node=spval
  allocate(this%LigthSatC4CarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%LigthSatC4CarboxyRate_node=spval
  allocate(this%NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%NutrientCtrlonC4Carboxy_node=spval
  allocate(this%CMassHCO3BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CMassHCO3BundleSheath_node=spval
  allocate(this%Vmax4RubiscoCarboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%Vmax4RubiscoCarboxy_node=spval
  allocate(this%ProteinCperm2LeafArea_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%ProteinCperm2LeafArea_node=spval
  allocate(this%CO2lmtRubiscoCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2lmtRubiscoCarboxyRate_node=spval
  allocate(this%Vmax4PEPCarboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%Vmax4PEPCarboxy_node=spval
  allocate(this%CO2lmtPEPCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2lmtPEPCarboxyRate_node=spval
  allocate(this%iPlantPhotosynthesisType(JP1));this%iPlantPhotosynthesisType=0
  allocate(this%Km4PEPCarboxy_pft(JP1));this%Km4PEPCarboxy_pft=spval
  allocate(this%O2L_pft(JP1));this%O2L_pft=spval
  allocate(this%CO2Solubility_pft(JP1));this%CO2Solubility_pft=spval
  allocate(this%CanopyGasCO2_pft(JP1));this%CanopyGasCO2_pft=spval
  allocate(this%ChillHours_pft(JP1));this%ChillHours_pft=spval
  allocate(this%SpecLeafChlAct_pft(JP1));this%SpecLeafChlAct_pft=spval
  allocate(this%LeafC3Chl2Protein_pft(JP1));this%LeafC3Chl2Protein_pft=spval
  allocate(this%LeafPEP2Protein_pft(JP1));this%LeafPEP2Protein_pft=spval
  allocate(this%LeafC4Chl2Protein_pft(JP1));this%LeafC4Chl2Protein_pft=spval
  allocate(this%LeafRubisco2Protein_pft(JP1));this%LeafRubisco2Protein_pft=spval
  allocate(this%VmaxPEPCarboxyRef_pft(JP1));this%VmaxPEPCarboxyRef_pft=spval
  allocate(this%VmaxRubOxyRef_pft(JP1));this%VmaxRubOxyRef_pft=spval
  allocate(this%VmaxSpecRubCarboxyRef_pft(JP1));this%VmaxSpecRubCarboxyRef_pft=spval
  allocate(this%XKCO2_pft(JP1));this%XKCO2_pft=spval
  allocate(this%XKO2_pft(JP1));this%XKO2_pft=spval
  allocate(this%RubiscoActivity_brch(MaxNumBranches,JP1));this%RubiscoActivity_brch=1._r8
  allocate(this%C4PhotosynDowreg_brch(MaxNumBranches,JP1));this%C4PhotosynDowreg_brch=spval
  allocate(this%aquCO2Intraleaf_pft(JP1));this%aquCO2Intraleaf_pft=spval
  allocate(this%Km4LeafaqCO2_pft(JP1));this%Km4LeafaqCO2_pft=spval
  allocate(this%Km4RubiscoCarboxy_pft(JP1));this%Km4RubiscoCarboxy_pft=spval
  allocate(this%LeafIntracellularCO2_pft(JP1));this%LeafIntracellularCO2_pft=spval
  allocate(this%O2I_pft(JP1));this%O2I_pft=spval
  allocate(this%RCS_pft(JP1));this%RCS_pft=spval
  allocate(this%CanopyCi2CaRatio_pft(JP1));this%CanopyCi2CaRatio_pft=spval
  allocate(this%H2OCuticleResist_pft(JP1));this%H2OCuticleResist_pft=spval

  end subroutine plt_photo_init
!------------------------------------------------------------------------

  subroutine plt_photo_destroy(this)
  class(plant_photosyns_type) :: this


  end subroutine plt_photo_destroy

!------------------------------------------------------------------------
  subroutine plt_pheno_init(this)
  implicit none
  class(plant_pheno_type) :: this


  allocate(this%TCChill4Seed_pft(JP1));this%TCChill4Seed_pft=spval
  allocate(this%TempOffset_pft(JP1));this%TempOffset_pft=spval
  allocate(this%MatureGroup_pft(JP1));this%MatureGroup_pft=spval
  allocate(this%PlantO2Stress_pft(JP1));this%PlantO2Stress_pft=spval
  allocate(this%NonstCMinConc2InitBranch_pft(JP1));this%NonstCMinConc2InitBranch_pft=spval
  allocate(this%MinNonstC2InitRoot_pft(JP1));this%MinNonstC2InitRoot_pft=spval
  allocate(this%fTCanopyGroth_pft(JP1));this%fTCanopyGroth_pft=spval
  allocate(this%TC4LeafOut_pft(JP1));this%TC4LeafOut_pft=spval
  allocate(this%TCGroth_pft(JP1));this%TCGroth_pft=spval
  allocate(this%TKGroth_pft(JP1));this%TKGroth_pft=spval
  allocate(this%TC4LeafOff_pft(JP1));this%TC4LeafOff_pft=spval
  allocate(this%HoursTooLowPsiCan_pft(JP1));this%HoursTooLowPsiCan_pft=spval
  allocate(this%LeafElmntRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafElmntRemobFlx_brch=spval
  allocate(this%PetioleChemElmRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetioleChemElmRemobFlx_brch=spval

  allocate(this%fTgrowRootP_vr(JZ1,JP1));this%fTgrowRootP_vr=spval
  allocate(this%GrainFillRate25C_pft(JP1));this%GrainFillRate25C_pft=spval
  allocate(this%ShootRootNonstElmConduts_pft(JP1));this%ShootRootNonstElmConduts_pft=spval
  allocate(this%SeedTempSens_pft(JP1));this%SeedTempSens_pft=spval
  allocate(this%HighTempLimitSeed_pft(JP1));this%HighTempLimitSeed_pft=spval
  allocate(this%PlantInitThermoAdaptZone_pft(JP1));this%PlantInitThermoAdaptZone_pft=spval
  allocate(this%rPlantThermoAdaptZone_pft(JP1));this%rPlantThermoAdaptZone_pft=0
  allocate(this%IsPlantActive_pft(JP1));this%IsPlantActive_pft=0
  allocate(this%iPlantState_pft(JP1));this%iPlantState_pft=0
  allocate(this%NetCumElmntFlx2Plant_pft(NumPlantChemElms,JP1));this%NetCumElmntFlx2Plant_pft=spval
  allocate(this%MatureGroup_brch(MaxNumBranches,JP1));this%MatureGroup_brch=spval
  allocate(this%iPlantShootState_pft(JP1));this%iPlantShootState_pft=0
  allocate(this%Hours2LeafOut_brch(MaxNumBranches,JP1));this%Hours2LeafOut_brch=spval
  allocate(this%HoursDoingRemob_brch(MaxNumBranches,JP1));this%HoursDoingRemob_brch=spval
  allocate(this%HourlyNodeNumNormByMatgrp_brch(MaxNumBranches,JP1));this%HourlyNodeNumNormByMatgrp_brch=spval
  allocate(this%dReproNodeNumNormByMatG_brch(MaxNumBranches,JP1));this%dReproNodeNumNormByMatG_brch=spval
  allocate(this%NodeNumNormByMatgrp_brch(MaxNumBranches,JP1));this%NodeNumNormByMatgrp_brch=spval
  allocate(this%ReprodNodeNumNormByMatrgrp_brch(MaxNumBranches,JP1));this%ReprodNodeNumNormByMatrgrp_brch=spval
  allocate(this%HourFailGrainFill_brch(MaxNumBranches,JP1));this%HourFailGrainFill_brch=spval
  allocate(this%iPlantPhenolType_pft(JP1));this%iPlantPhenolType_pft=0
  allocate(this%iPlantPhenolPattern_pft(JP1));this%iPlantPhenolPattern_pft=0
  allocate(this%iPlantTurnoverPattern_pft(JP1));this%iPlantTurnoverPattern_pft=0
  allocate(this%iPlantRootState_pft(JP1));this%iPlantRootState_pft=0
  allocate(this%iPlantDevelopPattern_pft(JP1));this%iPlantDevelopPattern_pft=0
  allocate(this%iPlantPhotoperiodType_pft(JP1));this%iPlantPhotoperiodType_pft=0
  allocate(this%doInitPlant_pft(JP1));this%doInitPlant_pft=0
  allocate(this%doReSeed_pft(JP1)); this%doReSeed_pft=.false.    
  allocate(this%iPlantRootProfile_pft(JP1));this%iPlantRootProfile_pft=0
  allocate(this%KHiestGroLeafNode_brch(MaxNumBranches,JP1));this%KHiestGroLeafNode_brch=0
  allocate(this%KLowestGroLeafNode_brch(MaxNumBranches,JP1));this%KLowestGroLeafNode_brch=0
  allocate(this%fRootGrowPSISense_pvr(jroots,JZ1,JP1)); this%fRootGrowPSISense_pvr=spval
  allocate(this%iPlantCalendar_brch(NumGrowthStages,MaxNumBranches,JP1));this%iPlantCalendar_brch=0
  allocate(this%TotalNodeNumNormByMatgrp_brch(MaxNumBranches,JP1));this%TotalNodeNumNormByMatgrp_brch=spval
  allocate(this%TotReproNodeNumNormByMatrgrp_brch(MaxNumBranches,JP1));this%TotReproNodeNumNormByMatrgrp_brch=spval
  allocate(this%LeafNumberAtFloralInit_brch(MaxNumBranches,JP1));this%LeafNumberAtFloralInit_brch=spval
  allocate(this%RefLeafAppearRate_pft(JP1));this%RefLeafAppearRate_pft=spval
  allocate(this%RefNodeInitRate_pft(JP1));this%RefNodeInitRate_pft=spval
  allocate(this%CriticPhotoPeriod_pft(JP1));this%CriticPhotoPeriod_pft=spval
  allocate(this%PhotoPeriodSens_pft(JP1));this%PhotoPeriodSens_pft=spval
  allocate(this%iPlantBranchState_brch(MaxNumBranches,JP1));this%iPlantBranchState_brch=0
  allocate(this%doRemobilization_brch(MaxNumBranches,JP1));this%doRemobilization_brch=0
  allocate(this%doPlantLeaveOff_brch(MaxNumBranches,JP1));this%doPlantLeaveOff_brch=0
  allocate(this%doPlantLeafOut_brch(MaxNumBranches,JP1));this%doPlantLeafOut_brch=0
  allocate(this%doInitLeafOut_brch(MaxNumBranches,JP1));this%doInitLeafOut_brch=0
  allocate(this%doSenescence_brch(MaxNumBranches,JP1));this%doSenescence_brch=0
  allocate(this%Prep4Literfall_brch(MaxNumBranches,JP1));this%Prep4Literfall_brch=0
  allocate(this%Hours4LiterfalAftMature_brch(MaxNumBranches,JP1));this%Hours4LiterfalAftMature_brch=0
  allocate(this%Hours4LenthenPhotoPeriod_brch(MaxNumBranches,JP1));this%Hours4LenthenPhotoPeriod_brch=spval
  allocate(this%Hours4ShortenPhotoPeriod_brch(MaxNumBranches,JP1));this%Hours4ShortenPhotoPeriod_brch=spval
  allocate(this%Hours4Leafout_brch(MaxNumBranches,JP1));this%Hours4Leafout_brch=spval
  allocate(this%HourReq4LeafOut_brch(NumCanopyLayers1,JP1));this%HourReq4LeafOut_brch=spval
  allocate(this%Hours4LeafOff_brch(MaxNumBranches,JP1));this%Hours4LeafOff_brch=spval
  allocate(this%HourReq4LeafOff_brch(NumCanopyLayers1,JP1));this%HourReq4LeafOff_brch=spval

  end subroutine plt_pheno_init
!------------------------------------------------------------------------

  subroutine plt_pheno_destroy(this)
  implicit none
  class(plant_pheno_type) :: this


  end subroutine plt_pheno_destroy
!------------------------------------------------------------------------

  subroutine plt_morph_init(this)
  implicit none
  class(plant_morph_type) :: this

  allocate(this%RootAreaPerPlant_pvr(jroots,JZ1,JP1));this%RootAreaPerPlant_pvr=spval
  allocate(this%RootLenDensPerPlant_pvr(jroots,JZ1,JP1));this%RootLenDensPerPlant_pvr=spval
  allocate(this%RootArea1stPP_pvr(jroots,JZ1,JP1)); this%RootArea1stPP_pvr=0._r8
  allocate(this%RootArea2ndPP_pvr(jroots,JZ1,JP1)); this%RootArea2ndPP_pvr=0._r8
  allocate(this%RootPoreVol_rpvr(jroots,JZ1,JP1));this%RootPoreVol_rpvr=spval
  allocate(this%RootVH2O_pvr(jroots,JZ1,JP1));this%RootVH2O_pvr=spval
  allocate(this%Root1stXNumL_rpvr(jroots,JZ1,JP1));this%Root1stXNumL_rpvr=spval
  allocate(this%Root2ndXNumL_rpvr(jroots,JZ1,JP1));this%Root2ndXNumL_rpvr=spval
  allocate(this%SeedCMass_pft(JP1));this%SeedCMass_pft=spval
  allocate(this%totRootLenDens_vr(JZ1));this%totRootLenDens_vr=spval
  allocate(this%RootBranchFreq_pft(JP1));this%RootBranchFreq_pft=spval
  allocate(this%ClumpFactorInit_pft(JP1));this%ClumpFactorInit_pft=spval
  allocate(this%ClumpFactorNow_pft(JP1));this%ClumpFactorNow_pft=spval
  allocate(this%HypoctoHeight_pft(JP1));this%HypoctoHeight_pft=spval
  allocate(this%RootPoreTortu4Gas_pft(jroots,JP1));this%RootPoreTortu4Gas_pft=spval
  allocate(this%rLen2WidthLeaf_pft(JP1));this%rLen2WidthLeaf_pft=spval
  allocate(this%MaxSeedNumPerSite_pft(JP1));this%MaxSeedNumPerSite_pft=spval
  allocate(this%GrothStalkMaxSeedSites_pft(JP1));this%GrothStalkMaxSeedSites_pft=spval
  allocate(this%SeedCMassMax_pft(JP1));this%SeedCMassMax_pft=spval
  allocate(this%Root1stMaxRadius1_pft(jroots,JP1));this%Root1stMaxRadius1_pft=spval
  allocate(this%Root2ndMaxRadius1_pft(jroots,JP1));this%Root2ndMaxRadius1_pft=spval
  allocate(this%RootRaidus_rpft(jroots,JP1));this%RootRaidus_rpft=spval
  allocate(this%Root1stRadius_pvr(jroots,JZ1,JP1));this%Root1stRadius_pvr=spval
  allocate(this%Root2ndRadius_rpvr(jroots,JZ1,JP1));this%Root2ndRadius_rpvr=spval
  allocate(this%Root1stMaxRadius_pft(jroots,JP1));this%Root1stMaxRadius_pft=spval
  allocate(this%Root2ndMaxRadius_pft(jroots,JP1));this%Root2ndMaxRadius_pft=spval

  allocate(this%Root1stDepz_pft(jroots,MaxNumRootAxes,JP1));this%Root1stDepz_pft=spval
  allocate(this%RootTotLenPerPlant_pvr(jroots,JZ1,JP1));this%RootTotLenPerPlant_pvr=spval
  allocate(this%RootLenPerPlant_pvr(jroots,JZ1,JP1));this%RootLenPerPlant_pvr=0._r8
  allocate(this%Root2ndEffLen4uptk_rpvr(jroots,JZ1,JP1));this%Root2ndEffLen4uptk_rpvr=spval
  allocate(this%Root1stSpecLen_pft(jroots,JP1));this%Root1stSpecLen_pft=spval
  allocate(this%Root2ndSpecLen_pft(jroots,JP1));this%Root2ndSpecLen_pft=spval
  allocate(this%Root1stLen_rpvr(jroots,JZ1,MaxNumRootAxes,JP1));this%Root1stLen_rpvr=spval
  allocate(this%Root2ndLen_rpvr(jroots,JZ1,MaxNumRootAxes,JP1));this%Root2ndLen_rpvr=spval
  allocate(this%Root2ndXNum_rpvr(jroots,JZ1,MaxNumRootAxes,JP1));this%Root2ndXNum_rpvr=spval
  allocate(this%iPlantNfixType_pft(JP1));this%iPlantNfixType_pft=0
  allocate(this%Myco_pft(JP1));this%Myco_pft=0
  allocate(this%CanopyHeight4WatUptake_pft(JP1));this%CanopyHeight4WatUptake_pft=spval
  allocate(this%KLeafNumber_brch(MaxNumBranches,JP1));this%KLeafNumber_brch=0
  allocate(this%NumOfLeaves_brch(MaxNumBranches,JP1));this%NumOfLeaves_brch=spval
  allocate(this%NGTopRootLayer_pft(JP1));this%NGTopRootLayer_pft=0;
  allocate(this%CanopyHeight_pft(JP1));this%CanopyHeight_pft=spval
  allocate(this%ShootNodeNumAtPlanting_pft(JP1));this%ShootNodeNumAtPlanting_pft=spval
  allocate(this%CanopyHeightZ_col(0:NumCanopyLayers1));this%CanopyHeightZ_col=spval
  allocate(this%CanopyStemArea_pft(JP1));this%CanopyStemArea_pft=spval
  allocate(this%CanopyLeafArea_pft(JP1));this%CanopyLeafArea_pft=spval
  allocate(this%MainBranchNum_pft(JP1));this%MainBranchNum_pft=0
  allocate(this%NMaxRootBotLayer_pft(JP1));this%NMaxRootBotLayer_pft=0
  allocate(this%NumPrimeRootAxes_pft(JP1));this%NumPrimeRootAxes_pft=0
  allocate(this%NumCogrowthNode_pft(JP1));this%NumCogrowthNode_pft=0
  allocate(this%BranchNumber_pft(JP1));this%BranchNumber_pft=0
  allocate(this%NumOfBranches_pft(JP1));this%NumOfBranches_pft=0
  allocate(this%NIXBotRootLayer_raxes(MaxNumRootAxes,JP1));this%NIXBotRootLayer_raxes=0
  allocate(this%PARTS_brch(NumOfPlantMorphUnits,MaxNumBranches,JP1));this%PARTS_brch=spval
  allocate(this%tlai_day_pft(JP1)); this%tlai_day_pft=spval
  allocate(this%tsai_day_pft(JP1)); this%tsai_day_pft=spval
  allocate(this%ShootNodeNum_brch(MaxNumBranches,JP1));this%ShootNodeNum_brch=spval
  allocate(this%NodeNum2InitFloral_brch(MaxNumBranches,JP1));this%NodeNum2InitFloral_brch=spval
  allocate(this%NodeNumberAtAnthesis_brch(MaxNumBranches,JP1));this%NodeNumberAtAnthesis_brch=spval
  allocate(this%SineBranchAngle_pft(JP1));this%SineBranchAngle_pft=spval
  allocate(this%LeafAreaZsec_brch(NumLeafZenithSectors1,NumCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafAreaZsec_brch=spval
  allocate(this%KMinNumLeaf4GroAlloc_brch(MaxNumBranches,JP1));this%KMinNumLeaf4GroAlloc_brch=0
  allocate(this%BranchNumerID_brch(MaxNumBranches,JP1));this%BranchNumerID_brch=0
  allocate(this%PotentialSeedSites_brch(MaxNumBranches,JP1));this%PotentialSeedSites_brch=spval
  allocate(this%CanPBranchHeight(MaxNumBranches,JP1));this%CanPBranchHeight=spval
  allocate(this%LeafAreaDying_brch(MaxNumBranches,JP1));this%LeafAreaDying_brch=spval
  allocate(this%LeafAreaLive_brch(MaxNumBranches,JP1));this%LeafAreaLive_brch=spval
  allocate(this%SinePetioleAngle_pft(JP1));this%SinePetioleAngle_pft=spval
  allocate(this%LeafAngleClass_pft(NumLeafZenithSectors1,JP1));this%LeafAngleClass_pft=spval
  allocate(this%CanopyStemAreaZ_pft(NumCanopyLayers1,JP1));this%CanopyStemAreaZ_pft=spval
  allocate(this%CanopyLeafAreaZ_pft(NumCanopyLayers1,JP1));this%CanopyLeafAreaZ_pft=spval
  allocate(this%LeafArea_node(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafArea_node=spval
  allocate(this%CanopySeedNum_pft(JP1));this%CanopySeedNum_pft=spval
  allocate(this%SeedDepth_pft(JP1));this%SeedDepth_pft=spval
  allocate(this%PlantinDepz_pft(JP1));this%PlantinDepz_pft=spval
  allocate(this%SeedMeanLen_pft(JP1));this%SeedMeanLen_pft=spval
  allocate(this%SeedVolumeMean_pft(JP1));this%SeedVolumeMean_pft=spval
  allocate(this%SeedAreaMean_pft(JP1));this%SeedAreaMean_pft=spval
  allocate(this%CanopyStemAareZ_col(NumCanopyLayers1));this%CanopyStemAareZ_col=spval
  allocate(this%PetoLen2Mass_pft(JP1));this%PetoLen2Mass_pft=spval
  allocate(this%NodeLenPergC_pft(JP1));this%NodeLenPergC_pft=spval
  allocate(this%SLA1_pft(JP1));this%SLA1_pft=spval
  allocate(this%CanopyLeafAareZ_col(NumCanopyLayers1));this%CanopyLeafAareZ_col=spval
  allocate(this%LeafStalkArea_pft(JP1));this%LeafStalkArea_pft=spval
  allocate(this%StalkNodeVertLength_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%StalkNodeVertLength_brch=spval
  allocate(this%PetoleLensNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%PetoleLensNode_brch=spval
  allocate(this%StalkNodeHeight_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%StalkNodeHeight_brch=spval
  allocate(this%StemAreaZsec_brch(NumLeafZenithSectors1,NumCanopyLayers1,MaxNumBranches,JP1));this%StemAreaZsec_brch=0._r8
  allocate(this%CanopyLeafArea_lnode(NumCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%CanopyLeafArea_lnode=0._r8
  allocate(this%CanopyStalkArea_lbrch(NumCanopyLayers1,MaxNumBranches,JP1));this%CanopyStalkArea_lbrch=spval
  allocate(this%TreeRingAveRadius_pft(JP1));this%TreeRingAveRadius_pft=spval
  allocate(this%MaxSoiL4Root_pft(JP1));this%MaxSoiL4Root_pft=0
  allocate(this%SeedSitesSet_brch(MaxNumBranches,JP1));this%SeedSitesSet_brch=spval
  allocate(this%ClumpFactor_pft(JP1));this%ClumpFactor_pft=spval
  allocate(this%RootVolPerMassC_pft(jroots,JP1));this%RootVolPerMassC_pft=spval
  allocate(this%RootPorosity_pft(jroots,JP1));this%RootPorosity_pft=spval
  allocate(this%Root2ndXSecArea_pft(jroots,JP1));this%Root2ndXSecArea_pft=spval
  allocate(this%Root1stXSecArea_pft(jroots,JP1));this%Root1stXSecArea_pft=spval
  allocate(this%RootRadialResist_pft(jroots,JP1));this%RootRadialResist_pft=spval
  allocate(this%RootAxialResist_pft(jroots,JP1));this%RootAxialResist_pft=spval
  allocate(this%iPlantGrainType_pft(JP1));this%iPlantGrainType_pft=0
  end subroutine plt_morph_init

!------------------------------------------------------------------------

  subroutine plt_morph_destroy(this)
  implicit none
  class(plant_morph_type) :: this

  end subroutine plt_morph_destroy
end module PlantAPIData
