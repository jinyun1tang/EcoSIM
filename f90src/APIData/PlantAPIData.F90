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

  integer, pointer :: NumGrowthStages       !number of growth stages
  integer, pointer :: MaxNumRootAxes          !maximum number of root layers
  integer, pointer :: MaxNumBranches         !maximum number of plant branches
  integer, pointer :: JP1         !number of plants
  integer, pointer :: NumOfSkyAzimuSects1        !number of sectors for the sky azimuth  [0,2*pi]
  integer, pointer :: jcplx       !number of organo-microbial complexes
  integer, pointer :: NumOfLeafAzimuthSectors1        !number of sectors for the leaf azimuth, [0,pi]
  integer, pointer :: NumOfCanopyLayers1         !number of canopy layers
  integer, pointer :: JZ1         !number of soil layers
  integer, pointer :: NumOfLeafZenithSectors1        !number of sectors for the leaf zenith [0,pi/2]
  integer, pointer :: MaxNodesPerBranch1      !number of canopy nodes
  integer, pointer :: jsken       !number of kinetic components in litter,  PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  integer, pointer :: NumLitterGroups     !number of litter groups nonstructural(0,*),
                          !     foliar(1,*),non-foliar(jroots,*),stalk(3,*),root(4,*), coarse woody (5,*)
  integer, pointer :: NumOfPlantMorphUnits        !number of organs involved in partition
  integer, pointer :: NumOfPlantLitrCmplxs
  integer, pointer :: jroots
!begin_data

  type, public :: plant_siteinfo_type
  real(r8) :: ALAT      !latitude	[degrees]
  real(r8) :: ATCA      !mean annual air temperature, [oC]
  real(r8) :: ALT       !altitude of grid cell, [m]
  real(r8) :: CCO2EI    !initial atmospheric CO2 concentration, [g m-3]
  real(r8) :: CO2EI     !initial atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: COXYE     !atmospheric O2 concentration, [g m-3]
  real(r8) :: SolarNoonHour_col    !time of solar noon, [h]
  real(r8) :: CO2E      !atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: DayLenthPrev      !daylength of previous day, [h]
  real(r8) :: DayLenthCurrent      !daylength, [h]
  real(r8) :: DayLenthMax      !maximum daylength, [h]
  real(r8) :: SoilSurfRoughnesst0_col        !initial soil surface roughness height, [m]
  real(r8) :: OXYE      !atmospheric O2 concentration, [umol mol-1]
  real(r8) :: PlantPopu_col       !total plant population, [d-2]
  real(r8) :: POROS1    !top layer soil porosity
  real(r8) :: WindSpeedAtm_col        !wind speed, [m h-1]
  real(r8), pointer :: PlantElemntStoreLandscape(:)  => null() !total plant element balance	g d-2
  real(r8) :: WindMesHeight        !wind speed measurement height, [m]
  real(r8) :: ZEROS2    !threshold zero
  real(r8) :: ZEROS     !threshold zero
  real(r8) :: ZERO      !threshold zero
  real(r8) :: ZERO2     !threshold zero
  integer :: KoppenClimZone     !Koppen climate zone
  integer :: NumActivePlants      !number of active PFT
  integer :: NL         !lowest soil layer number
  integer :: NP0        !intitial number of plant species
  integer :: MaxNumRootLays         !maximum root layer number
  integer :: NP         !number of plant species
  integer :: NU         !soil surface layer number
  integer :: DazCurrYear       !number of days in current year
  integer :: iYearCurrent       !current year

  real(r8) :: QH2OLoss_lnds    !total subsurface water flux	m3 d-2
  real(r8), pointer :: AtmGasc(:)  => null()  !atmospheric gas concentrations
  character(len=16), pointer :: DATAP(:) => null()   !parameter file name
  CHARACTER(len=16), pointer :: DATA(:)  => null()   !pft file
  real(r8), pointer :: AREA3(:)   => null()    !soil cross section area (vertical plan defined by its normal direction)
  real(r8), pointer :: ElmBalanceCum_pft(:,:)      => null()    !plant element balance, [g d-2]
  real(r8), pointer :: CumSoilThickness_vr(:)         => null()    !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8), pointer :: PPI_pft(:)                      => null()    !initial plant population, [m-2]
  real(r8), pointer :: PPatSeeding_pft(:)                      => null()    !plant population at seeding, [m-2]
  real(r8), pointer :: PPX_pft(:)                      => null()    !plant population, [m-2]
  logical , pointer :: flag_pft_active(:)          => null()
  real(r8), pointer :: PlantPopulation_pft(:)      => null()    !plant population, [d-2]
  real(r8), pointer :: DPTHZ_vr(:)   => null()    !depth to middle of soil layer from  surface of grid cell [m]
  real(r8), pointer :: FracSoiAsMicP_vr(:)    => null()    !micropore fraction
  real(r8), pointer :: DLYR3(:)   => null()    !vertical thickness of soil layer [m]
  real(r8), pointer :: VLWatMicPM_vr(:,:) => null()    !soil micropore water content, [m3 d-2]
  real(r8), pointer :: VLsoiAirPM(:,:) => null()    !soil air content, [m3 d-2]
  real(r8), pointer :: TortMicPM_vr(:,:)  => null()    !micropore soil tortuosity, []
  real(r8), pointer :: FILM(:,:)  => null()    !soil water film thickness , [m]
  contains
    procedure, public :: Init =>  plt_site_Init
    procedure, public :: Destroy => plt_site_destroy
  end type plant_siteinfo_type

  type, public :: plant_photosyns_type
  real(r8), pointer :: SpecChloryfilAct_pft(:)              => null()  !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: LeafC3ChlorofilConc_pft(:)           => null()  !leaf C3 chlorophyll content, [gC gC-1]
  real(r8), pointer :: FracLeafProtinAsPEPCarboxyl_pft(:)   => null()  !leaf PEP carboxylase content, [gC gC-1]
  real(r8), pointer :: LeafC4ChlorofilConc_pft(:)           => null()  !leaf C4 chlorophyll content, [gC gC-1]
  real(r8), pointer :: LeafRuBPConc_pft(:)                  => null()  !leaf rubisco content, [gC gC-1]
  real(r8), pointer :: VmaxPEPCarboxyRef_pft(:)  => null()  !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: VmaxRubOxyRef_pft(:)   => null()  !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: VmaxRubCarboxyRef_pft(:)   => null()  !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: XKCO2(:)  => null()  !Km for rubisco carboxylase activity, [uM]
  real(r8), pointer :: XKO2(:)   => null()  !Km for rubisco oxygenase activity, [uM]
  real(r8), pointer :: RubiscoActivity_brch(:,:) => null()   !branch down-regulation of CO2 fixation, [-]
  real(r8), pointer :: C4PhotosynDowreg_brch(:,:)=> null()   !down-regulation of C4 photosynthesis, [-]
  real(r8), pointer :: aquCO2Intraleaf_pft(:)   => null()   !leaf aqueous CO2 concentration, [uM]
  real(r8), pointer :: O2I(:)    => null()   !leaf gaseous O2 concentration, [umol m-3]
  real(r8), pointer :: LeafIntracellularCO2_pft(:)   => null()   !leaf gaseous CO2 concentration, [umol m-3]
  real(r8), pointer :: Km4RubiscoCarboxy_pft(:) => null()   !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8), pointer :: Km4LeafaqCO2_pft(:) => null()   !leaf aqueous CO2 Km no O2, [uM]
  real(r8), pointer :: RCS(:)    => null()   !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8), pointer :: CanPCi2CaRatio(:)   => null()   !Ci:Ca ratio, [-]
  real(r8), pointer :: H2OCuticleResist_pft(:)   => null()   !maximum stomatal resistance to vapor, [s h-1]
  real(r8), pointer :: ChillHours_pft(:)  => null()   !chilling effect on CO2 fixation, [-]
  real(r8), pointer :: CO2Solubility_pft(:)   => null()   !leaf CO2 solubility, [uM /umol mol-1]
  real(r8), pointer :: CanopyGasCO2_pft(:)   => null()   !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8), pointer :: O2L(:)    => null()   !leaf aqueous O2 concentration, [uM]
  real(r8), pointer :: Km4PEPCarboxy_pft(:) => null()   !Km for PEP carboxylase activity, [uM]
  integer,  pointer :: iPlantPhotosynthesisType(:)  => null()   !plant photosynthetic type (C3 or C4)
  real(r8), pointer :: MinCanPStomaResistH2O_pft(:)   => null()   !canopy minimum stomatal resistance, [s m-1]
  real(r8), pointer :: CuticleResist_pft(:)   => null()   !maximum stomatal resistance to vapor, [s m-1]
  real(r8), pointer :: CanPStomaResistH2O_pft(:)     => null()   !canopy stomatal resistance, [h m-1]
  real(r8), pointer :: CanopyBndlResist_pft(:)     => null()   !canopy boundary layer resistance, [h m-1]
  real(r8), pointer :: LeafO2Solubility_pft(:)    => null()   !leaf O2 solubility, [uM /umol mol-1]
  real(r8), pointer :: CO2CuticleResist_pft(:)   => null()   !maximum stomatal resistance to CO2, [s h-1]
  real(r8), pointer :: DiffCO2Atmos2Intracel_pft(:)   => null()   !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8), pointer :: LeafAUnshaded_zsec(:,:,:,:,:)  => null()   !leaf irradiated surface area, [m2 d-2]
  real(r8), pointer :: CPOOL3_node(:,:,:)     => null()   !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8), pointer :: CPOOL4_node(:,:,:)     => null()   !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: CMassCO2BundleSheath_node(:,:,:)       => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: CO2CompenPoint_node(:,:,:)      => null()   !CO2 compensation point, [uM]
  real(r8), pointer :: RubiscoCarboxyEff_node(:,:,:)       => null()   !carboxylation efficiency, [umol umol-1]
  real(r8), pointer :: C4CarboxyEff_node(:,:,:)      => null()   !C4 carboxylation efficiency, [umol umol-1]
  real(r8), pointer :: LigthSatCarboxyRate_node(:,:,:)      => null()   !maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: LigthSatC4CarboxyRate_node(:,:,:)      => null()   !maximum  light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: NutrientCtrlonC4Carboxy_node(:,:,:)      => null()   !down-regulation of C4 photosynthesis, [-]
  real(r8), pointer :: CMassHCO3BundleSheath_node(:,:,:)       => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: Vmax4RubiscoCarboxy_pft(:,:,:)      => null()   !maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: CO2lmtRubiscoCarboxyRate_node(:,:,:)       => null()   !carboxylation rate, [umol m-2 s-1]
  real(r8), pointer :: Vmax4PEPCarboxy_pft(:,:,:)      => null()   !maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: CO2lmtPEPCarboxyRate_node(:,:,:)      => null()   !C4 carboxylation rate, [umol m-2 s-1]
  real(r8), pointer :: AirConc_pft(:)           => null()   !total gas concentration, [mol m-3]

  contains
    procedure, public :: Init    =>  plt_photo_init
    procedure, public :: Destroy => plt_photo_destroy
  end type plant_photosyns_type

  type, public ::  plant_radiation_type
  real(r8) :: TotSineSkyAngles_grd     !sine of sky angles
  real(r8) :: SoilAlbedo      !soil albedo
  real(r8) :: SurfAlbedo_col      !Surface albedo
  real(r8) :: RadPARSolarBeam_col      !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8) :: RadSWDiffus_col     !diffuse shortwave radiation, [W m-2]
  real(r8) :: RadSWDirect_col      !direct shortwave radiation, [W m-2]
  real(r8) :: RadPARDiffus_col      !diffuse PAR, [umol m-2 s-1]
  real(r8) :: RadPARDirect_col     !direct PAR, [umol m-2 s-1]
  real(r8) :: Eco_NetRad_col       !ecosystem net radiation, [MJ d-2 h-1]
  real(r8) :: RadSWSolarBeam_col      !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8) :: FracSWRad2Grnd_col     !fraction of radiation intercepted by ground surface, [-]
  real(r8) :: RadSWGrnd_col      !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8) :: SineGrndSlope_col      !sine of slope, [-]
  real(r8) :: GroundSurfAzimuth_col     !azimuth of slope, [-]
  real(r8) :: CosineGrndSlope_col      !cosine of slope, [-]
  real(r8) :: LWRadGrnd    !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8) :: LWRadSky       !sky longwave radiation , [MJ d-2 h-1]
  real(r8) :: SineSunInclAnglNxtHour_col     !sine of solar angle next hour, [-]
  real(r8) :: SineSunInclAngle_col      !sine of solar angle, [-]
  real(r8), pointer :: RadSWLeafAlbedo_pft(:)     => null() !canopy shortwave albedo , [-]
  real(r8), pointer :: CanopyPARalbedo_pft(:)     => null() !canopy PAR albedo , [-]
  real(r8), pointer :: TAU_DirRadTransm(:)     => null() !fraction of radiation intercepted by canopy layer, [-]
  real(r8), pointer :: TAU_RadThru(:)     => null() !fraction of radiation transmitted by canopy layer, [-]
  real(r8), pointer :: LWRadCanopy_pft(:)    => null() !canopy longwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: RadSWbyCanopy_pft(:)     => null() !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: OMEGX(:,:,:)=> null() !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation
  real(r8), pointer :: OMEGAG(:)   => null() !sine of solar beam on leaf surface, [-]
  real(r8), pointer :: OMEGA(:,:,:)=> null() !sine of indirect sky radiation on leaf surface
  real(r8), pointer :: SineLeafAngle(:)     => null() !sine of leaf angle
  real(r8), pointer :: CosineLeafAngle(:)     => null() !cosine of leaf angle
  integer,  pointer :: iScatteringDiffus(:,:,:)=> null() !flag for calculating backscattering of radiation in canopy
  real(r8), pointer :: RadNet2Canopy_pft(:)     => null() !canopy net radiation , [MJ d-2 h-1]
  real(r8), pointer :: LeafSWabsorpty_pft(:)     => null() !canopy shortwave absorptivity , [-]
  real(r8), pointer :: LeafPARabsorpty_pft(:)     => null() !canopy PAR absorptivity
  real(r8), pointer :: RadSWLeafTransmis_pft(:)     => null() !canopy shortwave transmissivity , [-]
  real(r8), pointer :: RadPARLeafTransmis_pft(:)     => null() !canopy PAR transmissivity , [-]
  real(r8), pointer :: RadPARbyCanopy_pft(:)     => null() !canopy absorbed PAR , [umol m-2 s-1]
  real(r8), pointer :: FracPARads2Canopy_pft(:)    => null() !fraction of incoming PAR absorbed by canopy, [-]
  real(r8), pointer :: RadPAR_zsec(:,:,:,:)=> null()     !direct incoming PAR, [umol m-2 s-1]
  real(r8), pointer :: RadDifPAR_zsec(:,:,:,:)=> null()  !diffuse incoming PAR, [umol m-2 s-1]
  contains
    procedure, public :: Init    => plt_rad_init
    procedure, public :: Destroy => plt_rad_destroy
  end type plant_radiation_type

  type, public :: plant_morph_type
  real(r8) :: LeafStalkArea_col                     !stalk area of combined, each PFT canopy
  real(r8) :: CanopyLeafArea_col                     !grid canopy leaf area, [m2 d-2]
  real(r8) :: StemArea_col                 !grid canopy stem area, [m2 d-2]
  real(r8) :: CanopyHeight_col                        !canopy height , [m]
  REAL(R8), pointer :: PARTS_brch(:,:,:)         => null()  !fraction of C allocated to each morph unit
  real(r8), pointer :: RootVolPerMassC_pft(:,:)       => null() !root volume:mass ratio, [m3 g-1]
  real(r8), pointer :: RootPorosity_pft(:,:)       => null() !root porosity, [m3 m-3]
  real(r8), pointer :: Root2ndXSecArea_pft(:,:)     => null() !root  cross-sectional area  secondary axes, [m2]
  real(r8), pointer :: Root1stXSecArea_pft(:,:)     => null() !root cross-sectional area primary axes, [m2]
  real(r8), pointer :: Root1stMaxRadius1_pft(:,:)     => null() !root diameter primary axes, [m]
  real(r8), pointer :: Root2ndMaxRadius1_pft(:,:)     => null() !root diameter secondary axes, [m]
  real(r8), pointer :: SeedCMass_pft(:)         => null() !grain size at seeding, [g]
  real(r8), pointer :: RootPoreTortu4Gas(:,:)      => null() !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8), pointer :: Root1stLen_rpvr(:,:,:,:)  => null() !root layer length primary axes, [m d-2]
  real(r8), pointer :: Root2ndLen_pvr(:,:,:,:)  => null() !root layer length secondary axes, [m d-2]
  real(r8), pointer :: RootLenPerPlant_pvr(:,:,:)    => null() !root layer length per plant, [m p-1]
  real(r8), pointer :: Root2ndAveLen_pvr(:,:,:)    => null() !root layer average length, [m]
  real(r8), pointer :: Root1stSpecLen_pft(:,:)     => null() !specific root length primary axes, [m g-1]
  real(r8), pointer :: Root2ndSpecLen_pft(:,:)     => null() !specific root length secondary axes, [m g-1]
  real(r8), pointer :: Root2ndXNum_rpvr(:,:,:,:)   => null() !root layer number secondary axes, [d-2]
  real(r8), pointer :: Root1stDepz_pft(:,:,:)    => null() !root layer depth, [m]
  real(r8), pointer :: ClumpFactorInit_pft(:)          => null() !initial clumping factor for self-shading in canopy layer, [-]
  real(r8), pointer :: ClumpFactorNow_pft(:)          => null() !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8), pointer :: RootBranchFreq_pft(:)         => null() !root brancing frequency, [m-1]
  real(r8), pointer :: HypoctoHeight_pft(:)        => null() !cotyledon height, [m]
  integer,  pointer :: iPlantNfixType_pft(:)        => null() !N2 fixation type
  integer,  pointer :: MY(:)           => null() !mycorrhizal type (no or yes)
  real(r8), pointer :: CanPHeight4WatUptake(:)        => null() !canopy height, [m]
  real(r8), pointer :: CanopyHeightZ_col(:)  => null() !canopy layer height , [m]
  real(r8), pointer :: CanopyStemArea_pft(:)        => null() !plant stem area, [m2 d-2]
  real(r8), pointer :: CanopyLeafArea_pft(:)        => null() !plant canopy leaf area, [m2 d-2]
  integer,  pointer :: MainBranchNum_pft(:)          => null() !number of main branch
  integer,  pointer :: MaxSoiL4Root_pft(:)           => null() !maximum soil layer number for all root axes
  integer,  pointer :: NIXBotRootLayer_pft(:)          => null() !maximum soil layer number for all root axes, [-]
  integer,  pointer :: NumRootAxes_pft(:)          => null() !root primary axis number
  integer,  pointer :: NumCogrothNode_pft(:)         => null() !number of concurrently growing nodes
  integer,  pointer :: BranchNumber_pft(:)          => null() !branch number
  integer,  pointer :: NumOfBranches_pft(:)          => null() !branch number
  integer,  pointer :: NIXBotRootLayer_rpft(:,:)       => null() !maximum soil layer number for root axes, [-]
  real(r8), pointer :: ShootNodeNum_brch(:,:)       => null() !shoot node number, [-]
  real(r8), pointer :: NodeNum2InitFloral_brch(:,:)      => null() !shoot node number at floral initiation, [-]
  real(r8), pointer :: NodeNumberAtAnthesis_brch(:,:)      => null() !shoot node number at anthesis, [-]
  real(r8), pointer :: SineBranchAngle_pft(:)        => null() !branching angle, [degree from horizontal]
  real(r8), pointer :: LeafAreaZsec_brch(:,:,:,:,:) => null() !leaf surface area, [m2 d-2]
  integer,  pointer :: KMinNumLeaf4GroAlloc_brch(:,:)     => null() !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,  pointer :: BranchNumber_brch(:,:)       => null() !branch number, [-]
  real(r8), pointer :: PotentialSeedSites_brch(:,:)      => null() !branch potential grain number, [d-2]
  real(r8), pointer :: SeedNumSet_brch(:,:)      => null() !branch grain number, [d-2]
  real(r8), pointer :: CanPBranchHeight(:,:)     => null() !branch height, [m]
  real(r8), pointer :: LeafAreaDying_brch(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: LeafAreaLive_brch(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: SinePetioleAngle_pft(:)        => null() !sheath angle, [degree from horizontal]
  real(r8), pointer :: CLASS(:,:)      => null() !fractionction of leaves in different angle classes, [-]
  real(r8), pointer :: CanopyStemAreaZ_pft(:,:)      => null() !plant canopy layer stem area, [m2 d-2]
  real(r8), pointer :: CanopyLeafAreaZ_pft(:,:)      => null() !canopy layer leaf area, [m2 d-2]
  real(r8), pointer :: LeafAreaNode_brch(:,:,:)    => null() !leaf area, [m2 d-2]
  real(r8), pointer :: CanopySeedNum_pft(:)         => null() !canopy grain number, [d-2]
  real(r8), pointer :: SeedDepth_pft(:)        => null() !seeding depth, [m]
  real(r8), pointer :: PlantinDepz_pft(:)       => null() !planting depth, [m]
  real(r8), pointer :: SeedMeanLen_pft(:)         => null() !seed length, [m]
  real(r8), pointer :: SeedVolumeMean_pft(:)         => null() !seed volume, [m3 ]
  real(r8), pointer :: SeedAreaMean_pft(:)         => null() !seed surface area, [m2]
  real(r8), pointer :: CanopyStemAareZ_col(:)        => null() !total stem area, [m2 d-2]
  real(r8), pointer :: PetoLen2Mass_pft(:)         => null() !petiole length:mass during growth, [m gC-1]
  real(r8), pointer :: NodeLenPergC(:)         => null() !internode length:mass during growth, [m gC-1]
  real(r8), pointer :: SLA1_pft(:)         => null() !leaf area:mass during growth, [m2 gC-1]
  real(r8), pointer :: CanopyLeafAareZ_col(:)        => null() !total leaf area, [m2 d-2]
  real(r8), pointer :: LeafStalkArea_pft(:)        => null() !plant leaf+stem/stalk area, [m2 d-2]
  real(r8), pointer :: InternodeHeightDying_brch(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: PetoleLensNode_brch(:,:,:)    => null() !sheath height, [m]
  real(r8), pointer :: LiveInterNodeHight_brch(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: StemAreaZsec_brch(:,:,:,:)  => null() !stem surface area, [m2 d-2]
  real(r8), pointer :: CanopyLeafArea_lpft(:,:,:,:)  => null() !layer/node/branch leaf area, [m2 d-2]
  real(r8), pointer :: CanopyStalkArea_lbrch(:,:,:)    => null() !plant canopy layer branch stem area, [m2 d-2]
  real(r8), pointer :: ClumpFactor_pft(:)           => null() !clumping factor for self-shading in canopy layer, [-]
  real(r8), pointer :: ShootNodeNumAtPlanting_pft(:)         => null() !number of nodes in seed, [-]
  real(r8), pointer :: CanopyHeight_pft(:)           => null() !canopy height, [m]
  integer,  pointer :: NGTopRootLayer_pft(:)           => null() !soil layer at planting depth, [-]
  integer,  pointer :: KLeafNumber_brch(:,:)      => null() !leaf number, [-]
  real(r8), pointer :: NumOfLeaves_brch(:,:)       => null() !leaf number, [-]
  integer , pointer :: iPlantGrainType_pft(:)        => null() !grain type (below or above-ground)
  real(r8), pointer :: MaxPotentSeedNumber_pft(:)         => null() !maximum grain node number per branch, [-]
  real(r8), pointer :: MaxSeedNumPerSite_pft(:)         => null() !maximum grain number per node , [-]
  real(r8), pointer :: rLen2WidthLeaf_pft(:)         => null() !leaf length:width ratio, [-]
  real(r8), pointer :: MaxSeedCMass_pft(:)         => null() !maximum grain size   , [g]
  real(r8), pointer :: Root1stRadius_pvr(:,:,:)    => null() !root layer diameter primary axes, [m ]
  real(r8), pointer :: Root2ndRadius_pvr(:,:,:)    => null() !root layer diameter secondary axes, [m ]
  real(r8), pointer :: RootRaidus_rpft(:,:)      => null() !root internal radius, [m]
  real(r8), pointer :: Root1stMaxRadius_pft(:,:)     => null() !maximum radius of primary roots, [m]
  real(r8), pointer :: Root2ndMaxRadius_pft(:,:)     => null() !maximum radius of secondary roots, [m]
  real(r8), pointer :: RootRadialResist_pft(:,:)       => null() !root radial resistivity, [MPa h m-2]
  real(r8), pointer :: RootAxialResist_pft(:,:)       => null() !root axial resistivity, [MPa h m-4]
  real(r8), pointer :: totRootLenDens_vr(:)        => null() !total root length density, [m m-3]
  real(r8), pointer :: Root1stXNumL_pvr(:,:,:)     => null() !root layer number primary axes, [d-2]
  real(r8), pointer :: Root2ndXNum_pvr(:,:,:)     => null() !root layer number axes, [d-2]
  real(r8), pointer :: RootLenDensPerPlant_pvr(:,:,:)    => null() !root layer length density, [m m-3]
  real(r8), pointer :: RootPoreVol_pvr(:,:,:)    => null() !root layer volume air, [m2 d-2]
  real(r8), pointer :: RootVH2O_pvr(:,:,:)    => null() !root layer volume water, [m2 d-2]
  real(r8), pointer :: RootAreaPerPlant_pvr(:,:,:)    => null() !root layer area per plant, [m p-1]
  contains
    procedure, public :: Init    => plt_morph_init
    procedure, public :: Destroy => plt_morph_destroy
  end type plant_morph_type

  type, public :: plant_pheno_type
  real(r8), pointer :: fTgrowRootP_vr(:,:)   => null()     !root layer temperature growth functiom, [-]
  real(r8), pointer :: ShutRutNonstElmntConducts_pft(:)    => null()     !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(r8), pointer :: GrainFillRate25C_pft(:)    => null()     !maximum rate of fill per grain, [g h-1]
  real(r8), pointer :: TempOffset_pft(:)    => null()     !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
  real(r8), pointer :: PlantO2Stress_pft(:)     => null()     !plant O2 stress indicator, []
  real(r8), pointer :: MinNonstC2InitBranch_pft(:)       => null()     !branch nonstructural C content required for new branch, [gC gC-1]
  real(r8), pointer :: MinNonstC2InitRoot_pft(:)       => null()     !threshold root nonstructural C content for initiating new root axis, [gC gC-1]
  real(r8), pointer :: LeafElmntRemobFlx_brch(:,:,:) => null()    !element translocated from leaf during senescence, [g d-2 h-1]
  real(r8), pointer :: PetioleChemElmRemobFlx_brch(:,:,:) => null()    !element translocated from sheath during senescence, [g d-2 h-1]
  real(r8), pointer :: TC4LeafOut_pft(:)     => null()     !threshold temperature for spring leafout/dehardening, [oC]
  real(r8), pointer :: TCGroth_pft(:)     => null()     !canopy growth temperature, [oC]
  real(r8), pointer :: TC4LeafOff_pft(:)     => null()     !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8), pointer :: TKGroth_pft(:)     => null()     !canopy growth temperature, [K]
  real(r8), pointer :: fTCanopyGroth_pft(:)    => null()     !canopy temperature growth function, [-]
  real(r8), pointer :: HoursTooLowPsiCan_pft(:)    => null()     !canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY, []
  real(r8), pointer :: TCChill4Seed_pft(:)     => null()     !temperature below which seed set is adversely affected, [oC]
  real(r8), pointer :: iPlantThermoAdaptZone_pft(:)    => null()     !plant thermal adaptation zone, [-]
  real(r8), pointer :: PlantInitThermoAdaptZone(:)   => null()     !initial plant thermal adaptation zone, [-]
  real(r8), pointer :: HighTempLimitSeed_pft(:)     => null()     !temperature above which seed set is adversely affected, [oC]
  real(r8), pointer :: SeedTempSens_pft(:)    => null()     !sensitivity to HTC (seeds oC-1 above HTC)
  integer,  pointer :: iPlantState_pft(:)    => null()     !flag for species death
  integer,  pointer :: IsPlantActive_pft(:)   => null()     !flag for living pft
  real(r8), pointer :: NetCumElmntFlx2Plant_pft(:,:) => null()     !effect of canopy element status on seed set , []
  real(r8), pointer :: RefNodeInitRate_pft(:)    => null()     !rate of node initiation, [h-1 at 25 oC]
  real(r8), pointer :: RefLeafAppearRate_pft(:)    => null()     !rate of leaf initiation, [h-1 at 25 oC]
  real(r8), pointer :: CriticPhotoPeriod_pft(:)     => null()     !critical daylength for phenological progress, [h]
  real(r8), pointer :: PhotoPeriodSens_pft(:)    => null()     !difference between current and critical daylengths used to calculate  phenological progress [h]
  integer,  pointer :: iPlantBranchState_brch(:,:)  => null()  !flag to detect branch death , [-]
  integer,  pointer :: doRemobilization_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doPlantLeaveOff_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doPlantLeafOut_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doInitLeafOut_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doSenescence_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: Prep4Literfall_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: Hours4LiterfalAftMature_brch(:,:)  => null()  !branch phenology flag, [h]
  integer,  pointer :: KHiestGroLeafNode_brch(:,:)  => null()  !leaf growth stage counter, [-]
  integer,  pointer :: iPlantPhenolType_pft(:)    => null()  !climate signal for phenological progress: none, temperature, water stress
  integer,  pointer :: iPlantTurnoverPattern_pft(:)    => null()  !phenologically-driven above-ground turnover: all, foliar only, none
  integer,  pointer :: iPlantShootState_pft(:)    => null()  !flag to detect canopy death
  integer,  pointer :: iPlantPhenolPattern_pft(:)    => null()  !plant growth habit: annual or perennial
  integer,  pointer :: iPlantRootState_pft(:)    => null()  !flag to detect root system death
  integer,  pointer :: iPlantDevelopPattern_pft(:)    => null()  !plant growth habit (determinate or indeterminate)
  integer,  pointer :: iPlantPhotoperiodType_pft(:)    => null()  !photoperiod type (neutral, long day, short day)
  integer,  pointer :: iPlantRootProfile_pft(:)    => null()  !plant growth type (vascular, non-vascular)
  integer,  pointer :: doInitPlant_pft(:)    => null()  !PFT initialization flag:0=no,1=yes
  integer,  pointer :: KLowestGroLeafNode_brch(:,:) => null()  !leaf growth stage counter, [-]
  integer,  pointer :: iPlantCalendar_brch(:,:,:) => null()  !plant growth stage, [-]
  real(r8), pointer :: TotalNodeNumNormByMatgrp_brch(:,:) => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: TotReproNodeNumNormByMatrgrp_brch(:,:) => null()  !normalized node number during reproductive growth stages , [-]
  real(r8), pointer :: LeafNumberAtFloralInit_brch(:,:)  => null()  !leaf number at floral initiation, [-]
  real(r8), pointer :: Hours4LenthenPhotoPeriod_brch(:,:)   => null()  !initial heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: Hours4ShortenPhotoPeriod_brch(:,:)   => null()  !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: Hours4Leafout_brch(:,:)   => null()  !heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: HourReq4LeafOut_brch(:,:)   => null()  !hours above threshold temperature required for spring leafout/dehardening, [-]
  real(r8), pointer :: Hours4LeafOff_brch(:,:)   => null()  !cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: HourReq4LeafOff_brch(:,:)   => null()  !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8), pointer :: Hours2LeafOut_brch(:,:)   => null()  !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8), pointer :: HoursDoingRemob_brch(:,:)   => null()  !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  real(r8), pointer :: HourlyNodeNumNormByMatgrp_brch(:,:) => null()  !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8), pointer :: dReproNodeNumNormByMatG_brch(:,:) => null()  !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8), pointer :: MatureGroup_brch(:,:)  => null()  !plant maturity group, [-]
  real(r8), pointer :: NodeNumNormByMatgrp_brch(:,:)  => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: ReprodNodeNumNormByMatrgrp_brch(:,:)  => null()  !normalized node number during reproductive growth stages, [-]
  real(r8), pointer :: HourFailGrainFill_brch(:,:)   => null()  !flag to detect physiological maturity from  grain fill , [-]
  real(r8), pointer :: MatureGroup_pft(:)   => null()  !acclimated plant maturity group, [-]
  contains
    procedure, public :: Init    =>  plt_pheno_init
    procedure, public :: Destroy =>  plt_pheno_destroy
  end type plant_pheno_type

  type, public :: plant_soilchem_type
  real(r8), pointer :: FracBulkSOMC_vr(:,:)  => null()  !fraction of total organic C in complex, [-]
  real(r8), pointer :: ElmAllocmat4Litr(:,:,:,:)=> null() !litter kinetic fraction, [-]
  real(r8), pointer :: TScal4Difsvity_vr(:)     => null()  !temperature effect on diffusivity
  real(r8), pointer :: THETPM(:,:) => null()  !soil air-filled porosity, [m3 m-3]
  real(r8), pointer :: DiffusivitySolutEff(:,:)   => null()  !coefficient for dissolution - volatilization, []
  real(r8), pointer :: SoilResit4RootPentrate_vr(:)     => null()  !soil hydraulic resistance, [MPa h m-2]
  real(r8), pointer :: SoiBulkDensity_vr(:)     => null()  !soil bulk density, [Mg m-3]
  real(r8), pointer :: trc_solcl_vr(:,:) => null() !aqueous tracer concentration [g m-3]
  real(r8), pointer :: trc_gascl_vr(:,:) => null() !gaseous tracer concentration [g m-3]

  real(r8), pointer :: CSoilOrgM_vr(:,:)    => null()  !soil organic C content [gC kg soil-1]

  real(r8), pointer :: HydroCondMicP4RootUptake_vr(:)     => null()  !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]

  real(r8), pointer :: GasDifc_vr(:,:)=> null()  !gaseous diffusivity [m2 h-1]
  real(r8), pointer :: SoluteDifusvty_vr(:,:)=> null()  !aqueous diffusivity [m2 h-1]

  real(r8), pointer :: trc_gasml_vr(:,:)=> null()!gas layer mass [g d-2]
  real(r8), pointer :: GasSolbility_vr(:,:) => null() !gas solubility, [m3 m-3]
  real(r8), pointer :: THETW_vr(:)          => null()  !volumetric water content [m3 m-3]
  real(r8), pointer :: THETY_vr(:)          => null()  !air-dry water content, [m3 m-3]
  real(r8), pointer :: VLSoilPoreMicP_vr(:) => null()  !volume of soil layer	m3 d-2
  real(r8), pointer :: trcs_VLN_vr(:,:)     => null()
  real(r8), pointer :: VLSoilMicP_vr(:) => null()  !total micropore volume in layer [m3 d-2]
  real(r8), pointer :: VLiceMicP_vr(:)  => null()  !soil micropore ice content   [m3 d-2]
  real(r8), pointer :: VLWatMicP_vr(:)  => null()  !soil micropore water content [m3 d-2]
  real(r8), pointer :: VLMicP_vr(:)     => null()  !total volume in micropores [m3 d-2]

  real(r8), pointer :: DOM_vr(:,:,:)     => null()  !dissolved organic C micropore	[gC d-2]
  real(r8), pointer :: trc_solml_vr(:,:) => null() !aqueous tracer [g d-2]
  contains
    procedure, public :: Init => plt_soilchem_init
    procedure, public :: Destroy => plt_soilchem_destroy
  end type plant_soilchem_type

  type, public :: plant_allometry_type
  real(r8), pointer :: CPRTS_pft(:)    => null()  !root P:C ratio x root growth yield, [-]
  real(r8), pointer :: CNRTS_pft(:)    => null()  !root N:C ratio x root growth yield, [-]
  real(r8), pointer :: NodulerNC_pft(:)     => null()  !nodule N:C ratio, [gN gC-1]
  real(r8), pointer :: NodulerPC_pft(:)     => null()  !nodule P:C ratio, [gP gC-1]
  real(r8), pointer :: RootrNC_pft(:)     => null()  !root N:C ratio, [gN gC-1]
  real(r8), pointer :: RootrPC_pft(:)     => null()  !root P:C ratio, [gP gC-1]
  real(r8), pointer :: rCNNonstRemob_pft(:)     => null()  !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: rCPNonstRemob_pft(:)     => null()  !C:P ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: NoduGrowthYield_pft(:)     => null()  !nodule growth yield, [g g-1]
  real(r8), pointer :: RootBiomGrosYld_pft(:)     => null()  !root growth yield, [g g-1]
  real(r8), pointer :: rPCEar_pft(:)    => null()  !ear P:C ratio, [gP gC-1]
  real(r8), pointer :: PetioleBiomGrowthYield(:)    => null()  !sheath growth yield, [g g-1]
  real(r8), pointer :: rPCHusk_pft(:)    => null()  !husk P:C ratio, [gP gC-1]
  real(r8), pointer :: StalkBiomGrowthYield(:)    => null()  !stalk growth yield, [gC gC-1]
  real(r8), pointer :: HuskBiomGrowthYield(:)    => null()  !husk growth yield, [gC gC-1]
  real(r8), pointer :: ReserveBiomGrowthYield(:)    => null()  !reserve growth yield, [gC gC-1]
  real(r8), pointer :: GrainBiomGrowthYield(:)     => null()  !grain growth yield, [gC gC-1]
  real(r8), pointer :: EarBiomGrowthYield(:)    => null()  !ear growth yield, [gC gC-1]
  real(r8), pointer :: rNCHusk_pft(:)    => null()  !husk N:C ratio, [gN gC-1]
  real(r8), pointer :: rNCReserve_pft(:)    => null()  !reserve N:C ratio, [gN gC-1]
  real(r8), pointer :: rNCEar_pft(:)    => null()  !ear N:C ratio, [gN gC-1]
  real(r8), pointer :: rPCReserve_pft(:)    => null()  !reserve P:C ratio, [gP gC-1]
  real(r8), pointer :: CPGR(:)     => null()  !grain P:C ratio, [gP gP-1]
  real(r8), pointer :: rNCStalk_pft(:)    => null()  !stalk N:C ratio, [gN gC-1]
  real(r8), pointer :: FracShootLeafElmAlloc2Litr(:,:) => null()  !woody element allocation
  real(r8), pointer :: FracShootStalkElmAlloc2Litr(:,:) => null()  !leaf element allocation
  real(r8), pointer :: FracRootElmAlloc2Litr(:,:) => null()  !C woody fraction in root
  real(r8), pointer :: FracRootStalkElmAlloc2Litr(:,:) => null()  !woody element allocation
  real(r8), pointer :: FracHour4LeafoffRemob(:)     => null()  !allocation parameter
  real(r8), pointer :: LeafBiomGrowthYield(:)     => null()  !leaf growth yield, [g g-1]
  real(r8), pointer :: CNGR(:)     => null()  !grain N:C ratio, [g g-1]
  real(r8), pointer :: CPLF(:)     => null()  !maximum leaf P:C ratio, [g g-1]
  real(r8), pointer :: CPSHE(:)    => null()  !sheath P:C ratio, [g g-1]
  real(r8), pointer :: CNSHE(:)    => null()  !sheath N:C ratio, [g g-1]
  real(r8), pointer :: rPCStalk_pft(:)    => null()  !stalk P:C ratio, [g g-1]
  real(r8), pointer :: CNLF(:)     => null()  !maximum leaf N:C ratio, [g g-1]
  real(r8), pointer :: GrainSeedBiomCMean_brch(:,:)  => null()  !maximum grain C during grain fill, [g d-2]
  real(r8), pointer :: FracGroth2Node_pft(:)     => null()  !parameter for allocation of growth to nodes, [-]
  real(r8), pointer ::RootFracRemobilizableBiom(:)    => null()  !fraction of remobilizable nonstructural biomass in root, [-]

  contains
    procedure, public :: Init => plt_allom_init
    procedure, public :: Destroy => plt_allom_destroy
  end type plant_allometry_type

  type, public :: plant_biom_type
  real(r8), pointer :: StomatalStress_pft(:)  => null()   !stomatal stress from root turgor [0-1]  
  real(r8), pointer :: RootMassElm_pvr(:,:,:)                 => null()
  real(r8), pointer :: StandingDeadStrutElms_col(:)           => null()    !total standing dead element,                [g d-2]
  real(r8), pointer :: ZERO4LeafVar_pft(:)                    => null()    !threshold zero for leaf calculation
  real(r8), pointer :: ZERO4Groth_pft(:)                      => null()    !threshold zero for p calculation
  real(r8), pointer :: RootNodulNonstElms_pvr(:,:,:)          => null()    !root  layer nonstructural element,          [g d-2]
  real(r8), pointer :: RootNodulStrutElms_pvr(:,:,:)          => null()    !root layer nodule element,                  [g d-2]
  real(r8), pointer :: CanopyLeafCLyr_pft(:,:)                => null()    !canopy layer leaf C,                        [g d-2]
  real(r8), pointer :: RootMycoNonstElms_pft(:,:,:)           => null()
  real(r8), pointer :: RootMyco2ndStrutElms_rpvr(:,:,:,:,:)   => null()    !root layer element secondary axes,          [g d-2]
  real(r8), pointer :: RootMyco1stStrutElms_rpvr(:,:,:,:,:)   => null()    !root layer element primary axes,            [g d-2]
  real(r8), pointer :: RootMyco1stElm_raxs(:,:,:,:)           => null()    !root C primary axes,                        [g d-2]
  real(r8), pointer :: StandDeadKCompElms_pft(:,:,:)          => null()    !standing dead element fraction,             [g d-2]
  real(r8), pointer :: CanopyNonstElmConc_pft(:,:)            => null()    !canopy nonstructural element concentration, [g d-2]
  real(r8), pointer :: CanopyNonstElms_pft(:,:)               => null()    !canopy nonstructural element concentration, [g d-2]
  real(r8), pointer :: CanopyNodulNonstElms_pft(:,:)          => null()    !canopy nodule nonstructural element,        [g d-2]
  real(r8), pointer :: NoduleNonstructCconc_pft(:)            => null()    !nodule nonstructural C,                     [gC d-2]
  real(r8), pointer :: RootMycoActiveBiomC_pvr(:,:,:)         => null()    !root layer structural C,                    [g d-2]
  real(r8), pointer ::  PopuRootMycoC_pvr(:,:,:)              => null()    !root layer C,                               [g d-2]
  real(r8), pointer :: RootProteinC_pvr(:,:,:)                => null()    !root layer protein C,                       [g d-2]
  real(r8), pointer :: RootProteinConc_pvr(:,:,:)             => null()    !root layer protein C concentration,         [g g-1]
  real(r8), pointer :: RootMycoNonstElms_rpvr(:,:,:,:)        => null()    !root  layer nonstructural element,          [g d-2]
  real(r8), pointer :: RootNonstructElmConc_pvr(:,:,:,:)      => null()    !root  layer nonstructural C concentration,  [g g-1]
  real(r8), pointer :: LeafPetoNonstElmConc_brch(:,:,:)       => null()    !branch nonstructural C concentration,       [g d-2]
  real(r8), pointer :: InternodeStrutElms_brch(:,:,:,:)       => null()    !internode C,                                [g d-2]
  real(r8), pointer :: LeafElmntNode_brch(:,:,:,:)            => null()    !leaf element,                               [g d-2]
  real(r8), pointer :: LeafProteinCNode_brch(:,:,:)           => null()    !layer leaf protein C,                       [g d-2]
  real(r8), pointer :: PetioleElmntNode_brch(:,:,:,:)         => null()  !sheath element,                               [g d-2]
  real(r8), pointer :: PetoleProteinCNode_brch(:,:,:)         => null()    !layer sheath protein C,                     [g d-2]
  real(r8), pointer :: LeafElmsByLayerNode_brch(:,:,:,:,:) => null()    !layer leaf element,                         [g d-2]
  real(r8), pointer :: tCanLeafC_cl(:)                        => null()  !total leaf mass,                              [gC d-2]
  real(r8), pointer :: StandingDeadInitC_pft(:)               => null()  !initial standing dead C,                      [g C m-2]
  real(r8), pointer :: RootNodulElms_pft(:,:)                 => null()
  real(r8), pointer :: RootElms_pft(:,:)                      => null()  !plant root element,                           [gE d-2]
  real(r8), pointer :: RootElmsBeg_pft(:,:)                   => null()  !plant root element,                           [gE d-2]
  real(r8), pointer :: StandDeadStrutElmsBeg_pft(:,:)         => null()  !standing dead element,                        [gE d-2]
  real(r8), pointer :: CanopyNodulElms_pft(:,:)               => null()
  real(r8), pointer :: RootStrutElms_pft(:,:)                 => null()  !plant root structural element,                [gC d-2]
  real(r8), pointer :: SeedCPlanted_pft(:)                    => null()  !plant stored nonstructural C at planting,     [gC d-2]
  real(r8), pointer :: SeasonalNonstElms_pft(:,:)             => null()  !plant stored nonstructural element,           [gC d-2]
  real(r8), pointer :: CanopyLeafShethC_pft(:)                => null()  !canopy leaf + sheath C,                       [g d-2]
  real(r8), pointer :: ShootStrutElms_pft(:,:)                => null()  !canopy shoot C,                               [g d-2]
  real(r8), pointer :: AvgCanopyBiomC2Graze_pft(:)            => null()  !landscape average canopy shoot C,             [g d-2]
  real(r8), pointer :: StandDeadStrutElms_pft(:,:)            => null()  !standing dead element,                        [g d-2]
  real(r8), pointer :: NodulStrutElms_pft(:,:)                => null()  !root total nodule mass,                       element [g d-2]
  real(r8), pointer :: CanopyNonstElms_brch(:,:,:)            => null()  !branch nonstructural element,                 [g d-2]
  real(r8), pointer :: ShootElms_brch(:,:,:)                  => null()      !branch shoot elemental biomass,           [g d-2]
  real(r8), pointer :: ShootC4NonstC_brch(:,:)                => null()  !branch shoot nonstrucal elelment [g d-2]
  real(r8), pointer :: CanopyNodulNonstElms_brch(:,:,:)       => null()  !branch nodule nonstructural element,          [g d-2]
  real(r8), pointer :: LeafPetolBiomassC_brch(:,:)            => null()  !plant branch leaf + sheath C,                 [g d-2]
  real(r8), pointer :: StalkRsrvElms_brch(:,:,:)              => null()  !branch reserve element,                       [g d-2]
  real(r8), pointer :: LeafStrutElms_brch(:,:,:)              => null()   !branch leaf element,                         [g d-2]
  real(r8), pointer :: CanopyNodulStrutElms_brch(:,:,:)       => null()   !branch nodule element,                       [g d-2]
  real(r8), pointer :: PetoleStrutElms_brch(:,:,:)            => null()   !branch sheath element,                       [g d-2]
  real(r8), pointer :: EarStrutElms_brch(:,:,:)               => null()   !branch ear C,                                [g d-2]
  real(r8), pointer :: HuskStrutElms_brch(:,:,:)              => null()   !branch husk element,                         [g d-2]
  real(r8), pointer :: GrainStrutElms_brch(:,:,:)             => null()   !branch grain element,                        [g d-2]
  real(r8), pointer :: StalkStrutElms_brch(:,:,:)             => null()   !branch stalk element,                        [g d-2]
  real(r8), pointer :: ShootStrutElms_brch(:,:,:)             => null()   !branch shoot C,                              [g d-2]
  real(r8), pointer :: LeafChemElmRemob_brch(:,:,:)           => null()   !branch leaf structural element,              [g d-2]
  real(r8), pointer :: PetioleChemElmRemob_brch(:,:,:)        => null()   !branch sheath structural element,            [g d-2]
  real(r8), pointer :: SenecStalkStrutElms_brch(:,:,:)        => null()   !branch stalk structural element,             [g d-2]
  real(r8), pointer :: StalkBiomassC_brch(:,:)                => null()   !branch active stalk C,                       [g d-2]
  real(r8), pointer :: StalkStrutElms_pft(:,:)                => null()   !canopy stalk element,                        [g d-2]
  real(r8), pointer :: ShootElms_pft(:,:)                     => null()
  real(r8), pointer :: ShootElmsBeg_pft(:,:)                  => null()
  real(r8), pointer :: CanopyStalkC_pft(:)                    => null()   !canopy active stalk C,                       [g d-2
  real(r8), pointer :: LeafStrutElms_pft(:,:)   => null()   !canopy leaf elements,   [g d-2]
  real(r8), pointer :: PetoleStrutElms_pft(:,:) => null()   !canopy sheath element,  [g d-2]
  real(r8), pointer :: StalkRsrvElms_pft(:,:)   => null()   !canopy reserve element, [g d-2]
  real(r8), pointer :: HuskStrutElms_pft(:,:)   => null()   !canopy husk element,    [g d-2]
  real(r8), pointer :: RootBiomCPerPlant_pft(:) => null()   !root C per plant,       [g p-1]
  real(r8), pointer :: GrainStrutElms_pft(:,:)  => null()   !canopy grain element,   [g d-2]
  real(r8), pointer :: EarStrutElms_pft(:,:)    => null()   !canopy ear element,     [g d-2]
  contains
    procedure, public :: Init => plt_biom_init
    procedure, public :: Destroy => plt_biom_destroy
  end type plant_biom_type

  type, public :: plant_ew_type
  real(r8) :: SnowDepth          !snowpack depth, [m]
  real(r8) :: VcumWatSnow_col        !water volume in snowpack, [m3 d-2]
  real(r8) :: VcumDrySnoWE_col       !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8) :: Canopy_Heat_Sens_col               !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VLHeatCapSnowMin_col    !minimum snowpack heat capacity [MJ d-2 K-1]
  real(r8) :: VLHeatCapSurfSnow_col    !snowpack heat capacity, [MJ m-3 K-1]
  real(r8) :: CanopyHeatStor_col    !total canopy heat content, [MJ  d-2]
  real(r8) :: LWRadCanG     !grid total canopy LW emission, [MJ d-2 h-1]
  real(r8) :: QvET_col    !total canopy evaporation + transpiration, [m3 d-2]
  real(r8) :: VapXAir2Canopy_col    !grid canopy evaporation, [m3 d-2]
  real(r8) :: HeatFlx2Canopy_col    !total canopy heat flux, [MJ  d-2]
  real(r8) :: H2OLoss_CumYr_col     !total subsurface water flux, [m3 d-2]
  real(r8) :: VPA       !vapor concentration, [m3 m-3]
  real(r8) :: TairK       !air temperature, [K]
  real(r8) :: CanWat_col    !total canopy water content stored with dry matter, [m3 d-2]
  real(r8) :: Eco_Heat_Latent_col       !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8) :: Canopy_Heat_Latent_col      !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VcumIceSnow_col     !ice volume in snowpack, [m3 d-2]
  real(r8) :: TKSnow       !snow temperature, [K]
  real(r8) :: CanH2OHeldVg    !canopy surface water content, [m3 d-2]
  real(r8) :: Eco_Heat_Sens_col       !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8) :: RoughHeight        !canopy surface roughness height, [m]
  real(r8) :: ZERO4PlantDisplace_col        !zero plane displacement height, [m]
  real(r8) :: AbvCanopyBndlResist_col       !isothermal boundary layer resistance, [h m-1]
  real(r8) :: RIB       !Richardson number for calculating boundary layer resistance, [-]
  real(r8) :: Eco_Heat_Grnd_col       !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8), pointer :: Transpiration_pft(:)            => null()    !canopy transpiration,                                         [m2 d-2 h-1]
  real(r8), pointer :: PSICanopyTurg_pft(:)            => null()    !plant canopy turgor water potential,                          [MPa]
  real(r8), pointer :: PrecIntcptByCanopy_pft(:)       => null()    !water flux into canopy,                                       [m3 d-2 h-1]
  real(r8), pointer :: PSICanopy_pft(:)                => null()    !canopy total water potential,                                 [Mpa]
  real(r8), pointer :: VapXAir2Canopy_pft(:)           => null()    !canopy evaporation,                                           [m2 d-2 h-1]
  real(r8), pointer :: HeatStorCanopy_pft(:)                 => null()    !canopy storage heat flux,                                     [MJ d-2 h-1]
  real(r8), pointer :: EvapTransHeat_pft(:)            => null()    !canopy latent heat flux,                                      [MJ d-2 h-1]
  real(r8), pointer :: RAZ(:)                          => null()    !canopy roughness height,                                      [m]
  real(r8), pointer :: TKS_vr(:)                       => null()    !mean annual soil temperature,                                 [K]
  real(r8), pointer :: PSICanPDailyMin(:)              => null()    !minimum daily canopy water potential,                         [MPa]
  real(r8), pointer :: TCelciusCanopy_pft(:)           => null()    !canopy temperature,                                           [oC]
  real(r8), pointer :: DeltaTKC_pft(:)                 => null()    !change in canopy temperature,                                 [K]
  real(r8), pointer :: ENGYX_pft(:)                    => null()    !canopy heat storage from previous time step,                  [MJ d-2]
  real(r8), pointer :: TKC(:)                          => null()    !canopy temperature,                                           [K]
  real(r8), pointer :: PSICanopyOsmo_pft(:)            => null()    !canopy osmotic water potential,                               [Mpa]
  real(r8), pointer :: CanOsmoPsi0pt_pft(:)            => null()    !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8), pointer :: HeatXAir2PCan_pft(:)                => null()    !canopy sensible heat flux,                                    [MJ d-2 h-1]
  real(r8), pointer :: TKCanopy_pft(:)                 => null()    !canopy temperature,                                           [K]
  real(r8), pointer :: PSIRoot_pvr(:,:,:)              => null() !root total water potential,                                      [Mpa]
  real(r8), pointer :: PSIRootOSMO_vr(:,:,:)           => null() !root osmotic water potential,                                    [Mpa]
  real(r8), pointer :: PSIRootTurg_vr(:,:,:)           => null() !root turgor water potential,                                     [Mpa]
  real(r8), pointer :: AllPlantRootH2OUptake_vr(:,:,:) => null() !root water uptake,                                               [m2 d-2 h-1]
  real(r8), pointer :: THeatRootUptake_vr(:)           => null()   !total root heat uptake,                                        [MJ d-2]
  real(r8), pointer :: TPlantRootH2OUptake_vr(:)    => null()   !total root water uptake,                                       [m3 d-2]
  real(r8), pointer :: WatByPCanopy_pft(:)             => null()  !canopy surface water content,                                   [m3 d-2]
  real(r8), pointer :: CanopyWater_pft(:)              => null()  !canopy water content,                                           [m3 d-2]
  real(r8), pointer :: VHeatCapCanP_pft(:)             => null()  !canopy heat capacity,                                           [MJ d-2 K-1]
  real(r8), pointer :: TotalSoilH2OPSIMPa_vr(:)        => null()  !soil micropore total water potential [MPa]
  real(r8), pointer :: ETCanopy_CumYr_pft(:)           => null()  !total transpiration,                                            [m H2O d-2]
  real(r8), pointer :: QdewCanopy_pft(:)               => null()
  contains
    procedure, public :: Init => plt_ew_init
    procedure, public :: Destroy=> plt_ew_destroy
  end type plant_ew_type

  type, public :: plant_disturb_type

  real(r8) :: XCORP     !factor for surface litter incorporation and soil mixing
  integer  :: IYTYP      !fertilizer release type from fertilizer input file
  real(r8) :: DCORP      !soil mixing fraction with tillage, [-]
  integer  :: iSoilDisturbType_col      !soil disturbance type, [-]
  real(r8), pointer :: EcoHavstElmnt_CumYr_col(:)   => null()  !ecosystem harvest element, [gC d-2]
  real(r8), pointer :: O2ByFire_CumYr_pft(:)    => null()  !plant O2 uptake from fire, [g d-2 ]
  integer,  pointer :: iHarvestDay_pft(:)    => null()  !alternate day of harvest, [-]
  real(r8), pointer :: FERT(:)     => null()  !fertilizer application, [g m-2]
  integer,  pointer :: iYearPlanting_pft(:)     => null()  !year of planting
  integer,  pointer :: iPlantingDay_pft(:)    => null()  !alternate day of planting
  integer,  pointer :: iPlantingYear_pft(:)     => null()  !alternate year of planting
  integer,  pointer :: iHarvestYear_pft(:)     => null()  !alternate year of harvest
  integer,  pointer :: iDayPlanting_pft(:)    => null()  !day of planting
  integer,  pointer :: iDayPlantHarvest_pft(:)    => null()  !day of harvest
  integer,  pointer :: iYearPlantHarvest_pft(:)     => null()  !year of harvest
  integer,  pointer :: iHarvstType_pft(:)    => null()  !type of harvest
  integer,  pointer :: jHarvst_pft(:)    => null()  !flag for stand replacing disturbance
  real(r8), pointer :: FracBiomHarvsted(:,:,:)=> null()  !harvest efficiency, [-]
  real(r8), pointer :: FracCanopyHeightCut_pft(:)     => null()  !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  real(r8), pointer :: THIN_pft(:)     => null()  !thinning of plant population, [-]
  real(r8), pointer :: CH4ByFire_CumYr_pft(:)    => null()  !plant CH4 emission from fire, [g d-2 ]
  real(r8), pointer :: CO2ByFire_CumYr_pft(:)    => null()  !plant CO2 emission from fire, [g d-2 ]
  real(r8), pointer :: N2ObyFire_CumYr_pft(:)    => null()  !plant N2O emission from fire, [g d-2 ]
  real(r8), pointer :: NH3byFire_CumYr_pft(:)    => null()  !plant NH3 emission from fire, [g d-2 ]
  real(r8), pointer :: PO4byFire_CumYr_pft(:)    => null()  !plant PO4 emission from fire, [g d-2 ]
  real(r8), pointer :: EcoHavstElmntCum_pft(:,:) => null()  !total plant element harvest, [gC d-2 ]
  real(r8), pointer :: EcoHavstElmnt_CumYr_pft(:,:)  => null()  !plant element harvest, [g d-2 ]
  contains
    procedure, public :: Init    =>  plt_disturb_init
    procedure, public :: Destroy => plt_disturb_destroy
  end type plant_disturb_type

  type, public :: plant_bgcrate_type
  real(r8) :: Eco_NBP_CumYr_col      !total NBP, [g d-2]

  real(r8) :: NetCO2Flx2Canopy_col     !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8) :: ECO_ER_col      !ecosystem respiration, [g d-2 h-1]
  real(r8) :: Eco_AutoR_CumYr_col      !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8) :: TH2GZ     !total root H2 flux, [g d-2]
  real(r8) :: Canopy_NEE_col     !total net CO2 fixation, [gC d-2]
  real(r8), pointer :: LitrFallStrutElms_col(:)        => null() !total LitrFall element,                       [g d-2 h-1]
  real(r8), pointer :: NetPrimProduct_pft(:)           => null()   !total net primary productivity,             [gC d-2]
  real(r8), pointer :: NH3Dep2Can_pft(:)               => null()   !canopy NH3 flux,                            [g d-2 h-1]
  real(r8), pointer :: tRootMycoExud2Soil_vr(:,:,:)    => null()  !total root element exchange,                 [g d-2 h-1]
  real(r8), pointer :: RootN2Fix_pvr(:,:)              => null()  !root N2 fixation,                            [gN d-2 h-1]
  real(r8), pointer :: CanopyRespC_CumYr_pft(:)              => null()  !total autotrophic respiration,               [gC d-2 ]
  real(r8), pointer :: LitrfalStrutElms_pvr(:,:,:,:,:) => null()  !plant LitrFall element,                      [g d-2 h-1]
  real(r8), pointer :: REcoO2DmndResp_vr(:)            => null()  !total root + microbial O2 uptake,            [g d-2 h-1]
  real(r8), pointer :: REcoNH4DmndBand_vr(:)           => null()   !total root + microbial NH4 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: REcoH1PO4DmndSoil_vr(:)         => null()   !HPO4 demand in non-band by all microbial,   root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: REcoH1PO4DmndBand_vr(:)         => null()   !HPO4 demand in band by all microbial,       root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: REcoNO3DmndSoil_vr(:)           => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: REcoNH4DmndSoil_vr(:)           => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: REcoNO3DmndBand_vr(:)           => null()   !total root + microbial NO3 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: REcoH2PO4DmndSoil_vr(:)         => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: REcoH2PO4DmndBand_vr(:)         => null()   !total root + microbial PO4 uptake band,     [gP d-2 h-1]
  real(r8), pointer :: RH2PO4EcoDmndSoilPrev_vr(:)     => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: RH2PO4EcoDmndBandPrev_vr(:)     => null()   !total root + microbial PO4 uptake band,     [gP d-2 h-1]
  real(r8), pointer :: RH1PO4EcoDmndSoilPrev_vr(:)     => null()   !HPO4 demand in non-band by all microbial,   root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: RH1PO4EcoDmndBandPrev_vr(:)     => null()   !HPO4 demand in band by all microbial,       root, myco populations, [gP d-2 h-1]
  real(r8), pointer :: RNO3EcoDmndSoilPrev_vr(:)       => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4EcoDmndSoilPrev_vr(:)       => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4EcoDmndBandPrev_vr(:)       => null()   !total root + microbial NH4 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: RNO3EcoDmndBandPrev_vr(:)       => null()   !total root + microbial NO3 uptake band,     [gN d-2 h-1]
  real(r8), pointer :: RO2GasXchangePrev_vr(:)         => null()   !net gaseous O2 flux,                        [g d-2 h-1]
  real(r8), pointer :: RCO2GasFlxPrev_vr(:)            => null()   !net gaseous CO2 flux,                       [g d-2 h-1]
  real(r8), pointer :: RO2AquaXchangePrev_vr(:)        => null()   !net aqueous O2 flux,                        [g d-2 h-1]
  real(r8), pointer :: RO2EcoDmndPrev_vr(:)             => null()   !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8), pointer :: tRootCO2Emis_vr(:)               => null()   !total root CO2 flux,              [gC d-2 h-1]
  real(r8), pointer :: tRO2MicrbUptk_vr(:)              => null()   !total root internal O2 flux,      [g d-2 h-1]
  real(r8), pointer :: LitrfalStrutElms_vr(:,:,:,:)     => null() !total LitrFall element,             [g d-2 h-1]
  real(r8), pointer :: REcoDOMProd_vr(:,:,:)            => null()  !net microbial DOC flux,            [gC d-2 h-1]
  real(r8), pointer :: CO2NetFix_pft(:)                 => null()  !canopy net CO2 exchange,           [gC d-2 h-1]
  real(r8), pointer :: GrossCO2Fix_pft(:)               => null()  !total gross CO2 fixation,          [gC d-2 ]
  real(r8), pointer :: LitrfalStrutElms_pft(:,:)        => null()  !plant element LitrFall,            [g d-2 h-1]
  real(r8), pointer :: RootGasLossDisturb_pft(:,:)      => null() !gaseous flux fron root disturbance, [g d-2 h-1]
  real(r8), pointer :: SurfLitrfalStrutElms_CumYr_pft(:,:) => null()  !total surface LitrFall element,    [g d-2]
  real(r8), pointer :: LitrfalStrutElms_CumYr_pft(:,:)     => null()  !total plant element LitrFall,      [g d-2 ]
  real(r8), pointer :: GrossResp_pft(:)                 => null()  !total plant respiration,           [gC d-2 ]
  real(r8), pointer :: CanopyGrosRCO2_pft(:)            => null()
  real(r8), pointer :: RootGrosRCO2_pft(:)              => null()
  real(r8), pointer :: NodulInfectElms_pft(:,:)         => null()
  real(r8), pointer :: NodulInfectElmsCum_pft(:,:)      => null()
  real(r8), pointer :: NH3Emis_CumYr_pft(:)                => null()  !total canopy NH3 flux,   [gN d-2 ]
  real(r8), pointer :: PlantN2Fix_CumYr_pft(:)             => null()  !total plant N2 fixation, [g d-2 ]

  contains
    procedure, public :: Init  => plt_bgcrate_init
    procedure, public :: Destroy  => plt_bgcrate_destroy
  end type plant_bgcrate_type

  type, public :: plant_rootbgc_type
  real(r8), pointer :: TRootGasLossDisturb_pft(:)        => null()  !total root gas content [g d-2]
  real(r8), pointer :: trcs_plant_uptake_vr(:,:)         => null()   !total root-soil solute flux non-band,                          [g d-2 h-1]
  real(r8), pointer :: RootMycoExudElms_pft(:,:)         => null()  !total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
  real(r8), pointer :: RootN2Fix_pft(:)                  => null()  !total root N2 fixation,                                         [g d-2 h-1]
  real(r8), pointer :: RootNO3Uptake_pft(:)              => null()  !total root uptake of NO3,                                       [g d-2 h-1]
  real(r8), pointer :: RootNH4Uptake_pft(:)              => null()  !total root uptake of NH4,                                       [g d-2 h-1]
  real(r8), pointer :: RootHPO4Uptake_pft(:)             => null()  !total root uptake of HPO4,                                      [g d-2 h-1]
  real(r8), pointer :: RootH2PO4Uptake_pft(:)            => null()  !total root uptake of PO4,                                       [g d-2 h-1]
  real(r8), pointer :: RootMycoExudElm_pvr(:,:,:,:,:)    => null()  !root uptake (+ve) - exudation (-ve) of DOE,                     [g d-2 h-1]
  real(r8), pointer :: PlantRootSoilElmNetX_pft(:,:)     => null()  !net root element uptake (+ve) - exudation (-ve),                [gC d-2 h-1]
  real(r8), pointer :: RO2UptkSoilM_vr(:,:)              => null()  !total O2 sink,                                                  [g d-2 t-1]
  real(r8), pointer :: ZERO4Uptk_pft(:)                  => null()  !threshold zero for uptake calculation
  real(r8), pointer :: CMinPO4Root_pft(:,:)              => null()  !minimum PO4 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxPO4Root_pft(:,:)              => null()  !maximum root PO4 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmPO4Root_pft(:,:)                => null()  !Km for root PO4 uptake,                                         [g m-3]
  real(r8), pointer :: CminNO3Root_pft(:,:)              => null()  !minimum NO3 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxNO3Root_pft(:,:)              => null()  !maximum root NO3 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmNO3Root_pft(:,:)                => null()  !Km for root NO3 uptake,                                         [g m-3]
  real(r8), pointer :: CMinNH4Root_pft(:,:)              => null()  !minimum NH4 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxNH4Root_pft(:,:)              => null()  !maximum root NH4 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmNH4Root_pft(:,:)                => null()  !Km for root NH4 uptake,                                         [g m-3]
  real(r8), pointer :: RootCO2Emis_pvr(:,:,:)            => null()  !aqueous CO2 flux from roots to root water,                      [g d-2 h-1]
  real(r8), pointer :: RootO2Uptk_pvr(:,:,:)             => null()  !aqueous O2 flux from roots to root water,                       [g d-2 h-1]
  real(r8), pointer :: RUPGasSol_vr(:,:,:,:)             => null()  !aqueous CO2 flux from roots to soil water,                      [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx__pvr(:,:,:,:)   => null()  !gaseous tracer flux through roots,                              [g d-2 h-1]
  real(r8), pointer :: trcg_Root_DisEvap_flx_vr(:,:,:,:) => null()  !dissolution (+ve) - volatilization (-ve) gas flux in roots,     [g d-2 h-1]
  real(r8), pointer :: RootO2Dmnd4Resp_pvr(:,:,:)        => null()  !root  O2 demand from respiration,                               [g d-2 h-1]
  real(r8), pointer :: RootNH4DmndSoil_pvr(:,:,:)        => null()  !root uptake of NH4 non-band unconstrained by NH4,               [g d-2 h-1]
  real(r8), pointer :: RootNH4DmndBand_pvr(:,:,:)        => null()  !root uptake of NO3 band unconstrained by NO3,                   [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndSoil_pvr(:,:,:)        => null()  !root uptake of NH4 band unconstrained by NH4,                   [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndBand_pvr(:,:,:)        => null()  !root uptake of NO3 non-band unconstrained by NO3,               [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndSoil_pvr(:,:,:)      => null()  !root uptake of H2PO4 non-band
  real(r8), pointer :: RootH2PO4DmndBand_pvr(:,:,:)      => null()  !root uptake of H2PO4 band
  real(r8), pointer :: RootH1PO4DmndSoil_pvr(:,:,:)      => null()  !HPO4 demand in non-band by each root population
  real(r8), pointer :: RootH1PO4DmndBand_pvr(:,:,:)      => null()  !HPO4 demand in band by each root population
  real(r8), pointer :: RAutoRootO2Limter_pvr(:,:,:)      => null()  !O2 constraint to root respiration,                              []
  real(r8), pointer :: trcg_rootml_pvr(:,:,:,:)          => null() !root gas content,                                                [g d-2]
  real(r8), pointer :: trcs_rootml_pvr(:,:,:,:)          => null() !root aqueous content,                                            [g d-2]
  real(r8), pointer :: NH3Dep2Can_brch(:,:)              => null()  !gaseous NH3 flux fron root disturbance band,                    [g d-2 h-1]
  real(r8), pointer :: RootNutUptake_pvr(:,:,:,:)        => null()  !root uptake of Nutrient band,                                   [g d-2 h-1]
  real(r8), pointer :: RootOUlmNutUptake_pvr(:,:,:,:)    => null()  !root uptake of NH4 band unconstrained by O2,                    [g d-2 h-1]
  real(r8), pointer :: RootCUlmNutUptake_pvr(:,:,:,:)    => null()  !root uptake of NH4 band unconstrained by root nonstructural C,  [g d-2 h-1]
  real(r8), pointer :: RootRespPotent_pvr(:,:,:)         => null()  !root respiration unconstrained by O2,                           [g d-2 h-1]
  real(r8), pointer :: RootCO2EmisPot_pvr(:,:,:)         => null()  !root CO2 efflux unconstrained by root nonstructural C,          [g d-2 h-1]
  real(r8), pointer :: RootCO2Autor_pvr(:,:,:)           => null()  !root respiration constrained by O2,                             [g d-2 h-1]
  real(r8), pointer :: PlantExudElm_CumYr_pft(:,:)      => null()  !total net root element uptake (+ve) - exudation (-ve),          [gC d-2 ]
  real(r8), pointer :: trcg_root_vr(:,:)                 => null()   !total root internal gas flux,                                  [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx_vr(:,:)         => null()   !total internal root gas flux,                                  [gC d-2 h-1]
  real(r8), pointer :: CO2FixCL_pft(:)                   => null()   !Rubisco-limited CO2 fixation
  real(r8), pointer :: CO2FixLL_pft(:)                   => null()   !Light-limited CO2 fixation
  real(r8), pointer :: RootUptk_N_CumYr_pft(:)           => null()
  real(r8), pointer :: RootUptk_P_CumYr_pft(:)           => null()
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
  allocate(this%trcs_plant_uptake_vr(ids_beg:ids_end,JZ1)); this%trcs_plant_uptake_vr=0._r8
  allocate(this%trcg_rootml_pvr(idg_beg:idg_end-1,2,JZ1,JP1));this%trcg_rootml_pvr=spval
  allocate(this%trcs_rootml_pvr(idg_beg:idg_end-1,2,JZ1,JP1));this%trcs_rootml_pvr=spval
  allocate(this%TRootGasLossDisturb_pft(idg_beg:idg_end-1));this%TRootGasLossDisturb_pft=spval
  allocate(this%RO2UptkSoilM_vr(60,0:JZ1)); this%RO2UptkSoilM_vr=spval
  allocate(this%RootMycoExudElm_pvr(NumPlantChemElms,2,1:jcplx,0:JZ1,JP1));this%RootMycoExudElm_pvr=spval
  allocate(this%PlantRootSoilElmNetX_pft(NumPlantChemElms,JP1)); this%PlantRootSoilElmNetX_pft=spval
  allocate(this%PlantExudElm_CumYr_pft(NumPlantChemElms,JP1));this%PlantExudElm_CumYr_pft=spval
  allocate(this%RootUptk_N_CumYr_pft(JP1)); this%RootUptk_N_CumYr_pft=spval
  allocate(this%RootUptk_P_CumYr_pft(JP1)); this%RootUptk_P_CumYr_pft=spval  
  allocate(this%RootMycoExudElms_pft(NumPlantChemElms,JP1));this%RootMycoExudElms_pft=spval
  allocate(this%RootN2Fix_pft(JP1)); this%RootN2Fix_pft=spval
  allocate(this%RootNO3Uptake_pft(JP1)); this%RootNO3Uptake_pft=spval
  allocate(this%RootNH4Uptake_pft(JP1)); this%RootNH4Uptake_pft=spval
  allocate(this%RootHPO4Uptake_pft(JP1)); this%RootHPO4Uptake_pft=spval
  allocate(this%RootH2PO4Uptake_pft(JP1)); this%RootH2PO4Uptake_pft=spval

  allocate(this%ZERO4Uptk_pft(JP1)); this%ZERO4Uptk_pft=spval
  allocate(this%RootRespPotent_pvr(jroots,JZ1,JP1)); this%RootRespPotent_pvr=spval
  allocate(this%RootCO2EmisPot_pvr(jroots,JZ1,JP1)); this%RootCO2EmisPot_pvr=spval
  allocate(this%RootCO2Autor_pvr(jroots,JZ1,JP1)); this%RootCO2Autor_pvr=spval
  allocate(this%RootNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1)); this%RootNutUptake_pvr=0._r8
  allocate(this%RootOUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1));this%RootOUlmNutUptake_pvr=spval
  allocate(this%RootCUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1));this%RootCUlmNutUptake_pvr=spval
  allocate(this%NH3Dep2Can_brch(MaxNumBranches,JP1));this%NH3Dep2Can_brch=spval
  allocate(this%CO2FixCL_pft(JP1)); this%CO2FixCL_pft=spval
  allocate(this%CO2FixLL_pft(JP1)); this%CO2FixLL_pft=spval
  allocate(this%trcg_air2root_flx_vr(idg_beg:idg_end-1,JZ1));this%trcg_air2root_flx_vr=spval
  allocate(this%trcg_root_vr(idg_beg:idg_end-1,JZ1));this%trcg_root_vr=spval

  allocate(this%trcg_air2root_flx__pvr(idg_beg:idg_end-1,2,JZ1,JP1));this%trcg_air2root_flx__pvr=spval
  allocate(this%trcg_Root_DisEvap_flx_vr(idg_beg:idg_end-1,2,JZ1,JP1));this%trcg_Root_DisEvap_flx_vr=spval
  allocate(this%RootO2Dmnd4Resp_pvr(jroots,JZ1,JP1));this%RootO2Dmnd4Resp_pvr=spval
  allocate(this%RootNH4DmndSoil_pvr(jroots,JZ1,JP1));this%RootNH4DmndSoil_pvr=spval
  allocate(this%RootNH4DmndBand_pvr(jroots,JZ1,JP1));this%RootNH4DmndBand_pvr=spval
  allocate(this%RootNO3DmndSoil_pvr(jroots,JZ1,JP1));this%RootNO3DmndSoil_pvr=spval
  allocate(this%RootNO3DmndBand_pvr(jroots,JZ1,JP1));this%RootNO3DmndBand_pvr=spval
  allocate(this%RootH2PO4DmndSoil_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndSoil_pvr=spval
  allocate(this%RootH2PO4DmndBand_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndBand_pvr=spval
  allocate(this%RootH1PO4DmndSoil_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndSoil_pvr=spval
  allocate(this%RootH1PO4DmndBand_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndBand_pvr=spval

  allocate(this%RAutoRootO2Limter_pvr(jroots,JZ1,JP1));this%RAutoRootO2Limter_pvr=spval
  allocate(this%CMinPO4Root_pft(jroots,JP1));this%CMinPO4Root_pft=spval
  allocate(this%VmaxPO4Root_pft(jroots,JP1));this%VmaxPO4Root_pft=spval
  allocate(this%KmPO4Root_pft(jroots,JP1));this%KmPO4Root_pft=spval
  allocate(this%CminNO3Root_pft(jroots,JP1));this%CminNO3Root_pft=spval
  allocate(this%VmaxNO3Root_pft(jroots,JP1));this%VmaxNO3Root_pft=spval
  allocate(this%KmNO3Root_pft(jroots,JP1));this%KmNO3Root_pft=spval
  allocate(this%CMinNH4Root_pft(jroots,JP1));this%CMinNH4Root_pft=spval
  allocate(this%VmaxNH4Root_pft(jroots,JP1));this%VmaxNH4Root_pft=spval
  allocate(this%KmNH4Root_pft(jroots,JP1));this%KmNH4Root_pft=spval
  allocate(this%RootCO2Emis_pvr(jroots,JZ1,JP1));this%RootCO2Emis_pvr=spval
  allocate(this%RootO2Uptk_pvr(jroots,JZ1,JP1));this%RootO2Uptk_pvr=spval
  allocate(this%RUPGasSol_vr(idg_beg:idg_end,jroots,JZ1,JP1));this%RUPGasSol_vr=spval
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
!  if(allocated(RO2UptkSoilM_vr))deallocate(RO2UptkSoilM_vr)
!  if(allocated(RootMycoExudElm_pvr))deallocate(RootMycoExudElm_pvr)
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
!  if(allocated(RAutoRootO2Limter_pvr))deallocate(RAutoRootO2Limter_pvr)
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
  allocate(this%AtmGasc(idg_beg:idg_end-1));this%AtmGasc=spval
  allocate(this%DATAP(JP1)); this%DATAP=''
  allocate(this%DATA(30)); this%DATA=''
  allocate(this%AREA3(0:JZ1));this%AREA3=spval
  allocate(this%DLYR3(0:JZ1)); this%DLYR3=spval
  allocate(this%ElmBalanceCum_pft(NumPlantChemElms,JP1));this%ElmBalanceCum_pft=spval
  allocate(this%CumSoilThickness_vr(0:JZ1));this%CumSoilThickness_vr=spval
  allocate(this%DPTHZ_vr(0:JZ1));this%DPTHZ_vr=spval
  allocate(this%PPI_pft(JP1));this%PPI_pft=spval
  allocate(this%PPatSeeding_pft(JP1));this%PPatSeeding_pft=spval
  allocate(this%PPX_pft(JP1));this%PPX_pft=spval
  allocate(this%PlantPopulation_pft(JP1));this%PlantPopulation_pft=spval
  allocate(this%flag_pft_active(JP1));  this%flag_pft_active=.false.
  allocate(this%VLWatMicPM_vr(60,0:JZ1));this%VLWatMicPM_vr=spval
  allocate(this%VLsoiAirPM(60,0:JZ1));this%VLsoiAirPM=spval
  allocate(this%TortMicPM_vr(60,0:JZ1));this%TortMicPM_vr=spval
  allocate(this%FILM(60,0:JZ1)); this%FILM=spval

  end subroutine plt_site_Init
!----------------------------------------------------------------------
  subroutine plt_site_destroy(this)
  implicit none
  class(plant_siteinfo_type) :: this

!  if(allocated(DLYR3))deallocate(DLYR3)
!  if(allocated(FracSoiAsMicP_vr))deallocate(FracSoiAsMicP_vr)
!  if(allocated(TortMicPM_vr))deallocate(TortMicPM_vr)
!  if(allocated(FILM))deallocate(FILM)
!  if(allocated(DATAP))deallocate(DATAP)
!  if(allocated(DATA))deallocate(DATA)
!  call Destroy(this%PlantElemntStoreLandscape)
!  if(allocated(AREA3))deallocate(AREA3)
!  if(allocated(ElmBalanceCum_pft))deallocate(ElmBalanceCum_pft)
!  if(allocated(BALP))deallocate(BALP)
!  if(allocated(CumSoilThickness_vr))deallocate(CumSoilThickness_vr)
!  if(allocated(DPTHZ))deallocate(DPTHZ)
!  if(allocated(PPI_pft)) deallocate(PPI_pft)
!  if(allocated(PPatSeeding_pft)) deallocate(PPatSeeding_pft)
!  if(allocated(PPX)) deallocate(PPX)
!  if(allocated(PP))deallocate(PP)
!  if(allocated(VLWatMicPM))deallocate(VLWatMicPM)
!  if(allocated(VLsoiAirPM))deallocate(VLsoiAirPM)

  end subroutine plt_site_destroy
!----------------------------------------------------------------------

  subroutine plt_bgcrate_init(this)
  implicit none
  class(plant_bgcrate_type) :: this

  allocate(this%tRootCO2Emis_vr(JZ1)); this%tRootCO2Emis_vr=spval
  allocate(this%tRO2MicrbUptk_vr(JZ1)); this%tRO2MicrbUptk_vr=spval
  allocate(this%RH2PO4EcoDmndSoilPrev_vr(0:JZ1)); this%RH2PO4EcoDmndSoilPrev_vr=spval
  allocate(this%RH2PO4EcoDmndBandPrev_vr(0:JZ1)); this%RH2PO4EcoDmndBandPrev_vr=spval
  allocate(this%RH1PO4EcoDmndSoilPrev_vr(0:JZ1)); this%RH1PO4EcoDmndSoilPrev_vr=spval
  allocate(this%RH1PO4EcoDmndBandPrev_vr(0:JZ1)); this%RH1PO4EcoDmndBandPrev_vr=spval
  allocate(this%RNO3EcoDmndSoilPrev_vr(0:JZ1)); this%RNO3EcoDmndSoilPrev_vr=spval
  allocate(this%RNH4EcoDmndSoilPrev_vr(0:JZ1)); this%RNH4EcoDmndSoilPrev_vr =spval
  allocate(this%RNH4EcoDmndBandPrev_vr(0:JZ1)); this%RNH4EcoDmndBandPrev_vr=spval
  allocate(this%RNO3EcoDmndBandPrev_vr(0:JZ1)); this%RNO3EcoDmndBandPrev_vr=spval
  allocate(this%RO2GasXchangePrev_vr(0:JZ1)); this%RO2GasXchangePrev_vr=spval
  allocate(this%RCO2GasFlxPrev_vr(0:JZ1)); this%RCO2GasFlxPrev_vr=spval
  allocate(this%RO2AquaXchangePrev_vr(0:JZ1)); this%RO2AquaXchangePrev_vr=spval
  allocate(this%RO2EcoDmndPrev_vr(0:JZ1)); this%RO2EcoDmndPrev_vr=spval
  allocate(this%LitrfalStrutElms_vr(NumPlantChemElms,jsken,NumOfPlantLitrCmplxs,0:JZ1));this%LitrfalStrutElms_vr=spval
  allocate(this%GrossCO2Fix_pft(JP1));this%GrossCO2Fix_pft=spval
  allocate(this%REcoDOMProd_vr(1:NumPlantChemElms,1:jcplx,0:JZ1));this%REcoDOMProd_vr=spval
  allocate(this%CO2NetFix_pft(JP1));this%CO2NetFix_pft=spval
  allocate(this%RootGasLossDisturb_pft(idg_beg:idg_end-1,JP1));this%RootGasLossDisturb_pft=spval
  allocate(this%GrossResp_pft(JP1));this%GrossResp_pft=spval
  allocate(this%CanopyGrosRCO2_pft(JP1));this%CanopyGrosRCO2_pft=spval
  allocate(this%RootGrosRCO2_pft(JP1));this%RootGrosRCO2_pft=spval
  allocate(this%PlantN2Fix_CumYr_pft(JP1));this%PlantN2Fix_CumYr_pft=spval
  allocate(this%NH3Emis_CumYr_pft(JP1));this%NH3Emis_CumYr_pft=spval
  allocate(this%NodulInfectElms_pft(NumPlantChemElms,JP1));this%NodulInfectElms_pft=spval
  allocate(this%NodulInfectElmsCum_pft(NumPlantChemElms,JP1));this%NodulInfectElmsCum_pft=spval
  allocate(this%SurfLitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP1));this%SurfLitrfalStrutElms_CumYr_pft=spval
  allocate(this%LitrFallStrutElms_col(NumPlantChemElms));this%LitrFallStrutElms_col=spval
  allocate(this%NetPrimProduct_pft(JP1));this%NetPrimProduct_pft=spval
  allocate(this%NH3Dep2Can_pft(JP1));this%NH3Dep2Can_pft=spval
  allocate(this%tRootMycoExud2Soil_vr(NumPlantChemElms,1:jcplx,JZ1));this%tRootMycoExud2Soil_vr=spval
  allocate(this%RootN2Fix_pvr(JZ1,JP1));this%RootN2Fix_pvr=spval
  allocate(this%CanopyRespC_CumYr_pft(JP1));this%CanopyRespC_CumYr_pft=spval
  allocate(this%REcoH1PO4DmndBand_vr(0:JZ1));this%REcoH1PO4DmndBand_vr=spval
  allocate(this%REcoNO3DmndSoil_vr(0:JZ1));this%REcoNO3DmndSoil_vr=spval
  allocate(this%REcoH2PO4DmndSoil_vr(0:JZ1));this%REcoH2PO4DmndSoil_vr=spval
  allocate(this%REcoNH4DmndSoil_vr(0:JZ1));this%REcoNH4DmndSoil_vr=spval
  allocate(this%REcoO2DmndResp_vr(0:JZ1));this%REcoO2DmndResp_vr=spval
  allocate(this%REcoH2PO4DmndBand_vr(0:JZ1));this%REcoH2PO4DmndBand_vr=spval
  allocate(this%REcoNO3DmndBand_vr(0:JZ1));this%REcoNO3DmndBand_vr=spval
  allocate(this%REcoNH4DmndBand_vr(0:JZ1));this%REcoNH4DmndBand_vr=spval
  allocate(this%REcoH1PO4DmndSoil_vr(0:JZ1));this%REcoH1PO4DmndSoil_vr=spval

  allocate(this%LitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP1));this%LitrfalStrutElms_CumYr_pft=spval
  allocate(this%LitrfalStrutElms_pft(NumPlantChemElms,JP1));this%LitrfalStrutElms_pft=spval
  allocate(this%LitrfalStrutElms_pvr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ1,JP1))
  this%LitrfalStrutElms_pvr=spval

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
  allocate(this%N2ObyFire_CumYr_pft(JP1));this%N2ObyFire_CumYr_pft=spval
  allocate(this%NH3byFire_CumYr_pft(JP1));this%NH3byFire_CumYr_pft=spval
  allocate(this%PO4byFire_CumYr_pft(JP1));this%PO4byFire_CumYr_pft=spval

  allocate(this%iHarvestDay_pft(JP1)); this%iHarvestDay_pft=0
  allocate(this%O2ByFire_CumYr_pft(JP1));this%O2ByFire_CumYr_pft=spval
  allocate(this%FracBiomHarvsted(1:2,1:4,JP1));this%FracBiomHarvsted=spval
  allocate(this%FracCanopyHeightCut_pft(JP1)); this%FracCanopyHeightCut_pft=spval
  allocate(this%iYearPlantHarvest_pft(JP1));this%iYearPlantHarvest_pft=0
  allocate(this%FERT(1:20));this%FERT=spval
  allocate(this%iYearPlanting_pft(JP1));this%iYearPlanting_pft=0
  allocate(this%iPlantingDay_pft(JP1));this%iPlantingDay_pft=0
  allocate(this%iPlantingYear_pft(JP1));this%iPlantingYear_pft=0
  allocate(this%iHarvestYear_pft(JP1));this%iHarvestYear_pft=0
  allocate(this%iDayPlanting_pft(JP1));this%iDayPlanting_pft=0
  allocate(this%iDayPlantHarvest_pft(JP1));this%iDayPlantHarvest_pft=0
  allocate(this%iHarvstType_pft(JP1));this%iHarvstType_pft=-1
  allocate(this%jHarvst_pft(JP1));this%jHarvst_pft=0

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
!  if(allocated(jHarvst_pft))deallocate(jHarvst_pft)

  end subroutine plt_disturb_destroy
!----------------------------------------------------------------------
  subroutine  plt_ew_init(this)

  implicit none
  class(plant_ew_type) :: this

  allocate(this%QdewCanopy_pft(JP1)); this%QdewCanopy_pft=spval
  allocate(this%ETCanopy_CumYr_pft(JP1));this%ETCanopy_CumYr_pft=spval
  allocate(this%TotalSoilH2OPSIMPa_vr(0:JZ1));this%TotalSoilH2OPSIMPa_vr=spval
  allocate(this%THeatRootUptake_vr(0:JZ1));this%THeatRootUptake_vr=spval
  allocate(this%TKCanopy_pft(JP1));this%TKCanopy_pft=spval
  allocate(this%HeatXAir2PCan_pft(JP1));this%HeatXAir2PCan_pft=spval
  allocate(this%PrecIntcptByCanopy_pft(JP1));this%PrecIntcptByCanopy_pft=spval
  allocate(this%PSICanopyTurg_pft(JP1));this%PSICanopyTurg_pft=spval
  allocate(this%PSICanopy_pft(JP1));this%PSICanopy_pft=spval
  allocate(this%VapXAir2Canopy_pft(JP1));this%VapXAir2Canopy_pft=spval
  allocate(this%HeatStorCanopy_pft(JP1));this%HeatStorCanopy_pft=spval
  allocate(this%EvapTransHeat_pft(JP1));this%EvapTransHeat_pft=spval
  allocate(this%WatByPCanopy_pft(JP1));this%WatByPCanopy_pft=spval
  allocate(this%VHeatCapCanP_pft(JP1));this%VHeatCapCanP_pft=spval
  allocate(this%CanopyWater_pft(JP1));this%CanopyWater_pft=spval
  allocate(this%PSIRoot_pvr(jroots,JZ1,JP1));this%PSIRoot_pvr=spval
  allocate(this%PSIRootOSMO_vr(jroots,JZ1,JP1));this%PSIRootOSMO_vr=spval
  allocate(this%PSIRootTurg_vr(jroots,JZ1,JP1));this%PSIRootTurg_vr=spval
  allocate(this%AllPlantRootH2OUptake_vr(jroots,JZ1,JP1));this%AllPlantRootH2OUptake_vr=spval
  allocate(this%TPlantRootH2OUptake_vr(0:JZ1));this%TPlantRootH2OUptake_vr=spval
  allocate(this%Transpiration_pft(JP1));this%Transpiration_pft=spval
  allocate(this%PSICanopyOsmo_pft(JP1));this%PSICanopyOsmo_pft=spval
  allocate(this%TKS_vr(0:JZ1));this%TKS_vr=spval
  allocate(this%CanOsmoPsi0pt_pft(JP1));this%CanOsmoPsi0pt_pft=spval
  allocate(this%RAZ(JP1));this%RAZ=spval
  allocate(this%DeltaTKC_pft(JP1));this%DeltaTKC_pft=spval
  allocate(this%TKC(JP1));this%TKC=spval
  allocate(this%ENGYX_pft(JP1));this%ENGYX_pft=spval
  allocate(this%TCelciusCanopy_pft(JP1));this%TCelciusCanopy_pft=spval
  allocate(this%PSICanPDailyMin(JP1));this%PSICanPDailyMin=spval

  end subroutine plt_ew_init
!----------------------------------------------------------------------

  subroutine plt_ew_destroy(this)
  implicit none
  class(plant_ew_type) :: this

!  if(allocated(ETCanopy_CumYr_pft))deallocate(ETCanopy_CumYr_pft)
!  if(allocated(PSIST))deallocate(PSIST)
!  if(allocated(THeatRootUptake_vr))deallocate(THeatRootUptake_vr)
!  if(allocated(TKCanopy_pft))deallocate(TKCanopy_pft)
!  if(allocated(HeatXAir2PCan_pft))deallocate(HeatXAir2PCan_pft)
!  if(allocated(PrecIntcptByCanopy_pft))deallocate(PrecIntcptByCanopy_pft)
!  if(allocated(PSICanopyTurg_pft))deallocate(PSICanopyTurg_pft)
!  if(allocated(PSICanopy_pft))deallocate(PSICanopy_pft)
!  if(allocated(VapXAir2Canopy_pft))deallocate(VapXAir2Canopy_pft)
!  if(allocated(HeatStorCanopy_pft))deallocate(HeatStorCanopy_pft)
!  if(allocated(EvapTransHeat_pft))deallocate(EvapTransHeat_pft)
!  if(allocated(WatByPCanopy_pft))deallocate(WatByPCanopy_pft)
!  if(allocated(VHeatCapCanP_pft))deallocate(VHeatCapCanP_pft)
!  if(allocated(CanopyWater_pft))deallocate(CanopyWater_pft)
!  if(allocated(PSIRoot_pvr))deallocate(PSIRoot_pvr)
!  if(allocated(PSIRootOSMO_vr))deallocate(PSIRootOSMO_vr)
!  if(allocated(PSIRootTurg_vr))deallocate(PSIRootTurg_vr)
!  if(allocated(AllPlantRootH2OUptake_vr))deallocate(AllPlantRootH2OUptake_vr)
!  if(allocated(TAllPlantRootH2OUptake_vr))deallocate(TAllPlantRootH2OUptake_vr)
!  if(allocated(Transpiration_pft))deallocate(Transpiration_pft)
!  if(allocated(PSICanopyOsmo_pft))deallocate(PSICanopyOsmo_pft)
!  if(allocated(TKS))deallocate(TKS)
!  if(allocated(PSICanPDailyMin))deallocate(PSICanPDailyMin)
!  if(allocated(RAZ))deallocate(RAZ)
!  if(allocated(CanOsmoPsi0pt_pft))deallocate(CanOsmoPsi0pt_pft)
!  if(allocated(DeltaTKC_pft))deallocate(DeltaTKC_pft)
!  if(allocated(TKC))deallocate(TKC)
!  if(allocated(ENGYX_pft))deallocate(ENGYX_pft)
!  if(allocated(TCelciusCanopy_pft))deallocate(TCelciusCanopy_pft)

  end subroutine plt_ew_destroy

!----------------------------------------------------------------------

  subroutine plt_allom_init(this)
  implicit none
  class(plant_allometry_type) :: this


  allocate(this%RootFracRemobilizableBiom(JP1));this%RootFracRemobilizableBiom=spval
  allocate(this%FracGroth2Node_pft(JP1));this%FracGroth2Node_pft=spval
  allocate(this%GrainSeedBiomCMean_brch(MaxNumBranches,JP1));this%GrainSeedBiomCMean_brch=spval
  allocate(this%NoduGrowthYield_pft(JP1));this%NoduGrowthYield_pft=spval
  allocate(this%RootBiomGrosYld_pft(JP1));this%RootBiomGrosYld_pft=spval
  allocate(this%RootrPC_pft(JP1));this%RootrPC_pft=spval
  allocate(this%rCNNonstRemob_pft(JP1));this%rCNNonstRemob_pft=spval
  allocate(this%rCPNonstRemob_pft(JP1));this%rCPNonstRemob_pft=spval
  allocate(this%CPRTS_pft(JP1));this%CPRTS_pft=spval
  allocate(this%CNRTS_pft(JP1));this%CNRTS_pft=spval
  allocate(this%NodulerNC_pft(JP1));this%NodulerNC_pft=spval
  allocate(this%NodulerPC_pft(JP1));this%NodulerPC_pft=spval
  allocate(this%RootrNC_pft(JP1));this%RootrNC_pft=spval
  allocate(this%CPLF(JP1));this%CPLF=spval
  allocate(this%CPSHE(JP1));this%CPSHE=spval
  allocate(this%CNLF(JP1));this%CNLF=spval
  allocate(this%CNSHE(JP1));this%CNSHE=spval
  allocate(this%CNGR(JP1));this%CNGR=spval
  allocate(this%rPCStalk_pft(JP1));this%rPCStalk_pft=spval
  allocate(this%rNCStalk_pft(JP1));this%rNCStalk_pft=spval
  allocate(this%CPGR(JP1));this%CPGR=spval
  allocate(this%rPCEar_pft(JP1));this%rPCEar_pft=spval
  allocate(this%rPCReserve_pft(JP1));this%rPCReserve_pft=spval
  allocate(this%rNCReserve_pft(JP1));this%rNCReserve_pft=spval
  allocate(this%rPCHusk_pft(JP1));this%rPCHusk_pft=spval
  allocate(this%FracHour4LeafoffRemob(0:5));this%FracHour4LeafoffRemob=spval
  allocate(this%FracRootStalkElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracRootStalkElmAlloc2Litr=spval
  allocate(this%FracShootStalkElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracShootStalkElmAlloc2Litr=spval
  allocate(this%FracRootElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracRootElmAlloc2Litr=spval
  allocate(this%FracShootLeafElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracShootLeafElmAlloc2Litr=spval

  allocate(this%PetioleBiomGrowthYield(JP1));this%PetioleBiomGrowthYield=spval
  allocate(this%HuskBiomGrowthYield(JP1));this%HuskBiomGrowthYield=spval
  allocate(this%StalkBiomGrowthYield(JP1));this%StalkBiomGrowthYield=spval
  allocate(this%ReserveBiomGrowthYield(JP1));this%ReserveBiomGrowthYield=spval
  allocate(this%EarBiomGrowthYield(JP1));this%EarBiomGrowthYield=spval
  allocate(this%GrainBiomGrowthYield(JP1));this%GrainBiomGrowthYield=spval
  allocate(this%rNCHusk_pft(JP1));this%rNCHusk_pft=spval
  allocate(this%rNCEar_pft(JP1));this%rNCEar_pft=spval
  allocate(this%LeafBiomGrowthYield(JP1));this%LeafBiomGrowthYield=spval

  end subroutine plt_allom_init
!----------------------------------------------------------------------

  subroutine plt_allom_destroy(this)
  implicit none

  class(plant_allometry_type) :: this

!  if(allocated(FracGroth2Node_pft))deallocate(FracGroth2Node_pft)
!  if(allocated(GrainSeedBiomCMean_brch))deallocate(GrainSeedBiomCMean_brch)
!  if(allocated(NoduGrowthYield_pft))deallocate(NoduGrowthYield_pft)
!  if(allocated(RootBiomGrosYld_pft))deallocate(RootBiomGrosYld_pft)
!  if(allocated(RootrPC_pft))deallocate(RootrPC_pft)
!  if(allocated(rCNNonstRemob_pft))deallocate(rCNNonstRemob_pft)
!  if(allocated(rCPNonstRemob_pft))deallocate(rCPNonstRemob_pft)
!  if(allocated(CPRTS_pft))deallocate(CPRTS_pft)
!  if(allocated(CNRTS_pft))deallocate(CNRTS_pft)
!  if(allocated(NodulerNC_pft))deallocate(NodulerNC_pft)
!  if(allocated(NodulerPC_pft))deallocate(NodulerPC_pft)
!  if(allocated(RootrNC_pft))deallocate(RootrNC_pft)
!  if(allocated(RootFracRemobilizableBiom))deallocate(RootFracRemobilizableBiom)
!  if(allocated(CPLF))deallocate(CPLF)
!  if(allocated(CPSHE))deallocate(CPSHE)
!  if(allocated(CNSHE))deallocate(CNSHE)
!  if(allocated(CNLF))deallocate(CNLF)
!  if(allocated(LeafBiomGrowthYield))deallocate(LeafBiomGrowthYield)
!  if(allocated(rPCStalk_pft))deallocate(rPCStalk_pft)
!  if(allocated(CNGR))deallocate(CNGR)
!  if(allocated(rNCStalk_pft))deallocate(rNCStalk_pft)
!  if(allocated(CPGR))deallocate(CPGR)
!  if(allocated(PetioleBiomGrowthYield))deallocate(PetioleBiomGrowthYield)
!  if(allocated(StalkBiomGrowthYield))deallocate(StalkBiomGrowthYield)
!  if(allocated(GrainBiomGrowthYield))deallocate(GrainBiomGrowthYield)
!  if(allocated(ReserveBiomGrowthYield))deallocate(ReserveBiomGrowthYield)
!  if(allocated(EarBiomGrowthYield))deallocate(EarBiomGrowthYield)
!  if(allocated(HuskBiomGrowthYield))deallocate(HuskBiomGrowthYield)
!  if(allocated(FracHour4LeafoffRemob))deallocate(FracHour4LeafoffRemob)
!  if(allocated(FracShootStalkElmAlloc2Litr))deallocate(FracShootStalkElmAlloc2Litr)
!  if(allocated(FracRootStalkElmAlloc2Litr))deallocate(FracRootStalkElmAlloc2Litr)
!  if(allocated(FracRootElmAlloc2Litr))deallocate(FracRootElmAlloc2Litr)
!  if(allocated(FracShootLeafElmAlloc2Litr))deallocate(FracShootLeafElmAlloc2Litr)
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
  allocate(this%RootNodulStrutElms_pvr(NumPlantChemElms,JZ1,JP1));this%RootNodulStrutElms_pvr=spval
  allocate(this%CanopyLeafCLyr_pft(NumOfCanopyLayers1,JP1));this%CanopyLeafCLyr_pft=spval
  allocate(this%RootNodulNonstElms_pvr(NumPlantChemElms,JZ1,JP1));this%RootNodulNonstElms_pvr=spval
  allocate(this%StandDeadKCompElms_pft(NumPlantChemElms,jsken,JP1));this%StandDeadKCompElms_pft=0._r8
  allocate(this%RootMyco2ndStrutElms_rpvr(NumPlantChemElms,jroots,JZ1,NumOfCanopyLayers1,JP1))
  this%RootMyco2ndStrutElms_rpvr=spval
  allocate(this%RootMycoNonstElms_pft(NumPlantChemElms,jroots,JP1));this%RootMycoNonstElms_pft=spval
  allocate(this%RootMyco1stStrutElms_rpvr(NumPlantChemElms,jroots,JZ1,NumOfCanopyLayers1,JP1))
  this%RootMyco1stStrutElms_rpvr=spval
  allocate(this%RootNodulElms_pft(NumPlantChemElms,JP1));this%RootNodulElms_pft=spval
  allocate(this%CanopyNonstElmConc_pft(NumPlantChemElms,JP1));this%CanopyNonstElmConc_pft=spval
  allocate(this%CanopyNonstElms_pft(NumPlantChemElms,JP1));this%CanopyNonstElms_pft=spval
  allocate(this%CanopyNodulNonstElms_pft(NumPlantChemElms,JP1));this%CanopyNodulNonstElms_pft=spval
  allocate(this%NoduleNonstructCconc_pft(JP1));this%NoduleNonstructCconc_pft=spval
  allocate(this%RootProteinConc_pvr(jroots,JZ1,JP1));this%RootProteinConc_pvr=spval
  allocate(this%RootProteinC_pvr(jroots,JZ1,JP1));this%RootProteinC_pvr=spval
  allocate(this%RootMycoActiveBiomC_pvr(jroots,JZ1,JP1));this%RootMycoActiveBiomC_pvr=spval
  allocate(this%RootMassElm_pvr(NumPlantChemElms,JZ1,JP1)); this%RootMassElm_pvr = 0._r8
  allocate(this%PopuRootMycoC_pvr(jroots,JZ1,JP1));this%PopuRootMycoC_pvr=spval
  allocate(this%RootMycoNonstElms_rpvr(NumPlantChemElms,jroots,JZ1,JP1));this%RootMycoNonstElms_rpvr=spval
  allocate(this%RootNonstructElmConc_pvr(NumPlantChemElms,jroots,JZ1,JP1));this%RootNonstructElmConc_pvr=spval
  allocate(this%CanopyNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNonstElms_brch=spval
  allocate(this%ShootElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%ShootElms_brch=spval
  allocate(this%ShootC4NonstC_brch(MaxNumBranches,JP1));this%ShootC4NonstC_brch=spval
  allocate(this%CanopyStalkC_pft(JP1));this%CanopyStalkC_pft=spval
  allocate(this%ShootElms_pft(NumPlantChemElms,JP1));this%ShootElms_pft=spval
  allocate(this%ShootElmsBeg_pft(NumPlantChemElms,JP1));this%ShootElmsBeg_pft=spval
  allocate(this%LeafProteinCNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafProteinCNode_brch=spval
  allocate(this%PetoleProteinCNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%PetoleProteinCNode_brch=spval
  allocate(this%InternodeStrutElms_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%InternodeStrutElms_brch=spval
  allocate(this%LeafElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafElmntNode_brch=spval
  allocate(this%PetioleElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%PetioleElmntNode_brch=spval
  allocate(this%LeafElmsByLayerNode_brch(NumPlantChemElms,NumOfCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafElmsByLayerNode_brch=spval
  allocate(this%StalkBiomassC_brch(MaxNumBranches,JP1));this%StalkBiomassC_brch=spval
  allocate(this%CanopyNodulNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNodulNonstElms_brch=spval
  allocate(this%LeafPetoNonstElmConc_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafPetoNonstElmConc_brch=spval
  allocate(this%RootStrutElms_pft(NumPlantChemElms,JP1));this%RootStrutElms_pft=spval
  allocate(this%CanopyNodulElms_pft(NumPlantChemElms,JP1));this%CanopyNodulElms_pft=spval
  allocate(this%StandDeadStrutElmsBeg_pft(NumPlantChemElms,JP1));this%StandDeadStrutElmsBeg_pft=spval
  allocate(this%tCanLeafC_cl(NumOfCanopyLayers1));this%tCanLeafC_cl=spval
  allocate(this%RootElms_pft(NumPlantChemElms,JP1));this%RootElms_pft=spval
  allocate(this%RootElmsBeg_pft(NumPlantChemElms,JP1));this%RootElmsBeg_pft=spval
  allocate(this%SeedCPlanted_pft(JP1));this%SeedCPlanted_pft=spval
  allocate(this%SeasonalNonstElms_pft(NumPlantChemElms,JP1));this%SeasonalNonstElms_pft=spval
  allocate(this%CanopyLeafShethC_pft(JP1));this%CanopyLeafShethC_pft=spval
  allocate(this%StandDeadStrutElms_pft(NumPlantChemElms,JP1));this%StandDeadStrutElms_pft=spval
  allocate(this%RootBiomCPerPlant_pft(JP1));this%RootBiomCPerPlant_pft=spval
  allocate(this%PetoleStrutElms_pft(NumPlantChemElms,JP1));this%PetoleStrutElms_pft=spval
  allocate(this%StalkStrutElms_pft(NumPlantChemElms,JP1));this%StalkStrutElms_pft=spval
  allocate(this%StalkRsrvElms_pft(NumPlantChemElms,JP1));this%StalkRsrvElms_pft=spval
  allocate(this%GrainStrutElms_pft(NumPlantChemElms,JP1));this%GrainStrutElms_pft=spval
  allocate(this%HuskStrutElms_pft(NumPlantChemElms,JP1));this%HuskStrutElms_pft=spval
  allocate(this%EarStrutElms_pft(NumPlantChemElms,JP1));this%EarStrutElms_pft=spval
  allocate(this%NodulStrutElms_pft(NumPlantChemElms,JP1));this%NodulStrutElms_pft=spval
  allocate(this%LeafChemElmRemob_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafChemElmRemob_brch=spval
  allocate(this%PetioleChemElmRemob_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetioleChemElmRemob_brch=spval
  allocate(this%ShootStrutElms_pft(NumPlantChemElms,JP1));this%ShootStrutElms_pft=spval
  allocate(this%AvgCanopyBiomC2Graze_pft(JP1));this%AvgCanopyBiomC2Graze_pft=spval
  allocate(this%LeafStrutElms_pft(NumPlantChemElms,JP1));this%LeafStrutElms_pft=spval
  allocate(this%StandingDeadInitC_pft(JP1));this%StandingDeadInitC_pft=spval
  allocate(this%LeafPetolBiomassC_brch(MaxNumBranches,JP1));this%LeafPetolBiomassC_brch=spval
  allocate(this%StalkRsrvElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%StalkRsrvElms_brch=spval
  allocate(this%LeafStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafStrutElms_brch=spval
  allocate(this%CanopyNodulStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNodulStrutElms_brch=spval
  allocate(this%PetoleStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetoleStrutElms_brch=spval
  allocate(this%EarStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%EarStrutElms_brch=spval
  allocate(this%HuskStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%HuskStrutElms_brch=spval
  allocate(this%GrainStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%GrainStrutElms_brch=spval
  allocate(this%StalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%StalkStrutElms_brch=spval
  allocate(this%ShootStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%ShootStrutElms_brch=spval
  allocate(this%SenecStalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%SenecStalkStrutElms_brch=spval
  allocate(this%RootMyco1stElm_raxs(NumPlantChemElms,jroots,MaxNumRootAxes,JP1));this%RootMyco1stElm_raxs=spval

  end subroutine plt_biom_init
!----------------------------------------------------------------------

  subroutine plt_biom_destroy(this)

  implicit none
  class(plant_biom_type) :: this

!  if(allocated(ZERO4LeafVar_pft))deallocate(ZERO4LeafVar_pft)
!  if(allocated(ZERO4Groth_pft))deallocate(ZERO4Groth_pft)
!  if(allocated(RootNodulNonstElms_pvr))deallocate(RootNodulNonstElms_pvr)
!  if(allocated(RootNodulStrutElms_pvr))deallocate(RootNodulStrutElms_pvr)
!  if(allocated(CanopyLeafCLyr_pft))deallocate(CanopyLeafCLyr_pft)
!  call destroy(RootMyco1stElm_raxs)
!  if(allocated(RootMyco1stStrutElms_rpvr))deallocate(RootMyco1stStrutElms_rpvr)
!  if(allocated(StandDeadKCompElms_pft))deallocate(StandDeadKCompElms_pft)
!  if(allocated(RootMyco2ndStrutElms_rpvr))deallocate(RootMyco2ndStrutElms_rpvr)
!  if(allocated(CanopyNonstElmConc_pft))deallocate(CanopyNonstElmConc_pft)
!  if(allocated(CanopyNonstElms_pft))deallocate(CanopyNonstElms_pft)
!  if(allocated(CanopyNodulNonstElms_pft))deallocate(CanopyNodulNonstElms_pft)
!  if(allocated(NoduleNonstructCconc_pft))deallocate(NoduleNonstructCconc_pft)
!  if(allocated(RootProteinConc_pvr))deallocate(RootProteinConc_pvr)
!  if(allocated(RootProteinC_pvr))deallocate(RootProteinC_pvr)
!  if(allocated(RootMycoActiveBiomC_pvr))deallocate(RootMycoActiveBiomC_pvr)
!  if(allocated( PopuRootMycoC_pvr))deallocate( PopuRootMycoC_pvr)
!  if(allocated(RootMycoNonstElms_rpvr))deallocate(RootMycoNonstElms_rpvr)
!  if(allocated(WVSTK))deallocate(WVSTK)
!  if(allocated(LeafElmntNode_brch))deallocate(LeafElmntNode_brch)
!  if(allocated(LeafProteinCNode_brch))deallocate(LeafProteinCNode_brch)
!  if(allocated(PetoleProteinCNode_brch))deallocate(PetoleProteinCNode_brch)
!  if(allocated(LeafElmsByLayerNode_brch))deallocate(LeafElmsByLayerNode_brch)
!  if(allocated(InternodeStrutElms_brch))deallocate(InternodeStrutElms_brch)
!  if(allocated(PetioleElmntNode_brch))deallocate(PetioleElmntNode_brch)
!  if(allocated(StalkBiomassC_brch))deallocate(StalkBiomassC_brch)
!  if(allocated(PPOOL))deallocate(PPOOL)
!  if(allocated(CanopyNodulNonstElms_brch))deallocate(CanopyNodulNonstElms_brch)
!  if(allocated(LeafPetoNonstElmConc_brch))deallocate(LeafPetoNonstElmConc_brch)
!  if(allocated(RootStrutElms_pft))deallocate(RootStrutElms_pft)
!  if(allocated(tCanLeafC_cl))deallocate(tCanLeafC_cl)
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
!  if(allocated(LeafPetolBiomassC_brch))deallocate(LeafPetolBiomassC_brch)
!  if(allocated(StalkRsrvElms_brch))deallocate(StalkRsrvElms_brch)
!  if(allocated(LeafStrutElms_brch))deallocate(LeafStrutElms_brch)
!  if(allocated(CanopyNodulStrutElms_brch))deallocate(CanopyNodulStrutElms_brch)
!  if(allocated(PetoleStrutElms_brch))deallocate(PetoleStrutElms_brch)
!  if(allocated(EarStrutElms_brch))deallocate(EarStrutElms_brch)
!  if(allocated(HuskStrutElms_brch))deallocate(HuskStrutElms_brch)
!  if(allocated(GrainStrutElms_brch))deallocate(GrainStrutElms_brch)
!  if(allocated(StalkStrutElms_brch))deallocate(StalkStrutElms_brch)
!  if(allocated(ShootStrutElms_brch))deallocate(ShootStrutElms_brch)
!  if(allocated(SenecStalkStrutElms_brch))deallocate(SenecStalkStrutElms_brch)
!  if(allocated(WTSTDI))deallocate(WTSTDI)
!  if(allocated(SeedCPlanted_pft))deallocate(SeedCPlanted_pft)
!  if(allocated(WTLS))deallocate(WTLS)
!  if(allocated(ShootStrutElms_pft))deallocate(ShootStrutElms_pft)
!  if(allocated(AvgCanopyBiomC2Graze_pft))deallocate(AvgCanopyBiomC2Graze_pft)
!  if(allocated(StandDeadStrutElms_pft)deallocate(StandDeadStrutElms_pft)
!  if(allocated(NodulStrutElms_pft))deallocate(NodulStrutElms_pft)
  end subroutine plt_biom_destroy


!----------------------------------------------------------------------

  subroutine  plt_soilchem_init(this)

  implicit none

  class(plant_soilchem_type) :: this

  allocate(this%FracBulkSOMC_vr(1:jcplx,0:JZ1));this%FracBulkSOMC_vr=spval
  allocate(this%ElmAllocmat4Litr(NumPlantChemElms,0:NumLitterGroups,jsken,JP1));this%ElmAllocmat4Litr=spval
  allocate(this%TScal4Difsvity_vr(0:JZ1));this%TScal4Difsvity_vr=spval
  allocate(this%THETPM(60,0:JZ1));this%THETPM=spval
  allocate(this%DiffusivitySolutEff(60,0:JZ1));this%DiffusivitySolutEff=spval
  allocate(this%VLSoilMicP_vr(0:JZ1));this%VLSoilMicP_vr=spval
  allocate(this%VLiceMicP_vr(0:JZ1));this%VLiceMicP_vr=spval
  allocate(this%VLWatMicP_vr(0:JZ1));this%VLWatMicP_vr=spval
  allocate(this%VLMicP_vr(0:JZ1));this%VLMicP_vr=spval
  allocate(this%trcs_VLN_vr(ids_nuts_beg:ids_nuts_end,0:JZ1));this%trcs_VLN_vr=spval
  allocate(this%DOM_vr(idom_beg:idom_end,1:jcplx,0:JZ1));this%DOM_vr=spval
  allocate(this%trc_solml_vr(ids_beg:ids_end,0:JZ1));this%trc_solml_vr=spval

  allocate(this%trc_gasml_vr(idg_beg:idg_end,0:JZ1));this%trc_gasml_vr=spval
  allocate(this%trc_gascl_vr(idg_beg:idg_end,0:JZ1));this%trc_gascl_vr=spval
  allocate(this%CSoilOrgM_vr(1:NumPlantChemElms,0:JZ1));this%CSoilOrgM_vr=spval

  allocate(this%trc_solcl_vr(ids_beg:ids_end,0:jZ1));this%trc_solcl_vr=spval

  allocate(this%VLSoilPoreMicP_vr(0:JZ1));this%VLSoilPoreMicP_vr=spval
  allocate(this%THETW_vr(0:JZ1));this%THETW_vr=spval
  allocate(this%THETY_vr(0:JZ1));this%THETY_vr=spval

  allocate(this%GasSolbility_vr(idg_beg:idg_end,0:JZ1));this%GasSolbility_vr=spval
  allocate(this%GasDifc_vr(idg_beg:idg_end,0:JZ1));this%GasDifc_vr=spval
  allocate(this%SoluteDifusvty_vr(ids_beg:ids_end,0:JZ1));this%SoluteDifusvty_vr=spval
  allocate(this%SoilResit4RootPentrate_vr(JZ1));this%SoilResit4RootPentrate_vr=spval
  allocate(this%SoiBulkDensity_vr(0:JZ1));this%SoiBulkDensity_vr=spval
  allocate(this%HydroCondMicP4RootUptake_vr(JZ1));this%HydroCondMicP4RootUptake_vr=spval

  end subroutine plt_soilchem_init
!----------------------------------------------------------------------

  subroutine plt_soilchem_destroy(this)
  implicit none
  class(plant_soilchem_type) :: this

!  if(allocated(FracBulkSOMC_vr))deallocate(FracBulkSOMC_vr)

!  if(allocated(ElmAllocmat4Litr))deallocate(ElmAllocmat4Litr)
!  if(allocated(TScal4Difsvity_vr))deallocate(TScal4Difsvity_vr)
!  if(allocated(THETPM))deallocate(THETPM)
!  if(allocated(DiffusivitySolutEff))deallocate(DiffusivitySolutEff)
!  if(allocated(ZVSGL))deallocate(ZVSGL)
!  if(allocated(O2GSolubility))deallocate(O2GSolubility)

!  if(allocated(HydroCondMicP4RootUptake_vr))deallocate(HydroCondMicP4RootUptake_vr)
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
!  if(allocated(SoiBulkDensity_vr))deallocate(SoiBulkDensity_vr)
!  if(allocated(CO2G))deallocate(CO2G)
!  if(allocated(CLSGL))deallocate(CLSGL)
!  if(allocated(CQSGL))deallocate(CQSGL)

  end subroutine plt_soilchem_destroy
!----------------------------------------------------------------------


  subroutine InitPlantAPIData()

  implicit none

  JZ1                        => pltpar%JZ1
  NumOfCanopyLayers1         => pltpar%NumOfCanopyLayers1
  JP1                        => pltpar%JP1
  NumOfLeafAzimuthSectors1   => pltpar%NumOfLeafAzimuthSectors
  NumOfSkyAzimuSects1        => pltpar%NumOfSkyAzimuSects1
  NumOfLeafZenithSectors1    => pltpar%NumOfLeafZenithSectors1
  MaxNodesPerBranch1         => pltpar%MaxNodesPerBranch1
  !the following variable should be consistent with the soil bgc model
  jcplx => pltpar%jcplx
  jsken  => pltpar%jsken
  NumLitterGroups=> pltpar%NumLitterGroups
  MaxNumBranches    => pltpar%MaxNumBranches
  MaxNumRootAxes    => pltpar%MaxNumRootAxes
  NumOfPlantMorphUnits   => pltpar%NumOfPlantMorphUnits
  NumOfPlantLitrCmplxs => pltpar%NumOfPlantLitrCmplxs
  NumGrowthStages => pltpar%NumGrowthStages
  jroots => pltpar%jroots

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

  allocate(this%RadPAR_zsec(NumOfLeafZenithSectors1,NumOfSkyAzimuSects1,NumOfCanopyLayers1,JP1))
  allocate(this%RadDifPAR_zsec(NumOfLeafZenithSectors1,NumOfSkyAzimuSects1,NumOfCanopyLayers1,JP1))
  allocate(this%RadSWLeafAlbedo_pft(JP1))
  allocate(this%CanopyPARalbedo_pft(JP1))
  allocate(this%TAU_DirRadTransm(NumOfCanopyLayers1+1))
  allocate(this%TAU_RadThru(NumOfCanopyLayers1+1))
  allocate(this%LWRadCanopy_pft(JP1))
  allocate(this%RadSWbyCanopy_pft(JP1))
  allocate(this%OMEGX(NumOfSkyAzimuSects1,NumOfLeafZenithSectors1,NumOfLeafAzimuthSectors1))
  allocate(this%OMEGAG(NumOfSkyAzimuSects1))
  allocate(this%OMEGA(NumOfSkyAzimuSects1,NumOfLeafZenithSectors1,NumOfLeafAzimuthSectors1))
  allocate(this%SineLeafAngle(NumOfLeafZenithSectors1))
  allocate(this%CosineLeafAngle(NumOfLeafZenithSectors1))
  allocate(this%iScatteringDiffus(NumOfSkyAzimuSects1,NumOfLeafZenithSectors1,NumOfLeafAzimuthSectors1))
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
  allocate(this%MinCanPStomaResistH2O_pft(JP1));this%MinCanPStomaResistH2O_pft=spval
  allocate(this%LeafO2Solubility_pft(JP1));this%LeafO2Solubility_pft=spval
  allocate(this%CanopyBndlResist_pft(JP1));this%CanopyBndlResist_pft=spval
  allocate(this%CanPStomaResistH2O_pft(JP1));this%CanPStomaResistH2O_pft=spval
  allocate(this%DiffCO2Atmos2Intracel_pft(JP1));this%DiffCO2Atmos2Intracel_pft=spval
  allocate(this%AirConc_pft(JP1));this%AirConc_pft=spval
  allocate(this%CO2CuticleResist_pft(JP1));this%CO2CuticleResist_pft=spval
  allocate(this%LeafAUnshaded_zsec(NumOfLeafZenithSectors1,NumOfCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafAUnshaded_zsec=spval
  allocate(this%CPOOL3_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CPOOL3_node=spval
  allocate(this%CPOOL4_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CPOOL4_node=spval
  allocate(this%CMassCO2BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CMassCO2BundleSheath_node=spval
  allocate(this%CO2CompenPoint_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2CompenPoint_node=spval
  allocate(this%RubiscoCarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%RubiscoCarboxyEff_node=spval
  allocate(this%C4CarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%C4CarboxyEff_node=spval
  allocate(this%LigthSatCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%LigthSatCarboxyRate_node=spval
  allocate(this%LigthSatC4CarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%LigthSatC4CarboxyRate_node=spval
  allocate(this%NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%NutrientCtrlonC4Carboxy_node=spval
  allocate(this%CMassHCO3BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CMassHCO3BundleSheath_node=spval

  allocate(this%Vmax4RubiscoCarboxy_pft(MaxNodesPerBranch1,MaxNumBranches,JP1));this%Vmax4RubiscoCarboxy_pft=spval
  allocate(this%CO2lmtRubiscoCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2lmtRubiscoCarboxyRate_node=spval
  allocate(this%Vmax4PEPCarboxy_pft(MaxNodesPerBranch1,MaxNumBranches,JP1));this%Vmax4PEPCarboxy_pft=spval
  allocate(this%CO2lmtPEPCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2lmtPEPCarboxyRate_node=spval
  allocate(this%iPlantPhotosynthesisType(JP1));this%iPlantPhotosynthesisType=0
  allocate(this%Km4PEPCarboxy_pft(JP1));this%Km4PEPCarboxy_pft=spval
  allocate(this%O2L(JP1));this%O2L=spval
  allocate(this%CO2Solubility_pft(JP1));this%CO2Solubility_pft=spval
  allocate(this%CanopyGasCO2_pft(JP1));this%CanopyGasCO2_pft=spval
  allocate(this%ChillHours_pft(JP1));this%ChillHours_pft=spval
  allocate(this%SpecChloryfilAct_pft(JP1));this%SpecChloryfilAct_pft=spval
  allocate(this%LeafC3ChlorofilConc_pft(JP1));this%LeafC3ChlorofilConc_pft=spval
  allocate(this%FracLeafProtinAsPEPCarboxyl_pft(JP1));this%FracLeafProtinAsPEPCarboxyl_pft=spval
  allocate(this%LeafC4ChlorofilConc_pft(JP1));this%LeafC4ChlorofilConc_pft=spval
  allocate(this%LeafRuBPConc_pft(JP1));this%LeafRuBPConc_pft=spval
  allocate(this%VmaxPEPCarboxyRef_pft(JP1));this%VmaxPEPCarboxyRef_pft=spval
  allocate(this%VmaxRubOxyRef_pft(JP1));this%VmaxRubOxyRef_pft=spval
  allocate(this%VmaxRubCarboxyRef_pft(JP1));this%VmaxRubCarboxyRef_pft=spval
  allocate(this%XKCO2(JP1));this%XKCO2=spval
  allocate(this%XKO2(JP1));this%XKO2=spval
  allocate(this%RubiscoActivity_brch(MaxNumBranches,JP1));this%RubiscoActivity_brch=spval
  allocate(this%C4PhotosynDowreg_brch(MaxNumBranches,JP1));this%C4PhotosynDowreg_brch=spval
  allocate(this%aquCO2Intraleaf_pft(JP1));this%aquCO2Intraleaf_pft=spval
  allocate(this%Km4LeafaqCO2_pft(JP1));this%Km4LeafaqCO2_pft=spval
  allocate(this%Km4RubiscoCarboxy_pft(JP1));this%Km4RubiscoCarboxy_pft=spval
  allocate(this%LeafIntracellularCO2_pft(JP1));this%LeafIntracellularCO2_pft=spval
  allocate(this%O2I(JP1));this%O2I=spval
  allocate(this%RCS(JP1));this%RCS=spval
  allocate(this%CanPCi2CaRatio(JP1));this%CanPCi2CaRatio=spval
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
  allocate(this%MinNonstC2InitBranch_pft(JP1));this%MinNonstC2InitBranch_pft=spval
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
  allocate(this%ShutRutNonstElmntConducts_pft(JP1));this%ShutRutNonstElmntConducts_pft=spval
  allocate(this%SeedTempSens_pft(JP1));this%SeedTempSens_pft=spval
  allocate(this%HighTempLimitSeed_pft(JP1));this%HighTempLimitSeed_pft=spval
  allocate(this%PlantInitThermoAdaptZone(JP1));this%PlantInitThermoAdaptZone=spval
  allocate(this%iPlantThermoAdaptZone_pft(JP1));this%iPlantThermoAdaptZone_pft=0
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
  allocate(this%iPlantRootProfile_pft(JP1));this%iPlantRootProfile_pft=0
  allocate(this%KHiestGroLeafNode_brch(MaxNumBranches,JP1));this%KHiestGroLeafNode_brch=0
  allocate(this%KLowestGroLeafNode_brch(MaxNumBranches,JP1));this%KLowestGroLeafNode_brch=0
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
  allocate(this%HourReq4LeafOut_brch(NumOfCanopyLayers1,JP1));this%HourReq4LeafOut_brch=spval
  allocate(this%Hours4LeafOff_brch(MaxNumBranches,JP1));this%Hours4LeafOff_brch=spval
  allocate(this%HourReq4LeafOff_brch(NumOfCanopyLayers1,JP1));this%HourReq4LeafOff_brch=spval

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
  allocate(this%RootPoreVol_pvr(jroots,JZ1,JP1));this%RootPoreVol_pvr=spval
  allocate(this%RootVH2O_pvr(jroots,JZ1,JP1));this%RootVH2O_pvr=spval
  allocate(this%Root1stXNumL_pvr(jroots,JZ1,JP1));this%Root1stXNumL_pvr=spval
  allocate(this%Root2ndXNum_pvr(jroots,JZ1,JP1));this%Root2ndXNum_pvr=spval
  allocate(this%SeedCMass_pft(JP1));this%SeedCMass_pft=spval
  allocate(this%totRootLenDens_vr(JZ1));this%totRootLenDens_vr=spval
  allocate(this%RootBranchFreq_pft(JP1));this%RootBranchFreq_pft=spval
  allocate(this%ClumpFactorInit_pft(JP1));this%ClumpFactorInit_pft=spval
  allocate(this%ClumpFactorNow_pft(JP1));this%ClumpFactorNow_pft=spval
  allocate(this%HypoctoHeight_pft(JP1));this%HypoctoHeight_pft=spval
  allocate(this%RootPoreTortu4Gas(jroots,JP1));this%RootPoreTortu4Gas=spval
  allocate(this%rLen2WidthLeaf_pft(JP1));this%rLen2WidthLeaf_pft=spval
  allocate(this%MaxSeedNumPerSite_pft(JP1));this%MaxSeedNumPerSite_pft=spval
  allocate(this%MaxPotentSeedNumber_pft(JP1));this%MaxPotentSeedNumber_pft=spval
  allocate(this%MaxSeedCMass_pft(JP1));this%MaxSeedCMass_pft=spval
  allocate(this%Root1stMaxRadius1_pft(jroots,JP1));this%Root1stMaxRadius1_pft=spval
  allocate(this%Root2ndMaxRadius1_pft(jroots,JP1));this%Root2ndMaxRadius1_pft=spval
  allocate(this%RootRaidus_rpft(jroots,JP1));this%RootRaidus_rpft=spval
  allocate(this%Root1stRadius_pvr(jroots,JZ1,JP1));this%Root1stRadius_pvr=spval
  allocate(this%Root2ndRadius_pvr(jroots,JZ1,JP1));this%Root2ndRadius_pvr=spval
  allocate(this%Root1stMaxRadius_pft(jroots,JP1));this%Root1stMaxRadius_pft=spval
  allocate(this%Root2ndMaxRadius_pft(jroots,JP1));this%Root2ndMaxRadius_pft=spval

  allocate(this%Root1stDepz_pft(jroots,NumOfCanopyLayers1,JP1));this%Root1stDepz_pft=spval
  allocate(this%RootLenPerPlant_pvr(jroots,JZ1,JP1));this%RootLenPerPlant_pvr=spval
  allocate(this%Root2ndAveLen_pvr(jroots,JZ1,JP1));this%Root2ndAveLen_pvr=spval
  allocate(this%Root1stSpecLen_pft(jroots,JP1));this%Root1stSpecLen_pft=spval
  allocate(this%Root2ndSpecLen_pft(jroots,JP1));this%Root2ndSpecLen_pft=spval
  allocate(this%Root1stLen_rpvr(jroots,JZ1,NumOfCanopyLayers1,JP1));this%Root1stLen_rpvr=spval
  allocate(this%Root2ndLen_pvr(jroots,JZ1,NumOfCanopyLayers1,JP1));this%Root2ndLen_pvr=spval
  allocate(this%Root2ndXNum_rpvr(jroots,JZ1,NumOfCanopyLayers1,JP1));this%Root2ndXNum_rpvr=spval
  allocate(this%iPlantNfixType_pft(JP1));this%iPlantNfixType_pft=0
  allocate(this%MY(JP1));this%MY=0
  allocate(this%CanPHeight4WatUptake(JP1));this%CanPHeight4WatUptake=spval
  allocate(this%KLeafNumber_brch(MaxNumBranches,JP1));this%KLeafNumber_brch=0
  allocate(this%NumOfLeaves_brch(MaxNumBranches,JP1));this%NumOfLeaves_brch=spval
  allocate(this%NGTopRootLayer_pft(JP1));this%NGTopRootLayer_pft=0;
  allocate(this%CanopyHeight_pft(JP1));this%CanopyHeight_pft=spval
  allocate(this%ShootNodeNumAtPlanting_pft(JP1));this%ShootNodeNumAtPlanting_pft=spval
  allocate(this%CanopyHeightZ_col(0:NumOfCanopyLayers1));this%CanopyHeightZ_col=spval
  allocate(this%CanopyStemArea_pft(JP1));this%CanopyStemArea_pft=spval
  allocate(this%CanopyLeafArea_pft(JP1));this%CanopyLeafArea_pft=spval
  allocate(this%MainBranchNum_pft(JP1));this%MainBranchNum_pft=0
  allocate(this%NIXBotRootLayer_pft(JP1));this%NIXBotRootLayer_pft=0
  allocate(this%NumRootAxes_pft(JP1));this%NumRootAxes_pft=0
  allocate(this%NumCogrothNode_pft(JP1));this%NumCogrothNode_pft=0
  allocate(this%BranchNumber_pft(JP1));this%BranchNumber_pft=0
  allocate(this%NumOfBranches_pft(JP1));this%NumOfBranches_pft=0
  allocate(this%NIXBotRootLayer_rpft(NumOfCanopyLayers1,JP1));this%NIXBotRootLayer_rpft=0
  allocate(this%PARTS_brch(NumOfPlantMorphUnits,MaxNumBranches,JP1));this%PARTS_brch=spval
  allocate(this%ShootNodeNum_brch(MaxNumBranches,JP1));this%ShootNodeNum_brch=spval
  allocate(this%NodeNum2InitFloral_brch(MaxNumBranches,JP1));this%NodeNum2InitFloral_brch=spval
  allocate(this%NodeNumberAtAnthesis_brch(MaxNumBranches,JP1));this%NodeNumberAtAnthesis_brch=spval
  allocate(this%SineBranchAngle_pft(JP1));this%SineBranchAngle_pft=spval
  allocate(this%LeafAreaZsec_brch(NumOfLeafZenithSectors1,NumOfCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafAreaZsec_brch=spval
  allocate(this%KMinNumLeaf4GroAlloc_brch(MaxNumBranches,JP1));this%KMinNumLeaf4GroAlloc_brch=0
  allocate(this%BranchNumber_brch(MaxNumBranches,JP1));this%BranchNumber_brch=0
  allocate(this%PotentialSeedSites_brch(MaxNumBranches,JP1));this%PotentialSeedSites_brch=spval
  allocate(this%CanPBranchHeight(MaxNumBranches,JP1));this%CanPBranchHeight=spval
  allocate(this%LeafAreaDying_brch(MaxNumBranches,JP1));this%LeafAreaDying_brch=spval
  allocate(this%LeafAreaLive_brch(MaxNumBranches,JP1));this%LeafAreaLive_brch=spval
  allocate(this%SinePetioleAngle_pft(JP1));this%SinePetioleAngle_pft=spval
  allocate(this%CLASS(NumOfLeafZenithSectors1,JP1));this%CLASS=spval
  allocate(this%CanopyStemAreaZ_pft(NumOfCanopyLayers1,JP1));this%CanopyStemAreaZ_pft=spval
  allocate(this%CanopyLeafAreaZ_pft(NumOfCanopyLayers1,JP1));this%CanopyLeafAreaZ_pft=spval
  allocate(this%LeafAreaNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafAreaNode_brch=spval
  allocate(this%CanopySeedNum_pft(JP1));this%CanopySeedNum_pft=spval
  allocate(this%SeedDepth_pft(JP1));this%SeedDepth_pft=spval
  allocate(this%PlantinDepz_pft(JP1));this%PlantinDepz_pft=spval
  allocate(this%SeedMeanLen_pft(JP1));this%SeedMeanLen_pft=spval
  allocate(this%SeedVolumeMean_pft(JP1));this%SeedVolumeMean_pft=spval
  allocate(this%SeedAreaMean_pft(JP1));this%SeedAreaMean_pft=spval
  allocate(this%CanopyStemAareZ_col(NumOfCanopyLayers1));this%CanopyStemAareZ_col=spval
  allocate(this%PetoLen2Mass_pft(JP1));this%PetoLen2Mass_pft=spval
  allocate(this%NodeLenPergC(JP1));this%NodeLenPergC=spval
  allocate(this%SLA1_pft(JP1));this%SLA1_pft=spval
  allocate(this%CanopyLeafAareZ_col(NumOfCanopyLayers1));this%CanopyLeafAareZ_col=spval
  allocate(this%LeafStalkArea_pft(JP1));this%LeafStalkArea_pft=spval
  allocate(this%InternodeHeightDying_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%InternodeHeightDying_brch=spval
  allocate(this%PetoleLensNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%PetoleLensNode_brch=spval
  allocate(this%LiveInterNodeHight_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%LiveInterNodeHight_brch=spval
  allocate(this%StemAreaZsec_brch(NumOfLeafZenithSectors1,NumOfCanopyLayers1,MaxNumBranches,JP1));this%StemAreaZsec_brch=0._r8
  allocate(this%CanopyLeafArea_lpft(NumOfCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%CanopyLeafArea_lpft=0._r8
  allocate(this%CanopyStalkArea_lbrch(NumOfCanopyLayers1,MaxNumBranches,JP1));this%CanopyStalkArea_lbrch=spval
  allocate(this%MaxSoiL4Root_pft(JP1));this%MaxSoiL4Root_pft=0
  allocate(this%SeedNumSet_brch(MaxNumBranches,JP1));this%SeedNumSet_brch=spval
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
