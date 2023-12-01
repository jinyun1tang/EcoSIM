module PlantAPIData
  use data_kind_mod, only : r8 => DAT_KIND_R8
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
  integer, pointer :: JRS          !maximum number of root layers
  integer, pointer :: MaxNumBranches         !maximum number of plant branches
  integer, pointer :: JP1         !number of plants
  integer, pointer :: JSA1        !number of sectors for the sky azimuth  [0,2*pi]
  integer, pointer :: jcplx       !number of organo-microbial CO2CompenPoint_nodeexes
  integer, pointer :: JLA1        !number of sectors for the leaf azimuth, [0,pi]
  integer, pointer :: NumOfCanopyLayers1         !number of canopy layers
  integer, pointer :: JZ1         !number of soil layers
  integer, pointer :: JLI1        !number of sectors for the leaf zenith [0,pi/2]
  integer, pointer :: MaxNodesPerBranch1      !number of canopy nodes
  integer, pointer :: jsken       !number of kinetic components in litter,  PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  integer, pointer :: NumLitterGroups     !number of litter groups nonstructural(0,*),
                          !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
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
  real(r8) :: ZNOON     !time of solar noon, [h]
  real(r8) :: CO2E      !atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: CCO2E     !atmospheric CO2 concentration, [g m-3]
  real(r8) :: CCH4E     !atmospheric CH4 concentration, [g m-3]
  real(r8) :: CZ2OE     !atmospheric N2O concentration, [g m-3]
  real(r8) :: CH2GE     !atmospheric H2 concentration, [g m-3]
  real(r8) :: CNH3E     !atmospheric NH3 concentration, [g m-3]
  real(r8) :: DayLenthPrev      !daylength of previous day, [h]
  real(r8) :: DayLenthCurrent      !daylength, [h]
  real(r8) :: DayLenthMax      !maximum daylength, [h]
  real(r8) :: SoiSurfRoughnesst0        !initial soil surface roughness height, [m]
  real(r8) :: OXYE      !atmospheric O2 concentration, [umol mol-1]
  real(r8) :: PPT       !total plant population, [d-2]
  real(r8) :: POROS1    !top layer soil porosity
  real(r8) :: WindSpeedAtm        !wind speed, [m h-1]
  real(r8), pointer :: PlantElemntStoreLandscape(:)  => null() !total plant element balance	g d-2
  real(r8) :: WindMesHeight        !wind speed measurement height, [m]
  real(r8) :: ZEROS2    !threshold zero
  real(r8) :: ZEROS     !threshold zero
  real(r8) :: ZERO      !threshold zero
  real(r8) :: ZERO2     !threshold zero
  integer :: IETYP      !Koppen climate zone
  integer :: NumActivePlants      !number of active PFT
  integer :: NL         !lowest soil layer number
  integer :: NP0        !intitial number of plant species
  integer :: MaxRootLayNum         !maximum root layer number
  integer :: NP         !number of plant species
  integer :: NU         !soil surface layer number
  integer :: LYRC       !number of days in current year
  integer :: iYearCurrent       !current year

  real(r8) :: VOLWOU    !total subsurface water flux	m3 d-2
  character(len=16), pointer :: DATAP(:) => null()   !parameter file name
  CHARACTER(len=16), pointer :: DATA(:)  => null()   !pft file
  real(r8), pointer :: AREA3(:)   => null()    !soil cross section area (vertical plan defined by its normal direction)
  integer,  pointer :: IDATA(:)   => null()    !time keeper
  real(r8), pointer :: ElmntBalanceCum_pft(:,:)  => null()    !plant element balance, [g d-2]
  real(r8), pointer :: CumSoilThickness(:)  => null()    !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8), pointer :: PPI(:)     => null()    !initial plant population, [m-2]
  real(r8), pointer :: PPZ(:)     => null()    !plant population at seeding, [m-2]
  real(r8), pointer :: PPX(:)     => null()    !plant population, [m-2]
  real(r8), pointer :: pftPlantPopulation(:)      => null()    !plant population, [d-2]
  real(r8), pointer :: DPTHZ(:)   => null()    !depth to middle of soil layer from  surface of grid cell [m]
  real(r8), pointer :: FracSoiAsMicP(:)    => null()    !micropore fraction
  real(r8), pointer :: DLYR3(:)   => null()    !vertical thickness of soil layer [m]
  real(r8), pointer :: VLWatMicPM(:,:) => null()    !soil micropore water content, [m3 d-2]
  real(r8), pointer :: VLsoiAirPM(:,:) => null()    !soil air content, [m3 d-2]
  real(r8), pointer :: TortMicPM(:,:)  => null()    !micropore soil tortuosity, []
  real(r8), pointer :: FILM(:,:)  => null()    !soil water film thickness , [m]
  contains
    procedure, public :: Init =>  plt_site_Init
    procedure, public :: Destroy => plt_site_destroy
  end type plant_siteinfo_type

  type, public :: plant_photosyns_type
  real(r8), pointer :: ETMX(:)   => null()  !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: CHL(:)    => null()  !leaf C3 chlorophyll content, [gC gC-1]
  real(r8), pointer :: PEPC(:)   => null()  !leaf PEP carboxylase content, [gC gC-1]
  real(r8), pointer :: CHL4(:)   => null()  !leaf C4 chlorophyll content, [gC gC-1]
  real(r8), pointer :: RUBP(:)   => null()  !leaf rubisco content, [gC gC-1]
  real(r8), pointer :: VCMX4(:)  => null()  !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: VOMX(:)   => null()  !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: VCMX(:)   => null()  !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: XKCO2(:)  => null()  !Km for rubisco carboxylase activity, [uM]
  real(r8), pointer :: XKO2(:)   => null()  !Km for rubisco oxygenase activity, [uM]
  real(r8), pointer :: RubiscoActivity_brpft(:,:) => null()   !branch down-regulation of CO2 fixation, [-]
  real(r8), pointer :: C4PhotosynDowreg_brch(:,:)=> null()   !down-regulation of C4 photosynthesis, [-]
  real(r8), pointer :: CO2L(:)   => null()   !leaf aqueous CO2 concentration, [uM]
  real(r8), pointer :: O2I(:)    => null()   !leaf gaseous O2 concentration, [umol m-3]
  real(r8), pointer :: LeafIntracellularCO2_pft(:)   => null()   !leaf gaseous CO2 concentration, [umol m-3]
  real(r8), pointer :: Km4RubiscoCarboxy_pft(:) => null()   !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8), pointer :: XKCO2L(:) => null()   !leaf aqueous CO2 Km no O2, [uM]
  real(r8), pointer :: RCS(:)    => null()   !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8), pointer :: CanPCi2CaRatio(:)   => null()   !Ci:Ca ratio, [-]
  real(r8), pointer :: MaxCanPStomaResistH2O(:)   => null()   !maximum stomatal resistance to vapor, [s h-1]
  real(r8), pointer :: CHILL(:)  => null()   !chilling effect on CO2 fixation, [-]
  real(r8), pointer :: CO2Solubility_pft(:)   => null()   !leaf CO2 solubility, [uM /umol mol-1]
  real(r8), pointer :: CanopyGasCO2_pft(:)   => null()   !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8), pointer :: O2L(:)    => null()   !leaf aqueous O2 concentration, [uM]
  real(r8), pointer :: Km4PEPCarboxy_pft(:) => null()   !Km for PEP carboxylase activity, [uM]
  integer,  pointer :: iPlantPhotosynthesisType(:)  => null()   !plant photosynthetic type (C3 or C4)
  real(r8), pointer :: MinCanPStomaResistH2O(:)   => null()   !canopy minimum stomatal resistance, [s m-1]
  real(r8), pointer :: RSMX(:)   => null()   !maximum stomatal resistance to vapor, [s m-1]
  real(r8), pointer :: CanPStomaResistH2O(:)     => null()   !canopy stomatal resistance, [h m-1]
  real(r8), pointer :: CanPbndlResist(:)     => null()   !canopy boundary layer resistance, [h m-1]
  real(r8), pointer :: SO2(:)    => null()   !leaf O2 solubility, [uM /umol mol-1]
  real(r8), pointer :: CO2CuticleResist_pft(:)   => null()   !maximum stomatal resistance to CO2, [s h-1]
  real(r8), pointer :: DiffCO2Atmos2Intracel_pft(:)   => null()   !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8), pointer :: LeafAUnshaded_seclyrnodbrpft(:,:,:,:,:)  => null()   !leaf irradiated surface area, [m2 d-2]
  real(r8), pointer :: CPOOL3(:,:,:)     => null()   !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8), pointer :: CPOOL4(:,:,:)     => null()   !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
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
  real(r8) :: TYSIN     !sine of sky angles
  real(r8) :: SoilAlbedo      !soil albedo
  real(r8) :: ALBX      !Surface albedo
  real(r8) :: RAP0      !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8) :: RADY      !diffuse shortwave radiation, [W m-2]
  real(r8) :: RADS      !direct shortwave radiation, [W m-2]
  real(r8) :: RAPY      !diffuse PAR, [umol m-2 s-1]
  real(r8) :: RAPS      !direct PAR, [umol m-2 s-1]
  real(r8) :: Eco_NetRad_col       !ecosystem net radiation, [MJ d-2 h-1]
  real(r8) :: RAD0      !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8) :: FracSWRad2Grnd     !fraction of radiation intercepted by ground surface, [-]
  real(r8) :: SWRadOnGrnd      !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8) :: GSIN      !sine of slope, [-]
  real(r8) :: GAZI      !azimuth of slope, [-]
  real(r8) :: GCOS      !cosine of slope, [-]
  real(r8) :: LWRadGrnd    !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8) :: LWRadSky       !sky longwave radiation , [MJ d-2 h-1]
  real(r8) :: SineSolarAngleNextHour     !sine of solar angle next hour, [-]
  real(r8) :: SineSolarAngle      !sine of solar angle, [-]
  real(r8), pointer :: CanopySWAlbedo_pft(:)     => null() !canopy shortwave albedo , [-]
  real(r8), pointer :: CanopyPARalbedo_pft(:)     => null() !canopy PAR albedo , [-]
  real(r8), pointer :: TAUS(:)     => null() !fraction of radiation intercepted by canopy layer, [-]
  real(r8), pointer :: TAU0(:)     => null() !fraction of radiation transmitted by canopy layer, [-]
  real(r8), pointer :: LWRadCanP(:)    => null() !canopy longwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: SWRadByCanP(:)     => null() !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: OMEGX(:,:,:)=> null() !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation
  real(r8), pointer :: OMEGAG(:)   => null() !sine of solar beam on leaf surface, [-]
  real(r8), pointer :: OMEGA(:,:,:)=> null() !sine of indirect sky radiation on leaf surface
  real(r8), pointer :: ZSIN(:)     => null() !sine of leaf angle
  real(r8), pointer :: ZCOS(:)     => null() !cosine of leaf angle
  integer,  pointer :: IALBY(:,:,:)=> null() !flag for calculating backscattering of radiation in canopy
  real(r8), pointer :: RadNet2CanP(:)     => null() !canopy net radiation , [MJ d-2 h-1]
  real(r8), pointer :: CanopySWabsorpty_pft(:)     => null() !canopy shortwave absorptivity , [-]
  real(r8), pointer :: CanopyPARabsorpty_pft(:)     => null() !canopy PAR absorptivity
  real(r8), pointer :: TAUR(:)     => null() !canopy shortwave transmissivity , [-]
  real(r8), pointer :: TAUP(:)     => null() !canopy PAR transmissivity , [-]
  real(r8), pointer :: PARbyCanopy_pft(:)     => null() !canopy absorbed PAR , [umol m-2 s-1]
  real(r8), pointer :: FracPARbyCanopy_pft(:)    => null() !fraction of incoming PAR absorbed by canopy, [-]
  real(r8), pointer :: PAR(:,:,:,:)=> null()     !direct incoming PAR, [umol m-2 s-1]
  real(r8), pointer :: PARDIF(:,:,:,:)=> null()  !diffuse incoming PAR, [umol m-2 s-1]
  contains
    procedure, public :: Init    => plt_rad_init
    procedure, public :: Destroy => plt_rad_destroy
  end type plant_radiation_type

  type, public :: plant_morph_type
  real(r8) :: CanopyArea_grid                     !stalk area of combined, each PFT canopy
  real(r8) :: CanopyLA_grd                     !grid canopy leaf area, [m2 d-2]
  real(r8) :: StemAreag                 !grid canopy stem area, [m2 d-2]
  real(r8) :: GridMaxCanopyHeight                        !canopy height , [m]
  real(r8), pointer :: RootVolPerMassC_pft(:,:)       => null() !root volume:mass ratio, [m3 g-1]
  real(r8), pointer :: RootPorosity(:,:)       => null() !root porosity, [m3 m-3]
  real(r8), pointer :: SecndRootXSecArea(:,:)     => null() !root  cross-sectional area  secondary axes, [m2]
  real(r8), pointer :: PrimRootXSecArea(:,:)     => null() !root cross-sectional area primary axes, [m2]
  real(r8), pointer :: MaxPrimRootRadius1(:,:)     => null() !root diameter primary axes, [m]
  real(r8), pointer :: MaxSecndRootRadius1(:,:)     => null() !root diameter secondary axes, [m]
  real(r8), pointer :: SeedCMass(:)         => null() !grain size at seeding, [g]
  real(r8), pointer :: RootPoreTortu4Gas(:,:)      => null() !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8), pointer :: PrimRootLen(:,:,:,:)  => null() !root layer length primary axes, [m d-2]
  real(r8), pointer :: SecndRootLen(:,:,:,:)  => null() !root layer length secondary axes, [m d-2]
  real(r8), pointer :: RootLenPerPopu_pvr(:,:,:)    => null() !root layer length per plant, [m p-1]
  real(r8), pointer :: AveSecndRootLen(:,:,:)    => null() !root layer average length, [m]
  real(r8), pointer :: PrimRootSpecLen(:,:)     => null() !specific root length primary axes, [m g-1]
  real(r8), pointer :: SecndRootSpecLen(:,:)     => null() !specific root length secondary axes, [m g-1]
  real(r8), pointer :: SecndRootXNum_rpvr(:,:,:,:)   => null() !root layer number secondary axes, [d-2]
  real(r8), pointer :: PrimRootDepth(:,:,:)    => null() !root layer depth, [m]
  real(r8), pointer :: ClumpFactort0(:)          => null() !initial clumping factor for self-shading in canopy layer, [-]
  real(r8), pointer :: ClumpFactort(:)          => null() !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8), pointer :: RootBranchFreq_pft(:)         => null() !root brancing frequency, [m-1]
  real(r8), pointer :: HypoctoylHeight(:)        => null() !cotyledon height, [m]
  integer,  pointer :: iPlantNfixType(:)        => null() !N2 fixation type
  integer,  pointer :: MY(:)           => null() !mycorrhizal type (no or yes)
  real(r8), pointer :: CanPHeight4WatUptake(:)        => null() !canopy height, [m]
  real(r8), pointer :: CanopyHeightz(:)  => null() !canopy layer height , [m]
  real(r8), pointer :: CanopyStemA_pft(:)        => null() !plant stem area, [m2 d-2]
  real(r8), pointer :: CanopyLeafA_pft(:)        => null() !plant canopy leaf area, [m2 d-2]
  integer,  pointer :: NumOfMainBranch_pft(:)          => null() !number of main branch
  integer,  pointer :: NI(:)           => null() !maximum soil layer number for all root axes
  integer,  pointer :: NIXBotRootLayer_pft(:)          => null() !maximum soil layer number for all root axes, [-]
  integer,  pointer :: NumRootAxes_pft(:)          => null() !root primary axis number
  integer,  pointer :: NumConCurrentGrowinNode(:)         => null() !number of concurrently growing nodes
  integer,  pointer :: BranchNumber_pft(:)          => null() !branch number
  integer,  pointer :: NumOfBranches_pft(:)          => null() !branch number
  integer,  pointer :: NIXBotRootLayer_rpft(:,:)       => null() !maximum soil layer number for root axes, [-]
  real(r8), pointer :: ShootNodeNumber(:,:)       => null() !shoot node number, [-]
  real(r8), pointer :: NodeNumberToInitFloral(:,:)      => null() !shoot node number at floral initiation, [-]
  real(r8), pointer :: NodeNumberAtAnthesis(:,:)      => null() !shoot node number at anthesis, [-]
  real(r8), pointer :: SineBranchAngle_pft(:)        => null() !branching angle, [degree from horizontal]
  real(r8), pointer :: LeafA_lyrnodbrchpft(:,:,:,:,:) => null() !leaf surface area, [m2 d-2]
  integer,  pointer :: KLEAFX(:,:)     => null() !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,  pointer :: BranchNumber_brch(:,:)       => null() !branch number, [-]
  real(r8), pointer :: PotentialSeedSites_brch(:,:)      => null() !branch potential grain number, [d-2]
  real(r8), pointer :: SeedNumberSet_brch(:,:)      => null() !branch grain number, [d-2]
  real(r8), pointer :: CanPBranchHeight(:,:)     => null() !branch height, [m]
  real(r8), pointer :: LeafAreaDying_brch(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: LeafAreaLive_brch(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: SinePetioleAngle_pft(:)        => null() !sheath angle, [degree from horizontal]
  real(r8), pointer :: CLASS(:,:)      => null() !fractionction of leaves in different angle classes, [-]
  real(r8), pointer :: CanopyStemApft_lyr(:,:)      => null() !plant canopy layer stem area, [m2 d-2]
  real(r8), pointer :: CanopyLeafApft_lyr(:,:)      => null() !canopy layer leaf area, [m2 d-2]
  real(r8), pointer :: LeafAreaNode_brch(:,:,:)    => null() !leaf area, [m2 d-2]
  real(r8), pointer :: CanopySeedNumber_pft(:)         => null() !canopy grain number, [d-2]
  real(r8), pointer :: SeedinDepth(:)        => null() !seeding depth, [m]
  real(r8), pointer :: PlantinDepth(:)       => null() !planting depth, [m]
  real(r8), pointer :: SeedLengthMean_pft(:)         => null() !seed length, [m]
  real(r8), pointer :: SeedVolumeMean_pft(:)         => null() !seed volume, [m3 ]
  real(r8), pointer :: SeedAreaMean_pft(:)         => null() !seed surface area, [m2]
  real(r8), pointer :: CanopyStemA_lyr(:)        => null() !total stem area, [m2 d-2]
  real(r8), pointer :: SSL1(:)         => null() !petiole length:mass during growth, [m gC-1]
  real(r8), pointer :: SNL1(:)         => null() !internode length:mass during growth, [m gC-1]
  real(r8), pointer :: SLA1(:)         => null() !leaf area:mass during growth, [m2 gC-1]
  real(r8), pointer :: CanopyLAgrid_lyr(:)        => null() !total leaf area, [m2 d-2]
  real(r8), pointer :: CanopyArea_pft(:)        => null() !plant leaf+stem/stalk area, [m2 d-2]
  real(r8), pointer :: InternodeHeightDying_brch(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: PetioleLengthNode_brch(:,:,:)    => null() !sheath height, [m]
  real(r8), pointer :: InternodeHeightLive_brch(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: StemA_lyrnodbrchpft(:,:,:,:)  => null() !stem surface area, [m2 d-2]
  real(r8), pointer :: CanopyLeafAreaByLayer_pft(:,:,:,:)  => null() !layer/node/branch leaf area, [m2 d-2]
  real(r8), pointer :: CanopyBranchStemApft_lyr(:,:,:)    => null() !plant canopy layer branch stem area, [m2 d-2]
  real(r8), pointer :: ClumpFactor(:)           => null() !clumping factor for self-shading in canopy layer, [-]
  real(r8), pointer :: XTLI(:)         => null() !number of nodes in seed, [-]
  real(r8), pointer :: CanopyHeight(:)           => null() !canopy height, [m]
  integer,  pointer :: NGTopRootLayer_pft(:)           => null() !soil layer at planting depth, [-]
  integer,  pointer :: KLeafNumber_brch(:,:)      => null() !leaf number, [-]
  real(r8), pointer :: NumOfLeaves_brch(:,:)       => null() !leaf number, [-]
  integer , pointer :: iPlantGrainType(:)        => null() !grain type (below or above-ground)
  real(r8), pointer :: MaxPotentSeedNumber_pft(:)         => null() !maximum grain node number per branch, [-]
  real(r8), pointer :: MaxSeedNumPerSite_pft(:)         => null() !maximum grain number per node , [-]
  real(r8), pointer :: WDLF(:)         => null() !leaf length:width ratio, [-]
  real(r8), pointer :: MaxSeedCMass(:)         => null() !maximum grain size   , [g]
  real(r8), pointer :: PrimRootRadius(:,:,:)    => null() !root layer diameter primary axes, [m ]
  real(r8), pointer :: SecndRootRadius(:,:,:)    => null() !root layer diameter secondary axes, [m ]
  real(r8), pointer :: RRADP(:,:)      => null() !root internal radius, [m]
  real(r8), pointer :: MaxPrimRootRadius(:,:)     => null() !maximum radius of primary roots, [m]
  real(r8), pointer :: MaxSecndRootRadius(:,:)     => null() !maximum radius of secondary roots, [m]
  real(r8), pointer :: RSRR(:,:)       => null() !root radial resistivity, [MPa h m-2]
  real(r8), pointer :: RSRA(:,:)       => null() !root axial resistivity, [MPa h m-4]
  real(r8), pointer :: RTDNT(:)        => null() !total root length density, [m m-3]
  real(r8), pointer :: PrimRootXNumL(:,:,:)     => null() !root layer number primary axes, [d-2]
  real(r8), pointer :: SecndRootXNum_pvr(:,:,:)     => null() !root layer number axes, [d-2]
  real(r8), pointer :: RootLenthDensPerPopu_pvr(:,:,:)    => null() !root layer length density, [m m-3]
  real(r8), pointer :: RootVolume_vr(:,:,:)    => null() !root layer volume air, [m2 d-2]
  real(r8), pointer :: RootVH2O_vr(:,:,:)    => null() !root layer volume water, [m2 d-2]
  real(r8), pointer :: RootAreaPerPlant_vr(:,:,:)    => null() !root layer area per plant, [m p-1]
  contains
    procedure, public :: Init    => plt_morph_init
    procedure, public :: Destroy => plt_morph_destroy
  end type plant_morph_type

  type, public :: plant_pheno_type
  real(r8), pointer :: fTgrowRootP(:,:)   => null()     !root layer temperature growth functiom, [-]
  real(r8), pointer :: ShutRutNonstructElmntConducts(:)    => null()     !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(r8), pointer :: GrainFillRateat25C_pft(:)    => null()     !maximum rate of fill per grain, [g h-1]
  real(r8), pointer :: OFFST(:)    => null()     !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
  real(r8), pointer :: PlantO2Stress(:)     => null()     !plant O2 stress indicator, []
  real(r8), pointer :: MinNonstructalC4InitBranch(:)       => null()     !branch nonstructural C content required for new branch, [gC gC-1]
  real(r8), pointer :: MinNonstructuralC4InitRoot(:)       => null()     !threshold root nonstructural C content for initiating new root axis, [gC gC-1]
  real(r8), pointer :: LeafElmntRemobFlx_brch(:,:,:) => null()    !element translocated from leaf during senescence, [g d-2 h-1]
  real(r8), pointer :: PetioleChemElmntRemobFlx_brch(:,:,:) => null()    !element translocated from sheath during senescence, [g d-2 h-1]
  real(r8), pointer :: TCelciusChill4Leaf(:)     => null()     !threshold temperature for spring leafout/dehardening, [oC]
  real(r8), pointer :: TCG(:)     => null()     !canopy growth temperature, [oC]
  real(r8), pointer :: TCelcius4LeafOffHarden(:)     => null()     !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8), pointer :: TKG(:)     => null()     !canopy growth temperature, [K]
  real(r8), pointer :: fTgrowCanP(:)    => null()     !canopy temperature growth function, [-]
  real(r8), pointer :: HoursCanopyPSITooLow(:)    => null()     !canopy plant water stress indicator, number of hours PSICanP < PSILY, []
  real(r8), pointer :: TCelciusChill4Seed(:)     => null()     !temperature below which seed set is adversely affected, [oC]
  real(r8), pointer :: iPlantThermoAdaptZone(:)    => null()     !plant thermal adaptation zone, [-]
  real(r8), pointer :: iPlantInitThermoAdaptZone(:)   => null()     !initial plant thermal adaptation zone, [-]
  real(r8), pointer :: HTC(:)     => null()     !temperature above which seed set is adversely affected, [oC]
  real(r8), pointer :: SSTX(:)    => null()     !sensitivity to HTC (seeds oC-1 above HTC)
  integer,  pointer :: iPlantState(:)    => null()     !flag for species death
  integer,  pointer :: IsPlantActive(:)   => null()     !flag for living pft
  real(r8), pointer :: RSETE(:,:) => null()     !effect of canopy element status on seed set , []
  real(r8), pointer :: RefNodeInitRate(:)    => null()     !rate of node initiation, [h-1 at 25 oC]
  real(r8), pointer :: RefLeafAppearRate(:)    => null()     !rate of leaf initiation, [h-1 at 25 oC]
  real(r8), pointer :: CriticalPhotoPeriod_pft(:)     => null()     !critical daylength for phenological progress, [h]
  real(r8), pointer :: PhotoPeriodSens_pft(:)    => null()     !difference between current and critical daylengths used to calculate  phenological progress [h]
  integer,  pointer :: iPlantBranchState(:,:)  => null()  !flag to detect branch death , [-]
  integer,  pointer :: doRemobilization_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doPlantLeaveOff_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doPlantLeafOut_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doInitLeafOut_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: doSenescence_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: Prep4Literfall_brch(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: Hours4LiterfalAftMature_brch(:,:)  => null()  !branch phenology flag, [h]
  integer,  pointer :: KLeafNodeNumber(:,:)  => null()  !leaf growth stage counter, [-]
  integer,  pointer :: iPlantPhenologyType(:)    => null()  !climate signal for phenological progress: none, temperature, water stress
  integer,  pointer :: iPlantTurnoverPattern(:)    => null()  !phenologically-driven above-ground turnover: all, foliar only, none
  integer,  pointer :: iPlantShootState(:)    => null()  !flag to detect canopy death
  integer,  pointer :: iPlantPhenologyPattern(:)    => null()  !plant growth habit: annual or perennial
  integer,  pointer :: iPlantRootState(:)    => null()  !flag to detect root system death
  integer,  pointer :: iPlantDevelopPattern(:)    => null()  !plant growth habit (determinate or indeterminate)
  integer,  pointer :: iPlantPhotoperiodType(:)    => null()  !photoperiod type (neutral, long day, short day)
  integer,  pointer :: iPlantMorphologyType(:)    => null()  !plant growth type (vascular, non-vascular)
  integer,  pointer :: doInitPlant(:)    => null()  !PFT initialization flag:0=no,1=yes
  integer,  pointer :: KLeafNumLowestGrowing_pft(:,:) => null()  !leaf growth stage counter, [-]
  integer,  pointer :: iPlantCalendar(:,:,:) => null()  !plant growth stage, [-]
  real(r8), pointer :: TotalNodeNumNormByMatgrp_brch(:,:) => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: TotalReprodNodeNumNormByMatrgrp_brch(:,:) => null()  !normalized node number during reproductive growth stages , [-]
  real(r8), pointer :: LeafNumberAtFloralInit_brch(:,:)  => null()  !leaf number at floral initiation, [-]
  real(r8), pointer :: Hours4LenthenPhotoPeriod(:,:)   => null()  !initial heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: Hours4ShortenPhotoPeriod(:,:)   => null()  !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: Hours4Leafout(:,:)   => null()  !heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: HourThreshold4LeafOut_brch(:,:)   => null()  !hours above threshold temperature required for spring leafout/dehardening, [-]
  real(r8), pointer :: Hours4LeafOff(:,:)   => null()  !cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: HourThreshold4LeafOff(:,:)   => null()  !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8), pointer :: HourCounter4LeafOut_brch(:,:)   => null()  !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8), pointer :: HoursDoingRemob_brch(:,:)   => null()  !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  real(r8), pointer :: HourlyNodeNumNormByMatgrp_brch(:,:) => null()  !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8), pointer :: HourReprodNodeNumNormByMatrgrp_brch(:,:) => null()  !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8), pointer :: MatureGroup_brch(:,:)  => null()  !plant maturity group, [-]
  real(r8), pointer :: NodeNumNormByMatgrp_brch(:,:)  => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: ReprodNodeNumNormByMatrgrp_brch(:,:)  => null()  !normalized node number during reproductive growth stages, [-]
  real(r8), pointer :: HourWithNoGrainFill_brch(:,:)   => null()  !flag to detect physiological maturity from  grain fill , [-]
  real(r8), pointer :: MatureGroup_pft(:)   => null()  !acclimated plant maturity group, [-]
  contains
    procedure, public :: Init    =>  plt_pheno_init
    procedure, public :: Destroy =>  plt_pheno_destroy
  end type plant_pheno_type

  type, public :: plant_soilchem_type
  real(r8), pointer :: FOSRH(:,:)  => null()  !fraction of total organic C in CO2CompenPoint_nodeex, [-]
  real(r8), pointer :: CFOPE(:,:,:,:)=> null() !litter kinetic fraction, [-]
  real(r8), pointer :: TFND(:)     => null()  !temperature effect on diffusivity
  real(r8), pointer :: THETPM(:,:) => null()  !soil air-filled porosity, [m3 m-3]
  real(r8), pointer :: DiffusivitySolutEff(:,:)   => null()  !coefficient for dissolution - volatilization, []
  real(r8), pointer :: SoilResit4RootPentration(:)     => null()  !soil hydraulic resistance, [MPa h m-2]
  real(r8), pointer :: SoiBulkDensity(:)     => null()  !soil bulk density, [Mg m-3]
  real(r8), pointer :: trc_solcl(:,:) => null() !aqueous tracer concentration [g m-3]
  real(r8), pointer :: trc_gascl(:,:) => null() !gaseous tracer concentration [g m-3]

  real(r8), pointer :: CORGC(:)    => null()  !soil organic C content [gC kg soil-1]

  real(r8), pointer :: HydroCondMicP4RootUptake(:)     => null()  !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]

  real(r8), pointer :: GasDifc(:,:)=> null()  !gaseous diffusivity [m2 h-1]
  real(r8), pointer :: SolDifc(:,:)=> null()  !aqueous diffusivity [m2 h-1]

  real(r8), pointer :: trc_gasml(:,:)=> null()!gas layer mass [g d-2]

  real(r8), pointer :: GasSolbility(:,:)=> null() !gas solubility, [m3 m-3]

  real(r8), pointer :: THETW(:)    => null()  !volumetric water content [m3 m-3]
  real(r8), pointer :: THETY(:)    => null()  !air-dry water content, [m3 m-3]
  real(r8), pointer :: VLSoilPoreMicP(:)     => null()  !volume of soil layer	m3 d-2
  real(r8), pointer :: trcs_VLN(:,:)=> null()

  real(r8), pointer :: VLSoilMicP(:)     => null()  !total micropore volume in layer [m3 d-2]
  real(r8), pointer :: VLiceMicP(:)     => null()  !soil micropore ice content   [m3 d-2]
  real(r8), pointer :: VLWatMicP(:)     => null()  !soil micropore water content [m3 d-2]
  real(r8), pointer :: VLMicP(:)     => null()  !total volume in micropores [m3 d-2]

  real(r8), pointer :: DOM(:,:,:)    => null()  !dissolved organic C micropore	[gC d-2]
  real(r8), pointer :: trc_solml(:,:)=> null() !aqueous tracer [g d-2]
  contains
    procedure, public :: Init => plt_soilchem_init
    procedure, public :: Destroy => plt_soilchem_destroy
  end type plant_soilchem_type

  type, public :: plant_allometry_type
  real(r8), pointer :: CPRTS(:)    => null()  !root P:C ratio x root growth yield, [-]
  real(r8), pointer :: CNRTS(:)    => null()  !root N:C ratio x root growth yield, [-]
  real(r8), pointer :: NodulerNC_pft(:)     => null()  !nodule N:C ratio, [gN gC-1]
  real(r8), pointer :: NodulerPC_pft(:)     => null()  !nodule P:C ratio, [gP gC-1]
  real(r8), pointer :: RootrNC_pft(:)     => null()  !root N:C ratio, [gN gC-1]
  real(r8), pointer :: RootrPC_pft(:)     => null()  !root P:C ratio, [gP gC-1]
  real(r8), pointer :: rCNNonstructRemob_pft(:)     => null()  !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: rCPNonstructRemob_pft(:)     => null()  !C:P ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: NoduGrowthYield_pft(:)     => null()  !nodule growth yield, [g g-1]
  real(r8), pointer :: RootBiomGrowthYield(:)     => null()  !root growth yield, [g g-1]
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
  real(r8), pointer :: FWODBE(:,:) => null()  !woody element allocation
  real(r8), pointer :: FWODLE(:,:) => null()  !leaf element allocation
  real(r8), pointer :: FWODRE(:,:) => null()  !C woody fraction in root
  real(r8), pointer :: FWOODE(:,:) => null()  !woody element allocation
  real(r8), pointer :: FVRN(:)     => null()  !allocation parameter
  real(r8), pointer :: LeafBiomGrowthYield(:)     => null()  !leaf growth yield, [g g-1]
  real(r8), pointer :: CNGR(:)     => null()  !grain N:C ratio, [g g-1]
  real(r8), pointer :: CPLF(:)     => null()  !maximum leaf P:C ratio, [g g-1]
  real(r8), pointer :: CPSHE(:)    => null()  !sheath P:C ratio, [g g-1]
  real(r8), pointer :: CNSHE(:)    => null()  !sheath N:C ratio, [g g-1]
  real(r8), pointer :: rPCStalk_pft(:)    => null()  !stalk P:C ratio, [g g-1]
  real(r8), pointer :: CNLF(:)     => null()  !maximum leaf N:C ratio, [g g-1]
  real(r8), pointer :: GrainSeedBiomCMean_brch(:,:)  => null()  !maximum grain C during grain fill, [g d-2]
  real(r8), pointer :: FNOD(:)     => null()  !parameter for allocation of growth to nodes, [-]
  real(r8), pointer ::RootFracRemobilizableBiom(:)    => null()  !fraction of remobilizable nonstructural biomass in root, [-]

  contains
    procedure, public :: Init => plt_allom_init
    procedure, public :: Destroy => plt_allom_destroy
  end type plant_allometry_type

  type, public :: plant_biom_type
  real(r8), pointer :: StandingDeadChemElmnt_col(:)     => null()    !total standing dead element, [g d-2]
  real(r8), pointer :: ZEROL(:)       => null()    !threshold zero for leaf calculation
  real(r8), pointer :: ZEROP(:)       => null()    !threshold zero for p calculation
  real(r8), pointer :: RootNoduleNonstructElmnt_vr(:,:,:)  => null()    !root  layer nonstructural element, [g d-2]
  real(r8), pointer :: RootNodueChemElmnt_pvr(:,:,:)  => null()    !root layer nodule element, [g d-2]
  real(r8), pointer :: CanopyLeafCpft_lyr(:,:)     => null()    !canopy layer leaf C, [g d-2]
  real(r8), pointer :: Root2ndStructChemElmnt_pvr(:,:,:,:,:) => null()    !root layer element secondary axes, [g d-2]
  real(r8), pointer :: Root1stStructChemElmnt_pvr(:,:,:,:,:) => null()    !root layer element primary axes, [g d-2]
  real(r8), pointer :: Root1stChemElmnt(:,:,:,:)   => null()    !root C primary axes, [g d-2]
  real(r8), pointer :: StandingDeadKCompChemElmnts_pft(:,:,:)  => null()    !standing dead element fraction, [g d-2]
  real(r8), pointer :: CanopyNonstructElementConc_pft(:,:)    => null()    !canopy nonstructural element concentration, [g d-2]
  real(r8), pointer :: CanopyNonstructElements_pft(:,:)    => null()    !canopy nonstructural element concentration, [g d-2]
  real(r8), pointer :: NoduleNonstructElmnt_pft(:,:)    => null()    !canopy nodule nonstructural element, [g d-2]
  real(r8), pointer :: NoduleNonstructCconc_pft(:)      => null()    !nodule nonstructural C, [gC d-2]
  real(r8), pointer :: RootStructBiomC_vr(:,:,:)   => null()    !root layer structural C, [g d-2]
  real(r8), pointer ::  PopuPlantRootC_vr(:,:,:)   => null()    !root layer C, [g d-2]
  real(r8), pointer :: RootProteinC_pvr(:,:,:)   => null()    !root layer protein C, [g d-2]
  real(r8), pointer :: RootProteinConc_pftvr(:,:,:)  => null()    !root layer protein C concentration, [g g-1]
  real(r8), pointer ::  RootMycoNonstructElmnt_vr(:,:,:,:)=> null()    !root  layer nonstructural element, [g d-2]
  real(r8), pointer :: RootNonstructElementConcpft_vr(:,:,:,:)  => null()    !root  layer nonstructural C concentration, [g g-1]
  real(r8), pointer :: LeafPetioNonstructElmntConc_brch(:,:,:)    => null()    !branch nonstructural C concentration, [g d-2]
  real(r8), pointer :: InternodeChemElmnt_brch(:,:,:,:)  => null()    !internode C, [g d-2]
  real(r8), pointer :: LeafChemElmntNode_brch(:,:,:,:)    => null()    !leaf element, [g d-2]
  real(r8), pointer :: LeafProteinCNode_brch(:,:,:)    => null()    !layer leaf protein C, [g d-2]
  real(r8), pointer :: PetioleElmntNode_brch(:,:,:,:)   => null()  !sheath element , [g d-2]
  real(r8), pointer :: PetioleProteinCNode_brch(:,:,:)   => null()    !layer sheath protein C, [g d-2]
  real(r8), pointer :: LeafChemElmntByLayer_pft(:,:,:,:,:) => null()    !layer leaf element, [g d-2]
  real(r8), pointer :: WGLFT(:)       => null()  !total leaf mass, [gC d-2]
  real(r8), pointer :: StandingDeadInitC_pft(:)      => null()  !initial standing dead C, [g C m-2]
  real(r8), pointer :: RootChemElmnts_pft(:,:)     => null()  !plant root element, [gC d-2]
  real(r8), pointer :: RootStructChemElmnt_pft(:,:)    => null()  !plant root structural element, [gC d-2]
  real(r8), pointer :: SeedCPlanted_pft(:)       => null()  !plant stored nonstructural C at planting, [gC d-2]
  real(r8), pointer :: NonstructalChemElmnts_pft(:,:)     => null()  !plant stored nonstructural element, [gC d-2]
  real(r8), pointer :: CanopyLeafShethC_pft(:)        => null()  !canopy leaf + sheath C, [g d-2]
  real(r8), pointer :: ShootChemElmnts_pft(:,:)    => null()  !canopy shoot C, [g d-2]
  real(r8), pointer :: WTSHTA(:)      => null()  !landscape average canopy shoot C, [g d-2]
  real(r8), pointer :: StandingDeadChemElmnts_pft(:,:)    => null()  !standing dead element, [g d-2]
  real(r8), pointer :: NoduleChemElmnts_pft(:,:)     => null()  !root total nodule mass, element [g d-2]
  real(r8), pointer :: NonstructElmnt_brch(:,:,:)   => null()  !branch nonstructural element, [g d-2]
  real(r8), pointer :: NoduleNonstructElmnt_brch(:,:,:)  => null()  !branch nodule nonstructural element, [g d-2]
  real(r8), pointer :: LeafPetioleBiomassC_brch(:,:)     => null()  !plant branch leaf + sheath C, [g d-2]
  real(r8), pointer :: ReserveChemElmnts_brch(:,:,:) => null()  !branch reserve element, [g d-2]
  real(r8), pointer :: LeafChemElmnts_brch(:,:,:)  => null()   !branch leaf element, [g d-2]
  real(r8), pointer :: CanopyNoduleChemElmnt_brch(:,:,:)  => null()   !branch nodule element, [g d-2]
  real(r8), pointer :: PetioleChemElmnts_brch(:,:,:) => null()   !branch sheath element , [g d-2]
  real(r8), pointer :: EarChemElmnts_brch(:,:,:) => null()   !branch ear C, [g d-2]
  real(r8), pointer :: HuskChemElmnts_brch(:,:,:) => null()   !branch husk element, [g d-2]
  real(r8), pointer :: GrainChemElmnts_brch(:,:,:)  => null()   !branch grain element, [g d-2]
  real(r8), pointer :: StalkChemElmnts_brch(:,:,:) => null()   !branch stalk element, [g d-2]
  real(r8), pointer :: ShootChemElmnt_brch(:,:,:) => null()   !branch shoot C, [g d-2]
  real(r8), pointer :: LeafChemElmntRemob_brch(:,:,:)  => null()   !branch leaf structural element, [g d-2]
  real(r8), pointer :: PetioleChemElmntRemob_brch(:,:,:) => null()   !branch sheath structural element, [g d-2]
  real(r8), pointer :: BranchStalkChemElmnts_pft_pft(:,:,:) => null()   !branch stalk structural element, [g d-2]
  real(r8), pointer :: StalkBiomassC_brch(:,:)    => null()   !branch active stalk C, [g d-2]
  real(r8), pointer :: StalkChemElmnts_pft(:,:)    => null()   !canopy stalk element, [g d-2]
  real(r8), pointer :: CanopyStalkC_pft(:)       => null()   !canopy active stalk C, [g d-2
  real(r8), pointer :: LeafChemElmnts_pft(:,:)     => null()   !canopy leaf elements, [g d-2]
  real(r8), pointer :: PetioleChemElmnts_pft(:,:)    => null()   !canopy sheath element , [g d-2]
  real(r8), pointer :: ReserveChemElmnts_pft(:,:)    => null()   !canopy reserve element, [g d-2]
  real(r8), pointer :: HuskChemElmnts_pft(:,:)    => null()   !canopy husk element, [g d-2]
  real(r8), pointer :: RootBiomCPerPlant_pft(:)       => null()   !root C per plant, [g p-1]
  real(r8), pointer :: GrainChemElmnts_pft(:,:)     => null()   !canopy grain element, [g d-2]
  real(r8), pointer :: EarChemElmnts_pft(:,:)    => null()   !canopy ear element, [g d-2]
  contains
    procedure, public :: Init => plt_biom_init
    procedure, public :: Destroy => plt_biom_destroy
  end type plant_biom_type

  type, public :: plant_ew_type
  real(r8) :: SnowDepth          !snowpack depth, [m]
  real(r8) :: VcumWatSnow        !water volume in snowpack, [m3 d-2]
  real(r8) :: VcumDrySnoWE       !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8) :: Canopy_Heat_Sens_col               !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VLHeatCapSnowMin    !minimum snowpack heat capacity [MJ d-2 K-1]
  real(r8) :: VLHeatCapSurfSnow    !snowpack heat capacity, [MJ m-3 K-1]
  real(r8) :: TENGYC    !total canopy heat content, [MJ  d-2]
  real(r8) :: LWRadCanG     !grid total canopy LW emission, [MJ d-2 h-1]
  real(r8) :: TEVAPP    !total canopy evaporation + transpiration, [m3 d-2]
  real(r8) :: VapXAir2CanG    !grid canopy evaporation, [m3 d-2]
  real(r8) :: THFLXC    !total canopy heat flux, [MJ  d-2]
  real(r8) :: UVOLO     !total subsurface water flux, [m3 d-2]
  real(r8) :: VPA       !vapor concentration, [m3 m-3]
  real(r8) :: TairK       !air temperature, [K]
  real(r8) :: CanWatg    !total canopy water content stored with dry matter, [m3 d-2]
  real(r8) :: Eco_Heat_Latent_col       !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8) :: Canopy_Heat_Latent_col      !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VcumIceSnow     !ice volume in snowpack, [m3 d-2]
  real(r8) :: TKSnow       !snow temperature, [K]
  real(r8) :: CanH2OHeldVg    !canopy surface water content, [m3 d-2]
  real(r8) :: Eco_Heat_Sens_col       !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8) :: RoughHeight        !canopy surface roughness height, [m]
  real(r8) :: ZeroPlanDisp        !zero plane displacement height, [m]
  real(r8) :: BndlResistAboveCanG       !isothermal boundary layer resistance, [h m-1]
  real(r8) :: RIB       !Richardson number for calculating boundary layer resistance, [-]
  real(r8) :: Eco_Heat_Grnd_col       !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8), pointer :: PTrans(:)     => null()    !canopy transpiration, [m2 d-2 h-1]
  real(r8), pointer :: PSICanPTurg(:)  => null()    !plant canopy turgor water potential, [MPa]
  real(r8), pointer :: PrecIntcptByCanP(:)   => null()    !water flux into canopy, [m3 d-2 h-1]
  real(r8), pointer :: PSICanP(:)  => null()    !canopy total water potential , [Mpa]
  real(r8), pointer :: VapXAir2PCan(:)  => null()    !canopy evaporation, [m2 d-2 h-1]
  real(r8), pointer :: HeatStorCanP(:)  => null()    !canopy storage heat flux, [MJ d-2 h-1]
  real(r8), pointer :: EvapTransHeatP(:)  => null()    !canopy latent heat flux, [MJ d-2 h-1]
  real(r8), pointer :: RAZ(:)    => null()    !canopy roughness height, [m]
  real(r8), pointer :: TKS(:)    => null()    !mean annual soil temperature, [K]
  real(r8), pointer :: PSICanPDailyMin(:)  => null()    !minimum daily canopy water potential, [MPa]
  real(r8), pointer :: TCelciusCanopy(:)    => null()    !canopy temperature, [oC]
  real(r8), pointer :: DTKC(:)   => null()    !change in canopy temperature, [K]
  real(r8), pointer :: ENGYX(:)  => null()    !canopy heat storage from previous time step, [MJ d-2]
  real(r8), pointer :: TKC(:)    => null()    !canopy temperature, [K]
  real(r8), pointer :: PSICanPOsmo(:)  => null()    !canopy osmotic water potential, [Mpa]
  real(r8), pointer :: OSMO(:)   => null()    !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8), pointer :: HeatXAir2PCan(:)  => null()    !canopy sensible heat flux, [MJ d-2 h-1]
  real(r8), pointer :: TKCZ(:)   => null()    !canopy temperature, [K]
  real(r8), pointer :: PSIRoot(:,:,:) => null() !root total water potential , [Mpa]
  real(r8), pointer :: PSIRootOSMO(:,:,:) => null() !root osmotic water potential , [Mpa]
  real(r8), pointer :: PSIRootTurg(:,:,:) => null() !root turgor water potential , [Mpa]
  real(r8), pointer :: AllPlantRootH2OUptake_vr(:,:,:) => null() !root water uptake, [m2 d-2 h-1]
  real(r8), pointer :: THeatRootUptake(:)     => null()   !total root heat uptake, [MJ d-2]
  real(r8), pointer :: GridPlantRootH2OUptake_vr(:)    => null()   !total root water uptake, [m3 d-2]
  real(r8), pointer :: WatByPCanopy(:)     => null()  !canopy surface water content, [m3 d-2]
  real(r8), pointer :: CanWatP(:)     => null()  !canopy water content, [m3 d-2]
  real(r8), pointer :: VHeatCapCanP(:)     => null()  !canopy heat capacity, [MJ d-2 K-1]
  real(r8), pointer :: TotalSoilH2OPSIMPa(:)     => null()  !soil micropore total water potential [MPa]
  real(r8), pointer :: ETCanP(:)     =>  null()  !total transpiration, [m H2O d-2]
  contains
    procedure, public :: Init => plt_ew_init
    procedure, public :: Destroy=> plt_ew_destroy
  end type plant_ew_type

  type, public :: plant_disturb_type

  real(r8) :: XCORP     !factor for surface litter incorporation and soil mixing
  integer  :: IYTYP      !fertilizer release type from fertilizer input file
  real(r8) :: DCORP      !soil mixing fraction with tillage, [-]
  integer  :: ITILL      !soil disturbance type, [-]
  real(r8), pointer :: XHVSTE(:)   => null()  !ecosystem harvest element, [gC d-2]
  real(r8), pointer :: VOXYF(:)    => null()  !plant O2 uptake from fire, [g d-2 ]
  integer,  pointer :: IDAYY(:)    => null()  !alternate day of harvest, [-]
  real(r8), pointer :: FERT(:)     => null()  !fertilizer application, [g m-2]
  integer,  pointer :: iYearPlanting(:)     => null()  !year of planting
  integer,  pointer :: IDAYX(:)    => null()  !alternate day of planting
  integer,  pointer :: IYRX(:)     => null()  !alternate year of planting
  integer,  pointer :: IYRY(:)     => null()  !alternate year of harvest
  integer,  pointer :: iDayPlanting(:)    => null()  !day of planting
  integer,  pointer :: iDayPlantHarvest(:)    => null()  !day of harvest
  integer,  pointer :: iYearPlantHarvest(:)     => null()  !year of harvest
  integer,  pointer :: IHVST(:)    => null()  !type of harvest
  integer,  pointer :: JHVST(:)    => null()  !flag for stand replacing disturbance
  real(r8), pointer :: EHVST(:,:,:)=> null()  !harvest efficiency, [-]
  real(r8), pointer :: HVST(:)     => null()  !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  real(r8), pointer :: THIN_pft(:)     => null()  !thinning of plant population, [-]
  real(r8), pointer :: CH4ByFire_pft(:)    => null()  !plant CH4 emission from fire, [g d-2 ]
  real(r8), pointer :: CO2ByFire_pft(:)    => null()  !plant CO2 emission from fire, [g d-2 ]
  real(r8), pointer :: VN2OF(:)    => null()  !plant N2O emission from fire, [g d-2 ]
  real(r8), pointer :: VNH3F(:)    => null()  !plant NH3 emission from fire, [g d-2 ]
  real(r8), pointer :: VPO4F(:)    => null()  !plant PO4 emission from fire, [g d-2 ]
  real(r8), pointer :: THVSTE(:,:) => null()  !total plant element harvest, [gC d-2 ]
  real(r8), pointer :: HVSTE(:,:)  => null()  !plant element harvest, [g d-2 ]
  contains
    procedure, public :: Init    =>  plt_disturb_init
    procedure, public :: Destroy => plt_disturb_destroy
  end type plant_disturb_type

  type, public :: plant_bgcrate_type
  real(r8) :: Eco_NBP_col      !total NBP, [g d-2]
  real(r8) :: Eco_GPP_col      !ecosystem GPP, [g d-2 h-1]

  real(r8) :: CNETX     !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8) :: ECO_ER_col      !ecosystem respiration, [g d-2 h-1]
  real(r8) :: Eco_AutoR_col      !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8) :: TH2GZ     !total root H2 flux, [g d-2]
  real(r8) :: Canopy_NEE_col     !total net CO2 fixation, [gC d-2]
  real(r8), pointer :: LitterFallChemElmnt_col(:) => null() !total litterfall element, [g d-2 h-1]
  real(r8), pointer :: NetPrimaryProductvity_pft(:)       => null()   !total net primary productivity, [gC d-2]
  real(r8), pointer :: RNH3C(:)      => null()   !canopy NH3 flux, [g d-2 h-1]
  real(r8), pointer :: TDFOME(:,:,:)   =>  null()  !total root element exchange, [g d-2 h-1]
  real(r8), pointer :: RootN2Fix_pvr(:,:)    =>  null()  !root N2 fixation, [gN d-2 h-1]
  real(r8), pointer :: CanopyPlusNoduRespC_pft(:)      =>  null()  !total autotrophic respiration, [gC d-2 ]
  real(r8), pointer :: LitterFallChemElmnt_pftvr(:,:,:,:,:) =>  null()  !plant litterfall element, [g d-2 h-1]
  real(r8), pointer :: ROXYX(:)      =>  null()  !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8), pointer :: RNHBX(:)      => null()   !total root + microbial NH4 uptake band, [gN d-2 h-1]
  real(r8), pointer :: RP14X(:)      => null()   !HPO4 demand in non-band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), pointer :: RP1BX(:)      => null()   !HPO4 demand in band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), pointer :: RNO3X(:)      => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4X(:)      => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RN3BX(:)      => null()   !total root + microbial NO3 uptake band, [gN d-2 h-1]
  real(r8), pointer :: RPO4X(:)      => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: RPOBX(:)      => null()   !total root + microbial PO4 uptake band, [gP d-2 h-1]
  real(r8), pointer :: RPO4Y(:)      => null()   !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), pointer :: RPOBY(:)      => null()   !total root + microbial PO4 uptake band, [gP d-2 h-1]
  real(r8), pointer :: RP14Y(:)      => null()   !HPO4 demand in non-band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), pointer :: RP1BY(:)      => null()   !HPO4 demand in band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), pointer :: RNO3Y(:)      => null()   !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNH4Y(:)      => null()   !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), pointer :: RNHBY(:)      => null()   !total root + microbial NH4 uptake band, [gN d-2 h-1]
  real(r8), pointer :: RN3BY(:)      => null()   !total root + microbial NO3 uptake band, [gN d-2 h-1]
  real(r8), pointer :: ROXYF(:)      => null()   !net gaseous O2 flux, [g d-2 h-1]
  real(r8), pointer :: RCO2F(:)      => null()   !net gaseous CO2 flux, [g d-2 h-1]
  real(r8), pointer :: ROXYL(:)      => null()   !net aqueous O2 flux, [g d-2 h-1]
  real(r8), pointer :: ROXYY(:)      => null()   !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8), pointer :: TCO2P(:)      => null()   !total root CO2 flux, [gC d-2 h-1]
  real(r8), pointer :: TUPOXP(:)     => null()   !total root internal O2 flux, [g d-2 h-1]
  real(r8), pointer :: LitrfalChemElemnts_vr(:,:,:,:)   => null() !total litterfall element, [g d-2 h-1]
  real(r8), pointer :: RDOM_micb_flx(:,:,:)    => null()  !net microbial DOC flux, [gC d-2 h-1]
  real(r8), pointer :: CO2NetFix_pft(:)       => null()  !canopy net CO2 exchange, [gC d-2 h-1]
  real(r8), pointer :: GrossCO2Fix_pft(:)      => null()  !total gross CO2 fixation, [gC d-2 ]
  real(r8), pointer :: LitterFallChemElmnt_pft(:,:)    => null()  !plant element litterfall, [g d-2 h-1]
  real(r8), pointer :: RootGasLoss_disturb(:,:)=> null() !gaseous flux fron root disturbance, [g d-2 h-1]
  real(r8), pointer :: SurfLitrfallChemElmnts_pft(:,:)    => null()  !total surface litterfall element, [g d-2]
  real(r8), pointer :: LitrfallChemElmnts_pft(:,:)    => null()  !total plant element litterfall , [g d-2 ]
  real(r8), pointer :: GrossResp_pft(:)      => null()  !total plant respiration, [gC d-2 ]

  real(r8), pointer :: NH3EmiCum_pft(:)      => null()  !total canopy NH3 flux, [gN d-2 ]
  real(r8), pointer :: PlantN2FixCum_pft(:)     => null()  !total plant N2 fixation, [g d-2 ]

  contains
    procedure, public :: Init  => plt_bgcrate_init
    procedure, public :: Destroy  => plt_bgcrate_destroy
  end type plant_bgcrate_type

  type, public :: plant_rootbgc_type
  real(r8), pointer :: TRootGasLoss_disturb(:)   => null()  !total root gas content [g d-2]
  real(r8), pointer :: trcs_plant_uptake_vr(:,:)     => null()   !total root-soil solute flux non-band, [g d-2 h-1]  
  real(r8), pointer :: RootExudChemElmnt_pft(:,:)       => null()  !total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
  real(r8), pointer :: RootN2Fix_pft(:)          => null()  !total root N2 fixation, [g d-2 h-1]
  real(r8), pointer :: RootNO3Uptake_pft(:)         => null()  !total root uptake of NO3, [g d-2 h-1]
  real(r8), pointer :: RootNH4Uptake_pft(:)         => null()  !total root uptake of NH4, [g d-2 h-1]
  real(r8), pointer :: RootHPO4Uptake_pft(:)         => null()  !total root uptake of HPO4, [g d-2 h-1]
  real(r8), pointer :: RootH2PO4Uptake_pft(:)         => null()  !total root uptake of PO4, [g d-2 h-1]
  real(r8), pointer :: RDFOME(:,:,:,:,:)  => null()  !root uptake (+ve) - exudation (-ve) of DOE, [g d-2 h-1]
  real(r8), pointer :: PlantRootSoilChemNetX_pft(:,:)      => null()  !net root element uptake (+ve) - exudation (-ve), [gC d-2 h-1]
  real(r8), pointer :: ROXSK(:,:)       => null()  !total O2 sink, [g d-2 t-1]
  real(r8), pointer :: ZEROQ(:)         => null()  !threshold zero for uptake calculation
  real(r8), pointer :: UPMNPO(:,:)      => null()  !minimum PO4 concentration for root NH4 uptake, [g m-3]
  real(r8), pointer :: UPMXPO(:,:)      => null()  !maximum root PO4 uptake rate, [g m-2 h-1]
  real(r8), pointer :: UPKMPO(:,:)      => null()  !Km for root PO4 uptake, [g m-3]
  real(r8), pointer :: UPMNZO(:,:)      => null()  !minimum NO3 concentration for root NH4 uptake, [g m-3]
  real(r8), pointer :: UPMXZO(:,:)      => null()  !maximum root NO3 uptake rate, [g m-2 h-1]
  real(r8), pointer :: UPKMZO(:,:)      => null()  !Km for root NO3 uptake, [g m-3]
  real(r8), pointer :: UPMNZH(:,:)      => null()  !minimum NH4 concentration for root NH4 uptake, [g m-3]
  real(r8), pointer :: UPMXZH(:,:)      => null()  !maximum root NH4 uptake rate, [g m-2 h-1]
  real(r8), pointer :: UPKMZH(:,:)      => null()  !Km for root NH4 uptake, [g m-3]
  real(r8), pointer :: RCO2P(:,:,:)     => null()  !aqueous CO2 flux from roots to root water , [g d-2 h-1]
  real(r8), pointer :: RUPOXP(:,:,:)    => null()  !aqueous O2 flux from roots to root water , [g d-2 h-1]
  real(r8), pointer :: RCO2S(:,:,:)     => null()  !aqueous CO2 flux from roots to soil water, [g d-2 h-1]
  real(r8), pointer :: RUPOXS(:,:,:)    => null()  !aqueous O2 flux from roots to soil water, [g d-2 h-1]
  real(r8), pointer :: RUPCHS(:,:,:)    => null()  !aqueous CH4 flux from roots to soil water, [g d-2 h-1]
  real(r8), pointer :: RUPN2S(:,:,:)    => null()  !aqueous N2O flux from roots to soil water, [g d-2 h-1]
  real(r8), pointer :: RUPN3S(:,:,:)    => null()  !aqueous NH3 flux from roots to soil water non-band, [g d-2 h-1]
  real(r8), pointer :: RUPN3B(:,:,:)    => null()  !aqueous NH3 flux from roots to soil water band, [g d-2 h-1]
  real(r8), pointer :: RUPHGS(:,:,:)    => null()  !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx_pft_vr(:,:,:,:)    => null()  !gaseous tracer flux through roots, [g d-2 h-1]
  real(r8), pointer :: trcg_Root_DisEvap_flx_vr(:,:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
  real(r8), pointer :: ROXYP(:,:,:)     => null()  !root  O2 demand from respiration, [g d-2 h-1]
  real(r8), pointer :: RUNNHP(:,:,:)    => null()  !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8), pointer :: RUNNBP(:,:,:)    => null()  !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8), pointer :: RUNNOP(:,:,:)    => null()  !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8), pointer :: RUNNXP(:,:,:)    => null()  !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8), pointer :: RUPP2P(:,:,:)    => null()  !root uptake of H2PO4 non-band
  real(r8), pointer :: RUPP2B(:,:,:)    => null()  !root uptake of H2PO4 band
  real(r8), pointer :: RUPP1P(:,:,:)    => null()  !HPO4 demand in non-band by each root population
  real(r8), pointer :: RUPP1B(:,:,:)    => null()  !HPO4 demand in band by each root population
  real(r8), pointer :: RootAutoRO2Limiter_pvr(:,:,:)       => null()  !O2 constraint to root respiration, []
  real(r8), pointer :: trcg_rootml(:,:,:,:)=> null() !root gas content, [g d-2]
  real(r8), pointer :: trcs_rootml(:,:,:,:)=> null() !root aqueous content, [g d-2]
  real(r8), pointer :: RNH3B(:,:)       => null()  !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  real(r8), pointer :: RUPNHB(:,:,:)    => null()  !root uptake of NH4 band, [g d-2 h-1]
  real(r8), pointer :: RootNH4Uptake_pvr(:,:,:)    => null()  !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8), pointer :: RootH2PO4Uptake_pvr(:,:,:)    => null()  !root uptake of PO4 non-band, [g d-2 h-1]
  real(r8), pointer :: RUPNOB(:,:,:)    => null()  !root uptake of NO3 band, [g d-2 h-1]
  real(r8), pointer :: RootNO3Uptake_pvr(:,:,:)    => null()  !root uptake of NO3 non-band, [g d-2 h-1]
  real(r8), pointer :: RUPH1B(:,:,:)    => null()  !root uptake of HPO4 band
  real(r8), pointer :: RootHPO4Uptake_pvr(:,:,:)    => null()  !root uptake HPO4 non-band
  real(r8), pointer :: RUPH2B(:,:,:)    => null()  !root uptake of PO4 band, [g d-2 h-1]
  real(r8), pointer :: RUONHB(:,:,:)    => null()  !root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RUONH4(:,:,:)    => null()  !root uptake of NH4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RUOH2P(:,:,:)    => null()  !root uptake of PO4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RUONOB(:,:,:)    => null()  !root uptake of NO3 band unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RUONO3(:,:,:)    => null()  !root uptake of NO3 non-band unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RUOH1B(:,:,:)    => null()  !root HPO4 uptake in band unlimited by O2
  real(r8), pointer :: RUOH1P(:,:,:)    => null()  !root HPO4 uptake in non-band unlimited by O2
  real(r8), pointer :: RUOH2B(:,:,:)    => null()  !root uptake of PO4 band unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RUCNHB(:,:,:)    => null()  !root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RUCNH4(:,:,:)    => null()  !root uptake of NH4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RUCH2P(:,:,:)    => null()  !root uptake of PO4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RUCNOB(:,:,:)    => null()  !root uptake of NO3 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RUCNO3(:,:,:)    => null()  !root uptake of NO3 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RUCH1B(:,:,:)    => null()  !root HPO4 uptake in band unlimited by nonstructural C
  real(r8), pointer :: RUCH1P(:,:,:)    => null()  !root HPO4 uptake in non-band unlimited by nonstructural C
  real(r8), pointer :: RUCH2B(:,:,:)    => null()  !root uptake of PO4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RootRespPotential_vr(:,:,:)     => null()  !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RCO2N(:,:,:)     => null()  !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RCO2A(:,:,:)     => null()  !root respiration constrained by O2, [g d-2 h-1]
  real(r8), pointer :: PlantExudChemElmntCum_pft(:,:)      => null()  !total net root element uptake (+ve) - exudation (-ve), [gC d-2 ]
  real(r8), pointer :: trcg_TLP(:,:)    => null()   !total root internal gas flux, [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx_vr(:,:)   => null()   !total internal root gas flux , [gC d-2 h-1]

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
  allocate(this%trcs_plant_uptake_vr(ids_beg:ids_end,JZ1))
  allocate(this%trcg_rootml(idg_beg:idg_end-1,2,JZ1,JP1));this%trcg_rootml=0._r8
  allocate(this%trcs_rootml(idg_beg:idg_end-1,2,JZ1,JP1));this%trcs_rootml=0._r8
  allocate(this%TRootGasLoss_disturb(idg_beg:idg_end-1));this%TRootGasLoss_disturb=0._r8
  allocate(this%ROXSK(60,0:JZ1))
  allocate(this%RDFOME(NumOfPlantChemElmnts,2,1:jcplx,0:JZ1,JP1))
  allocate(this%PlantRootSoilChemNetX_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%PlantExudChemElmntCum_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%RootExudChemElmnt_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%RootN2Fix_pft(JP1))
  allocate(this%RootNO3Uptake_pft(JP1))
  allocate(this%RootNH4Uptake_pft(JP1))
  allocate(this%RootHPO4Uptake_pft(JP1))
  allocate(this%RootH2PO4Uptake_pft(JP1))

  allocate(this%ZEROQ(JP1))
  allocate(this%RootRespPotential_vr(2,JZ1,JP1))
  allocate(this%RCO2N(2,JZ1,JP1))
  allocate(this%RCO2A(2,JZ1,JP1))
  allocate(this%RUPNHB(2,JZ1,JP1))
  allocate(this%RootNH4Uptake_pvr(2,JZ1,JP1))
  allocate(this%RootH2PO4Uptake_pvr(2,JZ1,JP1))
  allocate(this%RUPNOB(2,JZ1,JP1))
  allocate(this%RootNO3Uptake_pvr(2,JZ1,JP1))
  allocate(this%RUPH1B(2,JZ1,JP1))
  allocate(this%RootHPO4Uptake_pvr(2,JZ1,JP1))
  allocate(this%RUPH2B(2,JZ1,JP1))
  allocate(this%RUONHB(2,JZ1,JP1))
  allocate(this%RUONH4(2,JZ1,JP1))
  allocate(this%RUOH2P(2,JZ1,JP1))
  allocate(this%RUONOB(2,JZ1,JP1))
  allocate(this%RUONO3(2,JZ1,JP1))
  allocate(this%RUOH1B(2,JZ1,JP1))
  allocate(this%RUOH1P(2,JZ1,JP1))
  allocate(this%RUOH2B(2,JZ1,JP1))
  allocate(this%RUCNHB(2,JZ1,JP1))
  allocate(this%RUCNH4(2,JZ1,JP1))
  allocate(this%RUCH2P(2,JZ1,JP1))
  allocate(this%RUCNOB(2,JZ1,JP1))
  allocate(this%RUCNO3(2,JZ1,JP1))
  allocate(this%RUCH1B(2,JZ1,JP1))
  allocate(this%RUCH1P(2,JZ1,JP1))
  allocate(this%RUCH2B(2,JZ1,JP1))
  allocate(this%RNH3B(MaxNumBranches,JP1))


  allocate(this%trcg_air2root_flx_vr(idg_beg:idg_end-1,JZ1))
  allocate(this%trcg_TLP(idg_beg:idg_end-1,JZ1))

  allocate(this%trcg_air2root_flx_pft_vr(idg_beg:idg_end-1,2,JZ1,JP1))
  allocate(this%trcg_Root_DisEvap_flx_vr(idg_beg:idg_end-1,2,JZ1,JP1))
  allocate(this%ROXYP(2,JZ1,JP1))
  allocate(this%RUNNHP(2,JZ1,JP1))
  allocate(this%RUNNBP(2,JZ1,JP1))
  allocate(this%RUNNOP(2,JZ1,JP1))
  allocate(this%RUNNXP(2,JZ1,JP1))
  allocate(this%RUPP2P(2,JZ1,JP1))
  allocate(this%RUPP2B(2,JZ1,JP1))
  allocate(this%RUPP1P(2,JZ1,JP1))
  allocate(this%RUPP1B(2,JZ1,JP1))

  allocate(this%RootAutoRO2Limiter_pvr(2,JZ1,JP1))
  allocate(this%UPMNPO(2,JP1))
  allocate(this%UPMXPO(2,JP1))
  allocate(this%UPKMPO(2,JP1))
  allocate(this%UPMNZO(2,JP1))
  allocate(this%UPMXZO(2,JP1))
  allocate(this%UPKMZO(2,JP1))
  allocate(this%UPMNZH(2,JP1))
  allocate(this%UPMXZH(2,JP1))
  allocate(this%UPKMZH(2,JP1))
  allocate(this%RCO2P(2,JZ1,JP1))
  allocate(this%RUPOXP(2,JZ1,JP1))
  allocate(this%RCO2S(2,JZ1,JP1))
  allocate(this%RUPOXS(2,JZ1,JP1))
  allocate(this%RUPCHS(2,JZ1,JP1))
  allocate(this%RUPN2S(2,JZ1,JP1))
  allocate(this%RUPN3S(2,JZ1,JP1))
  allocate(this%RUPN3B(2,JZ1,JP1))
  allocate(this%RUPHGS(2,JZ1,JP1))
  end subroutine plt_rootbgc_init
!----------------------------------------------------------------------

  subroutine plt_rootbgc_destroy(this)

  implicit none
  class(plant_rootbgc_type) :: this

!  call destroy(this%trcg_rootml)
!  call destroy(this%trcs_rootml)
!  if(allocated(CO2P))deallocate(CO2P)
!  if(allocated(CO2A))deallocate(CO2A)
!  if(allocated(H2GP))deallocate(H2GP)
!  if(allocated(H2GA))deallocate(H2GA)
!  if(allocated(OXYP))deallocate(OXYP)
!  if(allocated(OXYA))deallocate(OXYA)
!  if(allocated(TRootN2Fix_pft))deallocate(TRootN2Fix_pft)
!  if(allocated(ROXSK))deallocate(ROXSK)
!  if(allocated(RDFOME))deallocate(RDFOME)
!  if(allocated(PlantRootSoilChemNetX_pft))deallocate(PlantRootSoilChemNetX_pft)
!  if(allocated(PlantExudChemElmntCum_pft))deallocate(PlantExudChemElmntCum_pft)
!  if(allocated(RootExudChemElmnt_pft))deallocate(RootExudChemElmnt_pft)
!  if(allocated(RootN2Fix_pft))deallocate(RootN2Fix_pft)
!  if(allocated(RootNO3Uptake_pft))deallocate(RootNO3Uptake_pft)
!  if(allocated(RootNH4Uptake_pft))deallocate(RootNH4Uptake_pft)
!  if(allocated(RootHPO4Uptake_pft))deallocate(RootHPO4Uptake_pft)
!  if(allocated(RootH2PO4Uptake_pft))deallocate(RootH2PO4Uptake_pft)
!  if(allocated(RNH3B))deallocate(RNH3B)
!  if(allocated(ZEROQ))deallocate(ZEROQ)
!  if(allocated(RootRespPotential_vr))deallocate(RootRespPotential_vr)
!  if(allocated(RCO2N))deallocate(RCO2N)
!  if(allocated(RCO2A))deallocate(RCO2A)
!  if(allocated(RCO2P))deallocate(RCO2P)
!  if(allocated(RUPOXP))deallocate(RUPOXP)
!  if(allocated(RCO2S))deallocate(RCO2S)
!  if(allocated(RUPOXS))deallocate(RUPOXS)
!  if(allocated(RUPCHS))deallocate(RUPCHS)
!  if(allocated(RUPN2S))deallocate(RUPN2S)
!  if(allocated(RUPN3S))deallocate(RUPN3S)
!  if(allocated(RUPN3B))deallocate(RUPN3B)
!  if(allocated(RUPHGS))deallocate(RUPHGS)
!  if(allocated(RCOFLA))deallocate(RCOFLA)
!  if(allocated(ROXFLA))deallocate(ROXFLA)
!  if(allocated(RCHFLA))deallocate(RCHFLA)

!  if(allocated(RUPNHB))deallocate(RUPNHB)
!  if(allocated(RootNH4Uptake_pvr))deallocate(RootNH4Uptake_pvr)
!  if(allocated(RootH2PO4Uptake_pvr))deallocate(RootH2PO4Uptake_pvr)
!  if(allocated(RUPNOB))deallocate(RUPNOB)
!  if(allocated(RootNO3Uptake_pvr))deallocate(RootNO3Uptake_pvr)
!  if(allocated(RUPH1B))deallocate(RUPH1B)
!  if(allocated(RootHPO4Uptake_pvr))deallocate(RootHPO4Uptake_pvr)
!  if(allocated(RUPH2B))deallocate(RUPH2B)
!  if(allocated(RUONHB))deallocate(RUONHB)
!  if(allocated(RUONH4))deallocate(RUONH4)
!  if(allocated(RUOH2P))deallocate(RUOH2P)
!  if(allocated(RUONOB))deallocate(RUONOB)
!  if(allocated(RUONO3))deallocate(RUONO3)
!  if(allocated(RUOH1B))deallocate(RUOH1B)
!  if(allocated(RUOH1P))deallocate(RUOH1P)
!  if(allocated(RUOH2B))deallocate(RUOH2B)
!  if(allocated(RUCNHB))deallocate(RUCNHB)
!  if(allocated(RUCNH4))deallocate(RUCNH4)
!  if(allocated(RUCH2P))deallocate(RUCH2P)
!  if(allocated(RUCNOB))deallocate(RUCNOB)
!  if(allocated(RUCNO3))deallocate(RUCNO3)
!  if(allocated(RUCH1B))deallocate(RUCH1B)
!  if(allocated(RUCH1P))deallocate(RUCH1P)
!  if(allocated(RUCH2B))deallocate(RUCH2B)



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
!  if(allocated(ROXYP))deallocate(ROXYP)
!  if(allocated(RUNNHP))deallocate(RUNNHP)
!  if(allocated(RUNNBP))deallocate(RUNNBP)
!  if(allocated(RUNNOP))deallocate(RUNNOP)
!  if(allocated(RUNNXP))deallocate(RUNNXP)
!  if(allocated(RUPP2P))deallocate(RUPP2P)
!  if(allocated(RUPP2B))deallocate(RUPP2B)
!  if(allocated(RUPP1P))deallocate(RUPP1P)
!  if(allocated(RUPP1B))deallocate(RUPP1B)
!  if(allocated(RootAutoRO2Limiter_pvr))deallocate(RootAutoRO2Limiter_pvr)
!  if(allocated(UPMNPO))deallocate(UPMNPO)
!  if(allocated(UPMXPO))deallocate(UPMXPO)
!  if(allocated(UPKMPO))deallocate(UPKMPO)
!  if(allocated(UPMNZO))deallocate(UPMNZO)
!  if(allocated(UPMXZO))deallocate(UPMXZO)
!  if(allocated(UPKMZO))deallocate(UPKMZO)
!  if(allocated(UPMNZH))deallocate(UPMNZH)
!  if(allocated(UPMXZH))deallocate(UPMXZH)
!  if(allocated(UPKMZH))deallocate(UPKMZH)
  end subroutine plt_rootbgc_destroy
!----------------------------------------------------------------------
  subroutine plt_site_Init(this)
  implicit none
  class(plant_siteinfo_type) :: this


  allocate(this%PlantElemntStoreLandscape(NumOfPlantChemElmnts))
  allocate(this%FracSoiAsMicP(0:JZ1))
  allocate(this%DATAP(JP1))
  allocate(this%DATA(30))
  allocate(this%AREA3(0:JZ1))
  allocate(this%IDATA(60))
  allocate(this%DLYR3(0:JZ1))
  allocate(this%ElmntBalanceCum_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%CumSoilThickness(0:JZ1))
  allocate(this%DPTHZ(0:JZ1))
  allocate(this%PPI(JP1))
  allocate(this%PPZ(JP1))
  allocate(this%PPX(JP1))
  allocate(this%pftPlantPopulation(JP1))
  allocate(this%VLWatMicPM(60,0:JZ1))
  allocate(this%VLsoiAirPM(60,0:JZ1))
  allocate(this%TortMicPM(60,0:JZ1))
  allocate(this%FILM(60,0:JZ1))

  end subroutine plt_site_Init
!----------------------------------------------------------------------
  subroutine plt_site_destroy(this)
  implicit none
  class(plant_siteinfo_type) :: this



!  if(allocated(DLYR3))deallocate(DLYR3)


!  if(allocated(FracSoiAsMicP))deallocate(FracSoiAsMicP)
!  if(allocated(TortMicPM))deallocate(TortMicPM)
!  if(allocated(FILM))deallocate(FILM)
!  if(allocated(DATAP))deallocate(DATAP)
!  if(allocated(DATA))deallocate(DATA)
!  call Destroy(this%PlantElemntStoreLandscape)
!  if(allocated(AREA3))deallocate(AREA3)
!  if(allocated(IDATA))deallocate(IDATA)
!  if(allocated(ElmntBalanceCum_pft))deallocate(ElmntBalanceCum_pft)
!  if(allocated(BALP))deallocate(BALP)
!  if(allocated(CumSoilThickness))deallocate(CumSoilThickness)
!  if(allocated(DPTHZ))deallocate(DPTHZ)
!  if(allocated(PPI)) deallocate(PPI)
!  if(allocated(PPZ)) deallocate(PPZ)
!  if(allocated(PPX)) deallocate(PPX)
!  if(allocated(PP))deallocate(PP)
!  if(allocated(VLWatMicPM))deallocate(VLWatMicPM)
!  if(allocated(VLsoiAirPM))deallocate(VLsoiAirPM)

  end subroutine plt_site_destroy
!----------------------------------------------------------------------

  subroutine plt_bgcrate_init(this)
  implicit none
  class(plant_bgcrate_type) :: this

  allocate(this%TCO2P(JZ1))
  allocate(this%TUPOXP(JZ1))
  allocate(this%RPO4Y(0:JZ1))
  allocate(this%RPOBY(0:JZ1))
  allocate(this%RP14Y(0:JZ1))
  allocate(this%RP1BY(0:JZ1))
  allocate(this%RNO3Y(0:JZ1))
  allocate(this%RNH4Y(0:JZ1))
  allocate(this%RNHBY(0:JZ1))
  allocate(this%RN3BY(0:JZ1))
  allocate(this%ROXYF(0:JZ1))
  allocate(this%RCO2F(0:JZ1))
  allocate(this%ROXYL(0:JZ1))
  allocate(this%ROXYY(0:JZ1))
  allocate(this%LitrfalChemElemnts_vr(NumOfPlantChemElmnts,jsken,NumOfPlantLitrCmplxs,0:JZ1))
  allocate(this%GrossCO2Fix_pft(JP1))
  allocate(this%RDOM_micb_flx(idom_beg:idom_end,1:jcplx,0:JZ1))
  allocate(this%CO2NetFix_pft(JP1))
  allocate(this%RootGasLoss_disturb(idg_beg:idg_end-1,JP1))
  allocate(this%GrossResp_pft(JP1))
  allocate(this%PlantN2FixCum_pft(JP1))
  allocate(this%NH3EmiCum_pft(JP1))
  allocate(this%SurfLitrfallChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%LitterFallChemElmnt_col(NumOfPlantChemElmnts))
  allocate(this%NetPrimaryProductvity_pft(JP1))
  allocate(this%RNH3C(JP1))
  allocate(this%TDFOME(NumOfPlantChemElmnts,1:jcplx,JZ1))
  allocate(this%RootN2Fix_pvr(JZ1,JP1))
  allocate(this%CanopyPlusNoduRespC_pft(JP1))
  allocate(this%RP1BX(0:JZ1))
  allocate(this%RNO3X(0:JZ1))
  allocate(this%RPO4X(0:JZ1))
  allocate(this%RNH4X(0:JZ1))
  allocate(this%ROXYX(0:JZ1))
  allocate(this%RPOBX(0:JZ1))
  allocate(this%RN3BX(0:JZ1))
  allocate(this%RNHBX(0:JZ1))
  allocate(this%RP14X(0:JZ1))

  allocate(this%LitterFallChemElmnt_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%LitrfallChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%LitterFallChemElmnt_pftvr(NumOfPlantChemElmnts,jsken,1:NumOfPlantLitrCmplxs,0:JZ1,JP1))


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

  allocate(this%THIN_pft(JP1))
  allocate(this%XHVSTE(NumOfPlantChemElmnts))
  allocate(this%THVSTE(NumOfPlantChemElmnts,JP1))
  allocate(this%HVSTE(NumOfPlantChemElmnts,JP1))
  allocate(this%CH4ByFire_pft(JP1))
  allocate(this%CO2ByFire_pft(JP1))
  allocate(this%VN2OF(JP1))
  allocate(this%VNH3F(JP1))
  allocate(this%VPO4F(JP1))

  allocate(this%IDAYY(JP1))
  allocate(this%VOXYF(JP1))
  allocate(this%EHVST(1:2,1:4,JP1))
  allocate(this%HVST(JP1))
  allocate(this%iYearPlantHarvest(JP1))
  allocate(this%FERT(1:20))
  allocate(this%iYearPlanting(JP1))
  allocate(this%IDAYX(JP1))
  allocate(this%IYRX(JP1))
  allocate(this%IYRY(JP1))
  allocate(this%iDayPlanting(JP1))
  allocate(this%iDayPlantHarvest(JP1))
  allocate(this%IHVST(JP1))
  allocate(this%JHVST(JP1))

  end subroutine plt_disturb_init
!----------------------------------------------------------------------
  subroutine plt_disturb_destroy(this)
  implicit none
  class(plant_disturb_type) :: this


!  if(allocated(THIN_pft))deallocate(THIN_pft)
!  if(allocated(HVSTE))deallocate(HVSTE)
!  if(allocated(THVSTE))deallocate(THVSTE)
!  if(allocated(CH4ByFire_pft))deallocate(CH4ByFire_pft)
!  if(allocated(CO2ByFire_pft))deallocate(CO2ByFire_pft)
!  if(allocated(VN2OF))deallocate(VN2OF)
!  if(allocated(VNH3F))deallocate(VNH3F)
!  if(allocated(VPO4F))deallocate(VPO4F)

!  if(allocated(IDAYY))deallocate(IDAYY)
!  if(allocated(VOXYF))deallocate(VOXYF)
!  if(allocated(EHVST))deallocate(EHVST)
!  if(allocated(HVST))deallocate(HVST)
!  if(allocated(iYearPlantHarvest))deallocate(iYearPlantHarvest)
!  if(allocated(iDayPlanting))deallocate(iDayPlanting)
!  if(allocated(iDayPlantHarvest))deallocate(iDayPlantHarvest)
!  if(allocated(IYRY))deallocate(IYRY)
!  if(allocated(IDAYX))deallocate(IDAYX)
!  if(allocated(iYearPlanting))deallocate(iYearPlanting)
!  if(allocated(FERT))deallocate(FERT)
!  if(allocated(IYRX))deallocate(IYRX)
!  if(allocated(IHVST))deallocate(IHVST)
!  if(allocated(JHVST))deallocate(JHVST)

  end subroutine plt_disturb_destroy
!----------------------------------------------------------------------
  subroutine  plt_ew_init(this)

  implicit none
  class(plant_ew_type) :: this

  allocate(this%ETCanP(JP1))
  allocate(this%TotalSoilH2OPSIMPa(0:JZ1))
  allocate(this%THeatRootUptake(0:JZ1))
  allocate(this%TKCZ(JP1))
  allocate(this%HeatXAir2PCan(JP1))
  allocate(this%PrecIntcptByCanP(JP1))
  allocate(this%PSICanPTurg(JP1))
  allocate(this%PSICanP(JP1))
  allocate(this%VapXAir2PCan(JP1))
  allocate(this%HeatStorCanP(JP1))
  allocate(this%EvapTransHeatP(JP1))
  allocate(this%WatByPCanopy(JP1))
  allocate(this%VHeatCapCanP(JP1))
  allocate(this%CanWatP(JP1))
  allocate(this%PSIRoot(jroots,JZ1,JP1))
  allocate(this%PSIRootOSMO(jroots,JZ1,JP1))
  allocate(this%PSIRootTurg(jroots,JZ1,JP1))
  allocate(this%AllPlantRootH2OUptake_vr(jroots,JZ1,JP1))
  allocate(this%GridPlantRootH2OUptake_vr(0:JZ1))
  allocate(this%PTrans(JP1))
  allocate(this%PSICanPOsmo(JP1))
  allocate(this%TKS(0:JZ1))
  allocate(this%OSMO(JP1))
  allocate(this%RAZ(JP1))
  allocate(this%DTKC(JP1))
  allocate(this%TKC(JP1))
  allocate(this%ENGYX(JP1))
  allocate(this%TCelciusCanopy(JP1))
  allocate(this%PSICanPDailyMin(JP1))

  end subroutine plt_ew_init
!----------------------------------------------------------------------

  subroutine plt_ew_destroy(this)
  implicit none
  class(plant_ew_type) :: this

!  if(allocated(ETCanP))deallocate(ETCanP)
!  if(allocated(PSIST))deallocate(PSIST)
!  if(allocated(THeatRootUptake))deallocate(THeatRootUptake)
!  if(allocated(TKCZ))deallocate(TKCZ)
!  if(allocated(HeatXAir2PCan))deallocate(HeatXAir2PCan)
!  if(allocated(PrecIntcptByCanP))deallocate(PrecIntcptByCanP)
!  if(allocated(PSICanPTurg))deallocate(PSICanPTurg)
!  if(allocated(PSICanP))deallocate(PSICanP)
!  if(allocated(VapXAir2PCan))deallocate(VapXAir2PCan)
!  if(allocated(HeatStorCanP))deallocate(HeatStorCanP)
!  if(allocated(EvapTransHeatP))deallocate(EvapTransHeatP)
!  if(allocated(VHeatCapCanP))deallocate(VHeatCapCanP)
!  if(allocated(WatByPCanopy))deallocate(WatByPCanopy)
!  if(allocated(CanWatP))deallocate(CanWatP)
!  if(allocated(PSIRoot))deallocate(PSIRoot)
!  if(allocated(PSIRootOSMO))deallocate(PSIRootOSMO)
!  if(allocated(PSIRootTurg))deallocate(PSIRootTurg)
!  if(allocated(AllPlantRootH2OUptake_vr))deallocate(AllPlantRootH2OUptake_vr)
!  if(allocated(TAllPlantRootH2OUptake_vr))deallocate(TAllPlantRootH2OUptake_vr)
!  if(allocated(PTrans))deallocate(PTrans)
!  if(allocated(PSICanPOsmo))deallocate(PSICanPOsmo)
!  if(allocated(TKS))deallocate(TKS)
!  if(allocated(PSICanPDailyMin))deallocate(PSICanPDailyMin)
!  if(allocated(RAZ))deallocate(RAZ)
!  if(allocated(OSMO))deallocate(OSMO)
!  if(allocated(DTKC))deallocate(DTKC)
!  if(allocated(TKC))deallocate(TKC)
!  if(allocated(ENGYX))deallocate(ENGYX)
!  if(allocated(TCelciusCanopy))deallocate(TCelciusCanopy)

  end subroutine plt_ew_destroy

!----------------------------------------------------------------------

  subroutine plt_allom_init(this)
  implicit none
  class(plant_allometry_type) :: this


  allocate(this%RootFracRemobilizableBiom(JP1))
  allocate(this%FNOD(JP1))
  allocate(this%GrainSeedBiomCMean_brch(MaxNumBranches,JP1))
  allocate(this%NoduGrowthYield_pft(JP1))
  allocate(this%RootBiomGrowthYield(JP1))
  allocate(this%RootrPC_pft(JP1))
  allocate(this%rCNNonstructRemob_pft(JP1))
  allocate(this%rCPNonstructRemob_pft(JP1))
  allocate(this%CPRTS(JP1))
  allocate(this%CNRTS(JP1))
  allocate(this%NodulerNC_pft(JP1))
  allocate(this%NodulerPC_pft(JP1))
  allocate(this%RootrNC_pft(JP1))
  allocate(this%CPLF(JP1))
  allocate(this%CPSHE(JP1))
  allocate(this%CNLF(JP1))
  allocate(this%CNSHE(JP1))
  allocate(this%CNGR(JP1))
  allocate(this%rPCStalk_pft(JP1))
  allocate(this%rNCStalk_pft(JP1))
  allocate(this%CPGR(JP1))
  allocate(this%rPCEar_pft(JP1))
  allocate(this%rPCReserve_pft(JP1))
  allocate(this%rNCReserve_pft(JP1))
  allocate(this%rPCHusk_pft(JP1))
  allocate(this%FVRN(0:5))
  allocate(this%FWOODE(NumOfPlantChemElmnts,NumOfPlantLitrCmplxs))
  allocate(this%FWODLE(NumOfPlantChemElmnts,NumOfPlantLitrCmplxs))
  allocate(this%FWODRE(NumOfPlantChemElmnts,NumOfPlantLitrCmplxs))
  allocate(this%FWODBE(NumOfPlantChemElmnts,NumOfPlantLitrCmplxs))

  allocate(this%PetioleBiomGrowthYield(JP1))
  allocate(this%HuskBiomGrowthYield(JP1))
  allocate(this%StalkBiomGrowthYield(JP1))
  allocate(this%ReserveBiomGrowthYield(JP1))
  allocate(this%EarBiomGrowthYield(JP1))
  allocate(this%GrainBiomGrowthYield(JP1))
  allocate(this%rNCHusk_pft(JP1))
  allocate(this%rNCEar_pft(JP1))
  allocate(this%LeafBiomGrowthYield(JP1))

  end subroutine plt_allom_init
!----------------------------------------------------------------------

  subroutine plt_allom_destroy(this)
  implicit none

  class(plant_allometry_type) :: this

!  if(allocated(FNOD))deallocate(FNOD)
!  if(allocated(GrainSeedBiomCMean_brch))deallocate(GrainSeedBiomCMean_brch)
!  if(allocated(NoduGrowthYield_pft))deallocate(NoduGrowthYield_pft)
!  if(allocated(RootBiomGrowthYield))deallocate(RootBiomGrowthYield)
!  if(allocated(RootrPC_pft))deallocate(RootrPC_pft)
!  if(allocated(rCNNonstructRemob_pft))deallocate(rCNNonstructRemob_pft)
!  if(allocated(rCPNonstructRemob_pft))deallocate(rCPNonstructRemob_pft)
!  if(allocated(CPRTS))deallocate(CPRTS)
!  if(allocated(CNRTS))deallocate(CNRTS)
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
!  if(allocated(FVRN))deallocate(FVRN)
!  if(allocated(FWODLE))deallocate(FWODLE)
!  if(allocated(FWOODE))deallocate(FWOODE)
!  if(allocated(FWODRE))deallocate(FWODRE)
!  if(allocated(FWODBE))deallocate(FWODBE)
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

  allocate(this%ZEROL(JP1))
  allocate(this%ZEROP(JP1))
  allocate(this%StandingDeadChemElmnt_col(NumOfPlantChemElmnts))
  allocate(this%RootNodueChemElmnt_pvr(NumOfPlantChemElmnts,JZ1,JP1))
  allocate(this%CanopyLeafCpft_lyr(NumOfCanopyLayers1,JP1))
  allocate(this%RootNoduleNonstructElmnt_vr(NumOfPlantChemElmnts,JZ1,JP1))
  allocate(this%StandingDeadKCompChemElmnts_pft(NumOfPlantChemElmnts,jsken,JP1))
  allocate(this%Root2ndStructChemElmnt_pvr(NumOfPlantChemElmnts,jroots,JZ1,NumOfCanopyLayers1,JP1))
  allocate(this%Root1stStructChemElmnt_pvr(NumOfPlantChemElmnts,jroots,JZ1,NumOfCanopyLayers1,JP1))
  allocate(this%CanopyNonstructElementConc_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%CanopyNonstructElements_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%NoduleNonstructElmnt_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%NoduleNonstructCconc_pft(JP1))
  allocate(this%RootProteinConc_pftvr(jroots,JZ1,JP1))
  allocate(this%RootProteinC_pvr(jroots,JZ1,JP1))
  allocate(this%RootStructBiomC_vr(jroots,JZ1,JP1))
  allocate(this% PopuPlantRootC_vr(jroots,JZ1,JP1))
  allocate(this%RootMycoNonstructElmnt_vr(NumOfPlantChemElmnts,jroots,JZ1,JP1))
  allocate(this%RootNonstructElementConcpft_vr(NumOfPlantChemElmnts,jroots,JZ1,JP1))
  allocate(this%NonstructElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%CanopyStalkC_pft(JP1))
  allocate(this%LeafProteinCNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%PetioleProteinCNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%InternodeChemElmnt_brch(NumOfPlantChemElmnts,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%LeafChemElmntNode_brch(NumOfPlantChemElmnts,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%PetioleElmntNode_brch(NumOfPlantChemElmnts,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%LeafChemElmntByLayer_pft(NumOfPlantChemElmnts,NumOfCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%StalkBiomassC_brch(MaxNumBranches,JP1))
  allocate(this%NoduleNonstructElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%LeafPetioNonstructElmntConc_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%RootStructChemElmnt_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%WGLFT(NumOfCanopyLayers1))
  allocate(this%RootChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%SeedCPlanted_pft(JP1))
  allocate(this%NonstructalChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%CanopyLeafShethC_pft(JP1))
  allocate(this%StandingDeadChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%RootBiomCPerPlant_pft(JP1))
  allocate(this%PetioleChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%StalkChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%ReserveChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%GrainChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%HuskChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%EarChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%NoduleChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%LeafChemElmntRemob_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%PetioleChemElmntRemob_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%ShootChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%WTSHTA(JP1))
  allocate(this%LeafChemElmnts_pft(NumOfPlantChemElmnts,JP1))
  allocate(this%StandingDeadInitC_pft(JP1))
  allocate(this%LeafPetioleBiomassC_brch(MaxNumBranches,JP1))
  allocate(this%ReserveChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%LeafChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%CanopyNoduleChemElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%PetioleChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%EarChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%HuskChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%GrainChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%StalkChemElmnts_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%ShootChemElmnt_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%BranchStalkChemElmnts_pft_pft(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%Root1stChemElmnt(NumOfPlantChemElmnts,jroots,JRS,JP1))

  end subroutine plt_biom_init
!----------------------------------------------------------------------

  subroutine plt_biom_destroy(this)

  implicit none
  class(plant_biom_type) :: this

!  if(allocated(ZEROL))deallocate(ZEROL)
!  if(allocated(ZEROP))deallocate(ZEROP)
!  if(allocated(RootNoduleNonstructElmnt_vr))deallocate(RootNoduleNonstructElmnt_vr)
!  if(allocated(RootNodueChemElmnt_pvr))deallocate(RootNodueChemElmnt_pvr)
!  if(allocated(CanopyLeafCpft_lyr))deallocate(CanopyLeafCpft_lyr)
!  call destroy(Root1stChemElmnt)
!  if(allocated(Root1stStructChemElmnt_pvr))deallocate(Root1stStructChemElmnt_pvr)
!  if(allocated(StandingDeadKCompChemElmnts_pft))deallocate(StandingDeadKCompChemElmnts_pft)
!  if(allocated(Root2ndStructChemElmnt_pvr))deallocate(Root2ndStructChemElmnt_pvr)
!  if(allocated(CanopyNonstructElementConc_pft))deallocate(CanopyNonstructElementConc_pft)
!  if(allocated(CanopyNonstructElements_pft))deallocate(CanopyNonstructElements_pft)
!  if(allocated(NoduleNonstructElmnt_pft))deallocate(NoduleNonstructElmnt_pft)
!  if(allocated(NoduleNonstructCconc_pft))deallocate(NoduleNonstructCconc_pft)
!  if(allocated(RootProteinConc_pftvr))deallocate(RootProteinConc_pftvr)
!  if(allocated(RootProteinC_pvr))deallocate(RootProteinC_pvr)
!  if(allocated(RootStructBiomC_vr))deallocate(RootStructBiomC_vr)
!  if(allocated( PopuPlantRootC_vr))deallocate( PopuPlantRootC_vr)
!  if(allocated(RootMycoNonstructElmnt_vr))deallocate(RootMycoNonstructElmnt_vr)
!  if(allocated(WVSTK))deallocate(WVSTK)
!  if(allocated(LeafChemElmntNode_brch))deallocate(LeafChemElmntNode_brch)
!  if(allocated(LeafProteinCNode_brch))deallocate(LeafProteinCNode_brch)
!  if(allocated(PetioleProteinCNode_brch))deallocate(PetioleProteinCNode_brch)
!  if(allocated(LeafChemElmntByLayer_pft))deallocate(LeafChemElmntByLayer_pft)
!  if(allocated(InternodeChemElmnt_brch))deallocate(InternodeChemElmnt_brch)
!  if(allocated(PetioleElmntNode_brch))deallocate(PetioleElmntNode_brch)
!  if(allocated(StalkBiomassC_brch))deallocate(StalkBiomassC_brch)
!  if(allocated(PPOOL))deallocate(PPOOL)
!  if(allocated(NoduleNonstructElmnt_brch))deallocate(NoduleNonstructElmnt_brch)
!  if(allocated(LeafPetioNonstructElmntConc_brch))deallocate(LeafPetioNonstructElmntConc_brch)
!  if(allocated(RootStructChemElmnt_pft))deallocate(RootStructChemElmnt_pft)
!  if(allocated(WGLFT))deallocate(WGLFT)
!  if(allocated(NonstructElmnt_brch))deallocate(NonstructElmnt_brch)
!  if(allocated(RootBiomCPerPlant_pft))deallocate(RootBiomCPerPlant_pft)
!  if(allocated(WTLF))deallocate(WTLF)
!  if(allocated(WTSHE))deallocate(WTSHE)
!  if(allocated(WTRSV))deallocate(WTRSV)
!  if(allocated(WTSTK))deallocate(WTSTK)
!  if(allocated(WTLFP))deallocate(WTLFP)
!  if(allocated(WTGR))deallocate(WTGR)
!  if(allocated(WTLFN))deallocate(WTLFN)
!  if(allocated(WTEAR))deallocate(WTEAR)
!  if(allocated(HuskChemElmnts_pft))deallocate(HuskChemElmnts_pft)
!  if(allocated(LeafChemElmntRemob_brch))deallocate(LeafChemElmntRemob_brch)
!  if(allocated(PetioleChemElmntRemob_brch))deallocate(PetioleChemElmntRemob_brch)
!  if(allocated(LeafPetioleBiomassC_brch))deallocate(LeafPetioleBiomassC_brch)
!  if(allocated(ReserveChemElmnts_brch))deallocate(ReserveChemElmnts_brch)
!  if(allocated(LeafChemElmnts_brch))deallocate(LeafChemElmnts_brch)
!  if(allocated(CanopyNoduleChemElmnt_brch))deallocate(CanopyNoduleChemElmnt_brch)
!  if(allocated(PetioleChemElmnts_brch))deallocate(PetioleChemElmnts_brch)
!  if(allocated(EarChemElmnts_brch))deallocate(EarChemElmnts_brch)
!  if(allocated(HuskChemElmnts_brch))deallocate(HuskChemElmnts_brch)
!  if(allocated(GrainChemElmnts_brch))deallocate(GrainChemElmnts_brch)
!  if(allocated(StalkChemElmnts_brch))deallocate(StalkChemElmnts_brch)
!  if(allocated(ShootChemElmnt_brch))deallocate(ShootChemElmnt_brch)
!  if(allocated(BranchStalkChemElmnts_pft_pft))deallocate(BranchStalkChemElmnts_pft_pft)
!  if(allocated(WTSTDI))deallocate(WTSTDI)
!  if(allocated(SeedCPlanted_pft))deallocate(SeedCPlanted_pft)
!  if(allocated(WTLS))deallocate(WTLS)
!  if(allocated(ShootChemElmnts_pft))deallocate(ShootChemElmnts_pft)
!  if(allocated(WTSHTA))deallocate(WTSHTA)
!  if(allocated(StandingDeadChemElmnts_pft)deallocate(StandingDeadChemElmnts_pft)
!  if(allocated(NoduleChemElmnts_pft))deallocate(NoduleChemElmnts_pft)
  end subroutine plt_biom_destroy


!----------------------------------------------------------------------

  subroutine  plt_soilchem_init(this)

  implicit none

  class(plant_soilchem_type) :: this

  allocate(this%FOSRH(1:jcplx,0:JZ1))
  allocate(this%CFOPE(NumOfPlantChemElmnts,0:NumLitterGroups,jsken,JP1))
  allocate(this%TFND(0:JZ1))
  allocate(this%THETPM(60,0:JZ1))
  allocate(this%DiffusivitySolutEff(60,0:JZ1))
  allocate(this%VLSoilMicP(0:JZ1))
  allocate(this%VLiceMicP(0:JZ1))
  allocate(this%VLWatMicP(0:JZ1))
  allocate(this%VLMicP(0:JZ1))

  allocate(this%trcs_VLN(ids_nuts_beg:ids_nuts_end,0:JZ1))
  allocate(this%DOM(idom_beg:idom_end,1:jcplx,0:JZ1))

  allocate(this%trc_solml(ids_beg:ids_end,0:JZ1))

  allocate(this%trc_gasml(idg_beg:idg_end,0:JZ1))
  allocate(this%trc_gascl(idg_beg:idg_end,0:JZ1))
  allocate(this%CORGC(0:JZ1))

  allocate(this%trc_solcl(ids_beg:ids_end,0:jZ1))

  allocate(this%VLSoilPoreMicP(0:JZ1))
  allocate(this%THETW(0:JZ1))
  allocate(this%THETY(0:JZ1))

  allocate(this%GasSolbility(idg_beg:idg_end,0:JZ1))
  allocate(this%GasDifc(idg_beg:idg_end,0:JZ1))
  allocate(this%SolDifc(ids_beg:ids_end,0:JZ1))
  allocate(this%SoilResit4RootPentration(JZ1))
  allocate(this%SoiBulkDensity(0:JZ1))
  allocate(this%HydroCondMicP4RootUptake(JZ1))

  end subroutine plt_soilchem_init
!----------------------------------------------------------------------

  subroutine plt_soilchem_destroy(this)
  implicit none
  class(plant_soilchem_type) :: this

!  if(allocated(FOSRH))deallocate(FOSRH)

!  if(allocated(CFOPE))deallocate(CFOPE)
!  if(allocated(TFND))deallocate(TFND)
!  if(allocated(THETPM))deallocate(THETPM)
!  if(allocated(DiffusivitySolutEff))deallocate(DiffusivitySolutEff)
!  if(allocated(ZVSGL))deallocate(ZVSGL)
!  if(allocated(SOXYL))deallocate(SOXYL)

!  if(allocated(HydroCondMicP4RootUptake))deallocate(HydroCondMicP4RootUptake)
!  if(allocated(CGSGL))deallocate(CGSGL)
!  if(allocated(CHSGL))deallocate(CHSGL)
!  if(allocated(HGSGL))deallocate(HGSGL)
!  if(allocated(OGSGL))deallocate(OGSGL)
!  if(allocated(SoilResit4RootPentration))deallocate(SoilResit4RootPentration)

!   call destroy(this%GasSolbility)
!  if(allocated(THETW))deallocate(THETW)
!  if(allocated(THETY))deallocate(THETY)
!  if(allocated(VLSoilPoreMicP))deallocate(VLSoilPoreMicP)
!  if(allocated(CCH4G))deallocate(CCH4G)
!  if(allocated(CZ2OG))deallocate(CZ2OG)
!  if(allocated(CNH3G))deallocate(CNH3G)
!  if(allocated(CH2GG))deallocate(CH2GG)
!  if(allocated(CORGC))deallocate(CORGC)
!  if(allocated(H2PO4))deallocate(H2PO4)
!  if(allocated(HLSGL))deallocate(HLSGL)

!   call destroy(this%trcs_VLN)
!  if(allocated(OLSGL))deallocate(OLSGL)
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
!  if(allocated(SoiBulkDensity))deallocate(SoiBulkDensity)
!  if(allocated(CO2G))deallocate(CO2G)
!  if(allocated(CLSGL))deallocate(CLSGL)
!  if(allocated(CQSGL))deallocate(CQSGL)

  end subroutine plt_soilchem_destroy
!----------------------------------------------------------------------


  subroutine InitPlantAPIData()

  implicit none

  JZ1    => pltpar%JZ1
  NumOfCanopyLayers1    => pltpar%NumOfCanopyLayers1
  JP1    => pltpar%JP1
  JLA1   => pltpar%JLA1
  JSA1   => pltpar%JSA1
  JLI1   => pltpar%JLI1
  MaxNodesPerBranch1 => pltpar%MaxNodesPerBranch1
  !the following variable should be consistent with the soil bgc model
  jcplx => pltpar%jcplx
  jsken  => pltpar%jsken
  NumLitterGroups=> pltpar%NumLitterGroups
  MaxNumBranches    => pltpar%MaxNumBranches
  JRS    => pltpar%JRS
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

  allocate(this%PAR(JLI1,JSA1,NumOfCanopyLayers1,JP1))
  allocate(this%PARDIF(JLI1,JSA1,NumOfCanopyLayers1,JP1))
  allocate(this%CanopySWAlbedo_pft(JP1))
  allocate(this%CanopyPARalbedo_pft(JP1))
  allocate(this%TAUS(NumOfCanopyLayers1+1))
  allocate(this%TAU0(NumOfCanopyLayers1+1))
  allocate(this%LWRadCanP(JP1))
  allocate(this%SWRadByCanP(JP1))
  allocate(this%OMEGX(JSA1,JLI1,JLA1))
  allocate(this%OMEGAG(JSA1))
  allocate(this%OMEGA(JSA1,JLI1,JLA1))
  allocate(this%ZSIN(JLI1))
  allocate(this%ZCOS(JLI1))
  allocate(this%IALBY(JSA1,JLI1,JLA1))
  allocate(this%RadNet2CanP(JP1))
  allocate(this%CanopySWabsorpty_pft(JP1))
  allocate(this%CanopyPARabsorpty_pft(JP1))
  allocate(this%TAUP(JP1))
  allocate(this%TAUR(JP1))
  allocate(this%PARbyCanopy_pft(JP1))
  allocate(this%FracPARbyCanopy_pft(JP1))
  end subroutine plt_rad_init
!------------------------------------------------------------------------
  subroutine plt_rad_destroy(this)
! DESCRIPTION
! deallocate memory for plant_radiation_type
  implicit none
  class(plant_radiation_type) :: this

!  if(allocated(PAR))deallocate(PAR)
!  if(allocated(PARDIF))deallocate(PARDIF)
!  if(associated(this%CanopySWAlbedo_pft))deallocate(this%CanopySWAlbedo_pft)
!  if(associated(this%CanopyPARalbedo_pft))deallocate(this%CanopyPARalbedo_pft)
!  if(associated(this%TAUS))deallocate(this%TAUS)
!  if(associated(this%TAU0))deallocate(this%TAU0)
!  if(associated(this%LWRadCanP))deallocate(this%LWRadCanP)
!  if(associated(this%SWRadByCanP))deallocate(this%SWRadByCanP)
!  if(associated(this%OMEGX))deallocate(this%OMEGX)
!  if(associated(this%OMEGAG))deallocate(this%OMEGAG)
!  if(associated(this%ZSIN))deallocate(this%ZSIN)
!  if(allocated(ZCOS))deallocate(ZCOS)
!  if(associated(this%IALBY))deallocate(this%IALBY)
!  if(associated(this%OMEGA))deallocate(this%OMEGA)
!  if(associated(this%RadNet2CanP))deallocate(this%RadNet2CanP)
!  if(allocated(CanopySWabsorpty_pft))deallocate(CanopySWabsorpty_pft)
!  if(allocated(CanopyPARabsorpty_pft))deallocate(CanopyPARabsorpty_pft)
!  if(allocated(TAUR))deallocate(TAUR)
!  if(allocated(TAUP))deallocate(TAUP)
!  if(allocated(FracPARbyCanopy_pft))deallocate(FracPARbyCanopy_pft)
!  if(allocated(RADP))deallocate(RADP)

  end subroutine plt_rad_destroy

!------------------------------------------------------------------------

  subroutine plt_photo_init(this)
  class(plant_photosyns_type) :: this

  allocate(this%RSMX(JP1))
  allocate(this%MinCanPStomaResistH2O(JP1))
  allocate(this%SO2(JP1))
  allocate(this%CanPbndlResist(JP1))
  allocate(this%CanPStomaResistH2O(JP1))
  allocate(this%DiffCO2Atmos2Intracel_pft(JP1))
  allocate(this%AirConc_pft(JP1))
  allocate(this%CO2CuticleResist_pft(JP1))
  allocate(this%LeafAUnshaded_seclyrnodbrpft(JLI1,NumOfCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CPOOL3(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CPOOL4(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CMassCO2BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CO2CompenPoint_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%RubiscoCarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%C4CarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%LigthSatCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%LigthSatC4CarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CMassHCO3BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1))

  allocate(this%Vmax4RubiscoCarboxy_pft(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CO2lmtRubiscoCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%Vmax4PEPCarboxy_pft(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CO2lmtPEPCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%iPlantPhotosynthesisType(JP1))
  allocate(this%Km4PEPCarboxy_pft(JP1))
  allocate(this%O2L(JP1))
  allocate(this%CO2Solubility_pft(JP1))
  allocate(this%CanopyGasCO2_pft(JP1))
  allocate(this%CHILL(JP1))
  allocate(this%ETMX(JP1))
  allocate(this%CHL(JP1))
  allocate(this%PEPC(JP1))
  allocate(this%CHL4(JP1))
  allocate(this%RUBP(JP1))
  allocate(this%VCMX4(JP1))
  allocate(this%VOMX(JP1))
  allocate(this%VCMX(JP1))
  allocate(this%XKCO2(JP1))
  allocate(this%XKO2(JP1))
  allocate(this%RubiscoActivity_brpft(MaxNumBranches,JP1))
  allocate(this%C4PhotosynDowreg_brch(MaxNumBranches,JP1))
  allocate(this%CO2L(JP1))
  allocate(this%XKCO2L(JP1))
  allocate(this%Km4RubiscoCarboxy_pft(JP1))
  allocate(this%LeafIntracellularCO2_pft(JP1))
  allocate(this%O2I(JP1))
  allocate(this%RCS(JP1))
  allocate(this%CanPCi2CaRatio(JP1))
  allocate(this%MaxCanPStomaResistH2O(JP1))

  end subroutine plt_photo_init
!------------------------------------------------------------------------

  subroutine plt_photo_destroy(this)
  class(plant_photosyns_type) :: this

!  if(allocated(RSMX))deallocate(RSMX)
!  if(allocated(MinCanPStomaResistH2O))deallocate(MinCanPStomaResistH2O)
!  if(allocated(SO2))deallocate(SO2)
!  if(allocated(CanPStomaResistH2O))deallocate(CanPStomaResistH2O)
!  if(allocated(CanPbndlResist))deallocate(CanPbndlResist)
!  if(allocated(DiffCO2Atmos2Intracel_pft))deallocate(DiffCO2Atmos2Intracel_pft)
!  if(allocated(AirConc_pft))deallocate(AirConc_pft)
!  if(allocated(CO2CuticleResist_pft))deallocate(CO2CuticleResist_pft)
!  if(allocated(LeafAUnshaded_seclyrnodbrpft))deallocate(LeafAUnshaded_seclyrnodbrpft)
!  if(allocated(CPOOL3))deallocate(CPOOL3)
!  if(allocated(CPOOL4))deallocate(CPOOL4)
!  if(allocated(CMassCO2BundleSheath_node))deallocate(CMassCO2BundleSheath_node)
!  if(allocated(CO2CompenPoint_node))deallocate(CO2CompenPoint_node)
!  if(allocated(RubiscoCarboxyEff_node))deallocate(RubiscoCarboxyEff_node)
!  if(allocated(C4CarboxyEff_node))deallocate(C4CarboxyEff_node)
!  if(allocated(LigthSatCarboxyRate_node))deallocate(LigthSatCarboxyRate_node)
!  if(allocated(LigthSatC4CarboxyRate_node))deallocate(LigthSatC4CarboxyRate_node)
!  if(allocated(NutrientCtrlonC4Carboxy_node))deallocate(NutrientCtrlonC4Carboxy_node)
!  if(allocated(CMassHCO3BundleSheath_node))deallocate(CMassHCO3BundleSheath_node)

!  if(allocated(Vmax4RubiscoCarboxy_pft))deallocate(Vmax4RubiscoCarboxy_pft)
!  if(allocated(CO2lmtRubiscoCarboxyRate_node))deallocate(CO2lmtRubiscoCarboxyRate_node)
!  if(allocated(Vmax4PEPCarboxy_pft))deallocate(Vmax4PEPCarboxy_pft)
!  if(allocated(CO2lmtPEPCarboxyRate_node))deallocate(CO2lmtPEPCarboxyRate_node)
!  if(allocated(iPlantPhotosynthesisType))deallocate(iPlantPhotosynthesisType)
!  if(allocated(Km4PEPCarboxy_pft))deallocate(Km4PEPCarboxy_pft)
!  if(allocated(O2L))deallocate(O2L)
!  if(allocated(CO2Solubility_pft))deallocate(CO2Solubility_pft)
!  if(allocated(CanopyGasCO2_pft))deallocate(CanopyGasCO2_pft)
!  if(allocated(CHILL))deallocate(CHILL)
!  if(allocated(ETMX))deallocate(ETMX)
!  if(allocated(CHL))deallocate(CHL)
!  if(allocated(PEPC))deallocate(PEPC)
!  if(allocated(CHL4))deallocate(CHL4)
!  if(allocated(RUBP))deallocate(RUBP)
!  if(allocated(VCMX4))deallocate(VCMX4)
!  if(allocated(VOMX))deallocate(VOMX)
!  if(allocated(VCMX))deallocate(VCMX)
!  if(allocated(XKCO2))deallocate(XKCO2)
!  if(allocated(XKO2))deallocate(XKO2)
!  if(allocated(RubiscoActivity_brpft))deallocate(RubiscoActivity_brpft)
!  if(allocated(C4PhotosynDowreg_brch))deallocate(C4PhotosynDowreg_brch)
!  if(allocated(CO2L))deallocate(CO2L)
!  if(allocated(XKCO2L))deallocate(XKCO2L)
!  if(allocated(Km4RubiscoCarboxy_pft))deallocate(Km4RubiscoCarboxy_pft)
!  if(allocated(LeafIntracellularCO2_pft))deallocate(LeafIntracellularCO2_pft)
!  if(allocated(O2I))deallocate(O2I)
!  if(allocated(RCS))deallocate(RCS)
!  if(allocated(CanPCi2CaRatio))deallocate(CanPCi2CaRatio)
!  if(allocated(MaxCanPStomaResistH2O)) deallocate(MaxCanPStomaResistH2O)
  end subroutine plt_photo_destroy

!------------------------------------------------------------------------
  subroutine plt_pheno_init(this)
  implicit none
  class(plant_pheno_type) :: this


  allocate(this%TCelciusChill4Seed(JP1))
  allocate(this%OFFST(JP1))
  allocate(this%MatureGroup_pft(JP1))
  allocate(this%PlantO2Stress(JP1))
  allocate(this%MinNonstructalC4InitBranch(JP1))
  allocate(this%MinNonstructuralC4InitRoot(JP1))
  allocate(this%fTgrowCanP(JP1))
  allocate(this%TCelciusChill4Leaf(JP1))
  allocate(this%TCG(JP1))
  allocate(this%TKG(JP1))
  allocate(this%TCelcius4LeafOffHarden(JP1))
  allocate(this%HoursCanopyPSITooLow(JP1))
  allocate(this%LeafElmntRemobFlx_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))
  allocate(this%PetioleChemElmntRemobFlx_brch(NumOfPlantChemElmnts,MaxNumBranches,JP1))

  allocate(this%fTgrowRootP(JZ1,JP1))
  allocate(this%GrainFillRateat25C_pft(JP1))
  allocate(this%ShutRutNonstructElmntConducts(JP1))
  allocate(this%SSTX(JP1))
  allocate(this%HTC(JP1))
  allocate(this%iPlantInitThermoAdaptZone(JP1))
  allocate(this%iPlantThermoAdaptZone(JP1))
  allocate(this%IsPlantActive(JP1))
  allocate(this%iPlantState(JP1))
  allocate(this%RSETE(NumOfPlantChemElmnts,JP1))
  allocate(this%MatureGroup_brch(MaxNumBranches,JP1))
  allocate(this%iPlantShootState(JP1))
  allocate(this%HourCounter4LeafOut_brch(MaxNumBranches,JP1))
  allocate(this%HoursDoingRemob_brch(MaxNumBranches,JP1))
  allocate(this%HourlyNodeNumNormByMatgrp_brch(MaxNumBranches,JP1))
  allocate(this%HourReprodNodeNumNormByMatrgrp_brch(MaxNumBranches,JP1))
  allocate(this%NodeNumNormByMatgrp_brch(MaxNumBranches,JP1))
  allocate(this%ReprodNodeNumNormByMatrgrp_brch(MaxNumBranches,JP1))
  allocate(this%HourWithNoGrainFill_brch(MaxNumBranches,JP1))
  allocate(this%iPlantPhenologyType(JP1))
  allocate(this%iPlantPhenologyPattern(JP1))
  allocate(this%iPlantTurnoverPattern(JP1))
  allocate(this%iPlantRootState(JP1))
  allocate(this%iPlantDevelopPattern(JP1))
  allocate(this%iPlantPhotoperiodType(JP1))
  allocate(this%doInitPlant(JP1))
  allocate(this%iPlantMorphologyType(JP1))
  allocate(this%KLeafNodeNumber(MaxNumBranches,JP1))
  allocate(this%KLeafNumLowestGrowing_pft(MaxNumBranches,JP1))
  allocate(this%iPlantCalendar(NumGrowthStages,MaxNumBranches,JP1))
  allocate(this%TotalNodeNumNormByMatgrp_brch(MaxNumBranches,JP1))
  allocate(this%TotalReprodNodeNumNormByMatrgrp_brch(MaxNumBranches,JP1))
  allocate(this%LeafNumberAtFloralInit_brch(MaxNumBranches,JP1))
  allocate(this%RefLeafAppearRate(JP1))
  allocate(this%RefNodeInitRate(JP1))
  allocate(this%CriticalPhotoPeriod_pft(JP1))
  allocate(this%PhotoPeriodSens_pft(JP1))
  allocate(this%iPlantBranchState(MaxNumBranches,JP1))
  allocate(this%doRemobilization_brch(MaxNumBranches,JP1))
  allocate(this%doPlantLeaveOff_brch(MaxNumBranches,JP1))
  allocate(this%doPlantLeafOut_brch(MaxNumBranches,JP1))
  allocate(this%doInitLeafOut_brch(MaxNumBranches,JP1))
  allocate(this%doSenescence_brch(MaxNumBranches,JP1))
  allocate(this%Prep4Literfall_brch(MaxNumBranches,JP1))
  allocate(this%Hours4LiterfalAftMature_brch(MaxNumBranches,JP1))
  allocate(this%Hours4LenthenPhotoPeriod(MaxNumBranches,JP1))
  allocate(this%Hours4ShortenPhotoPeriod(MaxNumBranches,JP1))
  allocate(this%Hours4Leafout(MaxNumBranches,JP1))
  allocate(this%HourThreshold4LeafOut_brch(NumOfCanopyLayers1,JP1))
  allocate(this%Hours4LeafOff(MaxNumBranches,JP1))
  allocate(this%HourThreshold4LeafOff(NumOfCanopyLayers1,JP1))


  end subroutine plt_pheno_init
!------------------------------------------------------------------------

  subroutine plt_pheno_destroy(this)
  implicit none
  class(plant_pheno_type) :: this

!  if(allocated(TCelciusChill4Seed))deallocate(TCelciusChill4Seed)
!  if(allocated(OFFST))deallocate(OFFST)
!  if(allocated(MatureGroup_pft))deallocate(MatureGroup_pft)
!  if(allocated(PlantO2Stress))deallocate(PlantO2Stress)
!  if(allocated(PB))deallocate(PB)
!  if(allocated(PR))deallocate(PR)
!  if(allocated(fTgrowCanP))deallocate(fTgrowCanP)
!  if(allocated(TCelciusChill4Leaf))deallocate(TCelciusChill4Leaf)
!  if(allocated(TCG))deallocate(TCG)
!  if(allocated(TKG))deallocate(TKG)
!  if(allocated(TCelcius4LeafOffHarden))deallocate(TCelcius4LeafOffHarden)
!  if(allocated(HoursCanopyPSITooLow))deallocate(HoursCanopyPSITooLow)
!  if(allocated(LeafElmntRemobFlx_brch))deallocate(LeafElmntRemobFlx_brch)
!  if(allocated(PetioleChemElmntRemobFlx_brch))deallocate(PetioleChemElmntRemobFlx_brch)
!  if(allocated(fTgrowRootP))deallocate(fTgrowRootP)
!  if(allocated(GrainFillRateat25C_pft))deallocate(GrainFillRateat25C_pft)
!  if(allocated(RootAreaPerPlant_vr))deallocate(RootAreaPerPlant_vr)
!  if(allocated(ShutRutNonstructElmntConducts))deallocate(ShutRutNonstructElmntConducts)
!  if(allocated(SSTX)) deallocate(SSTX)
!  if(allocated(HTC))deallocate(HTC)
!  if(allocated(iPlantInitThermoAdaptZone))deallocate(iPlantInitThermoAdaptZone)
!  if(allocated(iPlantThermoAdaptZone))deallocate(iPlantThermoAdaptZone)
!  if(allocated(IsPlantActive))deallocate(IsPlantActive)
!  if(allocated(IDTH))deallocate(IDTH)
!  if(allocated(RSETE))deallocate(RSETE)
!  if(allocated(GROUP))deallocate(GROUP)
!  if(allocated(iPlantShootState))deallocate(iPlantShootState)
!  if(allocated(HourCounter4LeafOut_brch))deallocate(HourCounter4LeafOut_brch)
!  if(allocated(HoursDoingRemob_brch))deallocate(HoursDoingRemob_brch)
!  if(allocated(HourlyNodeNumNormByMatgrp_brch))deallocate(HourlyNodeNumNormByMatgrp_brch)
!  if(allocated(HourReprodNodeNumNormByMatrgrp_brch))deallocate(HourReprodNodeNumNormByMatrgrp_brch)
!  if(allocated(NodeNumNormByMatgrp_brch))deallocate(NodeNumNormByMatgrp_brch)
!  if(allocated(ReprodNodeNumNormByMatrgrp_brch))deallocate(ReprodNodeNumNormByMatrgrp_brch)
!  if(allocated(HourWithNoGrainFill_brch))deallocate(HourWithNoGrainFill_brch)
!  if(allocated(iPlantPhenologyType))deallocate(iPlantPhenologyType)
!  if(allocated(iPlantTurnoverPattern))deallocate(iPlantTurnoverPattern)
!  if(allocated(iPlantPhenologyPattern))deallocate(iPlantPhenologyPattern)
!  if(allocated(IDTHR))deallocate(IDTHR)
!  if(allocated(iPlantDevelopPattern))deallocate(iPlantDevelopPattern)
!  if(allocated(iPlantPhotoperiodType))deallocate(iPlantPhotoperiodType)
!  if(allocated(iPlantMorphologyType))deallocate(iPlantMorphologyType)
!  if(allocated(doInitPlant))deallocate(doInitPlant)
!  if(allocated(KLeafNumLowestGrowing_pft))deallocate(KLeafNumLowestGrowing_pft)
!  if(allocated(TotalNodeNumNormByMatgrp_brch))deallocate(TotalNodeNumNormByMatgrp_brch)
!  if(allocated(TotalReprodNodeNumNormByMatrgrp_brch))deallocate(TotalReprodNodeNumNormByMatrgrp_brch)
!  if(allocated(LeafNumberAtFloralInit_brch))deallocate(LeafNumberAtFloralInit_brch)
!  if(allocated(KLeafNodeNumber))deallocate(KLeafNodeNumber)
!  if(allocated(CriticalPhotoPeriod_pft))deallocate(CriticalPhotoPeriod_pft)
!  if(allocated(XRLA))deallocate(XRLA)
!  if(allocated(XRNI))deallocate(XRNI)
!  if(allocated(PhotoPeriodSens_pft))deallocate(PhotoPeriodSens_pft)
!  if(allocated(iPlantBranchState))deallocate(iPlantBranchState)
!  if(allocated(doRemobilization_brch))deallocate(doRemobilization_brch)
!  if(allocated(doPlantLeaveOff_brch))deallocate(doPlantLeaveOff_brch)
!  if(allocated(doPlantLeafOut_brch))deallocate(doPlantLeafOut_brch)
!  if(allocated(doInitLeafOut_brch))deallocate(doInitLeafOut_brch)
!  if(allocated(doSenescence_brch))deallocate(doSenescence_brch)
!  if(allocated(Prep4Literfall_brch))deallocate(Prep4Literfall_brch)
!  if(allocated(Hours4LiterfalAftMature_brch))deallocate(Hours4LiterfalAftMature_brch)
!  if(allocated(Hours4LenthenPhotoPeriod))deallocate(Hours4LenthenPhotoPeriod)
!  if(allocated(Hours4ShortenPhotoPeriod))deallocate(Hours4ShortenPhotoPeriod)
!  if(allocated(Hours4Leafout))deallocate(Hours4Leafout)
!  if(allocated(VRNL))deallocate(VRNL)
!  if(allocated(Hours4LeafOff))deallocate(Hours4LeafOff)
!  if(allocated(VRNX))deallocate(VRNX)
!  if(allocated(IDAY))deallocate(IDAY)

  end subroutine plt_pheno_destroy
!------------------------------------------------------------------------

  subroutine plt_morph_init(this)
  implicit none
  class(plant_morph_type) :: this

  allocate(this%RootAreaPerPlant_vr(jroots,JZ1,JP1))
  allocate(this%RootLenthDensPerPopu_pvr(jroots,JZ1,JP1))
  allocate(this%RootVolume_vr(jroots,JZ1,JP1))
  allocate(this%RootVH2O_vr(jroots,JZ1,JP1))
  allocate(this%PrimRootXNumL(jroots,JZ1,JP1))
  allocate(this%SecndRootXNum_pvr(jroots,JZ1,JP1))
  allocate(this%SeedCMass(JP1))
  allocate(this%RTDNT(JZ1))
  allocate(this%RootBranchFreq_pft(JP1))
  allocate(this%ClumpFactort0(JP1))
  allocate(this%ClumpFactort(JP1))
  allocate(this%HypoctoylHeight(JP1))
  allocate(this%RootPoreTortu4Gas(jroots,JP1))
  allocate(this%WDLF(JP1))
  allocate(this%MaxSeedNumPerSite_pft(JP1))
  allocate(this%MaxPotentSeedNumber_pft(JP1))
  allocate(this%MaxSeedCMass(JP1))
  allocate(this%MaxPrimRootRadius1(jroots,JP1))
  allocate(this%MaxSecndRootRadius1(jroots,JP1))
  allocate(this%RRADP(jroots,JP1))
  allocate(this%PrimRootRadius(jroots,JZ1,JP1))
  allocate(this%SecndRootRadius(jroots,JZ1,JP1))
  allocate(this%MaxPrimRootRadius(jroots,JP1))
  allocate(this%MaxSecndRootRadius(jroots,JP1))

  allocate(this%PrimRootDepth(jroots,NumOfCanopyLayers1,JP1))
  allocate(this%RootLenPerPopu_pvr(jroots,JZ1,JP1))
  allocate(this%AveSecndRootLen(jroots,JZ1,JP1))
  allocate(this%PrimRootSpecLen(jroots,JP1))
  allocate(this%SecndRootSpecLen(jroots,JP1))
  allocate(this%PrimRootLen(jroots,JZ1,NumOfCanopyLayers1,JP1))
  allocate(this%SecndRootLen(jroots,JZ1,NumOfCanopyLayers1,JP1))
  allocate(this%SecndRootXNum_rpvr(jroots,JZ1,NumOfCanopyLayers1,JP1))
  allocate(this%iPlantNfixType(JP1))
  allocate(this%MY(JP1))
  allocate(this%CanPHeight4WatUptake(JP1))
  allocate(this%KLeafNumber_brch(MaxNumBranches,JP1))
  allocate(this%NumOfLeaves_brch(MaxNumBranches,JP1))
  allocate(this%NGTopRootLayer_pft(JP1))
  allocate(this%CanopyHeight(JP1))
  allocate(this%XTLI(JP1))
  allocate(this%CanopyHeightz(0:NumOfCanopyLayers1))
  allocate(this%CanopyStemA_pft(JP1))
  allocate(this%CanopyLeafA_pft(JP1))
  allocate(this%NumOfMainBranch_pft(JP1))
  allocate(this%NIXBotRootLayer_pft(JP1))
  allocate(this%NumRootAxes_pft(JP1))
  allocate(this%NumConCurrentGrowinNode(JP1))
  allocate(this%BranchNumber_pft(JP1))
  allocate(this%NumOfBranches_pft(JP1))
  allocate(this%NIXBotRootLayer_rpft(NumOfCanopyLayers1,JP1))
  allocate(this%ShootNodeNumber(MaxNumBranches,JP1))
  allocate(this%NodeNumberToInitFloral(MaxNumBranches,JP1))
  allocate(this%NodeNumberAtAnthesis(MaxNumBranches,JP1))
  allocate(this%SineBranchAngle_pft(JP1))
  allocate(this%LeafA_lyrnodbrchpft(JLI1,NumOfCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%KLEAFX(MaxNumBranches,JP1))
  allocate(this%BranchNumber_brch(MaxNumBranches,JP1))
  allocate(this%PotentialSeedSites_brch(MaxNumBranches,JP1))
  allocate(this%CanPBranchHeight(MaxNumBranches,JP1))
  allocate(this%LeafAreaDying_brch(MaxNumBranches,JP1))
  allocate(this%LeafAreaLive_brch(MaxNumBranches,JP1))
  allocate(this%SinePetioleAngle_pft(JP1))
  allocate(this%CLASS(JLI1,JP1))
  allocate(this%CanopyStemApft_lyr(NumOfCanopyLayers1,JP1))
  allocate(this%CanopyLeafApft_lyr(NumOfCanopyLayers1,JP1))
  allocate(this%LeafAreaNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CanopySeedNumber_pft(JP1))
  allocate(this%SeedinDepth(JP1))
  allocate(this%PlantinDepth(JP1))
  allocate(this%SeedLengthMean_pft(JP1))
  allocate(this%SeedVolumeMean_pft(JP1))
  allocate(this%SeedAreaMean_pft(JP1))
  allocate(this%CanopyStemA_lyr(NumOfCanopyLayers1))
  allocate(this%SSL1(JP1))
  allocate(this%SNL1(JP1))
  allocate(this%SLA1(JP1))
  allocate(this%CanopyLAgrid_lyr(NumOfCanopyLayers1))
  allocate(this%CanopyArea_pft(JP1))
  allocate(this%InternodeHeightDying_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%PetioleLengthNode_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%InternodeHeightLive_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%StemA_lyrnodbrchpft(JLI1,NumOfCanopyLayers1,MaxNumBranches,JP1))
  allocate(this%CanopyLeafAreaByLayer_pft(NumOfCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  allocate(this%CanopyBranchStemApft_lyr(NumOfCanopyLayers1,MaxNumBranches,JP1))
  allocate(this%NI(JP1))
  allocate(this%SeedNumberSet_brch(MaxNumBranches,JP1))
  allocate(this%ClumpFactor(JP1))
  allocate(this%RootVolPerMassC_pft(jroots,JP1))
  allocate(this%RootPorosity(jroots,JP1))
  allocate(this%SecndRootXSecArea(jroots,JP1))
  allocate(this%PrimRootXSecArea(jroots,JP1))
  allocate(this%RSRR(jroots,JP1))
  allocate(this%RSRA(jroots,JP1))
  allocate(this%iPlantGrainType(JP1))
  end subroutine plt_morph_init

!------------------------------------------------------------------------

  subroutine plt_morph_destroy(this)
  implicit none
  class(plant_morph_type) :: this

  end subroutine plt_morph_destroy
end module PlantAPIData
