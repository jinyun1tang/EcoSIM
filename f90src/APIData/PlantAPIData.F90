module PlantAPIData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ElmIDMod
  use abortutils, only : destroy
  use EcoSiMParDataMod, only : pltpar
  use TracerIDMod
implicit none
  save
  character(len=*),private, parameter :: mod_filename = __FILE__
  public

  integer, pointer :: jpstgs       !number of growth stages
  integer, pointer :: JRS          !maximum number of root layers
  integer, pointer :: JBR         !maximum number of plant branches
  integer, pointer :: JP1         !number of plants
  integer, pointer :: JSA1        !number of sectors for the sky azimuth  [0,2*pi]
  integer, pointer :: jcplx       !number of organo-microbial complexes
  integer, pointer :: JLA1        !number of sectors for the leaf azimuth, [0,pi]
  integer, pointer :: JC1         !number of canopy layers
  integer, pointer :: JZ1         !number of soil layers
  integer, pointer :: JLI1        !number of sectors for the leaf zenith [0,pi/2]
  integer, pointer :: JNODS1      !number of canopy nodes
  integer, pointer :: jsken       !number of kinetic components in litter,  PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  integer, pointer :: Jlitgrp     !number of litter groups nonstructural(0,*),
                          !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
  integer, pointer :: JPRT        !number of organs involved in partition
  integer, pointer :: n_pltlitrk
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
  real(r8) :: DYLX      !daylength of previous day, [h]
  real(r8) :: DYLN      !daylength, [h]
  real(r8) :: DYLM      !maximum daylength, [h]
  real(r8) :: ZS        !initial soil surface roughness height, [m]
  real(r8) :: OXYE      !atmospheric O2 concentration, [umol mol-1]
  real(r8) :: PPT       !total plant population, [d-2]
  real(r8) :: POROS1    !top layer soil porosity
  real(r8) :: UA        !wind speed, [m h-1]
  real(r8), pointer :: TBALE(:)  => null() !total plant element balance	g d-2
  real(r8) :: Z0        !wind speed measurement height, [m]
  real(r8) :: ZEROS2    !threshold zero
  real(r8) :: ZEROS     !threshold zero
  real(r8) :: ZERO      !threshold zero
  real(r8) :: ZERO2     !threshold zero
  integer :: IETYP      !Koppen climate zone
  integer :: IFLGT      !number of active PFT
  integer :: NL         !lowest soil layer number
  integer :: NP0        !intitial number of plant species
  integer :: NJ         !maximum root layer number
  integer :: NP         !number of plant species
  integer :: NU         !soil surface layer number
  integer :: LYRC       !number of days in current year
  integer :: IYRC       !current year

  real(r8) :: VOLWOU    !total subsurface water flux	m3 d-2
  character(len=16), pointer :: DATAP(:) => null()   !parameter file name
  CHARACTER(len=16), pointer :: DATA(:)  => null()   !pft file
  real(r8), pointer :: AREA3(:)   => null()    !soil cross section area (vertical plan defined by its normal direction)
  integer,  pointer :: IDATA(:)   => null()    !time keeper
  real(r8), pointer :: BALE(:,:)  => null()    !plant element balance, [g d-2]
  real(r8), pointer :: CDPTHZ(:)  => null()    !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8), pointer :: PPI(:)     => null()    !initial plant population, [m-2]
  real(r8), pointer :: PPZ(:)     => null()    !plant population at seeding, [m-2]
  real(r8), pointer :: PPX(:)     => null()    !plant population, [m-2]
  real(r8), pointer :: pftPlantPopulation(:)      => null()    !plant population, [d-2]
  real(r8), pointer :: DPTHZ(:)   => null()    !depth to middle of soil layer from  surface of grid cell [m]
  real(r8), pointer :: FMPR(:)    => null()    !micropore fraction
  real(r8), pointer :: DLYR3(:)   => null()    !vertical thickness of soil layer [m]
  real(r8), pointer :: VOLWM(:,:) => null()    !soil micropore water content, [m3 d-2]
  real(r8), pointer :: VOLPM(:,:) => null()    !soil air content, [m3 d-2]
  real(r8), pointer :: TORT(:,:)  => null()    !soil tortuosity, []
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
  real(r8), pointer :: FDBK(:,:) => null()   !branch down-regulation of CO2 fixation, [-]
  real(r8), pointer :: FDBKX(:,:)=> null()   !down-regulation of C4 photosynthesis, [-]
  real(r8), pointer :: CO2L(:)   => null()   !leaf aqueous CO2 concentration, [uM]
  real(r8), pointer :: O2I(:)    => null()   !leaf gaseous O2 concentration, [umol m-3]
  real(r8), pointer :: CO2I(:)   => null()   !leaf gaseous CO2 concentration, [umol m-3]
  real(r8), pointer :: XKCO2O(:) => null()   !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8), pointer :: XKCO2L(:) => null()   !leaf aqueous CO2 Km no O2, [uM]
  real(r8), pointer :: RCS(:)    => null()   !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8), pointer :: FCO2(:)   => null()   !Ci:Ca ratio, [-]
  real(r8), pointer :: RSMH(:)   => null()   !maximum stomatal resistance to vapor, [s h-1]
  real(r8), pointer :: CHILL(:)  => null()   !chilling effect on CO2 fixation, [-]
  real(r8), pointer :: SCO2(:)   => null()   !leaf CO2 solubility, [uM /umol mol-1]
  real(r8), pointer :: CO2Q(:)   => null()   !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8), pointer :: O2L(:)    => null()   !leaf aqueous O2 concentration, [uM]
  real(r8), pointer :: XKCO24(:) => null()   !Km for PEP carboxylase activity, [uM]
  integer,  pointer :: ICTYP(:)  => null()   !plant photosynthetic type (C3 or C4)
  real(r8), pointer :: RSMN(:)   => null()   !canopy minimum stomatal resistance, [s m-1]
  real(r8), pointer :: RSMX(:)   => null()   !maximum stomatal resistance to vapor, [s m-1]
  real(r8), pointer :: RC(:)     => null()   !canopy stomatal resistance, [h m-1]
  real(r8), pointer :: RA(:)     => null()   !canopy boundary layer resistance, [h m-1]
  real(r8), pointer :: SO2(:)    => null()   !leaf O2 solubility, [uM /umol mol-1]
  real(r8), pointer :: RCMX(:)   => null()   !maximum stomatal resistance to CO2, [s h-1]
  real(r8), pointer :: DCO2(:)   => null()   !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8), pointer :: SURFX(:,:,:,:,:)  => null()   !leaf irradiated surface area, [m2 d-2]
  real(r8), pointer :: CPOOL3(:,:,:)     => null()   !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8), pointer :: CPOOL4(:,:,:)     => null()   !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: CO2B(:,:,:)       => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: COMPL(:,:,:)      => null()   !CO2 compensation point, [uM]
  real(r8), pointer :: CBXN(:,:,:)       => null()   !carboxylation efficiency, [umol umol-1]
  real(r8), pointer :: CBXN4(:,:,:)      => null()   !C4 carboxylation efficiency, [umol umol-1]
  real(r8), pointer :: ETGRO(:,:,:)      => null()   !maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: ETGR4(:,:,:)      => null()   !maximum  light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: FDBK4(:,:,:)      => null()   !down-regulation of C4 photosynthesis, [-]
  real(r8), pointer :: HCOB(:,:,:)       => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: VCGRO(:,:,:)      => null()   !maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: VGRO(:,:,:)       => null()   !carboxylation rate, [umol m-2 s-1]
  real(r8), pointer :: VCGR4(:,:,:)      => null()   !maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), pointer :: VGRO4(:,:,:)      => null()   !C4 carboxylation rate, [umol m-2 s-1]
  real(r8), pointer :: FMOL(:)           => null()   !total gas concentration, [mol m-3]

  contains
    procedure, public :: Init    =>  plt_photo_init
    procedure, public :: Destroy => plt_photo_destroy
  end type plant_photosyns_type

  type, public ::  plant_radiation_type
  real(r8) :: TYSIN     !sine of sky angles
  real(r8) :: ALBS      !snowpack albedo
  real(r8) :: ALBX      !Surface albedo
  real(r8) :: RAP0      !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8) :: RADY      !diffuse shortwave radiation, [W m-2]
  real(r8) :: RADS      !direct shortwave radiation, [W m-2]
  real(r8) :: RAPY      !diffuse PAR, [umol m-2 s-1]
  real(r8) :: RAPS      !direct PAR, [umol m-2 s-1]
  real(r8) :: TRN       !ecosystem net radiation, [MJ d-2 h-1]
  real(r8) :: RAD0      !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8) :: FRADG     !fraction of radiation intercepted by ground surface, [-]
  real(r8) :: RADG      !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8) :: GSIN      !sine of slope, [-]
  real(r8) :: GAZI      !azimuth of slope, [-]
  real(r8) :: GCOS      !cosine of slope, [-]
  real(r8) :: THRMGX    !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8) :: THS       !sky longwave radiation , [MJ d-2 h-1]
  real(r8) :: SSINN     !sine of solar angle next hour, [-]
  real(r8) :: SSIN      !sine of solar angle, [-]
  real(r8), pointer :: ALBR(:)     => null() !canopy shortwave albedo , [-]
  real(r8), pointer :: ALBP(:)     => null() !canopy PAR albedo , [-]
  real(r8), pointer :: TAUS(:)     => null() !fraction of radiation intercepted by canopy layer, [-]
  real(r8), pointer :: TAU0(:)     => null() !fraction of radiation transmitted by canopy layer, [-]
  real(r8), pointer :: THRM1(:)    => null() !canopy longwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: RADC(:)     => null() !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: OMEGX(:,:,:)=> null() !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation
  real(r8), pointer :: OMEGAG(:)   => null() !sine of solar beam on leaf surface, [-]
  real(r8), pointer :: OMEGA(:,:,:)=> null() !sine of indirect sky radiation on leaf surface
  real(r8), pointer :: ZSIN(:)     => null() !sine of leaf angle
  real(r8), pointer :: ZCOS(:)     => null() !cosine of leaf angle
  integer,  pointer :: IALBY(:,:,:)=> null() !flag for calculating backscattering of radiation in canopy
  real(r8), pointer :: RAD1(:)     => null() !canopy net radiation , [MJ d-2 h-1]
  real(r8), pointer :: ABSR(:)     => null() !canopy shortwave absorptivity , [-]
  real(r8), pointer :: ABSP(:)     => null() !canopy PAR absorptivity
  real(r8), pointer :: TAUR(:)     => null() !canopy shortwave transmissivity , [-]
  real(r8), pointer :: TAUP(:)     => null() !canopy PAR transmissivity , [-]
  real(r8), pointer :: RADP(:)     => null() !canopy absorbed PAR , [umol m-2 s-1]
  real(r8), pointer :: FRADP(:)    => null() !fraction of incoming PAR absorbed by canopy, [-]
  real(r8), pointer :: PAR(:,:,:,:)=> null()     !direct incoming PAR, [umol m-2 s-1]
  real(r8), pointer :: PARDIF(:,:,:,:)=> null()  !diffuse incoming PAR, [umol m-2 s-1]
  contains
    procedure, public :: Init    => plt_rad_init
    procedure, public :: Destroy => plt_rad_destroy
  end type plant_radiation_type

  type, public :: plant_morph_type
  real(r8) :: ARLSS                     !stalk area of combined, each PFT canopy
  real(r8) :: ARLFC                     !total canopy leaf area, [m2 d-2]
  real(r8) :: ARSTC                     !total canopy stem area, [m2 d-2]
  real(r8) :: GridMaxCanopyHeight                        !canopy height , [m]
  real(r8), pointer :: DMVL(:,:)       => null() !root volume:mass ratio, [m3 g-1]
  real(r8), pointer :: PORT(:,:)       => null() !root porosity, [m3 m-3]
  real(r8), pointer :: RTAR2X(:,:)     => null() !root  cross-sectional area  secondary axes, [m2]
  real(r8), pointer :: RTAR1X(:,:)     => null() !root cross-sectional area primary axes, [m2]
  real(r8), pointer :: RRAD1X(:,:)     => null() !root diameter primary axes, [m]
  real(r8), pointer :: RRAD2X(:,:)     => null() !root diameter secondary axes, [m]
  real(r8), pointer :: GRDM(:)         => null() !grain size at seeding, [g]
  real(r8), pointer :: PORTX(:,:)      => null() !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8), pointer :: RTLG1(:,:,:,:)  => null() !root layer length primary axes, [m d-2]
  real(r8), pointer :: RTLG2(:,:,:,:)  => null() !root layer length secondary axes, [m d-2]
  real(r8), pointer :: RTLGP(:,:,:)    => null() !root layer length per plant, [m p-1]
  real(r8), pointer :: RTLGA(:,:,:)    => null() !root layer average length, [m]
  real(r8), pointer :: RTLG1X(:,:)     => null() !specific root length primary axes, [m g-1]
  real(r8), pointer :: RTLG2X(:,:)     => null() !specific root length secondary axes, [m g-1]
  real(r8), pointer :: RTN2(:,:,:,:)   => null() !root layer number secondary axes, [d-2]
  real(r8), pointer :: RTDP1(:,:,:)    => null() !root layer depth, [m]
  real(r8), pointer :: CFI(:)          => null() !initial clumping factor for self-shading in canopy layer, [-]
  real(r8), pointer :: CFX(:)          => null() !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8), pointer :: RTFQ(:)         => null() !root brancing frequency, [m-1]
  real(r8), pointer :: HTCTL(:)        => null() !cotyledon height, [m]
  integer,  pointer :: INTYP(:)        => null() !N2 fixation type
  integer,  pointer :: MY(:)           => null() !mycorrhizal type (no or yes)
  real(r8), pointer :: HTSTZ(:)        => null() !canopy height, [m]
  real(r8), pointer :: ZL(:)           => null() !canopy layer height , [m]
  real(r8), pointer :: ARSTP(:)        => null() !plant stem area, [m2 d-2]
  real(r8), pointer :: ARLFP(:)        => null() !plant leaf area, [m2 d-2]
  integer,  pointer :: NB1(:)          => null() !number of main branch
  integer,  pointer :: NI(:)           => null() !maximum soil layer number for all root axes
  integer,  pointer :: NIX(:)          => null() !maximum soil layer number for all root axes, [-]
  integer,  pointer :: NRT(:)          => null() !root primary axis number
  integer,  pointer :: NNOD(:)         => null() !number of concurrently growing nodes
  integer,  pointer :: NBT(:)          => null() !branch number
  integer,  pointer :: NBR(:)          => null() !branch number
  integer,  pointer :: NINR(:,:)       => null() !maximum soil layer number for root axes, [-]
  real(r8), pointer :: PSTG(:,:)       => null() !shoot node number, [-]
  real(r8), pointer :: PSTGI(:,:)      => null() !shoot node number at floral initiation, [-]
  real(r8), pointer :: PSTGF(:,:)      => null() !shoot node number at anthesis, [-]
  real(r8), pointer :: ANGBR(:)        => null() !branching angle, [degree from horizontal]
  real(r8), pointer :: SURF(:,:,:,:,:) => null() !leaf surface area, [m2 d-2]
  integer,  pointer :: KLEAFX(:,:)     => null() !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,  pointer :: NBTB(:,:)       => null() !branch number, [-]
  real(r8), pointer :: GRNXB(:,:)      => null() !branch potential grain number, [d-2]
  real(r8), pointer :: GRNOB(:,:)      => null() !branch grain number, [d-2]
  real(r8), pointer :: HTSHEX(:,:)     => null() !branch height, [m]
  real(r8), pointer :: ARLFZ(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: ARLFB(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: ANGSH(:)        => null() !sheath angle, [degree from horizontal]
  real(r8), pointer :: CLASS(:,:)      => null() !fractionction of leaves in different angle classes, [-]
  real(r8), pointer :: ARSTV(:,:)      => null() !canopy layer stem area, [m2 d-2]
  real(r8), pointer :: ARLFV(:,:)      => null() !canopy layer leaf area, [m2 d-2]
  real(r8), pointer :: ARLF1(:,:,:)    => null() !leaf area, [m2 d-2]
  real(r8), pointer :: GRNO(:)         => null() !canopy grain number, [d-2]
  real(r8), pointer :: SDPTH(:)        => null() !seeding depth, [m]
  real(r8), pointer :: SDPTHI(:)       => null() !planting depth, [m]
  real(r8), pointer :: SDLG(:)         => null() !seed length, [m]
  real(r8), pointer :: SDVL(:)         => null() !seed volume, [m3 ]
  real(r8), pointer :: SDAR(:)         => null() !seed surface area, [m2]
  real(r8), pointer :: ARSTT(:)        => null() !total stem area, [m2 d-2]
  real(r8), pointer :: SSL1(:)         => null() !petiole length:mass during growth, [m gC-1]
  real(r8), pointer :: SNL1(:)         => null() !internode length:mass during growth, [m gC-1]
  real(r8), pointer :: SLA1(:)         => null() !leaf area:mass during growth, [m2 gC-1]
  real(r8), pointer :: ARLFT(:)        => null() !total leaf area, [m2 d-2]
  real(r8), pointer :: ARLFS(:)        => null() !plant leaf area, [m2 d-2]
  real(r8), pointer :: HTNODX(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: HTSHE(:,:,:)    => null() !sheath height, [m]
  real(r8), pointer :: HTNODE(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: SURFB(:,:,:,:)  => null() !stem surface area, [m2 d-2]
  real(r8), pointer :: ARLFL(:,:,:,:)  => null() !layer leaf area, [m2 d-2]
  real(r8), pointer :: ARSTK(:,:,:)    => null() !stem layer area, [m2 d-2]
  real(r8), pointer :: CF(:)           => null() !clumping factor for self-shading in canopy layer, [-]
  real(r8), pointer :: XTLI(:)         => null() !number of nodes in seed, [-]
  real(r8), pointer :: CanopyHeight(:)           => null() !canopy height, [m]
  integer,  pointer :: NG(:)           => null() !soil layer at planting depth, [-]
  integer,  pointer :: KLEAF(:,:)      => null() !leaf number, [-]
  real(r8), pointer :: VSTG(:,:)       => null() !leaf number, [-]
  integer , pointer :: IRTYP(:)        => null() !grain type (below or above-ground)
  real(r8), pointer :: STMX(:)         => null() !maximum grain node number per branch, [-]
  real(r8), pointer :: SDMX(:)         => null() !maximum grain number per node , [-]
  real(r8), pointer :: WDLF(:)         => null() !leaf length:width ratio, [-]
  real(r8), pointer :: GRMX(:)         => null() !maximum grain size   , [g]
  real(r8), pointer :: RRAD1(:,:,:)    => null() !root layer diameter primary axes, [m ]
  real(r8), pointer :: RRAD2(:,:,:)    => null() !root layer diameter secondary axes, [m ]
  real(r8), pointer :: RRADP(:,:)      => null() !root internal radius, [m]
  real(r8), pointer :: RRAD1M(:,:)     => null() !maximum radius of primary roots, [m]
  real(r8), pointer :: RRAD2M(:,:)     => null() !maximum radius of secondary roots, [m]
  real(r8), pointer :: RSRR(:,:)       => null() !root radial resistivity, [MPa h m-2]
  real(r8), pointer :: RSRA(:,:)       => null() !root axial resistivity, [MPa h m-4]
  real(r8), pointer :: RTDNT(:)        => null() !total root length density, [m m-3]
  real(r8), pointer :: RTN1(:,:,:)     => null() !root layer number primary axes, [d-2]
  real(r8), pointer :: RTNL(:,:,:)     => null() !root layer number axes, [d-2]
  real(r8), pointer :: RTDNP(:,:,:)    => null() !root layer length density, [m m-3]
  real(r8), pointer :: RTVLP(:,:,:)    => null() !root layer volume air, [m2 d-2]
  real(r8), pointer :: RTVLW(:,:,:)    => null() !root layer volume water, [m2 d-2]
  real(r8), pointer :: RTARP(:,:,:)    => null() !root layer area per plant, [m p-1]
  contains
    procedure, public :: Init    => plt_morph_init
    procedure, public :: Destroy => plt_morph_destroy
  end type plant_morph_type

  type, public :: plant_pheno_type
  real(r8), pointer :: TFN4(:,:)   => null()     !root layer temperature growth functiom, [-]
  real(r8), pointer :: PTSHT(:)    => null()     !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(r8), pointer :: GFILL(:)    => null()     !maximum rate of fill per grain, [g h-1]
  real(r8), pointer :: OFFST(:)    => null()     !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
  real(r8), pointer :: PlantO2Stress(:)     => null()     !plant O2 stress indicator, []
  real(r8), pointer :: PB(:)       => null()     !branch nonstructural C content required for new branch, [gC gC-1]
  real(r8), pointer :: PR(:)       => null()     !threshold root nonstructural C content for initiating new root axis, [gC gC-1]
  real(r8), pointer :: RCELX(:,:,:) => null()    !element translocated from leaf during senescence, [g d-2 h-1]
  real(r8), pointer :: RCESX(:,:,:) => null()    !element translocated from sheath during senescence, [g d-2 h-1]
  real(r8), pointer :: TCZ(:)     => null()     !threshold temperature for spring leafout/dehardening, [oC]
  real(r8), pointer :: TCG(:)     => null()     !canopy growth temperature, [oC]
  real(r8), pointer :: TCX(:)     => null()     !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8), pointer :: TKG(:)     => null()     !canopy growth temperature, [K]
  real(r8), pointer :: TFN3(:)    => null()     !canopy temperature growth function, [-]
  real(r8), pointer :: WSTR(:)    => null()     !canopy plant water stress indicator, number of hours PSILT < PSILY, []
  real(r8), pointer :: CTC(:)     => null()     !temperature below which seed set is adversely affected, [oC]
  real(r8), pointer :: ZTYP(:)    => null()     !plant thermal adaptation zone, [-]
  real(r8), pointer :: ZTYPI(:)   => null()     !initial plant thermal adaptation zone, [-]
  real(r8), pointer :: HTC(:)     => null()     !temperature above which seed set is adversely affected, [oC]
  real(r8), pointer :: SSTX(:)    => null()     !sensitivity to HTC (seeds oC-1 above HTC)
  integer,  pointer :: IDTH(:)    => null()     !flag for species death
  integer,  pointer :: IFLGC(:)   => null()     !flag for living pft
  real(r8), pointer :: RSETE(:,:) => null()     !effect of canopy element status on seed set , []
  real(r8), pointer :: XRNI(:)    => null()     !rate of node initiation, [h-1 at 25 oC]
  real(r8), pointer :: XRLA(:)    => null()     !rate of leaf initiation, [h-1 at 25 oC]
  real(r8), pointer :: XDL(:)     => null()     !critical daylength for phenological progress, [h]
  real(r8), pointer :: XPPD(:)    => null()     !difference between current and critical daylengths used to calculate  phenological progress [h]
  integer,  pointer :: IDTHB(:,:)  => null()  !flag to detect branch death , [-]
  integer,  pointer :: IFLGP(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGF(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGE(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGA(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGG(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGR(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGQ(:,:)  => null()  !branch phenology flag, [h]
  integer,  pointer :: KVSTG(:,:)  => null()  !leaf growth stage counter, [-]
  integer,  pointer :: IWTYP(:)    => null()  !climate signal for phenological progress: none, temperature, water stress
  integer,  pointer :: IBTYP(:)    => null()  !phenologically-driven above-ground turnover: all, foliar only, none
  integer,  pointer :: IDTHP(:)    => null()  !flag to detect canopy death
  integer,  pointer :: ISTYP(:)    => null()  !plant growth habit: annual or perennial
  integer,  pointer :: IDTHR(:)    => null()  !flag to detect root system death
  integer,  pointer :: IDTYP(:)    => null()  !plant growth habit (determinate or indeterminate)
  integer,  pointer :: IPTYP(:)    => null()  !photoperiod type (neutral, long day, short day)
  integer,  pointer :: IGTYP(:)    => null()  !plant growth type (vascular, non-vascular)
  integer,  pointer :: IFLGI(:)    => null()  !PFT initialization flag:0=no,1=yes
  integer,  pointer :: KVSTGN(:,:) => null()  !leaf growth stage counter, [-]
  integer,  pointer :: IDAY(:,:,:) => null()  !plant growth stage, [-]
  real(r8), pointer :: TGSTGI(:,:) => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: TGSTGF(:,:) => null()  !normalized node number during reproductive growth stages , [-]
  real(r8), pointer :: VSTGX(:,:)  => null()  !leaf number at floral initiation, [-]
  real(r8), pointer :: VRNY(:,:)   => null()  !initial heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: VRNZ(:,:)   => null()  !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: VRNS(:,:)   => null()  !heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: VRNL(:,:)   => null()  !hours above threshold temperature required for spring leafout/dehardening, [-]
  real(r8), pointer :: VRNF(:,:)   => null()  !cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: VRNX(:,:)   => null()  !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8), pointer :: ATRP(:,:)   => null()  !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8), pointer :: FLGZ(:,:)   => null()  !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  real(r8), pointer :: DGSTGI(:,:) => null()  !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8), pointer :: DGSTGF(:,:) => null()  !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8), pointer :: GROUP(:,:)  => null()  !plant maturity group, [-]
  real(r8), pointer :: GSTGI(:,:)  => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: GSTGF(:,:)  => null()  !normalized node number during reproductive growth stages, [-]
  real(r8), pointer :: FLG4(:,:)   => null()  !flag to detect physiological maturity from  grain fill , [-]
  real(r8), pointer :: GROUPI(:)   => null()  !acclimated plant maturity group, [-]
  contains
    procedure, public :: Init    =>  plt_pheno_init
    procedure, public :: Destroy =>  plt_pheno_destroy
  end type plant_pheno_type

  type, public :: plant_soilchem_type
  real(r8), pointer :: FOSRH(:,:)  => null()  !fraction of total organic C in complex, [-]
  real(r8), pointer :: CFOPE(:,:,:,:)=> null() !litter kinetic fraction, [-]
  real(r8), pointer :: TFND(:)     => null()  !temperature effect on diffusivity
  real(r8), pointer :: THETPM(:,:) => null()  !soil air-filled porosity, [m3 m-3]
  real(r8), pointer :: DFGS(:,:)   => null()  !coefficient for dissolution - volatilization, []
  real(r8), pointer :: RSCS(:)     => null()  !soil hydraulic resistance, [MPa h m-2]
  real(r8), pointer :: BKDS(:)     => null()  !soil bulk density, [Mg m-3]
  real(r8), pointer :: trc_solcl(:,:) => null() !aqueous tracer concentration [g m-3]
  real(r8), pointer :: trc_gascl(:,:) => null() !gaseous tracer concentration [g m-3]

  real(r8), pointer :: CORGC(:)    => null()  !soil organic C content [gC kg soil-1]

  real(r8), pointer :: CNDU(:)     => null()  !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]

  real(r8), pointer :: GasDifc(:,:)=> null()  !gaseous diffusivity [m2 h-1]
  real(r8), pointer :: SolDifc(:,:)=> null()  !aqueous diffusivity [m2 h-1]

  real(r8), pointer :: trc_gasml(:,:)=> null()!gas layer mass [g d-2]

  real(r8), pointer :: GSolbility(:,:)=> null() !gas solubility, [m3 m-3]

  real(r8), pointer :: THETW(:)    => null()  !volumetric water content [m3 m-3]
  real(r8), pointer :: THETY(:)    => null()  !air-dry water content, [m3 m-3]
  real(r8), pointer :: VOLX(:)     => null()  !volume of soil layer	m3 d-2
  real(r8), pointer :: trcs_VLN(:,:)=> null()

  real(r8), pointer :: VOLY(:)     => null()  !total micropore volume [m3 d-2]
  real(r8), pointer :: VOLI(:)     => null()  !soil micropore ice content   [m3 d-2]
  real(r8), pointer :: VOLW(:)     => null()  !soil micropore water content [m3 d-2]
  real(r8), pointer :: VOLA(:)     => null()  !total volume in micropores [m3 d-2]

  real(r8), pointer :: OQC(:,:)    => null()  !dissolved organic C micropore	[gC d-2]
  real(r8), pointer :: OQN(:,:)    => null()  !dissolved organic N micropore	[gN d-2]
  real(r8), pointer :: OQP(:,:)    => null()  !dissolved organic P micropore	[gP d-2]
  real(r8), pointer :: trc_solml(:,:)=> null() !aqueous tracer [g d-2]
  contains
    procedure, public :: Init => plt_soilchem_init
    procedure, public :: Destroy => plt_soilchem_destroy
  end type plant_soilchem_type

  type, public :: plant_allometry_type
  real(r8), pointer :: CPRTS(:)    => null()  !root P:C ratio x root growth yield, [-]
  real(r8), pointer :: CNRTS(:)    => null()  !root N:C ratio x root growth yield, [-]
  real(r8), pointer :: CNND(:)     => null()  !nodule N:C ratio, [gN gC-1]
  real(r8), pointer :: CPND(:)     => null()  !nodule P:C ratio, [gP gC-1]
  real(r8), pointer :: CNRT(:)     => null()  !root N:C ratio, [gN gC-1]
  real(r8), pointer :: CPRT(:)     => null()  !root P:C ratio, [gP gC-1]
  real(r8), pointer :: CNWS(:)     => null()  !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: CPWS(:)     => null()  !C:P ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: DMND(:)     => null()  !nodule growth yield, [g g-1]
  real(r8), pointer :: DMRT(:)     => null()  !root growth yield, [g g-1]
  real(r8), pointer :: CPEAR(:)    => null()  !ear P:C ratio, [gP gC-1]
  real(r8), pointer :: DMSHE(:)    => null()  !sheath growth yield, [g g-1]
  real(r8), pointer :: CPHSK(:)    => null()  !husk P:C ratio, [gP gC-1]
  real(r8), pointer :: DMSTK(:)    => null()  !stalk growth yield, [gC gC-1]
  real(r8), pointer :: DMHSK(:)    => null()  !husk growth yield, [gC gC-1]
  real(r8), pointer :: DMRSV(:)    => null()  !reserve growth yield, [gC gC-1]
  real(r8), pointer :: DMGR(:)     => null()  !grain growth yield, [gC gC-1]
  real(r8), pointer :: DMEAR(:)    => null()  !ear growth yield, [gC gC-1]
  real(r8), pointer :: CNHSK(:)    => null()  !husk N:C ratio, [gN gC-1]
  real(r8), pointer :: CNRSV(:)    => null()  !reserve N:C ratio, [gN gC-1]
  real(r8), pointer :: CNEAR(:)    => null()  !ear N:C ratio, [gN gC-1]
  real(r8), pointer :: CPRSV(:)    => null()  !reserve P:C ratio, [gP gC-1]
  real(r8), pointer :: CPGR(:)     => null()  !grain P:C ratio, [gP gP-1]
  real(r8), pointer :: CNSTK(:)    => null()  !stalk N:C ratio, [gN gC-1]
  real(r8), pointer :: FWODBE(:,:) => null()  !woody element allocation
  real(r8), pointer :: FWODLE(:,:) => null()  !leaf element allocation
  real(r8), pointer :: FWODRE(:,:) => null()  !C woody fraction in root
  real(r8), pointer :: FWOODE(:,:) => null()  !woody element allocation
  real(r8), pointer :: FVRN(:)     => null()  !allocation parameter
  real(r8), pointer :: DMLF(:)     => null()  !leaf growth yield, [g g-1]
  real(r8), pointer :: CNGR(:)     => null()  !grain N:C ratio, [g g-1]
  real(r8), pointer :: CPLF(:)     => null()  !maximum leaf P:C ratio, [g g-1]
  real(r8), pointer :: CPSHE(:)    => null()  !sheath P:C ratio, [g g-1]
  real(r8), pointer :: CNSHE(:)    => null()  !sheath N:C ratio, [g g-1]
  real(r8), pointer :: CPSTK(:)    => null()  !stalk P:C ratio, [g g-1]
  real(r8), pointer :: CNLF(:)     => null()  !maximum leaf N:C ratio, [g g-1]
  real(r8), pointer :: GRWTB(:,:)  => null()  !maximum grain C during grain fill, [g d-2]
  real(r8), pointer :: FNOD(:)     => null()  !parameter for allocation of growth to nodes, [-]
  real(r8), pointer :: CWSRT(:)    => null()  !fraction of remobilizable nonstructural biomass in root, [-]

  contains
    procedure, public :: Init => plt_allom_init
    procedure, public :: Destroy => plt_allom_destroy
  end type plant_allometry_type

  type, public :: plant_biom_type
  real(r8), pointer :: WTSTGET(:)     => null()    !total standing dead element, [g d-2]
  real(r8), pointer :: ZEROL(:)       => null()    !threshold zero for leaf calculation
  real(r8), pointer :: ZEROP(:)       => null()    !threshold zero for p calculation
  real(r8), pointer :: EPOOLN(:,:,:)  => null()    !root  layer nonstructural element, [g d-2]
  real(r8), pointer :: WTNDLE(:,:,:)  => null()    !root layer nodule element, [g d-2]
  real(r8), pointer :: WGLFV(:,:)     => null()    !canopy layer leaf C, [g d-2]
  real(r8), pointer :: WTRT2E(:,:,:,:,:) => null()    !root layer element secondary axes, [g d-2]
  real(r8), pointer :: WTRT1E(:,:,:,:,:) => null()    !root layer element primary axes, [g d-2]
  real(r8), pointer :: RTWT1E(:,:,:,:)   => null()    !root C primary axes, [g d-2]
  real(r8), pointer :: WTSTDE(:,:,:)  => null()    !standing dead element fraction, [g d-2]
  real(r8), pointer :: CEPOLP(:,:)    => null()    !canopy nonstructural element concentration, [g d-2]
  real(r8), pointer :: EPOOLP(:,:)    => null()    !canopy nonstructural element concentration, [g d-2]
  real(r8), pointer :: EPOLNP(:,:)    => null()    !canopy nodule nonstructural element, [g d-2]
  real(r8), pointer :: CCPLNP(:)      => null()    !nodule nonstructural C, [gC d-2]
  real(r8), pointer :: WTRTL(:,:,:)   => null()    !root layer structural C, [g d-2]
  real(r8), pointer :: WTRTD(:,:,:)   => null()    !root layer C, [g d-2]
  real(r8), pointer :: WSRTL(:,:,:)   => null()    !root layer protein C, [g d-2]
  real(r8), pointer :: CWSRTL(:,:,:)  => null()    !root layer protein C concentration, [g g-1]
  real(r8), pointer :: EPOOLR(:,:,:,:)=> null()    !root  layer nonstructural element, [g d-2]
  real(r8), pointer :: CEPOLR(:,:,:,:)  => null()    !root  layer nonstructural C concentration, [g g-1]
  real(r8), pointer :: CEPOLB(:,:,:)    => null()    !branch nonstructural C concentration, [g d-2]
  real(r8), pointer :: WGNODE(:,:,:,:)  => null()    !internode C, [g d-2]
  real(r8), pointer :: WGLFE(:,:,:,:)    => null()    !leaf element, [g d-2]
  real(r8), pointer :: WSLF(:,:,:)    => null()    !layer leaf protein C, [g d-2]
  real(r8), pointer :: WGSHE(:,:,:,:)   => null()  !sheath element , [g d-2]
  real(r8), pointer :: WSSHE(:,:,:)   => null()    !layer sheath protein C, [g d-2]
  real(r8), pointer :: WGLFLE(:,:,:,:,:) => null()    !layer leaf element, [g d-2]
  real(r8), pointer :: WGLFT(:)       => null()  !total leaf mass, [gC d-2]
  real(r8), pointer :: WTSTDI(:)      => null()  !initial standing dead C, [g C m-2]
  real(r8), pointer :: WTRTE(:,:)     => null()  !plant root element, [gC d-2]
  real(r8), pointer :: WTRTSE(:,:)    => null()  !plant root structural element, [gC d-2]
  real(r8), pointer :: WTRVX(:)       => null()  !plant stored nonstructural C at planting, [gC d-2]
  real(r8), pointer :: WTRVE(:,:)     => null()  !plant stored nonstructural element, [gC d-2]
  real(r8), pointer :: WTLS(:)        => null()  !canopy leaf + sheath C, [g d-2]
  real(r8), pointer :: WTSHTE(:,:)    => null()  !canopy shoot C, [g d-2]
  real(r8), pointer :: WTSHTA(:)      => null()  !landscape average canopy shoot C, [g d-2]
  real(r8), pointer :: WTSTGE(:,:)    => null()  !standing dead element, [g d-2]
  real(r8), pointer :: WTNDE(:,:)     => null()  !root total nodule mass, element [g d-2]
  real(r8), pointer :: EPOOL(:,:,:)   => null()  !branch nonstructural element, [g d-2]
  real(r8), pointer :: EPOLNB(:,:,:)  => null()  !branch nodule nonstructural element, [g d-2]
  real(r8), pointer :: WTLSB(:,:)     => null()  !branch leaf + sheath C, [g d-2]
  real(r8), pointer :: WTRSVBE(:,:,:) => null()  !branch reserve element, [g d-2]
  real(r8), pointer :: WTLFBE(:,:,:)  => null()   !branch leaf element, [g d-2]
  real(r8), pointer :: WTNDBE(:,:,:)  => null()   !branch nodule element, [g d-2]
  real(r8), pointer :: WTSHEBE(:,:,:) => null()   !branch sheath element , [g d-2]
  real(r8), pointer :: WTEARBE(:,:,:) => null()   !branch ear C, [g d-2]
  real(r8), pointer :: WTHSKBE(:,:,:) => null()   !branch husk element, [g d-2]
  real(r8), pointer :: WTGRBE(:,:,:)  => null()   !branch grain element, [g d-2]
  real(r8), pointer :: WTSTKBE(:,:,:) => null()   !branch stalk element, [g d-2]
  real(r8), pointer :: WTSHTBE(:,:,:) => null()   !branch shoot C, [g d-2]
  real(r8), pointer :: WGLFEX(:,:,:)  => null()   !branch leaf structural element, [g d-2]
  real(r8), pointer :: WGSHEXE(:,:,:) => null()   !branch sheath structural element, [g d-2]
  real(r8), pointer :: WTSTXBE(:,:,:) => null()   !branch stalk structural element, [g d-2]
  real(r8), pointer :: WVSTKB(:,:)    => null()   !branch active stalk C, [g d-2]
  real(r8), pointer :: WTSTKE(:,:)    => null()   !canopy stalk element, [g d-2]
  real(r8), pointer :: WVSTK(:)       => null()   !canopy active stalk C, [g d-2
  real(r8), pointer :: WTLFE(:,:)     => null()   !canopy leaf elements, [g d-2]
  real(r8), pointer :: WTSHEE(:,:)    => null()   !canopy sheath element , [g d-2]
  real(r8), pointer :: WTRSVE(:,:)    => null()   !canopy reserve element, [g d-2]
  real(r8), pointer :: WTHSKE(:,:)    => null()   !canopy husk element, [g d-2]
  real(r8), pointer :: WTRTA(:)       => null()   !root C per plant, [g p-1]
  real(r8), pointer :: WTGRE(:,:)     => null()   !canopy grain element, [g d-2]
  real(r8), pointer :: WTEARE(:,:)    => null()   !canopy ear element, [g d-2]
  contains
    procedure, public :: Init => plt_biom_init
    procedure, public :: Destroy => plt_biom_destroy
  end type plant_biom_type

  type, public :: plant_ew_type
  real(r8) :: DPTHS     !snowpack depth, [m]
  real(r8) :: VOLWS     !water volume in snowpack, [m3 d-2]
  real(r8) :: VOLSS     !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8) :: TSHC      !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VHCPWX    !snowpack heat capacity from previous time step, [MJ d-2 K-1]
  real(r8) :: VHCPW1    !snowpack heat capacity, [MJ m-3 K-1]
  real(r8) :: TENGYC    !total canopy heat content, [MJ  d-2]
  real(r8) :: THRMC     !total canopy LW emission, [MJ d-2 h-1]
  real(r8) :: TEVAPP    !total canopy evaporation + transpiration, [m3 d-2]
  real(r8) :: TEVAPC    !total canopy evaporation, [m3 d-2]
  real(r8) :: THFLXC    !total canopy heat flux, [MJ  d-2]
  real(r8) :: UVOLO     !total subsurface water flux, [m3 d-2]
  real(r8) :: VPA       !vapor concentration, [m3 m-3]
  real(r8) :: TairK       !air temperature, [K]
  real(r8) :: TVOLWP    !total canopy water content, [m3 d-2]
  real(r8) :: TLE       !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8) :: TLEC      !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VOLIS     !ice volume in snowpack, [m3 d-2]
  real(r8) :: TKW       !snow temperature, [K]
  real(r8) :: TVOLWC    !canopy surface water content, [m3 d-2]
  real(r8) :: TSH       !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8) :: ZR        !canopy surface roughness height, [m]
  real(r8) :: ZD        !zero plane displacement height, [m]
  real(r8) :: RAB       !isothermal boundary layer resistance, [h m-1]
  real(r8) :: RIB       !Richardson number for calculating boundary layer resistance, [-]
  real(r8) :: TGH       !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8), pointer :: EP(:)     => null()    !canopy transpiration, [m2 d-2 h-1]
  real(r8), pointer :: PSILG(:)  => null()    !canopy turgor water potential, [MPa]
  real(r8), pointer :: FLWC(:)   => null()    !water flux into canopy, [m3 d-2 h-1]
  real(r8), pointer :: PSILT(:)  => null()    !canopy total water potential , [Mpa]
  real(r8), pointer :: EVAPC(:)  => null()    !canopy evaporation, [m2 d-2 h-1]
  real(r8), pointer :: HFLXC(:)  => null()    !canopy storage heat flux, [MJ d-2 h-1]
  real(r8), pointer :: EFLXC(:)  => null()    !canopy latent heat flux, [MJ d-2 h-1]
  real(r8), pointer :: RAZ(:)    => null()    !canopy roughness height, [m]
  real(r8), pointer :: TKS(:)    => null()    !mean annual soil temperature, [K]
  real(r8), pointer :: PSILZ(:)  => null()    !minimum daily canopy water potential, [MPa]
  real(r8), pointer :: TCC(:)    => null()    !canopy temperature, [oC]
  real(r8), pointer :: DTKC(:)   => null()    !change in canopy temperature, [K]
  real(r8), pointer :: ENGYX(:)  => null()    !canopy heat storage from previous time step, [MJ d-2]
  real(r8), pointer :: TKC(:)    => null()    !canopy temperature, [K]
  real(r8), pointer :: PSILO(:)  => null()    !canopy osmotic water potential, [Mpa]
  real(r8), pointer :: OSMO(:)   => null()    !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8), pointer :: SFLXC(:)  => null()    !canopy sensible heat flux, [MJ d-2 h-1]
  real(r8), pointer :: TKCZ(:)   => null()    !canopy temperature, [K]
  real(r8), pointer :: PSIRT(:,:,:) => null() !root total water potential , [Mpa]
  real(r8), pointer :: PSIRO(:,:,:) => null() !root osmotic water potential , [Mpa]
  real(r8), pointer :: PSIRG(:,:,:) => null() !root turgor water potential , [Mpa]
  real(r8), pointer :: PopPlantRootH2OUptake_vr(:,:,:) => null() !root water uptake, [m2 d-2 h-1]
  real(r8), pointer :: TUPHT(:)     => null()   !total root heat uptake, [MJ d-2]
  real(r8), pointer :: GridPlantRootH2OUptake_vr(:)    => null()   !total root water uptake, [m3 d-2]
  real(r8), pointer :: VOLWC(:)     => null()  !canopy surface water content, [m3 d-2]
  real(r8), pointer :: VOLWP(:)     => null()  !canopy water content, [m3 d-2]
  real(r8), pointer :: VHCPC(:)     => null()  !canopy heat capacity, [MJ d-2 K-1]
  real(r8), pointer :: TotalSoilH2OPSIMPa(:)     => null()  !soil micropore total water potential [MPa]
  real(r8), pointer :: CTRAN(:)     =>  null()  !total transpiration, [m H2O d-2]
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
  integer,  pointer :: IYR0(:)     => null()  !year of planting
  integer,  pointer :: IDAYX(:)    => null()  !alternate day of planting
  integer,  pointer :: IYRX(:)     => null()  !alternate year of planting
  integer,  pointer :: IYRY(:)     => null()  !alternate year of harvest
  integer,  pointer :: IDAY0(:)    => null()  !day of planting
  integer,  pointer :: IDAYH(:)    => null()  !day of harvest
  integer,  pointer :: IYRH(:)     => null()  !year of harvest
  integer,  pointer :: IHVST(:)    => null()  !type of harvest
  integer,  pointer :: JHVST(:)    => null()  !flag for stand replacing disturbance
  real(r8), pointer :: EHVST(:,:,:)=> null()  !harvest efficiency, [-]
  real(r8), pointer :: HVST(:)     => null()  !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  real(r8), pointer :: THIN(:)     => null()  !thinning of plant population, [-]
  real(r8), pointer :: VCH4F(:)    => null()  !plant CH4 emission from fire, [g d-2 ]
  real(r8), pointer :: VCO2F(:)    => null()  !plant CO2 emission from fire, [g d-2 ]
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
  real(r8) :: TNBP      !total NBP, [g d-2]
  real(r8) :: TGPP      !ecosystem GPP, [g d-2 h-1]

  real(r8) :: CNETX     !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8) :: RECO      !ecosystem respiration, [g d-2 h-1]
  real(r8) :: TRAU      !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8) :: TH2GZ     !total root H2 flux, [g d-2]
  real(r8) :: TCCAN     !total net CO2 fixation, [gC d-2]
  real(r8), pointer :: ZESNC(:) => null() !total litterfall element, [g d-2 h-1]
  real(r8), pointer :: ZNPP(:)       => null()   !total net primary productivity, [gC d-2]
  real(r8), pointer :: RNH3C(:)      => null()   !canopy NH3 flux, [g d-2 h-1]
  real(r8), pointer :: TDFOME(:,:,:)   =>  null()  !total root element exchange, [g d-2 h-1]
  real(r8), pointer :: RUPNF(:,:)    =>  null()  !root N2 fixation, [gN d-2 h-1]
  real(r8), pointer :: TCO2A(:)      =>  null()  !total autotrophic respiration, [gC d-2 ]
  real(r8), pointer :: ESNC(:,:,:,:,:) =>  null()  !plant litterfall element, [g d-2 h-1]
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
  real(r8), pointer :: TUPNO3(:)     => null()   !total root-soil NO3 flux non-band, [gN d-2 h-1]
  real(r8), pointer :: TUPH2B(:)     => null()   !total root-soil PO4 flux band, [gP d-2 h-1]
  real(r8), pointer :: TUPH1B(:)     => null()   !soil-root exch of HPO4 in band [gP d-2 h-1]
  real(r8), pointer :: TUPN2S(:)     => null()   !total root-soil N2O flux, [gN d-2 h-1]
  real(r8), pointer :: TCO2S(:)      => null()   !total root-soil CO2 flux, [gC d-2 h-1]
  real(r8), pointer :: TCO2P(:)      => null()   !total root CO2 flux, [gC d-2 h-1]
  real(r8), pointer :: TUPOXP(:)     => null()   !total root internal O2 flux, [g d-2 h-1]
  real(r8), pointer :: TUPOXS(:)     => null()   !total root-soil O2 flux, [g d-2 h-1]
  real(r8), pointer :: TUPHGS(:)     => null()   !total root-soil H2 flux, [g d-2 h-1]
  real(r8), pointer :: TUPCHS(:)     => null()   !total root-soil CH4 flux, [gC d-2 h-1]
  real(r8), pointer :: TUPN3B(:)     => null()   !total root-soil NH3 flux band, [gN d-2 h-1]
  real(r8), pointer :: TUPN3S(:)     => null()   !total root-soil NH3 flux non-band, [gN d-2 h-1]
  real(r8), pointer :: TUPH2P(:)     => null()   !total root-soil PO4 flux non-band, [gP d-2 h-1]
  real(r8), pointer :: TUPH1P(:)     => null()   !soil-root exch of HPO4 in non-band, [gP d-2 h-1]
  real(r8), pointer :: TUPNH4(:)     => null()   !total root-soil NH4 flux non-band, [gN d-2 h-1]
  real(r8), pointer :: TUPNHB(:)     => null()   !total root-soil NH4 flux band, [gN d-2 h-1]
  real(r8), pointer :: TUPNOB(:)     => null()   !total root-soil NO3 flux band, [gN d-2 h-1]
  real(r8), pointer :: ESNT(:,:,:,:)   => null() !total litterfall element, [g d-2 h-1]
  real(r8), pointer :: XOQCS(:,:)    => null()  !net microbial DOC flux, [gC d-2 h-1]
  real(r8), pointer :: XOQNS(:,:)    => null()  !net microbial DON flux, [gN d-2 h-1]
  real(r8), pointer :: XOQPS(:,:)    => null()  !net microbial DOP flux, [gP d-2 h-1]
  real(r8), pointer :: CNET(:)       => null()  !canopy net CO2 exchange, [gC d-2 h-1]
  real(r8), pointer :: CARBN(:)      => null()  !total gross CO2 fixation, [gC d-2 ]
  real(r8), pointer :: HESNC(:,:)    => null()  !plant element litterfall, [g d-2 h-1]
  real(r8), pointer :: RFGas_root(:,:)=> null() !gaseous flux fron root disturbance, [g d-2 h-1]
  real(r8), pointer :: TESN0(:,:)    => null()  !total surface litterfall element, [g d-2]
  real(r8), pointer :: TESNC(:,:)    => null()  !total plant element litterfall , [g d-2 ]
  real(r8), pointer :: TCO2T(:)      => null()  !total plant respiration, [gC d-2 ]

  real(r8), pointer :: TNH3C(:)      => null()  !total canopy NH3 flux, [gN d-2 ]
  real(r8), pointer :: TZUPFX(:)     => null()  !total plant N2 fixation, [g d-2 ]

  contains
    procedure, public :: Init  => plt_bgcrate_init
    procedure, public :: Destroy  => plt_bgcrate_destroy
  end type plant_bgcrate_type

  type, public :: plant_rootbgc_type
  real(r8), pointer :: TRFGas_root(:)   => null()  !total root gas content [g d-2]
  real(r8), pointer :: UPOME(:,:)       => null()  !total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
  real(r8), pointer :: UPNF(:)          => null()  !total root N2 fixation, [g d-2 h-1]
  real(r8), pointer :: UPNO3(:)         => null()  !total root uptake of NO3, [g d-2 h-1]
  real(r8), pointer :: UPNH4(:)         => null()  !total root uptake of NH4, [g d-2 h-1]
  real(r8), pointer :: UPH1P(:)         => null()  !total root uptake of HPO4, [g d-2 h-1]
  real(r8), pointer :: UPH2P(:)         => null()  !total root uptake of PO4, [g d-2 h-1]
  real(r8), pointer :: RDFOME(:,:,:,:,:)  => null()  !root uptake (+ve) - exudation (-ve) of DOE, [g d-2 h-1]
  real(r8), pointer :: HEUPTK(:,:)      => null()  !net root element uptake (+ve) - exudation (-ve), [gC d-2 h-1]
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
  real(r8), pointer :: trcg_RFLA(:,:,:,:)    => null()  !gaseous tracer flux through roots, [g d-2 h-1]
  real(r8), pointer :: trcg_RDFA(:,:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
  real(r8), pointer :: ROXYP(:,:,:)     => null()  !root  O2 demand from respiration, [g d-2 h-1]
  real(r8), pointer :: RUNNHP(:,:,:)    => null()  !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8), pointer :: RUNNBP(:,:,:)    => null()  !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8), pointer :: RUNNOP(:,:,:)    => null()  !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8), pointer :: RUNNXP(:,:,:)    => null()  !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8), pointer :: RUPP2P(:,:,:)    => null()  !root uptake of H2PO4 non-band
  real(r8), pointer :: RUPP2B(:,:,:)    => null()  !root uptake of H2PO4 band
  real(r8), pointer :: RUPP1P(:,:,:)    => null()  !HPO4 demand in non-band by each root population
  real(r8), pointer :: RUPP1B(:,:,:)    => null()  !HPO4 demand in band by each root population
  real(r8), pointer :: WFR(:,:,:)       => null()  !O2 constraint to root respiration, []
  real(r8), pointer :: trcg_rootml(:,:,:,:)=> null() !root gas content, [g d-2]
  real(r8), pointer :: trcs_rootml(:,:,:,:)=> null() !root aqueous content, [g d-2]
  real(r8), pointer :: RNH3B(:,:)       => null()  !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  real(r8), pointer :: RUPNHB(:,:,:)    => null()  !root uptake of NH4 band, [g d-2 h-1]
  real(r8), pointer :: RUPNH4(:,:,:)    => null()  !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8), pointer :: RUPH2P(:,:,:)    => null()  !root uptake of PO4 non-band, [g d-2 h-1]
  real(r8), pointer :: RUPNOB(:,:,:)    => null()  !root uptake of NO3 band, [g d-2 h-1]
  real(r8), pointer :: RUPNO3(:,:,:)    => null()  !root uptake of NO3 non-band, [g d-2 h-1]
  real(r8), pointer :: RUPH1B(:,:,:)    => null()  !root uptake of HPO4 band
  real(r8), pointer :: RUPH1P(:,:,:)    => null()  !root uptake HPO4 non-band
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
  real(r8), pointer :: RCO2M(:,:,:)     => null()  !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8), pointer :: RCO2N(:,:,:)     => null()  !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), pointer :: RCO2A(:,:,:)     => null()  !root respiration constrained by O2, [g d-2 h-1]
  real(r8), pointer :: TEUPTK(:,:)      => null()  !total net root element uptake (+ve) - exudation (-ve), [gC d-2 ]
  real(r8), pointer :: trcg_TLP(:,:)    => null()   !total root internal gas flux, [g d-2 h-1]
  real(r8), pointer :: trcg_TFLA(:,:)   => null()   !total internal root gas flux , [gC d-2 h-1]
  real(r8), pointer :: TUPNF(:)         => null()   !total root N2 fixation, [g d-2 h-1]

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

  allocate(this%trcg_rootml(idg_beg:idg_end-1,2,JZ1,JP1));this%trcg_rootml=0._r8
  allocate(this%trcs_rootml(idg_beg:idg_end-1,2,JZ1,JP1));this%trcs_rootml=0._r8
  allocate(this%TRFGas_root(idg_beg:idg_end-1));this%TRFGas_root=0._r8
  allocate(this%TUPNF(JZ1))
  allocate(this%ROXSK(60,0:JZ1))
  allocate(this%RDFOME(npelms,2,1:jcplx,0:JZ1,JP1))
  allocate(this%HEUPTK(npelms,JP1))
  allocate(this%TEUPTK(npelms,JP1))
  allocate(this%UPOME(npelms,JP1))
  allocate(this%UPNF(JP1))
  allocate(this%UPNO3(JP1))
  allocate(this%UPNH4(JP1))
  allocate(this%UPH1P(JP1))
  allocate(this%UPH2P(JP1))

  allocate(this%ZEROQ(JP1))
  allocate(this%RCO2M(2,JZ1,JP1))
  allocate(this%RCO2N(2,JZ1,JP1))
  allocate(this%RCO2A(2,JZ1,JP1))
  allocate(this%RUPNHB(2,JZ1,JP1))
  allocate(this%RUPNH4(2,JZ1,JP1))
  allocate(this%RUPH2P(2,JZ1,JP1))
  allocate(this%RUPNOB(2,JZ1,JP1))
  allocate(this%RUPNO3(2,JZ1,JP1))
  allocate(this%RUPH1B(2,JZ1,JP1))
  allocate(this%RUPH1P(2,JZ1,JP1))
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
  allocate(this%RNH3B(JBR,JP1))


  allocate(this%trcg_TFLA(idg_beg:idg_end-1,JZ1))
  allocate(this%trcg_TLP(idg_beg:idg_end-1,JZ1))

  allocate(this%trcg_RFLA(idg_beg:idg_end-1,2,JZ1,JP1))
  allocate(this%trcg_RDFA(idg_beg:idg_end-1,2,JZ1,JP1))
  allocate(this%ROXYP(2,JZ1,JP1))
  allocate(this%RUNNHP(2,JZ1,JP1))
  allocate(this%RUNNBP(2,JZ1,JP1))
  allocate(this%RUNNOP(2,JZ1,JP1))
  allocate(this%RUNNXP(2,JZ1,JP1))
  allocate(this%RUPP2P(2,JZ1,JP1))
  allocate(this%RUPP2B(2,JZ1,JP1))
  allocate(this%RUPP1P(2,JZ1,JP1))
  allocate(this%RUPP1B(2,JZ1,JP1))

  allocate(this%WFR(2,JZ1,JP1))
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
!  if(allocated(TUPNF))deallocate(TUPNF)
!  if(allocated(ROXSK))deallocate(ROXSK)
!  if(allocated(RDFOME))deallocate(RDFOME)
!  if(allocated(HEUPTK))deallocate(HEUPTK)
!  if(allocated(TEUPTK))deallocate(TEUPTK)
!  if(allocated(UPOME))deallocate(UPOME)
!  if(allocated(UPNF))deallocate(UPNF)
!  if(allocated(UPNO3))deallocate(UPNO3)
!  if(allocated(UPNH4))deallocate(UPNH4)
!  if(allocated(UPH1P))deallocate(UPH1P)
!  if(allocated(UPH2P))deallocate(UPH2P)
!  if(allocated(RNH3B))deallocate(RNH3B)
!  if(allocated(ZEROQ))deallocate(ZEROQ)
!  if(allocated(RCO2M))deallocate(RCO2M)
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
!  if(allocated(RUPNH4))deallocate(RUPNH4)
!  if(allocated(RUPH2P))deallocate(RUPH2P)
!  if(allocated(RUPNOB))deallocate(RUPNOB)
!  if(allocated(RUPNO3))deallocate(RUPNO3)
!  if(allocated(RUPH1B))deallocate(RUPH1B)
!  if(allocated(RUPH1P))deallocate(RUPH1P)
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
!  if(allocated(WFR))deallocate(WFR)
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


  allocate(this%TBALE(npelms))
  allocate(this%FMPR(0:JZ1))
  allocate(this%DATAP(JP1))
  allocate(this%DATA(30))
  allocate(this%AREA3(0:JZ1))
  allocate(this%IDATA(60))
  allocate(this%DLYR3(0:JZ1))
  allocate(this%BALE(npelms,JP1))
  allocate(this%CDPTHZ(0:JZ1))
  allocate(this%DPTHZ(0:JZ1))
  allocate(this%PPI(JP1))
  allocate(this%PPZ(JP1))
  allocate(this%PPX(JP1))
  allocate(this%pftPlantPopulation(JP1))
  allocate(this%VOLWM(60,0:JZ1))
  allocate(this%VOLPM(60,0:JZ1))
  allocate(this%TORT(60,0:JZ1))
  allocate(this%FILM(60,0:JZ1))

  end subroutine plt_site_Init
!----------------------------------------------------------------------
  subroutine plt_site_destroy(this)
  implicit none
  class(plant_siteinfo_type) :: this



!  if(allocated(DLYR3))deallocate(DLYR3)


!  if(allocated(FMPR))deallocate(FMPR)
!  if(allocated(TORT))deallocate(TORT)
!  if(allocated(FILM))deallocate(FILM)
!  if(allocated(DATAP))deallocate(DATAP)
!  if(allocated(DATA))deallocate(DATA)
!  call Destroy(this%TBALE)
!  if(allocated(AREA3))deallocate(AREA3)
!  if(allocated(IDATA))deallocate(IDATA)
!  if(allocated(BALE))deallocate(BALE)
!  if(allocated(BALP))deallocate(BALP)
!  if(allocated(CDPTHZ))deallocate(CDPTHZ)
!  if(allocated(DPTHZ))deallocate(DPTHZ)
!  if(allocated(PPI)) deallocate(PPI)
!  if(allocated(PPZ)) deallocate(PPZ)
!  if(allocated(PPX)) deallocate(PPX)
!  if(allocated(PP))deallocate(PP)
!  if(allocated(VOLWM))deallocate(VOLWM)
!  if(allocated(VOLPM))deallocate(VOLPM)

  end subroutine plt_site_destroy
!----------------------------------------------------------------------

  subroutine plt_bgcrate_init(this)
  implicit none
  class(plant_bgcrate_type) :: this


  allocate(this%TCO2S(JZ1))
  allocate(this%TCO2P(JZ1))
  allocate(this%TUPN3B(JZ1))
  allocate(this%TUPH1B(JZ1))
  allocate(this%TUPN3S(JZ1))
  allocate(this%TUPHGS(JZ1))
  allocate(this%TUPOXP(JZ1))
  allocate(this%TUPOXS(JZ1))
  allocate(this%TUPN2S(JZ1))
  allocate(this%TUPCHS(JZ1))
  allocate(this%TUPNH4(JZ1))
  allocate(this%TUPNO3(JZ1))
  allocate(this%TUPH2P(JZ1))
  allocate(this%TUPH2B(JZ1))
  allocate(this%TUPH1P(JZ1))

  allocate(this%TUPNHB(JZ1))
  allocate(this%TUPNOB(JZ1))

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
  allocate(this%ESNT(npelms,jsken,n_pltlitrk,0:JZ1))
  allocate(this%CARBN(JP1))
  allocate(this%XOQCS(1:jcplx,0:JZ1))
  allocate(this%XOQNS(1:jcplx,0:JZ1))
  allocate(this%XOQPS(1:jcplx,0:JZ1))
  allocate(this%CNET(JP1))
  allocate(this%RFGas_root(idg_beg:idg_end-1,JP1))
  allocate(this%TCO2T(JP1))
  allocate(this%TZUPFX(JP1))
  allocate(this%TNH3C(JP1))
  allocate(this%TESN0(npelms,JP1))
  allocate(this%ZESNC(npelms))
  allocate(this%ZNPP(JP1))
  allocate(this%RNH3C(JP1))
  allocate(this%TDFOME(npelms,1:jcplx,JZ1))
  allocate(this%RUPNF(JZ1,JP1))
  allocate(this%TCO2A(JP1))
  allocate(this%RP1BX(0:JZ1))
  allocate(this%RNO3X(0:JZ1))
  allocate(this%RPO4X(0:JZ1))
  allocate(this%RNH4X(0:JZ1))
  allocate(this%ROXYX(0:JZ1))
  allocate(this%RPOBX(0:JZ1))
  allocate(this%RN3BX(0:JZ1))
  allocate(this%RNHBX(0:JZ1))
  allocate(this%RP14X(0:JZ1))

  allocate(this%HESNC(npelms,JP1))
  allocate(this%TESNC(npelms,JP1))
  allocate(this%ESNC(npelms,jsken,1:n_pltlitrk,0:JZ1,JP1))


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

  allocate(this%THIN(JP1))
  allocate(this%XHVSTE(npelms))
  allocate(this%THVSTE(npelms,JP1))
  allocate(this%HVSTE(npelms,JP1))
  allocate(this%VCH4F(JP1))
  allocate(this%VCO2F(JP1))
  allocate(this%VN2OF(JP1))
  allocate(this%VNH3F(JP1))
  allocate(this%VPO4F(JP1))

  allocate(this%IDAYY(JP1))
  allocate(this%VOXYF(JP1))
  allocate(this%EHVST(1:2,1:4,JP1))
  allocate(this%HVST(JP1))
  allocate(this%IYRH(JP1))
  allocate(this%FERT(1:20))
  allocate(this%IYR0(JP1))
  allocate(this%IDAYX(JP1))
  allocate(this%IYRX(JP1))
  allocate(this%IYRY(JP1))
  allocate(this%IDAY0(JP1))
  allocate(this%IDAYH(JP1))
  allocate(this%IHVST(JP1))
  allocate(this%JHVST(JP1))

  end subroutine plt_disturb_init
!----------------------------------------------------------------------
  subroutine plt_disturb_destroy(this)
  implicit none
  class(plant_disturb_type) :: this


!  if(allocated(THIN))deallocate(THIN)
!  if(allocated(HVSTE))deallocate(HVSTE)
!  if(allocated(THVSTE))deallocate(THVSTE)
!  if(allocated(VCH4F))deallocate(VCH4F)
!  if(allocated(VCO2F))deallocate(VCO2F)
!  if(allocated(VN2OF))deallocate(VN2OF)
!  if(allocated(VNH3F))deallocate(VNH3F)
!  if(allocated(VPO4F))deallocate(VPO4F)

!  if(allocated(IDAYY))deallocate(IDAYY)
!  if(allocated(VOXYF))deallocate(VOXYF)
!  if(allocated(EHVST))deallocate(EHVST)
!  if(allocated(HVST))deallocate(HVST)
!  if(allocated(IYRH))deallocate(IYRH)
!  if(allocated(IDAY0))deallocate(IDAY0)
!  if(allocated(IDAYH))deallocate(IDAYH)
!  if(allocated(IYRY))deallocate(IYRY)
!  if(allocated(IDAYX))deallocate(IDAYX)
!  if(allocated(IYR0))deallocate(IYR0)
!  if(allocated(FERT))deallocate(FERT)
!  if(allocated(IYRX))deallocate(IYRX)
!  if(allocated(IHVST))deallocate(IHVST)
!  if(allocated(JHVST))deallocate(JHVST)

  end subroutine plt_disturb_destroy
!----------------------------------------------------------------------
  subroutine  plt_ew_init(this)

  implicit none
  class(plant_ew_type) :: this

  allocate(this%CTRAN(JP1))
  allocate(this%TotalSoilH2OPSIMPa(0:JZ1))
  allocate(this%TUPHT(0:JZ1))
  allocate(this%TKCZ(JP1))
  allocate(this%SFLXC(JP1))
  allocate(this%FLWC(JP1))
  allocate(this%PSILG(JP1))
  allocate(this%PSILT(JP1))
  allocate(this%EVAPC(JP1))
  allocate(this%HFLXC(JP1))
  allocate(this%EFLXC(JP1))
  allocate(this%VOLWC(JP1))
  allocate(this%VHCPC(JP1))
  allocate(this%VOLWP(JP1))
  allocate(this%PSIRT(jroots,JZ1,JP1))
  allocate(this%PSIRO(jroots,JZ1,JP1))
  allocate(this%PSIRG(jroots,JZ1,JP1))
  allocate(this%PopPlantRootH2OUptake_vr(jroots,JZ1,JP1))
  allocate(this%GridPlantRootH2OUptake_vr(0:JZ1))
  allocate(this%EP(JP1))
  allocate(this%PSILO(JP1))
  allocate(this%TKS(0:JZ1))
  allocate(this%OSMO(JP1))
  allocate(this%RAZ(JP1))
  allocate(this%DTKC(JP1))
  allocate(this%TKC(JP1))
  allocate(this%ENGYX(JP1))
  allocate(this%TCC(JP1))
  allocate(this%PSILZ(JP1))

  end subroutine plt_ew_init
!----------------------------------------------------------------------

  subroutine plt_ew_destroy(this)
  implicit none
  class(plant_ew_type) :: this

!  if(allocated(CTRAN))deallocate(CTRAN)
!  if(allocated(PSIST))deallocate(PSIST)
!  if(allocated(TUPHT))deallocate(TUPHT)
!  if(allocated(TKCZ))deallocate(TKCZ)
!  if(allocated(SFLXC))deallocate(SFLXC)
!  if(allocated(FLWC))deallocate(FLWC)
!  if(allocated(PSILG))deallocate(PSILG)
!  if(allocated(PSILT))deallocate(PSILT)
!  if(allocated(EVAPC))deallocate(EVAPC)
!  if(allocated(HFLXC))deallocate(HFLXC)
!  if(allocated(EFLXC))deallocate(EFLXC)
!  if(allocated(VHCPC))deallocate(VHCPC)
!  if(allocated(VOLWC))deallocate(VOLWC)
!  if(allocated(VOLWP))deallocate(VOLWP)
!  if(allocated(PSIRT))deallocate(PSIRT)
!  if(allocated(PSIRO))deallocate(PSIRO)
!  if(allocated(PSIRG))deallocate(PSIRG)
!  if(allocated(PopPlantRootH2OUptake_vr))deallocate(PopPlantRootH2OUptake_vr)
!  if(allocated(TPopPlantRootH2OUptake_vr))deallocate(TPopPlantRootH2OUptake_vr)
!  if(allocated(EP))deallocate(EP)
!  if(allocated(PSILO))deallocate(PSILO)
!  if(allocated(TKS))deallocate(TKS)
!  if(allocated(PSILZ))deallocate(PSILZ)
!  if(allocated(RAZ))deallocate(RAZ)
!  if(allocated(OSMO))deallocate(OSMO)
!  if(allocated(DTKC))deallocate(DTKC)
!  if(allocated(TKC))deallocate(TKC)
!  if(allocated(ENGYX))deallocate(ENGYX)
!  if(allocated(TCC))deallocate(TCC)

  end subroutine plt_ew_destroy

!----------------------------------------------------------------------

  subroutine plt_allom_init(this)
  implicit none
  class(plant_allometry_type) :: this


  allocate(this%CWSRT(JP1))
  allocate(this%FNOD(JP1))
  allocate(this%GRWTB(JBR,JP1))
  allocate(this%DMND(JP1))
  allocate(this%DMRT(JP1))
  allocate(this%CPRT(JP1))
  allocate(this%CNWS(JP1))
  allocate(this%CPWS(JP1))
  allocate(this%CPRTS(JP1))
  allocate(this%CNRTS(JP1))
  allocate(this%CNND(JP1))
  allocate(this%CPND(JP1))
  allocate(this%CNRT(JP1))
  allocate(this%CPLF(JP1))
  allocate(this%CPSHE(JP1))
  allocate(this%CNLF(JP1))
  allocate(this%CNSHE(JP1))
  allocate(this%CNGR(JP1))
  allocate(this%CPSTK(JP1))
  allocate(this%CNSTK(JP1))
  allocate(this%CPGR(JP1))
  allocate(this%CPEAR(JP1))
  allocate(this%CPRSV(JP1))
  allocate(this%CNRSV(JP1))
  allocate(this%CPHSK(JP1))
  allocate(this%FVRN(0:5))
  allocate(this%FWOODE(npelms,n_pltlitrk))
  allocate(this%FWODLE(npelms,n_pltlitrk))
  allocate(this%FWODRE(npelms,n_pltlitrk))
  allocate(this%FWODBE(npelms,n_pltlitrk))

  allocate(this%DMSHE(JP1))
  allocate(this%DMHSK(JP1))
  allocate(this%DMSTK(JP1))
  allocate(this%DMRSV(JP1))
  allocate(this%DMEAR(JP1))
  allocate(this%DMGR(JP1))
  allocate(this%CNHSK(JP1))
  allocate(this%CNEAR(JP1))
  allocate(this%DMLF(JP1))

  end subroutine plt_allom_init
!----------------------------------------------------------------------

  subroutine plt_allom_destroy(this)
  implicit none

  class(plant_allometry_type) :: this

!  if(allocated(FNOD))deallocate(FNOD)
!  if(allocated(GRWTB))deallocate(GRWTB)
!  if(allocated(DMND))deallocate(DMND)
!  if(allocated(DMRT))deallocate(DMRT)
!  if(allocated(CPRT))deallocate(CPRT)
!  if(allocated(CNWS))deallocate(CNWS)
!  if(allocated(CPWS))deallocate(CPWS)
!  if(allocated(CPRTS))deallocate(CPRTS)
!  if(allocated(CNRTS))deallocate(CNRTS)
!  if(allocated(CNND))deallocate(CNND)
!  if(allocated(CPND))deallocate(CPND)
!  if(allocated(CNRT))deallocate(CNRT)
!  if(allocated(CWSRT))deallocate(CWSRT)
!  if(allocated(CPLF))deallocate(CPLF)
!  if(allocated(CPSHE))deallocate(CPSHE)
!  if(allocated(CNSHE))deallocate(CNSHE)
!  if(allocated(CNLF))deallocate(CNLF)
!  if(allocated(DMLF))deallocate(DMLF)
!  if(allocated(CPSTK))deallocate(CPSTK)
!  if(allocated(CNGR))deallocate(CNGR)
!  if(allocated(CNSTK))deallocate(CNSTK)
!  if(allocated(CPGR))deallocate(CPGR)
!  if(allocated(DMSHE))deallocate(DMSHE)
!  if(allocated(DMSTK))deallocate(DMSTK)
!  if(allocated(DMGR))deallocate(DMGR)
!  if(allocated(DMRSV))deallocate(DMRSV)
!  if(allocated(DMEAR))deallocate(DMEAR)
!  if(allocated(DMHSK))deallocate(DMHSK)
!  if(allocated(FVRN))deallocate(FVRN)
!  if(allocated(FWODLE))deallocate(FWODLE)
!  if(allocated(FWOODE))deallocate(FWOODE)
!  if(allocated(FWODRE))deallocate(FWODRE)
!  if(allocated(FWODBE))deallocate(FWODBE)
!  if(allocated(CPEAR))deallocate(CPEAR)
!  if(allocated(CPHSK))deallocate(CPHSK)
!  if(allocated(CNHSK))deallocate(CNHSK)
!  if(allocated(CNEAR))deallocate(CNEAR)
!  if(allocated(CNRSV))deallocate(CNRSV)
!  if(allocated(CPRSV))deallocate(CPRSV)

  end subroutine plt_allom_destroy
!----------------------------------------------------------------------
  subroutine plt_biom_init(this)
  implicit none
  class(plant_biom_type) :: this

  allocate(this%ZEROL(JP1))
  allocate(this%ZEROP(JP1))
  allocate(this%WTSTGET(npelms))
  allocate(this%WTNDLE(npelms,JZ1,JP1))
  allocate(this%WGLFV(JC1,JP1))
  allocate(this%EPOOLN(npelms,JZ1,JP1))
  allocate(this%WTSTDE(npelms,jsken,JP1))
  allocate(this%WTRT2E(npelms,jroots,JZ1,JC1,JP1))
  allocate(this%WTRT1E(npelms,jroots,JZ1,JC1,JP1))
  allocate(this%CEPOLP(npelms,JP1))
  allocate(this%EPOOLP(npelms,JP1))
  allocate(this%EPOLNP(npelms,JP1))
  allocate(this%CCPLNP(JP1))
  allocate(this%CWSRTL(jroots,JZ1,JP1))
  allocate(this%WSRTL(jroots,JZ1,JP1))
  allocate(this%WTRTL(jroots,JZ1,JP1))
  allocate(this%WTRTD(jroots,JZ1,JP1))
  allocate(this%EPOOLR(npelms,jroots,JZ1,JP1))
  allocate(this%CEPOLR(npelms,jroots,JZ1,JP1))
  allocate(this%EPOOL(npelms,JBR,JP1))
  allocate(this%WVSTK(JP1))
  allocate(this%WSLF(0:JNODS1,JBR,JP1))
  allocate(this%WSSHE(0:JNODS1,JBR,JP1))
  allocate(this%WGNODE(npelms,0:JNODS1,JBR,JP1))
  allocate(this%WGLFE(npelms,0:JNODS1,JBR,JP1))
  allocate(this%WGSHE(npelms,0:JNODS1,JBR,JP1))
  allocate(this%WGLFLE(npelms,JC1,0:JNODS1,JBR,JP1))
  allocate(this%WVSTKB(JBR,JP1))
  allocate(this%EPOLNB(npelms,JBR,JP1))
  allocate(this%CEPOLB(npelms,JBR,JP1))
  allocate(this%WTRTSE(npelms,JP1))
  allocate(this%WGLFT(JC1))
  allocate(this%WTRTE(npelms,JP1))
  allocate(this%WTRVX(JP1))
  allocate(this%WTRVE(npelms,JP1))
  allocate(this%WTLS(JP1))
  allocate(this%WTSTGE(npelms,JP1))
  allocate(this%WTRTA(JP1))
  allocate(this%WTSHEE(npelms,JP1))
  allocate(this%WTSTKE(npelms,JP1))
  allocate(this%WTRSVE(npelms,JP1))
  allocate(this%WTGRE(npelms,JP1))
  allocate(this%WTHSKE(npelms,JP1))
  allocate(this%WTEARE(npelms,JP1))
  allocate(this%WTNDE(npelms,JP1))
  allocate(this%WGLFEX(npelms,JBR,JP1))
  allocate(this%WGSHEXE(npelms,JBR,JP1))
  allocate(this%WTSHTE(npelms,JP1))
  allocate(this%WTSHTA(JP1))
  allocate(this%WTLFE(npelms,JP1))
  allocate(this%WTSTDI(JP1))
  allocate(this%WTLSB(JBR,JP1))
  allocate(this%WTRSVBE(npelms,JBR,JP1))
  allocate(this%WTLFBE(npelms,JBR,JP1))
  allocate(this%WTNDBE(npelms,JBR,JP1))
  allocate(this%WTSHEBE(npelms,JBR,JP1))
  allocate(this%WTEARBE(npelms,JBR,JP1))
  allocate(this%WTHSKBE(npelms,JBR,JP1))
  allocate(this%WTGRBE(npelms,JBR,JP1))
  allocate(this%WTSTKBE(npelms,JBR,JP1))
  allocate(this%WTSHTBE(npelms,JBR,JP1))
  allocate(this%WTSTXBE(npelms,JBR,JP1))
  allocate(this%RTWT1E(npelms,jroots,JRS,JP1))

  end subroutine plt_biom_init
!----------------------------------------------------------------------

  subroutine plt_biom_destroy(this)

  implicit none
  class(plant_biom_type) :: this

!  if(allocated(ZEROL))deallocate(ZEROL)
!  if(allocated(ZEROP))deallocate(ZEROP)
!  if(allocated(EPOOLN))deallocate(EPOOLN)
!  if(allocated(WTNDLE))deallocate(WTNDLE)
!  if(allocated(WGLFV))deallocate(WGLFV)
!  call destroy(RTWT1E)
!  if(allocated(WTRT1E))deallocate(WTRT1E)
!  if(allocated(WTSTDE))deallocate(WTSTDE)
!  if(allocated(WTRT2E))deallocate(WTRT2E)
!  if(allocated(CEPOLP))deallocate(CEPOLP)
!  if(allocated(EPOOLP))deallocate(EPOOLP)
!  if(allocated(EPOLNP))deallocate(EPOLNP)
!  if(allocated(CCPLNP))deallocate(CCPLNP)
!  if(allocated(CWSRTL))deallocate(CWSRTL)
!  if(allocated(WSRTL))deallocate(WSRTL)
!  if(allocated(WTRTL))deallocate(WTRTL)
!  if(allocated(WTRTD))deallocate(WTRTD)
!  if(allocated(EPOOLR))deallocate(EPOOLR)
!  if(allocated(WVSTK))deallocate(WVSTK)
!  if(allocated(WGLFE))deallocate(WGLFE)
!  if(allocated(WSLF))deallocate(WSLF)
!  if(allocated(WSSHE))deallocate(WSSHE)
!  if(allocated(WGLFLE))deallocate(WGLFLE)
!  if(allocated(WGNODE))deallocate(WGNODE)
!  if(allocated(WGSHE))deallocate(WGSHE)
!  if(allocated(WVSTKB))deallocate(WVSTKB)
!  if(allocated(PPOOL))deallocate(PPOOL)
!  if(allocated(EPOLNB))deallocate(EPOLNB)
!  if(allocated(CEPOLB))deallocate(CEPOLB)
!  if(allocated(WTRTSE))deallocate(WTRTSE)
!  if(allocated(WGLFT))deallocate(WGLFT)
!  if(allocated(EPOOL))deallocate(EPOOL)
!  if(allocated(WTRTA))deallocate(WTRTA)
!  if(allocated(WTLF))deallocate(WTLF)
!  if(allocated(WTSHE))deallocate(WTSHE)
!  if(allocated(WTRSV))deallocate(WTRSV)
!  if(allocated(WTSTK))deallocate(WTSTK)
!  if(allocated(WTLFP))deallocate(WTLFP)
!  if(allocated(WTGR))deallocate(WTGR)
!  if(allocated(WTLFN))deallocate(WTLFN)
!  if(allocated(WTEAR))deallocate(WTEAR)
!  if(allocated(WTHSKE))deallocate(WTHSKE)
!  if(allocated(WGLFEX))deallocate(WGLFEX)
!  if(allocated(WGSHEXE))deallocate(WGSHEXE)
!  if(allocated(WTLSB))deallocate(WTLSB)
!  if(allocated(WTRSVBE))deallocate(WTRSVBE)
!  if(allocated(WTLFBE))deallocate(WTLFBE)
!  if(allocated(WTNDBE))deallocate(WTNDBE)
!  if(allocated(WTSHEBE))deallocate(WTSHEBE)
!  if(allocated(WTEARBE))deallocate(WTEARBE)
!  if(allocated(WTHSKBE))deallocate(WTHSKBE)
!  if(allocated(WTGRBE))deallocate(WTGRBE)
!  if(allocated(WTSTKBE))deallocate(WTSTKBE)
!  if(allocated(WTSHTBE))deallocate(WTSHTBE)
!  if(allocated(WTSTXBE))deallocate(WTSTXBE)
!  if(allocated(WTSTDI))deallocate(WTSTDI)
!  if(allocated(WTRVX))deallocate(WTRVX)
!  if(allocated(WTLS))deallocate(WTLS)
!  if(allocated(WTSHTE))deallocate(WTSHTE)
!  if(allocated(WTSHTA))deallocate(WTSHTA)
!  if(allocated(WTSTGE)deallocate(WTSTGE)
!  if(allocated(WTNDE))deallocate(WTNDE)
  end subroutine plt_biom_destroy


!----------------------------------------------------------------------

  subroutine  plt_soilchem_init(this)

  implicit none

  class(plant_soilchem_type) :: this

  allocate(this%FOSRH(1:jcplx,0:JZ1))
  allocate(this%CFOPE(npelms,0:Jlitgrp,jsken,JP1))
  allocate(this%TFND(0:JZ1))
  allocate(this%THETPM(60,0:JZ1))
  allocate(this%DFGS(60,0:JZ1))
  allocate(this%VOLY(0:JZ1))
  allocate(this%VOLI(0:JZ1))
  allocate(this%VOLW(0:JZ1))
  allocate(this%VOLA(0:JZ1))

  allocate(this%trcs_VLN(ids_nuts_beg:ids_nuts_end,0:JZ1))
  allocate(this%OQC(1:jcplx,0:JZ1))
  allocate(this%OQN(1:jcplx,0:JZ1))
  allocate(this%OQP(1:jcplx,0:JZ1))

  allocate(this%trc_solml(ids_beg:ids_end,0:JZ1))

  allocate(this%trc_gasml(idg_beg:idg_end,0:JZ1))
  allocate(this%trc_gascl(idg_beg:idg_end,0:JZ1))
  allocate(this%CORGC(0:JZ1))

  allocate(this%trc_solcl(ids_beg:ids_end,0:jZ1))

  allocate(this%VOLX(0:JZ1))
  allocate(this%THETW(0:JZ1))
  allocate(this%THETY(0:JZ1))

  allocate(this%GSolbility(idg_beg:idg_end,0:JZ1))
  allocate(this%GasDifc(idg_beg:idg_end,0:JZ1))
  allocate(this%SolDifc(ids_beg:ids_end,0:JZ1))
  allocate(this%RSCS(JZ1))
  allocate(this%BKDS(0:JZ1))
  allocate(this%CNDU(JZ1))

  end subroutine plt_soilchem_init
!----------------------------------------------------------------------

  subroutine plt_soilchem_destroy(this)
  implicit none
  class(plant_soilchem_type) :: this

!  if(allocated(FOSRH))deallocate(FOSRH)

!  if(allocated(CFOPE))deallocate(CFOPE)
!  if(allocated(TFND))deallocate(TFND)
!  if(allocated(THETPM))deallocate(THETPM)
!  if(allocated(DFGS))deallocate(DFGS)
!  if(allocated(ZVSGL))deallocate(ZVSGL)
!  if(allocated(SOXYL))deallocate(SOXYL)

!  if(allocated(CNDU))deallocate(CNDU)
!  if(allocated(CGSGL))deallocate(CGSGL)
!  if(allocated(CHSGL))deallocate(CHSGL)
!  if(allocated(HGSGL))deallocate(HGSGL)
!  if(allocated(OGSGL))deallocate(OGSGL)
!  if(allocated(RSCS))deallocate(RSCS)

!   call destroy(this%GSolbility)
!  if(allocated(THETW))deallocate(THETW)
!  if(allocated(THETY))deallocate(THETY)
!  if(allocated(VOLX))deallocate(VOLX)
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
!  if(allocated(VOLY))deallocate(VOLY)
!  if(allocated(VOLI))deallocate(VOLI)
!  if(allocated(VOLW))deallocate(VOLW)
!  if(allocated(VOLA))deallocate(VOLA)
!  if(allocated(ZOSGL))deallocate(ZOSGL)

!  if(allocated(OQC))deallocate(OQC)
!  if(allocated(OQN))deallocate(OQN)
!  if(allocated(OQP))deallocate(OQP)

!  if(allocated(ZNSGL))deallocate(ZNSGL)

!  if(allocated(Z2SGL))deallocate(Z2SGL)
!  if(allocated(ZHSGL))deallocate(ZHSGL)
!  if(allocated(BKDS))deallocate(BKDS)
!  if(allocated(CO2G))deallocate(CO2G)
!  if(allocated(CLSGL))deallocate(CLSGL)
!  if(allocated(CQSGL))deallocate(CQSGL)

  end subroutine plt_soilchem_destroy
!----------------------------------------------------------------------


  subroutine InitPlantAPIData()

  implicit none

  JZ1    => pltpar%JZ1
  JC1    => pltpar%JC1
  JP1    => pltpar%JP1
  JLA1   => pltpar%JLA1
  JSA1   => pltpar%JSA1
  JLI1   => pltpar%JLI1
  JNODS1 => pltpar%JNODS1
  !the following variable should be consistent with the soil bgc model
  jcplx => pltpar%jcplx
  jsken  => pltpar%jsken
  Jlitgrp=> pltpar%Jlitgrp
  JBR    => pltpar%JBR
  JRS    => pltpar%JRS
  JPRT   => pltpar%JPRT
  n_pltlitrk => pltpar%n_pltlitrk
  jpstgs => pltpar%jpstgs
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

  allocate(this%PAR(JLI1,JSA1,JC1,JP1))
  allocate(this%PARDIF(JLI1,JSA1,JC1,JP1))
  allocate(this%ALBR(JP1))
  allocate(this%ALBP(JP1))
  allocate(this%TAUS(JC1+1))
  allocate(this%TAU0(JC1+1))
  allocate(this%THRM1(JP1))
  allocate(this%RADC(JP1))
  allocate(this%OMEGX(JSA1,JLI1,JLA1))
  allocate(this%OMEGAG(JSA1))
  allocate(this%OMEGA(JSA1,JLI1,JLA1))
  allocate(this%ZSIN(JLI1))
  allocate(this%ZCOS(JLI1))
  allocate(this%IALBY(JSA1,JLI1,JLA1))
  allocate(this%RAD1(JP1))
  allocate(this%ABSR(JP1))
  allocate(this%ABSP(JP1))
  allocate(this%TAUP(JP1))
  allocate(this%TAUR(JP1))
  allocate(this%RADP(JP1))
  allocate(this%FRADP(JP1))
  end subroutine plt_rad_init
!------------------------------------------------------------------------
  subroutine plt_rad_destroy(this)
! DESCRIPTION
! deallocate memory for plant_radiation_type
  implicit none
  class(plant_radiation_type) :: this

!  if(allocated(PAR))deallocate(PAR)
!  if(allocated(PARDIF))deallocate(PARDIF)
!  if(associated(this%ALBR))deallocate(this%ALBR)
!  if(associated(this%ALBP))deallocate(this%ALBP)
!  if(associated(this%TAUS))deallocate(this%TAUS)
!  if(associated(this%TAU0))deallocate(this%TAU0)
!  if(associated(this%THRM1))deallocate(this%THRM1)
!  if(associated(this%RADC))deallocate(this%RADC)
!  if(associated(this%OMEGX))deallocate(this%OMEGX)
!  if(associated(this%OMEGAG))deallocate(this%OMEGAG)
!  if(associated(this%ZSIN))deallocate(this%ZSIN)
!  if(allocated(ZCOS))deallocate(ZCOS)
!  if(associated(this%IALBY))deallocate(this%IALBY)
!  if(associated(this%OMEGA))deallocate(this%OMEGA)
!  if(associated(this%RAD1))deallocate(this%RAD1)
!  if(allocated(ABSR))deallocate(ABSR)
!  if(allocated(ABSP))deallocate(ABSP)
!  if(allocated(TAUR))deallocate(TAUR)
!  if(allocated(TAUP))deallocate(TAUP)
!  if(allocated(FRADP))deallocate(FRADP)
!  if(allocated(RADP))deallocate(RADP)

  end subroutine plt_rad_destroy

!------------------------------------------------------------------------

  subroutine plt_photo_init(this)
  class(plant_photosyns_type) :: this

  allocate(this%RSMX(JP1))
  allocate(this%RSMN(JP1))
  allocate(this%SO2(JP1))
  allocate(this%RA(JP1))
  allocate(this%RC(JP1))
  allocate(this%DCO2(JP1))
  allocate(this%FMOL(JP1))
  allocate(this%RCMX(JP1))
  allocate(this%SURFX(JLI1,JC1,JNODS1,JBR,JP1))
  allocate(this%CPOOL3(JNODS1,JBR,JP1))
  allocate(this%CPOOL4(JNODS1,JBR,JP1))
  allocate(this%CO2B(JNODS1,JBR,JP1))
  allocate(this%COMPL(JNODS1,JBR,JP1))
  allocate(this%CBXN(JNODS1,JBR,JP1))
  allocate(this%CBXN4(JNODS1,JBR,JP1))
  allocate(this%ETGRO(JNODS1,JBR,JP1))
  allocate(this%ETGR4(JNODS1,JBR,JP1))
  allocate(this%FDBK4(JNODS1,JBR,JP1))
  allocate(this%HCOB(JNODS1,JBR,JP1))

  allocate(this%VCGRO(JNODS1,JBR,JP1))
  allocate(this%VGRO(JNODS1,JBR,JP1))
  allocate(this%VCGR4(JNODS1,JBR,JP1))
  allocate(this%VGRO4(JNODS1,JBR,JP1))
  allocate(this%ICTYP(JP1))
  allocate(this%XKCO24(JP1))
  allocate(this%O2L(JP1))
  allocate(this%SCO2(JP1))
  allocate(this%CO2Q(JP1))
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
  allocate(this%FDBK(JBR,JP1))
  allocate(this%FDBKX(JBR,JP1))
  allocate(this%CO2L(JP1))
  allocate(this%XKCO2L(JP1))
  allocate(this%XKCO2O(JP1))
  allocate(this%CO2I(JP1))
  allocate(this%O2I(JP1))
  allocate(this%RCS(JP1))
  allocate(this%FCO2(JP1))
  allocate(this%RSMH(JP1))

  end subroutine plt_photo_init
!------------------------------------------------------------------------

  subroutine plt_photo_destroy(this)
  class(plant_photosyns_type) :: this

!  if(allocated(RSMX))deallocate(RSMX)
!  if(allocated(RSMN))deallocate(RSMN)
!  if(allocated(SO2))deallocate(SO2)
!  if(allocated(RC))deallocate(RC)
!  if(allocated(RA))deallocate(RA)
!  if(allocated(DCO2))deallocate(DCO2)
!  if(allocated(FMOL))deallocate(FMOL)
!  if(allocated(RCMX))deallocate(RCMX)
!  if(allocated(SURFX))deallocate(SURFX)
!  if(allocated(CPOOL3))deallocate(CPOOL3)
!  if(allocated(CPOOL4))deallocate(CPOOL4)
!  if(allocated(CO2B))deallocate(CO2B)
!  if(allocated(COMPL))deallocate(COMPL)
!  if(allocated(CBXN))deallocate(CBXN)
!  if(allocated(CBXN4))deallocate(CBXN4)
!  if(allocated(ETGRO))deallocate(ETGRO)
!  if(allocated(ETGR4))deallocate(ETGR4)
!  if(allocated(FDBK4))deallocate(FDBK4)
!  if(allocated(HCOB))deallocate(HCOB)

!  if(allocated(VCGRO))deallocate(VCGRO)
!  if(allocated(VGRO))deallocate(VGRO)
!  if(allocated(VCGR4))deallocate(VCGR4)
!  if(allocated(VGRO4))deallocate(VGRO4)
!  if(allocated(ICTYP))deallocate(ICTYP)
!  if(allocated(XKCO24))deallocate(XKCO24)
!  if(allocated(O2L))deallocate(O2L)
!  if(allocated(SCO2))deallocate(SCO2)
!  if(allocated(CO2Q))deallocate(CO2Q)
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
!  if(allocated(FDBK))deallocate(FDBK)
!  if(allocated(FDBKX))deallocate(FDBKX)
!  if(allocated(CO2L))deallocate(CO2L)
!  if(allocated(XKCO2L))deallocate(XKCO2L)
!  if(allocated(XKCO2O))deallocate(XKCO2O)
!  if(allocated(CO2I))deallocate(CO2I)
!  if(allocated(O2I))deallocate(O2I)
!  if(allocated(RCS))deallocate(RCS)
!  if(allocated(FCO2))deallocate(FCO2)
!  if(allocated(RSMH)) deallocate(RSMH)
  end subroutine plt_photo_destroy

!------------------------------------------------------------------------
  subroutine plt_pheno_init(this)
  implicit none
  class(plant_pheno_type) :: this


  allocate(this%CTC(JP1))
  allocate(this%OFFST(JP1))
  allocate(this%GROUPI(JP1))
  allocate(this%PlantO2Stress(JP1))
  allocate(this%PB(JP1))
  allocate(this%PR(JP1))
  allocate(this%TFN3(JP1))
  allocate(this%TCZ(JP1))
  allocate(this%TCG(JP1))
  allocate(this%TKG(JP1))
  allocate(this%TCX(JP1))
  allocate(this%WSTR(JP1))
  allocate(this%RCELX(npelms,JBR,JP1))
  allocate(this%RCESX(npelms,JBR,JP1))

  allocate(this%TFN4(JZ1,JP1))
  allocate(this%GFILL(JP1))
  allocate(this%PTSHT(JP1))
  allocate(this%SSTX(JP1))
  allocate(this%HTC(JP1))
  allocate(this%ZTYPI(JP1))
  allocate(this%ZTYP(JP1))
  allocate(this%IFLGC(JP1))
  allocate(this%IDTH(JP1))
  allocate(this%RSETE(npelms,JP1))
  allocate(this%GROUP(JBR,JP1))
  allocate(this%IDTHP(JP1))
  allocate(this%ATRP(JBR,JP1))
  allocate(this%FLGZ(JBR,JP1))
  allocate(this%DGSTGI(JBR,JP1))
  allocate(this%DGSTGF(JBR,JP1))
  allocate(this%GSTGI(JBR,JP1))
  allocate(this%GSTGF(JBR,JP1))
  allocate(this%FLG4(JBR,JP1))
  allocate(this%IWTYP(JP1))
  allocate(this%ISTYP(JP1))
  allocate(this%IBTYP(JP1))
  allocate(this%IDTHR(JP1))
  allocate(this%IDTYP(JP1))
  allocate(this%IPTYP(JP1))
  allocate(this%IFLGI(JP1))
  allocate(this%IGTYP(JP1))
  allocate(this%KVSTG(JBR,JP1))
  allocate(this%KVSTGN(JBR,JP1))
  allocate(this%IDAY(jpstgs,JBR,JP1))
  allocate(this%TGSTGI(JBR,JP1))
  allocate(this%TGSTGF(JBR,JP1))
  allocate(this%VSTGX(JBR,JP1))
  allocate(this%XRLA(JP1))
  allocate(this%XRNI(JP1))
  allocate(this%XDL(JP1))
  allocate(this%XPPD(JP1))
  allocate(this%IDTHB(JBR,JP1))
  allocate(this%IFLGP(JBR,JP1))
  allocate(this%IFLGF(JBR,JP1))
  allocate(this%IFLGE(JBR,JP1))
  allocate(this%IFLGA(JBR,JP1))
  allocate(this%IFLGG(JBR,JP1))
  allocate(this%IFLGR(JBR,JP1))
  allocate(this%IFLGQ(JBR,JP1))
  allocate(this%VRNY(JBR,JP1))
  allocate(this%VRNZ(JBR,JP1))
  allocate(this%VRNS(JBR,JP1))
  allocate(this%VRNL(JC1,JP1))
  allocate(this%VRNF(JBR,JP1))
  allocate(this%VRNX(JC1,JP1))


  end subroutine plt_pheno_init
!------------------------------------------------------------------------

  subroutine plt_pheno_destroy(this)
  implicit none
  class(plant_pheno_type) :: this

!  if(allocated(CTC))deallocate(CTC)
!  if(allocated(OFFST))deallocate(OFFST)
!  if(allocated(GROUPI))deallocate(GROUPI)
!  if(allocated(PlantO2Stress))deallocate(PlantO2Stress)
!  if(allocated(PB))deallocate(PB)
!  if(allocated(PR))deallocate(PR)
!  if(allocated(TFN3))deallocate(TFN3)
!  if(allocated(TCZ))deallocate(TCZ)
!  if(allocated(TCG))deallocate(TCG)
!  if(allocated(TKG))deallocate(TKG)
!  if(allocated(TCX))deallocate(TCX)
!  if(allocated(WSTR))deallocate(WSTR)
!  if(allocated(RCELX))deallocate(RCELX)
!  if(allocated(RCESX))deallocate(RCESX)
!  if(allocated(TFN4))deallocate(TFN4)
!  if(allocated(GFILL))deallocate(GFILL)
!  if(allocated(RTARP))deallocate(RTARP)
!  if(allocated(PTSHT))deallocate(PTSHT)
!  if(allocated(SSTX)) deallocate(SSTX)
!  if(allocated(HTC))deallocate(HTC)
!  if(allocated(ZTYPI))deallocate(ZTYPI)
!  if(allocated(ZTYP))deallocate(ZTYP)
!  if(allocated(IFLGC))deallocate(IFLGC)
!  if(allocated(IDTH))deallocate(IDTH)
!  if(allocated(RSETE))deallocate(RSETE)
!  if(allocated(GROUP))deallocate(GROUP)
!  if(allocated(IDTHP))deallocate(IDTHP)
!  if(allocated(ATRP))deallocate(ATRP)
!  if(allocated(FLGZ))deallocate(FLGZ)
!  if(allocated(DGSTGI))deallocate(DGSTGI)
!  if(allocated(DGSTGF))deallocate(DGSTGF)
!  if(allocated(GSTGI))deallocate(GSTGI)
!  if(allocated(GSTGF))deallocate(GSTGF)
!  if(allocated(FLG4))deallocate(FLG4)
!  if(allocated(IWTYP))deallocate(IWTYP)
!  if(allocated(IBTYP))deallocate(IBTYP)
!  if(allocated(ISTYP))deallocate(ISTYP)
!  if(allocated(IDTHR))deallocate(IDTHR)
!  if(allocated(IDTYP))deallocate(IDTYP)
!  if(allocated(IPTYP))deallocate(IPTYP)
!  if(allocated(IGTYP))deallocate(IGTYP)
!  if(allocated(IFLGI))deallocate(IFLGI)
!  if(allocated(KVSTGN))deallocate(KVSTGN)
!  if(allocated(TGSTGI))deallocate(TGSTGI)
!  if(allocated(TGSTGF))deallocate(TGSTGF)
!  if(allocated(VSTGX))deallocate(VSTGX)
!  if(allocated(KVSTG))deallocate(KVSTG)
!  if(allocated(XDL))deallocate(XDL)
!  if(allocated(XRLA))deallocate(XRLA)
!  if(allocated(XRNI))deallocate(XRNI)
!  if(allocated(XPPD))deallocate(XPPD)
!  if(allocated(IDTHB))deallocate(IDTHB)
!  if(allocated(IFLGP))deallocate(IFLGP)
!  if(allocated(IFLGF))deallocate(IFLGF)
!  if(allocated(IFLGE))deallocate(IFLGE)
!  if(allocated(IFLGA))deallocate(IFLGA)
!  if(allocated(IFLGG))deallocate(IFLGG)
!  if(allocated(IFLGR))deallocate(IFLGR)
!  if(allocated(IFLGQ))deallocate(IFLGQ)
!  if(allocated(VRNY))deallocate(VRNY)
!  if(allocated(VRNZ))deallocate(VRNZ)
!  if(allocated(VRNS))deallocate(VRNS)
!  if(allocated(VRNL))deallocate(VRNL)
!  if(allocated(VRNF))deallocate(VRNF)
!  if(allocated(VRNX))deallocate(VRNX)
!  if(allocated(IDAY))deallocate(IDAY)

  end subroutine plt_pheno_destroy
!------------------------------------------------------------------------

  subroutine plt_morph_init(this)
  implicit none
  class(plant_morph_type) :: this

  allocate(this%RTARP(jroots,JZ1,JP1))
  allocate(this%RTDNP(jroots,JZ1,JP1))
  allocate(this%RTVLP(jroots,JZ1,JP1))
  allocate(this%RTVLW(jroots,JZ1,JP1))
  allocate(this%RTN1(jroots,JZ1,JP1))
  allocate(this%RTNL(jroots,JZ1,JP1))
  allocate(this%GRDM(JP1))
  allocate(this%RTDNT(JZ1))
  allocate(this%RTFQ(JP1))
  allocate(this%CFI(JP1))
  allocate(this%CFX(JP1))
  allocate(this%HTCTL(JP1))
  allocate(this%PORTX(jroots,JP1))
  allocate(this%WDLF(JP1))
  allocate(this%SDMX(JP1))
  allocate(this%STMX(JP1))
  allocate(this%GRMX(JP1))
  allocate(this%RRAD1X(jroots,JP1))
  allocate(this%RRAD2X(jroots,JP1))
  allocate(this%RRADP(jroots,JP1))
  allocate(this%RRAD1(jroots,JZ1,JP1))
  allocate(this%RRAD2(jroots,JZ1,JP1))
  allocate(this%RRAD1M(jroots,JP1))
  allocate(this%RRAD2M(jroots,JP1))

  allocate(this%RTDP1(jroots,JC1,JP1))
  allocate(this%RTLGP(jroots,JZ1,JP1))
  allocate(this%RTLGA(jroots,JZ1,JP1))
  allocate(this%RTLG1X(jroots,JP1))
  allocate(this%RTLG2X(jroots,JP1))
  allocate(this%RTLG1(jroots,JZ1,JC1,JP1))
  allocate(this%RTLG2(jroots,JZ1,JC1,JP1))
  allocate(this%RTN2(jroots,JZ1,JC1,JP1))
  allocate(this%INTYP(JP1))
  allocate(this%MY(JP1))
  allocate(this%HTSTZ(JP1))
  allocate(this%KLEAF(JBR,JP1))
  allocate(this%VSTG(JBR,JP1))
  allocate(this%NG(JP1))
  allocate(this%CanopyHeight(JP1))
  allocate(this%XTLI(JP1))
  allocate(this%ZL(0:JC1))
  allocate(this%ARSTP(JP1))
  allocate(this%ARLFP(JP1))
  allocate(this%NB1(JP1))
  allocate(this%NIX(JP1))
  allocate(this%NRT(JP1))
  allocate(this%NNOD(JP1))
  allocate(this%NBT(JP1))
  allocate(this%NBR(JP1))
  allocate(this%NINR(JC1,JP1))
  allocate(this%PSTG(JBR,JP1))
  allocate(this%PSTGI(JBR,JP1))
  allocate(this%PSTGF(JBR,JP1))
  allocate(this%ANGBR(JP1))
  allocate(this%SURF(JLI1,JC1,JNODS1,JBR,JP1))
  allocate(this%KLEAFX(JBR,JP1))
  allocate(this%NBTB(JBR,JP1))
  allocate(this%GRNXB(JBR,JP1))
  allocate(this%HTSHEX(JBR,JP1))
  allocate(this%ARLFZ(JBR,JP1))
  allocate(this%ARLFB(JBR,JP1))
  allocate(this%ANGSH(JP1))
  allocate(this%CLASS(JLI1,JP1))
  allocate(this%ARSTV(JC1,JP1))
  allocate(this%ARLFV(JC1,JP1))
  allocate(this%ARLF1(0:JNODS1,JBR,JP1))
  allocate(this%GRNO(JP1))
  allocate(this%SDPTH(JP1))
  allocate(this%SDPTHI(JP1))
  allocate(this%SDLG(JP1))
  allocate(this%SDVL(JP1))
  allocate(this%SDAR(JP1))
  allocate(this%ARSTT(JC1))
  allocate(this%SSL1(JP1))
  allocate(this%SNL1(JP1))
  allocate(this%SLA1(JP1))
  allocate(this%ARLFT(JC1))
  allocate(this%ARLFS(JP1))
  allocate(this%HTNODX(0:JNODS1,JBR,JP1))
  allocate(this%HTSHE(0:JNODS1,JBR,JP1))
  allocate(this%HTNODE(0:JNODS1,JBR,JP1))
  allocate(this%SURFB(JLI1,JC1,JBR,JP1))
  allocate(this%ARLFL(JC1,0:JNODS1,JBR,JP1))
  allocate(this%ARSTK(JC1,JBR,JP1))
  allocate(this%NI(JP1))
  allocate(this%GRNOB(JBR,JP1))
  allocate(this%CF(JP1))
  allocate(this%DMVL(jroots,JP1))
  allocate(this%PORT(jroots,JP1))
  allocate(this%RTAR2X(jroots,JP1))
  allocate(this%RTAR1X(jroots,JP1))
  allocate(this%RSRR(jroots,JP1))
  allocate(this%RSRA(jroots,JP1))
  allocate(this%IRTYP(JP1))
  end subroutine plt_morph_init

!------------------------------------------------------------------------

  subroutine plt_morph_destroy(this)
  implicit none
  class(plant_morph_type) :: this

  end subroutine plt_morph_destroy
end module PlantAPIData
