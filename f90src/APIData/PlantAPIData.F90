module PlantAPIData
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use ElmIDMod
implicit none
  save
  character(len=*),private, parameter :: mod_filename = __FILE__
  public

! grid configuration
  integer  :: JP1         !number of plants
  integer  :: JSA1        !number of sectors for the sky azimuth  [0,2*pi]
  integer  :: jcplx11     !number of organo-microbial complexes
  integer  :: JLA1        !number of sectors for the leaf azimuth, [0,pi]
  integer  :: JC1         !number of canopy layers
  integer  :: JZ1         !number of soil layers
  integer  :: JLI1        !number of sectors for the leaf zenith [0,pi/2]
  integer  :: JNODS1      !number of canopy nodes
  integer  :: jsken       !number of kinetic components in litter,  PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  integer  :: Jlitgrp     !number of litter groups nonstructural(0,*),
                          !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
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
  real(r8), pointer :: PP(:)      => null()    !plant population, [d-2]
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
  real(r8) :: ZT                        !canopy height , [m]
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
  real(r8), pointer :: ZC(:)           => null() !canopy height, [m]
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
  real(r8), pointer :: OSTR(:)     => null()     !plant O2 stress indicator, []
  real(r8), pointer :: PB(:)       => null()     !branch nonstructural C content required for new branch, [gC gC-1]
  real(r8), pointer :: PR(:)       => null()     !threshold root nonstructural C content for initiating new root axis, [gC gC-1]
  real(r8), pointer :: RCCLX(:,:) => null()     !C translocated from leaf during senescence, [g d-2 h-1]
  real(r8), pointer :: RCZLX(:,:) => null()     !N translocated from leaf during senescence, [g d-2 h-1]
  real(r8), pointer :: RCPLX(:,:) => null()     !P translocated from leaf during senescence, [g d-2 h-1]
  real(r8), pointer :: RCCSX(:,:) => null()     !C translocated from sheath during senescence, [g d-2 h-1]
  real(r8), pointer :: RCZSX(:,:) => null()     !N translocated from sheath during senescence, [g d-2 h-1]
  real(r8), pointer :: RCPSX(:,:) => null()     !P translocated from sheath during senescence, [g d-2 h-1]
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
  real(r8), pointer :: CFOPC(:,:,:)=> null()  !litter kinetic fraction, [-]
  real(r8), pointer :: CFOPN(:,:,:)=> null()  !litterfall kinetic N fraction, [-]
  real(r8), pointer :: CFOPP(:,:,:)=> null()  !litter P kinetic fraction, [-]
  real(r8), pointer :: TFND(:)     => null()  !temperature effect on diffusivity
  real(r8), pointer :: THETPM(:,:) => null()  !soil air-filled porosity, [m3 m-3]
  real(r8), pointer :: DFGS(:,:)   => null()  !coefficient for dissolution - volatilization, []
  real(r8), pointer :: RSCS(:)     => null()  !soil hydraulic resistance, [MPa h m-2]
  real(r8), pointer :: ZHSGL(:)    => null()  !aqueous NH4 diffusivity, [m2 h-1]
  real(r8), pointer :: BKDS(:)     => null()  !soil bulk density, [Mg m-3]
  real(r8), pointer :: CH2P4(:)    => null()  !aqueous PO4 concentration non-band	[gP m-3]
  real(r8), pointer :: CH1P4B(:)   => null()  !aqueous H1PO4 concentration band [gP m-3]
  real(r8), pointer :: CH2P4B(:)   => null()  !aqueous PO4 concentration band	[gP m-3]
  real(r8), pointer :: CNO3S(:)    => null()  !NO3 concentration non-band micropore	[gN m-3]
  real(r8), pointer :: CH1P4(:)    => null()  !aqueous H1PO4 concentration non-band [gP m-3]
  real(r8), pointer :: CNO3B(:)    => null()  !NO3 concentration band micropore	[gN m-3]
  real(r8), pointer :: CNH4S(:)    => null()  !NH4 concentration non-band micropore	[gN m-3]
  real(r8), pointer :: CNH4B(:)    => null()  !NH4 concentration band micropore	[gN m-3]
  real(r8), pointer :: CO2G(:)     => null()  !gaseous CO2	[g d-2]
  real(r8), pointer :: CO2S(:)     => null()  !aqueous CO2  micropore	[gC d-2]
  real(r8), pointer :: CH4S(:)     => null()  !aqueous CH4  micropore	[gC d-2]
  real(r8), pointer :: CCH4S(:)    => null()  !aqueous CH4 concentration micropore	[gC m-3]
  real(r8), pointer :: CZ2OS(:)    => null()  !aqueous N2O concentration micropore	[gN m-3]
  real(r8), pointer :: CCH4G(:)    => null()  !gaseous CH4 concentration	[gC m-3]
  real(r8), pointer :: CNH3S(:)    => null()  !NH3 concentration non-band micropore	[gN m-3]
  real(r8), pointer :: CNH3B(:)    => null()  !NH3 concentration band micropore	[gN m-3]
  real(r8), pointer :: CH2GS(:)    => null()  !aqueous H2 concentration	[g m-3]
  real(r8), pointer :: HLSGL(:)    => null()  !aqueous H2 diffusivity, [m2 h-1]
  real(r8), pointer :: CLSGL(:)    => null()  !aqueous CO2 diffusivity	[m2 h-1]
  real(r8), pointer :: CQSGL(:)    => null()  !aqueous CH4 diffusivity	[m2 h-1]
  real(r8), pointer :: CZ2OG(:)    => null()  !gaseous N2O concentration	[gN m-3]
  real(r8), pointer :: CNH3G(:)    => null()  !gaseous NH3 concentration	[gN m-3]
  real(r8), pointer :: CH2GG(:)    => null()  !gaseous H2 concentration	[g m-3]
  real(r8), pointer :: CORGC(:)    => null()  !soil organic C content [gC kg soil-1]
  real(r8), pointer :: CPO4S(:)    => null()  !PO4 concentration non-band micropore	[g m-3]
  real(r8), pointer :: CNDU(:)     => null()  !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8), pointer :: CGSGL(:)    => null()  !gaseous CO2 diffusivity	[m2 h-1]
  real(r8), pointer :: CHSGL(:)    => null()  !gaseous CH4 diffusivity	[m2 h-1]
  real(r8), pointer :: HGSGL(:)    => null()  !gaseous H2 diffusivity  [m2 h-1]
  real(r8), pointer :: OGSGL(:)    => null()  !gaseous O2 diffusivity	[m2 h-1]
  real(r8), pointer :: Z2SGL(:)    => null()  !gaseous N2O diffusivity, [m2 h-1]
  real(r8), pointer :: ZOSGL(:)    => null()  !aqueous NO3 diffusivity, [m2 h-1]
  real(r8), pointer :: H1PO4(:)    => null()  !soil aqueous HPO4 content micropore non-band, [gP d-2]
  real(r8), pointer :: H2PO4(:)    => null()  !PO4 non-band micropore, [gP d-2]
  real(r8), pointer :: H1POB(:)    => null()  !soil aqueous HPO4 content micropore band, [gP d-2]
  real(r8), pointer :: H2POB(:)    => null()  !PO4 band micropore, [gP d-2]
  real(r8), pointer :: H2GS(:)     => null()  !aqueous H2 	[g d-2]
  real(r8), pointer :: OXYG(:)     => null()  !gaseous O2 	[g d-2]
  real(r8), pointer :: OXYS(:)     => null()  !aqueous O2  micropore	[g d-2]
  real(r8), pointer :: OLSGL(:)    => null()  !aqueous CO2 diffusivity	[m2 h-1]
  real(r8), pointer :: POSGL(:)    => null()  !aqueous PO4 diffusivity, [m2 h-1]
  real(r8), pointer :: SCO2L(:)    => null()  !solubility of CO2, [m3 m-3]
  real(r8), pointer :: SOXYL(:)    => null()  !solubility of O2, [m3 m-3]
  real(r8), pointer :: SCH4L(:)    => null()  !solubility of CH4, [m3 m-3]
  real(r8), pointer :: SN2OL(:)    => null()  !solubility of N2O, [m3 m-3]
  real(r8), pointer :: SNH3L(:)    => null()  !solubility of NH3, [m3 m-3]
  real(r8), pointer :: SH2GL(:)    => null()  !solubility of H2, [m3 m-3]
  real(r8), pointer :: THETW(:)    => null()  !volumetric water content [m3 m-3]
  real(r8), pointer :: THETY(:)    => null()  !air-dry water content, [m3 m-3]
  real(r8), pointer :: VOLX(:)     => null()  !volume of soil layer	m3 d-2
  real(r8), pointer :: VLPOB(:)    => null()  !PO4 band volume fracrion, [0-1]
  real(r8), pointer :: VLNO3(:)    => null()  !NO3 non-band volume fracrion, []
  real(r8), pointer :: VLPO4(:)    => null()  !PO4 non-band volume fracrion, []
  real(r8), pointer :: VLNOB(:)    => null()  !NO3 band volume fraction, []
  real(r8), pointer :: VLNH4(:)    => null()  !NH4 non-band volume fraction, []
  real(r8), pointer :: VLNHB(:)    => null()  !NH4 band volume fraction, []
  real(r8), pointer :: VOLY(:)     => null()  !total micropore volume [m3 d-2]
  real(r8), pointer :: VOLI(:)     => null()  !soil micropore ice content   [m3 d-2]
  real(r8), pointer :: VOLW(:)     => null()  !soil micropore water content [m3 d-2]
  real(r8), pointer :: VOLA(:)     => null()  !total volume in micropores [m3 d-2]
  real(r8), pointer :: ZNO3S(:)    => null()  !NO3 non-band micropore, [gN d-2]
  real(r8), pointer :: ZNO3B(:)    => null()  !NO3 band micropore, [Ng d-2]
  real(r8), pointer :: ZNSGL(:)    => null()  !aqueous NH3 diffusivity, [m2 h-1]
  real(r8), pointer :: ZVSGL(:)    => null()  !aqueous N2O diffusivity, [m2 h-1]
  real(r8), pointer :: ZNH4S(:)    => null()  !NH4 non-band micropore, [gN d-2]
  real(r8), pointer :: ZNH4B(:)    => null()  !NH4 band micropore, [gN d-2]
  real(r8), pointer :: Z2OS(:)     => null()  !aqueous N2O micropore, [gN d-2]
  real(r8), pointer :: ZNH3S(:)    => null()  !NH3 non-band micropore, [gN d-2]
  real(r8), pointer :: ZNH3B(:)    => null()  !NH3 band micropore, [g d-2]
  real(r8), pointer :: OQC(:,:)    => null()  !dissolved organic C micropore	[gC d-2]
  real(r8), pointer :: OQN(:,:)    => null()  !dissolved organic N micropore	[gN d-2]
  real(r8), pointer :: OQP(:,:)    => null()  !dissolved organic P micropore	[gP d-2]

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
  real(r8), pointer :: FWODB(:)    => null()  !woody C allocation
  real(r8), pointer :: FWODLN(:)   => null()  !leaf N allocation
  real(r8), pointer :: FWODLP(:)   => null()  !leaf P allocation
  real(r8), pointer :: FWODSP(:)   => null()  !P woody fraction in petiole
  real(r8), pointer :: FWODSN(:)   => null()  !N woody fraction in petiole
  real(r8), pointer :: FWODRN(:)   => null()  !N woody fraction in root
  real(r8), pointer :: FWODRP(:)   => null()  !P woody fraction in root
  real(r8), pointer :: FWODR(:)    => null()  !C woody fraction in root
  real(r8), pointer :: FWOOD(:)    => null()  !woody C allocation
  real(r8), pointer :: FVRN(:)     => null()  !allocation parameter
  real(r8), pointer :: FWOODN(:)   => null()  !woody N allocation
  real(r8), pointer :: FWOODP(:)   => null()  !woody P allocation
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
  real(r8) :: WTSTGT                               !total standing dead C, [g d-2]
  real(r8), pointer :: ZEROL(:)       => null()    !threshold zero for leaf calculation
  real(r8), pointer :: ZEROP(:)       => null()    !threshold zero for p calculation
  real(r8), pointer :: CPOOLN(:,:)    => null()    !root  layer nonstructural N, [g d-2]
  real(r8), pointer :: PPOOLN(:,:)    => null()    !nodule layer nonstructural P, [g d-2]
  real(r8), pointer :: WTNDLN(:,:)    => null()    !root layer nodule N, [g d-2]
  real(r8), pointer :: WTNDL(:,:)     => null()    !root layer nodule mass, [g d-2]
  real(r8), pointer :: WTNDLP(:,:)    => null()    !root layer nodule P, [g d-2]
  real(r8), pointer :: ZPOOLN(:,:)    => null()    !root nodule nonstructural N, [g d-2]
  real(r8), pointer :: WGLFV(:,:)     => null()    !canopy layer leaf C, [g d-2]
  real(r8), pointer :: WTRT2(:,:,:,:) => null()    !root layer C secondary axes, [g d-2]
  real(r8), pointer :: WTRT1(:,:,:,:) => null()    !root layer C primary axes, [g d-2]
  real(r8), pointer :: WTRT2N(:,:,:,:)=> null()    !root layer N secondary axes, [g d-2]
  real(r8), pointer :: WTRT1N(:,:,:,:)=> null()    !root layer N primary axes, [g d-2]
  real(r8), pointer :: WTRT2P(:,:,:,:)=> null()    !root layer P secondary axes, [g d-2]
  real(r8), pointer :: WTRT1P(:,:,:,:)=> null()    !root layer P primary axes, [g d-2]
  real(r8), pointer :: RTWT1(:,:,:)   => null()    !root C primary axes, [g d-2]
  real(r8), pointer :: RTWT1N(:,:,:)  => null()    !root N primary axes, [g d-2]
  real(r8), pointer :: RTWT1P(:,:,:)  => null()    !root P primary axes, [g d-2]
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
  real(r8), pointer :: CCPOLR(:,:,:)  => null()    !root  layer nonstructural C concentration, [g g-1]
  real(r8), pointer :: CZPOLR(:,:,:)  => null()    !root layer nonstructural N concentration, [g g-1]
  real(r8), pointer :: CPPOLR(:,:,:)  => null()    !root layer nonstructural P concentration, [g g-1]
  real(r8), pointer :: CEPOLB(:,:,:)    => null()    !branch nonstructural C concentration, [g d-2]
  real(r8), pointer :: WGNODE(:,:,:)  => null()    !internode C, [g d-2]
  real(r8), pointer :: WGNODN(:,:,:)  => null()    !internode N, [g d-2]
  real(r8), pointer :: WGNODP(:,:,:)  => null()    !nodule P, [g d-2]
  real(r8), pointer :: WGLF(:,:,:)    => null()    !leaf C, [g d-2]
  real(r8), pointer :: WGLFN(:,:,:)   => null()    !leaf N, [g d-2]
  real(r8), pointer :: WGLFP(:,:,:)   => null()    !leaf P, [g d-2]
  real(r8), pointer :: WSLF(:,:,:)    => null()    !layer leaf protein C, [g d-2]
  real(r8), pointer :: WGSHE(:,:,:)   => null()    !sheath C , [g d-2]
  real(r8), pointer :: WGSHN(:,:,:)   => null()    !sheath N, [g d-2]
  real(r8), pointer :: WGSHP(:,:,:)   => null()    !sheath P, [g d-2]
  real(r8), pointer :: WSSHE(:,:,:)   => null()    !layer sheath protein C, [g d-2]
  real(r8), pointer :: WGLFL(:,:,:,:) => null()    !layer leaf C, [g d-2]
  real(r8), pointer :: WGLFLN(:,:,:,:)=> null()    !layer leaf N, [g d-2]
  real(r8), pointer :: WGLFLP(:,:,:,:)=> null()    !leaf layer P, [g d-2]
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
  real(r8), pointer :: WTNDB(:,:)     => null()   !branch nodule C, [g d-2]
  real(r8), pointer :: WTSHEBE(:,:,:) => null()   !branch sheath element , [g d-2]
  real(r8), pointer :: WTEARBE(:,:,:) => null()   !branch ear C, [g d-2]
  real(r8), pointer :: WTHSKBE(:,:,:) => null()   !branch husk element, [g d-2]
  real(r8), pointer :: WTNDBN(:,:)    => null()   !branch nodule N, [g d-2]
  real(r8), pointer :: WTNDBP(:,:)    => null()   !branch nodule P, [g d-2]
  real(r8), pointer :: WTGRBE(:,:,:)  => null()   !branch grain element, [g d-2]
  real(r8), pointer :: WTSTKBE(:,:,:) => null()   !branch stalk element, [g d-2]
  real(r8), pointer :: WTSHTBE(:,:,:)    => null()   !branch shoot C, [g d-2]
  real(r8), pointer :: WGLFX(:,:)     => null()   !branch leaf structural C, [g d-2]
  real(r8), pointer :: WGLFNX(:,:)    => null()   !branch leaf structural N, [g d-2]
  real(r8), pointer :: WGLFPX(:,:)    => null()   !branch leaf structural P, [g d-2]
  real(r8), pointer :: WGSHEX(:,:)    => null()   !branch sheath structural C, [g d-2]
  real(r8), pointer :: WGSHNX(:,:)    => null()   !branch sheath structural N, [g d-2]
  real(r8), pointer :: WGSHPX(:,:)    => null()   !branch sheath structural P, [g d-2]
  real(r8), pointer :: WTSTXB(:,:)    => null()   !branch stalk structural C, [g d-2]
  real(r8), pointer :: WTSTXN(:,:)    => null()   !branch stalk structural N, [g d-2]
  real(r8), pointer :: WTSTXP(:,:)    => null()   !branch stalk structural P, [g d-2]
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
  real(r8) :: TKA       !air temperature, [K]
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
  real(r8), pointer :: UPWTR(:,:,:) => null() !root water uptake, [m2 d-2 h-1]
  real(r8), pointer :: TUPHT(:)     => null()   !total root heat uptake, [MJ d-2]
  real(r8), pointer :: TUPWTR(:)    => null()   !total root water uptake, [m3 d-2]
  real(r8), pointer :: VOLWC(:)     => null()  !canopy surface water content, [m3 d-2]
  real(r8), pointer :: VOLWP(:)     => null()  !canopy water content, [m3 d-2]
  real(r8), pointer :: VHCPC(:)     => null()  !canopy heat capacity, [MJ d-2 K-1]
  real(r8), pointer :: PSIST(:)     => null()  !soil micropore total water potential [MPa]
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
  real(r8), pointer :: TDFOMC(:,:)   =>  null()  !total root C exchange, [gC d-2 h-1]
  real(r8), pointer :: TDFOMN(:,:)   =>  null()  !total root N exchange, [gP d-2 h-1]
  real(r8), pointer :: TDFOMP(:,:)   =>  null()  !total root P exchange, [gP d-2 h-1]
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
  real(r8), pointer :: THGFLA(:)     => null()   !total root-atmosphere H2 flux, [g d-2 h-1]
  real(r8), pointer :: TLH2GP(:)     => null()   !total root-soil H2 flux, [g d-2 h-1]
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
  real(r8), pointer :: CSNT(:,:,:)   => null()   !total litterfall C, [gC d-2 h-1]
  real(r8), pointer :: ZSNT(:,:,:)   => null()   !total litterfall N, [gN d-2 h-1]
  real(r8), pointer :: PSNT(:,:,:)   => null()   !total litterfall P, [gP d-2 h-1]
  real(r8), pointer :: XOQCS(:,:)    => null()  !net microbial DOC flux, [gC d-2 h-1]
  real(r8), pointer :: XOQNS(:,:)    => null()  !net microbial DON flux, [gN d-2 h-1]
  real(r8), pointer :: XOQPS(:,:)    => null()  !net microbial DOP flux, [gP d-2 h-1]
  real(r8), pointer :: CNET(:)       => null()  !canopy net CO2 exchange, [gC d-2 h-1]
  real(r8), pointer :: CARBN(:)      => null()  !total gross CO2 fixation, [gC d-2 ]
  real(r8), pointer :: HESNC(:,:)    => null()  !plant element litterfall, [g d-2 h-1]
  real(r8), pointer :: RCO2Z(:)      => null()  !gaseous CO2 flux fron root disturbance, [gC d-2 h-1]
  real(r8), pointer :: ROXYZ(:)      => null()  !gaseous O2 flux fron root disturbance, [g d-2 h-1]
  real(r8), pointer :: RCH4Z(:)      => null()  !gaseous CH4 flux fron root disturbance, [g d-2 h-1]
  real(r8), pointer :: RN2OZ(:)      => null()  !gaseous N2O flux fron root disturbance, [g d-2 h-1]
  real(r8), pointer :: RNH3Z(:)      => null()  !gaseous NH3 flux fron root disturbance non-band, [g d-2 h-1]
  real(r8), pointer :: RH2GZ(:)      => null()  !gaseous H2 flux fron root disturbance, [g d-2 h-1]
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
  real(r8) :: TCO2Z     !total root CO2 content, [gC d-2]
  real(r8) :: TCH4Z     !total root CH4 content, [gC d-2]
  real(r8) :: TN2OZ     !total root N2O content, [g d-2]
  real(r8) :: TOXYZ     !total root O2 content, [g d-2]
  real(r8) :: TNH3Z     !total root NH3 content, [g d-2]
  real(r8), pointer :: UPOMC(:)         => null()  !total root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8), pointer :: UPOMN(:)         => null()  !total root uptake (+ve) - exudation (-ve) of DON, [g d-2 h-1]
  real(r8), pointer :: UPOMP(:)         => null()  !total root uptake (+ve) - exudation (-ve) of DOP, [g d-2 h-1]
  real(r8), pointer :: UPNF(:)          => null()  !total root N2 fixation, [g d-2 h-1]
  real(r8), pointer :: UPNO3(:)         => null()  !total root uptake of NO3, [g d-2 h-1]
  real(r8), pointer :: UPNH4(:)         => null()  !total root uptake of NH4, [g d-2 h-1]
  real(r8), pointer :: UPH1P(:)         => null()  !total root uptake of HPO4, [g d-2 h-1]
  real(r8), pointer :: UPH2P(:)         => null()  !total root uptake of PO4, [g d-2 h-1]
  real(r8), pointer :: RDFOMC(:,:,:,:)  => null()  !root uptake (+ve) - exudation (-ve) of DOC, [gC d-2 h-1]
  real(r8), pointer :: RDFOMN(:,:,:,:)  => null()  !root uptake (+ve) - exudation (-ve) of DON, [gN d-2 h-1]
  real(r8), pointer :: RDFOMP(:,:,:,:)  => null()  !root uptake (+ve) - exudation (-ve) of DOP, [gP d-2 h-1]
  real(r8), pointer :: HCUPTK(:)        => null()  !net root C uptake (+ve) - exudation (-ve), [gC d-2 h-1]
  real(r8), pointer :: HZUPTK(:)        => null()  !net root N uptake (+ve) - exudation (-ve), [gN d-2 h-1]
  real(r8), pointer :: HPUPTK(:)        => null()  !net root P uptake (+ve) - exudation (-ve), [gP d-2 h-1]
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
  real(r8), pointer :: RCOFLA(:,:,:)    => null()  !gaseous CO2 flux through roots, [g d-2 h-1]
  real(r8), pointer :: ROXFLA(:,:,:)    => null()  !gaseous O2 flux through roots, [g d-2 h-1]
  real(r8), pointer :: RCHFLA(:,:,:)    => null()  !gaseous CH4 flux through roots, [g d-2 h-1]
  real(r8), pointer :: RN2FLA(:,:,:)    => null()  !gaseous N2O flux through roots, [g d-2 h-1]
  real(r8), pointer :: RNHFLA(:,:,:)    => null()  !gaseous NH3 flux through roots, [g d-2 h-1]
  real(r8), pointer :: RHGFLA(:,:,:)    => null()  !gaseous H2 flux through roots, [g d-2 h-1]
  real(r8), pointer :: RCODFA(:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) CO2 flux in roots, [g d-2 h-1]
  real(r8), pointer :: ROXDFA(:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) O2 flux in roots, [g d-2 h-1]
  real(r8), pointer :: RCHDFA(:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) CH4 flux in roots, [g d-2 h-1]
  real(r8), pointer :: RN2DFA(:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) N2O flux in roots, [g d-2 h-1]
  real(r8), pointer :: RNHDFA(:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) NH3 flux in roots, [g d-2 h-1]
  real(r8), pointer :: RHGDFA(:,:,:)    => null()  !dissolution (+ve) - volatilization (-ve) H2 flux in roots, [g d-2 h-1]
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
  real(r8), pointer :: CO2P(:,:,:)      => null()  !root aqueous CO2 content, [g d-2 ]
  real(r8), pointer :: CO2A(:,:,:)      => null()  !root gaseous CO2 content, [g d-2 ]
  real(r8), pointer :: CH4P(:,:,:)      => null()  !root aqueous CH4 content, [g d-2 ]
  real(r8), pointer :: CH4A(:,:,:)      => null()  !root gaseous CH4 content, [g d-2 ]
  real(r8), pointer :: H2GP(:,:,:)      => null()  !aqueous H2 content of roots, [g d-2]
  real(r8), pointer :: H2GA(:,:,:)      => null()  !gaseous H2 content of roots, [g d-2]
  real(r8), pointer :: OXYP(:,:,:)      => null()  !root aqueous O2 content, [g d-2 ]
  real(r8), pointer :: OXYA(:,:,:)      => null()  !root gaseous O2 content, [g d-2 ]
  real(r8), pointer :: Z2OP(:,:,:)      => null()  !root aqueous N2O content, [g d-2 ]
  real(r8), pointer :: Z2OA(:,:,:)      => null()  !root gaseous N2O content, [g d-2 ]
  real(r8), pointer :: ZH3P(:,:,:)      => null()  !root aqueous NH3 content, [g d-2 ]
  real(r8), pointer :: ZH3A(:,:,:)      => null()  !root gaseous NH3 content, [g d-2 ]
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
  real(r8), pointer :: TLOXYP(:)        => null()   !total root internal O2 flux, [g d-2 h-1]
  real(r8), pointer :: TLN2OP(:)        => null()   !total root internal N2O flux, [gN d-2 h-1]
  real(r8), pointer :: TCOFLA(:)        => null()   !total internal root CO2 flux , [gC d-2 h-1]
  real(r8), pointer :: TNHFLA(:)        => null()   !total internal root NH3 flux , [gN d-2 h-1]
  real(r8), pointer :: TOXFLA(:)        => null()   !total internal root O2 flux , [g d-2 h-1]
  real(r8), pointer :: TN2FLA(:)        => null()   !total internal root N2O flux , [gN d-2 h-1]
  real(r8), pointer :: TCHFLA(:)        => null()   !total internal root CH4 flux , [gC d-2 h-1]
  real(r8), pointer :: TLCH4P(:)        => null()   !total root internal CH4 flux, [gC d-2 h-1]
  real(r8), pointer :: TLCO2P(:)        => null()   !total root internal CO2 flux, [gC d-2 h-1]
  real(r8), pointer :: TLNH3P(:)        => null()   !total root internal NH3 flux, [gN d-2 h-1]
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


  allocate(this%CO2P(2,JZ1,JP1))
  allocate(this%CO2A(2,JZ1,JP1))
  allocate(this%CH4P(2,JZ1,JP1))
  allocate(this%CH4A(2,JZ1,JP1))
  allocate(this%H2GP(2,JZ1,JP1))
  allocate(this%H2GA(2,JZ1,JP1))
  allocate(this%OXYP(2,JZ1,JP1))
  allocate(this%OXYA(2,JZ1,JP1))
  allocate(this%TUPNF(JZ1))
  allocate(this%ROXSK(60,0:JZ1))
  allocate(this%RDFOMC(2,0:jcplx11,0:JZ1,JP1))
  allocate(this%RDFOMN(2,0:jcplx11,0:JZ1,JP1))
  allocate(this%RDFOMP(2,0:jcplx11,0:JZ1,JP1))
  allocate(this%HCUPTK(JP1))
  allocate(this%HZUPTK(JP1))
  allocate(this%HPUPTK(JP1))
  allocate(this%TEUPTK(npelms,JP1))
  allocate(this%UPOMC(JP1))
  allocate(this%UPNF(JP1))
  allocate(this%UPNO3(JP1))
  allocate(this%UPNH4(JP1))
  allocate(this%UPOMN(JP1))
  allocate(this%UPH1P(JP1))
  allocate(this%UPH2P(JP1))
  allocate(this%UPOMP(JP1))

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
  allocate(this%RNH3B(JC1,JP1))


  allocate(this%TNHFLA(JZ1))
  allocate(this%TCOFLA(JZ1))
  allocate(this%TN2FLA(JZ1))
  allocate(this%TOXFLA(JZ1))
  allocate(this%TCHFLA(JZ1))
  allocate(this%TLCH4P(JZ1))
  allocate(this%TLOXYP(JZ1))
  allocate(this%TLNH3P(JZ1))
  allocate(this%TLCO2P(JZ1))
  allocate(this%TLN2OP(JZ1))


  allocate(this%RN2FLA(2,JZ1,JP1))
  allocate(this%RNHFLA(2,JZ1,JP1))
  allocate(this%RHGFLA(2,JZ1,JP1))
  allocate(this%RCODFA(2,JZ1,JP1))
  allocate(this%ROXDFA(2,JZ1,JP1))
  allocate(this%RCHDFA(2,JZ1,JP1))
  allocate(this%RN2DFA(2,JZ1,JP1))
  allocate(this%RHGDFA(2,JZ1,JP1))
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
  allocate(this%Z2OP(2,JZ1,JP1))
  allocate(this%Z2OA(2,JZ1,JP1))
  allocate(this%ZH3P(2,JZ1,JP1))
  allocate(this%ZH3A(2,JZ1,JP1))
  allocate(this%UPMNPO(2,JP1))
  allocate(this%UPMXPO(2,JP1))
  allocate(this%UPKMPO(2,JP1))
  allocate(this%UPMNZO(2,JP1))
  allocate(this%UPMXZO(2,JP1))
  allocate(this%UPKMZO(2,JP1))
  allocate(this%UPMNZH(2,JP1))
  allocate(this%UPMXZH(2,JP1))
  allocate(this%UPKMZH(2,JP1))
  allocate(this%RNHDFA(2,JZ1,JP1))
  allocate(this%RCO2P(2,JZ1,JP1))
  allocate(this%RUPOXP(2,JZ1,JP1))
  allocate(this%RCO2S(2,JZ1,JP1))
  allocate(this%RUPOXS(2,JZ1,JP1))
  allocate(this%RUPCHS(2,JZ1,JP1))
  allocate(this%RUPN2S(2,JZ1,JP1))
  allocate(this%RUPN3S(2,JZ1,JP1))
  allocate(this%RUPN3B(2,JZ1,JP1))
  allocate(this%RUPHGS(2,JZ1,JP1))
  allocate(this%RCOFLA(2,JZ1,JP1))
  allocate(this%ROXFLA(2,JZ1,JP1))
  allocate(this%RCHFLA(2,JZ1,JP1))
  end subroutine plt_rootbgc_init
!----------------------------------------------------------------------

  subroutine plt_rootbgc_destroy(this)

  implicit none
  class(plant_rootbgc_type) :: this

!  if(allocated(CO2P))deallocate(CO2P)
!  if(allocated(CO2A))deallocate(CO2A)
!  if(allocated(CH4P))deallocate(CH4P)
!  if(allocated(CH4A))deallocate(CH4A)
!  if(allocated(H2GP))deallocate(H2GP)
!  if(allocated(H2GA))deallocate(H2GA)
!  if(allocated(OXYP))deallocate(OXYP)
!  if(allocated(OXYA))deallocate(OXYA)
!  if(allocated(TUPNF))deallocate(TUPNF)
!  if(allocated(ROXSK))deallocate(ROXSK)
!  if(allocated(RDFOMC))deallocate(RDFOMC)
!  if(allocated(RDFOMN))deallocate(RDFOMN)
!  if(allocated(RDFOMP))deallocate(RDFOMP)
!  if(allocated(HCUPTK))deallocate(HCUPTK)
!  if(allocated(HZUPTK))deallocate(HZUPTK)
!  if(allocated(HPUPTK))deallocate(HPUPTK)
!  if(allocated(TEUPTK))deallocate(TEUPTK)
!  if(allocated(UPOMC))deallocate(UPOMC)
!  if(allocated(UPNF))deallocate(UPNF)
!  if(allocated(UPNO3))deallocate(UPNO3)
!  if(allocated(UPNH4))deallocate(UPNH4)
!  if(allocated(UPOMN))deallocate(UPOMN)
!  if(allocated(UPH1P))deallocate(UPH1P)
!  if(allocated(UPH2P))deallocate(UPH2P)
!  if(allocated(UPOMP))deallocate(UPOMP)
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
!  if(allocated(Z2OP))deallocate(Z2OP)
!  if(allocated(Z2OA))deallocate(Z2OA)
!  if(allocated(ZH3P))deallocate(ZH3P)
!  if(allocated(ZH3A))deallocate(ZH3A)
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
  allocate(this%PP(JP1))
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
  allocate(this%THGFLA(JZ1))
  allocate(this%TLH2GP(JZ1))
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

  allocate(this%RN2OZ(JP1))
  allocate(this%RNH3Z(JP1))
  allocate(this%RH2GZ(JP1))
  allocate(this%CSNT(jcplx11,0:1,0:JZ1))
  allocate(this%ZSNT(jcplx11,0:1,0:JZ1))
  allocate(this%PSNT(jcplx11,0:1,0:JZ1))
  allocate(this%CARBN(JP1))
  allocate(this%XOQCS(0:jcplx11,0:JZ1))
  allocate(this%XOQNS(0:jcplx11,0:JZ1))
  allocate(this%XOQPS(0:jcplx11,0:JZ1))
  allocate(this%CNET(JP1))
  allocate(this%RCO2Z(JP1))
  allocate(this%ROXYZ(JP1))
  allocate(this%RCH4Z(JP1))
  allocate(this%TCO2T(JP1))
  allocate(this%TZUPFX(JP1))
  allocate(this%TNH3C(JP1))
  allocate(this%TESN0(npelms,JP1))
  allocate(this%ZESNC(npelms))
  allocate(this%ZNPP(JP1))
  allocate(this%RNH3C(JP1))
  allocate(this%TDFOMC(0:jcplx11,JZ1))
  allocate(this%TDFOMN(0:jcplx11,JZ1))
  allocate(this%TDFOMP(0:jcplx11,JZ1))
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
  allocate(this%ESNC(jsken,npelms,0:1,0:JZ1,JP1))


  end subroutine plt_bgcrate_init
!----------------------------------------------------------------------
  subroutine plt_bgcrate_destroy(this)

  implicit none
  class(plant_bgcrate_type) :: this


!  if(allocated(TCO2S))deallocate(TCO2S)
!  if(allocated(TCO2P))deallocate(TCO2P)
!  if(allocated(THGFLA))deallocate(THGFLA)
!  if(allocated(TLH2GP))deallocate(TLH2GP)
!  if(allocated(TUPNH4))deallocate(TUPNH4)
!  if(allocated(TUPN2S))deallocate(TUPN2S)
!  if(allocated(TUPOXP))deallocate(TUPOXP)
!  if(allocated(TUPOXS))deallocate(TUPOXS)

!  if(allocated(TUPHGS))deallocate(TUPHGS)
!  if(allocated(TUPCHS))deallocate(TUPCHS)
!  if(allocated(TUPN3S))deallocate(TUPN3S)
!  if(allocated(TUPN3B))deallocate(TUPN3B)
!  if(allocated(TUPH1B))deallocate(TUPH1B)
!  if(allocated(TUPNO3))deallocate(TUPNO3)
!  if(allocated(TUPH1P))deallocate(TUPH1P)
!  if(allocated(TUPNOB))deallocate(TUPNOB)
!  if(allocated(TUPNHB))deallocate(TUPNHB)
!  if(allocated(TUPH2P))deallocate(TUPH2P)
!  if(allocated(TUPH2B))deallocate(TUPH2B)
!  if(allocated(RPO4Y))deallocate(RPO4Y)
!  if(allocated(RPOBY))deallocate(RPOBY)
!  if(allocated(RP14Y))deallocate(RP14Y)
!  if(allocated(RP1BY))deallocate(RP1BY)
!  if(allocated(RNO3Y))deallocate(RNO3Y)
!  if(allocated(RNH4Y))deallocate(RNH4Y)
!  if(allocated(RNHBY))deallocate(RNHBY)
!  if(allocated(RN3BY))deallocate(RN3BY)
!  if(allocated(ROXYF))deallocate(ROXYF)
!  if(allocated(RCO2F))deallocate(RCO2F)
!  if(allocated(ROXYL))deallocate(ROXYL)
!  if(allocated(ROXYY))deallocate(ROXYY)
!  if(allocated(RN2OZ))deallocate(RN2OZ)
!  if(allocated(RNH3Z))deallocate(RNH3Z)
!  if(allocated(RH2GZ))deallocate(RH2GZ)
!  if(allocated(CSNT))deallocate(CSNT)
!  if(allocated(ZSNT))deallocate(ZSNT)
!  if(allocated(PSNT))deallocate(PSNT)
!  if(allocated(CARBN))deallocate(CARBN)
!  if(allocated(XOQCS))deallocate(XOQCS)
!  if(allocated(XOQNS))deallocate(XOQNS)
!  if(allocated(XOQPS))deallocate(XOQPS)
!  if(allocated(CNET))deallocate(CNET)
!  if(allocated(RCO2Z))deallocate(RCO2Z)
!  if(allocated(ROXYZ))deallocate(ROXYZ)
!  if(allocated(RCH4Z))deallocate(RCH4Z)
!  if(allocated(TCO2T))deallocate(TCO2T)
!  if(allocated(TZUPFX))deallocate(TZUPFX)
!  if(allocated(TNH3C))deallocate(TNH3C)
!  if(allocated(TESN0))deallocate(TESN0)
!  if(allocated(ZNPP))deallocate(ZNPP)
!  call Destroy(ZESNC)
!  if(allocated(RNH3C))deallocate(RNH3C)
!  if(allocated(HESNC))deallocate(HESNC)
!  if(allocated(TESNC))deallocate(TESNC)
!  if(allocated(ESNC))deallocate(ESNC)
!  if(allocated(TCO2A))deallocate(TCO2A)
!  if(allocated(RUPNF))deallocate(RUPNF)
!  if(allocated(TDFOMP))deallocate(TDFOMP)
!  if(allocated(TDFOMN))deallocate(TDFOMN)
!  if(allocated(TDFOMC))deallocate(TDFOMC)
!  if(allocated(RP1BX))deallocate(RP1BX)
!  if(allocated(RNO3X))deallocate(RNO3X)
!  if(allocated(RP14X))deallocate(RP14X)
!  if(allocated(RNH4X))deallocate(RNH4X)
!  if(allocated(ROXYX))deallocate(ROXYX)
!  if(allocated(RPO4X))deallocate(RPO4X)
!  if(allocated(RNHBX))deallocate(RNHBX)
!  if(allocated(RN3BX))deallocate(RN3BX)
!  if(allocated(RPOBX))deallocate(RPOBX)
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
  allocate(this%PSIST(0:JZ1))
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
  allocate(this%PSIRT(2,JZ1,JP1))
  allocate(this%PSIRO(2,JZ1,JP1))
  allocate(this%PSIRG(2,JZ1,JP1))
  allocate(this%UPWTR(2,JZ1,JP1))
  allocate(this%TUPWTR(0:JZ1))
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
!  if(allocated(UPWTR))deallocate(UPWTR)
!  if(allocated(TUPWTR))deallocate(TUPWTR)
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
  allocate(this%GRWTB(JC1,JP1))
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
  allocate(this%FWOOD(0:1))
  allocate(this%FWODLP(0:1))
  allocate(this%FWODLN(0:1))
  allocate(this%FWODR(0:1))
  allocate(this%FWODSP(0:1))
  allocate(this%FWODB(0:1))
  allocate(this%FWODSN(0:1))
  allocate(this%FWODRN(0:1))
  allocate(this%FWODRP(0:1))
  allocate(this%FWOODN(0:1))
  allocate(this%FWOODP(0:1))

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
!  if(allocated(FWODLP))deallocate(FWODLP)
!  if(allocated(FWODLN))deallocate(FWODLN)
!  if(allocated(FWOOD))deallocate(FWOOD)
!  if(allocated(FWOODP))deallocate(FWOODP)
!  if(allocated(FWOODN))deallocate(FWOODN)
!  if(allocated(FWODR))deallocate(FWODR)
!  if(allocated(FWODRN))deallocate(FWODRN)
!  if(allocated(FWODSN))deallocate(FWODSN)
!  if(allocated(FWODSP))deallocate(FWODSP)
!  if(allocated(FWODB))deallocate(FWODB)
!  if(allocated(FWODRP))deallocate(FWODRP)
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

  allocate(this%WTNDLN(JZ1,JP1))
  allocate(this%WTNDL(JZ1,JP1))
  allocate(this%WTNDLP(JZ1,JP1))
  allocate(this%ZPOOLN(JZ1,JP1))
  allocate(this%WGLFV(JC1,JP1))
  allocate(this%CPOOLN(JZ1,JP1))
  allocate(this%PPOOLN(JZ1,JP1))
  allocate(this%WTSTDE(jsken,npelms,JP1))
  allocate(this%WTRT2(2,JZ1,JC1,JP1))
  allocate(this%WTRT1(2,JZ1,JC1,JP1))
  allocate(this%WTRT2N(2,JZ1,JC1,JP1))
  allocate(this%WTRT1N(2,JZ1,JC1,JP1))
  allocate(this%WTRT2P(2,JZ1,JC1,JP1))
  allocate(this%WTRT1P(2,JZ1,JC1,JP1))
  allocate(this%CEPOLP(npelms,JP1))
  allocate(this%EPOOLP(npelms,JP1))
  allocate(this%EPOLNP(npelms,JP1))
  allocate(this%CCPLNP(JP1))
  allocate(this%CWSRTL(2,JZ1,JP1))
  allocate(this%WSRTL(2,JZ1,JP1))
  allocate(this%WTRTL(2,JZ1,JP1))
  allocate(this%WTRTD(2,JZ1,JP1))
  allocate(this%EPOOLR(npelms,2,JZ1,JP1))
  allocate(this%CCPOLR(2,JZ1,JP1))
  allocate(this%CZPOLR(2,JZ1,JP1))
  allocate(this%CPPOLR(2,JZ1,JP1))
  allocate(this%EPOOL(JC1,npelms,JP1))
  allocate(this%WVSTK(JP1))
  allocate(this%WSLF(0:JNODS1,JC1,JP1))
  allocate(this%WSSHE(0:JNODS1,JC1,JP1))
  allocate(this%WGNODE(0:JNODS1,JC1,JP1))
  allocate(this%WGNODN(0:JNODS1,JC1,JP1))
  allocate(this%WGNODP(0:JNODS1,JC1,JP1))
  allocate(this%WGLF(0:JNODS1,JC1,JP1))
  allocate(this%WGLFN(0:JNODS1,JC1,JP1))
  allocate(this%WGLFP(0:JNODS1,JC1,JP1))
  allocate(this%WGSHE(0:JNODS1,JC1,JP1))
  allocate(this%WGSHN(0:JNODS1,JC1,JP1))
  allocate(this%WGSHP(0:JNODS1,JC1,JP1))
  allocate(this%WGLFL(JC1,0:JNODS1,JC1,JP1))
  allocate(this%WGLFLN(JC1,0:JNODS1,JC1,JP1))
  allocate(this%WGLFLP(JC1,0:JNODS1,JC1,JP1))
  allocate(this%WVSTKB(JC1,JP1))
  allocate(this%EPOLNB(JC1,npelms,JP1))
  allocate(this%CEPOLB(JC1,npelms,JP1))
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
  allocate(this%WGLFX(JC1,JP1))
  allocate(this%WGLFNX(JC1,JP1))
  allocate(this%WGLFPX(JC1,JP1))
  allocate(this%WGSHEX(JC1,JP1))
  allocate(this%WGSHNX(JC1,JP1))
  allocate(this%WGSHPX(JC1,JP1))
  allocate(this%WTSHTE(npelms,JP1))
  allocate(this%WTSHTA(JP1))
  allocate(this%WTLFE(npelms,JP1))
  allocate(this%WTSTDI(JP1))
  allocate(this%WTLSB(JC1,JP1))
  allocate(this%WTRSVBE(JC1,npelms,JP1))
  allocate(this%WTLFBE(JC1,npelms,JP1))
  allocate(this%WTNDB(JC1,JP1))
  allocate(this%WTSHEBE(JC1,npelms,JP1))
  allocate(this%WTEARBE(JC1,npelms,JP1))
  allocate(this%WTHSKBE(JC1,npelms,JP1))
  allocate(this%WTNDBN(JC1,JP1))
  allocate(this%WTNDBP(JC1,JP1))
  allocate(this%WTGRBE(JC1,npelms,JP1))
  allocate(this%WTSTKBE(JC1,npelms,JP1))
  allocate(this%WTSHTBE(JC1,npelms,JP1))
  allocate(this%WTSTXB(JC1,JP1))
  allocate(this%WTSTXN(JC1,JP1))
  allocate(this%WTSTXP(JC1,JP1))
  allocate(this%RTWT1(2,JC1,JP1))
  allocate(this%RTWT1N(2,JC1,JP1))
  allocate(this%RTWT1P(2,JC1,JP1))

  end subroutine plt_biom_init
!----------------------------------------------------------------------

  subroutine plt_biom_destroy(this)

  implicit none
  class(plant_biom_type) :: this

!  if(allocated(ZEROL))deallocate(ZEROL)
!  if(allocated(ZEROP))deallocate(ZEROP)
!  if(allocated(CPOOLN))deallocate(CPOOLN)
!  if(allocated(PPOOLN))deallocate(PPOOLN)
!  if(allocated(WTNDLN))deallocate(WTNDLN)
!  if(allocated(WTNDL))deallocate(WTNDL)
!  if(allocated(WTNDLP))deallocate(WTNDLP)
!  if(allocated(ZPOOLN))deallocate(ZPOOLN)
!  if(allocated(WGLFV))deallocate(WGLFV)
!  if(allocated(RTWT1))deallocate(RTWT1)
!  if(allocated(RTWT1N))deallocate(RTWT1N)
!  if(allocated(RTWT1P))deallocate(RTWT1P)
!  if(allocated(WTRT1))deallocate(WTRT1)
!  if(allocated(WTSTDE))deallocate(WTSTDE)
!  if(allocated(WTRT2))deallocate(WTRT2)
!  if(allocated(WTRT2N))deallocate(WTRT2N)
!  if(allocated(WTRT1N))deallocate(WTRT1N)
!  if(allocated(WTRT2P))deallocate(WTRT2P)
!  if(allocated(WTRT1P))deallocate(WTRT1P)
!  if(allocated(CEPOLP))deallocate(CEPOLP)
!  if(allocated(EPOOLP))deallocate(EPOOLP)
!  if(allocated(EPOLNP))deallocate(EPOLNP)
!  if(allocated(CCPLNP))deallocate(CCPLNP)
!  if(allocated(CWSRTL))deallocate(CWSRTL)
!  if(allocated(WSRTL))deallocate(WSRTL)
!  if(allocated(WTRTL))deallocate(WTRTL)
!  if(allocated(WTRTD))deallocate(WTRTD)
!  if(allocated(EPOOLR))deallocate(EPOOLR)
!  if(allocated(CCPOLR))deallocate(CCPOLR)
!  if(allocated(CZPOLR))deallocate(CZPOLR)
!  if(allocated(CPPOLR))deallocate(CPPOLR)
!  if(allocated(WVSTK))deallocate(WVSTK)
!  if(allocated(WGLF))deallocate(WGLF)
!  if(allocated(WGLFN))deallocate(WGLFN)
!  if(allocated(WGLFP))deallocate(WGLFP)
!  if(allocated(WSLF))deallocate(WSLF)
!  if(allocated(WSSHE))deallocate(WSSHE)
!  if(allocated(WGLFL))deallocate(WGLFL)
!  if(allocated(WGLFLN))deallocate(WGLFLN)
!  if(allocated(WGLFLP))deallocate(WGLFLP)
!  if(allocated(WGNODE))deallocate(WGNODE)
!  if(allocated(WGNODN))deallocate(WGNODN)
!  if(allocated(WGNODP))deallocate(WGNODP)
!  if(allocated(WGSHE))deallocate(WGSHE)
!  if(allocated(WGSHN))deallocate(WGSHN)
!  if(allocated(WGSHP))deallocate(WGSHP)
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
!  if(allocated(WGLFX))deallocate(WGLFX)
!  if(allocated(WGLFNX))deallocate(WGLFNX)
!  if(allocated(WGLFPX))deallocate(WGLFPX)
!  if(allocated(WGSHEX))deallocate(WGSHEX)
!  if(allocated(WGSHNX))deallocate(WGSHNX)
!  if(allocated(WGSHPX))deallocate(WGSHPX)
!  if(allocated(WTLSB))deallocate(WTLSB)
!  if(allocated(WTRSVBE))deallocate(WTRSVBE)
!  if(allocated(WTLFBE))deallocate(WTLFBE)
!  if(allocated(WTNDB))deallocate(WTNDB)
!  if(allocated(WTSHEBE))deallocate(WTSHEBE)
!  if(allocated(WTEARBE))deallocate(WTEARBE)
!  if(allocated(WTHSKBE))deallocate(WTHSKBE)
!  if(allocated(WTNDBN))deallocate(WTNDBN)
!  if(allocated(WTNDBP))deallocate(WTNDBP)
!  if(allocated(WTGRBE))deallocate(WTGRBE)
!  if(allocated(WTSTKBE))deallocate(WTSTKBE)
!  if(allocated(WTSHTBE))deallocate(WTSHTBE)
!  if(allocated(WTSTXB))deallocate(WTSTXB)
!  if(allocated(WTSTXN))deallocate(WTSTXN)
!  if(allocated(WTSTXP))deallocate(WTSTXP)
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

  allocate(this%FOSRH(0:jcplx11,0:JZ1))
  allocate(this%CFOPC(0:Jlitgrp,jsken,JP1))
  allocate(this%CFOPN(0:Jlitgrp,jsken,JP1))
  allocate(this%CFOPP(0:Jlitgrp,jsken,JP1))
  allocate(this%TFND(0:JZ1))
  allocate(this%THETPM(60,0:JZ1))
  allocate(this%DFGS(60,0:JZ1))
  allocate(this%VLPOB(0:JZ1))
  allocate(this%VLNO3(0:JZ1))
  allocate(this%VLPO4(0:JZ1))
  allocate(this%VOLY(0:JZ1))
  allocate(this%VOLI(0:JZ1))
  allocate(this%VOLW(0:JZ1))
  allocate(this%VOLA(0:JZ1))
  allocate(this%VLNOB(0:JZ1))
  allocate(this%VLNH4(0:JZ1))
  allocate(this%VLNHB(0:JZ1))
  allocate(this%ZOSGL(0:JZ1))
  allocate(this%OQC(0:jcplx11,0:JZ1))
  allocate(this%OQN(0:jcplx11,0:JZ1))
  allocate(this%OQP(0:jcplx11,0:JZ1))
  allocate(this%ZNO3S(0:JZ1))
  allocate(this%ZNO3B(0:JZ1))
  allocate(this%ZNSGL(0:JZ1))
  allocate(this%ZNH4S(0:JZ1))
  allocate(this%ZNH4B(0:JZ1))
  allocate(this%Z2OS(0:JZ1))
  allocate(this%ZNH3S(0:JZ1))
  allocate(this%ZNH3B(0:JZ1))
  allocate(this%H1PO4(0:JZ1))
  allocate(this%H2PO4(0:JZ1))
  allocate(this%H1POB(0:JZ1))
  allocate(this%H2POB(0:JZ1))
  allocate(this%H2GS(0:JZ1))
  allocate(this%HLSGL(0:JZ1))
  allocate(this%OXYG(0:JZ1))
  allocate(this%OXYS(0:JZ1))
  allocate(this%OLSGL(0:JZ1))
  allocate(this%POSGL(0:JZ1))
  allocate(this%CCH4G(0:JZ1))
  allocate(this%CZ2OG(0:JZ1))
  allocate(this%CNH3G(0:JZ1))
  allocate(this%CH2GG(0:JZ1))
  allocate(this%CORGC(0:JZ1))
  allocate(this%Z2SGL(JZ1))
  allocate(this%ZHSGL(JZ1))
  allocate(this%CH2P4(0:JZ1))
  allocate(this%CH1P4B(0:JZ1))
  allocate(this%CH2P4B(0:JZ1))
  allocate(this%CNO3S(0:JZ1))
  allocate(this%CH1P4(0:JZ1))
  allocate(this%CNO3B(0:JZ1))
  allocate(this%CNH4S(0:JZ1))
  allocate(this%CNH4B(0:JZ1))
  allocate(this%CO2G(0:JZ1))
  allocate(this%CO2S(0:JZ1))
  allocate(this%CH4S(0:JZ1))
  allocate(this%CCH4S(0:JZ1))
  allocate(this%CZ2OS(0:JZ1))
  allocate(this%CNH3S(0:JZ1))
  allocate(this%CNH3B(0:JZ1))
  allocate(this%CH2GS(0:JZ1))
  allocate(this%CLSGL(0:JZ1))
  allocate(this%CQSGL(0:JZ1))
  allocate(this%VOLX(0:JZ1))
  allocate(this%THETW(0:JZ1))
  allocate(this%THETY(0:JZ1))
  allocate(this%SCO2L(0:JZ1))
  allocate(this%SOXYL(0:JZ1))
  allocate(this%SCH4L(0:JZ1))
  allocate(this%SN2OL(0:JZ1))
  allocate(this%SNH3L(0:JZ1))
  allocate(this%SH2GL(0:JZ1))
  allocate(this%CGSGL(JZ1))
  allocate(this%CHSGL(JZ1))
  allocate(this%HGSGL(JZ1))
  allocate(this%OGSGL(JZ1))
  allocate(this%RSCS(JZ1))
  allocate(this%BKDS(0:JZ1))
  allocate(this%CNDU(JZ1))
  allocate(this%CPO4S(JZ1))
  allocate(this%ZVSGL(0:JZ1))
  end subroutine plt_soilchem_init
!----------------------------------------------------------------------

  subroutine plt_soilchem_destroy(this)
  implicit none
  class(plant_soilchem_type) :: this

!  if(allocated(FOSRH))deallocate(FOSRH)

!  if(allocated(CFOPC))deallocate(CFOPC)
!  if(allocated(CFOPN))deallocate(CFOPN)
!  if(allocated(CFOPP))deallocate(CFOPP)
!  if(allocated(TFND))deallocate(TFND)
!  if(allocated(THETPM))deallocate(THETPM)
!  if(allocated(DFGS))deallocate(DFGS)
!  if(allocated(ZVSGL))deallocate(ZVSGL)
!  if(allocated(SOXYL))deallocate(SOXYL)
!  if(allocated(CPO4S))deallocate(CPO4S)
!  if(allocated(CNDU))deallocate(CNDU)
!  if(allocated(CGSGL))deallocate(CGSGL)
!  if(allocated(CHSGL))deallocate(CHSGL)
!  if(allocated(HGSGL))deallocate(HGSGL)
!  if(allocated(OGSGL))deallocate(OGSGL)
!  if(allocated(RSCS))deallocate(RSCS)
!  if(allocated(SCO2L))deallocate(SCO2L)
!  if(allocated(SCH4L))deallocate(SCH4L)
!  if(allocated(SN2OL))deallocate(SN2OL)
!  if(allocated(SNH3L))deallocate(SNH3L)
!  if(allocated(SH2GL))deallocate(SH2GL)
!  if(allocated(THETW))deallocate(THETW)
!  if(allocated(THETY))deallocate(THETY)
!  if(allocated(VOLX))deallocate(VOLX)
!  if(allocated(CCH4G))deallocate(CCH4G)
!  if(allocated(CZ2OG))deallocate(CZ2OG)
!  if(allocated(CNH3G))deallocate(CNH3G)
!  if(allocated(CH2GG))deallocate(CH2GG)
!  if(allocated(CORGC))deallocate(CORGC)
!  if(allocated(H1PO4))deallocate(H1PO4)
!  if(allocated(H2PO4))deallocate(H2PO4)
!  if(allocated(H1POB))deallocate(H1POB)
!  if(allocated(H2POB))deallocate(H2POB)
!  if(allocated(H2GS))deallocate(H2GS)
!  if(allocated(HLSGL))deallocate(HLSGL)
!  if(allocated(OXYG))deallocate(OXYG)
!  if(allocated(OXYS))deallocate(OXYS)
!  if(allocated(OLSGL))deallocate(OLSGL)
!  if(allocated(POSGL))deallocate(POSGL)
!  if(allocated(VLPOB))deallocate(VLPOB)
!  if(allocated(VLNO3))deallocate(VLNO3)
!  if(allocated(VLPO4))deallocate(VLPO4)
!  if(allocated(VOLY))deallocate(VOLY)
!  if(allocated(VOLI))deallocate(VOLI)
!  if(allocated(VOLW))deallocate(VOLW)
!  if(allocated(VOLA))deallocate(VOLA)
!  if(allocated(VLNOB))deallocate(VLNOB)
!  if(allocated(VLNH4))deallocate(VLNH4)
!  if(allocated(VLNHB))deallocate(VLNHB)
!  if(allocated(ZOSGL))deallocate(ZOSGL)

!  if(allocated(OQC))deallocate(OQC)
!  if(allocated(OQN))deallocate(OQN)
!  if(allocated(OQP))deallocate(OQP)
!  if(allocated(ZNO3S))deallocate(ZNO3S)
!  if(allocated(ZNO3B))deallocate(ZNO3B)
!  if(allocated(ZNSGL))deallocate(ZNSGL)
!  if(allocated(ZNH4S))deallocate(ZNH4S)
!  if(allocated(ZNH4B))deallocate(ZNH4B)
!  if(allocated(Z2OS))deallocate(Z2OS)
!  if(allocated(ZNH3S))deallocate(ZNH3S)
!  if(allocated(ZNH3B))deallocate(ZNH3B)
!  if(allocated(Z2SGL))deallocate(Z2SGL)
!  if(allocated(ZHSGL))deallocate(ZHSGL)
!  if(allocated(BKDS))deallocate(BKDS)
!  if(allocated(CH2P4))deallocate(CH2P4)
!  if(allocated(CH1P4B))deallocate(CH1P4B)
!  if(allocated(CH2P4B))deallocate(CH2P4B)
!  if(allocated(CNO3S))deallocate(CNO3S)
!  if(allocated(CH1P4))deallocate(CH1P4)
!  if(allocated(CNO3B))deallocate(CNO3B)
!  if(allocated(CNH4S))deallocate(CNH4S)
!  if(allocated(CNH4B))deallocate(CNH4B)
!  if(allocated(CO2G))deallocate(CO2G)
!  if(allocated(CO2S))deallocate(CO2S)
!  if(allocated(CH4S))deallocate(CH4S)
!  if(allocated(CCH4S))deallocate(CCH4S)
!  if(allocated(CZ2OS))deallocate(CZ2OS)
!  if(allocated(CNH3S))deallocate(CNH3S)
!  if(allocated(CNH3B))deallocate(CNH3B)
!  if(allocated(CH2GS))deallocate(CH2GS)
!  if(allocated(CLSGL))deallocate(CLSGL)
!  if(allocated(CQSGL))deallocate(CQSGL)

  end subroutine plt_soilchem_destroy
!----------------------------------------------------------------------


  subroutine InitPlantAPIData(JZ,JC,JP,JSA,jcplx1,JLI,JLA,JNODS)

  implicit none
  integer, intent(in) :: JZ,JC,JP,JSA,JCplx1,JLI,JNODS,JLA

  JZ1=JZ
  JC1=JC
  JP1=JP
  JLA1=JLA
  JSA1=JSA
  JLI1   =JLI
  JNODS1 =JNODS
  !the following variable should be consistent with the soil bgc model
  jcplx11=jcplx1
  jsken  =4
  Jlitgrp=5

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
  allocate(this%SURFX(JLI1,JC1,JNODS1,JC1,JP1))
  allocate(this%CPOOL3(JNODS1,JC1,JP1))
  allocate(this%CPOOL4(JNODS1,JC1,JP1))
  allocate(this%CO2B(JNODS1,JC1,JP1))
  allocate(this%COMPL(JNODS1,JC1,JP1))
  allocate(this%CBXN(JNODS1,JC1,JP1))
  allocate(this%CBXN4(JNODS1,JC1,JP1))
  allocate(this%ETGRO(JNODS1,JC1,JP1))
  allocate(this%ETGR4(JNODS1,JC1,JP1))
  allocate(this%FDBK4(JNODS1,JC1,JP1))
  allocate(this%HCOB(JNODS1,JC1,JP1))

  allocate(this%VCGRO(JNODS1,JC1,JP1))
  allocate(this%VGRO(JNODS1,JC1,JP1))
  allocate(this%VCGR4(JNODS1,JC1,JP1))
  allocate(this%VGRO4(JNODS1,JC1,JP1))
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
  allocate(this%FDBK(JC1,JP1))
  allocate(this%FDBKX(JC1,JP1))
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
  allocate(this%OSTR(JP1))
  allocate(this%PB(JP1))
  allocate(this%PR(JP1))
  allocate(this%TFN3(JP1))
  allocate(this%TCZ(JP1))
  allocate(this%TCG(JP1))
  allocate(this%TKG(JP1))
  allocate(this%TCX(JP1))
  allocate(this%WSTR(JP1))
  allocate(this%RCCLX(JC1,JP1))
  allocate(this%RCZLX(JC1,JP1))
  allocate(this%RCPLX(JC1,JP1))
  allocate(this%RCCSX(JC1,JP1))
  allocate(this%RCZSX(JC1,JP1))
  allocate(this%RCPSX(JC1,JP1))

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
  allocate(this%GROUP(JC1,JP1))
  allocate(this%IDTHP(JP1))
  allocate(this%ATRP(JC1,JP1))
  allocate(this%FLGZ(JC1,JP1))
  allocate(this%DGSTGI(JC1,JP1))
  allocate(this%DGSTGF(JC1,JP1))
  allocate(this%GSTGI(JC1,JP1))
  allocate(this%GSTGF(JC1,JP1))
  allocate(this%FLG4(JC1,JP1))
  allocate(this%IWTYP(JP1))
  allocate(this%ISTYP(JP1))
  allocate(this%IBTYP(JP1))
  allocate(this%IDTHR(JP1))
  allocate(this%IDTYP(JP1))
  allocate(this%IPTYP(JP1))
  allocate(this%IFLGI(JP1))
  allocate(this%IGTYP(JP1))
  allocate(this%KVSTG(JC1,JP1))
  allocate(this%KVSTGN(JC1,JP1))
  allocate(this%IDAY(10,JC1,JP1))
  allocate(this%TGSTGI(JC1,JP1))
  allocate(this%TGSTGF(JC1,JP1))
  allocate(this%VSTGX(JC1,JP1))
  allocate(this%XRLA(JP1))
  allocate(this%XRNI(JP1))
  allocate(this%XDL(JP1))
  allocate(this%XPPD(JP1))
  allocate(this%IDTHB(JC1,JP1))
  allocate(this%IFLGP(JC1,JP1))
  allocate(this%IFLGF(JC1,JP1))
  allocate(this%IFLGE(JC1,JP1))
  allocate(this%IFLGA(JC1,JP1))
  allocate(this%IFLGG(JC1,JP1))
  allocate(this%IFLGR(JC1,JP1))
  allocate(this%IFLGQ(JC1,JP1))
  allocate(this%VRNY(JC1,JP1))
  allocate(this%VRNZ(JC1,JP1))
  allocate(this%VRNS(JC1,JP1))
  allocate(this%VRNL(JC1,JP1))
  allocate(this%VRNF(JC1,JP1))
  allocate(this%VRNX(JC1,JP1))


  end subroutine plt_pheno_init
!------------------------------------------------------------------------

  subroutine plt_pheno_destroy(this)
  implicit none
  class(plant_pheno_type) :: this

!  if(allocated(CTC))deallocate(CTC)
!  if(allocated(OFFST))deallocate(OFFST)
!  if(allocated(GROUPI))deallocate(GROUPI)
!  if(allocated(OSTR))deallocate(OSTR)
!  if(allocated(PB))deallocate(PB)
!  if(allocated(PR))deallocate(PR)
!  if(allocated(TFN3))deallocate(TFN3)
!  if(allocated(TCZ))deallocate(TCZ)
!  if(allocated(TCG))deallocate(TCG)
!  if(allocated(TKG))deallocate(TKG)
!  if(allocated(TCX))deallocate(TCX)
!  if(allocated(WSTR))deallocate(WSTR)
!  if(allocated(RCCLX))deallocate(RCCLX)
!  if(allocated(RCZLX))deallocate(RCZLX)
!  if(allocated(RCPLX))deallocate(RCPLX)
!  if(allocated(RCCSX))deallocate(RCCSX)
!  if(allocated(RCZSX))deallocate(RCZSX)
!  if(allocated(RCPSX))deallocate(RCPSX)
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

  allocate(this%RTARP(2,JZ1,JP1))
  allocate(this%RTDNP(2,JZ1,JP1))
  allocate(this%RTVLP(2,JZ1,JP1))
  allocate(this%RTVLW(2,JZ1,JP1))
  allocate(this%RTN1(2,JZ1,JP1))
  allocate(this%RTNL(2,JZ1,JP1))
  allocate(this%GRDM(JP1))
  allocate(this%RTDNT(JZ1))
  allocate(this%RTFQ(JP1))
  allocate(this%CFI(JP1))
  allocate(this%CFX(JP1))
  allocate(this%HTCTL(JP1))
  allocate(this%PORTX(2,JP1))
  allocate(this%WDLF(JP1))
  allocate(this%SDMX(JP1))
  allocate(this%STMX(JP1))
  allocate(this%GRMX(JP1))
  allocate(this%RRAD1X(2,JP1))
  allocate(this%RRAD2X(2,JP1))
  allocate(this%RRADP(2,JP1))
  allocate(this%RRAD1(2,JZ1,JP1))
  allocate(this%RRAD2(2,JZ1,JP1))
  allocate(this%RRAD1M(2,JP1))
  allocate(this%RRAD2M(2,JP1))

  allocate(this%RTDP1(2,JC1,JP1))
  allocate(this%RTLGP(2,JZ1,JP1))
  allocate(this%RTLGA(2,JZ1,JP1))
  allocate(this%RTLG1X(2,JP1))
  allocate(this%RTLG2X(2,JP1))
  allocate(this%RTLG1(2,JZ1,JC1,JP1))
  allocate(this%RTLG2(2,JZ1,JC1,JP1))
  allocate(this%RTN2(2,JZ1,JC1,JP1))
  allocate(this%INTYP(JP1))
  allocate(this%MY(JP1))
  allocate(this%HTSTZ(JP1))
  allocate(this%KLEAF(JC1,JP1))
  allocate(this%VSTG(JC1,JP1))
  allocate(this%NG(JP1))
  allocate(this%ZC(JP1))
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
  allocate(this%PSTG(JC1,JP1))
  allocate(this%PSTGI(JC1,JP1))
  allocate(this%PSTGF(JC1,JP1))
  allocate(this%ANGBR(JP1))
  allocate(this%SURF(JLI1,JC1,JNODS1,JC1,JP1))
  allocate(this%KLEAFX(JC1,JP1))
  allocate(this%NBTB(JC1,JP1))
  allocate(this%GRNXB(JC1,JP1))
  allocate(this%HTSHEX(JC1,JP1))
  allocate(this%ARLFZ(JC1,JP1))
  allocate(this%ARLFB(JC1,JP1))
  allocate(this%ANGSH(JP1))
  allocate(this%CLASS(JLI1,JP1))
  allocate(this%ARSTV(JC1,JP1))
  allocate(this%ARLFV(JC1,JP1))
  allocate(this%ARLF1(0:JNODS1,JC1,JP1))
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
  allocate(this%HTNODX(0:JNODS1,JC1,JP1))
  allocate(this%HTSHE(0:JNODS1,JC1,JP1))
  allocate(this%HTNODE(0:JNODS1,JC1,JP1))
  allocate(this%SURFB(JLI1,JC1,JC1,JP1))
  allocate(this%ARLFL(JC1,0:JNODS1,JC1,JP1))
  allocate(this%ARSTK(JC1,JC1,JP1))
  allocate(this%NI(JP1))
  allocate(this%GRNOB(JC1,JP1))
  allocate(this%CF(JP1))
  allocate(this%DMVL(2,JP1))
  allocate(this%PORT(2,JP1))
  allocate(this%RTAR2X(2,JP1))
  allocate(this%RTAR1X(2,JP1))
  allocate(this%RSRR(2,JP1))
  allocate(this%RSRA(2,JP1))
  allocate(this%IRTYP(JP1))
  end subroutine plt_morph_init

!------------------------------------------------------------------------

  subroutine plt_morph_destroy(this)
  implicit none
  class(plant_morph_type) :: this

!  if(allocated(RTDNP))deallocate(RTDNP)
!  if(allocated(RTVLP))deallocate(RTVLP)
!  if(allocated(RTVLW))deallocate(RTVLW)
!  if(allocated(RTN1))deallocate(RTN1)
!  if(allocated(RTNL))deallocate(RTNL)
!  if(allocated(GRDM))deallocate(GRDM)
!  if(allocated(RTDNT))deallocate(RTDNT)
!  if(allocated(RTFQ))deallocate(RTFQ)
!  if(allocated(CFI)) deallocate(CFI)
!  if(allocated(CFX))deallocate(CFX)
!  if(allocated(HTCTL))deallocate(HTCTL)
!  if(allocated(PORTX))deallocate(PORTX)
!  if(allocated(WDLF))deallocate(WDLF)
!  if(allocated(SDMX))deallocate(SDMX)
!  if(allocated(STMX))deallocate(STMX)
!  if(allocated(IRTYP))deallocate(IRTYP)
!  if(allocated(GRMX))deallocate(GRMX)
!  if(allocated(RSRR))deallocate(RSRR)
!  if(allocated(RSRA))deallocate(RSRA)
!  if(allocated(RRAD1X))deallocate(RRAD1X)
!  if(allocated(RRAD2X))deallocate(RRAD2X)
!  if(allocated(RRADP))deallocate(RRADP)
!  if(allocated(RRAD1))deallocate(RRAD1)
!  if(allocated(RRAD2))deallocate(RRAD2)
!  if(allocated(RRAD1M))deallocate(RRAD1M)
!  if(allocated(RRAD2M))deallocate(RRAD2M)

!  if(allocated(RTLGP))deallocate(RTLGP)
!  if(allocated(RTLGA))deallocate(RTLGA)
!  if(allocated(RTLG2X))deallocate(RTLG2X)
!  if(allocated(RTLG1X))deallocate(RTLG1X)
!  if(allocated(RTDP1))deallocate(RTDP1)
!  if(allocated(RTAR2X))deallocate(RTAR2X)
!  if(allocated(RTAR1X))deallocate(RTAR1X)
!  if(allocated(PORT))deallocate(PORT)
!  if(allocated(DMVL))deallocate(DMVL)
!  if(allocated(RTLG1))deallocate(RTLG1)
!  if(allocated(RTLG2))deallocate(RTLG2)
!  if(allocated(RTN2))deallocate(RTN2)
!  if(allocated(INTYP))deallocate(INTYP)
!  if(allocated(MY)) deallocate(MY)
!  if(allocated(HTSTZ))deallocate(HTSTZ)
!  if(allocated(KLEAF))deallocate(KLEAF)
!  if(allocated(NG)) deallocate(NG)
!  if(allocated(ZC))deallocate(ZC)
!  if(allocated(XTLI))deallocate(XTLI)
!  if(allocated(ZL))deallocate(ZL)
!  if(allocated(ARSTP))deallocate(ARSTP)
!  if(allocated(ARLFP))deallocate(ARLFP)
!  if(allocated(NB1)) deallocate(NB1)
!  if(allocated(NIX))deallocate(NIX)
!  if(allocated(NRT)) deallocate(NRT)
!  if(allocated(NNOD))deallocate(NNOD)
!  if(allocated(NBT)) deallocate(NBT)
!  if(allocated(NBR)) deallocate(NBR)
!  if(allocated(NINR))deallocate(NINR)
!  if(allocated(PSTG))deallocate(PSTG)
!  if(allocated(PSTGI))deallocate(PSTGI)
!  if(allocated(PSTGF))deallocate(PSTGF)
!  if(allocated(ANGBR))deallocate(ANGBR)
!  if(allocated(SURF))deallocate(SURF)
!  if(allocated(KLEAFX))deallocate(KLEAFX)
!  if(allocated(NBTB))deallocate(NBTB)
!  if(allocated(GRNXB))deallocate(GRNXB)
!  if(allocated(HTSHEX))deallocate(HTSHEX)
!  if(allocated(ARLFZ))deallocate(ARLFZ)
!  if(allocated(ARLFB))deallocate(ARLFB)
!  if(allocated(ANGSH))deallocate(ANGSH)
!  if(allocated(CLASS))deallocate(CLASS)
!  if(allocated(ARSTV))deallocate(ARSTV)
!  if(allocated(ARLFV))deallocate(ARLFV)
!  if(allocated(ARLF1))deallocate(ARLF1)
!  if(allocated(GRNO))deallocate(GRNO)
!  if(allocated(SDPTH))deallocate(SDPTH)
!  if(allocated(SDPTHI))deallocate(SDPTHI)
!  if(allocated(SDLG))deallocate(SDLG)
!  if(allocated(SDVL))deallocate(SDVL)
!  if(allocated(SDAR))deallocate(SDAR)
!  if(allocated(ARSTT))deallocate(ARSTT)
!  if(allocated(SSL1))deallocate(SSL1)
!  if(allocated(SNL1))deallocate(SNL1)
!  if(allocated(SLA1))deallocate(SLA1)
!  if(allocated(ARLFT))deallocate(ARLFT)
!  if(allocated(ARLFS))deallocate(ARLFS)
!  if(allocated(HTNODX))deallocate(HTNODX)
!  if(allocated(HTSHE))deallocate(HTSHE)
!  if(allocated(HTNODE))deallocate(HTNODE)
!  if(allocated(SURFB))deallocate(SURFB)
!  if(allocated(ARLFL))deallocate(ARLFL)
!  if(allocated(ARSTK))deallocate(ARSTK)
!  if(allocated(NI)) deallocate(NI)
!  if(allocated(GRNOB))deallocate(GRNOB)
!  if(allocated(CF))  deallocate(CF)
!  if(allocated(VSTG))deallocate(VSTG)
  end subroutine plt_morph_destroy
end module PlantAPIData
