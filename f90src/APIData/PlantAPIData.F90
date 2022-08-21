module PlantAPIData
  use data_kind_mod, only : r8 => SHR_KIND_R8
implicit none
  save
  public

  real(r8) :: ALATs1      !latitude	[degrees]
  real(r8) :: ATCAs1      !mean annual air temperature, [oC]
  real(r8) :: ALTs1       !altitude of grid cell, [m]

! grid configuration
  integer  :: JP1         !number of plants
  integer  :: JSA1        !number of sectors for the sky azimuth  [0,2*pi]
  integer  :: jcplx11     !number of organo-microbial complexes
  integer  :: JLA1        !number of sectors for the leaf azimuth, [0,pi]
  integer  :: JC1         !number of canopy layers
  integer  :: JZ1         !number of soil layers
  integer  :: JLI1        !number of sectors for the leaf zenith [0,pi/2]
  integer  :: JNODS1      !number of canopy nodes
!begin_data

  real(r8) :: VOLWSs1     !water volume in snowpack, [m3 d-2]
  real(r8) :: CCO2EIs1    !initial atmospheric CO2 concentration, [g m-3]
  real(r8) :: CO2EIs1     !initial atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: COXYEs1     !atmospheric O2 concentration, [g m-3]
  real(r8) :: RABs1       !isothermal boundary layer resistance, [h m-1]
  real(r8) :: VOLSSs1     !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8) :: TSHCs1      !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: CNETXs1     !total net canopy CO2 exchange, [g d-2 h-1]
  real(r8) :: ZNOONs1     !time of solar noon, [h]
  real(r8) :: CO2Es1      !atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: CCO2Es1     !atmospheric CO2 concentration, [g m-3]
  real(r8) :: CCH4Es1     !atmospheric CH4 concentration, [g m-3]
  real(r8) :: CZ2OEs1     !atmospheric N2O concentration, [g m-3]
  real(r8) :: CH2GEs1     !atmospheric H2 concentration, [g m-3]
  real(r8) :: CNH3Es1     !atmospheric NH3 concentration, [g m-3]
  real(r8) :: DYLXs1      !daylength of previous day, [h]

  real(r8) :: DYLNs1      !daylength, [h]
  real(r8) :: DPTHSs1     !snowpack depth, [m]
  real(r8) :: DYLMs1      !maximum daylength, [h]
  real(r8) :: ZSs1        !initial soil surface roughness height, [m]
  real(r8) :: OXYEs1      !atmospheric O2 concentration, [umol mol-1]
  real(r8) :: PPTs1       !total plant population, [d-2]
  real(r8) :: RIBs1       !Richardson number for calculating boundary layer resistance, [-]
  real(r8) :: RECOs1      !ecosystem respiration, [g d-2 h-1]
  real(r8) :: POROS1s1    !top layer soil porosity
  real(r8) :: UAs1        !wind speed, [m h-1]
  real(r8) :: TGHs1       !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8) :: TBALCs1     !total plant C balance	gC d-2
  real(r8) :: TBALNs1     !total plant N balance	gN d-2
  real(r8) :: TBALPs1     !total plant P balance	gP d-2
  real(r8) :: Z0s1        !wind speed measurement height, [m]
  real(r8) :: VHCPWXs1    !snowpack heat capacity from previous time step, [MJ d-2 K-1]
  real(r8) :: VHCPW1s1    !snowpack heat capacity, [MJ m-3 K-1]
  real(r8) :: TCCANs1     !total net CO2 fixation, [gC d-2]
  real(r8) :: TCH4Zs1     !total root CH4 content, [gC d-2]
  real(r8) :: TENGYCs1    !total canopy heat content, [MJ  d-2]
  real(r8) :: THRMCs1     !total canopy LW emission, [MJ d-2 h-1]
  real(r8) :: TEVAPPs1    !total canopy evaporation + transpiration, [m3 d-2]
  real(r8) :: TCO2Zs1     !total root CO2 content, [gC d-2]
  real(r8) :: TEVAPCs1    !total canopy evaporation, [m3 d-2]
  real(r8) :: TH2GZs1     !total root H2 flux, [g d-2]
  real(r8) :: THFLXCs1    !total canopy heat flux, [MJ  d-2]

  real(r8) :: TN2OZs1     !total root N2O content, [g d-2]
  real(r8) :: TLEs1       !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8) :: TOXYZs1     !total root O2 content, [g d-2]
  real(r8) :: TNH3Zs1     !total root NH3 content, [g d-2]
  real(r8) :: TVOLWPs1    !total canopy water content, [m3 d-2]
  real(r8) :: ZZSNCs1     !total litterfall N, [g d-2 h-1]
  real(r8) :: ZPSNCs1     !total litterfall P, [g d-2 h-1]
  real(r8) :: ZCSNCs1     !total litterfall C, [g d-2 h-1]
  real(r8) :: VOLISs1     !ice volume in snowpack, [m3 d-2]
  real(r8) :: TKWs1       !snow temperature, [K]
  real(r8) :: TVOLWCs1    !canopy surface water content, [m3 d-2]
  real(r8) :: TSHs1       !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8) :: TNBPs1      !total NBP, [g d-2]
  real(r8) :: TGPPs1      !ecosystem GPP, [g d-2 h-1]
  real(r8) :: TKAs1       !air temperature, [K]
  real(r8) :: TLECs1      !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: TRAUs1      !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8) :: UVOLOs1     !total subsurface water flux, [m3 d-2]
  real(r8) :: VPAs1       !vapor concentration, [m3 m-3]
  real(r8) :: XHVSTCs1    !ecosystem harvest C, [gC d-2]
  real(r8) :: XHVSTNs1    !ecosystem harvest N, [gN d-2]
  real(r8) :: XHVSTPs1    !ecosystem harvest P, [gP d-2]
  real(r8) :: XCORPs1     !factor for surface litter incorporation and soil mixing
  real(r8) :: ZEROS2s1    !threshold zero
  real(r8) :: ZEROSs1     !threshold zero
  real(r8) :: ZEROs1      !threshold zero
  real(r8) :: ZERO2s1     !threshold zero
  real(r8) :: ZRs1        !canopy surface roughness height, [m]
  real(r8) :: ZDs1        !zero plane displacement height, [m]
  integer :: IYTYPs1      !fertilizer release type from fertilizer input file
  integer :: IETYPs1      !Koppen climate zone
  integer :: IFLGTs1      !number of active PFT
  integer :: NLs1         !lowest soil layer number
  integer :: NP0s1        !intitial number of plant species
  integer :: NJs1         !maximum root layer number
  integer :: NPs1         !number of plant species
  integer :: NUs1         !soil surface layer number
  integer :: LYRCs1       !number of days in current year
  integer :: IYRCs1       !current year
  real(r8) :: DCORPs1     !soil mixing fraction with tillage, [-]
  real(r8) :: VOLWOUs1    !total subsurface water flux	m3 d-2
  integer :: ITILLs1      !soil disturbance type, [-]
  character(len=16), allocatable :: DATAPs1(:) !parameter file name
  CHARACTER(len=16), allocatable :: DATAs1(:)  !pft file
  real(r8), allocatable :: AREA3s1(:)    !soil cross section area (vertical plan defined by its normal direction)

  integer,  allocatable :: IYR0s1(:)     !year of planting
  integer,  allocatable :: IDAYXs1(:)    !!alternate day of planting
  integer,  allocatable :: IYRXs1(:)     !alternate year of planting
  integer,  allocatable :: IYRYs1(:)     !alternate year of harvest
  integer,  allocatable :: INTYPs1(:)    !N2 fixation type
  integer,  allocatable :: IDAY0s1(:)    !day of planting
  integer,  allocatable :: IDAYHs1(:)    !day of harvest
  integer,  allocatable :: IYRHs1(:)     !year of harvest

  integer,  allocatable :: IHVSTs1(:)    !type of harvest
  integer,  allocatable :: JHVSTs1(:)    !flag for stand replacing disturbance
  integer,  allocatable :: ICTYPs1(:)    !plant photosynthetic type (C3 or C4)
  integer,  allocatable :: IDATAs1(:)    !time keeper
  real(r8), allocatable :: FERTs1(:)     !fertilizer application, [g m-2]
  real(r8), allocatable :: TDFOMCs1(:,:) !total root C exchange, [gC d-2 h-1]
  real(r8), allocatable :: TDFOMNs1(:,:) !total root N exchange, [gP d-2 h-1]
  real(r8), allocatable :: TDFOMPs1(:,:) !total root P exchange, [gP d-2 h-1]
  real(r8), allocatable :: RUPNFs1(:,:)  !root N2 fixation, [gN d-2 h-1]

  real(r8), allocatable :: TCCs1(:)      !canopy temperature, [oC]
  real(r8), allocatable :: DTKCs1(:)     !change in canopy temperature, [K]
  real(r8), allocatable :: ENGYXs1(:)    !canopy heat storage from previous time step, [MJ d-2]
  real(r8), allocatable :: TKCs1(:)      !canopy temperature, [K]
  real(r8), allocatable :: PSILOs1(:)    !canopy osmotic water potential, [Mpa]
  real(r8), allocatable :: OSMOs1(:)     !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8), allocatable :: TCO2As1(:)    !total autotrophic respiration, [gC d-2 ]
  real(r8), allocatable :: CTRANs1(:)    !total transpiration, [m H2O d-2]
  real(r8), allocatable :: ZTYPs1(:)     !plant thermal adaptation zone, [-]
  real(r8), allocatable :: ZTYPIs1(:)    !initial plant thermal adaptation zone, [-]
  real(r8), allocatable :: HTCs1(:)      !temperature above which seed set is adversely affected, [oC]
  real(r8), allocatable :: SSTXs1(:)     !sensitivity to HTC (seeds oC-1 above HTC)

  real(r8), allocatable :: CDPTHZs1(:)   !depth to bottom of soil layer from  surface of grid cell [m]
  real(r8), allocatable :: PPIs1(:)      !initial plant population, [m-2]
  real(r8), allocatable :: PPZs1(:)      !plant population at seeding, [m-2]
  real(r8), allocatable :: PPXs1(:)      !plant population, [m-2]
  real(r8), allocatable :: PPs1(:)       !plant population, [d-2]
  real(r8), allocatable :: DPTHZs1(:)    !depth to middle of soil layer from  surface of grid cell [m]
  real(r8), allocatable :: CFIs1(:)      !initial clumping factor for self-shading in canopy layer, [-]
  real(r8), allocatable :: FMPRs1(:)     !micropore fraction

  real(r8), allocatable :: ROXYXs1(:)    !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8), allocatable :: RNHBXs1(:)    !total root + microbial NH4 uptake band, [gN d-2 h-1]
  real(r8), allocatable :: RP14Xs1(:)    !HPO4 demand in non-band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), allocatable :: RP1BXs1(:)    !HPO4 demand in band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), allocatable :: RNO3Xs1(:)    !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), allocatable :: RNH4Xs1(:)    !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), allocatable :: RN3BXs1(:)    !total root + microbial NO3 uptake band, [gN d-2 h-1]
  real(r8), allocatable :: RPO4Xs1(:)    !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), allocatable :: RPOBXs1(:)    !total root + microbial PO4 uptake band, [gP d-2 h-1]
  real(r8), allocatable :: RPO4Ys1(:)    !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  real(r8), allocatable :: RPOBYs1(:)    !total root + microbial PO4 uptake band, [gP d-2 h-1]
  real(r8), allocatable :: RP14Ys1(:)    !HPO4 demand in non-band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), allocatable :: RP1BYs1(:)    !HPO4 demand in band by all microbial,root,myco populations, [gP d-2 h-1]
  real(r8), allocatable :: RNO3Ys1(:)    !total root + microbial NO3 uptake non-band, [gN d-2 h-1]
  real(r8), allocatable :: RNH4Ys1(:)    !total root + microbial NH4 uptake non-band, [gN d-2 h-1]
  real(r8), allocatable :: RNHBYs1(:)    !total root + microbial NH4 uptake band, [gN d-2 h-1]
  real(r8), allocatable :: RN3BYs1(:)    !total root + microbial NO3 uptake band, [gN d-2 h-1]
  real(r8), allocatable :: ROXYFs1(:)    !net gaseous O2 flux, [g d-2 h-1]
  real(r8), allocatable :: RCO2Fs1(:)    !net gaseous CO2 flux, [g d-2 h-1]
  real(r8), allocatable :: ROXYLs1(:)    !net aqueous O2 flux, [g d-2 h-1]
  real(r8), allocatable :: ROXYYs1(:)    !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8), allocatable :: PSILZs1(:)    !minimum daily canopy water potential, [MPa]
  real(r8), allocatable :: RSMNs1(:)     !canopy minimum stomatal resistance, [s m-1]
  real(r8), allocatable :: DLYR3s1(:)    !vertical thickness of soil layer [m]
  real(r8), allocatable :: SFLXCs1(:)    !canopy sensible heat flux, [MJ d-2 h-1]
  real(r8), allocatable :: RCs1(:)       !canopy stomatal resistance, [h m-1]
  real(r8), allocatable :: TKCZs1(:)     !canopy temperature, [K]
  real(r8), allocatable :: RAs1(:)       !canopy boundary layer resistance, [h m-1]
  real(r8), allocatable :: SO2s1(:)      !leaf O2 solubility, [uM /umol mol-1]

  real(r8), allocatable :: TKSs1(:)      !mean annual soil temperature, [K]
  real(r8), allocatable :: TFNDs1(:)     !temperature effect on diffusivity

  real(r8), allocatable :: TUPNO3s1(:)   !total root-soil NO3 flux non-band, [gN d-2 h-1]
  real(r8), allocatable :: TUPH2Bs1(:)   !total root-soil PO4 flux band, [gP d-2 h-1]
  real(r8), allocatable :: TUPH1Bs1(:)   !soil-root exch of HPO4 in band [gP d-2 h-1]
  real(r8), allocatable :: TUPN2Ss1(:)   !total root-soil N2O flux, [gN d-2 h-1]
  real(r8), allocatable :: TCO2Ss1(:)    !total root-soil CO2 flux, [gC d-2 h-1]
  real(r8), allocatable :: TCO2Ps1(:)    !total root CO2 flux, [gC d-2 h-1]
  real(r8), allocatable :: TUPOXPs1(:)   !total root internal O2 flux, [g d-2 h-1]
  real(r8), allocatable :: THGFLAs1(:)   !total root-atmosphere H2 flux, [g d-2 h-1]
  real(r8), allocatable :: TLH2GPs1(:)   !total root-soil H2 flux, [g d-2 h-1]
  real(r8), allocatable :: TCHFLAs1(:)   !total internal root CH4 flux , [gC d-2 h-1]
  real(r8), allocatable :: TLCH4Ps1(:)   !total root internal CH4 flux, [gC d-2 h-1]
  real(r8), allocatable :: TLCO2Ps1(:)   !total root internal CO2 flux, [gC d-2 h-1]
  real(r8), allocatable :: TLNH3Ps1(:)   !total root internal NH3 flux, [gN d-2 h-1]
  real(r8), allocatable :: TUPHTs1(:)    !total root heat uptake, [MJ d-2]
  real(r8), allocatable :: TUPWTRs1(:)   !total root water uptake, [m3 d-2]
  real(r8), allocatable :: RTDNTs1(:)    !total root length density, [m m-3]
  real(r8), allocatable :: TUPNFs1(:)    !total root N2 fixation, [g d-2 h-1]
  real(r8), allocatable :: TLOXYPs1(:)   !total root internal O2 flux, [g d-2 h-1]
  real(r8), allocatable :: TLN2OPs1(:)   !total root internal N2O flux, [gN d-2 h-1]
  real(r8), allocatable :: TCOFLAs1(:)   !total internal root CO2 flux , [gC d-2 h-1]
  real(r8), allocatable :: TNHFLAs1(:)   !total internal root NH3 flux , [gN d-2 h-1]
  real(r8), allocatable :: TOXFLAs1(:)   !total internal root O2 flux , [g d-2 h-1]
  real(r8), allocatable :: TN2FLAs1(:)   !total internal root N2O flux , [gN d-2 h-1]
  real(r8), allocatable :: TUPOXSs1(:)   !total root-soil O2 flux, [g d-2 h-1]
  real(r8), allocatable :: TUPHGSs1(:)   !total root-soil H2 flux, [g d-2 h-1]
  real(r8), allocatable :: TUPCHSs1(:)   !total root-soil CH4 flux, [gC d-2 h-1]
  real(r8), allocatable :: TUPN3Bs1(:)   !total root-soil NH3 flux band, [gN d-2 h-1]
  real(r8), allocatable :: TUPN3Ss1(:)   !total root-soil NH3 flux non-band, [gN d-2 h-1]
  real(r8), allocatable :: TUPH2Ps1(:)   !total root-soil PO4 flux non-band, [gP d-2 h-1]
  real(r8), allocatable :: TUPH1Ps1(:)   !soil-root exch of HPO4 in non-band, [gP d-2 h-1]
  real(r8), allocatable :: TUPNH4s1(:)   !total root-soil NH4 flux non-band, [gN d-2 h-1]
  real(r8), allocatable :: TUPNHBs1(:)   !total root-soil NH4 flux band, [gN d-2 h-1]
  real(r8), allocatable :: TUPNOBs1(:)   !total root-soil NO3 flux band, [gN d-2 h-1]
  real(r8), allocatable :: CSNTs1(:,:,:)    !total litterfall C, [gC d-2 h-1]
  real(r8), allocatable :: ZSNTs1(:,:,:)    !total litterfall N, [gN d-2 h-1]
  real(r8), allocatable :: PSNTs1(:,:,:)    !total litterfall P, [gP d-2 h-1]
  real(r8), allocatable :: RDFOMCs1(:,:,:,:)  !root uptake (+ve) - exudation (-ve) of DOC, [gC d-2 h-1]
  real(r8), allocatable :: RDFOMNs1(:,:,:,:)  !root uptake (+ve) - exudation (-ve) of DON, [gN d-2 h-1]
  real(r8), allocatable :: RDFOMPs1(:,:,:,:)  !root uptake (+ve) - exudation (-ve) of DOP, [gP d-2 h-1]
  real(r8), allocatable :: FOSRHs1(:,:)   !fraction of total organic C in complex, [-]

  real(r8), allocatable :: XOQCSs1(:,:)   !net microbial DOC flux, [gC d-2 h-1]
  real(r8), allocatable :: XOQNSs1(:,:)   !net microbial DON flux, [gN d-2 h-1]
  real(r8), allocatable :: XOQPSs1(:,:)   !net microbial DOP flux, [gP d-2 h-1]
  real(r8), allocatable :: FLWCs1(:)      !water flux into canopy, [m3 d-2 h-1]
  real(r8), allocatable :: EHVSTs1(:,:,:) !harvest efficiency, [-]
  real(r8), allocatable :: HVSTs1(:)      !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  real(r8), allocatable :: THINs1(:)      !thinning of plant population, [-]
  real(r8), allocatable :: BALCs1(:)      !plant C balance, [gC d-2]
  real(r8), allocatable :: BALNs1(:)      !plant N balance, [gN d-2]
  real(r8), allocatable :: BALPs1(:)      !plant P balance, [gP d-2]

  real(r8), allocatable :: CTCs1(:)       !temperature below which seed set is adversely affected, [oC]
  real(r8), allocatable :: CNETs1(:)      !canopy net CO2 exchange, [gC d-2 h-1]
  real(r8), allocatable :: CARBNs1(:)     !total gross CO2 fixation, [gC d-2 ]
  real(r8), allocatable :: RSMXs1(:)      !maximum stomatal resistance to vapor, [s m-1]

  real(r8), allocatable :: CFXs1(:)       !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8), allocatable :: CWSRTs1(:)     !fraction of remobilizable nonstructural biomass in root, [-]
  real(r8), allocatable :: DCO2s1(:)      !gaesous CO2 concentration difference across stomates, [umol m-3]
  real(r8), allocatable :: FMOLs1(:)      !total gas concentration, [mol m-3]
  real(r8), allocatable :: FNODs1(:)      !parameter for allocation of growth to nodes, [-]
  real(r8), allocatable :: GROUPIs1(:)    !acclimated plant maturity group, [-]
  real(r8), allocatable :: HTCTLs1(:)     !cotyledon height, [m]
  real(r8), allocatable :: HCSNCs1(:)     !plant C litterfall, [gC d-2 h-1]
  real(r8), allocatable :: HZSNCs1(:)     !plant N litterfall, [gN d-2 h-1]
  real(r8), allocatable :: HPSNCs1(:)     !plant P litterfall, [gP d-2 h-1]
  real(r8), allocatable :: HVSTCs1(:)     !plant C harvest, [gC d-2 ]
  real(r8), allocatable :: HVSTNs1(:)     !plant N harvest, [gN d-2 ]
  real(r8), allocatable :: HVSTPs1(:)     !plant P harvest, [gP d-2 ]
  real(r8), allocatable :: HCUPTKs1(:)    !net root C uptake (+ve) - exudation (-ve), [gC d-2 h-1]
  real(r8), allocatable :: HZUPTKs1(:)    !net root N uptake (+ve) - exudation (-ve), [gN d-2 h-1]
  real(r8), allocatable :: HPUPTKs1(:)    !net root P uptake (+ve) - exudation (-ve), [gP d-2 h-1]

  real(r8), allocatable :: OFFSTs1(:)     !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
  real(r8), allocatable :: OSTRs1(:)      !plant O2 stress indicator, []
  real(r8), allocatable :: PSILTs1(:)     !canopy total water potential , [Mpa]
  real(r8), allocatable :: PBs1(:)        !branch nonstructural C content required for new branch, [gC gC-1]
  real(r8), allocatable :: PRs1(:)        !threshold root nonstructural C content for initiating new root axis, [gC gC-1]
  real(r8), allocatable :: PSILGs1(:)     !canopy turgor water potential, [MPa]
  real(r8), allocatable :: RCMXs1(:)      !maximum stomatal resistance to CO2, [s h-1]
  real(r8), allocatable :: RCO2Zs1(:)     !gaseous CO2 flux fron root disturbance, [gC d-2 h-1]
  real(r8), allocatable :: ROXYZs1(:)     !gaseous O2 flux fron root disturbance, [g d-2 h-1]
  real(r8), allocatable :: RCH4Zs1(:)     !gaseous CH4 flux fron root disturbance, [g d-2 h-1]
  real(r8), allocatable :: RN2OZs1(:)     !gaseous N2O flux fron root disturbance, [g d-2 h-1]
  real(r8), allocatable :: RNH3Zs1(:)     !gaseous NH3 flux fron root disturbance non-band, [g d-2 h-1]
  real(r8), allocatable :: RH2GZs1(:)     !gaseous H2 flux fron root disturbance, [g d-2 h-1]
  real(r8), allocatable :: RTFQs1(:)      !root brancing frequency, [m-1]
  real(r8), allocatable :: RAZs1(:)       !canopy roughness height, [m]

  real(r8), allocatable :: EVAPCs1(:)     !canopy evaporation, [m2 d-2 h-1]
  real(r8), allocatable :: HFLXCs1(:)     !canopy storage heat flux, [MJ d-2 h-1]
  real(r8), allocatable :: EFLXCs1(:)     !canopy latent heat flux, [MJ d-2 h-1]
  real(r8), allocatable :: TCZs1(:)       !threshold temperature for spring leafout/dehardening, [oC]
  real(r8), allocatable :: TCGs1(:)       !canopy growth temperature, [oC]
  real(r8), allocatable :: TCXs1(:)       !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8), allocatable :: TKGs1(:)       !canopy growth temperature, [K]
  real(r8), allocatable :: TFN3s1(:)      !canopy temperature growth function, [-]
  real(r8), allocatable :: TCSN0s1(:)     !total surface litterfall C, [g d-2]
  real(r8), allocatable :: TZSN0s1(:)     !total surface litterfall N, [g d-2]
  real(r8), allocatable :: TPSN0s1(:)     !total surface litterfall P, [g d-2]
  real(r8), allocatable :: TCSNCs1(:)     !total plant C litterfall , [g d-2 ]
  real(r8), allocatable :: TZSNCs1(:)     !total plant N litterfall , [g d-2 ]
  real(r8), allocatable :: TPSNCs1(:)     !total plant P litterfall , [g d-2 ]
  real(r8), allocatable :: TCO2Ts1(:)     !total plant respiration, [gC d-2 ]
  real(r8), allocatable :: TCUPTKs1(:)    !total net root C uptake (+ve) - exudation (-ve), [gC d-2 ]
  real(r8), allocatable :: THVSTCs1(:)    !total plant C harvest, [gC d-2 ]
  real(r8), allocatable :: TNH3Cs1(:)     !total canopy NH3 flux, [gN d-2 ]
  real(r8), allocatable :: TZUPTKs1(:)    !total net root N uptake (+ve) - exudation (-ve), [gN d-2 ]
  real(r8), allocatable :: THVSTNs1(:)    !total plant N harvest, [g d-2 ]
  real(r8), allocatable :: TZUPFXs1(:)    !total plant N2 fixation, [g d-2 ]
  real(r8), allocatable :: TPUPTKs1(:)    !total net root P uptake (+ve) - exudation (-ve), [g d-2 ]
  real(r8), allocatable :: THVSTPs1(:)    !total plant P harvest, [g d-2 ]
  real(r8), allocatable :: UPOMCs1(:)     !total root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8), allocatable :: UPOMNs1(:)     !total root uptake (+ve) - exudation (-ve) of DON, [g d-2 h-1]
  real(r8), allocatable :: UPOMPs1(:)     !total root uptake (+ve) - exudation (-ve) of DOP, [g d-2 h-1]
  real(r8), allocatable :: UPNFs1(:)      !total root N2 fixation, [g d-2 h-1]
  real(r8), allocatable :: UPNO3s1(:)     !total root uptake of NO3, [g d-2 h-1]
  real(r8), allocatable :: UPNH4s1(:)     !total root uptake of NH4, [g d-2 h-1]
  real(r8), allocatable :: UPH1Ps1(:)     !total root uptake of HPO4, [g d-2 h-1]
  real(r8), allocatable :: UPH2Ps1(:)     !total root uptake of PO4, [g d-2 h-1]
  real(r8), allocatable :: VOLWCs1(:)     !canopy surface water content, [m3 d-2]
  real(r8), allocatable :: VOLWPs1(:)     !canopy water content, [m3 d-2]
  real(r8), allocatable :: VHCPCs1(:)     !canopy heat capacity, [MJ d-2 K-1]
  real(r8), allocatable :: VCH4Fs1(:)     !plant CH4 emission from fire, [g d-2 ]
  real(r8), allocatable :: VCO2Fs1(:)     !plant CO2 emission from fire, [g d-2 ]
  real(r8), allocatable :: VN2OFs1(:)     !plant N2O emission from fire, [g d-2 ]
  real(r8), allocatable :: VNH3Fs1(:)     !plant NH3 emission from fire, [g d-2 ]
  real(r8), allocatable :: VPO4Fs1(:)     !plant PO4 emission from fire, [g d-2 ]
  real(r8), allocatable :: WSTRs1(:)      !canopy plant water stress indicator, number of hours PSILT < PSILY, []


  real(r8), allocatable :: ZEROLs1(:)     !threshold zero for leaf calculation
  real(r8), allocatable :: ZEROPs1(:)     !threshold zero for p calculation
  real(r8), allocatable :: ZEROQs1(:)     !threshold zero for uptake calculation
  real(r8), allocatable :: ZNPPs1(:)      !total net primary productivity, [gC d-2]

  real(r8), allocatable :: GRWTBs1(:,:)   !maximum grain C during grain fill, [g d-2]

  integer,  allocatable :: IDAYYs1(:)     !alternate day of harvest, [-]

  real(r8), allocatable :: RCCLXs1(:,:)   !C translocated from leaf during senescence, [g d-2 h-1]
  real(r8), allocatable :: RCZLXs1(:,:)   !N translocated from leaf during senescence, [g d-2 h-1]
  real(r8), allocatable :: RCPLXs1(:,:)   !P translocated from leaf during senescence, [g d-2 h-1]
  real(r8), allocatable :: RCCSXs1(:,:)   !C translocated from sheath during senescence, [g d-2 h-1]
  real(r8), allocatable :: RCZSXs1(:,:)   !N translocated from sheath during senescence, [g d-2 h-1]
  real(r8), allocatable :: RCPSXs1(:,:)   !P translocated from sheath during senescence, [g d-2 h-1]
  real(r8), allocatable :: RNH3Bs1(:,:)   !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]

  real(r8), allocatable :: SURFXs1(:,:,:,:,:)  !leaf irradiated surface area, [m2 d-2]
  real(r8), allocatable :: CPOOL3s1(:,:,:)    !minimum sink strength for nonstructural C transfer, [g d-2]
  real(r8), allocatable :: CPOOL4s1(:,:,:)    !leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
  real(r8), allocatable :: CO2Bs1(:,:,:)      !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), allocatable :: COMPLs1(:,:,:)     !CO2 compensation point, [uM]
  real(r8), allocatable :: CBXNs1(:,:,:)      !carboxylation efficiency, [umol umol-1]
  real(r8), allocatable :: CBXN4s1(:,:,:)     !C4 carboxylation efficiency, [umol umol-1]
  real(r8), allocatable :: ETGROs1(:,:,:)     !maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), allocatable :: ETGR4s1(:,:,:)     !maximum  light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), allocatable :: FDBK4s1(:,:,:)     !down-regulation of C4 photosynthesis, [-]
  real(r8), allocatable :: HCOBs1(:,:,:)      !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), allocatable :: VCGROs1(:,:,:)     !maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), allocatable :: VGROs1(:,:,:)      !carboxylation rate, [umol m-2 s-1]
  real(r8), allocatable :: VCGR4s1(:,:,:)     !maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  real(r8), allocatable :: VGRO4s1(:,:,:)     !C4 carboxylation rate, [umol m-2 s-1]

  real(r8), allocatable :: CO2Ps1(:,:,:)      !root aqueous CO2 content, [g d-2 ]
  real(r8), allocatable :: CO2As1(:,:,:)      !root gaseous CO2 content, [g d-2 ]
  real(r8), allocatable :: CH4Ps1(:,:,:)      !root aqueous CH4 content, [g d-2 ]
  real(r8), allocatable :: CH4As1(:,:,:)      !root gaseous CH4 content, [g d-2 ]
  real(r8), allocatable :: H2GPs1(:,:,:)      !aqueous H2 content of roots, [g d-2]
  real(r8), allocatable :: H2GAs1(:,:,:)      !gaseous H2 content of roots, [g d-2]
  real(r8), allocatable :: OXYPs1(:,:,:)      !root aqueous O2 content, [g d-2 ]
  real(r8), allocatable :: OXYAs1(:,:,:)      !root gaseous O2 content, [g d-2 ]
  real(r8), allocatable :: PSIRTs1(:,:,:)     !root total water potential , [Mpa]
  real(r8), allocatable :: PSIROs1(:,:,:)     !root osmotic water potential , [Mpa]
  real(r8), allocatable :: PSIRGs1(:,:,:)     !root turgor water potential , [Mpa]
  real(r8), allocatable :: RUPNHBs1(:,:,:)    !root uptake of NH4 band, [g d-2 h-1]
  real(r8), allocatable :: RUPNH4s1(:,:,:)    !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8), allocatable :: RUPH2Ps1(:,:,:)    !root uptake of PO4 non-band, [g d-2 h-1]
  real(r8), allocatable :: RUPNOBs1(:,:,:)    !root uptake of NO3 band, [g d-2 h-1]
  real(r8), allocatable :: RUPNO3s1(:,:,:)    !root uptake of NO3 non-band, [g d-2 h-1]
  real(r8), allocatable :: RUPH1Bs1(:,:,:)    !root uptake of HPO4 band
  real(r8), allocatable :: RUPH1Ps1(:,:,:)    !root uptake HPO4 non-band
  real(r8), allocatable :: RUPH2Bs1(:,:,:)    !root uptake of PO4 band, [g d-2 h-1]
  real(r8), allocatable :: RUONHBs1(:,:,:)    !root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RUONH4s1(:,:,:)    !root uptake of NH4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RUOH2Ps1(:,:,:)    !root uptake of PO4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RUONOBs1(:,:,:)    !root uptake of NO3 band unconstrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RUONO3s1(:,:,:)    !root uptake of NO3 non-band unconstrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RUOH1Bs1(:,:,:)    !root HPO4 uptake in band unlimited by O2
  real(r8), allocatable :: RUOH1Ps1(:,:,:)    !root HPO4 uptake in non-band unlimited by O2
  real(r8), allocatable :: RUOH2Bs1(:,:,:)    !root uptake of PO4 band unconstrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RUCNHBs1(:,:,:)    !root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), allocatable :: RUCNH4s1(:,:,:)    !root uptake of NH4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), allocatable :: RUCH2Ps1(:,:,:)    !root uptake of PO4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), allocatable :: RUCNOBs1(:,:,:)    !root uptake of NO3 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), allocatable :: RUCNO3s1(:,:,:)    !root uptake of NO3 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), allocatable :: RUCH1Bs1(:,:,:)    !root HPO4 uptake in band unlimited by nonstructural C
  real(r8), allocatable :: RUCH1Ps1(:,:,:)    !root HPO4 uptake in non-band unlimited by nonstructural C
  real(r8), allocatable :: RUCH2Bs1(:,:,:)    !root uptake of PO4 band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), allocatable :: RCO2Ms1(:,:,:)     !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RCO2Ns1(:,:,:)     !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8), allocatable :: RCO2As1(:,:,:)     !root respiration constrained by O2, [g d-2 h-1]
  real(r8), allocatable :: RTN1s1(:,:,:)      !root layer number primary axes, [d-2]
  real(r8), allocatable :: RTNLs1(:,:,:)      !root layer number axes, [d-2]
  real(r8), allocatable :: RTLGPs1(:,:,:)     !root layer length per plant, [m p-1]
  real(r8), allocatable :: RTDNPs1(:,:,:)     !root layer length density, [m m-3]
  real(r8), allocatable :: RTVLPs1(:,:,:)     !root layer volume air, [m2 d-2]
  real(r8), allocatable :: RTVLWs1(:,:,:)     !root layer volume water, [m2 d-2]
  real(r8), allocatable :: RRAD1s1(:,:,:)     !root layer diameter primary axes, [m ]
  real(r8), allocatable :: RRAD2s1(:,:,:)     !root layer diameter secondary axes, [m ]
  real(r8), allocatable :: RTARPs1(:,:,:)     !root layer area per plant, [m p-1]
  real(r8), allocatable :: RTLGAs1(:,:,:)     !root layer average length, [m]
  real(r8), allocatable :: RCO2Ps1(:,:,:)     !aqueous CO2 flux from roots to root water , [g d-2 h-1]
  real(r8), allocatable :: RUPOXPs1(:,:,:)    !aqueous O2 flux from roots to root water , [g d-2 h-1]
  real(r8), allocatable :: RCO2Ss1(:,:,:)     !aqueous CO2 flux from roots to soil water, [g d-2 h-1]
  real(r8), allocatable :: RUPOXSs1(:,:,:)    !aqueous O2 flux from roots to soil water, [g d-2 h-1]
  real(r8), allocatable :: RUPCHSs1(:,:,:)    !aqueous CH4 flux from roots to soil water, [g d-2 h-1]
  real(r8), allocatable :: RUPN2Ss1(:,:,:)    !aqueous N2O flux from roots to soil water, [g d-2 h-1]
  real(r8), allocatable :: RUPN3Ss1(:,:,:)    !aqueous NH3 flux from roots to soil water non-band, [g d-2 h-1]
  real(r8), allocatable :: RUPN3Bs1(:,:,:)    !aqueous NH3 flux from roots to soil water band, [g d-2 h-1]
  real(r8), allocatable :: RUPHGSs1(:,:,:)    !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8), allocatable :: RCOFLAs1(:,:,:)    !gaseous CO2 flux through roots, [g d-2 h-1]
  real(r8), allocatable :: ROXFLAs1(:,:,:)    !gaseous O2 flux through roots, [g d-2 h-1]
  real(r8), allocatable :: RCHFLAs1(:,:,:)    !gaseous CH4 flux through roots, [g d-2 h-1]
  real(r8), allocatable :: RN2FLAs1(:,:,:)    !gaseous N2O flux through roots, [g d-2 h-1]
  real(r8), allocatable :: RNHFLAs1(:,:,:)    !gaseous NH3 flux through roots, [g d-2 h-1]
  real(r8), allocatable :: RHGFLAs1(:,:,:)    !gaseous H2 flux through roots, [g d-2 h-1]
  real(r8), allocatable :: RCODFAs1(:,:,:)    !dissolution (+ve) - volatilization (-ve) CO2 flux in roots, [g d-2 h-1]
  real(r8), allocatable :: ROXDFAs1(:,:,:)    !dissolution (+ve) - volatilization (-ve) O2 flux in roots, [g d-2 h-1]
  real(r8), allocatable :: RCHDFAs1(:,:,:)    !dissolution (+ve) - volatilization (-ve) CH4 flux in roots, [g d-2 h-1]
  real(r8), allocatable :: RN2DFAs1(:,:,:)    !dissolution (+ve) - volatilization (-ve) N2O flux in roots, [g d-2 h-1]
  real(r8), allocatable :: RNHDFAs1(:,:,:)    !dissolution (+ve) - volatilization (-ve) NH3 flux in roots, [g d-2 h-1]
  real(r8), allocatable :: RHGDFAs1(:,:,:)    !dissolution (+ve) - volatilization (-ve) H2 flux in roots, [g d-2 h-1]
  real(r8), allocatable :: ROXYPs1(:,:,:)     !root  O2 demand from respiration, [g d-2 h-1]
  real(r8), allocatable :: RUNNHPs1(:,:,:)    !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8), allocatable :: RUNNBPs1(:,:,:)    !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8), allocatable :: RUNNOPs1(:,:,:)    !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8), allocatable :: RUNNXPs1(:,:,:)    !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8), allocatable :: RUPP2Ps1(:,:,:)    !root uptake of H2PO4 non-band
  real(r8), allocatable :: RUPP2Bs1(:,:,:)    !root uptake of H2PO4 band
  real(r8), allocatable :: RUPP1Ps1(:,:,:)    !HPO4 demand in non-band by each root population
  real(r8), allocatable :: RUPP1Bs1(:,:,:)    !HPO4 demand in band by each root population
  real(r8), allocatable :: UPWTRs1(:,:,:)     !root water uptake, [m2 d-2 h-1]
  real(r8), allocatable :: WFRs1(:,:,:)       !O2 constraint to root respiration, []

  real(r8), allocatable :: Z2OPs1(:,:,:)      !root aqueous N2O content, [g d-2 ]
  real(r8), allocatable :: Z2OAs1(:,:,:)      !root gaseous N2O content, [g d-2 ]
  real(r8), allocatable :: ZH3Ps1(:,:,:)      !root aqueous NH3 content, [g d-2 ]
  real(r8), allocatable :: ZH3As1(:,:,:)      !root gaseous NH3 content, [g d-2 ]
  real(r8), allocatable :: CPOOLNs1(:,:)      !root  layer nonstructural N, [g d-2]
  real(r8), allocatable :: PPOOLNs1(:,:)      !nodule layer nonstructural P, [g d-2]
  real(r8), allocatable :: TFN4s1(:,:)        !root layer temperature growth functiom, [-]
  real(r8), allocatable :: WTNDLNs1(:,:)      !root layer nodule N, [g d-2]
  real(r8), allocatable :: WTNDLs1(:,:)       !root layer nodule mass, [g d-2]
  real(r8), allocatable :: WTNDLPs1(:,:)      !root layer nodule P, [g d-2]
  real(r8), allocatable :: ZPOOLNs1(:,:)      !root nodule nonstructural N, [g d-2]

  real(r8), allocatable :: WGLFVs1(:,:)       !canopy layer leaf C, [g d-2]
  real(r8), allocatable :: DMVLs1(:,:)        !root volume:mass ratio, [m3 g-1]
  real(r8), allocatable :: PORTs1(:,:)        !root porosity, [m3 m-3]
  real(r8), allocatable :: PORTXs1(:,:)       !power function of root porosity used to calculate root gaseous diffusivity, [-]
  real(r8), allocatable :: RTAR2Xs1(:,:)      !root  cross-sectional area  secondary axes, [m2]
  real(r8), allocatable :: RTAR1Xs1(:,:)      !root cross-sectional area primary axes, [m2]
  real(r8), allocatable :: RRAD1Xs1(:,:)      !root diameter primary axes, [m]
  real(r8), allocatable :: RRAD2Xs1(:,:)      !root diameter secondary axes, [m]
  real(r8), allocatable :: GRDMs1(:)          !grain size at seeding, [g]
  real(r8), allocatable :: VOXYFs1(:)         !plant O2 uptake from fire, [g d-2 ]

  real(r8), allocatable :: PTSHTs1(:)         !shoot-root rate constant for nonstructural C exchange, [h-1]
  real(r8), allocatable :: EPs1(:)            !canopy transpiration, [m2 d-2 h-1]
  real(r8), allocatable :: WDLFs1(:)          !leaf length:width ratio, [-]
  real(r8), allocatable :: GRMXs1(:)          !maximum grain size   , [g]
  real(r8), allocatable :: GFILLs1(:)         !maximum rate of fill per grain, [g h-1]
  integer , allocatable :: IRTYPs1(:)         !grain type (below or above-ground)
  real(r8), allocatable :: STMXs1(:)          !maximum grain node number per branch, [-]
  real(r8), allocatable :: SDMXs1(:)          !maximum grain number per node , [-]

  real(r8), allocatable :: RNH3Cs1(:)         !canopy NH3 flux, [g d-2 h-1]


  real(r8), allocatable :: RTLG1Xs1(:,:)      !specific root length primary axes, [m g-1]
  real(r8), allocatable :: RRADPs1(:,:)       !root internal radius, [m]
  real(r8), allocatable :: RTLG2Xs1(:,:)      !specific root length secondary axes, [m g-1]
  real(r8), allocatable :: RSRRs1(:,:)        !root radial resistivity, [MPa h m-2]
  real(r8), allocatable :: RSRAs1(:,:)        !root axial resistivity, [MPa h m-4]
  real(r8), allocatable :: RRAD1Ms1(:,:)      !maximum radius of primary roots, [m]
  real(r8), allocatable :: RRAD2Ms1(:,:)      !maximum radius of secondary roots, [m]
  real(r8), allocatable :: UPMNPOs1(:,:)      !minimum PO4 concentration for root NH4 uptake, [g m-3]
  real(r8), allocatable :: UPMXPOs1(:,:)      !maximum root PO4 uptake rate, [g m-2 h-1]
  real(r8), allocatable :: UPKMPOs1(:,:)      !Km for root PO4 uptake, [g m-3]
  real(r8), allocatable :: UPMNZOs1(:,:)      !minimum NO3 concentration for root NH4 uptake, [g m-3]
  real(r8), allocatable :: UPMXZOs1(:,:)      !maximum root NO3 uptake rate, [g m-2 h-1]
  real(r8), allocatable :: UPKMZOs1(:,:)      !Km for root NO3 uptake, [g m-3]
  real(r8), allocatable :: UPMNZHs1(:,:)      !minimum NH4 concentration for root NH4 uptake, [g m-3]
  real(r8), allocatable :: UPMXZHs1(:,:)      !maximum root NH4 uptake rate, [g m-2 h-1]
  real(r8), allocatable :: UPKMZHs1(:,:)      !Km for root NH4 uptake, [g m-3]
  real(r8), allocatable :: RTLG1s1(:,:,:,:)   !root layer length primary axes, [m d-2]
  real(r8), allocatable :: RTLG2s1(:,:,:,:)   !root layer length secondary axes, [m d-2]
  real(r8), allocatable :: RTN2s1(:,:,:,:)    !root layer number secondary axes, [d-2]
  real(r8), allocatable :: WTRT2s1(:,:,:,:)   !root layer C secondary axes, [g d-2]
  real(r8), allocatable :: WTRT1s1(:,:,:,:)   !root layer C primary axes, [g d-2]
  real(r8), allocatable :: WTRT2Ns1(:,:,:,:)  !root layer N secondary axes, [g d-2]
  real(r8), allocatable :: WTRT1Ns1(:,:,:,:)  !root layer N primary axes, [g d-2]
  real(r8), allocatable :: WTRT2Ps1(:,:,:,:)  !root layer P secondary axes, [g d-2]
  real(r8), allocatable :: WTRT1Ps1(:,:,:,:)  !root layer P primary axes, [g d-2]
  real(r8), allocatable :: RTDP1s1(:,:,:)     !root layer depth, [m]
  real(r8), allocatable :: RTWT1s1(:,:,:)     !root C primary axes, [g d-2]
  real(r8), allocatable :: RTWT1Ns1(:,:,:)    !root N primary axes, [g d-2]
  real(r8), allocatable :: RTWT1Ps1(:,:,:)    !root P primary axes, [g d-2]
  real(r8), allocatable :: WTSTDGs1(:,:)      !standing dead C fraction, [g d-2]
  real(r8), allocatable :: WTSTDNs1(:,:)      !standing dead N fraction, [g d-2]
  real(r8), allocatable :: WTSTDPs1(:,:)      !standing dead P fraction, [g d-2]
  real(r8), allocatable :: CSNCs1(:,:,:,:)    !plant litterfall C, [g d-2 h-1]
  real(r8), allocatable :: PSNCs1(:,:,:,:)    !litterfall P flux, [g d-2 h-1]
  real(r8), allocatable :: ZSNCs1(:,:,:,:)    !total litterfall N, [g d-2 h-1]
  real(r8), allocatable :: CFOPCs1(:,:,:)     !litter kinetic fraction, [-]
  real(r8), allocatable :: CFOPNs1(:,:,:)     !litterfall kinetic N fraction, [-]
  real(r8), allocatable :: CFOPPs1(:,:,:)     !litter P kinetic fraction, [-]
  real(r8), allocatable :: PARs1(:,:,:,:)     !direct incoming PAR, [umol m-2 s-1]
  real(r8), allocatable :: PARDIFs1(:,:,:,:)  !diffuse incoming PAR, [umol m-2 s-1]
  real(r8), allocatable :: VOLWMs1(:,:)       !soil micropore water content, [m3 d-2]
  real(r8), allocatable :: VOLPMs1(:,:)       !soil air content, [m3 d-2]
  real(r8), allocatable :: TORTs1(:,:)        !soil tortuosity, []
  real(r8), allocatable :: FILMs1(:,:)        !soil water film thickness , [m]
  real(r8), allocatable :: ROXSKs1(:,:)       !total O2 sink, [g d-2 t-1]

  type, public :: plant_photosyns_type
  real(r8), pointer :: ETMXs1(:)   => null()  !cholorophyll activity , [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: CHLs1(:)    => null()  !leaf C3 chlorophyll content, [gC gC-1]
  real(r8), pointer :: PEPCs1(:)   => null()  !leaf PEP carboxylase content, [gC gC-1]
  real(r8), pointer :: CHL4s1(:)   => null()  !leaf C4 chlorophyll content, [gC gC-1]
  real(r8), pointer :: RUBPs1(:)   => null()  !leaf rubisco content, [gC gC-1]
  real(r8), pointer :: VCMX4s1(:)  => null()  !PEP carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: VOMXs1(:)   => null()  !rubisco oxygenase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: VCMXs1(:)   => null()  !rubisco carboxylase activity, [umol g-1 h-1 at 25 oC]
  real(r8), pointer :: XKCO2s1(:)  => null()  !Km for rubisco carboxylase activity, [uM]
  real(r8), pointer :: XKO2s1(:)   => null()  !Km for rubisco oxygenase activity, [uM]
  real(r8), pointer :: FDBKs1(:,:) => null()   !branch down-regulation of CO2 fixation, [-]
  real(r8), pointer :: FDBKXs1(:,:)=> null()   !down-regulation of C4 photosynthesis, [-]
  real(r8), pointer :: CO2Ls1(:)   => null()   !leaf aqueous CO2 concentration, [uM]
  real(r8), pointer :: O2Is1(:)    => null()   !leaf gaseous O2 concentration, [umol m-3]
  real(r8), pointer :: CO2Is1(:)   => null()   !leaf gaseous CO2 concentration, [umol m-3]
  real(r8), pointer :: XKCO2Os1(:) => null()   !leaf aqueous CO2 Km ambient O2, [uM]
  real(r8), pointer :: XKCO2Ls1(:) => null()   !leaf aqueous CO2 Km no O2, [uM]
  real(r8), pointer :: RCSs1(:)    => null()   !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8), pointer :: FCO2s1(:)   => null()   !Ci:Ca ratio, [-]
  real(r8), pointer :: RSMHs1(:)   => null()   !maximum stomatal resistance to vapor, [s h-1]
  real(r8), pointer :: CHILLs1(:)  => null()   !chilling effect on CO2 fixation, [-]
  real(r8), pointer :: SCO2s1(:)   => null()   !leaf CO2 solubility, [uM /umol mol-1]
  real(r8), pointer :: CO2Qs1(:)   => null()   !canopy gaesous CO2 concentration , [umol mol-1]
  real(r8), pointer :: O2Ls1(:)    => null()   !leaf aqueous O2 concentration, [uM]
  real(r8), pointer :: XKCO24s1(:) => null()   !Km for PEP carboxylase activity, [uM]

  contains
    procedure, public :: Init    =>  plt_photo_init
    procedure, public :: Destroy => plt_photo_destroy
  end type plant_photosyns_type

  type, public ::  plant_radiation_type
  real(r8) :: TYSINs1     !sine of sky angles
  real(r8) :: ALBSs1      !snowpack albedo
  real(r8) :: ALBXs1      !Surface albedo
  real(r8) :: RAP0s1      !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8) :: RADYs1      !diffuse shortwave radiation, [W m-2]
  real(r8) :: RADSs1      !direct shortwave radiation, [W m-2]
  real(r8) :: RAPYs1      !diffuse PAR, [umol m-2 s-1]
  real(r8) :: RAPSs1      !direct PAR, [umol m-2 s-1]
  real(r8) :: TRNs1       !ecosystem net radiation, [MJ d-2 h-1]
  real(r8) :: RAD0s1      !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8) :: FRADGs1     !fraction of radiation intercepted by ground surface, [-]
  real(r8) :: RADGs1      !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8) :: GSINs1      !sine of slope, [-]
  real(r8) :: GAZIs1      !azimuth of slope, [-]
  real(r8) :: GCOSs1      !cosine of slope, [-]
  real(r8) :: THRMGXs1    !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8) :: THSs1       !sky longwave radiation , [MJ d-2 h-1]
  real(r8) :: SSINNs1     !sine of solar angle next hour, [-]
  real(r8) :: SSINs1      !sine of solar angle, [-]
  real(r8), pointer :: ALBRs1(:)     => null() !canopy shortwave albedo , [-]
  real(r8), pointer :: ALBPs1(:)     => null() !canopy PAR albedo , [-]
  real(r8), pointer :: TAUSs1(:)     => null() !fraction of radiation intercepted by canopy layer, [-]
  real(r8), pointer :: TAU0s1(:)     => null() !fraction of radiation transmitted by canopy layer, [-]
  real(r8), pointer :: THRM1s1(:)    => null() !canopy longwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: RADCs1(:)     => null() !canopy absorbed shortwave radiation , [MJ d-2 h-1]
  real(r8), pointer :: OMEGXs1(:,:,:)=> null() !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation
  real(r8), pointer :: OMEGAGs1(:)   => null() !sine of solar beam on leaf surface, [-]
  real(r8), pointer :: OMEGAs1(:,:,:)=> null() !sine of indirect sky radiation on leaf surface
  real(r8), pointer :: ZSINs1(:)     => null() !sine of leaf angle
  real(r8), pointer :: ZCOSs1(:)     => null() !cosine of leaf angle
  integer,  pointer :: IALBYs1(:,:,:)=> null() !flag for calculating backscattering of radiation in canopy
  real(r8), pointer :: RAD1s1(:)     => null() !canopy net radiation , [MJ d-2 h-1]
  real(r8), pointer :: ABSRs1(:)     => null() !canopy shortwave absorptivity , [-]
  real(r8), pointer :: ABSPs1(:)     => null() !canopy PAR absorptivity
  real(r8), pointer :: TAURs1(:)     => null() !canopy shortwave transmissivity , [-]
  real(r8), pointer :: TAUPs1(:)     => null() !canopy PAR transmissivity , [-]
  real(r8), pointer :: RADPs1(:)     => null() !canopy absorbed PAR , [umol m-2 s-1]
  real(r8), pointer :: FRADPs1(:)    => null() !fraction of incoming PAR absorbed by canopy, [-]
  contains
    procedure, public :: Init    => plt_rad_init
    procedure, public :: Destroy => plt_rad_destroy
  end type plant_radiation_type

  type, public :: plant_morph_type
  real(r8) :: ARLSSs1                     !stalk area of combined, each PFT canopy
  real(r8) :: ARLFCs1                     !total canopy leaf area, [m2 d-2]
  real(r8) :: ARSTCs1                     !total canopy stem area, [m2 d-2]
  real(r8) :: ZTs1                        !canopy height , [m]
  integer,  pointer :: MYs1(:)           => null() !mycorrhizal type (no or yes)
  real(r8), pointer :: HTSTZs1(:)        => null() !canopy height, [m]
  real(r8), pointer :: ZLs1(:)           => null() !canopy layer height , [m]
  real(r8), pointer :: ARSTPs1(:)        => null() !plant stem area, [m2 d-2]
  real(r8), pointer :: ARLFPs1(:)        => null() !plant leaf area, [m2 d-2]
  integer,  pointer :: NB1s1(:)          => null() !number of main branch
  integer,  pointer :: NIs1(:)           => null() !maximum soil layer number for all root axes
  integer,  pointer :: NIXs1(:)          => null() !maximum soil layer number for all root axes, [-]
  integer,  pointer :: NRTs1(:)          => null() !root primary axis number
  integer,  pointer :: NNODs1(:)         => null() !number of concurrently growing nodes
  integer,  pointer :: NBTs1(:)          => null() !branch number
  integer,  pointer :: NBRs1(:)          => null() !branch number
  integer,  pointer :: NINRs1(:,:)       => null() !maximum soil layer number for root axes, [-]
  real(r8), pointer :: PSTGs1(:,:)       => null() !shoot node number, [-]
  real(r8), pointer :: PSTGIs1(:,:)      => null() !shoot node number at floral initiation, [-]
  real(r8), pointer :: PSTGFs1(:,:)      => null() !shoot node number at anthesis, [-]
  real(r8), pointer :: ANGBRs1(:)        => null() !branching angle, [degree from horizontal]
  real(r8), pointer :: SURFs1(:,:,:,:,:) => null() !leaf surface area, [m2 d-2]
  integer,  pointer :: KLEAFXs1(:,:)     => null() !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer,  pointer :: NBTBs1(:,:)       => null() !branch number, [-]
  real(r8), pointer :: GRNXBs1(:,:)      => null() !branch potential grain number, [d-2]
  real(r8), pointer :: GRNOBs1(:,:)      => null() !branch grain number, [d-2]
  real(r8), pointer :: HTSHEXs1(:,:)     => null() !branch height, [m]
  real(r8), pointer :: ARLFZs1(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: ARLFBs1(:,:)      => null() !branch leaf area, [m2 d-2]
  real(r8), pointer :: ANGSHs1(:)        => null() !sheath angle, [degree from horizontal]
  real(r8), pointer :: CLASSs1(:,:)      => null() !fractionction of leaves in different angle classes, [-]
  real(r8), pointer :: ARSTVs1(:,:)      => null() !canopy layer stem area, [m2 d-2]
  real(r8), pointer :: ARLFVs1(:,:)      => null() !canopy layer leaf area, [m2 d-2]
  real(r8), pointer :: ARLF1s1(:,:,:)    => null() !leaf area, [m2 d-2]
  real(r8), pointer :: GRNOs1(:)         => null() !canopy grain number, [d-2]
  real(r8), pointer :: SDPTHs1(:)        => null() !seeding depth, [m]
  real(r8), pointer :: SDPTHIs1(:)       => null() !planting depth, [m]
  real(r8), pointer :: SDLGs1(:)         => null() !seed length, [m]
  real(r8), pointer :: SDVLs1(:)         => null() !seed volume, [m3 ]
  real(r8), pointer :: SDARs1(:)         => null() !seed surface area, [m2]
  real(r8), pointer :: ARSTTs1(:)        => null() !total stem area, [m2 d-2]
  real(r8), pointer :: SSL1s1(:)         => null() !petiole length:mass during growth, [m gC-1]
  real(r8), pointer :: SNL1s1(:)         => null() !internode length:mass during growth, [m gC-1]
  real(r8), pointer :: SLA1s1(:)         => null() !leaf area:mass during growth, [m2 gC-1]
  real(r8), pointer :: ARLFTs1(:)        => null() !total leaf area, [m2 d-2]
  real(r8), pointer :: ARLFSs1(:)        => null() !plant leaf area, [m2 d-2]
  real(r8), pointer :: HTNODXs1(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: HTSHEs1(:,:,:)    => null() !sheath height, [m]
  real(r8), pointer :: HTNODEs1(:,:,:)   => null() !internode height, [m]
  real(r8), pointer :: SURFBs1(:,:,:,:)  => null() !stem surface area, [m2 d-2]
  real(r8), pointer :: ARLFLs1(:,:,:,:)  => null() !layer leaf area, [m2 d-2]
  real(r8), pointer :: ARSTKs1(:,:,:)    => null() !stem layer area, [m2 d-2]
  real(r8), pointer :: CFs1(:)           => null() !clumping factor for self-shading in canopy layer, [-]
  real(r8), pointer :: XTLIs1(:)         => null() !number of nodes in seed, [-]
  real(r8), pointer :: ZCs1(:)           => null() !canopy height, [m]
  integer,  pointer :: NGs1(:)           => null() !soil layer at planting depth, [-]
  integer,  pointer :: KLEAFs1(:,:)      => null() !leaf number, [-]
  real(r8), pointer :: VSTGs1(:,:)       => null() !leaf number, [-]
  contains
    procedure, public :: Init    => plt_morph_init
    procedure, public :: Destroy => plt_morph_destroy
  end type plant_morph_type

  type, public :: plant_pheno_type
  integer,  pointer :: IDTHs1(:)  => null()     !flag for species death
  integer,  pointer :: IFLGCs1(:) => null()     !flag for living pft
  real(r8), pointer :: RSETCs1(:) => null()     !effect of canopy C status on seed set , []
  real(r8), pointer :: RSETNs1(:) => null()     !effect of canopy N status on seed set , []
  real(r8), pointer :: RSETPs1(:) => null()     !effect of canopy P status on seed set , []
  real(r8), pointer :: XRNIs1(:)  => null()     !rate of node initiation, [h-1 at 25 oC]
  real(r8), pointer :: XRLAs1(:)  => null()     !rate of leaf initiation, [h-1 at 25 oC]
  real(r8), pointer :: XDLs1(:)   => null()     !critical daylength for phenological progress, [h]
  real(r8), pointer :: XPPDs1(:)  => null()     !difference between current and critical daylengths used to calculate  phenological progress [h]
  integer,  pointer :: IDTHBs1(:,:)  => null()  !flag to detect branch death , [-]
  integer,  pointer :: IFLGPs1(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGFs1(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGEs1(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGAs1(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGGs1(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGRs1(:,:)  => null()  !branch phenology flag, [-]
  integer,  pointer :: IFLGQs1(:,:)  => null()  !branch phenology flag, [h]
  integer,  pointer :: KVSTGs1(:,:)  => null()  !leaf growth stage counter, [-]
  integer,  pointer :: IWTYPs1(:)    => null()  !climate signal for phenological progress: none, temperature, water stress
  integer,  pointer :: IBTYPs1(:)    => null()  !phenologically-driven above-ground turnover: all, foliar only, none
  integer,  pointer :: IDTHPs1(:)    => null()  !flag to detect canopy death
  integer,  pointer :: ISTYPs1(:)    => null()  !plant growth habit: annual or perennial
  integer,  pointer :: IDTHRs1(:)    => null()  !flag to detect root system death
  integer,  pointer :: IDTYPs1(:)    => null()  !plant growth habit (determinate or indeterminate)
  integer,  pointer :: IPTYPs1(:)    => null()  !photoperiod type (neutral, long day, short day)
  integer,  pointer :: IGTYPs1(:)    => null()  !plant growth type (vascular, non-vascular)
  integer,  pointer :: IFLGIs1(:)    => null()  !PFT initialization flag:0=no,1=yes
  integer,  pointer :: KVSTGNs1(:,:) => null()  !leaf growth stage counter, [-]
  integer,  pointer :: IDAYs1(:,:,:) => null()  !plant growth stage, [-]
  real(r8), pointer :: TGSTGIs1(:,:) => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: TGSTGFs1(:,:) => null()  !normalized node number during reproductive growth stages , [-]
  real(r8), pointer :: VSTGXs1(:,:)  => null()  !leaf number at floral initiation, [-]
  real(r8), pointer :: VRNYs1(:,:)   => null()  !initial heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: VRNZs1(:,:)   => null()  !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: VRNSs1(:,:)   => null()  !heat requirement for spring leafout/dehardening, [h]
  real(r8), pointer :: VRNLs1(:,:)   => null()  !hours above threshold temperature required for spring leafout/dehardening, [-]
  real(r8), pointer :: VRNFs1(:,:)   => null()  !cold requirement for autumn leafoff/hardening, [h]
  real(r8), pointer :: VRNXs1(:,:)   => null()  !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8), pointer :: ATRPs1(:,:)   => null()  !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8), pointer :: FLGZs1(:,:)   => null()  !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  real(r8), pointer :: DGSTGIs1(:,:) => null()  !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8), pointer :: DGSTGFs1(:,:) => null()  !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8), pointer :: GROUPs1(:,:)  => null()  !plant maturity group, [-]
  real(r8), pointer :: GSTGIs1(:,:)  => null()  !normalized node number during vegetative growth stages , [-]
  real(r8), pointer :: GSTGFs1(:,:)  => null()  !normalized node number during reproductive growth stages, [-]
  real(r8), pointer :: FLG4s1(:,:)   => null()  !flag to detect physiological maturity from  grain fill , [-]

  contains
    procedure, public :: Init    =>  plt_pheno_init
    procedure, public :: Destroy =>  plt_pheno_destroy
  end type plant_pheno_type

  type, public :: plant_soilchem_type
  real(r8), pointer :: THETPMs1(:,:) => null()  !soil air-filled porosity, [m3 m-3]
  real(r8), pointer :: DFGSs1(:,:)   => null()  !coefficient for dissolution - volatilization, []
  real(r8), pointer :: RSCSs1(:)     => null()  !soil hydraulic resistance, [MPa h m-2]
  real(r8), pointer :: ZHSGLs1(:)    => null()  !aqueous NH4 diffusivity, [m2 h-1]
  real(r8), pointer :: BKDSs1(:)     => null()  !soil bulk density, [Mg m-3]
  real(r8), pointer :: CH2P4s1(:)    => null()  !aqueous PO4 concentration non-band	[gP m-3]
  real(r8), pointer :: CH1P4Bs1(:)   => null()  !aqueous H1PO4 concentration band [gP m-3]
  real(r8), pointer :: CH2P4Bs1(:)   => null()  !aqueous PO4 concentration band	[gP m-3]
  real(r8), pointer :: CNO3Ss1(:)    => null()  !NO3 concentration non-band micropore	[gN m-3]
  real(r8), pointer :: CH1P4s1(:)    => null()  !aqueous H1PO4 concentration non-band [gP m-3]
  real(r8), pointer :: CNO3Bs1(:)    => null()  !NO3 concentration band micropore	[gN m-3]
  real(r8), pointer :: CNH4Ss1(:)    => null()  !NH4 concentration non-band micropore	[gN m-3]
  real(r8), pointer :: CNH4Bs1(:)    => null()  !NH4 concentration band micropore	[gN m-3]
  real(r8), pointer :: CO2Gs1(:)     => null()  !gaseous CO2	[g d-2]
  real(r8), pointer :: CO2Ss1(:)     => null()  !aqueous CO2  micropore	[gC d-2]
  real(r8), pointer :: CH4Ss1(:)     => null()  !aqueous CH4  micropore	[gC d-2]
  real(r8), pointer :: CCH4Ss1(:)    => null()  !aqueous CH4 concentration micropore	[gC m-3]
  real(r8), pointer :: CZ2OSs1(:)    => null()  !aqueous N2O concentration micropore	[gN m-3]
  real(r8), pointer :: CCH4Gs1(:)    => null()  !gaseous CH4 concentration	[gC m-3]
  real(r8), pointer :: CNH3Ss1(:)    => null()  !NH3 concentration non-band micropore	[gN m-3]
  real(r8), pointer :: CNH3Bs1(:)    => null()  !NH3 concentration band micropore	[gN m-3]
  real(r8), pointer :: CH2GSs1(:)    => null()  !aqueous H2 concentration	[g m-3]
  real(r8), pointer :: HLSGLs1(:)    => null()  !aqueous H2 diffusivity, [m2 h-1]
  real(r8), pointer :: CLSGLs1(:)    => null()  !aqueous CO2 diffusivity	[m2 h-1]
  real(r8), pointer :: CQSGLs1(:)    => null()  !aqueous CH4 diffusivity	[m2 h-1]
  real(r8), pointer :: CZ2OGs1(:)    => null()  !gaseous N2O concentration	[gN m-3]
  real(r8), pointer :: CNH3Gs1(:)    => null()  !gaseous NH3 concentration	[gN m-3]
  real(r8), pointer :: CH2GGs1(:)    => null()  !gaseous H2 concentration	[g m-3]
  real(r8), pointer :: CORGCs1(:)    => null()  !soil organic C content [gC kg soil-1]
  real(r8), pointer :: CPO4Ss1(:)    => null()  !PO4 concentration non-band micropore	[g m-3]
  real(r8), pointer :: CNDUs1(:)     => null()  !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8), pointer :: CGSGLs1(:)    => null()  !gaseous CO2 diffusivity	[m2 h-1]
  real(r8), pointer :: CHSGLs1(:)    => null()  !gaseous CH4 diffusivity	[m2 h-1]
  real(r8), pointer :: HGSGLs1(:)    => null()  !gaseous H2 diffusivity  [m2 h-1]
  real(r8), pointer :: OGSGLs1(:)    => null()  !gaseous O2 diffusivity	[m2 h-1]
  real(r8), pointer :: Z2SGLs1(:)    => null()  !gaseous N2O diffusivity, [m2 h-1]
  real(r8), pointer :: ZOSGLs1(:)    => null()  !aqueous NO3 diffusivity, [m2 h-1]
  real(r8), pointer :: H1PO4s1(:)    => null()  !soil aqueous HPO4 content micropore non-band, [gP d-2]
  real(r8), pointer :: H2PO4s1(:)    => null()  !PO4 non-band micropore, [gP d-2]
  real(r8), pointer :: H1POBs1(:)    => null()  !soil aqueous HPO4 content micropore band, [gP d-2]
  real(r8), pointer :: H2POBs1(:)    => null()  !PO4 band micropore, [gP d-2]
  real(r8), pointer :: H2GSs1(:)     => null()  !aqueous H2 	[g d-2]
  real(r8), pointer :: OXYGs1(:)     => null()  !gaseous O2 	[g d-2]
  real(r8), pointer :: OXYSs1(:)     => null()  !aqueous O2  micropore	[g d-2]
  real(r8), pointer :: OLSGLs1(:)    => null()  !aqueous CO2 diffusivity	[m2 h-1]
  real(r8), pointer :: POSGLs1(:)    => null()  !aqueous PO4 diffusivity, [m2 h-1]
  real(r8), pointer :: PSISTs1(:)    => null()  !soil micropore total water potential [MPa]
  real(r8), pointer :: SCO2Ls1(:)    => null()  !solubility of CO2, [m3 m-3]
  real(r8), pointer :: SOXYLs1(:)    => null()  !solubility of O2, [m3 m-3]
  real(r8), pointer :: SCH4Ls1(:)    => null()  !solubility of CH4, [m3 m-3]
  real(r8), pointer :: SN2OLs1(:)    => null()  !solubility of N2O, [m3 m-3]
  real(r8), pointer :: SNH3Ls1(:)    => null()  !solubility of NH3, [m3 m-3]
  real(r8), pointer :: SH2GLs1(:)    => null()  !solubility of H2, [m3 m-3]
  real(r8), pointer :: THETWs1(:)    => null()  !volumetric water content [m3 m-3]
  real(r8), pointer :: THETYs1(:)    => null()  !air-dry water content, [m3 m-3]
  real(r8), pointer :: VOLXs1(:)     => null()  !volume of soil layer	m3 d-2
  real(r8), pointer :: VLPOBs1(:)    => null()  !PO4 band volume fracrion, [0-1]
  real(r8), pointer :: VLNO3s1(:)    => null()  !NO3 non-band volume fracrion, []
  real(r8), pointer :: VLPO4s1(:)    => null()  !PO4 non-band volume fracrion, []
  real(r8), pointer :: VLNOBs1(:)    => null()  !NO3 band volume fraction, []
  real(r8), pointer :: VLNH4s1(:)    => null()  !NH4 non-band volume fraction, []
  real(r8), pointer :: VLNHBs1(:)    => null()  !NH4 band volume fraction, []
  real(r8), pointer :: VOLYs1(:)     => null()  !total micropore volume [m3 d-2]
  real(r8), pointer :: VOLIs1(:)     => null()  !soil micropore ice content   [m3 d-2]
  real(r8), pointer :: VOLWs1(:)     => null()  !soil micropore water content [m3 d-2]
  real(r8), pointer :: VOLAs1(:)     => null()  !total volume in micropores [m3 d-2]
  real(r8), pointer :: ZNO3Ss1(:)    => null()  !NO3 non-band micropore, [gN d-2]
  real(r8), pointer :: ZNO3Bs1(:)    => null()  !NO3 band micropore, [Ng d-2]
  real(r8), pointer :: ZNSGLs1(:)    => null()  !aqueous NH3 diffusivity, [m2 h-1]
  real(r8), pointer :: ZVSGLs1(:)    => null()  !aqueous N2O diffusivity, [m2 h-1]
  real(r8), pointer :: ZNH4Ss1(:)    => null()  !NH4 non-band micropore, [gN d-2]
  real(r8), pointer :: ZNH4Bs1(:)    => null()  !NH4 band micropore, [gN d-2]
  real(r8), pointer :: Z2OSs1(:)     => null()  !aqueous N2O micropore, [gN d-2]
  real(r8), pointer :: ZNH3Ss1(:)    => null()  !NH3 non-band micropore, [gN d-2]
  real(r8), pointer :: ZNH3Bs1(:)    => null()  !NH3 band micropore, [g d-2]
  real(r8), pointer :: OQCs1(:,:)    => null()  !dissolved organic C micropore	[gC d-2]
  real(r8), pointer :: OQNs1(:,:)    => null()  !dissolved organic N micropore	[gN d-2]
  real(r8), pointer :: OQPs1(:,:)    => null()  !dissolved organic P micropore	[gP d-2]

  contains
    procedure, public :: Init => plt_soilchem_init
    procedure, public :: Destroy => plt_soilchem_destroy
  end type plant_soilchem_type

  type, public :: plant_allometry_type
  real(r8), pointer :: CPRTSs1(:)    => null()  !root P:C ratio x root growth yield, [-]
  real(r8), pointer :: CNRTSs1(:)    => null()  !root N:C ratio x root growth yield, [-]
  real(r8), pointer :: CNNDs1(:)     => null()  !nodule N:C ratio, [gN gC-1]
  real(r8), pointer :: CPNDs1(:)     => null()  !nodule P:C ratio, [gP gC-1]
  real(r8), pointer :: CNRTs1(:)     => null()  !root N:C ratio, [gN gC-1]
  real(r8), pointer :: CPRTs1(:)     => null()  !root P:C ratio, [gP gC-1]
  real(r8), pointer :: CNWSs1(:)     => null()  !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: CPWSs1(:)     => null()  !C:P ratio in remobilizable nonstructural biomass, [-]
  real(r8), pointer :: DMNDs1(:)     => null()  !nodule growth yield, [g g-1]
  real(r8), pointer :: DMRTs1(:)     => null()  !root growth yield, [g g-1]
  real(r8), pointer :: CPEARs1(:)    => null()  !ear P:C ratio, [gP gC-1]
  real(r8), pointer :: DMSHEs1(:)    => null()  !sheath growth yield, [g g-1]
  real(r8), pointer :: CPHSKs1(:)    => null()  !husk P:C ratio, [gP gC-1]
  real(r8), pointer :: DMSTKs1(:)    => null()  !stalk growth yield, [gC gC-1]
  real(r8), pointer :: DMHSKs1(:)    => null()  !husk growth yield, [gC gC-1]
  real(r8), pointer :: DMRSVs1(:)    => null()  !reserve growth yield, [gC gC-1]
  real(r8), pointer :: DMGRs1(:)     => null()  !grain growth yield, [gC gC-1]
  real(r8), pointer :: DMEARs1(:)    => null()  !ear growth yield, [gC gC-1]
  real(r8), pointer :: CNHSKs1(:)    => null()  !husk N:C ratio, [gN gC-1]
  real(r8), pointer :: CNRSVs1(:)    => null()  !reserve N:C ratio, [gN gC-1]
  real(r8), pointer :: CNEARs1(:)    => null()  !ear N:C ratio, [gN gC-1]
  real(r8), pointer :: CPRSVs1(:)    => null()  !reserve P:C ratio, [gP gC-1]
  real(r8), pointer :: CPGRs1(:)     => null()  !grain P:C ratio, [gP gP-1]
  real(r8), pointer :: CNSTKs1(:)    => null()  !stalk N:C ratio, [gN gC-1]
  real(r8), pointer :: FWODBs1(:)    => null()  !woody C allocation
  real(r8), pointer :: FWODLNs1(:)   => null()  !leaf N allocation
  real(r8), pointer :: FWODLPs1(:)   => null()  !leaf P allocation
  real(r8), pointer :: FWODSPs1(:)   => null()  !P woody fraction in petiole
  real(r8), pointer :: FWODSNs1(:)   => null()  !N woody fraction in petiole
  real(r8), pointer :: FWODRNs1(:)   => null()  !N woody fraction in root
  real(r8), pointer :: FWODRPs1(:)   => null()  !P woody fraction in root
  real(r8), pointer :: FWODRs1(:)    => null()  !C woody fraction in root
  real(r8), pointer :: FWOODs1(:)    => null()  !woody C allocation
  real(r8), pointer :: FVRNs1(:)     => null()  !allocation parameter
  real(r8), pointer :: FWOODNs1(:)   => null()  !woody N allocation
  real(r8), pointer :: FWOODPs1(:)   => null()  !woody P allocation
  real(r8), pointer :: DMLFs1(:)     => null()  !leaf growth yield, [g g-1]
  real(r8), pointer :: CNGRs1(:)     => null()  !grain N:C ratio, [g g-1]
  real(r8), pointer :: CPLFs1(:)     => null()  !maximum leaf P:C ratio, [g g-1]
  real(r8), pointer :: CPSHEs1(:)    => null()  !sheath P:C ratio, [g g-1]
  real(r8), pointer :: CNSHEs1(:)    => null()  !sheath N:C ratio, [g g-1]
  real(r8), pointer :: CPSTKs1(:)    => null()  !stalk P:C ratio, [g g-1]
  real(r8), pointer :: CNLFs1(:)     => null()  !maximum leaf N:C ratio, [g g-1]
  contains
    procedure, public :: Init => plt_allom_init
    procedure, public :: Destroy => plt_allom_destroy
  end type plant_allometry_type

  type, public :: plant_biom_type
  real(r8) :: WTSTGTs1                               !total standing dead C, [g d-2]
  real(r8), pointer :: PPOOLPs1(:)      => null()    !canopy nonstructural P, [gP d-2]
  real(r8), pointer :: PPOLNPs1(:)      => null()    !canopy nonstructural P concentration, [gP gC-1]
  real(r8), pointer :: CCPOLPs1(:)      => null()    !canopy nonstructural C concentration, [gC d-2]
  real(r8), pointer :: CPOOLPs1(:)      => null()    !canopy nonstructural P concentration, [gP d-2]
  real(r8), pointer :: CPOLNPs1(:)      => null()    !canopy nodule nonstructural P, [gP d-2]
  real(r8), pointer :: CCPLNPs1(:)      => null()    !nodule nonstructural C, [gC d-2]
  real(r8), pointer :: CZPOLPs1(:)      => null()    !canopy nonstructural N concentration, [g g-1]
  real(r8), pointer :: CPPOLPs1(:)      => null()    !canopy nonstructural P concentration, [g g-1]
  real(r8), pointer :: PPOOLRs1(:,:,:)  => null()    !root layer nonstructural P, [g d-2]
  real(r8), pointer :: WTRTLs1(:,:,:)   => null()    !root layer structural C, [g d-2]
  real(r8), pointer :: WTRTDs1(:,:,:)   => null()    !root layer C, [g d-2]
  real(r8), pointer :: WSRTLs1(:,:,:)   => null()    !root layer protein C, [g d-2]
  real(r8), pointer :: CWSRTLs1(:,:,:)  => null()    !root layer protein C concentration, [g g-1]
  real(r8), pointer :: ZPOOLRs1(:,:,:)  => null()    !root layer nonstructural N, [g d-2]
  real(r8), pointer :: CPOOLRs1(:,:,:)  => null()    !root  layer nonstructural C, [g d-2]
  real(r8), pointer :: CCPOLRs1(:,:,:)  => null()    !root  layer nonstructural C concentration, [g g-1]
  real(r8), pointer :: CZPOLRs1(:,:,:)  => null()    !root layer nonstructural N concentration, [g g-1]
  real(r8), pointer :: CPPOLRs1(:,:,:)  => null()    !root layer nonstructural P concentration, [g g-1]
  real(r8), pointer :: ZPOLNPs1(:)      => null()    !canopy nonstructural P concentration, [g g-1]
  real(r8), pointer :: CCPOLBs1(:,:)    => null()    !branch nonstructural C concentration, [g d-2]
  real(r8), pointer :: CZPOLBs1(:,:)    => null()    !branch nonstructural N concentration, [g g-1]
  real(r8), pointer :: CPPOLBs1(:,:)    => null()    !branch nonstructural P concentration, [g g-1]
  real(r8), pointer :: PPOLNBs1(:,:)    => null()    !branch nonstructural P concentration, [g g-1]
  real(r8), pointer :: WGNODEs1(:,:,:)  => null()    !internode C, [g d-2]
  real(r8), pointer :: WGNODNs1(:,:,:)  => null()    !internode N, [g d-2]
  real(r8), pointer :: WGNODPs1(:,:,:)  => null()    !nodule P, [g d-2]
  real(r8), pointer :: WGLFs1(:,:,:)    => null()    !leaf C, [g d-2]
  real(r8), pointer :: WGLFNs1(:,:,:)   => null()    !leaf N, [g d-2]
  real(r8), pointer :: WGLFPs1(:,:,:)   => null()    !leaf P, [g d-2]
  real(r8), pointer :: WSLFs1(:,:,:)    => null()    !layer leaf protein C, [g d-2]
  real(r8), pointer :: WGSHEs1(:,:,:)   => null()    !sheath C , [g d-2]
  real(r8), pointer :: WGSHNs1(:,:,:)   => null()    !sheath N, [g d-2]
  real(r8), pointer :: WGSHPs1(:,:,:)   => null()    !sheath P, [g d-2]
  real(r8), pointer :: WSSHEs1(:,:,:)   => null()    !layer sheath protein C, [g d-2]
  real(r8), pointer :: WGLFLs1(:,:,:,:) => null()    !layer leaf C, [g d-2]
  real(r8), pointer :: WGLFLNs1(:,:,:,:)=> null()    !layer leaf N, [g d-2]
  real(r8), pointer :: WGLFLPs1(:,:,:,:)=> null()    !leaf layer P, [g d-2]
  real(r8), pointer :: ZPOLNBs1(:,:) => null()  !branch nonstructural N concentration, [g g-1]
  real(r8), pointer :: WGLFTs1(:)    => null()  !total leaf mass, [gC d-2]
  real(r8), pointer :: WTSTDIs1(:)   => null()  !initial standing dead C, [g C m-2]
  real(r8), pointer :: WTRTs1(:)     => null()  !plant root C, [gC d-2]
  real(r8), pointer :: WTRTSs1(:)    => null()  !plant root structural C, [gC d-2]
  real(r8), pointer :: WTRTSNs1(:)   => null()  !plant root structural N, [gN d-2]
  real(r8), pointer :: WTRTSPs1(:)   => null()  !plant root structural P, [gP d-2]
  real(r8), pointer :: WTRVXs1(:)    => null()  !plant stored nonstructural C at planting, [gC d-2]
  real(r8), pointer :: WTRVCs1(:)    => null()  !plant stored nonstructural C, [gC d-2]
  real(r8), pointer :: WTLSs1(:)     => null()  !canopy leaf + sheath C, [g d-2]
  real(r8), pointer :: WTSHTs1(:)    => null()  !canopy shoot C, [g d-2]
  real(r8), pointer :: WTSHTAs1(:)   => null()  !landscape average canopy shoot C, [g d-2]
  real(r8), pointer :: WTSTGs1(:)    => null()  !standing dead C, [g d-2]
  real(r8), pointer :: WTSTGNs1(:)   => null()  !standing dead N, [g d-2]
  real(r8), pointer :: WTSTGPs1(:)   => null()  !standing dead P, [g d-2]
  real(r8), pointer :: WTNDs1(:)     => null()  !root total nodule mass, [g d-2]
  real(r8), pointer :: WTNDPs1(:)    => null()  !total nodule P, [g d-2]
  real(r8), pointer :: WTRTts1(:)    => null()  !plant root C, [g d-2]
  real(r8), pointer :: WTRTNs1(:)    => null()  !total root N , [g d-2]
  real(r8), pointer :: WTRTPs1(:)    => null()  !root total P, [g d-2]
  real(r8), pointer :: WTNDNs1(:)    => null()  !total canopy nodule N, [g d-2]
  real(r8), pointer :: WTSHNs1(:)    => null()  !canopy  N, [g d-2]
  real(r8), pointer :: WTSHPs1(:)    => null()  !canopy total P, [g d-2]
  real(r8), pointer :: WTRVNs1(:)    => null()  !plant stored nonstructural N, [g d-2]
  real(r8), pointer :: WTRVPs1(:)    => null()  !plant stored nonstructural P, [g d-2]
  real(r8), pointer :: ZPOOLPs1(:)   => null()  !canopy  nonstructural N, [gN d-2]
  real(r8), pointer :: PPOOLs1(:,:)  => null()  !branch nonstructural P, [g d-2]
  real(r8), pointer :: CPOOLs1(:,:)  => null()  !branch nonstructural C, [g d-2]
  real(r8), pointer :: ZPOOLs1(:,:)  => null()  !branch  nonstructural N, [g d-2]
  real(r8), pointer :: CPOLNBs1(:,:) => null()  !branch nodule nonstructural C, [g d-2]
  real(r8), pointer :: WTLSBs1(:,:)  => null()  !branch leaf + sheath C, [g d-2]
  real(r8), pointer :: WTRSVBs1(:,:) => null()  !branch reserve C, [g d-2]
  real(r8), pointer :: WTLFBs1(:,:)  => null()   !branch leaf C, [g d-2]
  real(r8), pointer :: WTNDBs1(:,:)  => null()   !branch nodule C, [g d-2]
  real(r8), pointer :: WTSHEBs1(:,:) => null()   !branch sheath C , [g d-2]
  real(r8), pointer :: WTEARBs1(:,:) => null()   !branch ear C, [g d-2]
  real(r8), pointer :: WTHSKBs1(:,:) => null()   !branch husk C, [g d-2]
  real(r8), pointer :: WTRSBNs1(:,:) => null()   !branch reserve N, [g d-2]
  real(r8), pointer :: WTLFBNs1(:,:) => null()   !branch leaf N, [g d-2]
  real(r8), pointer :: WTNDBNs1(:,:) => null()   !branch nodule N, [g d-2]
  real(r8), pointer :: WTSHBNs1(:,:) => null()   !branch sheath N, [g d-2]
  real(r8), pointer :: WTEABNs1(:,:) => null()   !branch ear N, [g d-2]
  real(r8), pointer :: WTHSBNs1(:,:) => null()   !branch husk N, [g d-2]
  real(r8), pointer :: WTRSBPs1(:,:) => null()   !branch reserve P, [g d-2]
  real(r8), pointer :: WTLFBPs1(:,:) => null()   !branch leaf P, [g d-2]
  real(r8), pointer :: WTNDBPs1(:,:) => null()   !branch nodule P, [g d-2]
  real(r8), pointer :: WTSHBPs1(:,:) => null()   !branch sheath P, [g d-2]
  real(r8), pointer :: WTEABPs1(:,:) => null()   !branch ear P, [g d-2]
  real(r8), pointer :: WTHSBPs1(:,:) => null()   !branch husk P, [g d-2]
  real(r8), pointer :: WTGRBs1(:,:)  => null()   !branch grain C, [g d-2]
  real(r8), pointer :: WTGRBNs1(:,:) => null()   !branch grain N, [g d-2]
  real(r8), pointer :: WTGRBPs1(:,:) => null()   !branch grain P, [g d-2]
  real(r8), pointer :: WTSTKBs1(:,:) => null()   !branch stalk C, [g d-2]
  real(r8), pointer :: WTSTBNs1(:,:) => null()   !branch stalk N, [g d-2]
  real(r8), pointer :: WTSTBPs1(:,:) => null()   !branch stalk P, [g d-2]
  real(r8), pointer :: WTSHTBs1(:,:) => null()   !branch shoot C, [g d-2]
  real(r8), pointer :: WTSHTNs1(:,:) => null()   !branch N, [g d-2]
  real(r8), pointer :: WTSHTPs1(:,:) => null()   !branch total P, [g d-2]
  real(r8), pointer :: WGLFXs1(:,:)  => null()   !branch leaf structural C, [g d-2]
  real(r8), pointer :: WGLFNXs1(:,:) => null()   !branch leaf structural N, [g d-2]
  real(r8), pointer :: WGLFPXs1(:,:) => null()   !branch leaf structural P, [g d-2]
  real(r8), pointer :: WGSHEXs1(:,:) => null()   !branch sheath structural C, [g d-2]
  real(r8), pointer :: WGSHNXs1(:,:) => null()   !branch sheath structural N, [g d-2]
  real(r8), pointer :: WGSHPXs1(:,:) => null()   !branch sheath structural P, [g d-2]
  real(r8), pointer :: WTSTXBs1(:,:) => null()   !branch stalk structural C, [g d-2]
  real(r8), pointer :: WTSTXNs1(:,:) => null()   !branch stalk structural N, [g d-2]
  real(r8), pointer :: WTSTXPs1(:,:) => null()   !branch stalk structural P, [g d-2]
  real(r8), pointer :: WVSTKBs1(:,:) => null()   !branch active stalk C, [g d-2]
  real(r8), pointer :: WTSTKs1(:)    => null()   !canopy stalk C, [g d-2]
  real(r8), pointer :: WTGRNNs1(:)   => null()   !canopy grain N, [g d-2]
  real(r8), pointer :: WTSHEPs1(:)   => null()   !canopy sheath P, [g d-2]
  real(r8), pointer :: WVSTKs1(:)    => null()   !canopy active stalk C, [g d-2
  real(r8), pointer :: WTLFPs1(:)    => null()   !canopy leaf P, [g d-2]
  real(r8), pointer :: WTLFs1(:)     => null()   !canopy leaf C, [g d-2]
  real(r8), pointer :: WTSHEs1(:)    => null()   !canopy sheath C , [g d-2]
  real(r8), pointer :: WTRSVs1(:)    => null()   !canopy reserve C, [g d-2]
  real(r8), pointer :: WTEARNs1(:)   => null()   !canopy ear N, [g d-2]
  real(r8), pointer :: WTLFNs1(:)    => null()   !canopy leaf N, [g d-2]
  real(r8), pointer :: WTRSVNs1(:)   => null()   !canopy reserve N, [g d-2]
  real(r8), pointer :: WTSHENs1(:)   => null()   !canopy sheath N, [g d-2]
  real(r8), pointer :: WTHSKNs1(:)   => null()   !canopy husk N, [g d-2]
  real(r8), pointer :: WTHSKs1(:)    => null()   !canopy husk C, [g d-2]
  real(r8), pointer :: WTRTAs1(:)    => null()   !root C per plant, [g p-1]
  real(r8), pointer :: WTGRs1(:)     => null()   !canopy grain C, [g d-2]
  real(r8), pointer :: WTEARs1(:)    => null()   !canopy ear C, [g d-2]
  real(r8), pointer :: WTSTKPs1(:)   => null()   !canopy stalk P, [g d-2]
  real(r8), pointer :: WTRSVPs1(:)   => null()   !canopy reserve P, [g d-2]
  real(r8), pointer :: WTEARPs1(:)   => null()   !canopy ear C, [g d-2]
  real(r8), pointer :: WTHSKPs1(:)   => null()   !canopy husk P, [g d-2]
  real(r8), pointer :: WTSTKNs1(:)   => null()   !canopy stalk N, [g d-2]
  real(r8), pointer :: WTGRNPs1(:)   => null()   !canopy grain P, [g d-2]
  contains
    procedure, public :: Init => plt_biom_init
    procedure, public :: Destroy => plt_biom_destroy
  end type plant_biom_type

  type(plant_allometry_type), public :: plt_allom           !plant allometric parameters
  type(plant_biom_type), public, target :: plt_biom         !plant biomass variables
  type(plant_soilchem_type), public, target :: plt_soilchem !soil bgc interface with plant root
  type(plant_pheno_type), public, target :: plt_pheno       !plant phenology
  type(plant_morph_type), public, target :: plt_morph       !plant morphology
  type(plant_radiation_type), public, target :: plt_rad     !plant radiation type
  type(plant_photosyns_type), public, target :: plt_photo   !plant photosynthesis type

  contains
!----------------------------------------------------------------------

  subroutine plt_allom_init(this)
  implicit none
  class(plant_allometry_type) :: this

  allocate(this%DMNDs1(JP1))
  allocate(this%DMRTs1(JP1))
  allocate(this%CPRTs1(JP1))
  allocate(this%CNWSs1(JP1))
  allocate(this%CPWSs1(JP1))
  allocate(this%CPRTSs1(JP1))
  allocate(this%CNRTSs1(JP1))
  allocate(this%CNNDs1(JP1))
  allocate(this%CPNDs1(JP1))
  allocate(this%CNRTs1(JP1))
  allocate(this%CPLFs1(JP1))
  allocate(this%CPSHEs1(JP1))
  allocate(this%CNLFs1(JP1))
  allocate(this%CNSHEs1(JP1))
  allocate(this%CNGRs1(JP1))
  allocate(this%CPSTKs1(JP1))
  allocate(this%CNSTKs1(JP1))
  allocate(this%CPGRs1(JP1))
  allocate(this%CPEARs1(JP1))
  allocate(this%CPRSVs1(JP1))
  allocate(this%CNRSVs1(JP1))
  allocate(this%CPHSKs1(JP1))
  allocate(this%FVRNs1(0:5))
  allocate(this%FWOODs1(0:1))
  allocate(this%FWODLPs1(0:1))
  allocate(this%FWODLNs1(0:1))
  allocate(this%FWODRs1(0:1))
  allocate(this%FWODSPs1(0:1))
  allocate(this%FWODBs1(0:1))
  allocate(this%FWODSNs1(0:1))
  allocate(this%FWODRNs1(0:1))
  allocate(this%FWODRPs1(0:1))
  allocate(this%FWOODNs1(0:1))
  allocate(this%FWOODPs1(0:1))

  allocate(this%DMSHEs1(JP1))
  allocate(this%DMHSKs1(JP1))
  allocate(this%DMSTKs1(JP1))
  allocate(this%DMRSVs1(JP1))
  allocate(this%DMEARs1(JP1))
  allocate(this%DMGRs1(JP1))
  allocate(this%CNHSKs1(JP1))
  allocate(this%CNEARs1(JP1))
  allocate(this%DMLFs1(JP1))

  end subroutine plt_allom_init
!----------------------------------------------------------------------

  subroutine plt_allom_destroy(this)
  implicit none

  class(plant_allometry_type) :: this

!  if(allocated(DMNDs1))deallocate(DMNDs1)
!  if(allocated(DMRTs1))deallocate(DMRTs1)
!  if(allocated(CPRTs1))deallocate(CPRTs1)
!  if(allocated(CNWSs1))deallocate(CNWSs1)
!  if(allocated(CPWSs1))deallocate(CPWSs1)
!  if(allocated(CPRTSs1))deallocate(CPRTSs1)
!  if(allocated(CNRTSs1))deallocate(CNRTSs1)
!  if(allocated(CNNDs1))deallocate(CNNDs1)
!  if(allocated(CPNDs1))deallocate(CPNDs1)
!  if(allocated(CNRTs1))deallocate(CNRTs1)
!  if(allocated(CWSRTs1))deallocate(CWSRTs1)
!  if(allocated(CPLFs1))deallocate(CPLFs1)
!  if(allocated(CPSHEs1))deallocate(CPSHEs1)
!  if(allocated(CNSHEs1))deallocate(CNSHEs1)
!  if(allocated(CNLFs1))deallocate(CNLFs1)
!  if(allocated(DMLFs1))deallocate(DMLFs1)
!  if(allocated(CPSTKs1))deallocate(CPSTKs1)
!  if(allocated(CNGRs1))deallocate(CNGRs1)
!  if(allocated(CNSTKs1))deallocate(CNSTKs1)
!  if(allocated(CPGRs1))deallocate(CPGRs1)
!  if(allocated(DMSHEs1))deallocate(DMSHEs1)
!  if(allocated(DMSTKs1))deallocate(DMSTKs1)
!  if(allocated(DMGRs1))deallocate(DMGRs1)
!  if(allocated(DMRSVs1))deallocate(DMRSVs1)
!  if(allocated(DMEARs1))deallocate(DMEARs1)
!  if(allocated(DMHSKs1))deallocate(DMHSKs1)
!  if(allocated(FVRNs1))deallocate(FVRNs1)
!  if(allocated(FWODLPs1))deallocate(FWODLPs1)
!  if(allocated(FWODLNs1))deallocate(FWODLNs1)
!  if(allocated(FWOODs1))deallocate(FWOODs1)
!  if(allocated(FWOODPs1))deallocate(FWOODPs1)
!  if(allocated(FWOODNs1))deallocate(FWOODNs1)
!  if(allocated(FWODRs1))deallocate(FWODRs1)
!  if(allocated(FWODRNs1))deallocate(FWODRNs1)
!  if(allocated(FWODSNs1))deallocate(FWODSNs1)
!  if(allocated(FWODSPs1))deallocate(FWODSPs1)
!  if(allocated(FWODBs1))deallocate(FWODBs1)
!  if(allocated(FWODRPs1))deallocate(FWODRPs1)
!  if(allocated(CPEARs1))deallocate(CPEARs1)
!  if(allocated(CPHSKs1))deallocate(CPHSKs1)
!  if(allocated(CNHSKs1))deallocate(CNHSKs1)
!  if(allocated(CNEARs1))deallocate(CNEARs1)
!  if(allocated(CNRSVs1))deallocate(CNRSVs1)
!  if(allocated(CPRSVs1))deallocate(CPRSVs1)

  end subroutine plt_allom_destroy
!----------------------------------------------------------------------
  subroutine plt_biom_init(this)
  implicit none
  class(plant_biom_type) :: this

  allocate(this%CCPOLPs1(JP1))
  allocate(this%CPOOLPs1(JP1))
  allocate(this%CPOLNPs1(JP1))
  allocate(this%CCPLNPs1(JP1))
  allocate(this%CZPOLPs1(JP1))
  allocate(this%CPPOLPs1(JP1))
  allocate(this%PPOOLRs1(2,JZ1,JP1))
  allocate(this%CWSRTLs1(2,JZ1,JP1))
  allocate(this%WSRTLs1(2,JZ1,JP1))
  allocate(this%ZPOOLRs1(2,JZ1,JP1))
  allocate(this%WTRTLs1(2,JZ1,JP1))
  allocate(this%WTRTDs1(2,JZ1,JP1))
  allocate(this%CPOOLRs1(2,JZ1,JP1))
  allocate(this%CCPOLRs1(2,JZ1,JP1))
  allocate(this%CZPOLRs1(2,JZ1,JP1))
  allocate(this%CPPOLRs1(2,JZ1,JP1))
  allocate(this%CPOOLs1(JC1,JP1))
  allocate(this%ZPOOLPs1(JP1))
  allocate(this%ZPOLNPs1(JP1))
  allocate(this%WVSTKs1(JP1))
  allocate(this%WSLFs1(0:JNODS1,JC1,JP1))
  allocate(this%WSSHEs1(0:JNODS1,JC1,JP1))
  allocate(this%WGNODEs1(0:JNODS1,JC1,JP1))
  allocate(this%WGNODNs1(0:JNODS1,JC1,JP1))
  allocate(this%WGNODPs1(0:JNODS1,JC1,JP1))
  allocate(this%WGLFs1(0:JNODS1,JC1,JP1))
  allocate(this%WGLFNs1(0:JNODS1,JC1,JP1))
  allocate(this%WGLFPs1(0:JNODS1,JC1,JP1))
  allocate(this%WGSHEs1(0:JNODS1,JC1,JP1))
  allocate(this%WGSHNs1(0:JNODS1,JC1,JP1))
  allocate(this%WGSHPs1(0:JNODS1,JC1,JP1))
  allocate(this%WGLFLs1(JC1,0:JNODS1,JC1,JP1))
  allocate(this%WGLFLNs1(JC1,0:JNODS1,JC1,JP1))
  allocate(this%WGLFLPs1(JC1,0:JNODS1,JC1,JP1))
  allocate(this%WVSTKBs1(JC1,JP1))
  allocate(this%ZPOOLs1(JC1,JP1))
  allocate(this%ZPOLNBs1(JC1,JP1))
  allocate(this%PPOOLs1(JC1,JP1))
  allocate(this%PPOLNBs1(JC1,JP1))
  allocate(this%CPOLNBs1(JC1,JP1))
  allocate(this%CCPOLBs1(JC1,JP1))
  allocate(this%CZPOLBs1(JC1,JP1))
  allocate(this%CPPOLBs1(JC1,JP1))
  allocate(this%WTRTSs1(JP1))
  allocate(this%WTRTSNs1(JP1))
  allocate(this%WTRTSPs1(JP1))
  allocate(this%WGLFTs1(JC1))
  allocate(this%WTRTs1(JP1))
  allocate(this%WTRVXs1(JP1))
  allocate(this%WTRVCs1(JP1))
  allocate(this%WTLSs1(JP1))
  allocate(this%WTSTGs1(JP1))
  allocate(this%WTSTKPs1(JP1))
  allocate(this%WTRTAs1(JP1))
  allocate(this%WTSHEs1(JP1))
  allocate(this%WTSTKs1(JP1))
  allocate(this%WTEARNs1(JP1))
  allocate(this%WTRSVs1(JP1))
  allocate(this%WTGRs1(JP1))
  allocate(this%WTLFNs1(JP1))
  allocate(this%WTSHENs1(JP1))
  allocate(this%WTRSVNs1(JP1))
  allocate(this%WTHSKNs1(JP1))
  allocate(this%WTHSKs1(JP1))
  allocate(this%WTEARs1(JP1))
  allocate(this%WTGRNNs1(JP1))
  allocate(this%WTRSVPs1(JP1))
  allocate(this%WTHSKPs1(JP1))
  allocate(this%WTSTKNs1(JP1))
  allocate(this%WTEARPs1(JP1))
  allocate(this%WTGRNPs1(JP1))
  allocate(this%WTSTGNs1(JP1))
  allocate(this%WTSTGPs1(JP1))
  allocate(this%WTNDs1(JP1))
  allocate(this%WTRTts1(JP1))
  allocate(this%WTNDNs1(JP1))
  allocate(this%WTRTNs1(JP1))
  allocate(this%WTSHNs1(JP1))
  allocate(this%WTRVNs1(JP1))
  allocate(this%WTNDPs1(JP1))
  allocate(this%WTRTPs1(JP1))
  allocate(this%WTSHPs1(JP1))
  allocate(this%WTRVPs1(JP1))
  allocate(this%WGLFXs1(JC1,JP1))
  allocate(this%WGLFNXs1(JC1,JP1))
  allocate(this%WGLFPXs1(JC1,JP1))
  allocate(this%WGSHEXs1(JC1,JP1))
  allocate(this%WGSHNXs1(JC1,JP1))
  allocate(this%WGSHPXs1(JC1,JP1))
  allocate(this%WTSHTs1(JP1))
  allocate(this%WTSHTAs1(JP1))
  allocate(this%WTSHEPs1(JP1))
  allocate(this%WTLFPs1(JP1))
  allocate(this%WTLFs1(JP1))
  allocate(this%WTSTDIs1(JP1))
  allocate(this%WTLSBs1(JC1,JP1))
  allocate(this%WTRSVBs1(JC1,JP1))
  allocate(this%WTLFBs1(JC1,JP1))
  allocate(this%WTNDBs1(JC1,JP1))
  allocate(this%WTSHEBs1(JC1,JP1))
  allocate(this%WTEARBs1(JC1,JP1))
  allocate(this%WTHSKBs1(JC1,JP1))
  allocate(this%WTRSBNs1(JC1,JP1))
  allocate(this%WTLFBNs1(JC1,JP1))
  allocate(this%WTNDBNs1(JC1,JP1))
  allocate(this%WTSHBNs1(JC1,JP1))
  allocate(this%WTEABNs1(JC1,JP1))
  allocate(this%WTHSBNs1(JC1,JP1))
  allocate(this%WTRSBPs1(JC1,JP1))
  allocate(this%WTLFBPs1(JC1,JP1))
  allocate(this%WTNDBPs1(JC1,JP1))
  allocate(this%WTSHBPs1(JC1,JP1))
  allocate(this%WTEABPs1(JC1,JP1))
  allocate(this%WTHSBPs1(JC1,JP1))
  allocate(this%WTGRBs1(JC1,JP1))
  allocate(this%WTGRBNs1(JC1,JP1))
  allocate(this%WTGRBPs1(JC1,JP1))
  allocate(this%WTSTKBs1(JC1,JP1))
  allocate(this%WTSTBNs1(JC1,JP1))
  allocate(this%WTSTBPs1(JC1,JP1))
  allocate(this%WTSHTBs1(JC1,JP1))
  allocate(this%WTSHTNs1(JC1,JP1))
  allocate(this%WTSHTPs1(JC1,JP1))
  allocate(this%WTSTXBs1(JC1,JP1))
  allocate(this%WTSTXNs1(JC1,JP1))
  allocate(this%WTSTXPs1(JC1,JP1))
  allocate(this%PPOOLPs1(JP1))
  allocate(this%PPOLNPs1(JP1))
  end subroutine plt_biom_init
!----------------------------------------------------------------------

  subroutine plt_biom_destroy(this)

  implicit none
  class(plant_biom_type) :: this


!  if(allocated(PPOOLPs1))deallocate(PPOOLPs1)
!  if(allocated(PPOLNPs1))deallocate(PPOLNPs1)
!  if(allocated(CCPOLPs1))deallocate(CCPOLPs1)
!  if(allocated(CPOOLPs1))deallocate(CPOOLPs1)
!  if(allocated(CPOLNPs1))deallocate(CPOLNPs1)
!  if(allocated(CCPLNPs1))deallocate(CCPLNPs1)
!  if(allocated(CZPOLPs1))deallocate(CZPOLPs1)
!  if(allocated(CPPOLPs1))deallocate(CPPOLPs1)
!  if(allocated(PPOOLRs1))deallocate(PPOOLRs1)
!  if(allocated(CWSRTLs1))deallocate(CWSRTLs1)
!  if(allocated(WSRTLs1))deallocate(WSRTLs1)
!  if(allocated(ZPOOLRs1))deallocate(ZPOOLRs1)
!  if(allocated(WTRTLs1))deallocate(WTRTLs1)
!  if(allocated(WTRTDs1))deallocate(WTRTDs1)
!  if(allocated(CPOOLRs1))deallocate(CPOOLRs1)
!  if(allocated(CCPOLRs1))deallocate(CCPOLRs1)
!  if(allocated(CZPOLRs1))deallocate(CZPOLRs1)
!  if(allocated(CPPOLRs1))deallocate(CPPOLRs1)
!  if(allocated(ZPOLNPs1))deallocate(ZPOLNPs1)
!  if(allocated(WVSTKs1))deallocate(WVSTKs1)
!  if(allocated(WGLFs1))deallocate(WGLFs1)
!  if(allocated(WGLFNs1))deallocate(WGLFNs1)
!  if(allocated(WGLFPs1))deallocate(WGLFPs1)
!  if(allocated(WSLFs1))deallocate(WSLFs1)
!  if(allocated(WSSHEs1))deallocate(WSSHEs1)
!  if(allocated(WGLFLs1))deallocate(WGLFLs1)
!  if(allocated(WGLFLNs1))deallocate(WGLFLNs1)
!  if(allocated(WGLFLPs1))deallocate(WGLFLPs1)
!  if(allocated(WGNODEs1))deallocate(WGNODEs1)
!  if(allocated(WGNODNs1))deallocate(WGNODNs1)
!  if(allocated(WGNODPs1))deallocate(WGNODPs1)
!  if(allocated(WGSHEs1))deallocate(WGSHEs1)
!  if(allocated(WGSHNs1))deallocate(WGSHNs1)
!  if(allocated(WGSHPs1))deallocate(WGSHPs1)
!  if(allocated(WVSTKBs1))deallocate(WVSTKBs1)
!  if(allocated(ZPOOLs1))deallocate(ZPOOLs1)
!  if(allocated(ZPOLNBs1))deallocate(ZPOLNBs1)
!  if(allocated(PPOOLs1))deallocate(PPOOLs1)
!  if(allocated(PPOLNBs1))deallocate(PPOLNBs1)
!  if(allocated(CPOLNBs1))deallocate(CPOLNBs1)
!  if(allocated(CCPOLBs1))deallocate(CCPOLBs1)
!  if(allocated(CZPOLBs1))deallocate(CZPOLBs1)
!  if(allocated(CPPOLBs1))deallocate(CPPOLBs1)
!  if(allocated(WTRTSs1))deallocate(WTRTSs1)
!  if(allocated(WTRTSNs1))deallocate(WTRTSNs1)
!  if(allocated(WTRTSPs1))deallocate(WTRTSPs1)
!  if(allocated(WTRTts1))deallocate(WTRTts1)
!  if(allocated(WGLFTs1))deallocate(WGLFTs1)
!  if(allocated(WTRTs1))deallocate(WTRTs1)
!  if(allocated(ZPOOLPs1))deallocate(ZPOOLPs1)
!  if(allocated(CPOOLs1))deallocate(CPOOLs1)
!  if(allocated(WTRTAs1))deallocate(WTRTAs1)
!  if(allocated(WTLFs1))deallocate(WTLFs1)
!  if(allocated(WTSHEs1))deallocate(WTSHEs1)
!  if(allocated(WTRSVs1))deallocate(WTRSVs1)
!  if(allocated(WTSTKs1))deallocate(WTSTKs1)
!  if(allocated(WTSHEPs1))deallocate(WTSHEPs1)
!  if(allocated(WTLFPs1))deallocate(WTLFPs1)
!  if(allocated(WTEARNs1))deallocate(WTEARNs1)
!  if(allocated(WTGRs1))deallocate(WTGRs1)
!  if(allocated(WTSHENs1))deallocate(WTSHENs1)
!  if(allocated(WTRSVNs1))deallocate(WTRSVNs1)
!  if(allocated(WTLFNs1))deallocate(WTLFNs1)
!  if(allocated(WTHSKNs1))deallocate(WTHSKNs1)
!  if(allocated(WTEARs1))deallocate(WTEARs1)
!  if(allocated(WTHSKs1))deallocate(WTHSKs1)
!  if(allocated(WTGRNNs1))deallocate(WTGRNNs1)
!  if(allocated(WTSTKPs1))deallocate(WTSTKPs1)
!  if(allocated(WTRSVPs1))deallocate(WTRSVPs1)
!  if(allocated(WTHSKPs1))deallocate(WTHSKPs1)
!  if(allocated(WTEARPs1))deallocate(WTEARPs1)
!  if(allocated(WTGRNPs1)) deallocate(WTGRNPs1)
!  if(allocated(WTSTKNs1)) deallocate(WTSTKNs1)
!  if(allocated(WGLFXs1))deallocate(WGLFXs1)
!  if(allocated(WGLFNXs1))deallocate(WGLFNXs1)
!  if(allocated(WGLFPXs1))deallocate(WGLFPXs1)
!  if(allocated(WGSHEXs1))deallocate(WGSHEXs1)
!  if(allocated(WGSHNXs1))deallocate(WGSHNXs1)
!  if(allocated(WGSHPXs1))deallocate(WGSHPXs1)
!  if(allocated(WTLSBs1))deallocate(WTLSBs1)
!  if(allocated(WTRSVBs1))deallocate(WTRSVBs1)
!  if(allocated(WTLFBs1))deallocate(WTLFBs1)
!  if(allocated(WTNDBs1))deallocate(WTNDBs1)
!  if(allocated(WTSHEBs1))deallocate(WTSHEBs1)
!  if(allocated(WTEARBs1))deallocate(WTEARBs1)
!  if(allocated(WTHSKBs1))deallocate(WTHSKBs1)
!  if(allocated(WTRSBNs1))deallocate(WTRSBNs1)
!  if(allocated(WTLFBNs1))deallocate(WTLFBNs1)
!  if(allocated(WTNDBNs1))deallocate(WTNDBNs1)
!  if(allocated(WTSHBNs1))deallocate(WTSHBNs1)
!  if(allocated(WTEABNs1))deallocate(WTEABNs1)
!  if(allocated(WTHSBNs1))deallocate(WTHSBNs1)
!  if(allocated(WTRSBPs1))deallocate(WTRSBPs1)
!  if(allocated(WTLFBPs1))deallocate(WTLFBPs1)
!  if(allocated(WTNDBPs1))deallocate(WTNDBPs1)
!  if(allocated(WTSHBPs1))deallocate(WTSHBPs1)
!  if(allocated(WTEABPs1))deallocate(WTEABPs1)
!  if(allocated(WTHSBPs1))deallocate(WTHSBPs1)
!  if(allocated(WTGRBs1))deallocate(WTGRBs1)
!  if(allocated(WTGRBNs1))deallocate(WTGRBNs1)
!  if(allocated(WTGRBPs1))deallocate(WTGRBPs1)
!  if(allocated(WTSTKBs1))deallocate(WTSTKBs1)
!  if(allocated(WTSTBNs1))deallocate(WTSTBNs1)
!  if(allocated(WTSTBPs1))deallocate(WTSTBPs1)
!  if(allocated(WTSHTBs1))deallocate(WTSHTBs1)
!  if(allocated(WTSHTNs1))deallocate(WTSHTNs1)
!  if(allocated(WTSHTPs1))deallocate(WTSHTPs1)
!  if(allocated(WTSTXBs1))deallocate(WTSTXBs1)
!  if(allocated(WTSTXNs1))deallocate(WTSTXNs1)
!  if(allocated(WTSTXPs1))deallocate(WTSTXPs1)
!  if(allocated(WTSTDIs1))deallocate(WTSTDIs1)
!  if(allocated(WTRVXs1))deallocate(WTRVXs1)
!  if(allocated(WTRVCs1))deallocate(WTRVCs1)
!  if(allocated(WTLSs1))deallocate(WTLSs1)
!  if(allocated(WTSHTs1))deallocate(WTSHTs1)
!  if(allocated(WTSHTAs1))deallocate(WTSHTAs1)
!  if(allocated(WTSTGs1))deallocate(WTSTGs1)
!  if(allocated(WTSTGNs1))deallocate(WTSTGNs1)
!  if(allocated(WTSTGPs1))deallocate(WTSTGPs1)
!  if(allocated(WTNDs1))deallocate(WTNDs1)
!  if(allocated(WTNDNs1))deallocate(WTNDNs1)
!  if(allocated(WTRTNs1))deallocate(WTRTNs1)
!  if(allocated(WTSHNs1))deallocate(WTSHNs1)
!  if(allocated(WTRVNs1))deallocate(WTRVNs1)
!  if(allocated(WTNDPs1))deallocate(WTNDPs1)
!  if(allocated(WTRTPs1))deallocate(WTRTPs1)
!  if(allocated(WTSHPs1))deallocate(WTSHPs1)
!  if(allocated(WTRVPs1))deallocate(WTRVPs1)
  end subroutine plt_biom_destroy


!----------------------------------------------------------------------

  subroutine  plt_soilchem_init(this)

  implicit none

  class(plant_soilchem_type) :: this

  allocate(this%THETPMs1(60,0:JZ1))
  allocate(this%DFGSs1(60,0:JZ1))
  allocate(this%VLPOBs1(0:JZ1))
  allocate(this%VLNO3s1(0:JZ1))
  allocate(this%VLPO4s1(0:JZ1))
  allocate(this%VOLYs1(0:JZ1))
  allocate(this%VOLIs1(0:JZ1))
  allocate(this%VOLWs1(0:JZ1))
  allocate(this%VOLAs1(0:JZ1))
  allocate(this%VLNOBs1(0:JZ1))
  allocate(this%VLNH4s1(0:JZ1))
  allocate(this%VLNHBs1(0:JZ1))
  allocate(this%ZOSGLs1(0:JZ1))
  allocate(this%OQCs1(0:jcplx11,0:JZ1))
  allocate(this%OQNs1(0:jcplx11,0:JZ1))
  allocate(this%OQPs1(0:jcplx11,0:JZ1))
  allocate(this%ZNO3Ss1(0:JZ1))
  allocate(this%ZNO3Bs1(0:JZ1))
  allocate(this%ZNSGLs1(0:JZ1))
  allocate(this%ZNH4Ss1(0:JZ1))
  allocate(this%ZNH4Bs1(0:JZ1))
  allocate(this%Z2OSs1(0:JZ1))
  allocate(this%ZNH3Ss1(0:JZ1))
  allocate(this%ZNH3Bs1(0:JZ1))
  allocate(this%H1PO4s1(0:JZ1))
  allocate(this%H2PO4s1(0:JZ1))
  allocate(this%H1POBs1(0:JZ1))
  allocate(this%H2POBs1(0:JZ1))
  allocate(this%H2GSs1(0:JZ1))
  allocate(this%HLSGLs1(0:JZ1))
  allocate(this%OXYGs1(0:JZ1))
  allocate(this%OXYSs1(0:JZ1))
  allocate(this%OLSGLs1(0:JZ1))
  allocate(this%POSGLs1(0:JZ1))
  allocate(this%CCH4Gs1(0:JZ1))
  allocate(this%CZ2OGs1(0:JZ1))
  allocate(this%CNH3Gs1(0:JZ1))
  allocate(this%CH2GGs1(0:JZ1))
  allocate(this%CORGCs1(0:JZ1))
  allocate(this%Z2SGLs1(JZ1))
  allocate(this%ZHSGLs1(JZ1))
  allocate(this%CH2P4s1(0:JZ1))
  allocate(this%CH1P4Bs1(0:JZ1))
  allocate(this%CH2P4Bs1(0:JZ1))
  allocate(this%CNO3Ss1(0:JZ1))
  allocate(this%CH1P4s1(0:JZ1))
  allocate(this%CNO3Bs1(0:JZ1))
  allocate(this%CNH4Ss1(0:JZ1))
  allocate(this%CNH4Bs1(0:JZ1))
  allocate(this%CO2Gs1(0:JZ1))
  allocate(this%CO2Ss1(0:JZ1))
  allocate(this%CH4Ss1(0:JZ1))
  allocate(this%CCH4Ss1(0:JZ1))
  allocate(this%CZ2OSs1(0:JZ1))
  allocate(this%CNH3Ss1(0:JZ1))
  allocate(this%CNH3Bs1(0:JZ1))
  allocate(this%CH2GSs1(0:JZ1))
  allocate(this%CLSGLs1(0:JZ1))
  allocate(this%CQSGLs1(0:JZ1))
  allocate(this%VOLXs1(0:JZ1))
  allocate(this%THETWs1(0:JZ1))
  allocate(this%THETYs1(0:JZ1))
  allocate(this%SCO2Ls1(0:JZ1))
  allocate(this%SOXYLs1(0:JZ1))
  allocate(this%SCH4Ls1(0:JZ1))
  allocate(this%SN2OLs1(0:JZ1))
  allocate(this%SNH3Ls1(0:JZ1))
  allocate(this%SH2GLs1(0:JZ1))
  allocate(this%CGSGLs1(JZ1))
  allocate(this%CHSGLs1(JZ1))
  allocate(this%HGSGLs1(JZ1))
  allocate(this%OGSGLs1(JZ1))
  allocate(this%RSCSs1(JZ1))
  allocate(this%BKDSs1(0:JZ1))
  allocate(this%CNDUs1(JZ1))
  allocate(this%CPO4Ss1(JZ1))
  allocate(this%PSISTs1(0:JZ1))
  allocate(this%ZVSGLs1(0:JZ1))
  end subroutine plt_soilchem_init
!----------------------------------------------------------------------

  subroutine plt_soilchem_destroy(this)
  implicit none
  class(plant_soilchem_type) :: this

!  if(allocated(THETPMs1))deallocate(THETPMs1)
!  if(allocated(DFGSs1))deallocate(DFGSs1)
!  if(allocated(ZVSGLs1))deallocate(ZVSGLs1)
!  if(allocated(PSISTs1))deallocate(PSISTs1)
!  if(allocated(SOXYLs1))deallocate(SOXYLs1)
!  if(allocated(CPO4Ss1))deallocate(CPO4Ss1)
!  if(allocated(CNDUs1))deallocate(CNDUs1)
!  if(allocated(CGSGLs1))deallocate(CGSGLs1)
!  if(allocated(CHSGLs1))deallocate(CHSGLs1)
!  if(allocated(HGSGLs1))deallocate(HGSGLs1)
!  if(allocated(OGSGLs1))deallocate(OGSGLs1)
!  if(allocated(RSCSs1))deallocate(RSCSs1)
!  if(allocated(SCO2Ls1))deallocate(SCO2Ls1)
!  if(allocated(SCH4Ls1))deallocate(SCH4Ls1)
!  if(allocated(SN2OLs1))deallocate(SN2OLs1)
!  if(allocated(SNH3Ls1))deallocate(SNH3Ls1)
!  if(allocated(SH2GLs1))deallocate(SH2GLs1)
!  if(allocated(THETWs1))deallocate(THETWs1)
!  if(allocated(THETYs1))deallocate(THETYs1)
!  if(allocated(VOLXs1))deallocate(VOLXs1)
!  if(allocated(CCH4Gs1))deallocate(CCH4Gs1)
!  if(allocated(CZ2OGs1))deallocate(CZ2OGs1)
!  if(allocated(CNH3Gs1))deallocate(CNH3Gs1)
!  if(allocated(CH2GGs1))deallocate(CH2GGs1)
!  if(allocated(CORGCs1))deallocate(CORGCs1)
!  if(allocated(H1PO4s1))deallocate(H1PO4s1)
!  if(allocated(H2PO4s1))deallocate(H2PO4s1)
!  if(allocated(H1POBs1))deallocate(H1POBs1)
!  if(allocated(H2POBs1))deallocate(H2POBs1)
!  if(allocated(H2GSs1))deallocate(H2GSs1)
!  if(allocated(HLSGLs1))deallocate(HLSGLs1)
!  if(allocated(OXYGs1))deallocate(OXYGs1)
!  if(allocated(OXYSs1))deallocate(OXYSs1)
!  if(allocated(OLSGLs1))deallocate(OLSGLs1)
!  if(allocated(POSGLs1))deallocate(POSGLs1)
!  if(allocated(VLPOBs1))deallocate(VLPOBs1)
!  if(allocated(VLNO3s1))deallocate(VLNO3s1)
!  if(allocated(VLPO4s1))deallocate(VLPO4s1)
!  if(allocated(VOLYs1))deallocate(VOLYs1)
!  if(allocated(VOLIs1))deallocate(VOLIs1)
!  if(allocated(VOLWs1))deallocate(VOLWs1)
!  if(allocated(VOLAs1))deallocate(VOLAs1)
!  if(allocated(VLNOBs1))deallocate(VLNOBs1)
!  if(allocated(VLNH4s1))deallocate(VLNH4s1)
!  if(allocated(VLNHBs1))deallocate(VLNHBs1)
!  if(allocated(ZOSGLs1))deallocate(ZOSGLs1)

!  if(allocated(OQCs1))deallocate(OQCs1)
!  if(allocated(OQNs1))deallocate(OQNs1)
!  if(allocated(OQPs1))deallocate(OQPs1)
!  if(allocated(ZNO3Ss1))deallocate(ZNO3Ss1)
!  if(allocated(ZNO3Bs1))deallocate(ZNO3Bs1)
!  if(allocated(ZNSGLs1))deallocate(ZNSGLs1)
!  if(allocated(ZNH4Ss1))deallocate(ZNH4Ss1)
!  if(allocated(ZNH4Bs1))deallocate(ZNH4Bs1)
!  if(allocated(Z2OSs1))deallocate(Z2OSs1)
!  if(allocated(ZNH3Ss1))deallocate(ZNH3Ss1)
!  if(allocated(ZNH3Bs1))deallocate(ZNH3Bs1)
!  if(allocated(Z2SGLs1))deallocate(Z2SGLs1)
!  if(allocated(ZHSGLs1))deallocate(ZHSGLs1)
!  if(allocated(BKDSs1))deallocate(BKDSs1)
!  if(allocated(CH2P4s1))deallocate(CH2P4s1)
!  if(allocated(CH1P4Bs1))deallocate(CH1P4Bs1)
!  if(allocated(CH2P4Bs1))deallocate(CH2P4Bs1)
!  if(allocated(CNO3Ss1))deallocate(CNO3Ss1)
!  if(allocated(CH1P4s1))deallocate(CH1P4s1)
!  if(allocated(CNO3Bs1))deallocate(CNO3Bs1)
!  if(allocated(CNH4Ss1))deallocate(CNH4Ss1)
!  if(allocated(CNH4Bs1))deallocate(CNH4Bs1)
!  if(allocated(CO2Gs1))deallocate(CO2Gs1)
!  if(allocated(CO2Ss1))deallocate(CO2Ss1)
!  if(allocated(CH4Ss1))deallocate(CH4Ss1)
!  if(allocated(CCH4Ss1))deallocate(CCH4Ss1)
!  if(allocated(CZ2OSs1))deallocate(CZ2OSs1)
!  if(allocated(CNH3Ss1))deallocate(CNH3Ss1)
!  if(allocated(CNH3Bs1))deallocate(CNH3Bs1)
!  if(allocated(CH2GSs1))deallocate(CH2GSs1)
!  if(allocated(CLSGLs1))deallocate(CLSGLs1)
!  if(allocated(CQSGLs1))deallocate(CQSGLs1)

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
  jcplx11=jcplx1
  JLI1=JLI
  JNODS1=JNODS

  call plt_pheno%Init()

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

  call plt_rad%Destroy()

  allocate(IYR0s1(JP1))
  allocate(INTYPs1(JP1))
  allocate(IYRXs1(JP1))
  allocate(IYRYs1(JP1))
  allocate(IDAYXs1(JP1))
  allocate(IDAY0s1(JP1))
  allocate(IDAYHs1(JP1))
  allocate(IYRHs1(JP1))
  allocate(IHVSTs1(JP1))
  allocate(JHVSTs1(JP1))

  allocate(IDATAs1(60))
  allocate(ICTYPs1(JP1))
  allocate(FERTs1(1:20))
  allocate(DLYR3s1(0:JZ1))


  allocate(CDPTHZs1(0:JZ1))
  allocate(DPTHZs1(0:JZ1))
  allocate(FMPRs1(0:JZ1))

  allocate(RPO4Ys1(0:JZ1))
  allocate(RPOBYs1(0:JZ1))
  allocate(RP14Ys1(0:JZ1))
  allocate(RP1BYs1(0:JZ1))
  allocate(RNO3Ys1(0:JZ1))
  allocate(RNH4Ys1(0:JZ1))
  allocate(RNHBYs1(0:JZ1))
  allocate(RN3BYs1(0:JZ1))
  allocate(ROXYFs1(0:JZ1))
  allocate(RCO2Fs1(0:JZ1))
  allocate(ROXYLs1(0:JZ1))
  allocate(ROXYYs1(0:JZ1))
  allocate(TKSs1(0:JZ1))
  allocate(TFNDs1(0:JZ1))
  allocate(AREA3s1(0:JZ1))

  allocate(RP1BXs1(0:JZ1))
  allocate(RNO3Xs1(0:JZ1))
  allocate(RPO4Xs1(0:JZ1))
  allocate(RNH4Xs1(0:JZ1))
  allocate(ROXYXs1(0:JZ1))
  allocate(RPOBXs1(0:JZ1))
  allocate(RN3BXs1(0:JZ1))
  allocate(RNHBXs1(0:JZ1))
  allocate(RP14Xs1(0:JZ1))
  allocate(TUPN3Bs1(JZ1))
  allocate(TUPH1Bs1(JZ1))
  allocate(TUPN3Ss1(JZ1))
  allocate(TUPHGSs1(JZ1))
  allocate(TUPOXPs1(JZ1))
  allocate(TNHFLAs1(JZ1))
  allocate(TCOFLAs1(JZ1))
  allocate(TN2FLAs1(JZ1))
  allocate(TOXFLAs1(JZ1))
  allocate(TCHFLAs1(JZ1))
  allocate(TLCH4Ps1(JZ1))
  allocate(TLOXYPs1(JZ1))
  allocate(TLNH3Ps1(JZ1))
  allocate(TLCO2Ps1(JZ1))
  allocate(TLN2OPs1(JZ1))
  allocate(TLH2GPs1(JZ1))
  allocate(THGFLAs1(JZ1))
  allocate(TUPOXSs1(JZ1))
  allocate(TCO2Ps1(JZ1))
  allocate(TUPN2Ss1(JZ1))
  allocate(TUPCHSs1(JZ1))
  allocate(TUPNH4s1(JZ1))
  allocate(TCO2Ss1(JZ1))
  allocate(TUPNO3s1(JZ1))
  allocate(TUPH2Ps1(JZ1))
  allocate(TUPH2Bs1(JZ1))
  allocate(TUPH1Ps1(JZ1))
  allocate(RTDNTs1(JZ1))
  allocate(TUPNFs1(JZ1))
  allocate(TUPNHBs1(JZ1))
  allocate(TUPNOBs1(JZ1))
  allocate(TUPHTs1(0:JZ1))
  allocate(TUPWTRs1(0:JZ1))
  allocate(TDFOMCs1(0:jcplx11,JZ1))
  allocate(TDFOMNs1(0:jcplx11,JZ1))
  allocate(TDFOMPs1(0:jcplx11,JZ1))
  allocate(CSNTs1(jcplx11,0:1,0:JZ1))
  allocate(ZSNTs1(jcplx11,0:1,0:JZ1))
  allocate(PSNTs1(jcplx11,0:1,0:JZ1))
  allocate(RDFOMCs1(2,0:jcplx11,0:JZ1,JP1))
  allocate(RDFOMNs1(2,0:jcplx11,0:JZ1,JP1))
  allocate(RDFOMPs1(2,0:jcplx11,0:JZ1,JP1))
  allocate(FOSRHs1(0:jcplx11,0:JZ1))
  allocate(XOQCSs1(0:jcplx11,0:JZ1))
  allocate(XOQNSs1(0:jcplx11,0:JZ1))
  allocate(XOQPSs1(0:jcplx11,0:JZ1))
  allocate(FLWCs1(JP1))
  allocate(EHVSTs1(1:2,1:4,JP1))
  allocate(HVSTs1(JP1))
  allocate(THINs1(JP1))
  allocate(BALCs1(JP1))
  allocate(BALNs1(JP1))
  allocate(BALPs1(JP1))
  allocate(PPIs1(JP1))
  allocate(PPZs1(JP1))
  allocate(PPXs1(JP1))
  allocate(SSTXs1(JP1))
  allocate(ZTYPs1(JP1))
  allocate(ZTYPIs1(JP1))
  allocate(CTCs1(JP1))
  allocate(CNETs1(JP1))
  allocate(CARBNs1(JP1))


  allocate(RSMXs1(JP1))
  allocate(CWSRTs1(JP1))
  allocate(CFIs1(JP1))
  allocate(CFXs1(JP1))
  allocate(DCO2s1(JP1))
  allocate(DATAPs1(JP1))
  allocate(DATAs1(30))

  allocate(FMOLs1(JP1))
  allocate(FNODs1(JP1))
  allocate(GROUPIs1(JP1))
  allocate(HTCTLs1(JP1))
  allocate(HCSNCs1(JP1))
  allocate(HZSNCs1(JP1))
  allocate(HPSNCs1(JP1))
  allocate(HVSTCs1(JP1))
  allocate(HVSTNs1(JP1))
  allocate(HVSTPs1(JP1))
  allocate(HCUPTKs1(JP1))
  allocate(HZUPTKs1(JP1))
  allocate(HPUPTKs1(JP1))
  allocate(OFFSTs1(JP1))
  allocate(OSTRs1(JP1))
  allocate(PPs1(JP1))
  allocate(PSILTs1(JP1))
  allocate(PBs1(JP1))
  allocate(PRs1(JP1))
  allocate(HTCs1(JP1))
  allocate(PSILGs1(JP1))
  allocate(RCMXs1(JP1))
  allocate(RCO2Zs1(JP1))
  allocate(ROXYZs1(JP1))
  allocate(RCH4Zs1(JP1))
  allocate(RN2OZs1(JP1))
  allocate(RNH3Zs1(JP1))
  allocate(RH2GZs1(JP1))

  allocate(RTFQs1(JP1))
  allocate(RAZs1(JP1))
  allocate(PSILZs1(JP1))
  allocate(RSMNs1(JP1))

  allocate(TKCZs1(JP1))
  allocate(SFLXCs1(JP1))
  allocate(SO2s1(JP1))
  allocate(RAs1(JP1))
  allocate(RCs1(JP1))
  allocate(EVAPCs1(JP1))
  allocate(HFLXCs1(JP1))
  allocate(EFLXCs1(JP1))
  allocate(TCZs1(JP1))
  allocate(TCGs1(JP1))
  allocate(TCXs1(JP1))
  allocate(TKGs1(JP1))
  allocate(TFN3s1(JP1))
  allocate(TCSN0s1(JP1))
  allocate(TZSN0s1(JP1))
  allocate(TPSN0s1(JP1))
  allocate(TCSNCs1(JP1))
  allocate(TZSNCs1(JP1))
  allocate(TPSNCs1(JP1))
  allocate(TCO2Ts1(JP1))
  allocate(TCUPTKs1(JP1))
  allocate(THVSTCs1(JP1))
  allocate(TNH3Cs1(JP1))
  allocate(TZUPTKs1(JP1))
  allocate(THVSTNs1(JP1))
  allocate(TZUPFXs1(JP1))
  allocate(TPUPTKs1(JP1))
  allocate(THVSTPs1(JP1))
  allocate(UPOMCs1(JP1))
  allocate(UPNFs1(JP1))
  allocate(UPNO3s1(JP1))
  allocate(UPNH4s1(JP1))
  allocate(UPOMNs1(JP1))
  allocate(UPH1Ps1(JP1))
  allocate(UPH2Ps1(JP1))
  allocate(UPOMPs1(JP1))
  allocate(VOLWCs1(JP1))
  allocate(VHCPCs1(JP1))
  allocate(VCH4Fs1(JP1))
  allocate(VCO2Fs1(JP1))
  allocate(VN2OFs1(JP1))
  allocate(VNH3Fs1(JP1))
  allocate(VPO4Fs1(JP1))
  allocate(VOLWPs1(JP1))
  allocate(WSTRs1(JP1))

  allocate(RNH3Cs1(JP1))
  allocate(VOXYFs1(JP1))
  allocate(WDLFs1(JP1))
  allocate(GFILLs1(JP1))
  allocate(STMXs1(JP1))
  allocate(SDMXs1(JP1))
  allocate(EPs1(JP1))
  allocate(IRTYPs1(JP1))
  allocate(GRMXs1(JP1))
  allocate(PTSHTs1(JP1))

  allocate(TCCs1(JP1))
  allocate(ENGYXs1(JP1))
  allocate(DTKCs1(JP1))
  allocate(TKCs1(JP1))
  allocate(OSMOs1(JP1))
  allocate(TCO2As1(JP1))
  allocate(PSILOs1(JP1))
  allocate(CTRANs1(JP1))

  allocate(GRDMs1(JP1))
  allocate(ZEROLs1(JP1))
  allocate(ZEROPs1(JP1))

  allocate(ZEROQs1(JP1))
  allocate(ZNPPs1(JP1))

  allocate(GRWTBs1(JC1,JP1))

  allocate(IDAYYs1(JP1))

  allocate(RCCLXs1(JC1,JP1))
  allocate(RCZLXs1(JC1,JP1))
  allocate(RCPLXs1(JC1,JP1))
  allocate(RCCSXs1(JC1,JP1))
  allocate(RCZSXs1(JC1,JP1))
  allocate(RCPSXs1(JC1,JP1))
  allocate(RNH3Bs1(JC1,JP1))

  allocate(SURFXs1(JLI1,JC1,JNODS1,JC1,JP1))
  allocate(CPOOL3s1(JNODS1,JC1,JP1))
  allocate(CPOOL4s1(JNODS1,JC1,JP1))
  allocate(CO2Bs1(JNODS1,JC1,JP1))
  allocate(COMPLs1(JNODS1,JC1,JP1))
  allocate(CBXNs1(JNODS1,JC1,JP1))
  allocate(CBXN4s1(JNODS1,JC1,JP1))
  allocate(ETGROs1(JNODS1,JC1,JP1))
  allocate(ETGR4s1(JNODS1,JC1,JP1))
  allocate(FDBK4s1(JNODS1,JC1,JP1))
  allocate(HCOBs1(JNODS1,JC1,JP1))
  allocate(VCGROs1(JNODS1,JC1,JP1))
  allocate(VGROs1(JNODS1,JC1,JP1))
  allocate(VCGR4s1(JNODS1,JC1,JP1))
  allocate(VGRO4s1(JNODS1,JC1,JP1))


  allocate(CO2Ps1(2,JZ1,JP1))
  allocate(CO2As1(2,JZ1,JP1))
  allocate(CH4Ps1(2,JZ1,JP1))
  allocate(CH4As1(2,JZ1,JP1))
  allocate(H2GPs1(2,JZ1,JP1))
  allocate(H2GAs1(2,JZ1,JP1))
  allocate(OXYPs1(2,JZ1,JP1))
  allocate(OXYAs1(2,JZ1,JP1))
  allocate(PSIRTs1(2,JZ1,JP1))
  allocate(PSIROs1(2,JZ1,JP1))
  allocate(PSIRGs1(2,JZ1,JP1))
  allocate(RUPNHBs1(2,JZ1,JP1))
  allocate(RUPNH4s1(2,JZ1,JP1))
  allocate(RUPH2Ps1(2,JZ1,JP1))
  allocate(RUPNOBs1(2,JZ1,JP1))
  allocate(RUPNO3s1(2,JZ1,JP1))
  allocate(RUPH1Bs1(2,JZ1,JP1))
  allocate(RUPH1Ps1(2,JZ1,JP1))
  allocate(RUPH2Bs1(2,JZ1,JP1))
  allocate(RUONHBs1(2,JZ1,JP1))
  allocate(RUONH4s1(2,JZ1,JP1))
  allocate(RUOH2Ps1(2,JZ1,JP1))
  allocate(RUONOBs1(2,JZ1,JP1))
  allocate(RUONO3s1(2,JZ1,JP1))
  allocate(RUOH1Bs1(2,JZ1,JP1))
  allocate(RUOH1Ps1(2,JZ1,JP1))
  allocate(RUOH2Bs1(2,JZ1,JP1))
  allocate(RUCNHBs1(2,JZ1,JP1))
  allocate(RUCNH4s1(2,JZ1,JP1))
  allocate(RUCH2Ps1(2,JZ1,JP1))
  allocate(RUCNOBs1(2,JZ1,JP1))
  allocate(RUCNO3s1(2,JZ1,JP1))
  allocate(RUCH1Bs1(2,JZ1,JP1))
  allocate(RUCH1Ps1(2,JZ1,JP1))
  allocate(RUCH2Bs1(2,JZ1,JP1))
  allocate(RCO2Ms1(2,JZ1,JP1))
  allocate(RCO2Ns1(2,JZ1,JP1))
  allocate(RCO2As1(2,JZ1,JP1))
  allocate(RTN1s1(2,JZ1,JP1))
  allocate(RTNLs1(2,JZ1,JP1))
  allocate(RTLGPs1(2,JZ1,JP1))
  allocate(RTDNPs1(2,JZ1,JP1))
  allocate(RTVLPs1(2,JZ1,JP1))
  allocate(RTVLWs1(2,JZ1,JP1))
  allocate(RRAD1s1(2,JZ1,JP1))
  allocate(RRAD2s1(2,JZ1,JP1))
  allocate(RTARPs1(2,JZ1,JP1))
  allocate(RTLGAs1(2,JZ1,JP1))
  allocate(RCO2Ps1(2,JZ1,JP1))
  allocate(RUPOXPs1(2,JZ1,JP1))
  allocate(RCO2Ss1(2,JZ1,JP1))
  allocate(RUPOXSs1(2,JZ1,JP1))
  allocate(RUPCHSs1(2,JZ1,JP1))
  allocate(RUPN2Ss1(2,JZ1,JP1))
  allocate(RUPN3Ss1(2,JZ1,JP1))
  allocate(RUPN3Bs1(2,JZ1,JP1))
  allocate(RUPHGSs1(2,JZ1,JP1))
  allocate(RCOFLAs1(2,JZ1,JP1))
  allocate(ROXFLAs1(2,JZ1,JP1))
  allocate(RCHFLAs1(2,JZ1,JP1))
  allocate(RN2FLAs1(2,JZ1,JP1))
  allocate(RNHFLAs1(2,JZ1,JP1))
  allocate(RHGFLAs1(2,JZ1,JP1))
  allocate(RCODFAs1(2,JZ1,JP1))
  allocate(ROXDFAs1(2,JZ1,JP1))
  allocate(RCHDFAs1(2,JZ1,JP1))
  allocate(RN2DFAs1(2,JZ1,JP1))
  allocate(RNHDFAs1(2,JZ1,JP1))
  allocate(RHGDFAs1(2,JZ1,JP1))
  allocate(ROXYPs1(2,JZ1,JP1))
  allocate(RUNNHPs1(2,JZ1,JP1))
  allocate(RUNNBPs1(2,JZ1,JP1))
  allocate(RUNNOPs1(2,JZ1,JP1))
  allocate(RUNNXPs1(2,JZ1,JP1))
  allocate(RUPP2Ps1(2,JZ1,JP1))
  allocate(RUPP2Bs1(2,JZ1,JP1))
  allocate(RUPP1Ps1(2,JZ1,JP1))
  allocate(RUPP1Bs1(2,JZ1,JP1))
  allocate(UPWTRs1(2,JZ1,JP1))
  allocate(WFRs1(2,JZ1,JP1))

  allocate(Z2OPs1(2,JZ1,JP1))
  allocate(Z2OAs1(2,JZ1,JP1))
  allocate(ZH3Ps1(2,JZ1,JP1))
  allocate(ZH3As1(2,JZ1,JP1))
  allocate(CPOOLNs1(JZ1,JP1))
  allocate(RUPNFs1(JZ1,JP1))
  allocate(PPOOLNs1(JZ1,JP1))
  allocate(TFN4s1(JZ1,JP1))
  allocate(WTNDLNs1(JZ1,JP1))
  allocate(WTNDLs1(JZ1,JP1))
  allocate(WTNDLPs1(JZ1,JP1))
  allocate(ZPOOLNs1(JZ1,JP1))
  allocate(WGLFVs1(JC1,JP1))
  allocate(DMVLs1(2,JP1))
  allocate(PORTs1(2,JP1))
  allocate(PORTXs1(2,JP1))
  allocate(RTAR2Xs1(2,JP1))
  allocate(RTAR1Xs1(2,JP1))
  allocate(RRAD1Xs1(2,JP1))
  allocate(RRAD2Xs1(2,JP1))
  allocate(RTLG1Xs1(2,JP1))
  allocate(RRADPs1(2,JP1))
  allocate(RTLG2Xs1(2,JP1))
  allocate(RSRRs1(2,JP1))
  allocate(RSRAs1(2,JP1))
  allocate(RRAD1Ms1(2,JP1))
  allocate(RRAD2Ms1(2,JP1))
  allocate(UPMNPOs1(2,JP1))
  allocate(UPMXPOs1(2,JP1))
  allocate(UPKMPOs1(2,JP1))
  allocate(UPMNZOs1(2,JP1))
  allocate(UPMXZOs1(2,JP1))
  allocate(UPKMZOs1(2,JP1))
  allocate(UPMNZHs1(2,JP1))
  allocate(UPMXZHs1(2,JP1))
  allocate(UPKMZHs1(2,JP1))
  allocate(RTLG1s1(2,JZ1,JC1,JP1))
  allocate(RTLG2s1(2,JZ1,JC1,JP1))
  allocate(RTN2s1(2,JZ1,JC1,JP1))
  allocate(WTRT2s1(2,JZ1,JC1,JP1))
  allocate(WTRT1s1(2,JZ1,JC1,JP1))
  allocate(WTRT2Ns1(2,JZ1,JC1,JP1))
  allocate(WTRT1Ns1(2,JZ1,JC1,JP1))
  allocate(WTRT2Ps1(2,JZ1,JC1,JP1))
  allocate(WTRT1Ps1(2,JZ1,JC1,JP1))
  allocate(RTDP1s1(2,JC1,JP1))
  allocate(RTWT1s1(2,JC1,JP1))
  allocate(RTWT1Ns1(2,JC1,JP1))
  allocate(RTWT1Ps1(2,JC1,JP1))
  allocate(WTSTDGs1(4,JP1))
  allocate(WTSTDNs1(4,JP1))
  allocate(WTSTDPs1(4,JP1))
  allocate(CSNCs1(4,0:1,0:JZ1,JP1))
  allocate(PSNCs1(4,0:1,0:JZ1,JP1))
  allocate(ZSNCs1(4,0:1,0:JZ1,JP1))
  allocate(CFOPCs1(0:5,4,JP1))
  allocate(CFOPNs1(0:5,4,JP1))
  allocate(CFOPPs1(0:5,4,JP1))
  allocate(PARs1(JLI1,JSA1,JC1,JP1))
  allocate(PARDIFs1(JLI1,JSA1,JC1,JP1))
  allocate(VOLWMs1(60,0:JZ1))
  allocate(VOLPMs1(60,0:JZ1))
  allocate(TORTs1(60,0:JZ1))
  allocate(FILMs1(60,0:JZ1))
  allocate(ROXSKs1(60,0:JZ1))

  end subroutine InitAllocate
!----------------------------------------------------------------------
  subroutine DestructPlantAPIData
  implicit none

  if(allocated(IYR0s1))deallocate(IYR0s1)
  if(allocated(INTYPs1))deallocate(INTYPs1)
  if(allocated(IDAYXs1))deallocate(IDAYXs1)
  if(allocated(IYRYs1))deallocate(IYRYs1)
  if(allocated(IYRXs1))deallocate(IYRXs1)

  if(allocated(IDAY0s1))deallocate(IDAY0s1)
  if(allocated(IDAYHs1))deallocate(IDAYHs1)
  if(allocated(IYRHs1))deallocate(IYRHs1)
  if(allocated(IHVSTs1))deallocate(IHVSTs1)

  if(allocated(JHVSTs1))deallocate(JHVSTs1)

  if(allocated(IDATAs1))deallocate(IDATAs1)

  if(allocated(ICTYPs1))deallocate(ICTYPs1)
  if(allocated(DATAs1))deallocate(DATAs1)

  if(allocated(FERTs1))deallocate(FERTs1)

  if(allocated(ZTYPs1))deallocate(ZTYPs1)
  if(allocated(ZTYPIs1))deallocate(ZTYPIs1)
  if(allocated(GRDMs1))deallocate(GRDMs1)
  if(allocated(HTCs1))deallocate(HTCs1)
  if(allocated(SSTXs1)) deallocate(SSTXs1)

  if(allocated(CDPTHZs1))deallocate(CDPTHZs1)

  if(allocated(CFIs1)) deallocate(CFIs1)
  if(allocated(PPIs1)) deallocate(PPIs1)
  if(allocated(PPZs1)) deallocate(PPZs1)

  if(allocated(PPXs1)) deallocate(PPXs1)
  if(allocated(DPTHZs1))deallocate(DPTHZs1)
  if(allocated(FMPRs1))deallocate(FMPRs1)

  if(allocated(RPO4Ys1))deallocate(RPO4Ys1)
  if(allocated(RPOBYs1))deallocate(RPOBYs1)
  if(allocated(RP14Ys1))deallocate(RP14Ys1)
  if(allocated(RP1BYs1))deallocate(RP1BYs1)
  if(allocated(RNO3Ys1))deallocate(RNO3Ys1)
  if(allocated(RNH4Ys1))deallocate(RNH4Ys1)
  if(allocated(RNHBYs1))deallocate(RNHBYs1)
  if(allocated(RN3BYs1))deallocate(RN3BYs1)
  if(allocated(ROXYFs1))deallocate(ROXYFs1)
  if(allocated(RCO2Fs1))deallocate(RCO2Fs1)
  if(allocated(ROXYLs1))deallocate(ROXYLs1)
  if(allocated(ROXYYs1))deallocate(ROXYYs1)

  if(allocated(DLYR3s1))deallocate(DLYR3s1)
  if(allocated(RSMNs1))deallocate(RSMNs1)


  if(allocated(TKCZs1))deallocate(TKCZs1)
  if(allocated(SO2s1))deallocate(SO2s1)
  if(allocated(RCs1))deallocate(RCs1)
  if(allocated(RAs1))deallocate(RAs1)
  if(allocated(SFLXCs1))deallocate(SFLXCs1)
  if(allocated(PSILZs1))deallocate(PSILZs1)

  if(allocated(TKSs1))deallocate(TKSs1)
  if(allocated(TFNDs1))deallocate(TFNDs1)

  if(allocated(AREA3s1))deallocate(AREA3s1)


  if(allocated(RP1BXs1))deallocate(RP1BXs1)
  if(allocated(RNO3Xs1))deallocate(RNO3Xs1)
  if(allocated(RP14Xs1))deallocate(RP14Xs1)
  if(allocated(RNH4Xs1))deallocate(RNH4Xs1)
  if(allocated(ROXYXs1))deallocate(ROXYXs1)
  if(allocated(RPO4Xs1))deallocate(RPO4Xs1)
  if(allocated(RNHBXs1))deallocate(RNHBXs1)
  if(allocated(RN3BXs1))deallocate(RN3BXs1)
  if(allocated(RPOBXs1))deallocate(RPOBXs1)
  if(allocated(TCO2Ss1))deallocate(TCO2Ss1)
  if(allocated(TUPNH4s1))deallocate(TUPNH4s1)
  if(allocated(TUPN2Ss1))deallocate(TUPN2Ss1)
  if(allocated(TUPOXPs1))deallocate(TUPOXPs1)
  if(allocated(TUPOXSs1))deallocate(TUPOXSs1)
  if(allocated(TCO2Ps1))deallocate(TCO2Ps1)
  if(allocated(THGFLAs1))deallocate(THGFLAs1)
  if(allocated(TN2FLAs1))deallocate(TN2FLAs1)
  if(allocated(TNHFLAs1))deallocate(TNHFLAs1)
  if(allocated(TCOFLAs1))deallocate(TCOFLAs1)
  if(allocated(TOXFLAs1))deallocate(TOXFLAs1)
  if(allocated(TCHFLAs1))deallocate(TCHFLAs1)
  if(allocated(TLCH4Ps1))deallocate(TLCH4Ps1)
  if(allocated(TLCO2Ps1))deallocate(TLCO2Ps1)
  if(allocated(TLNH3Ps1))deallocate(TLNH3Ps1)
  if(allocated(TUPHTs1))deallocate(TUPHTs1)
  if(allocated(RTDNTs1))deallocate(RTDNTs1)
  if(allocated(TUPNFs1))deallocate(TUPNFs1)
  if(allocated(TDFOMPs1))deallocate(TDFOMPs1)
  if(allocated(TDFOMNs1))deallocate(TDFOMNs1)
  if(allocated(TDFOMCs1))deallocate(TDFOMCs1)
  if(allocated(TUPWTRs1))deallocate(TUPWTRs1)
  if(allocated(TLOXYPs1))deallocate(TLOXYPs1)
  if(allocated(TLN2OPs1))deallocate(TLN2OPs1)
  if(allocated(TUPHGSs1))deallocate(TUPHGSs1)
  if(allocated(TLH2GPs1))deallocate(TLH2GPs1)
  if(allocated(TUPCHSs1))deallocate(TUPCHSs1)
  if(allocated(TUPN3Ss1))deallocate(TUPN3Ss1)
  if(allocated(TUPN3Bs1))deallocate(TUPN3Bs1)
  if(allocated(TUPH1Bs1))deallocate(TUPH1Bs1)
  if(allocated(TUPNO3s1))deallocate(TUPNO3s1)
  if(allocated(TUPH1Ps1))deallocate(TUPH1Ps1)
  if(allocated(TUPNOBs1))deallocate(TUPNOBs1)
  if(allocated(TUPNHBs1))deallocate(TUPNHBs1)
  if(allocated(TUPH2Ps1))deallocate(TUPH2Ps1)
  if(allocated(TUPH2Bs1))deallocate(TUPH2Bs1)
  if(allocated(RDFOMCs1))deallocate(RDFOMCs1)
  if(allocated(CSNTs1))deallocate(CSNTs1)
  if(allocated(ZSNTs1))deallocate(ZSNTs1)
  if(allocated(PSNTs1))deallocate(PSNTs1)
  if(allocated(RDFOMNs1))deallocate(RDFOMNs1)
  if(allocated(RDFOMPs1))deallocate(RDFOMPs1)
  if(allocated(FOSRHs1))deallocate(FOSRHs1)


  if(allocated(XOQCSs1))deallocate(XOQCSs1)
  if(allocated(XOQNSs1))deallocate(XOQNSs1)
  if(allocated(XOQPSs1))deallocate(XOQPSs1)
  if(allocated(FLWCs1))deallocate(FLWCs1)
  if(allocated(EHVSTs1))deallocate(EHVSTs1)
  if(allocated(HVSTs1))deallocate(HVSTs1)
  if(allocated(THINs1))deallocate(THINs1)
  if(allocated(BALCs1))deallocate(BALCs1)
  if(allocated(BALNs1))deallocate(BALNs1)
  if(allocated(BALPs1))deallocate(BALPs1)


  if(allocated(CTCs1))deallocate(CTCs1)
  if(allocated(CNETs1))deallocate(CNETs1)

  if(allocated(CARBNs1))deallocate(CARBNs1)

  if(allocated(RSMXs1))deallocate(RSMXs1)

  if(allocated(CFXs1))deallocate(CFXs1)
  if(allocated(DCO2s1))deallocate(DCO2s1)
  if(allocated(DATAPs1))deallocate(DATAPs1)

  if(allocated(FMOLs1))deallocate(FMOLs1)
  if(allocated(FNODs1))deallocate(FNODs1)
  if(allocated(GROUPIs1))deallocate(GROUPIs1)
  if(allocated(HTCTLs1))deallocate(HTCTLs1)
  if(allocated(HCSNCs1))deallocate(HCSNCs1)
  if(allocated(HZSNCs1))deallocate(HZSNCs1)
  if(allocated(HPSNCs1))deallocate(HPSNCs1)
  if(allocated(HVSTCs1))deallocate(HVSTCs1)
  if(allocated(HVSTNs1))deallocate(HVSTNs1)
  if(allocated(HVSTPs1))deallocate(HVSTPs1)
  if(allocated(HCUPTKs1))deallocate(HCUPTKs1)
  if(allocated(HZUPTKs1))deallocate(HZUPTKs1)
  if(allocated(HPUPTKs1))deallocate(HPUPTKs1)

  if(allocated(OFFSTs1))deallocate(OFFSTs1)
  if(allocated(OSTRs1))deallocate(OSTRs1)
  if(allocated(PPs1))deallocate(PPs1)
  if(allocated(PSILTs1))deallocate(PSILTs1)
  if(allocated(PBs1))deallocate(PBs1)
  if(allocated(PRs1))deallocate(PRs1)


  if(allocated(PSILGs1))deallocate(PSILGs1)
  if(allocated(RCMXs1))deallocate(RCMXs1)
  if(allocated(RCO2Zs1))deallocate(RCO2Zs1)
  if(allocated(ROXYZs1))deallocate(ROXYZs1)
  if(allocated(RCH4Zs1))deallocate(RCH4Zs1)
  if(allocated(RN2OZs1))deallocate(RN2OZs1)
  if(allocated(RNH3Zs1))deallocate(RNH3Zs1)
  if(allocated(RH2GZs1))deallocate(RH2GZs1)

  if(allocated(RTFQs1))deallocate(RTFQs1)
  if(allocated(RAZs1))deallocate(RAZs1)

  if(allocated(EVAPCs1))deallocate(EVAPCs1)
  if(allocated(HFLXCs1))deallocate(HFLXCs1)
  if(allocated(EFLXCs1))deallocate(EFLXCs1)

  if(allocated(TCZs1))deallocate(TCZs1)
  if(allocated(TCGs1))deallocate(TCGs1)
  if(allocated(TCXs1))deallocate(TCXs1)
  if(allocated(TKGs1))deallocate(TKGs1)
  if(allocated(TFN3s1))deallocate(TFN3s1)
  if(allocated(TCSN0s1))deallocate(TCSN0s1)
  if(allocated(TZSN0s1))deallocate(TZSN0s1)
  if(allocated(TPSN0s1))deallocate(TPSN0s1)
  if(allocated(TCSNCs1))deallocate(TCSNCs1)
  if(allocated(TZSNCs1))deallocate(TZSNCs1)
  if(allocated(TPSNCs1))deallocate(TPSNCs1)
  if(allocated(TCO2Ts1))deallocate(TCO2Ts1)
  if(allocated(TCUPTKs1))deallocate(TCUPTKs1)
  if(allocated(THVSTCs1))deallocate(THVSTCs1)
  if(allocated(TNH3Cs1))deallocate(TNH3Cs1)
  if(allocated(TZUPTKs1))deallocate(TZUPTKs1)
  if(allocated(THVSTNs1))deallocate(THVSTNs1)
  if(allocated(TZUPFXs1))deallocate(TZUPFXs1)
  if(allocated(TPUPTKs1))deallocate(TPUPTKs1)
  if(allocated(THVSTPs1))deallocate(THVSTPs1)
  if(allocated(UPOMCs1))deallocate(UPOMCs1)
  if(allocated(UPNFs1))deallocate(UPNFs1)
  if(allocated(UPNO3s1))deallocate(UPNO3s1)
  if(allocated(UPNH4s1))deallocate(UPNH4s1)
  if(allocated(UPOMNs1))deallocate(UPOMNs1)
  if(allocated(UPH1Ps1))deallocate(UPH1Ps1)
  if(allocated(UPH2Ps1))deallocate(UPH2Ps1)
  if(allocated(UPOMPs1))deallocate(UPOMPs1)
  if(allocated(VOLWCs1))deallocate(VOLWCs1)
  if(allocated(VHCPCs1))deallocate(VHCPCs1)
  if(allocated(VCH4Fs1))deallocate(VCH4Fs1)
  if(allocated(VCO2Fs1))deallocate(VCO2Fs1)
  if(allocated(VN2OFs1))deallocate(VN2OFs1)
  if(allocated(VNH3Fs1))deallocate(VNH3Fs1)
  if(allocated(VPO4Fs1))deallocate(VPO4Fs1)
  if(allocated(VOLWPs1))deallocate(VOLWPs1)
  if(allocated(WSTRs1))deallocate(WSTRs1)

  if(allocated(ZEROLs1))deallocate(ZEROLs1)
  if(allocated(ZEROPs1))deallocate(ZEROPs1)
  if(allocated(ZEROQs1))deallocate(ZEROQs1)

  if(allocated(ZNPPs1))deallocate(ZNPPs1)

  if(allocated(GRWTBs1))deallocate(GRWTBs1)

  if(allocated(IDAYYs1))deallocate(IDAYYs1)

  if(allocated(RCCLXs1))deallocate(RCCLXs1)
  if(allocated(RCZLXs1))deallocate(RCZLXs1)
  if(allocated(RCPLXs1))deallocate(RCPLXs1)
  if(allocated(RCCSXs1))deallocate(RCCSXs1)
  if(allocated(RCZSXs1))deallocate(RCZSXs1)
  if(allocated(RCPSXs1))deallocate(RCPSXs1)
  if(allocated(RNH3Bs1))deallocate(RNH3Bs1)

  if(allocated(SURFXs1))deallocate(SURFXs1)
  if(allocated(CPOOL3s1))deallocate(CPOOL3s1)
  if(allocated(CPOOL4s1))deallocate(CPOOL4s1)
  if(allocated(CO2Bs1))deallocate(CO2Bs1)
  if(allocated(COMPLs1))deallocate(COMPLs1)
  if(allocated(CBXNs1))deallocate(CBXNs1)
  if(allocated(CBXN4s1))deallocate(CBXN4s1)
  if(allocated(ETGROs1))deallocate(ETGROs1)
  if(allocated(ETGR4s1))deallocate(ETGR4s1)
  if(allocated(FDBK4s1))deallocate(FDBK4s1)
  if(allocated(HCOBs1))deallocate(HCOBs1)
  if(allocated(VCGROs1))deallocate(VCGROs1)
  if(allocated(VGROs1))deallocate(VGROs1)
  if(allocated(VCGR4s1))deallocate(VCGR4s1)
  if(allocated(VGRO4s1))deallocate(VGRO4s1)

  if(allocated(OSMOs1))deallocate(OSMOs1)
  if(allocated(TCCs1))deallocate(TCCs1)
  if(allocated(DTKCs1))deallocate(DTKCs1)
  if(allocated(TKCs1))deallocate(TKCs1)
  if(allocated(WTRT1s1))deallocate(WTRT1s1)

  if(allocated(ENGYXs1))deallocate(ENGYXs1)
  if(allocated(PSILOs1))deallocate(PSILOs1)

  if(allocated(CTRANs1))deallocate(CTRANs1)

  if(allocated(TCO2As1))deallocate(TCO2As1)

  if(allocated(CO2Ps1))deallocate(CO2Ps1)
  if(allocated(CO2As1))deallocate(CO2As1)
  if(allocated(CH4Ps1))deallocate(CH4Ps1)
  if(allocated(CH4As1))deallocate(CH4As1)
  if(allocated(H2GPs1))deallocate(H2GPs1)
  if(allocated(H2GAs1))deallocate(H2GAs1)
  if(allocated(OXYPs1))deallocate(OXYPs1)
  if(allocated(OXYAs1))deallocate(OXYAs1)
  if(allocated(PSIRTs1))deallocate(PSIRTs1)
  if(allocated(PSIROs1))deallocate(PSIROs1)
  if(allocated(PSIRGs1))deallocate(PSIRGs1)

  if(allocated(RUPNHBs1))deallocate(RUPNHBs1)
  if(allocated(RUPNH4s1))deallocate(RUPNH4s1)
  if(allocated(RUPH2Ps1))deallocate(RUPH2Ps1)
  if(allocated(RUPNOBs1))deallocate(RUPNOBs1)
  if(allocated(RUPNO3s1))deallocate(RUPNO3s1)
  if(allocated(RUPH1Bs1))deallocate(RUPH1Bs1)
  if(allocated(RUPH1Ps1))deallocate(RUPH1Ps1)
  if(allocated(RUPH2Bs1))deallocate(RUPH2Bs1)
  if(allocated(RUONHBs1))deallocate(RUONHBs1)
  if(allocated(RUONH4s1))deallocate(RUONH4s1)
  if(allocated(RUOH2Ps1))deallocate(RUOH2Ps1)
  if(allocated(RUONOBs1))deallocate(RUONOBs1)
  if(allocated(RUONO3s1))deallocate(RUONO3s1)
  if(allocated(RUOH1Bs1))deallocate(RUOH1Bs1)
  if(allocated(RUOH1Ps1))deallocate(RUOH1Ps1)
  if(allocated(RUOH2Bs1))deallocate(RUOH2Bs1)
  if(allocated(RUCNHBs1))deallocate(RUCNHBs1)
  if(allocated(RUCNH4s1))deallocate(RUCNH4s1)
  if(allocated(RUCH2Ps1))deallocate(RUCH2Ps1)
  if(allocated(RUCNOBs1))deallocate(RUCNOBs1)
  if(allocated(RUCNO3s1))deallocate(RUCNO3s1)
  if(allocated(RUCH1Bs1))deallocate(RUCH1Bs1)
  if(allocated(RUCH1Ps1))deallocate(RUCH1Ps1)
  if(allocated(RUCH2Bs1))deallocate(RUCH2Bs1)
  if(allocated(RCO2Ms1))deallocate(RCO2Ms1)
  if(allocated(RCO2Ns1))deallocate(RCO2Ns1)
  if(allocated(RCO2As1))deallocate(RCO2As1)
  if(allocated(RTN1s1))deallocate(RTN1s1)
  if(allocated(RTNLs1))deallocate(RTNLs1)
  if(allocated(RTLGPs1))deallocate(RTLGPs1)
  if(allocated(RTDNPs1))deallocate(RTDNPs1)
  if(allocated(RTVLPs1))deallocate(RTVLPs1)
  if(allocated(RTVLWs1))deallocate(RTVLWs1)
  if(allocated(RRAD1s1))deallocate(RRAD1s1)
  if(allocated(RRAD2s1))deallocate(RRAD2s1)
  if(allocated(RTARPs1))deallocate(RTARPs1)
  if(allocated(RTLGAs1))deallocate(RTLGAs1)
  if(allocated(RCO2Ps1))deallocate(RCO2Ps1)
  if(allocated(RUPOXPs1))deallocate(RUPOXPs1)
  if(allocated(RCO2Ss1))deallocate(RCO2Ss1)
  if(allocated(RUPOXSs1))deallocate(RUPOXSs1)
  if(allocated(RUPCHSs1))deallocate(RUPCHSs1)
  if(allocated(RUPN2Ss1))deallocate(RUPN2Ss1)
  if(allocated(RUPN3Ss1))deallocate(RUPN3Ss1)
  if(allocated(RUPN3Bs1))deallocate(RUPN3Bs1)
  if(allocated(RUPHGSs1))deallocate(RUPHGSs1)
  if(allocated(RCOFLAs1))deallocate(RCOFLAs1)
  if(allocated(ROXFLAs1))deallocate(ROXFLAs1)
  if(allocated(RCHFLAs1))deallocate(RCHFLAs1)
  if(allocated(RN2FLAs1))deallocate(RN2FLAs1)
  if(allocated(RNHFLAs1))deallocate(RNHFLAs1)
  if(allocated(RHGFLAs1))deallocate(RHGFLAs1)
  if(allocated(RCODFAs1))deallocate(RCODFAs1)
  if(allocated(ROXDFAs1))deallocate(ROXDFAs1)
  if(allocated(RCHDFAs1))deallocate(RCHDFAs1)
  if(allocated(RN2DFAs1))deallocate(RN2DFAs1)
  if(allocated(RNHDFAs1))deallocate(RNHDFAs1)
  if(allocated(RHGDFAs1))deallocate(RHGDFAs1)
  if(allocated(ROXYPs1))deallocate(ROXYPs1)
  if(allocated(RUNNHPs1))deallocate(RUNNHPs1)
  if(allocated(RUNNBPs1))deallocate(RUNNBPs1)
  if(allocated(RUNNOPs1))deallocate(RUNNOPs1)
  if(allocated(RUNNXPs1))deallocate(RUNNXPs1)
  if(allocated(RUPP2Ps1))deallocate(RUPP2Ps1)
  if(allocated(RUPP2Bs1))deallocate(RUPP2Bs1)
  if(allocated(RUPP1Ps1))deallocate(RUPP1Ps1)
  if(allocated(RUPP1Bs1))deallocate(RUPP1Bs1)
  if(allocated(UPWTRs1))deallocate(UPWTRs1)
  if(allocated(WFRs1))deallocate(WFRs1)

  if(allocated(Z2OPs1))deallocate(Z2OPs1)
  if(allocated(Z2OAs1))deallocate(Z2OAs1)
  if(allocated(ZH3Ps1))deallocate(ZH3Ps1)
  if(allocated(ZH3As1))deallocate(ZH3As1)
  if(allocated(CPOOLNs1))deallocate(CPOOLNs1)
  if(allocated(RUPNFs1))deallocate(RUPNFs1)
  if(allocated(PPOOLNs1))deallocate(PPOOLNs1)
  if(allocated(TFN4s1))deallocate(TFN4s1)
  if(allocated(WTNDLNs1))deallocate(WTNDLNs1)
  if(allocated(WTNDLs1))deallocate(WTNDLs1)
  if(allocated(WTNDLPs1))deallocate(WTNDLPs1)
  if(allocated(ZPOOLNs1))deallocate(ZPOOLNs1)

  if(allocated(WGLFVs1))deallocate(WGLFVs1)
  if(allocated(DMVLs1))deallocate(DMVLs1)
  if(allocated(PORTs1))deallocate(PORTs1)
  if(allocated(PORTXs1))deallocate(PORTXs1)
  if(allocated(RTAR2Xs1))deallocate(RTAR2Xs1)
  if(allocated(RTAR1Xs1))deallocate(RTAR1Xs1)
  if(allocated(RRAD1Xs1))deallocate(RRAD1Xs1)
  if(allocated(RRAD2Xs1))deallocate(RRAD2Xs1)

  if(allocated(VOXYFs1))deallocate(VOXYFs1)

  if(allocated(EPs1))deallocate(EPs1)
  if(allocated(WDLFs1))deallocate(WDLFs1)
  if(allocated(GFILLs1))deallocate(GFILLs1)
  if(allocated(IRTYPs1))deallocate(IRTYPs1)
  if(allocated(STMXs1))deallocate(STMXs1)
  if(allocated(SDMXs1))deallocate(SDMXs1)
  if(allocated(GRMXs1))deallocate(GRMXs1)

  if(allocated(PTSHTs1))deallocate(PTSHTs1)

  if(allocated(RNH3Cs1))deallocate(RNH3Cs1)

  if(allocated(RTLG1Xs1))deallocate(RTLG1Xs1)
  if(allocated(RRADPs1))deallocate(RRADPs1)
  if(allocated(RTLG2Xs1))deallocate(RTLG2Xs1)
  if(allocated(RSRRs1))deallocate(RSRRs1)
  if(allocated(RSRAs1))deallocate(RSRAs1)
  if(allocated(RRAD1Ms1))deallocate(RRAD1Ms1)
  if(allocated(RRAD2Ms1))deallocate(RRAD2Ms1)
  if(allocated(UPMNPOs1))deallocate(UPMNPOs1)
  if(allocated(UPMXPOs1))deallocate(UPMXPOs1)
  if(allocated(UPKMPOs1))deallocate(UPKMPOs1)
  if(allocated(UPMNZOs1))deallocate(UPMNZOs1)
  if(allocated(UPMXZOs1))deallocate(UPMXZOs1)
  if(allocated(UPKMZOs1))deallocate(UPKMZOs1)
  if(allocated(UPMNZHs1))deallocate(UPMNZHs1)
  if(allocated(UPMXZHs1))deallocate(UPMXZHs1)
  if(allocated(UPKMZHs1))deallocate(UPKMZHs1)
  if(allocated(RTLG1s1))deallocate(RTLG1s1)
  if(allocated(RTLG2s1))deallocate(RTLG2s1)
  if(allocated(RTN2s1))deallocate(RTN2s1)
  if(allocated(WTRT2s1))deallocate(WTRT2s1)
  if(allocated(WTRT2Ns1))deallocate(WTRT2Ns1)
  if(allocated(WTRT1Ns1))deallocate(WTRT1Ns1)
  if(allocated(WTRT2Ps1))deallocate(WTRT2Ps1)
  if(allocated(WTRT1Ps1))deallocate(WTRT1Ps1)
  if(allocated(RTDP1s1))deallocate(RTDP1s1)
  if(allocated(RTWT1s1))deallocate(RTWT1s1)
  if(allocated(RTWT1Ns1))deallocate(RTWT1Ns1)
  if(allocated(RTWT1Ps1))deallocate(RTWT1Ps1)
  if(allocated(WTSTDGs1))deallocate(WTSTDGs1)
  if(allocated(WTSTDNs1))deallocate(WTSTDNs1)
  if(allocated(WTSTDPs1))deallocate(WTSTDPs1)
  if(allocated(CSNCs1))deallocate(CSNCs1)
  if(allocated(PSNCs1))deallocate(PSNCs1)
  if(allocated(ZSNCs1))deallocate(ZSNCs1)
  if(allocated(CFOPCs1))deallocate(CFOPCs1)
  if(allocated(CFOPNs1))deallocate(CFOPNs1)
  if(allocated(CFOPPs1))deallocate(CFOPPs1)
  if(allocated(PARs1))deallocate(PARs1)
  if(allocated(PARDIFs1))deallocate(PARDIFs1)
  if(allocated(VOLWMs1))deallocate(VOLWMs1)
  if(allocated(VOLPMs1))deallocate(VOLPMs1)
  if(allocated(TORTs1))deallocate(TORTs1)
  if(allocated(FILMs1))deallocate(FILMs1)
  if(allocated(ROXSKs1))deallocate(ROXSKs1)

  end subroutine DestructPlantAPIData

!------------------------------------------------------------------------
  subroutine plt_rad_init(this)
! DESCRIPTION
! initialize data type for plant_radiation_type
  implicit none
  class(plant_radiation_type) :: this

  allocate(this%ALBRs1(JP1))
  allocate(this%ALBPs1(JP1))
  allocate(this%TAUSs1(JC1+1))
  allocate(this%TAU0s1(JC1+1))
  allocate(this%THRM1s1(JP1))
  allocate(this%RADCs1(JP1))
  allocate(this%OMEGXs1(JSA1,JLI1,JLA1))
  allocate(this%OMEGAGs1(JSA1))
  allocate(this%OMEGAs1(JSA1,JLI1,JLA1))
  allocate(this%ZSINs1(JLI1))
  allocate(this%ZCOSs1(JLI1))
  allocate(this%IALBYs1(JSA1,JLI1,JLA1))
  allocate(this%RAD1s1(JP1))
  allocate(this%ABSRs1(JP1))
  allocate(this%ABSPs1(JP1))
  allocate(this%TAUPs1(JP1))
  allocate(this%TAURs1(JP1))
  allocate(this%RADPs1(JP1))
  allocate(this%FRADPs1(JP1))
  end subroutine plt_rad_init
!------------------------------------------------------------------------
  subroutine plt_rad_destroy(this)
! DESCRIPTION
! deallocate memory for plant_radiation_type
  implicit none
  class(plant_radiation_type) :: this

!  if(associated(this%ALBRs1))deallocate(this%ALBRs1)
!  if(associated(this%ALBPs1))deallocate(this%ALBPs1)
!  if(associated(this%TAUSs1))deallocate(this%TAUSs1)
!  if(associated(this%TAU0s1))deallocate(this%TAU0s1)
!  if(associated(this%THRM1s1))deallocate(this%THRM1s1)
!  if(associated(this%RADCs1))deallocate(this%RADCs1)
!  if(associated(this%OMEGXs1))deallocate(this%OMEGXs1)
!  if(associated(this%OMEGAGs1))deallocate(this%OMEGAGs1)
!  if(associated(this%ZSINs1))deallocate(this%ZSINs1)
!  if(allocated(ZCOSs1))deallocate(ZCOSs1)
!  if(associated(this%IALBYs1))deallocate(this%IALBYs1)
!  if(associated(this%OMEGAs1))deallocate(this%OMEGAs1)
!  if(associated(this%RAD1s1))deallocate(this%RAD1s1)
!  if(allocated(ABSRs1))deallocate(ABSRs1)
!  if(allocated(ABSPs1))deallocate(ABSPs1)
!  if(allocated(TAURs1))deallocate(TAURs1)
!  if(allocated(TAUPs1))deallocate(TAUPs1)
!  if(allocated(FRADPs1))deallocate(FRADPs1)
!  if(allocated(RADPs1))deallocate(RADPs1)

  end subroutine plt_rad_destroy

!------------------------------------------------------------------------

  subroutine plt_photo_init(this)
  class(plant_photosyns_type) :: this

  allocate(this%XKCO24s1(JP1))
  allocate(this%O2Ls1(JP1))
  allocate(this%SCO2s1(JP1))
  allocate(this%CO2Qs1(JP1))
  allocate(this%CHILLs1(JP1))
  allocate(this%ETMXs1(JP1))
  allocate(this%CHLs1(JP1))
  allocate(this%PEPCs1(JP1))
  allocate(this%CHL4s1(JP1))
  allocate(this%RUBPs1(JP1))
  allocate(this%VCMX4s1(JP1))
  allocate(this%VOMXs1(JP1))
  allocate(this%VCMXs1(JP1))
  allocate(this%XKCO2s1(JP1))
  allocate(this%XKO2s1(JP1))
  allocate(this%FDBKs1(JC1,JP1))
  allocate(this%FDBKXs1(JC1,JP1))
  allocate(this%CO2Ls1(JP1))
  allocate(this%XKCO2Ls1(JP1))
  allocate(this%XKCO2Os1(JP1))
  allocate(this%CO2Is1(JP1))
  allocate(this%O2Is1(JP1))
  allocate(this%RCSs1(JP1))
  allocate(this%FCO2s1(JP1))
  allocate(this%RSMHs1(JP1))

  end subroutine plt_photo_init
!------------------------------------------------------------------------

  subroutine plt_photo_destroy(this)
  class(plant_photosyns_type) :: this

!  if(allocated(XKCO24s1))deallocate(XKCO24s1)
!  if(allocated(O2Ls1))deallocate(O2Ls1)
!  if(allocated(SCO2s1))deallocate(SCO2s1)
!  if(allocated(CO2Qs1))deallocate(CO2Qs1)
!  if(allocated(CHILLs1))deallocate(CHILLs1)
!  if(allocated(ETMXs1))deallocate(ETMXs1)
!  if(allocated(CHLs1))deallocate(CHLs1)
!  if(allocated(PEPCs1))deallocate(PEPCs1)
!  if(allocated(CHL4s1))deallocate(CHL4s1)
!  if(allocated(RUBPs1))deallocate(RUBPs1)
!  if(allocated(VCMX4s1))deallocate(VCMX4s1)
!  if(allocated(VOMXs1))deallocate(VOMXs1)
!  if(allocated(VCMXs1))deallocate(VCMXs1)
!  if(allocated(XKCO2s1))deallocate(XKCO2s1)
!  if(allocated(XKO2s1))deallocate(XKO2s1)
!  if(allocated(FDBKs1))deallocate(FDBKs1)
!  if(allocated(FDBKXs1))deallocate(FDBKXs1)
!  if(allocated(CO2Ls1))deallocate(CO2Ls1)
!  if(allocated(XKCO2Ls1))deallocate(XKCO2Ls1)
!  if(allocated(XKCO2Os1))deallocate(XKCO2Os1)
!  if(allocated(CO2Is1))deallocate(CO2Is1)
!  if(allocated(O2Is1))deallocate(O2Is1)
!  if(allocated(RCSs1))deallocate(RCSs1)
!  if(allocated(FCO2s1))deallocate(FCO2s1)
!  if(allocated(RSMHs1)) deallocate(RSMHs1)
  end subroutine plt_photo_destroy

!------------------------------------------------------------------------
  subroutine plt_pheno_init(this)
  implicit none
  class(plant_pheno_type) :: this

  allocate(this%IFLGCs1(JP1))
  allocate(this%IDTHs1(JP1))
  allocate(this%RSETCs1(JP1))
  allocate(this%RSETNs1(JP1))
  allocate(this%RSETPs1(JP1))
  allocate(this%GROUPs1(JC1,JP1))
  allocate(this%IDTHPs1(JP1))
  allocate(this%ATRPs1(JC1,JP1))
  allocate(this%FLGZs1(JC1,JP1))
  allocate(this%DGSTGIs1(JC1,JP1))
  allocate(this%DGSTGFs1(JC1,JP1))
  allocate(this%GSTGIs1(JC1,JP1))
  allocate(this%GSTGFs1(JC1,JP1))
  allocate(this%FLG4s1(JC1,JP1))
  allocate(this%IWTYPs1(JP1))
  allocate(this%ISTYPs1(JP1))
  allocate(this%IBTYPs1(JP1))
  allocate(this%IDTHRs1(JP1))
  allocate(this%IDTYPs1(JP1))
  allocate(this%IPTYPs1(JP1))
  allocate(this%IFLGIs1(JP1))
  allocate(this%IGTYPs1(JP1))
  allocate(this%KVSTGs1(JC1,JP1))
  allocate(this%KVSTGNs1(JC1,JP1))
  allocate(this%IDAYs1(10,JC1,JP1))
  allocate(this%TGSTGIs1(JC1,JP1))
  allocate(this%TGSTGFs1(JC1,JP1))
  allocate(this%VSTGXs1(JC1,JP1))
  allocate(this%XRLAs1(JP1))
  allocate(this%XRNIs1(JP1))
  allocate(this%XDLs1(JP1))
  allocate(this%XPPDs1(JP1))
  allocate(this%IDTHBs1(JC1,JP1))
  allocate(this%IFLGPs1(JC1,JP1))
  allocate(this%IFLGFs1(JC1,JP1))
  allocate(this%IFLGEs1(JC1,JP1))
  allocate(this%IFLGAs1(JC1,JP1))
  allocate(this%IFLGGs1(JC1,JP1))
  allocate(this%IFLGRs1(JC1,JP1))
  allocate(this%IFLGQs1(JC1,JP1))
  allocate(this%VRNYs1(JC1,JP1))
  allocate(this%VRNZs1(JC1,JP1))
  allocate(this%VRNSs1(JC1,JP1))
  allocate(this%VRNLs1(JC1,JP1))
  allocate(this%VRNFs1(JC1,JP1))
  allocate(this%VRNXs1(JC1,JP1))

  end subroutine plt_pheno_init
!------------------------------------------------------------------------

  subroutine plt_pheno_destroy(this)
  implicit none
  class(plant_pheno_type) :: this

!  if(allocated(IFLGCs1))deallocate(IFLGCs1)
!  if(allocated(IDTHs1))deallocate(IDTHs1)
!  if(allocated(RSETCs1))deallocate(RSETCs1)
!  if(allocated(RSETNs1))deallocate(RSETNs1)
!  if(allocated(RSETPs1))deallocate(RSETPs1)
!  if(allocated(GROUPs1))deallocate(GROUPs1)
!  if(allocated(IDTHPs1))deallocate(IDTHPs1)
!  if(allocated(ATRPs1))deallocate(ATRPs1)
!  if(allocated(FLGZs1))deallocate(FLGZs1)
!  if(allocated(DGSTGIs1))deallocate(DGSTGIs1)
!  if(allocated(DGSTGFs1))deallocate(DGSTGFs1)
!  if(allocated(GSTGIs1))deallocate(GSTGIs1)
!  if(allocated(GSTGFs1))deallocate(GSTGFs1)
!  if(allocated(FLG4s1))deallocate(FLG4s1)
!  if(allocated(IWTYPs1))deallocate(IWTYPs1)
!  if(allocated(IBTYPs1))deallocate(IBTYPs1)
!  if(allocated(ISTYPs1))deallocate(ISTYPs1)
!  if(allocated(IDTHRs1))deallocate(IDTHRs1)
!  if(allocated(IDTYPs1))deallocate(IDTYPs1)
!  if(allocated(IPTYPs1))deallocate(IPTYPs1)
!  if(allocated(IGTYPs1))deallocate(IGTYPs1)
!  if(allocated(IFLGIs1))deallocate(IFLGIs1)
!  if(allocated(KVSTGNs1))deallocate(KVSTGNs1)
!  if(allocated(TGSTGIs1))deallocate(TGSTGIs1)
!  if(allocated(TGSTGFs1))deallocate(TGSTGFs1)
!  if(allocated(VSTGXs1))deallocate(VSTGXs1)
!  if(allocated(KVSTGs1))deallocate(KVSTGs1)
!  if(allocated(XDLs1))deallocate(XDLs1)
!  if(allocated(XRLAs1))deallocate(XRLAs1)
!  if(allocated(XRNIs1))deallocate(XRNIs1)
!  if(allocated(XPPDs1))deallocate(XPPDs1)
!  if(allocated(IDTHBs1))deallocate(IDTHBs1)
!  if(allocated(IFLGPs1))deallocate(IFLGPs1)
!  if(allocated(IFLGFs1))deallocate(IFLGFs1)
!  if(allocated(IFLGEs1))deallocate(IFLGEs1)
!  if(allocated(IFLGAs1))deallocate(IFLGAs1)
!  if(allocated(IFLGGs1))deallocate(IFLGGs1)
!  if(allocated(IFLGRs1))deallocate(IFLGRs1)
!  if(allocated(IFLGQs1))deallocate(IFLGQs1)
!  if(allocated(VRNYs1))deallocate(VRNYs1)
!  if(allocated(VRNZs1))deallocate(VRNZs1)
!  if(allocated(VRNSs1))deallocate(VRNSs1)
!  if(allocated(VRNLs1))deallocate(VRNLs1)
!  if(allocated(VRNFs1))deallocate(VRNFs1)
!  if(allocated(VRNXs1))deallocate(VRNXs1)
!  if(allocated(IDAYs1))deallocate(IDAYs1)

  end subroutine plt_pheno_destroy
!------------------------------------------------------------------------

  subroutine plt_morph_init(this)
  implicit none
  class(plant_morph_type) :: this


  allocate(this%MYs1(JP1))
  allocate(this%HTSTZs1(JP1))
  allocate(this%KLEAFs1(JC1,JP1))
  allocate(this%VSTGs1(JC1,JP1))
  allocate(this%NGs1(JP1))
  allocate(this%ZCs1(JP1))
  allocate(this%XTLIs1(JP1))
  allocate(this%ZLs1(0:JC1))
  allocate(this%ARSTPs1(JP1))
  allocate(this%ARLFPs1(JP1))
  allocate(this%NB1s1(JP1))
  allocate(this%NIXs1(JP1))
  allocate(this%NRTs1(JP1))
  allocate(this%NNODs1(JP1))
  allocate(this%NBTs1(JP1))
  allocate(this%NBRs1(JP1))
  allocate(this%NINRs1(JC1,JP1))
  allocate(this%PSTGs1(JC1,JP1))
  allocate(this%PSTGIs1(JC1,JP1))
  allocate(this%PSTGFs1(JC1,JP1))
  allocate(this%ANGBRs1(JP1))
  allocate(this%SURFs1(JLI1,JC1,JNODS1,JC1,JP1))
  allocate(this%KLEAFXs1(JC1,JP1))
  allocate(this%NBTBs1(JC1,JP1))
  allocate(this%GRNXBs1(JC1,JP1))
  allocate(this%HTSHEXs1(JC1,JP1))
  allocate(this%ARLFZs1(JC1,JP1))
  allocate(this%ARLFBs1(JC1,JP1))
  allocate(this%ANGSHs1(JP1))
  allocate(this%CLASSs1(JLI1,JP1))
  allocate(this%ARSTVs1(JC1,JP1))
  allocate(this%ARLFVs1(JC1,JP1))
  allocate(this%ARLF1s1(0:JNODS1,JC1,JP1))
  allocate(this%GRNOs1(JP1))
  allocate(this%SDPTHs1(JP1))
  allocate(this%SDPTHIs1(JP1))
  allocate(this%SDLGs1(JP1))
  allocate(this%SDVLs1(JP1))
  allocate(this%SDARs1(JP1))
  allocate(this%ARSTTs1(JC1))
  allocate(this%SSL1s1(JP1))
  allocate(this%SNL1s1(JP1))
  allocate(this%SLA1s1(JP1))
  allocate(this%ARLFTs1(JC1))
  allocate(this%ARLFSs1(JP1))
  allocate(this%HTNODXs1(0:JNODS1,JC1,JP1))
  allocate(this%HTSHEs1(0:JNODS1,JC1,JP1))
  allocate(this%HTNODEs1(0:JNODS1,JC1,JP1))
  allocate(this%SURFBs1(JLI1,JC1,JC1,JP1))
  allocate(this%ARLFLs1(JC1,0:JNODS1,JC1,JP1))
  allocate(this%ARSTKs1(JC1,JC1,JP1))
  allocate(this%NIs1(JP1))
  allocate(this%GRNOBs1(JC1,JP1))
  allocate(this%CFs1(JP1))

  end subroutine plt_morph_init

!------------------------------------------------------------------------

  subroutine plt_morph_destroy(this)
  implicit none
  class(plant_morph_type) :: this

!  if(allocated(MYs1)) deallocate(MYs1)
!  if(allocated(HTSTZs1))deallocate(HTSTZs1)
!  if(allocated(KLEAFs1))deallocate(KLEAFs1)
!  if(allocated(NGs1)) deallocate(NGs1)
!  if(allocated(ZCs1))deallocate(ZCs1)
!  if(allocated(XTLIs1))deallocate(XTLIs1)
!  if(allocated(ZLs1))deallocate(ZLs1)
!  if(allocated(ARSTPs1))deallocate(ARSTPs1)
!  if(allocated(ARLFPs1))deallocate(ARLFPs1)
!  if(allocated(NB1s1)) deallocate(NB1s1)
!  if(allocated(NIXs1))deallocate(NIXs1)
!  if(allocated(NRTs1)) deallocate(NRTs1)
!  if(allocated(NNODs1))deallocate(NNODs1)
!  if(allocated(NBTs1)) deallocate(NBTs1)
!  if(allocated(NBRs1)) deallocate(NBRs1)
!  if(allocated(NINRs1))deallocate(NINRs1)
!  if(allocated(PSTGs1))deallocate(PSTGs1)
!  if(allocated(PSTGIs1))deallocate(PSTGIs1)
!  if(allocated(PSTGFs1))deallocate(PSTGFs1)
!  if(allocated(ANGBRs1))deallocate(ANGBRs1)
!  if(allocated(SURFs1))deallocate(SURFs1)
!  if(allocated(KLEAFXs1))deallocate(KLEAFXs1)
!  if(allocated(NBTBs1))deallocate(NBTBs1)
!  if(allocated(GRNXBs1))deallocate(GRNXBs1)
!  if(allocated(HTSHEXs1))deallocate(HTSHEXs1)
!  if(allocated(ARLFZs1))deallocate(ARLFZs1)
!  if(allocated(ARLFBs1))deallocate(ARLFBs1)
!  if(allocated(ANGSHs1))deallocate(ANGSHs1)
!  if(allocated(CLASSs1))deallocate(CLASSs1)
!  if(allocated(ARSTVs1))deallocate(ARSTVs1)
!  if(allocated(ARLFVs1))deallocate(ARLFVs1)
!  if(allocated(ARLF1s1))deallocate(ARLF1s1)
!  if(allocated(GRNOs1))deallocate(GRNOs1)
!  if(allocated(SDPTHs1))deallocate(SDPTHs1)
!  if(allocated(SDPTHIs1))deallocate(SDPTHIs1)
!  if(allocated(SDLGs1))deallocate(SDLGs1)
!  if(allocated(SDVLs1))deallocate(SDVLs1)
!  if(allocated(SDARs1))deallocate(SDARs1)
!  if(allocated(ARSTTs1))deallocate(ARSTTs1)
!  if(allocated(SSL1s1))deallocate(SSL1s1)
!  if(allocated(SNL1s1))deallocate(SNL1s1)
!  if(allocated(SLA1s1))deallocate(SLA1s1)
!  if(allocated(ARLFTs1))deallocate(ARLFTs1)
!  if(allocated(ARLFSs1))deallocate(ARLFSs1)
!  if(allocated(HTNODXs1))deallocate(HTNODXs1)
!  if(allocated(HTSHEs1))deallocate(HTSHEs1)
!  if(allocated(HTNODEs1))deallocate(HTNODEs1)
!  if(allocated(SURFBs1))deallocate(SURFBs1)
!  if(allocated(ARLFLs1))deallocate(ARLFLs1)
!  if(allocated(ARSTKs1))deallocate(ARSTKs1)
!  if(allocated(NIs1)) deallocate(NIs1)
!  if(allocated(GRNOBs1))deallocate(GRNOBs1)
!  if(allocated(CFs1))  deallocate(CFs1)
!  if(allocated(VSTGs1))deallocate(VSTGs1)
  end subroutine plt_morph_destroy
end module PlantAPIData
