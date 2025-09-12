module PlantBGCPars

! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts, only : jroots
  use EcoSIMConfig, only : jcplxc
  implicit none
  public
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
!
! 
  real(r8) :: FracHour4LeafoffRemob(0:5)          !allocation parameter, [-]  
  real(r8) :: PART1X                              !minimum fraction of growth allocated to leaf, [-]
  real(r8) :: PART2X                              !minimum fraction of growth allocated to petiole, [-]
  real(r8) :: VMXC                                !rate constant for nonstructural C oxidation in respiration, [h-1]
  real(r8) :: RSpecLiterFall                      !rate constant for LitrFall at end of growing season, [h-1]
  real(r8) :: Hours4PhyslMature                   !number of hours with no grain filling required for physilogical maturity, [h]
  real(r8) :: Hours4FullSenes                     !number of hours until full senescence after physl maturity, [h]
  real(r8) :: XFRX                                !maximum storage C content for remobiln from stalk,root reserves, [gC]
  real(r8) :: XFRY                                !rate const for remobiln to storage from stalk,root reserves, [h-1]
  real(r8) :: FSNK                                !min ratio of branch or mycorrhizae to root for calculating C transfer, [-]
  real(r8) :: FXFS                                !rate constant for remobilization of stalk C,N,P, [h-1]
  real(r8) :: FMYC                                !rate constant for root-mycorrhizal C,N,P exchange, [h-1]
  real(r8) :: CNKI                                !nonstructural N inhibition constant on growth, [g N, g-1 C]
  real(r8) :: CPKI                                !nonstructural P inhibition constant on growth, [g P g-1 C]
  real(r8) :: RmSpecPlant                         !specific maintenance respiration rate, [g C g-1 N h-1]
  real(r8) :: TurgPSIMin4OrganExtens              !minimum tugor pressure for plant organ expansion,extension, [MPa]
  real(r8) :: RCMN                                !minimum stomatal resistance to CO2, [s m-1]
  real(r8) :: RTDPX                               !distance behind growing point for secondary roots, [m]
  real(r8) :: Root2ndTipLen4uptk                  !effective root tip length four resource uptake, not whole root hair 4 uptake, [m]
  real(r8) :: EMODR                               !root modulus of elasticity, [MPa]
  real(r8) :: QNTM                                !quantum efficiency, [umol e- umol-1 PAR]
  real(r8) :: CURV                                !shape parameter for e- transport response to PAR,[-]

  real(r8) :: CNKI_rubisco                        !nonstruct N inhibition constant on rubisco (g N g-1 C),1.0E+02_r8
  real(r8) :: CPKI_rubisco                        !nonstruct P inhibition constant on rubisco (g P g-1 C),1.0E+03_r8
  real(r8) :: RSMY_stomaCO2                       !minimum stomatal resistance for CO2 uptake (h m-1)
  real(r8) :: C4KI_pepcarboxy                     !nonstructural C inhibition constant on PEP carboxylase (uM)
  real(r8) :: Hours4ConiferSpringDeharden         !hours to full dehardening of conifers in spring (h)

  real(r8) :: ELEC3                               !e- requirement for CO2 fixn by rubisco,        [umol e- umol CO2]
  real(r8) :: ELEC4                               !e- requirement for CO2 fixn by PEP carboxylase,[umol e- umol CO2]
  real(r8) :: CO2KI                               ! Ki for C3 leakage from bundle sheath to mesophyll in C4, [uM]
  real(r8) :: FCMassCO2BundleSheath_node          !partition decarboxylation to CO2 in C4, [-]
  real(r8) :: FCMassHCO3BundleSheath_node         !partition leakage to HCO3 in C4, [-]
  real(r8) :: COMP4                               !C4 CO2 compensation point, [uM] 
  real(r8) :: FWCLeaf                             !leaf water content, (g H2O (gC leaf)-1)
  real(r8) :: FWCBundlSheath                      !leaf water content in bundle sheath, in C4 CO2 fixiation, [m3 H2O (gC)-1]
  real(r8) :: FWCMesophyll                        !leaf water content in mesophyll in C4 CO2 fixation, [m3 H2O (gC)-1]
  real(r8) :: ZPLFM                               !min N:C,P:C in leaves relative to max values from PFT file, [-]

  real(r8) :: ZPGRM                               !min N:C,P:C in grain relative to max values from PFT file,[-]
  real(r8) :: FSTK                                !fraction of stalk area contributing to water,heat flow
  real(r8) :: ZSTX                                !maximum stalk inner radius for tranpsiration, [m]
  real(r8) :: StalkMassDensity                    !stalk density, [MgC m-3]
  real(r8) :: SpecStalkVolume                     !specific volume (m3 gC-1)
  real(r8) :: FRTX                                !Fraction used to calculate woody faction of stalk,root,[-]
  real(r8) :: SETC                                !Km for nonstructural C concn on seed set, [g g-1]
  real(r8) :: SETN                                !Km for nonstructural N concn on seed set, [g g-1]
  real(r8) :: SETP                                !Km for nonstructural P concn on seed set, [g g-1]
  real(r8) :: SLA2                                !parameter for calculating leaf area expansion, [-]
  real(r8) :: SSL2                                !parameter for calculating petiole extension,[-]
  real(r8) :: SNL2                                !parameter for calculating stalk extension, [-]
  real(r8) :: CNMX                                !maximum C:N ratio for nonstructural N transfer, [-]
  real(r8) :: CPMX                                !maximum C:P ratio for nonstructural P transfer, [-]
  real(r8) :: CNMN                                !minimum C:N ratio for nonstructural N transfer, [-]
  real(r8) :: CPMN                                !minimum C:P ratio for nonstructural P transfer, [-]
  real(r8) :: EN2F                                !N fixation yield from C oxidation, [g N g-1 C]
  real(r8) :: VMXO                                !specific respiration rate by bacterial N2 fixers, [g g-1 h-1]
  real(r8) :: SPNDLK                              !half saturation parameter for nodule maintenance respiration, [gC]
  real(r8) :: SPNDL                               !specific decomposition rate by bacterial N2 fixers, [h-1]
  real(r8) :: CCNGR                               !parameter to calculate nonstructural C,N,P exchange,[-]
  real(r8) :: CCNGB                               !parameter to calculate nonstructural C,N,P exchange, [-]
  real(r8) :: NodulBiomCAtInfection               !initial bacterial mass at infection, [gC]
  real(r8) :: CZKM                                !Km for nonstructural Nuptake by bacteria, [gN] 
  real(r8) :: CPKM                                !Km for nonstructural P uptake by bacteria, [gP]
  real(r8) :: RCCZR                               !minimum fractions for root C recycling,[-]
  real(r8) :: RCCYR                               !maximum fractions for root C recycling,[-]
  real(r8) :: RCCXR                               !maximum fractions for root N recycling,[-]
  real(r8) :: RCCQR                               !maximum fractions for root P recycling,[-]
  real(r8) :: RCCZN                               !minimum fractions for bacteria C recycling,[-]
  real(r8) :: RCCYN                               !maximum fractions for bacteria C recycling,[-]
  real(r8) :: RCCXN                               !maximum fractions for bacteria N recycling,[-]
  real(r8) :: RCCQN                               !maximum fractions for bacteria P recycling,[-]

  integer :: HoursReq4LiterfalAftMature           !required hours after physl maturity until start of LitrFall, [h]
  REAL(R8) :: FRSV(0:3)                           !rate constant for remobiln of storage chemical element during leafout, [h-1]
  real(r8) :: FXFY(0:1)                           !rate constant for leaf-reserve nonstructural C exchange, [h-1]
  real(r8) :: FXFZ(0:1)                           !rate constant for leaf-reserve nonstructural N,P exchange, [h-1]
  real(r8) :: RateK4ShootSeaStoreNonstEXfer(0:3)  !rate constant for leaf-storage nonstructural chemical element exchange, [h-1]
  real(r8) :: RateK4RootSeaStorNonstEXfer(0:3)    !rate constant for root-storage nonstructural chemical element exchange, [h-1]
  real(r8) :: FXRT(0:1)                           !root partitioning of storage C during leafout,[-]
  real(r8) :: FXSH(0:1)                           !shoot partitioning of storage C during leafout,[-]
  real(r8) :: FXRN(6)                             !rate constant for plant-bacteria nonstructl C,N,P exchange,[h-1]
  REAL(R8) :: RCCX(0:3)                           !maximum fractions for shoot N recycling,[-]
  real(r8) :: RCCQ(0:3)                           !maximum fractions for shoot P recycling,[-]
  REAL(R8) :: RCCZ(0:3)                           !minimum fractions for shoot C recycling,[-]
  real(r8) :: RCCY(0:3)                           !maximum fractions for shoot C recycling,[-]
  real(r8) :: Hours4SenesAftMature(0:3)           !number of hours after physiol maturity required for senescence,[h]
  real(r8) :: HourReq2InitSStor4LeafOut(0:1)      !number of hours required to initiate remobilization of storage C for leafout, [h]
  real(r8) :: GVMX(0:1)                           !specific oxidation rate of nonstructural C during leafout at 25 C, [h]
  real(r8) :: RTSK(0:3)                           !relative primary root sink strength 0.25=shallow,4.0=deep root profile,[-]

  !terminate [label for varaible parsing]
  integer, parameter :: ibackward=1  !index for backward scattering in canopy radiation
  integer, parameter :: iforward =2  !index for forward scattering in canopy radiation

  real(r8) :: CURV2                               !2xCURV, [-]
  real(r8) :: CURV4                               !4XCURV, [-]
  real(r8) :: ZPLFD                               !1-ZPLFM, [-]
  real(r8) :: ZPGRD                               !1-ZPGRM, [-]

  type, public :: plant_bgc_par_type
  integer :: inonstruct                  !group id of plant nonstructural litter
  integer :: ifoliar                     !group id of plant foliar litter
  integer :: inonfoliar                  !group id of plant non-foliar litter group
  integer :: istalk                      !group id of plant stalk litter group
  integer :: iroot                       !group id of plant root litter
  integer :: icwood                      !group id of coarse woody litter
  integer :: NumGrowthStages             !number of growth stages
  integer :: MaxNumRootAxes              !maximum number of root axes
  integer :: JP1                         !number of plants
  integer :: MaxNumBranches              !maximum number of branches
  integer :: NumOfSkyAzimuthSects1       !number of sectors for the sky azimuth  [0,2*pi]
  integer :: jcplx                       !number of organo-microbial complexes
  integer :: NumOfLeafAzimuthSectors     !number of sectors for the leaf azimuth, [0,pi]
  integer :: NumCanopyLayers1          !number of canopy layers
  integer :: JZ1                         !number of soil layers
  integer :: NumLeafZenithSectors1     !number of sectors for the leaf zenith [0,pi/2]
  integer :: MaxNodesPerBranch1          !maximum number of canopy nodes, 25
  integer :: jsken                       !number of kinetic components in litter
  integer :: NumLitterGroups             !number of litter groups nonstructural(0,*)                         
  integer :: NumOfPlantMorphUnits        !number of plant organs
  integer :: NumOfPlantLitrCmplxs        !number of plant litter microbial-om complexes
  integer :: iprotein                    !kinetic id of litter component as protein
  integer :: icarbhyro                   !kinetic id of litter component as carbonhydrate
  integer :: icellulos                   !kinetic id of litter component as cellulose
  integer :: ilignin                     !kinetic id of litter component as lignin
  integer :: k_woody_litr                !woody litter complex id
  integer :: k_fine_litr                 !fine litter complex id
  integer :: jroots                      !number of root groups, plant root + myco types
  end type plant_bgc_par_type

  contains

!----------------------------------------------------------------------------------------------------
  subroutine InitPlantTraitTable(pltpar,NumGrowthStages,MaxNumRootAxes)
  use PlantTraitTableMod, only : AllocPlantTraitTable
  implicit none
  type(plant_bgc_par_type)  :: pltpar
  integer, intent(out) :: NumGrowthStages
  integer, intent(out) :: MaxNumRootAxes
  integer :: npft,nkopenclms,npfts_tab

  call InitVegPars(pltpar,npft,nkopenclms,npfts_tab)
  
  NumGrowthStages = pltpar%NumGrowthStages
  MaxNumRootAxes  = pltpar%MaxNumRootAxes

  call AllocPlantTraitTable(jroots,npft,nkopenclms,npfts_tab)

  end subroutine InitPlantTraitTable

!----------------------------------------------------------------------------------------------------
  subroutine InitVegPars(pltpar,npft,nkopenclms,npfts_tab)
  use EcoSIMCtrlMod, only : pft_file_in,pft_nfid
  use abortutils, only : endrun
  use fileUtil, only : file_exists
  use ncdio_pio
  implicit none
  type(plant_bgc_par_type)  :: pltpar
  integer, intent(out) :: npft        !total pft types, exclude koppen climate code
  integer, intent(out) :: nkopenclms  !
  integer, intent(out) :: npfts_tab   !total pft records, pft_short name + numerical koppen climate code

  if (len_trim(pft_file_in) == 0)then
    write(*,*) "Setting PFTs to one"
    npfts_tab=1
  else
    if(.not. file_exists(trim(pft_file_in)))then
      call endrun(msg='Fail to locate plant trait file '//trim(pft_file_in)//' in ' &
        //mod_filename,line=__LINE__)
    else
      npfts_tab  = get_dim_len(pft_file_in, 'npfts');
      npft       = get_dim_len(pft_file_in, 'npft')
      nkopenclms = get_dim_len(pft_file_in,'nkopenclms')
    endif
  endif  
  pltpar%inonstruct = 0
  pltpar%ifoliar    = 1
  pltpar%inonfoliar = 2
  pltpar%istalk     = 3
  pltpar%iroot      = 4
  pltpar%icwood     = 5
  pltpar%jcplx      = jcplxc
  pltpar%jroots=jroots
  FracHour4LeafoffRemob =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)
  PART1X                      = 0.05_r8
  PART2X                      = 0.02_r8
  VMXC                        = 0.015_r8
  RSpecLiterFall              = 2.884E-03_r8
  Hours4PhyslMature           = 168.0_r8
  Hours4FullSenes             = 240.0_r8
  XFRX                        = 2.5E-02_r8
  XFRY                        = 2.5E-03_r8
  FSNK                        = 0.05_r8
  FXFS                        = 1.0_r8
  FMYC                        = 0.1_r8
  CNKI                        = 1.0E-01_r8
  CPKI                        = 1.0E-02_r8
  RmSpecPlant                 = 0.010_r8
  TurgPSIMin4OrganExtens      = 0.1_r8
  RCMN                        = 1.560E+01_r8
  RTDPX                       = 0.00_r8
  Root2ndTipLen4uptk          = 2.0E-03_r8
  EMODR                       = 5.0_r8
  QNTM                        = 0.45_r8
  CURV                        = 0.70_r8
  CURV2                       = 2.0_r8*CURV
  CURV4                       = 4.0_r8*CURV
  ELEC3                       = 4.5_r8
  ELEC4                       = 3.0_r8
  CO2KI                       = 1.0E+03_r8
  FCMassCO2BundleSheath_node  = 0.02_r8
  FCMassHCO3BundleSheath_node = 1.0_r8-FCMassCO2BundleSheath_node
  COMP4                       = 0.5_r8
  FWCLeaf                        = 6.0_r8
  FWCBundlSheath                         = 0.2_r8*FWCLeaf
  FWCMesophyll                         = 0.8_r8*FWCLeaf
  ZPLFM                       = 0.33_r8
  ZPLFD                       = 1.0_r8-ZPLFM
  ZPGRM                       = 0.75_r8
  ZPGRD                       = 1.0_r8-ZPGRM

  CNKI_rubisco                = 1.0E-02_r8
  CPKI_rubisco                = 1.0E-03_r8
  RSMY_stomaCO2               = 2.78E-03_r8
  C4KI_pepcarboxy             = 5.0E+06_r8
  Hours4ConiferSpringDeharden = 276.9_r8


  FSTK                  = 0.05_r8
  ZSTX                  = 1.0E-03_r8
  StalkMassDensity      = 0.225_r8
  SpecStalkVolume       = 1.0E-06_r8/StalkMassDensity
  FRTX                  = 1.0_r8/(1.0_r8-(1.0_r8-FSTK)**2)
  SETC                  = 1.0E-02_r8
  SETN                  = 1.0E-03_r8
  SETP                  = 1.0E-04_r8
  SLA2                  = -0.33_r8
  SSL2                  = -0.50_r8
  SNL2                  = -0.67_r8
  CNMX                  = 0.20_r8
  CPMX                  = 0.020_r8
  CNMN                  = 0.050_r8
  CPMN                  = 0.005_r8
  EN2F                  = 0.20_r8
  VMXO                  = 0.125_r8
  SPNDLK                = 0.01_r8
  SPNDL                 = 5.0E-04_r8
  CCNGR                 = 2.5E-01_r8
  CCNGB                 = 6.0E-04_r8
  NodulBiomCatInfection = 1.0E-03_r8
  CZKM                  = 2.5E-03_r8
  CPKM                  = 2.5E-04_r8
  RCCZR                 = 0.056_r8
  RCCYR                 = 0.167_r8
  RCCXR                 = 0.833_r8
  RCCQR                 = 0.833_r8
  RCCZN                 = 0.167_r8
  RCCYN                 = 0.833_r8
  RCCXN                 = 0.833_r8
  RCCQN                 = 0.833_r8

  HoursReq4LiterfalAftMature = 960


  RCCZ=real((/0.167,0.167,0.167,0.056/),r8)
  RCCY=real((/0.333,0.333,0.167,0.333/),r8)
  RCCX=real((/0.417,0.833,0.833,0.833/),r8)
  RCCQ=real((/0.417,0.833,0.833,0.833/),r8)

  RTSK=real((/0.25,1.0,4.0,4.0/),r8)
  FXRN=real((/0.25,0.125,0.0625,0.225,0.075,0.025/),r8)
  RateK4ShootSeaStoreNonstEXfer=real((/1.0E-02,1.0E-02,1.0E-05,5.0E-05/),r8)
  RateK4RootSeaStorNonstEXfer=real((/1.0E-02,1.0E-02,1.0E-05,5.0E-05/),r8)
  FXSH=real((/0.50,0.75/),r8);FXRT=real((/0.50,0.25/),r8)
  FRSV=real((/0.025,0.025,0.001,0.001/),r8)
  FXFY=real((/0.025,0.005/),r8);FXFZ=real((/0.25,0.05/),r8)
  Hours4SenesAftMature=real((/360.0,1440.0,720.0,720.0/),r8)
  HourReq2InitSStor4LeafOut=real((/68.96,276.9/),r8);GVMX=real((/0.010,0.0025/),r8)
  
  end subroutine InitVegPars

end module PlantBGCPars
