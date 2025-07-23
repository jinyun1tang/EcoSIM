module NitroPars
!
! DESCRIPTION:
! code defining parameters for nitro
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst
  implicit none
  public
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  save
! SUBSTRATE DECOMPOSITION BY MICROBIAL POPULATIONS

  real(r8) :: ORAD         !microbial radius, [m]
  real(r8) :: BIOS         !microbial density, [n m-3]
  real(r8) :: BIOA         !microbial surface area, [m2 m-3]
  real(r8) :: DCKI         !inhibition of decomposition by microbial concentration, [g C m-3]
  real(r8) :: RCCX         !maximum remobilization of microbial N, [-]
  real(r8) :: RCCQ         !maximum P recycling fractions , [-]
  real(r8) :: RCCZ         !minimum fractions for bacteria C recycling,[-]
  real(r8) :: RCCY         !maximum remobilization of microbial P, [-]
  real(r8) :: FPRIM        !fraction of nonstructural transferred with priming, [-] 
  real(r8) :: FPRIMM       !fraction of microbial C,N,P transferred with priming, [-]
  real(r8) :: OMGR         !rate constant for transferring nonstructural to structural microbial C,[h-1]
  real(r8) :: OQKI         !DOC product inhibition constant for decomposition, [g C m-3]
  real(r8) :: H2KI         !H2 product inhibition for methanogenesis, [g H m-3]
  real(r8) :: OAKI         !acetate product inhibition constant for decomposition, [g C m-3]
  real(r8) :: COMKI        !Km to slow microbial decomposition with low microbial C, [g micr C g-1 subs C]
  real(r8) :: COMKM        !Km to slow microbial maintenance respiration with low microbial C, [g micr C g-1 subs C]
  real(r8) :: CKC          !controls C remobilization of microbial C, [g C g-1 C]
  real(r8) :: FOSCZ0       !rate for mixing surface litter,  [h-1]
  real(r8) :: FOSCZL       !rate for mixing subsurface litter, [h-1]
  real(r8) :: FMN          !minimum ratio of total biological demand for any substrate by any microbial population,[-]
  real(r8) :: DCKM0        !Km for SOC decomposition, [g C g-1 soil]
  real(r8) :: DCKML        !Km for SOC decomposition, [g C g-1 soil] 
  real(r8) :: VMXO         !specific oxidation rates for all bacteria, [g C g-1C h-1]
  real(r8) :: VMXF         !specific oxidation rates for all fungi, [g C g-1C h-1] 
  real(r8) :: VMXM         !specific oxidation rates for acetotrophic methanogens,[g C g-1C h-1] 
  real(r8) :: VMXH         !specific oxidation rates for ammonia oxidizers, [g  g-1C h-1] 
  real(r8) :: VMXN         !specific oxidation rates for nitrite oxidizers, [g  g-1C h-1] 
  real(r8) :: VMX4         !specific oxidation rates for methanotrophs, [g   g-1C h-1] 
  real(r8) :: VMXC         !specific oxidation rates for hydrogenotrophic methanogens,[g   g-1C h-1] 
  real(r8) :: OQKM         !Km for DOC uptake by heterotrophs bacteria and fungi, [g C m-3]
  real(r8) :: OQKA         !Km for acetate uptake by heterotrophic fermenters, [g C m-3]
  real(r8) :: OQKAM        !Km for acetate uptake  by acetotrophic methanogens,[g C m-3]
  real(r8) :: CCKM         !Km for CO2 uptake,  [g C m-3]
  real(r8) :: CCK4         !Km for CH4 uptake, [g C m-3]
  real(r8) :: ZHKM         !Km for NH4 uptake by nitrifiers, [gN m-3]
  real(r8) :: ZNKM         !Km for NO2 uptake by nitrifiers, [gN m-3]
  real(r8) :: Z3KM         !Km for NO3 uptake by denitrifiers,[gN m-3]
  real(r8) :: Z2KM         !Km for NO2 uptake by denitrifiers, [gN m-3]
  real(r8) :: Z1KM         !Km for N2O uptake by denitrifiers, [gN m-3]
  real(r8) :: Z4MX         !maximum uptake rate for NH4 uptake kinetics by all microbial functional groups,[g N m-2 h-1]
  real(r8) :: Z4KU         !Km for NH4 uptake kinetics by all microbial functional groups,[g N m-3]
  real(r8) :: Z4MN         !minimum concentration for NH4 uptake kinetics by all microbial functional groups,[g N m-3]
  real(r8) :: ZOMX         !maximum uptake rate for NO3 uptake kinetics by all microbial functional groups, [g N m-2 h-1] 
  real(r8) :: ZOKU         !Km for NO3 uptake kinetics by all microbial functional groups, [g N m-3]
  real(r8) :: ZOMN         !minimum concentration for NO3 uptake kinetics by all microbial functional groups, [g N m-3]
  real(r8) :: HPMX         !maximum rate for H2PO4 uptake kinetics by all microbial functional groups, [g P m-2 h-1]
  real(r8) :: HPKU         !Km for H2PO4 uptake kinetics by all microbial functional groups,[g P m-3]
  real(r8) :: HPMN         !Minimum concentration for H2PO4 uptake kinetics by all microbial functional groups,[g P m-3]
  real(r8) :: ZFKM         !Km for N2 uptake by diazotrophs, [g N m-3]
  real(r8) :: H2KM         !Km for H2 uptake by hydrogenotrophic methanogens, [g H m-3]
  real(r8) :: ECNH         !efficiency of CO2 conversion to biomass by ammonia oxidizers,[-]
  real(r8) :: ECNO         !efficiency of CO2 conversion to biomass by nitrite oxidizers,[-]
  real(r8) :: ECHO         !efficiency of CO2 conversion to biomass by methane oxidizers,[-]
  real(r8) :: eQNO3toOxy   !N2:O2 ratios for e- transfers to NO3 by denitrifiers,[-]
  real(r8) :: eQNO2toOxy   !N2:O2 ratios for e- transfers to NO2 by denitrifiers,[-]
  real(r8) :: eQN2OtoOxy   !N2:O2 ratios for e- transfers to N2O by denitrifiers,[-]
  real(r8) :: RNFNI        !parameter for nitrification inhibition,[-]
  real(r8) :: ZHKI         !inhibition of nitrification inhibition by NH3, [g N m-3]
  real(r8) :: VMKI         !product inhibn for NOx reduction by denitrifiers, [g N m-3]
  real(r8) :: VHKI         !product inhibn for NH3 oxidation by nitrifiers, [g N m-3]
  real(r8) :: OXKA         !Km for O2 uptake by nitrifiers, [g O m-3]
  real(r8) :: EOMC         !energy requirements for microbial growth of aerobic bacteria, [kJ g-1 C]
  real(r8) :: EOMD         !energy requirements for microbial growth of denitrifiers, [kJ g-1 C]
  real(r8) :: EOMG         !energy requirements for microbial growth of fungi, [kJ g-1 C]
  real(r8) :: EOMF         !energy requirements for microbial growth of fermenters, [kJ g-1 C]
  real(r8) :: EOMH         !energy requirements for microbial growth of methanogens, [kJ g-1 C]
  real(r8) :: EOMN         !energy requirements for microbial growth of diazotrophs, [kJ g-1 C]
  real(r8) :: GO2X           !free energy yields of redox reactions for DOC-CO2, [kJ g-1 C]
  real(r8) :: GH4X           !free energy yields of redox reactions for CO2-CH4, [kJ g-1 C]
  real(r8) :: GCHX           !free energy yields of redox reactions for DOC-acetate, [kJ g-1 C]
  real(r8) :: GO2A           !free energy yields of redox reactions for acetate-CO2, [kJ g-1 C]
  real(r8) :: GC4X           !free energy yields of redox reactions for acetate-CH4, [kJ g-1 C]
  real(r8) :: GCOX           !free energy yields of redox reactions for CO2-CH4, [kJ g-1 C]
  real(r8) :: GNOX           !free energy yields of redox reactions for NO3-NO2, NO2-N2O,N2O-N2, [kJ g-1 N]
  real(r8) :: GN2X           !free energy yields of redox reactions for N2-NH3, [kJ g-1 N]
  real(r8) :: EN2X           !growth respiration efficiency for aerobic N2 fixation, [-]
  real(r8) :: EN2Y           !growth respiration efficiency for anaerobic N2 fixation, [-]
  real(r8) :: EO2X           !growth respiration efficiency for aerobic bacteria (DOC), [-]
  real(r8) :: EH4X           !growth respiration efficiency for fermenters, [-]
  real(r8) :: EO2G           !growth respiration efficiency for fungi, [-]
  real(r8) :: EO2D           !growth respiration efficiency for denitrifiers (aerobic), [-]
  real(r8) :: ENFX           !growth respiration efficiency for diazotrophs, [-]
  real(r8) :: ENOX           !growth respiration efficiency for denitrifiers (anaerobic), [-]
  real(r8) :: EO2A           !growth respiration efficiency for aerobic bacteria (acetate), [-]
  real(r8) :: TSORP          !sorption rate constant for OHC, [h-1]
  real(r8) :: HSORP          !sorption rate coefficient for OHC, [-]
  real(r8) :: SPOHC          !specific decomposition rate constant for adsorbed SOC, [g subs. C g-1 micr. C]
  real(r8) :: SPOHA          !specific decomposition rate constant for adsorbed acetate, [g subs. C g-1 micr. C]
  real(r8) :: RMOM           !specific maintenance respiration, [g C g-1 N h-1]
  real(r8) :: SPORC(2)       !specific decomposition rate constant microbial residue,  [g C g-1 N h-1]
  real(r8) :: SPOMC(2)       !specific decomposition rate constant microbial biomass,  [g C g-1 N h-1]
  real(r8) :: EN2F(7)        !N fixation yield from C oxidation, [g N g-1 C]
  real(r8) :: EFIRE(2,21:22) !partition coefficient for N loss as NH3 and P loss as PO4 during combustion, [g gC-1]
  contains

  subroutine initNitroPars
  implicit none

  ORAD   = ppmc
  BIOS   = ppmc/(4.19_r8*ORAD**3)
  BIOA   = BIOS*12.57_r8*ORAD**2
  DCKI   = 2.5_r8
  RCCX   = 0.833_r8
  RCCQ   = 0.833_r8
  RCCZ   = 0.167_r8
  RCCY   = 0.833_r8
  FPRIM  = 5.0E-02_r8
  FPRIMM = ppmc
  OMGR   = 2.5E-01_r8  !
  OQKI   = 1.2E+03_r8
  H2KI   = 1.0_r8
  OAKI   = 12.0_r8
  COMKI  = 1.0E-03_r8
  COMKM  = 1.0E-04_r8
  CKC    = 1.0E-03_r8
  FOSCZ0 = 2.0E-02_r8
  FOSCZL = 2.0E-06_r8
  FMN    = 1.0E-03_r8
  DCKM0  = 1.0E+03_r8
  DCKML  = 1.0E+03_r8
  VMXO   = 0.125_r8
  VMXF   = 0.125_r8
  VMXM   = 0.125_r8  !acetoclastic methanogenesis
  VMXH       = 0.375_r8
  VMXN       = 0.25_r8
  VMX4       = 0.375_r8
  VMXC       = 0.125_r8  !hydrogenotrophic methanogenesis
  OQKM       = 1.2E+01_r8
  OQKA       = 1.2E+01_r8
  OQKAM      = 1.2E+01_r8
  CCKM       = 0.15_r8
  CCK4       = 1.2E-03_r8
  ZHKM       = 1.4_r8
  ZNKM       = 1.4_r8
  Z3KM       = 1.4_r8
  Z2KM       = 1.4_r8
  Z1KM       = 0.014_r8
  Z4MX       = 5.0E-03_r8
  Z4KU       = 0.40_r8
  Z4MN       = 0.0125_r8
  ZOMX       = 5.0E-03_r8
  ZOKU       = 0.35_r8
  ZOMN       = 0.03_r8
  HPMX       = 1.0E-03_r8
  HPKU       = 0.075_r8
  HPMN       = 0.002_r8
  ZFKM       = 0.14_r8
  H2KM       = 0.01_r8
  ECNH       = 0.30_r8
  ECNO       = 0.10_r8
  ECHO       = 0.75_r8
  eQNO3toOxy = 0.857_r8
  eQNO2toOxy = 0.857_r8
  eQN2OtoOxy = 0.429_r8
  RNFNI      = 2.0E-04_r8
  ZHKI       = 7.0E+03_r8
  VMKI       = 0.25_r8
  VHKI       = 15.0_r8
  OXKA       = 0.16_r8
  EOMC       = 25.0_r8
  EOMD       = 37.5_r8
  EOMG       = 37.5_r8
  EOMF       = 75.0_r8
  EOMH       = 25.0_r8
  EOMN       = 75.0_r8
  GO2X       = 37.5_r8
  GH4X       = 66.5_r8
  GCHX       = 4.50_r8
  GO2A       = GO2X-GCHX
  GC4X       = 3.00_r8
  GCOX       = 11.00_r8
  GNOX       = 10.0_r8
  GN2X       = 187.5_r8
  EN2X       = GO2X/GN2X
  EN2Y       = GCHX/GN2X
  EO2X       = 1.0_r8/(1.0_r8+GO2X/EOMC)
  EH4X  = 1.0_r8/(1.0_r8+GH4X/EOMC)
  EO2G  = 1.0_r8/(1.0_r8+GO2X/EOMG)
  EO2D  = 1.0_r8/(1.0_r8+GO2X/EOMD)
  ENFX  = 1.0_r8/(1.0_r8+GO2X/EOMN)
  ENOX  = 1.0_r8/(1.0_r8+GNOX/EOMC)
  EO2A  = 1.0_r8/(1.0_r8+GO2A/EOMC)
  TSORP = 0.5_r8
  HSORP = 1.0_r8
  SPOHC = 0.25_r8
  SPOHA = 0.25_r8
  RMOM  = 0.010_r8

  SPORC = (/7.5_r8,1.5_r8/)                               !hydrolysis of microbial residue
  SPOMC = (/1.0E-02_r8,0.1E-02_r8/)                       !basal mortality rates
  EN2F  = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,EN2X,EN2Y/)
  EFIRE = reshape((/1.0_r8,1.0_r8,0.917_r8,0.167_r8/),shape(EFIRE))

  end subroutine initNitroPars

end module NitroPars
