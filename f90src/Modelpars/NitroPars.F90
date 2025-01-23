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
!
! ORAD=microbial radius (m), BIOS=microbial density (n m-3)
! BIOA=microbial surface area (m2 m-3), DCKI=inhibition of
! decomposition by microbial concentration (g C m-3)
! RCCX=maximum remobilization of microbial N (-)
! RCCY=maximum remobilization of microbial P (-)
!     RCCZ, RCCY = minimum, maximum remobilization of microbial C (-)
!     FPRIM, FPRIMM=fraction of nonstructural, microbial C,N,P
!     transferred with priming (-), OMGR=rate constant for
!     transferring nonstructural to structural microbial C (h-1)
!     OQKI=DOC product inhibition constant for decomposition (g C m-3)
!     H2KI=H2 product inhibition for methanogenesis (g H m-3)
!     COMKI, COMKM= Km to slow microbial decomposition, maintenance
!     respiration with low microbial C (g micr C g-1 subs C)
!     CKC=controls C remobilization of microbial C (g C g-1 C)
!     FOSCZ0, FOSCZL=rate constants for mixing surface (0) and
!     subsurface (L) litter (h-1),FMN=minimum ratio of total
!     biological demand for any substrate by any microbial population
!     DCKM0, DCKML=Km for SOC decomposition (g C g-1 soil)
!


!
!     SPECIFIC RESPIRATION RATES, M-M UPTAKE CONSTANTS,
!     STOICHIOMETRIC CONSTANTS FOR MICROBIAL REDOX REACTIONS
!
!     VMX*=specific oxidation rates (g C g-1C h-1_
!        O=all bacteria, F=fungi, M=acetotrophic methanogens
!        H=ammonia oxidizers, N=nitrite oxidizers, 4=methanotrophs
!        C=hydrogenotrophic methanogens
!     OQK*=Km for DOC uptake by heterotrophs (g C m-3)
!        M=all bacteria and fungi, A=acetate by fermenters
!        AM=acetate by acetotrophic methanogens
!     CCKM=Km for CO2 uptake, CCK4=Km for CH4 uptake (g C m-3)
!     Z*KM=Km for N uptake (g N m-3)
!        H=NH4 by nitrifiers, N=NO2 by nitrifiers
!        3=NO3 by denitrifiers, 2=NO2 by denitrifiers
!        1=N2O uptake by denitrifiers
!     Z4*=NH4 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
!       MX=maximum uptake rate, KU=Km, MN= minimum concentration
!     ZO*=NO3 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
!       MX=maximum uptake rate, KU=Km, MN= minimum concentration
!     HP*=H2PO4 uptake kinetics by all MFTs(g P m-2 h-1, g P m-3)
!       MX=maximum uptake rate, KU=Km, MN= minimum concentration
!     ZFKM=Km for N2 uptake by diazotrophs (g N m-3)
!     H2KM=Km for H2 uptake by hydrogenotrophic methanogens (g H m-3)
!     ECNH=efficiency CO2 conversion to biomass by ammonia oxidizers
!     ECNO=efficiency CO2 conversion to biomass by nitrite oxidizers
!     ECHO=efficiency CO2 conversion to biomass by methane oxidizers
!     eQNO3toOxy,eQNO2toOxy,eQN2OtoOxy=N2:O2 ratios for e- transfers to NO3, NO2 and N2O
!     by denitrifiers, RNFNI=parameter for nitrification inhibition
!     ZHKI=inhibition of nitrification inhibition by NH3 (g N m-3)
!     VMKI=product inhibn for NOx reduction by denitrifiers(g N m-3)
!     VHKI=product inhibn for NH3 oxidation by nitrifiers (g N m-3)
!     OXKA=Km for O2 uptake by nitrifiers(g O m-3)
!
!     ENERGY REQUIREMENTS FOR MICROBIAL GROWTH AND
!     ENERGY YIELDS FROM REDUCTION OF O2, OC, CH4, NO3, N2
!
!     EOM*=energy requirements for microbial growth (kJ g-1 C)
!        C=aerobic bacteria, D=denitrifiers, G=fungi, F=fermenters
!        H=methanogens, N=diazotrophs
!     G*=free energy yields of redox reactions (kJ g-1 C or N)
!        O2X=DOC-CO2, H4X=CO2-CH4, CHX=DOC-acetate, O2A=acetate-CO2
!        C4X=acetate-CH4, COX=CO2-CH4, NOX=NO3-NO2,NO2-N2O,N2O-N2
!        N2X=N2-NH3
!     E*=growth respiration efficiency (-)(growth yield=1.0-E*)
!        N2X=aerobic N2 fixation, N2Y=anaerobic N2 fixation
!        O2X=aerobic bacteria (DOC), H4X=fermenters, O2G=fungi
!        O2D=denitrifiers (aerobic), NFX=diazotrophs
!        NOX= denitrifiers (anaerobic),O2A=aerobic bacteria (acetate)
!
!
!     SORPTION COEFFICIENTS
!
!
!     SPECIFIC DECOMPOSITION RATES
!
!     DOSA=rate constant for litter colonization by heterotrophs (g C g-1 C)
!     SP*= specific decomposition rate constant (g subs. C g-1 micr. C)
!       OHC=adsorbed SOC, OHA=adsorbed acetate, OSC=SOC
!       (K=0,M=1,4 woody litter, K=1,M=1,4 non-woody litter,
!       K=2,M=1,4 manure, K=3,M=1,1 POC, K=4,M=1,2 humus)
!       ORC (M=1,2) microbial residue, OMC (M=1,2) microbial biomass
!     RMOM=specific maintenance respiration (g C g-1 N h-1)
!

  real(r8) :: ORAD
  real(r8) :: BIOS
  real(r8) :: BIOA
  real(r8) :: DCKI
  real(r8) :: RCCX
  real(r8) :: RCCQ
  real(r8) :: RCCZ
  real(r8) :: RCCY
  real(r8) :: FPRIM
  real(r8) :: FPRIMM
  real(r8) :: OMGR
  real(r8) :: OQKI
  real(r8) :: H2KI
  real(r8) :: OAKI
  real(r8) :: COMKI
  real(r8) :: COMKM
  real(r8) :: CKC
  real(r8) :: FOSCZ0   !rate for mixing surface litter  [1/h]
  real(r8) :: FOSCZL   !rate for mixing subsurface litter [1/h]
  real(r8) :: FMN
  real(r8) :: DCKM0
  real(r8) :: DCKML
  real(r8) :: VMXO
  real(r8) :: VMXF
  real(r8) :: VMXM
  real(r8) :: VMXH
  real(r8) :: VMXN
  real(r8) :: VMX4
  real(r8) :: VMXC
  real(r8) :: OQKM
  real(r8) :: OQKA
  real(r8) :: OQKAM
  real(r8) :: CCKM
  real(r8) :: CCK4
  real(r8) :: ZHKM
  real(r8) :: ZNKM
  real(r8) :: Z3KM
  real(r8) :: Z2KM
  real(r8) :: Z1KM
  real(r8) :: Z4MX
  real(r8) :: Z4KU
  real(r8) :: Z4MN
  real(r8) :: ZOMX
  real(r8) :: ZOKU
  real(r8) :: ZOMN
  real(r8) :: HPMX
  real(r8) :: HPKU
  real(r8) :: HPMN
  real(r8) :: ZFKM
  real(r8) :: H2KM
  real(r8) :: ECNH
  real(r8) :: ECNO
  real(r8) :: ECHO
  real(r8) :: eQNO3toOxy
  real(r8) :: eQNO2toOxy
  real(r8) :: eQN2OtoOxy
  real(r8) :: RNFNI
  real(r8) :: ZHKI
  real(r8) :: VMKI
  real(r8) :: VHKI
  real(r8) :: OXKA
  real(r8) :: EOMC
  real(r8) :: EOMD
  real(r8) :: EOMG
  real(r8) :: EOMF
  real(r8) :: EOMH
  real(r8) :: EOMN
  real(r8) :: GO2X
  real(r8) :: GH4X
  real(r8) :: GCHX
  real(r8) :: GO2A
  real(r8) :: GC4X
  real(r8) :: GCOX
  real(r8) :: GNOX
  real(r8) :: GN2X
  real(r8) :: EN2X
  real(r8) :: EN2Y
  real(r8) :: EO2X
  real(r8) :: EH4X
  real(r8) :: EO2G
  real(r8) :: EO2D
  real(r8) :: ENFX
  real(r8) :: ENOX
  real(r8) :: EO2A
  real(r8) :: TSORP
  real(r8) :: HSORP
  real(r8) :: SPOHC
  real(r8) :: SPOHA
  real(r8) :: RMOM
  real(r8) :: SPORC(2)
  real(r8) :: SPOMC(2)
  real(r8) :: EN2F(7)
  real(r8) :: EFIRE(2,21:22)
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
  OMGR   = 2.5E-01_r8
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
