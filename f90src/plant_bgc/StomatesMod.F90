  module stomatesMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use minimathmod
  use PlantAPIData
  implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  real(r8), PARAMETER :: QNTM=0.45_r8       !quantum efficiency (umol e- umol-1 PAR)
  real(r8), parameter :: CURV=0.70_r8       !shape parameter for e- transport response to PAR
  real(r8), PARAMETER :: CURV2=2.0_r8*CURV  !
  real(r8), PARAMETER :: CURV4=4.0_r8*CURV
  real(r8), PARAMETER :: ELEC3=4.5_r8       !e- requirement for CO2 fixed by rubisco (umol e- umol CO2)
  real(r8), PARAMETER :: ELEC4=3.0_r8       !e- requirement for CO2 fixed by PEP carboxylase (umol e- umol CO2)
  real(r8), PARAMETER :: CNKI=1.0E+02_r8    !nonstruct N inhibition constant on rubisco (g N g-1 C)
  real(r8), PARAMETER :: CPKI=1.0E+03_r8    !nonstruct P inhibition constant on rubisco (g P g-1 C)
  real(r8), PARAMETER :: RSMY=2.78E-03_r8   !minimum stomatal resistance for CO2 uptake (h m-1)
  real(r8), PARAMETER :: ATRPZ=276.9_r8     !hours to full dehardening of conifers in spring (h)
  real(r8), PARAMETER :: COMP4=0.5_r8       !C4 CO2 compensation point (uM)
  real(r8), PARAMETER :: FDML=6.0_r8        !leaf water content (g H2O g-1 C)
  real(r8), PARAMETER :: FBS=0.2_r8*FDML    !leaf water content in bundle sheath in C4 CO2 fixation
  real(r8), PARAMETER :: FMP=0.8_r8*FDML    !leaf water content in mesophyll in C4 CO2 fixationn
  real(r8), PARAMETER :: C4KI=5.0E+06_r8    !nonstructural C inhibition constant on PEP carboxylase (uM)
  real(r8), parameter :: Hours2KillAnuals(0:5)=real((/336.0,672.0,672.0,672.0,672.0,672.0/),r8)  !number of hours with no grain fill to terminate annuals

  public :: stomates
  contains

  subroutine stomates(I,J,NZ)
!
!     THIS subroutine CALCULATES CANOPY STOMATAL RESISTANCE AT MAXIMUM
!     CANOPY TURGOR FOR USE IN ENERGY BALANCE EQUATIONS IN 'UPTAKE'
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NZ

  integer :: K,L,M,NB,N
  REAL(R8):: CanopyBndlResist_pft4CO2
  real(r8):: RI
!     begin_execution
  associate(                            &
    RIB      =>  plt_ew%RIB       , &
    RAZ      =>  plt_ew%RAZ       , &
    TairK      =>  plt_ew%TairK       , &
    TKCanopy_pft     =>  plt_ew%TKCanopy_pft      , &
    CO2E     =>  plt_site%CO2E    , &
    CanopyGasCO2_pft     =>  plt_photo%CanopyGasCO2_pft   , &
    CanopyLeafArea_pft    =>  plt_morph%CanopyLeafArea_pft  , &
    ZEROP    =>  plt_biom%ZEROP   , &
    CNETX    =>  plt_bgcr%CNETX   , &
    SineSolarAngle     =>  plt_rad%SineSolarAngle     , &
    AirConc_pft     =>  plt_photo%AirConc_pft   , &
    MinCanPStomaResistH2O_pft     =>  plt_photo%MinCanPStomaResistH2O_pft   , &
    CanPCi2CaRatio     =>  plt_photo%CanPCi2CaRatio   , &
    MaxCanPStomaResistH2O_pft     =>  plt_photo%MaxCanPStomaResistH2O_pft   , &
    LeafIntracellularCO2_pft    =>  plt_photo%LeafIntracellularCO2_pft    &
  )
!
!     CANOPY TEMPERATURE + OFFSET FOR THERMAL ADAPTATION FROM 'READQ'
!
!     CANOPY BOUNDARY LAYER RESISTANCE
!
!     RI=Richardson's number
!     RIB=canopy isothermal Richardson's number
!     TairK,TKCanopy_pft=air,canopy temperature
!     RAZ=canopy isothermal boundary later resistance
!     CanopyBndlResist_pft4CO2=canopy boundary layer resistance to CO2, h/m
!     AirConc_pft=number of moles of air per m3
!
  RI=RichardsonNumber(RIB,TairK,TKCanopy_pft(NZ))

  CanopyBndlResist_pft4CO2=1.34*AMAX1(5.56E-03,RAZ(NZ)/(1.0-10.0*RI))

  !assuming pressure is one atmosphere
  AirConc_pft(NZ)=GetMolAirPerm3(TKCanopy_pft(NZ))
!
!     CANOPY CO2 CONCENTRATION FROM CO2 INFLUXES AND EFFLUXES
!
!     CanopyGasCO2_pft,CO2E=CO2 concentrations in canopy air,atmosphere, umol mol-1 (ppmv)
!     CNETX=net CO2 flux in canopy air from soil,plants, g d-2 h-1
! assuming steady state, canopy CO2 concentration is computed with mass balance. 
! how 8.33E+04 is determined. 
  CanopyGasCO2_pft(NZ)=CO2E-8.33E+04_r8*CNETX*CanopyBndlResist_pft4CO2/AirConc_pft(NZ)
  CanopyGasCO2_pft(NZ)=AMIN1(CO2E+200.0_r8,AZMAX1(CO2E-200.0_r8,CanopyGasCO2_pft(NZ)))
!
!     MESOPHYLL CO2 CONCENTRATION FROM CI:CA RATIO ENTERED IN 'READQ'
!
!     LeafIntracellularCO2_pft=intercellular CO2 concentration
!     CanPCi2CaRatio=intercellular:atmospheric CO2 concn ratio from PFT file, parameter
!     SineSolarAngle=sine of solar angle
!     CanopyLeafArea_pft=PFT leaf area
!
  LeafIntracellularCO2_pft(NZ)=CanPCi2CaRatio(NZ)*CanopyGasCO2_pft(NZ)

  IF(SineSolarAngle.GT.0.0.AND.CanopyLeafArea_pft(NZ).GT.ZEROP(NZ))THEN
!
    call PhotoActivePFT(NZ)
  ELSE
!
    MinCanPStomaResistH2O_pft(NZ)=MaxCanPStomaResistH2O_pft(NZ)
  ENDIF

  RETURN
  end associate
  END subroutine stomates
!------------------------------------------------------------------------------------------

  subroutine C3ShadedLeaves(K,N,M,L,NB,NZ,CH2O)
  implicit none
  integer, intent(in) :: K,N,M,L,NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8) :: ETLF,EGRO,PARX,PARJ
  real(r8) :: VL
!     begin_execution
  associate(                        &
    PARDiffus_zsec => plt_rad%PARDiffus_zsec  , &
    PARDirect_zsec   => plt_rad%PARDirect_zsec    , &
    CO2lmtRubiscoCarboxyRate_node  => plt_photo%CO2lmtRubiscoCarboxyRate_node , &
    LeafAUnshaded_seclyrnodbrpft  => plt_photo%LeafAUnshaded_seclyrnodbrpft , &
    RubiscoActivity_brpft   => plt_photo%RubiscoActivity_brpft  , &
    RubiscoCarboxyEff_node  => plt_photo%RubiscoCarboxyEff_node , &
    LigthSatCarboxyRate_node  => plt_photo%LigthSatCarboxyRate_node , &
    TAU0   => plt_rad%TAU0      &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PARDiffus_zsec=diffuse PAR flux
!     ETGR=light saturated e- transport rate
!     ETLF=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO=light-limited rubisco carboxylation rate
!
  PARX=QNTM*PARDiffus_zsec(N,M,L,NZ)
  PARJ=PARX+LigthSatCarboxyRate_node(K,NB,NZ)
  ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatCarboxyRate_node(K,NB,NZ)))/CURV2
  EGRO=ETLF*RubiscoCarboxyEff_node(K,NB,NZ)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     CH2O=total rubisco carboxylation rate
!     RubiscoActivity_brpft=N,P feedback inhibition on C3 CO2 fixation
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
  VL=AMIN1(CO2lmtRubiscoCarboxyRate_node(K,NB,NZ),EGRO)*RubiscoActivity_brpft(NB,NZ)
  CH2O=CH2O+VL*LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ)*TAU0(L+1)
  end associate
  end subroutine C3ShadedLeaves
!------------------------------------------------------------------------------------------

  subroutine C3SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
  implicit none
  integer, intent(in) :: K,N,M,L,NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8) :: ETLF,EGRO,PARX,PARJ
  real(r8) :: VL
!     begin_execution
  associate(                            &
    PARDirect_zsec    => plt_rad%PARDirect_zsec       , &
    CO2lmtRubiscoCarboxyRate_node   => plt_photo%CO2lmtRubiscoCarboxyRate_node    , &
    RubiscoCarboxyEff_node   => plt_photo%RubiscoCarboxyEff_node    , &
    LigthSatCarboxyRate_node   => plt_photo%LigthSatCarboxyRate_node    , &
    LeafAUnshaded_seclyrnodbrpft   => plt_photo%LeafAUnshaded_seclyrnodbrpft    , &
    RubiscoActivity_brpft    => plt_photo%RubiscoActivity_brpft     , &
    TAUS    =>  plt_rad%TAUS        &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PAR=direct PAR flux
!     LigthSatCarboxyRate_node=light saturated e- transport rate
!     ETLF=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO=light-limited rubisco carboxylation rate
!
  PARX=QNTM*PARDirect_zsec(N,M,L,NZ)
  PARJ=PARX+LigthSatCarboxyRate_node(K,NB,NZ)
  ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatCarboxyRate_node(K,NB,NZ)))/CURV2
  EGRO=ETLF*RubiscoCarboxyEff_node(K,NB,NZ)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     RubiscoActivity_brpft=N,P feedback inhibition on C3 CO2 fixation
!     CH2O=total rubisco carboxylation rate
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
  VL=AMIN1(CO2lmtRubiscoCarboxyRate_node(K,NB,NZ),EGRO)*RubiscoActivity_brpft(NB,NZ)
  CH2O=CH2O+VL*LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ)*TAUS(L+1)
  end associate
  end subroutine C3SunlitLeaves
!------------------------------------------------------------------------------------------

  subroutine C3PhotosynsCanopyLayerL(L,K,NB,NZ,CH2O)
  implicit none
  integer, intent(in) :: L,K,NB,NZ
  real(r8), intent(inout) :: CH2O
  integer :: N,M
!     begin_execution
  associate(                           &
    ZEROP  => plt_biom%ZEROP     , &
    LeafAUnshaded_seclyrnodbrpft  => plt_photo%LeafAUnshaded_seclyrnodbrpft    , &
    PARDirect_zsec   => plt_rad%PARDirect_zsec       , &
    PARDiffus_zsec => plt_rad%PARDiffus_zsec       &
  )
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO N=1,NumOfLeafZenithSectors1
    DO M=1,NumOfSkyAzimuthSectors1
      IF(LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ).GT.ZEROP(NZ))THEN
!
!     SUNLIT LEAVES
!
        IF(PARDirect_zsec(N,M,L,NZ).GT.0.0_r8)THEN
          call C3SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
!
!     SHADED LEAVES
        IF(PARDiffus_zsec(N,M,L,NZ).GT.0.0_r8)THEN
          call C3ShadedLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
!
      ENDIF
    ENDDO
  ENDDO
  end associate
  end subroutine C3PhotosynsCanopyLayerL
!------------------------------------------------------------------------------------------

  subroutine C3Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8), intent(in) :: TFN1,TFN2,TFNE,WSDN,XKO2L
  integer :: L
  real(r8) :: ETDN,VCDN
  real(r8) :: VOGRO
!     begin_execution
  associate(                          &
    CanopyLeafAreaByLayer_pft   =>  plt_morph%CanopyLeafAreaByLayer_pft , &
    ZEROP   =>  plt_biom%ZEROP  , &
    O2L     =>  plt_photo%O2L   , &
    aquCO2Intraleaf_pft    =>  plt_photo%aquCO2Intraleaf_pft  , &
    CHL     =>  plt_photo%CHL   , &
    Vmax4RubiscoCarboxy_pft   =>  plt_photo%Vmax4RubiscoCarboxy_pft , &
    Km4LeafaqCO2_pft  =>  plt_photo%Km4LeafaqCO2_pft, &
    Km4RubiscoCarboxy_pft  =>  plt_photo%Km4RubiscoCarboxy_pft, &
    ETMX    =>  plt_photo%ETMX  , &
    LigthSatCarboxyRate_node   =>  plt_photo%LigthSatCarboxyRate_node , &
    CO2lmtRubiscoCarboxyRate_node   =>  plt_photo%CO2lmtRubiscoCarboxyRate_node , &
    VCMX    =>  plt_photo%VCMX  , &
    RubiscoCarboxyEff_node   =>  plt_photo%RubiscoCarboxyEff_node , &
    VOMX    =>  plt_photo%VOMX  , &
    CO2CompenPoint_node   =>  plt_photo%CO2CompenPoint_node , &
    RUBP    =>  plt_photo%RUBP    &
  )
!
!     SURFICIAL DENSITY OF RUBISCO AND ITS CHLOROPHYLL
!
!     VCDN=surficial density of rubisco in mesophyll
!     ETDN=surficial density of chlorophyll in esophyll
!     RUBP=fraction of leaf protein in rubisco
!     CHL=fraction of leaf protein in mesophyll chlorophyll
!     WSDN=leaf protein surficial density
!
  VCDN=RUBP(NZ)*WSDN
  ETDN=CHL(NZ)*WSDN
!
!     CO2-LIMITED C3 CARBOXYLATION RATES
!
!     Vmax4RubiscoCarboxy_pft=rubisco carboxylation rate unlimited by CO2
!     VCMX=specific rubisco carboxylation activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN=surficial density of rubisco in mesophyll
!     VOGRO=rubisco oxygenation rate
!     TFN2=temperature function for oxygenation
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     aquCO2Intraleaf_pft,O2L=intercellular CO2,O2 concentrations (uM)
!     Km4LeafaqCO2_pft,Km4RubiscoCarboxy_pft=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!     CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2
!
  Vmax4RubiscoCarboxy_pft(K,NB,NZ)=VCMX(NZ)*TFN1*VCDN
  VOGRO=VOMX(NZ)*TFN2*VCDN
  CO2CompenPoint_node(K,NB,NZ)=0.5_r8*O2L(NZ)*VOGRO*Km4LeafaqCO2_pft(NZ)/(Vmax4RubiscoCarboxy_pft(K,NB,NZ)*XKO2L)
  CO2lmtRubiscoCarboxyRate_node(K,NB,NZ)=AZMAX1(Vmax4RubiscoCarboxy_pft(K,NB,NZ)*(aquCO2Intraleaf_pft(NZ)-CO2CompenPoint_node(K,NB,NZ)) &
    /(aquCO2Intraleaf_pft(NZ)+Km4RubiscoCarboxy_pft(NZ)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     LigthSatCarboxyRate_node=light-limited rubisco carboxylation rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN=surficial density of chlorophyll in mesophyll
!     RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
  LigthSatCarboxyRate_node(K,NB,NZ)=ETMX(NZ)*TFNE*ETDN
  RubiscoCarboxyEff_node(K,NB,NZ)=AZMAX1((aquCO2Intraleaf_pft(NZ)-CO2CompenPoint_node(K,NB,NZ))/(ELEC3*aquCO2Intraleaf_pft(NZ) &
    +10.5_r8*CO2CompenPoint_node(K,NB,NZ)))
!
!     FOR EACH CANOPY LAYER
!
!     CanopyLeafAreaByLayer_pft=leaf area
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!
  DO L=NumOfCanopyLayers1,1,-1
    IF(CanopyLeafAreaByLayer_pft(L,K,NB,NZ).GT.ZEROP(NZ))THEN
      call C3PhotosynsCanopyLayerL(L,K,NB,NZ,CH2O)
    ENDIF
  ENDDO
  end associate
  end subroutine C3Photosynthesis
!------------------------------------------------------------------------------------------

  subroutine C4Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8), intent(in) :: TFN1,TFN2,TFNE,XKO2L,WSDN
  integer :: L
  real(r8) :: CC4M
  real(r8) :: CCBS,ETDN4
  real(r8) :: ETDN
  real(r8) :: VCDN4,VCDN
  real(r8) :: VOGRO
!     begin_execution
  associate(                      &
    LeafChemElmntNode_brch   =>  plt_biom%LeafChemElmntNode_brch  , &
    ZEROP   =>  plt_biom%ZEROP  , &
    CanopyLeafAreaByLayer_pft   =>  plt_morph%CanopyLeafAreaByLayer_pft , &
    C4PhotosynDowreg_brch   =>  plt_photo%C4PhotosynDowreg_brch , &
    CHL4    =>  plt_photo%CHL4  , &
    O2L     =>  plt_photo%O2L   , &
    aquCO2Intraleaf_pft    =>  plt_photo%aquCO2Intraleaf_pft  , &
    ETMX    =>  plt_photo%ETMX  , &
    VCMX4   =>  plt_photo%VCMX4 , &
    RUBP    =>  plt_photo%RUBP  , &
    Km4PEPCarboxy_pft  =>  plt_photo%Km4PEPCarboxy_pft, &
    CMassCO2BundleSheath_node    =>  plt_photo%CMassCO2BundleSheath_node  , &
    CPOOL4  =>  plt_photo%CPOOL4, &
    Km4LeafaqCO2_pft  =>  plt_photo%Km4LeafaqCO2_pft, &
    NutrientCtrlonC4Carboxy_node   =>  plt_photo%NutrientCtrlonC4Carboxy_node , &
    Km4RubiscoCarboxy_pft  =>  plt_photo%Km4RubiscoCarboxy_pft, &
    CHL     =>  plt_photo%CHL   , &
    RubiscoCarboxyEff_node   =>  plt_photo%RubiscoCarboxyEff_node , &
    VCMX    =>  plt_photo%VCMX  , &
    LigthSatCarboxyRate_node   =>  plt_photo%LigthSatCarboxyRate_node , &
    CO2lmtRubiscoCarboxyRate_node   =>  plt_photo%CO2lmtRubiscoCarboxyRate_node , &
    Vmax4RubiscoCarboxy_pft   =>  plt_photo%Vmax4RubiscoCarboxy_pft , &
    CO2CompenPoint_node   =>  plt_photo%CO2CompenPoint_node , &
    Vmax4PEPCarboxy_pft   =>  plt_photo%Vmax4PEPCarboxy_pft , &
    CO2lmtPEPCarboxyRate_node  =>  plt_photo%CO2lmtPEPCarboxyRate_node, &
    VOMX    =>  plt_photo%VOMX  , &
    C4CarboxyEff_node   =>  plt_photo%C4CarboxyEff_node , &
    LigthSatC4CarboxyRate_node   =>  plt_photo%LigthSatC4CarboxyRate_node , &
    PEPC    =>  plt_photo%PEPC    &
  )
!
!     FEEDBACK ON C4 CARBOXYLATION FROM C4 NON-STRUCTURAL C
!
!     CC4M,CCBS=C4 nonstruct C concn in mesophyll,bundle sheath (uM)
!     CPOOL4,CMassCO2BundleSheath_node=C4 nonstructural C mass in mesophyll,bundle sheath
!     WGLF=leaf C mass
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!     NutrientCtrlonC4Carboxy_node=N,P feedback inhibition on C4 CO2 fixation
!
  CC4M=AZMAX1(0.021E+09_r8*CPOOL4(K,NB,NZ)/(LeafChemElmntNode_brch(ielmc,K,NB,NZ)*FMP))
  CCBS=AZMAX1(0.083E+09_r8*CMassCO2BundleSheath_node(K,NB,NZ)/(LeafChemElmntNode_brch(ielmc,K,NB,NZ)*FBS))
  NutrientCtrlonC4Carboxy_node(K,NB,NZ)=1.0_r8/(1.0_r8+CC4M/C4KI)
  NutrientCtrlonC4Carboxy_node(K,NB,NZ)=NutrientCtrlonC4Carboxy_node(K,NB,NZ)*C4PhotosynDowreg_brch(NB,NZ)
!
!     SURFICIAL DENSITY OF PEPC AND ITS CHLOROPHYLL
!
!     VCDN4=surficial density of PEP carboxylase in mesophyll
!     ETDN4=surficial density of chlorophyll in mesophyll
!     PEPC=fraction of leaf protein in PEP carboxylase
!     CHL4=fraction of leaf protein in mesophyll chlorophyll
!     WSDN=leaf protein surficial density
!
  VCDN4=PEPC(NZ)*WSDN
  ETDN4=CHL4(NZ)*WSDN
!
!     CO2-LIMITED C4 CARBOXYLATION RATES
!
!     Vmax4PEPCarboxy_pft,CO2lmtPEPCarboxyRate_node=PEP carboxylation rate unlimited,limited by CO2
!     VCMX4=specific PEP carboxylase activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN4=surficial density of PEP carboxylase in mesophyll
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     Km4PEPCarboxy_pft=Km for VCMX4 from PFT file (uM)
!
  Vmax4PEPCarboxy_pft(K,NB,NZ)=VCMX4(NZ)*TFN1*VCDN4
  CO2lmtPEPCarboxyRate_node(K,NB,NZ)=AZMAX1(Vmax4PEPCarboxy_pft(K,NB,NZ)*(aquCO2Intraleaf_pft(NZ)-COMP4)/(aquCO2Intraleaf_pft(NZ)+Km4PEPCarboxy_pft(NZ)))
!
!     C4 ELECTRON TRANSFER RATES
!
!     LigthSatC4CarboxyRate_node=light saturated e- transport rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN4=surficial density of chlorophyll in mesophyll
!     C4CarboxyEff_node=PEP caboxylation efficiency
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!
  LigthSatC4CarboxyRate_node(K,NB,NZ)=ETMX(NZ)*TFNE*ETDN4
  C4CarboxyEff_node(K,NB,NZ)=AZMAX1((aquCO2Intraleaf_pft(NZ)-COMP4)/(ELEC4*aquCO2Intraleaf_pft(NZ)+10.5_r8*COMP4))
!
!     FOR EACH CANOPY LAYER
!
!     CanopyLeafAreaByLayer_pft=leaf area
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!
  DO L=NumOfCanopyLayers1,1,-1
    IF(CanopyLeafAreaByLayer_pft(L,K,NB,NZ).GT.ZEROP(NZ))THEN
      call C4PhotosynsCanopyLayerL(L,K,NB,NZ,CH2O)
    ENDIF
  ENDDO
!
!     VARIABLES FOR C3 PHOTOSYNTHESIS DRIVEN BY C4
!
!     VCDN=surficial density of rubisco in bundle sheath
!     ETDN=surficial density of chlorophyll in bundle sheath
!     RUBP=fraction of leaf protein in rubisco
!     CHL=fraction of leaf protein in bundle sheath chlorophyll
!     WSDN=leaf protein surficial density
!
  VCDN=RUBP(NZ)*WSDN
  ETDN=CHL(NZ)*WSDN
!
!     CO2-LIMITED C3 CARBOXYLATION RATES
!
!     Vmax4RubiscoCarboxy_pft=rubisco carboxylation rate unlimited by CO2
!     VCMX=specific rubisco carboxylation activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN=surficial density of rubisco in bundle sheath
!     VOGRO=rubisco oxygenation rate
!     TFN2=temperature function for oxygenation
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     aquCO2Intraleaf_pft,O2L=intercellular CO2,O2 concentrations (uM)
!     Km4LeafaqCO2_pft,Km4RubiscoCarboxy_pft=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!     CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2
!     CCBS=C4 nonstruct C concn in bundle sheath (uM)
!
  Vmax4RubiscoCarboxy_pft(K,NB,NZ)=VCMX(NZ)*TFN1*VCDN
  VOGRO=VOMX(NZ)*TFN2*VCDN
  CO2CompenPoint_node(K,NB,NZ)=0.5_r8*O2L(NZ)*VOGRO*Km4LeafaqCO2_pft(NZ)/(Vmax4RubiscoCarboxy_pft(K,NB,NZ)*XKO2L)
  CO2lmtRubiscoCarboxyRate_node(K,NB,NZ)=AZMAX1(Vmax4RubiscoCarboxy_pft(K,NB,NZ)*(CCBS-CO2CompenPoint_node(K,NB,NZ))/(CCBS+Km4RubiscoCarboxy_pft(NZ)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     LigthSatCarboxyRate_node=light-limited rubisco carboxylation rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN=surficial density of chlorophyll in bundle sheath
!     RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
  LigthSatCarboxyRate_node(K,NB,NZ)=ETMX(NZ)*TFNE*ETDN
  RubiscoCarboxyEff_node(K,NB,NZ)=AZMAX1((CCBS-CO2CompenPoint_node(K,NB,NZ))/(ELEC3*CCBS+10.5*CO2CompenPoint_node(K,NB,NZ)))
  end associate
  end subroutine C4Photosynthesis
!------------------------------------------------------------------------------------------

  subroutine C4PhotosynsCanopyLayerL(L,K,NB,NZ,CH2O)
  implicit none
  integer, intent(in) :: L,K,NB,NZ
  real(r8), intent(inout) :: CH2O
  integer :: M,N
!     begin_execution
  associate(                        &
    LeafAUnshaded_seclyrnodbrpft  => plt_photo%LeafAUnshaded_seclyrnodbrpft , &
    ZEROP  =>  plt_biom%ZEROP , &
    PARDirect_zsec   => plt_rad%PARDirect_zsec    , &
    PARDiffus_zsec => plt_rad%PARDiffus_zsec    &
  )
!
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO  N=1,NumOfLeafZenithSectors1
    DO  M=1,NumOfSkyAzimuthSectors1
      IF(LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ).GT.ZEROP(NZ))THEN
!
!     SUNLIT LEAVES
        IF(PARDirect_zsec(N,M,L,NZ).GT.0.0_r8)THEN
          call C4SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
!
!     SHADED LEAVES
        IF(PARDiffus_zsec(N,M,L,NZ).GT.0.0_r8)THEN
          call C4ShadedLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  end associate
  end subroutine C4PhotosynsCanopyLayerL
!------------------------------------------------------------------------------------------

  subroutine C4ShadedLeaves(K,N,M,L,NB,NZ,CH2O)

  implicit none
  integer, intent(in) :: K,N,M,L,NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8) :: ETLF4,EGRO4,PARX,PARJ
  real(r8) :: VL
!     begin_execution
  associate(                          &
    LeafAUnshaded_seclyrnodbrpft  => plt_photo%LeafAUnshaded_seclyrnodbrpft   , &
    NutrientCtrlonC4Carboxy_node  => plt_photo%NutrientCtrlonC4Carboxy_node   , &
    C4CarboxyEff_node  => plt_photo%C4CarboxyEff_node   , &
    CO2lmtPEPCarboxyRate_node => plt_photo%CO2lmtPEPCarboxyRate_node  , &
    LigthSatC4CarboxyRate_node  => plt_photo%LigthSatC4CarboxyRate_node   , &
    PARDirect_zsec   => plt_rad%PARDirect_zsec      , &
    PARDiffus_zsec => plt_rad%PARDiffus_zsec    , &
    TAU0   => plt_rad%TAU0        &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PARDiffus_zsec=diffuse PAR flux
!     LigthSatC4CarboxyRate_node=light saturated e- transport rate
!     ETLF4=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO4=light-limited PEP carboxylation rate
!
  PARX=QNTM*PARDiffus_zsec(N,M,L,NZ)
  PARJ=PARX+LigthSatC4CarboxyRate_node(K,NB,NZ)
  ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatC4CarboxyRate_node(K,NB,NZ)))/CURV2
  EGRO4=ETLF4*C4CarboxyEff_node(K,NB,NZ)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     CO2lmtPEPCarboxyRate_node=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     CH2O=total PEP carboxylation rate
!     NutrientCtrlonC4Carboxy_node=N,P feedback inhibition on C4 CO2 fixation
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
  VL=AMIN1(CO2lmtPEPCarboxyRate_node(K,NB,NZ),EGRO4)*NutrientCtrlonC4Carboxy_node(K,NB,NZ)
  CH2O=CH2O+VL*LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ)*TAU0(L+1)
  end associate
  end subroutine C4ShadedLeaves
!------------------------------------------------------------------------------------------

  subroutine C4SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
  implicit none
  integer, intent(in) :: K,N,M,L,NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8) :: ETLF4,EGRO4,PARX,PARJ
  real(r8) :: VL
!     begin_execution
  associate(                       &
  LeafAUnshaded_seclyrnodbrpft  =>  plt_photo%LeafAUnshaded_seclyrnodbrpft , &
  C4CarboxyEff_node  =>  plt_photo%C4CarboxyEff_node , &
  LigthSatC4CarboxyRate_node  =>  plt_photo%LigthSatC4CarboxyRate_node , &
  NutrientCtrlonC4Carboxy_node  =>  plt_photo%NutrientCtrlonC4Carboxy_node , &
  CO2lmtPEPCarboxyRate_node =>  plt_photo%CO2lmtPEPCarboxyRate_node, &
  PARDirect_zsec   =>  plt_rad%PARDirect_zsec    , &
  TAUS   =>  plt_rad%TAUS      &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PAR=direct PAR flux
!     LigthSatC4CarboxyRate_node=light saturated e- transport rate
!     ETLF4=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO4=light-limited PEP carboxylation rate
!
  PARX=QNTM*PARDirect_zsec(N,M,L,NZ)
  PARJ=PARX+LigthSatC4CarboxyRate_node(K,NB,NZ)
  ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatC4CarboxyRate_node(K,NB,NZ)))/CURV2
  EGRO4=ETLF4*C4CarboxyEff_node(K,NB,NZ)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     CO2lmtPEPCarboxyRate_node=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     NutrientCtrlonC4Carboxy_node=N,P feedback inhibition on C4 CO2 fixation
!     CH2O=total PEP carboxylation rate
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
  VL=AMIN1(CO2lmtPEPCarboxyRate_node(K,NB,NZ),EGRO4)*NutrientCtrlonC4Carboxy_node(K,NB,NZ)
  CH2O=CH2O+VL*LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ)*TAUS(L+1)
  end associate
  end subroutine C4SunlitLeaves
!------------------------------------------------------------------------------------------

  subroutine LivingBranch(NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
  implicit none
  integer, intent(in):: NB,NZ
  real(r8),intent(inout) :: CH2O
  real(r8), intent(in) :: TFN1,TFN2,TFNE,XKO2L
  integer :: K
  real(r8) :: WSDN
!     begin_execution
  associate(                            &
    iPlantPhotosynthesisType   => plt_photo%iPlantPhotosynthesisType  , &
    Vmax4PEPCarboxy_pft    => plt_photo%Vmax4PEPCarboxy_pft   , &
    Vmax4RubiscoCarboxy_pft    => plt_photo%Vmax4RubiscoCarboxy_pft   , &
    ZERO     => plt_site%ZERO     , &
    LeafChemElmntNode_brch    => plt_biom%LeafChemElmntNode_brch    , &
    ZEROP    => plt_biom%ZEROP    , &
    LeafProteinCNode_brch     => plt_biom%LeafProteinCNode_brch     , &
    LeafAreaNode_brch    => plt_morph%LeafAreaNode_brch     &
  )
  DO K=1,MaxNodesPerBranch1
    IF(LeafAreaNode_brch(K,NB,NZ).GT.ZEROP(NZ).AND.LeafChemElmntNode_brch(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
      WSDN=LeafProteinCNode_brch(K,NB,NZ)/LeafAreaNode_brch(K,NB,NZ)
    ELSE
      WSDN=0.0_r8
    ENDIF

    IF(WSDN.GT.ZERO)THEN
!
!     iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4 from PFT file
!
      IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN
!     C4 PHOTOSYNTHESIS
        call C4Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
      ELSE IF(iPlantPhotosynthesisType(NZ).EQ.ic3_photo)then
!     C3 PHOTOSYNTHESIS
        call C3Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
      ENDIF
!
    ELSE
      Vmax4PEPCarboxy_pft(K,NB,NZ)=0.0_r8
      Vmax4RubiscoCarboxy_pft(K,NB,NZ)=0.0_r8
    ENDIF
  ENDDO
  end associate
  end subroutine LivingBranch
!------------------------------------------------------------------------------------------

  subroutine PhenoActiveBranch(NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8), intent(in) :: TFN1,TFN2,TFNE,XKO2L
!     begin_execution
  associate(                          &
    LeafPetioNonstructElmntConc_brch =>  plt_biom%LeafPetioNonstructElmntConc_brch  , &
    iPlantPhenologyType_pft =>  plt_pheno%iPlantPhenologyType_pft , &
    iPlantTurnoverPattern_pft =>  plt_pheno%iPlantTurnoverPattern_pft , &
    HourWithNoGrainFill_brch   =>  plt_pheno%HourWithNoGrainFill_brch   , &
    HourCounter4LeafOut_brch   =>  plt_pheno%HourCounter4LeafOut_brch   , &
    iPlantPhenologyPattern_pft =>  plt_pheno%iPlantPhenologyPattern_pft , &
    iPlantBranchState_brch  =>  plt_pheno%iPlantBranchState_brch  , &
    ZERO   =>  plt_site%ZERO    , &
    C4PhotosynDowreg_brch  => plt_photo%C4PhotosynDowreg_brch   , &
    RubiscoActivity_brpft   => plt_photo%RubiscoActivity_brpft      &
  )
!
!     FEEDBACK ON C3 CARBOXYLATION FROM NON-STRUCTURAL C:N:P
!
!     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
!     RubiscoActivity_brpft=N,P feedback inhibition on C3 CO2 fixation
!     CNKI,CPKI=nonstructural N,P inhibition constant on rubisco
!
  IF(LeafPetioNonstructElmntConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    RubiscoActivity_brpft(NB,NZ)=AMIN1(LeafPetioNonstructElmntConc_brch(ielmn,NB,NZ)/(LeafPetioNonstructElmntConc_brch(ielmn,NB,NZ)+LeafPetioNonstructElmntConc_brch(ielmc,NB,NZ)/CNKI) &
      ,LeafPetioNonstructElmntConc_brch(ielmp,NB,NZ)/(LeafPetioNonstructElmntConc_brch(ielmp,NB,NZ)+LeafPetioNonstructElmntConc_brch(ielmc,NB,NZ)/CPKI))
  ELSE
    RubiscoActivity_brpft(NB,NZ)=1.0_r8
  ENDIF
!
!     CHILLING
!
!     CHILL=accumulated chilling hours used to limit CO2 fixn
!
!     RubiscoActivity_brpft(NB,NZ)=RubiscoActivity_brpft(NB,NZ)/(1.0_r8+0.25_r8*CHILL(NZ))
!
!     DEHARDENING OF EVERGREENS IN SPRING
!
!     ATRP=hours above threshold temperature for dehardening since leafout
!     ATRPZ=hours to full dehardening of conifers in spring
! deciduous
  IF(iPlantPhenologyType_pft(NZ).NE.iphenotyp_evgreen.AND.iPlantTurnoverPattern_pft(NZ).GE.2)THEN
    RubiscoActivity_brpft(NB,NZ)=RubiscoActivity_brpft(NB,NZ)*AZMAX1(AMIN1(1.0_r8 &
      ,HourCounter4LeafOut_brch(NB,NZ)/(0.9_r8*ATRPZ)))
  ENDIF
!
!     TERMINATION OF ANNUALS
!
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     HourWithNoGrainFill_brch=number of hours with no grain fill after start of grain fill
!     Hours2KillAnuals=number of hours with no grain fill to terminate annuals
!
  IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.HourWithNoGrainFill_brch(NB,NZ).GT.0.0_r8)THEN
    C4PhotosynDowreg_brch(NB,NZ)=AZMAX1(1.0_r8-HourWithNoGrainFill_brch(NB,NZ)/Hours2KillAnuals(iPlantPhenologyType_pft(NZ)))
  ELSE
    C4PhotosynDowreg_brch(NB,NZ)=1.0_r8
  ENDIF
  RubiscoActivity_brpft(NB,NZ)=RubiscoActivity_brpft(NB,NZ)*C4PhotosynDowreg_brch(NB,NZ)
!
!     FOR EACH NODE
!
!     iPlantBranchState_brch=branch life flag:0=living,1=dead
!     LeafAreaNode_brch,WGLF,LeafProteinCNode_brch=leaf area,C mass,protein mass
!     WSDN=leaf protein surficial density
!
  IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
    call LivingBranch(NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
  ENDIF
  end associate
  end subroutine PhenoActiveBranch
!------------------------------------------------------------------------------------------

  subroutine PrepPhotosynthesis(NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(out) :: CH2O
  real(r8), intent(out) :: TFN1,TFN2,TFNE,XKO2L
  real(r8) :: ACTV,RTK
  real(r8) :: STK,TCCZ
  real(r8) :: TKCO
!     begin_execution
  associate(                          &
    LeafPetioNonstructElmntConc_brch =>  plt_biom%LeafPetioNonstructElmntConc_brch  , &
    TKCanopy_pft   =>  plt_ew%TKCanopy_pft      , &
    OFFST  =>  plt_pheno%OFFST  , &
    XKO2   =>  plt_photo%XKO2   , &
    CO2Solubility_pft  =>  plt_photo%CO2Solubility_pft  , &
    CanopyGasCO2_pft   =>  plt_photo%CanopyGasCO2_pft   , &
    O2L    =>  plt_photo%O2L    , &
    AirConc_pft   =>  plt_photo%AirConc_pft   , &
    DiffCO2Atmos2Intracel_pft  =>  plt_photo%DiffCO2Atmos2Intracel_pft  , &
    aquCO2Intraleaf_pft   =>  plt_photo%aquCO2Intraleaf_pft   , &
    LeafO2Solubility_pft    =>  plt_photo%LeafO2Solubility_pft    , &
    LeafIntracellularCO2_pft  =>  plt_photo%LeafIntracellularCO2_pft  , &
    O2I    =>  plt_photo%O2I    , &
    MaxCanPStomaResistH2O_pft   =>  plt_photo%MaxCanPStomaResistH2O_pft   , &
    Km4RubiscoCarboxy_pft =>  plt_photo%Km4RubiscoCarboxy_pft , &
    Km4LeafaqCO2_pft =>  plt_photo%Km4LeafaqCO2_pft , &
    XKCO2  =>  plt_photo%XKCO2    &

  )
!
!     CO2 AND O2 AQUEOUS SOLUBILITY
!
!     TCCZ=canopy temperature
!     SCO2,LeafO2Solubility_pft=solubility of CO2,O2 (uM/(umol mol-1))
!     aquCO2Intraleaf_pft,O2L=intercellular CO2,O2 concentrations (uM)
!     DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!
  TCCZ=TKCanopy_pft(NZ)-TC2K
  CO2Solubility_pft(NZ)=EXP(-2.621_r8-0.0317_r8*TCCZ)
  LeafO2Solubility_pft(NZ)=EXP(-6.175_r8-0.0211_r8*TCCZ)
  aquCO2Intraleaf_pft(NZ)=LeafIntracellularCO2_pft(NZ)*CO2Solubility_pft(NZ)
  O2L(NZ)=O2I(NZ)*LeafO2Solubility_pft(NZ)
  DiffCO2Atmos2Intracel_pft(NZ)=AirConc_pft(NZ)*(CanopyGasCO2_pft(NZ)-LeafIntracellularCO2_pft(NZ))
!
!     ARRHENIUS FUNCTIONS FOR CARBOXYLATION AND OXYGENATION
!
!     TKCanopy_pft,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
!     OFFST=shift in Arrhenius curve for thermal adaptation
!     TFN1,TFN2,TFNE=temperature function for carboxylation,
!     oxygenation,e- transport (25 oC =1)
!     8.313,710.0=gas constant,enthalpy
!     197500,222500=energy of high,low temp inactivn(KJ mol-1)
!     65000,60000,43000=activation energy for carboxylation,
!     oxygenation,e- transport
!
  CH2O=0.0_r8
  TKCO=TKCanopy_pft(NZ)+OFFST(NZ)
  RTK=RGAS*TKCO
  STK=710.0_r8*TKCO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  TFN1=EXP(26.237_r8-65000._r8/RTK)/ACTV
  TFN2=EXP(24.220_r8-60000._r8/RTK)/ACTV
  TFNE=EXP(17.362_r8-43000._r8/RTK)/ACTV
!
!     M-M CONSTANT FOR CARBOXYLATION FROM 'READQ' ADJUSTED FOR TEMPERATURE
!
!     Km4LeafaqCO2_pft,Km4RubiscoCarboxy_pft=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!
  Km4LeafaqCO2_pft(NZ)=XKCO2(NZ)*EXP(16.136_r8-40000._r8/RTK)
  XKO2L=XKO2(NZ)*EXP(8.067_r8-20000._r8/RTK)
  Km4RubiscoCarboxy_pft(NZ)=Km4LeafaqCO2_pft(NZ)*(1.0_r8+O2L(NZ)/XKO2L)
  end associate
  end subroutine PrepPhotosynthesis
!------------------------------------------------------------------------------------------

  subroutine PhotoActivePFT(NZ)
  implicit none
  integer, intent(in) :: NZ

  integer :: NB,K
  real(r8) :: CH2O
  real(r8) :: RSX
  real(r8) :: TFN1,TFN2,TFNE,XKO2L
  real(r8), parameter :: secsperhour=3600.0_r8

!     begin_execution
  associate(                          &
    HourThreshold4LeafOff_brch  =>  plt_pheno%HourThreshold4LeafOff_brch  , &
    Hours4LeafOff_brch   =>  plt_pheno%Hours4LeafOff_brch   , &
    HourThreshold4LeafOut_brch  =>  plt_pheno%HourThreshold4LeafOut_brch  , &
    Hours4Leafout_brch   =>  plt_pheno%Hours4Leafout_brch   , &
    iPlantPhenologyType_pft =>  plt_pheno%iPlantPhenologyType_pft , &
    ZEROP  =>  plt_biom%ZEROP   , &
    NU     =>  plt_site%NU      , &
    AREA3  =>  plt_site%AREA3   , &
    Vmax4RubiscoCarboxy_pft  =>  plt_photo%Vmax4RubiscoCarboxy_pft  , &
    RubiscoActivity_brpft   =>  plt_photo%RubiscoActivity_brpft   , &
    C4PhotosynDowreg_brch  =>  plt_photo%C4PhotosynDowreg_brch  , &
    Km4RubiscoCarboxy_pft =>  plt_photo%Km4RubiscoCarboxy_pft , &
    Vmax4PEPCarboxy_pft  =>  plt_photo%Vmax4PEPCarboxy_pft  , &
    MinCanPStomaResistH2O_pft   =>  plt_photo%MinCanPStomaResistH2O_pft   , &
    DiffCO2Atmos2Intracel_pft  =>  plt_photo%DiffCO2Atmos2Intracel_pft  , &
    MaxCanPStomaResistH2O_pft   =>  plt_photo%MaxCanPStomaResistH2O_pft   , &
    FracRadPARbyCanopy_pft  =>  plt_rad%FracRadPARbyCanopy_pft    , &
    NumOfBranches_pft    =>  plt_morph%NumOfBranches_pft      &
  )

  call PrepPhotosynthesis(NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
!
!     FOR EACH BRANCH
!
  DO NB=1,NumOfBranches_pft(NZ)
!
!     FEEDBACK ON CO2 FIXATION
!
!     iPlantPhenologyType_pft=phenology type from PFT file
!     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!
    IF(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_evgreen.OR.Hours4Leafout_brch(NB,NZ).GE.HourThreshold4LeafOut_brch(NB,NZ).OR.Hours4LeafOff_brch(NB,NZ).LT.HourThreshold4LeafOff_brch(NB,NZ))THEN

      call PhenoActiveBranch(NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
    ELSE
      RubiscoActivity_brpft(NB,NZ)=0.0_r8
      C4PhotosynDowreg_brch(NB,NZ)=1.0_r8
      DO K=1,MaxNodesPerBranch1
        Vmax4PEPCarboxy_pft(K,NB,NZ)=0.0_r8
        Vmax4RubiscoCarboxy_pft(K,NB,NZ)=0.0_r8
      ENDDO
    ENDIF
  ENDDO
!
!     MINIMUM CANOPY STOMATAL RESISTANCE FROM CO2 CONCENTRATION
!     DIFFERENCE DIVIDED BY TOTAL CO2 FIXATION
!
!     RSX,MinCanPStomaResistH2O_pft=minimum canopy stomatal resistance to CO2,H2O (h m-1)
!     CH2O=total PEP(C4) or rubisco(C3) carboxylation rate
!     FracRadPARbyCanopy_pft=fraction of radiation received by each PFT canopy
!     DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!     AREA=area of grid cell
!     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
! hourly time step
  IF(CH2O.GT.ZEROP(NZ))THEN
    RSX=FracRadPARbyCanopy_pft(NZ)*DiffCO2Atmos2Intracel_pft(NZ)*AREA3(NU)/(CH2O*secsperhour)
  ELSE
    RSX=MaxCanPStomaResistH2O_pft(NZ)*1.56_r8
  ENDIF
  MinCanPStomaResistH2O_pft(NZ)=AMIN1(MaxCanPStomaResistH2O_pft(NZ),AMAX1(RSMY,RSX*0.641_r8))
  end associate
  end subroutine PhotoActivePFT

end module stomatesMod
