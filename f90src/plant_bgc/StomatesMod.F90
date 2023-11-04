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
  REAL(R8):: CanPbndlResist4CO2
  real(r8):: RI
!     begin_execution
  associate(                            &
    RIB      =>  plt_ew%RIB       , &
    RAZ      =>  plt_ew%RAZ       , &
    TairK      =>  plt_ew%TairK       , &
    TKCZ     =>  plt_ew%TKCZ      , &
    CO2E     =>  plt_site%CO2E    , &
    CO2Q     =>  plt_photo%CO2Q   , &
    CanopyLeafA_pft    =>  plt_morph%CanopyLeafA_pft  , &
    ZEROP    =>  plt_biom%ZEROP   , &
    CNETX    =>  plt_bgcr%CNETX   , &
    SSIN     =>  plt_rad%SSIN     , &
    FMOL     =>  plt_photo%FMOL   , &
    MinCanPStomaResistH2O     =>  plt_photo%MinCanPStomaResistH2O   , &
    CanPCi2CaRatio     =>  plt_photo%CanPCi2CaRatio   , &
    MaxCanPStomaResistH2O     =>  plt_photo%MaxCanPStomaResistH2O   , &
    CO2I     =>  plt_photo%CO2I     &
  )
!
!     CANOPY TEMPERATURE + OFFSET FOR THERMAL ADAPTATION FROM 'READQ'
!
!     CANOPY BOUNDARY LAYER RESISTANCE
!
!     RI=Richardson's number
!     RIB=canopy isothermal Richardson's number
!     TairK,TKCZ=air,canopy temperature
!     RAZ=canopy isothermal boundary later resistance
!     CanPbndlResist4CO2=canopy boundary layer resistance to CO2, h/m
!     FMOL=number of moles of air per m3
!
  RI=RichardsonNumber(RIB,TairK,TKCZ(NZ))

  CanPbndlResist4CO2=1.34*AMAX1(5.56E-03,RAZ(NZ)/(1.0-10.0*RI))

  !assuming pressure is one atmosphere
  FMOL(NZ)=GetMolAirPerm3(TKCZ(NZ))
!
!     CANOPY CO2 CONCENTRATION FROM CO2 INFLUXES AND EFFLUXES
!
!     CO2Q,CO2E=CO2 concentrations in canopy air,atmosphere, umol mol-1 (ppmv)
!     CNETX=net CO2 flux in canopy air from soil,plants, g d-2 h-1
! assuming steady state, canopy CO2 concentration is computed with mass balance. 
! how 8.33E+04 is determined. 
  CO2Q(NZ)=CO2E-8.33E+04_r8*CNETX*CanPbndlResist4CO2/FMOL(NZ)
  CO2Q(NZ)=AMIN1(CO2E+200.0_r8,AZMAX1(CO2E-200.0_r8,CO2Q(NZ)))
!
!     MESOPHYLL CO2 CONCENTRATION FROM CI:CA RATIO ENTERED IN 'READQ'
!
!     CO2I=intercellular CO2 concentration
!     CanPCi2CaRatio=intercellular:atmospheric CO2 concn ratio from PFT file, parameter
!     SSIN=sine of solar angle
!     CanopyLeafA_pft=PFT leaf area
!
  CO2I(NZ)=CanPCi2CaRatio(NZ)*CO2Q(NZ)

  IF(SSIN.GT.0.0.AND.CanopyLeafA_pft(NZ).GT.ZEROP(NZ))THEN
!
    call PhotoActivePFT(NZ)
  ELSE
!
    MinCanPStomaResistH2O(NZ)=MaxCanPStomaResistH2O(NZ)
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
    PARDIF => plt_rad%PARDIF  , &
    PAR    => plt_rad%PAR     , &
    VGRO   => plt_photo%VGRO  , &
    LeafAUnshaded_seclyrnodbrpft  => plt_photo%LeafAUnshaded_seclyrnodbrpft , &
    RubiscoActivity_brpft   => plt_photo%RubiscoActivity_brpft  , &
    CBXN   => plt_photo%CBXN  , &
    ETGRO  => plt_photo%ETGRO , &
    TAU0   => plt_rad%TAU0      &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PARDIF=diffuse PAR flux
!     ETGR=light saturated e- transport rate
!     ETLF=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO=light-limited rubisco carboxylation rate
!
  PARX=QNTM*PARDIF(N,M,L,NZ)
  PARJ=PARX+ETGRO(K,NB,NZ)
  ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ)))/CURV2
  EGRO=ETLF*CBXN(K,NB,NZ)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     VGRO=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     CH2O=total rubisco carboxylation rate
!     RubiscoActivity_brpft=N,P feedback inhibition on C3 CO2 fixation
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
  VL=AMIN1(VGRO(K,NB,NZ),EGRO)*RubiscoActivity_brpft(NB,NZ)
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
    PAR     => plt_rad%PAR        , &
    VGRO    => plt_photo%VGRO     , &
    CBXN    => plt_photo%CBXN     , &
    ETGRO   => plt_photo%ETGRO    , &
    LeafAUnshaded_seclyrnodbrpft   => plt_photo%LeafAUnshaded_seclyrnodbrpft    , &
    RubiscoActivity_brpft    => plt_photo%RubiscoActivity_brpft     , &
    TAUS    =>  plt_rad%TAUS        &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PAR=direct PAR flux
!     ETGRO=light saturated e- transport rate
!     ETLF=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO=light-limited rubisco carboxylation rate
!
  PARX=QNTM*PAR(N,M,L,NZ)
  PARJ=PARX+ETGRO(K,NB,NZ)
  ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ)))/CURV2
  EGRO=ETLF*CBXN(K,NB,NZ)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     VGRO=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     RubiscoActivity_brpft=N,P feedback inhibition on C3 CO2 fixation
!     CH2O=total rubisco carboxylation rate
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
  VL=AMIN1(VGRO(K,NB,NZ),EGRO)*RubiscoActivity_brpft(NB,NZ)
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
    PAR    => plt_rad%PAR        , &
    PARDIF => plt_rad%PARDIF       &
  )
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO N=1,JLI1
    DO M=1,JSA1
      IF(LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ).GT.ZEROP(NZ))THEN
!
!     SUNLIT LEAVES
!

        IF(PAR(N,M,L,NZ).GT.0.0_r8)THEN
          call C3SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
!
!     SHADED LEAVES
        IF(PARDIF(N,M,L,NZ).GT.0.0_r8)THEN
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
    CanPLNBLA   =>  plt_morph%CanPLNBLA , &
    ZEROP   =>  plt_biom%ZEROP  , &
    O2L     =>  plt_photo%O2L   , &
    CO2L    =>  plt_photo%CO2L  , &
    CHL     =>  plt_photo%CHL   , &
    VCGRO   =>  plt_photo%VCGRO , &
    XKCO2L  =>  plt_photo%XKCO2L, &
    XKCO2O  =>  plt_photo%XKCO2O, &
    ETMX    =>  plt_photo%ETMX  , &
    ETGRO   =>  plt_photo%ETGRO , &
    VGRO    =>  plt_photo%VGRO  , &
    VCMX    =>  plt_photo%VCMX  , &
    CBXN    =>  plt_photo%CBXN  , &
    VOMX    =>  plt_photo%VOMX  , &
    COMPL   =>  plt_photo%COMPL , &
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
!     VCGRO=rubisco carboxylation rate unlimited by CO2
!     VCMX=specific rubisco carboxylation activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN=surficial density of rubisco in mesophyll
!     VOGRO=rubisco oxygenation rate
!     TFN2=temperature function for oxygenation
!     COMPL=C3 CO2 compensation point (uM)
!     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
!     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!     VGRO=rubisco carboxylation rate limited by CO2
!
  VCGRO(K,NB,NZ)=VCMX(NZ)*TFN1*VCDN
  VOGRO=VOMX(NZ)*TFN2*VCDN
  COMPL(K,NB,NZ)=0.5_r8*O2L(NZ)*VOGRO*XKCO2L(NZ)/(VCGRO(K,NB,NZ)*XKO2L)
  VGRO(K,NB,NZ)=AZMAX1(VCGRO(K,NB,NZ)*(CO2L(NZ)-COMPL(K,NB,NZ)) &
    /(CO2L(NZ)+XKCO2O(NZ)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     ETGRO=light-limited rubisco carboxylation rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN=surficial density of chlorophyll in mesophyll
!     CBXN=rubisco caboxylation efficiency
!     CO2L=intercellular CO2 concentrations (uM)
!     COMPL=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
  ETGRO(K,NB,NZ)=ETMX(NZ)*TFNE*ETDN
  CBXN(K,NB,NZ)=AZMAX1((CO2L(NZ)-COMPL(K,NB,NZ))/(ELEC3*CO2L(NZ) &
    +10.5_r8*COMPL(K,NB,NZ)))
!
!     FOR EACH CANOPY LAYER
!
!     CanPLNBLA=leaf area
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!
  DO L=NumOfCanopyLayers1,1,-1
    IF(CanPLNBLA(L,K,NB,NZ).GT.ZEROP(NZ))THEN
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
    WGLFE   =>  plt_biom%WGLFE  , &
    ZEROP   =>  plt_biom%ZEROP  , &
    CanPLNBLA   =>  plt_morph%CanPLNBLA , &
    FDBKX   =>  plt_photo%FDBKX , &
    CHL4    =>  plt_photo%CHL4  , &
    O2L     =>  plt_photo%O2L   , &
    CO2L    =>  plt_photo%CO2L  , &
    ETMX    =>  plt_photo%ETMX  , &
    VCMX4   =>  plt_photo%VCMX4 , &
    RUBP    =>  plt_photo%RUBP  , &
    XKCO24  =>  plt_photo%XKCO24, &
    CO2B    =>  plt_photo%CO2B  , &
    CPOOL4  =>  plt_photo%CPOOL4, &
    XKCO2L  =>  plt_photo%XKCO2L, &
    FDBK4   =>  plt_photo%FDBK4 , &
    XKCO2O  =>  plt_photo%XKCO2O, &
    CHL     =>  plt_photo%CHL   , &
    CBXN    =>  plt_photo%CBXN  , &
    VCMX    =>  plt_photo%VCMX  , &
    ETGRO   =>  plt_photo%ETGRO , &
    VGRO    =>  plt_photo%VGRO  , &
    VCGRO   =>  plt_photo%VCGRO , &
    COMPL   =>  plt_photo%COMPL , &
    VCGR4   =>  plt_photo%VCGR4 , &
    VGRO4   =>  plt_photo%VGRO4 , &
    VOMX    =>  plt_photo%VOMX  , &
    CBXN4   =>  plt_photo%CBXN4 , &
    ETGR4   =>  plt_photo%ETGR4 , &
    PEPC    =>  plt_photo%PEPC    &
  )
!
!     FEEDBACK ON C4 CARBOXYLATION FROM C4 NON-STRUCTURAL C
!
!     CC4M,CCBS=C4 nonstruct C concn in mesophyll,bundle sheath (uM)
!     CPOOL4,CO2B=C4 nonstructural C mass in mesophyll,bundle sheath
!     WGLF=leaf C mass
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
  CC4M=AZMAX1(0.021E+09_r8*CPOOL4(K,NB,NZ)/(WGLFE(ielmc,K,NB,NZ)*FMP))
  CCBS=AZMAX1(0.083E+09_r8*CO2B(K,NB,NZ)/(WGLFE(ielmc,K,NB,NZ)*FBS))
  FDBK4(K,NB,NZ)=1.0_r8/(1.0_r8+CC4M/C4KI)
  FDBK4(K,NB,NZ)=FDBK4(K,NB,NZ)*FDBKX(NB,NZ)
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
!     VCGR4,VGRO4=PEP carboxylation rate unlimited,limited by CO2
!     VCMX4=specific PEP carboxylase activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN4=surficial density of PEP carboxylase in mesophyll
!     CO2L=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     XKCO24=Km for VCMX4 from PFT file (uM)
!
  VCGR4(K,NB,NZ)=VCMX4(NZ)*TFN1*VCDN4
  VGRO4(K,NB,NZ)=AZMAX1(VCGR4(K,NB,NZ)*(CO2L(NZ)-COMP4)/(CO2L(NZ)+XKCO24(NZ)))
!
!     C4 ELECTRON TRANSFER RATES
!
!     ETGR4=light saturated e- transport rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN4=surficial density of chlorophyll in mesophyll
!     CBXN4=PEP caboxylation efficiency
!     CO2L=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!
  ETGR4(K,NB,NZ)=ETMX(NZ)*TFNE*ETDN4
  CBXN4(K,NB,NZ)=AZMAX1((CO2L(NZ)-COMP4)/(ELEC4*CO2L(NZ)+10.5_r8*COMP4))
!
!     FOR EACH CANOPY LAYER
!
!     CanPLNBLA=leaf area
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!
  DO L=NumOfCanopyLayers1,1,-1
    IF(CanPLNBLA(L,K,NB,NZ).GT.ZEROP(NZ))THEN
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
!     VCGRO=rubisco carboxylation rate unlimited by CO2
!     VCMX=specific rubisco carboxylation activity from PFT file
!     TFN1=temperature function for carboxylation
!     VCDN=surficial density of rubisco in bundle sheath
!     VOGRO=rubisco oxygenation rate
!     TFN2=temperature function for oxygenation
!     COMPL=C3 CO2 compensation point (uM)
!     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
!     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!     VGRO=rubisco carboxylation rate limited by CO2
!     CCBS=C4 nonstruct C concn in bundle sheath (uM)
!
  VCGRO(K,NB,NZ)=VCMX(NZ)*TFN1*VCDN
  VOGRO=VOMX(NZ)*TFN2*VCDN
  COMPL(K,NB,NZ)=0.5_r8*O2L(NZ)*VOGRO*XKCO2L(NZ)/(VCGRO(K,NB,NZ)*XKO2L)
  VGRO(K,NB,NZ)=AZMAX1(VCGRO(K,NB,NZ)*(CCBS-COMPL(K,NB,NZ))/(CCBS+XKCO2O(NZ)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     ETGRO=light-limited rubisco carboxylation rate
!     ETMX=specific chlorophyll activity from PFT file
!     TFNE=temperature function for e- transport
!     ETDN=surficial density of chlorophyll in bundle sheath
!     CBXN=rubisco caboxylation efficiency
!     CO2L=intercellular CO2 concentrations (uM)
!     COMPL=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
  ETGRO(K,NB,NZ)=ETMX(NZ)*TFNE*ETDN
  CBXN(K,NB,NZ)=AZMAX1((CCBS-COMPL(K,NB,NZ))/(ELEC3*CCBS+10.5*COMPL(K,NB,NZ)))
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
    PAR    => plt_rad%PAR     , &
    PARDIF => plt_rad%PARDIF    &
  )
!
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO  N=1,JLI1
    DO  M=1,JSA1
      IF(LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ).GT.ZEROP(NZ))THEN
!
!     SUNLIT LEAVES
        IF(PAR(N,M,L,NZ).GT.0.0_r8)THEN
          call C4SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
!
!     SHADED LEAVES
        IF(PARDIF(N,M,L,NZ).GT.0.0_r8)THEN
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
    FDBK4  => plt_photo%FDBK4   , &
    CBXN4  => plt_photo%CBXN4   , &
    VGRO4  => plt_photo%VGRO4   , &
    ETGR4  => plt_photo%ETGR4   , &
    PAR    => plt_rad%PAR       , &
    PARDIF => plt_rad%PARDIF    , &
    TAU0   => plt_rad%TAU0        &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PARDIF=diffuse PAR flux
!     ETGR4=light saturated e- transport rate
!     ETLF4=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO4=light-limited PEP carboxylation rate
!
  PARX=QNTM*PARDIF(N,M,L,NZ)
  PARJ=PARX+ETGR4(K,NB,NZ)
  ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ)))/CURV2
  EGRO4=ETLF4*CBXN4(K,NB,NZ)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     VGRO4=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     CH2O=total PEP carboxylation rate
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
  VL=AMIN1(VGRO4(K,NB,NZ),EGRO4)*FDBK4(K,NB,NZ)
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
  CBXN4  =>  plt_photo%CBXN4 , &
  ETGR4  =>  plt_photo%ETGR4 , &
  FDBK4  =>  plt_photo%FDBK4 , &
  VGRO4  =>  plt_photo%VGRO4 , &
  PAR    =>  plt_rad%PAR     , &
  TAUS   =>  plt_rad%TAUS      &
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     PAR=direct PAR flux
!     ETGR4=light saturated e- transport rate
!     ETLF4=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO4=light-limited PEP carboxylation rate
!
  PARX=QNTM*PAR(N,M,L,NZ)
  PARJ=PARX+ETGR4(K,NB,NZ)
  ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ)))/CURV2
  EGRO4=ETLF4*CBXN4(K,NB,NZ)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     VGRO4=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!     CH2O=total PEP carboxylation rate
!     LeafAUnshaded_seclyrnodbrpft=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
  VL=AMIN1(VGRO4(K,NB,NZ),EGRO4)*FDBK4(K,NB,NZ)
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
    ICTYP    => plt_photo%ICTYP   , &
    VCGR4    => plt_photo%VCGR4   , &
    VCGRO    => plt_photo%VCGRO   , &
    ZERO     => plt_site%ZERO     , &
    WGLFE    => plt_biom%WGLFE    , &
    ZEROP    => plt_biom%ZEROP    , &
    WSLF     => plt_biom%WSLF     , &
    ARLF1    => plt_morph%ARLF1     &
  )
  DO K=1,MaxCanopyNodes1
    IF(ARLF1(K,NB,NZ).GT.ZEROP(NZ).AND.WGLFE(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
      WSDN=WSLF(K,NB,NZ)/ARLF1(K,NB,NZ)
    ELSE
      WSDN=0.0_r8
    ENDIF

    IF(WSDN.GT.ZERO)THEN
!
!     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!
      IF(ICTYP(NZ).EQ.ic4_photo)THEN
!     C4 PHOTOSYNTHESIS
        call C4Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
      ELSE
!     C3 PHOTOSYNTHESIS
        call C3Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
      ENDIF
!
    ELSE
      VCGR4(K,NB,NZ)=0.0_r8
      VCGRO(K,NB,NZ)=0.0_r8
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
    CEPOLB =>  plt_biom%CEPOLB  , &
    IWTYP  =>  plt_pheno%IWTYP  , &
    IBTYP  =>  plt_pheno%IBTYP  , &
    FLG4   =>  plt_pheno%FLG4   , &
    HourCounter4LeafOut_brch   =>  plt_pheno%HourCounter4LeafOut_brch   , &
    ISTYP  =>  plt_pheno%ISTYP  , &
    iPlantBranchState  =>  plt_pheno%iPlantBranchState  , &
    ZERO   =>  plt_site%ZERO    , &
    FDBKX  => plt_photo%FDBKX   , &
    RubiscoActivity_brpft   => plt_photo%RubiscoActivity_brpft      &
  )
!
!     FEEDBACK ON C3 CARBOXYLATION FROM NON-STRUCTURAL C:N:P
!
!     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
!     RubiscoActivity_brpft=N,P feedback inhibition on C3 CO2 fixation
!     CNKI,CPKI=nonstructural N,P inhibition constant on rubisco
!
  IF(CEPOLB(ielmc,NB,NZ).GT.ZERO)THEN
    RubiscoActivity_brpft(NB,NZ)=AMIN1(CEPOLB(ielmn,NB,NZ)/(CEPOLB(ielmn,NB,NZ)+CEPOLB(ielmc,NB,NZ)/CNKI) &
      ,CEPOLB(ielmp,NB,NZ)/(CEPOLB(ielmp,NB,NZ)+CEPOLB(ielmc,NB,NZ)/CPKI))
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
  IF(IWTYP(NZ).NE.0.AND.IBTYP(NZ).GE.2)THEN
    RubiscoActivity_brpft(NB,NZ)=RubiscoActivity_brpft(NB,NZ)*AZMAX1(AMIN1(1.0_r8 &
      ,HourCounter4LeafOut_brch(NB,NZ)/(0.9_r8*ATRPZ)))
  ENDIF
!
!     TERMINATION OF ANNUALS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     FLG4=number of hours with no grain fill after start of grain fill
!     Hours2KillAnuals=number of hours with no grain fill to terminate annuals
!
  IF(ISTYP(NZ).EQ.iplt_annual.AND.FLG4(NB,NZ).GT.0.0_r8)THEN
    FDBKX(NB,NZ)=AZMAX1(1.0_r8-FLG4(NB,NZ)/Hours2KillAnuals(IWTYP(NZ)))
  ELSE
    FDBKX(NB,NZ)=1.0_r8
  ENDIF
  RubiscoActivity_brpft(NB,NZ)=RubiscoActivity_brpft(NB,NZ)*FDBKX(NB,NZ)
!
!     FOR EACH NODE
!
!     iPlantBranchState=branch life flag:0=living,1=dead
!     ARLF,WGLF,WSLF=leaf area,C mass,protein mass
!     WSDN=leaf protein surficial density
!
  IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
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
    CEPOLB =>  plt_biom%CEPOLB  , &
    TKCZ   =>  plt_ew%TKCZ      , &
    OFFST  =>  plt_pheno%OFFST  , &
    XKO2   =>  plt_photo%XKO2   , &
    SCO2   =>  plt_photo%SCO2   , &
    CO2Q   =>  plt_photo%CO2Q   , &
    O2L    =>  plt_photo%O2L    , &
    FMOL   =>  plt_photo%FMOL   , &
    DCO2   =>  plt_photo%DCO2   , &
    CO2L   =>  plt_photo%CO2L   , &
    SO2    =>  plt_photo%SO2    , &
    CO2I   =>  plt_photo%CO2I   , &
    O2I    =>  plt_photo%O2I    , &
    MaxCanPStomaResistH2O   =>  plt_photo%MaxCanPStomaResistH2O   , &
    XKCO2O =>  plt_photo%XKCO2O , &
    XKCO2L =>  plt_photo%XKCO2L , &
    XKCO2  =>  plt_photo%XKCO2    &

  )
!
!     CO2 AND O2 AQUEOUS SOLUBILITY
!
!     TCCZ=canopy temperature
!     SCO2,SO2=solubility of CO2,O2 (uM/(umol mol-1))
!     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
!     DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!
  TCCZ=TKCZ(NZ)-TC2K
  SCO2(NZ)=EXP(-2.621_r8-0.0317_r8*TCCZ)
  SO2(NZ)=EXP(-6.175_r8-0.0211_r8*TCCZ)
  CO2L(NZ)=CO2I(NZ)*SCO2(NZ)
  O2L(NZ)=O2I(NZ)*SO2(NZ)
  DCO2(NZ)=FMOL(NZ)*(CO2Q(NZ)-CO2I(NZ))
!
!     ARRHENIUS FUNCTIONS FOR CARBOXYLATION AND OXYGENATION
!
!     TKCZ,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
!     OFFST=shift in Arrhenius curve for thermal adaptation
!     TFN1,TFN2,TFNE=temperature function for carboxylation,
!     oxygenation,e- transport (25 oC =1)
!     8.313,710.0=gas constant,enthalpy
!     197500,222500=energy of high,low temp inactivn(KJ mol-1)
!     65000,60000,43000=activation energy for carboxylation,
!     oxygenation,e- transport
!
  CH2O=0.0_r8
  TKCO=TKCZ(NZ)+OFFST(NZ)
  RTK=RGAS*TKCO
  STK=710.0_r8*TKCO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  TFN1=EXP(26.237_r8-65000._r8/RTK)/ACTV
  TFN2=EXP(24.220_r8-60000._r8/RTK)/ACTV
  TFNE=EXP(17.362_r8-43000._r8/RTK)/ACTV
!
!     M-M CONSTANT FOR CARBOXYLATION FROM 'READQ' ADJUSTED FOR TEMPERATURE
!
!     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
!     XKO2L=Km for rubisco oxygenation
!
  XKCO2L(NZ)=XKCO2(NZ)*EXP(16.136_r8-40000._r8/RTK)
  XKO2L=XKO2(NZ)*EXP(8.067_r8-20000._r8/RTK)
  XKCO2O(NZ)=XKCO2L(NZ)*(1.0_r8+O2L(NZ)/XKO2L)
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
    HourThreshold4LeafOff  =>  plt_pheno%HourThreshold4LeafOff  , &
    Hours4LeafOff   =>  plt_pheno%Hours4LeafOff   , &
    HourThreshold4LeafOut  =>  plt_pheno%HourThreshold4LeafOut  , &
    VRNS   =>  plt_pheno%VRNS   , &
    IWTYP  =>  plt_pheno%IWTYP  , &
    ZEROP  =>  plt_biom%ZEROP   , &
    NU     =>  plt_site%NU      , &
    AREA3  =>  plt_site%AREA3   , &
    VCGRO  =>  plt_photo%VCGRO  , &
    RubiscoActivity_brpft   =>  plt_photo%RubiscoActivity_brpft   , &
    FDBKX  =>  plt_photo%FDBKX  , &
    XKCO2O =>  plt_photo%XKCO2O , &
    VCGR4  =>  plt_photo%VCGR4  , &
    MinCanPStomaResistH2O   =>  plt_photo%MinCanPStomaResistH2O   , &
    DCO2   =>  plt_photo%DCO2   , &
    MaxCanPStomaResistH2O   =>  plt_photo%MaxCanPStomaResistH2O   , &
    FracPARByCanP  =>  plt_rad%FracPARByCanP    , &
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
!     IWTYP=phenology type from PFT file
!     VRNS,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!
    IF(IWTYP(NZ).EQ.0.OR.VRNS(NB,NZ).GE.HourThreshold4LeafOut(NB,NZ).OR.Hours4LeafOff(NB,NZ).LT.HourThreshold4LeafOff(NB,NZ))THEN

      call PhenoActiveBranch(NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
    ELSE
      RubiscoActivity_brpft(NB,NZ)=0.0_r8
      FDBKX(NB,NZ)=1.0_r8
      DO K=1,MaxCanopyNodes1
        VCGR4(K,NB,NZ)=0.0_r8
        VCGRO(K,NB,NZ)=0.0_r8
      ENDDO
    ENDIF
  ENDDO
!
!     MINIMUM CANOPY STOMATAL RESISTANCE FROM CO2 CONCENTRATION
!     DIFFERENCE DIVIDED BY TOTAL CO2 FIXATION
!
!     RSX,MinCanPStomaResistH2O=minimum canopy stomatal resistance to CO2,H2O (h m-1)
!     CH2O=total PEP(C4) or rubisco(C3) carboxylation rate
!     FracPARByCanP=fraction of radiation received by each PFT canopy
!     DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!     AREA=area of grid cell
!     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
! hourly time step
  IF(CH2O.GT.ZEROP(NZ))THEN
    RSX=FracPARByCanP(NZ)*DCO2(NZ)*AREA3(NU)/(CH2O*secsperhour)
  ELSE
    RSX=MaxCanPStomaResistH2O(NZ)*1.56_r8
  ENDIF
  MinCanPStomaResistH2O(NZ)=AMIN1(MaxCanPStomaResistH2O(NZ),AMAX1(RSMY,RSX*0.641_r8))
  end associate
  end subroutine PhotoActivePFT

end module stomatesMod
