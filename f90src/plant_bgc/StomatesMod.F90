  module stomatesMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use PlantAPIData
  implicit none

  private
  character(len=*), parameter :: mod_filename = __FILE__

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
  real(r8), parameter :: FLG4Y(0:5)=real((/336.0,672.0,672.0,672.0,672.0,672.0/),r8)  !number of hours with no grain fill to terminate annuals

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
  REAL(R8):: RACL
  real(r8):: RI
!     begin_execution
  associate(                            &
    RIBs1      =>  plt_ew%RIBs1       , &
    RAZs1      =>  plt_ew%RAZs1       , &
    TKAs1      =>  plt_ew%TKAs1       , &
    CO2Qs1     =>  plt_photo%CO2Qs1   , &
    ARLFPs1    =>  plt_morph%ARLFPs1  , &
    SSINs1     =>  plt_rad%SSINs1     , &
    FCO2s1     =>  plt_photo%FCO2s1   , &
    RSMHs1     =>  plt_photo%RSMHs1   , &
    CO2Is1     =>  plt_photo%CO2Is1     &
  )
!
!     CANOPY TEMPERATURE + OFFSET FOR THERMAL ADAPTATION FROM 'READQ'
!
!     CANOPY BOUNDARY LAYER RESISTANCE
!
!     RI=Richardson's number
!     RIB=canopy isothermal Richardson�s number
!     TKA,TKCZ=air,canopy temperature
!     RAZ=canopy isothermal boundary later resistance
!     RACL=canopy boundary layer resistance to CO2
!     FMOL=number of moles of air per m3
!
  RI=AMAX1(-0.3,AMIN1(0.075,RIBs1*(TKAs1-TKCZs1(NZ))))
  RACL=1.34*AMAX1(5.56E-03,RAZs1(NZ)/(1.0-10.0*RI))
  FMOLs1(NZ)=1.2194E+04/TKCZs1(NZ)
!
!     CANOPY CO2 CONCENTRATION FROM CO2 INFLUXES AND EFFLUXES
!
!     CO2Q,CO2E=CO2 concentrations in canopy air,atmosphere, umol mol-1 (ppmv)
!     CNETX=net CO2 flux in canopy air from soil,plants, g d-2 h-1
!
  CO2Qs1(NZ)=CO2Es1-8.33E+04_r8*CNETXs1*RACL/FMOLs1(NZ)
  CO2Qs1(NZ)=AMIN1(CO2Es1+200.0_r8,AMAX1(0.0_r8,CO2Es1-200.0_r8,CO2Qs1(NZ)))
!
!     MESOPHYLL CO2 CONCENTRATION FROM CI:CA RATIO ENTERED IN 'READQ'
!
!     CO2I=intercellular CO2 concentration
!     FCO2=intercellular:atmospheric CO2 concn ratio from PFT file, parameter
!     SSIN=sine of solar angle
!     ARLFP=PFT leaf area
!
  CO2Is1(NZ)=FCO2s1(NZ)*CO2Qs1(NZ)

  IF(SSINs1.GT.0.0.AND.ARLFPs1(NZ).GT.ZEROPs1(NZ))THEN
!
    call PhotoActivePFT(NZ)
  ELSE
!
    RSMNs1(NZ)=RSMHs1(NZ)
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
    PARDIFs1 => plt_rad%PARDIFs1  , &
    PARs1    => plt_rad%PARs1     , &
    FDBKs1   => plt_photo%FDBKs1  , &
    TAU0s1   => plt_rad%TAU0s1      &
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
  PARX=QNTM*PARDIFs1(N,M,L,NZ)
  PARJ=PARX+ETGROs1(K,NB,NZ)
  ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
  EGRO=ETLF*CBXNs1(K,NB,NZ)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     VGRO=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     CH2O=total rubisco carboxylation rate
!     FDBK=N,P feedback inhibition on C3 CO2 fixation
!     SURFX=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
  VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*FDBKs1(NB,NZ)
  CH2O=CH2O+VL*SURFXs1(N,L,K,NB,NZ)*TAU0s1(L+1)
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
    PARs1     => plt_rad%PARs1        , &
    FDBKs1    => plt_photo%FDBKs1     , &
    TAUSs1    =>  plt_rad%TAUSs1        &
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
  PARX=QNTM*PARs1(N,M,L,NZ)
  PARJ=PARX+ETGROs1(K,NB,NZ)
  ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
  EGRO=ETLF*CBXNs1(K,NB,NZ)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     VGRO=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     FDBK=N,P feedback inhibition on C3 CO2 fixation
!     CH2O=total rubisco carboxylation rate
!     SURFX=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
  VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*FDBKs1(NB,NZ)
  CH2O=CH2O+VL*SURFXs1(N,L,K,NB,NZ)*TAUSs1(L+1)
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
    PARs1    => plt_rad%PARs1        , &
    PARDIFs1 => plt_rad%PARDIFs1       &
  )
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO N=1,JLI1
    DO M=1,JSA1
      IF(SURFXs1(N,L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!     SUNLIT LEAVES
!
        IF(PARs1(N,M,L,NZ).GT.0.0_r8)THEN
          call C3SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
!
!     SHADED LEAVES
        IF(PARDIFs1(N,M,L,NZ).GT.0.0_r8)THEN
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
    ARLFLs1   =>  plt_morph%ARLFLs1 , &
    O2Ls1     =>  plt_photo%O2Ls1   , &
    CO2Ls1    =>  plt_photo%CO2Ls1  , &
    CHLs1     =>  plt_photo%CHLs1   , &
    XKCO2Ls1  =>  plt_photo%XKCO2Ls1, &
    XKCO2Os1  =>  plt_photo%XKCO2Os1, &
    ETMXs1    =>  plt_photo%ETMXs1  , &
    VCMXs1    =>  plt_photo%VCMXs1  , &
    VOMXs1    =>  plt_photo%VOMXs1  , &
    RUBPs1    =>  plt_photo%RUBPs1    &
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
  VCDN=RUBPs1(NZ)*WSDN
  ETDN=CHLs1(NZ)*WSDN
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
  VCGROs1(K,NB,NZ)=VCMXs1(NZ)*TFN1*VCDN
  VOGRO=VOMXs1(NZ)*TFN2*VCDN
  COMPLs1(K,NB,NZ)=0.5_r8*O2Ls1(NZ)*VOGRO*XKCO2Ls1(NZ)/(VCGROs1(K,NB,NZ)*XKO2L)
  VGROs1(K,NB,NZ)=AMAX1(0.0_r8,VCGROs1(K,NB,NZ)*(CO2Ls1(NZ)-COMPLs1(K,NB,NZ)) &
    /(CO2Ls1(NZ)+XKCO2Os1(NZ)))
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
  ETGROs1(K,NB,NZ)=ETMXs1(NZ)*TFNE*ETDN
  CBXNs1(K,NB,NZ)=AMAX1(0.0_r8,(CO2Ls1(NZ)-COMPLs1(K,NB,NZ))/(ELEC3*CO2Ls1(NZ) &
    +10.5_r8*COMPLs1(K,NB,NZ)))
!
!     FOR EACH CANOPY LAYER
!
!     ARLFL=leaf area
!     SURFX=unself-shaded leaf surface area
!
  DO L=JC1,1,-1
    IF(ARLFLs1(L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
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
  associate(                          &
    WGLFs1    =>  plt_biom%WGLFs1   , &
    ARLFLs1   =>  plt_morph%ARLFLs1 , &
    FDBKXs1   =>  plt_photo%FDBKXs1 , &
    CHL4s1    =>  plt_photo%CHL4s1  , &
    O2Ls1     =>  plt_photo%O2Ls1   , &
    CO2Ls1    =>  plt_photo%CO2Ls1  , &
    ETMXs1    =>  plt_photo%ETMXs1  , &
    VCMX4s1   =>  plt_photo%VCMX4s1 , &
    RUBPs1    =>  plt_photo%RUBPs1  , &
    XKCO24s1  =>  plt_photo%XKCO24s1, &
    XKCO2Ls1  =>  plt_photo%XKCO2Ls1, &
    XKCO2Os1  =>  plt_photo%XKCO2Os1, &
    CHLs1     =>  plt_photo%CHLs1   , &
    VCMXs1    =>  plt_photo%VCMXs1  , &
    VOMXs1    =>  plt_photo%VOMXs1  , &
    PEPCs1    =>  plt_photo%PEPCs1    &
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
  CC4M=AMAX1(0.0_r8,0.021E+09_r8*CPOOL4s1(K,NB,NZ)/(WGLFs1(K,NB,NZ)*FMP))
  CCBS=AMAX1(0.0_r8,0.083E+09_r8*CO2Bs1(K,NB,NZ)/(WGLFs1(K,NB,NZ)*FBS))
  FDBK4s1(K,NB,NZ)=1.0_r8/(1.0_r8+CC4M/C4KI)
  FDBK4s1(K,NB,NZ)=FDBK4s1(K,NB,NZ)*FDBKXs1(NB,NZ)
!
!     SURFICIAL DENSITY OF PEPC AND ITS CHLOROPHYLL
!
!     VCDN4=surficial density of PEP carboxylase in mesophyll
!     ETDN4=surficial density of chlorophyll in mesophyll
!     PEPC=fraction of leaf protein in PEP carboxylase
!     CHL4=fraction of leaf protein in mesophyll chlorophyll
!     WSDN=leaf protein surficial density
!
  VCDN4=PEPCs1(NZ)*WSDN
  ETDN4=CHL4s1(NZ)*WSDN
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
  VCGR4s1(K,NB,NZ)=VCMX4s1(NZ)*TFN1*VCDN4
  VGRO4s1(K,NB,NZ)=AMAX1(0.0_r8,VCGR4s1(K,NB,NZ)*(CO2Ls1(NZ)-COMP4)/(CO2Ls1(NZ)+XKCO24s1(NZ)))
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
  ETGR4s1(K,NB,NZ)=ETMXs1(NZ)*TFNE*ETDN4
  CBXN4s1(K,NB,NZ)=AMAX1(0.0_r8,(CO2Ls1(NZ)-COMP4)/(ELEC4*CO2Ls1(NZ)+10.5_r8*COMP4))
!
!     FOR EACH CANOPY LAYER
!
!     ARLFL=leaf area
!     SURFX=unself-shaded leaf surface area
!
  DO L=JC1,1,-1
    IF(ARLFLs1(L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
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
  VCDN=RUBPs1(NZ)*WSDN
  ETDN=CHLs1(NZ)*WSDN
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
  VCGROs1(K,NB,NZ)=VCMXs1(NZ)*TFN1*VCDN
  VOGRO=VOMXs1(NZ)*TFN2*VCDN
  COMPLs1(K,NB,NZ)=0.5_r8*O2Ls1(NZ)*VOGRO*XKCO2Ls1(NZ)/(VCGROs1(K,NB,NZ)*XKO2L)
  VGROs1(K,NB,NZ)=AMAX1(0.0_r8,VCGROs1(K,NB,NZ)*(CCBS-COMPLs1(K,NB,NZ))/(CCBS+XKCO2Os1(NZ)))
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
  ETGROs1(K,NB,NZ)=ETMXs1(NZ)*TFNE*ETDN
  CBXNs1(K,NB,NZ)=AMAX1(0.0_r8,(CCBS-COMPLs1(K,NB,NZ))/(ELEC3*CCBS+10.5*COMPLs1(K,NB,NZ)))
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
    PARs1    => plt_rad%PARs1     , &
    PARDIFs1 => plt_rad%PARDIFs1    &
  )
!
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO  N=1,JLI1
    DO  M=1,JSA1
      IF(SURFXs1(N,L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!     SUNLIT LEAVES
        IF(PARs1(N,M,L,NZ).GT.0.0_r8)THEN
          call C4SunlitLeaves(K,N,M,L,NB,NZ,CH2O)
        ENDIF
!
!     SHADED LEAVES
        IF(PARDIFs1(N,M,L,NZ).GT.0.0_r8)THEN
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
    PARs1    => plt_rad%PARs1       , &
    PARDIFs1 => plt_rad%PARDIFs1    , &
    TAU0s1   => plt_rad%TAU0s1        &
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
  PARX=QNTM*PARDIFs1(N,M,L,NZ)
  PARJ=PARX+ETGR4s1(K,NB,NZ)
  ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4s1(K,NB,NZ)))/CURV2
  EGRO4=ETLF4*CBXN4s1(K,NB,NZ)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     VGRO4=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     CH2O=total PEP carboxylation rate
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!     SURFX=unself-shaded leaf surface area
!     TAU0=fraction of diffuse radiation transmitted from layer above
!
  VL=AMIN1(VGRO4s1(K,NB,NZ),EGRO4)*FDBK4s1(K,NB,NZ)
  CH2O=CH2O+VL*SURFXs1(N,L,K,NB,NZ)*TAU0s1(L+1)
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
  associate(                      &
  PARs1   =>  plt_rad%PARs1     , &
  TAUSs1  =>  plt_rad%TAUSs1      &
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
  PARX=QNTM*PARs1(N,M,L,NZ)
  PARJ=PARX+ETGR4s1(K,NB,NZ)
  ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4s1(K,NB,NZ)))/CURV2
  EGRO4=ETLF4*CBXN4s1(K,NB,NZ)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     VGRO4=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     FDBK4=N,P feedback inhibition on C4 CO2 fixation
!     CH2O=total PEP carboxylation rate
!     SURFX=unself-shaded leaf surface area
!     TAUS=fraction of direct radiation transmitted from layer above
!
  VL=AMIN1(VGRO4s1(K,NB,NZ),EGRO4)*FDBK4s1(K,NB,NZ)
  CH2O=CH2O+VL*SURFXs1(N,L,K,NB,NZ)*TAUSs1(L+1)
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
    WGLFs1     => plt_biom%WGLFs1     , &
    WSLFs1     => plt_biom%WSLFs1     , &
    ARLF1s1    => plt_morph%ARLF1s1     &
  )
  DO K=1,JNODS1
    IF(ARLF1s1(K,NB,NZ).GT.ZEROPs1(NZ).AND.WGLFs1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
      WSDN=WSLFs1(K,NB,NZ)/ARLF1s1(K,NB,NZ)
    ELSE
      WSDN=0.0
    ENDIF

    IF(WSDN.GT.ZEROs1)THEN
!
!     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!
      IF(ICTYPs1(NZ).EQ.4)THEN
!     C4 PHOTOSYNTHESIS
        call C4Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
      ELSE
!     C3 PHOTOSYNTHESIS
        call C3Photosynthesis(K,NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L,WSDN)
      ENDIF
!
    ELSE
      VCGR4s1(K,NB,NZ)=0.0
      VCGROs1(K,NB,NZ)=0.0
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
    CCPOLBs1 =>  plt_biom%CCPOLBs1  , &
    CZPOLBs1 =>  plt_biom%CZPOLBs1  , &
    CPPOLBs1 =>  plt_biom%CPPOLBs1  , &
    IWTYPs1  =>  plt_pheno%IWTYPs1  , &
    IBTYPs1  =>  plt_pheno%IBTYPs1  , &
    FLG4s1   =>  plt_pheno%FLG4s1   , &
    ATRPs1   =>  plt_pheno%ATRPs1   , &
    ISTYPs1  =>  plt_pheno%ISTYPs1  , &
    IDTHBs1  =>  plt_pheno%IDTHBs1  , &
    FDBKXs1   => plt_photo%FDBKXs1  , &
    FDBKs1   => plt_photo%FDBKs1      &
  )
!
!     FEEDBACK ON C3 CARBOXYLATION FROM NON-STRUCTURAL C:N:P
!
!     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
!     FDBK=N,P feedback inhibition on C3 CO2 fixation
!     CNKI,CPKI=nonstructural N,P inhibition constant on rubisco
!
  IF(CCPOLBs1(NB,NZ).GT.ZEROs1)THEN
    FDBKs1(NB,NZ)=AMIN1(CZPOLBs1(NB,NZ)/(CZPOLBs1(NB,NZ)+CCPOLBs1(NB,NZ)/CNKI) &
      ,CPPOLBs1(NB,NZ)/(CPPOLBs1(NB,NZ)+CCPOLBs1(NB,NZ)/CPKI))
  ELSE
    FDBKs1(NB,NZ)=1.0
  ENDIF
!
!     CHILLING
!
!     CHILL=accumulated chilling hours used to limit CO2 fixn
!
!     FDBKs1(NB,NZ)=FDBKs1(NB,NZ)/(1.0_r8+0.25_r8*CHILLs1(NZ))
!
!     DEHARDENING OF EVERGREENS IN SPRING
!
!     ATRP=hours above threshold temperature for dehardening since leafout
!     ATRPZ=hours to full dehardening of conifers in spring
!
  IF(IWTYPs1(NZ).NE.0.AND.IBTYPs1(NZ).GE.2)THEN
    FDBKs1(NB,NZ)=FDBKs1(NB,NZ)*AMAX1(0.0_r8,AMIN1(1.0_r8 &
      ,ATRPs1(NB,NZ)/(0.9_r8*ATRPZ)))
  ENDIF
!
!     TERMINATION OF ANNUALS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     FLG4=number of hours with no grain fill after start of grain fill
!     FLG4Y=number of hours with no grain fill to terminate annuals
!
  IF(ISTYPs1(NZ).EQ.0.AND.FLG4s1(NB,NZ).GT.0.0)THEN
    FDBKXs1(NB,NZ)=AMAX1(0.0_r8,1.0_r8-FLG4s1(NB,NZ)/FLG4Y(IWTYPs1(NZ)))
  ELSE
    FDBKXs1(NB,NZ)=1.0
  ENDIF
  FDBKs1(NB,NZ)=FDBKs1(NB,NZ)*FDBKXs1(NB,NZ)
!
!     FOR EACH NODE
!
!     IDTHB=branch life flag:0=living,1=dead
!     ARLF,WGLF,WSLF=leaf area,C mass,protein mass
!     WSDN=leaf protein surficial density
!
  IF(IDTHBs1(NB,NZ).EQ.0)THEN
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
    CCPOLBs1 =>  plt_biom%CCPOLBs1  , &
    XKO2s1   =>  plt_photo%XKO2s1   , &
    SCO2s1   =>  plt_photo%SCO2s1   , &
    CO2Qs1   =>  plt_photo%CO2Qs1   , &
    O2Ls1    =>  plt_photo%O2Ls1    , &
    CO2Ls1   =>  plt_photo%CO2Ls1   , &
    CO2Is1   =>  plt_photo%CO2Is1   , &
    O2Is1    =>  plt_photo%O2Is1    , &
    RSMHs1   =>  plt_photo%RSMHs1   , &
    XKCO2Os1 =>  plt_photo%XKCO2Os1 , &
    XKCO2Ls1 =>  plt_photo%XKCO2Ls1 , &
    XKCO2s1  =>  plt_photo%XKCO2s1    &

  )
!
!     CO2 AND O2 AQUEOUS SOLUBILITY
!
!     TCCZ=canopy temperature
!     SCO2,SO2=solubility of CO2,O2 (uM/(umol mol-1))
!     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
!     DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!
  TCCZ=TKCZs1(NZ)-TC2K
  SCO2s1(NZ)=EXP(-2.621_r8-0.0317_r8*TCCZ)
  SO2s1(NZ)=EXP(-6.175_r8-0.0211_r8*TCCZ)
  CO2Ls1(NZ)=CO2Is1(NZ)*SCO2s1(NZ)
  O2Ls1(NZ)=O2Is1(NZ)*SO2s1(NZ)
  DCO2s1(NZ)=FMOLs1(NZ)*(CO2Qs1(NZ)-CO2Is1(NZ))
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
  TKCO=TKCZs1(NZ)+OFFSTs1(NZ)
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
  XKCO2Ls1(NZ)=XKCO2s1(NZ)*EXP(16.136_r8-40000._r8/RTK)
  XKO2L=XKO2s1(NZ)*EXP(8.067_r8-20000._r8/RTK)
  XKCO2Os1(NZ)=XKCO2Ls1(NZ)*(1.0_r8+O2Ls1(NZ)/XKO2L)
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
  associate(                         &
    VRNXs1   =>  plt_pheno%VRNXs1  , &
    VRNFs1   =>  plt_pheno%VRNFs1  , &
    VRNLs1   =>  plt_pheno%VRNLs1  , &
    VRNSs1   =>  plt_pheno%VRNSs1  , &
    IWTYPs1  =>  plt_pheno%IWTYPs1  , &
    FDBKs1   => plt_photo%FDBKs1   , &
    FDBKXs1  => plt_photo%FDBKXs1  , &
    XKCO2Os1 => plt_photo%XKCO2Os1 , &
    RSMHs1   => plt_photo%RSMHs1   , &
    FRADPs1  => plt_rad%FRADPs1    , &
    NBRs1    => plt_morph%NBRs1      &
  )

  call PrepPhotosynthesis(NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
!
!     FOR EACH BRANCH
!
  DO NB=1,NBRs1(NZ)
!
!     FEEDBACK ON CO2 FIXATION
!
!     IWTYP=phenology type from PFT file
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!
    IF(IWTYPs1(NZ).EQ.0.OR.VRNSs1(NB,NZ).GE.VRNLs1(NB,NZ) &
      .OR.VRNFs1(NB,NZ).LT.VRNXs1(NB,NZ))THEN

      call PhenoActiveBranch(NB,NZ,CH2O,TFN1,TFN2,TFNE,XKO2L)
    ELSE
      FDBKs1(NB,NZ)=0.0_r8
      FDBKXs1(NB,NZ)=1.0_r8
      DO K=1,JNODS1
        VCGR4s1(K,NB,NZ)=0.0_r8
        VCGROs1(K,NB,NZ)=0.0_r8
      ENDDO
    ENDIF
  ENDDO
!
!     MINIMUM CANOPY STOMATAL RESISTANCE FROM CO2 CONCENTRATION
!     DIFFERENCE DIVIDED BY TOTAL CO2 FIXATION
!
!     RSX,RSMN=minimum canopy stomatal resistance to CO2,H2O (h m-1)
!     CH2O=total PEP(C4) or rubisco(C3) carboxylation rate
!     FRADP=fraction of radiation received by each PFT canopy
!     DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!     AREA=area of grid cell
!     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
! hourly time step
  IF(CH2O.GT.ZEROPs1(NZ))THEN
    RSX=FRADPs1(NZ)*DCO2s1(NZ)*AREA3s1(NUs1)/(CH2O*secsperhour)
  ELSE
    RSX=RSMHs1(NZ)*1.56_r8
  ENDIF
  RSMNs1(NZ)=AMIN1(RSMHs1(NZ),AMAX1(RSMY,RSX*0.641_r8))
  end associate
  end subroutine PhotoActivePFT

end module stomatesMod
