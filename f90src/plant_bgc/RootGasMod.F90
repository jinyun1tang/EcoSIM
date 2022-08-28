module RootGasMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StomatesMod   , only : stomates
  use minimathmod  , only : safe_adb,vapsat,test_aneb
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use PlantAPIData
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = __FILE__
  public :: RootSoilGasExchange
  contains
!------------------------------------------------------------------------

  subroutine RootSoilGasExchange(N,L,NZ,RRADL,FPQ,FRTDPX,RTARR,&
    UPWTRH,FOXYX,RUPOXT)

  implicit none
  integer , intent(in) :: N,L,NZ
  real(r8), intent(in) :: RRADL(2,JZ1),FPQ(2,JZ1,JP1),FRTDPX(JZ1,JP1)
  real(r8), intent(in) :: RTARR(2,JZ1),UPWTRH,FOXYX
  real(r8), intent(out):: RUPOXT
  integer :: M,MX
  real(r8) :: B,C
  real(r8) :: CO2A1,CO2P1,CO2G1,CO2S1,CH4A1
  real(r8) :: CH4P1,CH4S1,CCH4S1,CCH4P1
  real(r8) :: CN2OS1,CN2OP1,CNH3S1,CNH3B1,CNH3P1,CH2GS1,CH2GP1
  real(r8) :: CGSGL1,CHSGL1,CLSGL1
  real(r8) :: CQSGL1,CH4G1,CCO2S1,COXYS1,CCO2A1,COXYA1,CCH4A1
  real(r8) :: CZ2OA1,CNH3A1,CH2GA1,CCO2P1,COXYP1,COXYR,CO2PX,CH4PX
  real(r8) :: DIFOP,DIFCL,DIFZL,DIFNL,DIFNB,DIFHL,DIFOX,DFGSP
  real(r8) :: DFCOA,DFOXA,DFCHA,DFN2A,DFNHA,DFHGA,DFGP,DIFOL
  real(r8) :: H2GA1,H2GP1,H2GS1,HGSGL1,HLSGL1,H2GG1,H2GPX
  real(r8) :: OXYA1,OXYP1,OXYG1,OXYS1,OGSGL1
  real(r8) :: OLSGL1,OLSGLP,OXYPX,RTVLWA
  real(r8) :: RTVLWB,ROXYFX,RCO2FX,ROXYLX,ROXDFQ,RCHDFQ,RN2DFQ
  real(r8) :: RNHDFQ,RHGDFQ,ROXDF1,RCHDF1,RN2DF1,RNHDF1,RHGDF1
  real(r8) :: RTCR1,RTCR2,RTCRA,RTARRX,RCO2PX,RRADS,RMFCOS,RMFOXS
  real(r8) :: RMFCHS,RMFN2S,RMFN3S,RMFN3B,RMFHGS,RUPOXR,RDFOXS
  real(r8) :: RDFOXP,RUPOSX,RUPOPX,RDFCOS,RDXCOS,RCO2SX,RDFCHS
  real(r8) :: RDXCHS,RUPCSX,RDFN2S,RDXN2S,RUPZSX,RDFN3S,RDXNHS
  real(r8) :: RUPNSX,RDFN3B,RDXNHB,RUPNBX,RDFHGS,RDXHGS,RUPHGX
  real(r8) :: RCODFQ,RUPOST,RNBDFQ,RUPNTX,RCODF1,RCOFL1,ROXFL1
  real(r8) :: RCHFL1,RN2FL1,RNHFL1,RHGFL1,THETW1,THETM
  real(r8) :: UPMXP
  real(r8) :: VOLWCA,VOLWOA
  real(r8) :: VOLWC4,VOLWZA,VOLWNA,VOLWH2,VOLWMO,VOLWMM,VOLPMM
  real(r8) :: VOLWSP,VOLWMA,VOLWMB,VOLWSA,VOLWSB,VOLWCO,VOLWOX
  real(r8) :: VOLWCH,VOLWN2,VOLWNH,VOLWNB,VOLWHG,VOLPNH,VOLPNB
  real(r8) :: X
  real(r8) :: Z2OA1,Z2OP1,Z2OS1,ZH3A1,ZH3P1,ZH3S1
  real(r8) :: ZH3B1,Z2SGL1,ZHSGL1,ZVSGL1,ZNSGL1,Z2OG1,ZH3G1
  real(r8) :: ZH3PA,ZH3PB,ZH3GA,ZH3GB,Z2OPX,ZH3PX

!     begin_execution
  associate(                           &
    WTRTSs1  =>  plt_biom%WTRTSs1    , &
    ZEROPs1  =>  plt_biom%ZEROPs1    , &
    PPs1     =>  plt_site%PPs1       , &
    DPTHZs1  =>  plt_site%DPTHZs1    , &
    CCH4Es1  =>  plt_site%CCH4Es1    , &
    CCO2Es1  =>  plt_site%CCO2Es1    , &
    CH2GEs1  =>  plt_site%CH2GEs1    , &
    CNH3Es1  =>  plt_site%CNH3Es1    , &
    COXYEs1  =>  plt_site%COXYEs1    , &
    CZ2OEs1  =>  plt_site%CZ2OEs1    , &
    ZEROSs1  =>  plt_site%ZEROSs1    , &
    ZEROs1   =>  plt_site%ZEROs1     , &
    VOLWMs1  =>  plt_site%VOLWMs1    , &
    VOLPMs1  =>  plt_site%VOLPMs1    , &
    TORTs1   =>  plt_site%TORTs1     , &
    FILMs1   =>  plt_site%FILMs1     , &
    ROXYFs1  =>  plt_bgcr%ROXYFs1    , &
    RCO2Fs1  =>  plt_bgcr%RCO2Fs1    , &
    ROXYLs1  =>  plt_bgcr%ROXYLs1    , &
    WFRs1    =>  plt_rbgc%WFRs1      , &
    ZEROQs1  =>  plt_rbgc%ZEROQs1    , &
    RCO2Ms1  =>  plt_rbgc%RCO2Ms1    , &
    ROXYPs1  =>  plt_rbgc%ROXYPs1    , &
    RUPOXSs1 =>  plt_rbgc%RUPOXSs1   , &
    ROXSKs1  =>  plt_rbgc%ROXSKs1    , &
    RCO2Ps1  =>  plt_rbgc%RCO2Ps1    , &
    RHGFLAs1 =>  plt_rbgc%RHGFLAs1   , &
    RNHDFAs1 =>  plt_rbgc%RNHDFAs1   , &
    RUPN3Bs1 =>  plt_rbgc%RUPN3Bs1   , &
    RHGDFAs1 =>  plt_rbgc%RHGDFAs1   , &
    ROXFLAs1 =>  plt_rbgc%ROXFLAs1   , &
    RCHFLAs1 =>  plt_rbgc%RCHFLAs1   , &
    RN2FLAs1 =>  plt_rbgc%RN2FLAs1   , &
    RCOFLAs1 =>  plt_rbgc%RCOFLAs1   , &
    RCODFAs1 =>  plt_rbgc%RCODFAs1   , &
    RNHFLAs1 =>  plt_rbgc%RNHFLAs1   , &
    RUPHGSs1 =>  plt_rbgc%RUPHGSs1   , &
    RCO2Ss1  =>  plt_rbgc%RCO2Ss1    , &
    ROXDFAs1 =>  plt_rbgc%ROXDFAs1   , &
    RCHDFAs1 =>  plt_rbgc%RCHDFAs1   , &
    RN2DFAs1 =>  plt_rbgc%RN2DFAs1   , &
    RUPCHSs1 =>  plt_rbgc%RUPCHSs1   , &
    RUPN2Ss1 =>  plt_rbgc%RUPN2Ss1   , &
    RUPN3Ss1 =>  plt_rbgc%RUPN3Ss1   , &
    RUPOXPs1 =>  plt_rbgc%RUPOXPs1   , &
    RCO2As1  =>  plt_rbgc%RCO2As1    , &
    CO2Ps1   =>  plt_rbgc%CO2Ps1     , &
    CO2As1   =>  plt_rbgc%CO2As1     , &
    OXYAs1   =>  plt_rbgc%OXYAs1     , &
    OXYPs1   =>  plt_rbgc%OXYPs1     , &
    CH4As1   =>  plt_rbgc%CH4As1     , &
    CH4Ps1   =>  plt_rbgc%CH4Ps1     , &
    Z2OAs1   =>  plt_rbgc%Z2OAs1     , &
    Z2OPs1   =>  plt_rbgc%Z2OPs1     , &
    ZH3As1   =>  plt_rbgc%ZH3As1     , &
    ZH3Ps1   =>  plt_rbgc%ZH3Ps1     , &
    H2GAs1   =>  plt_rbgc%H2GAs1     , &
    H2GPs1   =>  plt_rbgc%H2GPs1     , &
    TFNDs1   =>  plt_soilchem%TFNDs1 , &
    CO2Gs1   =>  plt_soilchem%CO2Gs1 , &
    CO2Ss1   =>  plt_soilchem%CO2Ss1 , &
    ZNH3Ss1  =>  plt_soilchem%ZNH3Ss1, &
    SN2OLs1  =>  plt_soilchem%SN2OLs1, &
    H2GSs1   =>  plt_soilchem%H2GSs1 , &
    CZ2OSs1  =>  plt_soilchem%CZ2OSs1, &
    CNH3Gs1  =>  plt_soilchem%CNH3Gs1, &
    CH2GGs1  =>  plt_soilchem%CH2GGs1, &
    CH2GSs1  =>  plt_soilchem%CH2GSs1, &
    CZ2OGs1  =>  plt_soilchem%CZ2OGs1, &
    VLNH4s1  =>  plt_soilchem%VLNH4s1, &
    HGSGLs1  =>  plt_soilchem%HGSGLs1, &
    SCH4Ls1  =>  plt_soilchem%SCH4Ls1, &
    CCH4Gs1  =>  plt_soilchem%CCH4Gs1, &
    CLSGLs1  =>  plt_soilchem%CLSGLs1, &
    THETYs1  =>  plt_soilchem%THETYs1, &
    OLSGLs1  =>  plt_soilchem%OLSGLs1, &
    SOXYLs1  =>  plt_soilchem%SOXYLs1, &
    SNH3Ls1  =>  plt_soilchem%SNH3Ls1, &
    SH2GLs1  =>  plt_soilchem%SH2GLs1, &
    VOLYs1   =>  plt_soilchem%VOLYs1 , &
    ZNSGLs1  =>  plt_soilchem%ZNSGLs1, &
    Z2SGLs1  =>  plt_soilchem%Z2SGLs1, &
    ZHSGLs1  =>  plt_soilchem%ZHSGLs1, &
    ZVSGLs1  =>  plt_soilchem%ZVSGLs1, &
    SCO2Ls1  =>  plt_soilchem%SCO2Ls1, &
    HLSGLs1  =>  plt_soilchem%HLSGLs1, &
    CQSGLs1  =>  plt_soilchem%CQSGLs1, &
    VLNHBs1  =>  plt_soilchem%VLNHBs1, &
    ZNH3Bs1  =>  plt_soilchem%ZNH3Bs1, &
    CGSGLs1  =>  plt_soilchem%CGSGLs1, &
    CNH3Bs1  =>  plt_soilchem%CNH3Bs1, &
    CHSGLs1  =>  plt_soilchem%CHSGLs1, &
    OGSGLs1  =>  plt_soilchem%OGSGLs1, &
    Z2OSs1   =>  plt_soilchem%Z2OSs1 , &
    CNH3Ss1  =>  plt_soilchem%CNH3Ss1, &
    THETPMs1 =>  plt_soilchem%THETPMs1,&
    OXYGs1   =>  plt_soilchem%OXYGs1 , &
    DFGSs1   =>  plt_soilchem%DFGSs1 , &
    CCH4Ss1  =>  plt_soilchem%CCH4Ss1, &
    OXYSs1   =>  plt_soilchem%OXYSs1 , &
    CH4Ss1   =>  plt_soilchem%CH4Ss1 , &
    IDAYs1   =>  plt_pheno%IDAYs1    , &
    PORTXs1  =>  plt_morph%PORTXs1   , &
    RRAD1s1  =>  plt_morph%RRAD1s1   , &
    RTLGAs1  =>  plt_morph%RTLGAs1   , &
    RTVLPs1  =>  plt_morph%RTVLPs1   , &
    RTLGPs1  =>  plt_morph%RTLGPs1   , &
    RTNLs1   =>  plt_morph%RTNLs1    , &
    RRAD2s1  =>  plt_morph%RRAD2s1   , &
    RRADPs1  =>  plt_morph%RRADPs1   , &
    RTVLWs1  =>  plt_morph%RTVLWs1   , &
    PORTs1   =>  plt_morph%PORTs1    , &
    RTN1s1   =>  plt_morph%RTN1s1    , &
    NGs1     =>  plt_morph%NGs1      , &
    NB1s1    =>  plt_morph%NB1s1       &
  )
  IF(RCO2Ms1(N,L,NZ).GT.ZEROPs1(NZ).AND.RTVLWs1(N,L,NZ).GT.ZEROPs1(NZ) &
    .AND.FOXYX.GT.ZEROQs1(NZ))THEN
!
!     INITIALIZE VARIABLES USED IN ROOT GAS EXCHANGE
!     (CO2, O2, CH4, N2, N2O, NH3, H2)
!
!     CO2A1,CO2P1,CO2G1,CO2S1=gaseous,aqueous CO2 in root,soil
!     OXYA1,OXYP1,OXYG1,OXYS1=gaseous,aqueous O2 in root,soil
!     CH4A1,CH4P1,CH4G1,CH4S1=gaseous,aqueous CH4 in root,soil
!     Z2OA1,Z2OP1,Z2OG1,Z2OS1=gaseous,aqueous N2O in root,soil
!     ZH3A1,ZH3P1,ZH3G1,ZH3S1=gaseous,aqueous NH3 in root,soil
!     H2GA1,H2GP1,H2GG1,H2GS1=gaseous,aqueous H2 in root,soil
!     CCH4S1,CCH4P1=aqueous CH4 concentration in soil,root
!     CN2OS1,CN2OP1=aqueous N2O concentration in soil,root
!     CNH3S1,CNH3B1,CNH3P1=aqueous NH3 concn in soil non-band,band,root
!     CH2GS1,CH2GP1=aqueous H2 concentration in soil,root
!     RTVLWA,RTVLWB=root aqueous volume in non-band,band
!     XNPG=time step of flux calculation
!     UPMXP=O2 demand per plant at time step of flux calculation
!     ROXYFX=net O2 gas flux at time step of flux calculation
!     RCO2FX=net CO2 gas flux at time step of flux calculation
!     ROXYLX=net O2 aqueous flux at time step of flux calculation
!
    CO2A1=AMAX1(ZEROPs1(NZ),CO2As1(N,L,NZ))
    CO2P1=AMAX1(ZEROPs1(NZ),CO2Ps1(N,L,NZ))
    CO2G1=AMAX1(ZEROPs1(NZ),CO2Gs1(L)*FPQ(N,L,NZ))
    CO2S1=AMAX1(ZEROPs1(NZ),CO2Ss1(L)*FPQ(N,L,NZ))
    OXYA1=AMAX1(ZEROPs1(NZ),OXYAs1(N,L,NZ))
    OXYP1=AMAX1(ZEROPs1(NZ),OXYPs1(N,L,NZ))
    OXYG1=AMAX1(ZEROPs1(NZ),OXYGs1(L)*FOXYX)
    OXYS1=OXYSs1(L)*FOXYX
    CH4A1=CH4As1(N,L,NZ)
    CH4P1=CH4Ps1(N,L,NZ)
    CH4S1=CH4Ss1(L)*FPQ(N,L,NZ)
    CCH4S1=CCH4Ss1(L)
    CCH4P1=AMAX1(0.0,CH4P1/RTVLWs1(N,L,NZ))
    Z2OA1=Z2OAs1(N,L,NZ)
    Z2OP1=Z2OPs1(N,L,NZ)
    Z2OS1=Z2OSs1(L)*FPQ(N,L,NZ)
    CN2OS1=CZ2OSs1(L)
    CN2OP1=AMAX1(0.0,Z2OP1/RTVLWs1(N,L,NZ))
    ZH3A1=ZH3As1(N,L,NZ)
    ZH3P1=ZH3Ps1(N,L,NZ)
    ZH3S1=ZNH3Ss1(L)*FPQ(N,L,NZ)
    ZH3B1=ZNH3Bs1(L)*FPQ(N,L,NZ)
    CNH3S1=CNH3Ss1(L)
    CNH3B1=CNH3Bs1(L)
    CNH3P1=AMAX1(0.0,ZH3P1/RTVLWs1(N,L,NZ))
    H2GA1=H2GAs1(N,L,NZ)
    H2GP1=H2GPs1(N,L,NZ)
    H2GS1=H2GSs1(L)*FPQ(N,L,NZ)
    CH2GS1=CH2GSs1(L)
    CH2GP1=AMAX1(0.0,H2GP1/RTVLWs1(N,L,NZ))
    RTVLWA=RTVLWs1(N,L,NZ)*VLNH4s1(L)
    RTVLWB=RTVLWs1(N,L,NZ)*VLNHBs1(L)
    UPMXP=ROXYPs1(N,L,NZ)*XNPG/PPs1(NZ)
    ROXYFX=ROXYFs1(L)*FOXYX*XNPG
    RCO2FX=RCO2Fs1(L)*FOXYX*XNPG
    ROXYLX=ROXYLs1(L)*FOXYX*XNPG
!
!     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
!
!     *SGL1=diffusivity
!     PORTX=tortuosity effect of root porosity on diffusivity
!     CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
!     CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
!
    CGSGL1=CGSGLs1(L)*XNPG*PORTXs1(N,NZ)
    OGSGL1=OGSGLs1(L)*XNPG*PORTXs1(N,NZ)
    CHSGL1=CHSGLs1(L)*XNPG*PORTXs1(N,NZ)
    Z2SGL1=Z2SGLs1(L)*XNPG*PORTXs1(N,NZ)
    ZHSGL1=ZHSGLs1(L)*XNPG*PORTXs1(N,NZ)
    HGSGL1=HGSGLs1(L)*XNPG*PORTXs1(N,NZ)
    CLSGL1=CLSGLs1(L)*XNPG*FOXYX
    OLSGL1=OLSGLs1(L)*XNPG*FOXYX
    CQSGL1=CQSGLs1(L)*XNPG*FOXYX
    ZVSGL1=ZVSGLs1(L)*XNPG*FOXYX
    ZNSGL1=ZNSGLs1(L)*XNPG*FOXYX
    HLSGL1=HLSGLs1(L)*XNPG*FOXYX
    OLSGLP=OLSGLs1(L)*XNPG
    ROXDFQ=0.0_r8
    RCHDFQ=0.0_r8
    RN2DFQ=0.0_r8
    RNHDFQ=0.0_r8
    RHGDFQ=0.0_r8
    ROXDF1=0.0_r8
    RCHDF1=0.0_r8
    RN2DF1=0.0_r8
    RNHDF1=0.0_r8
    RHGDF1=0.0_r8
!
!     ROOT CONDUCTANCE TO GAS TRANSFER
!
!     WTRTS=total root,myco mass
!     FRTDPX=fraction of each soil layer with primary root
!     RTCR1,RTCR2,RTCRA=cross-sectional area/length of
!     primary,secondary,total root,myco system
!     RTN1,RTNL=number of root,myco primary,secondary axes
!     RRAD1,RRAD2=primary,secondary root radius
!     DPTHZ=depth of primary root from surface
!     RTLGA=average secondary root length
!
    IF(WTRTSs1(NZ).GT.ZEROPs1(NZ).AND.FRTDPX(L,NZ).GT.ZEROs1)THEN
      RTCR1=AMAX1(PPs1(NZ),RTN1s1(N,L,NZ)) &
        *PICON*RRAD1s1(N,L,NZ)**2/DPTHZs1(L)
      RTCR2=(RTNLs1(N,L,NZ)*PICON*RRAD2s1(N,L,NZ)**2 &
        /RTLGAs1(N,L,NZ))/FRTDPX(L,NZ)
      IF(RTCR2.GT.RTCR1)THEN
        RTCRA=RTCR1*RTCR2/(RTCR1+RTCR2)
      ELSE
        RTCRA=RTCR1
      ENDIF
    ELSE
      RTCRA=0.0_r8
    ENDIF
!
!     VARIABLES USED TO CALCULATE ROOT GAS TRANSFER
!     BETWEEN AQUEOUS AND GASEOUS PHASES
!
!     RTLGP=root,myco length per plant
!     IDAYs1(1,=emergence date
!     RTARR=root surface area/radius for uptake
!     RRADP=path length for radial diffusion within root
!     DIFOP=aqueous diffusivity of O2 within root
!     RTVLW=root,myco aqueous volume
!     S*L=solubility of gas in water from hour1.f:
!     CO2=CO2,OXY=O2,CH4=CH4,N2O=N2O,NH3=NH3,H2G=H2
!     DF*A=root-atmosphere gas conductance
!     DFGP=rate const for equilibrn of gas concn in gaseous-aqueous phases
!     RCO2PX=root CO2 gas flux at time step for gas flux calculations
!     RCO2A=root CO2 flux from grosub.f
!
    IF(N.EQ.1.AND.IDAYs1(1,NB1s1(NZ),NZ).GT.0 &
      .AND.RTLGPs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
      RTARRX=RTARR(N,L)/RRADPs1(N,NZ)
      DIFOP=OLSGLP*RTARRX
      VOLWCA=RTVLWs1(N,L,NZ)*SCO2Ls1(L)
      VOLWOA=RTVLWs1(N,L,NZ)*SOXYLs1(L)
      VOLWC4=RTVLWs1(N,L,NZ)*SCH4Ls1(L)
      VOLWZA=RTVLWs1(N,L,NZ)*SN2OLs1(L)
      VOLWNA=RTVLWs1(N,L,NZ)*SNH3Ls1(L)
      VOLWH2=RTVLWs1(N,L,NZ)*SH2GLs1(L)
      DFCOA=CGSGL1*RTCRA
      DFOXA=OGSGL1*RTCRA
      DFCHA=CHSGL1*RTCRA
      DFN2A=Z2SGL1*RTCRA
      DFNHA=ZHSGL1*RTCRA
      DFHGA=HGSGL1*RTCRA
    ELSE
      RTARRX=0.0_r8
      DIFOP=0.0_r8
      VOLWCA=0.0_r8
      VOLWOA=0.0_r8
      VOLWC4=0.0_r8
      VOLWZA=0.0_r8
      VOLWNA=0.0_r8
      VOLWH2=0.0_r8
      DFCOA=0.0_r8
      DFOXA=0.0_r8
      DFCHA=0.0_r8
      DFN2A=0.0_r8
      DFNHA=0.0_r8
      DFHGA=0.0_r8
    ENDIF
    DFGP=AMIN1(1.0,XNPD*SQRT(PORTs1(N,NZ))*TFNDs1(L))
    RCO2PX=-RCO2As1(N,L,NZ)*XNPG
!
!     SOLVE FOR GAS EXCHANGE IN SOIL AND ROOTS DURING ROOT UPTAKE
!     AT SMALLER TIME STEP NPH
!
    DO 99 M=1,NPH
!
!     AQUEOUS GAS DIFFUSIVITY THROUGH SOIL WATER TO ROOT
!
!     gas code:CO2=CO2,OXY=O2,CH4=CH4,Z2O=N2O,NH3=NH3 non-band,
!     NHB=NH3 band,H2G=H2
!     VOLWMM,VOLPMM=soil micropore water,air volume
!     FOXYX=root fraction of total O2 demand from previous hour
!     FPQ=PFT fraction of biome root mass
!     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
!     VOLX=soil volume excluding rock,macropores
!     THETW1=soil water concentration
!     TORT=soil tortuosity
!     FILM=soil water film thickness
!     RRADL=root radius
!     RRADS=path length for radial diffusion from soil to root
!     DIF*=aqueous diffusivity from soil to root:OL=O2,CL=CH4
!     ZL=N2O,NL=NH3 non-band,NB=NH4 band,HL=H2
!     C*G=soil gaseous concentration
!     VOLW*,VOLP*=VOLWMM,VOLPMM*gas solubility
!
      VOLWMO=VOLWMs1(M,L)*FOXYX
      VOLWMM=VOLWMs1(M,L)*FPQ(N,L,NZ)
      VOLPMM=VOLPMs1(M,L)*FPQ(N,L,NZ)
      VOLWSP=RTVLWs1(N,L,NZ)+VOLWMM
      VOLWMA=VOLWMM*VLNH4s1(L)
      VOLWMB=VOLWMM*VLNHBs1(L)
      VOLWSA=RTVLWA+VOLWMA
      VOLWSB=RTVLWB+VOLWMB
      THETW1=AMAX1(0.0,VOLWMs1(M,L)/VOLYs1(L))
      IF(THETW1.GT.THETYs1(L).AND.FPQ(N,L,NZ).GT.ZEROQs1(NZ))THEN
        THETM=TORTs1(M,L)*THETW1
        RRADS=LOG((FILMs1(M,L)+RRADL(N,L))/RRADL(N,L))
        RTARRX=RTARR(N,L)/RRADS
        DIFOL=THETM*OLSGL1*RTARRX
        DIFCL=THETM*CQSGL1*RTARRX
        DIFZL=THETM*ZVSGL1*RTARRX
        DIFNL=THETM*ZNSGL1*RTARRX*VLNH4s1(L)
        DIFNB=THETM*ZNSGL1*RTARRX*VLNHBs1(L)
        DIFHL=THETM*HLSGL1*RTARRX
        CH4G1=CCH4Gs1(L)*VOLPMM
        Z2OG1=CZ2OGs1(L)*VOLPMM
        ZH3G1=CNH3Gs1(L)*VOLPMM
        H2GG1=CH2GGs1(L)*VOLPMM
        VOLWCO=VOLWMM*SCO2Ls1(L)
        VOLWOX=VOLWMM*SOXYLs1(L)
        VOLWCH=VOLWMM*SCH4Ls1(L)
        VOLWN2=VOLWMM*SN2OLs1(L)
        VOLWNH=VOLWMM*SNH3Ls1(L)*VLNH4s1(L)
        VOLWNB=VOLWMM*SNH3Ls1(L)*VLNHBs1(L)
        VOLWHG=VOLWMM*SH2GLs1(L)
        VOLPNH=VOLPMM*VLNH4s1(L)
        VOLPNB=VOLPMM*VLNHBs1(L)
!
!     MASS FLOW OF GAS FROM SOIL TO ROOT AT SHORTER TIME STEP NPT
!
!     C*S1=soil aqueous concentration non-band
!     C*B1=soil aqueous concentration band
!     C*A1=root gaseous concentration
!     C*P1=root aqueous concentration
!     ROXYLX=soil net O2 aqueous flux
!     VOLWMM=micropore water volume
!     RTVLW,RTVLP=root aqueous,gaseous volume
!     RMF*=soil convective solute flux:COS=CO2,OXS=O2,CHS=CH4,
!     N2S=N2O,NHS=NH3 non-band,NHB=NH3 band,HGS=H2
!     UPWTRH=water uptake
!
        DO 90 MX=1,NPT
          OXYS1=OXYS1+ROXYLX
          CCO2S1=AMAX1(0.0,CO2S1/VOLWMM)
          COXYS1=AMIN1(COXYEs1*SOXYLs1(L),AMAX1(0.0,OXYS1/VOLWMO))
          CCH4S1=AMAX1(0.0,CH4S1/VOLWMM)
          CN2OS1=AMAX1(0.0,Z2OS1/VOLWMM)
          CNH3S1=AMAX1(0.0,ZH3S1/VOLWMM)
          CNH3B1=AMAX1(0.0,ZH3B1/VOLWMM)
          CH2GS1=AMAX1(0.0,H2GS1/VOLWMM)
          IF(RTVLPs1(N,L,NZ).GT.ZEROs1)THEN
            CCO2A1=AMAX1(0.0,CO2A1/RTVLPs1(N,L,NZ))
            COXYA1=AMAX1(0.0,OXYA1/RTVLPs1(N,L,NZ))
            CCH4A1=AMAX1(0.0,CH4A1/RTVLPs1(N,L,NZ))
            CZ2OA1=AMAX1(0.0,Z2OA1/RTVLPs1(N,L,NZ))
            CNH3A1=AMAX1(0.0,ZH3A1/RTVLPs1(N,L,NZ))
            CH2GA1=AMAX1(0.0,H2GA1/RTVLPs1(N,L,NZ))
          ELSE
            CCO2A1=0.0_r8
            COXYA1=0.0_r8
            CCH4A1=0.0_r8
            CZ2OA1=0.0_r8
            CNH3A1=0.0_r8
            CH2GA1=0.0_r8
          ENDIF
          CCO2P1=AMAX1(0.0,CO2P1/RTVLWs1(N,L,NZ))
          COXYP1=AMIN1(COXYEs1*SOXYLs1(L),AMAX1(0.0,OXYP1/RTVLWs1(N,L,NZ)))
          CCH4P1=AMAX1(0.0,CH4P1/RTVLWs1(N,L,NZ))
          CN2OP1=AMAX1(0.0,Z2OP1/RTVLWs1(N,L,NZ))
          CNH3P1=AMAX1(0.0,ZH3P1/RTVLWs1(N,L,NZ))
          CH2GP1=AMAX1(0.0,H2GP1/RTVLWs1(N,L,NZ))
          DIFOX=DIFOL+DIFOP
          RMFCOS=UPWTRH*CCO2S1
          RMFOXS=UPWTRH*COXYS1
          RMFCHS=UPWTRH*CCH4S1
          RMFN2S=UPWTRH*CN2OS1
          RMFN3S=UPWTRH*CNH3S1*VLNH4s1(L)
          RMFN3B=UPWTRH*CNH3B1*VLNHBs1(L)
          RMFHGS=UPWTRH*CH2GS1
!
!     SOLUTION FOR MASS FLOW + DIFFUSION OF O2 IN AQUEOUS PHASES OF
!     SOIL AND ROOT = ACTIVE UPTAKE OF O2 BY ROOT
!
!     DIFOL=O2 aqueous diffusivity from soil to root
!     UPWTRH=water uptake
!     DIFOP=aqueous diffusivity of O2 within root
!     COXYS1,COXYP1=soil,root aqueous O2 concentration
!     UPMXP=O2 demand per plant
!     RUPOXR=root O2 uptake per plant
!     COXYR=aqueous O2 concentration at root surface
!     RDFOXS,RDFOXP=aqueous O2 diffusion per plant:soil-root,within root
!
          X=(DIFOL+UPWTRH)*COXYS1+DIFOP*COXYP1
          IF(X.GT.ZEROs1.AND.OXYS1.GT.ZEROPs1(NZ))THEN
            B=-UPMXP-DIFOX*OXKM-X
            C=X*UPMXP
            RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
            COXYR=(X-RUPOXR)/DIFOX
            RDFOXS=RMFOXS+DIFOL*(COXYS1-COXYR)
            RDFOXP=DIFOP*(COXYP1-COXYR)
          ELSE
            X=DIFOP*COXYP1
            IF(X.GT.ZEROs1.AND.OXYP1.GT.ZEROPs1(NZ))THEN
              B=-UPMXP-DIFOP*OXKM-X
              C=X*UPMXP
              RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
              COXYR=(X-RUPOXR)/DIFOP
              RDFOXS=0.0_r8
              RDFOXP=DIFOP*(COXYP1-COXYR)
            ELSE
              RUPOXR=0.0_r8
              COXYR=0.0_r8
              RDFOXS=0.0_r8
              RDFOXP=0.0_r8
            ENDIF
          ENDIF
!
!     MASS FLOW + DIFFUSIVE EXCHANGE OF OTHER GASES
!     BETWEEN ROOT AND SOIL, CONSTRAINED BY COMPETITION
!     WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!     RUPOSX,RUPOPX=aqueous O2 uptake from soil,root
!     PP=PFT population
!     RDFCOS,RCO2SX=aqueous CO2 soil-root diffusion,root uptake
!     RDFCHS,RUPCSX=aqueous CH4 soil-root diffusion,root uptake
!     RDFN2S,RUPZSX=aqueous N2O soil-root diffusion,root uptake
!     RDFN3S,RUPNSX=aqueous NH3 soil-root diffusion,root uptake:non-band
!     RDFN3B,RUPNBX=aqueous NH3 soil-root diffusion,root uptake:band
!     RDFHGS,RUPHGX=aqueous H2 soil-root diffusion,root uptake
!     RMF*=soil convective solute flux
!     DIF*=aqueous diffusivity from soil to root
!     C*S1=soil aqueous concentration non-band
!     C*B1=soil aqueous concentration band
!     C*P1=root aqueous concentration
!

          RUPOSX=RDFOXS*PPs1(NZ)
          RUPOPX=RDFOXP*PPs1(NZ)
          RDFCOS=RMFCOS+DIFCL*(CCO2S1-CCO2P1)
          RDXCOS=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),CO2S1) &
            -VOLWMM*AMAX1(ZEROPs1(NZ),CO2P1))/VOLWSP
          IF(RDFCOS.GT.0.0)THEN
            RCO2SX=AMIN1(AMAX1(0.0,RDXCOS),RDFCOS*PPs1(NZ))
          ELSE
            RCO2SX=AMAX1(AMIN1(0.0,RDXCOS),RDFCOS*PPs1(NZ))
          ENDIF
          IF(N.EQ.1)THEN
            RDFCHS=RMFCHS+DIFCL*(CCH4S1-CCH4P1)
            RDXCHS=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),CH4S1) &
              -VOLWMM*AMAX1(ZEROPs1(NZ),CH4P1))/VOLWSP
            IF(RDFCHS.GT.0.0)THEN
              RUPCSX=AMIN1(AMAX1(0.0,RDXCHS),RDFCHS*PPs1(NZ))
            ELSE
              RUPCSX=AMAX1(AMIN1(0.0,RDXCHS),RDFCHS*PPs1(NZ))
            ENDIF
            RDFN2S=RMFN2S+DIFZL*(CN2OS1-CN2OP1)
            RDXN2S=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),Z2OS1) &
              -VOLWMM*AMAX1(ZEROPs1(NZ),Z2OP1))/VOLWSP
            IF(RDFN2S.GT.0.0)THEN
              RUPZSX=AMIN1(AMAX1(0.0,RDXN2S),RDFN2S*PPs1(NZ))
            ELSE
              RUPZSX=AMAX1(AMIN1(0.0,RDXN2S),RDFN2S*PPs1(NZ))
            ENDIF
            RDFN3S=RMFN3S+DIFNL*(CNH3S1-CNH3P1)
            IF(VOLWSA.GT.ZEROPs1(NZ))THEN
              ZH3PA=ZH3P1*VLNH4s1(L)
              RDXNHS=(RTVLWA*AMAX1(ZEROPs1(NZ),ZH3S1) &
                -VOLWMA*AMAX1(ZEROPs1(NZ),ZH3PA))/VOLWSA
            ELSE
              RDXNHS=0.0_r8
            ENDIF
            IF(RDFN3S.GT.0.0)THEN
              RUPNSX=AMIN1(AMAX1(0.0,RDXNHS),RDFN3S*PPs1(NZ))
            ELSE
              RUPNSX=AMAX1(AMIN1(0.0,RDXNHS),RDFN3S*PPs1(NZ))
            ENDIF
            RDFN3B=RMFN3B+DIFNB*(CNH3B1-CNH3P1)
            IF(VOLWSB.GT.ZEROPs1(NZ))THEN
              ZH3PB=ZH3P1*VLNHBs1(L)
              RDXNHB=(RTVLWB*AMAX1(ZEROPs1(NZ),ZH3B1) &
                -VOLWMB*AMAX1(ZEROPs1(NZ),ZH3PB))/VOLWSB
            ELSE
              RDXNHB=0.0_r8
            ENDIF
            IF(RDFN3B.GT.0.0)THEN
              RUPNBX=AMIN1(AMAX1(0.0,RDXNHB),RDFN3B*PPs1(NZ))
            ELSE
              RUPNBX=AMAX1(AMIN1(0.0,RDXNHB),RDFN3B*PPs1(NZ))
            ENDIF
            RDFHGS=RMFHGS+DIFHL*(CH2GS1-CH2GP1)
            RDXHGS=(RTVLWs1(N,L,NZ)*AMAX1(ZEROPs1(NZ),H2GS1) &
              -VOLWMM*AMAX1(ZEROPs1(NZ),H2GP1))/VOLWSP
            IF(RDFHGS.GT.0.0)THEN
              RUPHGX=AMIN1(AMAX1(0.0,RDXHGS),RDFHGS*PPs1(NZ))
            ELSE
              RUPHGX=AMAX1(AMIN1(0.0,RDXHGS),RDFHGS*PPs1(NZ))
            ENDIF
          ELSE
            RUPCSX=0.0_r8
            RUPZSX=0.0_r8
            RUPNSX=0.0_r8
            RUPNBX=0.0_r8
            RUPHGX=0.0_r8
          ENDIF
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN SOIL
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENTS
!     FROM 'WATSUB'
!
!     THETPM,THETX=air-filled porosity,minimum THETPM
!     R*DFQ=soil gas exchange between gaseous-aqueous phases
!     DFGS=rate constant for soil gas exchange from watsub.f
!     CO2G1,CO2S1=gaseous,aqueous CO2 in soil
!     OXYG1,OXYS1=gaseous,aqueous O2 in soil
!     CH4G1,CH4S1=gaseous,aqueous CH4 in soil
!     Z2OG1,Z2OS1=gaseous,aqueous N2O in soil
!     ZH3G1,ZH3S1,ZH3B1=gaseous,aqueous NH3 in soil non-band,band
!     H2GG1,H2GS1=gaseous,aqueous H2 in soil
!     RUPOSX=root aqueous O2 uptake
!     ROXYLX=root net O2 aqueous flux
!     RCO2SX=root aqueous CO2 uptake
!     RUPCSX=root aqueous CH4 uptake
!     RUPZSX=root aqueous N2O uptake
!     RUPNSX=root aqueous NH3 uptake non-band
!     RUPNBX=root aqueous NH3 uptake band
!     RUPHGX=root aqueous H2 uptake
!     VOLWMM,VOLPMM=soil micropore water,air volume
!     VOLW*=VOLWMM*gas solubility
!
          IF(THETPMs1(M,L).GT.THETX)THEN
            DFGSP=FPQ(N,L,NZ)*DFGSs1(M,L)
            RCODFQ=DFGSP*(AMAX1(ZEROPs1(NZ),CO2G1)*VOLWCO &
              -(AMAX1(ZEROSs1,CO2S1)-RCO2SX)*VOLPMM)/(VOLWCO+VOLPMM)
            RUPOST=RUPOSX-ROXYLX
            ROXDFQ=DFGSP*(AMAX1(ZEROPs1(NZ),OXYG1)*VOLWOX &
              -(AMAX1(ZEROSs1,OXYS1)-RUPOST)*VOLPMM)/(VOLWOX+VOLPMM)
            IF(N.EQ.1)THEN
              RCHDFQ=DFGSP*(AMAX1(ZEROPs1(NZ),CH4G1)*VOLWCH &
                -(AMAX1(ZEROSs1,CH4S1)-RUPCSX)*VOLPMM)/(VOLWCH+VOLPMM)
              RN2DFQ=DFGSP*(AMAX1(ZEROPs1(NZ),Z2OG1)*VOLWN2 &
                -(AMAX1(ZEROSs1,Z2OS1)-RUPZSX)*VOLPMM)/(VOLWN2+VOLPMM)
              IF(VOLWNH+VOLPNH.GT.ZEROPs1(NZ))THEN
                ZH3GA=ZH3G1*VLNH4s1(L)
                RNHDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROPs1(NZ),ZH3GA)*VOLWNH &
                  -(AMAX1(ZEROSs1,ZH3S1)-RUPNSX)*VOLPNH)/(VOLWNH+VOLPNH)))
              ELSE
                RNHDFQ=0.0_r8
              ENDIF
              IF(VOLWNB+VOLPNB.GT.ZEROPs1(NZ))THEN
                ZH3GB=ZH3G1*VLNHBs1(L)
                RNBDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROPs1(NZ),ZH3GB)*VOLWNB &
                  -(AMAX1(ZEROSs1,ZH3B1)-RUPNBX)*VOLPNB)/(VOLWNB+VOLPNB)))
              ELSE
                RNBDFQ=0.0_r8
              ENDIF
              RHGDFQ=DFGSP*(AMAX1(ZEROPs1(NZ),H2GG1)*VOLWHG &
                -(AMAX1(ZEROSs1,H2GS1)-RUPHGX)*VOLPMM)/(VOLWHG+VOLPMM)
            ELSE
              RCHDFQ=0.0_r8
              RN2DFQ=0.0_r8
              RNHDFQ=0.0_r8
              RNBDFQ=0.0_r8
              RHGDFQ=0.0_r8
            ENDIF
          ELSE
            RCODFQ=0.0_r8
            ROXDFQ=0.0_r8
            RCHDFQ=0.0_r8
            RN2DFQ=0.0_r8
            RNHDFQ=0.0_r8
            RNBDFQ=0.0_r8
            RHGDFQ=0.0_r8
          ENDIF
!
!     UPDATE GASEOUS, AQUEOUS GAS CONTENTS AND CONCENTRATIONS
!     FROM GASEOUS-AQUEOUS EXCHANGE, SOIL GAS TRANSFERS
!
          OXYG1=OXYG1-ROXDFQ+ROXYFX
          OXYS1=OXYS1+ROXDFQ-RUPOSX
          CO2G1=CO2G1-RCODFQ+RCO2FX
          CO2S1=CO2S1+RCODFQ-RCO2SX
          CH4S1=CH4S1+RCHDFQ-RUPCSX
          Z2OS1=Z2OS1+RN2DFQ-RUPZSX
          ZH3S1=ZH3S1+RNHDFQ-RUPNSX
          ZH3B1=ZH3B1+RNBDFQ-RUPNBX
          H2GS1=H2GS1+RHGDFQ-RUPHGX
!
!     GAS TRANSFER THROUGH ROOTS
!
          IF(N.EQ.1.AND.RTVLPs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
            RUPNTX=RUPNSX+RUPNBX
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
!
!     R*DF1=root gas exchange between gaseous-aqueous phases
!     R*FL1=root gas exchange with atmosphere
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     CO2A1,CO2P1=gaseous,aqueous CO2 in root
!     OXYA1,OXYP1=gaseous,aqueous O2 in root
!     CH4A1,CH4P1=gaseous,aqueous CH4 in root
!     Z2OA1,Z2OP1=gaseous,aqueous N2O in root
!     ZH3A1,ZH3P1=gaseous,aqueous NH3 in root
!     H2GA1,H2GP1=gaseous,aqueous H2 in root
!     RTVLW,RTVLP=root aqueous,gaseous volume
!     VOLW*=RTVLW*gas solubility
!     C*E,C*A1=atmosphere,root gas concentration
!     DF*A=root-atmosphere gas conductance
!
            CO2PX=CO2P1+RCO2PX
            RCODF1=AMAX1(-CO2PX,DFGP*(AMAX1(ZEROPs1(NZ),CO2A1)*VOLWCA &
              -CO2PX*RTVLPs1(N,L,NZ))/(VOLWCA+RTVLPs1(N,L,NZ)))
            OXYPX=OXYP1-RUPOPX
            ROXDF1=AMAX1(-OXYPX,DFGP*(AMAX1(ZEROPs1(NZ),OXYA1)*VOLWOA &
              -OXYPX*RTVLPs1(N,L,NZ))/(VOLWOA+RTVLPs1(N,L,NZ)))
            CH4PX=CH4P1+RUPCSX
            RCHDF1=AMAX1(-CH4PX,DFGP*(AMAX1(ZEROPs1(NZ),CH4A1)*VOLWC4 &
              -CH4PX*RTVLPs1(N,L,NZ))/(VOLWC4+RTVLPs1(N,L,NZ)))
            Z2OPX=Z2OP1+RUPZSX
            RN2DF1=AMAX1(-Z2OPX,DFGP*(AMAX1(ZEROPs1(NZ),Z2OA1)*VOLWZA &
              -Z2OPX*RTVLPs1(N,L,NZ))/(VOLWZA+RTVLPs1(N,L,NZ)))
            ZH3PX=ZH3P1+RUPNTX
            RNHDF1=AMAX1(-ZH3PX,DFGP*(AMAX1(ZEROPs1(NZ),ZH3A1)*VOLWNA &
              -ZH3PX*RTVLPs1(N,L,NZ))/(VOLWNA+RTVLPs1(N,L,NZ)))
            H2GPX=H2GP1+RUPHGX
            RHGDF1=AMAX1(-H2GPX,DFGP*(AMAX1(ZEROPs1(NZ),H2GA1)*VOLWH2 &
              -H2GPX*RTVLPs1(N,L,NZ))/(VOLWH2+RTVLPs1(N,L,NZ)))
            RCOFL1=AMIN1(DFCOA,RTVLPs1(N,L,NZ))*(CCO2Es1-CCO2A1)
            ROXFL1=AMIN1(DFOXA,RTVLPs1(N,L,NZ))*(COXYEs1-COXYA1)
            RCHFL1=AMIN1(DFCHA,RTVLPs1(N,L,NZ))*(CCH4Es1-CCH4A1)
            RN2FL1=AMIN1(DFN2A,RTVLPs1(N,L,NZ))*(CZ2OEs1-CZ2OA1)
            RNHFL1=AMIN1(DFNHA,RTVLPs1(N,L,NZ))*(CNH3Es1-CNH3A1)
            RHGFL1=AMIN1(DFHGA,RTVLPs1(N,L,NZ))*(CH2GEs1-CH2GA1)
          ELSE
            RCODF1=0.0_r8
            ROXDF1=0.0_r8
            RCHDF1=0.0_r8
            RN2DF1=0.0_r8
            RNHDF1=0.0_r8
            RHGDF1=0.0_r8
            RCOFL1=0.0_r8
            ROXFL1=0.0_r8
            RCHFL1=0.0_r8
            RN2FL1=0.0_r8
            RNHFL1=0.0_r8
            RHGFL1=0.0_r8
          ENDIF
!
!     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS
!     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
!
          CO2A1=CO2A1-RCODF1+RCOFL1
          OXYA1=OXYA1-ROXDF1+ROXFL1
          CH4A1=CH4A1-RCHDF1+RCHFL1
          Z2OA1=Z2OA1-RN2DF1+RN2FL1
          ZH3A1=ZH3A1-RNHDF1+RNHFL1
          H2GA1=H2GA1-RHGDF1+RHGFL1
          CO2P1=CO2P1+RCODF1+RCO2SX+RCO2PX
          OXYP1=OXYP1+ROXDF1-RUPOPX
          CH4P1=CH4P1+RCHDF1+RUPCSX
          Z2OP1=Z2OP1+RN2DF1+RUPZSX
          ZH3P1=ZH3P1+RNHDF1+RUPNSX+RUPNBX
          H2GP1=H2GP1+RHGDF1+RUPHGX
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2S=soil-root CO2 exchange
!     RUPOXS=soil-root O2 exchange
!     RUPCHS=soil-root CH4 exchange
!     RUPN2S=soil-root N2O exchange
!     RUPN3S=soil-root NH3 exchange non-band
!     RUPN3B=soil-root NH3 exchange band
!     RUPHGS=soil-root H2 exchange
!
          RCO2Ss1(N,L,NZ)=RCO2Ss1(N,L,NZ)+RCO2SX
          RUPOXSs1(N,L,NZ)=RUPOXSs1(N,L,NZ)+RUPOSX
          RUPCHSs1(N,L,NZ)=RUPCHSs1(N,L,NZ)+RUPCSX
          RUPN2Ss1(N,L,NZ)=RUPN2Ss1(N,L,NZ)+RUPZSX
          RUPN3Ss1(N,L,NZ)=RUPN3Ss1(N,L,NZ)+RUPNSX
          RUPN3Bs1(N,L,NZ)=RUPN3Bs1(N,L,NZ)+RUPNBX
          RUPHGSs1(N,L,NZ)=RUPHGSs1(N,L,NZ)+RUPHGX
!
!     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE
!
!     R*DFA=root aqueous-gaseous CO2 exchange
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!
          RCODFAs1(N,L,NZ)=RCODFAs1(N,L,NZ)+RCODF1
          ROXDFAs1(N,L,NZ)=ROXDFAs1(N,L,NZ)+ROXDF1
          RCHDFAs1(N,L,NZ)=RCHDFAs1(N,L,NZ)+RCHDF1
          RN2DFAs1(N,L,NZ)=RN2DFAs1(N,L,NZ)+RN2DF1
          RNHDFAs1(N,L,NZ)=RNHDFAs1(N,L,NZ)+RNHDF1
          RHGDFAs1(N,L,NZ)=RHGDFAs1(N,L,NZ)+RHGDF1
          RCOFLAs1(N,L,NZ)=RCOFLAs1(N,L,NZ)+RCOFL1
          ROXFLAs1(N,L,NZ)=ROXFLAs1(N,L,NZ)+ROXFL1
          RCHFLAs1(N,L,NZ)=RCHFLAs1(N,L,NZ)+RCHFL1
          RN2FLAs1(N,L,NZ)=RN2FLAs1(N,L,NZ)+RN2FL1
          RNHFLAs1(N,L,NZ)=RNHFLAs1(N,L,NZ)+RNHFL1
          RHGFLAs1(N,L,NZ)=RHGFLAs1(N,L,NZ)+RHGFL1
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2P=root CO2 emission into root
!     RUPOXP=root O2 uptake from root
!     ROXSK=total O2 uptake from soil by all microbial,root popns
!
          RCO2Ps1(N,L,NZ)=RCO2Ps1(N,L,NZ)+RCO2PX+RCO2SX
          RUPOXPs1(N,L,NZ)=RUPOXPs1(N,L,NZ)+RUPOPX
          ROXSKs1(M,L)=ROXSKs1(M,L)+RUPOSX

90      CONTINUE
      ENDIF
99  CONTINUE
!
!     O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO
!     OF ROOT O2 UPTAKE 'RUPOXT' TO ROOT O2 DEMAND 'ROXYP'
!
!     RUPOXT=O2 uptake from soil+root by each root,myco population
!     ROXYP=O2 demand by each root,myco population
!     WFR=constraint by O2 consumption on all root processes
!     imposed by O2 uptake
!
    RUPOXT=RUPOXPs1(N,L,NZ)+RUPOXSs1(N,L,NZ)
    WFRs1(N,L,NZ)=AMIN1(1.0,AMAX1(0.0 &
      ,RUPOXT/ROXYPs1(N,L,NZ)))
  ELSE
    RUPOXT=0.0_r8
    IF(L.GT.NGs1(NZ))THEN
      WFRs1(N,L,NZ)=WFRs1(N,L-1,NZ)
    ELSE
      WFRs1(N,L,NZ)=1.0
    ENDIF
  ENDIF
  end associate
  end subroutine RootSoilGasExchange

end module RootGasMod
