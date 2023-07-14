module RootGasMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use StomatesMod   , only : stomates
  use minimathmod  , only : safe_adb,vapsat,AZMAX1,AZMIN1
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use PlantAPIData
  use TracerIDMod
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = __FILE__
  public :: RootSoilGasExchange
  contains
!------------------------------------------------------------------------

  subroutine RootSoilGasExchange(N,L,NZ,RRADL,FPQ,FRTDPX,RTARR,&
    UPWTRH,FOXYX,PopPlantO2Uptake_vr)

  implicit none
  integer , intent(in) :: N,L,NZ
  real(r8), intent(in) :: RRADL(2,JZ1),FPQ(2,JZ1,JP1),FRTDPX(JZ1,JP1)
  real(r8), intent(in) :: RTARR(2,JZ1),UPWTRH,FOXYX
  real(r8), intent(out):: PopPlantO2Uptake_vr
  integer :: M,MX
  real(r8) :: B,C
  real(r8) :: trcg_gmas(idg_beg:idg_end-1)
  real(r8) :: CO2P1,CO2G1,CO2S1
  real(r8) :: CH4P1,CH4S1,CCH4S1,CCH4P1
  real(r8) :: CN2OS1,CN2OP1,CNH3S1,CNH3B1,CNH3P1,CH2GS1,CH2GP1
  real(r8) :: CGSGL1,CHSGL1,CLSGL1
  real(r8) :: CQSGL1,CH4G1,CCO2S1,COXYS1,trcg_gcon(idg_beg:idg_end-1)
  real(r8) :: CCO2P1,COXYP1,COXYR,CO2PX,CH4PX
  real(r8) :: DIFOP,DIFCL,DIFZL,DIFNL,DIFNB,DIFHL,DIFOX,DFGSP
  real(r8) :: DFCOA,DFOXA,DFCHA,DFN2A,DFNHA,DFHGA,DFGP,DIFOL
  real(r8) :: H2GP1,H2GS1,HGSGL1,HLSGL1,H2GG1,H2GPX
  real(r8) :: OXYP1,OXYG1,OXYS1,OGSGL1
  real(r8) :: OLSGL1,OLSGLP,OXYPX,RTVLWA
  real(r8) :: RTVLWB,ROXYFX,RCO2FX,ROXYLX,ROXDFQ,RCHDFQ,RN2DFQ
  real(r8) :: RNHDFQ,RHGDFQ
  real(r8) :: RTCR1,RTCR2,RTCRA,RTARRX,RCO2PX,RRADS,RMFCOS,RMFOXS
  real(r8) :: RMFCHS,RMFN2S,RMFN3S,RMFN3B,RMFHGS,RUPOXR,RDFOXS
  real(r8) :: RDFOXP,RUPOSX,RUPOPX,RDFCOS,RDXCOS,RCO2SX,RDFCHS
  real(r8) :: RDXCHS,RUPCSX,RDFN2S,RDXN2S,RUPZSX,RDFN3S,RDXNHS
  real(r8) :: RUPNSX,RDFN3B,RDXNHB,RUPNBX,RDFHGS,RDXHGS,RUPHGX
  real(r8) :: RCODFQ,RUPOST,RNBDFQ,RUPNTX
  real(r8) :: trcg_RDF1(idg_beg:idg_end-1),trcg_RFL1(idg_beg:idg_end-1)
  real(r8) :: THETW1,THETM
  real(r8) :: UPMXP
  real(r8) :: VOLWG(idg_beg:idg_end-1),VOLWMO,VOLWMM,VOLPMM
  real(r8) :: VOLWSP,VOLWMA,VOLWMB,VOLWSA,VOLWSB,VOLWCO,VOLWOX
  real(r8) :: VOLWCH,VOLWN2,VOLWNH,VOLWNB,VOLWHG,VOLPNH,VOLPNB
  real(r8) :: X
  real(r8) :: Z2OP1,Z2OS1,ZH3P1,ZH3S1
  real(r8) :: ZH3B1,Z2SGL1,ZHSGL1,ZVSGL1,ZNSGL1,Z2OG1,ZH3G1
  real(r8) :: ZH3PA,ZH3PB,ZH3GA,ZH3GB,Z2OPX,ZH3PX
  integer  :: NTG
!     begin_execution
  associate(                       &
    WTRTSE =>  plt_biom%WTRTSE   , &
    ZEROP  =>  plt_biom%ZEROP    , &
    pftPlantPopulation     =>  plt_site%pftPlantPopulation       , &
    DPTHZ  =>  plt_site%DPTHZ    , &
    CCH4E  =>  plt_site%CCH4E    , &
    CCO2E  =>  plt_site%CCO2E    , &
    CH2GE  =>  plt_site%CH2GE    , &
    CNH3E  =>  plt_site%CNH3E    , &
    COXYE  =>  plt_site%COXYE    , &
    CZ2OE  =>  plt_site%CZ2OE    , &
    ZEROS  =>  plt_site%ZEROS    , &
    ZERO   =>  plt_site%ZERO     , &
    VOLWM  =>  plt_site%VOLWM    , &
    VOLPM  =>  plt_site%VOLPM    , &
    TORT   =>  plt_site%TORT     , &
    FILM   =>  plt_site%FILM     , &
    ROXYF  =>  plt_bgcr%ROXYF    , &
    RCO2F  =>  plt_bgcr%RCO2F    , &
    ROXYL  =>  plt_bgcr%ROXYL    , &
    WFR    =>  plt_rbgc%WFR      , &
    ZEROQ  =>  plt_rbgc%ZEROQ    , &
    RCO2M  =>  plt_rbgc%RCO2M    , &
    ROXYP  =>  plt_rbgc%ROXYP    , &
    RUPOXS =>  plt_rbgc%RUPOXS   , &
    ROXSK  =>  plt_rbgc%ROXSK    , &
    RCO2P  =>  plt_rbgc%RCO2P    , &
    trcg_RFLA =>  plt_rbgc%trcg_RFLA   , &
    trcg_RDFA =>  plt_rbgc%trcg_RDFA   , &
    RUPN3B =>  plt_rbgc%RUPN3B   , &
    RUPHGS =>  plt_rbgc%RUPHGS   , &
    RCO2S  =>  plt_rbgc%RCO2S    , &
    RUPCHS =>  plt_rbgc%RUPCHS   , &
    RUPN2S =>  plt_rbgc%RUPN2S   , &
    RUPN3S =>  plt_rbgc%RUPN3S   , &
    RUPOXP =>  plt_rbgc%RUPOXP   , &
    RCO2A  =>  plt_rbgc%RCO2A    , &
    trcg_rootml   =>  plt_rbgc%trcg_rootml     , &
    trcs_rootml   =>  plt_rbgc%trcs_rootml     , &
    TFND   =>  plt_soilchem%TFND , &
    trcs_VLN  =>  plt_soilchem%trcs_VLN, &
    trc_solml  =>  plt_soilchem%trc_solml, &
    trc_gascl  =>  plt_soilchem%trc_gascl, &
    GasDifc=>  plt_soilchem%GasDifc,&
    SolDifc  =>  plt_soilchem%SolDifc, &
    GSolbility=> plt_soilchem%GSolbility,&
    trc_solcl  =>  plt_soilchem%trc_solcl, &
    THETY  =>  plt_soilchem%THETY, &
    VOLY   =>  plt_soilchem%VOLY , &
    THETPM =>  plt_soilchem%THETPM,&
    trc_gasml=> plt_soilchem%trc_gasml,&
    DFGS   =>  plt_soilchem%DFGS , &
    IDAY   =>  plt_pheno%IDAY    , &
    PORTX  =>  plt_morph%PORTX   , &
    RRAD1  =>  plt_morph%RRAD1   , &
    RTLGA  =>  plt_morph%RTLGA   , &
    RTVLP  =>  plt_morph%RTVLP   , &
    RTLGP  =>  plt_morph%RTLGP   , &
    RTNL   =>  plt_morph%RTNL    , &
    RRAD2  =>  plt_morph%RRAD2   , &
    RRADP  =>  plt_morph%RRADP   , &
    RTVLW  =>  plt_morph%RTVLW   , &
    PORT   =>  plt_morph%PORT    , &
    RTN1   =>  plt_morph%RTN1    , &
    NG     =>  plt_morph%NG      , &
    NB1    =>  plt_morph%NB1       &
  )
  IF(RCO2M(N,L,NZ).GT.ZEROP(NZ).AND.RTVLW(N,L,NZ).GT.ZEROP(NZ) &
    .AND.FOXYX.GT.ZEROQ(NZ))THEN
!
!     INITIALIZE VARIABLES USED IN ROOT GAS EXCHANGE
!     (CO2, O2, CH4, N2, N2O, NH3, H2)
!
!     trcg_gmas(idg_CO2),CO2P1,CO2G1,CO2S1=gaseous,aqueous CO2 in root,soil
!     trcg_gmas(idg_O2),OXYP1,OXYG1,OXYS1=gaseous,aqueous O2 in root,soil
!     trcg_gmas(idg_CH4),CH4P1,CH4G1,CH4S1=gaseous,aqueous CH4 in root,soil
!     trcg_gmas(idg_N2O),Z2OP1,Z2OG1,Z2OS1=gaseous,aqueous N2O in root,soil
!     trcg_gmas(idg_NH3),ZH3P1,ZH3G1,ZH3S1=gaseous,aqueous NH3 in root,soil
!     trcg_gmas(idg_H2),H2GP1,H2GG1,H2GS1=gaseous,aqueous H2 in root,soil
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
    trcg_gmas(idg_CO2)=AMAX1(ZEROP(NZ),trcg_rootml(idg_CO2,N,L,NZ))
    CO2P1=AMAX1(ZEROP(NZ),trcs_rootml(idg_CO2,N,L,NZ))
    CO2G1=AMAX1(ZEROP(NZ),trc_gasml(idg_CO2,L)*FPQ(N,L,NZ))
    CO2S1=AMAX1(ZEROP(NZ),trc_solml(idg_CO2,L)*FPQ(N,L,NZ))

    trcg_gmas(idg_O2)=AMAX1(ZEROP(NZ),trcg_rootml(idg_O2,N,L,NZ))
    OXYP1=AMAX1(ZEROP(NZ),trcs_rootml(idg_O2,N,L,NZ))
    OXYG1=AMAX1(ZEROP(NZ),trc_gasml(idg_O2,L)*FOXYX)
    OXYS1=trc_solml(idg_O2,L)*FOXYX

    trcg_gmas(idg_CH4)=trcg_rootml(idg_CH4,N,L,NZ)
    CH4P1=trcs_rootml(idg_CH4,N,L,NZ)
    CH4S1=trc_solml(idg_CH4,L)*FPQ(N,L,NZ)
    CCH4S1=trc_solcl(idg_CH4,L)
    CCH4P1=AZMAX1(CH4P1/RTVLW(N,L,NZ))

    trcg_gmas(idg_N2O)=trcg_rootml(idg_N2O,N,L,NZ)
    Z2OP1=trcs_rootml(idg_N2O,N,L,NZ)
    Z2OS1=trc_solml(idg_N2O,L)*FPQ(N,L,NZ)
    CN2OS1=trc_solcl(idg_N2O,L)
    CN2OP1=AZMAX1(Z2OP1/RTVLW(N,L,NZ))

    trcg_gmas(idg_NH3)=trcg_rootml(idg_NH3,N,L,NZ)
    ZH3P1=trcs_rootml(idg_NH3,N,L,NZ)
    ZH3S1=trc_solml(idg_NH3,L)*FPQ(N,L,NZ)
    ZH3B1=trc_solml(idg_NH3B,L)*FPQ(N,L,NZ)
    CNH3S1=trc_solcl(idg_NH3,L)
    CNH3B1=trc_solcl(idg_NH3B,L)
    CNH3P1=AZMAX1(ZH3P1/RTVLW(N,L,NZ))

    trcg_gmas(idg_H2)=trcg_rootml(idg_H2,N,L,NZ)
    H2GP1=trcs_rootml(idg_H2,N,L,NZ)
    H2GS1=trc_solml(idg_H2,L)*FPQ(N,L,NZ)
    CH2GS1=trc_solcl(idg_H2,L)
    CH2GP1=AZMAX1(H2GP1/RTVLW(N,L,NZ))

    RTVLWA=RTVLW(N,L,NZ)*trcs_VLN(ids_NH4,L)
    RTVLWB=RTVLW(N,L,NZ)*trcs_VLN(ids_NH4B,L)
    UPMXP=ROXYP(N,L,NZ)*XNPG/pftPlantPopulation(NZ)
    ROXYFX=ROXYF(L)*FOXYX*XNPG
    RCO2FX=RCO2F(L)*FOXYX*XNPG
    ROXYLX=ROXYL(L)*FOXYX*XNPG
!
!     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
!
!     *SGL1=diffusivity
!     PORTX=tortuosity effect of root porosity on diffusivity
!     CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
!     CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
!
    CGSGL1=GasDifc(idg_CO2,L)*XNPG*PORTX(N,NZ)
    OGSGL1=GasDifc(idg_O2,L)*XNPG*PORTX(N,NZ)
    CHSGL1=GasDifc(idg_CH4,L)*XNPG*PORTX(N,NZ)
    Z2SGL1=GasDifc(idg_N2O,L)*XNPG*PORTX(N,NZ)
    ZHSGL1=GasDifc(idg_NH3,L)*XNPG*PORTX(N,NZ)
    HGSGL1=GasDifc(idg_H2,L)*XNPG*PORTX(N,NZ)

    CLSGL1=SolDifc(idg_CO2,L)*XNPG*FOXYX
    OLSGL1=SolDifc(idg_O2,L)*XNPG*FOXYX
    CQSGL1=SolDifc(idg_CH4,L)*XNPG*FOXYX
    ZVSGL1=SolDifc(idg_N2O,L)*XNPG*FOXYX
    ZNSGL1=SolDifc(idg_NH3,L)*XNPG*FOXYX
    HLSGL1=SolDifc(idg_H2,L)*XNPG*FOXYX
    OLSGLP=SolDifc(idg_O2,L)*XNPG

    ROXDFQ=0.0_r8
    RCHDFQ=0.0_r8
    RN2DFQ=0.0_r8
    RNHDFQ=0.0_r8
    RHGDFQ=0.0_r8
    trcg_RDF1(idg_beg:idg_end-1)=0.0_r8
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
    IF(WTRTSE(ielmc,NZ).GT.ZEROP(NZ).AND.FRTDPX(L,NZ).GT.ZERO)THEN
      RTCR1=AMAX1(pftPlantPopulation(NZ),RTN1(N,L,NZ)) &
        *PICON*RRAD1(N,L,NZ)**2/DPTHZ(L)
      RTCR2=(RTNL(N,L,NZ)*PICON*RRAD2(N,L,NZ)**2 &
        /RTLGA(N,L,NZ))/FRTDPX(L,NZ)
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
!     IDAY(1,=emergence date
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
    IF(N.EQ.1.AND.IDAY(1,NB1(NZ),NZ).GT.0.AND.RTLGP(N,L,NZ).GT.ZEROP(NZ))THEN
      RTARRX=RTARR(N,L)/RRADP(N,NZ)
      DIFOP=OLSGLP*RTARRX

      DO NTG=idg_beg,idg_end-1
        VOLWG(NTG)=RTVLW(N,L,NZ)*GSolbility(NTG,L)
      ENDDO

      DFCOA=CGSGL1*RTCRA
      DFOXA=OGSGL1*RTCRA
      DFCHA=CHSGL1*RTCRA
      DFN2A=Z2SGL1*RTCRA
      DFNHA=ZHSGL1*RTCRA
      DFHGA=HGSGL1*RTCRA
    ELSE
      RTARRX=0.0_r8
      DIFOP=0.0_r8
      VOLWG(idg_beg:idg_end-1)=0.0_r8
      DFCOA=0.0_r8
      DFOXA=0.0_r8
      DFCHA=0.0_r8
      DFN2A=0.0_r8
      DFNHA=0.0_r8
      DFHGA=0.0_r8
    ENDIF
    DFGP=AMIN1(1.0,XNPD*SQRT(PORT(N,NZ))*TFND(L))
    RCO2PX=-RCO2A(N,L,NZ)*XNPG
!
!     SOLVE FOR GAS EXCHANGE IN SOIL AND ROOTS DURING ROOT UPTAKE
!     AT SMALLER TIME STEP NPH
!
    D99: DO M=1,NPH
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
      VOLWMO=VOLWM(M,L)*FOXYX
      VOLWMM=VOLWM(M,L)*FPQ(N,L,NZ)
      VOLPMM=VOLPM(M,L)*FPQ(N,L,NZ)
      VOLWSP=RTVLW(N,L,NZ)+VOLWMM
      VOLWMA=VOLWMM*trcs_VLN(ids_NH4,L)
      VOLWMB=VOLWMM*trcs_VLN(ids_NH4B,L)
      VOLWSA=RTVLWA+VOLWMA
      VOLWSB=RTVLWB+VOLWMB
      THETW1=AZMAX1(VOLWM(M,L)/VOLY(L))
      IF(THETW1.GT.THETY(L).AND.FPQ(N,L,NZ).GT.ZEROQ(NZ))THEN
        THETM=TORT(M,L)*THETW1
        RRADS=LOG((FILM(M,L)+RRADL(N,L))/RRADL(N,L))
        RTARRX=RTARR(N,L)/RRADS
        DIFOL=THETM*OLSGL1*RTARRX
        DIFCL=THETM*CQSGL1*RTARRX
        DIFZL=THETM*ZVSGL1*RTARRX
        DIFNL=THETM*ZNSGL1*RTARRX*trcs_VLN(ids_NH4,L)
        DIFNB=THETM*ZNSGL1*RTARRX*trcs_VLN(ids_NH4B,L)
        DIFHL=THETM*HLSGL1*RTARRX

        CH4G1=trc_gascl(idg_CH4,L)*VOLPMM
        Z2OG1=trc_gascl(idg_N2O,L)*VOLPMM
        ZH3G1=trc_gascl(idg_NH3,L)*VOLPMM
        H2GG1=trc_gascl(idg_H2,L)*VOLPMM

        VOLWCO=VOLWMM*GSolbility(idg_CO2,L)
        VOLWOX=VOLWMM*GSolbility(idg_O2,L)
        VOLWCH=VOLWMM*GSolbility(idg_CH4,L)
        VOLWN2=VOLWMM*GSolbility(idg_N2O,L)
        VOLWNH=VOLWMM*GSolbility(idg_NH3,L)*trcs_VLN(ids_NH4,L)
        VOLWNB=VOLWMM*GSolbility(idg_NH3,L)*trcs_VLN(ids_NH4B,L)
        VOLWHG=VOLWMM*GSolbility(idg_H2,L)
        VOLPNH=VOLPMM*trcs_VLN(ids_NH4,L)
        VOLPNB=VOLPMM*trcs_VLN(ids_NH4B,L)
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
        D90: DO MX=1,NPT
          OXYS1=OXYS1+ROXYLX
          CCO2S1=AZMAX1(CO2S1/VOLWMM)
          COXYS1=AMIN1(COXYE*GSolbility(idg_O2,L),AZMAX1(OXYS1/VOLWMO))
          CCH4S1=AZMAX1(CH4S1/VOLWMM)
          CN2OS1=AZMAX1(Z2OS1/VOLWMM)
          CNH3S1=AZMAX1(ZH3S1/VOLWMM)
          CNH3B1=AZMAX1(ZH3B1/VOLWMM)
          CH2GS1=AZMAX1(H2GS1/VOLWMM)
          IF(RTVLP(N,L,NZ).GT.ZERO)THEN
            DO NTG=idg_beg,idg_end-1
              trcg_gcon(NTG)=AZMAX1(trcg_gmas(NTG)/RTVLP(N,L,NZ))
            ENDDO
          ELSE
            trcg_gcon(idg_beg:idg_end-1)=0.0_r8
          ENDIF
          CCO2P1=AZMAX1(CO2P1/RTVLW(N,L,NZ))
          COXYP1=AMIN1(COXYE*GSolbility(idg_O2,L),AZMAX1(OXYP1/RTVLW(N,L,NZ)))
          CCH4P1=AZMAX1(CH4P1/RTVLW(N,L,NZ))
          CN2OP1=AZMAX1(Z2OP1/RTVLW(N,L,NZ))
          CNH3P1=AZMAX1(ZH3P1/RTVLW(N,L,NZ))
          CH2GP1=AZMAX1(H2GP1/RTVLW(N,L,NZ))

          DIFOX=DIFOL+DIFOP
          RMFCOS=UPWTRH*CCO2S1
          RMFOXS=UPWTRH*COXYS1
          RMFCHS=UPWTRH*CCH4S1
          RMFN2S=UPWTRH*CN2OS1
          RMFN3S=UPWTRH*CNH3S1*trcs_VLN(ids_NH4,L)
          RMFN3B=UPWTRH*CNH3B1*trcs_VLN(ids_NH4B,L)
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
          IF(X.GT.ZERO.AND.OXYS1.GT.ZEROP(NZ))THEN
            B=-UPMXP-DIFOX*OXKM-X
            C=X*UPMXP
            RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
            COXYR=(X-RUPOXR)/DIFOX
            RDFOXS=RMFOXS+DIFOL*(COXYS1-COXYR)
            RDFOXP=DIFOP*(COXYP1-COXYR)
          ELSE
            X=DIFOP*COXYP1
            IF(X.GT.ZERO.AND.OXYP1.GT.ZEROP(NZ))THEN
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

          RUPOSX=RDFOXS*pftPlantPopulation(NZ)
          RUPOPX=RDFOXP*pftPlantPopulation(NZ)
          RDFCOS=RMFCOS+DIFCL*(CCO2S1-CCO2P1)
          RDXCOS=(RTVLW(N,L,NZ)*AMAX1(ZEROP(NZ),CO2S1) &
            -VOLWMM*AMAX1(ZEROP(NZ),CO2P1))/VOLWSP
          IF(RDFCOS.GT.0.0)THEN
            RCO2SX=AMIN1(AZMAX1(RDXCOS),RDFCOS*pftPlantPopulation(NZ))
          ELSE
            RCO2SX=AMAX1(AZMIN1(RDXCOS),RDFCOS*pftPlantPopulation(NZ))
          ENDIF
          IF(N.EQ.1)THEN
            RDFCHS=RMFCHS+DIFCL*(CCH4S1-CCH4P1)
            RDXCHS=(RTVLW(N,L,NZ)*AMAX1(ZEROP(NZ),CH4S1) &
              -VOLWMM*AMAX1(ZEROP(NZ),CH4P1))/VOLWSP
            IF(RDFCHS.GT.0.0)THEN
              RUPCSX=AMIN1(AZMAX1(RDXCHS),RDFCHS*pftPlantPopulation(NZ))
            ELSE
              RUPCSX=AMAX1(AZMIN1(RDXCHS),RDFCHS*pftPlantPopulation(NZ))
            ENDIF
            RDFN2S=RMFN2S+DIFZL*(CN2OS1-CN2OP1)
            RDXN2S=(RTVLW(N,L,NZ)*AMAX1(ZEROP(NZ),Z2OS1) &
              -VOLWMM*AMAX1(ZEROP(NZ),Z2OP1))/VOLWSP
            IF(RDFN2S.GT.0.0)THEN
              RUPZSX=AMIN1(AZMAX1(RDXN2S),RDFN2S*pftPlantPopulation(NZ))
            ELSE
              RUPZSX=AMAX1(AZMIN1(RDXN2S),RDFN2S*pftPlantPopulation(NZ))
            ENDIF
            RDFN3S=RMFN3S+DIFNL*(CNH3S1-CNH3P1)
            IF(VOLWSA.GT.ZEROP(NZ))THEN
              ZH3PA=ZH3P1*trcs_VLN(ids_NH4,L)
              RDXNHS=(RTVLWA*AMAX1(ZEROP(NZ),ZH3S1) &
                -VOLWMA*AMAX1(ZEROP(NZ),ZH3PA))/VOLWSA
            ELSE
              RDXNHS=0.0_r8
            ENDIF
            IF(RDFN3S.GT.0.0)THEN
              RUPNSX=AMIN1(AZMAX1(RDXNHS),RDFN3S*pftPlantPopulation(NZ))
            ELSE
              RUPNSX=AMAX1(AZMIN1(RDXNHS),RDFN3S*pftPlantPopulation(NZ))
            ENDIF
            RDFN3B=RMFN3B+DIFNB*(CNH3B1-CNH3P1)
            IF(VOLWSB.GT.ZEROP(NZ))THEN
              ZH3PB=ZH3P1*trcs_VLN(ids_NH4B,L)
              RDXNHB=(RTVLWB*AMAX1(ZEROP(NZ),ZH3B1) &
                -VOLWMB*AMAX1(ZEROP(NZ),ZH3PB))/VOLWSB
            ELSE
              RDXNHB=0.0_r8
            ENDIF
            IF(RDFN3B.GT.0.0)THEN
              RUPNBX=AMIN1(AZMAX1(RDXNHB),RDFN3B*pftPlantPopulation(NZ))
            ELSE
              RUPNBX=AMAX1(AZMIN1(RDXNHB),RDFN3B*pftPlantPopulation(NZ))
            ENDIF
            RDFHGS=RMFHGS+DIFHL*(CH2GS1-CH2GP1)
            RDXHGS=(RTVLW(N,L,NZ)*AMAX1(ZEROP(NZ),H2GS1) &
              -VOLWMM*AMAX1(ZEROP(NZ),H2GP1))/VOLWSP
            IF(RDFHGS.GT.0.0)THEN
              RUPHGX=AMIN1(AZMAX1(RDXHGS),RDFHGS*pftPlantPopulation(NZ))
            ELSE
              RUPHGX=AMAX1(AZMIN1(RDXHGS),RDFHGS*pftPlantPopulation(NZ))
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
          IF(THETPM(M,L).GT.THETX)THEN
            DFGSP=FPQ(N,L,NZ)*DFGS(M,L)
            RCODFQ=DFGSP*(AMAX1(ZEROP(NZ),CO2G1)*VOLWCO &
              -(AMAX1(ZEROS,CO2S1)-RCO2SX)*VOLPMM)/(VOLWCO+VOLPMM)
            RUPOST=RUPOSX-ROXYLX
            ROXDFQ=DFGSP*(AMAX1(ZEROP(NZ),OXYG1)*VOLWOX &
              -(AMAX1(ZEROS,OXYS1)-RUPOST)*VOLPMM)/(VOLWOX+VOLPMM)
            IF(N.EQ.1)THEN
              RCHDFQ=DFGSP*(AMAX1(ZEROP(NZ),CH4G1)*VOLWCH &
                -(AMAX1(ZEROS,CH4S1)-RUPCSX)*VOLPMM)/(VOLWCH+VOLPMM)
              RN2DFQ=DFGSP*(AMAX1(ZEROP(NZ),Z2OG1)*VOLWN2 &
                -(AMAX1(ZEROS,Z2OS1)-RUPZSX)*VOLPMM)/(VOLWN2+VOLPMM)
              IF(VOLWNH+VOLPNH.GT.ZEROP(NZ))THEN
                ZH3GA=ZH3G1*trcs_VLN(ids_NH4,L)
                RNHDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROP(NZ),ZH3GA)*VOLWNH &
                  -(AMAX1(ZEROS,ZH3S1)-RUPNSX)*VOLPNH)/(VOLWNH+VOLPNH)))
              ELSE
                RNHDFQ=0.0_r8
              ENDIF
              IF(VOLWNB+VOLPNB.GT.ZEROP(NZ))THEN
                ZH3GB=ZH3G1*trcs_VLN(ids_NH4B,L)
                RNBDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX &
                  ,DFGSP*(AMAX1(ZEROP(NZ),ZH3GB)*VOLWNB &
                  -(AMAX1(ZEROS,ZH3B1)-RUPNBX)*VOLPNB)/(VOLWNB+VOLPNB)))
              ELSE
                RNBDFQ=0.0_r8
              ENDIF
              RHGDFQ=DFGSP*(AMAX1(ZEROP(NZ),H2GG1)*VOLWHG &
                -(AMAX1(ZEROS,H2GS1)-RUPHGX)*VOLPMM)/(VOLWHG+VOLPMM)
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
          IF(N.EQ.1.AND.RTVLP(N,L,NZ).GT.ZEROP(NZ))THEN
            RUPNTX=RUPNSX+RUPNBX
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
!
!     R*DF1=root gas exchange between gaseous-aqueous phases
!     R*FL1=root gas exchange with atmosphere
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     trcg_gmas(idg_CO2),CO2P1=gaseous,aqueous CO2 in root
!     trcg_gmas(idg_O2),OXYP1=gaseous,aqueous O2 in root
!     trcg_gmas(idg_CH4),CH4P1=gaseous,aqueous CH4 in root
!     trcg_gmas(idg_N2O),Z2OP1=gaseous,aqueous N2O in root
!     trcg_gmas(idg_NH3),ZH3P1=gaseous,aqueous NH3 in root
!     trcg_gmas(idg_H2),H2GP1=gaseous,aqueous H2 in root
!     RTVLW,RTVLP=root aqueous,gaseous volume
!     VOLW*=RTVLW*gas solubility
!     C*E,C*A1=atmosphere,root gas concentration
!     DF*A=root-atmosphere gas conductance
!
            CO2PX=CO2P1+RCO2PX
            trcg_RDF1(idg_CO2)=AMAX1(-CO2PX,DFGP*(AMAX1(ZEROP(NZ),trcg_gmas(idg_CO2))*VOLWG(idg_CO2) &
              -CO2PX*RTVLP(N,L,NZ))/(VOLWG(idg_CO2)+RTVLP(N,L,NZ)))
            OXYPX=OXYP1-RUPOPX
            trcg_RDF1(idg_O2)=AMAX1(-OXYPX,DFGP*(AMAX1(ZEROP(NZ),trcg_gmas(idg_O2))*VOLWG(idg_O2) &
              -OXYPX*RTVLP(N,L,NZ))/(VOLWG(idg_O2)+RTVLP(N,L,NZ)))
            CH4PX=CH4P1+RUPCSX
            trcg_RDF1(idg_CH4)=AMAX1(-CH4PX,DFGP*(AMAX1(ZEROP(NZ),trcg_gmas(idg_CH4))*VOLWG(idg_CH4) &
              -CH4PX*RTVLP(N,L,NZ))/(VOLWG(idg_CH4)+RTVLP(N,L,NZ)))
            Z2OPX=Z2OP1+RUPZSX
            trcg_RDF1(idg_N2O)=AMAX1(-Z2OPX,DFGP*(AMAX1(ZEROP(NZ),trcg_gmas(idg_N2O))*VOLWG(idg_N2O) &
              -Z2OPX*RTVLP(N,L,NZ))/(VOLWG(idg_N2O)+RTVLP(N,L,NZ)))
            ZH3PX=ZH3P1+RUPNTX
            trcg_RDF1(idg_NH3)=AMAX1(-ZH3PX,DFGP*(AMAX1(ZEROP(NZ),trcg_gmas(idg_NH3))*VOLWG(idg_NH3) &
              -ZH3PX*RTVLP(N,L,NZ))/(VOLWG(idg_NH3)+RTVLP(N,L,NZ)))
            H2GPX=H2GP1+RUPHGX
            trcg_RDF1(idg_H2)=AMAX1(-H2GPX,DFGP*(AMAX1(ZEROP(NZ),trcg_gmas(idg_H2))*VOLWG(idg_H2) &
              -H2GPX*RTVLP(N,L,NZ))/(VOLWG(idg_H2)+RTVLP(N,L,NZ)))

            trcg_RFL1(idg_CO2)=AMIN1(DFCOA,RTVLP(N,L,NZ))*(CCO2E-trcg_gcon(idg_CO2))
            trcg_RFL1(idg_O2)=AMIN1(DFOXA,RTVLP(N,L,NZ))*(COXYE-trcg_gcon(idg_O2))
            trcg_RFL1(idg_CH4)=AMIN1(DFCHA,RTVLP(N,L,NZ))*(CCH4E-trcg_gcon(idg_CH4))
            trcg_RFL1(idg_N2O)=AMIN1(DFN2A,RTVLP(N,L,NZ))*(CZ2OE-trcg_gcon(idg_N2O))
            trcg_RFL1(idg_NH3)=AMIN1(DFNHA,RTVLP(N,L,NZ))*(CNH3E-trcg_gcon(idg_NH3))
            trcg_RFL1(idg_H2)=AMIN1(DFHGA,RTVLP(N,L,NZ))*(CH2GE-trcg_gcon(idg_H2))
          ELSE
            trcg_RDF1(idg_beg:idg_end-1)=0.0_r8
            trcg_RFL1(idg_beg:idg_end-1)=0.0_r8
          ENDIF
!
!     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS
!     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
!
          DO NTG=idg_beg,idg_end-1
            trcg_gmas(NTG)=trcg_gmas(NTG)-trcg_RDF1(NTG)+trcg_RFL1(NTG)
          ENDDO

          CO2P1=CO2P1+trcg_RDF1(idg_CO2)+RCO2SX+RCO2PX
          OXYP1=OXYP1+trcg_RDF1(idg_O2)-RUPOPX
          CH4P1=CH4P1+trcg_RDF1(idg_CH4)+RUPCSX
          Z2OP1=Z2OP1+trcg_RDF1(idg_N2O)+RUPZSX
          ZH3P1=ZH3P1+trcg_RDF1(idg_NH3)+RUPNSX+RUPNBX
          H2GP1=H2GP1+trcg_RDF1(idg_H2)+RUPHGX
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
          RCO2S(N,L,NZ)=RCO2S(N,L,NZ)+RCO2SX
          RUPOXS(N,L,NZ)=RUPOXS(N,L,NZ)+RUPOSX
          RUPCHS(N,L,NZ)=RUPCHS(N,L,NZ)+RUPCSX
          RUPN2S(N,L,NZ)=RUPN2S(N,L,NZ)+RUPZSX
          RUPN3S(N,L,NZ)=RUPN3S(N,L,NZ)+RUPNSX
          RUPN3B(N,L,NZ)=RUPN3B(N,L,NZ)+RUPNBX
          RUPHGS(N,L,NZ)=RUPHGS(N,L,NZ)+RUPHGX
!
!     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE
!
!     R*DFA=root aqueous-gaseous CO2 exchange
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!
          DO NTG=idg_beg,idg_end-1
            trcg_RDFA(NTG,N,L,NZ)=trcg_RDFA(NTG,N,L,NZ)+trcg_RDF1(NTG)
            trcg_RFLA(NTG,N,L,NZ)=trcg_RFLA(NTG,N,L,NZ)+trcg_RFL1(NTG)
          ENDDO
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2P=root CO2 emission into root
!     RUPOXP=root O2 uptake from root
!     ROXSK=total O2 uptake from soil by all microbial,root popns
!
          RCO2P(N,L,NZ)=RCO2P(N,L,NZ)+RCO2PX+RCO2SX
          RUPOXP(N,L,NZ)=RUPOXP(N,L,NZ)+RUPOPX
          ROXSK(M,L)=ROXSK(M,L)+RUPOSX

        ENDDO D90
      ENDIF
    ENDDO D99
!
!     O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO
!     OF ROOT O2 UPTAKE 'RUPOXT' TO ROOT O2 DEMAND 'ROXYP'
!
!     RUPOXT=O2 uptake from soil+root by each root,myco population
!     ROXYP=O2 demand by each root,myco population
!     WFR=constraint by O2 consumption on all root processes
!     imposed by O2 uptake
!
    PopPlantO2Uptake_vr=RUPOXP(N,L,NZ)+RUPOXS(N,L,NZ)
    WFR(N,L,NZ)=AMIN1(1.0_r8,AZMAX1(PopPlantO2Uptake_vr/ROXYP(N,L,NZ)))
  ELSE
    PopPlantO2Uptake_vr=0.0_r8
    IF(L.GT.NG(NZ))THEN
      WFR(N,L,NZ)=WFR(N,L-1,NZ)
    ELSE
      WFR(N,L,NZ)=1.0
    ENDIF
  ENDIF
  end associate
  end subroutine RootSoilGasExchange

end module RootGasMod
