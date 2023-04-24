module RootMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod  , only : test_aeqb,safe_adb,AZMAX1,AZMIN1
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  use NoduleBGCMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: RootBGCModel
  contains

  subroutine RootBGCModel(I,J,NZ,IFLGZ,ICHK1,IDTHRN,NRX,PTRT,TFN6,CNRTW,CPRTW,XRTN1)

  implicit none
  integer, intent(in) :: I,J,NZ,IFLGZ
  integer, intent(inout) :: ICHK1(2,JZ1)
  integer, intent(inout)  :: NRX(2,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),CNRTW,CPRTW,XRTN1
  real(r8), intent(in):: PTRT

  integer :: IDTHRN
  real(r8) :: RTVL
  real(r8) :: WFNGR(2,JZ1)
  real(r8) :: RLNT(2,JZ1)
  real(r8) :: RTSK1(2,JZ1,10)
  real(r8) :: RTSK2(2,JZ1,10)
  real(r8) :: RTNT(2)

  associate(                                &
    WTRVE    =>   plt_biom%WTRVE      , &
    ZEROL    =>   plt_biom%ZEROL      , &
    DLYR3    =>   plt_site%DLYR3      , &
    PP       =>   plt_site%PP         , &
    ZERO     =>   plt_site%ZERO       , &
    ISTYP    =>   plt_pheno%ISTYP     , &
    IDTHR    =>   plt_pheno%IDTHR     , &
    IDTHP    =>   plt_pheno%IDTHP     , &
    SDLG     =>   plt_morph%SDLG      , &
    RTVLW    =>   plt_morph%RTVLW     , &
    RTARP    =>   plt_morph%RTARP     , &
    RTVLP    =>   plt_morph%RTVLP     , &
    RTDNP    =>   plt_morph%RTDNP     , &
    RTLGP    =>   plt_morph%RTLGP     , &
    NG       =>   plt_morph%NG        , &
    PORT     =>   plt_morph%PORT      , &
    SDVL     =>   plt_morph%SDVL      , &
    SDAR     =>   plt_morph%SDAR      , &
    NRT      =>  plt_morph%NRT          &
  )
!     ROOT GROWTH
!
  call RootBiochemistry(I,J,NZ,ICHK1,IDTHRN,NRX,TFN6,CNRTW,CPRTW,&
      XRTN1,WFNGR,RLNT,RTSK1,RTSK2,RTNT)
!
!     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
!     GERMINATION)
!
  RTLGP(ipltroot,NG(NZ),NZ)=RTLGP(ipltroot,NG(NZ),NZ)+SDLG(NZ)
  IF(DLYR3(NG(NZ)).GT.ZERO)THEN
    RTDNP(ipltroot,NG(NZ),NZ)=RTLGP(ipltroot,NG(NZ),NZ)/DLYR3(NG(NZ))
  ELSE
    RTDNP(ipltroot,NG(NZ),NZ)=0._r8
  ENDIF
  RTVL=RTVLP(ipltroot,NG(NZ),NZ)+RTVLW(ipltroot,NG(NZ),NZ)+SDVL(NZ)*PP(NZ)
  RTVLP(ipltroot,NG(NZ),NZ)=PORT(ipltroot,NZ)*RTVL
  RTVLW(ipltroot,NG(NZ),NZ)=(1.0_r8-PORT(ipltroot,NZ))*RTVL
  RTARP(ipltroot,NG(NZ),NZ)=RTARP(ipltroot,NG(NZ),NZ)+SDAR(NZ)
  IF(IDTHRN.EQ.NRT(NZ).OR.(WTRVE(ielmc,NZ).LE.ZEROL(NZ).AND.ISTYP(NZ).NE.iplt_annual))THEN
    IDTHR(NZ)=idead
    IDTHP(NZ)=idead
  ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
  call RootNoduleBiomchemistry(I,J,NZ,TFN6,WFNGR)  

  call NonstructlBiomTransfer(I,J,NZ,PTRT,RLNT,RTSK1,RTSK2,RTNT,IFLGZ)
  end associate
  end subroutine RootBGCModel

!------------------------------------------------------------------------------------------

  subroutine RootBiochemistry(I,J,NZ,ICHK1,IDTHRN,NRX,TFN6,CNRTW,CPRTW,XRTN1,WFNGR,RLNT,RTSK1,RTSK2,RTNT)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: ICHK1(2,JZ1)
  integer, intent(out) :: IDTHRN
  INTEGER, Intent(inout) :: NRX(2,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),CNRTW,CPRTW,XRTN1
  real(r8), intent(out) :: WFNGR(2,JZ1)
  REAL(R8), INTENT(OUT)  :: RLNT(2,JZ1)
  real(r8), INTENT(OUT) :: RTNT(2)
  real(r8), INTENT(OUT) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10)
  integer :: LL,LZ,L1,L,K,lx,M,NR,N,NTG
  real(r8) :: WFNR
  real(r8) :: WFNRG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLD
  real(r8) :: CPOOLX
  real(r8) :: DMRTD
  real(r8) :: FWTRT
  real(r8) :: RSCS2
  real(r8) :: RTLGL,RTLGZ
  real(r8) :: RTLGX
  real(r8) :: RTLGT
  real(r8) :: RTVL
  real(r8) :: RTAR
  real(r8) :: WTRTTX
  real(r8) :: WTRVCX
  real(r8) :: WTRTX,WTRTZ
  real(r8) :: WTRTLX
  real(r8) :: WTRTTT
  real(r8) :: WTRTT
  real(r8) :: XFRC

!     begin_execution
  associate(                            &
    EPOOLR   =>   plt_biom%EPOOLR     , &
    WTRTL    =>   plt_biom%WTRTL      , &
    WTRTE    =>   plt_biom%WTRTE      , &
    ZEROP    =>   plt_biom%ZEROP      , &
    WTRVE    =>   plt_biom%WTRVE      , &
    FWODRE   =>   plt_allom%FWODRE    , &
    DMRT     =>   plt_allom%DMRT      , &
    IGTYP    =>   plt_pheno%IGTYP     , &
    PSIRT    =>   plt_ew%PSIRT        , &
    PSIRG    =>   plt_ew%PSIRG        , &
    RFGas_root    =>   plt_bgcr%RFGas_root    , &
    trcg_rootml     =>   plt_rbgc%trcg_rootml       , &
    trcs_rootml  => plt_rbgc%trcs_rootml, &
    RSCS     =>   plt_soilchem%RSCS   , &
    VOLX     =>   plt_soilchem%VOLX   , &
    NU       =>   plt_site%NU         , &
    ZERO     =>   plt_site%ZERO       , &
    PP       =>   plt_site%PP         , &
    ZEROS2   =>   plt_site%ZEROS2     , &
    DLYR3    =>   plt_site%DLYR3      , &
    NL       =>   plt_site%NL         , &
    k_fine_litr=> pltpar%k_fine_litr  , &
    RTVLP    =>   plt_morph%RTVLP     , &
    RTDNP    =>   plt_morph%RTDNP     , &
    RTARP    =>   plt_morph%RTARP     , &
    RTVLW    =>   plt_morph%RTVLW     , &
    RRAD1X   =>   plt_morph%RRAD1X    , &
    RRAD2X   =>   plt_morph%RRAD2X    , &
    RRAD1M   =>   plt_morph%RRAD1M    , &
    RRAD1    =>   plt_morph%RRAD1     , &
    RTLGP    =>   plt_morph%RTLGP     , &
    RTLGA    =>   plt_morph%RTLGA     , &
    RRAD2    =>   plt_morph%RRAD2     , &
    RTAR2X   =>   plt_morph%RTAR2X    , &
    RRAD2M   =>   plt_morph%RRAD2M    , &
    RTNL     =>   plt_morph%RTNL      , &
    RTAR1X   =>   plt_morph%RTAR1X    , &
    PORT     =>   plt_morph%PORT      , &
    NI       =>   plt_morph%NI        , &
    DMVL     =>   plt_morph%DMVL      , &
    SDLG     =>   plt_morph%SDLG      , &
    MY       =>   plt_morph%MY        , &
    NG       =>   plt_morph%NG        , &
    NIX      =>   plt_morph%NIX         &
  )

  NIX(NZ)=NG(NZ)
  IDTHRN=0
!
  call SummarizeRootSink(NZ,XRTN1,RLNT,RTSK1,RTSK2,RTNT)
!
!     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
!
  D5010: DO N=1,MY(NZ)
    D5000: DO L=NU,NI(NZ)
!
!     IDENTIFY NEXT LOWER ROOT LAYER
!
!     VOLX=soil layer volume excluding macropore, rocks
!
      IF(VOLX(L).GT.ZEROS2)THEN
        D5003: DO LZ=L+1,NL
          IF(VOLX(LZ).GT.ZEROS2.OR.LZ.EQ.NL)THEN
            L1=LZ
            EXIT
          ENDIF
        ENDDO D5003
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     RSCS,RSCS2=soil resistance to secondary root penetration (MPa)
!     RRAD2=secondary root radius
!     WFNR=water function for root extension
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     WFNGR,WFNRG=growth,respiration function of root water potential
!     PSIRT,PSIRG=root total,turgor water potential
!     DMRT=root growth yield
!
        RSCS2=RSCS(L)*RRAD2(N,L,NZ)/1.0E-03_r8
        WFNR=AMIN1(1.0_r8,AZMAX1(PSIRG(N,L,NZ)-PSILM-RSCS2))
        IF(IGTYP(NZ).EQ.0)THEN
          WFNGR(N,L)=EXP(0.05*PSIRT(N,L,NZ))
          WFNRG=WFNR**0.10_r8
        ELSE
          WFNGR(N,L)=EXP(0.10_r8*PSIRT(N,L,NZ))
          WFNRG=WFNR**0.25_r8
        ENDIF
        DMRTD=1.0_r8-DMRT(NZ)
!
!     FOR EACH ROOT AXIS
!
        call GrowRootAxes(N,L,L1,NZ,NRX,WFNGR,ICHK1,WFNR,WFNRG,TFN6,XRTN1,DMRTD,&
          RLNT,RTSK1,RTSK2,CNRTW,CPRTW,RTLGZ,WTRTX,WTRTZ,RTLGL)

!
!     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
!     SEASONAL STORAGE POOL IS DEPLETED
!
!     WTRTL,WTRT=total root C mass
!     WTRVC=storage C
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!     CPOOLR=non-structural C mass in root
!
        IF(L.LE.NIX(NZ))THEN
          IF(WTRTL(N,L,NZ).GT.ZEROP(NZ).AND.WTRTE(ielmc,NZ).GT.ZEROP(NZ) &
            .AND.WTRVE(ielmc,NZ).LT.XFRX*WTRTE(ielmc,NZ))THEN
            FWTRT=WTRTL(N,L,NZ)/WTRTE(ielmc,NZ)
            WTRTLX=WTRTL(N,L,NZ)
            WTRTTX=WTRTE(ielmc,NZ)*FWTRT
            WTRTTT=WTRTLX+WTRTTX
            CPOOLX=AZMAX1(EPOOLR(ielmc,N,L,NZ))
            WTRVCX=AZMAX1(WTRVE(ielmc,NZ)*FWTRT)
            CPOOLD=(WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
            XFRC=AZMIN1(XFRY*CPOOLD)
            EPOOLR(ielmc,N,L,NZ)=EPOOLR(ielmc,N,L,NZ)+XFRC
            WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)-XFRC
          ENDIF
        ENDIF
!
!     ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
!     TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
!
!     RTLGZ=total primary root length
!     WTRTZ=total primary root C mass
!     RTLGL=total secondary root length
!     WTRTX=total secondary root C mass
!     RTLGT=total root length
!     WTRTT=total root C mass
!     FWOOD=C woody fraction in root:0=woody,1=non-woody
!     PP=PFT population
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVL,RTVLW,RTVLP=root or myco total,aqueous,gaseous volume
!     RRAD1,RRAD2=primary,secondary root radius
!     RTARP=root surface area per plant
!     RTLGA=average secondary root length
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!
        IF(N.EQ.ipltroot)THEN
          RTLGZ=RTLGZ*FWODRE(ielmc,k_fine_litr)
          RTLGL=RTLGL*FWODRE(ielmc,k_fine_litr)
        ENDIF
        RTLGX=RTLGZ*PP(NZ)
        RTLGT=RTLGL+RTLGX
        WTRTT=WTRTX+WTRTZ
        IF(RTLGT.GT.ZEROP(NZ).AND.WTRTT.GT.ZEROP(NZ) &
          .AND.PP(NZ).GT.ZEROP(NZ))THEN
          RTLGP(N,L,NZ)=RTLGT/PP(NZ)
          IF(DLYR3(L).GT.ZERO)THEN
            RTDNP(N,L,NZ)=RTLGP(N,L,NZ)/DLYR3(L)
          ELSE
            RTDNP(N,L,NZ)=0._r8
          ENDIF
          RTVL=AMAX1(RTAR1X(N,NZ)*RTLGX+RTAR2X(N,NZ)*RTLGL &
            ,WTRTT*DMVL(N,NZ)*PSIRG(N,L,NZ))
          RTVLP(N,L,NZ)=PORT(N,NZ)*RTVL
          RTVLW(N,L,NZ)=(1.0_r8-PORT(N,NZ))*RTVL
!primary roots
          RRAD1(N,L,NZ)=AMAX1(RRAD1X(N,NZ),(1.0_r8+PSIRT(N,L,NZ)/EMODR)*RRAD1M(N,NZ))
!secondary roots
          RRAD2(N,L,NZ)=AMAX1(RRAD2X(N,NZ),(1.0_r8+PSIRT(N,L,NZ)/EMODR)*RRAD2M(N,NZ))
          RTAR=PICON2s*(RRAD1(N,L,NZ)*RTLGX+RRAD2(N,L,NZ)*RTLGL)
          IF(RTNL(N,L,NZ).GT.ZEROP(NZ))THEN
            RTLGA(N,L,NZ)=AMAX1(RTLGAX,RTLGL/RTNL(N,L,NZ))
          ELSE
            RTLGA(N,L,NZ)=RTLGAX
          ENDIF
          RTARP(N,L,NZ)=RTAR/PP(NZ)
!     IF(N.EQ.1)THEN
!     RTARP(N,L,NZ)=RTARP(N,L,NZ)*RTLGAX/RTLGA(N,L,NZ)
!     ENDIF
        ELSE
          RTLGP(N,L,NZ)=0._r8
          RTDNP(N,L,NZ)=0._r8
          RTVLP(N,L,NZ)=0._r8
          RTVLW(N,L,NZ)=0._r8
          RRAD1(N,L,NZ)=RRAD1M(N,NZ)
          RRAD2(N,L,NZ)=RRAD2M(N,NZ)
          RTARP(N,L,NZ)=0._r8
          RTLGA(N,L,NZ)=RTLGAX
          DO NTG=idg_beg,idg_end-1
            RFGas_root(NTG,NZ)=RFGas_root(NTG,NZ)-(trcg_rootml(NTG,N,L,NZ)+trcs_rootml(NTG,N,L,NZ))
          ENDDO
          trcg_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
          trcs_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
        ENDIF
      ENDIF
    ENDDO D5000
  ENDDO D5010
  end associate
  end subroutine RootBiochemistry

!------------------------------------------------------------------------------------------

  subroutine GrowRootAxes(N,L,L1,NZ,NRX,WFNGR,ICHK1,WFNR,WFNRG,TFN6,XRTN1,DMRTD,RLNT,&
    RTSK1,RTSK2,CNRTW,CPRTW,RTLGL,RTLGZ,WTRTX,WTRTZ)
  implicit none
  INTEGER, INTENT(IN) :: N,L,L1,NZ
  integer, intent(inout) :: NRX(2,JZ1)
  real(r8), intent(in) :: TFN6(JZ1),XRTN1
  real(r8), intent(in) :: WFNGR(2,JZ1)
  real(r8), intent(in) :: DMRTD
  REAL(R8), INTENT(IN) :: RLNT(2,JZ1)
  real(r8), intent(in) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10),CNRTW,CPRTW
  integer, intent(inout) :: ICHK1(2,JZ1)
  real(r8), intent(inout):: WFNR,WFNRG
  real(r8), intent(out) :: RTLGL,RTLGZ,WTRTX,WTRTZ
  real(r8) :: CNRDA,CNRDM
  real(r8) :: CNPG
  real(r8) :: CCC,CNC,CPC
  real(r8) :: FSNC1,FSNC2
  real(r8) :: ZADD1,ZADD1M,ZADD2,ZADD2M
  real(r8) :: ZPOOLB
  real(r8) :: CGRORM
  real(r8) :: CGROR
  real(r8) :: DMRTR
  real(r8) :: FNP
  real(r8) :: FRTN
  real(r8) :: FRCO2
  real(r8) :: FSNCM
  real(r8) :: FSNCP
  real(r8) :: GRTWGM
  real(r8) :: GRTLGL
  real(r8) :: GRTWTG
  real(r8) :: GRTWTLE(npelms)
  real(r8) :: GRTWTM
  real(r8) :: PPOOLB
  real(r8) :: PADD2,PADD1
  real(r8) :: RCO2X,RCO2Y
  real(r8) :: RCO2G,RCO2T
  real(r8) :: RCO2XM
  real(r8) :: RCO2YM
  real(r8) :: RCO2GM
  real(r8) :: RCO2TM
  real(r8) :: RMNCR,RCO2RM,RCO2R
  real(r8) :: RCER(npelms)
  real(r8) :: RTN2X,RTN2Y
  real(r8) :: RTDP1X,RSCS1
  REAL(R8) :: SNCR,SNCRM
  real(r8) :: TFRCO2
  real(r8) :: RCCC,RCCN,RCCP
  integer :: NR,M,LL,LX,NE

!begin_execution
  associate(                          &
    RTWT1E  =>  plt_biom%RTWT1E     , &
    WTRT2E  =>  plt_biom%WTRT2E     , &
    CEPOLR  =>  plt_biom%CEPOLR     , &
    WTRT1E  =>  plt_biom%WTRT1E     , &
    EPOOLR  =>  plt_biom%EPOOLR     , &
    WSRTL   =>  plt_biom%WSRTL      , &
    WTRTL   =>  plt_biom%WTRTL      , &
    ZEROP   =>  plt_biom%ZEROP      , &
    CDPTHZ  =>  plt_site%CDPTHZ     , &
    RCO2A   =>  plt_rbgc%RCO2A      , &
    RCO2N   =>  plt_rbgc%RCO2N      , &
    RCO2M   =>  plt_rbgc%RCO2M      , &
    WFR     =>  plt_rbgc%WFR        , &
    ESNC    =>  plt_bgcr%ESNC       , &
    CNWS    =>  plt_allom%CNWS      , &
    CPWS    =>  plt_allom%CPWS      , &
    FWODRE  =>  plt_allom%FWODRE    , &
    CNRTS   =>  plt_allom%CNRTS     , &
    CPRTS   =>  plt_allom%CPRTS     , &
    DMRT    =>  plt_allom%DMRT      , &
    k_woody_litr=> pltpar%k_woody_litr,&
    k_fine_litr=> pltpar%k_fine_litr, &
    icwood  =>  pltpar%icwood       , &
    iroot   =>  pltpar%iroot        , &
    instruct=>  pltpar%instruct     , &
    IGTYP   =>  plt_pheno%IGTYP     , &
    IWTYP   =>  plt_pheno%IWTYP     , &
    TFN4    =>  plt_pheno%TFN4      , &
    IDAY    =>  plt_pheno%IDAY      , &
    BKDS    =>  plt_soilchem%BKDS   , &
    CFOPE   =>  plt_soilchem%CFOPE  , &
    RSCS    =>  plt_soilchem%RSCS   , &
    DLYR3   =>  plt_site%DLYR3      , &
    ZERO    =>  plt_site%ZERO       , &
    NJ      =>  plt_site%NJ         , &
    PSIRG   =>  plt_ew%PSIRG        , &
    RTNL    =>  plt_morph%RTNL      , &
    GRMX    =>  plt_morph%GRMX      , &
    RTDP1   =>  plt_morph%RTDP1     , &
    RTN1    =>  plt_morph%RTN1      , &
    NG      =>  plt_morph%NG        , &
    NIX     =>  plt_morph%NIX       , &
    RTN2    =>  plt_morph%RTN2      , &
    NRT     =>  plt_morph%NRT       , &
    RRAD1   =>  plt_morph%RRAD1     , &
    RTFQ    =>  plt_morph%RTFQ      , &
    SDPTH   =>  plt_morph%SDPTH     , &
    RTLG1   =>  plt_morph%RTLG1     , &
    RTLG2X  =>  plt_morph%RTLG2X    , &
    RTLG2   =>  plt_morph%RTLG2     , &
    NINR    =>  plt_morph%NINR      , &
    NB1     =>  plt_morph%NB1       , &
    FDBKX   =>  plt_photo%FDBKX       &
  )
  RTLGL=0._r8
  RTLGZ=0._r8
  WTRTX=0._r8
  WTRTZ=0._r8
  D5050: DO NR=1,NRT(NZ)
!
!     SECONDARY ROOT EXTENSION
!
    IF(L.LE.NINR(NR,NZ).AND.NRX(N,NR).EQ.0)THEN
!
!     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
!     TO CURRENT AXIS
!
!     RTSK2=total secondary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of secondary root sink strength in axis
!
      IF(RLNT(N,L).GT.ZEROP(NZ))THEN
        FRTN=RTSK2(N,L,NR)/RLNT(N,L)
      ELSE
        FRTN=1.0_r8
      ENDIF
!
!     N,P CONSTRAINT ON SECONDARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
      IF(CEPOLR(ielmc,N,L,NZ).GT.ZERO)THEN
        CNPG=AMIN1(CEPOLR(ielmn,N,L,NZ)/(CEPOLR(ielmn,N,L,NZ) &
          +CEPOLR(ielmc,N,L,NZ)*CNKI),CEPOLR(ielmp,N,L,NZ) &
          /(CEPOLR(ielmp,N,L,NZ)+CEPOLR(ielmc,N,L,NZ)*CPKI))
      ELSE
        CNPG=1.0_r8
      ENDIF
!
!     SECONDARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT2N=secondary root N mass
!     TFN6=temperature function for root maintenance respiration
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WFNGR=growth function of root water potential
!
      RMNCR=AZMAX1(RMPLT*WTRT2E(ielmn,N,L,NR,NZ))*TFN6(L)
      IF(IGTYP(NZ).EQ.0.OR.IWTYP(NZ).EQ.2)THEN
        RMNCR=RMNCR*WFNGR(N,L)
      ENDIF
!
!     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of secondary root sink strength in axis
!     CPOOL=non-structural C mass
!     TFN4=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     FDBKX=termination feedback inhibition on C3 CO2
!     WFNGR=growth function of root water potential
!
      RCO2RM=AZMAX1(VMXC*FRTN*EPOOLR(ielmc,N,L,NZ) &
        *TFN4(L,NZ))*CNPG*FDBKX(NB1(NZ),NZ)*WFNGR(N,L)
!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     WFR=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RCO2R=RCO2RM*WFR(N,L,NZ)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AZMAX1(RCO2XM)*WFNRG
      RCO2Y=AZMAX1(RCO2X)*WFNRG
!
!     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
      DMRTR=DMRTD*FRTN
      ZPOOLB=AZMAX1(EPOOLR(ielmn,N,L,NZ))
      PPOOLB=AZMAX1(EPOOLR(ielmp,N,L,NZ))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ),PPOOLB*DMRTR/CPRTS(NZ))
      IF(RCO2YM.GT.0.0_r8)THEN
        RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
        RCO2GM=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
        RCO2G=AMIN1(RCO2Y,FNP*WFR(N,L,NZ))
      ELSE
        RCO2G=0._r8
      ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,GRTWTG=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD2M,ZADD2,PADD2=nonstructural N,P unlimited,limited by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
      CGRORM=RCO2GM/DMRTD
      CGROR=RCO2G/DMRTD
      GRTWGM=CGRORM*DMRT(NZ)
      GRTWTG=CGROR*DMRT(NZ)
      ZADD2M=AZMAX1(GRTWGM*CNRTW)
      ZADD2=AZMAX1(AMIN1(FRTN*EPOOLR(ielmn,N,L,NZ),GRTWTG*CNRTW))
      PADD2=AZMAX1(AMIN1(FRTN*EPOOLR(ielmp,N,L,NZ),GRTWTG*CPRTW))
      CNRDM=AZMAX1(1.70_r8*ZADD2M)
      CNRDA=AZMAX1(1.70_r8*ZADD2)
!
!     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAY(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
      IF(IDAY(1,NB1(NZ),NZ).NE.0.AND.CEPOLR(ielmc,N,L,NZ).GT.ZERO)THEN
        CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(CEPOLR(ielmn,N,L,NZ),CEPOLR(ielmn,N,L,NZ) &
          +CEPOLR(ielmc,N,L,NZ)*CNKI) &
          ,safe_adb(CEPOLR(ielmp,N,L,NZ),CEPOLR(ielmp,N,L,NZ) &
          +CEPOLR(ielmc,N,L,NZ)*CPKI)))
        CNC=AZMAX1(AMIN1(1.0_r8 &
          ,safe_adb(CEPOLR(ielmc,N,L,NZ),CEPOLR(ielmc,N,L,NZ)+CEPOLR(ielmn,N,L,NZ)/CNKI)))
        CPC=AZMAX1(AMIN1(1.0_r8 &
          ,safe_adb(CEPOLR(ielmc,N,L,NZ),CEPOLR(ielmc,N,L,NZ)+CEPOLR(ielmp,N,L,NZ)/CPKI)))
      ELSE
        CCC=0._r8
        CNC=0._r8
        CPC=0._r8
      ENDIF
      RCCC=RCCZR+CCC*RCCYR
      RCCN=CNC*RCCXR
      RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P FROM SECONDARY ROOT DURING
!     REMOBILIZATION DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     WFR=constraint by O2 consumption on all root processes
!     RCER(ielmc),RCER(ielmn),RCER(ielmp)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC2=fraction of secondary root C to be remobilized
!
      IF(-RCO2XM.GT.0.0_r8)THEN
        IF(-RCO2XM.LT.WTRT2E(ielmc,N,L,NR,NZ)*RCCC)THEN
          SNCRM=-RCO2XM
        ELSE
          SNCRM=AZMAX1(WTRT2E(ielmc,N,L,NR,NZ)*RCCC)
        ENDIF
      ELSE
        SNCRM=0._r8
      ENDIF
      IF(-RCO2X.GT.0.0_r8)THEN
        IF(-RCO2X.LT.WTRT2E(ielmc,N,L,NR,NZ)*RCCC)THEN
          SNCR=-RCO2X
        ELSE
          SNCR=AZMAX1(WTRT2E(ielmc,N,L,NR,NZ)*RCCC)*WFR(N,L,NZ)
        ENDIF
      ELSE
        SNCR=0._r8
      ENDIF
      IF(SNCR.GT.0.0.AND.WTRT2E(ielmc,N,L,NR,NZ).GT.ZEROP(NZ))THEN
        RCER(ielmc)=RCCC*WTRT2E(ielmc,N,L,NR,NZ)
        RCER(ielmn)=WTRT2E(ielmn,N,L,NR,NZ)*(RCCN+(1.0_r8-RCCN)*RCER(ielmc)/WTRT2E(ielmc,N,L,NR,NZ))
        RCER(ielmp)=WTRT2E(ielmp,N,L,NR,NZ)*(RCCP+(1.0_r8-RCCP)*RCER(ielmc)/WTRT2E(ielmc,N,L,NR,NZ))
        IF(RCER(ielmc).GT.ZEROP(NZ))THEN
          FSNC2=AZMAX1(AMIN1(1.0_r8,SNCR/RCER(ielmc)))
        ELSE
          FSNC2=1.0_r8
        ENDIF
      ELSE
        RCER(1:npelms)=0._r8
        FSNC2=0._r8
      ENDIF
!
!     SECONDARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC2=fraction of secondary root C to be remobilized
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RCER(ielmc),RCER(ielmn),RCER(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      DO NE=1,npelms
        D6350: DO M=1,jsken
          ESNC(NE,M,k_woody_litr,L,NZ)=ESNC(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
            *FSNC2*(WTRT2E(NE,N,L,NR,NZ)-RCER(NE))*FWODRE(NE,k_woody_litr)

          ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
            *FSNC2*(WTRT2E(NE,N,L,NR,NZ)-RCER(NE))*FWODRE(NE,k_fine_litr)
        ENDDO D6350
      ENDDO
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC2=fraction of secondary root C to be remobilized
!     RCER(ielmc),RCER(ielmn),RCER(ielmp)=remobilization of C,N,P from senescing root
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      EPOOLR(ielmc,N,L,NZ)=EPOOLR(ielmc,N,L,NZ)-AMIN1(RMNCR,RCO2R) &
        -CGROR-CNRDA-SNCR+FSNC2*RCER(ielmc)
      EPOOLR(ielmn,N,L,NZ)=EPOOLR(ielmn,N,L,NZ)-ZADD2+FSNC2*RCER(ielmn)
      EPOOLR(ielmp,N,L,NZ)=EPOOLR(ielmp,N,L,NZ)-PADD2+FSNC2*RCER(ielmp)
!
!     TOTAL SECONDARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
      RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
      RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
      RCO2M(N,L,NZ)=RCO2M(N,L,NZ)+RCO2TM
      RCO2N(N,L,NZ)=RCO2N(N,L,NZ)+RCO2T
      RCO2A(N,L,NZ)=RCO2A(N,L,NZ)-RCO2T
!
!     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     GRTLGL=secondary root length extension
!     GRTWTG=secondary root C growth ltd by O2
!     RTLG2X=specific secondary root length from startq.f
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     FSNC2=fraction of secondary root C to be remobilized
!     RTLG2=secondary root length
!     GRTWTLE(ielmc),GRTWTLE(ielmn),GRTWTLE(ielmp)=net root C,N,P growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      GRTLGL=GRTWTG*RTLG2X(N,NZ)*WFNR*FWODRE(ielmc,k_fine_litr) &
        -FSNC2*RTLG2(N,L,NR,NZ)
      GRTWTLE(ielmc)=GRTWTG-FSNC2*WTRT2E(ielmc,N,L,NR,NZ)
      GRTWTLE(ielmn)=ZADD2-FSNC2*WTRT2E(ielmn,N,L,NR,NZ)
      GRTWTLE(ielmp)=PADD2-FSNC2*WTRT2E(ielmp,N,L,NR,NZ)
!
!     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
!     AND AXIS NUMBER
!
!     RTLG2=secondary root length
!     GRTLGL=secondary root length extension
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     GRTWTLE(ielmc),GRTWTLE(ielmn),GRTWTLE(ielmp)=net root C,N,P growth
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTFQ=root branching frequency from PFT file
!     RTN2,RTNL=number of secondary root axes
!
      RTLG2(N,L,NR,NZ)=RTLG2(N,L,NR,NZ)+GRTLGL
      DO NE=1,npelms
        WTRT2E(NE,N,L,NR,NZ)=WTRT2E(NE,N,L,NR,NZ)+GRTWTLE(NE)
      ENDDO
      WSRTL(N,L,NZ)=WSRTL(N,L,NZ)+AMIN1(CNWS(NZ)*WTRT2E(ielmn,N,L,NR,NZ) &
        ,CPWS(NZ)*WTRT2E(ielmp,N,L,NR,NZ))
      RTLGL=RTLGL+RTLG2(N,L,NR,NZ)
      WTRTX=WTRTX+WTRT2E(ielmc,N,L,NR,NZ)
      RTN2X=RTFQ(NZ)*XRTN1
      RTN2Y=RTFQ(NZ)*RTN2X
      RTN2(N,L,NR,NZ)=(RTN2X+RTN2Y)*DLYR3(L)
      RTNL(N,L,NZ)=RTNL(N,L,NZ)+RTN2(N,L,NR,NZ)
!
!     PRIMARY ROOT EXTENSION
!
!     BKDS=soil bulk density
!     RTDP1,RTDP1X=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     ICHKL=flag for identifying layer with primary root tip
!     RTN1=number of primary root axes
!     XRTN1=multiplier for number of primary root axes
!
      IF(N.EQ.ipltroot)THEN
        IF(BKDS(L).GT.ZERO)THEN
          RTDP1X=RTDP1(N,NR,NZ)-CDPTHZ(0)
        ELSE
          RTDP1X=RTDP1(N,NR,NZ)
        ENDIF
        IF(RTDP1X.GT.CDPTHZ(L-1).AND.ICHK1(N,NR).EQ.0)THEN
            RTN1(N,L,NZ)=RTN1(N,L,NZ)+XRTN1
            IF(RTDP1X.LE.CDPTHZ(L).OR.L.EQ.NJ)THEN
              ICHK1(N,NR)=1
!
!     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
!     ATTRIBUTED TO CURRENT AXIS
!
!     RTSK1=primary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of primary root sink strength in axis
!
              IF(RLNT(N,L).GT.ZEROP(NZ))THEN
                FRTN=RTSK1(N,L,NR)/RLNT(N,L)
              ELSE
                FRTN=1.0_r8
              ENDIF
!
!     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
!     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
!
!     RSCS,RSCS1=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
!     WFNR=water function for root extension
!     WFNRG=respiration function of root water potential
!
              RSCS1=RSCS(L)*RRAD1(N,L,NZ)/1.0E-03_r8
              WFNR=AMIN1(1.0_r8,AZMAX1(PSIRG(N,L,NZ)-PSILM-RSCS1))
              IF(IGTYP(NZ).EQ.0)THEN
                WFNRG=WFNR**0.10_r8
              ELSE
                WFNRG=WFNR**0.25_r8
              ENDIF
!
!     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
              IF(CEPOLR(ielmc,N,L,NZ).GT.ZERO)THEN
                CNPG=AMIN1(CEPOLR(ielmn,N,L,NZ)/(CEPOLR(ielmn,N,L,NZ) &
                  +CEPOLR(ielmc,N,L,NZ)*CNKI),CEPOLR(ielmp,N,L,NZ) &
                  /(CEPOLR(ielmp,N,L,NZ)+CEPOLR(ielmc,N,L,NZ)*CPKI))
              ELSE
                CNPG=1.0_r8
              ENDIF
!
!     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     ROOT STRUCTURAL N
!
!     RMNCR=root maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     WTRT1N=primary root N mass
!     TFN6=temperature function for root maintenance respiration
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WFNGR=growth function of root water potential
!
              RMNCR=AZMAX1(RMPLT*RTWT1E(ielmn,N,NR,NZ))*TFN6(L)
              IF(IGTYP(NZ).EQ.0.OR.IWTYP(NZ).EQ.2)THEN
                RMNCR=RMNCR*WFNGR(N,L)
              ENDIF
!
!     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
!     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
!
!     RCO2RM=respiration from non-structural C unlimited by O2
!     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of primary root sink strength in axis
!     CPOOL=non-structural C mass
!     TFN4=temperature function for root growth
!     CNPG=N,P constraint on respiration
!     FDBKX=termination feedback inhibition on C3 CO2
!     WFNGR=growth function of root water potential
!
              RCO2RM=AZMAX1(VMXC*FRTN*EPOOLR(ielmc,N,L,NZ) &
                *TFN4(L,NZ))*CNPG*FDBKX(NB1(NZ),NZ)*WFNGR(N,L)
              IF(RTDP1X.GE.CDPTHZ(NJ))THEN
                RCO2RM=AMIN1(RMNCR,RCO2RM)
              ENDIF
!
!     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     WFR=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
              RCO2R=RCO2RM*WFR(N,L,NZ)
              RCO2XM=RCO2RM-RMNCR
              RCO2X=RCO2R-RMNCR
              RCO2YM=AZMAX1(RCO2XM)*WFNRG
              RCO2Y=AZMAX1(RCO2X)*WFNRG
!
!     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
!     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
!
!     FRTN=fraction of secondary root sink strength in axis
!     ZPOOLR,PPOOLR=non-structural N,P mass in root
!     CNRTS,CPRTS=N,P root growth yield
!     FNP=growth respiration limited by non-structural N,P
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!
              DMRTR=DMRTD*FRTN
              ZPOOLB=AZMAX1(EPOOLR(ielmn,N,L,NZ))
              PPOOLB=AZMAX1(EPOOLR(ielmp,N,L,NZ))
              FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ),PPOOLB*DMRTR/CPRTS(NZ))
              IF(RCO2YM.GT.0.0_r8)THEN
                RCO2GM=AMIN1(RCO2YM,FNP)
              ELSE
                RCO2GM=0._r8
              ENDIF
              IF(RCO2Y.GT.0.0_r8)THEN
                RCO2G=AMIN1(RCO2Y,FNP*WFR(N,L,NZ))
              ELSE
                RCO2G=0._r8
              ENDIF
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGRORM,CGROR=total non-structural C used in growth and respn unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     DMRTD=root C respiration vs nonstructural C consumption
!     GRTWGM,GRTWTG=root C growth unltd,ltd by O2
!     DMRT=root growth yield
!     ZADD1M,ZADD1,PADD1=nonstructural N,P unltd,ltd by O2 used in growth
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!
              CGRORM=RCO2GM/DMRTD
              CGROR=RCO2G/DMRTD
              GRTWGM=CGRORM*DMRT(NZ)
              GRTWTG=CGROR*DMRT(NZ)
              ZADD1M=AZMAX1(GRTWGM*CNRTW)
              ZADD1=AZMAX1(AMIN1(FRTN*EPOOLR(ielmn,N,L,NZ),GRTWTG*CNRTW))
              PADD1=AZMAX1(AMIN1(FRTN*EPOOLR(ielmp,N,L,NZ),GRTWTG*CPRTW))
              CNRDM=AZMAX1(1.70_r8*ZADD1M)
              CNRDA=AZMAX1(1.70_r8*ZADD1)

              call PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RMNCR=root maintenance respiration
!     RCO2R=respiration from non-structural C limited by O2
!     CGROR=total non-structural C used in growth and respn ltd by O2
!     CNRDA=respiration for N assimilation unltd,ltd by O2
!     SNCR=excess maintenance respiration ltd by O2
!     FSNC1=fraction of primary root C to be remobilized
!     RCER(ielmc),RCER(ielmn),RCER(ielmp)=remobilization of C,N,P from senescing root
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!
              EPOOLR(ielmc,N,L,NZ)=EPOOLR(ielmc,N,L,NZ)-AMIN1(RMNCR,RCO2R)-CGROR-CNRDA-SNCR+FSNC1*RCER(ielmc)
              EPOOLR(ielmn,N,L,NZ)=EPOOLR(ielmn,N,L,NZ)-ZADD1+FSNC1*RCER(ielmn)
              EPOOLR(ielmp,N,L,NZ)=EPOOLR(ielmp,N,L,NZ)-PADD1+FSNC1*RCER(ielmp)
!
!     TOTAL PRIMARY ROOT RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!     RMNCR=root maintenance respiration
!     RCO2RM,RCO2R=respiration from non-structural C unltd,ltd by O2
!     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     CNRDM,CNRDA=respiration for N assimilation unltd,ltd by O2
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!
              RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
              RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
!
!     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
!     THROUGH WHICH PRIMARY ROOTS GROW
!
!     RTDP1=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     RTLG1=primary root length
!     SDPTH=seeding depth
!     FRCO2=fraction of primary root respiration attributed to layer
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!     RCO2TM,RCO2T=total C respiration unltd,ltd by O2
!
              IF(RTDP1(N,NR,NZ).GT.CDPTHZ(NG(NZ)))THEN
                TFRCO2=0._r8
                D5100: DO LL=NG(NZ),NINR(NR,NZ)
                  IF(LL.LT.NINR(NR,NZ))THEN
                    FRCO2=AMIN1(1.0_r8,RTLG1(N,LL,NR,NZ)/(RTDP1(N,NR,NZ)-SDPTH(NZ)))
                  ELSE
                    FRCO2=1.0_r8-TFRCO2
                  ENDIF
                  TFRCO2=TFRCO2+FRCO2
                  RCO2M(N,LL,NZ)=RCO2M(N,LL,NZ)+RCO2TM*FRCO2
                  RCO2N(N,LL,NZ)=RCO2N(N,LL,NZ)+RCO2T*FRCO2
                  RCO2A(N,LL,NZ)=RCO2A(N,LL,NZ)-RCO2T*FRCO2
                ENDDO D5100
              ELSE
                RCO2M(N,L,NZ)=RCO2M(N,L,NZ)+RCO2TM
                RCO2N(N,L,NZ)=RCO2N(N,L,NZ)+RCO2T
                RCO2A(N,L,NZ)=RCO2A(N,L,NZ)-RCO2T
              ENDIF
!
!     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
!     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
!     HAVE DISAPPEARED
!
!     GRTWTG=primary root C growth ltd by O2
!     GRTWTLE(ielmc),GRTWTLE(ielmn),GRTWTLE(ielmp)=net primary root C,N,P growth
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RTLG2=secondary root length
!
              GRTWTLE(ielmc)=GRTWTG-FSNC1*RTWT1E(ielmc,N,NR,NZ)
              GRTWTLE(ielmn)=ZADD1-FSNC1*RTWT1E(ielmn,N,NR,NZ)
              GRTWTLE(ielmp)=PADD1-FSNC1*RTWT1E(ielmp,N,NR,NZ)
              IF(GRTWTLE(ielmc).LT.0.0_r8)THEN
                LX=MAX(1,L-1)
                D5105: DO LL=L,LX,-1
                  GRTWTM=GRTWTLE(ielmc)
                  DO NE=1,npelms
                    IF(GRTWTLE(NE).LT.0.0_r8)THEN
                      IF(GRTWTLE(NE).GT.-WTRT2E(NE,N,LL,NR,NZ))THEN
                        if(NE==ielmc)RTLG2(N,LL,NR,NZ)=RTLG2(N,LL,NR,NZ)+GRTWTLE(NE) &
                          *RTLG2(N,LL,NR,NZ)/WTRT2E(NE,N,LL,NR,NZ)

                        WTRT2E(NE,N,LL,NR,NZ)=WTRT2E(NE,N,LL,NR,NZ)+GRTWTLE(NE)
                        GRTWTLE(NE)=0._r8
                      ELSE
                        if(NE==ielmc)RTLG2(N,LL,NR,NZ)=0._r8
                        GRTWTLE(NE)=GRTWTLE(NE)+WTRT2E(NE,N,LL,NR,NZ)
                        WTRT2E(NE,N,LL,NR,NZ)=0._r8
                      ENDIF
                    ENDIF
                  ENDDO
!
!     CONCURRENT LOSS OF MYCORRHIZAE AND NODULES WITH LOSS
!     OF SECONDARY ROOTS
!
!     GRTWTM=negative primary root C growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     FSNCM,FSNCP=fraction of mycorrhizal structural,nonstructural C to be remobilized
!     WTRTL=active root C mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTRT2,WTRT2N,WTRT2P=mycorrhizal C,N,P mass
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     RTLG2=mycorrhizal length
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
!
                  IF(GRTWTM.LT.0.0_r8)THEN
                    IF(WTRT2E(ielmc,ipltroot,LL,NR,NZ).GT.ZEROP(NZ))THEN
                      FSNCM=AMIN1(1.0_r8,ABS(GRTWTM)/WTRT2E(ielmc,ipltroot,LL,NR,NZ))
                    ELSE
                      FSNCM=1.0_r8
                    ENDIF
                    IF(WTRTL(ipltroot,LL,NZ).GT.ZEROP(NZ))THEN
                      FSNCP=AMIN1(1.0_r8,ABS(GRTWTM)/WTRTL(ipltroot,LL,NZ))
                    ELSE
                      FSNCP=1.0_r8
                    ENDIF

                    DO NE=1,npelms
                      D6450: DO M=1,jsken
                        ESNC(NE,M,k_woody_litr,LL,NZ)=ESNC(NE,M,k_woody_litr,LL,NZ)+CFOPE(NE,icwood,M,NZ) &
                          *FSNCM*AZMAX1(WTRT2E(NE,imycorrhz,LL,NR,NZ))*FWODRE(NE,k_woody_litr)

                        ESNC(NE,M,k_fine_litr,LL,NZ)=ESNC(NE,M,k_fine_litr,LL,NZ)+CFOPE(NE,iroot,M,NZ) &
                          *FSNCM*AZMAX1(WTRT2E(NE,imycorrhz,LL,NR,NZ))*FWODRE(NE,k_fine_litr)

                        ESNC(NE,M,k_fine_litr,LL,NZ)=ESNC(NE,M,k_fine_litr,LL,NZ)+CFOPE(NE,instruct,M,NZ) &
                          *FSNCP*AZMAX1(EPOOLR(NE,imycorrhz,LL,NZ))
                      ENDDO D6450
                      WTRT2E(NE,imycorrhz,LL,NR,NZ)=AZMAX1(WTRT2E(NE,imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)

                      EPOOLR(NE,imycorrhz,LL,NZ)=AZMAX1(EPOOLR(NE,imycorrhz,LL,NZ))*(1.0_r8-FSNCP)

                    ENDDO
                    RTLG2(imycorrhz,LL,NR,NZ)=AZMAX1(RTLG2(imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)
                  ENDIF
                ENDDO D5105
              ENDIF
!
              call PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,GRTWTG,GRTWTLE,&
                GRTLGL,RTLGZ,WTRTZ)
            ENDIF
!
!
            IF(L.EQ.NINR(NR,NZ))THEN
              call WithdrawPrimRoot(L,NR,NZ,N,RTDP1X,RLNT,RTSK1,RTSK2,XRTN1)
            ENDIF

!
!     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
!
            IF(WTRT1E(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
              EPOOLR(ielmc,N,L,NZ)=EPOOLR(ielmc,N,L,NZ)+WTRT1E(ielmc,N,L,NR,NZ)
              WTRT1E(ielmc,N,L,NR,NZ)=0._r8
            ENDIF
            IF(WTRT2E(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
              EPOOLR(ielmc,N,L,NZ)=EPOOLR(ielmc,N,L,NZ)+WTRT2E(ielmc,N,L,NR,NZ)
              WTRT2E(ielmc,N,L,NR,NZ)=0._r8
            ENDIF
!
!     TOTAL PRIMARY ROOT LENGTH AND MASS
!
!     RTLGZ=total primary root length
!     WTRTZ=total primary root C mass
!     RTLG1=primary root length in soil layer
!     WTRT1=primary root C mass in soil layer
!     NINR=deepest root layer
!
            RTLGZ=RTLGZ+RTLG1(N,L,NR,NZ)
            WTRTZ=WTRTZ+WTRT1E(ielmc,N,L,NR,NZ)
            NINR(NR,NZ)=MIN(NINR(NR,NZ),NJ)            
            IF(L.EQ.NINR(NR,NZ))NRX(N,NR)=1
          ENDIF
        ENDIF
        RTLGZ=RTLGZ+RTLG1(N,L,NR,NZ)
        WTRTZ=WTRTZ+WTRT1E(ielmc,N,L,NR,NZ)
!     ENDIF
      ENDIF
      NIX(NZ)=MAX(NIX(NZ),NINR(NR,NZ))
  ENDDO D5050
  end associate
  end subroutine GrowRootAxes

!------------------------------------------------------------------------------------------

  subroutine PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

  implicit none
  integer, intent(in) :: N,L,NZ,NR
  real(r8), intent(in) :: RCO2X,RCO2XM
  real(r8), intent(out) :: FSNC1
  integer :: M,NE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: RCER(npelms)
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SNCR,SNCRM
! begin_execution
  associate(                          &
    RTWT1E  =>  plt_biom%RTWT1E     , &
    CEPOLR  =>  plt_biom%CEPOLR     , &
    ZEROP   =>  plt_biom%ZEROP      , &
    ZERO    =>  plt_site%ZERO       , &
    icwood  =>  pltpar%icwood       , &
    k_woody_litr=> pltpar%k_woody_litr,&
    k_fine_litr=> pltpar%k_fine_litr, &
    iroot   =>  pltpar%iroot        , &
    FWODRE  =>  plt_allom%FWODRE    , &
    CFOPE   =>  plt_soilchem%CFOPE  , &
    WFR     =>  plt_rbgc%WFR        , &
    ESNC    =>  plt_bgcr%ESNC       , &
    IDAY    =>  plt_pheno%IDAY      , &
    NB1     =>  plt_morph%NB1         &
  )
!
!     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAY(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
  IF(IDAY(1,NB1(NZ),NZ).NE.0.AND.CEPOLR(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(CEPOLR(ielmn,N,L,NZ),CEPOLR(ielmn,N,L,NZ)+CEPOLR(ielmc,N,L,NZ)*CNKI) &
      ,safe_adb(CEPOLR(ielmp,N,L,NZ),CEPOLR(ielmp,N,L,NZ)+CEPOLR(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(CEPOLR(ielmc,N,L,NZ),CEPOLR(ielmc,N,L,NZ)+CEPOLR(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(CEPOLR(ielmc,N,L,NZ),CEPOLR(ielmc,N,L,NZ)+CEPOLR(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC=0._r8
    CNC=0._r8
    CPC=0._r8
  ENDIF
  RCCC=RCCZR+CCC*RCCYR
  RCCN=CNC*RCCXR
  RCCP=CPC*RCCQR
!
!     RECOVERY OF REMOBILIZABLE N,P DURING PRIMARY ROOT REMOBILIZATION
!     DEPENDS ON ROOT NON-STRUCTURAL C:N:P
!
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     WFR=constraint by O2 consumption on all root processes
!     RCER(ielmc),RCER(ielmn),RCER(ielmp)=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC1=fraction of primary root C to be remobilized
!
  IF(-RCO2XM.GT.0.0_r8)THEN
    IF(-RCO2XM.LT.RTWT1E(ielmc,N,NR,NZ)*RCCC)THEN
      SNCRM=-RCO2XM
    ELSE
      SNCRM=AZMAX1(RTWT1E(ielmc,N,NR,NZ)*RCCC)
    ENDIF
  ELSE
    SNCRM=0._r8
  ENDIF
  IF(-RCO2X.GT.0.0_r8)THEN
    IF(-RCO2X.LT.RTWT1E(ielmc,N,NR,NZ)*RCCC)THEN
      SNCR=-RCO2X
    ELSE
      SNCR=AZMAX1(RTWT1E(ielmc,N,NR,NZ)*RCCC)*WFR(N,L,NZ)
    ENDIF
  ELSE
    SNCR=0._r8
  ENDIF
  IF(SNCR.GT.0.0_r8.AND.RTWT1E(ielmc,N,NR,NZ).GT.ZEROP(NZ))THEN
    RCER(ielmc)=RCCC*RTWT1E(ielmc,N,NR,NZ)
    RCER(ielmn)=RTWT1E(ielmn,N,NR,NZ)*(RCCN+(1.0_r8-RCCN)*RCER(ielmc)/RTWT1E(ielmc,N,NR,NZ))
    RCER(ielmp)=RTWT1E(ielmp,N,NR,NZ)*(RCCP+(1.0_r8-RCCP)*RCER(ielmc)/RTWT1E(ielmc,N,NR,NZ))
    IF(RCER(ielmc).GT.ZEROP(NZ))THEN
      FSNC1=AZMAX1(AMIN1(1.0_r8,SNCR/RCER(ielmc)))
    ELSE
      FSNC1=1.0_r8
    ENDIF
  ELSE
    RCER(1:npelms)=0._r8
    FSNC1=0._r8
  ENDIF
!
!     PRIMARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=literfall C,N,P
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RCER(ielmc),RCER(ielmn),RCER(ielmp)=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
  D6355: DO M=1,jsken
    DO NE=1,npelms    
      ESNC(NE,M,k_woody_litr,L,NZ)=ESNC(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
        *FSNC1*(RTWT1E(NE,N,NR,NZ)-RCER(NE))*FWODRE(NE,k_woody_litr)

      ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
        *FSNC1*(RTWT1E(NE,N,NR,NZ)-RCER(NE))*FWODRE(NE,k_fine_litr)
    ENDDO    
  ENDDO D6355
  
  end associate
  end subroutine PrimRootRemobilization

!------------------------------------------------------------------------------------------
  subroutine PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,GRTWTG,GRTWTLE,GRTLGL,RTLGZ,WTRTZ)
  implicit none
  integer, intent(in) :: L,L1,N,NR,NZ
  real(r8), intent(in):: WFNR,FRTN,GRTWTG,GRTWTLE(npelms)
  real(r8), intent(inout) :: RTLGZ,WTRTZ
  real(r8), intent(out):: GRTLGL
  real(r8) :: FGROL,FGROZ
  integer :: NE
  REAL(R8) :: XFRE(npelms)
! begin_execution
  associate(                             &
    WTRT1E    =>  plt_biom%WTRT1E      , &
    RTWT1E    =>  plt_biom%RTWT1E      , &
    WSRTL     =>  plt_biom%WSRTL       , &
    WTRTD     =>  plt_biom%WTRTD       , &
    EPOOLR    =>  plt_biom%EPOOLR      , &
    ZEROP     =>  plt_biom%ZEROP       , &
    CNWS      =>  plt_allom%CNWS       , &
    CPWS      =>  plt_allom%CPWS       , &
    FWODRE    =>  plt_allom%FWODRE     , &
    CDPTHZ    =>  plt_site%CDPTHZ      , &
    DLYR3     =>  plt_site%DLYR3       , &
    NJ        =>  plt_site%NJ          , &
    PP        =>  plt_site%PP          , &
    PSIRG     =>  plt_ew%PSIRG         , &
    PSIRO     =>  plt_ew%PSIRO         , &
    PSIRT     =>  plt_ew%PSIRT         , &
    k_woody_litr=> pltpar%k_woody_litr , &
    k_fine_litr=> pltpar%k_fine_litr   , &
    RTLG1X    =>  plt_morph%RTLG1X     , &
    RRAD1     =>  plt_morph%RRAD1      , &
    RTLG1     =>  plt_morph%RTLG1      , &
    RTDP1     =>  plt_morph%RTDP1      , &
    NG        =>  plt_morph%NG         , &
    NINR      =>  plt_morph%NINR       , &
    SDPTH     =>  plt_morph%SDPTH        &
  )
!     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     GRTLGL=primary root length extension
!     GRTWTG=primary root C growth ltd by O2
!     RTLG1X=specific primary root length from startq.f
!     PP=PFT population
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     GRTWTLE(ielmc),GRTWTLE(ielmn),GRTWTLE(ielmp)=net primary root C,N,P growth
!     RTDP1=primary root depth from soil surface
!     SDPTH=seeding depth
!     FSNC1=fraction of primary root C to be remobilized
!     RTLG1=primary root length
!     GRTWTLE(ielmc),GRTWTLE(ielmn),GRTWTLE(ielmp)=net root C,N,P growth
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     DLYR=soil layer thickness
!
  IF(GRTWTLE(ielmc).LT.0.0.AND.RTWT1E(ielmc,N,NR,NZ).GT.ZEROP(NZ))THEN
    GRTLGL=GRTWTG*RTLG1X(N,NZ)/PP(NZ)*WFNR*FWODRE(ielmc,k_fine_litr) &
      +GRTWTLE(ielmc)*(RTDP1(N,NR,NZ)-SDPTH(NZ))/RTWT1E(ielmc,N,NR,NZ)
  ELSE
    GRTLGL=GRTWTG*RTLG1X(N,NZ)/PP(NZ)*WFNR*FWODRE(ielmc,k_fine_litr)
  ENDIF
  IF(L.LT.NJ)THEN
    GRTLGL=AMIN1(DLYR3(L1),GRTLGL)
  ENDIF
!
!     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
!     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
!     BOUNDARY OF CURRENT LAYER
!
!     GRTLGL=primary root length extension
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!
  IF(GRTLGL.GT.ZEROP(NZ).AND.L.LT.NJ)THEN
    FGROL=AZMAX1(AMIN1(1.0_r8,(CDPTHZ(L)-RTDP1(N,NR,NZ))/GRTLGL))
    IF(FGROL.LT.1.0_r8)FGROL=0._r8
    FGROZ=AZMAX1(1.0_r8-FGROL)
  ELSE
    FGROL=1.0_r8
    FGROZ=0._r8
  ENDIF
!
!     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
!     AND AXIS NUMBER
!
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     GRTWTLE(ielmc),GRTWTLE(ielmn),GRTWTLE(ielmp)=net root C,N,P growth
!     GRTLGL=primary root length extension
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTLG1=primary root length
!
  RTDP1(N,NR,NZ)=RTDP1(N,NR,NZ)+GRTLGL

  DO NE=1,npelms
    RTWT1E(NE,N,NR,NZ)=RTWT1E(NE,N,NR,NZ)+GRTWTLE(NE)
    WTRT1E(NE,N,L,NR,NZ)=WTRT1E(NE,N,L,NR,NZ)+GRTWTLE(NE)*FGROL
  ENDDO
  WSRTL(N,L,NZ)=WSRTL(N,L,NZ)+AMIN1(CNWS(NZ)*WTRT1E(ielmn,N,L,NR,NZ) &
    ,CPWS(NZ)*WTRT1E(ielmp,N,L,NR,NZ))
  RTLG1(N,L,NR,NZ)=RTLG1(N,L,NR,NZ)+GRTLGL*FGROL
!
!     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
!     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
!     OF CURRENT SOIL LAYER
!
!     FGROZ=fraction of GRTLGL in next lower soil layer
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     GRTWTLE(ielmc),GRTWTLE(ielmn),GRTWTLE(ielmp)=net root C,N,P growth
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     WTRTD=root C mass
!     RTLG1=primary root length
!     GRTLGL=primary root length extension
!     FRTN=fraction of primary root sink strength in axis
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     PSIRT,PSIRG,PSIRO=root total,turgor,osmotic water potential
!     NINR=deepest root layer
!
  IF(FGROZ.GT.0.0)THEN
    DO NE=1,npelms
      WTRT1E(NE,N,L1,NR,NZ)=WTRT1E(NE,N,L1,NR,NZ)+GRTWTLE(NE)*FGROZ
    ENDDO
    WSRTL(N,L1,NZ)=WSRTL(N,L1,NZ)+AMIN1(CNWS(NZ)*WTRT1E(ielmn,N,L1,NR,NZ) &
      ,CPWS(NZ)*WTRT1E(ielmp,N,L1,NR,NZ))
    WTRTD(N,L1,NZ)=WTRTD(N,L1,NZ)+WTRT1E(ielmc,N,L1,NR,NZ)
    RTLG1(N,L1,NR,NZ)=RTLG1(N,L1,NR,NZ)+GRTLGL*FGROZ
    RRAD1(N,L1,NZ)=RRAD1(N,L,NZ)
    RTLGZ=RTLGZ+RTLG1(N,L1,NR,NZ)
    WTRTZ=WTRTZ+WTRT1E(ielmc,N,L1,NR,NZ)

    DO NE=1,npelms
      XFRE(NE)=FRTN*EPOOLR(NE,N,L,NZ)
      EPOOLR(NE,N,L,NZ)=EPOOLR(NE,N,L,NZ)-XFRE(NE)
      EPOOLR(NE,N,L1,NZ)=EPOOLR(NE,N,L1,NZ)+XFRE(NE)
    ENDDO
    PSIRT(N,L1,NZ)=PSIRT(N,L,NZ)
    PSIRO(N,L1,NZ)=PSIRO(N,L,NZ)
    PSIRG(N,L1,NZ)=PSIRG(N,L,NZ)
    NINR(NR,NZ)=MAX(NG(NZ),L+1)
  ENDIF
  end associate
  end subroutine PrimRootExtension
!------------------------------------------------------------------------------------------

  subroutine WithdrawPrimRoot(L,NR,NZ,N,RTDP1X,RLNT,RTSK1,RTSK2,XRTN1)
  implicit none
  integer, intent(in) :: L,NR,NZ,N
  real(r8), intent(in):: RTDP1X
  real(r8), INTENT(IN) :: RLNT(2,JZ1)
  real(r8),intent(in) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10),XRTN1
  integer :: LL,NN,NE,NTG
  real(r8) :: XFRD,XFRW,FRTN
  real(r8) :: XFRE(npelms)

! begin_execution
  associate(                             &
    WTRT1E   =>  plt_biom%WTRT1E   , &
    WTRT2E   =>  plt_biom%WTRT2E   , &
    EPOOLR   =>  plt_biom%EPOOLR   , &
    WSRTL    =>  plt_biom%WSRTL    , &
    WTRTD    =>  plt_biom%WTRTD    , &
    WTNDLE   =>  plt_biom%WTNDLE   , &
    ZEROP    =>  plt_biom%ZEROP    , &
    EPOOLN   =>  plt_biom%EPOOLN   , &
    RFGas_root    =>  plt_bgcr%RFGas_root    , &
    CDPTHZ   =>  plt_site%CDPTHZ   , &
    ZEROS2   =>  plt_site%ZEROS2   , &
    DLYR3    =>  plt_site%DLYR3    , &
    trcg_rootml     =>  plt_rbgc%trcg_rootml , &
    trcs_rootml => plt_rbgc%trcs_rootml, &
    VOLX     =>  plt_soilchem%VOLX , &
    RTN1     =>  plt_morph%RTN1    , &
    RTNL     =>  plt_morph%RTNL    , &
    RTLG1    =>  plt_morph%RTLG1   , &
    RTN2     =>  plt_morph%RTN2    , &
    RTDP1    =>  plt_morph%RTDP1   , &
    NG       =>  plt_morph%NG      , &
    INTYP    =>  plt_morph%INTYP   , &
    MY       =>  plt_morph%MY      , &
    SDPTH    =>  plt_morph%SDPTH   , &
    NINR     =>  plt_morph%NINR      &
  )
!     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
!     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
!     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
!     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
!
!     NINR=deepest root layer
!     VOLX=soil layer volume excluding macropore, rocks
!     RTDP1X=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!     FRTN=fraction of primary root sink strength in axis
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTLG1=primary root length
!     WSRTL=root protein C mass
!     WTRTD=root C mass
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root

  D5115: DO LL=L,NG(NZ)+1,-1
    IF(VOLX(LL-1).GT.ZEROS2.AND.(RTDP1X.LT.CDPTHZ(LL-1).OR.RTDP1X.LT.SDPTH(NZ)))THEN
      IF(RLNT(N,LL).GT.ZEROP(NZ))THEN
        FRTN=(RTSK1(N,LL,NR)+RTSK2(N,LL,NR))/RLNT(N,LL)
      ELSE
        FRTN=1.0_r8
      ENDIF
      D5110: DO NN=1,MY(NZ)
        RTLG1(NN,LL-1,NR,NZ)=RTLG1(NN,LL-1,NR,NZ)+RTLG1(NN,LL,NR,NZ)
        RTLG1(NN,LL,NR,NZ)=0._r8
        DO NE=1,npelms
          WTRT1E(NE,NN,LL-1,NR,NZ)=WTRT1E(NE,NN,LL-1,NR,NZ)+WTRT1E(NE,NN,LL,NR,NZ)

          WTRT2E(NE,NN,LL-1,NR,NZ)=WTRT2E(NE,NN,LL-1,NR,NZ)+WTRT2E(NE,NN,LL,NR,NZ)

          WTRT1E(NE,NN,LL,NR,NZ)=0._r8

          WTRT2E(NE,NN,LL,NR,NZ)=0._r8

          XFRE(NE)=FRTN*EPOOLR(NE,NN,LL,NZ)

          EPOOLR(NE,NN,LL,NZ)=EPOOLR(NE,NN,LL,NZ)-XFRE(NE)

          EPOOLR(NE,NN,LL-1,NZ)=EPOOLR(NE,NN,LL-1,NZ)+XFRE(NE)
        ENDDO
        XFRW=FRTN*WSRTL(NN,L,NZ)
        XFRD=FRTN*WTRTD(NN,LL,NZ)

        WSRTL(NN,LL,NZ)=WSRTL(NN,LL,NZ)-XFRW
        WTRTD(NN,LL,NZ)=WTRTD(NN,LL,NZ)-XFRD
        WSRTL(NN,LL-1,NZ)=WSRTL(NN,LL-1,NZ)+XFRW
        WTRTD(NN,LL-1,NZ)=WTRTD(NN,LL-1,NZ)+XFRD
!
!     WITHDRAW GASES IN PRIMARY ROOTS
!
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
!     FRTN=fraction of primary root sink strength in axis
!
        DO NTG=idg_beg,idg_end-1
          RFGas_root(NTG,NZ)=RFGas_root(NTG,NZ)-FRTN &
            *(trcg_rootml(idg_CO2,NN,LL,NZ)+trcs_rootml(idg_CO2,NN,LL,NZ))
          trcg_rootml(NTG,NN,LL,NZ)=(1.0_r8-FRTN)*trcg_rootml(NTG,NN,LL,NZ)
          trcs_rootml(NTG,NN,LL,NZ)=(1.0_r8-FRTN)*trcs_rootml(NTG,NN,LL,NZ)
        ENDDO
      ENDDO D5110
!
!     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
!
!     RTN2,RTNL=number of secondary root axes
!     RTN1=number of primary root axes
!     RTLG1=primary root length
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!
      RTNL(N,LL,NZ)=RTNL(N,LL,NZ)-RTN2(N,LL,NR,NZ)
      RTNL(N,LL-1,NZ)=RTNL(N,LL-1,NZ)+RTN2(N,LL,NR,NZ)
      RTN2(N,LL,NR,NZ)=0._r8
      RTN1(N,LL,NZ)=RTN1(N,LL,NZ)-XRTN1
      IF(LL-1.GT.NG(NZ))THEN
        RTLG1(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CDPTHZ(LL-1)-RTDP1(N,NR,NZ))
      ELSE
        RTLG1(N,LL-1,NR,NZ)=DLYR3(LL-1)-(CDPTHZ(LL-1)-RTDP1(N,NR,NZ)) &
          -(SDPTH(NZ)-CDPTHZ(LL-2))
      ENDIF
!
!     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     FRTN=fraction of primary root sink strength in axis
!     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
!
      IF(INTYP(NZ).GE.1.AND.INTYP(NZ).LE.3)THEN
        DO NE=1,npelms
          XFRE(NE)=FRTN*WTNDLE(NE,LL,NZ)
          WTNDLE(NE,LL,NZ)=WTNDLE(NE,LL,NZ)-XFRE(NE)
          WTNDLE(NE,LL-1,NZ)=WTNDLE(NE,LL-1,NZ)+XFRE(NE)

          XFRE(NE)=FRTN*EPOOLN(NE,LL,NZ)
          EPOOLN(NE,LL,NZ)=EPOOLN(NE,LL,NZ)-XFRE(NE)
          EPOOLN(NE,LL-1,NZ)=EPOOLN(NE,LL-1,NZ)+XFRE(NE)
        ENDDO
      ENDIF
      NINR(NR,NZ)=MAX(NG(NZ),LL-1)
    ELSE
      EXIT
    ENDIF
  ENDDO D5115
  end associate
  end subroutine WithdrawPrimRoot
!------------------------------------------------------------------------------------------

  subroutine NonstructlBiomTransfer(I,J,NZ,PTRT,RLNT,RTSK1,RTSK2,RTNT,IFLGZ)
  implicit none
  integer, intent(in) :: I,J,NZ,IFLGZ
  real(r8), intent(in):: PTRT
  real(r8), INTENT(IN) :: RLNT(2,JZ1)
  real(r8),intent(in) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10)
  real(r8),intent(in) :: RTNT(2)
  integer :: L,NB,N,NR,NE
  real(r8) :: ZPOOLS,ZPOOLT
  real(r8) :: ZPOOLB
  real(r8) :: ZPOOLD,EPOOLD
  real(r8) :: XFRPX,XFRCX,XFRNX
  real(r8) :: WTLSBZ(JC1)
  real(r8) :: CPOOLZ(JC1),ZPOOLZ(JC1),PPOOLZ(JC1)
  REAL(R8) :: FWTR(JZ1),FWTB(JP1)
  real(r8) :: CPOOLT
  real(r8) :: CNL,CPL,CPOOLD
  real(r8) :: CPOOLB,CPOOLS
  real(r8) :: FWTC
  real(r8) :: FWTS
  real(r8) :: PPOOLB
  real(r8) :: PPOOLD
  real(r8) :: PPOOLT
  real(r8) :: PTSHTR
  real(r8) :: PPOOLS
  real(r8) :: WTPLTT
  real(r8) :: WTRTLX
  real(r8) :: WTRTD1
  real(r8) :: WTSTKT
  real(r8) :: WTRSVT,WTRSNT,WTRSPT
  real(r8) :: WTRSVD,WTRSND,WTRSPD
  real(r8) :: WTRTD2,WTLSBX,WTLSBB
  real(r8) :: WTRTLR
  real(r8) :: XFRE(npelms)
!     begin_execution
  associate(                                &
    FWODBE     =>   plt_allom%FWODBE  , &
    FWODRE     =>   plt_allom%FWODRE  , &
    CEPOLR     =>   plt_biom%CEPOLR   , &
    RTWT1E     =>   plt_biom%RTWT1E   , &
    WTRT1E     =>   plt_biom%WTRT1E   , &
    WTRT2E     =>   plt_biom%WTRT2E   , &
    WTRTL      =>   plt_biom%WTRTL    , &
    EPOOLR     =>   plt_biom%EPOOLR   , &
    WTLSB      =>   plt_biom%WTLSB    , &
    EPOOL      =>   plt_biom%EPOOL    , &
    WTRTD      =>   plt_biom%WTRTD    , &
    WVSTKB     =>   plt_biom%WVSTKB   , &
    WTRSVBE    =>   plt_biom%WTRSVBE  , &
    WTRVE      =>   plt_biom%WTRVE    , &
    WTLS       =>   plt_biom%WTLS     , &
    WTRTE      =>   plt_biom%WTRTE    , &
    ZEROL      =>   plt_biom%ZEROL    , &
    ZEROP      =>   plt_biom%ZEROP    , &
    IBTYP      =>   plt_pheno%IBTYP   , &
    IGTYP      =>   plt_pheno%IGTYP   , &
    IDTHB      =>   plt_pheno%IDTHB   , &
    ISTYP      =>   plt_pheno%ISTYP   , &
    PTSHT      =>   plt_pheno%PTSHT   , &
    IDAY       =>   plt_pheno%IDAY    , &
    ATRP       =>   plt_pheno%ATRP    , &
    RCO2A      =>   plt_rbgc%RCO2A    , &
    TCO2T      =>   plt_bgcr%TCO2T    , &
    RECO       =>   plt_bgcr%RECO     , &
    TRAU       =>   plt_bgcr%TRAU     , &
    NU         =>   plt_site%NU       , &
    ZERO       =>   plt_site%ZERO     , &
    k_woody_litr=> pltpar%k_woody_litr,&
    k_fine_litr=> pltpar%k_fine_litr, &
    NIX        =>   plt_morph%NIX     , &
    NINR       =>   plt_morph%NINR    , &
    RRAD2      =>   plt_morph%RRAD2   , &
    NI         =>   plt_morph%NI      , &
    MY         =>   plt_morph%MY      , &
    NRT        =>   plt_morph%NRT     , &
    NBR        =>   plt_morph%NBR       &
  )
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH LEAVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!     WHEN SEASONAL STORAGE C IS NOT BEING MOBILIZED
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     ATRP=hourly leafout counter
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     WTLSB=leaf+petiole mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!
  IF(NBR(NZ).GT.1)THEN
    WTPLTT=0._r8
    CPOOLT=0._r8
    ZPOOLT=0._r8
    PPOOLT=0._r8
    D300: DO NB=1,NBR(NZ)
      IF(IDTHB(NB,NZ).EQ.ialive)THEN
        IF(ATRP(NB,NZ).GT.ATRPX(ISTYP(NZ)))THEN
          WTLSBZ(NB)=AZMAX1(WTLSB(NB,NZ))
          CPOOLZ(NB)=AZMAX1(EPOOL(ielmc,NB,NZ))
          ZPOOLZ(NB)=AZMAX1(EPOOL(ielmn,NB,NZ))
          PPOOLZ(NB)=AZMAX1(EPOOL(ielmp,NB,NZ))
          WTPLTT=WTPLTT+WTLSBZ(NB)
          CPOOLT=CPOOLT+CPOOLZ(NB)
          ZPOOLT=ZPOOLT+ZPOOLZ(NB)
          PPOOLT=PPOOLT+PPOOLZ(NB)
        ENDIF
      ENDIF
    ENDDO D300
    D305: DO NB=1,NBR(NZ)
      IF(IDTHB(NB,NZ).EQ.ialive)THEN
        IF(ATRP(NB,NZ).GT.ATRPX(ISTYP(NZ)))THEN
          IF(WTPLTT.GT.ZEROP(NZ).AND.CPOOLT.GT.ZEROP(NZ))THEN
            CPOOLD=CPOOLT*WTLSBZ(NB)-CPOOLZ(NB)*WTPLTT
            ZPOOLD=ZPOOLT*CPOOLZ(NB)-ZPOOLZ(NB)*CPOOLT
            PPOOLD=PPOOLT*CPOOLZ(NB)-PPOOLZ(NB)*CPOOLT
            XFRE(ielmc)=0.01_r8*CPOOLD/WTPLTT
            XFRE(ielmn)=0.01_r8*ZPOOLD/CPOOLT
            XFRE(ielmp)=0.01_r8*PPOOLD/CPOOLT
            DO NE=1,npelms
              EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)+XFRE(NE)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO D305
  ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH STALK RESERVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     WVSTKB=stalk sapwood mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     IDAY(7,=start of grain filling and setting max seed size
!
  IF(NBR(NZ).GT.1)THEN
    WTSTKT=0._r8
    WTRSVT=0._r8
    WTRSNT=0._r8
    WTRSPT=0._r8
    D330: DO NB=1,NBR(NZ)
      IF(IDTHB(NB,NZ).EQ.ialive)THEN
        IF(IDAY(7,NB,NZ).NE.0)THEN
          WTSTKT=WTSTKT+WVSTKB(NB,NZ)
          WTRSVT=WTRSVT+WTRSVBE(ielmc,NB,NZ)
          WTRSNT=WTRSNT+WTRSVBE(ielmn,NB,NZ)
          WTRSPT=WTRSPT+WTRSVBE(ielmp,NB,NZ)
        ENDIF
      ENDIF
    ENDDO D330
    IF(WTSTKT.GT.ZEROP(NZ).AND.WTRSVT.GT.ZEROP(NZ))THEN
      D335: DO NB=1,NBR(NZ)
        IF(IDTHB(NB,NZ).EQ.ialive)THEN
          IF(IDAY(7,NB,NZ).NE.0)THEN
            WTRSVD=WTRSVT*WVSTKB(NB,NZ)-WTRSVBE(ielmc,NB,NZ)*WTSTKT
            XFRE(ielmc)=0.1_r8*WTRSVD/WTSTKT
            WTRSVBE(ielmc,NB,NZ)=WTRSVBE(ielmc,NB,NZ)+XFRE(ielmc)
            WTRSND=WTRSNT*WTRSVBE(ielmc,NB,NZ)-WTRSVBE(ielmn,NB,NZ)*WTRSVT
            XFRE(ielmn)=0.1_r8*WTRSND/WTRSVT
            WTRSVBE(ielmn,NB,NZ)=WTRSVBE(ielmn,NB,NZ)+XFRE(ielmn)
            WTRSPD=WTRSPT*WTRSVBE(ielmc,NB,NZ)-WTRSVBE(ielmp,NB,NZ)*WTRSVT
            XFRE(ielmp)=0.1_r8*WTRSPD/WTRSVT
            WTRSVBE(ielmp,NB,NZ)=WTRSVBE(ielmp,NB,NZ)+XFRE(ielmp)
          ENDIF
        ENDIF
      ENDDO D335
    ENDIF
  ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BWTWEEN ROOT AND MYCORRHIZAE
!     IN EACH ROOTED SOIL LAYER FROM NON-STRUCTURAL C,N,P
!     CONCENTRATION DIFFERENCES
!
!     MY=mycorrhizal:1=no,2=yes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in 1:root,2:mycorrhizae
!     WTRTD=1:root,2:mycorrhizal C mass
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     FMYC=rate constant for root-mycorrhizal C,N,P exchange (h-1)
!
  IF(MY(NZ).EQ.mycorarbu)THEN
    D425: DO L=NU,NIX(NZ)
      IF(EPOOLR(ielmc,ipltroot,L,NZ).GT.ZEROP(NZ).AND.WTRTD(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
!root
        WTRTD1=WTRTD(ipltroot,L,NZ)
        WTRTD2=AMIN1(WTRTD(ipltroot,L,NZ),AMAX1(FSNK &
          *WTRTD(ipltroot,L,NZ),WTRTD(imycorrhz,L,NZ)))
        WTPLTT=WTRTD1+WTRTD2
        IF(WTPLTT.GT.ZEROP(NZ))THEN
          CPOOLD=(EPOOLR(ielmc,ipltroot,L,NZ)*WTRTD2-EPOOLR(ielmc,imycorrhz,L,NZ)*WTRTD1)/WTPLTT
          XFRE(ielmc)=FMYC*CPOOLD
          EPOOLR(ielmc,ipltroot,L,NZ)=EPOOLR(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
          EPOOLR(ielmc,imycorrhz,L,NZ)=EPOOLR(ielmc,imycorrhz,L,NZ)+XFRE(ielmc)
          CPOOLT=EPOOLR(ielmc,ipltroot,L,NZ)+EPOOLR(ielmc,imycorrhz,L,NZ)

          IF(CPOOLT.GT.ZEROP(NZ))THEN
            DO NE=2,npelms
              EPOOLD=(EPOOLR(NE,ipltroot,L,NZ)*EPOOLR(ielmc,imycorrhz,L,NZ) &
                -EPOOLR(NE,imycorrhz,L,NZ)*EPOOLR(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(NE)=FMYC*EPOOLD
              EPOOLR(NE,ipltroot,L,NZ)=EPOOLR(NE,ipltroot,L,NZ)-XFRE(NE)
              EPOOLR(NE,imycorrhz,L,NZ)=EPOOLR(NE,imycorrhz,L,NZ)+XFRE(NE)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO D425
  ENDIF
!
!     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
!     IN PERENNIALS
!
  IF(IFLGZ.EQ.1.AND.ISTYP(NZ).NE.iplt_annual)THEN
    D5545: DO N=1,MY(NZ)
      D5550: DO L=NU,NI(NZ)
        IF(CEPOLR(ielmc,N,L,NZ).GT.ZERO)THEN
          CNL=CEPOLR(ielmc,N,L,NZ)/(CEPOLR(ielmc,N,L,NZ)+CEPOLR(ielmn,N,L,NZ)/CNKI)
          CPL=CEPOLR(ielmc,N,L,NZ)/(CEPOLR(ielmc,N,L,NZ)+CEPOLR(ielmp,N,L,NZ)/CPKI)
        ELSE
          CNL=0._r8
          CPL=0._r8
        ENDIF
        XFRCX=FXFR(IBTYP(NZ))*AZMAX1(EPOOLR(ielmc,N,L,NZ))
        XFRNX=FXFR(IBTYP(NZ))*AZMAX1(EPOOLR(ielmn,N,L,NZ))*(1.0_r8+CNL)
        XFRPX=FXFR(IBTYP(NZ))*AZMAX1(EPOOLR(ielmp,N,L,NZ))*(1.0_r8+CPL)
        XFRE(ielmc)=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
        XFRE(ielmn)=AMIN1(XFRNX,XFRE(ielmc)*CNMX,XFRPX*CNMX/CPMN*0.5_r8)
        XFRE(ielmp)=AMIN1(XFRPX,XFRE(ielmc)*CPMX,XFRNX*CPMX/CNMN*0.5_r8)
        DO NE=1,npelms
          EPOOLR(NE,N,L,NZ)=EPOOLR(NE,N,L,NZ)-XFRE(NE)
          WTRVE(NE,NZ)=WTRVE(NE,NZ)+XFRE(NE)
        ENDDO

      ENDDO D5550
    ENDDO D5545
  ENDIF
!
!     ROOT AND NODULE TOTALS
!
!     WTRTL,WTRTD=active,actual root C mass
!     WTRT1,WTRT2=primary,secondary root C mass in soil layer
!     TCO2T=total PFT respiration
!     RCO2A=total root respiration
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
  D5445: DO N=1,MY(NZ)
    D5450: DO L=NU,NI(NZ)
      WTRTL(N,L,NZ)=0._r8
      WTRTD(N,L,NZ)=0._r8
      D5460: DO NR=1,NRT(NZ)
        WTRTL(N,L,NZ)=WTRTL(N,L,NZ)+WTRT2E(ielmc,N,L,NR,NZ)
        WTRTD(N,L,NZ)=WTRTD(N,L,NZ)+WTRT2E(ielmc,N,L,NR,NZ)+WTRT1E(ielmc,N,L,NR,NZ)
      ENDDO D5460
      TCO2T(NZ)=TCO2T(NZ)+RCO2A(N,L,NZ)
      RECO=RECO+RCO2A(N,L,NZ)
      TRAU=TRAU+RCO2A(N,L,NZ)
    ENDDO D5450

    DO  NR=1,NRT(NZ)
      WTRTL(N,NINR(NR,NZ),NZ)=WTRTL(N,NINR(NR,NZ),NZ)+RTWT1E(ielmc,N,NR,NZ)
    ENDDO
  ENDDO D5445
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
!
!     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
!     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     WTLS,WTRT=total PFT leaf+petiole,root C mass
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     RLNT,RTNT=root layer,root system sink strength
!
!     IF(ISTYP(NZ).EQ.iplt_preanu)THEN
  IF(WTLS(NZ).GT.ZEROP(NZ))THEN
    FWTC=AMIN1(1.0_r8,0.667_r8*WTRTE(ielmc,NZ)/WTLS(NZ))
  ELSE
    FWTC=1.0_r8
  ENDIF
  IF(WTRTE(ielmc,NZ).GT.ZEROP(NZ))THEN
    FWTS=AMIN1(1.0_r8,WTLS(NZ)/(0.667_r8*WTRTE(ielmc,NZ)))
  ELSE
    FWTS=1.0_r8
  ENDIF
!     ELSE
!     FWTC=1.0_r8
!     FWTS=1.0_r8
!     ENDIF
  D290: DO L=NU,NI(NZ)
    IF(RTNT(1).GT.ZEROP(NZ))THEN
      FWTR(L)=AZMAX1(RLNT(ipltroot,L)/RTNT(1))
    ELSE
      FWTR(L)=1.0_r8
    ENDIF
  ENDDO D290
!     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
!     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
!
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!
  WTLS(NZ)=0._r8
  D309: DO NB=1,NBR(NZ)
    WTLS(NZ)=WTLS(NZ)+WTLSB(NB,NZ)
  ENDDO D309
!
!     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
!     OF TOTAL SINK STRENGTH OF THE CANOPY
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAY(8,=end date for setting final seed number
!     FWTB=branch sink weighting factor
!     PTSHT=rate constant for equilibrating shoot-root nonstructural C concn from PFT file
!     PTRT=allocation to leaf+petiole used to modify PTSHT in annuals
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     FWODB=C woody fraction in branch:0=woody,1=non-woody
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  D310: DO NB=1,NBR(NZ)
    IF(IDTHB(NB,NZ).EQ.ialive)THEN
      IF(WTLS(NZ).GT.ZEROP(NZ))THEN
        FWTB(NB)=AZMAX1(WTLSB(NB,NZ)/WTLS(NZ))
      ELSE
        FWTB(NB)=1.0_r8
      ENDIF
      IF(ISTYP(NZ).EQ.iplt_annual)THEN
        PTSHTR=PTSHT(NZ)*PTRT**0.167_r8
      ELSE
        PTSHTR=PTSHT(NZ)
      ENDIF
      D415: DO L=NU,NI(NZ)
        WTLSBX=WTLSB(NB,NZ)*FWODBE(ielmc,k_fine_litr)*FWTR(L)*FWTC
        WTRTLX=WTRTL(ipltroot,L,NZ)*FWODRE(ielmc,k_fine_litr)*FWTB(NB)*FWTS
        WTLSBB=AZMAX1(WTLSBX,FSNK*WTRTLX)
        WTRTLR=AZMAX1(WTRTLX,FSNK*WTLSBX)
        WTPLTT=WTLSBB+WTRTLR
        IF(WTPLTT.GT.ZEROP(NZ))THEN
          CPOOLB=AZMAX1(EPOOL(ielmc,NB,NZ)*FWTR(L))
          CPOOLS=AZMAX1(EPOOLR(ielmc,ipltroot,L,NZ)*FWTB(NB))
          CPOOLD=(CPOOLB*WTRTLR-CPOOLS*WTLSBB)/WTPLTT
          XFRE(ielmc)=PTSHTR*CPOOLD
          CPOOLT=CPOOLS+CPOOLB
          IF(CPOOLT.GT.ZEROP(NZ))THEN
            ZPOOLB=AZMAX1(EPOOL(ielmn,NB,NZ)*FWTR(L))
            ZPOOLS=AZMAX1(EPOOLR(ielmn,ipltroot,L,NZ)*FWTB(NB))
            ZPOOLD=(ZPOOLB*CPOOLS-ZPOOLS*CPOOLB)/CPOOLT
            XFRE(ielmn)=PTSHTR*ZPOOLD
            PPOOLB=AZMAX1(EPOOL(ielmp,NB,NZ)*FWTR(L))
            PPOOLS=AZMAX1(EPOOLR(ielmp,ipltroot,L,NZ)*FWTB(NB))
            PPOOLD=(PPOOLB*CPOOLS-PPOOLS*CPOOLB)/CPOOLT
            XFRE(ielmp)=PTSHTR*PPOOLD
          ELSE
            XFRE(ielmn)=0._r8
            XFRE(ielmp)=0._r8
          ENDIF
          DO NE=1,npelms
            EPOOL(NE,NB,NZ)=EPOOL(NE,NB,NZ)-XFRE(NE)
            EPOOLR(NE,ipltroot,L,NZ)=EPOOLR(NE,ipltroot,L,NZ)+XFRE(NE)
          ENDDO

        ENDIF
      ENDDO D415
    ENDIF
  ENDDO D310
  end associate
  end subroutine NonstructlBiomTransfer

!------------------------------------------------------------------------------------------
  subroutine SummarizeRootSink(NZ,XRTN1,RLNT,RTSK1,RTSK2,RTNT)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: XRTN1
  real(r8),INTENT(OUT) :: RLNT(2,JZ1)
  real(r8),intent(out) :: RTSK1(2,JZ1,10),RTSK2(2,JZ1,10)
  real(r8),INTENT(OUT) :: RTNT(2)
  integer :: N,L,K,NR,NE
  REAL(R8) :: RTDPL(10,JZ1)
  real(r8) :: CUPRL,CUPRO,CUPRC
  real(r8) :: RTDPP,RTDPS,RTSKP
  real(r8) :: RTSKS

  associate(                              &
    EPOOLR     =>   plt_biom%EPOOLR    , &
    ZEROP      =>   plt_biom%ZEROP     , &
    IGTYP      =>   plt_pheno%IGTYP    , &
    VOLX       =>   plt_soilchem%VOLX  , &
    ZEROS2     =>   plt_site%ZEROS2    , &
    NU         =>   plt_site%NU        , &
    CDPTHZ     =>   plt_site%CDPTHZ    , &
    ZERO       =>   plt_site%ZERO      , &
    DLYR3      =>   plt_site%DLYR3     , &
    RUPH1B     =>   plt_rbgc%RUPH1B    , &
    RUPH2P     =>   plt_rbgc%RUPH2P    , &
    RUPH2B     =>   plt_rbgc%RUPH2B    , &
    RUOH2B     =>   plt_rbgc%RUOH2B    , &
    RUOH1P     =>   plt_rbgc%RUOH1P    , &
    RUCH1B     =>   plt_rbgc%RUCH1B    , &
    RUCH2B     =>   plt_rbgc%RUCH2B    , &
    RUONH4     =>   plt_rbgc%RUONH4    , &
    RUCH1P     =>   plt_rbgc%RUCH1P    , &
    RUCH2P     =>   plt_rbgc%RUCH2P    , &
    RUCNOB     =>   plt_rbgc%RUCNOB    , &
    RUCNO3     =>   plt_rbgc%RUCNO3    , &
    RDFOME     =>   plt_rbgc%RDFOME    , &
    RUPNOB     =>   plt_rbgc%RUPNOB    , &
    RUOH2P     =>   plt_rbgc%RUOH2P    , &
    RUONOB     =>   plt_rbgc%RUONOB    , &
    RUONHB     =>   plt_rbgc%RUONHB    , &
    RUONO3     =>   plt_rbgc%RUONO3    , &
    RUCNHB     =>   plt_rbgc%RUCNHB    , &
    RUPNH4     =>   plt_rbgc%RUPNH4    , &
    RUPNO3     =>   plt_rbgc%RUPNO3    , &
    RUPNHB     =>   plt_rbgc%RUPNHB    , &
    RCO2N      =>   plt_rbgc%RCO2N     , &
    RUPH1P     =>   plt_rbgc%RUPH1P    , &
    RUCNH4     =>   plt_rbgc%RUCNH4    , &
    RUOH1B     =>   plt_rbgc%RUOH1B    , &
    RCO2M      =>   plt_rbgc%RCO2M     , &
    RCO2A      =>   plt_rbgc%RCO2A     , &
    HTSTZ      =>   plt_morph%HTSTZ    , &
    MY         =>   plt_morph%MY       , &
    RRAD1      =>   plt_morph%RRAD1    , &
    RTDP1      =>   plt_morph%RTDP1    , &
    HTCTL      =>   plt_morph%HTCTL    , &
    RRAD2      =>   plt_morph%RRAD2    , &
    RTN2       =>   plt_morph%RTN2     , &
    RTLGA      =>   plt_morph%RTLGA    , &
    SDPTH      =>   plt_morph%SDPTH    , &
    NI         =>   plt_morph%NI       , &
    NRT        =>   plt_morph%NRT        &
  )

!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER

  RLNT=0._R8
  RTSK1=0._r8
  RTSK2=0._r8
  RTNT=0._r8
  D4995: DO N=1,MY(NZ)
    D4990: DO L=NU,NI(NZ)
!
!     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
!     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
!
!     VOLX=soil layer volume excluding macropore, rocks
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!     RUONH4,RUONHB,RUON03,RUONOB=uptake from non-band,band of NH4,NO3 unlimited by O2
!     RUOH2P,RUOH2B,RUOH1P,RUOH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by O2
!     RUCNH4,RUCNHB,RUCN03,RUCNOB=uptake from non-band,band of NH4,NO3 unlimited by nonstructural C
!     RUCH2P,RUCH2B,RUCH1P,RUCH1B=uptake from non-band,band of H2PO4,HPO4 unlimited by nonstructural C
!
      IF(VOLX(L).GT.ZEROS2)THEN
        CUPRL=0.86_r8*(RUPNH4(N,L,NZ)+RUPNHB(N,L,NZ) &
          +RUPNO3(N,L,NZ)+RUPNOB(N,L,NZ)+RUPH2P(N,L,NZ) &
          +RUPH2B(N,L,NZ)+RUPH1P(N,L,NZ)+RUPH1B(N,L,NZ))
        CUPRO=0.86_r8*(RUONH4(N,L,NZ)+RUONHB(N,L,NZ) &
          +RUONO3(N,L,NZ)+RUONOB(N,L,NZ)+RUOH2P(N,L,NZ) &
          +RUOH2B(N,L,NZ)+RUOH1P(N,L,NZ)+RUOH1B(N,L,NZ))
        CUPRC=0.86_r8*(RUCNH4(N,L,NZ)+RUCNHB(N,L,NZ) &
          +RUCNO3(N,L,NZ)+RUCNOB(N,L,NZ)+RUCH2P(N,L,NZ) &
          +RUCH2B(N,L,NZ)+RUCH1P(N,L,NZ)+RUCH1B(N,L,NZ))
!
!     ACCUMULATE RESPIRATION IN FLUX ARRAYS
!
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     CPOOLR=non-structural C mass in root
!
        RCO2M(N,L,NZ)=RCO2M(N,L,NZ)+CUPRO
        RCO2N(N,L,NZ)=RCO2N(N,L,NZ)+CUPRC
        RCO2A(N,L,NZ)=RCO2A(N,L,NZ)-CUPRL
        EPOOLR(ielmc,N,L,NZ)=EPOOLR(ielmc,N,L,NZ)-CUPRL
!
!     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
!     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
!     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!
        D195: DO K=1,jcplx
          DO NE=1,npelms
            EPOOLR(NE,N,L,NZ)=EPOOLR(NE,N,L,NZ)+RDFOME(NE,N,K,L,NZ)
          ENDDO
        ENDDO D195
        EPOOLR(ielmn,N,L,NZ)=EPOOLR(ielmn,N,L,NZ)+RUPNH4(N,L,NZ)+RUPNHB(N,L,NZ) &
          +RUPNO3(N,L,NZ)+RUPNOB(N,L,NZ)
        EPOOLR(ielmp,N,L,NZ)=EPOOLR(ielmp,N,L,NZ)+RUPH2P(N,L,NZ)+RUPH2B(N,L,NZ) &
          +RUPH1P(N,L,NZ)+RUPH1B(N,L,NZ)
!
!     GROWTH OF EACH ROOT AXIS
!
        D4985: DO NR=1,NRT(NZ)
!
!     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
!
!     RTDP1=primary root depth from soil surface
!     RTDPP=primary root depth from canopy
!     CDPTHZ=depth from soil surface to layer bottom
!     HTSTZ=canopy height for water uptake
!     RTSK=relative primary root sink strength
!     RTSK1=primary root sink strength
!     XRTN1=number of primary root axes
!     RRAD1,RRAD2=primary, secondary root radius
!     RTNT,RLNT=total root sink strength
!
          IF(N.EQ.1)THEN
            IF(RTDP1(N,NR,NZ).GT.CDPTHZ(L-1))THEN
              IF(RTDP1(N,NR,NZ).LE.CDPTHZ(L))THEN
                RTDPP=RTDP1(N,NR,NZ)+HTSTZ(NZ)
                RTSK1(N,L,NR)=RTSK(IGTYP(NZ))*XRTN1*RRAD1(N,L,NZ)**2._r8/RTDPP
                RTNT(N)=RTNT(N)+RTSK1(N,L,NR)
                RLNT(N,L)=RLNT(N,L)+RTSK1(N,L,NR)
              ENDIF
            ENDIF
          ENDIF
!
!     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
!     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
!
!     RTDPL=depth of primary root axis in layer
!     RTDP1=primary root depth from soil surface
!     CDPTHZ=depth from soil surface to layer bottom
!     RTDPX=distance behind growing point for secondary roots
!     DLYR=layer thickness
!     SDPTH=seeding depth
!     HTCTL=hypocotyledon height
!     HTSTZ=canopy height for water uptake
!     RTDPS=secondary root depth from canopy
!     RTSKP,RTSKS=primary,secondary root sink strength
!     RTN2=number of secondary root axes
!     RTSK2=total secondary root sink strength
!     RTLGA=average secondary root length
!     RTNT,RLNT=total root sink strength
!
          IF(N.EQ.1)THEN
            RTDPL(NR,L)=AZMAX1(RTDP1(ipltroot,NR,NZ)-CDPTHZ(L-1)-RTDPX)
            RTDPL(NR,L)=AZMAX1(AMIN1(DLYR3(L),RTDPL(NR,L)) &
              -AZMAX1(SDPTH(NZ)-CDPTHZ(L-1)-HTCTL(NZ)))
            RTDPS=AMAX1(SDPTH(NZ),CDPTHZ(L-1))+0.5*RTDPL(NR,L)+HTSTZ(NZ)
            IF(RTDPS.GT.ZERO)THEN
              RTSKP=XRTN1*RRAD1(N,L,NZ)**2._r8/RTDPS
              RTSKS=safe_adb(RTN2(N,L,NR,NZ)*RRAD2(N,L,NZ)**2._r8,RTLGA(N,L,NZ))
              IF(RTSKP+RTSKS.GT.ZEROP(NZ))THEN
                RTSK2(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
              ELSE
                RTSK2(N,L,NR)=0._r8
              ENDIF
            ELSE
              RTSK2(N,L,NR)=0._r8
            ENDIF
          ELSE
            RTSK2(N,L,NR)=safe_adb(RTN2(N,L,NR,NZ)*RRAD2(N,L,NZ)**2._r8,RTLGA(N,L,NZ))
          ENDIF
          RTNT(N)=RTNT(N)+RTSK2(N,L,NR)
          RLNT(N,L)=RLNT(N,L)+RTSK2(N,L,NR)
        ENDDO D4985
      ENDIF
    ENDDO D4990
  ENDDO D4995
  end associate
  end subroutine SummarizeRootSink

end module RootMod
