module RootMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod  , only : test_aeqb,safe_adb
  use EcosimConst
  use GrosubPars
  use PlantAPIData
implicit none
  private

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
    SDLGs1     =>   plt_morph%SDLGs1      , &
    SDVLs1     =>   plt_morph%SDVLs1      , &
    SDARs1     =>   plt_morph%SDARs1      , &
    NRTs1      =>  plt_morph%NRTs1          &
  )
!     ROOT GROWTH
!
  call RootBiochemistry(I,J,NZ,ICHK1,IDTHRN,NRX,TFN6,CNRTW,CPRTW,&
      XRTN1,WFNGR,RLNT,RTSK1,RTSK2,RTNT)
!
!     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
!     GERMINATION)
!
  RTLGPs1(1,NGs1(NZ),NZ)=RTLGPs1(1,NGs1(NZ),NZ)+SDLGs1(NZ)
  IF(DLYR3s1(NGs1(NZ)).GT.ZEROs1)THEN
    RTDNPs1(1,NGs1(NZ),NZ)=RTLGPs1(1,NGs1(NZ),NZ)/DLYR3s1(NGs1(NZ))
  ELSE
    RTDNPs1(1,NGs1(NZ),NZ)=0._r8
  ENDIF
  RTVL=RTVLPs1(1,NGs1(NZ),NZ)+RTVLWs1(1,NGs1(NZ),NZ)+SDVLs1(NZ)*PPs1(NZ)
  RTVLPs1(1,NGs1(NZ),NZ)=PORTs1(1,NZ)*RTVL
  RTVLWs1(1,NGs1(NZ),NZ)=(1.0_r8-PORTs1(1,NZ))*RTVL
  RTARPs1(1,NGs1(NZ),NZ)=RTARPs1(1,NGs1(NZ),NZ)+SDARs1(NZ)
  IF(IDTHRN.EQ.NRTs1(NZ).OR.(WTRVCs1(NZ).LE.ZEROLs1(NZ).AND.ISTYPs1(NZ).NE.0))THEN
    IDTHRs1(NZ)=1
    IDTHPs1(NZ)=1
  ENDIF
!
!     ROOT N2 FIXATION (RHIZOBIA)
!
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
  integer :: LL,LZ,L1,L,K,lx,M,NR,N
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
  associate(                                &
    NIs1       =>   plt_morph%NIs1        , &
    SDLGs1     =>   plt_morph%SDLGs1      , &
    NIXs1      =>   plt_morph%NIXs1         &
  )

  NIXs1(NZ)=NGs1(NZ)
  IDTHRN=0
!
  call SummarizeRootSink(NZ,XRTN1,RLNT,RTSK1,RTSK2,RTNT)
!
!     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
!
  DO 5010 N=1,MYs1(NZ)
    DO 5000 L=NUs1,NIs1(NZ)
!
!     IDENTIFY NEXT LOWER ROOT LAYER
!
!     VOLX=soil layer volume excluding macropore, rocks
!
      IF(VOLXs1(L).GT.ZEROS2s1)THEN
!     WRITE(*,4994)'5004',I,J,NZ,N,L,NIs1(NZ)
!    2,NLs1,VOLXs1(L),CDPTHZs1(L-1)
        DO 5003 LZ=L+1,NLs1
!     WRITE(*,4994)'5003',I,J,NZ,N,L,LZ
!    2,LZ,VOLXs1(L),CDPTHZs1(LZ)
          IF(VOLXs1(LZ).GT.ZEROS2s1 &
            .OR.LZ.EQ.NLs1)THEN
            L1=LZ
            EXIT
          ENDIF
5003    CONTINUE
!     WRITE(*,4994)'5005',I,J,NZ,N,L,LZ
!    2,L1,VOLXs1(L),CDPTHZs1(L1)
!4994  FORMAT(A8,7I4,12E12.4)
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
        RSCS2=RSCSs1(L)*RRAD2s1(N,L,NZ)/1.0E-03
        WFNR=AMIN1(1.0,AMAX1(0.0,PSIRGs1(N,L,NZ)-PSILM-RSCS2))
        IF(IGTYPs1(NZ).EQ.0)THEN
          WFNGR(N,L)=EXP(0.05*PSIRTs1(N,L,NZ))
          WFNRG=WFNR**0.10
        ELSE
          WFNGR(N,L)=EXP(0.10*PSIRTs1(N,L,NZ))
          WFNRG=WFNR**0.25
        ENDIF
        DMRTD=1.0_r8-DMRTs1(NZ)
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
        IF(L.LE.NIXs1(NZ))THEN
          IF(WTRTLs1(N,L,NZ).GT.ZEROPs1(NZ) &
            .AND.WTRTs1(NZ).GT.ZEROPs1(NZ) &
            .AND.WTRVCs1(NZ).LT.XFRX*WTRTs1(NZ))THEN
            FWTRT=WTRTLs1(N,L,NZ)/WTRTs1(NZ)
            WTRTLX=WTRTLs1(N,L,NZ)
            WTRTTX=WTRTs1(NZ)*FWTRT
            WTRTTT=WTRTLX+WTRTTX
            CPOOLX=AMAX1(0.0,CPOOLRs1(N,L,NZ))
            WTRVCX=AMAX1(0.0,WTRVCs1(NZ)*FWTRT)
            CPOOLD=(WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
            XFRC=AMIN1(0.0,XFRY*CPOOLD)
            CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+XFRC
            WTRVCs1(NZ)=WTRVCs1(NZ)-XFRC
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
        IF(N.EQ.1)THEN
          RTLGZ=RTLGZ*FWODRs1(1)
          RTLGL=RTLGL*FWODRs1(1)
        ENDIF
        RTLGX=RTLGZ*PPs1(NZ)
        RTLGT=RTLGL+RTLGX
        WTRTT=WTRTX+WTRTZ
        IF(RTLGT.GT.ZEROPs1(NZ).AND.WTRTT.GT.ZEROPs1(NZ) &
          .AND.PPs1(NZ).GT.ZEROPs1(NZ))THEN
          RTLGPs1(N,L,NZ)=RTLGT/PPs1(NZ)
          IF(DLYR3s1(L).GT.ZEROs1)THEN
            RTDNPs1(N,L,NZ)=RTLGPs1(N,L,NZ)/DLYR3s1(L)
          ELSE
            RTDNPs1(N,L,NZ)=0._r8
          ENDIF
          RTVL=AMAX1(RTAR1Xs1(N,NZ)*RTLGX+RTAR2Xs1(N,NZ)*RTLGL &
            ,WTRTT*DMVLs1(N,NZ)*PSIRGs1(N,L,NZ))
          RTVLPs1(N,L,NZ)=PORTs1(N,NZ)*RTVL
          RTVLWs1(N,L,NZ)=(1.0_r8-PORTs1(N,NZ))*RTVL
          RRAD1s1(N,L,NZ)=AMAX1(RRAD1Xs1(N,NZ) &
            ,(1.0+PSIRTs1(N,L,NZ)/EMODR)*RRAD1Ms1(N,NZ))
          RRAD2s1(N,L,NZ)=AMAX1(RRAD2Xs1(N,NZ) &
            ,(1.0+PSIRTs1(N,L,NZ)/EMODR)*RRAD2Ms1(N,NZ))
          RTAR=6.283*RRAD1s1(N,L,NZ)*RTLGX &
            +6.283*RRAD2s1(N,L,NZ)*RTLGL
          IF(RTNLs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
            RTLGAs1(N,L,NZ)=AMAX1(RTLGAX,RTLGL/RTNLs1(N,L,NZ))
          ELSE
            RTLGAs1(N,L,NZ)=RTLGAX
          ENDIF
          RTARPs1(N,L,NZ)=RTAR/PPs1(NZ)
!     IF(N.EQ.1)THEN
!     RTARPs1(N,L,NZ)=RTARPs1(N,L,NZ)*RTLGAX/RTLGAs1(N,L,NZ)
!     ENDIF
        ELSE
          RTLGPs1(N,L,NZ)=0._r8
          RTDNPs1(N,L,NZ)=0._r8
          RTVLPs1(N,L,NZ)=0._r8
          RTVLWs1(N,L,NZ)=0._r8
          RRAD1s1(N,L,NZ)=RRAD1Ms1(N,NZ)
          RRAD2s1(N,L,NZ)=RRAD2Ms1(N,NZ)
          RTARPs1(N,L,NZ)=0._r8
          RTLGAs1(N,L,NZ)=RTLGAX
          RCO2Zs1(NZ)=RCO2Zs1(NZ)-(CO2As1(N,L,NZ)+CO2Ps1(N,L,NZ))
          ROXYZs1(NZ)=ROXYZs1(NZ)-(OXYAs1(N,L,NZ)+OXYPs1(N,L,NZ))
          RCH4Zs1(NZ)=RCH4Zs1(NZ)-(CH4As1(N,L,NZ)+CH4Ps1(N,L,NZ))
          RN2OZs1(NZ)=RN2OZs1(NZ)-(Z2OAs1(N,L,NZ)+Z2OPs1(N,L,NZ))
          RNH3Zs1(NZ)=RNH3Zs1(NZ)-(ZH3As1(N,L,NZ)+ZH3Ps1(N,L,NZ))
          RH2GZs1(NZ)=RH2GZs1(NZ)-(H2GAs1(N,L,NZ)+H2GPs1(N,L,NZ))
          CO2As1(N,L,NZ)=0._r8
          OXYAs1(N,L,NZ)=0._r8
          CH4As1(N,L,NZ)=0._r8
          Z2OAs1(N,L,NZ)=0._r8
          ZH3As1(N,L,NZ)=0._r8
          H2GAs1(N,L,NZ)=0._r8
          CO2Ps1(N,L,NZ)=0._r8
          OXYPs1(N,L,NZ)=0._r8
          CH4Ps1(N,L,NZ)=0._r8
          Z2OPs1(N,L,NZ)=0._r8
          ZH3Ps1(N,L,NZ)=0._r8
          H2GPs1(N,L,NZ)=0._r8
        ENDIF
      ENDIF
5000  CONTINUE
5010  CONTINUE
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
  real(r8) :: GRTWTP,GRTWTN,GRTWTL
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
  real(r8) :: RCCR,RCZR,RCPR
  real(r8) :: RTN2X,RTN2Y
  real(r8) :: RTDP1X,RSCS1
  REAL(R8) :: SNCR,SNCRM
  real(r8) :: TFRCO2
  real(r8) :: RCCC,RCCN,RCCP
  integer :: NR,M,LL,LX

!begin_execution
  associate(                              &
    NIXs1     =>  plt_morph%NIXs1       , &
    NRTs1     =>  plt_morph%NRTs1       , &
    SDPTHs1   =>  plt_morph%SDPTHs1     , &
    NINRs1    =>  plt_morph%NINRs1      , &
    NB1s1     =>  plt_morph%NB1s1       , &
    FDBKXs1   =>  plt_photo%FDBKXs1       &
  )
  RTLGL=0._r8
  RTLGZ=0._r8
  WTRTX=0._r8
  WTRTZ=0._r8
  DO 5050 NR=1,NRTs1(NZ)
!
!     SECONDARY ROOT EXTENSION
!
    IF(L.LE.NINRs1(NR,NZ).AND.NRX(N,NR).EQ.0)THEN
!
!     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
!     TO CURRENT AXIS
!
!     RTSK2=total secondary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of secondary root sink strength in axis
!
      IF(RLNT(N,L).GT.ZEROPs1(NZ))THEN
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
      IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
        CNPG=AMIN1(CZPOLRs1(N,L,NZ)/(CZPOLRs1(N,L,NZ) &
          +CCPOLRs1(N,L,NZ)*CNKI),CPPOLRs1(N,L,NZ) &
          /(CPPOLRs1(N,L,NZ)+CCPOLRs1(N,L,NZ)*CPKI))
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
      RMNCR=AMAX1(0.0,RMPLT*WTRT2Ns1(N,L,NR,NZ))*TFN6(L)
      IF(IGTYPs1(NZ).EQ.0.OR.IWTYPs1(NZ).EQ.2)THEN
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
      RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLRs1(N,L,NZ) &
        *TFN4s1(L,NZ))*CNPG*FDBKXs1(NB1s1(NZ),NZ)*WFNGR(N,L)
!
!     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCO2R=respiration from non-structural C limited by O2
!     WFR=constraint by O2 consumption on all root processes
!     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
!     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
!     WFNRG=respiration function of root water potential
!
      RCO2R=RCO2RM*WFRs1(N,L,NZ)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AMAX1(0.0,RCO2XM)*WFNRG
      RCO2Y=AMAX1(0.0,RCO2X)*WFNRG
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
      ZPOOLB=AMAX1(0.0,ZPOOLRs1(N,L,NZ))
      PPOOLB=AMAX1(0.0,PPOOLRs1(N,L,NZ))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTSs1(NZ) &
        ,PPOOLB*DMRTR/CPRTSs1(NZ))
      IF(RCO2YM.GT.0.0)THEN
        RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
        RCO2GM=0._r8
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
        RCO2G=AMIN1(RCO2Y,FNP*WFRs1(N,L,NZ))
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
      GRTWGM=CGRORM*DMRTs1(NZ)
      GRTWTG=CGROR*DMRTs1(NZ)
      ZADD2M=AMAX1(0.0,GRTWGM*CNRTW)
      ZADD2=AMAX1(0.0,AMIN1(FRTN*ZPOOLRs1(N,L,NZ),GRTWTG*CNRTW))
      PADD2=AMAX1(0.0,AMIN1(FRTN*PPOOLRs1(N,L,NZ),GRTWTG*CPRTW))
      CNRDM=AMAX1(0.0,1.70*ZADD2M)
      CNRDA=AMAX1(0.0,1.70*ZADD2)
!
!     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAYs1(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
      IF(IDAYs1(1,NB1s1(NZ),NZ).NE.0.AND.CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
        CCC=AMAX1(0.0_r8,AMIN1(1.0_r8,safe_adb(CZPOLRs1(N,L,NZ),CZPOLRs1(N,L,NZ) &
          +CCPOLRs1(N,L,NZ)*CNKI) &
          ,safe_adb(CPPOLRs1(N,L,NZ),CPPOLRs1(N,L,NZ) &
          +CCPOLRs1(N,L,NZ)*CPKI)))
        CNC=AMAX1(0.0_r8,AMIN1(1.0_r8 &
          ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)/CNKI)))
        CPC=AMAX1(0.0_r8,AMIN1(1.0_r8 &
          ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)/CPKI)))
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
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC2=fraction of secondary root C to be remobilized
!
      IF(-RCO2XM.GT.0.0)THEN
        IF(-RCO2XM.LT.WTRT2s1(N,L,NR,NZ)*RCCC)THEN
          SNCRM=-RCO2XM
        ELSE
          SNCRM=AMAX1(0.0,WTRT2s1(N,L,NR,NZ)*RCCC)
        ENDIF
      ELSE
        SNCRM=0._r8
      ENDIF
      IF(-RCO2X.GT.0.0)THEN
        IF(-RCO2X.LT.WTRT2s1(N,L,NR,NZ)*RCCC)THEN
          SNCR=-RCO2X
        ELSE
          SNCR=AMAX1(0.0,WTRT2s1(N,L,NR,NZ)*RCCC)*WFRs1(N,L,NZ)
        ENDIF
      ELSE
        SNCR=0._r8
      ENDIF
      IF(SNCR.GT.0.0.AND.WTRT2s1(N,L,NR,NZ).GT.ZEROPs1(NZ))THEN
        RCCR=RCCC*WTRT2s1(N,L,NR,NZ)
        RCZR=WTRT2Ns1(N,L,NR,NZ)*(RCCN+(1.0_r8-RCCN)*RCCR/WTRT2s1(N,L,NR,NZ))
        RCPR=WTRT2Ps1(N,L,NR,NZ)*(RCCP+(1.0_r8-RCCP)*RCCR/WTRT2s1(N,L,NR,NZ))
        IF(RCCR.GT.ZEROPs1(NZ))THEN
          FSNC2=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
        ELSE
          FSNC2=1.0_r8
        ENDIF
      ELSE
        RCCR=0._r8
        RCZR=0._r8
        RCPR=0._r8
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
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
      DO 6350 M=1,4
        CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
          *FSNC2*(WTRT2s1(N,L,NR,NZ)-RCCR)*FWODRs1(0)
        ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
          *FSNC2*(WTRT2Ns1(N,L,NR,NZ)-RCZR)*FWODRNs1(0)
        PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
          *FSNC2*(WTRT2Ps1(N,L,NR,NZ)-RCPR)*FWODRPs1(0)
        CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
          *FSNC2*(WTRT2s1(N,L,NR,NZ)-RCCR)*FWODRs1(1)
        ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
          *FSNC2*(WTRT2Ns1(N,L,NR,NZ)-RCZR)*FWODRNs1(1)
        PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
          *FSNC2*(WTRT2Ps1(N,L,NR,NZ)-RCPR)*FWODRPs1(1)
6350  CONTINUE
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
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-AMIN1(RMNCR,RCO2R) &
        -CGROR-CNRDA-SNCR+FSNC2*RCCR
      ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-ZADD2+FSNC2*RCZR
      PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-PADD2+FSNC2*RCPR
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
      RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)+RCO2TM
      RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)+RCO2T
      RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)-RCO2T
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
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
!
      GRTLGL=GRTWTG*RTLG2Xs1(N,NZ)*WFNR*FWODRs1(1) &
        -FSNC2*RTLG2s1(N,L,NR,NZ)
      GRTWTL=GRTWTG-FSNC2*WTRT2s1(N,L,NR,NZ)
      GRTWTN=ZADD2-FSNC2*WTRT2Ns1(N,L,NR,NZ)
      GRTWTP=PADD2-FSNC2*WTRT2Ps1(N,L,NR,NZ)
!
!     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
!     AND AXIS NUMBER
!
!     RTLG2=secondary root length
!     GRTLGL=secondary root length extension
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTFQ=root branching frequency from PFT file
!     RTN2,RTNL=number of secondary root axes
!
      RTLG2s1(N,L,NR,NZ)=RTLG2s1(N,L,NR,NZ)+GRTLGL
      WTRT2s1(N,L,NR,NZ)=WTRT2s1(N,L,NR,NZ)+GRTWTL
      WTRT2Ns1(N,L,NR,NZ)=WTRT2Ns1(N,L,NR,NZ)+GRTWTN
      WTRT2Ps1(N,L,NR,NZ)=WTRT2Ps1(N,L,NR,NZ)+GRTWTP
      WSRTLs1(N,L,NZ)=WSRTLs1(N,L,NZ) &
        +AMIN1(CNWSs1(NZ)*WTRT2Ns1(N,L,NR,NZ) &
        ,CPWSs1(NZ)*WTRT2Ps1(N,L,NR,NZ))
      RTLGL=RTLGL+RTLG2s1(N,L,NR,NZ)
      WTRTX=WTRTX+WTRT2s1(N,L,NR,NZ)
      RTN2X=RTFQs1(NZ)*XRTN1
      RTN2Y=RTFQs1(NZ)*RTN2X
      RTN2s1(N,L,NR,NZ)=(RTN2X+RTN2Y)*DLYR3s1(L)
      RTNLs1(N,L,NZ)=RTNLs1(N,L,NZ)+RTN2s1(N,L,NR,NZ)
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
      IF(N.EQ.1)THEN
        IF(BKDSs1(L).GT.ZEROs1)THEN
          RTDP1X=RTDP1s1(N,NR,NZ)-CDPTHZs1(0)
        ELSE
          RTDP1X=RTDP1s1(N,NR,NZ)
        ENDIF
        IF(RTDP1X.GT.CDPTHZs1(L-1).AND.ICHK1(N,NR).EQ.0)THEN
            RTN1s1(N,L,NZ)=RTN1s1(N,L,NZ)+XRTN1
            IF(RTDP1X.LE.CDPTHZs1(L).OR.L.EQ.NJs1)THEN
              ICHK1(N,NR)=1
!
!     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
!     ATTRIBUTED TO CURRENT AXIS
!
!     RTSK1=primary root sink strength
!     RLNT=total root sink strength
!     FRTN=fraction of primary root sink strength in axis
!
              IF(RLNT(N,L).GT.ZEROPs1(NZ))THEN
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
              RSCS1=RSCSs1(L)*RRAD1s1(N,L,NZ)/1.0E-03
              WFNR=AMIN1(1.0,AMAX1(0.0,PSIRGs1(N,L,NZ)-PSILM-RSCS1))
              IF(IGTYPs1(NZ).EQ.0)THEN
                WFNRG=WFNR**0.10
              ELSE
                WFNRG=WFNR**0.25
              ENDIF
!
!     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
!     NON-STRUCTURAL C:N:P
!
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNPG=N,P constraint on growth respiration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
              IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
                CNPG=AMIN1(CZPOLRs1(N,L,NZ)/(CZPOLRs1(N,L,NZ) &
                  +CCPOLRs1(N,L,NZ)*CNKI),CPPOLRs1(N,L,NZ) &
                  /(CPPOLRs1(N,L,NZ)+CCPOLRs1(N,L,NZ)*CPKI))
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
              RMNCR=AMAX1(0.0,RMPLT*RTWT1Ns1(N,NR,NZ))*TFN6(L)
              IF(IGTYPs1(NZ).EQ.0.OR.IWTYPs1(NZ).EQ.2)THEN
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
              RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLRs1(N,L,NZ) &
                *TFN4s1(L,NZ))*CNPG*FDBKXs1(NB1s1(NZ),NZ)*WFNGR(N,L)
              IF(RTDP1X.GE.CDPTHZs1(NJs1))THEN
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
              RCO2R=RCO2RM*WFRs1(N,L,NZ)
              RCO2XM=RCO2RM-RMNCR
              RCO2X=RCO2R-RMNCR
              RCO2YM=AMAX1(0.0,RCO2XM)*WFNRG
              RCO2Y=AMAX1(0.0,RCO2X)*WFNRG
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
              ZPOOLB=AMAX1(0.0,ZPOOLRs1(N,L,NZ))
              PPOOLB=AMAX1(0.0,PPOOLRs1(N,L,NZ))
              FNP=AMIN1(ZPOOLB*DMRTR/CNRTSs1(NZ),PPOOLB*DMRTR/CPRTSs1(NZ))
              IF(RCO2YM.GT.0.0)THEN
                RCO2GM=AMIN1(RCO2YM,FNP)
              ELSE
                RCO2GM=0._r8
              ENDIF
              IF(RCO2Y.GT.0.0)THEN
                RCO2G=AMIN1(RCO2Y,FNP*WFRs1(N,L,NZ))
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
              GRTWGM=CGRORM*DMRTs1(NZ)
              GRTWTG=CGROR*DMRTs1(NZ)
              ZADD1M=AMAX1(0.0,GRTWGM*CNRTW)
              ZADD1=AMAX1(0.0,AMIN1(FRTN*ZPOOLRs1(N,L,NZ),GRTWTG*CNRTW))
              PADD1=AMAX1(0.0,AMIN1(FRTN*PPOOLRs1(N,L,NZ),GRTWTG*CPRTW))
              CNRDM=AMAX1(0.0,1.70*ZADD1M)
              CNRDA=AMAX1(0.0,1.70*ZADD1)

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
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!
              CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-AMIN1(RMNCR,RCO2R)-CGROR-CNRDA-SNCR+FSNC1*RCCR
              ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-ZADD1+FSNC1*RCZR
              PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-PADD1+FSNC1*RCPR
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
              IF(RTDP1s1(N,NR,NZ).GT.CDPTHZs1(NGs1(NZ)))THEN
                TFRCO2=0._r8
                DO 5100 LL=NGs1(NZ),NINRs1(NR,NZ)
                  IF(LL.LT.NINRs1(NR,NZ))THEN
                    FRCO2=AMIN1(1.0,RTLG1s1(N,LL,NR,NZ) &
                      /(RTDP1s1(N,NR,NZ)-SDPTHs1(NZ)))
                  ELSE
                    FRCO2=1.0_r8-TFRCO2
                  ENDIF
                  TFRCO2=TFRCO2+FRCO2
                  RCO2Ms1(N,LL,NZ)=RCO2Ms1(N,LL,NZ)+RCO2TM*FRCO2
                  RCO2Ns1(N,LL,NZ)=RCO2Ns1(N,LL,NZ)+RCO2T*FRCO2
                  RCO2As1(N,LL,NZ)=RCO2As1(N,LL,NZ)-RCO2T*FRCO2
5100            CONTINUE
              ELSE
                RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)+RCO2TM
                RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)+RCO2T
                RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)-RCO2T
              ENDIF
!
!     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
!     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
!     HAVE DISAPPEARED
!
!     GRTWTG=primary root C growth ltd by O2
!     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
!     FSNC1=fraction of primary root C to be remobilized
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
!     RTLG2=secondary root length
!
              GRTWTL=GRTWTG-FSNC1*RTWT1s1(N,NR,NZ)
              GRTWTN=ZADD1-FSNC1*RTWT1Ns1(N,NR,NZ)
              GRTWTP=PADD1-FSNC1*RTWT1Ps1(N,NR,NZ)
              IF(GRTWTL.LT.0.0)THEN
                LX=MAX(1,L-1)
                DO 5105 LL=L,LX,-1
                  GRTWTM=GRTWTL
                  IF(GRTWTL.LT.0.0)THEN
                    IF(GRTWTL.GT.-WTRT2s1(N,LL,NR,NZ))THEN
                      RTLG2s1(N,LL,NR,NZ)=RTLG2s1(N,LL,NR,NZ)+GRTWTL &
                        *RTLG2s1(N,LL,NR,NZ)/WTRT2s1(N,LL,NR,NZ)
                      WTRT2s1(N,LL,NR,NZ)=WTRT2s1(N,LL,NR,NZ)+GRTWTL
                      GRTWTL=0._r8
                    ELSE
                      GRTWTL=GRTWTL+WTRT2s1(N,LL,NR,NZ)
                      RTLG2s1(N,LL,NR,NZ)=0._r8
                      WTRT2s1(N,LL,NR,NZ)=0._r8
                    ENDIF
                  ENDIF
                  IF(GRTWTN.LT.0.0)THEN
                    IF(GRTWTN.GT.-WTRT2Ns1(N,LL,NR,NZ))THEN
                      WTRT2Ns1(N,LL,NR,NZ)=WTRT2Ns1(N,LL,NR,NZ)+GRTWTN
                      GRTWTN=0._r8
                    ELSE
                      GRTWTN=GRTWTN+WTRT2Ns1(N,LL,NR,NZ)
                      WTRT2Ns1(N,LL,NR,NZ)=0._r8
                    ENDIF
                  ENDIF
                  IF(GRTWTP.LT.0.0)THEN
                    IF(GRTWTP.GT.-WTRT2Ps1(N,LL,NR,NZ))THEN
                      WTRT2Ps1(N,LL,NR,NZ)=WTRT2Ps1(N,LL,NR,NZ)+GRTWTP
                      GRTWTP=0._r8
                    ELSE
                      GRTWTP=GRTWTP+WTRT2Ps1(N,LL,NR,NZ)
                      WTRT2Ps1(N,LL,NR,NZ)=0._r8
                    ENDIF
                  ENDIF
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
                  IF(GRTWTM.LT.0.0)THEN
                    IF(WTRT2s1(1,LL,NR,NZ).GT.ZEROPs1(NZ))THEN
                      FSNCM=AMIN1(1.0,ABS(GRTWTM)/WTRT2s1(1,LL,NR,NZ))
                    ELSE
                      FSNCM=1.0_r8
                    ENDIF
                    IF(WTRTLs1(1,LL,NZ).GT.ZEROPs1(NZ))THEN
                      FSNCP=AMIN1(1.0,ABS(GRTWTM)/WTRTLs1(1,LL,NZ))
                    ELSE
                      FSNCP=1.0_r8
                    ENDIF
                    DO 6450 M=1,4
                      CSNCs1(M,0,LL,NZ)=CSNCs1(M,0,LL,NZ)+CFOPCs1(5,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2s1(2,LL,NR,NZ))*FWODRs1(0)
                      ZSNCs1(M,0,LL,NZ)=ZSNCs1(M,0,LL,NZ)+CFOPNs1(5,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ns1(2,LL,NR,NZ))*FWODRNs1(0)
                      PSNCs1(M,0,LL,NZ)=PSNCs1(M,0,LL,NZ)+CFOPPs1(5,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ps1(2,LL,NR,NZ))*FWODRPs1(0)
                      CSNCs1(M,1,LL,NZ)=CSNCs1(M,1,LL,NZ)+CFOPCs1(4,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2s1(2,LL,NR,NZ))*FWODRs1(1)
                      ZSNCs1(M,1,LL,NZ)=ZSNCs1(M,1,LL,NZ)+CFOPNs1(4,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ns1(2,LL,NR,NZ))*FWODRNs1(1)
                      PSNCs1(M,1,LL,NZ)=PSNCs1(M,1,LL,NZ)+CFOPPs1(4,M,NZ) &
                        *FSNCM*AMAX1(0.0,WTRT2Ps1(2,LL,NR,NZ))*FWODRPs1(1)
                      CSNCs1(M,1,LL,NZ)=CSNCs1(M,1,LL,NZ)+CFOPCs1(0,M,NZ) &
                        *FSNCP*AMAX1(0.0,CPOOLRs1(2,LL,NZ))
                      ZSNCs1(M,1,LL,NZ)=ZSNCs1(M,1,LL,NZ)+CFOPNs1(0,M,NZ) &
                        *FSNCP*AMAX1(0.0,ZPOOLRs1(2,LL,NZ))
                      PSNCs1(M,1,LL,NZ)=PSNCs1(M,1,LL,NZ)+CFOPPs1(0,M,NZ) &
                        *FSNCP*AMAX1(0.0,PPOOLRs1(2,LL,NZ))
6450                CONTINUE
                    RTLG2s1(2,LL,NR,NZ)=AMAX1(0.0,RTLG2s1(2,LL,NR,NZ))*(1.0_r8-FSNCM)
                    WTRT2s1(2,LL,NR,NZ)=AMAX1(0.0,WTRT2s1(2,LL,NR,NZ))*(1.0_r8-FSNCM)
                    WTRT2Ns1(2,LL,NR,NZ)=AMAX1(0.0,WTRT2Ns1(2,LL,NR,NZ))*(1.0_r8-FSNCM)
                    WTRT2Ps1(2,LL,NR,NZ)=AMAX1(0.0,WTRT2Ps1(2,LL,NR,NZ))*(1.0_r8-FSNCM)
                    CPOOLRs1(2,LL,NZ)=AMAX1(0.0,CPOOLRs1(2,LL,NZ))*(1.0_r8-FSNCP)
                    ZPOOLRs1(2,LL,NZ)=AMAX1(0.0,ZPOOLRs1(2,LL,NZ))*(1.0_r8-FSNCP)
                    PPOOLRs1(2,LL,NZ)=AMAX1(0.0,PPOOLRs1(2,LL,NZ))*(1.0_r8-FSNCP)
                  ENDIF
5105            CONTINUE
              ENDIF
!
              call PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,GRTWTG,GRTWTL,GRTWTN,GRTWTP,&
                GRTLGL,RTLGZ,WTRTZ)
            ENDIF
!
!
            IF(L.EQ.NINRs1(NR,NZ))THEN
              call WithdrawPrimRoot(L,NR,NZ,N,RTDP1X,RLNT,RTSK1,RTSK2,XRTN1)
            ENDIF

!
!     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
!
            IF(WTRT1s1(N,L,NR,NZ).LT.0.0)THEN
              CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+WTRT1s1(N,L,NR,NZ)
              WTRT1s1(N,L,NR,NZ)=0._r8
            ENDIF
            IF(WTRT2s1(N,L,NR,NZ).LT.0.0)THEN
              CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+WTRT2s1(N,L,NR,NZ)
              WTRT2s1(N,L,NR,NZ)=0._r8
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
            RTLGZ=RTLGZ+RTLG1s1(N,L,NR,NZ)
            WTRTZ=WTRTZ+WTRT1s1(N,L,NR,NZ)
            NINRs1(NR,NZ)=MIN(NINRs1(NR,NZ),NJs1)
            IF(L.EQ.NINRs1(NR,NZ))NRX(N,NR)=1
          ENDIF
        ENDIF
        RTLGZ=RTLGZ+RTLG1s1(N,L,NR,NZ)
        WTRTZ=WTRTZ+WTRT1s1(N,L,NR,NZ)
!     ENDIF
      ENDIF
      NIXs1(NZ)=MAX(NIXs1(NZ),NINRs1(NR,NZ))
5050  CONTINUE
  end associate
  end subroutine GrowRootAxes
!------------------------------------------------------------------------------------------

  subroutine RootNoduleBiomchemistry(I,J,NZ,TFN6,WFNGR)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: TFN6(JZ1)
  real(r8), intent(in) :: WFNGR(2,JZ1)
  integer :: L,M
  real(r8) :: ZADDN
  real(r8) :: ZPOOLD
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLT
  real(r8) :: CPOOLD
  real(r8) :: CCPOLN,CZPOLN
  real(r8) :: CPPOLN
  real(r8) :: CGNDL,CPOOLNX
  real(r8) :: CCNDLR
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: GRNDG
  real(r8) :: PPOOLD
  real(r8) :: PADDN
  real(r8) :: RCO2T
  real(r8) :: RCO2TM
  real(r8) :: RCNDL,RMNDL,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,RGN2F
  real(r8) :: RXNDLC,RXNDLN,RXNDLP
  real(r8) :: RDNDLC,RDNDLN,RDNDLP
  real(r8) :: RCNDLC,RCNDLN,RCNDLP
  real(r8) :: RGNDG
  real(r8) :: RXNSNC,RXNSNN,RXNSNP
  real(r8) :: RDNSNC,RDNSNN,RDNSNP
  real(r8) :: RCNSNC,RCNSNN,RCNSNP
  real(r8) :: RCNDLM,RXNDLM,RGNDLM
  real(r8) :: RSNDLM
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTRTD1,WTNDL1,WTRTDT
  real(r8) :: XFRC,XFRN,XFRP
  real(r8) :: RCCC,RCCN,RCCP
!     begin_execution
  associate(                               &
    NIXs1      =>   plt_morph%NIXs1        &
  )
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
  IF(INTYPs1(NZ).GE.1.AND.INTYPs1(NZ).LE.3)THEN
    DO 5400 L=NUs1,NIXs1(NZ)
      IF(WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
!
!     INITIAL INFECTION
!
        IF(WTNDLs1(L,NZ).LE.0.0)THEN
          WTNDLs1(L,NZ)=WTNDLs1(L,NZ)+WTNDI*AREA3s1(NUs1)
          WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)+WTNDI*AREA3s1(NUs1)*CNNDs1(NZ)
          WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)+WTNDI*AREA3s1(NUs1)*CPNDs1(NZ)
        ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
        IF(WTNDLs1(L,NZ).GT.ZEROPs1(NZ))THEN
          CCPOLN=AMAX1(0.0,CPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
          CZPOLN=AMAX1(0.0,ZPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
          CPPOLN=AMAX1(0.0,PPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
        ELSE
          CCPOLN=1.0_r8
          CZPOLN=1.0_r8
          CPPOLN=1.0_r8
        ENDIF
        IF(CCPOLN.GT.ZEROs1)THEN
          CCC=AMAX1(0.0,AMIN1(1.0,safe_adb(CZPOLN,CZPOLN+CCPOLN*CNKI) &
            ,safe_adb(CPPOLN,CPPOLN+CCPOLN*CPKI)))
!          if(curday==73)write(*,*)CCPOLN,CCPOLN,CZPOLN,CNKI
          CNC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CZPOLN/CNKI)))
          CPC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CPPOLN/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        IF(WTNDLs1(L,NZ).GT.ZEROPs1(NZ))THEN
          FCNPF=AMIN1(1.0 &
            ,SQRT(WTNDLNs1(L,NZ)/(WTNDLs1(L,NZ)*CNNDs1(NZ))) &
            ,SQRT(WTNDLPs1(L,NZ)/(WTNDLs1(L,NZ)*CPNDs1(NZ))))
        ELSE
          FCNPF=1.0_r8
        ENDIF
        SPNDLI=CCPOLN/(CCPOLN+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RCNDLM=respiration from non-structural C unltd by O2
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDL=bacterial C mass
!     TFN4=temperature function for root growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNGR=growth function of root water potential
!
        RCNDLM=AMAX1(0.0,AMIN1(CPOOLNs1(L,NZ) &
          ,VMXO*WTNDLs1(L,NZ))*FCNPF*TFN4s1(L,NZ)*WFNGR(1,L))
        CPOOLNX=CPOOLNs1(L,NZ)
!
!     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCNDL=respiration from non-structural C ltd by O2
!     WFR=constraint by O2 consumption on all root processes
!
        RCNDL=RCNDLM*WFRs1(1,L,NZ)
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN6=temperature function for root maintenance respiration
!     WTNDLN=bacterial N mass
!
        RMNDL=AMAX1(0.0,RMPLT*TFN6(L)*WTNDLNs1(L,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDLM,RXNDL=difference between non-structural C respn and mntc respn unltd,ltd by O2
!     RGNDLM,RGNDL=growth respiration unlimited by N,P and unltd,ltd by O2
!     RSNDLM,RSNDL=excess maintenance respiration unltd,ltd by O2
!
        RXNDLM=RCNDLM-RMNDL
        RXNDL=RCNDL-RMNDL
        RGNDLM=AMAX1(0.0,RXNDLM)
        RGNDL=AMAX1(0.0,RXNDL)
        RSNDLM=AMAX1(0.0,-RXNDLM)
        RSNDL=AMAX1(0.0,-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDL,WTNDLN=bacterial C,N mass
!     CNND=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     RGN2F=respiration for N2 fixation
!     RUPNF,UPNF=layer,total root N2 fixation
!
        RGN2P=AMAX1(0.0,WTNDLs1(L,NZ)*CNNDs1(NZ)-WTNDLNs1(L,NZ))/EN2F
        IF(RGNDL.GT.ZEROPs1(NZ))THEN
          RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
        ELSE
          RGN2F=0._r8
        ENDIF
        RUPNFs1(L,NZ)=RGN2F*EN2F
        UPNFs1(NZ)=UPNFs1(NZ)+RUPNFs1(L,NZ)
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     WTRTD=root C mass
!     CCNDLR=bacteria:root ratio
!     RDNDLX=effect of CCNDLR on bacteria decomposition rate
!     CCNKR=Km for bacterial vs root mass in decomposition
!     SPNDX=specific bacterial decomposition rate at current CCNDLR
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!
        RCCC=RCCZN+CCC*RCCYN
        RCCN=CNC*RCCXN
        RCCP=CPC*RCCQN
        SPNDX=SPNDL*SQRT(TFN4s1(L,NZ)*WFNGR(1,L))
        RXNDLC=SPNDX*WTNDLs1(L,NZ)
        RXNDLN=SPNDX*WTNDLNs1(L,NZ)
        RXNDLP=SPNDX*WTNDLPs1(L,NZ)
        RDNDLC=RXNDLC*(1.0_r8-RCCC)
        RDNDLN=RXNDLN*(1.0_r8-RCCC)*(1.0_r8-RCCN)
        RDNDLP=RXNDLP*(1.0_r8-RCCC)*(1.0_r8-RCCP)
        RCNDLC=RXNDLC-RDNDLC
        RCNDLN=RXNDLN-RDNDLN
        RCNDLP=RXNDLP-RDNDLP
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLC=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     RGN2F=respiration for N2 fixation
!     GRNDG=bacterial growth
!     DMND=bacterial growth yield
!     RGNDG=bacterial respiration for growth and N2 fixation
!     ZADDN,PADDN=nonstructural N,P used in growth
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
        CGNDL=AMIN1(CPOOLNs1(L,NZ)-AMIN1(RMNDL,RCNDL) &
          -RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0_r8-DMNDs1(NZ)))
        GRNDG=CGNDL*DMNDs1(NZ)
        RGNDG=RGN2F+CGNDL*(1.0_r8-DMNDs1(NZ))
        ZADDN=AMAX1(0.0,AMIN1(ZPOOLNs1(L,NZ),GRNDG*CNNDs1(NZ)))*CZPOLN/(CZPOLN+CZKM)
        PADDN=AMAX1(0.0,AMIN1(PPOOLNs1(L,NZ),GRNDG*CPNDs1(NZ)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!
        IF(RSNDL.GT.0.0.AND.WTNDLs1(L,NZ).GT.ZEROPs1(NZ).AND.RCCC.GT.ZEROs1)THEN
          RXNSNC=RSNDL/RCCC
          RXNSNN=RXNSNC*WTNDLNs1(L,NZ)/WTNDLs1(L,NZ)
          RXNSNP=RXNSNC*WTNDLPs1(L,NZ)/WTNDLs1(L,NZ)
          RDNSNC=RXNSNC*(1.0_r8-RCCC)
          RDNSNN=RXNSNN*(1.0_r8-RCCC)*(1.0_r8-RCCN)
          RDNSNP=RXNSNP*(1.0_r8-RCCC)*(1.0_r8-RCCP)
          RCNSNC=RXNSNC-RDNSNC
          RCNSNN=RXNSNN-RDNSNN
          RCNSNP=RXNSNP-RDNSNP
        ELSE
          RXNSNC=0._r8
          RXNSNN=0._r8
          RXNSNP=0._r8
          RDNSNC=0._r8
          RDNSNN=0._r8
          RDNSNP=0._r8
          RCNSNC=0._r8
          RCNSNN=0._r8
          RCNSNP=0._r8
        ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNC=bacterial C senescence to recycling
!     RCO2A=total root respiration
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
        RCO2TM=AMIN1(RMNDL,RCNDLM)+RGNDLM+RCNSNC
        RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
        RCO2Ms1(1,L,NZ)=RCO2Ms1(1,L,NZ)+RCO2TM
        RCO2Ns1(1,L,NZ)=RCO2Ns1(1,L,NZ)+RCO2T
        RCO2As1(1,L,NZ)=RCO2As1(1,L,NZ)-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!
        DO 6370 M=1,4
          CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ)*(RDNDLC+RDNSNC)
          ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ)*(RDNDLN+RDNSNN)
          PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ)*(RDNDLP+RDNSNP)
6370    CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNF=root N2 fixation
!
        CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLC
        ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)-ZADDN+RCNDLN+RCNSNN+RUPNFs1(L,NZ)
        PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)-PADDN+RCNDLP+RCNSNP
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
        WTNDLs1(L,NZ)=WTNDLs1(L,NZ)+GRNDG-RXNDLC-RXNSNC
        WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)+ZADDN-RXNDLN-RXNSNN
        WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)+PADDN-RXNDLP-RXNSNP
!     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,2122)'NODGR',I,J,NZ,L,RCNDLM,RCNDL,RMNDL,RGNDL,RGN2P
!    2,RGN2F,CGNDL,RSNDL,GRNDG,ZADDN,PADDN,RCCC,RCCN,RCCP
!    8,RDNDLC,RDNDLN,RDNDLP,RCNDLC,RDNDLX,WFRs1(1,L,NZ)
!    3,WTNDLs1(L,NZ),WTNDLNs1(L,NZ),WTNDLPs1(L,NZ)
!    2,CPOOLNs1(L,NZ),ZPOOLNs1(L,NZ),PPOOLNs1(L,NZ)
!    5,FCNPF,TFN4s1(L,NZ),WFNGR(1,L),PSIRTs1(1,L,NZ)
!    5,CCPOLN,CZPOLN,CPPOLN,CPOOLNX
!    6,VMXO*WTNDLs1(L,NZ)*TFN4s1(L,NZ)*FCNPF*WFNGR(1,L)
!2122  FORMAT(A8,4I4,60E14.6)
!     ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOLR,ZPOOLR,PPOOLR=root non-structural C,N,P mass
!     WTRTD=root C mass
!     WTNDL=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGR=parameter to calculate nonstructural C,N,P exchange
!     CCNDLR=bacteria:root ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(CPOOLRs1(1,L,NZ).GT.ZEROPs1(NZ) &
          .AND.WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
          CCNDLR=WTNDLs1(L,NZ)/WTRTDs1(1,L,NZ)
          WTRTD1=WTRTDs1(1,L,NZ)
          WTNDL1=AMIN1(WTRTDs1(1,L,NZ) &
          ,AMAX1(WTNDI*AREA3s1(NUs1),WTNDLs1(L,NZ)))
          WTRTDT=WTRTD1+WTNDL1
          IF(WTRTDT.GT.ZEROPs1(NZ))THEN
            FXRNX=FXRN(INTYPs1(NZ))/(1.0+CCNDLR/CCNGR)
!    2/(1.0+CCNDLR/(CCNGR*FXRN(INTYPs1(NZ))))
            CPOOLD=(CPOOLRs1(1,L,NZ)*WTNDL1-CPOOLNs1(L,NZ)*WTRTD1)/WTRTDT
            XFRC=FXRNX*CPOOLD
            CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)-XFRC
            CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)+XFRC
            CPOOLT=CPOOLRs1(1,L,NZ)+CPOOLNs1(L,NZ)
            IF(CPOOLT.GT.ZEROPs1(NZ))THEN
              ZPOOLD=(ZPOOLRs1(1,L,NZ)*CPOOLNs1(L,NZ) &
                -ZPOOLNs1(L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              XFRN=FXRNX*ZPOOLD
              PPOOLD=(PPOOLRs1(1,L,NZ)*CPOOLNs1(L,NZ) &
                -PPOOLNs1(L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              XFRP=FXRNX*PPOOLD
              ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)-XFRN
              PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)-XFRP
              ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)+XFRN
              PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)+XFRP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
5400  CONTINUE
  ENDIF
  end associate
  end subroutine RootNoduleBiomchemistry
!------------------------------------------------------------------------------------------

  subroutine PrimRootRemobilization(N,L,NZ,NR,FSNC1,RCO2X,RCO2XM)

  implicit none
  integer, intent(in) :: N,L,NZ,NR
  real(r8), intent(in) :: RCO2X,RCO2XM
  real(r8), intent(out) :: FSNC1
  integer :: M
  real(r8) :: CCC,CNC,CPC
  real(r8) :: RCCR,RCZR,RCPR
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SNCR,SNCRM
! begin_execution
  associate(                              &
    NB1s1     =>  plt_morph%NB1s1         &
  )
!
!     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
!     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LITTERFALL
!
!     IDAYs1(1,=emergence date
!     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZR,RCCYR=min,max fractions for root C recycling
!     RCCXR,RCCQR=max fractions for root N,P recycling
!
  IF(IDAYs1(1,NB1s1(NZ),NZ).NE.0 &
    .AND.CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
    CCC=AMAX1(0.0,AMIN1(1.0 &
      ,safe_adb(CZPOLRs1(N,L,NZ),CZPOLRs1(N,L,NZ) &
      +CCPOLRs1(N,L,NZ)*CNKI) &
      ,safe_adb(CPPOLRs1(N,L,NZ),CPPOLRs1(N,L,NZ) &
      +CCPOLRs1(N,L,NZ)*CPKI)))
    CNC=AMAX1(0.0,AMIN1(1.0 &
      ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ) &
      +CZPOLRs1(N,L,NZ)/CNKI)))
    CPC=AMAX1(0.0,AMIN1(1.0 &
      ,safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ) &
      +CPPOLRs1(N,L,NZ)/CPKI)))
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
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     FSNC1=fraction of primary root C to be remobilized
!
  IF(-RCO2XM.GT.0.0)THEN
    IF(-RCO2XM.LT.RTWT1s1(N,NR,NZ)*RCCC)THEN
      SNCRM=-RCO2XM
    ELSE
      SNCRM=AMAX1(0.0,RTWT1s1(N,NR,NZ)*RCCC)
    ENDIF
  ELSE
    SNCRM=0._r8
  ENDIF
  IF(-RCO2X.GT.0.0)THEN
    IF(-RCO2X.LT.RTWT1s1(N,NR,NZ)*RCCC)THEN
      SNCR=-RCO2X
    ELSE
      SNCR=AMAX1(0.0,RTWT1s1(N,NR,NZ)*RCCC)*WFRs1(N,L,NZ)
    ENDIF
  ELSE
    SNCR=0._r8
  ENDIF
  IF(SNCR.GT.0.0.AND.RTWT1s1(N,NR,NZ).GT.ZEROPs1(NZ))THEN
    RCCR=RCCC*RTWT1s1(N,NR,NZ)
    RCZR=RTWT1Ns1(N,NR,NZ)*(RCCN+(1.0_r8-RCCN)*RCCR/RTWT1s1(N,NR,NZ))
    RCPR=RTWT1Ps1(N,NR,NZ)*(RCCP+(1.0_r8-RCCP)*RCCR/RTWT1s1(N,NR,NZ))
    IF(RCCR.GT.ZEROPs1(NZ))THEN
      FSNC1=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
    ELSE
      FSNC1=1.0_r8
    ENDIF
  ELSE
    RCCR=0._r8
    RCZR=0._r8
    RCPR=0._r8
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
!     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!
  DO 6355 M=1,4
    CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
      *FSNC1*(RTWT1s1(N,NR,NZ)-RCCR)*FWODRs1(0)
    ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
      *FSNC1*(RTWT1Ns1(N,NR,NZ)-RCZR)*FWODRNs1(0)
    PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
      *FSNC1*(RTWT1Ps1(N,NR,NZ)-RCPR)*FWODRPs1(0)
    CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
      *FSNC1*(RTWT1s1(N,NR,NZ)-RCCR)*FWODRs1(1)
    ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
      *FSNC1*(RTWT1Ns1(N,NR,NZ)-RCZR)*FWODRNs1(1)
    PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
      *FSNC1*(RTWT1Ps1(N,NR,NZ)-RCPR)*FWODRPs1(1)
6355  CONTINUE
  end associate
  end subroutine PrimRootRemobilization

!------------------------------------------------------------------------------------------
  subroutine PrimRootExtension(L,L1,N,NR,NZ,WFNR,FRTN,GRTWTG,GRTWTL,GRTWTN,GRTWTP,GRTLGL,RTLGZ,WTRTZ)
  implicit none
  integer, intent(in) :: L,L1,N,NR,NZ
  real(r8), intent(in):: WFNR,FRTN,GRTWTG,GRTWTL,GRTWTN,GRTWTP
  real(r8), intent(inout) :: RTLGZ,WTRTZ
  real(r8), intent(out):: GRTLGL
  real(r8) :: FGROL,FGROZ
  REAL(R8) :: XFRC,XFRN,XFRP
! begin_execution
  associate(                                 &
    NINRs1    =>  plt_morph%NINRs1         , &
    SDPTHs1     =>  plt_morph%SDPTHs1        &
  )
!     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
!
!     GRTLGL=primary root length extension
!     GRTWTG=primary root C growth ltd by O2
!     RTLG1X=specific primary root length from startq.f
!     PP=PFT population
!     WFNR=water function for root extension
!     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
!     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
!     RTDP1=primary root depth from soil surface
!     SDPTH=seeding depth
!     FSNC1=fraction of primary root C to be remobilized
!     RTLG1=primary root length
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     DLYR=soil layer thickness
!
  IF(GRTWTL.LT.0.0.AND.RTWT1s1(N,NR,NZ) &
    .GT.ZEROPs1(NZ))THEN
    GRTLGL=GRTWTG*RTLG1Xs1(N,NZ)/PPs1(NZ)*WFNR*FWODRs1(1) &
      +GRTWTL*(RTDP1s1(N,NR,NZ)-SDPTHs1(NZ)) &
      /RTWT1s1(N,NR,NZ)
  ELSE
    GRTLGL=GRTWTG*RTLG1Xs1(N,NZ)/PPs1(NZ)*WFNR*FWODRs1(1)
  ENDIF
  IF(L.LT.NJs1)THEN
    GRTLGL=AMIN1(DLYR3s1(L1),GRTLGL)
  ENDIF
!
!     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
!     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
!     BOUNDARY OF CURRENT LAYER
!
!     GRTLGL=primary root length extension
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!
  IF(GRTLGL.GT.ZEROPs1(NZ).AND.L.LT.NJs1)THEN
    FGROL=AMAX1(0.0,AMIN1(1.0,(CDPTHZs1(L) &
      -RTDP1s1(N,NR,NZ))/GRTLGL))
    IF(FGROL.LT.1.0)FGROL=0._r8
    FGROZ=AMAX1(0.0,1.0_r8-FGROL)
  ELSE
    FGROL=1.0_r8
    FGROZ=0._r8
  ENDIF
!
!     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
!     AND AXIS NUMBER
!
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
!     GRTLGL=primary root length extension
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
!     WSRTL=total root protein C mass
!     CNWS,CPWS=protein:N,protein:P ratios from startq.f
!     RTLG1=primary root length
!
  RTWT1s1(N,NR,NZ)=RTWT1s1(N,NR,NZ)+GRTWTL
  RTWT1Ns1(N,NR,NZ)=RTWT1Ns1(N,NR,NZ)+GRTWTN
  RTWT1Ps1(N,NR,NZ)=RTWT1Ps1(N,NR,NZ)+GRTWTP
  RTDP1s1(N,NR,NZ)=RTDP1s1(N,NR,NZ)+GRTLGL
  WTRT1s1(N,L,NR,NZ)=WTRT1s1(N,L,NR,NZ)+GRTWTL*FGROL
  WTRT1Ns1(N,L,NR,NZ)=WTRT1Ns1(N,L,NR,NZ)+GRTWTN*FGROL
  WTRT1Ps1(N,L,NR,NZ)=WTRT1Ps1(N,L,NR,NZ)+GRTWTP*FGROL
  WSRTLs1(N,L,NZ)=WSRTLs1(N,L,NZ) &
    +AMIN1(CNWSs1(NZ)*WTRT1Ns1(N,L,NR,NZ) &
    ,CPWSs1(NZ)*WTRT1Ps1(N,L,NR,NZ))
  RTLG1s1(N,L,NR,NZ)=RTLG1s1(N,L,NR,NZ)+GRTLGL*FGROL
!
!     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
!     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
!     OF CURRENT SOIL LAYER
!
!     FGROZ=fraction of GRTLGL in next lower soil layer
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
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
    WTRT1s1(N,L1,NR,NZ)=WTRT1s1(N,L1,NR,NZ)+GRTWTL*FGROZ
    WTRT1Ns1(N,L1,NR,NZ)=WTRT1Ns1(N,L1,NR,NZ)+GRTWTN*FGROZ
    WTRT1Ps1(N,L1,NR,NZ)=WTRT1Ps1(N,L1,NR,NZ)+GRTWTP*FGROZ
    WSRTLs1(N,L1,NZ)=WSRTLs1(N,L1,NZ)+AMIN1(CNWSs1(NZ)*WTRT1Ns1(N,L1,NR,NZ) &
      ,CPWSs1(NZ)*WTRT1Ps1(N,L1,NR,NZ))
    WTRTDs1(N,L1,NZ)=WTRTDs1(N,L1,NZ)+WTRT1s1(N,L1,NR,NZ)
    RTLG1s1(N,L1,NR,NZ)=RTLG1s1(N,L1,NR,NZ)+GRTLGL*FGROZ
    RRAD1s1(N,L1,NZ)=RRAD1s1(N,L,NZ)
    RTLGZ=RTLGZ+RTLG1s1(N,L1,NR,NZ)
    WTRTZ=WTRTZ+WTRT1s1(N,L1,NR,NZ)
    XFRC=FRTN*CPOOLRs1(N,L,NZ)
    XFRN=FRTN*ZPOOLRs1(N,L,NZ)
    XFRP=FRTN*PPOOLRs1(N,L,NZ)
    CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-XFRC
    ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-XFRN
    PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-XFRP
    CPOOLRs1(N,L1,NZ)=CPOOLRs1(N,L1,NZ)+XFRC
    ZPOOLRs1(N,L1,NZ)=ZPOOLRs1(N,L1,NZ)+XFRN
    PPOOLRs1(N,L1,NZ)=PPOOLRs1(N,L1,NZ)+XFRP
    PSIRTs1(N,L1,NZ)=PSIRTs1(N,L,NZ)
    PSIROs1(N,L1,NZ)=PSIROs1(N,L,NZ)
    PSIRGs1(N,L1,NZ)=PSIRGs1(N,L,NZ)
    NINRs1(NR,NZ)=MAX(NGs1(NZ),L+1)
!     IF(NZ.EQ.2)THEN
!     WRITE(*,9877)'INFIL',I,J,NZ,NR,L,N,NINRs1(NR,NZ)
!    2,FRTN,WTRTDs1(N,L1,NZ),CPOOLRs1(N,L1,NZ)
!    2,FGROZ,RTDP1s1(N,NR,NZ),GRTLGL,CDPTHZs1(L)
!     ENDIF
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
  integer :: LL,NN
  real(r8) :: XFRD,XFRW,FRTN
  real(r8) :: XFRC,XFRN,XFRP

! begin_execution
  associate(                             &
    SDPTHs1    =>  plt_morph%SDPTHs1   , &
    NINRs1     =>  plt_morph%NINRs1      &
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

  DO 5115 LL=L,NGs1(NZ)+1,-1
    IF(VOLXs1(LL-1).GT.ZEROS2s1 &
      .AND.(RTDP1X.LT.CDPTHZs1(LL-1) &
      .OR.RTDP1X.LT.SDPTHs1(NZ)))THEN
      IF(RLNT(N,LL).GT.ZEROPs1(NZ))THEN
        FRTN=(RTSK1(N,LL,NR)+RTSK2(N,LL,NR))/RLNT(N,LL)
      ELSE
        FRTN=1.0_r8
      ENDIF
      DO 5110 NN=1,MYs1(NZ)
        WTRT1s1(NN,LL-1,NR,NZ)=WTRT1s1(NN,LL-1,NR,NZ)+WTRT1s1(NN,LL,NR,NZ)
        WTRT1Ns1(NN,LL-1,NR,NZ)=WTRT1Ns1(NN,LL-1,NR,NZ)+WTRT1Ns1(NN,LL,NR,NZ)
        WTRT1Ps1(NN,LL-1,NR,NZ)=WTRT1Ps1(NN,LL-1,NR,NZ)+WTRT1Ps1(NN,LL,NR,NZ)
        WTRT2s1(NN,LL-1,NR,NZ)=WTRT2s1(NN,LL-1,NR,NZ)+WTRT2s1(NN,LL,NR,NZ)
        WTRT2Ns1(NN,LL-1,NR,NZ)=WTRT2Ns1(NN,LL-1,NR,NZ)+WTRT2Ns1(NN,LL,NR,NZ)
        WTRT2Ps1(NN,LL-1,NR,NZ)=WTRT2Ps1(NN,LL-1,NR,NZ)+WTRT2Ps1(NN,LL,NR,NZ)
        RTLG1s1(NN,LL-1,NR,NZ)=RTLG1s1(NN,LL-1,NR,NZ)+RTLG1s1(NN,LL,NR,NZ)
        WTRT1s1(NN,LL,NR,NZ)=0._r8
        WTRT1Ns1(NN,LL,NR,NZ)=0._r8
        WTRT1Ps1(NN,LL,NR,NZ)=0._r8
        WTRT2s1(NN,LL,NR,NZ)=0._r8
        WTRT2Ns1(NN,LL,NR,NZ)=0._r8
        WTRT2Ps1(NN,LL,NR,NZ)=0._r8
        RTLG1s1(NN,LL,NR,NZ)=0._r8
        XFRC=FRTN*CPOOLRs1(NN,LL,NZ)
        XFRN=FRTN*ZPOOLRs1(NN,LL,NZ)
        XFRP=FRTN*PPOOLRs1(NN,LL,NZ)
        XFRW=FRTN*WSRTLs1(NN,L,NZ)
        XFRD=FRTN*WTRTDs1(NN,LL,NZ)
        CPOOLRs1(NN,LL,NZ)=CPOOLRs1(NN,LL,NZ)-XFRC
        ZPOOLRs1(NN,LL,NZ)=ZPOOLRs1(NN,LL,NZ)-XFRN
        PPOOLRs1(NN,LL,NZ)=PPOOLRs1(NN,LL,NZ)-XFRP
        WSRTLs1(NN,LL,NZ)=WSRTLs1(NN,LL,NZ)-XFRW
        WTRTDs1(NN,LL,NZ)=WTRTDs1(NN,LL,NZ)-XFRD
        CPOOLRs1(NN,LL-1,NZ)=CPOOLRs1(NN,LL-1,NZ)+XFRC
        ZPOOLRs1(NN,LL-1,NZ)=ZPOOLRs1(NN,LL-1,NZ)+XFRN
        PPOOLRs1(NN,LL-1,NZ)=PPOOLRs1(NN,LL-1,NZ)+XFRP
        WSRTLs1(NN,LL-1,NZ)=WSRTLs1(NN,LL-1,NZ)+XFRW
        WTRTDs1(NN,LL-1,NZ)=WTRTDs1(NN,LL-1,NZ)+XFRD
!
!     WITHDRAW GASES IN PRIMARY ROOTS
!
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
!     FRTN=fraction of primary root sink strength in axis
!
        RCO2Zs1(NZ)=RCO2Zs1(NZ)-FRTN*(CO2As1(NN,LL,NZ)+CO2Ps1(NN,LL,NZ))
        ROXYZs1(NZ)=ROXYZs1(NZ)-FRTN*(OXYAs1(NN,LL,NZ)+OXYPs1(NN,LL,NZ))
        RCH4Zs1(NZ)=RCH4Zs1(NZ)-FRTN*(CH4As1(NN,LL,NZ)+CH4Ps1(NN,LL,NZ))
        RN2OZs1(NZ)=RN2OZs1(NZ)-FRTN*(Z2OAs1(NN,LL,NZ)+Z2OPs1(NN,LL,NZ))
        RNH3Zs1(NZ)=RNH3Zs1(NZ)-FRTN*(ZH3As1(NN,LL,NZ)+ZH3Ps1(NN,LL,NZ))
        RH2GZs1(NZ)=RH2GZs1(NZ)-FRTN*(H2GAs1(NN,LL,NZ)+H2GPs1(NN,LL,NZ))
        CO2As1(NN,LL,NZ)=(1.0_r8-FRTN)*CO2As1(NN,LL,NZ)
        OXYAs1(NN,LL,NZ)=(1.0_r8-FRTN)*OXYAs1(NN,LL,NZ)
        CH4As1(NN,LL,NZ)=(1.0_r8-FRTN)*CH4As1(NN,LL,NZ)
        Z2OAs1(NN,LL,NZ)=(1.0_r8-FRTN)*Z2OAs1(NN,LL,NZ)
        ZH3As1(NN,LL,NZ)=(1.0_r8-FRTN)*ZH3As1(NN,LL,NZ)
        H2GAs1(NN,LL,NZ)=(1.0_r8-FRTN)*H2GAs1(NN,LL,NZ)
        CO2Ps1(NN,LL,NZ)=(1.0_r8-FRTN)*CO2Ps1(NN,LL,NZ)
        OXYPs1(NN,LL,NZ)=(1.0_r8-FRTN)*OXYPs1(NN,LL,NZ)
        CH4Ps1(NN,LL,NZ)=(1.0_r8-FRTN)*CH4Ps1(NN,LL,NZ)
        Z2OPs1(NN,LL,NZ)=(1.0_r8-FRTN)*Z2OPs1(NN,LL,NZ)
        ZH3Ps1(NN,LL,NZ)=(1.0_r8-FRTN)*ZH3Ps1(NN,LL,NZ)
        H2GPs1(NN,LL,NZ)=(1.0_r8-FRTN)*H2GPs1(NN,LL,NZ)

5110  CONTINUE
!
!     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
!
!     RTN2,RTNL=number of secondary root axes
!     RTN1=number of primary root axes
!     RTLG1=primary root length
!     CDPTHZ=depth from soil surface to layer bottom
!     SDPTH=seeding depth
!
      RTNLs1(N,LL,NZ)=RTNLs1(N,LL,NZ)-RTN2s1(N,LL,NR,NZ)
      RTNLs1(N,LL-1,NZ)=RTNLs1(N,LL-1,NZ)+RTN2s1(N,LL,NR,NZ)
      RTN2s1(N,LL,NR,NZ)=0._r8
      RTN1s1(N,LL,NZ)=RTN1s1(N,LL,NZ)-XRTN1
      IF(LL-1.GT.NGs1(NZ))THEN
        RTLG1s1(N,LL-1,NR,NZ)=DLYR3s1(LL-1)-(CDPTHZs1(LL-1)-RTDP1s1(N,NR,NZ))
      ELSE
        RTLG1s1(N,LL-1,NR,NZ)=DLYR3s1(LL-1)-(CDPTHZs1(LL-1)-RTDP1s1(N,NR,NZ)) &
          -(SDPTHs1(NZ)-CDPTHZs1(LL-2))
      ENDIF
!
!     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     FRTN=fraction of primary root sink strength in axis
!     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
!
      IF(INTYPs1(NZ).GE.1.AND.INTYPs1(NZ).LE.3)THEN
        XFRC=FRTN*WTNDLs1(LL,NZ)
        XFRN=FRTN*WTNDLNs1(LL,NZ)
        XFRP=FRTN*WTNDLPs1(LL,NZ)
        WTNDLs1(LL,NZ)=WTNDLs1(LL,NZ)-XFRC
        WTNDLNs1(LL,NZ)=WTNDLNs1(LL,NZ)-XFRN
        WTNDLPs1(LL,NZ)=WTNDLPs1(LL,NZ)-XFRP
        WTNDLs1(LL-1,NZ)=WTNDLs1(LL-1,NZ)+XFRC
        WTNDLNs1(LL-1,NZ)=WTNDLNs1(LL-1,NZ)+XFRN
        WTNDLPs1(LL-1,NZ)=WTNDLPs1(LL-1,NZ)+XFRP
        XFRC=FRTN*CPOOLNs1(LL,NZ)
        XFRN=FRTN*ZPOOLNs1(LL,NZ)
        XFRP=FRTN*PPOOLNs1(LL,NZ)
        CPOOLNs1(LL,NZ)=CPOOLNs1(LL,NZ)-XFRC
        ZPOOLNs1(LL,NZ)=ZPOOLNs1(LL,NZ)-XFRN
        PPOOLNs1(LL,NZ)=PPOOLNs1(LL,NZ)-XFRP
        CPOOLNs1(LL-1,NZ)=CPOOLNs1(LL-1,NZ)+XFRC
        ZPOOLNs1(LL-1,NZ)=ZPOOLNs1(LL-1,NZ)+XFRN
        PPOOLNs1(LL-1,NZ)=PPOOLNs1(LL-1,NZ)+XFRP
!     WRITE(*,9868)'WITHDRN',I,J,NZ,NR,LL,NN,NINRs1(NR,NZ)
!    2,WTNDLs1(LL,NZ),CPOOLNs1(LL,NZ),RTDP1s1(N,NR,NZ)
      ENDIF
      NINRs1(NR,NZ)=MAX(NGs1(NZ),LL-1)
    ELSE
      EXIT
    ENDIF
5115  CONTINUE
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
  integer :: L,NB,N,NR
  real(r8) :: ZPOOLS,ZPOOLT
  real(r8) :: ZPOOLB
  real(r8) :: ZPOOLD
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
  real(r8) :: XFRC,XFRN,XFRP
!     begin_execution
  associate(                                &
    NIXs1        =>   plt_morph%NIXs1     , &
    NINRs1       =>   plt_morph%NINRs1    , &
    NIs1         =>   plt_morph%NIs1      , &
    NRTs1        =>   plt_morph%NRTs1     , &
    NBRs1        =>   plt_morph%NBRs1       &
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
  IF(NBRs1(NZ).GT.1)THEN
    WTPLTT=0._r8
    CPOOLT=0._r8
    ZPOOLT=0._r8
    PPOOLT=0._r8
    DO 300 NB=1,NBRs1(NZ)
      IF(IDTHBs1(NB,NZ).EQ.0)THEN
        IF(ATRPs1(NB,NZ).GT.ATRPX(ISTYPs1(NZ)))THEN
          WTLSBZ(NB)=AMAX1(0.0,WTLSBs1(NB,NZ))
          CPOOLZ(NB)=AMAX1(0.0,CPOOLs1(NB,NZ))
          ZPOOLZ(NB)=AMAX1(0.0,ZPOOLs1(NB,NZ))
          PPOOLZ(NB)=AMAX1(0.0,PPOOLs1(NB,NZ))
          WTPLTT=WTPLTT+WTLSBZ(NB)
          CPOOLT=CPOOLT+CPOOLZ(NB)
          ZPOOLT=ZPOOLT+ZPOOLZ(NB)
          PPOOLT=PPOOLT+PPOOLZ(NB)
        ENDIF
      ENDIF
300 CONTINUE
    DO 305 NB=1,NBRs1(NZ)
      IF(IDTHBs1(NB,NZ).EQ.0)THEN
        IF(ATRPs1(NB,NZ).GT.ATRPX(ISTYPs1(NZ)))THEN
          IF(WTPLTT.GT.ZEROPs1(NZ) &
            .AND.CPOOLT.GT.ZEROPs1(NZ))THEN
            CPOOLD=CPOOLT*WTLSBZ(NB)-CPOOLZ(NB)*WTPLTT
            ZPOOLD=ZPOOLT*CPOOLZ(NB)-ZPOOLZ(NB)*CPOOLT
            PPOOLD=PPOOLT*CPOOLZ(NB)-PPOOLZ(NB)*CPOOLT
            XFRC=0.01*CPOOLD/WTPLTT
            XFRN=0.01*ZPOOLD/CPOOLT
            XFRP=0.01*PPOOLD/CPOOLT
            CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)+XFRC
            ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)+XFRN
            PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)+XFRP
          ENDIF
        ENDIF
      ENDIF
305 CONTINUE
  ENDIF
!
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH STALK RESERVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     WVSTKB=stalk sapwood mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     IDAYs1(7,=start of grain filling and setting max seed size
!
  IF(NBRs1(NZ).GT.1)THEN
    WTSTKT=0._r8
    WTRSVT=0._r8
    WTRSNT=0._r8
    WTRSPT=0._r8
    DO 330 NB=1,NBRs1(NZ)
      IF(IDTHBs1(NB,NZ).EQ.0)THEN
        IF(IDAYs1(7,NB,NZ).NE.0)THEN
          WTSTKT=WTSTKT+WVSTKBs1(NB,NZ)
          WTRSVT=WTRSVT+WTRSVBs1(NB,NZ)
          WTRSNT=WTRSNT+WTRSBNs1(NB,NZ)
          WTRSPT=WTRSPT+WTRSBPs1(NB,NZ)
        ENDIF
      ENDIF
330 CONTINUE
    IF(WTSTKT.GT.ZEROPs1(NZ) &
      .AND.WTRSVT.GT.ZEROPs1(NZ))THEN
      DO 335 NB=1,NBRs1(NZ)
        IF(IDTHBs1(NB,NZ).EQ.0)THEN
          IF(IDAYs1(7,NB,NZ).NE.0)THEN
            WTRSVD=WTRSVT*WVSTKBs1(NB,NZ)-WTRSVBs1(NB,NZ)*WTSTKT
            XFRC=0.1*WTRSVD/WTSTKT
            WTRSVBs1(NB,NZ)=WTRSVBs1(NB,NZ)+XFRC
            WTRSND=WTRSNT*WTRSVBs1(NB,NZ)-WTRSBNs1(NB,NZ)*WTRSVT
            XFRN=0.1*WTRSND/WTRSVT
            WTRSBNs1(NB,NZ)=WTRSBNs1(NB,NZ)+XFRN
            WTRSPD=WTRSPT*WTRSVBs1(NB,NZ)-WTRSBPs1(NB,NZ)*WTRSVT
            XFRP=0.1*WTRSPD/WTRSVT
            WTRSBPs1(NB,NZ)=WTRSBPs1(NB,NZ)+XFRP
          ENDIF
        ENDIF
335   CONTINUE
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
  IF(MYs1(NZ).EQ.2)THEN
    DO 425 L=NUs1,NIXs1(NZ)
      IF(CPOOLRs1(1,L,NZ).GT.ZEROPs1(NZ) &
        .AND.WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
        WTRTD1=WTRTDs1(1,L,NZ)
        WTRTD2=AMIN1(WTRTDs1(1,L,NZ),AMAX1(FSNK &
          *WTRTDs1(1,L,NZ),WTRTDs1(2,L,NZ)))
        WTPLTT=WTRTD1+WTRTD2
        IF(WTPLTT.GT.ZEROPs1(NZ))THEN
          CPOOLD=(CPOOLRs1(1,L,NZ)*WTRTD2 &
            -CPOOLRs1(2,L,NZ)*WTRTD1)/WTPLTT
          XFRC=FMYC*CPOOLD
          CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)-XFRC
          CPOOLRs1(2,L,NZ)=CPOOLRs1(2,L,NZ)+XFRC
          CPOOLT=CPOOLRs1(1,L,NZ)+CPOOLRs1(2,L,NZ)
          IF(CPOOLT.GT.ZEROPs1(NZ))THEN
            ZPOOLD=(ZPOOLRs1(1,L,NZ)*CPOOLRs1(2,L,NZ) &
              -ZPOOLRs1(2,L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
            XFRN=FMYC*ZPOOLD
            PPOOLD=(PPOOLRs1(1,L,NZ)*CPOOLRs1(2,L,NZ) &
              -PPOOLRs1(2,L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
            XFRP=FMYC*PPOOLD
            ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)-XFRN
            ZPOOLRs1(2,L,NZ)=ZPOOLRs1(2,L,NZ)+XFRN
            PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)-XFRP
            PPOOLRs1(2,L,NZ)=PPOOLRs1(2,L,NZ)+XFRP
    !     IF(L.EQ.NIXs1(NZ))THEN
    !     WRITE(*,9873)'MYCO',I,J,NZ,L,XFRC,XFRN,XFRP
!    2,CPOOLRs1(1,L,NZ),WTRTDs1(1,L,NZ)
!    3,CPOOLRs1(2,L,NZ),WTRTD2
!    3,WTPLTT,ZPOOLRs1(1,L,NZ),ZPOOLRs1(2,L,NZ)
!    4,PPOOLRs1(1,L,NZ),PPOOLRs1(2,L,NZ),CPOOLT
!9873  FORMAT(A8,4I4,20E24.16)
!     ENDIF
          ENDIF
        ENDIF
      ENDIF
425 CONTINUE
  ENDIF
!
!     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
!     IN PERENNIALS
!
  IF(IFLGZ.EQ.1.AND.ISTYPs1(NZ).NE.0)THEN
    DO 5545 N=1,MYs1(NZ)
      DO 5550 L=NUs1,NIs1(NZ)
        IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
          CNL=CCPOLRs1(N,L,NZ)/(CCPOLRs1(N,L,NZ) &
            +CZPOLRs1(N,L,NZ)/CNKI)
          CPL=CCPOLRs1(N,L,NZ)/(CCPOLRs1(N,L,NZ) &
            +CPPOLRs1(N,L,NZ)/CPKI)
        ELSE
          CNL=0._r8
          CPL=0._r8
        ENDIF
        XFRCX=FXFR(IGTYPs1(NZ)) &
          *AMAX1(0.0,CPOOLRs1(N,L,NZ))
        XFRNX=FXFR(IGTYPs1(NZ)) &
          *AMAX1(0.0,ZPOOLRs1(N,L,NZ))*(1.0+CNL)
        XFRPX=FXFR(IGTYPs1(NZ)) &
          *AMAX1(0.0,PPOOLRs1(N,L,NZ))*(1.0+CPL)
        XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
        XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN*0.5)
        XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN*0.5)
        CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-XFRC
        WTRVCs1(NZ)=WTRVCs1(NZ)+XFRC
        ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)-XFRN
        WTRVNs1(NZ)=WTRVNs1(NZ)+XFRN
        PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)-XFRP
        WTRVPs1(NZ)=WTRVPs1(NZ)+XFRP
5550  CONTINUE
5545  CONTINUE
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
  DO 5445 N=1,MYs1(NZ)
    DO 5450 L=NUs1,NIs1(NZ)
      WTRTLs1(N,L,NZ)=0._r8
      WTRTDs1(N,L,NZ)=0._r8
      DO 5460 NR=1,NRTs1(NZ)
        WTRTLs1(N,L,NZ)=WTRTLs1(N,L,NZ)+WTRT2s1(N,L,NR,NZ)
        WTRTDs1(N,L,NZ)=WTRTDs1(N,L,NZ)+WTRT2s1(N,L,NR,NZ) &
          +WTRT1s1(N,L,NR,NZ)
5460  CONTINUE
      TCO2Ts1(NZ)=TCO2Ts1(NZ)+RCO2As1(N,L,NZ)
      RECOs1=RECOs1+RCO2As1(N,L,NZ)
      TRAUs1=TRAUs1+RCO2As1(N,L,NZ)
5450  CONTINUE
    DO 5470 NR=1,NRTs1(NZ)
      WTRTLs1(N,NINRs1(NR,NZ),NZ)=WTRTLs1(N,NINRs1(NR,NZ),NZ) &
        +RTWT1s1(N,NR,NZ)
5470  CONTINUE
5445  CONTINUE
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
!     IF(ISTYPs1(NZ).EQ.1)THEN
  IF(WTLSs1(NZ).GT.ZEROPs1(NZ))THEN
    FWTC=AMIN1(1.0,0.667*WTRTs1(NZ)/WTLSs1(NZ))
  ELSE
    FWTC=1.0_r8
  ENDIF
  IF(WTRTs1(NZ).GT.ZEROPs1(NZ))THEN
    FWTS=AMIN1(1.0,WTLSs1(NZ)/(0.667*WTRTs1(NZ)))
  ELSE
    FWTS=1.0_r8
  ENDIF
!     ELSE
!     FWTC=1.0_r8
!     FWTS=1.0_r8
!     ENDIF
  DO 290 L=NUs1,NIs1(NZ)
    IF(RTNT(1).GT.ZEROPs1(NZ))THEN
      FWTR(L)=AMAX1(0.0,RLNT(1,L)/RTNT(1))
    ELSE
      FWTR(L)=1.0_r8
    ENDIF
290 CONTINUE
!     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
!     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
!
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!
  WTLSs1(NZ)=0._r8
  DO 309 NB=1,NBRs1(NZ)
    WTLSs1(NZ)=WTLSs1(NZ)+WTLSBs1(NB,NZ)
309 CONTINUE
!
!     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
!     OF TOTAL SINK STRENGTH OF THE CANOPY
!
!     IDTHB=branch living flag: 0=alive,1=dead
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAYs1(8,=end date for setting final seed number
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
  DO 310 NB=1,NBRs1(NZ)
    IF(IDTHBs1(NB,NZ).EQ.0)THEN
      IF(WTLSs1(NZ).GT.ZEROPs1(NZ))THEN
        FWTB(NB)=AMAX1(0.0,WTLSBs1(NB,NZ)/WTLSs1(NZ))
      ELSE
        FWTB(NB)=1.0_r8
      ENDIF
      IF(ISTYPs1(NZ).EQ.0)THEN
        PTSHTR=PTSHTs1(NZ)*PTRT**0.167
      ELSE
        PTSHTR=PTSHTs1(NZ)
      ENDIF
      DO 415 L=NUs1,NIs1(NZ)
        WTLSBX=WTLSBs1(NB,NZ)*FWODBs1(1)*FWTR(L)*FWTC
        WTRTLX=WTRTLs1(1,L,NZ)*FWODRs1(1)*FWTB(NB)*FWTS
        WTLSBB=AMAX1(0.0,WTLSBX,FSNK*WTRTLX)
        WTRTLR=AMAX1(0.0,WTRTLX,FSNK*WTLSBX)
        WTPLTT=WTLSBB+WTRTLR
        IF(WTPLTT.GT.ZEROPs1(NZ))THEN
          CPOOLB=AMAX1(0.0,CPOOLs1(NB,NZ)*FWTR(L))
          CPOOLS=AMAX1(0.0,CPOOLRs1(1,L,NZ)*FWTB(NB))
          CPOOLD=(CPOOLB*WTRTLR-CPOOLS*WTLSBB)/WTPLTT
          XFRC=PTSHTR*CPOOLD
          CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)-XFRC
          CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)+XFRC
          CPOOLT=CPOOLS+CPOOLB
          IF(CPOOLT.GT.ZEROPs1(NZ))THEN
            ZPOOLB=AMAX1(0.0,ZPOOLs1(NB,NZ)*FWTR(L))
            ZPOOLS=AMAX1(0.0,ZPOOLRs1(1,L,NZ)*FWTB(NB))
            ZPOOLD=(ZPOOLB*CPOOLS-ZPOOLS*CPOOLB)/CPOOLT
            XFRN=PTSHTR*ZPOOLD
            PPOOLB=AMAX1(0.0,PPOOLs1(NB,NZ)*FWTR(L))
            PPOOLS=AMAX1(0.0,PPOOLRs1(1,L,NZ)*FWTB(NB))
            PPOOLD=(PPOOLB*CPOOLS-PPOOLS*CPOOLB)/CPOOLT
            XFRP=PTSHTR*PPOOLD
          ELSE
            XFRN=0._r8
            XFRP=0._r8
          ENDIF
          ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)-XFRN
          ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)+XFRN
          PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)-XFRP
          PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)+XFRP

        ENDIF
415   CONTINUE
    ENDIF
310   CONTINUE
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
  integer :: N,L,K,NR
  REAL(R8) :: RTDPL(10,JZ1)
  real(r8) :: CUPRL,CUPRO,CUPRC
  real(r8) :: RTDPP,RTDPS,RTSKP
  real(r8) :: RTSKS

  associate(                                 &
    SDPTHs1      =>   plt_morph%SDPTHs1    , &
    NIs1         =>   plt_morph%NIs1       , &
    NRTs1        =>   plt_morph%NRTs1        &
  )

!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER

  RLNT=0._R8
  RTSK1=0._r8
  RTSK2=0._r8
  RTNT=0._r8
  DO 4995 N=1,MYs1(NZ)
    DO 4990 L=NUs1,NIs1(NZ)
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
      IF(VOLXs1(L).GT.ZEROS2s1)THEN
        CUPRL=0.86*(RUPNH4s1(N,L,NZ)+RUPNHBs1(N,L,NZ) &
          +RUPNO3s1(N,L,NZ)+RUPNOBs1(N,L,NZ)+RUPH2Ps1(N,L,NZ) &
          +RUPH2Bs1(N,L,NZ)+RUPH1Ps1(N,L,NZ)+RUPH1Bs1(N,L,NZ))
        CUPRO=0.86*(RUONH4s1(N,L,NZ)+RUONHBs1(N,L,NZ) &
          +RUONO3s1(N,L,NZ)+RUONOBs1(N,L,NZ)+RUOH2Ps1(N,L,NZ) &
          +RUOH2Bs1(N,L,NZ)+RUOH1Ps1(N,L,NZ)+RUOH1Bs1(N,L,NZ))
        CUPRC=0.86*(RUCNH4s1(N,L,NZ)+RUCNHBs1(N,L,NZ) &
          +RUCNO3s1(N,L,NZ)+RUCNOBs1(N,L,NZ)+RUCH2Ps1(N,L,NZ) &
          +RUCH2Bs1(N,L,NZ)+RUCH1Ps1(N,L,NZ)+RUCH1Bs1(N,L,NZ))
!
!     ACCUMULATE RESPIRATION IN FLUX ARRAYS
!
!     RCO2A=total root respiration
!     RCO2M,RCO2N=RCO2A unltd by O2,nonstructural C
!     CUPRL=C respiration for nutrient uptake
!     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
!     CPOOLR=non-structural C mass in root
!
        RCO2Ms1(N,L,NZ)=RCO2Ms1(N,L,NZ)+CUPRO
        RCO2Ns1(N,L,NZ)=RCO2Ns1(N,L,NZ)+CUPRC
        RCO2As1(N,L,NZ)=RCO2As1(N,L,NZ)-CUPRL
        CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)-CUPRL
!
!     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
!     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
!
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
!     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
!     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
!
        DO 195 K=0,jcplx11
          CPOOLRs1(N,L,NZ)=CPOOLRs1(N,L,NZ)+RDFOMCs1(N,K,L,NZ)
          ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ)+RDFOMNs1(N,K,L,NZ)
          PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ)+RDFOMPs1(N,K,L,NZ)
195     CONTINUE
        ZPOOLRs1(N,L,NZ)=ZPOOLRs1(N,L,NZ) &
          +RUPNH4s1(N,L,NZ)+RUPNHBs1(N,L,NZ) &
          +RUPNO3s1(N,L,NZ)+RUPNOBs1(N,L,NZ)
        PPOOLRs1(N,L,NZ)=PPOOLRs1(N,L,NZ) &
          +RUPH2Ps1(N,L,NZ)+RUPH2Bs1(N,L,NZ) &
          +RUPH1Ps1(N,L,NZ)+RUPH1Bs1(N,L,NZ)
!
!     GROWTH OF EACH ROOT AXIS
!
        DO 4985 NR=1,NRTs1(NZ)
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
            IF(RTDP1s1(N,NR,NZ).GT.CDPTHZs1(L-1))THEN
              IF(RTDP1s1(N,NR,NZ).LE.CDPTHZs1(L))THEN
                RTDPP=RTDP1s1(N,NR,NZ)+HTSTZs1(NZ)
                RTSK1(N,L,NR)=RTSK(IGTYPs1(NZ))*XRTN1*RRAD1s1(N,L,NZ)**2/RTDPP
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
            RTDPL(NR,L)=AMAX1(0.0,RTDP1s1(1,NR,NZ)-CDPTHZs1(L-1)-RTDPX)
            RTDPL(NR,L)=AMAX1(0.0,AMIN1(DLYR3s1(L),RTDPL(NR,L)) &
              -AMAX1(0.0,SDPTHs1(NZ)-CDPTHZs1(L-1)-HTCTLs1(NZ)))
            RTDPS=AMAX1(SDPTHs1(NZ),CDPTHZs1(L-1))+0.5*RTDPL(NR,L)+HTSTZs1(NZ)
            IF(RTDPS.GT.ZEROs1)THEN
              RTSKP=XRTN1*RRAD1s1(N,L,NZ)**2/RTDPS
              RTSKS=safe_adb(RTN2s1(N,L,NR,NZ)*RRAD2s1(N,L,NZ)**2,RTLGAs1(N,L,NZ))
              IF(RTSKP+RTSKS.GT.ZEROPs1(NZ))THEN
                RTSK2(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
              ELSE
                RTSK2(N,L,NR)=0._r8
              ENDIF
            ELSE
              RTSK2(N,L,NR)=0._r8
            ENDIF
          ELSE
            RTSK2(N,L,NR)=safe_adb(RTN2s1(N,L,NR,NZ)*RRAD2s1(N,L,NZ)**2,RTLGAs1(N,L,NZ))
          ENDIF
          RTNT(N)=RTNT(N)+RTSK2(N,L,NR)
          RLNT(N,L)=RLNT(N,L)+RTSK2(N,L,NR)
4985    CONTINUE
      ENDIF
4990  CONTINUE
4995  CONTINUE
  end associate
  end subroutine SummarizeRootSink

end module RootMod
