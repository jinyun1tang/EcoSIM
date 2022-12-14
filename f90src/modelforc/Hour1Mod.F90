module Hour1Mod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod  , only : test_aeqb,AZMAX1,AZMIN1
  use abortutils   , only : endrun, print_info
  use ATSUtilsMod
  use TracerPropMod
  use TracerIDMod
  use EcoSimConst
  use MiniFuncMod
  use SoilHydroParaMod
  use PlantAPI, only : PlantCanopyRadsModel
  use MicrobialDataType
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use CanopyRadDataType
  use EcoSIMSolverPar
  use ChemTracerParsMod
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use LandSurfDataType
  use PlantTraitDataType
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use CanopyDataType
  use EcoSimSumDataType
  use RootDataType
  use EcosimBGCFluxType
  use AqueChemDatatype
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use SedimentDataType
  use PlantDataRateType
  use GridDataType
  use EcoSIMConfig, only : jcplx1 => jcplx1c,jcplx=>jcplxc,nlbiomcp=>nlbiomcpc
  use EcoSIMConfig, only : ndbiomcp=>ndbiomcpc,jsken=>jskenc,NFGs=>NFGsc,do_instequil
  use EcoSiMParDataMod, only : micpar,pltpar
  implicit none

  private

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  real(r8) :: BKDSX  !     BKDSX=maximm soil bulk density
  real(r8) :: THETPW !     THETPW=minimum air-filled porosity for saturation (m3 m-3)
  real(r8) :: THETWP
  real(r8) :: XVOLWC(0:3)
  real(r8), pointer :: THETRX(:)
!
!     XVOLWC=foliar water retention capacity (m3 m-2)
!     THETRX=litter water retention capacity (m3 g C-1)

  public :: hour1
  public :: InitHour1
  contains

  subroutine InitHour1(n_litrsfk)

  implicit none
  integer, intent(in) :: n_litrsfk

  allocate(THETRX(1:n_litrsfk))
  BKDSX=1.89_r8

  THETPW=0.01_r8
  THETWP=1.0_r8-THETPW

  XVOLWC=real((/5.0E-04,2.5E-04,2.5E-04,2.5E-04/),r8)
  THETRX=real((/4.0E-06,8.0E-06,8.0E-06/),r8)

  end subroutine InitHour1
!------------------------------------------------------------------------------------------

  SUBROUTINE hour1(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE REINITIALIZES HOURLY VARIABLES USED IN OTHER
!     SUBROUTINES
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: L,NX,NY
  real(r8) :: THETPZ(JZ,JY,JX)
  real(r8) :: DPTH0(JY,JX)
  real(r8) :: FVOLR
  real(r8) :: VOLWCX
  real(r8) :: VOLWRX0,VOLR0
  real(r8) :: XJ
  integer :: NZ,NR,K
!     execution begins here

  XJ=J
  DOY=I-1+XJ/24
!
  call ResetLndscapeAccumlators()

  call SetAtmosTracerConcentration(I,NHW,NHE,NVN,NVS)
!
!     RESET FLUX ARRAYS USED IN OTHER SUBROUTINES
!
  call ResetFluxArrays(I,NHW,NHE,NVN,NVS)
!
!     IF SALT FLAG SET
!
  IF(ISALTG.NE.0)THEN
    call ResetSaltModelArrays(NHW,NHE,NVN,NVS)
  ENDIF
!
!     RESET SOIL PROPERTIES AND PEDOTRANSFER FUNCTIONS
!     FOLLOWING ANY SOIL DISTURBANCE
!
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      IF(J.EQ.1)THEN
        IFLGT(NY,NX)=0
        DO  NZ=1,NP(NY,NX)
          PSILZ(NZ,NY,NX)=0.0_r8
        ENDDO
      ENDIF
!
!     HYDROLOGICAL PRPOERTIES OR SURFACE LITTER
!
!     VOLWRX=liter water holding capacity
!     VOLR=dry litter volume
!     POROS0,FC,WP=litter porosity,field capacity,wilting point
!
      VOLWRX0=0._r8
      VOLR0=0._r8
      DO K=1,micpar%n_litrsfk
        VOLWRX0=VOLWRX0+THETRX(K)*RC0(K,NY,NX)
        VOLR0=VOLR0+RC0(K,NY,NX)/BKRS(K)
      ENDDO

      VOLWRX(NY,NX)=AZMAX1(VOLWRX0)
      VOLR(NY,NX)=AZMAX1(VOLR0*ppmc)

      IF(VOLR(NY,NX).GT.ZEROS(NY,NX))THEN
        FVOLR=VOLWRX(NY,NX)/VOLR(NY,NX)
      ELSE
        FVOLR=THETRX(micpar%k_fine_litr)/BKRS(micpar%k_fine_litr)
      ENDIF
      POROS0(NY,NX)=FVOLR
      FC(0,NY,NX)=0.500_r8*FVOLR
      WP(0,NY,NX)=0.125_r8*FVOLR
      PSL(0,NY,NX)=LOG(POROS0(NY,NX))
      FCL(0,NY,NX)=LOG(FC(0,NY,NX))
      WPL(0,NY,NX)=LOG(WP(0,NY,NX))
      PSD(0,NY,NX)=PSL(0,NY,NX)-FCL(0,NY,NX)
      FCD(0,NY,NX)=FCL(0,NY,NX)-WPL(0,NY,NX)
      SRP(0,NY,NX)=1.00_r8
!
!     RESET SURFACE LITTER PHYSICAL PROPERTIES (DENSITY, TEXTURE)
!     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)

      call SetLiterSoilPropAftDisturbance(I,J,NY,NX)
!
!
!     PARAMETERS FOR COHESION, EROSIVITY, AND ROUGHNESS OF SURFACE SOIL USED
!     FOR SURFACE WATER AND SEDIMENT TRANSPORT IN 'EROSION'

      call SetSurfaceProperty4SedErosion(NY,NX)
!
!     'RESET HOURLY ACCUMULATORS'
      call SetHourlyAccumulators(NY,NX)
!
!     RESET ARRAYS TO TRANSFER MATERIALS WITHIN SOILS
!     AND BETWEEN SOILS AND PLANTS
!
      call SetArrays4PlantSoilTransfer(NY,NX)
!
!     IF SOC FLAG IS SET
!

      IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
        call UpdateTotalSOC(NY,NX)
      ENDIF

      call ZeroHourlyArrays(NY,NX)

      call GetChemicalConcsInSoil(NY,NX,THETPZ)

      call GetSoluteConcentrations(NY,NX)

      call Prep4PlantMicrobeUptake(NY,NX)

      call CalGasSolubility(NY,NX)

      call GetSoilHydraulicVars(NY,NX)

!     CALCULATE ACTIVE LAYER DEPTH
      call DiagActiveLayerDepth(NY,NX)

!     OUTPUT FOR WATER TABLE DEPTH
      call GetOutput4WaterTableDepth(NY,NX,THETPZ)

      call GetSurfResidualProperties(NY,NX,DPTH0)

      call SetTracerPropertyInLiterAir(NY,NX)

      if(do_instequil)call ForceGasAquaEquil(NY,NX)
!
!      call CanopyConditionModel(I,J,NY,NX,DPTH0)

      call PlantCanopyRadsModel(I,J,NY,NX,DPTH0(NY,NX))
!
!     RESET HOURLY INDICATORS
!
      THRMCX(NY,NX)=THRMC(NY,NX)
      THRMGX(NY,NX)=THRMG(NY,NX)
      CNETX(NY,NX)=TCNET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      THRMC(NY,NX)=0.0_r8
      THRMG(NY,NX)=0.0_r8
      TLEX(NY,NX)=TLEC(NY,NX)
      TSHX(NY,NX)=TSHC(NY,NX)
      TLEC(NY,NX)=0.0_r8
      TSHC(NY,NX)=0.0_r8
      TRN(NY,NX)=0.0_r8
      TLE(NY,NX)=0.0_r8
      TSH(NY,NX)=0.0_r8
      TGH(NY,NX)=0.0_r8
      TCCAN(NY,NX)=0.0_r8
      TCNET(NY,NX)=0.0_r8
      RECO(NY,NX)=0.0_r8
!
!     CANOPY RETENTION OF PRECIPITATION
!
!     XVOLWC=foliar surface water retention capacity
!     ARLFP,ARSTP=leaf,stalk area of PFT
!     FLWC,TFLWC=water retention of PFT,combined canopy
!     PRECA=precipitation+irrigation
!     FRADP=fraction of radiation received by each PFT canopy
!     VOLWC=canopy surface water retention
!
      DO  NZ=1,NP(NY,NX)
        VOLWCX=XVOLWC(IGTYP(NZ,NY,NX))*(ARLFP(NZ,NY,NX)+ARSTP(NZ,NY,NX))
        FLWC(NZ,NY,NX)=AZMAX1(AMIN1(PRECA(NY,NX)*FRADP(NZ,NY,NX) &
          ,VOLWCX-VOLWC(NZ,NY,NX)))
        TFLWCI(NY,NX)=TFLWCI(NY,NX)+PRECA(NY,NX)*FRADP(NZ,NY,NX)
        TFLWC(NY,NX)=TFLWC(NY,NX)+FLWC(NZ,NY,NX)
!
!     NUMBERS OF TOP AND BOTTOM ROOTED SOIL LAYERS
!
!     NG=number of uppermost rooted layer
!     NINR=number of lowest rooted layer
!
        NG(NZ,NY,NX)=MAX(NG(NZ,NY,NX),NU(NY,NX))
        NIX(NZ,NY,NX)=MAX(NIX(NZ,NY,NX),NU(NY,NX))
        DO  NR=1,JC
          NINR(NR,NZ,NY,NX)=MAX(NINR(NR,NZ,NY,NX),NU(NY,NX))
        ENDDO
      ENDDO
!
!     WRITE SW AND PAR ALBEDO

    ENDDO
  ENDDO
!
!     FERTILIZER APPLICATIONS OCCUR AT SOLAR NOON
  call ApplyFertilizerAtNoon(I,J,NHW,NHE,NVN,NVS)

  END subroutine hour1
!------------------------------------------------------------------------------------------
  subroutine ResetLndscapeAccumlators()
!     RESET HOURLY SOIL ACCUMULATORS FOR WATER, HEAT, GASES, SOLUTES
!
  implicit none
  VOLWSO=0.0_r8
  HEATSO=0.0_r8
  OXYGSO=0.0_r8
  TLH2G=0.0_r8
  TSEDSO=0.0_r8
  TLRSDC=0.0_r8
  TLORGC=0.0_r8
  TLCO2G=0.0_r8
  TLRSDN=0.0_r8
  TLORGN=0.0_r8
  TLN2G=0.0_r8
  TLRSDP=0.0_r8
  TLORGP=0.0_r8
  TLNH4=0.0_r8
  TLNO3=0.0_r8
  TLPO4=0.0_r8
  TION=0.0_r8
  TBALE=0.0_r8
  end subroutine ResetLndscapeAccumlators
!------------------------------------------------------------------------------------------

  subroutine SetAtmosTracerConcentration(I,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS

  integer :: NX, NY
!     begin_execution
!
!     CONCENTRATIONS OF CO2, CH4, O2, N2, N2O, NH3, H2 IN ATMOSPHERE,
!     PRECIPITATION AND IRRIGATION FROM MIXING RATIOS READ IN 'READS'
!
!     C*E,C*R,C*Q=atmospheric,precipitation,irrigation solute concentrations
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
! CSTRR: surface irrigation ion strength, [g m-3]
  DO NX=NHW,NHE
    DO NY=NVN,NVS

      AtmGgms(idg_CO2,NY,NX)=CO2E(NY,NX)*5.36E-04_r8*TREF/TKA(NY,NX)  !gC/m3
      AtmGgms(idg_CH4,NY,NX)=CH4E(NY,NX)*5.36E-04_r8*TREF/TKA(NY,NX)  !gC/m3
      AtmGgms(idg_O2 ,NY,NX)=OXYE(NY,NX)*1.43E-03_r8*TREF/TKA(NY,NX)  !gO/m3
      AtmGgms(idg_N2 ,NY,NX)=Z2GE(NY,NX)*1.25E-03_r8*TREF/TKA(NY,NX)  !gN/m3
      AtmGgms(idg_N2O,NY,NX)=Z2OE(NY,NX)*1.25E-03_r8*TREF/TKA(NY,NX)  !gN/m3
      AtmGgms(idg_NH3,NY,NX)=ZNH3E(NY,NX)*6.25E-04_r8*TREF/TKA(NY,NX) !gN/m3
      AtmGgms(idg_H2 ,NY,NX)=H2GE(NY,NX)*8.92E-05_r8*TREF/TKA(NY,NX)  !gH/m3
      AtmGgms(idg_NH3B,NY,NX)=ZNH3E(NY,NX)*6.25E-04_r8*TREF/TKA(NY,NX) !gN/m3

      CCOR(NY,NX)=AtmGgms(idg_CO2,NY,NX)*gas_solubility(idg_CO2,TCA(NY,NX)) &
         /(EXP(ACTCG(idg_CO2)*CSTRR(NY,NX)))
      CCHR(NY,NX)=AtmGgms(idg_CH4,NY,NX)*gas_solubility(idg_CH4,TCA(NY,NX)) &
        /(EXP(ACTCG(idg_CH4)*CSTRR(NY,NX)))
      COXR(NY,NX)=AtmGgms(idg_O2,NY,NX)*gas_solubility(idg_O2, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_O2)*CSTRR(NY,NX)))
      CNNR(NY,NX)=AtmGgms(idg_N2,NY,NX)*gas_solubility(idg_N2, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_N2)*CSTRR(NY,NX)))
      CN2R(NY,NX)=AtmGgms(idg_N2O,NY,NX)*gas_solubility(idg_N2O, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_N2O)*CSTRR(NY,NX)))

      CCOQ(NY,NX)=AtmGgms(idg_CO2,NY,NX)*gas_solubility(idg_CO2, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_CO2)*CSTRQ(I,NY,NX)))
      CCHQ(NY,NX)=AtmGgms(idg_CH4,NY,NX)*gas_solubility(idg_CH4, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_CH4)*CSTRQ(I,NY,NX)))
      COXQ(NY,NX)=AtmGgms(idg_O2,NY,NX)*gas_solubility(idg_O2, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_O2)*CSTRQ(I,NY,NX)))
      CNNQ(NY,NX)=AtmGgms(idg_N2,NY,NX)*gas_solubility(idg_N2, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_N2)*CSTRQ(I,NY,NX)))
      CN2Q(NY,NX)=AtmGgms(idg_N2O,NY,NX)*gas_solubility(idg_N2O, TCA(NY,NX)) &
        /(EXP(ACTCG(idg_N2O)*CSTRQ(I,NY,NX)))

    ENDDO
  ENDDO
  end subroutine SetAtmosTracerConcentration
!------------------------------------------------------------------------------------------

  subroutine ResetFluxArrays(I,NHW,NHE,NVN,NVS)
!
  use EcoSIMConfig, only : column_mode
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS

  integer :: L,N,NX,NY,K,NN,NO,M,NGL
  integer :: extragrid
!     begin_execution
  extragrid=1
  if(column_mode)extragrid=0

  DO  NX=NHW,NHE+extragrid
    DO  NY=NVN,NVS+extragrid
!
!     WATER,SNOW,SOLUTE RUNOFF
!

      QR(1:2,1:2,NY,NX)=0.0_r8
      HQR(1:2,1:2,NY,NX)=0.0_r8

      XOCQRS(1:jcplx,1:2,1:2,NY,NX)=0.0_r8
      XONQRS(1:jcplx,1:2,1:2,NY,NX)=0.0_r8
      XOPQRS(1:jcplx,1:2,1:2,NY,NX)=0.0_r8
      XOAQRS(1:jcplx,1:2,1:2,NY,NX)=0.0_r8

      trcg_XRS(idg_beg:idg_end-1,1:2,1:2,NY,NX)=0.0_r8
      trcn_XRS(ids_nut_beg:ids_nuts_end,1:2,1:2,NY,NX)=0.0_r8

      QS(1:2,NY,NX)=0.0_r8
      QW(1:2,NY,NX)=0.0_r8
      QI(1:2,NY,NX)=0.0_r8
      HQS(1:2,NY,NX)=0.0_r8
      XCOQSS(1:2,NY,NX)=0.0_r8
      XCHQSS(1:2,NY,NX)=0.0_r8
      XOXQSS(1:2,NY,NX)=0.0_r8
      XNGQSS(1:2,NY,NX)=0.0_r8
      XN2QSS(1:2,NY,NX)=0.0_r8
      XN4QSS(1:2,NY,NX)=0.0_r8
      XN3QSS(1:2,NY,NX)=0.0_r8
      XNOQSS(1:2,NY,NX)=0.0_r8
      XP1QSS(1:2,NY,NX)=0.0_r8
      XP4QSS(1:2,NY,NX)=0.0_r8
!
!
!     GAS AND SOLUTE FLUXES
!
      DO  L=0,NL(NY,NX)+1

        trcs_XFLS(idg_beg:idg_end-1,1:3,L,NY,NX)=0.0_r8

        trcs_XFLS(ids_nut_beg:ids_nuts_end,1:3,L,NY,NX)=0.0_r8

        XOCFLS(1:jcplx,1:3,L,NY,NX)=0.0_r8
        XONFLS(1:jcplx,1:3,L,NY,NX)=0.0_r8
        XOPFLS(1:jcplx,1:3,L,NY,NX)=0.0_r8
        XOAFLS(1:jcplx,1:3,L,NY,NX)=0.0_r8
      ENDDO
!
!     BAND AND MACROPORE FLUXES
!
      DO L=1,NL(NY,NX)+1
        FLW(1:3,L,NY,NX)=0.0_r8
        FLWX(1:3,L,NY,NX)=0.0_r8
        FLWH(1:3,L,NY,NX)=0.0_r8
        HFLW(1:3,L,NY,NX)=0.0_r8

        trcs_XFLS(ids_beg:ids_end,1:3,L,NY,NX)=0.0_r8
        R3GasADTFlx(idg_beg:idg_end,1:3,L,NY,NX)=0._r8

        XOCFHS(1:jcplx,1:3,L,NY,NX)=0.0_r8
        XONFHS(1:jcplx,1:3,L,NY,NX)=0.0_r8
        XOPFHS(1:jcplx,1:3,L,NY,NX)=0.0_r8
        XOAFHS(1:jcplx,1:3,L,NY,NX)=0.0_r8

      ENDDO
    ENDDO
  ENDDO

!     IF EROSION FLAG SET
!
  IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
    DO NX=NHW,NHE+extragrid
      DO NY=NVN,NVS+extragrid
        XSEDER(1:2,1:2,NY,NX)=0.0_r8
        XSANER(1:2,1:2,NY,NX)=0.0_r8
        XSILER(1:2,1:2,NY,NX)=0.0_r8
        XCLAER(1:2,1:2,NY,NX)=0.0_r8
        XCECER(1:2,1:2,NY,NX)=0.0_r8
        XAECER(1:2,1:2,NY,NX)=0.0_r8
        XNH4ER(1:2,1:2,NY,NX)=0.0_r8
        XNH3ER(1:2,1:2,NY,NX)=0.0_r8
        XNHUER(1:2,1:2,NY,NX)=0.0_r8
        XNO3ER(1:2,1:2,NY,NX)=0.0_r8
        XNH4EB(1:2,1:2,NY,NX)=0.0_r8
        XNH3EB(1:2,1:2,NY,NX)=0.0_r8
        XNHUEB(1:2,1:2,NY,NX)=0.0_r8
        XNO3EB(1:2,1:2,NY,NX)=0.0_r8
        XN4ER(1:2,1:2,NY,NX)=0.0_r8
        XNBER(1:2,1:2,NY,NX)=0.0_r8
        XHYER(1:2,1:2,NY,NX)=0.0_r8
        XALER(1:2,1:2,NY,NX)=0.0_r8
        XCAER(1:2,1:2,NY,NX)=0.0_r8
        XMGER(1:2,1:2,NY,NX)=0.0_r8
        XNAER(1:2,1:2,NY,NX)=0.0_r8
        XKAER(1:2,1:2,NY,NX)=0.0_r8
        XHCER(1:2,1:2,NY,NX)=0.0_r8
        XAL2ER(1:2,1:2,NY,NX)=0.0_r8
        XOH0ER(1:2,1:2,NY,NX)=0.0_r8
        XOH1ER(1:2,1:2,NY,NX)=0.0_r8
        XOH2ER(1:2,1:2,NY,NX)=0.0_r8
        XH1PER(1:2,1:2,NY,NX)=0.0_r8
        XH2PER(1:2,1:2,NY,NX)=0.0_r8
        XOH0EB(1:2,1:2,NY,NX)=0.0_r8
        XOH1EB(1:2,1:2,NY,NX)=0.0_r8
        XOH2EB(1:2,1:2,NY,NX)=0.0_r8
        XH1PEB(1:2,1:2,NY,NX)=0.0_r8
        XH2PEB(1:2,1:2,NY,NX)=0.0_r8
        PALOER(1:2,1:2,NY,NX)=0.0_r8
        PFEOER(1:2,1:2,NY,NX)=0.0_r8
        PCACER(1:2,1:2,NY,NX)=0.0_r8
        PCASER(1:2,1:2,NY,NX)=0.0_r8
        PALPER(1:2,1:2,NY,NX)=0.0_r8
        PFEPER(1:2,1:2,NY,NX)=0.0_r8
        PCPDER(1:2,1:2,NY,NX)=0.0_r8
        PCPHER(1:2,1:2,NY,NX)=0.0_r8
        PCPMER(1:2,1:2,NY,NX)=0.0_r8
        PALPEB(1:2,1:2,NY,NX)=0.0_r8
        PFEPEB(1:2,1:2,NY,NX)=0.0_r8
        PCPDEB(1:2,1:2,NY,NX)=0.0_r8
        PCPHEB(1:2,1:2,NY,NX)=0.0_r8
        PCPMEB(1:2,1:2,NY,NX)=0.0_r8
        OMCER(:,:,1:2,1:2,NY,NX)=0.0_r8
        OMNER(:,:,1:2,1:2,NY,NX)=0.0_r8
        OMPER(:,:,1:2,1:2,NY,NX)=0.0_r8

        OMCERff(:,1:2,1:2,NY,NX)=0.0_r8
        OMNERff(:,1:2,1:2,NY,NX)=0.0_r8
        OMPERff(:,1:2,1:2,NY,NX)=0.0_r8

        ORCER(:,:,1:2,1:2,NY,NX)=0.0_r8
        ORNER(:,:,1:2,1:2,NY,NX)=0.0_r8
        ORPER(:,:,1:2,1:2,NY,NX)=0.0_r8
        OHCER(:,1:2,1:2,NY,NX)=0.0_r8
        OHNER(:,1:2,1:2,NY,NX)=0.0_r8
        OHPER(:,1:2,1:2,NY,NX)=0.0_r8
        OSCER(:,:,1:2,1:2,NY,NX)=0.0_r8
        OSAER(:,:,1:2,1:2,NY,NX)=0.0_r8
        OSNER(:,:,1:2,1:2,NY,NX)=0.0_r8
        OSPER(:,:,1:2,1:2,NY,NX)=0.0_r8
      ENDDO
    ENDDO
  ENDIF

  end subroutine ResetFluxArrays
!------------------------------------------------------------------------------------------

  subroutine ResetSaltModelArrays(NHW,NHE,NVN,NVS)
!
  use EcoSIMConfig, only : column_mode
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: N,NX,NY,L,NN,NSA
  integer :: extragrid
!     begin_execution

  extragrid=1
  if(column_mode)extragrid=0
  DO  NX=NHW,NHE+extragrid
    DO  NY=NVN,NVS+extragrid

      XQRAL(1:2,1:2,NY,NX)=0.0_r8
      XQRFE(1:2,1:2,NY,NX)=0.0_r8
      XQRHY(1:2,1:2,NY,NX)=0.0_r8
      XQRCA(1:2,1:2,NY,NX)=0.0_r8
      XQRMG(1:2,1:2,NY,NX)=0.0_r8
      XQRNA(1:2,1:2,NY,NX)=0.0_r8
      XQRKA(1:2,1:2,NY,NX)=0.0_r8
      XQROH(1:2,1:2,NY,NX)=0.0_r8
      XQRSO(1:2,1:2,NY,NX)=0.0_r8
      XQRCL(1:2,1:2,NY,NX)=0.0_r8
      XQRC3(1:2,1:2,NY,NX)=0.0_r8
      XQRHC(1:2,1:2,NY,NX)=0.0_r8
      XQRAL1(1:2,1:2,NY,NX)=0.0_r8
      XQRAL2(1:2,1:2,NY,NX)=0.0_r8
      XQRAL3(1:2,1:2,NY,NX)=0.0_r8
      XQRAL4(1:2,1:2,NY,NX)=0.0_r8
      XQRALS(1:2,1:2,NY,NX)=0.0_r8
      XQRFE1(1:2,1:2,NY,NX)=0.0_r8
      XQRFE2(1:2,1:2,NY,NX)=0.0_r8
      XQRFE3(1:2,1:2,NY,NX)=0.0_r8
      XQRFE4(1:2,1:2,NY,NX)=0.0_r8
      XQRFES(1:2,1:2,NY,NX)=0.0_r8
      XQRCAO(1:2,1:2,NY,NX)=0.0_r8
      XQRCAC(1:2,1:2,NY,NX)=0.0_r8
      XQRCAH(1:2,1:2,NY,NX)=0.0_r8
      XQRCAS(1:2,1:2,NY,NX)=0.0_r8
      XQRMGO(1:2,1:2,NY,NX)=0.0_r8
      XQRMGC(1:2,1:2,NY,NX)=0.0_r8
      XQRMGH(1:2,1:2,NY,NX)=0.0_r8
      XQRMGS(1:2,1:2,NY,NX)=0.0_r8
      XQRNAC(1:2,1:2,NY,NX)=0.0_r8
      XQRNAS(1:2,1:2,NY,NX)=0.0_r8
      XQRKAS(1:2,1:2,NY,NX)=0.0_r8
      XQRH0P(1:2,1:2,NY,NX)=0.0_r8
      XQRH3P(1:2,1:2,NY,NX)=0.0_r8
      XQRF1P(1:2,1:2,NY,NX)=0.0_r8
      XQRF2P(1:2,1:2,NY,NX)=0.0_r8
      XQRC0P(1:2,1:2,NY,NX)=0.0_r8
      XQRC1P(1:2,1:2,NY,NX)=0.0_r8
      XQRC2P(1:2,1:2,NY,NX)=0.0_r8
      XQRM1P(1:2,1:2,NY,NX)=0.0_r8

      trcsa_XQS(idsa_beg:idsa_end,1:2,NY,NX)=0.0_r8

      DO  L=1,NL(NY,NX)+1
        DO NSA=idsa_beg,idsab_end
          trcsa_XFLS(NSA,1:3,L,NY,NX)=0.0_r8
          trcsa_XFHS(NSA,1:3,L,NY,NX)=0.0_r8
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine ResetSaltModelArrays
!------------------------------------------------------------------------------------------

  subroutine SetSoilPropertyAftDisturbance(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: VMINL,VSAND
  real(r8) :: PSISK(0:100)
  real(r8) :: THETK(100)
  real(r8) :: CORGCM
  real(r8) :: PTDS
  real(r8) :: SUM2,SUM1
  real(r8) :: VORGC
  real(r8) :: XK,YK
  integer :: L,K,N,M

  !     begin_execution
  D9975: DO L=NUI(NY,NX),NLI(NY,NX)
    !
    !     AREA,DLYR=lateral(1,2), vertical(3) area,thickness of soil layer
    !     VOLT,VOLX,VOLY=layer volume including,excluding rock,macropores
    !
    IF(BKDS(L,NY,NX).LE.ZERO.AND.DLYR(3,L,NY,NX).LE.ZERO2)THEN
      VOLW(L,NY,NX)=0.0_r8
      VOLI(L,NY,NX)=0.0_r8
    ENDIF
    AREA(1,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(2,L,NY,NX)
    AREA(2,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(1,L,NY,NX)
    VOLT(L,NY,NX)=AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)

    VOLX(L,NY,NX)=AMAX1(VOLT(L,NY,NX)*FMPR(L,NY,NX),1.e-8_r8)
    IF(BKDS(L,NY,NX).LE.ZERO)THEN
      VOLY(L,NY,NX)=VOLX(L,NY,NX)
    ENDIF
    !
    !     BKVL=soil mass
    !     C*=concentration,ORGC=SOC,SAND=sand,SILT=silt,CLAY=clay
    !     PTDS=particle density
    !     PTDSNU=particle density of surface layer for use in erosion.f
    !     POROS=porosity used in diffusivity
    !     VOLA,VOLW,VOLI,VOLP=total,water-,ice-,air-filled micropore volume
    !     VOLAH,VOLWH,VOLIH,VOLPH=total,water-,ice-,air-filled macropore volume
    !     EHUM=fraction of microbial decomposition product allocated to humus
    !     EPOC=fraction of SOC decomposition product allocated to POC
    !     SRP=parameter for deviation from linear log-log water retention
    !     PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
    !     PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
    !     FC,WP=water contents at field capacity,wilting point
    !     FCL,WPL=log FC,WP
    !     FCD,PSD=FCL-WPL,log(POROS)-FCL
    !
    BKVL(L,NY,NX)=BKDS(L,NY,NX)*VOLX(L,NY,NX)

    IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CORGC(L,NY,NX)=AMIN1(0.55E+06,ORGC(L,NY,NX)/BKVL(L,NY,NX))
      CSAND(L,NY,NX)=SAND(L,NY,NX)/BKVL(L,NY,NX)
      CSILT(L,NY,NX)=SILT(L,NY,NX)/BKVL(L,NY,NX)
      CCLAY(L,NY,NX)=CLAY(L,NY,NX)/BKVL(L,NY,NX)
    ELSE
      CORGC(L,NY,NX)=0.0_r8
      CSAND(L,NY,NX)=0.0_r8
      CSILT(L,NY,NX)=0.0_r8
      CCLAY(L,NY,NX)=0.0_r8
    ENDIF
    IF(BKVL(L,NY,NX).GT.ZERO)THEN
      CORGCM=AZMAX1(AMIN1(1.0_r8,MWC2Soil*CORGC(L,NY,NX)))
      PTDS=1.30_r8*CORGCM+2.66_r8*(1.0_r8-CORGCM)
      IF(L.EQ.NU(NY,NX))THEN
!surface layer
        POROS(L,NY,NX)=AMAX1(POROS(L,NY,NX),1.0_r8-(BKDS(L,NY,NX)/PTDS))
      ELSE
!  avoid float exception
        if(BKDS(L,NY,NX)>=PTDS)then
          POROS(L,NY,NX)=1.0_r8
        else
          POROS(L,NY,NX)=1.0_r8-(BKDS(L,NY,NX)/PTDS)
        endif
      ENDIF
    ELSE
      PTDS=0.0_r8
      POROS(L,NY,NX)=1.0_r8
    ENDIF
    !     VOLA(L,NY,NX)=AMAX1(POROS(L,NY,NX)*VOLY(L,NY,NX)
    !    2,VOLW(L,NY,NX)+VOLI(L,NY,NX))
    !     VOLAH(L,NY,NX)=AMAX1(FHOL(L,NY,NX)*VOLT(L,NY,NX)
    !    2,VOLWH(L,NY,NX)+VOLIH(L,NY,NX))
    VOLA(L,NY,NX)=POROS(L,NY,NX)*VOLY(L,NY,NX)
    VOLAH(L,NY,NX)=FHOL(L,NY,NX)*VOLT(L,NY,NX)
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLP(L,NY,NX)=AZMAX1(VOLA(L,NY,NX)-VOLW(L,NY,NX) &
        -VOLI(L,NY,NX))+AZMAX1(VOLAH(L,NY,NX)-VOLWH(L,NY,NX) &
        -VOLIH(L,NY,NX))
    ELSE
      VOLP(L,NY,NX)=0.0_r8
    ENDIF
    EHUM(L,NY,NX)=0.200_r8+0.333_r8*AMIN1(0.5_r8,CCLAY(L,NY,NX))

    EPOC(L,NY,NX)=1.0_r8
    call SoilHydroProperty(L,NY,NX,I,J)
!
!     SOIL HEAT CAPACITY AND THERMAL CONDUCTIVITY OF SOLID PHASE
!     FROM SOC AND TEXTURE
!
!     VORGC,VMINL,VSAND=volume fractions of SOC,mineral,sand
!     STC,DTC=weighted thermal conductivity of soil solid component
!
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VORGC=CORGCM*BKDS(L,NY,NX)/PTDS
      VMINL=(CSILT(L,NY,NX)+CCLAY(L,NY,NX))*BKDS(L,NY,NX)/PTDS
      VSAND=CSAND(L,NY,NX)*BKDS(L,NY,NX)/PTDS
      STC(L,NY,NX)=(1.253_r8*VORGC*9.050E-04_r8+0.514_r8*VMINL*1.056E-02_r8 &
        +0.386_r8*VSAND*2.112E-02_r8)*FMPR(L,NY,NX) &
        +0.514_r8*ROCK(L,NY,NX)*1.056E-02_r8
      DTC(L,NY,NX)=(1.253_r8*VORGC+0.514_r8*VMINL+0.386_r8*VSAND) &
        *FMPR(L,NY,NX)+0.514_r8*ROCK(L,NY,NX)
    ELSE
      STC(L,NY,NX)=0.0_r8
      DTC(L,NY,NX)=0.0_r8
    ENDIF
  ENDDO D9975
  end subroutine SetSoilPropertyAftDisturbance
!------------------------------------------------------------------------------------------

  subroutine ResetSurfResidualProperty(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: TVOLG0,TVOLWI
  real(r8) :: THETWR
  real(r8) :: VOLIRZ,VOLWRZ
  real(r8) :: XVOLW0,XVOLI0
!     begin_execution
!
!     FCR=litter water content at -0.01 MPa
!     THETY=litter hygroscopic water content
!
  CORGC(0,NY,NX)=0.55E+06_r8
!
!     SOIL SURFACE WATER STORAGE CAPACITY
!
!     IDTBL=water table flag from site file
!     DTBLX,DTBLZ=current,initial natural water table depth
!     DTBLY,DTBLD=current,initial artificial water table depth
!     ZS,ZW=soil,water surface roughness
!     VOLWD=soil surface water retention capacity
!     VOLWG=VOLWD accounting for above-ground water table
!     EHUM=fraction of microbial decompn product allocated to surface humus
!     EPOC=fraction of SOC decomposition product allocated to surface POC
!
  IF(IDTBL(NY,NX).LE.1.OR.IDTBL(NY,NX).EQ.3)THEN
    DTBLX(NY,NX)=DTBLZ(NY,NX)
  ELSEIF(IDTBL(NY,NX).EQ.2.OR.IDTBL(NY,NX).EQ.4)THEN
    DTBLX(NY,NX)=DTBLZ(NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
  ENDIF
  IF(IDTBL(NY,NX).EQ.3.OR.IDTBL(NY,NX).EQ.4)THEN
    DTBLY(NY,NX)=DTBLD(NY,NX)
  ENDIF

  IF(BKDS(NU(NY,NX),NY,NX).GT.ZERO)THEN
    ZS(NY,NX)=0.020_r8
  ELSE
    ZS(NY,NX)=ZW
  ENDIF
  VOLWD(NY,NX)=AMAX1(0.001_r8,0.112_r8*ZS(NY,NX)+3.10_r8*ZS(NY,NX)**2_r8 &
    -0.012_r8*ZS(NY,NX)*SLOPE(0,NY,NX))*AREA(3,NU(NY,NX),NY,NX)
  VOLWG(NY,NX)=AMAX1(VOLWD(NY,NX) &
    ,-(DTBLX(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))*AREA(3,NU(NY,NX),NY,NX))

  DPTH(NU(NY,NX),NY,NX)=CDPTH(NU(NY,NX),NY,NX)-0.5_r8*DLYR(3,NU(NY,NX),NY,NX)
  IF(BKVL(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)/BKVL(NU(NY,NX),NY,NX)
    CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)/BKVL(NU(NY,NX),NY,NX)
    CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)/BKVL(NU(NY,NX),NY,NX)
  ELSE
    CCLAY(NU(NY,NX),NY,NX)=0.0_r8
    CSILT(NU(NY,NX),NY,NX)=0.0_r8
    CSAND(NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  EHUM(0,NY,NX)=0.200_r8+0.333_r8*AMIN1(0.5_r8,CCLAY(NU(NY,NX),NY,NX))
  EPOC(0,NY,NX)=0.150_r8
  end subroutine ResetSurfResidualProperty
!------------------------------------------------------------------------------------------

  subroutine SetLiterSoilPropAftDisturbance(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  real(r8) :: PSISK(0:100),THETK(100)
  REAL(R8) :: SUM2,SUM1
  real(r8) :: XK,YK
  integer :: K,M
!     begin_execution
!     IFLGS=disturbance flag
!     BKDS,BKDSI=current,initial bulk density
!

  IF(IFLGS(NY,NX).NE.0)THEN

    call LitterHydroproperty(NY,NX)
!
!   'RESET SOIL PHYSICAL PROPERTIES (DENSITY, TEXTURE)'
!     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
!
    call SetSoilPropertyAftDisturbance(I,J,NY,NX)
!
!   'SURFACE RESIDUE PROPERTIES'
    call ResetSurfResidualProperty(NY,NX)
!
!     IFLGS=reset disturbance flag
!
    IFLGS(NY,NX)=0
  ENDIF

  end subroutine SetLiterSoilPropAftDisturbance
!------------------------------------------------------------------------------------------

  subroutine SetHourlyAccumulators(NY,NX)
!     implicit none
  integer, intent(in) :: NX,NY

  integer :: L
!     begin_execution

  UCO2S(NY,NX)=0.0_r8
  TOMT(NY,NX)=0.0_r8
  TONT(NY,NX)=0.0_r8
  TOPT(NY,NX)=0.0_r8
  UVOLW(NY,NX)=0.0_r8
  URSDC(NY,NX)=0.0_r8
  UORGC(NY,NX)=0.0_r8
  URSDN(NY,NX)=0.0_r8
  UORGN(NY,NX)=0.0_r8
  URSDP(NY,NX)=0.0_r8
  UORGP(NY,NX)=0.0_r8
  UNH4(NY,NX)=0.0_r8
  UNO3(NY,NX)=0.0_r8
  UPO4(NY,NX)=0.0_r8
  UPP4(NY,NX)=0.0_r8
  UION(NY,NX)=0.0_r8
  HVOLO(NY,NX)=0.0_r8
  HCO2G(NY,NX)=0.0_r8
  HCH4G(NY,NX)=0.0_r8
  HOXYG(NY,NX)=0.0_r8
  HN2GG(NY,NX)=0.0_r8
  HN2OG(NY,NX)=0.0_r8
  HNH3G(NY,NX)=0.0_r8
  FLWR(NY,NX)=0.0_r8
  HFLWR(NY,NX)=0.0_r8
  THAWR(NY,NX)=0.0_r8
  HTHAWR(NY,NX)=0.0_r8
  HEATI(NY,NX)=0.0_r8
  HEATS(NY,NX)=0.0_r8
  HEATE(NY,NX)=0.0_r8
  HEATV(NY,NX)=0.0_r8
  HEATH(NY,NX)=0.0_r8
  TEVAPG(NY,NX)=0.0_r8


  GasSfAtmFlx(idg_beg:idg_end,NY,NX)=0._r8
  XCODFR(NY,NX)=0.0_r8
  XCHDFR(NY,NX)=0.0_r8
  XOXDFR(NY,NX)=0.0_r8
  XNGDFR(NY,NX)=0.0_r8
  XN2DFR(NY,NX)=0.0_r8
  XN3DFR(NY,NX)=0.0_r8
  XHGDFR(NY,NX)=0.0_r8
  TVOLWP(NY,NX)=0.0_r8
  TVOLWC(NY,NX)=0.0_r8
  TFLWCI(NY,NX)=0.0_r8
  TFLWC(NY,NX)=0.0_r8
  TEVAPP(NY,NX)=0.0_r8
  TEVAPC(NY,NX)=0.0_r8
  THFLXC(NY,NX)=0.0_r8
  TENGYC(NY,NX)=0.0_r8

  TRFGas_root(idg_beg:idg_end-1,NY,NX)=0.0_r8
  ZESNC(NY,NX,:)=0.0_r8
  WTSTGT(NY,NX)=0.0_r8
  PPT(NY,NX)=0.0_r8
! zero arrays in the snow layers
  FLSW(1:JS,NY,NX)=0.0_r8
  FLSWH(1:JS,NY,NX)=0.0_r8
  HFLSW(1:JS,NY,NX)=0.0_r8
  FLSWR(1:JS,NY,NX)=0.0_r8
  HFLSWR(1:JS,NY,NX)=0.0_r8
  XFLWS(1:JS,NY,NX)=0.0_r8
  XFLWW(1:JS,NY,NX)=0.0_r8
  XFLWI(1:JS,NY,NX)=0.0_r8
  XHFLWW(1:JS,NY,NX)=0.0_r8
  XWFLXS(1:JS,NY,NX)=0.0_r8
  XWFLXI(1:JS,NY,NX)=0.0_r8
  XTHAWW(1:JS,NY,NX)=0.0_r8

  trcg_XBLS(idg_beg:idg_end-1,1:JS,NY,NX)=0.0_r8

  trcn_XBLS(ids_nut_beg:ids_nuts_end,1:JS,NY,NX)=0.0_r8

  IF(ISALTG.NE.0)THEN
    trcsa_XBLS(idsa_beg:idsa_end,1:JS,NY,NX)=0.0_r8
  ENDIF
  end subroutine SetHourlyAccumulators
!------------------------------------------------------------------------------------------

  subroutine SetArrays4PlantSoilTransfer(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

!     begin_execution

  ESNT(1:jsken,1:npelms,1:pltpar%n_pltlitrk,0:NL(NY,NX),NY,NX)=0.0_r8

  XOQCS(1:jcplx,0:NL(NY,NX),NY,NX)=0.0_r8
  XOQNS(1:jcplx,0:NL(NY,NX),NY,NX)=0.0_r8
  XOQPS(1:jcplx,0:NL(NY,NX),NY,NX)=0.0_r8
  XOQAS(1:jcplx,0:NL(NY,NX),NY,NX)=0.0_r8

  XZHYS(0:NL(NY,NX),NY,NX)=0.0_r8
  TRN4S(0:NL(NY,NX),NY,NX)=0.0_r8
  TRN3S(0:NL(NY,NX),NY,NX)=0.0_r8
  TRN3G(0:NL(NY,NX),NY,NX)=0.0_r8
  TRNO3(0:NL(NY,NX),NY,NX)=0.0_r8
  TRNO2(0:NL(NY,NX),NY,NX)=0.0_r8
  TRH1P(0:NL(NY,NX),NY,NX)=0.0_r8
  TRH2P(0:NL(NY,NX),NY,NX)=0.0_r8

  trcx_TR(idx_NH4,0:NL(NY,NX),NY,NX)=0.0_r8
  trcx_TR(idx_AEC+1:idx_anion_soil_end,0:NL(NY,NX),NY,NX)=0.0_r8

  trcp_TR(idsp_AlPO4,0:NL(NY,NX),NY,NX)=0.0_r8
  trcp_TR(idsp_FePO4,0:NL(NY,NX),NY,NX)=0.0_r8
  trcp_TR(idsp_CaHPO4,0:NL(NY,NX),NY,NX)=0.0_r8
  trcp_TR(idsp_HA,0:NL(NY,NX),NY,NX)=0.0_r8
  trcp_TR(idsp_CaH2PO4,0:NL(NY,NX),NY,NX)=0.0_r8

  TUPWTR(0:NL(NY,NX),NY,NX)=0.0_r8
  TUPHT(0:NL(NY,NX),NY,NX)=0.0_r8

  GasDisFlx(idg_beg:idg_end,0:NL(NY,NX),NY,NX)=0.0_r8

  TDFOMC(1:jcplx,NU(NY,NX):NL(NY,NX),NY,NX)=0.0_r8
  TDFOMN(1:jcplx,NU(NY,NX):NL(NY,NX),NY,NX)=0.0_r8
  TDFOMP(1:jcplx,NU(NY,NX):NL(NY,NX),NY,NX)=0.0_r8
  ROXSK(1:NPH,NU(NY,NX):NL(NY,NX),NY,NX)=0.0_r8
  end subroutine SetArrays4PlantSoilTransfer
!------------------------------------------------------------------------------------------

  subroutine GetOutput4WaterTableDepth(NY,NX,THETPZ)
  implicit none
  integer, intent(in) :: NY,NX

  real(r8), intent(in) :: THETPZ(JZ,JY,JX)
  real(r8) :: PSIS1
  real(r8) :: THETW1
  real(r8) :: THETWM
  real(r8) :: THETPX
  integer :: LL,L
  integer :: IFLGY
!     begin_execution

  IFLGY=0
  DO L=NUI(NY,NX),NLI(NY,NX)
!     IDTBL=water table flag from site file
!     THETPZ,THETPW=current,minimum air-filled, porosity for water table
!     DPTH,DTBLX=depth of soil layer midpoint, water table
!     PSIS1=water potential in hydraulic equilibrium with layer below
!     THETW1,THETWP=water content at PSIS1,minimum SWC for water table
!     DPTHT=water table depth
!
  IF(IDTBL(NY,NX).NE.0)THEN
    IF(IFLGY.EQ.0)THEN
      IF(THETPZ(L,NY,NX).LT.THETPW.OR.L.EQ.NL(NY,NX))THEN
        IFLGY=1
        IF(DPTH(L,NY,NX).LT.DTBLX(NY,NX))THEN
          D5705: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
            IF(THETPZ(LL,NY,NX).GE.THETPW.AND.LL.NE.NL(NY,NX))THEN
              IFLGY=0
              exit
            ELSEIF(DPTH(LL,NY,NX).GE.DTBLX(NY,NX))THEN
              exit
            ENDIF
          ENDDO D5705
        ENDIF
        IF(IFLGY.EQ.1)THEN
          IF(THETPZ(L,NY,NX).GE.THETPW.AND.L.NE.NL(NY,NX))THEN
            PSIS1=PSISM(L+1,NY,NX)-0.0098*(DPTH(L+1,NY,NX)-DPTH(L,NY,NX))
            THETWM=THETWP*POROS(L,NY,NX)
            THETW1=AMIN1(THETWM,EXP((PSIMS(NY,NX)-LOG(-PSIS1)) &
              *PSD(L,NY,NX)/PSISD(NY,NX)+PSL(L,NY,NX)))
            IF(THETWM.GT.THETW1)THEN
              THETPX=AMIN1(1.0,AZMAX1((THETWM-THETW(L,NY,NX))/(THETWM-THETW1)))
              DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)*(1.0-THETPX)
            ELSE
              DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
            ENDIF
          ELSEIF(L.GT.NU(NY,NX))THEN
            PSIS1=PSISM(L,NY,NX)-0.0098_r8*(DPTH(L,NY,NX)-DPTH(L-1,NY,NX))
            THETWM=THETWP*POROS(L-1,NY,NX)
            THETW1=AMIN1(THETWM,EXP((PSIMS(NY,NX)-LOG(-PSIS1)) &
              *PSD(L-1,NY,NX)/PSISD(NY,NX)+PSL(L-1,NY,NX)))
            IF(THETWM.GT.THETW1)THEN
              THETPX=AMIN1(1.0,AZMAX1((THETWM-THETW(L-1,NY,NX))/(THETWM-THETW1)))
              DPTHT(NY,NX)=CDPTH(L-1,NY,NX)-DLYR(3,L-1,NY,NX)*(1.0-THETPX)
            ELSE
              DPTHT(NY,NX)=CDPTH(L-1,NY,NX)-DLYR(3,L-1,NY,NX)
            ENDIF
          ELSE
            DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  ENDDO
  end subroutine GetOutput4WaterTableDepth
!------------------------------------------------------------------------------------------

  subroutine SetSurfaceProperty4SedErosion(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: BKVLNX
  real(r8) :: CORGM
  real(r8) :: COHS
  real(r8) :: D50
  real(r8) :: VISCWL
  REAL(R8) :: ZD50
  !     begin_execution
  !
  !     DETS=soil detachability from rainfall impact
  !     D50=average particle size
  !     CER,XER=parameters for runoff transport capacity
  !     ZD50=particle size effect on surface roughness
  !     VLS=hourly sinking rate
  !     COHS=soil cohesion
  !     DETE=soil detachability
  !     ZM=surface roughness used in runoff velocity calculation in watsub.f
  !
  BKVLNU(NY,NX)=AZMAX1(BKVLNM(NY,NX)+MWC2Soil*ORGC(NU(NY,NX),NY,NX))
  BKVLNX=SAND(NU(NY,NX),NY,NX)+SILT(NU(NY,NX),NY,NX) &
    +CLAY(NU(NY,NX),NY,NX)+1.82E-06*ORGC(NU(NY,NX),NY,NX)
  IF(BKVLNX.GT.ZEROS(NY,NX))THEN
    CORGM=MWC2Soil*ORGC(NU(NY,NX),NY,NX)/BKVLNX
    CORGC(NU(NY,NX),NY,NX)=0.55E+06_r8*CORGM
    CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)/BKVLNX
    CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)/BKVLNX
    CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)/BKVLNX
  ELSE
    CORGM=0.0_r8
    CORGC(NU(NY,NX),NY,NX)=0.0_r8
    CSAND(NU(NY,NX),NY,NX)=0.0_r8
    CSILT(NU(NY,NX),NY,NX)=1.0
    CCLAY(NU(NY,NX),NY,NX)=0.0_r8
  ENDIF
  IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
    D50=1.0_r8*CCLAY(NU(NY,NX),NY,NX)+10.0_r8*CSILT(NU(NY,NX),NY,NX) &
      +100.0_r8*CSAND(NU(NY,NX),NY,NX)+100.0_r8*CORGM
    ZD50=0.041*(ppmc*D50)**0.167_r8
    ZM(NY,NX)=ZS(NY,NX)+ZD50+1.0_r8*VOLR(NY,NX)/AREA(3,0,NY,NX)
    CER(NY,NX)=((D50+5.0_r8)/0.32_r8)**(-0.6_r8)
    XER(NY,NX)=((D50+5.0_r8)/300.0_r8)**0.25_r8
    DETS(NY,NX)=ppmc*(1.0_r8+2.0_r8*(1.0_r8-CSILT(NU(NY,NX),NY,NX)-CORGM))
    COHS=2.0_r8+10.0_r8*(CCLAY(NU(NY,NX),NY,NX)+CORGM) &
      +5.0_r8*(1.0_r8-EXP(-2.0E-06_r8*RTDNT(NU(NY,NX),NY,NX)))
    DETE(NY,NX)=0.79_r8*EXP(-0.85_r8*AMAX1(1.0_r8,COHS))
    PTDSNU(NY,NX)=1.30_r8*CORGM+2.66_r8*(1.0_r8-CORGM)
    VISCWL=VISCW*EXP(0.533_r8-0.0267_r8*TCS(0,NY,NX))
    VLS(NY,NX)=3.6E+03_r8*9.8_r8*(PTDSNU(NY,NX)-1.0_r8) &
      *(ppmc*D50)**2/(18.0_r8*VISCWL)
  ENDIF
  end subroutine SetSurfaceProperty4SedErosion
!------------------------------------------------------------------------------------------

  subroutine UpdateTotalSOC(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: OC
  integer :: K,M,N,NGL,L
  !     begin_execution
  !
  !     TOTAL SOC FOR CALCULATING CHANGES IN SOC CALCULATED IN NITRO.F
  !
  !     OMC=microbial biomass, ORC=microbial residue
  !     OQC,OQCH=DOC in micropores,macropores
  !     OQA,OQAH=acetate in micropores,macropores
  !     OHC,OHA=adsorbed SOC,acetate
  !     OSC=SOC(K=0:woody litter, K=1:non-woody litter,
  !     K=2:manure, K=3:POC, K=4:humus)
  !
  OC=0.0_r8

  DO L=0,NL(NY,NX)
!  add heterotrophic complexs
    OC=OC+sum(OMC(1:nlbiomcp,1:NMICBSO,1:jcplx,L,NY,NX))

!  add autotrophic complex
    OC=OC+sum(OMCff(1:nlbiomcp,1:NMICBSA,L,NY,NX))
!  add microbial residue
    OC=OC+SUM(ORC(1:ndbiomcp,1:jcplx,L,NY,NX))
!  add dissolved/sorbed OM and acetate
    OC=OC+SUM(OQC(1:jcplx,L,NY,NX))+SUM(OQCH(1:jcplx,L,NY,NX)) &
         +SUM(OHC(1:jcplx,L,NY,NX))+SUM(OQA(1:jcplx,L,NY,NX)) &
         +SUM(OQAH(1:jcplx,L,NY,NX))+SUM(OHA(1:jcplx,L,NY,NX))
!  add OM complexes
    OC=OC+SUM(OSC(1:jsken,1:jcplx,L,NY,NX))
!
    ORGCX(L,NY,NX)=OC
  ENDDO
  end subroutine UpdateTotalSOC
!------------------------------------------------------------------------------------------

  subroutine DiagActiveLayerDepth(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: VOLIT,VOLAT,VOLITL,VOLWTL,VOLATL
  integer :: LL,L
  integer :: ICHKA
  logical :: goto5701

!     begin_execution
  ICHKA = 0
  DO L=NUI(NY,NX),NLI(NY,NX)
!
!     VOLI,VOLIH=ice volume in micropores,macropores
!     VOLW,VOLWH=water volume in micropores,macropores
!     VOLA,VOLAH=total volume in micropores,macropores
!     DPTHA=active layer depth
!     CDPTH,DLYR=depth to bottom,thickness of soil layer
!
  IF(ICHKA.EQ.0)THEN
    VOLIT=VOLI(L,NY,NX)+VOLIH(L,NY,NX)
    VOLAT=VOLA(L,NY,NX)+VOLAH(L,NY,NX)
    IF(VOLAT.GT.ZEROS2(NY,NX).AND.VOLIT.GT.0.01*VOLAT)THEN
      D5700: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
        VOLITL=VOLI(LL,NY,NX)+VOLIH(LL,NY,NX)
        VOLWTL=VOLW(LL,NY,NX)+VOLWH(LL,NY,NX)
        VOLATL=VOLA(LL,NY,NX)+VOLAH(LL,NY,NX)
        goto5701=(VOLATL.GT.ZEROS2(NY,NX).AND.VOLITL.LT.0.01_r8*VOLATL)
        if(goto5701)exit
      ENDDO D5700
      if(.not. goto5701)then
        IF(VOLAT.GT.ZEROS2(NY,NX))THEN
          DPTHA(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)*AMIN1(1.0_r8,VOLIT/VOLAT)
        ELSE
          DPTHA(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
        ENDIF
        ICHKA=1
      else
        DPTHA(NY,NX)=9999.0_r8
      endif
    ENDIF
  ENDIF
  ENDDO
  end subroutine DiagActiveLayerDepth
!------------------------------------------------------------------------------------------

  subroutine SetTracerPropertyInLiterAir(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: TFACL
  real(r8) :: TFACR
  real(r8) :: TFACA
  real(r8) :: TFACW
  integer :: K,L,NTG
!     begin_execution
!
!     LITTER GAS CONCENTRATIOS
!
!     C*G=soil gas gaseous concentration
!     *E=atmospheric concentration
!     TKS,TCS=litter temperature (K,C)
!     S*L=gas solubility
!     C*S=soil gas aqueous concentration
!
  trc_gascl(idg_CO2,0,NY,NX)=CO2E(NY,NX)*5.36E-04_r8*TREF/TKS(0,NY,NX)
  trc_gascl(idg_CH4,0,NY,NX)=CH4E(NY,NX)*5.36E-04_r8*TREF/TKS(0,NY,NX)
  trc_gascl(idg_O2,0,NY,NX)=OXYE(NY,NX)*1.43E-03_r8*TREF/TKS(0,NY,NX)
  trc_gascl(idg_N2,0,NY,NX)=Z2GE(NY,NX)*1.25E-03_r8*TREF/TKS(0,NY,NX)
  trc_gascl(idg_N2O,0,NY,NX)=Z2OE(NY,NX)*1.25E-03_r8*TREF/TKS(0,NY,NX)
  trc_gascl(idg_NH3,0,NY,NX)=ZNH3E(NY,NX)*6.25E-04_r8*TREF/TKS(0,NY,NX)
  trc_gascl(idg_H2,0,NY,NX)=H2GE(NY,NX)*8.92E-05*TREF/TKS(0,NY,NX)

! initialize all band nutrients to zero
  trc_solcl(ids_nutb_beg:ids_nutb_end,0,NY,NX)=0.0_r8
  IF(VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
! exclude NH3B,
    DO NTG=idg_beg,idg_end-1
      trc_solcl(NTG,0,NY,NX)=AZMAX1(trc_solml(NTG,0,NY,NX)/VOLW(0,NY,NX))
    ENDDO
  ELSE
    trc_solcl(idg_beg:idg_end-1,0,NY,NX)=0.0_r8
  ENDIF
!
!     TFACL=temperature effect on diffusivity
!     *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
!     *SG PARAMETER statement above
!
  TFACL=TEFAQUDIF(TKS(0,NY,NX))
  TFND(0,NY,NX)=TFACL

  SolDifc(idg_CO2,0,NY,NX)=CLSG*TFACL
  SolDifc(idg_CH4,0,NY,NX)=CQSG*TFACL
  SolDifc(idg_O2,0,NY,NX)=OLSG*TFACL
  SolDifc(idg_N2,0,NY,NX)=ZLSG*TFACL
  SolDifc(idg_NH3,0,NY,NX)=ZNSG*TFACL
  SolDifc(idg_H2,0,NY,NX)=HLSG*TFACL
  SolDifc(idg_N2O,0,NY,NX)=ZVSG*TFACL

  SolDifc(ids_NO3,0,NY,NX)=ZOSG*TFACL
  SolDifc(ids_H1PO4,0,NY,NX)=POSG*TFACL
  SolDifc(ids_NH4,0,NY,NX)   =SolDifc(idg_NH3,0,NY,NX)
  SolDifc(ids_NH4B,0,NY,NX)  =SolDifc(ids_NH4,0,NY,NX)
  SolDifc(idg_NH3B,0,NY,NX)  =SolDifc(idg_NH3,0,NY,NX)
  SolDifc(ids_NO3B,0,NY,NX)  =SolDifc(ids_NO3,0,NY,NX)
  SolDifc(ids_NO2,0,NY,NX)   =SolDifc(ids_NO3,0,NY,NX)
  SolDifc(ids_NO2B,0,NY,NX)  =SolDifc(ids_NO2,0,NY,NX)
  SolDifc(ids_H2PO4,0,NY,NX) =SolDifc(ids_H1PO4,0,NY,NX)
  SolDifc(ids_H1PO4B,0,NY,NX)=SolDifc(ids_H1PO4,0,NY,NX)
  SolDifc(ids_H2PO4B,0,NY,NX)=SolDifc(ids_H2PO4,0,NY,NX)

  OCSGL(0,NY,NX)=OCSG*TFACL
  ONSGL(0,NY,NX)=ONSG*TFACL
  OPSGL(0,NY,NX)=OPSG*TFACL
  OASGL(0,NY,NX)=OASG*TFACL
!
!     R*Y,R*X=total substrate uptake from previous,current hour
!     used in nitro.f, uptake.f
!
  ROXYY(0,NY,NX)=ROXYX(0,NY,NX)
  RNH4Y(0,NY,NX)=RNH4X(0,NY,NX)
  RNO3Y(0,NY,NX)=RNO3X(0,NY,NX)
  RNO2Y(0,NY,NX)=RNO2X(0,NY,NX)
  RN2OY(0,NY,NX)=RN2OX(0,NY,NX)
  RP14Y(0,NY,NX)=RP14X(0,NY,NX)
  RPO4Y(0,NY,NX)=RPO4X(0,NY,NX)
  ROXYX(0,NY,NX)=0.0_r8
  RNH4X(0,NY,NX)=0.0_r8
  RNO3X(0,NY,NX)=0.0_r8
  RNO2X(0,NY,NX)=0.0_r8
  RN2OX(0,NY,NX)=0.0_r8
  RP14X(0,NY,NX)=0.0_r8
  RPO4X(0,NY,NX)=0.0_r8
  D5055: DO K=1,jcplx
    ROQCY(K,0,NY,NX)=ROQCX(K,0,NY,NX)
    ROQAY(K,0,NY,NX)=ROQAX(K,0,NY,NX)
    ROQCX(K,0,NY,NX)=0.0_r8
    ROQAX(K,0,NY,NX)=0.0_r8
  ENDDO D5055
!
!     WGSGA,WGSGR,WGSGW=vapor diffusivity in air,litter,snowpack
!
  TFACA=TEFGASDIF(TKA(NY,NX))
  WGSGA(NY,NX)=WGSG*TFACA
  TFACR=TEFGASDIF(TKS(0,NY,NX))
  WGSGR(NY,NX)=WGSG*TFACR
  D5060: DO  L=1,JS
    TFACW=TEFGASDIF(TKW(L,NY,NX))
    WGSGW(L,NY,NX)=WGSG*TFACW
  ENDDO D5060
  end subroutine SetTracerPropertyInLiterAir
!------------------------------------------------------------------------------------------

  subroutine GetSoluteConcentrations(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
!     begin_execution
!     CALCULATE SOIL CONCENTRATIONS OF NH4, NH3, NO3, PO4
!     IN BAND AND NON-BAND ZONES
!
!     C*S=solute concentration in non-band
!     CH1P4,CH2P4=HPO4,H2PO4 concentration in non-band
!     Z*P=P ion pair amounts in non-band (see solute.f)
!     VLNH4,VLNO3,VLPO4=fraction of soil volume in NH4,NO3,PO4 non-band
!
  DO L=NUI(NY,NX),NLI(NY,NX)
    IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN

      IF(trcs_VLN(ids_NH4,L,NY,NX).GT.ZERO)THEN
        trc_solcl(ids_NH4,L,NY,NX)=AZMAX1(trc_solml(ids_NH4,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)))
        trc_solcl(idg_NH3,L,NY,NX)=AZMAX1(trc_solml(idg_NH3,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(idg_NH3,L,NY,NX)))
      ELSE
        trc_solcl(ids_NH4,L,NY,NX)=0.0_r8
        trc_solcl(idg_NH3,L,NY,NX)=0.0_r8
      ENDIF
      IF(trcs_VLN(ids_NO3,L,NY,NX).GT.ZERO)THEN
        trc_solcl(ids_NO3,L,NY,NX)=AZMAX1(trc_solml(ids_NO3,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)))
        trc_solcl(ids_NO2,L,NY,NX)=AZMAX1(trc_solml(ids_NO2,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_NO2,L,NY,NX)))
      ELSE
        trc_solcl(ids_NO3,L,NY,NX)=0.0_r8
        trc_solcl(ids_NO2,L,NY,NX)=0.0_r8
      ENDIF

      IF(trcs_VLN(ids_H1PO4,L,NY,NX).GT.ZERO)THEN
        trc_solcl(ids_H1PO4,L,NY,NX)=AZMAX1(trc_solml(ids_H1PO4,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)))
        trc_solcl(ids_H2PO4,L,NY,NX)=AZMAX1(trc_solml(ids_H2PO4,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_H2PO4,L,NY,NX)))

        CPO4S(L,NY,NX)=AZMAX1(((trcsa_solml(idsa_H0PO4,L,NY,NX)+trcsa_solml(idsa_H3PO4,L,NY,NX) &
          +trcsa_solml(idsa_FeHPO4,L,NY,NX)+trcsa_solml(idsa_FeH2PO4,L,NY,NX)+trcsa_solml(idsa_CaPO4,L,NY,NX) &
          +trcsa_solml(idsa_CaHPO4,L,NY,NX)+trcsa_solml(idsa_CaH2PO4,L,NY,NX)+trcsa_solml(idsa_MgHPO4,L,NY,NX))*patomw &
          +trc_solml(ids_H1PO4,L,NY,NX)+trc_solml(ids_H2PO4,L,NY,NX))/(VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)))
      ELSE
        trc_solcl(ids_H1PO4,L,NY,NX)=0.0_r8
        trc_solcl(ids_H2PO4,L,NY,NX)=0.0_r8
        CPO4S(L,NY,NX)=0.0_r8
      ENDIF
!
!     C*B=solute concentration in band
!     CH1PB,CH2PB=HPO4,H2PO4 concentration in band
!     Z*B=P ion pair amounts in band (see solute.f)
!     VLNHB,VLNOB,VLPOB=fraction of soil volume in NH4,NO3,PO4 band
!
      IF(trcs_VLN(ids_NH4B,L,NY,NX).GT.ZERO)THEN
        trc_solcl(ids_NH4B,L,NY,NX)=AZMAX1(trc_solml(ids_NH4B,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)))
        trc_solcl(idg_NH3B,L,NY,NX)=AZMAX1(trc_solml(idg_NH3B,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(idg_NH3B,L,NY,NX)))
      ELSE
        trc_solcl(ids_NH4B,L,NY,NX)=0.0_r8
        trc_solcl(idg_NH3B,L,NY,NX)=0.0_r8
      ENDIF

      IF(trcs_VLN(ids_NO3B,L,NY,NX).GT.ZERO)THEN
        trc_solcl(ids_NO3B,L,NY,NX)=AZMAX1(trc_solml(ids_NO3B,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)))
        trc_solcl(ids_NO2B,L,NY,NX)=AZMAX1(trc_solml(ids_NO2B,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_NO2B,L,NY,NX)))
      ELSE
        trc_solcl(ids_NO3B,L,NY,NX)=0.0_r8
        trc_solcl(ids_NO2B,L,NY,NX)=0.0_r8
      ENDIF

      IF(trcs_VLN(ids_H1PO4B,L,NY,NX).GT.ZERO)THEN
        trc_solcl(ids_H1PO4B,L,NY,NX)=AZMAX1(trc_solml(ids_H1PO4B,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)))
        trc_solcl(ids_H2PO4B,L,NY,NX)=AZMAX1(trc_solml(ids_H2PO4B,L,NY,NX)/(VOLW(L,NY,NX)*trcs_VLN(ids_H2PO4B,L,NY,NX)))

        CPO4B(L,NY,NX)=AZMAX1(((trcsa_solml(idsa_H0PO4B,L,NY,NX)+trcsa_solml(idsa_H3PO4B,L,NY,NX) &
          +trcsa_solml(idsa_FeHPO4B,L,NY,NX)+trcsa_solml(idsa_FeH2PO4B,L,NY,NX)+trcsa_solml(idsa_CaPO4B,L,NY,NX) &
          +trcsa_solml(idsa_CaHPO4B,L,NY,NX)+trcsa_solml(idsa_CaH2PO4B,L,NY,NX)+trcsa_solml(idsa_MgHPO4B,L,NY,NX))*patomw &
          +trc_solml(ids_H1PO4B,L,NY,NX)+trc_solml(ids_H2PO4B,L,NY,NX))/(VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)))
      ELSE
        trc_solcl(ids_H1PO4B,L,NY,NX)=0.0_r8
        trc_solcl(ids_H2PO4B,L,NY,NX)=0.0_r8
        CPO4B(L,NY,NX)=0.0_r8
      ENDIF

    ELSE
      trc_solcl(ids_nuts_beg:ids_nuts_end,L,NY,NX)=0.0_r8
      CPO4S(L,NY,NX)=0.0_r8
      CPO4B(L,NY,NX)=0.0_r8
    ENDIF
  ENDDO
  end subroutine GetSoluteConcentrations
!------------------------------------------------------------------------------------------

  subroutine Prep4PlantMicrobeUptake(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: TFACL,TFACG
  real(r8) :: FH2O
  real(r8) :: ZC3,ZA3,ZC2,ZA2,ZC1,ZA1
  REAL(R8) :: ZN
  real(r8) :: ZION1
  integer :: K,L

! begin_execution
  DO L=NUI(NY,NX),NLI(NY,NX)
!
! PREPARE ARRAYS FOR TOTAL O2 UPTAKE AND NH4,NO3.NO2,N2O,HPO4,H2PO4
! UPTAKE IN NON-BAND,BAND AND DOC,DON,DOP,ACETATE UPTAKE
!
! R*Y,R*X=total substrate uptake from previous,current hour
! used in nitro.f, uptake.f
!
    ROXYY(L,NY,NX)=ROXYX(L,NY,NX)
    RNH4Y(L,NY,NX)=RNH4X(L,NY,NX)
    RNO3Y(L,NY,NX)=RNO3X(L,NY,NX)
    RNO2Y(L,NY,NX)=RNO2X(L,NY,NX)
    RN2OY(L,NY,NX)=RN2OX(L,NY,NX)
    RP14Y(L,NY,NX)=RP14X(L,NY,NX)
    RPO4Y(L,NY,NX)=RPO4X(L,NY,NX)
    RNHBY(L,NY,NX)=RNHBX(L,NY,NX)
    RN3BY(L,NY,NX)=RN3BX(L,NY,NX)
    RN2BY(L,NY,NX)=RN2BX(L,NY,NX)
    RP1BY(L,NY,NX)=RP1BX(L,NY,NX)
    RPOBY(L,NY,NX)=RPOBX(L,NY,NX)
    ROXYX(L,NY,NX)=0.0_r8
    RNH4X(L,NY,NX)=0.0_r8
    RNO3X(L,NY,NX)=0.0_r8
    RNO2X(L,NY,NX)=0.0_r8
    RN2OX(L,NY,NX)=0.0_r8
    RP14X(L,NY,NX)=0.0_r8
    RPO4X(L,NY,NX)=0.0_r8
    RNHBX(L,NY,NX)=0.0_r8
    RN3BX(L,NY,NX)=0.0_r8
    RN2BX(L,NY,NX)=0.0_r8
    RP1BX(L,NY,NX)=0.0_r8
    RPOBX(L,NY,NX)=0.0_r8

    D5050: DO K=1,jcplx
      ROQCY(K,L,NY,NX)=ROQCX(K,L,NY,NX)
      ROQAY(K,L,NY,NX)=ROQAX(K,L,NY,NX)
      ROQCX(K,L,NY,NX)=0.0_r8
      ROQAX(K,L,NY,NX)=0.0_r8
    ENDDO D5050
!
! DIFFUSIVITY
!
! TFACG,TFACL=temperature effects on gaseous,aqueous diffusivity
!
! *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
! *SG PARAMETER statement above
!
    TFACG=TEFGASDIF(TKS(L,NY,NX))
    TFACL=TEFAQUDIF(TKS(L,NY,NX))
    TFND(L,NY,NX)=TFACL

    GasDifc(idg_CO2,L,NY,NX)=CGSG*TFACG
    GasDifc(idg_CH4,L,NY,NX)=CHSG*TFACG
    GasDifc(idg_O2,L,NY,NX)=OGSG*TFACG
    GasDifc(idg_N2,L,NY,NX)=ZGSG*TFACG
    GasDifc(idg_N2O,L,NY,NX)=Z2SG*TFACG
    GasDifc(idg_NH3,L,NY,NX)=ZHSG*TFACG
    GasDifc(idg_H2,L,NY,NX)=HGSG*TFACG
    GasDifc(idg_NH3B,L,NY,NX)=ZHSG*TFACG

    SolDifc(idg_CO2,L,NY,NX)=CLSG*TFACL
    SolDifc(idg_CH4,L,NY,NX)=CQSG*TFACL
    SolDifc(idg_O2,L,NY,NX)=OLSG*TFACL
    SolDifc(idg_N2,L,NY,NX)=ZLSG*TFACL
    SolDifc(idg_NH3,L,NY,NX)=ZNSG*TFACL
    SolDifc(idg_H2,L,NY,NX)=HLSG*TFACL
    SolDifc(idg_N2O,L,NY,NX)=ZVSG*TFACL
    SolDifc(idg_NH3B,L,NY,NX)=SolDifc(idg_NH3,L,NY,NX)
    SolDifc(ids_NO3,L,NY,NX)=ZOSG*TFACL
    SolDifc(ids_H1PO4,L,NY,NX)=POSG*TFACL

    SolDifc(ids_NH4,L,NY,NX)   =SolDifc(idg_NH3,L,NY,NX)
    SolDifc(ids_NH4B,L,NY,NX)  =SolDifc(ids_NH4,L,NY,NX)
    SolDifc(ids_NO3B,L,NY,NX)  =SolDifc(ids_NO3,L,NY,NX)
    SolDifc(ids_NO2,L,NY,NX)   =SolDifc(ids_NO3,L,NY,NX)
    SolDifc(ids_NO2B,L,NY,NX)  =SolDifc(ids_NO2,L,NY,NX)
    SolDifc(ids_H2PO4,L,NY,NX) =SolDifc(ids_H1PO4,L,NY,NX)
    SolDifc(ids_H1PO4B,L,NY,NX)=SolDifc(ids_H1PO4,L,NY,NX)
    SolDifc(ids_H2PO4B,L,NY,NX)=SolDifc(ids_H2PO4,L,NY,NX)

    OCSGL(L,NY,NX)=OCSG*TFACL
    ONSGL(L,NY,NX)=ONSG*TFACL
    OPSGL(L,NY,NX)=OPSG*TFACL
    OASGL(L,NY,NX)=OASG*TFACL
    WGSGL(L,NY,NX)=WGSG*TFACG

    IF(ISALTG.NE.0)THEN
      ALSGL(L,NY,NX)=ALSG*TFACL
      FESGL(L,NY,NX)=FESG*TFACL
      HYSGL(L,NY,NX)=HYSG*TFACL
      CASGL(L,NY,NX)=CASG*TFACL
      GMSGL(L,NY,NX)=GMSG*TFACL
      ANSGL(L,NY,NX)=ANSG*TFACL
      AKSGL(L,NY,NX)=AKSG*TFACL
      OHSGL(L,NY,NX)=OHSG*TFACL
      C3SGL(L,NY,NX)=C3SG*TFACL
      HCSGL(L,NY,NX)=HCSG*TFACL
      SOSGL(L,NY,NX)=SOSG*TFACL
      CLSXL(L,NY,NX)=CLSX*TFACL
!
!   TOTAL ION CONCENTRATION
!
!   ZC3,ZA3,ZC2,ZA2,ZC1,ZA1=total tri-,di-,univalent cations C,anions A
!   CSTR,CION=ion strength, total ion concentration
!
      ZC3=trcsa_solml(idsa_Al,L,NY,NX)+trcsa_solml(idsa_Fe,L,NY,NX)
      ZA3=trcsa_solml(idsa_H0PO4,L,NY,NX)+trcsa_solml(idsa_H0PO4B,L,NY,NX)
      ZC2=trcsa_solml(idsa_Ca,L,NY,NX)+trcsa_solml(idsa_Mg,L,NY,NX) &
        +trcsa_solml(idsa_AlOH,L,NY,NX)+trcsa_solml(idsa_FeOH,L,NY,NX) &
        +trcsa_solml(idsa_FeH2PO4,L,NY,NX)+trcsa_solml(idsa_FeH2PO4B,L,NY,NX)
      ZA2=trcsa_solml(idsa_SO4,L,NY,NX)+trcsa_solml(idsa_CO3,L,NY,NX)+trc_solml(ids_H1PO4,L,NY,NX)+trc_solml(ids_H1PO4B,L,NY,NX)

      ZC1=(trc_solml(ids_NH4,L,NY,NX)+trc_solml(ids_NH4B,L,NY,NX))/natomw+trcsa_solml(idsa_Hp,L,NY,NX) &
        +trcsa_solml(idsa_Na,L,NY,NX)+trcsa_solml(idsa_K,L,NY,NX) &
        +trcsa_solml(idsa_AlOH2,L,NY,NX)+trcsa_solml(idsa_FeOH2,L,NY,NX) &
        +trcsa_solml(idsa_AlSO4,L,NY,NX)+trcsa_solml(idsa_FeSO4,L,NY,NX) &
        +trcsa_solml(idsa_CaOH2,L,NY,NX)+trcsa_solml(idsa_CaHCO3,L,NY,NX) &
        +trcsa_solml(idsa_MgOH2,L,NY,NX)+trcsa_solml(idsa_MgHCO3,L,NY,NX)&
        +trcsa_solml(idsa_FeHPO4,L,NY,NX)+trcsa_solml(idsa_FeHPO4B,L,NY,NX) &
        +trcsa_solml(idsa_CaH2PO4,L,NY,NX)+trcsa_solml(idsa_CaH2PO4B,L,NY,NX)

      ZA1=(trc_solml(ids_NO3,L,NY,NX)+trc_solml(ids_NO3B,L,NY,NX))/natomw+trcsa_solml(idsa_OH,L,NY,NX) &
        +trcsa_solml(idsa_HCO3,L,NY,NX)+trcsa_solml(idsa_Cl,L,NY,NX) &
        +trcsa_solml(idsa_AlOH4,L,NY,NX)+trcsa_solml(idsa_FeOH4,L,NY,NX) &
        +trcsa_solml(idsa_NaCO3,L,NY,NX)+trcsa_solml(idsa_NaSO4,L,NY,NX) &
        +trcsa_solml(idsa_KSO4,L,NY,NX)+(trc_solml(ids_H2PO4,L,NY,NX) &
        +trc_solml(ids_H2PO4B,L,NY,NX))/patomw+trcsa_solml(idsa_CaPO4,L,NY,NX) &
        +trcsa_solml(idsa_CaPO4B,L,NY,NX)

      ZN=trc_solml(idg_CO2,L,NY,NX)/catomw+trc_solml(idg_CH4,L,NY,NX)/catomw+trc_solml(idg_O2,L,NY,NX)/32.0 &
        +(trc_solml(idg_N2,L,NY,NX)+trc_solml(idg_N2O,L,NY,NX)+trc_solml(idg_NH3,L,NY,NX)+trc_solml(idg_NH3B,L,NY,NX))/natomw &
        +trcsa_solml(idsa_AlOH3,L,NY,NX)+trcsa_solml(idsa_FeOH3,L,NY,NX)+trcsa_solml(idsa_CaCO3,L,NY,NX)+trcsa_solml(idsa_CaSO4,L,NY,NX) &
        +trcsa_solml(idsa_MgCO3,L,NY,NX)+trcsa_solml(idsa_MgSO4,L,NY,NX)+trcsa_solml(idsa_H3PO4,L,NY,NX)+trcsa_solml(idsa_CaHPO4,L,NY,NX) &
        +trcsa_solml(idsa_MgHPO4,L,NY,NX)+trcsa_solml(idsa_H3PO4B,L,NY,NX)+trcsa_solml(idsa_CaHPO4B,L,NY,NX)+trcsa_solml(idsa_MgHPO4B,L,NY,NX)

      ZION1=ABS(3.0_r8*(ZC3-ZA3)+2.0_r8*(ZC2-ZA2)+ZC1-ZA1)
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        CSTR(L,NY,NX)=AZMAX1(0.5E-03_r8*(9.0_r8*(ZC3+ZA3)+4.0_r8*(ZC2+ZA2) &
          +ZC1+ZA1+ZION1)/VOLW(L,NY,NX))
        CION(L,NY,NX)=AZMAX1((ZC3+ZA3+ZC2+ZA2+ZC1+ZA1+ZN)/VOLW(L,NY,NX))
      ELSE
        CSTR(L,NY,NX)=0.0_r8
        CION(L,NY,NX)=0.0_r8
      ENDIF
    ENDIF
!
! OSTWALD COEFFICIENTS FOR CO2, CH4, O2, N2, N2O, NH3 AND H2
! SOLUBILITY IN WATER
!

  ENDDO
  end subroutine Prep4PlantMicrobeUptake
!------------------------------------------------------------------------------------------

  subroutine GetSurfResidualProperties(NY,NX,DPTH0)

  implicit none
  integer, intent(in) :: NY,NX
  real(r8),intent(out) :: DPTH0(JY,JX)
  real(r8) :: TVOLG0,TVOLWI,THETWR
  real(r8) :: VOLWRZ
  real(r8) :: VOLIRZ
  real(r8) :: XVOLW0
  real(r8) :: XVOLI0
  integer  :: NTG,NTN
! begin_execution
! PHYSICAL PROPERTIES, AND WATER, GAS, AND MINERAL CONTENTS
! OF SURFACE RESIDUE
!
! VOLWRX=liter water holding capacity
! THETRX,RC0=specific WHC,mass of woody(0),fine(1),manure(2) litter
! VOLR=dry litter volume
! BKRS=dry bulk density of woody(0),fine(1),manure(2) litter
! TVOLG0=excess litter water+ice
! VOLT,VOLX=wet litter volume
! BKVL=litter mass
! VOLW,VOLI,VOLA,VOLP=litter water,ice,porosity,air volume
! THETW,THETI,THETA,THETP=litter water,ice,porosity,air concentration
! POROS=litter porosity
! THETW0,THETI0,DPTH0=litter excess water,ice,water+ice depth
! DLYR=litter thickness
! PSISM,PSISE=litter matric,saturation water potential
!
  TVOLG0=AZMAX1(VOLW(0,NY,NX)+VOLI(0,NY,NX)-VOLWRX(NY,NX))
  VOLT(0,NY,NX)=TVOLG0+VOLR(NY,NX)
  IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VOLX(0,NY,NX)=VOLT(0,NY,NX)
    BKVL(0,NY,NX)=MWC2Soil*ORGC(0,NY,NX)
    VOLA(0,NY,NX)=AZMAX1(VOLR(NY,NX)-BKVL(0,NY,NX)/1.30_r8)
    VOLP(0,NY,NX)=AZMAX1(VOLA(0,NY,NX)-VOLW(0,NY,NX)-VOLI(0,NY,NX))
    IF(VOLR(NY,NX).GT.ZEROS(NY,NX))THEN
      POROS(0,NY,NX)=VOLA(0,NY,NX)/VOLR(NY,NX)
      THETW(0,NY,NX)=AZMAX1(AMIN1(1.0_r8,VOLW(0,NY,NX)/VOLR(NY,NX)))
      THETI(0,NY,NX)=AZMAX1(AMIN1(1.0_r8,VOLI(0,NY,NX)/VOLR(NY,NX)))
      THETP(0,NY,NX)=AZMAX1(AMIN1(1.0_r8,VOLP(0,NY,NX)/VOLR(NY,NX)))
    ELSE
      POROS(0,NY,NX)=1.0_r8
      THETW(0,NY,NX)=0.0_r8
      THETI(0,NY,NX)=0.0_r8
      THETP(0,NY,NX)=0.0_r8
    ENDIF
    TVOLWI=VOLW(0,NY,NX)+VOLI(0,NY,NX)
    IF(TVOLWI.GT.ZEROS(NY,NX))THEN
      VOLWRZ=VOLW(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      VOLIRZ=VOLI(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      XVOLW0=AZMAX1(VOLW(0,NY,NX)-VOLWRZ)/AREA(3,NU(NY,NX),NY,NX)
      XVOLI0=AZMAX1(VOLI(0,NY,NX)-VOLIRZ)/AREA(3,NU(NY,NX),NY,NX)
    ELSE
      XVOLW0=0.0_r8
      XVOLI0=0.0_r8
    ENDIF
    DPTH0(NY,NX)=XVOLW0+XVOLI0
    DLYR(3,0,NY,NX)=VOLX(0,NY,NX)/AREA(3,0,NY,NX)
    IF(VOLR(NY,NX).GT.ZEROS(NY,NX).AND.VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW(0,NY,NX))/VOLR(NY,NX)
      IF(THETWR.LT.FC(0,NY,NX))THEN
        PSISM(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)+((FCL(0,NY,NX)-LOG(THETWR)) &
          /FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS(0,NY,NX))THEN
        PSISM(0,NY,NX)=-EXP(PSIMS(NY,NX)+(((PSL(0,NY,NX)-LOG(THETWR)) &
          /PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
        PSISM(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      PSISO(0,NY,NX)=0.0_r8
      PSISH(0,NY,NX)=0.0098_r8*(ALT(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX) &
        +0.5_r8*DLYR(3,0,NY,NX))
      PSIST(0,NY,NX)=AZMIN1(PSISM(0,NY,NX)+PSISO(0,NY,NX)+PSISH(0,NY,NX))
!
!     LITTER NH4,NH3,NO3,NO2,HPO4,H2PO4 CONCENTRATIONS
!
!     C*=litter solute concentrations
!
      DO NTN=ids_nut_beg,ids_nuts_end
        trc_solcl(NTN,0,NY,NX)=AZMAX1(trc_solml(NTN,0,NY,NX)/VOLW(0,NY,NX))
      ENDDO

    ELSE
      PSISM(0,NY,NX)=PSISM(NU(NY,NX),NY,NX)
      trc_solcl(ids_nut_beg:ids_nuts_end,0,NY,NX)=0.0_r8
    ENDIF
  ELSE
    VOLX(0,NY,NX)=0.0_r8
    BKVL(0,NY,NX)=0.0_r8
    VOLA(0,NY,NX)=0.0_r8
    VOLP(0,NY,NX)=0.0_r8
    POROS(0,NY,NX)=1.0
    DLYR(3,0,NY,NX)=0.0_r8
    THETW(0,NY,NX)=0.0_r8
    THETI(0,NY,NX)=0.0_r8
    THETP(0,NY,NX)=1.0
    VOLWRX(NY,NX)=0.0_r8
    PSISM(0,NY,NX)=PSISM(NU(NY,NX),NY,NX)

    trc_solcl(ids_nut_beg:ids_nuts_end,0,NY,NX)=0.0_r8
    trc_solcl(idg_beg:idg_end-1,0,NY,NX)=0.0_r8
  ENDIF
  end subroutine GetSurfResidualProperties

!------------------------------------------------------------------------------------------

  subroutine ApplyFertilizerAtNoon(I,J,NHW,NHE,NVN,NVS)
!
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS

  integer :: NX,NY
  real(r8) :: OFC(2),OFN(2),OFP(2)
  integer :: LFDPTH = 0
!     begin_execution

  D8990: DO NX=NHW,NHE
    D8995: DO NY=NVN,NVS
      IF(J.EQ.INT(ZNOON(NY,NX)))THEN

        call ApplyMineralFertilizer(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
!
!     SOIL LAYER NUMBER IN WHICH PLANT OR ANIMAL RESIDUES ARE APPLIED
!
        call ApplyPlantAnimalResidue(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
!
!     FERTILIZER UREA, NITRIFICATION INHIBITORS
        call ApplyUreaNitrifierInhibitor(I,J,NY,NX,LFDPTH)

      ENDIF

    ENDDO D8995
  ENDDO D8990
  end subroutine ApplyFertilizerAtNoon
!------------------------------------------------------------------------------------------

  subroutine ApplyUreaNitrifierInhibitor(I,J,NY,NX,LFDPTH)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(in) :: LFDPTH
  integer :: L
!     begin_execution
!
!     IYTYP=fertilizer release type from fertilizer input file
!     FERT=fertilizer type from fertilizer input file
!     IUTYP=urea hydrolysis inhibitor type (1=no,2=yes)
!     ZNHU0,ZNHUI=initial,current urea hydrolysis inhibition activity
!     ZNFN0,ZNFNI=initial,current nitrification inhibition activity
!
  IF(FERT(3,I,NY,NX).GT.0.0_r8.OR.FERT(7,I,NY,NX).GT.0.0_r8)THEN
    IF(IYTYP(0,I,NY,NX).EQ.0)THEN
      IUTYP(NY,NX)=0
    ELSEIF(IYTYP(0,I,NY,NX).EQ.1.OR.IYTYP(0,I,NY,NX).EQ.3)THEN
      IUTYP(NY,NX)=1
    ELSE
      !urea hydrolysis is on
      IUTYP(NY,NX)=2
    ENDIF
    D9964: DO L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN
        ZNHU0(L,NY,NX)=1.0_r8
        ZNHUI(L,NY,NX)=1.0_r8
      ELSE
        ZNHU0(L,NY,NX)=0.0_r8
        ZNHUI(L,NY,NX)=0.0_r8
      ENDIF
    ENDDO D9964
  ENDIF
  IF(IYTYP(0,I,NY,NX).EQ.3.OR.IYTYP(0,I,NY,NX).EQ.4)THEN
    D9965: DO L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN
        ZNFN0(L,NY,NX)=1.0_r8
        ZNFNI(L,NY,NX)=1.0_r8
      ELSE
        ZNFN0(L,NY,NX)=0.0_r8
        ZNFNI(L,NY,NX)=0.0_r8
      ENDIF
    ENDDO D9965
  ENDIF
  end subroutine ApplyUreaNitrifierInhibitor
!------------------------------------------------------------------------------------------

  subroutine ApplyPlantAnimalResidue(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(inout) :: LFDPTH
  real(r8), intent(in) :: OFC(2),OFN(2),OFP(2)
  real(r8) :: CNOF(4),CPOF(4)
  real(r8) :: CORGCX
  real(r8) :: CNOFT
  real(r8) :: CPOFT
  real(r8) :: FDPTHM
  real(r8) :: FRNT,FRPT
  REAL(R8) :: OSCI,OSNI,OSPI
  REAL(R8) :: OSCX,OSNX,OSPX
  REAL(R8) :: OMC1,OMN1,OMP1
  real(r8) :: OQC1,OQN1,OQP1
  real(r8) :: OSC1,OSN1,OSP1
  REAL(R8) :: RNT,RPT
  integer  :: L,K,M,N,NN,NGL
  real(r8) :: tglds
  real(r8) :: OMC1g,OMN1g,OMP1g
!     begin_execution
  associate(                           &
    k_fine_litr => micpar%k_fine_litr, &
    k_manure    => micpar%k_manure   , &
    ilignin     => micpar%ilignin    , &
    icellulos   => micpar%icellulos  , &
    icarbhyro   => micpar%icarbhyro  , &
    iprotein    => micpar%iprotein     &
  )
!     LFDPTH=layer number
!
  IF(OFC(1)+OFC(2).GT.0.0_r8)THEN
    DO  L=0,JZ
      FDPTHM=FDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      IF(FDPTHM.LE.0.0_r8)THEN
        LFDPTH=0
        exit
      ELSEIF(CDPTH(L,NY,NX).GE.FDPTHM)THEN
        LFDPTH=L
        exit
      ENDIF
    ENDDO
!
!     ALLOCATION OF PLANT RESIDUE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     CFOSC=fraction of litter allocated to protein(1)
!     soluble CH2O(2), cellulose(3) and lignin(4)
!     ITYPE=litter type entered in fertilizer input file
!
!     MAIZE
!
    IF(IYTYP(1,I,NY,NX).EQ.1)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.080_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.245_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.613_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.062_r8
!
!     WHEAT
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.2)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.125_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.171_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.560_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.144_r8
!
!     SOYBEAN
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.3)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.138_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.426_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.316_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.120_r8
!
!     OLD STRAW
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.4)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.075_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.125_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.550_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.250_r8
!
!     STRAW
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.5)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.036_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.044_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.767_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.153_r8
!
!     COMPOST
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.6)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.143_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.015_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.640_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.202_r8
!
!     GREEN MANURE
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.7)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.202_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.013_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.560_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.225_r8
!
!     SIMPLE SUBSTRATE
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.10)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.000_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=1.000_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.000_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.000_r8
    ELSE
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.075_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.125_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.550_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.250_r8
    ENDIF
!
!     ALLOCATION OF ANIMAL MANURE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     RUMINANT
!
    IF(IYTYP(2,I,NY,NX).EQ.1)THEN
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.036_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.044_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.630_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.290_r8
!
!     NON-RUMINANT
!
    ELSEIF(IYTYP(2,I,NY,NX).EQ.2)THEN
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.138_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.401_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.316_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.145_r8
!
!     GRAZING
!
    ELSEIF(IYTYP(2,I,NY,NX).EQ.3)THEN
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.036_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.044_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.630_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.290_r8
!
!     OTHER
!
    ELSE
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.138_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.401_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.316_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.145_r8
    ENDIF
!
!     DISTRIBUTE RESIDUE APPLICATION AMONG COMPONENTS OF RESIDUE COMPLEX
!
!     OFC,OFN,OFP=litter C,N,P application from fertilizer file
!

    D2965: DO K=1,2
      OSCI=OFC(K)*AREA(3,LFDPTH,NY,NX)
      OSNI=OFN(K)*AREA(3,LFDPTH,NY,NX)
      OSPI=OFP(K)*AREA(3,LFDPTH,NY,NX)
      IF(BKVL(LFDPTH,NY,NX).GT.ZEROS(NY,NX))THEN
        CORGCX=OSCI/BKVL(LFDPTH,NY,NX)
      ELSE
        CORGCX=0.55E+06_r8
      ENDIF
      OSCX=0.0_r8
      OSNX=0.0_r8
      OSPX=0.0_r8
!
!     BIOMASSES OF MICROBIAL POPULATIONS IN RESIDUE
!
!     OMC,OMN,OMP=microbial biomass in litter application
!     OMCI=microbial biomass content in litter
!     OMCF,OMCA=hetero,autotrophic biomass composition in litter
!
      D2960: DO N=1,NFGs
        tglds=JGnfo(N)-JGnfo(N)+1
        D2961: DO M=1,nlbiomcp
          OMC1=AZMAX1(AMIN1(OSCI*micpar%OMCI(M,K)*micpar%OMCF(N),OSCI-OSCX))
          OMN1=AZMAX1(AMIN1(OMC1*micpar%CNOMCa(M,N,K),OSNI-OSNX))
          OMP1=AZMAX1(AMIN1(OMC1*micpar%CPOMCa(M,N,K),OSPI-OSPX))
          DO NGL=JGnio(N),JGnfo(N)
            OMC1g=OMC1/tglds
            OMN1g=OMN1/tglds
            OMP1g=OMP1/tglds
            OMC(M,NGL,K,LFDPTH,NY,NX)=OMC(M,NGL,K,LFDPTH,NY,NX)+OMC1g
            OMN(M,NGL,K,LFDPTH,NY,NX)=OMN(M,NGL,K,LFDPTH,NY,NX)+OMN1g
            OMP(M,NGL,K,LFDPTH,NY,NX)=OMP(M,NGL,K,LFDPTH,NY,NX)+OMP1g
          ENDDO
          OSCX=OSCX+OMC1
          OSNX=OSNX+OMN1
          OSPX=OSPX+OMP1
          D2962: DO NN=1,NFGs
            tglds=JGnfA(N)-JGniA(N)+1
            DO NGL=JGniA(NN),JGnfA(NN)
              OMC1g=OMC1/tglds
              OMN1g=OMN1/tglds
              OMP1g=OMP1/tglds
              OMCff(M,NGL,LFDPTH,NY,NX)=OMCff(M,NGL,LFDPTH,NY,NX)+OMC1g*micpar%OMCA(NN)
              OMNff(M,NGL,LFDPTH,NY,NX)=OMNff(M,NGL,LFDPTH,NY,NX)+OMN1g*micpar%OMCA(NN)
              OMPff(M,NGL,LFDPTH,NY,NX)=OMPff(M,NGL,LFDPTH,NY,NX)+OMP1g*micpar%OMCA(NN)
            ENDDO
            OSCX=OSCX+OMC1*micpar%OMCA(NN)
            OSNX=OSNX+OMN1*micpar%OMCA(NN)
            OSPX=OSPX+OMP1*micpar%OMCA(NN)
          ENDDO D2962
        ENDDO D2961
      ENDDO D2960
!
!     DOC, DON AND DOP IN RESIDUE
!
!     OQC,OQN,OQP=DOC,DON,DOP in litter
!
      OQC1=AMIN1(0.1_r8*OSCX,OSCI-OSCX)
      OQN1=AMIN1(0.1_r8*OSNX,OSNI-OSNX)
      OQP1=AMIN1(0.1_r8*OSPX,OSPI-OSPX)
      OQC(K,LFDPTH,NY,NX)=OQC(K,LFDPTH,NY,NX)+OQC1
      OQN(K,LFDPTH,NY,NX)=OQN(K,LFDPTH,NY,NX)+OQN1
      OQP(K,LFDPTH,NY,NX)=OQP(K,LFDPTH,NY,NX)+OQP1
!
!     REMAINDER DISTRIBUTED TO RESIDUE FRACTIONS
!
!     OSC,OSN,OSP,OSA=SOC,SON,SOP,colonized SOC in litter
!     VOLT=litter volume
!     UORGF,UFERTN,UFERTP=accumulated litter C,N,P application
!     TNBP=accumulated net biome productivity
!
      OSCX=OSCX+OQC1
      OSNX=OSNX+OQN1
      OSPX=OSPX+OQP1
      CNOFT=0.0_r8
      CPOFT=0.0_r8
      IF(OSCI-OSCX.GT.ZEROS(NY,NX))THEN
        RNT=0.0_r8
        RPT=0.0_r8
        D965: DO M=1,jsken
          RNT=RNT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*micpar%CNOFC(M,K)
          RPT=RPT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*micpar%CPOFC(M,K)
        ENDDO D965
        FRNT=(OSNI-OSNX)/RNT
        FRPT=(OSPI-OSPX)/RPT
        D970: DO M=1,jsken
          CNOF(M)=micpar%CNOFC(M,K)*FRNT
          CPOF(M)=micpar%CPOFC(M,K)*FRPT
          CNOFT=CNOFT+CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)
          CPOFT=CPOFT+CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)
        ENDDO D970
      ELSE
        D975: DO M=1,jsken
          CNOF(M)=0.0_r8
          CPOF(M)=0.0_r8
        ENDDO D975
      ENDIF
      D2970: DO M=1,jsken
        OSC1=CFOSC(M,K,LFDPTH,NY,NX)*(OSCI-OSCX)
        IF(CNOFT.GT.ZERO)THEN
          OSN1=CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)/CNOFT*(OSNI-OSNX)
        ELSE
          OSN1=0.0_r8
        ENDIF
        IF(CPOFT.GT.ZERO)THEN
          OSP1=CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)/CPOFT*(OSPI-OSPX)
        ELSE
          OSP1=0.0_r8
        ENDIF
        OSC(M,K,LFDPTH,NY,NX)=OSC(M,K,LFDPTH,NY,NX)+OSC1
        OSA(M,K,LFDPTH,NY,NX)=OSA(M,K,LFDPTH,NY,NX)+OSC1*micpar%OMCI(1,K)
        OSN(M,K,LFDPTH,NY,NX)=OSN(M,K,LFDPTH,NY,NX)+OSN1
        OSP(M,K,LFDPTH,NY,NX)=OSP(M,K,LFDPTH,NY,NX)+OSP1
        IF(LFDPTH.EQ.0)THEN
          VOLT(LFDPTH,NY,NX)=VOLT(LFDPTH,NY,NX)+OSC1*ppmc/BKRS(micpar%k_fine_litr)
        ENDIF
      ENDDO D2970
      TORGF=TORGF+OSCI
      TORGN=TORGN+OSNI
      TORGP=TORGP+OSPI
      UORGF(NY,NX)=UORGF(NY,NX)+OSCI
      UFERTN(NY,NX)=UFERTN(NY,NX)+OSNI
      UFERTP(NY,NX)=UFERTP(NY,NX)+OSPI
      IF(IYTYP(2,I,NY,NX).LT.3)THEN
        TNBP(NY,NX)=TNBP(NY,NX)+OSCI
      ENDIF
    ENDDO D2965
  ENDIF
  end associate
  end subroutine ApplyPlantAnimalResidue
!------------------------------------------------------------------------------------------

  subroutine ApplyMineralFertilizer(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(out) :: OFC(2),OFN(2),OFP(2)
  integer, intent(out) :: LFDPTH
  real(r8) :: BAREF
  real(r8) :: CVRDF
  real(r8) :: CAC
  real(r8) :: CAS
  real(r8) :: CACX
  real(r8) :: CASX
  real(r8) :: FDPTHF
  real(r8) :: H0PO4T,H1PO4T,H2PO4T,H3PO4T
  real(r8) :: PMA,PMB,PHA
  real(r8) :: PALPOT,PFEPOT
  real(r8) :: PCAPDT,PCAPHT,PCAPMT
  real(r8) :: PMAX,PMBX,PHAX
  REAL(R8) :: XN4T
  real(r8) :: XOH0T,XOH1T,XOH2T,XH1PT,XH2PT
  real(r8) :: Z4A,Z3A,ZUA,ZOA,Z4B,Z3B
  REAL(R8) :: ZUB,ZOB
  real(r8) :: ZNH4T,ZNH3T,ZNO3T,ZNO2T
  real(r8) :: ZFE1PT,ZFE2PT
  real(r8) :: ZCA0PT,ZCA1PT,ZCA2PT
  real(r8) :: ZMG1PT,Z4AX,Z3AX,ZUAX,ZOAX
  real(r8) :: Z4BX,Z3BX,ZUBX,ZOBX
  integer :: L

!     begin_execution
!
!     NH4,NH3,UREA,NO3 FERTILIZER APPLICATION
!
!     *A,*B=broadcast,banded
!     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
!
  Z4A=FERT(1,I,NY,NX)
  Z3A=FERT(2,I,NY,NX)
  ZUA=FERT(3,I,NY,NX)
  ZOA=FERT(4,I,NY,NX)
  Z4B=FERT(5,I,NY,NX)
  Z3B=FERT(6,I,NY,NX)
  ZUB=FERT(7,I,NY,NX)
  ZOB=FERT(8,I,NY,NX)
!
!     MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE
!
!     PM*,PH*=Ca(H2PO4)2,apatite
!
  PMA=FERT(9,I,NY,NX)
  PMB=FERT(10,I,NY,NX)
  PHA=FERT(11,I,NY,NX)
!
!     LIME AND GYPSUM
!
!     CAC,CAS=CaCO3,CaSO4
!
  CAC=FERT(12,I,NY,NX)
  CAS=FERT(13,I,NY,NX)
!
!     PLANT(1) AND ANIMAL(2) RESIDUE C, N AND P
!
  OFC(1)=FERT(14,I,NY,NX)
  OFN(1)=FERT(15,I,NY,NX)
  OFP(1)=FERT(16,I,NY,NX)
  OFC(2)=FERT(17,I,NY,NX)
  OFN(2)=FERT(18,I,NY,NX)
  OFP(2)=FERT(19,I,NY,NX)
!
!     SOIL LAYER NUMBER AT DEPTH OF FERTILIZER APPLICATION
!
!     LFDPTH=layer number
!     CVRDF=fraction of fertilizer applied to surface litter
!
  IF(Z4A+Z3A+ZUA+ZOA+Z4B+Z3B+ZUB+ZOB+PMA+PMB+PHA+CAC+CAS.GT.0.0_r8)THEN
    FDPTHF=FDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
    IF(FDPTHF.LE.0.0_r8.AND.test_aeqb(Z4B+Z3B+ZUB+ZOB+PMB,0._r8))THEN
      LFDPTH=0
      CVRDF=1.0_r8-EXP(-0.8E-02_r8*(ORGC(0,NY,NX)/AREA(3,0,NY,NX)))
    ELSE
      D65: DO L=NUI(NY,NX),JZ
        IF(CDPTH(L,NY,NX).GE.FDPTHF)THEN
          LFDPTH=L
          CVRDF=1.0_r8
          exit
        ENDIF
      ENDDO D65
    ENDIF
    BAREF=1.0_r8-CVRDF
!
!     RESET WIDTH AND DEPTH OF NH4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWN=width of NH4 band row
!     DPNHB,WDNHB=depth,width of NH4 band
!     VLNHB,VLNH4=soil volume in NH4 band,non-band
!
    IF((Z4B+Z3B+ZUB.GT.0.0_r8).OR.((trc_solml(ids_NH4B,LFDPTH,NY,NX).GT.0.0_r8 &
      .OR.trc_solml(idg_NH3B,LFDPTH,NY,NX).GT.0.0_r8).AND.IFNHB(NY,NX).EQ.0))THEN
      IFNHB(NY,NX)=1
      ROWN(NY,NX)=ROWI(I,NY,NX)
      D50: DO L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          DPNHB(L,NY,NX)=DLYR(3,L,NY,NX)
          WDNHB(L,NY,NX)=0.0_r8
        ELSEIF(L.EQ.LFDPTH)THEN
          DPNHB(L,NY,NX)=AMAX1(0.025_r8,FDPTHF-CDPTH(L-1,NY,NX))
          WDNHB(L,NY,NX)=AMIN1(0.025_r8,ROWN(NY,NX))
        ELSE
          DPNHB(L,NY,NX)=0.0_r8
          WDNHB(L,NY,NX)=0.0_r8
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN(ids_NH4B,L,NY,NX)=AMIN1(0.999_r8,WDNHB(L,NY,NX)/ROWN(NY,NX) &
            *DPNHB(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          trcs_VLN(ids_NH4B,L,NY,NX)=0.0_r8
        ENDIF
        trcs_VLN(ids_NH4,L,NY,NX)=1.0_r8-trcs_VLN(ids_NH4B,L,NY,NX)
        trcs_VLN(idg_NH3B,L,NY,NX)=trcs_VLN(ids_NH4B,L,NY,NX)
        trcs_VLN(idg_NH3,L,NY,NX)=trcs_VLN(ids_NH4,L,NY,NX)
        ZNH4T=trc_solml(ids_NH4,L,NY,NX)+trc_solml(ids_NH4B,L,NY,NX)
        ZNH3T=trc_solml(idg_NH3,L,NY,NX)+trc_solml(idg_NH3B,L,NY,NX)
        XN4T=trcx_solml(idx_NH4,L,NY,NX)+trcx_solml(idx_NH4B,L,NY,NX)
        trc_solml(ids_NH4,L,NY,NX)=ZNH4T*trcs_VLN(ids_NH4,L,NY,NX)
        trc_solml(idg_NH3,L,NY,NX)=ZNH3T*trcs_VLN(idg_NH3,L,NY,NX)
        trc_solml(ids_NH4B,L,NY,NX)=ZNH4T*trcs_VLN(ids_NH4B,L,NY,NX)
        trc_solml(idg_NH3B,L,NY,NX)=ZNH3T*trcs_VLN(idg_NH3B,L,NY,NX)
        trcx_solml(idx_NH4,L,NY,NX)=XN4T*trcs_VLN(ids_NH4,L,NY,NX)
        trcx_solml(idx_NH4B,L,NY,NX)=XN4T*trcs_VLN(ids_NH4B,L,NY,NX)
      ENDDO D50
      DPNH4(NY,NX)=DPNHB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF NO3 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWO=width of NO3 band row
!     DPNOB,WDNOB=depth,width of NO3 band
!     VLNOB,VLNO3=soil volume in NO3 band,non-band
!
    IF((Z4B+Z3B+ZUB+ZOB.GT.0.0_r8).OR.((trc_solml(ids_NO3B,LFDPTH,NY,NX).GT.0.0_r8 &
      .OR.trc_solml(ids_NO2B,LFDPTH,NY,NX).GT.0.0_r8).AND.IFNOB(NY,NX).EQ.0))THEN
      IFNOB(NY,NX)=1
      ROWO(NY,NX)=ROWI(I,NY,NX)
      D45: DO L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          DPNOB(L,NY,NX)=DLYR(3,L,NY,NX)
          WDNOB(L,NY,NX)=0.0_r8
        ELSEIF(L.EQ.LFDPTH)THEN
          DPNOB(L,NY,NX)=AMAX1(0.01_r8,FDPTHF-CDPTH(L-1,NY,NX))
          WDNOB(L,NY,NX)=AMIN1(0.01_r8,ROWO(NY,NX))
        ELSE
          DPNOB(L,NY,NX)=0.0_r8
          WDNOB(L,NY,NX)=0.0_r8
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN(ids_NO3B,L,NY,NX)=AMIN1(0.999_r8,WDNOB(L,NY,NX)/ROWO(NY,NX) &
            *DPNOB(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          trcs_VLN(ids_NO3B,L,NY,NX)=0.0_r8
        ENDIF

        trcs_VLN(ids_NO3,L,NY,NX)=1.0_r8-trcs_VLN(ids_NO3B,L,NY,NX)
        trcs_VLN(ids_NO2B,L,NY,NX)=trcs_VLN(ids_NO3B,L,NY,NX)
        trcs_VLN(ids_NO2,L,NY,NX)=trcs_VLN(ids_NO3,L,NY,NX)
        ZNO3T=trc_solml(ids_NO3,L,NY,NX)+trc_solml(ids_NO3B,L,NY,NX)
        ZNO2T=trc_solml(ids_NO2,L,NY,NX)+trc_solml(ids_NO2B,L,NY,NX)

        trc_solml(ids_NO3,L,NY,NX)=ZNO3T*trcs_VLN(ids_NO3,L,NY,NX)
        trc_solml(ids_NO2,L,NY,NX)=ZNO2T*trcs_VLN(ids_NO2,L,NY,NX)
        trc_solml(ids_NO3B,L,NY,NX)=ZNO3T*trcs_VLN(ids_NO3B,L,NY,NX)
        trc_solml(ids_NO2B,L,NY,NX)=ZNO2T*trcs_VLN(ids_NO2B,L,NY,NX)
      ENDDO D45
      DPNO3(NY,NX)=DPNOB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF PO4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWP=width of H2PO4 band row
!     DPPOB,WDPOB=depth,width of H2PO4 band
!     VLPOB,VLPO4=soil volume in H2PO4 band,non-band
!
    IF((PMB.GT.0.0).OR.(trc_solml(ids_H2PO4B,LFDPTH,NY,NX).GT.0.0_r8.AND.IFPOB(NY,NX).EQ.0))THEN
      IFPOB(NY,NX)=1
      ROWP(NY,NX)=ROWI(I,NY,NX)
      DO  L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          DPPOB(L,NY,NX)=DLYR(3,L,NY,NX)
          WDPOB(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
        ELSEIF(L.EQ.LFDPTH)THEN
          DPPOB(L,NY,NX)=AMAX1(0.01,FDPTHF-CDPTH(L-1,NY,NX))
          WDPOB(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
        ELSE
          DPPOB(L,NY,NX)=0.0_r8
          WDPOB(L,NY,NX)=0.0_r8
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN(ids_H1PO4B,L,NY,NX)=AMIN1(0.999,WDPOB(L,NY,NX)/ROWP(NY,NX) &
          *DPPOB(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          trcs_VLN(ids_H1PO4B,L,NY,NX)=0.0_r8
        ENDIF
        trcs_VLN(ids_H1PO4,L,NY,NX)=1.0-trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcs_VLN(ids_H2PO4B,L,NY,NX)=trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcs_VLN(ids_H2PO4,L,NY,NX)=trcs_VLN(ids_H1PO4,L,NY,NX)

        H0PO4T=trcsa_solml(idsa_H0PO4,L,NY,NX)+trcsa_solml(idsa_H0PO4B,L,NY,NX)
        H1PO4T=trc_solml(ids_H1PO4,L,NY,NX)+trc_solml(ids_H1PO4B,L,NY,NX)
        H2PO4T=trc_solml(ids_H2PO4,L,NY,NX)+trc_solml(ids_H2PO4B,L,NY,NX)
        H3PO4T=trcsa_solml(idsa_H3PO4,L,NY,NX)+trcsa_solml(idsa_H3PO4B,L,NY,NX)
        ZFE1PT=trcsa_solml(idsa_FeHPO4,L,NY,NX)+trcsa_solml(idsa_FeHPO4B,L,NY,NX)
        ZFE2PT=trcsa_solml(idsa_FeH2PO4,L,NY,NX)+trcsa_solml(idsa_FeH2PO4B,L,NY,NX)
        ZCA0PT=trcsa_solml(idsa_CaPO4,L,NY,NX)+trcsa_solml(idsa_CaPO4B,L,NY,NX)
        ZCA1PT=trcsa_solml(idsa_CaHPO4,L,NY,NX)+trcsa_solml(idsa_CaHPO4B,L,NY,NX)
        ZCA2PT=trcsa_solml(idsa_CaH2PO4,L,NY,NX)+trcsa_solml(idsa_CaH2PO4B,L,NY,NX)
        ZMG1PT=trcsa_solml(idsa_MgHPO4,L,NY,NX)+trcsa_solml(idsa_MgHPO4B,L,NY,NX)

        XOH0T=trcx_solml(idx_OHe,L,NY,NX)+trcx_solml(idx_OHeB,L,NY,NX)
        XOH1T=trcx_solml(idx_OH,L,NY,NX)+trcx_solml(idx_OHB,L,NY,NX)
        XOH2T=trcx_solml(idx_OHp,L,NY,NX)+trcx_solml(idx_OHpB,L,NY,NX)
        XH1PT=trcx_solml(idx_HPO4,L,NY,NX)+trcx_solml(idx_HPO4B,L,NY,NX)
        XH2PT=trcx_solml(idx_H2PO4,L,NY,NX)+trcx_solml(idx_H2PO4B,L,NY,NX)
        PALPOT=trcp_salml(idsp_AlPO4,L,NY,NX)+trcp_salml(idsp_AlPO4B,L,NY,NX)
        PFEPOT=trcp_salml(idsp_FePO4,L,NY,NX)+trcp_salml(idsp_FePO4B,L,NY,NX)
        PCAPDT=trcp_salml(idsp_CaHPO4,L,NY,NX)+trcp_salml(idsp_CaHPO4B,L,NY,NX)
        PCAPHT=trcp_salml(idsp_HA,L,NY,NX)+trcp_salml(idsp_HAB,L,NY,NX)
        PCAPMT=trcp_salml(idsp_CaH2PO4,L,NY,NX)+trcp_salml(idsp_CaH2PO4B,L,NY,NX)

        trcsa_solml(idsa_H0PO4,L,NY,NX)=H0PO4T*trcs_VLN(ids_H1PO4,L,NY,NX)
        trc_solml(ids_H1PO4,L,NY,NX)=H1PO4T*trcs_VLN(ids_H1PO4,L,NY,NX)
        trc_solml(ids_H2PO4,L,NY,NX)=H2PO4T*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_H3PO4,L,NY,NX)=H3PO4T*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_FeHPO4,L,NY,NX)=ZFE1PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_FeH2PO4,L,NY,NX)=ZFE2PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_CaPO4,L,NY,NX)=ZCA0PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_CaHPO4,L,NY,NX)=ZCA1PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_CaH2PO4,L,NY,NX)=ZCA2PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_MgHPO4,L,NY,NX)=ZMG1PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcsa_solml(idsa_H0PO4B,L,NY,NX)=H0PO4T*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trc_solml(ids_H1PO4B,L,NY,NX)=H1PO4T*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trc_solml(ids_H2PO4B,L,NY,NX)=H2PO4T*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcsa_solml(idsa_H3PO4B,L,NY,NX)=H3PO4T*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcsa_solml(idsa_FeHPO4B,L,NY,NX)=ZFE1PT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcsa_solml(idsa_FeH2PO4B,L,NY,NX)=ZFE2PT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcsa_solml(idsa_CaPO4B,L,NY,NX)=ZCA0PT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcsa_solml(idsa_CaHPO4B,L,NY,NX)=ZCA1PT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcsa_solml(idsa_CaH2PO4B,L,NY,NX)=ZCA2PT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcsa_solml(idsa_MgHPO4B,L,NY,NX)=ZMG1PT*trcs_VLN(ids_H1PO4B,L,NY,NX)

        trcx_solml(idx_OHe,L,NY,NX)=XOH0T*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcx_solml(idx_OH,L,NY,NX)=XOH1T*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcx_solml(idx_OHp,L,NY,NX)=XOH2T*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcx_solml(idx_HPO4,L,NY,NX)=XH1PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcx_solml(idx_H2PO4,L,NY,NX)=XH2PT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcx_solml(idx_OHeB,L,NY,NX)=XOH0T*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcx_solml(idx_OHB,L,NY,NX)=XOH1T*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcx_solml(idx_OHpB,L,NY,NX)=XOH2T*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcx_solml(idx_HPO4B,L,NY,NX)=XH1PT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcx_solml(idx_H2PO4B,L,NY,NX)=XH2PT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcp_salml(idsp_AlPO4,L,NY,NX)=PALPOT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcp_salml(idsp_FePO4,L,NY,NX)=PFEPOT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcp_salml(idsp_CaHPO4,L,NY,NX)=PCAPDT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcp_salml(idsp_HA,L,NY,NX)=PCAPHT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcp_salml(idsp_CaH2PO4,L,NY,NX)=PCAPMT*trcs_VLN(ids_H1PO4,L,NY,NX)
        trcp_salml(idsp_AlPO4B,L,NY,NX)=PALPOT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcp_salml(idsp_FePO4B,L,NY,NX)=PFEPOT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcp_salml(idsp_CaHPO4B,L,NY,NX)=PCAPDT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcp_salml(idsp_HAB,L,NY,NX)=PCAPHT*trcs_VLN(ids_H1PO4B,L,NY,NX)
        trcp_salml(idsp_CaH2PO4B,L,NY,NX)=PCAPMT*trcs_VLN(ids_H1PO4B,L,NY,NX)
      ENDDO
      DPPO4(NY,NX)=DPPOB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
    ENDIF
!
!     UPDATE STATE VARIABLES FOR BROADCAST AND BANDED FERTILIZER
!     NH4, NH3, UREA, NO3, PO4, LIME AND GYPSUM IN SOIL
!     AND CONVERT FROM G TO MOLE
!
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=bdcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH2PO4,CaHPO4,apatite in non-band
!     PCAPMB,PCAPDB,PCAPHB=concn of precip CaH2PO4,CaHPO4,apatite in band
!     PCACO,PCASO=precipitated CaCO3,CaSO4
!
    Z4AX=Z4A*AREA(3,LFDPTH,NY,NX)/natomw
    Z3AX=Z3A*AREA(3,LFDPTH,NY,NX)/natomw
    ZUAX=ZUA*AREA(3,LFDPTH,NY,NX)/natomw
    ZOAX=ZOA*AREA(3,LFDPTH,NY,NX)/natomw
    Z4BX=Z4B*AREA(3,LFDPTH,NY,NX)/natomw
    Z3BX=Z3B*AREA(3,LFDPTH,NY,NX)/natomw
    ZUBX=ZUB*AREA(3,LFDPTH,NY,NX)/natomw
    ZOBX=ZOB*AREA(3,LFDPTH,NY,NX)/natomw
    PMAX=PMA*AREA(3,LFDPTH,NY,NX)/(2.0_r8*patomw)
    PMBX=PMB*AREA(3,LFDPTH,NY,NX)/(2.0_r8*patomw)
    PHAX=PHA*AREA(3,LFDPTH,NY,NX)/(3.0_r8*patomw)
    CACX=CAC*AREA(3,LFDPTH,NY,NX)/40.0_r8
    CASX=CAS*AREA(3,LFDPTH,NY,NX)/40.0_r8

    FertN_soil(ifert_nh4,LFDPTH,NY,NX)=FertN_soil(ifert_nh4,LFDPTH,NY,NX)+Z4AX*CVRDF
    FertN_soil(ifert_urea,LFDPTH,NY,NX)=FertN_soil(ifert_urea,LFDPTH,NY,NX)+ZUAX*CVRDF
    FertN_soil(ifert_no3,LFDPTH,NY,NX)=FertN_soil(ifert_no3,LFDPTH,NY,NX)+ZOAX*CVRDF

    FertN_band(ifert_nh4_band,LFDPTH,NY,NX)=FertN_band(ifert_nh4_band,LFDPTH,NY,NX)+Z4BX*CVRDF
    FertN_band(ifert_urea_band,LFDPTH,NY,NX)=FertN_band(ifert_urea_band,LFDPTH,NY,NX)+ZUBX*CVRDF
    FertN_band(ifert_no3_band,LFDPTH,NY,NX)=FertN_band(ifert_no3_band,LFDPTH,NY,NX)+ZOBX*CVRDF

    trcp_salml(idsp_CaH2PO4,LFDPTH,NY,NX)=trcp_salml(idsp_CaH2PO4,LFDPTH,NY,NX)+PMAX*trcs_VLN(ids_H1PO4,LFDPTH,NY,NX)*CVRDF
    trcp_salml(idsp_CaH2PO4B,LFDPTH,NY,NX)=trcp_salml(idsp_CaH2PO4B,LFDPTH,NY,NX)+PMAX*trcs_VLN(ids_H1PO4B,LFDPTH,NY,NX)*CVRDF+PMBX*CVRDF
    trcp_salml(idsp_HA,LFDPTH,NY,NX)=trcp_salml(idsp_HA,LFDPTH,NY,NX)+PHAX*trcs_VLN(ids_H1PO4,LFDPTH,NY,NX)*CVRDF
    trcp_salml(idsp_HAB,LFDPTH,NY,NX)=trcp_salml(idsp_HAB,LFDPTH,NY,NX)+PHAX*trcs_VLN(ids_H1PO4B,LFDPTH,NY,NX)*CVRDF
    IF(LFDPTH.EQ.0)THEN
      FertN_soil(ifert_nh4,NU(NY,NX),NY,NX)=FertN_soil(ifert_nh4,NU(NY,NX),NY,NX)+Z4AX*BAREF
      FertN_soil(ifert_nh3,NU(NY,NX),NY,NX)=FertN_soil(ifert_nh3,NU(NY,NX),NY,NX)+Z3AX
      FertN_soil(ifert_urea,NU(NY,NX),NY,NX)=FertN_soil(ifert_urea,NU(NY,NX),NY,NX)+ZUAX*BAREF
      FertN_soil(ifert_no3,NU(NY,NX),NY,NX)=FertN_soil(ifert_no3,NU(NY,NX),NY,NX)+ZOAX*BAREF

      FertN_band(ifert_nh4_band,NU(NY,NX),NY,NX)=FertN_band(ifert_nh4_band,NU(NY,NX),NY,NX)+Z4BX*BAREF
      FertN_band(ifert_nh3_band,NU(NY,NX),NY,NX)=FertN_band(ifert_nh3_band,NU(NY,NX),NY,NX)+Z3BX
      FertN_band(ifert_urea_band,NU(NY,NX),NY,NX)=FertN_band(ifert_urea_band,NU(NY,NX),NY,NX)+ZUBX*BAREF
      FertN_band(ifert_no3_band,NU(NY,NX),NY,NX)=FertN_band(ifert_no3_band,NU(NY,NX),NY,NX)+ZOBX*BAREF

      trcp_salml(idsp_CaH2PO4,NU(NY,NX),NY,NX)=trcp_salml(idsp_CaH2PO4,NU(NY,NX),NY,NX)+PMAX*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)*BAREF
      trcp_salml(idsp_CaH2PO4B,NU(NY,NX),NY,NX)=trcp_salml(idsp_CaH2PO4B,NU(NY,NX),NY,NX)+PMAX*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)*BAREF+PMBX*BAREF
      trcp_salml(idsp_HA,NU(NY,NX),NY,NX)=trcp_salml(idsp_HA,NU(NY,NX),NY,NX)+PHAX*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)*BAREF
      trcp_salml(idsp_HAB,NU(NY,NX),NY,NX)=trcp_salml(idsp_HAB,NU(NY,NX),NY,NX)+PHAX*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)*BAREF
    ELSE
      FertN_soil(ifert_nh3,LFDPTH,NY,NX)=FertN_soil(ifert_nh3,LFDPTH,NY,NX)+Z3AX*CVRDF
      FertN_band(ifert_nh3_band,LFDPTH,NY,NX)=FertN_band(ifert_nh3_band,LFDPTH,NY,NX)+Z3BX*CVRDF
    ENDIF
    trcp_salml(idsp_CaCO3,NU(NY,NX),NY,NX)=trcp_salml(idsp_CaCO3,NU(NY,NX),NY,NX)+CACX
    trcp_salml(idsp_CaSO4,NU(NY,NX),NY,NX)=trcp_salml(idsp_CaSO4,NU(NY,NX),NY,NX)+CASX
    TZIN=TZIN+natomw*(Z4AX+Z3AX+ZUAX+ZOAX+Z4BX+Z3BX+ZUBX+ZOBX)
    TPIN=TPIN+62.0*(PMAX+PMBX)+93.0*PHAX
    TIONIN=TIONIN+2.0*(CACX+CASX)
    UFERTN(NY,NX)=UFERTN(NY,NX)+natomw*(Z4AX+Z4BX+Z3AX+Z3BX+ZUAX+ZUBX+ZOAX+ZOBX)
    UFERTP(NY,NX)=UFERTP(NY,NX)+62.0_r8*(PMAX+PMBX)+93.0_r8*PHAX
  ENDIF
  end subroutine ApplyMineralFertilizer
!------------------------------------------------------------------------------------------

  subroutine GetChemicalConcsInSoil(NY,NX,THETPZ)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(out) :: THETPZ(JZ,JY,JX)
  integer :: L,NTG
!     begin_execution

!     CALCULATE SOIL CONCENTRATIONS OF SOLUTES, GASES
!
!     THETW,THETI,THETP=soil micropore water,ice,air concentration
!     THETPZ=soil micropore+macropore air concn for output
!

  DO L=NUI(NY,NX),NLI(NY,NX)

    IF(VOLX(L,NY,NX).LE.ZEROS(NY,NX))THEN
      THETW(L,NY,NX)=POROS(L,NY,NX)
      THETI(L,NY,NX)=0.0_r8
      THETP(L,NY,NX)=0.0_r8
    ELSE
      THETW(L,NY,NX)=AZMAX1(AMIN1(POROS(L,NY,NX),VOLW(L,NY,NX)/VOLY(L,NY,NX)))
      THETI(L,NY,NX)=AZMAX1(AMIN1(POROS(L,NY,NX),VOLI(L,NY,NX)/VOLY(L,NY,NX)))
      THETP(L,NY,NX)=AZMAX1(VOLP(L,NY,NX)/VOLY(L,NY,NX))
    ENDIF
    THETPZ(L,NY,NX)=AZMAX1(POROS(L,NY,NX)-THETW(L,NY,NX)-THETI(L,NY,NX))
!
!     GAS CONCENTRATIONS
!
!     C*G=soil gas gaseous concentration
!     C*S=soil gas aqueous concentration
!
    IF(THETP(L,NY,NX).GT.THETX)THEN
      DO NTG=idg_beg,idg_end-1
        trc_gascl(NTG,L,NY,NX)=AZMAX1(trc_gasml(NTG,L,NY,NX)/VOLP(L,NY,NX))
      ENDDO
    ELSE
      trc_gascl(idg_beg:idg_end-1,L,NY,NX)=0.0_r8
    ENDIF

    IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DO NTG=idg_beg,idg_end-1
        trc_solcl(NTG,L,NY,NX)=AZMAX1(trc_solml(NTG,L,NY,NX)/VOLW(L,NY,NX))
      ENDDO
    ELSE
      trc_solcl(idg_beg:idg_end-1,L,NY,NX)=0.0_r8
    ENDIF
!
!     CORGC=SOC concentration
!
    IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CORGC(L,NY,NX)=AMIN1(0.55E+06_r8,ORGC(L,NY,NX)/BKVL(L,NY,NX))
    ELSE
      CORGC(L,NY,NX)=0.0_r8
    ENDIF
  ENDDO
  end subroutine GetChemicalConcsInSoil
!------------------------------------------------------------------------------------------

  subroutine ZeroHourlyArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,L,NTSA
!     begin_execution
  DO L=NUI(NY,NX),NLI(NY,NX)
    FINH(L,NY,NX)=0.0_r8
    TCO2S(L,NY,NX)=0.0_r8
    TCO2P(L,NY,NX)=0.0_r8
    TCOFLA(L,NY,NX)=0.0_r8
    TCHFLA(L,NY,NX)=0.0_r8
    TLCO2P(L,NY,NX)=0.0_r8
    TUPOXP(L,NY,NX)=0.0_r8
    TUPOXS(L,NY,NX)=0.0_r8
    TUPCHS(L,NY,NX)=0.0_r8
    TUPN2S(L,NY,NX)=0.0_r8
    TUPN3S(L,NY,NX)=0.0_r8
    TUPN3B(L,NY,NX)=0.0_r8
    TUPHGS(L,NY,NX)=0.0_r8
    TOXFLA(L,NY,NX)=0.0_r8
    TCHFLA(L,NY,NX)=0.0_r8
    TN2FLA(L,NY,NX)=0.0_r8
    TNHFLA(L,NY,NX)=0.0_r8
    THGFLA(L,NY,NX)=0.0_r8
    TLOXYP(L,NY,NX)=0.0_r8
    TLCH4P(L,NY,NX)=0.0_r8
    TLN2OP(L,NY,NX)=0.0_r8
    TLNH3P(L,NY,NX)=0.0_r8
    TLH2GP(L,NY,NX)=0.0_r8
    TUPNH4(L,NY,NX)=0.0_r8
    TUPNO3(L,NY,NX)=0.0_r8
    TUPH2P(L,NY,NX)=0.0_r8
    TUPH1P(L,NY,NX)=0.0_r8
    TUPNHB(L,NY,NX)=0.0_r8
    TUPNOB(L,NY,NX)=0.0_r8
    TUPH2B(L,NY,NX)=0.0_r8
    TUPH1B(L,NY,NX)=0.0_r8
    TUPNF(L,NY,NX)=0.0_r8
    TRN4B(L,NY,NX)=0.0_r8
    TRN3B(L,NY,NX)=0.0_r8
    TRNOB(L,NY,NX)=0.0_r8
    TRN2B(L,NY,NX)=0.0_r8
    TRH1B(L,NY,NX)=0.0_r8
    TRH2B(L,NY,NX)=0.0_r8
    TRCO2(L,NY,NX)=0.0_r8
    TBCO2(L,NY,NX)=0.0_r8


    trcx_TR(idx_NH4B,L,NY,NX)=0.0_r8
    trcx_TR(idx_OHeB:idx_end,L,NY,NX)=0.0_r8

    TRXHY(L,NY,NX)=0.0_r8
    TRXAL(L,NY,NX)=0.0_r8
    TRXFE(L,NY,NX)=0.0_r8
    TRXCA(L,NY,NX)=0.0_r8
    TRXMG(L,NY,NX)=0.0_r8
    TRXNA(L,NY,NX)=0.0_r8
    TRXKA(L,NY,NX)=0.0_r8
    TRXHC(L,NY,NX)=0.0_r8
    TRXAL2(L,NY,NX)=0.0_r8
    TRXFE2(L,NY,NX)=0.0_r8

    trcp_TR(idsp_AlOH3,L,NY,NX)=0.0_r8
    trcp_TR(idsp_FeOH3,L,NY,NX)=0.0_r8
    trcp_TR(idsp_CaCO3,L,NY,NX)=0.0_r8
    trcp_TR(idsp_CaSO4,L,NY,NX)=0.0_r8
    trcp_TR(idsp_AlPO4B,L,NY,NX)=0.0_r8
    trcp_TR(idsp_FePO4B,L,NY,NX)=0.0_r8
    trcp_TR(idsp_CaHPO4B,L,NY,NX)=0.0_r8
    trcp_TR(idsp_HAB,L,NY,NX)=0.0_r8
    trcp_TR(idsp_CaH2PO4B,L,NY,NX)=0.0_r8

    trcs_XFXS(ids_beg:ids_end,L,NY,NX)=0.0_r8

    DO NTSA=idsa_beg,idsab_end
      trcsa_TR(NTSA,L,NY,NX)=0.0_r8
      trcsa_XFXS(NTSA,L,NY,NX)=0.0_r8
    ENDDO

    DO  K=1,jcplx
      XOCFXS(K,L,NY,NX)=0.0_r8
      XONFXS(K,L,NY,NX)=0.0_r8
      XOPFXS(K,L,NY,NX)=0.0_r8
      XOAFXS(K,L,NY,NX)=0.0_r8
    ENDDO
    THAW(L,NY,NX)=0.0_r8
    THAWH(L,NY,NX)=0.0_r8
    HTHAW(L,NY,NX)=0.0_r8
    trcg_XBLL(idg_beg:idg_end,L,NY,NX)=0.0_r8
    RTDNT(L,NY,NX)=0.0_r8
  ENDDO
  end subroutine ZeroHourlyArrays

!------------------------------------------------------------------------------------------

  subroutine CalGasSolubility(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer  :: L,NTG
  real(r8) :: FH2O

  L=0
  DO NTG=idg_beg,idg_end-1
    GSolbility(NTG,L,NY,NX)=gas_solubility(NTG,TCS(L,NY,NX))
  ENDDO

  DO  L=1,NL(NY,NX)+1
! S*L=solubility of gas in water
! TCS=soil temperature (oC)
! 5.56E+04_r8 := mole H2O / m3
    FH2O=5.56E+04_r8/(5.56E+04_r8+CION(L,NY,NX))
    DO NTG=idg_beg,idg_end-1
      GSolbility(NTG,L,NY,NX)=gas_solubility(NTG,TCS(L,NY,NX))/(EXP(ACTCG(NTG)*CSTR(L,NY,NX)))*FH2O
    ENDDO
  ENDDO
  end subroutine CalGasSolubility
end module Hour1Mod
