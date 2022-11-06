module StarteMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb
  use SOMDataType
  use TracerPropMod
  use TracerIDMod
  use ChemTranspDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use FertilizerDataType
  use EcosimConst
  use SnowDataType
  use SoilBGCDataType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use AqueChemDatatype
  use GridDataType
  use EcoSIMConfig
  use InitSoluteMod
  use SoluteParMod
  use SoluteChemDataType, only : solutedtype
  use EcoSiMParDataMod, only : micpar
  use ChemTracerParsMod
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  integer :: I,K,L,MM,M,NX,NY,NR1,NP2,NP3

  public :: starte
  contains

  SUBROUTINE starte(NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE INITIALIZES ALL SOIL CHEMISTRY VARIABLES
!
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  type(solutedtype)  :: solutevar
  REAL(R8) :: BKVLX
!     begin_execution
!
!     INITIALIZE CATION AND ANION CONCENTRATIONS
!     IN PRECIPITATION (K=1), IRRIGATION (K=2) AND SOIL (K=3)
!     FROM WEATHER, IRRIGATION AND SOIL FILES IN 'READS'
!
  DO   NX=NHW,NHE
    DO  NY=NVN,NVS
      solutevar%CCO2M=CCO2EI(NY,NX)/12.0_r8
      solutevar%CCH4M=AtmGgms(idg_CH4,NY,NX)/12.0_r8
      solutevar%COXYM=AtmGgms(idg_O2,NY,NX)/32.0_r8
      solutevar%CZ2GM=AtmGgms(idg_N2,NY,NX)/14.0_r8
      solutevar%CZ2OM=AtmGgms(idg_N2O,NY,NX)/14.0_r8
      solutevar%ATCA  = ATCA(NY,NX)
      solutevar%ZEROS = ZEROS(NY,NX)

      DO I=1,366
        DO  L=NU(NY,NX),NL(NY,NX)
          DO K=micpar%k_fine_litr,micpar%k_POM
            BKVLX=0._r8
!
            IF(K.EQ.micpar%k_manure.AND.L.EQ.1)THEN
!     INITIALIZE IRRIGATION WATER
              solutevar%CN4Z=CN4Q(I,NY,NX)
              solutevar%CNOZ=CNOQ(I,NY,NX)
              solutevar%CPOZ=CPOQ(I,NY,NX)
              solutevar%CALZ=CALQ(I,NY,NX)
              solutevar%CFEZ=CFEQ(I,NY,NX)
              solutevar%CCAZ=CCAQ(I,NY,NX)
              solutevar%CMGZ=CMGQ(I,NY,NX)
              solutevar%CNAZ=CNAQ(I,NY,NX)
              solutevar%CKAZ=CKAQ(I,NY,NX)
              solutevar%CSOZ=CSOQ(I,NY,NX)
              solutevar%CCLZ=CCLQ(I,NY,NX)
              solutevar%CHY1=10.0_r8**(-(PHQ(I,NY,NX)-3.0_r8))
              solutevar%COH1=DPH2O/solutevar%CHY1
            ELSE
              IF(I.EQ.1)then
                IF(K.EQ.micpar%k_fine_litr.AND.L.EQ.1)THEN
!     INITIALIZE RAINFALL
                  solutevar%CHY1=10.0_r8**(-(PHR(NY,NX)-3.0_r8))
                  solutevar%COH1=DPH2O/solutevar%CHY1
                  solutevar%CN4Z=CN4R(NY,NX)
                  solutevar%CNOZ=CNOR(NY,NX)
                  solutevar%CPOZ=CPOR(NY,NX)
                  solutevar%CALZ=CALR(NY,NX)
                  solutevar%CFEZ=CFER(NY,NX)
                  solutevar%CCAZ=CCAR(NY,NX)
                  solutevar%CMGZ=CMGR(NY,NX)
                  solutevar%CNAZ=CNAR(NY,NX)
                  solutevar%CKAZ=CKAR(NY,NX)
                  solutevar%CSOZ=CSOR(NY,NX)
                  solutevar%CCLZ=CCLR(NY,NX)
!               ELSEIF(K.EQ.micpar%k_POM.AND.(.not.is_restart_run).AND.is_first_year)THEN
!
!     INITIALIZE SOIL WATER

                  IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
                    BKVLX=BKVL(L,NY,NX)
                  ELSE
                    BKVLX=VOLA(L,NY,NX)
                  ENDIF
                  XCEC(L,NY,NX)=AMAX1(CNH4(L,NY,NX),CEC(L,NY,NX))*BKVLX
                  XAEC(L,NY,NX)=AMAX1(CPO4(L,NY,NX),AEC(L,NY,NX))*BKVLX
                  solutevar%CN4X=CNH4(L,NY,NX)
                  solutevar%CALX=CAL(L,NY,NX)
                  solutevar%CFEX=CFE(L,NY,NX)
                  solutevar%CCAX=CCA(L,NY,NX)
                  solutevar%CMGX=CMG(L,NY,NX)
                  solutevar%CNAX=CNA(L,NY,NX)
                  solutevar%CKAX=CKA(L,NY,NX)
                  solutevar%CSOX=CSO4(L,NY,NX)
                  solutevar%CCLX=CCL(L,NY,NX)
                  solutevar%CALPOX=CALPO(L,NY,NX)
                  solutevar%CFEPOX=CFEPO(L,NY,NX)
                  solutevar%CCAPDX=CCAPD(L,NY,NX)
                  solutevar%CCAPHX=CCAPH(L,NY,NX)
                  solutevar%CALOHX=CALOH(L,NY,NX)
                  solutevar%CFEOHX=CFEOH(L,NY,NX)
                  solutevar%CCACOX=CCACO(L,NY,NX)
                  solutevar%CCASOX=CCASO(L,NY,NX)
                  solutevar%CHY1=10.0_r8**(-(PH(L,NY,NX)-3.0_r8))
                  solutevar%CNOZ=CNO3(L,NY,NX)
                  solutevar%CPOZ=CPO4(L,NY,NX)
                  solutevar%COH1=DPH2O/solutevar%CHY1
                  solutevar%CN4Z=solutevar%CN4X
                  solutevar%CALZ=solutevar%CALX
                  solutevar%CFEZ=solutevar%CFEX
                  solutevar%CCAZ=solutevar%CCAX
                  solutevar%CMGZ=solutevar%CMGX
                  solutevar%CNAZ=solutevar%CNAX
                  solutevar%CKAZ=solutevar%CKAX
                  solutevar%CSOZ=solutevar%CSOX
                  solutevar%CCLZ=solutevar%CCLX
                ELSE
                  cycle
                ENDIF
              ELSE
                cycle
              ENDIF
            ENDIF

            solutevar%XAEC  = XAEC(L,NY,NX)
            solutevar%CEC   = CEC(L,NY,NX)
            solutevar%ORGC  = ORGC(L,NY,NX)
            solutevar%VLPO4 = VLPO4(L,NY,NX)
            solutevar%XCEC  = XCEC(L,NY,NX)
            solutevar%GKC4  = GKC4(L,NY,NX)
            solutevar%GKCA  = GKCA(L,NY,NX)
            solutevar%GKCH  = GKCH(L,NY,NX)
            solutevar%GKCK  = GKCK(L,NY,NX)
            solutevar%GKCN  = GKCN(L,NY,NX)
            solutevar%GKCM  = GKCM(L,NY,NX)
            solutevar%VOLW = VOLW(L,NY,NX)

            call InitSoluteModel(K,BKVLX,ISALTG,solutevar)

            call SoluteConcentrations(K,I,L,NY,NX,solutevar)
          ENDDO
        enddo
      ENDDO

      call InitialState(NY,NX)

    ENDDO
  ENDDO

  END subroutine starte
!------------------------------------------------------------------------------------------

  subroutine SoluteConcentrations(K,I,L,NY,NX,solutevar)

  implicit none
  integer, intent(in) :: K, I, L, NY, NX
  type(solutedtype), intent(in) :: solutevar
!
!     SOLUTE CONCENTRATIONS IN PRECIPITATION
!
  IF(K.EQ.micpar%k_fine_litr.AND.L.EQ.1.AND.I.EQ.1)THEN
    CCOR(NY,NX)=solutevar%CCO21
    CCHR(NY,NX)=solutevar%CCH41
    COXR(NY,NX)=solutevar%COXY1
    CNNR(NY,NX)=solutevar%CZ2G1
    CN2R(NY,NX)=solutevar%CZ2O1
    CN4R(NY,NX)=solutevar%CN41
    CN3R(NY,NX)=solutevar%CN31
    CALR(NY,NX)=solutevar%CAL1
    CFER(NY,NX)=solutevar%CFE1
    CHYR(NY,NX)=solutevar%CHY1
    CCAR(NY,NX)=solutevar%CCA1
    CMGR(NY,NX)=solutevar%CMG1
    CNAR(NY,NX)=solutevar%CNA1
    CKAR(NY,NX)=solutevar%CKA1
    COHR(NY,NX)=solutevar%COH1
    CSOR(NY,NX)=solutevar%CSO41
    CCLR(NY,NX)=solutevar%CCL1
    CC3R(NY,NX)=solutevar%CCO31
    CHCR(NY,NX)=solutevar%CHCO31
    CAL1R(NY,NX)=solutevar%CALO1
    CAL2R(NY,NX)=solutevar%CALO2
    CAL3R(NY,NX)=solutevar%CALO3
    CAL4R(NY,NX)=solutevar%CALO4
    CALSR(NY,NX)=solutevar%CALS1
    CFE1R(NY,NX)=solutevar%CFEO1
    CFE2R(NY,NX)=solutevar%CFEO2
    CFE3R(NY,NX)=solutevar%CFEO3
    CFE4R(NY,NX)=solutevar%CFEO4
    CFESR(NY,NX)=solutevar%CFES1
    CCAOR(NY,NX)=solutevar%CCAO1
    CCACR(NY,NX)=solutevar%CCAC1
    CCAHR(NY,NX)=solutevar%CCAH1
    CCASR(NY,NX)=solutevar%CCAS1
    CMGOR(NY,NX)=solutevar%CMGO1
    CMGCR(NY,NX)=solutevar%CMGC1
    CMGHR(NY,NX)=solutevar%CMGH1
    CMGSR(NY,NX)=solutevar%CMGS1
    CNACR(NY,NX)=solutevar%CNAC1
    CNASR(NY,NX)=solutevar%CNAS1
    CKASR(NY,NX)=solutevar%CKAS1
    CH0PR(NY,NX)=solutevar%CH0P1
    CH1PR(NY,NX)=solutevar%CH1P1
    CPOR(NY,NX)=solutevar%CH2P1
    CH3PR(NY,NX)=solutevar%CH3P1
    CF1PR(NY,NX)=solutevar%CF1P1
    CF2PR(NY,NX)=solutevar%CF2P1
    CC0PR(NY,NX)=solutevar%CC0P1
    CC1PR(NY,NX)=solutevar%CC1P1
    CC2PR(NY,NX)=solutevar%CC2P1
    CM1PR(NY,NX)=solutevar%CM1P1
    CSTRR(NY,NX)=solutevar%CSTR1
!
!     SOLUTE CONCENTRATIONS IN IRRIGATION
!
  ELSEIF(K.EQ.micpar%k_manure.AND.L.EQ.1)THEN
    CCOQ(NY,NX)=solutevar%CCO21
    CCHQ(NY,NX)=solutevar%CCH41
    COXQ(NY,NX)=solutevar%COXY1
    CNNQ(NY,NX)=solutevar%CZ2G1
    CN2Q(NY,NX)=solutevar%CZ2O1
    CN4Q(I,NY,NX)=solutevar%CN41
    CN3Q(I,NY,NX)=solutevar%CN31
    CALQ(I,NY,NX)=solutevar%CAL1
    CFEQ(I,NY,NX)=solutevar%CFE1
    CHYQ(I,NY,NX)=solutevar%CHY1
    CCAQ(I,NY,NX)=solutevar%CCA1
    CMGQ(I,NY,NX)=solutevar%CMG1
    CNAQ(I,NY,NX)=solutevar%CNA1
    CKAQ(I,NY,NX)=solutevar%CKA1
    COHQ(I,NY,NX)=solutevar%COH1
    CSOQ(I,NY,NX)=solutevar%CSO41
    CCLQ(I,NY,NX)=solutevar%CCL1
    CC3Q(I,NY,NX)=solutevar%CCO31
    CHCQ(I,NY,NX)=solutevar%CHCO31
    CAL1Q(I,NY,NX)=solutevar%CALO1
    CAL2Q(I,NY,NX)=solutevar%CALO2
    CAL3Q(I,NY,NX)=solutevar%CALO3
    CAL4Q(I,NY,NX)=solutevar%CALO4
    CALSQ(I,NY,NX)=solutevar%CALS1
    CFE1Q(I,NY,NX)=solutevar%CFEO1
    CFE2Q(I,NY,NX)=solutevar%CFEO2
    CFE3Q(I,NY,NX)=solutevar%CFEO3
    CFE4Q(I,NY,NX)=solutevar%CFEO4
    CFESQ(I,NY,NX)=solutevar%CFES1
    CCAOQ(I,NY,NX)=solutevar%CCAO1
    CCACQ(I,NY,NX)=solutevar%CCAC1
    CCAHQ(I,NY,NX)=solutevar%CCAH1
    CCASQ(I,NY,NX)=solutevar%CCAS1
    CMGOQ(I,NY,NX)=solutevar%CMGO1
    CMGCQ(I,NY,NX)=solutevar%CMGC1
    CMGHQ(I,NY,NX)=solutevar%CMGH1
    CMGSQ(I,NY,NX)=solutevar%CMGS1
    CNACQ(I,NY,NX)=solutevar%CNAC1
    CNASQ(I,NY,NX)=solutevar%CNAS1
    CKASQ(I,NY,NX)=solutevar%CKAS1
    CH0PQ(I,NY,NX)=solutevar%CH0P1
    CH1PQ(I,NY,NX)=solutevar%CH1P1
    CPOQ(I,NY,NX)=solutevar%CH2P1
    CH3PQ(I,NY,NX)=solutevar%CH3P1
    CF1PQ(I,NY,NX)=solutevar%CF1P1
    CF2PQ(I,NY,NX)=solutevar%CF2P1
    CC0PQ(I,NY,NX)=solutevar%CC0P1
    CC1PQ(I,NY,NX)=solutevar%CC1P1
    CC2PQ(I,NY,NX)=solutevar%CC2P1
    CM1PQ(I,NY,NX)=solutevar%CM1P1
    CSTRQ(I,NY,NX)=solutevar%CSTR1
!
!     SOLUTE CONCENTRATIONS IN SOIL
! for the POM complex, on the first day in the first year
! U means surface irrigation
  ELSEIF(K.EQ.micpar%k_POM.AND.I.EQ.1.AND.(.not.is_restart_run).AND.is_first_year)THEN
    CCOU=solutevar%CCO21
    CCHU=solutevar%CCH41
    COXU=0._r8
    CNNU=solutevar%CZ2G1
    CN2U=solutevar%CZ2O1
    CN4U(L,NY,NX)=solutevar%CN41
    CN3U(L,NY,NX)=solutevar%CN31
    CNOU(L,NY,NX)=solutevar%CNOX
    CALU(L,NY,NX)=solutevar%CAL1
    CFEU(L,NY,NX)=solutevar%CFE1
    CHYU(L,NY,NX)=solutevar%CHY1
    CCAU(L,NY,NX)=solutevar%CCA1
    CMGU(L,NY,NX)=solutevar%CMG1
    CNAU(L,NY,NX)=solutevar%CNA1
    CKAU(L,NY,NX)=solutevar%CKA1
    COHU(L,NY,NX)=solutevar%COH1
    CSOU(L,NY,NX)=solutevar%CSO41
    CCLU(L,NY,NX)=solutevar%CCL1
    CC3U(L,NY,NX)=solutevar%CCO31
    CHCU(L,NY,NX)=solutevar%CHCO31
    CAL1U(L,NY,NX)=solutevar%CALO1
    CAL2U(L,NY,NX)=solutevar%CALO2
    CAL3U(L,NY,NX)=solutevar%CALO3
    CAL4U(L,NY,NX)=solutevar%CALO4
    CALSU(L,NY,NX)=solutevar%CALS1
    CFE1U(L,NY,NX)=solutevar%CFEO1
    CFE2U(L,NY,NX)=solutevar%CFEO2
    CFE3U(L,NY,NX)=solutevar%CFEO3
    CFE4U(L,NY,NX)=solutevar%CFEO4
    CFESU(L,NY,NX)=solutevar%CFES1
    CCAOU(L,NY,NX)=solutevar%CCAO1
    CCACU(L,NY,NX)=solutevar%CCAC1
    CCAHU(L,NY,NX)=solutevar%CCAH1
    CCASU(L,NY,NX)=solutevar%CCAS1
    CMGOU(L,NY,NX)=solutevar%CMGO1
    CMGCU(L,NY,NX)=solutevar%CMGC1
    CMGHU(L,NY,NX)=solutevar%CMGH1
    CMGSU(L,NY,NX)=solutevar%CMGS1
    CNACU(L,NY,NX)=solutevar%CNAC1
    CNASU(L,NY,NX)=solutevar%CNAS1
    CKASU(L,NY,NX)=solutevar%CKAS1
    CH0PU(L,NY,NX)=solutevar%CH0P1
    CH1PU(L,NY,NX)=solutevar%CH1P1
    CH2PU(L,NY,NX)=solutevar%CH2P1
    CH3PU(L,NY,NX)=solutevar%CH3P1
    CF1PU(L,NY,NX)=solutevar%CF1P1
    CF2PU(L,NY,NX)=solutevar%CF2P1
    CC0PU(L,NY,NX)=solutevar%CC0P1
    CC1PU(L,NY,NX)=solutevar%CC1P1
    CC2PU(L,NY,NX)=solutevar%CC2P1
    CM1PU(L,NY,NX)=solutevar%CM1P1
!
!     INITIAL STATE VARIABLES FOR GAS IN SOIL
!
    CO2G(L,NY,NX)=CCO2EI(NY,NX)*VOLP(L,NY,NX)
    CH4G(L,NY,NX)=AtmGgms(idg_CH4,NY,NX)*VOLP(L,NY,NX)
    OXYG(L,NY,NX)=AtmGgms(idg_O2,NY,NX)*VOLP(L,NY,NX)
    Z2GG(L,NY,NX)=AtmGgms(idg_N2,NY,NX)*VOLP(L,NY,NX)
    Z2OG(L,NY,NX)=AtmGgms(idg_N2O,NY,NX)*VOLP(L,NY,NX)
    ZNH3G(L,NY,NX)=AtmGgms(idg_NH3,NY,NX)*VOLP(L,NY,NX)
    H2GG(L,NY,NX)=AtmGgms(idg_H2,NY,NX)*VOLP(L,NY,NX)
    IF(CDPTH(L-1,NY,NX).LT.DTBLZ(NY,NX))THEN
      OXYS(L,NY,NX)=AtmGgms(idg_O2,NY,NX)*gas_solubility(idg_O2, ATCA(NY,NX)) &
        /(EXP(AOXYX*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    ELSE
      OXYS(L,NY,NX)=0._r8
    ENDIF
    CO2S(L,NY,NX)=CCO2EI(NY,NX)*gas_solubility(idg_CO2, ATCA(NY,NX)) &
      /(EXP(ACO2X*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    CH4S(L,NY,NX)=AtmGgms(idg_CH4,NY,NX)*gas_solubility(idg_CH4, ATCA(NY,NX)) &
      /(EXP(ACH4X*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    Z2GS(L,NY,NX)=AtmGgms(idg_N2,NY,NX)*gas_solubility(idg_N2, ATCA(NY,NX)) &
      /(EXP(AN2GX*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    Z2OS(L,NY,NX)=AtmGgms(idg_N2O,NY,NX)*gas_solubility(idg_N2O, ATCA(NY,NX)) &
      /(EXP(AN2OX*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    H2GS(L,NY,NX)=AtmGgms(idg_H2,NY,NX)*gas_solubility(idg_H2, ATCA(NY,NX)) &
      /(EXP(AH2GX*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR MINERAL N AND P IN SOIL
!
    ZNH4S(L,NY,NX)=CN4U(L,NY,NX)*VOLW(L,NY,NX)*VLNH4(L,NY,NX)*14.0
    ZNH4B(L,NY,NX)=CN4U(L,NY,NX)*VOLW(L,NY,NX)*VLNHB(L,NY,NX)*14.0
    ZNH3S(L,NY,NX)=CN3U(L,NY,NX)*VOLW(L,NY,NX)*VLNH4(L,NY,NX)*14.0
    ZNH3B(L,NY,NX)=CN3U(L,NY,NX)*VOLW(L,NY,NX)*VLNHB(L,NY,NX)*14.0
    ZNO3S(L,NY,NX)=CNOU(L,NY,NX)*VOLW(L,NY,NX)*VLNO3(L,NY,NX)*14.0
    ZNO3B(L,NY,NX)=CNOU(L,NY,NX)*VOLW(L,NY,NX)*VLNOB(L,NY,NX)*14.0
    H2PO4(L,NY,NX)=CH2PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)*31.0
    H2POB(L,NY,NX)=CH2PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)*31.0
    H1PO4(L,NY,NX)=CH1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)*31.0
    H1POB(L,NY,NX)=CH1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)*31.0
    ZNO2S(L,NY,NX)=0._r8
    ZNO2B(L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS, ANIONS AND ION PAIRS IN SOIL
!
    ZAL(L,NY,NX)=CALU(L,NY,NX)*VOLW(L,NY,NX)
    ZFE(L,NY,NX)=CFEU(L,NY,NX)*VOLW(L,NY,NX)
    ZHY(L,NY,NX)=CHYU(L,NY,NX)*VOLW(L,NY,NX)
    ZCA(L,NY,NX)=CCAU(L,NY,NX)*VOLW(L,NY,NX)
    ZMG(L,NY,NX)=CMGU(L,NY,NX)*VOLW(L,NY,NX)
    ZNA(L,NY,NX)=CNAU(L,NY,NX)*VOLW(L,NY,NX)
    ZKA(L,NY,NX)=CKAU(L,NY,NX)*VOLW(L,NY,NX)
    ZOH(L,NY,NX)=COHU(L,NY,NX)*VOLW(L,NY,NX)
    ZSO4(L,NY,NX)=CSOU(L,NY,NX)*VOLW(L,NY,NX)
    ZCL(L,NY,NX)=CCLU(L,NY,NX)*VOLW(L,NY,NX)
    ZCO3(L,NY,NX)=CC3U(L,NY,NX)*VOLW(L,NY,NX)
    ZHCO3(L,NY,NX)=CHCU(L,NY,NX)*VOLW(L,NY,NX)
    ZALOH1(L,NY,NX)=CAL1U(L,NY,NX)*VOLW(L,NY,NX)
    ZALOH2(L,NY,NX)=CAL2U(L,NY,NX)*VOLW(L,NY,NX)
    ZALOH3(L,NY,NX)=CAL3U(L,NY,NX)*VOLW(L,NY,NX)
    ZALOH4(L,NY,NX)=CAL4U(L,NY,NX)*VOLW(L,NY,NX)
    ZALS(L,NY,NX)=CALSU(L,NY,NX)*VOLW(L,NY,NX)
    ZFEOH1(L,NY,NX)=CFE1U(L,NY,NX)*VOLW(L,NY,NX)
    ZFEOH2(L,NY,NX)=CFE2U(L,NY,NX)*VOLW(L,NY,NX)
    ZFEOH3(L,NY,NX)=CFE3U(L,NY,NX)*VOLW(L,NY,NX)
    ZFEOH4(L,NY,NX)=CFE4U(L,NY,NX)*VOLW(L,NY,NX)
    ZFES(L,NY,NX)=CFESU(L,NY,NX)*VOLW(L,NY,NX)
    ZCAO(L,NY,NX)=CCAOU(L,NY,NX)*VOLW(L,NY,NX)
    ZCAC(L,NY,NX)=CCACU(L,NY,NX)*VOLW(L,NY,NX)
    ZCAH(L,NY,NX)=CCAHU(L,NY,NX)*VOLW(L,NY,NX)
    ZCAS(L,NY,NX)=CCASU(L,NY,NX)*VOLW(L,NY,NX)
    ZMGO(L,NY,NX)=CMGOU(L,NY,NX)*VOLW(L,NY,NX)
    ZMGC(L,NY,NX)=CMGCU(L,NY,NX)*VOLW(L,NY,NX)
    ZMGH(L,NY,NX)=CMGHU(L,NY,NX)*VOLW(L,NY,NX)
    ZMGS(L,NY,NX)=CMGSU(L,NY,NX)*VOLW(L,NY,NX)
    ZNAC(L,NY,NX)=CNACU(L,NY,NX)*VOLW(L,NY,NX)
    ZNAS(L,NY,NX)=CNASU(L,NY,NX)*VOLW(L,NY,NX)
    ZKAS(L,NY,NX)=CKASU(L,NY,NX)*VOLW(L,NY,NX)
    H0PO4(L,NY,NX)=CH0PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    H3PO4(L,NY,NX)=CH3PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    ZFE1P(L,NY,NX)=CF1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    ZFE2P(L,NY,NX)=CF2PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    ZCA0P(L,NY,NX)=CC0PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    ZCA1P(L,NY,NX)=CC1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    ZCA2P(L,NY,NX)=CC2PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    ZMG1P(L,NY,NX)=CM1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
    H0POB(L,NY,NX)=CH0PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    H3POB(L,NY,NX)=CH3PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    ZFE1PB(L,NY,NX)=CF1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    ZFE2PB(L,NY,NX)=CF2PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    ZCA0PB(L,NY,NX)=CC0PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    ZCA1PB(L,NY,NX)=CC1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    ZCA2PB(L,NY,NX)=CC2PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
    ZMG1PB(L,NY,NX)=CM1PU(L,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR ALL MATERIAL IN SOIL MACROPORES
!
    CO2SH(L,NY,NX)=0._r8
    CH4SH(L,NY,NX)=0._r8
    OXYSH(L,NY,NX)=0._r8
    Z2GSH(L,NY,NX)=0._r8
    Z2OSH(L,NY,NX)=0._r8
    H2GSH(L,NY,NX)=0._r8
    ZNH4SH(L,NY,NX)=0._r8
    ZNH4BH(L,NY,NX)=0._r8
    ZNH3SH(L,NY,NX)=0._r8
    ZNH3BH(L,NY,NX)=0._r8
    ZNO3SH(L,NY,NX)=0._r8
    ZNO3BH(L,NY,NX)=0._r8
    ZNO2SH(L,NY,NX)=0._r8
    ZNO2BH(L,NY,NX)=0._r8
    H2PO4H(L,NY,NX)=0._r8
    H2POBH(L,NY,NX)=0._r8
    H1PO4H(L,NY,NX)=0._r8
    H1POBH(L,NY,NX)=0._r8
    ZALH(L,NY,NX)=0._r8
    ZFEH(L,NY,NX)=0._r8
    ZHYH(L,NY,NX)=0._r8
    ZCCH(L,NY,NX)=0._r8
    ZMAH(L,NY,NX)=0._r8
    ZNAH(L,NY,NX)=0._r8
    ZKAH(L,NY,NX)=0._r8
    ZOHH(L,NY,NX)=0._r8
    ZSO4H(L,NY,NX)=0._r8
    ZCLH(L,NY,NX)=0._r8
    ZCO3H(L,NY,NX)=0._r8
    ZHCO3H(L,NY,NX)=0._r8
    ZALO1H(L,NY,NX)=0._r8
    ZALO2H(L,NY,NX)=0._r8
    ZALO3H(L,NY,NX)=0._r8
    ZALO4H(L,NY,NX)=0._r8
    ZALSH(L,NY,NX)=0._r8
    ZFEO1H(L,NY,NX)=0._r8
    ZFEO2H(L,NY,NX)=0._r8
    ZFEO3H(L,NY,NX)=0._r8
    ZFEO4H(L,NY,NX)=0._r8
    ZFESH(L,NY,NX)=0._r8
    ZCAOH(L,NY,NX)=0._r8
    ZCACH(L,NY,NX)=0._r8
    ZCAHH(L,NY,NX)=0._r8
    ZCASH(L,NY,NX)=0._r8
    ZMGOH(L,NY,NX)=0._r8
    ZMGCH(L,NY,NX)=0._r8
    ZMGHH(L,NY,NX)=0._r8
    ZMGSH(L,NY,NX)=0._r8
    ZNACH(L,NY,NX)=0._r8
    ZNASH(L,NY,NX)=0._r8
    ZKASH(L,NY,NX)=0._r8
    H0PO4H(L,NY,NX)=0._r8
    H1PO4H(L,NY,NX)=0._r8
    H3PO4H(L,NY,NX)=0._r8
    ZFE1PH(L,NY,NX)=0._r8
    ZFE2PH(L,NY,NX)=0._r8
    ZCA0PH(L,NY,NX)=0._r8
    ZCA1PH(L,NY,NX)=0._r8
    ZCA2PH(L,NY,NX)=0._r8
    ZMG1PH(L,NY,NX)=0._r8
    H0POBH(L,NY,NX)=0._r8
    H1POBH(L,NY,NX)=0._r8
    H3POBH(L,NY,NX)=0._r8
    ZFE1BH(L,NY,NX)=0._r8
    ZFE2BH(L,NY,NX)=0._r8
    ZCA0BH(L,NY,NX)=0._r8
    ZCA1BH(L,NY,NX)=0._r8
    ZCA2BH(L,NY,NX)=0._r8
    ZMG1BH(L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR EXCHANGEABLE CATIONS AND ANIONS
!
    XN4(L,NY,NX)=solutevar%XN41*BKVL(L,NY,NX)*VLNH4(L,NY,NX)
    XNB(L,NY,NX)=solutevar%XN41*BKVL(L,NY,NX)*VLNHB(L,NY,NX)
    XHY(L,NY,NX)=solutevar%XHY1*BKVL(L,NY,NX)
    XAL(L,NY,NX)=solutevar%XAL1*BKVL(L,NY,NX)
    XFE(L,NY,NX)=solutevar%XFE1*BKVL(L,NY,NX)
    XCA(L,NY,NX)=solutevar%XCA1*BKVL(L,NY,NX)
    XMG(L,NY,NX)=solutevar%XMG1*BKVL(L,NY,NX)
    XNA(L,NY,NX)=solutevar%XNA1*BKVL(L,NY,NX)
    XKA(L,NY,NX)=solutevar%XKA1*BKVL(L,NY,NX)
    XHC(L,NY,NX)=solutevar%XHC1*BKVL(L,NY,NX)
    XALO2(L,NY,NX)=solutevar%XALO21*BKVL(L,NY,NX)
    XFEO2(L,NY,NX)=solutevar%XFEO21*BKVL(L,NY,NX)
    XOH0(L,NY,NX)=solutevar%XOH01*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    XOH1(L,NY,NX)=solutevar%XOH11*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    XOH2(L,NY,NX)=solutevar%XOH21*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    XH1P(L,NY,NX)=solutevar%XH1P1*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    XH2P(L,NY,NX)=solutevar%XH2P1*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    XOH0B(L,NY,NX)=solutevar%XOH01*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    XOH1B(L,NY,NX)=solutevar%XOH11*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    XOH2B(L,NY,NX)=solutevar%XOH21*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    XH1PB(L,NY,NX)=solutevar%XH1P1*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    XH2PB(L,NY,NX)=solutevar%XH2P1*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR PRECIPITATES
!
    PALOH(L,NY,NX)=solutevar%PALOH1*BKVL(L,NY,NX)
    PFEOH(L,NY,NX)=solutevar%PFEOH1*BKVL(L,NY,NX)
    PCACO(L,NY,NX)=solutevar%PCACO1*BKVL(L,NY,NX)
    PCASO(L,NY,NX)=solutevar%PCASO1*BKVL(L,NY,NX)
    PALPO(L,NY,NX)=solutevar%PALPO1*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    PFEPO(L,NY,NX)=solutevar%PFEPO1*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    PCAPD(L,NY,NX)=solutevar%PCAPD1*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    PCAPH(L,NY,NX)=solutevar%PCAPH1*BKVL(L,NY,NX)*VLPO4(L,NY,NX)
    PCAPM(L,NY,NX)=0._r8
    PALPB(L,NY,NX)=solutevar%PALPO1*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    PFEPB(L,NY,NX)=solutevar%PFEPO1*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    PCPDB(L,NY,NX)=solutevar%PCAPD1*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    PCPHB(L,NY,NX)=solutevar%PCAPH1*BKVL(L,NY,NX)*VLPOB(L,NY,NX)
    PCPMB(L,NY,NX)=0._r8
    ECND(L,NY,NX)=0._r8
    CSTR(L,NY,NX)=0._r8
    CION(L,NY,NX)=0._r8
    ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+0.5*OSN(1,2,L,NY,NX)
    OSN(1,2,L,NY,NX)=OSN(1,2,L,NY,NX)-0.5*OSN(1,2,L,NY,NX)
  ENDIF

  end subroutine SoluteConcentrations
!------------------------------------------------------------------------------------------

  subroutine InitialState(NY,NX)

  implicit none
  integer, intent(in) :: NY, NX
  real(R8) :: VOLWW
!
!     INITIAL STATE VARIABLES FOR MINERALS IN SURFACE RESIDUE
!
  IF(.not.is_restart_run.AND.is_first_year)THEN
    ZNH4S(0,NY,NX)=0._r8
    ZNH3S(0,NY,NX)=0._r8
    ZNO3S(0,NY,NX)=0._r8
    ZNO2S(0,NY,NX)=0._r8
    H2PO4(0,NY,NX)=0._r8
    H1PO4(0,NY,NX)=0._r8
    ZNH4B(0,NY,NX)=0._r8
    ZNH3B(0,NY,NX)=0._r8
    ZNO3B(0,NY,NX)=0._r8
    ZNO2B(0,NY,NX)=0._r8
    H2POB(0,NY,NX)=0._r8
    H1POB(0,NY,NX)=0._r8
    XN4(0,NY,NX)=0._r8
    XNB(0,NY,NX)=0._r8
    XOH0(0,NY,NX)=0._r8
    XOH1(0,NY,NX)=0._r8
    XOH2(0,NY,NX)=0._r8
    XH1P(0,NY,NX)=0._r8
    XH2P(0,NY,NX)=0._r8
    XNB(0,NY,NX)=0._r8
    XOH0B(0,NY,NX)=0._r8
    XOH1B(0,NY,NX)=0._r8
    XOH2B(0,NY,NX)=0._r8
    XH1PB(0,NY,NX)=0._r8
    XH2PB(0,NY,NX)=0._r8
    PALPO(0,NY,NX)=0._r8
    PFEPO(0,NY,NX)=0._r8
    PCAPD(0,NY,NX)=0._r8
    PCAPH(0,NY,NX)=0._r8
    PCAPM(0,NY,NX)=0._r8
    PALPB(0,NY,NX)=0._r8
    PFEPB(0,NY,NX)=0._r8
    PCPDB(0,NY,NX)=0._r8
    PCPHB(0,NY,NX)=0._r8
    PCPMB(0,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR MINERAL N AND P IN SNOWPACK
!
    D9985: DO L=1,JS
      IF(VHCPW(L,NY,NX).GT.VHCPWX(NY,NX))THEN
        VOLWW=VOLWSL(L,NY,NX)+VOLSSL(L,NY,NX)+VOLISL(L,NY,NX)*DENSI
        CO2W(L,NY,NX)=VOLWW*CCOR(NY,NX)
        CH4W(L,NY,NX)=VOLWW*CCHR(NY,NX)
        OXYW(L,NY,NX)=VOLWW*COXR(NY,NX)
        ZNGW(L,NY,NX)=VOLWW*CNNR(NY,NX)
        ZN2W(L,NY,NX)=VOLWW*CN2R(NY,NX)
        ZN4W(L,NY,NX)=VOLWW*CN4R(NY,NX)*natomw
        ZN3W(L,NY,NX)=VOLWW*CN3R(NY,NX)*natomw
        ZNOW(L,NY,NX)=VOLWW*CNOR(NY,NX)*natomw
        Z1PW(L,NY,NX)=VOLWW*CH1PR(NY,NX)*patomw
        ZHPW(L,NY,NX)=VOLWW*CPOR(NY,NX)*patomw
!
!     INITIAL STATE VARIABLES FOR CATIONS AND ANIONS IN SNOWPACK
!
        IF(ISALTG.NE.0)THEN
          ZALW(L,NY,NX)=VOLWW*CALR(NY,NX)
          ZFEW(L,NY,NX)=VOLWW*CFER(NY,NX)
          ZHYW(L,NY,NX)=VOLWW*CHYR(NY,NX)
          ZCAW(L,NY,NX)=VOLWW*CCAR(NY,NX)
          ZMGW(L,NY,NX)=VOLWW*CMGR(NY,NX)
          ZNAW(L,NY,NX)=VOLWW*CNAR(NY,NX)
          ZKAW(L,NY,NX)=VOLWW*CKAR(NY,NX)
          ZOHW(L,NY,NX)=VOLWW*COHR(NY,NX)
          ZSO4W(L,NY,NX)=VOLWW*CSOR(NY,NX)
          ZCLW(L,NY,NX)=VOLWW*CCLR(NY,NX)
          ZCO3W(L,NY,NX)=VOLWW*CC3R(NY,NX)
          ZHCO3W(L,NY,NX)=VOLWW*CHCR(NY,NX)
          ZALH1W(L,NY,NX)=VOLWW*CAL1R(NY,NX)
          ZALH2W(L,NY,NX)=VOLWW*CAL2R(NY,NX)
          ZALH3W(L,NY,NX)=VOLWW*CAL3R(NY,NX)
          ZALH4W(L,NY,NX)=VOLWW*CAL4R(NY,NX)
          ZALSW(L,NY,NX)=VOLWW*CALSR(NY,NX)
          ZFEH1W(L,NY,NX)=VOLWW*CFE1R(NY,NX)
          ZFEH2W(L,NY,NX)=VOLWW*CFE2R(NY,NX)
          ZFEH3W(L,NY,NX)=VOLWW*CFE3R(NY,NX)
          ZFEH4W(L,NY,NX)=VOLWW*CFE4R(NY,NX)
          ZFESW(L,NY,NX)=VOLWW*CFESR(NY,NX)
          ZCAOW(L,NY,NX)=VOLWW*CCAOR(NY,NX)
          ZCACW(L,NY,NX)=VOLWW*CCACR(NY,NX)
          ZCAHW(L,NY,NX)=VOLWW*CCAHR(NY,NX)
          ZCASW(L,NY,NX)=VOLWW*CCASR(NY,NX)
          ZMGOW(L,NY,NX)=VOLWW*CMGOR(NY,NX)
          ZMGCW(L,NY,NX)=VOLWW*CMGCR(NY,NX)
          ZMGHW(L,NY,NX)=VOLWW*CMGHR(NY,NX)
          ZMGSW(L,NY,NX)=VOLWW*CMGSR(NY,NX)
          ZNACW(L,NY,NX)=VOLWW*CNACR(NY,NX)
          ZNASW(L,NY,NX)=VOLWW*CNASR(NY,NX)
          ZKASW(L,NY,NX)=VOLWW*CKASR(NY,NX)
          H0PO4W(L,NY,NX)=VOLWW*CH0PR(NY,NX)
          H3PO4W(L,NY,NX)=VOLWW*CH3PR(NY,NX)
          ZFE1PW(L,NY,NX)=VOLWW*CF1PR(NY,NX)
          ZFE2PW(L,NY,NX)=VOLWW*CF2PR(NY,NX)
          ZCA0PW(L,NY,NX)=VOLWW*CC0PR(NY,NX)
          ZCA1PW(L,NY,NX)=VOLWW*CC1PR(NY,NX)
          ZCA2PW(L,NY,NX)=VOLWW*CC2PR(NY,NX)
          ZMG1PW(L,NY,NX)=VOLWW*CM1PR(NY,NX)
        ENDIF
      ELSE
        CO2W(L,NY,NX)=0._r8
        CH4W(L,NY,NX)=0._r8
        OXYW(L,NY,NX)=0._r8
        ZNGW(L,NY,NX)=0._r8
        ZN2W(L,NY,NX)=0._r8
        ZN4W(L,NY,NX)=0._r8
        ZN3W(L,NY,NX)=0._r8
        ZNOW(L,NY,NX)=0._r8
        Z1PW(L,NY,NX)=0._r8
        ZHPW(L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS AND ANIONS IN SNOWPACK
!
        IF(ISALTG.NE.0)THEN
          ZALW(L,NY,NX)=0._r8
          ZFEW(L,NY,NX)=0._r8
          ZHYW(L,NY,NX)=0._r8
          ZCAW(L,NY,NX)=0._r8
          ZMGW(L,NY,NX)=0._r8
          ZNAW(L,NY,NX)=0._r8
          ZKAW(L,NY,NX)=0._r8
          ZOHW(L,NY,NX)=0._r8
          ZSO4W(L,NY,NX)=0._r8
          ZCLW(L,NY,NX)=0._r8
          ZCO3W(L,NY,NX)=0._r8
          ZHCO3W(L,NY,NX)=0._r8
          ZALH1W(L,NY,NX)=0._r8
          ZALH2W(L,NY,NX)=0._r8
          ZALH3W(L,NY,NX)=0._r8
          ZALH4W(L,NY,NX)=0._r8
          ZALSW(L,NY,NX)=0._r8
          ZFEH1W(L,NY,NX)=0._r8
          ZFEH2W(L,NY,NX)=0._r8
          ZFEH3W(L,NY,NX)=0._r8
          ZFEH4W(L,NY,NX)=0._r8
          ZFESW(L,NY,NX)=0._r8
          ZCAOW(L,NY,NX)=0._r8
          ZCACW(L,NY,NX)=0._r8
          ZCAHW(L,NY,NX)=0._r8
          ZCASW(L,NY,NX)=0._r8
          ZMGOW(L,NY,NX)=0._r8
          ZMGCW(L,NY,NX)=0._r8
          ZMGHW(L,NY,NX)=0._r8
          ZMGSW(L,NY,NX)=0._r8
          ZNACW(L,NY,NX)=0._r8
          ZNASW(L,NY,NX)=0._r8
          ZKASW(L,NY,NX)=0._r8
          H0PO4W(L,NY,NX)=0._r8
          H3PO4W(L,NY,NX)=0._r8
          ZFE1PW(L,NY,NX)=0._r8
          ZFE2PW(L,NY,NX)=0._r8
          ZCA0PW(L,NY,NX)=0._r8
          ZCA1PW(L,NY,NX)=0._r8
          ZCA2PW(L,NY,NX)=0._r8
          ZMG1PW(L,NY,NX)=0._r8
        ENDIF
      ENDIF
    ENDDO D9985
  ENDIF
  end subroutine InitialState

end module StarteMod
