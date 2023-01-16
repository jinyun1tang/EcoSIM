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
      solutevar%CCO2M=CCO2EI(NY,NX)/catomw
      solutevar%CCH4M=AtmGgms(idg_CH4,NY,NX)/catomw
      solutevar%COXYM=AtmGgms(idg_O2,NY,NX)/32.0_r8
      solutevar%CZ2GM=AtmGgms(idg_N2,NY,NX)/natomw
      solutevar%CZ2OM=AtmGgms(idg_N2O,NY,NX)/natomw
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
                  trcx_solml(idx_CEC,L,NY,NX)=AMAX1(CNH4(L,NY,NX),CEC(L,NY,NX))*BKVLX
                  trcx_solml(idx_AEC,L,NY,NX)=AMAX1(CPO4(L,NY,NX),AEC(L,NY,NX))*BKVLX
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

            solutevar%XAEC  = trcx_solml(idx_AEC,L,NY,NX)
            solutevar%CEC   = CEC(L,NY,NX)
            solutevar%ORGC  = ORGC(L,NY,NX)
            solutevar%VLPO4 = trcs_VLN(ids_H1PO4,L,NY,NX)
            solutevar%XCEC  = trcx_solml(idx_CEC,L,NY,NX)
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
  associate(                       &
    iprotein  => micpar%iprotein,  &
    k_manure  => micpar%k_manure   &
  )
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
    trcn_irrig(ids_NH4,L,NY,NX)=solutevar%CN41
    trcn_irrig(idg_NH3,L,NY,NX)=solutevar%CN31
    trcn_irrig(ids_NO3,L,NY,NX)=solutevar%CNOX
    trcn_irrig(ids_NH4B,L,NY,NX)=trcn_irrig(ids_NH4,L,NY,NX)
    trcn_irrig(idg_NH3B,L,NY,NX)=trcn_irrig(idg_NH3,L,NY,NX)
    trcn_irrig(ids_NO3B,L,NY,NX)=trcn_irrig(ids_NO3,L,NY,NX)

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
    trcn_irrig(ids_H1PO4,L,NY,NX)=solutevar%CH1P1
    trcn_irrig(ids_H2PO4,L,NY,NX)=solutevar%CH2P1
    trcn_irrig(ids_H1PO4B,L,NY,NX)=trcn_irrig(ids_H1PO4,L,NY,NX)
    trcn_irrig(ids_H2PO4B,L,NY,NX)=trcn_irrig(ids_H2PO4,L,NY,NX)
    CH3PU(L,NY,NX)=solutevar%CH3P1
    CF1PU(L,NY,NX)=solutevar%CF1P1
    CF2PU(L,NY,NX)=solutevar%CF2P1
    CC0PU(L,NY,NX)=solutevar%CC0P1
    CC1PU(L,NY,NX)=solutevar%CC1P1
    CC2PU(L,NY,NX)=solutevar%CC2P1
    CM1PU(L,NY,NX)=solutevar%CM1P1
!
!   INITIAL STATE VARIABLES FOR GAS IN SOIL
!   CCO2EI is set to the first year, because AtmGgms(idg_CO2,x)
!   varies year by year, while other tracer gases are fixed year by year,
!   this is not quite right for CH4, and N2O. However, the current implementation
!   make sure the inexact restart run works. When exact restart is used, trc_gasml
!   and trc_solml will be read from restart file, so that the following inconsistent
!   use between CO2 and other gas tracers can be avoided.
!   Comment by Jinyun Tang, Nov 11, 2022
!
    trc_gasml(idg_CO2,L,NY,NX)=CCO2EI(NY,NX)*VOLP(L,NY,NX)
    trc_gasml(idg_CH4,L,NY,NX)=AtmGgms(idg_CH4,NY,NX)*VOLP(L,NY,NX)
    trc_gasml(idg_O2,L,NY,NX)=AtmGgms(idg_O2,NY,NX)*VOLP(L,NY,NX)
    trc_gasml(idg_N2,L,NY,NX)=AtmGgms(idg_N2,NY,NX)*VOLP(L,NY,NX)
    trc_gasml(idg_N2O,L,NY,NX)=AtmGgms(idg_N2O,NY,NX)*VOLP(L,NY,NX)
    trc_gasml(idg_NH3,L,NY,NX)=AtmGgms(idg_NH3,NY,NX)*VOLP(L,NY,NX)
    trc_gasml(idg_H2,L,NY,NX)=AtmGgms(idg_H2,NY,NX)*VOLP(L,NY,NX)

!   DTBLZ: external water table depth
    IF(CDPTH(L-1,NY,NX).LT.DTBLZ(NY,NX))THEN
! above water table
      trc_solml(idg_O2,L,NY,NX)=AtmGgms(idg_O2,NY,NX)*gas_solubility(idg_O2, ATCA(NY,NX)) &
        /(EXP(ACTCG(idg_O2)*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    ELSE
!below water table
      trc_solml(idg_O2,L,NY,NX)=0._r8
    ENDIF

    trc_solml(idg_CO2,L,NY,NX)=CCO2EI(NY,NX)*gas_solubility(idg_CO2, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_CO2)*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    trc_solml(idg_CH4,L,NY,NX)=AtmGgms(idg_CH4,NY,NX)*gas_solubility(idg_CH4, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_CH4)*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    trc_solml(idg_N2,L,NY,NX)=AtmGgms(idg_N2,NY,NX)*gas_solubility(idg_N2, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_N2)*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    trc_solml(idg_N2O,L,NY,NX)=AtmGgms(idg_N2O,NY,NX)*gas_solubility(idg_N2O, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_N2O)*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
    trc_solml(idg_H2,L,NY,NX)=AtmGgms(idg_H2,NY,NX)*gas_solubility(idg_H2, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_H2)*solutevar%CSTR1))*solutevar%FH2O*VOLW(L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR MINERAL N AND P IN SOIL
!
    trc_solml(ids_NH4,L,NY,NX)=trcn_irrig(ids_NH4,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    trc_solml(ids_NH4B,L,NY,NX)=trcn_irrig(ids_NH4B,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    trc_solml(idg_NH3,L,NY,NX)=trcn_irrig(idg_NH3,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    trc_solml(idg_NH3B,L,NY,NX)=trcn_irrig(idg_NH3B,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    trc_solml(ids_NO3,L,NY,NX)=trcn_irrig(ids_NO3,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)*natomw
    trc_solml(ids_NO3B,L,NY,NX)=trcn_irrig(ids_NO3B,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)*natomw

    trc_solml(ids_H2PO4,L,NY,NX)=trcn_irrig(ids_H2PO4,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    trc_solml(ids_H2PO4B,L,NY,NX)=trcn_irrig(ids_H2PO4B,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw
    trc_solml(ids_H1PO4,L,NY,NX)=trcn_irrig(ids_H1PO4,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    trc_solml(ids_H1PO4B,L,NY,NX)=trcn_irrig(ids_H1PO4B,L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw

    trc_solml(ids_NO2,L,NY,NX)=0._r8
    trc_solml(ids_NO2B,L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS, ANIONS AND ION PAIRS IN SOIL
!
    trcsa_solml(idsa_Al,L,NY,NX)=CALU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_Fe,L,NY,NX)=CFEU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_Hp,L,NY,NX)=CHYU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_Ca,L,NY,NX)=CCAU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_Mg,L,NY,NX)=CMGU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_Na,L,NY,NX)=CNAU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_K,L,NY,NX)=CKAU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_OH,L,NY,NX)=COHU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_SO4,L,NY,NX)=CSOU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_Cl,L,NY,NX)=CCLU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_CO3,L,NY,NX)=CC3U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_HCO3,L,NY,NX)=CHCU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_AlOH,L,NY,NX)=CAL1U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_AlOH2,L,NY,NX)=CAL2U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_AlOH3,L,NY,NX)=CAL3U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_AlOH4,L,NY,NX)=CAL4U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_AlSO4,L,NY,NX)=CALSU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_FeOH,L,NY,NX)=CFE1U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_FeOH2,L,NY,NX)=CFE2U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_FeOH3,L,NY,NX)=CFE3U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_FeOH4,L,NY,NX)=CFE4U(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_FeSO4,L,NY,NX)=CFESU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_CaOH2,L,NY,NX)=CCAOU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_CaCO3,L,NY,NX)=CCACU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_CaHCO3,L,NY,NX)=CCAHU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_CaSO4,L,NY,NX)=CCASU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_MgOH2,L,NY,NX)=CMGOU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_MgCO3,L,NY,NX)=CMGCU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_MgHCO3,L,NY,NX)=CMGHU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_MgSO4,L,NY,NX)=CMGSU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_NaCO3,L,NY,NX)=CNACU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_NaSO4,L,NY,NX)=CNASU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_KSO4,L,NY,NX)=CKASU(L,NY,NX)*VOLW(L,NY,NX)
    trcsa_solml(idsa_H0PO4,L,NY,NX)=CH0PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_H3PO4,L,NY,NX)=CH3PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_FeHPO4,L,NY,NX)=CF1PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_FeH2PO4,L,NY,NX)=CF2PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_CaPO4,L,NY,NX)=CC0PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_CaHPO4,L,NY,NX)=CC1PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_CaH2PO4,L,NY,NX)=CC2PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_MgHPO4,L,NY,NX)=CM1PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcsa_solml(idsa_H0PO4B,L,NY,NX)=CH0PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_solml(idsa_H3PO4B,L,NY,NX)=CH3PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_solml(idsa_FeHPO4B,L,NY,NX)=CF1PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_solml(idsa_FeH2PO4B,L,NY,NX)=CF2PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_solml(idsa_CaPO4B,L,NY,NX)=CC0PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_solml(idsa_CaHPO4B,L,NY,NX)=CC1PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_solml(idsa_CaH2PO4B,L,NY,NX)=CC2PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcsa_solml(idsa_MgHPO4B,L,NY,NX)=CM1PU(L,NY,NX)*VOLW(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR ALL MATERIAL IN SOIL MACROPORES
!
    trc_soHml(ids_beg:ids_end,L,NY,NX)=0._r8
    trcsa_soHml(idsa_beg:idsab_end,L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR EXCHANGEABLE CATIONS AND ANIONS
!
    trcx_solml(idx_NH4,L,NY,NX)=solutevar%XN41*BKVL(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)
    trcx_solml(idx_NH4B,L,NY,NX)=solutevar%XN41*BKVL(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)
    trcx_solml(idx_Hp,L,NY,NX)=solutevar%XHY1*BKVL(L,NY,NX)
    trcx_solml(idx_Al,L,NY,NX)=solutevar%XAL1*BKVL(L,NY,NX)
    trcx_solml(idx_Fe,L,NY,NX)=solutevar%XFE1*BKVL(L,NY,NX)
    trcx_solml(idx_Ca,L,NY,NX)=solutevar%XCA1*BKVL(L,NY,NX)
    trcx_solml(idx_Mg,L,NY,NX)=solutevar%XMG1*BKVL(L,NY,NX)
    trcx_solml(idx_Na,L,NY,NX)=solutevar%XNA1*BKVL(L,NY,NX)
    trcx_solml(idx_K,L,NY,NX)=solutevar%XKA1*BKVL(L,NY,NX)
    trcx_solml(idx_COOH,L,NY,NX)=solutevar%XHC1*BKVL(L,NY,NX)
    trcx_solml(idx_AlOH2,L,NY,NX)=solutevar%XALO21*BKVL(L,NY,NX)
    trcx_solml(idx_FeOH2,L,NY,NX)=solutevar%XFEO21*BKVL(L,NY,NX)
    trcx_solml(idx_OHe,L,NY,NX)=solutevar%XOH01*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_OH,L,NY,NX)=solutevar%XOH11*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_OHp,L,NY,NX)=solutevar%XOH21*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_HPO4,L,NY,NX)=solutevar%XH1P1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_H2PO4,L,NY,NX)=solutevar%XH2P1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_OHeB,L,NY,NX)=solutevar%XOH01*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_OHB,L,NY,NX)=solutevar%XOH11*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_OHpB,L,NY,NX)=solutevar%XOH21*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_HPO4B,L,NY,NX)=solutevar%XH1P1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_H2PO4B,L,NY,NX)=solutevar%XH2P1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR PRECIPITATES
!
    trcp_salml(idsp_AlOH3,L,NY,NX)=solutevar%PALOH1*BKVL(L,NY,NX)
    trcp_salml(idsp_FeOH3,L,NY,NX)=solutevar%PFEOH1*BKVL(L,NY,NX)
    trcp_salml(idsp_CaCO3,L,NY,NX)=solutevar%PCACO1*BKVL(L,NY,NX)
    trcp_salml(idsp_CaSO4,L,NY,NX)=solutevar%PCASO1*BKVL(L,NY,NX)
    trcp_salml(idsp_AlPO4,L,NY,NX)=solutevar%PALPO1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_FePO4,L,NY,NX)=solutevar%PFEPO1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_CaHPO4,L,NY,NX)=solutevar%PCAPD1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_HA,L,NY,NX)=solutevar%PCAPH1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_CaH2PO4,L,NY,NX)=0._r8
    trcp_salml(idsp_AlPO4B,L,NY,NX)=solutevar%PALPO1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcp_salml(idsp_FePO4B,L,NY,NX)=solutevar%PFEPO1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcp_salml(idsp_CaHPO4B,L,NY,NX)=solutevar%PCAPD1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcp_salml(idsp_HAB,L,NY,NX)=solutevar%PCAPH1*BKVL(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcp_salml(idsp_CaH2PO4B,L,NY,NX)=0._r8
    ECND(L,NY,NX)=0._r8
    CSTR(L,NY,NX)=0._r8
    CION(L,NY,NX)=0._r8
! the following line is quite interesting, Jinyun Tang, Nov 17, 2022
    trc_solml(ids_NH4,L,NY,NX)=trc_solml(ids_NH4,L,NY,NX)+0.5_r8*OSN(iprotein,k_manure,L,NY,NX)
    OSN(iprotein,k_manure,L,NY,NX)=OSN(iprotein,k_manure,L,NY,NX)-0.5_r8*OSN(iprotein,k_manure,L,NY,NX)
  ENDIF
  end associate
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
    trc_solml(ids_nuts_beg:ids_nuts_end,0,NY,NX)=0._r8

    trcx_solml(idx_NH4,0,NY,NX)=0._r8
    trcx_solml(idx_NH4B,0,NY,NX)=0._r8
    trcx_solml(idx_OHe,0,NY,NX)=0._r8
    trcx_solml(idx_OH,0,NY,NX)=0._r8
    trcx_solml(idx_OHp,0,NY,NX)=0._r8
    trcx_solml(idx_HPO4,0,NY,NX)=0._r8
    trcx_solml(idx_H2PO4,0,NY,NX)=0._r8
    trcx_solml(idx_NH4B,0,NY,NX)=0._r8
    trcx_solml(idx_OHeB,0,NY,NX)=0._r8
    trcx_solml(idx_OHB,0,NY,NX)=0._r8
    trcx_solml(idx_OHpB,0,NY,NX)=0._r8
    trcx_solml(idx_HPO4B,0,NY,NX)=0._r8
    trcx_solml(idx_H2PO4B,0,NY,NX)=0._r8

    trcp_salml(idsp_AlPO4,0,NY,NX)=0._r8
    trcp_salml(idsp_FePO4,0,NY,NX)=0._r8
    trcp_salml(idsp_CaHPO4,0,NY,NX)=0._r8
    trcp_salml(idsp_HA,0,NY,NX)=0._r8
    trcp_salml(idsp_CaH2PO4,0,NY,NX)=0._r8
    trcp_salml(idsp_AlPO4B,0,NY,NX)=0._r8
    trcp_salml(idsp_FePO4B,0,NY,NX)=0._r8
    trcp_salml(idsp_CaHPO4B,0,NY,NX)=0._r8
    trcp_salml(idsp_HAB,0,NY,NX)=0._r8
    trcp_salml(idsp_CaH2PO4B,0,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR MINERAL N AND P IN SNOWPACK
!
    D9985: DO L=1,JS
      IF(VHCPW(L,NY,NX).GT.VHCPWX(NY,NX))THEN
        VOLWW=VOLWSL(L,NY,NX)+VOLSSL(L,NY,NX)+VOLISL(L,NY,NX)*DENSI
        trcg_solsml(idg_CO2,L,NY,NX)=VOLWW*CCOR(NY,NX)
        trcg_solsml(idg_CH4,L,NY,NX)=VOLWW*CCHR(NY,NX)
        trcg_solsml(idg_O2,L,NY,NX)=VOLWW*COXR(NY,NX)
        trcg_solsml(idg_N2,L,NY,NX)=VOLWW*CNNR(NY,NX)
        trcg_solsml(idg_N2O,L,NY,NX)=VOLWW*CN2R(NY,NX)
        trcg_solsml(idg_NH3,L,NY,NX)=VOLWW*CN3R(NY,NX)*natomw

        trcn_solsml(ids_NH4,L,NY,NX)=VOLWW*CN4R(NY,NX)*natomw
        trcn_solsml(ids_NO3,L,NY,NX)=VOLWW*CNOR(NY,NX)*natomw
        trcn_solsml(ids_H1PO4,L,NY,NX)=VOLWW*CH1PR(NY,NX)*patomw
        trcn_solsml(ids_H2PO4,L,NY,NX)=VOLWW*CPOR(NY,NX)*patomw
!
!     INITIAL STATE VARIABLES FOR CATIONS AND ANIONS IN SNOWPACK
!
        IF(salt_model)THEN
          trcs_solsml(idsa_Al,L,NY,NX)=VOLWW*CALR(NY,NX)
          trcs_solsml(idsa_Fe,L,NY,NX)=VOLWW*CFER(NY,NX)
          trcs_solsml(idsa_Hp,L,NY,NX)=VOLWW*CHYR(NY,NX)
          trcs_solsml(idsa_Ca,L,NY,NX)=VOLWW*CCAR(NY,NX)
          trcs_solsml(idsa_Mg,L,NY,NX)=VOLWW*CMGR(NY,NX)
          trcs_solsml(idsa_Na,L,NY,NX)=VOLWW*CNAR(NY,NX)
          trcs_solsml(idsa_K,L,NY,NX)=VOLWW*CKAR(NY,NX)
          trcs_solsml(idsa_OH,L,NY,NX)=VOLWW*COHR(NY,NX)
          trcs_solsml(idsa_SO4,L,NY,NX)=VOLWW*CSOR(NY,NX)
          trcs_solsml(idsa_Cl,L,NY,NX)=VOLWW*CCLR(NY,NX)
          trcs_solsml(idsa_CO3,L,NY,NX)=VOLWW*CC3R(NY,NX)
          trcs_solsml(idsa_HCO3,L,NY,NX)=VOLWW*CHCR(NY,NX)
          trcs_solsml(idsa_AlOH,L,NY,NX)=VOLWW*CAL1R(NY,NX)
          trcs_solsml(idsa_AlOH2,L,NY,NX)=VOLWW*CAL2R(NY,NX)
          trcs_solsml(idsa_AlOH3,L,NY,NX)=VOLWW*CAL3R(NY,NX)
          trcs_solsml(idsa_AlOH4,L,NY,NX)=VOLWW*CAL4R(NY,NX)
          trcs_solsml(idsa_AlSO4,L,NY,NX)=VOLWW*CALSR(NY,NX)
          trcs_solsml(idsa_FeOH,L,NY,NX)=VOLWW*CFE1R(NY,NX)
          trcs_solsml(idsa_FeOH2,L,NY,NX)=VOLWW*CFE2R(NY,NX)
          trcs_solsml(idsa_FeOH3,L,NY,NX)=VOLWW*CFE3R(NY,NX)
          trcs_solsml(idsa_FeOH4,L,NY,NX)=VOLWW*CFE4R(NY,NX)
          trcs_solsml(idsa_FeSO4,L,NY,NX)=VOLWW*CFESR(NY,NX)
          trcs_solsml(idsa_CaOH2,L,NY,NX)=VOLWW*CCAOR(NY,NX)
          trcs_solsml(idsa_CaCO3,L,NY,NX)=VOLWW*CCACR(NY,NX)
          trcs_solsml(idsa_CaHCO3,L,NY,NX)=VOLWW*CCAHR(NY,NX)
          trcs_solsml(idsa_CaSO4,L,NY,NX)=VOLWW*CCASR(NY,NX)
          trcs_solsml(idsa_MgOH2,L,NY,NX)=VOLWW*CMGOR(NY,NX)
          trcs_solsml(idsa_MgCO3,L,NY,NX)=VOLWW*CMGCR(NY,NX)
          trcs_solsml(idsa_MgHCO3,L,NY,NX)=VOLWW*CMGHR(NY,NX)
          trcs_solsml(idsa_MgSO4,L,NY,NX)=VOLWW*CMGSR(NY,NX)
          trcs_solsml(idsa_NaCO3,L,NY,NX)=VOLWW*CNACR(NY,NX)
          trcs_solsml(idsa_NaSO4,L,NY,NX)=VOLWW*CNASR(NY,NX)
          trcs_solsml(idsa_KSO4,L,NY,NX)=VOLWW*CKASR(NY,NX)
          trcs_solsml(idsa_H0PO4,L,NY,NX)=VOLWW*CH0PR(NY,NX)
          trcs_solsml(idsa_H3PO4,L,NY,NX)=VOLWW*CH3PR(NY,NX)
          trcs_solsml(idsa_FeHPO4,L,NY,NX)=VOLWW*CF1PR(NY,NX)
          trcs_solsml(idsa_FeH2PO4,L,NY,NX)=VOLWW*CF2PR(NY,NX)
          trcs_solsml(idsa_CaPO4,L,NY,NX)=VOLWW*CC0PR(NY,NX)
          trcs_solsml(idsa_CaHPO4,L,NY,NX)=VOLWW*CC1PR(NY,NX)
          trcs_solsml(idsa_CaH2PO4,L,NY,NX)=VOLWW*CC2PR(NY,NX)
          trcs_solsml(idsa_MgHPO4,L,NY,NX)=VOLWW*CM1PR(NY,NX)
        ENDIF
      ELSE
        trcg_solsml(idg_beg:idg_end-1,L,NY,NX)=0._r8
        trcn_solsml(ids_nut_beg:ids_nuts_end,L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS AND ANIONS IN SNOWPACK
!
        IF(salt_model)THEN
          trcs_solsml(idsa_beg:idsa_end,L,NY,NX)=0._r8
        ENDIF
      ENDIF
    ENDDO D9985
  ENDIF
  end subroutine InitialState

end module StarteMod
