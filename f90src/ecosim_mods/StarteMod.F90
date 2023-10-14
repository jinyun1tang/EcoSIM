module StarteMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
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
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  integer :: I,K,L,MM,M,NX,NY,NR1,NP2,NP3

  public :: starte
  contains

  SUBROUTINE starte(NHW,NHE,NVN,NVS)
!
! DESCRIPTION:
!     THIS SUBROUTINE INITIALIZES ALL SOIL CHEMISTRY VARIABLES
! The top layer are initialized every year to accomadate changes in
! boundary conditions (irrigation, rainfall, manure application) every year
! other layers are done only in the first year.

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  type(solutedtype)  :: solutevar
  REAL(R8) :: BulkSoilMass
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
            BulkSoilMass=0._r8
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
              solutevar%H_1p_conc=10.0_r8**(-(PHQ(I,NY,NX)-3.0_r8))
              solutevar%OH_1e_conc=DPH2O/solutevar%H_1p_conc
            ELSE
              IF(I.EQ.1)then
                IF(K.EQ.micpar%k_fine_litr.AND.L.EQ.1)THEN
!     INITIALIZE RAINFALL, top layer
                  solutevar%H_1p_conc=10.0_r8**(-(PHR(NY,NX)-3.0_r8))
                  solutevar%OH_1e_conc=DPH2O/solutevar%H_1p_conc
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
!               ELSEIF(K.EQ.micpar%k_POM.AND.(.not.is_restart()).AND.is_first_year)THEN
!
!     INITIALIZE SOIL WATER

                  IF(SoilMicPMassLayer(L,NY,NX).GT.ZEROS(NY,NX))THEN
                    BulkSoilMass=SoilMicPMassLayer(L,NY,NX)
                  ELSE
                    !this is for water
                    BulkSoilMass=VLMicP(L,NY,NX)
                  ENDIF
                  trcx_solml(idx_CEC,L,NY,NX)=AMAX1(CNH4(L,NY,NX),CEC(L,NY,NX))*BulkSoilMass
                  trcx_solml(idx_AEC,L,NY,NX)=AMAX1(CPO4(L,NY,NX),AEC(L,NY,NX))*BulkSoilMass
                  solutevar%CN4X=CNH4(L,NY,NX)
                  solutevar%CALX=CAL(L,NY,NX)
                  solutevar%CFEX=CFE(L,NY,NX)
                  solutevar%CaX_conc=CCA(L,NY,NX)
                  solutevar%MgX_conc=CMG(L,NY,NX)
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
                  solutevar%H_1p_conc=10.0_r8**(-(PH(L,NY,NX)-3.0_r8))
                  solutevar%CNOZ=CNO3(L,NY,NX)
                  solutevar%CPOZ=CPO4(L,NY,NX)
                  solutevar%OH_1e_conc=DPH2O/solutevar%H_1p_conc
                  solutevar%CN4Z=solutevar%CN4X
                  solutevar%CALZ=solutevar%CALX
                  solutevar%CFEZ=solutevar%CFEX
                  solutevar%CCAZ=solutevar%CaX_conc
                  solutevar%CMGZ=solutevar%MgX_conc
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
            solutevar%VLWatMicP = VLWatMicP(L,NY,NX)

            call InitSoluteModel(K,BulkSoilMass,ISALTG,solutevar)

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
  !litter pool, top soil layer
    CCOR(NY,NX)=solutevar%H2CO3_aqu_conc
    CCHR(NY,NX)=solutevar%CH4_aqu_conc
    COXR(NY,NX)=solutevar%O2_aqu_conc
    CNNR(NY,NX)=solutevar%N2_aqu_conc
    CN2R(NY,NX)=solutevar%N2O_aqu_conc
    CN4R(NY,NX)=solutevar%NH4_1p_conc
    CN3R(NY,NX)=solutevar%NH3_aqu_conc
    CALR(NY,NX)=solutevar%Al_3p_conc
    CFER(NY,NX)=solutevar%Fe_3p_conc
    CHYR(NY,NX)=solutevar%H_1p_conc
    CCAR(NY,NX)=solutevar%Ca_2p_conc
    CMGR(NY,NX)=solutevar%Mg_2p_conc
    CNAR(NY,NX)=solutevar%Na_1p_conc
    CKAR(NY,NX)=solutevar%K_1p_conc
    COHR(NY,NX)=solutevar%OH_1e_conc
    CSOR(NY,NX)=solutevar%SO4_2e_conc
    CCLR(NY,NX)=solutevar%Cl_e_conc
    CC3R(NY,NX)=solutevar%CO3_2e_conc
    CHCR(NY,NX)=solutevar%HCO3_e_conc
    CAL1R(NY,NX)=solutevar%AlOH_2p_conc
    CAL2R(NY,NX)=solutevar%AlOH2_p_conc
    CAL3R(NY,NX)=solutevar%AlOH3_conc
    CAL4R(NY,NX)=solutevar%AlOH4_1e_conc
    CALSR(NY,NX)=solutevar%AlSO4_1p_conc
    CFE1R(NY,NX)=solutevar%FeOH_2p_conc
    CFE2R(NY,NX)=solutevar%FeO2H2_p_conc
    CFE3R(NY,NX)=solutevar%FeO3H3_conc
    CFE4R(NY,NX)=solutevar%FeO4H4_1e_conc
    CFESR(NY,NX)=solutevar%FeSO4_1p_conc
    CCAOR(NY,NX)=solutevar%CaO2H2_conc
    CCACR(NY,NX)=solutevar%CaCO3_conc
    CCAHR(NY,NX)=solutevar%CaHCO3_1p_conc
    CCASR(NY,NX)=solutevar%CaSO4_conc
    CMGOR(NY,NX)=solutevar%MgOH_1p_conc
    CMGCR(NY,NX)=solutevar%MgCO3_conc
    CMGHR(NY,NX)=solutevar%MgHCO3_1p_conc
    CMGSR(NY,NX)=solutevar%MgSO4_conc
    CNACR(NY,NX)=solutevar%NaCO3_1e_conc
    CNASR(NY,NX)=solutevar%NaSO4_1e_conc
    CKASR(NY,NX)=solutevar%KSO4_1e_conc
    CH0PR(NY,NX)=solutevar%H0PO4_conc
    CH1PR(NY,NX)=solutevar%H1PO4_2e_conc
    CPOR(NY,NX)=solutevar%H2PO4_1e_conc
    CH3PR(NY,NX)=solutevar%H3PO4_conc
    CF1PR(NY,NX)=solutevar%FeHPO4_conc
    CF2PR(NY,NX)=solutevar%FeH2PO4_conc
    CC0PR(NY,NX)=solutevar%CaPO4_1e_con
    CC1PR(NY,NX)=solutevar%CaHPO4_conc
    CC2PR(NY,NX)=solutevar%CaH2PO4_1p_conc
    CM1PR(NY,NX)=solutevar%MgHPO4_conc
    CSTRR(NY,NX)=solutevar%CSTR1
!
!     SOLUTE CONCENTRATIONS IN IRRIGATION
!
  ELSEIF(K.EQ.micpar%k_manure.AND.L.EQ.1)THEN
  ! manure, top layer
    CCOQ(NY,NX)=solutevar%H2CO3_aqu_conc
    CCHQ(NY,NX)=solutevar%CH4_aqu_conc
    COXQ(NY,NX)=solutevar%O2_aqu_conc
    CNNQ(NY,NX)=solutevar%N2_aqu_conc
    CN2Q(NY,NX)=solutevar%N2O_aqu_conc
    CN4Q(I,NY,NX)=solutevar%NH4_1p_conc
    CN3Q(I,NY,NX)=solutevar%NH3_aqu_conc
    CALQ(I,NY,NX)=solutevar%Al_3p_conc
    CFEQ(I,NY,NX)=solutevar%Fe_3p_conc
    CHYQ(I,NY,NX)=solutevar%H_1p_conc
    CCAQ(I,NY,NX)=solutevar%Ca_2p_conc
    CMGQ(I,NY,NX)=solutevar%Mg_2p_conc
    CNAQ(I,NY,NX)=solutevar%Na_1p_conc
    CKAQ(I,NY,NX)=solutevar%K_1p_conc
    COHQ(I,NY,NX)=solutevar%OH_1e_conc
    CSOQ(I,NY,NX)=solutevar%SO4_2e_conc
    CCLQ(I,NY,NX)=solutevar%Cl_e_conc
    CC3Q(I,NY,NX)=solutevar%CO3_2e_conc
    CHCQ(I,NY,NX)=solutevar%HCO3_e_conc
    CAL1Q(I,NY,NX)=solutevar%AlOH_2p_conc
    CAL2Q(I,NY,NX)=solutevar%AlOH2_p_conc
    CAL3Q(I,NY,NX)=solutevar%AlOH3_conc
    CAL4Q(I,NY,NX)=solutevar%AlOH4_1e_conc
    CALSQ(I,NY,NX)=solutevar%AlSO4_1p_conc
    CFE1Q(I,NY,NX)=solutevar%FeOH_2p_conc
    CFE2Q(I,NY,NX)=solutevar%FeO2H2_p_conc
    CFE3Q(I,NY,NX)=solutevar%FeO3H3_conc
    CFE4Q(I,NY,NX)=solutevar%FeO4H4_1e_conc
    CFESQ(I,NY,NX)=solutevar%FeSO4_1p_conc
    CCAOQ(I,NY,NX)=solutevar%CaO2H2_conc
    CCACQ(I,NY,NX)=solutevar%CaCO3_conc
    CCAHQ(I,NY,NX)=solutevar%CaHCO3_1p_conc
    CCASQ(I,NY,NX)=solutevar%CaSO4_conc
    CMGOQ(I,NY,NX)=solutevar%MgOH_1p_conc
    CMGCQ(I,NY,NX)=solutevar%MgCO3_conc
    CMGHQ(I,NY,NX)=solutevar%MgHCO3_1p_conc
    CMGSQ(I,NY,NX)=solutevar%MgSO4_conc
    CNACQ(I,NY,NX)=solutevar%NaCO3_1e_conc
    CNASQ(I,NY,NX)=solutevar%NaSO4_1e_conc
    CKASQ(I,NY,NX)=solutevar%KSO4_1e_conc
    CH0PQ(I,NY,NX)=solutevar%H0PO4_conc
    CH1PQ(I,NY,NX)=solutevar%H1PO4_2e_conc
    CPOQ(I,NY,NX)=solutevar%H2PO4_1e_conc
    CH3PQ(I,NY,NX)=solutevar%H3PO4_conc
    CF1PQ(I,NY,NX)=solutevar%FeHPO4_conc
    CF2PQ(I,NY,NX)=solutevar%FeH2PO4_conc
    CC0PQ(I,NY,NX)=solutevar%CaPO4_1e_con
    CC1PQ(I,NY,NX)=solutevar%CaHPO4_conc
    CC2PQ(I,NY,NX)=solutevar%CaH2PO4_1p_conc
    CM1PQ(I,NY,NX)=solutevar%MgHPO4_conc
    CSTRQ(I,NY,NX)=solutevar%CSTR1
!
!     SOLUTE CONCENTRATIONS IN SOIL
! for the POM complex, on the first day in the first year
! U means surface irrigation
  ELSEIF(K.EQ.micpar%k_POM.AND.I.EQ.1.AND.(.not.is_restart()).AND.is_first_year)THEN
    CCOU=solutevar%H2CO3_aqu_conc
    CCHU=solutevar%CH4_aqu_conc
    COXU=0._r8
    CNNU=solutevar%N2_aqu_conc
    CN2U=solutevar%N2O_aqu_conc
    trcn_irrig(ids_NH4,L,NY,NX)=solutevar%NH4_1p_conc
    trcn_irrig(idg_NH3,L,NY,NX)=solutevar%NH3_aqu_conc
    trcn_irrig(ids_NO3,L,NY,NX)=solutevar%CNOX
    trcn_irrig(ids_NH4B,L,NY,NX)=trcn_irrig(ids_NH4,L,NY,NX)
    trcn_irrig(idg_NH3B,L,NY,NX)=trcn_irrig(idg_NH3,L,NY,NX)
    trcn_irrig(ids_NO3B,L,NY,NX)=trcn_irrig(ids_NO3,L,NY,NX)

    CALU(L,NY,NX)=solutevar%Al_3p_conc
    CFEU(L,NY,NX)=solutevar%Fe_3p_conc
    CHYU(L,NY,NX)=solutevar%H_1p_conc
    CCAU(L,NY,NX)=solutevar%Ca_2p_conc
    CMGU(L,NY,NX)=solutevar%Mg_2p_conc
    CNAU(L,NY,NX)=solutevar%Na_1p_conc
    CKAU(L,NY,NX)=solutevar%K_1p_conc
    COHU(L,NY,NX)=solutevar%OH_1e_conc
    CSOU(L,NY,NX)=solutevar%SO4_2e_conc
    CCLU(L,NY,NX)=solutevar%Cl_e_conc
    CC3U(L,NY,NX)=solutevar%CO3_2e_conc
    CHCU(L,NY,NX)=solutevar%HCO3_e_conc
    CAL1U(L,NY,NX)=solutevar%AlOH_2p_conc
    CAL2U(L,NY,NX)=solutevar%AlOH2_p_conc
    CAL3U(L,NY,NX)=solutevar%AlOH3_conc
    CAL4U(L,NY,NX)=solutevar%AlOH4_1e_conc
    CALSU(L,NY,NX)=solutevar%AlSO4_1p_conc
    CFE1U(L,NY,NX)=solutevar%FeOH_2p_conc
    CFE2U(L,NY,NX)=solutevar%FeO2H2_p_conc
    CFE3U(L,NY,NX)=solutevar%FeO3H3_conc
    CFE4U(L,NY,NX)=solutevar%FeO4H4_1e_conc
    CFESU(L,NY,NX)=solutevar%FeSO4_1p_conc
    CCAOU(L,NY,NX)=solutevar%CaO2H2_conc
    CCACU(L,NY,NX)=solutevar%CaCO3_conc
    CCAHU(L,NY,NX)=solutevar%CaHCO3_1p_conc
    CCASU(L,NY,NX)=solutevar%CaSO4_conc
    CMGOU(L,NY,NX)=solutevar%MgOH_1p_conc
    CMGCU(L,NY,NX)=solutevar%MgCO3_conc
    CMGHU(L,NY,NX)=solutevar%MgHCO3_1p_conc
    CMGSU(L,NY,NX)=solutevar%MgSO4_conc
    CNACU(L,NY,NX)=solutevar%NaCO3_1e_conc
    CNASU(L,NY,NX)=solutevar%NaSO4_1e_conc
    CKASU(L,NY,NX)=solutevar%KSO4_1e_conc
    CH0PU(L,NY,NX)=solutevar%H0PO4_conc
    trcn_irrig(ids_H1PO4,L,NY,NX)=solutevar%H1PO4_2e_conc
    trcn_irrig(ids_H2PO4,L,NY,NX)=solutevar%H2PO4_1e_conc
    trcn_irrig(ids_H1PO4B,L,NY,NX)=trcn_irrig(ids_H1PO4,L,NY,NX)
    trcn_irrig(ids_H2PO4B,L,NY,NX)=trcn_irrig(ids_H2PO4,L,NY,NX)
    CH3PU(L,NY,NX)=solutevar%H3PO4_conc
    CF1PU(L,NY,NX)=solutevar%FeHPO4_conc
    CF2PU(L,NY,NX)=solutevar%FeH2PO4_conc
    CC0PU(L,NY,NX)=solutevar%CaPO4_1e_con
    CC1PU(L,NY,NX)=solutevar%CaHPO4_conc
    CC2PU(L,NY,NX)=solutevar%CaH2PO4_1p_conc
    CM1PU(L,NY,NX)=solutevar%MgHPO4_conc
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
    trc_gasml(idg_CO2,L,NY,NX)=CCO2EI(NY,NX)*VLsoiAirP(L,NY,NX)
    trc_gasml(idg_CH4,L,NY,NX)=AtmGgms(idg_CH4,NY,NX)*VLsoiAirP(L,NY,NX)
    trc_gasml(idg_O2,L,NY,NX)=AtmGgms(idg_O2,NY,NX)*VLsoiAirP(L,NY,NX)
    trc_gasml(idg_N2,L,NY,NX)=AtmGgms(idg_N2,NY,NX)*VLsoiAirP(L,NY,NX)
    trc_gasml(idg_N2O,L,NY,NX)=AtmGgms(idg_N2O,NY,NX)*VLsoiAirP(L,NY,NX)
    trc_gasml(idg_NH3,L,NY,NX)=AtmGgms(idg_NH3,NY,NX)*VLsoiAirP(L,NY,NX)
    trc_gasml(idg_H2,L,NY,NX)=AtmGgms(idg_H2,NY,NX)*VLsoiAirP(L,NY,NX)

!   ExtWaterTablet0: external water table depth
    IF(CumDepth2LayerBottom(L-1,NY,NX).LT.ExtWaterTablet0(NY,NX))THEN
      ! above water table
      trc_solml(idg_O2,L,NY,NX)=AtmGgms(idg_O2,NY,NX)*gas_solubility(idg_O2, ATCA(NY,NX)) &
        /(EXP(ACTCG(idg_O2)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP(L,NY,NX)
    ELSE
      !below water table
      trc_solml(idg_O2,L,NY,NX)=0._r8
    ENDIF

    trc_solml(idg_CO2,L,NY,NX)=CCO2EI(NY,NX)*gas_solubility(idg_CO2, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_CO2)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP(L,NY,NX)
    trc_solml(idg_CH4,L,NY,NX)=AtmGgms(idg_CH4,NY,NX)*gas_solubility(idg_CH4, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_CH4)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP(L,NY,NX)
    trc_solml(idg_N2,L,NY,NX)=AtmGgms(idg_N2,NY,NX)*gas_solubility(idg_N2, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_N2)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP(L,NY,NX)
    trc_solml(idg_N2O,L,NY,NX)=AtmGgms(idg_N2O,NY,NX)*gas_solubility(idg_N2O, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_N2O)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP(L,NY,NX)
    trc_solml(idg_H2,L,NY,NX)=AtmGgms(idg_H2,NY,NX)*gas_solubility(idg_H2, ATCA(NY,NX)) &
      /(EXP(ACTCG(idg_H2)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP(L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR MINERAL N AND P IN SOIL
!
    trc_solml(ids_NH4,L,NY,NX)=trcn_irrig(ids_NH4,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    trc_solml(ids_NH4B,L,NY,NX)=trcn_irrig(ids_NH4B,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    trc_solml(idg_NH3,L,NY,NX)=trcn_irrig(idg_NH3,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    trc_solml(idg_NH3B,L,NY,NX)=trcn_irrig(idg_NH3B,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    trc_solml(ids_NO3,L,NY,NX)=trcn_irrig(ids_NO3,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)*natomw
    trc_solml(ids_NO3B,L,NY,NX)=trcn_irrig(ids_NO3B,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)*natomw

    trc_solml(ids_H2PO4,L,NY,NX)=trcn_irrig(ids_H2PO4,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    trc_solml(ids_H2PO4B,L,NY,NX)=trcn_irrig(ids_H2PO4B,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw
    trc_solml(ids_H1PO4,L,NY,NX)=trcn_irrig(ids_H1PO4,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    trc_solml(ids_H1PO4B,L,NY,NX)=trcn_irrig(ids_H1PO4B,L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw

    trc_solml(ids_NO2,L,NY,NX)=0._r8
    trc_solml(ids_NO2B,L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS, ANIONS AND ION PAIRS IN SOIL
!
    trcSalt_solml(idsalt_Al,L,NY,NX)=CALU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_Fe,L,NY,NX)=CFEU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_Hp,L,NY,NX)=CHYU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_Ca,L,NY,NX)=CCAU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_Mg,L,NY,NX)=CMGU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_Na,L,NY,NX)=CNAU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_K,L,NY,NX)=CKAU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_OH,L,NY,NX)=COHU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_SO4,L,NY,NX)=CSOU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_Cl,L,NY,NX)=CCLU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_CO3,L,NY,NX)=CC3U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_HCO3,L,NY,NX)=CHCU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_AlOH,L,NY,NX)=CAL1U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_AlOH2,L,NY,NX)=CAL2U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_AlOH3,L,NY,NX)=CAL3U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_AlOH4,L,NY,NX)=CAL4U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_AlSO4,L,NY,NX)=CALSU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_FeOH,L,NY,NX)=CFE1U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_FeOH2,L,NY,NX)=CFE2U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_FeOH3,L,NY,NX)=CFE3U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_FeOH4,L,NY,NX)=CFE4U(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_FeSO4,L,NY,NX)=CFESU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_CaOH2,L,NY,NX)=CCAOU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_CaCO3,L,NY,NX)=CCACU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_CaHCO3,L,NY,NX)=CCAHU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_CaSO4,L,NY,NX)=CCASU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_MgOH2,L,NY,NX)=CMGOU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_MgCO3,L,NY,NX)=CMGCU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_MgHCO3,L,NY,NX)=CMGHU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_MgSO4,L,NY,NX)=CMGSU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_NaCO3,L,NY,NX)=CNACU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_NaSO4,L,NY,NX)=CNASU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_KSO4,L,NY,NX)=CKASU(L,NY,NX)*VLWatMicP(L,NY,NX)
    trcSalt_solml(idsalt_H0PO4,L,NY,NX)=CH0PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_H3PO4,L,NY,NX)=CH3PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_FeHPO4,L,NY,NX)=CF1PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_FeH2PO4,L,NY,NX)=CF2PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_CaPO4,L,NY,NX)=CC0PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_CaHPO4,L,NY,NX)=CC1PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_CaH2PO4,L,NY,NX)=CC2PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_MgHPO4,L,NY,NX)=CM1PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcSalt_solml(idsalt_H0PO4B,L,NY,NX)=CH0PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcSalt_solml(idsalt_H3PO4B,L,NY,NX)=CH3PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcSalt_solml(idsalt_FeHPO4B,L,NY,NX)=CF1PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcSalt_solml(idsalt_FeH2PO4B,L,NY,NX)=CF2PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcSalt_solml(idsalt_CaPO4B,L,NY,NX)=CC0PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcSalt_solml(idsalt_CaHPO4B,L,NY,NX)=CC1PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcSalt_solml(idsalt_CaH2PO4B,L,NY,NX)=CC2PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcSalt_solml(idsalt_MgHPO4B,L,NY,NX)=CM1PU(L,NY,NX)*VLWatMicP(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR ALL MATERIAL IN SOIL MACROPORES
!
    trc_soHml(ids_beg:ids_end,L,NY,NX)=0._r8
    trcSalt_soHml(idsalt_beg:idsaltb_end,L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR EXCHANGEABLE CATIONS AND ANIONS
!
    trcx_solml(idx_NH4,L,NY,NX)=solutevar%XNH4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)
    trcx_solml(idx_NH4B,L,NY,NX)=solutevar%XNH4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)
    trcx_solml(idx_Hp,L,NY,NX)=solutevar%XHY1*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_Al,L,NY,NX)=solutevar%XAl_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_Fe,L,NY,NX)=solutevar%XFe_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_Ca,L,NY,NX)=solutevar%XCa_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_Mg,L,NY,NX)=solutevar%XMg_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_Na,L,NY,NX)=solutevar%XNa_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_K,L,NY,NX)=solutevar%XK_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_COOH,L,NY,NX)=solutevar%XHC1*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_AlOH2,L,NY,NX)=solutevar%XAlO2H2_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_FeOH2,L,NY,NX)=solutevar%XFeO2H2_conc*SoilMicPMassLayer(L,NY,NX)
    trcx_solml(idx_OHe,L,NY,NX)=solutevar%XOH_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_OH,L,NY,NX)=solutevar%XROH1_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_OHp,L,NY,NX)=solutevar%XROH2_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_HPO4,L,NY,NX)=solutevar%XHPO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_H2PO4,L,NY,NX)=solutevar%XH2PO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcx_solml(idx_OHeB,L,NY,NX)=solutevar%XOH_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_OHB,L,NY,NX)=solutevar%XROH1_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_OHpB,L,NY,NX)=solutevar%XROH2_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_HPO4B,L,NY,NX)=solutevar%XHPO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcx_solml(idx_H2PO4B,L,NY,NX)=solutevar%XH2PO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR PRECIPITATES
!
    trcp_salml(idsp_AlOH3,L,NY,NX)=solutevar%Precp_AlO3H3_conc*SoilMicPMassLayer(L,NY,NX)
    trcp_salml(idsp_FeOH3,L,NY,NX)=solutevar%Precp_FeO3H3_conc*SoilMicPMassLayer(L,NY,NX)
    trcp_salml(idsp_CaCO3,L,NY,NX)=solutevar%Precp_CaCO3_conc*SoilMicPMassLayer(L,NY,NX)
    trcp_salml(idsp_CaSO4,L,NY,NX)=solutevar%Precp_CaSO4_conc*SoilMicPMassLayer(L,NY,NX)
    trcp_salml(idsp_AlPO4,L,NY,NX)=solutevar%Precp_AlPO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_FePO4,L,NY,NX)=solutevar%Precp_FePO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_CaHPO4,L,NY,NX)=solutevar%Precp_CaHPO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_HA,L,NY,NX)=solutevar%Precp_Ca5P3O12O3H3_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)
    trcp_salml(idsp_CaH2PO4,L,NY,NX)=0._r8
    trcp_salml(idsp_AlPO4B,L,NY,NX)=solutevar%Precp_AlPO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcp_salml(idsp_FePO4B,L,NY,NX)=solutevar%Precp_FePO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcp_salml(idsp_CaHPO4B,L,NY,NX)=solutevar%Precp_CaHPO4_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
    trcp_salml(idsp_HAB,L,NY,NX)=solutevar%Precp_Ca5P3O12O3H3_conc*SoilMicPMassLayer(L,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)
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
  IF(.not.is_restart().AND.is_first_year)THEN
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
      IF(VLHeatCapSnow(L,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
        VOLWW=VLWatSnow(L,NY,NX)+VLDrySnoWE(L,NY,NX)+VLIceSnow(L,NY,NX)*DENSICE
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
          trcs_solsml(idsalt_Al,L,NY,NX)=VOLWW*CALR(NY,NX)
          trcs_solsml(idsalt_Fe,L,NY,NX)=VOLWW*CFER(NY,NX)
          trcs_solsml(idsalt_Hp,L,NY,NX)=VOLWW*CHYR(NY,NX)
          trcs_solsml(idsalt_Ca,L,NY,NX)=VOLWW*CCAR(NY,NX)
          trcs_solsml(idsalt_Mg,L,NY,NX)=VOLWW*CMGR(NY,NX)
          trcs_solsml(idsalt_Na,L,NY,NX)=VOLWW*CNAR(NY,NX)
          trcs_solsml(idsalt_K,L,NY,NX)=VOLWW*CKAR(NY,NX)
          trcs_solsml(idsalt_OH,L,NY,NX)=VOLWW*COHR(NY,NX)
          trcs_solsml(idsalt_SO4,L,NY,NX)=VOLWW*CSOR(NY,NX)
          trcs_solsml(idsalt_Cl,L,NY,NX)=VOLWW*CCLR(NY,NX)
          trcs_solsml(idsalt_CO3,L,NY,NX)=VOLWW*CC3R(NY,NX)
          trcs_solsml(idsalt_HCO3,L,NY,NX)=VOLWW*CHCR(NY,NX)
          trcs_solsml(idsalt_AlOH,L,NY,NX)=VOLWW*CAL1R(NY,NX)
          trcs_solsml(idsalt_AlOH2,L,NY,NX)=VOLWW*CAL2R(NY,NX)
          trcs_solsml(idsalt_AlOH3,L,NY,NX)=VOLWW*CAL3R(NY,NX)
          trcs_solsml(idsalt_AlOH4,L,NY,NX)=VOLWW*CAL4R(NY,NX)
          trcs_solsml(idsalt_AlSO4,L,NY,NX)=VOLWW*CALSR(NY,NX)
          trcs_solsml(idsalt_FeOH,L,NY,NX)=VOLWW*CFE1R(NY,NX)
          trcs_solsml(idsalt_FeOH2,L,NY,NX)=VOLWW*CFE2R(NY,NX)
          trcs_solsml(idsalt_FeOH3,L,NY,NX)=VOLWW*CFE3R(NY,NX)
          trcs_solsml(idsalt_FeOH4,L,NY,NX)=VOLWW*CFE4R(NY,NX)
          trcs_solsml(idsalt_FeSO4,L,NY,NX)=VOLWW*CFESR(NY,NX)
          trcs_solsml(idsalt_CaOH2,L,NY,NX)=VOLWW*CCAOR(NY,NX)
          trcs_solsml(idsalt_CaCO3,L,NY,NX)=VOLWW*CCACR(NY,NX)
          trcs_solsml(idsalt_CaHCO3,L,NY,NX)=VOLWW*CCAHR(NY,NX)
          trcs_solsml(idsalt_CaSO4,L,NY,NX)=VOLWW*CCASR(NY,NX)
          trcs_solsml(idsalt_MgOH2,L,NY,NX)=VOLWW*CMGOR(NY,NX)
          trcs_solsml(idsalt_MgCO3,L,NY,NX)=VOLWW*CMGCR(NY,NX)
          trcs_solsml(idsalt_MgHCO3,L,NY,NX)=VOLWW*CMGHR(NY,NX)
          trcs_solsml(idsalt_MgSO4,L,NY,NX)=VOLWW*CMGSR(NY,NX)
          trcs_solsml(idsalt_NaCO3,L,NY,NX)=VOLWW*CNACR(NY,NX)
          trcs_solsml(idsalt_NaSO4,L,NY,NX)=VOLWW*CNASR(NY,NX)
          trcs_solsml(idsalt_KSO4,L,NY,NX)=VOLWW*CKASR(NY,NX)
          trcs_solsml(idsalt_H0PO4,L,NY,NX)=VOLWW*CH0PR(NY,NX)
          trcs_solsml(idsalt_H3PO4,L,NY,NX)=VOLWW*CH3PR(NY,NX)
          trcs_solsml(idsalt_FeHPO4,L,NY,NX)=VOLWW*CF1PR(NY,NX)
          trcs_solsml(idsalt_FeH2PO4,L,NY,NX)=VOLWW*CF2PR(NY,NX)
          trcs_solsml(idsalt_CaPO4,L,NY,NX)=VOLWW*CC0PR(NY,NX)
          trcs_solsml(idsalt_CaHPO4,L,NY,NX)=VOLWW*CC1PR(NY,NX)
          trcs_solsml(idsalt_CaH2PO4,L,NY,NX)=VOLWW*CC2PR(NY,NX)
          trcs_solsml(idsalt_MgHPO4,L,NY,NX)=VOLWW*CM1PR(NY,NX)
        ENDIF
      ELSE
        trcg_solsml(idg_beg:idg_end-1,L,NY,NX)=0._r8
        trcn_solsml(ids_nut_beg:ids_nuts_end,L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS AND ANIONS IN SNOWPACK
!
        IF(salt_model)THEN
          trcs_solsml(idsalt_beg:idsalt_end,L,NY,NX)=0._r8
        ENDIF
      ENDIF
    ENDDO D9985
  ENDIF
  end subroutine InitialState

end module StarteMod
