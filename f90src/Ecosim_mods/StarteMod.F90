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
  integer :: NY,NX,L,I,K
!     begin_execution
!
!     INITIALIZE CATION AND ANION CONCENTRATIONS
!     IN PRECIPITATION (K=1), IRRIGATION (K=2) AND SOIL (K=3)
!     FROM WEATHER, IRRIGATION AND SOIL FILES IN 'READS'
!
  DO   NX=NHW,NHE
    DO  NY=NVN,NVS
      solutevar%CCO2M = CCO2EI(NY,NX)/catomw
      solutevar%CCH4M = AtmGasCgperm3(idg_CH4,NY,NX)/catomw
      solutevar%COXYM = AtmGasCgperm3(idg_O2,NY,NX)/32.0_r8
      solutevar%CZ2GM = AtmGasCgperm3(idg_N2,NY,NX)/natomw
      solutevar%CZ2OM = AtmGasCgperm3(idg_N2O,NY,NX)/natomw
      solutevar%ATCA  = ATCA(NY,NX)
      solutevar%ZEROS = ZEROS(NY,NX)

      DO I=1,366
        DO  L=NU(NY,NX),NL(NY,NX)
          D2000: DO K=micpar%k_fine_litr,micpar%k_POM
            BulkSoilMass=0._r8
!
            IF(K.EQ.micpar%k_manure.AND.L.EQ.1)THEN
!     INITIALIZE IRRIGATION WATER
              solutevar%CN4Z = NH4_irrig_mole_conc(I,NY,NX)
              solutevar%CNOZ = NO3_irrig_mole_conc(I,NY,NX)
              solutevar%CPOZ = H2PO4_irrig_mole_conc(I,NY,NX)
              if(salt_model)then
                solutevar%Al_mole_conc = trcsalt_irrig_mole_conc_col(idsalt_Al,I,NY,NX)
                solutevar%CFEZ = trcsalt_irrig_mole_conc_col(idsalt_Fe,I,NY,NX)
                solutevar%CCAZ = trcsalt_irrig_mole_conc_col(idsalt_Ca,I,NY,NX)
                solutevar%CMGZ = trcsalt_irrig_mole_conc_col(idsalt_Mg,I,NY,NX)
                solutevar%CNAZ = trcsalt_irrig_mole_conc_col(idsalt_Na,I,NY,NX)
                solutevar%CKAZ = trcsalt_irrig_mole_conc_col(idsalt_K,I,NY,NX)
                solutevar%CSOZ = trcsalt_irrig_mole_conc_col(idsalt_SO4,I,NY,NX)
                solutevar%CCLZ = trcsalt_irrig_mole_conc_col(idsalt_Cl,I,NY,NX)
              endif

              solutevar%H_1p_aqua_mole_conc  = 10.0_r8**(-(PHQ(I,NY,NX)-3.0_r8))
              solutevar%OH_1e_aqua_mole_conc = DPH2O/solutevar%H_1p_aqua_mole_conc
            ELSE
              IF(I.EQ.1)then
                IF(K.EQ.micpar%k_fine_litr.AND.L.EQ.1)THEN
!     INITIALIZE RAINFALL, top layer
                  solutevar%H_1p_aqua_mole_conc  = 10.0_r8**(-(PHR(NY,NX)-3.0_r8))
                  solutevar%OH_1e_aqua_mole_conc = DPH2O/solutevar%H_1p_aqua_mole_conc
                  solutevar%CN4Z       = NH4_rain_mole_conc(NY,NX)
                  solutevar%CNOZ       = NO3_rain_mole_conc(NY,NX)
                  solutevar%CPOZ       = H2PO4_rain_mole_conc(NY,NX)
                  if(salt_model)then
                    solutevar%Al_mole_conc = trcsalt_rain_mole_conc_col(idsalt_Al,NY,NX)
                    solutevar%CFEZ = trcsalt_rain_mole_conc_col(idsalt_Fe,NY,NX)
                    solutevar%CCAZ = trcsalt_rain_mole_conc_col(idsalt_Ca,NY,NX)
                    solutevar%CMGZ = trcsalt_rain_mole_conc_col(idsalt_Mg,NY,NX)
                    solutevar%CNAZ = trcsalt_rain_mole_conc_col(idsalt_Na,NY,NX)
                    solutevar%CKAZ = trcsalt_rain_mole_conc_col(idsalt_K,NY,NX)
                    solutevar%CSOZ = trcsalt_rain_mole_conc_col(idsalt_SO4,NY,NX)
                    solutevar%CCLZ = trcsalt_rain_mole_conc_col(idsalt_Cl,NY,NX)
                  endif
!               ELSEIF(K.EQ.micpar%k_POM.AND.(.not.is_restart()).AND.is_first_year)THEN
!
!     INITIALIZE SOIL WATER

                  IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
                    BulkSoilMass=VLSoilMicPMass_vr(L,NY,NX)
                  ELSE
                    !this is for water
                    BulkSoilMass=VLMicP_vr(L,NY,NX)
                  ENDIF
                  trcx_solml_vr(idx_CEC,L,NY,NX) = AMAX1(CNH4_vr(L,NY,NX),CEC_vr(L,NY,NX))*BulkSoilMass
                  trcx_solml_vr(idx_AEC,L,NY,NX) = AMAX1(CPO4_vr(L,NY,NX),AEC_vr(L,NY,NX))*BulkSoilMass
                  solutevar%CN4X                 = CNH4_vr(L,NY,NX)
                  solutevar%CALX                 = CAL_vr(L,NY,NX)
                  solutevar%CFEX                 = CFE_vr(L,NY,NX)
                  solutevar%CaX_conc             = CCA_vr(L,NY,NX)
                  solutevar%MgX_conc             = CMG_vr(L,NY,NX)
                  solutevar%CNAX                 = CNA_vr(L,NY,NX)
                  solutevar%CKAX                 = CKA_vr(L,NY,NX)
                  solutevar%CSOX                 = CSO4_vr(L,NY,NX)
                  solutevar%CCLX                 = CCL_vr(L,NY,NX)
                  solutevar%CALPOX               = CALPO_vr(L,NY,NX)
                  solutevar%CFEPOX               = CFEPO_vr(L,NY,NX)
                  solutevar%CCAPDX               = CCAPD_vr(L,NY,NX)
                  solutevar%CCAPHX               = CCAPH_vr(L,NY,NX)
                  solutevar%CALOHX               = CALOH_vr(L,NY,NX)
                  solutevar%CFEOHX               = CFEOH_vr(L,NY,NX)
                  solutevar%CCACOX               = CCACO_vr(L,NY,NX)
                  solutevar%CCASOX               = CCASO_vr(L,NY,NX)
                  solutevar%H_1p_aqua_mole_conc  = 10.0_r8**(-(PH_vr(L,NY,NX)-3.0_r8))
                  solutevar%CNOZ                 = CNO3_vr(L,NY,NX)
                  solutevar%CPOZ                 = CPO4_vr(L,NY,NX)
                  solutevar%OH_1e_aqua_mole_conc = DPH2O/solutevar%H_1p_aqua_mole_conc
                  solutevar%CN4Z                 = solutevar%CN4X
                  solutevar%Al_mole_conc                 = solutevar%CALX
                  solutevar%CFEZ                 = solutevar%CFEX
                  solutevar%CCAZ                 = solutevar%CaX_conc
                  solutevar%CMGZ                 = solutevar%MgX_conc
                  solutevar%CNAZ                 = solutevar%CNAX
                  solutevar%CKAZ                 = solutevar%CKAX
                  solutevar%CSOZ                 = solutevar%CSOX
                  solutevar%CCLZ                 = solutevar%CCLX
                ELSE
                  cycle
                ENDIF
              ELSE
                cycle
              ENDIF
            ENDIF

            solutevar%XAEC  = trcx_solml_vr(idx_AEC,L,NY,NX)
            solutevar%CEC   = CEC_vr(L,NY,NX)
            solutevar%ORGC  = SoilOrgM_vr(ielmc,L,NY,NX)
            solutevar%VLPO4 = trcs_VLN_vr(ids_H1PO4,L,NY,NX)
            solutevar%XCEC  = trcx_solml_vr(idx_CEC,L,NY,NX)
            solutevar%GKC4  = GKC4_vr(L,NY,NX)
            solutevar%GKCA  = GKCA_vr(L,NY,NX)
            solutevar%GKCH  = GKCH_vr(L,NY,NX)
            solutevar%GKCK  = GKCK_vr(L,NY,NX)
            solutevar%GKCN  = GKCN_vr(L,NY,NX)
            solutevar%GKCM  = GKCM_vr(L,NY,NX)
            solutevar%VLWatMicP = VLWatMicP_vr(L,NY,NX)

            call InitSoluteModel(K,BulkSoilMass,solutevar)

            call InitSoluteConcs(K,I,L,NY,NX,solutevar)
          ENDDO D2000
        enddo
      ENDDO

      call InitialState(NY,NX)

    ENDDO
  ENDDO

  END subroutine starte
!------------------------------------------------------------------------------------------

  subroutine InitSoluteConcs(K,I,L,NY,NX,solutevar)

  implicit none
  integer, intent(in) :: K, I, L, NY, NX
  type(solutedtype), intent(in) :: solutevar
  integer :: nsalts,ids,idg
!
!     SOLUTE CONCENTRATIONS IN PRECIPITATION
!
  associate(                       &
    iprotein  => micpar%iprotein,  &
    k_manure  => micpar%k_manure   &
  )
  IF(K.EQ.micpar%k_fine_litr.AND.L.EQ.1.AND.I.EQ.1)THEN
  !litter pool, top soil layer
    trcg_rain_mole_conc_col(idg_CO2,NY,NX) = solutevar%H2CO3_aqua_mole_conc
    trcg_rain_mole_conc_col(idg_CH4,NY,NX) = solutevar%CH4_aqua_mole_conc
    trcg_rain_mole_conc_col(idg_O2,NY,NX)  = solutevar%O2_aqua_mole_conc
    trcg_rain_mole_conc_col(idg_N2,NY,NX)  = solutevar%N2_aqua_mole_conc
    trcg_rain_mole_conc_col(idg_N2O,NY,NX) = solutevar%N2O_aqua_mole_conc
    trcg_rain_mole_conc_col(idg_NH3,NY,NX) = solutevar%NH3_aqua_mole_conc
    NH4_rain_mole_conc(NY,NX)              = solutevar%NH4_1p_aqua_mole_conc    
    if(salt_model)then
      trcsalt_rain_mole_conc_col(idsalt_Al,NY,NX)       = solutevar%Al_3p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_Fe,NY,NX)       = solutevar%Fe_3p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_Hp,NY,NX)       = solutevar%H_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_Ca,NY,NX)       = solutevar%Ca_2p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_Mg,NY,NX)       = solutevar%Mg_2p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_Na,NY,NX)       = solutevar%Na_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_K,NY,NX)        = solutevar%K_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_OH,NY,NX)       = solutevar%OH_1e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_SO4,NY,NX)      = solutevar%SO4_2e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_Cl,NY,NX)       = solutevar%Cl_e_conc
      trcsalt_rain_mole_conc_col(idsalt_CO3,NY,NX)      = solutevar%CO3_2e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_HCO3,NY,NX)     = solutevar%HCO3_e_conc
      trcsalt_rain_mole_conc_col(idsalt_AlOH,NY,NX)     = solutevar%AlOH_2p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_AlOH2,NY,NX)    = solutevar%AlO2H2_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_AlOH3,NY,NX)    = solutevar%AlO3H3_conc
      trcsalt_rain_mole_conc_col(idsalt_AlOH4,NY,NX)    = solutevar%AlO4H4_1e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_AlSO4,NY,NX)    = solutevar%AlSO4_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_FeOH,NY,NX)     = solutevar%FeOH_2p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_FeOH2,NY,NX)    = solutevar%FeO2H2_p_conc
      trcsalt_rain_mole_conc_col(idsalt_FeOH3,NY,NX)    = solutevar%FeO3H3_conc
      trcsalt_rain_mole_conc_col(idsalt_FeOH4,NY,NX)    = solutevar%FeO4H4_1e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_FeSO4,NY,NX)    = solutevar%FeSO4_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_CaOH,NY,NX)     = solutevar%CaO2H2_conc
      trcsalt_rain_mole_conc_col(idsalt_CaCO3,NY,NX)    = solutevar%CaCO3_conc
      trcsalt_rain_mole_conc_col(idsalt_CaHCO3,NY,NX)   = solutevar%CaHCO3_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_CaSO4,NY,NX)    = solutevar%CaSO4_conc
      trcsalt_rain_mole_conc_col(idsalt_MgOH2,NY,NX)    = solutevar%MgOH_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_MgCO3,NY,NX)    = solutevar%MgCO3_conc
      trcsalt_rain_mole_conc_col(idsalt_MgHCO3,NY,NX)   = solutevar%MgHCO3_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_MgSO4,NY,NX)    = solutevar%MgSO4_conc
      trcsalt_rain_mole_conc_col(idsalt_NaCO3,NY,NX)    = solutevar%NaCO3_1e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_NaSO4,NY,NX)    = solutevar%NaSO4_1e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_KSO4,NY,NX)     = solutevar%KSO4_1e_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_H0PO4,NY,NX)    = solutevar%H0PO4_3e_conc
      trcsalt_rain_mole_conc_col(idsalt_H3PO4,NY,NX)    = solutevar%H3PO4_conc
      trcsalt_rain_mole_conc_col(idsalt_FeHPO4,NY,NX)   = solutevar%FeHPO4_p_conc
      trcsalt_rain_mole_conc_col(idsalt_FeH2PO4,NY,NX)  = solutevar%FeH2PO4_2p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_CaPO4,NY,NX)    = solutevar%CaPO4_1e_con
      trcsalt_rain_mole_conc_col(idsalt_CaHPO4,NY,NX)   = solutevar%CaHPO4_conc
      trcsalt_rain_mole_conc_col(idsalt_CaH4P2O8,NY,NX) = solutevar%CaH4P2O8_1p_aqua_mole_conc
      trcsalt_rain_mole_conc_col(idsalt_MgHPO4,NY,NX)   = solutevar%MgHPO4_conc
    endif
    HPO4_rain_mole_conc(NY,NX)      = solutevar%H1PO4_2e_aqua_mole_conc
    H2PO4_rain_mole_conc(NY,NX)     = solutevar%H2PO4_1e_aqua_mole_conc
    SurfIrrig_IonStrenth_col(NY,NX) = solutevar%CSTR1
!
!     SOLUTE CONCENTRATIONS IN IRRIGATION
!
  ELSEIF(K.EQ.micpar%k_manure .AND. L.EQ.1)THEN
  ! manure, top layer
    trcg_irrig_mole_conc_col(idg_CO2,NY,NX) = solutevar%H2CO3_aqua_mole_conc
    trcg_irrig_mole_conc_col(idg_CH4,NY,NX) = solutevar%CH4_aqua_mole_conc
    trcg_irrig_mole_conc_col(idg_O2,NY,NX)  = solutevar%O2_aqua_mole_conc
    trcg_irrig_mole_conc_col(idg_N2,NY,NX)  = solutevar%N2_aqua_mole_conc
    trcg_irrig_mole_conc_col(idg_N2O,NY,NX) = solutevar%N2O_aqua_mole_conc
    trcg_irrig_mole_conc_col(idg_NH3,NY,NX) = solutevar%NH3_aqua_mole_conc
    NH4_irrig_mole_conc(I,NY,NX)            = solutevar%NH4_1p_aqua_mole_conc
    HPO4_irrig_mole_conc(I,NY,NX)           = solutevar%H1PO4_2e_aqua_mole_conc
    H2PO4_irrig_mole_conc(I,NY,NX)          = solutevar%H2PO4_1e_aqua_mole_conc
    CSTRQ(I,NY,NX)                          = solutevar%CSTR1
    if(salt_model)then
      trcsalt_irrig_mole_conc_col(idsalt_Al,I,NY,NX)       = solutevar%Al_3p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_Fe,I,NY,NX)       = solutevar%Fe_3p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_Hp,I,NY,NX)       = solutevar%H_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_Ca,I,NY,NX)       = solutevar%Ca_2p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_Mg,I,NY,NX)       = solutevar%Mg_2p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_Na,I,NY,NX)       = solutevar%Na_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_K,I,NY,NX)        = solutevar%K_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_OH,I,NY,NX)       = solutevar%OH_1e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_SO4,I,NY,NX)      = solutevar%SO4_2e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_Cl,I,NY,NX)       = solutevar%Cl_e_conc
      trcsalt_irrig_mole_conc_col(idsalt_CO3,I,NY,NX)      = solutevar%CO3_2e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_HCO3,I,NY,NX)     = solutevar%HCO3_e_conc
      trcsalt_irrig_mole_conc_col(idsalt_AlOH,I,NY,NX)     = solutevar%AlOH_2p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_AlOH2,I,NY,NX)    = solutevar%AlO2H2_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_AlOH3,I,NY,NX)    = solutevar%AlO3H3_conc
      trcsalt_irrig_mole_conc_col(idsalt_AlOH4,I,NY,NX)    = solutevar%AlO4H4_1e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_AlSO4,I,NY,NX)    = solutevar%AlSO4_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_FeOH,I,NY,NX)     = solutevar%FeOH_2p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_FeOH2,I,NY,NX)    = solutevar%FeO2H2_p_conc
      trcsalt_irrig_mole_conc_col(idsalt_FeOH3,I,NY,NX)    = solutevar%FeO3H3_conc
      trcsalt_irrig_mole_conc_col(idsalt_FeOH4,I,NY,NX)    = solutevar%FeO4H4_1e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_FeSO4,I,NY,NX)    = solutevar%FeSO4_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_CaOH,I,NY,NX)     = solutevar%CaO2H2_conc
      trcsalt_irrig_mole_conc_col(idsalt_CaCO3,I,NY,NX)    = solutevar%CaCO3_conc
      trcsalt_irrig_mole_conc_col(idsalt_CaHCO3,I,NY,NX)   = solutevar%CaHCO3_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_CaSO4,I,NY,NX)    = solutevar%CaSO4_conc
      trcsalt_irrig_mole_conc_col(idsalt_MgOH2,I,NY,NX)    = solutevar%MgOH_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_MgCO3,I,NY,NX)    = solutevar%MgCO3_conc
      trcsalt_irrig_mole_conc_col(idsalt_MgHCO3,I,NY,NX)   = solutevar%MgHCO3_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_MgSO4,I,NY,NX)    = solutevar%MgSO4_conc
      trcsalt_irrig_mole_conc_col(idsalt_NaCO3,I,NY,NX)    = solutevar%NaCO3_1e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_NaSO4,I,NY,NX)    = solutevar%NaSO4_1e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_KSO4,I,NY,NX)     = solutevar%KSO4_1e_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_H0PO4,I,NY,NX)    = solutevar%H0PO4_3e_conc
      trcsalt_irrig_mole_conc_col(idsalt_H3PO4,I,NY,NX)    = solutevar%H3PO4_conc
      trcsalt_irrig_mole_conc_col(idsalt_FeHPO4,I,NY,NX)   = solutevar%FeHPO4_p_conc
      trcsalt_irrig_mole_conc_col(idsalt_FeH2PO4,I,NY,NX)  = solutevar%FeH2PO4_2p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_CaPO4,I,NY,NX)    = solutevar%CaPO4_1e_con
      trcsalt_irrig_mole_conc_col(idsalt_CaHPO4,I,NY,NX)   = solutevar%CaHPO4_conc
      trcsalt_irrig_mole_conc_col(idsalt_CaH4P2O8,I,NY,NX) = solutevar%CaH4P2O8_1p_aqua_mole_conc
      trcsalt_irrig_mole_conc_col(idsalt_MgHPO4,I,NY,NX)   = solutevar%MgHPO4_conc
    endif
!
!     SOLUTE CONCENTRATIONS IN SOIL
! for the POM complex, on the first day in the first year
! U means surface irrigation
! first year and not a restart run
  ELSEIF(K.EQ.micpar%k_POM.AND.I.EQ.1.AND.(.not.is_restart()).AND.is_first_year)THEN
    CCOU                           = solutevar%H2CO3_aqua_mole_conc*catomw
    CCHU                           = solutevar%CH4_aqua_mole_conc*catomw
    COXU                           = 0._r8
    CNNU                           = solutevar%N2_aqua_mole_conc
    CN2U                           = solutevar%N2O_aqua_mole_conc
    trcn_irrig_vr(ids_NH4,L,NY,NX)    = solutevar%NH4_1p_aqua_mole_conc
    trcn_irrig_vr(idg_NH3,L,NY,NX)    = solutevar%NH3_aqua_mole_conc
    trcn_irrig_vr(ids_NO3,L,NY,NX)    = solutevar%CNOX
    trcn_irrig_vr(ids_NH4B,L,NY,NX)   = trcn_irrig_vr(ids_NH4,L,NY,NX)
    trcn_irrig_vr(idg_NH3B,L,NY,NX)   = trcn_irrig_vr(idg_NH3,L,NY,NX)
    trcn_irrig_vr(ids_NO3B,L,NY,NX)   = trcn_irrig_vr(ids_NO3,L,NY,NX)
    trcn_irrig_vr(ids_H1PO4,L,NY,NX)  = solutevar%H1PO4_2e_aqua_mole_conc
    trcn_irrig_vr(ids_H2PO4,L,NY,NX)  = solutevar%H2PO4_1e_aqua_mole_conc
    trcn_irrig_vr(ids_H1PO4B,L,NY,NX) = trcn_irrig_vr(ids_H1PO4,L,NY,NX)
    trcn_irrig_vr(ids_H2PO4B,L,NY,NX) = trcn_irrig_vr(ids_H2PO4,L,NY,NX)

    if(salt_model)then
      trcSalt_irrig_vr(idsalt_Al,L,NY,NX)       = solutevar%Al_3p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_Fe,L,NY,NX)       = solutevar%Fe_3p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_Hp,L,NY,NX)       = solutevar%H_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_Ca,L,NY,NX)       = solutevar%Ca_2p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_Mg,L,NY,NX)       = solutevar%Mg_2p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_Na,L,NY,NX)       = solutevar%Na_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_K,L,NY,NX)        = solutevar%K_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_OH,L,NY,NX)       = solutevar%OH_1e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_SO4,L,NY,NX)      = solutevar%SO4_2e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_Cl,L,NY,NX)       = solutevar%Cl_e_conc
      trcSalt_irrig_vr(idsalt_CO3,L,NY,NX)      = solutevar%CO3_2e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_HCO3,L,NY,NX)     = solutevar%HCO3_e_conc
      trcSalt_irrig_vr(idsalt_AlOH,L,NY,NX)     = solutevar%AlOH_2p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_AlOH2,L,NY,NX)    = solutevar%AlO2H2_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_AlOH3,L,NY,NX)    = solutevar%AlO3H3_conc
      trcSalt_irrig_vr(idsalt_AlOH4,L,NY,NX)    = solutevar%AlO4H4_1e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_AlSO4,L,NY,NX)    = solutevar%AlSO4_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_FeOH,L,NY,NX)     = solutevar%FeOH_2p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_FeOH2,L,NY,NX)    = solutevar%FeO2H2_p_conc
      trcSalt_irrig_vr(idsalt_FeOH3,L,NY,NX)    = solutevar%FeO3H3_conc
      trcSalt_irrig_vr(idsalt_FeOH4,L,NY,NX)    = solutevar%FeO4H4_1e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_FeSO4,L,NY,NX)    = solutevar%FeSO4_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_CaOH,L,NY,NX)     = solutevar%CaO2H2_conc
      trcSalt_irrig_vr(idsalt_CaCO3,L,NY,NX)    = solutevar%CaCO3_conc
      trcSalt_irrig_vr(idsalt_CaHCO3,L,NY,NX)   = solutevar%CaHCO3_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_CaSO4,L,NY,NX)    = solutevar%CaSO4_conc
      trcSalt_irrig_vr(idsalt_MgOH2,L,NY,NX)    = solutevar%MgOH_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_MgCO3,L,NY,NX)    = solutevar%MgCO3_conc
      trcSalt_irrig_vr(idsalt_MgHCO3,L,NY,NX)   = solutevar%MgHCO3_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_MgSO4,L,NY,NX)    = solutevar%MgSO4_conc
      trcSalt_irrig_vr(idsalt_NaCO3,L,NY,NX)    = solutevar%NaCO3_1e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_NaSO4,L,NY,NX)    = solutevar%NaSO4_1e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_KSO4,L,NY,NX)     = solutevar%KSO4_1e_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_H0PO4,L,NY,NX)    = solutevar%H0PO4_3e_conc
      trcSalt_irrig_vr(idsalt_H3PO4,L,NY,NX)    = solutevar%H3PO4_conc
      trcSalt_irrig_vr(idsalt_FeHPO4,L,NY,NX)   = solutevar%FeHPO4_p_conc
      trcSalt_irrig_vr(idsalt_FeH2PO4,L,NY,NX)  = solutevar%FeH2PO4_2p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_CaPO4,L,NY,NX)    = solutevar%CaPO4_1e_con
      trcSalt_irrig_vr(idsalt_CaHPO4,L,NY,NX)   = solutevar%CaHPO4_conc
      trcSalt_irrig_vr(idsalt_CaH4P2O8,L,NY,NX) = solutevar%CaH4P2O8_1p_aqua_mole_conc
      trcSalt_irrig_vr(idsalt_MgHPO4,L,NY,NX)   = solutevar%MgHPO4_conc
    endif
!
!   INITIAL STATE VARIABLES FOR GAS IN SOIL
!   CCO2EI is set to the first year, because AtmGasCgperm3(idg_CO2,x)
!   varies year by year, while other tracer gases are fixed year by year,
!   this is not quite right for CH4, and N2O. However, the current implementation
!   make sure the inexact restart run works. When exact restart is used, trcg_gasml_vr
!   and trcs_solml_vr will be read from restart file, so that the following inconsistent
!   use between CO2 and other gas tracers can be avoided.
!   Comment by Jinyun Tang, Nov 11, 2022
!
    trcg_gasml_vr(idg_CO2,L,NY,NX) = CCO2EI(NY,NX)*VLsoiAirP_vr(L,NY,NX)
    DO idg=idg_beg,idg_NH3
      trcg_gasml_vr(idg,L,NY,NX) = AtmGasCgperm3(idg,NY,NX)*VLsoiAirP_vr(L,NY,NX)
    ENDDO

!   ExtWaterTablet0_col: external water table depth
    IF(CumDepz2LayBottom_vr(L-1,NY,NX).LT.ExtWaterTablet0_col(NY,NX))THEN
      ! above water table
      trcs_solml_vr(idg_O2,L,NY,NX)=AtmGasCgperm3(idg_O2,NY,NX)*gas_solubility(idg_O2, ATCA(NY,NX)) &
        /(EXP(GasSechenovConst(idg_O2)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP_vr(L,NY,NX)
    ELSE
      !below water table
      trcs_solml_vr(idg_O2,L,NY,NX)=0._r8
    ENDIF

    
    DO idg=idg_beg,idg_NH3-1  
      trcs_solml_vr(idg,L,NY,NX)=AtmGasCgperm3(idg,NY,NX)*gas_solubility(idg, ATCA(NY,NX)) &
        /(EXP(GasSechenovConst(idg)*solutevar%CSTR1))*solutevar%FH2O*VLWatMicP_vr(L,NY,NX)
    ENDDO  
!
!     INITIAL STATE VARIABLES FOR MINERAL N AND P IN SOIL
!
    trcs_solml_vr(ids_NH4,L,NY,NX)  = trcn_irrig_vr(ids_NH4,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
    trcs_solml_vr(ids_NH4B,L,NY,NX) = trcn_irrig_vr(ids_NH4B,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
    trcs_solml_vr(idg_NH3,L,NY,NX)  = trcn_irrig_vr(idg_NH3,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
    trcs_solml_vr(idg_NH3B,L,NY,NX) = trcn_irrig_vr(idg_NH3B,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
    trcs_solml_vr(ids_NO3,L,NY,NX)  = trcn_irrig_vr(ids_NO3,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)*natomw
    trcs_solml_vr(ids_NO3B,L,NY,NX) = trcn_irrig_vr(ids_NO3B,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)*natomw

    trcs_solml_vr(ids_H2PO4,L,NY,NX)  = trcn_irrig_vr(ids_H2PO4,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
    trcs_solml_vr(ids_H2PO4B,L,NY,NX) = trcn_irrig_vr(ids_H2PO4B,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw
    trcs_solml_vr(ids_H1PO4,L,NY,NX)  = trcn_irrig_vr(ids_H1PO4,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
    trcs_solml_vr(ids_H1PO4B,L,NY,NX) = trcn_irrig_vr(ids_H1PO4B,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw

    trcs_solml_vr(ids_NO2,L,NY,NX)  = 0._r8
    trcs_solml_vr(ids_NO2B,L,NY,NX) = 0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS, ANIONS AND ION PAIRS IN SOIL
!
    if(salt_model)then
      do nsalts=idsalt_beg,idsalt_KSO4
        trcSalt_solml_vr(nsalts,L,NY,NX)=trcSalt_irrig_vr(nsalts,L,NY,NX)*VLWatMicP_vr(L,NY,NX)
      ENDDO

      DO nsalts=idsalt_H0PO4,idsalt_MgHPO4
        ids=nsalts-idsalt_H0PO4+idsalt_H0PO4B
        trcSalt_solml_vr(nsalts,L,NY,NX)=trcSalt_irrig_vr(nsalts,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcSalt_solml_vr(ids,L,NY,NX)=trcSalt_irrig_vr(nsalts,L,NY,NX)*VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
      ENDDO
      trcSalt_soHml_vr(idsalt_beg:idsaltb_end,L,NY,NX)=0._r8
    endif

!
!     INITIAL STATE VARIABLES FOR ALL MATERIAL IN SOIL MACROPORES
!
    trcs_soHml_vr(ids_beg:ids_end,L,NY,NX)=0._r8
!
!     INITIAL STATE VARIABLES FOR EXCHANGEABLE CATIONS AND ANIONS
!
    trcx_solml_vr(idx_NH4,L,NY,NX)    = solutevar%XNH4_mole_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
    trcx_solml_vr(idx_NH4B,L,NY,NX)   = solutevar%XNH4_mole_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
    trcx_solml_vr(idx_Hp,L,NY,NX)     = solutevar%XHY1*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_Al,L,NY,NX)     = solutevar%XAl_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_Fe,L,NY,NX)     = solutevar%XFe_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_Ca,L,NY,NX)     = solutevar%XCa_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_Mg,L,NY,NX)     = solutevar%XMg_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_Na,L,NY,NX)     = solutevar%XNa_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_K,L,NY,NX)      = solutevar%XK_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_COOH,L,NY,NX)   = solutevar%XHC1*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_AlOH2,L,NY,NX)  = solutevar%XAlO2H2_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_FeOH2,L,NY,NX)  = solutevar%XFeO2H2_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcx_solml_vr(idx_OHe,L,NY,NX)    = solutevar%XOH_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcx_solml_vr(idx_OH,L,NY,NX)     = solutevar%XROH1_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcx_solml_vr(idx_OHp,L,NY,NX)    = solutevar%XROH2_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcx_solml_vr(idx_HPO4,L,NY,NX)   = solutevar%XHPO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcx_solml_vr(idx_H2PO4,L,NY,NX)  = solutevar%XH2PO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcx_solml_vr(idx_OHeB,L,NY,NX)   = solutevar%XOH_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcx_solml_vr(idx_OHB,L,NY,NX)    = solutevar%XROH1_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcx_solml_vr(idx_OHpB,L,NY,NX)   = solutevar%XROH2_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcx_solml_vr(idx_HPO4B,L,NY,NX)  = solutevar%XHPO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcx_solml_vr(idx_H2PO4B,L,NY,NX) = solutevar%XH2PO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
!
!     INITIAL STATE VARIABLES FOR PRECIPITATES
!
    trcp_saltpml_vr(idsp_AlOH3,L,NY,NX)     = solutevar%Precp_AlO3H3_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcp_saltpml_vr(idsp_FeOH3,L,NY,NX)     = solutevar%Precp_FeO3H3_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcp_saltpml_vr(idsp_CaCO3,L,NY,NX)     = solutevar%Precp_CaCO3_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcp_saltpml_vr(idsp_CaSO4,L,NY,NX)     = solutevar%Precp_CaSO4_conc*VLSoilMicPMass_vr(L,NY,NX)
    trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)     = solutevar%Precp_AlPO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcp_saltpml_vr(idsp_FePO4,L,NY,NX)     = solutevar%Precp_FePO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)    = solutevar%Precp_CaHPO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcp_saltpml_vr(idsp_HA,L,NY,NX)        = solutevar%Precp_Ca5P3O12O3H3_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)  = 0._r8
    trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)    = solutevar%Precp_AlPO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)    = solutevar%Precp_FePO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)   = solutevar%Precp_CaHPO4_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcp_saltpml_vr(idsp_HAB,L,NY,NX)       = solutevar%Precp_Ca5P3O12O3H3_conc*VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX) = 0._r8
    ElectricConductivity_vr(L,NY,NX)        = 0._r8
    SolutesIonStrenth_vr(L,NY,NX)           = 0._r8
    SolutesIonConc_vr(L,NY,NX)              = 0._r8
! the following line is quite interesting, Jinyun Tang, Nov 17, 2022
    trcs_solml_vr(ids_NH4,L,NY,NX)               = trcs_solml_vr(ids_NH4,L,NY,NX)+0.5_r8*SolidOM_vr(ielmn,iprotein,k_manure,L,NY,NX)
    SolidOM_vr(ielmn,iprotein,k_manure,L,NY,NX) = SolidOM_vr(ielmn,iprotein,k_manure,L,NY,NX)-0.5_r8*SolidOM_vr(ielmn,iprotein,k_manure,L,NY,NX)
  ENDIF
  end associate
  end subroutine InitSoluteConcs
!------------------------------------------------------------------------------------------

  subroutine InitialState(NY,NX)

  implicit none
  integer, intent(in) :: NY, NX
  real(R8) :: VOLWW
  integer :: nsalts
  integer :: L,idg
!
!     INITIAL STATE VARIABLES FOR MINERALS IN SURFACE RESIDUE
! not restart run and first year
  IF(.not.is_restart().AND.is_first_year)THEN
    trcs_solml_vr(ids_nuts_beg:ids_nuts_end,0,NY,NX)=0._r8

    trcx_solml_vr(idx_NH4,0,NY,NX)    = 0._r8
    trcx_solml_vr(idx_NH4B,0,NY,NX)   = 0._r8
    trcx_solml_vr(idx_OHe,0,NY,NX)    = 0._r8
    trcx_solml_vr(idx_OH,0,NY,NX)     = 0._r8
    trcx_solml_vr(idx_OHp,0,NY,NX)    = 0._r8
    trcx_solml_vr(idx_HPO4,0,NY,NX)   = 0._r8
    trcx_solml_vr(idx_H2PO4,0,NY,NX)  = 0._r8
    trcx_solml_vr(idx_NH4B,0,NY,NX)   = 0._r8
    trcx_solml_vr(idx_OHeB,0,NY,NX)   = 0._r8
    trcx_solml_vr(idx_OHB,0,NY,NX)    = 0._r8
    trcx_solml_vr(idx_OHpB,0,NY,NX)   = 0._r8
    trcx_solml_vr(idx_HPO4B,0,NY,NX)  = 0._r8
    trcx_solml_vr(idx_H2PO4B,0,NY,NX) = 0._r8

    trcp_saltpml_vr(idsp_AlPO4,0,NY,NX)     = 0._r8
    trcp_saltpml_vr(idsp_FePO4,0,NY,NX)     = 0._r8
    trcp_saltpml_vr(idsp_CaHPO4,0,NY,NX)    = 0._r8
    trcp_saltpml_vr(idsp_HA,0,NY,NX)        = 0._r8
    trcp_saltpml_vr(idsp_CaH4P2O8,0,NY,NX)  = 0._r8
    trcp_saltpml_vr(idsp_AlPO4B,0,NY,NX)    = 0._r8
    trcp_saltpml_vr(idsp_FePO4B,0,NY,NX)    = 0._r8
    trcp_saltpml_vr(idsp_CaHPO4B,0,NY,NX)   = 0._r8
    trcp_saltpml_vr(idsp_HAB,0,NY,NX)       = 0._r8
    trcp_saltpml_vr(idsp_CaH4P2O8B,0,NY,NX) = 0._r8
!
!     INITIAL STATE VARIABLES FOR MINERAL N AND P IN SNOWPACK
!
    D9985: DO L=1,JS
      IF(VLHeatCapSnow_snvr(L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
        VOLWW                             = VLWatSnow_snvr(L,NY,NX)+VLDrySnoWE_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE

        DO idg=idg_beg,idg_NH3
          trcg_solsml_snvr(idg,L,NY,NX) = VOLWW*trcg_rain_mole_conc_col(idg,NY,NX)*MolecularWeight(idg)
        ENDDO

        trcn_solsml_snvr(ids_NH4,L,NY,NX)   = VOLWW*NH4_rain_mole_conc(NY,NX)*natomw
        trcn_solsml_snvr(ids_NO3,L,NY,NX)   = VOLWW*NO3_rain_mole_conc(NY,NX)*natomw
        trcn_solsml_snvr(ids_H1PO4,L,NY,NX) = VOLWW*HPO4_rain_mole_conc(NY,NX)*patomw
        trcn_solsml_snvr(ids_H2PO4,L,NY,NX) = VOLWW*H2PO4_rain_mole_conc(NY,NX)*patomw
!
!     INITIAL STATE VARIABLES FOR CATIONS AND ANIONS IN SNOWPACK
!
        IF(salt_model)THEN
          do nsalts=idsalt_beg,idsalt_end
            trc_Saltml_snvr(nsalts,L,NY,NX)=VOLWW*trcsalt_rain_mole_conc_col(nsalts,NY,NX)
          enddo
        ENDIF
      ELSE
        trcg_solsml_snvr(idg_beg:idg_NH3,L,NY,NX)        = 0._r8
        trcn_solsml_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = 0._r8
!
!     INITIAL STATE VARIABLES FOR CATIONS AND ANIONS IN SNOWPACK
!
        IF(salt_model)THEN
          trc_Saltml_snvr(idsalt_beg:idsalt_end,L,NY,NX)=0._r8
        ENDIF
      ENDIF
    ENDDO D9985
  ENDIF
  end subroutine InitialState

end module StarteMod
