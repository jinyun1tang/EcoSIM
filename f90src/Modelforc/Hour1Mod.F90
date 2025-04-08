module Hour1Mod
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: GravAcceleration=>DAT_CONST_G
  use minimathmod,    only: isclose, AZMAX1, AZMIN1
  use abortutils,     only: endrun,  print_info
  use fileUtil,      only: iulog
  use PlantMod,      only: PlantCanopyRadsModel
  use FertilizerMod, only: ApplyFertilizerAtNoon
  use EcoSIMConfig     , only : jcplx=>jcplxc,nlbiomcp=>NumLiveMicrbCompts
  use EcoSIMConfig     , only : ndbiomcp=>NumDeadMicrbCompts,jsken=>jskenc
  use EcoSIMConfig     , only : NumMicbFunGrupsPerCmplx=>NumMicbFunGrupsPerCmplx,do_instequil
  use EcoSiMParDataMod,  only : micpar, pltpar
  use SoilBGCNLayMod,    only : sumORGMLayL
  use PlantMgmtDataType, only : NP
  use BalancesMod,       only : SummarizeTracerMass,BegCheckBalances
  use DebugToolMod
  use ATSUtilsMod
  use TracerPropMod
  use TracerIDMod
  use EcoSimConst
  use EcoSIMCtrlMod
  use MiniFuncMod
  use SoilHydroParaMod
  use MicrobialDataType
  use SOMDataType
  use ChemTranspDataType
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
  implicit none

  private

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: hour1
  public :: InitHour1

  real(r8) :: FoliarWatRetcap(0:3)
  real(r8), pointer :: THETRX(:)
  real(r8), parameter :: mGravAccelerat=1.e-3_r8*GravAcceleration  !gravitational constant devided by 1000.
!
!     FoliarWatRetcap=foliar water retention capacity (m3 m-2)
!     THETRX=litter water retention capacity (m3 g C-1)

  contains

  subroutine InitHour1(NumOfLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  allocate(THETRX(1:NumOfLitrCmplxs))

  FoliarWatRetcap=real((/5.0E-04,2.5E-04,2.5E-04,2.5E-04/),r8)
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

  character(len=*), parameter :: subname='hour1'
  integer :: L,NX,NY
  real(r8) :: THETPZ_vr(JZ)   !air-filled soil pore
  real(r8) :: DPTH0           !water+ice thickness in litter

  integer :: NZ,NR,K
  logical :: dosum

!     execution begins here
  call PrintInfo('beg '//subname)
!
  if(lverb)write(*,*)'ResetLndscapeAccumlators'
  call ResetLndscapeAccumlators()

  if(lverb)write(*,*)'SetAtmsTracerConc'
  call SetAtmsTracerConc(I,J,NHW,NHE,NVN,NVS)
!
!     RESET FLUX ARRAYS USED IN OTHER SUBROUTINES
!
  if(lverb)write(*,*)'ResetFluxArrays'
  call ResetFluxArrays(I,NHW,NHE,NVN,NVS)
!
!     IF SALT FLAG SET
!
  IF(salt_model)THEN
    if(lverb)write(*,*)'ResetSaltModelArrays'
    call ResetSaltModelArrays(NHW,NHE,NVN,NVS)
  ENDIF
!
!     RESET SOIL PROPERTIES AND PEDOTRANSFER FUNCTIONS
!     FOLLOWING ANY SOIL DISTURBANCE
!
  if(lverb)write(*,*)'set atms gas conc'
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS  

      IF(J.EQ.1)THEN
        NumActivePlants(NY,NX)=0
        DO  NZ=1,NP(NY,NX)
          PSICanPDailyMin(NZ,NY,NX)=0._r8
        ENDDO
      ENDIF
    ENDDO  
  ENDDO

!     HYDROLOGICAL PRPOERTIES OR SURFACE LITTER
  if(lverb)write(*,*)'UpdateLiterPropertz'
  call UpdateLiterPropertz(NHW,NHE,NVN,NVS)
!
!     RESET SURFACE LITTER PHYSICAL PROPERTIES (DENSITY, TEXTURE)
!     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
  if(lverb)write(*,*)'SetLiterSoilPropAftDisturb'
  call SetLiterSoilPropAftDisturb(I,J,NHW,NHE,NVN,NVS,dosum)

  if(dosum) call SummarizeTracerMass(I,J,NHW,NHE,NVN,NVS)  

  call BegCheckBalances(I,J,NHW,NHE,NVN,NVS)

  if(lverb)write(*,*)'SetSurfaceProp4SedErosion'
  call SetSurfaceProp4SedErosion(NHW,NHE,NVN,NVS)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
!
!
!     PARAMETERS FOR COHESION, EROSIVITY, AND ROUGHNESS OF SURFACE SOIL USED
!     FOR SURFACE WATER AND SEDIMENT TRANSPORT IN 'EROSION'
!
      if(lverb)write(*,*)'RESET HOURLY ACCUMULATORS'
      call SetHourlyDiagnostics(I,J,NY,NX)
!
!     RESET ARRAYS TO TRANSFER MATERIALS WITHIN SOILS
!     AND BETWEEN SOILS AND PLANTS
!
      call SetArrays4PlantSoilTransfer(NY,NX)
!
!     IF SOC FLAG IS SET
!
      IF(iErosionMode.EQ.ieros_frzthawsom .OR. iErosionMode.EQ.ieros_frzthawsomeros)THEN
        call UpdateTotalSOC(NY,NX)
      ENDIF

      call ZeroHourlyArrays(NY,NX)

      call GetChemicalConcsInSoil(I,J,NY,NX,THETPZ_vr)

      call GetSoluteConcentrations(NY,NX)

      call Prep4PlantMicrobeUptake(NY,NX)

      call CalGasSolubility(NY,NX)

      call GetSoilHydraulicVars(I,J,NY,NX)

!     CALCULATE ACTIVE LAYER DEPTH
      call DiagActiveLayerDepth(NY,NX)

!     OUTPUT FOR WATER TABLE DEPTH
      call DiagWaterTBLDepz(NY,NX,THETPZ_vr)

      call GetSurfResidualProperties(I,J,NY,NX,DPTH0)

      call SetTracerPropertyInLiterAir(NY,NX)

      if(do_instequil)call ForceGasAquaEquil(NY,NX)
!
      call PlantCanopyRadsModel(I,J,NY,NX,DPTH0)
!
      if(lverb)write(*,*)'RESET HOURLY INDICATORS'
!
      LWRadCanGPrev_col(NY,NX)         = LWRadCanG(NY,NX)
      LWRadGrnd(NY,NX)                 = LWRadBySurf_col(NY,NX)
      NetCO2Flx2Canopy_col(NY,NX)      = Eco_NEE_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      LWRadCanG(NY,NX)                 = 0._r8
      LWRadBySurf_col(NY,NX)           = 0._r8
      TLEX_col(NY,NX)                  = Air_Heat_Latent_store_col(NY,NX)
      TSHX_col(NY,NX)                  = Air_Heat_Sens_store_col(NY,NX)

      Air_Heat_Latent_store_col(NY,NX) = 0._r8
      Air_Heat_Sens_store_col(NY,NX)   = 0._r8
      Eco_NetRad_col(NY,NX)            = 0._r8
      Eco_Heat_Latent_col(NY,NX)       = 0._r8
      Eco_Heat_Sens_col(NY,NX)         = 0._r8
      Eco_Heat_GrndSurf_col(NY,NX)     = 0._r8
      Canopy_NEE_col(NY,NX)            = 0._r8
      Eco_NEE_col(NY,NX)               = 0._r8
      ECO_ER_col(NY,NX)                = 0._r8

      DO  NZ=1,NP(NY,NX)
!
!     NUMBERS OF TOP AND BOTTOM ROOTED SOIL LAYERS
!
!     NG=number of uppermost rooted layer
!     NIXBotRootLayer_rpft=number of lowest rooted layer
!
        NGTopRootLayer_pft(NZ,NY,NX)  = MAX(NGTopRootLayer_pft(NZ,NY,NX),NU(NY,NX))
        NIXBotRootLayer_pft(NZ,NY,NX) = MAX(NIXBotRootLayer_pft(NZ,NY,NX),NU(NY,NX))
        DO  NR=1,NumOfCanopyLayers
          NIXBotRootLayer_rpft(NR,NZ,NY,NX)=MAX(NIXBotRootLayer_rpft(NR,NZ,NY,NX),NU(NY,NX))
        ENDDO
      ENDDO

      if(lverb)write(*,*)'CanopyInterceptPrecp'
      CALL CanopyInterceptPrecp(NY,NX)
!
!     WRITE SW AND PAR ALBEDO

    ENDDO
  ENDDO
!
  if(lverb)write(*,*)'ApplyFertilizerAtNoon'
!     FERTILIZER APPLICATIONS OCCUR AT SOLAR NOON
  call ApplyFertilizerAtNoon(I,J,NHW,NHE,NVN,NVS)

  call PrintInfo('end '//subname)
  END subroutine hour1
!------------------------------------------------------------------------------------------

  subroutine CanopyInterceptPrecp(NY,NX)
  !
  !DESCRIPTION
  !precipitation intercepation by canopy
  implicit none
  integer, intent(in) :: NY,NX
  integer :: NZ
  real(r8) :: CanopyWatHeldCap  !maximum precipitation holding capacity by canopy (leaf+stem) [m3 H2O]
  real(r8) :: prec2canopy_pft   !precipiation onto canopy [m H2O/h]
!
!     CANOPY RETENTION OF PRECIPITATION
!
!     FoliarWatRetcap=foliar surface water retention capacity
!     CanopyLeafArea_pft,CanopyStemArea_pft=leaf,stalk area of PFT
!     FLWC,TFLWC=water retention of PFT,combined canopy
!     PRECA=precipitation+irrigation
!     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
!     VOLWC=canopy surface water retention
!
!     Warning: No snofall intercepation is considered at the moment.

  DO  NZ=1,NP(NY,NX)
    CanopyWatHeldCap                 = FoliarWatRetcap(iPlantRootProfile_pft(NZ,NY,NX)) &
      *(CanopyLeafArea_pft(NZ,NY,NX)+CanopyStemArea_pft(NZ,NY,NX))
      
    prec2canopy_pft                  = PrecRainAndIrrig_col(NY,NX)*FracPARads2Canopy_pft(NZ,NY,NX)
    PrecIntcptByCanopy_pft(NZ,NY,NX) = AZMAX1(AMIN1(prec2canopy_pft,CanopyWatHeldCap-WatHeldOnCanopy_pft(NZ,NY,NX)))
    Prec2Canopy_col(NY,NX)           = Prec2Canopy_col(NY,NX)+prec2canopy_pft
    PrecIntceptByCanopy_col(NY,NX)   = PrecIntceptByCanopy_col(NY,NX)+PrecIntcptByCanopy_pft(NZ,NY,NX)
  ENDDO
  RainPrecThrufall_col(NY,NX) = PrecRainAndIrrig_col(NY,NX)-PrecIntceptByCanopy_col(NY,NX)

  end subroutine CanopyInterceptPrecp

!------------------------------------------------------------------------------------------
  subroutine ResetLndscapeAccumlators()
!     RESET HOURLY SOIL ACCUMULATORS FOR WATER, HEAT, GASES, SOLUTES
!
  implicit none
  WatMassStore_lnd     = 0._r8
  HeatStore_lnd        = 0._r8
  TSoilO2G_lnd         = 0._r8
  TSoilH2G_lnd         = 0._r8
  TSEDSO               = 0._r8
  LitRMStoreLndscap(:) = 0._r8

  POMHumStoreLndscap(:)        = 0._r8
  TGasC_lnd                    = 0._r8
  TGasN_lnd                    = 0._r8
  TDisolNH4_lnd                = 0._r8
  tNO3_lnd                     = 0._r8
  TDisolPi_lnd                 = 0._r8
  TION                         = 0._r8
  PlantElemntStoreLandscape(:) = 0._r8
  end subroutine ResetLndscapeAccumlators
!------------------------------------------------------------------------------------------

  subroutine SetAtmsTracerConc(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS

  integer :: NX, NY, idg
  real(r8) :: tPBOT,tmp
!     begin_execution
!
!     CONCENTRATIONS OF CO2, CH4, O2, N2, N2O, NH3, H2 IN ATMOSPHERE,
!     PRECIPITATION AND IRRIGATION FROM MIXING RATIOS READ IN 'READS'
!
!
! SurfIrrig_IonStrenth_col: surface irrigation ion strength, [g m-3]
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      !obtain the mass density using the ideal gas law, taking TREF as reference tempeature
      tPBOT                        = PBOT_col(NY,NX)/1.01325E+02_r8
      tmp                          = Tref/TairKClimMean(NY,NX)*tPBOT
      CCO2EI(NY,NX)                = CO2EI(NY,NX)*5.36E-04_r8*tmp
      
      AtmGasCgperm3(idg_CO2,NY,NX) = CO2E_col(NY,NX)*5.36E-04_r8*tmp !gC/m3
      AtmGasCgperm3(idg_CH4,NY,NX) = CH4E_col(NY,NX)*5.36E-04_r8*tmp !gC/m3
      AtmGasCgperm3(idg_O2,NY,NX)  = OXYE_col(NY,NX)*1.43E-03_r8*tmp !gO/m3
      AtmGasCgperm3(idg_N2,NY,NX)  = Z2GE_col(NY,NX)*1.25E-03_r8*tmp !gN/m3
      AtmGasCgperm3(idg_N2O,NY,NX) = Z2OE_col(NY,NX)*1.25E-03_r8*tmp !gN/m3 
      AtmGasCgperm3(idg_AR,NY,NX)  = ARGE_col(NY,NX)*1.78E-02_r8*tmp  !gAr/m3      
      AtmGasCgperm3(idg_NH3,NY,NX) = ZNH3E_col(NY,NX)*6.25E-04_r8*tmp !gN/m3
      AtmGasCgperm3(idg_H2,NY,NX)  = H2GE_col(NY,NX)*8.92E-05_r8*tmp  !gN/m3

      DO idg=idg_beg,idg_NH3
        trcg_rain_mole_conc_col(idg,NY,NX) = AtmGasCgperm3(idg,NY,NX)*gas_solubility(idg,TCA_col(NY,NX)) &
           /(EXP(GasSechenovConst(idg)*SurfIrrig_IonStrenth_col(NY,NX))*MolecularWeight(idg))
        trcg_irrig_mole_conc_col(idg,NY,NX) = AtmGasCgperm3(idg,NY,NX)*gas_solubility(idg, TCA_col(NY,NX)) &
          /(EXP(GasSechenovConst(idg)*CSTRQ(I,NY,NX))*MolecularWeight(idg))
      ENDDO
      GDD_col(NY,NX) = GDD_col(NY,NX)+TCA_col(NY,NX)/24._r8
    ENDDO
  ENDDO
  end subroutine SetAtmsTracerConc
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

  XGridSurfRunoff_2DH(1:2,1:2,:,:)       = 0._r8
  HeatXGridBySurfRunoff_2DH(1:2,1:2,:,:) = 0._r8

  DO  NX=NHW,NHE+extragrid
    DO  NY=NVN,NVS+extragrid
!
!     WATER,SNOW,SOLUTE RUNOFF
!
      QWatIntLaterFlow_col(NY,NX)            = 0._r8
      QCanopyWat2Dist_col(NY,NX)             = 0._r8
      HeatCanopy2Dist_col(NY,NX)             = 0._r8
      HydroSufDOCFlx_col(NY,NX)              = 0._r8
      HydroSubsDOCFlx_col(NY,NX)             = 0._r8
      HydroSufDICFlx_col(NY,NX)              = 0._r8
      HydroSubsDICFlx_col(NY,NX)             = 0._r8
      HydroSubsDONFlx_col(NY,NX)             = 0._r8
      HydroSubsDINFlx_col(NY,NX)             = 0._r8
      HydroSubsDOPFlx_col(NY,NX)             = 0._r8
      HydroSubsDIPFlx_col(NY,NX)             = 0._r8
      SurfRunoffWatFluxM_2DH(:,NY,NX)        = 0._r8

      DOM_FloXSurRunoff_2D(idom_beg:idom_end,1:jcplx,1:2,1:2,NY,NX)=0._r8

      trcg_FloXSurRunoff_2D(idg_beg:idg_NH3,1:2,1:2,NY,NX)=0._r8
      trcn_FloXSurRunoff_2D(ids_nut_beg:ids_nuts_end,1:2,1:2,NY,NX)=0._r8

      DrySnoBySnoRedistrib_2DH(1:2,NY,NX)             = 0._r8
      WatBySnowRedistrib_2DH(1:2,NY,NX)            = 0._r8
      IceBySnowRedistrib_2DH(1:2,NY,NX)            = 0._r8
      HeatBySnowRedistrib_2DH(1:2,NY,NX)           = 0._r8
      trcg_FloXSnow_2DH(idg_beg:idg_NH3,1:2,NY,NX) = 0._r8

      trcn_FloXSnow_2DH(ids_NH4,1:2,NY,NX)   = 0._r8
      trcn_FloXSnow_2DH(ids_NO3,1:2,NY,NX)   = 0._r8
      trcn_FloXSnow_2DH(ids_H1PO4,1:2,NY,NX) = 0._r8
      trcn_FloXSnow_2DH(ids_H2PO4,1:2,NY,NX) = 0._r8
!
!
!     GAS AND SOLUTE FLUXES
!
      DO  L=0,NL(NY,NX)+1

        trcs_TransptMicP_3D(ids_beg:ids_end,1:3,L,NY,NX)=0._r8

        DOM_MicpTransp_3D(idom_beg:idom_end,1:jcplx,1:3,L,NY,NX)=0._r8
      ENDDO
!
!     BAND AND MACROPORE FLUXES
!
      trcs_deadroot2soil_col(idg_beg:idg_NH3,NY,NX) = 0._r8
      DO L=1,NL(NY,NX)+1
        WaterFlowSoiMicP_3D(1:3,L,NY,NX)  = 0._r8
        WaterFlowSoiMicPX_3D(1:3,L,NY,NX) = 0._r8
        WaterFlowSoiMacP_3D(1:3,L,NY,NX)  = 0._r8
        HeatFlow2Soil_3D(1:3,L,NY,NX)     = 0._r8
        RootCO2Ar2Soil_vr(L,NY,NX)        = 0._r8
        RootCO2Ar2Root_vr(L,NY,NX)        = 0._r8
        trcs_deadroot2soil_vr(idg_beg:idg_NH3,L,NY,NX)                 = 0._r8
        trcs_TransptMicP_3D(ids_beg:ids_end,1:3,L,NY,NX)              = 0._r8
        Gas_AdvDif_Flx_3D(idg_beg:idg_end,1:3,L,NY,NX)                = 0._r8
        DOM_Macp_Transp_flx_3D(idom_beg:idom_end,1:jcplx,1:3,L,NY,NX) = 0._r8

      ENDDO
    ENDDO
  ENDDO

!     IF EROSION FLAG SET
!
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    DO NX=NHW,NHE+extragrid
      DO NY=NVN,NVS+extragrid
        cumSed_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XSand_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XSilt_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XClay_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XNH4Soil_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XNH3Soil_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XUreaSoil_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XNO3Soil_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XNH4Band_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XNH3Band_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XUreaBand_Eros_2D(1:2,1:2,NY,NX)=0._r8
        XNO3Band_Eros_2D(1:2,1:2,NY,NX)=0._r8

        trcx_Eros_2D(idx_beg:idx_end,1:2,1:2,NY,NX)  = 0._r8
        trcp_Eros_2D(idsp_beg:idsp_end,1:2,1:2,NY,NX) = 0._r8

        OMEERhetr(:,:,:,1:2,1:2,NY,NX) = 0._r8
        OMEERauto(:,:,1:2,1:2,NY,NX)   = 0._r8

        OMBioResdu_Eros_2D(1:NumPlantChemElms,:,:,1:2,1:2,NY,NX) = 0._r8
        SolidOMAct_Eros_2D(:,:,1:2,1:2,NY,NX)                    = 0._r8
        SorbedOM_Eros_2D(1:NumPlantChemElms,:,1:2,1:2,NY,NX)   = 0._r8
        SolidOM_Eros_2D(1:NumPlantChemElms,:,:,1:2,1:2,NY,NX) = 0._r8
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

      trcSalt_FloXSurRunoff_2D(idsalt_beg:idsalt_end,1:2,1:2,NY,NX) = 0._r8
      trcSalt_FloXSnow_2DH(idsalt_beg:idsalt_end,1:2,NY,NX)             = 0._r8

      DO  L=1,NL(NY,NX)+1
        DO NSA=idsalt_beg,idsaltb_end
          trcSalt_TransptMicP_3D(NSA,1:3,L,NY,NX) = 0._r8
          trcSalt_TransptMacP_3D(NSA,1:3,L,NY,NX)      = 0._r8
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine ResetSaltModelArrays
!------------------------------------------------------------------------------------------

  subroutine SetSoilPropertyAftDisturb(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: VMINL,VSAND
  real(r8) :: PSISK(0:100)
  real(r8) :: THETK(100)
  real(r8) :: CORGCM
  real(r8) :: ParticleDens
  real(r8) :: SUM2,SUM1
  real(r8) :: VORGC
  real(r8) :: XK,YK
  integer :: L,K,N,M

  !     begin_execution
  !write(*,*) "In SetSoilPropertyAftDisturbance: "
  D9975: DO L=NUI(NY,NX),NLI(NY,NX)
    !
    !     AREA,DLYR=lateral(1,2), vertical(3) area,thickness of soil layer
    !     VOLT,VLSoilPoreMicP_vr,VLSoilMicP=layer volume including,excluding rock,macropores
    !removed ponding layer
    IF(SoilBulkDensity_vr(L,NY,NX).LE.ZERO .AND. DLYR_3D(3,L,NY,NX).LE.ZERO2)THEN
      VLWatMicP_vr(L,NY,NX) = 0._r8
      VLiceMicP_vr(L,NY,NX) = 0._r8
    ENDIF
    AREA(1,L,NY,NX)        = DLYR_3D(3,L,NY,NX)*DLYR_3D(2,L,NY,NX)
    AREA(2,L,NY,NX)        = DLYR_3D(3,L,NY,NX)*DLYR_3D(1,L,NY,NX)
    VGeomLayer_vr(L,NY,NX) = AREA(3,L,NY,NX)*DLYR_3D(3,L,NY,NX)

    VLSoilPoreMicP_vr(L,NY,NX) = AMAX1(VGeomLayer_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX),1.e-8_r8)
    IF(SoilBulkDensity_vr(L,NY,NX).LE.ZERO)THEN
      VLSoilMicP_vr(L,NY,NX)=VLSoilPoreMicP_vr(L,NY,NX)
    ENDIF
    !
    !     BKVL=soil mass
    !     C*=concentration,ORGC=SOC,SAND=sand,SILT=silt,CLAY=clay
    !     ParticleDens=particle density
    !     PrtcleDensitySurfLay_col=particle density of surface layer for use in erosion.f
    !     POROS=porosity used in diffusivity
    !     VOLA,VOLW,VOLI,VOLP=total,water-,ice-,air-filled micropore volume
    !     VOLAH,VOLWH,VOLIH,VOLPH=total,water-,ice-,air-filled macropore volume
    !     EHUM=fraction of microbial decomposition product allocated to humus
    !     EPOC=fraction of SOC decomposition product allocated to POC
    !     SRP=parameter for deviation from linear log-log water retention
    !     PSIMX,PSIMN,LOGPSIAtSat=log water potential at FC,WP,POROS
    !     PSISD,PSIMD=PSIMX-LOGPSIAtSat,PSIMN-PSIMX
    !     FC,WP=water contents at field capacity,wilting point
    !     FCL,LOGWiltPoint=log FC,WP
    !     FCD,PSD=FCL-LOGWiltPoint,log(POROS)-FCL
    !
    VLSoilMicPMass_vr(L,NY,NX)=SoilBulkDensity_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)

    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CSoilOrgM_vr(ielmc,L,NY,NX) = AMIN1(orgcden,SoilOrgM_vr(ielmc,L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX))
      CSAND_vr(L,NY,NX)              = SAND(L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
      CSILT(L,NY,NX)              = SILT(L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
      CCLAY_vr(L,NY,NX)              = CLAY(L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
    ELSE
      CSoilOrgM_vr(ielmc,L,NY,NX)=0._r8
      CSAND_vr(L,NY,NX) = 0._r8
      CSILT(L,NY,NX) = 0._r8
      CCLAY_vr(L,NY,NX) = 0._r8
    ENDIF
    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZERO)THEN
      CORGCM       = AZMAX1(AMIN1(1.0_r8,MWC2Soil*CSoilOrgM_vr(ielmc,L,NY,NX)))
      ParticleDens = 1.30_r8*CORGCM+2.66_r8*(1.0_r8-CORGCM)
      IF(L.EQ.NU(NY,NX))THEN
!surface layer
        POROS_vr(L,NY,NX)=AMAX1(POROS_vr(L,NY,NX),1.0_r8-(SoilBulkDensity_vr(L,NY,NX)/ParticleDens))
      ELSE
!  avoid float exception
        if(SoilBulkDensity_vr(L,NY,NX)>=ParticleDens)then
          POROS_vr(L,NY,NX)=1.0_r8
        else
          POROS_vr(L,NY,NX)=1.0_r8-(SoilBulkDensity_vr(L,NY,NX)/ParticleDens)
        endif
      ENDIF
    ELSE
      ParticleDens      = 0._r8
      POROS_vr(L,NY,NX) = 1.0_r8
    ENDIF
    !     VLMicP_vr(L,NY,NX)=AMAX1(POROS_vr(L,NY,NX)*VLSoilMicP_vr(L,NY,NX)
    !    2,VLWatMicP_vr(L,NY,NX)+VLiceMicP_vr(L,NY,NX))
    !     VLMacP_vr(L,NY,NX)=AMAX1(SoilFracAsMacP_vr(L,NY,NX)*VGeomLayer_vr(L,NY,NX)
    !    2,VLWatMacP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
    VLMicP_vr(L,NY,NX) = POROS_vr(L,NY,NX)*VLSoilMicP_vr(L,NY,NX)
    VLMacP_vr(L,NY,NX) = SoilFracAsMacP_vr(L,NY,NX)*VGeomLayer_vr(L,NY,NX)
    IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      VLsoiAirP_vr(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)) &
        +AZMAX1(VLMacP_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX))
    ELSE
      VLsoiAirP_vr(L,NY,NX)=0._r8
    ENDIF
    EHUM(L,NY,NX) = 0.200_r8+0.333_r8*AMIN1(0.5_r8,CCLAY_vr(L,NY,NX))
    EPOC(L,NY,NX) = 1.0_r8

    call SoilHydroProperty(I,J,L,NY,NX)
!
!     SOIL HEAT CAPACITY AND THERMAL CONDUCTIVITY OF SOLID PHASE
!     FROM SOC AND TEXTURE
!
!     VORGC,VMINL,VSAND=volume fractions of SOC,mineral,sand
!     STC,DTC=weighted thermal conductivity of soil solid component
!
    if(lverb)write(*,*)'setthermcond'
    IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      VORGC                        = CORGCM*SoilBulkDensity_vr(L,NY,NX)/ParticleDens
      VMINL                        = (CSILT(L,NY,NX)+CCLAY_vr(L,NY,NX))*SoilBulkDensity_vr(L,NY,NX)/ParticleDens
      VSAND                        = CSAND_vr(L,NY,NX)*SoilBulkDensity_vr(L,NY,NX)/ParticleDens
      NumerSolidThermCond_vr(L,NY,NX) = (1.253_r8*VORGC*9.050E-04_r8+0.514_r8*VMINL*1.056E-02_r8 &
        +0.386_r8*VSAND*2.112E-02_r8)*FracSoiAsMicP_vr(L,NY,NX) &
        +0.514_r8*ROCK_vr(L,NY,NX)*1.056E-02_r8
      DenomSolidThermCond_vr(L,NY,NX)=(1.253_r8*VORGC+0.514_r8*VMINL+0.386_r8*VSAND) &
        *FracSoiAsMicP_vr(L,NY,NX)+0.514_r8*ROCK_vr(L,NY,NX)
    ELSE
      NumerSolidThermCond_vr(L,NY,NX) = 0._r8
      DenomSolidThermCond_vr(L,NY,NX) = 0._r8
    ENDIF
  ENDDO D9975

  end subroutine SetSoilPropertyAftDisturb
!------------------------------------------------------------------------------------------

  subroutine ResetSurfResidualProperty(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: VxcessWatLitR,TVOLWI
  real(r8) :: ThetaWLitR
  real(r8) :: VOLIRZ,VWatLitrZ
  real(r8) :: XVOLW0,XVOLI0
!     begin_execution
!
!     FCR=litter water content at -0.01 MPa
!     THETY=litter hygroscopic water content
!
  CSoilOrgM_vr(ielmc,0,NY,NX)=orgcden
!
!     SOIL SURFACE WATER STORAGE CAPACITY
!
!     IDWaterTable=water table flag from site file
!     ExtWaterTable,ExtWaterTablet0_col=current,initial natural water table depth
!     TileWaterTable_col,DTBLD=current,initial artificial water table depth
!     SoilSurfRoughnesst0_col,ZW=soil,water surface roughness
!     VLWatHeldCapSurf_col=soil surface water retention capacity
!     VWatStoreCapSurf_col=VLWatHeldCapSurf_col accounting for above-ground water table
!     EHUM=fraction of microbial decompn product allocated to surface humus
!     EPOC=fraction of SOC decomposition product allocated to surface POC
!
  !Stationary water table
  IF(IDWaterTable_col(NY,NX).LE.1 .OR. IDWaterTable_col(NY,NX).EQ.3)THEN
    ExtWaterTable_col(NY,NX)=ExtWaterTablet0_col(NY,NX)
  !Mobile water table  
  ELSEIF(IDWaterTable_col(NY,NX).EQ.2 .OR. IDWaterTable_col(NY,NX).EQ.4)THEN
    ExtWaterTable_col(NY,NX)=ExtWaterTablet0_col(NY,NX)+CumDepz2LayBottom_vr(NU(NY,NX)-1,NY,NX)
  ENDIF

  IF(IDWaterTable_col(NY,NX).EQ.3.OR.IDWaterTable_col(NY,NX).EQ.4)THEN
    TileWaterTable_col(NY,NX)=DTBLD(NY,NX)
  ENDIF

  IF(SoilBulkDensity_vr(NU(NY,NX),NY,NX).GT.ZERO)THEN
    SoilSurfRoughnesst0_col(NY,NX)=0.020_r8
  ELSE
    SoilSurfRoughnesst0_col(NY,NX)=ZW
  ENDIF
  VLWatHeldCapSurf_col(NY,NX)=AMAX1(0.001_r8,0.112_r8*SoilSurfRoughnesst0_col(NY,NX)+&
    3.10_r8*SoilSurfRoughnesst0_col(NY,NX)**2._r8 &
    -0.012_r8*SoilSurfRoughnesst0_col(NY,NX)*SLOPE(0,NY,NX))*AREA(3,NU(NY,NX),NY,NX)

  VWatStoreCapSurf_col(NY,NX)=AMAX1(VLWatHeldCapSurf_col(NY,NX),-(ExtWaterTable_col(NY,NX)-&
    CumDepz2LayBottom_vr(NU(NY,NX)-1,NY,NX))*AREA(3,NU(NY,NX),NY,NX))

  SoilDepthMidLay_vr(NU(NY,NX),NY,NX)=CumDepz2LayBottom_vr(NU(NY,NX),NY,NX)-0.5_r8*DLYR_3D(3,NU(NY,NX),NY,NX)
  IF(VLSoilMicPMass_vr(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    CCLAY_vr(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)/VLSoilMicPMass_vr(NU(NY,NX),NY,NX)
    CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)/VLSoilMicPMass_vr(NU(NY,NX),NY,NX)
    CSAND_vr(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)/VLSoilMicPMass_vr(NU(NY,NX),NY,NX)
  ELSE
    CCLAY_vr(NU(NY,NX),NY,NX)=0._r8
    CSILT(NU(NY,NX),NY,NX)=0._r8
    CSAND_vr(NU(NY,NX),NY,NX)=0._r8
  ENDIF
  EHUM(0,NY,NX)=0.200_r8+0.333_r8*AMIN1(0.5_r8,CCLAY_vr(NU(NY,NX),NY,NX))
  EPOC(0,NY,NX)=0.150_r8
  end subroutine ResetSurfResidualProperty
!------------------------------------------------------------------------------------------

  subroutine SetLiterSoilPropAftDisturb(I,J,NHW,NHE,NVN,NVS,dosum)
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  logical :: dosum
  integer :: NY,NX

  real(r8) :: PSISK(0:100),THETK(100)
  REAL(R8) :: SUM2,SUM1
  real(r8) :: XK,YK
  integer :: K,M
!     begin_execution
!     iResetSoilProf_col=disturbance flag
!     SoilBulkDensity_vr,BKDSI=current,initial bulk density
!
  dosum=.false.
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      IF(iResetSoilProf_col(NY,NX).NE.ifalse)THEN

        if(lverb)write(*,*)'LitterHydroproperty'
        call LitterHydroproperty(NY,NX)
    !
    !   'RESET SOIL PHYSICAL PROPERTIES (DENSITY, TEXTURE)'
    !     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
    !
        if(lverb)write(*,*)'SetSoilPropertyAftDisturb'
        call SetSoilPropertyAftDisturb(I,J,NY,NX)
    !
    !   'SURFACE RESIDUE PROPERTIES'
        if(lverb)write(*,*)'ResetSurfResidualProperty'
        call ResetSurfResidualProperty(NY,NX)
    !
    !     iResetSoilProf_col=reset disturbance flag
    !
        iResetSoilProf_col(NY,NX) = ifalse
        dosum                     = .true.
      ENDIF
    ENDDO  
  ENDDO
  end subroutine SetLiterSoilPropAftDisturb
!------------------------------------------------------------------------------------------

  subroutine SetHourlyDiagnostics(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NX,NY

  integer :: L
!     begin_execution
  RootCO2AutorPrev_col(NY,NX)                = RootCO2Autor_col(NY,NX)

  trcg_snowMassloss_col(:,NY,NX)             = 0._r8
  trcn_snowMassloss_col(:,NY,NX)             = 0._r8
  trcs_RMicbUptake_col(:,NY,NX)              = 0._r8
  RGasNetProd_col(idg_beg:idg_NH3,NY,NX)     = 0._r8
  Soil_Gas_pressure_vr(:,NY,NX)              = 0._r8
  Gas_WetDeposition_col(:,NY,NX)             = 0._r8
  RootCO2Autor_col(NY,NX)                    = 0._r8
  QIceInflx_vr(:,NY,NX)                      = 0._r8
  QIceInflx_col(NY,NX)                       = 0._r8
  Rain2ExposedSurf_col(NY,NX)                = 0._r8
  Rain2Soil_col(NY,NX)                       = 0._r8
  Rain2LitR_col(NY,NX)                       = 0._r8
  RainPrec2Sno_col(NY,NX)                    = 0._r8
  QSnoIceXfer2Soil_col(NY,NX)                = 0._r8
  QSnoWatXfer2Soil_col(NY,NX)                = 0._r8
  QSnoHeatXfer2Soil_col(NY,NX)               = 0._r8
  QSnowH2Oloss_col(NY,NX)                    = 0._r8
  PrecHeat2Snow_col(NY,NX)                   = 0._r8
  Prec2Snow_col(NY,NX)                       = 0._r8
  ECO_HR_CO2_vr(:,NY,NX)                     = 0._r8
  ECO_HR_CH4_vr(:,NY,NX)                     = 0._r8
  ECO_HR_CO2_col(NY,NX)                      = 0._r8
  ECO_HR_CH4_col(NY,NX)                      = 0._r8
  Eco_RadSW_col(NY,NX)                       = 0._r8
  RootCO2Autor_vr(:,NY,NX)                   = 0._r8
  tRDIM2DOM_col(1:NumPlantChemElms,NY,NX)    = 0._r8
  QRunSurf_col(NY,NX)                        = 0._r8
  HeatRunSurf_col(NY,NX)                     = 0._r8
  Qinflx2Soil_col(NY,NX)                     = 0._r8
  HeatFlx2Grnd_col(NY,NX)                    = 0._r8
  DIC_mass_col(NY,NX)                        = 0._r8
  tMicBiome_col(1:NumPlantChemElms,NY,NX)    = 0._r8
  tSoilOrgM_col(1:NumPlantChemElms,NY,NX)    = 0._r8
  WatMass_col(NY,NX)                         = 0._r8
  HeatStore_col(NY,NX)                       = 0._r8
  tLitrOM_col(1:NumPlantChemElms,NY,NX)      = 0._r8
  tHumOM_col(1:NumPlantChemElms,NY,NX)       = 0._r8
  tNH4_col(NY,NX)                            = 0._r8
  tNO3_col(NY,NX)                            = 0._r8
  tHxPO4_col(NY,NX)                          = 0._r8
  tXPO4_col(NY,NX)                           = 0._r8
  UION_col(NY,NX)                            = 0._r8
  QDischar_col(NY,NX)                        = 0._r8
  PrecHeat_col(NY,NX)                        = 0._r8
  QDrain_col(NY,NX)                          = 0._r8
  HeatDrain_col(NY,NX)                       = 0._r8
  trcg_AquaADV_Snow2Litr_flx(:,NY,NX)        = 0._r8
  trcn_AquaADV_Snow2Litr_flx(:,NY,NX)        = 0._r8
  trcn_AquaADV_Snow2Soil_flx(:,NY,NX)        = 0._r8
  trcn_AquaADV_Snow2Band_flx(:,NY,NX)        = 0._r8
  trcg_AquaADV_Snow2Soil_flx(:,NY,NX)        = 0._r8
  RootCO2Ar2Soil_col(NY,NX)                  = 0._r8
  RootCO2Ar2Root_col(NY,NX)                  = 0._r8
  GasHydroLossFlx_col(idg_beg:idg_end,NY,NX) = 0._r8
  SurfGasEmisFlx_col(idg_beg:idg_NH3,NY,NX)          = 0._r8
  GasDiff2Surf_flx_col(idg_beg:idg_NH3,NY,NX)        = 0._r8
  WatFLo2LitR_col(NY,NX)                             = 0._r8
  HeatFLoByWat2LitR_col(NY,NX)                       = 0._r8
  TLitrIceFlxThaw_col(NY,NX)                         = 0._r8
  TLitrIceHeatFlxFrez_col(NY,NX)                     = 0._r8
  HeatByRad2Surf_col(NY,NX)                          = 0._r8
  HeatSensAir2Surf_col(NY,NX)                        = 0._r8
  HeatEvapAir2Surf_col(NY,NX)                        = 0._r8
  HeatSensVapAir2Surf_col(NY,NX)                     = 0._r8
  HeatNet2Surf_col(NY,NX)                            = 0._r8
  VapXAir2GSurf_col(NY,NX)                           = 0._r8
  THeatSnowThaw_col(NY,NX)                           = 0._r8
  THeatSoiThaw_col(NY,NX)                            = 0._r8
  trcs_TransptMacP_3D(:,:,:,:,:)                     = 0._r8
  trcg_DisolEvap_Atm2Soil_flx(idg_beg:idg_end,NY,NX) = 0._r8
  trcg_DisolEvap_Atm2Litr_flx(idg_beg:idg_NH3,NY,NX) = 0._r8

  TPlantRootH2OUptake_col(NY,NX)                   = 0._r8
  CanopyWat_col(NY,NX)                             = 0._r8
  WatHeldOnCanopy_col(NY,NX)                       = 0._r8
  Prec2Canopy_col(NY,NX)                           = 0._r8
  PrecIntceptByCanopy_col(NY,NX)                   = 0._r8
  QVegET_col(NY,NX)                                = 0._r8
  VapXAir2Canopy_col(NY,NX)                        = 0._r8
  HeatFlx2Canopy_col(NY,NX)                        = 0._r8
  CanopyHeatStor_col(NY,NX)                        = 0._r8
  
  TRootGasLossDisturb_col(idg_beg:idg_NH3,NY,NX) = 0._r8
  LitrFallStrutElms_col(:,NY,NX)                      = 0._r8
  StandingDeadStrutElms_col(1:NumPlantChemElms,NY,NX) = 0._r8
  PlantPopu_col(NY,NX)                                = 0._r8
! zero arrays in the snow layers
  WatConvSno2MicP_snvr(1:JS,NY,NX)   = 0._r8
  WatConvSno2MacP_snvr(1:JS,NY,NX)   = 0._r8
  HeatConvSno2Soi_snvr(1:JS,NY,NX)   = 0._r8
  WatConvSno2LitR_snvr(1:JS,NY,NX)   = 0._r8
  HeatConvSno2LitR_snvr(1:JS,NY,NX)  = 0._r8
  SnoXfer2SnoLay_snvr(1:JS,NY,NX)    = 0._r8
  WatXfer2SnoLay_snvr(1:JS,NY,NX)    = 0._r8
  IceXfer2SnoLay_snvr(1:JS,NY,NX)    = 0._r8
  HeatXfer2SnoLay_snvr(1:JS,NY,NX)   = 0._r8
  XPhaseChangeHeatL_snvr(1:JS,NY,NX) = 0._r8
  HeatSource_vr(:,NY,NX)             = 0._r8
  MoistSensDecomp_vr(:,NY,NX)        = 0._r8
  trcg_AquaAdv_flx_snvr(idg_beg:idg_NH3,1:JS,NY,NX) = 0._r8
  trcn_AquaAdv_flx_snvr(ids_nut_beg:ids_nuts_end,1:JS,NY,NX)      = 0._r8
  
  IF(salt_model)THEN
    trcSalt_AquaAdv_flx_snvr(idsalt_beg:idsalt_end,1:JS,NY,NX)=0._r8
    trcSalt_AquaADV_Snow2Litr_flx(:,NY,NX) = 0._r8
    trcSalt_AquaADV_Snow2Soil_flx(:,NY,NX) = 0._r8
    trcSalt_snowMassloss_col(:,NY,NX)      = 0._r8
  ENDIF
  end subroutine SetHourlyDiagnostics
!------------------------------------------------------------------------------------------

  subroutine SetArrays4PlantSoilTransfer(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

!     begin_execution

  LitrfalStrutElms_vr(1:NumPlantChemElms,1:jsken,1:pltpar%NumOfPlantLitrCmplxs,0:NL(NY,NX),NY,NX) = 0._r8
  HeatSource_col(NY,NX)                                                       = 0._r8
  REcoDOMProd_vr(idom_beg:idom_end,1:jcplx,0:NL(NY,NX),NY,NX)                 = 0._r8
  RProd_Hp_vr(0:NL(NY,NX),NY,NX)                                              = 0._r8
  trcn_GeoChem_soil_vr(ids_nut_beg:ids_nuts_end,0:NL(NY,NX),NY,NX)            = 0._r8
  TRChem_sol_NH3_soil_vr(0:NL(NY,NX),NY,NX)                                   = 0._r8
  TRChem_gas_NH3_geochem_vr(0:NL(NY,NX),NY,NX)                                = 0._r8
  trcx_TRSoilChem_vr(idx_beg:idx_end,0:NL(NY,NX),NY,NX)                       = 0._r8
  trcp_RChem_soil_vr(idsp_psoi_beg:idsp_psoi_end,0:NL(NY,NX),NY,NX)           = 0._r8
  TPlantRootH2OLoss_vr(0:NL(NY,NX),NY,NX)                                     = 0._r8
  THeatLossRoot2Soil_vr(0:NL(NY,NX),NY,NX)                                    = 0._r8
  Gas_Disol_Flx_vr(idg_beg:idg_end,0:NL(NY,NX),NY,NX)                         = 0._r8
  tRootMycoExud2Soil_vr(1:NumPlantChemElms,1:jcplx,NU(NY,NX):NL(NY,NX),NY,NX) = 0._r8
  RO2UptkSoilM_vr(1:NPH,NU(NY,NX):NL(NY,NX),NY,NX)                            = 0._r8
  RainLitr_col(NY,NX)                                                         = 0._r8
  end subroutine SetArrays4PlantSoilTransfer
!------------------------------------------------------------------------------------------

  subroutine DiagWaterTBLDepz(NY,NX,THETPZ_vr)
  !
  !Description:
  !Diagnose water table depth in the soil column
  implicit none
  integer, intent(in) :: NY,NX

  real(r8), intent(in) :: THETPZ_vr(JZ)  !air-filled porosity in layer
  real(r8) :: PSIEquil                   !equilibrium matric potential
  real(r8) :: THETW1
  real(r8) :: THETWM
  real(r8) :: THETPX
  real(r8) :: THETPW !     THETPW=minimum air-filled porosity for saturation (m3 m-3)  
  real(r8) :: THETWP  
  integer :: LL,L
  logical :: FoundWaterTable
!     begin_execution

  THETPW          = 0.01_r8
  THETWP          = 1.0_r8-THETPW
  FoundWaterTable = .false.

  DO L=NUI(NY,NX),NLI(NY,NX)
!     IDWaterTable=water table flag from site file
!     THETPZ,THETPW=current,minimum air-filled, porosity for water table
!     DPTH,ExtWaterTable=depth of soil layer midpoint, water table
!     PSIEquil=water potential in hydraulic equilibrium with layer below
!     THETW1,THETWP=water content at PSIEquil,minimum SWC for water table
!     DepzIntWTBL_col=water table depth
!
    IF(IDWaterTable_col(NY,NX).NE.0)THEN
      IF(FoundWaterTable)exit

      IF(THETPZ_vr(L).LT.THETPW .OR. L.EQ.NL(NY,NX))THEN
        FoundWaterTable=.true.
        IF(SoilDepthMidLay_vr(L,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN   !above external water table
          D5705: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
            IF(THETPZ_vr(LL).GE.THETPW .AND. LL.NE.NL(NY,NX))THEN
              !air-filled pore greater minimum, i.e. not saturated
              FoundWaterTable=.false.
              exit
            ELSE IF(SoilDepthMidLay_vr(LL,NY,NX).GE.ExtWaterTable_col(NY,NX))THEN
              !current layer is lower than external water table
              exit
            ENDIF
          END DO D5705
        ENDIF

          !THETPW=saturation criterion for water table identification
          IF(FoundWaterTable)THEN

            IF(THETPZ_vr(L).GE.THETPW .AND. L.NE.NL(NY,NX))THEN !saturated and inside the hydrologically active zone
              !not bottom layer, saturated
              !PSIeqv in saturated layer
              PSIEquil = PSISoilMatricP_vr(L+1,NY,NX)-mGravAccelerat*(SoilDepthMidLay_vr(L+1,NY,NX)-SoilDepthMidLay_vr(L,NY,NX))
              THETWM   = THETWP*POROS_vr(L,NY,NX)
              THETW1   = AMIN1(THETWM,EXP((LOGPSIAtSat(NY,NX)-LOG(-PSIEquil)) &
                *PSD_vr(L,NY,NX)/LOGPSIMXD(NY,NX)+LOGPOROS_vr(L,NY,NX)))

              IF(THETWM.GT.THETW1)THEN
                THETPX                 = AMIN1(1.0_r8,AZMAX1((THETWM-THETW_vr(L,NY,NX))/(THETWM-THETW1)))
                DepzIntWTBL_col(NY,NX) = CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)*(1.0_r8-THETPX)
              ELSE
                DepzIntWTBL_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)
              ENDIF
            ELSEIF(L.GT.NU(NY,NX))THEN
              !not bottom layer, and not topsoil layer, partially saturated
              PSIEquil = PSISoilMatricP_vr(L,NY,NX)-mGravAccelerat*(SoilDepthMidLay_vr(L,NY,NX)-SoilDepthMidLay_vr(L-1,NY,NX))
              THETWM   = THETWP*POROS_vr(L-1,NY,NX)
              THETW1   = AMIN1(THETWM,EXP((LOGPSIAtSat(NY,NX)-LOG(-PSIEquil)) &
                *PSD_vr(L-1,NY,NX)/LOGPSIMXD(NY,NX)+LOGPOROS_vr(L-1,NY,NX)))
              IF(THETWM.GT.THETW1)THEN
                THETPX                 = AMIN1(1.0_r8,AZMAX1((THETWM-THETW_vr(L-1,NY,NX))/(THETWM-THETW1)))
                DepzIntWTBL_col(NY,NX) = CumDepz2LayBottom_vr(L-1,NY,NX)-DLYR_3D(3,L-1,NY,NX)*(1.0_r8-THETPX)
              ELSE
                DepzIntWTBL_col(NY,NX)=CumDepz2LayBottom_vr(L-1,NY,NX)-DLYR_3D(3,L-1,NY,NX)
              ENDIF
            ELSE
              DepzIntWTBL_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)
            ENDIF
          ENDIF
        ENDIF
      
    ENDIF
  END DO
  end subroutine DiagWaterTBLDepz
!------------------------------------------------------------------------------------------

  subroutine SetSurfaceProp4SedErosion(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX
  real(r8) :: BKVLNX
  real(r8) :: CORGM
  real(r8) :: COHS
  real(r8) :: D50
  real(r8) :: VISCWL
  REAL(R8) :: ZD50
  !     begin_execution
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      !
      !     SoilDetachability4Erosion1=soil detachability from rainfall impact
      !     D50=average particle size
      !     CER,XER=parameters for runoff transport capacity
      !     ZD50=particle size effect on surface roughness
      !     VLS=hourly sinking rate
      !     COHS=soil cohesion
      !     SoilDetachability4Erosion2=soil detachability
      !     ZM=surface roughness used in runoff velocity calculation in watsub.f
      !
      SoilMicPMassLayerMX(NY,NX)=AZMAX1(SoilMicPMassLayerMn(NY,NX)+&
        MWC2Soil*SoilOrgM_vr(ielmc,NU(NY,NX),NY,NX))
      BKVLNX=SAND(NU(NY,NX),NY,NX)+SILT(NU(NY,NX),NY,NX) &
        +CLAY(NU(NY,NX),NY,NX)+1.82E-06*SoilOrgM_vr(ielmc,NU(NY,NX),NY,NX)
      IF(BKVLNX.GT.ZEROS(NY,NX))THEN
        CORGM                               = MWC2Soil*SoilOrgM_vr(ielmc,NU(NY,NX),NY,NX)/BKVLNX
        CSoilOrgM_vr(ielmc,NU(NY,NX),NY,NX) = orgcden*CORGM
        CSAND_vr(NU(NY,NX),NY,NX)              = SAND(NU(NY,NX),NY,NX)/BKVLNX
        CSILT(NU(NY,NX),NY,NX)              = SILT(NU(NY,NX),NY,NX)/BKVLNX
        CCLAY_vr(NU(NY,NX),NY,NX)              = CLAY(NU(NY,NX),NY,NX)/BKVLNX
      ELSE
        CORGM                               = 0._r8
        CSoilOrgM_vr(ielmc,NU(NY,NX),NY,NX) = 0._r8
        CSAND_vr(NU(NY,NX),NY,NX)              = 0._r8
        CSILT(NU(NY,NX),NY,NX)              = 1.0
        CCLAY_vr(NU(NY,NX),NY,NX)              = 0._r8
      ENDIF
      
      D50=1.0_r8*CCLAY_vr(NU(NY,NX),NY,NX)+10._r8*CSILT(NU(NY,NX),NY,NX) &
        +100._r8*CSAND_vr(NU(NY,NX),NY,NX)+100._r8*CORGM
      ZD50                    = 0.041*(ppmc*D50)**0.167_r8
      SoiSurfRoughness(NY,NX) = SoilSurfRoughnesst0_col(NY,NX)+ZD50+1.0_r8*VLitR_col(NY,NX)/AREA(3,0,NY,NX)
      IF(iErosionMode.EQ.ieros_frzthawsom .OR. iErosionMode.EQ.ieros_frzthawsomeros)THEN        
        CER(NY,NX)              = ((D50+5.0_r8)/0.32_r8)**(-0.6_r8)
        XER(NY,NX)              = ((D50+5.0_r8)/300._r8)**0.25_r8
        print*,'SoiSurfRoughness',SoiSurfRoughness(NY,NX)
        SoilDetachability4Erosion1(NY,NX)=ppmc*(1.0_r8+2.0_r8*(1.0_r8-CSILT(NU(NY,NX),NY,NX)-CORGM))
        COHS=2.0_r8+10._r8*(CCLAY_vr(NU(NY,NX),NY,NX)+CORGM) &
          +5.0_r8*(1.0_r8-EXP(-2.0E-06_r8*totRootLenDens_vr(NU(NY,NX),NY,NX)))
        SoilDetachability4Erosion2(NY,NX)=0.79_r8*EXP(-0.85_r8*AMAX1(1.0_r8,COHS))

        PrtcleDensitySurfLay_col(NY,NX) = 1.30_r8*CORGM+2.66_r8*(1.0_r8-CORGM)
        VISCWL                        = VISCW*EXP(0.533_r8-0.0267_r8*TCS_vr(0,NY,NX))
        VLS_col(NY,NX)                    = 3.6E+03_r8*9.8_r8*(PrtcleDensitySurfLay_col(NY,NX)-1.0_r8) &
          *(ppmc*D50)**2/(18.0_r8*VISCWL)
      ENDIF
    ENDDO
  ENDDO
  end subroutine SetSurfaceProp4SedErosion
!------------------------------------------------------------------------------------------

  subroutine UpdateTotalSOC(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: OC
  integer :: K,M,N,NGL,L,NB,NC
  real(r8) :: orgm(1:NumPlantChemElms)
  !     begin_execution
  !
  DO L=0,NL(NY,NX)
    call sumORGMLayL(L,NY,NX,ORGM,.true.)
    ORGCX_vr(L,NY,NX)=ORGM(ielmc)
  ENDDO
  end subroutine UpdateTotalSOC
!------------------------------------------------------------------------------------------

  subroutine DiagActiveLayerDepth(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: VLiceTot,VLPoreTot,VLiceTotL,VOLWTL,VLPoreTotL
  integer :: LL,L
  logical :: FoundActiveLayer
  logical :: goto5701

!     begin_execution
  FoundActiveLayer = .false.
  DO L=NUI(NY,NX),NLI(NY,NX)
!
!     VOLI,VOLIH=ice volume in micropores,macropores
!     VOLW,VOLWH=water volume in micropores,macropores
!     VOLA,VOLAH=total volume in micropores,macropores
!     ActiveLayDepZ_col=active layer depth
!     CDPTH,DLYR=depth to bottom,thickness of soil layer
!
    IF(FoundActiveLayer)exit
    !current layer
    VLiceTot=VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)
    VLPoreTot=VLMicP_vr(L,NY,NX)+VLMacP_vr(L,NY,NX)

    IF(VLPoreTot.GT.ZEROS2(NY,NX).AND.VLiceTot.GT.0.01*VLPoreTot)THEN
      !significant ice and pore
      D5700: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
        VLiceTotL  = VLiceMicP_vr(LL,NY,NX)+VLiceMacP_vr(LL,NY,NX)
        VOLWTL     = VLWatMicP_vr(LL,NY,NX)+VLWatMacP_vr(LL,NY,NX)
        VLPoreTotL = VLMicP_vr(LL,NY,NX)+VLMacP_vr(LL,NY,NX)
        !defined as significant increase in ice content
        goto5701=(VLPoreTotL.GT.ZEROS2(NY,NX).AND.VLiceTotL.LT.ZERO2*VLPoreTotL)
        if(goto5701)exit
      ENDDO D5700

      if(.not. goto5701)then
        IF(VLPoreTot.GT.ZEROS2(NY,NX))THEN
          ActiveLayDepZ_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)*AMIN1(1.0_r8,VLiceTot/VLPoreTot)
        ELSE
          ActiveLayDepZ_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)
        ENDIF
        FoundActiveLayer=.true.
      else
        !not active layer
        ActiveLayDepZ_col(NY,NX)=9999.0_r8
      endif
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
  integer :: K,L,idg
!     begin_execution
!
!     LITTER GAS CONCENTRATIOS
!
!     C*G=soil gas gaseous concentration
!     *E=atmospheric concentration
!     TKS_vr,TCS=litter temperature (K,C)
!     S*L=gas solubility
!     C*S=soil gas aqueous concentration
!
  trcg_gascl_vr(idg_CO2,0,NY,NX) = CO2E_col(NY,NX)*5.36E-04_r8*TREF/TKS_vr(0,NY,NX)
  trcg_gascl_vr(idg_CH4,0,NY,NX) = CH4E_col(NY,NX)*5.36E-04_r8*TREF/TKS_vr(0,NY,NX)
  trcg_gascl_vr(idg_O2,0,NY,NX)  = OXYE_col(NY,NX)*1.43E-03_r8*TREF/TKS_vr(0,NY,NX)
  trcg_gascl_vr(idg_N2,0,NY,NX)  = Z2GE_col(NY,NX)*1.25E-03_r8*TREF/TKS_vr(0,NY,NX)
  trcg_gascl_vr(idg_N2O,0,NY,NX) = Z2OE_col(NY,NX)*1.25E-03_r8*TREF/TKS_vr(0,NY,NX)
  trcg_gascl_vr(idg_NH3,0,NY,NX) = ZNH3E_col(NY,NX)*6.25E-04_r8*TREF/TKS_vr(0,NY,NX)
  trcg_gascl_vr(idg_H2,0,NY,NX)  = H2GE_col(NY,NX)*8.92E-05*TREF/TKS_vr(0,NY,NX)

! initialize all band nutrients to zero
  trc_solcl_vr(ids_nutb_beg:ids_nutb_end,0,NY,NX)=0._r8
  IF(VLWatMicP_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
! exclude NH3B,
    DO idg=idg_beg,idg_NH3
      trc_solcl_vr(idg,0,NY,NX)=AZMAX1(trcs_solml_vr(idg,0,NY,NX)/VLWatMicP_vr(0,NY,NX))
    ENDDO
  ELSE
    trc_solcl_vr(idg_beg:idg_NH3,0,NY,NX)=0._r8
  ENDIF
!
!     TFACL=temperature effect on diffusivity
!     *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
!     *SG PARAMETER statement above
!
  TFACL                      = TEFAQUDIF(TKS_vr(0,NY,NX))
  TScal4Difsvity_vr(0,NY,NX) = TFACL

  SoluteDifusvty_vr(idg_CO2,0,NY,NX)   = CLSG*TFACL
  SoluteDifusvty_vr(idg_CH4,0,NY,NX)   = CQSG*TFACL
  SoluteDifusvty_vr(idg_O2,0,NY,NX)    = OLSG*TFACL
  SoluteDifusvty_vr(idg_N2,0,NY,NX)    = ZLSG*TFACL
  SoluteDifusvty_vr(idg_NH3,0,NY,NX)   = ZNSG*TFACL
  SoluteDifusvty_vr(idg_H2,0,NY,NX)    = HLSG*TFACL
  SoluteDifusvty_vr(idg_N2O,0,NY,NX)   = ZVSG*TFACL
  SoluteDifusvty_vr(ids_NO3,0,NY,NX)   = ZOSG*TFACL
  SoluteDifusvty_vr(ids_H1PO4,0,NY,NX) = POSG*TFACL
  
  SoluteDifusvty_vr(ids_NH4,0,NY,NX)   =SoluteDifusvty_vr(idg_NH3,0,NY,NX)
  SoluteDifusvty_vr(ids_NH4B,0,NY,NX)  =SoluteDifusvty_vr(ids_NH4,0,NY,NX)
  SoluteDifusvty_vr(idg_NH3B,0,NY,NX)  =SoluteDifusvty_vr(idg_NH3,0,NY,NX)
  SoluteDifusvty_vr(ids_NO3B,0,NY,NX)  =SoluteDifusvty_vr(ids_NO3,0,NY,NX)
  SoluteDifusvty_vr(ids_NO2,0,NY,NX)   =SoluteDifusvty_vr(ids_NO3,0,NY,NX)
  SoluteDifusvty_vr(ids_NO2B,0,NY,NX)  =SoluteDifusvty_vr(ids_NO2,0,NY,NX)
  SoluteDifusvty_vr(ids_H2PO4,0,NY,NX) =SoluteDifusvty_vr(ids_H1PO4,0,NY,NX)
  SoluteDifusvty_vr(ids_H1PO4B,0,NY,NX)=SoluteDifusvty_vr(ids_H1PO4,0,NY,NX)
  SoluteDifusvty_vr(ids_H2PO4B,0,NY,NX)=SoluteDifusvty_vr(ids_H2PO4,0,NY,NX)

  DOMdiffusivity_vr(idom_doc,0,NY,NX)     = OCSG*TFACL
  DOMdiffusivity_vr(idom_don,0,NY,NX)     = ONSG*TFACL
  DOMdiffusivity_vr(idom_dop,0,NY,NX)     = OPSG*TFACL
  DOMdiffusivity_vr(idom_acetate,0,NY,NX) = OASG*TFACL
!
!     R*Y,R*X=total substrate uptake from previous,current hour
!     used in nitro.f, uptake.f
!
  RO2EcoDmndPrev_vr(0,NY,NX)        = REcoO2DmndResp_vr(0,NY,NX)
  RNH4EcoDmndSoilPrev_vr(0,NY,NX)   = REcoNH4DmndSoil_vr(0,NY,NX)
  RNO3EcoDmndSoilPrev_vr(0,NY,NX)   = REcoNO3DmndSoil_vr(0,NY,NX)
  RNO2EcoUptkSoilPrev_vr(0,NY,NX)   = RNO2EcoUptkSoil_vr(0,NY,NX)
  RN2OEcoUptkSoilPrev_vr(0,NY,NX)   = RN2OEcoUptkSoil_vr(0,NY,NX)
  RH1PO4EcoDmndSoilPrev_vr(0,NY,NX) = REcoH1PO4DmndSoil_vr(0,NY,NX)
  RH2PO4EcoDmndSoilPrev_vr(0,NY,NX) = REcoH2PO4DmndSoil_vr(0,NY,NX)
  REcoO2DmndResp_vr(0,NY,NX)        = 0._r8
  REcoNH4DmndSoil_vr(0,NY,NX)       = 0._r8
  REcoNO3DmndSoil_vr(0,NY,NX)       = 0._r8
  RNO2EcoUptkSoil_vr(0,NY,NX)       = 0._r8
  RN2OEcoUptkSoil_vr(0,NY,NX)       = 0._r8
  REcoH1PO4DmndSoil_vr(0,NY,NX)     = 0._r8
  REcoH2PO4DmndSoil_vr(0,NY,NX)     = 0._r8
  D5055: DO K=1,jcplx
    RDOMEcoDmndPrev_vr(K,0,NY,NX)     = RDOMEcoDmndK_vr(K,0,NY,NX)
    RAcetateEcoDmndPrev_vr(K,0,NY,NX) = RAcetateEcoDmndK_vr(K,0,NY,NX)
    RDOMEcoDmndK_vr(K,0,NY,NX)        = 0._r8
    RAcetateEcoDmndK_vr(K,0,NY,NX)    = 0._r8
  ENDDO D5055
!
!     WVapDifusvityAir_col,VaporDiffusivityLitR_col,WGSGW=vapor diffusivity in air,litter,snowpack
!
  TFACA                           = TEFGASDIF(TairK_col(NY,NX))
  WVapDifusvityAir_col(NY,NX)     = WGSG*TFACA
  TFACR                           = TEFGASDIF(TKS_vr(0,NY,NX))
  VaporDiffusivityLitR_col(NY,NX) = WGSG*TFACR
  D5060: DO  L=1,JS
    TFACW                   = TEFGASDIF(TKSnow_snvr(L,NY,NX))
    H2OVapDifscSno(L,NY,NX) = WGSG*TFACW
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
    IF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

      IF(trcs_VLN_vr(ids_NH4,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_NH4,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_NH4,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)))
        trc_solcl_vr(idg_NH3,L,NY,NX)=AZMAX1(trcs_solml_vr(idg_NH3,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(idg_NH3,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NH4,L,NY,NX)=0._r8
        trc_solcl_vr(idg_NH3,L,NY,NX)=0._r8
      ENDIF
      IF(trcs_VLN_vr(ids_NO3,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_NO3,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_NO3,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)))
        trc_solcl_vr(ids_NO2,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_NO2,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO2,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NO3,L,NY,NX)=0._r8
        trc_solcl_vr(ids_NO2,L,NY,NX)=0._r8
      ENDIF

      IF(trcs_VLN_vr(ids_H1PO4,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_H1PO4,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_H1PO4,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)))
        trc_solcl_vr(ids_H2PO4,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_H2PO4,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H2PO4,L,NY,NX)))

        CPO4S_vr(L,NY,NX)=trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4,L,NY,NX)

        if(salt_model)then
          CPO4S_vr(L,NY,NX)=CPO4S_vr(L,NY,NX)+(trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX))*patomw
        endif  

        CPO4S_vr(L,NY,NX)=AZMAX1(CPO4S_vr(L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)))

      ELSE
        trc_solcl_vr(ids_H1PO4,L,NY,NX)=0._r8
        trc_solcl_vr(ids_H2PO4,L,NY,NX)=0._r8
        CPO4S_vr(L,NY,NX)=0._r8
      ENDIF
!
!     C*B=solute concentration in band
!     CH1PB,CH2PB=HPO4,H2PO4 concentration in band
!     Z*B=P ion pair amounts in band (see solute.f)
!     VLNHB,VLNOB,VLPOB=fraction of soil volume in NH4,NO3,PO4 band
!
      IF(trcs_VLN_vr(ids_NH4B,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_NH4B,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_NH4B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)))
        trc_solcl_vr(idg_NH3B,L,NY,NX)=AZMAX1(trcs_solml_vr(idg_NH3B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(idg_NH3B,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NH4B,L,NY,NX)=0._r8
        trc_solcl_vr(idg_NH3B,L,NY,NX)=0._r8
      ENDIF

      IF(trcs_VLN_vr(ids_NO3B,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_NO3B,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_NO3B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)))
        trc_solcl_vr(ids_NO2B,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_NO2B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO2B,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NO3B,L,NY,NX)=0._r8
        trc_solcl_vr(ids_NO2B,L,NY,NX)=0._r8
      ENDIF
      IF(trcs_VLN_vr(ids_H1PO4B,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_H1PO4B,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_H1PO4B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)))
        trc_solcl_vr(ids_H2PO4B,L,NY,NX)=AZMAX1(trcs_solml_vr(ids_H2PO4B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H2PO4B,L,NY,NX)))

        CPO4B_vr(L,NY,NX)=trcs_solml_vr(ids_H1PO4B,L,NY,NX)+trcs_solml_vr(ids_H2PO4B,L,NY,NX)
        if(salt_model)then
          CPO4B_vr(L,NY,NX)=CPO4B_vr(L,NY,NX)+(trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX))*patomw 
        endif  
        CPO4B_vr(L,NY,NX)=AZMAX1(CPO4B_vr(L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_H1PO4B,L,NY,NX) = 0._r8
        trc_solcl_vr(ids_H2PO4B,L,NY,NX) = 0._r8
        CPO4B_vr(L,NY,NX)                   = 0._r8
      ENDIF

    ELSE
      trc_solcl_vr(ids_nuts_beg:ids_nuts_end,L,NY,NX) = 0._r8
      CPO4S_vr(L,NY,NX)                                  = 0._r8
      CPO4B_vr(L,NY,NX)                                  = 0._r8
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
    RO2EcoDmndPrev_vr(L,NY,NX)        = REcoO2DmndResp_vr(L,NY,NX)
    RNH4EcoDmndSoilPrev_vr(L,NY,NX)   = REcoNH4DmndSoil_vr(L,NY,NX)
    RNO3EcoDmndSoilPrev_vr(L,NY,NX)   = REcoNO3DmndSoil_vr(L,NY,NX)
    RNO2EcoUptkSoilPrev_vr(L,NY,NX)   = RNO2EcoUptkSoil_vr(L,NY,NX)
    RN2OEcoUptkSoilPrev_vr(L,NY,NX)   = RN2OEcoUptkSoil_vr(L,NY,NX)
    RH1PO4EcoDmndSoilPrev_vr(L,NY,NX) = REcoH1PO4DmndSoil_vr(L,NY,NX)
    RH2PO4EcoDmndSoilPrev_vr(L,NY,NX) = REcoH2PO4DmndSoil_vr(L,NY,NX)
    RNH4EcoDmndBandPrev_vr(L,NY,NX)   = REcoNH4DmndBand_vr(L,NY,NX)
    RNO3EcoDmndBandPrev_vr(L,NY,NX)   = REcoNO3DmndBand_vr(L,NY,NX)
    RNO2EcoUptkBandPrev_vr(L,NY,NX)   = RNO2EcoUptkBand_vr(L,NY,NX)
    RH1PO4EcoDmndBandPrev_vr(L,NY,NX) = REcoH1PO4DmndBand_vr(L,NY,NX)
    RH2PO4EcoDmndBandPrev_vr(L,NY,NX) = REcoH2PO4DmndBand_vr(L,NY,NX)
    REcoO2DmndResp_vr(L,NY,NX)        = 0._r8
    REcoNH4DmndSoil_vr(L,NY,NX)       = 0._r8
    REcoNO3DmndSoil_vr(L,NY,NX)       = 0._r8
    RNO2EcoUptkSoil_vr(L,NY,NX)       = 0._r8
    RN2OEcoUptkSoil_vr(L,NY,NX)       = 0._r8
    REcoH1PO4DmndSoil_vr(L,NY,NX)     = 0._r8
    REcoH2PO4DmndSoil_vr(L,NY,NX)     = 0._r8
    REcoNH4DmndBand_vr(L,NY,NX)       = 0._r8
    REcoNO3DmndBand_vr(L,NY,NX)       = 0._r8
    RNO2EcoUptkBand_vr(L,NY,NX)       = 0._r8
    REcoH1PO4DmndBand_vr(L,NY,NX)     = 0._r8
    REcoH2PO4DmndBand_vr(L,NY,NX)     = 0._r8

    D5050: DO K=1,jcplx
      RDOMEcoDmndPrev_vr(K,L,NY,NX)     = RDOMEcoDmndK_vr(K,L,NY,NX)
      RAcetateEcoDmndPrev_vr(K,L,NY,NX) = RAcetateEcoDmndK_vr(K,L,NY,NX)
      RDOMEcoDmndK_vr(K,L,NY,NX)        = 0._r8
      RAcetateEcoDmndK_vr(K,L,NY,NX)    = 0._r8
    ENDDO D5050
!
! DIFFUSIVITY
!
! TFACG,TFACL=temperature effects on gaseous,aqueous diffusivity
!
! *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
! *SG PARAMETER statement above
!
    TFACG                      = TEFGASDIF(TKS_vr(L,NY,NX))
    TFACL                      = TEFAQUDIF(TKS_vr(L,NY,NX))
    TScal4Difsvity_vr(L,NY,NX) = TFACL

    GasDifc_vr(idg_CO2,L,NY,NX)  = CGSG*TFACG
    GasDifc_vr(idg_CH4,L,NY,NX)  = CHSG*TFACG
    GasDifc_vr(idg_O2,L,NY,NX)   = OGSG*TFACG
    GasDifc_vr(idg_N2,L,NY,NX)   = ZGSG*TFACG
    GasDifc_vr(idg_N2O,L,NY,NX)  = Z2SG*TFACG
    GasDifc_vr(idg_NH3,L,NY,NX)  = ZHSG*TFACG
    GasDifc_vr(idg_H2,L,NY,NX)   = HGSG*TFACG
    GasDifc_vr(idg_NH3B,L,NY,NX) = ZHSG*TFACG
    GasDifc_vr(idg_AR,L,NY,NX)   = ARSG*TFACG

    SoluteDifusvty_vr(idg_CO2,L,NY,NX)   = CLSG*TFACL
    SoluteDifusvty_vr(idg_CH4,L,NY,NX)   = CQSG*TFACL
    SoluteDifusvty_vr(idg_O2,L,NY,NX)    = OLSG*TFACL
    SoluteDifusvty_vr(idg_N2,L,NY,NX)    = ZLSG*TFACL
    SoluteDifusvty_vr(idg_NH3,L,NY,NX)   = ZNSG*TFACL
    SoluteDifusvty_vr(idg_H2,L,NY,NX)    = HLSG*TFACL
    SoluteDifusvty_vr(idg_N2O,L,NY,NX)   = ZVSG*TFACL
    SoluteDifusvty_vr(idg_AR,L,NY,NX)    = ARSL*TFACL
    SoluteDifusvty_vr(idg_NH3B,L,NY,NX)  = SoluteDifusvty_vr(idg_NH3,L,NY,NX)
    SoluteDifusvty_vr(ids_NO3,L,NY,NX)   = ZOSG*TFACL
    SoluteDifusvty_vr(ids_H1PO4,L,NY,NX) = POSG*TFACL

    SoluteDifusvty_vr(ids_NH4,L,NY,NX)   =SoluteDifusvty_vr(idg_NH3,L,NY,NX)
    SoluteDifusvty_vr(ids_NH4B,L,NY,NX)  =SoluteDifusvty_vr(ids_NH4,L,NY,NX)
    SoluteDifusvty_vr(ids_NO3B,L,NY,NX)  =SoluteDifusvty_vr(ids_NO3,L,NY,NX)
    SoluteDifusvty_vr(ids_NO2,L,NY,NX)   =SoluteDifusvty_vr(ids_NO3,L,NY,NX)
    SoluteDifusvty_vr(ids_NO2B,L,NY,NX)  =SoluteDifusvty_vr(ids_NO2,L,NY,NX)
    SoluteDifusvty_vr(ids_H2PO4,L,NY,NX) =SoluteDifusvty_vr(ids_H1PO4,L,NY,NX)
    SoluteDifusvty_vr(ids_H1PO4B,L,NY,NX)=SoluteDifusvty_vr(ids_H1PO4,L,NY,NX)
    SoluteDifusvty_vr(ids_H2PO4B,L,NY,NX)=SoluteDifusvty_vr(ids_H2PO4,L,NY,NX)

    DOMdiffusivity_vr(idom_doc,L,NY,NX)     = OCSG*TFACL
    DOMdiffusivity_vr(idom_don,L,NY,NX)     = ONSG*TFACL
    DOMdiffusivity_vr(idom_dop,L,NY,NX)     = OPSG*TFACL
    DOMdiffusivity_vr(idom_acetate,L,NY,NX) = OASG*TFACL
    WVapDifusvitySoil_vr(L,NY,NX)           = WGSG*TFACG

    IF(salt_model)THEN
      AquaIonDifusivty_vr(idsalt_Al,L,NY,NX)   = ALSG*TFACL
      AquaIonDifusivty_vr(idsalt_Fe,L,NY,NX)   = FESG*TFACL
      AquaIonDifusivty_vr(idsalt_Hp,L,NY,NX)   = HYSG*TFACL
      AquaIonDifusivty_vr(idsalt_Ca,L,NY,NX)   = CASG*TFACL
      AquaIonDifusivty_vr(idsalt_Mg,L,NY,NX)   = GMSG*TFACL
      AquaIonDifusivty_vr(idsalt_Na,L,NY,NX)   = ANSG*TFACL
      AquaIonDifusivty_vr(idsalt_K,L,NY,NX)    = AKSG*TFACL
      AquaIonDifusivty_vr(idsalt_OH,L,NY,NX)   = OHSG*TFACL
      AquaIonDifusivty_vr(idsalt_CO3,L,NY,NX)  = C3SG*TFACL
      AquaIonDifusivty_vr(idsalt_HCO3,L,NY,NX) = HCSG*TFACL
      AquaIonDifusivty_vr(idsalt_SO4,L,NY,NX)  = SOSG*TFACL
      AquaIonDifusivty_vr(idsalt_Cl,L,NY,NX)   = CLSX*TFACL
!
!   TOTAL ION CONCENTRATION
!
!   ZC3,ZA3,ZC2,ZA2,ZC1,ZA1=total tri-,di-,univalent cations C,anions A
!   CSTR,SolutesIonConc_vr=ion strength, total ion concentration
!
      ZC3=trcSalt_solml_vr(idsalt_Al,L,NY,NX)+trcSalt_solml_vr(idsalt_Fe,L,NY,NX)
      ZA3=trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)
      ZC2=trcSalt_solml_vr(idsalt_Ca,L,NY,NX)+trcSalt_solml_vr(idsalt_Mg,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_AlOH,L,NY,NX)+trcSalt_solml_vr(idsalt_FeOH,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)
      ZA2=trcSalt_solml_vr(idsalt_SO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CO3,L,NY,NX)+trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H1PO4B,L,NY,NX)

      ZC1=(trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(ids_NH4B,L,NY,NX))/natomw+trcSalt_solml_vr(idsalt_Hp,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_Na,L,NY,NX)+trcSalt_solml_vr(idsalt_K,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_AlOH2,L,NY,NX)+trcSalt_solml_vr(idsalt_FeOH2,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_AlSO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeSO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_CaOH,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHCO3,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_MgOH2,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHCO3,L,NY,NX)&
        +trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)

      ZA1=(trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO3B,L,NY,NX))/natomw+trcSalt_solml_vr(idsalt_OH,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_HCO3,L,NY,NX)+trcSalt_solml_vr(idsalt_Cl,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_AlOH4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeOH4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_NaCO3,L,NY,NX)+trcSalt_solml_vr(idsalt_NaSO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_KSO4,L,NY,NX)+(trcs_solml_vr(ids_H2PO4,L,NY,NX) &
        +trcs_solml_vr(ids_H2PO4B,L,NY,NX))/patomw+trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)

      ZN=trcs_solml_vr(idg_CO2,L,NY,NX)/catomw+trcs_solml_vr(idg_CH4,L,NY,NX)/catomw+trcs_solml_vr(idg_O2,L,NY,NX)/32.0 &
        +(trcs_solml_vr(idg_N2,L,NY,NX)+trcs_solml_vr(idg_N2O,L,NY,NX)+trcs_solml_vr(idg_NH3,L,NY,NX)+trcs_solml_vr(idg_NH3B,L,NY,NX))/natomw &
        +trcSalt_solml_vr(idsalt_AlOH3,L,NY,NX)+trcSalt_solml_vr(idsalt_FeOH3,L,NY,NX)+trcSalt_solml_vr(idsalt_CaCO3,L,NY,NX)+trcSalt_solml_vr(idsalt_CaSO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_MgCO3,L,NY,NX)+trcSalt_solml_vr(idsalt_MgSO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)

      ZION1=ABS(3.0_r8*(ZC3-ZA3)+2.0_r8*(ZC2-ZA2)+ZC1-ZA1)
      IF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        SolutesIonStrenth_vr(L,NY,NX)=AZMAX1(0.5E-03_r8*(9.0_r8*(ZC3+ZA3)+4.0_r8*(ZC2+ZA2) &
          +ZC1+ZA1+ZION1)/VLWatMicP_vr(L,NY,NX))
        SolutesIonConc_vr(L,NY,NX)=AZMAX1((ZC3+ZA3+ZC2+ZA2+ZC1+ZA1+ZN)/VLWatMicP_vr(L,NY,NX))
      ELSE
        SolutesIonStrenth_vr(L,NY,NX)=0._r8
        SolutesIonConc_vr(L,NY,NX)=0._r8
      ENDIF
    ENDIF
!
! OSTWALD COEFFICIENTS FOR CO2, CH4, O2, N2, N2O, NH3 AND H2
! SOLUBILITY IN WATER
!

  ENDDO
  end subroutine Prep4PlantMicrobeUptake
!------------------------------------------------------------------------------------------

  subroutine GetSurfResidualProperties(I,J,NY,NX,DPTH0)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8),intent(out) :: DPTH0    !water+ice thickness in litter layer
  real(r8) :: VxcessWatLitR,TVOLWI,ThetaWLitR
  real(r8) :: VWatLitrZ
  real(r8) :: VOLIRZ
  real(r8) :: XVOLW0
  real(r8) :: XVOLI0
  integer  :: idg,NTN
! begin_execution
! PHYSICAL PROPERTIES, AND WATER, GAS, AND MINERAL CONTENTS
! OF SURFACE RESIDUE
!
! VWatLitRHoldCapcity=liter water holding capacity
! THETRX,RC0=specific WHC,mass of woody(0),fine(1),manure(2) litter
! VLitR=dry litter volume
! BulkDensLitR=dry bulk density of woody(0),fine(1),manure(2) litter
! VxcessWatLitR=excess litter water+ice
! VOLT,VLSoilPoreMicP_vr=wet litter volume
! BKVL=litter mass
! VOLW,VOLI,VOLA,VOLP=litter water,ice,porosity,air volume
! THETW,THETI,THETA,THETP=litter water,ice,porosity,air concentration
! POROS=litter porosity
! THETW0,THETI0,DPTH0=litter excess water,ice,water+ice depth
! DLYR=litter thickness
! PSISM,PSISE=litter matric,saturation water potential
!
  !write(*,*) "In GetSurfResidualProperties: "
  VxcessWatLitR          = AZMAX1(VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)-VWatLitRHoldCapcity_col(NY,NX))
  VGeomLayer_vr(0,NY,NX) = VxcessWatLitR+VLitR_col(NY,NX)
  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VLSoilPoreMicP_vr(0,NY,NX) = VGeomLayer_vr(0,NY,NX)
    VLSoilMicPMass_vr(0,NY,NX) = MWC2Soil*SoilOrgM_vr(ielmc,0,NY,NX)
    VLMicP_vr(0,NY,NX)         = AZMAX1(VLitR_col(NY,NX)-VLSoilMicPMass_vr(0,NY,NX)/1.30_r8)
    VLsoiAirP_vr(0,NY,NX)      = AZMAX1(VLMicP_vr(0,NY,NX)-VLWatMicP_vr(0,NY,NX)-VLiceMicP_vr(0,NY,NX))
    IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX))THEN
      POROS_vr(0,NY,NX)    = VLMicP_vr(0,NY,NX)/VLitR_col(NY,NX)
      THETW_vr(0,NY,NX)    = AZMAX1(AMIN1(1.0_r8,VLWatMicP_vr(0,NY,NX)/VLitR_col(NY,NX)))
      THETI_vr(0,NY,NX)    = AZMAX1(AMIN1(1.0_r8,VLiceMicP_vr(0,NY,NX)/VLitR_col(NY,NX)))
      ThetaAir_vr(0,NY,NX) = AZMAX1(AMIN1(1.0_r8,VLsoiAirP_vr(0,NY,NX)/VLitR_col(NY,NX)))
    ELSE
      POROS_vr(0,NY,NX)    = 1.0_r8
      THETW_vr(0,NY,NX)    = 0._r8
      THETI_vr(0,NY,NX)    = 0._r8
      ThetaAir_vr(0,NY,NX) = 0._r8
    ENDIF
    TVOLWI=VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)
    IF(TVOLWI.GT.ZEROS(NY,NX))THEN
      VWatLitrZ = VLWatMicP_vr(0,NY,NX)/TVOLWI*VWatLitRHoldCapcity_col(NY,NX)    !water helding capacity by litter
      VOLIRZ    = VLiceMicP_vr(0,NY,NX)/TVOLWI*VWatLitRHoldCapcity_col(NY,NX)    !ice helding capacity by litter
      XVOLW0    = AZMAX1(VLWatMicP_vr(0,NY,NX)-VWatLitrZ)/AREA(3,NU(NY,NX),NY,NX)
      XVOLI0    = AZMAX1(VLiceMicP_vr(0,NY,NX)-VOLIRZ)/AREA(3,NU(NY,NX),NY,NX)
    ELSE
      XVOLW0=0._r8
      XVOLI0=0._r8
    ENDIF
    DPTH0    = XVOLW0+XVOLI0
    DLYR_3D(3,0,NY,NX) = VLSoilPoreMicP_vr(0,NY,NX)/AREA(3,0,NY,NX)

    IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX) .AND. VLWatMicP_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VLWatMicP_vr(0,NY,NX))/VLitR_col(NY,NX)
      IF(ThetaWLitR.LT.FieldCapacity_vr(0,NY,NX))THEN
        PSISoilMatricP_vr(0,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX)+((LOGFldCapacity_vr(0,NY,NX)-LOG(ThetaWLitR)) &
          /FCD_vr(0,NY,NX)*LOGPSIMND(NY,NX))))
      ELSEIF(ThetaWLitR.LT.POROS_vr(0,NY,NX))THEN
        PSISoilMatricP_vr(0,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(0,NY,NX)-LOG(ThetaWLitR)) &
          /PSD_vr(0,NY,NX))**SRP_vr(0,NY,NX)*LOGPSIMXD(NY,NX)))
      ELSE
        PSISoilMatricP_vr(0,NY,NX)=PSISE_vr(0,NY,NX)
      ENDIF
      PSISoilOsmotic_vr(0,NY,NX)     = 0._r8
      PSIGrav_vr(0,NY,NX)            = mGravAccelerat*(ALT(NY,NX)-CumDepz2LayBottom_vr(NU(NY,NX)-1,NY,NX)+0.5_r8*DLYR_3D(3,0,NY,NX))
      ElvAdjstedSoilH2OPSIMPa_vr(0,NY,NX) = AZMIN1(PSISoilMatricP_vr(0,NY,NX)+PSISoilOsmotic_vr(0,NY,NX)+PSIGrav_vr(0,NY,NX))
!
!     LITTER NH4,NH3,NO3,NO2,HPO4,H2PO4 CONCENTRATIONS
!
!     C*=litter solute concentrations
!
      DO NTN=ids_nut_beg,ids_nuts_end
        trc_solcl_vr(NTN,0,NY,NX)=AZMAX1(trcs_solml_vr(NTN,0,NY,NX)/VLWatMicP_vr(0,NY,NX))
      ENDDO

    ELSE
      PSISoilMatricP_vr(0,NY,NX)                     = PSISoilMatricP_vr(NU(NY,NX),NY,NX)
      trc_solcl_vr(ids_nut_beg:ids_nuts_end,0,NY,NX) = 0._r8
    ENDIF
  ELSE
    VLSoilPoreMicP_vr(0,NY,NX)                     = 0._r8
    VLSoilMicPMass_vr(0,NY,NX)                     = 0._r8
    VLMicP_vr(0,NY,NX)                             = 0._r8
    VLsoiAirP_vr(0,NY,NX)                          = 0._r8
    POROS_vr(0,NY,NX)                              = 1.0_r8
    DLYR_3D(3,0,NY,NX)                             = 0._r8
    THETW_vr(0,NY,NX)                              = 0._r8
    THETI_vr(0,NY,NX)                              = 0._r8
    ThetaAir_vr(0,NY,NX)                           = 1.0_r8
    VWatLitRHoldCapcity_col(NY,NX)                 = 0._r8
    PSISoilMatricP_vr(0,NY,NX)                     = PSISoilMatricP_vr(NU(NY,NX),NY,NX)
    trc_solcl_vr(ids_nut_beg:ids_nuts_end,0,NY,NX) = 0._r8
    trc_solcl_vr(idg_beg:idg_NH3,0,NY,NX)        = 0._r8
  ENDIF
  end subroutine GetSurfResidualProperties


!------------------------------------------------------------------------------------------

  subroutine GetChemicalConcsInSoil(I,J,NY,NX,THETPZ_vr)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: NY,NX
  real(r8), intent(out) :: THETPZ_vr(JZ)  !air-filled soil pore, m3/d2
  integer :: L,idg
!     begin_execution

!     CALCULATE SOIL CONCENTRATIONS OF SOLUTES, GASES
!
!     THETW,THETI,THETP=soil micropore water,ice,air concentration
!     THETPZ=soil micropore+macropore air concn for output
!

  DO L=NUI(NY,NX),NLI(NY,NX)

    IF(VLSoilPoreMicP_vr(L,NY,NX).LE.ZEROS(NY,NX))THEN
      THETW_vr(L,NY,NX)    = POROS_vr(L,NY,NX)
      THETI_vr(L,NY,NX)    = 0._r8
      ThetaAir_vr(L,NY,NX) = 0._r8
    ELSE
      THETW_vr(L,NY,NX)    = AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLWatMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)))
      THETI_vr(L,NY,NX)    = AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLiceMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)))
      ThetaAir_vr(L,NY,NX) = AZMAX1(VLsoiAirP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX))
    ENDIF
    THETPZ_vr(L)=AZMAX1(POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX)-THETI_vr(L,NY,NX))
!
!     GAS CONCENTRATIONS
!
    IF(ThetaAir_vr(L,NY,NX).GT.AirFillPore_Min)THEN
      DO idg=idg_beg,idg_NH3
        trcg_gascl_vr(idg,L,NY,NX)=AZMAX1(trcg_gasml_vr(idg,L,NY,NX)/VLsoiAirP_vr(L,NY,NX))
      ENDDO
    ELSE
      trcg_gascl_vr(idg_beg:idg_NH3,L,NY,NX)=0._r8
    ENDIF

    IF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DO idg=idg_beg,idg_NH3
        trc_solcl_vr(idg,L,NY,NX)=AZMAX1(trcs_solml_vr(idg,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
      ENDDO
    ELSE
      trc_solcl_vr(idg_beg:idg_NH3,L,NY,NX)=0._r8
    ENDIF
!    print*,L,trcs_solml_vr(idg_CH4,L,NY,NX),VLWatMicP_vr(L,NY,NX)

!
!     CORGC=SOC concentration
!
    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CSoilOrgM_vr(ielmc,L,NY,NX)=AMIN1(orgcden,SoilOrgM_vr(ielmc,L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX))
    ELSE
      CSoilOrgM_vr(ielmc,L,NY,NX)=0._r8
    ENDIF
  ENDDO
  end subroutine GetChemicalConcsInSoil
!------------------------------------------------------------------------------------------

  subroutine ZeroHourlyArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,L,NTSA
!     begin_execution

  RootN2Fix_col(NY,NX)        = 0._r8
  RUptkRootO2_col(NY,NX)      = 0._r8
  RootCO2Emis2Root_col(NY,NX) = 0._r8
  trcs_plant_uptake_col(ids_beg:ids_end,NY,NX) = 0._r8  
  DO L=NUI(NY,NX),NLI(NY,NX)
    RootN2Fix_vr(L,NY,NX)                         = 0._r8
    FWatExMacP2MicP_vr(L,NY,NX)                   = 0._r8
    trcs_plant_uptake_vr(ids_beg:ids_end,L,NY,NX) = 0._r8
    RootCO2Emis2Root_vr(L,NY,NX)                 = 0._r8
    RUptkRootO2_vr(L,NY,NX)                       = 0._r8
    trcg_air2root_flx_vr(idg_beg:idg_NH3,L,NY,NX) = 0._r8

    trcn_RChem_band_soil_vr(ids_NH4B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(ids_NO3B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(ids_NO2B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(ids_H1PO4B,L,NY,NX) = 0._r8
    trcn_RChem_band_soil_vr(ids_H2PO4B,L,NY,NX) = 0._r8
    TProd_CO2_geochem_soil_vr(L,NY,NX)          = 0._r8
    Txchem_CO2_vr(L,NY,NX)                      = 0._r8

    trcx_TRSoilChem_vr(idx_NH4B,L,NY,NX)=0._r8
    trcx_TRSoilChem_vr(idx_OHeB:idx_end,L,NY,NX)=0._r8

    TRChem_H_p_sorbed_soil_vr(L,NY,NX)       = 0._r8
    TRChem_Al_sorbed_soil_vr(L,NY,NX)        = 0._r8
    TRChem_Fe_sorbed_soil_vr(L,NY,NX)     = 0._r8
    TRChem_Ca_sorbed_soil_vr(L,NY,NX)        = 0._r8
    TRChem_Mg_sorbed_soil_vr(L,NY,NX)        = 0._r8
    TRChem_Na_sorbed_soil_vr(L,NY,NX)        = 0._r8
    TRChem_K_sorbed_soil_vr(L,NY,NX)         = 0._r8
    TRChem_HCO3_sorbed_soil_vr(L,NY,NX)      = 0._r8
    TRChem_AlO2H2_sorbed_soil_vr(L,NY,NX)    = 0._r8
    TRChem_FeO2H2_sorbed_soil_vr(L,NY,NX) = 0._r8

    trcp_RChem_soil_vr(idsp_beg:idsp_psoi_beg-1,L,NY,NX)=0._r8

    trcp_RChem_soil_vr(idsp_beg_band:idsp_end,L,NY,NX)=0._r8

    trcs_Mac2MicPore_flx_vr(ids_beg:ids_end,L,NY,NX)=0._r8

    DO NTSA=idsalt_beg,idsaltb_end
      trcSalt_RGeoChem_flx_vr(NTSA,L,NY,NX)=0._r8
      trcSalt_Mac2MicPore_flx_vr(NTSA,L,NY,NX)=0._r8
    ENDDO

    DO  K=1,jcplx
      DOM_Mac2MicPore_flx_vr(idom_beg:idom_end,K,L,NY,NX)=0._r8
    ENDDO
    TLIceThawMicP_vr(L,NY,NX)=0._r8
    TLIceThawMacP_vr(L,NY,NX)=0._r8
    TLPhaseChangeHeat2Soi_vr(L,NY,NX)=0._r8
    trcg_ebu_flx_vr(idg_beg:idg_end,L,NY,NX)=0._r8
    totRootLenDens_vr(L,NY,NX)=0._r8

  ENDDO
  trcg_ebu_flx_col(idg_beg:idg_NH3,NY,NX)=0._r8
  trcg_air2root_flx_col(idg_beg:idg_NH3,NY,NX)=0._r8
  end subroutine ZeroHourlyArrays

!------------------------------------------------------------------------------------------

  subroutine CalGasSolubility(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer  :: L,idg
  real(r8) :: FH2O

  L=0
  DO idg=idg_beg,idg_NH3
    GasSolbility_vr(idg,L,NY,NX)=gas_solubility(idg,TCS_vr(L,NY,NX))
  ENDDO

  GasSolbility_vr(idg_NH3B,L,NY,NX)=GasSolbility_vr(idg_NH3,L,NY,NX)
 
  DO  L=1,NL(NY,NX)+1
    ! S*L=solubility of gas in water
    ! TCS=soil temperature (oC)
    ! 5.56E+04_r8 := mole H2O / m3
    FH2O=5.56E+04_r8/(5.56E+04_r8+SolutesIonConc_vr(L,NY,NX))
    DO idg=idg_beg,idg_NH3
      GasSolbility_vr(idg,L,NY,NX)=gas_solubility(idg,TCS_vr(L,NY,NX))*EXP(-GasSechenovConst(idg)*SolutesIonStrenth_vr(L,NY,NX))*FH2O
    ENDDO
    GasSolbility_vr(idg_NH3B,L,NY,NX)=GasSolbility_vr(idg_NH3,L,NY,NX)
  ENDDO
  end subroutine CalGasSolubility
!------------------------------------------------------------------------------------------
  subroutine UpdateLiterPropertz(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS  
  real(r8) :: FVLitR  !litter porosity
  integer :: NY,NX
!
!     VWatLitRHoldCapcity=liter water holding capacity
!     VLitR=dry litter volume
!     POROS0,FC,WP=litter porosity,field capacity,wilting point
!
  real(r8) :: VWatLitRHoldCapcity0,VLitR0
  integer  :: K

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      VWatLitRHoldCapcity0 = 0._r8
      VLitR0               = 0._r8
      DO K=1, micpar%NumOfLitrCmplxs
        VWatLitRHoldCapcity0 = VWatLitRHoldCapcity0+THETRX(K)*RC0(K,NY,NX)
        VLitR0               = VLitR0+RC0(K,NY,NX)/BulkDensLitR(K)
      ENDDO

      VWatLitRHoldCapcity_col(NY,NX) = AZMAX1(VWatLitRHoldCapcity0)
      VLitR_col(NY,NX)               = AZMAX1(VLitR0*ppmc)

      IF(AMIN1(VLitR_col(NY,NX),VWatLitRHoldCapcity_col(NY,NX)).GT.ZEROS(NY,NX))THEN
        FVLitR=VWatLitRHoldCapcity_col(NY,NX)/VLitR_col(NY,NX)
      ELSE
        FVLitR=THETRX(micpar%k_fine_litr)/BulkDensLitR(micpar%k_fine_litr)
      ENDIF      
      POROS0_col(NY,NX)          = FVLitR      
      FieldCapacity_vr(0,NY,NX)  = 0.500_r8*FVLitR
      WiltPoint_vr(0,NY,NX)      = 0.125_r8*FVLitR      
      LOGPOROS_vr(0,NY,NX)       = LOG(POROS0_col(NY,NX))
      LOGFldCapacity_vr(0,NY,NX) = LOG(FieldCapacity_vr(0,NY,NX))
      LOGWiltPoint_vr(0,NY,NX)   = LOG(WiltPoint_vr(0,NY,NX))
      PSD_vr(0,NY,NX)            = LOGPOROS_vr(0,NY,NX)-LOGFldCapacity_vr(0,NY,NX)
      FCD_vr(0,NY,NX)            = LOGFldCapacity_vr(0,NY,NX)-LOGWiltPoint_vr(0,NY,NX)
      SRP_vr(0,NY,NX)            = 1.00_r8
    enddo
  enddo    
  end subroutine UpdateLiterPropertz


end module Hour1Mod
