module Hour1Mod
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: GravAcceleration=>DAT_CONST_G
  use minimathmod,    only: isclose, AZMAX1, AZMIN1
  use abortutils,     only: endrun,  print_info
  use fileUtil,       only: iulog
  use PlantMod,       only: PlantCanopyRadsModel
  use EcoSIMConfig, only : jcplx=>jcplxc,nlbiomcp=>NumLiveMicrbCompts
  use EcoSIMConfig, only : ndbiomcp=>NumDeadMicrbCompts,jsken=>jskenc
  use EcoSIMConfig, only : NumMicbFunGrupsPerCmplx=>NumMicbFunGrupsPerCmplx,do_instequil
  use EcoSiMParDataMod, only : micpar,pltpar
  use SoilBGCNLayMod, only : sumORGMLayL
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
  integer :: L,NX,NY
  real(r8) :: THETPZ_vr(JZ,JY,JX)
  real(r8) :: DPTH0(JY,JX)

  real(r8) :: tPBOT,tmp
  integer :: NZ,NR,K
!     execution begins here
!  write(111,*)'xxxxxx',I,J,CanopyLeafArea_lpft(1,25,1,1,1,5)
!
  if(lverb)write(*,*)'ResetLndscapeAccumlators'
  call ResetLndscapeAccumlators()

  if(lverb)write(*,*)'SetAtmosTracerConcentration'
  call SetAtmosTracerConcentration(I,NHW,NHE,NVN,NVS)
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

      tPBOT                        = PBOT_col(NY,NX)/1.01325E+02_r8
      tmp                          = Tref/TairKClimMean(NY,NX)*tPBOT
      CCO2EI(NY,NX)                = CO2EI(NY,NX)*5.36E-04_r8*tmp
      AtmGasCgperm3(idg_CO2,NY,NX) = CO2E_col(NY,NX)*5.36E-04_r8*tmp
      AtmGasCgperm3(idg_CH4,NY,NX) = CH4E_col(NY,NX)*5.36E-04_r8*tmp
      AtmGasCgperm3(idg_O2,NY,NX)  = OXYE(NY,NX)*1.43E-03_r8*tmp
      AtmGasCgperm3(idg_N2,NY,NX)  = Z2GE(NY,NX)*1.25E-03_r8*tmp
      AtmGasCgperm3(idg_N2O,NY,NX) = Z2OE(NY,NX)*1.25E-03_r8*tmp
      AtmGasCgperm3(idg_NH3,NY,NX) = ZNH3E_col(NY,NX)*6.25E-04_r8*tmp
      AtmGasCgperm3(idg_H2,NY,NX)  = H2GE(NY,NX)*8.92E-05_r8*tmp

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
!
!     RESET SURFACE LITTER PHYSICAL PROPERTIES (DENSITY, TEXTURE)
!     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
  if(lverb)write(*,*)'SetLiterSoilPropAftDisturb'
  call SetLiterSoilPropAftDisturb(I,J,NHW,NHE,NVN,NVS)

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
      call SetHourlyAccumulators(NY,NX)
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

      call GetChemicalConcsInSoil(NY,NX,THETPZ_vr)

      call GetSoluteConcentrations(NY,NX)

      call Prep4PlantMicrobeUptake(NY,NX)

      call CalGasSolubility(NY,NX)

      call GetSoilHydraulicVars(NY,NX)

!     CALCULATE ACTIVE LAYER DEPTH
      call DiagActiveLayerDepth(NY,NX)

!     OUTPUT FOR WATER TABLE DEPTH
      call GetOutput4WaterTableDepth(NY,NX,THETPZ_vr)

      call GetSurfResidualProperties(I,J,NY,NX,DPTH0)

      call SetTracerPropertyInLiterAir(NY,NX)

      if(do_instequil)call ForceGasAquaEquil(NY,NX)
!
      call PlantCanopyRadsModel(I,J,NY,NX,DPTH0(NY,NX))
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
      HeatPrec_col(NY,NX)              = 0._r8

      DO  NZ=1,NP(NY,NX)
!
!     NUMBERS OF TOP AND BOTTOM ROOTED SOIL LAYERS
!
!     NG=number of uppermost rooted layer
!     NIXBotRootLayer_rpft=number of lowest rooted layer
!
        NGTopRootLayer_pft(NZ,NY,NX) =MAX(NGTopRootLayer_pft(NZ,NY,NX),NU(NY,NX))
        NIXBotRootLayer_pft(NZ,NY,NX)=MAX(NIXBotRootLayer_pft(NZ,NY,NX),NU(NY,NX))
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

  call SummarizeTracers(NHW,NHE,NVN,NVS)

  END subroutine hour1
!------------------------------------------------------------------------------------------

  subroutine CanopyInterceptPrecp(NY,NX)
  !
  !DESCRIPTION
  !precipitation intercepation by canopy
  implicit none
  integer, intent(in) :: NY,NX
  integer :: NZ
  real(r8) :: VOLWCX  !maximum precipitation holding capacity by canopy (leaf+stem) [m H2O/h]
  real(r8) :: prec2canopy_pft  !precipiation onto canopy [m H2O/h]
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
  DO  NZ=1,NP(NY,NX)
    VOLWCX                           = FoliarWatRetcap(iPlantRootProfile_pft(NZ,NY,NX)) &
      *(CanopyLeafArea_pft(NZ,NY,NX)+CanopyStemArea_pft(NZ,NY,NX))
    prec2canopy_pft                  = PrecRainAndIrrig_col(NY,NX)*FracPARads2Canopy_pft(NZ,NY,NX)
    PrecIntcptByCanopy_pft(NZ,NY,NX) = AZMAX1(AMIN1(prec2canopy_pft,VOLWCX-WatByPCanopy_pft(NZ,NY,NX)))
    Prec2Canopy_col(NY,NX)           = Prec2Canopy_col(NY,NX)+prec2canopy_pft
    PrecIntceptByCanopy_col(NY,NX)   = PrecIntceptByCanopy_col(NY,NX)+PrecIntcptByCanopy_pft(NZ,NY,NX)
  ENDDO

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

      AtmGasCgperm3(idg_CO2,NY,NX)  = CO2E_col(NY,NX)*5.36E-04_r8*TREF/TairK_col(NY,NX)  !gC/m3
      AtmGasCgperm3(idg_CH4,NY,NX)  = CH4E_col(NY,NX)*5.36E-04_r8*TREF/TairK_col(NY,NX)  !gC/m3
      AtmGasCgperm3(idg_O2 ,NY,NX)  = OXYE(NY,NX)*1.43E-03_r8*TREF/TairK_col(NY,NX)  !gO/m3
      AtmGasCgperm3(idg_N2 ,NY,NX)  = Z2GE(NY,NX)*1.25E-03_r8*TREF/TairK_col(NY,NX)  !gN/m3
      AtmGasCgperm3(idg_N2O,NY,NX)  = Z2OE(NY,NX)*1.25E-03_r8*TREF/TairK_col(NY,NX)  !gN/m3
      AtmGasCgperm3(idg_NH3,NY,NX)  = ZNH3E_col(NY,NX)*6.25E-04_r8*TREF/TairK_col(NY,NX) !gN/m3
      AtmGasCgperm3(idg_H2 ,NY,NX)  = H2GE(NY,NX)*8.92E-05_r8*TREF/TairK_col(NY,NX)  !gH/m3
      AtmGasCgperm3(idg_NH3B,NY,NX) = ZNH3E_col(NY,NX)*6.25E-04_r8*TREF/TairK_col(NY,NX) !gN/m3

      trcVolatile_rain_conc(idg_CO2,NY,NX) = AtmGasCgperm3(idg_CO2,NY,NX)*gas_solubility(idg_CO2,TCA_col(NY,NX)) &
         /(EXP(ACTCG(idg_CO2)*CSTRR(NY,NX)))
      trcVolatile_rain_conc(idg_CH4,NY,NX) = AtmGasCgperm3(idg_CH4,NY,NX)*gas_solubility(idg_CH4,TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_CH4)*CSTRR(NY,NX)))
      trcVolatile_rain_conc(idg_O2,NY,NX) = AtmGasCgperm3(idg_O2,NY,NX)*gas_solubility(idg_O2, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_O2)*CSTRR(NY,NX)))
      trcVolatile_rain_conc(idg_N2,NY,NX) = AtmGasCgperm3(idg_N2,NY,NX)*gas_solubility(idg_N2, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_N2)*CSTRR(NY,NX)))
      trcVolatile_rain_conc(idg_N2O,NY,NX) = AtmGasCgperm3(idg_N2O,NY,NX)*gas_solubility(idg_N2O, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_N2O)*CSTRR(NY,NX)))

      trcVolatile_irrig_conc(idg_CO2,NY,NX) = AtmGasCgperm3(idg_CO2,NY,NX)*gas_solubility(idg_CO2, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_CO2)*CSTRQ(I,NY,NX)))
      trcVolatile_irrig_conc(idg_CH4,NY,NX) = AtmGasCgperm3(idg_CH4,NY,NX)*gas_solubility(idg_CH4, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_CH4)*CSTRQ(I,NY,NX)))
      trcVolatile_irrig_conc(idg_O2,NY,NX) = AtmGasCgperm3(idg_O2,NY,NX)*gas_solubility(idg_O2, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_O2)*CSTRQ(I,NY,NX)))
      trcVolatile_irrig_conc(idg_N2,NY,NX) = AtmGasCgperm3(idg_N2,NY,NX)*gas_solubility(idg_N2, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_N2)*CSTRQ(I,NY,NX)))
      trcVolatile_irrig_conc(idg_N2O,NY,NX) = AtmGasCgperm3(idg_N2O,NY,NX)*gas_solubility(idg_N2O, TCA_col(NY,NX)) &
        /(EXP(ACTCG(idg_N2O)*CSTRQ(I,NY,NX)))
      GDD_col(NY,NX) = GDD_col(NY,NX)+TCA_col(NY,NX)/24._r8
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

  XGridSurfRunoff_2DH(1:2,1:2,:,:)       = 0._r8
  HeatXGridBySurfRunoff_2DH(1:2,1:2,:,:) = 0._r8

  DO  NX=NHW,NHE+extragrid
    DO  NY=NVN,NVS+extragrid
!
!     WATER,SNOW,SOLUTE RUNOFF
!
      HydroSufDOCFlx_col(NY,NX)            = 0._r8
      HydroSubsDOCFlx_col(NY,NX)           = 0._r8
      HydroSufDICFlx_col(NY,NX)            = 0._r8
      HydroSubsDICFlx_col(NY,NX)           = 0._r8
      HydroSubsDONFlx_col(NY,NX)           = 0._r8
      HydroSubsDINFlx_col(NY,NX)           = 0._r8
      HydroSubsDOPFlx_col(NY,NX)           = 0._r8
      HydroSubsDIPFlx_col(NY,NX)           = 0._r8
      WatFlux4ErosionM_2DH(:,NY,NX)        = 0._r8

      dom_2DFloXSurRunoff(idom_beg:idom_end,1:jcplx,1:2,1:2,NY,NX)=0._r8

      trcg_FloXSurRunoff_2D(idg_beg:idg_end-1,1:2,1:2,NY,NX)=0._r8
      trcn_FloXSurRunoff_2D(ids_nut_beg:ids_nuts_end,1:2,1:2,NY,NX)=0._r8

      DrysnoBySnowRedistrib(1:2,NY,NX)             = 0._r8
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
      DO L=1,NL(NY,NX)+1
        WaterFlowSoiMicP_3D(1:3,L,NY,NX) = 0._r8
        WaterFlowSoiMicPX(1:3,L,NY,NX)   = 0._r8
        WaterFlowMacP_3D(1:3,L,NY,NX)    = 0._r8
        HeatFlow2Soil_3D(1:3,L,NY,NX)    = 0._r8

        trcs_TransptMicP_3D(ids_beg:ids_end,1:3,L,NY,NX)=0._r8
        Gas_3DAdvDif_Flx_vr(idg_beg:idg_end,1:3,L,NY,NX)=0._r8

        DOM_3DMacp_Transp_flx(idom_beg:idom_end,1:jcplx,1:3,L,NY,NX)=0._r8

      ENDDO
    ENDDO
  ENDDO

!     IF EROSION FLAG SET
!
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    DO NX=NHW,NHE+extragrid
      DO NY=NVN,NVS+extragrid
        cumSedErosion(1:2,1:2,NY,NX)=0._r8
        XSANER(1:2,1:2,NY,NX)=0._r8
        XSILER(1:2,1:2,NY,NX)=0._r8
        XCLAER(1:2,1:2,NY,NX)=0._r8
        XNH4ER(1:2,1:2,NY,NX)=0._r8
        XNH3ER(1:2,1:2,NY,NX)=0._r8
        XNHUER(1:2,1:2,NY,NX)=0._r8
        XNO3ER(1:2,1:2,NY,NX)=0._r8
        XNH4EB(1:2,1:2,NY,NX)=0._r8
        XNH3EB(1:2,1:2,NY,NX)=0._r8
        XNHUEB(1:2,1:2,NY,NX)=0._r8
        XNO3EB(1:2,1:2,NY,NX)=0._r8

        trcx_XER(idx_beg:idx_end,1:2,1:2,NY,NX)=0._r8
        trcp_ER(idsp_beg:idsp_end,1:2,1:2,NY,NX)=0._r8

        OMEERhetr(:,:,:,1:2,1:2,NY,NX)=0._r8
        OMEERauto(:,:,1:2,1:2,NY,NX)=0._r8

        ORMER(1:NumPlantChemElms,:,:,1:2,1:2,NY,NX)=0._r8
        OSAER(:,:,1:2,1:2,NY,NX)=0._r8
        OHMER(1:NumPlantChemElms,:,1:2,1:2,NY,NX)=0._r8
        OSMER(1:NumPlantChemElms,:,:,1:2,1:2,NY,NX)=0._r8
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

      trc_salt_rof_bounds(idsalt_beg:idsalt_end,1:2,1:2,NY,NX)=0._r8
      trcSalt_XQS(idsalt_beg:idsalt_end,1:2,NY,NX)=0._r8

      DO  L=1,NL(NY,NX)+1
        DO NSA=idsalt_beg,idsaltb_end
          trcSalt3DFlo2Cell(NSA,1:3,L,NY,NX)=0._r8
          trcSalt_XFHS(NSA,1:3,L,NY,NX)=0._r8
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
    !
    IF(SoiBulkDensity_vr(L,NY,NX).LE.ZERO.AND.DLYR(3,L,NY,NX).LE.ZERO2)THEN
      VLWatMicP_vr(L,NY,NX)=0._r8
      VLiceMicP_vr(L,NY,NX)=0._r8
    ENDIF
    AREA(1,L,NY,NX)        = DLYR(3,L,NY,NX)*DLYR(2,L,NY,NX)
    AREA(2,L,NY,NX)        = DLYR(3,L,NY,NX)*DLYR(1,L,NY,NX)
    VGeomLayer_vr(L,NY,NX) = AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)

    VLSoilPoreMicP_vr(L,NY,NX) = AMAX1(VGeomLayer_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX),1.e-8_r8)
    IF(SoiBulkDensity_vr(L,NY,NX).LE.ZERO)THEN
      VLSoilMicP_vr(L,NY,NX)=VLSoilPoreMicP_vr(L,NY,NX)
    ENDIF
    !
    !     BKVL=soil mass
    !     C*=concentration,ORGC=SOC,SAND=sand,SILT=silt,CLAY=clay
    !     ParticleDens=particle density
    !     ParticleDensitySurfLay=particle density of surface layer for use in erosion.f
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
    VLSoilMicPMass_vr(L,NY,NX)=SoiBulkDensity_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)

    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CSoilOrgM_vr(ielmc,L,NY,NX) = AMIN1(orgcden,SoilOrgM_vr(ielmc,L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX))
      CSAND(L,NY,NX)              = SAND(L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
      CSILT(L,NY,NX)              = SILT(L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
      CCLAY(L,NY,NX)              = CLAY(L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
    ELSE
      CSoilOrgM_vr(ielmc,L,NY,NX)=0._r8
      CSAND(L,NY,NX)=0._r8
      CSILT(L,NY,NX)=0._r8
      CCLAY(L,NY,NX)=0._r8
    ENDIF
    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZERO)THEN
      CORGCM=AZMAX1(AMIN1(1.0_r8,MWC2Soil*CSoilOrgM_vr(ielmc,L,NY,NX)))
      ParticleDens=1.30_r8*CORGCM+2.66_r8*(1.0_r8-CORGCM)
      IF(L.EQ.NU(NY,NX))THEN
!surface layer
        POROS_vr(L,NY,NX)=AMAX1(POROS_vr(L,NY,NX),1.0_r8-(SoiBulkDensity_vr(L,NY,NX)/ParticleDens))
      ELSE
!  avoid float exception
        if(SoiBulkDensity_vr(L,NY,NX)>=ParticleDens)then
          POROS_vr(L,NY,NX)=1.0_r8
        else
          POROS_vr(L,NY,NX)=1.0_r8-(SoiBulkDensity_vr(L,NY,NX)/ParticleDens)
        endif
      ENDIF
    ELSE
      ParticleDens=0._r8
      POROS_vr(L,NY,NX)=1.0_r8
    ENDIF
    !     VLMicP_vr(L,NY,NX)=AMAX1(POROS_vr(L,NY,NX)*VLSoilMicP_vr(L,NY,NX)
    !    2,VLWatMicP_vr(L,NY,NX)+VLiceMicP_vr(L,NY,NX))
    !     VLMacP_vr(L,NY,NX)=AMAX1(SoilFracAsMacP_vr(L,NY,NX)*VGeomLayer_vr(L,NY,NX)
    !    2,VLWatMacP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
    VLMicP_vr(L,NY,NX) = POROS_vr(L,NY,NX)*VLSoilMicP_vr(L,NY,NX)
    VLMacP_vr(L,NY,NX) = SoilFracAsMacP_vr(L,NY,NX)*VGeomLayer_vr(L,NY,NX)
    IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      VLsoiAirP_vr(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)) &
        +AZMAX1(VLMacP_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX))
    ELSE
      VLsoiAirP_vr(L,NY,NX)=0._r8
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
    if(lverb)write(*,*)'setthermcond'
    IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      VORGC=CORGCM*SoiBulkDensity_vr(L,NY,NX)/ParticleDens
      VMINL=(CSILT(L,NY,NX)+CCLAY(L,NY,NX))*SoiBulkDensity_vr(L,NY,NX)/ParticleDens
      VSAND=CSAND(L,NY,NX)*SoiBulkDensity_vr(L,NY,NX)/ParticleDens
      NumerSolidThermCond(L,NY,NX)=(1.253_r8*VORGC*9.050E-04_r8+0.514_r8*VMINL*1.056E-02_r8 &
        +0.386_r8*VSAND*2.112E-02_r8)*FracSoiAsMicP_vr(L,NY,NX) &
        +0.514_r8*ROCK_vr(L,NY,NX)*1.056E-02_r8
      DenomSolidThermCond(L,NY,NX)=(1.253_r8*VORGC+0.514_r8*VMINL+0.386_r8*VSAND) &
        *FracSoiAsMicP_vr(L,NY,NX)+0.514_r8*ROCK_vr(L,NY,NX)
    ELSE
      NumerSolidThermCond(L,NY,NX)=0._r8
      DenomSolidThermCond(L,NY,NX)=0._r8
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
!     ExtWaterTable,ExtWaterTablet0=current,initial natural water table depth
!     DTBLY,DTBLD=current,initial artificial water table depth
!     SoilSurfRoughnesst0_col,ZW=soil,water surface roughness
!     VLWatheldCapSurf_col=soil surface water retention capacity
!     VWatStoreCapSurf_col=VLWatheldCapSurf_col accounting for above-ground water table
!     EHUM=fraction of microbial decompn product allocated to surface humus
!     EPOC=fraction of SOC decomposition product allocated to surface POC
!
  IF(IDWaterTable(NY,NX).LE.1 .OR. IDWaterTable(NY,NX).EQ.3)THEN
    ExtWaterTable_col(NY,NX)=ExtWaterTablet0(NY,NX)
  ELSEIF(IDWaterTable(NY,NX).EQ.2.OR.IDWaterTable(NY,NX).EQ.4)THEN
    ExtWaterTable_col(NY,NX)=ExtWaterTablet0(NY,NX)+CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
  ENDIF

  IF(IDWaterTable(NY,NX).EQ.3.OR.IDWaterTable(NY,NX).EQ.4)THEN
    DTBLY(NY,NX)=DTBLD(NY,NX)
  ENDIF

  IF(SoiBulkDensity_vr(NU(NY,NX),NY,NX).GT.ZERO)THEN
    SoilSurfRoughnesst0_col(NY,NX)=0.020_r8
  ELSE
    SoilSurfRoughnesst0_col(NY,NX)=ZW
  ENDIF
  VLWatheldCapSurf_col(NY,NX)=AMAX1(0.001_r8,0.112_r8*SoilSurfRoughnesst0_col(NY,NX)+&
    3.10_r8*SoilSurfRoughnesst0_col(NY,NX)**2._r8 &
    -0.012_r8*SoilSurfRoughnesst0_col(NY,NX)*SLOPE(0,NY,NX))*AREA(3,NU(NY,NX),NY,NX)

  VWatStoreCapSurf_col(NY,NX)=AMAX1(VLWatheldCapSurf_col(NY,NX),-(ExtWaterTable_col(NY,NX)-&
    CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX))*AREA(3,NU(NY,NX),NY,NX))

  SoiDepthMidLay_vr(NU(NY,NX),NY,NX)=CumDepz2LayerBot_vr(NU(NY,NX),NY,NX)-0.5_r8*DLYR(3,NU(NY,NX),NY,NX)
  IF(VLSoilMicPMass_vr(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)/VLSoilMicPMass_vr(NU(NY,NX),NY,NX)
    CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)/VLSoilMicPMass_vr(NU(NY,NX),NY,NX)
    CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)/VLSoilMicPMass_vr(NU(NY,NX),NY,NX)
  ELSE
    CCLAY(NU(NY,NX),NY,NX)=0._r8
    CSILT(NU(NY,NX),NY,NX)=0._r8
    CSAND(NU(NY,NX),NY,NX)=0._r8
  ENDIF
  EHUM(0,NY,NX)=0.200_r8+0.333_r8*AMIN1(0.5_r8,CCLAY(NU(NY,NX),NY,NX))
  EPOC(0,NY,NX)=0.150_r8
  end subroutine ResetSurfResidualProperty
!------------------------------------------------------------------------------------------

  subroutine SetLiterSoilPropAftDisturb(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX

  real(r8) :: PSISK(0:100),THETK(100)
  REAL(R8) :: SUM2,SUM1
  real(r8) :: XK,YK
  integer :: K,M
!     begin_execution
!     iResetSoilProf_col=disturbance flag
!     SoiBulkDensity_vr,BKDSI=current,initial bulk density
!

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
        iResetSoilProf_col(NY,NX)=ifalse
      ENDIF
    ENDDO  
  ENDDO
  end subroutine SetLiterSoilPropAftDisturb
!------------------------------------------------------------------------------------------

  subroutine SetHourlyAccumulators(NY,NX)
!     implicit none
  integer, intent(in) :: NX,NY

  integer :: L
!     begin_execution

  ECO_HR_CO2_col(NY,NX)                   = 0._r8
  ECO_HR_CH4_col(NY,NX)                   = 0._r8
  Eco_RadSW_col(NY,NX)                    = 0._r8
  RootCO2Autor_vr(:,NY,NX)                = 0._r8
  tRDIM2DOM_col(1:NumPlantChemElms,NY,NX) = 0._r8
  QRunSurf_col(NY,NX)                     = 0._r8
  HeatRunSurf_col(NY,NX)                  = 0._r8
  tHeatUptk_col(NY,NX)                    = 0._r8
  Qinflx2Soil_col(NY,NX)                  = 0._r8
  HeatFlx2Grnd_col(NY,NX)                 = 0._r8
  DIC_mass_col(NY,NX)                     = 0._r8
  tMicBiome_col(1:NumPlantChemElms,NY,NX) = 0._r8
  tSoilOrgM_col(1:NumPlantChemElms,NY,NX) = 0._r8
  WatMass_col(NY,NX)                      = 0._r8
  HeatStore_col(NY,NX)                    = 0._r8
  tLitrOM_col(1:NumPlantChemElms,NY,NX)   = 0._r8
  tHumOM_col(1:NumPlantChemElms,NY,NX)    = 0._r8
  tNH4_col(NY,NX)                         = 0._r8
  tNO3_col(NY,NX)                         = 0._r8
  tHxPO4_col(NY,NX)                       = 0._r8
  tXPO4_col(NY,NX)                        = 0._r8
  UION(NY,NX)                             = 0._r8
  QDischar_col(NY,NX)                     = 0._r8
  SurfGasFlx_col(idg_beg:idg_NH3,NY,NX)   = 0._r8
  WatFLo2LitR_col(NY,NX)                      = 0._r8
  HeatFLo2LitrByWat(NY,NX)                = 0._r8
  TLitrIceFlxThaw_col(NY,NX)              = 0._r8
  TLitrIceHeatFlxFrez_col(NY,NX)          = 0._r8
  HeatByRad2Surf_col(NY,NX)               = 0._r8
  HeatSensAir2Surf_col(NY,NX)             = 0._r8
  HeatEvapAir2Surf_col(NY,NX)             = 0._r8
  HeatSensVapAir2Surf_col(NY,NX)          = 0._r8
  HeatNet2Surf_col(NY,NX)                 = 0._r8
  VapXAir2GSurf_col(NY,NX)                = 0._r8

  trcs_TransptMacP_3D(:,:,:,:,:) = 0._r8
  Gas_Flx_atmDif2soil_col(idg_beg:idg_end,NY,NX) = 0._r8
  trcg_surf_disevap_flx(idg_beg:idg_end-1,NY,NX) = 0._r8

  CanWat_col(NY,NX)              = 0._r8
  CanH2OHeldVg_col(NY,NX)        = 0._r8
  Prec2Canopy_col(NY,NX)         = 0._r8
  PrecIntceptByCanopy_col(NY,NX) = 0._r8
  QvET_col(NY,NX)                = 0._r8
  VapXAir2Canopy_col(NY,NX)      = 0._r8
  HeatFlx2Canopy_col(NY,NX)      = 0._r8
  CanopyHeatStor_col(NY,NX)      = 0._r8

  TRootGasLossDisturb_pft(idg_beg:idg_end-1,NY,NX)    = 0._r8
  LitrFallStrutElms_col(:,NY,NX)                      = 0._r8
  StandingDeadStrutElms_col(1:NumPlantChemElms,NY,NX) = 0._r8
  PlantPopu_col(NY,NX)                                = 0._r8
! zero arrays in the snow layers
  WatConvSno2MicP_snvr(1:JS,NY,NX)                         = 0._r8
  WatConvSno2MacP_snvr(1:JS,NY,NX)                         = 0._r8
  HeatConvSno2Soi_snvr(1:JS,NY,NX)                         = 0._r8
  WatConvSno2LitR_snvr(1:JS,NY,NX)                         = 0._r8
  HeatConvSno2LitR_snvr(1:JS,NY,NX)                        = 0._r8
  SnoXfer2SnoLay_snvr(1:JS,NY,NX)                          = 0._r8
  WatXfer2SnoLay_snvr(1:JS,NY,NX)                          = 0._r8
  IceXfer2SnoLay_snvr(1:JS,NY,NX)                          = 0._r8
  HeatXfer2SnoLay_snvr(1:JS,NY,NX)                         = 0._r8
  XPhaseChangeHeatL_snvr(1:JS,NY,NX)                       = 0._r8
  HeatSource_vr(:,NY,NX)                                   = 0._r8
  trcVolatile_Xbndl_flx_snvr(idg_beg:idg_end-1,1:JS,NY,NX) = 0._r8
  trcn_Xbndl_flx(ids_nut_beg:ids_nuts_end,1:JS,NY,NX)      = 0._r8
  IF(salt_model)THEN
    trcSaltFlo2SnowLay(idsalt_beg:idsalt_end,1:JS,NY,NX)=0._r8
  ENDIF
  end subroutine SetHourlyAccumulators
!------------------------------------------------------------------------------------------

  subroutine SetArrays4PlantSoilTransfer(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

!     begin_execution

  LitrfalStrutElms_vr(1:NumPlantChemElms,1:jsken,1:pltpar%NumOfPlantLitrCmplxs,0:NL(NY,NX),NY,NX) = 0._r8

  REcoDOMProd_vr(idom_beg:idom_end,1:jcplx,0:NL(NY,NX),NY,NX)                 = 0._r8
  XZHYS(0:NL(NY,NX),NY,NX)                                                    = 0._r8
  trcn_RChem_soil_vr(ids_nut_beg:ids_nuts_end,0:NL(NY,NX),NY,NX)              = 0._r8
  TR_NH3_soil_vr(0:NL(NY,NX),NY,NX)                                           = 0._r8
  TR_NH3_geochem_vr(0:NL(NY,NX),NY,NX)                                        = 0._r8
  trcx_TRSoilChem_vr(idx_beg:idx_end,0:NL(NY,NX),NY,NX)                       = 0._r8
  trcp_RChem_soil(idsp_psoi_beg:idsp_psoi_end,0:NL(NY,NX),NY,NX)              = 0._r8
  TPlantRootH2OUptake_vr(0:NL(NY,NX),NY,NX)                                   = 0._r8
  THeatRootUptake_vr(0:NL(NY,NX),NY,NX)                                       = 0._r8
  Gas_Disol_Flx_vr(idg_beg:idg_end,0:NL(NY,NX),NY,NX)                         = 0._r8
  tRootMycoExud2Soil_vr(1:NumPlantChemElms,1:jcplx,NU(NY,NX):NL(NY,NX),NY,NX) = 0._r8
  RO2UptkSoilM_vr(1:NPH,NU(NY,NX):NL(NY,NX),NY,NX)                            = 0._r8
  end subroutine SetArrays4PlantSoilTransfer
!------------------------------------------------------------------------------------------

  subroutine GetOutput4WaterTableDepth(NY,NX,THETPZ_vr)
  implicit none
  integer, intent(in) :: NY,NX

  real(r8), intent(in) :: THETPZ_vr(JZ,JY,JX)
  real(r8) :: PSIEquil    !equilibrium matric potential
  real(r8) :: THETW1
  real(r8) :: THETWM
  real(r8) :: THETPX
  real(r8) :: THETPW !     THETPW=minimum air-filled porosity for saturation (m3 m-3)  
  real(r8) :: THETWP  
  integer :: LL,L
  logical :: FoundWaterTable
!     begin_execution
  THETPW=0.01_r8
  THETWP=1.0_r8-THETPW

  FoundWaterTable=.false.
  DO L=NUI(NY,NX),NLI(NY,NX)
!     IDWaterTable=water table flag from site file
!     THETPZ,THETPW=current,minimum air-filled, porosity for water table
!     DPTH,ExtWaterTable=depth of soil layer midpoint, water table
!     PSIEquil=water potential in hydraulic equilibrium with layer below
!     THETW1,THETWP=water content at PSIEquil,minimum SWC for water table
!     DepthInternalWTBL=water table depth
!
    IF(IDWaterTable(NY,NX).NE.0)THEN
      IF(.not.FoundWaterTable)THEN
        IF(THETPZ_vr(L,NY,NX).LT.THETPW.OR.L.EQ.NL(NY,NX))THEN
          FoundWaterTable=.true.
          IF(SoiDepthMidLay_vr(L,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN
            D5705: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
              IF(THETPZ_vr(LL,NY,NX).GE.THETPW.AND.LL.NE.NL(NY,NX))THEN
                !air-filled pore greater minimum, i.e. not saturated
                FoundWaterTable=.false.
                exit
              ELSE IF(SoiDepthMidLay_vr(LL,NY,NX).GE.ExtWaterTable_col(NY,NX))THEN
                !current layer is lower than external water table
                exit
              ENDIF
            END DO D5705
          ENDIF

          !THETPW=saturation criterion for water table identification
          IF(FoundWaterTable)THEN
            IF(THETPZ_vr(L,NY,NX).GE.THETPW.AND.L.NE.NL(NY,NX))THEN
              !not bottom layer, saturated
              !PSIeqv
              PSIEquil=PSISoilMatricP_vr(L+1,NY,NX)-mGravAccelerat*(SoiDepthMidLay_vr(L+1,NY,NX)-SoiDepthMidLay_vr(L,NY,NX))

              THETWM=THETWP*POROS_vr(L,NY,NX)
              THETW1=AMIN1(THETWM,EXP((LOGPSIAtSat(NY,NX)-LOG(-PSIEquil)) &
                *PSD(L,NY,NX)/LOGPSIMXD(NY,NX)+LOGPOROS_vr(L,NY,NX)))

              IF(THETWM.GT.THETW1)THEN
                THETPX=AMIN1(1.0_r8,AZMAX1((THETWM-THETW_vr(L,NY,NX))/(THETWM-THETW1)))
                DepthInternalWTBL(NY,NX)=CumDepz2LayerBot_vr(L,NY,NX)-DLYR(3,L,NY,NX)*(1.0_r8-THETPX)
              ELSE
                DepthInternalWTBL(NY,NX)=CumDepz2LayerBot_vr(L,NY,NX)-DLYR(3,L,NY,NX)
              ENDIF
            ELSEIF(L.GT.NU(NY,NX))THEN
              !not bottom layer, and not topsoil layer, partially saturated
              PSIEquil = PSISoilMatricP_vr(L,NY,NX)-mGravAccelerat*(SoiDepthMidLay_vr(L,NY,NX)-SoiDepthMidLay_vr(L-1,NY,NX))
              THETWM = THETWP*POROS_vr(L-1,NY,NX)
              THETW1 = AMIN1(THETWM,EXP((LOGPSIAtSat(NY,NX)-LOG(-PSIEquil)) &
                *PSD(L-1,NY,NX)/LOGPSIMXD(NY,NX)+LOGPOROS_vr(L-1,NY,NX)))
              IF(THETWM.GT.THETW1)THEN
                THETPX=AMIN1(1.0_r8,AZMAX1((THETWM-THETW_vr(L-1,NY,NX))/(THETWM-THETW1)))
                DepthInternalWTBL(NY,NX)=CumDepz2LayerBot_vr(L-1,NY,NX)-DLYR(3,L-1,NY,NX)*(1.0_r8-THETPX)
              ELSE
                DepthInternalWTBL(NY,NX)=CumDepz2LayerBot_vr(L-1,NY,NX)-DLYR(3,L-1,NY,NX)
              ENDIF
            ELSE
              DepthInternalWTBL(NY,NX)=CumDepz2LayerBot_vr(L,NY,NX)-DLYR(3,L,NY,NX)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  END DO
  end subroutine GetOutput4WaterTableDepth
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
        CORGM=MWC2Soil*SoilOrgM_vr(ielmc,NU(NY,NX),NY,NX)/BKVLNX
        CSoilOrgM_vr(ielmc,NU(NY,NX),NY,NX)=orgcden*CORGM
        CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)/BKVLNX
        CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)/BKVLNX
        CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)/BKVLNX
      ELSE
        CORGM=0._r8
        CSoilOrgM_vr(ielmc,NU(NY,NX),NY,NX)=0._r8
        CSAND(NU(NY,NX),NY,NX)=0._r8
        CSILT(NU(NY,NX),NY,NX)=1.0
        CCLAY(NU(NY,NX),NY,NX)=0._r8
      ENDIF
      
      IF(iErosionMode.EQ.ieros_frzthawsom.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
        D50=1.0_r8*CCLAY(NU(NY,NX),NY,NX)+10._r8*CSILT(NU(NY,NX),NY,NX) &
          +100._r8*CSAND(NU(NY,NX),NY,NX)+100._r8*CORGM
        ZD50=0.041*(ppmc*D50)**0.167_r8
        SoiSurfRoughness(NY,NX)=SoilSurfRoughnesst0_col(NY,NX)+ZD50+1.0_r8*VLitR_col(NY,NX)/AREA(3,0,NY,NX)
        CER(NY,NX)=((D50+5.0_r8)/0.32_r8)**(-0.6_r8)
        XER(NY,NX)=((D50+5.0_r8)/300._r8)**0.25_r8

        SoilDetachability4Erosion1(NY,NX)=ppmc*(1.0_r8+2.0_r8*(1.0_r8-CSILT(NU(NY,NX),NY,NX)-CORGM))
        COHS=2.0_r8+10._r8*(CCLAY(NU(NY,NX),NY,NX)+CORGM) &
          +5.0_r8*(1.0_r8-EXP(-2.0E-06_r8*totRootLenDens_vr(NU(NY,NX),NY,NX)))
        SoilDetachability4Erosion2(NY,NX)=0.79_r8*EXP(-0.85_r8*AMAX1(1.0_r8,COHS))

        ParticleDensitySurfLay(NY,NX)=1.30_r8*CORGM+2.66_r8*(1.0_r8-CORGM)
        VISCWL=VISCW*EXP(0.533_r8-0.0267_r8*TCS(0,NY,NX))
        VLS(NY,NX)=3.6E+03_r8*9.8_r8*(ParticleDensitySurfLay(NY,NX)-1.0_r8) &
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
!     ActiveLayDepth=active layer depth
!     CDPTH,DLYR=depth to bottom,thickness of soil layer
!
    IF(FoundActiveLayer)exit
    !current layer
    VLiceTot=VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)
    VLPoreTot=VLMicP_vr(L,NY,NX)+VLMacP_vr(L,NY,NX)

    IF(VLPoreTot.GT.ZEROS2(NY,NX).AND.VLiceTot.GT.0.01*VLPoreTot)THEN
      !significant ice and pore
      D5700: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
        VLiceTotL=VLiceMicP_vr(LL,NY,NX)+VLiceMacP_vr(LL,NY,NX)
        VOLWTL=VLWatMicP_vr(LL,NY,NX)+VLWatMacP_vr(LL,NY,NX)
        VLPoreTotL=VLMicP_vr(LL,NY,NX)+VLMacP_vr(LL,NY,NX)
        !defined as significant increase in ice content
        goto5701=(VLPoreTotL.GT.ZEROS2(NY,NX).AND.VLiceTotL.LT.ZERO2*VLPoreTotL)
        if(goto5701)exit
      ENDDO D5700

      if(.not. goto5701)then
        IF(VLPoreTot.GT.ZEROS2(NY,NX))THEN
          ActiveLayDepth(NY,NX)=CumDepz2LayerBot_vr(L,NY,NX)-DLYR(3,L,NY,NX)*AMIN1(1.0_r8,VLiceTot/VLPoreTot)
        ELSE
          ActiveLayDepth(NY,NX)=CumDepz2LayerBot_vr(L,NY,NX)-DLYR(3,L,NY,NX)
        ENDIF
        FoundActiveLayer=.true.
      else
        !not active layer
        ActiveLayDepth(NY,NX)=9999.0_r8
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
  integer :: K,L,NTG
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
  trc_gascl_vr(idg_CO2,0,NY,NX) = CO2E_col(NY,NX)*5.36E-04_r8*TREF/TKS_vr(0,NY,NX)
  trc_gascl_vr(idg_CH4,0,NY,NX) = CH4E_col(NY,NX)*5.36E-04_r8*TREF/TKS_vr(0,NY,NX)
  trc_gascl_vr(idg_O2,0,NY,NX)  = OXYE(NY,NX)*1.43E-03_r8*TREF/TKS_vr(0,NY,NX)
  trc_gascl_vr(idg_N2,0,NY,NX)  = Z2GE(NY,NX)*1.25E-03_r8*TREF/TKS_vr(0,NY,NX)
  trc_gascl_vr(idg_N2O,0,NY,NX) = Z2OE(NY,NX)*1.25E-03_r8*TREF/TKS_vr(0,NY,NX)
  trc_gascl_vr(idg_NH3,0,NY,NX) = ZNH3E_col(NY,NX)*6.25E-04_r8*TREF/TKS_vr(0,NY,NX)
  trc_gascl_vr(idg_H2,0,NY,NX)  = H2GE(NY,NX)*8.92E-05*TREF/TKS_vr(0,NY,NX)

! initialize all band nutrients to zero
  trc_solcl_vr(ids_nutb_beg:ids_nutb_end,0,NY,NX)=0._r8
  IF(VLWatMicP_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
! exclude NH3B,
    DO NTG=idg_beg,idg_end-1
      trc_solcl_vr(NTG,0,NY,NX)=AZMAX1(trc_solml_vr(NTG,0,NY,NX)/VLWatMicP_vr(0,NY,NX))
    ENDDO
  ELSE
    trc_solcl_vr(idg_beg:idg_end-1,0,NY,NX)=0._r8
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
        trc_solcl_vr(ids_NH4,L,NY,NX)=AZMAX1(trc_solml_vr(ids_NH4,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)))
        trc_solcl_vr(idg_NH3,L,NY,NX)=AZMAX1(trc_solml_vr(idg_NH3,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(idg_NH3,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NH4,L,NY,NX)=0._r8
        trc_solcl_vr(idg_NH3,L,NY,NX)=0._r8
      ENDIF
      IF(trcs_VLN_vr(ids_NO3,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_NO3,L,NY,NX)=AZMAX1(trc_solml_vr(ids_NO3,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)))
        trc_solcl_vr(ids_NO2,L,NY,NX)=AZMAX1(trc_solml_vr(ids_NO2,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO2,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NO3,L,NY,NX)=0._r8
        trc_solcl_vr(ids_NO2,L,NY,NX)=0._r8
      ENDIF

      IF(trcs_VLN_vr(ids_H1PO4,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_H1PO4,L,NY,NX)=AZMAX1(trc_solml_vr(ids_H1PO4,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)))
        trc_solcl_vr(ids_H2PO4,L,NY,NX)=AZMAX1(trc_solml_vr(ids_H2PO4,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H2PO4,L,NY,NX)))

        CPO4S(L,NY,NX)=trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_solml_vr(ids_H2PO4,L,NY,NX)

        if(salt_model)then
          CPO4S(L,NY,NX)=CPO4S(L,NY,NX)+(trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX))*patomw
        endif  

        CPO4S(L,NY,NX)=AZMAX1(CPO4S(L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)))

      ELSE
        trc_solcl_vr(ids_H1PO4,L,NY,NX)=0._r8
        trc_solcl_vr(ids_H2PO4,L,NY,NX)=0._r8
        CPO4S(L,NY,NX)=0._r8
      ENDIF
!
!     C*B=solute concentration in band
!     CH1PB,CH2PB=HPO4,H2PO4 concentration in band
!     Z*B=P ion pair amounts in band (see solute.f)
!     VLNHB,VLNOB,VLPOB=fraction of soil volume in NH4,NO3,PO4 band
!
      IF(trcs_VLN_vr(ids_NH4B,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_NH4B,L,NY,NX)=AZMAX1(trc_solml_vr(ids_NH4B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)))
        trc_solcl_vr(idg_NH3B,L,NY,NX)=AZMAX1(trc_solml_vr(idg_NH3B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(idg_NH3B,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NH4B,L,NY,NX)=0._r8
        trc_solcl_vr(idg_NH3B,L,NY,NX)=0._r8
      ENDIF

      IF(trcs_VLN_vr(ids_NO3B,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_NO3B,L,NY,NX)=AZMAX1(trc_solml_vr(ids_NO3B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)))
        trc_solcl_vr(ids_NO2B,L,NY,NX)=AZMAX1(trc_solml_vr(ids_NO2B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_NO2B,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_NO3B,L,NY,NX)=0._r8
        trc_solcl_vr(ids_NO2B,L,NY,NX)=0._r8
      ENDIF
      IF(trcs_VLN_vr(ids_H1PO4B,L,NY,NX).GT.ZERO)THEN
        trc_solcl_vr(ids_H1PO4B,L,NY,NX)=AZMAX1(trc_solml_vr(ids_H1PO4B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)))
        trc_solcl_vr(ids_H2PO4B,L,NY,NX)=AZMAX1(trc_solml_vr(ids_H2PO4B,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H2PO4B,L,NY,NX)))

        CPO4B(L,NY,NX)=trc_solml_vr(ids_H1PO4B,L,NY,NX)+trc_solml_vr(ids_H2PO4B,L,NY,NX)
        if(salt_model)then
          CPO4B(L,NY,NX)=CPO4B(L,NY,NX)+(trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX) &
            +trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX))*patomw 
        endif  
        CPO4B(L,NY,NX)=AZMAX1(CPO4B(L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)))
      ELSE
        trc_solcl_vr(ids_H1PO4B,L,NY,NX) = 0._r8
        trc_solcl_vr(ids_H2PO4B,L,NY,NX) = 0._r8
        CPO4B(L,NY,NX)                   = 0._r8
      ENDIF

    ELSE
      trc_solcl_vr(ids_nuts_beg:ids_nuts_end,L,NY,NX) = 0._r8
      CPO4S(L,NY,NX)                                  = 0._r8
      CPO4B(L,NY,NX)                                  = 0._r8
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

    SoluteDifusvty_vr(idg_CO2,L,NY,NX)   = CLSG*TFACL
    SoluteDifusvty_vr(idg_CH4,L,NY,NX)   = CQSG*TFACL
    SoluteDifusvty_vr(idg_O2,L,NY,NX)    = OLSG*TFACL
    SoluteDifusvty_vr(idg_N2,L,NY,NX)    = ZLSG*TFACL
    SoluteDifusvty_vr(idg_NH3,L,NY,NX)   = ZNSG*TFACL
    SoluteDifusvty_vr(idg_H2,L,NY,NX)    = HLSG*TFACL
    SoluteDifusvty_vr(idg_N2O,L,NY,NX)   = ZVSG*TFACL
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
      ZA2=trcSalt_solml_vr(idsalt_SO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CO3,L,NY,NX)+trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_solml_vr(ids_H1PO4B,L,NY,NX)

      ZC1=(trc_solml_vr(ids_NH4,L,NY,NX)+trc_solml_vr(ids_NH4B,L,NY,NX))/natomw+trcSalt_solml_vr(idsalt_Hp,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_Na,L,NY,NX)+trcSalt_solml_vr(idsalt_K,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_AlOH2,L,NY,NX)+trcSalt_solml_vr(idsalt_FeOH2,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_AlSO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeSO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_CaOH,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHCO3,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_MgOH2,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHCO3,L,NY,NX)&
        +trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)

      ZA1=(trc_solml_vr(ids_NO3,L,NY,NX)+trc_solml_vr(ids_NO3B,L,NY,NX))/natomw+trcSalt_solml_vr(idsalt_OH,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_HCO3,L,NY,NX)+trcSalt_solml_vr(idsalt_Cl,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_AlOH4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeOH4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_NaCO3,L,NY,NX)+trcSalt_solml_vr(idsalt_NaSO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_KSO4,L,NY,NX)+(trc_solml_vr(ids_H2PO4,L,NY,NX) &
        +trc_solml_vr(ids_H2PO4B,L,NY,NX))/patomw+trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX) &
        +trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)

      ZN=trc_solml_vr(idg_CO2,L,NY,NX)/catomw+trc_solml_vr(idg_CH4,L,NY,NX)/catomw+trc_solml_vr(idg_O2,L,NY,NX)/32.0 &
        +(trc_solml_vr(idg_N2,L,NY,NX)+trc_solml_vr(idg_N2O,L,NY,NX)+trc_solml_vr(idg_NH3,L,NY,NX)+trc_solml_vr(idg_NH3B,L,NY,NX))/natomw &
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
  real(r8),intent(out) :: DPTH0(JY,JX)
  real(r8) :: VxcessWatLitR,TVOLWI,ThetaWLitR
  real(r8) :: VWatLitrZ
  real(r8) :: VOLIRZ
  real(r8) :: XVOLW0
  real(r8) :: XVOLI0
  integer  :: NTG,NTN
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
      VWatLitrZ = VLWatMicP_vr(0,NY,NX)/TVOLWI*VWatLitRHoldCapcity_col(NY,NX)
      VOLIRZ    = VLiceMicP_vr(0,NY,NX)/TVOLWI*VWatLitRHoldCapcity_col(NY,NX)
      XVOLW0    = AZMAX1(VLWatMicP_vr(0,NY,NX)-VWatLitrZ)/AREA(3,NU(NY,NX),NY,NX)
      XVOLI0    = AZMAX1(VLiceMicP_vr(0,NY,NX)-VOLIRZ)/AREA(3,NU(NY,NX),NY,NX)
    ELSE
      XVOLW0=0._r8
      XVOLI0=0._r8
    ENDIF
    DPTH0(NY,NX)    = XVOLW0+XVOLI0
    DLYR(3,0,NY,NX) = VLSoilPoreMicP_vr(0,NY,NX)/AREA(3,0,NY,NX)

    IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX) .AND. VLWatMicP_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VLWatMicP_vr(0,NY,NX))/VLitR_col(NY,NX)
      IF(ThetaWLitR.LT.FieldCapacity_vr(0,NY,NX))THEN
        PSISoilMatricP_vr(0,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX)+((LOGFldCapacity_vr(0,NY,NX)-LOG(ThetaWLitR)) &
          /FCD(0,NY,NX)*LOGPSIMND(NY,NX))))
      ELSEIF(ThetaWLitR.LT.POROS_vr(0,NY,NX))THEN
        PSISoilMatricP_vr(0,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(0,NY,NX)-LOG(ThetaWLitR)) &
          /PSD(0,NY,NX))**SRP(0,NY,NX)*LOGPSIMXD(NY,NX)))
      ELSE
        PSISoilMatricP_vr(0,NY,NX)=PSISE_vr(0,NY,NX)
      ENDIF
      PSISoilOsmotic_vr(0,NY,NX)     = 0._r8
      PSIGrav_vr(0,NY,NX)            = mGravAccelerat*(ALT(NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)+0.5_r8*DLYR(3,0,NY,NX))
      TotalSoilH2OPSIMPa_vr(0,NY,NX) = AZMIN1(PSISoilMatricP_vr(0,NY,NX)+PSISoilOsmotic_vr(0,NY,NX)+PSIGrav_vr(0,NY,NX))
!
!     LITTER NH4,NH3,NO3,NO2,HPO4,H2PO4 CONCENTRATIONS
!
!     C*=litter solute concentrations
!
      DO NTN=ids_nut_beg,ids_nuts_end
        trc_solcl_vr(NTN,0,NY,NX)=AZMAX1(trc_solml_vr(NTN,0,NY,NX)/VLWatMicP_vr(0,NY,NX))
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
    DLYR(3,0,NY,NX)                                = 0._r8
    THETW_vr(0,NY,NX)                              = 0._r8
    THETI_vr(0,NY,NX)                              = 0._r8
    ThetaAir_vr(0,NY,NX)                           = 1.0_r8
    VWatLitRHoldCapcity_col(NY,NX)                 = 0._r8
    PSISoilMatricP_vr(0,NY,NX)                     = PSISoilMatricP_vr(NU(NY,NX),NY,NX)
    trc_solcl_vr(ids_nut_beg:ids_nuts_end,0,NY,NX) = 0._r8
    trc_solcl_vr(idg_beg:idg_end-1,0,NY,NX)        = 0._r8
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
      IF(J.EQ.INT(SolarNoonHour_col(NY,NX)))THEN
        if(lverb)write(iulog,*)'ApplyMineralFertilizer'
        call ApplyMineralFertilizer(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
!
!     SOIL LAYER NUMBER IN WHICH PLANT OR ANIMAL RESIDUES ARE APPLIED
!
       if(lverb)write(iulog,*) 'ApplyManure'
        call ApplyManure(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
!
!     FERTILIZER UREA, NITRIFICATION INHIBITORS
        if(lverb)write(iulog,*)'ApplyUreaNitrifierInhibitor'
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
! urea application
  IF(FERT(ifert_urea,I,NY,NX).GT.0._r8 .OR. FERT(ifert_urea_band,I,NY,NX).GT.0._r8)THEN  
    IF(IYTYP(0,I,NY,NX).EQ.0)THEN
      IUTYP(NY,NX)=0
    ELSEIF(IYTYP(0,I,NY,NX).EQ.1 .OR. IYTYP(0,I,NY,NX).EQ.3)THEN
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
        ZNHU0(L,NY,NX)=0._r8
        ZNHUI(L,NY,NX)=0._r8
      ENDIF
    ENDDO D9964
  ENDIF

  IF(IYTYP(0,I,NY,NX).EQ.3 .OR. IYTYP(0,I,NY,NX).EQ.4)THEN
    D9965: DO L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN
        ZNFN0(L,NY,NX)=1.0_r8
        ZNFNI(L,NY,NX)=1.0_r8
      ELSE
        ZNFN0(L,NY,NX)=0._r8
        ZNFNI(L,NY,NX)=0._r8
      ENDIF
    ENDDO D9965
  ENDIF
  end subroutine ApplyUreaNitrifierInhibitor
!------------------------------------------------------------------------------------------

  subroutine ApplyManure(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
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
  integer  :: L,K,M,N,NN,NGL,MID
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
  IF(OFC(1)+OFC(2).GT.0._r8)THEN
    DO  L=0,JZ
      FDPTHM=FDPTH(I,NY,NX)+CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
      IF(FDPTHM.LE.0._r8)THEN
        LFDPTH=0
        exit
      ELSEIF(CumDepz2LayerBot_vr(L,NY,NX).GE.FDPTHM)THEN
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
!     DISTRIBUTE RESIDUE APPLICATION AMONG COMPONENTS OF RESIDUE complex
!
!     OFC,OFN,OFP=litter C,N,P application from fertilizer file
!

    D2965: DO K=1,2
      OSCI=OFC(K)*AREA(3,LFDPTH,NY,NX)
      OSNI=OFN(K)*AREA(3,LFDPTH,NY,NX)
      OSPI=OFP(K)*AREA(3,LFDPTH,NY,NX)
      IF(VLSoilMicPMass_vr(LFDPTH,NY,NX).GT.ZEROS(NY,NX))THEN
        CORGCX=OSCI/VLSoilMicPMass_vr(LFDPTH,NY,NX)
      ELSE
        CORGCX=orgcden
      ENDIF
      OSCX=0._r8
      OSNX=0._r8
      OSPX=0._r8
!
!     BIOMASSES OF MICROBIAL POPULATIONS IN RESIDUE
!
!     OMC,OMN,OMP=microbial biomass in litter application
!     OMCI=microbial biomass content in litter
!     OMCF,OMCA=hetero,autotrophic biomass composition in litter
!
      D2960: DO N=1,NumMicbFunGrupsPerCmplx
        tglds=JGnfo(N)-JGnfo(N)+1
        D2961: DO M=1,nlbiomcp
          OMC1=AZMAX1(AMIN1(OSCI*micpar%OMCI(M,K)*micpar%OMCF(N),OSCI-OSCX))
          OMN1=AZMAX1(AMIN1(OMC1*micpar%rNCOMCa(M,N,K),OSNI-OSNX))
          OMP1=AZMAX1(AMIN1(OMC1*micpar%rPCOMCa(M,N,K),OSPI-OSPX))
          DO NGL=JGnio(N),JGnfo(N)
            MID=micpar%get_micb_id(M,NGL)
            OMC1g=OMC1/tglds
            OMN1g=OMN1/tglds
            OMP1g=OMP1/tglds
            mBiomeHeter_vr(ielmc,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmc,MID,K,LFDPTH,NY,NX)+OMC1g
            mBiomeHeter_vr(ielmn,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmn,MID,K,LFDPTH,NY,NX)+OMN1g
            mBiomeHeter_vr(ielmp,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmp,MID,K,LFDPTH,NY,NX)+OMP1g
          ENDDO
          OSCX=OSCX+OMC1
          OSNX=OSNX+OMN1
          OSPX=OSPX+OMP1
          D2962: DO NN=1,NumMicbFunGrupsPerCmplx
            tglds=JGnfA(N)-JGniA(N)+1
            DO NGL=JGniA(NN),JGnfA(NN)
              MID=micpar%get_micb_id(M,NGL)
              OMC1g=OMC1/tglds
              OMN1g=OMN1/tglds
              OMP1g=OMP1/tglds
              mBiomeAutor_vr(ielmc,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmc,MID,LFDPTH,NY,NX)+OMC1g*micpar%OMCA(NN)
              mBiomeAutor_vr(ielmn,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmn,MID,LFDPTH,NY,NX)+OMN1g*micpar%OMCA(NN)
              mBiomeAutor_vr(ielmp,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmp,MID,LFDPTH,NY,NX)+OMP1g*micpar%OMCA(NN)
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

      DOM_vr(idom_doc,K,LFDPTH,NY,NX)=DOM_vr(idom_doc,K,LFDPTH,NY,NX)+OQC1
      DOM_vr(idom_don,K,LFDPTH,NY,NX)=DOM_vr(idom_don,K,LFDPTH,NY,NX)+OQN1
      DOM_vr(idom_dop,K,LFDPTH,NY,NX)=DOM_vr(idom_dop,K,LFDPTH,NY,NX)+OQP1
!
!     REMAINDER DISTRIBUTED TO RESIDUE FRACTIONS
!
!     OSC,OSN,OSP,OSA=SOC,SON,SOP,colonized SOC in litter
!     VOLT=litter volume
!     AmendCFlx_CumYr_col,FertNFlx_CumYr_col,FerPFlx_CumYr_col=accumulated litter C,N,P application
!     Eco_NBP_CumYr_col=accumulated net biome productivity
!
      OSCX=OSCX+OQC1
      OSNX=OSNX+OQN1
      OSPX=OSPX+OQP1
      CNOFT=0._r8
      CPOFT=0._r8
      IF(OSCI-OSCX.GT.ZEROS(NY,NX))THEN
        RNT=0._r8
        RPT=0._r8
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
          CNOF(M)=0._r8
          CPOF(M)=0._r8
        ENDDO D975
      ENDIF
      D2970: DO M=1,jsken
        OSC1=CFOSC(M,K,LFDPTH,NY,NX)*(OSCI-OSCX)
        IF(CNOFT.GT.ZERO)THEN
          OSN1=CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)/CNOFT*(OSNI-OSNX)
        ELSE
          OSN1=0._r8
        ENDIF
        IF(CPOFT.GT.ZERO)THEN
          OSP1=CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)/CPOFT*(OSPI-OSPX)
        ELSE
          OSP1=0._r8
        ENDIF
        SolidOM_vr(ielmc,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmc,M,K,LFDPTH,NY,NX)+OSC1
        SolidOM_vr(ielmn,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmn,M,K,LFDPTH,NY,NX)+OSN1
        SolidOM_vr(ielmp,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmp,M,K,LFDPTH,NY,NX)+OSP1
        SolidOMAct_vr(M,K,LFDPTH,NY,NX)=SolidOMAct_vr(M,K,LFDPTH,NY,NX)+OSC1*micpar%OMCI(1,K)
        IF(LFDPTH.EQ.0)THEN
          VGeomLayer_vr(LFDPTH,NY,NX)=VGeomLayer_vr(LFDPTH,NY,NX)+OSC1*ppmc/BulkDensLitR(micpar%k_fine_litr)
        ENDIF
      ENDDO D2970
      TORGF=TORGF+OSCI
      TORGN=TORGN+OSNI
      TORGP=TORGP+OSPI
      AmendCFlx_CumYr_col(NY,NX)=AmendCFlx_CumYr_col(NY,NX)+OSCI
      FertNFlx_CumYr_col(NY,NX)=FertNFlx_CumYr_col(NY,NX)+OSNI
      FerPFlx_CumYr_col(NY,NX)=FerPFlx_CumYr_col(NY,NX)+OSPI
      IF(IYTYP(2,I,NY,NX).LT.3)THEN
        Eco_NBP_CumYr_col(NY,NX)=Eco_NBP_CumYr_col(NY,NX)+OSCI
      ENDIF
    ENDDO D2965
  ENDIF
  end associate
  end subroutine ApplyManure
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
  IF(Z4A+Z3A+ZUA+ZOA+Z4B+Z3B+ZUB+ZOB+PMA+PMB+PHA+CAC+CAS.GT.0._r8)THEN
    FDPTHF=FDPTH(I,NY,NX)+CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
    IF(FDPTHF.LE.0._r8.AND.isclose(Z4B+Z3B+ZUB+ZOB+PMB,0._r8))THEN
      LFDPTH=0
      CVRDF=1.0_r8-EXP(-0.8E-02_r8*(SoilOrgM_vr(ielmc,0,NY,NX)/AREA(3,0,NY,NX)))
    ELSE
      D65: DO L=NUI(NY,NX),JZ
        IF(CumDepz2LayerBot_vr(L,NY,NX).GE.FDPTHF)THEN
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
    IF((Z4B+Z3B+ZUB.GT.0._r8).OR.((trc_solml_vr(ids_NH4B,LFDPTH,NY,NX).GT.0._r8 &
      .OR.trc_solml_vr(idg_NH3B,LFDPTH,NY,NX).GT.0._r8).AND.IFNHB(NY,NX).EQ.0))THEN
      IFNHB(NY,NX)=1
      ROWN(NY,NX)=ROWI(I,NY,NX)
      D50: DO L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessNH4_vr(L,NY,NX)=DLYR(3,L,NY,NX)
          BandWidthNH4_vr(L,NY,NX)=0._r8
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessNH4_vr(L,NY,NX)=AMAX1(0.025_r8,FDPTHF-CumDepz2LayerBot_vr(L-1,NY,NX))
          BandWidthNH4_vr(L,NY,NX)=AMIN1(0.025_r8,ROWN(NY,NX))
        ELSE
          BandThicknessNH4_vr(L,NY,NX)=0._r8
          BandWidthNH4_vr(L,NY,NX)=0._r8
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_NH4B,L,NY,NX)=AMIN1(0.999_r8,BandWidthNH4_vr(L,NY,NX)/ROWN(NY,NX) &
            *BandThicknessNH4_vr(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_NH4B,L,NY,NX)=0._r8
        ENDIF
        trcs_VLN_vr(ids_NH4,L,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_VLN_vr(idg_NH3B,L,NY,NX)=trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_VLN_vr(idg_NH3,L,NY,NX)=trcs_VLN_vr(ids_NH4,L,NY,NX)
        ZNH4T=trc_solml_vr(ids_NH4,L,NY,NX)+trc_solml_vr(ids_NH4B,L,NY,NX)
        ZNH3T=trc_solml_vr(idg_NH3,L,NY,NX)+trc_solml_vr(idg_NH3B,L,NY,NX)
        XN4T=trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX)
        trc_solml_vr(ids_NH4,L,NY,NX)=ZNH4T*trcs_VLN_vr(ids_NH4,L,NY,NX)
        trc_solml_vr(idg_NH3,L,NY,NX)=ZNH3T*trcs_VLN_vr(idg_NH3,L,NY,NX)
        trc_solml_vr(ids_NH4B,L,NY,NX)=ZNH4T*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trc_solml_vr(idg_NH3B,L,NY,NX)=ZNH3T*trcs_VLN_vr(idg_NH3B,L,NY,NX)
        trcx_solml_vr(idx_NH4,L,NY,NX)=XN4T*trcs_VLN_vr(ids_NH4,L,NY,NX)
        trcx_solml_vr(idx_NH4B,L,NY,NX)=XN4T*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      ENDDO D50
      BandDepthNH4_col(NY,NX)=BandThicknessNH4_vr(LFDPTH,NY,NX)+CumDepz2LayerBot_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF NO3 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWO=width of NO3 band row
!     DPNOB,WDNOB=depth,width of NO3 band
!     VLNOB,VLNO3=soil volume in NO3 band,non-band
!
    IF((Z4B+Z3B+ZUB+ZOB.GT.0._r8).OR.((trc_solml_vr(ids_NO3B,LFDPTH,NY,NX).GT.0._r8 &
      .OR.trc_solml_vr(ids_NO2B,LFDPTH,NY,NX).GT.0._r8).AND.IFNOB(NY,NX).EQ.0))THEN
      IFNOB(NY,NX)=1
      ROWO(NY,NX)=ROWI(I,NY,NX)
      D45: DO L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessNO3_vr(L,NY,NX)=DLYR(3,L,NY,NX)
          BandWidthNO3_vr(L,NY,NX)=0._r8
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessNO3_vr(L,NY,NX)=AMAX1(0.01_r8,FDPTHF-CumDepz2LayerBot_vr(L-1,NY,NX))
          BandWidthNO3_vr(L,NY,NX)=AMIN1(0.01_r8,ROWO(NY,NX))
        ELSE
          BandThicknessNO3_vr(L,NY,NX)=0._r8
          BandWidthNO3_vr(L,NY,NX)=0._r8
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_NO3B,L,NY,NX)=AMIN1(0.999_r8,BandWidthNO3_vr(L,NY,NX)/ROWO(NY,NX) &
            *BandThicknessNO3_vr(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_NO3B,L,NY,NX)=0._r8
        ENDIF

        trcs_VLN_vr(ids_NO3,L,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_VLN_vr(ids_NO2B,L,NY,NX)=trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_VLN_vr(ids_NO2,L,NY,NX)=trcs_VLN_vr(ids_NO3,L,NY,NX)
        ZNO3T=trc_solml_vr(ids_NO3,L,NY,NX)+trc_solml_vr(ids_NO3B,L,NY,NX)
        ZNO2T=trc_solml_vr(ids_NO2,L,NY,NX)+trc_solml_vr(ids_NO2B,L,NY,NX)

        trc_solml_vr(ids_NO3,L,NY,NX)=ZNO3T*trcs_VLN_vr(ids_NO3,L,NY,NX)
        trc_solml_vr(ids_NO2,L,NY,NX)=ZNO2T*trcs_VLN_vr(ids_NO2,L,NY,NX)
        trc_solml_vr(ids_NO3B,L,NY,NX)=ZNO3T*trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trc_solml_vr(ids_NO2B,L,NY,NX)=ZNO2T*trcs_VLN_vr(ids_NO2B,L,NY,NX)
      ENDDO D45
      BandDepthNO3_col(NY,NX)=BandThicknessNO3_vr(LFDPTH,NY,NX)+CumDepz2LayerBot_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF PO4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWP=width of H2PO4 band row
!     DPPOB,WDPOB=depth,width of H2PO4 band
!     VLPOB,VLPO4=soil volume in H2PO4 band,non-band
!
    IF((PMB.GT.0.0).OR.(trc_solml_vr(ids_H2PO4B,LFDPTH,NY,NX).GT.0._r8.AND.IFPOB(NY,NX).EQ.0))THEN
      IFPOB(NY,NX)=1
      ROWP(NY,NX)=ROWI(I,NY,NX)
      DO  L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessPO4_vr(L,NY,NX)=DLYR(3,L,NY,NX)
          BandWidthPO4_vr(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessPO4_vr(L,NY,NX)=AMAX1(0.01,FDPTHF-CumDepz2LayerBot_vr(L-1,NY,NX))
          BandWidthPO4_vr(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
        ELSE
          BandThicknessPO4_vr(L,NY,NX)=0._r8
          BandWidthPO4_vr(L,NY,NX)=0._r8
        ENDIF
        IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=AMIN1(0.999_r8,BandWidthPO4_vr(L,NY,NX)/ROWP(NY,NX) &
          *BandThicknessPO4_vr(L,NY,NX)/DLYR(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=0._r8
        ENDIF
        trcs_VLN_vr(ids_H1PO4,L,NY,NX)=1.0-trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_VLN_vr(ids_H2PO4B,L,NY,NX)=trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_VLN_vr(ids_H2PO4,L,NY,NX)=trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        H1PO4T=trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_solml_vr(ids_H1PO4B,L,NY,NX)
        H2PO4T=trc_solml_vr(ids_H2PO4,L,NY,NX)+trc_solml_vr(ids_H2PO4B,L,NY,NX)

        XOH0T=trcx_solml_vr(idx_OHe,L,NY,NX)+trcx_solml_vr(idx_OHeB,L,NY,NX)
        XOH1T=trcx_solml_vr(idx_OH,L,NY,NX)+trcx_solml_vr(idx_OHB,L,NY,NX)
        XOH2T=trcx_solml_vr(idx_OHp,L,NY,NX)+trcx_solml_vr(idx_OHpB,L,NY,NX)
        XH1PT=trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_HPO4B,L,NY,NX)
        XH2PT=trcx_solml_vr(idx_H2PO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX)
        PALPOT=trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)+trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)
        PFEPOT=trcp_saltpml_vr(idsp_FePO4,L,NY,NX)+trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)
        PCAPDT=trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)+trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)
        PCAPHT=trcp_saltpml_vr(idsp_HA,L,NY,NX)+trcp_saltpml_vr(idsp_HAB,L,NY,NX)
        PCAPMT=trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)+trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX)

        trc_solml_vr(ids_H1PO4,L,NY,NX)=H1PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trc_solml_vr(ids_H2PO4,L,NY,NX)=H2PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trc_solml_vr(ids_H1PO4B,L,NY,NX)=H1PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trc_solml_vr(ids_H2PO4B,L,NY,NX)=H2PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

        trcx_solml_vr(idx_OHe,L,NY,NX)=XOH0T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OH,L,NY,NX)=XOH1T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OHp,L,NY,NX)=XOH2T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_HPO4,L,NY,NX)=XH1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_H2PO4,L,NY,NX)=XH2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OHeB,L,NY,NX)=XOH0T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_OHB,L,NY,NX)=XOH1T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_OHpB,L,NY,NX)=XOH2T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_HPO4B,L,NY,NX)=XH1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_H2PO4B,L,NY,NX)=XH2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)=PALPOT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_FePO4,L,NY,NX)=PFEPOT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)=PCAPDT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_HA,L,NY,NX)=PCAPHT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)=PCAPMT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)=PALPOT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)=PFEPOT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)=PCAPDT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_HAB,L,NY,NX)=PCAPHT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX)=PCAPMT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

        if(salt_model)then
          H0PO4T=trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)
          H3PO4T=trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)
          ZFE1PT=trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)
          ZFE2PT=trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)
          ZCA0PT=trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)
          ZCA1PT=trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)
          ZCA2PT=trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)
          ZMG1PT=trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)
          
          trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)=H0PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)=H3PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)=ZFE1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)=ZFE2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)=ZCA0PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)=ZCA1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)=ZCA2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)=ZMG1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)=H0PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

          trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)=H3PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)=ZFE1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)=ZFE2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)=ZCA0PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)=ZCA1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)=ZCA2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)=ZMG1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        endif

      ENDDO
      BandDepthPO4_col(NY,NX)=BandThicknessPO4_vr(LFDPTH,NY,NX)+CumDepz2LayerBot_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     UPDATE STATE VARIABLES FOR BROADCAST AND BANDED FERTILIZER
!     NH4, NH3, UREA, NO3, PO4, LIME AND GYPSUM IN SOIL
!     AND CONVERT FROM G TO MOLE
!
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=bdcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH4P2O8,CaHPO4,apatite in non-band
!     PCAPMB,PCAPDB,PCAPHB=concn of precip CaH4P2O8,CaHPO4,apatite in band
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
    CACX=CAC*AREA(3,LFDPTH,NY,NX)/40._r8
    CASX=CAS*AREA(3,LFDPTH,NY,NX)/40._r8

    FertN_soil_vr(ifert_nh4,LFDPTH,NY,NX)=FertN_soil_vr(ifert_nh4,LFDPTH,NY,NX)+Z4AX*CVRDF
    FertN_soil_vr(ifert_urea,LFDPTH,NY,NX)=FertN_soil_vr(ifert_urea,LFDPTH,NY,NX)+ZUAX*CVRDF
    FertN_soil_vr(ifert_no3,LFDPTH,NY,NX)=FertN_soil_vr(ifert_no3,LFDPTH,NY,NX)+ZOAX*CVRDF

    FertN_Band_vr(ifert_nh4_band,LFDPTH,NY,NX)=FertN_Band_vr(ifert_nh4_band,LFDPTH,NY,NX)+Z4BX*CVRDF
    FertN_Band_vr(ifert_urea_band,LFDPTH,NY,NX)=FertN_Band_vr(ifert_urea_band,LFDPTH,NY,NX)+ZUBX*CVRDF
    FertN_Band_vr(ifert_no3_band,LFDPTH,NY,NX)=FertN_Band_vr(ifert_no3_band,LFDPTH,NY,NX)+ZOBX*CVRDF

    trcp_saltpml_vr(idsp_CaH4P2O8,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8,LFDPTH,NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4,LFDPTH,NY,NX)*CVRDF
    trcp_saltpml_vr(idsp_CaH4P2O8B,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8B,LFDPTH,NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4B,LFDPTH,NY,NX)*CVRDF+PMBX*CVRDF
    trcp_saltpml_vr(idsp_HA,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_HA,LFDPTH,NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4,LFDPTH,NY,NX)*CVRDF
    trcp_saltpml_vr(idsp_HAB,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_HAB,LFDPTH,NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4B,LFDPTH,NY,NX)*CVRDF
    IF(LFDPTH.EQ.0)THEN
      FertN_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)+Z4AX*BAREF
      FertN_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)+Z3AX
      FertN_soil_vr(ifert_urea,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_urea,NU(NY,NX),NY,NX)+ZUAX*BAREF
      FertN_soil_vr(ifert_no3,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_no3,NU(NY,NX),NY,NX)+ZOAX*BAREF

      FertN_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)+Z4BX*BAREF
      FertN_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)+Z3BX
      FertN_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX)+ZUBX*BAREF
      FertN_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)+ZOBX*BAREF

      trcp_saltpml_vr(idsp_CaH4P2O8,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8,NU(NY,NX),NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)*BAREF
      trcp_saltpml_vr(idsp_CaH4P2O8B,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8B,NU(NY,NX),NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)*BAREF+PMBX*BAREF
      trcp_saltpml_vr(idsp_HA,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_HA,NU(NY,NX),NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)*BAREF
      trcp_saltpml_vr(idsp_HAB,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_HAB,NU(NY,NX),NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)*BAREF
    ELSE
      FertN_soil_vr(ifert_nh3,LFDPTH,NY,NX)=FertN_soil_vr(ifert_nh3,LFDPTH,NY,NX)+Z3AX*CVRDF
      FertN_Band_vr(ifert_nh3_band,LFDPTH,NY,NX)=FertN_Band_vr(ifert_nh3_band,LFDPTH,NY,NX)+Z3BX*CVRDF
    ENDIF
    trcp_saltpml_vr(idsp_CaCO3,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaCO3,NU(NY,NX),NY,NX)+CACX
    trcp_saltpml_vr(idsp_CaSO4,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaSO4,NU(NY,NX),NY,NX)+CASX
    TZIN=TZIN+natomw*(Z4AX+Z3AX+ZUAX+ZOAX+Z4BX+Z3BX+ZUBX+ZOBX)
    TPIN=TPIN+62.0_r8*(PMAX+PMBX)+93.0_r8*PHAX
    TIONIN=TIONIN+2.0_r8*(CACX+CASX)
    FertNFlx_CumYr_col(NY,NX)=FertNFlx_CumYr_col(NY,NX)+natomw*(Z4AX+Z4BX+Z3AX+Z3BX+ZUAX+ZUBX+ZOAX+ZOBX)
    FerPFlx_CumYr_col(NY,NX)=FerPFlx_CumYr_col(NY,NX)+62.0_r8*(PMAX+PMBX)+93.0_r8*PHAX
  ENDIF
  end subroutine ApplyMineralFertilizer
!------------------------------------------------------------------------------------------

  subroutine GetChemicalConcsInSoil(NY,NX,THETPZ_vr)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(out) :: THETPZ_vr(JZ,JY,JX)
  integer :: L,NTG
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
    THETPZ_vr(L,NY,NX)=AZMAX1(POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX)-THETI_vr(L,NY,NX))
!
!     GAS CONCENTRATIONS
!
!     C*G=soil gas gaseous concentration
!     C*S=soil gas aqueous concentration
!
    IF(ThetaAir_vr(L,NY,NX).GT.THETX)THEN
      DO NTG=idg_beg,idg_end-1
        trc_gascl_vr(NTG,L,NY,NX)=AZMAX1(trc_gasml_vr(NTG,L,NY,NX)/VLsoiAirP_vr(L,NY,NX))
      ENDDO
    ELSE
      trc_gascl_vr(idg_beg:idg_end-1,L,NY,NX)=0._r8
    ENDIF

    IF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DO NTG=idg_beg,idg_end-1
        trc_solcl_vr(NTG,L,NY,NX)=AZMAX1(trc_solml_vr(NTG,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
      ENDDO
    ELSE
      trc_solcl_vr(idg_beg:idg_end-1,L,NY,NX)=0._r8
    ENDIF
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

  DO L=NUI(NY,NX),NLI(NY,NX)
    FWatExMacP2MicP(L,NY,NX)                        = 0._r8
    trcs_plant_uptake_vr(ids_beg:ids_end,L,NY,NX)   = 0._r8
    tRootCO2Emis_vr(L,NY,NX)                        = 0._r8
    trcg_root_vr(idg_beg:idg_end-1,L,NY,NX)         = 0._r8
    tRO2MicrbUptk_vr(L,NY,NX)                       = 0._r8
    trcg_air2root_flx_vr(idg_beg:idg_end-1,L,NY,NX) = 0._r8

    trcn_RChem_band_soil_vr(ids_NH4B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(ids_NO3B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(ids_NO2B,L,NY,NX)   = 0._r8
    trcn_RChem_band_soil_vr(ids_H1PO4B,L,NY,NX) = 0._r8
    trcn_RChem_band_soil_vr(ids_H2PO4B,L,NY,NX) = 0._r8
    TR_CO2_gchem_soil_vr(L,NY,NX)                 = 0._r8
    Txchem_CO2_vr(L,NY,NX)                      = 0._r8

    trcx_TRSoilChem_vr(idx_NH4B,L,NY,NX)=0._r8
    trcx_TRSoilChem_vr(idx_OHeB:idx_end,L,NY,NX)=0._r8

    TR_H_p_sorbed_soil(L,NY,NX)       = 0._r8
    TR_Al_sorbed_soil(L,NY,NX)        = 0._r8
    TR_Fe_sorbed_soil_vr(L,NY,NX)     = 0._r8
    TR_Ca_sorbed_soil(L,NY,NX)        = 0._r8
    TR_Mg_sorbed_soil(L,NY,NX)        = 0._r8
    TR_Na_sorbed_soil(L,NY,NX)        = 0._r8
    TR_K_sorbed_soil(L,NY,NX)         = 0._r8
    TR_HCO3_sorbed_soil(L,NY,NX)      = 0._r8
    TR_AlO2H2_sorbed_soil(L,NY,NX)    = 0._r8
    TR_FeO2H2_sorbed_soil_vr(L,NY,NX) = 0._r8

    trcp_RChem_soil(idsp_beg:idsp_psoi_beg-1,L,NY,NX)=0._r8

    trcp_RChem_soil(idsp_beg_band:idsp_end,L,NY,NX)=0._r8

    trcs_Mac2MicXfer_vr(ids_beg:ids_end,L,NY,NX)=0._r8

    DO NTSA=idsalt_beg,idsaltb_end
      trcSalt_TR(NTSA,L,NY,NX)=0._r8
      trcSalt_XFXS(NTSA,L,NY,NX)=0._r8
    ENDDO

    DO  K=1,jcplx
      DOM_PoreTranspFlx(idom_beg:idom_end,K,L,NY,NX)=0._r8
    ENDDO
    TLIceThawMicP(L,NY,NX)=0._r8
    TLIceThawMacP(L,NY,NX)=0._r8
    TLPhaseChangeHeat2Soi(L,NY,NX)=0._r8
    trcg_ebu_flx_vr(idg_beg:idg_end,L,NY,NX)=0._r8
    totRootLenDens_vr(L,NY,NX)=0._r8

  ENDDO
  trcg_ebu_flx_col(idg_beg:idg_NH3,NY,NX)=0._r8
  trcg_pltroot_flx_col(idg_beg:idg_NH3,NY,NX)=0._r8
  end subroutine ZeroHourlyArrays

!------------------------------------------------------------------------------------------

  subroutine CalGasSolubility(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer  :: L,NTG
  real(r8) :: FH2O

  L=0
  DO NTG=idg_beg,idg_end-1
    GasSolbility_vr(NTG,L,NY,NX)=gas_solubility(NTG,TCS(L,NY,NX))
  ENDDO

  DO  L=1,NL(NY,NX)+1
! S*L=solubility of gas in water
! TCS=soil temperature (oC)
! 5.56E+04_r8 := mole H2O / m3
    FH2O=5.56E+04_r8/(5.56E+04_r8+SolutesIonConc_vr(L,NY,NX))
    DO NTG=idg_beg,idg_end-1
      GasSolbility_vr(NTG,L,NY,NX)=gas_solubility(NTG,TCS(L,NY,NX))*EXP(-ACTCG(NTG)*SolutesIonStrenth_vr(L,NY,NX))*FH2O
    ENDDO
  ENDDO
  end subroutine CalGasSolubility
!------------------------------------------------------------------------------------------
  subroutine UpdateLiterPropertz(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS  
  real(r8) :: FVLitR
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

      IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX))THEN
        FVLitR=VWatLitRHoldCapcity_col(NY,NX)/VLitR_col(NY,NX)
      ELSE
        FVLitR=THETRX(micpar%k_fine_litr)/BulkDensLitR(micpar%k_fine_litr)
      ENDIF
      POROS0(NY,NX)              = FVLitR
      FieldCapacity_vr(0,NY,NX)  = 0.500_r8*FVLitR
      WiltPoint_vr(0,NY,NX)      = 0.125_r8*FVLitR
      LOGPOROS_vr(0,NY,NX)       = LOG(POROS0(NY,NX))
      LOGFldCapacity_vr(0,NY,NX) = LOG(FieldCapacity_vr(0,NY,NX))
      LOGWiltPoint_vr(0,NY,NX)   = LOG(WiltPoint_vr(0,NY,NX))
      PSD(0,NY,NX)               = LOGPOROS_vr(0,NY,NX)-LOGFldCapacity_vr(0,NY,NX)
      FCD(0,NY,NX)               = LOGFldCapacity_vr(0,NY,NX)-LOGWiltPoint_vr(0,NY,NX)
      SRP(0,NY,NX)               = 1.00_r8
    enddo
  enddo    
  end subroutine UpdateLiterPropertz

!------------------------------------------------------------------------------------------
  subroutine SummarizeTracers(NHW,NHE,NVN,NVS)
  !
  !Description
  !sum up mass of tracers
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS  
  integer :: NY,NX
  integer :: idg,L

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      trcVolatileMass_col(idg_beg:idg_end,NY,NX)=0._r8
      DO L=NUI(NY,NX),NLI(NY,NX)
        DO idg=idg_beg,idg_end
          trcVolatileMass_col(idg,NY,NX)=trcVolatileMass_col(idg,NY,NX)+trc_gasml_vr(idg,L,NY,NX) + &
            trc_solml_vr(idg,L,NY,NX) + trc_soHml_vr(idg,L,NY,NX)
        ENDDO
      ENDDO
      DO idg=idg_beg,idg_NH3
        trcVolatileMass_col(idg,NY,NX)=trcVolatileMass_col(idg,NY,NX)+trc_solml_vr(idg,0,NY,NX)
        DO L=1,JS
          trcVolatileMass_col(idg,NY,NX)=trcVolatileMass_col(idg,NY,NX)+trcg_solsml_snvr(idg,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine SummarizeTracers
end module Hour1Mod
