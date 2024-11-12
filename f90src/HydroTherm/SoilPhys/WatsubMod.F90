module WatsubMod
!!
! Description:
! Do water and enerby balance calculation.
! The module diagnoses the mass and energy fluxes associated
! with soil/snow water (vapor, liquid and ice) and energy, and update
! them in redistmod.F90

  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: GravAcceleration=>DAT_CONST_G
  use abortutils,     only: endrun,           print_info
  use minimathmod,    only: isclose,          isclose,  safe_adb, vapsat, AZMAX1, AZMIN1, AZMAX1t
  use ElmIDMod,       only: iEastWestDirection,         iNorthSouthDirection, iVerticalDirection
  use SurfPhysData,   only: InitSurfPhysData, DestructSurfPhysData
  use SnowBalanceMod, only : DebugSnowPrint
  use EcoSIMCtrlMod,  only: lverb,snowRedist_model  
  use RootDataType
  use EcosimConst
  use MiniFuncMod
  use EcoSIMSolverPar
  use SOMDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilWaterDataType
  use SoilHeatDatatype
  use EcoSIMCtrlDataType
  use LandSurfDataType
  use ClimForcDataType
  use FertilizerDataType
  use SnowDataType
  use PlantTraitDataType
  use SurfLitterDataType
  use SurfSoilDataType
  use SurfSoilDataType
  use CanopyDataType
  use ChemTranspDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use IrrigationDataType
  use GridDataType
  use WatsubDataMod
  use PhysPars    
  use SnowPhysMod
  use HydroThermData
  use SurfPhysMod
  use SoilPhysParaMod
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8), parameter :: mGravAccelerat=1.e-3_r8*GravAcceleration  !gravitational constant devided by 1000.
  integer :: curday,curhour

  public :: watsub,InitWatsub
  public :: DestroyWatsub
  contains
!------------------------------------------------------------------------------------------
  subroutine InitWatsub

  implicit none

  call InitSurfPhysData  

  call InitWatSubData

  end subroutine InitWatsub
!------------------------------------------------------------------------------------------
  subroutine DestroyWatsub
  implicit none
  call DestructSurfPhysData
  end subroutine DestroyWatsub
!------------------------------------------------------------------------------------------

  SUBROUTINE watsub(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CACULATES ENERGY BALANCES OF SNOW, RESIDUE
!     AND SOIL SURFACES, FREEZING, THAWING, AND HEAT AND WATER
!     TRANSFER THROUGH SOIL PROFILES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: K0,K1
  integer :: KL,L,L2,LL,M,MM,M1,M2,M3,M4,M5,M6,NX,NY
  integer :: N,N1,N2,N3,N4,N5,N6,NN,N4B,N5B,NUX
  real(r8):: ResistanceLitRLay(JY,JX),HeatFluxAir2Soi(JY,JX)
  REAL(R8):: KSatRedusByRainKinetEnergyS(JY,JX)

  REAL(R8) :: TopLayWatVol(JY,JX)
  real(r8) :: Qinfl2MicP(JY,JX)
  real(r8) :: Hinfl2Soil(JY,JX)
! begin_execution
!

  call LocalCopySoilVars(I,J,NHW,NHE,NVN,NVS)

  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  call InitSoilHydrauics(NHW,NHE,NVN,NVS)

! DYNAMIC LOOP FOR FLUX CALCULATIONS
  ! iterate for NPH times
!  TLPhaseChangeHeat2Soi1s(:,:,:)=0._r8
  D3320: DO M=1,NPH

    call FWDCopyTopLayerWatVolMit(NHW,NHE,NVN,NVS,TopLayWatVol)

!   run surface energy balance model, uses ResistanceLitRLay, update top layer soil moisture
    call RunSurfacePhysModel(I,J,M,NHE,NHW,NVS,NVN,ResistanceLitRLay,KSatRedusByRainKinetEnergyS,&
      TopLayWatVol,HeatFluxAir2Soi,Qinfl2MicP,Hinfl2Soil)

!   prepare update for the other soil layers
    call CopySoilWatVolMit(I,J,M,NHW,NHE,NVN,NVS,TopLayWatVol)
        
    call Subsurface3DFlowMit(I,J,M,NHW,NHE,NVN,NVS,KSatRedusByRainKinetEnergyS,HeatFluxAir2Soi)

    call LateralWatHeatExchMit(I,J,M,NHW,NHE,NVN,NVS,KSatRedusByRainKinetEnergyS)

!   update states and fluxes
    DO NX=NHW,NHE
      DO  NY=NVN,NVS
        HeatFlx2Grnd_col(NY,NX) = HeatFlx2Grnd_col(NY,NX)+Hinfl2Soil(NY,NX)
        Qinflx2Soil_col(NY,NX)  = Qinflx2Soil_col(NY,NX)+Qinfl2MicP(NY,NX)        
      ENDDO
    ENDDO  

    if(snowRedist_model)call AccumulateSnowRedisFlux(I,J,M,NHW,NHE,NVN,NVS)
    

    call UpdateSoilMoistTemp(I,J,M,NHW,NHE,NVN,NVS)

    IF(M.NE.NPH)THEN
      call UpdateSurfaceAtM(I,J,M,NHW,NHE,NVN,NVS)

      call UpdateStateFluxAtM(M,NHW,NHE,NVN,NVS)
    ELSE
      call UpdateFluxAtExit(I,J,NHW,NHE,NVN,NVS)
    ENDIF

  ENDDO D3320

  END subroutine watsub

!------------------------------------------------------------------------------------------  

  subroutine FWDCopyTopLayerWatVolMit(NHW,NHE,NVN,NVS,TopLayWatVol)
  !
  !make a copy of toplayer water content
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(out) :: TopLayWatVol(JY,JX)
  integer :: NY,NX

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      TopLayWatVol(NY,NX)= VLWatMicP1_vr(NUM(NY,NX),NY,NX)
    ENDDO
  ENDDO
  end subroutine FWDCopyTopLayerWatVolMit

!------------------------------------------------------------------------------------------  

  subroutine CopySoilWatVolMit(I,J,M,NHW,NHE,NVN,NVS,TopLayWatVol)

  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: TopLayWatVol
  integer :: L,NY,NX
  
  !VOLW2 will be updated in soil surface model

  DO NX=NHW,NHE
    DO  NY=NVN,NVS  
      DO L=NUM(NY,NX)+1,NL(NY,NX)        
        VLWatMicP2_vr(L,NY,NX)=VLWatMicP1_vr(L,NY,NX)
      ENDDO  
      VLWatMicP2_vr(NUM(NY,NX),NY,NX)=TopLayWatVol(NY,NX)
    ENDDO
  ENDDO
  end subroutine CopySoilWatVolMit


!------------------------------------------------------------------------------------------  
  subroutine LocalCopySoilVars(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX

  integer :: L,LyrIrrig
  real(r8) :: VLTSoiPore,rsatMacP

  DX995: DO NX=NHW,NHE
    DX990: DO NY=NVN,NVS

    !make a local copy of the upper boundary index
      NUM(NY,NX)=NU(NY,NX)

    ! CDPTH=depth to bottom of soil layer
    ! WDPTH,LyrIrrig=depth,layer of subsurface irrigation

      !identify the layer where irrigation is applied
      D65: DO L=NUM(NY,NX),NL(NY,NX)
        IF(CumDepz2LayerBot_vr(L,NY,NX).GE.WDPTH(I,NY,NX))THEN
          LyrIrrig=L
          exit
        ENDIF
      ENDDO D65

      D30: DO L=NUM(NY,NX),NL(NY,NX)
    !
    !   ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
    !   FOR USE AT INTERNAL TIME STEP IN30    con SOIL LAYERS
    !
    !   PSISM1,PSISM=matric water potential
    !   VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of micropores
    !   VLWatMicPX1=VLWatMicP1 accounting for wetting front
    !   VOLAH*,VOLWH*,VOLIH*,VOLPH*=pore,water,ice,air macropores
    !   BKDS=bulk density
    !   CCLAY=clay concentration
    !   FVOLAH=parameter for clay effect on macropore volume
    !   VLSoilPoreMicP_vr,VOLT=soil,total volumes
    !   WP=wilting point
    !   THETW*,THETI*,THETP*=water,ice,air-filled porosity
    !   VHeatCapacity1_vr,VHCM=volumetric heat capacities of total volume, solid
    !   VLHeatCapacityA,VLHeatCapacityB=volumetric heat capacities of micropore,macropore
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   input to intercept the soil moisture model
        PSISM1_vr(L,NY,NX)      = PSISoilMatricP_vr(L,NY,NX)
        VLMicP1_vr(L,NY,NX)     = AZMAX1(VLMicP_vr(L,NY,NX))
        VLWatMicP1_vr(L,NY,NX)  = AZMAX1(VLWatMicP_vr(L,NY,NX))
        VLWatMicPX1_vr(L,NY,NX) = AZMAX1(VLWatMicPX_vr(L,NY,NX))
        VLiceMicP1_vr(L,NY,NX)  = AZMAX1(VLiceMicP_vr(L,NY,NX))
        VLWatMacP1_vr(L,NY,NX)  = AZMAX1(VLWatMacP_vr(L,NY,NX))
        VLiceMacP1_vr(L,NY,NX)  = AZMAX1(VLiceMacP_vr(L,NY,NX))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
          VLairMicP_vr(L,NY,NX)  = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
          VLairMicP1_vr(L,NY,NX) = AZMAX1(VLairMicP_vr(L,NY,NX))
        ELSE
          VLairMicP_vr(L,NY,NX)  = 0.0_r8
          VLairMicP1_vr(L,NY,NX) = 0.0_r8
        ENDIF

        !FVOLAH accounts for clay swelling effect due to change in micropore water, but it is set to zero
        VLMacP1_vr(L,NY,NX)=AZMAX1(VLMacP_vr(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX) &
          *(safe_adb(VLWatMicP1_vr(L,NY,NX),VLSoilMicP_vr(L,NY,NX))-WiltPoint_vr(L,NY,NX))*VGeomLayer_vr(L,NY,NX))

        IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
          VLairMacP_vr(L,NY,NX)  = VLMacP1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
          VLairMacP1_vr(L,NY,NX) = AZMAX1(VLairMacP_vr(L,NY,NX))
        ELSE
          VLairMacP_vr(L,NY,NX)  = 0.0_r8
          VLairMacP1_vr(L,NY,NX) = 0.0_r8
        ENDIF        
        VLWatMicPM_vr(1,L,NY,NX) = VLWatMicP1_vr(L,NY,NX)
        VLWatMacPM(1,L,NY,NX)    = VLWatMacP1_vr(L,NY,NX)
        VLsoiAirPM(1,L,NY,NX)    = VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX)+&
          THETPI*(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))
        
        VLTSoiPore = VLSoilMicP_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)
        IF(VLTSoiPore.GT.ZEROS2(NY,NX))THEN
          !fraction as water
          FracSoiPAsWat_vr(L,NY,NX)=AZMAX1t((VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX))/VLTSoiPore)
          !fraction as ice
          FracSoiPAsIce_vr(L,NY,NX)=AZMAX1t((VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))/VLTSoiPore)
          !fraction as air
          FracSoiPAsAir_vr(L,NY,NX)=AZMAX1t((VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX))/VLTSoiPore)
        ELSE
          FracSoiPAsWat_vr(L,NY,NX)=POROS_vr(L,NY,NX)
          FracSoiPAsIce_vr(L,NY,NX)=0.0_r8
          FracSoiPAsAir_vr(L,NY,NX)=0.0_r8
        ENDIF
        THETPM(1,L,NY,NX)=FracSoiPAsAir_vr(L,NY,NX)
        IF(VLMicP1_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          FracSoilAsAirt(L,NY,NX)=AZMAX1((VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX)) &
            /(VLMicP1_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)))
        ELSE
          FracSoilAsAirt(L,NY,NX)=0.0_r8
        ENDIF
        !VHeatCapacity1_vr=total heat capacity
        !VLHeatCapacityA=heat capcity without macropore water/ice
        !VLHeatCapacityB=heat capacity for macropore water/ice
        VLHeatCapacityA_vr(L,NY,NX)=VHeatCapacitySoilM_vr(L,NY,NX)+cpw*VLWatMicP1_vr(L,NY,NX) &
          +cpi*VLiceMicP1_vr(L,NY,NX)
        if(VLHeatCapacityA_vr(L,NY,NX)>0._r8)VLHeatCapacityA_vr(L,NY,NX)=VLHeatCapacityA_vr(L,NY,NX)+cpo*RootMassElm_vr(ielmc,L,NY,NX)      
        VLHeatCapacityB_vr(L,NY,NX) = cpw*VLWatMacP1_vr(L,NY,NX)+cpi*VLiceMacP1_vr(L,NY,NX)
        VHeatCapacity1_vr(L,NY,NX)  = VLHeatCapacityA_vr(L,NY,NX)+VLHeatCapacityB_vr(L,NY,NX)
    !
    !   MACROPOROSITY
    !
    !   FMAC,SoilFracAsMicP=macropore,micropore volume fractions
    !   CNDH*=macropore hydraulic conductivity
    !   TKS_vr,TK1=soil temperature
    !   FLU,HeatIrrigation=subsurface water,convective heat fluxes
    !   AREAU,AREAD=fractions of layer below natural,artifl water table
    !
        IF(VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          rsatMacP                    = VLMacP1_vr(L,NY,NX)/VLMacP_vr(L,NY,NX)
          SoilFracAsMacP1_vr(L,NY,NX) = SoilFracAsMacP_vr(L,NY,NX)*rsatMacP
          HydroCondMacP1_vr(L,NY,NX)  = HydroCondMacP_vr(L,NY,NX)*rsatMacP**2._r8
        ELSE
          SoilFracAsMacP1_vr(L,NY,NX) = 0.0_r8
          HydroCondMacP1_vr(L,NY,NX)  = 0.0_r8
        ENDIF
        SoilFracAsMicP_vr(L,NY,NX) = 1.0_r8-SoilFracAsMacP1_vr(L,NY,NX)
        TKSoi1_vr(L,NY,NX)         = TKS_vr(L,NY,NX)

        !LyrIrrig=layer number where irrigation is applied
        IF(L.EQ.LyrIrrig)THEN
          FWatIrrigate2MicP_vr(L,NY,NX)  = IrrigSubsurf_col(NY,NX)
          HeatIrrigation(L,NY,NX)        = cpw*TairK_col(NY,NX)*IrrigSubsurf_col(NY,NX)
          FWatIrrigate2MicP1_vr(L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*dts_HeatWatTP
          HeatIrrigation1(L,NY,NX)       = HeatIrrigation(L,NY,NX)*dts_HeatWatTP
        ELSE
          FWatIrrigate2MicP_vr(L,NY,NX)  = 0.0_r8
          HeatIrrigation(L,NY,NX)        = 0.0_r8
          FWatIrrigate2MicP1_vr(L,NY,NX) = 0.0_r8
          HeatIrrigation1(L,NY,NX)       = 0.0_r8
        ENDIF
        IF(CumDepz2LayerBot_vr(L,NY,NX).GE.ExtWaterTable_col(NY,NX))THEN
          AREAU(L,NY,NX)=AMIN1(1.0_r8,AZMAX1(safe_adb(CumDepz2LayerBot_vr(L,NY,NX)-ExtWaterTable_col(NY,NX),DLYR(3,L,NY,NX))))
        ELSE
          AREAU(L,NY,NX)=0.0_r8
        ENDIF
        IF(CumDepz2LayerBot_vr(L,NY,NX).GE.DTBLY(NY,NX))THEN
          AreaUnderWaterTBL(L,NY,NX)=AMIN1(1.0_r8,AZMAX1(safe_adb(CumDepz2LayerBot_vr(L,NY,NX)-DTBLY(NY,NX),DLYR(3,L,NY,NX))))
        ELSE
          AreaUnderWaterTBL(L,NY,NX)=0.0_r8
        ENDIF
      ENDDO D30
    ENDDO DX990
  ENDDO DX995
  
  end subroutine LocalCopySoilVars

!------------------------------------------------------------------------------------------

  subroutine InitSoilHydrauics(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: N,N1,N2,N3,N4,N5,N6
  integer :: NY,NX,L
!     begin_execution
!
!     INITIALIZE SOIL HYDRAULIC PARAMETERS IN LOCAL ARRAYS
!     FOR LATER USE IN WATER TRANSFER ALGORITHMS
!
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      D35: DO L=NUM(NY,NX),NL(NY,NX)
        D40: DO N=FlowDirIndicator(NY,NX),3
          N1=NX;N2=NY;N3=L
! in the EW direction
          IF(N.EQ.iEastWestDirection)THEN
            IF(NX.EQ.NHE)THEN
              cycle
            ELSE
              N4=NX+1
              N5=NY
              N6=L
            ENDIF
! in the NS direction
          ELSEIF(N.EQ.iNorthSouthDirection)THEN
            IF(NY.EQ.NVS)THEN
              cycle
            ELSE
              N4=NX
              N5=NY+1
              N6=L
            ENDIF
! in the vertical direction
          ELSEIF(N.EQ.iVerticalDirection)THEN
            IF(L.EQ.NL(NY,NX))THEN
              cycle
            ELSE
              N4=NX
              N5=NY
              N6=L+1
            ENDIF
          ENDIF
!
    !     MACROPORE CONDUCTIVITY FROM 'HOUR1' AND GRAVITATIONAL
    !     GRADIENT USED TO CALCULATE MACROPORE FLOW FOR USE BELOW
    !
    !     HydroCondMacP1=macropore hydraulic conductivity
    !     AVCNHL=macropore hydraulic conductance
    !     DLYR=layer depth
    !
          IF(HydroCondMacP1_vr(N3,N2,N1).GT.ZERO .AND. HydroCondMacP1_vr(N6,N5,N4).GT.ZERO)THEN
            !when both src and dest have macro pores
            AVCNHL(N,N6,N5,N4)=2.0_r8*HydroCondMacP1_vr(N3,N2,N1)*HydroCondMacP1_vr(N6,N5,N4) &
              /(HydroCondMacP1_vr(N3,N2,N1)*DLYR(N,N6,N5,N4)+HydroCondMacP1_vr(N6,N5,N4) &
              *DLYR(N,N3,N2,N1))
          ELSE
            AVCNHL(N,N6,N5,N4)=0.0_r8
          ENDIF
        ENDDO D40
      ENDDO D35
    ENDDO
  ENDDO
  end subroutine InitSoilHydrauics


!------------------------------------------------------------------------------------------

  subroutine Subsurface3DFlowMit(I,J,M,NHW,NHE,NVN,NVS,KSatRedusByRainKinetEnergy,HeatFluxAir2Soi)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in)  :: M,NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: KSatRedusByRainKinetEnergy(:,:)
  real(r8), dimension(:,:),intent(in) :: HeatFluxAir2Soi(:,:)
  integer :: N,N1,N2,N3,N4,N5,N6,L,LL,K1,KL,NY,NX
  real(r8) :: WTHET1,FCDX,FCLX,FCX
  real(r8) :: PSISV1,TKY,PSDX
  real(r8) :: WPLX,WPX,WTHET2

  real(r8) :: ConvectVapFlux,ConvectHeatFluxMicP,PSISVL
  real(r8) :: TKLX
  real(r8) :: ConvectiveHeatFlxMacP,HeatByWatFlowMicP,HFLWS,THETA1,THETAL  
  logical  :: LInvalidMacP     !disable macropore?
  !     begin_execution
  !
  !     WATER AND ENERGY TRANSFER THROUGH SOIL PROFILE
  !
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      call InitSoil3DModelMit(M,NY,NX)

      LInvalidMacP=.false.
      D4400: DO L=1,NL(NY,NX)
        N1=NX;N2=NY;N3=L
    !
    !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
    !
        D4320: DO N=FlowDirIndicator(N2,N1),3
          IF(N.EQ.iEastWestDirection)THEN
            !west-east direction
            !west to east
            IF(NX.EQ.NHE)THEN
              !east boundary
              cycle
            ELSE
              N4=NX+1
              N5=NY
              N6=L
    !
          !     ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
          !
          !     IF(N2.EQ.2.AND.(N1.EQ.2.OR.N1.EQ.3).AND.L.LE.15)THEN
          !     CYCLE
          !     ENDIF
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN
            !north-south direction
            !north to south
            IF(NY.EQ.NVS)THEN
              !south boundary
              cycle
            ELSE
              N4=NX
              N5=NY+1
              N6=L
    !
          !     ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
          !
          !     IF(N1.EQ.3.AND.(N2.EQ.1.OR.N2.EQ.2).AND.L.LE.15)THEN
          !     CYCLE
          !     ENDIF
          !
          !     END ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
          !
            ENDIF
          ELSEIF(N.EQ.iVerticalDirection)THEN
            !vertical direction
            IF(L.EQ.NL(NY,NX))THEN
              cycle
            ELSE
              !destination is layer below
              N4=NX
              N5=NY
              N6=L+1
            ENDIF
          ENDIF
!
!     SKIP NON-EXISTENT DESTINATION SOIL LAYERS
          ! identified by soil volume
          D1100: DO LL=N6,NL(NY,NX)
            IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
              N6=LL
              exit
            ENDIF
          ENDDO D1100
          ! need to better understand the following line
          IF(N3.EQ.NU(N2,N1))N6X(N2,N1)=N6
      !
      !     POROSITIES 'THETP*', WATER CONTENTS 'THETA*', AND POTENTIALS
      !     'PSIS*' FOR EACH GRID CELL
      !
      !     THETA1,THETAL=micropore water concn in source,destination cells
      !     THETY=hygroscopic water concentration
      !     POROS=soil porosity
      !     VLSoilPoreMicP_vrI=soil volume excluding rock, macropore
      !
          IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4).AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
              !within the calculation domain
    !          if(I>=132 .and. j==19)print*,' source layer'
              call CalcSoilWatPotential(NY,NX,N1,N2,N3,PSISoilMatricPtmp_vr(N3,N2,N1),THETA1)
   !           if(I>=132 .and. j==19)print*,'dest'
              call CalcSoilWatPotential(NY,NX,N4,N5,N6,PSISoilMatricPtmp_vr(N6,N5,N4),THETAL)
              !
              !     ACCOUNT FOR WETTING FRONTS WHEN CALCULATING WATER CONTENTS,
              !     MATRIC WATER POTENTIALS AND HYDRAULIC CONDUCTIVITIES USED
              !     IN WATER FLUX CALCULATIONS
              !
              !     HydCondSrc,CNDL=hydraulic conductivities in source,destination cells
              !     KSatRedusByRainKinetEnergy=reduction in soil surface Ksat from rainfall energy impact
              !     PSISM1=soil matric potential
              !     VLWatMicPX1=VLWatMicP1 accounting for wetting front
              !
  !            if(I>=132 .and. j==19)print*,'DARCY FLOW IF BOTH CELLS ARE SATURATED'
              !     (CURRENT WATER POTENTIAL > AIR ENTRY WATER POTENTIAL)
              !
              call MicporeDarcyFlow(NY,NX,N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,&
                KSatRedusByRainKinetEnergy(NY,NX),HeatByWatFlowMicP,PSISV1,PSISVL)          

          !
              call MacporeFLow(NY,NX,M,N,N1,N2,N3,N4,N5,N6,ConvectiveHeatFlxMacP,LInvalidMacP)

              call WaterVaporFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectVapFlux,&
                ConvectHeatFluxMicP)

          !
          !     FLWL=total water+vapor flux to destination
          !     WatXChange2WatTableX=total unsaturated water+vapor flux to destination
          !     HeatByWatFlowMicP=total convective heat flux from water+vapor flux
          !
              WatXChange2WatTable(N,N6,N5,N4)  = WatXChange2WatTable(N,N6,N5,N4)+ConvectVapFlux
              WatXChange2WatTableX(N,N6,N5,N4) = WatXChange2WatTableX(N,N6,N5,N4)+ConvectVapFlux
              HeatByWatFlowMicP                = HeatByWatFlowMicP+ConvectHeatFluxMicP
              HeatFlow2Soili(N,N6,N5,N4)       = HeatByWatFlowMicP+ConvectiveHeatFlxMacP
          !
              call Solve4Heat(I,J,N,NY,NX,N1,N2,N3,N4,N5,N6,ConvectHeatFluxMicP,HeatFluxAir2Soi(NY,NX))

          !
          !     TOTAL WATER, VAPOR AND HEAT FLUXES
          !
          !     FLW,FLWX,FLWH=total water flux through micropores,macropores
          !     HFLW=total heat flux
          !     WaterFlow2MicPM=water flux used for solute flux calculations in TranspNoSalt.f
          !
              WaterFlowSoiMicP_3D(N,N6,N5,N4) = WaterFlowSoiMicP_3D(N,N6,N5,N4)+WatXChange2WatTable(N,N6,N5,N4)
              WaterFlowSoiMicPX(N,N6,N5,N4)   = WaterFlowSoiMicPX(N,N6,N5,N4)+WatXChange2WatTableX(N,N6,N5,N4)
              WaterFlowMacP_3D(N,N6,N5,N4)    = WaterFlowMacP_3D(N,N6,N5,N4)+ConvWaterFlowMacP_3D(N,N6,N5,N4)
              HeatFlow2Soil_3D(N,N6,N5,N4)    = HeatFlow2Soil_3D(N,N6,N5,N4)+HeatFlow2Soili(N,N6,N5,N4)
              WaterFlow2MicPM(M,N,N6,N5,N4)   = WatXChange2WatTable(N,N6,N5,N4)

              IF(N.EQ.iVerticalDirection)THEN
            !
            !     WATER FILM THICKNESS FOR CALCULATING GAS EXCHANGE IN TranspNoSalt.F
            !
                FILM(M,N6,N5,N4)=FilmThickness(PSISoilMatricPtmp_vr(N6,N5,N4))
              ENDIF
            ELSEIF(N.NE.iVerticalDirection)THEN
              WatXChange2WatTable(N,N6,N5,N4)  = 0.0_r8
              WatXChange2WatTableX(N,N6,N5,N4) = 0.0_r8
              ConvWaterFlowMacP_3D(N,N6,N5,N4) = 0.0_r8
              HeatFlow2Soili(N,N6,N5,N4)       = 0.0_r8
              WaterFlow2MicPM(M,N,N6,N5,N4)    = 0.0_r8
              WaterFlow2MacPM(M,N,N6,N5,N4)    = 0.0_r8
            ENDIF

          ELSE
            IF(N.EQ.iVerticalDirection)THEN
              WatXChange2WatTable(N,N3,N2,N1)  = 0.0_r8
              WatXChange2WatTableX(N,N3,N2,N1) = 0.0_r8
              ConvWaterFlowMacP_3D(N,N3,N2,N1) = 0.0_r8
              HeatFlow2Soili(N,N3,N2,N1)       = 0.0_r8
              WaterFlow2MacPM(M,N,N3,N2,N1)    = 0.0_r8
              WaterFlow2MacPM(M,N,N3,N2,N1)    = 0.0_r8
            ELSE
              WatXChange2WatTable(N,N6,N5,N4)  = 0.0_r8
              WatXChange2WatTableX(N,N6,N5,N4) = 0.0_r8
              ConvWaterFlowMacP_3D(N,N6,N5,N4) = 0.0_r8
              HeatFlow2Soili(N,N6,N5,N4)       = 0.0_r8
              WaterFlow2MicPM(M,N,N6,N5,N4)    = 0.0_r8
              WaterFlow2MacPM(M,N,N6,N5,N4)    = 0.0_r8
            ENDIF
          ENDIF
        ENDDO D4320
      ENDDO D4400
    ENDDO
  ENDDO  
  end subroutine Subsurface3DFlowMit
!------------------------------------------------------------------------------------------

  subroutine LateralWatHeatExchMit(I,J,M,NHW,NHE,NVN,NVS,KSatRedusByRainKinetEnergyS)
  !
  !Description
  ! boundary flow involes exchange with external water table, and through tile drainage
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: KSatRedusByRainKinetEnergyS(JY,JX)
  integer :: NY,NX
  integer :: L,LL
  integer :: N,NN,N1,N2,N3,N4,N5,N4B,N5B,N6
  integer :: M1,M2,M3,M4,M5,M6,K1,KL
  integer :: IFLGDH,IFLGU,IFLGUH,IFLGD
  real(r8) :: HydcondSrc,VLWatMicP1X
  real(r8) :: FLWT
  real(r8) :: RechargSurf,RechargSubSurf,RechargRateWTBL
  real(r8) :: DPTHH,CNDL
  real(r8) :: FINHX,THETAX  
  real(r8) :: AirfMicP,VOLPX2,AirfMacP
  real(r8) :: XN,THETA1
  real(r8) :: VOLP1X,VLWatMacP1X,VOLPH1X

!     begin_execution
!     AirfMicP,AirfMacP=air-filled porosity in micropores,macropores
  D9595: DO  NX=NHW,NHE
    D9590: DO  NY=NVN,NVS
      D9585: DO L=NUM(NY,NX),NL(NY,NX)
        AirfMicP = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
        VOLPX2   = AirfMicP
        AirfMacP = VLMacP1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
!
        call Config4WaterTableDrain(L,NY,NX,IFLGU,IFLGUH,DPTHH)

!
!       'IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO TILE DRAIN'
        call Config4TileDrainage(L,NY,NX,IFLGD,IFLGDH,DPTHH)
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
!     N3,N2,N1=L,NY,NX of source grid cell
!     M6,M5,M4=L,NY,NX of destination grid cell
!       
        N1=NX;N2=NY;N3=L
!
!     LOCATE EXTERNAL BOUNDARIES
!
        D9580: DO N=FlowDirIndicator(NY,NX),3
          D9575: DO NN=1,2
            IF(N.EQ.iEastWestDirection)THEN
              ! along the W-E direction
              N4=NX+1
              N5=NY
              N4B=NX-1
              N5B=NY
              N6=L
              IF(NN.EQ.1)THEN
                !eastern boundary
                IF(NX.EQ.NHE)THEN                
                  M1              = NX
                  M2              = NY
                  M3              = L
                  M4              = NX+1
                  M5              = NY
                  M6              = L
                  XN              = -1.0_r8   !going out
                  RechargSurf     = RechargEastSurf(M2,M1)
                  RechargSubSurf  = RechargEastSubSurf(M2,M1)
                  RechargRateWTBL = RechargRateEastWTBL(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                !west boundary
                IF(NX.EQ.NHW)THEN
                  M1              = NX+1
                  M2              = NY
                  M3              = L
                  M4              = NX
                  M5              = NY
                  M6              = L
                  XN              = 1.0_r8    !coming in
                  RechargSurf     = RechargWestSurf(M5,M4)
                  RechargSubSurf  = RechargWestSubSurf(M5,M4)
                  RechargRateWTBL = RechargRateWestWTBL(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              ! along the N-S direction
              N4=NX
              N5=NY+1
              N4B=NX
              N5B=NY-1
              N6=L
              IF(NN.EQ.1)THEN
                !south boundary
                IF(NY.EQ.NVS)THEN
                  M1              = NX
                  M2              = NY
                  M3              = L
                  M4              = NX
                  M5              = NY+1
                  M6              = L
                  XN              = -1.0_r8    !going out
                  RechargSurf     = RechargSouthSurf(M2,M1)
                  RechargSubSurf  = RechargSouthSubSurf(M2,M1)
                  RechargRateWTBL = RechargRateSouthWTBL(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                !north boundary
                IF(NY.EQ.NVN)THEN
                  M1              = NX
                  M2              = NY+1
                  M3              = L
                  M4              = NX
                  M5              = NY
                  M6              = L
                  XN              = 1.0_r8   !coming in
                  RechargSurf     = RechargNorthSurf(M5,M4)
                  RechargSubSurf  = RechargNorthSubSurf(M5,M4)
                  RechargRateWTBL = RechargRateNorthWTBL(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN
! in the vertical direction
              N4=NX
              N5=NY
              N6=L+1
              IF(NN.EQ.1)THEN
                !bottom
                IF(L.EQ.NL(NY,NX))THEN
                  M1              = NX
                  M2              = NY
                  M3              = L
                  M4              = NX
                  M5              = NY
                  M6              = L+1
                  XN              = -1.0_r8    !going out
                  RechargSubSurf  = RCHGD(M2,M1)
                  RechargRateWTBL = 1.0_r8
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                cycle
              ENDIF
            ENDIF
!
!     BOUNDARY SURFACE RUNOFF DEPENDING ON ASPECT, SLOPE
!     VELOCITY, HYDRAULIC RADIUS AND SURFACE WATER STORAGE
!
!     CDPTH,CumSoilDeptht0=current,initial surface elevation
!     BKDS=bulk density
!     XGridRunoffFlag,RCHQ*=runoff boundary flags

!           top soil layer and surface soil layer is active, litter layer is lower than its initial thickness 
!           or the grid is a soil
!           Not in vertical direction, 
            IF(L.EQ.NUM(N2,N1) .AND. N.NE.iVerticalDirection &
              .AND. (CumDepz2LayerBot_vr(NU(N2,N1)-1,N2,N1).LE.CumSoilDeptht0(N2,N1) &
              .OR. SoiBulkDensity_vr(NUI(N2,N1),N2,N1).GT.ZERO))THEN
!  NO runoff
              IF(.not.XGridRunoffFlag(NN,N,N2,N1).OR.isclose(RechargSurf,0._r8).OR. &
                ABS(WatFlux4ErosionM_2DH(M,N2,N1)).LT.ZEROS(N2,N1))THEN
                WatFlx2LitRByRunoff(N,NN,M5,M4)  = 0.0_r8
                HeatFlx2LitRByRunoff(N,NN,M5,M4) = 0.0_r8
              ELSE
                call SurfaceRunoff(M,N,NN,N1,N2,M4,M5,RechargSurf,XN)
!
        !     BOUNDARY SNOW FLUX
        !
        !     DrySnoFlxBySnowRedistribut,WatFlxBySnowRedistribut,IceFlxBySnowRedistribut=snow,water,ice transfer
        !     HeatFlxBySnowRedistribut=convective heat transfer from snow,water,ice transfer
        !     QS,QW,QI=cumulative hourly snow,water,ice transfer
        !     HQS=cumulative hourly convective heat transfer from snow,water,ice transfer
!
                IF(NN.EQ.1)THEN
                  DrySnoFlxBySnowRedistribut(N,M5,M4)  = 0.0_r8
                  WatFlxBySnowRedistribut(N,M5,M4)     = 0.0_r8
                  IceFlxBySnowRedistribut(N,M5,M4)     = 0.0_r8
                  HeatFlxBySnowRedistribut(N,M5,M4)    = 0.0_r8
                  DrySnoFlxBySnoRedistM_2DH(M,N,M5,M4) = DrySnoFlxBySnowRedistribut(N,M5,M4)
                ENDIF
              ENDIF
            ELSE
              IF(N.NE.iVerticalDirection)THEN
                WatFlx2LitRByRunoff(N,NN,M5,M4)  = 0.0_r8
                HeatFlx2LitRByRunoff(N,NN,M5,M4) = 0.0_r8
              ENDIF
            ENDIF
!
          ! BOUNDARY SUBSURFACE WATER AND HEAT TRANSFER DEPENDING
          ! ON LEVEL OF WATER TABLE
!
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
              IF(FlowDirIndicator(N2,N1).NE.iVerticalDirection .OR. N.EQ.iVerticalDirection)THEN
              !including lateral connection or woking on vertical direction
!
              ! IF NO WATER TABLE
              !
              ! IDWaterTable=water table flag
              ! THETA1,THETAX=water content ahead,behind wetting front
              ! K1,KL=pore water class ahead,behind wetting front
              ! HydcondSrc,CNDL=hydraulic conductivity ahead,behind wetting front
              ! KSatRedusByRainKinetEnergy=reduction in soil surface Ksat from rainfall energy impact
              ! FLWL,WatXChange2WatTableX=lower boundary micropore water flux
              ! ConvWaterFlowMacP_3D=lower boundary macropore water flux
              ! HFLWL=convective heat from lower boundary water flux
              ! XH,XN,dts_HeatWatTP=rate constant,direction indicator,time step
              ! SLOPE=sin(vertical slope)=1
              ! RCHG*=boundary flags
!
                IF(IDWaterTable(N2,N1).EQ.0 .OR. N.EQ.iVerticalDirection)THEN              
                  !involve no water table or vertical direction
                  THETA1 = AMAX1(THETY_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),&
                    safe_adb(VLWatMicP1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
                  THETAX = AMAX1(THETY_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),&
                    safe_adb(VLWatMicPX1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
                  K1     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-THETA1)/POROS_vr(N3,N2,N1))+1))
                  KL     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-THETAX)/POROS_vr(N3,N2,N1))+1))

                  IF(N3.EQ.NUM(NY,NX))THEN
                   HydcondSrc=HydroCond_3D(N,K1,N3,N2,N1)*KSatRedusByRainKinetEnergyS(NY,NX)
                  ELSE
                   HydcondSrc=HydroCond_3D(N,K1,N3,N2,N1)
                  ENDIF
                
                  WatXChange2WatTable(N,M6,M5,M4)=AMIN1(VLWatMicP1_vr(N3,N2,N1)*dts_wat, &
                    XN*mGravAccelerat*(-ABS(SLOPE(N,N2,N1)))*HydcondSrc*AREA(3,N3,N2,N1)) &
                    *RechargSubSurf*RechargRateWTBL*dts_HeatWatTP

                  WatXChange2WatTableX(N,M6,M5,M4) = WatXChange2WatTable(N,M6,M5,M4)
                  ConvWaterFlowMacP_3D(N,M6,M5,M4) = AMIN1(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
                    ,XN*mGravAccelerat*(-ABS(SLOPE(N,N2,N1)))*HydroCondMacP1_vr(N3,N2,N1)*AREA(3,N3,N2,N1)) &
                    *RechargSubSurf*RechargRateWTBL*dts_HeatWatTP
                  HeatFlow2Soili(N,M6,M5,M4)=cpw*TKSoi1_vr(N3,N2,N1)*(WatXChange2WatTable(N,M6,M5,M4) &
                    +ConvWaterFlowMacP_3D(N,M6,M5,M4))

                ELSE
!
                  CALL WaterTBLDrain(N,N1,N2,N3,M4,M5,M6,IFLGU,IFLGUH,RechargSubSurf,RechargRateWTBL,DPTHH,XN)

                  call TileDrain(N,N1,N2,N3,M4,M5,M6,IFLGD,IFLGDH,RechargRateWTBL,RechargSubSurf,DPTHH,XN)

                  call SubSufXChangeWithExtWatTable(NY,NX,N,N1,N2,N3,M4,M5,M6,DPTHH,RechargSubSurf,&
                    RechargRateWTBL,XN,AirfMicP,VOLPX2,AirfMacP)

                ENDIF
!
        !     SUBSURFACE HEAT SOURCE/SINK
        !
        !     HFLWL=heat flux across lower boundary
        !     TK1=lower boundary soil temperature
        !     TKSD=deep source/sink temperature from geothermal flux
        !     TCNDG=thermal conductivity below lower boundary
        !     SoilHeatSrcDepth,CDPTH=depth of thermal sink/source, lower boundary
        !     KoppenClimZone=Koppen climate zone
                IF(N.EQ.iVerticalDirection .AND. KoppenClimZone_col(N2,N1).NE.-2)THEN
                  HeatFlow2Soili(N,M6,M5,M4)=HeatFlow2Soili(N,M6,M5,M4)+(TKSoi1_vr(N3,N2,N1)-TKSD(N2,N1))* &
                    TCNDG/(SoilHeatSrcDepth(N2,N1)-CumDepz2LayerBot_vr(N3,N2,N1)) &
                    *AREA(N,N3,N2,N1)*dts_HeatWatTP
                ENDIF
                WaterFlowSoiMicP_3D(N,M6,M5,M4) = WaterFlowSoiMicP_3D(N,M6,M5,M4)+WatXChange2WatTable(N,M6,M5,M4)
                WaterFlowSoiMicPX(N,M6,M5,M4)   = WaterFlowSoiMicPX(N,M6,M5,M4)+WatXChange2WatTableX(N,M6,M5,M4)
                WaterFlowMacP_3D(N,M6,M5,M4)    = WaterFlowMacP_3D(N,M6,M5,M4)+ConvWaterFlowMacP_3D(N,M6,M5,M4)
                HeatFlow2Soil_3D(N,M6,M5,M4)    = HeatFlow2Soil_3D(N,M6,M5,M4)+HeatFlow2Soili(N,M6,M5,M4)
                WaterFlow2MicPM(M,N,M6,M5,M4)   = WatXChange2WatTable(N,M6,M5,M4)
                WaterFlow2MacPM(M,N,M6,M5,M4)   = ConvWaterFlowMacP_3D(N,M6,M5,M4)
              ENDIF
            ELSE
              WatXChange2WatTable(N,M6,M5,M4)  = 0.0_r8
              WatXChange2WatTableX(N,M6,M5,M4) = 0.0_r8
              ConvWaterFlowMacP_3D(N,M6,M5,M4) = 0.0_r8
              HeatFlow2Soili(N,M6,M5,M4)       = 0.0_r8
              WaterFlow2MicPM(M,N,M6,M5,M4)    = 0.0_r8
              WaterFlow2MacPM(M,N,M6,M5,M4)    = 0.0_r8
            ENDIF
          ENDDO D9575
    !
    !     NET WATER AND HEAT FLUXES IN RUNOFF AND SNOW DRIFT
    !
    !     TQR1,THQR1=net runoff,convective heat from runoff
    !     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
    !     WatFlx2LitRByRunoff,HeatFlx2LitRByRunoff=runoff, convective heat from runoff
    !     DrySnoFlxBySnowRedistribut,WatFlxBySnowRedistribut,IceFlxBySnowRedistribut=snow,water,ice transfer
    !     HeatFlxBySnowRedistribut=convective heat transfer from snow,water,ice transfer
!        if(I>=132 .and. j==19)print*,'top layer snow redistribution'!
          IF(L.EQ.NUM(N2,N1).AND.N.NE.iVerticalDirection)THEN
            if(snowRedist_model)call SumSnowDriftByRunoff(M,N,N1,N2,N4,N5,N4B,N5B)
          ENDIF
!
        !     NET WATER AND HEAT FLUXES THROUGH SOIL AND SNOWPACK
        !
        !     TFLWL,THFLWL=net water micropore,macropore flux
        !     THFLWL=net convective+conductive heat flux
        !     FLWL =micropore water,heat flux
        !     ConvWaterFlowMacP_3D=macropore water,heat flux
        !     HFLWL=soil heat flux
!
          IF(FlowDirIndicator(N2,N1).NE.iVerticalDirection.OR.N.EQ.iVerticalDirection)THEN
            D1200: DO LL=N6,NL(N5,N4)
              IF(VLSoilPoreMicP_vr(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
                N6=LL
                exit
              ENDIF
            ENDDO D1200
            !exchange with water table if micropore is non-zero
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
              TWatCharge2MicP(N3,N2,N1)          = TWatCharge2MicP(N3,N2,N1)+WatXChange2WatTable(N,N3,N2,N1) &
                -WatXChange2WatTable(N,N6,N5,N4)
              TWatXChange2WatTableX(N3,N2,N1)    = TWatXChange2WatTableX(N3,N2,N1)+WatXChange2WatTableX(N,N3,N2,N1) &
                -WatXChange2WatTableX(N,N6,N5,N4)
              TConvWaterFlowMacP_3D_vr(N3,N2,N1) = TConvWaterFlowMacP_3D_vr(N3,N2,N1)+ConvWaterFlowMacP_3D(N,N3,N2,N1) &
                -ConvWaterFlowMacP_3D(N,N6,N5,N4)
              THeatFlow2Soili_vr(N3,N2,N1)       = THeatFlow2Soili_vr(N3,N2,N1)+HeatFlow2Soili(N,N3,N2,N1) &
                -HeatFlow2Soili(N,N6,N5,N4)
            ELSE
              TWatCharge2MicP(N3,N2,N1)          = 0.0_r8
              TWatXChange2WatTableX(N3,N2,N1)    = 0.0_r8
              TConvWaterFlowMacP_3D_vr(N3,N2,N1) = 0.0_r8
              THeatFlow2Soili_vr(N3,N2,N1)       = 0.0_r8
            ENDIF
          ENDIF
        ENDDO D9580
!
!     INFILTRATION OF WATER FROM MACROPORES INTO MICROPORES
!
!     VLWatMacP1=macropore volume
!     FINHX,FINHL=macro-micropore transfer unltd,ltd by water,air volume
!     FWatExMacP2MicPM=macro-micropore transfer for use in TranspNoSalt.f
!     HydroCond_3D=hydraulic conductivity
!     PSISE,PSISoilAirEntry=air entry,matric water potentials
!     PHOL,MacPRadius=path length between,radius of macropores from hour1.f
!     dts_HeatWatTP=time step
!     VLWatMicP1X,VOLP1X=current micropore water,air volume
!     VLWatMacP1X,VOLPH1X=current macropore water,air volume
!
        IF(VLWatMacP1_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
          !south-north direction, is it true?
          FINHX = TwoPiCON*HydroCond_3D(2,1,N3,N2,N1)*AREA(3,N3,N2,N1) &
            *(PSISE_vr(N3,N2,N1)-PSISoilMatricPtmp_vr(N3,N2,N1)) &
            /LOG(PathLenMacP(N3,N2,N1)/MacPRadius(N3,N2,N1))*dts_HeatWatTP
          VLWatMicP1X = VLWatMicP1_vr(N3,N2,N1)+TWatCharge2MicP(N3,N2,N1)+FWatIrrigate2MicP1_vr(N3,N2,N1)
          VOLP1X      = AZMAX1(VLMicP1_vr(N3,N2,N1)-VLWatMicP1X-VLiceMicP1_vr(N3,N2,N1))
          VLWatMacP1X = VLWatMacP1_vr(N3,N2,N1)+TConvWaterFlowMacP_3D_vr(N3,N2,N1)
          VOLPH1X     = AZMAX1(VLMacP1_vr(N3,N2,N1)-VLWatMacP1X-VLiceMacP1_vr(N3,N2,N1))

          IF(FINHX.GT.0.0_r8)THEN
            FWatExMacP2MicPi(N3,N2,N1)=AZMAX1(AMIN1(FINHX,VLWatMacP1X,VOLP1X))
          ELSE
            FWatExMacP2MicPi(N3,N2,N1)=AZMIN1(AMAX1(FINHX,-VOLPH1X,-VLWatMicP1X))
          ENDIF
          FWatExMacP2MicPM(M,N3,N2,N1) = FWatExMacP2MicPi(N3,N2,N1)
          FWatExMacP2MicP(N3,N2,N1)    = FWatExMacP2MicP(N3,N2,N1)+FWatExMacP2MicPi(N3,N2,N1)
        ELSE
          FWatExMacP2MicPi(N3,N2,N1)   = 0.0_r8
          FWatExMacP2MicPM(M,N3,N2,N1) = 0.0_r8
        ENDIF

        call FreezeThawMit(NY,NX,L,N1,N2,N3)
!
!     DISSIPATE WETTING FRONT
!
!     VLWatMicP1=soil micropore water content
!     VLWatMicPX1=soil micropore water content behind wetting front
!     FLWVL=water flux from wetted to drier soil
!
!        TLPhaseChangeHeat2Soi1s(L,NY,NX)=TLPhaseChangeHeat2Soi1s(L,NY,NX)+TLPhaseChangeHeat2Soi1(L,NY,NX)
      ENDDO D9585
    ENDDO D9590
  ENDDO D9595

  end subroutine LateralWatHeatExchMit
!------------------------------------------------------------------------------------------
  subroutine UpdateSoilMoistTemp(I,J,M,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX,L
  real(r8) :: tk1pres,tk1l
  real(r8) :: ENGY1,VLTSoiPore,VHXX
  real(r8) :: TKXX,VLWMicPre,VLHeatCapacityPre

! begin_execution
  D11: DO NX=NHW,NHE
    D12: DO NY=NVN,NVS
      D13: DO L=NUM(NY,NX),NL(NY,NX)

        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VLWMicPre              = VLWatMicP1_vr(L,NY,NX)
          VLWatMicP1_vr(L,NY,NX) = VLWatMicP1_vr(L,NY,NX)+TWatCharge2MicP(L,NY,NX)+FWatExMacP2MicPi(L,NY,NX) &
            +TMLiceThawMicP(L,NY,NX)+FWatIrrigate2MicP1_vr(L,NY,NX)
          
          VLWatMicPX1_vr(L,NY,NX) = VLWatMicPX1_vr(L,NY,NX)+TWatXChange2WatTableX(L,NY,NX)+FWatExMacP2MicPi(L,NY,NX) &
             +TMLiceThawMicP(L,NY,NX)+FWatIrrigate2MicP1_vr(L,NY,NX)
          VLWatMicPX1_vr(L,NY,NX)   = AMIN1(VLWatMicP1_vr(L,NY,NX),VLWatMicPX1_vr(L,NY,NX))
          VLiceMicP1_vr(L,NY,NX) = VLiceMicP1_vr(L,NY,NX)-TMLiceThawMicP(L,NY,NX)/DENSICE

          VLWatMacP1_vr(L,NY,NX) = VLWatMacP1_vr(L,NY,NX)+TConvWaterFlowMacP_3D_vr(L,NY,NX)-FWatExMacP2MicPi(L,NY,NX) &
            +TMLiceThawMacP(L,NY,NX)
          VLiceMacP1_vr(L,NY,NX) = VLiceMacP1_vr(L,NY,NX)-TMLiceThawMacP(L,NY,NX)/DENSICE

          IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
            ! air-filled space
            VLairMicP_vr(L,NY,NX)  = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
            VLairMicP1_vr(L,NY,NX) = AZMAX1(VLairMicP_vr(L,NY,NX))
            VLairMacP_vr(L,NY,NX)  = VLMacP1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
            VLairMacP1_vr(L,NY,NX) = AZMAX1(VLairMacP_vr(L,NY,NX))
            VLMacP1_vr(L,NY,NX)    = AZMAX1(VLMacP_vr(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX) &
              *(safe_adb(VLWatMicP1_vr(L,NY,NX),VLSoilMicP_vr(L,NY,NX))-WiltPoint_vr(L,NY,NX))*VGeomLayer_vr(L,NY,NX))
          ELSE
            VLairMicP_vr(L,NY,NX)  = 0.0_r8
            VLairMicP1_vr(L,NY,NX) = 0.0_r8
            VLairMacP_vr(L,NY,NX)  = 0.0_r8
            VLairMacP1_vr(L,NY,NX) = 0.0_r8
            VLMicP1_vr(L,NY,NX)    = VLWatMicP1_vr(L,NY,NX)+VLiceMicP1_vr(L,NY,NX)
            VLMacP1_vr(L,NY,NX)    = 0.0_r8
          ENDIF

          VLTSoiPore                  = VLSoilMicP_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)
          FracSoiPAsWat_vr(L,NY,NX)   = AZMAX1t((VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX))/VLTSoiPore)
          FracSoiPAsIce_vr(L,NY,NX)   = AZMAX1t((VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))/VLTSoiPore)
          FracSoiPAsAir_vr(L,NY,NX)   = AZMAX1t((VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX))/VLTSoiPore)
          
          IF(VLMicP1_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            FracSoilAsAirt(L,NY,NX)=AZMAX1((VLairMicP1_vr(L,NY,NX)+ &
              VLairMacP1_vr(L,NY,NX))/(VLMicP1_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)))
          ELSE
            FracSoilAsAirt(L,NY,NX)=0.0_r8
          ENDIF
          IF(VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            SoilFracAsMacP1_vr(L,NY,NX) = SoilFracAsMacP_vr(L,NY,NX)*VLMacP1_vr(L,NY,NX)/VLMacP_vr(L,NY,NX)
            HydroCondMacP1_vr(L,NY,NX)  = HydroCondMacP_vr(L,NY,NX)*(VLMacP1_vr(L,NY,NX)/VLMacP_vr(L,NY,NX))**2
          ELSE
            SoilFracAsMacP1_vr(L,NY,NX)   = 0.0_r8
            HydroCondMacP1_vr(L,NY,NX) = 0.0_r8
          ENDIF
          SoilFracAsMicP_vr(L,NY,NX) = 1.0_r8-SoilFracAsMacP1_vr(L,NY,NX)
          TKXX                    = TKSoi1_vr(L,NY,NX)
          VHXX                    = VHeatCapacity1_vr(L,NY,NX)
          ENGY1                   = VHeatCapacity1_vr(L,NY,NX)*TKSoi1_vr(L,NY,NX)

          if(TKSoi1_vr(L,NY,NX)>400._r8.or.TKSoi1_vr(L,NY,NX)<100._r8)then
            write(*,*)'======'
            write(*,*)'VLWatMicP1_vr(L,NY,NX)=',VLWMicPre,VLWatMicP1_vr(L,NY,NX),L,VLairMicP_vr(L,NY,NX),&
              VLiceMicP1_vr(L,NY,NX)
            write(*,*)TWatCharge2MicP(L,NY,NX),FWatExMacP2MicPi(L,NY,NX), &
              TMLiceThawMicP(L,NY,NX),FWatIrrigate2MicP1_vr(L,NY,NX)
            write(*,*)'M, L=',M,L,NY,NX,NUM(NY,NX),VGeomLayer_vr(L,NY,NX),ZEROS2(NY,NX),VLWatMicPX1_vr(L,NY,NX)
            write(*,*)'SoiBulkDensity_vr(L,NY,NX)=',SoiBulkDensity_vr(L,NY,NX),TKS_vr(L,NY,NX),ZEROS(NY,NX)
            write(*,*)'VHeatCapacity1_vr(L,NY,NX),TKSoi1_vr(L,NY,NX),TKXX',L,VHeatCapacity1_vr(L,NY,NX),TKSoi1_vr(L,NY,NX),TKXX
            write(*,*)VLMicP1_vr(L,NY,NX),VLMacP1_vr(L,NY,NX)
            if(TKSoi1_vr(L,NY,NX)>1.e3_r8.or.TKSoi1_vr(L,NY,NX)<0._r8)call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          
          VLHeatCapacityPre           = VHeatCapacity1_vr(L,NY,NX)
          VLHeatCapacityA_vr(L,NY,NX) = VHeatCapacitySoilM_vr(L,NY,NX)+cpw*VLWatMicP1_vr(L,NY,NX)+cpi*VLiceMicP1_vr(L,NY,NX)
          if(VLHeatCapacityA_vr(L,NY,NX)>0._r8)VLHeatCapacityA_vr(L,NY,NX)=VLHeatCapacityA_vr(L,NY,NX)+cpo*RootMassElm_vr(ielmc,L,NY,NX)
          VLHeatCapacityB_vr(L,NY,NX) = cpw*VLWatMacP1_vr(L,NY,NX)+cpi*VLiceMacP1_vr(L,NY,NX)
          VHeatCapacity1_vr(L,NY,NX)  = VLHeatCapacityA_vr(L,NY,NX)+VLHeatCapacityB_vr(L,NY,NX)

         IF(VHeatCapacity1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

            tk1l            = TKSoi1_vr(L,NY,NX)       
            TKSoi1_vr(L,NY,NX) = (ENGY1+THeatFlow2Soili_vr(L,NY,NX)+TLPhaseChangeHeat2Soi1(L,NY,NX)+&
              HeatIrrigation1(L,NY,NX))/VHeatCapacity1_vr(L,NY,NX)
            
            if(TKSoi1_vr(L,NY,NX)>400.0_r8)then
              write(*,*)'TKSoi1_vr(L,NY,NX)',M,L,TKSoi1_vr(L,NY,NX),tk1l,VLHeatCapacityPre,VHeatCapacity1_vr(L,NY,NX)
              write(*,*)ENGY1/VHeatCapacity1_vr(L,NY,NX),THeatFlow2Soili_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX),&
               TLPhaseChangeHeat2Soi1(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX),HeatIrrigation1(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX)
              write(*,*)VLWatMicPM_vr(M,L,NY,NX),VLWatMicPM_vr(M+1,L,NY,NX),THeatFlow2Soili_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX)
              call endrun()
            endif
          ELSEIF(L.EQ.1)THEN
            TKSoi1_vr(L,NY,NX)=TairK_col(NY,NX)
          ELSE
            TKSoi1_vr(L,NY,NX)=TKSoi1_vr(L-1,NY,NX)
          ENDIF
        ENDIF
      ENDDO D13
    ENDDO D12
  ENDDO D11

  end subroutine UpdateSoilMoistTemp

!------------------------------------------------------------------------------------------

  subroutine UpdateStateFluxAtM(M,NHW,NHE,NVN,NVS)
  !
  !Description
  ! Early exit of the watsub solver
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX
  integer :: L,NUX,LL,Ls
  
! begin_execution
  D9795: DO NX=NHW,NHE
    D9790: DO NY=NVN,NVS

      !
      ! SOIL LAYER WATER, ICE AND TEMPERATURE
      !
      D9785: DO L=NUM(NY,NX),NL(NY,NX)
        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

          !record intermediate variables for bgc calculation
          VLWatMicPM_vr(M+1,L,NY,NX) = VLWatMicP1_vr(L,NY,NX)
          VLWatMacPM(M+1,L,NY,NX)    = VLWatMacP1_vr(L,NY,NX)
          VLsoiAirPM(M+1,L,NY,NX)    = VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX) &
            +THETPI*(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))

          !change in soil air volume
          ReductVLsoiAirPM(M,L,NY,NX) = VLsoiAirPM(M,L,NY,NX)-VLsoiAirPM(M+1,L,NY,NX)

          THETPM(M+1,L,NY,NX)         = FracSoiPAsAir_vr(L,NY,NX)
 !
        ELSE
          !layer L disappears
          VLWatMicPM_vr(M+1,L,NY,NX)  = 0.0_r8
          VLWatMacPM(M+1,L,NY,NX)     = 0.0_r8
          VLsoiAirPM(M+1,L,NY,NX)     = 0.0_r8
          ReductVLsoiAirPM(M,L,NY,NX) = VLsoiAirPM(M,L,NY,NX)
          THETPM(M+1,L,NY,NX)         = 0.0_r8
        ENDIF

      ENDDO D9785

      !
      !       RESET SURFACE LAYER NUMBER AND TRANSFER ALL WATER TO SOIL SURFACE LAYER
      !       IF LAKE SURFACE LAYER IS LOST TO EVAPORATION
      !
      !       NUM=new surface layer number after CO2CompenPoint_nodeete lake evaporation
      !       LakeSurfFlowMicP,LakeSurfFlowMacP,LakeSurfHeatFlux=lake surface water flux, heat flux if lake surface disappears
!
      IF(SoiBulkDensity_vr(NUM(NY,NX),NY,NX).LE.ZERO .AND. VHeatCapacity1_vr(NUM(NY,NX),NY,NX).LE.VHCPNX(NY,NX))THEN
        !the soil/water profile moves down        
        NUX=NUM(NY,NX)
        D9970: DO  LL=NUX+1,NL(NY,NX)
          IF(VLSoilPoreMicP_vr(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
            NUM(NY,NX)                     = LL
            H2OFlow2TopSoiMicP_col(NY,NX)  = WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX)
            H2OFlow2TopSoiMicPX_col(NY,NX) = WaterFlowSoiMicPX(3,NUM(NY,NX),NY,NX)
            H2OFlow2TopSoiMacP_col(NY,NX)  = WaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)
            HeatFlow2TopSoi_col(NY,NX)     = HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX)
            exit
          ENDIF
        ENDDO D9970
      ENDIF
    ENDDO D9790
  ENDDO D9795

  end subroutine UpdateStateFluxAtM

!------------------------------------------------------------------------------------------

  subroutine UpdateFluxAtExit(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
! begin_execution
  D9695: DO NX=NHW,NHE
    D9690: DO NY=NVN,NVS
      IF(NUM(NY,NX).EQ.NU(NY,NX))THEN
        !vertical flux, incoming
        !top layer is the same
        LakeSurfFlowMicP_col(NY,NX) = WaterFlowSoiMicP_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfFlowMicPX(NY,NX)    = WaterFlowSoiMicPX(3,N6X(NY,NX),NY,NX)
        LakeSurfFlowMacP_col(NY,NX) = WaterFlowMacP_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfHeatFlux_col(NY,NX) = HeatFlow2Soil_3D(3,N6X(NY,NX),NY,NX)
      ELSE
        !the top soil/water layer has changed
        LakeSurfFlowMicP_col(NY,NX) = H2OFlow2TopSoiMicP_col(NY,NX)
        LakeSurfFlowMicPX(NY,NX)    = H2OFlow2TopSoiMicPX_col(NY,NX)
        LakeSurfFlowMacP_col(NY,NX)     = H2OFlow2TopSoiMacP_col(NY,NX)
        LakeSurfHeatFlux_col(NY,NX) = HeatFlow2TopSoi_col(NY,NX)
      ENDIF
!      call writeSurfDiagnosis(I,J,NY,NX)
    ENDDO D9690
  ENDDO D9695
  
  end subroutine UpdateFluxAtExit

!--------------------------------------------------------------------------

  subroutine Config4TileDrainage(L,NY,NX,IFLGD,IFLGDH,DPTHH)

  implicit none
  integer , intent(in) :: L,NY,NX
  integer, intent(out) :: IFLGD,IFLGDH
  real(r8),intent(out) :: DPTHH

  integer :: LL
  real(r8) :: DTBLYX
!
!     IDWaterTable=water table flag
!     DPTH,DTBLY=depth to layer midpoint, artificial water table
!     PSISM1,PSISE=soil,air entry matric potential
!     DTBLYX=equilibrium water potential with artificial water table
!     IFLGD=micropore discharge flag to artificial water table
!
  IF(IDWaterTable(NY,NX).GE.3.AND.SoiDepthMidLay_vr(L,NY,NX).LT.DTBLY(NY,NX))THEN
    IF(PSISM1_vr(L,NY,NX).GT.mGravAccelerat*(SoiDepthMidLay_vr(L,NY,NX)-DTBLY(NY,NX)))THEN
      IFLGD=0
      IF(L.LT.NL(NY,NX))THEN
        D9568: DO  LL=L+1,NL(NY,NX)
          DTBLYX=DTBLY(NY,NX)+PSISE_vr(LL,NY,NX)/mGravAccelerat
          IF(SoiDepthMidLay_vr(LL,NY,NX).LT.DTBLYX)THEN
            IF((PSISM1_vr(LL,NY,NX).LE.mGravAccelerat*(SoiDepthMidLay_vr(LL,NY,NX)-DTBLYX) &
              .AND.L.NE.NL(NY,NX)).OR.SoiDepthMidLay_vr(LL,NY,NX).GT.ActiveLayDepth(NY,NX))THEN
              IFLGD=1
            ENDIF
          ENDIF
        ENDDO D9568
      ENDIF
    ELSE
      IFLGD=1
    ENDIF
  ELSE
    IFLGD=1
  ENDIF
!
!     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO TILE DRAIN
!
!     VLMacP1,VLWatMacP1,VLiceMacP1=macropore volume,water,ice content
!     CDPTH=depth to layer bottom
!     DLYR=layer thickness
!     IFLGDH=macropore discharge flag to artificial water table
!
  IF(VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    DPTHH=CumDepz2LayerBot_vr(L,NY,NX)-(VLWatMacP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX)) &
      /VLMacP1_vr(L,NY,NX)*DLYR(3,L,NY,NX)
  ELSE
    DPTHH=CumDepz2LayerBot_vr(L,NY,NX)
  ENDIF

  IF(IDWaterTable(NY,NX).GE.3 .AND. DPTHH.LT.DTBLY(NY,NX) .AND. VLWatMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
! artificial water table, e.g. tile drainage
    IFLGDH=0
    IF(L.LT.NL(NY,NX))THEN
      D9569: DO  LL=L+1,NL(NY,NX)
        IF(SoiDepthMidLay_vr(LL,NY,NX).LT.DTBLY(NY,NX))THEN
! the layer is above tile drain water table
          IF(VLMacP1_vr(LL,NY,NX).LE.ZEROS(NY,NX))THEN
! no macropore flow
            IFLGDH=1
          ENDIF
        ENDIF
      ENDDO D9569
    ENDIF
  ELSE
    IFLGDH=1
  ENDIF
  end subroutine Config4TileDrainage
!------------------------------------------------------------------

  subroutine Config4WaterTableDrain(L,NY,NX,IFLGU,IFLGUH,DPTHH)

  implicit none
  integer , intent(in) :: L,NY,NX
  integer , intent(out):: IFLGU,IFLGUH
  real(r8), intent(out):: DPTHH

  real(r8) :: ExtWaterTableEquil
  integer :: LL

!     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO WATER TABLE
!
!     IDWaterTable=water table flag
!     DPTH,ExtWaterTable=depth to layer midpoint,natural water table
!     PSISM1_vr(<0),PSISE_vr(<0)=matric,air entry water potential [MPa], from water height to MPa, 1.e-6*9.8*1.e3*H
!     ExtWaterTableEquil=equilibrium water potential with natural water table
!     ActiveLayDepth=active layer depth
!     IFLGU=micropore discharge flag to natural water table
!  total water potential psi+rho*grav*h

  IF(IDWaterTable(NY,NX).NE.0.AND.SoiDepthMidLay_vr(L,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN
    !the layer mid-depth is lower than water table
    IF(PSISM1_vr(L,NY,NX).GT.mGravAccelerat*(SoiDepthMidLay_vr(L,NY,NX)-ExtWaterTable_col(NY,NX)))THEN
      IFLGU=0
      D9565: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
        !water level difference
        ExtWaterTableEquil=ExtWaterTable_col(NY,NX)+PSISE_vr(LL,NY,NX)/mGravAccelerat  
        IF(SoiDepthMidLay_vr(LL,NY,NX).LT.ExtWaterTableEquil)THEN
          IF((PSISM1_vr(LL,NY,NX).LE.mGravAccelerat*(SoiDepthMidLay_vr(LL,NY,NX)-ExtWaterTableEquil) &
            .AND.L.NE.NL(NY,NX)).OR.SoiDepthMidLay_vr(LL,NY,NX).GT.ActiveLayDepth(NY,NX))THEN
            IFLGU=1
          ENDIF
        ENDIF
      ENDDO D9565
    ELSE
      IFLGU=1
    ENDIF
  ELSE
    IFLGU=1
  ENDIF
!
!     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO WATER TABLE
!
!     VLMacP1,VLWatMacP1,VLiceMacP1=macropore volume,water,ice content
!     DPTHH depth to layer macropore water
!     CDPTH=depth to layer bottom
!     DLYR=layer thickness
!     IFLGUH=macropore discharge flag to natural water table
!
  IF(VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    DPTHH=CumDepz2LayerBot_vr(L,NY,NX)-(VLWatMacP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX)) &
      /VLMacP1_vr(L,NY,NX)*DLYR(3,L,NY,NX)
  ELSE
    DPTHH=CumDepz2LayerBot_vr(L,NY,NX)
  ENDIF
  !
  IF(IDWaterTable(NY,NX).NE.0 .AND. DPTHH.LT.ExtWaterTable_col(NY,NX) &
    .AND. VLWatMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    !active water table
    IFLGUH=0
!     DO 9566 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
!     IF(SoiDepthMidLay_vr(LL,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN
!     IF(VLMacP1_vr(LL,NY,NX).LE.ZEROS(NY,NX))THEN
!     IFLGUH=1
!     ENDIF
!     ENDIFx
!9566  CONTINUE
  ELSE
    IFLGUH=1
  ENDIF
  end subroutine Config4WaterTableDrain

!------------------------------------------------------------------------------------------

  subroutine InitSoil3DModelMit(M,NY,NX)

  implicit none
  integer, intent(in) :: NY,NX,M
  integer :: L
  real(r8) :: THETWT,TScal4Aquadifsvity,THETWA
  real(r8) :: VLSoiPorAV,VLWatSoi,scalar
  real(r8) :: THETWH,Z3S

  D9885: DO L=NUM(NY,NX),NL(NY,NX)
    TMLiceThawMicP(L,NY,NX)           = 0.0_r8
    TMLiceThawMacP(L,NY,NX)           = 0.0_r8
    TLPhaseChangeHeat2Soi1(L,NY,NX)   = 0.0_r8
    TWatCharge2MicP(L,NY,NX)          = 0.0_r8
    TWatXChange2WatTableX(L,NY,NX)    = 0.0_r8
    TConvWaterFlowMacP_3D_vr(L,NY,NX) = 0.0_r8
    THeatFlow2Soili_vr(L,NY,NX)       = 0.0_r8
!
!   GAS EXCHANGE COEFFICIENTS SOIL LAYERS
!
!   VOLA1,VLiceMicP1,VLWatMicP1=total,ice-,water-filled microporosity
!   VLMacP1,VLiceMacP1,VLWatMacP1=total,ice-,water-filled macroporosity
!   VLsoiAirPM=air-filled porosity
!   TScal4Aquadifsvity=temperature effect on gas diffusivity
!   DiffusivitySolutEff=rate constant for air-water gas exchange
!   Z1S,Z2SW,Z2SD,Z3SX=parameters for soil air-water gas transfers
!   XNPD=time step for gas transfer calculations
!   TORT,TortMacPM=tortuosity for aqueous diffn in micropores,macropres
!
    VLWatSoi   = VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)
    VLSoiPorAV = VLMicP1_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
    IF(VLSoiPorAV.GT.ZEROS2(NY,NX) .AND. VLsoiAirPM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWA                         = AZMAX1(AMIN1(1.0_r8,VLWatSoi/VLSoiPorAV))
      TScal4Aquadifsvity             = TEFAQUDIF(TKSoi1_vr(0,NY,NX))
      Z3S                            = FieldCapacity_vr(L,NY,NX)/POROS_vr(L,NY,NX)
      scalar                         = TScal4Aquadifsvity*XNPD
      DiffusivitySolutEff(M,L,NY,NX) = fDiffusivitySolutEff(scalar,THETWA,Z3S)
    ELSE
      DiffusivitySolutEff(M,L,NY,NX)=0.0_r8
    ENDIF

    IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      THETWT=safe_adb(VLWatMicPM_vr(M,L,NY,NX),VLSoilMicP_vr(L,NY,NX))
      TortMicPM_vr(M,L,NY,NX)=TortMicporew(THETWT)*(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
    ELSE
!   standing water has tortuosity 0.7?
      TortMicPM_vr(M,L,NY,NX)=0.7_r8
    ENDIF
    IF(VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWH=VLWatMacPM(M,L,NY,NX)/VLMacP1_vr(L,NY,NX)
      TortMacPM(M,L,NY,NX)=TortMacporew(THETWH)*SoilFracAsMacP_vr(L,NY,NX)
    ELSE
      TortMacPM(M,L,NY,NX)=0.0_r8
    ENDIF
  ENDDO D9885
  end subroutine InitSoil3DModelMit


!------------------------------------------------------------------------------------------
  subroutine MicporeDarcyFlow(NY,NX,N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,KSatRedusByRainKinetEnergy,&
    HeatByWatFlowMicP,PSISV1,PSISVL)          
  implicit none
  integer, intent(in)  :: NY,NX,N
  integer, intent(in)  :: N1,N2,N3  !source grid
  integer, intent(in)  :: N4,N5,N6  !destination grid
  real(r8), intent(in) :: THETA1,THETAL,KSatRedusByRainKinetEnergy
  real(r8), intent(out):: HeatByWatFlowMicP,PSISV1,PSISVL
  real(r8) :: AVE_CONDUCTANCE,HydCondSrc,HydCondDest
  real(r8) :: PredDarcyFlowMax,PSISTL,PSIST1,THETW1
  real(r8) :: WatDarcyFlowMicP,HeatByDarcyFlowMicP,PtWatDarcyFlux,PredDarcyFlow,PSISTMP,THETWL
  integer :: K1,KL


  !     THETW1,THETWL=water concentrations in source,destination cells
  !air entry pressure is the pressure when large-pores started to be air-filled, when soil matric pressure is less
  !than air entry pressure, it is unsaturated
  IF(PSISoilMatricPtmp_vr(N3,N2,N1).GT.PSISoilAirEntry(N3,N2,N1).AND. &
    PSISoilMatricPtmp_vr(N6,N5,N4).GT.PSISoilAirEntry(N6,N5,N4))THEN
    !both are saturated
    THETW1              = THETA1
    THETWL              = THETAL
    K1                  = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-THETW1)/POROS_vr(N3,N2,N1))+1))
    KL                  = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N6,N5,N4)-THETWL)/POROS_vr(N6,N5,N4))+1))
    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    PSISM1_vr(N6,N5,N4) = PSISoilMatricPtmp_vr(N6,N5,N4)
    !
    !     GREEN-AMPT FLOW IF ONE LAYER IS SATURATED
    !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POENTIAL)
    !
    !     GREEN-AMPT FLOW IF SOURCE CELL SATURATED
    !ThetaSat_vr=micropore soil water content
  ELSEIF(PSISoilMatricPtmp_vr(N3,N2,N1).GT.PSISoilAirEntry(N3,N2,N1))THEN
    !source grid is saturated
    THETW1              = THETA1
    THETWL              = AMAX1(THETY_vr(N6,N5,N4),AMIN1(POROS_vr(N6,N5,N4),&
      safe_adb(VLWatMicPX1_vr(N6,N5,N4),VLSoilMicP_vr(N6,N5,N4))))
    K1                  = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-THETW1)/POROS_vr(N3,N2,N1))+1))
    KL                  = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N6,N5,N4) &
      -AMIN1(ThetaSat_vr(N6,N5,N4),THETWL))/POROS_vr(N6,N5,N4))+1))
    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    
    IF(VLSoilMicPMass_vr(N6,N5,N4).GT.ZEROS(NY,NX))THEN
      !dest layer is pore active
      IF(THETWL.LT.FieldCapacity_vr(N6,N5,N4))THEN
        !less than field capacity
        PSISM1_vr(N6,N5,N4)=AMAX1(PSIHY,-EXP(LOGPSIFLD(N5,N4)+((LOGFldCapacity_vr(N6,N5,N4)-LOG(THETWL)) &
          /FCD(N6,N5,N4)*LOGPSIMND(N5,N4))))
      ELSEIF(THETWL.LT.POROS_vr(N6,N5,N4)-DTHETW)THEN
        PSISM1_vr(N6,N5,N4)=-EXP(LOGPSIAtSat(N5,N4)+(((LOGPOROS_vr(N6,N5,N4)-LOG(THETWL)) &
          /PSD(N6,N5,N4))**SRP(N6,N5,N4)*LOGPSIMXD(N5,N4)))
      ELSE
        !saturated
        THETWL=POROS_vr(N6,N5,N4)
        PSISM1_vr(N6,N5,N4)=PSISE_vr(N6,N5,N4)
      ENDIF
    ELSE
      THETWL=POROS_vr(N6,N5,N4)
      PSISM1_vr(N6,N5,N4)=PSISE_vr(N6,N5,N4)
    ENDIF
    !
    !     GREEN-AMPT FLOW IF ADJACENT CELL SATURATED
!
  ELSEIF(PSISoilMatricPtmp_vr(N6,N5,N4).GT.PSISoilAirEntry(N6,N5,N4))THEN
    !source grid is saturated, dest grid is not
    THETW1 = AMAX1(THETY_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicPX1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
    THETWL = THETAL
    K1     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-AMIN1(ThetaSat_vr(N3,N2,N1),THETW1))/POROS_vr(N3,N2,N1))+1))
    KL     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N6,N5,N4)-THETWL)/POROS_vr(N6,N5,N4))+1))

    IF(VLSoilMicPMass_vr(N3,N2,N1).GT.ZEROS(NY,NX))THEN
      IF(THETW1.LT.FieldCapacity_vr(N3,N2,N1))THEN
        PSISTMP             = -EXP(LOGPSIFLD(N2,N1)+(LOGFldCapacity_vr(N3,N2,N1)-LOG(THETW1)) &
          /FCD(N3,N2,N1)*LOGPSIMND(N2,N1))
        PSISM1_vr(N3,N2,N1) = AMAX1(PSIHY,PSISTMP)
      ELSEIF(THETW1.LT.POROS_vr(N3,N2,N1)-DTHETW)THEN
        PSISM1_vr(N3,N2,N1)=-EXP(LOGPSIAtSat(N2,N1)+(((LOGPOROS_vr(N3,N2,N1)-LOG(THETW1)) &
          /PSD(N3,N2,N1))**SRP(N3,N2,N1)*LOGPSIMXD(N2,N1)))
      ELSE
        THETW1              = POROS_vr(N3,N2,N1)
        PSISM1_vr(N3,N2,N1) = PSISE_vr(N3,N2,N1)
      ENDIF
    ELSE
      THETW1=POROS_vr(N3,N2,N1)
      PSISM1_vr(N3,N2,N1)=PSISE_vr(N3,N2,N1)
    ENDIF
    !
    !     RICHARDS FLOW IF NEITHER CELL IS SATURATED
    !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POTENTIAL)
    !
  ELSE
    THETW1              = THETA1
    THETWL              = THETAL
    K1                  = MAX(1,MIN(100,INT(100.0*(POROS_vr(N3,N2,N1)-THETW1)/POROS_vr(N3,N2,N1))+1))
    KL                  = MAX(1,MIN(100,INT(100.0*(POROS_vr(N6,N5,N4)-THETWL)/POROS_vr(N6,N5,N4))+1))
    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    PSISM1_vr(N6,N5,N4) = PSISoilMatricPtmp_vr(N6,N5,N4)
  ENDIF

  !
  !     HYDRAULIC CONUCTIVITY
  !
  !     HydCondSrc,HydCondDest=hydraulic conductivity of source,destination layer
  !     HydroCond_3D=lateral(1,2),vertical(3) micropore hydraulic conductivity
  !
  IF(N3.EQ.NUM(NY,NX))THEN
    !surface soil
    HydCondSrc=HydroCond_3D(N,K1,N3,N2,N1)*KSatRedusByRainKinetEnergy
  ELSE
    HydCondSrc=HydroCond_3D(N,K1,N3,N2,N1)
  ENDIF
  HydCondDest=HydroCond_3D(N,KL,N6,N5,N4)
  !
  !     TOTAL SOIL WATER POTENTIAL = MATRIC, GRAVIMETRIC + OSMOTIC
  !
  !     PSISM1,PSISH,PSISO=soil matric,gravitational,osmotic potentials
  !src
  PSIST1=PSISM1_vr(N3,N2,N1)+PSIGrav_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  PSISV1=PSISM1_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  !dest
  PSISTL=PSISM1_vr(N6,N5,N4)+PSIGrav_vr(N6,N5,N4)+PSISoilOsmotic_vr(N6,N5,N4)
  PSISVL=PSISM1_vr(N6,N5,N4)+PSISoilOsmotic_vr(N6,N5,N4)
  !
  !     HYDRAULIC CONDUCTIVITY FROM CURRENT WATER CONTENT
  !     AND LOOKUP ARRAY GENERATED IN 'HOUR1'
  !
  !     HydCondSrc,HydCondDest=hydraulic conductivities in source,destination cells
  !     KSatRedusByRainKinetEnergy=reduction in soil surface Ksat from rainfall energy impact
  !     AVE_CONDUCTANCE=source-destination hydraulic conductance
  !     DLYR=layer thickness
  !
  IF(HydCondSrc.GT.ZERO.AND.HydCondDest.GT.ZERO)THEN
    AVE_CONDUCTANCE=2.0_r8*HydCondSrc*HydCondDest/(HydCondSrc*DLYR(N,N6,N5,N4)+HydCondDest*DLYR(N,N3,N2,N1))
  ELSE
    AVE_CONDUCTANCE=0.0_r8
  ENDIF
  !
  !     WATER FLUX FROM WATER POTENTIALS, HYDRAULIC CONDUCTIVITY
  !     CONSTRAINED BY WATER POTENTIAL GRADIENT, COUPLED WITH
  !     CONVECTIVE HEAT FLUX FROM WATER FLUX
  !
  !     PtWatDarcyFlux,WatDarcyFlowMicP=micropore water flux unlimited,limited by source water
  !     dts_HeatWatTP=time step of flux calculations
  !     VOLW2,VOLP1=water,air contents of source,destination micropores
  !     HeatByWatFlowMicP=convective heat flux from micropore water flux
  !     VLairMicP=excess water+ice relative to porosity
  !
  PtWatDarcyFlux=AVE_CONDUCTANCE*(PSIST1-PSISTL)*AREA(N,N3,N2,N1)*dts_HeatWatTP

  IF(PtWatDarcyFlux.GE.0.0_r8)THEN
    !flow from src to dest
    IF(THETW1.GT.ThetaSat_vr(N3,N2,N1))THEN
      !water more than air-entry saturation
      PredDarcyFlow=PtWatDarcyFlux+AMIN1((THETW1-ThetaSat_vr(N3,N2,N1))*VLSoilMicP_vr(N3,N2,N1),&
        AZMAX1((ThetaSat_vr(N6,N5,N4)-THETWL)*VLSoilMicP_vr(N6,N5,N4)))*dts_wat
    ELSE
      PredDarcyFlow=PtWatDarcyFlux
    ENDIF

    WatDarcyFlowMicP=AZMAX1(AMIN1(PredDarcyFlow,VLWatMicP2_vr(N3,N2,N1)*dts_wat,VLairMicP1_vr(N6,N5,N4)*dts_wat))
    PredDarcyFlowMax=AZMAX1(AMIN1(PtWatDarcyFlux,VLWatMicP2_vr(N3,N2,N1)*dts_wat,VLairMicP1_vr(N6,N5,N4)*dts_wat))
    !     FLQL1=(THETW1-ThetaSat_vr(N3,N2,N1))*VLSoilMicP_vr(N3,N2,N1)
    !     FLQL2=(ThetaSat_vr(N6,N5,N4)-THETWL)*VLSoilMicP_vr(N6,N5,N4)
    !     FLQL3=PtWatDarcyFlux+AMIN1(FLQL1,AZMAX1(FLQL2))*dts_wat
    !     FLQL4=AZMAX1(AMIN1(FLQL3,VLairMicP1_vr(N6,N5,N4)*dts_wat))
  ELSE
    !flow from dest to src
    IF(THETWL.GT.ThetaSat_vr(N6,N5,N4))THEN
      PredDarcyFlow=PtWatDarcyFlux+AMAX1((ThetaSat_vr(N6,N5,N4)-THETWL)*VLSoilMicP_vr(N6,N5,N4), &
        AZMIN1((THETW1-ThetaSat_vr(N3,N2,N1))*VLSoilMicP_vr(N3,N2,N1)))*dts_wat
    ELSE
      PredDarcyFlow=PtWatDarcyFlux
    ENDIF

    WatDarcyFlowMicP=AZMIN1(AMAX1(PredDarcyFlow,-VLWatMicP2_vr(N6,N5,N4)*dts_wat,-VLairMicP1_vr(N3,N2,N1)*dts_wat))
    PredDarcyFlowMax=AZMIN1(AMAX1(PtWatDarcyFlux,-VLWatMicP2_vr(N6,N5,N4)*dts_wat,-VLairMicP1_vr(N3,N2,N1)*dts_wat))
    !     FLQL1=(ThetaSat_vr(N6,N5,N4)-THETWL)*VLSoilMicP_vr(N6,N5,N4)
    !     FLQL2=(THETW1-ThetaSat_vr(N3,N2,N1))*VLSoilMicP_vr(N3,N2,N1)
    !     FLQL3=PtWatDarcyFlux+AMAX1(FLQL1,AZMIN1(FLQL2))*dts_wat
    !     FLQL4=AZMIN1(AMAX1(FLQL3,-VLairMicP1_vr(N3,N2,N1)*dts_wat))
  ENDIF

  IF(N.EQ.iVerticalDirection .AND. VLairMicP_vr(N6,N5,N4).LT.0.0_r8)THEN
    !flow in the vertical direction, destination grid is saturated
    WatDarcyFlowMicP = WatDarcyFlowMicP+AZMIN1(AMAX1(-VLWatMicP2_vr(N6,N5,N4)*dts_wat,VLairMicP_vr(N6,N5,N4)))
    PredDarcyFlowMax = PredDarcyFlowMax+AZMIN1(AMAX1(-VLWatMicP2_vr(N6,N5,N4)*dts_wat,VLairMicP_vr(N6,N5,N4)))
  ENDIF

  IF(WatDarcyFlowMicP.GT.0.0_r8)THEN
    HeatByDarcyFlowMicP=cpw*TKSoi1_vr(N3,N2,N1)*WatDarcyFlowMicP
  ELSE
    HeatByDarcyFlowMicP=cpw*TKSoi1_vr(N6,N5,N4)*WatDarcyFlowMicP
  ENDIF

  VLWatMicP2_vr(N3,N2,N1)          = VLWatMicP2_vr(N3,N2,N1)-WatDarcyFlowMicP
  VLWatMicP2_vr(N6,N5,N4)          = VLWatMicP2_vr(N6,N5,N4)+WatDarcyFlowMicP
  WatXChange2WatTable(N,N6,N5,N4)  = WatDarcyFlowMicP
  WatXChange2WatTableX(N,N6,N5,N4) = PredDarcyFlowMax
  HeatByWatFlowMicP                = HeatByDarcyFlowMicP
  end subroutine MicporeDarcyFlow
!------------------------------------------------------------------------------------------
  subroutine MacporeFLow(NY,NX,M,N,N1,N2,N3,N4,N5,N6,ConvectiveHeatFlxMacP,LInvalidMacP)
  implicit none
  integer, intent(in) :: NY,NX,M,N
  integer, intent(in) :: N1,N2,N3   !source grid
  integer, intent(in) :: N4,N5,N6   !destination grid
  real(r8), intent(out) :: ConvectiveHeatFlxMacP   ! >0, flux into dst, 
  logical, intent(inout) :: LInvalidMacP   !==1, dest grid invalid
  real(r8) :: FLWHX,PSISH1,PSISHL
  !     PSISH1,PSISHL=macropore total water potl in source,destination
  !     DLYR=layer thickness
  !     VLWatMacP1,VOLPH1=macropore water,air content

  IF(VLMacP1_vr(N3,N2,N1).GT.ZEROS2(N2,N1).AND.VLMacP1_vr(N6,N5,N4).GT.ZEROS2(N5,N4).AND.(.not.LInvalidMacP))THEN
    PSISH1=PSIGrav_vr(N3,N2,N1)+mGravAccelerat*DLYR(3,N3,N2,N1)*(AMIN1(1.0_r8,&
      AZMAX1(VLWatMacP1_vr(N3,N2,N1)/VLMacP1_vr(N3,N2,N1)))-0.5_r8)
    PSISHL=PSIGrav_vr(N6,N5,N4)+mGravAccelerat*DLYR(3,N6,N5,N4)*(AMIN1(1.0_r8,&
      AZMAX1(VLWatMacP1_vr(N6,N5,N4)/VLMacP1_vr(N6,N5,N4)))-0.5_r8)
    !
    !     MACROPORE FLOW IF GRAVITATIONAL GRADIENT IS POSITIVE
    !     AND MACROPORE POROSITY EXISTS IN ADJACENT CELL
    !
    !     FLWHX,ConvWaterFlowMacP_3D=macropore water flux unltd,ltd by source water
    !     dts_HeatWatTP=time step of flux calculations
    !     VOLW2,VOLP1=water,air contents of source,destination micropores
    !     ConvectiveHeatFlxMacP=convective heat flux from micropore water flux
    !
    FLWHX=AVCNHL(N,N6,N5,N4)*(PSISH1-PSISHL)*AREA(N,N3,N2,N1)*dts_HeatWatTP
    IF(N.NE.iVerticalDirection)THEN
      !horizontal direction
      IF(PSISH1.GT.PSISHL)THEN
        !flow from source to dest grid
        ConvWaterFlowMacP_3D(N,N6,N5,N4)=AZMAX1(AMIN1(AMIN1(VLWatMacP1_vr(N3,N2,N1),&
          VLairMacP1_vr(N6,N5,N4))*dts_wat,FLWHX))
      ELSEIF(PSISH1.LT.PSISHL)THEN
        !flow from dest grid to source
        ConvWaterFlowMacP_3D(N,N6,N5,N4)=AZMIN1(AMAX1(AMAX1(-VLWatMacP1_vr(N6,N5,N4),&
          -VLairMacP1_vr(N3,N2,N1))*dts_wat,FLWHX))
      ELSE
        ConvWaterFlowMacP_3D(N,N6,N5,N4)=0.0_r8
      ENDIF
    ELSE
      !vertical direction
      ConvWaterFlowMacP_3D(N,N6,N5,N4)=AZMAX1(AMIN1(AMIN1(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
        +ConvWaterFlowMacP_3D(N,N3,N2,N1),VLairMacP1_vr(N6,N5,N4)*dts_wat),FLWHX))
    ENDIF

    IF(N.EQ.iVerticalDirection)THEN
      ConvWaterFlowMacP_3D(N,N6,N5,N4)=ConvWaterFlowMacP_3D(N,N6,N5,N4)+AZMIN1(VLairMacP_vr(N6,N5,N4))
    ENDIF
    WaterFlow2MacPM(M,N,N6,N5,N4)=ConvWaterFlowMacP_3D(N,N6,N5,N4)
  ELSE
    ConvWaterFlowMacP_3D(N,N6,N5,N4)=0.0_r8
    WaterFlow2MacPM(M,N,N6,N5,N4)=0.0_r8
    IF(VLairMacP1_vr(N6,N5,N4).LE.0.0_r8)LInvalidMacP=.true.
  ENDIF

  IF(ConvWaterFlowMacP_3D(N,N6,N5,N4).GT.0.0_r8)THEN
    ConvectiveHeatFlxMacP=cpw*TKSoi1_vr(N3,N2,N1)*ConvWaterFlowMacP_3D(N,N6,N5,N4)
  ELSE
    ConvectiveHeatFlxMacP=cpw*TKSoi1_vr(N6,N5,N4)*ConvWaterFlowMacP_3D(N,N6,N5,N4)
  ENDIF
  end subroutine MacporeFLow  
!------------------------------------------------------------------------------------------

  subroutine WaterVaporFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectVapFlux,ConvectHeatFluxMicP)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), intent(in) :: PSISV1,PSISVL
  real(r8), intent(out) :: ConvectVapFlux,ConvectHeatFluxMicP
  real(r8) :: TK11,TK12,VP1,VPL,VPY,CNV1,CNVL
  REAL(R8) :: ATCNVL,PotentialVaporFlux,MaxVaporFlux
  !     VAPOR PRESSURE AND DIFFUSIVITY IN EACH GRID CELL
  !
  
  !     THETPM,THETX=current, minimum air-filled porosity
  !     TK11,TK12=interim soil temperature in source,destination
  !     VP1,VPL=vapor concentration in source,destination
  !     PSISV1,PSISVL=matric+osmotic water potl in source,destination
  !     CNV1,CNVL=vapor conductivities of source, destination
  !     POROS,POROQ=porosity, tortuosity
  !     WVapDifusvitySoil_vr=vapor diffusivity
  !     ATCNVL=source,destination vapor conductance
  !     DLYR=soil layer depth
  !     PotentialVaporFlux,MaxVaporFlux=vapor flux unlimited,limited by vapor
  !     VPY=equilibrium vapor concentration
  !     dts_wat=time step for flux calculations
  !     ConvectVapFlux,ConvectiveHeatFlux=vapor flux and its convective heat flux
!
  IF(THETPM(M,N3,N2,N1).GT.THETX.AND.THETPM(M,N6,N5,N4).GT.THETX)THEN
    TK11 = TKSoi1_vr(N3,N2,N1)
    TK12 = TKSoi1_vr(N6,N5,N4)

    VP1    = vapsat(TK11)*EXP(18.0_r8*PSISV1/(RGASC*TK11))
    VPL    = vapsat(TK12)*EXP(18.0_r8*PSISVL/(RGASC*TK12))
    CNV1   = WVapDifusvitySoil_vr(N3,N2,N1)*THETPM(M,N3,N2,N1)*POROQ*THETPM(M,N3,N2,N1)/POROS_vr(N3,N2,N1)
    CNVL   = WVapDifusvitySoil_vr(N6,N5,N4)*THETPM(M,N6,N5,N4)*POROQ*THETPM(M,N6,N5,N4)/POROS_vr(N6,N5,N4)
    ATCNVL = 2.0_r8*CNV1*CNVL/(CNV1*DLYR(N,N6,N5,N4)+CNVL*DLYR(N,N3,N2,N1))
    !
    !     VAPOR FLUX FROM VAPOR PRESSURE AND DIFFUSIVITY,
    !     AND CONVECTIVE HEAT FLUX FROM VAPOR FLUX
!
    PotentialVaporFlux = ATCNVL*(VP1-VPL)*AREA(N,N3,N2,N1)*dts_HeatWatTP
    VPY                = (VP1*VLsoiAirPM(M,N3,N2,N1)+VPL*VLsoiAirPM(M,N6,N5,N4))/(VLsoiAirPM(M,N3,N2,N1)&
      +VLsoiAirPM(M,N6,N5,N4))
    MaxVaporFlux = (VP1-VPY)*VLsoiAirPM(M,N3,N2,N1)*dts_wat
    !out from source grid 
    IF(PotentialVaporFlux.GE.0.0_r8)THEN
      ConvectVapFlux      = AZMAX1(AMIN1(PotentialVaporFlux,MaxVaporFlux))
      ConvectHeatFluxMicP = (cpw*TKSoi1_vr(N3,N2,N1)+EvapLHTC)*ConvectVapFlux
    !out from dest grid  
    ELSE
      ConvectVapFlux      = AZMIN1(AMAX1(PotentialVaporFlux,MaxVaporFlux))
      ConvectHeatFluxMicP = (cpw*TKSoi1_vr(N6,N5,N4)+EvapLHTC)*ConvectVapFlux
    ENDIF
  ELSE
    ConvectVapFlux      = 0.0_r8
    ConvectHeatFluxMicP = 0.0_r8
  ENDIF
  end subroutine WaterVaporFlow  
!------------------------------------------------------------------------------------------

  subroutine Solve4Heat(I,J,N,NY,NX,N1,N2,N3,N4,N5,N6,ConvectHeatFluxMicP,HeatFluxAir2Soi)
  implicit none
  integer , intent(in) :: I,J
  integer , intent(in) :: N,NY,NX
  integer , intent(in) :: N1,N2,N3  !source
  integer , intent(in) :: N4,N5,N6  !dest
  real(r8), intent(in) :: ConvectHeatFluxMicP,HeatFluxAir2Soi
  real(r8) :: TK1X,TKLX,TKY,HFLWC,HFLWX,HeatCondSoi
  real(r8) :: ATCNDL,DTKX
  real(r8) :: ThermCondDst,ThermCondSrc
  !

  !     THERMAL CONDUCTIVITY IN EACH GRID CELL
  !
  !     DTH*,RYL*,DNU*,TRB*,XNUS*=turbulence effects on thermal conductivity
  !     THETWX,FracSoiPAsAir=water,air concentration
  !     TCNDW*,TCNDA*=thermal conductivity of water,air
  !     ThermCondSrc,TCNDL=soil thermal conductivity in source,destination
  !     WTHET*=multiplier for air concn in thermal conductivity
  !     ATCNDL=source-destination thermal conductance
  !
  DTKX=ABS(TKSoi1_vr(N3,N2,N1)-TKSoi1_vr(N6,N5,N4))*ppmc

  call CalcSoilThermConductivity(N1,N2,N3,DTKX,ThermCondSrc)

  call CalcSoilThermConductivity(N4,N5,N6,DTKX,ThermCondDst)

  ATCNDL=safe_adb(2.0_r8*ThermCondSrc*ThermCondDst,(ThermCondSrc*DLYR(N,N6,N5,N4)+ThermCondDst*DLYR(N,N3,N2,N1)))

  !     HEAT FLOW FROM THERMAL CONDUCTIVITY AND TEMPERATURE GRADIENT
  !
  !     VHeatCapacity1_vr,VHCPW=volumetric heat capacity of soil,snowpack
  !     TK1X,TKLX=interim temperatures of source,destination
  !     ConvectiveHeatFlux,HeatFluxAir2Soi=convective heat from soil vapor flux
  !     HeatFluxAir2Soi=storage heat flux from snowpack
  !     TKY=equilibrium source-destination temperature
  !     HFLWC,HFLWX=source-destination heat flux unltd,ltd by heat
  !     ATCNDL=source-destination thermal conductance
  !     HeatCondSoi=source-destination conductive heat flux
  !     HFLWL=total conductive+convective source-destination heat flux
  !
  
  IF(VHeatCapacity1_vr(N3,N2,N1).GT.VHCPNX(NY,NX))THEN
    IF(N3.EQ.NUM(NY,NX) .AND. VLHeatCapSnow_snvr(1,N2,N1).LE.VLHeatCapSnowMin_col(N2,N1))THEN
      !surface layer, not significant snowpack
      TK1X=TKSoi1_vr(N3,N2,N1)-(ConvectHeatFluxMicP-HeatFluxAir2Soi)/VHeatCapacity1_vr(N3,N2,N1)
      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoi1_vr(N3,N2,N1)-ConvectHeatFluxMicP/VHeatCapacity1_vr(N3,N2,N1)',&
          TKSoi1_vr(N3,N2,N1),ConvectHeatFluxMicP,HeatFluxAir2Soi,VHeatCapacity1_vr(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    ELSE
      TK1X=TKSoi1_vr(N3,N2,N1)-ConvectHeatFluxMicP/VHeatCapacity1_vr(N3,N2,N1)
      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoi1_vr(N3,N2,N1)-ConvectHeatFluxMicP/VHeatCapacity1_vr(N3,N2,N1)',&
          TKSoi1_vr(N3,N2,N1),ConvectHeatFluxMicP,VHeatCapacity1_vr(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    ENDIF
  ELSE
    TK1X=TKSoi1_vr(N3,N2,N1)
  ENDIF

  !if(I>=132 .and. J==19)print*,'tklx'
  IF(VHeatCapacity1_vr(N6,N5,N4).GT.ZEROS(NY,NX))THEN
    TKLX=TKSoi1_vr(N6,N5,N4)+ConvectHeatFluxMicP/VHeatCapacity1_vr(N6,N5,N4)
  ELSE
    TKLX=TKSoi1_vr(N6,N5,N4)
  ENDIF
  
  if(VHeatCapacity1_vr(N3,N2,N1)+VHeatCapacity1_vr(N6,N5,N4)>0._r8)then
    !temperature at the interface between (N3,N2,N1) and (N6,N5,N4)
    TKY=(VHeatCapacity1_vr(N3,N2,N1)*TK1X+VHeatCapacity1_vr(N6,N5,N4)*TKLX)/(VHeatCapacity1_vr(N3,N2,N1)&
      +VHeatCapacity1_vr(N6,N5,N4))
  ELSE
    TKY=(TK1X+TKLX)/2._r8
  endif 

  !if(I>=132 .and. J==19)print*,'aftky'
  HFLWX=(TK1X-TKY)*VHeatCapacity1_vr(N3,N2,N1)*dts_wat
  HFLWC=ATCNDL*(TK1X-TKLX)*AREA(N,N3,N2,N1)*dts_HeatWatTP
  IF(HFLWC.GE.0.0_r8)THEN
    !from (N3,N2,N1) into (N6,N5,N4)
    HeatCondSoi=AZMAX1(AMIN1(HFLWX,HFLWC))
  ELSE
    !from (N6,N5,N4) into (N3,N2,N1)
    HeatCondSoi=AZMIN1(AMAX1(HFLWX,HFLWC))
  ENDIF
  HeatFlow2Soili(N,N6,N5,N4)=HeatFlow2Soili(N,N6,N5,N4)+HeatCondSoi

  end subroutine Solve4Heat  
!------------------------------------------------------------------------------------------
  subroutine WaterTBLDrain(N,N1,N2,N3,M4,M5,M6,IFLGU,IFLGUH,RechargSubSurf,RechargRateWTBL,DPTHH,XN)
  implicit none
  integer, intent(in) :: N,N1,N2,N3,M4,M5,M6
  integer, intent(in) :: IFLGU,IFLGUH
  real(r8), intent(in):: RechargSubSurf,RechargRateWTBL,DPTHH,XN
  real(r8) :: FLWT,MacPoreDischarg,FLWTH,PSISWD,PSISWT,PSISWTH

!     MICROPORE DISCHARGE ABOVE WATER TABLE
!
!     IFLGU=micropore discharge flag to natural water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     WaterTBLSlope=water table slope
!     PSISWT=water potential driving micropore discharge
!     PSISoilMatricPtmp_vr,PSISO=matric,osmotic water potential
!     DPTH,ExtWaterTable=depth to layer midpoint,natural water table
!     DepthInternalWTBL=depth to internal water table
!     FLWL=micropore discharge to natural water table
!     HFLWL=convective heat from discharge to natural water table
!     HydroCond_3D=saturated hydraulic conductivity
!     AREAU=fraction of layer below natural water table
!
  IF(IFLGU.EQ.0 .AND. (.not.isclose(RechargRateWTBL,0._r8)))THEN
    PSISWD = XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope(N2,N1))
    PSISWT = AZMIN1(-PSISoilMatricPtmp_vr(N3,N2,N1)-0.03_r8*PSISoilOsmotic_vr(N3,N2,N1) &
      +mGravAccelerat*(SoiDepthMidLay_vr(N3,N2,N1)-ExtWaterTable_col(N2,N1)) &
      -mGravAccelerat*AZMAX1(SoiDepthMidLay_vr(N3,N2,N1)-DepthInternalWTBL(N2,N1)))
    IF(PSISWT.LT.0.0_r8)PSISWT=PSISWT-PSISWD
    FLWT=PSISWT*HydroCond_3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*(1.0_r8-AREAU(N3,N2,N1))/(RechargSubSurf+1.0_r8) &
      *RechargRateWTBL*dts_HeatWatTP
    WatXChange2WatTable(N,M6,M5,M4)  = XN*FLWT
    WatXChange2WatTableX(N,M6,M5,M4) = XN*FLWT
    HeatFlow2Soili(N,M6,M5,M4)       = cpw*TKSoi1_vr(N3,N2,N1)*XN*FLWT
!    if(M6==1)then
!      write(*,*)'3HeatFlow2Soili(N,M6,M5,M4)=',HeatFlow2Soili(N,M6,M5,M4)
!    endif
  ELSE
    WatXChange2WatTable(N,M6,M5,M4)  = 0.0_r8
    WatXChange2WatTableX(N,M6,M5,M4) = 0.0_r8
    HeatFlow2Soili(N,M6,M5,M4)       = 0.0_r8
  ENDIF
!
!     MACROPORE DISCHARGE ABOVE WATER TABLE
!
!     IFLGUH=macropore discharge flag to natural water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     WaterTBLSlope=water table slope
!     PSISWTH=water potential driving macropore discharge
!     PSISO=osmotic water potential
!     DPTHH,ExtWaterTable=depth to layer macropore water,natural water table
!     DepthInternalWTBL=depth to internal water table
!     FLWTH,MacPoreDischarg=macropore discharge unltd,ltd by macropore water
!     HydroCondMacP1=macropore hydraulic conductivity
!     ConvWaterFlowMacP_3D=macropore discharge to natural water table
!     HFLWL=convective heat from discharge to natural water table
!     HydroCond_3D=saturated hydraulic conductivity
  !     AREAU=fraction of layer below natural water table
!
  IF(IFLGUH.EQ.0 .AND. (.not.isclose(RechargRateWTBL,0._r8)).AND.VLMacP1_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    PSISWD  = XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope(N2,N1))
    PSISWTH = -0.03_r8*PSISoilOsmotic_vr(N3,N2,N1)+mGravAccelerat*(DPTHH-ExtWaterTable_col(N2,N1)) &
      -mGravAccelerat*AZMAX1(DPTHH-DepthInternalWTBL(N2,N1))
    IF(PSISWTH.LT.0.0_r8)PSISWTH=PSISWTH-PSISWD
    FLWTH=PSISWTH*HydroCondMacP1_vr(N3,N2,N1)*AREA(N,N3,N2,N1) &
      *(1.0_r8-AREAU(N3,N2,N1))/(RechargSubSurf+1.0)*RechargRateWTBL*dts_HeatWatTP

    MacPoreDischarg=AMAX1(FLWTH,AZMIN1(-(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
      +ConvWaterFlowMacP_3D(3,N3,N2,N1)-ConvWaterFlowMacP_3D(3,N3+1,N2,N1))))
    ConvWaterFlowMacP_3D(N,M6,M5,M4)=XN*MacPoreDischarg
    HeatFlow2Soili(N,M6,M5,M4)=HeatFlow2Soili(N,M6,M5,M4)+cpw*TKSoi1_vr(N3,N2,N1)*XN*MacPoreDischarg
  ELSE
    ConvWaterFlowMacP_3D(N,M6,M5,M4)=0.0_r8
  ENDIF
  end subroutine WaterTBLDrain  
!------------------------------------------------------------------------------------------
  subroutine TileDrain(N,N1,N2,N3,M4,M5,M6,IFLGD,IFLGDH,RechargRateWTBL,RechargSubSurf,DPTHH,XN)
  implicit none
  integer, intent(in) :: N,N1,N2,N3,M4,M5,M6,IFLGD,IFLGDH
  real(r8),intent(in) :: RechargRateWTBL,RechargSubSurf,DPTHH,XN
  real(r8) :: FLWT,FLWTH,MacPoreDischarg,PSISWD,PSISWT,PSISWTH
!
  !     MICROPORE DISCHARGE ABOVE TILE DRAIN
  !
  !     IFLGD=micropore discharge flag to artificial water table
  !     PSISWD=water potential from water table slope
  !     XN,RCHG*=direction indicator,boundary flag
  !     SLOPE=sin(lateral slope)
  !     DLYR=layer width
  !     WaterTBLSlope=water table slope
  !     PSISWT=water potential driving micropore discharge
  !     PSISoilAirEntry,PSISO=matric,osmotic water potential
  !     DPTH,DTBLY=depth to layer midpoint,artificial water table
  !     DepthInternalWTBL=depth to internal water table
  !     FLWL=micropore discharge to natural+artificial water table
  !     HFLWL=convective heat from dischg to natural+artifl water table
  !     HydroCond_3D=saturated hydraulic conductivity
  !     AreaUnderWaterTBL=fraction of layer below artificial water table
!
  IF(IFLGD.EQ.0.AND.(.not.isclose(RechargRateWTBL,0._r8)))THEN
    PSISWD=XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope(N2,N1))
    PSISWT=AZMIN1(-PSISoilMatricPtmp_vr(N3,N2,N1)-0.03_r8*PSISoilOsmotic_vr(N3,N2,N1) &
      +mGravAccelerat*(SoiDepthMidLay_vr(N3,N2,N1)-DTBLY(N2,N1)) &
      -mGravAccelerat*AZMAX1(SoiDepthMidLay_vr(N3,N2,N1)-DepthInternalWTBL(N2,N1)))

    IF(PSISWT.LT.0.0_r8)PSISWT=PSISWT-PSISWD
    FLWT = PSISWT*HydroCond_3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*(1.0_r8-AreaUnderWaterTBL(N3,N2,N1))&
      /(RechargSubSurf+1.0)*RechargRateWTBL*dts_HeatWatTP
    WatXChange2WatTable(N,M6,M5,M4)  = WatXChange2WatTable(N,M6,M5,M4)+XN*FLWT
    WatXChange2WatTableX(N,M6,M5,M4) = WatXChange2WatTableX(N,M6,M5,M4)+XN*FLWT
    HeatFlow2Soili(N,M6,M5,M4)       = HeatFlow2Soili(N,M6,M5,M4)+cpw*TKSoi1_vr(N3,N2,N1)*XN*FLWT
  ENDIF
!
!     MACROPORE DISCHARGE ABOVE TILE DRAIN
!
!     IFLGDH=macropore discharge flag to artificial water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     WaterTBLSlope=water table slope
!     PSISWTH=water potential driving macropore discharge
!     PSISO=osmotic water potential
!     DPTHH,DTBLY=depth to layer macropore water,artificl water table
!     DepthInternalWTBL=depth to internal water table
!     FLWTH,MacPoreDischarg=macropore discharge unltd,ltd by macropore water
!     HydroCondMacP1=macropore hydraulic conductivity
!     ConvWaterFlowMacP_3D=macropore discharge to artificial water table
!     HFLWL=convective heat from discharge to artificial water table
!     HydroCond_3D=saturated hydraulic conductivity
!     AreaUnderWaterTBL=fraction of layer below artificial water table
!
  IF(IFLGDH.EQ.0.AND.(.not.isclose(RechargRateWTBL,0._r8)).AND.VLMacP1_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    PSISWD  = XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope(N2,N1))
    PSISWTH = -0.03_r8*PSISoilOsmotic_vr(N3,N2,N1)+mGravAccelerat*(DPTHH-DTBLY(N2,N1)) &
      -mGravAccelerat*AZMAX1(DPTHH-DepthInternalWTBL(N2,N1))
    IF(PSISWTH.LT.0.0_r8)PSISWTH = PSISWTH-PSISWD

    FLWTH = PSISWTH*HydroCondMacP1_vr(N3,N2,N1)*AREA(N,N3,N2,N1)*(1.0_r8-AreaUnderWaterTBL(N3,N2,N1)) &
      /(RechargSubSurf+1.0_r8)*RechargRateWTBL*dts_HeatWatTP
    MacPoreDischarg=AMAX1(FLWTH,AZMIN1(-(VLWatMacP1_vr(N3,N2,N1)*dts_wat+ConvWaterFlowMacP_3D(3,N3,N2,N1) &
      -ConvWaterFlowMacP_3D(3,N3+1,N2,N1))))
    ConvWaterFlowMacP_3D(N,M6,M5,M4) = ConvWaterFlowMacP_3D(N,M6,M5,M4)+XN*MacPoreDischarg
    HeatFlow2Soili(N,M6,M5,M4)       = HeatFlow2Soili(N,M6,M5,M4)+cpw*TKSoi1_vr(N3,N2,N1)*XN*MacPoreDischarg
  ENDIF
  end subroutine TileDrain
!------------------------------------------------------------------------------------------

  SUBROUTINE SubSufXChangeWithExtWatTable(NY,NX,N,N1,N2,N3,M4,M5,M6,DPTHH,RechargSubSurf,&
    RechargRateWTBL,XN,AirfMicP,VOLPX2,AirfMacP)
  !
  !subsurface recharge to soil micropore and macropores from external water table
  !it considers the existence of frozen layers
  implicit none
  integer, intent(in) :: NY,NX,N,N1,N2,N3,M4,M5,M6
  real(r8),intent(in) :: RechargRateWTBL,XN,DPTHH,RechargSubSurf
  real(r8),intent(inout):: AirfMicP,VOLPX2,AirfMacP
  real(r8) :: FLWU,FLWUL,FLWUX,FLWUH,FLWUHL  
  real(r8) :: PSISUT,PSISUTH,PSISWD

  !     MICROPORE RECHARGE BELOW WATER TABLE
  !
  !     ActiveLayDepth=active layer depth
  !     AirfMicP=air volume
  !     PSISWD=water potential from water table slope
  !     XN,RCHG*=direction indicator,boundary flag
  !     SLOPE=sin(lateral slope)
  !     DLYR=layer width
  !     WaterTBLSlope=water table slope
  !     PSISUT=water potential driving micropore recharge
  !     PSISoilAirEntry,PSISO=matric,osmotic water potential
  !     DPTH,ExtWaterTable=depth to layer midpoint,natural water table
  !     FLWU,FLWUL=micropore recharge unltd,ltd by micropore air volume
  !     FLWL=micropore recharge from natural water table
  !     HFLWL=convective heat from recharge from natural water table
  !     HydroCond_3D=saturated hydraulic conductivity
  !     AREAU=fraction of layer below natural water table
  !     VLairMicP=air filled 
      IF(SoiDepthMidLay_vr(N3,N2,N1).GE.ExtWaterTable_col(N2,N1)     &
        .AND.ActiveLayDepth(N2,N1).GT.ExtWaterTable_col(N2,N1)   &
        .AND.SoiDepthMidLay_vr(N3,N2,N1).LT.ActiveLayDepth(N2,N1) &
        .AND.(AirfMicP.GT.ZEROS2(N2,N1).OR.SoiBulkDensity_vr(N3,N2,N1).LE.ZERO) &
        .AND.VLairMicP_vr(N3,N2,N1).GT.0.0_r8 &
        .AND.(.not.isclose(RechargRateWTBL,0._r8)))THEN
        PSISWD = XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope(N2,N1))
        PSISUT = AZMAX1(-PSISoilMatricPtmp_vr(N3,N2,N1)-0.03_r8*PSISoilOsmotic_vr(N3,N2,N1)+&
          mGravAccelerat*(SoiDepthMidLay_vr(N3,N2,N1)-ExtWaterTable_col(N2,N1)))
        IF(PSISUT.GT.0.0_r8)PSISUT=PSISUT+PSISWD
        FLWU=PSISUT*HydroCond_3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*&
          AREAU(N3,N2,N1)/(RechargSubSurf+1.0)*RechargRateWTBL*dts_HeatWatTP
        !within a time step, the incoming flow cannot exceed avaiable air-filled pores   
        IF(SoiBulkDensity_vr(N3,N2,N1).GT.ZERO)THEN
          FLWUL=AMIN1(FLWU,AirfMicP)
          FLWUX=AMIN1(FLWU,VOLPX2)
        ELSE
          FLWUL=FLWU
          FLWUX=FLWU
        ENDIF
        WatXChange2WatTable(N,M6,M5,M4)  = WatXChange2WatTable(N,M6,M5,M4)+XN*FLWUL
        WatXChange2WatTableX(N,M6,M5,M4) = WatXChange2WatTableX(N,M6,M5,M4)+XN*FLWUX
        HeatFlow2Soili(N,M6,M5,M4)       = HeatFlow2Soili(N,M6,M5,M4)+cpw*TKSoi1_vr(N3,N2,N1)*XN*FLWUL
!        if(M6==1)then
!          write(*,*)'5HeatFlow2Soili(N,M6,M5,M4)=',HeatFlow2Soili(N,M6,M5,M4)
!        endif
        AirfMicP = AirfMicP-XN*WatXChange2WatTable(N,M6,M5,M4)
        VOLPX2   = VOLPX2-XN*WatXChange2WatTableX(N,M6,M5,M4)
      ENDIF
!
      !     MACROPORE RECHARGE BELOW WATER TABLE
      !
      !     PSISWD=water potential from water table slope
      !     XN,RCHG*=direction indicator,boundary flag
      !     SLOPE=sin(lateral slope)
      !     DLYR=layer width
      !     WaterTBLSlope=water table slope
      !     PSISUTH=water potential driving macropore recharge
      !     PSISO=osmotic water potential
      !     DPTHH,ExtWaterTable=depth to layer macropore water,natural water table
      !     DepthInternalWTBL=depth to internal water table
      !     HydroCondMacP1=macropore hydraulic conductivity
      !     FLWUH,FLWUHL=macropore recharge unltd,ltd by macropore air volume
      !     ConvWaterFlowMacP_3D=macropore discharge to natural water table
      !     HFLWL=convective heat from discharge to natural water table
      !     HydroCond_3D=saturated hydraulic conductivity
      !     AREAU=fraction of layer below natural water table
!
      IF(DPTHH.GT.ExtWaterTable_col(N2,N1)                         & !deeper than water table
        .AND.ActiveLayDepth(N2,N1).GT.ExtWaterTable_col(N2,N1)     & !active layer below water table
        .AND.SoiDepthMidLay_vr(N3,N2,N1).LT.ActiveLayDepth(N2,N1) & !midlayer depth above active water layer
        .AND.AirfMacP.GT.ZEROS2(NY,NX)                         & !macropore has air-filled fraction
        .AND.(.not.isclose(RechargRateWTBL,0.0_r8)))THEN         !recharge is on
        PSISWD                    = XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope(N2,N1))
        PSISUTH                   = -0.03_r8*PSISoilOsmotic_vr(N3,N2,N1)+mGravAccelerat*(DPTHH-ExtWaterTable_col(N2,N1))
        IF(PSISUTH.GT.0.0_r8)PSISUTH = PSISUTH+PSISWD
        FLWUH                     = PSISUTH*HydroCondMacP1_vr(N3,N2,N1)*AREA(N,N3,N2,N1)*&
        AREAU(N3,N2,N1)/(RechargSubSurf+1.0_r8)*RechargRateWTBL*dts_HeatWatTP
        !is *dts_wat needed below?  
        FLWUHL                           = AMIN1(FLWUH,AirfMacP*dts_wat)
        ConvWaterFlowMacP_3D(N,M6,M5,M4) = ConvWaterFlowMacP_3D(N,M6,M5,M4)+XN*FLWUHL
        HeatFlow2Soili(N,M6,M5,M4)       = HeatFlow2Soili(N,M6,M5,M4)+cpw*TKSoi1_vr(N3,N2,N1)*XN*FLWUHL
        AirfMacP                         = AirfMacP-XN*ConvWaterFlowMacP_3D(N,M6,M5,M4)
      ENDIF
    
  end SUBROUTINE SubSufXChangeWithExtWatTable              
!------------------------------------------------------------------------------------------

  subroutine FreezeThawMit(NY,NX,L,N1,N2,N3)
  !
  !Description
  !applying freezethaw in cell (N3,N2,N1)
  implicit none
  integer, intent(in) :: NY,NX,L,N1,N2,N3
  real(r8) :: VLHeatCapacityBX,PSISMX,TFREEZ,VLWatMicP1X
  real(r8) :: VLWatMacP1X,VLHeatCapacityAX,TK1App,MacPIceHeatFlxFrezPt,MacPIceHeatFlxFrez
  real(r8) :: ENGY1,MicPIceHeatFlxFrez,MicPIceHeatFlxFrezPt,VLHeatCapacityX
!
!     FREEZE-THAW IN SOIL LAYER MICROPORE FROM NET CHANGE IN SOIL
!     LAYER HEAT STORAGE
!
!     TFREEZ=micropore freezing temperature
!     PSISoilAirEntry,PSISO=micropore matric,osmotic potential
!     VLWatMicP1*,VLiceMicP1=micropore water,ice volume
!     VLWatMacP1*,VLiceMacP1=macropore water,ice volume
!     VLHeatCapacityX,VLHeatCapacityAX,VLHeatCapacityBX=total soil,micropore,macropore heat capacity
!     VHCM=soil solid volumetric heat capacity
!     TK1*=soil temperature
!     THFLWL=total soil conductive, convective heat flux
!     HeatIrrigation1=subsurface convective heat flux
!     MicPIceHeatFlxFrezPt,MicPIceHeatFlxFrez=latent heat from micro freeze-thaw unltd,ltd by water,ice
!     FIceThawMicP=soil water flux from micropore freeze-thaw
!     MacPIceHeatFlxFrezPt,SoiPLIceHeatFlxFrez=latent heat from macro freeze-thaw unltd,ltd by water,ice
!     FIceThawMacP=soil water flux from macropore freeze-thaw
!
  PSISMX          = PSISoilMatricPtmp_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  TFREEZ          = -9.0959E+04_r8/(PSISMX-LtHeatIceMelt)
  VLWatMicP1X     = VLWatMicP1_vr(N3,N2,N1)+TWatCharge2MicP(N3,N2,N1)+FWatExMacP2MicPi(N3,N2,N1)&
    +FWatIrrigate2MicP1_vr(N3,N2,N1)
  VLWatMacP1X     = VLWatMacP1_vr(N3,N2,N1)+TConvWaterFlowMacP_3D_vr(N3,N2,N1)-FWatExMacP2MicPi(N3,N2,N1)
  ENGY1           = VHeatCapacity1_vr(N3,N2,N1)*TKSoi1_vr(N3,N2,N1)
  VLHeatCapacityX = VHeatCapacitySoilM_vr(N3,N2,N1)+cpw*(VLWatMicP1X+VLWatMacP1X) &
    +cpi*(VLiceMicP1_vr(N3,N2,N1)+VLiceMacP1_vr(N3,N2,N1))
  if(VLHeatCapacityX>0._r8)VLHeatCapacityX=VLHeatCapacityX+cpo*RootMassElm_vr(ielmc,N3,N2,N1)

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    TK1App=(ENGY1+THeatFlow2Soili_vr(N3,N2,N1)+HeatIrrigation1(N3,N2,N1))/VLHeatCapacityX
    IF((TK1App.LT.TFREEZ .AND. VLWatMicP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)) &
      .OR.(TK1App.GT.TFREEZ .AND. VLiceMicP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)))THEN
      VLHeatCapacityAX     = VHeatCapacitySoilM_vr(N3,N2,N1)+cpw*VLWatMicP1X+cpi*VLiceMicP1_vr(N3,N2,N1)
      if(VLHeatCapacityAX>0._r8)VLHeatCapacityAX=VLHeatCapacityAX+cpo*RootMassElm_vr(ielmc,N3,N2,N1)
      MicPIceHeatFlxFrezPt = VLHeatCapacityAX*(TFREEZ-TK1App)/((1.0_r8+6.2913E-03_r8*TFREEZ)*(1.0_r8-0.10_r8*PSISMX))*dts_wat
      IF(MicPIceHeatFlxFrezPt.LT.0.0_r8)THEN
        !thaw
        MicPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1_vr(N3,N2,N1)*dts_wat,MicPIceHeatFlxFrezPt)
      ELSE
        MicPIceHeatFlxFrez=AMIN1(LtHeatIceMelt*VLWatMicP1X*dts_wat,MicPIceHeatFlxFrezPt)
      ENDIF
      FIceThawMicP(N3,N2,N1)=-MicPIceHeatFlxFrez/LtHeatIceMelt
    ELSE
      MicPIceHeatFlxFrez=0.0_r8
      FIceThawMicP(N3,N2,N1)=0.0_r8
    ENDIF
  ELSE
    MicPIceHeatFlxFrez=0.0_r8
    FIceThawMicP(N3,N2,N1)=0.0_r8
  ENDIF
!
!     FREEZE-THAW IN SOIL LAYER MACROPORE FROM NET CHANGE IN SOIL
!     LAYER HEAT STORAGE
!  TFice=frozen temperature
  IF((TK1App.LT.TFice.AND.VLWatMacP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)) &
    .OR.(TK1App.GT.TFice.AND.VLiceMacP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)))THEN
    !there is freeze-thaw 
    VLHeatCapacityBX=cpw*VLWatMacP1X+cpi*VLiceMacP1_vr(L,NY,NX)
    MacPIceHeatFlxFrezPt=VLHeatCapacityBX*(TFREEZ-TK1App)/((1.0_r8+6.2913E-03_r8*TFREEZ) &
        *(1.0_r8-0.10_r8*PSISMX))*dts_wat

    IF(MacPIceHeatFlxFrezPt.LT.0.0_r8)THEN
      !ice thaw
      MacPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMacP1_vr(N3,N2,N1)*dts_wat,MacPIceHeatFlxFrezPt)
    ELSE
      !water freeze
      MacPIceHeatFlxFrez=AMIN1(LtHeatIceMelt*VLWatMacP1X*dts_wat,MacPIceHeatFlxFrezPt)
    ENDIF
    FIceThawMacP(N3,N2,N1)=-MacPIceHeatFlxFrez/LtHeatIceMelt
  ELSE
    MacPIceHeatFlxFrez     = 0.0_r8
    FIceThawMacP(N3,N2,N1) = 0.0_r8
  ENDIF
  SoiPLIceHeatFlxFrez(N3,N2,N1)=MicPIceHeatFlxFrez+MacPIceHeatFlxFrez

  !
  !     TOTAL AND HOURLY ACCUMULATED FREEZE-THAW FLUXES
  !
  !     THAW,TLIceThawMacP=hourly accumulated freeze-thaw flux in micropores,macropores
  !     HTHAW=hourly accumulated freeze-thaw latent heat flux
  !     TWFLXL,TMLiceThawMacP=total accumulated freeze-thaw in micropores,macropores
  !     TLPhaseChangeHeat2Soi1=total latent heat flux from melting

  ! Summarize fluxes for the current iteration
  TLIceThawMicP(N3,N2,N1)          = TLIceThawMicP(N3,N2,N1)+FIceThawMicP(N3,N2,N1)
  TLIceThawMacP(N3,N2,N1)          = TLIceThawMacP(N3,N2,N1)+FIceThawMacP(N3,N2,N1)
  TLPhaseChangeHeat2Soi(N3,N2,N1)  = TLPhaseChangeHeat2Soi(N3,N2,N1)+SoiPLIceHeatFlxFrez(N3,N2,N1)
  TMLiceThawMicP(N3,N2,N1)         = TMLiceThawMicP(N3,N2,N1)+FIceThawMicP(N3,N2,N1)
  TMLiceThawMacP(N3,N2,N1)         = TMLiceThawMacP(N3,N2,N1)+FIceThawMacP(N3,N2,N1)
  TLPhaseChangeHeat2Soi1(N3,N2,N1) = TLPhaseChangeHeat2Soi1(N3,N2,N1)+SoiPLIceHeatFlxFrez(N3,N2,N1)
  end subroutine FreezeThawMit        

end module WatsubMod
