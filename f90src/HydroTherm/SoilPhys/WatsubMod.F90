module WatsubMod
!!
! Description:
! Do water and enerby balance calculation.
! The module diagnoses the mass and energy fluxes associated
! with soil/snow water (vapor, liquid and ice) and energy, and update
! them in redistmod.F90

  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: GravAcceleration=>DAT_CONST_G
  use abortutils,     only: endrun,   print_info
  use ElmIDMod,       only: iEastWestDirection,   iNorthSouthDirection, iVerticalDirection
  use SurfPhysData,   only: InitSurfPhysData, DestructSurfPhysData
  use PerturbationMod, only : check_Soil_Warming,is_warming_layerL  
  use DebugToolMod
  use EcoSIMCtrlMod  
  use minimathmod      
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
  real(r8), parameter :: tiny_wat=1.e-12_r8   !threshold water flux
  integer :: curday,curhour
  logical :: do_warming

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
  REAL(R8):: RainEkReducedKsat(JY,JX)

  REAL(R8) :: TopLayWatVol_col(JY,JX)
  real(r8) :: Qinfl2MicP_col(JY,JX)
  real(r8) :: HeatInfl2Soil(JY,JX)
  real(r8) :: Qinfl2MacP_col(JY,JX)
  real(r8) :: twatmass0(JY,JX)
  integer  :: NUX0(JY,JX)
!  real(r8) :: dtime
! begin_execution

  call PrintInfo('beg watsub')
  if(fixWaterLevel)HeatAdv_scal=0._r8

  do_warming=check_Soil_Warming(iYearCurrent,I)
! this enters the EcoSIM interaction with the soil-water-temperature module
  call LocalCopySoilVars(I,J,NHW,NHE,NVN,NVS,twatmass0)
! 
  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)
!
  call InitSoilHydrauics(NHW,NHE,NVN,NVS)

  !if(lverb)  write(*,*)'DYNAMIC LOOP FOR FLUX CALCULATIONS'
!  dtime=0._r8
  ! iterate for NPH times
  D3320: DO M=1,NPH
    DO NX=NHW,NHE
      DO  NY=NVN,NVS
        NUX0(NY,NX)                  = NUM(NY,NX)
        QDischarM_col(NY,NX)         = 0._r8
        QDrainM_col(NY,NX)           = 0._r8
        Qinflx2SoilM_col(NY,NX)      = 0._r8
        QWatIntLaterFlowM_col(NY,NX) = 0._r8
      ENDDO
    ENDDO
!    dtime=dtime+dts_HeatWatTP
    call FWDCopyTopLayerWatVolitM(NHW,NHE,NVN,NVS,TopLayWatVol_col)
    
    if(lverb)write(*,*)'run surface energy balance model, uses ResistanceLitRLay, update top layer soil moisture'
    call RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistanceLitRLay,RainEkReducedKsat,&
      TopLayWatVol_col,HeatFluxAir2Soi,Qinfl2MicP_col,HeatInfl2Soil,Qinfl2MacP_col)
    
    if(lverb)write(*,*)'prepare update for the other soil layers'
    call CopySoilWatVolIterateM(I,J,M,NHW,NHE,NVN,NVS,TopLayWatVol_col)

    if(lverb)write(*,*)'subsurface 3D flow'    
    call Subsurface3DInternalFlowM(I,J,M,NHW,NHE,NVN,NVS,RainEkReducedKsat,HeatFluxAir2Soi)

    if(lverb)write(*,*)'LateralWatHeatExchIterateM'    
    call LateralWatHeatExchIterateM(I,J,M,NHW,NHE,NVN,NVS,RainEkReducedKsat)

!   update states and fluxes
    DO NX=NHW,NHE
      DO  NY=NVN,NVS
        HeatFlx2Grnd_col(NY,NX) = HeatFlx2Grnd_col(NY,NX)+HeatInfl2Soil(NY,NX)
        Qinflx2Soil_col(NY,NX)  = Qinflx2Soil_col(NY,NX)+Qinfl2MicP_col(NY,NX)+Qinfl2MacP_col(NY,NX)
        Qinflx2SoilM_col(NY,NX) = Qinflx2SoilM_col(NY,NX)+Qinfl2MicP_col(NY,NX)+Qinfl2MacP_col(NY,NX)
      ENDDO
    ENDDO  

    if(snowRedist_model)call AccumulateSnowRedisFluxM(I,J,M,NHW,NHE,NVN,NVS)

    call UpdateSoilMoistTemp(I,J,M,NHW,NHE,NVN,NVS,twatmass0)

    call AggregateSurfRunoffFluxM(I,J,M,NHW,NHE,NVN,NVS)

    call UpdateSurfaceAtM(I,J,M,NHW,NHE,NVN,NVS)

    IF(M.NE.NPH)THEN
      call UpdateStateFluxAtM(I,J,M,NHW,NHE,NVN,NVS)
    ELSE
      call UpdateFluxAtExit(I,J,NHW,NHE,NVN,NVS)
    ENDIF

    call checkMassBalance(I,J,M,NHW,NHE,NVN,NVS,twatmass0,NUX0)

  ENDDO D3320
  call PrintInfo('end watsub')

  END subroutine watsub

!------------------------------------------------------------------------------------------  

  subroutine checkMassBalance(I,J,M,NHW,NHE,NVN,NVS,twatmass0,NUX0)

  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(inout) :: twatmass0(JY,JX)
  integer, intent(in)    :: NUX0(JY,JX)
  real(r8)  :: twatmass1(JY,JX)

  real(r8) :: dwat,dwat0
  integer :: NY,NX,L

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
     twatmass1(NY,NX)=0._r8
     dwat0=0._r8
     DO L=NUX0(NY,NX),NUM(NY,NX)-1
       dwat0=dwat0 + VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
     ENDDO
      D131: DO L=NUM(NY,NX),NL(NY,NX)
        twatmass1(NY,NX)=twatmass1(NY,NX)+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
      ENDDO D131
      if(fixWaterLevel)then
        dwat=twatmass0(NY,NX)-twatmass1(NY,NX)
      else
        dwat=twatmass0(NY,NX)-twatmass1(NY,NX)+Qinflx2Soil_col(NY,NX)+QWatIntLaterFlow_col(NY,NX)-QDischar_col(NY,NX)-QDrain_col(NY,NX)
      endif
!      dwat=twatmass0(NY,NX)-twatmass1(NY,NX)+Qinflx2SoilM_col(NY,NX)-QDischarM_col(NY,NX)-QDrainM_col(NY,NX)+QWatIntLaterFlowM_col(NY,NX)
!      if(I==141 .and. J>=13)then
!        write(211,*)I+J/24.,NY,NX,M,'wat',NUM(NY,NX),dwat,twatmass0(NY,NX),twatmass1(NY,NX),Qinflx2Soil_col(NY,NX),QWatIntLaterFlow_col(NY,NX),&
!          QDischar_col(NY,NX),QDrain_col(NY,NX)
!        write(211,*)I+J/24.,NY,NX,M,'wat',NUM(NY,NX),dwat,twatmass0(NY,NX),twatmass1(NY,NX),Qinflx2SoilM_col(NY,NX),QWatIntLaterFlowM_col(NY,NX),&
!          QDischarM_col(NY,NX),QDrainM_col(NY,NX)          
!        if(NX==1)then
!          write(311,*)I+J/24.,NY,NX,M,'wat',dwat,twatmass0(NY,NX),twatmass1(NY,NX),Qinflx2SoilM_col(NY,NX),QDischarM_col(NY,NX),&
!            QDrainM_col(NY,NX),QWatIntLaterFlowM_col(NY,NX)  
!          write(401,*)I+J/24.,NY,NX,M,'chk', (VLWatMicP1_vr(L,NY,NX),L=NUM(NY,NX),NL(NY,NX)) 
!        else
!          write(312,*)I+J/24.,NY,NX,M,'wat',dwat,twatmass0(NY,NX),twatmass1(NY,NX),Qinflx2SoilM_col(NY,NX),QDischarM_col(NY,NX),&
!            QDrainM_col(NY,NX),QWatIntLaterFlowM_col(NY,NX)  
!        endif  
!        write(213,*)NY,NX,NUM(NY,NX),NU(NY,NX),VLWatMicP_vr(1,NY,NX),VLWatMicP1_vr(1,NY,NX),H2OFlow2TopSoiMicP_col(NY,NX)  
!      endif
!       if(I>=156)write(211,*)I*1000+J,NY,NX,M,'wat',NUM(NY,NX),dwat,twatmass0(NY,NX),twatmass1(NY,NX),Qinflx2Soil_col(NY,NX),QWatIntLaterFlow_col(NY,NX),&
!          QDischar_col(NY,NX),QDrain_col(NY,NX),dwat0

!      if(I==358 .and. J==12)write(211,*)I+J/24.,NY,NX,M,'wat',NUX0(NY,NX),NUM(NY,NX),dwat,twatmass0(NY,NX),twatmass1(NY,NX),Qinflx2Soil_col(NY,NX),QWatIntLaterFlow_col(NY,NX),&
!          QDischar_col(NY,NX),QDrain_col(NY,NX)
      if(abs(dwat)>1.e-4_r8)then
        call endrun('soil H2O error test failure in '//trim(mod_filename)//' at line',__LINE__)
      endif

!      twatmass0(NY,NX)=twatmass1(NY,NX)
    ENDDO
  ENDDO  

  end subroutine checkMassBalance
!------------------------------------------------------------------------------------------  

  subroutine FWDCopyTopLayerWatVolitM(NHW,NHE,NVN,NVS,TopLayWatVol_col)
  !
  !make a copy of toplayer water content
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(out) :: TopLayWatVol_col(JY,JX)
  integer :: NY,NX

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      TopLayWatVol_col(NY,NX)= VLWatMicP1_vr(NUM(NY,NX),NY,NX)
    ENDDO
  ENDDO
  end subroutine FWDCopyTopLayerWatVolitM

!------------------------------------------------------------------------------------------  

  subroutine CopySoilWatVolIterateM(I,J,M,NHW,NHE,NVN,NVS,TopLayWatVol_col)

  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: TopLayWatVol_col
  integer :: L,NY,NX
  
  !VOLW2 will be updated in soil surface model

  DO NX=NHW,NHE
    DO  NY=NVN,NVS  
      DO L=NUM(NY,NX)+1,NL(NY,NX)        
        VLWatMicP2_vr(L,NY,NX)=VLWatMicP1_vr(L,NY,NX)
      ENDDO  
      VLWatMicP2_vr(NUM(NY,NX),NY,NX)=TopLayWatVol_col(NY,NX)
    ENDDO
  ENDDO
  end subroutine CopySoilWatVolIterateM


!------------------------------------------------------------------------------------------  
  subroutine LocalCopySoilVars(I,J,NHW,NHE,NVN,NVS,twatmass0)
  !
  !coupling interface between EcoSIM and the soil physics module 
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(out):: twatmass0(JY,JX)
  integer :: NY,NX

  integer :: L,LyrIrrig
  real(r8) :: VLTSoiPore,rsatMacP

  DX995: DO NX=NHW,NHE
    DX990: DO NY=NVN,NVS
      twatmass0(NY,NX)=0._r8
    !make a local copy of the upper boundary index
      NUM(NY,NX)=NU(NY,NX)

    ! CDPTH=depth to bottom of soil layer
    ! WDPTH,LyrIrrig=depth,layer of subsurface irrigation

      !identify the layer where irrigation is applied
      D65: DO L=NUM(NY,NX),NL(NY,NX)
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.WDPTH(I,NY,NX))THEN
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
    !   CCLAY_vr=clay concentration
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
        twatmass0(NY,NX)        = twatmass0(NY,NX)+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
          VLairMicP_vr(L,NY,NX)  = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
          VLairMicP1_vr(L,NY,NX) = AZMAX1(VLairMicP_vr(L,NY,NX))
        ELSE
          VLairMicP_vr(L,NY,NX)  = 0.0_r8
          VLairMicP1_vr(L,NY,NX) = 0.0_r8
        ENDIF

        !FVOLAH accounts for clay swelling effect due to change in micropore water, but it is set to zero
        VLMacP1_vr(L,NY,NX)=AZMAX1(VLMacP_vr(L,NY,NX)-FVOLAH*CCLAY_vr(L,NY,NX) &
          *(safe_adb(VLWatMicP1_vr(L,NY,NX),VLSoilMicP_vr(L,NY,NX))-WiltPoint_vr(L,NY,NX))*VGeomLayer_vr(L,NY,NX))

        IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
          VLairMacP_vr(L,NY,NX)  = VLMacP1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
          VLairMacP1_vr(L,NY,NX) = AZMAX1(VLairMacP_vr(L,NY,NX))
        ELSE
          VLairMacP_vr(L,NY,NX)  = 0.0_r8
          VLairMacP1_vr(L,NY,NX) = 0.0_r8
        ENDIF        
        VLWatMicPM_vr(1,L,NY,NX) = VLWatMicP1_vr(L,NY,NX)
        VLWatMacPM_vr(1,L,NY,NX)    = VLWatMacP1_vr(L,NY,NX)
        VLsoiAirPM_vr(1,L,NY,NX) = VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX)+&
          THETPI*(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))
        
        VLTSoiPore = VLSoilMicP_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)
        IF(VLTSoiPore.GT.ZEROS2(NY,NX))THEN
          !fraction as water
          FracSoiPAsWat_vr(L,NY,NX)=AZMAX1t((VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX))/VLTSoiPore)
          !fraction as ice
          FracSoiPAsIce_vr(L,NY,NX)=AZMAX1t((VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))/VLTSoiPore)
          !fraction as air
          AirFilledSoilPore_vr(L,NY,NX)=AZMAX1t((VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX))/VLTSoiPore)
        ELSE
          FracSoiPAsWat_vr(L,NY,NX)=POROS_vr(L,NY,NX)
          FracSoiPAsIce_vr(L,NY,NX)=0.0_r8
          AirFilledSoilPore_vr(L,NY,NX)=0.0_r8
        ENDIF
        AirFilledSoilPoreM_vr(1,L,NY,NX)=AirFilledSoilPore_vr(L,NY,NX)
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
        if(plantOM4Heat .and. VLHeatCapacityA_vr(L,NY,NX)>0._r8)then
          VLHeatCapacityA_vr(L,NY,NX)=VLHeatCapacityA_vr(L,NY,NX)+cpo*RootMassElm_vr(ielmc,L,NY,NX)      
        endif  
        VLHeatCapacityB_vr(L,NY,NX) = cpw*VLWatMacP1_vr(L,NY,NX)+cpi*VLiceMacP1_vr(L,NY,NX)
        VHeatCapacity1_vr(L,NY,NX)  = VLHeatCapacityA_vr(L,NY,NX)+VLHeatCapacityB_vr(L,NY,NX)
    !
    !   MACROPOROSITY
    !
    !   FMAC,SoilFracAsMicP=macropore,micropore volume fractions
    !   CNDH*=macropore hydraulic conductivity
    !   TKS_vr,TK1=soil temperature
    !   FLU,HeatIrrigation=subsurface water,convective heat fluxes
    !   FracLayVolBelowExtWTBL_vr,AREAD=fractions of layer below natural,artifl water table
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
        TKSoil1_vr(L,NY,NX)         = TKS_vr(L,NY,NX)

        !LyrIrrig=layer number where irrigation is applied
        IF(L.EQ.LyrIrrig .and. .not.fixWaterLevel)THEN
          FWatIrrigate2MicP_vr(L,NY,NX)  = IrrigSubsurf_col(NY,NX)
          HeatIrrigation_vr(L,NY,NX)     = cpw*TairK_col(NY,NX)*IrrigSubsurf_col(NY,NX)
          FWatIrrigate2MicP1_vr(L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*dts_HeatWatTP
          HeatIrrigation1_vr(L,NY,NX)    = HeatIrrigation_vr(L,NY,NX)*dts_HeatWatTP
        ELSE
          FWatIrrigate2MicP_vr(L,NY,NX)  = 0.0_r8
          HeatIrrigation_vr(L,NY,NX)     = 0.0_r8
          FWatIrrigate2MicP1_vr(L,NY,NX) = 0.0_r8
          HeatIrrigation1_vr(L,NY,NX)    = 0.0_r8
        ENDIF
        !lower than natural water table
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.ExtWaterTable_col(NY,NX))THEN
          FracLayVolBelowExtWTBL_vr(L,NY,NX)=AMIN1(1.0_r8, &
            AZMAX1(safe_adb(CumDepz2LayBottom_vr(L,NY,NX)-ExtWaterTable_col(NY,NX),DLYR_3D(3,L,NY,NX))))
        ELSE
          FracLayVolBelowExtWTBL_vr(L,NY,NX)=0.0_r8
        ENDIF
        !
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.TileWaterTable_col(NY,NX))THEN
          FracLayVolBelowTileWTBL_vr(L,NY,NX)=AMIN1(1.0_r8, &
            AZMAX1(safe_adb(CumDepz2LayBottom_vr(L,NY,NX)-TileWaterTable_col(NY,NX),DLYR_3D(3,L,NY,NX))))
        ELSE
          FracLayVolBelowTileWTBL_vr(L,NY,NX)=0.0_r8
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
    !     AVCNHL_3D=macropore hydraulic conductance
    !     DLYR=layer depth
    !
          IF(HydroCondMacP1_vr(N3,N2,N1).GT.ZERO .AND. HydroCondMacP1_vr(N6,N5,N4).GT.ZERO)THEN
            !when both src and dest have macro pores
            AVCNHL_3D(N,N6,N5,N4)=2.0_r8*HydroCondMacP1_vr(N3,N2,N1)*HydroCondMacP1_vr(N6,N5,N4) &
              /(HydroCondMacP1_vr(N3,N2,N1)*DLYR_3D(N,N6,N5,N4)+HydroCondMacP1_vr(N6,N5,N4) &
              *DLYR_3D(N,N3,N2,N1))
          ELSE
            AVCNHL_3D(N,N6,N5,N4)=0.0_r8
          ENDIF
        ENDDO D40
      ENDDO D35
    ENDDO
  ENDDO
  end subroutine InitSoilHydrauics

!------------------------------------------------------------------------------------------

  subroutine Subsurface3DInternalFlowM(I,J,M,NHW,NHE,NVN,NVS,KSatRedusByRainKinetEnergy,HeatFluxAir2Soi)
  !
  !Description
  !Internal 3D flow
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in)  :: M,NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: KSatRedusByRainKinetEnergy(:,:)
  real(r8), dimension(:,:),intent(in) :: HeatFluxAir2Soi(:,:)
  integer :: N,N1,N2,N3,N4,N5,N6,L,LL,K1,KL,NY,NX
  real(r8) :: WTHET1,FCDX,FCLX,FCX
  real(r8) :: PSISV1,TKY,PSDX
  real(r8) :: WPLX,WPX,WTHET2

  real(r8) :: ConvectVapFlux,HeatByConvectVapFlux,PSISVL
  real(r8) :: PredDarcyFlowMax
  real(r8) :: TKLX
  real(r8) :: WatDarcyFlowMicP,HeatByDarcyFlowMicP,WaterMacpFlow
  real(r8) :: HeatByFlowMacP,HeatByWatFlowMicP,HFLWS,THETA1,THETAL  
  logical  :: LInvalidMacP     !disable macropore?
  !     begin_execution
  !
  !     WATER AND ENERGY TRANSFER THROUGH SOIL PROFILE
  !
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !

  call InitSoil3DModelIterateM(I,J,M,NHW,NHE,NVN,NVS)
  
  !If waterlevel is fixed, like shallow lake, without considering the change
  !of water depth due to hydrological fluxes 
  
  DO NX=NHW,NHE
    DO NY=NVN,NVS

      LInvalidMacP=.false.
      D4400: DO L=1,NL(NY,NX)
        N1=NX;N2=NY;N3=L
    !
    !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
    !
        D4320: DO N=FlowDirIndicator(N2,N1),3
          IF(N.EQ.iEastWestDirection)THEN            
            !west to east
            IF(NX.EQ.NHE)THEN
              !east boundary
              cycle
            ELSE
              N4=NX+1
              N5=NY
              N6=L  
          !     ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
          !
          !     IF(N2.EQ.2.AND.(N1.EQ.2.OR.N1.EQ.3).AND.L.LE.15)THEN
          !     CYCLE
          !     ENDIF
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN
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
            !bottom layer
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
!         write(*,*)'SKIP NON-EXISTENT DESTINATION SOIL LAYERS'
          ! identified by soil volume
          D1100: DO LL=N6,NL(NY,NX)
            IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
              N6=LL
              exit
            ENDIF
          ENDDO D1100
          ! source grid at surface, make a copy of the destination layer number
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
            !Both source and destination grids are legitimate
            !Source top layer is matched to destination top layer 
            IF(N3.GE.NUM(N2,N1) .AND. N6.GE.NUM(N5,N4) .AND. N3.LE.NL(N2,N1) .AND. N6.LE.NL(N5,N4))THEN

              call CalcSoilWatPotential(NY,NX,N1,N2,N3,PSISoilMatricPtmp_vr(N3,N2,N1),THETA1)

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
              !     (CURRENT WATER POTENTIAL > AIR ENTRY WATER POTENTIAL)
              
              call MicporeDarcyFlow(NY,NX,N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,&
                KSatRedusByRainKinetEnergy(NY,NX),PredDarcyFlowMax,WatDarcyFlowMicP,&
                HeatByDarcyFlowMicP,PSISV1,PSISVL)          
                
              call MacporeFLow(NY,NX,M,N,N1,N2,N3,N4,N5,N6,WaterMacpFlow,HeatByFlowMacP,LInvalidMacP)

              call WaterVaporFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectVapFlux,&
                HeatByConvectVapFlux)

          !
          !     FLWL=total water+vapor flux to destination
          !     WaterFlow2MicptX_3D=total unsaturated water+vapor flux to destination
          !     HeatByWatFlowMicP=total convective heat flux from water+vapor flux
          !
              WaterFlow2Micpt_3D(N,N6,N5,N4)  = WatDarcyFlowMicP+ConvectVapFlux
              WaterFlow2MicptX_3D(N,N6,N5,N4) = PredDarcyFlowMax+ConvectVapFlux
              WaterFlow2Macpt_3D(N,N6,N5,N4)  = WaterMacpFlow              
              HeatFlow2Soili_3D(N,N6,N5,N4)   = HeatByDarcyFlowMicP+HeatByConvectVapFlux+HeatByFlowMacP
!              if(N6==16)then
!              write(411,*)I+J/24.,M,'heatdarcy',HeatByDarcyFlowMicP,HeatByConvectVapFlux,HeatByFlowMacP,N6,N5,N4,N
!              endif
              WaterFlow2Micptl_3D(N,N6,N5,N4) = WaterFlow2Micpt_3D(N,N6,N5,N4)
              WaterFlow2Macptl_3D(N,N6,N5,N4) = WaterFlow2Macpt_3D(N,N6,N5,N4)
          !   compute heat flux by conduction
              call Solve4HeatByConduction(I,J,M,N,NY,NX,N1,N2,N3,N4,N5,N6,HeatByConvectVapFlux,HeatFluxAir2Soi(NY,NX))

          !
          !     TOTAL WATER, VAPOR AND HEAT FLUXES
          !
          !     FLW,FLWX,FLWH=total water flux through micropores,macropores
          !     HFLW=total heat flux
          !     WaterFlow2MicPM_3D=water flux used for solute flux calculations in TranspNoSalt.f
          !
              if(.not.fixWaterLevel)then
                WaterFlowSoiMicP_3D(N,N6,N5,N4)  = WaterFlowSoiMicP_3D(N,N6,N5,N4)+WaterFlow2Micpt_3D(N,N6,N5,N4)
                WaterFlowSoiMicPX_3D(N,N6,N5,N4) = WaterFlowSoiMicPX_3D(N,N6,N5,N4)+WaterFlow2MicptX_3D(N,N6,N5,N4)
                WaterFlowSoiMacP_3D(N,N6,N5,N4)  = WaterFlowSoiMacP_3D(N,N6,N5,N4)+WaterFlow2Macpt_3D(N,N6,N5,N4)
              endif
              HeatFlow2Soil_3D(N,N6,N5,N4)     = HeatFlow2Soil_3D(N,N6,N5,N4)+HeatFlow2Soili_3D(N,N6,N5,N4)
              if(N.NE.iVerticalDirection)then
                QWatIntLaterFlow_col(N2,N1)  = QWatIntLaterFlow_col(N2,N1)-WaterFlow2Micpt_3D(N,N6,N5,N4)-WaterFlow2Macpt_3D(N,N6,N5,N4)
                QWatIntLaterFlowM_col(N2,N1) = QWatIntLaterFlowM_col(N2,N1)-WaterFlow2Micpt_3D(N,N6,N5,N4)-WaterFlow2Macpt_3D(N,N6,N5,N4)
              endif
              WaterFlow2MicPM_3D(M,N,N6,N5,N4) = WaterFlow2Micpt_3D(N,N6,N5,N4)
              WaterFlow2MacPM_3D(M,N,N6,N5,N4) = WaterFlow2Macpt_3D(N,N6,N5,N4)                            
              IF(N.EQ.iVerticalDirection)THEN
            !
            !     WATER FILM THICKNESS FOR CALCULATING GAS EXCHANGE IN TranspNoSalt.F
            !
                FILMM_vr(M,N6,N5,N4)=FilmThickness(PSISoilMatricPtmp_vr(N6,N5,N4))
              ENDIF
            ELSEIF(N.NE.iVerticalDirection)THEN
              WaterFlow2Micpt_3D(N,N6,N5,N4)   = 0.0_r8
              WaterFlow2MicptX_3D(N,N6,N5,N4)  = 0.0_r8
              WaterFlow2Macpt_3D(N,N6,N5,N4)   = 0.0_r8
              HeatFlow2Soili_3D(N,N6,N5,N4)    = 0.0_r8
              WaterFlow2MicPM_3D(M,N,N6,N5,N4) = 0.0_r8
              WaterFlow2MacPM_3D(M,N,N6,N5,N4) = 0.0_r8
            ENDIF

          ELSE
            IF(N.EQ.iVerticalDirection)THEN
              WaterFlow2Micpt_3D(N,N3,N2,N1)   = 0.0_r8
              WaterFlow2MicptX_3D(N,N3,N2,N1)  = 0.0_r8
              WaterFlow2Macpt_3D(N,N3,N2,N1)   = 0.0_r8
              HeatFlow2Soili_3D(N,N3,N2,N1)    = 0.0_r8
              WaterFlow2MacPM_3D(M,N,N3,N2,N1) = 0.0_r8
              WaterFlow2MacPM_3D(M,N,N3,N2,N1) = 0.0_r8
            ELSE
              WaterFlow2Micpt_3D(N,N6,N5,N4)   = 0.0_r8
              WaterFlow2MicptX_3D(N,N6,N5,N4)  = 0.0_r8
              WaterFlow2Macpt_3D(N,N6,N5,N4)   = 0.0_r8
              HeatFlow2Soili_3D(N,N6,N5,N4)    = 0.0_r8
              WaterFlow2MicPM_3D(M,N,N6,N5,N4) = 0.0_r8
              WaterFlow2MacPM_3D(M,N,N6,N5,N4) = 0.0_r8
            ENDIF
          ENDIF
        ENDDO D4320
      ENDDO D4400
    ENDDO
  ENDDO  
  !only do this for 

  DO N=FlowDirIndicator(N2,N1),2
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        N1=NX;N2=NY
        DO L=NU(N2,N1),NL(NY,NX)
          N3=L    
          QWatIntLaterFlow_col(N2,N1)  = QWatIntLaterFlow_col(N2,N1)+WaterFlow2Micptl_3D(N,N3,N2,N1)+WaterFlow2Macptl_3D(N,N3,N2,N1)
          QWatIntLaterFlowM_col(N2,N1) = QWatIntLaterFlowM_col(N2,N1)+WaterFlow2Micptl_3D(N,N3,N2,N1)+WaterFlow2Macptl_3D(N,N3,N2,N1)
        ENDDO        
      ENDDO  
    ENDDO
  ENDDO

  end subroutine Subsurface3DInternalFlowM
!------------------------------------------------------------------------------------------

  subroutine LateralWatHeatExchIterateM(I,J,M,NHW,NHE,NVN,NVS,RainEkReducedKsat)
  !
  !Description
  ! boundary flow involes exchange with external water table, and through tile drainage
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: RainEkReducedKsat(JY,JX)
  integer :: NY,NX
  integer :: L,LL
  integer :: N,NN,N1,N2,N3
  integer :: N4,N5,N4B,N5B,N6           !grid indices for snow redistribution
  integer :: M1,M2,M3,M4,M5,M6,K1,KL
  logical :: DoMacPDischarg2Tile
  logical :: DoMicPDischarg2ExtWTBL
  logical :: DoMacPDischarg2ExtWTBL
  logical :: DoMicPDischarg2Tile
  real(r8) :: VLWatMicP1X
  real(r8) :: FLWT
  real(r8) :: RechargSurf,RechargSubSurf,RechargRateWTBL
  real(r8) :: DPTHH,CNDL
  real(r8) :: FINHX,THETAX  
  real(r8) :: AirfMicP,VOLPX2,AirfMacP
  real(r8) :: XN,THETA1,HeatFlx
  real(r8) :: VOLP1X,VLWatMacP1X,VOLPH1X  
  character(len=20) :: str_dir
  logical :: donot_drain,lZeroRunoff

!     begin_execution
!     AirfMicP,AirfMacP=air-filled porosity in micropores,macropores
  D9595: DO  NX=NHW,NHE
    D9590: DO  NY=NVN,NVS
!      write(110,*)I+J/24.,M,FlowDirIndicator(NY,NX)
      D9585: DO L=NUM(NY,NX),NL(NY,NX)

        AirfMicP = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
        VOLPX2   = AirfMicP
        AirfMacP = VLMacP1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
        DPTHH    = get_DPTHH(L,NY,NX)

        ! configure water table info
        call Config4WaterTableDrain(L,NY,NX,DPTHH,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL)

        call Config4TileDrainage(L,NY,NX,DPTHH,DoMicPDischarg2Tile,DoMacPDischarg2Tile)
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
!     N3,N2,N1 for source grid cell
!     M6,M5,M4 for destination grid cell
!       
        N1=NX;N2=NY;N3=L
!
!     LOCATE EXTERNAL BOUNDARIES
!
        D9580: DO N=FlowDirIndicator(NY,NX),3
          D9575: DO NN=1,2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            !Locate the boundary 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
            str_dir='none'       
            IF(N.EQ.iEastWestDirection)THEN
              ! along the W-E direction              
              N4  = NX+1; N5  = NY  !eastern
              N4B = NX-1; N5B = NY  !western
              N6  = L
              IF(NN.EQ.iOutflow)THEN
                !eastern boundary  |---|->
                IF(NX.EQ.NHE)THEN                
                  M1 = NX; M2  = NY; M3 = L  !source
                  M4 = NX+1;M5 = NY;M6  = L  !target
                  XN = -1.0_r8               !going out of eastern boundary
                  RechargSurf     = RechargEastSurf(M2,M1)
                  RechargSubSurf  = RechargEastSubSurf(M2,M1)
                  RechargRateWTBL = RechargRateEastWTBL(M2,M1)
                  str_dir='east'
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN
                !west boundary   -|-> |
                IF(NX.EQ.NHW)THEN                  
                  M1 = NX+1;M2 = NY; M3 = L  !source                  
                  M4 = NX; M5 = NY; M6 = L   !target
                  XN = 1.0_r8                !coming in from western boundary
                  RechargSurf     = RechargWestSurf(M5,M4)
                  RechargSubSurf  = RechargWestSubSurf(M5,M4)
                  RechargRateWTBL = RechargRateWestWTBL(M5,M4)
                  str_dir='west'
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              ! along the N-S direction              
              N4  = NX; N5  = NY+1 !south
              N4B = NX; N5B = NY-1 !north
              N6  = L
              IF(NN.EQ.iOutflow)THEN     !-----
                !south boundary     \|/
                IF(NY.EQ.NVS)THEN !-----
                  M1 = NX; M2 = NY; M3   = L   !source
                  M4 = NX; M5 = NY+1; M6 = L   !target
                  XN              = -1.0_r8    !going out of south boundary
                  RechargSurf     = RechargSouthSurf(M2,M1)
                  RechargSubSurf  = RechargSouthSubSurf(M2,M1)
                  RechargRateWTBL = RechargRateSouthWTBL(M2,M1)
                  str_dir='south'
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN  !\|/
                !north boundary     ----
                IF(NY.EQ.NVN)THEN  !----
                  M1 = NX; M2 = NY+1; M3 = L !source grid 
                  M4 = NX;M5  = NY; M6   = L !target grid
                  XN          = 1.0_r8       !coming in from north boundary
                  RechargSurf     = RechargNorthSurf(M5,M4)
                  RechargSubSurf  = RechargNorthSubSurf(M5,M4)
                  RechargRateWTBL = RechargRateNorthWTBL(M5,M4)
                  str_dir='north'
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN
! in the vertical direction
              N4 = NX
              N5 = NY
              N6 = L+1
              IF(NN.EQ.iOutflow)THEN
                !bottom
                IF(L.EQ.NL(NY,NX))THEN
                  M1 = NX; M2 = NY; M3  = L   !source
                  M4 = NX; M5 = NY; M6 = L+1  !target
                  XN              = -1.0_r8    !going out lower boundary
                  RechargSubSurf  = RechargBottom_col(M2,M1)
                  RechargRateWTBL = 1.0_r8
                  str_dir='vert'
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN
                cycle
              ENDIF
            ENDIF
!
!     BOUNDARY SURFACE RUNOFF DEPENDING ON ASPECT, SLOPE
!     VELOCITY, HYDRAULIC RADIUS AND SURFACE WATER STORAGE
!
!     CDPTH,CumLitRDepzInit_col=current,initial surface elevation
!     BKDS=bulk density
!     XGridRunoffFlag,RCHQ*=runoff boundary flags

!           top soil layer and surface soil layer is active, litter layer is lower than its initial thickness 
!           or the grid is a soil
!           surface lateral flow 
            IF(L.EQ.NUM(N2,N1) .AND. N.NE.iVerticalDirection                                & ! lateral flow at surface
              .AND. (CumDepz2LayBottom_vr(NU(N2,N1)-1,N2,N1).LE.CumLitRDepzInit_col(N2,N1)  & ! in the soil
              .OR. SoilBulkDensity_vr(NUI(N2,N1),N2,N1).GT.ZERO))THEN                         ! it is soil

              !  NO runoff
              lZeroRunoff=.not.XGridRunoffFlag(NN,N,N2,N1)   &              !Runoff flag is off
                .OR. isclose(RechargSurf,0._r8)              &              !Allowed runoff rate is zero
                .OR. ABS(SurfRunoffWatFluxM_2DH(M,N2,N1)).LT.ZEROS(N2,N1)   !Runoff from infiltration partition is inisignificant

               !do runoff
              IF(.not.LZeroRunoff)then
                call XBoundSurfaceRunoff(I,J,M,N,NN,N1,N2,M4,M5,RechargSurf,XN)
                !
                !     BOUNDARY SNOW FLUX
                !
                !  DrySnoFlxBySnowRedistribut,WatFlxBySnowRedistribut,IceFlxBySnowRedistribut=snow,water,ice transfer
                !  HeatFlxBySnowRedistribut=convective heat transfer from snow,water,ice transfer

                ! Zero out for eastern and southern
                IF(NN.EQ.iOutflow)THEN
                  DrySnoFlxBySnowRedistribut(N,M5,M4)  = 0.0_r8
                  WatFlxBySnowRedistribut(N,M5,M4)     = 0.0_r8
                  IceFlxBySnowRedistribut(N,M5,M4)     = 0.0_r8
                  HeatFlxBySnowRedistribut(N,M5,M4)    = 0.0_r8
                  DrySnoFlxBySnoRedistM_2DH(M,N,M5,M4) = DrySnoFlxBySnowRedistribut(N,M5,M4)
                ENDIF
              ENDIF
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! BOUNDARY SUBSURFACE WATER AND HEAT TRANSFER DEPENDING
          ! ON LEVEL OF WATER TABLE
!
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
              IF(FlowDirIndicator(N2,N1).NE.iVerticalDirection .OR. N.EQ.iVerticalDirection)THEN
                !including lateral connection or woking on vertical direction

                ! IF NO WATER TABLE

                IF(IDWaterTable_col(N2,N1).EQ.0 .OR. N.EQ.iVerticalDirection)THEN    
                  !vertical flow or lateral drainage
                  call BoundaryDrain(I,J,N,N1,N2,N3,M4,M5,M6,RainEkReducedKsat(NY,NX),XN,RechargSubSurf,RechargRateWTBL)

                  !lateral flow
                ELSE

                  WaterFlow2Micpt_3D(N,M6,M5,M4)  = 0.0_r8
                  WaterFlow2MicptX_3D(N,M6,M5,M4) = 0.0_r8
                  WaterFlow2Macpt_3D(N,M6,M5,M4)  = 0.0_r8
                  HeatFlow2Soili_3D(N,M6,M5,M4)   = 0.0_r8

                  if(DoMicPDischarg2ExtWTBL.or.DoMacPDischarg2ExtWTBL)then
                    !discharge to external water table
                    CALL DischargeOverWaterTBL(I,J,N,N1,N2,N3,M4,M5,M6,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL,&
                      RechargSubSurf,RechargRateWTBL,DPTHH,XN,ExtWaterTable_col(N2,N1),FracLayVolBelowExtWTBL_vr(N3,N2,N1))
                  endif

                  if(DoMicPDischarg2Tile .or. DoMacPDischarg2Tile)then
                    !tile drainage
                    CALL DischargeOverWaterTBL(I,J,N,N1,N2,N3,M4,M5,M6,DoMicPDischarg2Tile,DoMacPDischarg2Tile,&
                      RechargSubSurf,RechargRateWTBL,DPTHH,XN,TileWaterTable_col(N2,N1),FracLayVolBelowTileWTBL_vr(N3,N2,N1))                      
                  endif
!                  write(214,*)I+J/24.,M,N3,N2,N1,M6,M5,M4,'chking str'//trim(str_dir)
!                  donot_drain=(I==267 .and. J==24 .and. M==1 .and. trim(str_dir)=='south')

                  call RechargeFromExtWaterTBL(I,J,M,N,N1,N2,N3,M4,M5,M6,DPTHH,RechargSubSurf,&
                    RechargRateWTBL,XN,AirfMicP,VOLPX2,AirfMacP,str_dir)!,donot_drain)

                ENDIF
!
        !     SUBSURFACE HEAT SOURCE/SINK
        !
        !     HFLWL=heat flux across lower boundary
        !     TK1=lower boundary soil temperature
        !     TKSD=deep source/sink temperature from geothermal flux
        !     TCNDG=thermal conductivity below lower boundary
        !     SoilHeatSrcDepth_col,CDPTH=depth of thermal sink/source, lower boundary
        !     KoppenClimZone=Koppen climate zone
                IF(N.EQ.iVerticalDirection .AND. KoppenClimZone_col(N2,N1).NE.-2)THEN
                !heat flux going out (<0)
                  HeatFlx=(TKSoil1_vr(N3,N2,N1)-TKSD(N2,N1))*TCNDG/(SoilHeatSrcDepth_col(N2,N1)-CumDepz2LayBottom_vr(N3,N2,N1)) &
                    *AREA(N,N3,N2,N1)*dts_HeatWatTP
                  HeatFlow2Soili_3D(N,M6,M5,M4) = HeatFlow2Soili_3D(N,M6,M5,M4)+heatFlx                  
                  HeatSource_col(N2,N1)         = HeatSource_col(N2,N1)-HeatFlx
                ENDIF

                if(.not.fixWaterLevel)then
                  WaterFlowSoiMicP_3D(N,M6,M5,M4)  = WaterFlowSoiMicP_3D(N,M6,M5,M4)+WaterFlow2Micpt_3D(N,M6,M5,M4)
                  WaterFlowSoiMicPX_3D(N,M6,M5,M4) = WaterFlowSoiMicPX_3D(N,M6,M5,M4)+WaterFlow2MicptX_3D(N,M6,M5,M4)
                  WaterFlowSoiMacP_3D(N,M6,M5,M4)  = WaterFlowSoiMacP_3D(N,M6,M5,M4)+WaterFlow2Macpt_3D(N,M6,M5,M4)
                  WaterFlow2MicPM_3D(M,N,M6,M5,M4) = WaterFlow2Micpt_3D(N,M6,M5,M4)
                  WaterFlow2MacPM_3D(M,N,M6,M5,M4) = WaterFlow2Macpt_3D(N,M6,M5,M4)
                endif
!                if(M6==16)then
!                write(411,*)I+J/24.,HeatFlow2Soili_3D(N,M6,M5,M4),heatFlx,M6,M6,M4
!                endif
                HeatFlow2Soil_3D(N,M6,M5,M4)     = HeatFlow2Soil_3D(N,M6,M5,M4)+HeatFlow2Soili_3D(N,M6,M5,M4)
              ENDIF
            ELSE
              WaterFlow2Micpt_3D(N,M6,M5,M4)   = 0.0_r8
              WaterFlow2MicptX_3D(N,M6,M5,M4)  = 0.0_r8
              WaterFlow2Macpt_3D(N,M6,M5,M4)   = 0.0_r8
              HeatFlow2Soili_3D(N,M6,M5,M4)    = 0.0_r8

              WaterFlow2MicPM_3D(M,N,M6,M5,M4) = 0.0_r8
              WaterFlow2MacPM_3D(M,N,M6,M5,M4) = 0.0_r8
            ENDIF
          ENDDO D9575  !loop NN
    !
    !     NET WATER AND HEAT FLUXES IN RUNOFF AND SNOW DRIFT
    !
    !     TQR1,THQR1=net runoff,convective heat from runoff
    !     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
    !     WatFlx2LitRByRunoff,HeatFlx2LitRByRunoff=runoff, convective heat from runoff
    !     DrySnoFlxBySnowRedistribut,WatFlxBySnowRedistribut,IceFlxBySnowRedistribut=snow,water,ice transfer
    !     HeatFlxBySnowRedistribut=convective heat transfer from snow,water,ice transfer

          IF(L.EQ.NUM(N2,N1) .AND. N.NE.iVerticalDirection)THEN
            if(snowRedist_model)call SumSnowDriftByRunoffM(M,N,N1,N2,N4,N5,N4B,N5B)
          ENDIF
!
        !     NET WATER AND HEAT FLUXES THROUGH SOIL AND SNOWPACK
        !
        !     TFLWL,THFLWL=net water micropore,macropore flux
        !     THFLWL=net convective+conductive heat flux
        !     FLWL =micropore water,heat flux
        !     WaterFlow2Macpt_3D=macropore water,heat flux
        !     HFLWL=soil heat flux
!
          IF(FlowDirIndicator(N2,N1).NE.iVerticalDirection .OR. N.EQ.iVerticalDirection)THEN
            D1200: DO LL=N6,NL(N5,N4)
              IF(VLSoilPoreMicP_vr(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
                N6=LL
                exit
              ENDIF
            ENDDO D1200
            !compute the net flux 
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
              if(.not.fixWaterLevel)then
                TWatFlow2MicP_3DM_vr(N3,N2,N1)     = TWatFlow2MicP_3DM_vr(N3,N2,N1)+WaterFlow2Micpt_3D(N,N3,N2,N1)-WaterFlow2Micpt_3D(N,N6,N5,N4)
                TWaterFlow2MicptX_3DM_vr(N3,N2,N1) = TWaterFlow2MicptX_3DM_vr(N3,N2,N1)+WaterFlow2MicptX_3D(N,N3,N2,N1) -WaterFlow2MicptX_3D(N,N6,N5,N4)
                TWaterFlow2Macpt_3DM_vr(N3,N2,N1)  = TWaterFlow2Macpt_3DM_vr(N3,N2,N1)+WaterFlow2Macpt_3D(N,N3,N2,N1)-WaterFlow2Macpt_3D(N,N6,N5,N4)
              endif
              THeatFlow2Soil_3DM_vr(N3,N2,N1)    = THeatFlow2Soil_3DM_vr(N3,N2,N1)+HeatFlow2Soili_3D(N,N3,N2,N1)-HeatFlow2Soili_3D(N,N6,N5,N4)

            ELSE
              TWatFlow2MicP_3DM_vr(N3,N2,N1)     = 0.0_r8
              TWaterFlow2MicptX_3DM_vr(N3,N2,N1) = 0.0_r8
              TWaterFlow2Macpt_3DM_vr(N3,N2,N1)  = 0.0_r8
              THeatFlow2Soil_3DM_vr(N3,N2,N1)    = 0.0_r8
            ENDIF
          ENDIF
        ENDDO D9580  !loop boundaries through 3 directions
!
!     INFILTRATION OF WATER FROM MACROPORES INTO MICROPORES
!
!     VLWatMacP1=macropore volume
!     FINHX,FINHL=macro-micropore transfer unltd,ltd by water,air volume
!     FWatExMacP2MicPM_vr=macro-micropore transfer for use in TranspNoSalt.f
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
          VLWatMicP1X = VLWatMicP1_vr(N3,N2,N1)+TWatFlow2MicP_3DM_vr(N3,N2,N1)+FWatIrrigate2MicP1_vr(N3,N2,N1)
          VOLP1X      = AZMAX1(VLMicP1_vr(N3,N2,N1)-VLWatMicP1X-VLiceMicP1_vr(N3,N2,N1))
          VLWatMacP1X = VLWatMacP1_vr(N3,N2,N1)+TWaterFlow2Macpt_3DM_vr(N3,N2,N1)
          VOLPH1X     = AZMAX1(VLMacP1_vr(N3,N2,N1)-VLWatMacP1X-VLiceMacP1_vr(N3,N2,N1))

          IF(FINHX.GT.0.0_r8)THEN
            FWatExMacP2MicPiM_vr(N3,N2,N1)=AZMAX1(AMIN1(FINHX,VLWatMacP1X,VOLP1X))
          ELSE
            FWatExMacP2MicPiM_vr(N3,N2,N1)=AZMIN1(AMAX1(FINHX,-VOLPH1X,-VLWatMicP1X))
          ENDIF
          FWatExMacP2MicPM_vr(M,N3,N2,N1) = FWatExMacP2MicPiM_vr(N3,N2,N1)
          FWatExMacP2MicP_vr(N3,N2,N1)    = FWatExMacP2MicP_vr(N3,N2,N1)+FWatExMacP2MicPiM_vr(N3,N2,N1)
        ELSE
          FWatExMacP2MicPiM_vr(N3,N2,N1)   = 0.0_r8
          FWatExMacP2MicPM_vr(M,N3,N2,N1) = 0.0_r8
        ENDIF

        call FreezeThawIterateM(I,J,M,NY,NX,L,N1,N2,N3)
!
!     DISSIPATE WETTING FRONT
!
!     VLWatMicP1=soil micropore water content
!     VLWatMicPX1=soil micropore water content behind wetting front
!     FLWVL=water flux from wetted to drier soil
!
      ENDDO D9585  !loop vertically through layers
    ENDDO D9590
  ENDDO D9595

  end subroutine LateralWatHeatExchIterateM
!------------------------------------------------------------------------------------------

  subroutine BoundaryDrain(I,J,N,N1,N2,N3,M4,M5,M6,RainEkReducedKsat,XN,RechargSubSurf,RechargRateWTBL)    
  !
  !Drainage across the lateral or bottom boundaries.      
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: N
  integer,  intent(in) :: N1,N2,N3  !source
  integer,  intent(in) :: M4,M5,M6  !dest
  real(r8), intent(in) :: RainEkReducedKsat 
  real(r8), intent(in) :: XN        !direction
  real(r8), intent(in) :: RechargSubSurf,RechargRateWTBL
  real(r8) :: THETA1,THETAX,HydcondSrc,RechargRate
  real(r8) :: watflx,heatflx
  integer :: K1,KL
!

! IDWaterTable=water table flag

  RechargRate=RechargSubSurf*RechargRateWTBL
  !involve no water table or in vertical direction
  THETA1 = AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicP1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
  THETAX = AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicPX1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
  K1     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-THETA1)/POROS_vr(N3,N2,N1))+1))
  KL     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-THETAX)/POROS_vr(N3,N2,N1))+1))

  !surface
  IF(N3.EQ.NUM(N2,N1))THEN
    HydcondSrc=HydroCond_3D(N,K1,N3,N2,N1)*RainEkReducedKsat
  ELSE
    HydcondSrc=HydroCond_3D(N,K1,N3,N2,N1)
  ENDIF

  WaterFlow2Micpt_3D(N,M6,M5,M4)=AMIN1(VLWatMicP1_vr(N3,N2,N1)*dts_wat, &
    XN*mGravAccelerat*(-ABS(SLOPE(N,N2,N1)))*HydcondSrc*AREA(3,N3,N2,N1)) &
    *RechargRate*dts_HeatWatTP

  if(abs(WaterFlow2Micpt_3D(N,M6,M5,M4)).LT.tiny_wat)then
    WaterFlow2Micpt_3D(N,M6,M5,M4) = 0._r8
  endif

  WaterFlow2MicptX_3D(N,M6,M5,M4) = WaterFlow2Micpt_3D(N,M6,M5,M4)  
  WaterFlow2Macpt_3D(N,M6,M5,M4)  = AMIN1(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
    ,XN*mGravAccelerat*(-ABS(SLOPE(N,N2,N1)))*HydroCondMacP1_vr(N3,N2,N1)*AREA(3,N3,N2,N1)) &
    *RechargRate*dts_HeatWatTP

  if(abs(WaterFlow2Macpt_3D(N,M6,M5,M4)).LT.tiny_wat)then
    WaterFlow2Macpt_3D(N,M6,M5,M4)=0._r8
  endif  

  watflx                        = WaterFlow2Micpt_3D(N,M6,M5,M4)+WaterFlow2Macpt_3D(N,M6,M5,M4)
  heatflx                       = cpw*TKSoil1_vr(N3,N2,N1)*watflx
  HeatFlow2Soili_3D(N,M6,M5,M4) = heatflx
  
  if(N3==M6 .and. N2==M5 .and. N1==M4)then
    QDrain_col(N2,N1)    = QDrain_col(N2,N1) - watflx
    QDrainM_col(N2,N1)   = QDrainM_col(N2,N1) - watflx
    HeatDrain_col(N2,N1) = HeatDrain_col(N2,N1)-heatflx
  else
    QDrain_col(N2,N1)    = QDrain_col(N2,N1) + watflx
    QDrainM_col(N2,N1)   = QDrainM_col(N2,N1) + watflx
    HeatDrain_col(N2,N1) = HeatDrain_col(N2,N1)+heatflx
  endif
  
  end subroutine BoundaryDrain

!------------------------------------------------------------------------------------------
  subroutine UpdateSoilMoistTemp(I,J,M,NHW,NHE,NVN,NVS,twatmass0)
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: twatmass0(JY,JX)
  integer :: NY,NX,L
  real(r8) :: tk1pres,tk1l
  real(r8) :: ENGY1,VLTSoiPore,VHXX
  real(r8) :: TKXX,VLWMicPre,VLHeatCapacityPre,dHeat
  integer :: it   !current model step in the simulation year
  real(r8) :: dwat,dwat1
  real(r8),parameter :: tau=100._r8
  it=(I-1)*24+J

! begin_execution
  D11: DO NX=NHW,NHE
    D12: DO NY=NVN,NVS
      D13: DO L = NUM(NY,NX), NL(NY,NX)

        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

          VLWMicPre              = VLWatMicP1_vr(L,NY,NX)
          VLWatMicP1_vr(L,NY,NX) = VLWatMicP1_vr(L,NY,NX)+TWatFlow2MicP_3DM_vr(L,NY,NX)+FWatExMacP2MicPiM_vr(L,NY,NX) &
            +TMLiceThawedMicP_vr(L,NY,NX)+FWatIrrigate2MicP1_vr(L,NY,NX)

          VLWatMicPX1_vr(L,NY,NX) = VLWatMicPX1_vr(L,NY,NX)+TWaterFlow2MicptX_3DM_vr(L,NY,NX)+FWatExMacP2MicPiM_vr(L,NY,NX) &
             +TMLiceThawedMicP_vr(L,NY,NX)+FWatIrrigate2MicP1_vr(L,NY,NX)
          VLWatMicPX1_vr(L,NY,NX) = AMIN1(VLWatMicP1_vr(L,NY,NX),VLWatMicPX1_vr(L,NY,NX))
          VLiceMicP1_vr(L,NY,NX)  = VLiceMicP1_vr(L,NY,NX)-TMLiceThawedMicP_vr(L,NY,NX)/DENSICE

          VLWatMacP1_vr(L,NY,NX) = VLWatMacP1_vr(L,NY,NX)+TWaterFlow2Macpt_3DM_vr(L,NY,NX)-FWatExMacP2MicPiM_vr(L,NY,NX) &
            +TMLiceThawedMacP_vr(L,NY,NX)
          VLiceMacP1_vr(L,NY,NX) = VLiceMacP1_vr(L,NY,NX)-TMLiceThawedMacP_vr(L,NY,NX)/DENSICE

          IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
            ! air-filled space
            VLairMicP_vr(L,NY,NX)  = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
            VLairMicP1_vr(L,NY,NX) = AZMAX1(VLairMicP_vr(L,NY,NX))
            VLairMacP_vr(L,NY,NX)  = VLMacP1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
            VLairMacP1_vr(L,NY,NX) = AZMAX1(VLairMacP_vr(L,NY,NX))
            VLMacP1_vr(L,NY,NX)    = AZMAX1(VLMacP_vr(L,NY,NX)-FVOLAH*CCLAY_vr(L,NY,NX) &
              *(safe_adb(VLWatMicP1_vr(L,NY,NX),VLSoilMicP_vr(L,NY,NX))-WiltPoint_vr(L,NY,NX))*VGeomLayer_vr(L,NY,NX))
          ELSE
            VLairMicP_vr(L,NY,NX)  = 0.0_r8
            VLairMicP1_vr(L,NY,NX) = 0.0_r8
            VLairMacP_vr(L,NY,NX)  = 0.0_r8
            VLairMacP1_vr(L,NY,NX) = 0.0_r8
            VLMicP1_vr(L,NY,NX)    = VLWatMicP1_vr(L,NY,NX)+VLiceMicP1_vr(L,NY,NX)
            VLMacP1_vr(L,NY,NX)    = 0.0_r8
          ENDIF
!          if(I>=253)then
!            write(211,*)I+J/24.,L,VLWatMicP1_vr(L,NY,NX)+VLiceMicP1_vr(L,NY,NX)*DENSICE-(VLWatMicP_vr(L,NY,NX)+VLiceMicP_vr(L,NY,NX)*DENSICE),&
!              TWatFlow2MicP_3DM_vr(L,NY,NX),FWatExMacP2MicPiM_vr(L,NY,NX), TMLiceThawedMicP_vr(L,NY,NX),FWatIrrigate2MicP1_vr(L,NY,NX)
!          endif
          dwat=dwat+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
          VLTSoiPore                    = VLSoilMicP_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)
          FracSoiPAsWat_vr(L,NY,NX)     = AZMAX1t((VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX))/VLTSoiPore)
          FracSoiPAsIce_vr(L,NY,NX)     = AZMAX1t((VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))/VLTSoiPore)
          AirFilledSoilPore_vr(L,NY,NX) = AZMAX1t((VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX))/VLTSoiPore)
          
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
            SoilFracAsMacP1_vr(L,NY,NX) = 0.0_r8
            HydroCondMacP1_vr(L,NY,NX)  = 0.0_r8
          ENDIF
          SoilFracAsMicP_vr(L,NY,NX) = 1.0_r8-SoilFracAsMacP1_vr(L,NY,NX)
          TKXX                       = TKSoil1_vr(L,NY,NX)
          VHXX                       = VHeatCapacity1_vr(L,NY,NX)
          ENGY1                      = VHXX*TKXX
          
          VLHeatCapacityA_vr(L,NY,NX) = VHeatCapacitySoilM_vr(L,NY,NX)+cpw*VLWatMicP1_vr(L,NY,NX)+cpi*VLiceMicP1_vr(L,NY,NX)
          if(plantOM4Heat .and. VLHeatCapacityA_vr(L,NY,NX)>0._r8)then
            VLHeatCapacityA_vr(L,NY,NX)=VLHeatCapacityA_vr(L,NY,NX)+cpo*RootMassElm_vr(ielmc,L,NY,NX)
          endif
          VLHeatCapacityB_vr(L,NY,NX) = cpw*VLWatMacP1_vr(L,NY,NX)+cpi*VLiceMacP1_vr(L,NY,NX)
          VHeatCapacity1_vr(L,NY,NX)  = VLHeatCapacityA_vr(L,NY,NX)+VLHeatCapacityB_vr(L,NY,NX)

         !if(I>=253)then
         !  write(211,*)L,VHeatCapacity1_vr(L,NY,NX),THeatFlow2Soil_3DM_vr(L,NY,NX),HeatIrrigation1_vr(L,NY,NX), &
         !     TLPhaseChangeHeat2Soi1_vr(L,NY,NX)
         !endif

         IF(VHeatCapacity1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            
            TKSoil1_vr(L,NY,NX) = (ENGY1+THeatFlow2Soil_3DM_vr(L,NY,NX)+HeatIrrigation1_vr(L,NY,NX)&
              +TLPhaseChangeHeat2Soi1_vr(L,NY,NX))/VHeatCapacity1_vr(L,NY,NX)
!            if(L>=14)then
!            write(300+L,*)L,TKSoil1_vr(L,NY,NX),TKXX,THeatFlow2Soil_3DM_vr(L,NY,NX),HeatIrrigation1_vr(L,NY,NX), &
!              TLPhaseChangeHeat2Soi1_vr(L,NY,NX),(THeatFlow2Soil_3DM_vr(L,NY,NX)+HeatIrrigation1_vr(L,NY,NX) &
!              +TLPhaseChangeHeat2Soi1_vr(L,NY,NX))/VHeatCapacity1_vr(L,NY,NX),TKSD(NY,NX),SoilHeatSrcDepth_col(NY,NX)
!            endif
            if(TKSoil1_vr(L,NY,NX)>400._r8.or.TKSoil1_vr(L,NY,NX)<100._r8)then
              write(*,*)'======'
              write(*,*)'VLWatMicP1_vr(L,NY,NX)=',VLWMicPre,VLWatMicP1_vr(L,NY,NX),L,VLairMicP_vr(L,NY,NX),&
                VLiceMicP1_vr(L,NY,NX)
              write(*,*)TWatFlow2MicP_3DM_vr(L,NY,NX),FWatExMacP2MicPiM_vr(L,NY,NX), &
                TMLiceThawedMicP_vr(L,NY,NX),FWatIrrigate2MicP1_vr(L,NY,NX)
              write(*,*)'M, L=',M,L,NY,NX,NUM(NY,NX),VGeomLayer_vr(L,NY,NX),ZEROS2(NY,NX),VLWatMicPX1_vr(L,NY,NX)
              write(*,*)'SoilBulkDensity_vr(L,NY,NX)=',SoilBulkDensity_vr(L,NY,NX),TKS_vr(L,NY,NX),ZEROS(NY,NX)
              write(*,*)'VHeatCapacity1_vr(L,NY,NX),TKSoil1_vr(L,NY,NX),TKXX',L,VHeatCapacity1_vr(L,NY,NX),TKSoil1_vr(L,NY,NX),TKXX
              write(*,*)VLMicP1_vr(L,NY,NX),VLMacP1_vr(L,NY,NX)
              if(TKSoil1_vr(L,NY,NX)>1.e3_r8.or.TKSoil1_vr(L,NY,NX)<0._r8)call endrun(trim(mod_filename)//' at line',__LINE__)
            endif

            ! do soil warming
            if(do_warming)then
              if(is_warming_layerL(L,NY,NX))then     
                tk1l                   = TKSoil1_vr(L,NY,NX)
                dHeat                  = (TKS_ref_vr(it,L,NY,NX)-tk1l)*VHeatCapacity1_vr(L,NY,NX)*dts_HeatWatTP
                HeatSource_vr(L,NY,NX) = HeatSource_vr(L,NY,NX)+(TKS_ref_vr(it,L,NY,NX)-tk1l)*VHeatCapacity1_vr(L,NY,NX)*dts_HeatWatTP
                TKSoil1_vr(L,NY,NX)    = tk1l+dHeat/VHeatCapacity1_vr(L,NY,NX)
              endif
            endif
            
            if(TKSoil1_vr(L,NY,NX)>400.0_r8)then
              write(*,*)'TKSoil1_vr(L,NY,NX)',I+J/24.,M,L,TKSoil1_vr(L,NY,NX),tk1l,VHXX,VHeatCapacity1_vr(L,NY,NX),TKS_vr(L,NY,NX)
              write(*,*)ENGY1/VHeatCapacity1_vr(L,NY,NX),THeatFlow2Soil_3DM_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX),&
               TLPhaseChangeHeat2Soi1_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX),HeatIrrigation1_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX)
              write(*,*)VLWatMicPM_vr(M,L,NY,NX),VLWatMicPM_vr(M+1,L,NY,NX),THeatFlow2Soil_3DM_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX)
              call endrun('crazy temperature in '//trim(mod_filename)//' at line',__LINE__)
            endif
          ELSEIF(L.EQ.1)THEN
            TKSoil1_vr(L,NY,NX)=TairK_col(NY,NX)
          ELSE
            TKSoil1_vr(L,NY,NX)=TKSoil1_vr(L-1,NY,NX)
          ENDIF
        ENDIF
      ENDDO D13
!      if(I>=317 .and. J>=18)then
!        if(NX==1)then
!          write(401,*)I+J/24.,NY,NX,M,'xdwat',twatmass0(NY,NX),dwat
!          write(401,*)I+J/24.,NY,NX,M,'wat',(VLWatMicP1_vr(L,NY,NX),L=NUM(NY,NX), NL(NY,NX))
!        else
!          write(402,*)I+J/24.,NY,NX,M,'xdwat',twatmass0(NY,NX),dwat
!        endif  
!      endif
    ENDDO D12
  ENDDO D11

  end subroutine UpdateSoilMoistTemp

!------------------------------------------------------------------------------------------

  subroutine UpdateStateFluxAtM(I,J,M,NHW,NHE,NVN,NVS)
  !
  !Description
  ! Early exit of the watsub solver
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX
  integer :: L,NUX,LL,Ls
  real(r8):: VHCPX,HFLX,dCPX,dwat

! begin_execution
  D9795: DO NX=NHW,NHE
    D9790: DO NY=NVN,NVS
      dwat=0._r8
      !
      ! SOIL LAYER WATER, ICE AND TEMPERATURE
      !
      D9785: DO L=NUM(NY,NX),NL(NY,NX)
        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

          !record intermediate variables for bgc calculation
          VLWatMicPM_vr(M+1,L,NY,NX) = VLWatMicP1_vr(L,NY,NX)
          VLWatMacPM_vr(M+1,L,NY,NX) = VLWatMacP1_vr(L,NY,NX)
          VLsoiAirPM_vr(M+1,L,NY,NX) = VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX)+THETPI*(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))

          !change in soil air volume
          ReductVLsoiAirPM_vr(M,L,NY,NX) = VLsoiAirPM_vr(M,L,NY,NX)-VLsoiAirPM_vr(M+1,L,NY,NX)
          AirFilledSoilPoreM_vr(M+1,L,NY,NX)         = AirFilledSoilPore_vr(L,NY,NX)
 !
        ELSE
          !layer L disappears
          VLWatMicPM_vr(M+1,L,NY,NX)  = 0.0_r8
          VLWatMacPM_vr(M+1,L,NY,NX)  = 0.0_r8
          VLsoiAirPM_vr(M+1,L,NY,NX)  = 0.0_r8
          ReductVLsoiAirPM_vr(M,L,NY,NX) = VLsoiAirPM_vr(M,L,NY,NX)
          AirFilledSoilPoreM_vr(M+1,L,NY,NX)         = 0.0_r8
        ENDIF
        dwat=dwat+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
      ENDDO D9785

      !
      !       RESET SURFACE LAYER NUMBER AND TRANSFER ALL WATER TO SOIL SURFACE LAYER
      !       IF LAKE SURFACE LAYER IS LOST TO EVAPORATION
      !
      !       NUM=new surface layer number after CO2CompenPoint_nodeete lake evaporation
      !       LakeSurfFlowMicP,LakeSurfFlowMacP,LakeSurfHeatFlux=lake surface water flux, heat flux if lake surface disappears
      ! NUM layer is pure water, and its heat capacity and thereby water is too little
      IF(SoilBulkDensity_vr(NUM(NY,NX),NY,NX).LE.ZERO .AND. VHeatCapacity1_vr(NUM(NY,NX),NY,NX).LE.VHCPNX_col(NY,NX))THEN
        !the soil/water profile moves down        
        NUX=NUM(NY,NX)        
!        if(I>=317 .and. J>=18)then
!          if(NX==1)then
!            write(401,*)'dwat0',dwat
!            write(401,*)'nux0',(VLWatMicP1_vr(L,NY,NX),L=NUX,NL(NY,NX))
!            write(401,*)'nux00',(VLWatMacP1_vr(L,NY,NX),L=NUX,NL(NY,NX))
!            write(401,*)'nuxii',((VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX)),L=NUX,NL(NY,NX))
!          endif  
!        endif             
        D9970: DO  LL=NUX+1,NL(NY,NX)
          IF(VLSoilPoreMicP_vr(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
            NUM(NY,NX)                     = LL
            dCPX = cpw*(VLWatMicP1_vr(NUX,NY,NX)+VLWatMacP1_vr(NUX,NY,NX))+cpi*(VLiceMicP1_vr(NUX,NY,NX))
            HFLX = dCPX*TKSoil1_vr(NUX,NY,NX)
            
            WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX)  = WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX) +VLWatMicP1_vr(NUX,NY,NX)
            WaterFlowSoiMicPX_3D(3,NUM(NY,NX),NY,NX) = WaterFlowSoiMicPX_3D(3,NUM(NY,NX),NY,NX)+VLWatMicP1_vr(NUX,NY,NX)
            WaterFlowSoiMacP_3D(3,NUM(NY,NX),NY,NX)  = WaterFlowSoiMacP_3D(3,NUM(NY,NX),NY,NX) +VLWatMacP1_vr(NUX,NY,NX)

            H2OFlow2TopSoiMicP_col(NY,NX)  = WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX) 
            H2OFlow2TopSoiMicPX_col(NY,NX) = WaterFlowSoiMicPX_3D(3,NUM(NY,NX),NY,NX)
            H2OFlow2TopSoiMacP_col(NY,NX)  = WaterFlowSoiMacP_3D(3,NUM(NY,NX),NY,NX) 
            
            HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX) = HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX) + HFLX
            HeatFlow2TopSoi_col(NY,NX)           = HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX)

            QIceInflx_vr(NUM(NY,NX),NY,NX)      = QIceInflx_vr(NUM(NY,NX),NY,NX)+VLiceMicP1_vr(NUX,NY,NX)
            VLWatMicP1_vr(NUM(NY,NX),NY,NX)     = VLWatMicP1_vr(NUM(NY,NX),NY,NX)+ VLWatMicP1_vr(NUX,NY,NX)
            VLWatMacP1_vr(NUM(NY,NX),NY,NX)     = VLWatMacP1_vr(NUM(NY,NX),NY,NX)+ VLWatMacP1_vr(NUX,NY,NX)
            VLiceMicP1_vr(NUM(NY,NX),NY,NX)     = VLiceMicP1_vr(NUM(NY,NX),NY,NX)+ VLiceMicP1_vr(NUX,NY,NX)
            VHCPX                               = VHeatCapacity1_vr(NUM(NY,NX),NY,NX)
            VHeatCapacity1_vr(NUM(NY,NX),NY,NX) = VHeatCapacity1_vr(NUM(NY,NX),NY,NX)+dCPX
            TKSoil1_vr(NUM(NY,NX),NY,NX)        = (VHCPX*TKSoil1_vr(NUM(NY,NX),NY,NX)+HFLX)/VHeatCapacity1_vr(NUM(NY,NX),NY,NX)
            exit
          ENDIF
        ENDDO D9970
!        if(I>=317 .and. J>=18)then
!          write(214,*)'nox',I+J/24.,M,VLWatMicP1_vr(NUX:NUM(NY,NX)-1,NY,NX)
!          write(214,*)'NUX',M,NUX,NUM(NY,NX),NU(NY,NX),VLWatMicP1_vr(NUX,NY,NX),&
!          VLWatMacP1_vr(NUX,NY,NX),WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX),TKSoil1_vr(NUX,NY,NX)
!          if(NX==1)then
!            write(401,*)'nux1',(VLWatMicP1_vr(L,NY,NX),L=NUX,NL(NY,NX))
!            write(401,*)'nux2',(VLiceMacP1_vr(L,NY,NX),L=NUX,NL(NY,NX))
!          endif  
!        endif  
        VLWatMicP1_vr(NUX,NY,NX) = 0._r8
        VLWatMacP1_vr(NUX,NY,NX) = 0._r8
        VLiceMicP1_vr(NUX,NY,NX) = 0._r8
        VLiceMacP1_vr(NUX,NY,NX) = 0._r8
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
  if(fixWaterLevel)return

  D9695: DO NX=NHW,NHE
    D9690: DO NY=NVN,NVS
      !there is not change in surface layer      
      IF(NUM(NY,NX).EQ.NU(NY,NX))THEN
        !vertical flux, incoming
        !top layer is the same
        LakeSurfFlowMicP_col(NY,NX)  = WaterFlowSoiMicP_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfFlowMicPX_col(NY,NX) = WaterFlowSoiMicPX_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfFlowMacP_col(NY,NX)  = WaterFlowSoiMacP_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfHeatFlux_col(NY,NX)  = HeatFlow2Soil_3D(3,N6X(NY,NX),NY,NX)
      !the top soil/water layer has changed, ponding water
      ELSE
        LakeSurfFlowMicP_col(NY,NX)  = H2OFlow2TopSoiMicP_col(NY,NX)
        LakeSurfFlowMicPX_col(NY,NX) = H2OFlow2TopSoiMicPX_col(NY,NX)
        LakeSurfFlowMacP_col(NY,NX)  = H2OFlow2TopSoiMacP_col(NY,NX)
        LakeSurfHeatFlux_col(NY,NX)  = HeatFlow2TopSoi_col(NY,NX)
      ENDIF
!      call writeSurfDiagnosis(I,J,NY,NX)
    ENDDO D9690
  ENDDO D9695
  
  end subroutine UpdateFluxAtExit

!--------------------------------------------------------------------------
  function get_DPTHH(L,NY,NX)result(DPTHH)
  !
  !get the depth to layer macropore water
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8) :: DPTHH

  IF(VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    DPTHH= CumDepz2LayBottom_vr(L,NY,NX)-(VLWatMacP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))/VLMacP1_vr(L,NY,NX)*DLYR_3D(3,L,NY,NX)
  ELSE
    DPTHH=CumDepz2LayBottom_vr(L,NY,NX)
  ENDIF

  end function get_DPTHH

!--------------------------------------------------------------------------
  subroutine Config4TileDrainage(L,NY,NX,DPTHH,DoMicPDischarg2Tile,DoMacPDischarg2Tile)

  implicit none
  integer , intent(in) :: L,NY,NX
  real(r8),intent(in) :: DPTHH   !depth to layer macropore water
  logical, intent(out) :: DoMicPDischarg2Tile
  logical, intent(out) :: DoMacPDischarg2Tile

  integer :: LL
  real(r8) :: TileWaterTableX
!
!     IDWaterTable=water table flag
!     DPTH,TileWaterTable_col=depth to layer midpoint, artificial water table
!     PSISM1,PSISE=soil,air entry matric potential
!     TileWaterTableX=equilibrium water potential with artificial water table
!     DoMicPDischarg2Tile=micropore discharge flag to artificial water table
! 
  IF(IDWaterTable_col(NY,NX).GE.3 .AND. SoilDepthMidLay_vr(L,NY,NX).LT.TileWaterTable_col(NY,NX))THEN
    IF(PSISM1_vr(L,NY,NX).GT.mGravAccelerat*(SoilDepthMidLay_vr(L,NY,NX)-TileWaterTable_col(NY,NX)))THEN
      DoMicPDischarg2Tile=.true.
      IF(L.LT.NL(NY,NX))THEN
        D9568: DO  LL=L+1,NL(NY,NX)
          TileWaterTableX=TileWaterTable_col(NY,NX)+PSISE_vr(LL,NY,NX)/mGravAccelerat
          IF(SoilDepthMidLay_vr(LL,NY,NX).LT.TileWaterTableX)THEN
            IF((PSISM1_vr(LL,NY,NX).LE.mGravAccelerat*(SoilDepthMidLay_vr(LL,NY,NX)-TileWaterTableX) &
              .AND.L.NE.NL(NY,NX)).OR.SoilDepthMidLay_vr(LL,NY,NX).GT.ActiveLayDepZ_col(NY,NX))THEN
              DoMicPDischarg2Tile=.false.
            ENDIF
          ENDIF
        ENDDO D9568
      ENDIF
    ELSE
      DoMicPDischarg2Tile=.false.
    ENDIF
  ELSE
    DoMicPDischarg2Tile=.false.
  ENDIF
!
!     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO TILE DRAIN
!
!     VLMacP1,VLWatMacP1,VLiceMacP1=macropore volume,water,ice content
!     CDPTH=depth to layer bottom
!     DLYR=layer thickness
!     DoMacPDischarg2Tile=macropore discharge flag to artificial water table
!
  IF(IDWaterTable_col(NY,NX).GE.3 .AND. DPTHH.LT.TileWaterTable_col(NY,NX) .AND. VLWatMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
! artificial water table, e.g. tile drainage
    DoMacPDischarg2Tile=.true.
    IF(L.LT.NL(NY,NX))THEN
      D9569: DO  LL=L+1,NL(NY,NX)
        IF(SoilDepthMidLay_vr(LL,NY,NX).LT.TileWaterTable_col(NY,NX))THEN
          ! the layer is above tile drain water table
          IF(VLMacP1_vr(LL,NY,NX).LE.ZEROS(NY,NX))THEN
            ! no macropore flow
            DoMacPDischarg2Tile=.false.
          ENDIF
        ENDIF
      ENDDO D9569
    ENDIF
  ELSE
    DoMacPDischarg2Tile=.false.
  ENDIF
  end subroutine Config4TileDrainage
!------------------------------------------------------------------

  subroutine Config4WaterTableDrain(L,NY,NX,DPTHH,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL)

  implicit none
  integer , intent(in) :: L,NY,NX
  real(r8), intent(in):: DPTHH  !depth to layer macropore water
  logical , intent(out):: DoMicPDischarg2ExtWTBL
  logical , intent(out):: DoMacPDischarg2ExtWTBL

  real(r8) :: ExtWaterTableEquil
  integer :: LL

!     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO WATER TABLE
!
!     IDWaterTable=water table flag
!     DPTH,ExtWaterTable=depth to layer midpoint,natural water table
!     PSISM1_vr(<0),PSISE_vr(<0)=matric,air entry water potential [MPa], from water height to MPa, 1.e-6*9.8*1.e3*H
!     ExtWaterTableEquil=equilibrium water potential with natural water table
!     ActiveLayDepZ_col=active layer depth
!     DoMicPDischarg2ExtWTBL=micropore discharge flag to natural water table
!  total water potential psi+rho*grav*h

  IF(IDWaterTable_col(NY,NX).NE.0 .AND. SoilDepthMidLay_vr(L,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN
    !the layer mid-depth is lower than water table
    IF(PSISM1_vr(L,NY,NX).GT.mGravAccelerat*(SoilDepthMidLay_vr(L,NY,NX)-ExtWaterTable_col(NY,NX)))THEN
      DoMicPDischarg2ExtWTBL=.true.
      D9565: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
        !water level difference
        ExtWaterTableEquil=ExtWaterTable_col(NY,NX)+PSISE_vr(LL,NY,NX)/mGravAccelerat  

        IF(SoilDepthMidLay_vr(LL,NY,NX).LT.ExtWaterTableEquil)THEN
          IF((PSISM1_vr(LL,NY,NX).LE.mGravAccelerat*(SoilDepthMidLay_vr(LL,NY,NX)-ExtWaterTableEquil) &
            .AND. L.NE.NL(NY,NX)) .OR. SoilDepthMidLay_vr(LL,NY,NX).GT.ActiveLayDepZ_col(NY,NX))THEN
            DoMicPDischarg2ExtWTBL=.false.
          ENDIF
        ENDIF
      ENDDO D9565
    ELSE
      DoMicPDischarg2ExtWTBL=.false.
    ENDIF
  ELSE
    DoMicPDischarg2ExtWTBL=.false.
  ENDIF
!
!     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO WATER TABLE
!
!     VLMacP1,VLWatMacP1,VLiceMacP1=macropore volume,water,ice content
!     CDPTH=depth to layer bottom
!     DLYR=layer thickness
!     DoMacPDischarg2ExtWTBL=macropore discharge flag to natural water table
!  
  IF(IDWaterTable_col(NY,NX).NE.0 .AND. DPTHH.LT.ExtWaterTable_col(NY,NX) &
    .AND. VLWatMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    !active water table
    DoMacPDischarg2ExtWTBL=.true.
!     DO 9566 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
!     IF(SoilDepthMidLay_vr(LL,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN
!     IF(VLMacP1_vr(LL,NY,NX).LE.ZEROS(NY,NX))THEN
!     DoMacPDischarg2ExtWTBL=1
!     ENDIF
!     ENDIFx
!9566  CONTINUE
  ELSE
    DoMacPDischarg2ExtWTBL=.false.
  ENDIF
  end subroutine Config4WaterTableDrain

!------------------------------------------------------------------------------------------

  subroutine InitSoil3DModelIterateM(I,J,M,NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX,L
  real(r8) :: THETWT,TScal4Aquadifsvity,THETWA
  real(r8) :: VLSoiPorAvail   !avaiable pore space for air + water, excluding ice
  real(r8) :: VLWatSoi,scalar
  real(r8) :: THETWH,Z3S

  WaterFlow2Micptl_3D=0._r8  
  WaterFlow2Macptl_3D=0._r8        

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      D9885: DO L=NUM(NY,NX),NL(NY,NX)
        TMLiceThawedMicP_vr(L,NY,NX)       = 0.0_r8
        TMLiceThawedMacP_vr(L,NY,NX)       = 0.0_r8
        TLPhaseChangeHeat2Soi1_vr(L,NY,NX) = 0.0_r8
        TWatFlow2MicP_3DM_vr(L,NY,NX)      = 0.0_r8
        TWaterFlow2MicptX_3DM_vr(L,NY,NX)  = 0.0_r8
        TWaterFlow2Macpt_3DM_vr(L,NY,NX)   = 0.0_r8
        THeatFlow2Soil_3DM_vr(L,NY,NX)     = 0.0_r8
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
    !   TORT,TortMacPM_vr=tortuosity for aqueous diffn in micropores,macropres
    !
        VLWatSoi      = VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)
        VLSoiPorAvail = VLMicP1_vr(L,NY,NX)+VLMacP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)

        !unsaturated soil
        IF(VLSoiPorAvail.GT.ZEROS2(NY,NX) .AND. VLsoiAirPM_vr(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
          THETWA                         = AZMAX1(AMIN1(1.0_r8,VLWatSoi/VLSoiPorAvail))
          TScal4Aquadifsvity             = TEFAQUDIF(TKSoil1_vr(L,NY,NX))
          Z3S                            = FieldCapacity_vr(L,NY,NX)/POROS_vr(L,NY,NX)
          scalar                         = TScal4Aquadifsvity*XNPD
          DiffusivitySolutEffM_vr(M,L,NY,NX) = fDiffusivitySolutEff(scalar,THETWA,Z3S)          
        ELSE
          DiffusivitySolutEffM_vr(M,L,NY,NX)=0.0_r8
        ENDIF

        !Pure water layer
        IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
          THETWT                  = safe_adb(VLWatMicPM_vr(M,L,NY,NX),VLSoilMicP_vr(L,NY,NX))
          TortMicPM_vr(M,L,NY,NX) = TortMicporew(THETWT)*(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
        !soil layer  
        ELSE
          TortMicPM_vr(M,L,NY,NX)=0.7_r8
        ENDIF
        !
        IF(VLMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          THETWH                  = VLWatMacPM_vr(M,L,NY,NX)/VLMacP1_vr(L,NY,NX)
          TortMacPM_vr(M,L,NY,NX) = TortMacporew(THETWH)*SoilFracAsMacP_vr(L,NY,NX)
        ELSE
          TortMacPM_vr(M,L,NY,NX)=0.0_r8
        ENDIF
      ENDDO D9885
    ENDDO
  ENDDO
  end subroutine InitSoil3DModelIterateM

!------------------------------------------------------------------------------------------
  subroutine MicporeDarcyFlow(NY,NX,N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,KSatRedusByRainKinetEnergy,&
    PredDarcyFlowMax,WatDarcyFlowMicP,HeatByDarcyFlowMicP,PSISV1,PSISVL)          
  implicit none
  integer, intent(in)  :: NY,NX,N
  integer, intent(in)  :: N1,N2,N3  !source grid
  integer, intent(in)  :: N4,N5,N6  !destination grid
  real(r8), intent(in) :: THETA1,THETAL,KSatRedusByRainKinetEnergy
  real(r8), intent(out):: PredDarcyFlowMax
  real(r8), intent(out):: WatDarcyFlowMicP
  real(r8), intent(out):: HeatByDarcyFlowMicP
  real(r8), intent(out):: PSISV1,PSISVL
  real(r8) :: AVE_CONDUCTANCE,HydCondSrc,HydCondDest
  real(r8) :: PSISTL,PSIST1,THETW1
  real(r8) :: PtWatDarcyFlux,PredDarcyFlow,PSISTMP,THETWL
  integer :: K1,KL

   PredDarcyFlowMax    = 0._r8
   WatDarcyFlowMicP    = 0._r8
   HeatByDarcyFlowMicP = 0._r8
   PSISV1              = 0._r8
   PSISVL              = 0._r8

  IF(DLYR_3D(N,N6,N5,N4).LE.0._r8 .OR. DLYR_3D(N,N3,N2,N1).LE.0._r8 .OR. fixWaterLevel)return

  !     THETW1,THETWL=water concentrations in source,destination cells
  !air entry pressure is the pressure when large-pores started to be air-filled, when soil matric pressure is less
  !than air entry pressure, it is unsaturated
  IF(PSISoilMatricPtmp_vr(N3,N2,N1).GT.PSISoilAirEntry(N3,N2,N1) .AND. &
    PSISoilMatricPtmp_vr(N6,N5,N4).GT.PSISoilAirEntry(N6,N5,N4))THEN
!    print*,'both are saturated',POROS_vr(N3,N2,N1),POROS_vr(N6,N5,N4)
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
!    print*,'source grid is saturated'
    THETW1              = THETA1
    THETWL              = AMAX1(SoilWatAirDry_vr(N6,N5,N4),AMIN1(POROS_vr(N6,N5,N4),safe_adb(VLWatMicPX1_vr(N6,N5,N4),VLSoilMicP_vr(N6,N5,N4))))
    K1                  = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-THETW1)/POROS_vr(N3,N2,N1))+1))
    KL                  = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N6,N5,N4)-AMIN1(ThetaSat_vr(N6,N5,N4),THETWL))/POROS_vr(N6,N5,N4))+1))
    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    
    IF(VLSoilMicPMass_vr(N6,N5,N4).GT.ZEROS(NY,NX))THEN
!      print*,'dest layer is pore active'
      IF(THETWL.LT.FieldCapacity_vr(N6,N5,N4))THEN
        !less than field capacity
        PSISM1_vr(N6,N5,N4)=AMAX1(PSIHY,-EXP(LOGPSIFLD(N5,N4)+((LOGFldCapacity_vr(N6,N5,N4)-LOG(THETWL)) &
          /FCD_vr(N6,N5,N4)*LOGPSIMND(N5,N4))))
      ELSEIF(THETWL.LT.POROS_vr(N6,N5,N4)-DTHETW)THEN
        PSISM1_vr(N6,N5,N4)=-EXP(LOGPSIAtSat(N5,N4)+(((LOGPOROS_vr(N6,N5,N4)-LOG(THETWL)) &
          /PSD_vr(N6,N5,N4))**SRP_vr(N6,N5,N4)*LOGPSIMXD(N5,N4)))
      ELSE
        !saturated
        THETWL              = POROS_vr(N6,N5,N4)
        PSISM1_vr(N6,N5,N4) = PSISE_vr(N6,N5,N4)
      ENDIF
    ELSE
      THETWL              = POROS_vr(N6,N5,N4)
      PSISM1_vr(N6,N5,N4) = PSISE_vr(N6,N5,N4)
    ENDIF
    !
    !     GREEN-AMPT FLOW IF ADJACENT CELL SATURATED
!
  ELSEIF(PSISoilMatricPtmp_vr(N6,N5,N4).GT.PSISoilAirEntry(N6,N5,N4))THEN
!    print*,'source grid is saturated, dest grid is not'
    THETW1 = AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicPX1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
    THETWL = THETAL
    K1     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N3,N2,N1)-AMIN1(ThetaSat_vr(N3,N2,N1),THETW1))/POROS_vr(N3,N2,N1))+1))
    KL     = MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(N6,N5,N4)-THETWL)/POROS_vr(N6,N5,N4))+1))

    IF(VLSoilMicPMass_vr(N3,N2,N1).GT.ZEROS(NY,NX))THEN
      IF(THETW1.LT.FieldCapacity_vr(N3,N2,N1))THEN
        PSISTMP             = -EXP(LOGPSIFLD(N2,N1)+(LOGFldCapacity_vr(N3,N2,N1)-LOG(THETW1)) &
          /FCD_vr(N3,N2,N1)*LOGPSIMND(N2,N1))
        PSISM1_vr(N3,N2,N1) = AMAX1(PSIHY,PSISTMP)
      ELSEIF(THETW1.LT.POROS_vr(N3,N2,N1)-DTHETW)THEN
        PSISM1_vr(N3,N2,N1)=-EXP(LOGPSIAtSat(N2,N1)+(((LOGPOROS_vr(N3,N2,N1)-LOG(THETW1)) &
          /PSD_vr(N3,N2,N1))**SRP_vr(N3,N2,N1)*LOGPSIMXD(N2,N1)))
      ELSE
        THETW1              = POROS_vr(N3,N2,N1)
        PSISM1_vr(N3,N2,N1) = PSISE_vr(N3,N2,N1)
      ENDIF
    ELSE
      THETW1              = POROS_vr(N3,N2,N1)
      PSISM1_vr(N3,N2,N1) = PSISE_vr(N3,N2,N1)
    ENDIF
    !
    !     RICHARDS FLOW IF NEITHER CELL IS SATURATED
    !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POTENTIAL)
    !
  ELSE
!    print*,'1508xxx'
    THETW1              = THETA1
    THETWL              = THETAL
    K1                  = MAX(1,MIN(100,INT(100.0*(POROS_vr(N3,N2,N1)-THETW1)/POROS_vr(N3,N2,N1))+1))
    KL                  = MAX(1,MIN(100,INT(100.0*(POROS_vr(N6,N5,N4)-THETWL)/POROS_vr(N6,N5,N4))+1))
    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    PSISM1_vr(N6,N5,N4) = PSISoilMatricPtmp_vr(N6,N5,N4)
  ENDIF

  !
!  print*,'HYDRAULIC CONUCTIVITY'
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
  PSIST1 = PSISM1_vr(N3,N2,N1)+PSIGrav_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  PSISV1 = PSISM1_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  !dest
  PSISTL = PSISM1_vr(N6,N5,N4)+PSIGrav_vr(N6,N5,N4)+PSISoilOsmotic_vr(N6,N5,N4)
  PSISVL = PSISM1_vr(N6,N5,N4)+PSISoilOsmotic_vr(N6,N5,N4)
  !
!  print*,'     HYDRAULIC CONDUCTIVITY FROM CURRENT WATER CONTENT'
  !     AND LOOKUP ARRAY GENERATED IN 'HOUR1'
  !
  !     HydCondSrc,HydCondDest=hydraulic conductivities in source,destination cells
  !     KSatRedusByRainKinetEnergy=reduction in soil surface Ksat from rainfall energy impact
  !     AVE_CONDUCTANCE=source-destination hydraulic conductance
  !     DLYR=layer thickness
  !
!  print*,HydCondSrc,DLYR_3D(N,N6,N5,N4),HydCondDest,DLYR_3D(N,N3,N2,N1)
  IF(HydCondSrc.GT.ZERO.AND.HydCondDest.GT.ZERO)THEN
    AVE_CONDUCTANCE=2.0_r8*HydCondSrc*HydCondDest/(HydCondSrc*DLYR_3D(N,N6,N5,N4)+HydCondDest*DLYR_3D(N,N3,N2,N1))
  ELSE
    AVE_CONDUCTANCE=0.0_r8
  ENDIF
  !
  !     WATER FLUX FROM WATER POTENTIALS, HYDRAULIC CONDUCTIVITY
  !     CONSTRAINED BY WATER POTENTIAL GRADIENT, COUPLED WITH
  !     CONVECTIVE HEAT FLUX FROM WATER FLUX
  !
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

    WatDarcyFlowMicP = AZMAX1(AMIN1(PredDarcyFlow,VLWatMicP2_vr(N3,N2,N1)*dts_wat,VLairMicP1_vr(N6,N5,N4)*dts_wat))
    PredDarcyFlowMax = AZMAX1(AMIN1(PtWatDarcyFlux,VLWatMicP2_vr(N3,N2,N1)*dts_wat,VLairMicP1_vr(N6,N5,N4)*dts_wat))
    ! FLQL1            = (THETW1-ThetaSat_vr(N3,N2,N1))*VLSoilMicP_vr(N3,N2,N1)
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
    HeatByDarcyFlowMicP=cpw*TKSoil1_vr(N3,N2,N1)*WatDarcyFlowMicP*HeatAdv_scal
  ELSE
    HeatByDarcyFlowMicP=cpw*TKSoil1_vr(N6,N5,N4)*WatDarcyFlowMicP*HeatAdv_scal
  ENDIF

  VLWatMicP2_vr(N3,N2,N1)         = VLWatMicP2_vr(N3,N2,N1)-WatDarcyFlowMicP
  VLWatMicP2_vr(N6,N5,N4)         = VLWatMicP2_vr(N6,N5,N4)+WatDarcyFlowMicP

  end subroutine MicporeDarcyFlow
!------------------------------------------------------------------------------------------
  subroutine MacporeFLow(NY,NX,M,N,N1,N2,N3,N4,N5,N6,WaterMacpFlow,HeatByFlowMacP,LInvalidMacP)
  implicit none
  integer, intent(in) :: NY,NX,M,N
  integer, intent(in) :: N1,N2,N3   !source grid
  integer, intent(in) :: N4,N5,N6   !destination grid
  real(r8), intent(out) :: HeatByFlowMacP   ! >0, flux into dst, 
  real(r8), intent(out) :: WaterMacpFlow
  logical, intent(inout) :: LInvalidMacP   !==1, dest grid invalid

  real(r8) :: FLWHX,PSISH1,PSISHL
  !     PSISH1,PSISHL=macropore total water potl in source,destination
  !     DLYR=layer thickness
  !     VLWatMacP1,VOLPH1=macropore water,air content

  HeatByFlowMacP = 0._r8
  WaterMacpFlow  = 0._r8
  LInvalidMacP   = .true.
  if(fixWaterLevel)return
  IF(VLMacP1_vr(N3,N2,N1).GT.ZEROS2(N2,N1) .AND. VLMacP1_vr(N6,N5,N4).GT.ZEROS2(N5,N4) .AND. (.not.LInvalidMacP))THEN
    PSISH1=PSIGrav_vr(N3,N2,N1)+mGravAccelerat*DLYR_3D(3,N3,N2,N1)*(AMIN1(1.0_r8,&
      AZMAX1(VLWatMacP1_vr(N3,N2,N1)/VLMacP1_vr(N3,N2,N1)))-0.5_r8)
    PSISHL=PSIGrav_vr(N6,N5,N4)+mGravAccelerat*DLYR_3D(3,N6,N5,N4)*(AMIN1(1.0_r8,&
      AZMAX1(VLWatMacP1_vr(N6,N5,N4)/VLMacP1_vr(N6,N5,N4)))-0.5_r8)
    !
    !     MACROPORE FLOW IF GRAVITATIONAL GRADIENT IS POSITIVE
    !     AND MACROPORE POROSITY EXISTS IN ADJACENT CELL
    !
    !     FLWHX,WaterFlow2Macpt_3D=macropore water flux unltd,ltd by source water
    !     dts_HeatWatTP=time step of flux calculations
    !     VOLW2,VOLP1=water,air contents of source,destination micropores
    !     HeatByFlowMacP=convective heat flux from micropore water flux
    !
    FLWHX=AVCNHL_3D(N,N6,N5,N4)*(PSISH1-PSISHL)*AREA(N,N3,N2,N1)*dts_HeatWatTP
    IF(N.NE.iVerticalDirection)THEN
      !horizontal direction
      IF(PSISH1.GT.PSISHL)THEN
        !flow from source to dest grid
        WaterMacpFlow=AZMAX1(AMIN1(AMIN1(VLWatMacP1_vr(N3,N2,N1),VLairMacP1_vr(N6,N5,N4))*dts_wat,FLWHX))
      ELSEIF(PSISH1.LT.PSISHL)THEN
        !flow from dest grid to source
        WaterMacpFlow=AZMIN1(AMAX1(AMAX1(-VLWatMacP1_vr(N6,N5,N4),-VLairMacP1_vr(N3,N2,N1))*dts_wat,FLWHX))
      ELSE
        WaterMacpFlow=0.0_r8
      ENDIF
    ELSE
      !vertical direction
      WaterMacpFlow=AZMAX1(AMIN1(AMIN1(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
        +WaterFlow2Macpt_3D(N,N3,N2,N1),VLairMacP1_vr(N6,N5,N4)*dts_wat),FLWHX))
    ENDIF

    IF(N.EQ.iVerticalDirection)THEN
      WaterMacpFlow=WaterMacpFlow+AZMIN1(VLairMacP_vr(N6,N5,N4))
    ENDIF

    if(abs(WaterMacpFlow).LT.tiny_wat)then
      WaterMacpFlow=0._r8
    endif
  ELSE
    WaterMacpFlow   = 0.0_r8
    IF(VLairMacP1_vr(N6,N5,N4).LE.0.0_r8)LInvalidMacP=.true.
  ENDIF

  IF(WaterMacpFlow.GT.0.0_r8)THEN
    HeatByFlowMacP=cpw*TKSoil1_vr(N3,N2,N1)*WaterMacpFlow
  ELSE
    HeatByFlowMacP=cpw*TKSoil1_vr(N6,N5,N4)*WaterMacpFlow
  ENDIF

  end subroutine MacporeFLow
!------------------------------------------------------------------------------------------

  subroutine WaterVaporFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectVapFlux,HeatByConvectVapFlux)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), intent(in) :: PSISV1,PSISVL
  real(r8), intent(out) :: ConvectVapFlux
  real(r8), intent(out) :: HeatByConvectVapFlux
  real(r8) :: TK11,TK12,VP1,VPL,VPY,CNV1,CNVL
  REAL(R8) :: ATCNVL,PotentialVaporFlux,MaxVaporFlux
  !     VAPOR PRESSURE AND DIFFUSIVITY IN EACH GRID CELL
  !
  
  !     THETPM,AirFillPore_Min=current, minimum air-filled porosity
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
  ConvectVapFlux       = 0._r8
  HeatByConvectVapFlux = 0._r8
  if(fixWaterLevel)return
  IF(AirFilledSoilPoreM_vr(M,N3,N2,N1).GT.AirFillPore_Min .AND. AirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min)THEN
    TK11   = TKSoil1_vr(N3,N2,N1)
    TK12   = TKSoil1_vr(N6,N5,N4)
    VP1    = vapsat(TK11)*EXP(18.0_r8*PSISV1/(RGASC*TK11))
    VPL    = vapsat(TK12)*EXP(18.0_r8*PSISVL/(RGASC*TK12))
    CNV1   = WVapDifusvitySoil_vr(N3,N2,N1)*AirFilledSoilPoreM_vr(M,N3,N2,N1)*POROQ*AirFilledSoilPoreM_vr(M,N3,N2,N1)/POROS_vr(N3,N2,N1)
    CNVL   = WVapDifusvitySoil_vr(N6,N5,N4)*AirFilledSoilPoreM_vr(M,N6,N5,N4)*POROQ*AirFilledSoilPoreM_vr(M,N6,N5,N4)/POROS_vr(N6,N5,N4)
    ATCNVL = 2.0_r8*CNV1*CNVL/(CNV1*DLYR_3D(N,N6,N5,N4)+CNVL*DLYR_3D(N,N3,N2,N1))
    !
    !     VAPOR FLUX FROM VAPOR PRESSURE AND DIFFUSIVITY,
    !     AND CONVECTIVE HEAT FLUX FROM VAPOR FLUX
!
    PotentialVaporFlux = ATCNVL*(VP1-VPL)*AREA(N,N3,N2,N1)*dts_HeatWatTP
    VPY          = (VP1*VLsoiAirPM_vr(M,N3,N2,N1)+VPL*VLsoiAirPM_vr(M,N6,N5,N4))/(VLsoiAirPM_vr(M,N3,N2,N1)+VLsoiAirPM_vr(M,N6,N5,N4))
    MaxVaporFlux = (VP1-VPY)*VLsoiAirPM_vr(M,N3,N2,N1)*dts_wat
    !out from source grid 
    IF(PotentialVaporFlux.GE.0.0_r8)THEN
      ConvectVapFlux       = AZMAX1d(AMIN1(PotentialVaporFlux,MaxVaporFlux),tiny_wat)
      HeatByConvectVapFlux = (cpw*TKSoil1_vr(N3,N2,N1)+EvapLHTC)*ConvectVapFlux
    !out from dest grid  
    ELSE
      ConvectVapFlux       = AZMIN1d(AMAX1(PotentialVaporFlux,MaxVaporFlux),tiny_wat)
      HeatByConvectVapFlux = (cpw*TKSoil1_vr(N6,N5,N4)+EvapLHTC)*ConvectVapFlux
    ENDIF
  ELSE
    ConvectVapFlux       = 0.0_r8
    HeatByConvectVapFlux = 0.0_r8
  ENDIF
  end subroutine WaterVaporFlow  
!------------------------------------------------------------------------------------------

  subroutine Solve4HeatByConduction(I,J,M,N,NY,NX,N1,N2,N3,N4,N5,N6,HeatByConvectVapFlux,HeatFluxAir2Soi)
  !
  !Description
  !Solve for heat flux due to conduction along temperature gradient
  implicit none
  integer , intent(in) :: I,J,M
  integer , intent(in) :: N,NY,NX
  integer , intent(in) :: N1,N2,N3  !source
  integer , intent(in) :: N4,N5,N6  !dest
  real(r8), intent(in) :: HeatByConvectVapFlux    !convective heat flux out of source grid into dest grid
  real(r8), intent(in) :: HeatFluxAir2Soi        !ground heat into surface layer
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
  DTKX=ABS(TKSoil1_vr(N3,N2,N1)-TKSoil1_vr(N6,N5,N4))*ppmc

  call CalcSoilThermConductivity(N1,N2,N3,DTKX,ThermCondSrc)

  call CalcSoilThermConductivity(N4,N5,N6,DTKX,ThermCondDst)

  ATCNDL=safe_adb(2.0_r8*ThermCondSrc*ThermCondDst,(ThermCondSrc*DLYR_3D(N,N6,N5,N4)+ThermCondDst*DLYR_3D(N,N3,N2,N1)))

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
  
  IF(VHeatCapacity1_vr(N3,N2,N1).GT.VHCPNX_col(NY,NX))THEN
    !surface layer and does not have significant snowpack
    IF(N3.EQ.NUM(NY,NX) .AND. VLHeatCapSnow_snvr(1,N2,N1).LE.VLHeatCapSnowMin_col(N2,N1))THEN
      TK1X=TKSoil1_vr(N3,N2,N1)-(HeatByConvectVapFlux-HeatFluxAir2Soi)/VHeatCapacity1_vr(N3,N2,N1)
      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoil1_vr(N3,N2,N1)-HeatByConvectVapFlux/VHeatCapacity1_vr(N3,N2,N1)',&
          TKSoil1_vr(N3,N2,N1),HeatByConvectVapFlux,HeatFluxAir2Soi,VHeatCapacity1_vr(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    ELSE
      TK1X=TKSoil1_vr(N3,N2,N1)-HeatByConvectVapFlux/VHeatCapacity1_vr(N3,N2,N1)
      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoil1_vr(N3,N2,N1)-HeatByConvectVapFlux/VHeatCapacity1_vr(N3,N2,N1)',&
          TKSoil1_vr(N3,N2,N1),HeatByConvectVapFlux,VHeatCapacity1_vr(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    ENDIF
  ELSE
    TK1X=TKSoil1_vr(N3,N2,N1)
  ENDIF

  IF(VHeatCapacity1_vr(N6,N5,N4).GT.ZEROS(NY,NX))THEN
    TKLX=TKSoil1_vr(N6,N5,N4)+HeatByConvectVapFlux/VHeatCapacity1_vr(N6,N5,N4)
  ELSE
    TKLX=TKSoil1_vr(N6,N5,N4)
  ENDIF
  
  if(VHeatCapacity1_vr(N3,N2,N1)+VHeatCapacity1_vr(N6,N5,N4)>0._r8)then
    !temperature at the interface between (N3,N2,N1) and (N6,N5,N4)
    TKY=(VHeatCapacity1_vr(N3,N2,N1)*TK1X+VHeatCapacity1_vr(N6,N5,N4)*TKLX)/(VHeatCapacity1_vr(N3,N2,N1)+VHeatCapacity1_vr(N6,N5,N4))
  ELSE
    TKY=(TK1X+TKLX)/2._r8
  endif 

  HFLWX = (TK1X-TKY)*VHeatCapacity1_vr(N3,N2,N1)*dts_wat
  HFLWC = ATCNDL*(TK1X-TKLX)*AREA(N,N3,N2,N1)*dts_HeatWatTP
  IF(HFLWC.GE.0.0_r8)THEN
    !from (N3,N2,N1) into (N6,N5,N4)
    HeatCondSoi=AZMAX1(AMIN1(HFLWX,HFLWC))
  ELSE
    !from (N6,N5,N4) into (N3,N2,N1)
    HeatCondSoi=AZMIN1(AMAX1(HFLWX,HFLWC))
  ENDIF
  HeatFlow2Soili_3D(N,N6,N5,N4)=HeatFlow2Soili_3D(N,N6,N5,N4)+HeatCondSoi
!  if(N6==16)then
!  write(411,*)I+J/24.,M,'heatcond',HeatFlow2Soili_3D(N,N6,N5,N4),N6,N5,N4,N
!  endif
  end subroutine Solve4HeatByConduction  
!------------------------------------------------------------------------------------------
  subroutine DischargeOverWaterTBL(I,J,N,N1,N2,N3,M4,M5,M6,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL,&
    RechargSubSurf,RechargRateWTBL,DPTHH,XN,TargetWaterTBL,FracVolBelowWTBL)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: N
  integer, intent(in) :: N1,N2,N3  !source grid
  integer, intent(in) :: M4,M5,M6  !dest grid
  logical, intent(in) :: DoMicPDischarg2ExtWTBL
  logical, intent(in) :: DoMacPDischarg2ExtWTBL
  real(r8), intent(in):: RechargSubSurf,RechargRateWTBL
  real(r8), intent(in):: XN                         !flow direction, <0 out of source grid
  real(r8), intent(in):: DPTHH                      !depth to layer macropore water 
  real(r8), intent(in):: TargetWaterTBL
  real(r8), intent(in) :: FracVolBelowWTBL
  real(r8) :: FLWT,MacPoreDischarg,FLWTH,PSISWD,PSISWT,PSISWTH
  real(r8) :: watflx, heatflx


!     MICROPORE DISCHARGE ABOVE WATER TABLE
!
!     DoMicPDischarg2ExtWTBL=micropore discharge flag to natural water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     WaterTBLSlope_col=water table slope
!     PSISWT=water potential driving micropore discharge
!     PSISoilMatricPtmp_vr,PSISO=matric,osmotic water potential
!     DPTH,ExtWaterTable=depth to layer midpoint,natural water table
!     DepzIntWTBL_col=depth to internal water table
!     FLWL=micropore discharge to natural water table
!     HFLWL=convective heat from discharge to natural water table
!     HydroCond_3D=saturated hydraulic conductivity
!     FracLayVolBelowExtWTBL_vr=fraction of layer below natural water table
!
  IF(DoMicPDischarg2ExtWTBL .AND. (.not.isclose(RechargRateWTBL,0._r8)))THEN
    PSISWD = XN*0.005_r8*SLOPE(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope_col(N2,N1))
    PSISWT = AZMIN1(-PSISoilMatricPtmp_vr(N3,N2,N1)-0.03_r8*PSISoilOsmotic_vr(N3,N2,N1) &
      +mGravAccelerat*(SoilDepthMidLay_vr(N3,N2,N1)-ExtWaterTable_col(N2,N1)) &
      -mGravAccelerat*AZMAX1(SoilDepthMidLay_vr(N3,N2,N1)-DepzIntWTBL_col(N2,N1)))

    IF(PSISWT.LT.0.0_r8)PSISWT=PSISWT-PSISWD
    FLWT=PSISWT*HydroCond_3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*(1.0_r8-FracVolBelowWTBL)/(RechargSubSurf+1.0_r8) &
      *RechargRateWTBL*dts_HeatWatTP
    watflx                          = XN*FLWT
    heatflx                         = cpw*TKSoil1_vr(N3,N2,N1)*watflx
    WaterFlow2Micpt_3D(N,M6,M5,M4)  = WaterFlow2Micpt_3D(N,M6,M5,M4)+watflx
    WaterFlow2MicptX_3D(N,M6,M5,M4) = WaterFlow2MicptX_3D(N,M6,M5,M4)+watflx
    HeatFlow2Soili_3D(N,M6,M5,M4)   = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx

    QDrain_col(N2,N1)               = QDrain_col(N2,N1) + watflx
    HeatDrain_col(N2,N1)            = HeatDrain_col(N2,N1)+heatflx
  ENDIF
!
!     MACROPORE DISCHARGE ABOVE WATER TABLE
!
!     DoMacPDischarg2ExtWTBL=macropore discharge flag to natural water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     WaterTBLSlope_col=water table slope
!     PSISWTH=water potential driving macropore discharge
!     PSISO=osmotic water potential
!     ExtWaterTable=natural water table
!     DepzIntWTBL_col=depth to internal water table
!     FLWTH,MacPoreDischarg=macropore discharge unltd,ltd by macropore water
!     HydroCondMacP1=macropore hydraulic conductivity
!     WaterFlow2Macpt_3D=macropore discharge to natural water table
!     HFLWL=convective heat from discharge to natural water table
!     HydroCond_3D=saturated hydraulic conductivity
  !     FracLayVolBelowExtWTBL_vr=fraction of layer below natural water table
!
  IF(DoMacPDischarg2ExtWTBL .AND. (.not.isclose(RechargRateWTBL,0._r8)) .AND. VLMacP1_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    PSISWD  = XN*0.005_r8*SLOPE(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope_col(N2,N1))
    PSISWTH = -0.03_r8*PSISoilOsmotic_vr(N3,N2,N1)+mGravAccelerat*(DPTHH-TargetWaterTBL) &
      -mGravAccelerat*AZMAX1(DPTHH-DepzIntWTBL_col(N2,N1))

    IF(PSISWTH.LT.0.0_r8)PSISWTH=PSISWTH-PSISWD
    FLWTH=PSISWTH*HydroCondMacP1_vr(N3,N2,N1)*AREA(N,N3,N2,N1) &
      *(1.0_r8-FracVolBelowWTBL)/(RechargSubSurf+1.0)*RechargRateWTBL*dts_HeatWatTP

    MacPoreDischarg=AMAX1(FLWTH,AZMIN1(-(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
      +WaterFlow2Macpt_3D(3,N3,N2,N1)-WaterFlow2Macpt_3D(3,N3+1,N2,N1))))
    watflx                         = XN*MacPoreDischarg
    heatflx                        = cpw*TKSoil1_vr(N3,N2,N1)*watflx
    WaterFlow2Macpt_3D(N,M6,M5,M4) = WaterFlow2Macpt_3D(N,M6,M5,M4)+watflx
    HeatFlow2Soili_3D(N,M6,M5,M4)  = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx

    QDrain_col(N2,N1)              = QDrain_col(N2,N1) + watflx
    HeatDrain_col(N2,N1)           = HeatDrain_col(N2,N1)+heatflx
  ENDIF
  end subroutine DischargeOverWaterTBL
!------------------------------------------------------------------------------------------

  SUBROUTINE RechargeFromExtWaterTBL(I,J,M,N,N1,N2,N3,M4,M5,M6,DPTHH,RechargSubSurf,&
    RechargRateWTBL,XN,AirfMicP,VOLPX2,AirfMacP,str_dir)!,donot_drain)
  !
  !subsurface recharge to soil micropore and macropores from external water table
  !it considers the existence of frozen layers
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: N
  integer, intent(in) :: N1,N2,N3   !source grid
  integer, intent(in) :: M4,M5,M6   !dest grid
  real(r8),intent(in) :: RechargRateWTBL,DPTHH,RechargSubSurf
  real(r8),intent(in) :: XN   !flow direction, > 0 (in), < 0 (out)
  real(r8),intent(inout):: AirfMicP   !air-filled space in grid (N3,N2,N1)
  real(r8),intent(inout) :: VOLPX2,AirfMacP
  character(len=*),intent(in) :: str_dir  
!  logical, intent(in) :: donot_drain

  character(len=*), parameter :: subname='RechargeFromExtWaterTBL'
  real(r8) :: FLWU,FLWUL,FLWUX,FLWUH,FLWUHL  
  real(r8) :: PSISUT,PSISUTH,PSISWD
  real(r8) :: heatflx,watflx,watflxu

  call PrintInfo('beg '//subname)
  !     MICROPORE RECHARGE BELOW WATER TABLE

  IF(SoilDepthMidLay_vr(N3,N2,N1).GE.ExtWaterTable_col(N2,N1)                      & !deeper than water table
    .AND. ActiveLayDepZ_col(N2,N1).GT.ExtWaterTable_col(N2,N1)                     & !active layer deepd enough
    .AND. SoilDepthMidLay_vr(N3,N2,N1).LT.ActiveLayDepZ_col(N2,N1)                 & !inside active layer
    .AND. (AirfMicP.GT.ZEROS2(N2,N1) .OR. SoilBulkDensity_vr(N3,N2,N1).LE.ZERO)    & !able to accept water
    .AND. VLairMicP_vr(N3,N2,N1).GT.0.0_r8                                         & !source grid is unsaturated
    .AND. (.not.isclose(RechargRateWTBL,0._r8)))THEN                                 !water exchange with watertable enabled

    PSISWD = XN*0.005_r8*SLOPE(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope_col(N2,N1))
    PSISUT = AZMAX1(-PSISoilMatricPtmp_vr(N3,N2,N1)-0.03_r8*PSISoilOsmotic_vr(N3,N2,N1)+&
      mGravAccelerat*(SoilDepthMidLay_vr(N3,N2,N1)-ExtWaterTable_col(N2,N1)))

    !outflow/drainage  (>0.)
    IF(PSISUT.GT.0.0_r8)PSISUT=PSISUT+PSISWD

    FLWU = PSISUT*HydroCond_3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*FracLayVolBelowExtWTBL_vr(N3,N2,N1) &
      /(RechargSubSurf+1.0)*RechargRateWTBL*dts_HeatWatTP

    !within a time step, the incoming flow cannot exceed avaiable air-filled pores   
    IF(SoilBulkDensity_vr(N3,N2,N1).GT.ZERO)THEN
      FLWUL = AMIN1(FLWU,AirfMicP)
      FLWUX = AMIN1(FLWU,VOLPX2)
    ELSE
      FLWUL = FLWU
      FLWUX = FLWU
    ENDIF
    watflx                 = XN*FLWUL
    watflxu                = XN*FLWUX
    heatflx                = cpw*TKSoil1_vr(N3,N2,N1)*watflx
!        if(donot_drain)then

!        else
      if(N3==M6 .and. N2==M5 .and. N1==M4)then
        QDischar_col(N2,N1)    = QDischar_col(N2,N1)-watflx
        QDischarM_col(N2,N1)   = QDischarM_col(N2,N1)-watflx
        HeatDischar_col(N2,N1) = HeatDischar_col(N2,N1)-heatflx
      else
        QDischar_col(N2,N1)    = QDischar_col(N2,N1)+watflx
        QDischarM_col(N2,N1)   = QDischarM_col(N2,N1)+watflx
        HeatDischar_col(N2,N1) = HeatDischar_col(N2,N1)+heatflx
      endif
      WaterFlow2Micpt_3D(N,M6,M5,M4)  = WaterFlow2Micpt_3D(N,M6,M5,M4)+watflx
      WaterFlow2MicptX_3D(N,M6,M5,M4) = WaterFlow2MicptX_3D(N,M6,M5,M4)+watflxu
      HeatFlow2Soili_3D(N,M6,M5,M4)   = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx
      AirfMicP                        = AirfMicP-XN*watflx
      VOLPX2                          = VOLPX2-XN*watflxu
!        endif

!        if(I>=267 .and. J>=23)write(214,*)I+J/24.,M,N3,N2,N1,M6,M5,M4,'dis1'//str_dir,watflx,XN,N !,donot_drain
  ENDIF
!
  !     MACROPORE RECHARGE BELOW WATER TABLE


  IF(DPTHH.GT.ExtWaterTable_col(N2,N1)                            & !deeper than water table
    .AND. ActiveLayDepZ_col(N2,N1).GT.ExtWaterTable_col(N2,N1)    & !active layer below water table
    .AND. SoilDepthMidLay_vr(N3,N2,N1).LT.ActiveLayDepZ_col(N2,N1) & !midlayer depth above active water layer
    .AND. AirfMacP.GT.ZEROS2(N2,N1)                               & !macropore has air-filled fraction
    .AND. (.not.isclose(RechargRateWTBL,0.0_r8)))THEN               !recharge is on

    PSISWD  = XN*0.005*SLOPE(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTBLSlope_col(N2,N1))
    PSISUTH = -0.03_r8*PSISoilOsmotic_vr(N3,N2,N1)+mGravAccelerat*(DPTHH-ExtWaterTable_col(N2,N1))

    !outflow/drainage
    IF(PSISUTH.GT.0.0_r8)PSISUTH = PSISUTH+PSISWD
    FLWUH  = PSISUTH*HydroCondMacP1_vr(N3,N2,N1)*AREA(N,N3,N2,N1)*FracLayVolBelowExtWTBL_vr(N3,N2,N1) &
      /(RechargSubSurf+1.0_r8)*RechargRateWTBL*dts_HeatWatTP
    !is *dts_wat needed below?  
    FLWUHL                         = AMIN1(FLWUH,AirfMacP)
    watflx                         = XN*FLWUHL
    heatflx                        = cpw*TKSoil1_vr(N3,N2,N1)*watflx
    WaterFlow2Macpt_3D(N,M6,M5,M4) = WaterFlow2Macpt_3D(N,M6,M5,M4)+watflx
    HeatFlow2Soili_3D(N,M6,M5,M4)  = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx
    AirfMacP                       = AirfMacP-XN*watflx

!        write(214,*)I+J/24.,'dis2',watflx,XN*watflx        
    if(N3==M6 .and. N2==M5 .and. N1==M4)then
      QDischar_col(N2,N1)    = QDischar_col(N2,N1)-watflx
      QDischarM_col(N2,N1)   = QDischarM_col(N2,N1)-watflx
      HeatDischar_col(N2,N1) = HeatDischar_col(N2,N1)-heatflx
    else
      QDischar_col(N2,N1)    = QDischar_col(N2,N1)+watflx
      QDischarM_col(N2,N1)   = QDischarM_col(N2,N1)+watflx
      HeatDischar_col(N2,N1) = HeatDischar_col(N2,N1)+heatflx
    endif
  ENDIF
  call PrintInfo('end '//subname)    
  end SUBROUTINE RechargeFromExtWaterTBL
!------------------------------------------------------------------------------------------

  subroutine FreezeThawIterateM(I,J,M,NY,NX,L,N1,N2,N3)
  !
  !Description
  !Applying freezethaw in cell (N3,N2,N1)
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NY,NX,L,N1,N2,N3
  real(r8) :: VLHeatCapacityBX
  real(r8) :: PSISMX
  real(r8) :: TFREEZ           !freezing point 
  real(r8) :: VLWatMicP1X      !micropore liquid water before freeze/thaw
  real(r8) :: VLWatMacP1X      !macropore liquid water before freeze/thaw
  real(r8) :: VLHeatCapacityAX !
  real(r8) :: TK1App           !apparent temperature before freeze/thaw
  real(r8) :: MacPIceHeatFlxFrezPt,MacPIceHeatFlxFrez
  real(r8) :: ENGY1
  real(r8) :: MicPIceHeatFlxFrez   !actual latent heat from micropore freeze-thaw
  real(r8) :: MicPIceHeatFlxFrezPt !potential heat flux supporting micropore freeze-thaw
  real(r8) :: VLHeatCapacityX  !heat capacity before freeze/thaw
  real(r8) :: dcpo             !heat capacity due to living root [MJ/K/d-2]
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
!     HeatIrrigation1_vr=subsurface convective heat flux
!     MicPIceHeatFlxFrezPt,MicPIceHeatFlxFrez=latent heat from micro freeze-thaw unltd,ltd by water,ice
!     FIceThawedMicP_vr=soil water flux from micropore freeze-thaw
!     MacPIceHeatFlxFrezPt,SoilWatFrezHeatRelease_vr=latent heat from macro freeze-thaw unltd,ltd by water,ice
!     FIceThawedMacP_vr=soil water flux from macropore freeze-thaw
!
  if(plantOM4Heat)then
    dcpo = cpo*RootMassElm_vr(ielmc,N3,N2,N1)
  else
    dcpo=0._r8
  endif
  PSISMX          = PSISoilMatricPtmp_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  !apply the Clausis-Clapeyron equation for depressed freezing temperature, 
  !Cary and Mayland (1972), 9.0959E+04_r8~=273.15*LtHeatIceMelt  
  TFREEZ          = -9.0959E+04_r8/(PSISMX-LtHeatIceMelt)
!  TFREEZ = TFICE
  !obtain the potential water
  VLWatMicP1X     = VLWatMicP1_vr(N3,N2,N1)+TWatFlow2MicP_3DM_vr(N3,N2,N1)+FWatExMacP2MicPiM_vr(N3,N2,N1)+FWatIrrigate2MicP1_vr(N3,N2,N1)
  VLWatMacP1X     = VLWatMacP1_vr(N3,N2,N1)+TWaterFlow2Macpt_3DM_vr(N3,N2,N1)-FWatExMacP2MicPiM_vr(N3,N2,N1)
  ENGY1           = VHeatCapacity1_vr(N3,N2,N1)*TKSoil1_vr(N3,N2,N1)
  VLHeatCapacityX = VHeatCapacitySoilM_vr(N3,N2,N1)+cpw*(VLWatMicP1X+VLWatMacP1X)+cpi*(VLiceMicP1_vr(N3,N2,N1)+VLiceMacP1_vr(N3,N2,N1))

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    VLHeatCapacityX=VLHeatCapacityX+dcpo
    TK1App=(ENGY1+THeatFlow2Soil_3DM_vr(N3,N2,N1)+HeatIrrigation1_vr(N3,N2,N1))/VLHeatCapacityX
!
!     FREEZE-THAW IN SOIL LAYER MACROPORE FROM NET CHANGE IN SOIL
!     LAYER HEAT STORAGE
!   TFice=frozen temperature
!   There is freeze-thaw in macropore, it is assumed there is no freezing temperature depression in 
!   macropores, therefore, water will first freeze in macropores.
    IF((TK1App.LT.TFICE .AND. VLWatMacP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)) &
      .OR.(TK1App.GT.TFICE .AND. VLiceMacP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)))THEN
      
      VLHeatCapacityBX     = cpw*VLWatMacP1X+cpi*VLiceMacP1_vr(L,NY,NX)
      !where is the following equation come from?
      MacPIceHeatFlxFrezPt = VLHeatCapacityBX*(TFICE-TK1App)/((1.0_r8+6.2913E-03_r8*TFICE))*dts_wat 
!          /(1.0_r8-0._r8*0.10_r8*PSISMX)

      !Ice thawed, absorb heat
      IF(MacPIceHeatFlxFrezPt.LT.0.0_r8)THEN        
        MacPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMacP1_vr(N3,N2,N1)*dts_wat,MacPIceHeatFlxFrezPt)        
      !Water frozen, release heat
      ELSE        
        MacPIceHeatFlxFrez=AMIN1(LtHeatIceMelt*VLWatMacP1X*dts_wat,MacPIceHeatFlxFrezPt)
      ENDIF
      FIceThawedMacP_vr(N3,N2,N1)=-MacPIceHeatFlxFrez/LtHeatIceMelt
    !No freeze-thaw in macropore 
    ELSE
      MacPIceHeatFlxFrez     = 0.0_r8
      FIceThawedMacP_vr(N3,N2,N1) = 0.0_r8
    ENDIF

    !freeze-thaw in micropore
    IF((TK1App.LT.TFREEZ .AND. VLWatMicP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)) &
      .OR.(TK1App.GT.TFREEZ .AND. VLiceMicP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)))THEN
      
      VLHeatCapacityAX     = VHeatCapacitySoilM_vr(N3,N2,N1)+cpw*VLWatMicP1X+cpi*VLiceMicP1_vr(N3,N2,N1)
      if(VLHeatCapacityAX>0._r8)VLHeatCapacityAX=VLHeatCapacityAX+dcpo
      MicPIceHeatFlxFrezPt = VLHeatCapacityAX*(TFREEZ-TK1App)/((1.0_r8+6.2913E-03_r8*TFREEZ))*dts_wat

      !fusion energy absorb (<0) in thaw
      IF(MicPIceHeatFlxFrezPt.LT.0.0_r8)THEN
        MicPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1_vr(N3,N2,N1)*dts_wat,MicPIceHeatFlxFrezPt)
      !fusion energy release (>0) in freeze
      ELSE
        MicPIceHeatFlxFrez=AMIN1(LtHeatIceMelt*VLWatMicP1X*dts_wat,MicPIceHeatFlxFrezPt)
      ENDIF
      FIceThawedMicP_vr(N3,N2,N1)=-MicPIceHeatFlxFrez/LtHeatIceMelt
    ELSE
      MicPIceHeatFlxFrez     = 0.0_r8
      FIceThawedMicP_vr(N3,N2,N1) = 0.0_r8
    ENDIF

  ELSE
    MicPIceHeatFlxFrez          = 0.0_r8
    FIceThawedMicP_vr(N3,N2,N1) = 0.0_r8
    MacPIceHeatFlxFrez          = 0.0_r8
    FIceThawedMacP_vr(N3,N2,N1) = 0.0_r8
  ENDIF
  !total heat released due to phase change 
  SoilWatFrezHeatRelease_vr(N3,N2,N1)=MicPIceHeatFlxFrez+MacPIceHeatFlxFrez !+ &
    !(cpw-cpi*DENSICE)*(FIceThawedMicP_vr(N3,N2,N1)+FIceThawedMacP_vr(N3,N2,N1))*TFREEZ
  
  !
  !     TOTAL AND HOURLY ACCUMULATED FREEZE-THAW FLUXES
  !
  !     THAW,TLIceThawMacP=hourly accumulated freeze-thaw flux in micropores,macropores
  !     HTHAW=hourly accumulated freeze-thaw latent heat flux
  !     TWFLXL,TMLiceThawMacP=total accumulated freeze-thaw in micropores,macropores
  !     TLPhaseChangeHeat2Soi1_vr=total latent heat flux from melting

  ! Summarize fluxes for the current iteration
  TLIceThawMicP_vr(N3,N2,N1)          = TLIceThawMicP_vr(N3,N2,N1)+FIceThawedMicP_vr(N3,N2,N1)
  TLIceThawMacP_vr(N3,N2,N1)          = TLIceThawMacP_vr(N3,N2,N1)+FIceThawedMacP_vr(N3,N2,N1)
  TLPhaseChangeHeat2Soi_vr(N3,N2,N1)  = TLPhaseChangeHeat2Soi_vr(N3,N2,N1)+SoilWatFrezHeatRelease_vr(N3,N2,N1)
  TMLiceThawedMicP_vr(N3,N2,N1)       = TMLiceThawedMicP_vr(N3,N2,N1)+FIceThawedMicP_vr(N3,N2,N1)
  TMLiceThawedMacP_vr(N3,N2,N1)       = TMLiceThawedMacP_vr(N3,N2,N1)+FIceThawedMacP_vr(N3,N2,N1)
  TLPhaseChangeHeat2Soi1_vr(N3,N2,N1) = TLPhaseChangeHeat2Soi1_vr(N3,N2,N1)+SoilWatFrezHeatRelease_vr(N3,N2,N1)
  end subroutine FreezeThawIterateM

end module WatsubMod
