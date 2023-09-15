module WatsubMod
!!
! Description:
! Do water and enerby balance calculation.
! The module diagnoses the mass and energy fluxes associated
! with soil/snow water (vapor, liquid and ice) and energy, and update
! them in redistmod.F90

  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : GravAcceleration=>DAT_CONST_G
  use abortutils   , only : endrun, print_info
  use minimathmod  , only : isclose, isclose,safe_adb,vapsat,AZMAX1,AZMIN1,AZMAX1t
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
  use SurfPhysData, only : InitSurfPhysData,DestructSurfPhysData  
  implicit none

  private

  character(len=*), parameter :: mod_filename = __FILE__
  integer, parameter :: iewstdir=1   !east-west direction
  integer, parameter :: insthdir=2   !north-south direction
  integer, parameter :: ivertdir=3   !vertical direction
  real(r8), parameter :: mGravAcceleration=1.e-3_r8*GravAcceleration  !gravitational constant devided by 1000.
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
  real(r8):: RAR1(JY,JX),HeatFlux2Ground(JY,JX)
  REAL(R8):: FKSATS(JY,JX)

  REAL(R8) :: TopLayWatVol(JY,JX)
! begin_execution
!
  curday=i;curhour=j

  call PrepWaterEnergyBalance(I,J,NHW,NHE,NVN,NVS,RAR1)

  call InitSoilHydrauics(NHW,NHE,NVN,NVS)

! DYNAMIC LOOP FOR FLUX CALCULATIONS
  ! iterate for NPH times
  D3320: DO M=1,NPH

    call FWDCopyTopLayerWatVolMit(NHW,NHE,NVN,NVS,TopLayWatVol)

    ! run surface energy balance model, uses RAR1
    call SurfacePhysModel(M,NHE,NHW,NVS,NVN,RAR1,FKSATS,HeatFlux2Ground,TopLayWatVol)

    call CopySoilWatVolMit(NHW,NHE,NVN,NVS,TopLayWatVol)
        
    ! do 3D water flow
    call Subsurface3DFlowMit(M,NHW,NHE,NVN,NVS,FKSATS,HeatFlux2Ground)

    call LateralWatHeatExchMit(M,NHW,NHE,NVN,NVS,FKSATS)
!
    IF(M.NE.NPH)THEN
!     intermediate iteration
      call UpdateStateFluxAtM(M,NHW,NHE,NVN,NVS)
    ELSE
!     last iteration
      call UpdateFluxAtExit(NHW,NHE,NVN,NVS)
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
      TopLayWatVol(NY,NX)= VOLW2(NUM(NY,NX),NY,NX)
    ENDDO
  ENDDO
  end subroutine FWDCopyTopLayerWatVolMit

!------------------------------------------------------------------------------------------  

  subroutine CopySoilWatVolMit(NHW,NHE,NVN,NVS,TopLayWatVol)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: TopLayWatVol
  integer :: L,NY,NX
  
  !VOLW2 will be updated in soil surface model

  DO NX=NHW,NHE
    DO  NY=NVN,NVS  
      DO L=NUM(NY,NX)+1,NL(NY,NX)        
        VOLW2(L,NY,NX)=VLWatMicP1(L,NY,NX)
      ENDDO  
      VOLW2(NUM(NY,NX),NY,NX)=TopLayWatVol(NY,NX)
    ENDDO
  ENDDO
  end subroutine CopySoilWatVolMit


!------------------------------------------------------------------------------------------  
  subroutine LocalCopySoilVars(I,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS
  integer :: NY,NX

  integer :: L,LWDPTH
  real(r8) :: VLTSoiPore

  DX995: DO NX=NHW,NHE
    DX990: DO NY=NVN,NVS

    ! CDPTH=depth to bottom of soil layer
    ! WDPTH,LWDPTH=depth,layer of subsurface irrigation

      !identify the layer where irrigation is applied
      D65: DO L=NUM(NY,NX),NL(NY,NX)
        IF(CumDepth2LayerBottom(L,NY,NX).GE.WDPTH(I,NY,NX))THEN
          LWDPTH=L
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
    !   VLSoilPoreMicP,VOLT=soil,total volumes
    !   WP=wilting point
    !   THETW*,THETI*,THETP*=water,ice,air-filled porosity
    !   VLHeatCapacity,VHCM=volumetric heat capacities of total volume, solid
    !   VLHeatCapacityA,VLHeatCapacityB=volumetric heat capacities of micropore,macropore
    !
        PSISM1(L,NY,NX)=PSISoilMatricP(L,NY,NX)
        VLMicP1(L,NY,NX)=VLMicP(L,NY,NX)
        VLWatMicP1(L,NY,NX)=VLWatMicP(L,NY,NX)
        VLWatMicPX1(L,NY,NX)=VLWatMicPX(L,NY,NX)
        VLiceMicP1(L,NY,NX)=VLiceMicP(L,NY,NX)
        VLWatMacP1(L,NY,NX)=VLWatMacP(L,NY,NX)
        VLiceMacP1(L,NY,NX)=VLiceMacP(L,NY,NX)
        
        IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
          VLairMicP(L,NY,NX)=VLMicP1(L,NY,NX)-VLWatMicP1(L,NY,NX)-VLiceMicP1(L,NY,NX)
          VLairMicPc(L,NY,NX)=AZMAX1(VLairMicP(L,NY,NX))
        ELSE
          VLairMicP(L,NY,NX)=0.0_r8
          VLairMicPc(L,NY,NX)=0.0_r8
        ENDIF

        !FVOLAH accounts for clay swelling effect due to change in micropore water, but it is set to zero
        VLMacP1(L,NY,NX)=AZMAX1(VLMacP(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX) &
          *(safe_adb(VLWatMicP1(L,NY,NX),VLSoilMicP(L,NY,NX))-WiltPoint(L,NY,NX))*VGeomLayer(L,NY,NX))

        IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
          VLairMacP(L,NY,NX)=VLMacP1(L,NY,NX)-VLWatMacP1(L,NY,NX)-VLiceMacP1(L,NY,NX)
          VLairMacPc(L,NY,NX)=AZMAX1(VLairMacP(L,NY,NX))
        ELSE
          VLairMacP(L,NY,NX)=0.0_r8
          VLairMacPc(L,NY,NX)=0.0_r8
        ENDIF        
        VLWatMicPM(1,L,NY,NX)=VLWatMicP1(L,NY,NX)
        VLWatMacPM(1,L,NY,NX)=VLWatMacP1(L,NY,NX)
        VLsoiAirPM(1,L,NY,NX)=VLairMicPc(L,NY,NX)+VLairMacPc(L,NY,NX)+THETPI*(VLiceMicP1(L,NY,NX)+VLiceMacP1(L,NY,NX))
        
        VLTSoiPore=VLSoilMicP(L,NY,NX)+VLMacP1(L,NY,NX)
        IF(VLTSoiPore.GT.ZEROS2(NY,NX))THEN
          !fraction as water
          FracSoiPAsWat(L,NY,NX)=AZMAX1t((VLWatMicP1(L,NY,NX)+VLWatMacP1(L,NY,NX))/VLTSoiPore)
          !fraction as ice
          FracSoiPAsIce(L,NY,NX)=AZMAX1t((VLiceMicP1(L,NY,NX)+VLiceMacP1(L,NY,NX))/VLTSoiPore)
          !fraction as air
          FracSoiPAsAir(L,NY,NX)=AZMAX1t((VLairMicPc(L,NY,NX)+VLairMacPc(L,NY,NX))/VLTSoiPore)
        ELSE
          FracSoiPAsWat(L,NY,NX)=POROS(L,NY,NX)
          FracSoiPAsIce(L,NY,NX)=0.0_r8
          FracSoiPAsAir(L,NY,NX)=0.0_r8
        ENDIF
        THETPM(1,L,NY,NX)=FracSoiPAsAir(L,NY,NX)
        IF(VLMicP1(L,NY,NX)+VLMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          FracSoilAsAirt(L,NY,NX)=AZMAX1((VLairMicPc(L,NY,NX)+VLairMacPc(L,NY,NX))/(VLMicP1(L,NY,NX)+VLMacP1(L,NY,NX)))
        ELSE
          FracSoilAsAirt(L,NY,NX)=0.0_r8
        ENDIF
        !VLHeatCapacity=total heat capacity
        !VLHeatCapacityA=heat capcity without macropore water/ice
        !VLHeatCapacityB=heat capacity for macropore water/ice
        VLHeatCapacityA(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+cpw*VLWatMicP1(L,NY,NX)+cpi*VLiceMicP1(L,NY,NX)    
        VLHeatCapacityB(L,NY,NX)=cpw*VLWatMacP1(L,NY,NX)+cpi*VLiceMacP1(L,NY,NX)      
        VLHeatCapacity(L,NY,NX)=VLHeatCapacityA(L,NY,NX)+VLHeatCapacityB(L,NY,NX)
    !
    !   MACROPOROSITY
    !
    !   FMAC,SoilFracAsMicP=macropore,micropore volume fractions
    !   CNDH*=macropore hydraulic conductivity
    !   TKS,TK1=soil temperature
    !   FLU,HeatIrrigation=subsurface water,convective heat fluxes
    !   AREAU,AREAD=fractions of layer below natural,artifl water table
    !
        IF(VLMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          SoilFracAsMacP1(L,NY,NX)=SoilFracAsMacP(L,NY,NX)*VLMacP1(L,NY,NX)/VLMacP(L,NY,NX)
          HydroCondMacP1(L,NY,NX)=HydroCondMacP(L,NY,NX)*(VLMacP1(L,NY,NX)/VLMacP(L,NY,NX))**2._r8
        ELSE
          SoilFracAsMacP1(L,NY,NX)=0.0_r8
          HydroCondMacP1(L,NY,NX)=0.0_r8
        ENDIF
        SoilFracAsMicP(L,NY,NX)=1.0_r8-SoilFracAsMacP1(L,NY,NX)
        TKSoi1(L,NY,NX)=TKS(L,NY,NX)
        if(TKSoi1(L,NY,NX)>1.e3_r8)then
          write(*,*)'TKS(L,NY,NX)',L,TKS(L,NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
        !LWDPTH=layer number where irrigation is applied
        IF(L.EQ.LWDPTH)THEN
          FWatIrrigate2MicP(L,NY,NX)=PRECU(NY,NX)
          HeatIrrigation(L,NY,NX)=cpw*TairK(NY,NX)*PRECU(NY,NX)
          FWatIrrigate2MicP1(L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*dts_HeatWatTP
          HeatIrrigation1(L,NY,NX)=HeatIrrigation(L,NY,NX)*dts_HeatWatTP
        ELSE
          FWatIrrigate2MicP(L,NY,NX)=0.0_r8
          HeatIrrigation(L,NY,NX)=0.0_r8
          FWatIrrigate2MicP1(L,NY,NX)=0.0_r8
          HeatIrrigation1(L,NY,NX)=0.0_r8
        ENDIF
        IF(CumDepth2LayerBottom(L,NY,NX).GE.DTBLX(NY,NX))THEN
          AREAU(L,NY,NX)=AMIN1(1.0_r8,AZMAX1(safe_adb(CumDepth2LayerBottom(L,NY,NX)-DTBLX(NY,NX),DLYR(3,L,NY,NX))))
        ELSE
          AREAU(L,NY,NX)=0.0_r8
        ENDIF
        IF(CumDepth2LayerBottom(L,NY,NX).GE.DTBLY(NY,NX))THEN
          AreaUnderWaterTBL(L,NY,NX)=AMIN1(1.0_r8,AZMAX1(safe_adb(CumDepth2LayerBottom(L,NY,NX)-DTBLY(NY,NX),DLYR(3,L,NY,NX))))
        ELSE
          AreaUnderWaterTBL(L,NY,NX)=0.0_r8
        ENDIF
      ENDDO D30
    ENDDO DX990
  ENDDO DX995

  end subroutine LocalCopySoilVars


!------------------------------------------------------------------------------------------

  subroutine PrepWaterEnergyBalance(I,J,NHW,NHE,NVN,NVS,RAR1)
  implicit none
  integer :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),dimension(:,:),intent(out):: RAR1
  integer :: NY,NX

! begin_execution

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
    !make a local copy of the upper boundary index
      NUM(NY,NX)=NU(NY,NX)
    ENDDO D9990
  ENDDO D9995

  call StageSurfStateVars(I,J,NHW,NHE,NVN,NVS,RAR1)

  call LocalCopySoilVars(I,NHW,NHE,NVN,NVS)

  end subroutine PrepWaterEnergyBalance
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
        D40: DO N=NCN(NY,NX),3
          N1=NX
          N2=NY
          N3=L
! in the EW direction
          IF(N.EQ.iewstdir)THEN
            IF(NX.EQ.NHE)THEN
              cycle
            ELSE
              N4=NX+1
              N5=NY
              N6=L
            ENDIF
! in the NS direction
          ELSEIF(N.EQ.insthdir)THEN
            IF(NY.EQ.NVS)THEN
              cycle
            ELSE
              N4=NX
              N5=NY+1
              N6=L
            ENDIF
! in the vertical direction
          ELSEIF(N.EQ.ivertdir)THEN
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
          IF(HydroCondMacP1(N3,N2,N1).GT.ZERO.AND.HydroCondMacP1(N6,N5,N4).GT.ZERO)THEN
            AVCNHL(N,N6,N5,N4)=2.0_r8*HydroCondMacP1(N3,N2,N1)*HydroCondMacP1(N6,N5,N4) &
              /(HydroCondMacP1(N3,N2,N1)*DLYR(N,N6,N5,N4)+HydroCondMacP1(N6,N5,N4) &
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

  subroutine Subsurface3DFlowMit(M,NHW,NHE,NVN,NVS,FKSAT,HeatFlux2Ground)
  implicit none
  integer, intent(in)  :: M,NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: FKSAT(:,:)
  real(r8), dimension(:,:),intent(in) :: HeatFlux2Ground(:,:)
  integer :: N,N1,N2,N3,N4,N5,N6,L,LL,K1,KL,NY,NX
  real(r8) :: DTKX
  real(r8) :: WTHET1,FCDX,FCLX,FCX
  real(r8) :: PSISV1,TKY,PSDX
  real(r8) :: WPLX,WPX,WTHET2
  real(r8) :: ThermCondDst,ThermCondSrc

  real(r8) :: AVE_CONDUCTANCE,ATCNDL,CNDL  
  real(r8) :: ConvectiveVaporFlux,ConvectiveHeatFluxMicP,PSISVL
  real(r8) :: TKLX
  real(r8) :: ConvectiveHeatFlxMacP,THeatDueWatFlow,HFLWS,THETA1,THETAL  
  integer  :: IFLGH
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

  IFLGH=0
  D4400: DO L=1,NL(NY,NX)
    N1=NX
    N2=NY
    N3=L
    !
    !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
    !
    D4320: DO N=NCN(N2,N1),3
      IF(N.EQ.iewstdir)THEN
        !west-east direction
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
      ELSEIF(N.EQ.insthdir)THEN
        !north-south direction
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
      ELSEIF(N.EQ.ivertdir)THEN
        !vertical direction
        IF(L.EQ.NL(NY,NX))THEN
          cycle
        ELSE
          N4=NX
          N5=NY
          N6=L+1
        ENDIF
      ENDIF
!
!     SKIP NON-EXISTENT DESTINATION SOIL LAYERS
!     identified by soil volume
      D1100: DO LL=N6,NL(NY,NX)
        IF(VLSoilPoreMicP(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
          N6=LL
          exit
        ENDIF
      ENDDO D1100

      IF(N3.EQ.NU(N2,N1))N6X(N2,N1)=N6
      !
      !     POROSITIES 'THETP*', WATER CONTENTS 'THETA*', AND POTENTIALS
      !     'PSIS*' FOR EACH GRID CELL
      !
      !     THETA1,THETAL=micropore water concn in source,destination cells
      !     THETY=hygroscopic water concentration
      !     POROS=soil porosity
      !     VLSoilPoreMicPI=soil volume excluding rock, macropore
      !
      IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
        IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4).AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN

          ! source layer
          call CalcSoilWatPotential(NY,NX,N1,N2,N3,PSISoilMatricPtmp(N3,N2,N1),THETA1)
          !dest
          call CalcSoilWatPotential(NY,NX,N4,N5,N6,PSISoilMatricPtmp(N6,N5,N4),THETAL)
          !
          !     ACCOUNT FOR WETTING FRONTS WHEN CALCULATING WATER CONTENTS,
          !     MATRIC WATER POTENTIALS AND HYDRAULIC CONDUCTIVITIES USED
          !     IN WATER FLUX CALCULATIONS
          !
          !     CND1,CNDL=hydraulic conductivities in source,destination cells
          !     FKSAT=reduction in soil surface Ksat from rainfall energy impact
          !     PSISM1=soil matric potential
          !     VLWatMicPX1=VLWatMicP1 accounting for wetting front
          !
          !     DARCY FLOW IF BOTH CELLS ARE SATURATED
          !     (CURRENT WATER POTENTIAL > AIR ENTRY WATER POTENTIAL)
          !
          call MicporeDarcyFlow(NY,NX,N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,&
            FKSAT(NY,NX),THeatDueWatFlow,PSISV1,PSISVL)          

!
          !     MACROPORE FLOW FROM POISEUILLE FLOW IF MACROPORES PRESENT
          !
          call MacporeFLow(NY,NX,M,N,N1,N2,N3,N4,N5,N6,ConvectiveHeatFlxMacP,IFLGH)

!         micropore flow
          call WaterVaporFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectiveVaporFlux,ConvectiveHeatFluxMicP)

          !
          !     FLWL=total water+vapor flux to destination
          !     WatXChange2WatTableX=total unsaturated water+vapor flux to destination
          !     THeatDueWatFlow=total convective heat flux from water+vapor flux
          !
          WatXChange2WatTable(N,N6,N5,N4)=WatXChange2WatTable(N,N6,N5,N4)+ConvectiveVaporFlux
          WatXChange2WatTableX(N,N6,N5,N4)=WatXChange2WatTableX(N,N6,N5,N4)+ConvectiveVaporFlux
          THeatDueWatFlow=THeatDueWatFlow+ConvectiveHeatFluxMicP
          HeatFlowi(N,N6,N5,N4)=THeatDueWatFlow+ConvectiveHeatFlxMacP
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
          DTKX=ABS(TKSoi1(N3,N2,N1)-TKSoi1(N6,N5,N4))*ppmc

          call CalcThermConduc(N1,N2,N3,DTKX,ThermCondSrc)

          call CalcThermConduc(N4,N5,N6,DTKX,ThermCondDst)

          ATCNDL=(2.0_r8*ThermCondSrc*ThermCondDst)/(ThermCondSrc*DLYR(N,N6,N5,N4)+ThermCondDst*DLYR(N,N3,N2,N1))

          call Solve4Heat(N,NY,NX,N1,N2,N3,N4,N5,N6,ATCNDL,ConvectiveHeatFluxMicP,HeatFlux2Ground(NY,NX))

          !
          !     TOTAL WATER, VAPOR AND HEAT FLUXES
          !
          !     FLW,FLWX,FLWH=total water flux through micropores,macropores
          !     HFLW=total heat flux
          !     FLWM=water flux used for solute flux calculations in trnsfr.f
          !
          WaterFlowSoiMicP(N,N6,N5,N4)=WaterFlowSoiMicP(N,N6,N5,N4)+WatXChange2WatTable(N,N6,N5,N4)
          FLWX(N,N6,N5,N4)=FLWX(N,N6,N5,N4)+WatXChange2WatTableX(N,N6,N5,N4)
          WaterFlowMacP(N,N6,N5,N4)=WaterFlowMacP(N,N6,N5,N4)+ConvectWaterFlowMacP(N,N6,N5,N4)
          HeatFlow(N,N6,N5,N4)=HeatFlow(N,N6,N5,N4)+HeatFlowi(N,N6,N5,N4)
          FLWM(M,N,N6,N5,N4)=WatXChange2WatTable(N,N6,N5,N4)

          IF(N.EQ.ivertdir)THEN
            !
            !     WATER FILM THICKNESS FOR CALCULATING GAS EXCHANGE IN TRNSFR.F
            !
            FILM(M,N6,N5,N4)=FilmThickness(PSISoilMatricPtmp(N6,N5,N4))
          ENDIF
        ELSEIF(N.NE.ivertdir)THEN
          WatXChange2WatTable(N,N6,N5,N4)=0.0_r8
          WatXChange2WatTableX(N,N6,N5,N4)=0.0_r8
          ConvectWaterFlowMacP(N,N6,N5,N4)=0.0_r8
          HeatFlowi(N,N6,N5,N4)=0.0_r8
          FLWM(M,N,N6,N5,N4)=0.0_r8
          WaterFlowMacPi(M,N,N6,N5,N4)=0.0_r8
        ENDIF

      ELSE
        IF(N.EQ.ivertdir)THEN
          WatXChange2WatTable(N,N3,N2,N1)=0.0_r8
          WatXChange2WatTableX(N,N3,N2,N1)=0.0_r8
          ConvectWaterFlowMacP(N,N3,N2,N1)=0.0_r8
          HeatFlowi(N,N3,N2,N1)=0.0_r8
          WaterFlowMacPi(M,N,N3,N2,N1)=0.0_r8
          WaterFlowMacPi(M,N,N3,N2,N1)=0.0_r8
        ELSE
          WatXChange2WatTable(N,N6,N5,N4)=0.0_r8
          WatXChange2WatTableX(N,N6,N5,N4)=0.0_r8
          ConvectWaterFlowMacP(N,N6,N5,N4)=0.0_r8
          HeatFlowi(N,N6,N5,N4)=0.0_r8
          FLWM(M,N,N6,N5,N4)=0.0_r8
          WaterFlowMacPi(M,N,N6,N5,N4)=0.0_r8
        ENDIF
      ENDIF
    ENDDO D4320
  ENDDO D4400
  ENDDO
  ENDDO  
  end subroutine Subsurface3DFlowMit
!------------------------------------------------------------------------------------------

  subroutine LateralWatHeatExchMit(M,NHW,NHE,NVN,NVS,FKSATS)
  !
  !Description
  ! boundary flow involes exchange with external water table, and through tile drainage
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: FKSATS(JY,JX)
  integer :: NY,NX
  integer :: L,LL
  integer :: N,NN,N1,N2,N3,N4,N5,N4B,N5B,N6
  integer :: M1,M2,M3,M4,M5,M6,K1,KL
  integer :: IFLGDH,IFLGU,IFLGUH,IFLGD
  real(r8) :: CND1,VLWatMicP1X
  real(r8) :: FLWT
  real(r8) :: RCHQF,RCHGFU,RCHGFT
  real(r8) :: DPTHH,CNDL
  real(r8) :: FINHX,THETAX  
  real(r8) :: VOLP2,VOLPX2,VOLPH2
  real(r8) :: XN,THETA1
  real(r8) :: VOLP1X,VLWatMacP1X,VOLPH1X

!     begin_execution
!     VOLP2,VOLPH2=air-filled porosity in micropores,macropores
  D9595: DO  NX=NHW,NHE
    D9590: DO  NY=NVN,NVS
      D9585: DO L=NUM(NY,NX),NL(NY,NX)
        VOLP2=VLMicP1(L,NY,NX)-VLWatMicP1(L,NY,NX)-VLiceMicP1(L,NY,NX)
        VOLPX2=VOLP2
        VOLPH2=VLMacP1(L,NY,NX)-VLWatMacP1(L,NY,NX)-VLiceMacP1(L,NY,NX)
!
        call Config4WaterTableDrain(L,NY,NX,IFLGU,IFLGUH,DPTHH)

!
!     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO TILE DRAIN
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
        D9580: DO N=NCN(NY,NX),3
          D9575: DO NN=1,2
            IF(N.EQ.iewstdir)THEN
! along the W-E direction
              N4=NX+1
              N5=NY
              N4B=NX-1
              N5B=NY
              N6=L
              IF(NN.EQ.1)THEN
                !east boundary
                IF(NX.EQ.NHE)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX+1
                  M5=NY
                  M6=L
                  XN=-1.0_r8
                  RCHQF=RCHQE(M2,M1)
                  RCHGFU=RCHGEU(M2,M1)
                  RCHGFT=RCHGET(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                !west boundary
                IF(NX.EQ.NHW)THEN
                  M1=NX+1
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0_r8
                  RCHQF=RCHQW(M5,M4)
                  RCHGFU=RCHGWU(M5,M4)
                  RCHGFT=RCHGWT(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.insthdir)THEN
! along the N-S direction
              N4=NX
              N5=NY+1
              N4B=NX
              N5B=NY-1
              N6=L
              IF(NN.EQ.1)THEN
                !south boundary
                IF(NY.EQ.NVS)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY+1
                  M6=L
                  XN=-1.0_r8
                  RCHQF=RCHQS(M2,M1)
                  RCHGFU=RCHGSU(M2,M1)
                  RCHGFT=RCHGST(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                !north boundary
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY+1
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0_r8
                  RCHQF=RCHQN(M5,M4)
                  RCHGFU=RCHGNU(M5,M4)
                  RCHGFT=RCHGNT(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.ivertdir)THEN
! in the vertical direction
              N4=NX
              N5=NY
              N6=L+1
              IF(NN.EQ.1)THEN
                !bottom
                IF(L.EQ.NL(NY,NX))THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L+1
                  XN=-1.0_r8
                  RCHGFU=RCHGD(M2,M1)
                  RCHGFT=1.0_r8
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
!     IRCHG,RCHQ*=runoff boundary flags
!           top soil layer and surface soil layer is active
            IF(L.EQ.NUM(N2,N1).AND.N.NE.ivertdir.AND. &
              (CumDepth2LayerBottom(NU(N2,N1)-1,N2,N1).LE.CumSoilDeptht0(N2,N1) &
              .OR.SoiBulkDensity(NUI(N2,N1),N2,N1).GT.ZERO))THEN
!not in vertical direction
              IF(IRCHG(NN,N,N2,N1).EQ.0.OR.isclose(RCHQF,0._r8).OR.ABS(QRM(M,N2,N1)).LT.ZEROS(N2,N1))THEN
                !no runoff
                QR1(N,NN,M5,M4)=0.0_r8
                HQR1(N,NN,M5,M4)=0.0_r8
              ELSE

                call SurfaceRunoff(M,N,NN,N1,N2,M4,M5,RCHQF,XN)
!
        !     BOUNDARY SNOW FLUX
        !
        !     QS1,QW1,QI1=snow,water,ice transfer
        !     HQS1=convective heat transfer from snow,water,ice transfer
        !     QS,QW,QI=cumulative hourly snow,water,ice transfer
        !     HQS=cumulative hourly convective heat transfer from snow,water,ice transfer
!
                IF(NN.EQ.1)THEN
                  QS1(N,M5,M4)=0.0_r8
                  QW1(N,M5,M4)=0.0_r8
                  QI1(N,M5,M4)=0.0_r8
                  HQS1(N,M5,M4)=0.0_r8
                  QSM(M,N,M5,M4)=QS1(N,M5,M4)
                ENDIF
              ENDIF
            ELSE
              IF(N.NE.ivertdir)THEN
                QR1(N,NN,M5,M4)=0.0_r8
                HQR1(N,NN,M5,M4)=0.0_r8
              ENDIF
            ENDIF
!
          ! BOUNDARY SUBSURFACE WATER AND HEAT TRANSFER DEPENDING
          ! ON LEVEL OF WATER TABLE
!
            IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
              IF(NCN(N2,N1).NE.ivertdir.OR.N.EQ.ivertdir)THEN
              !including lateral connection or woking on vertical direction
!
              ! IF NO WATER TABLE
              !
              ! IDTBL=water table flag
              ! THETA1,THETAX=water content ahead,behind wetting front
              ! K1,KL=pore water class ahead,behind wetting front
              ! CND1,CNDL=hydraulic conductivity ahead,behind wetting front
              ! FKSAT=reduction in soil surface Ksat from rainfall energy impact
              ! FLWL,WatXChange2WatTableX=lower boundary micropore water flux
              ! ConvectWaterFlowMacP=lower boundary macropore water flux
              ! HFLWL=convective heat from lower boundary water flux
              ! XH,XN,dts_HeatWatTP=rate constant,direction indicator,time step
              ! SLOPE=sin(vertical slope)=1
              ! RCHG*=boundary flags
!
              IF(IDTBL(N2,N1).EQ.0.OR.N.EQ.ivertdir)THEN
                THETA1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1),safe_adb(VLWatMicP1(N3,N2,N1),VLSoilMicP(N3,N2,N1))))
                THETAX=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1),safe_adb(VLWatMicPX1(N3,N2,N1),VLSoilMicP(N3,N2,N1))))
                K1=MAX(1,MIN(100,INT(100.0_r8*(POROS(N3,N2,N1)-THETA1)/POROS(N3,N2,N1))+1))
                KL=MAX(1,MIN(100,INT(100.0_r8*(POROS(N3,N2,N1)-THETAX)/POROS(N3,N2,N1))+1))
                IF(N3.EQ.NUM(NY,NX))THEN
                  CND1=HydroCond3D(N,K1,N3,N2,N1)*FKSATS(NY,NX)
                ELSE
                  CND1=HydroCond3D(N,K1,N3,N2,N1)
                ENDIF
                
                WatXChange2WatTable(N,M6,M5,M4)=AMIN1(VLWatMicP1(N3,N2,N1)*dts_wat, &
                  XN*mGravAcceleration*(-ABS(SLOPE(N,N2,N1)))*CND1*AREA(3,N3,N2,N1)) &
                  *RCHGFU*RCHGFT*dts_HeatWatTP

                if(abs(WatXChange2WatTable(N,M6,M5,M4))>1.e20)then
                  write(*,*)'VLWatMicP1(N3,N2,N1)*dts_wat=',VLWatMicP1(N3,N2,N1),dts_wat
                  write(*,*)'XN=',XN,CND1
                  write(*,*)'RCHGFU*RCHGFT*dts_HeatWatTP=',RCHGFU,RCHGFT,dts_HeatWatTP
                  write(*,*)'at line',__LINE__
                  call endrun(trim(mod_filename)//'at line',__LINE__)
                endif
                WatXChange2WatTableX(N,M6,M5,M4)=WatXChange2WatTable(N,M6,M5,M4)
                ConvectWaterFlowMacP(N,M6,M5,M4)=AMIN1(VLWatMacP1(L,NY,NX)*dts_wat &
                 ,XN*mGravAcceleration*(-ABS(SLOPE(N,N2,N1)))*HydroCondMacP1(L,NY,NX)*AREA(3,N3,N2,N1)) &
                 *RCHGFU*RCHGFT*dts_HeatWatTP
                HeatFlowi(N,M6,M5,M4)=cpw*TKSoi1(N3,N2,N1)*(WatXChange2WatTable(N,M6,M5,M4)+ConvectWaterFlowMacP(N,M6,M5,M4))
              ELSE
!
                CALL WaterTBLDrain(N,N1,N2,N3,M4,M5,M6,IFLGU,IFLGUH,RCHGFU,RCHGFT,DPTHH,XN)

                call TileDrain(N,N1,N2,N3,M4,M5,M6,IFLGD,IFLGDH,RCHGFT,RCHGFU,DPTHH,XN)

                call SubSufRecharge(NY,NX,N,N1,N2,N3,M4,M5,M6,DPTHH,RCHGFU,RCHGFT,XN,VOLP2,VOLPX2,VOLPH2)

              ENDIF
!
        !     SUBSURFACE HEAT SOURCE/SINK
        !
        !     HFLWL=heat flux across lower boundary
        !     TK1=lower boundary soil temperature
        !     TKSD=deep source/sink temperature from geothermal flux
        !     TCNDG=thermal conductivity below lower boundary
        !     SoilHeatSrcDepth,CDPTH=depth of thermal sink/source, lower boundary
        !     IETYP=Koppen climate zone
              IF(N.EQ.ivertdir.AND.IETYP(N2,N1).NE.-2)THEN
                HeatFlowi(N,M6,M5,M4)=HeatFlowi(N,M6,M5,M4)+(TKSoi1(N3,N2,N1)-TKSD(N2,N1))* &
                  TCNDG/(SoilHeatSrcDepth(N2,N1)-CumDepth2LayerBottom(N3,N2,N1)) &
                  *AREA(N,N3,N2,N1)*dts_HeatWatTP
              ENDIF
              WaterFlowSoiMicP(N,M6,M5,M4)=WaterFlowSoiMicP(N,M6,M5,M4)+WatXChange2WatTable(N,M6,M5,M4)
              FLWX(N,M6,M5,M4)=FLWX(N,M6,M5,M4)+WatXChange2WatTableX(N,M6,M5,M4)
              WaterFlowMacP(N,M6,M5,M4)=WaterFlowMacP(N,M6,M5,M4)+ConvectWaterFlowMacP(N,M6,M5,M4)
              HeatFlow(N,M6,M5,M4)=HeatFlow(N,M6,M5,M4)+HeatFlowi(N,M6,M5,M4)
              FLWM(M,N,M6,M5,M4)=WatXChange2WatTable(N,M6,M5,M4)
              WaterFlowMacPi(M,N,M6,M5,M4)=ConvectWaterFlowMacP(N,M6,M5,M4)
            ENDIF
          ELSE
            WatXChange2WatTable(N,M6,M5,M4)=0.0_r8
            WatXChange2WatTableX(N,M6,M5,M4)=0.0_r8
            ConvectWaterFlowMacP(N,M6,M5,M4)=0.0_r8
            HeatFlowi(N,M6,M5,M4)=0.0_r8
            FLWM(M,N,M6,M5,M4)=0.0_r8
            WaterFlowMacPi(M,N,M6,M5,M4)=0.0_r8
          ENDIF
        ENDDO D9575
    !
    !     NET WATER AND HEAT FLUXES IN RUNOFF AND SNOW DRIFT
    !
    !     TQR1,THQR1=net runoff,convective heat from runoff
    !     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
    !     QR1,HQR1=runoff, convective heat from runoff
    !     QS1,QW1,QI1=snow,water,ice transfer
    !     HQS1=convective heat transfer from snow,water,ice transfer
!
        IF(L.EQ.NUM(N2,N1).AND.N.NE.ivertdir)THEN
          !top layer snow redistribution
          call SumSnowRoffDrift(M,N,N1,N2,N4,N5,N4B,N5B)
        ENDIF
!
        !     NET WATER AND HEAT FLUXES THROUGH SOIL AND SNOWPACK
        !
        !     TFLWL,THFLWL=net water micropore,macropore flux
        !     THFLWL=net convective+conductive heat flux
        !     FLWL =micropore water,heat flux
        !     ConvectWaterFlowMacP=macropore water,heat flux
        !     HFLWL=soil heat flux
!
        IF(NCN(N2,N1).NE.ivertdir.OR.N.EQ.ivertdir)THEN
          D1200: DO LL=N6,NL(N5,N4)
            IF(VLSoilPoreMicP(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
              N6=LL
              exit
            ENDIF
          ENDDO D1200
          !exchange with water table if micropore is non-zero
          IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            TWatCharge2MicP(N3,N2,N1)=TWatCharge2MicP(N3,N2,N1)+WatXChange2WatTable(N,N3,N2,N1)-WatXChange2WatTable(N,N6,N5,N4)
            TWatXChange2WatTableX(N3,N2,N1)=TWatXChange2WatTableX(N3,N2,N1)+WatXChange2WatTableX(N,N3,N2,N1)&
              -WatXChange2WatTableX(N,N6,N5,N4)
            TConvectWaterFlowMacP(N3,N2,N1)=TConvectWaterFlowMacP(N3,N2,N1)+ConvectWaterFlowMacP(N,N3,N2,N1) &
              -ConvectWaterFlowMacP(N,N6,N5,N4)
            THeatFlowi(N3,N2,N1)=THeatFlowi(N3,N2,N1)+HeatFlowi(N,N3,N2,N1)-HeatFlowi(N,N6,N5,N4)
          ELSE
            TWatCharge2MicP(N3,N2,N1)=0.0_r8
            TWatXChange2WatTableX(N3,N2,N1)=0.0_r8
            TConvectWaterFlowMacP(N3,N2,N1)=0.0_r8
            THeatFlowi(N3,N2,N1)=0.0_r8
          ENDIF
        ENDIF
      ENDDO D9580
!
!     INFILTRATION OF WATER FROM MACROPORES INTO MICROPORES
!
!     VLWatMacP1=macropore volume
!     FINHX,FINHL=macro-micropore transfer unltd,ltd by water,air volume
!     FINHM=macro-micropore transfer for use in trnsfr.f
!     HydroCond3D=hydraulic conductivity
!     PSISE,PSISoilAirEntry=air entry,matric water potentials
!     PHOL,MacPRadius=path length between,radius of macropores from hour1.f
!     dts_HeatWatTP=time step
!     VLWatMicP1X,VOLP1X=current micropore water,air volume
!     VLWatMacP1X,VOLPH1X=current macropore water,air volume
!
      IF(VLWatMacP1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
        !south-north direction, is it true?
        FINHX=TwoPiCON*HydroCond3D(2,1,N3,N2,N1)*AREA(3,N3,N2,N1) &
          *(PSISE(N3,N2,N1)-PSISoilMatricPtmp(N3,N2,N1)) &
          /LOG(PathLenMacP(N3,N2,N1)/MacPRadius(N3,N2,N1))*dts_HeatWatTP
        VLWatMicP1X=VLWatMicP1(N3,N2,N1)+TWatCharge2MicP(N3,N2,N1)+FWatIrrigate2MicP1(N3,N2,N1)
        VOLP1X=AZMAX1(VLMicP1(N3,N2,N1)-VLWatMicP1X-VLiceMicP1(N3,N2,N1))
        VLWatMacP1X=VLWatMacP1(N3,N2,N1)+TConvectWaterFlowMacP(N3,N2,N1)
        VOLPH1X=AZMAX1(VLMacP1(N3,N2,N1)-VLWatMacP1X-VLiceMacP1(N3,N2,N1))
        IF(FINHX.GT.0.0_r8)THEN
          FWatExMacP2MicP(N3,N2,N1)=AZMAX1(AMIN1(FINHX,VLWatMacP1X,VOLP1X))
        ELSE
          FWatExMacP2MicP(N3,N2,N1)=AZMIN1(AMAX1(FINHX,-VOLPH1X,-VLWatMicP1X))
        ENDIF
        FINHM(M,N3,N2,N1)=FWatExMacP2MicP(N3,N2,N1)
        FINH(N3,N2,N1)=FINH(N3,N2,N1)+FWatExMacP2MicP(N3,N2,N1)
      ELSE
        FWatExMacP2MicP(N3,N2,N1)=0.0_r8
        FINHM(M,N3,N2,N1)=0.0_r8
      ENDIF

      call FreezeThawMit(NY,NX,L,N1,N2,N3)
!
!     DISSIPATE WETTING FRONT
!
!     VLWatMicP1=soil micropore water content
!     VLWatMicPX1=soil micropore water content behind wetting front
!     FLWVL=water flux from wetted to drier soil
!
      ENDDO D9585
    ENDDO D9590
  ENDDO D9595

  end subroutine LateralWatHeatExchMit
!------------------------------------------------------------------------------------------

  subroutine UpdateStateFluxAtM(M,NHW,NHE,NVN,NVS)
  !
  !Description
  ! Early exit of the watsub solver
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX
  integer :: L,NUX,LL,Ls

  real(r8) :: tk1pres,tk1l
  real(r8) :: ENGY1,VLTSoiPore,VHXX
  real(r8) :: TKXX
  
! begin_execution
  D9795: DO NX=NHW,NHE
    D9790: DO NY=NVN,NVS

      call UpdateSurfaceAtM(M,NY,NX)

      !
      ! SOIL LAYER WATER, ICE AND TEMPERATURE
      !
      ! VLWatMicP1,VLiceMicP1=micropore water,ice volume
      ! VLWatMicPX1=micropore water volume behind wetting front
      ! VLWatMacP1,VLiceMacP1=macropore water,ice volume
      ! TFLWL=net water flux
      ! FINHL=micropore-macropore flux
      ! TMLiceThawMicP,TMLiceThawMacP=total accumulated freeze-thaw in micropores,macropores
      ! FLU1=subsurface water input
      ! DENSICE=ice density
      ! VOLA1,VLMacP1=micropore,macropore volume
      ! VOLP1,VOLPH1=micropore,macropore air volume
      ! VOLWM,VLWatMacPM,VLsoiAirPM,FLPM=micropore,macropore water volume, air volume
      ! and change in air volume for use in trnsfr.f
      ! THETWX,FracSoiPAsIce,FracSoiPAsAir,FracSoilAsAirt=bulk water,ice,air concn,air-filled porosity
      ! THETPM=air concentration for use in trnsfr.f
      ! FMAC,SoilFracAsMicP=macropore,micropore fraction
      ! HydroCondMacP1=maropore hydraulic conductivity
      ! VLHeatCapacity,VHCM=volumetric heat capacities of total volume, solid
      ! VLHeatCapacityA,VLHeatCapacityB=volumetric heat capacities of soil+micropore,macropore
      ! TK1=soil temperature
      !

      D9785: DO L=NUM(NY,NX),NL(NY,NX)
        IF(VGeomLayer(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VLWatMicP1(L,NY,NX)=VLWatMicP1(L,NY,NX)+TWatCharge2MicP(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
            +TMLiceThawMicP(L,NY,NX)+FWatIrrigate2MicP1(L,NY,NX)

          if(abs(VLWatMicP1(L,NY,NX))>1.e20_r8)then
            write(*,*)'VLWatMicP1(L,NY,NX)=',VLWatMicP1(L,NY,NX),L
            write(*,*)'TFLWL=',TWatCharge2MicP(L,NY,NX)
            write(*,*)'FINHL=',FWatExMacP2MicP(L,NY,NX)
            write(*,*)'TMLiceThawMicP=',TMLiceThawMicP(L,NY,NX)+FWatIrrigate2MicP1(L,NY,NX)
            write(*,*)'at line',__LINE__
            call endrun(trim(mod_filename)//'at line',__LINE__)
          endif
          VLWatMicPX1(L,NY,NX)=VLWatMicPX1(L,NY,NX)+TWatXChange2WatTableX(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
             +TMLiceThawMicP(L,NY,NX)+FWatIrrigate2MicP1(L,NY,NX)
          VLWatMicPX1(L,NY,NX)=AMIN1(VLWatMicP1(L,NY,NX),VLWatMicPX1(L,NY,NX))
          VLiceMicP1(L,NY,NX)=VLiceMicP1(L,NY,NX)-TMLiceThawMicP(L,NY,NX)/DENSICE
          VLWatMacP1(L,NY,NX)=VLWatMacP1(L,NY,NX)+TConvectWaterFlowMacP(L,NY,NX)-FWatExMacP2MicP(L,NY,NX)+TMLiceThawMacP(L,NY,NX)
          VLiceMacP1(L,NY,NX)=VLiceMacP1(L,NY,NX)-TMLiceThawMacP(L,NY,NX)/DENSICE

          IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
 ! air-filled space
            VLairMicP(L,NY,NX)=VLMicP1(L,NY,NX)-VLWatMicP1(L,NY,NX)-VLiceMicP1(L,NY,NX)
            VLairMicPc(L,NY,NX)=AZMAX1(VLairMicP(L,NY,NX))
            VLairMacP(L,NY,NX)=VLMacP1(L,NY,NX)-VLWatMacP1(L,NY,NX)-VLiceMacP1(L,NY,NX)
            VLairMacPc(L,NY,NX)=AZMAX1(VLairMacP(L,NY,NX))
            VLMacP1(L,NY,NX)=AZMAX1(VLMacP(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX) &
              *(safe_adb(VLWatMicP1(L,NY,NX),VLSoilMicP(L,NY,NX))-WiltPoint(L,NY,NX))*VGeomLayer(L,NY,NX))
          ELSE
            VLairMicP(L,NY,NX)=0.0_r8
            VLairMicPc(L,NY,NX)=0.0_r8
            VLairMacP(L,NY,NX)=0.0_r8
            VLairMacPc(L,NY,NX)=0.0_r8
            VLMicP1(L,NY,NX)=VLWatMicP1(L,NY,NX)+VLiceMicP1(L,NY,NX)
            VLMacP1(L,NY,NX)=0.0_r8
          ENDIF

          !record intermediate variables for bgc calculation
          VLWatMicPM(M+1,L,NY,NX)=VLWatMicP1(L,NY,NX)
          VLWatMacPM(M+1,L,NY,NX)=VLWatMacP1(L,NY,NX)
          VLsoiAirPM(M+1,L,NY,NX)=VLairMicPc(L,NY,NX)+VLairMacPc(L,NY,NX) &
            +THETPI*(VLiceMicP1(L,NY,NX)+VLiceMacP1(L,NY,NX))
          !change in soil air volume
          FLPM(M,L,NY,NX)=VLsoiAirPM(M,L,NY,NX)-VLsoiAirPM(M+1,L,NY,NX)
          VLTSoiPore=VLSoilMicP(L,NY,NX)+VLMacP1(L,NY,NX)
          FracSoiPAsWat(L,NY,NX)=AZMAX1t((VLWatMicP1(L,NY,NX)+VLWatMacP1(L,NY,NX))/VLTSoiPore)
          FracSoiPAsIce(L,NY,NX)=AZMAX1t((VLiceMicP1(L,NY,NX)+VLiceMacP1(L,NY,NX))/VLTSoiPore)
          FracSoiPAsAir(L,NY,NX)=AZMAX1t((VLairMicPc(L,NY,NX)+VLairMacPc(L,NY,NX))/VLTSoiPore)
          THETPM(M+1,L,NY,NX)=FracSoiPAsAir(L,NY,NX)
          IF(VLMicP1(L,NY,NX)+VLMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            FracSoilAsAirt(L,NY,NX)=AZMAX1((VLairMicPc(L,NY,NX)+VLairMacPc(L,NY,NX))/(VLMicP1(L,NY,NX)+VLMacP1(L,NY,NX)))
          ELSE
            FracSoilAsAirt(L,NY,NX)=0.0_r8
          ENDIF
          IF(VLMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            SoilFracAsMacP1(L,NY,NX)=SoilFracAsMacP(L,NY,NX)*VLMacP1(L,NY,NX)/VLMacP(L,NY,NX)
            HydroCondMacP1(L,NY,NX)=HydroCondMacP(L,NY,NX)*(VLMacP1(L,NY,NX)/VLMacP(L,NY,NX))**2
          ELSE
            SoilFracAsMacP1(L,NY,NX)=0.0_r8
            HydroCondMacP1(L,NY,NX)=0.0_r8
          ENDIF
          SoilFracAsMicP(L,NY,NX)=1.0_r8-SoilFracAsMacP1(L,NY,NX)
          TKXX=TKSoi1(L,NY,NX)
          VHXX=VLHeatCapacity(L,NY,NX)
          ENGY1=VLHeatCapacity(L,NY,NX)*TKSoi1(L,NY,NX)
          if(TKSoi1(L,NY,NX)>1.e3_r8.or.TKSoi1(L,NY,NX)<0._r8)then
            write(*,*)'L=',L,NY,NX,NUM(NY,NX)
            write(*,*)'SoiBulkDensity(L,NY,NX)=',SoiBulkDensity(L,NY,NX)
            write(*,*)'VLHeatCapacity(L,NY,NX),TKSoi1(L,NY,NX)',L,VLHeatCapacity(L,NY,NX),TKSoi1(L,NY,NX)
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          VLHeatCapacityA(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+cpw*VLWatMicP1(L,NY,NX)+cpi*VLiceMicP1(L,NY,NX)
          VLHeatCapacityB(L,NY,NX)=cpw*VLWatMacP1(L,NY,NX)+cpi*VLiceMacP1(L,NY,NX)
          VLHeatCapacity(L,NY,NX)=VLHeatCapacityA(L,NY,NX)+VLHeatCapacityB(L,NY,NX)
          !
          !         BEGIN ARTIFICIAL SOIL WARMING
          !
          !         THFLWL=THFLWL incremented for soil warming
          !         TKSZ=temperature used to calculate additional heat flux for warming
          !
          !         IF(NX.EQ.3.AND.NY.EQ.2.AND.L.GT.NUM(NY,NX)
          !           3.AND.L.LE.17.AND.I.GE.152.AND.I.LE.304)THEN
          !           THeatFlowi(L,NY,NX)=THeatFlowi(L,NY,NX)
          !             2+(TKSZ(I,J,L)-TKSoi1(L,NY,NX))*VLHeatCapacity(L,NY,NX)*dts_HeatWatTP
          !             WRITE(*,3379)'TKSZ',I,J,M,NX,NY,L,TKSZ(I,J,L)
          !               2,TKSoi1(L,NY,NX),VLHeatCapacity(L,NY,NX),THeatFlowi(L,NY,NX)
          !3379  FORMAT(A8,6I4,12E12.4)
          !           ENDIF
          !
          !           END ARTIFICIAL SOIL WARMING
!
          IF(VLHeatCapacity(L,NY,NX).GT.ZEROS(NY,NX))THEN
            tk1l=TKSoi1(L,NY,NX)
            TKSoi1(L,NY,NX)=(ENGY1+THeatFlowi(L,NY,NX)+TLSoiPIceHeatFlxFrez1(L,NY,NX)+&
              HeatIrrigation1(L,NY,NX))/VLHeatCapacity(L,NY,NX)
            if(abs(TKSoi1(L,NY,NX)/tk1l-1._r8)>0.025_r8)then
              TKSoi1(L,NY,NX)=tk1l
            endif
          ELSEIF(L.EQ.1)THEN
            TKSoi1(L,NY,NX)=TairK(NY,NX)
          ELSE
            TKSoi1(L,NY,NX)=TKSoi1(L-1,NY,NX)
          ENDIF
        ELSE
          VLWatMicPM(M+1,L,NY,NX)=0.0_r8
          VLWatMacPM(M+1,L,NY,NX)=0.0_r8
          VLsoiAirPM(M+1,L,NY,NX)=0.0_r8
          FLPM(M,L,NY,NX)=VLsoiAirPM(M,L,NY,NX)
          THETPM(M+1,L,NY,NX)=0.0_r8
        ENDIF

      ENDDO D9785
      !
      !       RESET SURFACE LAYER NUMBER AND TRANSFER ALL WATER TO SOIL SURFACE LAYER
      !       IF LAKE SURFACE LAYER IS LOST TO EVAPORATION
      !
      !       NUM=new surface layer number after complete lake evaporation
      !       LakeSurfFlow,FLWHNU,LakeSurfHeatFlux=lake surface water flux, heat flux if lake surface disappears
!
      IF(SoiBulkDensity(NUM(NY,NX),NY,NX).LE.ZERO.AND.VLHeatCapacity(NUM(NY,NX),NY,NX).LE.VHCPNX(NY,NX))THEN
        NUX=NUM(NY,NX)
        DO  LL=NUX+1,NL(NY,NX)
          IF(VLSoilPoreMicP(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
            NUM(NY,NX)=LL
            FLWNX(NY,NX)=WaterFlowSoiMicP(3,NUM(NY,NX),NY,NX)
            FLWXNX(NY,NX)=FLWX(3,NUM(NY,NX),NY,NX)
            FLWHNX(NY,NX)=WaterFlowMacP(3,NUM(NY,NX),NY,NX)
            HFLWNX(NY,NX)=HeatFlow(3,NUM(NY,NX),NY,NX)
            exit
          ENDIF
        ENDDO
      ENDIF
    ENDDO D9790
  ENDDO D9795

  end subroutine UpdateStateFluxAtM

!------------------------------------------------------------------------------------------

  subroutine UpdateFluxAtExit(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
! begin_execution
  D9695: DO NX=NHW,NHE
    D9690: DO NY=NVN,NVS
      IF(NUM(NY,NX).EQ.NU(NY,NX))THEN
        LakeSurfFlow(NY,NX)=WaterFlowSoiMicP(3,N6X(NY,NX),NY,NX)
        FLWXNU(NY,NX)=FLWX(3,N6X(NY,NX),NY,NX)
        FLWHNU(NY,NX)=WaterFlowMacP(3,N6X(NY,NX),NY,NX)
        LakeSurfHeatFlux(NY,NX)=HeatFlow(3,N6X(NY,NX),NY,NX)
      ELSE
        LakeSurfFlow(NY,NX)=FLWNX(NY,NX)
        FLWXNU(NY,NX)=FLWXNX(NY,NX)
        FLWHNU(NY,NX)=FLWHNX(NY,NX)
        LakeSurfHeatFlux(NY,NX)=HFLWNX(NY,NX)
      ENDIF
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
!     IDTBL=water table flag
!     DPTH,DTBLY=depth to layer midpoint, artificial water table
!     PSISM1,PSISE=soil,air entry matric potential
!     DTBLYX=equilibrium water potential with artificial water table
!     IFLGD=micropore discharge flag to artificial water table
!
  IF(IDTBL(NY,NX).GE.3.AND.DPTH(L,NY,NX).LT.DTBLY(NY,NX))THEN
    IF(PSISM1(L,NY,NX).GT.mGravAcceleration*(DPTH(L,NY,NX)-DTBLY(NY,NX)))THEN
      IFLGD=0
      IF(L.LT.NL(NY,NX))THEN
        D9568: DO  LL=L+1,NL(NY,NX)
          DTBLYX=DTBLY(NY,NX)+PSISE(LL,NY,NX)/mGravAcceleration
          IF(DPTH(LL,NY,NX).LT.DTBLYX)THEN
            IF((PSISM1(LL,NY,NX).LE.mGravAcceleration*(DPTH(LL,NY,NX)-DTBLYX) &
              .AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
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
  IF(VLMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    DPTHH=CumDepth2LayerBottom(L,NY,NX)-(VLWatMacP1(L,NY,NX)+VLiceMacP1(L,NY,NX))/VLMacP1(L,NY,NX)*DLYR(3,L,NY,NX)
  ELSE
    DPTHH=CumDepth2LayerBottom(L,NY,NX)
  ENDIF

  IF(IDTBL(NY,NX).GE.3.AND.DPTHH.LT.DTBLY(NY,NX).AND.VLWatMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
! artificial water table, e.g. tile drainage
    IFLGDH=0
    IF(L.LT.NL(NY,NX))THEN
      D9569: DO  LL=L+1,NL(NY,NX)
        IF(DPTH(LL,NY,NX).LT.DTBLY(NY,NX))THEN
! the layer is above tile drain water table
          IF(VLMacP1(LL,NY,NX).LE.ZEROS(NY,NX))THEN
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

  real(r8) :: DTBLXX
  integer :: LL

!     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO WATER TABLE
!
!     IDTBL=water table flag
!     DPTH,DTBLX=depth to layer midpoint,natural water table
!     PSISM1(<0),PSISE(<0)=matric,air entry water potential [MPa], from water height to MPa, 1.e-6*9.8*1.e3*H
!     DTBLXX=equilibrium water potential with natural water table
!     DPTHA=active layer depth
!     IFLGU=micropore discharge flag to natural water table
!  total water potential psi+rho*grav*h

  IF(IDTBL(NY,NX).NE.0.AND.DPTH(L,NY,NX).LT.DTBLX(NY,NX))THEN
    !the layer mid-depth is lower than water table
    IF(PSISM1(L,NY,NX).GT.mGravAcceleration*(DPTH(L,NY,NX)-DTBLX(NY,NX)))THEN
      IFLGU=0
      D9565: DO LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
        !water level difference
        DTBLXX=DTBLX(NY,NX)+PSISE(LL,NY,NX)/mGravAcceleration  
        IF(DPTH(LL,NY,NX).LT.DTBLXX)THEN
          IF((PSISM1(LL,NY,NX).LE.mGravAcceleration*(DPTH(LL,NY,NX)-DTBLXX) &
            .AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
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
  IF(VLMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    DPTHH=CumDepth2LayerBottom(L,NY,NX)-(VLWatMacP1(L,NY,NX)+VLiceMacP1(L,NY,NX))/VLMacP1(L,NY,NX)*DLYR(3,L,NY,NX)
  ELSE
    DPTHH=CumDepth2LayerBottom(L,NY,NX)
  ENDIF
  IF(IDTBL(NY,NX).NE.0.AND.DPTHH.LT.DTBLX(NY,NX).AND.VLWatMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    !active water table
    IFLGUH=0
!     DO 9566 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
!     IF(DPTH(LL,NY,NX).LT.DTBLX(NY,NX))THEN
!     IF(VLMacP1(LL,NY,NX).LE.ZEROS(NY,NX))THEN
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
  real(r8) :: THETWT,TFND1,THETWA
  real(r8) :: VLSoiPorAV,VLWatSoi,scalar
  real(r8) :: THETWH,Z3S

  D9885: DO L=NUM(NY,NX),NL(NY,NX)
    TMLiceThawMicP(L,NY,NX)=0.0_r8
    TMLiceThawMacP(L,NY,NX)=0.0_r8
    TLSoiPIceHeatFlxFrez1(L,NY,NX)=0.0_r8
    TWatCharge2MicP(L,NY,NX)=0.0_r8
    TWatXChange2WatTableX(L,NY,NX)=0.0_r8
    TConvectWaterFlowMacP(L,NY,NX)=0.0_r8
    THeatFlowi(L,NY,NX)=0.0_r8

    !VOLW2 will be updated in soil surface model
    !VOLW2(L,NY,NX)=VLWatMicP1(L,NY,NX)
!
!   GAS EXCHANGE COEFFICIENTS SOIL LAYERS
!
!   VOLA1,VLiceMicP1,VLWatMicP1=total,ice-,water-filled microporosity
!   VLMacP1,VLiceMacP1,VLWatMacP1=total,ice-,water-filled macroporosity
!   VLsoiAirPM=air-filled porosity
!   TFND1=temperature effect on gas diffusivity
!   DFGS=rate constant for air-water gas exchange
!   Z1S,Z2SW,Z2SD,Z3SX=parameters for soil air-water gas transfers
!   XNPD=time step for gas transfer calculations
!   TORT,TortMacPM=tortuosity for aqueous diffn in micropores,macropres
!
    VLWatSoi=VLWatMicP1(L,NY,NX)+VLWatMacP1(L,NY,NX)
    VLSoiPorAV=VLMicP1(L,NY,NX)+VLMacP1(L,NY,NX)-VLiceMicP1(L,NY,NX)-VLiceMacP1(L,NY,NX)
    IF(VLSoiPorAV.GT.ZEROS2(NY,NX).AND.VLsoiAirPM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWA=AZMAX1(AMIN1(1.0_r8,VLWatSoi/VLSoiPorAV))
      TFND1=TEFAQUDIF(TKSoi1(0,NY,NX))
      Z3S=FieldCapacity(L,NY,NX)/POROS(L,NY,NX)
      scalar=TFND1*XNPD
      DFGS(M,L,NY,NX)=fDFGS(scalar,THETWA,Z3S)
    ELSE
      DFGS(M,L,NY,NX)=0.0_r8
    ENDIF

    IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
      THETWT=safe_adb(VLWatMicPM(M,L,NY,NX),VLSoilMicP(L,NY,NX))
      TortMicPM(M,L,NY,NX)=TortMicporew(THETWT)*(1.0_r8-SoilFracAsMacP(L,NY,NX))
    ELSE
!   standing water has tortuosity 0.7?
      TortMicPM(M,L,NY,NX)=0.7_r8
    ENDIF
    IF(VLMacP1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWH=VLWatMacPM(M,L,NY,NX)/VLMacP1(L,NY,NX)
      TortMacPM(M,L,NY,NX)=TortMacporew(THETWH)*SoilFracAsMacP(L,NY,NX)
    ELSE
      TortMacPM(M,L,NY,NX)=0.0_r8
    ENDIF
  ENDDO D9885
  end subroutine InitSoil3DModelMit


!------------------------------------------------------------------------------------------
  subroutine MicporeDarcyFlow(NY,NX,N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,FKSAT,THeatDueWatFlow,PSISV1,PSISVL)          
  implicit none
  integer, intent(in)  :: NY,NX,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: THETA1,THETAL,FKSAT
  real(r8), intent(out):: THeatDueWatFlow,PSISV1,PSISVL
  real(r8) :: AVE_CONDUCTANCE,CND1,CNDL,FLQ2,PSISTL,PSIST1,THETW1
  real(r8) :: FLQL,HWFLQL,FLQX,FLQZ,PSISTMP,THETWL
  integer :: K1,KL


  !     THETW1,THETWL=water concentrations in source,destination cells

  IF(PSISoilMatricPtmp(N3,N2,N1).GT.PSISoilAirEntry(N3,N2,N1).AND.PSISoilMatricPtmp(N6,N5,N4).GT.PSISoilAirEntry(N6,N5,N4))THEN
    THETW1=THETA1
    THETWL=THETAL
    K1=MAX(1,MIN(100,INT(100.0_r8*(POROS(N3,N2,N1)-THETW1)/POROS(N3,N2,N1))+1))
    KL=MAX(1,MIN(100,INT(100.0_r8*(POROS(N6,N5,N4)-THETWL)/POROS(N6,N5,N4))+1))
    PSISM1(N3,N2,N1)=PSISoilMatricPtmp(N3,N2,N1)
    PSISM1(N6,N5,N4)=PSISoilMatricPtmp(N6,N5,N4)
    !
    !     GREEN-AMPT FLOW IF ONE LAYER IS SATURATED
    !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POENTIAL)
    !
    !     GREEN-AMPT FLOW IF SOURCE CELL SATURATED
    !THETS=micropore soil water content
  ELSEIF(PSISoilMatricPtmp(N3,N2,N1).GT.PSISoilAirEntry(N3,N2,N1))THEN
    THETW1=THETA1
    THETWL=AMAX1(THETY(N6,N5,N4),AMIN1(POROS(N6,N5,N4),safe_adb(VLWatMicPX1(N6,N5,N4),VLSoilMicP(N6,N5,N4))))
    K1=MAX(1,MIN(100,INT(100.0_r8*(POROS(N3,N2,N1)-THETW1)/POROS(N3,N2,N1))+1))
    KL=MAX(1,MIN(100,INT(100.0_r8*(POROS(N6,N5,N4)-AMIN1(THETS(N6,N5,N4),THETWL))/POROS(N6,N5,N4))+1))
    PSISM1(N3,N2,N1)=PSISoilMatricPtmp(N3,N2,N1)
    
    IF(SoilMicPMassLayer(N6,N5,N4).GT.ZEROS(NY,NX))THEN
      IF(THETWL.LT.FieldCapacity(N6,N5,N4))THEN
        PSISM1(N6,N5,N4)=AMAX1(PSIHY,-EXP(LOGPSIFLD(N5,N4)+((LOGFldCapacity(N6,N5,N4)-LOG(THETWL)) &
          /FCD(N6,N5,N4)*LOGPSIMND(N5,N4))))
      ELSEIF(THETWL.LT.POROS(N6,N5,N4)-DTHETW)THEN
        PSISM1(N6,N5,N4)=-EXP(LOGPSIAtSat(N5,N4)+(((LOGPOROS(N6,N5,N4)-LOG(THETWL)) &
          /PSD(N6,N5,N4))**SRP(N6,N5,N4)*LOGPSIMXD(N5,N4)))
      ELSE
        THETWL=POROS(N6,N5,N4)
        PSISM1(N6,N5,N4)=PSISE(N6,N5,N4)
      ENDIF
    ELSE
      THETWL=POROS(N6,N5,N4)
      PSISM1(N6,N5,N4)=PSISE(N6,N5,N4)
    ENDIF
    !
    !     GREEN-AMPT FLOW IF ADJACENT CELL SATURATED
!
  ELSEIF(PSISoilMatricPtmp(N6,N5,N4).GT.PSISoilAirEntry(N6,N5,N4))THEN
    THETW1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1),safe_adb(VLWatMicPX1(N3,N2,N1),VLSoilMicP(N3,N2,N1))))
    THETWL=THETAL
    K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)-AMIN1(THETS(N3,N2,N1),THETW1))/POROS(N3,N2,N1))+1))
    KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)-THETWL)/POROS(N6,N5,N4))+1))
    IF(SoilMicPMassLayer(N3,N2,N1).GT.ZEROS(NY,NX))THEN
      IF(THETW1.LT.FieldCapacity(N3,N2,N1))THEN
        PSISTMP=-EXP(LOGPSIFLD(N2,N1)+(LOGFldCapacity(N3,N2,N1)-LOG(THETW1))/FCD(N3,N2,N1)*LOGPSIMND(N2,N1))
        PSISM1(N3,N2,N1)=AMAX1(PSIHY,PSISTMP)
      ELSEIF(THETW1.LT.POROS(N3,N2,N1)-DTHETW)THEN
        PSISM1(N3,N2,N1)=-EXP(LOGPSIAtSat(N2,N1)+(((LOGPOROS(N3,N2,N1)-LOG(THETW1)) &
          /PSD(N3,N2,N1))**SRP(N3,N2,N1)*LOGPSIMXD(N2,N1)))
      ELSE
        THETW1=POROS(N3,N2,N1)
        PSISM1(N3,N2,N1)=PSISE(N3,N2,N1)
      ENDIF
    ELSE
      THETW1=POROS(N3,N2,N1)
      PSISM1(N3,N2,N1)=PSISE(N3,N2,N1)
    ENDIF
    !
    !     RICHARDS FLOW IF NEITHER CELL IS SATURATED
    !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POTENTIAL)
    !
  ELSE
    THETW1=THETA1
    THETWL=THETAL
    K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)-THETW1)/POROS(N3,N2,N1))+1))
    KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)-THETWL)/POROS(N6,N5,N4))+1))
    PSISM1(N3,N2,N1)=PSISoilMatricPtmp(N3,N2,N1)
    PSISM1(N6,N5,N4)=PSISoilMatricPtmp(N6,N5,N4)
  ENDIF

  !
  !     HYDRAULIC CONUCTIVITY
  !
  !     CND1,CNDL=hydraulic conductivity of source,destination layer
  !     HydroCond3D=lateral(1,2),vertical(3) micropore hydraulic conductivity
  !
  IF(N3.EQ.NUM(NY,NX))THEN
    CND1=HydroCond3D(N,K1,N3,N2,N1)*FKSAT
  ELSE
    CND1=HydroCond3D(N,K1,N3,N2,N1)
  ENDIF
  CNDL=HydroCond3D(N,KL,N6,N5,N4)
  !
  !     TOTAL SOIL WATER POTENTIAL = MATRIC, GRAVIMETRIC + OSMOTIC
  !
  !     PSISM1,PSISH,PSISO=soil matric,gravitational,osmotic potentials
!
  PSIST1=PSISM1(N3,N2,N1)+PSISH(N3,N2,N1)+PSISoilOsmotic(N3,N2,N1)
  PSISTL=PSISM1(N6,N5,N4)+PSISH(N6,N5,N4)+PSISoilOsmotic(N6,N5,N4)
  PSISV1=PSISM1(N3,N2,N1)+PSISoilOsmotic(N3,N2,N1)
  PSISVL=PSISM1(N6,N5,N4)+PSISoilOsmotic(N6,N5,N4)
  !
  !     HYDRAULIC CONDUCTIVITY FROM CURRENT WATER CONTENT
  !     AND LOOKUP ARRAY GENERATED IN 'HOUR1'
  !
  !     CND1,CNDL=hydraulic conductivities in source,destination cells
  !     FKSAT=reduction in soil surface Ksat from rainfall energy impact
  !     AVE_CONDUCTANCE=source-destination hydraulic conductance
  !     DLYR=layer thickness
  !
  IF(CND1.GT.ZERO.AND.CNDL.GT.ZERO)THEN
    AVE_CONDUCTANCE=2.0_r8*CND1*CNDL/(CND1*DLYR(N,N6,N5,N4)+CNDL*DLYR(N,N3,N2,N1))
  ELSE
    AVE_CONDUCTANCE=0.0_r8
  ENDIF
  !
  !     WATER FLUX FROM WATER POTENTIALS, HYDRAULIC CONDUCTIVITY
  !     CONSTRAINED BY WATER POTENTIAL GRADIENT, COUPLED WITH
  !     CONVECTIVE HEAT FLUX FROM WATER FLUX
  !
  !     FLQX,FLQL=micropore water flux unlimited,limited by source water
  !     dts_HeatWatTP=time step of flux calculations
  !     VOLW2,VOLP1=water,air contents of source,destination micropores
  !     THeatDueWatFlow=convective heat flux from micropore water flux
  !     VLairMicP=excess water+ice relative to porosity
  !
  FLQX=AVE_CONDUCTANCE*(PSIST1-PSISTL)*AREA(N,N3,N2,N1)*dts_HeatWatTP
  IF(FLQX.GE.0.0_r8)THEN
    IF(THETW1.GT.THETS(N3,N2,N1))THEN
      FLQZ=FLQX+AMIN1((THETW1-THETS(N3,N2,N1))*VLSoilMicP(N3,N2,N1),AZMAX1((THETS(N6,N5,N4)-THETWL) &
        *VLSoilMicP(N6,N5,N4)))*dts_wat
    ELSE
      FLQZ=FLQX
    ENDIF
    FLQL=AZMAX1(AMIN1(FLQZ,VOLW2(N3,N2,N1)*dts_wat,VLairMicPc(N6,N5,N4)*dts_wat))
    FLQ2=AZMAX1(AMIN1(FLQX,VOLW2(N3,N2,N1)*dts_wat,VLairMicPc(N6,N5,N4)*dts_wat))
    !     FLQL1=(THETW1-THETS(N3,N2,N1))*VLSoilMicP(N3,N2,N1)
    !     FLQL2=(THETS(N6,N5,N4)-THETWL)*VLSoilMicP(N6,N5,N4)
    !     FLQL3=FLQX+AMIN1(FLQL1,AZMAX1(FLQL2))*dts_wat
    !     FLQL4=AZMAX1(AMIN1(FLQL3,VLairMicPc(N6,N5,N4)*dts_wat))
  ELSE
    IF(THETWL.GT.THETS(N6,N5,N4))THEN
      FLQZ=FLQX+AMAX1((THETS(N6,N5,N4)-THETWL)*VLSoilMicP(N6,N5,N4), &
        AZMIN1((THETW1-THETS(N3,N2,N1))*VLSoilMicP(N3,N2,N1)))*dts_wat
    ELSE
      FLQZ=FLQX
    ENDIF
    FLQL=AZMIN1(AMAX1(FLQZ,-VOLW2(N6,N5,N4)*dts_wat,-VLairMicPc(N3,N2,N1)*dts_wat))
    FLQ2=AZMIN1(AMAX1(FLQX,-VOLW2(N6,N5,N4)*dts_wat,-VLairMicPc(N3,N2,N1)*dts_wat))
    !     FLQL1=(THETS(N6,N5,N4)-THETWL)*VLSoilMicP(N6,N5,N4)
    !     FLQL2=(THETW1-THETS(N3,N2,N1))*VLSoilMicP(N3,N2,N1)
    !     FLQL3=FLQX+AMAX1(FLQL1,AZMIN1(FLQL2))*dts_wat
    !     FLQL4=AZMIN1(AMAX1(FLQL3,-VLairMicPc(N3,N2,N1)*dts_wat))
  ENDIF
  IF(N.EQ.ivertdir.AND.VLairMicP(N6,N5,N4).LT.0.0_r8)THEN
    FLQL=FLQL+AZMIN1(AMAX1(-VOLW2(N6,N5,N4)*dts_wat,VLairMicP(N6,N5,N4)))
    FLQ2=FLQ2+AZMIN1(AMAX1(-VOLW2(N6,N5,N4)*dts_wat,VLairMicP(N6,N5,N4)))
  ENDIF
  IF(FLQL.GT.0.0_r8)THEN
    HWFLQL=cpw*TKSoi1(N3,N2,N1)*FLQL
  ELSE
    HWFLQL=cpw*TKSoi1(N6,N5,N4)*FLQL
  ENDIF

  VOLW2(N3,N2,N1)=VOLW2(N3,N2,N1)-FLQL
  VOLW2(N6,N5,N4)=VOLW2(N6,N5,N4)+FLQL

  WatXChange2WatTable(N,N6,N5,N4)=FLQL
  WatXChange2WatTableX(N,N6,N5,N4)=FLQ2

  THeatDueWatFlow=HWFLQL
  end subroutine MicporeDarcyFlow
!------------------------------------------------------------------------------------------
  subroutine MacporeFLow(NY,NX,M,N,N1,N2,N3,N4,N5,N6,ConvectiveHeatFlxMacP,IFLGH)
  implicit none
  integer, intent(in) :: NY,NX,M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(out) :: ConvectiveHeatFlxMacP
  integer, intent(inout) :: IFLGH
  real(r8) :: FLWHX,PSISH1,PSISHL
  !     PSISH1,PSISHL=macropore total water potl in source,destination
  !     DLYR=layer thickness
  !     VLWatMacP1,VOLPH1=macropore water,air content

  IF(VLMacP1(N3,N2,N1).GT.ZEROS2(N2,N1).AND.VLMacP1(N6,N5,N4).GT.ZEROS2(N5,N4).AND.IFLGH.EQ.0)THEN
    PSISH1=PSISH(N3,N2,N1)+mGravAcceleration*DLYR(3,N3,N2,N1)*(AMIN1(1.0_r8,AZMAX1(VLWatMacP1(N3,N2,N1)/VLMacP1(N3,N2,N1)))-0.5_r8)
    PSISHL=PSISH(N6,N5,N4)+mGravAcceleration*DLYR(3,N6,N5,N4)*(AMIN1(1.0_r8,AZMAX1(VLWatMacP1(N6,N5,N4)/VLMacP1(N6,N5,N4)))-0.5_r8)
    !
    !     MACROPORE FLOW IF GRAVITATIONAL GRADIENT IS POSITIVE
    !     AND MACROPORE POROSITY EXISTS IN ADJACENT CELL
    !
    !     FLWHX,ConvectWaterFlowMacP=macropore water flux unltd,ltd by source water
    !     dts_HeatWatTP=time step of flux calculations
    !     VOLW2,VOLP1=water,air contents of source,destination micropores
    !     ConvectiveHeatFlxMacP=convective heat flux from micropore water flux
    !
    FLWHX=AVCNHL(N,N6,N5,N4)*(PSISH1-PSISHL)*AREA(N,N3,N2,N1)*dts_HeatWatTP
    IF(N.NE.ivertdir)THEN
      IF(PSISH1.GT.PSISHL)THEN
        ConvectWaterFlowMacP(N,N6,N5,N4)=AZMAX1(AMIN1(AMIN1(VLWatMacP1(N3,N2,N1),VLairMacPc(N6,N5,N4))*dts_wat,FLWHX))
      ELSEIF(PSISH1.LT.PSISHL)THEN
        ConvectWaterFlowMacP(N,N6,N5,N4)=AZMIN1(AMAX1(AMAX1(-VLWatMacP1(N6,N5,N4),-VLairMacPc(N3,N2,N1))*dts_wat,FLWHX))
      ELSE
        ConvectWaterFlowMacP(N,N6,N5,N4)=0.0_r8
      ENDIF
    ELSE
      ConvectWaterFlowMacP(N,N6,N5,N4)=AZMAX1(AMIN1(AMIN1(VLWatMacP1(N3,N2,N1)*dts_wat &
        +ConvectWaterFlowMacP(N,N3,N2,N1),VLairMacPc(N6,N5,N4)*dts_wat),FLWHX))
    ENDIF
    IF(N.EQ.ivertdir)THEN
      ConvectWaterFlowMacP(N,N6,N5,N4)=ConvectWaterFlowMacP(N,N6,N5,N4)+AZMIN1(VLairMacP(N6,N5,N4))
    ENDIF
    WaterFlowMacPi(M,N,N6,N5,N4)=ConvectWaterFlowMacP(N,N6,N5,N4)
  ELSE
    ConvectWaterFlowMacP(N,N6,N5,N4)=0.0_r8
    WaterFlowMacPi(M,N,N6,N5,N4)=0.0_r8
    IF(VLairMacPc(N6,N5,N4).LE.0.0_r8)IFLGH=1
  ENDIF
  IF(ConvectWaterFlowMacP(N,N6,N5,N4).GT.0.0_r8)THEN
    ConvectiveHeatFlxMacP=cpw*TKSoi1(N3,N2,N1)*ConvectWaterFlowMacP(N,N6,N5,N4)
  ELSE
    ConvectiveHeatFlxMacP=cpw*TKSoi1(N6,N5,N4)*ConvectWaterFlowMacP(N,N6,N5,N4)
  ENDIF
  end subroutine MacporeFLow  
!------------------------------------------------------------------------------------------

  subroutine WaterVaporFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectiveVaporFlux,ConvectiveHeatFluxMicP)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), intent(in) :: PSISV1,PSISVL
  real(r8), intent(out) :: ConvectiveVaporFlux,ConvectiveHeatFluxMicP
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
  !     WGSGL=vapor diffusivity
  !     ATCNVL=source,destination vapor conductance
  !     DLYR=soil layer depth
  !     PotentialVaporFlux,MaxVaporFlux=vapor flux unlimited,limited by vapor
  !     VPY=equilibrium vapor concentration
  !     dts_wat=time step for flux calculations
  !     ConvectiveVaporFlux,ConvectiveHeatFlux=vapor flux and its convective heat flux
!
  IF(THETPM(M,N3,N2,N1).GT.THETX.AND.THETPM(M,N6,N5,N4).GT.THETX)THEN
    TK11=TKSoi1(N3,N2,N1)
    TK12=TKSoi1(N6,N5,N4)

    VP1=vapsat(TK11)*EXP(18.0_r8*PSISV1/(RGAS*TK11))
    VPL=vapsat(TK12)*EXP(18.0_r8*PSISVL/(RGAS*TK12))
    CNV1=WGSGL(N3,N2,N1)*THETPM(M,N3,N2,N1)*POROQ*THETPM(M,N3,N2,N1)/POROS(N3,N2,N1)
    CNVL=WGSGL(N6,N5,N4)*THETPM(M,N6,N5,N4)*POROQ*THETPM(M,N6,N5,N4)/POROS(N6,N5,N4)
    ATCNVL=2.0_r8*CNV1*CNVL/(CNV1*DLYR(N,N6,N5,N4)+CNVL*DLYR(N,N3,N2,N1))
    !
    !     VAPOR FLUX FROM VAPOR PRESSURE AND DIFFUSIVITY,
    !     AND CONVECTIVE HEAT FLUX FROM VAPOR FLUX
!
    PotentialVaporFlux=ATCNVL*(VP1-VPL)*AREA(N,N3,N2,N1)*dts_HeatWatTP
    VPY=(VP1*VLsoiAirPM(M,N3,N2,N1)+VPL*VLsoiAirPM(M,N6,N5,N4))/(VLsoiAirPM(M,N3,N2,N1)+VLsoiAirPM(M,N6,N5,N4))
    MaxVaporFlux=(VP1-VPY)*VLsoiAirPM(M,N3,N2,N1)*dts_wat
    IF(PotentialVaporFlux.GE.0.0_r8)THEN
      ConvectiveVaporFlux=AZMAX1(AMIN1(PotentialVaporFlux,MaxVaporFlux))
      ConvectiveHeatFluxMicP=(cpw*TKSoi1(N3,N2,N1)+EvapLHTC)*ConvectiveVaporFlux
    ELSE
      ConvectiveVaporFlux=AZMIN1(AMAX1(PotentialVaporFlux,MaxVaporFlux))
      ConvectiveHeatFluxMicP=(cpw*TKSoi1(N6,N5,N4)+EvapLHTC)*ConvectiveVaporFlux
    ENDIF
    if((ConvectiveHeatFluxMicP>1.e10_r8 .or.curday>=282).and. max(N3,N6)<=3)then
      print*,'high ConvectiveHeatFluxMicP',ConvectiveHeatFluxMicP,VP1,VPL,VPY,ATCNVL
      print*,'TK',TK11,TKS(N3,N2,N1),TK12,TKS(N6,N5,N4)
      print*,N3,N2,N1,N6,N5,N4
    endif
  ELSE
    ConvectiveVaporFlux=0.0_r8
    ConvectiveHeatFluxMicP=0.0_r8
  ENDIF
  end subroutine WaterVaporFlow  
!------------------------------------------------------------------------------------------

  subroutine Solve4Heat(N,NY,NX,N1,N2,N3,N4,N5,N6,ATCNDL,ConvectiveHeatFluxMicP,HeatFlux2Ground)
  implicit none
  integer , intent(in) :: N,NY,NX,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: ATCNDL,ConvectiveHeatFluxMicP,HeatFlux2Ground
  real(r8) :: TK1X,TKLX,TKY,HFLWC,HFLWX,HFLWSX
  !
  !     HEAT FLOW FROM THERMAL CONDUCTIVITY AND TEMPERATURE GRADIENT
  !
  !     VLHeatCapacity,VHCPW=volumetric heat capacity of soil,snowpack
  !     TK1X,TKLX=interim temperatures of source,destination
  !     ConvectiveHeatFlux,HeatFlux2Ground=convective heat from soil vapor flux
  !     HeatFlux2Ground=storage heat flux from snowpack
  !     TKY=equilibrium source-destination temperature
  !     HFLWC,HFLWX=source-destination heat flux unltd,ltd by heat
  !     ATCNDL=source-destination thermal conductance
  !     HFLWSX=source-destination conductive heat flux
  !     HFLWL=total conductive+convective source-destination heat flux
  !
  IF(VLHeatCapacity(N3,N2,N1).GT.VHCPNX(NY,NX))THEN
    IF(N3.EQ.NUM(NY,NX).AND.VLHeatCapSnow(1,N2,N1).LE.VLHeatCapSnowMN(N2,N1))THEN
      !surface layer, not significant snowpack
      TK1X=TKSoi1(N3,N2,N1)-(ConvectiveHeatFluxMicP-HeatFlux2Ground)/VLHeatCapacity(N3,N2,N1)
      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoi1(N3,N2,N1)-ConvectiveHeatFluxMicP/VLHeatCapacity(N3,N2,N1)',&
          TKSoi1(N3,N2,N1),ConvectiveHeatFluxMicP,HeatFlux2Ground,VLHeatCapacity(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    ELSE
      TK1X=TKSoi1(N3,N2,N1)-ConvectiveHeatFluxMicP/VLHeatCapacity(N3,N2,N1)
      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoi1(N3,N2,N1)-ConvectiveHeatFluxMicP/VLHeatCapacity(N3,N2,N1)',&
          TKSoi1(N3,N2,N1),ConvectiveHeatFluxMicP,VLHeatCapacity(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    ENDIF
  ELSE
    TK1X=TKSoi1(N3,N2,N1)
  ENDIF

  IF(VLHeatCapacity(N6,N5,N4).GT.ZEROS(NY,NX))THEN
    TKLX=TKSoi1(N6,N5,N4)+ConvectiveHeatFluxMicP/VLHeatCapacity(N6,N5,N4)
  ELSE
    TKLX=TKSoi1(N6,N5,N4)
  ENDIF
  
  if(VLHeatCapacity(N3,N2,N1)+VLHeatCapacity(N6,N5,N4)>0._r8)then
    TKY=(VLHeatCapacity(N3,N2,N1)*TK1X+VLHeatCapacity(N6,N5,N4)*TKLX)/(VLHeatCapacity(N3,N2,N1)+VLHeatCapacity(N6,N5,N4))
  ELSE
    TKY=(TK1X+TKLX)/2._r8
  endif 
          !
  HFLWX=(TK1X-TKY)*VLHeatCapacity(N3,N2,N1)*dts_wat
  HFLWC=ATCNDL*(TK1X-TKLX)*AREA(N,N3,N2,N1)*dts_HeatWatTP
  IF(HFLWC.GE.0.0_r8)THEN
    HFLWSX=AZMAX1(AMIN1(HFLWX,HFLWC))
  ELSE
    HFLWSX=AZMIN1(AMAX1(HFLWX,HFLWC))
  ENDIF
  HeatFlowi(N,N6,N5,N4)=HeatFlowi(N,N6,N5,N4)+HFLWSX
  end subroutine Solve4Heat  
!------------------------------------------------------------------------------------------
  subroutine WaterTBLDrain(N,N1,N2,N3,M4,M5,M6,IFLGU,IFLGUH,RCHGFU,RCHGFT,DPTHH,XN)
  implicit none
  integer, intent(in) :: N,N1,N2,N3,M4,M5,M6,IFLGU,IFLGUH
  real(r8), intent(in):: RCHGFU,RCHGFT,DPTHH,XN
  real(r8) :: FLWT,FLWTHL,FLWTH,PSISWD,PSISWT,PSISWTH
!     MICROPORE DISCHARGE ABOVE WATER TABLE
!
!     IFLGU=micropore discharge flag to natural water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     DTBLG=water table slope
!     PSISWT=water potential driving micropore discharge
!     PSISoilMatricPtmp,PSISO=matric,osmotic water potential
!     DPTH,DTBLX=depth to layer midpoint,natural water table
!     DPTHT=depth to internal water table
!     FLWL=micropore discharge to natural water table
!     HFLWL=convective heat from discharge to natural water table
!     HydroCond3D=saturated hydraulic conductivity
!     AREAU=fraction of layer below natural water table
!
  IF(IFLGU.EQ.0.AND.(.not.isclose(RCHGFT,0._r8)))THEN
    PSISWD=XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-DTBLG(N2,N1))
    PSISWT=AZMIN1(-PSISoilMatricPtmp(N3,N2,N1)-0.03_r8*PSISoilOsmotic(N3,N2,N1) &
      +mGravAcceleration*(DPTH(N3,N2,N1)-DTBLX(N2,N1)) &
      -mGravAcceleration*AZMAX1(DPTH(N3,N2,N1)-DPTHT(N2,N1)))
    IF(PSISWT.LT.0.0_r8)PSISWT=PSISWT-PSISWD
    FLWT=PSISWT*HydroCond3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*(1.0_r8-AREAU(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*dts_HeatWatTP
    WatXChange2WatTable(N,M6,M5,M4)=XN*FLWT
    WatXChange2WatTableX(N,M6,M5,M4)=XN*FLWT
    HeatFlowi(N,M6,M5,M4)=cpw*TKSoi1(N3,N2,N1)*XN*FLWT
  ELSE
    WatXChange2WatTable(N,M6,M5,M4)=0.0_r8
    WatXChange2WatTableX(N,M6,M5,M4)=0.0_r8
    HeatFlowi(N,M6,M5,M4)=0.0_r8
  ENDIF
!
!     MACROPORE DISCHARGE ABOVE WATER TABLE
!
!     IFLGUH=macropore discharge flag to natural water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     DTBLG=water table slope
!     PSISWTH=water potential driving macropore discharge
!     PSISO=osmotic water potential
!     DPTHH,DTBLX=depth to layer macropore water,natural water table
!     DPTHT=depth to internal water table
!     FLWTH,FLWTHL=macropore discharge unltd,ltd by macropore water
!     HydroCondMacP1=macropore hydraulic conductivity
!     ConvectWaterFlowMacP=macropore discharge to natural water table
!     HFLWL=convective heat from discharge to natural water table
!     HydroCond3D=saturated hydraulic conductivity
  !     AREAU=fraction of layer below natural water table
!
  IF(IFLGUH.EQ.0.AND.(.not.isclose(RCHGFT,0._r8)).AND.VLMacP1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-DTBLG(N2,N1))
    PSISWTH=-0.03*PSISoilOsmotic(N3,N2,N1)+mGravAcceleration*(DPTHH-DTBLX(N2,N1)) &
      -mGravAcceleration*AZMAX1(DPTHH-DPTHT(N2,N1))
    IF(PSISWTH.LT.0.0_r8)PSISWTH=PSISWTH-PSISWD
    FLWTH=PSISWTH*HydroCondMacP1(N3,N2,N1)*AREA(N,N3,N2,N1) &
      *(1.0_r8-AREAU(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*dts_HeatWatTP
    FLWTHL=AMAX1(FLWTH,AZMIN1(-(VLWatMacP1(N3,N2,N1)*dts_wat &
      +ConvectWaterFlowMacP(3,N3,N2,N1)-ConvectWaterFlowMacP(3,N3+1,N2,N1))))
    ConvectWaterFlowMacP(N,M6,M5,M4)=XN*FLWTHL
    HeatFlowi(N,M6,M5,M4)=HeatFlowi(N,M6,M5,M4)+cpw*TKSoi1(N3,N2,N1)*XN*FLWTHL
  ELSE
    ConvectWaterFlowMacP(N,M6,M5,M4)=0.0_r8
  ENDIF
  end subroutine WaterTBLDrain  
!------------------------------------------------------------------------------------------
  subroutine TileDrain(N,N1,N2,N3,M4,M5,M6,IFLGD,IFLGDH,RCHGFT,RCHGFU,DPTHH,XN)
  implicit none
  integer, intent(in) :: N,N1,N2,N3,M4,M5,M6,IFLGD,IFLGDH
  real(r8),intent(in) :: RCHGFT,RCHGFU,DPTHH,XN
  real(r8) :: FLWT,FLWTH,FLWTHL,PSISWD,PSISWT,PSISWTH
!
  !     MICROPORE DISCHARGE ABOVE TILE DRAIN
  !
  !     IFLGD=micropore discharge flag to artificial water table
  !     PSISWD=water potential from water table slope
  !     XN,RCHG*=direction indicator,boundary flag
  !     SLOPE=sin(lateral slope)
  !     DLYR=layer width
  !     DTBLG=water table slope
  !     PSISWT=water potential driving micropore discharge
  !     PSISoilAirEntry,PSISO=matric,osmotic water potential
  !     DPTH,DTBLY=depth to layer midpoint,artificial water table
  !     DPTHT=depth to internal water table
  !     FLWL=micropore discharge to natural+artificial water table
  !     HFLWL=convective heat from dischg to natural+artifl water table
  !     HydroCond3D=saturated hydraulic conductivity
  !     AreaUnderWaterTBL=fraction of layer below artificial water table
!
  IF(IFLGD.EQ.0.AND.(.not.isclose(RCHGFT,0._r8)))THEN
    PSISWD=XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-DTBLG(N2,N1))
    PSISWT=AZMIN1(-PSISoilMatricPtmp(N3,N2,N1)-0.03_r8*PSISoilOsmotic(N3,N2,N1) &
      +mGravAcceleration*(DPTH(N3,N2,N1)-DTBLY(N2,N1))-mGravAcceleration*AZMAX1(DPTH(N3,N2,N1)-DPTHT(N2,N1)))
    IF(PSISWT.LT.0.0_r8)PSISWT=PSISWT-PSISWD
    FLWT=PSISWT*HydroCond3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*(1.0_r8-AreaUnderWaterTBL(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*dts_HeatWatTP
    WatXChange2WatTable(N,M6,M5,M4)=WatXChange2WatTable(N,M6,M5,M4)+XN*FLWT
    WatXChange2WatTableX(N,M6,M5,M4)=WatXChange2WatTableX(N,M6,M5,M4)+XN*FLWT
    HeatFlowi(N,M6,M5,M4)=HeatFlowi(N,M6,M5,M4)+cpw*TKSoi1(N3,N2,N1)*XN*FLWT
  ENDIF
!
!     MACROPORE DISCHARGE ABOVE TILE DRAIN
!
!     IFLGDH=macropore discharge flag to artificial water table
!     PSISWD=water potential from water table slope
!     XN,RCHG*=direction indicator,boundary flag
!     SLOPE=sin(lateral slope)
!     DLYR=layer width
!     DTBLG=water table slope
!     PSISWTH=water potential driving macropore discharge
!     PSISO=osmotic water potential
!     DPTHH,DTBLY=depth to layer macropore water,artificl water table
!     DPTHT=depth to internal water table
!     FLWTH,FLWTHL=macropore discharge unltd,ltd by macropore water
!     HydroCondMacP1=macropore hydraulic conductivity
!     ConvectWaterFlowMacP=macropore discharge to artificial water table
!     HFLWL=convective heat from discharge to artificial water table
!     HydroCond3D=saturated hydraulic conductivity
!     AreaUnderWaterTBL=fraction of layer below artificial water table
!
  IF(IFLGDH.EQ.0.AND.(.not.isclose(RCHGFT,0._r8)).AND.VLMacP1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    PSISWD=XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-DTBLG(N2,N1))
    PSISWTH=-0.03_r8*PSISoilOsmotic(N3,N2,N1)+mGravAcceleration*(DPTHH-DTBLY(N2,N1))-mGravAcceleration*AZMAX1(DPTHH-DPTHT(N2,N1))
    IF(PSISWTH.LT.0.0_r8)PSISWTH=PSISWTH-PSISWD
    FLWTH=PSISWTH*HydroCondMacP1(N3,N2,N1)*AREA(N,N3,N2,N1)*(1.0_r8-AreaUnderWaterTBL(N3,N2,N1))/(RCHGFU+1.0_r8)*RCHGFT*dts_HeatWatTP
    FLWTHL=AMAX1(FLWTH,AZMIN1(-(VLWatMacP1(N3,N2,N1)*dts_wat+ConvectWaterFlowMacP(3,N3,N2,N1)-ConvectWaterFlowMacP(3,N3+1,N2,N1))))
    ConvectWaterFlowMacP(N,M6,M5,M4)=ConvectWaterFlowMacP(N,M6,M5,M4)+XN*FLWTHL
    HeatFlowi(N,M6,M5,M4)=HeatFlowi(N,M6,M5,M4)+cpw*TKSoi1(N3,N2,N1)*XN*FLWTHL
  ENDIF
  end subroutine TileDrain
!------------------------------------------------------------------------------------------

  SUBROUTINE SubSufRecharge(NY,NX,N,N1,N2,N3,M4,M5,M6,DPTHH,RCHGFU,RCHGFT,XN,VOLP2,VOLPX2,VOLPH2)
  implicit none
  integer, intent(in) :: NY,NX,N,N1,N2,N3,M4,M5,M6
  real(r8),intent(in) :: RCHGFT,XN,DPTHH,RCHGFU
  real(r8),intent(inout):: VOLP2,VOLPX2,VOLPH2
  real(r8) :: FLWU,FLWUL,FLWUX,FLWUH,FLWUHL  
  real(r8) :: PSISUT,PSISUTH,PSISWD

  !     MICROPORE RECHARGE BELOW WATER TABLE
  !
  !     DPTHA=active layer depth
  !     VOLP2=air volume
  !     PSISWD=water potential from water table slope
  !     XN,RCHG*=direction indicator,boundary flag
  !     SLOPE=sin(lateral slope)
  !     DLYR=layer width
  !     DTBLG=water table slope
  !     PSISUT=water potential driving micropore recharge
  !     PSISoilAirEntry,PSISO=matric,osmotic water potential
  !     DPTH,DTBLX=depth to layer midpoint,natural water table
  !     FLWU,FLWUL=micropore recharge unltd,ltd by micropore air volume
  !     FLWL=micropore recharge from natural water table
  !     HFLWL=convective heat from recharge from natural water table
  !     HydroCond3D=saturated hydraulic conductivity
  !     AREAU=fraction of layer below natural water table
  !     VLairMicP=air filled 
      IF(DPTH(N3,N2,N1).GE.DTBLX(N2,N1)     &
        .AND.DPTHA(N2,N1).GT.DTBLX(N2,N1)   &
        .AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1) &
        .AND.(VOLP2.GT.ZEROS2(N2,N1).OR.SoiBulkDensity(N3,N2,N1).LE.ZERO) &
        .AND.VLairMicP(N3,N2,N1).GT.0.0_r8 &
        .AND.(.not.isclose(RCHGFT,0._r8)))THEN
        PSISWD=XN*0.005_r8*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-DTBLG(N2,N1))
        PSISUT=AZMAX1(-PSISoilMatricPtmp(N3,N2,N1)-0.03_r8*PSISoilOsmotic(N3,N2,N1)+mGravAcceleration*(DPTH(N3,N2,N1)-DTBLX(N2,N1)))
        IF(PSISUT.GT.0.0_r8)PSISUT=PSISUT+PSISWD
        FLWU=PSISUT*HydroCond3D(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)*AREAU(N3,N2,N1)/(RCHGFU+1.0)*RCHGFT*dts_HeatWatTP
        IF(SoiBulkDensity(N3,N2,N1).GT.ZERO)THEN
          FLWUL=AMIN1(FLWU,VOLP2)
          FLWUX=AMIN1(FLWU,VOLPX2)
        ELSE
          FLWUL=FLWU
          FLWUX=FLWU
        ENDIF
        WatXChange2WatTable(N,M6,M5,M4)=WatXChange2WatTable(N,M6,M5,M4)+XN*FLWUL
        WatXChange2WatTableX(N,M6,M5,M4)=WatXChange2WatTableX(N,M6,M5,M4)+XN*FLWUX
        HeatFlowi(N,M6,M5,M4)=HeatFlowi(N,M6,M5,M4)+cpw*TKSoi1(N3,N2,N1)*XN*FLWUL
        VOLP2=VOLP2-XN*WatXChange2WatTable(N,M6,M5,M4)
        VOLPX2=VOLPX2-XN*WatXChange2WatTableX(N,M6,M5,M4)
      ENDIF
!
      !     MACROPORE RECHARGE BELOW WATER TABLE
      !
      !     PSISWD=water potential from water table slope
      !     XN,RCHG*=direction indicator,boundary flag
      !     SLOPE=sin(lateral slope)
      !     DLYR=layer width
      !     DTBLG=water table slope
      !     PSISUTH=water potential driving macropore recharge
      !     PSISO=osmotic water potential
      !     DPTHH,DTBLX=depth to layer macropore water,natural water table
      !     DPTHT=depth to internal water table
      !     HydroCondMacP1=macropore hydraulic conductivity
      !     FLWUH,FLWUHL=macropore recharge unltd,ltd by macropore air volume
      !     ConvectWaterFlowMacP=macropore discharge to natural water table
      !     HFLWL=convective heat from discharge to natural water table
      !     HydroCond3D=saturated hydraulic conductivity
      !     AREAU=fraction of layer below natural water table
!
      IF(DPTHH.GT.DTBLX(N2,N1)                & !deeper than water table
        .AND.DPTHA(N2,N1).GT.DTBLX(N2,N1)     & !active layer below water table
        .AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1)   & !midlayer depth above active water layer
        .AND.VOLPH2.GT.ZEROS2(NY,NX)          & !macropore has air-filled fraction
        .AND.(.not.isclose(RCHGFT,0.0_r8)))THEN      !recharge is on
        PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)*(1.0_r8-DTBLG(N2,N1))
        PSISUTH=-0.03_r8*PSISoilOsmotic(N3,N2,N1)+mGravAcceleration*(DPTHH-DTBLX(N2,N1))
        IF(PSISUTH.GT.0.0_r8)PSISUTH=PSISUTH+PSISWD
        FLWUH=PSISUTH*HydroCondMacP1(N3,N2,N1)*AREA(N,N3,N2,N1)*AREAU(N3,N2,N1)/(RCHGFU+1.0_r8)*RCHGFT*dts_HeatWatTP
        FLWUHL=AMIN1(FLWUH,VOLPH2*dts_wat)
        ConvectWaterFlowMacP(N,M6,M5,M4)=ConvectWaterFlowMacP(N,M6,M5,M4)+XN*FLWUHL
        HeatFlowi(N,M6,M5,M4)=HeatFlowi(N,M6,M5,M4)+cpw*TKSoi1(N3,N2,N1)*XN*FLWUHL
        VOLPH2=VOLPH2-XN*ConvectWaterFlowMacP(N,M6,M5,M4)
      ENDIF
    
  end SUBROUTINE SubSufRecharge              
!------------------------------------------------------------------------------------------

  subroutine FreezeThawMit(NY,NX,L,N1,N2,N3)
  !
  !Description
  !applying freezethaw in cell (N3,N2,N1)
  implicit none
  integer, intent(in) :: NY,NX,L,N1,N2,N3
  real(r8) :: VLHeatCapacityBX,PSISMX,TFREEZ,VLWatMicP1X
  real(r8) :: VLWatMacP1X,VLHeatCapacityAX,TK1X,MacPIceHeatFlxFrezPt,MacPIceHeatFlxFrez
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
  PSISMX=PSISoilMatricPtmp(N3,N2,N1)+PSISoilOsmotic(N3,N2,N1)
  TFREEZ=-9.0959E+04_r8/(PSISMX-LtHeatIceMelt)
  VLWatMicP1X=VLWatMicP1(N3,N2,N1)+TWatCharge2MicP(N3,N2,N1)+FWatExMacP2MicP(N3,N2,N1)+FWatIrrigate2MicP1(N3,N2,N1)
  VLWatMacP1X=VLWatMacP1(N3,N2,N1)+TConvectWaterFlowMacP(N3,N2,N1)-FWatExMacP2MicP(N3,N2,N1)
  ENGY1=VLHeatCapacity(N3,N2,N1)*TKSoi1(N3,N2,N1)
  VLHeatCapacityX=VHeatCapacitySoilM(N3,N2,N1)+cpw*(VLWatMicP1X+VLWatMacP1X) &
    +cpi*(VLiceMicP1(N3,N2,N1)+VLiceMacP1(N3,N2,N1))

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    TK1X=(ENGY1+THeatFlowi(N3,N2,N1)+HeatIrrigation1(N3,N2,N1))/VLHeatCapacityX
    IF((TK1X.LT.TFREEZ.AND.VLWatMicP1(N3,N2,N1).GT.ZERO*VGeomLayer(N3,N2,N1)) &
      .OR.(TK1X.GT.TFREEZ.AND.VLiceMicP1(N3,N2,N1).GT.ZERO*VGeomLayer(N3,N2,N1)))THEN
      VLHeatCapacityAX=VHeatCapacitySoilM(N3,N2,N1)+cpw*VLWatMicP1X+cpi*VLiceMicP1(N3,N2,N1)
      MicPIceHeatFlxFrezPt=VLHeatCapacityAX*(TFREEZ-TK1X)/((1.0_r8+6.2913E-03_r8*TFREEZ)*(1.0_r8-0.10_r8*PSISMX))*dts_wat
      IF(MicPIceHeatFlxFrezPt.LT.0.0_r8)THEN
        !thaw
        MicPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1(N3,N2,N1)*dts_wat,MicPIceHeatFlxFrezPt)
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
  IF((TK1X.LT.TFice.AND.VLWatMacP1(N3,N2,N1).GT.ZERO*VGeomLayer(N3,N2,N1)) &
    .OR.(TK1X.GT.TFice.AND.VLiceMacP1(N3,N2,N1).GT.ZERO*VGeomLayer(N3,N2,N1)))THEN
    !there is freeze-thaw 
    VLHeatCapacityBX=cpw*VLWatMacP1X+cpi*VLiceMacP1(L,NY,NX)
    MacPIceHeatFlxFrezPt=VLHeatCapacityBX*(TFREEZ-TK1X)/((1.0_r8+6.2913E-03_r8*TFREEZ) &
        *(1.0_r8-0.10_r8*PSISMX))*dts_wat
    IF(MacPIceHeatFlxFrezPt.LT.0.0_r8)THEN
      !ice thaw
      MacPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMacP1(N3,N2,N1)*dts_wat,MacPIceHeatFlxFrezPt)
    ELSE
      !water freeze
      MacPIceHeatFlxFrez=AMIN1(LtHeatIceMelt*VLWatMacP1X*dts_wat,MacPIceHeatFlxFrezPt)
    ENDIF
    FIceThawMacP(N3,N2,N1)=-MacPIceHeatFlxFrez/LtHeatIceMelt
  ELSE
    MacPIceHeatFlxFrez=0.0_r8
    FIceThawMacP(N3,N2,N1)=0.0_r8
  ENDIF
  SoiPLIceHeatFlxFrez(N3,N2,N1)=MicPIceHeatFlxFrez+MacPIceHeatFlxFrez

  !
  !     TOTAL AND HOURLY ACCUMULATED FREEZE-THAW FLUXES
  !
  !     THAW,TLIceThawMacP=hourly accumulated freeze-thaw flux in micropores,macropores
  !     HTHAW=hourly accumulated freeze-thaw latent heat flux
  !     TWFLXL,TMLiceThawMacP=total accumulated freeze-thaw in micropores,macropores
  !     TLSoiPIceHeatFlxFrez1=total latent heat flux from melting

  ! Summarize fluxes for the current iteration
  TLIceThawMicP(N3,N2,N1)=TLIceThawMicP(N3,N2,N1)+FIceThawMicP(N3,N2,N1)
  TLIceThawMacP(N3,N2,N1)=TLIceThawMacP(N3,N2,N1)+FIceThawMacP(N3,N2,N1)
  TLSoiPIceHeatFlxFrez(N3,N2,N1)=TLSoiPIceHeatFlxFrez(N3,N2,N1)+SoiPLIceHeatFlxFrez(N3,N2,N1)
  TMLiceThawMicP(N3,N2,N1)=TMLiceThawMicP(N3,N2,N1)+FIceThawMicP(N3,N2,N1)
  TMLiceThawMacP(N3,N2,N1)=TMLiceThawMacP(N3,N2,N1)+FIceThawMacP(N3,N2,N1)
  TLSoiPIceHeatFlxFrez1(N3,N2,N1)=TLSoiPIceHeatFlxFrez1(N3,N2,N1)+SoiPLIceHeatFlxFrez(N3,N2,N1)
  end subroutine FreezeThawMit        

end module WatsubMod
