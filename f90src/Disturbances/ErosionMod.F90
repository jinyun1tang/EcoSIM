module ErosionMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use EcoSIMConfig, only : nlbiomcp => NumLiveMicrbCompts
  use EcoSIMConfig, only : ndbiomcp=> NumDeadMicrbCompts
  use EcoSIMConfig, only : jcplx1=> jcplxcm1, NumMicbFunGrupsPerCmplx => NumMicbFunGrupsPerCmplx,jcplx=>jcplxc
  use EcoSIMConfig, only : column_mode
  use EcoSiMParDataMod    , only : micpar  
  use EcoSIMCtrlMod, only : iErosionMode
  use MicrobialDataType
  use SOMDataType
  use EcoSIMSolverPar
  use FertilizerDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use LandSurfDataType
  use SurfSoilDataType
  use ChemTranspDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  USE SedimentDataType
  use GridDataType
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  real(r8), PARAMETER :: FSINK=0.01_r8

  real(r8), allocatable :: SedErosionM(:,:,:,:)
  real(r8), allocatable :: TERSED_col(:,:)
  real(r8), allocatable :: RDTSED_col(:,:)
  real(r8), allocatable :: FVOLIM(:,:)
  real(r8), allocatable :: FVLWatMicPM(:,:)
  real(r8), allocatable :: FracVol4Erosion(:,:)
  real(r8), allocatable :: BaseErosionRate(:,:)

  public :: erosion
  public :: InitErosion
  public :: DestructErosion
  contains

  subroutine InitErosion()

  implicit none

  allocate(SedErosionM(2,2,JV,JH))
  allocate(TERSED_col(JY,JX))
  allocate(RDTSED_col(JY,JX))
  allocate(FVOLIM(JY,JX))
  allocate(FVLWatMicPM(JY,JX))
  allocate(FracVol4Erosion(JY,JX))
  allocate(BaseErosionRate(JY,JX))

  end subroutine InitErosion
!------------------------------------------------------------------------------------------

  SUBROUTINE erosion(I,J,NHW,NHE,NVN,NVS)
      !
!     THIS SUBROUTINE CALCULATES DETACHMENT AND OVERLAND TRANSPORT
!     OF SURFACE SEDIMENT FROM PRECIPITATION IN WEATHER FILE AND
!     FROM RUNOFF IN 'WATSUB'
      implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: M

!     execution begins here
!
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    DO M=1,NPH
      call SedimentDetachmentM(M,NHW,NHE,NVN,NVS)

      call SedimentTransportM(M,NHW,NHE,NVN,NVS)
    ENDDO
  !
  !     INTERNAL SEDIMENT FLUXES
  !
    call LateralXGridSedmentFlux(NHW, NHE,NVN,NVS)
  !
  !     EXTERNAL BOUNDARY SEDIMENT FLUXES
  !
    if(.not.column_mode)call XBoundErosionFlux(NHW,NHE,NVN,NVS)
  ENDIF
  END subroutine erosion
!---------------------------------------------------------------------------------------------------

  subroutine SedimentDetachmentM(M,NHW,NHE,NVN,NVS)
!     INTERNAL TIME STEP AT WHICH SEDIMENT DETACHMENT AND TRANSPORT
!     IS CALCULATED. DETACHMENT IS THE SUM OF THAT BY RAINFALL AND
!     OVERLAND FLOW
!
  implicit none

  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8) :: DETW,SEDX,CSEDD,CSEDX
  real(r8) :: DEPI,DETI,SoilDetachRate,STPR
  integer :: NY,NX


  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      TERSED_col(NY,NX)=0._r8
      RDTSED_col(NY,NX)=0._r8
      FVOLIM(NY,NX)=AMIN1(1.0_r8,AZMAX1(XVLiceMicPM(M,NY,NX)/VWatStoreCapSurf_col(NY,NX)))
      FVLWatMicPM(NY,NX)=AMIN1(1.0_r8,AZMAX1(XVLMobileWatMicPM(M,NY,NX)/VWatStoreCapSurf_col(NY,NX)))
      FracVol4Erosion(NY,NX)=(1.0_r8-FVOLIM(NY,NX))*FVLWatMicPM(NY,NX)
!
!     DETACHMENT BY RAINFALL WHEN SURFACE WATER IS PRESENT
!
      IF(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.ZERO.AND.EnergyImpact4Erosion_colM(M,NY,NX).GT.0.0_r8 &
        .AND.XVLMobileWatMicPM(M,NY,NX).GT.ZEROS(NY,NX))THEN
!
!     DETACHMENT OF SEDIMENT FROM SURFACE SOIL DEPENDS ON RAINFALL
!     KINETIC ENERGY AND FROM DETACHMENT COEFFICIENT IN 'HOUR1'
!     ATTENUATED BY DEPTH OF SURFACE WATER
!
        DETW=SoilDetachability4Erosion1(NY,NX)*(1.0_r8+2.0_r8*VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)/VLMicP_vr(NU_col(NY,NX),NY,NX))
        SoilDetachRate=AMIN1(VLSoilMicPMass_vr(NU_col(NY,NX),NY,NX)*dts_wat &
          ,DETW*EnergyImpact4Erosion_colM(M,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
          *FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*FracSurfSnoFree_col(NY,NX)*(1.0-FVOLIM(NY,NX)))
        RDTSED_col(NY,NX)=RDTSED_col(NY,NX)+SoilDetachRate
      ENDIF
!
!     DEPOSITION OF SEDIMENT TO SOIL SURFACE FROM IMMOBILE SURFACE WATER
!

      IF(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.ZERO.AND.FracVol4Erosion(NY,NX).GT.ZERO)THEN

        SEDX=SED_col(NY,NX)+RDTSED_col(NY,NX)
        IF(XVLMobileWaterLitRM(M,NY,NX).LE.VWatStoreCapSurf_col(NY,NX))THEN
          IF(SEDX.GT.ZEROS(NY,NX))THEN
            CSEDD=AZMAX1(SEDX/XVLMobileWatMicPM(M,NY,NX))
        
            DEPI=AMAX1(-SEDX,VLS_col(NY,NX)*(0.0_r8-CSEDD)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
              *FracVol4Erosion(NY,NX)*FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*dts_HeatWatTP)
            RDTSED_col(NY,NX)=RDTSED_col(NY,NX)+DEPI
          ENDIF
        ELSE
    !     DETACHMENT IN SURFACE WATER FROM OVERLAND WATER
    !     VELOCITY FROM 'WATSUB' USED TO CALCULATE STREAM POWER,
    !     AND FROM SEDIMENT TRANSPORT CAPACITY VS. CURRENT SEDIMENT
    !     CONCENTRATION IN SURFACE WATER, MODIFIED BY SOIL COHESION
    !     FROM 'HOUR1'
    !
          STPR=1.0E+02_r8*RunoffVelocityM_col(M,NY,NX)*ABS(SLOPE_col(0,NY,NX))
          CSEDX=PtclDensitySurfLay_col(NY,NX)*CER_col(NY,NX)*AZMAX1(STPR-0.4_r8)**XER_col(NY,NX)
          CSEDD=AZMAX1(SEDX/XVLMobileWatMicPM(M,NY,NX))

          IF(CSEDX.GT.CSEDD)THEN
            DETI=AMIN1(VLSoilMicPMass_vr(NU_col(NY,NX),NY,NX)*dts_wat &
              ,SoilDetachability4Erosion2(NY,NX)*(CSEDX-CSEDD)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
              *FracVol4Erosion(NY,NX)*FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*dts_HeatWatTP)
          ELSE
            IF(SEDX.GT.ZEROS(NY,NX))THEN
              DETI=AMAX1(-SEDX,VLS_col(NY,NX)*(CSEDX-CSEDD)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
                *FracVol4Erosion(NY,NX)*FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*dts_HeatWatTP)
            ELSE
              DETI=0._r8
            ENDIF
          ENDIF
          RDTSED_col(NY,NX)=RDTSED_col(NY,NX)+DETI
        ENDIF
      ENDIF
!
!
!     TRANSPORT OF SEDIMENT IN OVERLAND FLOW FROM SEDIMENT
!     CONCENTRATION TIMES OVERLAND WATER FLUX FROM 'WATSUB'
!

    ENDDO
  ENDDO

  end subroutine SedimentDetachmentM
!------------------------------------------------------------------------------------------
  subroutine OverLandFlowSedTransp(M,NY,NX,NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: M,NY,NX,NHW,NHE,NVN,NVS
  integer :: N1,N2,N,NN,N4,N5,N4B,N5B
  real(r8) :: SEDX,CSEDE,FERM

  N1=NX
  N2=NY
!     SEDX=SED_col(N2,N1)
  IF(SurfRunoffPotentM_col(M,N2,N1).LE.0.0_r8.OR.SoilBulkDensity_vr(NU_col(N2,N1),N2,N1).LE.ZERO)THEN
    BaseErosionRate(N2,N1)=0._r8
  ELSE
    IF(XVLMobileWatMicPM(M,N2,N1).GT.ZEROS2(N2,N1))THEN
      SEDX=SED_col(N2,N1)+RDTSED_col(N2,N1)
      CSEDE=AZMAX1(SEDX/XVLMobileWatMicPM(M,N2,N1))
       BaseErosionRate(N2,N1)=AMIN1(SEDX,CSEDE*SurfRunoffPotentM_col(M,N2,N1)*(1.0_r8-FVOLIM(N2,N1)))
    ELSE
      BaseErosionRate(N2,N1)=0._r8
    ENDIF
  ENDIF
!
!     LOCATE INTERNAL BOUNDARIES
!
  IF(BaseErosionRate(N2,N1).GT.0.0_r8)THEN
    DO  N=1,2
      DO  NN=1,2
        IF(N.EQ.iWestEastDirection)THEN
          !skip boundariers        
          IF(NX.EQ.NHE .AND. NN.EQ.iFront .OR. NX.EQ.NHW .AND. NN.EQ.iBehind)THEN
            cycle
          !internal grids  
          ELSE            
            N4=NX+1;N5=NY    !dest in the east        
            N4B=NX-1;N5B=NY  !dest in the west
          ENDIF
        ELSEIF(N.EQ.iNorthSouthDirection)THEN
          !skip boundariers
          IF(NY.EQ.NVS .AND. NN.EQ.iFront .OR. NY.EQ.NVN .AND. NN.EQ.iBehind)THEN
            cycle
          ELSE            
            N4=NX;N5=NY+1   !dest in the south             
            N4B=NX;N5B=NY-1 !dest in the north
          ENDIF
        ENDIF

        IF(BaseErosionRate(N2,N1).GT.ZEROS(N2,N1))THEN
          IF(NN.EQ.iFront)THEN  !            
            FERM                      = QflxSurfRunoffM_2DH(M,N,iBehind,N5,N4)/SurfRunoffPotentM_col(M,N2,N1)
            SedErosionM(N,iBehind,N5,N4)    = BaseErosionRate(N2,N1)*FERM
            cumSed_Eros_2D(N,iBehind,N5,N4) = cumSed_Eros_2D(N,iBehind,N5,N4)+SedErosionM(N,iBehind,N5,N4)
          ELSE
            SedErosionM(N,iBehind,N5,N4)=0._r8
          ENDIF

          IF(NN.EQ.iBehind)THEN
            IF(N4B.GT.0.AND.N5B.GT.0)THEN
              !well-defined dest grid
              FERM                             = QflxSurfRunoffM_2DH(M,N,iFront,N5B,N4B)/SurfRunoffPotentM_col(M,N2,N1)
              SedErosionM(N,iFront,N5B,N4B)    = BaseErosionRate(N2,N1)*FERM
              cumSed_Eros_2D(N,iFront,N5B,N4B) = cumSed_Eros_2D(N,iFront,N5B,N4B)+SedErosionM(N,iFront,N5B,N4B)
            ELSE
              SedErosionM(N,iFront,N5B,N4B)=0._r8
            ENDIF
          ENDIF
        ELSE
          SedErosionM(N,iBehind,N5,N4)=0._r8
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            SedErosionM(N,iFront,N5B,N4B)=0._r8
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  end subroutine OverLandFlowSedTransp
!------------------------------------------------------------------------------------------
  subroutine SedimentTransportM(M,NHW,NHE,NVN,NVS)
!     INTERNAL TIME STEP AT WHICH SEDIMENT DETACHMENT AND TRANSPORT
!     IS CALCULATED. DETACHMENT IS THE SUM OF THAT BY RAINFALL AND
!     OVERLAND FLOW
!
  implicit none

  integer, intent(in) :: M,NHW,NHE,NVN,NVS

  real(r8) :: CSEDE,SEDX
  integer :: NGL
  integer :: N1,N2,NY,NX
!
!     BOUNDARY SEDIMENT FLUXES
!
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      N1=NX;N2=NY
      call OverLandFlowSedTransp(M,NY,NX,NHW,NHE,NVN,NVS)
      
      IF(SurfRunoffPotentM_col(M,N2,N1).LE.0.0.OR.SoilBulkDensity_vr(NU_col(N2,N1),N2,N1).LE.ZERO)THEN
        BaseErosionRate(N2,N1)=0._r8
      ELSE
        IF(XVLMobileWatMicPM(M,N2,N1).GT.ZEROS2(N2,N1))THEN
          SEDX  = SED_col(NY,NX)+RDTSED_col(NY,NX)
          CSEDE = AZMAX1(SEDX/XVLMobileWatMicPM(M,N2,N1))
          BaseErosionRate(N2,N1)=AMIN1(SEDX,CSEDE*SurfRunoffPotentM_col(M,N2,N1))
        ELSE
          BaseErosionRate(N2,N1)=0._r8
        ENDIF
      ENDIF
!
      if(.not. column_mode)call XBoundSedTranspM(M,NY,NX,NHW,NHE,NVN,NVS,N1,N2)
!
!     UPDATE STATE VARIABLES FOR SEDIMENT TRANSPORT
!
      SED_col(NY,NX)=SED_col(NY,NX)+TERSED_col(NY,NX)+RDTSED_col(NY,NX)

    ENDDO
  ENDDO
  end subroutine SedimentTransportM
!------------------------------------------------------------------------------------------
  subroutine XBoundSedTranspM(M,NY,NX,NHW,NHE,NVN,NVS,N1,N2)

  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX,NHW,NHE,NVN,NVS
!     LOCATE EXTERNAL BOUNDARIES
!
  integer :: N,NN,N4,N5,N4B,N5B
  integer :: M1,M2,M4,M5
  real(r8) :: RCHQF,FERM,XN

  D9580: DO  N=1,2
    D9575: DO  NN=1,2
      IF(N.EQ.iWestEastDirection)THEN
        N4=NX+1;N5=NY
        N4B=NX-1;N5B=NY
        IF(NN.EQ.iFront)THEN 
          IF(NX.EQ.NHE)THEN  !eastern boundary
            M1=NX;M2=NY
            M4=NX+1;M5=NY
            XN    = -1.0_r8
            RCHQF = RechargEastSurf_col(M2,M1)
          ELSE
            cycle
          ENDIF
        ELSEIF(NN.EQ.iBehind)THEN
          IF(NX.EQ.NHW)THEN  !western boundary
            M1=NX;M2=NY
            M4=NX;M5=NY
            XN    = 1.0_r8
            RCHQF = RechargWestSurf_col(M5,M4)
          ELSE
            cycle
          ENDIF
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        N4  = NX;N5  = NY+1
        N4B = NX;N5B = NY-1
        IF(NN.EQ.iFront)THEN          
          IF(NY.EQ.NVS)THEN !southern boundary
            M1 = NX;M2 = NY
            M4 = NX;M5 = NY+1
            XN = -1.0_r8
            RCHQF=RechargSouthSurf_col(M2,M1)
          ELSE
            cycle
          ENDIF
        ELSEIF(NN.EQ.iBehind)THEN          
          IF(NY.EQ.NVN)THEN !north boundary
            M1 = NX;M2 = NY
            M4 = NX;M5 = NY
            XN = 1.0_r8
            RCHQF=RechargNorthSurf_col(M5,M4)
          ELSE
            cycle
          ENDIF
        ENDIF
      ENDIF
!
!     SEDIMENT TRANSPORT ACROSS BOUNDARY FROM BOUNDARY RUNOFF
!     IN 'WATSUB' TIMES BOUNDARY SEDIMENT CONCENTRATION IN
!     SURFACE WATER
!
      IF(.not.XGridRunoffFlag_2DH(NN,N,N2,N1) .OR. isclose(RCHQF,0._r8) &
        .OR. BaseErosionRate(N2,N1).LE.ZEROS(N2,N1))THEN
        SedErosionM(N,NN,M5,M4)=0._r8
      ELSE
        IF(SurfRunoffPotentM_col(M,N2,N1).GT.ZEROS(N2,N1))THEN
          IF((NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
            .OR. (NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
            FERM                       = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/SurfRunoffPotentM_col(M,N2,N1)
            SedErosionM(N,NN,M5,M4)    = BaseErosionRate(N2,N1)*FERM
            cumSed_Eros_2D(N,NN,M5,M4) = cumSed_Eros_2D(N,NN,M5,M4)+SedErosionM(N,NN,M5,M4)
          ELSEIF((NN.EQ.iBehind.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
            .OR.(NN.EQ.iFront.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
            SedErosionM(N,NN,M5,M4)=0._r8
          ELSE
            SedErosionM(N,NN,M5,M4)=0._r8
          ENDIF
        ENDIF
      ENDIF
    ENDDO D9575
!
!     TOTAL SEDIMENT FLUXES
!
    D1202: DO  NN=1,2
      TERSED_col(N2,N1)=TERSED_col(N2,N1)+SedErosionM(N,NN,N2,N1)
      !legitimate front grid
      IF(IFLBM_2DH(M,N,NN,N5,N4).EQ.0)THEN
        TERSED_col(N2,N1)=TERSED_col(N2,N1)-SedErosionM(N,NN,N5,N4)
      ENDIF
      !legtimate behind grid
      IF(N4B.GT.0.AND.N5B.GT.0 .AND. NN.EQ.iFront)THEN
        TERSED_col(N2,N1)=TERSED_col(N2,N1)-SedErosionM(N,NN,N5B,N4B)
      ENDIF
    ENDDO D1202
  ENDDO D9580
  end subroutine XBoundSedTranspM
!------------------------------------------------------------------------------------------

  subroutine LateralXGridSedmentFlux(NHW, NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NGL,NTX,NTP,NE,idom
  integer :: K,NN,N,NO,M,NY,NX
  integer :: N1,N2,N4,N5,N4B,N5B,MID
  real(r8) :: FSEDER

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
        N1=NX
        N2=NY
        DO  N=1,2
          DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              IF(NX.EQ.NHE .AND. NN.EQ.iFront .OR. NX.EQ.NHW .AND. NN.EQ.iBehind)THEN
                cycle
              ELSE
                N4  = NX+1
                N5  = NY
                N4B = NX-1
                N5B = NY
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              IF(NY.EQ.NVS .AND. NN.EQ.iFront .OR. NY.EQ.NVN .AND. NN.EQ.iBehind)THEN
                cycle
              ELSE
                N4=NX
                N5=NY+1
                N4B=NX
                N5B=NY-1
              ENDIF
            ENDIF
!
!     FLUXES OF ALL SOLID MATERIALS IN SEDIMENT ARE CALCULATED
!     FROM VALUES OF THEIR CURRENT STATE VARIABLES MULTIPLIED
!     BY THE FRACTION OF THE TOTAL SURFACE LAYER MASS THAT IS
!     TRANSPORTED IN SEDIMENT
!
!     SOIL MINERALS
!
            IF(NN.EQ.iFront)THEN
              FSEDER = AMIN1(1.0,cumSed_Eros_2D(N,iBehind,N5,N4)/SoilMicPMassLayerMX(N2,N1))

              XSand_Eros_2D(N,iBehind,N5,N4) = FSEDER*SAND_vr(NU_col(N2,N1),N2,N1)
              XSilt_Eros_2D(N,iBehind,N5,N4) = FSEDER*SILT_vr(NU_col(N2,N1),N2,N1)
              XClay_Eros_2D(N,iBehind,N5,N4) = FSEDER*CLAY_vr(NU_col(N2,N1),N2,N1)
!
!     FERTILIZER POOLS
!
!     *ER=sediment flux from erosion.f
!
              XNH4Soil_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh4,NU_col(N2,N1),N2,N1)
              XNH3Soil_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh3,NU_col(N2,N1),N2,N1)
              XUreaSoil_Eros_molN_2D(N,iBehind,N5,N4) = FSEDER*FertN_mole_soil_vr(ifert_N_urea,NU_col(N2,N1),N2,N1)
              XNO3Soil_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_soil_vr(ifert_N_no3,NU_col(N2,N1),N2,N1)
              XPO4Soil_Eros_molP_2D(N,iBehind,N5,N4) =  FSEDER*FertP_mole_soil_vr(NU_col(N2,N1),N2,N1)

              XNH4Band_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(N2,N1),N2,N1)
              XNH3Band_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(N2,N1),N2,N1)
              XUreaBand_Eros_molN_2D(N,iBehind,N5,N4) = FSEDER*FertN_mole_Band_vr(ifert_N_urea_band,NU_col(N2,N1),N2,N1)
              XNO3Band_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_Band_vr(ifert_N_no3_band,NU_col(N2,N1),N2,N1)
              XPO4Band_Eros_molP_2D(N,iBehind,N5,N4)  = FSEDER*FertP_mole_band_vr(NU_col(N2,N1),N2,N1)
!
!     EXCHANGEABLE CATIONS AND ANIONS
!     sediment code
!
              DO NTX=idx_beg,idx_end
                trcx_Eros_mol_2D(NTX,N,iBehind,N5,N4)=FSEDER*trcx_solml_vr(NTX,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     PRECIPITATES
!
!     sediment code
!
              DO NTP=idsp_beg,idsp_end
                trcp_Eros_mol_2D(NTP,N,iBehind,N5,N4)   =FSEDER*trcp_saltpml_vr(NTP,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     ORGANIC MATTER
!
              DO  K=1,jcplx
                DO NO=1,NumMicbFunGrupsPerCmplx
                  DO NGL=JGniH(NO),JGnfH(NO)
                    DO M=1,nlbiomcp
                      MID=micpar%get_micb_id(M,NGL)
                      DO NE=1,NumPlantChemElms
                        OMEERhetr_2D(NE,MID,K,N,iBehind,N5,N4)=FSEDER*mBiomeHeter_vr(NE,MID,K,NU_col(N2,N1),N2,N1)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO NO=1,NumMicbFunGrupsPerCmplx
                DO NGL=JGniA(NO),JGnfA(NO)
                  DO M=1,nlbiomcp
                    MID=micpar%get_micb_id(M,NGL)
                    DO NE=1,NumPlantChemElms
                      OMEERauto_2D(NE,MID,N,iBehind,N5,N4)=FSEDER*mBiomeAutor_vr(NE,MID,NU_col(N2,N1),N2,N1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO  K=1,jcplx
                DO  M=1,ndbiomcp
                  DO NE=1,NumPlantChemElms
                    OMBioResdu_Eros_2D(NE,M,K,N,iBehind,N5,N4)=FSEDER*OMBioResdu_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                ENDDO
                DO idom=idom_beg,idom_end
                  SorbedOM_Eros_2D(idom,K,N,iBehind,N5,N4)=FSEDER*SorbedOM_vr(idom,K,NU_col(N2,N1),N2,N1)
                ENDDO
                DO  M=1,jsken
                  DO NE=1,NumPlantChemElms
                    SolidOM_Eros_2D(NE,M,K,N,iBehind,N5,N4)=FSEDER*SolidOM_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                  SolidOMAct_Eros_2D(M,K,N,iBehind,N5,N4)=FSEDER*SolidOMAct_vr(M,K,NU_col(N2,N1),N2,N1)
                ENDDO
              ENDDO
            ELSE
              XSand_Eros_2D(N,iBehind,N5,N4)=0._r8
              XSilt_Eros_2D(N,iBehind,N5,N4)=0._r8
              XClay_Eros_2D(N,iBehind,N5,N4)=0._r8
!
!     FERTILIZER POOLS
!
              XNH4Soil_Eros_molN_2D(N,iBehind,N5,N4)  = 0._r8
              XNH3Soil_Eros_molN_2D(N,iBehind,N5,N4)  = 0._r8
              XUreaSoil_Eros_molN_2D(N,iBehind,N5,N4) = 0._r8
              XNO3Soil_Eros_molN_2D(N,iBehind,N5,N4)  = 0._r8
              XPO4Soil_Eros_molP_2D(N,iBehind,N5,N4)  = 0._r8
              XNH4Band_Eros_molN_2D(N,iBehind,N5,N4)  = 0._r8
              XNH3Band_Eros_molN_2D(N,iBehind,N5,N4)  = 0._r8
              XUreaBand_Eros_molN_2D(N,iBehind,N5,N4) = 0._r8
              XNO3Band_Eros_molN_2D(N,iBehind,N5,N4)  = 0._r8
              XPO4Band_Eros_molP_2D(N,iBehind,N5,N4)  = 0._r8
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
              trcx_Eros_mol_2D(idx_beg:idx_end,N,iBehind,N5,N4)=0._r8
!
!     PRECIPITATES
!
              trcp_Eros_mol_2D(idsp_beg:idsp_end,N,iBehind,N5,N4)=0._r8
!
!     ORGANIC MATTER
!
              DO  K=1,jcplx
                DO  NO=1,NumMicbFunGrupsPerCmplx
                  DO NGL=JGniH(NO),JGnfH(NO)
                    DO  M=1,nlbiomcp
                      MID=micpar%get_micb_id(M,NGL)
                      DO NE=1,NumPlantChemElms
                        OMEERhetr_2D(NE,MID,K,N,iBehind,N5,N4)=0._r8
                      ENDDO
                    enddo
                  ENDDO
                enddo
              ENDDO

              DO  NO=1,NumMicbFunGrupsPerCmplx
                DO NGL=JGniA(NO),JGnfA(NO)
                  DO  M=1,nlbiomcp
                    MID=micpar%get_micb_id(M,NGL)
                    DO NE=1,NumPlantChemElms                    
                      OMEERauto_2D(NE,MID,N,iBehind,N5,N4)=0._r8
                    ENDDO
                  enddo
                ENDDO
              enddo

              DO  K=1,jcplx
                DO  M=1,ndbiomcp
                  DO NE=1,NumPlantChemElms
                    OMBioResdu_Eros_2D(NE,M,K,N,iBehind,N5,N4)=0._r8
                  ENDDO
                ENDDO
                DO idom=idom_beg,idom_end
                  SorbedOM_Eros_2D(idom,K,N,iBehind,N5,N4)=0._r8
                ENDDO
                DO  M=1,jsken
                  SolidOMAct_Eros_2D(M,K,N,iBehind,N5,N4)=0._r8       
                  DO NE=1,NumPlantChemElms         
                    SolidOM_Eros_2D(NE,M,K,N,iBehind,N5,N4)=0._r8
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
            IF(NN.EQ.iBehind)THEN
              IF(N4B.GT.0.AND.N5B.GT.0)THEN
                FSEDER=AMIN1(1.0_r8,cumSed_Eros_2D(N,iFront,N5B,N4B)/SoilMicPMassLayerMX(N2,N1))
                XSand_Eros_2D(N,iFront,N5B,N4B) = FSEDER*SAND_vr(NU_col(N2,N1),N2,N1)
                XSilt_Eros_2D(N,iFront,N5B,N4B) = FSEDER*SILT_vr(NU_col(N2,N1),N2,N1)
                XClay_Eros_2D(N,iFront,N5B,N4B) = FSEDER*CLAY_vr(NU_col(N2,N1),N2,N1)
!
!     FERTILIZER POOLS
!
                XNH4Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh4,NU_col(N2,N1),N2,N1)
                XNH3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh3,NU_col(N2,N1),N2,N1)
                XUreaSoil_Eros_molN_2D(N,iFront,N5B,N4B) = FSEDER*FertN_mole_soil_vr(ifert_N_urea,NU_col(N2,N1),N2,N1)
                XNO3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_soil_vr(ifert_N_no3,NU_col(N2,N1),N2,N1)
                XPO4Soil_Eros_molP_2D(N,iFront,N5B,N4B)  = FSEDER*FertP_mole_soil_vr(NU_col(N2,N1),N2,N1)
                XNH4Band_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(N2,N1),N2,N1)
                XNH3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(N2,N1),N2,N1)
                XUreaBand_Eros_molN_2D(N,iFront,N5B,N4B) = FSEDER*FertN_mole_Band_vr(ifert_N_urea_band,NU_col(N2,N1),N2,N1)
                XNO3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_Band_vr(ifert_N_no3_band,NU_col(N2,N1),N2,N1)
                XPO4Band_Eros_molP_2D(N,iFront,N5B,N4B)  = FSEDER*FertP_mole_band_vr(NU_col(N2,N1),N2,N1)

!
!     EXCHANGEABLE CATIONS AND ANIONS
!
!     sediment code
!
                DO NTX=idx_beg,idx_end
                  trcx_Eros_mol_2D(NTX,N,iFront,N5B,N4B)=FSEDER*trcx_solml_vr(NTX,NU_col(N2,N1),N2,N1)
                ENDDO
!
!     PRECIPITATES
!
!     sediment code
                DO NTP=idsp_beg,idsp_end
                  trcp_Eros_mol_2D(NTP,N,iFront,N5B,N4B)=FSEDER*trcp_saltpml_vr(NTP,NU_col(N2,N1),N2,N1)
                ENDDO
!
!     ORGANIC MATTER
!
                DO  K=1,jcplx
                  DO  NO=1,NumMicbFunGrupsPerCmplx
                    DO NGL=JGniH(NO),JGnfH(NO)
                      DO  M=1,nlbiomcp
                        MID=micpar%get_micb_id(M,NGL)      
                        DO NE=1,NumPlantChemElms                                                 
                          OMEERhetr_2D(NE,MID,K,N,iFront,N5B,N4B)=FSEDER*mBiomeHeter_vr(NE,MID,K,NU_col(N2,N1),N2,N1)
                        ENDDO
                      enddo
                    enddo
                  ENDDO
                ENDDO
                DO  NO=1,NumMicbFunGrupsPerCmplx
                  DO NGL=JGniA(NO),JGnfA(NO)
                    DO  M=1,nlbiomcp
                      MID=micpar%get_micb_id(M,NGL)
                      DO NE=1,NumPlantChemElms         
                        OMEERauto_2D(NE,MID,N,iFront,N5B,N4B)=FSEDER*mBiomeAutor_vr(NE,MID,NU_col(N2,N1),N2,N1)
                      ENDDO
                    enddo
                  enddo
                ENDDO

                DO  K=1,jcplx
                  DO  M=1,ndbiomcp
                    DO NE=1,NumPlantChemElms   
                      OMBioResdu_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=FSEDER*OMBioResdu_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                    ENDDO
                  ENDDO
                  DO idom=idom_beg,idom_end
                    SorbedOM_Eros_2D(idom,K,N,iFront,N5B,N4B)=FSEDER*SorbedOM_vr(idom,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                  DO  M=1,jsken
                    SolidOMAct_Eros_2D(M,K,N,iFront,N5B,N4B)=FSEDER*SolidOMAct_vr(M,K,NU_col(N2,N1),N2,N1)
                    DO NE=1,NumPlantChemElms                      
                      SolidOM_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=FSEDER*SolidOM_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                    ENDDO  
                  ENDDO
                ENDDO
              ELSE
                XSand_Eros_2D(N,iFront,N5B,N4B) = 0._r8
                XSilt_Eros_2D(N,iFront,N5B,N4B) = 0._r8
                XClay_Eros_2D(N,iFront,N5B,N4B) = 0._r8
!
!     FERTILIZER POOLS
!
                XNH4Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XNH3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XUreaSoil_Eros_molN_2D(N,iFront,N5B,N4B) = 0._r8
                XNO3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XPO4Soil_Eros_molP_2D(N,iFront,N5B,N4B)  = 0._r8
                XNH4Band_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XNH3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XUreaBand_Eros_molN_2D(N,iFront,N5B,N4B) = 0._r8
                XNO3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XPO4Band_Eros_molP_2D(N,iFront,N5B,N4B)= 0._r8
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
                trcx_Eros_mol_2D(idx_beg:idx_end,N,iFront,N5B,N4B)=0._r8
!
!     PRECIPITATES
!
                trcp_Eros_mol_2D(idsp_beg:idsp_end,N,iFront,N5B,N4B)=0._r8
!
!     ORGANIC MATTER
!
                DO  K=1,jcplx
                  DO  NO=1,NumMicbFunGrupsPerCmplx
                    DO NGL=JGniH(NO),JGnfH(NO)
                      DO  M=1,nlbiomcp
                        MID=micpar%get_micb_id(M,NGL)
                        DO NE=1,NumPlantChemElms  
                          OMEERhetr_2D(NE,MID,K,N,iFront,N5B,N4B)=0._r8
                        ENDDO  
                      enddo
                    ENDDO
                  enddo
                ENDDO

                DO  NO=1,NumMicbFunGrupsPerCmplx
                  DO NGL=JGniA(NO),JGnfA(NO)
                    DO  M=1,nlbiomcp
                      MID=micpar%get_micb_id(M,NGL)
                      DO NE=1,NumPlantChemElms                        
                        OMEERauto_2D(NE,MID,N,iFront,N5B,N4B)=0._r8
                      ENDDO  
                    enddo
                  ENDDO
                enddo

                DO  K=1,jcplx
                  DO  M=1,ndbiomcp
                    DO NE=1,NumPlantChemElms                    
                      OMBioResdu_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=0._r8
                    ENDDO  
                  ENDDO
                  DO idom=idom_beg,idom_end
                    SorbedOM_Eros_2D(idom,K,N,iFront,N5B,N4B)=0._r8
                  ENDDO
                  DO  M=1,jsken
                    SolidOMAct_Eros_2D(M,K,N,iFront,N5B,N4B)=0._r8
                    DO NE=1,NumPlantChemElms       
                      SolidOM_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=0._r8
                    ENDDO  
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
    ENDDO
  ENDDO
  end subroutine LateralXGridSedmentFlux
!----------------------------------------------------------------------------------------------

  subroutine XBoundErosionFlux(NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) :: NHW,NHE,NVN,NVS

  real(r8) :: RCHQF,FSEDER
  integer :: NGL
  integer :: NY,NX,NE,idom
  integer :: N,NN,K,NO,M,NTX,NTP
  integer :: N1,N2,N4,N5
  integer :: M1,M2,M4,M5,MID
  real(r8) :: XN

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      IF(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.ZERO)THEN
        N1=NX
        N2=NY
        D8980: DO  N=1,2
          D8975: DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              N4=NX+1
              N5=NY
              IF(NN.EQ.iFront)THEN
                IF(NX.EQ.NHE)THEN
                  M1    = NX
                  M2    = NY
                  M4    = NX+1
                  M5    = NY
                  XN    = -1.0_r8
                  RCHQF = RechargEastSurf_col(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                IF(NX.EQ.NHW)THEN
                  M1 = NX
                  M2 = NY
                  M4 = NX
                  M5 = NY
                  XN = 1.0_r8
                  RCHQF=RechargWestSurf_col(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              N4=NX
              N5=NY+1
              IF(NN.EQ.iFront)THEN
                IF(NY.EQ.NVS)THEN
                  M1 = NX
                  M2 = NY
                  M4 = NX
                  M5 = NY+1
                  XN = -1.0_r8
                  RCHQF=RechargSouthSurf_col(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY
                  M4=NX
                  M5=NY
                  XN=1.0_r8
                  RCHQF=RechargNorthSurf_col(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ENDIF
            IF(.not.XGridRunoffFlag_2DH(NN,N,N2,N1).OR.isclose(RCHQF,0._r8) &
              .OR. ABS(cumSed_Eros_2D(N,NN,M5,M4)).LE.ZEROS(N2,N1))THEN
              XSand_Eros_2D(N,NN,M5,M4)=0._r8
              XSilt_Eros_2D(N,NN,M5,M4)=0._r8
              XClay_Eros_2D(N,NN,M5,M4)=0._r8
!
!     FERTILIZER POOLS
!
              XNH4Soil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNH3Soil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XUreaSoil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNO3Soil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XPO4Soil_Eros_molP_2D(N,NN,M5,M4)=0._r8
              XNH4Band_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNH3Band_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XUreaBand_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNO3Band_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XPO4Band_Eros_molP_2D(N,NN,M5,M4)=0._r8
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
              trcx_Eros_mol_2D(idx_beg:idx_end,N,NN,M5,M4)=0._r8
!
!     PRECIPITATES
!
              trcp_Eros_mol_2D(idsp_beg:idsp_end,N,NN,M5,M4)=0._r8
!
!     ORGANIC MATTER
!
              DO  K=1,jcplx
                DO  NO=1,NumMicbFunGrupsPerCmplx
                  DO NGL=JGniH(NO),JGnfH(NO)
                    DO  M=1,nlbiomcp
                      MID=micpar%get_micb_id(M,NGL)
                      DO NE=1,NumPlantChemElms
                        OMEERhetr_2D(NE,MID,K,N,NN,M5,M4)=0._r8
                      ENDDO
                    enddo
                  ENDDO
                enddo
              enddo
              DO  NO=1,NumMicbFunGrupsPerCmplx
                DO NGL=JGniA(NO),JGnfA(NO)
                  DO  M=1,nlbiomcp
                    MID=micpar%get_micb_id(M,NGL)
                    DO NE=1,NumPlantChemElms       
                      OMEERauto_2D(NE,MID,N,NN,M5,M4)=0._r8
                    ENDDO  
                  enddo
                ENDDO
              enddo
              DO  K=1,jcplx
                DO  M=1,ndbiomcp
                  DO NE=1,NumPlantChemElms                       
                    OMBioResdu_Eros_2D(NE,M,K,N,NN,M5,M4)=0._r8
                  ENDDO  
                ENDDO
                DO idom=idom_beg,idom_end
                  SorbedOM_Eros_2D(idom,K,N,NN,M5,M4)=0._r8
                ENDDO
                DO  M=1,jsken
                  SolidOMAct_Eros_2D(M,K,N,NN,M5,M4)=0._r8
                  DO NE=1,NumPlantChemElms                       
                    SolidOM_Eros_2D(NE,M,K,N,NN,M5,M4)=0._r8
                  ENDDO
                ENDDO
              ENDDO
!
!     CALCULATE FRACTION OF SURACE MATERIAL ERODED
!
            ELSE
              FSEDER=AMIN1(1.0_r8,cumSed_Eros_2D(N,NN,N5,N4)/SoilMicPMassLayerMX(N2,N1))
!
!     SOIL MINERALS
!
              XSand_Eros_2D(N,NN,M5,M4)=FSEDER*SAND_vr(NU_col(N2,N1),N2,N1)
              XSilt_Eros_2D(N,NN,M5,M4)=FSEDER*SILT_vr(NU_col(N2,N1),N2,N1)
              XClay_Eros_2D(N,NN,M5,M4)=FSEDER*CLAY_vr(NU_col(N2,N1),N2,N1)
!
!     FERTILIZER POOLS
!
              XNH4Soil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_nh4,NU_col(N2,N1),N2,N1)
              XNH3Soil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_nh3,NU_col(N2,N1),N2,N1)
              XUreaSoil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_urea,NU_col(N2,N1),N2,N1)
              XNO3Soil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_no3,NU_col(N2,N1),N2,N1)
              XPO4Soil_Eros_molP_2D(N,NN,M5,M4)=FSEDER*FertP_mole_soil_vr(NU_col(N2,N1),N2,N1)

              XNH4Band_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(N2,N1),N2,N1)
              XNH3Band_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(N2,N1),N2,N1)
              XUreaBand_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_urea_band,NU_col(N2,N1),N2,N1)
              XNO3Band_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_no3_band,NU_col(N2,N1),N2,N1)
              XPO4Band_Eros_molP_2D(N,NN,M5,M4)=FSEDER*FertP_mole_band_vr(NU_col(N2,N1),N2,N1)
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
              DO NTX=idx_beg,idx_end
                trcx_Eros_mol_2D(NTX,N,NN,M5,M4)=FSEDER*trcx_solml_vr(NTX,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     PRECIPITATES
!
              DO NTP=idsp_beg,idsp_end
                trcp_Eros_mol_2D(NTP,N,NN,M5,M4)=FSEDER*trcp_saltpml_vr(NTP,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     ORGANIC MATTER
!
              DO  K=1,jcplx
                DO NO=1,NumMicbFunGrupsPerCmplx
                  DO NGL=JGniH(NO),JGnfH(NO)
                    DO M=1,nlbiomcp
                      MID=micpar%get_micb_id(M,NGL)                    
                      DO NE=1,NumPlantChemElms                       
                        OMEERhetr_2D(NE,MID,K,N,NN,M5,M4)=FSEDER*mBiomeHeter_vr(NE,MID,K,NU_col(N2,N1),N2,N1)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
              DO NO=1,NumMicbFunGrupsPerCmplx
                DO NGL=JGniA(NO),JGnfH(NO)
                  DO M=1,nlbiomcp
                    MID=micpar%get_micb_id(M,NGL)
                    DO NE=1,NumPlantChemElms                       
                      OMEERauto_2D(NE,MID,N,NN,M5,M4)=FSEDER*mBiomeAutor_vr(NE,MID,NU_col(N2,N1),N2,N1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO  K=1,jcplx
                DO  M=1,ndbiomcp
                  DO NE=1,NumPlantChemElms                       
                    OMBioResdu_Eros_2D(NE,M,K,N,NN,M5,M4)=FSEDER*OMBioResdu_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO  
                ENDDO
                DO idom=idom_beg,idom_end
                  SorbedOM_Eros_2D(idom,K,N,NN,M5,M4)=FSEDER*SorbedOM_vr(idom,K,NU_col(N2,N1),N2,N1)
                ENDDO
                DO  M=1,jsken
                  SolidOMAct_Eros_2D(M,K,N,NN,M5,M4)=FSEDER*SolidOMAct_vr(M,K,NU_col(N2,N1),N2,N1)
                  DO NE=1,NumPlantChemElms                       
                    SolidOM_Eros_2D(NE,M,K,N,NN,M5,M4)=FSEDER*SolidOM_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO D8975
        ENDDO D8980
      ENDIF
    ENDDO
  ENDDO

  end subroutine XBoundErosionFlux

!------------------------------------------------------------------------------------------
  subroutine DestructErosion()
  use abortutils, only : destroy

  implicit none

  call destroy(SedErosionM)
  call destroy(TERSED_col)
  call destroy(RDTSED_col)
  call destroy(FVOLIM)
  call destroy(FVLWatMicPM)
  call destroy(FracVol4Erosion)
  call destroy(BaseErosionRate)

  end subroutine DestructErosion
  end module ErosionMod
