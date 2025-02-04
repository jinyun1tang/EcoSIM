module SoilHydroParaMod

  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar
  use minimathmod,      only: isclose, AZMAX1, AZMIN1
  use DebugToolMod
  use EcoSimConst
  use SoilWaterDataType
  use SoilPropertyDataType
  use SoilHeatDataType
  use GridDataType
  use SoilPhysDataType
  use AqueChemDatatype
  use EcoSIMCtrlDataType
  use SOMDataType
  use FlagDataType
  use LandSurfDataType
  use EcoSIMHistMod
  use EcoSIMConfig
  use SurfLitterDataType
  use EcoSIMCtrlMod
  use PhysPars

implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8), parameter :: FORGW=0.25E+06_r8 !threshold for  C concentration in organic soil 	g Mg-1
  real(r8), parameter :: mGravAccelerat=1.e-3_r8*GravAcceleration  !gravitational constant devided by 1000.
  real(r8), parameter :: tiny_val=1.e-6_r8   !minimum value 
  public :: GetSoilHydraulicVars
  public :: SoilHydroProperty
  public :: LitterHydroproperty
  PUBLIC :: ComputeSoilHydroPars
contains
!------------------------------------------------------------------------------------------

  subroutine GetSoilHydraulicVars(I,J,NY,NX)
  !
  !DESCRIPTIONS
  !compute hydraulic properties
  !called in hour1.F90 before doing hydrology
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  character(len=*), parameter :: subname='GetSoilHydraulicVars'
  REAL(R8) :: FCX,FCLX
  real(r8) :: FCDX
  real(r8) :: PSDX
  real(r8) :: THETW1
  real(r8) :: WPX,WPLX
  integer :: K,L

  ! begin_execution
  call PrintInfo('beg '//subname)

  DO L=NUI(NY,NX),NLI(NY,NX)
  ! WATER POTENTIALS
  !
  ! FC,WP=water contents at field capacity,wilting point,saturation
  ! PSISM,PSISE=matric,saturation water potential
  ! SRP=parameter for deviation from linear log-log water retention
  ! FC,WP=water contents at field capacity,wilting point from soil file
  ! FCL,WPL=log FC,WP
  ! FCD,PSD=FCL-WPL,log(POROS)-FCL
  ! FCI,WPI=FC,WP of ice
  ! THETIX=ice concentration

    !significant soil mass
    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX) .AND. VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      THETW1=AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLWatMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)),tiny_val)
      !water content less than field capacity
      IF(THETW1.LT.FieldCapacity_vr(L,NY,NX))THEN
        PSISoilMatricP_vr(L,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX) &
          +((LOGFldCapacity_vr(L,NY,NX)-LOG(THETW1))/FCD_vr(L,NY,NX)*LOGPSIMND(NY,NX))))
      !water content exceeds field capacity but less than sturation    
      ELSE IF(THETW1.LT.POROS_vr(L,NY,NX)-DTHETW)THEN
        PSISoilMatricP_vr(L,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(L,NY,NX)-LOG(THETW1)) &
          /PSD_vr(L,NY,NX))**SRP_vr(L,NY,NX)*LOGPSIMXD(NY,NX)))
      !saturated    
      ELSE
        PSISoilMatricP_vr(L,NY,NX)=PSISE_vr(L,NY,NX)
      ENDIF
      !ice 
    ELSE IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX) .and. THETI_vr(L,NY,NX)>ZEROS2(NY,NX))THEN
      FCX  = FCI*THETI_vr(L,NY,NX)
      WPX  = WPI*THETI_vr(L,NY,NX)
      FCLX = LOG(FCX)
      WPLX = LOG(WPX)
      PSDX = LOGPOROS_vr(L,NY,NX)-FCLX
      FCDX = FCLX-WPLX
      IF(THETW_vr(L,NY,NX).LT.FCX)THEN
        PSISoilMatricP_vr(L,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX) &
          +((FCLX-LOG(THETW_vr(L,NY,NX)))/FCDX*LOGPSIMND(NY,NX))))
      ELSE IF(THETW_vr(L,NY,NX).LT.POROS_vr(L,NY,NX)-DTHETW)THEN
        PSISoilMatricP_vr(L,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(L,NY,NX)-LOG(THETW_vr(L,NY,NX))) &
          /PSDX)*LOGPSIMXD(NY,NX)))
      ELSE
        PSISoilMatricP_vr(L,NY,NX)=PSISE_vr(L,NY,NX)
      ENDIF
    ELSE
      PSISoilMatricP_vr(L,NY,NX)=PSISE_vr(L,NY,NX)
    ENDIF
!
!     SOIL OSMOTIC, GRAVIMETRIC AND MATRIC WATER POTENTIALS
!
!     PSISM,PSISO,PSIGrav_vr,PSIST=matric,osmotic,gravimetric,total water potential
!
    PSISoilOsmotic_vr(L,NY,NX)          = -RGASC*1.E-6_r8*TKS_vr(L,NY,NX)*SolutesIonConc_vr(L,NY,NX)
    PSIGrav_vr(L,NY,NX)                 = mGravAccelerat*(ALT(NY,NX)-SoilDepthMidLay_vr(L,NY,NX))
    ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX) = AZMIN1(PSISoilMatricP_vr(L,NY,NX)+PSISoilOsmotic_vr(L,NY,NX)+PSIGrav_vr(L,NY,NX))
!    write(212,*)I+J/24.,L,ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX),PSISoilMatricP_vr(L,NY,NX)
!
!     SOIL RESISTANCE TO ROOT PENETRATION
!
!     SoilResit4RootPentrate_vr=soil resistance to root penetration (MPa)
!
!     IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
!     CCLAY_vrT=CCLAY_vr(L,NY,NX)*1.0E+02
!     CORGCT=CSoilOrgM_vr(ielmc,L,NY,NX)*1.0E-04
!     CC=EXP(-3.6733-0.1447*CCLAY_vrT+0.7653*CORGCT)
!     DD=-0.4805-0.1239*CCLAY_vrT+0.2080*CORGCT
!     EE=3.8521+0.0963*CCLAY_vrT
!     SoilResit4RootPentrate_vr(L,NY,NX)=CC*THETW_vr(L,NY,NX)**DD*SoilBulkDensity_vr(L,NY,NX)**EE
!     ELSE
    SoilResit4RootPentrate_vr(L,NY,NX)=0.0_r8
!     ENDIF
!
!     SOIL HYDRAULIC CONDUCTIVITIES FROM AMBIENT SOIL WATER CONTENTS
!
!     HydroCondMicP4RootUptake_vr=soil hydraulic conductivity for root uptake
!
    K=MAX(1,MIN(100,INT(100.0_r8*(POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX))/POROS_vr(L,NY,NX))+1))
    HydroCondMicP4RootUptake_vr(L,NY,NX)=0.5_r8*(HydroCond_3D(1,K,L,NY,NX)+HydroCond_3D(3,K,L,NY,NX))
    
  END DO
  call PrintInfo('end '//subname)
  end subroutine GetSoilHydraulicVars


!------------------------------------------------------------------------------------------

  subroutine SoilHydroProperty(I,J,L,NY,NX)
  !
  !Set up soil hydraulic property
  implicit none
  integer, intent(in) :: L,NY,NX
  integer, intent(in) :: I,J
  real(r8) :: THETF
  real(r8) :: H2OSOIatK(100),PSISK(0:100)

  integer :: K,M,N
  real(r8) :: XK,YK,SUM1,SUM2
  real(r8) :: VISCWL

  IF(CSoilOrgM_vr(ielmc,L,NY,NX).GT.FORGC)THEN
    SRP_vr(L,NY,NX)=0.25_r8
  ELSE IF(CSoilOrgM_vr(ielmc,L,NY,NX).GT.0.5_r8*FORGC)THEN
    SRP_vr(L,NY,NX)=0.33_r8
  ELSE
    SRP_vr(L,NY,NX)=1.00_r8
  ENDIF
  
  if(lverb)write(*,*)'SoilHydroProperty::setshape',POROS_vr(L,NY,NX),cold_run()
! double check cold_run() setup
  LOGPOROS_vr(L,NY,NX)=LOG(POROS_vr(L,NY,NX))

  IF((ISOIL(isoi_fc,L,NY,NX).EQ.isoi_set .AND. ISOIL(isoi_wp,L,NY,NX).EQ.isoi_set) .OR. (.not.cold_run()))THEN
  ! read from check point file or if soil properties are set with soil file    
    LOGFldCapacity_vr(L,NY,NX) = LOG(FieldCapacity_vr(L,NY,NX))
    LOGWiltPoint_vr(L,NY,NX)   = LOG(WiltPoint_vr(L,NY,NX))
    PSD_vr(L,NY,NX)               = LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX)
    FCD_vr(L,NY,NX)               = LOGFldCapacity_vr(L,NY,NX)-LOGWiltPoint_vr(L,NY,NX)
  ELSE
    !
    !     DEFAULT SOIL HYDROLOGIC PPTYS (FIELD CAPACITY, WILTING POINT)
    !     IF ACTUAL VALUES WERE NOT INPUT TO THE SOIL FILE
    !
    !     THW,THI=initial soil water,ice content from soil file
    !
    IF(cold_run())THEN
      call SetColdRunSoilStates(I,J,L,NY,NX)
    ENDIF
  ENDIF

  if(lverb)write(*,*)'finish soilp set'
  VLsoiAirP_vr(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)) &
    +AZMAX1(VLMacP_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX))

  IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    !ratio of total air-filled pore to micropore
    ThetaAir_vr(L,NY,NX)=VLsoiAirP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)
  ELSE
    ThetaAir_vr(L,NY,NX)=0.0_r8
  ENDIF

  IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
    SoilWatAirDry_vr(L,NY,NX)=EXP((LOGPSIFLD(NY,NX)-LOG(-PSIHY))*FCD_vr(L,NY,NX)/LOGPSIMND(NY,NX)+LOGFldCapacity_vr(L,NY,NX))
  ELSE
    SoilWatAirDry_vr(L,NY,NX)=ZERO2
  ENDIF
      !
      !     SATURATED HYDRAULIC CONDUCTIVITY FROM SWC AT SATURATION VS.
      !     -0.033 MPA (MINERAL SOILS) IF NOT ENTERED IN SOIL FILE IN 'READS'
      !
      !     SCNV,SCNH=vertical,lateral saturated hydraulic conductivity
  !
  IF(ISOIL(isoi_scnv,L,NY,NX).EQ.isoi_unset)THEN
    !computing vertical saturated hydraulic conductivity
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
        *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD(NY,NX)+LOGPOROS_vr(L,NY,NX)))
      SatHydroCondVert_vr(L,NY,NX)=1.54_r8*((POROS_vr(L,NY,NX)-THETF)/THETF)**2
    ELSE
      SatHydroCondVert_vr(L,NY,NX) = 0.10_r8+75.0_r8*1.0E-15_r8**SoilBulkDensity_vr(L,NY,NX)
      SatHydroCondVert_vr(L,NY,NX) = SatHydroCondVert_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
    ENDIF
  ENDIF

  IF(ISOIL(isoi_scnh,L,NY,NX).EQ.isoi_unset)THEN
    !computing horizontal saturated hydraulic conductivity
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
        *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD(NY,NX)+LOGPOROS_vr(L,NY,NX)))
      SatHydroCondHrzn_vr(L,NY,NX)=1.54_r8*((POROS_vr(L,NY,NX)-THETF)/THETF)**2._r8
    ELSE
      SatHydroCondHrzn_vr(L,NY,NX) = 0.10_r8+75.0_r8*1.0E-15_r8**SoilBulkDensity_vr(L,NY,NX)
      SatHydroCondHrzn_vr(L,NY,NX) = SatHydroCondHrzn_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
    ENDIF
  ENDIF

  !
  !     HYDRAULIC CONDUCTIVITY FUNCTION FROM KSAT AND SOIL WATER RELEASE CURVE
  !
  !     THETK,PSISK=micropore class water content,potential
  !     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
  !
  !     IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
  SUM2=0.0_r8
  DO  K=1,100
    XK           = K-1
    H2OSOIatK(K) = POROS_vr(L,NY,NX)-(XK/100.0_r8*POROS_vr(L,NY,NX))
    IF(H2OSOIatK(K).LT.FieldCapacity_vr(L,NY,NX))THEN
      PSISK(K)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX)+((LOGFldCapacity_vr(L,NY,NX)-LOG(H2OSOIatK(K))) &
        /FCD_vr(L,NY,NX)*LOGPSIMND(NY,NX))))
    ELSEIF(H2OSOIatK(K).LT.POROS_vr(L,NY,NX)-DTHETW)THEN
      !almost saturated
      PSISK(K)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(L,NY,NX)-LOG(H2OSOIatK(K))) &
        /PSD_vr(L,NY,NX))**SRP_vr(L,NY,NX)*LOGPSIMXD(NY,NX)))
    ELSE
      !fully saturated
      PSISK(K)=PSISE_vr(L,NY,NX)
    ENDIF
    SUM2=SUM2+(2*K-1)/(PSISK(K)**2_r8)
  ENDDO

  DO  K=1,100
    SUM1 = 0.0_r8
    XK   = K-1
    YK   = ((100.0_r8-XK)/100.0_r8)**1.33_r8
    DO M    = K, 100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2_r8)
    ENDDO

    DO  N=1,3
      IF(N.EQ.3)THEN
        !vertical
        HydroCond_3D(N,K,L,NY,NX)=SatHydroCondVert_vr(L,NY,NX)*YK*SUM1/SUM2
        IF(K.GT.1.AND.PSISK(K).LT.PSISoilAirEntry(L,NY,NX))then
          if(PSISK(K-1).GE.PSISoilAirEntry(L,NY,NX))THEN
            !moisture at air-entry saturation
            ThetaSat_vr(L,NY,NX)=H2OSOIatK(K)
          endif
        ENDIF
      ELSE
        !horizontal
        HydroCond_3D(N,K,L,NY,NX)=SatHydroCondHrzn_vr(L,NY,NX)*YK*SUM1/SUM2
      ENDIF
    ENDDO
  ENDDO

!     SOIL MACROPORE DIMENSIONS AND CONDUCTIVITY FROM MACROPORE FRACTION
!     ENTERED IN 'READS'
!
!     PathLenMacP,MacPNumLayer,MacPRadius=path length between, number,radius of macropores
!     CNDH=macropore hydraulic conductivity
!
  MacPRadius(L,NY,NX)   = 0.5E-03_r8
  MacPNumLayer(L,NY,NX) = INT(VLMacP_vr(L,NY,NX)/(PICON*MacPRadius(L,NY,NX)**2._r8*VGeomLayert0_vr(L,NY,NX)))

  IF(MacPNumLayer(L,NY,NX).GT.0.0_r8)THEN
    PathLenMacP(L,NY,NX)=1.0_r8/(SQRT(PICON*MacPNumLayer(L,NY,NX)))
  ELSE
    PathLenMacP(L,NY,NX)=1.0_r8
  ENDIF
  VISCWL                    = VISCW*EXP(0.533_r8-0.0267_r8*TCS_vr(L,NY,NX))
  HydroCondMacP_vr(L,NY,NX) = 3.6E+03_r8*PICON*MacPNumLayer(L,NY,NX)*MacPRadius(L,NY,NX)**4._r8/(8.0_r8*VISCWL)
!  write(*,*)L,SatHydroCondVert_vr(L,NY,NX),HydroCondMacP_vr(L,NY,NX)
  end subroutine SoilHydroProperty

!------------------------------------------------------------------------------------------
  subroutine SetColdRunSoilStates(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  
! restart is defined as simulation starting from a previous run
  IF(ISOIL(isoi_fc,L,NY,NX).EQ.isoi_unset .OR. ISOIL(isoi_wp,L,NY,NX).EQ.isoi_unset)THEN
    !calculating FC or WP
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      FieldCapacity_vr(L,NY,NX)=0.2576_r8-0.20_r8*CSAND_vr(L,NY,NX) &
            +0.36_r8*CCLAY_vr(L,NY,NX)+0.60E-06*CSoilOrgM_vr(ielmc,L,NY,NX)
    ELSE
      IF(SoilBulkDensity_vr(L,NY,NX).LT.0.075_r8)THEN
        FieldCapacity_vr(L,NY,NX)=0.27_r8
      ELSEIF(SoilBulkDensity_vr(L,NY,NX).LT.0.195_r8)THEN
        FieldCapacity_vr(L,NY,NX)=0.62_r8
      ELSE
        FieldCapacity_vr(L,NY,NX)=0.71_r8
      ENDIF
    ENDIF
    FieldCapacity_vr(L,NY,NX) = FieldCapacity_vr(L,NY,NX)/(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
    FieldCapacity_vr(L,NY,NX) = AMIN1(0.75_r8*POROS_vr(L,NY,NX),FieldCapacity_vr(L,NY,NX))
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      WiltPoint_vr(L,NY,NX)=0.0260_r8+0.50_r8*CCLAY_vr(L,NY,NX)+0.32E-06_r8*CSoilOrgM_vr(ielmc,L,NY,NX)
    ELSE
      IF(SoilBulkDensity_vr(L,NY,NX).LT.0.075_r8)THEN
        WiltPoint_vr(L,NY,NX)=0.04_r8
      ELSEIF(SoilBulkDensity_vr(L,NY,NX).LT.0.195_r8)THEN
        WiltPoint_vr(L,NY,NX)=0.15_r8
      ELSE
        WiltPoint_vr(L,NY,NX)=0.22_r8
      ENDIF
    ENDIF
    WiltPoint_vr(L,NY,NX) = WiltPoint_vr(L,NY,NX)/(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
    WiltPoint_vr(L,NY,NX) = AMIN1(0.75_r8*FieldCapacity_vr(L,NY,NX),WiltPoint_vr(L,NY,NX))
  ENDIF
  LOGFldCapacity_vr(L,NY,NX) = LOG(FieldCapacity_vr(L,NY,NX))
  LOGWiltPoint_vr(L,NY,NX)   = LOG(WiltPoint_vr(L,NY,NX))
  PSD_vr(L,NY,NX)               = LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX)
  FCD_vr(L,NY,NX)               = LOGFldCapacity_vr(L,NY,NX)-LOGWiltPoint_vr(L,NY,NX)

!   IBEGIN:   start date of model run

    !IDATA(9): start year of model run
  IF(I.EQ.IBEGIN .AND. J.EQ.1 .AND. is_first_year)THEN
    !first time step at the beginning year
    !THW=initial soil water content
    !DPTH=depth to middle of soil layer [m]
    !ExtWaterTablet0_col=external water table depth, [m]
    IF(THW_vr(L,NY,NX).GT.1.0_r8 .OR. SoilDepthMidLay_vr(L,NY,NX).GE.ExtWaterTablet0_col(NY,NX))THEN
      !below the water table, thus it is saturated
      THETW_vr(L,NY,NX)=POROS_vr(L,NY,NX)
    ELSEIF(isclose(THW_vr(L,NY,NX),1._r8))THEN
      !at field capacity
      THETW_vr(L,NY,NX)=FieldCapacity_vr(L,NY,NX)
    ELSEIF(isclose(THW_vr(L,NY,NX),0._r8))THEN
      !at wilting point
      THETW_vr(L,NY,NX)=WiltPoint_vr(L,NY,NX)
    ELSEIF(THW_vr(L,NY,NX).LT.0.0_r8)THEN
      !CO2CompenPoint_nodeetely dry
      THETW_vr(L,NY,NX)=0.0_r8
    ENDIF

    IF(THI_vr(L,NY,NX).GT.1.0_r8.OR.SoilDepthMidLay_vr(L,NY,NX).GE.ExtWaterTablet0_col(NY,NX))THEN
      THETI_vr(L,NY,NX)=AZMAX1(AMIN1(POROS_vr(L,NY,NX),POROS_vr(L,NY,NX)-THW_vr(L,NY,NX)))
    ELSEIF(isclose(THI_vr(L,NY,NX),1._r8))THEN
      THETI_vr(L,NY,NX)=AZMAX1(AMIN1(FieldCapacity_vr(L,NY,NX),POROS_vr(L,NY,NX)-THW_vr(L,NY,NX)))
    ELSEIF(isclose(THI_vr(L,NY,NX),0._r8))THEN
      THETI_vr(L,NY,NX)=AZMAX1(AMIN1(WiltPoint_vr(L,NY,NX),POROS_vr(L,NY,NX)-THW_vr(L,NY,NX)))
    ELSEIF(THI_vr(L,NY,NX).LT.0.0_r8)THEN
      THETI_vr(L,NY,NX)=0.0_r8
    ENDIF

  !in a cold run, set it
    VLWatMicP_vr(L,NY,NX)     = THETW_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
    VLWatMicPX_vr(L,NY,NX)    = VLWatMicP_vr(L,NY,NX)
    VLWatMacP_vr(L,NY,NX)     = THETW_vr(L,NY,NX)*VLMacP_vr(L,NY,NX)
    VLiceMicP_vr(L,NY,NX)     = THETI_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
    VLiceMacP_vr(L,NY,NX)     = THETI_vr(L,NY,NX)*VLMacP_vr(L,NY,NX)
    ThetaH2OZ_vr(L,NY,NX)     = THETW_vr(L,NY,NX)
    ThetaICEZ_vr(L,NY,NX)     = THETI_vr(L,NY,NX)
    VHeatCapacity_vr(L,NY,NX) = VHeatCapacitySoilM_vr(L,NY,NX)+Cpw*(VLWatMicP_vr(L,NY,NX) &
      +VLWatMacP_vr(L,NY,NX))+Cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
  ENDIF
  end subroutine SetColdRunSoilStates
!------------------------------------------------------------------------------------------  
  subroutine LitterHydroproperty(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer  :: K,M
  real(r8) :: XK,YK
  real(r8) :: SUM1,SUM2
  integer, parameter :: n100=100
  real(r8) :: H2OSOIatK(n100),PSISK(0:n100)

  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    SoilBulkDensity_vr(0,NY,NX)=VLSoilMicPMass_vr(0,NY,NX)/VGeomLayer_vr(0,NY,NX)
  ELSE
    SoilBulkDensity_vr(0,NY,NX)=BulkDensLitR(micpar%k_fine_litr)
  ENDIF
  SoilWatAirDry_vr(0,NY,NX)=EXP((LOGPSIFLD(NY,NX)-LOG(-PSIHY))*FCD_vr(0,NY,NX)/LOGPSIMND(NY,NX)+LOGFldCapacity_vr(0,NY,NX))
  SUM2             =0.0_r8
  D1220: DO  K=1,n100
    XK           = K-1
    H2OSOIatK(K) = POROS0_col(NY,NX)-(XK/n100*POROS0_col(NY,NX))
    IF(H2OSOIatK(K).LT.FieldCapacity_vr(0,NY,NX))THEN
      PSISK(K)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX)+((LOGFldCapacity_vr(0,NY,NX)-LOG(H2OSOIatK(K))) &
          /FCD_vr(0,NY,NX)*LOGPSIMND(NY,NX))))
    ELSEIF(H2OSOIatK(K).LT.POROS0_col(NY,NX))THEN
      PSISK(K)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(0,NY,NX)-LOG(H2OSOIatK(K))) &
          /PSD_vr(0,NY,NX))**SRP_vr(0,NY,NX)*LOGPSIMXD(NY,NX)))
    ELSE
      PSISK(K)=PSISE_vr(0,NY,NX)
    ENDIF
    SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
  ENDDO D1220

  D1235: DO  K=1,n100
    SUM1  = 0.0_r8
    XK    = K-1
    YK    = ((n100-XK)/n100)**1.33_r8
    D1230: DO M = K, n100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2._r8)
    ENDDO D1230
    HydroCond_3D(3,K,0,NY,NX) = SatHydroCondVert_vr(0,NY,NX)*YK*SUM1/SUM2
    HydroCond_3D(1,K,0,NY,NX) = 0.0_r8
    HydroCond_3D(2,K,0,NY,NX) = 0.0_r8
    IF(K.GT.1.AND.(PSISK(K).LT.PSISoilAirEntry(0,NY,NX)))THEN
      IF(PSISK(K-1).GE.PSISoilAirEntry(0,NY,NX))THEN
        !moisture at air-entry saturation
        ThetaSat_vr(0,NY,NX)=H2OSOIatK(K)
      ENDIF
    ENDIF
  ENDDO D1235
  end subroutine LitterHydroproperty

!------------------------------------------------------------------------------------------

  subroutine ComputeSoilHydroPars(NY,NX,NU,NM)
  implicit none
  integer, intent(in) :: NY,NX,NU,NM
  integer :: L
  DO  L=NU,NM
    IF(FieldCapacity_vr(L,NY,NX).LT.0.0_r8)THEN
      !computing field capacity
      ISOIL(isoi_fc,L,NY,NX)  = isoi_unset
      PSIAtFldCapacity(NY,NX) = -0.033_r8
    ELSE
      ISOIL(isoi_fc,L,NY,NX)=isoi_set
    ENDIF

    IF(WiltPoint_vr(L,NY,NX).LT.0.0_r8)THEN
      !computing wilting point
      ISOIL(isoi_wp,L,NY,NX) = isoi_unset
      PSIAtWiltPoint(NY,NX)  = -1.5_r8
    ELSE
      ISOIL(isoi_wp,L,NY,NX)=isoi_set
    ENDIF
    IF(SatHydroCondVert_vr(L,NY,NX).LT.0.0_r8)THEN
      !soil vertical saturated hydraulic conductivity
      ISOIL(isoi_scnv,L,NY,NX)=isoi_unset
    ELSE
      ISOIL(isoi_scnv,L,NY,NX)=isoi_set
    ENDIF
    IF(SatHydroCondHrzn_vr(L,NY,NX).LT.0.0_r8)THEN
      !soil horizontal saturated hydraulic conductivity
      ISOIL(isoi_scnh,L,NY,NX)=isoi_unset
    ELSE
      ISOIL(isoi_scnh,L,NY,NX)=isoi_set
    ENDIF
  ENDDO
  END subroutine ComputeSoilHydroPars


end module SoilHydroParaMod
