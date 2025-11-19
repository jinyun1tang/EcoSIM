module SoilHydroParaMod

  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar
  use minimathmod,      only: isclose, AZMAX1, AZMIN1, Viscosity_H2O
  Use SoilPhysParaMod,  only: ComputePsiMCM,ComputePSIPond,getMoistK
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
  real(r8), parameter :: FORGW=0.25E+06_r8   !threshold for  C concentration in organic soil 	gC Mg-1
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

  real(r8) :: THETW1           !water-filled soil micropores ![m3 water /(m3 micropore soil)]
  real(r8) :: WPX,WPLX,CC,DD,EE,CCLAY_pct,CORGCT
  integer :: K,L

  ! begin_execution
  call PrintInfo('beg '//subname)

  DO L=NUI_col(NY,NX),NLI_col(NY,NX)
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
      
      call ComputePsiMCM(L,NY,NX,THETW1,PSISoilMatricP_vr(L,NY,NX))

      !not significant soil mass, but have significant ice, for (partially) frozen ponding water 
    ELSE IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

      call ComputePSIPond(L,NY,NX,THETI_vr(L,NY,NX),THETW_vr(L,NY,NX),THETW1,PSISoilMatricP_vr(L,NY,NX))

    ELSE
      PSISoilMatricP_vr(L,NY,NX)=PSISE_vr(L,NY,NX)
    ENDIF
    !
    !     SOIL OSMOTIC, GRAVIMETRIC AND MATRIC WATER POTENTIALS
    !
    !     PSISM,PSISO,PSIGrav_vr,PSIST=matric,osmotic,gravimetric,total water potential
    !
    PSISoilOsmotic_vr(L,NY,NX)          = -RGASC*TKS_vr(L,NY,NX)*SolutesIonConc_vr(L,NY,NX)*1.E-6_r8  !MPa
    PSIGrav_vr(L,NY,NX)                 = mGravAccelerat*(ALT_col(NY,NX)-SoilDepthMidLay_vr(L,NY,NX))
    ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX) = AZMIN1(PSISoilMatricP_vr(L,NY,NX)+PSISoilOsmotic_vr(L,NY,NX)+PSIGrav_vr(L,NY,NX))
    !
    !     SOIL RESISTANCE TO ROOT PENETRATION
    !
    !     SoilResist4RootPentrate_vr=soil resistance to root penetration [MPa]
    !
    IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      !Eq.(8) in Grant (1993), Simulation model of soil compaction and root growth I. Model structure.
      !fitted using data from Rickman et al. (1992)
      if(THETW_vr(L,NY,NX)>0.01_r8)THEN
        CCLAY_pct = CCLAY_vr(L,NY,NX)*1.0E+02_r8
        CORGCT    = CSoilOrgM_vr(ielmc,L,NY,NX)*1.0E-04_r8
        CC        = EXP(-3.6733_r8-0.1447_r8*CCLAY_pct+0.7653_r8*CORGCT)
        DD        = -0.4805_r8-0.1239_r8*CCLAY_pct+0.2080_r8*CORGCT
        EE        = 3.8521_r8+0.0963_r8*CCLAY_pct
        SoilResist4RootPentrate_vr(L,NY,NX) = CC*THETW_vr(L,NY,NX)**DD*SoilBulkDensity_vr(L,NY,NX)**EE
      ELSE
        SoilResist4RootPentrate_vr(L,NY,NX)=  300._r8
      ENDIF
    ELSE
       SoilResist4RootPentrate_vr(L,NY,NX)=0.0_r8
    ENDIF
    !
    !     SOIL HYDRAULIC CONDUCTIVITIES FROM AMBIENT SOIL WATER CONTENTS
    !
    !     HYCDMicP4RootUptake_vr=soil hydraulic conductivity for root uptake
    !
    K=getMoistK(THETW_vr(L,NY,NX),POROS_vr(L,NY,NX))
    HYCDMicP4RootUptake_vr(L,NY,NX)=0.5_r8*(HydroCond_3D(1,K,L,NY,NX)+HydroCond_3D(3,K,L,NY,NX))
    
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
  real(r8) :: H2OSOIK(100),PSISK(0:100)
  character(len=*), parameter :: subname='SoilHydroProperty'

  integer :: K,M,N,N100
  real(r8) :: XK,YK,SUM1,SUM2
  real(r8) :: VISCWL

  call PrintInfo('beg '//subname)
  !gC/Mg soil
  IF(CSoilOrgM_vr(ielmc,L,NY,NX).GT.FORGC)THEN
    !organic soil
    SRP_vr(L,NY,NX)=0.25_r8
  ELSE IF(CSoilOrgM_vr(ielmc,L,NY,NX).GT.0.5_r8*FORGC)THEN
    SRP_vr(L,NY,NX)=0.33_r8
  ELSE
    SRP_vr(L,NY,NX)=1.00_r8
  ENDIF
  
  ! double check cold_run() setup
  LOGPOROS_vr(L,NY,NX)=LOG(POROS_vr(L,NY,NX))

  IF((ISOIL_vr(isoi_fc,L,NY,NX).EQ.isoi_set .AND. ISOIL_vr(isoi_wp,L,NY,NX).EQ.isoi_set) .OR. (.not.cold_run()))THEN
    ! read from check point file or if soil properties are set with soil file    
    LOGFldCapacity_vr(L,NY,NX) = LOG(FieldCapacity_vr(L,NY,NX))
    LOGWiltPoint_vr(L,NY,NX)   = LOG(WiltPoint_vr(L,NY,NX))
    PSD_vr(L,NY,NX)            = LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX)
    FCD_vr(L,NY,NX)            = LOGFldCapacity_vr(L,NY,NX)-LOGWiltPoint_vr(L,NY,NX)
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

  VLsoiAirP_vr(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)) &
    +AZMAX1(VLMacP_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX))

  IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    !ratio of total air-filled pore to micropore
    ThetaAir_vr(L,NY,NX)=VLsoiAirP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)
  ELSE
    ThetaAir_vr(L,NY,NX)=0.0_r8
  ENDIF

  IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
    !PSIHY: soil matric potential at dessication
    SoilWatAirDry_vr(L,NY,NX)=EXP((LOGPSIFLD_col(NY,NX)-LOG(-PSIHY))*FCD_vr(L,NY,NX)/LOGPSIMND_col(NY,NX)+LOGFldCapacity_vr(L,NY,NX))
  ELSE
    SoilWatAirDry_vr(L,NY,NX)=ZERO2
  ENDIF
      !
      !     SATURATED HYDRAULIC CONDUCTIVITY FROM SWC AT SATURATION VS.
      !     -0.033 MPA (MINERAL SOILS) IF NOT ENTERED IN SOIL FILE IN 'READS'
      !
      !     SCNV,SCNH=vertical,lateral saturated hydraulic conductivity
  !
  IF(ISOIL_vr(isoi_scnv,L,NY,NX).EQ.isoi_unset)THEN
    !computing vertical saturated hydraulic conductivity
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
        *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))
      SatHydroCondVert_vr(L,NY,NX)=1.54_r8*((POROS_vr(L,NY,NX)-THETF)/THETF)**2
    ELSE
      SatHydroCondVert_vr(L,NY,NX) = 0.10_r8+75.0_r8*1.0E-15_r8**SoilBulkDensity_vr(L,NY,NX)
      SatHydroCondVert_vr(L,NY,NX) = SatHydroCondVert_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
    ENDIF
  ENDIF

  IF(ISOIL_vr(isoi_scnh,L,NY,NX).EQ.isoi_unset)THEN
    !computing horizontal saturated hydraulic conductivity
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
        *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))
      SatHydroCondHrzn_vr(L,NY,NX)=1.54_r8*((POROS_vr(L,NY,NX)-THETF)/THETF)**2
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
  !using Eq.(1) and II in Table 1 from Green and Corey, 1971
  SUM2=0.0_r8;N100=100
  !integrate of 1/Psi(x) over [0,1] 
  DO  K=1,n100
    !air-filled relative pore volume
    XK  = (K-1._r8)/n100
    !moisture content
    H2OSOIK(K) = POROS_vr(L,NY,NX)*(1._r8-XK)
    call ComputePsiMCM(L,NY,NX,H2OSOIK(K),PSISK(K))
    SUM2=SUM2+(2._r8*K-1)/(PSISK(K)**2)
  ENDDO

  !progress from dry to saturation
  !completely dry at XK==1
  DO  K=1,N100
    
    !air-filled relative pore volume
    XK   = (K-1._r8)/N100
    !saturation at K
    YK   = (1._r8-XK)**1.33_r8
    SUM1 = 0.0_r8
    DO M    = K, N100
      SUM1=SUM1+(2*(M-K)+1._r8)/(PSISK(M)**2)
    ENDDO

    DO  N=1,3
      IF(N.EQ.3)THEN
        !vertical
        HydroCond_3D(N,K,L,NY,NX)=SatHydroCondVert_vr(L,NY,NX)*YK*SUM1/SUM2

        !obtain water content at saturation
        IF(K.GT.1 .AND. PSISK(K).LT.PSISoilAirEntry_vr(L,NY,NX))then
          if(PSISK(K-1).GE.PSISoilAirEntry_vr(L,NY,NX))THEN
            !moisture at air-entry saturation
            ThetaSat_vr(L,NY,NX)=H2OSOIK(K)
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
  !VGeomLayert0_vr: geometric volume of the layer 
  !VLMacP_vr: volume of macropores in the layer
  !Assuming soil macropores are approximated as (cylydrical) pipes, and has layer thickness as pipe length
  !Eq.(3.3)-(3.4) from Dimitrov et al. (2010)
  
  MacPoreRadius_vr(L,NY,NX)   = 0.5E-03_r8
  MacPoreNumbers_vr(L,NY,NX) = INT(VLMacP_vr(L,NY,NX)/(PICON*MacPoreRadius_vr(L,NY,NX)**2*VGeomLayert0_vr(L,NY,NX)))

  IF(MacPoreNumbers_vr(L,NY,NX).GT.0.0_r8)THEN
    PathLenMacPore_vr(L,NY,NX)=1.0_r8/(SQRT(PICON*MacPoreNumbers_vr(L,NY,NX)))
  ELSE
    PathLenMacPore_vr(L,NY,NX)=1.0_r8
  ENDIF
  !water viscosity as a function of temperature
  VISCWL  = Viscosity_H2O(TCS_vr(L,NY,NX))
  !Apply the Poiseuille's Law
  HydroCondMacP_vr(L,NY,NX) = 3.6E+03_r8*PICON*MacPoreNumbers_vr(L,NY,NX)*MacPoreRadius_vr(L,NY,NX)**4/(8.0_r8*VISCWL)

  call PrintInfo('end '//subname)
  end subroutine SoilHydroProperty

!------------------------------------------------------------------------------------------
  subroutine SetColdRunSoilStates(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  
! restart is defined as simulation starting from a previous run
  IF(ISOIL_vr(isoi_fc,L,NY,NX).EQ.isoi_unset .OR. ISOIL_vr(isoi_wp,L,NY,NX).EQ.isoi_unset)THEN
    !calculating FC or WP
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      FieldCapacity_vr(L,NY,NX)=0.2576_r8-0.20_r8*CSAND_vr(L,NY,NX) &
            +0.36_r8*CCLAY_vr(L,NY,NX)+0.60E-06*CSoilOrgM_vr(ielmc,L,NY,NX)
    ELSE
      !organic soil
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
      !organic soil
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
  PSD_vr(L,NY,NX)            = LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX)
  FCD_vr(L,NY,NX)            = LOGFldCapacity_vr(L,NY,NX)-LOGWiltPoint_vr(L,NY,NX)

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
  real(r8) :: H2OSOIK(n100),PSISK(0:n100)

  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    SoilBulkDensity_vr(0,NY,NX)=VLSoilMicPMass_vr(0,NY,NX)/VGeomLayer_vr(0,NY,NX)
  ELSE
    SoilBulkDensity_vr(0,NY,NX)=BulkDensLitR(micpar%k_fine_litr)
  ENDIF
  SoilWatAirDry_vr(0,NY,NX)=EXP((LOGPSIFLD_col(NY,NX)-LOG(-PSIHY))*FCD_vr(0,NY,NX)/LOGPSIMND_col(NY,NX)+LOGFldCapacity_vr(0,NY,NX))
  SUM2             =0.0_r8
  D1220: DO  K=1,n100
    XK         = (K-1._r8)/n100
    H2OSOIK(K) = POROS0_col(NY,NX)*(1._r8-XK)
    call ComputePsiMCM(0,NY,NX,H2OSOIK(K),PSISK(K))
    SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
  ENDDO D1220

  D1235: DO  K=1,n100
    SUM1  = 0.0_r8
    XK    = (K-1._r8)/N100
    !air-filled pore
    YK    = (1._r8-XK)**1.33_r8
    D1230: DO M = K, n100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2)
    ENDDO D1230
    HydroCond_3D(3,K,0,NY,NX) = SatHydroCondVert_vr(0,NY,NX)*YK*SUM1/SUM2
    HydroCond_3D(1,K,0,NY,NX) = 0.0_r8
    HydroCond_3D(2,K,0,NY,NX) = 0.0_r8
    IF(K.GT.1.AND.(PSISK(K).LT.PSISoilAirEntry_vr(0,NY,NX)))THEN
      IF(PSISK(K-1).GE.PSISoilAirEntry_vr(0,NY,NX))THEN
        !moisture at air-entry saturation
        ThetaSat_vr(0,NY,NX)=H2OSOIK(K)
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
      ISOIL_vr(isoi_fc,L,NY,NX)  = isoi_unset
      PSIAtFldCapacity_col(NY,NX) = -0.033_r8
    ELSE
      ISOIL_vr(isoi_fc,L,NY,NX)=isoi_set
    ENDIF

    IF(WiltPoint_vr(L,NY,NX).LT.0.0_r8)THEN
      !computing wilting point
      ISOIL_vr(isoi_wp,L,NY,NX) = isoi_unset
      PSIAtWiltPoint_col(NY,NX)  = -1.5_r8
    ELSE
      ISOIL_vr(isoi_wp,L,NY,NX)=isoi_set
    ENDIF
    IF(SatHydroCondVert_vr(L,NY,NX).LT.0.0_r8)THEN
      !soil vertical saturated hydraulic conductivity
      ISOIL_vr(isoi_scnv,L,NY,NX)=isoi_unset
    ELSE
      ISOIL_vr(isoi_scnv,L,NY,NX)=isoi_set
    ENDIF
    IF(SatHydroCondHrzn_vr(L,NY,NX).LT.0.0_r8)THEN
      !soil horizontal saturated hydraulic conductivity
      ISOIL_vr(isoi_scnh,L,NY,NX)=isoi_unset
    ELSE
      ISOIL_vr(isoi_scnh,L,NY,NX)=isoi_set
    ENDIF
  ENDDO
  END subroutine ComputeSoilHydroPars


end module SoilHydroParaMod
