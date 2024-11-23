module SurfPhysMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : endrun
  use GridDataType
  use HydroThermData
  use CanopyDataType
  use SoilPropertyDataType
  use ClimForcDataType
  use SurfSoilDataType
  use LandSurfDataType
  use SoilWaterDataType
  use SoilHeatDataType
  use SurfLitterDataType
  use MiniFuncMod
  use minimathmod
  use SurfPhysData
  use EcoSIMCtrlDataType 
  use SnowDataType
  use SOMDataType
  USE ChemTranspDataType
  use HydroThermData
  use SnowPhysMod
  use PlantTraitDataType
  use SoilPhysDataType
  use SoilBGCDataType
  use IrrigationDataType
  use SnowPhysData
  use PhysPars
  use SurfLitterPhysMod
  use EcoSIMSolverPar
  use EcoSIMCtrlMod
  use SoilPhysParaMod
  use SnowBalanceMod
implicit none
  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  real(r8), parameter :: tiny_wat=1.e-13_r8

  !surface model
  public :: StageSurfacePhysModel
  public :: RunSurfacePhysModel
  public :: SurfaceEnergyModel
  public :: AggregateSurfRunoffFlux  
  public :: writeSurfDiagnosis
  !
  public :: SurfaceRunoff
  public :: UpdateSurfaceAtM

! FEnergyImpact4Erosion=rate constant for restoring surface Ksat
! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
! RACX,LitRSurfResistance=minimum boundary layer resistances of canopy,litter (h m-1)

  real(r8), parameter :: FEnergyImpact4Erosion = 1.0E-03_r8
  real(r8), parameter :: SatHydroCondLitR      = 25.0_r8      !saturated hydraulic conductivity of surface litter
  real(r8), parameter :: LitRSurfResistance    = 0.0139_r8    !minimum boundary layer resistances of litter [h m-1]
  real(r8), parameter :: SoilEmisivity         = 0.97_r8      !soil emissivity
  real(r8), parameter :: SnowEmisivity         = 0.97_r8      !snowpack emissivity
  real(r8), parameter :: SurfLitREmisivity     = 0.97_r8      !surfce litter emissivity
  real(r8), parameter :: RACX                  = 0.0139_r8    !total canopy boundary later resistance h/m

  real(r8) :: PrecHeat2LitR2,Prec2LitR2  
  real(r8) :: HeatSensVapAir2Soi
  real(r8) :: HeatSensAir2Grnd
  real(r8) :: Radnet2LitGrnd   !net radiation onto ground [MJ]
  real(r8) :: LatentHeatEvapAir2Grnd,NetWatFlx2SoiMacP
  real(r8) :: CumHeatSensAir2LitR
  real(r8) :: cumHeatSensAir2Soil,NetWatFlx2LitR
  real(r8) :: WatNetFlo2TopSoiMicP
  real(r8) :: PrecHeat2SoiNet,PrecNet2SoiMicP,PrecAir2LitR
  real(r8) :: PrecHeatAir2LitR,PrecNet2SoiMacP  
contains


  subroutine StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  use SnowPhysMod, only : CopySnowStates
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8),dimension(:,:),intent(OUT) :: ResistanceLitRLay(JY,JX)
  character(len=*), parameter :: subn=trim(mod_filename)//'::StageSurfacePhysModel'
  integer :: NY,NX

  watflw =0._r8;waticefl=0._r8

  if(ATS_cpl_mode) then 
    DO NX=NHW,NHE
      DO NY=NVN,NHE  
         NUM(NY,NX)=1 
      enddo
    enddo
  endif

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
    !reset accumulators to zero
    TEvapXAir2Toplay_col(NY,NX)       = 0._r8
    TEvapXAir2LitR_col(NY,NX)         = 0._r8
    TEvapXAir2Snow_col(NY,NX)         = 0._r8
    TXGridSurfRunoff_2DH(NY,NX)       = 0.0_r8
    THeatXGridBySurfRunoff_2DH(NY,NX) = 0.0_r8

    !make a local copy of the upper boundary index
!
!     ADJUST SURFACE ELEVATION USED IN RUNOFF FOR FREEZE-THAW, EROSION
!     AND SOC
!
!     Altitude_grid,ALT=current,initial elevation of ground surface
!     CumDepz2LayerBot_vr(NUM(NY,NX)-1),=depth of ground surface
!     EnergyImpact4Erosion=cumulative rainfall energy impact on soil surface
!
      Altitude_grid(NY,NX)        = ALT(NY,NX)-CumDepz2LayerBot_vr(NUM(NY,NX)-1,NY,NX)
      EnergyImpact4Erosion(NY,NX) = EnergyImpact4Erosion(NY,NX)*(1.0_r8-FEnergyImpact4Erosion)

      call CopySnowStates(I,J,NY,NX)
      
      call CopySurfaceVars(I,J,NY,NX)
      
      call PartionSurfaceFraction(NY,NX)

      call PartitionPrecip(I,J,NY,NX)
      
      call SurfaceResistances(I,J,NY,NX,ResistanceLitRLay(NY,NX))

      call SurfaceRadiation(I,J,NY,NX)

      call SetCanopyProperty(NY,NX)

    ENDDO D9990
  ENDDO D9995
  end subroutine StageSurfacePhysModel

!------------------------------------------------------------------------------------------  

  subroutine CopySurfaceVars(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8) :: VWatLitrZ,TVWatIceLitR,VOLIRZ
!
! SET INITIAL SOIL VALUES
!
! LitrIceFlxThaw,LitrIceHeatFlxFrez=initialize surface litter freeze,thaw,latent heat
  LitrIceFlxThaw_col(NY,NX)     = 0.0_r8
  LitrIceHeatFlxFrez_col(NY,NX) = 0.0_r8

!
! ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
!     FOR USE AT INTERNAL TIME STEP IN SURFACE LITTER
!
!     LWRadBySurf_col=longwave emission from litter surface
!     VHeatCapacity1_vr=volumetric heat capacity of litter
!     VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of litter
!     VWatLitRHoldCapcity=maximum water retention by litter
!     XVOLT,XVOLW=free surface water+ice,water
!     VHeatCapLitR=min heat capacity for litter water,heat fluxes
!     VLitR=litter volume
!     THETW*,THETI*,THETP*=water,ice,air concentrations
!     PSISM*=litter matric water potential
!
  LWRadBySurf_col(NY,NX)        = 0.0_r8
  VHeatCapacity1_vr(0,NY,NX)    = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
  VLPoreLitR_col(NY,NX)         = VLMicP_vr(0,NY,NX)
  VLWatMicP1_vr(0,NY,NX)        = AZMAX1(VLWatMicP_vr(0,NY,NX))
  VLiceMicP1_vr(0,NY,NX)        = AZMAX1(VLiceMicP_vr(0,NY,NX))
  VLairMicP1_vr(0,NY,NX)        = AZMAX1(VLPoreLitR_col(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))
  VLWatMicPM_vr(1,0,NY,NX)      = VLWatMicP1_vr(0,NY,NX)
  VLsoiAirPM(1,0,NY,NX)         = VLairMicP1_vr(0,NY,NX)
  TVWatIceLitR                  = VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX)
  XVLMobileWaterLitR_col(NY,NX) = AZMAX1(TVWatIceLitR-VWatLitRHoldCapcity_col(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ               = VLWatMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    VOLIRZ                  = VLiceMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    XVLMobileWatMicP(NY,NX) = AZMAX1(VLWatMicP1_vr(0,NY,NX)-VWatLitrZ)
    XVLiceMicP_col(NY,NX)   = AZMAX1(VLiceMicP1_vr(0,NY,NX)-VOLIRZ)
  ELSE
    XVLMobileWatMicP(NY,NX)=0.0_r8
    XVLiceMicP_col(NY,NX)=0.0_r8
  ENDIF
  
  XVLMobileWaterLitRM(1,NY,NX) = XVLMobileWaterLitR_col(NY,NX)
  XVLMobileWatMicPM(1,NY,NX)   = XVLMobileWatMicP(NY,NX)
  XVLiceMicPM(1,NY,NX)         = XVLiceMicP_col(NY,NX)
  
  IF(VLitR_col(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat_vr(0,NY,NX)     = AZMAX1t(VLWatMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoiPAsIce_vr(0,NY,NX)     = AZMAX1t(VLiceMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoilPoreAsAir_vr(0,NY,NX) = AZMAX1t(VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)) * &
      AZMAX1t((1.0_r8-XVLMobileWaterLitR_col(NY,NX)/VLWatheldCapSurf_col(NY,NX)))
  ELSE
    FracSoiPAsWat_vr(0,NY,NX)     = 0.0_r8
    FracSoiPAsIce_vr(0,NY,NX)     = 0.0_r8
    FracSoilPoreAsAir_vr(0,NY,NX) = 1.0
  ENDIF
!  write(145,*)I+J/24.,VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX),VLWatMicP1_vr(0,NY,NX),VLiceMicP1_vr(0,NY,NX),&
!    VLairMicP1_vr(0,NY,NX),VLitR_col(NY,NX)
  THETPM(1,0,NY,NX)   = FracSoilPoreAsAir_vr(0,NY,NX)
  PSISM1_vr(0,NY,NX)  = PSISoilMatricP_vr(0,NY,NX)
  TKSoil1_vr(0,NY,NX) = TKS_vr(0,NY,NX)

  end subroutine CopySurfaceVars

!------------------------------------------------------------------------------------------

  subroutine PartionSurfaceFraction(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

!     SNOW AND RESIDUE COVERAGE OF SOIL SURFACE
!     FSNW,FSNX=fractions of snow,snow-free cover
!     SnowDepth=snowpack depth
!     MinSnowDepth=minimum snowpack depth for full cover
!     BARE,CVRD=fractions of soil,litter cover

  FracSurfAsSnow_col(NY,NX)  = AMIN1(1.0_r8,SQRT((SnowDepth_col(NY,NX)/MinSnowDepth)))
  FracSurfSnoFree_col(NY,NX) = 1.0_r8-FracSurfAsSnow_col(NY,NX)
  !if there is heat-wise significant litter layer
  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    FracSurfBareSoil_col(NY,NX)=AMIN1(1.0_r8,AZMAX1(EXP(-0.8E-02_r8*(SoilOrgM_vr(ielmc,0,NY,NX)/AREA(3,0,NY,NX)))))
  ELSE
    FracSurfBareSoil_col(NY,NX)=1.0_r8
  ENDIF

  FracSurfByLitR_col(NY,NX)=1.0_r8-FracSurfBareSoil_col(NY,NX)
  end subroutine PartionSurfaceFraction

!------------------------------------------------------------------------------------------
  subroutine SetCanopyProperty(NY,NX)      
  
  implicit none
  integer, intent(in) :: NY,NX

  real(r8), parameter :: SpecHeatCapAir=1.25E-03_r8    !heat capacity of 1m3 air, [MJ/(m3 K)]

  !TLEX=total latent heat flux x boundary layer resistance, [MJ m-1]
  !TSHX=total sensible heat flux x boundary layer resistance, [MJ m-1]
  !VPQ_col=vapor pressure in canopy air, 
  !TKQ=temperature in canopy air, Kelvin

  !write(*,*) "For VPQ_col and TKQ calc:"
  !write(*,*) "VPA(NY,NX) = ", VPA(NY,NX)
  !write(*,*) "TLEX_col(NY,NX) = ", TLEX_col(NY,NX)
  !write(*,*) "EvapLHTC = ", EvapLHTC
  !write(*,*) "NY = ", NY
  !write(*,*) "NX = ", NX
  !write(*,*) "NUM(NY,NX) = ", NUM(NY,NX)
  !write(*,*) "AREA(3,NUM(NY,NX),NY,NX) = ", AREA(3,NUM(NY,NX),NY,NX)
  !write(*,*) "TairK_col(NY,NX) = ", TairK_col(NY,NX)
  !write(*,*) "TSHX_col(NY,NX) = ", TSHX_col(NY,NX)
  !write(*,*) "SpecHeatCapAir = ", SpecHeatCapAir

  VPQ_col(NY,NX) = VPA(NY,NX)-TLEX_col(NY,NX)/(EvapLHTC*AREA(3,NUM(NY,NX),NY,NX))
  TKQ_col(NY,NX) = TairK_col(NY,NX)-TSHX_col(NY,NX)/(SpecHeatCapAir*AREA(3,NUM(NY,NX),NY,NX))

  end subroutine SetCanopyProperty
!------------------------------------------------------------------------------------------

  subroutine SurfaceRadiation(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: THRYX   !long-wave radiation incident on ground [MJ]
  real(r8) :: RADGX   !short-wave radiation incident on ground, [MJ]
  real(r8) :: Stefboltz_area
!
!     INITIALIZE PARAMETERS, FLUXES FOR ENERGY EXCHANGE
!     AT SNOW, RESIDUE AND SOIL SURFACES
!
!     RADGX=shortwave radiation at ground surface
!     RadSWonSno,RadSWonSoil_col,RadSWonLitR_col= shortwave radn at snowpack,soil,litter, [MJ]
!     FracSWRad2Grnd_col=fraction of shortwave radiation at ground surface
!     FSNW,FSNX=fractions of snow,snow-free cover
!     BARE,CVRD=fractions of soil,litter cover
!     XNPS=internal time step for fluxes through snowpack
!     THRYX=longwave radiation at ground surface
!     LWRad2Snow_col,LWRad2Soil_col,LWRad2LitR_col=longwave radn incident at snowpack,soil,litter
!     LWEmscefSnow_col,LWEmscefSoil_col,LWEmscefLitR_col=longwave radn emitted by snowpack,soil,litter
!     SnowEmisivity,SoilEmisivity,SurfLitREmisivity=emissivity of snowpack,soil,litter surfaces
!     THS=sky longwave radiation
!     LWRadCanGPrev_col=longwave radiation emitted by canopy

  RADGX                  = RadSWGrnd_col(NY,NX)*dts_HeatWatTP
  RadSWonSno_col(NY,NX)  = RADGX*FracSurfAsSnow_col(NY,NX)*XNPS
  RadSWonSoil_col(NY,NX) = RADGX*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)
  RadSWonLitR_col(NY,NX) = RADGX*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR  

  THRYX                 = (LWRadSky_col(NY,NX)*FracSWRad2Grnd_col(NY,NX)+LWRadCanGPrev_col(NY,NX))*dts_HeatWatTP
  LWRad2Snow_col(NY,NX) = THRYX*FracSurfAsSnow_col(NY,NX)*XNPS
  LWRad2Soil_col(NY,NX) = THRYX*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)
  LWRad2LitR_col(NY,NX) = THRYX*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR  

  ! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
  !stefboltz_const is stefan-boltzman constant converted into [MJ /(m^2 K^4 h)]
  Stefboltz_area          = stefboltz_const*AREA(3,NUM(NY,NX),NY,NX)
  LWEmscefSnow_col(NY,NX) = SnowEmisivity*Stefboltz_area*FracSurfAsSnow_col(NY,NX)*dts_HeatWatTP
  LWEmscefSoil_col(NY,NX) = SoilEmisivity*Stefboltz_area*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)*dts_HeatWatTP
  LWEmscefLitR_col(NY,NX) = SurfLitREmisivity*Stefboltz_area*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dts_HeatWatTP
!
  end subroutine SurfaceRadiation
!------------------------------------------------------------------------------------------
  subroutine SurfaceResistances(I,J,NY,NX,ResistanceLitRLay)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8), intent(out):: ResistanceLitRLay     !litter layer resistance [h/m]
  real(r8) :: FracSoiPAsAir0,DFVR,TFACR  
  real(r8) :: PAREX  !area scaled conductances for latent heat fluxes, m^2/h
  real(r8) :: PARSX  !area scaled conductances for sensible heat fluxes, m^2/h
  real(r8) :: RAS    !Snow layer resistance, m
  real(r8) :: ALFZ,WindSpeedGrnd
!
!     BOUNDARY LAYER CONDUCTANCES FOR EXPORT TO TranspNoSalt.F
!

!     AERODYNAMIC RESISTANCE OF CANOPY TO SNOW/RESIDUE/SOIL
!     SURFACE ENERGY EXCHANGE WITH ATMOSPHERE
!
!     ALFZ=parameter for canopy effect on windspeed
!     FracSWRad2Grnd_col=fraction of shortwave radiation at ground surface
!     AbvCanopyBndlResist_col,RAC=isothermal blr above canopy, canopy blr, [h/m]
!     ZT,SoilSurfRoughnesst0_col=canopy, surface roughness heights
!     UA,WindSpeedGrnd=windspeeds above,below canopy
!     VPQ_col,VPA=vapor pressure within,above canopy
!     TKQ,TairK=temperature within,above canopy
!     TLEX,TSHX=net latent,sensible heat fluxes x blrs from prev hour
!     VAP=latent heat of evaporation
!     1.25E-03=heat capacity of air
!     AREA=surface area of grid cell
!
  ALFZ=2.0_r8*(1.0_r8-FracSWRad2Grnd_col(NY,NX))
  IF(AbvCanopyBndlResist_col(NY,NX).GT.ZERO .AND. CanopyHeight_col(NY,NX).GT.SoilSurfRoughnesst0_col(NY,NX) &
    .AND. ALFZ.GT.ZERO)THEN
    !If there are plants
    CanopyBndlResist_col(NY,NX)=AMIN1(RACX,AZMAX1(CanopyHeight_col(NY,NX)*EXP(ALFZ) &
      /(ALFZ/AbvCanopyBndlResist_col(NY,NX))*AZMAX1(EXP(-ALFZ*SoilSurfRoughnesst0_col(NY,NX)/CanopyHeight_col(NY,NX)) &
      -EXP(-ALFZ*(ZERO4PlantDisplace_col(NY,NX)+RoughHeight_col(NY,NX))/CanopyHeight_col(NY,NX)))))
    WindSpeedGrnd=WindSpeedAtm_col(NY,NX)*EXP(-ALFZ)
  ELSE
    CanopyBndlResist_col(NY,NX) = 0.0_r8
    WindSpeedGrnd               = WindSpeedAtm_col(NY,NX)
  ENDIF
!
!     AERODYNAMIC RESISTANCE OF SNOWPACK, RESIDUE AND SOIL
!     SURFACES TO ENERGY EXCHANGE WITH ATMOSPHERE
!     Soil Sci. Soc. Am. J. 48:25-32
!
!     RAR=porosity-unlimited litter blr
!     DLYRR=litter depth
!     VaporDiffusivityLitR_col=vapor diffusivity in litter
!     RAG,RAGW,RAGR=isothermal blrs at ground,snowpack,litter surfaces
!     LitRSurfResistance=boundary layer resistance (blr) of litter surface
!     FracSoiPAsAir*=air-filled porosity of litter
!     DFVR=porosity limitation to diffusion through litter
!     POROQ=litter tortuosity
!     ResistanceLitRLay=porosity-limited litter blr [h/m]
!     PAREX,PARSX=
!     PAREW,PARSW=conductances for snowpack latent,sensible heatfluxes
!     PAREG,PARSG=conductances for soil latent,sensible heat fluxes
!     PARER,PARSR=conductances for litter latent,sensible heat fluxes
!     XNPR=internal time step for fluxes through litter
!     dts_HeatWatTP=time step for heat calculation in watsub 1/h 

! Note: VaporDiffusivityLitR_col is computed in Hour1 which is not used in the
!       ATS coupler. However the calculation is simple so I'm just reproducing it
!       here. This depends on the temperature of the liter (TKS) which is the same
!       as TKSoi1 at the surface
  TFACR                           = TEFGASDIF(TKS_vr(0,NY,NX))
  VaporDiffusivityLitR_col(NY,NX) = TFACR*7.70E-02_r8
  VapDiffusResistanceLitR(NY,NX)  = DLYRR_COL(NY,NX)/VaporDiffusivityLitR_col(NY,NX)

  ResistAreodynOverSoil_col(NY,NX) = CanopyBndlResist_col(NY,NX)+AbvCanopyBndlResist_col(NY,NX)
  ResistAreodynOverSnow_col(NY,NX) = ResistAreodynOverSoil_col(NY,NX)
  ResistAreodynOverLitr_col(NY,NX) = ResistAreodynOverSoil_col(NY,NX)+LitRSurfResistance
  RARG(NY,NX)                      = ResistAreodynOverLitr_col(NY,NX)
  FracSoiPAsAir0                   = AMAX1(ZERO2,FracSoilPoreAsAir_vr(0,NY,NX))
  DFVR                             = FracSoiPAsAir0*POROQ*FracSoiPAsAir0/POROS_vr(0,NY,NX)
  ResistanceLitRLay                = ResistAreodynOverSoil_col(NY,NX)+VapDiffusResistanceLitR(NY,NX)/DFVR
  PAREX                            = AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP               !conductance for latent heat flux, m^2 h
  ! assuming volumetric heat capacity of air is 1.25E-3_r8 MJ/(m^3 K)
  PARSX                            = 1.25E-03_r8*AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP   !conductance for sensible heat flux, MJ h /(m K)
  
  AScaledCdWOverSnow_col(NY,NX) = PAREX*FracSurfAsSnow_col(NY,NX)*XNPS
  AScaledCdHOverSnow_col(NY,NX) = PARSX*FracSurfAsSnow_col(NY,NX)*XNPS
  AScaledCdWOverSoil_col(NY,NX) = PAREX*FracSurfSnoFree_col(NY,NX)
  AScaledCdWOverLitr_col(NY,NX) = PAREX*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR
  AScaledCdHOverSoil_col(NY,NX) = PARSX*FracSurfSnoFree_col(NY,NX)
  AScaledCdHOverLitr_col(NY,NX) = PARSX*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR
!  write(143,*)I+J/24.,ResistAreodynOverSoil_col(NY,NX),FracSoilPoreAsAir_vr(0,NY,NX),VapDiffusResistanceLitR(NY,NX),DFVR

!     PARR=boundary layer conductance above litter,soil surfaces, (m^2/h)/(h/m)
!
  RAS         = SnowBNDResistance(NY,NX)
  PARR(NY,NX) = AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP/(ResistAreodynOverLitr_col(NY,NX)+RAS)   !this includes snow layer resistance
  end subroutine SurfaceResistances

  !------------------------------------------------------------------------------------------

  subroutine SoilSRFEnerbyBalance(M,I,J,NY,NX,PSISV1,LWRadGrnd,ResistanceLitRLay,TopLayWatVol,&
    VapXAir2TopLay,HeatFluxAir2Soi)
  implicit none
  integer, intent(in) :: M,I,J,NY,NX
  real(r8), intent(out):: PSISV1      !surface soil suction pressure [MPa]
  real(r8), intent(out) :: LWRadGrnd  !outgoing long wave radiation from soil surface [MJ]
  real(r8), intent(inout) :: ResistanceLitRLay
  real(r8), intent(inout) :: TopLayWatVol
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi   !total heat flux into soil after surface energy budget [MJ]
  real(r8) :: RAa            ! [h/m]
  real(r8) :: VaporSoi1      ! soil vapor pressure [ton water/m3]
  real(r8) :: CdSoiEvap      ! scaled conductance for evaporation   [m^3]
  real(r8) :: CdSoiHSens     ! scaled conductance for sensible heat [MJ/K]
  real(r8) :: RAGX           ! resistance [h/m]
  real(r8) :: tRadIncid      !total incoming radiation to soil or liter surface, short + long [MJ]
  real(r8) :: RI,WPLX,WPX,FracSoiPAsAir0
  real(r8) :: VLWatGrnd,VLIceGrnd,DFVR
  real(r8) :: TKX1,THETA1S
  real(r8) :: tHeatAir2Grnd   !residual heat flux into soil from incoming radiation minus sensible and latent heat [MJ]
  real(r8) :: AlbedoGrnd
  real(r8) :: RadSWbySoil     !shortwave radiation absorbed by exposed soil [MJ]
! begin_execution
!
! PHYSICAL AND HYDRAULIC PROPERTIES OF SOIL SURFACE INCLUDING
! AIR AND WATER-FILLED POROSITY, AND WATER POTENTIAL USED IN
! FLUX CALCULATIONS
!
! THETY=current,hygroscopic water concentration
! POROS=porosity
! VLWatMicP1,VLSoilPoreMicP_vrI=volume of micropore water
! VLSoilPoreMicP_vrI=soil volume less macropore,rock
! FC,WP=water contents at field capacity,wilting point
! FCL,WPL=log FC,WP
! FCD,PSD=FCL-WPL,log(POROS)-FCL
! PSISM1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
! LOGPSIMX,LOGPSIMN,LOGPSIAtSat=log water potential at FC,WP,POROS
! PSISD,PSIMD=LOGPSIMX-LOGPSIAtSat,LOGPSIMN-LOGPSIMX
! SRP=parameter for deviation from linear log-log water retention
! function from hour1.f
! PSISO=osmotic potential
! BKVL=bulk soil mass for a given layer 

  call CalcSoilWatPotential(NY,NX,NX,NY,NUM(NY,NX),PSISM1_vr(NUM(NY,NX),NY,NX),THETA1S)

  PSISV1=PSISM1_vr(NUM(NY,NX),NY,NX)+PSISoilOsmotic_vr(NUM(NY,NX),NY,NX)
!  write(109,*)I+J/24.,M,PSISM1_vr(NUM(NY,NX),NY,NX),PSISoilOsmotic_vr(NUM(NY,NX),NY,NX),&
!    VLWatMicP1_vr(NUM(NY,NX),NY,NX),VLSoilMicP_vr(NUM(NY,NX),NY,NX)
!
! IF(NX.EQ.4.AND.NY.EQ.5)THEN
!     WRITE(*,3232)'PSISV1',I,J,M,NX,NY,NUM(NY,NX),PSISV1
!    2,PSISM1_vr(NUM(NY,NX),NY,NX),PSISoilOsmotic_vr(NUM(NY,NX),NY,NX)
!    3,FracSoiPAsWat_vr(NUM(NY,NX),NY,NX),THETW1,POROS_vr(NUM(NY,NX),NY,NX)
!    4,PSL(NUM(NY,NX),NY,NX),LOG(THETW1),PSD(NUM(NY,NX),NY,NX)
!    5,VLWatMicP1_vr(NUM(NY,NX),NY,NX),VLSoilMicP_vr(NUM(NY,NX),NY,NX)
!    5,VLSoilPoreMicP_vr(NUM(NY,NX),NY,NX)
!    5,SRP(NUM(NY,NX),NY,NX)
!3232  FORMAT(A8,6I4,20E14.6)
! ENDIF
!
! SOIL SURFACE ALBEDO, NET RADIATION
!
! VLWatMicP1,VLiceMicP1=water,ice volume in micopores
! VLWatMacP1,VLiceMacP1=water,ice volume in macopores
! AlbedoGrnd,SoilAlbedo=albedo of ground surface,soil
! BKVL=soil mass
! RadSWonSoil_col,LWRad2Soil_col,Radnet2LitGrnd=incoming shortwave,longwave,net radiation
! LWRadGrnd,LWEmscefSoil_col=emitted longwave radiation, emissivity
! TK1=soil temperature
! albedo of water and ice are set to 0.06, and 0.30 respectively
  VLWatGrnd=VLWatMicP1_vr(NUM(NY,NX),NY,NX)+VLWatMacP1_vr(NUM(NY,NX),NY,NX)
  VLIceGrnd=VLiceMicP1_vr(NUM(NY,NX),NY,NX)+VLiceMacP1_vr(NUM(NY,NX),NY,NX)

  IF((VLWatGrnd+VLIceGrnd).GT.ZEROS2(NY,NX))THEN
    !top soil has water or ice
    !ice albedo seems too low.
    !write(*,*) "Albedo recomputed"
    AlbedoGrnd=(SoilAlbedo_col(NY,NX)*VLSoilMicPMass_vr(NUM(NY,NX),NY,NX)+0.06_r8*VLWatGrnd &
      +0.30_r8*VLIceGrnd)/(VLSoilMicPMass_vr(NUM(NY,NX),NY,NX)+VLWatGrnd+VLIceGrnd)
  ELSE
    !write(*,*) "Albedo from soil"
    AlbedoGrnd=SoilAlbedo_col(NY,NX)
  ENDIF
  !absorbed radiation
  !Radnet2LitGrnd=net radiation, after taking out outgoing surface layer radiation  
  !LWRadGrnd=emitted longwave radiation  
  RadSWbySoil          = (1.0_r8-AlbedoGrnd)*RadSWonSoil_col(NY,NX)
  tRadIncid            = RadSWbySoil+LWRad2Soil_col(NY,NX)
  LWRadGrnd            = LWEmscefSoil_col(NY,NX)*TKSoil1_vr(NUM(NY,NX),NY,NX)**4._r8
  Radnet2LitGrnd       = tRadIncid-LWRadGrnd
  Eco_RadSW_col(NY,NX) = Eco_RadSW_col(NY,NX) + RadSWbySoil
!
! AERODYNAMIC RESISTANCE ABOVE SOIL SURFACE INCLUDING
! RESISTANCE IMPOSED BY PLANT CANOPY
!
! FracSoiPAsAir*=air-filled porosity of soil
! DFVR=porosity limitation to diffusion through soil
! POROQ=soil tortuosity
! ResistanceLitRLay=porosity-limited litter blr
! RAGZ=combined soil+litter blr
! RI=Richardsons number
! RIB=isothermal RI
! TKQ,TK1=canopy air,soil temperature
! RAGZ,RAa=soil+litter blr
! ResistBndlSurf_col=isothermal blr at ground surface
!
  FracSoiPAsAir0            = AMAX1(ZERO,FracSoilPoreAsAir_vr(0,NY,NX))
  DFVR                      = FracSoiPAsAir0*POROQ*FracSoiPAsAir0/POROS_vr(0,NY,NX)
  ResistanceLitRLay         = ResistAreodynOverSoil_col(NY,NX)+VapDiffusResistanceLitR(NY,NX)/DFVR
  RI                        = RichardsonNumber(RIB(NY,NX),TKQ_col(NY,NX),TKSoil1_vr(NUM(NY,NX),NY,NX))
  RAGX                      = AMAX1(RAM,0.8_r8*ResistBndlSurf_col(NY,NX),AMIN1(1.2_r8*ResistBndlSurf_col(NY,NX), ResistanceLitRLay/(1.0_r8-10.0_r8*RI)))
  ResistBndlSurf_col(NY,NX) = RAGX
  RAa                       = ResistAreodynOverLitr_col(NY,NX)+ResistBndlSurf_col(NY,NX)
!
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!
! CdSoiEvap,CdSoiHSens=conductance for latent,sensible heat fluxes over soil
! PAREG,PARSG=conductances for latent,sensible heat fluxes
! RZ=minimum surface resistance
! VaporSoi1,VPQ_col=vapor pressure at soil surface, canopy air
! VapXAir2TopLay=negative of evaporation, or water vapor transfer from canopy air to top layer soil/lake
! LatentHeatEvapAir2Grnd=latent heat flux
! XH=rate constant
! VOLW2=soil water volume
! VAP=latent heat of evaporation
! HeatSensVapAir2Soi=convective heat of evaporation flux
!
  CdSoiEvap  = AScaledCdWOverSoil_col(NY,NX)/(RAa+RZ)
  CdSoiHSens = AScaledCdHOverSoil_col(NY,NX)/RAa
  TKX1       = TKSoil1_vr(NUM(NY,NX),NY,NX)
  if(TKX1<0._r8)then
  write(*,*)'TKX1=',TKX1
  call endrun("Negative temperature in "//trim(mod_filename),__LINE__)
  endif
  VaporSoi1=vapsat(TKX1)*EXP(18.0_r8*PSISV1/(RGASC*TKX1))

  !evaporation, no more than what is available, g H2O
  VapXAir2TopLay=AMAX1(CdSoiEvap*(VPQ_col(NY,NX)-VaporSoi1),-AZMAX1(TopLayWatVol*dts_wat))   

  !latent heat > 0 into soil/ground
  LatentHeatEvapAir2Grnd=VapXAir2TopLay*EvapLHTC
  IF(VapXAir2TopLay.LT.0.0_r8)THEN
    !evaporation (<0 into atmosphere)
    HeatSensVapAir2Soi=VapXAir2TopLay*cpw*TKSoil1_vr(NUM(NY,NX),NY,NX)
  ELSE
    !condensation (>0 into ground)
    HeatSensVapAir2Soi=VapXAir2TopLay*cpw*TKQ_col(NY,NX)
  ENDIF
  !take away water from evaporation
  TopLayWatVol=TopLayWatVol+VapXAir2TopLay
!
! SOLVE FOR SOIL SURFACE TEMPERATURE AT WHICH ENERGY
! BALANCE OCCURS, SOLVE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
! HeatSensAir2Grnd,LatentHeatEvapAir2Grnd,Radnet2LitGrnd=sensible,latent heat fluxes, net radiation
! HeatSensVapAir2Soi=convective heat flux from LatentHeatEvapAir2Grnd, > 0 into ground
! HeatFluxAir2Soi=storage heat flux
!
  HeatSensAir2Grnd = CdSoiHSens*(TKQ_col(NY,NX)-TKSoil1_vr(NUM(NY,NX),NY,NX))
  tHeatAir2Grnd    = Radnet2LitGrnd+LatentHeatEvapAir2Grnd+HeatSensAir2Grnd !net energy into soil, subtracting latent heat and sensible heat
  HeatFluxAir2Soi  = tHeatAir2Grnd+HeatSensVapAir2Soi !total heat plus convective heat > 0 to ground

  !write(*,*) "printing heat flux vars: "
  !write(*,*) "CdSoiEvap:", CdSoiEvap
  !write(*,*) "CdSoiHSens:", CdSoiHSens
  !write(*,*) "TKX1:", TKX1
  !write(*,*) "VaporSoi1:", VaporSoi1
  !write(*,*) "VapXAir2TopLay:", VapXAir2TopLay
  !write(*,*) "LatentHeatEvapAir2Grnd:", LatentHeatEvapAir2Grnd
  !write(*,*) "HeatSensVapAir2Soi:", HeatSensVapAir2Soi
  !write(*,*) "TopLayWatVol:", TopLayWatVol
  !write(*,*) "HeatSensAir2Grnd:", HeatSensAir2Grnd
  !write(*,*) "tHeatAir2Grnd:", tHeatAir2Grnd
  !write(*,*) "HeatFluxAir2Soi:", HeatFluxAir2Soi 
  !write(*,*) "AScaledCdWOverSoil_col(NY,NX):", AScaledCdWOverSoil_col(NY,NX)
  !write(*,*) "RAa:", RAa
  !write(*,*) "RZ:", RZ
  !write(*,*) "AScaledCdHOverSoil_col(NY,NX):", AScaledCdHOverSoil_col(NY,NX)
  !write(*,*) "VPQ_col(NY,NX):", VPQ_col(NY,NX)
  !write(*,*) "PSISV1:", PSISV1
  !write(*,*) "RGASC:", RGASC
  !write(*,*) "EvapLHTC:", EvapLHTC
  !write(*,*) "cpw:", cpw
  !write(*,*) "TKQ_col(NY,NX):", TKQ_col(NY,NX)
  !write(*,*) "Radnet2LitGrnd:", Radnet2LitGrnd
  !write(*,*) "dts_wat:", dts_wat
  !write(*,*) "NUM(NY,NX):", NUM(NY,NX)

  end subroutine SoilSRFEnerbyBalance

!------------------------------------------------------------------------------------------

  subroutine ExposedSoilFlux(M,I,J,NY,NX,CumNetWatFlow2LitR,CumNetHeatFlow2LitR,&
    CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil,ResistanceLitRLay,&
    TopLayWatVol,VapXAir2TopLay,HeatFluxAir2Soi,NetWatFlx2SoiMicP)
  implicit none
  integer, intent(in) :: M,I,J   !soil heat-flow iteration id
  integer, intent(in) :: NY,NX
  real(r8), intent(in)  :: CumNetWatFlow2LitR,CumNetHeatFlow2LitR
  real(r8), intent(in) :: CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil
  real(r8) ,intent(inout) :: ResistanceLitRLay
  real(r8), intent(inout) :: TopLayWatVol   !water content in topsoil layer 
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi    !heat flux from air to top soil [MJ]
  real(r8), intent(out) :: NetWatFlx2SoiMicP  !
  real(r8) :: PSISV1
  real(r8) :: LWRadGrnd
  real(r8) :: HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,EvapLitR2Soi1,HeatFluxAir2LitR

! begin_execution
! Watch out for L, is its value defined?
  call SoilSRFEnerbyBalance(M,I,J,NY,NX,PSISV1,LWRadGrnd,ResistanceLitRLay,TopLayWatVol,&
    VapXAir2TopLay,HeatFluxAir2Soi)
!
  call SurfLitREnergyBalance(I,J,M,NY,NX,PSISV1,Prec2LitR2,PrecHeat2LitR2,CumNetWatFlow2LitR,&
    CumNetHeatFlow2LitR,CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil,&
    HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,EvapLitR2Soi1,HeatFluxAir2LitR)
!
  call SumAftEnergyBalance(NY,NX,LWRadGrnd,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,HeatFluxAir2LitR,HeatFluxAir2Soi,NetWatFlx2SoiMicP)

  end subroutine ExposedSoilFlux

!------------------------------------------------------------------------------------------

  subroutine AtmLandSurfExchange(I,J,M,NY,NX,ResistanceLitRLay,TopLayWatVol,LatentHeatAir2Sno,&
    HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,VapXAir2TopLay,HeatFluxAir2Soi1)
  implicit none
  integer, intent(in) :: M           !soil heat-flow iteration id
  integer, intent(in) :: NY,NX,I,J
  real(r8), intent(inout) :: ResistanceLitRLay
  real(r8), intent(inout) :: TopLayWatVol  
  real(r8), intent(out) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvapAir2Snow
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi1   !MJ/d2
  real(r8) :: cumNetHeatFlow2Soil
  real(r8) :: CumWatFlx2SoiMicP,CumWatFlx2SoiMacP
  real(r8) :: CumNetHeatFlow2LitR
  real(r8) :: CumWatXFlx2SoiMicP
  real(r8) :: CumNetWatFlow2LitR
  real(r8) :: HeatNetFlx2Snow
  real(r8) :: NetWatFlx2SoiMicP
  integer  :: L
! begin_execution
  VapXAir2TopLay       = 0._r8
  CumWatXFlx2SoiMicP   = 0._r8
  CumNetWatFlow2LitR  = 0._r8
  CumNetHeatFlow2LitR = 0._r8
  CumWatFlx2SoiMicP    = 0._r8
  CumWatFlx2SoiMacP    = 0._r8
  HeatSensEvapAir2Snow = 0._r8
  HeatNetFlx2Snow      = 0._r8
  LatentHeatAir2Sno    = 0._r8
  HeatSensAir2Snow     = 0._r8
  Radnet2Snow          = 0._r8
  cumNetHeatFlow2Soil   = 0._r8
  HeatFluxAir2Soi1     = 0._r8

  !solve if there is significant snow layer 
  IF(VLSnowHeatCapM_snvr(M,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
!    print*,'SolveSnowpack'
!   VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
    call SolveSnowpack(I,J,M,NY,NX,LatentHeatAir2Sno,Radnet2Snow,HeatSensEvapAir2Snow,HeatSensAir2Snow,&
      HeatNetFlx2Snow,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP,CumNetWatFlow2LitR,&
      CumNetHeatFlow2LitR,cumNetHeatFlow2Soil)
  ENDIF

  !AGGREGATE fluxes from solving the snow model
  WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)   = CumWatFlx2SoiMicP
  WaterFlow2MicptX_3D(3,NUM(NY,NX),NY,NX)  = CumWatXFlx2SoiMicP
  ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX) = CumWatFlx2SoiMacP
  HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)    = cumNetHeatFlow2Soil
  WatFLow2LitR_col(NY,NX)                  = CumNetWatFlow2LitR
  HeatFLoByWat2LitRi_col(NY,NX)            = CumNetHeatFlow2LitR
!
! ENERGY EXCHANGE AT SOIL SURFACE IF EXPOSED UNDER SNOWPACK
! FSNW,FSNX=fractions of snow,snow-free cover
  IF(FracSurfSnoFree_col(NY,NX).GT.0.0_r8 .AND. (SoiBulkDensity_vr(NUM(NY,NX),NY,NX).GT.ZERO.OR. &
    VHeatCapacity1_vr(NUM(NY,NX),NY,NX).GT.VHCPNX(NY,NX)))THEN
    !Ground partically covered by snow, focus on litter-soil interaction     
    call ExposedSoilFlux(M,I,J,NY,NX,CumNetWatFlow2LitR,CumNetHeatFlow2LitR,&
      CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil,ResistanceLitRLay,&
      TopLayWatVol,VapXAir2TopLay,HeatFluxAir2Soi1,NetWatFlx2SoiMicP)
  ELSE
!   ground is fully snow covered, thus no flux from soil & litter
    call SnowCoveredTopSoilFlux(NY,NX,NetWatFlx2SoiMicP)
  ENDIF

!
! AGGREGATE RESIDUE AND SOIL SURFACE FLUXES BENEATH SNOW
! AND ATMOSPHERE
!
! WaterFlow2Micpt_3D,WaterFlow2MicptX_3D=total water flux into soil micropores
! FLWHL=total water flux into soil macropores
! HFLWL=total heat flux into soil
! FLWRL,WaterFlow2MicptX_3D=total water flux into litter
! HFLWRL=total heat flux into litter
! FLWV*=total internal vapor flux in soil
!
  WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)   = WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)+WatNetFlo2TopSoiMicP
  WaterFlow2MicptX_3D(3,NUM(NY,NX),NY,NX)  = WaterFlow2MicptX_3D(3,NUM(NY,NX),NY,NX)+NetWatFlx2SoiMicP
  ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX) = ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)+NetWatFlx2SoiMacP
  HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)    = HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)+cumHeatSensAir2Soil
  !fluxes to litter
  WatFLow2LitR_col(NY,NX)                  = WatFLow2LitR_col(NY,NX)+NetWatFlx2LitR
  HeatFLoByWat2LitRi_col(NY,NX)            = HeatFLoByWat2LitRi_col(NY,NX)+CumHeatSensAir2LitR

  end subroutine AtmLandSurfExchange
!------------------------------------------------------------------------------------------

  subroutine SnowCoveredTopSoilFlux(NY,NX,NetWatFlx2SoiMicP)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8),intent(out):: NetWatFlx2SoiMicP
! begin_execution
  Radnet2LitGrnd         = 0.0_r8           !net radiation into soil
  LatentHeatEvapAir2Grnd = 0.0_r8   !latent heat flux from air and topsoil
  HeatSensVapAir2Soi     = 0.0_r8       !convective heat flux from air and topsoil
  HeatSensAir2Grnd       = 0.0_r8         !sensible heat flux from air to topsoil

  WatNetFlo2TopSoiMicP    = 0.0_r8     !total water flux from air into top soil
  NetWatFlx2SoiMicP       = 0.0_r8        !total water flux from air into soil micropores
  NetWatFlx2SoiMacP       = 0.0_r8        !total water flux from air into macropores
  cumHeatSensAir2Soil     = 0.0_r8      !total water associated heat flux from air into soil
  NetWatFlx2LitR          = 0.0_r8            !total water flux from air to litter
  CumHeatSensAir2LitR     = 0.0_r8       !total water associated heat flux from air to litter
  VapXAir2LitR_col(NY,NX) = 0.0_r8  !evaporative flux from air into litter
  end subroutine SnowCoveredTopSoilFlux  
!------------------------------------------------------------------------------------------

  subroutine SurfLitrSoilWaterExchange(I,J,M,NY,NX,KSatReductByRainKineticEnergy)
  implicit none
  integer, intent(in) :: I,J,M,NY,NX
  real(r8),intent(in) :: KSatReductByRainKineticEnergy

  real(r8) :: THETW1,ThetaWLitR,PSIST1
  real(r8) :: PSIST0,HeatFlxLitR2Soi,FLQZ,DarcyFlxLitR2Soil
  real(r8) :: WatDarcyFloLitR2SoiMicP,FLQ2,CNDR,DarcyCondLitR2Soil
  real(r8) :: HeatFlowLitR2MacP,WatFlowLitR2MacP    
  real(r8) :: CND1
  integer :: K0,K1
  real(r8) :: TKRt,VHCapLitR,ENGYR,VLWatLitR

! begin_execution
! CNDR,HCNDLitr_col=current,saturated litter hydraulic conductivity
! PSISE,PSISM1_vr(0,=air entry,current litter water potential
! VLWatMicP1_vr(0,VWatLitRHoldCapcity=current,maximum litter water volume
! CND1,HydroCond_3D=soil hydraulic conductivity
! KSatReductByRainKineticEnergy=reduction in soil surface Ksat from rainfall energy impact
! K1=soil relative water-filled porosity
! THETWX,POROS=soil water content,porosity
! DarcyCondLitR2Soil=litter-soil hydraulic conductance
! DLYRR,DLYR=litter,soil thicknesses
! WatDarcyFloLitR2SoiMicP,FLQL=litter-soil water flux unltd,ltd by water
! dts_HeatWatTP=time step of flux calculations
! CVRD=fraction of litter cover
! PSISM1_vr(NUM=soil matric water potential
! VLWatMicP1_vr(NUM=soil water volume
! HFLQL=convective heat from litter-soil water flux
! FLWL,HFLWL=micropore water,heat flux
! WatFLow2LitR,HFLWRL=total litter water,heat flux
! FLWRM=litter-soil water flux for solute transfer in TranspNoSalt.f
! CND1,CNDL=hydraulic conductivity of source,destination layer
! HydroCond_3D=lateral(1,2),vertical(3) micropore hydraulic conductivity
!

  !check litter temperature
  VLWatLitR=AMIN1(VLWatMicP1_vr(0,NY,NX)+WatFLow2LitR_col(NY,NX),VLWatMicP1_vr(0,NY,NX))

!  if(I==178 .and. j>=1)then
!  ENGYR     = TKSoil1_vr(0,NY,NX)*VHeatCapacity1_vr(0,NY,NX)
!  VHCapLitR = VHeatCapacity1_vr(0,NY,NX)+cpw*WatFLow2LitR_col(NY,NX)
!  TKRt      = (ENGYR+HeatFLoByWat2LitRi_col(NY,NX))/VHCapLitR
!  write(112,*)I,J,M,'afsnowup1',TKRt-273.15,VHCapLitR
!  endif

  IF(SoiBulkDensity_vr(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    !top layer is soil
    IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS2(NY,NX))THEN
      !surface litter holds water
      ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VLWatLitR)/VLitR_col(NY,NX)
    ELSE
      ThetaWLitR=POROS0(NY,NX)
    ENDIF
    THETW1=AMAX1(THETY_vr(NUM(NY,NX),NY,NX),AMIN1(POROS_vr(NUM(NY,NX),NY,NX) &
      ,safe_adb(VLWatMicP1_vr(NUM(NY,NX),NY,NX),VLSoilMicP_vr(NUM(NY,NX),NY,NX))))
    !K0 litter layer  
    !K1 topsoil layer    
    !DarcyFlxLitR2Soil = water flux from litter layer into the topsoil    
    K0                 = MAX(1,MIN(100,INT(100.0_r8*(AZMAX1(POROS0(NY,NX)-ThetaWLitR))/POROS0(NY,NX))+1))
    K1                 = MAX(1,MIN(100,INT(100.0_r8*(AZMAX1(POROS_vr(NUM(NY,NX),NY,NX)-THETW1))/POROS_vr(NUM(NY,NX),NY,NX))+1))
    CNDR               = HydroCond_3D(3,K0,0,NY,NX)
    CND1               = HydroCond_3D(3,K1,NUM(NY,NX),NY,NX)*KSatReductByRainKineticEnergy
    DarcyCondLitR2Soil = 2.0_r8*CNDR*CND1/(CNDR*DLYR(3,NUM(NY,NX),NY,NX)+CND1*DLYRR_COL(NY,NX))
    PSIST0             = PSISM1_vr(0,NY,NX)+PSIGrav_vr(0,NY,NX)+PSISoilOsmotic_vr(0,NY,NX)
    PSIST1             = PSISM1_vr(NUM(NY,NX),NY,NX)+PSIGrav_vr(NUM(NY,NX),NY,NX)+PSISoilOsmotic_vr(NUM(NY,NX),NY,NX)
    DarcyFlxLitR2Soil  = DarcyCondLitR2Soil*(PSIST0-PSIST1)*AREA(3,NUM(NY,NX),NY,NX)*FracEffAsLitR_col(NY,NX)*dts_HeatWatTP

    IF(DarcyFlxLitR2Soil.GE.0.0_r8)THEN
      !flow from liter to soil
      !more water than saturation
      !note dts_wat:=AMIN1(1.0_r8,20.0_r8*dts_HeatWatTP)
      IF(ThetaWLitR.GT.ThetaSat_vr(0,NY,NX))THEN
        FLQZ=DarcyFlxLitR2Soil+AMIN1((ThetaWLitR-ThetaSat_vr(0,NY,NX))*VLitR_col(NY,NX),&
          AZMAX1((ThetaSat_vr(NUM(NY,NX),NY,NX)-THETW1)*VLSoilMicP_vr(NUM(NY,NX),NY,NX)))*dts_wat
      ELSE
        FLQZ=DarcyFlxLitR2Soil
      ENDIF
      WatDarcyFloLitR2SoiMicP = AZMAX1(AMIN1(FLQZ,VLWatLitR*dts_wat,VLairMicP1_vr(NUM(NY,NX),NY,NX)))
      FLQ2                    = AZMAX1(AMIN1(DarcyFlxLitR2Soil,VLWatLitR*dts_wat,VLairMicP1_vr(NUM(NY,NX),NY,NX)))
    ELSE
      !from soil to liter
      IF(THETW1.GT.ThetaSat_vr(NUM(NY,NX),NY,NX))THEN
        FLQZ=DarcyFlxLitR2Soil+AMAX1((ThetaSat_vr(NUM(NY,NX),NY,NX)-THETW1)*VLSoilMicP_vr(NUM(NY,NX),NY,NX),&
          AZMIN1((ThetaWLitR-ThetaSat_vr(0,NY,NX))*VLitR_col(NY,NX)))*dts_wat
      ELSE
        FLQZ=DarcyFlxLitR2Soil
      ENDIF
      WatDarcyFloLitR2SoiMicP = AZMIN1(AMAX1(FLQZ,-VLWatMicP1_vr(NUM(NY,NX),NY,NX)*dts_wat,-VLairMicP1_vr(0,NY,NX)))
      FLQ2                    = AZMIN1(AMAX1(DarcyFlxLitR2Soil,-VLWatMicP1_vr(NUM(NY,NX),NY,NX)*dts_wat,-VLairMicP1_vr(0,NY,NX)))
    ENDIF

    !top soil layer has negative air volume 
    IF(VLairMicP_vr(NUM(NY,NX),NY,NX).LT.0.0_r8)THEN
      WatDarcyFloLitR2SoiMicP = WatDarcyFloLitR2SoiMicP+AZMIN1(AMAX1(-VLWatMicP1_vr(NUM(NY,NX),NY,NX)*dts_wat,VLairMicP_vr(NUM(NY,NX),NY,NX)))
      FLQ2                    = FLQ2+AZMIN1(AMAX1(-VLWatMicP1_vr(NUM(NY,NX),NY,NX)*dts_wat,VLairMicP_vr(NUM(NY,NX),NY,NX)))
    ENDIF

    IF(WatDarcyFloLitR2SoiMicP.GT.0.0_r8)THEN
      !litter layer to soil
      HeatFlxLitR2Soi=cpw*TKSoil1_vr(0,NY,NX)*WatDarcyFloLitR2SoiMicP
    ELSE
      !soil to litter layer
      HeatFlxLitR2Soi=cpw*TKSoil1_vr(NUM(NY,NX),NY,NX)*WatDarcyFloLitR2SoiMicP
    ENDIF

    WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX) = WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)  = HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)+HeatFlxLitR2Soi
    WatFLow2LitR_col(NY,NX)                = WatFLow2LitR_col(NY,NX)-WatDarcyFloLitR2SoiMicP
    HeatFLoByWat2LitRi_col(NY,NX)          = HeatFLoByWat2LitRi_col(NY,NX)-HeatFlxLitR2Soi
    WatFLo2LitrM(M,NY,NX)                  = WatDarcyFloLitR2SoiMicP
    !soil to litter
  ELSE
    !top layer is water
    WatDarcyFloLitR2SoiMicP                = XVLMobileWatMicP(NY,NX)*dts_wat
    HeatFlxLitR2Soi                        = cpw*TKSoil1_vr(0,NY,NX)*WatDarcyFloLitR2SoiMicP
    WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX) = WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)  = HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)+HeatFlxLitR2Soi
    WatFLow2LitR_col(NY,NX)                = WatFLow2LitR_col(NY,NX)-WatDarcyFloLitR2SoiMicP
    HeatFLoByWat2LitRi_col(NY,NX)          = HeatFLoByWat2LitRi_col(NY,NX)-HeatFlxLitR2Soi
    WatFLo2LitrM(M,NY,NX)                  = WatDarcyFloLitR2SoiMicP
  ENDIF

!     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
!     OF THE LITTER IS EXCEEDED
!
!     VOLPH1=air-filled macroporosity
!     FINHR,HFINHR=water,convective heat from litter to macropores
!     VLWatMicP1_vr(0,VWatLitRHoldCapcity=current,maximum litter water volume
!     WaterFlow2Micpt_3D,HFLWL=micropore water,heat flux
!     WatFLow2LitR,HFLWRL=total litter water,heat flux
!
  IF(VLairMacP1_vr(NUM(NY,NX),NY,NX).GT.0.0_r8 .AND. XVLMobileWatMicP(NY,NX).GT.0.0_r8)THEN
    WatFlowLitR2MacP                         = AMIN1(XVLMobileWatMicP(NY,NX)*dts_wat,VLairMacP1_vr(NUM(NY,NX),NY,NX))
    HeatFlowLitR2MacP                        = WatFlowLitR2MacP*cpw*TKSoil1_vr(0,NY,NX)
    ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX) = ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)+WatFlowLitR2MacP
    HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)    = HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)+HeatFlowLitR2MacP
!    if(HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)>10._r8)then
!      write(*,*)'HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)x=',HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX),HeatFlowLitR2MacP
!    endif
    WatFLow2LitR_col(NY,NX)       = WatFLow2LitR_col(NY,NX)-WatFlowLitR2MacP
    HeatFLoByWat2LitRi_col(NY,NX) = HeatFLoByWat2LitRi_col(NY,NX)-HeatFlowLitR2MacP
  ENDIF

  end subroutine SurfLitrSoilWaterExchange
!------------------------------------------------------------------------------------------

  subroutine InfilSRFRoffPartition(I,J,M,NY,NX)
  !
  !Partition surface water into infiltration and runoff
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NY,NX
  integer :: N1,N2
  real(r8) :: LitrIceHeatFlxFrezPt,TK1X,ENGYR,VLWatMicP1X,VLHeatCapacityX
  real(r8) :: TFREEZ,TFLX,VWatExcess

  real(r8) :: WatExcessDetph
  real(r8) :: Q   !flux of mobile water (>0)
  real(r8) :: HydraulicRadius,CrossSectVelocity

!     begin_execution
!
!     FREEZE-THAW IN RESIDUE SURFACE FROM NET CHANGE IN RESIDUE
!     SURFACE HEAT STORAGE
!
!     TFREEZ=litter freezing temperature
!     PSISM1=litter water potential
!     VLWatMicP1*,VLiceMicP1=litter water,ice volume
!     VLHeatCapacity*=litter volumetric heat capacity
!     TK1*=litter temperature
!     ORGC_vr=litter organic C
!     HFLWRL=total litter conductive, convective heat flux
!     LitrIceHeatFlxFrezPt,TFLX=unltd,ltd latent heat from freeze-thaw
!     LitrIceHeatFlxFrez,LitrIceFlxThaw=litter water,latent heat flux from freeze-thaw
! using Clausius-Clapeyron equation for freezing temeprature calculation
  TFREEZ          = -9.0959E+04_r8/(PSISM1_vr(0,NY,NX)-LtHeatIceMelt)
  VLWatMicP1X     = AZMAX1(VLWatMicP1_vr(0,NY,NX)+WatFLow2LitR_col(NY,NX))
  ENGYR           = VHeatCapacity1_vr(0,NY,NX)*TKSoil1_vr(0,NY,NX)
  VLHeatCapacityX = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP1X+cpi*VLiceMicP1_vr(0,NY,NX)

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    TK1X=(ENGYR+HeatFLoByWat2LitRi_col(NY,NX))/VLHeatCapacityX
  ELSE
    TK1X=TKSoil1_vr(0,NY,NX)
  ENDIF

  !do freeze thaw of litter water
  IF((TK1X.LT.TFREEZ .AND. VLWatMicP1_vr(0,NY,NX).GT.ZERO*VGeomLayer_vr(0,NY,NX)) &
    .OR. (TK1X.GT.TFREEZ .AND. VLiceMicP1_vr(0,NY,NX).GT.ZERO*VGeomLayer_vr(0,NY,NX)))THEN
    LitrIceHeatFlxFrezPt = VHeatCapacity1_vr(0,NY,NX)*(TFREEZ-TK1X) &
      /((1.0_r8+TFREEZ*6.2913E-03_r8)*(1.0_r8-0.10_r8*PSISM1_vr(0,NY,NX)))*dts_wat

    IF(LitrIceHeatFlxFrezPt.LT.0.0_r8)THEN
      !ice thaw, current temperature above freezing TFREEZ<TK1X
      TFLX=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1_vr(0,NY,NX)*dts_wat,LitrIceHeatFlxFrezPt)
    ELSE
      !water freeze (TFLX > 0._r8)
      TFLX=AMIN1(LtHeatIceMelt*VLWatMicP1X*dts_wat,LitrIceHeatFlxFrezPt)
    ENDIF
    LitrIceHeatFlxFrez_col(NY,NX) = TFLX
    LitrIceFlxThaw_col(NY,NX)     = -TFLX/LtHeatIceMelt
  ELSE
    LitrIceFlxThaw_col(NY,NX)     = 0.0_r8
    LitrIceHeatFlxFrez_col(NY,NX) = 0.0_r8
  ENDIF
!
!     THICKNESS OF WATER FILMS IN LITTER AND SOIL SURFACE
!     FROM WATER POTENTIALS FOR GAS EXCHANGE IN TranspNoSalt.F
!
  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    FILM(M,0,NY,NX)=FilmThickness(PSISM1_vr(0,NY,NX), is_top_layer=.true.)
  ELSE
    FILM(M,0,NY,NX)=1.0E-03_r8
  ENDIF
  FILM(M,NUM(NY,NX),NY,NX)=FilmThickness(PSISM1_vr(NUM(NY,NX),NY,NX),is_top_layer=.true.)
!
!     OVERLAND FLOW WHEN WATER STORAGE CAPACITY
!     OF THE SOIL SURFACE PLUS MACROPORES IS EXCEEDED
!
  N1=NX;N2=NY
!
!     SURFACE WATER FLUX
!
!     N2,N1=NY,NX of source grid cell
!     XVOLT,XVOLW=excess water+ice,water in source grid cell
!     VWatStoreCapSurf_col=ground surface water retention capacity
!     VWatExcess=ponded water volume above surface retention capacity
!     D,R,S,V=depth,perimeter,slope,velocity of runoff
!     DIST=distance between source,destination
!     ZM=surface roughness height for runoff
!     Q=runoff from Mannings equation
!     QRM,QRV=runoff,velocity for erosion, solute transfer
! there is mobile water
  IF(XVLMobileWaterLitR_col(N2,N1).GT.VWatStoreCapSurf_col(N2,N1))THEN
    VWatExcess                    = XVLMobileWaterLitR_col(N2,N1)-VWatStoreCapSurf_col(N2,N1)
    WatExcessDetph                = VWatExcess/AREA(3,0,N2,N1)
    HydraulicRadius               = WatExcessDetph/2.828_r8
    CrossSectVelocity             = GaucklerManningVelocity(HydraulicRadius,SLOPE(0,N2,N1))/SoiSurfRoughness(N2,N1)  ![1/s]
    Q                             = CrossSectVelocity*WatExcessDetph*AREA(3,NUM(N2,N1),N2,N1)*3.6E+03_r8*dts_HeatWatTP  ![kg/h/d2]
    VLWatMicP1X                   = AZMAX1(VLWatMicP1_vr(0,N2,N1)+LitrIceFlxThaw_col(N2,N1))
    RunoffVelocity(M,N2,N1)       = CrossSectVelocity
    WatFlux4ErosionM_2DH(M,N2,N1) = AMIN1(Q,VWatExcess*dts_wat,VLWatMicP1X*dts_wat) &
      *XVLMobileWatMicP(N2,N1)/XVLMobileWaterLitR_col(N2,N1)
  ELSE
    RunoffVelocity(M,N2,N1)       = 0.0_r8
    WatFlux4ErosionM_2DH(M,N2,N1) = 0.0_r8
  ENDIF
  end subroutine InfilSRFRoffPartition
!------------------------------------------------------------------------------------------

  subroutine AccumWaterVaporHeatFluxes(I,J,M,NY,NX,LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,&
    Radnet2Snow,VapXAir2TopLay)
  implicit none
  integer , intent(in) :: I,J  
  integer , intent(in) :: M,NY,NX
  real(r8), intent(in) :: LatentHeatAir2Sno,Radnet2Snow
  real(r8), intent(in) :: HeatSensAir2Snow,HeatSensEvapAir2Snow
  real(r8), intent(in) :: VapXAir2TopLay
  real(r8) :: WatFLo2LitrPrev,VLWatLitR, dHeatLitR
! begin_execution
! HOURLY-ACCUMULATED WATER, VAPOR AND HEAT FLUXES THROUGH
! SURFACE RESIDUE AND SOIL SURFACE
!
! TLitrIceFlxThaw,TLitrIceHeatFlxFrez=litter water,heat fluxes from freeze-thaw
! FLW,FLWH,HFLW=soil micropore,macropore,heat fluxes
! FLWR,HFLWR=litter water,heat fluxes
! FLSW,WatConvSno2MacP=water from snowpack to soil micropores,macropores
! HeatConvSno2Soi=convective heat from snowpack to soil
! WatConvSno2LitR=water flux from snowpack to litter
! HeatConvSno2LitR=convective heat flux from snowpack to litter
! HEATI,HeatEvapAir2Surf_col,HeatSensAir2Surf_col,HEATG=net radiation,latent,sensible,storage heat
! VapXAir2GSurf=total evaporation
! FLWM,WaterFlow2MacPM_3D=water flux into soil micropore,macropore for use in TranspNoSalt.f
! VLWatMicPX1=VLWatMicP1 accounting for wetting front
!
  VLWatLitR=VLWatMicP1_vr(0,NY,NX)+WatFLow2LitR_col(NY,NX)+LitrIceFlxThaw_col(NY,NX)
  if(VLWatLitR<0._r8)then
    dHeatLitR                     = safe_adb(VLWatLitR-tiny_wat,WatFLow2LitR_col(NY,NX))*HeatFLoByWat2LitRi_col(NY,NX)
    WatFLow2LitR_col(NY,NX)       = WatFLow2LitR_col(NY,NX)-VLWatLitR+tiny_wat
    HeatFLoByWat2LitRi_col(NY,NX) = HeatFLoByWat2LitRi_col(NY,NX)-dHeatLitR
  endif
  WatFLo2LitrPrev                         = WatFLo2LitR_col(NY,NX)
  TLitrIceFlxThaw_col(NY,NX)              = TLitrIceFlxThaw_col(NY,NX)+LitrIceFlxThaw_col(NY,NX)
  TLitrIceHeatFlxFrez_col(NY,NX)          = TLitrIceHeatFlxFrez_col(NY,NX)+LitrIceHeatFlxFrez_col(NY,NX)
  WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX) = WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX)+WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)
  WaterFlowSoiMicPX(3,NUM(NY,NX),NY,NX)   = WaterFlowSoiMicPX(3,NUM(NY,NX),NY,NX)+WaterFlow2MicptX_3D(3,NUM(NY,NX),NY,NX)
  WaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)    = WaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)+ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)
  HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX)    = HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX)+HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)
  WatFLo2LitR_col(NY,NX)                  = WatFLo2LitR_col(NY,NX)+WatFLow2LitR_col(NY,NX)
  HeatFLo2LitrByWat(NY,NX)                = HeatFLo2LitrByWat(NY,NX)+HeatFLoByWat2LitRi_col(NY,NX)
  HeatByRad2Surf_col(NY,NX)               = HeatByRad2Surf_col(NY,NX)+Radnet2LitGrnd+Radnet2Snow
  HeatSensAir2Surf_col(NY,NX)             = HeatSensAir2Surf_col(NY,NX)+HeatSensAir2Grnd+HeatSensAir2Snow
  HeatEvapAir2Surf_col(NY,NX)             = HeatEvapAir2Surf_col(NY,NX)+LatentHeatEvapAir2Grnd+LatentHeatAir2Sno
  HeatSensVapAir2Surf_col(NY,NX)          = HeatSensVapAir2Surf_col(NY,NX)+HeatSensVapAir2Soi+HeatSensEvapAir2Snow
  HeatNet2Surf_col(NY,NX)                 = HeatNet2Surf_col(NY,NX)+Radnet2LitGrnd+Radnet2Snow &
    +HeatSensAir2Grnd+HeatSensAir2Snow+LatentHeatEvapAir2Grnd+LatentHeatAir2Sno &
    +HeatSensVapAir2Soi+HeatSensEvapAir2Snow


  !EVAPG=negative evaporation from ground/top soil layer
  !EVAPR=evaporation from litter layer   
  !EVAPSN=evaporation from snow, sublimation+evaporation
  VapXAir2GSurf_col(NY,NX)                 = VapXAir2GSurf_col(NY,NX)+VapXAir2TopLay+VapXAir2LitR_col(NY,NX)+VapXAir2Sno(NY,NX)   !>0 into ground
  WaterFlow2MicPM_3D(M,3,NUM(NY,NX),NY,NX) = WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)
  WaterFlow2MacPM_3D(M,3,NUM(NY,NX),NY,NX) = ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)
  TEvapXAir2Toplay_col(NY,NX)              = TEvapXAir2Toplay_col(NY,NX)+VapXAir2TopLay
  TEvapXAir2LitR_col(NY,NX)                = TevapXAir2LitR_col(NY,NX)+VapXAir2LitR_col(NY,NX)
  TEvapXAir2Snow_col(NY,NX)                = TEvapXAir2Snow_col(NY,NX)+VapXAir2Sno(NY,NX)
  end subroutine AccumWaterVaporHeatFluxes

!------------------------------------------------------------------------------------------

  subroutine InitSurfModel(I,J,M,NY,NX,ResistanceLitRLay,KSatReductByRainKineticEnergy)
  implicit none
  integer, intent(in) :: M   !soil heat-flow iteration id
  integer, intent(in) :: NY,NX,I,J
  real(r8),intent(in) :: ResistanceLitRLay
  real(r8),intent(out):: KSatReductByRainKineticEnergy
  integer  :: L  
  real(r8) :: scalar,THETWT,HFLQR1,FLQRS
  real(r8) :: FLQRH,VOLAT0,ENGYD
  real(r8) :: ENGYB,RAS,TScal4Aquadifsvity,THETWA
  real(r8) :: MobileH2OSurf   !mm H2O exceed surface water holding capacity

! begin_execution
! INITIALIZE NET SURFACE FLUX ACCUMULATORS
!
! TQS1,TQW1,TQI1=net water and snowpack snow,water,ice runoff
! cumHeatFlx2LitRByRunoff_col,THQS1=net convective heat from surface water and snow runoff
! FracAsExposedSoil_col,FracEffAsLitR_col=fractions of soil,litter cover including free water+ice
! ResistBndlSurf_col= boundary layer resistance at soil surface
! CondGasXSnowM_col=boundary layer conductance above soil surface
!
  call ZeroSnowFlux(NY,NX)

  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    FracAsExposedSoil_col(NY,NX)=AZMAX1(FracSurfBareSoil_col(NY,NX)-AMIN1(1.0_r8,&
      AZMAX1(XVLMobileWaterLitR_col(NY,NX)/VLWatheldCapSurf_col(NY,NX))))
  ELSE
    FracAsExposedSoil_col(NY,NX)=1.0_r8
  ENDIF
  FracEffAsLitR_col(NY,NX)   = 1.0_r8-FracAsExposedSoil_col(NY,NX)
  ResistBndlSurf_col(NY,NX)  = 1.0_r8/(FracAsExposedSoil_col(NY,NX)/ResistAreodynOverLitr_col(NY,NX)+FracEffAsLitR_col(NY,NX)/ResistanceLitRLay)
  RAS                        = SnowBNDResistance(NY,NX)
  CondGasXSnowM_col(M,NY,NX) = AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP/(ResistBndlSurf_col(NY,NX)+RAS)  !m^2 h/(h/m) = m3
!  write(141,*)I+J/24.,M,CondGasXSnowM_col(M,NY,NX),ResistBndlSurf_col(NY,NX),RAS,ResistAreodynOverLitr_col(NY,NX),FracEffAsLitR_col(NY,NX),ResistanceLitRLay
!
! REDISTRIBUTE INCOMING PRECIPITATION
! BETWEEN RESIDUE AND SOIL SURFACE
!
! BKDS=bulk density
! FLQRS,FLQRH=water flux from soil micropores,macropores to litter
! Prec2SoiMicP1,Prec2SoiMacP1,Prec2LitR1=rain+irrigation to micropores,macropores,litter
! VOLP1,VOLPH1=air-filled microporosity,macroporosity
! HFLQR1=convective heat flux from soil to litter
! PrecAir2LitR,PrecHeatAir2LitR=total water flux, convective heat flux to litter
! PrecNet2SoiMicP,PrecNet2SoiMacP=total water flux to soil micropores, macropores
! PrecHeat2SoiNet=total convective heat flux to soil micropores, macropores
! XNPR=time step for litter water,heat flux calculations
!
  IF(SoiBulkDensity_vr(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    !soil bulk density significant
    FLQRS            = AZMAX1(Prec2SoiMicP1(NY,NX)-VLairMicP1_vr(NUM(NY,NX),NY,NX))
    FLQRH            = AZMAX1(Prec2SoiMacP1(NY,NX)-VLairMacP1_vr(NUM(NY,NX),NY,NX))
    HFLQR1           = cpw*TairK_col(NY,NX)*(FLQRS+FLQRH)
    PrecAir2LitR     = Prec2LitR1(NY,NX)+FLQRS+FLQRH
    PrecHeatAir2LitR = PrecHeat2LitR1(NY,NX)+HFLQR1
    PrecNet2SoiMicP  = Prec2SoiMicP1(NY,NX)-FLQRS
    PrecNet2SoiMacP  = Prec2SoiMacP1(NY,NX)-FLQRH
    PrecHeat2SoiNet  = PrecHeat2SoiMicP1(NY,NX)-HFLQR1
  ELSE
    PrecAir2LitR     = Prec2LitR1(NY,NX)
    PrecHeatAir2LitR = PrecHeat2LitR1(NY,NX)
    PrecNet2SoiMicP  = Prec2SoiMicP1(NY,NX)
    PrecNet2SoiMacP  = Prec2SoiMacP1(NY,NX)
    PrecHeat2SoiNet  = PrecHeat2SoiMicP1(NY,NX)
  ENDIF
  Prec2LitR2     = PrecAir2LitR*XNPR
  PrecHeat2LitR2 = PrecHeatAir2LitR*XNPR
!
! WATER GAS EXCHANGE COEFFICIENTS IN SURFACE LITTER
!
! VOLA1,VLiceMicP1,VLWatMicP1,VsoiPM=total,ice-,water-,air-filled porosity
! TScal4Aquadifsvity=temperature effect on gas diffusivity
! DiffusivitySolutEff=rate constant for air-water gas exchange
! Z1R,Z2RW,Z2RD,Z3RX=parameters for litter air-water gas transfers
! XNPD=time step for gas transfer calculations, it is tunable parameter
! TORT=tortuosity for aqueous diffusivity
! VOLAT0=ice-excluded porosity in litter

  VOLAT0=VLPoreLitR_col(NY,NX)-VLiceMicP1_vr(0,NY,NX)
  IF(VOLAT0.GT.ZEROS2(NY,NX).AND.VLsoiAirPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    !litter layer is not saturated
    THETWA                         = AZMAX1(AMIN1(1.0_r8,VLWatMicP1_vr(0,NY,NX)/VOLAT0))
    TScal4Aquadifsvity             = TEFAQUDIF(TKSoil1_vr(0,NY,NX))
    scalar                         = TScal4Aquadifsvity*XNPD
    DiffusivitySolutEff(M,0,NY,NX) = fDiffusivitySolutEff(scalar,THETWA,0.0_r8,is_litter=.true.)
  ELSE
    !litter layer saturated
    DiffusivitySolutEff(M,0,NY,NX)=0.0_r8
  ENDIF
! VWatLitRHoldCapcity=surface litter water holding capacity, [m3 d-2]
  IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS(NY,NX))THEN
    !litter is able to hold water
    THETWT=AMIN1(1.0_r8,VLWatMicP_vr(0,NY,NX)/VWatLitRHoldCapcity_col(NY,NX))
  ELSE
    THETWT=1.0
  ENDIF
  !TORT=tortuosity in litter (treated as micropore)
  TortMicPM_vr(M,0,NY,NX)=TortMicporeW(THETWT)
!
! KINETIC ENERGY OF DIRECT RAINFALL AND THROUGHFALL
!
! PRECD,PRECB=direct,indirect precipn+irrign at soil surface
! ENGYD,ENGYB=energy impact of direct,indirect precipn+irrign at soil surface
! VWatStoreCapSurf_col=ground surface water retention capacity
! XVOLW=free surface water
! ZT=canopy height
! EnergyImpact4ErosionM=total energy impact for use in erosion.f
! EnergyImpact4Erosion=cumulative rainfall energy impact on soil surface
! KSatReductByRainKineticEnergy=reduction in soil surface Ksat from rainfall energy impact
! Note: A good reference for the following formula and alternatives
! is "Rainfall intensity-kinetic energy relationships for soil loss prediction",
! Kinnell, 1981

  IF(PrecDirect2Grnd_col(NY,NX).GT.ZERO)THEN
    ENGYD=AZMAX1(8.95_r8+8.44_r8*LOG(PRECM_col(NY,NX)))
  ELSE
    ENGYD=0.0_r8
  ENDIF
  IF(PrecIndirect2Grnd_col(NY,NX).GT.ZERO)THEN
    ENGYB=AZMAX1(15.8_r8*SQRT(AMIN1(2.5_r8,CanopyHeight_col(NY,NX)))-5.87_r8)
  ELSE
    ENGYB=0.0_r8
  ENDIF

  IF(ENGYD+ENGYB.GT.ZERO)THEN
    MobileH2OSurf                  = 1.0E+03_r8*AZMAX1(XVLMobileWaterLitR_col(NY,NX)-VWatStoreCapSurf_col(NY,NX))/AREA(3,NUM(NY,NX),NY,NX)
    EnergyImpact4ErosionM(M,NY,NX) = (ENGYD*PrecDirect2Grnd_col(NY,NX)+ENGYB*PrecIndirect2Grnd_col(NY,NX))*EXP(-2.0_r8*MobileH2OSurf) &
      *FracSurfBareSoil_col(NY,NX)*dts_HeatWatTP
    EnergyImpact4Erosion(NY,NX)=EnergyImpact4Erosion(NY,NX)+EnergyImpact4ErosionM(M,NY,NX)
  ELSE
    EnergyImpact4ErosionM(M,NY,NX)=0.0_r8
  ENDIF
  KSatReductByRainKineticEnergy=EXP(-2.0E-03_r8*(CSILT(NUM(NY,NX),NY,NX)+CCLAY(NUM(NY,NX),NY,NX)) &
    *EnergyImpact4Erosion(NY,NX))

!
!  SNOWPACK FLUX ACCUMULATORS
!
   call InitSnowAccums(I,J,M,NY,NX)

  call PrepIterSnowLayer(I,J,M,NY,NX)

!
  end subroutine InitSurfModel
!------------------------------------------------------------------------------------------

  subroutine SurfaceRunoff(I,J,M,N,NN,N1,N2,M4,M5,RCHQF,XN)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,N,NN
  integer, intent(in) :: N1,N2  !source grid
  integer, intent(in) :: M4,M5  !dest grid
  real(r8),intent(in) :: RCHQF  !water flux scalar
  real(r8),intent(in) :: XN     !flow direction
  real(r8) :: ALT1     !elevation at the source grid
  real(r8) :: ALT2     !elevation in the dest center
  real(r8) :: DPTHW1,DPTHW2
  real(r8) :: VX
  !
  ! SURFACE BOUNDARY WATER FLUX
  !
  ! DPTHW1,DPTHW2=surface water depth of source,destination
  ! ALT1,ALT2=elevation of source,destination
  ! XVOLT=excess surface water+ice
  ! VWatStoreCapSurf_col=ground surface water retention capacity
  ! ExtWaterTable=natural water table depth
  ! QR1,HeatFlx2LitRByRunoff=runoff, convective heat from runoff
  ! QR,HQR=hourly-accumulated runoff, convective heat from runoff
  ! QRM,QRV=runoff,velocity for erosion, solute transfer
  ! XN=direction
  !
  ! RUNOFF
  !
  DPTHW1 = XVLMobileWaterLitR_col(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
  DPTHW2 = VWatStoreCapSurf_col(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)

  ALT1 = Altitude_grid(N2,N1)+DPTHW1
  ALT2 = Altitude_grid(N2,N1)+DPTHW2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)

  !grid elevation is higher than outside the grid, and in grid water layer higher than external water table
  !depth is counting downward
  IF(ALT1.GT.ALT2 .AND. CumDepz2LayerBot_vr(NU(N2,N1)-1,N2,N1)-DPTHW1.LT.ExtWaterTable_col(N2,N1))THEN
    !grid has more water than water table, 
    !out of grid (N2,N1), WatFlux4ErosionM is computed from surface physics model
    WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)   = -XN*WatFlux4ErosionM_2DH(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)  = cpw*TKSoil1_vr(0,N2,N1)*WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    XGridSurfRunoff_2DH(N,NN,M5,M4)       = XGridSurfRunoff_2DH(N,NN,M5,M4)+WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    HeatXGridBySurfRunoff_2DH(N,NN,M5,M4) = HeatXGridBySurfRunoff_2DH(N,NN,M5,M4)+HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
! RUNON
! water table in higher than grid surface (accouting for minimum water )
  ELSEIF(CumDepz2LayerBot_vr(NU(N2,N1)-1,N2,N1)-DPTHW1.GT.ExtWaterTable_col(N2,N1))THEN
    !elevation difference
    VX                                    = AZMIN1((ExtWaterTable_col(N2,N1)-CumDepz2LayerBot_vr(NU(N2,N1)-1,N2,N1)+DPTHW1)*AREA(3,NUM(N2,N1),N2,N1))
    WatFlux4ErosionM_2DH(M,N2,N1)         = VX*dts_wat
    RunoffVelocity(M,N2,N1)               = 0.0_r8
    WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)   = -XN*WatFlux4ErosionM_2DH(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)  = cpw*TKSoil1_vr(0,N2,N1)*WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    XGridSurfRunoff_2DH(N,NN,M5,M4)       = XGridSurfRunoff_2DH(N,NN,M5,M4)+WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    HeatXGridBySurfRunoff_2DH(N,NN,M5,M4) = HeatXGridBySurfRunoff_2DH(N,NN,M5,M4)+HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
  ELSE
    WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)  = 0.0_r8
    HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4) = 0.0_r8
  ENDIF
  QflxSurfRunoffM_2DH(M,N,NN,M5,M4) = WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
  IFLBM(M,N,NN,M5,M4)               = 0

  end subroutine SurfaceRunoff
!------------------------------------------------------------------------------------------

  subroutine PartitionPrecip(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  real(r8) :: SnoFall,HeatByPrec
  real(r8) :: PrecThrufall2SoiMicP,PrecThrufall2SoiMacP,Rain4ToSno
  real(r8) :: PrecThrufall2LitR,PrecThrufall2Soil,PrecHeat2LitR,PrecHeat2Soil
  real(r8) :: PrecThruFall  
  real(r8), parameter :: m2mm=1.e3_r8
!     PRECA=precipitation+irrigation
!     PRECD,PRECB=direct,indirect precipn+irrign at soil surface
!     Prec2Canopy_col=net ice transfer to canopy, updated in hour1


  !convert water flux from m/hour to mm/hour
  PRECM_col(NY,NX)             = m2mm*PrecRainAndIrrig_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  PrecDirect2Grnd_col(NY,NX)   = m2mm*(PrecRainAndIrrig_col(NY,NX)-Prec2Canopy_col(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  PrecIndirect2Grnd_col(NY,NX) = m2mm*(Prec2Canopy_col(NY,NX)-PrecIntceptByCanopy_col(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
!
!     RESIDUE WATER ABSORPTION CAPACITY
!
!     HCNDLitr_col=litter saturated hydraulic conductivity
!     DLYRR=litter depth
!
!  HCNDLitr_col(NY,NX) = SatHydroCondLitR
  DLYRR_COL(NY,NX)    = AMAX1(2.5E-03_r8,DLYR(3,0,NY,NX))
!
!     DISTRIBUTION OF PRECIPITATION AND ITS HEAT AMONG SURFACE
!     RESIDUE, SOIL SURFACE, AND MACROPORES
!
!     PRECA,SnoFalPrec=rainfall+irrigation,snowfall (water equiv)
!     Rain4ToSno=rainfall to snowpack
!     SnoFall=snowfall to snowpack
!     HeatByPrec=convective heat flux to snowpack
!     PrecThruFall=precip to litter+soil surfaces
!     PrecThrufall2Soil,PrecThrufall2LitR=precip to soil,litter surfaces
!     PrecHeat2Soil,PrecHeat2LitR=convective heat flux to soil,litter surfaces
!     PrecThrufall2SoiMicP,PrecThrufall2SoiMacP=precip to soil micropores,macropores
!     TFLWC=canopy intercepted precipitation
!     FSNW=fraction of snow cover
!     PrecIntceptByCanopy_col=precipitation intercepted by plant canopy

  IF(PrecRainAndIrrig_col(NY,NX).GT.0.0_r8 .OR. SnoFalPrec_col(NY,NX).GT.0.0_r8)THEN
  ! there is precipitation
    Rain4ToSno           = (PrecRainAndIrrig_col(NY,NX)-PrecIntceptByCanopy_col(NY,NX))*FracSurfAsSnow_col(NY,NX)
    SnoFall              = SnoFalPrec_col(NY,NX)                                                        !snowfall
    HeatByPrec           = cps*TairK_col(NY,NX)*SnoFall+cpw*TairK_col(NY,NX)*Rain4ToSno                  !incoming heat flux from precipitations to snow-covered surface
    PrecThruFall         = (PrecRainAndIrrig_col(NY,NX)-PrecIntceptByCanopy_col(NY,NX))*FracSurfSnoFree_col(NY,NX)       !incoming precipitation to snow-free surface
    PrecThrufall2LitR    = PrecThruFall*FracSurfByLitR_col(NY,NX)                             !water flux to snow-free coverd by litter
    PrecHeat2LitR        = cpw*TairK_col(NY,NX)*PrecThrufall2LitR                                 !heat flux to snow-free surface covered by litter
    PrecThrufall2Soil    = PrecThruFall*FracSurfBareSoil_col(NY,NX)                          !heat flux to snow-free surface not covered by litter
    PrecHeat2Soil        = cpw*TairK_col(NY,NX)*PrecThrufall2Soil
    PrecThrufall2SoiMicP = PrecThrufall2Soil*SoilFracAsMicP_vr(NUM(NY,NX),NY,NX)          !water flux to micropore
    PrecThrufall2SoiMacP = PrecThrufall2Soil*SoilFracAsMacP1_vr(NUM(NY,NX),NY,NX)         !water flux to macropore
  ELSE
  ! no precipitation
    Rain4ToSno           = -PrecIntceptByCanopy_col(NY,NX)*FracSurfAsSnow_col(NY,NX)                   !
    SnoFall              = 0.0_r8
    HeatByPrec           = cpw*TairK_col(NY,NX)*Rain4ToSno
    PrecThruFall         = -PrecIntceptByCanopy_col(NY,NX)*FracSurfSnoFree_col(NY,NX)
    PrecThrufall2LitR    = PrecThruFall*FracSurfByLitR_col(NY,NX)
    PrecHeat2LitR        = cpw*TairK_col(NY,NX)*PrecThrufall2LitR
    PrecThrufall2Soil    = PrecThruFall*FracSurfBareSoil_col(NY,NX)
    PrecHeat2Soil        = cpw*TairK_col(NY,NX)*PrecThrufall2Soil
    PrecThrufall2SoiMicP = PrecThrufall2Soil*SoilFracAsMicP_vr(NUM(NY,NX),NY,NX)
    PrecThrufall2SoiMacP = PrecThrufall2Soil*SoilFracAsMacP1_vr(NUM(NY,NX),NY,NX)
  ENDIF
  HeatPrec_col(NY,NX)=HeatPrec_col(NY,NX)+(HeatByPrec+PrecHeat2Soil)*dts_HeatWatTP
!
!     PRECIP ON SNOW ARRAYS EXPORTED TO TranspNoSalt.F, TranspSalt.F
!     FOR SOLUTE FLUX CALCULATIONS
!
!     SnoFalPrec,RainFalPrec_col,PrecAtm_col,PRECI=snow,rain,snow+rain,irrigation
!     VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
!     Rain2LitRSurf_col,Irrig2LitRSurf=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!
  IF(SnoFalPrec_col(NY,NX).GT.0.0_r8 .OR. (RainFalPrec_col(NY,NX).GT.0.0_r8 &
    .AND. VLHeatCapSnow_snvr(1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)))THEN
    !there is precipitation, there is significant snow layer
    Rain2LitRSurf_col(NY,NX)=0.0_r8
    Irrig2LitRSurf(NY,NX)=0.0_r8
    Rain2SoilSurf_col(NY,NX)=PrecAtm_col(NY,NX)
    Irrig2SoilSurf(NY,NX)=IrrigSurface_col(NY,NX)
  ELSEIF((PrecAtm_col(NY,NX).GT.0.0.OR.IrrigSurface_col(NY,NX).GT.0.0_r8) &
    .AND.VLHeatCapSnow_snvr(1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX))THEN
    !there is insignificant snow layer
    Rain2LitRSurf_col(NY,NX)=PrecThrufall2LitR*PrecAtm_col(NY,NX)/(PrecAtm_col(NY,NX)+IrrigSurface_col(NY,NX))
    Irrig2LitRSurf(NY,NX)=PrecThrufall2LitR*IrrigSurface_col(NY,NX)/(PrecAtm_col(NY,NX)+IrrigSurface_col(NY,NX))
    Rain2SoilSurf_col(NY,NX)=PrecAtm_col(NY,NX)-Rain2LitRSurf_col(NY,NX)
    Irrig2SoilSurf(NY,NX)=IrrigSurface_col(NY,NX)-Irrig2LitRSurf(NY,NX)
  ELSE
    !no precipitation
    Rain2LitRSurf_col(NY,NX)=0.0_r8
    Irrig2LitRSurf(NY,NX)=0.0_r8
    Rain2SoilSurf_col(NY,NX)=0.0_r8
    Irrig2SoilSurf(NY,NX)=0.0_r8
  ENDIF
!
!     GATHER PRECIPITATION AND MELTWATER FLUXES AND THEIR HEATS
!     AMONG ATMOSPHERE, SNOWPACK, RESIDUE AND SOIL SURFACES
!     INTO LOCAL ARRAYS FOR USE IN MASS AND ENERGY EXCHANGE
!     ALGORITHMS
!
!     dts_HeatWatTP=internal time step for fluxes through soil profile
!
!     FLW0S,Ice2Snowt,Rain2Snowt=snow,ice,water input to snowpack
!     HeatFall2Snowt=convective heat flux to snowpack
!     Prec2SoiMicP1,Prec2SoiMacP1,Prec2LitR1=rain+irrigation to micropores,macropores,litter
!     PrecHeat2SoiMicP1,HWFLY1=convective heat flux to soil,litter surfaces
!
  SnowFallt(NY,NX)      = SnoFall*dts_HeatWatTP
  Ice2Snowt(NY,NX)      = 0.0_r8
  Rain2Snowt(NY,NX)     = Rain4ToSno*dts_HeatWatTP
  HeatFall2Snowt(NY,NX) = HeatByPrec*dts_HeatWatTP

  Prec2SoiMicP1(NY,NX) = PrecThrufall2SoiMicP*dts_HeatWatTP
  Prec2SoiMacP1(NY,NX) = PrecThrufall2SoiMacP*dts_HeatWatTP

  Prec2LitR1(NY,NX)        = PrecThrufall2LitR*dts_HeatWatTP
  PrecHeat2SoiMicP1(NY,NX) = PrecHeat2Soil*dts_HeatWatTP
  PrecHeat2LitR1(NY,NX)    = PrecHeat2LitR*dts_HeatWatTP

  end subroutine PartitionPrecip
!------------------------------------------------------------------------------------------  

  subroutine UpdateSurfaceAtM(I,J,M,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: M,I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer  :: NY,NX

! begin_execution
  D9795: DO NX=NHW,NHE
    D9790: DO NY=NVN,NVS

      !update snow
      call UpdateSnowAtM(I,J,M,NY,NX)
      
      call UpdateLitRAftRunoff(I,J,M,NY,NX)
      
    ENDDO D9790
  ENDDO D9795

  end subroutine UpdateSurfaceAtM
!------------------------------------------------------------------------------------------

  subroutine SumAftEnergyBalance(NY,NX,LWRadGrnd,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,HeatFluxAir2LitR,HeatFluxAir2Soi,NetWatFlx2SoiMicP)
  implicit none
  integer, intent(in)  :: NY,NX
  real(r8), intent(in) :: LWRadGrnd            !>0 into atmosphere, long wave radiation from ground surface
  real(r8), intent(in) :: VapXAir2TopLay
  real(r8), intent(in) :: HeatSensLitR2Soi1
  real(r8), intent(in) :: HeatSensVapLitR2Soi1
  real(r8), intent(in) :: EvapLitR2Soi1
  real(r8), intent(in) :: HeatFluxAir2LitR
  real(r8), intent(in) :: HeatFluxAir2Soi
  real(r8), intent(out):: NetWatFlx2SoiMicP
  real(r8) :: FLWVLS  
! begin_execution
!
! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
! FOR LATER UPDATES TO STATE VARIABLES
!
! WatNetFlo2TopSoiMicP,NetWatFlx2SoiMacP=water flux from atm to soil micropores,macropores
! cumHeatSensAir2Soil=convective heat flux from atm to soil
! NetWatFlx2LitR=water flux from atm to litter
! CumHeatSensAir2LitR=convective heat flux from atm to litter
! FLWVLS=water flux within soil accounting for wetting front
!
  WatNetFlo2TopSoiMicP = PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1
  NetWatFlx2SoiMicP    = PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1
  NetWatFlx2SoiMacP    = PrecNet2SoiMacP
  cumHeatSensAir2Soil  = PrecHeat2SoiNet+HeatFluxAir2Soi+HeatSensVapLitR2Soi1+HeatSensLitR2Soi1
  !total water and heat fluxes to litter
  NetWatFlx2LitR       = PrecAir2LitR+VapXAir2LitR_col(NY,NX)-EvapLitR2Soi1
  CumHeatSensAir2LitR  = PrecHeatAir2LitR+HeatFluxAir2LitR-HeatSensVapLitR2Soi1-HeatSensLitR2Soi1
  FLWVLS               = (VLWatMicP1_vr(NUM(NY,NX),NY,NX)-VLWatMicPX1_vr(NUM(NY,NX),NY,NX))*dts_HeatWatTP
!
! GENERATE NEW SNOWPACK
!
! XFLWS,WatXfer2SnoLay,IceXfer2SnoLay=hourly snow,water,ice transfer
! SnowFallt,Rain2Snowt,Ice2Snowt=snow,water,ice input to snowpack
! HeatXfer2SnoLay=hourly convective heat flux from snow,water,ice transfer
! HeatFall2Snowt=convective heat flux from snow,water,ice to snowpack
!
  IF(VLHeatCapSnow_snvr(1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) .AND. SnowFallt(NY,NX).GT.ZEROS(NY,NX))THEN
    SnoXfer2SnoLay_snvr(1,NY,NX)  = SnoXfer2SnoLay_snvr(1,NY,NX)+SnowFallt(NY,NX)
    WatXfer2SnoLay_snvr(1,NY,NX)  = WatXfer2SnoLay_snvr(1,NY,NX)+Rain2Snowt(NY,NX)
    IceXfer2SnoLay_snvr(1,NY,NX)  = IceXfer2SnoLay_snvr(1,NY,NX)+Ice2Snowt(NY,NX)
    HeatXfer2SnoLay_snvr(1,NY,NX) = HeatXfer2SnoLay_snvr(1,NY,NX)+HeatFall2Snowt(NY,NX)
  ENDIF
  !LWRadBySurf_col=longwave emission from litter and surface soil into atmosphere
  LWRadBySurf_col(NY,NX)=LWRadBySurf_col(NY,NX)+LWRadGrnd
  end subroutine SumAftEnergyBalance
!------------------------------------------------------------------------------------------
  subroutine RunSurfacePhysModel(I,J,M,NHE,NHW,NVS,NVN,ResistanceLitRLay,KSatReductByRainKineticEnergy,&
    TopLayWatVol,HeatFluxAir2Soi,Qinfl2MicP,Hinfl2Soil)
  !
  !run surface energy/water model for iteration M  
  implicit none
  integer, intent(in) :: I,J !day, hour
  integer, intent(in) :: M   !soil heat-water iteration id
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8), dimension(:,:),intent(inout) :: ResistanceLitRLay(JY,JX)
  REAL(R8), dimension(:,:),INTENT(OUT) :: KSatReductByRainKineticEnergy
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol(JY,JX)
  real(r8), dimension(:,:),intent(out) :: HeatFluxAir2Soi(JY,JX)
  real(r8), dimension(:,:),optional,intent(out) :: Qinfl2MicP(JY,JX)
  real(r8), dimension(:,:),optional,intent(out) :: Hinfl2Soil(JY,JX)    !heat flow into soil [MJ]
  real(r8) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvapAir2Snow,VapXAir2TopLay
  integer :: N1,N2,NX,NY,L

  WatFlx2LitRByRunoff_2DH(:,:,:,:)=0._r8
  HeatFlx2LitRByRunoff_2DH(:,:,:,:)=0._r8

  D9895: DO  NX=NHW,NHE
    D9890: DO  NY=NVN,NVS

      call SurfaceEnergyModel(I,J,M,NX,NY,ResistanceLitRLay(NY,NX),KSatReductByRainKineticEnergy(NY,NX),&
        HeatFluxAir2Soi(NY,NX),LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,&
        TopLayWatVol(NY,NX),VapXAir2TopLay)

    ! CAPILLARY EXCHANGE OF WATER BETWEEN SOIL SURFACE AND RESIDUE
      call SurfLitrSoilWaterExchange(I,J,M,NY,NX,KSatReductByRainKineticEnergy(NY,NX))

      call InfilSRFRoffPartition(I,J,M,NY,NX)
    !
      if(.not.ATS_cpl_mode)call LateralGridsHdryoExch(I,J,M,NY,NX,NHE,NHW,NVS,NVN)

      if(snowRedist_model)call SnowRedistribution(M,NY,NX,NHE,NHW,NVS,NVN)

    ! In ATS coupled mode we do not run the full Redist so we put the snow
    ! models here instead
      if (ATS_cpl_mode) then
        call SnowMassUpdate(I,J,NY,NX)
        call SnowpackLayering(I,J,NY,NX)
      end if

      call AccumWaterVaporHeatFluxes(I,J,M,NY,NX,LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,&
        Radnet2Snow,VapXAir2TopLay)

      call UpdateLitRB4RunoffM(I,J,M,NY,NX)
      if(present(Qinfl2MicP))Qinfl2MicP(NY,NX)=WaterFlow2Micpt_3D(3,NUM(NY,NX),NY,NX)
      if(present(Hinfl2Soil))Hinfl2Soil(NY,NX)=HeatFlow2Soili_3D(3,NUM(NY,NX),NY,NX)

    ENDDO D9890
  ENDDO D9895

  end subroutine RunSurfacePhysModel

!------------------------------------------------------------------------------------------
  subroutine SurfaceEnergyModel(I,J,M,NX,NY,ResistanceLitRLay,KSatReductByRainKineticEnergy,&
    HeatFluxAir2Soi1,LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,&
    TopLayWatVol,VapXAir2TopLay)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M     !soil heat-flow iteration
  integer, intent(in) :: NX,NY
  real(r8),intent(inout) :: ResistanceLitRLay
  REAL(R8), INTENT(OUT) :: KSatReductByRainKineticEnergy,HeatFluxAir2Soi1
  real(r8), intent(out) :: Radnet2Snow,LatentHeatAir2Sno,HeatSensAir2Snow,HeatSensEvapAir2Snow
  real(r8), intent(inout) :: TopLayWatVol
  real(r8), intent(out) :: VapXAir2TopLay
  integer :: N1,N2,L

  !ResistanceLitRLay is input
  call InitSurfModel(I,J,M,NY,NX,ResistanceLitRLay,KSatReductByRainKineticEnergy)

! updates ResistanceLitRLay
  call AtmLandSurfExchange(I,J,M,NY,NX,ResistanceLitRLay,TopLayWatVol,LatentHeatAir2Sno,&
    HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,VapXAir2TopLay,HeatFluxAir2Soi1)

  !update snow pack before doing snow redistribution to avoid negative mass values  
  call UpdateSnowPack1(I,J,M,NY,NX)

  end subroutine SurfaceEnergyModel

!------------------------------------------------------------------------------------------
  pure function GaucklerManningVelocity(HydraulicRadius, slope)result(CrossSectVelocity)

  !https://en.wikipedia.org/wiki/Manning_formula
  !V=k/n*Rh^(2/3)*S^(1/2)
  implicit none
  real(r8), intent(in) :: HydraulicRadius
  real(r8), intent(in) :: slope   !stream slope, or hydraulic gradient, or channel bed slope
  
  real(r8) :: CrossSectVelocity    !m/s

  CrossSectVelocity=HydraulicRadius**0.67_r8*SQRT(SLOPE)

  end function GaucklerManningVelocity

!------------------------------------------------------------------------------------------

  subroutine writeSurfDiagnosis(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,j
  integer, intent(in) :: NY,NX

  write(106,*)I+J/24.,TEvapXAir2Toplay_col(NY,NX)*1000._r8/AREA(3,NU(NY,NX),NY,NX), &
    TEvapXAir2LitR_col(NY,NX)*1000._r8/AREA(3,NU(NY,NX),NY,NX), & 
    TEvapXAir2Snow_col(NY,NX)*1000._r8/AREA(3,NU(NY,NX),NY,NX) 

  end subroutine writeSurfDiagnosis

!------------------------------------------------------------------------------------------

  subroutine AggregateSurfRunoffFlux(I,J,M,NHW,NHE,NVN,NVS)
  implicit none
  
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS

  integer :: N1,N2   !source grid, with which the lateral exchange is computed 
  integer :: N,NN,N4,N5,N4B,N5B
  integer :: NY,NX
  real(r8) :: VLWatLitR,VLicelitR
  real(r8) :: Heatflxlitr
!     begin_execution

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      cumWatFlx2LitRByRunoff_col(NY,NX)      = 0.0_r8
      cumHeatFlx2LitRByRunoff_col(NY,NX)     = 0.0_r8

  VLWatLitR  = AZMAX1(VLWatMicP_vr(0,NY,NX)+WatFLo2LitR_col(NY,NX)+TLitrIceFlxThaw_col(NY,NX)+TXGridSurfRunoff_2DH(NY,NX))
  VLicelitR  = AZMAX1(VLiceMicP_vr(0,NY,NX)-TLitrIceFlxThaw_col(NY,NX)/DENSICE)

      N1=NX;N2=NY  
      DO  N=1,2
        IF(N.EQ.iEastWestDirection)THEN
          !exchange in the x direction, west-east
          N4  = NX+1   !east
          N5  = NY
          N4B = NX-1  !west
          N5B = NY
        ELSEIF(N.EQ.iNorthSouthDirection)THEN
          !exchange in the y direction, north-south
          N4  = NX
          N5  = NY+1    !south
          N4B = NX
          N5B = NY-1   !north
        ENDIF

        DO  NN=1,2        
          !coming in
          cumWatFlx2LitRByRunoff_col(N2,N1)  = cumWatFlx2LitRByRunoff_col(N2,N1)+WatFlx2LitRByRunoff_2DH(N,NN,N2,N1)
          cumHeatFlx2LitRByRunoff_col(N2,N1) = cumHeatFlx2LitRByRunoff_col(N2,N1)+HeatFlx2LitRByRunoff_2DH(N,NN,N2,N1)
          
          VLWatLitR=VLWatMicP1_vr(0,NY,NX)+cumWatFlx2LitRByRunoff_col(N2,N1)
          if(VLWatLitR<0._r8)then
            VLWatLitR                            = VLWatLitR-tiny_wat
            Heatflxlitr                          = safe_adb(VLWatLitR,WatFlx2LitRByRunoff_2DH(N,NN,N2,N1))*HeatFlx2LitRByRunoff_2DH(N,NN,N2,N1)
            WatFlx2LitRByRunoff_2DH(N,NN,N2,N1)  = WatFlx2LitRByRunoff_2DH(N,NN,N2,N1)-VLWatLitR
            HeatFlx2LitRByRunoff_2DH(N,NN,N2,N1) = HeatFlx2LitRByRunoff_2DH(N,NN,N2,N1)-Heatflxlitr
            cumWatFlx2LitRByRunoff_col(N2,N1)    = cumWatFlx2LitRByRunoff_col(N2,N1)-VLWatLitR
            cumHeatFlx2LitRByRunoff_col(N2,N1)   = cumHeatFlx2LitRByRunoff_col(N2,N1)-Heatflxlitr
          endif

          !going out
          IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
            !there is outflow in east or south
            cumWatFlx2LitRByRunoff_col(N2,N1)  = cumWatFlx2LitRByRunoff_col(N2,N1)-WatFlx2LitRByRunoff_2DH(N,NN,N5,N4)
            cumHeatFlx2LitRByRunoff_col(N2,N1) = cumHeatFlx2LitRByRunoff_col(N2,N1)-HeatFlx2LitRByRunoff_2DH(N,NN,N5,N4)

            VLWatLitR=VLWatMicP1_vr(0,NY,NX)+cumWatFlx2LitRByRunoff_col(N2,N1)
            if(VLWatLitR<0._r8)then
              VLWatLitR                            = VLWatLitR-tiny_wat            
              Heatflxlitr                          = safe_adb(VLWatLitR,WatFlx2LitRByRunoff_2DH(N,NN,N5,N4))*HeatFlx2LitRByRunoff_2DH(N,NN,N5,N4)
              WatFlx2LitRByRunoff_2DH(N,NN,N5,N4)  = WatFlx2LitRByRunoff_2DH(N,NN,N5,N4)+VLWatLitR
              HeatFlx2LitRByRunoff_2DH(N,NN,N5,N4) = HeatFlx2LitRByRunoff_2DH(N,NN,N5,N4)+Heatflxlitr
              cumWatFlx2LitRByRunoff_col(N2,N1)    = cumWatFlx2LitRByRunoff_col(N2,N1)-VLWatLitR
              cumHeatFlx2LitRByRunoff_col(N2,N1)   = cumHeatFlx2LitRByRunoff_col(N2,N1)-Heatflxlitr
              if(NN.EQ.2)then
                XGridSurfRunoff_2DH(N,2,N5,N4)       = XGridSurfRunoff_2DH(N,2,N5,N4)+VLWatLitR
                HeatXGridBySurfRunoff_2DH(N,2,N5,N4) = HeatXGridBySurfRunoff_2DH(N,2,N5,N4)+Heatflxlitr
              endif
            endif
          ENDIF

          !going out
          IF(N4B.GT.0 .AND. N5B.GT.0 .AND. NN.EQ.1)THEN
            !there is outflow in west and north
            cumWatFlx2LitRByRunoff_col(N2,N1)  = cumWatFlx2LitRByRunoff_col(N2,N1)-WatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B)
            cumHeatFlx2LitRByRunoff_col(N2,N1) = cumHeatFlx2LitRByRunoff_col(N2,N1)-HeatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B)

            VLWatLitR=VLWatMicP1_vr(0,NY,NX)+cumWatFlx2LitRByRunoff_col(N2,N1)
            if(VLWatLitR<0._r8)then
              VLWatLitR                              = VLWatLitR-tiny_wat
              Heatflxlitr                            = safe_adb(VLWatLitR,WatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B))*HeatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B)
              WatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B)  = WatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B)+VLWatLitR
              HeatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B) = HeatFlx2LitRByRunoff_2DH(N,NN,N5B,N4B)+Heatflxlitr
              cumWatFlx2LitRByRunoff_col(N2,N1)      = cumWatFlx2LitRByRunoff_col(N2,N1)-VLWatLitR
              cumHeatFlx2LitRByRunoff_col(N2,N1)     = cumHeatFlx2LitRByRunoff_col(N2,N1)-Heatflxlitr
              XGridSurfRunoff_2DH(N,1,N5B,N4B)       = XGridSurfRunoff_2DH(N,1,N5B,N4B)+VLWatLitR
              HeatXGridBySurfRunoff_2DH(N,1,N5B,N4B) = HeatXGridBySurfRunoff_2DH(N,1,N5B,N4B)+Heatflxlitr

            endif

          ENDIF

        ENDDO
      ENDDO

      TXGridSurfRunoff_2DH(N2,N1)       = TXGridSurfRunoff_2DH(N2,N1)+cumWatFlx2LitRByRunoff_col(N2,N1)
      THeatXGridBySurfRunoff_2DH(N2,N1) = THeatXGridBySurfRunoff_2DH(N2,N1)+cumHeatFlx2LitRByRunoff_col(N2,N1)

      if(M.EQ.NPH)then
        VLWatMicP1_vr(0,NY,NX)   = VLWatMicP1_vr(0,NY,NX) +cumWatFlx2LitRByRunoff_col(NY,NX)
        VLairMicP1_vr(0,NY,NX)   = AZMAX1(VLPoreLitR_col(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))
        VLWatMicPM_vr(M,0,NY,NX) = VLWatMicP1_vr(0,NY,NX)
        VLsoiAirPM(M,0,NY,NX)    = VLairMicP1_vr(0,NY,NX)
      ELSE

        VLWatLitR  = VLWatMicP_vr(0,NY,NX)+WatFLo2LitR_col(NY,NX)+TLitrIceFlxThaw_col(NY,NX)+TXGridSurfRunoff_2DH(NY,NX)
        VLicelitR  = VLiceMicP_vr(0,NY,NX)-TLitrIceFlxThaw_col(NY,NX)/DENSICE

        if(abs(VLWatLitR-VLWatMicP1_vr(0,NY,NX)-cumWatFlx2LitRByRunoff_col(N2,N1))>tiny_wat)then
          write(*,*)(I*1000+J)*10+M,VLWatLitR,VLWatMicP1_vr(0,NY,NX)+cumWatFlx2LitRByRunoff_col(N2,N1),VLWatMicP_vr(0,NY,NX),VLWatMicP1_vr(0,NY,NX)
          write(*,*)WatFLo2LitR_col(NY,NX),TLitrIceFlxThaw_col(NY,NX),TXGridSurfRunoff_2DH(NY,NX)
          call endrun(trim(mod_filename)//'at line',__LINE__)
        endif

      endif

    ENDDO
  ENDDO    

  end subroutine AggregateSurfRunoffFlux  
end module SurfPhysMod
