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

  !surface model
  public :: StageSurfacePhysModel
  public :: RunSurfacePhysModel
  public :: SurfaceEnergyModel
  public :: writeSurfDiagnosis
  !
  public :: SurfaceRunoff
  public :: UpdateSurfaceAtM

! SatHydroCondLitR=saturated hydraulic conductivity of surface litter
! FEnergyImpact4Erosion=rate constant for restoring surface Ksat
! LitRSurfResistance=minimum boundary layer resistances of litter (h m-1)
! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
! RACX,LitRSurfResistance=minimum boundary layer resistances of canopy,litter (h m-1)

  real(r8), parameter :: FEnergyImpact4Erosion=1.0E-03_r8
  real(r8), parameter :: SatHydroCondLitR=25.0_r8
  real(r8), parameter :: LitRSurfResistance=0.0139_r8
  real(r8), parameter :: SoilEmisivity=0.97_r8        !soil emissivity
  real(r8), parameter :: SnowEmisivity=0.97_r8        !snowpack emissivity
  real(r8), parameter :: SurfLitREmisivity=0.97_r8    !surfce litter emissivity
  real(r8), parameter :: RACX=0.0139_r8    !total canopy boundary later resistance h/m  

  real(r8) :: PrecHeat2LitR2,Prec2LitR2  
  real(r8) :: HeatSensVapAir2Soi
  real(r8) :: HeatSensAir2Grnd,Radnet2LitGrnd
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
    TEvapXAir2Toplay_col(NY,NX)=0._r8
    TEvapXAir2LitR_col(NY,NX) = 0._r8
    TEvapXAir2Snow_col(NY,NX) = 0._r8

    !make a local copy of the upper boundary index
!
!     ADJUST SURFACE ELEVATION USED IN RUNOFF FOR FREEZE-THAW, EROSION
!     AND SOC
!
!     Altitude_grid,ALT=current,initial elevation of ground surface
!     CumDepz2LayerBot_vr(NUM(NY,NX)-1),=depth of ground surface
!     EnergyImpact4Erosion=cumulative rainfall energy impact on soil surface
!
      Altitude_grid(NY,NX)=ALT(NY,NX)-CumDepz2LayerBot_vr(NUM(NY,NX)-1,NY,NX)
      EnergyImpact4Erosion(NY,NX)=EnergyImpact4Erosion(NY,NX)*(1.0_r8-FEnergyImpact4Erosion)

      call CopySnowStates(I,J,NY,NX)
      
      call CopySurfaceVars(I,J,NY,NX)
      
      call PartionSurfaceFraction(NY,NX)

      call PartitionPrecip(I,J,NY,NX)
      
      call SurfaceRadiation(I,J,NY,NX)

      call SurfaceResistances(I,J,NY,NX,ResistanceLitRLay)

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
  LitrIceFlxThaw(NY,NX)     = 0.0_r8
  LitrIceHeatFlxFrez(NY,NX) = 0.0_r8

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
  LWRadBySurf_col(NY,NX)     = 0.0_r8
  VHeatCapacity1_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
  VLPoreLitR(NY,NX)          = VLMicP_vr(0,NY,NX)
  VLWatMicP1_vr(0,NY,NX)     = AZMAX1(VLWatMicP_vr(0,NY,NX))
  VLiceMicP1_vr(0,NY,NX)     = AZMAX1(VLiceMicP_vr(0,NY,NX))
  VLairMicP1_vr(0,NY,NX)     = AZMAX1(VLPoreLitR(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))
  VLWatMicPM_vr(1,0,NY,NX)   = VLWatMicP1_vr(0,NY,NX)
  VLsoiAirPM(1,0,NY,NX)      = VLairMicP1_vr(0,NY,NX)
  TVWatIceLitR               = VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX)
  XVLMobileWaterLitR_col(NY,NX)  = AZMAX1(TVWatIceLitR-VWatLitRHoldCapcity_col(NY,NX))
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
    FracSoiPAsWat_vr(0,NY,NX)=AZMAX1t(VLWatMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoiPAsIce_vr(0,NY,NX)=AZMAX1t(VLiceMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    
    FracSoiPAsAir_vr(0,NY,NX)=AZMAX1t(VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))* &
      AZMAX1t((1.0_r8-XVLMobileWaterLitR_col(NY,NX)/MaxVLWatByLitR_col(NY,NX)))
  ELSE
    FracSoiPAsWat_vr(0,NY,NX)=0.0_r8
    FracSoiPAsIce_vr(0,NY,NX)=0.0_r8
    FracSoiPAsAir_vr(0,NY,NX)=1.0
  ENDIF
  
  THETPM(1,0,NY,NX)  = FracSoiPAsAir_vr(0,NY,NX)
  PSISM1_vr(0,NY,NX) = PSISoilMatricP_vr(0,NY,NX)
  TKSoi1_vr(0,NY,NX)    = TKS_vr(0,NY,NX)

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

  FracSurfAsSnow(NY,NX)=AMIN1(1.0_r8,SQRT((SnowDepth_col(NY,NX)/MinSnowDepth)))
  FracSurfSnoFree(NY,NX)=1.0_r8-FracSurfAsSnow(NY,NX)
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

  real(r8), parameter :: SensHeatCondctance=1.25E-03_r8

  !TLEX=total latent heat flux x boundary layer resistance, [MJ m-1]
  !TSHX=total sensible heat flux x boundary layer resistance, [MJ m-1]
  !VPQ_col=vapor pressure in canopy air, 
  !TKQ=temperature in canopy air, Kelvin

  !write(*,*) "For VPQ_col and TKQ calc:"
  !write(*,*) "VPA(NY,NX) = ", VPA(NY,NX)
  !write(*,*) "TLEX(NY,NX) = ", TLEX(NY,NX)
  !write(*,*) "EvapLHTC = ", EvapLHTC
  !write(*,*) "NY = ", NY
  !write(*,*) "NX = ", NX
  !write(*,*) "NUM(NY,NX) = ", NUM(NY,NX)
  !write(*,*) "AREA(3,NUM(NY,NX),NY,NX) = ", AREA(3,NUM(NY,NX),NY,NX)
  !write(*,*) "TairK_col(NY,NX) = ", TairK_col(NY,NX)
  !write(*,*) "TSHX(NY,NX) = ", TSHX(NY,NX)
  !write(*,*) "SensHeatCondctance = ", SensHeatCondctance

  VPQ_col(NY,NX)=VPA(NY,NX)-TLEX(NY,NX)/(EvapLHTC*AREA(3,NUM(NY,NX),NY,NX))
  TKQ(NY,NX)=TairK_col(NY,NX)-TSHX(NY,NX)/(SensHeatCondctance*AREA(3,NUM(NY,NX),NY,NX))

  end subroutine SetCanopyProperty
!------------------------------------------------------------------------------------------

  subroutine SurfaceRadiation(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: THRYX,RADGX
!
!     INITIALIZE PARAMETERS, FLUXES FOR ENERGY EXCHANGE
!     AT SNOW, RESIDUE AND SOIL SURFACES
!
!     RADGX=shortwave radiation at ground surface
!     RadSWonSno,RadSWonSoi,RadSWonLitR= shortwave radn at snowpack,soil,litter
!     FracSWRad2Grnd_col=fraction of shortwave radiation at ground surface
!     FSNW,FSNX=fractions of snow,snow-free cover
!     BARE,CVRD=fractions of soil,litter cover
!     XNPS=internal time step for fluxes through snowpack
!     THRYX=longwave radiation at ground surface
!     LWRad2Snow,LWRad2Grnd,LWRad2LitR=longwave radn incident at snowpack,soil,litter
!     LWEmscefSnow_col,LWEmscefSoil_col,LWEmscefLitR_col=longwave radn emitted by snowpack,soil,litter
!     SnowEmisivity,SoilEmisivity,SurfLitREmisivity=emissivity of snowpack,soil,litter surfaces
!     THS=sky longwave radiation
!     LWRadCanGPrev=longwave radiation emitted by canopy

  RADGX=RadSWGrnd_col(NY,NX)*dts_HeatWatTP
  RadSWonSno(NY,NX)=RADGX*FracSurfAsSnow(NY,NX)*XNPS
  RadSWonSoi(NY,NX)=RADGX*FracSurfSnoFree(NY,NX)*FracSurfBareSoil_col(NY,NX)      
  RadSWonLitR(NY,NX)=RADGX*FracSurfSnoFree(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR

  THRYX=(LWRadSky(NY,NX)*FracSWRad2Grnd_col(NY,NX)+LWRadCanGPrev(NY,NX))*dts_HeatWatTP
!  if(I<=1 .or. I>=365)then
!  print*,'THRYX',LWRadSky(NY,NX),FracSWRad2Grnd_col(NY,NX),LWRadCanGPrev(NY,NX)
!  print*,THRYX,FracSurfAsSnow(NY,NX),XNPS
!  endif
  LWRad2Snow(NY,NX)=THRYX*FracSurfAsSnow(NY,NX)*XNPS
  LWRad2Grnd(NY,NX)=THRYX*FracSurfSnoFree(NY,NX)*FracSurfBareSoil_col(NY,NX)
  LWRad2LitR(NY,NX)=THRYX*FracSurfSnoFree(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR
  ! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
  !stefboltz_const is stefan-boltzman constant converted into MJ m-2 K-4/per hour
  LWEmscefSnow_col(NY,NX)=SnowEmisivity*stefboltz_const*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX)*dts_HeatWatTP
  LWEmscefSoil_col(NY,NX)=SoilEmisivity*stefboltz_const*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree(NY,NX)*&
    FracSurfBareSoil_col(NY,NX)*dts_HeatWatTP
  LWEmscefLitR_col(NY,NX)=SurfLitREmisivity*stefboltz_const*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree(NY,NX)*&
    FracSurfByLitR_col(NY,NX)*dts_HeatWatTP
!  print*,'frac',FracSurfByLitR_col(NY,NX),FracSurfSnoFree(NY,NX)
!

  end subroutine SurfaceRadiation
!------------------------------------------------------------------------------------------
  subroutine SurfaceResistances(I,J,NY,NX,ResistanceLitRLay)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8), dimension(:,:), intent(out):: ResistanceLitRLay(JY,JX)
  real(r8) :: FracSoiPAsAir0,DFVR,TFACR  
  real(r8) :: PAREX,PARSX,RAS
  real(r8) :: ALFZ,WindSpeedGrnd
!
!     BOUNDARY LAYER CONDUCTANCES FOR EXPORT TO TranspNoSalt.F
!

!     AERODYNAMIC RESISTANCE OF CANOPY TO SNOW/RESIDUE/SOIL
!     SURFACE ENERGY EXCHANGE WITH ATMOSPHERE
!
!     ALFZ=parameter for canopy effect on windspeed
!     FracSWRad2Grnd_col=fraction of shortwave radiation at ground surface
!     AbvCanopyBndlResist_col,RAC=isothermal blr above canopy, canopy blr
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
    BndlResistCanopy_col(NY,NX)=AMIN1(RACX,AZMAX1(CanopyHeight_col(NY,NX)*EXP(ALFZ) &
      /(ALFZ/AbvCanopyBndlResist_col(NY,NX))*AZMAX1(EXP(-ALFZ*SoilSurfRoughnesst0_col(NY,NX)/CanopyHeight_col(NY,NX)) &
      -EXP(-ALFZ*(ZERO4PlantDisplace_col(NY,NX)+RoughHeight_col(NY,NX))/CanopyHeight_col(NY,NX)))))
    WindSpeedGrnd=WindSpeedAtm_col(NY,NX)*EXP(-ALFZ)
  ELSE
    BndlResistCanopy_col(NY,NX) = 0.0_r8
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
!     ResistanceLitRLay=porosity-limited litter blr
!     PAREX,PARSX=conductances for latent,sensible heat fluxes
!     PAREW,PARSW=conductances for snowpack latent,sensible heatfluxes
!     PAREG,PARSG=conductances for soil latent,sensible heat fluxes
!     PARER,PARSR=conductances for litter latent,sensible heat fluxes
!     XNPR=internal time step for fluxes through litter

! Note: VaporDiffusivityLitR_col is computed in Hour1 which is not used in the
!       ATS coupler. However the calculation is simple so I'm just reproducing it
!       here. This depends on the temperature of the liter (TKS) which is the same
!       as TKSoi1 at the surface
  TFACR                           = TEFGASDIF(TKS_vr(0,NY,NX))
  VaporDiffusivityLitR_col(NY,NX) = TFACR*7.70E-02_r8
  VapDiffusResistanceLitR(NY,NX)  = DLYRR_COL(NY,NX)/VaporDiffusivityLitR_col(NY,NX)

  RAG(NY,NX)               = BndlResistCanopy_col(NY,NX)+AbvCanopyBndlResist_col(NY,NX)
  RAGW(NY,NX)              = RAG(NY,NX)
  RAGR(NY,NX)              = RAG(NY,NX)+LitRSurfResistance
  RARG(NY,NX)              = RAGR(NY,NX)
  FracSoiPAsAir0           = AMAX1(ZERO2,FracSoiPAsAir_vr(0,NY,NX))
  DFVR                     = FracSoiPAsAir0*POROQ*FracSoiPAsAir0/POROS_vr(0,NY,NX)
  ResistanceLitRLay(NY,NX) = RAG(NY,NX)+VapDiffusResistanceLitR(NY,NX)/DFVR
  PAREX                    = AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP               !conductance for latent heat flux
  PARSX                    = 1.25E-03_r8*AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP   !conductance for sensible heat flux
  
  PAREW(NY,NX)=PAREX*FracSurfAsSnow(NY,NX)*XNPS
  PARSW(NY,NX)=PARSX*FracSurfAsSnow(NY,NX)*XNPS
  PAREG(NY,NX)=PAREX*FracSurfSnoFree(NY,NX)
  PARER(NY,NX)=PAREX*FracSurfSnoFree(NY,NX)*XNPR*FracSurfByLitR_col(NY,NX)
  PARSG(NY,NX)=PARSX*FracSurfSnoFree(NY,NX)
  PARSR(NY,NX)=PARSX*FracSurfSnoFree(NY,NX)*XNPR*FracSurfByLitR_col(NY,NX)

!     PARR=boundary layer conductance above litter,soil surfaces
!
  RAS=SnowBNDResistance(NY,NX)

  PARR(NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP/(RAGR(NY,NX)+RAS)   !this includes snow layer resistance
  end subroutine SurfaceResistances

  !------------------------------------------------------------------------------------------

  subroutine SoilSRFEnerbyBalance(M,I,J,NY,NX,PSISV1,LWRadGrnd,ResistanceLitRLay,TopLayWatVol,&
    VapXAir2TopLay,HeatFluxAir2Soi)
  implicit none
  integer, intent(in) :: M,I,J,NY,NX
  real(r8), intent(out):: PSISV1,LWRadGrnd
  real(r8),dimension(:,:), intent(inout) :: ResistanceLitRLay(JY,JX)
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol(JY,JX)
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi
  real(r8) :: RAa,VaporSoi1,CdSoiEvap,CdSoiHSens,RAGX
  real(r8) :: tRadIncid    !total incoming radiation to soil or liter surface, short + long
  real(r8) :: RI,WPLX,WPX,FracSoiPAsAir0
  real(r8) :: VLWatGrnd,VLIceGrnd,DFVR
  real(r8) :: TKX1,THETA1S
  real(r8) :: tHeatAir2Grnd,AlbedoGrnd
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
! RadSWonSoi,LWRad2Grnd,Radnet2LitGrnd=incoming shortwave,longwave,net radiation
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
  tRadIncid=(1.0_r8-AlbedoGrnd)*RadSWonSoi(NY,NX)+LWRad2Grnd(NY,NX)
  LWRadGrnd=LWEmscefSoil_col(NY,NX)*TKSoi1_vr(NUM(NY,NX),NY,NX)**4._r8
  Radnet2LitGrnd=tRadIncid-LWRadGrnd
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
! RAGS=isothermal blr at ground surface
!
  FracSoiPAsAir0           = AMAX1(ZERO,FracSoiPAsAir_vr(0,NY,NX))
  DFVR                     = FracSoiPAsAir0*POROQ*FracSoiPAsAir0/POROS_vr(0,NY,NX)
  ResistanceLitRLay(NY,NX) = RAG(NY,NX)+VapDiffusResistanceLitR(NY,NX)/DFVR
  RI                       = RichardsonNumber(RIB(NY,NX),TKQ(NY,NX),TKSoi1_vr(NUM(NY,NX),NY,NX))
  RAGX                     = AMAX1(RAM,0.8_r8*RAGS(NY,NX),AMIN1(1.2_r8*RAGS(NY,NX),&
    ResistanceLitRLay(NY,NX)/(1.0_r8-10.0_r8*RI)))
  RAGS(NY,NX)              = RAGX
  RAa                      = RAGR(NY,NX)+RAGS(NY,NX)
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
  CdSoiEvap  = PAREG(NY,NX)/(RAa+RZ)
  CdSoiHSens = PARSG(NY,NX)/RAa
  TKX1       = TKSoi1_vr(NUM(NY,NX),NY,NX)
  if(TKX1<0._r8)then
  write(*,*)'TKX1=',TKX1
  call endrun("Negative temperature in "//trim(mod_filename),__LINE__)
  endif
  VaporSoi1=vapsat(TKX1)*EXP(18.0_r8*PSISV1/(RGASC*TKX1))

  !evaporation, no more than what is available, g H2O
  VapXAir2TopLay=AMAX1(CdSoiEvap*(VPQ_col(NY,NX)-VaporSoi1),-AZMAX1(TopLayWatVol(NY,NX)*dts_wat))   

  !latent heat > 0 into soil/ground
  LatentHeatEvapAir2Grnd=VapXAir2TopLay*EvapLHTC
  IF(VapXAir2TopLay.LT.0.0_r8)THEN
    !evaporation (<0 into atmosphere)
    HeatSensVapAir2Soi=VapXAir2TopLay*cpw*TKSoi1_vr(NUM(NY,NX),NY,NX)
  ELSE
    !condensation (>0 into ground)
    HeatSensVapAir2Soi=VapXAir2TopLay*cpw*TKQ(NY,NX)
  ENDIF
  !take away water from evaporation
  TopLayWatVol(NY,NX)=TopLayWatVol(NY,NX)+VapXAir2TopLay
!
! SOLVE FOR SOIL SURFACE TEMPERATURE AT WHICH ENERGY
! BALANCE OCCURS, SOLVE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
! HeatSensAir2Grnd,LatentHeatEvapAir2Grnd,Radnet2LitGrnd=sensible,latent heat fluxes, net radiation
! HeatSensVapAir2Soi=convective heat flux from LatentHeatEvapAir2Grnd, > 0 into ground
! HeatFluxAir2Soi=storage heat flux
!
  HeatSensAir2Grnd=CdSoiHSens*(TKQ(NY,NX)-TKSoi1_vr(NUM(NY,NX),NY,NX))
  !net energy into soil, subtracting latent heat and sensible heat
  tHeatAir2Grnd=Radnet2LitGrnd+LatentHeatEvapAir2Grnd+HeatSensAir2Grnd
  !total heat plus convective heat > 0 to ground
  HeatFluxAir2Soi=tHeatAir2Grnd+HeatSensVapAir2Soi

  !write(*,*) "printing heat flux vars: "
  !write(*,*) "CdSoiEvap:", CdSoiEvap
  !write(*,*) "CdSoiHSens:", CdSoiHSens
  !write(*,*) "TKX1:", TKX1
  !write(*,*) "VaporSoi1:", VaporSoi1
  !write(*,*) "VapXAir2TopLay:", VapXAir2TopLay
  !write(*,*) "LatentHeatEvapAir2Grnd:", LatentHeatEvapAir2Grnd
  !write(*,*) "HeatSensVapAir2Soi:", HeatSensVapAir2Soi
  !write(*,*) "TopLayWatVol(NY,NX):", TopLayWatVol(NY,NX)
  !write(*,*) "HeatSensAir2Grnd:", HeatSensAir2Grnd
  !write(*,*) "tHeatAir2Grnd:", tHeatAir2Grnd
  !write(*,*) "HeatFluxAir2Soi:", HeatFluxAir2Soi 
  !write(*,*) "PAREG(NY,NX):", PAREG(NY,NX)
  !write(*,*) "RAa:", RAa
  !write(*,*) "RZ:", RZ
  !write(*,*) "PARSG(NY,NX):", PARSG(NY,NX)
  !write(*,*) "VPQ_col(NY,NX):", VPQ_col(NY,NX)
  !write(*,*) "PSISV1:", PSISV1
  !write(*,*) "RGASC:", RGASC
  !write(*,*) "EvapLHTC:", EvapLHTC
  !write(*,*) "cpw:", cpw
  !write(*,*) "TKQ(NY,NX):", TKQ(NY,NX)
  !write(*,*) "Radnet2LitGrnd:", Radnet2LitGrnd
  !write(*,*) "dts_wat:", dts_wat
  !write(*,*) "NUM(NY,NX):", NUM(NY,NX)

  end subroutine SoilSRFEnerbyBalance

!------------------------------------------------------------------------------------------

  subroutine ExposedSoilFlux(M,I,J,NY,NX,ResistanceLitRLay,TopLayWatVol,VapXAir2TopLay,HeatFluxAir2Soi,&
    NetWatFlx2SoiMicP)
  implicit none
  integer, intent(in) :: M,I,J   !soil heat-flow iteration id
  integer, intent(in) :: NY,NX
  real(r8),dimension(:,:),intent(inout) :: ResistanceLitRLay(JY,JX)
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol(JY,JX)  
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi
  real(r8), intent(out) :: NetWatFlx2SoiMicP
  real(r8) :: PSISV1
  real(r8) :: LWRadGrnd
  real(r8) :: HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR

! begin_execution
! Watch out for L, is its value defined?
  call SoilSRFEnerbyBalance(M,I,J,NY,NX,PSISV1,LWRadGrnd,ResistanceLitRLay,TopLayWatVol,&
    VapXAir2TopLay,HeatFluxAir2Soi)
!
  call SRFLitterEnergyBalance(M,NY,NX,PSISV1,Prec2LitR2,PrecHeat2LitR2,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR)
!
  call SumAftEnergyBalance(NY,NX,LWRadGrnd,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR,HeatFluxAir2Soi,NetWatFlx2SoiMicP)

  end subroutine ExposedSoilFlux

!------------------------------------------------------------------------------------------

  subroutine AtmLandSurfExchange(I,J,M,NY,NX,ResistanceLitRLay,TopLayWatVol,LatentHeatAir2Sno,&
    HeatSensEvap,HeatSensAir2Snow,Radnet2Snow,VapXAir2TopLay,HeatFluxAir2Soi1)
  implicit none
  integer, intent(in) :: M           !soil heat-flow iteration id
  integer, intent(in) :: NY,NX,I,J
  real(r8), dimension(:,:), intent(inout) :: ResistanceLitRLay(JY,JX)
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol(JY,JX)  
  real(r8), intent(out) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvap
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi1   !MJ/d2
  real(r8) :: cumHeatFlowSno2Soi
  real(r8) :: CumWatFlx2SoiMicP,CumWatFlx2SoiMacP
  real(r8) :: CumHeatFlow2LitR
  real(r8) :: CumWatXFlx2SoiMicP
  real(r8) :: CumWatFlow2LitR
  real(r8) :: HeatNetFlx2Snow
  real(r8) :: NetWatFlx2SoiMicP
  integer  :: L
! begin_execution
  VapXAir2TopLay     = 0._r8
  CumWatXFlx2SoiMicP = 0._r8
  CumWatFLow2LitR    = 0._r8
  CumHeatFlow2LitR   = 0._r8
  CumWatFlx2SoiMicP  = 0._r8
  CumWatFlx2SoiMacP  = 0._r8
  HeatSensEvap       = 0._r8
  HeatNetFlx2Snow    = 0._r8
  LatentHeatAir2Sno  = 0._r8
  HeatSensAir2Snow   = 0._r8
  Radnet2Snow        = 0._r8
  cumHeatFlowSno2Soi = 0._r8
  HeatFluxAir2Soi1   = 0._r8

  !solve if there is significant snow layer 
  IF(VLSnowHeatCapM_snvr(M,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
!    print*,'SolveSnowpack'
!   VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
    call SolveSnowpack(I,J,M,NY,NX,LatentHeatAir2Sno,Radnet2Snow,HeatSensEvap,HeatSensAir2Snow,&
      HeatNetFlx2Snow,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP,CumWatFlow2LitR,&
      CumHeatFlow2LitR,cumHeatFlowSno2Soi)
  ENDIF

!
! ENERGY EXCHANGE AT SOIL SURFACE IF EXPOSED UNDER SNOWPACK
! FSNW,FSNX=fractions of snow,snow-free cover
  IF(FracSurfSnoFree(NY,NX).GT.0.0_r8 .AND. (SoiBulkDensity_vr(NUM(NY,NX),NY,NX).GT.ZERO.OR. &
    VHeatCapacity1_vr(NUM(NY,NX),NY,NX).GT.VHCPNX(NY,NX)))THEN
    !Ground partically covered by snow
     
    call ExposedSoilFlux(M,I,J,NY,NX,ResistanceLitRLay,TopLayWatVol,VapXAir2TopLay,HeatFluxAir2Soi1,&
      NetWatFlx2SoiMicP)

  ELSE
!   ground is fully snow covered, thus no flux from soil & litter
    call SnowCoveredTopSoilFlux(NY,NX,NetWatFlx2SoiMicP)
  ENDIF

!
! AGGREGATE RESIDUE AND SOIL SURFACE FLUXES BENEATH SNOW
! AND ATMOSPHERE
!
! WatXChange2WatTable,WatXChange2WatTableX=total water flux into soil micropores
! FLWHL=total water flux into soil macropores
! HFLWL=total heat flux into soil
! FLWRL,WatXChange2WatTableX=total water flux into litter
! HFLWRL=total heat flux into litter
! FLWV*=total internal vapor flux in soil
!
  WatXChange2WatTable(3,NUM(NY,NX),NY,NX)  = CumWatFlx2SoiMicP+WatNetFlo2TopSoiMicP
  WatXChange2WatTableX(3,NUM(NY,NX),NY,NX) = CumWatXFlx2SoiMicP+NetWatFlx2SoiMicP
  ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX) = CumWatFlx2SoiMacP+NetWatFlx2SoiMacP
  HeatFlow2Soili(3,NUM(NY,NX),NY,NX)       = cumHeatFlowSno2Soi+cumHeatSensAir2Soil
!  if(NUM(NY,NX)==1.and.HeatFlow2Soili(3,NUM(NY,NX),NY,NX)>10._r8)then
!    write(*,*)'atlmdn',cumHeatFlowSno2Soi,cumHeatSensAir2Soil
!  endif
  WatFLow2LitR_col(NY,NX)   = CumWatFlow2LitR+NetWatFlx2LitR
  HeatFLoByWat2LitRi_col(NY,NX) = CumHeatFlow2LitR+CumHeatSensAir2LitR

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

  WatNetFlo2TopSoiMicP = 0.0_r8     !total water flux from air into top soil
  NetWatFlx2SoiMicP    = 0.0_r8        !total water flux from air into soil micropores
  NetWatFlx2SoiMacP    = 0.0_r8        !total water flux from air into macropores
  cumHeatSensAir2Soil  = 0.0_r8      !total water associated heat flux from air into soil
  NetWatFlx2LitR       = 0.0_r8            !total water flux from air to litter
  CumHeatSensAir2LitR  = 0.0_r8       !total water associated heat flux from air to litter
  VapXAir2LitR(NY,NX)  = 0.0_r8  !evaporative flux from air into litter
  end subroutine SnowCoveredTopSoilFlux  
!------------------------------------------------------------------------------------------

  subroutine SurfLitrSoilWaterExchange(I,J,M,NY,NX,KSatReductByRainKineticEnergy)
  implicit none
  integer, intent(in) :: I,J,M,NY,NX
  real(r8),intent(in) :: KSatReductByRainKineticEnergy

  real(r8) :: THETW1,ThetaWLitR,PSIST1
  real(r8) :: PSIST0,HeatFlxLitR2Soi,FLQZ,DarcyFlxLitR2Soil
  real(r8) :: WatDarcyFloLitR2SoiMicP,FLQ2,CNDR,AVCNDR
  real(r8) :: HeatFlowLitR2MacP,WatFlowLitR2MacP    
  real(r8) :: CND1
  integer :: K0,K1

! begin_execution
! CNDR,HCNDR=current,saturated litter hydraulic conductivity
! PSISE,PSISM1_vr(0,=air entry,current litter water potential
! VLWatMicP1_vr(0,VWatLitRHoldCapcity=current,maximum litter water volume
! CND1,HydroCond_3D=soil hydraulic conductivity
! KSatReductByRainKineticEnergy=reduction in soil surface Ksat from rainfall energy impact
! K1=soil relative water-filled porosity
! THETWX,POROS=soil water content,porosity
! AVCNDR=litter-soil hydraulic conductance
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
  IF(SoiBulkDensity_vr(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    !top layer is soil
    IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS2(NY,NX))THEN
      !surface litter holds water
      ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VLWatMicP1_vr(0,NY,NX))/VLitR_col(NY,NX)
    ELSE
      ThetaWLitR=POROS0(NY,NX)
    ENDIF
    THETW1=AMAX1(THETY_vr(NUM(NY,NX),NY,NX),AMIN1(POROS_vr(NUM(NY,NX),NY,NX) &
      ,safe_adb(VLWatMicP1_vr(NUM(NY,NX),NY,NX),VLSoilMicP_vr(NUM(NY,NX),NY,NX))))
    !K0 litter layer  
    !K1 topsoil layer    
    !DarcyFlxLitR2Soil = water flux from litter layer into the topsoil    
    K0                = MAX(1,MIN(100,INT(100.0*(AZMAX1(POROS0(NY,NX)-ThetaWLitR))/POROS0(NY,NX))+1))
    K1                = MAX(1,MIN(100,INT(100.0*(AZMAX1(POROS_vr(NUM(NY,NX),NY,NX)-THETW1))/POROS_vr(NUM(NY,NX),NY,NX))+1))
    CNDR              = HydroCond_3D(3,K0,0,NY,NX)
    CND1              = HydroCond_3D(3,K1,NUM(NY,NX),NY,NX)*KSatReductByRainKineticEnergy
    AVCNDR            = 2.0_r8*CNDR*CND1/(CNDR*DLYR(3,NUM(NY,NX),NY,NX)+CND1*DLYRR_COL(NY,NX))
    PSIST0            = PSISM1_vr(0,NY,NX)+PSIGrav_vr(0,NY,NX)+PSISoilOsmotic_vr(0,NY,NX)
    PSIST1            = PSISM1_vr(NUM(NY,NX),NY,NX)+PSIGrav_vr(NUM(NY,NX),NY,NX)+PSISoilOsmotic_vr(NUM(NY,NX),NY,NX)
    DarcyFlxLitR2Soil = AVCNDR*(PSIST0-PSIST1)*AREA(3,NUM(NY,NX),NY,NX)*CVRDW(NY,NX)*dts_HeatWatTP

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
      WatDarcyFloLitR2SoiMicP=AZMAX1(AMIN1(FLQZ,VLWatMicP1_vr(0,NY,NX)*dts_wat,VLairMicP1_vr(NUM(NY,NX),NY,NX)))
      FLQ2=AZMAX1(AMIN1(DarcyFlxLitR2Soil,VLWatMicP1_vr(0,NY,NX)*dts_wat,VLairMicP1_vr(NUM(NY,NX),NY,NX)))
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

    IF(VLairMicP_vr(NUM(NY,NX),NY,NX).LT.0.0_r8)THEN
      WatDarcyFloLitR2SoiMicP=WatDarcyFloLitR2SoiMicP+AZMIN1(AMAX1(-VLWatMicP1_vr(NUM(NY,NX),NY,NX)*dts_wat,&
        VLairMicP_vr(NUM(NY,NX),NY,NX)))
      FLQ2=FLQ2+AZMIN1(AMAX1(-VLWatMicP1_vr(NUM(NY,NX),NY,NX)*dts_wat,VLairMicP_vr(NUM(NY,NX),NY,NX)))
    ENDIF

    IF(WatDarcyFloLitR2SoiMicP.GT.0.0_r8)THEN
      !litter layer to soil
      HeatFlxLitR2Soi=cpw*TKSoi1_vr(0,NY,NX)*WatDarcyFloLitR2SoiMicP
    ELSE
      !soil to litter layer
      HeatFlxLitR2Soi=cpw*TKSoi1_vr(NUM(NY,NX),NY,NX)*WatDarcyFloLitR2SoiMicP
    ENDIF

    WatXChange2WatTable(3,NUM(NY,NX),NY,NX)=WatXChange2WatTable(3,NUM(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    if(abs(WatXChange2WatTable(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qWatXChange2WatTable(3,NUM(NY,NX),NY,NX)=',WatDarcyFloLitR2SoiMicP
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HeatFlow2Soili(3,NUM(NY,NX),NY,NX)=HeatFlow2Soili(3,NUM(NY,NX),NY,NX)+HeatFlxLitR2Soi
!    if(HeatFlow2Soili(3,NUM(NY,NX),NY,NX)>10._r8)then
!      write(*,*)'HeatFlow2Soili(3,NUM(NY,NX),NY,NX)',HeatFlow2Soili(3,NUM(NY,NX),NY,NX),HeatFlxLitR2Soi
!    endif
    WatFLow2LitR_col(NY,NX)   = WatFLow2LitR_col(NY,NX)-WatDarcyFloLitR2SoiMicP
    HeatFLoByWat2LitRi_col(NY,NX) = HeatFLoByWat2LitRi_col(NY,NX)-HeatFlxLitR2Soi
    WatFLo2LitrM(M,NY,NX)     = WatDarcyFloLitR2SoiMicP
!    write(*,*)I+J/24.,'if1',HeatFlxLitR2Soi,WatDarcyFloLitR2SoiMicP,&
!      safe_adb(HeatFlxLitR2Soi,(cpw*WatDarcyFloLitR2SoiMicP))
  ELSE
    !top layer is water
    WatDarcyFloLitR2SoiMicP                 = XVLMobileWatMicP(NY,NX)*dts_wat
    HeatFlxLitR2Soi                         = cpw*TKSoi1_vr(0,NY,NX)*WatDarcyFloLitR2SoiMicP
    WatXChange2WatTable(3,NUM(NY,NX),NY,NX) = WatXChange2WatTable(3,NUM(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    if(abs(WatXChange2WatTable(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qrWatXChange2WatTable(3,NUM(NY,NX),NY,NX)=',XVLMobileWatMicP(NY,NX),dts_wat
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HeatFlow2Soili(3,NUM(NY,NX),NY,NX)=HeatFlow2Soili(3,NUM(NY,NX),NY,NX)+HeatFlxLitR2Soi
!    if(HeatFlow2Soili(3,NUM(NY,NX),NY,NX)>10._r8)then
!      write(*,*)'HeatFlow2Soili(3,NUM(NY,NX),NY,NX)=',HeatFlow2Soili(3,NUM(NY,NX),NY,NX),HeatFlxLitR2Soi
!    endif
    WatFLow2LitR_col(NY,NX)   = WatFLow2LitR_col(NY,NX)-WatDarcyFloLitR2SoiMicP
    HeatFLoByWat2LitRi_col(NY,NX) = HeatFLoByWat2LitRi_col(NY,NX)-HeatFlxLitR2Soi
    WatFLo2LitrM(M,NY,NX)     = WatDarcyFloLitR2SoiMicP
  ENDIF

!     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
!     OF THE LITTER IS EXCEEDED
!
!     VOLPH1=air-filled macroporosity
!     FINHR,HFINHR=water,convective heat from litter to macropores
!     VLWatMicP1_vr(0,VWatLitRHoldCapcity=current,maximum litter water volume
!     WatXChange2WatTable,HFLWL=micropore water,heat flux
!     WatFLow2LitR,HFLWRL=total litter water,heat flux
!
  IF(VLairMacP1_vr(NUM(NY,NX),NY,NX).GT.0.0_r8 .AND.XVLMobileWatMicP(NY,NX).GT.0.0_r8)THEN
    WatFlowLitR2MacP                         = AMIN1(XVLMobileWatMicP(NY,NX)*dts_wat,VLairMacP1_vr(NUM(NY,NX),NY,NX))
    HeatFlowLitR2MacP                        = WatFlowLitR2MacP*cpw*TKSoi1_vr(0,NY,NX)
    ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX) = ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)+WatFlowLitR2MacP
    HeatFlow2Soili(3,NUM(NY,NX),NY,NX)       = HeatFlow2Soili(3,NUM(NY,NX),NY,NX)+HeatFlowLitR2MacP
!    if(HeatFlow2Soili(3,NUM(NY,NX),NY,NX)>10._r8)then
!      write(*,*)'HeatFlow2Soili(3,NUM(NY,NX),NY,NX)x=',HeatFlow2Soili(3,NUM(NY,NX),NY,NX),HeatFlowLitR2MacP
!    endif
    WatFLow2LitR_col(NY,NX)   = WatFLow2LitR_col(NY,NX)-WatFlowLitR2MacP
    HeatFLoByWat2LitRi_col(NY,NX) = HeatFLoByWat2LitRi_col(NY,NX)-HeatFlowLitR2MacP
!    write(*,*)I+J/24.,'surif',HeatFlowLitR2MacP,WatFlowLitR2MacP,safe_adb(HeatFlowLitR2MacP,WatFlowLitR2MacP*cpw)
  ENDIF
!  if(M==11)write(*,*)'SurfLitrSoilWaterExchange macP',M,NY,NX,HeatFLoByWat2LitRi_col(NY,NX)
  end subroutine SurfLitrSoilWaterExchange
!------------------------------------------------------------------------------------------

  subroutine InfilSRFRoffPartition(I,J,M,NY,NX,N1,N2)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NY,NX
  integer, intent(out):: N1,N2
  real(r8) :: LitrIceHeatFlxFrezPt,TK1X,ENGYR,VLWatMicP1X,VLHeatCapacityX
  real(r8) :: TFREEZ,TFLX,VWatExcess

  real(r8) :: WatExcessDetph,Q,HydraulicRadius,CrossSectVelocity

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
!
  TFREEZ          = -9.0959E+04_r8/(PSISM1_vr(0,NY,NX)-LtHeatIceMelt)
  VLWatMicP1X     = AZMAX1(VLWatMicP1_vr(0,NY,NX)+WatFLow2LitR_col(NY,NX))
  ENGYR           = VHeatCapacity1_vr(0,NY,NX)*TKSoi1_vr(0,NY,NX)
  VLHeatCapacityX = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP1X+cpi*VLiceMicP1_vr(0,NY,NX)

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    TK1X=(ENGYR+HeatFLoByWat2LitRi_col(NY,NX))/VLHeatCapacityX
  ELSE
    TK1X=TKSoi1_vr(0,NY,NX)
  ENDIF
  IF((TK1X.LT.TFREEZ .AND. VLWatMicP1_vr(0,NY,NX).GT.ZERO*VGeomLayer_vr(0,NY,NX)) &
    .OR.(TK1X.GT.TFREEZ .AND. VLiceMicP1_vr(0,NY,NX).GT.ZERO*VGeomLayer_vr(0,NY,NX)))THEN
    LitrIceHeatFlxFrezPt=VHeatCapacity1_vr(0,NY,NX)*(TFREEZ-TK1X) &
      /((1.0_r8+TFREEZ*6.2913E-03_r8)*(1.0_r8-0.10_r8*PSISM1_vr(0,NY,NX)))*dts_wat
    IF(LitrIceHeatFlxFrezPt.LT.0.0_r8)THEN
      !ice thaw
      TFLX=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1_vr(0,NY,NX)*dts_wat,LitrIceHeatFlxFrezPt)
    ELSE
      !water freeze
      TFLX=AMIN1(LtHeatIceMelt*VLWatMicP1X*dts_wat,LitrIceHeatFlxFrezPt)
    ENDIF
    LitrIceHeatFlxFrez(NY,NX) = TFLX
    LitrIceFlxThaw(NY,NX)     = -TFLX/LtHeatIceMelt
  ELSE
    LitrIceFlxThaw(NY,NX)     = 0.0_r8
    LitrIceHeatFlxFrez(NY,NX) = 0.0_r8
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
!     VWatStoreCapSurf=ground surface water retention capacity
!     VWatExcess=ponded water volume above surface retention capacity
!     D,R,S,V=depth,perimeter,slope,velocity of runoff
!     DIST=distance between source,destination
!     ZM=surface roughness height for runoff
!     Q=runoff from Mannings equation
!     QRM,QRV=runoff,velocity for erosion, solute transfer
!
  IF(XVLMobileWaterLitR_col(N2,N1).GT.VWatStoreCapSurf(N2,N1))THEN
    VWatExcess                    = XVLMobileWaterLitR_col(N2,N1)-VWatStoreCapSurf(N2,N1)
    WatExcessDetph                = VWatExcess/AREA(3,0,N2,N1)
    HydraulicRadius               = WatExcessDetph/2.828_r8
    CrossSectVelocity             = GaucklerManningVelocity(HydraulicRadius,SLOPE(0,N2,N1))/SoiSurfRoughness(N2,N1)  !1/s
    Q                             = CrossSectVelocity*WatExcessDetph*AREA(3,NUM(N2,N1),N2,N1)*3.6E+03_r8*dts_HeatWatTP
    VLWatMicP1X                   = AZMAX1(VLWatMicP1_vr(0,N2,N1)+LitrIceFlxThaw(N2,N1))
    RunoffVelocity(M,N2,N1)       = CrossSectVelocity
    WatFlux4ErosionM_2DH(M,N2,N1) = AMIN1(Q,VWatExcess*dts_wat,VLWatMicP1X*dts_wat) &
      *XVLMobileWatMicP(N2,N1)/XVLMobileWaterLitR_col(N2,N1)
  ELSE
    RunoffVelocity(M,N2,N1)       = 0.0_r8
    WatFlux4ErosionM_2DH(M,N2,N1) = 0.0_r8
  ENDIF
  end subroutine InfilSRFRoffPartition
!------------------------------------------------------------------------------------------

  subroutine LateralGridsHdryoExch(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
  !
  !between grid horizontal water flow
  implicit none
  integer, intent(in) :: M,NY,NX,NHE,NHW,NVS,NVN
  integer, intent(in) :: N1,N2   !source grid, with which the lateral exchange is computed 
  integer :: N,NN,N4,N5,N4B,N5B
  real(r8) :: ALT1,ALT2,ALTB,QRQ1
  integer, parameter :: idirew=1
  integer, parameter :: idirns=2
!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO  N=1,2
    DO  NN=1,2
      IF(N.EQ.iEastWestDirection)THEN
        !east-west
        IF((NX.EQ.NHE .AND. NN.EQ.1) .OR. (NX.EQ.NHW .AND. NN.EQ.2))THEN
          !at the boundary
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N4B=NX-1
          N5B=NY
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        !south-north
        IF((NY.EQ.NVS .AND. NN.EQ.1) .OR. (NY.EQ.NVN .AND. NN.EQ.2))THEN
          !at the boundary
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N4B=NX
          N5B=NY-1
        ENDIF
      ENDIF
!
!     ELEVATION OF EACH PAIR OF ADJACENT GRID CELLS
!
!     XVOLT,XVOLW=excess water+ice,water in destination grid cell
!     ALT1,ALT2=elevation of source,destination
!     QRQ1=equilibrium runoff
!     WatFlx2LitRByRunoff,HeatFlx2LitRByRunoff=runoff, convective heat from runoff
!     QR,HQR=hourly-accumulated runoff, convective heat from runoff
!     WatFlux4ErosionM=runoff water flux
      IF(WatFlux4ErosionM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
        ! there is runoff
        ! source grid elevation
        ALT1=Altitude_grid(N2,N1)+XVLMobileWaterLitR_col(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
!
!     EAST OR SOUTH RUNOFF
!
        IF(NN.EQ.1)THEN
          !destination grid elevation
          ALT2=Altitude_grid(N5,N4)+XVLMobileWaterLitR_col(N5,N4)/AREA(3,NU(N5,N4),N5,N4)
          IF(ALT1.GT.ALT2)THEN
            !source grid into dest grid
            QRQ1=AZMAX1(((ALT1-ALT2)*AREA(3,NUM(N2,N1),N2,N1) &
              *AREA(3,NU(N5,N4),N5,N4)-XVLMobileWaterLitR_col(N5,N4)*AREA(3,NUM(N2,N1),N2,N1) &
              +XVLMobileWaterLitR_col(N2,N1)*AREA(3,NU(N5,N4),N5,N4)) &
              /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5,N4),N5,N4)))
            WatFlx2LitRByRunoff(N,2,N5,N4)=AMIN1(QRQ1,WatFlux4ErosionM_2DH(M,N2,N1))*FSLOPE(N,N2,N1)
            HeatFlx2LitRByRunoff(N,2,N5,N4)=cpw*TKSoi1_vr(0,N2,N1)*WatFlx2LitRByRunoff(N,2,N5,N4)
            XGridSurfRunoff_2DH(N,2,N5,N4)=XGridSurfRunoff_2DH(N,2,N5,N4)+WatFlx2LitRByRunoff(N,2,N5,N4)
            HeatXGridBySurfRunoff_2DH(N,2,N5,N4)=HeatXGridBySurfRunoff_2DH(N,2,N5,N4)+HeatFlx2LitRByRunoff(N,2,N5,N4)
            QflxSurfRunoffM(M,N,2,N5,N4)=WatFlx2LitRByRunoff(N,2,N5,N4)
            IFLBM(M,N,2,N5,N4)=0
!            if(N2==1 .AND. N1==1 .and. M>=26)write(192,*)'LateralGridsHdryoExch',M,N,2,N5,N4,WatFlx2LitRByRunoff(N,2,N5,N4),QRQ1
          ELSE
            WatFlx2LitRByRunoff(N,2,N5,N4)=0.0_r8
            HeatFlx2LitRByRunoff(N,2,N5,N4)=0.0_r8
            QflxSurfRunoffM(M,N,2,N5,N4)=0.0_r8
            IFLBM(M,N,2,N5,N4)=1
          ENDIF
        ENDIF
!
!     WEST OR NORTH RUNOFF
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            ALTB=Altitude_grid(N5B,N4B)+XVLMobileWaterLitR_col(N5B,N4B)/AREA(3,NU(N5,N4B),N5B,N4B)
            IF(ALT1.GT.ALTB)THEN
              QRQ1=AZMAX1(((ALT1-ALTB)*AREA(3,NUM(N2,N1),N2,N1) &
                *AREA(3,NU(N5B,N4B),N5B,N4B)-XVLMobileWaterLitR_col(N5B,N4B) &
                *AREA(3,NUM(N2,N1),N2,N1) &
                +XVLMobileWaterLitR_col(N2,N1)*AREA(3,NU(N5B,N4B),N5B,N4B)) &
                /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5B,N4B),N5B,N4B)))
              WatFlx2LitRByRunoff(N,1,N5B,N4B)=AMIN1(QRQ1,WatFlux4ErosionM_2DH(M,N2,N1))*FSLOPE(N,N2,N1)
              HeatFlx2LitRByRunoff(N,1,N5B,N4B)=cpw*TKSoi1_vr(0,N2,N1)*WatFlx2LitRByRunoff(N,1,N5B,N4B)
              XGridSurfRunoff_2DH(N,1,N5B,N4B)=XGridSurfRunoff_2DH(N,1,N5B,N4B)+WatFlx2LitRByRunoff(N,1,N5B,N4B)
              HeatXGridBySurfRunoff_2DH(N,1,N5B,N4B)=HeatXGridBySurfRunoff_2DH(N,1,N5B,N4B)+HeatFlx2LitRByRunoff(N,1,N5B,N4B)
              QflxSurfRunoffM(M,N,1,N5B,N4B)=WatFlx2LitRByRunoff(N,1,N5B,N4B)
              IFLBM(M,N,1,N5B,N4B)=1
            ELSE
              WatFlx2LitRByRunoff(N,1,N5B,N4B)=0.0_r8
              HeatFlx2LitRByRunoff(N,1,N5B,N4B)=0.0_r8
              QflxSurfRunoffM(M,N,1,N5B,N4B)=0.0_r8
              IFLBM(M,N,1,N5B,N4B)=0
            ENDIF
          ENDIF
        ENDIF
      ELSE
        !there is no runoff
        WatFlx2LitRByRunoff(N,2,N5,N4)=0.0_r8
        HeatFlx2LitRByRunoff(N,2,N5,N4)=0.0_r8
        QflxSurfRunoffM(M,N,2,N5,N4)=0.0_r8
        IFLBM(M,N,2,N5,N4)=0
        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          WatFlx2LitRByRunoff(N,1,N5B,N4B)=0.0_r8
          HeatFlx2LitRByRunoff(N,1,N5B,N4B)=0.0_r8
          QflxSurfRunoffM(M,N,1,N5B,N4B)=0.0_r8
          IFLBM(M,N,1,N5B,N4B)=0
        ENDIF
      ENDIF

    ENDDO
  ENDDO

  end subroutine LateralGridsHdryoExch
!------------------------------------------------------------------------------------------

  subroutine AccumWaterVaporHeatFluxes(M,NY,NX,LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,&
    Radnet2Snow,VapXAir2TopLay)
  implicit none
  integer , intent(in) :: M,NY,NX
  real(r8), intent(in) :: LatentHeatAir2Sno,Radnet2Snow
  real(r8), intent(in) :: HeatSensAir2Snow,HeatSensEvap
  real(r8), intent(in) :: VapXAir2TopLay
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
! FLWM,WaterFlow2MacPM=water flux into soil micropore,macropore for use in TranspNoSalt.f
! VLWatMicPX1=VLWatMicP1 accounting for wetting front
!
  TLitrIceFlxThaw(NY,NX)=TLitrIceFlxThaw(NY,NX)+LitrIceFlxThaw(NY,NX)
  TLitrIceHeatFlxFrez(NY,NX)=TLitrIceHeatFlxFrez(NY,NX)+LitrIceHeatFlxFrez(NY,NX)
  WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX)=WaterFlowSoiMicP_3D(3,NUM(NY,NX),NY,NX) &
    +WatXChange2WatTable(3,NUM(NY,NX),NY,NX)
  WaterFlowSoiMicPX(3,NUM(NY,NX),NY,NX)=WaterFlowSoiMicPX(3,NUM(NY,NX),NY,NX) &
    +WatXChange2WatTableX(3,NUM(NY,NX),NY,NX)
  WaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)=WaterFlowMacP_3D(3,NUM(NY,NX),NY,NX) &
    +ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)
  HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX)=HeatFlow2Soil_3D(3,NUM(NY,NX),NY,NX) &
    +HeatFlow2Soili(3,NUM(NY,NX),NY,NX)
  WatFLo2Litr(NY,NX)       = WatFLo2Litr(NY,NX)+WatFLow2LitR_col(NY,NX)
  HeatFLo2LitrByWat(NY,NX) = HeatFLo2LitrByWat(NY,NX)+HeatFLoByWat2LitRi_col(NY,NX)

  HeatByRadiation_col(NY,NX)     = HeatByRadiation_col(NY,NX)+Radnet2LitGrnd+Radnet2Snow
  HeatSensAir2Surf_col(NY,NX)    = HeatSensAir2Surf_col(NY,NX)+HeatSensAir2Grnd+HeatSensAir2Snow
  HeatEvapAir2Surf_col(NY,NX)    = HeatEvapAir2Surf_col(NY,NX)+LatentHeatEvapAir2Grnd+LatentHeatAir2Sno
  HeatSensVapAir2Surf_col(NY,NX) = HeatSensVapAir2Surf_col(NY,NX)+HeatSensVapAir2Soi+HeatSensEvap
  !HeatNet2Surf_col > 0. into grnd
  HeatNet2Surf_col(NY,NX)=HeatNet2Surf_col(NY,NX)+Radnet2LitGrnd+Radnet2Snow &
    +HeatSensAir2Grnd+HeatSensAir2Snow+LatentHeatEvapAir2Grnd+LatentHeatAir2Sno &
    +HeatSensVapAir2Soi+HeatSensEvap

  !EVAPG=negative evaporation from ground/top soil layer
  !EVAPR=evaporation from litter layer   
  !EVAPSN=evaporation from snow, sublimation+evaporation
  VapXAir2GSurf_col(NY,NX)              = VapXAir2GSurf_col(NY,NX)+VapXAir2TopLay+VapXAir2LitR(NY,NX)+VapXAir2Sno(NY,NX)   !>0 into ground
  WaterFlow2MicPM(M,3,NUM(NY,NX),NY,NX) = WatXChange2WatTable(3,NUM(NY,NX),NY,NX)
  WaterFlow2MacPM(M,3,NUM(NY,NX),NY,NX) = ConvWaterFlowMacP_3D(3,NUM(NY,NX),NY,NX)
  TEvapXAir2Toplay_col(NY,NX)           = TEvapXAir2Toplay_col(NY,NX)+VapXAir2TopLay
  TEvapXAir2LitR_col(NY,NX)             = TevapXAir2LitR_col(NY,NX)+VapXAir2LitR(NY,NX)
  TEvapXAir2Snow_col(NY,NX)             = TEvapXAir2Snow_col(NY,NX)+VapXAir2Sno(NY,NX)
  end subroutine AccumWaterVaporHeatFluxes

!------------------------------------------------------------------------------------------

  subroutine InitSurfModel(I,J,M,NY,NX,ResistanceLitRLay,KSatReductByRainKineticEnergy)
  implicit none
  integer, intent(in) :: M   !soil heat-flow iteration id
  integer, intent(in) :: NY,NX,I,J
  real(r8),dimension(:,:),intent(in) :: ResistanceLitRLay(JY,JX)
  real(r8),intent(out):: KSatReductByRainKineticEnergy
  integer :: L  
  real(r8) :: scalar,THETWT,HFLQR1,FLQRS
  real(r8) :: FLQRH,VOLAT0,ENGYD
  real(r8) :: ENGYB,RAS,TScal4Aquadifsvity,THETWA
  real(r8) :: HV

! begin_execution
! INITIALIZE NET SURFACE FLUX ACCUMULATORS
!
! TQS1,TQW1,TQI1=net water and snowpack snow,water,ice runoff
! cumHeatFlx2LitRByRunoff,THQS1=net convective heat from surface water and snow runoff
! BAREW,CVRDW=fractions of soil,litter cover including free water+ice
! RAGS= boundary layer resistance at soil surface
! PARG=boundary layer conductance above soil surface
!
  call ZeroSnowFlux(NY,NX)

  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    BAREW(NY,NX)=AZMAX1(FracSurfBareSoil_col(NY,NX)-AMIN1(1.0_r8,&
      AZMAX1(XVLMobileWaterLitR_col(NY,NX)/MaxVLWatByLitR_col(NY,NX))))
  ELSE
    BAREW(NY,NX)=1.0_r8
  ENDIF
  CVRDW(NY,NX)=1.0_r8-BAREW(NY,NX)
  RAGS(NY,NX)=1.0_r8/(BAREW(NY,NX)/RAGR(NY,NX)+CVRDW(NY,NX)/ResistanceLitRLay(NY,NX))
  RAS=SnowBNDResistance(NY,NX)
  PARG(M,NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP/(RAGS(NY,NX)+RAS)
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
  Prec2LitR2=PrecAir2LitR*XNPR
  PrecHeat2LitR2=PrecHeatAir2LitR*XNPR
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

  VOLAT0=VLPoreLitR(NY,NX)-VLiceMicP1_vr(0,NY,NX)
  IF(VOLAT0.GT.ZEROS2(NY,NX).AND.VLsoiAirPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    !litter layer is not saturated
    THETWA                         = AZMAX1(AMIN1(1.0_r8,VLWatMicP1_vr(0,NY,NX)/VOLAT0))
    TScal4Aquadifsvity             = TEFAQUDIF(TKSoi1_vr(0,NY,NX))
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
! VWatStoreCapSurf=ground surface water retention capacity
! XVOLW=free surface water
! ZT=canopy height
! EnergyImpact4ErosionM=total energy impact for use in erosion.f
! EnergyImpact4Erosion=cumulative rainfall energy impact on soil surface
! KSatReductByRainKineticEnergy=reduction in soil surface Ksat from rainfall energy impact
! Note: A good reference for the following formula and alternatives
! is "Rainfall intensity-kinetic energy relationships for soil loss prediction",
! Kinnell, 1981

  IF(PRECD_col(NY,NX).GT.ZERO)THEN
    ENGYD=AZMAX1(8.95_r8+8.44_r8*LOG(PRECM_col(NY,NX)))
  ELSE
    ENGYD=0.0_r8
  ENDIF
  IF(PRECB_col(NY,NX).GT.ZERO)THEN
    ENGYB=AZMAX1(15.8_r8*SQRT(AMIN1(2.5_r8,CanopyHeight_col(NY,NX)))-5.87_r8)
  ELSE
    ENGYB=0.0_r8
  ENDIF

  IF(ENGYD+ENGYB.GT.ZERO)THEN
    HV=1.0E+03_r8*AZMAX1(XVLMobileWaterLitR_col(NY,NX)-VWatStoreCapSurf(NY,NX))/AREA(3,NUM(NY,NX),NY,NX)
    EnergyImpact4ErosionM(M,NY,NX)=(ENGYD*PRECD_col(NY,NX)+ENGYB*PRECB_col(NY,NX))*EXP(-2.0_r8*HV) &
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

  subroutine SurfaceRunoff(M,N,NN,N1,N2,M4,M5,RCHQF,XN)
  implicit none
  integer, intent(in) :: M,N,NN
  integer, intent(in) :: N1,N2  !source grid
  integer, intent(in) :: M4,M5  !dest grid
  real(r8),intent(in) :: RCHQF  !water flux scalar
  real(r8),intent(in) :: XN     !flow direction
  real(r8) :: ALT1,ALT2,DPTHW1,DPTHW2
  real(r8) :: VX
  !
  ! SURFACE BOUNDARY WATER FLUX
  !
  ! DPTHW1,DPTHW2=surface water depth of source,destination
  ! ALT1,ALT2=elevation of source,destination
  ! XVOLT=excess surface water+ice
  ! VWatStoreCapSurf=ground surface water retention capacity
  ! ExtWaterTable=natural water table depth
  ! QR1,HeatFlx2LitRByRunoff=runoff, convective heat from runoff
  ! QR,HQR=hourly-accumulated runoff, convective heat from runoff
  ! QRM,QRV=runoff,velocity for erosion, solute transfer
  ! XN=direction
  !
  ! RUNOFF
  !
  DPTHW1=XVLMobileWaterLitR_col(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
  DPTHW2=VWatStoreCapSurf(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
  !elevation at the source center
  ALT1=Altitude_grid(N2,N1)+DPTHW1      
  !elevation in the dest center
  ALT2=Altitude_grid(N2,N1)+DPTHW2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)  
  !grid elevation is higher than outside the grid, and in grid water layer higher than external water table
  IF(ALT1.GT.ALT2.AND.CumDepz2LayerBot_vr(NU(N2,N1)-1,N2,N1)-DPTHW1.LT.ExtWaterTable_col(N2,N1))THEN
    !out of grid (N2,N1), WatFlux4ErosionM is computed from surface physics model
    WatFlx2LitRByRunoff(N,NN,M5,M4)   = -XN*WatFlux4ErosionM_2DH(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HeatFlx2LitRByRunoff(N,NN,M5,M4)  = cpw*TKSoi1_vr(0,N2,N1)*WatFlx2LitRByRunoff(N,NN,M5,M4)
    XGridSurfRunoff_2DH(N,NN,M5,M4)  = XGridSurfRunoff_2DH(N,NN,M5,M4)+WatFlx2LitRByRunoff(N,NN,M5,M4)
    HeatXGridBySurfRunoff_2DH(N,NN,M5,M4) = HeatXGridBySurfRunoff_2DH(N,NN,M5,M4)+HeatFlx2LitRByRunoff(N,NN,M5,M4)
! RUNON
! water table in higher than grid surface (accouting for minimum water )
  ELSEIF(CumDepz2LayerBot_vr(NU(N2,N1)-1,N2,N1)-DPTHW1.GT.ExtWaterTable_col(N2,N1))THEN
    !elevation difference
    VX        = AZMIN1((ExtWaterTable_col(N2,N1)-CumDepz2LayerBot_vr(NU(N2,N1)-1,N2,N1)+DPTHW1)*AREA(3,NUM(N2,N1),N2,N1))
    WatFlux4ErosionM_2DH(M,N2,N1)     = VX*dts_wat
    RunoffVelocity(M,N2,N1)           = 0.0_r8
    WatFlx2LitRByRunoff(N,NN,M5,M4)   = -XN*WatFlux4ErosionM_2DH(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HeatFlx2LitRByRunoff(N,NN,M5,M4)  = cpw*TKSoi1_vr(0,N2,N1)*WatFlx2LitRByRunoff(N,NN,M5,M4)
    XGridSurfRunoff_2DH(N,NN,M5,M4)  = XGridSurfRunoff_2DH(N,NN,M5,M4)+WatFlx2LitRByRunoff(N,NN,M5,M4)
    HeatXGridBySurfRunoff_2DH(N,NN,M5,M4) = HeatXGridBySurfRunoff_2DH(N,NN,M5,M4)+HeatFlx2LitRByRunoff(N,NN,M5,M4)
  ELSE
    WatFlx2LitRByRunoff(N,NN,M5,M4)=0.0_r8
    HeatFlx2LitRByRunoff(N,NN,M5,M4)=0.0_r8
  ENDIF
  QflxSurfRunoffM(M,N,NN,M5,M4)=WatFlx2LitRByRunoff(N,NN,M5,M4)
  IFLBM(M,N,NN,M5,M4)=0
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
!     TFLWCI=net ice transfer to canopy, updated in hour1


  !convert water flux from m/hour to mm/hour
  PRECM_col(NY,NX)=m2mm*PrecRainAndIrrig_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  PRECD_col(NY,NX)=m2mm*(PrecRainAndIrrig_col(NY,NX)-TFLWCI(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  PRECB_col(NY,NX)=m2mm*(TFLWCI(NY,NX)-PrecIntceptByCanopy_col(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
!
!     RESIDUE WATER ABSORPTION CAPACITY
!
!     HCNDR=litter saturated hydraulic conductivity
!     DLYRR=litter depth
!
  HCNDR(NY,NX)=SatHydroCondLitR
  DLYRR_COL(NY,NX)=AMAX1(2.5E-03_r8,DLYR(3,0,NY,NX))
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

  IF(PrecRainAndIrrig_col(NY,NX).GT.0.0_r8 .OR. SnoFalPrec_col(NY,NX).GT.0.0_r8)THEN
  ! there is precipitation
    Rain4ToSno=(PrecRainAndIrrig_col(NY,NX)-PrecIntceptByCanopy_col(NY,NX))*FracSurfAsSnow(NY,NX)
    SnoFall=SnoFalPrec_col(NY,NX)                                                        !snowfall
    HeatByPrec=cps*TairK_col(NY,NX)*SnoFall+cpw*TairK_col(NY,NX)*Rain4ToSno                  !incoming heat flux from precipitations to snow-covered surface
    PrecThruFall=(PrecRainAndIrrig_col(NY,NX)-PrecIntceptByCanopy_col(NY,NX))*FracSurfSnoFree(NY,NX)       !incoming precipitation to snow-free surface
    PrecThrufall2LitR=PrecThruFall*FracSurfByLitR_col(NY,NX)                             !water flux to snow-free coverd by litter
    PrecHeat2LitR=cpw*TairK_col(NY,NX)*PrecThrufall2LitR                                 !heat flux to snow-free surface covered by litter
    PrecThrufall2Soil=PrecThruFall*FracSurfBareSoil_col(NY,NX)                          !heat flux to snow-free surface not covered by litter
    PrecHeat2Soil=cpw*TairK_col(NY,NX)*PrecThrufall2Soil
    PrecThrufall2SoiMicP=PrecThrufall2Soil*SoilFracAsMicP_vr(NUM(NY,NX),NY,NX)          !water flux to micropore
    PrecThrufall2SoiMacP=PrecThrufall2Soil*SoilFracAsMacP1_vr(NUM(NY,NX),NY,NX)         !water flux to macropore
  ELSE
  ! no precipitation
    Rain4ToSno=-PrecIntceptByCanopy_col(NY,NX)*FracSurfAsSnow(NY,NX)                   !
    SnoFall=0.0_r8
    HeatByPrec=cpw*TairK_col(NY,NX)*Rain4ToSno
    PrecThruFall=-PrecIntceptByCanopy_col(NY,NX)*FracSurfSnoFree(NY,NX)
    PrecThrufall2LitR=PrecThruFall*FracSurfByLitR_col(NY,NX)
    PrecHeat2LitR=cpw*TairK_col(NY,NX)*PrecThrufall2LitR
    PrecThrufall2Soil=PrecThruFall*FracSurfBareSoil_col(NY,NX)
    PrecHeat2Soil=cpw*TairK_col(NY,NX)*PrecThrufall2Soil
    PrecThrufall2SoiMicP=PrecThrufall2Soil*SoilFracAsMicP_vr(NUM(NY,NX),NY,NX)
    PrecThrufall2SoiMacP=PrecThrufall2Soil*SoilFracAsMacP1_vr(NUM(NY,NX),NY,NX)
  ENDIF
!
!     PRECIP ON SNOW ARRAYS EXPORTED TO TranspNoSalt.F, TranspSalt.F
!     FOR SOLUTE FLUX CALCULATIONS
!
!     SnoFalPrec,RainFalPrec,PrecAtm_col,PRECI=snow,rain,snow+rain,irrigation
!     VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
!     Rain2LitRSurf_col,Irrig2LitRSurf=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!
  IF(SnoFalPrec_col(NY,NX).GT.0.0_r8 .OR. (RainFalPrec(NY,NX).GT.0.0_r8 &
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
  SnowFallt(NY,NX)=SnoFall*dts_HeatWatTP
  Ice2Snowt(NY,NX)=0.0_r8
  Rain2Snowt(NY,NX)=Rain4ToSno*dts_HeatWatTP
  HeatFall2Snowt(NY,NX)=HeatByPrec*dts_HeatWatTP

  Prec2SoiMicP1(NY,NX)=PrecThrufall2SoiMicP*dts_HeatWatTP
  Prec2SoiMacP1(NY,NX)=PrecThrufall2SoiMacP*dts_HeatWatTP

  Prec2LitR1(NY,NX)=PrecThrufall2LitR*dts_HeatWatTP
  PrecHeat2SoiMicP1(NY,NX)=PrecHeat2Soil*dts_HeatWatTP
  PrecHeat2LitR1(NY,NX)=PrecHeat2LitR*dts_HeatWatTP

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

      call UpdateLitRAftRunoff(M,NY,NX)
      
    ENDDO D9790
  ENDDO D9795

  end subroutine UpdateSurfaceAtM
!------------------------------------------------------------------------------------------

  subroutine SumAftEnergyBalance(NY,NX,LWRadGrnd,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR,HeatFluxAir2Soi,NetWatFlx2SoiMicP)
  implicit none
  integer, intent(in)  :: NY,NX
  real(r8), intent(in) :: LWRadGrnd            !>0 into atmosphere, long wave radiation from ground surface
  real(r8), intent(in) :: VapXAir2TopLay
  real(r8), intent(in) :: HeatSensLitR2Soi1
  real(r8), intent(in) :: HeatSensVapLitR2Soi1
  real(r8), intent(in) :: EvapLitR2Soi1
  real(r8), intent(in) :: TotHeatAir2LitR
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
  WatNetFlo2TopSoiMicP=PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1
  if(abs(WatNetFlo2TopSoiMicP)>1.e20_r8)then
    write(*,*)'PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1',PrecNet2SoiMicP,VapXAir2TopLay,EvapLitR2Soi1
    call endrun(trim(mod_filename)//' at line',__LINE__)
  endif

  NetWatFlx2SoiMicP   = PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1
  NetWatFlx2SoiMacP   = PrecNet2SoiMacP
  cumHeatSensAir2Soil = PrecHeat2SoiNet+HeatFluxAir2Soi+HeatSensVapLitR2Soi1+HeatSensLitR2Soi1
  NetWatFlx2LitR      = PrecAir2LitR+VapXAir2LitR(NY,NX)-EvapLitR2Soi1
  CumHeatSensAir2LitR = PrecHeatAir2LitR+TotHeatAir2LitR-HeatSensVapLitR2Soi1-HeatSensLitR2Soi1
  FLWVLS              = (VLWatMicP1_vr(NUM(NY,NX),NY,NX)-VLWatMicPX1_vr(NUM(NY,NX),NY,NX))*dts_HeatWatTP
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
  real(r8), dimension(:,:),optional,intent(out) :: Hinfl2Soil(JY,JX)
  real(r8) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvap,VapXAir2TopLay
  integer :: N1,N2,NX,NY,L


  D9895: DO  NX=NHW,NHE
    D9890: DO  NY=NVN,NVS
!      write(*,*)'RunSurfacePhysModel',NY,NX,'M=',M,TKS_vr(0,NY,NX)

      call SurfaceEnergyModel(I,J,M,NX,NY,ResistanceLitRLay,KSatReductByRainKineticEnergy(NY,NX),&
        HeatFluxAir2Soi(NY,NX),LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,Radnet2Snow,&
        TopLayWatVol,VapXAir2TopLay)

    ! CAPILLARY EXCHANGE OF WATER BETWEEN SOIL SURFACE AND RESIDUE
      call SurfLitrSoilWaterExchange(I,J,M,NY,NX,KSatReductByRainKineticEnergy(NY,NX))

      call InfilSRFRoffPartition(I,J,M,NY,NX,N1,N2)
    !
      if(.not.ATS_cpl_mode)call LateralGridsHdryoExch(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)

      if(snowRedist_model)call SnowRedistribution(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)

    ! In ATS coupled mode we do not run the full Redist so we put the snow
    ! models here instead
      if (ATS_cpl_mode) then
        call SnowMassUpdate(I,J,NY,NX)
        call SnowpackLayering(I,J,NY,NX)
      end if

      call AccumWaterVaporHeatFluxes(M,NY,NX,LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,&
        Radnet2Snow,VapXAir2TopLay)

       call UpdateLitRB4RunoffM(I,J,M,NY,NX)
       if(present(Qinfl2MicP))Qinfl2MicP(NY,NX)=WatXChange2WatTable(3,NUM(NY,NX),NY,NX)
       if(present(Hinfl2Soil))Hinfl2Soil(NY,NX)=HeatFlow2Soili(3,NUM(NY,NX),NY,NX)

    ENDDO D9890
  ENDDO D9895

  end subroutine RunSurfacePhysModel

!------------------------------------------------------------------------------------------
  subroutine SurfaceEnergyModel(I,J,M,NX,NY,ResistanceLitRLay,KSatReductByRainKineticEnergy,&
    HeatFluxAir2Soi1,LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,Radnet2Snow,&
    TopLayWatVol,VapXAir2TopLay)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M     !soil heat-flow iteration
  integer, intent(in) :: NX,NY
  real(r8), dimension(:,:),intent(inout) :: ResistanceLitRLay(JY,JX)
  REAL(R8),INTENT(OUT) :: KSatReductByRainKineticEnergy,HeatFluxAir2Soi1
  real(r8), intent(out) :: Radnet2Snow,LatentHeatAir2Sno,HeatSensAir2Snow,HeatSensEvap
  real(r8),dimension(:,:),intent(inout) :: TopLayWatVol(JY,JX)
  real(r8), intent(out) :: VapXAir2TopLay
  integer :: N1,N2,L

  !ResistanceLitRLay is input
  call InitSurfModel(I,J,M,NY,NX,ResistanceLitRLay,KSatReductByRainKineticEnergy)

! updates ResistanceLitRLay
  call AtmLandSurfExchange(I,J,M,NY,NX,ResistanceLitRLay,TopLayWatVol,LatentHeatAir2Sno,&
    HeatSensEvap,HeatSensAir2Snow,Radnet2Snow,VapXAir2TopLay,HeatFluxAir2Soi1)

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
end module SurfPhysMod
