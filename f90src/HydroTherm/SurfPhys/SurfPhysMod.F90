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
  use SoilPhysParaMod
implicit none
  private
  character(len=*), parameter :: mod_filename=&
  __FILE__

  !surface model
  public :: StageSurfStateVars
  public :: SurfacePhysModel
  public :: SurfaceEnergyModel
  !
  public :: SurfaceRunoff
  public :: UpdateSurfaceAtM

! SatHydroCondLitR=saturated hydraulic conductivity of surface litter
! FENGYP=rate constant for restoring surface Ksat
! RARX=minimum boundary layer resistances of litter (h m-1)
! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
! RACX,RARX=minimum boundary layer resistances of canopy,litter (h m-1)

  real(r8), parameter :: FENGYP=1.0E-03_r8
  real(r8), parameter :: SatHydroCondLitR=25.0_r8
  real(r8), parameter :: RARX=0.0139_r8
  real(r8), parameter :: SoilEmisivity=0.97_r8        !soil emissivity
  real(r8), parameter :: SnowEmisivity=0.97_r8        !snowpack emissivity
  real(r8), parameter :: SurfLitREmisivity=0.97_r8    !surfce litter emissivity
  real(r8), parameter :: RACX=0.0139_r8    !total canopy boundary later resistance h/m  

  real(r8) :: PrecHeat2LitR2,Prec2LitR2  
  real(r8) :: HeatSensVapAir2Soi
  real(r8) :: HeatSensAir2Grnd,Radnet2LitGrnd,HeatFluxAir2Soi
  real(r8) :: LatentHeatEvapAir2Grnd,NetWatFlx2SoiMacP
  real(r8) :: CumHeatSensAir2LitR
  real(r8) :: cumHeatSensAir2Soil,NetWatFlx2LitR
  real(r8) :: WatNetFlo2TopSoiMicP,NetWatFlx2SoiMicP
  real(r8) :: PrecHeat2SoiNet,PrecNet2SoiMicP,PrecAir2LitR
  real(r8) :: PrecHeatAir2LitR,PrecNet2SoiMacP  
contains


  subroutine StageSurfStateVars(I,J,NHW,NHE,NVN,NVS,RAR1)

  use SnowPhysMod, only : CopySnowStates
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8),dimension(:,:),intent(OUT) :: RAR1
  character(len=*), parameter :: subn=trim(mod_filename)//'::StageSurfStateVars'
  integer :: NY,NX

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
    !make a local copy of the upper boundary index
!
!     ADJUST SURFACE ELEVATION USED IN RUNOFF FOR FREEZE-THAW, EROSION
!     AND SOC
!
!     ALTG,ALT=current,initial elevation of ground surface
!     CumDepth2LayerBottom(NUM(NY,NX)-1,=depth of ground surface
!     ENGYP=cumulative rainfall energy impact on soil surface
!
      ALTG(NY,NX)=ALT(NY,NX)-CumDepth2LayerBottom(NUM(NY,NX)-1,NY,NX)
      ENGYP(NY,NX)=ENGYP(NY,NX)*(1.0_r8-FENGYP)

      call CopySnowStates(NY,NX)

      call CopySurfaceVars(NY,NX)
!
      call PartionSurfaceFraction(NY,NX)

      call PartitionPrecip(NY,NX)

      call SurfaceRadiation(NY,NX)

      call SurfaceResistances(NY,NX,RAR1)

      call SetCanopyProperty(NY,NX)

    ENDDO D9990
  ENDDO D9995
  end subroutine StageSurfStateVars

!------------------------------------------------------------------------------------------  

  subroutine CopySurfaceVars(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: VWatLitrZ,TVWatIceLitR,VOLIRZ
!
! SET INITIAL SOIL VALUES
!
! LitrIceFlxThaw,LitrIceHeatFlxFrez=initialize surface litter freeze,thaw,latent heat
  LitrIceFlxThaw(NY,NX)=0.0_r8
  LitrIceHeatFlxFrez(NY,NX)=0.0_r8

!
! ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
!     FOR USE AT INTERNAL TIME STEP IN SURFACE LITTER
!
!     LWRadBySurf=longwave emission from litter surface
!     VLHeatCapacity=volumetric heat capacity of litter
!     VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of litter
!     VWatLitrX=maximum water retention by litter
!     XVOLT,XVOLW=free surface water+ice,water
!     VHeatCapLitR=min heat capacity for litter water,heat fluxes
!     VLitR=litter volume
!     THETW*,THETI*,THETP*=water,ice,air concentrations
!     PSISM*=litter matric water potential
!
  LWRadBySurf(NY,NX)=0.0_r8
  VLHeatCapacity(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VLWatMicP(0,NY,NX)+cpi*VLiceMicP(0,NY,NX)
  VLPoreLitR(NY,NX)=VLMicP(0,NY,NX)
  VLWatMicP1(0,NY,NX)=AZMAX1(VLWatMicP(0,NY,NX))
  VLiceMicP1(0,NY,NX)=AZMAX1(VLiceMicP(0,NY,NX))
  VLairMicP1(0,NY,NX)=AZMAX1(VLPoreLitR(NY,NX)-VLWatMicP1(0,NY,NX)-VLiceMicP1(0,NY,NX))
  VLWatMicPM(1,0,NY,NX)=VLWatMicP1(0,NY,NX)
  VLsoiAirPM(1,0,NY,NX)=VLairMicP1(0,NY,NX)
  TVWatIceLitR=VLWatMicP1(0,NY,NX)+VLiceMicP1(0,NY,NX)
  XVGeomLayer(NY,NX)=AZMAX1(TVWatIceLitR-VWatLitrX(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ=VLWatMicP1(0,NY,NX)/TVWatIceLitR*VWatLitrX(NY,NX)
    VOLIRZ=VLiceMicP1(0,NY,NX)/TVWatIceLitR*VWatLitrX(NY,NX)
    XVLWatMicP(NY,NX)=AZMAX1(VLWatMicP1(0,NY,NX)-VWatLitrZ)
    XVLiceMicP(NY,NX)=AZMAX1(VLiceMicP1(0,NY,NX)-VOLIRZ)
  ELSE
    XVLWatMicP(NY,NX)=0.0_r8
    XVLiceMicP(NY,NX)=0.0_r8
  ENDIF
  XVOLTM(1,NY,NX)=XVGeomLayer(NY,NX)
  XVLWatMicPM(1,NY,NX)=XVLWatMicP(NY,NX)
  XVLiceMicPM(1,NY,NX)=XVLiceMicP(NY,NX)
  IF(VLitR(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat(0,NY,NX)=AZMAX1t(VLWatMicP1(0,NY,NX)/VLitR(NY,NX))
    FracSoiPAsIce(0,NY,NX)=AZMAX1t(VLiceMicP1(0,NY,NX)/VLitR(NY,NX))
    FracSoiPAsAir(0,NY,NX)=AZMAX1t(VLairMicP1(0,NY,NX)/VLitR(NY,NX))*AZMAX1t((1.0_r8-XVGeomLayer(NY,NX)/VOLWD(NY,NX)))
  ELSE
    FracSoiPAsWat(0,NY,NX)=0.0_r8
    FracSoiPAsIce(0,NY,NX)=0.0_r8
    FracSoiPAsAir(0,NY,NX)=1.0
  ENDIF
  THETPM(1,0,NY,NX)=FracSoiPAsAir(0,NY,NX)
  PSISM1(0,NY,NX)=PSISoilMatricP(0,NY,NX)
  TKSoi1(0,NY,NX)=TKS(0,NY,NX)

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
  FracSurfAsSnow(NY,NX)=AMIN1(1.0_r8,SQRT((SnowDepth(NY,NX)/MinSnowDepth)))
  FracSurfSnoFree(NY,NX)=1.0_r8-FracSurfAsSnow(NY,NX)
  !if there is heat-wise significant litter layer
  IF(VLHeatCapacity(0,NY,NX).GT.VHeatCapLitR(NY,NX))THEN
    FracSurfAsBareSoi(NY,NX)=AMIN1(1.0_r8,AZMAX1(EXP(-0.8E-02_r8*(ORGC(0,NY,NX)/AREA(3,0,NY,NX)))))
  ELSE
    FracSurfAsBareSoi(NY,NX)=1.0_r8
  ENDIF

  FracSurfByLitR(NY,NX)=1.0_r8-FracSurfAsBareSoi(NY,NX)
  end subroutine PartionSurfaceFraction

!------------------------------------------------------------------------------------------
  subroutine SetCanopyProperty(NY,NX)      
  
  implicit none
  integer, intent(in) :: NY,NX

  real(r8), parameter :: SensHeatCondctance=1.25E-03_r8

  !TLEX=total latent heat flux x boundary layer resistance, [MJ m-1]
  !TSHX=total sensible heat flux x boundary layer resistance, [MJ m-1]
  !VPQ=vapor pressure in canopy air, 
  !TKQ=temperature in canopy air, Kelvin

  VPQ(NY,NX)=VPA(NY,NX)-TLEX(NY,NX)/(EvapLHTC*AREA(3,NUM(NY,NX),NY,NX))
  TKQ(NY,NX)=TairK(NY,NX)-TSHX(NY,NX)/(SensHeatCondctance*AREA(3,NUM(NY,NX),NY,NX))

  end subroutine SetCanopyProperty
!------------------------------------------------------------------------------------------

  subroutine SurfaceRadiation(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: THRYX,RADGX
!
!     INITIALIZE PARAMETERS, FLUXES FOR ENERGY EXCHANGE
!     AT SNOW, RESIDUE AND SOIL SURFACES
!
!     RADGX=shortwave radiation at ground surface
!     RADXW,RADXG,RADXR= shortwave radn at snowpack,soil,litter
!     FRADG=fraction of shortwave radiation at ground surface
!     FSNW,FSNX=fractions of snow,snow-free cover
!     BARE,CVRD=fractions of soil,litter cover
!     XNPS=internal time step for fluxes through snowpack
!     THRYX=longwave radiation at ground surface
!     LWRad2Snow,LWRad2Grnd,LWRad2LitR=longwave radn incident at snowpack,soil,litter
!     THRMW,THRMS,THRMR=longwave radn emitted by snowpack,soil,litter
!     SnowEmisivity,SoilEmisivity,SurfLitREmisivity=emissivity of snowpack,soil,litter surfaces
!     THS=sky longwave radiation
!     LWRadCanGPrev=longwave radiation emitted by canopy

  RADGX=RADG(NY,NX)*dts_HeatWatTP
  RADXW(NY,NX)=RADGX*FracSurfAsSnow(NY,NX)*XNPS
  RADXG(NY,NX)=RADGX*FracSurfSnoFree(NY,NX)*FracSurfAsBareSoi(NY,NX)      
  RADXR(NY,NX)=RADGX*FracSurfSnoFree(NY,NX)*FracSurfByLitR(NY,NX)*XNPR

  THRYX=(LWRadSky(NY,NX)*FRADG(NY,NX)+LWRadCanGPrev(NY,NX))*dts_HeatWatTP
  LWRad2Snow(NY,NX)=THRYX*FracSurfAsSnow(NY,NX)*XNPS
  LWRad2Grnd(NY,NX)=THRYX*FracSurfSnoFree(NY,NX)*FracSurfAsBareSoi(NY,NX)
  LWRad2LitR(NY,NX)=THRYX*FracSurfSnoFree(NY,NX)*FracSurfByLitR(NY,NX)*XNPR
  ! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
  !what is 2.04E-10_r8 it is stefan-boltzman constant converted into MJ m-2 K-4/per hour
  THRMW(NY,NX)=SnowEmisivity*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX)*dts_snohttp
  THRMS(NY,NX)=SoilEmisivity*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree(NY,NX)*FracSurfAsBareSoi(NY,NX)*dts_HeatWatTP
  THRMR(NY,NX)=SurfLitREmisivity*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree(NY,NX)*FracSurfByLitR(NY,NX)*dts_litrhtwtp
!
  end subroutine SurfaceRadiation
!------------------------------------------------------------------------------------------
  subroutine SurfaceResistances(NY,NX,RAR1)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), dimension(:,:), intent(out):: RAR1
  real(r8) :: FracSoiPAsAir0,DFVR  
  real(r8) :: PAREX,PARSX,RAS
  real(r8) :: ALFZ,UAG
!
!     BOUNDARY LAYER CONDUCTANCES FOR EXPORT TO TRNSFR.F
!

!     AERODYNAMIC RESISTANCE OF CANOPY TO SNOW/RESIDUE/SOIL
!     SURFACE ENERGY EXCHANGE WITH ATMOSPHERE
!
!     ALFZ=parameter for canopy effect on windspeed
!     FRADG=fraction of shortwave radiation at ground surface
!     BndlResistAboveCanG,RAC=isothermal blr above canopy, canopy blr
!     ZT,SoiSurfRoughnesst0=canopy, surface roughness heights
!     UA,UAG=windspeeds above,below canopy
!     VPQ,VPA=vapor pressure within,above canopy
!     TKQ,TairK=temperature within,above canopy
!     TLEX,TSHX=net latent,sensible heat fluxes x blrs from prev hour
!     VAP=latent heat of evaporation
!     1.25E-03=heat capacity of air
!     AREA=surface area of grid cell
!
  ALFZ=2.0_r8*(1.0_r8-FRADG(NY,NX))
  IF(BndlResistAboveCanG(NY,NX).GT.ZERO.AND.GridMaxCanopyHeight(NY,NX).GT.SoiSurfRoughnesst0(NY,NX).AND.ALFZ.GT.ZERO)THEN
    !if there is plant
    BndlResistCanG(NY,NX)=AMIN1(RACX,AZMAX1(GridMaxCanopyHeight(NY,NX)*EXP(ALFZ) &
      /(ALFZ/BndlResistAboveCanG(NY,NX))*AZMAX1(EXP(-ALFZ*SoiSurfRoughnesst0(NY,NX)/GridMaxCanopyHeight(NY,NX)) &
      -EXP(-ALFZ*(ZeroPlanDisp(NY,NX)+RoughHeight(NY,NX))/GridMaxCanopyHeight(NY,NX)))))
    UAG=UA(NY,NX)*EXP(-ALFZ)
  ELSE
    BndlResistCanG(NY,NX)=0.0_r8
    UAG=UA(NY,NX)
  ENDIF
!
!     AERODYNAMIC RESISTANCE OF SNOWPACK, RESIDUE AND SOIL
!     SURFACES TO ENERGY EXCHANGE WITH ATMOSPHERE
!     Soil Sci. Soc. Am. J. 48:25-32
!
!     RAR=porosity-unlimited litter blr
!     DLYRR=litter depth
!     WGSGR=vapor diffusivity in litter
!     RAG,RAGW,RAGR=isothermal blrs at ground,snowpack,litter surfaces
!     RARX=boundary layer resistance (blr) of litter surface
!     FracSoiPAsAir*=air-filled porosity of litter
!     DFVR=porosity limitation to diffusion through litter
!     POROQ=litter tortuosity
!     RAR1=porosity-limited litter blr
!     PAREX,PARSX=conductances for latent,sensible heat fluxes
!     PAREW,PARSW=conductances for snowpack latent,sensible heatfluxes
!     PAREG,PARSG=conductances for soil latent,sensible heat fluxes
!     PARER,PARSR=conductances for litter latent,sensible heat fluxes
!     XNPR=internal time step for fluxes through litter
!
  RAR(NY,NX)=DLYRR(NY,NX)/WGSGR(NY,NX)
  RAG(NY,NX)=BndlResistCanG(NY,NX)+BndlResistAboveCanG(NY,NX)
  RAGW(NY,NX)=RAG(NY,NX)
  RAGR(NY,NX)=RAG(NY,NX)+RARX
  RARG(NY,NX)=RAGR(NY,NX)
  FracSoiPAsAir0=AMAX1(ZERO2,FracSoiPAsAir(0,NY,NX))
  DFVR=FracSoiPAsAir0*POROQ*FracSoiPAsAir0/POROS(0,NY,NX)
  RAR1(NY,NX)=RAG(NY,NX)+RAR(NY,NX)/DFVR
  PAREX=AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP               !conductance for latent heat flux
  PARSX=1.25E-03_r8*AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP   !conductance for sensible heat flux
  
  PAREW(NY,NX)=PAREX*FracSurfAsSnow(NY,NX)*XNPS
  PARSW(NY,NX)=PARSX*FracSurfAsSnow(NY,NX)*XNPS
  PAREG(NY,NX)=PAREX*FracSurfSnoFree(NY,NX)
  PARER(NY,NX)=PAREX*FracSurfSnoFree(NY,NX)*XNPR*FracSurfByLitR(NY,NX)
  PARSG(NY,NX)=PARSX*FracSurfSnoFree(NY,NX)
  PARSR(NY,NX)=PARSX*FracSurfSnoFree(NY,NX)*XNPR*FracSurfByLitR(NY,NX)

!     PARR=boundary layer conductance above litter,soil surfaces
!
  RAS=SnowBNDResistance(NY,NX)

  PARR(NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*dts_HeatWatTP/(RAGR(NY,NX)+RAS)   !this includes snow layer resistance
  end subroutine SurfaceResistances

  !------------------------------------------------------------------------------------------

  subroutine SoilSRFEnerbyBalance(M,NY,NX,PSISV1,LWRadGrnd,RAR1,TopLayWatVol,VapXAir2TopLay)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(out):: PSISV1,LWRadGrnd
  real(r8),dimension(:,:), intent(inout) :: RAR1
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8) :: RAa,VaporSoi1,PARE,PARS,RAGX,RFLX0,RI,WPLX,WPX,FracSoiPAsAir0
  real(r8) :: VLWatGrnd,VLIceGrnd,DFVR
  real(r8) :: TKX1
  real(r8) :: HFLX0,AlbedoGrnd
! begin_execution
!
! PHYSICAL AND HYDRAULIC PROPERTIES OF SOIL SURFACE INCLUDING
! AIR AND WATER-FILLED POROSITY, AND WATER POTENTIAL USED IN
! FLUX CALCULATIONS
!
! THETY=current,hygroscopic water concentration
! POROS=porosity
! VLWatMicP1,VLSoilPoreMicPI=volume of micropore water
! VLSoilPoreMicPI=soil volume less macropore,rock
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

  call CalcSoilWatPotential(NY,NX,NX,NY,NUM(NY,NX),PSISM1(NUM(NY,NX),NY,NX))

  PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISoilOsmotic(NUM(NY,NX),NY,NX)
!
! IF(NX.EQ.4.AND.NY.EQ.5)THEN
!     WRITE(*,3232)'PSISV1',I,J,M,NX,NY,NUM(NY,NX),PSISV1
!    2,PSISM1(NUM(NY,NX),NY,NX),PSISoilOsmotic(NUM(NY,NX),NY,NX)
!    3,FracSoiPAsWat(NUM(NY,NX),NY,NX),THETW1,POROS(NUM(NY,NX),NY,NX)
!    4,PSL(NUM(NY,NX),NY,NX),LOG(THETW1),PSD(NUM(NY,NX),NY,NX)
!    5,VLWatMicP1(NUM(NY,NX),NY,NX),VLSoilMicP(NUM(NY,NX),NY,NX)
!    5,VLSoilPoreMicP(NUM(NY,NX),NY,NX)
!    5,SRP(NUM(NY,NX),NY,NX)
!3232  FORMAT(A8,6I4,20E14.6)
! ENDIF
!
! SOIL SURFACE ALBEDO, NET RADIATION
!
! VLWatMicP1,VLiceMicP1=water,ice volume in micopores
! VLWatMacP1,VLiceMacP1=water,ice volume in macopores
! AlbedoGrnd,ALBS=albedo of ground surface,soil
! BKVL=soil mass
! RADXG,LWRad2Grnd,Radnet2LitGrnd=incoming shortwave,longwave,net radiation
! LWRadGrnd,THRMS=emitted longwave radiation, emissivity
! TK1=soil temperature
! albedo of water and ice are set to 0.06, and 0.30 respectively
  VLWatGrnd=VLWatMicP1(NUM(NY,NX),NY,NX)+VLWatMacP1(NUM(NY,NX),NY,NX)
  VLIceGrnd=VLiceMicP1(NUM(NY,NX),NY,NX)+VLiceMacP1(NUM(NY,NX),NY,NX)

  IF((VLWatGrnd+VLIceGrnd).GT.ZEROS2(NY,NX))THEN
    !top soil has water or ice
    !ice albedo seems too low.
    AlbedoGrnd=(ALBS(NY,NX)*SoilMicPMassLayer(NUM(NY,NX),NY,NX)+0.06_r8*VLWatGrnd &
      +0.30_r8*VLIceGrnd)/(SoilMicPMassLayer(NUM(NY,NX),NY,NX)+VLWatGrnd+VLIceGrnd)
  ELSE
    AlbedoGrnd=ALBS(NY,NX)
  ENDIF
  !absorbed radiation
  !Radnet2LitGrnd=net radiation, after taking out outgoing surface layer radiation  
  !LWRadGrnd=emitted longwave radiation  
  RFLX0=(1.0_r8-AlbedoGrnd)*RADXG(NY,NX)+LWRad2Grnd(NY,NX)
  LWRadGrnd=THRMS(NY,NX)*TKSoi1(NUM(NY,NX),NY,NX)**4._r8
  Radnet2LitGrnd=RFLX0-LWRadGrnd
!
! AERODYNAMIC RESISTANCE ABOVE SOIL SURFACE INCLUDING
! RESISTANCE IMPOSED BY PLANT CANOPY
!
! FracSoiPAsAir*=air-filled porosity of soil
! DFVR=porosity limitation to diffusion through soil
! POROQ=soil tortuosity
! RAR1=porosity-limited litter blr
! RAGZ=combined soil+litter blr
! RI=Richardsons number
! RIB=isothermal RI
! TKQ,TK1=canopy air,soil temperature
! RAGZ,RAa=soil+litter blr
! RAGS=isothermal blr at ground surface
!
  FracSoiPAsAir0=AMAX1(ZERO,FracSoiPAsAir(0,NY,NX))
  DFVR=FracSoiPAsAir0*POROQ*FracSoiPAsAir0/POROS(0,NY,NX)
  RAR1(NY,NX)=RAG(NY,NX)+RAR(NY,NX)/DFVR
  RI=RichardsonNumber(RIB(NY,NX),TKQ(NY,NX),TKSoi1(NUM(NY,NX),NY,NX))
  RAGX=AMAX1(RAM,0.8_r8*RAGS(NY,NX),AMIN1(1.2_r8*RAGS(NY,NX),RAR1(NY,NX)/(1.0_r8-10.0_r8*RI)))
  RAGS(NY,NX)=RAGX
  RAa=RAGR(NY,NX)+RAGS(NY,NX)
! IF(I.EQ.63.AND.NX.EQ.1)THEN
!     WRITE(*,7776)'RAGX',I,J,M,NX,NY,RAGZ,FracSurfAsBareSoi(NY,NX),RAG(NY,NX)
!    2,CVRDW(NY,NX),RAR1,RI,RIB(NY,NX),TKQ(NY,NX),TKSoi1(NUM(NY,NX),NY,NX)
!    3,TKSoi1(0,NY,NX),RAGX,RAM,RAGS(NY,NX),RA
!    4,RAR(NY,NX),DFVR,FracSoiPAsAir0,POROQ,FracSoiPAsAir(0,NY,NX)
!    5,DLYRR(NY,NX),WGSGR(NY,NX)
!7776  FORMAT(A8,5I4,30E12.4)
! ENDIF
!
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!
! PARE,PARS=blcs for latent,sensible heat fluxes over soil
! PAREG,PARSG=conductances for latent,sensible heat fluxes
! RZ=minimum surface resistance
! VaporSoi1,VPQ=vapor pressure at soil surface, canopy air
! VapXAir2TopLay=negative of evaporation, or water vapor transfer from canopy air to top layer soil/lake
! LatentHeatEvapAir2Grnd=latent heat flux
! XH=rate constant
! VOLW2=soil water volume
! VAP=latent heat of evaporation
! HeatSensVapAir2Soi=convective heat of evaporation flux
!
  PARE=PAREG(NY,NX)/(RAa+RZ)
  PARS=PARSG(NY,NX)/RAa

  TKX1=TKSoi1(NUM(NY,NX),NY,NX)
  VaporSoi1=vapsat(TKX1)*EXP(18.0*PSISV1/(RGAS*TKX1))

  !evaporation, no more than what is available, g H2O
  VapXAir2TopLay=AMAX1(PARE*(VPQ(NY,NX)-VaporSoi1),-AZMAX1(TopLayWatVol(NY,NX)*dts_wat))
  !latent heat
  LatentHeatEvapAir2Grnd=VapXAir2TopLay*EvapLHTC
  IF(VapXAir2TopLay.LT.0.0_r8)THEN
    !evaporation
    HeatSensVapAir2Soi=VapXAir2TopLay*cpw*TKSoi1(NUM(NY,NX),NY,NX)
  ELSE
    !condensation
    HeatSensVapAir2Soi=VapXAir2TopLay*cpw*TKQ(NY,NX)
  ENDIF
  !take away water from evaporation
  TopLayWatVol(NY,NX)=TopLayWatVol(NY,NX)+VapXAir2TopLay
!
! SOLVE FOR SOIL SURFACE TEMPERATURE AT WHICH ENERGY
! BALANCE OCCURS, SOLVE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
! HeatSensAir2Grnd,LatentHeatEvapAir2Grnd,Radnet2LitGrnd=sensible,latent heat fluxes, net radiation
! HeatSensVapAir2Soi=convective heat flux from LatentHeatEvapAir2Grnd
! HeatFlux2Ground=storage heat flux
!
  HeatSensAir2Grnd=PARS*(TKQ(NY,NX)-TKSoi1(NUM(NY,NX),NY,NX))
  !net energy into soil, subtracting latent heat and sensible heat
  HFLX0=Radnet2LitGrnd+LatentHeatEvapAir2Grnd+HeatSensAir2Grnd
  !total heat plus convective heat 
  HeatFluxAir2Soi=HFLX0+HeatSensVapAir2Soi

  end subroutine SoilSRFEnerbyBalance

!------------------------------------------------------------------------------------------

  subroutine ExposedSoilFlux(M,NY,NX,RAR1,TopLayWatVol,VapXAir2TopLay)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),dimension(:,:),intent(inout) :: RAR1
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol  
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8) :: PSISV1
  real(r8) :: LWRadGrnd
  real(r8) :: HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR

! begin_execution
! Watch out for L, is its value defined?
  call SoilSRFEnerbyBalance(M,NY,NX,PSISV1,LWRadGrnd,RAR1,TopLayWatVol,VapXAir2TopLay)
!
  call SRFLitterEnergyBalance(M,NY,NX,PSISV1,Prec2LitR2,PrecHeat2LitR2,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR)
!
  call SumAftEnergyBalance(NY,NX,LWRadGrnd,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR)

  end subroutine ExposedSoilFlux

!------------------------------------------------------------------------------------------

  subroutine AtmLandSurfExchange(M,NY,NX,RAR1,TopLayWatVol,LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,&
    Radnet2Snow,VapXAir2TopLay)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), dimension(:,:), intent(inout) :: RAR1
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol  
  real(r8), intent(out) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvap
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8) :: cumHeatFlowSno2Soi
  real(r8) :: CumWatFlx2SoiMicP,CumWatFlx2SoiMacP
  real(r8) :: CumHeatFlow2LitR
  real(r8) :: CumWatXFlx2SoiMicP
  real(r8) :: CumWatFlow2LitR
  real(r8) :: HeatNetFlx2Snow

! begin_execution
  VapXAir2TopLay=0._r8
  !solve if there is significant snow layer 
  
  IF(VLSnowHeatCapM(M,1,NY,NX).GT.VLHeatCapSnowMN(NY,NX))THEN
!   VHCPW,VLHeatCapSnowMN=current, minimum snowpack heat capacities
    call SolveSnowpack(M,NY,NX,LatentHeatAir2Sno,Radnet2Snow,HeatSensEvap,HeatSensAir2Snow,HeatNetFlx2Snow,&
      CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP,CumWatFlow2LitR,&
      CumHeatFlow2LitR,cumHeatFlowSno2Soi)
  ENDIF
!
! ENERGY EXCHANGE AT SOIL SURFACE IF EXPOSED UNDER SNOWPACK
! FSNW,FSNX=fractions of snow,snow-free cover
  IF(FracSurfSnoFree(NY,NX).GT.0.0_r8.AND.(SoiBulkDensity(NUM(NY,NX),NY,NX).GT.ZERO.OR. &
    VLHeatCapacity(NUM(NY,NX),NY,NX).GT.VHCPNX(NY,NX)))THEN
    !Ground partically covered by snow
     
    call ExposedSoilFlux(M,NY,NX,RAR1,TopLayWatVol,VapXAir2TopLay)

  ELSE
!   ground is fully snow covered, thus no flux from soil & litter
    call SnowCoveredTopSoilFlux(NY,NX)
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
  WatXChange2WatTable(3,NUM(NY,NX),NY,NX)=CumWatFlx2SoiMicP+WatNetFlo2TopSoiMicP
  WatXChange2WatTableX(3,NUM(NY,NX),NY,NX)=CumWatXFlx2SoiMicP+NetWatFlx2SoiMicP
  ConvectWaterFlowMacP(3,NUM(NY,NX),NY,NX)=CumWatFlx2SoiMacP+NetWatFlx2SoiMacP
  HeatFlow2Soili(3,NUM(NY,NX),NY,NX)=cumHeatFlowSno2Soi+cumHeatSensAir2Soil
  WatFLow2LitR(NY,NX)=CumWatFlow2LitR+NetWatFlx2LitR
  HeatFLoByWat2LitRi(NY,NX)=CumHeatFlow2LitR+CumHeatSensAir2LitR

!  write(*,*)'AtmLandSurfExchange  init',M,NY,NX,HeatFLoByWat2LitRi(NY,NX),CumHeatFlow2LitR,CumHeatSensAir2LitR
  end subroutine AtmLandSurfExchange  
!------------------------------------------------------------------------------------------

  subroutine SnowCoveredTopSoilFlux(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

! begin_execution
  Radnet2LitGrnd=0.0_r8   !net radiation into soil
  LatentHeatEvapAir2Grnd=0.0_r8   !latent heat flux from air and topsoil
  HeatSensVapAir2Soi=0.0_r8   !convective heat flux from air and topsoil
  HeatSensAir2Grnd=0.0_r8   !sensible heat flux from air to topsoil
  HeatFluxAir2Soi=0.0_r8   !net heat flux into from air into topsoil

  WatNetFlo2TopSoiMicP=0.0_r8   !total water flux from air into soil
  NetWatFlx2SoiMicP=0.0_r8
  NetWatFlx2SoiMacP=0.0_r8        !total water flux from air into macropore
  cumHeatSensAir2Soil=0.0_r8        !total water associated heat flux from air into soil
  NetWatFlx2LitR=0.0_r8        !total water flux from air to litter
  CumHeatSensAir2LitR=0.0_r8       !total water associated heat flux from air to litter
  VapXAir2LitR(NY,NX)=0.0_r8  !evaporative flux from air into litter
  end subroutine SnowCoveredTopSoilFlux  
!------------------------------------------------------------------------------------------

  subroutine SurfLitrSoilWaterExchange(M,NY,NX,FKSAT)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FKSAT

  real(r8) :: THETW1,THETWR,PSIST1
  real(r8) :: PSIST0,HeatFlxLitR2Soi,FLQZ,FLQX
  real(r8) :: WatDarcyFloLitR2SoiMicP,FLQ2,CNDR,AVCNDR
  real(r8) :: HeatFlowLitR2MacP,WatFlowLitR2MacP    
  real(r8) :: CND1
  integer :: K0,K1

! begin_execution
! CNDR,HCNDR=current,saturated litter hydraulic conductivity
! PSISE,PSISM1(0,=air entry,current litter water potential
! VLWatMicP1(0,VWatLitrX=current,maximum litter water volume
! CND1,HydroCond3D=soil hydraulic conductivity
! FKSAT=reduction in soil surface Ksat from rainfall energy impact
! K1=soil relative water-filled porosity
! THETWX,POROS=soil water content,porosity
! AVCNDR=litter-soil hydraulic conductance
! DLYRR,DLYR=litter,soil thicknesses
! WatDarcyFloLitR2SoiMicP,FLQL=litter-soil water flux unltd,ltd by water
! dts_HeatWatTP=time step of flux calculations
! CVRD=fraction of litter cover
! PSISM1(NUM=soil matric water potential
! VLWatMicP1(NUM=soil water volume
! HFLQL=convective heat from litter-soil water flux
! FLWL,HFLWL=micropore water,heat flux
! WatFLow2LitR,HFLWRL=total litter water,heat flux
! FLWRM=litter-soil water flux for solute transfer in trnsfr.f
! CND1,CNDL=hydraulic conductivity of source,destination layer
! HydroCond3D=lateral(1,2),vertical(3) micropore hydraulic conductivity
!
  IF(SoiBulkDensity(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    IF(VWatLitrX(NY,NX).GT.ZEROS2(NY,NX))THEN
      !surface litter holds water
      THETWR=AMIN1(VWatLitrX(NY,NX),VLWatMicP1(0,NY,NX))/VLitR(NY,NX)
    ELSE
      THETWR=POROS0(NY,NX)
    ENDIF
    THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX),AMIN1(POROS(NUM(NY,NX),NY,NX) &
      ,safe_adb(VLWatMicP1(NUM(NY,NX),NY,NX),VLSoilMicP(NUM(NY,NX),NY,NX))))
    !litter layer  
    K0=MAX(1,MIN(100,INT(100.0*(AZMAX1(POROS0(NY,NX)-THETWR))/POROS0(NY,NX))+1))
    !topsoil layer
    K1=MAX(1,MIN(100,INT(100.0*(AZMAX1(POROS(NUM(NY,NX),NY,NX)-THETW1))/POROS(NUM(NY,NX),NY,NX))+1))
    CNDR=HydroCond3D(3,K0,0,NY,NX)
    CND1=HydroCond3D(3,K1,NUM(NY,NX),NY,NX)*FKSAT
    AVCNDR=2.0_r8*CNDR*CND1/(CNDR*DLYR(3,NUM(NY,NX),NY,NX)+CND1*DLYRR(NY,NX))
    PSIST0=PSISM1(0,NY,NX)+PSIGrav(0,NY,NX)+PSISoilOsmotic(0,NY,NX)
    PSIST1=PSISM1(NUM(NY,NX),NY,NX)+PSIGrav(NUM(NY,NX),NY,NX)+PSISoilOsmotic(NUM(NY,NX),NY,NX)
    !FLQX=water flux from litter layer into the topsoil
    FLQX=AVCNDR*(PSIST0-PSIST1)*AREA(3,NUM(NY,NX),NY,NX)*CVRDW(NY,NX)*dts_HeatWatTP

    IF(FLQX.GE.0.0_r8)THEN
      !flow from liter to soil
      IF(THETWR.GT.THETS(0,NY,NX))THEN
        FLQZ=FLQX+AMIN1((THETWR-THETS(0,NY,NX))*VLitR(NY,NX),&
          AZMAX1((THETS(NUM(NY,NX),NY,NX)-THETW1)*VLSoilMicP(NUM(NY,NX),NY,NX)))*dts_wat
      ELSE
        FLQZ=FLQX
      ENDIF
      WatDarcyFloLitR2SoiMicP=AZMAX1(AMIN1(FLQZ,VLWatMicP1(0,NY,NX)*dts_wat,VLairMicP1(NUM(NY,NX),NY,NX)))
      FLQ2=AZMAX1(AMIN1(FLQX,VLWatMicP1(0,NY,NX)*dts_wat,VLairMicP1(NUM(NY,NX),NY,NX)))
    ELSE
      !from soil to liter
      IF(THETW1.GT.THETS(NUM(NY,NX),NY,NX))THEN
        FLQZ=FLQX+AMAX1((THETS(NUM(NY,NX),NY,NX)-THETW1)*VLSoilMicP(NUM(NY,NX),NY,NX),&
          AZMIN1((THETWR-THETS(0,NY,NX))*VLitR(NY,NX)))*dts_wat
      ELSE
        FLQZ=FLQX
      ENDIF
      WatDarcyFloLitR2SoiMicP=AZMIN1(AMAX1(FLQZ,-VLWatMicP1(NUM(NY,NX),NY,NX)*dts_wat,-VLairMicP1(0,NY,NX)))
      FLQ2=AZMIN1(AMAX1(FLQX,-VLWatMicP1(NUM(NY,NX),NY,NX)*dts_wat,-VLairMicP1(0,NY,NX)))
    ENDIF

    IF(VLairMicP(NUM(NY,NX),NY,NX).LT.0.0_r8)THEN
      WatDarcyFloLitR2SoiMicP=WatDarcyFloLitR2SoiMicP+AZMIN1(AMAX1(-VLWatMicP1(NUM(NY,NX),NY,NX)*dts_wat,&
        VLairMicP(NUM(NY,NX),NY,NX)))
      FLQ2=FLQ2+AZMIN1(AMAX1(-VLWatMicP1(NUM(NY,NX),NY,NX)*dts_wat,VLairMicP(NUM(NY,NX),NY,NX)))
    ENDIF

    IF(WatDarcyFloLitR2SoiMicP.GT.0.0_r8)THEN
      HeatFlxLitR2Soi=cpw*TKSoi1(0,NY,NX)*WatDarcyFloLitR2SoiMicP
    ELSE
      HeatFlxLitR2Soi=cpw*TKSoi1(NUM(NY,NX),NY,NX)*WatDarcyFloLitR2SoiMicP
    ENDIF
    WatXChange2WatTable(3,NUM(NY,NX),NY,NX)=WatXChange2WatTable(3,NUM(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    if(abs(WatXChange2WatTable(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qWatXChange2WatTable(3,NUM(NY,NX),NY,NX)=',WatDarcyFloLitR2SoiMicP
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HeatFlow2Soili(3,NUM(NY,NX),NY,NX)=HeatFlow2Soili(3,NUM(NY,NX),NY,NX)+HeatFlxLitR2Soi
    WatFLow2LitR(NY,NX)=WatFLow2LitR(NY,NX)-WatDarcyFloLitR2SoiMicP
    HeatFLoByWat2LitRi(NY,NX)=HeatFLoByWat2LitRi(NY,NX)-HeatFlxLitR2Soi
    WatFLo2LitrM(M,NY,NX)=WatDarcyFloLitR2SoiMicP

  ELSE
    WatDarcyFloLitR2SoiMicP=XVLWatMicP(NY,NX)*dts_wat
    HeatFlxLitR2Soi=cpw*TKSoi1(0,NY,NX)*WatDarcyFloLitR2SoiMicP
    WatXChange2WatTable(3,NUM(NY,NX),NY,NX)=WatXChange2WatTable(3,NUM(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    if(abs(WatXChange2WatTable(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qrWatXChange2WatTable(3,NUM(NY,NX),NY,NX)=',XVLWatMicP(NY,NX),dts_wat
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HeatFlow2Soili(3,NUM(NY,NX),NY,NX)=HeatFlow2Soili(3,NUM(NY,NX),NY,NX)+HeatFlxLitR2Soi
    WatFLow2LitR(NY,NX)=WatFLow2LitR(NY,NX)-WatDarcyFloLitR2SoiMicP
    HeatFLoByWat2LitRi(NY,NX)=HeatFLoByWat2LitRi(NY,NX)-HeatFlxLitR2Soi

    WatFLo2LitrM(M,NY,NX)=WatDarcyFloLitR2SoiMicP

  ENDIF
!  write(*,*)'SurfLitrSoilWaterExchange micP',M,NY,NX,HeatFLoByWat2LitRi(NY,NX)

!     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
!     OF THE LITTER IS EXCEEDED
!
!     VOLPH1=air-filled macroporosity
!     FINHR,HFINHR=water,convective heat from litter to macropores
!     VLWatMicP1(0,VWatLitrX=current,maximum litter water volume
!     WatXChange2WatTable,HFLWL=micropore water,heat flux
!     WatFLow2LitR,HFLWRL=total litter water,heat flux
!
  IF(VLairMacP1(NUM(NY,NX),NY,NX).GT.0.0_r8 .AND.XVLWatMicP(NY,NX).GT.0.0_r8)THEN
    WatFlowLitR2MacP=AMIN1(XVLWatMicP(NY,NX)*dts_wat,VLairMacP1(NUM(NY,NX),NY,NX))
    HeatFlowLitR2MacP=WatFlowLitR2MacP*cpw*TKSoi1(0,NY,NX)
    ConvectWaterFlowMacP(3,NUM(NY,NX),NY,NX)=ConvectWaterFlowMacP(3,NUM(NY,NX),NY,NX)+WatFlowLitR2MacP
    HeatFlow2Soili(3,NUM(NY,NX),NY,NX)=HeatFlow2Soili(3,NUM(NY,NX),NY,NX)+HeatFlowLitR2MacP
    WatFLow2LitR(NY,NX)=WatFLow2LitR(NY,NX)-WatFlowLitR2MacP
    HeatFLoByWat2LitRi(NY,NX)=HeatFLoByWat2LitRi(NY,NX)-HeatFlowLitR2MacP
  ENDIF
!  write(*,*)'SurfLitrSoilWaterExchange macP',M,NY,NX,HeatFLoByWat2LitRi(NY,NX)
  end subroutine SurfLitrSoilWaterExchange
!------------------------------------------------------------------------------------------

  subroutine InfilSRFRoffPartition(M,NY,NX,N1,N2)
  implicit none
  integer, intent(in) :: M,NY,NX
  integer, intent(out):: N1,N2
  real(r8) :: LitrIceHeatFlxFrezPt,TK1X,ENGYR,VLWatMicP1X,VLHeatCapacityX
  real(r8) :: TFREEZ,TFLX,VX

  real(r8) :: D,Q,R,V

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
!     ORGC=litter organic C
!     HFLWRL=total litter conductive, convective heat flux
!     LitrIceHeatFlxFrezPt,TFLX=unltd,ltd latent heat from freeze-thaw
!     LitrIceHeatFlxFrez,LitrIceFlxThaw=litter water,latent heat flux from freeze-thaw
!
  TFREEZ=-9.0959E+04_r8/(PSISM1(0,NY,NX)-LtHeatIceMelt)
  VLWatMicP1X=AZMAX1(VLWatMicP1(0,NY,NX)+WatFLow2LitR(NY,NX))
  ENGYR=VLHeatCapacity(0,NY,NX)*TKSoi1(0,NY,NX)
  VLHeatCapacityX=cpo*ORGC(0,NY,NX)+cpw*VLWatMicP1X+cpi*VLiceMicP1(0,NY,NX)

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    TK1X=(ENGYR+HeatFLoByWat2LitRi(NY,NX))/VLHeatCapacityX
  ELSE
    TK1X=TKSoi1(0,NY,NX)
  ENDIF

  IF((TK1X.LT.TFREEZ.AND.VLWatMicP1(0,NY,NX).GT.ZERO*VGeomLayer(0,NY,NX)) &
    .OR.(TK1X.GT.TFREEZ.AND.VLiceMicP1(0,NY,NX).GT.ZERO*VGeomLayer(0,NY,NX)))THEN
    LitrIceHeatFlxFrezPt=VLHeatCapacity(0,NY,NX)*(TFREEZ-TK1X) &
      /((1.0_r8+TFREEZ*6.2913E-03_r8)*(1.0_r8-0.10_r8*PSISM1(0,NY,NX)))*dts_wat
    IF(LitrIceHeatFlxFrezPt.LT.0.0_r8)THEN
      !ice thaw
      TFLX=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1(0,NY,NX)*dts_wat,LitrIceHeatFlxFrezPt)
    ELSE
      !water freeze
      TFLX=AMIN1(LtHeatIceMelt*VLWatMicP1X*dts_wat,LitrIceHeatFlxFrezPt)
    ENDIF
    LitrIceHeatFlxFrez(NY,NX)=TFLX
    LitrIceFlxThaw(NY,NX)=-TFLX/LtHeatIceMelt
  ELSE
    LitrIceFlxThaw(NY,NX)=0.0_r8
    LitrIceHeatFlxFrez(NY,NX)=0.0_r8
  ENDIF
!
!     THICKNESS OF WATER FILMS IN LITTER AND SOIL SURFACE
!     FROM WATER POTENTIALS FOR GAS EXCHANGE IN TRNSFR.F
!
  IF(VLHeatCapacity(0,NY,NX).GT.VHeatCapLitR(NY,NX))THEN
    FILM(M,0,NY,NX)=FilmThickness(PSISM1(0,NY,NX), is_top_layer=.true.)
  ELSE
    FILM(M,0,NY,NX)=1.0E-03_r8
  ENDIF
  FILM(M,NUM(NY,NX),NY,NX)=FilmThickness(PSISM1(NUM(NY,NX),NY,NX),is_top_layer=.true.)
!
!     OVERLAND FLOW WHEN WATER STORAGE CAPACITY
!     OF THE SOIL SURFACE PLUS MACROPORES IS EXCEEDED
!
  N1=NX
  N2=NY
!
!     SURFACE WATER FLUX
!
!     N2,N1=NY,NX of source grid cell
!     XVOLT,XVOLW=excess water+ice,water in source grid cell
!     VOLWG=ground surface water retention capacity
!     VX=ponded water volume above surface retention capacity
!     D,R,S,V=depth,perimeter,slope,velocity of runoff
!     DIST=distance between source,destination
!     ZM=surface roughness height for runoff
!     Q=runoff from Mannings equation
!     QRM,QRV=runoff,velocity for erosion, solute transfer
!
  IF(XVGeomLayer(N2,N1).GT.VOLWG(N2,N1))THEN
    VX=XVGeomLayer(N2,N1)-VOLWG(N2,N1)
    D=VX/AREA(3,0,N2,N1)
    R=D/2.828_r8
    V=R**0.67_r8*safe_adb(SQRT(SLOPE(0,N2,N1)),SoiSurfRoughness(N2,N1))
    Q=V*D*AREA(3,NUM(N2,N1),N2,N1)*3.6E+03_r8*dts_HeatWatTP
    VLWatMicP1X=AZMAX1(VLWatMicP1(0,N2,N1)+LitrIceFlxThaw(N2,N1))
    QRM(M,N2,N1)=AMIN1(Q,VX*dts_wat,VLWatMicP1X*dts_wat)*XVLWatMicP(N2,N1)/XVGeomLayer(N2,N1)
    QRV(M,N2,N1)=V
  ELSE
    QRM(M,N2,N1)=0.0_r8
    QRV(M,N2,N1)=0.0_r8
  ENDIF
  end subroutine InfilSRFRoffPartition
!------------------------------------------------------------------------------------------

  subroutine LateralHydroExchange(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
  implicit none
  integer, intent(in) :: M,NY,NX,NHE,NHW,NVS,NVN
  integer, intent(in) :: N1,N2
  integer :: N,NN,N4,N5,N4B,N5B
  real(r8) :: ALT1,ALT2,ALTB,QRQ1
  integer, parameter :: idirew=1
  integer, parameter :: idirns=2
!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO  N=1,2
    DO  NN=1,2
      IF(N.EQ.idirew)THEN
        !east-west
        IF(NX.EQ.NHE.AND.NN.EQ.1.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
          !at the boundary
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N4B=NX-1
          N5B=NY
        ENDIF
      ELSEIF(N.EQ.idirns)THEN
        !south-north
        IF(NY.EQ.NVS.AND.NN.EQ.1.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
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
!     QR1,HQR1=runoff, convective heat from runoff
!     QR,HQR=hourly-accumulated runoff, convective heat from runoff
!     QRM=runoff water flux
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        ! there is runoff
        ! source grid elevation
        ALT1=ALTG(N2,N1)+XVGeomLayer(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
!
!     EAST OR SOUTH RUNOFF
!
        IF(NN.EQ.1)THEN
          !destination grid elevation
          ALT2=ALTG(N5,N4)+XVGeomLayer(N5,N4)/AREA(3,NU(N5,N4),N5,N4)
          IF(ALT1.GT.ALT2)THEN
            QRQ1=AZMAX1(((ALT1-ALT2)*AREA(3,NUM(N2,N1),N2,N1) &
              *AREA(3,NU(N5,N4),N5,N4)-XVGeomLayer(N5,N4)*AREA(3,NUM(N2,N1),N2,N1) &
              +XVGeomLayer(N2,N1)*AREA(3,NU(N5,N4),N5,N4)) &
              /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5,N4),N5,N4)))
            QR1(N,2,N5,N4)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
            HQR1(N,2,N5,N4)=cpw*TKSoi1(0,N2,N1)*QR1(N,2,N5,N4)
            QR(N,2,N5,N4)=QR(N,2,N5,N4)+QR1(N,2,N5,N4)
            HQR(N,2,N5,N4)=HQR(N,2,N5,N4)+HQR1(N,2,N5,N4)
            QRMN(M,N,2,N5,N4)=QR1(N,2,N5,N4)
            IFLBM(M,N,2,N5,N4)=0
          ELSE
            QR1(N,2,N5,N4)=0.0_r8
            HQR1(N,2,N5,N4)=0.0_r8
            QRMN(M,N,2,N5,N4)=0.0_r8
            IFLBM(M,N,2,N5,N4)=1
          ENDIF
        ENDIF
!
!     WEST OR NORTH RUNOFF
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            ALTB=ALTG(N5B,N4B)+XVGeomLayer(N5B,N4B)/AREA(3,NU(N5,N4B),N5B,N4B)
            IF(ALT1.GT.ALTB)THEN
              QRQ1=AZMAX1(((ALT1-ALTB)*AREA(3,NUM(N2,N1),N2,N1) &
                *AREA(3,NU(N5B,N4B),N5B,N4B)-XVGeomLayer(N5B,N4B) &
                *AREA(3,NUM(N2,N1),N2,N1) &
                +XVGeomLayer(N2,N1)*AREA(3,NU(N5B,N4B),N5B,N4B)) &
                /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5B,N4B),N5B,N4B)))
              QR1(N,1,N5B,N4B)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
              HQR1(N,1,N5B,N4B)=cpw*TKSoi1(0,N2,N1)*QR1(N,1,N5B,N4B)
              QR(N,1,N5B,N4B)=QR(N,1,N5B,N4B)+QR1(N,1,N5B,N4B)
              HQR(N,1,N5B,N4B)=HQR(N,1,N5B,N4B)+HQR1(N,1,N5B,N4B)
              QRMN(M,N,1,N5B,N4B)=QR1(N,1,N5B,N4B)
              IFLBM(M,N,1,N5B,N4B)=1
            ELSE
              QR1(N,1,N5B,N4B)=0.0_r8
              HQR1(N,1,N5B,N4B)=0.0_r8
              QRMN(M,N,1,N5B,N4B)=0.0_r8
              IFLBM(M,N,1,N5B,N4B)=0
            ENDIF
          ENDIF
        ENDIF
      ELSE
        !there is no runoff
        QR1(N,2,N5,N4)=0.0_r8
        HQR1(N,2,N5,N4)=0.0_r8
        QRMN(M,N,2,N5,N4)=0.0_r8
        IFLBM(M,N,2,N5,N4)=0
        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          QR1(N,1,N5B,N4B)=0.0_r8
          HQR1(N,1,N5B,N4B)=0.0_r8
          QRMN(M,N,1,N5B,N4B)=0.0_r8
          IFLBM(M,N,1,N5B,N4B)=0
        ENDIF
      ENDIF

    ENDDO
  ENDDO

  call SnowRedistrub(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)

  end subroutine LateralHydroExchange
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
! HEATI,HeatEvapAir2Surf,HeatSensAir2Surf,HEATG=net radiation,latent,sensible,storage heat
! VapXAir2GSurf=total evaporation
! FLWM,WaterFlow2MacPM=water flux into soil micropore,macropore for use in trnsfr.f
! VLWatMicPX1=VLWatMicP1 accounting for wetting front
!
  TLitrIceFlxThaw(NY,NX)=TLitrIceFlxThaw(NY,NX)+LitrIceFlxThaw(NY,NX)
  TLitrIceHeatFlxFrez(NY,NX)=TLitrIceHeatFlxFrez(NY,NX)+LitrIceHeatFlxFrez(NY,NX)
  WaterFlowSoiMicP(3,NUM(NY,NX),NY,NX)=WaterFlowSoiMicP(3,NUM(NY,NX),NY,NX)+WatXChange2WatTable(3,NUM(NY,NX),NY,NX)
  WaterFlowSoiMicPX(3,NUM(NY,NX),NY,NX)=WaterFlowSoiMicPX(3,NUM(NY,NX),NY,NX)+WatXChange2WatTableX(3,NUM(NY,NX),NY,NX)
  WaterFlowMacP(3,NUM(NY,NX),NY,NX)=WaterFlowMacP(3,NUM(NY,NX),NY,NX)+ConvectWaterFlowMacP(3,NUM(NY,NX),NY,NX)
  HeatFlow2Soil(3,NUM(NY,NX),NY,NX)=HeatFlow2Soil(3,NUM(NY,NX),NY,NX)+HeatFlow2Soili(3,NUM(NY,NX),NY,NX)
  WatFLo2Litr(NY,NX)=WatFLo2Litr(NY,NX)+WatFLow2LitR(NY,NX)
  HeatFLo2LitrByWat(NY,NX)=HeatFLo2LitrByWat(NY,NX)+HeatFLoByWat2LitRi(NY,NX)

  HeatRadiation(NY,NX)=HeatRadiation(NY,NX)+Radnet2LitGrnd+Radnet2Snow
  HeatSensAir2Surf(NY,NX)=HeatSensAir2Surf(NY,NX)+HeatSensAir2Grnd+HeatSensAir2Snow
  HeatEvapAir2Surf(NY,NX)=HeatEvapAir2Surf(NY,NX)+LatentHeatEvapAir2Grnd+LatentHeatAir2Sno
  HeatSensVapAir2Surf(NY,NX)=HeatSensVapAir2Surf(NY,NX)+HeatSensVapAir2Soi+HeatSensEvap
  HeatNet2Surf(NY,NX)=HeatNet2Surf(NY,NX)+Radnet2LitGrnd+Radnet2Snow &
    +HeatSensAir2Grnd+HeatSensAir2Snow+LatentHeatEvapAir2Grnd &
    +LatentHeatAir2Sno+HeatSensVapAir2Soi+HeatSensEvap

  !EVAPG=negative evaporation from ground/top soil layer
  !EVAPR=evaporation from litter layer   
  !EVAPSN=evaporation from snow, sublimation+evaporation
  VapXAir2GSurf(NY,NX)=VapXAir2GSurf(NY,NX)+VapXAir2TopLay+VapXAir2LitR(NY,NX)+VapXAir2Sno(NY,NX)
  WaterFlow2MicPM(M,3,NUM(NY,NX),NY,NX)=WatXChange2WatTable(3,NUM(NY,NX),NY,NX)
  WaterFlow2MacPM(M,3,NUM(NY,NX),NY,NX)=ConvectWaterFlowMacP(3,NUM(NY,NX),NY,NX)
  end subroutine AccumWaterVaporHeatFluxes

!------------------------------------------------------------------------------------------

  subroutine InitSurfModel(M,NY,NX,RAR1,FKSAT)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),dimension(:,:),intent(in) :: RAR1
  real(r8),intent(out):: FKSAT
  integer :: L  
  real(r8) :: scalar,THETWT,HFLQR1,FLQRS
  real(r8) :: FLQRH,VOLAT0,ENGYD
  real(r8) :: ENGYB,RAS,TFND1,THETWA
  real(r8) :: HV

! begin_execution
! INITIALIZE NET SURFACE FLUX ACCUMULATORS
!
! TQR1,TQS1,TQW1,TQI1=net water and snowpack snow,water,ice runoff
! THQR1,THQS1=net convective heat from surface water and snow runoff
! BAREW,CVRDW=fractions of soil,litter cover including free water+ice
! RAGS= boundary layer resistance at soil surface
! PARG=boundary layer conductance above soil surface
!
  call ZeroSnowFlux(NY,NX)

  IF(VLHeatCapacity(0,NY,NX).GT.VHeatCapLitR(NY,NX))THEN
    BAREW(NY,NX)=AZMAX1(FracSurfAsBareSoi(NY,NX)-AMIN1(1.0_r8,AZMAX1(XVGeomLayer(NY,NX)/VOLWD(NY,NX))))
  ELSE
    BAREW(NY,NX)=1.0_r8
  ENDIF
  CVRDW(NY,NX)=1.0_r8-BAREW(NY,NX)
  RAGS(NY,NX)=1.0_r8/(BAREW(NY,NX)/RAGR(NY,NX)+CVRDW(NY,NX)/RAR1(NY,NX))
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
  IF(SoiBulkDensity(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    !soil bulk density significant
    FLQRS=AZMAX1(Prec2SoiMicP1(NY,NX)-VLairMicP1(NUM(NY,NX),NY,NX))
    FLQRH=AZMAX1(Prec2SoiMacP1(NY,NX)-VLairMacP1(NUM(NY,NX),NY,NX))
    HFLQR1=cpw*TairK(NY,NX)*(FLQRS+FLQRH)
    PrecAir2LitR=Prec2LitR1(NY,NX)+FLQRS+FLQRH
    PrecHeatAir2LitR=PrecHeat2LitR1(NY,NX)+HFLQR1
    PrecNet2SoiMicP=Prec2SoiMicP1(NY,NX)-FLQRS
    PrecNet2SoiMacP=Prec2SoiMacP1(NY,NX)-FLQRH
    PrecHeat2SoiNet=PrecHeat2SoiMicP1(NY,NX)-HFLQR1
  ELSE
    PrecAir2LitR=Prec2LitR1(NY,NX)
    PrecHeatAir2LitR=PrecHeat2LitR1(NY,NX)
    PrecNet2SoiMicP=Prec2SoiMicP1(NY,NX)
    PrecNet2SoiMacP=Prec2SoiMacP1(NY,NX)
    PrecHeat2SoiNet=PrecHeat2SoiMicP1(NY,NX)
  ENDIF
  Prec2LitR2=PrecAir2LitR*XNPR
  PrecHeat2LitR2=PrecHeatAir2LitR*XNPR
!
! WATER GAS EXCHANGE COEFFICIENTS IN SURFACE LITTER
!
! VOLA1,VLiceMicP1,VLWatMicP1,VsoiPM=total,ice-,water-,air-filled porosity
! TFND1=temperature effect on gas diffusivity
! DFGS=rate constant for air-water gas exchange
! Z1R,Z2RW,Z2RD,Z3RX=parameters for litter air-water gas transfers
! XNPD=time step for gas transfer calculations, it is tunable parameter
! TORT=tortuosity for aqueous diffusivity
! VOLAT0=ice-excluded porosity in litter

  VOLAT0=VLPoreLitR(NY,NX)-VLiceMicP1(0,NY,NX)
  IF(VOLAT0.GT.ZEROS2(NY,NX).AND.VLsoiAirPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    !litter layer is not saturated
    THETWA=AZMAX1(AMIN1(1.0_r8,VLWatMicP1(0,NY,NX)/VOLAT0))
    TFND1=TEFAQUDIF(TKSoi1(0,NY,NX))
    scalar=TFND1*XNPD
    DFGS(M,0,NY,NX)=fDFGS(scalar,THETWA,0.0_r8,is_litter=.true.)
  ELSE
    !litter layer saturated
    DFGS(M,0,NY,NX)=0.0_r8
  ENDIF
! VWatLitrX=surface litter water holding capacity, [m3 d-2]
  IF(VWatLitrX(NY,NX).GT.ZEROS(NY,NX))THEN
    !litter is able to hold water
    THETWT=AMIN1(1.0_r8,VLWatMicP(0,NY,NX)/VWatLitrX(NY,NX))
  ELSE
    THETWT=1.0
  ENDIF
  !TORT=tortuosity in litter (treated as micropore)
  TortMicPM(M,0,NY,NX)=TortMicporeW(THETWT)
!
! KINETIC ENERGY OF DIRECT RAINFALL AND THROUGHFALL
!
! PRECD,PRECB=direct,indirect precipn+irrign at soil surface
! ENGYD,ENGYB=energy impact of direct,indirect precipn+irrign at soil surface
! VOLWG=ground surface water retention capacity
! XVOLW=free surface water
! ZT=canopy height
! ENGYPM=total energy impact for use in erosion.f
! ENGYP=cumulative rainfall energy impact on soil surface
! FKSAT=reduction in soil surface Ksat from rainfall energy impact
! Note: A good reference for the following formula and alternatives
! is "Rainfall intensity-kinetic energy relationships for soil loss prediction",
! Kinnell, 1981
  IF(PRECD(NY,NX).GT.ZERO)THEN
    ENGYD=AZMAX1(8.95_r8+8.44_r8*LOG(PRECM(NY,NX)))
  ELSE
    ENGYD=0.0_r8
  ENDIF
  IF(PRECB(NY,NX).GT.ZERO)THEN
    ENGYB=AZMAX1(15.8_r8*SQRT(AMIN1(2.5_r8,GridMaxCanopyHeight(NY,NX)))-5.87_r8)
  ELSE
    ENGYB=0.0_r8
  ENDIF

  IF(ENGYD+ENGYB.GT.ZERO)THEN
    HV=1.0E+03_r8*AZMAX1(XVGeomLayer(NY,NX)-VOLWG(NY,NX))/AREA(3,NUM(NY,NX),NY,NX)
    ENGYPM(M,NY,NX)=(ENGYD*PRECD(NY,NX)+ENGYB*PRECB(NY,NX))*EXP(-2.0_r8*HV)*FracSurfAsBareSoi(NY,NX)*dts_HeatWatTP
    ENGYP(NY,NX)=ENGYP(NY,NX)+ENGYPM(M,NY,NX)
  ELSE
    ENGYPM(M,NY,NX)=0.0_r8
  ENDIF
  FKSAT=EXP(-2.0E-03_r8*(CSILT(NUM(NY,NX),NY,NX)+CCLAY(NUM(NY,NX),NY,NX))*ENGYP(NY,NX))

!
!  SNOWPACK FLUX ACCUMULATORS
!
   call InitSnowAccums(NY,NX)
!
! SURFACE FLUX ACCUMULATORS
!
! TWFLXL,TMLiceThawMacP=total freeze-thaw in micropores,macropores
! TLPhaseChangeHeat2Soi=total latent heat from freeze-thaw
! TFLWL,TFLWHL=net water flux in micropores,macropores
! THFLWL=net heat flux
!
!
! ENERGY EXCHANGE VARIABLES AT SNOW SURFACE IF PRESENT
!
! Radnet2Snow,LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,HeatNetFlx2Snow=netradn,latent,convective,sensible
! and storage heat fluxes
! CumWatFlx2SoiMacP=water from snowpack to soil micropores,macropores
! cumHeatFlowSno2Soi=conv heat from snowpack to soil micropores,macropores
!
! EVAPS,EVAPW=evaporation from soil,snowpack surfaces
! WatFlowSno2LitRM,WatFlowSno2MicPM,WatFlowSno2MacPM=water into litter,soil micropores,micropores for use in trnsfr.f
!
  VapXAir2Sno(NY,NX)=0._r8;EVAPS(NY,NX)=0.0_r8;EVAPW(NY,NX)=0.0_r8

  WatFlowSno2LitRM(M,NY,NX)=0.0_r8
  WatFlowSno2MicPM(M,NY,NX)=0.0_r8
  WatFlowSno2MacPM(M,NY,NX)=0.0_r8
!
  call PrepIterSnowLayer(M,NY,NX)
!
  end subroutine InitSurfModel  
!------------------------------------------------------------------------------------------

  subroutine SurfaceRunoff(M,N,NN,N1,N2,M4,M5,RCHQF,XN)
  implicit none
  integer, intent(in) :: M,N,NN,N1,N2,M4,M5
  real(r8), intent(in):: RCHQF,XN
  real(r8) :: ALT1,ALT2,DPTHW1,DPTHW2
  real(r8) :: VX
  !
  ! SURFACE BOUNDARY WATER FLUX
  !
  ! DPTHW1,DPTHW2=surface water depth of source,destination
  ! ALT1,ALT2=elevation of source,destination
  ! XVOLT=excess surface water+ice
  ! VOLWG=ground surface water retention capacity
  ! ExtWaterTable=natural water table depth
  ! QR1,HQR1=runoff, convective heat from runoff
  ! QR,HQR=hourly-accumulated runoff, convective heat from runoff
  ! QRM,QRV=runoff,velocity for erosion, solute transfer
  ! XN=direction
  !
  ! RUNOFF
  !
  DPTHW1=XVGeomLayer(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
  DPTHW2=VOLWG(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
  ALT1=ALTG(N2,N1)+DPTHW1
  ALT2=ALTG(N2,N1)+DPTHW2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)
  IF(ALT1.GT.ALT2.AND.CumDepth2LayerBottom(NU(N2,N1)-1,N2,N1)-DPTHW1.LT.ExtWaterTable(N2,N1))THEN
    QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HQR1(N,NN,M5,M4)=cpw*TKSoi1(0,N2,N1)*QR1(N,NN,M5,M4)
    QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
    HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
!
! RUNON
!
  ELSEIF(CumDepth2LayerBottom(NU(N2,N1)-1,N2,N1)-DPTHW1.GT.ExtWaterTable(N2,N1))THEN
    VX=AZMIN1((ExtWaterTable(N2,N1)-CumDepth2LayerBottom(NU(N2,N1)-1,N2,N1)+DPTHW1)*AREA(3,NUM(N2,N1),N2,N1))
    QRM(M,N2,N1)=VX*dts_wat
    QRV(M,N2,N1)=0.0_r8
    QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HQR1(N,NN,M5,M4)=cpw*TKSoi1(0,N2,N1)*QR1(N,NN,M5,M4)
    QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
    HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)

  ELSE
    QR1(N,NN,M5,M4)=0.0_r8
    HQR1(N,NN,M5,M4)=0.0_r8
  ENDIF
  QRMN(M,N,NN,M5,M4)=QR1(N,NN,M5,M4)
  IFLBM(M,N,NN,M5,M4)=0
  end subroutine SurfaceRunoff
!------------------------------------------------------------------------------------------

  subroutine PartitionPrecip(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: SnoFall,HeatByPrec
  real(r8) :: PrecThrufall2SoiMicP,PrecThrufall2SoiMacP,Rain4ToSno
  real(r8) :: PrecThrufall2LitR,PrecThrufall2Soil,PrecHeat2LitR,PrecHeat2Soil
  real(r8) :: PrecThruFall  
!     PRECA=precipitation+irrigation
!     PRECD,PRECB=direct,indirect precipn+irrign at soil surface
!     TFLWCI=net ice transfer to canopy, updated in hour1


  !convert water flux from m/hour to mm/hour
  PRECM(NY,NX)=1.0E+03_r8*PRECA(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  PRECD(NY,NX)=1.0E+03_r8*(PRECA(NY,NX)-TFLWCI(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  PRECB(NY,NX)=1.0E+03_r8*(TFLWCI(NY,NX)-PrecIntcptByCanG(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
!
!     RESIDUE WATER ABSORPTION CAPACITY
!
!     HCNDR=litter saturated hydraulic conductivity
!     DLYRR=litter depth
!
  HCNDR(NY,NX)=SatHydroCondLitR
  DLYRR(NY,NX)=AMAX1(2.5E-03_r8,DLYR(3,0,NY,NX))
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

  IF(PRECA(NY,NX).GT.0.0_r8.OR.SnoFalPrec(NY,NX).GT.0.0_r8)THEN
  ! there is precipitation
    Rain4ToSno=(PRECA(NY,NX)-PrecIntcptByCanG(NY,NX))*FracSurfAsSnow(NY,NX)
    SnoFall=SnoFalPrec(NY,NX)                                                        !snowfall
    HeatByPrec=cps*TairK(NY,NX)*SnoFall+cpw*TairK(NY,NX)*Rain4ToSno                  !incoming heat flux from precipitations to snow-covered surface
    PrecThruFall=(PRECA(NY,NX)-PrecIntcptByCanG(NY,NX))*FracSurfSnoFree(NY,NX)       !incoming precipitation to snow-free surface
    PrecThrufall2LitR=PrecThruFall*FracSurfByLitR(NY,NX)                             !water flux to snow-free coverd by litter
    PrecHeat2LitR=cpw*TairK(NY,NX)*PrecThrufall2LitR                                 !heat flux to snow-free surface covered by litter
    PrecThrufall2Soil=PrecThruFall*FracSurfAsBareSoi(NY,NX)                          !heat flux to snow-free surface not covered by litter
    PrecHeat2Soil=cpw*TairK(NY,NX)*PrecThrufall2Soil
    PrecThrufall2SoiMicP=PrecThrufall2Soil*SoilFracAsMicP(NUM(NY,NX),NY,NX)          !water flux to micropore
    PrecThrufall2SoiMacP=PrecThrufall2Soil*SoilFracAsMacP1(NUM(NY,NX),NY,NX)         !water flux to macropore
  ELSE
  ! no precipitation
    Rain4ToSno=-PrecIntcptByCanG(NY,NX)*FracSurfAsSnow(NY,NX)                   !
    SnoFall=0.0_r8
    HeatByPrec=cpw*TairK(NY,NX)*Rain4ToSno
    PrecThruFall=-PrecIntcptByCanG(NY,NX)*FracSurfSnoFree(NY,NX)
    PrecThrufall2LitR=PrecThruFall*FracSurfByLitR(NY,NX)
    PrecHeat2LitR=cpw*TairK(NY,NX)*PrecThrufall2LitR
    PrecThrufall2Soil=PrecThruFall*FracSurfAsBareSoi(NY,NX)
    PrecHeat2Soil=cpw*TairK(NY,NX)*PrecThrufall2Soil
    PrecThrufall2SoiMicP=PrecThrufall2Soil*SoilFracAsMicP(NUM(NY,NX),NY,NX)
    PrecThrufall2SoiMacP=PrecThrufall2Soil*SoilFracAsMacP1(NUM(NY,NX),NY,NX)
  ENDIF
!
!     PRECIP ON SNOW ARRAYS EXPORTED TO TRNSFR.F, TRNSFRS.F
!     FOR SOLUTE FLUX CALCULATIONS
!
!     SnoFalPrec,RainFalPrec,PRECQ,PRECI=snow,rain,snow+rain,irrigation
!     VHCPW,VLHeatCapSnowMN=current, minimum snowpack heat capacities
!     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!
  IF(SnoFalPrec(NY,NX).GT.0.0_r8.OR.(RainFalPrec(NY,NX).GT.0.0_r8 &
    .AND.VLHeatCapSnow(1,NY,NX).GT.VLHeatCapSnowMN(NY,NX)))THEN
    !there is precipitation, there is significant snow layer
    FLQRQ(NY,NX)=0.0_r8
    FLQRI(NY,NX)=0.0_r8
    FLQGQ(NY,NX)=PRECQ(NY,NX)
    FLQGI(NY,NX)=PRECI(NY,NX)
  ELSEIF((PRECQ(NY,NX).GT.0.0.OR.PRECI(NY,NX).GT.0.0_r8) &
    .AND.VLHeatCapSnow(1,NY,NX).LE.VLHeatCapSnowMN(NY,NX))THEN
    !there is insignificant snow layer
    FLQRQ(NY,NX)=PrecThrufall2LitR*PRECQ(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
    FLQRI(NY,NX)=PrecThrufall2LitR*PRECI(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
    FLQGQ(NY,NX)=PRECQ(NY,NX)-FLQRQ(NY,NX)
    FLQGI(NY,NX)=PRECI(NY,NX)-FLQRI(NY,NX)
  ELSE
    !no precipitation
    FLQRQ(NY,NX)=0.0_r8
    FLQRI(NY,NX)=0.0_r8
    FLQGQ(NY,NX)=0.0_r8
    FLQGI(NY,NX)=0.0_r8
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

  subroutine UpdateSurfaceAtM(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VOLIRZ,ENGYR,VLHeatCapLitRPre
  real(r8) :: TK0Prev,TVWatIceLitR,VWatLitrZ

  !update snow
  call UpdateSnowAtM(M,NY,NX)

  ! SURFACE RESIDUE WATER AND TEMPERATURE
  !
  ! XVOLT,XVOLW=free water+ice,water in litter layer
  ! VOLWM,VsoiPM=surface water,air content for use in trnsfr.f
  ! VWatLitrX=maximum water retention by litter
  ! VLHeatCapacity=volumetric heat capacity of litter
  ! VOLA1,VLWatMicP1,VLiceMicP1,VOLP1=pore,water,ice,air volumes of litter
  ! VWatLitrX=maximum water retention by litter
  ! LitrIceHeatFlxFrez,LitrIceFlxThaw=litter water,latent heat flux from freeze-thaw
  ! VLitR=dry litter volume
  ! THETWX,FracSoiPAsIce,FracSoiPAsAir=water,ice,air concentrations
  ! VLHeatCapacity=volumetric heat capacity of litter
  ! TK1=litter temperature
  ! HFLWRL,LitrIceHeatFlxFrez,THQR1=litter total cond+conv,latent,runoff heat flux

  VLWatMicP1(0,NY,NX)=AZMAX1(VLWatMicP1(0,NY,NX)+WatFLow2LitR(NY,NX)+LitrIceFlxThaw(NY,NX)+TQR1(NY,NX))
  VLiceMicP1(0,NY,NX)=AZMAX1(VLiceMicP1(0,NY,NX)-LitrIceFlxThaw(NY,NX)/DENSICE)
  VLairMicP1(0,NY,NX)=AZMAX1(VLPoreLitR(NY,NX)-VLWatMicP1(0,NY,NX)-VLiceMicP1(0,NY,NX))
  VLWatMicPM(M+1,0,NY,NX)=VLWatMicP1(0,NY,NX)
  VLsoiAirPM(M+1,0,NY,NX)=VLairMicP1(0,NY,NX)
  TVWatIceLitR=VLWatMicP1(0,NY,NX)+VLiceMicP1(0,NY,NX)
  XVGeomLayer(NY,NX)=AZMAX1(TVWatIceLitR-VWatLitrX(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ=VLWatMicP1(0,NY,NX)/TVWatIceLitR*VWatLitrX(NY,NX)
    VOLIRZ=VLiceMicP1(0,NY,NX)/TVWatIceLitR*VWatLitrX(NY,NX)
    XVLWatMicP(NY,NX)=AZMAX1(VLWatMicP1(0,NY,NX)-VWatLitrZ)
    XVLiceMicP(NY,NX)=AZMAX1(VLiceMicP1(0,NY,NX)-VOLIRZ)
  ELSE
    XVLWatMicP(NY,NX)=0.0_r8
    XVLiceMicP(NY,NX)=0.0_r8
  ENDIF
  XVOLTM(M+1,NY,NX)=XVGeomLayer(NY,NX)
  XVLWatMicPM(M+1,NY,NX)=XVLWatMicP(NY,NX)
  XVLiceMicPM(M+1,NY,NX)=XVLiceMicP(NY,NX)
  IF(VLitR(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat(0,NY,NX)=AZMAX1t(VLWatMicP1(0,NY,NX)/VLitR(NY,NX))
    FracSoiPAsIce(0,NY,NX)=AZMAX1t(VLiceMicP1(0,NY,NX)/VLitR(NY,NX))
    FracSoiPAsAir(0,NY,NX)=AZMAX1t(VLairMicP1(0,NY,NX)/VLitR(NY,NX))*AZMAX1t((1.0_r8-XVGeomLayer(NY,NX)/VOLWD(NY,NX)))
  ELSE
    FracSoiPAsWat(0,NY,NX)=0.0_r8
    FracSoiPAsIce(0,NY,NX)=0.0_r8
    FracSoiPAsAir(0,NY,NX)=1.0_r8
  ENDIF
  THETPM(M+1,0,NY,NX)=FracSoiPAsAir(0,NY,NX)
  VLHeatCapLitRPre=VLHeatCapacity(0,NY,NX)              !heat capacity
  TK0Prev=TKSoi1(0,NY,NX)                                 !residual temperature
  ENGYR=VLHeatCapacity(0,NY,NX)*TKSoi1(0,NY,NX)  !initial energy content
  VLHeatCapacity(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VLWatMicP1(0,NY,NX)+cpi*VLiceMicP1(0,NY,NX)  !update heat capacity

  IF(VLHeatCapacity(0,NY,NX).GT.VHeatCapLitR(NY,NX))THEN
    TKSoi1(0,NY,NX)=(ENGYR+HeatFLoByWat2LitRi(NY,NX)+LitrIceHeatFlxFrez(NY,NX)+THQR1(NY,NX))/VLHeatCapacity(0,NY,NX)
!    write(*,*)'UpdateSurfaceAtM=',M,NY,NX,TKSoi1(0,NY,NX),TK0Prev,HeatFLoByWat2LitRi(NY,NX)
!    write(*,*)ENGYR/VLHeatCapacity(0,NY,NX),HeatFLoByWat2LitRi(NY,NX)/VLHeatCapacity(0,NY,NX),&
!      LitrIceHeatFlxFrez(NY,NX)/VLHeatCapacity(0,NY,NX),THQR1(NY,NX)/VLHeatCapacity(0,NY,NX)
!    IF(ABS(VLHeatCapacity(0,NY,NX)/VLHeatCapLitRPre-1._r8)>0.025_r8.or. &
!      abs(TKSoi1(0,NY,NX)/TK0Prev-1._r8)>0.025_r8)THEN
!      TKSoi1(0,NY,NX)=TKSoi1(NUM(NY,NX),NY,NX)
!    ENDIF
  ELSE
    TKSoi1(0,NY,NX)=TKSoi1(NUM(NY,NX),NY,NX)
  ENDIF
  end subroutine UpdateSurfaceAtM
!------------------------------------------------------------------------------------------

  subroutine SumAftEnergyBalance(NY,NX,LWRadGrnd,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,TotHeatAir2LitR)
  implicit none
  integer, intent(in)  :: NY,NX
  real(r8), intent(in) :: LWRadGrnd
  real(r8), intent(in) :: VapXAir2TopLay
  real(r8), intent(in) :: HeatSensLitR2Soi1
  real(r8), intent(in) :: HeatSensVapLitR2Soi1
  real(r8), intent(in) :: EvapLitR2Soi1
  real(r8), intent(in) :: TotHeatAir2LitR
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

  NetWatFlx2SoiMicP=PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1
  NetWatFlx2SoiMacP=PrecNet2SoiMacP
  cumHeatSensAir2Soil=PrecHeat2SoiNet+HeatFluxAir2Soi+HeatSensVapLitR2Soi1+HeatSensLitR2Soi1
  NetWatFlx2LitR=PrecAir2LitR+VapXAir2LitR(NY,NX)-EvapLitR2Soi1
  CumHeatSensAir2LitR=PrecHeatAir2LitR+TotHeatAir2LitR-HeatSensVapLitR2Soi1-HeatSensLitR2Soi1
!  write(*,*)'CumHeatSensAir2LitR ',PrecHeatAir2LitR,TotHeatAir2LitR,HeatSensVapLitR2Soi1,HeatSensLitR2Soi1
  FLWVLS=(VLWatMicP1(NUM(NY,NX),NY,NX)-VLWatMicPX1(NUM(NY,NX),NY,NX))*dts_HeatWatTP
!
! GENERATE NEW SNOWPACK
!
! XFLWS,WatXfer2SnoLay,IceXfer2SnoLay=hourly snow,water,ice transfer
! SnowFallt,Rain2Snowt,Ice2Snowt=snow,water,ice input to snowpack
! HeatXfer2SnoLay=hourly convective heat flux from snow,water,ice transfer
! HeatFall2Snowt=convective heat flux from snow,water,ice to snowpack
!
  IF(VLHeatCapSnow(1,NY,NX).LE.VLHeatCapSnowMN(NY,NX).AND.SnowFallt(NY,NX).GT.ZEROS(NY,NX))THEN
    SnoXfer2SnoLay(1,NY,NX)=SnoXfer2SnoLay(1,NY,NX)+SnowFallt(NY,NX)
    WatXfer2SnoLay(1,NY,NX)=WatXfer2SnoLay(1,NY,NX)+Rain2Snowt(NY,NX)
    IceXfer2SnoLay(1,NY,NX)=IceXfer2SnoLay(1,NY,NX)+Ice2Snowt(NY,NX)
    HeatXfer2SnoLay(1,NY,NX)=HeatXfer2SnoLay(1,NY,NX)+HeatFall2Snowt(NY,NX)
!     WRITE(*,4422)'INIT',I,J,SnowFallt(NY,NX),Rain2Snowt(NY,NX)
!    3,Ice2Snowt(NY,NX),HeatFall2Snowt(NY,NX),SnoXfer2SnoLay(1,NY,NX),WatXfer2SnoLay(1,NY,NX)
!    2,IceXfer2SnoLay(1,NY,NX),HeatXfer2SnoLay(1,NY,NX),HeatFlow2Soili(3,NUM(NY,NX),NY,NX)
!    3,HeatFlow2Soil(3,NUM(NY,NX),NY,NX),FracSurfSnoFree(NY,NX),VLHeatCapacity(NUM(NY,NX),NY,NX)
!    4*TKSoi1(NUM(NY,NX),NY,NX),HeatFLoByWat2LitRi(NY,NX),HeatFLo2LitrByWat(NY,NX)
!    5,VLHeatCapacity(0,NY,NX)*TKSoi1(0,NY,NX),HeatNet2Surf(NY,NX),Radnet2LitGrnd,Radnet2Snow
!    2,HeatSensAir2Grnd,HeatSensAir2Snow,LatentHeatEvapAir2Grnd
!      ,LatentHeatAir2Sno,HeatSensVapAir2Soi,HeatSensEvap
  ENDIF
  !LWRadBySurf=longwave emission from litter and surface soil
  LWRadBySurf(NY,NX)=LWRadBySurf(NY,NX)+LWRadGrnd
  end subroutine SumAftEnergyBalance
!------------------------------------------------------------------------------------------
  subroutine SurfacePhysModel(M,NHE,NHW,NVS,NVN,RAR1,FKSAT,HeatFlux2Ground,TopLayWatVol)
  implicit none
  integer, intent(in) :: M,NHE,NHW,NVS,NVN
  real(r8), dimension(:,:),intent(inout) :: RAR1
  REAL(R8), dimension(:,:),INTENT(OUT) :: FKSAT
  real(r8), dimension(:,:),intent(out) :: HeatFlux2Ground
  real(r8),dimension(:,:),intent(inout) :: TopLayWatVol
  real(r8) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvap,VapXAir2TopLay
  integer :: N1,N2,NX,NY


  D9895: DO  NX=NHW,NHE
    D9890: DO  NY=NVN,NVS

!      write(*,*)'SurfacePhysModel',NY,NX,'M=',M,TKS(0,NY,NX)

      call SurfaceEnergyModel(M,NX,NY,RAR1,FKSAT(NY,NX),HeatFlux2Ground(NY,NX),LatentHeatAir2Sno,&
        HeatSensEvap,HeatSensAir2Snow,Radnet2Snow,TopLayWatVol,VapXAir2TopLay)

!      write(*,*)'TXKR SurfaceEnergyModel MM=',M,TKSoi1(0,NY,NX)

    ! CAPILLARY EXCHANGE OF WATER BETWEEN SOIL SURFACE AND RESIDUE
      call SurfLitrSoilWaterExchange(M,NY,NX,FKSAT(NY,NX))

      call InfilSRFRoffPartition(M,NY,NX,N1,N2)
    !
      call LateralHydroExchange(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
    !
      call AccumWaterVaporHeatFluxes(M,NY,NX,LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,&
        Radnet2Snow,VapXAir2TopLay)
!      write(*,*)'end SurfacePhysModel'
    ENDDO D9890
  ENDDO D9895

  end subroutine SurfacePhysModel

!------------------------------------------------------------------------------------------
  subroutine SurfaceEnergyModel(M,NX,NY,RAR1,FKSAT,HeatFlux2Ground1,LatentHeatAir2Sno,&
    HeatSensEvap,HeatSensAir2Snow,Radnet2Snow,TopLayWatVol,VapXAir2TopLay)
  implicit none
  integer, intent(in) :: M,NX,NY
  real(r8), dimension(:,:),intent(inout) :: RAR1
  REAL(R8),INTENT(OUT) :: FKSAT,HeatFlux2Ground1
  real(r8), intent(out) :: Radnet2Snow,LatentHeatAir2Sno,HeatSensAir2Snow,HeatSensEvap
  real(r8),dimension(:,:),intent(inout) :: TopLayWatVol
  real(r8), intent(out) :: VapXAir2TopLay
  integer :: N1,N2

  !RAR1 is input
  call InitSurfModel(M,NY,NX,RAR1,FKSAT)

! updates RAR1
  call AtmLandSurfExchange(M,NY,NX,RAR1,TopLayWatVol,LatentHeatAir2Sno,HeatSensEvap,HeatSensAir2Snow,&
    Radnet2Snow,VapXAir2TopLay)

  HeatFlux2Ground1=HeatFluxAir2Soi
  
  end subroutine SurfaceEnergyModel



end module SurfPhysMod