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
  use EcoSIMSolverPar
  use SoilPhysParaMod
implicit none
  private
  character(len=*), parameter :: mod_filename=__FILE__

  !surface model
  public :: StageSurfStateVars
  public :: SurfacePhysModel
  public :: SurfaceEnergyModel
  !
  public :: SurfaceRunoff
  public :: UpdateSurfaceAtM
  public :: CalcThermConduc


! HCNDRR=saturated hydraulic conductivity of surface litter
! FENGYP=rate constant for restoring surface Ksat
! RARX=minimum boundary layer resistances of litter (h m-1)
! EMMS,EMMW,EMMR=emissivities of surface soil, snow and litter
! RACX,RARX=minimum boundary layer resistances of canopy,litter (h m-1)

  real(r8), parameter :: FENGYP=1.0E-03_r8
  real(r8), parameter :: HCNDRR=25.0_r8
  real(r8), parameter :: RARX=0.0139_r8
  real(r8), parameter :: EMMS=0.97_r8    !soil emissivity
  real(r8), parameter :: EMMW=0.97_r8    !snowpack emissivity
  real(r8), parameter :: EMMR=0.97_r8    !surfce litter emissivity
  real(r8), parameter :: RACX=0.0139_r8    !total canopy boundary later resistance h/m  

  real(r8) :: VOLW12,VHCPR2,VHCP12,VFLXR
  real(r8) :: FLV1,TKS1,TKR1,THRMZ,SFLXR
  real(r8) :: HFLXR,HWFLM2,RFLXR,FLYM2  
  real(r8) :: EFLXR,HWFLV1,HFLCR1,VFLXG
  real(r8) :: SFLXG,RFLXG,HeatFlux2Ground
  real(r8) :: EFLXG,VFLXW,FLWHLG,SFLXW
  real(r8) :: RFLXW,HFLXW,HFLWRLW,HFLWRLG
  real(r8) :: HFLWLG,HFLWLW,FLWRLG,FLWRLW
  real(r8) :: FLWLG,FLWLXG,FLWLXW,FLWLW
  real(r8) :: FLWHLW,EFLXW,HWFLQM,FLQM,FLYM
  real(r8) :: HWFLYM,FLHM  
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
  real(r8) :: VOLWRZ,TVOLWI,VOLIRZ
!
! SET INITIAL SOIL VALUES
!
! WFLXR,TFLXR=initialize surface litter freeze,thaw,latent heat
  WFLXR(NY,NX)=0.0_r8
  TFLXR(NY,NX)=0.0_r8

!
! ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
!     FOR USE AT INTERNAL TIME STEP IN SURFACE LITTER
!
!     THRMG=longwave emission from litter surface
!     VHCP1=volumetric heat capacity of litter
!     VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of litter
!     VOLWRX=maximum water retention by litter
!     XVOLT,XVOLW=free surface water+ice,water
!     VHCPRX=min heat capacity for litter water,heat fluxes
!     VOLR=litter volume
!     THETW*,THETI*,THETP*=water,ice,air concentrations
!     PSISM*=litter matric water potential
!
  THRMG(NY,NX)=0.0_r8
  VHCP1(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW(0,NY,NX)+cpi*VOLI(0,NY,NX)
  VOLA10(NY,NX)=VOLA(0,NY,NX)
  VOLW1(0,NY,NX)=AZMAX1(VOLW(0,NY,NX))
  VOLI1(0,NY,NX)=AZMAX1(VOLI(0,NY,NX))
  VOLP1(0,NY,NX)=AZMAX1(VOLA10(NY,NX)-VOLW1(0,NY,NX)-VOLI1(0,NY,NX))
  VOLWM(1,0,NY,NX)=VOLW1(0,NY,NX)
  VOLPM(1,0,NY,NX)=VOLP1(0,NY,NX)
  TVOLWI=VOLW1(0,NY,NX)+VOLI1(0,NY,NX)
  XVOLT(NY,NX)=AZMAX1(TVOLWI-VOLWRX(NY,NX))
  IF(TVOLWI.GT.ZEROS(NY,NX))THEN
    VOLWRZ=VOLW1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
    VOLIRZ=VOLI1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
    XVOLW(NY,NX)=AZMAX1(VOLW1(0,NY,NX)-VOLWRZ)
    XVOLI(NY,NX)=AZMAX1(VOLI1(0,NY,NX)-VOLIRZ)
  ELSE
    XVOLW(NY,NX)=0.0_r8
    XVOLI(NY,NX)=0.0_r8
  ENDIF
  XVOLTM(1,NY,NX)=XVOLT(NY,NX)
  XVOLWM(1,NY,NX)=XVOLW(NY,NX)
  XVOLIM(1,NY,NX)=XVOLI(NY,NX)
  IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
    THETWX(0,NY,NX)=AZMAX1t(VOLW1(0,NY,NX)/VOLR(NY,NX))
    THETIX(0,NY,NX)=AZMAX1t(VOLI1(0,NY,NX)/VOLR(NY,NX))
    THETPX(0,NY,NX)=AZMAX1t(VOLP1(0,NY,NX)/VOLR(NY,NX))*AZMAX1t((1.0_r8-XVOLT(NY,NX)/VOLWD(NY,NX)))
  ELSE
    THETWX(0,NY,NX)=0.0_r8
    THETIX(0,NY,NX)=0.0_r8
    THETPX(0,NY,NX)=1.0
  ENDIF
  THETPM(1,0,NY,NX)=THETPX(0,NY,NX)
  PSISM1(0,NY,NX)=PSISM(0,NY,NX)
  TK1(0,NY,NX)=TKS(0,NY,NX)

  end subroutine CopySurfaceVars

!------------------------------------------------------------------------------------------

  subroutine PartionSurfaceFraction(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

!     SNOW AND RESIDUE COVERAGE OF SOIL SURFACE
!     FSNW,FSNX=fractions of snow,snow-free cover
!     DPTHS=snowpack depth
!     DPTHSX=minimum snowpack depth for full cover
!     BARE,CVRD=fractions of soil,litter cover
  FSNW(NY,NX)=AMIN1(1.0_r8,SQRT((DPTHS(NY,NX)/DPTHSX)))
  FSNX(NY,NX)=1.0_r8-FSNW(NY,NX)
  !if there is heat-wise significant litter layer
  IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
    BARE(NY,NX)=AMIN1(1.0_r8,AZMAX1(EXP(-0.8E-02_r8*(ORGC(0,NY,NX)/AREA(3,0,NY,NX)))))
  ELSE
    BARE(NY,NX)=1.0_r8
  ENDIF

  CVRD(NY,NX)=1.0_r8-BARE(NY,NX)
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

  VPQ(NY,NX)=VPA(NY,NX)-TLEX(NY,NX)/(VAP*AREA(3,NUM(NY,NX),NY,NX))
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
!     THRYW,THRYG,THRYR=longwave radn incident at snowpack,soil,litter
!     THRMW,THRMS,THRMR=longwave radn emitted by snowpack,soil,litter
!     EMMW,EMMS,EMMR=emissivity of snowpack,soil,litter surfaces
!     THS=sky longwave radiation
!     THRMCX=longwave radiation emitted by canopy

  RADGX=RADG(NY,NX)*XNPH
  RADXW(NY,NX)=RADGX*FSNW(NY,NX)*XNPS
  RADXG(NY,NX)=RADGX*FSNX(NY,NX)*BARE(NY,NX)      
  RADXR(NY,NX)=RADGX*FSNX(NY,NX)*CVRD(NY,NX)*XNPR

  THRYX=(THS(NY,NX)*FRADG(NY,NX)+THRMCX(NY,NX))*XNPH
  THRYW(NY,NX)=THRYX*FSNW(NY,NX)*XNPS
  THRYG(NY,NX)=THRYX*FSNX(NY,NX)*BARE(NY,NX)
  THRYR(NY,NX)=THRYX*FSNX(NY,NX)*CVRD(NY,NX)*XNPR
  ! EMMS,EMMW,EMMR=emissivities of surface soil, snow and litter
  !what is 2.04E-10_r8 it is stefan-boltzman constant converted into MJ m-2 K-4/per hour
  THRMW(NY,NX)=EMMW*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPY
  THRMS(NY,NX)=EMMS*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*BARE(NY,NX)*XNPH
  THRMR(NY,NX)=EMMR*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*CVRD(NY,NX)*XNPZ
!
  end subroutine SurfaceRadiation
!------------------------------------------------------------------------------------------
  subroutine SurfaceResistances(NY,NX,RAR1)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), dimension(:,:), intent(out):: RAR1
  real(r8) :: THETPX0,DFVR  
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
!     RAB,RAC=isothermal blr above canopy, canopy blr
!     ZT,ZS=canopy, surface roughness heights
!     UA,UAG=windspeeds above,below canopy
!     VPQ,VPA=vapor pressure within,above canopy
!     TKQ,TairK=temperature within,above canopy
!     TLEX,TSHX=net latent,sensible heat fluxes x blrs from prev hour
!     VAP=latent heat of evaporation
!     1.25E-03=heat capacity of air
!     AREA=surface area of grid cell
!
  ALFZ=2.0_r8*(1.0_r8-FRADG(NY,NX))
  IF(RAB(NY,NX).GT.ZERO.AND.GridMaxCanopyHeight(NY,NX).GT.ZS(NY,NX).AND.ALFZ.GT.ZERO)THEN
    RAC(NY,NX)=AMIN1(RACX,AZMAX1(GridMaxCanopyHeight(NY,NX)*EXP(ALFZ) &
      /(ALFZ/RAB(NY,NX))*AZMAX1(EXP(-ALFZ*ZS(NY,NX)/GridMaxCanopyHeight(NY,NX)) &
      -EXP(-ALFZ*(ZD(NY,NX)+ZR(NY,NX))/GridMaxCanopyHeight(NY,NX)))))
    UAG=UA(NY,NX)*EXP(-ALFZ)
  ELSE
    RAC(NY,NX)=0.0_r8
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
!     THETPX*=air-filled porosity of litter
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
  RAG(NY,NX)=RAC(NY,NX)+RAB(NY,NX)
  RAGW(NY,NX)=RAG(NY,NX)
  RAGR(NY,NX)=RAG(NY,NX)+RARX
  RARG(NY,NX)=RAGR(NY,NX)
  THETPX0=AMAX1(ZERO2,THETPX(0,NY,NX))
  DFVR=THETPX0*POROQ*THETPX0/POROS(0,NY,NX)
  RAR1(NY,NX)=RAG(NY,NX)+RAR(NY,NX)/DFVR
  PAREX=AREA(3,NUM(NY,NX),NY,NX)*XNPH               !conductance for latent heat flux
  PARSX=1.25E-03_r8*AREA(3,NUM(NY,NX),NY,NX)*XNPH   !conductance for sensible heat flux
  
  PAREW(NY,NX)=PAREX*FSNW(NY,NX)*XNPS
  PARSW(NY,NX)=PARSX*FSNW(NY,NX)*XNPS
  PAREG(NY,NX)=PAREX*FSNX(NY,NX)
  PARER(NY,NX)=PAREX*FSNX(NY,NX)*XNPR*CVRD(NY,NX)
  PARSG(NY,NX)=PARSX*FSNX(NY,NX)
  PARSR(NY,NX)=PARSX*FSNX(NY,NX)*XNPR*CVRD(NY,NX)

!     PARR=boundary layer conductance above litter,soil surfaces
!
  RAS=SnowBNDResistance(NY,NX)

  PARR(NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*XNPH/(RAGR(NY,NX)+RAS)   !this includes snow layer resistance
  end subroutine SurfaceResistances
!------------------------------------------------------------------------------------------

  subroutine SurfLitterIterate(M,NY,NX,VOLWR2,ATCNDR,ATCNVR,PSISV1,RFLX0)

  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: ATCNDR,ATCNVR,PSISV1,RFLX0
  real(r8),intent(inout) :: VOLWR2
  integer  :: NN
  real(r8) :: tk1pre,VHCPRXX
  real(r8) :: RAa,VP1,VPR,VPY
  real(r8) :: TKXR,TK1X,TKY
  real(r8) :: FLVC,FLVX,ENGYR
  real(r8) :: HFLWC,HFLWX,EVAPR2,EFLXR2,VFLXR2
  real(r8) :: HFLX02,THRMZ2,THETWR,SFLXR2,FLV2,HWFLV2
  real(r8) :: PARE,PARS,RAGX,RI,RFLXR2,HFLXR2,HFLCR2

! begin_execution

  D5000: DO NN=1,NPR
    IF(VHCPR2.GT.VHCPRX(NY,NX))THEN
      !
      ! AERODYNAMIC RESISTANCE ABOVE RESIDUE INCLUDING
      ! RESISTANCE IMPOSED BY PLANT CANOPY
      !
      ! RI=Richardsons number
      ! RIB=isothermal RI
      ! TKQ,TKR1=canopy air,litter temperature
      ! RZ=surface resistance to evaporation, given as a prescribed parameter
      ! RAGX,RA=litter blr
      ! RAG,RAGR=isothermal blr at ground surface
      ! PARE,PARS=blcs for litter latent,sensible heat fluxes
      !
      RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB(NY,NX)*(TKQ(NY,NX)-TKR1)))
      RAGX=AMAX1(RAM,0.8_r8*RAGR(NY,NX),AMIN1(1.2_r8*RAGR(NY,NX),RARG(NY,NX)/(1.0_r8-10.0_r8*RI)))
      RAGR(NY,NX)=RAGX
      RAa=RAGX
      PARE=PARER(NY,NX)/(RAa+RZ)
      PARS=PARSR(NY,NX)/RAa
!
!     NET RADIATION AT RESIDUE SURFACE
!
!     THRMZ2=longwave radiation emitted by litter
!     RFLXR2=litter net radiation,after taking out outgoing radiation from litter 
!     THETWR=litter water content
!     VOLWRX=maximum water retention by litter
!     PSISM1=litter matric water potential
!
      THRMZ2=THRMR(NY,NX)*TKR1**4._r8
      RFLXR2=RFLX0-THRMZ2
      IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
        THETWR=AMIN1(VOLWRX(NY,NX),VOLWR2)/VOLR(NY,NX)
      ELSE
        THETWR=POROS0(NY,NX)
      ENDIF

      IF(VOLR(NY,NX).GT.ZEROS(NY,NX).AND.VOLW1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
        THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
        IF(THETWR.LT.FC(0,NY,NX))THEN
          PSISM1(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
            +((FCL(0,NY,NX)-LOG(THETWR))/FCD(0,NY,NX)*PSIMD(NY,NX))))
        ELSEIF(THETWR.LT.POROS0(NY,NX))THEN
          PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX)+(((PSL(0,NY,NX)-LOG(THETWR)) &
            /PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
        ELSE
          THETWR=POROS0(NY,NX)
          PSISM1(0,NY,NX)=PSISE(0,NY,NX)
        ENDIF
      ELSE
        THETWR=POROS0(NY,NX)
        PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
!
!     VAPOR FLUX AT RESIDUE SURFACE
!
!     VPR,VP1,VPQ=vapor pressure in litter,soil,canopy air, m3/m-3
!     TKS1=soil temperature
!     EVAPR2=litter evaporation
!     EFLXR2=litter latent heat flux
!     VAP=latent heat of evaporation
!     VFLXR2=convective heat of evaporation flux
!
      VPR=vapsat(TKR1)*EXP(18.0_r8*PSISM1(0,NY,NX)/(RGAS*TKR1))   !in litter
      if(abs(VPR)>1.e20_r8)then
        write(*,*)'TKR1=',TKR1,TK1(0,NY,NX),TK1(NUM(NY,NX),NY,NX)
        write(*,*)'PSISM1(0,NY,NX)=',PSISM1(0,NY,NX)
        call endrun(trim(mod_filename)//'at line',__LINE__)
      endif

      VP1=vapsat(TKS1)*EXP(18.0_r8*PSISV1/(RGAS*TKS1))    !in soil
      EVAPR2=AMAX1(-AZMAX1(VOLWR2*XNPX),PARE*(VPQ(NY,NX)-VPR))

      EFLXR2=EVAPR2*VAP             !energy flux
      VFLXR2=EVAPR2*cpw*TKR1        !energy flux
      !
      ! SOLVE FOR RESIDUE TO SOIL SURFACE HEAT FLUXES
      !
      ! FLVC,FLVX=vapor unconstrained,vapor constrained vapor flux
      ! XNPZ=time step for litter flux calculations
      ! VPY=equilibrium vapor concentration
      ! VOLPM=litter,soil air filled porosity
      ! FLV2=litter soil vapor flux
      ! HWFLV2=convective heat of litter soil vapor flux
      ! TKXR,TK1X=interim calculation of litter,soil temperatures
      ! TKY=equilibrium litter-soil temperature
      ! HFLWC,HFLWX=litter-soil heat flux unlimited,limited by heat
      ! HFLCR2=litter-soil heat flux
      ! THETPM: air-filled porosity
      IF(THETPM(M,0,NY,NX).GT.THETX.AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
        FLVC=ATCNVR*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*CVRD(NY,NX)*XNPZ
        if(abs(FLVC)>1.e20_r8)then
          write(*,*)'ATCNVR*(VPR-VP1)=',ATCNVR,VPR,VP1
          write(*,*)'FSNX(NY,NX),CVRD(NY,NX)=',FSNX(NY,NX),CVRD(NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
        VPY=(VPR*VOLPM(M,0,NY,NX)+VP1*VOLPM(M,NUM(NY,NX),NY,NX)) &
          /(VOLPM(M,0,NY,NX)+VOLPM(M,NUM(NY,NX),NY,NX))
        FLVX=(VPR-VPY)*VOLPM(M,0,NY,NX)*XNPB
        if(abs(FLVX)>1.e20_r8)then
          write(*,*)'(VPR-VPY)*VOLPM(M,0,NY,NX)',VPR,VPY,VOLPM(M,0,NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
        IF(FLVC.GE.0.0_r8)THEN
          !from litter to soil
          FLV2=AZMAX1(AMIN1(FLVC,FLVX))
          if(abs(FLV2)>1.e20_r8)then
            write(*,*)'FLVC,FLVX',FLVC,FLVX
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          HWFLV2=(cpw*TKR1+VAP)*FLV2
        ELSE
          !from soil to litter
          FLV2=AZMIN1(AMAX1(FLVC,FLVX))
          if(abs(FLV2)>1.e20_r8)then
            write(*,*)'FLVC,FLVX',FLVC,FLVX
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          HWFLV2=(cpw*TKS1+VAP)*FLV2
        ENDIF
      ELSE
        FLV2=0.0_r8
        HWFLV2=0.0_r8
      ENDIF
      TKXR=TKR1-HWFLV2/VHCPR2         !update litter layer temperature
      TK1X=TKS1+HWFLV2/VHCP12         !update soil layer temperature
      TKY=(TKXR*VHCPR2+TK1X*VHCP12)/(VHCPR2+VHCP12)   !equilibrium temperature
      HFLWX=(TKXR-TKY)*VHCPR2*XNPB    !sensible heat flux
      HFLWC=ATCNDR*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*CVRD(NY,NX)*XNPZ
      IF(HFLWC.GE.0.0_r8)THEN
        HFLCR2=AZMAX1(AMIN1(HFLWX,HFLWC))
      ELSE
        HFLCR2=AZMIN1(AMAX1(HFLWX,HFLWC))
      ENDIF
!
!     SOLVE FOR RESIDUE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
!     SFLXR2,RFLXR2,EFLXR2=litter sensible,net radn,latent heat fluxes
!     HFLX02,HFLXR2=storage,total litter heat flux
!
      SFLXR2=PARS*(TKQ(NY,NX)-TKR1)    !sensible heat flux between canopy air and litter surface
      HFLX02=RFLXR2+EFLXR2+SFLXR2      !
      HFLXR2=HFLX02+VFLXR2
!
!     AGGREGATE WATER AND ENERGY FLUXES FROM TIME STEP FOR LITTER
!     CALCULATIONS TO THAT FOR SOIL PROFILE
!
      EVAPR(NY,NX)=EVAPR(NY,NX)+EVAPR2
      RFLXR=RFLXR+RFLXR2
      EFLXR=EFLXR+EFLXR2
      VFLXR=VFLXR+VFLXR2
      SFLXR=SFLXR+SFLXR2
      HFLXR=HFLXR+HFLXR2
      ! write(*,*)'.........HFLXR',HFLXR,HFLXR2,TKQ(NY,NX),TairK(NY,NX),TKR1
      ! write(*,*)'HFLX02+VFLXR2,NN',HFLX02,VFLXR2,NN
      ! write(*,*)'RFLXR2+EFLXR2+SFLXR2',RFLXR2,EFLXR2,SFLXR2
      ! write(*,*)'PARS*(TKQ(NY,NX)-TKR1)',PARS,TKQ(NY,NX),TKR1
      ! call print_info(trim(mod_filename)//' at line',__LINE__)
      FLV1=FLV1+FLV2
      HWFLV1=HWFLV1+HWFLV2
      HFLCR1=HFLCR1+HFLCR2
      THRMZ=THRMZ+THRMZ2
    ELSE
      !not significant litter layer heat capacity
      EVAPR2=0.0_r8
      RFLXR2=0.0_r8
      EFLXR2=0.0_r8
      VFLXR2=0.0_r8
      SFLXR2=0.0_r8
      HFLXR2=0.0_r8
      FLV2=0.0_r8
      HWFLV2=0.0_r8
      HFLCR2=0.0_r8
      THRMZ2=0.0_r8
    ENDIF
    !VOLWR2,VOLW12=water in litter, topsoil layer 
    VOLWR2=VOLWR2+FLYM2+EVAPR2-FLV2
    VOLW12=VOLW12+FLV2
    ENGYR=VHCPR2*TKR1
    VHCPRXX=VHCPR2
    
    ! VHCPR2: heat capacity, kJ/kg/Kelvin
    VHCPR2=cpo*ORGC(0,NY,NX)+cpw*VOLWR2+cpi*VOLI1(0,NY,NX)
    VHCP12=VHCP12+cpw*FLV2
    tk1pre=TKR1
    TKR1=(ENGYR+HFLXR2+HWFLM2-HWFLV2-HFLCR2)/VHCPR2
    IF(ABS(VHCPRXX/VHCPR2-1._r8)>0.025_r8.or.abs(TKR1/tk1pre-1._r8)>0.025_r8)then
      TKR1=TK1(0,NY,NX)
    endif
    TKS1=TKS1+(HWFLV2+HFLCR2)/VHCP12

  ENDDO D5000
  end subroutine SurfLitterIterate
!------------------------------------------------------------------------------------------
  subroutine SRFLitterEnergyBalance(M,NY,NX,PSISV1)
!!
! Description:
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(in) :: PSISV1
  real(r8) :: VOLWR2,ATCNDR,ATCNVR,RFLX0
  real(r8) :: WTHET0,CNVR,CNV1
  real(r8) :: THETRR,TCNDR,XNUSW0,XNUSA0
  real(r8) :: TCNDA0,RYLXW0,RYLNW0
  real(r8) :: TCNDW0,TCND1,RYLXA0
  real(r8) :: RYLNA0,DTKX,DTHW0,DTHA0
  real(r8) :: ALBL  
! begin_execution
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!

  EVAPR(NY,NX)=0.0_r8
  RFLXR=0.0_r8
  EFLXR=0.0_r8
  VFLXR=0.0_r8
  SFLXR=0.0_r8
  HFLXR=0.0_r8
  FLV1=0.0_r8
  HWFLV1=0.0_r8
  HFLCR1=0.0_r8
  THRMZ=0.0_r8
  IF(VHCP1(0,NY,NX).LE.VHCPRX(NY,NX))THEN
    TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
    RETURN
  ENDIF  
!
! NET RADIATION AT RESIDUE SURFACE
!
! ALBL=litter albedo
! BKVL=litter mass
! VOLW1,VOLI1=water,ice volume in litter
! RADXR,THRYR=incoming shortwave,longwave radiation
! TKR1,TKS1=litter,soil temperature
! VOLWR2,VOLW12=litter,soil water volume
! VHCPR2,VHCP12=litter,soil heat capacity
!
! litter layer
! albedo
  ALBL=(0.20_r8*BKVL(0,NY,NX)+0.06_r8*VOLW1(0,NY,NX)+0.30_r8 &
    *VOLI1(0,NY,NX))/(BKVL(0,NY,NX)+VOLW1(0,NY,NX)+VOLI1(0,NY,NX))

  !radiation incident on litter layer  
  RFLX0=(1.0_r8-ALBL)*RADXR(NY,NX)+THRYR(NY,NX)  
  !kelvin, initial litter layer temperature
  TKR1=TK1(0,NY,NX)                              
  !volumetric water content
  VOLWR2=VOLW1(0,NY,NX)                          
  !heat capacity, initialized with residual layer
  VHCPR2=VHCP1(0,NY,NX)                          

! top soil layer
  TKS1=TK1(NUM(NY,NX),NY,NX)
  VOLW12=VOLW1(NUM(NY,NX),NY,NX)
  VHCP12=VHCP1(NUM(NY,NX),NY,NX)
  !
  ! THERMAL CONDUCTIVITY BETWEEN SURFACE RESIDUE AND SOIL SURFACE
  !
  ! CNVR,CNV1=litter,soil vapor conductivity
  ! THETPM=litter air concentration
  ! POROS,POROQ=litter porosity, tortuosity
  ! WGSGR,WGSGL=litter,soil vapor diffusivity
  ! CVRD=litter cover fraction
  ! ATCNVR=litter-soil vapor conductance
  ! DLYRR,DLYR=litter,soil depths
  ! THETRR=dry litter concentration
  ! DTH*,RYL*,DNU*,TRB*=turbulence effects on thermal conductivity
  ! WTHET0,WTHET1=multiplier for air concn in thermal conductivity
  ! TCNDW*,TCNDA*=thermal conductivity of water,air
  ! TCNDR,TCND1=litter,soil thermal conductivity
  ! ATCNDR=litter-soil thermal conductance
  !
  CNVR=WGSGR(NY,NX)*THETPM(M,0,NY,NX)*POROQ*THETPM(M,0,NY,NX)/POROS(0,NY,NX)
  CNV1=WGSGL(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ &
    *THETPM(M,NUM(NY,NX),NY,NX)/POROS(NUM(NY,NX),NY,NX)

  IF(CVRD(NY,NX).GT.ZERO)THEN
    !there is litter layer
    IF(CNVR.GT.ZERO.AND.CNV1.GT.ZERO)THEN
      ATCNVR=2.0_r8*CNVR*CNV1/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV1*DLYRR(NY,NX))
    ELSE
      !below is a numerical hack
      ATCNVR=2.0_r8*CNVR/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))*CVRD(NY,NX)
    ENDIF
  ELSE
    ATCNVR=0.0_r8
  ENDIF
  !THETRR=litter in relative volume
  !THETPY=relative volume as air  
  THETRR=AZMAX1(1.0_r8-THETPX(0,NY,NX)-THETWX(0,NY,NX)-THETIX(0,NY,NX))

  DTKX=ABS(TK1(0,NY,NX)-TK1(NUM(NY,NX),NY,NX))*ppmc

  DTHW0=AZMAX1(THETWX(0,NY,NX)-TRBW)**3._r8
  DTHA0=AZMAX1(THETPX(0,NY,NX)-TRBA)**3._r8
  RYLXW0=DTKX*DTHW0
  RYLXA0=DTKX*DTHA0
  RYLNW0=AMIN1(1.0E+04_r8,RYLXW*RYLXW0)
  RYLNA0=AMIN1(1.0E+04_r8,RYLXA*RYLXA0)
  XNUSW0=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW0**0.25_r8/DNUSW)
  XNUSA0=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA0**0.25_r8/DNUSA)
  TCNDW0=2.067E-03_r8*XNUSW0
  TCNDA0=9.050E-05_r8*XNUSA0  
  WTHET0=1.467_r8-0.467_r8*THETPY(0,NY,NX)
  TCNDR=(0.779_r8*THETRR*9.050E-04_r8+0.622_r8*THETWX(0,NY,NX)*TCNDW0 &
    +0.380_r8*THETIX(0,NY,NX)*7.844E-03_r8+WTHET0*THETPX(0,NY,NX)*TCNDA0) &
    /(0.779_r8*THETRR+0.622_r8*THETWX(0,NY,NX)+0.380_r8*THETIX(0,NY,NX)+WTHET0*THETPX(0,NY,NX))

  call CalcThermConduc(NX,NY,NUM(NY,NX),DTKX,TCND1)

  ATCNDR=2.0_r8*TCNDR*TCND1/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX)+TCND1*DLYRR(NY,NX))
!
! SMALLER TIME STEP FOR SOLVING SURFACE RESIDUE ENERGY EXCHANGE
!
  call SurfLitterIterate(M,NY,NX,VOLWR2,ATCNDR,ATCNVR,PSISV1,RFLX0)

  end subroutine SRFLitterEnergyBalance

  !------------------------------------------------------------------------------------------

  subroutine SoilSRFEnerbyBalance(M,NY,NX,PSISV1,THRMA,RAR1,TopLayerWaterVolume)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(out):: PSISV1,THRMA
  real(r8),dimension(:,:), intent(inout) :: RAR1
  real(r8), dimension(:,:),intent(inout) :: TopLayerWaterVolume
  real(r8) :: RAa,VP1,PARE,PARS,RAGX,RFLX0,RI,WPLX,WPX,THETPX0
  real(r8) :: VOLWXG,VOLIXG,DFVR
  real(r8) :: TKX1
  real(r8) :: HFLX0,ALBG
! begin_execution
!
! PHYSICAL AND HYDRAULIC PROPERTIES OF SOIL SURFACE INCLUDING
! AIR AND WATER-FILLED POROSITY, AND WATER POTENTIAL USED IN
! FLUX CALCULATIONS
!
! THETY=current,hygroscopic water concentration
! POROS=porosity
! VOLW1,VOLXI=volume of micropore water
! VOLXI=soil volume less macropore,rock
! FC,WP=water contents at field capacity,wilting point
! FCL,WPL=log FC,WP
! FCD,PSD=FCL-WPL,log(POROS)-FCL
! PSISM1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
! PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
! PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
! SRP=parameter for deviation from linear log-log water retention
! function from hour1.f
! PSISO=osmotic potential
! BKVL=bulk soil mass for a given layer 

  call CalcSoilWaterPotential(NY,NX,NX,NY,NUM(NY,NX),PSISM1(NUM(NY,NX),NY,NX))

  PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
!
! IF(NX.EQ.4.AND.NY.EQ.5)THEN
!     WRITE(*,3232)'PSISV1',I,J,M,NX,NY,NUM(NY,NX),PSISV1
!    2,PSISM1(NUM(NY,NX),NY,NX),PSISO(NUM(NY,NX),NY,NX)
!    3,THETWX(NUM(NY,NX),NY,NX),THETW1,POROS(NUM(NY,NX),NY,NX)
!    4,PSL(NUM(NY,NX),NY,NX),LOG(THETW1),PSD(NUM(NY,NX),NY,NX)
!    5,VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX)
!    5,VOLX(NUM(NY,NX),NY,NX)
!    5,SRP(NUM(NY,NX),NY,NX)
!3232  FORMAT(A8,6I4,20E14.6)
! ENDIF
!
! SOIL SURFACE ALBEDO, NET RADIATION
!
! VOLW1,VOLI1=water,ice volume in micopores
! VOLWH1,VOLIH1=water,ice volume in macopores
! ALBG,ALBS=albedo of ground surface,soil
! BKVL=soil mass
! RADXG,THRYG,RFLXG=incoming shortwave,longwave,net radiation
! THRMA,THRMS=emitted longwave radiation, emissivity
! TK1=soil temperature
! albedo of water and ice are set to 0.06, and 0.30 respectively
  VOLWXG=VOLW1(NUM(NY,NX),NY,NX)+VOLWH1(NUM(NY,NX),NY,NX)
  VOLIXG=VOLI1(NUM(NY,NX),NY,NX)+VOLIH1(NUM(NY,NX),NY,NX)

  IF(VOLWXG+VOLIXG.GT.ZEROS2(NY,NX))THEN
    !top soil has water or ice
    !confirm if VOLI is defined as water equivalent, ice albedo seems too low.
    ALBG=(ALBS(NY,NX)*BKVL(NUM(NY,NX),NY,NX)+0.06_r8*VOLWXG &
      +0.30_r8*VOLIXG)/(BKVL(NUM(NY,NX),NY,NX)+VOLWXG+VOLIXG)
  ELSE
    ALBG=ALBS(NY,NX)
  ENDIF
  !absorbed radiation
  !RFLXG=net radiation, after taking out outgoing surface layer radiation  
  !THRMA=emitted longwave radiation  
  RFLX0=(1.0_r8-ALBG)*RADXG(NY,NX)+THRYG(NY,NX)
  THRMA=THRMS(NY,NX)*TK1(NUM(NY,NX),NY,NX)**4._r8
  RFLXG=RFLX0-THRMA
!
! AERODYNAMIC RESISTANCE ABOVE SOIL SURFACE INCLUDING
! RESISTANCE IMPOSED BY PLANT CANOPY
!
! THETPX*=air-filled porosity of soil
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
  THETPX0=AMAX1(ZERO,THETPX(0,NY,NX))
  DFVR=THETPX0*POROQ*THETPX0/POROS(0,NY,NX)
  RAR1(NY,NX)=RAG(NY,NX)+RAR(NY,NX)/DFVR
  RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB(NY,NX)*(TKQ(NY,NX)-TK1(NUM(NY,NX),NY,NX))))
  RAGX=AMAX1(RAM,0.8_r8*RAGS(NY,NX),AMIN1(1.2_r8*RAGS(NY,NX),RAR1(NY,NX)/(1.0_r8-10.0_r8*RI)))
  RAGS(NY,NX)=RAGX
  RAa=RAGR(NY,NX)+RAGS(NY,NX)
! IF(I.EQ.63.AND.NX.EQ.1)THEN
!     WRITE(*,7776)'RAGX',I,J,M,NX,NY,RAGZ,BARE(NY,NX),RAG(NY,NX)
!    2,CVRDW(NY,NX),RAR1,RI,RIB(NY,NX),TKQ(NY,NX),TK1(NUM(NY,NX),NY,NX)
!    3,TK1(0,NY,NX),RAGX,RAM,RAGS(NY,NX),RA
!    4,RAR(NY,NX),DFVR,THETPX0,POROQ,THETPX(0,NY,NX)
!    5,DLYRR(NY,NX),WGSGR(NY,NX)
!7776  FORMAT(A8,5I4,30E12.4)
! ENDIF
!
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!
! PARE,PARS=blcs for latent,sensible heat fluxes over soil
! PAREG,PARSG=conductances for latent,sensible heat fluxes
! RZ=minimum surface resistance
! VP1,VPQ=vapor pressure at soil surface, canopy air
! EVAPG=evaporation
! EFLXG=latent heat flux
! XH=rate constant
! VOLW2=soil water volume
! VAP=latent heat of evaporation
! VFLXG=convective heat of evaporation flux
!
  PARE=PAREG(NY,NX)/(RAa+RZ)
  PARS=PARSG(NY,NX)/RAa

  TKX1=TK1(NUM(NY,NX),NY,NX)
  VP1=vapsat(TKX1)*EXP(18.0*PSISV1/(RGAS*TKX1))

  !evaporation, no more than what is available, g H2O
  EVAPG(NY,NX)=AMAX1(PARE*(VPQ(NY,NX)-VP1),-AZMAX1(TopLayerWaterVolume(NY,NX)*XNPX))
  !latent heat
  EFLXG=EVAPG(NY,NX)*VAP
  IF(EVAPG(NY,NX).LT.0.0_r8)THEN
    !evaporation
    VFLXG=EVAPG(NY,NX)*cpw*TK1(NUM(NY,NX),NY,NX)
  ELSE
    !condensation
    VFLXG=EVAPG(NY,NX)*cpw*TKQ(NY,NX)
  ENDIF
  !take away water from evaporation
  TopLayerWaterVolume(NY,NX)=TopLayerWaterVolume(NY,NX)+EVAPG(NY,NX)
!
! SOLVE FOR SOIL SURFACE TEMPERATURE AT WHICH ENERGY
! BALANCE OCCURS, SOLVE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
! SFLXG,EFLXG,RFLXG=sensible,latent heat fluxes, net radiation
! VFLXG=convective heat flux from EFLXG
! HeatFlux2Ground=storage heat flux
!
  SFLXG=PARS*(TKQ(NY,NX)-TK1(NUM(NY,NX),NY,NX))
  !net energy into soil, subtracting latent heat and sensible heat
  HFLX0=RFLXG+EFLXG+SFLXG
  !total heat plus convective heat 
  HeatFlux2Ground=HFLX0+VFLXG

  end subroutine SoilSRFEnerbyBalance

!------------------------------------------------------------------------------------------

  subroutine ExposedSoilFlux(M,NY,NX,RAR1,TopLayerWaterVolume)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),dimension(:,:),intent(inout) :: RAR1
  real(r8), dimension(:,:),intent(inout) :: TopLayerWaterVolume  
  real(r8) :: PSISV1
  real(r8) :: THRMA

! begin_execution
! Watch out for L, is its value defined?
  call SoilSRFEnerbyBalance(M,NY,NX,PSISV1,THRMA,RAR1,TopLayerWaterVolume)
!
  call SRFLitterEnergyBalance(M,NY,NX,PSISV1)
!
  call SumAftEnergyBalance(NY,NX,THRMA)

  end subroutine ExposedSoilFlux

!------------------------------------------------------------------------------------------

  subroutine AtmLandSurfExchange(M,NY,NX,RAR1,TopLayerWaterVolume)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), dimension(:,:), intent(inout) :: RAR1
  real(r8), dimension(:,:),intent(inout) :: TopLayerWaterVolume  
! begin_execution

  !solve if there is significant snow layer 
  
  IF(VHCPWM(M,1,NY,NX).GT.VHCPWX(NY,NX))THEN
!   VHCPW,VHCPWX=current, minimum snowpack heat capacities
    call SolveSnowpack(M,NY,NX,EFLXW,RFLXW,VFLXW,SFLXW,HFLXW,&
      FLWHLW,FLWLW,FLWLXW,FLWRLW,HFLWRLW,HFLWLW)
  ENDIF
!
! ENERGY EXCHANGE AT SOIL SURFACE IF EXPOSED UNDER SNOWPACK
! FSNW,FSNX=fractions of snow,snow-free cover
  IF(FSNX(NY,NX).GT.0.0_r8.AND.(BKDS(NUM(NY,NX),NY,NX).GT.ZERO.OR.VHCP1(NUM(NY,NX),NY,NX).GT.VHCPNX(NY,NX)))THEN
   !
    call ExposedSoilFlux(M,NY,NX,RAR1,TopLayerWaterVolume)
!
  ELSE
!   ground is fully snow covered, thus no flux from soil & litter
    call SnowCoveredTopSoilFlux(NY,NX)
  ENDIF
!
! AGGREGATE RESIDUE AND SOIL SURFACE FLUXES BENEATH SNOW
! AND ATMOSPHERE
!
! FLWL,FLWLX=total water flux into soil micropores
! FLWHL=total water flux into soil macropores
! HFLWL=total heat flux into soil
! FLWRL,FLWLX=total water flux into litter
! HFLWRL=total heat flux into litter
! FLWV*=total internal vapor flux in soil
!
  FLWL(3,NUM(NY,NX),NY,NX)=FLWLW+FLWLG
  FLWLX(3,NUM(NY,NX),NY,NX)=FLWLXW+FLWLXG
  FLWHL(3,NUM(NY,NX),NY,NX)=FLWHLW+FLWHLG
  HFLWL(3,NUM(NY,NX),NY,NX)=HFLWLW+HFLWLG
  FLWRL(NY,NX)=FLWRLW+FLWRLG
  HFLWRL(NY,NX)=HFLWRLW+HFLWRLG

  end subroutine AtmLandSurfExchange  
!------------------------------------------------------------------------------------------

  subroutine SnowCoveredTopSoilFlux(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

! begin_execution
  RFLXG=0.0_r8   !net radiation into soil
  EFLXG=0.0_r8   !latent heat flux from air and topsoil
  VFLXG=0.0_r8   !convective heat flux from air and topsoil
  SFLXG=0.0_r8   !sensible heat flux from air to topsoil
  HeatFlux2Ground=0.0_r8   !net heat flux into from air into topsoil
  RFLXR=0.0_r8   !net radiation flux into litter
  EFLXR=0.0_r8   !latent heat flux from air and litter
  VFLXR=0.0_r8   !convective heat flux from air to litter
  SFLXR=0.0_r8   !sensible heat flux from air to litter
  HFLXR=0.0_r8   !net heat flux from air to litter
  FLWLG=0.0_r8   !total water flux from air into soil
  FLWLXG=0.0_r8
  FLWHLG=0.0_r8        !total water flux from air into macropore
  HFLWLG=0.0_r8        !total water associated heat flux from air into soil
  FLWRLG=0.0_r8        !total water flux from air to litter
  HFLWRLG=0.0_r8       !total water associated heat flux from air to litter
  EVAPG(NY,NX)=0.0_r8  !evaporative flux from air into soil
  EVAPR(NY,NX)=0.0_r8  !evaporative flux from air into litter
  end subroutine SnowCoveredTopSoilFlux  
!------------------------------------------------------------------------------------------

  subroutine SurfSoilResidueWaterCapillExch(M,NY,NX,FKSAT)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FKSAT

  real(r8) :: THETW1,THETWR,PSIST1
  real(r8) :: PSIST0,HFLQR,FLQZ,FLQX
  real(r8) :: FLQR,FLQ2,CNDR,AVCNDR
  real(r8) :: CND1
  integer :: K0,K1

! begin_execution
! CNDR,HCNDR=current,saturated litter hydraulic conductivity
! PSISE,PSISM1(0,=air entry,current litter water potential
! VOLW1(0,VOLWRX=current,maximum litter water volume
! CND1,HCND=soil hydraulic conductivity
! FKSAT=reduction in soil surface Ksat from rainfall energy impact
! K1=soil relative water-filled porosity
! THETWX,POROS=soil water content,porosity
! AVCNDR=litter-soil hydraulic conductance
! DLYRR,DLYR=litter,soil thicknesses
! FLQR,FLQL=litter-soil water flux unltd,ltd by water
! XNPH=time step of flux calculations
! CVRD=fraction of litter cover
! PSISM1(NUM=soil matric water potential
! VOLW1(NUM=soil water volume
! HFLQL=convective heat from litter-soil water flux
! FLWL,HFLWL=micropore water,heat flux
! FLWRL,HFLWRL=total litter water,heat flux
! FLWRM=litter-soil water flux for solute transfer in trnsfr.f
! CND1,CNDL=hydraulic conductivity of source,destination layer
! HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
!
  IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      !surface litter holds water
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
    ELSE
      THETWR=POROS0(NY,NX)
    ENDIF
    THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX),AMIN1(POROS(NUM(NY,NX),NY,NX) &
      ,safe_adb(VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX))))
    !litter layer  
    K0=MAX(1,MIN(100,INT(100.0*(AZMAX1(POROS0(NY,NX)-THETWR))/POROS0(NY,NX))+1))
    !topsoil layer
    K1=MAX(1,MIN(100,INT(100.0*(AZMAX1(POROS(NUM(NY,NX),NY,NX)-THETW1))/POROS(NUM(NY,NX),NY,NX))+1))
    CNDR=HCND(3,K0,0,NY,NX)
    CND1=HCND(3,K1,NUM(NY,NX),NY,NX)*FKSAT
    AVCNDR=2.0_r8*CNDR*CND1/(CNDR*DLYR(3,NUM(NY,NX),NY,NX)+CND1*DLYRR(NY,NX))
    PSIST0=PSISM1(0,NY,NX)+PSISH(0,NY,NX)+PSISO(0,NY,NX)
    PSIST1=PSISM1(NUM(NY,NX),NY,NX)+PSISH(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
    !FLQX=water flux from litter layer into the topsoil
    FLQX=AVCNDR*(PSIST0-PSIST1)*AREA(3,NUM(NY,NX),NY,NX)*CVRDW(NY,NX)*XNPH

    IF(FLQX.GE.0.0_r8)THEN
      IF(THETWR.GT.THETS(0,NY,NX))THEN
        FLQZ=FLQX+AMIN1((THETWR-THETS(0,NY,NX))*VOLR(NY,NX),&
          AZMAX1((THETS(NUM(NY,NX),NY,NX)-THETW1)*VOLY(NUM(NY,NX),NY,NX)))*XNPX
      ELSE
        FLQZ=FLQX
      ENDIF
      FLQR=AZMAX1(AMIN1(FLQZ,VOLW1(0,NY,NX)*XNPX,VOLP1(NUM(NY,NX),NY,NX)))
      FLQ2=AZMAX1(AMIN1(FLQX,VOLW1(0,NY,NX)*XNPX,VOLP1(NUM(NY,NX),NY,NX)))
    ELSE
      IF(THETW1.GT.THETS(NUM(NY,NX),NY,NX))THEN
        FLQZ=FLQX+AMAX1((THETS(NUM(NY,NX),NY,NX)-THETW1)*VOLY(NUM(NY,NX),NY,NX),&
          AZMIN1((THETWR-THETS(0,NY,NX))*VOLR(NY,NX)))*XNPX
      ELSE
        FLQZ=FLQX
      ENDIF
      FLQR=AZMIN1(AMAX1(FLQZ,-VOLW1(NUM(NY,NX),NY,NX)*XNPX,-VOLP1(0,NY,NX)))
      FLQ2=AZMIN1(AMAX1(FLQX,-VOLW1(NUM(NY,NX),NY,NX)*XNPX,-VOLP1(0,NY,NX)))
    ENDIF

    IF(VOLP1Z(NUM(NY,NX),NY,NX).LT.0.0_r8)THEN
      FLQR=FLQR+AZMIN1(AMAX1(-VOLW1(NUM(NY,NX),NY,NX)*XNPX,VOLP1Z(NUM(NY,NX),NY,NX)))
      FLQ2=FLQ2+AZMIN1(AMAX1(-VOLW1(NUM(NY,NX),NY,NX)*XNPX,VOLP1Z(NUM(NY,NX),NY,NX)))
    ENDIF

    IF(FLQR.GT.0.0_r8)THEN
      HFLQR=cpw*TK1(0,NY,NX)*FLQR
    ELSE
      HFLQR=cpw*TK1(NUM(NY,NX),NY,NX)*FLQR
    ENDIF
    FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
    if(abs(FLWL(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qFLWL(3,NUM(NY,NX),NY,NX)=',FLQR
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
    FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
    HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR
    FLWRM(M,NY,NX)=FLQR

  ELSE
    FLQR=XVOLW(NY,NX)*XNPX
    HFLQR=cpw*TK1(0,NY,NX)*FLQR
    FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
    if(abs(FLWL(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qrFLWL(3,NUM(NY,NX),NY,NX)=',XVOLW(NY,NX),XNPX
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
    FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
    HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR

    FLWRM(M,NY,NX)=FLQR

  ENDIF

  end subroutine SurfSoilResidueWaterCapillExch
!------------------------------------------------------------------------------------------

  subroutine InfilSRFRoffPartition(M,NY,NX,N1,N2)
  implicit none
  integer, intent(in) :: M,NY,NX
  integer, intent(out):: N1,N2
  real(r8) :: TFLX1,TK1X,ENGYR,VOLW1X,VHCP1X
  real(r8) :: TFREEZ,TFLX,VX
  real(r8) :: HFLQHR,FLQHR  
  real(r8) :: D,Q,R,V

!     begin_execution
!     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
!     OF THE LITTER IS EXCEEDED
!
!     VOLPH1=air-filled macroporosity
!     FINHR,HFINHR=water,convective heat from litter to macropores
!     VOLW1(0,VOLWRX=current,maximum litter water volume
!     FLWL,HFLWL=micropore water,heat flux
!     FLWRL,HFLWRL=total litter water,heat flux
!
  IF(VOLPH1(NUM(NY,NX),NY,NX).GT.0.0_r8 .AND.XVOLW(NY,NX).GT.0.0_r8)THEN
    FLQHR=AMIN1(XVOLW(NY,NX)*XNPX,VOLPH1(NUM(NY,NX),NY,NX))
    HFLQHR=FLQHR*cpw*TK1(0,NY,NX)
    FLWHL(3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)+FLQHR
    HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQHR
    FLWRL(NY,NX)=FLWRL(NY,NX)-FLQHR
    HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQHR
  ENDIF
!
!     FREEZE-THAW IN RESIDUE SURFACE FROM NET CHANGE IN RESIDUE
!     SURFACE HEAT STORAGE
!
!     TFREEZ=litter freezing temperature
!     PSISM1=litter water potential
!     VOLW1*,VOLI1=litter water,ice volume
!     VHCP1*=litter volumetric heat capacity
!     TK1*=litter temperature
!     ORGC=litter organic C
!     HFLWRL=total litter conductive, convective heat flux
!     TFLX1,TFLX=unltd,ltd latent heat from freeze-thaw
!     TFLXR,WFLXR=litter water,latent heat flux from freeze-thaw
!
  TFREEZ=-9.0959E+04_r8/(PSISM1(0,NY,NX)-333.0_r8)
  VOLW1X=AZMAX1(VOLW1(0,NY,NX)+FLWRL(NY,NX))
  ENGYR=VHCP1(0,NY,NX)*TK1(0,NY,NX)
  VHCP1X=cpo*ORGC(0,NY,NX)+cpw*VOLW1X+cpi*VOLI1(0,NY,NX)
  IF(VHCP1X.GT.ZEROS(NY,NX))THEN
    TK1X=(ENGYR+HFLWRL(NY,NX))/VHCP1X
  ELSE
    TK1X=TK1(0,NY,NX)
  ENDIF
  IF((TK1X.LT.TFREEZ.AND.VOLW1(0,NY,NX).GT.ZERO*VOLT(0,NY,NX)) &
    .OR.(TK1X.GT.TFREEZ.AND.VOLI1(0,NY,NX).GT.ZERO*VOLT(0,NY,NX)))THEN
    TFLX1=VHCP1(0,NY,NX)*(TFREEZ-TK1X) &
      /((1.0_r8+TFREEZ*6.2913E-03_r8)*(1.0_r8-0.10_r8*PSISM1(0,NY,NX)))*XNPX
    IF(TFLX1.LT.0.0_r8)THEN
      TFLX=AMAX1(-333.0_r8*DENSI*VOLI1(0,NY,NX)*XNPX,TFLX1)
    ELSE
      TFLX=AMIN1(333.0_r8*VOLW1X*XNPX,TFLX1)
    ENDIF
    TFLXR(NY,NX)=TFLX
    WFLXR(NY,NX)=-TFLX/333.0_r8
  ELSE
    WFLXR(NY,NX)=0.0_r8
    TFLXR(NY,NX)=0.0_r8
  ENDIF
!
!     THICKNESS OF WATER FILMS IN LITTER AND SOIL SURFACE
!     FROM WATER POTENTIALS FOR GAS EXCHANGE IN TRNSFR.F
!
  IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
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
  IF(XVOLT(N2,N1).GT.VOLWG(N2,N1))THEN
    VX=XVOLT(N2,N1)-VOLWG(N2,N1)
    D=VX/AREA(3,0,N2,N1)
    R=D/2.828_r8
    V=R**0.67_r8*safe_adb(SQRT(SLOPE(0,N2,N1)),ZM(N2,N1))
    Q=V*D*AREA(3,NUM(N2,N1),N2,N1)*3.6E+03_r8*XNPH
    VOLW1X=AZMAX1(VOLW1(0,N2,N1)+WFLXR(N2,N1))
    QRM(M,N2,N1)=AMIN1(Q,VX*XNPX,VOLW1X*XNPX)*XVOLW(N2,N1)/XVOLT(N2,N1)
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
        ALT1=ALTG(N2,N1)+XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
!
!     EAST OR SOUTH RUNOFF
!
        IF(NN.EQ.1)THEN
          !destination grid elevation
          ALT2=ALTG(N5,N4)+XVOLT(N5,N4)/AREA(3,NU(N5,N4),N5,N4)
          IF(ALT1.GT.ALT2)THEN
            QRQ1=AZMAX1(((ALT1-ALT2)*AREA(3,NUM(N2,N1),N2,N1) &
              *AREA(3,NU(N5,N4),N5,N4)-XVOLT(N5,N4)*AREA(3,NUM(N2,N1),N2,N1) &
              +XVOLT(N2,N1)*AREA(3,NU(N5,N4),N5,N4)) &
              /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5,N4),N5,N4)))
            QR1(N,2,N5,N4)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
            HQR1(N,2,N5,N4)=cpw*TK1(0,N2,N1)*QR1(N,2,N5,N4)
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
            ALTB=ALTG(N5B,N4B)+XVOLT(N5B,N4B)/AREA(3,NU(N5,N4B),N5B,N4B)
            IF(ALT1.GT.ALTB)THEN
              QRQ1=AZMAX1(((ALT1-ALTB)*AREA(3,NUM(N2,N1),N2,N1) &
                *AREA(3,NU(N5B,N4B),N5B,N4B)-XVOLT(N5B,N4B) &
                *AREA(3,NUM(N2,N1),N2,N1) &
                +XVOLT(N2,N1)*AREA(3,NU(N5B,N4B),N5B,N4B)) &
                /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5B,N4B),N5B,N4B)))
              QR1(N,1,N5B,N4B)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
              HQR1(N,1,N5B,N4B)=cpw*TK1(0,N2,N1)*QR1(N,1,N5B,N4B)
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

  subroutine AccumWaterVaporHeatFluxes(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
! begin_execution
! HOURLY-ACCUMULATED WATER, VAPOR AND HEAT FLUXES THROUGH
! SURFACE RESIDUE AND SOIL SURFACE
!
! THAWR,HTHAWR=litter water,heat fluxes from freeze-thaw
! FLW,FLWH,HFLW=soil micropore,macropore,heat fluxes
! FLWR,HFLWR=litter water,heat fluxes
! FLSW,FLSWH=water from snowpack to soil micropores,macropores
! HFLSW=convective heat from snowpack to soil
! FLSWR=water flux from snowpack to litter
! HFLSWR=convective heat flux from snowpack to litter
! HEATI,HEATE,HEATS,HEATG=net radiation,latent,sensible,storage heat
! TEVAPG=total evaporation
! FLWM,FLWHM=water flux into soil micropore,macropore for use in trnsfr.f
! VOLWX1=VOLW1 accounting for wetting front
!
  THAWR(NY,NX)=THAWR(NY,NX)+WFLXR(NY,NX)
  HTHAWR(NY,NX)=HTHAWR(NY,NX)+TFLXR(NY,NX)
  FLW(3,NUM(NY,NX),NY,NX)=FLW(3,NUM(NY,NX),NY,NX)+FLWL(3,NUM(NY,NX),NY,NX)
  FLWX(3,NUM(NY,NX),NY,NX)=FLWX(3,NUM(NY,NX),NY,NX)+FLWLX(3,NUM(NY,NX),NY,NX)
  FLWH(3,NUM(NY,NX),NY,NX)=FLWH(3,NUM(NY,NX),NY,NX)+FLWHL(3,NUM(NY,NX),NY,NX)
  HFLW(3,NUM(NY,NX),NY,NX)=HFLW(3,NUM(NY,NX),NY,NX)+HFLWL(3,NUM(NY,NX),NY,NX)
  FLWR(NY,NX)=FLWR(NY,NX)+FLWRL(NY,NX)
  HFLWR(NY,NX)=HFLWR(NY,NX)+HFLWRL(NY,NX)
  HEATI(NY,NX)=HEATI(NY,NX)+RFLXG+RFLXR+RFLXW
  HEATS(NY,NX)=HEATS(NY,NX)+SFLXG+SFLXR+SFLXW
  HEATE(NY,NX)=HEATE(NY,NX)+EFLXG+EFLXR+EFLXW
  HEATV(NY,NX)=HEATV(NY,NX)+VFLXG+VFLXR+VFLXW
  HEATH(NY,NX)=HEATH(NY,NX)+RFLXG+RFLXR+RFLXW &
    +SFLXG+SFLXR+SFLXW+EFLXG+EFLXR+EFLXW+VFLXG+VFLXR+VFLXW

  !EVAPG=evaporation from ground/top soil layer
  !EVAPR=evaporation from litter layer   
  !EVAPSN=evaporation from snow, sublimation+evaporation
  TEVAPG(NY,NX)=TEVAPG(NY,NX)+EVAPG(NY,NX)+EVAPR(NY,NX)+EVAPSN(NY,NX)
  FLWM(M,3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)
  FLWHM(M,3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)
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

  IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
    BAREW(NY,NX)=AZMAX1(BARE(NY,NX)-AMIN1(1.0_r8,AZMAX1(XVOLT(NY,NX)/VOLWD(NY,NX))))
  ELSE
    BAREW(NY,NX)=1.0_r8
  ENDIF
  CVRDW(NY,NX)=1.0_r8-BAREW(NY,NX)
  RAGS(NY,NX)=1.0_r8/(BAREW(NY,NX)/RAGR(NY,NX)+CVRDW(NY,NX)/RAR1(NY,NX))
  RAS=SnowBNDResistance(NY,NX)
  PARG(M,NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*XNPH/(RAGS(NY,NX)+RAS)
!
! REDISTRIBUTE INCOMING PRECIPITATION
! BETWEEN RESIDUE AND SOIL SURFACE
!
! BKDS=bulk density
! FLQRS,FLQRH=water flux from soil micropores,macropores to litter
! FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
! VOLP1,VOLPH1=air-filled microporosity,macroporosity
! HFLQR1=convective heat flux from soil to litter
! FLYM,HWFLYM=total water flux, convective heat flux to litter
! FLQM,FLHM=total water flux to soil micropores, macropores
! HWFLQM=total convective heat flux to soil micropores, macropores
! XNPR=time step for litter water,heat flux calculations
!
  IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    FLQRS=AZMAX1(FLQ1(NY,NX)-VOLP1(NUM(NY,NX),NY,NX))
    FLQRH=AZMAX1(FLH1(NY,NX)-VOLPH1(NUM(NY,NX),NY,NX))
    HFLQR1=cpw*TairK(NY,NX)*(FLQRS+FLQRH)
    FLYM=FLY1(NY,NX)+FLQRS+FLQRH
    HWFLYM=HWFLY1(NY,NX)+HFLQR1
    FLQM=FLQ1(NY,NX)-FLQRS
    FLHM=FLH1(NY,NX)-FLQRH
    HWFLQM=HWFLQ1(NY,NX)-HFLQR1
  ELSE
    FLYM=FLY1(NY,NX)
    HWFLYM=HWFLY1(NY,NX)
    FLQM=FLQ1(NY,NX)
    FLHM=FLH1(NY,NX)
    HWFLQM=HWFLQ1(NY,NX)
  ENDIF
  FLYM2=FLYM*XNPR
  HWFLM2=HWFLYM*XNPR
!
! WATER GAS EXCHANGE COEFFICIENTS IN SURFACE LITTER
!
! VOLA1,VOLI1,VOLW1,VOLPM=total,ice-,water-,air-filled porosity
! TFND1=temperature effect on gas diffusivity
! DFGS=rate constant for air-water gas exchange
! Z1R,Z2RW,Z2RD,Z3RX=parameters for litter air-water gas transfers
! XNPD=time step for gas transfer calculations, it is tunable parameter
! TORT=tortuosity for aqueous diffusivity
! VOLAT0=ice-excluded porosity in litter

  VOLAT0=VOLA10(NY,NX)-VOLI1(0,NY,NX)
  IF(VOLAT0.GT.ZEROS2(NY,NX).AND.VOLPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    !litter layer is not saturated
    THETWA=AZMAX1(AMIN1(1.0_r8,VOLW1(0,NY,NX)/VOLAT0))
    TFND1=TEFAQUDIF(TK1(0,NY,NX))
    scalar=TFND1*XNPD
    DFGS(M,0,NY,NX)=fDFGS(scalar,THETWA,0.0_r8,is_litter=.true.)
  ELSE
    !litter layer saturated
    DFGS(M,0,NY,NX)=0.0_r8
  ENDIF
! VOLWRX=surface litter water holding capacity, [m3 d-2]
  IF(VOLWRX(NY,NX).GT.ZEROS(NY,NX))THEN
    !litter is able to hold water
    THETWT=AMIN1(1.0_r8,VOLW(0,NY,NX)/VOLWRX(NY,NX))
  ELSE
    THETWT=1.0
  ENDIF
  !TORT=tortuosity in litter (treated as micropore)
  TORT(M,0,NY,NX)=TortMicporeW(THETWT)
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
    HV=1.0E+03_r8*AZMAX1(XVOLT(NY,NX)-VOLWG(NY,NX))/AREA(3,NUM(NY,NX),NY,NX)
    ENGYPM(M,NY,NX)=(ENGYD*PRECD(NY,NX)+ENGYB*PRECB(NY,NX))*EXP(-2.0_r8*HV)*BARE(NY,NX)*XNPH
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
! TWFLXL,TWFLXH=total freeze-thaw in micropores,macropores
! TTFLXL=total latent heat from freeze-thaw
! TFLWL,TFLWHL=net water flux in micropores,macropores
! THFLWL=net heat flux
!
!
! ENERGY EXCHANGE VARIABLES AT SNOW SURFACE IF PRESENT
!
! RFLXW,EFLXW,VFLXW,SFLXW,HFLXW=netradn,latent,convective,sensible
! and storage heat fluxes
! FLWLW,FLWHLW=water from snowpack to soil micropores,macropores
! HFLWLW=conv heat from snowpack to soil micropores,macropores
! FLWRLW=water flux from snowpack to litter
! HFLWRLW=convective heat flux from snowpack to litter
!
  RFLXW=0.0_r8
  EFLXW=0.0_r8
  VFLXW=0.0_r8
  SFLXW=0.0_r8
  HFLXW=0.0_r8
  FLWLW=0.0_r8
  FLWLXW=0.0_r8
  FLWHLW=0.0_r8
  HFLWLW=0.0_r8
  FLWRLW=0.0_r8
  HFLWRLW=0.0_r8
!
! EVAPS,EVAPW=evaporation from soil,snowpack surfaces
! FLQRM,FLQSM,FLQHM=water into litter,soil micropores,micropores for use in trnsfr.f
!
  EVAPSN(NY,NX)=0._r8;EVAPS(NY,NX)=0.0_r8;EVAPW(NY,NX)=0.0_r8

  FLQRM(M,NY,NX)=0.0_r8
  FLQSM(M,NY,NX)=0.0_r8
  FLQHM(M,NY,NX)=0.0_r8
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
  ! DTBLX=natural water table depth
  ! QR1,HQR1=runoff, convective heat from runoff
  ! QR,HQR=hourly-accumulated runoff, convective heat from runoff
  ! QRM,QRV=runoff,velocity for erosion, solute transfer
  ! XN=direction
  !
  ! RUNOFF
  !
  DPTHW1=XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
  DPTHW2=VOLWG(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
  ALT1=ALTG(N2,N1)+DPTHW1
  ALT2=ALTG(N2,N1)+DPTHW2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)
  IF(ALT1.GT.ALT2.AND.CumDepth2LayerBottom(NU(N2,N1)-1,N2,N1)-DPTHW1.LT.DTBLX(N2,N1))THEN
    QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HQR1(N,NN,M5,M4)=cpw*TK1(0,N2,N1)*QR1(N,NN,M5,M4)
    QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
    HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
!
! RUNON
!
  ELSEIF(CumDepth2LayerBottom(NU(N2,N1)-1,N2,N1)-DPTHW1.GT.DTBLX(N2,N1))THEN
    VX=AZMIN1((DTBLX(N2,N1)-CumDepth2LayerBottom(NU(N2,N1)-1,N2,N1)+DPTHW1)*AREA(3,NUM(N2,N1),N2,N1))
    QRM(M,N2,N1)=VX*XNPX
    QRV(M,N2,N1)=0.0_r8
    QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
    HQR1(N,NN,M5,M4)=cpw*TK1(0,N2,N1)*QR1(N,NN,M5,M4)
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

  real(r8) :: FLWSW,HFLWSW
  real(r8) :: FLWQAS,FLWQAH,FLWQW
  real(r8) :: FLWQBX,FLWQAX,HFLWQB,HFLWQA
  real(r8) :: FLWQB  
!     PRECA=precipitation+irrigation
!     PRECD,PRECB=direct,indirect precipn+irrign at soil surface
!     TFLWCI=net ice transfer to canopy, updated in hour1


  !convert water flux from m/hour to mm/hour
  PRECM(NY,NX)=1.0E+03_r8*PRECA(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  PRECD(NY,NX)=1.0E+03_r8*(PRECA(NY,NX)-TFLWCI(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  PRECB(NY,NX)=1.0E+03_r8*(TFLWCI(NY,NX)-TFLWC(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
!
!     RESIDUE WATER ABSORPTION CAPACITY
!
!     HCNDR=litter saturated hydraulic conductivity
!     DLYRR=litter depth
!
  HCNDR(NY,NX)=HCNDRR
  DLYRR(NY,NX)=AMAX1(2.5E-03_r8,DLYR(3,0,NY,NX))
!
!     DISTRIBUTION OF PRECIPITATION AND ITS HEAT AMONG SURFACE
!     RESIDUE, SOIL SURFACE, AND MACROPORES
!
!     PRECA,PRECW=rainfall+irrigation,snowfall (water equiv)
!     FLWQW=rainfall to snowpack
!     FLWSW=snowfall to snowpack
!     HFLWSW=convective heat flux to snowpack
!     FLWQB=precip to litter+soil surfaces
!     FLWQAX,FLWQBX=precip to soil,litter surfaces
!     HFLWQA,HFLWQB=convective heat flux to soil,litter surfaces
!     FLWQAS,FLWQAH=precip to soil micropores,macropores
!     TFLWC=canopy intercepted precipitation
!     FSNW=fraction of snow cover

  IF(PRECA(NY,NX).GT.0.0_r8.OR.PRECW(NY,NX).GT.0.0_r8)THEN
  ! there is precipitation
    FLWQW=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNW(NY,NX)
    FLWSW=PRECW(NY,NX)                                !snowfall
    HFLWSW=cps*TairK(NY,NX)*FLWSW+cpw*TairK(NY,NX)*FLWQW  !incoming heat flux from precipitations to snow-covered surface
    FLWQB=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNX(NY,NX)     !incoming precipitation to snow-free surface
    FLWQBX=FLWQB*CVRD(NY,NX)                          !water flux to snow-free coverd by litter
    HFLWQB=cpw*TairK(NY,NX)*FLWQBX                      !heat flux to snow-free surface covered by litter
    FLWQAX=FLWQB*BARE(NY,NX)                          !heat flux to snow-free surface not covered by litter
    HFLWQA=cpw*TairK(NY,NX)*FLWQAX
    FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)              !water flux to micropore
    FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)              !water flux to macropore
  ELSE
  ! no precipitation
    FLWQW=-TFLWC(NY,NX)*FSNW(NY,NX)                   !
    FLWSW=0.0_r8
    HFLWSW=cpw*TairK(NY,NX)*FLWQW
    FLWQB=-TFLWC(NY,NX)*FSNX(NY,NX)
    FLWQBX=FLWQB*CVRD(NY,NX)
    HFLWQB=cpw*TairK(NY,NX)*FLWQBX
    FLWQAX=FLWQB*BARE(NY,NX)
    HFLWQA=cpw*TairK(NY,NX)*FLWQAX
    FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)
    FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)
  ENDIF
!
!     PRECIP ON SNOW ARRAYS EXPORTED TO TRNSFR.F, TRNSFRS.F
!     FOR SOLUTE FLUX CALCULATIONS
!
!     PRECW,PRECR,PRECQ,PRECI=snow,rain,snow+rain,irrigation
!     VHCPW,VHCPWX=current, minimum snowpack heat capacities
!     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!
  IF(PRECW(NY,NX).GT.0.0_r8.OR.(PRECR(NY,NX).GT.0.0_r8 .AND.VHCPW(1,NY,NX).GT.VHCPWX(NY,NX)))THEN
    !there is precipitation, there is significant snow layer
    FLQRQ(NY,NX)=0.0_r8
    FLQRI(NY,NX)=0.0_r8
    FLQGQ(NY,NX)=PRECQ(NY,NX)
    FLQGI(NY,NX)=PRECI(NY,NX)
  ELSEIF((PRECQ(NY,NX).GT.0.0.OR.PRECI(NY,NX).GT.0.0_r8).AND.VHCPW(1,NY,NX).LE.VHCPWX(NY,NX))THEN
    !there is insignificant snow layer
    FLQRQ(NY,NX)=FLWQBX*PRECQ(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
    FLQRI(NY,NX)=FLWQBX*PRECI(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
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
!     XNPH=internal time step for fluxes through soil profile
!
!     FLW0S,FLQ0I,FLQ0W=snow,ice,water input to snowpack
!     HWFLQ0=convective heat flux to snowpack
!     FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
!     HWFLQ1,HWFLY1=convective heat flux to soil,litter surfaces
!
  FLQ0S(NY,NX)=FLWSW*XNPH
  FLQ0I(NY,NX)=0.0_r8
  FLQ0W(NY,NX)=FLWQW*XNPH
  HWFLQ0(NY,NX)=HFLWSW*XNPH

  FLQ1(NY,NX)=FLWQAS*XNPH
  FLH1(NY,NX)=FLWQAH*XNPH

  FLY1(NY,NX)=FLWQBX*XNPH
  HWFLQ1(NY,NX)=HFLWQA*XNPH
  HWFLY1(NY,NX)=HFLWQB*XNPH

  end subroutine PartitionPrecip
!------------------------------------------------------------------------------------------  

  subroutine UpdateSurfaceAtM(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VOLIRZ,ENGYR,VHCPXX
  real(r8) :: TK0XX,tk1pres,TVOLWI,VOLWRZ

  !update snow
  call UpdateSnowAtM(M,NY,NX)

  ! SURFACE RESIDUE WATER AND TEMPERATURE
  !
  ! XVOLT,XVOLW=free water+ice,water in litter layer
  ! VOLWM,VOLPM=surface water,air content for use in trnsfr.f
  ! VOLWRX=maximum water retention by litter
  ! VHCP1=volumetric heat capacity of litter
  ! VOLA1,VOLW1,VOLI1,VOLP1=pore,water,ice,air volumes of litter
  ! VOLWRX=maximum water retention by litter
  ! TFLXR,WFLXR=litter water,latent heat flux from freeze-thaw
  ! VOLR=dry litter volume
  ! THETWX,THETIX,THETPX=water,ice,air concentrations
  ! VHCP1=volumetric heat capacity of litter
  ! TK1=litter temperature
  ! HFLWRL,TFLXR,THQR1=litter total cond+conv,latent,runoff heat flux

  VOLW1(0,NY,NX)=AZMAX1(VOLW1(0,NY,NX)+FLWRL(NY,NX)+WFLXR(NY,NX)+TQR1(NY,NX))
  VOLI1(0,NY,NX)=AZMAX1(VOLI1(0,NY,NX)-WFLXR(NY,NX)/DENSI)
  VOLP1(0,NY,NX)=AZMAX1(VOLA10(NY,NX)-VOLW1(0,NY,NX)-VOLI1(0,NY,NX))
  VOLWM(M+1,0,NY,NX)=VOLW1(0,NY,NX)
  VOLPM(M+1,0,NY,NX)=VOLP1(0,NY,NX)
  TVOLWI=VOLW1(0,NY,NX)+VOLI1(0,NY,NX)
  XVOLT(NY,NX)=AZMAX1(TVOLWI-VOLWRX(NY,NX))
  IF(TVOLWI.GT.ZEROS(NY,NX))THEN
    VOLWRZ=VOLW1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
    VOLIRZ=VOLI1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
    XVOLW(NY,NX)=AZMAX1(VOLW1(0,NY,NX)-VOLWRZ)
    XVOLI(NY,NX)=AZMAX1(VOLI1(0,NY,NX)-VOLIRZ)
  ELSE
    XVOLW(NY,NX)=0.0_r8
    XVOLI(NY,NX)=0.0_r8
  ENDIF
  XVOLTM(M+1,NY,NX)=XVOLT(NY,NX)
  XVOLWM(M+1,NY,NX)=XVOLW(NY,NX)
  XVOLIM(M+1,NY,NX)=XVOLI(NY,NX)
  IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
    THETWX(0,NY,NX)=AZMAX1t(VOLW1(0,NY,NX)/VOLR(NY,NX))
    THETIX(0,NY,NX)=AZMAX1t(VOLI1(0,NY,NX)/VOLR(NY,NX))
    THETPX(0,NY,NX)=AZMAX1t(VOLP1(0,NY,NX)/VOLR(NY,NX))*AZMAX1t((1.0_r8-XVOLT(NY,NX)/VOLWD(NY,NX)))
  ELSE
    THETWX(0,NY,NX)=0.0_r8
    THETIX(0,NY,NX)=0.0_r8
    THETPX(0,NY,NX)=1.0_r8
  ENDIF
  THETPM(M+1,0,NY,NX)=THETPX(0,NY,NX)
  VHCPXX=VHCP1(0,NY,NX)              !heat capacity
  TK0XX=TK1(0,NY,NX)                 !residual temperature
  ENGYR=VHCP1(0,NY,NX)*TK1(0,NY,NX)  !initial energy content
  VHCP1(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW1(0,NY,NX)+cpi*VOLI1(0,NY,NX)  !update heat capacity
  IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
    tk1pres=TK1(0,NY,NX)
    TK1(0,NY,NX)=(ENGYR+HFLWRL(NY,NX)+TFLXR(NY,NX)+THQR1(NY,NX))/VHCP1(0,NY,NX)
    IF(ABS(VHCP1(0,NY,NX)/VHCPXX-1._r8)>0.025_r8.or. &
      abs(TK1(0,NY,NX)/tk1pres-1._r8)>0.025_r8)THEN
      TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
    ENDIF

  ELSE
    TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
  ENDIF
  end subroutine UpdateSurfaceAtM
!------------------------------------------------------------------------------------------

  subroutine SumAftEnergyBalance(NY,NX,THRMA)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in):: THRMA
  real(r8) :: FLWVLS  
! begin_execution
!
! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
! FOR LATER UPDATES TO STATE VARIABLES
!
! FLWLG,FLWHLG=water flux from atm to soil micropores,macropores
! HFLWLG=convective heat flux from atm to soil
! FLWRLG=water flux from atm to litter
! HFLWRLG=convective heat flux from atm to litter
! FLWVLS=water flux within soil accounting for wetting front
!
  FLWLG=FLQM+EVAPG(NY,NX)+FLV1
  if(abs(FLWLG)>1.e20_r8)then
    write(*,*)'FLQM+EVAPG(NY,NX)+FLV1',FLQM,EVAPG(NY,NX),FLV1
    call endrun(trim(mod_filename)//' at line',__LINE__)
  endif
  FLWLXG=FLQM+EVAPG(NY,NX)+FLV1
  FLWHLG=FLHM
  HFLWLG=HWFLQM+HeatFlux2Ground+HWFLV1+HFLCR1
  FLWRLG=FLYM+EVAPR(NY,NX)-FLV1
  HFLWRLG=HWFLYM+HFLXR-HWFLV1-HFLCR1

  FLWVLS=(VOLW1(NUM(NY,NX),NY,NX)-VOLWX1(NUM(NY,NX),NY,NX))*XNPH
!
! GENERATE NEW SNOWPACK
!
! XFLWS,XFLWW,XFLWI=hourly snow,water,ice transfer
! FLQ0S,FLQ0W,FLQ0I=snow,water,ice input to snowpack
! XHFLWW=hourly convective heat flux from snow,water,ice transfer
! HWFLQ0=convective heat flux from snow,water,ice to snowpack
!
  IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX).AND.FLQ0S(NY,NX).GT.ZEROS(NY,NX))THEN
    XFLWS(1,NY,NX)=XFLWS(1,NY,NX)+FLQ0S(NY,NX)
    XFLWW(1,NY,NX)=XFLWW(1,NY,NX)+FLQ0W(NY,NX)
    XFLWI(1,NY,NX)=XFLWI(1,NY,NX)+FLQ0I(NY,NX)
    XHFLWW(1,NY,NX)=XHFLWW(1,NY,NX)+HWFLQ0(NY,NX)
!     WRITE(*,4422)'INIT',I,J,FLQ0S(NY,NX),FLQ0W(NY,NX)
!    3,FLQ0I(NY,NX),HWFLQ0(NY,NX),XFLWS(1,NY,NX),XFLWW(1,NY,NX)
!    2,XFLWI(1,NY,NX),XHFLWW(1,NY,NX),HFLWL(3,NUM(NY,NX),NY,NX)
!    3,HFLW(3,NUM(NY,NX),NY,NX),FSNX(NY,NX),VHCP1(NUM(NY,NX),NY,NX)
!    4*TK1(NUM(NY,NX),NY,NX),HFLWRL(NY,NX),HFLWR(NY,NX)
!    5,VHCP1(0,NY,NX)*TK1(0,NY,NX),HEATH(NY,NX),RFLXG,RFLXR,RFLXW
!    2,SFLXG,SFLXR,SFLXW,EFLXG,EFLXR,EFLXW,VFLXG,VFLXR,VFLXW
  ENDIF
  !THRMG=longwave emission from litter and surface soil
  THRMG(NY,NX)=THRMG(NY,NX)+THRMA+THRMZ
  end subroutine SumAftEnergyBalance
!------------------------------------------------------------------------------------------
  subroutine SurfacePhysModel(M,NHE,NHW,NVS,NVN,RAR1,FKSAT,HeatFlux2Ground,TopLayerWaterVolume)
  implicit none
  integer, intent(in) :: M,NHE,NHW,NVS,NVN
  real(r8), dimension(:,:),intent(inout) :: RAR1
  REAL(R8), dimension(:,:),INTENT(OUT) :: FKSAT
  real(r8), dimension(:,:),intent(out) :: HeatFlux2Ground
  real(r8),dimension(:,:),intent(inout) :: TopLayerWaterVolume
  integer :: N1,N2,NX,NY

  D9895: DO  NX=NHW,NHE
    D9890: DO  NY=NVN,NVS
      call SurfaceEnergyModel(M,NX,NY,RAR1,FKSAT(NY,NX),HeatFlux2Ground(NY,NX),TopLayerWaterVolume)
      
    ! CAPILLARY EXCHANGE OF WATER BETWEEN SOIL SURFACE AND RESIDUE
      call SurfSoilResidueWaterCapillExch(M,NY,NX,FKSAT(NY,NX))

      call InfilSRFRoffPartition(M,NY,NX,N1,N2)
    !
      call LateralHydroExchange(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
    !
      call AccumWaterVaporHeatFluxes(M,NY,NX)
    ENDDO D9890
  ENDDO D9895

  end subroutine SurfacePhysModel

!------------------------------------------------------------------------------------------
  subroutine SurfaceEnergyModel(M,NX,NY,RAR1,FKSAT,HeatFlux2Ground1,TopLayerWaterVolume)
  implicit none
  integer, intent(in) :: M,NX,NY
  real(r8), dimension(:,:),intent(inout) :: RAR1
  REAL(R8),INTENT(OUT) :: FKSAT,HeatFlux2Ground1
  real(r8),dimension(:,:),intent(inout) :: TopLayerWaterVolume
  integer :: N1,N2

  !RAR1 is input
  call InitSurfModel(M,NY,NX,RAR1,FKSAT)

! updates RAR1
  call AtmLandSurfExchange(M,NY,NX,RAR1,TopLayerWaterVolume)
  HeatFlux2Ground1=HeatFlux2Ground

  end subroutine SurfaceEnergyModel

!------------------------------------------------------------------------------------------
  subroutine CalcThermConduc(N1,N2,N3,DTKX,TCND1)
  implicit none
  integer , intent(in) :: N1,N2,N3
  real(r8), intent(in) :: DTKX
  real(r8), intent(out):: TCND1

  real(r8) :: DTHW1,DTHA1,RYLXW1,RYLXA1,RYLNW1,RYLNA1
  REAL(R8) :: XNUSW1,XNUSA1,TCNDW1,TCNDA1,WTHET1

  IF(BKDS(N3,N2,N1).GT.ZERO.OR.THETWX(N3,N2,N1)+THETIX(N3,N2,N1).GT.ZERO)THEN
    !it is a soil layer or pure water layer
    DTHW1=AZMAX1(THETWX(N3,N2,N1)-TRBW)**3._r8
    DTHA1=AZMAX1(THETPX(N3,N2,N1)-TRBA)**3._r8
    RYLXW1=DTKX*DTHW1
    RYLXA1=DTKX*DTHA1
    RYLNW1=AMIN1(1.0E+04_r8,RYLXW*RYLXW1)
    RYLNA1=AMIN1(1.0E+04_r8,RYLXA*RYLXA1)
    XNUSW1=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW1**0.25_r8/DNUSW)
    XNUSA1=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA1**0.25_r8/DNUSA)
    TCNDW1=2.067E-03_r8*XNUSW1
    TCNDA1=9.050E-05_r8*XNUSA1
    WTHET1=1.467_r8-0.467_r8*THETPY(N3,N2,N1)
    TCND1=(STC(N3,N2,N1)+THETWX(N3,N2,N1)*TCNDW1 &
      +0.611_r8*THETIX(N3,N2,N1)*7.844E-03_r8 &
      +WTHET1*THETPX(N3,N2,N1)*TCNDA1) &
      /(DTC(N3,N2,N1)+THETWX(N3,N2,N1)+0.611_r8*THETIX(N3,N2,N1) &
      +WTHET1*THETPX(N3,N2,N1))
  ELSE
    TCND1=0.0_r8
  ENDIF
  end subroutine CalcThermConduc  

end module SurfPhysMod