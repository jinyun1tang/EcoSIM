module SnowPhysMod
!
! Description
! the snow model
! required input
!

! codes for snow physics
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils   , only : endrun   
  use SnowDataType
  use GridDataType
  use MiniMathMod
  use ClimForcDataType
  use EcoSimConst
  use SnowPhysData
  use EcoSIMSolverPar
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use HydroThermData
  use SurfLitterDataType
  use SurfSoilDataType  
  use SoilHeatDataType
  use SoilPropertyDataType
  use SoilPhysDataType
  use ChemTranspDataType
  use LandSurfDataType
  use PhysPars  
implicit none
  private
  character(len=*), parameter :: mod_filename=__FILE__
  
  public :: InitSnowLayers
  public :: SnowRedistribution
  public :: SnowBNDResistance
  public :: ZeroSnowFlux  
  public :: PrepIterSnowLayer
  public :: InitSnowAccums
  public :: CopySnowStates
  public :: SolveSnowpack
  public :: SumSnowDriftByRunoff
  public :: UpdateSnowAtM
  public :: UpdateSnowPack1
contains
!------------------------------------------------------------------------------------------
  subroutine InitSnowLayers(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: DLYRSI
  real(r8) :: VOLSWI
  real(r8), parameter :: cumSnowDepthI(JS)=(/0.05_r8,0.15_r8,0.30_r8,0.60_r8,1.00_r8/)
  !the maximum snow layer is 1.0 m.
  integer :: L
! begin_execution
!
! cumSnowDepth=depth to bottom
! DENSICE=ice density, 0.92 Mg m-3
! NewSnowDens=snow density (Mg m-3)
! VOLSS,VcumWatSnow,VOLIS,VOLS=snow,water,ice,total snowpack volume(m3)

! cumSnowDepthI=depth to bottom of snowpack layers (m), i.e. from current layere to surface
! DLYRS=snowpack layer thickness (m)
! VOLSSL,VcumWatSnowL,VOLISL,VLSnoDWI=snow,water,ice,total layer volume(m3), water equivalent snow 
! DENSS=layer density (Mg m-3)
! TKW,TCSnow=later temperature K,oC
! VHCPW=layer volumetric heat capacity (MJ m-3 K-1)
! SnowDepth=total snow height in the column

  cumSnowDepth(0,NY,NX)=0.0_r8
  NewSnowDens(NY,NX)=0.10_r8
  VcumDrySnoWE(NY,NX)=SnowDepth(NY,NX)*NewSnowDens(NY,NX)*DH(NY,NX)*DV(NY,NX)
  VcumWatSnow(NY,NX)=0.0_r8
  VcumIceSnow(NY,NX)=0.0_r8
  VcumSnoDWI(NY,NX)=VcumDrySnoWE(NY,NX)/NewSnowDens(NY,NX)+VcumWatSnow(NY,NX)+VcumIceSnow(NY,NX)

!  VOLSWI=0.0_r8

  !build the snow profile, topdown
  D9580: DO L=1,JS
    IF(L.EQ.1)THEN
      !bottom snow layer
      DLYRSI=cumSnowDepthI(L)
      SnowLayerThick(L,NY,NX)=AMIN1(DLYRSI,SnowDepth(NY,NX))
    ELSE
      DLYRSI=cumSnowDepthI(L)-cumSnowDepthI(L-1)
      SnowLayerThick(L,NY,NX)=AMIN1(DLYRSI,AZMAX1(SnowDepth(NY,NX)-cumSnowDepthI(L-1)))
    ENDIF
    VLDrySnoWE(L,NY,NX)=SnowLayerThick(L,NY,NX)*NewSnowDens(NY,NX)*DH(NY,NX)*DV(NY,NX)
    VLWatSnow(L,NY,NX)=0.0_r8
    VLIceSnow(L,NY,NX)=0.0_r8

!    IF(L.EQ.1)THEN
!      VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX)+VLIceSnow(L,NY,NX)*DENSICE)
!    ELSE
!      VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE(L-1,NY,NX)+VLWatSnow(L-1,NY,NX) &
!        +VLIceSnow(L-1,NY,NX)*DENSICE+VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX) &
!        +VLIceSnow(L,NY,NX)*DENSICE)
!    ENDIF

    SnoDensL(L,NY,NX)=NewSnowDens(NY,NX)
    VLSnoDWI(L,NY,NX)=VLDrySnoWE(L,NY,NX)/SnoDensL(L,NY,NX)+VLWatSnow(L,NY,NX)+VLIceSnow(L,NY,NX)
    VLSnowt0(L,NY,NX)=DLYRSI*DH(NY,NX)*DV(NY,NX)      !it is a non-zero number, potential/maximum volume
    cumSnowDepth(L,NY,NX)=cumSnowDepth(L-1,NY,NX)+SnowLayerThick(L,NY,NX)
    TKSnow(L,NY,NX)=AMIN1(Tref,TairKClimMean(NY,NX))
    TCSnow(L,NY,NX)=AZMIN1(ATCA(NY,NX))
    VLHeatCapSnow(L,NY,NX)=cps*VLDrySnoWE(L,NY,NX)+cpw*VLWatSnow(L,NY,NX)+cpi*VLIceSnow(L,NY,NX)
  ENDDO D9580

!
!     VLHeatCapSnowMin,=minimum heat capacities for solving
!      snowpack water and heat fluxes
!
  VLHeatCapSnowMin(NY,NX)=VLHeatCapSnoMin*AREA(3,NU(NY,NX),NY,NX)

  end subroutine InitSnowLayers

!------------------------------------------------------------------------------------------

  subroutine InitSnowAccums(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
!  TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
!  CumHeat2SnowLay=convective heat fluxes of snow,water,ice in snowpack
!

  D9875: DO L=1,JS
    CumSno2SnowLay(L,NY,NX)=0.0_r8
    CumWat2SnowLay(L,NY,NX)=0.0_r8
    CumIce2SnowLay(L,NY,NX)=0.0_r8
    CumHeat2SnowLay(L,NY,NX)=0.0_r8
    XSnowThawMassL(L,NY,NX)=0._r8
    XIceThawMassL(L,NY,NX)=0._r8    
  ENDDO D9875

  end subroutine InitSnowAccums


!------------------------------------------------------------------------------------------

  subroutine SnowPackIterationM(M,NY,NX,TotWatXFlx2SoiMicP,TotHeatFlow2Soi,WatFlowSno2MacP,&
    TotWatFlow2LitrByWat,TotHeatFlow2LitrByWat,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,&
    CumWatXFlx2SoiMicP,CumWatFlow2LitR,CumHeatFlow2LitR,cumHeatFlowSno2Soi)
  implicit none
  integer,  intent(in) :: M,NY,NX  
  real(r8), intent(out) :: TotWatXFlx2SoiMicP,TotHeatFlow2Soi,WatFlowSno2MacP
  real(r8), intent(out) :: TotWatFlow2LitrByWat,TotHeatFlow2LitrByWat
  real(r8), intent(inout) :: CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP
  real(r8), intent(inout) :: CumWatFlow2LitR,CumHeatFlow2LitR,cumHeatFlowSno2Soi
  real(r8) :: TCND1W,TCNDR
  real(r8) :: ATCNDW,VapCondSnoWeited,TCNDS
  real(r8) :: H2OVapFlx,H2OVapFlxMax,CumVapFlxLitr2Soi
  real(r8) :: CumVapFlxSno2Litr,CumHeatConvFlxSno2Litr
  real(r8) :: VapConvFlxInSnow,HeatbyVapConvInSnow
  real(r8) :: cumHeatConvFlxLitr2Soi1,CumHeatCndFlxSno2Litr
  real(r8) :: HeatCnduct,cumHeatCndFlxLitr2Soi,HeatCnductMax
  real(r8) :: PSISV1,VLairSno1,VLairSno2,VPY
  real(r8) :: TKY,FracAsAirSno1,FracAsAirInSno
  real(r8) :: FCDX,LOGFCX,FCX,VapFlxSno2Soi1,HeatConvFlxSno2Soi1
  real(r8) :: VapCond1,VapCond2,CNVR,DENSW1,DENSW2
  real(r8) :: WatVapFloInSnow,HeatByWatVapFloInSnow,WatFlowSno2Soil,WatFlowSno2MicP,PtWatFlowSno2Soi
  real(r8) :: WatFloInSno,HeatFlxByWatFloInSno,WatFlowSno2LitR
  real(r8) :: HeatFlowSno2LitrByWat,WatFloInSnoMax,HeatFlowSno2SoiByWat
  real(r8) :: HeatCndFlxSno2Soi,HeatCndFlxInSno,PSDX,TCND2W,THETRR
  real(r8) :: TK0X,TK1X,VapSnoSrc,VapSnoDest,WPX
  integer :: L,L2,ICHKL

  ! begin_execution
  ! PHYSICAL AND HYDRAULIC PROPERTIES OF SNOWPACK INCLUDING
  ! AIR AND WATER-FILLED POROSITY, WATER POTENTIAL OF UNDERLYING
  ! SOIL SURFACE USED IN FLUX CALCULATIONS
  !
  ! VHCPW,VLHeatCapSnowMin=current, minimum snowpack heat capacities
  ! VLDrySnoWE0M,VOLI0M,VOLW0M,VLSnoDWI1=snow,ice,water,total snowpack volume
  ! DENSS,DENSICE,NewSnowDens=snow,ice,minimum snow density
  ! AREA=area of grid cell
  ! VLairSno1=snowpack air volume
  ! FracAsAirSno1=snowpack air concentration
  ! VapCond1=snowpack vapor conductivity
  ! VapSnoSrc=snowpack vapor concentration
  ! TKSnow1=snowpack temperature
  ! WGSGW=vapor diffusivity
  ! DENSW1=snowpack density
  ! VHCPWMM=previous snowpack heat capacity

  ICHKL=0
  !loop from surface to bottom
  D9880: DO L=1,JS

    IF(VLHeatCapSnowM1(L,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
      ! active snow layer
      VLSnoDWI1(L,NY,NX)=VLDrySnoWE0M(L,NY,NX)/SnoDensL(L,NY,NX)+VLWatSnow0M(L,NY,NX) &
        +VLIceSnow0M(L,NY,NX)
      SnowLayerThick0(L,NY,NX)=VLSnoDWI1(L,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
      VLairSno1=AZMAX1(VLSnoDWI1(L,NY,NX)-VLDrySnoWE0M(L,NY,NX)-VLIceSnow0M(L,NY,NX) &
        -VLWatSnow0M(L,NY,NX))
      FracAsAirSno1=AMAX1(THETPI,VLairSno1/VLSnoDWI1(L,NY,NX))
      VapCond1=FracAsAirSno1**2.0_r8*H2OVapDifscSno(L,NY,NX)
      VapSnoSrc=vapsat(TKSnow1(L,NY,NX))  

      IF(VLSnoDWI1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        !maximum snow density is 0.6 g/cm3
        DENSW1=AMIN1(0.6_r8,(VLDrySnoWE0M(L,NY,NX)+VLWatSnow0M(L,NY,NX)+&
          VLIceSnow0M(L,NY,NX)*DENSICE)/VLSnoDWI1(L,NY,NX))
      ELSE
        DENSW1=NewSnowDens(NY,NX)
      ENDIF
      !
      ! SNOW THERMAL CONDUCTIVITY FROM J GLACIOL 43:26-41
      !
      ! TCND1W=snow thermal conductivity
      ! WatFloInSnoMax=porosity-unconstrained snow water flux
      !
      TCND1W=0.0036_r8*10._r8**(2.650_r8*DENSW1-1.652_r8)
      !
      ! DISCHARGE OF MELTWATER AND ITS HEAT FROM SNOWPACK LAYER
      ! TO LOWER SNOWPACK LAYER
      ! the value 0.05 below is arbitrary, 
      WatFloInSnoMax=AZMAX1(AZMAX1(VLWatSnow0M(L,NY,NX))-0.05_r8*AZMAX1(VLDrySnoWE0M(L,NY,NX)))*dts_sno
      
!      if(M>=28 .and. NY==3 .and. L==1)THEN
!        write(*,*)'WatFloInSnoMax',WatFloInSnoMax,VLWatSnow0M(L,NY,NX)
!      endif
      !
      ! WATER AND HEAT FLUXES IN SNOWPACK
      !
      ! SnowLayerThick0=snow layer thickness
      ! WatFloInSno=porosity-constrained snow water flux
      ! HeatFlxByWatFloInSno=convective heat flux from water flux
      ! id of next/target snow layer
      ! maximum JS snow layers 
      ! next layer is L2
      L2=MIN(JS,L+1)

      !not bottom layer, and it is heat significant
      IF(L.LT.JS.AND.VLHeatCapSnowM1(L2,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
        !if L==JS-1, L2==JS, so top layer is treated here.

        VLSnoDWI1(L2,NY,NX)=VLDrySnoWE0M(L2,NY,NX)/SnoDensL(L2,NY,NX)+VLWatSnow0M(L2,NY,NX) &
          +VLIceSnow0M(L2,NY,NX)
        SnowLayerThick0(L2,NY,NX)=VLSnoDWI1(L2,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
        VLairSno2=VLSnoDWI1(L2,NY,NX)-VLDrySnoWE0M(L2,NY,NX)-VLIceSnow0M(L2,NY,NX)-VLWatSnow0M(L2,NY,NX)        
        FracAsAirInSno=AMAX1(THETPI,VLairSno2/VLSnoDWI1(L2,NY,NX))
        !water flows only in air-filled pores
        WatFloInSno=AMIN1(FracAsAirInSno,WatFloInSnoMax)
        HeatFlxByWatFloInSno=cpw*TKSnow1(L,NY,NX)*WatFloInSno
        !
        ! VAPOR FLUX IN SNOWPACK
        !
        ! VLairSno1,VLairSno2=air-filled volumes of source, destination layers
        ! L2=destination layer
        ! VapCond1,VapCond2=vapor conductivities of source, destination layers
        ! VapSnoSrc,VapSnoDest=vapor concentrations of source, destination layers
        ! TKSnow1=soil temperature
        ! VapCondSnoWeited=snow vapor conductance
        ! SnowLayerThick0=snow layer thickness
        ! H2OVapFlx,H2OVapFlxMax=vapor-unconstrained,vapor-constrained vapor flux
        ! VapConvFlxInSnow,HeatbyVapConvInSnow=vapor flux and its convective heat flux
        !
        IF(VLairSno1.GT.ZEROS2(NY,NX).AND.VLairSno2.GT.ZEROS2(NY,NX))THEN
          !both layer has air-filled pores
          VapCond2=FracAsAirInSno**2.0_r8*H2OVapDifscSno(L2,NY,NX)
          VapSnoDest=vapsat(TKSnow1(L2,NY,NX))
          VapCondSnoWeited=2.0_r8*VapCond1*VapCond2/(VapCond1*SnowLayerThick0(L2,NY,NX) &
            +VapCond2*SnowLayerThick0(L,NY,NX))
          H2OVapFlx=VapCondSnoWeited*(VapSnoSrc-VapSnoDest)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX)*dts_snohttp
          VPY=(VapSnoSrc*VLairSno1+VapSnoDest*VLairSno2)/(VLairSno1+VLairSno2)
          !H2OVapFlxMax>0, moves out from layer L
          H2OVapFlxMax=(VapSnoSrc-VPY)*VLairSno1*dts_sno
          IF(H2OVapFlx.GE.0.0_r8)THEN
            VapConvFlxInSnow=AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax,VLWatSnow0M(L,NY,NX)*dts_wat))
            HeatbyVapConvInSnow=(cpw*TKSnow1(L,NY,NX)+EvapLHTC)*VapConvFlxInSnow
          ELSE
            VapConvFlxInSnow=AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax,-VLWatSnow0M(L2,NY,NX)*dts_wat))
            HeatbyVapConvInSnow=(cpw*TKSnow1(L2,NY,NX)+EvapLHTC)*VapConvFlxInSnow
          ENDIF
        ELSE
          VapConvFlxInSnow=0.0_r8
          HeatbyVapConvInSnow=0.0_r8
        ENDIF
        !
        ! HEAT FLUX IN SNOWPACK
        !
        ! DENSW2,TCNDW2=density,thermal conductivity in destination layer
        ! ATCNDW=thermal conductance
        ! SnowLayerThick0=layer thickness
        ! TKY=equilibrium temperature
        ! HeatCnductMax,HeatCnduct=heat-constrained,heat-unconstrained heat fluxes
        ! VHCPWMM,TKSnow1=volumetric heat capacity,temperature
        ! dts_wat=time step for flux calculations
        ! FSNW=snow cover fraction
        ! dts_snohttp=time step for snowpack flux calculations
        ! HeatCndFlxInSno=snowpack heat flux
        ! FLW0S,FLQ0I,FLQ0W=snow,ice,water fluxes through snowpack
        ! HFLW0W=convective heat flux snow,water,ice fluxes
        !
        IF(VLSnoDWI1(L2,NY,NX).GT.ZEROS2(NY,NX))THEN
          DENSW2=AMIN1(0.6_r8,(VLDrySnoWE0M(L2,NY,NX)+VLWatSnow0M(L2,NY,NX) &
            +VLIceSnow0M(L2,NY,NX)*DENSICE)/VLSnoDWI1(L2,NY,NX))
        ELSE
          !initialize new snow layer
          DENSW2=NewSnowDens(NY,NX)  
        ENDIF
        TCND2W=0.0036_r8*10._r8**(2.650_r8*DENSW2-1.652_r8)
        ATCNDW=2.0_r8*TCND1W*TCND2W/(TCND1W*SnowLayerThick0(L2,NY,NX)+TCND2W*SnowLayerThick0(L,NY,NX))
        TKY=(TKSnow1(L,NY,NX)*VLHeatCapSnowM1(L,NY,NX)+TKSnow1(L2,NY,NX) &
          *VLHeatCapSnowM1(L2,NY,NX))/(VLHeatCapSnowM1(L,NY,NX)+VLHeatCapSnowM1(L2,NY,NX))
        HeatCnductMax=(TKSnow1(L,NY,NX)-TKY)*VLHeatCapSnowM1(L,NY,NX)*dts_sno
        HeatCnduct=ATCNDW*(TKSnow1(L,NY,NX)-TKSnow1(L2,NY,NX))*AREA(3,NUM(NY,NX),NY,NX) &
          *FracSurfAsSnow(NY,NX)*dts_snohttp

        IF(HeatCnduct.GE.0.0_r8)THEN
          HeatCndFlxInSno=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
        ELSE
          HeatCndFlxInSno=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
        ENDIF

        WatVapFloInSnow=WatFloInSno+VapConvFlxInSnow
        HeatByWatVapFloInSnow=HeatFlxByWatFloInSno+HeatbyVapConvInSnow+HeatCndFlxInSno
        SnoX2SnoLay(L2,NY,NX)=0.0_r8
        WatX2SnoLay(L2,NY,NX)=WatVapFloInSnow
        IceX2SnoLay(L2,NY,NX)=0.0_r8
        HeatX2SnoLay(L2,NY,NX)=HeatByWatVapFloInSnow
        WatFlowInSnowM(M,L2,NY,NX)=WatFlowInSnowM(M,L2,NY,NX)+WatFloInSno
        !
        ! DISCHARGE OF MELTWATER AND ITS HEAT FROM LOWEST SNOWPACK LAYER
        ! TO RESIDUE, SURFACE SOIL MICROPORES AND MACROPORES
        !
        ! WatFloInSnoMax,WatFlowSno2LitR=porosity-unconstrained water flux to soil,litter
        ! PtWatFlowSno2Soi,WatFlowSno2MicP,WatFlowSno2MacP=water flux to soil surface,
        ! micropores,macropores
        ! VOLP1,VOLPH1=air volumes of soil micropores,macropores
        ! SoilFracAsMacP1,FGRD=macropore,micropore volume fractions
        ! HeatFlowSno2SoiByWat,HeatFlowSno2LitrByWat=convective heat fluxes to soil,litter

        ! THETY=hygroscopic water concentration
        ! POROS=soil porosity
        ! FC,WP,FCL,WPL=field capacity,wilting point, log(FC),log(WP)
        ! FCI,WPI=FC,WP of ice
        ! FracSoiPAsIce=ice concentration
        ! BKVL=bulk density x volume of soil layer
        ! VLitR=surface litter volume
        ! FracSoiPAsWat=relative pore fraciton filled by water
      ELSE
        !dealing with residual snow
        !L==JS, or the layer L2 has insignificant snow, so layer L is directly interacting with litter and soil surface        
        IF(ICHKL.EQ.0)THEN

          !interaction with respect to litter layer and topsoil 
          !potential flow to soil
          PtWatFlowSno2Soi=WatFloInSnoMax*FracSurfAsBareSoi(NY,NX)
          !actual flow to micro- and macropores
          WatFlowSno2MicP=AMIN1(VLairMicP1(NUM(NY,NX),NY,NX)*dts_wat,PtWatFlowSno2Soi*SoilFracAsMicP(NUM(NY,NX),NY,NX))
          WatFlowSno2MacP=AMIN1(VLairMacP1(NUM(NY,NX),NY,NX)*dts_wat,PtWatFlowSno2Soi*SoilFracAsMacP1(NUM(NY,NX),NY,NX))
          WatFlowSno2Soil=WatFlowSno2MicP+WatFlowSno2MacP
          HeatFlowSno2SoiByWat=cpw*TKSnow1(L,NY,NX)*WatFlowSno2Soil
          WatFlowSno2LitR=WatFloInSnoMax-WatFlowSno2Soil
          HeatFlowSno2LitrByWat=cpw*TKSnow1(L,NY,NX)*WatFlowSno2LitR

          call SnowTopSoilExch(M,L,NY,NX,VapCond1,VapSnoSrc,VLairSno1,TCND1W,&
            VapFlxSno2Soi1,HeatConvFlxSno2Soi1,HeatCndFlxSno2Soi,VapCond2,PSISV1,TCNDS)

          !
          ! HEAT FLUX AMONG SNOWPACK, SURFACE RESIDUE AND SURFACE SOIL
          !
          ! CumVapFlxSno2Litr=snowpack-litter vapor flux
          ! CumHeatConvFlxSno2Litr,CumHeatCndFlxSno2Litr=snowpack-litter convective,conductive heat fluxes
          ! VapFlxSno2Soi1=snowpack-soil vapor flux
          ! HeatConvFlxSno2Soi1,HeatCndFlxSno2Soi=snowpack-soil convective,conductive heat fluxes
          ! VLHeatCapacity,VHeatCapLitR=current,minimum litter heat capacities
          ! TK0X,TKXR,TK1X=snowpack,litter,soil temperatures
          ! CNVR,VapCond1,VapCond2=litter,snowpack,soil vapor conductivity
          ! THETP*,FracSoiPAsWat,FracSoiPAsIce=litter air,water,ice concentration
          ! POROS,POROQ=litter porosity, tortuosity
          ! CVRD=litter cover fraction
          ! VaporDiffusivityLitR=litter vapor diffusivity
          ! AvgVaporCondctSnowLitR,AvgVaporCondctSoilLitR=snowpack-litter,litter-soil vapor conductance
          ! DLYRR,SnowLayerThick0,DLYR=litter,snowpack,soil depths
          ! THETRR=dry litter concentration
          ! TCNDR,TCND1W,TCNDS=litter,snowpack,soil thermal conductivity
          ! AvgThermCondctSnoLitR,AvgThermCondctSoilLitR=snow-litter,litter-soil thermal conductance
          !
          CumVapFlxSno2Litr=0.0_r8
          CumHeatConvFlxSno2Litr=0.0_r8
          CumHeatCndFlxSno2Litr=0.0_r8
          CumVapFlxLitr2Soi=0.0_r8
          cumHeatConvFlxLitr2Soi1=0.0_r8
          cumHeatCndFlxLitr2Soi=0.0_r8

          !surface litter layer is active
          IF(VLHeatCapacity(0,NY,NX).GT.VHeatCapLitR(NY,NX))THEN
            call SnowSurLitterExch(M,L,NY,NX,VapCond1,VapCond2,TCND1W,VLairSno1,PSISV1,TCNDS,&
              CumVapFlxSno2Litr,CumHeatConvFlxSno2Litr,CumHeatCndFlxSno2Litr,CumVapFlxLitr2Soi,&
              cumHeatConvFlxLitr2Soi1,cumHeatCndFlxLitr2Soi)
          ENDIF
          !
          ! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
          ! FOR LATER UPDATES TO STATE VARIABLES
          !
          ! TotWatXFlx2SoiMicP,CumWatFlx2SoiMicP=total,accumulated water flux to soil micropores
          ! CumWatXFlx2SoiMicP,CumWatFlx2SoiMacP=total,accumd snow-soil micropore,macropore water
          ! TotHeatFlow2Soi,cumHeatFlowSno2Soi=total,accumulated snow+litter heat flux to soil
          ! TotWatFlow2LitrByWat,CumWatFlow2LitR=total,accumulated snow+soil water flux to litter
          ! TotHeatFlow2LitrByWat,CumHeatFlow2LitR=total,accumulated snow+soil heat flux to litter
          ! WatFlowSno2LitRM,WatFlowSno2MicPM,WatFlowSno2MacPM=total water flux to litter,soil micropore,macropore
          ! FLSW,WatConvSno2MacP,WatConvSno2LitR=water flux from lowest snow layer to soil macropore,micropore,litter
          ! HeatConvSno2Soi,HeatConvSno2LitR=heat flux from lowest snow layer to soil,litter
!
          TotWatXFlx2SoiMicP=WatFlowSno2MicP+VapFlxSno2Soi1+CumVapFlxLitr2Soi
          CumWatFlx2SoiMicP=CumWatFlx2SoiMicP+TotWatXFlx2SoiMicP
          if(abs(CumWatFlx2SoiMicP)>1.e20_r8)then
            write(*,*)'CumWatFlx2SoiMicP=',WatFlowSno2MicP,VapFlxSno2Soi1,CumVapFlxLitr2Soi
            write(*,*)'at line',__LINE__
            call endrun(trim(mod_filename)//'at line',__LINE__)
          endif
          CumWatXFlx2SoiMicP=CumWatXFlx2SoiMicP+WatFlowSno2MicP
          CumWatFlx2SoiMacP=CumWatFlx2SoiMacP+WatFlowSno2MacP
          TotHeatFlow2Soi=HeatFlowSno2SoiByWat+HeatConvFlxSno2Soi1+HeatCndFlxSno2Soi &
            +cumHeatConvFlxLitr2Soi1+cumHeatCndFlxLitr2Soi
          cumHeatFlowSno2Soi=cumHeatFlowSno2Soi+TotHeatFlow2Soi
          TotWatFlow2LitrByWat=WatFlowSno2LitR+CumVapFlxSno2Litr-CumVapFlxLitr2Soi
!          if(M>=28 .and. NY==3 .and. L==1)THEN
!            write(*,*)'TotWatFlow2LitrByWat',TotWatFlow2LitrByWat,WatFlowSno2LitR,CumVapFlxSno2Litr,CumVapFlxLitr2Soi
!          endif
          CumWatFlow2LitR=CumWatFlow2LitR+TotWatFlow2LitrByWat
          TotHeatFlow2LitrByWat=HeatFlowSno2LitrByWat+CumHeatConvFlxSno2Litr+CumHeatCndFlxSno2Litr &
            -cumHeatConvFlxLitr2Soi1-cumHeatCndFlxLitr2Soi
          CumHeatFlow2LitR=CumHeatFlow2LitR+TotHeatFlow2LitrByWat

          WatFlowSno2LitRM(M,NY,NX)=WatFlowSno2LitRM(M,NY,NX)+WatFlowSno2LitR
          WatFlowSno2MicPM(M,NY,NX)=WatFlowSno2MicPM(M,NY,NX)+WatFlowSno2MicP
          WatFlowSno2MacPM(M,NY,NX)=WatFlowSno2MacPM(M,NY,NX)+WatFlowSno2MacP
          FLSW(L,NY,NX)=FLSW(L,NY,NX)+TotWatXFlx2SoiMicP
          WatConvSno2MacP(L,NY,NX)=WatConvSno2MacP(L,NY,NX)+WatFlowSno2MacP
          HeatConvSno2Soi(L,NY,NX)=HeatConvSno2Soi(L,NY,NX)+TotHeatFlow2Soi
          WatConvSno2LitR(L,NY,NX)=WatConvSno2LitR(L,NY,NX)+TotWatFlow2LitrByWat
          HeatConvSno2LitR(L,NY,NX)=HeatConvSno2LitR(L,NY,NX)+TotHeatFlow2LitrByWat

          ICHKL=1
        ENDIF
      ENDIF
    ENDIF
  ENDDO D9880
  end subroutine SnowPackIterationM

!------------------------------------------------------------------------------------------

  subroutine SolveSnowpack(M,NY,NX,LatentHeatAir2Sno,Radnet2Snow,HeatSensEvap,HeatSensAir2Snow,&
    HeatNetFlx2Snow,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP,CumWatFlow2LitR,&
    CumHeatFlow2LitR,cumHeatFlowSno2Soi)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(inout) :: HeatSensEvap
  real(r8), intent(inout) :: CumWatFlow2LitR
  real(r8), intent(out) :: HeatNetFlx2Snow
  real(r8), intent(out) :: LatentHeatAir2Sno
  real(r8), intent(out) :: Radnet2Snow
  real(r8), intent(out) :: HeatSensAir2Snow
  real(r8), intent(out) :: CumHeatFlow2LitR
  real(r8), intent(out) :: cumHeatFlowSno2Soi
  real(r8), intent(out) :: CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP
  integer :: MM,L,L2

  real(r8) :: TotWatXFlx2SoiMicP,TotHeatFlow2Soi  
  real(r8) :: ENGY0,FVOLI0,FVOLS0
  real(r8) :: TKApp,TKXR,TK1X

  real(r8) :: WatFlowSno2MacP,TotWatFlow2LitrByWat,TotHeatFlow2LitrByWat
  real(r8) :: NetIce2LayL,NetSno2LayL,NetWat2LayL,HeatByFrezThaw,TFLX1

  real(r8) :: NetHeat2LayL,TotSnowLMass,VLHeatCapSnowMX
  real(r8) :: VOLI0X,VOLS0X,VOLW0X,IceThawMass,SnowThawMass
  !     begin_execution
  !     HEAT AND VAPOR FLUXES BETWEEN SNOWPACK AND ATMOSPHERE
  !
  !     VHCPWM=volumetric heat capacity of snowpack
  !     NPS=number of cycles for solving snowpack heat and water fluxes
  !     SnowAlbedo=snowpack albedo
  !     VLDrySnoWE0M,VOLI0M,VOLW0M=snow,ice,water volumes
  !     RFLX0=net radiation input
  !     RADXW=shortwave radiation at snowpack surface
  !     LWRad2Snow=longwave radn incident at snowpack surface
  !     LWRadSno1=longwave radn emitted by snowpack surface
  !     TKSnow1=snowpack surface temperature
  !     RadNet2Sno2=net radiation

  CumWatFlx2SoiMacP =0._r8
  cumHeatFlowSno2Soi=0._r8
  CumWatFlx2SoiMicP=0._r8
  CumHeatFlow2LitR=0.0_r8
  CumWatXFlx2SoiMicP=0._r8
  HeatSensAir2Snow=0.0_r8
  Radnet2Snow=0.0_r8
  LatentHeatAir2Sno=0.0_r8
  HeatNetFlx2Snow=0.0_r8

  D3000: DO MM=1,NPS

     if(TKSoi1(0,NY,NX)<100._r8 .or. TKSoi1(0,NY,NX)>400._r8)write(*,*)'TXKR MM=',MM,TKSoi1(0,NY,NX)
!    write(*,*)'CumHeatFlow2LitR MM=',MM,NY,NX,CumHeatFlow2LitR
    call SnowAtmosExchange(M,NY,NX,LatentHeatAir2Sno,HeatSensEvap,HeatNetFlx2Snow,Radnet2Snow,HeatSensAir2Snow)

    call SnowPackIterationM(M,NY,NX,TotWatXFlx2SoiMicP,TotHeatFlow2Soi,WatFlowSno2MacP,&
      TotWatFlow2LitrByWat,TotHeatFlow2LitrByWat,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,&
      CumWatXFlx2SoiMicP,CumWatFlow2LitR,CumHeatFlow2LitR,cumHeatFlowSno2Soi)
!
!     ACCUMULATE SNOWPACK FLUXES TO LONGER TIME STEP FOR
!     LITTER, SOIL FLUX CALCULATIONS
!
!     LWRadBySurf=total longwave emission
!     XFLWS,WatXfer2SnoLay,IceXfer2SnoLay=hourly accumulated snow,water,ice transfer
!     HeatXfer2SnoLay=hourly convective heat flux from snow,water,ice transfer
!     NetSno2LayL,NetWat2LayL,NetIce2LayL=net snow,water,ice transfer
!     NetHeat2LayL=convective heat flux from net snow,water,ice transfer
!     TFLWS,TFLWW,TFLWI=accumulated net snow,water,ice transfer
!     CumHeat2SnowLay=convective heat flux from accumd snow,water,ice transfer
!
    D9860: DO L=1,JS
      SnoXfer2SnoLay(L,NY,NX)=SnoXfer2SnoLay(L,NY,NX)+SnoX2SnoLay(L,NY,NX)
      WatXfer2SnoLay(L,NY,NX)=WatXfer2SnoLay(L,NY,NX)+WatX2SnoLay(L,NY,NX)
      IceXfer2SnoLay(L,NY,NX)=IceXfer2SnoLay(L,NY,NX)+IceX2SnoLay(L,NY,NX)
      HeatXfer2SnoLay(L,NY,NX)=HeatXfer2SnoLay(L,NY,NX)+HeatX2SnoLay(L,NY,NX)
      !next layer
      L2=MIN(JS,L+1)
!
!     IF WITHIN SNOWPACK
!
      IF(L.LT.JS.AND.VLHeatCapSnowM1(L2,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
        NetSno2LayL=SnoX2SnoLay(L,NY,NX)-SnoX2SnoLay(L2,NY,NX)
        NetWat2LayL=WatX2SnoLay(L,NY,NX)-WatX2SnoLay(L2,NY,NX)
        NetIce2LayL=IceX2SnoLay(L,NY,NX)-IceX2SnoLay(L2,NY,NX)
        NetHeat2LayL=HeatX2SnoLay(L,NY,NX)-HeatX2SnoLay(L2,NY,NX)
        CumSno2SnowLay(L,NY,NX)=CumSno2SnowLay(L,NY,NX)+NetSno2LayL
        CumWat2SnowLay(L,NY,NX)=CumWat2SnowLay(L,NY,NX)+NetWat2LayL
        CumIce2SnowLay(L,NY,NX)=CumIce2SnowLay(L,NY,NX)+NetIce2LayL
        CumHeat2SnowLay(L,NY,NX)=CumHeat2SnowLay(L,NY,NX)+NetHeat2LayL
!
!     IF AT BOTTOM OF SNOWPACK
!
      ELSEIF(VLHeatCapSnowM1(L,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
        NetSno2LayL=SnoX2SnoLay(L,NY,NX)
        NetWat2LayL=WatX2SnoLay(L,NY,NX)-TotWatFlow2LitrByWat-TotWatXFlx2SoiMicP-WatFlowSno2MacP
        NetIce2LayL=IceX2SnoLay(L,NY,NX)
        NetHeat2LayL=HeatX2SnoLay(L,NY,NX)-TotHeatFlow2LitrByWat-TotHeatFlow2Soi
        CumSno2SnowLay(L,NY,NX)=CumSno2SnowLay(L,NY,NX)+NetSno2LayL
        CumWat2SnowLay(L,NY,NX)=CumWat2SnowLay(L,NY,NX)+NetWat2LayL
        CumIce2SnowLay(L,NY,NX)=CumIce2SnowLay(L,NY,NX)+NetIce2LayL
        CumHeat2SnowLay(L,NY,NX)=CumHeat2SnowLay(L,NY,NX)+NetHeat2LayL
!        if(M>=28.and.NY==3.and.L==1)then
!          write(*,*)'NetWat2LayL',WatX2SnoLay(L,NY,NX),TotWatFlow2LitrByWat,TotWatXFlx2SoiMicP,WatFlowSno2MacP
!        endif
      ELSE
        NetSno2LayL=0.0_r8
        NetWat2LayL=0.0_r8
        NetIce2LayL=0.0_r8
        NetHeat2LayL=0.0_r8
      ENDIF
!
!     FREEZE-THAW IN SNOWPACK FROM NET CHANGE IN SNOWPACK
!     HEAT STORAGE
!
!     VLDrySnoWE0M,VOLW0M,VOLI0M=snow,water,ice volume
!     VHCPWMM,VLHeatCapSnowMX,VLHeatCapSnowMin=previous,current,minimum heat capacity
!     NetHeat2LayL=net conductive+convective heat flux
!     TFLX1=unconstrained latent heat flux from freeze-thaw
!     FVOLS0,FVOLI0=fractions of total water in water,ice
!     HeatByFrezThaw=source-limited latent heat flux from freeze-thaw
!     SnowThawMass,IceThawMass=freeze-thaw changes in water,ice
!     SnowThawMassL,IceThawMassL=accumulated freeze-thaw
!     PhaseChangeHeatL=accumulated latent heat flux from freeze-thaw
!     XSnowThawMassL,XIceThawMassL=hourly accumulated freeze-thaw
!     XPhaseChangeHeatL=hourly accumulated latent heat flux from freeze-thaw
!
      VOLS0X=AZMAX1(VLDrySnoWE0M(L,NY,NX))
      VOLW0X=AZMAX1(VLWatSnow0M(L,NY,NX))
      VOLI0X=AZMAX1(VLIceSnow0M(L,NY,NX))
      ENGY0=VLHeatCapSnowM1(L,NY,NX)*TKSnow1(L,NY,NX)
      VLHeatCapSnowMX=cps*VOLS0X+cpw*VOLW0X+cpi*VOLI0X

      IF(VLHeatCapSnowMX.GT.VLHeatCapSnowMin(NY,NX))THEN
        !apparent temperature before freeze-thaw
        TKApp=(ENGY0+NetHeat2LayL)/VLHeatCapSnowMX
        IF((TKApp.LT.TFice.AND.VOLW0X.GT.ZERO*VcumSnoDWI(NY,NX)) &
          .OR.(TKApp.GT.TFice.AND.VOLI0X+VOLS0X.GT.ZERO*VcumSnoDWI(NY,NX)))THEN
          !freeze-thaw condition met
          TFLX1=VLHeatCapSnowMX*(TFice-TKApp)/2.7185_r8*dts_wat

          IF(TFLX1.LT.0.0_r8)THEN
            !thaw
            TotSnowLMass=VOLS0X+VOLI0X*DENSICE
            IF(TotSnowLMass.GT.ZEROS2(NY,NX))THEN
              FVOLS0=VOLS0X/TotSnowLMass
              FVOLI0=VOLI0X*DENSICE/TotSnowLMass
            ELSE
              FVOLS0=0.0_r8
              FVOLI0=0.0_r8
            ENDIF
            HeatByFrezThaw=AMAX1(-333.0_r8*TotSnowLMass*dts_wat,TFLX1)
            SnowThawMass=-HeatByFrezThaw*FVOLS0/333.0_r8
            IceThawMass=-HeatByFrezThaw*FVOLI0/333.0_r8
          ELSE
            !freeze
            FVOLS0=0.0_r8
            FVOLI0=0.0_r8
            HeatByFrezThaw=AMIN1(333.0_r8*VOLW0X*dts_wat,TFLX1)
            SnowThawMass=0.0_r8
            IceThawMass=-HeatByFrezThaw/333.0_r8
          ENDIF
        ELSE
          TFLX1=0.0_r8
          FVOLS0=0.0_r8
          FVOLI0=0.0_r8
          HeatByFrezThaw=0.0_r8
          SnowThawMass=0.0_r8
          IceThawMass=0.0_r8
        ENDIF
        SnowThawMassL(L,NY,NX)=SnowThawMassL(L,NY,NX)+SnowThawMass
        IceThawMassL(L,NY,NX)=IceThawMassL(L,NY,NX)+IceThawMass
        PhaseChangeHeatL(L,NY,NX)=PhaseChangeHeatL(L,NY,NX)+HeatByFrezThaw
        XSnowThawMassL(L,NY,NX)=XSnowThawMassL(L,NY,NX)+SnowThawMass
        XIceThawMassL(L,NY,NX)=XIceThawMassL(L,NY,NX)+IceThawMass
        XPhaseChangeHeatL(L,NY,NX)=XPhaseChangeHeatL(L,NY,NX)+HeatByFrezThaw
      ELSE
        HeatByFrezThaw=0.0_r8
        SnowThawMass=0.0_r8
        IceThawMass=0.0_r8
      ENDIF
!
!     INTERNAL SNOWPACK SNOW, WATER, ICE, TEMPERATURE
!
!     VLDrySnoWE0M,VOLW0M,VOLI0M=snow water eqv,water,ice volume
!     NetSno2LayL,NetWat2LayL,NetIce2LayL=net snow,water,ice transfer
!     NetHeat2LayL=conductive+convective heat from snow,water,ice transfer
!     SnowThawMass,IceThawMass=freeze-thaw changes in water,ice
!     HeatByFrezThaw=source-limited latent heat flux from freeze-thaw
!     DENSICE=ice density
!     TKSnow1,TairK=snowpack,air temperature
!     VHCPWMM,VLHeatCapSnowMin=snowpack, minimum heat capacity
!
!      if(M>=28.and.NY==3.and.L==1)then
!        write(*,*)'SolveSnowpackbbbf  ',MM,VLDrySnoWE0M(L,NY,NX),VLWatSnow0M(L,NY,NX),VLIceSnow0M(L,NY,NX)
!      endif

      VLDrySnoWE0M(L,NY,NX)=VLDrySnoWE0M(L,NY,NX)+NetSno2LayL-SnowThawMass
      VLWatSnow0M(L,NY,NX)=VLWatSnow0M(L,NY,NX)+NetWat2LayL+SnowThawMass+IceThawMass      
      VLIceSnow0M(L,NY,NX)=VLIceSnow0M(L,NY,NX)-IceThawMass/DENSICE      
!      if(M>=28.and.NY==3.and.L==1)then
!        write(*,*)'SolveSnowpackaf  ',MM,VLDrySnoWE0M(L,NY,NX),VLWatSnow0M(L,NY,NX),VLIceSnow0M(L,NY,NX)
!        if(VLWatSnow0M(L,NY,NX)<0._r8)then
!          write(*,*)'flx',NetWat2LayL,SnowThawMass,IceThawMass      
!          call endrun(trim(mod_filename)//' at line',__LINE__)        
!        endif
!      endif
      ENGY0=VLHeatCapSnowM1(L,NY,NX)*TKSnow1(L,NY,NX)
      VLHeatCapSnowM1(L,NY,NX)=cps*VLDrySnoWE0M(L,NY,NX)+cpw*VLWatSnow0M(L,NY,NX)+cpi*VLIceSnow0M(L,NY,NX)
      IF(VLHeatCapSnowM1(L,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
        TKSnow1(L,NY,NX)=(ENGY0+NetHeat2LayL+HeatByFrezThaw)/VLHeatCapSnowM1(L,NY,NX)
      ELSEIF(L.EQ.1)THEN
        TKSnow1(L,NY,NX)=TairK(NY,NX)
      ELSE
        TKSnow1(L,NY,NX)=TKSnow1(L-1,NY,NX)
      ENDIF

    ENDDO D9860
  ENDDO D3000
  end subroutine SolveSnowpack  
!------------------------------------------------------------------------------------------
  subroutine SnowAtmosExchange(M,NY,NX,LatentHeatAir2Sno,HeatSensEvap,HeatNetFlx2Snow,Radnet2Snow,HeatSensAir2Snow)
  implicit none  
  integer, intent(in) :: M,NY,NX
  real(r8), intent(inout) :: LatentHeatAir2Sno,HeatSensEvap,HeatNetFlx2Snow,Radnet2Snow,HeatSensAir2Snow
  real(r8) :: SnowAlbedo,RFLX0,RI
  real(r8) :: LWRadSno1,RadNet2Sno2
  real(r8) :: RAGX
  real(r8):: Raa,PARE,PARS
  real(r8) :: VPSno0,EVAPW2,EVAPX2,HeatSensAir2Sno2
  real(r8) :: LatentHeatAir2Sno2,EvapSublimation2,MaxVapXAir2Sno
  real(r8) :: HeatNetFlx2Sno1,HeatNetFlx2Sno2,HeatSensAir2SnoByEvap2
  real(r8) :: NetHeatAir2Snow,FLW0W2
  real(r8) :: HeatSnofall2Snow,FLW0S2,FLW0I2
  real(r8) :: FLQ0I2,FLQ0S2,FLQ0W2

  SnowAlbedo=(0.85_r8*VLDrySnoWE0M(1,NY,NX)+0.30_r8*VLIceSnow0M(1,NY,NX)+0.06_r8*VLWatSnow0M(1,NY,NX)) &
    /(VLDrySnoWE0M(1,NY,NX)+VLIceSnow0M(1,NY,NX)+VLWatSnow0M(1,NY,NX))

  RFLX0=(1.0_r8-SnowAlbedo)*RADXW(NY,NX)+LWRad2Snow(NY,NX)    !incoming radiation, short + longwave
  LWRadSno1=THRMW(NY,NX)*TKSnow1(1,NY,NX)**4._r8         !emitting longwave radiation,
  RadNet2Sno2=RFLX0-LWRadSno1                            !net radiation
  !
  !     AERODYNAMIC RESISTANCE ABOVE SNOWPACK INCLUDING
  !     RESISTANCE IMPOSED BY PLANT CANOPY
  !
  !     RI=Richardsons number
  !     RIB=isothermal RI
  !     TKQ=canopy air temperature
  !     RAGX,RA=snowpack blr
  !     RAG,RAGW=isothermal blrs at ground,snowpack surfaces
  !
  RI=RichardsonNumber(RIB(NY,NX),TKQ(NY,NX),TKSnow1(1,NY,NX))
  RAGX=AMAX1(RAM,0.8_r8*RAGW(NY,NX),AMIN1(1.2_r8*RAGW(NY,NX),RAG(NY,NX)/(1.0_r8-10.0_r8*RI)))
  RAGW(NY,NX)=RAGX
  RAa=RAGX
  !
  ! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
  !
  !     PARE,PARS=blcs for snowpack latent,sensible heat fluxes
  !     PAREW,PARSW=conductances for latent,sensible heat fluxes
  !     RZ=surface resistance
  !     VPSno0,VPQ=vapor pressure at snowpack surface, canopy air
  !     MaxVapXAir2Sno,EVAPW2,EvapSublimation2=evaporation total, water,snow
  !     XNPS=1/NPS
  !     LatentHeatAir2Sno2=latent heat flux
  !     VAP,VAPS=latent heat of evaporation,sublimation
  !     HeatSensAir2SnoByEvap2=convective heat of evaporation flux
  !
  PARE=PAREW(NY,NX)/(RAa+RZ)
  PARS=PARSW(NY,NX)/RAa    
  VPSno0=vapsat(TKSnow1(1,NY,NX))

  MaxVapXAir2Sno=PARE*(VPQ(NY,NX)-VPSno0)  
  !first the loss is evaporation from snow held water
  EVAPW2=AMAX1(MaxVapXAir2Sno,-AZMAX1(VLWatSnow0M(1,NY,NX)*dts_sno))    
  !then the loss is sublimation from dry snow
  EVAPX2=AZMIN1(MaxVapXAir2Sno-EVAPW2)
  !then the loss is from dry snow
  EvapSublimation2=AMAX1(EVAPX2,-AZMAX1(VLDrySnoWE0M(1,NY,NX)*dts_sno))
  LatentHeatAir2Sno2=EVAPW2*EvapLHTC+EvapSublimation2*SublmHTC

  IF(MaxVapXAir2Sno.LT.0.0_r8)THEN
    !snow is losing water/heat
    HeatSensAir2SnoByEvap2=(EVAPW2*cpw+EvapSublimation2*cps)*TKSnow1(1,NY,NX)
  ELSE
    !snow is gaining water/heat, condensation/deposition
    HeatSensAir2SnoByEvap2=(EVAPW2*cpw+EvapSublimation2*cps)*TKQ(NY,NX)
  ENDIF
!
!     SOLVE FOR SNOWPACK SURFACE TEMPERATURE AT WHICH ENERGY
!     BALANCE OCCURS, SOLVE AND ACCUMULATE LATENT, SENSIBLE
!     STORAGE HEAT FLUXES AND EVAPORATION
!
!     HeatSensAir2Sno2,LatentHeatAir2Sno2,RadNet2Sno2=sensible,latent heat fluxes, net radiation
!     HeatSensAir2SnoByEvap2=convective heat flux from LatentHeatAir2Sno2
!     HeatNetFlx2Sno1=storage heat flux
!     FLQ0S2,FLQ0W2,FLQ0I2=snow,water,ice input to snowpack
!     HeatSnofall2Snow=convective heat from snow,water,ice input to snowpack
!
  HeatSensAir2Sno2=PARS*(TKQ(NY,NX)-TKSnow1(1,NY,NX))
  !occasionally, RadNet2Sno2 and HeatSensAir2Sno2 go to infinity
  HeatNetFlx2Sno1=RadNet2Sno2+LatentHeatAir2Sno2+HeatSensAir2Sno2   
  HeatNetFlx2Sno2=HeatNetFlx2Sno1+HeatSensAir2SnoByEvap2
  Radnet2Snow=Radnet2Snow+RadNet2Sno2
  LatentHeatAir2Sno=LatentHeatAir2Sno+LatentHeatAir2Sno2
  HeatSensEvap=HeatSensEvap+HeatSensAir2SnoByEvap2
  HeatSensAir2Snow=HeatSensAir2Snow+HeatSensAir2Sno2
  HeatNetFlx2Snow=HeatNetFlx2Snow+HeatNetFlx2Sno2

  EVAPS(NY,NX)=EVAPS(NY,NX)+EvapSublimation2
  EVAPW(NY,NX)=EVAPW(NY,NX)+EVAPW2
  VapXAir2Sno(NY,NX)=VapXAir2Sno(NY,NX)+EvapSublimation2+EVAPW2
  FLQ0S2=SnowFallt(NY,NX)*XNPS
  FLQ0W2=Rain2Snowt(NY,NX)*XNPS
  FLQ0I2=Ice2Snowt(NY,NX)*XNPS
  HeatSnofall2Snow=HeatFall2Snowt(NY,NX)*XNPS
  FLW0S2=FLQ0S2+EvapSublimation2
  FLW0W2=FLQ0W2+EVAPW2
  FLW0I2=FLQ0I2
  NetHeatAir2Snow=HeatSnofall2Snow+HeatNetFlx2Sno2
  SnoX2SnoLay(1,NY,NX)=FLW0S2
  WatX2SnoLay(1,NY,NX)=FLW0W2
  IceX2SnoLay(1,NY,NX)=FLW0I2
  HeatX2SnoLay(1,NY,NX)=NetHeatAir2Snow
  WatFlowInSnowM(M,1,NY,NX)=WatFlowInSnowM(M,1,NY,NX)+FLQ0S2+FLQ0I2+FLQ0W2
  LWRadBySurf(NY,NX)=LWRadBySurf(NY,NX)+LWRadSno1
!     IF(NX.EQ.3.AND.NY.EQ.3)THEN
!     WRITE(*,7759)'EVAP',I,J,M,MM,FLW0S2
!    2,FLQ0S2,EvapSublimation2,FLW0W2,FLQ0W2
!    3,FracSurfAsSnow(NY,NX),FLW0I2,FLQ0I2,RadNet2Sno2,LatentHeatAir2Sno2
!    4,HeatSensAir2Sno2,HeatSensAir2SnoByEvap2,RA,MaxVapXAir2Sno,EVAPX2,VPQ(NY,NX),VPSno0
!    5,VLWatSnow0M(1,NY,NX),VLDrySnoWE0M(1,NY,NX),VLIceSnow0M(1,NY,NX)
!    6,HeatX2SnoLay(1,NY,NX),NetHeatAir2Snow,HWFLQ02,HeatNetFlx2Sno2,RadNet2Sno2,LatentHeatAir2Sno2
!    7,HeatSensAir2Sno2,HeatSensAir2SnoByEvap2,TKSnow1(1,NY,NX),TKQ(NY,NX)
!    8,PARE,RA,RZ,EvapSublimation2,EVAPW2,MaxVapXAir2Sno
!7759  FORMAT(A8,4I4,40E14.6)
!     ENDIF
!
  end subroutine SnowAtmosExchange  
!------------------------------------------------------------------------------------------
  subroutine PrepIterSnowLayer(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  integer :: L
! FLUX VARIABLES IN SNOWPACK
!
! PhaseChangeHeatL=latent heat from freeze-thaw
! SnowThawMassL,IceThawMassL=freeze-thaw between snow,ice and water
! FLW0S,FLW0W,FLW0I=snow,water,ice fluxes
! HFLW0W=convective heat flux from snow,water,ice fluxes
! WatFlowInSnowM=snowpack water flux
! VLDrySnoWE0M,VOLW0M,VOLI0M=snow,water,ice contents
! VHCPWMM=volumetric heat capacity
! TKSnow1=snow temperature
!
  D9765: DO L=1,JS
    PhaseChangeHeatL(L,NY,NX)=0.0_r8
    SnowThawMassL(L,NY,NX)=0.0_r8
    IceThawMassL(L,NY,NX)=0.0_r8
    SnoX2SnoLay(L,NY,NX)=0.0_r8
    WatX2SnoLay(L,NY,NX)=0.0_r8
    IceX2SnoLay(L,NY,NX)=0.0_r8
    HeatX2SnoLay(L,NY,NX)=0.0_r8
    WatFlowInSnowM(M,L,NY,NX)=0.0_r8
    VLDrySnoWE0M(L,NY,NX)=VLDrySnoWE0(L,NY,NX)
    VLWatSnow0M(L,NY,NX)=VLWatSnow0(L,NY,NX)
    VLIceSnow0M(L,NY,NX)=VLIceSnow0(L,NY,NX)
    VLHeatCapSnowM1(L,NY,NX)=VLSnowHeatCapM(M,L,NY,NX)    
    TKSnow1(L,NY,NX)=TKSnow0(L,NY,NX)
  ENDDO D9765
  end subroutine PrepIterSnowLayer  
!------------------------------------------------------------------------------------------

  subroutine UpdateSnowPack1(M,NY,NX)

  implicit NONE
  integer, intent(in) :: M,NY,NX
  
  integer :: L
  real(r8) :: ENGY0
  real(r8) :: dHPhaseChange,VLDrySnoWEtmp

!
!     VOLS0,VOLW0,VOLI0=snow,water,ice volumes in snowpack
!     TFLWS,TFLWW=net snow,water flux
!     SnowThawMassL,IceThawMassL=snow-water,ice-water freeze-thaw flux
!     DENSICE=ice density
!     VHCPWM=snowpack volumetric heat capacity
!     TK0=snowpack temperature
!     CumHeat2SnowLay=total snowpack conductive+convective heat flux
!     PhaseChangeHeatL=snowpack latent heat flux from freeze-thaw
!
!  if(M>=28 .and. NY==3)write(*,*)'UpdateSnowAtM0000 NY,NX ',NY,NX,M,VLDrySnoWE0(1,NY,NX),&
!    VLWatSnow0(1,NY,NX),VLIceSnow0(1,NY,NX)
  D9780: DO L=1,JS
    VLDrySnoWE0(L,NY,NX)=VLDrySnoWE0(L,NY,NX)+CumSno2SnowLay(L,NY,NX)-SnowThawMassL(L,NY,NX)
    VLWatSnow0(L,NY,NX)=VLWatSnow0(L,NY,NX)+CumWat2SnowLay(L,NY,NX)+SnowThawMassL(L,NY,NX)+IceThawMassL(L,NY,NX)
    VLIceSnow0(L,NY,NX)=VLIceSnow0(L,NY,NX)-IceThawMassL(L,NY,NX)/DENSICE
    ENGY0=VLSnowHeatCapM(M,L,NY,NX)*TKSnow0(L,NY,NX)
    VLSnowHeatCapM(M+1,L,NY,NX)=cps*VLDrySnoWE0(L,NY,NX)+cpw*VLWatSnow0(L,NY,NX)+cpi*VLIceSnow0(L,NY,NX)
    IF(VLSnowHeatCapM(M+1,L,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
      TKSnow0(L,NY,NX)=(ENGY0+CumHeat2SnowLay(L,NY,NX)+PhaseChangeHeatL(L,NY,NX))/VLSnowHeatCapM(M+1,L,NY,NX)
    ELSEIF(L.EQ.1)THEN
      TKSnow0(L,NY,NX)=TairK(NY,NX)
    ELSE
      TKSnow0(L,NY,NX)=TKSnow0(L-1,NY,NX)
    ENDIF
  ENDDO D9780
!  if(M>=28 .and. NY==3)write(*,*)'UpdateSnowAtM1111 NY,NX ',NY,NX,M,VLWatSnow0(1,NY,NX)
  if(VLWatSnow0(1,NY,NX)<0._r8)then
    L=1
!    write(*,*)VLDrySnoWE0(1,NY,NX),VLWatSnow0(1,NY,NX),VLIceSnow0(1,NY,NX)
!    write(*,*)CumWat2SnowLay(1,NY,NX),SnowThawMassL(1,NY,NX),IceThawMassL(1,NY,NX)
    !make up the negative water from dry snow
    VLDrySnoWEtmp=VLDrySnoWE0(L,NY,NX)+VLWatSnow0(L,NY,NX)
    !this is equivalent to thaw, which releases latent heat
    if(VLDrySnoWEtmp>0._r8)then
      VLDrySnoWE0(L,NY,NX)=VLDrySnoWEtmp
      ENGY0=VLSnowHeatCapM(M,L,NY,NX)*TKSnow0(L,NY,NX)
      VLSnowHeatCapM(M+1,L,NY,NX)=cps*VLDrySnoWE0(L,NY,NX)+cpw*VLWatSnow0(L,NY,NX)+cpi*VLIceSnow0(L,NY,NX)
      dHPhaseChange=-333.0_r8*VLWatSnow0(L,NY,NX)
      TKSnow0(L,NY,NX)=(ENGY0+dHPhaseChange)/VLSnowHeatCapM(M+1,L,NY,NX)
      PhaseChangeHeatL(L,NY,NX)=PhaseChangeHeatL(L,NY,NX)+dHPhaseChange
      VLWatSnow0(L,NY,NX)=0._r8      
    else
      call endrun(trim(mod_filename)//' at line',__LINE__)    
    endif

  endif

  end subroutine UpdateSnowPack1   

!------------------------------------------------------------------------------------------

  subroutine UpdateSnowAtM(M,NY,NX)
  implicit none    
  integer, intent(in) :: M,NY,NX
  
  integer :: L
  real(r8) :: tk1pres
  real(r8) :: ENGY0,ENGY1
  real(r8) :: FLWI,FLWW,FLWS
  real(r8) :: HFLWS
!     SNOWPACK WATER, ICE, SNOW AND TEMPERATURE

!  call UpdateSnowPack1(M,NY,NX)

  !      if(curday>=176)then
  !        write(*,*)'line',__LINE__,'tk1',TKSoi1(8,ny,nx),TKSoi1(9,ny,nx),M
  !      endif
  !
  !     SNOW RUNOFF
  !
  !  cumDrySnoFlxByRedistribut,cumWatFlxBySnowRedistribut,cumIceFlxBySnowRedistribut,
  !  cumHeatFlxBySnowRedistribut=net snow,water,ice, heat from snowpack runoff
  !
  

  VLDrySnoWE0(1,NY,NX)=VLDrySnoWE0(1,NY,NX)+cumDrySnoFlxByRedistribut(NY,NX)
  VLWatSnow0(1,NY,NX)=VLWatSnow0(1,NY,NX)+cumWatFlxBySnowRedistribut(NY,NX)
  VLIceSnow0(1,NY,NX)=VLIceSnow0(1,NY,NX)+cumIceFlxBySnowRedistribut(NY,NX)
  ENGY0=VLSnowHeatCapM(M+1,1,NY,NX)*TKSnow0(1,NY,NX)
  VLSnowHeatCapM(M+1,1,NY,NX)=cps*VLDrySnoWE0(1,NY,NX)+cpw*VLWatSnow0(1,NY,NX)+cpi*VLIceSnow0(1,NY,NX)
  IF(VLSnowHeatCapM(M+1,1,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
    TKSnow0(1,NY,NX)=(ENGY0+cumHeatFlxBySnowRedistribut(NY,NX))/VLSnowHeatCapM(M+1,1,NY,NX)
  ELSE
    TKSnow0(1,NY,NX)=TairK(NY,NX)
  ENDIF

  
!  if(M>=28 .and. NY==3)write(*,*)'UpdateSnowAtM NY,NX ',NY,NX,M,VLWatSnow0(1,NY,NX),cumWatFlxBySnowRedistribut(NY,NX)
  if(VLWatSnow0(1,NY,NX)<0._r8)then    
    write(*,*)'UpdateSnowAtM NY,NX ',NY,NX,M,VLWatSnow0(1,NY,NX),cumWatFlxBySnowRedistribut(NY,NX)
    call endrun(trim(mod_filename)//' at line',__LINE__)
  endif
    !
    !     IF SNOWPACK DISAPPEARS ALL MATERIAL,HEAT TRANSFERRED TO SOIL SURFACE
    !
    !     VLHeatCapacity,VLHeatCapacityA,VLHeatCapacityP=total soil,soil+micropore,macropore heat capacity
    !     TK1=soil surface temperature, why not to litter layer
    !
  IF(VLHeatCapSnow(1,NY,NX).LE.VLHeatCapSnowMin(NY,NX).AND.TairK(NY,NX).GT.TFice)THEN
    FLWS=VLDrySnoWE0(1,NY,NX)
    FLWW=VLWatSnow0(1,NY,NX)
    FLWI=VLIceSnow0(1,NY,NX)
    HFLWS=(cpw*FLWW+cps*FLWS+cpi*FLWI)*TKSnow0(1,NY,NX)
    VLDrySnoWE0(1,NY,NX)=VLDrySnoWE0(1,NY,NX)-FLWS
    VLWatSnow0(1,NY,NX)=VLWatSnow0(1,NY,NX)-FLWW
    VLIceSnow0(1,NY,NX)=VLIceSnow0(1,NY,NX)-FLWI
 
    !add to litter layer?
    VLWatMicP1(NUM(NY,NX),NY,NX)=VLWatMicP1(NUM(NY,NX),NY,NX)+FLWW
    VLiceMicP1(NUM(NY,NX),NY,NX)=VLiceMicP1(NUM(NY,NX),NY,NX)+FLWI+FLWS/DENSICE
    
    ENGY1=VLHeatCapacity(NUM(NY,NX),NY,NX)*TKSoi1(NUM(NY,NX),NY,NX)
    VLHeatCapacityA(NUM(NY,NX),NY,NX)=VHeatCapacitySoilM(NUM(NY,NX),NY,NX) &
      +cpw*VLWatMicP1(NUM(NY,NX),NY,NX)+cpi*VLiceMicP1(NUM(NY,NX),NY,NX)
    VLHeatCapacityB(NUM(NY,NX),NY,NX)=cpw*VLWatMacP1(NUM(NY,NX),NY,NX) &
      +cpi*VLiceMacP1(NUM(NY,NX),NY,NX)
    VLHeatCapacity(NUM(NY,NX),NY,NX)=VLHeatCapacityA(NUM(NY,NX),NY,NX)+VLHeatCapacityB(NUM(NY,NX),NY,NX)

    IF(VLHeatCapacity(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    ! topsoil layer is there
      tk1pres=TKSoi1(NUM(NY,NX),NY,NX)
      TKSoi1(NUM(NY,NX),NY,NX)=(ENGY1+HFLWS)/VLHeatCapacity(NUM(NY,NX),NY,NX)
!      if(abs(tk1pres/TKSoi1(NUM(NY,NX),NY,NX)-1._r8)>0.025_r8)then
!        TKSoi1(NUM(NY,NX),NY,NX)=TairK(NY,NX)
!      endif
    ELSE
      TKSoi1(NUM(NY,NX),NY,NX)=TairK(NY,NX)
    ENDIF
  ENDIF
!  if(M>=28 .and. NY==3)write(*,*)'UpdateSnowAtM222 NY,NX ',NY,NX,M,VLWatSnow0(1,NY,NX)  
  end subroutine UpdateSnowAtM

!------------------------------------------------------------------------------------------
  subroutine SnowRedistribution(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
!
! SNOW redistribution
! currently, it does consider wind effect
  implicit none
  integer, intent(in) :: M,NY,NX,NHE,NHW,NVS,NVN
  integer, intent(in) :: N1,N2

  integer :: N,NN,N4,N5,N4B,N5B
  real(r8) :: ALTS1,ALTS2
  real(r8) :: QSX,SnowDepthLateralGradient
  integer, parameter :: idirew=1
  integer, parameter :: idirns=2
!
!     SNOW REDISTRIBUTION FROM SNOWPACK
!
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     ALTS1,ALTS2=elevation of source,destination snowpack surfaces
!     SnowDepthLateralGradient,DIST=slope,distance between source,destination
!     QSX=transfer fraction
!     DrySnoFlxBySnowRedistribut,WatFlxBySnowRedistribut,IceFlxBySnowRedistribut=snow,water,ice transfer
!     HeatFlxBySnowRedistribut=convective heat transfer from snow,water,ice transfer
!     VOLS0,VOLW0,VOLI0=snow,water,ice volume
!     MinSnowDepth=minimum snowpack depth for full cover
!     QS,WatBySnowRedistribution,IceBySnowRedistribution=hourly-accumulated snow,water,ice transfer
!     HeatBySnowRedistribution=hourly-accumd convective heat from snow,water,ice transfer
!     DrySnoFlxBySnowRedistributM=snow transfer for solute flux calculation

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

      ALTS1=Altitude_grid(N2,N1)+SnowDepth(N2,N1)

      IF(NN.EQ.1)THEN
        !east or south        
        ALTS2=Altitude_grid(N5,N4)+SnowDepth(N5,N4)
        
        SnowDepthLateralGradient=(ALTS1-ALTS2)/DIST(N,NU(N5,N4),N5,N4)
        QSX=SnowDepthLateralGradient/AMAX1(1.0_r8,DIST(N,NU(N5,N4),N5,N4))*dts_HeatWatTP

        IF(SnowDepthLateralGradient.GT.0.0_r8.AND.SnowDepth(N2,N1).GT.MinSnowDepth)THEN
          DrySnoFlxBySnowRedistribut(N,N5,N4)=QSX*AZMAX1(VLDrySnoWE0(1,N2,N1))
          WatFlxBySnowRedistribut(N,N5,N4)=QSX*AZMAX1(VLWatSnow0(1,N2,N1))
          IceFlxBySnowRedistribut(N,N5,N4)=QSX*AZMAX1(VLIceSnow0(1,N2,N1))
          HeatFlxBySnowRedistribut(N,N5,N4)=TKSnow0(1,N2,N1)*(cps*DrySnoFlxBySnowRedistribut(N,N5,N4) &
            +cpw*WatFlxBySnowRedistribut(N,N5,N4)+cpi*IceFlxBySnowRedistribut(N,N5,N4))
        ELSEIF(SnowDepthLateralGradient.LT.0.0_r8.AND.SnowDepth(N5,N4).GT.MinSnowDepth)THEN
          DrySnoFlxBySnowRedistribut(N,N5,N4)=QSX*AZMAX1(VLDrySnoWE0(1,N5,N4))
          WatFlxBySnowRedistribut(N,N5,N4)=QSX*AZMAX1(VLWatSnow0(1,N5,N4))
          IceFlxBySnowRedistribut(N,N5,N4)=QSX*AZMAX1(VLIceSnow0(1,N5,N4))
          HeatFlxBySnowRedistribut(N,N5,N4)=TKSnow0(1,N5,N4)*(cps*DrySnoFlxBySnowRedistribut(N,N5,N4) &
            +cpw*WatFlxBySnowRedistribut(N,N5,N4)+cpi*IceFlxBySnowRedistribut(N,N5,N4))
        ELSE
          DrySnoFlxBySnowRedistribut(N,N5,N4)=0.0_r8
          WatFlxBySnowRedistribut(N,N5,N4)=0.0_r8
          IceFlxBySnowRedistribut(N,N5,N4)=0.0_r8
          HeatFlxBySnowRedistribut(N,N5,N4)=0.0_r8
        ENDIF
        DrysnoBySnowRedistribution(N,N5,N4)=DrysnoBySnowRedistribution(N,N5,N4)+DrySnoFlxBySnowRedistribut(N,N5,N4)
        WatBySnowRedistribution(N,N5,N4)=WatBySnowRedistribution(N,N5,N4)+WatFlxBySnowRedistribut(N,N5,N4)
        IceBySnowRedistribution(N,N5,N4)=IceBySnowRedistribution(N,N5,N4)+IceFlxBySnowRedistribut(N,N5,N4)
        HeatBySnowRedistribution(N,N5,N4)=HeatBySnowRedistribution(N,N5,N4)+HeatFlxBySnowRedistribut(N,N5,N4)
        DrySnoFlxBySnowRedistributM(M,N,N5,N4)=DrySnoFlxBySnowRedistribut(N,N5,N4)
      ENDIF

      !add west and south
      IF(NN.EQ.2)THEN

      ENDIF
    ENDDO
  ENDDO    
  end subroutine SnowRedistribution

!------------------------------------------------------------------------------------------
  subroutine SumSnowDriftByRunoff(M,N,N1,N2,N4,N5,N4B,N5B)
  use SoilWaterDataType, only : IFLBM
  implicit none 
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2
  integer, intent(in) :: N4,N5    !forward dest grid
  integer, intent(in) :: N4B,N5B  !backward dest grid
  integer :: NN

  D1202: DO NN=1,2
    cumWatFlx2LitRByRunoff(N2,N1)=cumWatFlx2LitRByRunoff(N2,N1)+WatFlx2LitRByRunoff(N,NN,N2,N1)
    cumHeatFlx2LitRByRunoff(N2,N1)=cumHeatFlx2LitRByRunoff(N2,N1)+HeatFlx2LitRByRunoff(N,NN,N2,N1)
    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      !there is runoff
      cumWatFlx2LitRByRunoff(N2,N1)=cumWatFlx2LitRByRunoff(N2,N1)-WatFlx2LitRByRunoff(N,NN,N5,N4)
      cumHeatFlx2LitRByRunoff(N2,N1)=cumHeatFlx2LitRByRunoff(N2,N1)-HeatFlx2LitRByRunoff(N,NN,N5,N4)
    ENDIF

    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      cumWatFlx2LitRByRunoff(N2,N1)=cumWatFlx2LitRByRunoff(N2,N1)-WatFlx2LitRByRunoff(N,NN,N5B,N4B)
      cumHeatFlx2LitRByRunoff(N2,N1)=cumHeatFlx2LitRByRunoff(N2,N1)-HeatFlx2LitRByRunoff(N,NN,N5B,N4B)
    ENDIF

    IF(M.EQ.NPH)THEN
      IFLBH(N,NN,N5,N4)=IFLBM(M,N,NN,N5,N4)
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        IFLBH(N,NN,N5B,N4B)=IFLBM(M,N,NN,N5B,N4B)
      ENDIF
    ENDIF
  ENDDO D1202

  cumDrySnoFlxByRedistribut(N2,N1)=cumDrySnoFlxByRedistribut(N2,N1)+DrySnoFlxBySnowRedistribut(N,N2,N1) &
    -DrySnoFlxBySnowRedistribut(N,N5,N4)
  cumWatFlxBySnowRedistribut(N2,N1)=cumWatFlxBySnowRedistribut(N2,N1)+WatFlxBySnowRedistribut(N,N2,N1) &
    -WatFlxBySnowRedistribut(N,N5,N4)
  cumIceFlxBySnowRedistribut(N2,N1)=cumIceFlxBySnowRedistribut(N2,N1)+IceFlxBySnowRedistribut(N,N2,N1) &
    -IceFlxBySnowRedistribut(N,N5,N4)
  cumHeatFlxBySnowRedistribut(N2,N1)=cumHeatFlxBySnowRedistribut(N2,N1)+HeatFlxBySnowRedistribut(N,N2,N1) &
    -HeatFlxBySnowRedistribut(N,N5,N4)
  END subroutine SumSnowDriftByRunoff
!------------------------------------------------------------------------------------------
  subroutine ZeroSnowFlux(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  cumWatFlx2LitRByRunoff(NY,NX)=0.0_r8
  cumHeatFlx2LitRByRunoff(NY,NX)=0.0_r8
  cumDrySnoFlxByRedistribut(NY,NX)=0.0_r8
  cumWatFlxBySnowRedistribut(NY,NX)=0.0_r8
  cumIceFlxBySnowRedistribut(NY,NX)=0.0_r8
  cumHeatFlxBySnowRedistribut(NY,NX)=0.0_r8

  end subroutine ZeroSnowFlux  
!------------------------------------------------------------------------------------------
  FUNCTION SnowBNDResistance(NY,NX)result(RAS)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: RAS
  real(r8) :: RASL,RASX
  real(r8) :: FracAsAirSno
  integer :: L

!     RAS,RASL=blrs of snowpack,snowpack layer
!     VOLS,VLSnoDWI1=volume of snowpack,snowpack layer
!     DLYRS=snowpack later depth
!     WGSGW=vapor diffusivity in snowpack
!     FracAsAirSno=snowpack air-filled porosity
!     THETPI=air content of ice
!     VOLS0,VOLI0,VOLW0,VLSnoDWI1=snow,ice,water,total volumes of snowpack

  RAS=0.0_r8
  IF(VcumSnoDWI(NY,NX).GT.ZEROS2(NY,NX))THEN
    D9775: DO L=1,JS
      IF(VLSnoDWI1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        RASX=SnowLayerThick(L,NY,NX)/H2OVapDifscSno(L,NY,NX)
        FracAsAirSno=AMAX1(THETPI,1.0_r8-(VLDrySnoWE0(L,NY,NX)+VLIceSnow0(L,NY,NX) &
            +VLWatSnow0(L,NY,NX))/VLSnoDWI1(L,NY,NX))
        RASL=RASX/AMAX1(ZERO,FracAsAirSno)**2.0_r8
        RAS=RAS+RASL
      ENDIF
    ENDDO D9775
  ENDIF
  END FUNCTION SnowBNDResistance      

!------------------------------------------------------------------------------------------

  subroutine SnowSurfLitRIterate(L,M,NY,NX,AvgThermCondctSnoLitR,AvgThermCondctSoilLitR,&
    AvgVaporCondctSoilLitR,AvgVaporCondctSnowLitR,PSISV1,VLairSno1,TK0X,TKXR,TK1X,CumVapFlxSno2Litr,&
    CumVapFlxLitr2Soi,cumHeatConvFlxLitr2Soi1,CumHeatConvFlxSno2Litr,CumHeatCndFlxSno2Litr,&
    cumHeatCndFlxLitr2Soi)
  implicit none
  integer, intent(in) :: L,M,NY,NX
  real(r8), intent(in) :: AvgThermCondctSnoLitR,AvgThermCondctSoilLitR,AvgVaporCondctSnowLitR
  real(r8), intent(in) :: AvgVaporCondctSoilLitR,PSISV1,VLairSno1
  real(r8), intent(inout) :: TK0X,TKXR,TK1X
  real(r8), intent(inout) :: CumVapFlxSno2Litr,CumVapFlxLitr2Soi
  real(r8), intent(inout) :: cumHeatConvFlxLitr2Soi1,CumHeatConvFlxSno2Litr
  real(r8), intent(inout) :: CumHeatCndFlxSno2Litr,cumHeatCndFlxLitr2Soi
  integer :: NN
  real(r8) :: VP1,VPR,VapSnow0,VPY
  real(r8) :: H2OVapFlx,VapConvFlxLitr2Soi
  real(r8) :: HeatConvFlxLitr2Soi
  real(r8) :: VapFlxSno2Litr,HeatConvFlxSno2Litr,HeatCndFlxSno2Litr
  real(r8) :: H2OVapFlxMax,HeatCnduct,HeatCndFlxLitr2Soi
  real(r8) :: HeatCnductMax,TKY

! begin_execution
  D4000: DO NN=1,NPR
    !
    ! VAPOR FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! VapSnow0,VPR,VPY=snowpack,litter, equilibrium vapor concentration
    ! TK0X,TKXR=snowpack,litter temperature
    ! PSISM1=litter matric water potential
    ! H2OVapFlx,H2OVapFlxMax=vapor-unconstrained,vapor-constrained vapor flux
    ! AREA=area of grid cell
    ! FSNW,CVRD=snow,litter cover fraction
    ! dts_litrvapht=time step for flux calculation
    ! FLVRSX=snow-litter vapor flux
    ! HFLVRSX=convective heat flux from snow-litter vapor flux
    !
    !residue vapor pressure
    VPR=vapsat(TKXR)*EXP(18.0_r8*PSISM1(0,NY,NX)/(RGAS*TKXR))
    if(abs(VPR)>1.e20_r8)then
      write(*,*)'TKXR=',TKXR,TKSoi1(0,NY,NX),TKSoi1(NUM(NY,NX),NY,NX),NN
      write(*,*)'PSISM1(0,NY,NX)=',PSISM1(0,NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif

    !AvgVaporCondctSnowLitR=snowpack-litter conductance
    !there is litter layer with sufficient moisture
    IF(VLairSno1.GT.ZEROS2(NY,NX).AND.THETPM(M,0,NY,NX).GT.THETX)THEN      
      VapSnow0=vapsat(TK0X) !snow vapor pressure, saturated
      !snow <-> residue vapor flux
      H2OVapFlx=AvgVaporCondctSnowLitR*(VapSnow0-VPR)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX) &
        *FracSurfByLitR(NY,NX)*dts_litrvapht 
      !volume weighted vapor pressure
      VPY=(VapSnow0*VLairSno1+VPR*VLsoiAirPM(M,0,NY,NX))/(VLairSno1+VLsoiAirPM(M,0,NY,NX))             
      H2OVapFlxMax=(VapSnow0-VPY)*VLairSno1*dt_watvap

      IF(H2OVapFlx.GE.0.0_r8)THEN
        !vapor flux from snow to litter
        VapFlxSno2Litr=AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax))       !water flux
        HeatConvFlxSno2Litr=(cpw*TK0X+EvapLHTC)*VapFlxSno2Litr         !enthalpy flux associated with water flux
      ELSE
        !vapor flux from litter to snow
        VapFlxSno2Litr=AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax))
        HeatConvFlxSno2Litr=(cpw*TKXR+EvapLHTC)*VapFlxSno2Litr
      ENDIF
    ELSE
      VapFlxSno2Litr=0.0_r8
      HeatConvFlxSno2Litr=0.0_r8
    ENDIF
    !
    ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! TKY=snow-litter equilibrium temperature
    ! HeatCnduct,HeatCnductMax=snow-litter heat flux unltd,ltd by heat
    ! HeatCndFlxSno2Litr=snow-litter heat flux
    ! VLSnowHeatCapMM= volumetric heat capacity in snow layer
    TKY=(TK0X*VLHeatCapSnowM1(L,NY,NX)+TKXR*VLHeatCapacity(0,NY,NX))/(VLHeatCapSnowM1(L,NY,NX)+VLHeatCapacity(0,NY,NX))
    HeatCnductMax=(TK0X-TKY)*VLHeatCapSnowM1(L,NY,NX)*dt_watvap
    HeatCnduct=AvgThermCondctSnoLitR*(TK0X-TKXR)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX)*FracSurfByLitR(NY,NX)*dts_litrvapht
    IF(HeatCnduct.GE.0.0_r8)THEN
      HeatCndFlxSno2Litr=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
    ELSE
      HeatCndFlxSno2Litr=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
    ENDIF
!
!     VAPOR FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
!
!     THETPM,VsoiPM=air-filled porosity,volume
!     VP1,VPY=soil,litter-soil equilibrium vapor concentration
!     TK1X=soil temperature
!     PSISV1=soil matric+osmotic water potentials
!     H2OVapFlx,H2OVapFlxMax=vapor-unconstrained,vapor-constrained vapor flux
!     VapConvFlxLitr2Soi=litter-soil vapor flux
!     HeatConvFlxLitr2Soi=convective heat of litter-soil vapor flux
!     TKXR,TK1X=interim calculation of litter,soil temperatures
!
    !both litter layer and topsoil are not-saturated
    IF(VLsoiAirPM(M,0,NY,NX).GT.ZEROS(NY,NX).AND.VLsoiAirPM(M,NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      VP1=vapsat(TK1X)*EXP(18.0_r8*PSISV1/(RGAS*TK1X))
      H2OVapFlx=AvgVaporCondctSoilLitR*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX)*FracSurfByLitR(NY,NX)*dts_litrvapht

      if(abs(H2OVapFlx)>1.e20_r8)then
        write(*,*)'AvgVaporCondctSoilLitR=',AvgVaporCondctSoilLitR,VPR,VP1
        write(*,*)'FracSurfAsSnow(NY,NX)*FracSurfByLitR(NY,NX)=',FracSurfAsSnow(NY,NX),FracSurfByLitR(NY,NX)
        write(*,*)'at line',__LINE__
        call endrun(trim(mod_filename)//'at line',__LINE__)
      endif
      VPY=(VPR*VLsoiAirPM(M,0,NY,NX)+VP1*VLsoiAirPM(M,NUM(NY,NX),NY,NX)) &
        /(VLsoiAirPM(M,0,NY,NX)+VLsoiAirPM(M,NUM(NY,NX),NY,NX))
      H2OVapFlxMax=(VPR-VPY)*VLsoiAirPM(M,0,NY,NX)*dt_watvap
      IF(H2OVapFlx.GE.0.0_r8)THEN
        VapConvFlxLitr2Soi=AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax,VLWatSnow0M(L,NY,NX)*XNPB))
        if(abs(VapConvFlxLitr2Soi)>1.0e20_r8)then
          write(*,*)'H2OVapFlx,H2OVapFlxMax,VLWatSnow0M(L,NY,NX)*XNPB=',H2OVapFlx,H2OVapFlxMax,VLWatSnow0M(L,NY,NX)*XNPB
          write(*,*)'at line',__LINE__
          call endrun(trim(mod_filename)//'at line',__LINE__)
        endif
        HeatConvFlxLitr2Soi=(cpw*TKXR+EvapLHTC)*VapConvFlxLitr2Soi  !enthalpy flux
      ELSE
        VapConvFlxLitr2Soi=AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax))
        if(abs(VapConvFlxLitr2Soi)>1.0e20_r8)then
          write(*,*)'H2OVapFlx,H2OVapFlxMax=',H2OVapFlx,H2OVapFlxMax
          write(*,*)'at line',__LINE__
          call endrun(trim(mod_filename)//'at line',__LINE__)
        endif
        HeatConvFlxLitr2Soi=(cpw*TK1X+EvapLHTC)*VapConvFlxLitr2Soi
      ENDIF
    ELSE
      VapConvFlxLitr2Soi=0.0_r8
      HeatConvFlxLitr2Soi=0.0_r8
    ENDIF
    !update litter layer temperature
    TKXR=TKXR-HeatConvFlxLitr2Soi/VLHeatCapacity(0,NY,NX)   
    if(TKXR<0._r8)then
      write(*,*)'unphysical litter temeprature',NN
      write(*,*)'TKXR',TKXR,HeatConvFlxLitr2Soi,VLHeatCapacity(0,NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    !update top soil layer temperature
    TK1X=TK1X+HeatConvFlxLitr2Soi/VLHeatCapacity(NUM(NY,NX),NY,NX)
!
!     HEAT FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
!
!     TKY=litter-soil equilibrium temperature
!     HeatCnduct,HeatCnductMax=litter-soil heat flux unltd,ltg by heat
!     HeatCndFlxLitr2Soi=litter-soil heat flux
!
    TKY=(TKXR*VLHeatCapacity(0,NY,NX)+TK1X*VLHeatCapacity(NUM(NY,NX),NY,NX)) &
      /(VLHeatCapacity(0,NY,NX)+VLHeatCapacity(NUM(NY,NX),NY,NX))
    HeatCnductMax=(TKXR-TKY)*VLHeatCapacity(0,NY,NX)*dt_watvap
    HeatCnduct=AvgThermCondctSoilLitR*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX) &
      *FracSurfByLitR(NY,NX)*dts_litrvapht
    
    IF(HeatCnduct.GE.0.0_r8)THEN
      HeatCndFlxLitr2Soi=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
    ELSE
      HeatCndFlxLitr2Soi=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
    ENDIF
!
!     ACCUMULATE SNOW-LITTER, LITTER-SOIL HEAT FLUXES
!     WITHIN LONGER TIME STEP FOR SNOWPACK FLUX CALCULATIONS
!
    CumVapFlxSno2Litr=CumVapFlxSno2Litr+VapFlxSno2Litr
    CumHeatConvFlxSno2Litr=CumHeatConvFlxSno2Litr+HeatConvFlxSno2Litr
    CumHeatCndFlxSno2Litr=CumHeatCndFlxSno2Litr+HeatCndFlxSno2Litr
    CumVapFlxLitr2Soi=CumVapFlxLitr2Soi+VapConvFlxLitr2Soi
    if(abs(CumVapFlxLitr2Soi)>1.0e20_r8)then
      write(*,*)'VapConvFlxLitr2Soi=',VapConvFlxLitr2Soi
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    cumHeatConvFlxLitr2Soi1=cumHeatConvFlxLitr2Soi1+HeatConvFlxLitr2Soi
    cumHeatCndFlxLitr2Soi=cumHeatCndFlxLitr2Soi+HeatCndFlxLitr2Soi
    TK0X=TK0X-HeatConvFlxSno2Litr/VLHeatCapSnowM1(L,NY,NX)
    TKXR=TKXR+(HeatConvFlxSno2Litr-HeatCndFlxLitr2Soi)/VLHeatCapacity(0,NY,NX)
    TK1X=TK1X+HeatCndFlxLitr2Soi/VLHeatCapacity(NUM(NY,NX),NY,NX)
    if(TKXR<0._r8)then
      write(*,*)'unphysical litter temeprature'
      write(*,*)'TKXR',TKXR,HeatConvFlxSno2Litr,HeatCndFlxLitr2Soi,VLHeatCapacity(0,NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif  

  ENDDO D4000
  
  end subroutine SnowSurfLitRIterate

!------------------------------------------------------------------------------------------
  subroutine CopySnowStates(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
! VOLS0,VOLSSL=snowpack snow content (water equivalent)
! VOLI0,VOLISSL=snowpack ice content
! VOLW0,VOLWSL=snowpack water content
! VLSnoDWI1,VLSnoDWI=snowpack volume
! SnowLayerThick0,DLYRS=snowpack depth
! VLSnowHeatCapM,VHCPW=snowpack heat capacity
! TK0,TKW=snowpack temperature
!
  D60: DO L=1,JS
    VLDrySnoWE0(L,NY,NX)=VLDrySnoWE(L,NY,NX)
    VLIceSnow0(L,NY,NX)=VLIceSnow(L,NY,NX)
    VLWatSnow0(L,NY,NX)=VLWatSnow(L,NY,NX)
    VLSnoDWI1(L,NY,NX)=VLSnoDWI(L,NY,NX)
    SnowLayerThick0(L,NY,NX)=SnowLayerThick(L,NY,NX)
    VLSnowHeatCapM(1,L,NY,NX)=VLHeatCapSnow(L,NY,NX)
    TKSnow0(L,NY,NX)=TKSnow(L,NY,NX)
  ENDDO D60  
  end subroutine CopySnowStates
!------------------------------------------------------------------------------------------

  subroutine UpdateSoilWaterPotential(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: THETW1,VPY,FCX,WPX
  real(r8) :: LOGFCX,LOGWPX,PSDX,FCDX

  THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX),AMIN1(POROS(NUM(NY,NX),NY,NX) &
    ,safe_adb(VLWatMicP1(NUM(NY,NX),NY,NX),VLSoilMicP(NUM(NY,NX),NY,NX))))
  IF(SoilMicPMassLayer(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    IF(THETW1.LT.FieldCapacity(NUM(NY,NX),NY,NX))THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX) &
        +((LOGFldCapacity(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
        /FCD(NUM(NY,NX),NY,NX)*LOGPSIMND(NY,NX))))
    ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(LOGPSIAtSat(NY,NX) &
        +(((LOGPOROS(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
        /PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*LOGPSIMXD(NY,NX)))
    ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
    ENDIF
  ELSEIF(VLSoilPoreMicP(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    FCX=FCI*FracSoiPAsIce(NUM(NY,NX),NY,NX)
    WPX=WPI*FracSoiPAsIce(NUM(NY,NX),NY,NX)
    LOGFCX=LOG(FCX)
    LOGWPX=LOG(WPX)
    PSDX=LOGPOROS(NUM(NY,NX),NY,NX)-LOGFCX
    FCDX=LOGFCX-LOGWPX
    IF(FracSoiPAsWat(NUM(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX) &
        +((LOGFCX-LOG(FracSoiPAsWat(NUM(NY,NX),NY,NX)))/FCDX*LOGPSIMND(NY,NX))))
    ELSEIF(FracSoiPAsWat(NUM(NY,NX),NY,NX).LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(LOGPSIAtSat(NY,NX) &
        +(((LOGPOROS(NUM(NY,NX),NY,NX)-LOG(FracSoiPAsWat(NUM(NY,NX),NY,NX)))/PSDX)*LOGPSIMXD(NY,NX)))
    ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
    ENDIF
  ELSE
    THETW1=POROS(NUM(NY,NX),NY,NX)
    PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
  ENDIF
  end subroutine UpdateSoilWaterPotential  
!------------------------------------------------------------------------------------------
  subroutine SnowSurLitterExch(M,L,NY,NX,VapCond1,VapCond2,TCND1W,VLairSno1,PSISV1,TCNDS,&
    CumVapFlxSno2Litr,CumHeatConvFlxSno2Litr,CumHeatCndFlxSno2Litr,CumVapFlxLitr2Soi,&
    cumHeatConvFlxLitr2Soi1,cumHeatCndFlxLitr2Soi )
  implicit none
  integer  , intent(in) :: M,L,NY,NX
  real(r8) , intent(in) :: VapCond1,VapCond2,TCND1W,VLairSno1,PSISV1,TCNDS
  real(r8) , intent(inout) :: CumVapFlxSno2Litr,CumHeatConvFlxSno2Litr
  real(r8) , intent(inout) :: CumHeatCndFlxSno2Litr,CumVapFlxLitr2Soi
  real(r8) , intent(inout) :: cumHeatConvFlxLitr2Soi1,cumHeatCndFlxLitr2Soi
  real(r8) :: TK0X,TKXR,TK1X
  real(r8) :: CNVR,AvgVaporCondctSnowLitR,AvgVaporCondctSoilLitR
  real(r8) :: THETRR,TCNDR,AvgThermCondctSnoLitR,AvgThermCondctSoilLitR
  real(r8) :: THETWR

  ! THETWR,THETW1=litter, soil water concentration
  ! VWatLitRHoldCapcity=litter water retention capacity
  ! PSISM1(0,PSISM1(NUM=litter,soil water potentials

  IF(VLitR(NY,NX).GT.ZEROS(NY,NX).AND.VLWatMicP1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    THETWR=AMIN1(VWatLitRHoldCapcity(NY,NX),VLWatMicP1(0,NY,NX))/VLitR(NY,NX)
    IF(THETWR.LT.FieldCapacity(0,NY,NX))THEN
      PSISM1(0,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX)+((LOGFldCapacity(0,NY,NX)-LOG(THETWR))/FCD(0,NY,NX) &
        *LOGPSIMND(NY,NX))))
    ELSEIF(THETWR.LT.POROS0(NY,NX))THEN
      PSISM1(0,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS(0,NY,NX)-LOG(THETWR))/PSD(0,NY,NX))**SRP(0,NY,NX) &
        *LOGPSIMXD(NY,NX)))
    ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
    ENDIF
  ELSE
    THETWR=POROS0(NY,NX)
    PSISM1(0,NY,NX)=PSISE(0,NY,NX)
  ENDIF


  TK0X=TKSnow1(L,NY,NX)
  TKXR=TKSoi1(0,NY,NX)
  TK1X=TKSoi1(NUM(NY,NX),NY,NX)
  CNVR=VaporDiffusivityLitR(NY,NX)*THETPM(M,0,NY,NX)*POROQ*THETPM(M,0,NY,NX)/POROS(0,NY,NX)
  if(TKXR<0._r8)then
    write(*,*)'SnowSurLitterExch negative M, L, TKR',M,L,TKXR
  endif
  IF(FracSurfByLitR(NY,NX).GT.ZERO)THEN
    IF(VapCond1.GT.ZERO.AND.CNVR.GT.ZERO)THEN
      AvgVaporCondctSnowLitR=2.0_r8*CNVR*VapCond1/(VapCond1*DLYRR(NY,NX)+CNVR*SnowLayerThick0(L,NY,NX))
    ELSE
      AvgVaporCondctSnowLitR=2.0_r8*VapCond1/(DLYRR(NY,NX)+SnowLayerThick0(L,NY,NX))
    ENDIF

    IF(CNVR.GT.ZERO.AND.VapCond2.GT.ZERO)THEN
      AvgVaporCondctSoilLitR=2.0_r8*CNVR*VapCond2/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+VapCond2*DLYRR(NY,NX))
    ELSE
      AvgVaporCondctSoilLitR=2.0_r8*VapCond2/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))
    ENDIF
    THETRR=AZMAX1(1.0_r8-FracSoiPAsAir(0,NY,NX)-FracSoiPAsWat(0,NY,NX)-FracSoiPAsIce(0,NY,NX))
    TCNDR=(0.779_r8*THETRR*9.050E-04_r8+0.622_r8*FracSoiPAsWat(0,NY,NX) &
      *2.067E-03_r8+0.380_r8*FracSoiPAsIce(0,NY,NX)*7.844E-03_r8+FracSoiPAsAir(0,NY,NX) &
      *9.050E-05_r8)/(0.779_r8*THETRR+0.622_r8*FracSoiPAsWat(0,NY,NX) &
      +0.380_r8*FracSoiPAsIce(0,NY,NX)+FracSoiPAsAir(0,NY,NX))

    IF(TCND1W.GT.ZERO.AND.TCNDR.GT.ZERO)THEN
      AvgThermCondctSnoLitR=2.0_r8*TCND1W*TCNDR/(TCND1W*DLYRR(NY,NX)+TCNDR*SnowLayerThick0(L,NY,NX))
    ELSE
      AvgThermCondctSnoLitR=0.0_r8
    ENDIF

    IF(TCNDR.GT.ZERO.AND.TCNDS.GT.ZERO)THEN
      AvgThermCondctSoilLitR=2.0_r8*TCNDR*TCNDS/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX)+TCNDS*DLYRR(NY,NX))
    ELSE
      AvgThermCondctSoilLitR=0.0_r8
    ENDIF
  ELSE
    AvgVaporCondctSnowLitR=0.0_r8
    AvgVaporCondctSoilLitR=0.0_r8
    AvgThermCondctSnoLitR=0.0_r8
    AvgThermCondctSoilLitR=0.0_r8
  ENDIF
  !
  ! SHORTER TIME STEP FOR SURFACE RESIDUE FLUX CALCULATIONS
  call SnowSurfLitRIterate(L,M,NY,NX,AvgThermCondctSnoLitR,AvgThermCondctSoilLitR,&
    AvgVaporCondctSoilLitR,AvgVaporCondctSnowLitR,PSISV1,VLairSno1,TK0X,TKXR,TK1X,CumVapFlxSno2Litr,&
    CumVapFlxLitr2Soi,cumHeatConvFlxLitr2Soi1,CumHeatConvFlxSno2Litr,CumHeatCndFlxSno2Litr,&
    cumHeatCndFlxLitr2Soi)

  end subroutine SnowSurLitterExch

!------------------------------------------------------------------------------------------

  subroutine SnowTopSoilExch(M,L,NY,NX,VapCond1,VapSnoSrc,VLairSno1,TCND1W,&
    VapFlxSno2Soi1,HeatConvFlxSno2Soi1,HeatCndFlxSno2Soi,VapCond2,PSISV1,TCNDS)
  implicit none
  integer , intent(in) :: M,L,NY,NX
  real(r8), intent(in) :: VapCond1,VapSnoSrc,VLairSno1,TCND1W
  real(r8), intent(out):: VapFlxSno2Soi1,HeatConvFlxSno2Soi1,HeatCndFlxSno2Soi,VapCond2,PSISV1,TCNDS
  real(r8) :: AvgThermCondctSoilLitR
  real(r8) :: WTHET2,VPY,VapSoiDest,AvgVaporCondctSoilLitR,H2OVapFlx,H2OVapFlxMax
  real(r8) :: HeatCnduct,HeatCnductMax,TKWX1,TKY

  call UpdateSoilWaterPotential(NY,NX)
  !micropore pressure, excluding gravitational pressure
  PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISoilOsmotic(NUM(NY,NX),NY,NX)
!
  ! VAPOR FLUX BETWEEN SNOWPACK AND SOIL SURFACE
  !
  ! VLairSno1,THETPM=air volume,concentration
  ! VapCond1,VapCond2=vapor conductances of source, destination layers
  ! VapSnoSrc,VapSoiDest=vapor concentrations of source, destination layers
  ! POROS,POROQ=porosity, tortuosity
  ! WGSGL=vapor diffusivity
  ! TKSnow1,TK1=snow,soil surface temperature
  ! PSISV1=soil matric+osmotic potential
  ! AvgVaporCondctSoilLitR=snow-soil vapor conductance
  ! DLYR=soil surface layer depth
  ! H2OVapFlx,H2OVapFlxMax=vapor flux unlimited,limited by vapor
  ! VPY=equilibrium vapor concentration
  ! dts_wat=time step for flux calculations
  ! VapFlxSno2Soi1,HeatConvFlxSno2Soi1=vapor flux and its convective heat flux
  !
  IF(VLairSno1.GT.ZEROS2(NY,NX).AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
    VapCond2=WGSGL(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ &
      *THETPM(M,NUM(NY,NX),NY,NX)/POROS(NUM(NY,NX),NY,NX)
    VapSoiDest=vapsat(TKSoi1(NUM(NY,NX),NY,NX))*EXP(18.0_r8*PSISV1/(RGAS*TKSoi1(NUM(NY,NX),NY,NX)))
    AvgVaporCondctSoilLitR=2.0_r8*VapCond1*VapCond2/(VapCond1*DLYR(3,NUM(NY,NX),NY,NX)+VapCond2*SnowLayerThick0(L,NY,NX))
    H2OVapFlx=AvgVaporCondctSoilLitR*(VapSnoSrc-VapSoiDest)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfAsSnow(NY,NX)&
      *FracSurfAsBareSoi(NY,NX)*dts_snohttp
    VPY=(VapSnoSrc*VLairSno1+VapSoiDest*VLsoiAirPM(M,NUM(NY,NX),NY,NX))/(VLairSno1+VLsoiAirPM(M,NUM(NY,NX),NY,NX))
    H2OVapFlxMax=(VapSnoSrc-VPY)*VLairSno1*dts_sno

    IF(H2OVapFlx.GE.0.0_r8)THEN
      !water flux goes into soil
      VapFlxSno2Soi1=AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax))
      !heat flux
      HeatConvFlxSno2Soi1=(cpw*TKSnow1(L,NY,NX)+EvapLHTC)*VapFlxSno2Soi1
    ELSE
      !water flux out of soil
      VapFlxSno2Soi1=AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax))
      HeatConvFlxSno2Soi1=(cpw*TKSoi1(NUM(NY,NX),NY,NX)+EvapLHTC)*VapFlxSno2Soi1
    ENDIF
  ELSE
    VapCond2=0.0_r8
    VapFlxSno2Soi1=0.0_r8
    HeatConvFlxSno2Soi1=0.0_r8
  ENDIF

  !
  ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE SOIL
  !
  ! WTHET2=multiplier for air concentration in thermal conductivity
  ! TCND1W,TCNDS=thermal conductivity of snowpack, soil surface
  ! STC,DTC=mineral component of thermal conductivity
  ! FracSoiPAsWat,FracSoiPAsIce,FracSoiPAsAir=soil surface water,ice,air concentrations
  ! BARE=soil surface fraction
  ! AvgThermCondctSoilLitR=snowpack-soil thermal conductance
  ! TKWX1=interim snowpack temperature
  ! TKY=equilibrium temperature
  ! HeatCnductMax,HeatCnduct=heat-constrained,heat-unconstrained heat fluxes
  ! dts_snohttp=time step for snowpack flux calculations
  ! HeatCndFlxSno2Soi=snowpack-soil heat flux
  ! FracSoilAsAirt=air-filled pores (including macropores and micropores)
  WTHET2=1.467_r8-0.467_r8*FracSoilAsAirt(NUM(NY,NX),NY,NX)
  TCNDS=(STC(NUM(NY,NX),NY,NX)+FracSoiPAsWat(NUM(NY,NX),NY,NX) &
    *2.067E-03_r8+0.611_r8*FracSoiPAsIce(NUM(NY,NX),NY,NX)*7.844E-03_r8 &
    +WTHET2*FracSoiPAsAir(NUM(NY,NX),NY,NX)*9.050E-05_r8) &
    /(DTC(NUM(NY,NX),NY,NX)+FracSoiPAsWat(NUM(NY,NX),NY,NX) &
    +0.611_r8*FracSoiPAsIce(NUM(NY,NX),NY,NX) &
    +WTHET2*FracSoiPAsAir(NUM(NY,NX),NY,NX))

  !the mean thermal conductivity
  IF(FracSurfAsBareSoi(NY,NX).GT.ZERO)THEN
    AvgThermCondctSoilLitR=2.0_r8*TCND1W*TCNDS/(TCND1W*DLYR(3,NUM(NY,NX),NY,NX)+TCNDS*SnowLayerThick0(L,NY,NX))
  ELSE
    AvgThermCondctSoilLitR=0.0_r8
  ENDIF

  TKWX1=TKSoi1(NUM(NY,NX),NY,NX)+HeatConvFlxSno2Soi1/VLHeatCapacity(NUM(NY,NX),NY,NX)
  TKY=(TKSnow1(L,NY,NX)*VLHeatCapSnowM1(L,NY,NX)+TKWX1*VLHeatCapacity(NUM(NY,NX),NY,NX)) &
    /(VLHeatCapSnowM1(L,NY,NX)+VLHeatCapacity(NUM(NY,NX),NY,NX))
  HeatCnductMax=(TKSnow1(L,NY,NX)-TKY)*VLHeatCapSnowM1(L,NY,NX)*dts_sno
  HeatCnduct=AvgThermCondctSoilLitR*(TKSnow1(L,NY,NX)-TKWX1)*AREA(3,NUM(NY,NX),NY,NX) &
    *FracSurfAsSnow(NY,NX)*FracSurfAsBareSoi(NY,NX)*dts_snohttp
  IF(HeatCnduct.GE.0.0_r8)THEN
    HeatCndFlxSno2Soi=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
  ELSE
    HeatCndFlxSno2Soi=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
  ENDIF
  end subroutine SnowTopSoilExch
end module SnowPhysMod
