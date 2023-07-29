module SnowPhysMod
!
! Description
! the snow model
! required input
!

! codes for snow physics
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : safe_adb,AZMAX1  
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
  public :: SnowRedistrub
  public :: SnowBNDResistance
  public :: ZeroSnowFlux  
  public :: PrepIterSnowLayer
  public :: InitSnowAccums
  public :: CopySnowStates
  public :: SolveSnowpack
  public :: SumSnowRoffDrift
  public :: UpdateSnowAtM
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
! DENSI=ice density, 0.92 Mg m-3
! DENS0=snow density (Mg m-3)
! VOLSS,VOLWS,VOLIS,VOLS=snow,water,ice,total snowpack volume(m3)

! cumSnowDepthI=depth to bottom of snowpack layers (m), i.e. from current layere to surface
! DLYRS=snowpack layer thickness (m)
! VOLSSL,VOLWSL,VOLISL,VOLSL=snow,water,ice,total layer volume(m3), water equivalent snow 
! DENSS=layer density (Mg m-3)
! TKW,TCW=later temperature K,oC
! VHCPW=layer volumetric heat capacity (MJ m-3 K-1)
! SnowDepth=total snow height in the column

  cumSnowDepth(0,NY,NX)=0.0_r8
  DENS0(NY,NX)=0.10_r8
  VOLSS(NY,NX)=SnowDepth(NY,NX)*DENS0(NY,NX)*DH(NY,NX)*DV(NY,NX)
  VOLWS(NY,NX)=0.0_r8
  VOLIS(NY,NX)=0.0_r8
  VOLS(NY,NX)=VOLSS(NY,NX)/DENS0(NY,NX)+VOLWS(NY,NX)+VOLIS(NY,NX)

!  VOLSWI=0.0_r8

  !build the snow profile, topdown
  D9580: DO L=1,JS
    IF(L.EQ.1)THEN
      !bottom snow layer
      DLYRSI=cumSnowDepthI(L)
      DLYRS(L,NY,NX)=AMIN1(DLYRSI,SnowDepth(NY,NX))
    ELSE
      DLYRSI=cumSnowDepthI(L)-cumSnowDepthI(L-1)
      DLYRS(L,NY,NX)=AMIN1(DLYRSI,AZMAX1(SnowDepth(NY,NX)-cumSnowDepthI(L-1)))
    ENDIF
    VOLSSL(L,NY,NX)=DLYRS(L,NY,NX)*DENS0(NY,NX)*DH(NY,NX)*DV(NY,NX)
    VOLWSL(L,NY,NX)=0.0_r8
    VOLISL(L,NY,NX)=0.0_r8

!    IF(L.EQ.1)THEN
!      VOLSWI=VOLSWI+0.5_r8*(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)*DENSI)
!    ELSE
!      VOLSWI=VOLSWI+0.5_r8*(VOLSSL(L-1,NY,NX)+VOLWSL(L-1,NY,NX) &
!        +VOLISL(L-1,NY,NX)*DENSI+VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
!        +VOLISL(L,NY,NX)*DENSI)
!    ENDIF

    DENSS(L,NY,NX)=DENS0(NY,NX)
    VOLSL(L,NY,NX)=VOLSSL(L,NY,NX)/DENSS(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)
    VOLSI(L,NY,NX)=DLYRSI*DH(NY,NX)*DV(NY,NX)      !it is a non-zero number, potential/maximum volume
    cumSnowDepth(L,NY,NX)=cumSnowDepth(L-1,NY,NX)+DLYRS(L,NY,NX)
    TKW(L,NY,NX)=AMIN1(Tref,TairKClimMean(NY,NX))
    TCW(L,NY,NX)=AZMIN1(ATCA(NY,NX))
    VHCPW(L,NY,NX)=cps*VOLSSL(L,NY,NX)+cpw*VOLWSL(L,NY,NX)+cpi*VOLISL(L,NY,NX)
  ENDDO D9580

!
!     VHCPWX,=minimum heat capacities for solving
!      snowpack water and heat fluxes
!
  VHCPWX(NY,NX)=VHCPWMin*AREA(3,NU(NY,NX),NY,NX)

  end subroutine InitSnowLayers

!------------------------------------------------------------------------------------------

  subroutine InitSnowAccums(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
!  TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
!  THFLWW=convective heat fluxes of snow,water,ice in snowpack
!

  D9875: DO L=1,JS
    TFLWS(L,NY,NX)=0.0_r8
    TFLWW(L,NY,NX)=0.0_r8
    TFLWI(L,NY,NX)=0.0_r8
    THFLWW(L,NY,NX)=0.0_r8
  ENDDO D9875

  end subroutine InitSnowAccums


!------------------------------------------------------------------------------------------

  subroutine SnowPackIteration(M,NY,NX,FLWLT,HFLWLT,FLWQGH,FLWRT,&
    HFLWRT,FLWHLW,FLWLW,FLWLXW,FLWRLW,HFLWRLW,HFLWLW)
  implicit none
  integer, intent(in) :: M,NY,NX  
  real(r8),intent(out) :: FLWLT,HFLWLT,FLWQGH,FLWRT,HFLWRT
  real(r8), intent(inout) :: FLWHLW,FLWLW,FLWLXW,FLWRLW,HFLWRLW,HFLWLW
  real(r8) :: TCND1W,TCNDR
  real(r8) :: ATCNDR,ATCNDS,ATCNDW,ATCNVW,ATCNVR,TCNDS
  real(r8) :: FLVC,FLVX,FLVR1
  real(r8) :: FLVSR,HFLVSR
  real(r8) :: FLVSS,HFLVSS
  real(r8) :: HFLVR1,HFLWSR
  real(r8) :: HFLWC,HFLWR1,HFLWX
  real(r8) :: PSISV1,VOLP01,VOLP02,VPY
  real(r8) :: TKY,THETP1,THETP2
  real(r8) :: FCDX,FCLX,FCX,FLVS1,HFLVS1
  real(r8) :: CNV1,CNV2,CNVR,DENSW1,DENSW2
  real(r8) :: FLW0T,HFLW0T,FLWQG,FLWQGS,FLWQGX
  real(r8) :: FLWQM,HFLWQM,FLWQR,HFLWQR,FLWQX,HFLWQG
  real(r8) :: HFLWS1,HFLWSS,PSDX,TCND2W,THETRR
  real(r8) :: TK0X,TK1X,VP1,VP2,WPLX,WPX
  integer :: L,L2,ICHKL

  ! begin_execution
  ! PHYSICAL AND HYDRAULIC PROPERTIES OF SNOWPACK INCLUDING
  ! AIR AND WATER-FILLED POROSITY, WATER POTENTIAL OF UNDERLYING
  ! SOIL SURFACE USED IN FLUX CALCULATIONS
  !
  ! VHCPW,VHCPWX=current, minimum snowpack heat capacities
  ! VOLS0M,VOLI0M,VOLW0M,VOLS1=snow,ice,water,total snowpack volume
  ! DENSS,DENSI,DENS0=snow,ice,minimum snow density
  ! AREA=area of grid cell
  ! VOLP01=snowpack air volume
  ! THETP1=snowpack air concentration
  ! CNV1=snowpack vapor conductivity
  ! VP1=snowpack vapor concentration
  ! TempK_snowpack=snowpack temperature
  ! WGSGW=vapor diffusivity
  ! DENSW1=snowpack density
  ! VHCPWMM=previous snowpack heat capacity

  ICHKL=0
  D9880: DO L=1,JS

    IF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
    ! active snow layer
      VOLS1(L,NY,NX)=VOLS0M(L,NY,NX)/DENSS(L,NY,NX)+VOLW0M(L,NY,NX)+VOLI0M(L,NY,NX)
      DLYRS0(L,NY,NX)=VOLS1(L,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
      VOLP01=AZMAX1(VOLS1(L,NY,NX)-VOLS0M(L,NY,NX)-VOLI0M(L,NY,NX)-VOLW0M(L,NY,NX))
      THETP1=AMAX1(THETPI,VOLP01/VOLS1(L,NY,NX))
      CNV1=THETP1**2.0*WGSGW(L,NY,NX)
      VP1=vapsat(TempK_snowpack(L,NY,NX))      
      IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        DENSW1=AMIN1(0.6_r8,(VOLS0M(L,NY,NX)+VOLW0M(L,NY,NX)+VOLI0M(L,NY,NX)*DENSI)/VOLS1(L,NY,NX))
      ELSE
        DENSW1=DENS0(NY,NX)
      ENDIF
      !
      ! SNOW THERMAL CONDUCTIVITY FROM J GLACIOL 43:26-41
      !
      ! TCND1W=snow thermal conductivity
      ! FLWQX=porosity-unconstrained snow water flux
      !
      TCND1W=0.0036_r8*10**(2.650_r8*DENSW1-1.652_r8)
      !
      ! DISCHARGE OF MELTWATER AND ITS HEAT FROM SNOWPACK LAYER
      ! TO LOWER SNOWPACK LAYER
      !
      FLWQX=AZMAX1(AZMAX1(VOLW0M(L,NY,NX))-0.05_r8*AZMAX1(VOLS0M(L,NY,NX)))*XNPA
      !
      ! WATER AND HEAT FLUXES IN SNOWPACK
      !
      ! DLYRS0=snow layer thickness
      ! FLWQM=porosity-constrained snow water flux
      ! HFLWQM=convective heat flux from water flux
      !
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWMM(L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        !if L==JS-1, L2==JS, so top layer is treated here.
        VOLS1(L2,NY,NX)=VOLS0M(L2,NY,NX)/DENSS(L2,NY,NX)+VOLW0M(L2,NY,NX)+VOLI0M(L2,NY,NX)
        DLYRS0(L2,NY,NX)=VOLS1(L2,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
        VOLP02=VOLS1(L2,NY,NX)-VOLS0M(L2,NY,NX)-VOLI0M(L2,NY,NX)-VOLW0M(L2,NY,NX)
        THETP2=AMAX1(THETPI,VOLP02/VOLS1(L2,NY,NX))
        FLWQM=AMIN1(THETP2,FLWQX)
        HFLWQM=cpw*TempK_snowpack(L,NY,NX)*FLWQM
        !
        ! VAPOR FLUX IN SNOWPACK
        !
        ! VOLP01,VOLP02=air-filled volumes of source, destination layers
        ! L2=destination layer
        ! CNV1,CNV2=vapor conductivities of source, destination layers
        ! VP1,VP2=vapor concentrations of source, destination layers
        ! TempK_snowpack=soil temperature
        ! ATCNVW=snow vapor conductance
        ! DLYRS0=snow layer thickness
        ! FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
        ! FLVSS,HFLVSS=vapor flux and its convective heat flux
        !
        IF(VOLP01.GT.ZEROS2(NY,NX).AND.VOLP02.GT.ZEROS2(NY,NX))THEN
          CNV2=THETP2**2.0_r8*WGSGW(L2,NY,NX)
          VP2=vapsat(TempK_snowpack(L2,NY,NX))
          ATCNVW=2.0_r8*CNV1*CNV2/(CNV1*DLYRS0(L2,NY,NX)+CNV2*DLYRS0(L,NY,NX))
          FLVC=ATCNVW*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPY
          VPY=(VP1*VOLP01+VP2*VOLP02)/(VOLP01+VOLP02)
          FLVX=(VP1-VPY)*VOLP01*XNPA
          IF(FLVC.GE.0.0_r8)THEN
            FLVSS=AZMAX1(AMIN1(FLVC,FLVX,VOLW0M(L,NY,NX)*XNPX))
            HFLVSS=(cpw*TempK_snowpack(L,NY,NX)+VAP)*FLVSS
          ELSE
            FLVSS=AZMIN1(AMAX1(FLVC,FLVX,-VOLW0M(L2,NY,NX)*XNPX))
            HFLVSS=(cpw*TempK_snowpack(L2,NY,NX)+VAP)*FLVSS
          ENDIF
        ELSE
          FLVSS=0.0_r8
          HFLVSS=0.0_r8
        ENDIF
        !
        ! HEAT FLUX IN SNOWPACK
        !
        ! DENSW2,TCNDW2=density,thermal conductivity in destination layer
        ! ATCNDW=thermal conductance
        ! DLYRS0=layer thickness
        ! TKY=equilibrium temperature
        ! HFLWX,HFLWC=heat-constrained,heat-unconstrained heat fluxes
        ! VHCPWMM,TempK_snowpack=volumetric heat capacity,temperature
        ! XNPX=time step for flux calculations
        ! FSNW=snow cover fraction
        ! XNPY=time step for snowpack flux calculations
        ! HFLWSS=snowpack heat flux
        ! FLW0S,FLQ0I,FLQ0W=snow,ice,water fluxes through snowpack
        ! HFLW0W=convective heat flux snow,water,ice fluxes
        !
        IF(VOLS1(L2,NY,NX).GT.ZEROS2(NY,NX))THEN
          DENSW2=AMIN1(0.6_r8,(VOLS0M(L2,NY,NX)+VOLW0M(L2,NY,NX) &
            +VOLI0M(L2,NY,NX)*DENSI)/VOLS1(L2,NY,NX))
        ELSE
          DENSW2=DENS0(NY,NX)
        ENDIF
        TCND2W=0.0036_r8*10**(2.650_r8*DENSW2-1.652_r8)
        ATCNDW=2.0_r8*TCND1W*TCND2W/(TCND1W*DLYRS0(L2,NY,NX)+TCND2W*DLYRS0(L,NY,NX))
        TKY=(TempK_snowpack(L,NY,NX)*VHCPWMM(L,NY,NX)+TempK_snowpack(L2,NY,NX) &
          *VHCPWMM(L2,NY,NX))/(VHCPWMM(L,NY,NX)+VHCPWMM(L2,NY,NX))
        HFLWX=(TempK_snowpack(L,NY,NX)-TKY)*VHCPWMM(L,NY,NX)*XNPA
        HFLWC=ATCNDW*(TempK_snowpack(L,NY,NX)-TempK_snowpack(L2,NY,NX))*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPY

        IF(HFLWC.GE.0.0_r8)THEN
          HFLWSS=AZMAX1(AMIN1(HFLWX,HFLWC))
        ELSE
          HFLWSS=AZMIN1(AMAX1(HFLWX,HFLWC))
        ENDIF
        FLW0T=FLWQM+FLVSS
        HFLW0T=HFLWQM+HFLVSS+HFLWSS
        FLW0S(L2,NY,NX)=0.0_r8
        FLW0W(L2,NY,NX)=FLW0T
        FLW0I(L2,NY,NX)=0.0_r8
        HFLW0W(L2,NY,NX)=HFLW0T
        FLQWM(M,L2,NY,NX)=FLQWM(M,L2,NY,NX)+FLWQM
        !
        ! DISCHARGE OF MELTWATER AND ITS HEAT FROM LOWEST SNOWPACK LAYER
        ! TO RESIDUE, SURFACE SOIL MICROPORES AND MACROPORES
        !
        ! FLWQX,FLWQR=porosity-unconstrained water flux to soil,litter
        ! FLWQGX,FLWQGS,FLWQGH=water flux to soil surface,
        ! micropores,macropores
        ! VOLP1,VOLPH1=air volumes of soil micropores,macropores
        ! FMAC,FGRD=macropore,micropore volume fractions
        ! HFLWQG,HFLWQR=convective heat fluxes to soil,litter

        ! THETY=hygroscopic water concentration
        ! POROS=soil porosity
        ! FC,WP,FCL,WPL=field capacity,wilting point, log(FC),log(WP)
        ! FCI,WPI=FC,WP of ice
        ! THETIX=ice concentration
        ! BKVL=bulk density x volume of soil layer
        ! VOLR=surface litter volume
        ! THETWX=relative pore fraciton filled by water
      ELSE
        !dealing with residual snow
        !L==JS, top layer
        IF(ICHKL.EQ.0)THEN

          !interaction with respect to litter layer and topsoil 
          FLWQGX=FLWQX*BARE(NY,NX)
          FLWQGS=AMIN1(VOLP1(NUM(NY,NX),NY,NX)*XNPX,FLWQGX*FGRD(NUM(NY,NX),NY,NX))
          FLWQGH=AMIN1(VOLPH1(NUM(NY,NX),NY,NX)*XNPX,FLWQGX*FMAC(NUM(NY,NX),NY,NX))
          FLWQG=FLWQGS+FLWQGH
          HFLWQG=cpw*TempK_snowpack(L,NY,NX)*FLWQG
          FLWQR=FLWQX-FLWQG
          HFLWQR=cpw*TempK_snowpack(L,NY,NX)*FLWQR

          call SnowTopSoilExch(M,L,NY,NX,CNV1,VP1,VOLP01,TCND1W,FLVS1,HFLVS1,HFLWS1,CNV2,PSISV1,TCNDS)

          !
          ! HEAT FLUX AMONG SNOWPACK, SURFACE RESIDUE AND SURFACE SOIL
          !
          ! FLVSR=snowpack-litter vapor flux
          ! HFLVSR,HFLWSR=snowpack-litter convective,conductive heat fluxes
          ! FLVS1=snowpack-soil vapor flux
          ! HFLVS1,HFLWS1=snowpack-soil convective,conductive heat fluxes
          ! VHCP1,VHCPRX=current,minimum litter heat capacities
          ! TK0X,TKXR,TK1X=snowpack,litter,soil temperatures
          ! CNVR,CNV1,CNV2=litter,snowpack,soil vapor conductivity
          ! THETP*,THETWX,THETIX=litter air,water,ice concentration
          ! POROS,POROQ=litter porosity, tortuosity
          ! CVRD=litter cover fraction
          ! WGSGR=litter vapor diffusivity
          ! ATCNVR,ATCNVS=snowpack-litter,litter-soil vapor conductance
          ! DLYRR,DLYRS0,DLYR=litter,snowpack,soil depths
          ! THETRR=dry litter concentration
          ! TCNDR,TCND1W,TCNDS=litter,snowpack,soil thermal conductivity
          ! ATCNDR,ATCNDS=snow-litter,litter-soil thermal conductance
          !
          FLVSR=0.0_r8
          HFLVSR=0.0_r8
          HFLWSR=0.0_r8
          FLVR1=0.0_r8
          HFLVR1=0.0_r8
          HFLWR1=0.0_r8

          !surface litter layer is active
          IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
            call SnowSurLitterExch(M,L,NY,NX,CNV1,CNV2,TCND1W,VOLP01,PSISV1,TCNDS,&
              FLVSR,HFLVSR,HFLWSR,FLVR1,HFLVR1,HFLWR1)
          ENDIF
          !
          ! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
          ! FOR LATER UPDATES TO STATE VARIABLES
          !
          ! FLWLT,FLWLW=total,accumulated water flux to soil micropores
          ! FLWLXW,FLWHLW=total,accumd snow-soil micropore,macropore water
          ! HFLWLT,HFLWLW=total,accumulated snow+litter heat flux to soil
          ! FLWRT,FLWRLW=total,accumulated snow+soil water flux to litter
          ! HFLWRT,HFLWRLW=total,accumulated snow+soil heat flux to litter
          ! FLQRM,FLQSM,FLQHM=total water flux to litter,soil micropore,macropore
          ! FLSW,FLSWH,FLSWR=water flux from lowest snow layer to soil macropore,micropore,litter
          ! HFLSW,HFLSWR=heat flux from lowest snow layer to soil,litter
!
          FLWLT=FLWQGS+FLVS1+FLVR1
          FLWLW=FLWLW+FLWLT
          if(abs(FLWLW)>1.e20_r8)then
            write(*,*)'FLWLW=',FLWQGS,FLVS1,FLVR1
            write(*,*)'at line',__LINE__
            call endrun(trim(mod_filename)//'at line',__LINE__)
          endif
          FLWLXW=FLWLXW+FLWQGS
          FLWHLW=FLWHLW+FLWQGH
          HFLWLT=HFLWQG+HFLVS1+HFLWS1+HFLVR1+HFLWR1
          HFLWLW=HFLWLW+HFLWLT
          FLWRT=FLWQR+FLVSR-FLVR1
          FLWRLW=FLWRLW+FLWRT
          HFLWRT=HFLWQR+HFLVSR+HFLWSR-HFLVR1-HFLWR1
          HFLWRLW=HFLWRLW+HFLWRT

          FLQRM(M,NY,NX)=FLQRM(M,NY,NX)+FLWQR
          FLQSM(M,NY,NX)=FLQSM(M,NY,NX)+FLWQGS
          FLQHM(M,NY,NX)=FLQHM(M,NY,NX)+FLWQGH
          FLSW(L,NY,NX)=FLSW(L,NY,NX)+FLWLT
          FLSWH(L,NY,NX)=FLSWH(L,NY,NX)+FLWQGH
          HFLSW(L,NY,NX)=HFLSW(L,NY,NX)+HFLWLT
          FLSWR(L,NY,NX)=FLSWR(L,NY,NX)+FLWRT
          HFLSWR(L,NY,NX)=HFLSWR(L,NY,NX)+HFLWRT

          ICHKL=1
        ENDIF
      ENDIF
    ENDIF
  ENDDO D9880
  end subroutine SnowPackIteration

!------------------------------------------------------------------------------------------

  subroutine SolveSnowpack(M,NY,NX,EFLXW,RFLXW,VFLXW,SFLXW,HFLXW,FLWHLW,&
    FLWLW,FLWLXW,FLWRLW,HFLWRLW,HFLWLW)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(inout) :: EFLXW,RFLXW,VFLXW,SFLXW,HFLXW,FLWHLW,FLWLW,FLWLXW
  real(r8), intent(inout) :: FLWRLW,HFLWRLW,HFLWLW
  integer :: MM,L,L2
  real(r8):: Raa,ALBW
  real(r8) :: FLWLT,HFLWLT  
  real(r8) :: ENGY0
  real(r8) :: TK0X,TKXR,TK1X
  real(r8) :: VP0
  real(r8) :: EFLXW2,EVAPS2,EVAPT2,EVAPW2,EVAPX2
  real(r8) :: FLW0W2,HFLW0W2,FLWQGH,FLWRT,HFLWRT
  real(r8) :: HWFLQ02,FLW0S2,FLW0I2,HFLXW2,RFLXW2
  real(r8) :: FLQ0I2,FLQ0S2,FLQ0W2,FVOLI0,FVOLS0
  real(r8) :: SFLXW2,TFLWIX,TFLWSX,TFLWWX,TFLX0X,TFLX1
  real(r8) :: HFLX02,PARE,PARS,RAGX,RFLX0,RI
  real(r8) :: THFLWWX,THRMX,TVOLWS,VFLXW2,VHCPWMX
  real(r8) :: VOLI0X,VOLS0X,VOLW0X,WFLXIX,WFLXSX
  !     begin_execution
  !     HEAT AND VAPOR FLUXES BETWEEN SNOWPACK AND ATMOSPHERE
  !
  !     VHCPWM=volumetric heat capacity of snowpack
  !     NPS=number of cycles for solving snowpack heat and water fluxes
  !     ALBW=snowpack albedo
  !     VOLS0M,VOLI0M,VOLW0M=snow,ice,water volumes
  !     RFLX0=net radiation input
  !     RADXW=shortwave radiation at snowpack surface
  !     THRYW=longwave radn incident at snowpack surface
  !     THRMX=longwave radn emitted by snowpack surface
  !     TempK_snowpack=snowpack surface temperature
  !     RFLXW2=net radiation

  D3000: DO MM=1,NPS

    ALBW=(0.85_r8*VOLS0M(1,NY,NX)+0.30_r8*VOLI0M(1,NY,NX)+0.06_r8*VOLW0M(1,NY,NX)) &
      /(VOLS0M(1,NY,NX)+VOLI0M(1,NY,NX)+VOLW0M(1,NY,NX))

    RFLX0=(1.0_r8-ALBW)*RADXW(NY,NX)+THRYW(NY,NX)    !incoming radiation, short + longwave
    THRMX=THRMW(NY,NX)*TempK_snowpack(1,NY,NX)**4           !emitting longwave radiation,
    RFLXW2=RFLX0-THRMX                            !net radiation
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
    RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB(NY,NX)*(TKQ(NY,NX)-TempK_snowpack(1,NY,NX))))
    RAGX=AMAX1(RAM,0.8_r8*RAGW(NY,NX),AMIN1(1.2_r8*RAGW(NY,NX),RAG(NY,NX)/(1.0_r8-10.0_r8*RI)))
    RAGW(NY,NX)=RAGX
    RAa=RAGX
    !
    ! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
    !
    !     PARE,PARS=blcs for snowpack latent,sensible heat fluxes
    !     PAREW,PARSW=conductances for latent,sensible heat fluxes
    !     RZ=surface resistance
    !     VP0,VPQ=vapor pressure at snowpack surface, canopy air
    !     EVAPT2,EVAPW2,EVAPS2=evaporation total, water,snow
    !     XNPS=1/NPS
    !     EFLXW2=latent heat flux
    !     VAP,VAPS=latent heat of evaporation,sublimation
    !     VFLXW2=convective heat of evaporation flux
    !
    PARE=PAREW(NY,NX)/(RAa+RZ)
    PARS=PARSW(NY,NX)/RAa    
    VP0=vapsat(TempK_snowpack(1,NY,NX))
    EVAPT2=PARE*(VPQ(NY,NX)-VP0)
    EVAPW2=AMAX1(EVAPT2,-AZMAX1(VOLW0M(1,NY,NX)*XNPA))
    EVAPX2=AZMIN1(EVAPT2-EVAPW2)
    EVAPS2=AMAX1(EVAPX2,-AZMAX1(VOLS0M(1,NY,NX)*XNPA))
    EFLXW2=EVAPW2*VAP+EVAPS2*VAPS

    IF(EVAPT2.LT.0.0_r8)THEN
      VFLXW2=(EVAPW2*cpw+EVAPS2*cps)*TempK_snowpack(1,NY,NX)
    ELSE
      VFLXW2=(EVAPW2*cpw+EVAPS2*cps)*TKQ(NY,NX)
    ENDIF
!
!     SOLVE FOR SNOWPACK SURFACE TEMPERATURE AT WHICH ENERGY
!     BALANCE OCCURS, SOLVE AND ACCUMULATE LATENT, SENSIBLE
!     STORAGE HEAT FLUXES AND EVAPORATION
!
!     SFLXW2,EFLXW2,RFLXW2=sensible,latent heat fluxes, net radiation
!     VFLXW2=convective heat flux from EFLXW2
!     HFLX02=storage heat flux
!     FLQ0S2,FLQ0W2,FLQ0I2=snow,water,ice input to snowpack
!     HWFLQ02=convective heat from snow,water,ice input to snowpack
!
    SFLXW2=PARS*(TKQ(NY,NX)-TempK_snowpack(1,NY,NX))
    !occasionally, RFLXW2 and SFLXW2 go to infinity
    HFLX02=RFLXW2+EFLXW2+SFLXW2   
    HFLXW2=HFLX02+VFLXW2
    RFLXW=RFLXW+RFLXW2
    EFLXW=EFLXW+EFLXW2
    VFLXW=VFLXW+VFLXW2
    SFLXW=SFLXW+SFLXW2
    HFLXW=HFLXW+HFLXW2
    EVAPS(NY,NX)=EVAPS(NY,NX)+EVAPS2
    EVAPW(NY,NX)=EVAPW(NY,NX)+EVAPW2
    EVAPSN(NY,NX)=EVAPSN(NY,NX)+EVAPS2+EVAPW2
    FLQ0S2=FLQ0S(NY,NX)*XNPS
    FLQ0W2=FLQ0W(NY,NX)*XNPS
    FLQ0I2=FLQ0I(NY,NX)*XNPS
    HWFLQ02=HWFLQ0(NY,NX)*XNPS
    FLW0S2=FLQ0S2+EVAPS2
    FLW0W2=FLQ0W2+EVAPW2
    FLW0I2=FLQ0I2
    HFLW0W2=HWFLQ02+HFLXW2
    FLW0S(1,NY,NX)=FLW0S2
    FLW0W(1,NY,NX)=FLW0W2
    FLW0I(1,NY,NX)=FLW0I2
    HFLW0W(1,NY,NX)=HFLW0W2
    FLQWM(M,1,NY,NX)=FLQWM(M,1,NY,NX)+FLQ0S2+FLQ0I2+FLQ0W2
!     IF(NX.EQ.3.AND.NY.EQ.3)THEN
!     WRITE(*,7759)'EVAP',I,J,M,MM,FLW0S2
!    2,FLQ0S2,EVAPS2,FLW0W2,FLQ0W2
!    3,FSNW(NY,NX),FLW0I2,FLQ0I2,RFLXW2,EFLXW2
!    4,SFLXW2,VFLXW2,RA,EVAPT2,EVAPX2,VPQ(NY,NX),VP0
!    5,VOLW0M(1,NY,NX),VOLS0M(1,NY,NX),VOLI0M(1,NY,NX)
!    6,HFLW0W(1,NY,NX),HFLW0W2,HWFLQ02,HFLXW2,RFLXW2,EFLXW2
!    7,SFLXW2,VFLXW2,TempK_snowpack(1,NY,NX),TKQ(NY,NX)
!    8,PARE,RA,RZ,EVAPS2,EVAPW2,EVAPT2
!7759  FORMAT(A8,4I4,40E14.6)
!     ENDIF
!
    call SnowPackIteration(M,NY,NX,FLWLT,HFLWLT,FLWQGH,FLWRT,HFLWRT,FLWHLW,FLWLW,&
      FLWLXW,FLWRLW,HFLWRLW,HFLWLW)
!
!     ACCUMULATE SNOWPACK FLUXES TO LONGER TIME STEP FOR
!     LITTER, SOIL FLUX CALCULATIONS
!
!     THRMG=total longwave emission
!     XFLWS,XFLWW,XFLWI=hourly accumulated snow,water,ice transfer
!     XHFLWW=hourly convective heat flux from snow,water,ice transfer
!     TFLWSX,TFLWWX,TFLWIX=net snow,water,ice transfer
!     THFLWWX=convective heat flux from net snow,water,ice transfer
!     TFLWS,TFLWW,TFLWI=accumulated net snow,water,ice transfer
!     THFLWW=convective heat flux from accumd snow,water,ice transfer
!
    THRMG(NY,NX)=THRMG(NY,NX)+THRMX
    D9860: DO L=1,JS
      XFLWS(L,NY,NX)=XFLWS(L,NY,NX)+FLW0S(L,NY,NX)
      XFLWW(L,NY,NX)=XFLWW(L,NY,NX)+FLW0W(L,NY,NX)
      XFLWI(L,NY,NX)=XFLWI(L,NY,NX)+FLW0I(L,NY,NX)
      XHFLWW(L,NY,NX)=XHFLWW(L,NY,NX)+HFLW0W(L,NY,NX)
      L2=MIN(JS,L+1)
!
!     IF WITHIN SNOWPACK
!
      IF(L.LT.JS.AND.VHCPWMM(L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        TFLWSX=FLW0S(L,NY,NX)-FLW0S(L2,NY,NX)
        TFLWWX=FLW0W(L,NY,NX)-FLW0W(L2,NY,NX)
        TFLWIX=FLW0I(L,NY,NX)-FLW0I(L2,NY,NX)
        THFLWWX=HFLW0W(L,NY,NX)-HFLW0W(L2,NY,NX)
        TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX
        TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX
        TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
        THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
!
!     IF AT BOTTOM OF SNOWPACK
!
      ELSEIF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
        TFLWSX=FLW0S(L,NY,NX)
        TFLWWX=FLW0W(L,NY,NX)-FLWRT-FLWLT-FLWQGH
        TFLWIX=FLW0I(L,NY,NX)
        THFLWWX=HFLW0W(L,NY,NX)-HFLWRT-HFLWLT
        TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX
        TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX
        TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
        THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
      ELSE
        TFLWSX=0.0_r8
        TFLWWX=0.0_r8
        TFLWIX=0.0_r8
        THFLWWX=0.0_r8
      ENDIF
!
!     FREEZE-THAW IN SNOWPACK FROM NET CHANGE IN SNOWPACK
!     HEAT STORAGE
!
!     VOLS0M,VOLW0M,VOLI0M=snow,water,ice volume
!     VHCPWMM,VHCPWMX,VHCPWX=previous,current,minimum heat capacity
!     THFLWWX=net conductive+convective heat flux
!     TFLX1=unconstrained latent heat flux from freeze-thaw
!     FVOLS0,FVOLI0=fractions of total water in water,ice
!     TFLX0X=source-limited latent heat flux from freeze-thaw
!     WFLXSX,WFLXIX=freeze-thaw changes in water,ice
!     WFLXS,WFLXI=accumulated freeze-thaw
!     TFLX0=accumulated latent heat flux from freeze-thaw
!     XWFLXS,XWFLXI=hourly accumulated freeze-thaw
!     XTHAWW=hourly accumulated latent heat flux from freeze-thaw
!
      VOLS0X=AZMAX1(VOLS0M(L,NY,NX))
      VOLW0X=AZMAX1(VOLW0M(L,NY,NX))
      VOLI0X=AZMAX1(VOLI0M(L,NY,NX))
      ENGY0=VHCPWMM(L,NY,NX)*TempK_snowpack(L,NY,NX)
      VHCPWMX=cps*VOLS0X+cpw*VOLW0X+cpi*VOLI0X
      IF(VHCPWMX.GT.VHCPWX(NY,NX))THEN
        TK0X=(ENGY0+THFLWWX)/VHCPWMX
        IF((TK0X.LT.TFice.AND.VOLW0X.GT.ZERO*VOLS(NY,NX)) &
          .OR.(TK0X.GT.TFice.AND.VOLI0X+VOLS0X.GT.ZERO*VOLS(NY,NX)))THEN
          TFLX1=VHCPWMX*(TFice-TK0X)/2.7185*XNPX
          IF(TFLX1.LT.0.0_r8)THEN
            TVOLWS=VOLS0X+VOLI0X*DENSI
            IF(TVOLWS.GT.ZEROS2(NY,NX))THEN
              FVOLS0=VOLS0X/TVOLWS
              FVOLI0=VOLI0X*DENSI/TVOLWS
            ELSE
              FVOLS0=0.0_r8
              FVOLI0=0.0_r8
            ENDIF
            TFLX0X=AMAX1(-333.0*TVOLWS*XNPX,TFLX1)
            WFLXSX=-TFLX0X*FVOLS0/333.0
            WFLXIX=-TFLX0X*FVOLI0/333.0
          ELSE
            FVOLS0=0.0_r8
            FVOLI0=0.0_r8
            TFLX0X=AMIN1(333.0*VOLW0X*XNPX,TFLX1)
            WFLXSX=0.0_r8
            WFLXIX=-TFLX0X/333.0
          ENDIF
        ELSE
          TFLX1=0.0_r8
          FVOLS0=0.0_r8
          FVOLI0=0.0_r8
          TFLX0X=0.0_r8
          WFLXSX=0.0_r8
          WFLXIX=0.0_r8
        ENDIF
        WFLXS(L,NY,NX)=WFLXS(L,NY,NX)+WFLXSX
        WFLXI(L,NY,NX)=WFLXI(L,NY,NX)+WFLXIX
        TFLX0(L,NY,NX)=TFLX0(L,NY,NX)+TFLX0X
        XWFLXS(L,NY,NX)=XWFLXS(L,NY,NX)+WFLXSX
        XWFLXI(L,NY,NX)=XWFLXI(L,NY,NX)+WFLXIX
        XTHAWW(L,NY,NX)=XTHAWW(L,NY,NX)+TFLX0X
      ELSE
        TFLX0X=0.0_r8
        WFLXSX=0.0_r8
        WFLXIX=0.0_r8
      ENDIF
!
!     INTERNAL SNOWPACK SNOW, WATER, ICE, TEMPERATURE
!
!     VOLS0M,VOLW0M,VOLI0M=snow water eqv,water,ice volume
!     TFLWSX,TFLWWX,TFLWIX=net snow,water,ice transfer
!     THFLWWX=conductive+convective heat from snow,water,ice transfer
!     WFLXSX,WFLXIX=freeze-thaw changes in water,ice
!     TFLX0X=source-limited latent heat flux from freeze-thaw
!     DENSI=ice density
!     TempK_snowpack,TairK=snowpack,air temperature
!     VHCPWMM,VHCPWX=snowpack, minimum heat capacity
!
      VOLS0M(L,NY,NX)=VOLS0M(L,NY,NX)+TFLWSX-WFLXSX
      VOLW0M(L,NY,NX)=VOLW0M(L,NY,NX)+TFLWWX+WFLXSX+WFLXIX
      VOLI0M(L,NY,NX)=VOLI0M(L,NY,NX)-WFLXIX/DENSI
      ENGY0=VHCPWMM(L,NY,NX)*TempK_snowpack(L,NY,NX)
      VHCPWMM(L,NY,NX)=cps*VOLS0M(L,NY,NX)+cpw*VOLW0M(L,NY,NX)+cpi*VOLI0M(L,NY,NX)
      IF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
        TempK_snowpack(L,NY,NX)=(ENGY0+THFLWWX+TFLX0X)/VHCPWMM(L,NY,NX)
      ELSEIF(L.EQ.1)THEN
        TempK_snowpack(L,NY,NX)=TairK(NY,NX)
      ELSE
        TempK_snowpack(L,NY,NX)=TempK_snowpack(L-1,NY,NX)
      ENDIF
    ENDDO D9860
  ENDDO D3000
  end subroutine SolveSnowpack  
!------------------------------------------------------------------------------------------
  subroutine PrepIterSnowLayer(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  integer :: L
! FLUX VARIABLES IN SNOWPACK
!
! TFLX0=latent heat from freeze-thaw
! WFLXS,WFLXI=freeze-thaw between snow,ice and water
! FLW0S,FLW0W,FLW0I=snow,water,ice fluxes
! HFLW0W=convective heat flux from snow,water,ice fluxes
! FLQWM=snowpack water flux
! VOLS0M,VOLW0M,VOLI0M=snow,water,ice contents
! VHCPWMM=volumetric heat capacity
! TempK_snowpack=snow temperature
!
  D9765: DO L=1,JS
    TFLX0(L,NY,NX)=0.0_r8
    WFLXS(L,NY,NX)=0.0_r8
    WFLXI(L,NY,NX)=0.0_r8
    FLW0S(L,NY,NX)=0.0_r8
    FLW0W(L,NY,NX)=0.0_r8
    FLW0I(L,NY,NX)=0.0_r8
    HFLW0W(L,NY,NX)=0.0_r8
    FLQWM(M,L,NY,NX)=0.0_r8
    VOLS0M(L,NY,NX)=VOLS0(L,NY,NX)
    VOLW0M(L,NY,NX)=VOLW0(L,NY,NX)
    VOLI0M(L,NY,NX)=VOLI0(L,NY,NX)
    VHCPWMM(L,NY,NX)=VHCPWM(M,L,NY,NX)    
    TempK_snowpack(L,NY,NX)=TK0(L,NY,NX)
  ENDDO D9765
  end subroutine PrepIterSnowLayer  
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
!
!     VOLS0,VOLW0,VOLI0=snow,water,ice volumes in snowpack
!     TFLWS,TFLWW=net snow,water flux
!     WFLXS,WFLXI=snow-water,ice-water freeze-thaw flux
!     DENSI=ice density
!     VHCPWM=snowpack volumetric heat capacity
!     TK0=snowpack temperature
!     THFLWW=total snowpack conductive+convective heat flux
!     TFLX0=snowpack latent heat flux from freeze-thaw
!
  D9780: DO L=1,JS
    VOLS0(L,NY,NX)=VOLS0(L,NY,NX)+TFLWS(L,NY,NX)-WFLXS(L,NY,NX)
    VOLW0(L,NY,NX)=VOLW0(L,NY,NX)+TFLWW(L,NY,NX)+WFLXS(L,NY,NX)+WFLXI(L,NY,NX)
    VOLI0(L,NY,NX)=VOLI0(L,NY,NX)-WFLXI(L,NY,NX)/DENSI
    ENGY0=VHCPWM(M,L,NY,NX)*TK0(L,NY,NX)
    VHCPWM(M+1,L,NY,NX)=cps*VOLS0(L,NY,NX)+cpw*VOLW0(L,NY,NX)+cpi*VOLI0(L,NY,NX)
    IF(VHCPWM(M+1,L,NY,NX).GT.VHCPWX(NY,NX))THEN
      TK0(L,NY,NX)=(ENGY0+THFLWW(L,NY,NX)+TFLX0(L,NY,NX))/VHCPWM(M+1,L,NY,NX)
    ELSEIF(L.EQ.1)THEN
      TK0(L,NY,NX)=TairK(NY,NX)
    ELSE
      TK0(L,NY,NX)=TK0(L-1,NY,NX)
    ENDIF
  ENDDO D9780
    !      if(curday>=176)then
    !        write(*,*)'line',__LINE__,'tk1',TK1(8,ny,nx),TK1(9,ny,nx),M
    !      endif
    !
    !     SNOW RUNOFF
    !
    !     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
    !
  VOLS0(1,NY,NX)=VOLS0(1,NY,NX)+TQS1(NY,NX)
  VOLW0(1,NY,NX)=VOLW0(1,NY,NX)+TQW1(NY,NX)
  VOLI0(1,NY,NX)=VOLI0(1,NY,NX)+TQI1(NY,NX)
  ENGY0=VHCPWM(M+1,1,NY,NX)*TK0(1,NY,NX)
  VHCPWM(M+1,1,NY,NX)=cps*VOLS0(1,NY,NX)+cpw*VOLW0(1,NY,NX)+cpi*VOLI0(1,NY,NX)
  IF(VHCPWM(M+1,1,NY,NX).GT.VHCPWX(NY,NX))THEN
    TK0(1,NY,NX)=(ENGY0+THQS1(NY,NX))/VHCPWM(M+1,1,NY,NX)
  ELSE
    TK0(1,NY,NX)=TairK(NY,NX)
  ENDIF
    !
    !     IF SNOWPACK DISAPPEARS ALL MATERIAL,HEAT TRANSFERRED TO SOIL SURFACE
    !
    !     VHCP1,VHCP1A,VHCP1P=total soil,soil+micropore,macropore heat capacity
    !     TK1=soil surface temperature, why not to litter layer
    !
  IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX).AND.TairK(NY,NX).GT.TFice)THEN
    FLWS=VOLS0(1,NY,NX)
    FLWW=VOLW0(1,NY,NX)
    FLWI=VOLI0(1,NY,NX)
    HFLWS=(cpw*FLWW+cps*FLWS+cpi*FLWI)*TK0(1,NY,NX)
    VOLS0(1,NY,NX)=VOLS0(1,NY,NX)-FLWS
    VOLW0(1,NY,NX)=VOLW0(1,NY,NX)-FLWW
    VOLI0(1,NY,NX)=VOLI0(1,NY,NX)-FLWI

    !add to litter layer?
    VOLW1(NUM(NY,NX),NY,NX)=VOLW1(NUM(NY,NX),NY,NX)+FLWW
    VOLI1(NUM(NY,NX),NY,NX)=VOLI1(NUM(NY,NX),NY,NX)+FLWI+FLWS/DENSI
    
    ENGY1=VHCP1(NUM(NY,NX),NY,NX)*TK1(NUM(NY,NX),NY,NX)
    VHCP1A(NUM(NY,NX),NY,NX)=VHCM(NUM(NY,NX),NY,NX) &
      +cpw*VOLW1(NUM(NY,NX),NY,NX)+cpi*VOLI1(NUM(NY,NX),NY,NX)
    VHCP1B(NUM(NY,NX),NY,NX)=cpw*VOLWH1(NUM(NY,NX),NY,NX) &
      +cpi*VOLIH1(NUM(NY,NX),NY,NX)
    VHCP1(NUM(NY,NX),NY,NX)=VHCP1A(NUM(NY,NX),NY,NX)+VHCP1B(NUM(NY,NX),NY,NX)

    IF(VHCP1(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    ! topsoil layer is there
      tk1pres=TK1(NUM(NY,NX),NY,NX)
      TK1(NUM(NY,NX),NY,NX)=(ENGY1+HFLWS)/VHCP1(NUM(NY,NX),NY,NX)
      if(abs(tk1pres/TK1(NUM(NY,NX),NY,NX)-1._r8)>0.025_r8)then
        TK1(NUM(NY,NX),NY,NX)=TairK(NY,NX)
      endif
    ELSE
      TK1(NUM(NY,NX),NY,NX)=TairK(NY,NX)
    ENDIF
  ENDIF
  end subroutine UpdateSnowAtM

!------------------------------------------------------------------------------------------
  subroutine SnowRedistrub(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
!
! SNOW redistribution
  implicit none
  integer, intent(in) :: M,NY,NX,NHE,NHW,NVS,NVN
  integer, intent(in) :: N1,N2

  integer :: N,NN,N4,N5,N4B,N5B
  real(r8) :: ALTS1,ALTS2
  real(r8) :: QSX,SS
  integer, parameter :: idirew=1
  integer, parameter :: idirns=2
!
!     SNOW REDISTRIBUTION FROM SNOWPACK
!
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     ALTS1,ALTS2=elevation of source,destination snowpack surfaces
!     SS,DIST=slope,distance between source,destination
!     QSX=transfer fraction
!     QS1,QW1,QI1=snow,water,ice transfer
!     HQS1=convective heat transfer from snow,water,ice transfer
!     VOLS0,VOLW0,VOLI0=snow,water,ice volume
!     MinSnowDepth=minimum snowpack depth for full cover
!     QS,QW,QI=hourly-accumulated snow,water,ice transfer
!     HQS=hourly-accumd convective heat from snow,water,ice transfer
!     QSM=snow transfer for solute flux calculation

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

      IF(NN.EQ.1)THEN
        !east or south
        ALTS1=ALTG(N2,N1)+SnowDepth(N2,N1)
        ALTS2=ALTG(N5,N4)+SnowDepth(N5,N4)
        SS=(ALTS1-ALTS2)/DIST(N,NU(N5,N4),N5,N4)
        QSX=SS/AMAX1(1.0_r8,DIST(N,NU(N5,N4),N5,N4))*XNPH
        IF(SS.GT.0.0_r8.AND.SnowDepth(N2,N1).GT.MinSnowDepth)THEN
          QS1(N,N5,N4)=QSX*VOLS0(1,N2,N1)
          QW1(N,N5,N4)=QSX*VOLW0(1,N2,N1)
          QI1(N,N5,N4)=QSX*VOLI0(1,N2,N1)
          HQS1(N,N5,N4)=TK0(1,N2,N1)*(cps*QS1(N,N5,N4)+cpw*QW1(N,N5,N4)+cpi*QI1(N,N5,N4))
        ELSEIF(SS.LT.0.0_r8.AND.SnowDepth(N5,N4).GT.MinSnowDepth)THEN
          QS1(N,N5,N4)=QSX*VOLS0(1,N5,N4)
          QW1(N,N5,N4)=QSX*VOLW0(1,N5,N4)
          QI1(N,N5,N4)=QSX*VOLI0(1,N5,N4)
          HQS1(N,N5,N4)=TK0(1,N5,N4)*(cps*QS1(N,N5,N4)+cpw*QW1(N,N5,N4)+cpi*QI1(N,N5,N4))
        ELSE
          QS1(N,N5,N4)=0.0_r8
          QW1(N,N5,N4)=0.0_r8
          QI1(N,N5,N4)=0.0_r8
          HQS1(N,N5,N4)=0.0_r8
        ENDIF
        QS(N,N5,N4)=QS(N,N5,N4)+QS1(N,N5,N4)
        QW(N,N5,N4)=QW(N,N5,N4)+QW1(N,N5,N4)
        QI(N,N5,N4)=QI(N,N5,N4)+QI1(N,N5,N4)
        HQS(N,N5,N4)=HQS(N,N5,N4)+HQS1(N,N5,N4)
        QSM(M,N,N5,N4)=QS1(N,N5,N4)
      ENDIF
      !add west and south
      IF(NN.EQ.2)THEN

      ENDIF
    ENDDO
  ENDDO    
  end subroutine SnowRedistrub

!------------------------------------------------------------------------------------------
  subroutine SumSnowRoffDrift(M,N,N1,N2,N4,N5,N4B,N5B)
  use SoilWaterDataType, only : IFLBM
  implicit none 
  integer, intent(in) :: M,N,N1,N2,N4,N5,N4B,N5B
  integer :: NN

  D1202: DO NN=1,2
    TQR1(N2,N1)=TQR1(N2,N1)+QR1(N,NN,N2,N1)
    THQR1(N2,N1)=THQR1(N2,N1)+HQR1(N,NN,N2,N1)
    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
    !there is runoff
      TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5,N4)
      THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5,N4)
    ENDIF

    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5B,N4B)
      THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5B,N4B)
    ENDIF

    IF(M.EQ.NPH)THEN
      IFLBH(N,NN,N5,N4)=IFLBM(M,N,NN,N5,N4)
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        IFLBH(N,NN,N5B,N4B)=IFLBM(M,N,NN,N5B,N4B)
      ENDIF
    ENDIF
  ENDDO D1202

  TQS1(N2,N1)=TQS1(N2,N1)+QS1(N,N2,N1)-QS1(N,N5,N4)
  TQW1(N2,N1)=TQW1(N2,N1)+QW1(N,N2,N1)-QW1(N,N5,N4)
  TQI1(N2,N1)=TQI1(N2,N1)+QI1(N,N2,N1)-QI1(N,N5,N4)
  THQS1(N2,N1)=THQS1(N2,N1)+HQS1(N,N2,N1)-HQS1(N,N5,N4)
  END subroutine SumSnowRoffDrift
!------------------------------------------------------------------------------------------
  subroutine ZeroSnowFlux(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  TQR1(NY,NX)=0.0_r8
  THQR1(NY,NX)=0.0_r8
  TQS1(NY,NX)=0.0_r8
  TQW1(NY,NX)=0.0_r8
  TQI1(NY,NX)=0.0_r8
  THQS1(NY,NX)=0.0_r8
  end subroutine ZeroSnowFlux  
!------------------------------------------------------------------------------------------
  FUNCTION SnowBNDResistance(NY,NX)result(RAS)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: RAS
  real(r8) :: RASL,RASX
  real(r8) :: THETPL
  integer :: L

!     RAS,RASL=blrs of snowpack,snowpack layer
!     VOLS,VOLS1=volume of snowpack,snowpack layer
!     DLYRS=snowpack later depth
!     WGSGW=vapor diffusivity in snowpack
!     THETPL=snowpack air-filled porosity
!     THETPI=air content of ice
!     VOLS0,VOLI0,VOLW0,VOLS1=snow,ice,water,total volumes of snowpack

  RAS=0.0_r8
  IF(VOLS(NY,NX).GT.ZEROS2(NY,NX))THEN
    D9775: DO L=1,JS
      IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        RASX=DLYRS(L,NY,NX)/WGSGW(L,NY,NX)
        THETPL=AMAX1(THETPI,1.0_r8-(VOLS0(L,NY,NX)+VOLI0(L,NY,NX) &
            +VOLW0(L,NY,NX))/VOLS1(L,NY,NX))
        RASL=RASX/AMAX1(ZERO,THETPL)**2.0_r8
        RAS=RAS+RASL
      ENDIF
    ENDDO D9775
  ENDIF
  END FUNCTION SnowBNDResistance      

!------------------------------------------------------------------------------------------

  subroutine SnowSurLitterIterate(L,M,NY,NX,ATCNDR,ATCNDS,ATCNVS,ATCNVR,&
    PSISV1,VOLP01,TK0X,TKXR,TK1X,FLVSR,HFLVR1,HFLVSR,HFLWSR,HFLWR1)
  implicit none
  integer, intent(in) :: L,M,NY,NX
  real(r8), intent(in) :: ATCNDR,ATCNDS,ATCNVR,ATCNVS,PSISV1,VOLP01
  real(r8), intent(inout) :: TK0X,TKXR,TK1X,FLVSR
  real(r8), intent(inout) :: HFLVR1,HFLVSR,HFLWSR,HFLWR1
  integer :: NN
  real(r8) :: VP1,VPR,VP0,VPY
  real(r8) :: FLVC,FLVR1,FLVR1X
  real(r8) :: HFLVR1X
  real(r8) :: FLVSRX,HFLVSRX,HFLWSRX
  real(r8) :: FLVX,HFLWC,HFLWR1X
  real(r8) :: HFLWX,TKY

! begin_execution
  D4000: DO NN=1,NPR
    !
    ! VAPOR FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! VP0,VPR,VPY=snowpack,litter, equilibrium vapor concentration
    ! TK0X,TKXR=snowpack,litter temperature
    ! PSISM1=litter matric water potential
    ! FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
    ! AREA=area of grid cell
    ! FSNW,CVRD=snow,litter cover fraction
    ! XNPQ=time step for flux calculation
    ! FLVRSX=snow-litter vapor flux
    ! HFLVRSX=convective heat flux from snow-litter vapor flux
    !
    !residue vapor pressure
    VPR=vapsat(TKXR)*EXP(18.0_r8*PSISM1(0,NY,NX)/(RGAS*TKXR))
    if(abs(VPR)>1.e20_r8)then
      write(*,*)'TKXR=',TKXR,TK1(0,NY,NX),TK1(NUM(NY,NX),NY,NX)
      write(*,*)'PSISM1(0,NY,NX)=',PSISM1(0,NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif

    !ATCNVR=snowpack-litter conductance
    !there is litter layer with sufficient moisture
    IF(VOLP01.GT.ZEROS2(NY,NX).AND.THETPM(M,0,NY,NX).GT.THETX)THEN      
      VP0=vapsat(TK0X) !snow vapor pressure, saturated
      !snow <-> residue vapor flux
      FLVC=ATCNVR*(VP0-VPR)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ 
      !volume weighted vapor pressure
      VPY=(VP0*VOLP01+VPR*VOLPM(M,0,NY,NX))/(VOLP01+VOLPM(M,0,NY,NX))             
      FLVX=(VP0-VPY)*VOLP01*XNPC

      IF(FLVC.GE.0.0_r8)THEN
        !vapor flux from snow to litter
        FLVSRX=AZMAX1(AMIN1(FLVC,FLVX))       !water flux
        HFLVSRX=(cpw*TK0X+VAP)*FLVSRX         !enthalpy flux associated with water flux
      ELSE
        !vapor flux from litter to snow
        FLVSRX=AZMIN1(AMAX1(FLVC,FLVX))
        HFLVSRX=(cpw*TKXR+VAP)*FLVSRX
      ENDIF
    ELSE
      FLVSRX=0.0_r8
      HFLVSRX=0.0_r8
    ENDIF
    !
    ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! TKY=snow-litter equilibrium temperature
    ! HFLWC,HFLWX=snow-litter heat flux unltd,ltd by heat
    ! HFLWSRX=snow-litter heat flux
    ! VHCPWMM= volumetric heat capacity in snow layer
    TKY=(TK0X*VHCPWMM(L,NY,NX)+TKXR*VHCP1(0,NY,NX))/(VHCPWMM(L,NY,NX)+VHCP1(0,NY,NX))
    HFLWX=(TK0X-TKY)*VHCPWMM(L,NY,NX)*XNPC
    HFLWC=ATCNDR*(TK0X-TKXR)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
    IF(HFLWC.GE.0.0_r8)THEN
      HFLWSRX=AZMAX1(AMIN1(HFLWX,HFLWC))
    ELSE
      HFLWSRX=AZMIN1(AMAX1(HFLWX,HFLWC))
    ENDIF
!
!     VAPOR FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
!
!     THETPM,VOLPM=air-filled porosity,volume
!     VP1,VPY=soil,litter-soil equilibrium vapor concentration
!     TK1X=soil temperature
!     PSISV1=soil matric+osmotic water potentials
!     FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
!     FLVR1X=litter-soil vapor flux
!     HFLVR1X=convective heat of litter-soil vapor flux
!     TKXR,TK1X=interim calculation of litter,soil temperatures
!
    !both litter layer and topsoil are not-saturated
    IF(VOLPM(M,0,NY,NX).GT.ZEROS(NY,NX).AND.VOLPM(M,NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      VP1=vapsat(TK1X)*EXP(18.0_r8*PSISV1/(RGAS*TK1X))
      FLVC=ATCNVS*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ

      if(abs(FLVC)>1.e20_r8)then
        write(*,*)'ATCNVS=',ATCNVS,VPR,VP1
        write(*,*)'FSNW(NY,NX)*CVRD(NY,NX)=',FSNW(NY,NX),CVRD(NY,NX)
        write(*,*)'at line',__LINE__
        call endrun(trim(mod_filename)//'at line',__LINE__)
      endif
      VPY=(VPR*VOLPM(M,0,NY,NX)+VP1*VOLPM(M,NUM(NY,NX),NY,NX)) &
        /(VOLPM(M,0,NY,NX)+VOLPM(M,NUM(NY,NX),NY,NX))
      FLVX=(VPR-VPY)*VOLPM(M,0,NY,NX)*XNPC
      IF(FLVC.GE.0.0_r8)THEN
        FLVR1X=AZMAX1(AMIN1(FLVC,FLVX,VOLW0M(L,NY,NX)*XNPB))
        if(abs(FLVR1X)>1.0e20_r8)then
          write(*,*)'FLVC,FLVX,VOLW0M(L,NY,NX)*XNPB=',FLVC,FLVX,VOLW0M(L,NY,NX)*XNPB
          write(*,*)'at line',__LINE__
          call endrun(trim(mod_filename)//'at line',__LINE__)
        endif
        HFLVR1X=(cpw*TKXR+VAP)*FLVR1X  !enthalpy flux
      ELSE
        FLVR1X=AZMIN1(AMAX1(FLVC,FLVX))
        if(abs(FLVR1X)>1.0e20_r8)then
          write(*,*)'FLVC,FLVX=',FLVC,FLVX
          write(*,*)'at line',__LINE__
          call endrun(trim(mod_filename)//'at line',__LINE__)
        endif
        HFLVR1X=(cpw*TK1X+VAP)*FLVR1X
      ENDIF
    ELSE
      FLVR1X=0.0_r8
      HFLVR1X=0.0_r8
    ENDIF
    !update litter layer temperature
    TKXR=TKXR-HFLVR1X/VHCP1(0,NY,NX)   
    !update top soil layer temperature
    TK1X=TK1X+HFLVR1X/VHCP1(NUM(NY,NX),NY,NX)
!
!     HEAT FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
!
!     TKY=litter-soil equilibrium temperature
!     HFLWC,HFLWX=litter-soil heat flux unltd,ltg by heat
!     HFLWR1X=litter-soil heat flux
!
    TKY=(TKXR*VHCP1(0,NY,NX)+TK1X*VHCP1(NUM(NY,NX),NY,NX))/(VHCP1(0,NY,NX)+VHCP1(NUM(NY,NX),NY,NX))
    HFLWX=(TKXR-TKY)*VHCP1(0,NY,NX)*XNPC
    HFLWC=ATCNDS*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
    IF(HFLWC.GE.0.0_r8)THEN
      HFLWR1X=AZMAX1(AMIN1(HFLWX,HFLWC))
    ELSE
      HFLWR1X=AZMIN1(AMAX1(HFLWX,HFLWC))
    ENDIF
!
!     ACCUMULATE SNOW-LITTER, LITTER-SOIL HEAT FLUXES
!     WITHIN LONGER TIME STEP FOR SNOWPACK FLUX CALCULATIONS
!
    FLVSR=FLVSR+FLVSRX
    HFLVSR=HFLVSR+HFLVSRX
    HFLWSR=HFLWSR+HFLWSRX
    FLVR1=FLVR1+FLVR1X
    if(abs(FLVR1)>1.0e20_r8)then
      write(*,*)'FLVR1X=',FLVR1X
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HFLVR1=HFLVR1+HFLVR1X
    HFLWR1=HFLWR1+HFLWR1X
    TK0X=TK0X-HFLVSRX/VHCPWMM(L,NY,NX)
    TKXR=TKXR+(HFLVSRX-HFLWR1X)/VHCP1(0,NY,NX)
    TK1X=TK1X+HFLWR1X/VHCP1(NUM(NY,NX),NY,NX)

  ENDDO D4000
  
  end subroutine SnowSurLitterIterate

!------------------------------------------------------------------------------------------
  subroutine CopySnowStates(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
! VOLS0,VOLSSL=snowpack snow content (water equivalent)
! VOLI0,VOLISSL=snowpack ice content
! VOLW0,VOLWSL=snowpack water content
! VOLS1,VOLSL=snowpack volume
! DLYRS0,DLYRS=snowpack depth
! VHCPWM,VHCPW=snowpack heat capacity
! TK0,TKW=snowpack temperature
!
  D60: DO L=1,JS
    VOLS0(L,NY,NX)=VOLSSL(L,NY,NX)
    VOLI0(L,NY,NX)=VOLISL(L,NY,NX)
    VOLW0(L,NY,NX)=VOLWSL(L,NY,NX)
    VOLS1(L,NY,NX)=VOLSL(L,NY,NX)
    DLYRS0(L,NY,NX)=DLYRS(L,NY,NX)
    VHCPWM(1,L,NY,NX)=VHCPW(L,NY,NX)
    TK0(L,NY,NX)=TKW(L,NY,NX)
  ENDDO D60  
  end subroutine CopySnowStates
!------------------------------------------------------------------------------------------

  subroutine UpdateSoilWaterPotential(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: THETW1,VPY,FCX,WPX
  real(r8) :: FCLX,WPLX,PSDX,FCDX

  THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX),AMIN1(POROS(NUM(NY,NX),NY,NX) &
    ,safe_adb(VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX))))
  IF(BKVL(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
        +((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
        /FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
    ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX) &
        +(((PSL(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
        /PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*PSISD(NY,NX)))
    ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
    ENDIF
  ELSEIF(VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    FCX=FCI*THETIX(NUM(NY,NX),NY,NX)
    WPX=WPI*THETIX(NUM(NY,NX),NY,NX)
    FCLX=LOG(FCX)
    WPLX=LOG(WPX)
    PSDX=PSL(NUM(NY,NX),NY,NX)-FCLX
    FCDX=FCLX-WPLX
    IF(THETWX(NUM(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
        +((FCLX-LOG(THETWX(NUM(NY,NX),NY,NX)))/FCDX*PSIMD(NY,NX))))
    ELSEIF(THETWX(NUM(NY,NX),NY,NX).LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX) &
        +(((PSL(NUM(NY,NX),NY,NX)-LOG(THETWX(NUM(NY,NX),NY,NX)))/PSDX)*PSISD(NY,NX)))
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
  subroutine SnowSurLitterExch(M,L,NY,NX,CNV1,CNV2,TCND1W,VOLP01,PSISV1,TCNDS,&
    FLVSR,HFLVSR,HFLWSR,FLVR1,HFLVR1,HFLWR1 )
  implicit none
  integer,  intent(in) :: M,L,NY,NX
  real(r8) , intent(in):: CNV1,CNV2,TCND1W,VOLP01,PSISV1,TCNDS
  real(r8) ,intent(inout) :: FLVSR,HFLVSR,HFLWSR,FLVR1,HFLVR1,HFLWR1
  real(r8) :: TK0X,TKXR,TK1X
  real(r8) :: CNVR,ATCNVR,ATCNVS
  real(r8) :: THETRR,TCNDR,ATCNDR,ATCNDS
  real(r8) :: THETWR

  ! THETWR,THETW1=litter, soil water concentration
  ! VOLWRX=litter water retention capacity
  ! PSISM1(0,PSISM1(NUM=litter,soil water potentials

  IF(VOLR(NY,NX).GT.ZEROS(NY,NX).AND.VOLW1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
    IF(THETWR.LT.FC(0,NY,NX))THEN
      PSISM1(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)+((FCL(0,NY,NX)-LOG(THETWR))/FCD(0,NY,NX)*PSIMD(NY,NX))))
    ELSEIF(THETWR.LT.POROS0(NY,NX))THEN
      PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX)+(((PSL(0,NY,NX)-LOG(THETWR))/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
    ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
    ENDIF
  ELSE
    THETWR=POROS0(NY,NX)
    PSISM1(0,NY,NX)=PSISE(0,NY,NX)
  ENDIF


  TK0X=TempK_snowpack(L,NY,NX)
  TKXR=TK1(0,NY,NX)
  TK1X=TK1(NUM(NY,NX),NY,NX)
  CNVR=WGSGR(NY,NX)*THETPM(M,0,NY,NX)*POROQ*THETPM(M,0,NY,NX)/POROS(0,NY,NX)

  IF(CVRD(NY,NX).GT.ZERO)THEN
    IF(CNV1.GT.ZERO.AND.CNVR.GT.ZERO)THEN
      ATCNVR=2.0_r8*CNVR*CNV1/(CNV1*DLYRR(NY,NX)+CNVR*DLYRS0(L,NY,NX))
    ELSE
      ATCNVR=2.0_r8*CNV1/(DLYRR(NY,NX)+DLYRS0(L,NY,NX))
    ENDIF

    IF(CNVR.GT.ZERO.AND.CNV2.GT.ZERO)THEN
      ATCNVS=2.0_r8*CNVR*CNV2/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRR(NY,NX))
    ELSE
      ATCNVS=2.0_r8*CNV2/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))
    ENDIF
    THETRR=AZMAX1(1.0_r8-THETPX(0,NY,NX)-THETWX(0,NY,NX)-THETIX(0,NY,NX))
    TCNDR=(0.779_r8*THETRR*9.050E-04_r8+0.622_r8*THETWX(0,NY,NX) &
      *2.067E-03_r8+0.380_r8*THETIX(0,NY,NX)*7.844E-03_r8+THETPX(0,NY,NX) &
      *9.050E-05_r8)/(0.779_r8*THETRR+0.622_r8*THETWX(0,NY,NX) &
      +0.380_r8*THETIX(0,NY,NX)+THETPX(0,NY,NX))

    IF(TCND1W.GT.ZERO.AND.TCNDR.GT.ZERO)THEN
      ATCNDR=2.0_r8*TCND1W*TCNDR/(TCND1W*DLYRR(NY,NX)+TCNDR*DLYRS0(L,NY,NX))
    ELSE
      ATCNDR=0.0_r8
    ENDIF

    IF(TCNDR.GT.ZERO.AND.TCNDS.GT.ZERO)THEN
      ATCNDS=2.0_r8*TCNDR*TCNDS/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX)+TCNDS*DLYRR(NY,NX))
    ELSE
      ATCNDS=0.0_r8
    ENDIF
  ELSE
    ATCNVR=0.0_r8
    ATCNVS=0.0_r8
    ATCNDR=0.0_r8
    ATCNDS=0.0_r8
  ENDIF
  !
  ! SHORTER TIME STEP FOR SURFACE RESIDUE FLUX CALCULATIONS
  call SnowSurLitterIterate(L,M,NY,NX,ATCNDR,ATCNDS,ATCNVS,ATCNVR,PSISV1,&
    VOLP01,TK0X,TKXR,TK1X,FLVSR,HFLVR1,HFLVSR,HFLWSR,HFLWR1)

  end subroutine SnowSurLitterExch

!------------------------------------------------------------------------------------------

  subroutine SnowTopSoilExch(M,L,NY,NX,CNV1,VP1,VOLP01,TCND1W,FLVS1,HFLVS1,HFLWS1,CNV2,PSISV1,TCNDS)
  implicit none
  integer , intent(in) :: M,L,NY,NX
  real(r8), intent(in) :: CNV1,VP1,VOLP01,TCND1W
  real(r8), intent(out):: FLVS1,HFLVS1,HFLWS1,CNV2,PSISV1,TCNDS
  real(r8) :: ATCNDS
  real(r8) :: WTHET2,VPY,VP2,ATCNVS,FLVC,FLVX
  real(r8) :: HFLWC,HFLWX,TKWX1,TKY

  call UpdateSoilWaterPotential(NY,NX)

  PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
!
  ! VAPOR FLUX BETWEEN SNOWPACK AND SOIL SURFACE
  !
  ! VOLP01,THETPM=air volume,concentration
  ! CNV1,CNV2=vapor conductances of source, destination layers
  ! VP1,VP2=vapor concentrations of source, destination layers
  ! POROS,POROQ=porosity, tortuosity
  ! WGSGL=vapor diffusivity
  ! TempK_snowpack,TK1=snow,soil surface temperature
  ! PSISV1=soil matric+osmotic potential
  ! ATCNVS=snow-soil vapor conductance
  ! DLYR=soil surface layer depth
  ! FLVC,FLVX=vapor flux unlimited,limited by vapor
  ! VPY=equilibrium vapor concentration
  ! XNPX=time step for flux calculations
  ! FLVS1,HFLVS1=vapor flux and its convective heat flux
  !
  IF(VOLP01.GT.ZEROS2(NY,NX).AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
    CNV2=WGSGL(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ &
      *THETPM(M,NUM(NY,NX),NY,NX)/POROS(NUM(NY,NX),NY,NX)
    VP2=vapsat(TK1(NUM(NY,NX),NY,NX))*EXP(18.0_r8*PSISV1/(RGAS*TK1(NUM(NY,NX),NY,NX)))
    ATCNVS=2.0_r8*CNV1*CNV2/(CNV1*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRS0(L,NY,NX))
    FLVC=ATCNVS*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*BARE(NY,NX)*XNPY
    VPY=(VP1*VOLP01+VP2*VOLPM(M,NUM(NY,NX),NY,NX))/(VOLP01+VOLPM(M,NUM(NY,NX),NY,NX))
    FLVX=(VP1-VPY)*VOLP01*XNPA

    IF(FLVC.GE.0.0_r8)THEN
      !water flux goes into soil
      FLVS1=AZMAX1(AMIN1(FLVC,FLVX))
      !heat flux
      HFLVS1=(cpw*TempK_snowpack(L,NY,NX)+VAP)*FLVS1
    ELSE
      !water flux out of soil
      FLVS1=AZMIN1(AMAX1(FLVC,FLVX))
      HFLVS1=(cpw*TK1(NUM(NY,NX),NY,NX)+VAP)*FLVS1
    ENDIF
  ELSE
    CNV2=0.0_r8
    FLVS1=0.0_r8
    HFLVS1=0.0_r8
  ENDIF

  !
  ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE SOIL
  !
  ! WTHET2=multiplier for air concentration in thermal conductivity
  ! TCND1W,TCNDS=thermal conductivity of snowpack, soil surface
  ! STC,DTC=mineral component of thermal conductivity
  ! THETWX,THETIX,THETPX=soil surface water,ice,air concentrations
  ! BARE=soil surface fraction
  ! ATCNDS=snowpack-soil thermal conductance
  ! TKWX1=interim snowpack temperature
  ! TKY=equilibrium temperature
  ! HFLWX,HFLWC=heat-constrained,heat-unconstrained heat fluxes
  ! XNPY=time step for snowpack flux calculations
  ! HFLWS1=snowpack-soil heat flux
  ! THETPY=air-filled pores (including macropores and micropores)
  WTHET2=1.467_r8-0.467_r8*THETPY(NUM(NY,NX),NY,NX)
  TCNDS=(STC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX) &
    *2.067E-03_r8+0.611_r8*THETIX(NUM(NY,NX),NY,NX)*7.844E-03_r8 &
    +WTHET2*THETPX(NUM(NY,NX),NY,NX)*9.050E-05_r8) &
    /(DTC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX) &
    +0.611_r8*THETIX(NUM(NY,NX),NY,NX) &
    +WTHET2*THETPX(NUM(NY,NX),NY,NX))

  !the mean thermal conductivity
  IF(BARE(NY,NX).GT.ZERO)THEN
    ATCNDS=2.0_r8*TCND1W*TCNDS/(TCND1W*DLYR(3,NUM(NY,NX),NY,NX)+TCNDS*DLYRS0(L,NY,NX))
  ELSE
    ATCNDS=0.0_r8
  ENDIF

  TKWX1=TK1(NUM(NY,NX),NY,NX)+HFLVS1/VHCP1(NUM(NY,NX),NY,NX)
  TKY=(TempK_snowpack(L,NY,NX)*VHCPWMM(L,NY,NX)+TKWX1*VHCP1(NUM(NY,NX),NY,NX))/(VHCPWMM(L,NY,NX)+VHCP1(NUM(NY,NX),NY,NX))
  HFLWX=(TempK_snowpack(L,NY,NX)-TKY)*VHCPWMM(L,NY,NX)*XNPA
  HFLWC=ATCNDS*(TempK_snowpack(L,NY,NX)-TKWX1)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*BARE(NY,NX)*XNPY
  IF(HFLWC.GE.0.0_r8)THEN
    HFLWS1=AZMAX1(AMIN1(HFLWX,HFLWC))
  ELSE
    HFLWS1=AZMIN1(AMAX1(HFLWX,HFLWC))
  ENDIF
  end subroutine SnowTopSoilExch
end module SnowPhysMod
