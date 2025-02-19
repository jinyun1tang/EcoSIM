module SurfLitterPhysMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only : etimer  
  use DebugToolMod
  use MiniMathMod
  use HydroThermData
  use SurfSoilDataType
  use GridDataType
  use SurfLitterDataType
  use ChemTranspDataType
  use SoilWaterDataType
  use SoilPropertyDataType
  use SurfPhysData
  use SoilPhysDataType
  use EcoSIMCtrlDataType
  use SOMDataType
  use LandSurfDataType
  use EcoSimConst
  use EcoSIMSolverPar
  use abortutils
  use PhysPars
  use UnitMod, only : units  
  use SoilPhysParaMod
  USE SoilHeatDataType
  USE ClimForcDataType
  use SnowPhysData
implicit none
  private
  character(len=*), parameter, private :: mod_filename=&
  __FILE__
  real(r8), parameter :: tiny_wat=1.e-13_r8

  public :: SurfLitREnergyBalanceM
  public :: UpdateLitRPhys
  public :: LateralGridsHdryoExch
  public :: UpdateLitRBe4RunoffM
  public :: UpdateLitRAftRunoff
  contains


!------------------------------------------------------------------------------------------
  subroutine SurfLitREnergyBalanceM(I,J,M,NY,NX,PSISV1,Prec2LitR2,RainHeat2LitR2,&
    CumNetWatFlow2LitR,CumNetHeatFlow2LitR,CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil,&
    HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,EvapLitR2Soi1,VapXAir2LitR,HeatFluxAir2LitR)
!!
! Description:
  implicit none
  integer , intent(in) :: I,J
  integer , intent(in) :: M       !soil heat-water iteration id
  integer , intent(in) :: NY,NX
  real(r8), intent(in) :: PSISV1
  real(r8), intent(in) :: Prec2LitR2       !precipitation flux to litter
  real(r8), intent(in) :: RainHeat2LitR2   !heat flux to litter by precipitation
  real(r8), intent(in) :: CumNetWatFlow2LitR
  real(r8), intent(in) :: CumNetHeatFlow2LitR
  real(r8), intent(in) :: CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil
  real(r8), intent(out) :: HeatSensLitR2Soi1
  real(r8), intent(out) :: HeatSensVapLitR2Soi1
  real(r8), intent(out) :: EvapLitR2Soi1
  real(r8), intent(out) :: VapXAir2LitR         !water vapor flux from canopy air to litr
  real(r8), intent(out) :: HeatFluxAir2LitR     !Heat flux from canopy air to litr
  
  character(len=*), parameter :: subname='SurfLitREnergyBalanceM'
  real(r8) :: HeatSensAir2LitR
  real(r8) :: HeatSensEvapAir2LitR
  real(r8) :: LWRadLitR  
  real(r8) :: VWatLitr2,CdTLit2Soil,CdVaporLit2Soil,Radt2LitR
  real(r8) :: CNVR,CNV1
  real(r8) :: TCNDR
  real(r8) :: TCND1
  real(r8) :: DTKX
  real(r8) :: LitRAlbedo  
  real(r8) :: Radnet2LitR         !net radiation on to litter [MJ]
  real(r8) :: RadSWByLitR         !shortwave radiation absorbed by litter [MJ]
  real(r8) :: LatentHeatAir2LitR


! begin_execution
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!
  call PrintInfo('beg '//subname)

  VapXAir2LitR         = 0.0_r8
  Radnet2LitR          = 0.0_r8
  LatentHeatAir2LitR   = 0.0_r8
  HeatSensEvapAir2LitR = 0.0_r8
  HeatSensAir2LitR     = 0.0_r8
  HeatFluxAir2LitR     = 0.0_r8

  EvapLitR2Soi1        = 0.0_r8
  HeatSensVapLitR2Soi1 = 0.0_r8
  HeatSensLitR2Soi1    = 0.0_r8
  LWRadLitR            = 0.0_r8

  IF(VHeatCapacity1_vr(0,NY,NX).LE.VHeatCapLitRMin_col(NY,NX))THEN
    TKSoil1_vr(0,NY,NX)=TKSoil1_vr(NUM(NY,NX),NY,NX)
  ENDIF  
!
! NET RADIATION AT RESIDUE SURFACE
!
! LitRAlbedo=litter albedo
! VLWatMicP1,VLiceMicP1=water,ice volume in litter
! RadSW2LitR_col,LWRad2LitR_col=incoming shortwave,longwave radiation

!
  !
  ! THERMAL CONDUCTIVITY BETWEEN SURFACE RESIDUE AND SOIL SURFACE
  !
  
  IF(FracSurfByLitR_col(NY,NX).GT.ZERO)THEN
  ! litter layer
  ! albedo
    LitRAlbedo=(0.20_r8*VLSoilMicPMass_vr(0,NY,NX)+0.06_r8*VLWatMicP1_vr(0,NY,NX)+0.30_r8 &
      *VLiceMicP1_vr(0,NY,NX))/(VLSoilMicPMass_vr(0,NY,NX)+VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX))

    !radiation incident on litter layer  
    RadSWByLitR          = (1.0_r8-LitRAlbedo)*RadSW2LitR_col(NY,NX)
    Radt2LitR            = RadSWByLitR+LWRad2LitR_col(NY,NX)
    Eco_RadSW_col(NY,NX) = Eco_RadSW_col(NY,NX) + RadSWByLitR

    CNVR = VaporDiffusivityLitR_col(NY,NX)*AirFilledSoilPoreM_vr(M,0,NY,NX)*POROQ*AirFilledSoilPoreM_vr(M,0,NY,NX)/POROS_vr(0,NY,NX)
    CNV1 = WVapDifusvitySoil_vr(NUM(NY,NX),NY,NX)*AirFilledSoilPoreM_vr(M,NUM(NY,NX),NY,NX)*POROQ &
      *AirFilledSoilPoreM_vr(M,NUM(NY,NX),NY,NX)/POROS_vr(NUM(NY,NX),NY,NX)
    !there is litter layer
    IF(CNVR.GT.ZERO.AND.CNV1.GT.ZERO)THEN
      CdVaporLit2Soil=2.0_r8*CNVR*CNV1/(CNVR*DLYR_3D(3,NUM(NY,NX),NY,NX)+CNV1*DLYRR_COL(NY,NX))
    ELSE
      !below is a numerical hack
      CdVaporLit2Soil=2.0_r8*CNVR/(DLYR_3D(3,NUM(NY,NX),NY,NX)+DLYRR_COL(NY,NX))*FracSurfByLitR_col(NY,NX)
    ENDIF

    DTKX=ABS(TKSoil1_vr(0,NY,NX)-TKSoil1_vr(NUM(NY,NX),NY,NX))*ppmc

    call CalcLitRThermConductivity(NY,NX,DTKX,TCNDR)

    call CalcSoilThermConductivity(NX,NY,NUM(NY,NX),DTKX,TCND1)

    CdTLit2Soil=2.0_r8*TCNDR*TCND1/(TCNDR*DLYR_3D(3,NUM(NY,NX),NY,NX)+TCND1*DLYRR_COL(NY,NX))
  !
  ! SMALLER TIME STEP FOR SOLVING SURFACE RESIDUE ENERGY EXCHANGE
  !
    call SurfLitterIterationM(I,J,M,NY,NX,CdTLit2Soil,CdVaporLit2Soil,PSISV1,Radt2LitR,Prec2LitR2,RainHeat2LitR2,&
      CumNetWatFlow2LitR,CumNetHeatFlow2LitR,CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil,&
      EvapLitR2Soi1,HeatSensAir2LitR,HeatSensEvapAir2LitR,HeatSensLitR2Soi1,&
      HeatSensVapLitR2Soi1,LatentHeatAir2LitR,LWRadLitR,Radnet2LitR,VapXAir2LitR,HeatFluxAir2LitR)

    HeatByRad2Surf_col(NY,NX)      = HeatByRad2Surf_col(NY,NX)+Radnet2LitR
    HeatNet2Surf_col(NY,NX)        = HeatNet2Surf_col(NY,NX)+Radnet2LitR+HeatSensEvapAir2LitR+HeatSensAir2LitR
    HeatSensAir2Surf_col(NY,NX)    = HeatSensAir2Surf_col(NY,NX)+HeatSensAir2LitR
    HeatEvapAir2Surf_col(NY,NX)    = HeatEvapAir2Surf_col(NY,NX)+LatentHeatAir2LitR
    HeatSensVapAir2Surf_col(NY,NX) = HeatSensVapAir2Surf_col(NY,NX)+HeatSensEvapAir2LitR
    LWRadBySurf_col(NY,NX)         = LWRadBySurf_col(NY,NX)+LWRadLitR
  ELSE
    CdVaporLit2Soil=0.0_r8
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine SurfLitREnergyBalanceM
!------------------------------------------------------------------------------------------

  subroutine CalcLitRThermConductivity(NY,NX,DTKX,TCNDR)

  use PhysPars
  implicit none
  integer , intent(in) :: NY,NX
  real(r8), intent(in) :: DTKX    !absolute temperautre gradient between litter and soil
  real(r8), intent(out) :: TCNDR
  real(r8) :: HeatDiffusByWat0,HeatDiffusByAir0
  real(r8) :: RYLXW0,RYLXA0,RYLNA0,RYLNW0
  real(r8) :: XNUSW0,XNUSA0
  real(r8) :: TCNDW0,TCNDA0,WTHET0,THETRR

  !THETRR=litter in relative volume
  !FracSoilAsAirt=relative volume as air  
  THETRR           = AZMAX1(1.0_r8-AirFilledSoilPore_vr(0,NY,NX)-FracSoiPAsWat_vr(0,NY,NX)-FracSoiPAsIce_vr(0,NY,NX))
  HeatDiffusByWat0 = AZMAX1(FracSoiPAsWat_vr(0,NY,NX)-TRBW)**3._r8
  HeatDiffusByAir0 = AZMAX1(AirFilledSoilPore_vr(0,NY,NX)-TRBA)**3._r8
  RYLXW0           = DTKX*HeatDiffusByWat0
  RYLXA0           = DTKX*HeatDiffusByAir0
  RYLNW0           = AMIN1(1.0E+04_r8,RYLXW*RYLXW0)
  RYLNA0           = AMIN1(1.0E+04_r8,RYLXA*RYLXA0)
  XNUSW0           = AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW0**0.25_r8/DNUSW)
  XNUSA0           = AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA0**0.25_r8/DNUSA)
  TCNDW0           = 2.067E-03_r8*XNUSW0
  TCNDA0           = 9.050E-05_r8*XNUSA0
  WTHET0           = 1.467_r8-0.467_r8*FracSoilAsAirt(0,NY,NX)
  TCNDR            = (0.779_r8*THETRR*9.050E-04_r8+0.622_r8*FracSoiPAsWat_vr(0,NY,NX)*TCNDW0 &
    +0.380_r8*FracSoiPAsIce_vr(0,NY,NX)*7.844E-03_r8+WTHET0*AirFilledSoilPore_vr(0,NY,NX)*TCNDA0) &
    /(0.779_r8*THETRR+0.622_r8*FracSoiPAsWat_vr(0,NY,NX)+0.380_r8*FracSoiPAsIce_vr(0,NY,NX)+WTHET0*AirFilledSoilPore_vr(0,NY,NX))
  end subroutine CalcLitRThermConductivity

!------------------------------------------------------------------------------------------

  subroutine SurfLitterIterationM(I,J,M,NY,NX,CdTLit2Soil,CdVaporLit2Soil,PSISV1,Radt2LitR,Prec2LitR2,RainHeat2LitR2,&
    CumNetWatFlow2LitR,CumNetHeatFlow2LitR,CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumNetHeatFlow2Soil,&
    EvapLitR2Soi1,HeatSensAir2LitR,HeatSensEvapAir2LitR,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,LatentHeatAir2LitR,LWRadLitR,Radnet2LitR,VapXAir2LitR,HeatFluxAir2LitR)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M     !soil heat-flow iteration id
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: CdTLit2Soil     !heat conductance between litter and soil [1/(m2 K h)]
  real(r8), intent(in) :: CdVaporLit2Soil     !water conductance between litter and soil [m/h]
  real(r8), intent(in) :: PSISV1     !surface soil suction pressure [MPa]
  real(r8), intent(in) :: Radt2LitR    !incoming radiation (from canopy air and sun) onto litter [MJ]
  real(r8), intent(in) :: Prec2LitR2                !precipitation water flux added to litter [m3H2O d-2]
  real(r8), intent(in) :: RainHeat2LitR2           !precipitation heat flux added to litter [MJ]
  real(r8), intent(in) :: CumNetWatFlow2LitR      !cumulative net water flux to litter
  real(r8), intent(in) :: CumNetHeatFlow2LitR     !cumulative net heat flux to litter
  real(r8), intent(in) :: CumWatFlx2SoiMicP
  real(r8), intent(in) :: CumWatFlx2SoiMacP
  real(r8), intent(in) :: cumNetHeatFlow2Soil
  real(r8), intent(inout) :: EvapLitR2Soi1         !cumulated evaporation flux from litter to soil [kg H2O]
  real(r8), intent(inout) :: HeatSensAir2LitR      !cumulated sensible heat flux from air to litter [MJ]
  real(r8), intent(inout) :: HeatSensEvapAir2LitR  !cumulated sensible heat associated with water flux [MJ]
  real(r8), intent(inout) :: HeatSensLitR2Soi1     !cumulated sensible heat flux from litter to soil [MJ]
  real(r8), intent(inout) :: HeatSensVapLitR2Soi1  !cumulated heat flux from litter to soil due to water vapor flux [MJ]
  real(r8), intent(inout) :: LatentHeatAir2LitR    !cumulated latent heat from air to litter [MJ]
  real(r8), intent(inout) :: LWRadLitR             !cumulated outgoing long wave radiation from litter [MJ]
  real(r8), intent(inout) :: Radnet2LitR           !cumulated net radiation on litter [MJ]
  real(r8), intent(inout) :: VapXAir2LitR
  real(r8), intent(inout) :: HeatFluxAir2LitR       !cumulated heat into litter [MJ]

  character(len=*), parameter :: subname='SurfLitterIterationM'
  integer  :: NN
  real(r8) :: tk1pre,VHCPRXX
  real(r8) :: RAa,VaporSoi1,VapLitR,VPY
  real(r8) :: TKXR,TK1X,TKY
  real(r8) :: FLVC,FLVX,ENGYR
  real(r8) :: HFLWC,HFLWX
  real(r8) :: EVAPR2                 !(>0)evporation flux from air to litter [ton H2O]
  real(r8) :: LatentHeatAir2LitR2    !latent heat flux from air to litter [MJ]
  real(r8) :: HeatSensEvapAir2LitR2  !sensible heat associated with water flux [MJ]
  real(r8) :: NetHeatAir2LitR2       !residual heat flux into litter by subtracting latent and sensible heat from radiation [MJ]
  real(r8) :: LWRadLitR2             !outgoing long wave radiation from litter [MJ]
  real(r8) :: ThetaWLitR             !volumetric water held in litter [m3 water/m3 litter]
  real(r8) :: HeatSensAir2LitR2      !sensible heat flux from air to litter [MJ]
  real(r8) :: EvapLitR2Soi2          !evaporation flux from litter to soil [kg H2O]
  real(r8) :: HeatSensVapLitR2Soi2   !advective heat by evaporation [MJ]
  real(r8) :: CdLitREvap             !conductance of evaporation [m3]
  real(r8) :: CdLitRHSens            !conductance of sensible heat [MJ/K]
  real(r8) :: RAGX,RI
  real(r8) :: Radnet2LitR2           !litter net radiation,after taking out outgoing radiation from litter [MJ]
  real(r8) :: HeatFluxAir2LitR2      !total residual heat into litter [MJ]
  real(r8) :: HeatSensLitR2Soi2      !sensible heat flux from litter to soil
  real(r8) :: VLHeatCapcityLitR2,VLHeatCapacity2
  real(r8) :: VLWatMicP12,VWatLitr2
  real(r8) :: TKR1,TKS1,dWLitR
  real(r8) :: ENGYS
  real(r8) :: dt_litrHeat            !time step for solving litter heat/water fluxes
  real(r8) :: dLWdTSoil,dLWSoil
  real(r8) :: tHeatLitR2Soil2

! begin_execution
  call PrintInfo('beg '//subname)

  VWatLitr2          = VLWatMicP1_vr(0,NY,NX)+CumNetWatFlow2LitR
  VLHeatCapcityLitR2 = VHeatCapacity1_vr(0,NY,NX)+cpw*CumNetWatFlow2LitR
  TKR1               = (TKSoil1_vr(0,NY,NX)*VHeatCapacity1_vr(0,NY,NX)+CumNetHeatFlow2LitR)/VLHeatCapcityLitR2

  VLHeatCapacity2    = VHeatCapacity1_vr(NUM(NY,NX),NY,NX)
  VLWatMicP12        = VLWatMicP1_vr(NUM(NY,NX),NY,NX)
  TKS1               = TKSoil1_vr(NUM(NY,NX),NY,NX)
  dLWdTSoil          = 0._r8

  !embedded iteration, local time step size
  dt_litrHeat=dts_HeatWatTP/real(NPR,kind=r8)      !time step for litter flux calculation

  D5000: DO NN=1,NPR
!    write(*,*)'HeatFluxAir2LitR  x',NN,HeatFluxAir2LitR
    IF(VLHeatCapcityLitR2.GT.VHeatCapLitRMin_col(NY,NX))THEN
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
      ! CdLitREvap,CdLitRHSens=conductance for litter latent,sensible heat fluxes
      !
      RI   = RichardsonNumber(RIB_col(NY,NX),TKQ_col(NY,NX),TKR1)
      RAGX = AMAX1(RAM,0.8_r8*ResistAreodynOverLitr_col(NY,NX),&
        AMIN1(1.2_r8*ResistAreodynOverLitr_col(NY,NX),RARG(NY,NX)/(1.0_r8-10.0_r8*RI)))
      ResistAreodynOverLitr_col(NY,NX) = RAGX
      RAa                              = RAGX
      CdLitREvap                       = AScaledCdWOverLitr_col(NY,NX)/(RAa+RZ)
      CdLitRHSens                      = AScaledCdHOverLitr_col(NY,NX)/RAa
!
!     NET RADIATION AT RESIDUE SURFACE
!
!     LWRadLitR2=longwave radiation emitted by litter
!     ThetaWLitR=litter water content
!     VWatLitRHoldCapcity=maximum water retention by litter
!     PSISM1=litter matric water potential
!
      LWRadLitR2   = LWEmscefLitR_col(NY,NX)*TKR1**4._r8/real(NPR,kind=r8)
      Radnet2LitR2 = Radt2LitR-LWRadLitR2

      IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS2(NY,NX))THEN
        ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VWatLitr2)/VLitR_col(NY,NX)
      ELSE
        ThetaWLitR=POROS0_col(NY,NX)
      ENDIF

      IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX) .AND. VLWatMicP1_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
        ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VWatLitr2)/VLitR_col(NY,NX)
        IF(ThetaWLitR.LT.FieldCapacity_vr(0,NY,NX))THEN
          PSISM1_vr(0,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD(NY,NX) &
            +((LOGFldCapacity_vr(0,NY,NX)-LOG(ThetaWLitR))/FCD_vr(0,NY,NX)*LOGPSIMND(NY,NX))))
        ELSEIF(ThetaWLitR.LT.POROS0_col(NY,NX))THEN
          PSISM1_vr(0,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(0,NY,NX)-LOG(ThetaWLitR)) &
            /PSD_vr(0,NY,NX))**SRP_vr(0,NY,NX)*LOGPSIMXD(NY,NX)))
        ELSE
          ThetaWLitR         = POROS0_col(NY,NX)
          PSISM1_vr(0,NY,NX) = PSISE_vr(0,NY,NX)
        ENDIF
      ELSE
        ThetaWLitR         = POROS0_col(NY,NX)
        PSISM1_vr(0,NY,NX) = PSISE_vr(0,NY,NX)
      ENDIF
!
!     VAPOR FLUX AT RESIDUE SURFACE
!
!     VapLitR,VaporSoi1,VPQ_col=vapor pressure in litter,soil,canopy air, m3/m-3
!     TKS1=soil temperature
!     EVAPR2=negative of litter evaporation,<0 into atmosphere
!     LatentHeatAir2LitR2=litter latent heat flux
!     VAP=latent heat of evaporation
!     HeatSensEvapAir2LitR2=convective heat of evaporation flux
!
!     in litter      
      VapLitR               = vapsat(TKR1)*EXP(18.0_r8*PSISM1_vr(0,NY,NX)/(RGASC*TKR1))
      VaporSoi1             = vapsat(TKS1)*EXP(18.0_r8*PSISV1/(RGASC*TKS1))    !vapor pressure in soil, ton H2O/m3
      EVAPR2                = AMAX1(-AZMAX1(VWatLitr2)*dts_wat,CdLitREvap*(VPQ_col(NY,NX)-VapLitR)) ![ton H2O]
      LatentHeatAir2LitR2   = EVAPR2*EvapLHTC          !latent energy flux,                             MJ/ton H2O * ton H2O = MJ
      HeatSensEvapAir2LitR2 = EVAPR2*cpw*TKR1          !mass energy flux
      !
      ! SOLVE FOR RESIDUE TO SOIL SURFACE HEAT FLUXES
      !
      ! FLVC,FLVX=vapor unconstrained,vapor constrained vapor flux
      ! dt_litrHeat=time step for litter flux calculations
      ! VPY=equilibrium vapor concentration
      ! VsoiPM=litter,soil air filled porosity
      ! EvapLitR2Soi2=litter soil vapor flux
      ! HeatSensVapLitR2Soi2=convective heat of litter soil vapor flux
      ! TKXR,TK1X=interim calculation of litter,soil temperatures
      ! TKY=equilibrium litter-soil temperature
      ! HFLWC,HFLWX=litter-soil heat flux unlimited,limited by heat
      ! HeatSensLitR2Soi2=litter-soil heat flux
      ! THETPM: air-filled porosity

      ! Water flux by evaporation
      IF(AirFilledSoilPoreM_vr(M,0,NY,NX).GT.THETX .AND. AirFilledSoilPoreM_vr(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
        FLVC = CdVaporLit2Soil*(VapLitR-VaporSoi1)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dt_litrHeat
        VPY  = (VapLitR*VLsoiAirPM_vr(M,0,NY,NX)+VaporSoi1*VLsoiAirPM_vr(M,NUM(NY,NX),NY,NX))/(VLsoiAirPM_vr(M,0,NY,NX)+VLsoiAirPM_vr(M,NUM(NY,NX),NY,NX))
        FLVX = (VapLitR-VPY)*VLsoiAirPM_vr(M,0,NY,NX)*XNPB

        IF(FLVC.GE.0.0_r8)THEN
          !from litter to soil
          EvapLitR2Soi2        = AZMAX1(AMIN1(FLVC,FLVX))
          HeatSensVapLitR2Soi2 = (cpw*TKR1*HeatAdv_scal+EvapLHTC)*EvapLitR2Soi2
        ELSE
          !from soil to litter
          EvapLitR2Soi2        = AZMIN1(AMAX1(FLVC,FLVX))
          HeatSensVapLitR2Soi2 = (cpw*TKS1*HeatAdv_scal+EvapLHTC)*EvapLitR2Soi2
        ENDIF
      ELSE
        EvapLitR2Soi2        = 0.0_r8
        HeatSensVapLitR2Soi2 = 0.0_r8
      ENDIF

      TKXR  = TKR1-HeatSensVapLitR2Soi2/VLHeatCapcityLitR2                                          !update litter layer temperature
      TK1X  = TKS1+HeatSensVapLitR2Soi2/VLHeatCapacity2                                             !update soil layer temperature
      TKY   = (TKXR*VLHeatCapcityLitR2+TK1X*VLHeatCapacity2)/(VLHeatCapcityLitR2+VLHeatCapacity2)   !equilibrium temperature
      HFLWX = (TKXR-TKY)*VLHeatCapcityLitR2*XNPB                                                    !sensible heat flux > 0 into soil
      HFLWC = CdTLit2Soil*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dt_litrHeat
      IF(HFLWC.GE.0.0_r8)THEN
        !heat from litter to soil
        HeatSensLitR2Soi2=AZMAX1(AMIN1(HFLWX,HFLWC))
      ELSE
        !heat from soil to litter
        HeatSensLitR2Soi2=AZMIN1(AMAX1(HFLWX,HFLWC))
      ENDIF
      tHeatLitR2Soil2=HeatSensLitR2Soi2+HeatSensVapLitR2Soi2
!
!     SOLVE FOR RESIDUE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
!     HeatSensAir2LitR2,Radnet2LitR2,LatentHeatAir2LitR2=litter sensible,net radn,latent heat fluxes
!     NetHeatAir2LitR2,HeatFluxAir2LitR2=storage,total litter heat flux
!
      HeatSensAir2LitR2 = CdLitRHSens*(TKQ_col(NY,NX)-TKR1)    !sensible heat flux between canopy air and litter surface
      NetHeatAir2LitR2  = Radnet2LitR2+LatentHeatAir2LitR2+HeatSensAir2LitR2      !
      HeatFluxAir2LitR2 = NetHeatAir2LitR2+HeatSensEvapAir2LitR2
!
!     AGGREGATE WATER AND ENERGY FLUXES FROM TIME STEP FOR LITTER
!     CALCULATIONS TO THAT FOR SOIL PROFILE
!
      VapXAir2LitR         = VapXAir2LitR+EVAPR2
      Radnet2LitR          = Radnet2LitR+Radnet2LitR2
      LatentHeatAir2LitR   = LatentHeatAir2LitR+LatentHeatAir2LitR2
      HeatSensEvapAir2LitR = HeatSensEvapAir2LitR+HeatSensEvapAir2LitR2
      HeatSensAir2LitR     = HeatSensAir2LitR+HeatSensAir2LitR2
      HeatFluxAir2LitR     = HeatFluxAir2LitR+HeatFluxAir2LitR2

      EvapLitR2Soi1        = EvapLitR2Soi1+EvapLitR2Soi2
      HeatSensVapLitR2Soi1 = HeatSensVapLitR2Soi1+HeatSensVapLitR2Soi2
      HeatSensLitR2Soi1    = HeatSensLitR2Soi1+HeatSensLitR2Soi2
      LWRadLitR            = LWRadLitR+LWRadLitR2
    ELSE
      !not significant litter layer heat capacity
      EVAPR2                = 0.0_r8
      Radnet2LitR2          = 0.0_r8
      LatentHeatAir2LitR2   = 0.0_r8
      HeatSensEvapAir2LitR2 = 0.0_r8
      HeatSensAir2LitR2     = 0.0_r8
      HeatFluxAir2LitR2     = 0.0_r8
      EvapLitR2Soi2         = 0.0_r8
      HeatSensVapLitR2Soi2  = 0.0_r8
      HeatSensLitR2Soi2     = 0.0_r8
      LWRadLitR2            = 0.0_r8
    ENDIF
    !VWatLitr2,VLWatMicP12=water in litter, topsoil layer 
    dWLitR      = Prec2LitR2+EVAPR2-EvapLitR2Soi2
    VWatLitr2   = VWatLitr2+dWLitR
    VLWatMicP12 = VLWatMicP12+EvapLitR2Soi2
    ENGYR       = VLHeatCapcityLitR2*TKR1
    VHCPRXX     = VLHeatCapcityLitR2
    ENGYS       = VLHeatCapacity2*TKS1
    ! VLHeatCapcityLitR2: heat capacity, kJ/kg/Kelvin
    VLHeatCapcityLitR2 = VLHeatCapcityLitR2+cpw*dWLitR
    VLHeatCapacity2    = VLHeatCapacity2+cpw*EvapLitR2Soi2
    dLWSoil            = dLWdTSoil*(TKS1-TKSoil1_vr(NUM(NY,NX),NY,NX))
    tk1pre             = TKR1
    TKR1               = (ENGYR+HeatFluxAir2LitR2+RainHeat2LitR2-tHeatLitR2Soil2)/VLHeatCapcityLitR2
    TKS1               = (ENGYS+tHeatLitR2Soil2+dLWSoil)/VLHeatCapacity2

  ENDDO D5000
  call PrintInfo('end '//subname)

  end subroutine SurfLitterIterationM

!------------------------------------------------------------------------------------------
  subroutine UpdateLitRPhys(I,J,NY,NX,dWat,dHeat,HEATIN_lnd)
  !
  !Description
  !Update Litter physical variables
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: NY,NX
  real(r8), intent(in) :: dWat   !water added due to litter fall, m3 H2O/d2/h
  real(r8), intent(in) :: dHeat  !heat added due to litter fall, MJ/d2/h
  real(r8), intent(inout) :: HEATIN_lnd

  real(r8) :: VHeatCapacityLitrX  !old litr heat capacity
  real(r8) :: VHeatCapacityLitR   !current litr heat capacity
  real(r8) :: dVHeatCapacityLitr  !change in heat capacity
  real(r8) :: tkspre,ENGYR,VLWatMicPr,VLiceMicPr
  real(r8) :: ENGYZ,HeatByLitrMassChange,HS

  IF(FracSurfByLitR_col(NY,NX).LE.ZERO)return
  ! CALCULATE SURFACE RESIDUE TEMPERATURE FROM ITS CHANGE
  ! IN HEAT STORAGE
  !
  VHeatCapacityLitrX = VHeatCapacity_vr(0,NY,NX)
  VHeatCapacityLitR  = cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)+cpo*SoilOrgM_vr(ielmc,0,NY,NX)
  VLWatMicPr         = VLWatMicP_vr(0,NY,NX)
  VLiceMicPr         = VLiceMicP_vr(0,NY,NX)

  !TairK: air temperature in kelvin, HeatByLitrMassChange represents increase heat in litr
  dVHeatCapacityLitr   = VHeatCapacityLitR-VHeatCapacityLitrX
  HeatByLitrMassChange = dVHeatCapacityLitr*TairK_col(NY,NX)
  ENGYZ                = VHeatCapacityLitrX*TKS_vr(0,NY,NX)

  !update water, ice content and heat capacity of residue
!  if(etimer%get_curr_yearAD()<=1981)then  
!  write(118,*)I+J/24.,VLWatMicP_vr(0,NY,NX),WatFLo2LitR_col(NY,NX),TLitrIceFlxThaw_col(NY,NX), &
!    WatInByRunoff,PrecRainAndIrrig_col(NY,NX),FracSurfByLitR_col(NY,NX)
!  endif

  VLWatMicP_vr(0,NY,NX)     = VLWatMicP_vr(0,NY,NX)+dWat
  VHeatCapacity_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)

  THeatSoiThaw_col(NY,NX)   = THeatSoiThaw_col(NY,NX)+TLitrIceHeatFlxFrez_col(NY,NX)
  IF(VHeatCapacity_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    !when there are still significant heat capacity of the residual layer
    tkspre          = TKS_vr(0,NY,NX)
    TKS_vr(0,NY,NX) = (ENGYZ+HeatByLitrMassChange +dHeat)/VHeatCapacity_vr(0,NY,NX)

    if(TKS_vr(0,NY,NX)<100._r8 .or. TKS_vr(0,NY,NX)>360._r8)then
      write(*,*)I,J,NY,NX,TKS_vr(0,NY,NX),tkspre
      write(*,*)'WatFLo2Litr =',WatFLo2LitR_col(NY,NX)
      write(*,*)'wat flo2litr icethaw runoff',VLWatMicPr,VLWatMicP_vr(0,NY,NX),WatFLo2LitR_col(NY,NX),TLitrIceFlxThaw_col(NY,NX)
      write(*,*)'ice',VLiceMicPr,VLiceMicP_vr(0,NY,NX)
      write(*,*)'engy',ENGYZ,HeatFLoByWat2LitR_col(NY,NX),TLitrIceHeatFlxFrez_col(NY,NX),HeatByLitrMassChange, &
        VHeatCapacity_vr(0,NY,NX)        
      write(*,*)'vhc',VHeatCapacityLitrX,VHeatCapacityLitR,dVHeatCapacityLitR,TairK_col(NY,NX),SoilOrgM_vr(ielmc,0,NY,NX)   
      write(*,*)'tengz',ENGYZ/VHeatCapacity_vr(0,NY,NX),HeatFLoByWat2LitR_col(NY,NX)/VHeatCapacity_vr(0,NY,NX),&
        TLitrIceHeatFlxFrez_col(NY,NX)/VHeatCapacity_vr(0,NY,NX),HeatByLitrMassChange/VHeatCapacity_vr(0,NY,NX)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif  
    HEATIN_lnd = HEATIN_lnd+HeatByLitrMassChange

  ELSE
    HEATIN_lnd      = HEATIN_lnd+HeatByLitrMassChange+(TKS_vr(NUM(NY,NX),NY,NX)-TKS_vr(0,NY,NX))*VHeatCapacity_vr(0,NY,NX)
    TKS_vr(0,NY,NX) = TKS_vr(NUM(NY,NX),NY,NX)
  ENDIF

  TCS_vr(0,NY,NX)  = units%Kelvin2Celcius(TKS_vr(0,NY,NX))
  HEATIN_lnd    = HEATIN_lnd+TLitrIceHeatFlxFrez_col(NY,NX)

  end subroutine UpdateLitRPhys

!------------------------------------------------------------------------------------------
  subroutine UpdateLitRBe4RunoffM(I,J,M,NY,NX)
  !
  !update litter physical properties by processes before surface runoff
  implicit none
  integer, intent(in) :: M,NY,NX,I,J

  character(len=*), parameter :: subname='UpdateLitRBe4RunoffM'
  real(r8) :: VOLIRZ,ENGYR,VLHeatCapLitRPre
  real(r8) :: TK0Prev,TVWatIceLitR,VWatLitrZ
  real(r8) :: VLWatMicP10,VLiceMicP10,VLWatLitR,VLicelitR


  ! SURFACE RESIDUE WATER AND TEMPERATURE
  !
  call PrintInfo('beg '//subname)

  VLWatMicP10                   = VLWatMicP1_vr(0,NY,NX)
  VLiceMicP10                   = VLiceMicP1_vr(0,NY,NX)

  VLWatMicP1_vr(0,NY,NX)        = (VLWatMicP1_vr(0,NY,NX)+WatFLo2LitRM_col(NY,NX)+LitrIceFlxThaw_col(NY,NX))
  VLiceMicP1_vr(0,NY,NX)        = (VLiceMicP1_vr(0,NY,NX)-LitrIceFlxThaw_col(NY,NX)/DENSICE)

  VLairMicP1_vr(0,NY,NX)     = AZMAX1(VLPoreLitR_col(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))

  VLWatMicPM_vr(M+1,0,NY,NX) = VLWatMicP1_vr(0,NY,NX)
  VLsoiAirPM_vr(M+1,0,NY,NX) = VLairMicP1_vr(0,NY,NX)

!  VLWatLitR  = VLWatMicP_vr(0,NY,NX)+TLitrIceFlxThaw_col(NY,NX)+WatFLo2LitR_col(NY,NX)+TXGridSurfRunoff_2DH(NY,NX)
!  VLicelitR  = VLiceMicP_vr(0,NY,NX)-TLitrIceFlxThaw_col(NY,NX)/DENSICE
  
  TVWatIceLitR                  = VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX)
  XVLMobileWaterLitR_col(NY,NX) = AZMAX1(TVWatIceLitR-VWatLitRHoldCapcity_col(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ               = VLWatMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    VOLIRZ                  = VLiceMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    XVLMobileWatMicP(NY,NX) = AZMAX1(VLWatMicP1_vr(0,NY,NX)-VWatLitrZ)
    XVLiceMicP_col(NY,NX)   = AZMAX1(VLiceMicP1_vr(0,NY,NX)-VOLIRZ)
  ELSE
    XVLMobileWatMicP(NY,NX) = 0.0_r8
    XVLiceMicP_col(NY,NX)   = 0.0_r8
  ENDIF

  XVLMobileWaterLitRM(M+1,NY,NX) = XVLMobileWaterLitR_col(NY,NX)
  XVLMobileWatMicPM(M+1,NY,NX)   = XVLMobileWatMicP(NY,NX)
  XVLiceMicPM(M+1,NY,NX)         = XVLiceMicP_col(NY,NX)
  IF(VLitR_col(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat_vr(0,NY,NX)     = AZMIN1(1._r8,AZMAX1t(VLWatMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)))
    FracSoiPAsIce_vr(0,NY,NX)     = AZMAX1t(VLiceMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    AirFilledSoilPore_vr(0,NY,NX) = AZMAX1t(VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)) &
      *AZMAX1t((1.0_r8-XVLMobileWaterLitR_col(NY,NX)/VLWatheldCapSurf_col(NY,NX)))
  ELSE
    FracSoiPAsWat_vr(0,NY,NX)     = 0.0_r8
    FracSoiPAsIce_vr(0,NY,NX)     = 0.0_r8
    AirFilledSoilPore_vr(0,NY,NX) = 1.0_r8
  ENDIF

  AirFilledSoilPoreM_vr(M+1,0,NY,NX)        = AirFilledSoilPore_vr(0,NY,NX)
  VLHeatCapLitRPre           = VHeatCapacity1_vr(0,NY,NX)                      !heat capacity
  TK0Prev                    = TKSoil1_vr(0,NY,NX)                             !residual temperature
  ENGYR                      = VHeatCapacity1_vr(0,NY,NX)*TKSoil1_vr(0,NY,NX)  !initial energy content
  VHeatCapacity1_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP1_vr(0,NY,NX)+cpi*VLiceMicP1_vr(0,NY,NX)  !update heat capacity

  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    TKSoil1_vr(0,NY,NX)=(ENGYR+HeatFLoByWat2LitRM_col(NY,NX)+LitrIceHeatFlxFrez_col(NY,NX))/VHeatCapacity1_vr(0,NY,NX)
    if(TKSoil1_vr(0,NY,NX)<200._r8 .or. abs(TKSoil1_vr(0,NY,NX)-TK0Prev)>60._r8)then
      write(*,*)'IJ, weird litter temp UpdateLitRBe4RunoffM=',I*1000+J,TKSoil1_vr(0,NY,NX),TK0Prev,TairK_col(NY,NX),TKSoil1_vr(NUM(NY,NX),NY,NX)
      write(*,*)'VLHeatcap',VHeatCapacity1_vr(0,NY,NX),VLHeatCapLitRPre
      write(*,*)'engy',ENGYR/VHeatCapacity1_vr(0,NY,NX),HeatFLoByWat2LitRM_col(NY,NX)/VHeatCapacity1_vr(0,NY,NX),&
        LitrIceHeatFlxFrez_col(NY,NX)/VHeatCapacity1_vr(0,NY,NX)
      write(*,*)'cpo',cpo*SoilOrgM_vr(ielmc,0,NY,NX),cpw*VLWatMicP1_vr(0,NY,NX),cpi*VLiceMicP1_vr(0,NY,NX),VLHeatCapLitRPre,VHeatCapLitRMin_col(NY,NX) 
      write(*,*)'cpw',cpw*VLWatMicP10,VLWatMicP10,VLiceMicP10,SoilOrgM_vr(ielmc,0,NY,NX),VLWatMicP1_vr(0,NY,NX),VLiceMicP1_vr(0,NY,NX)
      write(*,*)'watflw',ENGYR,WatFLo2LitRM_col(NY,NX),HeatFLoByWat2LitRM_col(NY,NX),HeatFLoByWat2LitRM_col(NY,NX)/(WatFLo2LitRM_col(NY,NX)*cpw)
      write(*,*)'vlwat',VLWatMicP10,WatFLo2LitRM_col(NY,NX),cumWatFlx2LitRByRunoff_col(NY,NX),cumHeatFlx2LitRByRunoff_col(NY,NX)
      write(*,*)I,J,M      
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
!    IF(ABS(VHeatCapacity1_vr(0,NY,NX)/VLHeatCapLitRPre-1._r8)>0.025_r8.or. &
!      abs(TKSoil1_vr(0,NY,NX)/TK0Prev-1._r8)>0.025_r8)THEN
!      TKSoil1_vr(0,NY,NX)=TKSoil1_vr(NUM(NY,NX),NY,NX)
!    ENDIF
  ELSE
    TKSoil1_vr(0,NY,NX)=TKSoil1_vr(NUM(NY,NX),NY,NX)
  ENDIF

!  watflw(NY, NX)  = watflw(NY,NX)+WatFLo2LitRM_col(NY,NX)
!  waticefl(NY,NX) = waticefl(NY,NX)+LitrIceFlxThaw_col(NY,NX)
  call PrintInfo('end '//subname)
  end subroutine UpdateLitRBe4RunoffM

!------------------------------------------------------------------------------------------
  subroutine UpdateLitRAftRunoff(I,J,M,NY,NX)
  !
  !Description
  !Update litter water after computing runoff
  implicit none
  integer, intent(in) :: I,J,M,NY,NX
  real(r8) :: VOLIRZ,ENGYR,VLHeatCapLitRPre
  real(r8) :: TK0Prev,TVWatIceLitR,VWatLitrZ
  real(r8) :: VLWatMicP10,VLiceMicP10
  integer :: K

  ! SURFACE RESIDUE WATER AND TEMPERATURE
  
  VLWatMicP10 = VLWatMicP1_vr(0,NY,NX)
  VLiceMicP10 = VLiceMicP1_vr(0,NY,NX)

  VLWatMicP1_vr(0,NY,NX)     = AZMAX1(VLWatMicP1_vr(0,NY,NX)+cumWatFlx2LitRByRunoff_col(NY,NX))
  VLairMicP1_vr(0,NY,NX)     = AZMAX1(VLPoreLitR_col(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))

  if(M.NE.NPH)then
    VLWatMicPM_vr(M+1,0,NY,NX) = VLWatMicP1_vr(0,NY,NX)
    VLsoiAirPM_vr(M+1,0,NY,NX) = VLairMicP1_vr(0,NY,NX)
  endif

  TVWatIceLitR                  = VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX)
  XVLMobileWaterLitR_col(NY,NX) = AZMAX1(TVWatIceLitR-VWatLitRHoldCapcity_col(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ               = VLWatMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    VOLIRZ                  = VLiceMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    XVLMobileWatMicP(NY,NX) = AZMAX1(VLWatMicP1_vr(0,NY,NX)-VWatLitrZ)
    XVLiceMicP_col(NY,NX)   = AZMAX1(VLiceMicP1_vr(0,NY,NX)-VOLIRZ)
  ELSE
    XVLMobileWatMicP(NY,NX) = 0.0_r8
    XVLiceMicP_col(NY,NX)   = 0.0_r8
  ENDIF

  if(M.NE.NPH)then
    XVLMobileWaterLitRM(M+1,NY,NX) = XVLMobileWaterLitR_col(NY,NX)
    XVLMobileWatMicPM(M+1,NY,NX)   = XVLMobileWatMicP(NY,NX)
    XVLiceMicPM(M+1,NY,NX)         = XVLiceMicP_col(NY,NX)
  endif

  IF(VLitR_col(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat_vr(0,NY,NX)     = AZMIN1(1._r8,AZMAX1t(VLWatMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)))
    FracSoiPAsIce_vr(0,NY,NX)     = AZMIN1(1._r8,AZMAX1t(VLiceMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)))
    AirFilledSoilPore_vr(0,NY,NX) = AZMAX1t(VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)) &
      *AZMAX1t((1.0_r8-XVLMobileWaterLitR_col(NY,NX)/VLWatheldCapSurf_col(NY,NX)))
  ELSE
    FracSoiPAsWat_vr(0,NY,NX)     = 0.0_r8
    FracSoiPAsIce_vr(0,NY,NX)     = 0.0_r8
    AirFilledSoilPore_vr(0,NY,NX) = 1.0_r8
  ENDIF

  if(M.NE.NPH)AirFilledSoilPoreM_vr(M+1,0,NY,NX) = AirFilledSoilPore_vr(0,NY,NX)
  VLHeatCapLitRPre           = VHeatCapacity1_vr(0,NY,NX)                !heat capacity
  TK0Prev                    = TKSoil1_vr(0,NY,NX)                                 !residual temperature
  ENGYR                      = VHeatCapacity1_vr(0,NY,NX)*TKSoil1_vr(0,NY,NX)  !initial energy content
  VHeatCapacity1_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP1_vr(0,NY,NX)+cpi*VLiceMicP1_vr(0,NY,NX)  !update heat capacity

  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    TKSoil1_vr(0,NY,NX)=(ENGYR+cumHeatFlx2LitRByRunoff_col(NY,NX))/VHeatCapacity1_vr(0,NY,NX)
!    write(*,'(50A)')('-',K=1,50)
!    write(*,*)I,J,M,NY,NX
!    write(*,*)'weird litter temp UpdateLitRAftRunoff  =',TKSoil1_vr(0,NY,NX),TK0Prev,VHeatCapacity1_vr(0,NY,NX),VLHeatCapLitRPre
!    write(*,*)'engyr',ENGYR/VHeatCapacity1_vr(0,NY,NX),cumHeatFlx2LitRByRunoff_col(NY,NX)/VHeatCapacity1_vr(0,NY,NX)
!    write(*,*)'cpo',cpo*SoilOrgM_vr(ielmc,0,NY,NX),cpw*VLWatMicP1_vr(0,NY,NX),cpi*VLiceMicP1_vr(0,NY,NX),VHeatCapLitRMin_col(NY,NX) 
!    write(*,*)'cpw',SoilOrgM_vr(ielmc,0,NY,NX),VLWatMicP10,VLiceMicP10,VLWatMicP1_vr(0,NY,NX),VLiceMicP1_vr(0,NY,NX)
!    write(*,*)'vlwat',VLWatMicP10,cumWatFlx2LitRByRunoff_col(NY,NX),cumHeatFlx2LitRByRunoff_col(NY,NX)
   
    if(TKSoil1_vr(0,NY,NX)<100._r8 .or. TKSoil1_vr(0,NY,NX)>400._r8)then
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
  ELSE
    TKSoil1_vr(0,NY,NX)=TKSoil1_vr(NUM(NY,NX),NY,NX)
  ENDIF

  TKS_vr(0,NY,NX)           = TKSoil1_vr(0,NY,NX)
  VLWatMicP_vr(0,NY,NX)     = VLWatMicP1_vr(0,NY,NX)
  VLiceMicP_vr(0,NY,NX)     = VLiceMicP1_vr(0,NY,NX)
  VHeatCapacity_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
  
  end subroutine UpdateLitRAftRunoff

!------------------------------------------------------------------------------------------

  subroutine LateralGridsHdryoExch(I,J,M,NY,NX,NHE,NHW,NVS,NVN)
  !
  !between grid horizontal water flow
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M,NY,NX,NHE,NHW,NVS,NVN
  integer :: N1,N2   !source grid, with which the lateral exchange is computed 
  integer :: N,NN,N4,N5,N4B,N5B
  real(r8) :: ALT1,ALT2,ALTB
  real(r8) :: QRQ1  !equilibrium outgoing water flux
  real(r8) :: VLWatLitR,VLWatLitR1
  integer, parameter :: idirew=1
  integer, parameter :: idirns=2

!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!

  N1=NX;N2=NY

  DO  N=1,2
    DO  NN=1,2
      IF(N.EQ.iEastWestDirection)THEN
        !east-west
        IF((NX.EQ.NHE .AND. NN.EQ.iOutflow) .OR. (NX.EQ.NHW .AND. NN.EQ.iInflow))THEN
          !at the eastern/western boundary 
          cycle
        ELSE
          N4  = NX+1   !right/east
          N5  = NY
          N4B = NX-1   !left/west
          N5B = NY
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        !south-north
        IF((NY.EQ.NVS .AND. NN.EQ.iOutflow) .OR. (NY.EQ.NVN .AND. NN.EQ.iInflow))THEN
          !at the boundary
          cycle
        ELSE
          N4  = NX
          N5  = NY+1  !upper/south
          N4B = NX
          N5B = NY-1  !lower/north
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
!     FSLOPE=partitions surface water flow in (N=1)EW,(N=2)NS direction
!     XVLMobileWaterLitR_col: mobile water
!     IFLBM=runoff direction (0 = E or S, 1 = W or N)

!     there is water moving
      IF(SurfRunoffWatFluxM_2DH(M,N2,N1).GT.ZEROS(N2,N1))THEN
        ! there is runoff
        ! source grid elevation
        ALT1=Altitude_grid(N2,N1)+XVLMobileWaterLitR_col(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
!
!     EAST OR SOUTH RUNOFF
!
        IF(NN.EQ.iOutflow)THEN
          !destination grid (N5,N4) elevation
          ALT2=Altitude_grid(N5,N4)+XVLMobileWaterLitR_col(N5,N4)/AREA(3,NU(N5,N4),N5,N4)
          IF(ALT1.GT.ALT2)THEN
            !out flow from source grid into dest grid, AZMAX1 ensures upstream flow,
            QRQ1=AZMAX1((ALT1-ALT2)*AREA(3,NUM(N2,N1),N2,N1)*AREA(3,NU(N5,N4),N5,N4) &  
              -XVLMobileWaterLitR_col(N5,N4)*AREA(3,NUM(N2,N1),N2,N1) &
              +XVLMobileWaterLitR_col(N2,N1)*AREA(3,NU(N5,N4),N5,N4)) &            
              /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5,N4),N5,N4))

            !assign flux to target
            WatFlx2LitRByRunoff_2DH(N,2,N5,N4)   = AMIN1(QRQ1,SurfRunoffWatFluxM_2DH(M,N2,N1))*FSLOPE(N,N2,N1)
            HeatFlx2LitRByRunoff_2DH(N,2,N5,N4)  = cpw*TKSoil1_vr(0,N2,N1)*WatFlx2LitRByRunoff_2DH(N,2,N5,N4)
            XGridSurfRunoff_2DH(N,2,N5,N4)       = XGridSurfRunoff_2DH(N,2,N5,N4)+WatFlx2LitRByRunoff_2DH(N,2,N5,N4)
            HeatXGridBySurfRunoff_2DH(N,2,N5,N4) = HeatXGridBySurfRunoff_2DH(N,2,N5,N4)+HeatFlx2LitRByRunoff_2DH(N,2,N5,N4)
            QflxSurfRunoffM_2DH(M,N,2,N5,N4)     = WatFlx2LitRByRunoff_2DH(N,2,N5,N4)
            IFLBM(M,N,2,N5,N4)                   = 0

          ELSE
            WatFlx2LitRByRunoff_2DH(N,2,N5,N4)  = 0.0_r8
            HeatFlx2LitRByRunoff_2DH(N,2,N5,N4) = 0.0_r8
            QflxSurfRunoffM_2DH(M,N,2,N5,N4)    = 0.0_r8
            IFLBM(M,N,2,N5,N4)                  = 0
          ENDIF
        ENDIF
!
!     WEST OR NORTH RUNOFF
!
        IF(NN.EQ.iInflow)THEN
          IF(N4B.GT.0 .AND. N5B.GT.0)THEN
            !destination grid (N5B,N4B) 
            ALTB = Altitude_grid(N5B,N4B)+XVLMobileWaterLitR_col(N5B,N4B)/AREA(3,NU(N5,N4B),N5B,N4B)
            
            IF(ALT1.GT.ALTB)THEN
              QRQ1 = AZMAX1((ALT1-ALTB)*AREA(3,NUM(N2,N1),N2,N1)*AREA(3,NU(N5B,N4B),N5B,N4B) &
                -XVLMobileWaterLitR_col(N5B,N4B)*AREA(3,NUM(N2,N1),N2,N1) &
                +XVLMobileWaterLitR_col(N2,N1)*AREA(3,NU(N5B,N4B),N5B,N4B)) &
                /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5B,N4B),N5B,N4B))
              !assign flux to target
              WatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)   = AMIN1(QRQ1,SurfRunoffWatFluxM_2DH(M,N2,N1))*FSLOPE(N,N2,N1)
              HeatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)  = cpw*TKSoil1_vr(0,N2,N1)*WatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)
              XGridSurfRunoff_2DH(N,1,N5B,N4B)       = XGridSurfRunoff_2DH(N,1,N5B,N4B)+WatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)
              HeatXGridBySurfRunoff_2DH(N,1,N5B,N4B) = HeatXGridBySurfRunoff_2DH(N,1,N5B,N4B)+HeatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)
              QflxSurfRunoffM_2DH(M,N,1,N5B,N4B)     = WatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)
              IFLBM(M,N,1,N5B,N4B)                   = 1

            ELSE
              WatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)  = 0.0_r8
              HeatFlx2LitRByRunoff_2DH(N,1,N5B,N4B) = 0.0_r8
              QflxSurfRunoffM_2DH(M,N,1,N5B,N4B)    = 0.0_r8
              IFLBM(M,N,1,N5B,N4B)                  = 1
            ENDIF
          ENDIF
        ENDIF
      ELSE
        !there is no runoff
        WatFlx2LitRByRunoff_2DH(N,2,N5,N4)  = 0.0_r8
        HeatFlx2LitRByRunoff_2DH(N,2,N5,N4) = 0.0_r8
        QflxSurfRunoffM_2DH(M,N,2,N5,N4)    = 0.0_r8
        IFLBM(M,N,2,N5,N4)                  = 0
        IF(N4B.GT.0 .AND. N5B.GT.0)THEN
          WatFlx2LitRByRunoff_2DH(N,1,N5B,N4B)  = 0.0_r8
          HeatFlx2LitRByRunoff_2DH(N,1,N5B,N4B) = 0.0_r8
          QflxSurfRunoffM_2DH(M,N,1,N5B,N4B)    = 0.0_r8
          IFLBM(M,N,1,N5B,N4B)                  = 0
        ENDIF
      ENDIF
      IF(M.EQ.NPH)THEN
        IFLBH(N,NN,N5,N4)=IFLBM(M,N,NN,N5,N4)
        IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.iOutflow)THEN
          IFLBH(N,NN,N5B,N4B)=IFLBM(M,N,NN,N5B,N4B)
        ENDIF
      ENDIF

    ENDDO  
  ENDDO
  
  end subroutine LateralGridsHdryoExch


end module SurfLitterPhysMod