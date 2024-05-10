module RedistMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : padr, print_info,endrun,destroy
  use minimathmod, only : safe_adb,AZMAX1
  use EcoSiMParDataMod, only : micpar
  use SurfLitterPhysMod, only : UpdateLitRPhys
  use nitrosMod  
  use ElmIDMod
  use EcosimConst
  use FlagDataType
  use GridDataType
  use ErosionBalMod
  use RunoffBalMod
  use TFlxTypeMod
  use SoilLayerDynMod
  use SOMDataType
  USE SoilPropertyDataType
  use SoilWaterDataType
  USE SurfLitterDataType
  USE SoilHeatDataType
  USE EcoSIMCtrlDataType
  use EcoSIMCtrlMod
  USE SoilBGCDataType
  USE EcosimBGCFluxType
  use PlantDataRateType
  USE CanopyDataType
  USE ClimForcDataType
  USE SurfSoilDataType
  USE SnowDataType
  USE LandSurfDataType
  USE FertilizerDataType
  USE AqueChemDatatype
  USE ChemTranspDataType
  USE IrrigationDataType
  USE RootDataType
  USE SoilPhysDataType
  use MicrobialDataType
  USE SedimentDataType
  USE EcoSimSumDataType
  USE LateralTranspMod
  use UnitMod, only : units
  use SnowBalanceMod
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  real(r8), pointer :: ThetaCX(:)



  integer :: curday, curhour
  public :: redist, InitRedist
  contains

  subroutine InitRedist
  implicit none

  allocate(ThetaCX(micpar%NumOfPlantLitrCmplxs))

  ThetaCX=(/8.0E-06_r8,8.0E-06_r8/)

  call InitTflxType()

  end subroutine InitRedist

!------------------------------------------------------------------------------------------

  SUBROUTINE redist(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE UPDATES SOIL STATE VARIABLES WITH WATER, HEAT,
!     C, N, P, SOLUTE FLUXES CALCULATED IN EARLIER SUBROUTINES
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L,LG
  real(r8) :: DORGC(JZ),DVLiceMicP(JZ)
  real(r8) :: TXCO2(JY,JX),DORGE(JY,JX)
  real(r8) :: VOLISO,VOLPT,VOLTT
  real(r8) :: TFLWT
!     execution begins here
  curday=I
  curhour=J
  VOLISO=0.0_r8
  TFLWT=0.0_r8
  VOLPT=0.0_r8
  VOLTT=0.0_r8
  
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      TXCO2(NY,NX)=0.0_r8
      DORGE(NY,NX)=0.0_r8
      WQRH(NY,NX)=0.0_r8
!
      call AddFlux2SurfaceResidue(NY,NX)
!
      call SinkChemicals(NY,NX)
!
!     RUNOFF AND SUBSURFACE BOUNDARY FLUXES
!
      call RunoffBal(I,J,NY,NX,NHW,NHE,NVN,NVS)
!
!     CHANGE EXTERNAL WATER TABLE DEPTH THROUGH DISTURBANCE
!
      call ModifyExWTBLByDisturbance(I,J,NY,NX)

      call LateralTranspt(I,J,NY,NX,LG)

      call SnowMassUpdate(NY,NX)

      call HandleSurfaceBoundary(I,NY,NX)
!
!     OVERLAND FLOW
      call OverlandFlow(NY,NX)

      call SoilErosion(NY,NX,DORGE)
!
      call ChemicalBySnowRedistribution(NY,NX)
!
      call DiagSnowChemMass(NY,NX)
!
      call LitterLayerSummary(NY,NX)
!
      call update_physVar_Profile(NY,NX,VOLISO,DVLiceMicP)    

      call UpdateChemInSoilLays(NY,NX,LG,DORGC,TXCO2,DORGE)
!
!     SNOWPACK LAYERING
      call SnowpackLayering(NY,NX)

      call RelayerSoilProfile(NY,NX,DORGC,DVLiceMicP)

      call UpdateOutputVars(I,J,NY,NX,TXCO2)
!
!     CHECK MATERIAL BALANCES
!
!      IF(I.EQ.365.AND.J.EQ.24)THEN
!        WRITE(19,2221)'ORGC',I,J,iYearCurrent,NX,NY,(SoilOrgM_vr(ielmc,L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
!        WRITE(20,2221)'ORGN',I,J,iYearCurrent,NX,NY,(SoilOrgM_vr(ielmn,L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
2221    FORMAT(A8,5I6,21E14.6)
!      ENDIF

    ENDDO D9990
  ENDDO D9995
  RETURN

  END subroutine redist

!------------------------------------------------------------------------------------------

  subroutine UpdateOutputVars(I,J,NY,NX,TXCO2)

  use TillageMixMod
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: TXCO2(JY,JX)   !what does TXCO2 mean, be careful?
  real(r8) :: VLSoilPoreMicP_vrX,VOLTX
  integer  :: L

  Eco_NetRad_col(NY,NX)=Eco_NetRad_col(NY,NX)+HeatByRadiation(NY,NX)
  Eco_Heat_Latent_col(NY,NX)=Eco_Heat_Latent_col(NY,NX)+HeatEvapAir2Surf(NY,NX)
  Eco_Heat_Sens_col(NY,NX)=Eco_Heat_Sens_col(NY,NX)+HeatSensAir2Surf(NY,NX)
  Eco_Heat_Grnd_col(NY,NX)=Eco_Heat_Grnd_col(NY,NX)-(HeatNet2Surf(NY,NX)-HeatSensVapAir2Surf(NY,NX))
  Canopy_Heat_Latent_col(NY,NX)=Canopy_Heat_Latent_col(NY,NX)+HeatEvapAir2Surf(NY,NX)*BndlResistCanG(NY,NX)
  Canopy_Heat_Sens_col(NY,NX)=Canopy_Heat_Sens_col(NY,NX)+HeatSensAir2Surf(NY,NX)*BndlResistCanG(NY,NX)
  Eco_NEE_col(NY,NX)=Canopy_NEE_col(NY,NX)+SurfGasFlx_col(idg_CO2,NY,NX)
  ECO_ER_col(NY,NX)=ECO_ER_col(NY,NX)+SurfGasFlx_col(idg_CO2,NY,NX)
  Eco_NPP_col(NY,NX)=Eco_GPP_col(NY,NX)+Eco_AutoR_col(NY,NX)
  Eco_NBP_col(NY,NX)=Eco_NBP_col(NY,NX)+Canopy_NEE_col(NY,NX) &
    +SurfGasFlx_col(idg_CO2,NY,NX)+SurfGasFlx_col(idg_CH4,NY,NX) &
    +TXCO2(NY,NX)-HydroSufDOCFlx_col(NY,NX)-HydroSufDICFlx_col(NY,NX)-HydroSubsDOCFlx_col(NY,NX)-HydroSubsDICFlx_col(NY,NX)
    
  IF(NU(NY,NX).GT.NUI(NY,NX))THEN  !the surface is lowered
    DO L=NUI(NY,NX),NU(NY,NX)-1
      IF(VLSoilPoreMicP_vr(L,NY,NX).LE.ZEROS2(NY,NX))THEN
        TKS(L,NY,NX)=TKS(NU(NY,NX),NY,NX)
        TCS(L,NY,NX)=units%Kelvin2Celcius(TKS(L,NY,NX))
      ENDIF
    ENDDO
  ENDIF
  !
  !     MIX ALL SOIL STATE VARIABLES AND INCORPORATE ALL SURFACE
  !     RESIDUE STATE VARIABLES WITHIN THE TILLAGE ZONE TO THE EXTENT
  !     ASSOCIATED IN 'DAY' WITH EACH TILLAGE EVENT ENTERED IN THE
  !     TILLAGE FILE
  call ApplyTillageMixing(I,J,NY,NX)

  !
  !     OUTPUT FOR SOIL WATER, ICE CONTENTS
  !
  THETWZ(0,NY,NX)=AZMAX1((VLWatMicP_vr(0,NY,NX)-VWatLitRHoldCapcity(NY,NX))/AREA(3,0,NY,NX))
  THETIZ(0,NY,NX)=AZMAX1((VLiceMicP(0,NY,NX)-VWatLitRHoldCapcity(NY,NX))/AREA(3,0,NY,NX))
  !THETWZ(0,NY,NX)=AZMAX1(AMIN1(1.0,VLWatMicP_vr(0,NY,NX)/VLitR(NY,NX)))
  !THETIZ(0,NY,NX)=AZMAX1(AMIN1(1.0,VLiceMicP(0,NY,NX)/VLitR(NY,NX)))
  
  D9945: DO L=NUI(NY,NX),NL(NY,NX)
    VLSoilPoreMicP_vrX=AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)*FracSoiAsMicP(L,NY,NX)
    VOLTX=VLSoilPoreMicP_vrX+VLMacP(L,NY,NX)
    THETWZ(L,NY,NX)=safe_adb(VLWatMicP_vr(L,NY,NX)+AMIN1(VLMacP(L,NY,NX),&
      VLWatMacP(L,NY,NX)),VOLTX)
    THETIZ(L,NY,NX)=safe_adb(VLiceMicP(L,NY,NX)+AMIN1(VLMacP(L,NY,NX) &
        ,VLiceMacP(L,NY,NX)),VOLTX)
  ENDDO D9945
  end subroutine UpdateOutputVars

!------------------------------------------------------------------------------------------

  subroutine DiagSnowChemMass(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: SSW,ENGYW,WS
  integer :: L,nsalts
  !
  !     UPDATE STATE VARIABLES WITH TOTAL FLUXES CALCULATED ABOVE
  !
  !     IF(J.EQ.24)THEN
  !
  !     SNOWPACK VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
  !     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
  !     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS
  !
  DO  L=1,JS
    WS=VLDrySnoWE(L,NY,NX)+VLWatSnow(L,NY,NX)+VLIceSnow(L,NY,NX)*DENSICE

    WaterStoreLandscape=WaterStoreLandscape+WS
    UVLWatMicP(NY,NX)=UVLWatMicP(NY,NX)+WS
    ENGYW=VLHeatCapSnow(L,NY,NX)*TKSnow(L,NY,NX)
    HeatStoreLandscape=HeatStoreLandscape+ENGYW
    TGasC_lnd=TGasC_lnd+trcg_solsml(idg_CO2,L,NY,NX)+trcg_solsml(idg_CH4,L,NY,NX)
    DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)+trcg_solsml(idg_CO2,L,NY,NX)+trcg_solsml(idg_CH4,L,NY,NX)
    OXYGSO=OXYGSO+trcg_solsml(idg_O2,L,NY,NX)
    TGasN_lnd=TGasN_lnd+trcg_solsml(idg_N2,L,NY,NX)+trcg_solsml(idg_N2O,L,NY,NX)
    TDisolNH4_lnd=TDisolNH4_lnd+trcn_solsml(ids_NH4,L,NY,NX)+trcg_solsml(idg_NH3,L,NY,NX)
    tNO3_lnd=tNO3_lnd+trcn_solsml(ids_NO3,L,NY,NX)
    TDisolPi_lnd=TDisolPi_lnd+trcn_solsml(ids_H1PO4,L,NY,NX)+trcn_solsml(ids_H2PO4,L,NY,NX)

    IF(salt_model)THEN
      SSW=0._r8
      do nsalts=idsalt_beg,idsalt_end
        SSW=SSW+trcs_solsml(nsalts,L,NY,NX)*trcSaltIonNumber(nsalts)
      ENDDO  
      TION=TION+SSW
    ENDIF
  ENDDO
  end subroutine DiagSnowChemMass
!------------------------------------------------------------------------------------------

  subroutine ModifyExWTBLByDisturbance(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: DCORPW
!     begin_execution
!     SolarNoonHour_col=hour of solar noon from weather file
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     DCORP=mixing intensity (fire) or depth (tillage,drainage) of disturbance
!     CumDepth2LayerBottom(NU=soil surface elevation
!     DTBLI,DTBLDI=depth of natural,artificial water table from readi.f
!     ExtWaterTable,ExtWaterTablet0=current,initial natural water table depth
!     DTBLY,DTBLD=current,initial artificial water table depth
!     IDWaterTable=water table flag from readi.f
!        :0=none
!        :1,2=natural stationary,mobile
!        :3,4=artificial stationary,mobile
!     FWatDischarge=hourly water loss through lateral and lower boundaries
!
  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.iSoilDisturbType_col(I,NY,NX).EQ.23)THEN
    ! drainage is on
    DCORPW=DCORP(I,NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
    DTBLI(NY,NX)=DCORPW
    ExtWaterTablet0(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX))
    ExtWaterTable(NY,NX)=ExtWaterTablet0(NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
  ENDIF

  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.iSoilDisturbType_col(I,NY,NX).EQ.24)THEN
    ! drainage in on
    DCORPW=DCORP(I,NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
    IF(IDWaterTable(NY,NX).EQ.1)THEN
      IDWaterTable(NY,NX)=3
    ELSEIF(IDWaterTable(NY,NX).EQ.2)THEN
      IDWaterTable(NY,NX)=4
    ENDIF
    DTBLDI(NY,NX)=DCORPW
    DTBLD(NY,NX)=AZMAX1(DTBLDI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX)))
    DTBLY(NY,NX)=DTBLD(NY,NX)
  ENDIF
!
!     SET DEPTH OF MOBILE EXTERNAL WATER TABLE
! switched on for change of water table due to discharge/drainage
! why 0.00167, time relaxization constant, ?
! 4 is mobile tile drainge.
  IF(IDWaterTable(NY,NX).EQ.2.OR.IDWaterTable(NY,NX).EQ.4)THEN
    ExtWaterTable(NY,NX)=ExtWaterTable(NY,NX)-FWatDischarge(NY,NX)/AREA(3,NU(NY,NX),NY,NX) &
      -0.00167_r8*(ExtWaterTable(NY,NX)-ExtWaterTablet0(NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX))
    ExtWaterTable(NY,NX)=ExtWaterTablet0(NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
  ENDIF
  
  IF(IDWaterTable(NY,NX).EQ.4)THEN
    DTBLY(NY,NX)=DTBLY(NY,NX)-FWatDischarge(NY,NX)/AREA(3,NU(NY,NX),NY,NX) &
      -0.00167_r8*(DTBLY(NY,NX)-DTBLD(NY,NX))
  ENDIF
  end subroutine ModifyExWTBLByDisturbance
!------------------------------------------------------------------------------------------

  subroutine HandleSurfaceBoundary(I,NY,NX)

  implicit none
  integer, intent(in) :: I,NY,NX

  integer :: L,K,LS,NTG,NTP,NTX,nsalts
  real(r8):: vhcp1s
  real(r8) :: CI,CH,CO,CX
  real(r8) :: OI,OO
  real(r8) :: HI,HO
  real(r8) :: PI,PXB
  real(r8) :: SIN,SGN,SIP,SNB
  real(r8) :: SPB,SNM0,SPM0,SIR,SII,SBU
  real(r8) :: WI,WO
  real(r8) :: ZSI,ZXB,ZGI
  real(r8) :: ZNGGIN,ZN2OIN,ZNH3IN
  integer :: idom,ids,idg
  ! begin_execution

  !
  call UpdateLitRPhys(NY,NX,TWat2GridBySurfRunoff(NY,NX),THeat2GridBySurfRunoff(NY,NX),&
    HeatStoreLandscape,HEATIN)

  !     UVLWatMicP(NY,NX)=UVLWatMicP(NY,NX)-VLWatMicP_vr(0,NY,NX)-VLiceMicP(0,NY,NX)*DENSICE
  !
  !     SURFACE BOUNDARY WATER FLUXES
  !
  WI=PrecAtm_col(NY,NX)+IrrigSurface(NY,NX)   !total incoming water flux=rain/snowfall + irrigation
  CRAIN=CRAIN+WI
  URAIN(NY,NX)=URAIN(NY,NX)+WI
  WO=VapXAir2GSurf(NY,NX)+TEVAPP(NY,NX) !total outgoing water flux
  CEVAP=CEVAP-WO
  UEVAP(NY,NX)=UEVAP(NY,NX)-WO
  EvapoTransp_col(NY,NX)=-WO
  VOLWOU=VOLWOU-IrrigSubsurf_col(NY,NX)
  FWatDischarge(NY,NX)=FWatDischarge(NY,NX)-IrrigSubsurf_col(NY,NX)
  AnualH2OLoss_col(NY,NX)=AnualH2OLoss_col(NY,NX)-IrrigSubsurf_col(NY,NX)
  UDRAIN(NY,NX)=UDRAIN(NY,NX)+WaterFlowSoiMicP(3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY HEAT FLUXES
  !
  HEATIN=HEATIN+cpw*TairK(NY,NX)*PrecRainAndSurfirrig(NY,NX)+cps*TairK(NY,NX)*SnoFalPrec(NY,NX)
  HEATIN=HEATIN+HeatNet2Surf(NY,NX)+THFLXC(NY,NX)
  D5150: DO L=1,JS
    HEATIN=HEATIN+XPhaseChangeHeatL(L,NY,NX)
  ENDDO D5150
  HEATOU=HEATOU-cpw*TairK(NY,NX)*IrrigSubsurf_col(NY,NX)
!
! SURFACE BOUNDARY CO2, CH4 AND DOC FLUXES
! XCODFS: surface - atmosphere CO2 dissolution (+ve) - volatilization (-ve)
! XCOFLG: gaseous CO2 flux, [g d-2 h-1]
! TCO2Z: total root CO2 content
! FLQGQ: precipitation flux into soil surface
! Rain2LitRSurf_col: precipitation flux into surface litter
! FLQGI: irrifation flux into soil surface
! Irrig2LitRSurf: irrigation flux into surface litter
! XCODFG: soil CO2 dissolution (+ve) - volatilization (-ve)
! XCODFR: soil surface CO2 dissolution (+ve) - volatilization
! UCO2G: total soil CO2 flux, [g d-2]
! HCO2G: hourly soil CO2 flux, [g d-2 h-1]
  CI=GasSfAtmFlx_col(idg_CO2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_CO2,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_CO2,NY,NX) &
    +(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*CO2_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*CO2_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_CO2,0,NY,NX)+trcg_surf_disevap_flx(idg_CO2,NY,NX)
  CH=GasSfAtmFlx_col(idg_CH4,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_CH4,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_CH4,NY,NX) &
    +(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*CH4_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*CH4_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_CH4,0,NY,NX)+trcg_surf_disevap_flx(idg_CH4,NY,NX)
  CO=-IrrigSubsurf_col(NY,NX)*CO2_irrig_conc(NY,NX)
  CX=-IrrigSubsurf_col(NY,NX)*CH4_irrig_conc(NY,NX)
  SurfGasFlx_col(idg_CO2,NY,NX)=SurfGasFlx_col(idg_CO2,NY,NX)+CI
  SurfGasFlx_col(idg_CH4,NY,NX)=SurfGasFlx_col(idg_CH4,NY,NX)+CH
  CO2GIN=CO2GIN+CI+CH
  TOMOU(ielmc)=TOMOU(ielmc)+CO+CX
  !
  !     SURFACE BOUNDARY O2 FLUXES
  !
  OI=GasSfAtmFlx_col(idg_O2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_O2,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_O2,NY,NX) &
    +(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*O2_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*O2_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_O2,0,NY,NX)+trcg_surf_disevap_flx(idg_O2,NY,NX)
  OXYGIN=OXYGIN+OI
  OO=trcg_RMicbTransf_vr(idg_O2,0,NY,NX)-IrrigSubsurf_col(NY,NX)*O2_irrig_conc(NY,NX)
  OXYGOU=OXYGOU+OO
  SurfGasFlx_col(idg_O2,NY,NX)=SurfGasFlx_col(idg_O2,NY,NX)+OI
  HI=GasSfAtmFlx_col(idg_H2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_H2,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_H2,NY,NX) &
    +Gas_Disol_Flx_vr(idg_H2,0,NY,NX)+trcg_surf_disevap_flx(idg_H2,NY,NX)
  H2GIN=H2GIN+HI
  HO=trcg_RMicbTransf_vr(idg_H2,0,NY,NX)
  H2GOU=H2GOU+HO
  !
  !     SURFACE BOUNDARY N2, N2O, NH3, NH4, NO3, AND DON FLUXES
  !
  ZSI=((Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX)) &
      *(NH4_rain_conc(NY,NX)+NH3_rain_conc(NY,NX)+NO3_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX)))*14.0
  ZXB=-IrrigSubsurf_col(NY,NX)*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX))-IrrigSubsurf_col(NY,NX) &
      *(NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX))*14.0
  TZIN=TZIN+ZSI
  TOMOU(ielmn)=TOMOU(ielmn)+ZXB
  ZGI=(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*(N2_rain_conc(NY,NX)+N2O_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX)) &
      +GasSfAtmFlx_col(idg_N2,NY,NX)+GasSfAtmFlx_col(idg_N2O,NY,NX)+GasSfAtmFlx_col(idg_NH3,NY,NX) &
      +GasSfAtmFlx_col(idg_NH3B,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX) &
      +Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX)+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
      +TRootGasLossDisturb_pft(idg_N2O,NY,NX)+TRootGasLossDisturb_pft(idg_NH3,NY,NX) &
      +Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,0,NY,NX) &
      +trcg_surf_disevap_flx(idg_N2,NY,NX)+trcg_surf_disevap_flx(idg_N2O,NY,NX) &
      +trcg_surf_disevap_flx(idg_NH3,NY,NX)
  ZN2GIN=ZN2GIN+ZGI
  ZDRAIN(NY,NX)=ZDRAIN(NY,NX)+trcs_3DTransp2MicP_vr(ids_NH4,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP_vr(idg_NH3,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP_vr(ids_NO3,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP_vr(ids_NO2,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP_vr(ids_NH4B,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP_vr(idg_NH3B,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP_vr(ids_NO3B,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP_vr(ids_NO2B,3,NK(NY,NX),NY,NX)

  ZNGGIN=GasSfAtmFlx_col(idg_N2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX)
  ZN2OIN=GasSfAtmFlx_col(idg_N2O,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX)+Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)
  ZNH3IN=GasSfAtmFlx_col(idg_NH3,NY,NX)+GasSfAtmFlx_col(idg_NH3B,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
    +Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)

  SurfGasFlx_col(idg_N2,NY,NX)=SurfGasFlx_col(idg_N2,NY,NX)+ZNGGIN
  SurfGasFlx_col(idg_N2O,NY,NX)=SurfGasFlx_col(idg_N2O,NY,NX)+ZN2OIN
  SurfGasFlx_col(idg_NH3,NY,NX)=SurfGasFlx_col(idg_NH3,NY,NX)+ZNH3IN
  SurfGasFlx_col(idg_N2,NY,NX)=SurfGasFlx_col(idg_N2,NY,NX)+Micb_N2Fixation_vr(0,NY,NX)
  SurfGasFlx_col(idg_H2,NY,NX)=SurfGasFlx_col(idg_H2,NY,NX)+HI
  !
  !     SURFACE BOUNDARY PO4 AND DOP FLUXES
  !
  PI=patomw*((Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX)) &
      *(H2PO4_rain_conc(NY,NX)+HPO4_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(H2PO4_irrig_conc(I,NY,NX)+HPO4_irrig_conc(I,NY,NX)))
  PXB=-patomw*IrrigSubsurf_col(NY,NX)*(H2PO4_irrig_conc(I,NY,NX)+HPO4_irrig_conc(I,NY,NX))
  TPIN=TPIN+PI
  TOMOU(ielmp)=TOMOU(ielmp)+PXB
  PDRAIN(NY,NX)=PDRAIN(NY,NX)+trcs_3DTransp2MicP_vr(ids_H2PO4,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP_vr(ids_H2PO4B,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP_vr(ids_H1PO4,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP_vr(ids_H1PO4B,3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY ION FLUXES
  !
  SIN=((Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX)) &
      *(2.0_r8*NH4_rain_conc(NY,NX)+NH3_rain_conc(NY,NX)+NO3_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(2.0*NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX)))
  SGN=(2.0_r8*(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*(N2_rain_conc(NY,NX)+N2O_rain_conc(NY,NX)) &
      +2.0_r8*(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX)) &
      +2.0_r8*(GasSfAtmFlx_col(idg_N2,NY,NX)+GasSfAtmFlx_col(idg_N2O,NY,NX))+GasSfAtmFlx_col(idg_NH3,NY,NX) &
      +GasSfAtmFlx_col(idg_NH3B,NY,NX)+2.0_r8*(Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX) &
      +Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX))+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
      +2.0_r8*TRootGasLossDisturb_pft(idg_N2O,NY,NX)+TRootGasLossDisturb_pft(idg_NH3,NY,NX) &
      +2.0_r8*(Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX))+Gas_Disol_Flx_vr(idg_NH3,0,NY,NX) &
      +2.0_r8*(trcg_surf_disevap_flx(idg_N2,NY,NX)+trcg_surf_disevap_flx(idg_N2O,NY,NX))+trcg_surf_disevap_flx(idg_NH3,NY,NX))/natomw
  SIP=((Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*(3.0_r8*H2PO4_rain_conc(NY,NX)+2.0_r8*HPO4_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(3.0_r8*H2PO4_irrig_conc(I,NY,NX)+2.0_r8*HPO4_irrig_conc(I,NY,NX)))
  SNB=-IrrigSubsurf_col(NY,NX)*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX))-IrrigSubsurf_col(NY,NX) &
      *(2.0_r8*NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX))
      SPB=-IrrigSubsurf_col(NY,NX)*(3.0*H2PO4_irrig_conc(I,NY,NX)+2.0*HPO4_irrig_conc(I,NY,NX))
  SNM0=(2.0_r8*RNutMicbTransf_vr(ids_NH4,0,NY,NX)+RNutMicbTransf_vr(ids_NO3,0,NY,NX)+RNutMicbTransf_vr(ids_NO2,0,NY,NX) &
      -2.0_r8*Micb_N2Fixation_vr(0,NY,NX))/natomw
  SPM0=(2.0_r8*RNutMicbTransf_vr(ids_H1PO4,0,NY,NX)+3.0_r8*RNutMicbTransf_vr(ids_H2PO4,0,NY,NX))/patomw
  !
  !     ACCUMULATE PLANT LitrFall FLUXES
  !
  Litrfall_lnds(ielmc)=Litrfall_lnds(ielmc)+LitrFallStrutElms_col(ielmc,NY,NX)
  Litrfall_lnds(ielmn)=Litrfall_lnds(ielmn)+LitrFallStrutElms_col(ielmn,NY,NX)
  Litrfall_lnds(ielmp)=Litrfall_lnds(ielmp)+LitrFallStrutElms_col(ielmp,NY,NX)

  LiterfalOrgM_col(ielmc,NY,NX)=LiterfalOrgM_col(ielmc,NY,NX)+LitrFallStrutElms_col(ielmc,NY,NX)
  LiterfalOrgM_col(ielmn,NY,NX)=LiterfalOrgM_col(ielmn,NY,NX)+LitrFallStrutElms_col(ielmn,NY,NX)
  LiterfalOrgM_col(ielmp,NY,NX)=LiterfalOrgM_col(ielmp,NY,NX)+LitrFallStrutElms_col(ielmp,NY,NX)
  !
  !     SURFACE BOUNDARY SALT FLUXES FROM RAINFALL AND SURFACE IRRIGATION
  !
  IF(salt_model)THEN
    SIR=0._r8
    SII=0._r8
    do nsalts=idsalt_beg,idsalt_end
      SIR=SIR+trcsalt_rain_conc(nsalts,NY,NX)*trcSaltIonNumber(nsalts)
      SII=SII+trcsalt_irrig_conc(idsalt_Al,I,NY,NX)*trcSaltIonNumber(nsalts)
    ENDDO  
    SIR=SIR*PrecAtm_col(NY,NX)
    SBU=-IrrigSubsurf_col(NY,NX)*SII
    SII=SII*IrrigSurface(NY,NX)
    TIONIN=TIONIN+SIR+SII
    !
    !     SUBSURFACE BOUNDARY SALT FLUXES FROM SUBSURFACE IRRIGATION
    !
    TIONOU=TIONOU+SBU
  ENDIF
  !
  ! GAS EXCHANGE FROM SURFACE VOLATILIZATION-DISSOLUTION
  !
  D9680: DO K=1,micpar%NumOfLitrCmplxs 
    do idom=idom_beg,idom_end   
      DOM_vr(idom,K,0,NY,NX)=DOM_vr(idom,K,0,NY,NX)+DOM_3DMicp_Transp_flx(idom,K,3,0,NY,NX)
    enddo
  ENDDO D9680

  do idg=idg_beg,idg_NH3-1
    trc_solml_vr(idg,0,NY,NX)=trc_solml_vr(idg,0,NY,NX)+trcg_surf_disevap_flx(idg,NY,NX) &
      +trcs_3DTransp2MicP_vr(idg,3,0,NY,NX)+Gas_Disol_Flx_vr(idg,0,NY,NX)-trcg_RMicbTransf_vr(idg,0,NY,NX)
  enddo

  trc_solml_vr(idg_N2,0,NY,NX)=trc_solml_vr(idg_N2,0,NY,NX)-Micb_N2Fixation_vr(0,NY,NX)

  trc_solml_vr(idg_NH3,0,NY,NX)=trc_solml_vr(idg_NH3,0,NY,NX)+trcg_surf_disevap_flx(idg_NH3,NY,NX) &
    +trcs_3DTransp2MicP_vr(idg_NH3,3,0,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)+TR_NH3_soil_vr(0,NY,NX)

  do ids=ids_nut_beg,ids_nuts_end
    trc_solml_vr(ids,0,NY,NX)=trc_solml_vr(ids,0,NY,NX)+trcs_3DTransp2MicP_vr(ids,3,0,NY,NX) &
      +RNutMicbTransf_vr(ids,0,NY,NX)+trcn_RChem_soil_vr(ids,0,NY,NX)
  enddo  

!update aqueous concentrations
  DO NTG=idg_beg,idg_end
    trc_solml_vr(NTG,NU(NY,NX),NY,NX)=trc_solml_vr(NTG,NU(NY,NX),NY,NX)+GasSfAtmFlx_col(NTG,NY,NX)
  ENDDO

  Eco_HR_col(NY,NX)=Eco_HR_col(NY,NX)+trcg_RMicbTransf_vr(idg_CO2,0,NY,NX)+trcg_RMicbTransf_vr(idg_CH4,0,NY,NX)
  SurfGasFlx_col(idg_N2,NY,NX)=SurfGasFlx_col(idg_N2,NY,NX)+trcg_RMicbTransf_vr(idg_N2,0,NY,NX)

  RO2GasXchangePrev_vr(0,NY,NX)=Gas_Disol_Flx_vr(idg_O2,0,NY,NX)
  RCO2GasFlxPrev_vr(0,NY,NX)=Gas_Disol_Flx_vr(idg_CO2,0,NY,NX)
  RCH4F(0,NY,NX)=Gas_Disol_Flx_vr(idg_CH4,0,NY,NX)

  !RO2AquaXchangePrev_vr:=soil surface O2 dissolution + aqueous O2 flux micropore
  RO2AquaXchangePrev_vr(0,NY,NX)=trcg_surf_disevap_flx(idg_O2,NY,NX)+trcs_3DTransp2MicP_vr(idg_O2,3,0,NY,NX) &
    -(Rain2LitRSurf_col(NY,NX)*O2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*O2_irrig_conc(NY,NX))

  RCH4PhysexchPrev_vr(0,NY,NX)=trcg_surf_disevap_flx(idg_CH4,NY,NX)+trcs_3DTransp2MicP_vr(idg_CH4,3,0,NY,NX) &
    -(Rain2LitRSurf_col(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*CH4_irrig_conc(NY,NX))
    
  RO2AquaXchangePrev_vr(NU(NY,NX),NY,NX)=RO2AquaXchangePrev_vr(NU(NY,NX),NY,NX)+GasSfAtmFlx_col(idg_O2,NY,NX)
  RCH4PhysexchPrev_vr(NU(NY,NX),NY,NX)=RCH4PhysexchPrev_vr(NU(NY,NX),NY,NX)+GasSfAtmFlx_col(idg_CH4,NY,NX)
  !
  !     SURFACE LITTER ION EXCHANGE AND PRECIPITATION
  !

  trcx_solml(idx_NH4,0,NY,NX)=trcx_solml(idx_NH4,0,NY,NX)+trcx_TRSoilChem_vr(idx_NH4,0,NY,NX)
  DO NTX=idx_AEC+1,idx_anion_soil_end
    trcx_solml(NTX,0,NY,NX)=trcx_solml(NTX,0,NY,NX)+trcx_TRSoilChem_vr(NTX,0,NY,NX)
  ENDDO

  DO NTP=idsp_psoi_beg,idsp_psoi_end
    trcp_salml(NTP,0,NY,NX)=trcp_salml(NTP,0,NY,NX)+trcp_RChem_soil(NTP,0,NY,NX)
  ENDDO
  !  !

  end subroutine HandleSurfaceBoundary
!------------------------------------------------------------------------------------------

  subroutine OverlandFlow(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,NTG,idom
  !     begin_execution
  !
  IF(ABS(TWat2GridBySurfRunoff(NY,NX)).GT.ZEROS(NY,NX))THEN
    !
    !     DOC, DON, DOP
    !
    D8570: DO K=1,micpar%NumOfLitrCmplxs    
      DO idom=idom_beg,idom_end
        DOM_vr(idom,K,0,NY,NX)=DOM_vr(idom,K,0,NY,NX)+TOMQRS(idom,K,NY,NX)
      ENDDO
    ENDDO D8570
!
    call OverlandSnowFlow(NY,NX)

  ENDIF
  end subroutine OverlandFlow
!------------------------------------------------------------------------------------------

  subroutine SoilErosion(NY,NX,DORGE)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(out) :: DORGE(JY,JX)
  integer :: K,NO,M,NGL,NTX,NTP,MID,NE,idom

  REAL(R8) :: DORGP

  ! begin_execution
  !
  ! INTERNAL SURFACE SEDIMENT TRANSPORT
  !
  IF((iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros) &
    .AND.ABS(tErosionSedmLoss(NY,NX)).GT.ZEROS(NY,NX))THEN
    TSED(NY,NX)=TSED(NY,NX)+tErosionSedmLoss(NY,NX)
    !
    !     SOIL MINERAL FRACTIONS
    !
    SAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)+TSANER(NY,NX)
    SILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)+TSILER(NY,NX)
    CLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)+TCLAER(NY,NX)
    !
    !     FERTILIZER POOLS
!
    FertN_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)+TNH4Eros_col(NY,NX)
    FertN_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)+TNH3Eros_col(NY,NX)
    FertN_soil_vr(ifert_urea,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_urea,NU(NY,NX),NY,NX)+TNUreaEros_col(NY,NX)
    FertN_soil_vr(ifert_no3,NU(NY,NX),NY,NX)=FertN_soil_vr(ifert_no3,NU(NY,NX),NY,NX)+TNO3Eros_col(NY,NX)
    FertN_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)+TNH4ErosBand_col(NY,NX)
    FertN_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)+TNH3ErosBand_col(NY,NX)
    FertN_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX)+TNUreaErosBand_col(NY,NX)
    FertN_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)=FertN_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)+TNO3ErosBand_col(NY,NX)
!
    !   EXCHANGEABLE CATIONS AND ANIONS
!
    DO NTX=idx_beg,idx_end
      trcx_solml(NTX,NU(NY,NX),NY,NX)=trcx_solml(NTX,NU(NY,NX),NY,NX)+trcx_TER(NTX,NY,NX)
    ENDDO

!
!     PRECIPITATES
!
    DO NTP=idsp_beg,idsp_end
      trcp_salml(NTP,NU(NY,NX),NY,NX)=trcp_salml(NTP,NU(NY,NX),NY,NX)+trcp_TER(NTP,NY,NX)
    ENDDO
!
!   ORGANIC CONSTITUENTS
!
    DORGP=0.0_r8
    D9280: DO K=1,jcplx
      DO  NO=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(NO),JGnfo(NO)
          DO  M=1,nlbiomcp
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              mBiomeHeter_vr(NE,MID,K,NU(NY,NX),NY,NX)=mBiomeHeter_vr(NE,MID,K,NU(NY,NX),NY,NX)+TOMEERhetr(NE,MID,K,NY,NX)
            ENDDO
            DORGE(NY,NX)=DORGE(NY,NX)+TOMEERhetr(ielmc,MID,K,NY,NX)
            DORGP=DORGP+TOMEERhetr(ielmp,MID,K,NY,NX)
          enddo
        enddo
      enddo
    ENDDO D9280

    DO  NO=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGniA(NO),JGnfA(NO)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor_vr(NE,MID,NU(NY,NX),NY,NX)=mBiomeAutor_vr(NE,MID,NU(NY,NX),NY,NX)+TOMEERauto(NE,MID,NY,NX)
          ENDDO
          DORGE(NY,NX)=DORGE(NY,NX)+TOMEERauto(ielmc,MID,NY,NX)
          DORGP=DORGP+TOMEERauto(ielmp,MID,NY,NX)
        enddo
      enddo
    enddo

    D9275: DO K=1,jcplx
      D9270: DO M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          OMBioResdu_vr(NE,M,K,NU(NY,NX),NY,NX)=OMBioResdu_vr(NE,M,K,NU(NY,NX),NY,NX)+TORMER(NE,M,K,NY,NX)
        ENDDO
        DORGE(NY,NX)=DORGE(NY,NX)+TORMER(ielmc,M,K,NY,NX)
        DORGP=DORGP+TORMER(ielmp,M,K,NY,NX)
      ENDDO D9270
      DO idom=idom_beg,idom_end
        SorbedOM_vr(idom,K,NU(NY,NX),NY,NX)=SorbedOM_vr(idom,K,NU(NY,NX),NY,NX)+TOHMER(idom,K,NY,NX)
      ENDDO
      DORGE(NY,NX)=DORGE(NY,NX)+TOHMER(idom_doc,K,NY,NX)+TOHMER(idom_acetate,K,NY,NX)
      DORGP=DORGP+TOHMER(idom_dop,K,NY,NX)
      D9265: DO M=1,jsken
        SolidOMAct_vr(M,K,NU(NY,NX),NY,NX)=SolidOMAct_vr(M,K,NU(NY,NX),NY,NX)+TOSAER(M,K,NY,NX)
        DO NE=1,NumPlantChemElms
          SolidOM_vr(NE,M,K,NU(NY,NX),NY,NX)=SolidOM_vr(NE,M,K,NU(NY,NX),NY,NX)+TOSMER(NE,M,K,NY,NX)
        ENDDO
        DORGE(NY,NX)=DORGE(NY,NX)+TOSMER(ielmc,M,K,NY,NX)
        DORGP=DORGP+TOSMER(ielmp,M,K,NY,NX)
      ENDDO D9265
    ENDDO D9275
  ENDIF
  end subroutine SoilErosion

!------------------------------------------------------------------------------------------

  subroutine LitterLayerSummary(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: K,N,M,NGL,NB,NE
  real(r8) :: PSS,HS,CS
  real(r8) :: litrOM(1:NumPlantChemElms),OS
  real(r8) :: POS,POX,POP
  real(r8) :: SSS,ZG,Z4S,WS,Z4X
  real(r8) :: Z4F,ZOS,ZOF

!     begin_execution
!     TOTAL C,N,P, SALTS IN SURFACE RESIDUE
!

  call sumSurfOMCK(NY,NX,RC0(:,NY,NX),RC0ff(NY,NX))

  call sumMicBiomLayL(0,NY,NX,tMicBiome_col(1:NumPlantChemElms,NY,NX))
  
  call sumLitrOMLayL(0,NY,NX,litrOM)

  SoilOrgM_vr(1:NumPlantChemElms,0,NY,NX)=litrOM

  DO NE=1,NumPlantChemElms
    tSoilOrgM_col(NE,NY,NX)=tSoilOrgM_col(NE,NY,NX)+SoilOrgM_vr(NE,0,NY,NX)
  enddo

  OMLitrC_vr(0,NY,NX)=litrOM(ielmc)

  DO NE=1,NumPlantChemElms
    LitRMStoreLndscap(NE)=LitRMStoreLndscap(NE)+litrOM(NE)
    tLitrOM_col(NE,NY,NX)=tLitrOM_col(NE,NY,NX)+litrOM(NE)
  ENDDO

  WS=CanH2OHeldVg(NY,NX)+CanWat_col(NY,NX)+VLWatMicP_vr(0,NY,NX)+VLiceMicP(0,NY,NX)*DENSICE
  WaterStoreLandscape=WaterStoreLandscape+WS
  UVLWatMicP(NY,NX)=UVLWatMicP(NY,NX)+WS
  HeatStoreLandscape=HeatStoreLandscape+TENGYC(NY,NX)
  CS=trc_solml_vr(idg_CO2,0,NY,NX)+trc_solml_vr(idg_CH4,0,NY,NX)
  TGasC_lnd=TGasC_lnd+CS
  DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)+CS
  HS=trc_solml_vr(idg_H2,0,NY,NX)
  TLH2G=TLH2G+HS
  OS=trc_solml_vr(idg_O2,0,NY,NX)
  OXYGSO=OXYGSO+OS
  ZG=trc_solml_vr(idg_N2,0,NY,NX)+trc_solml_vr(idg_N2O,0,NY,NX)
  TGasN_lnd=TGasN_lnd+ZG
  Z4S=trc_solml_vr(ids_NH4,0,NY,NX)+trc_solml_vr(idg_NH3,0,NY,NX)
  Z4X=natomw*trcx_solml(idx_NH4,0,NY,NX)
  Z4F=natomw*(FertN_soil_vr(ifert_nh4,0,NY,NX)+FertN_soil_vr(ifert_urea,0,NY,NX) &
    +FertN_soil_vr(ifert_nh3,0,NY,NX))
  TDisolNH4_lnd=TDisolNH4_lnd+Z4S+Z4X+Z4F
  tNH4_col(NY,NX)=tNH4_col(NY,NX)+Z4S+Z4X

  ZOS=trc_solml_vr(ids_NO3,0,NY,NX)+trc_solml_vr(ids_NO2,0,NY,NX)
  ZOF=natomw*FertN_soil_vr(ifert_no3,0,NY,NX)
  tNO3_lnd=tNO3_lnd+ZOS+ZOF
  tNO3_col(NY,NX)=tNO3_col(NY,NX)+ZOS

  POS=trc_solml_vr(ids_H1PO4,0,NY,NX)+trc_solml_vr(ids_H2PO4,0,NY,NX)
  POX=patomw*(trcx_solml(idx_HPO4,0,NY,NX)+trcx_solml(idx_H2PO4,0,NY,NX))
  POP=patomw*(trcp_salml(idsp_AlPO4,0,NY,NX)+trcp_salml(idsp_FePO4,0,NY,NX) &
    +trcp_salml(idsp_CaHPO4,0,NY,NX))+2._r8*patomw*trcp_salml(idsp_CaH4P2O8,0,NY,NX) &
    +3._r8*patomw*trcp_salml(idsp_HA,0,NY,NX)
  TDisolPi_lnd=TDisolPi_lnd+POS+POX+POP
  tHxPO4_col(NY,NX)=tHxPO4_col(NY,NX)+POX
  tXPO4_col(NY,NX)=tXPO4_col(NY,NX)+POP

  IF(salt_model)call DiagSurfLitRLayerSalt(NY,NX,TDisolPi_lnd)

  end subroutine LitterLayerSummary
!------------------------------------------------------------------------------------------

  subroutine DiagSurfLitRLayerSalt(NY,NX,TDisolPi_lnd)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(inout) :: TDisolPi_lnd
  real(r8) :: SSS,PSS
  INTEGER :: nsalts

  DO nsalts=idsalt_beg,idsalt_end
    trcSalt_solml(nsalts,0,NY,NX)=trcSalt_solml(nsalts,0,NY,NX)+trcSalt3DFlo2Cell(nsalts,3,0,NY,NX)
  ENDDO

  PSS=patomw*(trcSalt_solml(idsalt_H0PO4,0,NY,NX)+trcSalt_solml(idsalt_H3PO4,0,NY,NX)&
    +trcSalt_solml(idsalt_FeHPO4,0,NY,NX) &
    +trcSalt_solml(idsalt_FeH2PO4,0,NY,NX)+trcSalt_solml(idsalt_CaPO4,0,NY,NX)&
    +trcSalt_solml(idsalt_CaHPO4,0,NY,NX) &
    +trcSalt_solml(idsalt_CaH4P2O8,0,NY,NX)+trcSalt_solml(idsalt_MgHPO4,0,NY,NX))
  TDisolPi_lnd=TDisolPi_lnd+PSS
  SSS=0._r8
  do nsalts=idsalt_beg,idsalt_end
    SSS=SSS+trcSalt_solml(nsalts,0,NY,NX)*trcSaltIonNumber(nsalts)
  enddo

  TION=TION+SSS
  UION(NY,NX)=UION(NY,NX)+SSS
  end subroutine DiagSurfLitRLayerSalt
!------------------------------------------------------------------------------------------
  subroutine update_physVar_Profile(NY,NX,VOLISO,DVLiceMicP)    
  !     WATER, ICE, HEAT, TEMPERATUR
  !
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(inout) :: VOLISO  
  REAL(R8),INTENT(OUT) :: DVLiceMicP(JZ)  !change in ice volume
  real(r8) :: TKS00,TKSX,TKSpre
  real(r8) :: ENGY
  real(r8) :: TVHeatCapacity
  real(r8) :: TVHeatCapacitySoilM,TVOLW,TVOLWH,TVOLI,TVOLIH,TENGY
  real(r8) :: VOLWXX,VOLIXX,VHeatCapacityX,WS
  real(r8) :: DVLWatMicP_vr(JZ,JY,JX)   !change in water volume
  integer :: L

  TVHeatCapacity=0.0_r8
  TVHeatCapacitySoilM=0.0_r8
  TVOLW=0.0_r8
  TVOLWH=0.0_r8
  TVOLI=0.0_r8
  TVOLIH=0.0_r8
  TENGY=0.0_r8
  DVLiceMicP=0._r8
  
  DO L=NU(NY,NX),NL(NY,NX)
    TKSX=TKS(L,NY,NX)
    VHeatCapacityX=VHeatCapacity(L,NY,NX)
    VOLWXX=VLWatMicP_vr(L,NY,NX)
    VOLIXX=VLiceMicP(L,NY,NX)
    !micropore
    VLWatMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)+TWatFlowCellMicP(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
      +WatIceThawMicP(L,NY,NX)+GridPlantRootH2OUptake_vr(L,NY,NX)+FWatIrrigate2MicP(L,NY,NX)
    VLWatMicPX(L,NY,NX)=VLWatMicPX(L,NY,NX)+TWatFlowCellMicPX(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
      +WatIceThawMicP(L,NY,NX)+GridPlantRootH2OUptake_vr(L,NY,NX)+FWatIrrigate2MicP(L,NY,NX)

    !do a numerical correction
    VLWatMicPX(L,NY,NX)=AMIN1(VLWatMicP_vr(L,NY,NX),VLWatMicPX(L,NY,NX)+0.01_r8*(VLWatMicP_vr(L,NY,NX)-VLWatMicPX(L,NY,NX)))
    VLiceMicP(L,NY,NX)=VLiceMicP(L,NY,NX)-WatIceThawMicP(L,NY,NX)/DENSICE

    !micropore
    VLWatMacP(L,NY,NX)=VLWatMacP(L,NY,NX)+TWaterFlowMacP(L,NY,NX)-FWatExMacP2MicP(L,NY,NX)+WatIceThawMacP(L,NY,NX)
    VLiceMacP(L,NY,NX)=VLiceMacP(L,NY,NX)-WatIceThawMacP(L,NY,NX)/DENSICE

    !volume change
    DVLWatMicP_vr(L,NY,NX)=VLWatMicP1(L,NY,NX)+VLWatMacP1(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLWatMacP(L,NY,NX)
    DVLiceMicP(L)=VLiceMicP1(L,NY,NX)+VLiceMacP1(L,NY,NX)-VLiceMicP(L,NY,NX)-VLiceMacP(L,NY,NX)

    !update water/ice-unfilled pores
    IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
      VLsoiAirP(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP(L,NY,NX) &
        +VLMacP(L,NY,NX)-VLWatMacP(L,NY,NX)-VLiceMacP(L,NY,NX))
    ELSE
      VLsoiAirP(L,NY,NX)=0.0_r8
!     VLMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)+VLiceMicP(L,NY,NX)
!    2+DVLWatMicP_vr(L,NY,NX)+DVLiceMicP(L,NY,NX)
!     VLSoilPoreMicP_vr(L,NY,NX)=VLMicP_vr(L,NY,NX)
!     VGeomLayer(L,NY,NX)=VLMicP_vr(L,NY,NX)
    ENDIF
    ENGY=VHeatCapacityX*TKSX
    VHeatCapacity(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP(L,NY,NX)) &
      +cpi*(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX))

    TVHeatCapacity=TVHeatCapacity+VHeatCapacity(L,NY,NX)
    TVHeatCapacitySoilM=TVHeatCapacitySoilM+VHeatCapacitySoilM(L,NY,NX)    
    TVOLW=TVOLW+VLWatMicP_vr(L,NY,NX)
    TVOLWH=TVOLWH+VLWatMacP(L,NY,NX)
    TVOLI=TVOLI+VLiceMicP(L,NY,NX)
    TVOLIH=TVOLIH+VLiceMacP(L,NY,NX)
    TENGY=TENGY+ENGY
    !
    !     ARTIFICIAL SOIL WARMING
    !
    !     IF(NX.EQ.3.AND.NY.EQ.2.AND.L.GT.NU(NY,NX)
    !    3.AND.L.LE.17.AND.I.GE.152.AND.I.LE.304)THEN
    !     THeatFlow2Soil(L,NY,NX)=THeatFlow2Soil(L,NY,NX)
    !    2+(TKSZ(I,J,L)-TKS(L,NY,NX))*VHeatCapacity(L,NY,NX)
    !     WRITE(*,3379)'TKSZ',I,J,NX,NY,L,TKSZ(I,J,L)
    !    2,TKS(L,NY,NX),VHeatCapacity(L,NY,NX),THeatFlow2Soil(L,NY,NX)
    !3379  FORMAT(A8,6I4,12E12.4)
    !     ENDIF
    !
    !     END ARTIFICIAL SOIL WARMING
    !
    TKSpre=TKS(L,NY,NX)
    IF(VHeatCapacity(L,NY,NX).GT.ZEROS(NY,NX) &
      .and.VHeatCapacity(L,NY,NX)/(VHeatCapacityX+VHeatCapacity(L,NY,NX))>0.05_r8)THEN
      TKS00=TKS(L,NY,NX)
      TKS(L,NY,NX)=(ENGY+THeatFlow2Soil(L,NY,NX)+THeatSoiThaw(L,NY,NX) &
        +THeatRootUptake(L,NY,NX)+HeatIrrigation(L,NY,NX))/VHeatCapacity(L,NY,NX)

!      if(curday>=285.and.L<=2)write(*,*)'rexL=',L,NY,NX,curhour,VHeatCapacityX,VHeatCapacity(L,NY,NX),&
!        SoiBulkDensity(L,NY,NX),NU(NY,NX)
      if(TKS(L,NY,NX)>4.e2)then
        write(*,*)'weird temperature in redist',L,NY,NX,TKSpre,TKS(L,NY,NX)
        write(*,*)'heatcap',VHeatCapacityX,VHeatCapacity(L,NY,NX),ZEROS(NY,NX),SoiBulkDensity(L,NY,NX)
        write(*,*)'itemized',ENGY/VHeatCapacity(L,NY,NX),&
          THeatFlow2Soil(L,NY,NX)/VHeatCapacity(L,NY,NX),&
          THeatSoiThaw(L,NY,NX)/VHeatCapacity(L,NY,NX), &
          THeatRootUptake(L,NY,NX)/VHeatCapacity(L,NY,NX),&
          HeatIrrigation(L,NY,NX)/VHeatCapacity(L,NY,NX)
        write(*,*)'wat',VLWatMicP_vr(L,NY,NX),VLWatMacP(L,NY,NX), &
          VLiceMicP(L,NY,NX),VLiceMacP(L,NY,NX)  
        write(*,*)'heat',ENGY,THeatFlow2Soil(L,NY,NX),VHeatCapacitySoilM(L,NY,NX)  
        call endrun()  
      endif

    ELSE
      TKS(L,NY,NX)=TKS(NUM(NY,NX),NY,NX)
    ENDIF
    TCS(L,NY,NX)=units%Kelvin2Celcius(TKS(L,NY,NX))
    WS=VLWatMicP_vr(L,NY,NX)+VLWatMacP(L,NY,NX)+(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX))*DENSICE
    WaterStoreLandscape=WaterStoreLandscape+WS
    VOLISO=VOLISO+VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX)
    UVLWatMicP(NY,NX)=UVLWatMicP(NY,NX)+WS
!    2-WiltPoint(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
    HeatStoreLandscape=HeatStoreLandscape+VHeatCapacity(L,NY,NX)*TKS(L,NY,NX)
  ENDDO
  end subroutine update_physVar_Profile
!------------------------------------------------------------------------------------------
  subroutine UpdateChemInSoilLays(NY,NX,LG,DORGC,TXCO2,DORGE)
  !

  implicit none
  integer,  intent(in) :: NY,NX,LG
  real(r8), intent(in) :: DORGE(JY,JX)
  real(r8), intent(inout) :: TXCO2(JY,JX)
  real(r8), intent(out) :: DORGC(JZ)

  integer  :: L,K,M,N,LL,NGL,NTX,NTP,NTG,NTS
  real(r8) :: HS,CS
  real(r8) :: CIB,CHB,OIB,COB
  real(r8) :: HGB,HOB,OS,OOB
  real(r8) :: POS,POX,POP
  real(r8) :: SNM,SPM,SSB,SD

  real(r8) :: WX,ZG,Z4S,Z4X,Z4F,ZOS,ZOF
  real(r8) :: ZGB,Z2B,ZHB
  integer  :: idg,idom,ids,NE

  !     begin_execution
  !     UPDATE SOIL LAYER VARIABLES WITH TOTAL FLUXES

  DORGC=0._r8
  D125: DO L=NU(NY,NX),NL(NY,NX)
    !

    SurfGasFlx_col(idg_H2,NY,NX)=SurfGasFlx_col(idg_H2,NY,NX)+Micb_N2Fixation_vr(L,NY,NX)

    !
    !     RESIDUE FROM PLANT LitrFall
!
    D8565: DO K=1,micpar%NumOfPlantLitrCmplxs
      DO  M=1,jsken
        SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)+LitrfalStrutElms_vr(ielmc,M,K,L,NY,NX)*micpar%OMCI(1,K)
        DO NE=1,NumPlantChemElms
          SolidOM_vr(NE,M,K,L,NY,NX)=SolidOM_vr(NE,M,K,L,NY,NX)+LitrfalStrutElms_vr(NE,M,K,L,NY,NX)
        ENDDO
      enddo
    ENDDO D8565
!
    !     DOC, DON, DOP FROM AQUEOUS TRANSPORT
!
    D8560: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_vr(idom,K,L,NY,NX)=DOM_vr(idom,K,L,NY,NX)+DOM_Transp2Micp_flx(idom,K,L,NY,NX) &
          +DOM_PoreTranspFlx(idom,K,L,NY,NX)
          
        DOM_MacP_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)+DOM_Transp2Macp_flx(idom,K,L,NY,NX) &
          -DOM_PoreTranspFlx(idom,K,L,NY,NX)

      enddo
    ENDDO D8560
    !
    !     DOC, DON, DOP FROM PLANT EXUDATION
    !
    D195: DO K=1,jcplx
      DO NE=1,NumPlantChemElms
        DOM_vr(NE,K,L,NY,NX)=DOM_vr(NE,K,L,NY,NX)+tRootMycoExud2Soil_vr(NE,K,L,NY,NX)
      ENDDO
    ENDDO D195
    !
    !     SOIL SOLUTES FROM AQUEOUS TRANSPORT, MICROBIAL AND ROOT
    !     EXCHANGE, EQUILIBRIUM REACTIONS, GAS EXCHANGE,
    !     MICROPORE-MACROPORE EXCHANGE,
    !
    trc_solml_vr(idg_N2,L,NY,NX)=trc_solml_vr(idg_N2,L,NY,NX) &
      -Micb_N2Fixation_vr(L,NY,NX) 

    trc_solml_vr(idg_CO2,L,NY,NX)=trc_solml_vr(idg_CO2,L,NY,NX) &
      +TR_CO2_aqu_soil_vr(L,NY,NX) 

    do idg=idg_beg,idg_NH3-1
    trc_solml_vr(idg,L,NY,NX)=trc_solml_vr(idg,L,NY,NX) &
      +trcs_Transp2MicP_vr(idg,L,NY,NX)+Gas_Disol_Flx_vr(idg,L,NY,NX) &
      -trcg_RMicbTransf_vr(idg,L,NY,NX)-trcs_plant_uptake_vr(idg,L,NY,NX) &
      +trcs_Irrig_vr(idg,L,NY,NX)+trcs_PoreTranspFlx_vr(idg,L,NY,NX) &
      +trcg_ebu_flx_vr(idg,L,NY,NX)
    enddo

    trc_solml_vr(idg_NH3,L,NY,NX)=trc_solml_vr(idg_NH3,L,NY,NX) &
      +TR_NH3_soil_vr(L,NY,NX) &
      +trcs_Transp2MicP_vr(idg_NH3,L,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,L,NY,NX) &
      -trcs_plant_uptake_vr(idg_NH3,L,NY,NX)&
      +trcs_Irrig_vr(idg_NH3,L,NY,NX) &
      +trcs_PoreTranspFlx_vr(idg_NH3,L,NY,NX)+trcg_ebu_flx_vr(idg_NH3,L,NY,NX)

    trc_solml_vr(idg_NH3B,L,NY,NX)=trc_solml_vr(idg_NH3B,L,NY,NX) &
      +Gas_Disol_Flx_vr(idg_NH3B,L,NY,NX)+trcg_ebu_flx_vr(idg_NH3B,L,NY,NX) &
      +trcs_Transp2MicP_vr(idg_NH3B,L,NY,NX) &
      +trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX)-trcs_plant_uptake_vr(idg_NH3B,L,NY,NX) &
      +trcs_Irrig_vr(idg_NH3B,L,NY,NX)+trcs_PoreTranspFlx_vr(idg_NH3B,L,NY,NX)

   DO ids=ids_nut_beg,ids_nuts_end
      trc_solml_vr(ids,L,NY,NX)=trc_solml_vr(ids,L,NY,NX) &
        +trcs_Transp2MicP_vr(ids,L,NY,NX)+RNutMicbTransf_vr(ids,L,NY,NX) &
        +trcn_RChem_soil_vr(ids,L,NY,NX)-trcs_plant_uptake_vr(ids,L,NY,NX) &
        +trcs_Irrig_vr(ids,L,NY,NX)+trcs_PoreTranspFlx_vr(ids,L,NY,NX)
    ENDDO
    
    do ids=ids_NH4B,ids_nutb_end
      trc_solml_vr(ids,L,NY,NX)=trc_solml_vr(ids,L,NY,NX) &
        +trcs_Transp2MicP_vr(ids,L,NY,NX)+RNutMicbTransf_vr(ids,L,NY,NX) &
        +trcn_RChem_band_soil_vr(ids,L,NY,NX)-trcs_plant_uptake_vr(ids,L,NY,NX) &
        +trcs_Irrig_vr(ids,L,NY,NX)+trcs_PoreTranspFlx_vr(ids,L,NY,NX)
    enddo


    Eco_HR_col(NY,NX)=Eco_HR_col(NY,NX)+trcg_RMicbTransf_vr(idg_CO2,L,NY,NX) &
      +trcg_RMicbTransf_vr(idg_CH4,L,NY,NX)
    SurfGasFlx_col(idg_N2,NY,NX)=SurfGasFlx_col(idg_N2,NY,NX)+trcg_RMicbTransf_vr(idg_N2,L,NY,NX)
    !
    !     EXCHANGEABLE CATIONS AND ANIONS FROM EXCHANGE REACTIONS
    !

    trcx_solml(idx_NH4,L,NY,NX)=trcx_solml(idx_NH4,L,NY,NX)+trcx_TRSoilChem_vr(idx_NH4,L,NY,NX)
    trcx_solml(idx_NH4B,L,NY,NX)=trcx_solml(idx_NH4B,L,NY,NX)+trcx_TRSoilChem_vr(idx_NH4B,L,NY,NX)

    DO NTX=idx_AEC+1,idx_end
      trcx_solml(NTX,L,NY,NX)=trcx_solml(NTX,L,NY,NX)+trcx_TRSoilChem_vr(NTX,L,NY,NX)
    ENDDO

    !
    !     PRECIPITATES FROM PRECIPITATION-DISSOLUTION REACTIONS
    !
    DO NTP=idsp_p_beg,idsp_p_end
      trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)+trcp_RChem_soil(NTP,L,NY,NX)
    ENDDO
    !
    !     MACROPORE SOLUTES FROM MACROPORE-MICROPORE EXCHANGE
    !
    DO NTS=ids_beg,ids_end
      trc_soHml_vr(NTS,L,NY,NX)=trc_soHml_vr(NTS,L,NY,NX)+trcs_Transp2MacP_vr(NTS,L,NY,NX) &
        -trcs_PoreTranspFlx_vr(NTS,L,NY,NX)
    ENDDO
    !
    !     GASES FROM VOLATILIZATION-DISSOLUTION AND GAS TRANSFER
!
    DO NTG=idg_beg,idg_NH3
      trc_gasml_vr(NTG,L,NY,NX)=trc_gasml_vr(NTG,L,NY,NX)+Gas_AdvDif_Flx_vr(NTG,L,NY,NX) &
        -Gas_Disol_Flx_vr(NTG,L,NY,NX)
    ENDDO

    trc_gasml_vr(idg_NH3,L,NY,NX)=trc_gasml_vr(idg_NH3,L,NY,NX)-Gas_Disol_Flx_vr(idg_NH3B,L,NY,NX) &
      +TRN3G(L,NY,NX)
    RO2GasXchangePrev_vr(L,NY,NX)=Gas_AdvDif_Flx_vr(idg_O2,L,NY,NX)
    RCO2GasFlxPrev_vr(L,NY,NX)=Gas_AdvDif_Flx_vr(idg_CO2,L,NY,NX)
    RCH4F(L,NY,NX)=Gas_AdvDif_Flx_vr(idg_CH4,L,NY,NX)
    RO2AquaXchangePrev_vr(L,NY,NX)=trcs_Transp2MicP_vr(idg_O2,L,NY,NX)+trcs_Irrig_vr(idg_O2,L,NY,NX) &
      +trcs_PoreTranspFlx_vr(idg_O2,L,NY,NX)+trcg_ebu_flx_vr(idg_O2,L,NY,NX)
    RCH4PhysexchPrev_vr(L,NY,NX)=trcs_Transp2MicP_vr(idg_CH4,L,NY,NX)+trcs_Irrig_vr(idg_CH4,L,NY,NX) &
      +trcs_PoreTranspFlx_vr(idg_CH4,L,NY,NX)+trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
    !
    !     GRID CELL BOUNDARY FLUXES FROM ROOT GAS TRANSFER
!   watch out the following code for changes
    HEATIN=HEATIN+THeatSoiThaw(L,NY,NX)+THeatRootUptake(L,NY,NX)
    CIB=trcg_air2root_flx_vr(idg_CO2,L,NY,NX)
    CHB=trcg_air2root_flx_vr(idg_CH4,L,NY,NX)
    OIB=trcg_air2root_flx_vr(idg_O2,L,NY,NX)
    HGB=trcg_air2root_flx_vr(idg_H2,L,NY,NX)
    ZGB=trcg_air2root_flx_vr(idg_N2,L,NY,NX)
    Z2B=trcg_air2root_flx_vr(idg_N2O,L,NY,NX)
    ZHB=trcg_air2root_flx_vr(idg_NH3,L,NY,NX)
!
!     GRID CELL BOUNDARY FLUXES BUBBLING
!
    IF(LG.EQ.0)THEN
      LL=0
      CIB=CIB+trcg_ebu_flx_vr(idg_CO2,L,NY,NX)
      CHB=CHB+trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
      OIB=OIB+trcg_ebu_flx_vr(idg_O2,L,NY,NX)
      ZGB=ZGB+trcg_ebu_flx_vr(idg_N2,L,NY,NX)
      Z2B=Z2B+trcg_ebu_flx_vr(idg_N2O,L,NY,NX)
      ZHB=ZHB+trcg_ebu_flx_vr(idg_NH3,L,NY,NX)+trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)
      HGB=HGB+trcg_ebu_flx_vr(idg_H2,L,NY,NX)
    ELSE
      LL=MIN(L,LG)
      DO NTG=idg_beg,idg_NH3
        trc_gasml_vr(NTG,LL,NY,NX)=trc_gasml_vr(NTG,LL,NY,NX)-trcg_ebu_flx_vr(NTG,L,NY,NX)
      ENDDO
      trc_gasml_vr(idg_NH3,LL,NY,NX)=trc_gasml_vr(idg_NH3,LL,NY,NX)-trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)

      IF(LG.LT.L)THEN
        TGasC_lnd=TGasC_lnd-trcg_ebu_flx_vr(idg_CO2,L,NY,NX)-trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
        DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)-trcg_ebu_flx_vr(idg_CO2,L,NY,NX)-trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
        OXYGSO=OXYGSO-trcg_ebu_flx_vr(idg_O2,L,NY,NX)
        TGasN_lnd=TGasN_lnd-trcg_ebu_flx_vr(idg_N2,L,NY,NX)-trcg_ebu_flx_vr(idg_N2O,L,NY,NX) &
          -trcg_ebu_flx_vr(idg_NH3,L,NY,NX)-trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)
        TLH2G=TLH2G-trcg_ebu_flx_vr(idg_H2,L,NY,NX)
      ENDIF
    ENDIF
    CO2GIN=CO2GIN+CIB+CHB
    COB=tRootCO2Emis_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_CO2,L,NY,NX)-TR_CO2_aqu_soil_vr(L,NY,NX)
    TOMOU(ielmc)=TOMOU(ielmc)+COB
    SurfGasFlx_col(idg_CO2,NY,NX)=SurfGasFlx_col(idg_CO2,NY,NX)+CIB
    SurfGasFlx_col(idg_CH4,NY,NX)=SurfGasFlx_col(idg_CH4,NY,NX)+CHB
    UCOP(NY,NX)=UCOP(NY,NX)+tRootCO2Emis_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_CO2,L,NY,NX)
    HydroSubsDICFlx_col(NY,NX)=HydroSubsDICFlx_col(NY,NX)-12.0*TBCO2(L,NY,NX)
    TXCO2(NY,NX)=TXCO2(NY,NX)+12.0*TBCO2(L,NY,NX)
    OXYGIN=OXYGIN+OIB
    OOB=trcg_RMicbTransf_vr(idg_O2,L,NY,NX)+tRO2MicrbUptk_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_O2,L,NY,NX)
    OXYGOU=OXYGOU+OOB
    H2GIN=H2GIN+HGB
    HOB=trcg_RMicbTransf_vr(idg_H2,L,NY,NX)+trcs_plant_uptake_vr(idg_H2,L,NY,NX)
    H2GOU=H2GOU+HOB
    ZN2GIN=ZN2GIN+ZGB+Z2B+ZHB
    SurfGasFlx_col(idg_O2,NY,NX)=SurfGasFlx_col(idg_O2,NY,NX)+OIB
    SurfGasFlx_col(idg_N2,NY,NX)=SurfGasFlx_col(idg_N2,NY,NX)+ZGB
    SurfGasFlx_col(idg_N2O,NY,NX)=SurfGasFlx_col(idg_N2O,NY,NX)+Z2B
    SurfGasFlx_col(idg_NH3,NY,NX)=SurfGasFlx_col(idg_NH3,NY,NX)+ZHB
    SurfGasFlx_col(idg_H2,NY,NX)=SurfGasFlx_col(idg_H2,NY,NX)+HGB
    !
    !     GRID CELL BOUNDARY FLUXES FROM EQUILIBRIUM REACTIONS
!
    SNM=(2.0_r8*(RNutMicbTransf_vr(ids_NH4,L,NY,NX)+RNutMicbTransf_vr(ids_NH4B,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_NH4,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_NH4B,L,NY,NX)-Micb_N2Fixation_vr(L,NY,NX))&
      -trcs_plant_uptake_vr(idg_NH3,L,NY,NX) &
      -trcs_plant_uptake_vr(idg_NH3B,L,NY,NX) &
      +RNutMicbTransf_vr(ids_NO3,L,NY,NX)+RNutMicbTransf_vr(ids_NO3B,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_NO3,L,NY,NX)-trcs_plant_uptake_vr(ids_NO3B,L,NY,NX) &
      +RNutMicbTransf_vr(ids_NO2,L,NY,NX)+RNutMicbTransf_vr(ids_NO2B,L,NY,NX))/natomw
    SPM=(2.0_r8*(RNutMicbTransf_vr(ids_H1PO4,L,NY,NX)+RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_H1PO4,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_H1PO4B,L,NY,NX))+3.0*(RNutMicbTransf_vr(ids_H2PO4,L,NY,NX) &
      +RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_H2PO4,L,NY,NX)-trcs_plant_uptake_vr(ids_H2PO4B,L,NY,NX)))/patomw
    SSB=TRH2O(L,NY,NX)+TBCO2(L,NY,NX)+XZHYS(L,NY,NX)+TBION(L,NY,NX)
    TIONOU=TIONOU-SSB
    !     HydroIonFlx_col(NY,NX)=HydroIonFlx_col(NY,NX)-SSB
    !     WRITE(20,3339)'SSB',I,J,L,SSB,TRH2O(L,NY,NX)
    !    2,TBCO2(L,NY,NX),XZHYS(L,NY,NX),TBION(L,NY,NX)
    !
    !     GAS AND SOLUTE EXCHANGE WITHIN GRID CELL ADDED TO ECOSYSTEM
    !     TOTALS FOR CALCULATING COMPETITION CONSTRAINTS ON MICROBIAL
    !     AND ROOT POPULATIONS
!
    call SumMicBGCFluxes(L,NY,NX)

    !
    !     GRID CELL VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
    !     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
    !     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS SINCE START OF RUN
    !
    !     IF(J.EQ.24)THEN
    SD=SAND(L,NY,NX)+SILT(L,NY,NX)+CLAY(L,NY,NX)
    TSEDSO=TSEDSO+SD
    CS=trc_gasml_vr(idg_CO2,L,NY,NX)+trc_solml_vr(idg_CO2,L,NY,NX) &
      +trc_soHml_vr(idg_CO2,L,NY,NX)+trcg_root_vr(idg_CO2,L,NY,NX) &
      +trc_gasml_vr(idg_CH4,L,NY,NX)+trc_solml_vr(idg_CH4,L,NY,NX) &
      +trc_soHml_vr(idg_CH4,L,NY,NX)+trcg_root_vr(idg_CH4,L,NY,NX)
    TGasC_lnd=TGasC_lnd+CS
    DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)+CS
    HS=trc_gasml_vr(idg_H2,L,NY,NX)+trc_solml_vr(idg_H2,L,NY,NX) &
      +trc_soHml_vr(idg_H2,L,NY,NX)+trcg_root_vr(idg_H2,L,NY,NX)
    TLH2G=TLH2G+HS

    OS=trc_gasml_vr(idg_O2,L,NY,NX)+trc_solml_vr(idg_O2,L,NY,NX) &
      +trc_soHml_vr(idg_O2,L,NY,NX)+trcg_root_vr(idg_O2,L,NY,NX)
    OXYGSO=OXYGSO+OS
!  should I add N2 to ZG below? Dec, 18, 2022.
    ZG=trc_gasml_vr(idg_N2,L,NY,NX)+trc_solml_vr(idg_N2,L,NY,NX) &
      +trc_soHml_vr(idg_N2,L,NY,NX)+trcg_root_vr(idg_N2O,L,NY,NX) &
      +trc_gasml_vr(idg_N2O,L,NY,NX)+trc_solml_vr(idg_N2O,L,NY,NX) &
      +trc_soHml_vr(idg_N2O,L,NY,NX)+trcg_root_vr(idg_NH3,L,NY,NX) &
      +trc_gasml_vr(idg_NH3,L,NY,NX)
    TGasN_lnd=TGasN_lnd+ZG

    Z4S=trc_solml_vr(ids_NH4,L,NY,NX)+trc_soHml_vr(ids_NH4,L,NY,NX) &
      +trc_solml_vr(ids_NH4B,L,NY,NX)+trc_soHml_vr(ids_NH4B,L,NY,NX)&
      +trc_solml_vr(idg_NH3,L,NY,NX)+trc_soHml_vr(idg_NH3,L,NY,NX) &
      +trc_solml_vr(idg_NH3B,L,NY,NX)+trc_soHml_vr(idg_NH3B,L,NY,NX)

    Z4X=natomw*(trcx_solml(idx_NH4,L,NY,NX)+trcx_solml(idx_NH4B,L,NY,NX))
    Z4F=natomw*(FertN_soil_vr(ifert_nh4,L,NY,NX)+FertN_soil_vr(ifert_urea,L,NY,NX) &
      +FertN_soil_vr(ifert_nh3,L,NY,NX)+FertN_Band_vr(ifert_nh4_band,L,NY,NX) &
      +FertN_Band_vr(ifert_urea_band,L,NY,NX)+FertN_Band_vr(ifert_nh3_band,L,NY,NX))
    TDisolNH4_lnd=TDisolNH4_lnd+Z4S+Z4X+Z4F
    tNH4_col(NY,NX)=tNH4_col(NY,NX)+Z4S+Z4X

    ZOS=trc_solml_vr(ids_NO3,L,NY,NX)+trc_soHml_vr(ids_NO3,L,NY,NX) &
      +trc_solml_vr(ids_NO3B,L,NY,NX)+trc_soHml_vr(ids_NO3B,L,NY,NX) &
      +trc_solml_vr(ids_NO2,L,NY,NX)+trc_soHml_vr(ids_NO2,L,NY,NX) &
      +trc_solml_vr(ids_NO2B,L,NY,NX)+trc_soHml_vr(ids_NO2B,L,NY,NX)
    ZOF=natomw*(FertN_soil_vr(ifert_no3,L,NY,NX)+FertN_soil_vr(ifert_no3,L,NY,NX))
    tNO3_lnd=tNO3_lnd+ZOS+ZOF
    tNO3_col(NY,NX)=tNO3_col(NY,NX)+ZOS

    POS=trc_solml_vr(ids_H2PO4,L,NY,NX)+trc_soHml_vr(ids_H2PO4,L,NY,NX)+trc_solml_vr(ids_H2PO4B,L,NY,NX) &
      +trc_soHml_vr(ids_H2PO4B,L,NY,NX)+trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_soHml_vr(ids_H1PO4,L,NY,NX) &
      +trc_solml_vr(ids_H1PO4B,L,NY,NX)+trc_soHml_vr(ids_H1PO4B,L,NY,NX)

    !exchangeable P
    POX=patomw*(trcx_solml(idx_HPO4,L,NY,NX)+trcx_solml(idx_H2PO4,L,NY,NX) &
      +trcx_solml(idx_HPO4B,L,NY,NX)+trcx_solml(idx_H2PO4B,L,NY,NX))

    !precipitated P
    POP=patomw*(trcp_salml(idsp_AlPO4,L,NY,NX)+trcp_salml(idsp_FePO4,L,NY,NX)&
      +trcp_salml(idsp_CaHPO4,L,NY,NX)+trcp_salml(idsp_AlPO4B,L,NY,NX)&
      +trcp_salml(idsp_FePO4B,L,NY,NX)+trcp_salml(idsp_CaHPO4B,L,NY,NX)) &
      +2._r8*patomw*(trcp_salml(idsp_CaH4P2O8,L,NY,NX)+trcp_salml(idsp_CaH4P2O8B,L,NY,NX)) &
      +3._r8*patomw*(trcp_salml(idsp_HA,L,NY,NX)+trcp_salml(idsp_HAB,L,NY,NX))

    TDisolPi_lnd=TDisolPi_lnd+POS+POX+POP
    tHxPO4_col(NY,NX)=tHxPO4_col(NY,NX)+POX
    tXPO4_col(NY,NX)=tXPO4_col(NY,NX)+POP
    !
    call SumOMStates(L,NY,NX,DORGC(L),DORGE(NY,NX))
!
!     TOTAL SALT IONS
!
    IF(salt_model)call UpdateSaltIonInSoilLayers(L,NY,NX,TDisolPi_lnd)

  ENDDO D125
  end subroutine UpdateChemInSoilLays

!------------------------------------------------------------------------------------------
  subroutine UpdateSaltIonInSoilLayers(L,NY,NX,TDisolPi_lnd)
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(inout) :: TDisolPi_lnd
  integer  :: NTP,nsalts
  real(r8) :: ECHY,ECOH,ECAL,ECFE,ECCA,ECMG,ECNA,ECKA,ECCO,ECHC
  real(r8) :: ECNO,ECSO,ECCL
  real(r8) :: PSS,SSS,SSH,SSF,SSX,SST,SSP

  DO nsalts=idsalt_beg,idsaltb_end
    trcSalt_solml(nsalts,L,NY,NX)=trcSalt_solml(nsalts,L,NY,NX)+trcSalt_TR(nsalts,L,NY,NX) &
      +trcSalt_Flo2MicP_vr(nsalts,L,NY,NX)+trcSalt_RFLU(nsalts,L,NY,NX)+trcSalt_XFXS(nsalts,L,NY,NX)

    trcSalt_soHml(nsalts,L,NY,NX)=trcSalt_soHml(nsalts,L,NY,NX)+trcSalt_Flo2MacP_vr(nsalts,L,NY,NX) &
      -trcSalt_XFXS(nsalts,L,NY,NX)
  ENDDO
  trcSalt_solml(idsalt_AlOH2,L,NY,NX)=trcSalt_solml(idsalt_AlOH2,L,NY,NX)-TR_AlO2H2_sorbed_soil(L,NY,NX)
  trcSalt_solml(idsalt_FeOH2,L,NY,NX)=trcSalt_solml(idsalt_FeOH2,L,NY,NX)-TR_FeO2H2_sorbed_soil(L,NY,NX)

  trcx_solml(idx_Hp,L,NY,NX)=trcx_solml(idx_Hp,L,NY,NX)+TR_H_p_sorbed_soil(L,NY,NX)
  trcx_solml(idx_Al,L,NY,NX)=trcx_solml(idx_Al,L,NY,NX)+TR_Al_sorbed_soil(L,NY,NX)
  trcx_solml(idx_Fe,L,NY,NX)=trcx_solml(idx_Fe,L,NY,NX)+TR_Fe_sorbed_soil(L,NY,NX)
  trcx_solml(idx_Ca,L,NY,NX)=trcx_solml(idx_Ca,L,NY,NX)+TR_Ca_sorbed_soil(L,NY,NX)
  trcx_solml(idx_Mg,L,NY,NX)=trcx_solml(idx_Mg,L,NY,NX)+TR_Mg_sorbed_soil(L,NY,NX)
  trcx_solml(idx_Na,L,NY,NX)=trcx_solml(idx_Na,L,NY,NX)+TR_Na_sorbed_soil(L,NY,NX)
  trcx_solml(idx_K,L,NY,NX)=trcx_solml(idx_K,L,NY,NX)+TR_K_sorbed_soil(L,NY,NX)
  trcx_solml(idx_COOH,L,NY,NX)=trcx_solml(idx_COOH,L,NY,NX)+TR_HCO3_sorbed_soil(L,NY,NX)
  trcx_solml(idx_AlOH2,L,NY,NX)=trcx_solml(idx_AlOH2,L,NY,NX)+TR_AlO2H2_sorbed_soil(L,NY,NX)
  trcx_solml(idx_FeOH2,L,NY,NX)=trcx_solml(idx_FeOH2,L,NY,NX)+TR_FeO2H2_sorbed_soil(L,NY,NX)

! all non-P precipitates
  DO NTP=idsp_beg,idsp_p_beg-1
    trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)+trcp_RChem_soil(NTP,L,NY,NX)
  ENDDO

  PSS=0._r8
  DO nsalts=idsalt_psoil_beg,idsalt_pband_end
    PSS=PSS+trcSalt_solml(nsalts,L,NY,NX)+trcSalt_soHml(nsalts,L,NY,NX)
  ENDDO
  PSS=PSS*patomw  

  TDisolPi_lnd=TDisolPi_lnd+PSS

  SSS=0._r8
  SSH=0._r8
  DO nsalts=idsalt_beg,idsaltb_end
    SSS=SSS+trcSalt_solml(nsalts,L,NY,NX)*trcSaltIonNumber(nsalts)
    SSH=SSH+trcSalt_soHml(idsalt_Al,L,NY,NX)*trcSaltIonNumber(nsalts)
  ENDDO

 
!     TOTAL FERILIZER,EXCHANGEABLE CATIONS AND ANIONS, PRECIPITATES
!
  SSF=FertN_soil_vr(ifert_nh3,L,NY,NX)+FertN_soil_vr(ifert_urea,L,NY,NX) &
    +FertN_soil_vr(ifert_no3,L,NY,NX) &
    +FertN_Band_vr(ifert_nh3_band,L,NY,NX)+FertN_Band_vr(ifert_urea_band,L,NY,NX) &
    +FertN_Band_vr(ifert_no3_band,L,NY,NX) &
    +2.0_r8*(FertN_soil_vr(ifert_nh4,L,NY,NX)+FertN_Band_vr(ifert_nh4_band,L,NY,NX))

  SSX=trcx_solml(idx_Hp,L,NY,NX)+trcx_solml(idx_Al,L,NY,NX) &
    +trcx_solml(idx_Fe,L,NY,NX)+trcx_solml(idx_Ca,L,NY,NX)+trcx_solml(idx_Mg,L,NY,NX) &
    +trcx_solml(idx_Na,L,NY,NX)+trcx_solml(idx_K,L,NY,NX)+trcx_solml(idx_COOH,L,NY,NX) &
    +trcx_solml(idx_OHe,L,NY,NX)+trcx_solml(idx_OHeB,L,NY,NX) &
    +2.0_r8*(trcx_solml(idx_NH4,L,NY,NX)+trcx_solml(idx_NH4B,L,NY,NX) &
    +trcx_solml(idx_OH,L,NY,NX)+trcx_solml(idx_OHB,L,NY,NX)) &
    +3.0_r8*(trcx_solml(idx_AlOH2,L,NY,NX)+trcx_solml(idx_FeOH2,L,NY,NX) &
    +trcx_solml(idx_OHp,L,NY,NX)+trcx_solml(idx_OHpB,L,NY,NX) &
    +trcx_solml(idx_HPO4,L,NY,NX)+trcx_solml(idx_HPO4B,L,NY,NX)) &
    +4.0_r8*(trcx_solml(idx_H2PO4,L,NY,NX)+trcx_solml(idx_H2PO4B,L,NY,NX))

  SSP=2.0_r8*(trcp_salml(idsp_CaCO3,L,NY,NX)+trcp_salml(idsp_CaSO4,L,NY,NX) &
    +trcp_salml(idsp_AlPO4,L,NY,NX)+trcp_salml(idsp_FePO4,L,NY,NX) &
    +trcp_salml(idsp_AlPO4B,L,NY,NX)+trcp_salml(idsp_FePO4B,L,NY,NX)) &
    +3.0_r8*(trcp_salml(idsp_CaHPO4,L,NY,NX)+trcp_salml(idsp_CaHPO4B,L,NY,NX)) &
    +4.0_r8*(trcp_salml(idsp_AlOH3,L,NY,NX)+trcp_salml(idsp_FeOH3,L,NY,NX)) &
    +7.0_r8*(trcp_salml(idsp_CaH4P2O8,L,NY,NX)+trcp_salml(idsp_CaH4P2O8B,L,NY,NX)) &
    +9.0_r8*(trcp_salml(idsp_HA,L,NY,NX)+trcp_salml(idsp_HAB,L,NY,NX))

  SST=SSS+SSH+SSF+SSX+SSP
  TION=TION+SST
  UION(NY,NX)=UION(NY,NX)+SST

!
!     SOIL ELECTRICAL CONDUCTIVITY
!
  IF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    ECHY=0.337_r8*AZMAX1(trcSalt_solml(idsalt_Hp,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECOH=0.192_r8*AZMAX1(trcSalt_solml(idsalt_OH,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECAL=0.056_r8*AZMAX1(trcSalt_solml(idsalt_Al,L,NY,NX)*3.0_r8/VLWatMicP_vr(L,NY,NX))
    ECFE=0.051_r8*AZMAX1(trcSalt_solml(idsalt_Fe,L,NY,NX)*3.0_r8/VLWatMicP_vr(L,NY,NX))
    ECCA=0.060_r8*AZMAX1(trcSalt_solml(idsalt_Ca,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECMG=0.053_r8*AZMAX1(trcSalt_solml(idsalt_Mg,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECNA=0.050_r8*AZMAX1(trcSalt_solml(idsalt_Na,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECKA=0.070_r8*AZMAX1(trcSalt_solml(idsalt_K,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECCO=0.072_r8*AZMAX1(trcSalt_solml(idsalt_CO3,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECHC=0.044_r8*AZMAX1(trcSalt_solml(idsalt_HCO3,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECSO=0.080_r8*AZMAX1(trcSalt_solml(idsalt_SO4,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECCL=0.076_r8*AZMAX1(trcSalt_solml(idsalt_Cl,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECNO=0.071_r8*AZMAX1(trc_solml_vr(ids_NO3,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*natomw))
    ECND(L,NY,NX)=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA &
      +ECCO+ECHC+ECSO+ECCL+ECNO

  ELSE
    ECND(L,NY,NX)=0.0_r8
  ENDIF
  end subroutine UpdateSaltIonInSoilLayers

!------------------------------------------------------------------------------------------
  subroutine SumOMStates(L,NY,NX,DORGCL,DORGEC)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in):: DORGEC
  real(r8), intent(out):: DORGCL

  real(r8) :: ORGM(NumPlantChemElms)
  real(r8) :: litrOM(NumPlantChemElms)
  real(r8) :: HumOM(NumPlantChemElms)
  
  integer :: NE
    !     TOTAL SOC,SON,SOP
    !
    !     OMC=microbial biomass, ORC=microbial residue
    !     OQC,OQCH=DOC in micropores,macropores
    !     OQA,OQAH=acetate in micropores,macropores
    !     OHC,OHA=adsorbed SOC,acetate
    !     OSC=SOC(K=0:woody litter, K=1:non-woody litter,
    !     K=2:manure, K=3:POC, K=4:humus)
!
  !sum up living microbes
  call sumMicBiomLayL(L,NY,NX,ORGM)

  DO NE=1,NumPlantChemElms
    tMicBiome_col(NE,NY,NX)=tMicBiome_col(NE,NY,NX)+ORGM(NE)
  ENDDO
    
  call sumORGMLayL(L,NY,NX,SoilOrgM_vr(1:NumPlantChemElms,L,NY,NX))
  
  DO NE=1,NumPlantChemElms
    tSoilOrgM_col(NE,NY,NX)=tSoilOrgM_col(NE,NY,NX)+SoilOrgM_vr(NE,L,NY,NX)
  ENDDO

  !sum up litter complexes
  call sumLitrOMLayL(L,NY,NX,litrOM)
  OMLitrC_vr(L,NY,NX)=litrOM(ielmc)

  IF(iErosionMode.EQ.ieros_frzthawsom .OR. iErosionMode.EQ.ieros_frzthawsomeros)THEN
! change in organic C
    DORGCL=ORGCX_vr(L,NY,NX)-SoilOrgM_vr(ielmc,L,NY,NX)
    IF(L.EQ.NU(NY,NX))THEN
      DORGCL=DORGCL+DORGEC
    ENDIF
  ELSE
    DORGCL=0.0_r8
  ENDIF
  !sum non-litter complexes
  call sumHumOMLayL(L,NY,NX,HumOM)

  DO NE=1,NumPlantChemElms
    LitRMStoreLndscap(NE)=LitRMStoreLndscap(NE)+litrOM(NE)
    tLitrOM_col(NE,NY,NX)=tLitrOM_col(NE,NY,NX)+litrOM(NE)
    
    POMHumStoreLndscap(NE)=POMHumStoreLndscap(NE)+HumOM(NE)
    tHumOM_col(NE,NY,NX)=tHumOM_col(NE,NY,NX)+HumOM(NE)
  ENDDO
  TSEDSO=TSEDSO+(litrOM(ielmc)+HumOM(ielmc))*ppmc
  end subroutine SumOMStates

!------------------------------------------------------------------------------------------

  subroutine AddFlux2SurfaceResidue(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: M,N,K,NGL,NE
  real(r8) :: HRAINR,RAINR

!     begin_execution
!     ADD ABOVE-GROUND LitrFall FROM EXTRACT.F TO SURFACE RESIDUE
!
!     OSC,OSN,OSP=SOC,SON,SOP
!     CSNT,ZSNT,PSNT=total C,N,P LitrFall
!     ORGC=total SOC
!     RAINR,HRAINR=water,heat in LitrFall
!     FLWR,HFLWR=water,heat flux into litter
!     HEATIN=cumulative net surface heat transfer
!
  DO   K=1,micpar%NumOfPlantLitrCmplxs
    DO  M=1,jsken
      SolidOMAct_vr(M,K,0,NY,NX)=SolidOMAct_vr(M,K,0,NY,NX)+LitrfalStrutElms_vr(ielmc,M,K,0,NY,NX)*micpar%OMCI(1,K)
      DO NE=1,NumPlantChemElms
        SolidOM_vr(NE,M,K,0,NY,NX)=SolidOM_vr(NE,M,K,0,NY,NX)+LitrfalStrutElms_vr(NE,M,K,0,NY,NX)
      ENDDO

      SoilOrgM_vr(ielmc,0,NY,NX)=SoilOrgM_vr(ielmc,0,NY,NX)+LitrfalStrutElms_vr(ielmc,M,K,0,NY,NX)
      RAINR=LitrfalStrutElms_vr(ielmc,M,K,0,NY,NX)*ThetaCX(K)
      HRAINR=RAINR*cpw*TairK(NY,NX)
      WatFLo2Litr(NY,NX)=WatFLo2Litr(NY,NX)+RAINR
      HeatFLo2LitrByWat(NY,NX)=HeatFLo2LitrByWat(NY,NX)+HRAINR

      CRAIN=CRAIN+RAINR
      HEATIN=HEATIN+HRAINR
    enddo
  ENDDO

  call SumSurfMicBGCFluxes(NY,NX)

  end subroutine AddFlux2SurfaceResidue

!------------------------------------------------------------------------------------------

  subroutine SumSurfMicBGCFluxes(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  INTEGER :: K,N,NGL
  !
  !     GAS AND SOLUTE EXCHANGE WITHIN SURFACE LITTER ADDED TO ECOSYSTEM
  !     TOTALS FOR CALCULATING COMPETITION CONSTRAINTS ON MICROBIAL
  !     AND ROOT POPULATIONS IN NITRO.F AND UPTAKE.F
  !
  !     REcoO2DmndResp_vr=O2 demand by all microbial,root,myco populations
  !     REcoNH4DmndSoil_vr=NH4 demand in non-band by all microbial,root,myco populations
  !     REcoNO3DmndSoil_vr=NO3 demand in non-band by all microbial,root,myco populations
  !     REcoH2PO4DmndSoil_vr=H2PO4 demand in non-band by all microbial,root,myco populations
  !     REcoH1PO4DmndSoil_vr=HPO4 demand in non-band by all microbial,root,myco populations
  !     REcoNH4DmndBand_vr=NH4 demand in band by all microbial,root,myco populations
  !     REcoNO3DmndBand_vr=NO3 demand in band by all microbial,root,myco populations
  !     REcoH2PO4DmndBand_vr=H2PO4 demand in band by all microbial,root,myco populations
  !     REcoH1PO4DmndBand_vr=HPO4 demand in band by all microbial,root,myco populations
  !     RO2DmndHetert=O2 demand from DOC,DOA oxidation
  !     RNO3ReduxDmndSoilHeter_vr=demand for NO3 reduction
  !     RNO2DmndReduxSoilHeter_vr=demand for NO2 oxidation
  !     RN2ODmndReduxHeter_vr=demand for N2O reduction
  !     RNH4DmndSoilHeter_vr,RNH4DmndLitrHeter_col=substrate-unlimited NH4 immobilization
  !     RNO3DmndSoilHeter_vr,RNO3DmndLitrHeter_col=substrate-unlimited NO3 immobilization
  !     RH2PO4DmndSoilHeter_vr,RH2PO4DmndLitrHeter_col=substrate-unlimited H2PO4 immobilization
  !     RH1PO4DmndSoilHeter_vr,RH1PO4DmndLitrHeter_col=substrate-unlimited HPO4 immobilization
  !     ROQCX,ROQAX=total DOC,DOA demand from DOC,DOA oxidation
  !     RDOCUptkHeter_vr,RAcetateUptkHeter_vr=DOC,DOA demand from DOC,DOA oxidation
  !     RNO2EcoUptkSoil_vr=total demand for NO2 reduction
  !     RNO2DmndSoilChemo_vr=demand for NO2 reduction
!
  DO K=1,jcplx
    IF(micpar%is_litter(K))THEN
      DO  N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          REcoO2DmndResp_vr(0,NY,NX)=REcoO2DmndResp_vr(0,NY,NX)+RO2DmndHetert(NGL,K,0,NY,NX)
          REcoNO3DmndSoil_vr(0,NY,NX)=REcoNO3DmndSoil_vr(0,NY,NX)+RNO3ReduxDmndSoilHeter_vr(NGL,K,0,NY,NX)
          RNO2EcoUptkSoil_vr(0,NY,NX)=RNO2EcoUptkSoil_vr(0,NY,NX)+RNO2DmndReduxSoilHeter_vr(NGL,K,0,NY,NX)
          RN2OEcoUptkSoil_vr(0,NY,NX)=RN2OEcoUptkSoil_vr(0,NY,NX)+RN2ODmndReduxHeter_vr(NGL,K,0,NY,NX)
          REcoNH4DmndSoil_vr(0,NY,NX)=REcoNH4DmndSoil_vr(0,NY,NX)+RNH4DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoNO3DmndSoil_vr(0,NY,NX)=REcoNO3DmndSoil_vr(0,NY,NX)+RNO3DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoH2PO4DmndSoil_vr(0,NY,NX)=REcoH2PO4DmndSoil_vr(0,NY,NX)+RH2PO4DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoH1PO4DmndSoil_vr(0,NY,NX)=REcoH1PO4DmndSoil_vr(0,NY,NX)+RH1PO4DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoNH4DmndSoil_vr(NU(NY,NX),NY,NX)=REcoNH4DmndSoil_vr(NU(NY,NX),NY,NX)+RNH4DmndLitrHeter_col(NGL,K,NY,NX)
          REcoNO3DmndSoil_vr(NU(NY,NX),NY,NX)=REcoNO3DmndSoil_vr(NU(NY,NX),NY,NX)+RNO3DmndLitrHeter_col(NGL,K,NY,NX)
          REcoH2PO4DmndSoil_vr(NU(NY,NX),NY,NX)=REcoH2PO4DmndSoil_vr(NU(NY,NX),NY,NX)+RH2PO4DmndLitrHeter_col(NGL,K,NY,NX)
          REcoH1PO4DmndSoil_vr(NU(NY,NX),NY,NX)=REcoH1PO4DmndSoil_vr(NU(NY,NX),NY,NX)+RH1PO4DmndLitrHeter_col(NGL,K,NY,NX)
          ROQCX(K,0,NY,NX)=ROQCX(K,0,NY,NX)+RDOCUptkHeter_vr(NGL,K,0,NY,NX)
          ROQAX(K,0,NY,NX)=ROQAX(K,0,NY,NX)+RAcetateUptkHeter_vr(NGL,K,0,NY,NX)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      REcoO2DmndResp_vr(0,NY,NX)=REcoO2DmndResp_vr(0,NY,NX)+RO2DmndAutort(NGL,0,NY,NX)
      REcoNH4DmndSoil_vr(0,NY,NX)=REcoNH4DmndSoil_vr(0,NY,NX)+RNH3OxidAutor(NGL,0,NY,NX)
      RNO2EcoUptkSoil_vr(0,NY,NX)=RNO2EcoUptkSoil_vr(0,NY,NX)+RNO2OxidAutor(NGL,0,NY,NX)
      RN2OEcoUptkSoil_vr(0,NY,NX)=RN2OEcoUptkSoil_vr(0,NY,NX)+RN2ODmndReduxAutor_vr(NGL,0,NY,NX)
      REcoNH4DmndSoil_vr(0,NY,NX)=REcoNH4DmndSoil_vr(0,NY,NX)+RNH4UptkSoilAutor_vr(NGL,0,NY,NX)
      REcoNO3DmndSoil_vr(0,NY,NX)=REcoNO3DmndSoil_vr(0,NY,NX)+RNO3UptkSoilAutor_vr(NGL,0,NY,NX)
      REcoH2PO4DmndSoil_vr(0,NY,NX)=REcoH2PO4DmndSoil_vr(0,NY,NX)+RH2PO4UptkSoilAutor_vr(NGL,0,NY,NX)
      REcoH1PO4DmndSoil_vr(0,NY,NX)=REcoH1PO4DmndSoil_vr(0,NY,NX)+RH1PO4UptkSoilAutor_vr(NGL,0,NY,NX)
      REcoNH4DmndSoil_vr(NU(NY,NX),NY,NX)=REcoNH4DmndSoil_vr(NU(NY,NX),NY,NX)+RNH4UptkLitrAutor_col(NGL,NY,NX)
      REcoNO3DmndSoil_vr(NU(NY,NX),NY,NX)=REcoNO3DmndSoil_vr(NU(NY,NX),NY,NX)+RNO3UptkLitrAutor_col(NGL,NY,NX)
      REcoH2PO4DmndSoil_vr(NU(NY,NX),NY,NX)=REcoH2PO4DmndSoil_vr(NU(NY,NX),NY,NX)+RH2PO4UptkLitrAutor_col(NGL,NY,NX)
      REcoH1PO4DmndSoil_vr(NU(NY,NX),NY,NX)=REcoH1PO4DmndSoil_vr(NU(NY,NX),NY,NX)+RH1PO4UptkLitrAutor_col(NGL,NY,NX)

    ENDDO
  ENDDO

  RNO2EcoUptkSoil_vr(0,NY,NX)=RNO2EcoUptkSoil_vr(0,NY,NX)+RNO2DmndSoilChemo_vr(0,NY,NX)
  RNO2EcoUptkBand_vr(0,NY,NX)=RNO2EcoUptkBand_vr(0,NY,NX)+RNO2DmndBandChemo_vr(0,NY,NX)
  end subroutine SumSurfMicBGCFluxes

!------------------------------------------------------------------------------------------


  subroutine SumMicBGCFluxes(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX

  integer :: K,N,NGL

  REcoO2DmndResp_vr(L,NY,NX)=REcoO2DmndResp_vr(L,NY,NX)+SUM(RO2DmndHetert(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoNH4DmndSoil_vr(L,NY,NX)=REcoNH4DmndSoil_vr(L,NY,NX)+SUM(RNH4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoNO3DmndSoil_vr(L,NY,NX)=REcoNO3DmndSoil_vr(L,NY,NX)+SUM(RNO3ReduxDmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)) &
    +SUM(RNO3DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  RNO2EcoUptkSoil_vr(L,NY,NX)=RNO2EcoUptkSoil_vr(L,NY,NX)+SUM(RNO2DmndReduxSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  RN2OEcoUptkSoil_vr(L,NY,NX)=RN2OEcoUptkSoil_vr(L,NY,NX)+SUM(RN2ODmndReduxHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoH2PO4DmndSoil_vr(L,NY,NX)=REcoH2PO4DmndSoil_vr(L,NY,NX)+SUM(RH2PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoH1PO4DmndSoil_vr(L,NY,NX)=REcoH1PO4DmndSoil_vr(L,NY,NX)+SUM(RH1PO4DmndSoilHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoNH4DmndBand_vr(L,NY,NX)=REcoNH4DmndBand_vr(L,NY,NX)+SUM(RNH4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoNO3DmndBand_vr(L,NY,NX)=REcoNO3DmndBand_vr(L,NY,NX)+SUM(RNO3ReduxDmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX)) &
    +SUM(RNO3DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  RNO2EcoUptkBand_vr(L,NY,NX)=RNO2EcoUptkBand_vr(L,NY,NX)+SUM(RNO2DmndReduxBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoH2PO4DmndBand_vr(L,NY,NX)=REcoH2PO4DmndBand_vr(L,NY,NX)+SUM(RH2PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  REcoH1PO4DmndBand_vr(L,NY,NX)=REcoH1PO4DmndBand_vr(L,NY,NX)+SUM(RH1PO4DmndBandHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  DO K=1,jcplx
    ROQCX(K,L,NY,NX)=ROQCX(K,L,NY,NX)+SUM(RDOCUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
    ROQAX(K,L,NY,NX)=ROQAX(K,L,NY,NX)+SUM(RAcetateUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  ENDDO
  REcoO2DmndResp_vr(L,NY,NX)=REcoO2DmndResp_vr(L,NY,NX)+SUM(RO2DmndAutort(1:NumHetetrMicCmplx,L,NY,NX))
  REcoNH4DmndSoil_vr(L,NY,NX)=REcoNH4DmndSoil_vr(L,NY,NX)+SUM(RNH3OxidAutor(1:NumHetetrMicCmplx,L,NY,NX)) &
    +SUM(RNH4UptkSoilAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  REcoNO3DmndSoil_vr(L,NY,NX)=REcoNO3DmndSoil_vr(L,NY,NX)+SUM(RNO3UptkSoilAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  RNO2EcoUptkSoil_vr(L,NY,NX)=RNO2EcoUptkSoil_vr(L,NY,NX)+SUM(RNO2OxidAutor(1:NumHetetrMicCmplx,L,NY,NX))
  RN2OEcoUptkSoil_vr(L,NY,NX)=RN2OEcoUptkSoil_vr(L,NY,NX)+SUM(RN2ODmndReduxAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  REcoH2PO4DmndSoil_vr(L,NY,NX)=REcoH2PO4DmndSoil_vr(L,NY,NX)+SUM(RH2PO4UptkSoilAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  REcoH1PO4DmndSoil_vr(L,NY,NX)=REcoH1PO4DmndSoil_vr(L,NY,NX)+SUM(RH1PO4UptkSoilAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  REcoNH4DmndBand_vr(L,NY,NX)=REcoNH4DmndBand_vr(L,NY,NX)+SUM(RNH3OxidAutorBand(1:NumHetetrMicCmplx,L,NY,NX)) &
    +SUM(RNH4UptkBandAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  REcoNO3DmndBand_vr(L,NY,NX)=REcoNO3DmndBand_vr(L,NY,NX)+SUM(RNO3UptkBandAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  RNO2EcoUptkBand_vr(L,NY,NX)=RNO2EcoUptkBand_vr(L,NY,NX)+SUM(RNO2OxidAutorBand(1:NumHetetrMicCmplx,L,NY,NX))
  REcoH2PO4DmndBand_vr(L,NY,NX)=REcoH2PO4DmndBand_vr(L,NY,NX)+SUM(RH2PO4UptkBandAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))
  REcoH1PO4DmndBand_vr(L,NY,NX)=REcoH1PO4DmndBand_vr(L,NY,NX)+SUM(RH1PO4UptkBandAutor_vr(1:NumHetetrMicCmplx,L,NY,NX))

  RNO2EcoUptkSoil_vr(L,NY,NX)=RNO2EcoUptkSoil_vr(L,NY,NX)+RNO2DmndSoilChemo_vr(L,NY,NX)
  RNO2EcoUptkBand_vr(L,NY,NX)=RNO2EcoUptkBand_vr(L,NY,NX)+RNO2DmndBandChemo_vr(L,NY,NX)

  end subroutine SumMicBGCFluxes
end module RedistMod
