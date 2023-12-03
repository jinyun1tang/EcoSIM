module RedistMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : padr, print_info,endrun,destroy
  use minimathmod, only : safe_adb,AZMAX1
  use EcoSiMParDataMod, only : micpar
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
  use SurfLitterPhysMod, only : UpdateLitRPhys
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
  real(r8) :: DORGC(JZ,JY,JX),DVLiceMicP(JZ,JY,JX)
  real(r8) :: TXCO2(JY,JX),DORGE(JY,JX)
  real(r8) :: UDVLiceMicP,UDLYXF
  real(r8) :: VOLISO,VOLPT,VOLTT
  real(r8) :: TFLWT
!     execution begins here
  curday=I
  curhour=J
  VOLISO=0.0_r8
  UDVLiceMicP=0.0_r8
  UDLYXF=0.0_r8
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
      call SumLitterLayerChemMass(NY,NX)
!
      call update_physVar_Profile(NY,NX,VOLISO,DVLiceMicP)    

      call UpdateChemInSoilLays(NY,NX,LG,DORGC,TXCO2,DORGE)
!
!     SNOWPACK LAYERING
      call SnowpackLayering(NY,NX)

      call RelayerSoilProfile(NY,NX,DORGC(:,NY,NX),DVLiceMicP(:,NY,NX),UDVLiceMicP,UDLYXF)

      call UpdateOutputVars(I,J,NY,NX,TXCO2)
!
!     CHECK MATERIAL BALANCES
!
!      IF(I.EQ.365.AND.J.EQ.24)THEN
!        WRITE(19,2221)'ORGC',I,J,iYearCurrent,NX,NY,(ORGC(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
!        WRITE(20,2221)'ORGN',I,J,iYearCurrent,NX,NY,(ORGN(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
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
  real(r8) :: VLSoilPoreMicPX,VOLTX
  integer  :: L

  Eco_NetRad_col(NY,NX)=Eco_NetRad_col(NY,NX)+HeatByRadiation(NY,NX)
  Eco_Heat_Latent_col(NY,NX)=Eco_Heat_Latent_col(NY,NX)+HeatEvapAir2Surf(NY,NX)
  Eco_Heat_Sens_col(NY,NX)=Eco_Heat_Sens_col(NY,NX)+HeatSensAir2Surf(NY,NX)
  Eco_Heat_Grnd_col(NY,NX)=Eco_Heat_Grnd_col(NY,NX)-(HeatNet2Surf(NY,NX)-HeatSensVapAir2Surf(NY,NX))
  Canopy_Heat_Latent_col(NY,NX)=Canopy_Heat_Latent_col(NY,NX)+HeatEvapAir2Surf(NY,NX)*BndlResistCanG(NY,NX)
  Canopy_Heat_Sens_col(NY,NX)=Canopy_Heat_Sens_col(NY,NX)+HeatSensAir2Surf(NY,NX)*BndlResistCanG(NY,NX)
  Eco_NEE_col(NY,NX)=Canopy_NEE_col(NY,NX)+SurfGasFlx(idg_CO2,NY,NX)
  ECO_ER_col(NY,NX)=ECO_ER_col(NY,NX)+SurfGasFlx(idg_CO2,NY,NX)
  Eco_NPP_col(NY,NX)=Eco_GPP_col(NY,NX)+Eco_AutoR_col(NY,NX)
  Eco_NBP_col(NY,NX)=Eco_NBP_col(NY,NX)+Canopy_NEE_col(NY,NX) &
    +SurfGasFlx(idg_CO2,NY,NX)+SurfGasFlx(idg_CH4,NY,NX) &
    +TXCO2(NY,NX)-HDOCQ(NY,NX)-HDICQ(NY,NX)-HDOCD(NY,NX)-HDICD(NY,NX)
    
  IF(NU(NY,NX).GT.NUI(NY,NX))THEN  !the surface is lowered
    DO L=NUI(NY,NX),NU(NY,NX)-1
      IF(VLSoilPoreMicP(L,NY,NX).LE.ZEROS2(NY,NX))THEN
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
  THETWZ(0,NY,NX)=AZMAX1((VLWatMicP(0,NY,NX)-VWatLitRHoldCapcity(NY,NX))/AREA(3,0,NY,NX))
  THETIZ(0,NY,NX)=AZMAX1((VLiceMicP(0,NY,NX)-VWatLitRHoldCapcity(NY,NX))/AREA(3,0,NY,NX))
  !THETWZ(0,NY,NX)=AZMAX1(AMIN1(1.0,VLWatMicP(0,NY,NX)/VLitR(NY,NX)))
  !THETIZ(0,NY,NX)=AZMAX1(AMIN1(1.0,VLiceMicP(0,NY,NX)/VLitR(NY,NX)))
  
  D9945: DO L=NUI(NY,NX),NL(NY,NX)
    VLSoilPoreMicPX=AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)*FracSoiAsMicP(L,NY,NX)
    VOLTX=VLSoilPoreMicPX+VLMacP(L,NY,NX)
    THETWZ(L,NY,NX)=safe_adb(VLWatMicP(L,NY,NX)+AMIN1(VLMacP(L,NY,NX),&
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
    TLCO2G=TLCO2G+trcg_solsml(idg_CO2,L,NY,NX)+trcg_solsml(idg_CH4,L,NY,NX)
    DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)+trcg_solsml(idg_CO2,L,NY,NX)+trcg_solsml(idg_CH4,L,NY,NX)
    OXYGSO=OXYGSO+trcg_solsml(idg_O2,L,NY,NX)
    TLN2G=TLN2G+trcg_solsml(idg_N2,L,NY,NX)+trcg_solsml(idg_N2O,L,NY,NX)
    TLNH4=TLNH4+trcn_solsml(ids_NH4,L,NY,NX)+trcg_solsml(idg_NH3,L,NY,NX)
    TLNO3=TLNO3+trcn_solsml(ids_NO3,L,NY,NX)
    TLPO4=TLPO4+trcn_solsml(ids_H1PO4,L,NY,NX)+trcn_solsml(ids_H2PO4,L,NY,NX)

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
  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.ITILL(I,NY,NX).EQ.23)THEN
    ! drainage is on
    DCORPW=DCORP(I,NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
    DTBLI(NY,NX)=DCORPW
    ExtWaterTablet0(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX))
    ExtWaterTable(NY,NX)=ExtWaterTablet0(NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
  ENDIF

  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.ITILL(I,NY,NX).EQ.24)THEN
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

  use ElmIDMod
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

  !     UVLWatMicP(NY,NX)=UVLWatMicP(NY,NX)-VLWatMicP(0,NY,NX)-VLiceMicP(0,NY,NX)*DENSICE
  !
  !     SURFACE BOUNDARY WATER FLUXES
  !
  WI=PrecAtm(NY,NX)+IrrigSurface(NY,NX)   !total incoming water flux=rain/snowfall + irrigation
  CRAIN=CRAIN+WI
  URAIN(NY,NX)=URAIN(NY,NX)+WI
  WO=VapXAir2GSurf(NY,NX)+TEVAPP(NY,NX) !total outgoing water flux
  CEVAP=CEVAP-WO
  UEVAP(NY,NX)=UEVAP(NY,NX)-WO
  VOLWOU=VOLWOU-IrrigSubsurf(NY,NX)
  FWatDischarge(NY,NX)=FWatDischarge(NY,NX)-IrrigSubsurf(NY,NX)
  UVOLO(NY,NX)=UVOLO(NY,NX)-IrrigSubsurf(NY,NX)
  UDRAIN(NY,NX)=UDRAIN(NY,NX)+WaterFlowSoiMicP(3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY HEAT FLUXES
  !
  HEATIN=HEATIN+cpw*TairK(NY,NX)*PrecRainAndSurfirrig(NY,NX)+cps*TairK(NY,NX)*SnoFalPrec(NY,NX)
  HEATIN=HEATIN+HeatNet2Surf(NY,NX)+THFLXC(NY,NX)
  D5150: DO L=1,JS
    HEATIN=HEATIN+XPhaseChangeHeatL(L,NY,NX)
  ENDDO D5150
  HEATOU=HEATOU-cpw*TairK(NY,NX)*IrrigSubsurf(NY,NX)
!
! SURFACE BOUNDARY CO2, CH4 AND DOC FLUXES
! XCODFS: surface - atmosphere CO2 dissolution (+ve) - volatilization (-ve)
! XCOFLG: gaseous CO2 flux, [g d-2 h-1]
! TCO2Z: total root CO2 content
! FLQGQ: precipitation flux into soil surface
! Rain2LitRSurf: precipitation flux into surface litter
! FLQGI: irrifation flux into soil surface
! Irrig2LitRSurf: irrigation flux into surface litter
! XCODFG: soil CO2 dissolution (+ve) - volatilization (-ve)
! XCODFR: soil surface CO2 dissolution (+ve) - volatilization
! UCO2G: total soil CO2 flux, [g d-2]
! HCO2G: hourly soil CO2 flux, [g d-2 h-1]
  CI=GasSfAtmFlx(idg_CO2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_CO2,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_CO2,NY,NX) &
    +(Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX))*CO2_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*CO2_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_CO2,0,NY,NX)+trcg_surf_disevap_flx(idg_CO2,NY,NX)
  CH=GasSfAtmFlx(idg_CH4,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_CH4,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_CH4,NY,NX) &
    +(Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX))*CH4_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*CH4_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_CH4,0,NY,NX)+trcg_surf_disevap_flx(idg_CH4,NY,NX)
  CO=-IrrigSubsurf(NY,NX)*CO2_irrig_conc(NY,NX)
  CX=-IrrigSubsurf(NY,NX)*CH4_irrig_conc(NY,NX)
  SurfGasFlx(idg_CO2,NY,NX)=SurfGasFlx(idg_CO2,NY,NX)+CI
  SurfGasFlx(idg_CH4,NY,NX)=SurfGasFlx(idg_CH4,NY,NX)+CH
  CO2GIN=CO2GIN+CI+CH
  TCOU=TCOU+CO+CX
  !
  !     SURFACE BOUNDARY O2 FLUXES
  !
  OI=GasSfAtmFlx(idg_O2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_O2,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_O2,NY,NX) &
    +(Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX))*O2_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*O2_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_O2,0,NY,NX)+trcg_surf_disevap_flx(idg_O2,NY,NX)
  OXYGIN=OXYGIN+OI
  OO=trcg_RMicbTransf_vr(idg_O2,0,NY,NX)-IrrigSubsurf(NY,NX)*O2_irrig_conc(NY,NX)
  OXYGOU=OXYGOU+OO
  SurfGasFlx(idg_O2,NY,NX)=SurfGasFlx(idg_O2,NY,NX)+OI
  HI=GasSfAtmFlx(idg_H2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_H2,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_H2,NY,NX) &
    +Gas_Disol_Flx_vr(idg_H2,0,NY,NX)+trcg_surf_disevap_flx(idg_H2,NY,NX)
  H2GIN=H2GIN+HI
  HO=trcg_RMicbTransf_vr(idg_H2,0,NY,NX)
  H2GOU=H2GOU+HO
  !
  !     SURFACE BOUNDARY N2, N2O, NH3, NH4, NO3, AND DON FLUXES
  !
  ZSI=((Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX)) &
      *(NH4_rain_conc(NY,NX)+NH3_rain_conc(NY,NX)+NO3_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX)))*14.0
  ZXB=-IrrigSubsurf(NY,NX)*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX))-IrrigSubsurf(NY,NX) &
      *(NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX))*14.0
  TZIN=TZIN+ZSI
  TZOU=TZOU+ZXB
  ZGI=(Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX))*(N2_rain_conc(NY,NX)+N2O_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX)) &
      +GasSfAtmFlx(idg_N2,NY,NX)+GasSfAtmFlx(idg_N2O,NY,NX)+GasSfAtmFlx(idg_NH3,NY,NX) &
      +GasSfAtmFlx(idg_NH3B,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX) &
      +Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX)+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
      +TRootGasLoss_disturb(idg_N2O,NY,NX)+TRootGasLoss_disturb(idg_NH3,NY,NX) &
      +Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,0,NY,NX) &
      +trcg_surf_disevap_flx(idg_N2,NY,NX)+trcg_surf_disevap_flx(idg_N2O,NY,NX) &
      +trcg_surf_disevap_flx(idg_NH3,NY,NX)
  ZN2GIN=ZN2GIN+ZGI
  ZDRAIN(NY,NX)=ZDRAIN(NY,NX)+trcs_3DTransp2MicP(ids_NH4,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP(idg_NH3,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP(ids_NO3,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP(ids_NO2,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP(ids_NH4B,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP(idg_NH3B,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP(ids_NO3B,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP(ids_NO2B,3,NK(NY,NX),NY,NX)

  ZNGGIN=GasSfAtmFlx(idg_N2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX)
  ZN2OIN=GasSfAtmFlx(idg_N2O,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX)+Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)
  ZNH3IN=GasSfAtmFlx(idg_NH3,NY,NX)+GasSfAtmFlx(idg_NH3B,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
    +Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)

  SurfGasFlx(idg_N2,NY,NX)=SurfGasFlx(idg_N2,NY,NX)+ZNGGIN
  SurfGasFlx(idg_N2O,NY,NX)=SurfGasFlx(idg_N2O,NY,NX)+ZN2OIN
  SurfGasFlx(idg_NH3,NY,NX)=SurfGasFlx(idg_NH3,NY,NX)+ZNH3IN
  SurfGasFlx(idg_N2,NY,NX)=SurfGasFlx(idg_N2,NY,NX)+Micb_N2Fixation_vr(0,NY,NX)
  SurfGasFlx(idg_H2,NY,NX)=SurfGasFlx(idg_H2,NY,NX)+HI
  !
  !     SURFACE BOUNDARY PO4 AND DOP FLUXES
  !
  PI=patomw*((Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX)) &
      *(H2PO4_rain_conc(NY,NX)+HPO4_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(H2PO4_irrig_conc(I,NY,NX)+HPO4_irrig_conc(I,NY,NX)))
  PXB=-patomw*IrrigSubsurf(NY,NX)*(H2PO4_irrig_conc(I,NY,NX)+HPO4_irrig_conc(I,NY,NX))
  TPIN=TPIN+PI
  TPOU=TPOU+PXB
  PDRAIN(NY,NX)=PDRAIN(NY,NX)+trcs_3DTransp2MicP(ids_H2PO4,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP(ids_H2PO4B,3,NK(NY,NX),NY,NX)+trcs_3DTransp2MicP(ids_H1PO4,3,NK(NY,NX),NY,NX) &
      +trcs_3DTransp2MicP(ids_H1PO4B,3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY ION FLUXES
  !
  SIN=((Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX)) &
      *(2.0_r8*NH4_rain_conc(NY,NX)+NH3_rain_conc(NY,NX)+NO3_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(2.0*NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX)))
  SGN=(2.0_r8*(Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX))*(N2_rain_conc(NY,NX)+N2O_rain_conc(NY,NX)) &
      +2.0_r8*(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX)) &
      +2.0_r8*(GasSfAtmFlx(idg_N2,NY,NX)+GasSfAtmFlx(idg_N2O,NY,NX))+GasSfAtmFlx(idg_NH3,NY,NX) &
      +GasSfAtmFlx(idg_NH3B,NY,NX)+2.0_r8*(Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX) &
      +Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX))+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
      +2.0_r8*TRootGasLoss_disturb(idg_N2O,NY,NX)+TRootGasLoss_disturb(idg_NH3,NY,NX) &
      +2.0_r8*(Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX))+Gas_Disol_Flx_vr(idg_NH3,0,NY,NX) &
      +2.0_r8*(trcg_surf_disevap_flx(idg_N2,NY,NX)+trcg_surf_disevap_flx(idg_N2O,NY,NX))+trcg_surf_disevap_flx(idg_NH3,NY,NX))/natomw
  SIP=((Rain2SoilSurf(NY,NX)+Rain2LitRSurf(NY,NX))*(3.0_r8*H2PO4_rain_conc(NY,NX)+2.0_r8*HPO4_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(3.0_r8*H2PO4_irrig_conc(I,NY,NX)+2.0_r8*HPO4_irrig_conc(I,NY,NX)))
  SNB=-IrrigSubsurf(NY,NX)*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX))-IrrigSubsurf(NY,NX) &
      *(2.0_r8*NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX))
      SPB=-IrrigSubsurf(NY,NX)*(3.0*H2PO4_irrig_conc(I,NY,NX)+2.0*HPO4_irrig_conc(I,NY,NX))
  SNM0=(2.0_r8*RNutMicbTransf_vr(ids_NH4,0,NY,NX)+RNutMicbTransf_vr(ids_NO3,0,NY,NX)+RNutMicbTransf_vr(ids_NO2,0,NY,NX) &
      -2.0_r8*Micb_N2Fixation_vr(0,NY,NX))/natomw
  SPM0=(2.0_r8*RNutMicbTransf_vr(ids_H1PO4,0,NY,NX)+3.0_r8*RNutMicbTransf_vr(ids_H2PO4,0,NY,NX))/patomw
  !
  !     ACCUMULATE PLANT LITTERFALL FLUXES
  !
  XCSN=XCSN+LitterFallChemElmnt_col(ielmc,NY,NX)
  XZSN=XZSN+LitterFallChemElmnt_col(ielmn,NY,NX)
  XPSN=XPSN+LitterFallChemElmnt_col(ielmp,NY,NX)
  UXCSN(NY,NX)=UXCSN(NY,NX)+LitterFallChemElmnt_col(ielmc,NY,NX)
  UXZSN(NY,NX)=UXZSN(NY,NX)+LitterFallChemElmnt_col(ielmn,NY,NX)
  UXPSN(NY,NX)=UXPSN(NY,NX)+LitterFallChemElmnt_col(ielmp,NY,NX)
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
    SIR=SIR*PrecAtm(NY,NX)
    SBU=-IrrigSubsurf(NY,NX)*SII
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
      DOM(idom,K,0,NY,NX)=DOM(idom,K,0,NY,NX)+DOM_3DMicp_Transp_flx(idom,K,3,0,NY,NX)
    enddo
  ENDDO D9680

  do idg=idg_beg,idg_NH3-1
    trc_solml(idg,0,NY,NX)=trc_solml(idg,0,NY,NX)+trcg_surf_disevap_flx(idg,NY,NX) &
      +trcs_3DTransp2MicP(idg,3,0,NY,NX)+Gas_Disol_Flx_vr(idg,0,NY,NX)-trcg_RMicbTransf_vr(idg,0,NY,NX)
  enddo

  trc_solml(idg_N2,0,NY,NX)=trc_solml(idg_N2,0,NY,NX)-Micb_N2Fixation_vr(0,NY,NX)

  trc_solml(idg_NH3,0,NY,NX)=trc_solml(idg_NH3,0,NY,NX)+trcg_surf_disevap_flx(idg_NH3,NY,NX) &
    +trcs_3DTransp2MicP(idg_NH3,3,0,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)+TR_NH3_soil_vr(0,NY,NX)

  do ids=ids_nut_beg,ids_nuts_end
    trc_solml(ids,0,NY,NX)=trc_solml(ids,0,NY,NX)+trcs_3DTransp2MicP(ids,3,0,NY,NX) &
      +RNutMicbTransf_vr(ids,0,NY,NX)+trcn_RChem_soil_vr(ids,0,NY,NX)
  enddo  

!update aqueous concentrations
  DO NTG=idg_beg,idg_end
    trc_solml(NTG,NU(NY,NX),NY,NX)=trc_solml(NTG,NU(NY,NX),NY,NX)+GasSfAtmFlx(NTG,NY,NX)
  ENDDO

  Eco_HR_col(NY,NX)=Eco_HR_col(NY,NX)+trcg_RMicbTransf_vr(idg_CO2,0,NY,NX)+trcg_RMicbTransf_vr(idg_CH4,0,NY,NX)
  SurfGasFlx(idg_N2,NY,NX)=SurfGasFlx(idg_N2,NY,NX)+trcg_RMicbTransf_vr(idg_N2,0,NY,NX)

  ROXYF(0,NY,NX)=Gas_Disol_Flx_vr(idg_O2,0,NY,NX)
  RCO2F(0,NY,NX)=Gas_Disol_Flx_vr(idg_CO2,0,NY,NX)
  RCH4F(0,NY,NX)=Gas_Disol_Flx_vr(idg_CH4,0,NY,NX)

  !ROXYL:=soil surface O2 dissolution + aqueous O2 flux micropore
  ROXYL(0,NY,NX)=trcg_surf_disevap_flx(idg_O2,NY,NX)+trcs_3DTransp2MicP(idg_O2,3,0,NY,NX) &
    -(Rain2LitRSurf(NY,NX)*O2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*O2_irrig_conc(NY,NX))

  RCH4L(0,NY,NX)=trcg_surf_disevap_flx(idg_CH4,NY,NX)+trcs_3DTransp2MicP(idg_CH4,3,0,NY,NX) &
    -(Rain2LitRSurf(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*CH4_irrig_conc(NY,NX))
    
  ROXYL(NU(NY,NX),NY,NX)=ROXYL(NU(NY,NX),NY,NX)+GasSfAtmFlx(idg_O2,NY,NX)
  RCH4L(NU(NY,NX),NY,NX)=RCH4L(NU(NY,NX),NY,NX)+GasSfAtmFlx(idg_CH4,NY,NX)
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

  integer :: K,NTG
  !     begin_execution
  !
  IF(ABS(TWat2GridBySurfRunoff(NY,NX)).GT.ZEROS(NY,NX))THEN
    !
    !     DOC, DON, DOP
    !
    D8570: DO K=1,micpar%NumOfLitrCmplxs
      DOM(idom_doc,K,0,NY,NX)=DOM(idom_doc,K,0,NY,NX)+TOCQRS(K,NY,NX)
      DOM(idom_don,K,0,NY,NX)=DOM(idom_don,K,0,NY,NX)+TONQRS(K,NY,NX)
      DOM(idom_dop,K,0,NY,NX)=DOM(idom_dop,K,0,NY,NX)+TOPQRS(K,NY,NX)
      DOM(idom_acetate,K,0,NY,NX)=DOM(idom_acetate,K,0,NY,NX)+TOAQRS(K,NY,NX)
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
  integer :: K,NO,M,NGL,NTX,NTP

  REAL(R8) :: DORGP

  ! begin_execution
  !
  ! INTERNAL SURFACE SEDIMENT TRANSPORT
  !
  IF((IERSNG.EQ.1.OR.IERSNG.EQ.3).AND.ABS(TSEDER(NY,NX)).GT.ZEROS(NY,NX))THEN
    TSED(NY,NX)=TSED(NY,NX)+TSEDER(NY,NX)
    !
    !     SOIL MINERAL FRACTIONS
    !
    SAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)+TSANER(NY,NX)
    SILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)+TSILER(NY,NX)
    CLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)+TCLAER(NY,NX)
    !
    !     FERTILIZER POOLS
!
    FertN_soil(ifert_nh4,NU(NY,NX),NY,NX)=FertN_soil(ifert_nh4,NU(NY,NX),NY,NX)+TNH4ER(NY,NX)
    FertN_soil(ifert_nh3,NU(NY,NX),NY,NX)=FertN_soil(ifert_nh3,NU(NY,NX),NY,NX)+TNH3ER(NY,NX)
    FertN_soil(ifert_urea,NU(NY,NX),NY,NX)=FertN_soil(ifert_urea,NU(NY,NX),NY,NX)+TNHUER(NY,NX)
    FertN_soil(ifert_no3,NU(NY,NX),NY,NX)=FertN_soil(ifert_no3,NU(NY,NX),NY,NX)+TNO3ER(NY,NX)
    FertN_band(ifert_nh4_band,NU(NY,NX),NY,NX)=FertN_band(ifert_nh4_band,NU(NY,NX),NY,NX)+TNH4EB(NY,NX)
    FertN_band(ifert_nh3_band,NU(NY,NX),NY,NX)=FertN_band(ifert_nh3_band,NU(NY,NX),NY,NX)+TNH3EB(NY,NX)
    FertN_band(ifert_urea_band,NU(NY,NX),NY,NX)=FertN_band(ifert_urea_band,NU(NY,NX),NY,NX)+TNHUEB(NY,NX)
    FertN_band(ifert_no3_band,NU(NY,NX),NY,NX)=FertN_band(ifert_no3_band,NU(NY,NX),NY,NX)+TNO3EB(NY,NX)
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
      DO  NO=1,NFGs
        DO NGL=JGnio(NO),JGnfo(NO)
          DO  M=1,nlbiomcp
            OMC(M,NGL,K,NU(NY,NX),NY,NX)=OMC(M,NGL,K,NU(NY,NX),NY,NX)+TOMCER(M,NGL,K,NY,NX)
            OMN(M,NGL,K,NU(NY,NX),NY,NX)=OMN(M,NGL,K,NU(NY,NX),NY,NX)+TOMNER(M,NGL,K,NY,NX)
            OMP(M,NGL,K,NU(NY,NX),NY,NX)=OMP(M,NGL,K,NU(NY,NX),NY,NX)+TOMPER(M,NGL,K,NY,NX)
            DORGE(NY,NX)=DORGE(NY,NX)+TOMCER(M,NGL,K,NY,NX)
            DORGP=DORGP+TOMPER(M,NGL,K,NY,NX)
          enddo
        enddo
      enddo
    ENDDO D9280

    DO  NO=1,NFGs
      DO NGL=JGniA(NO),JGnfA(NO)
        DO  M=1,nlbiomcp
          OMCff(M,NGL,NU(NY,NX),NY,NX)=OMCff(M,NGL,NU(NY,NX),NY,NX)+TOMCERff(M,NGL,NY,NX)
          OMNff(M,NGL,NU(NY,NX),NY,NX)=OMNff(M,NGL,NU(NY,NX),NY,NX)+TOMNERff(M,NGL,NY,NX)
          OMPff(M,NGL,NU(NY,NX),NY,NX)=OMPff(M,NGL,NU(NY,NX),NY,NX)+TOMPERff(M,NGL,NY,NX)
          DORGE(NY,NX)=DORGE(NY,NX)+TOMCERff(M,NGL,NY,NX)
          DORGP=DORGP+TOMPER(M,NGL,K,NY,NX)
        enddo
      enddo
    enddo

    D9275: DO K=1,jcplx
      D9270: DO M=1,ndbiomcp
        ORC(M,K,NU(NY,NX),NY,NX)=ORC(M,K,NU(NY,NX),NY,NX)+TORCER(M,K,NY,NX)
        ORN(M,K,NU(NY,NX),NY,NX)=ORN(M,K,NU(NY,NX),NY,NX)+TORNER(M,K,NY,NX)
        ORP(M,K,NU(NY,NX),NY,NX)=ORP(M,K,NU(NY,NX),NY,NX)+TORPER(M,K,NY,NX)
        DORGE(NY,NX)=DORGE(NY,NX)+TORCER(M,K,NY,NX)
        DORGP=DORGP+TORPER(M,K,NY,NX)
      ENDDO D9270
      OHC(K,NU(NY,NX),NY,NX)=OHC(K,NU(NY,NX),NY,NX)+TOHCER(K,NY,NX)
      OHN(K,NU(NY,NX),NY,NX)=OHN(K,NU(NY,NX),NY,NX)+TOHNER(K,NY,NX)
      OHP(K,NU(NY,NX),NY,NX)=OHP(K,NU(NY,NX),NY,NX)+TOHPER(K,NY,NX)
      OHA(K,NU(NY,NX),NY,NX)=OHA(K,NU(NY,NX),NY,NX)+TOHAER(K,NY,NX)
      DORGE(NY,NX)=DORGE(NY,NX)+TOHCER(K,NY,NX)+TOHAER(K,NY,NX)
      DORGP=DORGP+TOHPER(K,NY,NX)
      D9265: DO M=1,jsken
        OSC(M,K,NU(NY,NX),NY,NX)=OSC(M,K,NU(NY,NX),NY,NX)+TOSCER(M,K,NY,NX)
        OSA(M,K,NU(NY,NX),NY,NX)=OSA(M,K,NU(NY,NX),NY,NX)+TOSAER(M,K,NY,NX)
        OSN(M,K,NU(NY,NX),NY,NX)=OSN(M,K,NU(NY,NX),NY,NX)+TOSNER(M,K,NY,NX)
        OSP(M,K,NU(NY,NX),NY,NX)=OSP(M,K,NU(NY,NX),NY,NX)+TOSPER(M,K,NY,NX)
        DORGE(NY,NX)=DORGE(NY,NX)+TOSCER(M,K,NY,NX)
        DORGP=DORGP+TOSPER(M,K,NY,NX)
      ENDDO D9265
    ENDDO D9275
  ENDIF
  end subroutine SoilErosion

!------------------------------------------------------------------------------------------

  subroutine SumLitterLayerChemMass(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: K,N,M,NGL
  real(r8) :: PSS,HS,CS
  real(r8) :: DC,DN,DP,OS
  real(r8) :: POS,POX,POP
  real(r8) :: SSS,ZG,Z4S,WS,Z4X
  real(r8) :: Z4F,ZOS,ZOF
  real(r8) :: tDC,tDN,tDP
  integer :: NumOfLitrCmplxs
!     begin_execution
!     TOTAL C,N,P, SALTS IN SURFACE RESIDUE
!

  NumOfLitrCmplxs   = micpar%NumOfLitrCmplxs

  DC=0.0_r8
  DN=0.0_r8
  DP=0.0_r8
  DO K=1,jcplx
    RC0(K,NY,NX)=0.0_r8
  ENDDO
  RC0ff(NY,NX)=0.0_r8

  OMCL(0,NY,NX)=0.0_r8
  OMNL(0,NY,NX)=0.0_r8

  DO K=1,NumOfLitrCmplxs
    !
    ! TOTAL heterotrophic MICROBIAL C,N,P
    !
    tDC=SUM(OMC(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,K,0,NY,NX))
    tDN=SUM(OMN(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,K,0,NY,NX))
    tDP=SUM(OMP(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,K,0,NY,NX))
    DC=DC+tDC
    DN=DN+tDN
    DP=DP+tDP
    RC0(K,NY,NX)=RC0(K,NY,NX)+tDC
    TOMT(NY,NX)=TOMT(NY,NX)+tDC
    TONT(NY,NX)=TONT(NY,NX)+tDN
    TOPT(NY,NX)=TOPT(NY,NX)+tDP
    OMCL(0,NY,NX)=OMCL(0,NY,NX)+tDC
    OMNL(0,NY,NX)=OMNL(0,NY,NX)+tDN
  ENDDO

  !
  ! TOTAL autotrophic MICROBIAL C,N,P
  !
  tDC=SUM(OMCff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,0,NY,NX))
  tDN=SUM(OMNff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,0,NY,NX))
  tDP=SUM(OMPff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,0,NY,NX))
  DC=DC+tDC
  DN=DN+tDN
  DP=DP+tDP

  RC0ff(NY,NX)=RC0ff(NY,NX)+tDC
  TOMT(NY,NX)=TOMT(NY,NX)+tDC
  TONT(NY,NX)=TONT(NY,NX)+tDN
  TOPT(NY,NX)=TOPT(NY,NX)+tDP
  OMCL(0,NY,NX)=OMCL(0,NY,NX)+tDC
  OMNL(0,NY,NX)=OMNL(0,NY,NX)+tDN

  !
  !     TOTAL MICROBIAL RESIDUE C,N,P
  !
  DC=DC+SUM(ORC(1:ndbiomcp,1:NumOfLitrCmplxs,0,NY,NX))
  DN=DN+SUM(ORN(1:ndbiomcp,1:NumOfLitrCmplxs,0,NY,NX))
  DP=DP+SUM(ORP(1:ndbiomcp,1:NumOfLitrCmplxs,0,NY,NX))
  DO K=1,NumOfLitrCmplxs
    RC0(K,NY,NX)=RC0(K,NY,NX)+SUM(ORC(1:ndbiomcp,K,0,NY,NX))
    RC0(K,NY,NX)=RC0(K,NY,NX)+DOM(idom_doc,K,0,NY,NX)+DOM_Macp(idom_doc,K,0,NY,NX) &
      +OHC(K,0,NY,NX)+DOM(idom_acetate,K,0,NY,NX)+DOM_Macp(idom_acetate,K,0,NY,NX)+OHA(K,0,NY,NX)
    RC0(K,NY,NX)=RC0(K,NY,NX)+SUM(OSC(1:jsken,K,0,NY,NX))
  ENDDO

!
!     TOTAL DOC, DON, DOP
!
  DC=DC+SUM(DOM(idom_doc,1:NumOfLitrCmplxs,0,NY,NX))+SUM(DOM_Macp(idom_doc,1:NumOfLitrCmplxs,0,NY,NX)) &
       +SUM(OHC(1:NumOfLitrCmplxs,0,NY,NX))+SUM(DOM(idom_acetate,1:NumOfLitrCmplxs,0,NY,NX)) &
       +SUM(DOM_Macp(idom_acetate,1:NumOfLitrCmplxs,0,NY,NX))+SUM(OHA(1:NumOfLitrCmplxs,0,NY,NX))

  DN=DN+SUM(DOM(idom_don,1:NumOfLitrCmplxs,0,NY,NX))+SUM(DOM_Macp(idom_don,1:NumOfLitrCmplxs,0,NY,NX)) &
       +SUM(OHN(1:NumOfLitrCmplxs,0,NY,NX))
  DP=DP+SUM(DOM(idom_dop,1:NumOfLitrCmplxs,0,NY,NX))+SUM(DOM_Macp(idom_dop,1:NumOfLitrCmplxs,0,NY,NX)) &
       +SUM(OHP(1:NumOfLitrCmplxs,0,NY,NX))
!
    !     TOTAL PLANT RESIDUE C,N,P
!
  DC=DC+SUM(OSC(1:jsken,1:NumOfLitrCmplxs,0,NY,NX))
  DN=DN+SUM(OSN(1:jsken,1:NumOfLitrCmplxs,0,NY,NX))
  DP=DP+SUM(OSP(1:jsken,1:NumOfLitrCmplxs,0,NY,NX))

  ORGC(0,NY,NX)=DC
  ORGN(0,NY,NX)=DN
  ORGR(0,NY,NX)=DC

  LitRCStoreLandscape=LitRCStoreLandscape+DC
  URSDC(NY,NX)=URSDC(NY,NX)+DC
  LitRNStoreLandscape=LitRNStoreLandscape+DN
  URSDN(NY,NX)=URSDN(NY,NX)+DN
  LitRPStoreLandscape=LitRPStoreLandscape+DP
  URSDP(NY,NX)=URSDP(NY,NX)+DP
  WS=CanH2OHeldVg(NY,NX)+CanWatg(NY,NX)+VLWatMicP(0,NY,NX)+VLiceMicP(0,NY,NX)*DENSICE
  WaterStoreLandscape=WaterStoreLandscape+WS
  UVLWatMicP(NY,NX)=UVLWatMicP(NY,NX)+WS
  HeatStoreLandscape=HeatStoreLandscape+TENGYC(NY,NX)
  CS=trc_solml(idg_CO2,0,NY,NX)+trc_solml(idg_CH4,0,NY,NX)
  TLCO2G=TLCO2G+CS
  DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)+CS
  HS=trc_solml(idg_H2,0,NY,NX)
  TLH2G=TLH2G+HS
  OS=trc_solml(idg_O2,0,NY,NX)
  OXYGSO=OXYGSO+OS
  ZG=trc_solml(idg_N2,0,NY,NX)+trc_solml(idg_N2O,0,NY,NX)
  TLN2G=TLN2G+ZG
  Z4S=trc_solml(ids_NH4,0,NY,NX)+trc_solml(idg_NH3,0,NY,NX)
  Z4X=natomw*trcx_solml(idx_NH4,0,NY,NX)
  Z4F=natomw*(FertN_soil(ifert_nh4,0,NY,NX)+FertN_soil(ifert_urea,0,NY,NX) &
    +FertN_soil(ifert_nh3,0,NY,NX))
  TLNH4=TLNH4+Z4S+Z4X+Z4F
  UNH4(NY,NX)=UNH4(NY,NX)+Z4S+Z4X

  ZOS=trc_solml(ids_NO3,0,NY,NX)+trc_solml(ids_NO2,0,NY,NX)
  ZOF=natomw*FertN_soil(ifert_no3,0,NY,NX)
  TLNO3=TLNO3+ZOS+ZOF
  UNO3(NY,NX)=UNO3(NY,NX)+ZOS
  POS=trc_solml(ids_H1PO4,0,NY,NX)+trc_solml(ids_H2PO4,0,NY,NX)
  POX=patomw*(trcx_solml(idx_HPO4,0,NY,NX)+trcx_solml(idx_H2PO4,0,NY,NX))
  POP=patomw*(trcp_salml(idsp_AlPO4,0,NY,NX)+trcp_salml(idsp_FePO4,0,NY,NX) &
    +trcp_salml(idsp_CaHPO4,0,NY,NX))+2._r8*patomw*trcp_salml(idsp_CaH4P2O8,0,NY,NX) &
    +3._r8*patomw*trcp_salml(idsp_HA,0,NY,NX)
  TLPO4=TLPO4+POS+POX+POP
  UPO4(NY,NX)=UPO4(NY,NX)+POX
  UPP4(NY,NX)=UPP4(NY,NX)+POP

  IF(salt_model)call DiagSurfLitRLayerSalt(NY,NX,TLPO4)

  end subroutine SumLitterLayerChemMass
!------------------------------------------------------------------------------------------

  subroutine DiagSurfLitRLayerSalt(NY,NX,TLPO4)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(inout) :: TLPO4
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
  TLPO4=TLPO4+PSS
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
  REAL(R8),INTENT(OUT) :: DVLiceMicP(JZ,JY,JX)  !change in ice volume
  real(r8) :: TKS00,TKSX
  real(r8) :: ENGY
  real(r8) :: TVHeatCapacity
  real(r8) :: TVHeatCapacitySoilM,TVOLW,TVOLWH,TVOLI,TVOLIH,TENGY
  real(r8) :: VOLWXX,VOLIXX,VHeatCapacityX,WS
  real(r8) :: DVLWatMicP(JZ,JY,JX)   !change in water volume
  integer :: L

  TVHeatCapacity=0.0_r8
  TVHeatCapacitySoilM=0.0_r8
  TVOLW=0.0_r8
  TVOLWH=0.0_r8
  TVOLI=0.0_r8
  TVOLIH=0.0_r8
  TENGY=0.0_r8
  
  DO L=NU(NY,NX),NL(NY,NX)

    TKSX=TKS(L,NY,NX)
    VHeatCapacityX=VHeatCapacity(L,NY,NX)
    VOLWXX=VLWatMicP(L,NY,NX)
    VOLIXX=VLiceMicP(L,NY,NX)
    !micropore
    VLWatMicP(L,NY,NX)=VLWatMicP(L,NY,NX)+TWatFlowCellMicP(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
      +WatIceThawMicP(L,NY,NX)+GridPlantRootH2OUptake_vr(L,NY,NX)+FWatIrrigate2MicP(L,NY,NX)
    VLWatMicPX(L,NY,NX)=VLWatMicPX(L,NY,NX)+TWatFlowCellMicPX(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
      +WatIceThawMicP(L,NY,NX)+GridPlantRootH2OUptake_vr(L,NY,NX)+FWatIrrigate2MicP(L,NY,NX)

    !do a numerical correction
    VLWatMicPX(L,NY,NX)=AMIN1(VLWatMicP(L,NY,NX),VLWatMicPX(L,NY,NX)+0.01_r8*(VLWatMicP(L,NY,NX)-VLWatMicPX(L,NY,NX)))
    VLiceMicP(L,NY,NX)=VLiceMicP(L,NY,NX)-WatIceThawMicP(L,NY,NX)/DENSICE

    !micropore
    VLWatMacP(L,NY,NX)=VLWatMacP(L,NY,NX)+TWaterFlowMacP(L,NY,NX)-FWatExMacP2MicP(L,NY,NX)+WatIceThawMacP(L,NY,NX)
    VLiceMacP(L,NY,NX)=VLiceMacP(L,NY,NX)-WatIceThawMacP(L,NY,NX)/DENSICE

    !volume change
    DVLWatMicP(L,NY,NX)=VLWatMicP1(L,NY,NX)+VLWatMacP1(L,NY,NX)-VLWatMicP(L,NY,NX)-VLWatMacP(L,NY,NX)
    DVLiceMicP(L,NY,NX)=VLiceMicP1(L,NY,NX)+VLiceMacP1(L,NY,NX)-VLiceMicP(L,NY,NX)-VLiceMacP(L,NY,NX)

    !update water/ice-unfilled pores
    IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
      VLsoiAirP(L,NY,NX)=AZMAX1(VLMicP(L,NY,NX)-VLWatMicP(L,NY,NX)-VLiceMicP(L,NY,NX) &
        +VLMacP(L,NY,NX)-VLWatMacP(L,NY,NX)-VLiceMacP(L,NY,NX))
    ELSE
      VLsoiAirP(L,NY,NX)=0.0_r8
!     VLMicP(L,NY,NX)=VLWatMicP(L,NY,NX)+VLiceMicP(L,NY,NX)
!    2+DVLWatMicP(L,NY,NX)+DVLiceMicP(L,NY,NX)
!     VLSoilPoreMicP(L,NY,NX)=VLMicP(L,NY,NX)
!     VGeomLayer(L,NY,NX)=VLMicP(L,NY,NX)
    ENDIF
    ENGY=VHeatCapacityX*TKSX
    VHeatCapacity(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+cpw*(VLWatMicP(L,NY,NX)+VLWatMacP(L,NY,NX)) &
      +cpi*(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX))

    TVHeatCapacity=TVHeatCapacity+VHeatCapacity(L,NY,NX)
    TVHeatCapacitySoilM=TVHeatCapacitySoilM+VHeatCapacitySoilM(L,NY,NX)    
    TVOLW=TVOLW+VLWatMicP(L,NY,NX)
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
    IF(VHeatCapacity(L,NY,NX).GT.ZEROS(NY,NX))THEN
      TKS00=TKS(L,NY,NX)
      TKS(L,NY,NX)=(ENGY+THeatFlow2Soil(L,NY,NX)+THeatSoiThaw(L,NY,NX) &
        +THeatRootUptake(L,NY,NX)+HeatIrrigation(L,NY,NX))/VHeatCapacity(L,NY,NX)

!      if(L==1.and.abs(TKS(L,NY,NX)/TKS00-1._r8)>0.025_r8)then
!        TKS(L,NY,NX)=TKS00
!      endif
    ELSE
      TKS(L,NY,NX)=TKS(NUM(NY,NX),NY,NX)
    ENDIF

    TCS(L,NY,NX)=units%Kelvin2Celcius(TKS(L,NY,NX))
    WS=VLWatMicP(L,NY,NX)+VLWatMacP(L,NY,NX)+(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX))*DENSICE
    WaterStoreLandscape=WaterStoreLandscape+WS
    VOLISO=VOLISO+VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX)
    UVLWatMicP(NY,NX)=UVLWatMicP(NY,NX)+WS
!    2-WiltPoint(L,NY,NX)*VLSoilPoreMicP(L,NY,NX)
    HeatStoreLandscape=HeatStoreLandscape+VHeatCapacity(L,NY,NX)*TKS(L,NY,NX)
  ENDDO
  end subroutine update_physVar_Profile
!------------------------------------------------------------------------------------------
  subroutine UpdateChemInSoilLays(NY,NX,LG,DORGC,TXCO2,DORGE)
  !
  use ElmIDMod
  implicit none
  integer,  intent(in) :: NY,NX,LG
  real(r8), intent(in) :: DORGE(JY,JX)
  real(r8), intent(inout) :: TXCO2(JY,JX)
  real(r8), intent(out) :: DORGC(JZ,JY,JX)

  integer  :: L,K,M,N,LL,NGL,NTX,NTP,NTG,NTS
  real(r8) :: HS,CS
  real(r8) :: CIB,CHB,OIB,COB
  real(r8) :: HGB,HOB,OS,OOB
  real(r8) :: POS,POX,POP
  real(r8) :: SNM,SPM,SSB,SD

  real(r8) :: WX,ZG,Z4S,Z4X,Z4F,ZOS,ZOF
  real(r8) :: ZGB,Z2B,ZHB
  integer  :: idg,idom,ids

  !     begin_execution
  !     UPDATE SOIL LAYER VARIABLES WITH TOTAL FLUXES

  D125: DO L=NU(NY,NX),NL(NY,NX)
    !

    SurfGasFlx(idg_H2,NY,NX)=SurfGasFlx(idg_H2,NY,NX)+Micb_N2Fixation_vr(L,NY,NX)

    !
    !     RESIDUE FROM PLANT LITTERFALL
!
    D8565: DO K=1,micpar%NumOfPlantLitrCmplxs
      DO  M=1,jsken
        OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)+LitrfalChemElemnts_vr(ielmc,M,K,L,NY,NX)
        OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)+LitrfalChemElemnts_vr(ielmc,M,K,L,NY,NX)*micpar%OMCI(1,K)
        OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+LitrfalChemElemnts_vr(ielmn,M,K,L,NY,NX)
        OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+LitrfalChemElemnts_vr(ielmp,M,K,L,NY,NX)
      enddo
    ENDDO D8565
!
    !     DOC, DON, DOP FROM AQUEOUS TRANSPORT
!
    D8560: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM(idom,K,L,NY,NX)=DOM(idom,K,L,NY,NX)+DOM_Transp2Micp_flx(idom,K,L,NY,NX)+DOM_PoreTranspFlx(idom,K,L,NY,NX)
        DOM(idom,K,L,NY,NX)=DOM(idom,K,L,NY,NX)+DOM_Transp2Micp_flx(idom,K,L,NY,NX)+DOM_PoreTranspFlx(idom,K,L,NY,NX)
        DOM(idom,K,L,NY,NX)=DOM(idom,K,L,NY,NX)+DOM_Transp2Micp_flx(idom,K,L,NY,NX)+DOM_PoreTranspFlx(idom,K,L,NY,NX)
        DOM(idom,K,L,NY,NX)=DOM(idom,K,L,NY,NX)+DOM_Transp2Micp_flx(idom,K,L,NY,NX)+DOM_PoreTranspFlx(idom,K,L,NY,NX)

        DOM_Macp(idom,K,L,NY,NX)=DOM_Macp(idom,K,L,NY,NX)+DOM_Transp2Macp_flx(idom,K,L,NY,NX)-DOM_PoreTranspFlx(idom,K,L,NY,NX)
        DOM_Macp(idom,K,L,NY,NX)=DOM_Macp(idom,K,L,NY,NX)+DOM_Transp2Macp_flx(idom,K,L,NY,NX)-DOM_PoreTranspFlx(idom,K,L,NY,NX)
        DOM_Macp(idom,K,L,NY,NX)=DOM_Macp(idom,K,L,NY,NX)+DOM_Transp2Macp_flx(idom,K,L,NY,NX)-DOM_PoreTranspFlx(idom,K,L,NY,NX)
        DOM_Macp(idom,K,L,NY,NX)=DOM_Macp(idom,K,L,NY,NX)+DOM_Transp2Macp_flx(idom,K,L,NY,NX)-DOM_PoreTranspFlx(idom,K,L,NY,NX)
      enddo
    ENDDO D8560
    !
    !     DOC, DON, DOP FROM PLANT EXUDATION
    !
    D195: DO K=1,jcplx
      DOM(idom_doc,K,L,NY,NX)=DOM(idom_doc,K,L,NY,NX)+TDFOME(ielmc,K,L,NY,NX)
      DOM(idom_don,K,L,NY,NX)=DOM(idom_don,K,L,NY,NX)+TDFOME(ielmn,K,L,NY,NX)
      DOM(idom_dop,K,L,NY,NX)=DOM(idom_dop,K,L,NY,NX)+TDFOME(ielmp,K,L,NY,NX)
    ENDDO D195
    !
    !     SOIL SOLUTES FROM AQUEOUS TRANSPORT, MICROBIAL AND ROOT
    !     EXCHANGE, EQUILIBRIUM REACTIONS, GAS EXCHANGE,
    !     MICROPORE-MACROPORE EXCHANGE,
    !
    trc_solml(idg_N2,L,NY,NX)=trc_solml(idg_N2,L,NY,NX) &
      -Micb_N2Fixation_vr(L,NY,NX) 

    trc_solml(idg_CO2,L,NY,NX)=trc_solml(idg_CO2,L,NY,NX) &
      +TR_CO2_aqu_soil_vr(L,NY,NX) 

    do idg=idg_beg,idg_NH3-1
    trc_solml(idg,L,NY,NX)=trc_solml(idg,L,NY,NX) &
      +trcs_Transp2MicP_vr(idg,L,NY,NX)+Gas_Disol_Flx_vr(idg,L,NY,NX) &
      -trcg_RMicbTransf_vr(idg,L,NY,NX)-trcs_plant_uptake_vr(idg,L,NY,NX) &
      +trcs_Irrig_vr(idg,L,NY,NX)+trcs_PoreTranspFlx_vr(idg,L,NY,NX) &
      +trcg_ebu_flx_vr(idg,L,NY,NX)
    enddo

    trc_solml(idg_NH3,L,NY,NX)=trc_solml(idg_NH3,L,NY,NX) &
      +TR_NH3_soil_vr(L,NY,NX) &
      +trcs_Transp2MicP_vr(idg_NH3,L,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,L,NY,NX) &
      -trcs_plant_uptake_vr(idg_NH3,L,NY,NX)&
      +trcs_Irrig_vr(idg_NH3,L,NY,NX) &
      +trcs_PoreTranspFlx_vr(idg_NH3,L,NY,NX)+trcg_ebu_flx_vr(idg_NH3,L,NY,NX)

    trc_solml(idg_NH3B,L,NY,NX)=trc_solml(idg_NH3B,L,NY,NX) &
      +Gas_Disol_Flx_vr(idg_NH3B,L,NY,NX)+trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)+ &
      +trcs_Transp2MicP_vr(idg_NH3B,L,NY,NX) &
      +trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX)-trcs_plant_uptake_vr(idg_NH3B,L,NY,NX) &
      +trcs_Irrig_vr(idg_NH3B,L,NY,NX)+trcs_PoreTranspFlx_vr(idg_NH3B,L,NY,NX)

   DO ids=ids_nut_beg,ids_nuts_end
      trc_solml(ids,L,NY,NX)=trc_solml(ids,L,NY,NX) &
        +trcs_Transp2MicP_vr(ids,L,NY,NX)+RNutMicbTransf_vr(ids,L,NY,NX) &
        +trcn_RChem_soil_vr(ids,L,NY,NX)-trcs_plant_uptake_vr(ids,L,NY,NX) &
        +trcs_Irrig_vr(ids,L,NY,NX)+trcs_PoreTranspFlx_vr(ids,L,NY,NX)
    ENDDO

    do ids=ids_NH4B,ids_nutb_end
      trc_solml(ids,L,NY,NX)=trc_solml(ids,L,NY,NX) &
        +trcs_Transp2MicP_vr(ids,L,NY,NX)+RNutMicbTransf_vr(ids,L,NY,NX) &
        +trcn_RChem_band_soil_vr(ids,L,NY,NX)-trcs_plant_uptake_vr(ids,L,NY,NX) &
        +trcs_Irrig_vr(ids,L,NY,NX)+trcs_PoreTranspFlx_vr(ids,L,NY,NX)
    enddo


    Eco_HR_col(NY,NX)=Eco_HR_col(NY,NX)+trcg_RMicbTransf_vr(idg_CO2,L,NY,NX)+trcg_RMicbTransf_vr(idg_CH4,L,NY,NX)
    SurfGasFlx(idg_N2,NY,NX)=SurfGasFlx(idg_N2,NY,NX)+trcg_RMicbTransf_vr(idg_N2,L,NY,NX)
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
      trc_soHml(NTS,L,NY,NX)=trc_soHml(NTS,L,NY,NX)+trcs_Transp2MacP_vr(NTS,L,NY,NX)-trcs_PoreTranspFlx_vr(NTS,L,NY,NX)
    ENDDO
    !
    !     GASES FROM VOLATILIZATION-DISSOLUTION AND GAS TRANSFER
!
    DO NTG=idg_beg,idg_NH3
      trc_gasml(NTG,L,NY,NX)=trc_gasml(NTG,L,NY,NX)+Gas_AdvDif_Flx_vr(NTG,L,NY,NX)-Gas_Disol_Flx_vr(NTG,L,NY,NX)
    ENDDO

    trc_gasml(idg_NH3,L,NY,NX)=trc_gasml(idg_NH3,L,NY,NX)-Gas_Disol_Flx_vr(idg_NH3B,L,NY,NX)+TRN3G(L,NY,NX)
    ROXYF(L,NY,NX)=Gas_AdvDif_Flx_vr(idg_O2,L,NY,NX)
    RCO2F(L,NY,NX)=Gas_AdvDif_Flx_vr(idg_CO2,L,NY,NX)
    RCH4F(L,NY,NX)=Gas_AdvDif_Flx_vr(idg_CH4,L,NY,NX)
    ROXYL(L,NY,NX)=trcs_Transp2MicP_vr(idg_O2,L,NY,NX)+trcs_Irrig_vr(idg_O2,L,NY,NX) &
      +trcs_PoreTranspFlx_vr(idg_O2,L,NY,NX)+trcg_ebu_flx_vr(idg_O2,L,NY,NX)
    RCH4L(L,NY,NX)=trcs_Transp2MicP_vr(idg_CH4,L,NY,NX)+trcs_Irrig_vr(idg_CH4,L,NY,NX) &
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
        trc_gasml(NTG,LL,NY,NX)=trc_gasml(NTG,LL,NY,NX)-trcg_ebu_flx_vr(NTG,L,NY,NX)
      ENDDO
      trc_gasml(idg_NH3,LL,NY,NX)=trc_gasml(idg_NH3,LL,NY,NX)-trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)

      IF(LG.LT.L)THEN
        TLCO2G=TLCO2G-trcg_ebu_flx_vr(idg_CO2,L,NY,NX)-trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
        DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)-trcg_ebu_flx_vr(idg_CO2,L,NY,NX)-trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
        OXYGSO=OXYGSO-trcg_ebu_flx_vr(idg_O2,L,NY,NX)
        TLN2G=TLN2G-trcg_ebu_flx_vr(idg_N2,L,NY,NX)-trcg_ebu_flx_vr(idg_N2O,L,NY,NX) &
          -trcg_ebu_flx_vr(idg_NH3,L,NY,NX)-trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)
        TLH2G=TLH2G-trcg_ebu_flx_vr(idg_H2,L,NY,NX)
      ENDIF
    ENDIF
    CO2GIN=CO2GIN+CIB+CHB
    COB=TCO2P(L,NY,NX)+trcs_plant_uptake_vr(idg_CO2,L,NY,NX)-TR_CO2_aqu_soil_vr(L,NY,NX)
    TCOU=TCOU+COB
    SurfGasFlx(idg_CO2,NY,NX)=SurfGasFlx(idg_CO2,NY,NX)+CIB
    SurfGasFlx(idg_CH4,NY,NX)=SurfGasFlx(idg_CH4,NY,NX)+CHB
    UCOP(NY,NX)=UCOP(NY,NX)+TCO2P(L,NY,NX)+trcs_plant_uptake_vr(idg_CO2,L,NY,NX)
    HDICD(NY,NX)=HDICD(NY,NX)-12.0*TBCO2(L,NY,NX)
    TXCO2(NY,NX)=TXCO2(NY,NX)+12.0*TBCO2(L,NY,NX)
    OXYGIN=OXYGIN+OIB
    OOB=trcg_RMicbTransf_vr(idg_O2,L,NY,NX)+TUPOXP(L,NY,NX)+trcs_plant_uptake_vr(idg_O2,L,NY,NX)
    OXYGOU=OXYGOU+OOB
    H2GIN=H2GIN+HGB
    HOB=trcg_RMicbTransf_vr(idg_H2,L,NY,NX)+trcs_plant_uptake_vr(idg_H2,L,NY,NX)
    H2GOU=H2GOU+HOB
    ZN2GIN=ZN2GIN+ZGB+Z2B+ZHB
    SurfGasFlx(idg_O2,NY,NX)=SurfGasFlx(idg_O2,NY,NX)+OIB
    SurfGasFlx(idg_N2,NY,NX)=SurfGasFlx(idg_N2,NY,NX)+ZGB
    SurfGasFlx(idg_N2O,NY,NX)=SurfGasFlx(idg_N2O,NY,NX)+Z2B
    SurfGasFlx(idg_NH3,NY,NX)=SurfGasFlx(idg_NH3,NY,NX)+ZHB
    SurfGasFlx(idg_H2,NY,NX)=SurfGasFlx(idg_H2,NY,NX)+HGB
    !
    !     GRID CELL BOUNDARY FLUXES FROM EQUILIBRIUM REACTIONS
!
    SNM=(2.0_r8*(RNutMicbTransf_vr(ids_NH4,L,NY,NX)+RNutMicbTransf_vr(ids_NH4B,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_NH4,L,NY,NX) &
      -trcs_plant_uptake_vr(ids_NH4B,L,NY,NX)-Micb_N2Fixation_vr(L,NY,NX))-trcs_plant_uptake_vr(idg_NH3,L,NY,NX) &
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
    !     UIONOU(NY,NX)=UIONOU(NY,NX)-SSB
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
    CS=trc_gasml(idg_CO2,L,NY,NX)+trc_solml(idg_CO2,L,NY,NX) &
      +trc_soHml(idg_CO2,L,NY,NX)+trcg_TLP(idg_CO2,L,NY,NX) &
      +trc_gasml(idg_CH4,L,NY,NX)+trc_solml(idg_CH4,L,NY,NX) &
      +trc_soHml(idg_CH4,L,NY,NX)+trcg_TLP(idg_CH4,L,NY,NX)
    TLCO2G=TLCO2G+CS
    DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)+CS
    HS=trc_gasml(idg_H2,L,NY,NX)+trc_solml(idg_H2,L,NY,NX) &
      +trc_soHml(idg_H2,L,NY,NX)+trcg_TLP(idg_H2,L,NY,NX)
    TLH2G=TLH2G+HS

    OS=trc_gasml(idg_O2,L,NY,NX)+trc_solml(idg_O2,L,NY,NX) &
      +trc_soHml(idg_O2,L,NY,NX)+trcg_TLP(idg_O2,L,NY,NX)
    OXYGSO=OXYGSO+OS
!  should I add N2 to ZG below? Dec, 18, 2022.
    ZG=trc_gasml(idg_N2,L,NY,NX)+trc_solml(idg_N2,L,NY,NX) &
      +trc_soHml(idg_N2,L,NY,NX)+trcg_TLP(idg_N2O,L,NY,NX) &
      +trc_gasml(idg_N2O,L,NY,NX)+trc_solml(idg_N2O,L,NY,NX) &
      +trc_soHml(idg_N2O,L,NY,NX)+trcg_TLP(idg_NH3,L,NY,NX) &
      +trc_gasml(idg_NH3,L,NY,NX)
    TLN2G=TLN2G+ZG
    Z4S=trc_solml(ids_NH4,L,NY,NX)+trc_soHml(ids_NH4,L,NY,NX) &
      +trc_solml(ids_NH4B,L,NY,NX)+trc_soHml(ids_NH4B,L,NY,NX)&
      +trc_solml(idg_NH3,L,NY,NX)+trc_soHml(idg_NH3,L,NY,NX) &
      +trc_solml(idg_NH3B,L,NY,NX)+trc_soHml(idg_NH3B,L,NY,NX)

    Z4X=natomw*(trcx_solml(idx_NH4,L,NY,NX)+trcx_solml(idx_NH4B,L,NY,NX))
    Z4F=natomw*(FertN_soil(ifert_nh4,L,NY,NX)+FertN_soil(ifert_urea,L,NY,NX) &
      +FertN_soil(ifert_nh3,L,NY,NX)+FertN_band(ifert_nh4_band,L,NY,NX) &
      +FertN_band(ifert_urea_band,L,NY,NX)+FertN_band(ifert_nh3_band,L,NY,NX))
    TLNH4=TLNH4+Z4S+Z4X+Z4F
    UNH4(NY,NX)=UNH4(NY,NX)+Z4S+Z4X

    ZOS=trc_solml(ids_NO3,L,NY,NX)+trc_soHml(ids_NO3,L,NY,NX)+trc_solml(ids_NO3B,L,NY,NX) &
      +trc_soHml(ids_NO3B,L,NY,NX)+trc_solml(ids_NO2,L,NY,NX)+trc_soHml(ids_NO2,L,NY,NX) &
      +trc_solml(ids_NO2B,L,NY,NX)+trc_soHml(ids_NO2B,L,NY,NX)
    ZOF=natomw*(FertN_soil(ifert_no3,L,NY,NX)+FertN_soil(ifert_no3,L,NY,NX))
    TLNO3=TLNO3+ZOS+ZOF
    UNO3(NY,NX)=UNO3(NY,NX)+ZOS

    POS=trc_solml(ids_H2PO4,L,NY,NX)+trc_soHml(ids_H2PO4,L,NY,NX)+trc_solml(ids_H2PO4B,L,NY,NX) &
      +trc_soHml(ids_H2PO4B,L,NY,NX)+trc_solml(ids_H1PO4,L,NY,NX)+trc_soHml(ids_H1PO4,L,NY,NX) &
      +trc_solml(ids_H1PO4B,L,NY,NX)+trc_soHml(ids_H1PO4B,L,NY,NX)

    !exchangeable P
    POX=patomw*(trcx_solml(idx_HPO4,L,NY,NX)+trcx_solml(idx_H2PO4,L,NY,NX) &
      +trcx_solml(idx_HPO4B,L,NY,NX)+trcx_solml(idx_H2PO4B,L,NY,NX))

    !precipitated P
    POP=patomw*(trcp_salml(idsp_AlPO4,L,NY,NX)+trcp_salml(idsp_FePO4,L,NY,NX)&
      +trcp_salml(idsp_CaHPO4,L,NY,NX)+trcp_salml(idsp_AlPO4B,L,NY,NX)&
      +trcp_salml(idsp_FePO4B,L,NY,NX)+trcp_salml(idsp_CaHPO4B,L,NY,NX)) &
      +2._r8*patomw*(trcp_salml(idsp_CaH4P2O8,L,NY,NX)+trcp_salml(idsp_CaH4P2O8B,L,NY,NX)) &
      +3._r8*patomw*(trcp_salml(idsp_HA,L,NY,NX)+trcp_salml(idsp_HAB,L,NY,NX))

    TLPO4=TLPO4+POS+POX+POP
    UPO4(NY,NX)=UPO4(NY,NX)+POX
    UPP4(NY,NX)=UPP4(NY,NX)+POP
    !
    call SumOMStates(L,NY,NX,DORGC(L,NY,NX),DORGE(NY,NX))
!
!     TOTAL SALT IONS
!
    IF(salt_model)call UpdateSaltIonInSoilLayers(L,NY,NX,TLPO4)

  ENDDO D125
  end subroutine UpdateChemInSoilLays

!------------------------------------------------------------------------------------------
  subroutine UpdateSaltIonInSoilLayers(L,NY,NX,TLPO4)
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(inout) :: TLPO4
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

  TLPO4=TLPO4+PSS

  SSS=0._r8
  SSH=0._r8
  DO nsalts=idsalt_beg,idsaltb_end
    SSS=SSS+trcSalt_solml(nsalts,L,NY,NX)*trcSaltIonNumber(nsalts)
    SSH=SSH+trcSalt_soHml(idsalt_Al,L,NY,NX)*trcSaltIonNumber(nsalts)
  ENDDO

 
!     TOTAL FERILIZER,EXCHANGEABLE CATIONS AND ANIONS, PRECIPITATES
!
  SSF=FertN_soil(ifert_nh3,L,NY,NX)+FertN_soil(ifert_urea,L,NY,NX) &
    +FertN_soil(ifert_no3,L,NY,NX) &
    +FertN_band(ifert_nh3_band,L,NY,NX)+FertN_band(ifert_urea_band,L,NY,NX) &
    +FertN_band(ifert_no3_band,L,NY,NX) &
    +2.0_r8*(FertN_soil(ifert_nh4,L,NY,NX)+FertN_band(ifert_nh4_band,L,NY,NX))

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
  IF(VLWatMicP(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    ECHY=0.337_r8*AZMAX1(trcSalt_solml(idsalt_Hp,L,NY,NX)/VLWatMicP(L,NY,NX))
    ECOH=0.192_r8*AZMAX1(trcSalt_solml(idsalt_OH,L,NY,NX)/VLWatMicP(L,NY,NX))
    ECAL=0.056_r8*AZMAX1(trcSalt_solml(idsalt_Al,L,NY,NX)*3.0_r8/VLWatMicP(L,NY,NX))
    ECFE=0.051_r8*AZMAX1(trcSalt_solml(idsalt_Fe,L,NY,NX)*3.0_r8/VLWatMicP(L,NY,NX))
    ECCA=0.060_r8*AZMAX1(trcSalt_solml(idsalt_Ca,L,NY,NX)*2.0_r8/VLWatMicP(L,NY,NX))
    ECMG=0.053_r8*AZMAX1(trcSalt_solml(idsalt_Mg,L,NY,NX)*2.0_r8/VLWatMicP(L,NY,NX))
    ECNA=0.050_r8*AZMAX1(trcSalt_solml(idsalt_Na,L,NY,NX)/VLWatMicP(L,NY,NX))
    ECKA=0.070_r8*AZMAX1(trcSalt_solml(idsalt_K,L,NY,NX)/VLWatMicP(L,NY,NX))
    ECCO=0.072_r8*AZMAX1(trcSalt_solml(idsalt_CO3,L,NY,NX)*2.0_r8/VLWatMicP(L,NY,NX))
    ECHC=0.044_r8*AZMAX1(trcSalt_solml(idsalt_HCO3,L,NY,NX)/VLWatMicP(L,NY,NX))
    ECSO=0.080_r8*AZMAX1(trcSalt_solml(idsalt_SO4,L,NY,NX)*2.0_r8/VLWatMicP(L,NY,NX))
    ECCL=0.076_r8*AZMAX1(trcSalt_solml(idsalt_Cl,L,NY,NX)/VLWatMicP(L,NY,NX))
    ECNO=0.071_r8*AZMAX1(trc_solml(ids_NO3,L,NY,NX)/(VLWatMicP(L,NY,NX)*natomw))
    ECND(L,NY,NX)=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA &
      +ECCO+ECHC+ECSO+ECCL+ECNO

  ELSE
    ECND(L,NY,NX)=0.0_r8
  ENDIF
  end subroutine UpdateSaltIonInSoilLayers

!------------------------------------------------------------------------------------------
  subroutine SumOMStates(L,NY,NX,DORGCC,DORGEC)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in):: DORGEC
  real(r8), intent(out):: DORGCC
  real(r8) :: DC,DN,DP,OC,ON,OP

  integer :: K,N,M,NGL
  integer :: NumOfLitrCmplxs
  real(r8) :: tDC,tDN,tDP
    !     TOTAL SOC,SON,SOP
    !
    !     OMC=microbial biomass, ORC=microbial residue
    !     OQC,OQCH=DOC in micropores,macropores
    !     OQA,OQAH=acetate in micropores,macropores
    !     OHC,OHA=adsorbed SOC,acetate
    !     OSC=SOC(K=0:woody litter, K=1:non-woody litter,
    !     K=2:manure, K=3:POC, K=4:humus)
!
  NumOfLitrCmplxs  = micpar%NumOfLitrCmplxs

  DC=0.0_r8
  DN=0.0_r8
  DP=0.0_r8
  OC=0.0_r8
  ON=0.0_r8
  OP=0.0_r8
  OMCL(L,NY,NX)=0.0_r8
  OMNL(L,NY,NX)=0.0_r8

! add living microbes
  DO K=1,jcplx
    tDC=SUM(OMC(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,K,L,NY,NX))
    tDN=SUM(OMN(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,K,L,NY,NX))
    tDP=SUM(OMP(1:nlbiomcp,1:NumOfMicrobs1HetertrophCmplx,K,L,NY,NX))
    IF(micpar%is_litter(K))THEN  
      !K=0,1,2: woody litr, nonwoody litr, and manure
      DC=DC+tDC
      DN=DN+tDN
      DP=DP+tDP

      TOMT(NY,NX)=TOMT(NY,NX)+tDC
      TONT(NY,NX)=TONT(NY,NX)+tDN
      TOPT(NY,NX)=TOPT(NY,NX)+tDP
      OMCL(L,NY,NX)=OMCL(L,NY,NX)+tDC
      OMNL(L,NY,NX)=OMNL(L,NY,NX)+tDN
    ELSE
      OC=OC+tDC
      ON=ON+tDN
      OP=OP+tDP

      TOMT(NY,NX)=TOMT(NY,NX)+tDC
      TONT(NY,NX)=TONT(NY,NX)+tDN
      TOPT(NY,NX)=TOPT(NY,NX)+tDP
      OMCL(L,NY,NX)=OMCL(L,NY,NX)+tDC
      OMNL(L,NY,NX)=OMNL(L,NY,NX)+tDN
    ENDIF
  ENDDO

! add autotrophs
  tDC=SUM(OMCff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX))
  tDN=SUM(OMNff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX))
  tDP=SUM(OMPff(1:nlbiomcp,1:NumOfMicrobsInAutotrophCmplx,L,NY,NX))
  OC=OC+tDC
  ON=ON+tDN
  OP=OP+tDP
  TOMT(NY,NX)=TOMT(NY,NX)+tDC
  TONT(NY,NX)=TONT(NY,NX)+tDN
  TOPT(NY,NX)=TOPT(NY,NX)+tDP
  OMCL(L,NY,NX)=OMCL(L,NY,NX)+tDC
  OMNL(L,NY,NX)=OMNL(L,NY,NX)+tDN

  DO K=1,jcplx
    IF(micpar%is_litter(K))THEN
! litter + manure
      DC=DC+SUM(ORC(1:ndbiomcp,K,L,NY,NX))
      DN=DN+SUM(ORN(1:ndbiomcp,K,L,NY,NX))
      DP=DP+SUM(ORP(1:ndbiomcp,K,L,NY,NX))
      DC=DC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHC(K,L,NY,NX) &
        +DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+DOM(idom_don,K,L,NY,NX)+DOM_Macp(idom_don,K,L,NY,NX)+OHN(K,L,NY,NX)
      DP=DP+DOM(idom_dop,K,L,NY,NX)+DOM_Macp(idom_dop,K,L,NY,NX)+OHP(K,L,NY,NX)

      DC=DC+SUM(OSC(1:jsken,K,L,NY,NX))
      DN=DN+SUM(OSN(1:jsken,K,L,NY,NX))
      DP=DP+SUM(OSP(1:jsken,K,L,NY,NX))

    ELSE
      OC=OC+SUM(ORC(1:ndbiomcp,K,L,NY,NX))
      ON=ON+SUM(ORN(1:ndbiomcp,K,L,NY,NX))
      OP=OP+SUM(ORP(1:ndbiomcp,K,L,NY,NX))

      OC=OC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHC(K,L,NY,NX) &
        +DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHA(K,L,NY,NX)
      ON=ON+DOM(idom_don,K,L,NY,NX)+DOM_Macp(idom_don,K,L,NY,NX)+OHN(K,L,NY,NX)
      OP=OP+DOM(idom_dop,K,L,NY,NX)+DOM_Macp(idom_dop,K,L,NY,NX)+OHP(K,L,NY,NX)

      OC=OC+SUM(OSC(1:jsken,K,L,NY,NX))
      ON=ON+SUM(OSN(1:jsken,K,L,NY,NX))
      OP=OP+SUM(OSP(1:jsken,K,L,NY,NX))
    ENDIF
  ENDDO
  !litter plus non-litter
  ORGC(L,NY,NX)=DC+OC
  ORGN(L,NY,NX)=DN+ON
  ORGR(L,NY,NX)=DC

  IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
! change in organic C
    DORGCC=ORGCX(L,NY,NX)-ORGC(L,NY,NX)
    IF(L.EQ.NU(NY,NX))THEN
      DORGCC=DORGCC+DORGEC
    ENDIF
  ELSE
    DORGCC=0.0_r8
  ENDIF

  LitRCStoreLandscape=LitRCStoreLandscape+DC
  URSDC(NY,NX)=URSDC(NY,NX)+DC
  LitRNStoreLandscape=LitRNStoreLandscape+DN
  URSDN(NY,NX)=URSDN(NY,NX)+DN
  LitRPStoreLandscape=LitRPStoreLandscape+DP
  URSDP(NY,NX)=URSDP(NY,NX)+DP
  PomHumCStoreLandscape=PomHumCStoreLandscape+OC
  UORGC(NY,NX)=UORGC(NY,NX)+OC
  PomHumNStoreLandscape=PomHumNStoreLandscape+ON
  UORGN(NY,NX)=UORGN(NY,NX)+ON
  PomHumPStoreLandscape=PomHumPStoreLandscape+OP
  UORGP(NY,NX)=UORGP(NY,NX)+OP
  TSEDSO=TSEDSO+(DC+OC)*ppmc
  end subroutine SumOMStates

!------------------------------------------------------------------------------------------

  subroutine AddFlux2SurfaceResidue(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: M,N,K,NGL
  real(r8) :: HRAINR,RAINR

!     begin_execution
!     ADD ABOVE-GROUND LITTERFALL FROM EXTRACT.F TO SURFACE RESIDUE
!
!     OSC,OSN,OSP=SOC,SON,SOP
!     CSNT,ZSNT,PSNT=total C,N,P litterfall
!     ORGC=total SOC
!     RAINR,HRAINR=water,heat in litterfall
!     FLWR,HFLWR=water,heat flux into litter
!     HEATIN=cumulative net surface heat transfer
!
  DO   K=1,micpar%NumOfPlantLitrCmplxs
    DO  M=1,jsken
      OSC(M,K,0,NY,NX)=OSC(M,K,0,NY,NX)+LitrfalChemElemnts_vr(ielmc,M,K,0,NY,NX)
      OSA(M,K,0,NY,NX)=OSA(M,K,0,NY,NX)+LitrfalChemElemnts_vr(ielmc,M,K,0,NY,NX)*micpar%OMCI(1,K)
      OSN(M,K,0,NY,NX)=OSN(M,K,0,NY,NX)+LitrfalChemElemnts_vr(ielmn,M,K,0,NY,NX)
      OSP(M,K,0,NY,NX)=OSP(M,K,0,NY,NX)+LitrfalChemElemnts_vr(ielmp,M,K,0,NY,NX)
      ORGC(0,NY,NX)=ORGC(0,NY,NX)+LitrfalChemElemnts_vr(ielmc,M,K,0,NY,NX)
      RAINR=LitrfalChemElemnts_vr(ielmc,M,K,0,NY,NX)*ThetaCX(K)
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
  !     ROXYX=O2 demand by all microbial,root,myco populations
  !     RNH4X=NH4 demand in non-band by all microbial,root,myco populations
  !     RNO3X=NO3 demand in non-band by all microbial,root,myco populations
  !     RPO4X=H2PO4 demand in non-band by all microbial,root,myco populations
  !     RP14X=HPO4 demand in non-band by all microbial,root,myco populations
  !     RNHBX=NH4 demand in band by all microbial,root,myco populations
  !     RN3BX=NO3 demand in band by all microbial,root,myco populations
  !     RPOBX=H2PO4 demand in band by all microbial,root,myco populations
  !     RP1BX=HPO4 demand in band by all microbial,root,myco populations
  !     ROXYS=O2 demand from DOC,DOA oxidation
  !     RVMX4=demand for NH4 oxidation
  !     RVMX3=demand for NO3 reduction
  !     RVMX2=demand for NO2 oxidation
  !     RVMX1=demand for N2O reduction
  !     RINHO,RINHOR=substrate-unlimited NH4 immobilization
  !     RINOO,RINOOR=substrate-unlimited NO3 immobilization
  !     RIPOO,RIPOOR=substrate-unlimited H2PO4 immobilization
  !     RIPO1,RIPO1R=substrate-unlimited HPO4 immobilization
  !     ROQCX,ROQAX=total DOC,DOA demand from DOC,DOA oxidation
  !     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation
  !     RNO2X=total demand for NO2 reduction
  !     RVMXC=demand for NO2 reduction
!
  DO K=1,jcplx
    IF(micpar%is_litter(K))THEN
      DO  N=1,NFGs
        DO NGL=JGnio(N),JGnfo(N)
          ROXYX(0,NY,NX)=ROXYX(0,NY,NX)+ROXYS(NGL,K,0,NY,NX)
          RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RVMX4(NGL,K,0,NY,NX)
          RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RVMX3(NGL,K,0,NY,NX)
          RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMX2(NGL,K,0,NY,NX)
          RN2OX(0,NY,NX)=RN2OX(0,NY,NX)+RVMX1(NGL,K,0,NY,NX)
          RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RINHO(NGL,K,0,NY,NX)
          RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RINOO(NGL,K,0,NY,NX)
          RPO4X(0,NY,NX)=RPO4X(0,NY,NX)+RIPOO(NGL,K,0,NY,NX)
          RP14X(0,NY,NX)=RP14X(0,NY,NX)+RIPO1(NGL,K,0,NY,NX)
          RNH4X(NU(NY,NX),NY,NX)=RNH4X(NU(NY,NX),NY,NX)+RINHOR(NGL,K,NY,NX)
          RNO3X(NU(NY,NX),NY,NX)=RNO3X(NU(NY,NX),NY,NX)+RINOOR(NGL,K,NY,NX)
          RPO4X(NU(NY,NX),NY,NX)=RPO4X(NU(NY,NX),NY,NX)+RIPOOR(NGL,K,NY,NX)
          RP14X(NU(NY,NX),NY,NX)=RP14X(NU(NY,NX),NY,NX)+RIPO1R(NGL,K,NY,NX)
          ROQCX(K,0,NY,NX)=ROQCX(K,0,NY,NX)+ROQCS(NGL,K,0,NY,NX)
          ROQAX(K,0,NY,NX)=ROQAX(K,0,NY,NX)+ROQAS(NGL,K,0,NY,NX)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DO  N=1,NFGs
    DO NGL=JGniA(N),JGnfA(N)
      ROXYX(0,NY,NX)=ROXYX(0,NY,NX)+ROXYSff(NGL,0,NY,NX)
      RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RVMX4ff(NGL,0,NY,NX)
      RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RVMX3ff(NGL,0,NY,NX)
      RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMX2ff(NGL,0,NY,NX)
      RN2OX(0,NY,NX)=RN2OX(0,NY,NX)+RVMX1ff(NGL,0,NY,NX)
      RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RINHOff(NGL,0,NY,NX)
      RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RINOOff(NGL,0,NY,NX)
      RPO4X(0,NY,NX)=RPO4X(0,NY,NX)+RIPOOff(NGL,0,NY,NX)
      RP14X(0,NY,NX)=RP14X(0,NY,NX)+RIPO1ff(NGL,0,NY,NX)
      RNH4X(NU(NY,NX),NY,NX)=RNH4X(NU(NY,NX),NY,NX)+RINHORff(NGL,NY,NX)
      RNO3X(NU(NY,NX),NY,NX)=RNO3X(NU(NY,NX),NY,NX)+RINOORff(NGL,NY,NX)
      RPO4X(NU(NY,NX),NY,NX)=RPO4X(NU(NY,NX),NY,NX)+RIPOORff(NGL,NY,NX)
      RP14X(NU(NY,NX),NY,NX)=RP14X(NU(NY,NX),NY,NX)+RIPO1Rff(NGL,NY,NX)

    ENDDO
  ENDDO

  RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMXC(0,NY,NX)
  RN2BX(0,NY,NX)=RN2BX(0,NY,NX)+RVMBC(0,NY,NX)
  end subroutine SumSurfMicBGCFluxes

!------------------------------------------------------------------------------------------


  subroutine SumMicBGCFluxes(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX

  integer :: K,N,NGL

  ROXYX(L,NY,NX)=ROXYX(L,NY,NX)+SUM(ROXYS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RNH4X(L,NY,NX)=RNH4X(L,NY,NX)+SUM(RVMX4(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)) &
    +SUM(RINHO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RNO3X(L,NY,NX)=RNO3X(L,NY,NX)+SUM(RVMX3(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)) &
    +SUM(RINOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+SUM(RVMX2(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RN2OX(L,NY,NX)=RN2OX(L,NY,NX)+SUM(RVMX1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RPO4X(L,NY,NX)=RPO4X(L,NY,NX)+SUM(RIPOO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RP14X(L,NY,NX)=RP14X(L,NY,NX)+SUM(RIPO1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RNHBX(L,NY,NX)=RNHBX(L,NY,NX)+SUM(RVMB4(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)) &
    +SUM(RINHB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RN3BX(L,NY,NX)=RN3BX(L,NY,NX)+SUM(RVMB3(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX)) &
    +SUM(RINOB(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+SUM(RVMB2(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RPOBX(L,NY,NX)=RPOBX(L,NY,NX)+SUM(RIPBO(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  RP1BX(L,NY,NX)=RP1BX(L,NY,NX)+SUM(RIPB1(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  DO K=1,jcplx
    ROQCX(K,L,NY,NX)=ROQCX(K,L,NY,NX)+SUM(ROQCS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
    ROQAX(K,L,NY,NX)=ROQAX(K,L,NY,NX)+SUM(ROQAS(1:NumOfMicrobs1HetertrophCmplx,1:jcplx,L,NY,NX))
  ENDDO
  ROXYX(L,NY,NX)=ROXYX(L,NY,NX)+SUM(ROXYSff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RNH4X(L,NY,NX)=RNH4X(L,NY,NX)+SUM(RVMX4ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX)) &
    +SUM(RINHOff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RNO3X(L,NY,NX)=RNO3X(L,NY,NX)+SUM(RVMX3ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX)) &
    +SUM(RINOOff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+SUM(RVMX2ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RN2OX(L,NY,NX)=RN2OX(L,NY,NX)+SUM(RVMX1ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RPO4X(L,NY,NX)=RPO4X(L,NY,NX)+SUM(RIPOOff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RP14X(L,NY,NX)=RP14X(L,NY,NX)+SUM(RIPO1ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RNHBX(L,NY,NX)=RNHBX(L,NY,NX)+SUM(RVMB4ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX)) &
    +SUM(RINHBff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RN3BX(L,NY,NX)=RN3BX(L,NY,NX)+SUM(RVMB3ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX)) &
    +SUM(RINOBff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+SUM(RVMB2ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RPOBX(L,NY,NX)=RPOBX(L,NY,NX)+SUM(RIPBOff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))
  RP1BX(L,NY,NX)=RP1BX(L,NY,NX)+SUM(RIPB1ff(1:NumOfMicrobs1HetertrophCmplx,L,NY,NX))

  RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+RVMXC(L,NY,NX)
  RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+RVMBC(L,NY,NX)

  end subroutine SumMicBGCFluxes
end module RedistMod
