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

  real(r8), pointer :: THETCX(:)



  integer :: curday, curhour
  public :: redist, InitRedist
  contains

  subroutine InitRedist
  implicit none

  allocate(THETCX(micpar%NumOfPlantLitrCmplxs))

  THETCX=(/8.0E-06_r8,8.0E-06_r8/)

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
      call UpdateSnowChemFlux(NY,NX)
!
      call CalcLitterLayerChemicalMass(NY,NX)
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
!        WRITE(19,2221)'ORGC',I,J,IYRC,NX,NY,(ORGC(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
!        WRITE(20,2221)'ORGN',I,J,IYRC,NX,NY,(ORGN(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
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
  real(r8), intent(in) :: TXCO2(JY,JX)
  real(r8) :: VLSoilPoreMicPX,VOLTX
  integer  :: L
  TRN(NY,NX)=TRN(NY,NX)+HeatByRadiation(NY,NX)
  TLE(NY,NX)=TLE(NY,NX)+HeatEvapAir2Surf(NY,NX)
  TSH(NY,NX)=TSH(NY,NX)+HeatSensAir2Surf(NY,NX)
  TGH(NY,NX)=TGH(NY,NX)-(HeatNet2Surf(NY,NX)-HeatSensVapAir2Surf(NY,NX))
  TLEC(NY,NX)=TLEC(NY,NX)+HeatEvapAir2Surf(NY,NX)*BndlResistCanG(NY,NX)
  TSHC(NY,NX)=TSHC(NY,NX)+HeatSensAir2Surf(NY,NX)*BndlResistCanG(NY,NX)
  TCNET(NY,NX)=TCCAN(NY,NX)+HCO2G(NY,NX)
  RECO(NY,NX)=RECO(NY,NX)+HCO2G(NY,NX)
  TCAN(NY,NX)=TCAN(NY,NX)+TCCAN(NY,NX)
  TNPP(NY,NX)=TGPP(NY,NX)+TRAU(NY,NX)
  TNBP(NY,NX)=TCAN(NY,NX)+UCO2G(NY,NX)+UCH4G(NY,NX) &
    -UDOCQ(NY,NX)-UDICQ(NY,NX)-UDOCD(NY,NX)-UDICD(NY,NX)+TXCO2(NY,NX)
    
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

  subroutine UpdateSnowChemFlux(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: SSW,ENGYW,WS
  integer :: L
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
    UCO2S(NY,NX)=UCO2S(NY,NX)+trcg_solsml(idg_CO2,L,NY,NX)+trcg_solsml(idg_CH4,L,NY,NX)
    OXYGSO=OXYGSO+trcg_solsml(idg_O2,L,NY,NX)
    TLN2G=TLN2G+trcg_solsml(idg_N2,L,NY,NX)+trcg_solsml(idg_N2O,L,NY,NX)
    TLNH4=TLNH4+trcn_solsml(ids_NH4,L,NY,NX)+trcg_solsml(idg_NH3,L,NY,NX)
    TLNO3=TLNO3+trcn_solsml(ids_NO3,L,NY,NX)
    TLPO4=TLPO4+trcn_solsml(ids_H1PO4,L,NY,NX)+trcn_solsml(ids_H2PO4,L,NY,NX)

    IF(salt_model)THEN
      SSW=trcs_solsml(idsalt_Al,L,NY,NX)+trcs_solsml(idsalt_Fe,L,NY,NX)+trcs_solsml(idsalt_Hp,L,NY,NX)+trcs_solsml(idsalt_Ca,L,NY,NX) &
        +trcs_solsml(idsalt_Mg,L,NY,NX)+trcs_solsml(idsalt_Na,L,NY,NX)+trcs_solsml(idsalt_K,L,NY,NX)+trcs_solsml(idsalt_OH,L,NY,NX) &
        +trcs_solsml(idsalt_SO4,L,NY,NX)+trcs_solsml(idsalt_Cl,L,NY,NX)+trcs_solsml(idsalt_CO3,L,NY,NX)+trcs_solsml(idsalt_H0PO4,L,NY,NX) &
        +2.0_r8*(trcs_solsml(idsalt_HCO3,L,NY,NX)+trcs_solsml(idsalt_AlOH,L,NY,NX) &
        +trcs_solsml(idsalt_AlSO4,L,NY,NX)+trcs_solsml(idsalt_FeOH,L,NY,NX)+trcs_solsml(idsalt_FeSO4,L,NY,NX)+trcs_solsml(idsalt_CaOH2,L,NY,NX) &
        +trcs_solsml(idsalt_CaCO3,L,NY,NX)+trcs_solsml(idsalt_CaSO4,L,NY,NX)+trcs_solsml(idsalt_MgOH2,L,NY,NX)+trcs_solsml(idsalt_MgCO3,L,NY,NX) &
        +trcs_solsml(idsalt_MgSO4,L,NY,NX)+trcs_solsml(idsalt_NaCO3,L,NY,NX)+trcs_solsml(idsalt_NaSO4,L,NY,NX)+trcs_solsml(idsalt_KSO4,L,NY,NX) &
        +trcs_solsml(idsalt_CaPO4,L,NY,NX)) &
        +3.0*(trcs_solsml(idsalt_AlOH2,L,NY,NX)+trcs_solsml(idsalt_FeOH2,L,NY,NX)+trcs_solsml(idsalt_CaHCO3,L,NY,NX) &
        +trcs_solsml(idsalt_MgHCO3,L,NY,NX)+trcs_solsml(idsalt_FeHPO4,L,NY,NX)+trcs_solsml(idsalt_CaHPO4,L,NY,NX) &
        +trcs_solsml(idsalt_MgHPO4,L,NY,NX)) &
        +4.0*(trcs_solsml(idsalt_AlOH3,L,NY,NX)+trcs_solsml(idsalt_FeOH3,L,NY,NX)+trcs_solsml(idsalt_H3PO4,L,NY,NX) &
        +trcs_solsml(idsalt_FeH2PO4,L,NY,NX)+trcs_solsml(idsalt_CaH2PO4,L,NY,NX)) &
        +5.0*(trcs_solsml(idsalt_AlOH4,L,NY,NX)+trcs_solsml(idsalt_FeOH4,L,NY,NX))
      TION=TION+SSW

    ENDIF
  ENDDO
  end subroutine UpdateSnowChemFlux
!------------------------------------------------------------------------------------------

  subroutine ModifyExWTBLByDisturbance(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: DCORPW
!     begin_execution
!     ZNOON=hour of solar noon from weather file
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
  IF(J.EQ.INT(ZNOON(NY,NX)).AND.ITILL(I,NY,NX).EQ.23)THEN
    ! drainage is on
    DCORPW=DCORP(I,NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
    DTBLI(NY,NX)=DCORPW
    ExtWaterTablet0(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX))
    ExtWaterTable(NY,NX)=ExtWaterTablet0(NY,NX)+CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
  ENDIF

  IF(J.EQ.INT(ZNOON(NY,NX)).AND.ITILL(I,NY,NX).EQ.24)THEN
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

  integer :: L,K,LS,NTG,NTP,NTX
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
  ! begin_execution

  !
  call UpdateLitRPhys(NY,NX,TWat2GridBySurfRunoff(NY,NX),THeat2GridBySurfRunoff(NY,NX),HeatStoreLandscape,HEATIN)

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
! FLQRQ: precipitation flux into surface litter
! FLQGI: irrifation flux into soil surface
! FLQRI: irrigation flux into surface litter
! XCODFG: soil CO2 dissolution (+ve) - volatilization (-ve)
! XCODFR: soil surface CO2 dissolution (+ve) - volatilization
! UCO2G: total soil CO2 flux, [g d-2]
! HCO2G: hourly soil CO2 flux, [g d-2 h-1]
  CI=GasSfAtmFlx(idg_CO2,NY,NX)+R3GasADTFlx(idg_CO2,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_CO2,NY,NX) &
    +(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX) &
    +(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX) &
    +GasDisFlx(idg_CO2,0,NY,NX)+trcg_XDFR(idg_CO2,NY,NX)
  CH=GasSfAtmFlx(idg_CH4,NY,NX)+R3GasADTFlx(idg_CH4,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_CH4,NY,NX) &
    +(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCHR(NY,NX) &
    +(FLQGI(NY,NX)+FLQRI(NY,NX))*CCHQ(NY,NX) &
    +GasDisFlx(idg_CH4,0,NY,NX)+trcg_XDFR(idg_CH4,NY,NX)
  CO=-IrrigSubsurf(NY,NX)*CCOQ(NY,NX)
  CX=-IrrigSubsurf(NY,NX)*CCHQ(NY,NX)
  UCO2G(NY,NX)=UCO2G(NY,NX)+CI
  HCO2G(NY,NX)=HCO2G(NY,NX)+CI
  UCH4G(NY,NX)=UCH4G(NY,NX)+CH
  HCH4G(NY,NX)=HCH4G(NY,NX)+CH
  CO2GIN=CO2GIN+CI+CH
  TCOU=TCOU+CO+CX
  !
  !     SURFACE BOUNDARY O2 FLUXES
  !
  OI=GasSfAtmFlx(idg_O2,NY,NX)+R3GasADTFlx(idg_O2,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_O2,NY,NX) &
    +(FLQGQ(NY,NX)+FLQRQ(NY,NX))*COXR(NY,NX) &
    +(FLQGI(NY,NX)+FLQRI(NY,NX))*COXQ(NY,NX) &
    +GasDisFlx(idg_O2,0,NY,NX)+trcg_XDFR(idg_O2,NY,NX)
  OXYGIN=OXYGIN+OI
  OO=RUPOXO(0,NY,NX)-IrrigSubsurf(NY,NX)*COXQ(NY,NX)
  OXYGOU=OXYGOU+OO
  UOXYG(NY,NX)=UOXYG(NY,NX)+OI
  HOXYG(NY,NX)=HOXYG(NY,NX)+OI
  HI=GasSfAtmFlx(idg_H2,NY,NX)+R3GasADTFlx(idg_H2,3,NU(NY,NX),NY,NX) &
    +TRootGasLoss_disturb(idg_H2,NY,NX) &
    +GasDisFlx(idg_H2,0,NY,NX)+trcg_XDFR(idg_H2,NY,NX)
  H2GIN=H2GIN+HI
  HO=RH2GO(0,NY,NX)
  H2GOU=H2GOU+HO
  !
  !     SURFACE BOUNDARY N2, N2O, NH3, NH4, NO3, AND DON FLUXES
  !
  ZSI=((FLQGQ(NY,NX)+FLQRQ(NY,NX)) &
      *(CN4R(NY,NX)+CN3R(NY,NX)+CNOR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX)) &
      *(CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX)))*14.0
  ZXB=-IrrigSubsurf(NY,NX)*(CNNQ(NY,NX)+CN2Q(NY,NX))-IrrigSubsurf(NY,NX) &
      *(CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX))*14.0
  TZIN=TZIN+ZSI
  TZOU=TZOU+ZXB
  ZGI=(FLQGQ(NY,NX)+FLQRQ(NY,NX))*(CNNR(NY,NX)+CN2R(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX))*(CNNQ(NY,NX)+CN2Q(NY,NX)) &
      +GasSfAtmFlx(idg_N2,NY,NX)+GasSfAtmFlx(idg_N2O,NY,NX)+GasSfAtmFlx(idg_NH3,NY,NX) &
      +GasSfAtmFlx(idg_NH3B,NY,NX)+R3GasADTFlx(idg_N2,3,NU(NY,NX),NY,NX) &
      +R3GasADTFlx(idg_N2O,3,NU(NY,NX),NY,NX)+R3GasADTFlx(idg_NH3,3,NU(NY,NX),NY,NX) &
      +TRootGasLoss_disturb(idg_N2O,NY,NX)+TRootGasLoss_disturb(idg_NH3,NY,NX) &
      +GasDisFlx(idg_N2O,0,NY,NX)+GasDisFlx(idg_N2,0,NY,NX)+GasDisFlx(idg_NH3,0,NY,NX) &
      +trcg_XDFR(idg_N2,NY,NX)+trcg_XDFR(idg_N2O,NY,NX)+trcg_XDFR(idg_NH3,NY,NX)
  ZN2GIN=ZN2GIN+ZGI
  ZDRAIN(NY,NX)=ZDRAIN(NY,NX)+trcs_XFLS(ids_NH4,3,NK(NY,NX),NY,NX) &
      +trcs_XFLS(idg_NH3,3,NK(NY,NX),NY,NX)+trcs_XFLS(ids_NO3,3,NK(NY,NX),NY,NX) &
      +trcs_XFLS(ids_NO2,3,NK(NY,NX),NY,NX)+trcs_XFLS(ids_NH4B,3,NK(NY,NX),NY,NX) &
      +trcs_XFLS(idg_NH3B,3,NK(NY,NX),NY,NX)+trcs_XFLS(ids_NO3B,3,NK(NY,NX),NY,NX) &
      +trcs_XFLS(ids_NO2B,3,NK(NY,NX),NY,NX)

  ZNGGIN=GasSfAtmFlx(idg_N2,NY,NX)+R3GasADTFlx(idg_N2,3,NU(NY,NX),NY,NX)+GasDisFlx(idg_N2,0,NY,NX)
  ZN2OIN=GasSfAtmFlx(idg_N2O,NY,NX)+R3GasADTFlx(idg_N2O,3,NU(NY,NX),NY,NX)+GasDisFlx(idg_N2O,0,NY,NX)
  ZNH3IN=GasSfAtmFlx(idg_NH3,NY,NX)+GasSfAtmFlx(idg_NH3B,NY,NX)+R3GasADTFlx(idg_NH3,3,NU(NY,NX),NY,NX) &
    +GasDisFlx(idg_NH3,0,NY,NX)
  ! UN2GG(NY,NX)=UN2GG(NY,NX)+ZNGGIN
  !     HN2GG(NY,NX)=HN2GG(NY,NX)+ZNGGIN
  UN2OG(NY,NX)=UN2OG(NY,NX)+ZN2OIN
  HN2OG(NY,NX)=HN2OG(NY,NX)+ZN2OIN
  UNH3G(NY,NX)=UNH3G(NY,NX)+ZNH3IN
  HNH3G(NY,NX)=HNH3G(NY,NX)+ZNH3IN
  UN2GS(NY,NX)=UN2GS(NY,NX)+XN2GS(0,NY,NX)
  UH2GG(NY,NX)=UH2GG(NY,NX)+HI
  !
  !     SURFACE BOUNDARY PO4 AND DOP FLUXES
  !
  PI=patomw*((FLQGQ(NY,NX)+FLQRQ(NY,NX)) &
      *(CPOR(NY,NX)+CH1PR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX)) &
      *(CPOQ(I,NY,NX)+CH1PQ(I,NY,NX)))
  PXB=-patomw*IrrigSubsurf(NY,NX)*(CPOQ(I,NY,NX)+CH1PQ(I,NY,NX))
  TPIN=TPIN+PI
  TPOU=TPOU+PXB
  PDRAIN(NY,NX)=PDRAIN(NY,NX)+trcs_XFLS(ids_H2PO4,3,NK(NY,NX),NY,NX) &
      +trcs_XFLS(ids_H2PO4B,3,NK(NY,NX),NY,NX)+trcs_XFLS(ids_H1PO4,3,NK(NY,NX),NY,NX) &
      +trcs_XFLS(ids_H1PO4B,3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY ION FLUXES
  !
  SIN=((FLQGQ(NY,NX)+FLQRQ(NY,NX)) &
      *(2.0_r8*CN4R(NY,NX)+CN3R(NY,NX)+CNOR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX)) &
      *(2.0*CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX)))
  SGN=(2.0_r8*(FLQGQ(NY,NX)+FLQRQ(NY,NX))*(CNNR(NY,NX)+CN2R(NY,NX)) &
      +2.0_r8*(FLQGI(NY,NX)+FLQRI(NY,NX))*(CNNQ(NY,NX)+CN2Q(NY,NX)) &
      +2.0_r8*(GasSfAtmFlx(idg_N2,NY,NX)+GasSfAtmFlx(idg_N2O,NY,NX))+GasSfAtmFlx(idg_NH3,NY,NX) &
      +GasSfAtmFlx(idg_NH3B,NY,NX)+2.0_r8*(R3GasADTFlx(idg_N2,3,NU(NY,NX),NY,NX) &
      +R3GasADTFlx(idg_N2O,3,NU(NY,NX),NY,NX))+R3GasADTFlx(idg_NH3,3,NU(NY,NX),NY,NX) &
      +2.0_r8*TRootGasLoss_disturb(idg_N2O,NY,NX)+TRootGasLoss_disturb(idg_NH3,NY,NX) &
      +2.0_r8*(GasDisFlx(idg_N2O,0,NY,NX)+GasDisFlx(idg_N2,0,NY,NX))+GasDisFlx(idg_NH3,0,NY,NX) &
      +2.0_r8*(trcg_XDFR(idg_N2,NY,NX)+trcg_XDFR(idg_N2O,NY,NX))+trcg_XDFR(idg_NH3,NY,NX))/natomw
  SIP=((FLQGQ(NY,NX)+FLQRQ(NY,NX))*(3.0_r8*CPOR(NY,NX)+2.0_r8*CH1PR(NY,NX)) &
      +(FLQGI(NY,NX)+FLQRI(NY,NX))*(3.0_r8*CPOQ(I,NY,NX)+2.0_r8*CH1PQ(I,NY,NX)))
  SNB=-IrrigSubsurf(NY,NX)*(CNNQ(NY,NX)+CN2Q(NY,NX))-IrrigSubsurf(NY,NX) &
      *(2.0_r8*CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX))
      SPB=-IrrigSubsurf(NY,NX)*(3.0*CPOQ(I,NY,NX)+2.0*CH1PQ(I,NY,NX))
  SNM0=(2.0_r8*XNH4S(0,NY,NX)+XNO3S(0,NY,NX)+XNO2S(0,NY,NX) &
      -2.0_r8*XN2GS(0,NY,NX))/natomw
  SPM0=(2.0_r8*XH1PS(0,NY,NX)+3.0_r8*XH2PS(0,NY,NX))/patomw
  !
  !     ACCUMULATE PLANT LITTERFALL FLUXES
  !
  XCSN=XCSN+ZESNC(ielmc,NY,NX)
  XZSN=XZSN+ZESNC(ielmn,NY,NX)
  XPSN=XPSN+ZESNC(ielmp,NY,NX)
  UXCSN(NY,NX)=UXCSN(NY,NX)+ZESNC(ielmc,NY,NX)
  UXZSN(NY,NX)=UXZSN(NY,NX)+ZESNC(ielmn,NY,NX)
  UXPSN(NY,NX)=UXPSN(NY,NX)+ZESNC(ielmp,NY,NX)
  !
  !     SURFACE BOUNDARY SALT FLUXES FROM RAINFALL AND SURFACE IRRIGATION
  !
  IF(salt_model)THEN
    SIR=PrecAtm(NY,NX)*(CALR(NY,NX)+CFER(NY,NX)+CHYR(NY,NX)+CCAR(NY,NX) &
      +CMGR(NY,NX)+CNAR(NY,NX)+CKAR(NY,NX)+COHR(NY,NX)+CSOR(NY,NX) &
      +CCLR(NY,NX)+CC3R(NY,NX)+CH0PR(NY,NX) &
      +2.0_r8*(CHCR(NY,NX)+CAL1R(NY,NX)+CALSR(NY,NX)+CFE1R(NY,NX) &
      +CFESR(NY,NX)+CCAOR(NY,NX)+CCACR(NY,NX)+CCASR(NY,NX)+CMGOR(NY,NX) &
      +CMGCR(NY,NX)+CMGSR(NY,NX)+CNACR(NY,NX)+CNASR(NY,NX) &
      +CKASR(NY,NX)+CC0PR(NY,NX)) &
      +3.0_r8*(CAL2R(NY,NX)+CFE2R(NY,NX)+CCAHR(NY,NX)+CMGHR(NY,NX) &
      +CF1PR(NY,NX)+CC1PR(NY,NX)+CM1PR(NY,NX)) &
      +4.0_r8*(CAL3R(NY,NX)+CFE3R(NY,NX)+CH3PR(NY,NX)+CF2PR(NY,NX) &
      +CC2PR(NY,NX)) &
      +5.0_r8*(CAL4R(NY,NX)+CFE4R(NY,NX)))
    SII=IrrigSurface(NY,NX)*(CALQ(I,NY,NX)+CFEQ(I,NY,NX)+CHYQ(I,NY,NX) &
      +CCAQ(I,NY,NX)+CMGQ(I,NY,NX)+CNAQ(I,NY,NX)+CKAQ(I,NY,NX) &
      +COHQ(I,NY,NX)+CSOQ(I,NY,NX)+CCLQ(I,NY,NX)+CC3Q(I,NY,NX) &
      +CH0PQ(I,NY,NX) &
      +2.0_r8*(CHCQ(I,NY,NX)+CAL1Q(I,NY,NX)+CALSQ(I,NY,NX) &
      +CFE1Q(I,NY,NX)+CFESQ(I,NY,NX)+CCAOQ(I,NY,NX)+CCACQ(I,NY,NX) &
      +CCASQ(I,NY,NX)+CMGOQ(I,NY,NX)+CMGCQ(I,NY,NX)+CMGSQ(I,NY,NX) &
      +CNACQ(I,NY,NX)+CNASQ(I,NY,NX)+CKASQ(I,NY,NX)+CC0PQ(I,NY,NX)) &
      +3.0_r8*(CAL2Q(I,NY,NX)+CFE2Q(I,NY,NX)+CCAHQ(I,NY,NX) &
      +CMGHQ(I,NY,NX)+CF1PQ(I,NY,NX)+CC1PQ(I,NY,NX)+CM1PQ(I,NY,NX)) &
      +4.0_r8*(CAL3Q(I,NY,NX)+CFE3Q(I,NY,NX) &
      +CH3PQ(I,NY,NX)+CF2PQ(I,NY,NX)+CC2PQ(I,NY,NX)) &
      +5.0_r8*(CAL4Q(I,NY,NX)+CFE4Q(I,NY,NX)))
    TIONIN=TIONIN+SIR+SII
    !
    !     SUBSURFACE BOUNDARY SALT FLUXES FROM SUBSURFACE IRRIGATION
    !
    SBU=-IrrigSubsurf(NY,NX)*(CALQ(I,NY,NX)+CFEQ(I,NY,NX)+CHYQ(I,NY,NX) &
      +CCAQ(I,NY,NX)+CMGQ(I,NY,NX)+CNAQ(I,NY,NX)+CKAQ(I,NY,NX) &
      +COHQ(I,NY,NX)+CSOQ(I,NY,NX)+CCLQ(I,NY,NX)+CC3Q(I,NY,NX) &
      +CH0PQ(I,NY,NX) &
      +2.0_r8*(CHCQ(I,NY,NX)+CAL1Q(I,NY,NX)+CALSQ(I,NY,NX) &
      +CFE1Q(I,NY,NX)+CFESQ(I,NY,NX)+CCAOQ(I,NY,NX)+CCACQ(I,NY,NX) &
      +CCASQ(I,NY,NX)+CMGOQ(I,NY,NX)+CMGCQ(I,NY,NX)+CMGSQ(I,NY,NX) &
      +CNACQ(I,NY,NX)+CNASQ(I,NY,NX)+CKASQ(I,NY,NX)+CC0PQ(I,NY,NX)) &
      +3.0_r8*(CAL2Q(I,NY,NX)+CFE2Q(I,NY,NX)+CCAHQ(I,NY,NX)+CMGHQ(I,NY,NX) &
      +CF1PQ(I,NY,NX)+CC1PQ(I,NY,NX)+CM1PQ(I,NY,NX)) &
      +4.0_r8*(CAL3Q(I,NY,NX)+CFE3Q(I,NY,NX) &
      +CH3PQ(I,NY,NX)+CF2PQ(I,NY,NX)+CC2PQ(I,NY,NX)) &
      +5.0_r8*(CAL4Q(I,NY,NX)+CFE4Q(I,NY,NX)))
    TIONOU=TIONOU+SBU
  ENDIF
  !
  ! GAS EXCHANGE FROM SURFACE VOLATILIZATION-DISSOLUTION
  !
  D9680: DO K=1,micpar%NumOfLitrCmplxs
    OQC(K,0,NY,NX)=OQC(K,0,NY,NX)+XOCFLS(K,3,0,NY,NX)
    OQN(K,0,NY,NX)=OQN(K,0,NY,NX)+XONFLS(K,3,0,NY,NX)
    OQP(K,0,NY,NX)=OQP(K,0,NY,NX)+XOPFLS(K,3,0,NY,NX)
    OQA(K,0,NY,NX)=OQA(K,0,NY,NX)+XOAFLS(K,3,0,NY,NX)
  ENDDO D9680

  trc_solml(idg_CO2,0,NY,NX)=trc_solml(idg_CO2,0,NY,NX)+trcg_XDFR(idg_CO2,NY,NX) &
    +trcs_XFLS(idg_CO2,3,0,NY,NX)+GasDisFlx(idg_CO2,0,NY,NX)-RCO2O(0,NY,NX)
  trc_solml(idg_CH4,0,NY,NX)=trc_solml(idg_CH4,0,NY,NX)+trcg_XDFR(idg_CH4,NY,NX) &
    +trcs_XFLS(idg_CH4,3,0,NY,NX)+GasDisFlx(idg_CH4,0,NY,NX)-RCH4O(0,NY,NX)
  trc_solml(idg_O2,0,NY,NX)=trc_solml(idg_O2,0,NY,NX)+trcg_XDFR(idg_O2,NY,NX) &
    +trcs_XFLS(idg_O2,3,0,NY,NX)+GasDisFlx(idg_O2,0,NY,NX)-RUPOXO(0,NY,NX)
  trc_solml(idg_N2,0,NY,NX)=trc_solml(idg_N2,0,NY,NX)+trcg_XDFR(idg_N2,NY,NX) &
    +trcs_XFLS(idg_N2,3,0,NY,NX)+GasDisFlx(idg_N2,0,NY,NX)-RN2G(0,NY,NX)-XN2GS(0,NY,NX)
  trc_solml(idg_N2O,0,NY,NX)=trc_solml(idg_N2O,0,NY,NX)+trcg_XDFR(idg_N2O,NY,NX) &
    +trcs_XFLS(idg_N2O,3,0,NY,NX)+GasDisFlx(idg_N2O,0,NY,NX)-RN2O(0,NY,NX)
  trc_solml(idg_H2,0,NY,NX)=trc_solml(idg_H2,0,NY,NX)+trcg_XDFR(idg_H2,NY,NX) &
    +trcs_XFLS(idg_H2,3,0,NY,NX)+GasDisFlx(idg_H2,0,NY,NX)-RH2GO(0,NY,NX)
  trc_solml(idg_NH3,0,NY,NX)=trc_solml(idg_NH3,0,NY,NX)+trcg_XDFR(idg_NH3,NY,NX) &
    +trcs_XFLS(idg_NH3,3,0,NY,NX)+GasDisFlx(idg_NH3,0,NY,NX)+TR_NH3_soil(0,NY,NX)

  trc_solml(ids_NH4,0,NY,NX)=trc_solml(ids_NH4,0,NY,NX)+trcs_XFLS(ids_NH4,3,0,NY,NX) &
    +XNH4S(0,NY,NX)+TR_NH4_soil(0,NY,NX)
  trc_solml(ids_NO3,0,NY,NX)=trc_solml(ids_NO3,0,NY,NX)+trcs_XFLS(ids_NO3,3,0,NY,NX) &
    +XNO3S(0,NY,NX)+TRNO3(0,NY,NX)
  trc_solml(ids_NO2,0,NY,NX)=trc_solml(ids_NO2,0,NY,NX)+trcs_XFLS(ids_NO2,3,0,NY,NX) &
    +XNO2S(0,NY,NX)

  trc_solml(ids_H1PO4,0,NY,NX)=trc_solml(ids_H1PO4,0,NY,NX)+trcs_XFLS(ids_H1PO4,3,0,NY,NX) &
    +XH1PS(0,NY,NX)+TR_H1PO4_soil(0,NY,NX)
  trc_solml(ids_H2PO4,0,NY,NX)=trc_solml(ids_H2PO4,0,NY,NX)+trcs_XFLS(ids_H2PO4,3,0,NY,NX) &
    +XH2PS(0,NY,NX)+TR_H2PO4_soil(0,NY,NX)

!update aqueous concentrations
  DO NTG=idg_beg,idg_end
    trc_solml(NTG,NU(NY,NX),NY,NX)=trc_solml(NTG,NU(NY,NX),NY,NX)+GasSfAtmFlx(NTG,NY,NX)
  ENDDO

  THRE(NY,NX)=THRE(NY,NX)+RCO2O(0,NY,NX)+RCH4O(0,NY,NX)
  UN2GG(NY,NX)=UN2GG(NY,NX)+RN2G(0,NY,NX)
  HN2GG(NY,NX)=HN2GG(NY,NX)+RN2G(0,NY,NX)

  ROXYF(0,NY,NX)=GasDisFlx(idg_O2,0,NY,NX)
  RCO2F(0,NY,NX)=GasDisFlx(idg_CO2,0,NY,NX)
  RCH4F(0,NY,NX)=GasDisFlx(idg_CH4,0,NY,NX)

  !ROXYL:=soil surface O2 dissolution + aqueous O2 flux micropore
  ROXYL(0,NY,NX)=trcg_XDFR(idg_O2,NY,NX)+trcs_XFLS(idg_O2,3,0,NY,NX) &
    -(FLQRQ(NY,NX)*COXR(NY,NX)+FLQRI(NY,NX)*COXQ(NY,NX))
  RCH4L(0,NY,NX)=trcg_XDFR(idg_CH4,NY,NX)+trcs_XFLS(idg_CH4,3,0,NY,NX) &
    -(FLQRQ(NY,NX)*CCHR(NY,NX)+FLQRI(NY,NX)*CCHQ(NY,NX))
  ROXYL(NU(NY,NX),NY,NX)=ROXYL(NU(NY,NX),NY,NX)+GasSfAtmFlx(idg_O2,NY,NX)
  RCH4L(NU(NY,NX),NY,NX)=RCH4L(NU(NY,NX),NY,NX)+GasSfAtmFlx(idg_CH4,NY,NX)
  !
  !     SURFACE LITTER ION EXCHANGE AND PRECIPITATION
  !

  trcx_solml(idx_NH4,0,NY,NX)=trcx_solml(idx_NH4,0,NY,NX)+trcx_TR(idx_NH4,0,NY,NX)
  DO NTX=idx_AEC+1,idx_anion_soil_end
    trcx_solml(NTX,0,NY,NX)=trcx_solml(NTX,0,NY,NX)+trcx_TR(NTX,0,NY,NX)
  ENDDO

  DO NTP=idsp_psoi_beg,idsp_psoi_end
    trcp_salml(NTP,0,NY,NX)=trcp_salml(NTP,0,NY,NX)+trcp_TR(NTP,0,NY,NX)
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
      OQC(K,0,NY,NX)=OQC(K,0,NY,NX)+TOCQRS(K,NY,NX)
      OQN(K,0,NY,NX)=OQN(K,0,NY,NX)+TONQRS(K,NY,NX)
      OQP(K,0,NY,NX)=OQP(K,0,NY,NX)+TOPQRS(K,NY,NX)
      OQA(K,0,NY,NX)=OQA(K,0,NY,NX)+TOAQRS(K,NY,NX)
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

  subroutine CalcLitterLayerChemicalMass(NY,NX)
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
    RC0(K,NY,NX)=RC0(K,NY,NX)+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX) &
      +OHC(K,0,NY,NX)+OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
    RC0(K,NY,NX)=RC0(K,NY,NX)+SUM(OSC(1:jsken,K,0,NY,NX))
  ENDDO

!
!     TOTAL DOC, DON, DOP
!
  DC=DC+SUM(OQC(1:NumOfLitrCmplxs,0,NY,NX))+SUM(OQCH(1:NumOfLitrCmplxs,0,NY,NX)) &
       +SUM(OHC(1:NumOfLitrCmplxs,0,NY,NX))+SUM(OQA(1:NumOfLitrCmplxs,0,NY,NX)) &
       +SUM(OQAH(1:NumOfLitrCmplxs,0,NY,NX))+SUM(OHA(1:NumOfLitrCmplxs,0,NY,NX))

  DN=DN+SUM(OQN(1:NumOfLitrCmplxs,0,NY,NX))+SUM(OQNH(1:NumOfLitrCmplxs,0,NY,NX)) &
       +SUM(OHN(1:NumOfLitrCmplxs,0,NY,NX))
  DP=DP+SUM(OQP(1:NumOfLitrCmplxs,0,NY,NX))+SUM(OQPH(1:NumOfLitrCmplxs,0,NY,NX)) &
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
  UCO2S(NY,NX)=UCO2S(NY,NX)+CS
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
    +trcp_salml(idsp_CaHPO4,0,NY,NX))+2._r8*patomw*trcp_salml(idsp_CaH2PO4,0,NY,NX) &
    +3._r8*patomw*trcp_salml(idsp_HA,0,NY,NX)
  TLPO4=TLPO4+POS+POX+POP
  UPO4(NY,NX)=UPO4(NY,NX)+POX
  UPP4(NY,NX)=UPP4(NY,NX)+POP

  IF(salt_model)call UpdateSurfaceLayerSalt(NY,NX,TLPO4)

  end subroutine CalcLitterLayerChemicalMass
!------------------------------------------------------------------------------------------

  subroutine UpdateSurfaceLayerSalt(NY,NX,TLPO4)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(inout) :: TLPO4
  real(r8) :: SSS,PSS
  INTEGER :: NTSA

  DO NTSA=idsalt_beg,idsalt_end
    trcSalt_solml(NTSA,0,NY,NX)=trcSalt_solml(NTSA,0,NY,NX)+trcSalt_XFLS(NTSA,3,0,NY,NX)
  ENDDO

  PSS=patomw*(trcSalt_solml(idsalt_H0PO4,0,NY,NX)+trcSalt_solml(idsalt_H3PO4,0,NY,NX)&
    +trcSalt_solml(idsalt_FeHPO4,0,NY,NX) &
    +trcSalt_solml(idsalt_FeH2PO4,0,NY,NX)+trcSalt_solml(idsalt_CaPO4,0,NY,NX)&
    +trcSalt_solml(idsalt_CaHPO4,0,NY,NX) &
    +trcSalt_solml(idsalt_CaH2PO4,0,NY,NX)+trcSalt_solml(idsalt_MgHPO4,0,NY,NX))
  TLPO4=TLPO4+PSS
  SSS=trcSalt_solml(idsalt_Al,0,NY,NX)+trcSalt_solml(idsalt_Fe,0,NY,NX) &
    +trcSalt_solml(idsalt_Hp,0,NY,NX)+trcSalt_solml(idsalt_Ca,0,NY,NX) &
    +trcSalt_solml(idsalt_Mg,0,NY,NX)+trcSalt_solml(idsalt_Na,0,NY,NX) &
    +trcSalt_solml(idsalt_K,0,NY,NX)+trcSalt_solml(idsalt_OH,0,NY,NX) &
    +trcSalt_solml(idsalt_SO4,0,NY,NX)+trcSalt_solml(idsalt_Cl,0,NY,NX) &
    +trcSalt_solml(idsalt_CO3,0,NY,NX)+trcSalt_solml(idsalt_H0PO4,0,NY,NX) &
    +2.0*(trcSalt_solml(idsalt_HCO3,0,NY,NX)+trcSalt_solml(idsalt_AlOH,0,NY,NX) &
    +trcSalt_solml(idsalt_AlSO4,0,NY,NX) &
    +trcSalt_solml(idsalt_FeOH,0,NY,NX)+trcSalt_solml(idsalt_FeSO4,0,NY,NX) &
    +trcSalt_solml(idsalt_CaOH2,0,NY,NX)+trcSalt_solml(idsalt_CaCO3,0,NY,NX) &
    +trcSalt_solml(idsalt_CaSO4,0,NY,NX)+trcSalt_solml(idsalt_MgOH2,0,NY,NX) &
    +trcSalt_solml(idsalt_MgCO3,0,NY,NX)+trcSalt_solml(idsalt_MgSO4,0,NY,NX) &
    +trcSalt_solml(idsalt_NaCO3,0,NY,NX)+trcSalt_solml(idsalt_NaSO4,0,NY,NX) &
    +trcSalt_solml(idsalt_KSO4,0,NY,NX)+trcSalt_solml(idsalt_CaPO4,0,NY,NX)) &
    +3.0*(trcSalt_solml(idsalt_AlOH2,0,NY,NX)+trcSalt_solml(idsalt_FeOH2,0,NY,NX) &
    +trcSalt_solml(idsalt_CaHCO3,0,NY,NX) &
    +trcSalt_solml(idsalt_MgHCO3,0,NY,NX)+trcSalt_solml(idsalt_FeHPO4,0,NY,NX) &
    +trcSalt_solml(idsalt_CaHPO4,0,NY,NX)+trcSalt_solml(idsalt_MgHPO4,0,NY,NX)) &
    +4.0*(trcSalt_solml(idsalt_AlOH3,0,NY,NX)+trcSalt_solml(idsalt_FeOH3,0,NY,NX) &
    +trcSalt_solml(idsalt_H3PO4,0,NY,NX) &
    +trcSalt_solml(idsalt_FeH2PO4,0,NY,NX)+trcSalt_solml(idsalt_CaH2PO4,0,NY,NX)) &
    +5.0*(trcSalt_solml(idsalt_AlOH4,0,NY,NX)+trcSalt_solml(idsalt_FeOH4,0,NY,NX))
  TION=TION+SSS
  UION(NY,NX)=UION(NY,NX)+SSS
  end subroutine UpdateSurfaceLayerSalt
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

  !     begin_execution
  !     UPDATE SOIL LAYER VARIABLES WITH TOTAL FLUXES

  D125: DO L=NU(NY,NX),NL(NY,NX)
    !

    UN2GS(NY,NX)=UN2GS(NY,NX)+XN2GS(L,NY,NX)

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
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+TOCFLS(K,L,NY,NX)+XOCFXS(K,L,NY,NX)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+TONFLS(K,L,NY,NX)+XONFXS(K,L,NY,NX)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+TOPFLS(K,L,NY,NX)+XOPFXS(K,L,NY,NX)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+TOAFLS(K,L,NY,NX)+XOAFXS(K,L,NY,NX)
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)+TOCFHS(K,L,NY,NX)-XOCFXS(K,L,NY,NX)
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)+TONFHS(K,L,NY,NX)-XONFXS(K,L,NY,NX)
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)+TOPFHS(K,L,NY,NX)-XOPFXS(K,L,NY,NX)
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)+TOAFHS(K,L,NY,NX)-XOAFXS(K,L,NY,NX)

    ENDDO D8560
    !
    !     DOC, DON, DOP FROM PLANT EXUDATION
    !
    D195: DO K=1,jcplx
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+TDFOME(ielmc,K,L,NY,NX)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+TDFOME(ielmn,K,L,NY,NX)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+TDFOME(ielmp,K,L,NY,NX)
    ENDDO D195
    !
    !     SOIL SOLUTES FROM AQUEOUS TRANSPORT, MICROBIAL AND ROOT
    !     EXCHANGE, EQUILIBRIUM REACTIONS, GAS EXCHANGE,
    !     MICROPORE-MACROPORE EXCHANGE,
    !
    trc_solml(idg_CO2,L,NY,NX)=trc_solml(idg_CO2,L,NY,NX) &
      +trcs_TFLS(idg_CO2,L,NY,NX)+GasDisFlx(idg_CO2,L,NY,NX) &
      -RCO2O(L,NY,NX)-TCO2S(L,NY,NX)+trcs_RFLU(idg_CO2,L,NY,NX) &
      +trcs_XFXS(idg_CO2,L,NY,NX) &
      +TR_CO2_aqu_soil(L,NY,NX)+trcg_XBLL(idg_CO2,L,NY,NX)
    trc_solml(idg_CH4,L,NY,NX)=trc_solml(idg_CH4,L,NY,NX) &
      +trcs_TFLS(idg_CH4,L,NY,NX)+GasDisFlx(idg_CH4,L,NY,NX) &
      -RCH4O(L,NY,NX)-TUPCHS(L,NY,NX)+trcs_RFLU(idg_CH4,L,NY,NX) &
      +trcs_XFXS(idg_CH4,L,NY,NX)+trcg_XBLL(idg_CH4,L,NY,NX)
    trc_solml(idg_O2,L,NY,NX)=trc_solml(idg_O2,L,NY,NX) &
      +trcs_TFLS(idg_O2,L,NY,NX)+GasDisFlx(idg_O2,L,NY,NX) &
      -RUPOXO(L,NY,NX)-TUPOXS(L,NY,NX)+trcs_RFLU(idg_O2,L,NY,NX) &
      +trcs_XFXS(idg_O2,L,NY,NX)+trcg_XBLL(idg_O2,L,NY,NX)

    trc_solml(idg_N2,L,NY,NX)=trc_solml(idg_N2,L,NY,NX) &
      +trcs_TFLS(idg_N2,L,NY,NX)+GasDisFlx(idg_N2,L,NY,NX) &
      -RN2G(L,NY,NX)-TUPNF(L,NY,NX)+trcs_RFLU(idg_N2,L,NY,NX)+trcs_XFXS(idg_N2,L,NY,NX) &
      -XN2GS(L,NY,NX)+trcg_XBLL(idg_N2,L,NY,NX)
    trc_solml(idg_N2O,L,NY,NX)=trc_solml(idg_N2O,L,NY,NX) &
      +trcs_TFLS(idg_N2O,L,NY,NX)+GasDisFlx(idg_N2O,L,NY,NX) &
      -RN2O(L,NY,NX)-TUPN2S(L,NY,NX)+trcs_RFLU(idg_N2O,L,NY,NX)+trcs_XFXS(idg_N2O,L,NY,NX) &
      +trcg_XBLL(idg_N2O,L,NY,NX)

    trc_solml(idg_H2,L,NY,NX)=trc_solml(idg_H2,L,NY,NX) &
      +trcs_TFLS(idg_H2,L,NY,NX)+GasDisFlx(idg_H2,L,NY,NX) &
      -RH2GO(L,NY,NX)-TUPHGS(L,NY,NX)+trcs_RFLU(idg_H2,L,NY,NX) &
      +trcs_XFXS(idg_H2,L,NY,NX)+trcg_XBLL(idg_H2,L,NY,NX)
    trc_solml(idg_NH3,L,NY,NX)=trc_solml(idg_NH3,L,NY,NX) &
      +trcs_TFLS(idg_NH3,L,NY,NX)+GasDisFlx(idg_NH3,L,NY,NX) &
      +TR_NH3_soil(L,NY,NX)-TUPN3S(L,NY,NX)+trcs_RFLU(idg_NH3,L,NY,NX) &
      +trcs_XFXS(idg_NH3,L,NY,NX)+trcg_XBLL(idg_NH3,L,NY,NX)
    trc_solml(ids_NH4,L,NY,NX)=trc_solml(ids_NH4,L,NY,NX) &
      +trcs_TFLS(ids_NH4,L,NY,NX)+XNH4S(L,NY,NX) &
      +TR_NH4_soil(L,NY,NX)-TUPNH4(L,NY,NX)+trcs_RFLU(ids_NH4,L,NY,NX) &
      +trcs_XFXS(ids_NH4,L,NY,NX)

    trc_solml(ids_NO3,L,NY,NX)=trc_solml(ids_NO3,L,NY,NX) &
      +trcs_TFLS(ids_NO3,L,NY,NX)+XNO3S(L,NY,NX) &
      +TRNO3(L,NY,NX)-TUPNO3(L,NY,NX)+trcs_RFLU(ids_NO3,L,NY,NX) &
      +trcs_XFXS(ids_NO3,L,NY,NX)
    trc_solml(ids_NO2,L,NY,NX)=trc_solml(ids_NO2,L,NY,NX) &
      +trcs_TFLS(ids_NO2,L,NY,NX)+XNO2S(L,NY,NX) &
      +TRNO2(L,NY,NX)+trcs_XFXS(ids_NO2,L,NY,NX)

    trc_solml(ids_H1PO4,L,NY,NX)=trc_solml(ids_H1PO4,L,NY,NX) &
      +trcs_TFLS(ids_H1PO4,L,NY,NX)+XH1PS(L,NY,NX) &
      +TR_H1PO4_soil(L,NY,NX)-TUPH1P(L,NY,NX)+trcs_RFLU(ids_H1PO4,L,NY,NX)+trcs_XFXS(ids_H1PO4,L,NY,NX)
    trc_solml(ids_H2PO4,L,NY,NX)=trc_solml(ids_H2PO4,L,NY,NX) &
      +trcs_TFLS(ids_H2PO4,L,NY,NX)+XH2PS(L,NY,NX) &
      +TR_H2PO4_soil(L,NY,NX)-TUPH2P(L,NY,NX)+trcs_RFLU(ids_H2PO4,L,NY,NX)+trcs_XFXS(ids_H2PO4,L,NY,NX)
    trc_solml(idg_NH3B,L,NY,NX)=trc_solml(idg_NH3B,L,NY,NX) &
      +trcs_TFLS(idg_NH3B,L,NY,NX)+GasDisFlx(idg_NH3B,L,NY,NX) &
      +TR_NH3_band_soil(L,NY,NX)-TUPN3B(L,NY,NX)+trcs_RFLU(idg_NH3B,L,NY,NX) &
      +trcs_XFXS(idg_NH3B,L,NY,NX)+trcg_XBLL(idg_NH3B,L,NY,NX)
    trc_solml(ids_NH4B,L,NY,NX)=trc_solml(ids_NH4B,L,NY,NX) &
      +trcs_TFLS(ids_NH4B,L,NY,NX)+XNH4B(L,NY,NX) &
      +TR_NH4_band_soil(L,NY,NX)-TUPNHB(L,NY,NX)+trcs_RFLU(ids_NH4B,L,NY,NX) &
      +trcs_XFXS(ids_NH4B,L,NY,NX)
    trc_solml(ids_NO3B,L,NY,NX)=trc_solml(ids_NO3B,L,NY,NX) &
      +trcs_TFLS(ids_NO3B,L,NY,NX)+XNO3B(L,NY,NX) &
      +TRNOB(L,NY,NX)-TUPNOB(L,NY,NX)+trcs_RFLU(ids_NO3B,L,NY,NX) &
      +trcs_XFXS(ids_NO3B,L,NY,NX)
    trc_solml(ids_NO2B,L,NY,NX)=trc_solml(ids_NO2B,L,NY,NX) &
      +trcs_TFLS(ids_NO2B,L,NY,NX)+XNO2B(L,NY,NX) &
      +TRN2B(L,NY,NX)+trcs_XFXS(ids_NO2B,L,NY,NX)
    trc_solml(ids_H1PO4B,L,NY,NX)=trc_solml(ids_H1PO4B,L,NY,NX) &
      +trcs_TFLS(ids_H1PO4B,L,NY,NX)+XH1BS(L,NY,NX) &
      +TR_H1PO4_band_soil(L,NY,NX)-TUPH1B(L,NY,NX)+trcs_RFLU(ids_H1PO4B,L,NY,NX) &
      +trcs_XFXS(ids_H1PO4B,L,NY,NX)
    trc_solml(ids_H2PO4B,L,NY,NX)=trc_solml(ids_H2PO4B,L,NY,NX) &
      +trcs_TFLS(ids_H2PO4B,L,NY,NX)+XH2BS(L,NY,NX) &
      +TR_H2PO4_band_soil(L,NY,NX) -TUPH2B(L,NY,NX)+trcs_RFLU(ids_H2PO4B,L,NY,NX) &
      +trcs_XFXS(ids_H2PO4B,L,NY,NX)

    THRE(NY,NX)=THRE(NY,NX)+RCO2O(L,NY,NX)+RCH4O(L,NY,NX)
    UN2GG(NY,NX)=UN2GG(NY,NX)+RN2G(L,NY,NX)
    HN2GG(NY,NX)=HN2GG(NY,NX)+RN2G(L,NY,NX)
    !
    !     EXCHANGEABLE CATIONS AND ANIONS FROM EXCHANGE REACTIONS
    !


    trcx_solml(idx_NH4,L,NY,NX)=trcx_solml(idx_NH4,L,NY,NX)+trcx_TR(idx_NH4,L,NY,NX)
    trcx_solml(idx_NH4B,L,NY,NX)=trcx_solml(idx_NH4B,L,NY,NX)+trcx_TR(idx_NH4B,L,NY,NX)

    DO NTX=idx_AEC+1,idx_end
      trcx_solml(NTX,L,NY,NX)=trcx_solml(NTX,L,NY,NX)+trcx_TR(NTX,L,NY,NX)
    ENDDO

    !
    !     PRECIPITATES FROM PRECIPITATION-DISSOLUTION REACTIONS
    !
    DO NTP=idsp_p_beg,idsp_p_end
      trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)+trcp_TR(NTP,L,NY,NX)
    ENDDO
    !
    !     MACROPORE SOLUTES FROM MACROPORE-MICROPORE EXCHANGE
    !
    DO NTS=ids_beg,ids_end
      trc_soHml(NTS,L,NY,NX)=trc_soHml(NTS,L,NY,NX)+trcs_TFHS(NTS,L,NY,NX)-trcs_XFXS(NTS,L,NY,NX)
    ENDDO
    !
    !     GASES FROM VOLATILIZATION-DISSOLUTION AND GAS TRANSFER
!
    DO NTG=idg_beg,idg_end-1
      trc_gasml(NTG,L,NY,NX)=trc_gasml(NTG,L,NY,NX)+RTGasADFlx(NTG,L,NY,NX)-GasDisFlx(NTG,L,NY,NX)
    ENDDO

    trc_gasml(idg_NH3,L,NY,NX)=trc_gasml(idg_NH3,L,NY,NX)-GasDisFlx(idg_NH3B,L,NY,NX)+TRN3G(L,NY,NX)
    ROXYF(L,NY,NX)=RTGasADFlx(idg_O2,L,NY,NX)
    RCO2F(L,NY,NX)=RTGasADFlx(idg_CO2,L,NY,NX)
    RCH4F(L,NY,NX)=RTGasADFlx(idg_CH4,L,NY,NX)
    ROXYL(L,NY,NX)=trcs_TFLS(idg_O2,L,NY,NX)+trcs_RFLU(idg_O2,L,NY,NX) &
      +trcs_XFXS(idg_O2,L,NY,NX)+trcg_XBLL(idg_O2,L,NY,NX)
    RCH4L(L,NY,NX)=trcs_TFLS(idg_CH4,L,NY,NX)+trcs_RFLU(idg_CH4,L,NY,NX) &
      +trcs_XFXS(idg_CH4,L,NY,NX)+trcg_XBLL(idg_CH4,L,NY,NX)
    !
    !     GRID CELL BOUNDARY FLUXES FROM ROOT GAS TRANSFER
!   watch out the following code for changes
    HEATIN=HEATIN+THeatSoiThaw(L,NY,NX)+THeatRootUptake(L,NY,NX)
    CIB=trcg_TFLA(idg_CO2,L,NY,NX)
    CHB=trcg_TFLA(idg_CH4,L,NY,NX)
    OIB=trcg_TFLA(idg_O2,L,NY,NX)
!    HGB=trcg_TFLA(idg_H2,L,NY,NX)
    HGB=0.0_r8
    ZGB=0.0_r8
    Z2B=trcg_TFLA(idg_N2O,L,NY,NX)
    ZHB=trcg_TFLA(idg_NH3,L,NY,NX)
!
!     GRID CELL BOUNDARY FLUXES BUBBLING
!
    IF(LG.EQ.0)THEN
      LL=0
      CIB=CIB+trcg_XBLL(idg_CO2,L,NY,NX)
      CHB=CHB+trcg_XBLL(idg_CH4,L,NY,NX)
      OIB=OIB+trcg_XBLL(idg_O2,L,NY,NX)
      ZGB=ZGB+trcg_XBLL(idg_N2,L,NY,NX)
      Z2B=Z2B+trcg_XBLL(idg_N2O,L,NY,NX)
      ZHB=ZHB+trcg_XBLL(idg_NH3,L,NY,NX)+trcg_XBLL(idg_NH3B,L,NY,NX)
      HGB=HGB+trcg_XBLL(idg_H2,L,NY,NX)
    ELSE
      LL=MIN(L,LG)
      DO NTG=idg_beg,idg_end-1
        trc_gasml(NTG,LL,NY,NX)=trc_gasml(NTG,LL,NY,NX)-trcg_XBLL(NTG,L,NY,NX)
      ENDDO
      trc_gasml(idg_NH3,LL,NY,NX)=trc_gasml(idg_NH3,LL,NY,NX)-trcg_XBLL(idg_NH3B,L,NY,NX)

      IF(LG.LT.L)THEN
        TLCO2G=TLCO2G-trcg_XBLL(idg_CO2,L,NY,NX)-trcg_XBLL(idg_CH4,L,NY,NX)
        UCO2S(NY,NX)=UCO2S(NY,NX)-trcg_XBLL(idg_CO2,L,NY,NX)-trcg_XBLL(idg_CH4,L,NY,NX)
        OXYGSO=OXYGSO-trcg_XBLL(idg_O2,L,NY,NX)
        TLN2G=TLN2G-trcg_XBLL(idg_N2,L,NY,NX)-trcg_XBLL(idg_N2O,L,NY,NX) &
          -trcg_XBLL(idg_NH3,L,NY,NX)-trcg_XBLL(idg_NH3B,L,NY,NX)
        TLH2G=TLH2G-trcg_XBLL(idg_H2,L,NY,NX)
      ENDIF
    ENDIF
    CO2GIN=CO2GIN+CIB+CHB
    COB=TCO2P(L,NY,NX)+TCO2S(L,NY,NX)-TR_CO2_aqu_soil(L,NY,NX)
    TCOU=TCOU+COB
    HCO2G(NY,NX)=HCO2G(NY,NX)+CIB
    UCO2G(NY,NX)=UCO2G(NY,NX)+CIB
    HCH4G(NY,NX)=HCH4G(NY,NX)+CHB
    UCH4G(NY,NX)=UCH4G(NY,NX)+CHB
    UCOP(NY,NX)=UCOP(NY,NX)+TCO2P(L,NY,NX)+TCO2S(L,NY,NX)
    UDICD(NY,NX)=UDICD(NY,NX)-12.0*TBCO2(L,NY,NX)
    TXCO2(NY,NX)=TXCO2(NY,NX)+12.0*TBCO2(L,NY,NX)
    OXYGIN=OXYGIN+OIB
    OOB=RUPOXO(L,NY,NX)+TUPOXP(L,NY,NX)+TUPOXS(L,NY,NX)
    OXYGOU=OXYGOU+OOB
    UOXYG(NY,NX)=UOXYG(NY,NX)+OIB
    HOXYG(NY,NX)=HOXYG(NY,NX)+OIB
    H2GIN=H2GIN+HGB
    HOB=RH2GO(L,NY,NX)+TUPHGS(L,NY,NX)
    H2GOU=H2GOU+HOB
    ZN2GIN=ZN2GIN+ZGB+Z2B+ZHB
!     UN2GG(NY,NX)=UN2GG(NY,NX)+ZGB
!     HN2GG(NY,NX)=HN2GG(NY,NX)+ZGB
    UN2OG(NY,NX)=UN2OG(NY,NX)+Z2B
    HN2OG(NY,NX)=HN2OG(NY,NX)+Z2B
    UNH3G(NY,NX)=UNH3G(NY,NX)+ZHB
    HNH3G(NY,NX)=HNH3G(NY,NX)+ZHB
    UH2GG(NY,NX)=UH2GG(NY,NX)+HGB
    !
    !     GRID CELL BOUNDARY FLUXES FROM EQUILIBRIUM REACTIONS
!
    SNM=(2.0_r8*(XNH4S(L,NY,NX)+XNH4B(L,NY,NX)-TUPNH4(L,NY,NX) &
      -TUPNHB(L,NY,NX)-XN2GS(L,NY,NX))-TUPN3S(L,NY,NX)-TUPN3B(L,NY,NX) &
      +XNO3S(L,NY,NX)+XNO3B(L,NY,NX)-TUPNO3(L,NY,NX)-TUPNOB(L,NY,NX) &
      +XNO2S(L,NY,NX)+XNO2B(L,NY,NX))/natomw
    SPM=(2.0_r8*(XH1PS(L,NY,NX)+XH1BS(L,NY,NX)-TUPH1P(L,NY,NX) &
      -TUPH1B(L,NY,NX))+3.0*(XH2PS(L,NY,NX)+XH2BS(L,NY,NX) &
      -TUPH2P(L,NY,NX)-TUPH2B(L,NY,NX)))/patomw
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
    UCO2S(NY,NX)=UCO2S(NY,NX)+CS
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

    POX=patomw*(trcx_solml(idx_HPO4,L,NY,NX)+trcx_solml(idx_H2PO4,L,NY,NX) &
      +trcx_solml(idx_HPO4B,L,NY,NX)+trcx_solml(idx_H2PO4B,L,NY,NX))
    POP=patomw*(trcp_salml(idsp_AlPO4,L,NY,NX)+trcp_salml(idsp_FePO4,L,NY,NX)&
      +trcp_salml(idsp_CaHPO4,L,NY,NX)+trcp_salml(idsp_AlPO4B,L,NY,NX)&
      +trcp_salml(idsp_FePO4B,L,NY,NX)+trcp_salml(idsp_CaHPO4B,L,NY,NX)) &
      +2._r8*patomw*(trcp_salml(idsp_CaH2PO4,L,NY,NX)+trcp_salml(idsp_CaH2PO4B,L,NY,NX)) &
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
  integer  :: NTP,NTSA
  real(r8) :: ECHY,ECOH,ECAL,ECFE,ECCA,ECMG,ECNA,ECKA,ECCO,ECHC
  real(r8) :: ECNO,ECSO,ECCL
  real(r8) :: PSS,SSS,SSH,SSF,SSX,SST,SSP

  DO NTSA=idsalt_beg,idsaltb_end
    trcSalt_solml(NTSA,L,NY,NX)=trcSalt_solml(NTSA,L,NY,NX)+trcSalt_TR(NTSA,L,NY,NX) &
      +trcSalt_TFLS(NTSA,L,NY,NX)+trcSalt_RFLU(NTSA,L,NY,NX)+trcSalt_XFXS(NTSA,L,NY,NX)

    trcSalt_soHml(NTSA,L,NY,NX)=trcSalt_soHml(NTSA,L,NY,NX)+trcSalt_TFHS(NTSA,L,NY,NX) &
      -trcSalt_XFXS(NTSA,L,NY,NX)
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
    trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)+trcp_TR(NTP,L,NY,NX)
  ENDDO

  PSS=patomw*(trcSalt_solml(idsalt_H0PO4,L,NY,NX)+trcSalt_solml(idsalt_H3PO4,L,NY,NX) &
    +trcSalt_solml(idsalt_FeHPO4,L,NY,NX) &
    +trcSalt_solml(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml(idsalt_CaPO4,L,NY,NX) &
    +trcSalt_solml(idsalt_CaHPO4,L,NY,NX) &
    +trcSalt_solml(idsalt_CaH2PO4,L,NY,NX)+trcSalt_solml(idsalt_MgHPO4,L,NY,NX) &
    +trcSalt_solml(idsalt_H0PO4B,L,NY,NX) &
    +trcSalt_solml(idsalt_H3PO4B,L,NY,NX)+trcSalt_solml(idsalt_FeHPO4B,L,NY,NX) &
    +trcSalt_solml(idsalt_FeH2PO4B,L,NY,NX) &
    +trcSalt_solml(idsalt_CaPO4B,L,NY,NX)+trcSalt_solml(idsalt_CaHPO4B,L,NY,NX) &
    +trcSalt_solml(idsalt_CaH2PO4B,L,NY,NX) &
    +trcSalt_solml(idsalt_MgHPO4B,L,NY,NX)+trcSalt_soHml(idsalt_H0PO4,L,NY,NX)&
    +trcSalt_soHml(idsalt_H3PO4,L,NY,NX)+trcSalt_soHml(idsalt_FeHPO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_FeH2PO4,L,NY,NX)+trcSalt_soHml(idsalt_CaPO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_CaHPO4,L,NY,NX)+trcSalt_soHml(idsalt_CaH2PO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_MgHPO4,L,NY,NX)+trcSalt_soHml(idsalt_H0PO4B,L,NY,NX) &
    +trcSalt_soHml(idsalt_H3PO4B,L,NY,NX)+trcSalt_soHml(idsalt_FeHPO4B,L,NY,NX) &
    +trcSalt_soHml(idsalt_FeH2PO4B,L,NY,NX)+trcSalt_soHml(idsalt_CaPO4B,L,NY,NX)&
    +trcSalt_soHml(idsalt_CaHPO4B,L,NY,NX)+trcSalt_soHml(idsalt_CaH2PO4B,L,NY,NX) &
    +trcSalt_soHml(idsalt_MgHPO4B,L,NY,NX))
  TLPO4=TLPO4+PSS

  SSS=trcSalt_solml(idsalt_Al,L,NY,NX)+trcSalt_solml(idsalt_Fe,L,NY,NX) &
    +trcSalt_solml(idsalt_Hp,L,NY,NX)+trcSalt_solml(idsalt_Ca,L,NY,NX) &
    +trcSalt_solml(idsalt_Mg,L,NY,NX)+trcSalt_solml(idsalt_Na,L,NY,NX) &
    +trcSalt_solml(idsalt_K,L,NY,NX)+trcSalt_solml(idsalt_OH,L,NY,NX) &
    +trcSalt_solml(idsalt_SO4,L,NY,NX)+trcSalt_solml(idsalt_Cl,L,NY,NX) &
    +trcSalt_solml(idsalt_CO3,L,NY,NX)+trcSalt_solml(idsalt_H0PO4,L,NY,NX) &
    +trcSalt_solml(idsalt_H0PO4B,L,NY,NX) &
    +2.0_r8*(trcSalt_solml(idsalt_HCO3,L,NY,NX)+trcSalt_solml(idsalt_AlOH,L,NY,NX) &
    +trcSalt_solml(idsalt_AlSO4,L,NY,NX)+trcSalt_solml(idsalt_FeOH,L,NY,NX)&
    +trcSalt_solml(idsalt_FeSO4,L,NY,NX)+trcSalt_solml(idsalt_CaOH2,L,NY,NX) &
    +trcSalt_solml(idsalt_CaCO3,L,NY,NX)+trcSalt_solml(idsalt_CaSO4,L,NY,NX)&
    +trcSalt_solml(idsalt_MgOH2,L,NY,NX)+trcSalt_solml(idsalt_MgCO3,L,NY,NX) &
    +trcSalt_solml(idsalt_MgSO4,L,NY,NX)+trcSalt_solml(idsalt_NaCO3,L,NY,NX)&
    +trcSalt_solml(idsalt_NaSO4,L,NY,NX)+trcSalt_solml(idsalt_KSO4,L,NY,NX) &
    +trcSalt_solml(idsalt_CaPO4,L,NY,NX)+trcSalt_solml(idsalt_CaPO4B,L,NY,NX)) &
    +3.0_r8*(trcSalt_solml(idsalt_AlOH2,L,NY,NX)+trcSalt_solml(idsalt_FeOH2,L,NY,NX)&
    +trcSalt_solml(idsalt_CaHCO3,L,NY,NX) &
    +trcSalt_solml(idsalt_MgHCO3,L,NY,NX)+trcSalt_solml(idsalt_FeHPO4,L,NY,NX)&
    +trcSalt_solml(idsalt_CaHPO4,L,NY,NX)+trcSalt_solml(idsalt_MgHPO4,L,NY,NX) &
    +trcSalt_solml(idsalt_FeHPO4B,L,NY,NX)+trcSalt_solml(idsalt_CaHPO4B,L,NY,NX)&
    +trcSalt_solml(idsalt_MgHPO4B,L,NY,NX)) &
    +4.0_r8*(trcSalt_solml(idsalt_AlOH3,L,NY,NX)+trcSalt_solml(idsalt_FeOH3,L,NY,NX) &
    +trcSalt_solml(idsalt_H3PO4,L,NY,NX) &
    +trcSalt_solml(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml(idsalt_CaH2PO4,L,NY,NX) &
    +trcSalt_solml(idsalt_H3PO4B,L,NY,NX)+trcSalt_solml(idsalt_FeH2PO4B,L,NY,NX) &
    +trcSalt_solml(idsalt_CaH2PO4B,L,NY,NX)) &
    +5.0_r8*(trcSalt_solml(idsalt_AlOH4,L,NY,NX)+trcSalt_solml(idsalt_FeOH4,L,NY,NX))

  SSH=trcSalt_soHml(idsalt_Al,L,NY,NX)+trcSalt_soHml(idsalt_Fe,L,NY,NX)&
    +trcSalt_soHml(idsalt_Hp,L,NY,NX)+trcSalt_soHml(idsalt_Ca,L,NY,NX) &
    +trcSalt_soHml(idsalt_Mg,L,NY,NX)+trcSalt_soHml(idsalt_Na,L,NY,NX) &
    +trcSalt_soHml(idsalt_K,L,NY,NX)+trcSalt_soHml(idsalt_OH,L,NY,NX)  &
    +trcSalt_soHml(idsalt_SO4,L,NY,NX)+trcSalt_soHml(idsalt_Cl,L,NY,NX) &
    +trcSalt_soHml(idsalt_CO3,L,NY,NX) +trcSalt_soHml(idsalt_H0PO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_H0PO4B,L,NY,NX) &
    +2.0_r8*(trcSalt_soHml(idsalt_HCO3,L,NY,NX)+trcSalt_soHml(idsalt_AlOH,L,NY,NX) &
    +trcSalt_soHml(idsalt_AlSO4,L,NY,NX)+trcSalt_soHml(idsalt_FeOH,L,NY,NX) &
    +trcSalt_soHml(idsalt_FeSO4,L,NY,NX)+trcSalt_soHml(idsalt_CaOH2,L,NY,NX) &
    +trcSalt_soHml(idsalt_CaCO3,L,NY,NX)+trcSalt_soHml(idsalt_CaSO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_MgOH2,L,NY,NX)+trcSalt_soHml(idsalt_MgCO3,L,NY,NX) &
    +trcSalt_soHml(idsalt_MgSO4,L,NY,NX)+trcSalt_soHml(idsalt_NaCO3,L,NY,NX) &
    +trcSalt_soHml(idsalt_NaSO4,L,NY,NX)+trcSalt_soHml(idsalt_KSO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_CaPO4,L,NY,NX)+trcSalt_soHml(idsalt_CaPO4B,L,NY,NX)) &
    +3.0_r8*(trcSalt_soHml(idsalt_AlOH2,L,NY,NX)+trcSalt_soHml(idsalt_FeOH2,L,NY,NX) &
    +trcSalt_soHml(idsalt_CaHCO3,L,NY,NX) &
    +trcSalt_soHml(idsalt_MgHCO3,L,NY,NX)+trcSalt_soHml(idsalt_FeHPO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_CaHPO4,L,NY,NX)+trcSalt_soHml(idsalt_MgHPO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_FeHPO4B,L,NY,NX)+trcSalt_soHml(idsalt_CaHPO4B,L,NY,NX) &
    +trcSalt_soHml(idsalt_MgHPO4B,L,NY,NX)) &
    +4.0_r8*(trcSalt_soHml(idsalt_AlOH3,L,NY,NX)+trcSalt_soHml(idsalt_FeOH3,L,NY,NX) &
    +trcSalt_soHml(idsalt_H3PO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_FeH2PO4,L,NY,NX)+trcSalt_soHml(idsalt_CaH2PO4,L,NY,NX) &
    +trcSalt_soHml(idsalt_H3PO4B,L,NY,NX) &
    +trcSalt_soHml(idsalt_FeH2PO4B,L,NY,NX)+trcSalt_soHml(idsalt_CaH2PO4B,L,NY,NX)) &
    +5.0_r8*(trcSalt_soHml(idsalt_AlOH4,L,NY,NX)+trcSalt_soHml(idsalt_FeOH4,L,NY,NX))
!
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
    +7.0_r8*(trcp_salml(idsp_CaH2PO4,L,NY,NX)+trcp_salml(idsp_CaH2PO4B,L,NY,NX)) &
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
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
        +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      DP=DP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)

      DC=DC+SUM(OSC(1:jsken,K,L,NY,NX))
      DN=DN+SUM(OSN(1:jsken,K,L,NY,NX))
      DP=DP+SUM(OSP(1:jsken,K,L,NY,NX))

    ELSE
      OC=OC+SUM(ORC(1:ndbiomcp,K,L,NY,NX))
      ON=ON+SUM(ORN(1:ndbiomcp,K,L,NY,NX))
      OP=OP+SUM(ORP(1:ndbiomcp,K,L,NY,NX))

      OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
        +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      ON=ON+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      OP=OP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)

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
      RAINR=LitrfalChemElemnts_vr(ielmc,M,K,0,NY,NX)*THETCX(K)
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
