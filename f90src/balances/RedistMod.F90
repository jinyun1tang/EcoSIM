module RedistMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use abortutils,        only: padr,     print_info, endrun, destroy
  use minimathmod,       only: safe_adb, AZMAX1, fixEXflux, fixnegmass
  use EcoSiMParDataMod,  only: micpar
  use SurfLitterPhysMod, only: UpdateLitRPhys
  use SoilBGCNLayMod  
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
  use SnowTransportMod
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
  real(r8) :: DORGC(JZ),DVLiceMicP_vr(JZ)
  real(r8) :: TXCO2(JY,JX),DORGE(JY,JX)
  real(r8) :: VOLISO,VOLPT,VOLTT
  real(r8) :: TFLWT,orgm(1:NumPlantChemElms)
!     execution begins here
  curday=I
  curhour=J
  VOLISO = 0.0_r8
  TFLWT  = 0.0_r8
  VOLPT  = 0.0_r8
  VOLTT  = 0.0_r8
  
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      TXCO2(NY,NX)        = 0.0_r8
      DORGE(NY,NX)        = 0.0_r8
      call AddFlux2SurfaceResidue(I,J,NY,NX)
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

      call XGridTranspt(I,J,NY,NX,LG)

      call SnowMassUpdate(I,J,NY,NX)

      call HandleSurfaceBoundary(I,J,NY,NX)

!     OVERLAND FLOW
      call OverlandFlow(I,J,NY,NX)

      call SoilErosion(NY,NX,DORGE)
!
      call DiagSnowChemMass(I,J,NY,NX,HeatStore_col(NY,NX))
!
      call LitterLayerSummary(NY,NX)

      call update_physVar_Profile(I,J,NY,NX,VOLISO,DVLiceMicP_vr)    

      call UpdateChemInSoilLays(I,J,NY,NX,LG,DORGC,TXCO2,DORGE)
!
!     SNOWPACK LAYERING

      call SnowpackLayering(I,J,NY,NX)

      call RelayerSoilProfile(NY,NX,DORGC,DVLiceMicP_vr)

      call UpdateOutputVars(I,J,NY,NX,TXCO2)

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
  real(r8) :: VLSoilPoreMicPX   !maximal soil micropore allowed
  integer  :: L
  if(lverb)write(*,*)'UpdateOutputVars'

  Eco_NetRad_col(NY,NX)         = Eco_NetRad_col(NY,NX)+HeatByRadiation_col(NY,NX)
  Eco_Heat_Latent_col(NY,NX)    = Eco_Heat_Latent_col(NY,NX)+HeatEvapAir2Surf_col(NY,NX)
  Eco_Heat_Sens_col(NY,NX)      = Eco_Heat_Sens_col(NY,NX)+HeatSensAir2Surf_col(NY,NX)
  Eco_Heat_Grnd_col(NY,NX)      = Eco_Heat_Grnd_col(NY,NX)-(HeatNet2Surf_col(NY,NX)-HeatSensVapAir2Surf_col(NY,NX))
  Canopy_Heat_Latent_col(NY,NX) = Canopy_Heat_Latent_col(NY,NX)+HeatEvapAir2Surf_col(NY,NX)*BndlResistCanopy_col(NY,NX)
  Canopy_Heat_Sens_col(NY,NX)   = Canopy_Heat_Sens_col(NY,NX)+HeatSensAir2Surf_col(NY,NX)*BndlResistCanopy_col(NY,NX)
  Eco_NEE_col(NY,NX)            = Canopy_NEE_col(NY,NX)+SurfGasFlx_col(idg_CO2,NY,NX)
  ECO_ER_col(NY,NX)             = ECO_ER_col(NY,NX)+SurfGasFlx_col(idg_CO2,NY,NX)
  Eco_NPP_CumYr_col(NY,NX)      = Eco_GPP_CumYr_col(NY,NX)+Eco_AutoR_CumYr_col(NY,NX)
  Eco_NBP_CumYr_col(NY,NX)      = Eco_NBP_CumYr_col(NY,NX)+Canopy_NEE_col(NY,NX) &
    +SurfGasFlx_col(idg_CO2,NY,NX)+SurfGasFlx_col(idg_CH4,NY,NX) &
    +TXCO2(NY,NX)-HydroSufDOCFlx_col(NY,NX)-HydroSufDICFlx_col(NY,NX)-HydroSubsDOCFlx_col(NY,NX)-HydroSubsDICFlx_col(NY,NX)
    
  IF(NU(NY,NX).GT.NUI(NY,NX))THEN  !the surface is lowered
    DO L=NUI(NY,NX),NU(NY,NX)-1
      IF(VLSoilPoreMicP_vr(L,NY,NX).LE.ZEROS2(NY,NX))THEN
        !set default values
        TKS_vr(L,NY,NX) = TKS_vr(NU(NY,NX),NY,NX)
        TCS(L,NY,NX)    = units%Kelvin2Celcius(TKS_vr(L,NY,NX))
      ENDIF
    ENDDO
  ENDIF
  !     MIX ALL SOIL STATE VARIABLES AND INCORPORATE ALL SURFACE
  !     RESIDUE STATE VARIABLES WITHIN THE TILLAGE ZONE TO THE EXTENT
  !     ASSOCIATED IN 'DAY' WITH EACH TILLAGE EVENT ENTERED IN THE
  !     TILLAGE FILE
  if(lverb)write(*,*)'tillage'
  call ApplyTillageMixing(I,J,NY,NX)

  !
  !     OUTPUT FOR SOIL WATER, ICE CONTENTS
  !
  ThetaH2OZ_vr(0,NY,NX)=AZMAX1((VLWatMicP_vr(0,NY,NX)-VWatLitRHoldCapcity_col(NY,NX))/AREA(3,0,NY,NX))
  ThetaICEZ_vr(0,NY,NX)=AZMAX1((VLiceMicP_vr(0,NY,NX)-VWatLitRHoldCapcity_col(NY,NX))/AREA(3,0,NY,NX))
  
  D9945: DO L=NUI(NY,NX),NL(NY,NX)
    VLWatMacP_vr(L,NY,NX) = AZMAX1(VLWatMacP_vr(L,NY,NX))
    VLiceMacP_vr(L,NY,NX) = AZMAX1(VLiceMacP_vr(L,NY,NX))
    VLSoilPoreMicPX       = AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
    VOLTX_vr(L,NY,NX)     = VLSoilPoreMicPX+VLMacP_vr(L,NY,NX)
    ThetaH2OZ_vr(L,NY,NX) = safe_adb(VLWatMicP_vr(L,NY,NX)+AMIN1(VLMacP_vr(L,NY,NX),VLWatMacP_vr(L,NY,NX)),VOLTX_vr(L,NY,NX))/POROS_vr(L,NY,NX)
    ThetaICEZ_vr(L,NY,NX) = safe_adb(VLiceMicP_vr(L,NY,NX)+AMIN1(VLMacP_vr(L,NY,NX),VLiceMacP_vr(L,NY,NX)),VOLTX_vr(L,NY,NX))/POROS_vr(L,NY,NX)
  ENDDO D9945

  end subroutine UpdateOutputVars


!------------------------------------------------------------------------------------------

  subroutine ModifyExWTBLByDisturbance(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8) :: DCORPW
!     begin_execution
!     SolarNoonHour_col=hour of solar noon from weather file
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     DCORP=mixing intensity (fire) or depth (tillage,drainage) of disturbance
!     CumDepz2LayerBot_vr(NU=soil surface elevation
!     DTBLI,WtblDepzTile_col=depth of natural,artificial water table from readi.f
!     ExtWaterTable,ExtWaterTablet0=current,initial natural water table depth
!     DTBLY,DTBLD=current,initial artificial water table depth
!     IDWaterTable=water table flag from readi.f
!        :0=none
!        :1,2=natural stationary,mobile
!        :3,4=artificial stationary,mobile
!     QDischar_col=hourly water loss through lateral and lower boundaries
!
  if(lverb)write(*,*)'ModifyExWTBLByDisturbance'
  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.iSoilDisturbType_col(I,NY,NX).EQ.23)THEN
    ! drainage is on
    DCORPW=DepzCorp_col(I,NY,NX)+CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
    NatWtblDepz_col(NY,NX)=DCORPW
    ExtWaterTablet0(NY,NX)=NatWtblDepz_col(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX))
    ExtWaterTable_col(NY,NX)=ExtWaterTablet0(NY,NX)+CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
  ENDIF

  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.iSoilDisturbType_col(I,NY,NX).EQ.24)THEN
    ! drainage in on
    DCORPW=DepzCorp_col(I,NY,NX)+CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
    IF(IDWaterTable(NY,NX).EQ.1)THEN
      IDWaterTable(NY,NX)=3
    ELSEIF(IDWaterTable(NY,NX).EQ.2)THEN
      IDWaterTable(NY,NX)=4
    ENDIF
    WtblDepzTile_col(NY,NX)=DCORPW
    DTBLD(NY,NX)=AZMAX1(WtblDepzTile_col(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX)))
    DTBLY(NY,NX)=DTBLD(NY,NX)
  ENDIF
!
!     SET DEPTH OF MOBILE EXTERNAL WATER TABLE
! switched on for change of water table due to discharge/drainage
! why 0.00167, time relaxization constant, ?
! 4 is mobile tile drainge.
  IF(IDWaterTable(NY,NX).EQ.2 .OR. IDWaterTable(NY,NX).EQ.4)THEN
    ExtWaterTable_col(NY,NX)=ExtWaterTable_col(NY,NX)-QDischar_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX) &
      -0.00167_r8*(ExtWaterTable_col(NY,NX)-ExtWaterTablet0(NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX))
    ExtWaterTable_col(NY,NX)=ExtWaterTablet0(NY,NX)+CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
  ENDIF
  
  IF(IDWaterTable(NY,NX).EQ.4)THEN
    DTBLY(NY,NX)=DTBLY(NY,NX)-QDischar_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX) &
      -0.00167_r8*(DTBLY(NY,NX)-DTBLD(NY,NX))
  ENDIF
  end subroutine ModifyExWTBLByDisturbance
!------------------------------------------------------------------------------------------

  subroutine HandleSurfaceBoundary(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX

  integer :: L,K,LS,NTG,NTP,NTX,nsalts,NE
  real(r8):: vhcp1s
  real(r8) :: CI,CH,CO,CX
  real(r8) :: OI,OO
  real(r8) :: HI,HO
  real(r8) :: PI,PXB
  real(r8) :: SIN,SGN,SIP,SNB
  real(r8) :: SPB,SNM0,SPM0,SIR,SII,SBU
  real(r8) :: WI,WO,dval,dflx,rval,pval
  real(r8) :: ZSI,ZXB,ZGI
  real(r8) :: ZNGGIN,ZN2OIN,ZNH3IN
  integer :: idom,ids,idg
  ! begin_execution

  if(lverb)write(*,*)'HandleSurfaceBoundary'

  call UpdateLitRPhys(I,J,NY,NX,TXGridSurfRunoff_2DH(NY,NX),THeatXGridBySurfRunoff_2DH(NY,NX),&
    HeatStore_lnd,HEATIN_lnd)

  do idg=idg_beg,idg_NH3-1
    trc_solml_vr(idg,0,NY,NX)=trc_solml_vr(idg,0,NY,NX)+trcg_surf_disevap_flx(idg,NY,NX) &
      +trcs_Transp2MicP_3D(idg,3,0,NY,NX)+Gas_Disol_Flx_vr(idg,0,NY,NX)-trcg_RMicbTransf_vr(idg,0,NY,NX)
    trc_solml_vr(idg,0,NY,NX)=fixnegmass(trc_solml_vr(idg,0,NY,NX))  
  enddo

  trc_solml_vr(idg_N2,0,NY,NX)=trc_solml_vr(idg_N2,0,NY,NX)-Micb_N2Fixation_vr(0,NY,NX)
  trc_solml_vr(idg_N2,0,NY,NX)=fixnegmass(trc_solml_vr(idg_N2,0,NY,NX))
  rval=trc_solml_vr(idg_NH3,0,NY,NX)
  
  dflx=trcg_surf_disevap_flx(idg_NH3,NY,NX)+trcs_Transp2MicP_3D(idg_NH3,3,0,NY,NX) &
    +Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)+TR_NH3_soil_vr(0,NY,NX)

  trc_solml_vr(idg_NH3,0,NY,NX)=rval+dflx  
  !negative value correction  
  
  if(trc_solml_vr(idg_NH3,0,NY,NX)<-1.e-7)then
    pval=abs(rval/dflx)
    trcg_surf_disevap_flx(idg_NH3,NY,NX)= trcg_surf_disevap_flx(idg_NH3,NY,NX)*pval
    trcs_Transp2MicP_3D(idg_NH3,3,0,NY,NX)=trcs_Transp2MicP_3D(idg_NH3,3,0,NY,NX)*pval
    Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)=Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)*pval
    TR_NH3_soil_vr(0,NY,NX)=TR_NH3_soil_vr(0,NY,NX)*pval
    trc_solml_vr(idg_NH3,0,NY,NX)=0._r8
  endif
  trc_solml_vr(idg_NH3,0,NY,NX)=fixnegmass(trc_solml_vr(idg_NH3,0,NY,NX))

  ! GAS EXCHANGE FROM SURFACE VOLATILIZATION-DISSOLUTION
  !

  do ids=ids_nut_beg,ids_nuts_end  
    trc_solml_vr(ids,0,NY,NX)=trc_solml_vr(ids,0,NY,NX)  &
      +RNutMicbTransf_vr(ids,0,NY,NX)+trcn_RChem_soil_vr(ids,0,NY,NX)

    trc_solml_vr(ids,0,NY,NX)=fixnegmass(trc_solml_vr(ids,0,NY,NX))  
    
    dval=-trcs_Transp2MicP_3D(ids,3,0,NY,NX)
    call fixEXflux(trc_solml_vr(ids,0,NY,NX),dval)
    trcs_Transp2MicP_3D(ids,3,0,NY,NX)=-dval
  enddo  

!update aqueous concentrations
  DO NTG=idg_beg,idg_end
    dval=-Gas_Flx_atmDif2soil_col(NTG,NY,NX)
    trc_solml_vr(NTG,NU(NY,NX),NY,NX)=trc_solml_vr(NTG,NU(NY,NX),NY,NX)+Gas_Flx_atmDif2soil_col(NTG,NY,NX)
    call fixEXflux(trc_solml_vr(NTG,NU(NY,NX),NY,NX),dval)
    Gas_Flx_atmDif2soil_col(NTG,NY,NX)=-dval
    if(trc_solml_vr(NTG,NU(NY,NX),NY,NX)<0._r8)then
      write(*,*)'redistntg ',trcs_names(NTG),trc_solml_vr(NTG,NU(NY,NX),NY,NX),Gas_Flx_atmDif2soil_col(NTG,NY,NX)
    endif
  ENDDO

  !
  !     SURFACE BOUNDARY WATER FLUXES
  !
  WI                       = PrecAtm_col(NY,NX)+IrrigSurface_col(NY,NX)   !total incoming water flux    = rain/snowfall + irrigation
  CRAIN                    = CRAIN+WI
  QRain_CumYr_col(NY,NX)   = WI
  WO                       = VapXAir2GSurf_col(NY,NX)+QvET_col(NY,NX)        !total outgoing water flux, > 0 into ground surface
  CEVAP                    = CEVAP-WO
  QEvap_CumYr_col(NY,NX)   = QEvap_CumYr_col(NY,NX)-WO         !>0 into atmosphere
  EvapoTransp_col(NY,NX)   = -WO
  QDischar_col(NY,NX)      = QDischar_col(NY,NX)-IrrigSubsurf_col(NY,NX)
  H2OLoss_CumYr_col(NY,NX) = H2OLoss_CumYr_col(NY,NX)-IrrigSubsurf_col(NY,NX)
  QDrain_col(NY,NX)        = WaterFlowSoiMicP_3D(3,NK(NY,NX),NY,NX)
  HeatDrain_col(NY,NX)     = HeatFlow2Soil_3D(3,NK(NY,NX),NY,NX)
  QH2OLoss_lnds            = QH2OLoss_lnds-IrrigSubsurf_col(NY,NX)
  !
  !     SURFACE BOUNDARY HEAT FLUXES
  !
  HEATIN_lnd=HEATIN_lnd+cpw*TairK_col(NY,NX)*PrecRainAndIrrig_col(NY,NX)+cps*TairK_col(NY,NX)*SnoFalPrec_col(NY,NX)
  HEATIN_lnd=HEATIN_lnd+HeatNet2Surf_col(NY,NX)+HeatFlx2Canopy_col(NY,NX)

  D5150: DO L=1,JS
    HEATIN_lnd=HEATIN_lnd+XPhaseChangeHeatL_snvr(L,NY,NX)
  ENDDO D5150
  HeatOut_lnds=HeatOut_lnds-cpw*TairK_col(NY,NX)*IrrigSubsurf_col(NY,NX)
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
  CI=Gas_Flx_atmDif2soil_col(idg_CO2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_CO2,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_CO2,NY,NX) &
    +(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*CO2_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*CO2_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_CO2,0,NY,NX)+trcg_surf_disevap_flx(idg_CO2,NY,NX)
  CH=Gas_Flx_atmDif2soil_col(idg_CH4,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_CH4,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_CH4,NY,NX) &
    +(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*CH4_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*CH4_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_CH4,0,NY,NX)+trcg_surf_disevap_flx(idg_CH4,NY,NX)
  CO                            = -IrrigSubsurf_col(NY,NX)*CO2_irrig_conc(NY,NX)
  CX                            = -IrrigSubsurf_col(NY,NX)*CH4_irrig_conc(NY,NX)
  SurfGasFlx_col(idg_CO2,NY,NX) = SurfGasFlx_col(idg_CO2,NY,NX)+CI
  SurfGasFlx_col(idg_CH4,NY,NX) = SurfGasFlx_col(idg_CH4,NY,NX)+CH
  SurfGas_CO2_lnd               = SurfGas_CO2_lnd+CI+CH
  TOMOU_lnds(ielmc)             = TOMOU_lnds(ielmc)+CO+CX
  !
  !     SURFACE BOUNDARY O2 FLUXES
  !
  OI=Gas_Flx_atmDif2soil_col(idg_O2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_O2,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_O2,NY,NX) &
    +(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*O2_rain_conc(NY,NX) &
    +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*O2_irrig_conc(NY,NX) &
    +Gas_Disol_Flx_vr(idg_O2,0,NY,NX)+trcg_surf_disevap_flx(idg_O2,NY,NX)
  SurfGas_O2_lnd               = SurfGas_O2_lnd+OI
  OO                           = trcg_RMicbTransf_vr(idg_O2,0,NY,NX)-IrrigSubsurf_col(NY,NX)*O2_irrig_conc(NY,NX)
  OXYGOU                       = OXYGOU+OO
  SurfGasFlx_col(idg_O2,NY,NX) = SurfGasFlx_col(idg_O2,NY,NX)+OI
  HI                           = Gas_Flx_atmDif2soil_col(idg_H2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_H2,3,NU(NY,NX),NY,NX) &
    +TRootGasLossDisturb_pft(idg_H2,NY,NX)+Gas_Disol_Flx_vr(idg_H2,0,NY,NX)+trcg_surf_disevap_flx(idg_H2,NY,NX)
  SurfGas_H2_lnd = SurfGas_H2_lnd+HI
  HO             = trcg_RMicbTransf_vr(idg_H2,0,NY,NX)
  H2GOU          = H2GOU+HO
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
  TOMOU_lnds(ielmn)=TOMOU_lnds(ielmn)+ZXB
  ZGI=(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*(N2_rain_conc(NY,NX)+N2O_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX)) &
      +Gas_Flx_atmDif2soil_col(idg_N2,NY,NX)+Gas_Flx_atmDif2soil_col(idg_N2O,NY,NX)+Gas_Flx_atmDif2soil_col(idg_NH3,NY,NX) &
      +Gas_Flx_atmDif2soil_col(idg_NH3B,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX) &
      +Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX)+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
      +TRootGasLossDisturb_pft(idg_N2O,NY,NX)+TRootGasLossDisturb_pft(idg_NH3,NY,NX) &
      +Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,0,NY,NX) &
      +trcg_surf_disevap_flx(idg_N2,NY,NX)+trcg_surf_disevap_flx(idg_N2O,NY,NX) &
      +trcg_surf_disevap_flx(idg_NH3,NY,NX)
  SurfGas_N2_lnd=SurfGas_N2_lnd+ZGI
  ZDRAIN(NY,NX)=ZDRAIN(NY,NX)+trcs_Transp2MicP_3D(ids_NH4,3,NK(NY,NX),NY,NX) &
      +trcs_Transp2MicP_3D(idg_NH3,3,NK(NY,NX),NY,NX)+trcs_Transp2MicP_3D(ids_NO3,3,NK(NY,NX),NY,NX) &
      +trcs_Transp2MicP_3D(ids_NO2,3,NK(NY,NX),NY,NX)+trcs_Transp2MicP_3D(ids_NH4B,3,NK(NY,NX),NY,NX) &
      +trcs_Transp2MicP_3D(idg_NH3B,3,NK(NY,NX),NY,NX)+trcs_Transp2MicP_3D(ids_NO3B,3,NK(NY,NX),NY,NX) &
      +trcs_Transp2MicP_3D(ids_NO2B,3,NK(NY,NX),NY,NX)

  ZNGGIN=Gas_Flx_atmDif2soil_col(idg_N2,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX)+Gas_Disol_Flx_vr(idg_N2,0,NY,NX)
  ZN2OIN=Gas_Flx_atmDif2soil_col(idg_N2O,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_N2O,3,NU(NY,NX),NY,NX)+Gas_Disol_Flx_vr(idg_N2O,0,NY,NX)
  ZNH3IN=Gas_Flx_atmDif2soil_col(idg_NH3,NY,NX)+Gas_Flx_atmDif2soil_col(idg_NH3B,NY,NX)+Gas_3DAdvDif_Flx_vr(idg_NH3,3,NU(NY,NX),NY,NX) &
    +Gas_Disol_Flx_vr(idg_NH3,0,NY,NX)

  SurfGasFlx_col(idg_N2,NY,NX)  = SurfGasFlx_col(idg_N2,NY,NX)+ZNGGIN
  SurfGasFlx_col(idg_N2O,NY,NX) = SurfGasFlx_col(idg_N2O,NY,NX)+ZN2OIN
  SurfGasFlx_col(idg_NH3,NY,NX) = SurfGasFlx_col(idg_NH3,NY,NX)+ZNH3IN
  SurfGasFlx_col(idg_N2,NY,NX)  = SurfGasFlx_col(idg_N2,NY,NX)+Micb_N2Fixation_vr(0,NY,NX)
  SurfGasFlx_col(idg_H2,NY,NX)  = SurfGasFlx_col(idg_H2,NY,NX)+HI
  !
  !     SURFACE BOUNDARY PO4 AND DOP FLUXES
  !
  PI=patomw*((Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX)) &
      *(H2PO4_rain_conc(NY,NX)+HPO4_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(H2PO4_irrig_conc(I,NY,NX)+HPO4_irrig_conc(I,NY,NX)))
  PXB=-patomw*IrrigSubsurf_col(NY,NX)*(H2PO4_irrig_conc(I,NY,NX)+HPO4_irrig_conc(I,NY,NX))
  TPIN=TPIN+PI
  TOMOU_lnds(ielmp)=TOMOU_lnds(ielmp)+PXB
  PDRAIN(NY,NX)=PDRAIN(NY,NX)+trcs_Transp2MicP_3D(ids_H2PO4,3,NK(NY,NX),NY,NX) &
      +trcs_Transp2MicP_3D(ids_H2PO4B,3,NK(NY,NX),NY,NX)+trcs_Transp2MicP_3D(ids_H1PO4,3,NK(NY,NX),NY,NX) &
      +trcs_Transp2MicP_3D(ids_H1PO4B,3,NK(NY,NX),NY,NX)
  !
  !     SURFACE BOUNDARY ION FLUXES
  !
  SIN=((Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX)) &
      *(2.0_r8*NH4_rain_conc(NY,NX)+NH3_rain_conc(NY,NX)+NO3_rain_conc(NY,NX)) &
      +(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX)) &
      *(2.0*NH4_irrig_conc(I,NY,NX)+NH3_irrig_conc(I,NY,NX)+NO3_irrig_conc(I,NY,NX)))
  SGN=(2.0_r8*(Rain2SoilSurf_col(NY,NX)+Rain2LitRSurf_col(NY,NX))*(N2_rain_conc(NY,NX)+N2O_rain_conc(NY,NX)) &
      +2.0_r8*(Irrig2SoilSurf(NY,NX)+Irrig2LitRSurf(NY,NX))*(N2_irrig_conc(NY,NX)+N2O_irrig_conc(NY,NX)) &
      +2.0_r8*(Gas_Flx_atmDif2soil_col(idg_N2,NY,NX)+Gas_Flx_atmDif2soil_col(idg_N2O,NY,NX))+Gas_Flx_atmDif2soil_col(idg_NH3,NY,NX) &
      +Gas_Flx_atmDif2soil_col(idg_NH3B,NY,NX)+2.0_r8*(Gas_3DAdvDif_Flx_vr(idg_N2,3,NU(NY,NX),NY,NX) &
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
  DO NE=1,NumPlantChemElms
    Litrfall_lnds(NE)=Litrfall_lnds(NE)+LitrFallStrutElms_col(NE,NY,NX)
    LiterfalOrgM_col(NE,NY,NX)=LiterfalOrgM_col(NE,NY,NX)+LitrFallStrutElms_col(NE,NY,NX)
  ENDDO
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
    SII=SII*IrrigSurface_col(NY,NX)
    TIONIN=TIONIN+SIR+SII    
    !
    !     SUBSURFACE BOUNDARY SALT FLUXES FROM SUBSURFACE IRRIGATION
    !
    TIONOU=TIONOU+SBU
  ENDIF
  !========================================================================
  !
  D9680: DO K=1,micpar%NumOfLitrCmplxs 
    do idom=idom_beg,idom_end   
      DOM_vr(idom,K,0,NY,NX)=DOM_vr(idom,K,0,NY,NX)+DOM_MicpTransp_3D(idom,K,3,0,NY,NX)
    enddo
  ENDDO D9680

  Eco_HR_CumYr_col(NY,NX)      = Eco_HR_CumYr_col(NY,NX)+trcg_RMicbTransf_vr(idg_CO2,0,NY,NX)+trcg_RMicbTransf_vr(idg_CH4,0,NY,NX)
  SurfGasFlx_col(idg_N2,NY,NX) = SurfGasFlx_col(idg_N2,NY,NX)+trcg_RMicbTransf_vr(idg_N2,0,NY,NX)

  RO2GasXchangePrev_vr(0,NY,NX) = Gas_Disol_Flx_vr(idg_O2,0,NY,NX)
  RCO2GasFlxPrev_vr(0,NY,NX)    = Gas_Disol_Flx_vr(idg_CO2,0,NY,NX)
  RCH4F(0,NY,NX)                = Gas_Disol_Flx_vr(idg_CH4,0,NY,NX)

  !RO2AquaXchangePrev_vr:=soil surface O2 dissolution + aqueous O2 flux micropore
  RO2AquaXchangePrev_vr(0,NY,NX)=trcg_surf_disevap_flx(idg_O2,NY,NX)+trcs_Transp2MicP_3D(idg_O2,3,0,NY,NX) &
    -(Rain2LitRSurf_col(NY,NX)*O2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*O2_irrig_conc(NY,NX))

  RCH4PhysexchPrev_vr(0,NY,NX)=trcg_surf_disevap_flx(idg_CH4,NY,NX)+trcs_Transp2MicP_3D(idg_CH4,3,0,NY,NX) &
    -(Rain2LitRSurf_col(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*CH4_irrig_conc(NY,NX))
    
  RO2AquaXchangePrev_vr(NU(NY,NX),NY,NX) = RO2AquaXchangePrev_vr(NU(NY,NX),NY,NX)+Gas_Flx_atmDif2soil_col(idg_O2,NY,NX)
  RCH4PhysexchPrev_vr(NU(NY,NX),NY,NX)   = RCH4PhysexchPrev_vr(NU(NY,NX),NY,NX)+Gas_Flx_atmDif2soil_col(idg_CH4,NY,NX)
  !
  !     SURFACE LITTER ION EXCHANGE AND PRECIPITATION
  !
  trcx_solml_vr(idx_NH4,0,NY,NX)=trcx_solml_vr(idx_NH4,0,NY,NX)+trcx_TRSoilChem_vr(idx_NH4,0,NY,NX)
  DO NTX=idx_AEC+1,idx_anion_soil_end
    trcx_solml_vr(NTX,0,NY,NX)=trcx_solml_vr(NTX,0,NY,NX)+trcx_TRSoilChem_vr(NTX,0,NY,NX)
  ENDDO

  DO NTP=idsp_psoi_beg,idsp_psoi_end
    trcp_saltpml_vr(NTP,0,NY,NX)=trcp_saltpml_vr(NTP,0,NY,NX)+trcp_RChem_soil(NTP,0,NY,NX)
  ENDDO
  !  !

  end subroutine HandleSurfaceBoundary
!------------------------------------------------------------------------------------------

  subroutine OverlandFlow(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX

  integer :: K,idom
  !     begin_execution
  !
  if(lverb)write(*,*)'OverlandFlow'
  IF(ABS(TXGridSurfRunoff_2DH(NY,NX)).GT.ZEROS(NY,NX))THEN
    !
    !     DOC, DON, DOP
    !
    D8570: DO K=1,micpar%NumOfLitrCmplxs    
      DO idom=idom_beg,idom_end
        DOM_vr(idom,K,0,NY,NX)=DOM_vr(idom,K,0,NY,NX)+TOMQRS(idom,K,NY,NX)
      ENDDO
    ENDDO D8570
!
    call OverlandFlowThruSnow(I,J,NY,NX)

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
  if(lverb)write(*,*)'SoilErosion'
  IF((iErosionMode.EQ.ieros_frzthaweros .OR. iErosionMode.EQ.ieros_frzthawsomeros) &
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
    FertN_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)       = FertN_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)+TNH4Eros_col(NY,NX)
    FertN_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)       = FertN_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)+TNH3Eros_col(NY,NX)
    FertN_soil_vr(ifert_urea,NU(NY,NX),NY,NX)      = FertN_soil_vr(ifert_urea,NU(NY,NX),NY,NX)+TNUreaEros_col(NY,NX)
    FertN_soil_vr(ifert_no3,NU(NY,NX),NY,NX)       = FertN_soil_vr(ifert_no3,NU(NY,NX),NY,NX)+TNO3Eros_col(NY,NX)
    FertN_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)  = FertN_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)+TNH4ErosBand_col(NY,NX)
    FertN_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)  = FertN_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)+TNH3ErosBand_col(NY,NX)
    FertN_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX) = FertN_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX)+TNUreaErosBand_col(NY,NX)
    FertN_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)  = FertN_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)+TNO3ErosBand_col(NY,NX)
!
    !   EXCHANGEABLE CATIONS AND ANIONS
!
    DO NTX=idx_beg,idx_end
      trcx_solml_vr(NTX,NU(NY,NX),NY,NX)=trcx_solml_vr(NTX,NU(NY,NX),NY,NX)+trcx_TER(NTX,NY,NX)
    ENDDO

!
!     PRECIPITATES
!
    DO NTP=idsp_beg,idsp_end
      trcp_saltpml_vr(NTP,NU(NY,NX),NY,NX)=trcp_saltpml_vr(NTP,NU(NY,NX),NY,NX)+trcp_TER(NTP,NY,NX)
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
  real(r8) :: OS
  real(r8) :: POS,POX,POP
  real(r8) :: SSS,ZG,Z4S,WS,Z4X
  real(r8) :: Z4F,ZOS,ZOF

  if(lverb)write(*,*)'LitterLayerSummary'
!     begin_execution
!     TOTAL C,N,P, SALTS IN SURFACE RESIDUE
!
  call sumSurfOMCK(NY,NX,RC0(:,NY,NX),RC0ff(NY,NX))

  call sumMicBiomLayL(0,NY,NX,tMicBiome_col(1:NumPlantChemElms,NY,NX))

  call sumLitrOMLayL(0,NY,NX,litrOM_vr(1:NumPlantChemElms,0,NY,NX))

  SoilOrgM_vr(1:NumPlantChemElms,0,NY,NX)=litrOM_vr(1:NumPlantChemElms,0,NY,NX)

  !update heat capacity for consistency
  VHeatCapacity_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
  DO NE=1,NumPlantChemElms
    tSoilOrgM_col(NE,NY,NX)=tSoilOrgM_col(NE,NY,NX)+SoilOrgM_vr(NE,0,NY,NX)
  enddo

  OMLitrC_vr(0,NY,NX)=litrOM_vr(ielmc,0,NY,NX)

  DO NE=1,NumPlantChemElms
    LitRMStoreLndscap(NE)=LitRMStoreLndscap(NE)+litrOM_vr(NE,0,NY,NX)
    tLitrOM_col(NE,NY,NX)=tLitrOM_col(NE,NY,NX)+litrOM_vr(NE,0,NY,NX)
  ENDDO
  !water mass 
  WS                   = CanH2OHeldVg_col(NY,NX)+CanWat_col(NY,NX)+VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)*DENSICE
  WatMassStore_lnd     = WatMassStore_lnd+WS
  WatMass_col(NY,NX)   = WatMass_col(NY,NX)+WS
  HeatStore_col(NY,NX) = HeatStore_col(NY,NX)+CanopyHeatStor_col(NY,NX)
  HeatStore_lnd        = HeatStore_lnd+CanopyHeatStor_col(NY,NX)
  CS                   = trc_solml_vr(idg_CO2,0,NY,NX)+trc_solml_vr(idg_CH4,0,NY,NX)
  TGasC_lnd           = TGasC_lnd+CS
  DIC_mass_col(NY,NX) = DIC_mass_col(NY,NX)+CS
  HS                  = trc_solml_vr(idg_H2,0,NY,NX)
  TSoilH2G_lnd        = TSoilH2G_lnd+HS
  OS                  = trc_solml_vr(idg_O2,0,NY,NX)
  TSoilO2G_lnd        = TSoilO2G_lnd+OS
  ZG                  = trc_solml_vr(idg_N2,0,NY,NX)+trc_solml_vr(idg_N2O,0,NY,NX)
  TGasN_lnd           = TGasN_lnd+ZG
  Z4S                 = trc_solml_vr(ids_NH4,0,NY,NX)+trc_solml_vr(idg_NH3,0,NY,NX)
  if(Z4S<0._r8)then
  write(*,*)'redistsurfnh4',NY,NX,trc_solml_vr(ids_NH4,0,NY,NX),trc_solml_vr(idg_NH3,0,NY,NX)
  stop
  endif
  Z4X=natomw*trcx_solml_vr(idx_NH4,0,NY,NX)
  Z4F=natomw*(FertN_soil_vr(ifert_nh4,0,NY,NX)+FertN_soil_vr(ifert_urea,0,NY,NX) &
    +FertN_soil_vr(ifert_nh3,0,NY,NX))
  TDisolNH4_lnd=TDisolNH4_lnd+Z4S+Z4X+Z4F
  tNH4_col(NY,NX)=tNH4_col(NY,NX)+Z4S+Z4X

  ZOS             = trc_solml_vr(ids_NO3,0,NY,NX)+trc_solml_vr(ids_NO2,0,NY,NX)
  ZOF             = natomw*FertN_soil_vr(ifert_no3,0,NY,NX)
  tNO3_lnd        = tNO3_lnd+ZOS+ZOF
  tNO3_col(NY,NX) = tNO3_col(NY,NX)+ZOS

  POS=trc_solml_vr(ids_H1PO4,0,NY,NX)+trc_solml_vr(ids_H2PO4,0,NY,NX)
  POX=patomw*(trcx_solml_vr(idx_HPO4,0,NY,NX)+trcx_solml_vr(idx_H2PO4,0,NY,NX))
  POP=patomw*(trcp_saltpml_vr(idsp_AlPO4,0,NY,NX)+trcp_saltpml_vr(idsp_FePO4,0,NY,NX) &
    +trcp_saltpml_vr(idsp_CaHPO4,0,NY,NX))+2._r8*patomw*trcp_saltpml_vr(idsp_CaH4P2O8,0,NY,NX) &
    +3._r8*patomw*trcp_saltpml_vr(idsp_HA,0,NY,NX)
  TDisolPi_lnd      = TDisolPi_lnd+POS+POX+POP
  tHxPO4_col(NY,NX) = tHxPO4_col(NY,NX)+POX
  tXPO4_col(NY,NX)  = tXPO4_col(NY,NX)+POP

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
    trcSalt_solml_vr(nsalts,0,NY,NX)=trcSalt_solml_vr(nsalts,0,NY,NX)+trcSalt3DFlo2Cell(nsalts,3,0,NY,NX)
  ENDDO

  PSS=patomw*(trcSalt_solml_vr(idsalt_H0PO4,0,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4,0,NY,NX)&
    +trcSalt_solml_vr(idsalt_FeHPO4,0,NY,NX) &
    +trcSalt_solml_vr(idsalt_FeH2PO4,0,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4,0,NY,NX)&
    +trcSalt_solml_vr(idsalt_CaHPO4,0,NY,NX) &
    +trcSalt_solml_vr(idsalt_CaH4P2O8,0,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4,0,NY,NX))
  TDisolPi_lnd=TDisolPi_lnd+PSS
  SSS=0._r8
  do nsalts=idsalt_beg,idsalt_end
    SSS=SSS+trcSalt_solml_vr(nsalts,0,NY,NX)*trcSaltIonNumber(nsalts)
  enddo

  TION=TION+SSS
  UION(NY,NX)=UION(NY,NX)+SSS
  end subroutine DiagSurfLitRLayerSalt
!------------------------------------------------------------------------------------------
  subroutine update_physVar_Profile(I,J,NY,NX,VOLISO,DVLiceMicP_vr)    
  !     WATER, ICE, HEAT, TEMPERATUR
  !
  use PerturbationMod, only : check_Soil_Warming,is_warming_layerL
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8), intent(inout) :: VOLISO  
  REAL(R8),INTENT(OUT) :: DVLiceMicP_vr(JZ)  !change in ice volume
  real(r8) :: TKS00,TKSX,TKSpre
  real(r8) :: ENGY
  real(r8) :: TVHeatCapacity
  real(r8) :: TVHeatCapacitySoilM,TVOLW,TVOLWH,TVOLI,TVOLIH,TENGY
  real(r8) :: VOLWXX,VOLIXX,VHeatCapacityX,WS,HS
  real(r8) :: DVLWatMicP_vr(JZ,JY,JX)   !change in water volume
  integer :: L,it
  logical :: do_warming

  do_warming=check_Soil_Warming(iYearCurrent,I)
  it=(I-1)*24+J

  if(lverb)write(*,*)'update_physVar_Profile'
  TVHeatCapacity      = 0.0_r8
  TVHeatCapacitySoilM = 0.0_r8
  TVOLW               = 0.0_r8
  TVOLWH              = 0.0_r8
  TVOLI               = 0.0_r8
  TVOLIH              = 0.0_r8
  TENGY               = 0.0_r8
  DVLiceMicP_vr       = 0._r8
  
  DO L=NU(NY,NX),NL(NY,NX)
    TKSX           = TKS_vr(L,NY,NX)
    VHeatCapacityX = VHeatCapacity_vr(L,NY,NX)
    VOLWXX         = VLWatMicP_vr(L,NY,NX)
    VOLIXX         = VLiceMicP_vr(L,NY,NX)
    !micropore
    VLWatMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)+TWatFlowCellMicP_vr(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
      +WatIceThawMicP_vr(L,NY,NX)+TPlantRootH2OUptake_vr(L,NY,NX)+FWatIrrigate2MicP_vr(L,NY,NX)
    VLWatMicPX_vr(L,NY,NX)=VLWatMicPX_vr(L,NY,NX)+TWatFlowCellMicPX_vr(L,NY,NX)+FWatExMacP2MicP(L,NY,NX) &
      +WatIceThawMicP_vr(L,NY,NX)+TPlantRootH2OUptake_vr(L,NY,NX)+FWatIrrigate2MicP_vr(L,NY,NX)

    !do a numerical correction
    VLWatMicPX_vr(L,NY,NX) = AMIN1(VLWatMicP_vr(L,NY,NX),VLWatMicPX_vr(L,NY,NX)+0.01_r8*(VLWatMicP_vr(L,NY,NX)-VLWatMicPX_vr(L,NY,NX)))
    VLiceMicP_vr(L,NY,NX)  = VLiceMicP_vr(L,NY,NX)-WatIceThawMicP_vr(L,NY,NX)/DENSICE

    !macropore
    VLWatMacP_vr(L,NY,NX)=VLWatMacP_vr(L,NY,NX)+TWaterFlowMacP_vr(L,NY,NX)-FWatExMacP2MicP(L,NY,NX)+WatIceThawMacP_vr(L,NY,NX)
    VLiceMacP_vr(L,NY,NX)=VLiceMacP_vr(L,NY,NX)-WatIceThawMacP_vr(L,NY,NX)/DENSICE

    !volume change
    DVLWatMicP_vr(L,NY,NX) = VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)
    DVLiceMicP_vr(L)       = VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX)

    !update water/ice-unfilled pores
    IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      VLsoiAirP_vr(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX) &
        +VLMacP_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX))
    ELSE
      VLsoiAirP_vr(L,NY,NX)=0.0_r8
!     VLMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)+VLiceMicP_vr(L,NY,NX)
!    2+DVLWatMicP_vr(L,NY,NX)+DVLiceMicP_vr(L,NY,NX)
!     VLSoilPoreMicP_vr(L,NY,NX)=VLMicP_vr(L,NY,NX)
!     VGeomLayer_vr(L,NY,NX)=VLMicP_vr(L,NY,NX)
    ENDIF
    ENGY=VHeatCapacityX*TKSX
    VHeatCapacity_vr(L,NY,NX)=VHeatCapacitySoilM_vr(L,NY,NX)+cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
      +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))

    TVHeatCapacity      = TVHeatCapacity+VHeatCapacity_vr(L,NY,NX)
    TVHeatCapacitySoilM = TVHeatCapacitySoilM+VHeatCapacitySoilM_vr(L,NY,NX)
    TVOLW               = TVOLW+VLWatMicP_vr(L,NY,NX)
    TVOLWH              = TVOLWH+VLWatMacP_vr(L,NY,NX)
    TVOLI               = TVOLI+VLiceMicP_vr(L,NY,NX)
    TVOLIH              = TVOLIH+VLiceMacP_vr(L,NY,NX)
    TENGY               = TENGY+ENGY
    !
    !     ARTIFICIAL SOIL WARMING
    !
    !     IF(NX.EQ.3.AND.NY.EQ.2.AND.L.GT.NU(NY,NX)
    !    3.AND.L.LE.17.AND.I.GE.152.AND.I.LE.304)THEN
    !     THeatFlow2Soil_vr(L,NY,NX)=THeatFlow2Soil_vr(L,NY,NX)
    !    2+(TKSZ(I,J,L)-TKS_vr(L,NY,NX))*VHeatCapacity_vr(L,NY,NX)
    !     WRITE(*,3379)'TKSZ',I,J,NX,NY,L,TKSZ(I,J,L)
    !    2,TKS_vr(L,NY,NX),VHeatCapacity_vr(L,NY,NX),THeatFlow2Soil_vr(L,NY,NX)
    !3379  FORMAT(A8,6I4,12E12.4)
    !     ENDIF
    !
    !     END ARTIFICIAL SOIL WARMING
    !
    TKSpre=TKS_vr(L,NY,NX)
    IF(VHeatCapacity_vr(L,NY,NX).GT.ZEROS(NY,NX) .and. VHeatCapacity_vr(L,NY,NX)/(VHeatCapacityX+VHeatCapacity_vr(L,NY,NX))>0.05_r8)THEN
      if(do_warming .and. is_warming_layerL(L,NY,NX))then     
        THeatFlow2Soil_vr(L,NY,NX) = THeatFlow2Soil_vr(L,NY,NX) + &
          (TKS_ref_vr(it,L,NY,NX)-TKS_vr(L,NY,NX))*VHeatCapacity_vr(L,NY,NX)
      endif   

      TKS00=TKS_vr(L,NY,NX)
      TKS_vr(L,NY,NX)=(ENGY+THeatFlow2Soil_vr(L,NY,NX)+THeatSoiThaw_vr(L,NY,NX) &
        +THeatRootUptake_vr(L,NY,NX)+HeatIrrigation(L,NY,NX))/VHeatCapacity_vr(L,NY,NX)
      tHeatUptk_col(NY,NX)=tHeatUptk_col(NY,NX)+THeatRootUptake_vr(L,NY,NX)      
!      if(curday>=285.and.L<=2)write(*,*)'rexL=',L,NY,NX,curhour,VHeatCapacityX,VHeatCapacity_vr(L,NY,NX),&
!        SoiBulkDensity_vr(L,NY,NX),NU(NY,NX)
      if(TKS_vr(L,NY,NX)>4.e2)then
        write(*,*)'weird temperature in redist',L,NY,NX,TKSpre,TKS_vr(L,NY,NX)
        write(*,*)'heatcap',VHeatCapacityX,VHeatCapacity_vr(L,NY,NX),ZEROS(NY,NX),SoiBulkDensity_vr(L,NY,NX)
        write(*,*)'itemized',ENGY/VHeatCapacity_vr(L,NY,NX),&
          THeatFlow2Soil_vr(L,NY,NX)/VHeatCapacity_vr(L,NY,NX),&
          THeatSoiThaw_vr(L,NY,NX)/VHeatCapacity_vr(L,NY,NX), &
          THeatRootUptake_vr(L,NY,NX)/VHeatCapacity_vr(L,NY,NX),&
          HeatIrrigation(L,NY,NX)/VHeatCapacity_vr(L,NY,NX)
        write(*,*)'wat',VLWatMicP_vr(L,NY,NX),VLWatMacP_vr(L,NY,NX), &
          VLiceMicP_vr(L,NY,NX),VLiceMacP_vr(L,NY,NX)  
        write(*,*)'heat',ENGY,THeatFlow2Soil_vr(L,NY,NX),VHeatCapacitySoilM_vr(L,NY,NX)  
        call endrun()  
      endif

    ELSE
      TKS_vr(L,NY,NX)=TKS_vr(NUM(NY,NX),NY,NX)
    ENDIF
    TCS(L,NY,NX)         = units%Kelvin2Celcius(TKS_vr(L,NY,NX))
    WS                   = VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)+(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE
    HS                   = VHeatCapacity_vr(L,NY,NX)*TKS_vr(L,NY,NX)
    WatMassStore_lnd     = WatMassStore_lnd+WS
    VOLISO               = VOLISO+VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)
    WatMass_col(NY,NX)   = WatMass_col(NY,NX)+WS
    HeatStore_col(NY,NX) = HeatStore_col(NY,NX)+HS
!    2-WiltPoint_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
    HeatStore_lnd=HeatStore_lnd+HS
  ENDDO
  end subroutine update_physVar_Profile
!------------------------------------------------------------------------------------------
  subroutine UpdateChemInSoilLays(I,J,NY,NX,LG,DORGC,TXCO2,DORGE)
  !

  implicit none
  integer,  intent(in) :: I,J
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
  real(r8) :: ZGB,Z2B,ZHB,dval,dval0,pval
  integer  :: idg,idom,ids,NE

  !     begin_execution
  !     UPDATE SOIL LAYER VARIABLES WITH TOTAL FLUXES

  if(lverb)write(*,*)'UpdateChemInSoilLays'
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
        DOM_vr(idom,K,L,NY,NX)=AZMAX1(DOM_vr(idom,K,L,NY,NX)+DOM_Transp2Micp_vr(idom,K,L,NY,NX) &
          +DOM_PoreTranspFlx(idom,K,L,NY,NX))
          
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
        if(DOM_vr(NE,K,L,NY,NX)<0._r8)print*,'933rootx',L,tRootMycoExud2Soil_vr(NE,K,L,NY,NX)
      ENDDO
    ENDDO D195

    !     MACROPORE SOLUTES FROM MACROPORE-MICROPORE EXCHANGE
    !
    if(SoilFracAsMacP_vr(L,NY,NX)>0._r8)then
      DO NTS=ids_beg,ids_end
        dval0=-trcs_Transp2MacP_vr(NTS,L,NY,NX)      
        if(dval0<0._r8)then
          trc_soHml_vr(NTS,L,NY,NX)=trc_soHml_vr(NTS,L,NY,NX)+trcs_Transp2MacP_vr(NTS,L,NY,NX) 
        else
          dval=dval0
          call fixEXflux(trc_soHml_vr(NTS,L,NY,NX),dval)      
          if(dval<dval0)then
            trcs_Transp2MacP_vr(NTS,L,NY,NX) =trcs_Transp2MacP_vr(NTS,L,NY,NX)*dval/dval0
          endif
        endif  
        dval0=trcs_PoreTranspFlx_vr(NTS,L,NY,NX)
        dval=dval0
        if(dval0<0._r8)then
          trc_soHml_vr(NTS,L,NY,NX)=trc_soHml_vr(NTS,L,NY,NX)-trcs_PoreTranspFlx_vr(NTS,L,NY,NX)
        else
          dval=dval0
          call fixEXflux(trc_soHml_vr(NTS,L,NY,NX),dval)      
          if(dval<dval0)then
            trcs_PoreTranspFlx_vr(NTS,L,NY,NX)=trcs_PoreTranspFlx_vr(NTS,L,NY,NX)*dval/dval0
          endif
        endif  
      ENDDO
    endif
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
      dval=trc_solml_vr(idg,L,NY,NX)
      trc_solml_vr(idg,L,NY,NX)=trc_solml_vr(idg,L,NY,NX) &
        +trcs_Transp2MicP_vr(idg,L,NY,NX)+Gas_Disol_Flx_vr(idg,L,NY,NX) &
        +trcs_Irrig_vr(idg,L,NY,NX)+trcs_PoreTranspFlx_vr(idg,L,NY,NX) &
        +trcg_ebu_flx_vr(idg,L,NY,NX)

       trc_solml_vr(idg,L,NY,NX)=fixnegmass(trc_solml_vr(idg,L,NY,NX))       
       if(trcg_RMicbTransf_vr(idg,L,NY,NX)<=0._r8)then
         !production
         trc_solml_vr(idg,L,NY,NX)=trc_solml_vr(idg,L,NY,NX)-trcg_RMicbTransf_vr(idg,L,NY,NX)
         call fixEXflux(trc_solml_vr(idg,L,NY,NX),trcs_plant_uptake_vr(idg,L,NY,NX))
       else
         dval0=trcs_plant_uptake_vr(idg,L,NY,NX) + trcg_RMicbTransf_vr(idg,L,NY,NX)
         dval=dval0
         call fixEXflux(trc_solml_vr(idg,L,NY,NX),dval)
         if(dval<dval0)then
           pval=dval/dval0
           trcs_plant_uptake_vr(idg,L,NY,NX) = trcs_plant_uptake_vr(idg,L,NY,NX)*pval
           trcg_RMicbTransf_vr(idg,L,NY,NX)  = trcg_RMicbTransf_vr(idg,L,NY,NX)*pval
         endif
       endif

      if(trc_solml_vr(idg,L,NY,NX)<0._r8)then
        print*,idg,'xxx',L,dval,trc_solml_vr(idg,L,NY,NX)
        print*,trcs_Transp2MicP_vr(idg,L,NY,NX)+Gas_Disol_Flx_vr(idg,L,NY,NX), &
          -trcg_RMicbTransf_vr(idg,L,NY,NX),-trcs_plant_uptake_vr(idg,L,NY,NX), &
          +trcs_Irrig_vr(idg,L,NY,NX)+trcs_PoreTranspFlx_vr(idg,L,NY,NX) &
          +trcg_ebu_flx_vr(idg,L,NY,NX)
      endif  
    enddo

    trc_solml_vr(idg_NH3,L,NY,NX)=trc_solml_vr(idg_NH3,L,NY,NX) &
      +TR_NH3_soil_vr(L,NY,NX) &
      +trcs_Transp2MicP_vr(idg_NH3,L,NY,NX)+Gas_Disol_Flx_vr(idg_NH3,L,NY,NX) &
      +trcs_Irrig_vr(idg_NH3,L,NY,NX) &
      +trcs_PoreTranspFlx_vr(idg_NH3,L,NY,NX)+trcg_ebu_flx_vr(idg_NH3,L,NY,NX)

    trc_solml_vr(idg_NH3,L,NY,NX)=fixnegmass(trc_solml_vr(idg_NH3,L,NY,NX))

    call fixEXflux(trc_solml_vr(idg_NH3,L,NY,NX),trcs_plant_uptake_vr(idg_NH3,L,NY,NX))

    trc_solml_vr(idg_NH3B,L,NY,NX)=trc_solml_vr(idg_NH3B,L,NY,NX) &
      +Gas_Disol_Flx_vr(idg_NH3B,L,NY,NX)+trcg_ebu_flx_vr(idg_NH3B,L,NY,NX) &
      +trcs_Transp2MicP_vr(idg_NH3B,L,NY,NX) &
      +trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX) &
      +trcs_Irrig_vr(idg_NH3B,L,NY,NX)+trcs_PoreTranspFlx_vr(idg_NH3B,L,NY,NX)
 
    trc_solml_vr(idg_NH3B,L,NY,NX)=fixnegmass(trc_solml_vr(idg_NH3B,L,NY,NX))
    call fixEXflux(trc_solml_vr(idg_NH3B,L,NY,NX),trcs_plant_uptake_vr(idg_NH3B,L,NY,NX))

   DO ids=ids_nut_beg,ids_nuts_end
      trc_solml_vr(ids,L,NY,NX)=trc_solml_vr(ids,L,NY,NX) &
        +trcs_Transp2MicP_vr(ids,L,NY,NX)  &
        +trcn_RChem_soil_vr(ids,L,NY,NX) &
        +trcs_Irrig_vr(ids,L,NY,NX)+trcs_PoreTranspFlx_vr(ids,L,NY,NX)     
      trc_solml_vr(ids,L,NY,NX)=fixnegmass(trc_solml_vr(ids,L,NY,NX))

      if(RNutMicbTransf_vr(ids,L,NY,NX)>=0._r8)then
        trc_solml_vr(ids,L,NY,NX)=trc_solml_vr(ids,L,NY,NX)+RNutMicbTransf_vr(ids,L,NY,NX)
        call fixEXflux(trc_solml_vr(ids,L,NY,NX),trcs_plant_uptake_vr(ids,L,NY,NX))
      else
        dval0=-RNutMicbTransf_vr(ids,L,NY,NX)+trcs_plant_uptake_vr(ids,L,NY,NX)
        dval=dval0
        call fixEXflux(trc_solml_vr(ids,L,NY,NX),dval)
        if(dval<dval0)then
          pval=dval/dval0
          RNutMicbTransf_vr(ids,L,NY,NX)    = RNutMicbTransf_vr(ids,L,NY,NX)*pval
          trcs_plant_uptake_vr(ids,L,NY,NX) = trcs_plant_uptake_vr(ids,L,NY,NX)*pval
        endif
      endif
    ENDDO
 
    do ids=ids_NH4B,ids_nutb_end
      trc_solml_vr(ids,L,NY,NX)=trc_solml_vr(ids,L,NY,NX) &
        +trcs_Transp2MicP_vr(ids,L,NY,NX) &
        +trcn_RChem_band_soil_vr(ids,L,NY,NX) &
        +trcs_Irrig_vr(ids,L,NY,NX)+trcs_PoreTranspFlx_vr(ids,L,NY,NX)
      trc_solml_vr(ids,L,NY,NX)=fixnegmass(trc_solml_vr(ids,L,NY,NX))
    
      if(RNutMicbTransf_vr(ids,L,NY,NX)>=0._r8)then
        trc_solml_vr(ids,L,NY,NX)=trc_solml_vr(ids,L,NY,NX)+RNutMicbTransf_vr(ids,L,NY,NX)
        call fixEXflux(trc_solml_vr(ids,L,NY,NX),trcs_plant_uptake_vr(ids,L,NY,NX))
      else
        dval0=-RNutMicbTransf_vr(ids,L,NY,NX)+trcs_plant_uptake_vr(ids,L,NY,NX)
        dval=dval0
        call fixEXflux(trc_solml_vr(ids,L,NY,NX),dval)
        if(dval<dval0)then
          pval=dval/dval0
          RNutMicbTransf_vr(ids,L,NY,NX)    = RNutMicbTransf_vr(ids,L,NY,NX)*pval
          trcs_plant_uptake_vr(ids,L,NY,NX) = trcs_plant_uptake_vr(ids,L,NY,NX)*pval
        endif
      endif  
    enddo

    Eco_HR_CumYr_col(NY,NX)=Eco_HR_CumYr_col(NY,NX)+trcg_RMicbTransf_vr(idg_CO2,L,NY,NX) &
      +trcg_RMicbTransf_vr(idg_CH4,L,NY,NX)
    SurfGasFlx_col(idg_N2,NY,NX)=SurfGasFlx_col(idg_N2,NY,NX)+trcg_RMicbTransf_vr(idg_N2,L,NY,NX)
    !
    !     EXCHANGEABLE CATIONS AND ANIONS FROM EXCHANGE REACTIONS
    !

    trcx_solml_vr(idx_NH4,L,NY,NX)=trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_TRSoilChem_vr(idx_NH4,L,NY,NX)
    trcx_solml_vr(idx_NH4B,L,NY,NX)=trcx_solml_vr(idx_NH4B,L,NY,NX)+trcx_TRSoilChem_vr(idx_NH4B,L,NY,NX)

    DO NTX=idx_AEC+1,idx_end
      trcx_solml_vr(NTX,L,NY,NX)=trcx_solml_vr(NTX,L,NY,NX)+trcx_TRSoilChem_vr(NTX,L,NY,NX)
    ENDDO

    !
    !     PRECIPITATES FROM PRECIPITATION-DISSOLUTION REACTIONS
    !
    DO NTP=idsp_p_beg,idsp_p_end
      trcp_saltpml_vr(NTP,L,NY,NX)=trcp_saltpml_vr(NTP,L,NY,NX)+trcp_RChem_soil(NTP,L,NY,NX)
    ENDDO
    !
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
    !  watch out the following code for changes
    HEATIN_lnd=HEATIN_lnd+THeatSoiThaw_vr(L,NY,NX)+THeatRootUptake_vr(L,NY,NX)
    CIB=trcg_air2root_flx_vr(idg_CO2,L,NY,NX)
    CHB=trcg_air2root_flx_vr(idg_CH4,L,NY,NX)
    OIB=trcg_air2root_flx_vr(idg_O2,L,NY,NX)
    HGB=trcg_air2root_flx_vr(idg_H2,L,NY,NX)
    ZGB=trcg_air2root_flx_vr(idg_N2,L,NY,NX)
    Z2B=trcg_air2root_flx_vr(idg_N2O,L,NY,NX)
    ZHB=trcg_air2root_flx_vr(idg_NH3,L,NY,NX)

    trcg_pltroot_flx_col(idg_beg:idg_NH3,NY,NX)=trcg_pltroot_flx_col(idg_beg:idg_NH3,NY,NX)+trcg_air2root_flx_vr(idg_beg:idg_NH3,L,NY,NX)

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

      trcg_ebu_flx_col(idg_beg:idg_NH3,NY,NX)=trcg_ebu_flx_vr(idg_beg:idg_NH3,L,NY,NX)
      trcg_ebu_flx_col(idg_NH3,NY,NX)=trcg_ebu_flx_col(idg_NH3,NY,NX)+trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)

    ELSE
      LL=MIN(L,LG)
      DO NTG=idg_beg,idg_NH3
        trc_gasml_vr(NTG,LL,NY,NX)=trc_gasml_vr(NTG,LL,NY,NX)-trcg_ebu_flx_vr(NTG,L,NY,NX)
      ENDDO
      trc_gasml_vr(idg_NH3,LL,NY,NX)=trc_gasml_vr(idg_NH3,LL,NY,NX)-trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)

      IF(LG.LT.L)THEN
        TGasC_lnd=TGasC_lnd-trcg_ebu_flx_vr(idg_CO2,L,NY,NX)-trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
        DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)-trcg_ebu_flx_vr(idg_CO2,L,NY,NX)-trcg_ebu_flx_vr(idg_CH4,L,NY,NX)
        TSoilO2G_lnd=TSoilO2G_lnd-trcg_ebu_flx_vr(idg_O2,L,NY,NX)
        TGasN_lnd=TGasN_lnd-trcg_ebu_flx_vr(idg_N2,L,NY,NX)-trcg_ebu_flx_vr(idg_N2O,L,NY,NX) &
          -trcg_ebu_flx_vr(idg_NH3,L,NY,NX)-trcg_ebu_flx_vr(idg_NH3B,L,NY,NX)
        TSoilH2G_lnd=TSoilH2G_lnd-trcg_ebu_flx_vr(idg_H2,L,NY,NX)
      ENDIF
    ENDIF
    SurfGas_CO2_lnd   = SurfGas_CO2_lnd+CIB+CHB
    COB               = tRootCO2Emis_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_CO2,L,NY,NX)-TR_CO2_aqu_soil_vr(L,NY,NX)
    TOMOU_lnds(ielmc) = TOMOU_lnds(ielmc)+COB

    RootResp_CumYr_col(NY,NX)  = RootResp_CumYr_col(NY,NX)+tRootCO2Emis_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_CO2,L,NY,NX)
    HydroSubsDICFlx_col(NY,NX) = HydroSubsDICFlx_col(NY,NX)-catomw*Txchem_CO2_vr(L,NY,NX)
    TXCO2(NY,NX)               = TXCO2(NY,NX)+catomw*Txchem_CO2_vr(L,NY,NX)
    SurfGas_O2_lnd             = SurfGas_O2_lnd+OIB
    OOB                        = trcg_RMicbTransf_vr(idg_O2,L,NY,NX)+tRO2MicrbUptk_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_O2,L,NY,NX)
    OXYGOU                     = OXYGOU+OOB
    SurfGas_H2_lnd             = SurfGas_H2_lnd+HGB
    HOB                        = trcg_RMicbTransf_vr(idg_H2,L,NY,NX)+trcs_plant_uptake_vr(idg_H2,L,NY,NX)
    H2GOU                      = H2GOU+HOB
    SurfGas_N2_lnd             = SurfGas_N2_lnd+ZGB+Z2B+ZHB
    
    SurfGasFlx_col(idg_CO2,NY,NX) = SurfGasFlx_col(idg_CO2,NY,NX)+CIB
    SurfGasFlx_col(idg_CH4,NY,NX) = SurfGasFlx_col(idg_CH4,NY,NX)+CHB
    SurfGasFlx_col(idg_O2,NY,NX)  = SurfGasFlx_col(idg_O2,NY,NX)+OIB
    SurfGasFlx_col(idg_N2,NY,NX)  = SurfGasFlx_col(idg_N2,NY,NX)+ZGB
    SurfGasFlx_col(idg_N2O,NY,NX) = SurfGasFlx_col(idg_N2O,NY,NX)+Z2B
    SurfGasFlx_col(idg_NH3,NY,NX) = SurfGasFlx_col(idg_NH3,NY,NX)+ZHB
    SurfGasFlx_col(idg_H2,NY,NX)  = SurfGasFlx_col(idg_H2,NY,NX)+HGB
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
    SSB=TRH2O(L,NY,NX)+Txchem_CO2_vr(L,NY,NX)+XZHYS(L,NY,NX)+TBION(L,NY,NX)
    TIONOU=TIONOU-SSB
    !     HydroIonFlx_CumYr_col(NY,NX)=HydroIonFlx_CumYr_col(NY,NX)-SSB
    !     WRITE(20,3339)'SSB',I,J,L,SSB,TRH2O(L,NY,NX)
    !    2,Txchem_CO2_vr(L,NY,NX),XZHYS(L,NY,NX),TBION(L,NY,NX)
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
      +trcg_root_vr(idg_CO2,L,NY,NX)+trc_gasml_vr(idg_CH4,L,NY,NX) &
      +trc_solml_vr(idg_CH4,L,NY,NX)+trcg_root_vr(idg_CH4,L,NY,NX) 

    TGasC_lnd=TGasC_lnd+CS
    DIC_mass_col(NY,NX)=DIC_mass_col(NY,NX)+CS
    HS=trc_gasml_vr(idg_H2,L,NY,NX)+trc_solml_vr(idg_H2,L,NY,NX) &
      +trcg_root_vr(idg_H2,L,NY,NX)


    OS=trc_gasml_vr(idg_O2,L,NY,NX)+trc_solml_vr(idg_O2,L,NY,NX) &
      +trcg_root_vr(idg_O2,L,NY,NX)


!  should I add N2 to ZG below? Dec, 18, 2022.
    ZG=trc_gasml_vr(idg_N2,L,NY,NX)+trc_solml_vr(idg_N2,L,NY,NX) &
      +trcg_root_vr(idg_N2O,L,NY,NX)+trc_gasml_vr(idg_N2O,L,NY,NX) &
      +trc_solml_vr(idg_N2O,L,NY,NX)+trcg_root_vr(idg_NH3,L,NY,NX) &
      +trc_gasml_vr(idg_NH3,L,NY,NX)  

    Z4S=trc_solml_vr(ids_NH4,L,NY,NX)+trc_solml_vr(ids_NH4B,L,NY,NX) &
      +trc_solml_vr(idg_NH3,L,NY,NX)+trc_solml_vr(idg_NH3B,L,NY,NX) 

    ZOS=trc_solml_vr(ids_NO3,L,NY,NX)+trc_solml_vr(ids_NO3B,L,NY,NX)&
      +trc_solml_vr(ids_NO2,L,NY,NX)+trc_solml_vr(ids_NO2B,L,NY,NX) 

    POS=trc_solml_vr(ids_H2PO4,L,NY,NX)+trc_solml_vr(ids_H2PO4B,L,NY,NX) &
      +trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_solml_vr(ids_H1PO4B,L,NY,NX)

    if(SoilFracAsMacP_vr(L,NY,NX)>0._r8)then
    CS=CS+trc_soHml_vr(idg_CO2,L,NY,NX)+trc_soHml_vr(idg_CH4,L,NY,NX)

    HS=HS+trc_soHml_vr(idg_H2,L,NY,NX)

    OS=OS+trc_soHml_vr(idg_O2,L,NY,NX)

    ZG=ZG+trc_soHml_vr(idg_N2,L,NY,NX)+trc_soHml_vr(idg_N2O,L,NY,NX)

    Z4S=Z4S+trc_soHml_vr(ids_NH4,L,NY,NX)+trc_soHml_vr(ids_NH4B,L,NY,NX)&
      +trc_soHml_vr(idg_NH3,L,NY,NX)+trc_soHml_vr(idg_NH3B,L,NY,NX)

    ZOS=ZOS+trc_soHml_vr(ids_NO3,L,NY,NX)+trc_soHml_vr(ids_NO3B,L,NY,NX) &
      +trc_soHml_vr(ids_NO2,L,NY,NX)+trc_soHml_vr(ids_NO2B,L,NY,NX)

    POS=POS+trc_soHml_vr(ids_H2PO4B,L,NY,NX)+trc_soHml_vr(ids_H2PO4,L,NY,NX) &
      +trc_soHml_vr(ids_H1PO4,L,NY,NX)+trc_soHml_vr(ids_H1PO4B,L,NY,NX)
    endif
    Z4X=natomw*(trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX))
    Z4F=natomw*(FertN_soil_vr(ifert_nh4,L,NY,NX)+FertN_soil_vr(ifert_urea,L,NY,NX) &
      +FertN_soil_vr(ifert_nh3,L,NY,NX)+FertN_Band_vr(ifert_nh4_band,L,NY,NX) &
      +FertN_Band_vr(ifert_urea_band,L,NY,NX)+FertN_Band_vr(ifert_nh3_band,L,NY,NX))

    TSoilH2G_lnd=TSoilH2G_lnd+HS
    TSoilO2G_lnd=TSoilO2G_lnd+OS
    TGasN_lnd=TGasN_lnd+ZG
    TDisolNH4_lnd=TDisolNH4_lnd+Z4S+Z4X+Z4F
    tNH4_col(NY,NX)=tNH4_col(NY,NX)+Z4S+Z4X


    ZOF=natomw*(FertN_soil_vr(ifert_no3,L,NY,NX)+FertN_soil_vr(ifert_no3,L,NY,NX))
    tNO3_lnd=tNO3_lnd+ZOS+ZOF
    tNO3_col(NY,NX)=tNO3_col(NY,NX)+ZOS

    !exchangeable P
    POX=patomw*(trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4,L,NY,NX) &
      +trcx_solml_vr(idx_HPO4B,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX))

    !precipitated P
    POP=patomw*(trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)+trcp_saltpml_vr(idsp_FePO4,L,NY,NX)&
      +trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)+trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)&
      +trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)+trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)) &
      +2._r8*patomw*(trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)+trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX)) &
      +3._r8*patomw*(trcp_saltpml_vr(idsp_HA,L,NY,NX)+trcp_saltpml_vr(idsp_HAB,L,NY,NX))

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
    trcSalt_solml_vr(nsalts,L,NY,NX)=trcSalt_solml_vr(nsalts,L,NY,NX)+trcSalt_TR(nsalts,L,NY,NX) &
      +trcSalt_Flo2MicP_vr(nsalts,L,NY,NX)+trcSalt_RFLU(nsalts,L,NY,NX)+trcSalt_XFXS(nsalts,L,NY,NX)

    trcSalt_soHml_vr(nsalts,L,NY,NX)=trcSalt_soHml_vr(nsalts,L,NY,NX)+trcSalt_Flo2MacP_vr(nsalts,L,NY,NX) &
      -trcSalt_XFXS(nsalts,L,NY,NX)
  ENDDO
  trcSalt_solml_vr(idsalt_AlOH2,L,NY,NX)=trcSalt_solml_vr(idsalt_AlOH2,L,NY,NX)-TR_AlO2H2_sorbed_soil(L,NY,NX)
  trcSalt_solml_vr(idsalt_FeOH2,L,NY,NX)=trcSalt_solml_vr(idsalt_FeOH2,L,NY,NX)-TR_FeO2H2_sorbed_soil_vr(L,NY,NX)

  trcx_solml_vr(idx_Hp,L,NY,NX)=trcx_solml_vr(idx_Hp,L,NY,NX)+TR_H_p_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_Al,L,NY,NX)=trcx_solml_vr(idx_Al,L,NY,NX)+TR_Al_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_Fe,L,NY,NX)=trcx_solml_vr(idx_Fe,L,NY,NX)+TR_Fe_sorbed_soil_vr(L,NY,NX)
  trcx_solml_vr(idx_Ca,L,NY,NX)=trcx_solml_vr(idx_Ca,L,NY,NX)+TR_Ca_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_Mg,L,NY,NX)=trcx_solml_vr(idx_Mg,L,NY,NX)+TR_Mg_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_Na,L,NY,NX)=trcx_solml_vr(idx_Na,L,NY,NX)+TR_Na_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_K,L,NY,NX)=trcx_solml_vr(idx_K,L,NY,NX)+TR_K_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_COOH,L,NY,NX)=trcx_solml_vr(idx_COOH,L,NY,NX)+TR_HCO3_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_AlOH2,L,NY,NX)=trcx_solml_vr(idx_AlOH2,L,NY,NX)+TR_AlO2H2_sorbed_soil(L,NY,NX)
  trcx_solml_vr(idx_FeOH2,L,NY,NX)=trcx_solml_vr(idx_FeOH2,L,NY,NX)+TR_FeO2H2_sorbed_soil_vr(L,NY,NX)

! all non-P precipitates
  DO NTP=idsp_beg,idsp_p_beg-1
    trcp_saltpml_vr(NTP,L,NY,NX)=trcp_saltpml_vr(NTP,L,NY,NX)+trcp_RChem_soil(NTP,L,NY,NX)
  ENDDO

  PSS=0._r8
  DO nsalts=idsalt_psoil_beg,idsalt_pband_end
    PSS=PSS+trcSalt_solml_vr(nsalts,L,NY,NX)+trcSalt_soHml_vr(nsalts,L,NY,NX)
  ENDDO
  PSS=PSS*patomw  

  TDisolPi_lnd=TDisolPi_lnd+PSS

  SSS=0._r8
  SSH=0._r8
  DO nsalts=idsalt_beg,idsaltb_end
    SSS=SSS+trcSalt_solml_vr(nsalts,L,NY,NX)*trcSaltIonNumber(nsalts)
    SSH=SSH+trcSalt_soHml_vr(idsalt_Al,L,NY,NX)*trcSaltIonNumber(nsalts)
  ENDDO

 
!     TOTAL FERILIZER,EXCHANGEABLE CATIONS AND ANIONS, PRECIPITATES
!
  SSF=FertN_soil_vr(ifert_nh3,L,NY,NX)+FertN_soil_vr(ifert_urea,L,NY,NX) &
    +FertN_soil_vr(ifert_no3,L,NY,NX) &
    +FertN_Band_vr(ifert_nh3_band,L,NY,NX)+FertN_Band_vr(ifert_urea_band,L,NY,NX) &
    +FertN_Band_vr(ifert_no3_band,L,NY,NX) &
    +2.0_r8*(FertN_soil_vr(ifert_nh4,L,NY,NX)+FertN_Band_vr(ifert_nh4_band,L,NY,NX))

  SSX=trcx_solml_vr(idx_Hp,L,NY,NX)+trcx_solml_vr(idx_Al,L,NY,NX) &
    +trcx_solml_vr(idx_Fe,L,NY,NX)+trcx_solml_vr(idx_Ca,L,NY,NX)+trcx_solml_vr(idx_Mg,L,NY,NX) &
    +trcx_solml_vr(idx_Na,L,NY,NX)+trcx_solml_vr(idx_K,L,NY,NX)+trcx_solml_vr(idx_COOH,L,NY,NX) &
    +trcx_solml_vr(idx_OHe,L,NY,NX)+trcx_solml_vr(idx_OHeB,L,NY,NX) &
    +2.0_r8*(trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX) &
    +trcx_solml_vr(idx_OH,L,NY,NX)+trcx_solml_vr(idx_OHB,L,NY,NX)) &
    +3.0_r8*(trcx_solml_vr(idx_AlOH2,L,NY,NX)+trcx_solml_vr(idx_FeOH2,L,NY,NX) &
    +trcx_solml_vr(idx_OHp,L,NY,NX)+trcx_solml_vr(idx_OHpB,L,NY,NX) &
    +trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_HPO4B,L,NY,NX)) &
    +4.0_r8*(trcx_solml_vr(idx_H2PO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX))

  SSP=2.0_r8*(trcp_saltpml_vr(idsp_CaCO3,L,NY,NX)+trcp_saltpml_vr(idsp_CaSO4,L,NY,NX) &
    +trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)+trcp_saltpml_vr(idsp_FePO4,L,NY,NX) &
    +trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)+trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)) &
    +3.0_r8*(trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)+trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)) &
    +4.0_r8*(trcp_saltpml_vr(idsp_AlOH3,L,NY,NX)+trcp_saltpml_vr(idsp_FeOH3,L,NY,NX)) &
    +7.0_r8*(trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)+trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX)) &
    +9.0_r8*(trcp_saltpml_vr(idsp_HA,L,NY,NX)+trcp_saltpml_vr(idsp_HAB,L,NY,NX))

  SST=SSS+SSH+SSF+SSX+SSP
  TION=TION+SST
  UION(NY,NX)=UION(NY,NX)+SST

!
!     SOIL ELECTRICAL CONDUCTIVITY
!
  IF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    ECHY=0.337_r8*AZMAX1(trcSalt_solml_vr(idsalt_Hp,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECOH=0.192_r8*AZMAX1(trcSalt_solml_vr(idsalt_OH,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECAL=0.056_r8*AZMAX1(trcSalt_solml_vr(idsalt_Al,L,NY,NX)*3.0_r8/VLWatMicP_vr(L,NY,NX))
    ECFE=0.051_r8*AZMAX1(trcSalt_solml_vr(idsalt_Fe,L,NY,NX)*3.0_r8/VLWatMicP_vr(L,NY,NX))
    ECCA=0.060_r8*AZMAX1(trcSalt_solml_vr(idsalt_Ca,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECMG=0.053_r8*AZMAX1(trcSalt_solml_vr(idsalt_Mg,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECNA=0.050_r8*AZMAX1(trcSalt_solml_vr(idsalt_Na,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECKA=0.070_r8*AZMAX1(trcSalt_solml_vr(idsalt_K,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECCO=0.072_r8*AZMAX1(trcSalt_solml_vr(idsalt_CO3,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECHC=0.044_r8*AZMAX1(trcSalt_solml_vr(idsalt_HCO3,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECSO=0.080_r8*AZMAX1(trcSalt_solml_vr(idsalt_SO4,L,NY,NX)*2.0_r8/VLWatMicP_vr(L,NY,NX))
    ECCL=0.076_r8*AZMAX1(trcSalt_solml_vr(idsalt_Cl,L,NY,NX)/VLWatMicP_vr(L,NY,NX))
    ECNO=0.071_r8*AZMAX1(trc_solml_vr(ids_NO3,L,NY,NX)/(VLWatMicP_vr(L,NY,NX)*natomw))
    ECND_vr(L,NY,NX)=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA &
      +ECCO+ECHC+ECSO+ECCL+ECNO

  ELSE
    ECND_vr(L,NY,NX)=0.0_r8
  ENDIF
  end subroutine UpdateSaltIonInSoilLayers

!------------------------------------------------------------------------------------------
  subroutine SumOMStates(L,NY,NX,DORGCL,DORGEC)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in):: DORGEC
  real(r8), intent(out):: DORGCL

  real(r8) :: ORGM(NumPlantChemElms)
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
  !total organic matter, litter + POM/humus, exclude living microbial biomass  
  call sumORGMLayL(L,NY,NX,SoilOrgM_vr(1:NumPlantChemElms,L,NY,NX),info='redist')
  
  DO NE=1,NumPlantChemElms
    if(SoilOrgM_vr(NE,L,NY,NX)<0._r8)write(*,*)'redist',NE,L,SoilOrgM_vr(1:NE,L,NY,NX)
    tSoilOrgM_col(NE,NY,NX)=tSoilOrgM_col(NE,NY,NX)+SoilOrgM_vr(NE,L,NY,NX)
  ENDDO

  !sum up litter complexes
  call sumLitrOMLayL(L,NY,NX,litrOM_vr(1:NumPlantChemElms,L,NY,NX))
  OMLitrC_vr(L,NY,NX)=litrOM_vr(ielmc,L,NY,NX)

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
    LitRMStoreLndscap(NE)=LitRMStoreLndscap(NE)+litrOM_vr(NE,L,NY,NX)
    tLitrOM_col(NE,NY,NX)=tLitrOM_col(NE,NY,NX)+litrOM_vr(NE,L,NY,NX)
    
    POMHumStoreLndscap(NE)=POMHumStoreLndscap(NE)+HumOM(NE)
    tHumOM_col(NE,NY,NX)=tHumOM_col(NE,NY,NX)+HumOM(NE)
  ENDDO
  TSEDSO=TSEDSO+(litrOM_vr(ielmc,L,NY,NX)+HumOM(ielmc))*ppmc
  end subroutine SumOMStates

!------------------------------------------------------------------------------------------

  subroutine AddFlux2SurfaceResidue(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX

  integer :: M,N,K,NGL,NE
  real(r8) :: HRAINR,RAINR
  real(r8) :: VLWatMicP1X,ENGYR,VLHeatCapacityX,TK1X
!     begin_execution
!     ADD ABOVE-GROUND LitrFall FROM EXTRACT.F TO SURFACE RESIDUE
!
!     OSC,OSN,OSP=SOC,SON,SOP
!     CSNT,ZSNT,PSNT=total C,N,P LitrFall
!     ORGC=total SOC
!     RAINR,HRAINR=water,heat in LitrFall
!     FLWR,HFLWR=water,heat flux into litter
!     HEATIN_lnd=cumulative net surface heat transfer
!

  DO   K=1,micpar%NumOfPlantLitrCmplxs
    DO  M=1,jsken
      SolidOMAct_vr(M,K,0,NY,NX)=SolidOMAct_vr(M,K,0,NY,NX)+LitrfalStrutElms_vr(ielmc,M,K,0,NY,NX)*micpar%OMCI(1,K)
      DO NE=1,NumPlantChemElms
        SolidOM_vr(NE,M,K,0,NY,NX)=SolidOM_vr(NE,M,K,0,NY,NX)+LitrfalStrutElms_vr(NE,M,K,0,NY,NX)
      ENDDO

      SoilOrgM_vr(ielmc,0,NY,NX) = SoilOrgM_vr(ielmc,0,NY,NX)+LitrfalStrutElms_vr(ielmc,M,K,0,NY,NX)
      RAINR                      = LitrfalStrutElms_vr(ielmc,M,K,0,NY,NX)*ThetaCX(K)
      HRAINR                     = RAINR*cpw*TairK_col(NY,NX)+LitrfalStrutElms_vr(ielmc,M,K,0,NY,NX)*cpo*TairK_col(NY,NX)
      WatFLo2Litr(NY,NX)         = WatFLo2Litr(NY,NX)+RAINR
      HeatFLo2LitrByWat(NY,NX)   = HeatFLo2LitrByWat(NY,NX)+HRAINR
      CRAIN                      = CRAIN+RAINR
      HEATIN_lnd                 = HEATIN_lnd+HRAINR
    enddo
  ENDDO

  call SumSurfMicBGCFluxes(I,J,NY,NX)

  end subroutine AddFlux2SurfaceResidue

!------------------------------------------------------------------------------------------

  subroutine SumSurfMicBGCFluxes(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
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
  !     RDOMEcoDmndK_vr,RAcetateEcoDmndK_vr=total DOC,DOA demand from DOC,DOA oxidation
  !     RDOCUptkHeter_vr,RAcetateUptkHeter_vr=DOC,DOA demand from DOC,DOA oxidation
  !     RNO2EcoUptkSoil_vr=total demand for NO2 reduction
  !     RNO2DmndSoilChemo_vr=demand for NO2 reduction
!
  DO K=1,jcplx
    IF(micpar%is_litter(K))THEN
      DO  N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          REcoO2DmndResp_vr(0,NY,NX)            = REcoO2DmndResp_vr(0,NY,NX)+RO2DmndHetert(NGL,K,0,NY,NX)
          REcoNO3DmndSoil_vr(0,NY,NX)           = REcoNO3DmndSoil_vr(0,NY,NX)+RNO3ReduxDmndSoilHeter_vr(NGL,K,0,NY,NX)
          RNO2EcoUptkSoil_vr(0,NY,NX)           = RNO2EcoUptkSoil_vr(0,NY,NX)+RNO2DmndReduxSoilHeter_vr(NGL,K,0,NY,NX)
          RN2OEcoUptkSoil_vr(0,NY,NX)           = RN2OEcoUptkSoil_vr(0,NY,NX)+RN2ODmndReduxHeter_vr(NGL,K,0,NY,NX)
          REcoNH4DmndSoil_vr(0,NY,NX)           = REcoNH4DmndSoil_vr(0,NY,NX)+RNH4DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoNO3DmndSoil_vr(0,NY,NX)           = REcoNO3DmndSoil_vr(0,NY,NX)+RNO3DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoH2PO4DmndSoil_vr(0,NY,NX)         = REcoH2PO4DmndSoil_vr(0,NY,NX)+RH2PO4DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoH1PO4DmndSoil_vr(0,NY,NX)         = REcoH1PO4DmndSoil_vr(0,NY,NX)+RH1PO4DmndSoilHeter_vr(NGL,K,0,NY,NX)
          REcoNH4DmndSoil_vr(NU(NY,NX),NY,NX)   = REcoNH4DmndSoil_vr(NU(NY,NX),NY,NX)+RNH4DmndLitrHeter_col(NGL,K,NY,NX)
          REcoNO3DmndSoil_vr(NU(NY,NX),NY,NX)   = REcoNO3DmndSoil_vr(NU(NY,NX),NY,NX)+RNO3DmndLitrHeter_col(NGL,K,NY,NX)
          REcoH2PO4DmndSoil_vr(NU(NY,NX),NY,NX) = REcoH2PO4DmndSoil_vr(NU(NY,NX),NY,NX)+RH2PO4DmndLitrHeter_col(NGL,K,NY,NX)
          REcoH1PO4DmndSoil_vr(NU(NY,NX),NY,NX) = REcoH1PO4DmndSoil_vr(NU(NY,NX),NY,NX)+RH1PO4DmndLitrHeter_col(NGL,K,NY,NX)
          RDOMEcoDmndK_vr(K,0,NY,NX)                      = RDOMEcoDmndK_vr(K,0,NY,NX)+RDOCUptkHeter_vr(NGL,K,0,NY,NX)
          RAcetateEcoDmndK_vr(K,0,NY,NX)                      = RAcetateEcoDmndK_vr(K,0,NY,NX)+RAcetateUptkHeter_vr(NGL,K,0,NY,NX)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DO  N=1,NumMicbFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      REcoO2DmndResp_vr(0,NY,NX)=REcoO2DmndResp_vr(0,NY,NX)+RO2DmndAutort_vr(NGL,0,NY,NX)
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
    RDOMEcoDmndK_vr(K,L,NY,NX)=RDOMEcoDmndK_vr(K,L,NY,NX)+SUM(RDOCUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
    RAcetateEcoDmndK_vr(K,L,NY,NX)=RAcetateEcoDmndK_vr(K,L,NY,NX)+SUM(RAcetateUptkHeter_vr(1:NumHetetrMicCmplx,1:jcplx,L,NY,NX))
  ENDDO
  REcoO2DmndResp_vr(L,NY,NX)=REcoO2DmndResp_vr(L,NY,NX)+SUM(RO2DmndAutort_vr(1:NumHetetrMicCmplx,L,NY,NX))
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
