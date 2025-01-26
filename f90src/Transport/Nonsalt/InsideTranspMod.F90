module InsideTranspMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: destroy
  use minimathmod,   only: safe_adb, AZMAX1, AZMIN1
  use EcoSIMCtrlMod, only: PrintInfo
  use GridConsts
  use EcoSIMSolverPar
  use TranspNoSaltDataMod
  use AqueChemDatatype
  use GridDataType
  use SoilBGCDataType
  use SurfSoilDataType
  use SurfLitterDataType
  use SnowDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use ClimForcDataType
  USE EcoSimConst
  USE SoilHeatDataType
  use SurfaceFluxMod
  use TracerIDMod
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__


  public :: ModelTracerHydroFluxMM
  public :: InitInsTp
  public :: DestructInsTp
  contains
!------------------------------------------------------------------------------------------
  subroutine InitInsTp()
  implicit none

  end subroutine InitInsTp
!------------------------------------------------------------------------------------------
  subroutine DestructInsTp
  implicit none
  

  end subroutine DestructInsTp
!------------------------------------------------------------------------------------------

  subroutine ModelTracerHydroFluxMM(I,J,M,MX,NHW, NHE, NVN, NVS,WaterFlow2Soil,RGasAtmDisol2LitrM,RGasAtmDisol2SoilM)
  !
  !Description
  !Hydrological tracer flux in iteration MM
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,MX
  integer, intent(in) :: NHW, NHE, NVN, NVS
  real(r8), intent(inout) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8), intent(out) :: RGasAtmDisol2LitrM(idg_beg:idg_NH3,JY,JX)  !atmospheric gas dissolves to litter solutes
  real(r8), intent(out) :: RGasAtmDisol2SoilM(idg_beg:idg_end,JY,JX)  !atmospheric gas dissolves to soil solutes
  character(len=*), parameter :: subname='ModelTracerHydroFluxMM' 
  integer :: NY,NX
  real(r8) :: FLWRM1
  real(r8) :: trcg_FloSno2LitR(idg_beg:idg_NH3)
  real(r8) :: trcn_FloSno2LitR(ids_nut_beg:ids_nuts_end)
  real(r8) :: RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)
  real(r8) :: RGas_Dif_Atm2Litr_FlxMM(idg_beg:idg_NH3)

  call PrintInfo('beg '//subname)

  RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)=0._r8

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      call ResetFluxAccumulatorsMM(I,J,M,NY,NX,MX)

      IF(M.NE.MX)THEN
        !this is a new iteration for moisture-temp update
 
        call SnowSoluteDischargeM(M,NY,NX,trcg_FloSno2LitR,trcn_FloSno2LitR)
!
!     SOLUTE FLUXES AT SOIL SURFACE FROM SURFACE WATER
!     CONTEnsolutes, WATER FLUXES 'WaterFlow2Soil' AND ATMOSPHERE BOUNDARY
!     LAYER RESISTANCES 'PARGM' FROM 'WATSUB'
!
        call SoluteFluxSurfaceM(I,J,M,NY,NX,NHE,NHW,NVS,NVN,&
          WaterFlow2Soil,trcg_FloSno2LitR,trcn_FloSno2LitR,RGas_Dif_Atm2Soil_FlxMM,&
          RGas_Dif_Atm2Litr_FlxMM,RGasAtmDisol2SoilM(:,NY,NX),RGasAtmDisol2LitrM(:,NY,NX))
!
      ENDIF
!
!     VOLATILIZATION-DISSOLUTION OF GASES IN RESIDUE AND SOIL SURFACE
!     LAYERS FROM GASEOUS CONCENTRATIONS VS. THEIR AQUEOUS
!     EQUIVALEnsolutes DEPENDING ON SOLUBILITY FROM 'HOUR1'
!     AND TRANSFER COEFFICIENT 'DiffusivitySolutEff' FROM 'WATSUB'
!
      call LitterGasVolatilDissolMM(I,J,M,NY,NX,RGas_Dif_Atm2Litr_FlxMM)
!
!     SURFACE GAS EXCHANGE FROM GAS DIFFUSIVITY THROUGH
!     SOIL SURFACE LAYER AND THROUGH ATMOSPHERE BOUNDARY
!     LAYER
!
      call SurfSoilFluxGasDifAdvMM(M,NY,NX,WaterFlow2Soil,RGas_Dif_Atm2Soil_FlxMM)
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
!
      call TracerExchXGridsMM(M,MX,NY,NX,NHE,NVS,WaterFlow2Soil,RGasAtmDisol2SoilM(:,NY,NX))

    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine ModelTracerHydroFluxMM
!------------------------------------------------------------------------------------------

  subroutine ResetFluxAccumulatorsMM(I,J,M,NY,NX,MX)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M, NY, NX,MX

  integer :: K,L,idg,nsol,idom
  real(r8) :: PARGM

  IF(M.NE.MX)THEN
    !
    !     GASEOUS BOUNDARY LAYER CONDUCTANCES
    !
    !     PARG=boundary layer conductance above soil surface from watsub.f
    !
    PARGM                   = CondGasXSnowM_col(M,NY,NX)*dt_GasCyc
    PARG_cef(idg_CO2,NY,NX) = PARGM*0.74_r8
    PARG_cef(idg_CH4,NY,NX) = PARGM*1.04_r8
    PARG_cef(idg_O2,NY,NX)  = PARGM*0.83_r8
    PARG_cef(idg_N2,NY,NX)  = PARGM*0.86_r8
    PARG_cef(idg_N2O,NY,NX) = PARGM*0.74_r8
    PARG_cef(idg_NH3,NY,NX) = PARGM*1.02_r8
    PARG_cef(idg_H2,NY,NX)  = PARGM*2.08_r8
!
!     RESET RUNOFF SOLUTE FLUX ACCUMULATORS
!
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in micropores
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in micropores
!     ZN3G=gaseous NH3
!

    DO  K=1,jcplx
      DOM_SurfRunoff_flxM(idom_beg:idom_end,K,NY,NX)=0.0_r8
    ENDDO

    trcg_SurfRunoff_flxM(idg_beg:idg_NH3,NY,NX)             = 0.0_r8
    trcn_SurfRunoff_flxM(ids_nut_beg:ids_nuts_end,NY,NX) = 0.0_r8
    trcg_SnowDrift_flxM(idg_beg:idg_NH3,NY,NX)                  = 0.0_r8
    trcn_SnowDrift_flxM(ids_nut_beg:ids_nuts_end,NY,NX)         = 0.0_r8

!   because NH3 is gas-aqua dual phase
    trc_solml2_vr(idg_NH3,0,NY,NX)=trc_solml2_vr(idg_NH3,0,NY,NX)-RBGCSinkSoluteM_vr(idg_NH3,0,NY,NX)
!   exclude nutrients in band
    DO nsol=ids_nut_beg,ids_nuts_end
      trc_solml2_vr(nsol,0,NY,NX)=trc_solml2_vr(nsol,0,NY,NX)-RBGCSinkSoluteM_vr(nsol,0,NY,NX)
    ENDDO
    RBGCSinkGasMM_vr(idg_O2,0,NY,NX)=RO2UptkSoilM_vr(M,0,NY,NX)*dt_GasCyc
  ENDIF

  DO idg=idg_beg,idg_NH3
    if(idg/=idg_NH3)then
      trc_solml2_vr(idg,0,NY,NX)=trc_solml2_vr(idg,0,NY,NX)-RBGCSinkGasMM_vr(idg,0,NY,NX)
    ELSE
      trc_gasml2_vr(idg_NH3,0,NY,NX)=trc_gasml2_vr(idg_NH3,0,NY,NX)-RBGCSinkGasMM_vr(idg_NH3,0,NY,NX)
    ENDIF
  ENDDO

!
!     INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
!
  IF(M.NE.MX)THEN
    DO  L=1,JS
      trcg_TBLS_snvr(idg_beg:idg_NH3,L,NY,NX)     = 0._r8
      trcn_TBLS(ids_nut_beg:ids_nuts_end,L,NY,NX) = 0._r8
    ENDDO

    DO L=NU(NY,NX),NL(NY,NX)    
      DO  K=1,jcplx
        DOM_Transp2Micp_vr(idom_beg:idom_end,K,L,NY,NX) = 0.0_r8
        DOM_Transp2Macp_flx(idom_beg:idom_end,K,L,NY,NX)= 0.0_r8
      ENDDO

      TR3MicPoreSolFlx_vr(ids_beg:ids_end,L,NY,NX) = 0.0_r8
      TR3MacPoreSolFlx_vr(ids_beg:ids_end,L,NY,NX) = 0._r8
!
!     ADD SOLUTE SINKS
!
!     R*SK2=total flux from nitro.f, uptake.f, solute.f
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band micropores
!     ZNH4B,ZNH3B,ZNO3B,ZNO2B,H1POB,H2POB=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band micropores
!
!     include NH3 and band nutrients
      DO nsol=ids_nuts_beg,ids_nuts_end
        trc_solml2_vr(nsol,L,NY,NX)=trc_solml2_vr(nsol,L,NY,NX)-RBGCSinkSoluteM_vr(nsol,L,NY,NX)
      ENDDO
      !add oxygen uptake here
      RBGCSinkGasMM_vr(idg_O2,L,NY,NX)=RO2UptkSoilM_vr(M,L,NY,NX)*dt_GasCyc
    ENDDO  
  ENDIF
!
!     INITIALIZE SOIL SOLUTE NET FLUX ACCUMULATORS
!
  DO L=NU(NY,NX),NL(NY,NX)
!
!     SOIL GAS FLUX ACCUMULATORS
    Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,L,NY,NX)=0._r8
    DO idg=idg_beg,idg_NH3
      if(idg/=idg_NH3)then
        trc_solml2_vr(idg,L,NY,NX)=trc_solml2_vr(idg,L,NY,NX)-RBGCSinkGasMM_vr(idg,L,NY,NX)
      else
        trc_gasml2_vr(idg_NH3,L,NY,NX)=trc_gasml2_vr(idg_NH3,L,NY,NX)-RBGCSinkGasMM_vr(idg_NH3,L,NY,NX)
      endif
    ENDDO

  ENDDO
  end subroutine ResetFluxAccumulatorsMM


!------------------------------------------------------------------------------------------

  subroutine TracerExchXGridsMM(M,MX,NY,NX,NHE,NVS,WaterFlow2Soil,RGasAtmDisol2SoilM)
!
! DESCRIPTION
! exchanges tracers within (gaseous vs aqueous phase) and between
! grid cells.
  implicit none

  integer, intent(in) :: M,MX, NY, NX, NHE, NVS
  real(r8), intent(inout) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8), intent(in) :: RGasAtmDisol2SoilM(idg_beg:idg_end)
  real(r8) :: VFLW
  real(r8) :: VOLH2A,VOLH2B
  real(r8) :: VLWatMacPS,VOLWT

  integer :: iFlagEbu,N,L,K,LL
  integer :: N1,N2,N3  !source grid index
  integer :: N4,N5,N6  !dest grid index

! begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  iFlagEbu=0
  D125: DO L=1,NL(NY,NX)
    !source
    N1=NX;N2=NY;N3=L    
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
    D120: DO N=FlowDirIndicator(N2,N1),3

      IF(N.EQ.iEastWestDirection)THEN
        !WEST-EAST
        IF(NX.EQ.NHE)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N6=L
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        !NORTH-SOUTH
        IF(NY.EQ.NVS)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N6=L
        ENDIF
      ELSEIF(N.EQ.iVerticalDirection)THEN
        !VERTICAL
        IF(L.EQ.NL(NY,NX))THEN
          cycle
        ELSE
          N4 = NX
          N5 = NY
          N6 = L+1
        ENDIF
      ENDIF

      DO LL=N6,NL(NY,NX)
        IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
          N6=LL
          exit
        ENDIF
      ENDDO
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
!     WATER CONTEnsolutes AND WATER FLUXES 'WaterFlow2Soil' FROM 'WATSUB'
!
!     VLSoilPoreMicP_vr,VLSoilMicP=soil volume excluding rock, macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     VLWatMicPM,VLWatMacPM=micropore,macropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!     ReductVLsoiAirPM=change in air volume
!     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
!
      IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
        IF(N3.GE.NUM(N2,N1) .AND. N6.GE.NUM(N5,N4) .AND. N3.LE.NL(N2,N1) .AND. N6.LE.NL(N5,N4))THEN

          IF(M.NE.MX)THEN
            VLWatMicPMA_vr(N6,N5,N4) = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
            VLWatMicPMB_vr(N6,N5,N4) = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
            VLWatMicPXA_vr(N6,N5,N4) = natomw*VLWatMicPMA_vr(N6,N5,N4)
            VLWatMicPXB_vr(N6,N5,N4)    = natomw*VLWatMicPMB_vr(N6,N5,N4)

            VLsoiAirPMA(N6,N5,N4)         = VLsoiAirPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
            VLsoiAirPMB(N6,N5,N4)         = VLsoiAirPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
            CumReductVLsoiAirPM(N6,N5,N4) = ReductVLsoiAirPM(M,N6,N5,N4)*dt_GasCyc
!
            WaterFlow2Soil(N,N6,N5,N4)=(WaterFlow2MicPM_3D(M,N,N6,N5,N4) &
              +WaterFlow2MacPM_3D(M,N,N6,N5,N4))*dt_GasCyc
!
            call SoluteAdvDifusTranspM(M,N,N1,N2,N3,N4,N5,N6)
!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
            IF(N.EQ.iVerticalDirection)THEN
              call MicMacPoresSoluteXchangeM(M,N,N1,N2,N3,N4,N5,N6)
            ENDIF
          ENDIF
!
!     GASEOUS TRANSPORT FROM GASEOUS DIFFUSIVITY AND CONCENTRATION
!     DIFFERENCES BETWEEN ADJACENT GRID CELLS
!
          call GasTransportMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)

        ELSEIF(N.NE.iVerticalDirection)THEN
          call ZeroTransport1MM(N,N1,N2,N3,N4,N5,N6)
        ENDIF
      ELSE
        call ZeroTransport2MM(N,N1,N2,N3,N4,N5,N6)
      ENDIF
    ENDDO D120
!
!     CHECK FOR BUBBLING IF THE SUM OF ALL GASEOUS EQUIVALENT
!     PARTIAL CONCENTRATIONS EXCEEDS ATMOSPHERIC PRESSURE
    if(M.NE.MX)call BubbleEffluxM(M,N1,N2,N3,NY,NX,RGasAtmDisol2SoilM,iFlagEbu)

  ENDDO D125
  end subroutine TracerExchXGridsMM

! ----------------------------------------------------------------------

   subroutine MicroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Adv2MicP_flx)
   implicit none
   integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
   real(r8),intent(out) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)

   real(r8) :: trcs_adv_flx(ids_beg:ids_end)
   integer :: K,nsol,idom
   real(r8) :: VFLW
   
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     WaterFlow2MicPM_3D=water flux through soil micropore from watsub.f
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     RFL*S=solute diffusive flux through micropore
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *S2,*B2=micropore solute content in non-band,band
!
  IF(WaterFlow2MicPM_3D(M,N,N6,N5,N4).GT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,N6,N5,N4)/VLWatMicPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2(idom,K,N3,N2,N1))
      enddo
    ENDDO

    DO nsol=ids_beg,ids_end
      trcs_adv_flx(nsol)=VFLW*AZMAX1(trc_solml2_vr(nsol,N3,N2,N1))
    ENDDO
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS TO CURRENT FROM
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
  ELSE
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MicPM_3D(M,N,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    D9815: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2(idom,K,N6,N5,N4))
      enddo
    ENDDO D9815
    DO nsol=ids_beg,ids_end
      trcs_adv_flx(nsol)=VFLW*AZMAX1(trc_solml2_vr(nsol,N6,N5,N4))
    ENDDO
  ENDIF

  DO nsol=ids_beg,ids_end
    R3PoreSolFlx_3D(nsol,N,N6,N5,N4)=trcs_adv_flx(nsol)
  ENDDO

  end subroutine MicroporeSoluteAdvectionM

! ----------------------------------------------------------------------
  subroutine MicroporeSoluteDiffusionM(M,N,N1,N2,N3,N4,N5,N6,THETW1_vr,DOM_Difus_Mac2Micp_flxM)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: THETW1_vr(JZ,JY,JX)
  real(r8), intent(out):: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)

  real(r8) :: VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: VLWatMicP2A,VLWatMicP2B,VLWatMicP3A,VLWatMicP3B,VLWatMicP4A,VLWatMicP4B
  real(r8) :: trcsolc1(ids_beg:ids_end)
  real(r8) :: trcsolc2(ids_beg:ids_end)
  real(r8) :: CDOM_MicP1(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MicP2(idom_beg:idom_end,1:jcplx)
  real(r8) :: SDifc(ids_beg:ids_end),SDifFlx(ids_beg:ids_end)
  real(r8) :: DISPN,DIFOM(idom_beg:idom_end)
  real(r8) :: DLYR1,DLYR2,TORTL
  integer  :: K,nsol,idg,idom

  IF(THETW1_vr(N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1) .AND. THETW1_vr(N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4) &
    .AND. VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) .AND. VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN

      VLWatMicP2A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
      VLWatMicP2B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
      VLWatMicP3A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3,N3,N2,N1)
      VLWatMicP3B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3B,N3,N2,N1)
      VLWatMicP4A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4,N3,N2,N1)
      VLWatMicP4B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4B,N3,N2,N1)

      VLWatMicPPA = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
      VLWatMicPPB = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
      VLWatMicPOA = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
      VLWatMicPOB = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)

!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!
!     MICROPORE CONCENTRATIONS FROM WATER-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
!
!     C*1,C*2=solute concentration in source,destination layer
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *S2,*B2=soil solute content in non-band,band
!
    D9810: DO K=1,jcplx
      do idom=idom_beg,idom_end
        !source
        CDOM_MicP1(idom,K)=AZMAX1(DOM_MicP2(idom,K,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
        !dest
        CDOM_MicP2(idom,K)=AZMAX1(DOM_MicP2(idom,K,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
      enddo
    ENDDO D9810

    DO idg=idg_beg,idg_end-2
      trcsolc1(idg)=AZMAX1(trc_solml2_vr(idg,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
    ENDDO

    IF(VLWatMicP4A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4) = AZMAX1(trc_solml2_vr(ids_NH4,N3,N2,N1)/VLWatMicP4A)
      trcsolc1(idg_NH3) = AZMAX1(trc_solml2_vr(idg_NH3,N3,N2,N1)/VLWatMicP4A)
    ELSE
      trcsolc1(ids_NH4)=0.0_r8
      trcsolc1(idg_NH3)=0.0_r8
    ENDIF
    IF(VLWatMicP3A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3) = AZMAX1(trc_solml2_vr(ids_NO3,N3,N2,N1)/VLWatMicP3A)
      trcsolc1(ids_NO2) = AZMAX1(trc_solml2_vr(ids_NO2,N3,N2,N1)/VLWatMicP3A)
    ELSE
      trcsolc1(ids_NO3)=0.0_r8
      trcsolc1(ids_NO2)=0.0_r8
    ENDIF
    IF(VLWatMicP2A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4) = AZMAX1(trc_solml2_vr(ids_H1PO4,N3,N2,N1)/VLWatMicP2A)
      trcsolc1(ids_H2PO4) = AZMAX1(trc_solml2_vr(ids_H2PO4,N3,N2,N1)/VLWatMicP2A)
    ELSE
      trcsolc1(ids_H1PO4)=0.0_r8
      trcsolc1(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VLWatMicP4B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4B) = AZMAX1(trc_solml2_vr(ids_NH4B,N3,N2,N1)/VLWatMicP4B)
      trcsolc1(idg_NH3B) = AZMAX1(trc_solml2_vr(idg_NH3B,N3,N2,N1)/VLWatMicP4B)
    ELSE
      trcsolc1(ids_NH4B)=0.0_r8
      trcsolc1(idg_NH3B)=0.0_r8
    ENDIF
    IF(VLWatMicP3B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3B) = AZMAX1(trc_solml2_vr(ids_NO3B,N3,N2,N1)/VLWatMicP3B)
      trcsolc1(ids_NO2B) = AZMAX1(trc_solml2_vr(ids_NO2B,N3,N2,N1)/VLWatMicP3B)
    ELSE
      trcsolc1(ids_NO3B) = trcsolc1(ids_NO3)
      trcsolc1(ids_NO2B) = trcsolc1(ids_NO2)
    ENDIF
    IF(VLWatMicP2B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4B) = AZMAX1(trc_solml2_vr(ids_H1PO4B,N3,N2,N1)/VLWatMicP2B)
      trcsolc1(ids_H2PO4B) = AZMAX1(trc_solml2_vr(ids_H2PO4B,N3,N2,N1)/VLWatMicP2B)
    ELSE
      trcsolc1(ids_H1PO4B)=trcsolc1(ids_H1PO4)
      trcsolc1(ids_H2PO4B)=trcsolc1(ids_H2PO4)
    ENDIF

    DO idg=idg_beg,idg_end-2
      trcsolc2(idg)=AZMAX1(trc_solml2_vr(idg,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
    ENDDO


    IF(VLWatMicPMA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      trcsolc2(idg_NH3) = AZMAX1(trc_solml2_vr(idg_NH3,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4))
      trcsolc2(ids_NH4) = AZMAX1(trc_solml2_vr(ids_NH4,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4))
    ELSE
      trcsolc2(idg_NH3) = 0.0_r8
      trcsolc2(ids_NH4) = 0.0_r8
    ENDIF
    IF(VLWatMicPOA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3) = AZMAX1(trc_solml2_vr(ids_NO3,N6,N5,N4)/VLWatMicPOA)
      trcsolc2(ids_NO2) = AZMAX1(trc_solml2_vr(ids_NO2,N6,N5,N4)/VLWatMicPOA)
    ELSE
      trcsolc2(ids_NO3) = 0.0_r8
      trcsolc2(ids_NO2) = 0.0_r8
    ENDIF
    IF(VLWatMicPPA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4) = AZMAX1(trc_solml2_vr(ids_H1PO4,N6,N5,N4)/VLWatMicPPA)
      trcsolc2(ids_H2PO4) = AZMAX1(trc_solml2_vr(ids_H2PO4,N6,N5,N4)/VLWatMicPPA)
    ELSE
      trcsolc2(ids_H1PO4) = 0.0_r8
      trcsolc2(ids_H2PO4) = 0.0_r8
    ENDIF
    IF(VLWatMicPMB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      trcsolc2(idg_NH3B) = AZMAX1(trc_solml2_vr(idg_NH3B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4))
      trcsolc2(ids_NH4B) = AZMAX1(trc_solml2_vr(ids_NH4B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4))
    ELSE
      trcsolc2(idg_NH3B) = trcsolc2(idg_NH3)
      trcsolc2(ids_NH4B) = trcsolc2(ids_NH4)
    ENDIF
    IF(VLWatMicPOB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3B) = AZMAX1(trc_solml2_vr(ids_NO3B,N6,N5,N4)/VLWatMicPOB)
      trcsolc2(ids_NO2B) = AZMAX1(trc_solml2_vr(ids_NO2B,N6,N5,N4)/VLWatMicPOB)
    ELSE
      trcsolc2(ids_NO3B) = trcsolc2(ids_NO3)
      trcsolc2(ids_NO2B) = trcsolc2(ids_NO2)
    ENDIF
    IF(VLWatMicPPB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4B) = AZMAX1(trc_solml2_vr(ids_H1PO4B,N6,N5,N4)/VLWatMicPPB)
      trcsolc2(ids_H2PO4B) = AZMAX1(trc_solml2_vr(ids_H2PO4B,N6,N5,N4)/VLWatMicPPB)
    ELSE
      trcsolc2(ids_H1PO4B) = trcsolc2(ids_H1PO4)
      trcsolc2(ids_H2PO4B) = trcsolc2(ids_H2PO4)
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!
!     DLYR=soil layer thickness
!     TORT=micropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     WaterFlow2MicPM_3D=water flux through soil micropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through micropore
!     *SGL2=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     XDPTH=cross-sectional area/distance between layers
!     C*1,C*2=micropore solute concentration in source,destination layer
!     DFV*=diffusive solute transfer through soil micropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMicPM_vr(M,N3,N2,N1)*DLYR1+TortMicPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN = DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MicPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DO idom=idom_beg,idom_end
      DIFOM(idom)=(DOMdiffusivity2_vr(idom,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

    DO nsol=ids_beg,ids_end
      SDifc(nsol)=(SoluteDifusivitytscaled_vr(nsol,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
    D9805: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Mac2Micp_flxM(idom,K)=DIFOM(idom)*(CDOM_MicP1(idom,K)-CDOM_MicP2(idom,K))        
      enddo
    ENDDO D9805

    DO idg=idg_beg,idg_end-2
      SDifFlx(idg)=SDifc(idg)*(trcsolc1(idg)-trcsolc2(idg))
    ENDDO

    DO nsol=ids_nuts_beg,ids_nuts_end
      SDifFlx(nsol)=SDifc(nsol)*(trcsolc1(nsol)-trcsolc2(nsol))*AMIN1(trcs_VLN_vr(nsol,N3,N2,N1),trcs_VLN_vr(nsol,N6,N5,N4))
    ENDDO
  ELSE
    D9905: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9905
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF

  DO nsol=ids_beg,ids_end
    R3PoreSolFlx_3D(nsol,N,N6,N5,N4)=R3PoreSolFlx_3D(nsol,N,N6,N5,N4)+SDifFlx(nsol)
  ENDDO

  end subroutine MicroporeSoluteDiffusionM

! ----------------------------------------------------------------------
  subroutine MacroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Adv2MacP_flxM)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(out):: DOM_Adv2MacP_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: trcs_RFH(ids_beg:ids_end)
  integer  :: K,idg,nsol,idom
  real(r8) :: VFLW


!     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
!

  IF(WaterFlow2MacPM_3D(M,N,N6,N5,N4).GT.0.0_r8)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN CURRENT
!     GRID CELL
!
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     VLWatMicPAH=macropore porosity
!     RFH*=solute diffusive flux through macropore
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *SH2,*BH2=macropore solute content in non-band,band
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    IF(VLWatMacPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
!
!     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE
!
    IF(N.EQ.3.AND.VLMacP_vr(N6,N5,N4).GT.VLWatMacPM_vr(M,N6,N5,N4))THEN
      D9800: DO K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1((DOM_MacP2(idom,K,N3,N2,N1) &
            -AZMIN1(DOM_Mac2MicPore_flxM_vr(idom,K,NU(N2,N1),N2,N1))))
        enddo
      ENDDO D9800

      DO idg=idg_beg,idg_end-2
        trcs_RFH(idg)=VFLW*AZMAX1((trc_soHml2_vr(idg,N3,N2,N1)-AZMIN1(trcs_Mac2MicPore_flxM_vr(idg,NU(N2,N1),N2,N1))))
      ENDDO

      DO nsol=ids_nuts_beg,ids_nuts_end
        trcs_RFH(nsol)=VFLW*AZMAX1((trc_soHml2_vr(nsol,N3,N2,N1) &
          -AZMIN1(trcs_Mac2MicPore_flxM_vr(nsol,NU(N2,N1),N2,N1)*trcs_VLN_vr(nsol,N3,N2,N1)))) &
          *trcs_VLN_vr(nsol,N6,N5,N4)
      ENDDO
!
!     OTHERWISE
!
    ELSE
      D9850: DO K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1(DOM_MacP2(idom,K,N3,N2,N1))
        enddo
      ENDDO D9850
!exclude NH3 and NH3B
      DO idg=idg_beg,idg_end-2
        trcs_RFH(idg)=VFLW*AZMAX1(trc_soHml2_vr(idg,N3,N2,N1))
      ENDDO

      DO nsol=ids_nuts_beg,ids_nuts_end
        trcs_RFH(nsol)=VFLW*AZMAX1(trc_soHml2_vr(nsol,N3,N2,N1))*trcs_VLN_vr(nsol,N6,N5,N4)
      ENDDO
    ENDIF
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN ADJACENT
!     GRID CELL
!
  ELSEIF(WaterFlow2MacPM_3D(M,N,N6,N5,N4).LT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    D9665: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1(DOM_MacP2(idom,K,N6,N5,N4))
      enddo
    ENDDO D9665

    DO idg=idg_beg,idg_end-2
      trcs_RFH(idg)=VFLW*AZMAX1(trc_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO nsol=ids_nuts_beg,ids_nuts_end
      trcs_RFH(nsol)=VFLW*AZMAX1(trc_soHml2_vr(nsol,N6,N5,N4))*trcs_VLN_vr(nsol,N6,N5,N4)
    ENDDO
  ELSE
!
!     NO MACROPORE FLUX
!
    D9795: DO K=1,jcplx
      DOM_Adv2MacP_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9795
    trcs_RFH(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO nsol=ids_beg,ids_end
    R3PoreSoHFlx_3D(nsol,N,N6,N5,N4)=trcs_RFH(nsol)
  ENDDO

  end subroutine MacroporeSoluteAdvectionM

! ----------------------------------------------------------------------
  subroutine MacroporeSoluteDispersionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Difus_Macp_flxM)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(out):: DOM_Difus_Macp_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: VOLH2A,VOLH2B,VOLH3A,VOLH3B,VOLH4A,VOLH4B
  real(r8) :: VOLHMA,VOLHMB,VOLHOA,VOLHOB,VOLHPA,VOLHPB
  real(r8) :: trcs_coH1(ids_beg:ids_end)
  real(r8) :: trcs_coH2(ids_beg:ids_end)
  real(r8) :: SDifc(ids_beg:ids_end),TORTL
  real(r8) :: DIFOM(idom_beg:idom_end)
  real(r8) :: DISPN,DLYR1,DLYR2
  real(r8) :: SDifHFlx(ids_beg:ids_end)
  integer  :: K,nsol,idg,idom
  real(r8) :: CDOM_MacP1(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MacP2(idom_beg:idom_end,1:jcplx)


!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MACROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     THETY=hygroscopic water content
!     VOLAH=total macropore volume
!
  IF(VLWatMacPM_vr(M,N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1)*VLMacP_vr(N3,N2,N1) &
    .AND. VLWatMacPM_vr(M,N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4)*VLMacP_vr(N6,N5,N4))THEN
!
!     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
!
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     *H2=macropore solute content
!     VLWatMacPM=macropore water content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
    VOLH4A=VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4,N3,N2,N1)
    VOLH4B=VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4B,N3,N2,N1)
    VOLH3A=VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3,N3,N2,N1)
    VOLH3B=VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3B,N3,N2,N1)
    VOLH2A=VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    VOLH2B=VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)

    VOLHOA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
    VOLHOB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)
    VOLHPA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    VOLHPB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)

    D9790: DO K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MacP1(idom,K)=AZMAX1(DOM_MacP2(idom,K,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
        CDOM_MacP2(idom,K)=AZMAX1(DOM_MacP2(idom,K,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
      enddo
    ENDDO D9790
    !exclude NH3 and NH3B
    DO idg=idg_beg,idg_end-2
      trcs_coH1(idg)=AZMAX1(trc_soHml2_vr(idg,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
    ENDDO

    IF(VOLH4A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4)=AZMAX1(trc_soHml2_vr(ids_NH4,N3,N2,N1)/VOLH4A)
      trcs_coH1(idg_NH3)=AZMAX1(trc_soHml2_vr(idg_NH3,N3,N2,N1)/VOLH4A)
    ELSE
      trcs_coH1(ids_NH4)=0.0_r8
      trcs_coH1(idg_NH3)=0.0_r8
    ENDIF
    IF(VOLH3A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3)=AZMAX1(trc_soHml2_vr(ids_NO3,N3,N2,N1)/VOLH3A)
      trcs_coH1(ids_NO2)=AZMAX1(trc_soHml2_vr(ids_NO2,N3,N2,N1)/VOLH3A)
    ELSE
      trcs_coH1(ids_NO3)=0.0_r8
      trcs_coH1(ids_NO2)=0.0_r8
    ENDIF
    IF(VOLH2A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4)=AZMAX1(trc_soHml2_vr(ids_H1PO4,N3,N2,N1)/VOLH2A)
      trcs_coH1(ids_H2PO4)=AZMAX1(trc_soHml2_vr(ids_H2PO4,N3,N2,N1)/VOLH2A)
    ELSE
      trcs_coH1(ids_H1PO4)=0.0_r8
      trcs_coH1(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VOLH4B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4B)=AZMAX1(trc_soHml2_vr(ids_NH4B,N3,N2,N1)/VOLH4B)
      trcs_coH1(idg_NH3B)=AZMAX1(trc_soHml2_vr(idg_NH3B,N3,N2,N1)/VOLH4B)
    ELSE
      trcs_coH1(ids_NH4B)=trcs_coH1(ids_NH4)
      trcs_coH1(idg_NH3B)=trcs_coH1(idg_NH3)
    ENDIF
    IF(VOLH3B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3B)=AZMAX1(trc_soHml2_vr(ids_NO3B,N3,N2,N1)/VOLH3B)
      trcs_coH1(ids_NO2B)=AZMAX1(trc_soHml2_vr(ids_NO2B,N3,N2,N1)/VOLH3B)
    ELSE
      trcs_coH1(ids_NO3B)=trcs_coH1(ids_NO3)
      trcs_coH1(ids_NO2B)=trcs_coH1(ids_NO2)
    ENDIF
    IF(VOLH2B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4B)=AZMAX1(trc_soHml2_vr(ids_H1PO4B,N3,N2,N1)/VOLH2B)
      trcs_coH1(ids_H2PO4B)=AZMAX1(trc_soHml2_vr(ids_H2PO4B,N3,N2,N1)/VOLH2B)
    ELSE
      trcs_coH1(ids_H1PO4B)=trcs_coH1(ids_H1PO4)
      trcs_coH1(ids_H2PO4B)=trcs_coH1(ids_H2PO4)
    ENDIF
   !excldue NH3 and NH3B
    DO idg=idg_beg,idg_end-2
      trcs_coH2(idg)=AZMAX1(trc_soHml2_vr(idg,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
    ENDDO

    VOLHMA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
    IF(VOLHMA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4)=AZMAX1(trc_soHml2_vr(ids_NH4,N6,N5,N4)/VOLHMA)
      trcs_coH2(idg_NH3)=AZMAX1(trc_soHml2_vr(idg_NH3,N6,N5,N4)/VOLHMA)
    ELSE
      trcs_coH2(ids_NH4)=0.0_r8
      trcs_coH2(idg_NH3)=0.0_r8
    ENDIF
    VOLHOA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
    IF(VOLHOA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3)=AZMAX1(trc_soHml2_vr(ids_NO3,N6,N5,N4)/VOLHOA)
      trcs_coH2(ids_NO2)=AZMAX1(trc_soHml2_vr(ids_NO2,N6,N5,N4)/VOLHOA)
    ELSE
      trcs_coH2(ids_NO3)=0.0_r8
      trcs_coH2(ids_NO2)=0.0_r8
    ENDIF
    VOLHPA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    IF(VOLHPA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4)=AZMAX1(trc_soHml2_vr(ids_H1PO4,N6,N5,N4)/VOLHPA)
      trcs_coH2(ids_H2PO4)=AZMAX1(trc_soHml2_vr(ids_H2PO4,N6,N5,N4)/VOLHPA)
    ELSE
      trcs_coH2(ids_H1PO4)=0.0_r8
      trcs_coH2(ids_H2PO4)=0.0_r8
    ENDIF
    VOLHMB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
    IF(VOLHMB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4B)=AZMAX1(trc_soHml2_vr(ids_NH4B,N6,N5,N4)/VOLHMB)
      trcs_coH2(idg_NH3B)=AZMAX1(trc_soHml2_vr(idg_NH3B,N6,N5,N4)/VOLHMB)
    ELSE
      trcs_coH2(ids_NH4B)=trcs_coH2(ids_NH4)
      trcs_coH2(idg_NH3B)=trcs_coH2(idg_NH3)
    ENDIF
    VOLHOB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)
    IF(VOLHOB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3B)=AZMAX1(trc_soHml2_vr(ids_NO3B,N6,N5,N4)/VOLHOB)
      trcs_coH2(ids_NO2B)=AZMAX1(trc_soHml2_vr(ids_NO2B,N6,N5,N4)/VOLHOB)
    ELSE
      trcs_coH2(ids_NO3B)=trcs_coH2(ids_NO3)
      trcs_coH2(ids_NO2B)=trcs_coH2(ids_NO2)
    ENDIF
    VOLHPB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    IF(VOLHPB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4B)=AZMAX1(trc_soHml2_vr(ids_H1PO4B,N6,N5,N4)/VOLHPB)
      trcs_coH2(ids_H2PO4B)=AZMAX1(trc_soHml2_vr(ids_H2PO4B,N6,N5,N4)/VOLHPB)
    ELSE
      trcs_coH2(ids_H1PO4B)=trcs_coH2(ids_H1PO4)
      trcs_coH2(ids_H2PO4B)=trcs_coH2(ids_H2PO4)
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
!
!     DLYR=soil layer thickness
!     TortMacPM_vr=macropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through macropore
!     *SGL2=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     XDPTH=cross-sectional area/distance between layers
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     DFH*=diffusive solute transfer through soil macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMacPM_vr(M,N3,N2,N1)*DLYR1+TortMacPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN = DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MacPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DO idom=idom_beg,idom_end
      DIFOM(idom)=(DOMdiffusivity2_vr(idom,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

    DO nsol=ids_beg,ids_end
      SDifc(nsol)=(SoluteDifusivitytscaled_vr(nsol,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MACROPORES
!
    D9785: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Macp_flxM(idom,K)=DIFOM(idom)*(CDOM_MacP1(idom,K)-CDOM_MacP2(idom,K))
      enddo
    ENDDO D9785
! exclude NH3 and NH3B
    DO idg=idg_beg,idg_end-2
      SDifHFlx(idg)=SDifc(idg)*(trcs_coH1(idg)-trcs_coH2(idg))
    ENDDO

    DO nsol=ids_nuts_beg,ids_end
      SDifHFlx(nsol)=SDifc(nsol)*(trcs_coH1(nsol)-trcs_coH2(nsol)) &
        *AMIN1(trcs_VLN_vr(nsol,N3,N2,N1),trcs_VLN_vr(nsol,N6,N5,N4))
    ENDDO
  ELSE
    D9780: DO K=1,jcplx
      DOM_Difus_Macp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9780
    SDifHFlx(ids_beg:ids_end)=0.0_r8
  ENDIF
  DO nsol = ids_beg,ids_end
    R3PoreSoHFlx_3D(nsol,N,N6,N5,N4)=R3PoreSoHFlx_3D(nsol,N,N6,N5,N4)+SDifHFlx(nsol)
  ENDDO
  end subroutine MacroporeSoluteDispersionM

! ----------------------------------------------------------------------
  subroutine SoluteAdvDifusTranspM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N3  !source grid
  integer, intent(in) :: N4,N5,N6  !dest grid

  real(r8)  :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)  
  real(r8) :: THETW1_vr(JZ,JY,JX)  !soil saturation 
  integer  :: K,nsol,idom
  real(r8) :: DOM_Adv2MacP_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: DOM_Difus_Macp_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)

  THETW1_vr(N3,N2,N1)=AZMAX1(safe_adb(VLWatMicPM_vr(M,N3,N2,N1),VLSoilMicP_vr(N3,N2,N1)))
  THETW1_vr(N6,N5,N4)=AZMAX1(safe_adb(VLWatMicPM_vr(M,N6,N5,N4),VLSoilMicP_vr(N6,N5,N4)))

!     SOLUTE TRANSPORT IN MICROPORES
!
  call MicroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Adv2MicP_flx)
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MICROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  call MicroporeSoluteDiffusionM(M,N,N1,N2,N3,N4,N5,N6,THETW1_vr,DOM_Difus_Mac2Micp_flxM)
!
!     SOLUTE TRANSPORT IN MACROPORES
!
  call MacroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Adv2MacP_flxM)
!
  call MacroporeSoluteDispersionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Difus_Macp_flxM)
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES
!
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,R*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,R*FHB=convective + diffusive solute flux through macropores in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!
  D9765: DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)=DOM_Adv2MicP_flx(idom,K)+DOM_Difus_Mac2Micp_flxM(idom,K)
      DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)=DOM_Adv2MacP_flxM(idom,K)+DOM_Difus_Macp_flxM(idom,K)
    enddo
  ENDDO D9765

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux through micropores
!     X*FLW,X*FLB= hourly convective + diffusive solute flux through micropores in non-band,band
!     X*FHS=hourly convective + diffusive solute flux through macropores
!     X*FHW,X*FHB= hourly convective + diffusive solute flux through macropores in non-band,band
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,X*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,X*FHB=convective + diffusive solute flux through macropores in non-band,band
!
  D9755: DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTransp_3D(idom,K,N,N6,N5,N4)=DOM_MicpTransp_3D(idom,K,N,N6,N5,N4) &
        +DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)
      DOM_3DMacp_Transp_flx(idom,K,N,N6,N5,N4)=DOM_3DMacp_Transp_flx(idom,K,N,N6,N5,N4) &
        +DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)
    enddo
  ENDDO D9755

  DO nsol=ids_beg,ids_end
    trcs_TransptMicP_3D(nsol,N,N6,N5,N4)=trcs_TransptMicP_3D(nsol,N,N6,N5,N4)+R3PoreSolFlx_3D(nsol,N,N6,N5,N4)
    trcs_TransptMacP_3D(nsol,N,N6,N5,N4)=trcs_TransptMacP_3D(nsol,N,N6,N5,N4)+R3PoreSoHFlx_3D(nsol,N,N6,N5,N4)
  ENDDO

  end subroutine SoluteAdvDifusTranspM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteXchangeM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  integer :: K,nsol,idg,idom

  do idg=idg_beg,idg_NH3-1
    trcg_VLWatMicP_vr(idg,N6,N5,N4)=VLWatMicPM_vr(M,N6,N5,N4)*GasSolbility_vr(idg,N6,N5,N4)
  enddo

  trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4)=VLWatMicPMA_vr(N6,N5,N4)*GasSolbility_vr(idg_NH3,N6,N5,N4)
  trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4)=VLWatMicPMB_vr(N6,N5,N4)*GasSolbility_vr(idg_NH3,N6,N5,N4)
!
!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND

  call MicMacPoresSoluteAdvExchange(M,N,N1,N2,N3,N4,N5,N6)

!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION
!     DIFFERENCES
  call MicMacPoresSoluteDifExchange(M,N,N1,N2,N3,N4,N5,N6)

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!
  D9945: DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Mac2MicPore_flx_vr(idom,K,N6,N5,N4)=DOM_Mac2MicPore_flx_vr(idom,K,N6,N5,N4)+DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)
    enddo
  ENDDO D9945

  DO nsol=ids_beg,ids_end
    trcs_Mac2MicPore_flx_vr(nsol,N6,N5,N4)=trcs_Mac2MicPore_flx_vr(nsol,N6,N5,N4)+trcs_Mac2MicPore_flxM_vr(nsol,N6,N5,N4)
  ENDDO

  end subroutine MicMacPoresSoluteXchangeM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteAdvExchange(M,N,N1,N2,N3,N4,N5,N6)
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)
  real(r8) :: VFLW    !fraction in flow
  real(r8) :: trcs_adv_flx(ids_beg:ids_end)
  integer :: K,idg,nsol,idom
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FWatExMacP2MicPM_vr=macro-micropore water transfer from watsub.f
!     VLWatMicPM,VLWatMacPM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FWatExMacP2MicPM_vr(M,N6,N5,N4).GT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF
    D9970: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MacP2(idom,K,N6,N5,N4))
      enddo
    ENDDO D9970

    DO idg=idg_beg,idg_end-2
      trcs_adv_flx(idg)=VFLW*AZMAX1(trc_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO nsol=ids_nuts_beg,ids_nuts_end
      trcs_adv_flx(nsol)=VFLW*AZMAX1(trc_soHml2_vr(nsol,N6,N5,N4))*trcs_VLN_vr(nsol,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FWatExMacP2MicPM_vr(M,N6,N5,N4).LT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    D9965: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2(idom,K,N6,N5,N4))
      enddo
    ENDDO D9965
!exclude NH3 and NH3B
    DO idg=idg_beg,idg_end-2
      trcs_adv_flx(idg)=VFLW*AZMAX1(trc_solml2_vr(idg,N6,N5,N4))
    ENDDO

    DO nsol=ids_nuts_beg,ids_nuts_end
      trcs_adv_flx(nsol)=VFLW*AZMAX1(trc_solml2_vr(nsol,N6,N5,N4))*trcs_VLN_vr(nsol,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    D9960: DO K=1,jcplx
      DOM_Adv2MicP_flx(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9960
    trcs_adv_flx(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO nsol=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(nsol,N6,N5,N4)=trcs_adv_flx(nsol)
  ENDDO
!
!     TOTAL CONVECTIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
!     R*FXS,R*FXB=total convective + diffusive solute flux between macro- and micropore in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     RFL*=convective flux between macro- and micropore
!     DFV*=diffusive solute flux between macro- and micropore
!

  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)=DOM_Adv2MicP_flx(idom,K)
    enddo
  enddo

  end subroutine MicMacPoresSoluteAdvExchange

! ----------------------------------------------------------------------

  subroutine MicMacPoresSoluteDifExchange(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: trcs_DFV(ids_beg:ids_end)
  real(r8) :: VLWatMacPS,VOLWT
  integer  :: K,nsol,idg,idom
  real(r8) :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)

!
!     VLWatMicPM,VLWatMacPM=micropore,macropore water-filled porosity from watsub.f
!     DFV*S,DFV*B=diffusive solute flux between macro- and micropore in non-band,band
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *2,*H2=solute content of micropores,macropores
!
  IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
    VLWatMacPS=AMIN1(XFRS*VGeomLayer_vr(N6,N5,N4),VLWatMacPM_vr(M,N6,N5,N4))
    VOLWT=VLWatMicPM_vr(M,N6,N5,N4)+VLWatMacPS

    D9955: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Mac2Micp_flxM(idom,K)=dts_HeatWatTP*(AZMAX1(DOM_MacP2(idom,K,N6,N5,N4)) &
          *VLWatMicPM_vr(M,N6,N5,N4)-AZMAX1(DOM_MicP2(idom,K,N6,N5,N4))*VLWatMacPS)/VOLWT
      enddo
    ENDDO D9955

    DO idg=idg_beg,idg_NH3-1
      trcs_DFV(idg)=dts_HeatWatTP*(AZMAX1(trc_soHml2_vr(idg,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trc_solml2_vr(idg,N6,N5,N4))*VLWatMacPS)/VOLWT
    ENDDO

    DO nsol=ids_nuts_beg,ids_nuts_end
      trcs_DFV(nsol)=dts_HeatWatTP*(AZMAX1(trc_soHml2_vr(nsol,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trc_solml2_vr(nsol,N6,N5,N4))*VLWatMacPS)/VOLWT &
        *trcs_VLN_vr(nsol,N6,N5,N4)
    ENDDO

  ELSE
    D9975: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9975

    trcs_DFV(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO nsol=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(nsol,N6,N5,N4)=trcs_Mac2MicPore_flxM_vr(nsol,N6,N5,N4)+trcs_DFV(nsol)
  ENDDO
!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
!     R*FXS,R*FXB=total convective + diffusive solute flux between macro- and micropore in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     RFL*=convective flux between macro- and micropore
!     DFV*=diffusive solute flux between macro- and micropore
!

  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)=DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)+DOM_Difus_Mac2Micp_flxM(idom,K)
    enddo
  enddo
  end subroutine MicMacPoresSoluteDifExchange

! ----------------------------------------------------------------------
  subroutine BubbleEffluxM(M,N1,N2,N3,NY,NX,RGasAtmDisol2SoilM,iFlagEbu)
  implicit none
  integer, intent(in) :: M,N1,N2,N3,NY,NX
  real(r8), intent(in) :: RGasAtmDisol2SoilM(idg_beg:idg_end)
  integer, intent(inout) :: iFlagEbu
  real(r8) :: THETW1
  real(r8) :: trcg_SLX(idg_beg:idg_end)
  real(r8) :: trcg_VOLG(idg_beg:idg_end)
  integer  :: idg
  real(r8) :: VTATM,VTGAS,DVTGAS
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     VLSoilMicP=micropore volume
!     iFlagEbu=bubbling flag:0=enabled,1=disabled
!     S*L=solubility of gas in water from hour1.f
!
  IF(N3.GE.NUM(N2,N1))THEN
    THETW1=AZMAX1(safe_adb(VLWatMicPM_vr(M,N3,N2,N1),VLSoilMicP_vr(N3,N2,N1)))
    IF(THETW1.GT.SoilWatAirDry_vr(N3,N2,N1).AND.iFlagEbu.EQ.0)THEN

      trcg_SLX(idg_CO2) =catomw*GasSolbility_vr(idg_CO2,N3,N2,N1)  !conver into carbon g C/mol
      trcg_SLX(idg_CH4) =catomw*GasSolbility_vr(idg_CH4,N3,N2,N1)
      trcg_SLX(idg_O2)  =32.0_r8*GasSolbility_vr(idg_O2,N3,N2,N1)
      trcg_SLX(idg_N2)  =28.0_r8*GasSolbility_vr(idg_N2,N3,N2,N1)
      trcg_SLX(idg_N2O) =28.0_r8*GasSolbility_vr(idg_N2O,N3,N2,N1)
      trcg_SLX(idg_NH3) =natomw*GasSolbility_vr(idg_NH3,N3,N2,N1)
      trcg_SLX(idg_H2)  =2.0_r8*GasSolbility_vr(idg_H2,N3,N2,N1)
      trcg_SLX(idg_NH3B)=trcg_SLX(idg_NH3)
!
!     GASEOUS EQUIVALENT PARTIAL CONCENTRATIONS
!
!     V*G2=molar gas concentration
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     R*DFS=gas exchange between atmosphere and soil surface water
!
      IF(N3.EQ.NU(N2,N1))THEN
        DO idg=idg_beg,idg_end
          trcg_VOLG(idg)=(trc_solml2_vr(idg,N3,N2,N1)+RGasAtmDisol2SoilM(idg))/trcg_SLX(idg)
        ENDDO
      ELSE
        DO idg=idg_beg,idg_end
          trcg_VOLG(idg)=trc_solml2_vr(idg,N3,N2,N1)/trcg_SLX(idg)
        ENDDO
      ENDIF
!
!     GASEOUS EQUIVALENT ATMOSPHERIC CONCENTRATION
!
!     VTATM=molar gas concentration at atmospheric pressure
!     VTGAS=total molar gas concentration
!
      VTATM=AZMAX1(1.2194E+04_r8*VLWatMicPM_vr(M,N3,N2,N1)/TKS_vr(N3,N2,N1))

      VTGAS=sum(trcg_VOLG(idg_beg:idg_end))
!
!     PROPORTIONAL REMOVAL OF EXCESS AQUEOUS GASES
!
!     R*BBL=bubble flux
!     gas code:*CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     V*G2=molar gas concentration
!
      IF(VTGAS.GT.VTATM)THEN
        DVTGAS=0.5_r8*(VTATM-VTGAS)  !<0, bubbling occurs
        DO idg=idg_beg,idg_end
          trcg_Ebu_vr(idg,N3,N2,N1)=AZMIN1(DVTGAS*trcg_VOLG(idg)/VTGAS)*trcg_SLX(idg)  
          IF(trcg_Ebu_vr(idg,N3,N2,N1)<1.e-10_r8)trcg_Ebu_vr(idg,N3,N2,N1)=0._r8
        ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*BBL=hourly bubble flux
!
        DO idg=idg_beg,idg_end
          trcg_ebu_flx_vr(idg,N3,N2,N1)=trcg_ebu_flx_vr(idg,N3,N2,N1)+trcg_Ebu_vr(idg,N3,N2,N1) !<0, ebullition
        ENDDO
      ELSE
        trcg_Ebu_vr(idg_beg:idg_end,N3,N2,N1)=0.0_r8
      ENDIF
    ELSE
      iFlagEbu=1
      trcg_Ebu_vr(idg_beg:idg_end,N3,N2,N1)=0.0_r8
    ENDIF

  ENDIF
  end subroutine BubbleEffluxM
! ----------------------------------------------------------------------

  subroutine GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: trc_gasc1(idg_beg:idg_end)
  real(r8) :: trc_gasc2(idg_beg:idg_end)

  real(r8) :: CNDC1
  real(r8) :: CNDC2
  real(r8) :: DFLG2,DFLGL
  integer :: idg

!     GASEOUS DIFFUSIVITIES
!
!     DFLG2,DFLGL=air-filled porosity effect on gaseous diffusivity in source,destination layer
!     POROQ=Penman Water Linear Reduction tortuosity from starts.f
!     POROS=total porosity
!     DLYR=soil layer thickness
!     D*G=gaseous diffusivity in soil
!     CND*1,CND*2=gaseous conductance in source,destination layer
!     *SGL2= gaseous diffusivity in air
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
  DFLG2=2.0_r8*AZMAX1(AirFilledSoilPoreM_vr(M,N3,N2,N1))*POROQ*AirFilledSoilPoreM_vr(M,N3,N2,N1)/POROS_vr(N3,N2,N1) &
    *AREA(N,N3,N2,N1)/DLYR_3D(N,N3,N2,N1)

  DFLGL=2.0_r8*AZMAX1(AirFilledSoilPoreM_vr(M,N6,N5,N4))*POROQ*AirFilledSoilPoreM_vr(M,N6,N5,N4)/POROS_vr(N6,N5,N4) &
    *AREA(N,N6,N5,N4)/DLYR_3D(N,N6,N5,N4)

!
!     GASOUS CONDUCTANCES
!
!     D*G=gaseous diffusivity in soil
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
  DO idg=idg_beg,idg_end
    CNDC1=DFLG2*GasDifctScaled_vr(idg,N3,N2,N1)
    CNDC2=DFLGL*GasDifctScaled_vr(idg,N6,N5,N4)
    DifuscG_vr(idg,N,N6,N5,N4)=(CNDC1*CNDC2)/(CNDC1+CNDC2)
  ENDDO
!
!     GASEOUS CONCENTRATIONS FROM AIR-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
!
!     C*G1,C*G2=gaseous concentration in source,destination layer
!     *G2=gaseous content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     VLsoiAirPM=air-filled porosity
!
!
!     DIFFUSIVE GAS TRANSFER DRIVEN BY GAS CONCENTRATIONS IN
!     ADJACENT GRID CELLS
!
!     DFV*G=diffusive gas flux
!     C*G1,C*G2=gaseous concentration in source,destination layer
!
! does not include band NH3
  DO idg=idg_beg,idg_NH3
    trc_gasc1(idg)               = AZMAX1(trc_gasml2_vr(idg,N3,N2,N1)/VLsoiAirPM_vr(M,N3,N2,N1))
    trc_gasc2(idg)               = AZMAX1(trc_gasml2_vr(idg,N6,N5,N4)/VLsoiAirPM_vr(M,N6,N5,N4))
    RGasADFlxMM_3D(idg,N,N6,N5,N4) = DifuscG_vr(idg,N,N6,N5,N4)*(trc_gasc1(idg)-trc_gasc2(idg))
  ENDDO

  end subroutine GasDiffusionMM


! ----------------------------------------------------------------------
  subroutine GasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)
  implicit none
  real(r8), intent(in) :: WaterFlow2Soil(3,JD,JV,JH)
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: RGasAdv
  real(r8) :: VFLW,FLQW
  integer :: idg
!
!     CONVECTIVE GAS TRANSFER DRIVEN BY SOIL WATER FLUXES
!     FROM 'WATSUB' AND GAS CONCENTRATIONS IN THE ADJACENT GRID CELLS
!     DEPENDING ON WATER FLUX DIRECTION
!
!     by assuming volume conservation, gases and water flow in opposite direction
!     WaterFlow2Soil=total water flux into soil micropore+macropore from watsub.f
!     VLsoiAirPM=air-filled porosity
!     RFL*G=convective gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     *G2=gaseous content
!
  FLQW=WaterFlow2Soil(N,N6,N5,N4)
  IF(FLQW.GT.0.0_r8)THEN
    IF(VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=-AZMAX1(AMIN1(VFLWX,FLQW/VLsoiAirPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO idg=idg_beg,idg_NH3
      RGasAdv                       = VFLW*AZMAX1(trc_gasml2_vr(idg,N6,N5,N4))
      RGasADFlxMM_3D(idg,N,N6,N5,N4) = RGasADFlxMM_3D(idg,N,N6,N5,N4)+RGasAdv
    ENDDO
  ELSE
    IF(VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=-AZMIN1(AMAX1(-VFLWX,FLQW/VLsoiAirPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO idg=idg_beg,idg_NH3
      RGasAdv                       = VFLW*AZMAX1(trc_gasml2_vr(idg,N3,N2,N1))
      RGasADFlxMM_3D(idg,N,N6,N5,N4) = RGasADFlxMM_3D(idg,N,N6,N5,N4)+RGasAdv
    ENDDO
  ENDIF

  end subroutine GasAdvectionMM

! ----------------------------------------------------------------------
  subroutine GasTransportMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)

  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: WaterFlow2Soil(3,JD,JV,JH)
  integer :: idg

!     AirFilledSoilPoreM_vr,VLsoiAirPM=air-filled porosity,volume from watsub.f

  IF(AirFilledSoilPoreM_vr(M,N3,N2,N1).GT.THETX .AND. AirFilledSoilPoreM_vr(M,N6,N5,N4).GT.THETX &
    .AND.VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
    .AND.VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN

!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE
    call GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

!     TOTAL SOIL GAS FLUX FROM CONVECTIVE FLUX
    call GasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLG=hourly total convective+diffusive gas flux
!
    DO idg=idg_beg,idg_NH3
      Gas_3DAdvDif_Flx_vr(idg,N,N6,N5,N4)=Gas_3DAdvDif_Flx_vr(idg,N,N6,N5,N4)+RGasADFlxMM_3D(idg,N,N6,N5,N4)
    ENDDO

  ELSE
    RGasADFlxMM_3D(idg_beg:idg_NH3,N,N6,N5,N4)=0.0_r8
  ENDIF

  call GasDissolutionMM(M,N,N1,N2,N3,N4,N5,N6)

  end subroutine GasTransportMM

! ----------------------------------------------------------------------
  subroutine GasDissolutionMM(M,N,N1,N2,N3,N4,N5,N6)
  !
  !Description
  !Compute gas dissolution 
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N3,N4,N5,N6

  integer :: idg
!
!     VOLATILIZATION-DISSOLUTION OF GASES IN SOIL
!     LAYER FROM GASEOUS CONCENTRATIONS VS. THEIR AQUEOUS
!     EQUIVALEnsolutes DEPENDING ON SOLUBILITY FROM 'HOUR1'
!     AND TRANSFER COEFFICIENT 'DiffusivitySolutEff' FROM 'WATSUB'
!
!     AirFilledSoilPoreM_vr,VLWatMicPPM=air-filled porosity,volume
!     R*DFG=water-air gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     DiffusivitySolutEff=rate constant for air-water gas exchange from watsub.f
!     *G2,*S2=gaseous,aqueous gas content
!     VLWatMicP*=equivalent aqueous volume for gas
!
  IF(N.EQ.iVerticalDirection)THEN
    IF(AirFilledSoilPoreM_vr(M,N6,N5,N4).GT.THETX)THEN
      do idg=idg_beg,idg_NH3-1
        RGas_Disol_FlxMM_vr(idg,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg,N6,N5,N4))*trcg_VLWatMicP_vr(idg,N6,N5,N4) &
          -trc_solml2_vr(idg,N6,N5,N4)*VLsoiAirPM_vr(M,N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg,N6,N5,N4)+VLsoiAirPM_vr(M,N6,N5,N4))
      enddo    

      IF(VLsoiAirPMA(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGas_Disol_FlxMM_vr(idg_NH3,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg_NH3,N6,N5,N4))*trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4) &
          -trc_solml2_vr(idg_NH3,N6,N5,N4)*VLsoiAirPMA(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4)+VLsoiAirPMA(N6,N5,N4))

      ELSE
        RGas_Disol_FlxMM_vr(idg_NH3,N6,N5,N4)=0.0_r8
      ENDIF

      IF(VLsoiAirPMB(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGas_Disol_FlxMM_vr(idg_NH3B,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
          (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg_NH3,N6,N5,N4))*trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4) &
          -trc_solml2_vr(idg_NH3B,N6,N5,N4)*VLsoiAirPMB(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4)+VLsoiAirPMB(N6,N5,N4))

      ELSE
        RGas_Disol_FlxMM_vr(idg_NH3B,N6,N5,N4)=0.0_r8
      ENDIF

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly water-air gas flux
!
      DO idg=idg_beg,idg_end
        Gas_Disol_Flx_vr(idg,N6,N5,N4)=Gas_Disol_Flx_vr(idg,N6,N5,N4)+RGas_Disol_FlxMM_vr(idg,N6,N5,N4)
      ENDDO
    ELSE
      RGas_Disol_FlxMM_vr(idg_beg:idg_end,N6,N5,N4)=0.0_r8
    ENDIF
  ENDIF
  end subroutine GasDissolutionMM

! ----------------------------------------------------------------------
  subroutine ZeroTransport1MM(N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,N6

  integer :: K

  DifuscG_vr(idg_beg:idg_end,N,N6,N5,N4) = 0._r8

  D9750: DO K=1,jcplx
    DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,N6,N5,N4)=0.0_r8
    DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,N6,N5,N4)=0.0_r8
  ENDDO D9750

  R3PoreSolFlx_3D(ids_beg:ids_end,N,N6,N5,N4)=0.0_r8

  R3PoreSoHFlx_3D(ids_beg:ids_end,N,N6,N5,N4)=0.0_r8

  RGasADFlxMM_3D(idg_beg:idg_NH3,N,N6,N5,N4)=0.0_r8
  end subroutine ZeroTransport1MM
! ----------------------------------------------------------------------

  subroutine ZeroTransport2MM(N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,N6

  integer :: K

  DifuscG_vr(idg_beg:idg_NH3,N,N3,N2,N1)=0.0_r8

  D9751: DO K=1,jcplx
    DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,N3,N2,N1) = 0.0_r8
    DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,N3,N2,N1) = 0.0_r8
  ENDDO D9751

  R3PoreSolFlx_3D(ids_beg:ids_end,N,N3,N2,N1) = 0.0_r8
  R3PoreSoHFlx_3D(ids_beg:ids_end,N,N3,N2,N1) = 0.0_r8
  RGasADFlxMM_3D(idg_beg:idg_NH3,N,N3,N2,N1)    = 0.0_r8
  end subroutine ZeroTransport2MM
end module InsideTranspMod
