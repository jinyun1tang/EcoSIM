module InsideTranspMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use abortutils,        only: destroy,  endrun
  use minimathmod,       only: safe_adb, AZMAX1, AZMIN1
  use TracerPropMod,     only: MolecularWeight
  use PlantDataRateType, only: trcs_deadroot2soil_vr
  use DebugToolMod  
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

  subroutine ModelTracerHydroFluxMM(I,J,M,MX,NHW, NHE, NVN, NVS,WaterFlow2SoilMM_3D,RGasAtmDisol2LitrM,RGasAtmDisol2SoilM)
  !
  !Description
  !Hydrological tracer flux in iteration MM
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,MX
  integer, intent(in) :: NHW, NHE, NVN, NVS
  real(r8), intent(inout) :: WaterFlow2SoilMM_3D(3,JD,JV,JH)
  real(r8), intent(out) :: RGasAtmDisol2LitrM(idg_beg:idg_NH3,JY,JX)  !atmospheric gas dissolves to litter solutes
  real(r8), intent(out) :: RGasAtmDisol2SoilM(idg_beg:idg_end,JY,JX)  !atmospheric gas dissolves to soil solutes
  
  character(len=*), parameter :: subname='ModelTracerHydroFluxMM' 
  integer :: NY,NX
  real(r8) :: FLWRM1
  real(r8) :: RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)
  real(r8) :: RGas_Dif_Atm2Litr_FlxMM(idg_beg:idg_NH3)

  call PrintInfo('beg '//subname)

  RGas_Dif_Atm2Soil_FlxMM(idg_beg:idg_end)=0._r8

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      call ResetFluxAccumulatorsMM(I,J,M,NY,NX,MX)

      call ApplyTracerBGCSinkMM(I,J,M,NY,NX,MX)

      IF(M.NE.MX)THEN
        !this is a new iteration for moisture-temp update
 
        call SnowSoluteDischargeM(I,J,M,NY,NX)
!
!     SOLUTE FLUXES AT SOIL SURFACE FROM SURFACE WATER
!     CONTEnsolutes, WATER FLUXES 'WaterFlow2SoilMM_3D' AND ATMOSPHERE BOUNDARY
!     LAYER RESISTANCES 'PARGM' FROM 'WATSUB'
!
        call SoluteFluxSurfaceM(I,J,M,NY,NX,NHE,NHW,NVS,NVN, WaterFlow2SoilMM_3D,RGas_Dif_Atm2Soil_FlxMM,&
          RGas_Dif_Atm2Litr_FlxMM,RGasAtmDisol2SoilM(:,NY,NX),RGasAtmDisol2LitrM(:,NY,NX))
!
      ENDIF
!
!     VOLATILIZATION-DISSOLUTION OF GASES IN RESIDUE AND SOIL SURFACE
!     LAYERS FROM GASEOUS CONCENTRATIONS VS. THEIR AQUEOUS
!
      call LitterGasVolatilDissolMM(I,J,M,NY,NX,RGas_Dif_Atm2Litr_FlxMM)
!
!     SURFACE GAS EXCHANGE FROM GAS DIFFUSIVITY THROUGH
!     SOIL SURFACE LAYER AND THROUGH ATMOSPHERE BOUNDARY
!     LAYER
!
      call SurfSoilFluxGasDifAdvMM(M,NY,NX,WaterFlow2SoilMM_3D,RGas_Dif_Atm2Soil_FlxMM)
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
!
      call TracerExchXGridsMM(I,J,M,MX,NY,NX,NHE,NVS,WaterFlow2SoilMM_3D,RGasAtmDisol2SoilM(:,NY,NX))

    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine ModelTracerHydroFluxMM
!------------------------------------------------------------------------------------------

  subroutine ResetFluxAccumulatorsMM(I,J,M,NY,NX,MX)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M, NY, NX,MX

  integer :: K,L,idg,ids,idom
  real(r8) :: PARGM

  IF(M.NE.MX)THEN
    !
    !     GASEOUS BOUNDARY LAYER CONDUCTANCES
    !
    !     PARG=boundary layer conductance above soil surface from watsub.f
    !
    PARGM                   = CondGasXSnowM_col(M,NY,NX)*dt_GasCyc
    PARGas_CefMM(idg_CO2,NY,NX) = PARGM*0.74_r8
    PARGas_CefMM(idg_CH4,NY,NX) = PARGM*1.04_r8
    PARGas_CefMM(idg_O2,NY,NX)  = PARGM*0.83_r8
    PARGas_CefMM(idg_N2,NY,NX)  = PARGM*0.86_r8
    PARGas_CefMM(idg_N2O,NY,NX) = PARGM*0.74_r8
    PARGas_CefMM(idg_NH3,NY,NX) = PARGM*1.02_r8
    PARGas_CefMM(idg_H2,NY,NX)  = PARGM*2.08_r8
    PARGas_CefMM(idg_AR,NY,NX)  = PARGM*0.72_r8    
!    
!     RESET RUNOFF SOLUTE FLUX ACCUMULATORS
!

    DO  K=1,jcplx
      DOM_SurfRunoff_flxM(idom_beg:idom_end,K,NY,NX)=0.0_r8
    ENDDO

    trcg_SurfRunoff_flx(idg_beg:idg_NH3,NY,NX)          = 0.0_r8
    trcn_SurfRunoff_flx(ids_nut_beg:ids_nuts_end,NY,NX) = 0.0_r8
    trcg_SnowDrift_flxM(idg_beg:idg_NH3,NY,NX)           = 0.0_r8
    trcn_SnowDrift_flxM(ids_nut_beg:ids_nuts_end,NY,NX)  = 0.0_r8
  ENDIF

!
  IF(M.NE.MX)THEN
!   INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
    DO  L=1,JS
      trcg_Aqua_flxM_snvr(idg_beg:idg_NH3,L,NY,NX)          = 0._r8
      trcn_Aqua_flxM_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = 0._r8
    ENDDO
    !fluxes in soil
    DO L=NU(NY,NX),NL(NY,NX)    
      DO  K=1,jcplx
        DOM_Transp2Micp_flxM_vr(idom_beg:idom_end,K,L,NY,NX) = 0.0_r8
        DOM_Transp2Macp_flxM_vr(idom_beg:idom_end,K,L,NY,NX)= 0.0_r8
      ENDDO
      trcs_Transp2Micp_flxM_vr(ids_beg:ids_end,L,NY,NX) = 0.0_r8
      trcs_Transp2Macp_flxM_vr(ids_beg:ids_end,L,NY,NX) = 0._r8
    ENDDO  
  ENDIF
!
!     SOIL GAS FLUX ACCUMULATORS
!
  DO L=NU(NY,NX),NL(NY,NX)

    DO idg=idg_beg,idg_NH3-1
      Gas_AdvDif_FlxMM_vr(idg,L,NY,NX) = 0._r8
    ENDDO    
    Gas_AdvDif_FlxMM_vr(idg_NH3,L,NY,NX) = 0._r8
    
  ENDDO
  end subroutine ResetFluxAccumulatorsMM
!------------------------------------------------------------------------------------------
  subroutine ApplyTracerBGCSinkMM(I,J,M,NY,NX,MX)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M, NY, NX,MX

  integer :: L,ids,idg
  logical :: lbegIterM   !flag for the begining of iteration M
  
  lbegIterM=M.NE.MX  
  IF(lbegIterM)THEN
    !first round of iteration M
    ! because NH3 is gas-aqua dual phase
    trcs_solml2_vr(idg_NH3,0,NY,NX)=trcs_solml2_vr(idg_NH3,0,NY,NX)-RBGCSinkSoluteM_vr(idg_NH3,0,NY,NX)

   ! exclude nutrients in band
    DO ids=ids_nut_beg,ids_nuts_end
      trcs_solml2_vr(ids,0,NY,NX)=trcs_solml2_vr(ids,0,NY,NX)-RBGCSinkSoluteM_vr(ids,0,NY,NX)
    ENDDO
    RBGCSinkGasMM_vr(idg_O2,0,NY,NX)=RO2UptkSoilM_vr(M,0,NY,NX)*dt_GasCyc

    ! ADD SOLUTE SINKS
    DO L=NU(NY,NX),NL(NY,NX)    
      ! include NH3 and band nutrients
      DO ids=ids_nuts_beg,ids_nuts_end
        trcs_solml2_vr(ids,L,NY,NX)=trcs_solml2_vr(ids,L,NY,NX)-RBGCSinkSoluteM_vr(ids,L,NY,NX)
      ENDDO
      !first round of iteration M
      !add oxygen uptake here
      RBGCSinkGasMM_vr(idg_O2,L,NY,NX)=RO2UptkSoilM_vr(M,L,NY,NX)*dt_GasCyc-trcs_deadroot2soil_vr(idg_O2,L,NY,NX)*dts_gas
    ENDDO
  ENDIF

  trc_gasml2_vr(idg_NH3,0,NY,NX)=AZMAX1(trc_gasml2_vr(idg_NH3,0,NY,NX)-RBGCSinkGasMM_vr(idg_NH3,0,NY,NX))  
  DO idg=idg_beg,idg_NH3-1
    trcs_solml2_vr(idg,0,NY,NX)=trcs_solml2_vr(idg,0,NY,NX)-RBGCSinkGasMM_vr(idg,0,NY,NX)
  ENDDO

  DO L=NU(NY,NX),NL(NY,NX)
    !gaseous NH3 is consumed through geochemistry
    trc_gasml2_vr(idg_NH3,L,NY,NX)  = AZMAX1(trc_gasml2_vr(idg_NH3,L,NY,NX)-RBGCSinkGasMM_vr(idg_NH3,L,NY,NX))    
    
    !other gases are taken from aqueous phases
    DO idg=idg_beg,idg_NH3-1
      trcs_solml2_vr(idg,L,NY,NX) = AZMAX1(trcs_solml2_vr(idg,L,NY,NX)-RBGCSinkGasMM_vr(idg,L,NY,NX))
    ENDDO

  ENDDO
  end subroutine ApplyTracerBGCSinkMM

!------------------------------------------------------------------------------------------

  subroutine TracerExchXGridsMM(I,J,M,MX,NY,NX,NHE,NVS,WaterFlow2SoilMM_3D,RGasAtmDisol2SoilM)
!
! DESCRIPTION
! exchanges tracers within (gaseous vs aqueous phase) and between
! grid cells.
  implicit none

  integer, intent(in) :: I,J,M,MX
  integer, intent(in) :: NY, NX, NHE, NVS
  real(r8), intent(inout) :: WaterFlow2SoilMM_3D(3,JD,JV,JH)
  real(r8), intent(in) :: RGasAtmDisol2SoilM(idg_beg:idg_end)

  character(len=*), parameter :: subname='TracerExchXGridsMM'
  real(r8) :: VFLW
  real(r8) :: VOLH2A,VOLH2B
  real(r8) :: VLWatMacPS,VOLWT

  integer :: iDisableEbu,N,L,K,LL
  integer :: N1,N2,N3  !source grid index
  integer :: N4,N5,N6  !dest grid index

! begin_execution
  call PrintInfo('beg '//subname)

!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  iDisableEbu=ifalse
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
          N4 = NX+1
          N5 = NY
          N6 = L
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        !NORTH-SOUTH
        IF(NY.EQ.NVS)THEN
          cycle
        ELSE
          N4 = NX
          N5 = NY+1
          N6 = L
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
!     WATER CONTEnsolutes AND WATER FLUXES 'WaterFlow2SoilMM_3D' FROM 'WATSUB'

!     ReductVLsoiAirPM=change in air volume
!     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
!
      IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
        IF(N3.GE.NUM(N2,N1) .AND. N6.GE.NUM(N5,N4) .AND. N3.LE.NL(N2,N1) .AND. N6.LE.NL(N5,N4))THEN

          IF(M.NE.MX)THEN
            VLWatMicPMA_vr(N6,N5,N4) = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
            VLWatMicPMB_vr(N6,N5,N4) = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
            VLWatMicPXA_vr(N6,N5,N4) = natomw*VLWatMicPMA_vr(N6,N5,N4)
            VLWatMicPXB_vr(N6,N5,N4) = natomw*VLWatMicPMB_vr(N6,N5,N4)

            VLsoiAirPMA_vr(N6,N5,N4)          = VLsoiAirPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
            VLsoiAirPMB_vr(N6,N5,N4)          = VLsoiAirPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
            CumReductVLsoiAirPMM_vr(N6,N5,N4) = ReductVLsoiAirPM_vr(M,N6,N5,N4)*dt_GasCyc
!
            WaterFlow2SoilMM_3D(N,N6,N5,N4)=(WaterFlow2MicPM_3D(M,N,N6,N5,N4) &
              +WaterFlow2MacPM_3D(M,N,N6,N5,N4))*dt_GasCyc
!
            call SoluteAdvDifusTranspM(I,J,M,N,N1,N2,N3,N4,N5,N6)
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
          call GasTransportMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2SoilMM_3D)

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
    if(M.NE.MX)call BubbleEffluxM(M,N1,N2,N3,NY,NX,RGasAtmDisol2SoilM,iDisableEbu)
  ENDDO D125
  call PrintInfo('end '//subname)
  end subroutine TracerExchXGridsMM

! ----------------------------------------------------------------------

   subroutine MicroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Adv2MicP_flx)
   implicit none
   integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
   real(r8),intent(out) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)

   real(r8) :: trcs_adv_flx(ids_beg:ids_end)
   integer :: K,ids,idom
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

    DO ids=ids_beg,ids_end
      trcs_adv_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N3,N2,N1))
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
    DO ids=ids_beg,ids_end
      trcs_adv_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))
    ENDDO
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)=trcs_adv_flx(ids)
    if(abs(trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4))>1.e10)then
      write(*,*)ids,N,N6,N5,N4,trcs_adv_flx(ids)
      write(*,*)VFLW,trcs_solml2_vr(ids,N6,N5,N4),trcs_solml2_vr(ids,N3,N2,N1),trcs_names(ids)
      write(*,*)trcs_solml_vr(ids,N6,N5,N4),trcs_solml_vr(ids,N3,N2,N1)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDDO

  end subroutine MicroporeSoluteAdvectionM

! ----------------------------------------------------------------------
  subroutine MicroporeSoluteDiffusionM(M,N,N1,N2,N3,N4,N5,N6,THETW1_vr,DOM_Difus_Mac2Micp_flxM)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: THETW1_vr(JZ,JY,JX)
  real(r8), intent(out):: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)

  character(len=*), parameter :: subname='MicroporeSoluteDiffusionM'
  real(r8) :: VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: VLWatMicP2A,VLWatMicP2B,VLWatMicP3A,VLWatMicP3B,VLWatMicP4A,VLWatMicP4B
  real(r8) :: trcsolc1(ids_beg:ids_end)
  real(r8) :: trcsolc2(ids_beg:ids_end)
  real(r8) :: CDOM_MicP1(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MicP2(idom_beg:idom_end,1:jcplx)
  real(r8) :: SDifc(ids_beg:ids_end),SDifFlx(ids_beg:ids_end)
  real(r8) :: DISPN,DIFOM(idom_beg:idom_end)
  real(r8) :: DLYR1,DLYR2,TORTL
  integer  :: K,ids,idg,idom

  call PrintInfo('beg '//subname)
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

    DO idg=idg_beg,idg_NH3-1
      trcsolc1(idg)=AZMAX1(trcs_solml2_vr(idg,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
    ENDDO

    IF(VLWatMicP4A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4,N3,N2,N1)/VLWatMicP4A),tracerSolc_max(ids_NH4))
      trcsolc1(idg_NH3) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3,N3,N2,N1)/VLWatMicP4A),tracerSolc_max(idg_NH3))

    ELSE
      trcsolc1(ids_NH4)=0.0_r8
      trcsolc1(idg_NH3)=0.0_r8
    ENDIF
    IF(VLWatMicP3A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3,N3,N2,N1)/VLWatMicP3A),tracerSolc_max(ids_NO3))
      trcsolc1(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,N3,N2,N1)/VLWatMicP3A)      
    ELSE
      trcsolc1(ids_NO3)=0.0_r8
      trcsolc1(ids_NO2)=0.0_r8
    ENDIF
    IF(VLWatMicP2A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,N3,N2,N1)/VLWatMicP2A)
      trcsolc1(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,N3,N2,N1)/VLWatMicP2A)
    ELSE
      trcsolc1(ids_H1PO4)=0.0_r8
      trcsolc1(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VLWatMicP4B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4B,N3,N2,N1)/VLWatMicP4B),tracerSolc_max(ids_NH4B))
      trcsolc1(idg_NH3B) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3B,N3,N2,N1)/VLWatMicP4B),tracerSolc_max(idg_NH3B))
    ELSE
      trcsolc1(ids_NH4B)=0.0_r8
      trcsolc1(idg_NH3B)=0.0_r8
    ENDIF
    IF(VLWatMicP3B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3B,N3,N2,N1)/VLWatMicP3B),tracerSolc_max(ids_NO3B))
      trcsolc1(ids_NO2B) = AZMAX1(trcs_solml2_vr(ids_NO2B,N3,N2,N1)/VLWatMicP3B)
    ELSE
      trcsolc1(ids_NO3B) = trcsolc1(ids_NO3)
      trcsolc1(ids_NO2B) = trcsolc1(ids_NO2)
    ENDIF
    IF(VLWatMicP2B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,N3,N2,N1)/VLWatMicP2B)
      trcsolc1(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,N3,N2,N1)/VLWatMicP2B)
    ELSE
      trcsolc1(ids_H1PO4B)=trcsolc1(ids_H1PO4)
      trcsolc1(ids_H2PO4B)=trcsolc1(ids_H2PO4)
    ENDIF

    DO idg=idg_beg,idg_NH3-1
      trcsolc2(idg)=AZMAX1(trcs_solml2_vr(idg,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
    ENDDO


    IF(VLWatMicPMA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      !gN/m3 H2O, the maximum solubility is about 9.7gN/100g water at 25oC
      !
      trcsolc2(idg_NH3) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4)),tracerSolc_max(idg_NH3))
      trcsolc2(ids_NH4) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4)),tracerSolc_max(ids_NH4))

    ELSE
      trcsolc2(idg_NH3) = 0.0_r8
      trcsolc2(ids_NH4) = 0.0_r8
    ENDIF
    IF(VLWatMicPOA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3) = AZMAX1(trcs_solml2_vr(ids_NO3,N6,N5,N4)/VLWatMicPOA)
      trcsolc2(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,N6,N5,N4)/VLWatMicPOA)
    ELSE
      trcsolc2(ids_NO3) = 0.0_r8
      trcsolc2(ids_NO2) = 0.0_r8
    ENDIF
    IF(VLWatMicPPA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,N6,N5,N4)/VLWatMicPPA)
      trcsolc2(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,N6,N5,N4)/VLWatMicPPA)
    ELSE
      trcsolc2(ids_H1PO4) = 0.0_r8
      trcsolc2(ids_H2PO4) = 0.0_r8
    ENDIF
    IF(VLWatMicPMB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      trcsolc2(idg_NH3B) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4)),tracerSolc_max(idg_NH3B))
      trcsolc2(ids_NH4B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4)),tracerSolc_max(ids_NH4B))
    ELSE
      trcsolc2(idg_NH3B) = trcsolc2(idg_NH3)
      trcsolc2(ids_NH4B) = trcsolc2(ids_NH4)
    ENDIF
    IF(VLWatMicPOB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3B,N6,N5,N4)/VLWatMicPOB),tracerSolc_max(ids_NO3B))
      trcsolc2(ids_NO2B) = AZMAX1(trcs_solml2_vr(ids_NO2B,N6,N5,N4)/VLWatMicPOB)
    ELSE
      trcsolc2(ids_NO3B) = trcsolc2(ids_NO3)
      trcsolc2(ids_NO2B) = trcsolc2(ids_NO2)
    ENDIF
    IF(VLWatMicPPB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,N6,N5,N4)/VLWatMicPPB)
      trcsolc2(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,N6,N5,N4)/VLWatMicPPB)
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
    DISPN = DISP_3D(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MicPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DO idom=idom_beg,idom_end
      DIFOM(idom)=(DOM_diffusivitytscaledM_vr(idom,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

    DO ids=ids_beg,ids_end
      SDifc(ids)=(SoluteDifusivitytscaledM_vr(ids,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
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

    DO idg=idg_beg,idg_NH3-1
      SDifFlx(idg)=SDifc(idg)*(trcsolc1(idg)-trcsolc2(idg))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      SDifFlx(ids)=SDifc(ids)*(trcsolc1(ids)-trcsolc2(ids))*AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))
      
    ENDDO
  ELSE
    D9905: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9905
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)=trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)+SDifFlx(ids)
    if(abs(trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4))>1.e10)then
      write(*,*)ids,N,N6,N5,N4,trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4),SDifFlx(ids),trcs_names(ids)   
      write(*,*)SDifc(ids),trcsolc1(ids),trcsolc2(ids),AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))   
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine MicroporeSoluteDiffusionM

! ----------------------------------------------------------------------
  subroutine MacroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6,DOM_Adv2MacP_flxM)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(out):: DOM_Adv2MacP_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: trcs_RFH(ids_beg:ids_end)
  integer  :: K,idg,ids,idom
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
        trcs_RFH(idg)=VFLW*AZMAX1((trcs_soHml2_vr(idg,N3,N2,N1)-AZMIN1(trcs_Mac2MicPore_flxM_vr(idg,NU(N2,N1),N2,N1))))
      ENDDO

      DO ids=ids_nuts_beg,ids_nuts_end
        trcs_RFH(ids)=VFLW*AZMAX1((trcs_soHml2_vr(ids,N3,N2,N1) &
          -AZMIN1(trcs_Mac2MicPore_flxM_vr(ids,NU(N2,N1),N2,N1)*trcs_VLN_vr(ids,N3,N2,N1)))) &
          *trcs_VLN_vr(ids,N6,N5,N4)
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
        trcs_RFH(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N3,N2,N1))
      ENDDO

      DO ids=ids_nuts_beg,ids_nuts_end
        trcs_RFH(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N3,N2,N1))*trcs_VLN_vr(ids,N6,N5,N4)
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
      trcs_RFH(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_RFH(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
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

  DO ids=ids_beg,ids_end
    trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)=trcs_RFH(ids)
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
  integer  :: K,ids,idg,idom
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
      trcs_coH1(idg)=AZMAX1(trcs_soHml2_vr(idg,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
    ENDDO

    IF(VOLH4A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4)=AZMAX1(trcs_soHml2_vr(ids_NH4,N3,N2,N1)/VOLH4A)
      trcs_coH1(idg_NH3)=AZMAX1(trcs_soHml2_vr(idg_NH3,N3,N2,N1)/VOLH4A)
    ELSE
      trcs_coH1(ids_NH4)=0.0_r8
      trcs_coH1(idg_NH3)=0.0_r8
    ENDIF
    IF(VOLH3A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3)=AZMAX1(trcs_soHml2_vr(ids_NO3,N3,N2,N1)/VOLH3A)
      trcs_coH1(ids_NO2)=AZMAX1(trcs_soHml2_vr(ids_NO2,N3,N2,N1)/VOLH3A)
    ELSE
      trcs_coH1(ids_NO3)=0.0_r8
      trcs_coH1(ids_NO2)=0.0_r8
    ENDIF
    IF(VOLH2A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4)=AZMAX1(trcs_soHml2_vr(ids_H1PO4,N3,N2,N1)/VOLH2A)
      trcs_coH1(ids_H2PO4)=AZMAX1(trcs_soHml2_vr(ids_H2PO4,N3,N2,N1)/VOLH2A)
    ELSE
      trcs_coH1(ids_H1PO4)=0.0_r8
      trcs_coH1(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VOLH4B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4B)=AZMAX1(trcs_soHml2_vr(ids_NH4B,N3,N2,N1)/VOLH4B)
      trcs_coH1(idg_NH3B)=AZMAX1(trcs_soHml2_vr(idg_NH3B,N3,N2,N1)/VOLH4B)
    ELSE
      trcs_coH1(ids_NH4B)=trcs_coH1(ids_NH4)
      trcs_coH1(idg_NH3B)=trcs_coH1(idg_NH3)
    ENDIF
    IF(VOLH3B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3B)=AZMAX1(trcs_soHml2_vr(ids_NO3B,N3,N2,N1)/VOLH3B)
      trcs_coH1(ids_NO2B)=AZMAX1(trcs_soHml2_vr(ids_NO2B,N3,N2,N1)/VOLH3B)
    ELSE
      trcs_coH1(ids_NO3B)=trcs_coH1(ids_NO3)
      trcs_coH1(ids_NO2B)=trcs_coH1(ids_NO2)
    ENDIF
    IF(VOLH2B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4B)=AZMAX1(trcs_soHml2_vr(ids_H1PO4B,N3,N2,N1)/VOLH2B)
      trcs_coH1(ids_H2PO4B)=AZMAX1(trcs_soHml2_vr(ids_H2PO4B,N3,N2,N1)/VOLH2B)
    ELSE
      trcs_coH1(ids_H1PO4B)=trcs_coH1(ids_H1PO4)
      trcs_coH1(ids_H2PO4B)=trcs_coH1(ids_H2PO4)
    ENDIF
   !excldue NH3 and NH3B
    DO idg=idg_beg,idg_end-2
      trcs_coH2(idg)=AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
    ENDDO

    VOLHMA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
    IF(VOLHMA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4)=AZMAX1(trcs_soHml2_vr(ids_NH4,N6,N5,N4)/VOLHMA)
      trcs_coH2(idg_NH3)=AZMAX1(trcs_soHml2_vr(idg_NH3,N6,N5,N4)/VOLHMA)
    ELSE
      trcs_coH2(ids_NH4)=0.0_r8
      trcs_coH2(idg_NH3)=0.0_r8
    ENDIF
    VOLHOA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
    IF(VOLHOA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3)=AZMAX1(trcs_soHml2_vr(ids_NO3,N6,N5,N4)/VOLHOA)
      trcs_coH2(ids_NO2)=AZMAX1(trcs_soHml2_vr(ids_NO2,N6,N5,N4)/VOLHOA)
    ELSE
      trcs_coH2(ids_NO3)=0.0_r8
      trcs_coH2(ids_NO2)=0.0_r8
    ENDIF
    VOLHPA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    IF(VOLHPA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4)=AZMAX1(trcs_soHml2_vr(ids_H1PO4,N6,N5,N4)/VOLHPA)
      trcs_coH2(ids_H2PO4)=AZMAX1(trcs_soHml2_vr(ids_H2PO4,N6,N5,N4)/VOLHPA)
    ELSE
      trcs_coH2(ids_H1PO4)=0.0_r8
      trcs_coH2(ids_H2PO4)=0.0_r8
    ENDIF
    VOLHMB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
    IF(VOLHMB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4B)=AZMAX1(trcs_soHml2_vr(ids_NH4B,N6,N5,N4)/VOLHMB)
      trcs_coH2(idg_NH3B)=AZMAX1(trcs_soHml2_vr(idg_NH3B,N6,N5,N4)/VOLHMB)
    ELSE
      trcs_coH2(ids_NH4B)=trcs_coH2(ids_NH4)
      trcs_coH2(idg_NH3B)=trcs_coH2(idg_NH3)
    ENDIF
    VOLHOB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)
    IF(VOLHOB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3B)=AZMAX1(trcs_soHml2_vr(ids_NO3B,N6,N5,N4)/VOLHOB)
      trcs_coH2(ids_NO2B)=AZMAX1(trcs_soHml2_vr(ids_NO2B,N6,N5,N4)/VOLHOB)
    ELSE
      trcs_coH2(ids_NO3B)=trcs_coH2(ids_NO3)
      trcs_coH2(ids_NO2B)=trcs_coH2(ids_NO2)
    ENDIF
    VOLHPB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    IF(VOLHPB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4B)=AZMAX1(trcs_soHml2_vr(ids_H1PO4B,N6,N5,N4)/VOLHPB)
      trcs_coH2(ids_H2PO4B)=AZMAX1(trcs_soHml2_vr(ids_H2PO4B,N6,N5,N4)/VOLHPB)
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
    DISPN = DISP_3D(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MacPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DO idom=idom_beg,idom_end
      DIFOM(idom)=(DOM_diffusivitytscaledM_vr(idom,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

    DO ids=ids_beg,ids_end
      SDifc(ids)=(SoluteDifusivitytscaledM_vr(ids,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
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

    DO ids=ids_nuts_beg,ids_end
      SDifHFlx(ids)=SDifc(ids)*(trcs_coH1(ids)-trcs_coH2(ids)) &
        *AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))
    ENDDO
  ELSE
    D9780: DO K=1,jcplx
      DOM_Difus_Macp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9780
    SDifHFlx(ids_beg:ids_end)=0.0_r8
  ENDIF
  DO ids = ids_beg,ids_end
    trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)=trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)+SDifHFlx(ids)
  ENDDO
  end subroutine MacroporeSoluteDispersionM

! ----------------------------------------------------------------------
  subroutine SoluteAdvDifusTranspM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: I,J,M,N
  integer, intent(in) :: N1,N2,N3  !source grid
  integer, intent(in) :: N4,N5,N6  !dest grid
  character(len=*), parameter :: subname='SoluteAdvDifusTranspM'
  real(r8)  :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)  
  real(r8) :: THETW1_vr(JZ,JY,JX)  !soil saturation 
  integer  :: K,ids,idom
  real(r8) :: DOM_Adv2MacP_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: DOM_Difus_Macp_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)

  call PrintInfo('beg '//subname)
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
!
  D9755: DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTransp_3D(idom,K,N,N6,N5,N4)=DOM_MicpTransp_3D(idom,K,N,N6,N5,N4) &
        +DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)
      DOM_Macp_Transp_flx_3D(idom,K,N,N6,N5,N4)=DOM_Macp_Transp_flx_3D(idom,K,N,N6,N5,N4) &
        +DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)
    enddo
  ENDDO D9755

  DO ids=ids_beg,ids_end
    trcs_TransptMicP_3D(ids,N,N6,N5,N4)=trcs_TransptMicP_3D(ids,N,N6,N5,N4)+trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)
    trcs_TransptMacP_3D(ids,N,N6,N5,N4)=trcs_TransptMacP_3D(ids,N,N6,N5,N4)+trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)
    if(abs(trcs_TransptMicP_3D(ids,N,N6,N5,N4))>1.e10)then
      write(*,*)(I*1000+J)*100+M,ids,N,N6,N5,N4,trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4),trcs_names(ids)       
      call endrun(trim(mod_filename)//' at line',__LINE__)         
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SoluteAdvDifusTranspM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteXchangeM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  integer :: K,ids,idg,idom

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

  DO ids=ids_beg,ids_end
    trcs_Mac2MicPore_flx_vr(ids,N6,N5,N4)=trcs_Mac2MicPore_flx_vr(ids,N6,N5,N4)+trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)
  ENDDO

  end subroutine MicMacPoresSoluteXchangeM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteAdvExchange(M,N,N1,N2,N3,N4,N5,N6)
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)
  real(r8) :: VFLW    !fraction in flow
  real(r8) :: trcs_adv_flx(ids_beg:ids_end)
  integer :: K,idg,ids,idom
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
      trcs_adv_flx(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_adv_flx(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
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
      trcs_adv_flx(idg)=VFLW*AZMAX1(trcs_solml2_vr(idg,N6,N5,N4))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_adv_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
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

  DO ids=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)=trcs_adv_flx(ids)
  ENDDO
!
!     TOTAL CONVECTIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
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
  integer  :: K,ids,idg,idom
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
      trcs_DFV(idg)=dts_HeatWatTP*(AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcs_solml2_vr(idg,N6,N5,N4))*VLWatMacPS)/VOLWT
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_DFV(ids)=dts_HeatWatTP*(AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))*VLWatMacPS)/VOLWT &
        *trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO

  ELSE
    D9975: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9975

    trcs_DFV(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)=trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)+trcs_DFV(ids)
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
  subroutine BubbleEffluxM(M,N1,N2,N3,NY,NX,RGasAtmDisol2SoilM,iDisableEbu)
  implicit none
  integer, intent(in) :: M,N1,N2,N3,NY,NX
  real(r8), intent(in) :: RGasAtmDisol2SoilM(idg_beg:idg_end)  !newly dissolved gas 
  integer, intent(inout) :: iDisableEbu                        !bubbling flag

  character(len=*), parameter :: subname='BubbleEffluxM'  
  real(r8) :: THETW1
  real(r8) :: GasMassSolubility(idg_beg:idg_end)
  real(r8) :: trcg_VOLG(idg_beg:idg_end)   !gas concentation corresponding to each volatile [mol d-2]
  integer  :: idg
  real(r8) :: VTATM,VTGAS,DVTGAS
  real(r8) :: FracEbu                      !fraction removed through ebullition
  real(r8) :: dPond
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     VLSoilMicP=micropore volume
!     S*L=solubility of gas in water from hour1.f
!
  call PrintInfo('beg '//subname)
  dPond=0._r8
  IF(N3.GE.NUM(N2,N1))THEN
    THETW1=AZMAX1(safe_adb(VLWatMicPM_vr(M,N3,N2,N1),VLSoilMicP_vr(N3,N2,N1)))
    IF(THETW1.GT.SoilWatAirDry_vr(N3,N2,N1) .AND. iDisableEbu.EQ.ifalse)THEN

      do idg=idg_beg,idg_NH3
        GasMassSolubility(idg) =MolecularWeight(idg)*GasSolbility_vr(idg,N3,N2,N1)  !conver into carbon g /mol
      enddo        
      GasMassSolubility(idg_NH3B)=GasMassSolubility(idg_NH3)
!
!     GASEOUS EQUIVALENT PARTIAL CONCENTRATIONS
!
      IF(N3.EQ.NU(N2,N1))THEN
        DO idg=idg_beg,idg_end
          trcg_VOLG(idg)=AZMAX1(trcs_solml2_vr(idg,N3,N2,N1)+RGasAtmDisol2SoilM(idg))/GasMassSolubility(idg)    !mol/d2 
        ENDDO
      ELSE
        DO idg=idg_beg,idg_end
          trcg_VOLG(idg)=AZMAX1(trcs_solml2_vr(idg,N3,N2,N1))/GasMassSolubility(idg)
        ENDDO
      ENDIF
!
!     GASEOUS EQUIVALENT ATMOSPHERIC CONCENTRATION
!
!     VTATM=molar gas concentration at atmospheric pressure
!     VTGAS=total molar gas concentration     
!    
      if(iPondFlag_col(N2,N1))then 
        if(SoilBulkDensity_vr(N3,N2,N1).LE.ZERO)then
          !still in the ponding water, \rho*g*h,  1.e3*10
          dPond=(CumDepz2LayBottom_vr(N3,N2,N1)-CumDepz2LayBottom_vr(NUM(N2,N1)-1,N2,N1))*10._r8
        else
          !in the soil
          dPond=(CumDepz2LayBottom_vr(iPondBotLev_col(N2,N1),N2,N1)-CumDepz2LayBottom_vr(NUM(N2,N1)-1,N2,N1))*10._r8
        endif  
      endif  
      VTATM = AZMAX1((PBOT_col(NY,NX)+dPond)*VLWatMicPM_vr(M,N3,N2,N1)/(RGasC*TKS_vr(N3,N2,N1)))*1.E3_r8  !mol gas/ d2
      VTGAS = sum(trcg_VOLG(idg_beg:idg_end))
!
!     PROPORTIONAL REMOVAL OF EXCESS AQUEOUS GASES
!
      IF(VTGAS.GT.VTATM)THEN
        DVTGAS  = 0.5_r8*(VTATM-VTGAS)         !<0, bubbling occurs
        FracEbu = AZMIN1(DVTGAS/VTGAS)
        DO idg=idg_beg,idg_end
          trcg_Ebu_flxM_vr(idg,N3,N2,N1)=FracEbu*trcg_VOLG(idg)*GasMassSolubility(idg)  
        ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
        DO idg=idg_beg,idg_end
          trcg_ebu_flx_vr(idg,N3,N2,N1)=trcg_ebu_flx_vr(idg,N3,N2,N1)+trcg_Ebu_flxM_vr(idg,N3,N2,N1) !<0, ebullition
        ENDDO
      ELSE
        trcg_Ebu_flxM_vr(idg_beg:idg_end,N3,N2,N1)=0.0_r8
      ENDIF
    ELSE
      iDisableEbu=itrue
      trcg_Ebu_flxM_vr(idg_beg:idg_end,N3,N2,N1)=0.0_r8
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine BubbleEffluxM
! ----------------------------------------------------------------------

  subroutine GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: trc_gasc1(idg_beg:idg_end)
  real(r8) :: trc_gasc2(idg_beg:idg_end)

  real(r8) :: CNDC1   !conductance in source grid
  real(r8) :: CNDC2   !conductance in destination grid
  real(r8) :: DFLG2   !air-filled porosity effect on gaseous diffusivity in source
  real(r8) :: DFLGL   !air-filled porosity effect on gaseous diffusivity in destination
  integer :: idg

!     GASEOUS DIFFUSIVITIES
!
!     DFLG2,DFLGL= in source,destination layer
!     POROQ=Penman Water Linear Reduction tortuosity from starts.f
!     POROS=total porosity
!     DLYR=soil layer thickness
!     D*G=gaseous diffusivity in soil
!     CND*1,CND*2=gaseous conductance in source,destination layer
!     *SGL2= gaseous diffusivity in air
!
  DFLG2=2.0_r8*AZMAX1(AirFilledSoilPoreM_vr(M,N3,N2,N1))*POROQ*AirFilledSoilPoreM_vr(M,N3,N2,N1)/POROS_vr(N3,N2,N1) &
    *AREA(N,N3,N2,N1)/DLYR_3D(N,N3,N2,N1)

  DFLGL=2.0_r8*AZMAX1(AirFilledSoilPoreM_vr(M,N6,N5,N4))*POROQ*AirFilledSoilPoreM_vr(M,N6,N5,N4)/POROS_vr(N6,N5,N4) &
    *AREA(N,N6,N5,N4)/DLYR_3D(N,N6,N5,N4)

!
!     GASOUS CONDUCTANCES
!
!
  DO idg=idg_beg,idg_NH3
    CNDC1=DFLG2*GasDifctScaledMM_vr(idg,N3,N2,N1)
    CNDC2=DFLGL*GasDifctScaledMM_vr(idg,N6,N5,N4)
    GasDifuscoefMM_3D(idg,N,N6,N5,N4)=(CNDC1*CNDC2)/(CNDC1+CNDC2)
  ENDDO
!
!     GASEOUS CONCENTRATIONS FROM AIR-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
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
    Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4) = GasDifuscoefMM_3D(idg,N,N6,N5,N4)*(trc_gasc1(idg)-trc_gasc2(idg))
  ENDDO

  end subroutine GasDiffusionMM


! ----------------------------------------------------------------------
  subroutine UpstreamGasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2SoilMM_3D)

  implicit none
  real(r8), intent(in) :: WaterFlow2SoilMM_3D(3,JD,JV,JH)
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
!     WaterFlow2SoilMM_3D=total water flux into soil micropore+macropore from watsub.f
!
  FLQW=WaterFlow2SoilMM_3D(N,N6,N5,N4)
  !flow out of dest grid
  IF(FLQW.GT.0.0_r8)THEN
    !dest grid is not saturated
    IF(VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN   
      VFLW=-AZMAX1(AMIN1(VFLWX,FLQW/VLsoiAirPM_vr(M,N6,N5,N4)))   !negative flow 
    !dest grid is aturated  
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO idg=idg_beg,idg_NH3
      RGasAdv                       = VFLW*AZMAX1(trc_gasml2_vr(idg,N6,N5,N4))
      Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4)+RGasAdv
    ENDDO
  !flow out of source grid  
  ELSE
    IF(VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=-AZMIN1(AMAX1(-VFLWX,FLQW/VLsoiAirPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO idg=idg_beg,idg_NH3
      RGasAdv                       = VFLW*AZMAX1(trc_gasml2_vr(idg,N3,N2,N1))
      Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4)+RGasAdv
    ENDDO
  ENDIF

  end subroutine UpstreamGasAdvectionMM

! ----------------------------------------------------------------------
  subroutine GasTransportMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2SoilMM_3D)

  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: WaterFlow2SoilMM_3D(3,JD,JV,JH)
  integer :: idg

!     AirFilledSoilPoreM_vr,VLsoiAirPM=air-filled porosity,volume from watsub.f

  IF(AirFilledSoilPoreM_vr(M,N3,N2,N1).GT.AirFillPore_Min .AND. AirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min &
    .AND. VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) &       !source grid is not saturated
    .AND. VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN    !dest grid is not saturated

!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE
    call GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

!     TOTAL SOIL GAS FLUX FROM CONVECTIVE FLUX
    call UpstreamGasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2SoilMM_3D)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
    DO idg=idg_beg,idg_NH3
      Gas_AdvDif_Flx_3D(idg,N,N6,N5,N4)=Gas_AdvDif_Flx_3D(idg,N,N6,N5,N4)+Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4)
    ENDDO

  ELSE
    Gas_AdvDif_FlxMM_3D(idg_beg:idg_NH3,N,N6,N5,N4)=0.0_r8
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
    IF(AirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min)THEN
      do idg=idg_beg,idg_NH3-1
        RGas_Disol_FlxMM_vr(idg,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg,N6,N5,N4))*trcg_VLWatMicP_vr(idg,N6,N5,N4) &
          -trcs_solml2_vr(idg,N6,N5,N4)*VLsoiAirPM_vr(M,N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg,N6,N5,N4)+VLsoiAirPM_vr(M,N6,N5,N4))
      enddo    

      IF(VLsoiAirPMA_vr(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGas_Disol_FlxMM_vr(idg_NH3,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg_NH3,N6,N5,N4))*trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4) &
          -trcs_solml2_vr(idg_NH3,N6,N5,N4)*VLsoiAirPMA_vr(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4)+VLsoiAirPMA_vr(N6,N5,N4))

      ELSE
        RGas_Disol_FlxMM_vr(idg_NH3,N6,N5,N4)=0.0_r8
      ENDIF

      IF(VLsoiAirPMB_vr(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGas_Disol_FlxMM_vr(idg_NH3B,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
          (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg_NH3,N6,N5,N4))*trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4) &
          -trcs_solml2_vr(idg_NH3B,N6,N5,N4)*VLsoiAirPMB_vr(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4)+VLsoiAirPMB_vr(N6,N5,N4))

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
      RGas_Disol_FlxMM_vr(idg_beg:idg_NH3,N6,N5,N4)=0.0_r8
    ENDIF
  ENDIF
  end subroutine GasDissolutionMM

! ----------------------------------------------------------------------
  subroutine ZeroTransport1MM(N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,N6

  integer :: K

  GasDifuscoefMM_3D(idg_beg:idg_end,N,N6,N5,N4) = 0._r8

  D9750: DO K=1,jcplx
    DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,N6,N5,N4)=0.0_r8
    DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,N6,N5,N4)=0.0_r8
  ENDDO D9750

  trcs_MicpTranspFlxM_3D(ids_beg:ids_end,N,N6,N5,N4)=0.0_r8

  trcs_MacpTranspFlxM_3D(ids_beg:ids_end,N,N6,N5,N4)=0.0_r8

  Gas_AdvDif_FlxMM_3D(idg_beg:idg_NH3,N,N6,N5,N4)=0.0_r8
  end subroutine ZeroTransport1MM
! ----------------------------------------------------------------------

  subroutine ZeroTransport2MM(N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,N6

  integer :: K

  GasDifuscoefMM_3D(idg_beg:idg_NH3,N,N3,N2,N1)=0.0_r8

  D9751: DO K=1,jcplx
    DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,N3,N2,N1) = 0.0_r8
    DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,N3,N2,N1) = 0.0_r8
  ENDDO D9751

  trcs_MicpTranspFlxM_3D(ids_beg:ids_end,N,N3,N2,N1) = 0.0_r8
  trcs_MacpTranspFlxM_3D(ids_beg:ids_end,N,N3,N2,N1) = 0.0_r8
  Gas_AdvDif_FlxMM_3D(idg_beg:idg_NH3,N,N3,N2,N1)         = 0.0_r8
  end subroutine ZeroTransport2MM
end module InsideTranspMod
