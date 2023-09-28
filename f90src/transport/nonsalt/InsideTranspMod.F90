module InsideTranspMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use minimathmod, only : safe_adb,AZMAX1,AZMIN1
  use GridConsts
  use EcoSIMSolverPar
  use TransfrDataMod
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

  real(r8),allocatable :: RFHOC(:)
  real(r8),allocatable :: RFHON(:)
  real(r8),allocatable :: RFHOP(:)
  real(r8),allocatable :: RFHOA(:)
  real(r8),allocatable :: COQCH1(:)
  real(r8),allocatable :: COQCH2(:)
  real(r8),allocatable :: COQNH1(:)
  real(r8),allocatable :: COQNH2(:)
  real(r8),allocatable :: COQPH1(:)
  real(r8),allocatable :: COQPH2(:)
  real(r8),allocatable :: COQAH1(:)
  real(r8),allocatable :: COQAH2(:)

  public :: ModelTracerHydroFlux
  public :: InitInsTp
  public :: DestructInsTp
  contains
!------------------------------------------------------------------------------------------
  subroutine InitInsTp()
  implicit none
  allocate(RFLOC(1:jcplx))
  allocate(RFLON(1:jcplx))
  allocate(RFLOP(1:jcplx))
  allocate(RFLOA(1:jcplx))
  allocate(RFHOC(1:jcplx))
  allocate(RFHON(1:jcplx))
  allocate(RFHOP(1:jcplx))
  allocate(RFHOA(1:jcplx))
  allocate(COQC1(1:jcplx))
  allocate(COQC2(1:jcplx))
  allocate(COQN1(1:jcplx))
  allocate(COQN2(1:jcplx))
  allocate(COQP1(1:jcplx))
  allocate(COQP2(1:jcplx))
  allocate(COQA1(1:jcplx))
  allocate(COQA2(1:jcplx))
  allocate(COQCH1(1:jcplx))
  allocate(COQCH2(1:jcplx))
  allocate(COQNH1(1:jcplx))
  allocate(COQNH2(1:jcplx))
  allocate(COQPH1(1:jcplx))
  allocate(COQPH2(1:jcplx))
  allocate(COQAH1(1:jcplx))
  allocate(COQAH2(1:jcplx))
  allocate(DFVOC(1:jcplx))
  allocate(DFVON(1:jcplx))
  allocate(DFVOP(1:jcplx))
  allocate(DFVOA(1:jcplx))
  allocate(DFHOC(1:jcplx))
  allocate(DFHON(1:jcplx))
  allocate(DFHOP(1:jcplx))
  allocate(DFHOA(1:jcplx))

  end subroutine InitInsTp
!------------------------------------------------------------------------------------------
  subroutine DestructInsTp
  implicit none
  call destroy(RFLOC)
  call destroy(RFLON)
  call destroy(RFLOP)
  call destroy(RFLOA)
  call destroy(RFHOC)
  call destroy(RFHON)
  call destroy(RFHOP)
  call destroy(RFHOA)
  call destroy(COQC1)
  call destroy(COQC2)
  call destroy(COQN1)
  call destroy(COQN2)
  call destroy(COQP1)
  call destroy(COQP2)
  call destroy(COQA1)
  call destroy(COQA2)
  call destroy(COQCH1)
  call destroy(COQCH2)
  call destroy(COQNH1)
  call destroy(COQNH2)
  call destroy(COQPH1)
  call destroy(COQPH2)
  call destroy(COQAH1)
  call destroy(COQAH2)
  call destroy(DFVOC)
  call destroy(DFVON)
  call destroy(DFVOP)
  call destroy(DFVOA)
  call destroy(DFHOC)
  call destroy(DFHON)
  call destroy(DFHOP)
  call destroy(DFHOA)

  end subroutine DestructInsTp
!------------------------------------------------------------------------------------------

  subroutine ModelTracerHydroFlux(M,MX, NHW, NHE, NVN, NVS,WaterFlow2Soil)
  implicit none

  integer, intent(in) :: M,MX, NHW, NHE, NVN, NVS
  real(r8), intent(inout) :: WaterFlow2Soil(3,JD,JV,JH)
  integer :: NY,NX
  real(r8) :: FLWRM1
  real(r8) :: trcg_RFLS0(idg_beg:idg_end-1)
  real(r8) :: trcn_RFLW0(ids_nut_beg:ids_nuts_end)
  real(r8) :: RDXS_gas(idg_beg:idg_end)

  RDXS_gas(idg_beg:idg_end)=0._r8
  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      call ResetFluxAccumulators(M,NY,NX,MX)

      IF(M.NE.MX)THEN
!
!     This IF statement is the next ~1700 lines so I'm leaving it in here
!
!     SOLUTE FLUXES FROM MELTING SNOWPACK TO
!     RESIDUE AND SOIL SURFACE FROM SNOWMELT IN 'WATSUB' AND
!     CONCENTRATIONS IN SNOWPACK
!
        call SoluteFluxSnowpackDisch(M,NY,NX,trcg_RFLS0,trcn_RFLW0)
!
!     SOLUTE FLUXES AT SOIL SURFACE FROM SURFACE WATER
!     CONTENTS, WATER FLUXES 'WaterFlow2Soil' AND ATMOSPHERE BOUNDARY
!     LAYER RESISTANCES 'PARGM' FROM 'WATSUB'
!
        call SoluteFluxSurface(M,NY,NX,NHE,NHW,NVS,NVN,&
          WaterFlow2Soil,trcg_RFLS0,trcn_RFLW0,RDXS_gas)
!
      ENDIF
!
!     VOLATILIZATION-DISSOLUTION OF GASES IN RESIDUE AND SOIL SURFACE
!     LAYERS FROM GASEOUS CONCENTRATIONS VS. THEIR AQUEOUS
!     EQUIVALENTS DEPENDING ON SOLUBILITY FROM 'HOUR1'
!     AND TRANSFER COEFFICIENT 'DFGS' FROM 'WATSUB'
!
      call LitterGasVolatilDissol(M,NY,NX)
!
!     SURFACE GAS EXCHANGE FROM GAS DIFFUSIVITY THROUGH
!     SOIL SURFACE LAYER AND THROUGH ATMOSPHERE BOUNDARY
!     LAYER
!
      call SurfSoilFluxGasDifAdv(M,NY,NX,WaterFlow2Soil,RDXS_gas)
!
!     SOIL SURFACE WATER-AIR GAS EXCHANGE
!
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
!
      call TracerExchInBetweenCells(M,MX,NY,NX,NHE,NVS,WaterFlow2Soil)

    ENDDO
  ENDDO
  end subroutine ModelTracerHydroFlux
!------------------------------------------------------------------------------------------

  subroutine ResetFluxAccumulators(M,NY,NX,MX)
  implicit none

  integer, intent(in) :: M, NY, NX,MX

  integer :: K,L,NTG,NTS
  real(r8) :: PARGM

  IF(M.NE.MX)THEN
!
!     GASEOUS BOUNDARY LAYER CONDUCTANCES
!
!     PARG=boundary layer conductance above soil surface from watsub.f
!
    PARGM=PARG(M,NY,NX)*XNPT
    PARG_cef(idg_CO2,NY,NX)=PARGM*0.74_r8
    PARG_cef(idg_CH4,NY,NX)=PARGM*1.04_r8
    PARG_cef(idg_O2,NY,NX)=PARGM*0.83_r8
    PARG_cef(idg_N2,NY,NX)=PARGM*0.86_r8
    PARG_cef(idg_N2O,NY,NX)=PARGM*0.74_r8
    PARG_cef(idg_NH3,NY,NX)=PARGM*1.02_r8
    PARG_cef(idg_H2,NY,NX)=PARGM*2.08_r8
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
      TQROC(K,NY,NX)=0.0_r8
      TQRON(K,NY,NX)=0.0_r8
      TQROP(K,NY,NX)=0.0_r8
      TQROA(K,NY,NX)=0.0_r8
      OQC2(K,0,NY,NX)=OQC2(K,0,NY,NX)-ROCSK2(K,0,NY,NX)
      OQN2(K,0,NY,NX)=OQN2(K,0,NY,NX)-RONSK2(K,0,NY,NX)
      OQP2(K,0,NY,NX)=OQP2(K,0,NY,NX)-ROPSK2(K,0,NY,NX)
      OQA2(K,0,NY,NX)=OQA2(K,0,NY,NX)-ROASK2(K,0,NY,NX)
    ENDDO

    trcg_TQR(idg_beg:idg_end-1,NY,NX)=0.0_r8
    trcn_TQR(ids_nut_beg:ids_nuts_end,NY,NX)=0.0_r8
    trcg_TQ(idg_beg:idg_end-1,NY,NX)=0.0_r8
    trcn_TQ(ids_nut_beg:ids_nuts_end,NY,NX)=0.0_r8

    TQSCOS(NY,NX)=0.0_r8

!   because NH3 is gas-aqua dual phase
    trc_solml2(idg_NH3,0,NY,NX)=trc_solml2(idg_NH3,0,NY,NX)-RBGCSinkS(idg_NH3,0,NY,NX)
!   exclude nutrients in band
    DO NTS=ids_nut_beg,ids_nuts_end
      trc_solml2(NTS,0,NY,NX)=trc_solml2(NTS,0,NY,NX)-RBGCSinkS(NTS,0,NY,NX)
    ENDDO
    RBGCSinkG(idg_O2,0,NY,NX)=ROXSK(M,0,NY,NX)*XNPT
  ENDIF

  DO NTG=idg_beg,idg_end-1
    if(NTG/=idg_NH3)then
      trc_solml2(NTG,0,NY,NX)=trc_solml2(NTG,0,NY,NX)-RBGCSinkG(NTG,0,NY,NX)
    ELSE
      trc_gasml2(idg_NH3,0,NY,NX)=trc_gasml2(idg_NH3,0,NY,NX)-RBGCSinkG(idg_NH3,0,NY,NX)
    ENDIF
  ENDDO
!
!     INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
!
  IF(M.NE.MX)THEN
    DO  L=1,JS
      trcg_TBLS(idg_beg:idg_end-1,L,NY,NX)=0._r8
      trcn_TBLS(ids_nut_beg:ids_nuts_end,L,NY,NX)=0._r8
    ENDDO
  ENDIF
!
!     INITIALIZE SOIL SOLUTE NET FLUX ACCUMULATORS
!
  DO L=NU(NY,NX),NL(NY,NX)
    IF(M.NE.MX)THEN
      DO  K=1,jcplx
        TOCFLS(K,L,NY,NX)=0.0_r8
        TONFLS(K,L,NY,NX)=0.0_r8
        TOPFLS(K,L,NY,NX)=0.0_r8
        TOAFLS(K,L,NY,NX)=0.0_r8
        TOCFHS(K,L,NY,NX)=0.0_r8
        TONFHS(K,L,NY,NX)=0.0_r8
        TOPFHS(K,L,NY,NX)=0.0_r8
        TOAFHS(K,L,NY,NX)=0.0_r8
        OQC2(K,L,NY,NX)=OQC2(K,L,NY,NX)-ROCSK2(K,L,NY,NX)
        OQN2(K,L,NY,NX)=OQN2(K,L,NY,NX)-RONSK2(K,L,NY,NX)
        OQP2(K,L,NY,NX)=OQP2(K,L,NY,NX)-ROPSK2(K,L,NY,NX)
        OQA2(K,L,NY,NX)=OQA2(K,L,NY,NX)-ROASK2(K,L,NY,NX)
      ENDDO
      R3PorTSolFlx(ids_beg:ids_end,L,NY,NX)=0.0_r8
      R3PorTSoHFlx(ids_beg:ids_end,L,NY,NX)=0._r8
!
!     ADD SOLUTE SINKS
!
!     R*SK2=total flux from nitro.f, uptake.f, solute.f
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band micropores
!     ZNH4B,ZNH3B,ZNO3B,ZNO2B,H1POB,H2POB=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band micropores
!
!     include NH3 and band nutrients
      DO NTS=ids_nuts_beg,ids_nuts_end
        trc_solml2(NTS,L,NY,NX)=trc_solml2(NTS,L,NY,NX)-RBGCSinkS(NTS,L,NY,NX)
      ENDDO

      RBGCSinkG(idg_O2,L,NY,NX)=ROXSK(M,L,NY,NX)*XNPT
    ENDIF
!
!     SOIL GAS FLUX ACCUMULATORS
!
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in micropores
!     ZN3G=gaseous NH3
!
    RTGasADFlx(idg_beg:idg_end-1,L,NY,NX)=0._r8
    DO NTG=idg_beg,idg_end-1
      if(NTG/=idg_NH3)then
        trc_solml2(NTG,L,NY,NX)=trc_solml2(NTG,L,NY,NX)-RBGCSinkG(NTG,L,NY,NX)
      else
        trc_gasml2(idg_NH3,L,NY,NX)=trc_gasml2(idg_NH3,L,NY,NX)-RBGCSinkG(NTG,L,NY,NX)
      endif
    ENDDO

  ENDDO
  end subroutine ResetFluxAccumulators


!------------------------------------------------------------------------------------------

  subroutine TracerExchInBetweenCells(M,MX,NY,NX,NHE,NVS,WaterFlow2Soil)
!
! DESCRIPTION
! exchanges tracers within (gaseous vs aqueous phase) and between
! grid cells.
  implicit none

  integer, intent(in) :: M,MX, NY, NX, NHE, NVS
  real(r8), intent(inout) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8) :: VFLW
  real(r8) :: VOLH2A,VOLH2B
  real(r8) :: VLWatMacPS,VOLWT

  integer :: IFLGB,N,L,K,LL
  integer :: N1,N2,N3,N4,N5,N6

! begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  IFLGB=0
  D125: DO L=1,NL(NY,NX)
    N1=NX;N2=NY;N3=L
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
    D120: DO N=FlowDirIndicator(N2,N1),3

      IF(N.EQ.1)THEN
        !WEST-EAST
        IF(NX.EQ.NHE)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N6=L
        ENDIF
      ELSEIF(N.EQ.2)THEN
        !NORTH-SOUTH
        IF(NY.EQ.NVS)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N6=L
        ENDIF
      ELSEIF(N.EQ.3)THEN
        !VERTICAL
        IF(L.EQ.NL(NY,NX))THEN
          cycle
        ELSE
          N4=NX
          N5=NY
          N6=L+1
        ENDIF
      ENDIF

      DO LL=N6,NL(NY,NX)
        IF(VLSoilPoreMicP(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
          N6=LL
          exit
        ENDIF
      ENDDO
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
!     WATER CONTENTS AND WATER FLUXES 'WaterFlow2Soil' FROM 'WATSUB'
!
!     VLSoilPoreMicP,VLSoilMicP=soil volume excluding rock, macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     VLWatMicPM,VLWatMacPM=micropore,macropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!     ReductVLsoiAirPM=change in air volume
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
      IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
        IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4) &
          .AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
          IF(M.NE.MX)THEN
            VLWatMicPMA(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*trcs_VLN(ids_NH4,N6,N5,N4)
            VLWatMicPMB(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*trcs_VLN(ids_NH4B,N6,N5,N4)
            VLWatMicPXA(N6,N5,N4)=natomw*VLWatMicPMA(N6,N5,N4)
            VLWatMicPXB(N6,N5,N4)=natomw*VLWatMicPMB(N6,N5,N4)

            VLsoiAirPMA(N6,N5,N4)=VLsoiAirPM(M,N6,N5,N4)*trcs_VLN(ids_NH4,N6,N5,N4)
            VLsoiAirPMB(N6,N5,N4)=VLsoiAirPM(M,N6,N5,N4)*trcs_VLN(ids_NH4B,N6,N5,N4)
            CumReductVLsoiAirPM(N6,N5,N4)=ReductVLsoiAirPM(M,N6,N5,N4)*XNPT
!
!     GASEOUS SOLUBILITIES
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     VLWatMicP*=equivalent aqueous volume for gas
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     S*L=solubility of gas in water from hour1.f
!     WaterFlow2Soil=total water flux into soil micropore+macropore from watsub.f
!
            WaterFlow2Soil(N,N6,N5,N4)=(WaterFlow2MicPM(M,N,N6,N5,N4)+WaterFlow2MacPM(M,N,N6,N5,N4))*XNPT
!
            call SoluteAdvDifTransport(M,N,N1,N2,N3,N4,N5,N6)

!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
            IF(N.EQ.3)THEN
              call MicMacPoresSoluteExchange(M,N,N1,N2,N3,N4,N5,N6)
            ENDIF
          ENDIF
!
!     GASEOUS TRANSPORT FROM GASEOUS DIFFUSIVITY AND CONCENTRATION
!     DIFFERENCES BETWEEN ADJACENT GRID CELLS
!
          call GaseousTransport(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)

        ELSEIF(N.NE.3)THEN
          call ZeroTransport1(N,N1,N2,N3,N4,N5,N6)

        ENDIF
      ELSE
        call ZeroTransport2(N,N1,N2,N3,N4,N5,N6)
      ENDIF
    ENDDO D120
!
!     CHECK FOR BUBBLING IF THE SUM OF ALL GASEOUS EQUIVALENT
!     PARTIAL CONCENTRATIONS EXCEEDS ATMOSPHERIC PRESSURE
    call BubbleEfflux(M,N1,N2,N3,NY,NX,MX,IFLGB)

  ENDDO D125
  end subroutine TracerExchInBetweenCells

! ----------------------------------------------------------------------

   subroutine SoluteAdvTranspMicropore(M,N,N1,N2,N3,N4,N5,N6)
   implicit none
   integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
   real(r8) :: trc_RFL(ids_beg:ids_end)

   integer :: K,NTS
   real(r8) :: VFLW
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     WaterFlow2MicPM=water flux through soil micropore from watsub.f
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     RFL*S=solute diffusive flux through micropore
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *S2,*B2=micropore solute content in non-band,band
!
  IF(WaterFlow2MicPM(M,N,N6,N5,N4).GT.0.0_r8)THEN
    IF(VLWatMicPM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MicPM(M,N,N6,N5,N4)/VLWatMicPM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,N3,N2,N1))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,N3,N2,N1))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,N3,N2,N1))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,N3,N2,N1))
    ENDDO

    DO NTS=ids_beg,ids_end
    trc_RFL(NTS)=VFLW*AZMAX1(trc_solml2(NTS,N3,N2,N1))
    ENDDO
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS TO CURRENT FROM
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
  ELSE
    IF(VLWatMicPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MicPM(M,N,N6,N5,N4)/VLWatMicPM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    D9815: DO K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,N6,N5,N4))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,N6,N5,N4))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,N6,N5,N4))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,N6,N5,N4))
    ENDDO D9815
    DO NTS=ids_beg,ids_end
      trc_RFL(NTS)=VFLW*AZMAX1(trc_solml2(NTS,N6,N5,N4))
    ENDDO
  ENDIF

  DO NTS=ids_beg,ids_end
    R3PoreSolFlx(NTS,N,N6,N5,N4)=trc_RFL(NTS)
  ENDDO

  end subroutine SoluteAdvTranspMicropore

! ----------------------------------------------------------------------
  subroutine SoluteDifTranspMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: THETW1(JZ,JY,JX)
  real(r8) :: VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: VLWatMicP2A,VLWatMicP2B,VLWatMicP3A,VLWatMicP3B,VLWatMicP4A,VLWatMicP4B
  real(r8) :: trcsolc1(ids_beg:ids_end)
  real(r8) :: trcsolc2(ids_beg:ids_end)

  real(r8) :: SDifc(ids_beg:ids_end),SDifFlx(ids_beg:ids_end)
  real(r8) :: DISPN,DIFOC,DIFON,DIFOP,DIFOA
  real(r8) :: DLYR1,DLYR2,TORTL
  integer  :: K,NTS,NTG

  IF(THETW1(N3,N2,N1).GT.THETY(N3,N2,N1).AND.THETW1(N6,N5,N4).GT.THETY(N6,N5,N4) &
    .AND.VLWatMicPM(M,N3,N2,N1).GT.ZEROS2(N2,N1).AND.VLWatMicPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN

      VLWatMicP2A=VLWatMicPM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4,N3,N2,N1)
      VLWatMicP2B=VLWatMicPM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      VLWatMicP3A=VLWatMicPM(M,N3,N2,N1)*trcs_VLN(ids_NO3,N3,N2,N1)
      VLWatMicP3B=VLWatMicPM(M,N3,N2,N1)*trcs_VLN(ids_NO3B,N3,N2,N1)
      VLWatMicP4A=VLWatMicPM(M,N3,N2,N1)*trcs_VLN(ids_NH4,N3,N2,N1)
      VLWatMicP4B=VLWatMicPM(M,N3,N2,N1)*trcs_VLN(ids_NH4B,N3,N2,N1)

      VLWatMicPPA=VLWatMicPM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4,N6,N5,N4)
      VLWatMicPPB=VLWatMicPM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      VLWatMicPOA=VLWatMicPM(M,N6,N5,N4)*trcs_VLN(ids_NO3,N6,N5,N4)
      VLWatMicPOB=VLWatMicPM(M,N6,N5,N4)*trcs_VLN(ids_NO3B,N6,N5,N4)

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
      COQC1(K)=AZMAX1(OQC2(K,N3,N2,N1)/VLWatMicPM(M,N3,N2,N1))
      COQN1(K)=AZMAX1(OQN2(K,N3,N2,N1)/VLWatMicPM(M,N3,N2,N1))
      COQP1(K)=AZMAX1(OQP2(K,N3,N2,N1)/VLWatMicPM(M,N3,N2,N1))
      COQA1(K)=AZMAX1(OQA2(K,N3,N2,N1)/VLWatMicPM(M,N3,N2,N1))

      COQC2(K)=AZMAX1(OQC2(K,N6,N5,N4)/VLWatMicPM(M,N6,N5,N4))
      COQN2(K)=AZMAX1(OQN2(K,N6,N5,N4)/VLWatMicPM(M,N6,N5,N4))
      COQP2(K)=AZMAX1(OQP2(K,N6,N5,N4)/VLWatMicPM(M,N6,N5,N4))
      COQA2(K)=AZMAX1(OQA2(K,N6,N5,N4)/VLWatMicPM(M,N6,N5,N4))
    ENDDO D9810

    DO NTG=idg_beg,idg_end-2
      trcsolc1(NTG)=AZMAX1(trc_solml2(NTG,N3,N2,N1)/VLWatMicPM(M,N3,N2,N1))
    ENDDO

    IF(VLWatMicP4A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4)=AZMAX1(trc_solml2(ids_NH4,N3,N2,N1)/VLWatMicP4A)
      trcsolc1(idg_NH3)=AZMAX1(trc_solml2(idg_NH3,N3,N2,N1)/VLWatMicP4A)
    ELSE
      trcsolc1(ids_NH4)=0.0_r8
      trcsolc1(idg_NH3)=0.0_r8
    ENDIF
    IF(VLWatMicP3A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3)=AZMAX1(trc_solml2(ids_NO3,N3,N2,N1)/VLWatMicP3A)
      trcsolc1(ids_NO2)=AZMAX1(trc_solml2(ids_NO2,N3,N2,N1)/VLWatMicP3A)
    ELSE
      trcsolc1(ids_NO3)=0.0_r8
      trcsolc1(ids_NO2)=0.0_r8
    ENDIF
    IF(VLWatMicP2A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4)=AZMAX1(trc_solml2(ids_H1PO4,N3,N2,N1)/VLWatMicP2A)
      trcsolc1(ids_H2PO4)=AZMAX1(trc_solml2(ids_H2PO4,N3,N2,N1)/VLWatMicP2A)
    ELSE
      trcsolc1(ids_H1PO4)=0.0_r8
      trcsolc1(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VLWatMicP4B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4B)=AZMAX1(trc_solml2(ids_NH4B,N3,N2,N1)/VLWatMicP4B)
      trcsolc1(idg_NH3B)=AZMAX1(trc_solml2(idg_NH3B,N3,N2,N1)/VLWatMicP4B)
    ELSE
      trcsolc1(ids_NH4B)=0.0_r8
      trcsolc1(idg_NH3B)=0.0_r8
    ENDIF
    IF(VLWatMicP3B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3B)=AZMAX1(trc_solml2(ids_NO3B,N3,N2,N1)/VLWatMicP3B)
      trcsolc1(ids_NO2B)=AZMAX1(trc_solml2(ids_NO2B,N3,N2,N1)/VLWatMicP3B)
    ELSE
      trcsolc1(ids_NO3B)=trcsolc1(ids_NO3)
      trcsolc1(ids_NO2B)=trcsolc1(ids_NO2)
    ENDIF
    IF(VLWatMicP2B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4B)=AZMAX1(trc_solml2(ids_H1PO4B,N3,N2,N1)/VLWatMicP2B)
      trcsolc1(ids_H2PO4B)=AZMAX1(trc_solml2(ids_H2PO4B,N3,N2,N1)/VLWatMicP2B)
    ELSE
      trcsolc1(ids_H1PO4B)=trcsolc1(ids_H1PO4)
      trcsolc1(ids_H2PO4B)=trcsolc1(ids_H2PO4)
    ENDIF

    DO NTG=idg_beg,idg_end-2
      trcsolc2(NTG)=AZMAX1(trc_solml2(NTG,N6,N5,N4)/VLWatMicPM(M,N6,N5,N4))
    ENDDO


    IF(VLWatMicPMA(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      trcsolc2(idg_NH3)=AZMAX1(trc_solml2(idg_NH3,N6,N5,N4)/VLWatMicPMA(N6,N5,N4))
      trcsolc2(ids_NH4)=AZMAX1(trc_solml2(ids_NH4,N6,N5,N4)/VLWatMicPMA(N6,N5,N4))
    ELSE
      trcsolc2(idg_NH3)=0.0_r8
      trcsolc2(ids_NH4)=0.0_r8
    ENDIF
    IF(VLWatMicPOA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3)=AZMAX1(trc_solml2(ids_NO3,N6,N5,N4)/VLWatMicPOA)
      trcsolc2(ids_NO2)=AZMAX1(trc_solml2(ids_NO2,N6,N5,N4)/VLWatMicPOA)
    ELSE
      trcsolc2(ids_NO3)=0.0_r8
      trcsolc2(ids_NO2)=0.0_r8
    ENDIF
    IF(VLWatMicPPA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4)=AZMAX1(trc_solml2(ids_H1PO4,N6,N5,N4)/VLWatMicPPA)
      trcsolc2(ids_H2PO4)=AZMAX1(trc_solml2(ids_H2PO4,N6,N5,N4)/VLWatMicPPA)
    ELSE
      trcsolc2(ids_H1PO4)=0.0_r8
      trcsolc2(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VLWatMicPMB(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      trcsolc2(idg_NH3B)=AZMAX1(trc_solml2(idg_NH3B,N6,N5,N4)/VLWatMicPMB(N6,N5,N4))
      trcsolc2(ids_NH4B)=AZMAX1(trc_solml2(ids_NH4B,N6,N5,N4)/VLWatMicPMB(N6,N5,N4))
    ELSE
      trcsolc2(idg_NH3B)=trcsolc2(idg_NH3)
      trcsolc2(ids_NH4B)=trcsolc2(ids_NH4)
    ENDIF
    IF(VLWatMicPOB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3B)=AZMAX1(trc_solml2(ids_NO3B,N6,N5,N4)/VLWatMicPOB)
      trcsolc2(ids_NO2B)=AZMAX1(trc_solml2(ids_NO2B,N6,N5,N4)/VLWatMicPOB)
    ELSE
      trcsolc2(ids_NO3B)=trcsolc2(ids_NO3)
      trcsolc2(ids_NO2B)=trcsolc2(ids_NO2)
    ENDIF
    IF(VLWatMicPPB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4B)=AZMAX1(trc_solml2(ids_H1PO4B,N6,N5,N4)/VLWatMicPPB)
      trcsolc2(ids_H2PO4B)=AZMAX1(trc_solml2(ids_H2PO4B,N6,N5,N4)/VLWatMicPPB)
    ELSE
      trcsolc2(ids_H1PO4B)=trcsolc2(ids_H1PO4)
      trcsolc2(ids_H2PO4B)=trcsolc2(ids_H2PO4)
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!
!     DLYR=soil layer thickness
!     TORT=micropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     WaterFlow2MicPM=water flux through soil micropore from watsub.f
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
    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TortMicPM(M,N3,N2,N1)*DLYR1+TortMicPM(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MicPM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DIFOC=(OCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFON=(ONSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOP=(OPSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOA=(OASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DO NTS=ids_beg,ids_end
      SDifc(NTS)=(SolDifcc(NTS,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
    D9805: DO K=1,jcplx
      DFVOC(K)=DIFOC*(COQC1(K)-COQC2(K))
      DFVON(K)=DIFON*(COQN1(K)-COQN2(K))
      DFVOP(K)=DIFOP*(COQP1(K)-COQP2(K))
      DFVOA(K)=DIFOA*(COQA1(K)-COQA2(K))
    ENDDO D9805

    DO NTG=idg_beg,idg_end-2
      SDifFlx(NTG)=SDifc(NTG)*(trcsolc1(NTG)-trcsolc2(NTG))
    ENDDO

    DO NTS=ids_nuts_beg,ids_nuts_end
      SDifFlx(NTS)=SDifc(NTS)*(trcsolc1(NTS)-trcsolc2(NTS)) &
        *AMIN1(trcs_VLN(NTS,N3,N2,N1),trcs_VLN(NTS,N6,N5,N4))
    ENDDO
  ELSE
    D9905: DO K=1,jcplx
      DFVOC(K)=0.0_r8
      DFVON(K)=0.0_r8
      DFVOP(K)=0.0_r8
      DFVOA(K)=0.0_r8
    ENDDO D9905
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF

  DO NTS=ids_beg,ids_end
    R3PoreSolFlx(NTS,N,N6,N5,N4)=R3PoreSolFlx(NTS,N,N6,N5,N4)+SDifFlx(NTS)
  ENDDO

  end subroutine SoluteDifTranspMicropore

! ----------------------------------------------------------------------
  subroutine SoluteAdvTranspMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: trcs_RFH(ids_beg:ids_end)
  integer  :: K,NTG,NTS
  real(r8) :: VFLW
!     WaterFlow2MacPM=water flux through soil macropore from watsub.f
!

  IF(WaterFlow2MacPM(M,N,N6,N5,N4).GT.0.0_r8)THEN
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
    IF(VLWatMacPM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MacPM(M,N,N6,N5,N4)/VLWatMacPM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
!
!     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE
!
    IF(N.EQ.3.AND.VLMacP(N6,N5,N4).GT.VLWatMacPM(M,N6,N5,N4))THEN
      D9800: DO K=1,jcplx
        RFHOC(K)=VFLW*AZMAX1((OQCH2(K,N3,N2,N1)-AZMIN1(ROCFXS(K,NU(N2,N1),N2,N1))))
        RFHON(K)=VFLW*AZMAX1((OQNH2(K,N3,N2,N1)-AZMIN1(RONFXS(K,NU(N2,N1),N2,N1))))
        RFHOP(K)=VFLW*AZMAX1((OQPH2(K,N3,N2,N1)-AZMIN1(ROPFXS(K,NU(N2,N1),N2,N1))))
        RFHOA(K)=VFLW*AZMAX1((OQAH2(K,N3,N2,N1)-AZMIN1(ROAFXS(K,NU(N2,N1),N2,N1))))
      ENDDO D9800

      DO NTG=idg_beg,idg_end-2
        trcs_RFH(NTG)=VFLW*AZMAX1((trc_soHml2(NTG,N3,N2,N1) &
          -AZMIN1(RporeSoXFlx(NTG,NU(N2,N1),N2,N1))))
      ENDDO

      DO NTS=ids_nuts_beg,ids_nuts_end
        trcs_RFH(NTS)=VFLW*AZMAX1((trc_soHml2(NTS,N3,N2,N1) &
          -AZMIN1(RporeSoXFlx(NTS,NU(N2,N1),N2,N1)*trcs_VLN(NTS,N3,N2,N1)))) &
          *trcs_VLN(NTS,N6,N5,N4)
      ENDDO
!
!     OTHERWISE
!
    ELSE
      D9850: DO K=1,jcplx
        RFHOC(K)=VFLW*AZMAX1(OQCH2(K,N3,N2,N1))
        RFHON(K)=VFLW*AZMAX1(OQNH2(K,N3,N2,N1))
        RFHOP(K)=VFLW*AZMAX1(OQPH2(K,N3,N2,N1))
        RFHOA(K)=VFLW*AZMAX1(OQAH2(K,N3,N2,N1))
      ENDDO D9850
!exclude NH3 and NH3B
      DO NTG=idg_beg,idg_end-2
        trcs_RFH(NTG)=VFLW*AZMAX1(trc_soHml2(NTG,N3,N2,N1))
      ENDDO

      DO NTS=ids_nuts_beg,ids_nuts_end
        trcs_RFH(NTS)=VFLW*AZMAX1(trc_soHml2(NTS,N3,N2,N1))*trcs_VLN(NTS,N6,N5,N4)
      ENDDO
    ENDIF
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN ADJACENT
!     GRID CELL
!
  ELSEIF(WaterFlow2MacPM(M,N,N6,N5,N4).LT.0.0_r8)THEN
    IF(VLWatMacPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MacPM(M,N,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    D9665: DO K=1,jcplx
      RFHOC(K)=VFLW*AZMAX1(OQCH2(K,N6,N5,N4))
      RFHON(K)=VFLW*AZMAX1(OQNH2(K,N6,N5,N4))
      RFHOP(K)=VFLW*AZMAX1(OQPH2(K,N6,N5,N4))
      RFHOA(K)=VFLW*AZMAX1(OQAH2(K,N6,N5,N4))
    ENDDO D9665

    DO NTG=idg_beg,idg_end-2
      trcs_RFH(NTG)=VFLW*AZMAX1(trc_soHml2(NTG,N6,N5,N4))
    ENDDO

    DO NTS=ids_nuts_beg,ids_nuts_end
      trcs_RFH(NTS)=VFLW*AZMAX1(trc_soHml2(NTS,N6,N5,N4))*trcs_VLN(NTS,N6,N5,N4)
    ENDDO
  ELSE
!
!     NO MACROPORE FLUX
!
    D9795: DO K=1,jcplx
      RFHOC(K)=0.0_r8
      RFHON(K)=0.0_r8
      RFHOP(K)=0.0_r8
      RFHOA(K)=0.0_r8
    ENDDO D9795
    trcs_RFH(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO NTS=ids_beg,ids_end
    R3PoreSoHFlx(NTS,N,N6,N5,N4)=trcs_RFH(NTS)
  ENDDO

  end subroutine SoluteAdvTranspMacropore

! ----------------------------------------------------------------------
  subroutine SoluteDifTranspMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: VOLH2A,VOLH2B,VOLH3A,VOLH3B,VOLH4A,VOLH4B
  real(r8) :: VOLHMA,VOLHMB,VOLHOA,VOLHOB,VOLHPA,VOLHPB
  real(r8) :: trcs_coH1(ids_beg:ids_end)
  real(r8) :: trcs_coH2(ids_beg:ids_end)
  real(r8) :: SDifc(ids_beg:ids_end),TORTL
  real(r8) :: DIFOC,DIFON,DIFOP,DIFOA
  real(r8) :: DISPN,DLYR1,DLYR2
  real(r8) :: SDifHFlx(ids_beg:ids_end)
  integer  :: K,NTS,NTG

!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MACROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     THETY=hygroscopic water content
!     VOLAH=total macropore volume
!
  IF(VLWatMacPM(M,N3,N2,N1).GT.THETY(N3,N2,N1)*VLMacP(N3,N2,N1) &
    .AND.VLWatMacPM(M,N6,N5,N4).GT.THETY(N6,N5,N4)*VLMacP(N6,N5,N4))THEN
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
    VOLH4A=VLWatMacPM(M,N3,N2,N1)*trcs_VLN(ids_NH4,N3,N2,N1)
    VOLH4B=VLWatMacPM(M,N3,N2,N1)*trcs_VLN(ids_NH4B,N3,N2,N1)
    VOLH3A=VLWatMacPM(M,N3,N2,N1)*trcs_VLN(ids_NO3,N3,N2,N1)
    VOLH3B=VLWatMacPM(M,N3,N2,N1)*trcs_VLN(ids_NO3B,N3,N2,N1)
    VOLH2A=VLWatMacPM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4,N3,N2,N1)
    VOLH2B=VLWatMacPM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4B,N3,N2,N1)

    VOLHOA=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_NO3,N6,N5,N4)
    VOLHOB=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_NO3B,N6,N5,N4)
    VOLHPA=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    VOLHPB=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4B,N6,N5,N4)

    D9790: DO K=1,jcplx
      COQCH1(K)=AZMAX1(OQCH2(K,N3,N2,N1)/VLWatMacPM(M,N3,N2,N1))
      COQNH1(K)=AZMAX1(OQNH2(K,N3,N2,N1)/VLWatMacPM(M,N3,N2,N1))
      COQPH1(K)=AZMAX1(OQPH2(K,N3,N2,N1)/VLWatMacPM(M,N3,N2,N1))
      COQAH1(K)=AZMAX1(OQAH2(K,N3,N2,N1)/VLWatMacPM(M,N3,N2,N1))
      COQCH2(K)=AZMAX1(OQCH2(K,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4))
      COQNH2(K)=AZMAX1(OQNH2(K,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4))
      COQPH2(K)=AZMAX1(OQPH2(K,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4))
      COQAH2(K)=AZMAX1(OQAH2(K,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4))
    ENDDO D9790
!exclude NH3 and NH3B
    DO NTG=idg_beg,idg_end-2
      trcs_coH1(NTG)=AZMAX1(trc_soHml2(NTG,N3,N2,N1)/VLWatMacPM(M,N3,N2,N1))
    ENDDO

    IF(VOLH4A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4)=AZMAX1(trc_soHml2(ids_NH4,N3,N2,N1)/VOLH4A)
      trcs_coH1(idg_NH3)=AZMAX1(trc_soHml2(idg_NH3,N3,N2,N1)/VOLH4A)
    ELSE
      trcs_coH1(ids_NH4)=0.0_r8
      trcs_coH1(idg_NH3)=0.0_r8
    ENDIF
    IF(VOLH3A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3)=AZMAX1(trc_soHml2(ids_NO3,N3,N2,N1)/VOLH3A)
      trcs_coH1(ids_NO2)=AZMAX1(trc_soHml2(ids_NO2,N3,N2,N1)/VOLH3A)
    ELSE
      trcs_coH1(ids_NO3)=0.0_r8
      trcs_coH1(ids_NO2)=0.0_r8
    ENDIF
    IF(VOLH2A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4)=AZMAX1(trc_soHml2(ids_H1PO4,N3,N2,N1)/VOLH2A)
      trcs_coH1(ids_H2PO4)=AZMAX1(trc_soHml2(ids_H2PO4,N3,N2,N1)/VOLH2A)
    ELSE
      trcs_coH1(ids_H1PO4)=0.0_r8
      trcs_coH1(ids_H2PO4)=0.0_r8
    ENDIF
    IF(VOLH4B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4B)=AZMAX1(trc_soHml2(ids_NH4B,N3,N2,N1)/VOLH4B)
      trcs_coH1(idg_NH3B)=AZMAX1(trc_soHml2(idg_NH3B,N3,N2,N1)/VOLH4B)
    ELSE
      trcs_coH1(ids_NH4B)=trcs_coH1(ids_NH4)
      trcs_coH1(idg_NH3B)=trcs_coH1(idg_NH3)
    ENDIF
    IF(VOLH3B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3B)=AZMAX1(trc_soHml2(ids_NO3B,N3,N2,N1)/VOLH3B)
      trcs_coH1(ids_NO2B)=AZMAX1(trc_soHml2(ids_NO2B,N3,N2,N1)/VOLH3B)
    ELSE
      trcs_coH1(ids_NO3B)=trcs_coH1(ids_NO3)
      trcs_coH1(ids_NO2B)=trcs_coH1(ids_NO2)
    ENDIF
    IF(VOLH2B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4B)=AZMAX1(trc_soHml2(ids_H1PO4B,N3,N2,N1)/VOLH2B)
      trcs_coH1(ids_H2PO4B)=AZMAX1(trc_soHml2(ids_H2PO4B,N3,N2,N1)/VOLH2B)
    ELSE
      trcs_coH1(ids_H1PO4B)=trcs_coH1(ids_H1PO4)
      trcs_coH1(ids_H2PO4B)=trcs_coH1(ids_H2PO4)
    ENDIF
!excldue NH3 and NH3B
    DO NTG=idg_beg,idg_end-2
      trcs_coH2(NTG)=AZMAX1(trc_soHml2(NTG,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4))
    ENDDO

    VOLHMA=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_NH4,N6,N5,N4)
    IF(VOLHMA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4)=AZMAX1(trc_soHml2(ids_NH4,N6,N5,N4)/VOLHMA)
      trcs_coH2(idg_NH3)=AZMAX1(trc_soHml2(idg_NH3,N6,N5,N4)/VOLHMA)
    ELSE
      trcs_coH2(ids_NH4)=0.0_r8
      trcs_coH2(idg_NH3)=0.0_r8
    ENDIF
    VOLHOA=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_NO3,N6,N5,N4)
    IF(VOLHOA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3)=AZMAX1(trc_soHml2(ids_NO3,N6,N5,N4)/VOLHOA)
      trcs_coH2(ids_NO2)=AZMAX1(trc_soHml2(ids_NO2,N6,N5,N4)/VOLHOA)
    ELSE
      trcs_coH2(ids_NO3)=0.0_r8
      trcs_coH2(ids_NO2)=0.0_r8
    ENDIF
    VOLHPA=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    IF(VOLHPA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4)=AZMAX1(trc_soHml2(ids_H1PO4,N6,N5,N4)/VOLHPA)
      trcs_coH2(ids_H2PO4)=AZMAX1(trc_soHml2(ids_H2PO4,N6,N5,N4)/VOLHPA)
    ELSE
      trcs_coH2(ids_H1PO4)=0.0_r8
      trcs_coH2(ids_H2PO4)=0.0_r8
    ENDIF
    VOLHMB=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_NH4B,N6,N5,N4)
    IF(VOLHMB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4B)=AZMAX1(trc_soHml2(ids_NH4B,N6,N5,N4)/VOLHMB)
      trcs_coH2(idg_NH3B)=AZMAX1(trc_soHml2(idg_NH3B,N6,N5,N4)/VOLHMB)
    ELSE
      trcs_coH2(ids_NH4B)=trcs_coH2(ids_NH4)
      trcs_coH2(idg_NH3B)=trcs_coH2(idg_NH3)
    ENDIF
    VOLHOB=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_NO3B,N6,N5,N4)
    IF(VOLHOB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3B)=AZMAX1(trc_soHml2(ids_NO3B,N6,N5,N4)/VOLHOB)
      trcs_coH2(ids_NO2B)=AZMAX1(trc_soHml2(ids_NO2B,N6,N5,N4)/VOLHOB)
    ELSE
      trcs_coH2(ids_NO3B)=trcs_coH2(ids_NO3)
      trcs_coH2(ids_NO2B)=trcs_coH2(ids_NO2)
    ENDIF
    VOLHPB=VLWatMacPM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    IF(VOLHPB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4B)=AZMAX1(trc_soHml2(ids_H1PO4B,N6,N5,N4)/VOLHPB)
      trcs_coH2(ids_H2PO4B)=AZMAX1(trc_soHml2(ids_H2PO4B,N6,N5,N4)/VOLHPB)
    ELSE
      trcs_coH2(ids_H1PO4B)=trcs_coH2(ids_H1PO4)
      trcs_coH2(ids_H2PO4B)=trcs_coH2(ids_H2PO4)
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
!
!     DLYR=soil layer thickness
!     TortMacPM=macropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     WaterFlow2MacPM=water flux through soil macropore from watsub.f
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
    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TortMacPM(M,N3,N2,N1)*DLYR1+TortMacPM(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MacPM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DIFOC=(OCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFON=(ONSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOP=(OPSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOA=(OASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DO NTS=ids_beg,ids_end
      SDifc(NTS)=(SolDifcc(NTS,N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    ENDDO

!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MACROPORES
!
    D9785: DO K=1,jcplx
      DFHOC(K)=DIFOC*(COQCH1(K)-COQCH2(K))
      DFHON(K)=DIFON*(COQNH1(K)-COQNH2(K))
      DFHOP(K)=DIFOP*(COQPH1(K)-COQPH2(K))
      DFHOA(K)=DIFOA*(COQAH1(K)-COQAH2(K))
    ENDDO D9785
! exclude NH3 and NH3B
    DO NTG=idg_beg,idg_end-2
      SDifHFlx(NTG)=SDifc(NTG)*(trcs_coH1(NTG)-trcs_coH2(NTG))
    ENDDO

    DO NTS=ids_nuts_beg,ids_end
      SDifHFlx(NTS)=SDifc(NTS)*(trcs_coH1(NTS)-trcs_coH2(NTS)) &
        *AMIN1(trcs_VLN(NTS,N3,N2,N1),trcs_VLN(NTS,N6,N5,N4))
    ENDDO
  ELSE
    D9780: DO K=1,jcplx
      DFHOC(K)=0.0_r8
      DFHON(K)=0.0_r8
      DFHOP(K)=0.0_r8
      DFHOA(K)=0.0_r8
    ENDDO D9780
    SDifHFlx(ids_beg:ids_end)=0.0_r8
  ENDIF
  DO NTS = ids_beg,ids_end
    R3PoreSoHFlx(NTS,N,N6,N5,N4)=R3PoreSoHFlx(NTS,N,N6,N5,N4)+SDifHFlx(NTS)
  ENDDO
  end subroutine SoluteDifTranspMacropore

! ----------------------------------------------------------------------
  subroutine SoluteAdvDifTransport(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: THETW1(JZ,JY,JX)
  integer :: K,NTS

  THETW1(N3,N2,N1)=AZMAX1(safe_adb(VLWatMicPM(M,N3,N2,N1),VLSoilMicP(N3,N2,N1)))
  THETW1(N6,N5,N4)=AZMAX1(safe_adb(VLWatMicPM(M,N6,N5,N4),VLSoilMicP(N6,N5,N4)))

!     SOLUTE TRANSPORT IN MICROPORES
!
  call SoluteAdvTranspMicropore(M,N,N1,N2,N3,N4,N5,N6)
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MICROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  call SoluteDifTranspMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
!
!     SOLUTE TRANSPORT IN MACROPORES
!
  call SoluteAdvTranspMacropore(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifTranspMacropore(M,N,N1,N2,N3,N4,N5,N6)
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
    ROCFLS(K,N,N6,N5,N4)=RFLOC(K)+DFVOC(K)
    RONFLS(K,N,N6,N5,N4)=RFLON(K)+DFVON(K)
    ROPFLS(K,N,N6,N5,N4)=RFLOP(K)+DFVOP(K)
    ROAFLS(K,N,N6,N5,N4)=RFLOA(K)+DFVOA(K)
    ROCFHS(K,N,N6,N5,N4)=RFHOC(K)+DFHOC(K)
    RONFHS(K,N,N6,N5,N4)=RFHON(K)+DFHON(K)
    ROPFHS(K,N,N6,N5,N4)=RFHOP(K)+DFHOP(K)
    ROAFHS(K,N,N6,N5,N4)=RFHOA(K)+DFHOA(K)
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
    XOCFLS(K,N,N6,N5,N4)=XOCFLS(K,N,N6,N5,N4)+ROCFLS(K,N,N6,N5,N4)
    XONFLS(K,N,N6,N5,N4)=XONFLS(K,N,N6,N5,N4)+RONFLS(K,N,N6,N5,N4)
    XOPFLS(K,N,N6,N5,N4)=XOPFLS(K,N,N6,N5,N4)+ROPFLS(K,N,N6,N5,N4)
    XOAFLS(K,N,N6,N5,N4)=XOAFLS(K,N,N6,N5,N4)+ROAFLS(K,N,N6,N5,N4)
    XOCFHS(K,N,N6,N5,N4)=XOCFHS(K,N,N6,N5,N4)+ROCFHS(K,N,N6,N5,N4)
    XONFHS(K,N,N6,N5,N4)=XONFHS(K,N,N6,N5,N4)+RONFHS(K,N,N6,N5,N4)
    XOPFHS(K,N,N6,N5,N4)=XOPFHS(K,N,N6,N5,N4)+ROPFHS(K,N,N6,N5,N4)
    XOAFHS(K,N,N6,N5,N4)=XOAFHS(K,N,N6,N5,N4)+ROAFHS(K,N,N6,N5,N4)
  ENDDO D9755

  DO NTS=ids_beg,ids_end
    trcs_XFLS(NTS,N,N6,N5,N4)=trcs_XFLS(NTS,N,N6,N5,N4)+R3PoreSolFlx(NTS,N,N6,N5,N4)
    trcs_XFHS(NTS,N,N6,N5,N4)=trcs_XFHS(NTS,N,N6,N5,N4)+R3PoreSoHFlx(NTS,N,N6,N5,N4)
  ENDDO

  end subroutine SoluteAdvDifTransport

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteExchange(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  integer :: K,NTS

  VLWatMicPCO(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*GasSolbility(idg_CO2,N6,N5,N4)
  VLWatMicPCH(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*GasSolbility(idg_CH4,N6,N5,N4)
  VLWatMicPOX(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*GasSolbility(idg_O2,N6,N5,N4)
  VLWatMicPNG(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*GasSolbility(idg_N2,N6,N5,N4)
  VLWatMicPN2(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*GasSolbility(idg_N2O,N6,N5,N4)
  VLWatMicPN3(N6,N5,N4)=VLWatMicPMA(N6,N5,N4)*GasSolbility(idg_NH3,N6,N5,N4)
  VLWatMicPHG(N6,N5,N4)=VLWatMicPM(M,N6,N5,N4)*GasSolbility(idg_H2,N6,N5,N4)

  VLWatMicPNB(N6,N5,N4)=VLWatMicPMB(N6,N5,N4)*GasSolbility(idg_NH3,N6,N5,N4)
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
!     X*FXS,X*FXB= hourly convective + diffusive solute flux between macro- and micropore in non-band,band
!     R*FXS,R*FXB=convective + diffusive solute flux between macro- and micropore in non-band,band
!
  D9945: DO K=1,jcplx
    XOCFXS(K,N6,N5,N4)=XOCFXS(K,N6,N5,N4)+ROCFXS(K,N6,N5,N4)
    XONFXS(K,N6,N5,N4)=XONFXS(K,N6,N5,N4)+RONFXS(K,N6,N5,N4)
    XOPFXS(K,N6,N5,N4)=XOPFXS(K,N6,N5,N4)+ROPFXS(K,N6,N5,N4)
    XOAFXS(K,N6,N5,N4)=XOAFXS(K,N6,N5,N4)+ROAFXS(K,N6,N5,N4)
  ENDDO D9945

  DO NTS=ids_beg,ids_end
    trcs_XFXS(NTS,N6,N5,N4)=trcs_XFXS(NTS,N6,N5,N4)+RporeSoXFlx(NTS,N6,N5,N4)
  ENDDO

  end subroutine MicMacPoresSoluteExchange

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteAdvExchange(M,N,N1,N2,N3,N4,N5,N6)
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: VFLW
  real(r8) :: trc_RFL(ids_beg:ids_end)
  integer :: K,NTG,NTS
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FWatExMacP2MicPM=macro-micropore water transfer from watsub.f
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
  IF(FWatExMacP2MicPM(M,N6,N5,N4).GT.0.0_r8)THEN
    IF(VLWatMacPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM(M,N6,N5,N4)/VLWatMacPM(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF
    D9970: DO K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQCH2(K,N6,N5,N4))
      RFLON(K)=VFLW*AZMAX1(OQNH2(K,N6,N5,N4))
      RFLOP(K)=VFLW*AZMAX1(OQPH2(K,N6,N5,N4))
      RFLOA(K)=VFLW*AZMAX1(OQAH2(K,N6,N5,N4))
    ENDDO D9970

    DO NTG=idg_beg,idg_end-2
      trc_RFL(NTG)=VFLW*AZMAX1(trc_soHml2(NTG,N6,N5,N4))
    ENDDO

    DO NTS=ids_nuts_beg,ids_nuts_end
      trc_RFL(NTS)=VFLW*AZMAX1(trc_soHml2(NTS,N6,N5,N4))*trcs_VLN(NTS,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FWatExMacP2MicPM(M,N6,N5,N4).LT.0.0_r8)THEN
    IF(VLWatMicPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM(M,N6,N5,N4)/VLWatMicPM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    D9965: DO K=1,jcplx
      RFLOC(K)=VFLW*AZMAX1(OQC2(K,N6,N5,N4))
      RFLON(K)=VFLW*AZMAX1(OQN2(K,N6,N5,N4))
      RFLOP(K)=VFLW*AZMAX1(OQP2(K,N6,N5,N4))
      RFLOA(K)=VFLW*AZMAX1(OQA2(K,N6,N5,N4))
    ENDDO D9965
!exclude NH3 and NH3B
    DO NTG=idg_beg,idg_end-2
      trc_RFL(NTG)=VFLW*AZMAX1(trc_solml2(NTG,N6,N5,N4))
    ENDDO

    DO NTS=ids_nuts_beg,ids_nuts_end
      trc_RFL(NTS)=VFLW*AZMAX1(trc_solml2(NTS,N6,N5,N4))*trcs_VLN(NTS,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    D9960: DO K=1,jcplx
      RFLOC(K)=0.0_r8
      RFLON(K)=0.0_r8
      RFLOP(K)=0.0_r8
      RFLOA(K)=0.0_r8
    ENDDO D9960
    trc_RFL(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO NTS=ids_beg,ids_end
    RporeSoXFlx(NTS,N6,N5,N4)=trc_RFL(NTS)
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
    ROCFXS(K,N6,N5,N4)=RFLOC(K)
    RONFXS(K,N6,N5,N4)=RFLON(K)
    ROPFXS(K,N6,N5,N4)=RFLOP(K)
    ROAFXS(K,N6,N5,N4)=RFLOA(K)
  enddo

  end subroutine MicMacPoresSoluteAdvExchange

! ----------------------------------------------------------------------

  subroutine MicMacPoresSoluteDifExchange(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: trcs_DFV(ids_beg:ids_end)
  real(r8) :: VLWatMacPS,VOLWT
  integer  :: K,NTS,NTG
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
  IF(VLWatMacPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
    VLWatMacPS=AMIN1(XFRS*VGeomLayer(N6,N5,N4),VLWatMacPM(M,N6,N5,N4))
    VOLWT=VLWatMicPM(M,N6,N5,N4)+VLWatMacPS

    D9955: DO K=1,jcplx
      DFVOC(K)=dts_HeatWatTP*(AZMAX1(OQCH2(K,N6,N5,N4))*VLWatMicPM(M,N6,N5,N4) &
        -AZMAX1(OQC2(K,N6,N5,N4))*VLWatMacPS)/VOLWT
      DFVON(K)=dts_HeatWatTP*(AZMAX1(OQNH2(K,N6,N5,N4))*VLWatMicPM(M,N6,N5,N4) &
        -AZMAX1(OQN2(K,N6,N5,N4))*VLWatMacPS)/VOLWT
      DFVOP(K)=dts_HeatWatTP*(AZMAX1(OQPH2(K,N6,N5,N4))*VLWatMicPM(M,N6,N5,N4) &
        -AZMAX1(OQP2(K,N6,N5,N4))*VLWatMacPS)/VOLWT
      DFVOA(K)=dts_HeatWatTP*(AZMAX1(OQAH2(K,N6,N5,N4))*VLWatMicPM(M,N6,N5,N4) &
        -AZMAX1(OQA2(K,N6,N5,N4))*VLWatMacPS)/VOLWT
    ENDDO D9955

    DO NTG=idg_beg,idg_end-2
      trcs_DFV(NTG)=dts_HeatWatTP*(AZMAX1(trc_soHml2(NTG,N6,N5,N4))*VLWatMicPM(M,N6,N5,N4) &
        -AZMAX1(trc_solml2(NTG,N6,N5,N4))*VLWatMacPS)/VOLWT
    ENDDO

    DO NTS=ids_nuts_beg,ids_nuts_end
      trcs_DFV(NTS)=dts_HeatWatTP*(AZMAX1(trc_soHml2(NTS,N6,N5,N4))*VLWatMicPM(M,N6,N5,N4) &
        -AZMAX1(trc_solml2(NTS,N6,N5,N4))*VLWatMacPS)/VOLWT &
        *trcs_VLN(NTS,N6,N5,N4)
    ENDDO

  ELSE
    D9975: DO K=1,jcplx
      DFVOC(K)=0.0_r8
      DFVON(K)=0.0_r8
      DFVOP(K)=0.0_r8
      DFVOA(K)=0.0_r8
    ENDDO D9975

    trcs_DFV(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO NTS=ids_beg,ids_end
    RporeSoXFlx(NTS,N6,N5,N4)=RporeSoXFlx(NTS,N6,N5,N4)+trcs_DFV(NTS)
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
    ROCFXS(K,N6,N5,N4)=ROCFXS(K,N6,N5,N4)+DFVOC(K)
    RONFXS(K,N6,N5,N4)=RONFXS(K,N6,N5,N4)+DFVON(K)
    ROPFXS(K,N6,N5,N4)=ROPFXS(K,N6,N5,N4)+DFVOP(K)
    ROAFXS(K,N6,N5,N4)=ROAFXS(K,N6,N5,N4)+DFVOA(K)
  enddo
  end subroutine MicMacPoresSoluteDifExchange

! ----------------------------------------------------------------------
  subroutine BubbleEfflux(M,N1,N2,N3,NY,NX,MX,IFLGB)
  implicit none
  integer, intent(in) :: M,N1,N2,N3,NY,NX,MX
  integer, intent(inout) :: IFLGB
  real(r8) :: THETW1
  real(r8) :: trcg_SLX(idg_beg:idg_end)
  real(r8) :: trcg_VOLG(idg_beg:idg_end)
  integer  :: NTG
  real(r8) :: VTATM,VTGAS,DVTGAS
!
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     VLSoilMicP=micropore volume
!     IFLGB=bubbling flag:0=enabled,1=disabled
!     S*L=solubility of gas in water from hour1.f
!
  IF(N3.GE.NUM(N2,N1).AND.M.NE.MX)THEN
    THETW1=AZMAX1(safe_adb(VLWatMicPM(M,N3,N2,N1),VLSoilMicP(N3,N2,N1)))
    IF(THETW1.GT.THETY(N3,N2,N1).AND.IFLGB.EQ.0)THEN

      trcg_SLX(idg_CO2) =catomw*GasSolbility(idg_CO2,N3,N2,N1)  !conver into carbon g C/mol
      trcg_SLX(idg_CH4) =catomw*GasSolbility(idg_CH4,N3,N2,N1)
      trcg_SLX(idg_O2)  =32.0_r8*GasSolbility(idg_O2,N3,N2,N1)
      trcg_SLX(idg_N2)  =28.0_r8*GasSolbility(idg_N2,N3,N2,N1)
      trcg_SLX(idg_N2O) =28.0_r8*GasSolbility(idg_N2O,N3,N2,N1)
      trcg_SLX(idg_NH3) =natomw*GasSolbility(idg_NH3,N3,N2,N1)
      trcg_SLX(idg_H2)  =2.0_r8*GasSolbility(idg_H2,N3,N2,N1)
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
        DO NTG=idg_beg,idg_end
          trcg_VOLG(NTG)=(trc_solml2(NTG,N3,N2,N1)+RGasSSVol(NTG,NY,NX))/trcg_SLX(NTG)
        ENDDO
      ELSE
        DO NTG=idg_beg,idg_end
          trcg_VOLG(NTG)=trc_solml2(NTG,N3,N2,N1)/trcg_SLX(NTG)
        ENDDO
      ENDIF
!
!     GASEOUS EQUIVALENT ATMOSPHERIC CONCENTRATION
!
!     VTATM=molar gas concentration at atmospheric pressure
!     VTGAS=total molar gas concentration
!
      VTATM=AZMAX1(1.2194E+04_r8*VLWatMicPM(M,N3,N2,N1)/TKS(N3,N2,N1))

      VTGAS=sum(trcg_VOLG(idg_beg:idg_end))
!
!     PROPORTIONAL REMOVAL OF EXCESS AQUEOUS GASES
!
!     R*BBL=bubble flux
!     gas code:*CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     V*G2=molar gas concentration
!
      IF(VTGAS.GT.VTATM)THEN
        DVTGAS=0.5_r8*(VTATM-VTGAS)
        DO NTG=idg_beg,idg_end
          trcg_BBL(NTG,N3,N2,N1)=AZMIN1(DVTGAS*trcg_VOLG(NTG)/VTGAS)*trcg_SLX(NTG)
        ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*BBL=hourly bubble flux
!
        DO NTG=idg_beg,idg_end
          trcg_XBLL(NTG,N3,N2,N1)=trcg_XBLL(NTG,N3,N2,N1)+trcg_BBL(NTG,N3,N2,N1)
        ENDDO
      ELSE
        trcg_BBL(idg_beg:idg_end,N3,N2,N1)=0.0_r8
      ENDIF
    ELSE
      IFLGB=1
      trcg_BBL(idg_beg:idg_end,N3,N2,N1)=0.0_r8
    ENDIF

  ENDIF
  end subroutine BubbleEfflux
! ----------------------------------------------------------------------

  subroutine GasDifTransport(M,N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: trc_gasc1(idg_beg:idg_end)
  real(r8) :: trc_gasc2(idg_beg:idg_end)

  real(r8) :: CNDC1
  real(r8) :: CNDC2
  real(r8) :: DFLG2,DFLGL
  real(r8) :: DFVGG(idg_beg:idg_end)
  integer :: NTG

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
  DFLG2=2.0_r8*AZMAX1(THETPM(M,N3,N2,N1))*POROQ*THETPM(M,N3,N2,N1)/POROS(N3,N2,N1) &
    *AREA(N,N3,N2,N1)/DLYR(N,N3,N2,N1)

  DFLGL=2.0_r8*AZMAX1(THETPM(M,N6,N5,N4))*POROQ*THETPM(M,N6,N5,N4)/POROS(N6,N5,N4) &
    *AREA(N,N6,N5,N4)/DLYR(N,N6,N5,N4)

!
!     GASOUS CONDUCTANCES
!
!     D*G=gaseous diffusivity in soil
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
  DO NTG=idg_beg,idg_end
    CNDC1=DFLG2*GasDifcc(NTG,N3,N2,N1)
    CNDC2=DFLGL*GasDifcc(NTG,N6,N5,N4)
    DifuscG(NTG,N,N6,N5,N4)=(CNDC1*CNDC2)/(CNDC1+CNDC2)
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
  DO NTG=idg_beg,idg_end-1
    trc_gasc1(NTG)=AZMAX1(trc_gasml2(NTG,N3,N2,N1)/VLsoiAirPM(M,N3,N2,N1))
    trc_gasc2(NTG)=AZMAX1(trc_gasml2(NTG,N6,N5,N4)/VLsoiAirPM(M,N6,N5,N4))
    DFVGG(NTG)=DifuscG(NTG,N,N6,N5,N4)*(trc_gasc1(NTG)-trc_gasc2(NTG))
    R3GasADFlx(NTG,N,N6,N5,N4)=DFVGG(NTG)
  ENDDO

  end subroutine GasDifTransport


! ----------------------------------------------------------------------
  subroutine GasAdvTransport(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)
  implicit none
  real(r8), intent(in) :: WaterFlow2Soil(3,JD,JV,JH)
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: RGasAdv
  real(r8) :: VFLW,FLQW
  integer :: NTG
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
    IF(VLsoiAirPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=-AZMAX1(AMIN1(VFLWX,FLQW/VLsoiAirPM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO NTG=idg_beg,idg_end-1
      RGasAdv=VFLW*AZMAX1(trc_gasml2(NTG,N6,N5,N4))
      R3GasADFlx(NTG,N,N6,N5,N4)=R3GasADFlx(NTG,N,N6,N5,N4)+RGasAdv
    ENDDO
  ELSE
    IF(VLsoiAirPM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=-AZMIN1(AMAX1(-VFLWX,FLQW/VLsoiAirPM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO NTG=idg_beg,idg_end-1
      RGasAdv=VFLW*AZMAX1(trc_gasml2(NTG,N3,N2,N1))
      R3GasADFlx(NTG,N,N6,N5,N4)=R3GasADFlx(NTG,N,N6,N5,N4)+RGasAdv
    ENDDO
  ENDIF

  end subroutine GasAdvTransport

! ----------------------------------------------------------------------
  subroutine GaseousTransport(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)

  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: WaterFlow2Soil(3,JD,JV,JH)
  integer :: NTG

!     THETPM,VLsoiAirPM=air-filled porosity,volume from watsub.f

  IF(THETPM(M,N3,N2,N1).GT.THETX.AND.THETPM(M,N6,N5,N4).GT.THETX &
    .AND.VLsoiAirPM(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
    .AND.VLsoiAirPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN

!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE
    call GasDifTransport(M,N,N1,N2,N3,N4,N5,N6)

!     TOTAL SOIL GAS FLUX FROM CONVECTIVE FLUX
    call GasAdvTransport(M,N,N1,N2,N3,N4,N5,N6,WaterFlow2Soil)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLG=hourly total convective+diffusive gas flux
!
    DO NTG=idg_beg,idg_end-1
      R3GasADTFlx(NTG,N,N6,N5,N4)=R3GasADTFlx(NTG,N,N6,N5,N4) &
        +R3GasADFlx(NTG,N,N6,N5,N4)
    ENDDO

  ELSE
    R3GasADFlx(idg_beg:idg_end-1,N,N6,N5,N4)=0.0_r8
  ENDIF
  call VolatilizationDissolution(M,N,N1,N2,N3,N4,N5,N6)
  end subroutine GaseousTransport

! ----------------------------------------------------------------------
  subroutine VolatilizationDissolution(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: CNH3B0,CNH4B0,CNH4S0
  real(r8) :: CNH3S0
  integer :: NTG
!
!     VOLATILIZATION-DISSOLUTION OF GASES IN SOIL
!     LAYER FROM GASEOUS CONCENTRATIONS VS. THEIR AQUEOUS
!     EQUIVALENTS DEPENDING ON SOLUBILITY FROM 'HOUR1'
!     AND TRANSFER COEFFICIENT 'DFGS' FROM 'WATSUB'
!
!     THETPM,VLWatMicPPM=air-filled porosity,volume
!     R*DFG=water-air gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     DFGS=rate constant for air-water gas exchange from watsub.f
!     *G2,*S2=gaseous,aqueous gas content
!     VLWatMicP*=equivalent aqueous volume for gas
!
  IF(N.EQ.3)THEN
    IF(THETPM(M,N6,N5,N4).GT.THETX)THEN
      RGasDSFlx(idg_CO2,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
        ,trc_gasml2(idg_CO2,N6,N5,N4))*VLWatMicPCO(N6,N5,N4) &
        -trc_solml2(idg_CO2,N6,N5,N4)*VLsoiAirPM(M,N6,N5,N4)) &
        /(VLWatMicPCO(N6,N5,N4)+VLsoiAirPM(M,N6,N5,N4))
      RGasDSFlx(idg_CH4,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
        ,trc_gasml2(idg_CH4,N6,N5,N4))*VLWatMicPCH(N6,N5,N4) &
        -trc_solml2(idg_CH4,N6,N5,N4)*VLsoiAirPM(M,N6,N5,N4)) &
        /(VLWatMicPCH(N6,N5,N4)+VLsoiAirPM(M,N6,N5,N4))
      RGasDSFlx(idg_O2,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
        ,trc_gasml2(idg_O2,N6,N5,N4))*VLWatMicPOX(N6,N5,N4) &
        -trc_solml2(idg_O2,N6,N5,N4)*VLsoiAirPM(M,N6,N5,N4)) &
        /(VLWatMicPOX(N6,N5,N4)+VLsoiAirPM(M,N6,N5,N4))
      RGasDSFlx(idg_N2,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
        ,trc_gasml2(idg_N2,N6,N5,N4))*VLWatMicPNG(N6,N5,N4) &
        -trc_solml2(idg_N2,N6,N5,N4)*VLsoiAirPM(M,N6,N5,N4)) &
        /(VLWatMicPNG(N6,N5,N4)+VLsoiAirPM(M,N6,N5,N4))
      RGasDSFlx(idg_N2O,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
        ,trc_gasml2(idg_N2O,N6,N5,N4))*VLWatMicPN2(N6,N5,N4) &
        -trc_solml2(idg_N2O,N6,N5,N4)*VLsoiAirPM(M,N6,N5,N4)) &
        /(VLWatMicPN2(N6,N5,N4)+VLsoiAirPM(M,N6,N5,N4))

      IF(VLsoiAirPMA(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXA(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGasDSFlx(idg_NH3,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
          ,trc_gasml2(idg_NH3,N6,N5,N4))*VLWatMicPN3(N6,N5,N4) &
          -trc_solml2(idg_NH3,N6,N5,N4)*VLsoiAirPMA(N6,N5,N4)) &
          /(VLWatMicPN3(N6,N5,N4)+VLsoiAirPMA(N6,N5,N4))
        CNH3S0=AZMAX1((trc_solml2(idg_NH3,N6,N5,N4) &
          +RGasDSFlx(idg_NH3,N6,N5,N4))/VLWatMicPXA(N6,N5,N4))
        CNH4S0=AZMAX1(trc_solml2(ids_NH4,N6,N5,N4))/VLWatMicPXA(N6,N5,N4)
      ELSE
        RGasDSFlx(idg_NH3,N6,N5,N4)=0.0_r8
      ENDIF

      RGasDSFlx(idg_H2,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
        ,trc_gasml2(idg_H2,N6,N5,N4))*VLWatMicPHG(N6,N5,N4) &
        -trc_solml2(idg_H2,N6,N5,N4)*VLsoiAirPM(M,N6,N5,N4)) &
        /(VLWatMicPHG(N6,N5,N4)+VLsoiAirPM(M,N6,N5,N4))


      IF(VLsoiAirPMB(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXB(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGasDSFlx(idg_NH3B,N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
          ,trc_gasml2(idg_NH3,N6,N5,N4))*VLWatMicPNB(N6,N5,N4) &
          -trc_solml2(idg_NH3B,N6,N5,N4)*VLsoiAirPMB(N6,N5,N4)) &
          /(VLWatMicPNB(N6,N5,N4)+VLsoiAirPMB(N6,N5,N4))
        CNH3B0=AZMAX1((trc_solml2(idg_NH3B,N6,N5,N4) &
          +RGasDSFlx(idg_NH3B,N6,N5,N4))/VLWatMicPXB(N6,N5,N4))
        CNH4B0=AZMAX1(trc_solml2(ids_NH4B,N6,N5,N4))/VLWatMicPXB(N6,N5,N4)
      ELSE
        RGasDSFlx(idg_NH3B,N6,N5,N4)=0.0_r8
      ENDIF

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly water-air gas flux
!
      DO NTG=idg_beg,idg_end
        GasDisFlx(NTG,N6,N5,N4)=GasDisFlx(NTG,N6,N5,N4)+RGasDSFlx(NTG,N6,N5,N4)
      ENDDO
    ELSE
      RGasDSFlx(idg_beg:idg_end,N6,N5,N4)=0.0_r8
    ENDIF
  ENDIF
  end subroutine VolatilizationDissolution

! ----------------------------------------------------------------------
  subroutine ZeroTransport1(N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,N6

  integer :: K

  DifuscG(idg_beg:idg_end,N,N6,N5,N4) = 0._r8

  D9750: DO K=1,jcplx
    ROCFLS(K,N,N6,N5,N4)=0.0_r8
    RONFLS(K,N,N6,N5,N4)=0.0_r8
    ROPFLS(K,N,N6,N5,N4)=0.0_r8
    ROAFLS(K,N,N6,N5,N4)=0.0_r8
    ROCFHS(K,N,N6,N5,N4)=0.0_r8
    RONFHS(K,N,N6,N5,N4)=0.0_r8
    ROPFHS(K,N,N6,N5,N4)=0.0_r8
    ROAFHS(K,N,N6,N5,N4)=0.0_r8
  ENDDO D9750

  R3PoreSolFlx(ids_beg:ids_end,N,N6,N5,N4)=0.0_r8

  R3PoreSoHFlx(ids_beg:ids_end,N,N6,N5,N4)=0.0_r8

  R3GasADFlx(idg_beg:idg_end-1,N,N6,N5,N4)=0.0_r8
  end subroutine ZeroTransport1
! ----------------------------------------------------------------------

  subroutine ZeroTransport2(N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,N6

  integer :: K

  DifuscG(idg_beg:idg_end-1,N,N3,N2,N1)=0.0_r8

  D9751: DO K=1,jcplx
    ROCFLS(K,N,N3,N2,N1)=0.0_r8
    RONFLS(K,N,N3,N2,N1)=0.0_r8
    ROPFLS(K,N,N3,N2,N1)=0.0_r8
    ROAFLS(K,N,N3,N2,N1)=0.0_r8
    ROCFHS(K,N,N3,N2,N1)=0.0_r8
    RONFHS(K,N,N3,N2,N1)=0.0_r8
    ROPFHS(K,N,N3,N2,N1)=0.0_r8
    ROAFHS(K,N,N3,N2,N1)=0.0_r8
  ENDDO D9751

  R3PoreSolFlx(ids_beg:ids_end,N,N3,N2,N1)=0.0_r8
  R3PoreSoHFlx(ids_beg:ids_end,N,N3,N2,N1)=0.0_r8
  R3GasADFlx(idg_beg:idg_end-1,N,N3,N2,N1)=0.0_r8
  end subroutine ZeroTransport2
end module InsideTranspMod
