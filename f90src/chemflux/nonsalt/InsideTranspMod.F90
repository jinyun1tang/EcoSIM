module InsideTranspMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : destroy
  use minimathmod, only : safe_adb
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
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  real(r8), parameter :: XFRS=0.05
  real(r8),allocatable :: RFLOC(:)
  real(r8),allocatable :: RFLON(:)
  real(r8),allocatable :: RFLOP(:)
  real(r8),allocatable :: RFLOA(:)
  real(r8),allocatable :: RFHOC(:)
  real(r8),allocatable :: RFHON(:)
  real(r8),allocatable :: RFHOP(:)
  real(r8),allocatable :: RFHOA(:)
  real(r8),allocatable :: COQC1(:)
  real(r8),allocatable :: COQC2(:)
  real(r8),allocatable :: COQN1(:)
  real(r8),allocatable :: COQN2(:)
  real(r8),allocatable :: COQP1(:)
  real(r8),allocatable :: COQP2(:)
  real(r8),allocatable :: COQA1(:)
  real(r8),allocatable :: COQA2(:)
  real(r8),allocatable :: COQCH1(:)
  real(r8),allocatable :: COQCH2(:)
  real(r8),allocatable :: COQNH1(:)
  real(r8),allocatable :: COQNH2(:)
  real(r8),allocatable :: COQPH1(:)
  real(r8),allocatable :: COQPH2(:)
  real(r8),allocatable :: COQAH1(:)
  real(r8),allocatable :: COQAH2(:)
  real(r8),allocatable :: DFVOC(:)
  real(r8),allocatable :: DFVON(:)
  real(r8),allocatable :: DFVOP(:)
  real(r8),allocatable :: DFVOA(:)
  real(r8),allocatable :: DFHOC(:)
  real(r8),allocatable :: DFHON(:)
  real(r8),allocatable :: DFHOP(:)
  real(r8),allocatable :: DFHOA(:)

  real(r8) :: ROXFLS1,ROXFLS0,RNOFLW1,RNGFLS0
  real(r8) :: RCOFLS0,RCHFLS0,RN2FLS0,RN4FLW0
  real(r8) :: RN3FLW0,RNOFLW0,RH1PFS0,RH2PFS0,RCOFLS1,RCHFLS1
  real(r8) :: RNGFLS1,RN2FLS1,RN4FLW1,RN3FLW1
  real(r8) :: RH1PFS1,RH2PFS1,RN4FLB1,RN3FLB1,RNOFLB1,RH1BFB1
  real(r8) :: RH2BFB1,RCODXR,RCHDXR,ROXDXR,RNGDXR,RN2DXR,RN3DXR
  real(r8) :: RHGDXR,RCODXS,RCHDXS,ROXDXS,RNGDXS,RN2DXS,RN3DXS
  real(r8) :: RNBDXS,RHGDXS,RFLCOS,RFLCHS,RFLOXS,RFLNGS,RFLN2S
  real(r8) :: RFLHGS,RFLNH4,RFLNH3,RFLNO3,RFLNO2,RFLP14,RFLPO4
  real(r8) :: RFLN4B,RFLN3B,RFLNOB,RFLN2B,RFLP1B,RFLPOB,RFLCOG
  real(r8) :: RFLCHG,RFLOXG,RFLNGG,RFLN2G,RFLN3G,RFLH2G,RFHCOS
  real(r8) :: RFHCHS,RFHOXS,RFHNGS,RFHN2S,RFHHGS,RFHNH4,RFHNH3
  real(r8) :: RFHNO3,RFHNO2,RFHP14,RFHPO4,RFHN4B,RFHN3B,RFHNOB
  real(r8) :: RFHN2B,RFHP1B,RFHPOB,RCHGFU,RCHGFT
  real(r8) :: CCH4S1,CCH4S2,CNO3B2
  real(r8) :: COXYS1,CZ2GS1,CZ2OS1,CH2GS1,CNH4S1
  real(r8) :: CNH3S1,CNO3S1,CNO2S1,CP14S1,CPO4S1
  real(r8) :: CCO2S1,CCO2S2
  real(r8) :: COXYS2,CZ2GS2,CZ2OS2,CH2GS2,CNH3S2,CNH4S2,CNO3S2
  real(r8) :: CNO2S2,CP14S2,CPO4S2
  real(r8) :: CNH3B2,CNH4B2,CNO2B2
  real(r8) :: CP14B2,CPO4B2
  real(r8) :: DFVCOS,DFVCHS,DFVOXS
  real(r8) :: DFVNGS,DFVN2S,DFVHGS,DFVNH4,DFVNH3,DFVNO3,DFVNO2
  real(r8) :: DFVP14,DFVPO4,DFVN4B,DFVN3B,DFVNOB,DFVN2B,DFVP1B
  real(r8) :: DFVPOB,DFVCOG,DFVCHG,DFVOXG,DFVNGG


  public :: ModelSoluteHydroFlux
  public :: InitInsTp
  public :: DestructInsTp
  contains
!------------------------------------------------------------------------------------------
  subroutine InitInsTp()
  implicit none
  allocate(RFLOC(0:jcplx1))
  allocate(RFLON(0:jcplx1))
  allocate(RFLOP(0:jcplx1))
  allocate(RFLOA(0:jcplx1))
  allocate(RFHOC(0:jcplx1))
  allocate(RFHON(0:jcplx1))
  allocate(RFHOP(0:jcplx1))
  allocate(RFHOA(0:jcplx1))
  allocate(COQC1(0:jcplx1))
  allocate(COQC2(0:jcplx1))
  allocate(COQN1(0:jcplx1))
  allocate(COQN2(0:jcplx1))
  allocate(COQP1(0:jcplx1))
  allocate(COQP2(0:jcplx1))
  allocate(COQA1(0:jcplx1))
  allocate(COQA2(0:jcplx1))
  allocate(COQCH1(0:jcplx1))
  allocate(COQCH2(0:jcplx1))
  allocate(COQNH1(0:jcplx1))
  allocate(COQNH2(0:jcplx1))
  allocate(COQPH1(0:jcplx1))
  allocate(COQPH2(0:jcplx1))
  allocate(COQAH1(0:jcplx1))
  allocate(COQAH2(0:jcplx1))
  allocate(DFVOC(0:jcplx1))
  allocate(DFVON(0:jcplx1))
  allocate(DFVOP(0:jcplx1))
  allocate(DFVOA(0:jcplx1))
  allocate(DFHOC(0:jcplx1))
  allocate(DFHON(0:jcplx1))
  allocate(DFHOP(0:jcplx1))
  allocate(DFHOA(0:jcplx1))

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

  subroutine ModelSoluteHydroFlux(M,MX, NHW, NHE, NVN, NVS)
  implicit none

  integer, intent(in) :: M,MX, NHW, NHE, NVN, NVS
  integer :: NY,NX
  real(r8) :: FLWRM1
  real(r8) :: FLQM(3,JD,JV,JH)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      call ResetandInitFluxAccumulators(M,NY,NX,MX)

      IF(M.NE.MX)THEN
!     This IF statement is the next ~1700 lines so I'm leaving it in here
!
!     SOLUTE FLUXES FROM MELTING SNOWPACK TO
!     RESIDUE AND SOIL SURFACE FROM SNOWMELT IN 'WATSUB' AND
!     CONCENTRATIONS IN SNOWPACK
!
        call SoluteFluxSnowpack(M,NY,NX)
!
!     SOLUTE FLUXES AT SOIL SURFACE FROM SURFACE WATER
!     CONTENTS, WATER FLUXES 'FLQM' AND ATMOSPHERE BOUNDARY
!     LAYER RESISTANCES 'PARGM' FROM 'WATSUB'
!
        call SoluteFluxSurface(M,NY,NX,FLQM)
!
!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
        call ConvectiveSurfaceSoluteFlux(M,NY,NX,FLWRM1)
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
        call DiffusiveFluxAtSoilSurface(M,NY,NX,FLWRM1)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES
!
        call TotalPoreFluxAdjacentCell(NY,NX)
!
!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     SURFACE LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
        call MacroMicroTransfer(M,NY,NX)
!
!     SOLUTE TRANSPORT FROM WATER OVERLAND FLOW
!     IN 'WATSUB' AND FROM SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE LAYER
!
        call OverlandFlowSnowdriftTransport(M,NY,NX,NHE,NHW,NVS,NVN)
      ENDIF
!
!     VOLATILIZATION-DISSOLUTION OF GASES IN RESIDUE AND SOIL SURFACE
!     LAYERS FROM GASEOUS CONCENTRATIONS VS. THEIR AQUEOUS
!     EQUIVALENTS DEPENDING ON SOLUBILITY FROM 'HOUR1'
!     AND TRANSFER COEFFICIENT 'DFGS' FROM 'WATSUB'
!
      call SurfaceGasVolatilDissol(M,NY,NX)
!
!     SURFACE GAS EXCHANGE FROM GAS DIFFUSIVITY THROUGH
!     SOIL SURFACE LAYER AND THROUGH ATMOSPHERE BOUNDARY
!     LAYER
!
      call GasDiffusionConvection(M,NY,NX,FLQM)
!
!     SOIL SURFACE WATER-AIR GAS EXCHANGE
!
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
!
      call SoluteFluxAdjacentCells(M,MX,NY,NX,NHE,NVS)

    ENDDO
  ENDDO
  end subroutine ModelSoluteHydroFlux
!------------------------------------------------------------------------------------------

  subroutine ResetandInitFluxAccumulators(M,NY,NX,MX)
  implicit none

  integer, intent(in) :: M, NY, NX,MX

  integer :: K,L
  real(r8) :: PARGM

  IF(M.NE.MX)THEN
!
!     GASEOUS BOUNDARY LAYER CONDUCTANCES
!
!     PARG=boundary layer conductance above soil surface from watsub.f
!
    PARGM=PARG(M,NY,NX)*XNPT
    PARGCO(NY,NX)=PARGM*0.74
    PARGCH(NY,NX)=PARGM*1.04
    PARGOX(NY,NX)=PARGM*0.83
    PARGNG(NY,NX)=PARGM*0.86
    PARGN2(NY,NX)=PARGM*0.74
    PARGN3(NY,NX)=PARGM*1.02
    PARGH2(NY,NX)=PARGM*2.08
!
!     RESET RUNOFF SOLUTE FLUX ACCUMULATORS
!
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in micropores
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in micropores
!     ZN3G=gaseous NH3
!

    DO  K=0,jcplx1
      TQROC(K,NY,NX)=0.0
      TQRON(K,NY,NX)=0.0
      TQROP(K,NY,NX)=0.0
      TQROA(K,NY,NX)=0.0
      OQC2(K,0,NY,NX)=OQC2(K,0,NY,NX)-ROCSK2(K,0,NY,NX)
      OQN2(K,0,NY,NX)=OQN2(K,0,NY,NX)-RONSK2(K,0,NY,NX)
      OQP2(K,0,NY,NX)=OQP2(K,0,NY,NX)-ROPSK2(K,0,NY,NX)
      OQA2(K,0,NY,NX)=OQA2(K,0,NY,NX)-ROASK2(K,0,NY,NX)
    ENDDO
    TQRCOS(NY,NX)=0.0
    TQRCHS(NY,NX)=0.0
    TQROXS(NY,NX)=0.0
    TQRNGS(NY,NX)=0.0
    TQRN2S(NY,NX)=0.0
    TQRHGS(NY,NX)=0.0
    TQRNH4(NY,NX)=0.0
    TQRNH3(NY,NX)=0.0
    TQRNO3(NY,NX)=0.0
    TQRNO2(NY,NX)=0.0
    TQRH1P(NY,NX)=0.0
    TQRH2P(NY,NX)=0.0
    TQSCOS(NY,NX)=0.0
    TQSCHS(NY,NX)=0.0
    TQSOXS(NY,NX)=0.0
    TQSNGS(NY,NX)=0.0
    TQSN2S(NY,NX)=0.0
    TQSNH4(NY,NX)=0.0
    TQSNH3(NY,NX)=0.0
    TQSNO3(NY,NX)=0.0
    TQSH1P(NY,NX)=0.0
    TQSH2P(NY,NX)=0.0
    ZNH4S2(0,NY,NX)=ZNH4S2(0,NY,NX)-RN4SK2(0,NY,NX)
    ZNH3S2(0,NY,NX)=ZNH3S2(0,NY,NX)-RN3SK2(0,NY,NX)
    ZNO3S2(0,NY,NX)=ZNO3S2(0,NY,NX)-RNOSK2(0,NY,NX)
    ZNO2S2(0,NY,NX)=ZNO2S2(0,NY,NX)-RNXSK2(0,NY,NX)
    H2PO42(0,NY,NX)=H2PO42(0,NY,NX)-RHPSK2(0,NY,NX)
    H1PO42(0,NY,NX)=H1PO42(0,NY,NX)-R1PSK2(0,NY,NX)
    ROXSK2(0,NY,NX)=ROXSK(M,0,NY,NX)*XNPT
  ENDIF
  CO2S2(0,NY,NX)=CO2S2(0,NY,NX)-RCOSK2(0,NY,NX)
  CH4S2(0,NY,NX)=CH4S2(0,NY,NX)-RCHSK2(0,NY,NX)
  OXYS2(0,NY,NX)=OXYS2(0,NY,NX)-ROXSK2(0,NY,NX)
  Z2GS2(0,NY,NX)=Z2GS2(0,NY,NX)-RNGSK2(0,NY,NX)
  Z2OS2(0,NY,NX)=Z2OS2(0,NY,NX)-RN2SK2(0,NY,NX)
  H2GS2(0,NY,NX)=H2GS2(0,NY,NX)-RHGSK2(0,NY,NX)
  ZN3G2(0,NY,NX)=ZN3G2(0,NY,NX)-RNHSK2(0,NY,NX)
!
!     INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
!
  IF(M.NE.MX)THEN
    DO  L=1,JS
      TCOBLS(L,NY,NX)=0.0
      TCHBLS(L,NY,NX)=0.0
      TOXBLS(L,NY,NX)=0.0
      TNGBLS(L,NY,NX)=0.0
      TN2BLS(L,NY,NX)=0.0
      TN4BLW(L,NY,NX)=0.0
      TN3BLW(L,NY,NX)=0.0
      TNOBLW(L,NY,NX)=0.0
      TH1PBS(L,NY,NX)=0.0
      TH2PBS(L,NY,NX)=0.0
    ENDDO
  ENDIF
!
!     INITIALIZE SOIL SOLUTE NET FLUX ACCUMULATORS
!
  DO L=NU(NY,NX),NL(NY,NX)
    IF(M.NE.MX)THEN
      DO  K=0,jcplx1
        TOCFLS(K,L,NY,NX)=0.0
        TONFLS(K,L,NY,NX)=0.0
        TOPFLS(K,L,NY,NX)=0.0
        TOAFLS(K,L,NY,NX)=0.0
        TOCFHS(K,L,NY,NX)=0.0
        TONFHS(K,L,NY,NX)=0.0
        TOPFHS(K,L,NY,NX)=0.0
        TOAFHS(K,L,NY,NX)=0.0
        OQC2(K,L,NY,NX)=OQC2(K,L,NY,NX)-ROCSK2(K,L,NY,NX)
        OQN2(K,L,NY,NX)=OQN2(K,L,NY,NX)-RONSK2(K,L,NY,NX)
        OQP2(K,L,NY,NX)=OQP2(K,L,NY,NX)-ROPSK2(K,L,NY,NX)
        OQA2(K,L,NY,NX)=OQA2(K,L,NY,NX)-ROASK2(K,L,NY,NX)
      ENDDO
      TCOFLS(L,NY,NX)=0.0
      TCHFLS(L,NY,NX)=0.0
      TOXFLS(L,NY,NX)=0.0
      TNGFLS(L,NY,NX)=0.0
      TN2FLS(L,NY,NX)=0.0
      THGFLS(L,NY,NX)=0.0
      TN4FLW(L,NY,NX)=0.0
      TN3FLW(L,NY,NX)=0.0
      TNOFLW(L,NY,NX)=0.0
      TNXFLS(L,NY,NX)=0.0
      TH1PFS(L,NY,NX)=0.0
      TH2PFS(L,NY,NX)=0.0
      TN4FLB(L,NY,NX)=0.0
      TN3FLB(L,NY,NX)=0.0
      TNOFLB(L,NY,NX)=0.0
      TNXFLB(L,NY,NX)=0.0
      TH1BFB(L,NY,NX)=0.0
      TH2BFB(L,NY,NX)=0.0
      TCOFHS(L,NY,NX)=0.0
      TCHFHS(L,NY,NX)=0.0
      TOXFHS(L,NY,NX)=0.0
      TNGFHS(L,NY,NX)=0.0
      TN2FHS(L,NY,NX)=0.0
      THGFHS(L,NY,NX)=0.0
      TN4FHW(L,NY,NX)=0.0
      TN3FHW(L,NY,NX)=0.0
      TNOFHW(L,NY,NX)=0.0
      TNXFHS(L,NY,NX)=0.0
      TH1PHS(L,NY,NX)=0.0
      TH2PHS(L,NY,NX)=0.0
      TN4FHB(L,NY,NX)=0.0
      TN3FHB(L,NY,NX)=0.0
      TNOFHB(L,NY,NX)=0.0
      TNXFHB(L,NY,NX)=0.0
      TH1BHB(L,NY,NX)=0.0
      TH2BHB(L,NY,NX)=0.0
!
!     ADD SOLUTE SINKS
!
!     R*SK2=total flux from nitro.f, uptake.f, solute.f
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band micropores
!     ZNH4B,ZNH3B,ZNO3B,ZNO2B,H1POB,H2POB=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band micropores
!
!
      ZNH4S2(L,NY,NX)=ZNH4S2(L,NY,NX)-RN4SK2(L,NY,NX)
      ZNH3S2(L,NY,NX)=ZNH3S2(L,NY,NX)-RN3SK2(L,NY,NX)
      ZNO3S2(L,NY,NX)=ZNO3S2(L,NY,NX)-RNOSK2(L,NY,NX)
      ZNO2S2(L,NY,NX)=ZNO2S2(L,NY,NX)-RNXSK2(L,NY,NX)
      H2PO42(L,NY,NX)=H2PO42(L,NY,NX)-RHPSK2(L,NY,NX)
      H1PO42(L,NY,NX)=H1PO42(L,NY,NX)-R1PSK2(L,NY,NX)
      ZNH4B2(L,NY,NX)=ZNH4B2(L,NY,NX)-R4BSK2(L,NY,NX)
      ZNH3B2(L,NY,NX)=ZNH3B2(L,NY,NX)-R3BSK2(L,NY,NX)
      ZNO3B2(L,NY,NX)=ZNO3B2(L,NY,NX)-RNBSK2(L,NY,NX)
      ZNO2B2(L,NY,NX)=ZNO2B2(L,NY,NX)-RNZSK2(L,NY,NX)
      H2POB2(L,NY,NX)=H2POB2(L,NY,NX)-RHBSK2(L,NY,NX)
      H1POB2(L,NY,NX)=H1POB2(L,NY,NX)-R1BSK2(L,NY,NX)
      ROXSK2(L,NY,NX)=ROXSK(M,L,NY,NX)*XNPT
    ENDIF
!
!     SOIL GAS FLUX ACCUMULATORS
!
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in micropores
!     ZN3G=gaseous NH3
!
    TCOFLG(L,NY,NX)=0.0
    TCHFLG(L,NY,NX)=0.0
    TOXFLG(L,NY,NX)=0.0
    TNGFLG(L,NY,NX)=0.0
    TN2FLG(L,NY,NX)=0.0
    TN3FLG(L,NY,NX)=0.0
    THGFLG(L,NY,NX)=0.0
    CO2S2(L,NY,NX)=CO2S2(L,NY,NX)-RCOSK2(L,NY,NX)
    CH4S2(L,NY,NX)=CH4S2(L,NY,NX)-RCHSK2(L,NY,NX)
    OXYS2(L,NY,NX)=OXYS2(L,NY,NX)-ROXSK2(L,NY,NX)
    Z2GS2(L,NY,NX)=Z2GS2(L,NY,NX)-RNGSK2(L,NY,NX)
    Z2OS2(L,NY,NX)=Z2OS2(L,NY,NX)-RN2SK2(L,NY,NX)
    H2GS2(L,NY,NX)=H2GS2(L,NY,NX)-RHGSK2(L,NY,NX)
    ZN3G2(L,NY,NX)=ZN3G2(L,NY,NX)-RNHSK2(L,NY,NX)
!     IF(I.EQ.105)THEN
!     WRITE(*,444)'CO2S1',I,J,NX,NY,L,M,MM,CO2S2(L,NY,NX)
!    2,RCOSK2(L,NY,NX),RCO2O(L,NY,NX),TCO2S(L,NY,NX)
!    3,TRCO2(L,NY,NX)
!     ENDIF
  ENDDO
  end subroutine ResetandInitFluxAccumulators
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSnowpack(M,NY,NX)

  implicit none

  integer, intent(in) :: M, NY, NX
  integer :: ICHKL,L,L2
  real(r8) :: VFLWS,VFLWW,VFLWR
  real(r8) :: VFLWNOB,VFLWNO3,VFLWNHB
  real(r8) :: VFLWNH4,VFLWPO4,VFLWPOB
!
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     VOLWSL=snowpack water content
!     FLQWM=snowpack water flux
!     R*BLS=solute flux in snowpack
!     X*BLS=hourly solute flux in snowpack
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
  ICHKL=0
  DO  L=1,JS
    IF(VHCPWM(M,L,NY,NX).GT.VHCPWX(NY,NX))THEN
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWM(M,L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AMAX1(0.0,AMIN1(1.0,FLQWM(M,L2,NY,NX)/VOLWSL(L,NY,NX)))
        ELSE
          VFLWW=1.0
        ENDIF
        RCOBLS(L2,NY,NX)=CO2W2(L,NY,NX)*VFLWW
        RCHBLS(L2,NY,NX)=CH4W2(L,NY,NX)*VFLWW
        ROXBLS(L2,NY,NX)=OXYW2(L,NY,NX)*VFLWW
        RNGBLS(L2,NY,NX)=ZNGW2(L,NY,NX)*VFLWW
        RN2BLS(L2,NY,NX)=ZN2W2(L,NY,NX)*VFLWW
        RN4BLW(L2,NY,NX)=ZN4W2(L,NY,NX)*VFLWW
        RN3BLW(L2,NY,NX)=ZN3W2(L,NY,NX)*VFLWW
        RNOBLW(L2,NY,NX)=ZNOW2(L,NY,NX)*VFLWW
        RH1PBS(L2,NY,NX)=Z1PW2(L,NY,NX)*VFLWW
        RH2PBS(L2,NY,NX)=ZHPW2(L,NY,NX)*VFLWW
        XCOBLS(L2,NY,NX)=XCOBLS(L2,NY,NX)+RCOBLS(L2,NY,NX)
        XCHBLS(L2,NY,NX)=XCHBLS(L2,NY,NX)+RCHBLS(L2,NY,NX)
        XOXBLS(L2,NY,NX)=XOXBLS(L2,NY,NX)+ROXBLS(L2,NY,NX)
        XNGBLS(L2,NY,NX)=XNGBLS(L2,NY,NX)+RNGBLS(L2,NY,NX)
        XN2BLS(L2,NY,NX)=XN2BLS(L2,NY,NX)+RN2BLS(L2,NY,NX)
        XN4BLW(L2,NY,NX)=XN4BLW(L2,NY,NX)+RN4BLW(L2,NY,NX)
        XN3BLW(L2,NY,NX)=XN3BLW(L2,NY,NX)+RN3BLW(L2,NY,NX)
        XNOBLW(L2,NY,NX)=XNOBLW(L2,NY,NX)+RNOBLW(L2,NY,NX)
        XH1PBS(L2,NY,NX)=XH1PBS(L2,NY,NX)+RH1PBS(L2,NY,NX)
        XH2PBS(L2,NY,NX)=XH2PBS(L2,NY,NX)+RH2PBS(L2,NY,NX)
      ELSE
        IF(L.LT.JS)THEN
          RCOBLS(L2,NY,NX)=0.0
          RCHBLS(L2,NY,NX)=0.0
          ROXBLS(L2,NY,NX)=0.0
          RNGBLS(L2,NY,NX)=0.0
          RN2BLS(L2,NY,NX)=0.0
          RN4BLW(L2,NY,NX)=0.0
          RN3BLW(L2,NY,NX)=0.0
          RNOBLW(L2,NY,NX)=0.0
          RH1PBS(L2,NY,NX)=0.0
          RH2PBS(L2,NY,NX)=0.0
        ENDIF
!
!     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
!
!     VOLWSL=snowpack water content
!     FLQRM,FLQSM,FLQHM=total water flux to litter,soil micropore,macropore
!     CVRD,BARE=litter cover fraction,1-CVRD
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     R*S0=solute flux to surface litter
!     R*S1,R*B1=solute flux to soil surface non-band,band
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
        IF(ICHKL.EQ.0)THEN
          IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AMAX1(0.0,AMIN1(1.0,FLQRM(M,NY,NX)/VOLWSL(L,NY,NX)))
            VFLWS=AMAX1(0.0,AMIN1(1.0,(FLQSM(M,NY,NX)+FLQHM(M,NY,NX))/VOLWSL(L,NY,NX)))
          ELSE
            VFLWR=CVRD(NY,NX)
            VFLWS=BARE(NY,NX)
          ENDIF
          VFLWNH4=VFLWS*VLNH4(NU(NY,NX),NY,NX)
          VFLWNHB=VFLWS*VLNHB(NU(NY,NX),NY,NX)
          VFLWNO3=VFLWS*VLNO3(NU(NY,NX),NY,NX)
          VFLWNOB=VFLWS*VLNOB(NU(NY,NX),NY,NX)
          VFLWPO4=VFLWS*VLPO4(NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*VLPOB(NU(NY,NX),NY,NX)
          RCOFLS0=CO2W2(L,NY,NX)*VFLWR
          RCHFLS0=CH4W2(L,NY,NX)*VFLWR
          ROXFLS0=OXYW2(L,NY,NX)*VFLWR
          RNGFLS0=ZNGW2(L,NY,NX)*VFLWR
          RN2FLS0=ZN2W2(L,NY,NX)*VFLWR
          RN4FLW0=ZN4W2(L,NY,NX)*VFLWR
          RN3FLW0=ZN3W2(L,NY,NX)*VFLWR
          RNOFLW0=ZNOW2(L,NY,NX)*VFLWR
          RH1PFS0=Z1PW2(L,NY,NX)*VFLWR
          RH2PFS0=ZHPW2(L,NY,NX)*VFLWR
          RCOFLS1=CO2W2(L,NY,NX)*VFLWS
          RCHFLS1=CH4W2(L,NY,NX)*VFLWS
          ROXFLS1=OXYW2(L,NY,NX)*VFLWS
          RNGFLS1=ZNGW2(L,NY,NX)*VFLWS
          RN2FLS1=ZN2W2(L,NY,NX)*VFLWS
          RN4FLW1=ZN4W2(L,NY,NX)*VFLWNH4
          RN3FLW1=ZN3W2(L,NY,NX)*VFLWNH4
          RNOFLW1=ZNOW2(L,NY,NX)*VFLWNO3
          RH1PFS1=Z1PW2(L,NY,NX)*VFLWPO4
          RH2PFS1=ZHPW2(L,NY,NX)*VFLWPO4
          RN4FLB1=ZN4W2(L,NY,NX)*VFLWNHB
          RN3FLB1=ZN3W2(L,NY,NX)*VFLWNHB
          RNOFLB1=ZNOW2(L,NY,NX)*VFLWNOB
          RH1BFB1=Z1PW2(L,NY,NX)*VFLWPOB
          RH2BFB1=ZHPW2(L,NY,NX)*VFLWPOB
          ICHKL=1
        ENDIF
      ENDIF
    ELSE
      RCOFLS0=0.0
      RCHFLS0=0.0
      ROXFLS0=0.0
      RNGFLS0=0.0
      RN2FLS0=0.0
      RN4FLW0=0.0
      RN3FLW0=0.0
      RNOFLW0=0.0
      RH1PFS0=0.0
      RH2PFS0=0.0
      RCOFLS1=0.0
      RCHFLS1=0.0
      ROXFLS1=0.0
      RNGFLS1=0.0
      RN2FLS1=0.0
      RN4FLW1=0.0
      RN3FLW1=0.0
      RNOFLW1=0.0
      RH1PFS1=0.0
      RH2PFS1=0.0
      RN4FLB1=0.0
      RN3FLB1=0.0
      RNOFLB1=0.0
      RH1BFB1=0.0
      RH2BFB1=0.0
    ENDIF
  ENDDO
  end subroutine SoluteFluxSnowpack
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSurface(M,NY,NX,FLQM)
  implicit none
  integer, intent(in) :: NY, NX,M
  real(r8),intent(out) :: FLQM(3,JD,JV,JH)
  integer :: K
  real(r8) :: VOLWPA,VOLWPB,VOLWOB
  real(r8) :: VOLWOA
  real(r8) :: TORT0,TORT1
  real(r8) :: DLYR0,DLYR1
  real(r8) :: CCO2GQ,CCH4GQ
  real(r8) :: COXYGQ,CZ2GGQ,CZ2OGQ,CZN3GQ,CH2GGQ,CZN3BQ
  real(r8) :: DFGSCH,DFGSCO,DFGSOX,DFGSNG,DFGSN2,DFGSN3,DFGSHL

!     VOLWM,VOLWHM,VOLPM,FLPM=micropore,macropore water volume, air volume and change in air volume
!     FLWM,FLWHM=water flux into soil micropore,macropore from watsub.f
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     FLVM,FLM=air,water flux in gas flux calculations
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
  VOLWMA(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*VLNH4(NU(NY,NX),NY,NX)
  VOLWMB(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX)*VLNHB(NU(NY,NX),NY,NX)
  VOLWXA(NU(NY,NX),NY,NX)=14.0*VOLWMA(NU(NY,NX),NY,NX)
  VOLWXB(NU(NY,NX),NY,NX)=14.0*VOLWMB(NU(NY,NX),NY,NX)
  VOLWOA=VOLWM(M,NU(NY,NX),NY,NX)*VLNO3(NU(NY,NX),NY,NX)
  VOLWOB=VOLWM(M,NU(NY,NX),NY,NX)*VLNOB(NU(NY,NX),NY,NX)
  VOLWPA=VOLWM(M,NU(NY,NX),NY,NX)*VLPO4(NU(NY,NX),NY,NX)
  VOLWPB=VOLWM(M,NU(NY,NX),NY,NX)*VLPOB(NU(NY,NX),NY,NX)
  VOLPMA(NU(NY,NX),NY,NX)=VOLPM(M,NU(NY,NX),NY,NX)*VLNH4(NU(NY,NX),NY,NX)
  VOLPMB(NU(NY,NX),NY,NX)=VOLPM(M,NU(NY,NX),NY,NX)*VLNHB(NU(NY,NX),NY,NX)
  FLVM(NU(NY,NX),NY,NX)=FLPM(M,NU(NY,NX),NY,NX)*XNPT
  FLQM(3,NU(NY,NX),NY,NX)=(FLWM(M,3,NU(NY,NX),NY,NX)+FLWHM(M,3,NU(NY,NX),NY,NX))*XNPT
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE RESIDUE
!
!     VOLT,DLYR,AREA=litter volume, thickness, area
!     VOLWM=micropore water-filled porosity from watsub.f
!     TORT=tortuosity from hour1.f
!     *SGL*=solute diffusivity
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     DFGS*=effective solute diffusivity
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in litter
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in litter
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in litter
!     C*1=solute concentration in litter
!
  IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX))
    TORT0=TORT(M,0,NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5*DLYR0)*CVRD(NY,NX)
    DFGSCO=CLSGL2(0,NY,NX)*TORT0
    DFGSCH=CQSGL2(0,NY,NX)*TORT0
    DFGSOX=OLSGL2(0,NY,NX)*TORT0
    DFGSNG=ZLSGL2(0,NY,NX)*TORT0
    DFGSN2=ZNSGL2(0,NY,NX)*TORT0
    DFGSN3=ZVSGL2(0,NY,NX)*TORT0
    DFGSHL=HLSGL2(0,NY,NX)*TORT0
    DO  K=0,jcplx1
      COQC1(K)=AMAX1(0.0,OQC2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
      COQN1(K)=AMAX1(0.0,OQN2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
      COQP1(K)=AMAX1(0.0,OQP2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
      COQA1(K)=AMAX1(0.0,OQA2(K,0,NY,NX)/VOLWM(M,0,NY,NX))
    ENDDO
    CCO2S1=AMAX1(0.0,CO2S2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CCH4S1=AMAX1(0.0,CH4S2(0,NY,NX)/VOLWM(M,0,NY,NX))
    COXYS1=AMAX1(0.0,OXYS2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CZ2GS1=AMAX1(0.0,Z2GS2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CZ2OS1=AMAX1(0.0,Z2OS2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CH2GS1=AMAX1(0.0,H2GS2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CNH4S1=AMAX1(0.0,ZNH4S2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CNH3S1=AMAX1(0.0,ZNH3S2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CNO3S1=AMAX1(0.0,ZNO3S2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CNO2S1=AMAX1(0.0,ZNO2S2(0,NY,NX)/VOLWM(M,0,NY,NX))
    CP14S1=AMAX1(0.0,H1PO42(0,NY,NX)/VOLWM(M,0,NY,NX))
    CPO4S1=AMAX1(0.0,H2PO42(0,NY,NX)/VOLWM(M,0,NY,NX))
!
!     EQUILIBRIUM CONCENTRATIONS AT RESIDUE SURFACE AT WHICH
!     AQUEOUS DIFFUSION THROUGH RESIDUE SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM AQUEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
!     PARR=boundary layer conductance above litter surface from watsub.f
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at litter surface
!     S*L=solubility of gas in water from hour1.f
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     DFGS*=effective solute diffusivity
!
    CCO2GQ=(PARR(NY,NX)*CCO2E(NY,NX)*SCO2L(0,NY,NX)+DFGSCO*CCO2S1) &
      /(DFGSCO+PARR(NY,NX))
    CCH4GQ=(PARR(NY,NX)*CCH4E(NY,NX)*SCH4L(0,NY,NX)+DFGSCH*CCH4S1) &
      /(DFGSCH+PARR(NY,NX))
    COXYGQ=(PARR(NY,NX)*COXYE(NY,NX)*SOXYL(0,NY,NX)+DFGSOX*COXYS1) &
      /(DFGSOX+PARR(NY,NX))
    CZ2GGQ=(PARR(NY,NX)*CZ2GE(NY,NX)*SN2GL(0,NY,NX)+DFGSNG*CZ2GS1) &
      /(DFGSNG+PARR(NY,NX))
    CZ2OGQ=(PARR(NY,NX)*CZ2OE(NY,NX)*SN2OL(0,NY,NX)+DFGSN2*CZ2OS1) &
      /(DFGSN2+PARR(NY,NX))
    CZN3GQ=(PARR(NY,NX)*CNH3E(NY,NX)*SNH3L(0,NY,NX)+DFGSN3*CNH3S1) &
      /(DFGSN3+PARR(NY,NX))
    CH2GGQ=(PARR(NY,NX)*CH2GE(NY,NX)*SH2GL(0,NY,NX)+DFGSHL*CH2GS1) &
      /(DFGSHL+PARR(NY,NX))
!
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND RESIDUE SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!
!     R*DFR,X*DFR=current,hourly gas exchange between atmosphere and surface litter water
!     gas code: *CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at litter surface
!     VOLWM=litter water volume
!     DFGS*=effective solute diffusivity
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     R*DXR=R*DFR for gas flux calculations
!
    RCODFR(NY,NX)=(CCO2GQ-CCO2S1)*AMIN1(VOLWM(M,0,NY,NX),DFGSCO)
    RCHDFR(NY,NX)=(CCH4GQ-CCH4S1)*AMIN1(VOLWM(M,0,NY,NX),DFGSCH)
    ROXDFR(NY,NX)=(COXYGQ-COXYS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSOX)
    RNGDFR(NY,NX)=(CZ2GGQ-CZ2GS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSNG)
    RN2DFR(NY,NX)=(CZ2OGQ-CZ2OS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSN2)
    RN3DFR(NY,NX)=(CZN3GQ-CNH3S1)*AMIN1(VOLWM(M,0,NY,NX),DFGSN3)
    RHGDFR(NY,NX)=(CH2GGQ-CH2GS1)*AMIN1(VOLWM(M,0,NY,NX),DFGSHL)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
    XCODFR(NY,NX)=XCODFR(NY,NX)+RCODFR(NY,NX)
    XCHDFR(NY,NX)=XCHDFR(NY,NX)+RCHDFR(NY,NX)
    XOXDFR(NY,NX)=XOXDFR(NY,NX)+ROXDFR(NY,NX)
    XNGDFR(NY,NX)=XNGDFR(NY,NX)+RNGDFR(NY,NX)
    XN2DFR(NY,NX)=XN2DFR(NY,NX)+RN2DFR(NY,NX)
    XN3DFR(NY,NX)=XN3DFR(NY,NX)+RN3DFR(NY,NX)
    XHGDFR(NY,NX)=XHGDFR(NY,NX)+RHGDFR(NY,NX)
!     IF(J.EQ.24)THEN
!     WRITE(*,1118)'RCODFR',I,J,NX,NY,M,MM,NU(NY,NX),RCODFR(NY,NX)
!    2,PARR(NY,NX),CO2S2(0,NY,NX),CCO2E(NY,NX),CCO2Q,DFGSCO
!    3,DCO21,CCO22,SCO2L(0,NY,NX),TORT1
!    4,DLYR(3,0,NY,NX),VOLWM(M,0,NY,NX)
!     WRITE(*,1118)'RCHDFR',I,J,NX,NY,M,MM,NU(NY,NX),RCHDFR(NY,NX)
!    2,CH4GQ,CH4S2(0,NY,NX),PARR(NY,NX),CCH4E(NY,NX)
!    3,VOLWCH(0,NY,NX),DFGSCH,TORT(M,0,NY,NX)
!     WRITE(*,1118)'ROXDFR',I,J,NX,NY,M,MM,NU(NY,NX),ROXDFR(NY,NX)
!    2,XOXDFR(NY,NX),COXYGQ,COXYS1,PARR(NY,NX),DFGSOX,COXYE(NY,NX)
!    3,TORT(M,0,NY,NX),VOLWM(M,0,NY,NX),VOLT(0,NY,NX),DLYR(3,0,NY,NX)
!1118  FORMAT(A8,7I4,30E12.4)
!     ENDIF
  ELSE
    RCODFR(NY,NX)=0.0
    RCHDFR(NY,NX)=0.0
    ROXDFR(NY,NX)=0.0
    RNGDFR(NY,NX)=0.0
    RN2DFR(NY,NX)=0.0
    RN3DFR(NY,NX)=0.0
    RHGDFR(NY,NX)=0.0
  ENDIF
  RCODXR=RCODFR(NY,NX)*XNPT
  RCHDXR=RCHDFR(NY,NX)*XNPT
  ROXDXR=ROXDFR(NY,NX)*XNPT
  RNGDXR=RNGDFR(NY,NX)*XNPT
  RN2DXR=RN2DFR(NY,NX)*XNPT
  RN3DXR=RN3DFR(NY,NX)*XNPT
  RHGDXR=RHGDFR(NY,NX)*XNPT
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE SOIL LAYER
!
!     VOLWM=micropore water-filled porosity from watsub.f
!     TORT=tortuosity from hour1.f
!     *SGL*=solute diffusivity
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     DFGS*=effective solute diffusivity
!     C*2=solute concentration in soil surface
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in soil
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in litter
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band micropores
!     ZNH4B,ZNH3B,ZNO3B,ZNO2B,H1POB,H2POB=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band micropores
!
  IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
    TORT1=TORT(M,NU(NY,NX),NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5*DLYR1)
    DFGSCO=CLSGL2(NU(NY,NX),NY,NX)*TORT1
    DFGSCH=CQSGL2(NU(NY,NX),NY,NX)*TORT1
    DFGSOX=OLSGL2(NU(NY,NX),NY,NX)*TORT1
    DFGSNG=ZLSGL2(NU(NY,NX),NY,NX)*TORT1
    DFGSN2=ZNSGL2(NU(NY,NX),NY,NX)*TORT1
    DFGSN3=ZVSGL2(NU(NY,NX),NY,NX)*TORT1
    DFGSHL=HLSGL2(NU(NY,NX),NY,NX)*TORT1
    DO  K=0,jcplx1
      COQC2(K)=AMAX1(0.0,OQC2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
      COQN2(K)=AMAX1(0.0,OQN2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
      COQP2(K)=AMAX1(0.0,OQP2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
      COQA2(K)=AMAX1(0.0,OQA2(K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    ENDDO
    CCO2S2=AMAX1(0.0,CO2S2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CCH4S2=AMAX1(0.0,CH4S2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    COXYS2=AMAX1(0.0,OXYS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CZ2GS2=AMAX1(0.0,Z2GS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CZ2OS2=AMAX1(0.0,Z2OS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    CH2GS2=AMAX1(0.0,H2GS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
    IF(VOLWMA(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      CNH3S2=AMAX1(0.0,ZNH3S2(NU(NY,NX),NY,NX)/VOLWMA(NU(NY,NX),NY,NX))
      CNH4S2=AMAX1(0.0,ZNH4S2(NU(NY,NX),NY,NX)/VOLWMA(NU(NY,NX),NY,NX))
    ELSE
      CNH3S2=0.0
      CNH4S2=0.0
    ENDIF
    IF(VOLWOA.GT.ZEROS2(NY,NX))THEN
      CNO3S2=AMAX1(0.0,ZNO3S2(NU(NY,NX),NY,NX)/VOLWOA)
      CNO2S2=AMAX1(0.0,ZNO2S2(NU(NY,NX),NY,NX)/VOLWOA)
    ELSE
      CNO3S2=0.0
      CNO2S2=0.0
    ENDIF
    IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
      CP14S2=AMAX1(0.0,H1PO42(NU(NY,NX),NY,NX)/VOLWPA)
      CPO4S2=AMAX1(0.0,H2PO42(NU(NY,NX),NY,NX)/VOLWPA)
    ELSE
      CP14S2=0.0
      CPO4S2=0.0
    ENDIF
    IF(VOLWMB(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      CNH3B2=AMAX1(0.0,ZNH3B2(NU(NY,NX),NY,NX)/VOLWMB(NU(NY,NX),NY,NX))
      CNH4B2=AMAX1(0.0,ZNH4B2(NU(NY,NX),NY,NX)/VOLWMB(NU(NY,NX),NY,NX))
    ELSE
      CNH3B2=CNH3S2
      CNH4B2=CNH4S2
    ENDIF
    IF(VOLWOB.GT.ZEROS2(NY,NX))THEN
      CNO3B2=AMAX1(0.0,ZNO3B2(NU(NY,NX),NY,NX)/VOLWOB)
      CNO2B2=AMAX1(0.0,ZNO2B2(NU(NY,NX),NY,NX)/VOLWOB)
    ELSE
      CNO3B2=CNO3S2
      CNO2B2=CNO2S2
    ENDIF
    IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
      CP14B2=AMAX1(0.0,H1POB2(NU(NY,NX),NY,NX)/VOLWPB)
      CPO4B2=AMAX1(0.0,H2POB2(NU(NY,NX),NY,NX)/VOLWPB)
    ELSE
      CP14B2=CP14S2
      CPO4B2=CPO4S2
    ENDIF
!
!     EQUILIBRIUM CONCENTRATIONS AT SOIL SURFACE AT WHICH
!     AQUEOUS DIFFUSION THROUGH SOIL SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM AQUEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
!     PARG=boundary layer conductance above soil surface from watsub.f
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at soil surface
!     S*L=solubility of gas in water from hour1.f
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     DFGS*=effective solute diffusivity
!
    CCO2GQ=(PARG(M,NY,NX)*CCO2E(NY,NX)*SCO2L(NU(NY,NX),NY,NX) &
      +DFGSCO*CCO2S2)/(DFGSCO+PARG(M,NY,NX))
    CCH4GQ=(PARG(M,NY,NX)*CCH4E(NY,NX)*SCH4L(NU(NY,NX),NY,NX) &
      +DFGSCH*CCH4S2)/(DFGSCH+PARG(M,NY,NX))
    COXYGQ=(PARG(M,NY,NX)*COXYE(NY,NX)*SOXYL(NU(NY,NX),NY,NX) &
      +DFGSOX*COXYS2)/(DFGSOX+PARG(M,NY,NX))
    CZ2GGQ=(PARG(M,NY,NX)*CZ2GE(NY,NX)*SN2GL(NU(NY,NX),NY,NX) &
      +DFGSNG*CZ2GS2)/(DFGSNG+PARG(M,NY,NX))
    CZ2OGQ=(PARG(M,NY,NX)*CZ2OE(NY,NX)*SN2OL(NU(NY,NX),NY,NX) &
      +DFGSN2*CZ2OS2)/(DFGSN2+PARG(M,NY,NX))
    CZN3GQ=(PARG(M,NY,NX)*CNH3E(NY,NX)*SNH3L(NU(NY,NX),NY,NX) &
      +DFGSN3*CNH3S2)/(DFGSN3+PARG(M,NY,NX))
    CZN3BQ=(PARG(M,NY,NX)*CNH3E(NY,NX)*SNH3L(NU(NY,NX),NY,NX) &
      +DFGSN3*CNH3B2)/(DFGSN3+PARG(M,NY,NX))
    CH2GGQ=(PARG(M,NY,NX)*CH2GE(NY,NX)*SH2GL(NU(NY,NX),NY,NX) &
      +DFGSHL*CH2GS2)/(DFGSHL+PARG(M,NY,NX))
!
!     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
!     BETWEEN ATMOSPHERIC AND SOIL SURFACE EQUILIBRIUM
!     CONCENTRATIONS
!
!     R*DFS,X*DFS=current,cumulative gas exchange between atmosphere and soil surface water
!     gas code:*CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     C*E=atmospheric gas concentration from hour1.f
!     C*Q=equilibrium gas concentration at soil surface
!     VOLWM=litter water volume
!     DFGS*=effective solute diffusivity
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     R*DXS=R*DFS for gas flux calculations
!
    RCODFS(NY,NX)=(CCO2GQ-CCO2S2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGSCO)
    RCHDFS(NY,NX)=(CCH4GQ-CCH4S2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGSCH)
    ROXDFS(NY,NX)=(COXYGQ-COXYS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGSOX)
    RNGDFS(NY,NX)=(CZ2GGQ-CZ2GS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGSNG)
    RN2DFS(NY,NX)=(CZ2OGQ-CZ2OS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGSN2)
    RN3DFS(NY,NX)=(CZN3GQ-CNH3S2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX)*VLNH4(NU(NY,NX),NY,NX),DFGSN3)
    RNBDFS(NY,NX)=(CZN3BQ-CNH3B2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX)*VLNHB(NU(NY,NX),NY,NX),DFGSN3)
    RHGDFS(NY,NX)=(CH2GGQ-CH2GS2)*AMIN1(VOLWM(M,NU(NY,NX),NY,NX),DFGSHL)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
    XCODFS(NY,NX)=XCODFS(NY,NX)+RCODFS(NY,NX)
    XCHDFS(NY,NX)=XCHDFS(NY,NX)+RCHDFS(NY,NX)
    XOXDFS(NY,NX)=XOXDFS(NY,NX)+ROXDFS(NY,NX)
    XNGDFS(NY,NX)=XNGDFS(NY,NX)+RNGDFS(NY,NX)
    XN2DFS(NY,NX)=XN2DFS(NY,NX)+RN2DFS(NY,NX)
    XN3DFS(NY,NX)=XN3DFS(NY,NX)+RN3DFS(NY,NX)
    XNBDFS(NY,NX)=XNBDFS(NY,NX)+RNBDFS(NY,NX)
    XHGDFS(NY,NX)=XHGDFS(NY,NX)+RHGDFS(NY,NX)
!     IF(J.EQ.24)THEN
!     WRITE(*,1118)'RCODFS',I,J,NX,NY,M,MM,NU(NY,NX),RCODFS(NY,NX)
!    2,XCODFS(NY,NX),CO2GQ,CO2S2X,CO2S2(NU(NY,NX),NY,NX),PARG(M,NY,NX)
!    3,CCO2E(NY,NX),VOLWCO(NU(NY,NX),NY,NX),DFGSCO
!    2,TORT(M,NU(NY,NX),NY,NX),CLSGL2(NU(NY,NX),NY,NX),TORT1
!    4,DLYR(3,NU(NY,NX),NY,NX),CO2S(NU(NY,NX),NY,NX)
!    5,RCOSK2(NU(NY,NX),NY,NX)
!     WRITE(*,1118)'RCHDFS',I,J,NX,NY,M,MM,NU(NY,NX),RCHDFS(NY,NX)
!    2,XCHDFS(NY,NX),CH4GQ,CH4S2X,CH4S2(NU(NY,NX),NY,NX),PARG(M,NY,NX)
!    3,CCH4E(NY,NX),VOLWCH(NU(NY,NX),NY,NX),DFGSCH,TORT(M,0,NY,NX)
!     WRITE(*,1118)'ROXDFS',I,J,NX,NY,M,MM,NU(NY,NX),ROXDFS(NY,NX)
!    2,XOXDFS(NY,NX),COXYGQ,COXYS2,OXYS2(NU(NY,NX),NY,NX),PARG(M,NY,NX)
!    3,COXYE(NY,NX),VOLWM(M,NU(NY,NX),NY,NX),DFGSOX,TORT(M,0,NY,NX)
!     WRITE(*,1118)'RNGDFS',I,J,NX,NY,M,MM,NU(NY,NX),RNGDFS(NY,NX)
!    2,XNGDFS(NY,NX),Z2GGQ,Z2GS2X,Z2GS2(NU(NY,NX),NY,NX),PARG(M,NY,NX)
!    3,CZ2GE(NY,NX),VOLWNG(NU(NY,NX),NY,NX),DFGSNG,TORT(M,0,NY,NX)
!     WRITE(*,1118)'RN2DFS',I,J,NX,NY,M,MM,NU(NY,NX),RN2DFS(NY,NX)
!    2,XN2DFS(NY,NX),Z2OGQ,Z2OS2X,Z2OS2(NU(NY,NX),NY,NX),PARG(M,NY,NX)
!    3,CZ2OE(NY,NX),VOLWN2(NU(NY,NX),NY,NX),DFGSN2,TORT(M,0,NY,NX)
!    4,VOLWM(M,NU(NY,NX),NY,NX),SN2OL(NU(NY,NX),NY,NX)
!    5,TKS(NU(NY,NX),NY,NX)
!     WRITE(*,1118)'RN3DFS',I,J,NX,NY,M,MM,NU(NY,NX),RN3DFS(NY,NX)
!    2,XN3DFS(NY,NX),ZN3GQ,ZN3S2X,PARG(M,NY,NX),CNH3E(NY,NX)
!    3,VOLWN3(NU(NY,NX),NY,NX),DFGSN3,ZNH3S2(NU(NY,NX),NY,NX)
!    4,ZN3G2(NU(NY,NX),NY,NX)
!     ENDIF
  ELSE
    RCODFS(NY,NX)=0.0
    RCHDFS(NY,NX)=0.0
    ROXDFS(NY,NX)=0.0
    RNGDFS(NY,NX)=0.0
    RN2DFS(NY,NX)=0.0
    RN3DFS(NY,NX)=0.0
    RNBDFS(NY,NX)=0.0
    RHGDFS(NY,NX)=0.0
  ENDIF
  RCODXS=RCODFS(NY,NX)*XNPT
  RCHDXS=RCHDFS(NY,NX)*XNPT
  ROXDXS=ROXDFS(NY,NX)*XNPT
  RNGDXS=RNGDFS(NY,NX)*XNPT
  RN2DXS=RN2DFS(NY,NX)*XNPT
  RN3DXS=RN3DFS(NY,NX)*XNPT
  RNBDXS=RNBDFS(NY,NX)*XNPT
  RHGDXS=RHGDFS(NY,NX)*XNPT
  end subroutine SoluteFluxSurface
!------------------------------------------------------------------------------------------

  subroutine ConvectiveSurfaceSoluteFlux(M,NY,NX,FLWRM1)
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8),intent(out) :: FLWRM1
  REAL(R8) :: VFLW
  integer :: K

  FLWRM1=FLWRM(M,NY,NX)
!
!     FLWRM=litter-soil water flux from watsub.f
!
!     IF WATER FLUX FROM 'WATSUB' IS FROM RESIDUE TO
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN RESIDUE
!
!     VOLWM=litter water volume
!     RFL*=soil-litter convective solute flux
!     *S2=litter solute content
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
  IF(FLWRM1.GT.0.0)THEN
    IF(VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AMAX1(0.0,AMIN1(VFLWX,FLWRM1/VOLWM(M,0,NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQC2(K,0,NY,NX))
      RFLON(K)=VFLW*AMAX1(0.0,OQN2(K,0,NY,NX))
      RFLOP(K)=VFLW*AMAX1(0.0,OQP2(K,0,NY,NX))
      RFLOA(K)=VFLW*AMAX1(0.0,OQA2(K,0,NY,NX))
    ENDDO
    RFLCOS=VFLW*AMAX1(0.0,CO2S2(0,NY,NX))
    RFLCHS=VFLW*AMAX1(0.0,CH4S2(0,NY,NX))
    RFLOXS=VFLW*AMAX1(0.0,OXYS2(0,NY,NX))
    RFLNGS=VFLW*AMAX1(0.0,Z2GS2(0,NY,NX))
    RFLN2S=VFLW*AMAX1(0.0,Z2OS2(0,NY,NX))
    RFLHGS=VFLW*AMAX1(0.0,H2GS2(0,NY,NX))
    RFLNH4=VFLW*AMAX1(0.0,ZNH4S2(0,NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNH3=VFLW*AMAX1(0.0,ZNH3S2(0,NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNO3=VFLW*AMAX1(0.0,ZNO3S2(0,NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLNO2=VFLW*AMAX1(0.0,ZNO2S2(0,NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLP14=VFLW*AMAX1(0.0,H1PO42(0,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLPO4=VFLW*AMAX1(0.0,H2PO42(0,NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLN4B=VFLW*AMAX1(0.0,ZNH4S2(0,NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLN3B=VFLW*AMAX1(0.0,ZNH3S2(0,NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLNOB=VFLW*AMAX1(0.0,ZNO3S2(0,NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLN2B=VFLW*AMAX1(0.0,ZNO2S2(0,NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLP1B=VFLW*AMAX1(0.0,H1PO42(0,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
    RFLPOB=VFLW*AMAX1(0.0,H2PO42(0,NY,NX))*VLPOB(NU(NY,NX),NY,NX)
!
!     IF WATER FLUX FROM 'WATSUB' IS TO RESIDUE FROM
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE
!
!     VOLWM=micropore water-filled porosity from watsub.f
!     RFL*=soil-litter convective solute flux
!     *S2=soil solute content
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
  ELSE
    IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AMIN1(0.0,AMAX1(-VFLWX,FLWRM1/VOLWM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQC2(K,NU(NY,NX),NY,NX))
      RFLON(K)=VFLW*AMAX1(0.0,OQN2(K,NU(NY,NX),NY,NX))
      RFLOP(K)=VFLW*AMAX1(0.0,OQP2(K,NU(NY,NX),NY,NX))
      RFLOA(K)=VFLW*AMAX1(0.0,OQA2(K,NU(NY,NX),NY,NX))
    ENDDO
    RFLCOS=VFLW*AMAX1(0.0,CO2S2(NU(NY,NX),NY,NX))
    RFLCHS=VFLW*AMAX1(0.0,CH4S2(NU(NY,NX),NY,NX))
    RFLOXS=VFLW*AMAX1(0.0,OXYS2(NU(NY,NX),NY,NX))
    RFLNGS=VFLW*AMAX1(0.0,Z2GS2(NU(NY,NX),NY,NX))
    RFLN2S=VFLW*AMAX1(0.0,Z2OS2(NU(NY,NX),NY,NX))
    RFLHGS=VFLW*AMAX1(0.0,H2GS2(NU(NY,NX),NY,NX))
    RFLNH4=VFLW*AMAX1(0.0,ZNH4S2(NU(NY,NX),NY,NX))
    RFLNH3=VFLW*AMAX1(0.0,ZNH3S2(NU(NY,NX),NY,NX))
    RFLNO3=VFLW*AMAX1(0.0,ZNO3S2(NU(NY,NX),NY,NX))
    RFLNO2=VFLW*AMAX1(0.0,ZNO2S2(NU(NY,NX),NY,NX))
    RFLP14=VFLW*AMAX1(0.0,H1PO42(NU(NY,NX),NY,NX))
    RFLPO4=VFLW*AMAX1(0.0,H2PO42(NU(NY,NX),NY,NX))
    RFLN4B=VFLW*AMAX1(0.0,ZNH4B2(NU(NY,NX),NY,NX))
    RFLN3B=VFLW*AMAX1(0.0,ZNH3B2(NU(NY,NX),NY,NX))
    RFLNOB=VFLW*AMAX1(0.0,ZNO3B2(NU(NY,NX),NY,NX))
    RFLN2B=VFLW*AMAX1(0.0,ZNO2B2(NU(NY,NX),NY,NX))
    RFLP1B=VFLW*AMAX1(0.0,H1POB2(NU(NY,NX),NY,NX))
    RFLPOB=VFLW*AMAX1(0.0,H2POB2(NU(NY,NX),NY,NX))
  ENDIF
  end subroutine ConvectiveSurfaceSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine DiffusiveFluxAtSoilSurface(M,NY,NX,FLWRM1)
  implicit none

  integer, intent(in) :: M, NY, NX
  real(r8), intent(in) :: FLWRM1
  real(r8) :: TORT0,TORT1
  real(r8) :: DLYR0,DLYR1
  real(r8) :: DIFCQ,DIFCS,DIFHG,DIFN2,DIFNG,DISPN
  real(r8) :: DIFOS,DIFOC,DIFON,DIFOP,DIFOA,DIFNH,DIFNO,DIFPO
  real(r8) :: DIFOC0,DIFON0,DIFOP0,DIFOA0,DIFNH0
  real(r8) :: DIFNO0,DIFPO0,DIFCS0,DIFCQ0,DIFOS0,DIFNG0,DIFN20
  real(r8) :: DIFPO1,DIFCS1,DIFCQ1,DIFOS1,DIFNG1,DIFN21,DIFHG1
  real(r8) :: DIFHG0,DIFOC1,DIFON1,DIFOP1,DIFOA1,DIFNH1,DIFNO1
  integer :: K
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VOLWM=micropore water-filled porosity from watsub.f
!
  IF((VOLT(0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
    .AND.(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN
!
!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE
!
!     DLYR0,DLYR1=litter, soil surface thickness
!     TORT=tortuosity from hour1.f
!     CVRD,BARE=litter cover fraction,1-CVRD
!     DISP=dispersivity parameter
!     FLWRM=litter-soil water flux from watsub.f
!     DIF*0,DIF*1=aqueous diffusivity-dispersivity in litter, soil surface
!     *SGL*=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!
    DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX))
    TORT0=TORT(M,0,NY,NX)/DLYR0*CVRD(NY,NX)
    DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
    TORT1=TORT(M,NU(NY,NX),NY,NX)/DLYR1
    DISPN=DISP(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))
    DIFOC0=(OCSGL2(0,NY,NX)*TORT0+DISPN)
    DIFON0=(ONSGL2(0,NY,NX)*TORT0+DISPN)
    DIFOP0=(OPSGL2(0,NY,NX)*TORT0+DISPN)
    DIFOA0=(OASGL2(0,NY,NX)*TORT0+DISPN)
    DIFNH0=(ZNSGL2(0,NY,NX)*TORT0+DISPN)
    DIFNO0=(ZOSGL2(0,NY,NX)*TORT0+DISPN)
    DIFPO0=(POSGL2(0,NY,NX)*TORT0+DISPN)
    DIFCS0=(CLSGL2(0,NY,NX)*TORT0+DISPN)
    DIFCQ0=(CQSGL2(0,NY,NX)*TORT0+DISPN)
    DIFOS0=(OLSGL2(0,NY,NX)*TORT0+DISPN)
    DIFNG0=(ZLSGL2(0,NY,NX)*TORT0+DISPN)
    DIFN20=(ZVSGL2(0,NY,NX)*TORT0+DISPN)
    DIFHG0=(HLSGL2(0,NY,NX)*TORT0+DISPN)
    DIFOC1=(OCSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFON1=(ONSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFOP1=(OPSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFOA1=(OASGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFNH1=(ZNSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFNO1=(ZOSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFPO1=(POSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFCS1=(CLSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFCQ1=(CQSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFOS1=(OLSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFNG1=(ZLSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFN21=(ZVSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFHG1=(HLSGL2(NU(NY,NX),NY,NX)*TORT1+DISPN)
    DIFOC=DIFOC0*DIFOC1/(DIFOC0+DIFOC1)*AREA(3,NU(NY,NX),NY,NX)
    DIFON=DIFON0*DIFON1/(DIFON0+DIFON1)*AREA(3,NU(NY,NX),NY,NX)
    DIFOP=DIFOP0*DIFOP1/(DIFOP0+DIFOP1)*AREA(3,NU(NY,NX),NY,NX)
    DIFOA=DIFOA0*DIFOA1/(DIFOA0+DIFOA1)*AREA(3,NU(NY,NX),NY,NX)
    DIFNH=DIFNH0*DIFNH1/(DIFNH0+DIFNH1)*AREA(3,NU(NY,NX),NY,NX)
    DIFNO=DIFNO0*DIFNO1/(DIFNO0+DIFNO1)*AREA(3,NU(NY,NX),NY,NX)
    DIFPO=DIFPO0*DIFPO1/(DIFPO0+DIFPO1)*AREA(3,NU(NY,NX),NY,NX)
    DIFCS=DIFCS0*DIFCS1/(DIFCS0+DIFCS1)*AREA(3,NU(NY,NX),NY,NX)
    DIFCQ=DIFCQ0*DIFCQ1/(DIFCQ0+DIFCQ1)*AREA(3,NU(NY,NX),NY,NX)
    DIFOS=DIFOS0*DIFOS1/(DIFOS0+DIFOS1)*AREA(3,NU(NY,NX),NY,NX)
    DIFNG=DIFNG0*DIFNG1/(DIFNG0+DIFNG1)*AREA(3,NU(NY,NX),NY,NX)
    DIFN2=DIFN20*DIFN21/(DIFN20+DIFN21)*AREA(3,NU(NY,NX),NY,NX)
    DIFHG=DIFHG0*DIFHG1/(DIFHG0+DIFHG1)*AREA(3,NU(NY,NX),NY,NX)
!     DFV*S,DFV*B=diffusive solute flux between litter and soil surface in non-band,band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     C*1,C*2=solute concentration in litter,soil surface
!
    DO  K=0,jcplx1
      DFVOC(K)=DIFOC*(COQC1(K)-COQC2(K))
      DFVON(K)=DIFON*(COQN1(K)-COQN2(K))
      DFVOP(K)=DIFOP*(COQP1(K)-COQP2(K))
      DFVOA(K)=DIFOA*(COQA1(K)-COQA2(K))
    ENDDO
    DFVCOS=DIFCS*(CCO2S1-CCO2S2)
    DFVCHS=DIFCQ*(CCH4S1-CCH4S2)
    DFVOXS=DIFOS*(COXYS1-COXYS2)
    DFVNGS=DIFNG*(CZ2GS1-CZ2GS2)
    DFVN2S=DIFN2*(CZ2OS1-CZ2OS2)
    DFVHGS=DIFHG*(CH2GS1-CH2GS2)
    DFVNH4=DIFNH*(CNH4S1-CNH4S2)*VLNH4(NU(NY,NX),NY,NX)
    DFVNH3=DIFNH*(CNH3S1-CNH3S2)*VLNH4(NU(NY,NX),NY,NX)
    DFVNO3=DIFNO*(CNO3S1-CNO3S2)*VLNO3(NU(NY,NX),NY,NX)
    DFVNO2=DIFNO*(CNO2S1-CNO2S2)*VLNO3(NU(NY,NX),NY,NX)
    DFVP14=DIFPO*(CP14S1-CP14S2)*VLPO4(NU(NY,NX),NY,NX)
    DFVPO4=DIFPO*(CPO4S1-CPO4S2)*VLPO4(NU(NY,NX),NY,NX)
    DFVN4B=DIFNH*(CNH4S1-CNH4B2)*VLNHB(NU(NY,NX),NY,NX)
    DFVN3B=DIFNH*(CNH3S1-CNH3B2)*VLNHB(NU(NY,NX),NY,NX)
    DFVNOB=DIFNO*(CNO3S1-CNO3B2)*VLNOB(NU(NY,NX),NY,NX)
    DFVN2B=DIFNO*(CNO2S1-CNO2B2)*VLNOB(NU(NY,NX),NY,NX)
    DFVP1B=DIFPO*(CP14S1-CP14B2)*VLPOB(NU(NY,NX),NY,NX)
    DFVPOB=DIFPO*(CPO4S1-CPO4B2)*VLPOB(NU(NY,NX),NY,NX)
  ELSE
    DO  K=0,jcplx1
      DFVOC(K)=0.0
      DFVON(K)=0.0
      DFVOP(K)=0.0
      DFVOA(K)=0.0
    ENDDO
    DFVCOS=0.0
    DFVCHS=0.0
    DFVOXS=0.0
    DFVNGS=0.0
    DFVN2S=0.0
    DFVHGS=0.0
    DFVNH4=0.0
    DFVNH3=0.0
    DFVNO3=0.0
    DFVNO2=0.0
    DFVP14=0.0
    DFVPO4=0.0
    DFVN4B=0.0
    DFVN3B=0.0
    DFVNOB=0.0
    DFVN2B=0.0
    DFVP1B=0.0
    DFVPOB=0.0
  ENDIF
  end subroutine DiffusiveFluxAtSoilSurface
!------------------------------------------------------------------------------------------

  subroutine TotalPoreFluxAdjacentCell(NY,NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K
!     R*FLS=convective + diffusive solute flux between litter, soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into soil in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     R*FL0,R*FL1=convective flux into surface litter, soil surface
!     RFL*=convective flux between surface litter and soil surface
!     DFV*=diffusive solute flux between litter and soil surface
!
  DO K=0,2
    ROCFLS(K,3,0,NY,NX)=ROCFL0(K,NY,NX)-RFLOC(K)-DFVOC(K)
    RONFLS(K,3,0,NY,NX)=RONFL0(K,NY,NX)-RFLON(K)-DFVON(K)
    ROPFLS(K,3,0,NY,NX)=ROPFL0(K,NY,NX)-RFLOP(K)-DFVOP(K)
    ROAFLS(K,3,0,NY,NX)=ROAFL0(K,NY,NX)-RFLOA(K)-DFVOA(K)
    ROCFLS(K,3,NU(NY,NX),NY,NX)=ROCFL1(K,NY,NX)+RFLOC(K)+DFVOC(K)
    RONFLS(K,3,NU(NY,NX),NY,NX)=RONFL1(K,NY,NX)+RFLON(K)+DFVON(K)
    ROPFLS(K,3,NU(NY,NX),NY,NX)=ROPFL1(K,NY,NX)+RFLOP(K)+DFVOP(K)
    ROAFLS(K,3,NU(NY,NX),NY,NX)=ROAFL1(K,NY,NX)+RFLOA(K)+DFVOA(K)
  ENDDO
  RCOFLS(3,0,NY,NX)=RCOFL0(NY,NX)+RCOFLS0-RFLCOS-DFVCOS
!     IF(I.GE.360)THEN
!     WRITE(*,8788)'RCOFLS0',I,J,M,MM,NX,NY,NU(NY,NX)
!    2,RCOFLS(3,0,NY,NX)
!    2,RCOFL0(NY,NX),RCOFLS0,RFLCOS,DFVCOS,VFLW,CO2S2(0,NY,NX)
!    3,DIFCS,CCO2S1,CCO2S2,DIFCS0,DIFCS1,CO2S2(NU(NY,NX),NY,NX)
!    4,VOLWM(M,NU(NY,NX),NY,NX)
!8788  FORMAT(A8,7I4,20E12.4)
!     ENDIF
  RCHFLS(3,0,NY,NX)=RCHFL0(NY,NX)+RCHFLS0-RFLCHS-DFVCHS
  ROXFLS(3,0,NY,NX)=ROXFL0(NY,NX)+ROXFLS0-RFLOXS-DFVOXS
  RNGFLS(3,0,NY,NX)=RNGFL0(NY,NX)+RNGFLS0-RFLNGS-DFVNGS
  RN2FLS(3,0,NY,NX)=RN2FL0(NY,NX)+RN2FLS0-RFLN2S-DFVN2S
  RHGFLS(3,0,NY,NX)=RHGFL0(NY,NX)-RFLHGS-DFVHGS
  RN4FLW(3,0,NY,NX)=RN4FL0(NY,NX)+RN4FLW0-RFLNH4-DFVNH4-RFLN4B-DFVN4B
  RN3FLW(3,0,NY,NX)=RN3FL0(NY,NX)+RN3FLW0-RFLNH3-DFVNH3-RFLN3B-DFVN3B
  RNOFLW(3,0,NY,NX)=RNOFL0(NY,NX)+RNOFLW0-RFLNO3-DFVNO3-RFLNOB-DFVNOB
  RNXFLS(3,0,NY,NX)=RNXFL0(NY,NX)-RFLNO2-DFVNO2-RFLN2B-DFVN2B
  RH1PFS(3,0,NY,NX)=RH1PF0(NY,NX)+RH1PFS0-RFLP14-DFVP14-RFLP1B-DFVP1B
  RH2PFS(3,0,NY,NX)=RH2PF0(NY,NX)+RH2PFS0-RFLPO4-DFVPO4-RFLPOB-DFVPOB
  RCOFLS(3,NU(NY,NX),NY,NX)=RCOFL1(NY,NX)+RCOFLS1+RFLCOS+DFVCOS
  RCHFLS(3,NU(NY,NX),NY,NX)=RCHFL1(NY,NX)+RCHFLS1+RFLCHS+DFVCHS
  ROXFLS(3,NU(NY,NX),NY,NX)=ROXFL1(NY,NX)+ROXFLS1+RFLOXS+DFVOXS
  RNGFLS(3,NU(NY,NX),NY,NX)=RNGFL1(NY,NX)+RNGFLS1+RFLNGS+DFVNGS
  RN2FLS(3,NU(NY,NX),NY,NX)=RN2FL1(NY,NX)+RN2FLS1+RFLN2S+DFVN2S
  RHGFLS(3,NU(NY,NX),NY,NX)=RHGFL1(NY,NX)+RFLHGS+DFVHGS
  RN4FLW(3,NU(NY,NX),NY,NX)=RN4FL1(NY,NX)+RN4FLW1+RFLNH4+DFVNH4
  RN3FLW(3,NU(NY,NX),NY,NX)=RN3FL1(NY,NX)+RN3FLW1+RFLNH3+DFVNH3
  RNOFLW(3,NU(NY,NX),NY,NX)=RNOFL1(NY,NX)+RNOFLW1+RFLNO3+DFVNO3
  RNXFLS(3,NU(NY,NX),NY,NX)=RNXFL1(NY,NX)+RFLNO2+DFVNO2
  RH1PFS(3,NU(NY,NX),NY,NX)=RH1PF1(NY,NX)+RH1PFS1+RFLP14+DFVP14
  RH2PFS(3,NU(NY,NX),NY,NX)=RH2PF1(NY,NX)+RH2PFS1+RFLPO4+DFVPO4
  RN4FLB(3,NU(NY,NX),NY,NX)=RN4FL2(NY,NX)+RN4FLB1+RFLN4B+DFVN4B
  RN3FLB(3,NU(NY,NX),NY,NX)=RN3FL2(NY,NX)+RN3FLB1+RFLN3B+DFVN3B
  RNOFLB(3,NU(NY,NX),NY,NX)=RNOFL2(NY,NX)+RNOFLB1+RFLNOB+DFVNOB
  RNXFLB(3,NU(NY,NX),NY,NX)=RNXFL2(NY,NX)+RFLN2B+DFVN2B
  RH1BFB(3,NU(NY,NX),NY,NX)=RH1BF2(NY,NX)+RH1BFB1+RFLP1B+DFVP1B
  RH2BFB(3,NU(NY,NX),NY,NX)=RH2BF2(NY,NX)+RH2BFB1+RFLPOB+DFVPOB
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  DO K=0,2
    XOCFLS(K,3,0,NY,NX)=XOCFLS(K,3,0,NY,NX)-RFLOC(K)-DFVOC(K)
    XONFLS(K,3,0,NY,NX)=XONFLS(K,3,0,NY,NX)-RFLON(K)-DFVON(K)
    XOPFLS(K,3,0,NY,NX)=XOPFLS(K,3,0,NY,NX)-RFLOP(K)-DFVOP(K)
    XOAFLS(K,3,0,NY,NX)=XOAFLS(K,3,0,NY,NX)-RFLOA(K)-DFVOA(K)
    XOCFLS(K,3,NU(NY,NX),NY,NX)=XOCFLS(K,3,NU(NY,NX),NY,NX)+RFLOC(K)+DFVOC(K)
    XONFLS(K,3,NU(NY,NX),NY,NX)=XONFLS(K,3,NU(NY,NX),NY,NX)+RFLON(K)+DFVON(K)
    XOPFLS(K,3,NU(NY,NX),NY,NX)=XOPFLS(K,3,NU(NY,NX),NY,NX)+RFLOP(K)+DFVOP(K)
    XOAFLS(K,3,NU(NY,NX),NY,NX)=XOAFLS(K,3,NU(NY,NX),NY,NX)+RFLOA(K)+DFVOA(K)
  ENDDO
  XCOFLS(3,0,NY,NX)=XCOFLS(3,0,NY,NX)+RCOFLS0-RFLCOS-DFVCOS
  XCHFLS(3,0,NY,NX)=XCHFLS(3,0,NY,NX)+RCHFLS0-RFLCHS-DFVCHS
  XOXFLS(3,0,NY,NX)=XOXFLS(3,0,NY,NX)+ROXFLS0-RFLOXS-DFVOXS
  XNGFLS(3,0,NY,NX)=XNGFLS(3,0,NY,NX)+RNGFLS0-RFLNGS-DFVNGS
  XN2FLS(3,0,NY,NX)=XN2FLS(3,0,NY,NX)+RN2FLS0-RFLN2S-DFVN2S
  XHGFLS(3,0,NY,NX)=XHGFLS(3,0,NY,NX)-RFLHGS-DFVHGS
  XN4FLW(3,0,NY,NX)=XN4FLW(3,0,NY,NX)+RN4FLW0-RFLNH4-DFVNH4-RFLN4B-DFVN4B
  XN3FLW(3,0,NY,NX)=XN3FLW(3,0,NY,NX)+RN3FLW0-RFLNH3-DFVNH3-RFLN3B-DFVN3B
  XNOFLW(3,0,NY,NX)=XNOFLW(3,0,NY,NX)+RNOFLW0-RFLNO3-DFVNO3-RFLNOB-DFVNOB
  XNXFLS(3,0,NY,NX)=XNXFLS(3,0,NY,NX)-RFLNO2-DFVNO2-RFLN2B-DFVN2B
  XH1PFS(3,0,NY,NX)=XH1PFS(3,0,NY,NX)+RH1PFS0-RFLP14-DFVP14-RFLP1B-DFVP1B
  XH2PFS(3,0,NY,NX)=XH2PFS(3,0,NY,NX)+RH2PFS0-RFLPO4-DFVPO4-RFLPOB-DFVPOB
  XCOFLS(3,NU(NY,NX),NY,NX)=XCOFLS(3,NU(NY,NX),NY,NX)+RCOFLS1+RFLCOS+DFVCOS
  XCHFLS(3,NU(NY,NX),NY,NX)=XCHFLS(3,NU(NY,NX),NY,NX)+RCHFLS1+RFLCHS+DFVCHS
  XOXFLS(3,NU(NY,NX),NY,NX)=XOXFLS(3,NU(NY,NX),NY,NX)+ROXFLS1+RFLOXS+DFVOXS
  XNGFLS(3,NU(NY,NX),NY,NX)=XNGFLS(3,NU(NY,NX),NY,NX)+RNGFLS1+RFLNGS+DFVNGS
  XN2FLS(3,NU(NY,NX),NY,NX)=XN2FLS(3,NU(NY,NX),NY,NX)+RN2FLS1+RFLN2S+DFVN2S
  XHGFLS(3,NU(NY,NX),NY,NX)=XHGFLS(3,NU(NY,NX),NY,NX)+RFLHGS+DFVHGS
  XN4FLW(3,NU(NY,NX),NY,NX)=XN4FLW(3,NU(NY,NX),NY,NX)+RN4FLW1+RFLNH4+DFVNH4
  XN3FLW(3,NU(NY,NX),NY,NX)=XN3FLW(3,NU(NY,NX),NY,NX)+RN3FLW1+RFLNH3+DFVNH3
  XNOFLW(3,NU(NY,NX),NY,NX)=XNOFLW(3,NU(NY,NX),NY,NX)+RNOFLW1+RFLNO3+DFVNO3
  XNXFLS(3,NU(NY,NX),NY,NX)=XNXFLS(3,NU(NY,NX),NY,NX)+RFLNO2+DFVNO2
  XH1PFS(3,NU(NY,NX),NY,NX)=XH1PFS(3,NU(NY,NX),NY,NX)+RH1PFS1+RFLP14+DFVP14
  XH2PFS(3,NU(NY,NX),NY,NX)=XH2PFS(3,NU(NY,NX),NY,NX)+RH2PFS1+RFLPO4+DFVPO4
  XN4FLB(3,NU(NY,NX),NY,NX)=XN4FLB(3,NU(NY,NX),NY,NX)+RN4FLB1+RFLN4B+DFVN4B
  XN3FLB(3,NU(NY,NX),NY,NX)=XN3FLB(3,NU(NY,NX),NY,NX)+RN3FLB1+RFLN3B+DFVN3B
  XNOFLB(3,NU(NY,NX),NY,NX)=XNOFLB(3,NU(NY,NX),NY,NX)+RNOFLB1+RFLNOB+DFVNOB
  XNXFLB(3,NU(NY,NX),NY,NX)=XNXFLB(3,NU(NY,NX),NY,NX)+RFLN2B+DFVN2B
  XH1BFB(3,NU(NY,NX),NY,NX)=XH1BFB(3,NU(NY,NX),NY,NX)+RH1BFB1+RFLP1B+DFVP1B
  XH2BFB(3,NU(NY,NX),NY,NX)=XH2BFB(3,NU(NY,NX),NY,NX)+RH2BFB1+RFLPOB+DFVPOB
  end subroutine TotalPoreFluxAdjacentCell
!------------------------------------------------------------------------------------------

  subroutine MacroMicroTransfer(M,NY,NX)
  implicit none

  integer, intent(in) :: M, NY, NX
  integer :: K
  real(r8) :: VFLW
  real(r8) :: VOLWHS,VOLWT

!
!     FINHM=macro-micropore water transfer from watsub.f
!     VOLWM,VOLWHM=micropore,macropore water volume
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
  IF(FINHM(M,NU(NY,NX),NY,NX).GT.0.0)THEN
    IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AMAX1(0.0,AMIN1(VFLWX,FINHM(M,NU(NY,NX),NY,NX) &
      /VOLWHM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQCH2(K,NU(NY,NX),NY,NX))
      RFLON(K)=VFLW*AMAX1(0.0,OQNH2(K,NU(NY,NX),NY,NX))
      RFLOP(K)=VFLW*AMAX1(0.0,OQPH2(K,NU(NY,NX),NY,NX))
      RFLOA(K)=VFLW*AMAX1(0.0,OQAH2(K,NU(NY,NX),NY,NX))
    ENDDO
    RFLCOS=VFLW*AMAX1(0.0,CO2SH2(NU(NY,NX),NY,NX))
    RFLCHS=VFLW*AMAX1(0.0,CH4SH2(NU(NY,NX),NY,NX))
    RFLOXS=VFLW*AMAX1(0.0,OXYSH2(NU(NY,NX),NY,NX))
    RFLNGS=VFLW*AMAX1(0.0,Z2GSH2(NU(NY,NX),NY,NX))
    RFLN2S=VFLW*AMAX1(0.0,Z2OSH2(NU(NY,NX),NY,NX))
    RFLHGS=VFLW*AMAX1(0.0,H2GSH2(NU(NY,NX),NY,NX))
    RFLNH4=VFLW*AMAX1(0.0,ZNH4H2(NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNH3=VFLW*AMAX1(0.0,ZNH3H2(NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNO3=VFLW*AMAX1(0.0,ZNO3H2(NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLNO2=VFLW*AMAX1(0.0,ZNO2H2(NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLP14=VFLW*AMAX1(0.0,H1P4H2(NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLPO4=VFLW*AMAX1(0.0,H2P4H2(NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLN4B=VFLW*AMAX1(0.0,ZN4BH2(NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLN3B=VFLW*AMAX1(0.0,ZN3BH2(NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLNOB=VFLW*AMAX1(0.0,ZNOBH2(NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLN2B=VFLW*AMAX1(0.0,ZN2BH2(NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLP1B=VFLW*AMAX1(0.0,H1PBH2(NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
    RFLPOB=VFLW*AMAX1(0.0,H2PBH2(NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FINHM(M,NU(NY,NX),NY,NX).LT.0.0)THEN
    IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AMIN1(0.0,AMAX1(-VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO  K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQC2(K,NU(NY,NX),NY,NX))
      RFLON(K)=VFLW*AMAX1(0.0,OQN2(K,NU(NY,NX),NY,NX))
      RFLOP(K)=VFLW*AMAX1(0.0,OQP2(K,NU(NY,NX),NY,NX))
      RFLOA(K)=VFLW*AMAX1(0.0,OQA2(K,NU(NY,NX),NY,NX))
    ENDDO
    RFLCOS=VFLW*AMAX1(0.0,CO2S2(NU(NY,NX),NY,NX))
    RFLCHS=VFLW*AMAX1(0.0,CH4S2(NU(NY,NX),NY,NX))
    RFLOXS=VFLW*AMAX1(0.0,OXYS2(NU(NY,NX),NY,NX))
    RFLNGS=VFLW*AMAX1(0.0,Z2GS2(NU(NY,NX),NY,NX))
    RFLN2S=VFLW*AMAX1(0.0,Z2OS2(NU(NY,NX),NY,NX))
    RFLHGS=VFLW*AMAX1(0.0,H2GS2(NU(NY,NX),NY,NX))
    RFLNH4=VFLW*AMAX1(0.0,ZNH4S2(NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNH3=VFLW*AMAX1(0.0,ZNH3S2(NU(NY,NX),NY,NX))*VLNH4(NU(NY,NX),NY,NX)
    RFLNO3=VFLW*AMAX1(0.0,ZNO3S2(NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLNO2=VFLW*AMAX1(0.0,ZNO2S2(NU(NY,NX),NY,NX))*VLNO3(NU(NY,NX),NY,NX)
    RFLP14=VFLW*AMAX1(0.0,H1PO42(NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLPO4=VFLW*AMAX1(0.0,H2PO42(NU(NY,NX),NY,NX))*VLPO4(NU(NY,NX),NY,NX)
    RFLN4B=VFLW*AMAX1(0.0,ZNH4B2(NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLN3B=VFLW*AMAX1(0.0,ZNH3B2(NU(NY,NX),NY,NX))*VLNHB(NU(NY,NX),NY,NX)
    RFLNOB=VFLW*AMAX1(0.0,ZNO3B2(NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLN2B=VFLW*AMAX1(0.0,ZNO2B2(NU(NY,NX),NY,NX))*VLNOB(NU(NY,NX),NY,NX)
    RFLP1B=VFLW*AMAX1(0.0,H1POB2(NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
    RFLPOB=VFLW*AMAX1(0.0,H2POB2(NU(NY,NX),NY,NX))*VLPOB(NU(NY,NX),NY,NX)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    DO  K=0,jcplx1
      RFLOC(K)=0.0
      RFLON(K)=0.0
      RFLOP(K)=0.0
      RFLOA(K)=0.0
    ENDDO
    RFLCOS=0.0
    RFLCHS=0.0
    RFLOXS=0.0
    RFLNGS=0.0
    RFLN2S=0.0
    RFLHGS=0.0
    RFLNH4=0.0
    RFLNH3=0.0
    RFLNO3=0.0
    RFLNO2=0.0
    RFLP14=0.0
    RFLPO4=0.0
    RFLN4B=0.0
    RFLN3B=0.0
    RFLNOB=0.0
    RFLN2B=0.0
    RFLP1B=0.0
    RFLPOB=0.0
  ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
!     VOLWM,VOLWHM=micropore,macropore water volume
!     XFRS*VOLT=maximum macropore volume for solute transfer
!     DFV*=diffusive macropore-micropore solute transfer
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *H2,*2=macropore,micropore solute content
!
  IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VOLWHS=AMIN1(XFRS*VOLT(NU(NY,NX),NY,NX),VOLWHM(M,NU(NY,NX),NY,NX))
    VOLWT=VOLWM(M,NU(NY,NX),NY,NX)+VOLWHS
    DO  K=0,jcplx1
      DFVOC(K)=XNPH*(AMAX1(0.0,OQCH2(K,NU(NY,NX),NY,NX))*VOLWM(M,NU(NY,NX),NY,NX) &
        -AMAX1(0.0,OQC2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
      DFVON(K)=XNPH*(AMAX1(0.0,OQNH2(K,NU(NY,NX),NY,NX)) &
        *VOLWM(M,NU(NY,NX),NY,NX) &
        -AMAX1(0.0,OQN2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
      DFVOP(K)=XNPH*(AMAX1(0.0,OQPH2(K,NU(NY,NX),NY,NX)) &
        *VOLWM(M,NU(NY,NX),NY,NX) &
        -AMAX1(0.0,OQP2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
      DFVOA(K)=XNPH*(AMAX1(0.0,OQAH2(K,NU(NY,NX),NY,NX)) &
        *VOLWM(M,NU(NY,NX),NY,NX) &
        -AMAX1(0.0,OQA2(K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    ENDDO
    DFVCOS=XNPH*(AMAX1(0.0,CO2SH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,CO2S2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    DFVCHS=XNPH*(AMAX1(0.0,CH4SH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,CH4S2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    DFVOXS=XNPH*(AMAX1(0.0,OXYSH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,OXYS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    DFVNGS=XNPH*(AMAX1(0.0,Z2GSH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,Z2GS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    DFVN2S=XNPH*(AMAX1(0.0,Z2OSH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,Z2OS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    DFVHGS=XNPH*(AMAX1(0.0,H2GSH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,H2GS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
    DFVNH4=XNPH*(AMAX1(0.0,ZNH4H2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNH4S2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNH4(NU(NY,NX),NY,NX)
    DFVNH3=XNPH*(AMAX1(0.0,ZNH3H2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNH3S2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNH4(NU(NY,NX),NY,NX)
    DFVNO3=XNPH*(AMAX1(0.0,ZNO3H2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNO3S2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNO3(NU(NY,NX),NY,NX)
    DFVNO2=XNPH*(AMAX1(0.0,ZNO2H2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNO2S2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNO3(NU(NY,NX),NY,NX)
    DFVP14=XNPH*(AMAX1(0.0,H1P4H2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,H1PO42(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPO4(NU(NY,NX),NY,NX)
    DFVPO4=XNPH*(AMAX1(0.0,H2P4H2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,H2PO42(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPO4(NU(NY,NX),NY,NX)
    DFVN4B=XNPH*(AMAX1(0.0,ZN4BH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNH4B2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNHB(NU(NY,NX),NY,NX)
    DFVN3B=XNPH*(AMAX1(0.0,ZN3BH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNH3B2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNHB(NU(NY,NX),NY,NX)
    DFVNOB=XNPH*(AMAX1(0.0,ZNOBH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNO3B2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNOB(NU(NY,NX),NY,NX)
    DFVN2B=XNPH*(AMAX1(0.0,ZN2BH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,ZNO2B2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLNOB(NU(NY,NX),NY,NX)
    DFVP1B=XNPH*(AMAX1(0.0,H1PBH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,H1POB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPOB(NU(NY,NX),NY,NX)
    DFVPOB=XNPH*(AMAX1(0.0,H2PBH2(NU(NY,NX),NY,NX)) &
      *VOLWM(M,NU(NY,NX),NY,NX) &
      -AMAX1(0.0,H2POB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT &
      *VLPOB(NU(NY,NX),NY,NX)
  ELSE
    DO  K=0,jcplx1
      DFVOC(K)=0.0
      DFVON(K)=0.0
      DFVOP(K)=0.0
      DFVOA(K)=0.0
    ENDDO
    DFVCOS=0.0
    DFVCHS=0.0
    DFVOXS=0.0
    DFVNGS=0.0
    DFVN2S=0.0
    DFVHGS=0.0
    DFVNH4=0.0
    DFVNH3=0.0
    DFVNO3=0.0
    DFVNO2=0.0
    DFVP14=0.0
    DFVPO4=0.0
    DFVN4B=0.0
    DFVN3B=0.0
    DFVNOB=0.0
    DFVN2B=0.0
    DFVP1B=0.0
    DFVPOB=0.0
  ENDIF
!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  DO  K=0,jcplx1
    ROCFXS(K,NU(NY,NX),NY,NX)=RFLOC(K)+DFVOC(K)
    RONFXS(K,NU(NY,NX),NY,NX)=RFLON(K)+DFVON(K)
    ROPFXS(K,NU(NY,NX),NY,NX)=RFLOP(K)+DFVOP(K)
    ROAFXS(K,NU(NY,NX),NY,NX)=RFLOA(K)+DFVOA(K)
  ENDDO
  RCOFXS(NU(NY,NX),NY,NX)=RFLCOS+DFVCOS
  RCHFXS(NU(NY,NX),NY,NX)=RFLCHS+DFVCHS
  ROXFXS(NU(NY,NX),NY,NX)=RFLOXS+DFVOXS
  RNGFXS(NU(NY,NX),NY,NX)=RFLNGS+DFVNGS
  RN2FXS(NU(NY,NX),NY,NX)=RFLN2S+DFVN2S
  RHGFXS(NU(NY,NX),NY,NX)=RFLHGS+DFVHGS
  RN4FXW(NU(NY,NX),NY,NX)=RFLNH4+DFVNH4
  RN3FXW(NU(NY,NX),NY,NX)=RFLNH3+DFVNH3
  RNOFXW(NU(NY,NX),NY,NX)=RFLNO3+DFVNO3
  RNXFXS(NU(NY,NX),NY,NX)=RFLNO2+DFVNO2
  RH1PXS(NU(NY,NX),NY,NX)=RFLP14+DFVP14
  RH2PXS(NU(NY,NX),NY,NX)=RFLPO4+DFVPO4
  RN4FXB(NU(NY,NX),NY,NX)=RFLN4B+DFVN4B
  RN3FXB(NU(NY,NX),NY,NX)=RFLN3B+DFVN3B
  RNOFXB(NU(NY,NX),NY,NX)=RFLNOB+DFVNOB
  RNXFXB(NU(NY,NX),NY,NX)=RFLN2B+DFVN2B
  RH1BXB(NU(NY,NX),NY,NX)=RFLP1B+DFVP1B
  RH2BXB(NU(NY,NX),NY,NX)=RFLPOB+DFVPOB
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FXS=hourly convective + diffusive solute flux between macropores and micropores
!     R*FXS=total convective + diffusive solute flux between macropores and micropores
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  DO  K=0,jcplx1
    XOCFXS(K,NU(NY,NX),NY,NX)=XOCFXS(K,NU(NY,NX),NY,NX) &
      +ROCFXS(K,NU(NY,NX),NY,NX)
    XONFXS(K,NU(NY,NX),NY,NX)=XONFXS(K,NU(NY,NX),NY,NX) &
      +RONFXS(K,NU(NY,NX),NY,NX)
    XOPFXS(K,NU(NY,NX),NY,NX)=XOPFXS(K,NU(NY,NX),NY,NX) &
      +ROPFXS(K,NU(NY,NX),NY,NX)
    XOAFXS(K,NU(NY,NX),NY,NX)=XOAFXS(K,NU(NY,NX),NY,NX) &
      +ROAFXS(K,NU(NY,NX),NY,NX)
  ENDDO
  XCOFXS(NU(NY,NX),NY,NX)=XCOFXS(NU(NY,NX),NY,NX) &
    +RCOFXS(NU(NY,NX),NY,NX)
  XCHFXS(NU(NY,NX),NY,NX)=XCHFXS(NU(NY,NX),NY,NX) &
    +RCHFXS(NU(NY,NX),NY,NX)
  XOXFXS(NU(NY,NX),NY,NX)=XOXFXS(NU(NY,NX),NY,NX) &
    +ROXFXS(NU(NY,NX),NY,NX)
  XNGFXS(NU(NY,NX),NY,NX)=XNGFXS(NU(NY,NX),NY,NX) &
    +RNGFXS(NU(NY,NX),NY,NX)
  XN2FXS(NU(NY,NX),NY,NX)=XN2FXS(NU(NY,NX),NY,NX) &
    +RN2FXS(NU(NY,NX),NY,NX)
  XHGFXS(NU(NY,NX),NY,NX)=XHGFXS(NU(NY,NX),NY,NX) &
    +RHGFXS(NU(NY,NX),NY,NX)
  XN4FXW(NU(NY,NX),NY,NX)=XN4FXW(NU(NY,NX),NY,NX) &
    +RN4FXW(NU(NY,NX),NY,NX)
  XN3FXW(NU(NY,NX),NY,NX)=XN3FXW(NU(NY,NX),NY,NX) &
    +RN3FXW(NU(NY,NX),NY,NX)
  XNOFXW(NU(NY,NX),NY,NX)=XNOFXW(NU(NY,NX),NY,NX) &
    +RNOFXW(NU(NY,NX),NY,NX)
  XNXFXS(NU(NY,NX),NY,NX)=XNXFXS(NU(NY,NX),NY,NX) &
    +RNXFXS(NU(NY,NX),NY,NX)
  XH1PXS(NU(NY,NX),NY,NX)=XH1PXS(NU(NY,NX),NY,NX) &
    +RH1PXS(NU(NY,NX),NY,NX)
  XH2PXS(NU(NY,NX),NY,NX)=XH2PXS(NU(NY,NX),NY,NX) &
    +RH2PXS(NU(NY,NX),NY,NX)
  XN4FXB(NU(NY,NX),NY,NX)=XN4FXB(NU(NY,NX),NY,NX) &
    +RN4FXB(NU(NY,NX),NY,NX)
  XN3FXB(NU(NY,NX),NY,NX)=XN3FXB(NU(NY,NX),NY,NX) &
    +RN3FXB(NU(NY,NX),NY,NX)
  XNOFXB(NU(NY,NX),NY,NX)=XNOFXB(NU(NY,NX),NY,NX) &
    +RNOFXB(NU(NY,NX),NY,NX)
  XNXFXB(NU(NY,NX),NY,NX)=XNXFXB(NU(NY,NX),NY,NX) &
    +RNXFXB(NU(NY,NX),NY,NX)
  XH1BXB(NU(NY,NX),NY,NX)=XH1BXB(NU(NY,NX),NY,NX) &
    +RH1BXB(NU(NY,NX),NY,NX)
  XH2BXB(NU(NY,NX),NY,NX)=XH2BXB(NU(NY,NX),NY,NX) &
    +RH2BXB(NU(NY,NX),NY,NX)
!     IF(I.EQ.133)THEN
!     WRITE(*,441)'RH1BXB',I,J,M,NX,NY,RH1BXB(NU(NY,NX),NY,NX)
!    2,RFLP1B,DFVP1B,XH1BXB(NU(NY,NX),NY,NX),H1PBH2(NU(NY,NX),NY,NX)
!    2,H1POB2(NU(NY,NX),NY,NX)
!441   FORMAT(A8,5I4,20E17.8)
!     ENDIF
!
  end subroutine MacroMicroTransfer
!------------------------------------------------------------------------------------------

  subroutine OverlandFlowSnowdriftTransport(M,NY,NX,NHE,NHW,NVS,NVN)
  implicit none

  integer, intent(in) :: M, NY, NX, NHE, NHW, NVS, NVN
  real(r8) :: FQRM,VFLW
  integer :: K,N,NN,N1,N2,N4,N5,N4B,N5B
!     QRM=runoff from watsub.f
!     RQR*=solute in runoff
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VOLWM=litter water volume from watsub.f
!     *S2=litter solute content
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     X*QRS=accumulated hourly solute in runoff
!
  N1=NX
  N2=NY
  IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
    IF(VOLWM(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,QRM(M,N2,N1)/VOLWM(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=0,jcplx1
      RQROC0(K,N2,N1)=VFLW*AMAX1(0.0,OQC2(K,0,N2,N1))
      RQRON0(K,N2,N1)=VFLW*AMAX1(0.0,OQN2(K,0,N2,N1))
      RQROP0(K,N2,N1)=VFLW*AMAX1(0.0,OQP2(K,0,N2,N1))
      RQROA0(K,N2,N1)=VFLW*AMAX1(0.0,OQA2(K,0,N2,N1))
    ENDDO
    RQRCOS0(N2,N1)=VFLW*AMAX1(0.0,CO2S2(0,N2,N1))
    RQRCHS0(N2,N1)=VFLW*AMAX1(0.0,CH4S2(0,N2,N1))
    RQROXS0(N2,N1)=VFLW*AMAX1(0.0,OXYS2(0,N2,N1))
    RQRNGS0(N2,N1)=VFLW*AMAX1(0.0,Z2GS2(0,N2,N1))
    RQRN2S0(N2,N1)=VFLW*AMAX1(0.0,Z2OS2(0,N2,N1))
    RQRHGS0(N2,N1)=VFLW*AMAX1(0.0,H2GS2(0,N2,N1))
    RQRNH40(N2,N1)=VFLW*AMAX1(0.0,ZNH4S2(0,N2,N1))
    RQRNH30(N2,N1)=VFLW*AMAX1(0.0,ZNH3S2(0,N2,N1))
    RQRNO30(N2,N1)=VFLW*AMAX1(0.0,ZNO3S2(0,N2,N1))
    RQRNO20(N2,N1)=VFLW*AMAX1(0.0,ZNO2S2(0,N2,N1))
    RQRH1P0(N2,N1)=VFLW*AMAX1(0.0,H1PO42(0,N2,N1))
    RQRH2P0(N2,N1)=VFLW*AMAX1(0.0,H2PO42(0,N2,N1))
  ELSE
    DO K=0,jcplx1
      RQROC0(K,N2,N1)=0.0
      RQRON0(K,N2,N1)=0.0
      RQROP0(K,N2,N1)=0.0
      RQROA0(K,N2,N1)=0.0
    ENDDO
    RQRCOS0(N2,N1)=0.0
    RQRCHS0(N2,N1)=0.0
    RQROXS0(N2,N1)=0.0
    RQRNGS0(N2,N1)=0.0
    RQRN2S0(N2,N1)=0.0
    RQRHGS0(N2,N1)=0.0
    RQRNH40(N2,N1)=0.0
    RQRNH30(N2,N1)=0.0
    RQRNO30(N2,N1)=0.0
    RQRNO20(N2,N1)=0.0
    RQRH1P0(N2,N1)=0.0
    RQRH2P0(N2,N1)=0.0
  ENDIF
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO 4310 N=1,2
    DO 4305 NN=1,2
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE.AND.NN.EQ.1.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N4B=NX-1
          N5B=NY
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS.AND.NN.EQ.1.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N4B=NX
          N5B=NY-1
        ENDIF
      ENDIF
!
!     QRM=runoff from watsub.f
!     RQR*=solute in runoff
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VOLWM=litter water volume from watsub.f
!     *S2=litter solute content
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        IF(NN.EQ.1)THEN
          FQRM=QRMN(M,N,2,N5,N4)/QRM(M,N2,N1)
          DO  K=0,jcplx1
            RQROC(K,N,2,N5,N4)=RQROC0(K,N2,N1)*FQRM
            RQRON(K,N,2,N5,N4)=RQRON0(K,N2,N1)*FQRM
            RQROP(K,N,2,N5,N4)=RQROP0(K,N2,N1)*FQRM
            RQROA(K,N,2,N5,N4)=RQROA0(K,N2,N1)*FQRM
          ENDDO
          RQRCOS(N,2,N5,N4)=RQRCOS0(N2,N1)*FQRM
          RQRCHS(N,2,N5,N4)=RQRCHS0(N2,N1)*FQRM
          RQROXS(N,2,N5,N4)=RQROXS0(N2,N1)*FQRM
          RQRNGS(N,2,N5,N4)=RQRNGS0(N2,N1)*FQRM
          RQRN2S(N,2,N5,N4)=RQRN2S0(N2,N1)*FQRM
          RQRHGS(N,2,N5,N4)=RQRHGS0(N2,N1)*FQRM
          RQRNH4(N,2,N5,N4)=RQRNH40(N2,N1)*FQRM
          RQRNH3(N,2,N5,N4)=RQRNH30(N2,N1)*FQRM
          RQRNO3(N,2,N5,N4)=RQRNO30(N2,N1)*FQRM
          RQRNO2(N,2,N5,N4)=RQRNO20(N2,N1)*FQRM
          RQRH1P(N,2,N5,N4)=RQRH1P0(N2,N1)*FQRM
          RQRH2P(N,2,N5,N4)=RQRH2P0(N2,N1)*FQRM
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
          DO  K=0,jcplx1
            XOCQRS(K,N,2,N5,N4)=XOCQRS(K,N,2,N5,N4)+RQROC(K,N,2,N5,N4)
            XONQRS(K,N,2,N5,N4)=XONQRS(K,N,2,N5,N4)+RQRON(K,N,2,N5,N4)
            XOPQRS(K,N,2,N5,N4)=XOPQRS(K,N,2,N5,N4)+RQROP(K,N,2,N5,N4)
            XOAQRS(K,N,2,N5,N4)=XOAQRS(K,N,2,N5,N4)+RQROA(K,N,2,N5,N4)
          ENDDO
          XCOQRS(N,2,N5,N4)=XCOQRS(N,2,N5,N4)+RQRCOS(N,2,N5,N4)
          XCHQRS(N,2,N5,N4)=XCHQRS(N,2,N5,N4)+RQRCHS(N,2,N5,N4)
          XOXQRS(N,2,N5,N4)=XOXQRS(N,2,N5,N4)+RQROXS(N,2,N5,N4)
          XNGQRS(N,2,N5,N4)=XNGQRS(N,2,N5,N4)+RQRNGS(N,2,N5,N4)
          XN2QRS(N,2,N5,N4)=XN2QRS(N,2,N5,N4)+RQRN2S(N,2,N5,N4)
          XHGQRS(N,2,N5,N4)=XHGQRS(N,2,N5,N4)+RQRHGS(N,2,N5,N4)
          XN4QRW(N,2,N5,N4)=XN4QRW(N,2,N5,N4)+RQRNH4(N,2,N5,N4)
          XN3QRW(N,2,N5,N4)=XN3QRW(N,2,N5,N4)+RQRNH3(N,2,N5,N4)
          XNOQRW(N,2,N5,N4)=XNOQRW(N,2,N5,N4)+RQRNO3(N,2,N5,N4)
          XNXQRS(N,2,N5,N4)=XNXQRS(N,2,N5,N4)+RQRNO2(N,2,N5,N4)
          XP1QRW(N,2,N5,N4)=XP1QRW(N,2,N5,N4)+RQRH1P(N,2,N5,N4)
          XP4QRW(N,2,N5,N4)=XP4QRW(N,2,N5,N4)+RQRH2P(N,2,N5,N4)
        ELSE
          DO K=0,jcplx1
            RQROC(K,N,2,N5,N4)=0.0
            RQRON(K,N,2,N5,N4)=0.0
            RQROP(K,N,2,N5,N4)=0.0
            RQROA(K,N,2,N5,N4)=0.0
          ENDDO
          RQRCOS(N,2,N5,N4)=0.0
          RQRCHS(N,2,N5,N4)=0.0
          RQROXS(N,2,N5,N4)=0.0
          RQRNGS(N,2,N5,N4)=0.0
          RQRN2S(N,2,N5,N4)=0.0
          RQRHGS(N,2,N5,N4)=0.0
          RQRNH4(N,2,N5,N4)=0.0
          RQRNH3(N,2,N5,N4)=0.0
          RQRNO3(N,2,N5,N4)=0.0
          RQRNO2(N,2,N5,N4)=0.0
          RQRH1P(N,2,N5,N4)=0.0
          RQRH2P(N,2,N5,N4)=0.0
        ENDIF
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            FQRM=QRMN(M,N,1,N5B,N4B)/QRM(M,N2,N1)
            DO  K=0,jcplx1
              RQROC(K,N,1,N5B,N4B)=RQROC0(K,N2,N1)*FQRM
              RQRON(K,N,1,N5B,N4B)=RQRON0(K,N2,N1)*FQRM
              RQROP(K,N,1,N5B,N4B)=RQROP0(K,N2,N1)*FQRM
              RQROA(K,N,1,N5B,N4B)=RQROA0(K,N2,N1)*FQRM
            ENDDO
            RQRCOS(N,1,N5B,N4B)=RQRCOS0(N2,N1)*FQRM
            RQRCHS(N,1,N5B,N4B)=RQRCHS0(N2,N1)*FQRM
            RQROXS(N,1,N5B,N4B)=RQROXS0(N2,N1)*FQRM
            RQRNGS(N,1,N5B,N4B)=RQRNGS0(N2,N1)*FQRM
            RQRN2S(N,1,N5B,N4B)=RQRN2S0(N2,N1)*FQRM
            RQRHGS(N,1,N5B,N4B)=RQRHGS0(N2,N1)*FQRM
            RQRNH4(N,1,N5B,N4B)=RQRNH40(N2,N1)*FQRM
            RQRNH3(N,1,N5B,N4B)=RQRNH30(N2,N1)*FQRM
            RQRNO3(N,1,N5B,N4B)=RQRNO30(N2,N1)*FQRM
            RQRNO2(N,1,N5B,N4B)=RQRNO20(N2,N1)*FQRM
            RQRH1P(N,1,N5B,N4B)=RQRH1P0(N2,N1)*FQRM
            RQRH2P(N,1,N5B,N4B)=RQRH2P0(N2,N1)*FQRM
            DO K=0,jcplx1
              XOCQRS(K,N,1,N5B,N4B)=XOCQRS(K,N,1,N5B,N4B)+RQROC(K,N,1,N5B,N4B)
              XONQRS(K,N,1,N5B,N4B)=XONQRS(K,N,1,N5B,N4B)+RQRON(K,N,1,N5B,N4B)
              XOPQRS(K,N,1,N5B,N4B)=XOPQRS(K,N,1,N5B,N4B)+RQROP(K,N,1,N5B,N4B)
              XOAQRS(K,N,1,N5B,N4B)=XOAQRS(K,N,1,N5B,N4B)+RQROA(K,N,1,N5B,N4B)
            ENDDO
            XCOQRS(N,1,N5B,N4B)=XCOQRS(N,1,N5B,N4B)+RQRCOS(N,1,N5B,N4B)
            XCHQRS(N,1,N5B,N4B)=XCHQRS(N,1,N5B,N4B)+RQRCHS(N,1,N5B,N4B)
            XOXQRS(N,1,N5B,N4B)=XOXQRS(N,1,N5B,N4B)+RQROXS(N,1,N5B,N4B)
            XNGQRS(N,1,N5B,N4B)=XNGQRS(N,1,N5B,N4B)+RQRNGS(N,1,N5B,N4B)
            XN2QRS(N,1,N5B,N4B)=XN2QRS(N,1,N5B,N4B)+RQRN2S(N,1,N5B,N4B)
            XHGQRS(N,1,N5B,N4B)=XHGQRS(N,1,N5B,N4B)+RQRHGS(N,1,N5B,N4B)
            XN4QRW(N,1,N5B,N4B)=XN4QRW(N,1,N5B,N4B)+RQRNH4(N,1,N5B,N4B)
            XN3QRW(N,1,N5B,N4B)=XN3QRW(N,1,N5B,N4B)+RQRNH3(N,1,N5B,N4B)
            XNOQRW(N,1,N5B,N4B)=XNOQRW(N,1,N5B,N4B)+RQRNO3(N,1,N5B,N4B)
            XNXQRS(N,1,N5B,N4B)=XNXQRS(N,1,N5B,N4B)+RQRNO2(N,1,N5B,N4B)
            XP1QRW(N,1,N5B,N4B)=XP1QRW(N,1,N5B,N4B)+RQRH1P(N,1,N5B,N4B)
            XP4QRW(N,1,N5B,N4B)=XP4QRW(N,1,N5B,N4B)+RQRH2P(N,1,N5B,N4B)
          ELSE
            DO  K=0,jcplx1
              RQROC(K,N,1,N5B,N4B)=0.0
              RQRON(K,N,1,N5B,N4B)=0.0
              RQROP(K,N,1,N5B,N4B)=0.0
              RQROA(K,N,1,N5B,N4B)=0.0
            ENDDO
            RQRCOS(N,1,N5B,N4B)=0.0
            RQRCHS(N,1,N5B,N4B)=0.0
            RQROXS(N,1,N5B,N4B)=0.0
            RQRNGS(N,1,N5B,N4B)=0.0
            RQRN2S(N,1,N5B,N4B)=0.0
            RQRHGS(N,1,N5B,N4B)=0.0
            RQRNH4(N,1,N5B,N4B)=0.0
            RQRNH3(N,1,N5B,N4B)=0.0
            RQRNO3(N,1,N5B,N4B)=0.0
            RQRNO2(N,1,N5B,N4B)=0.0
            RQRH1P(N,1,N5B,N4B)=0.0
            RQRH2P(N,1,N5B,N4B)=0.0
          ENDIF
        ENDIF
      ELSE
        DO K=0,jcplx1
          RQROC(K,N,2,N5,N4)=0.0
          RQRON(K,N,2,N5,N4)=0.0
          RQROP(K,N,2,N5,N4)=0.0
          RQROA(K,N,2,N5,N4)=0.0
        ENDDO
        RQRCOS(N,2,N5,N4)=0.0
        RQRCHS(N,2,N5,N4)=0.0
        RQROXS(N,2,N5,N4)=0.0
        RQRNGS(N,2,N5,N4)=0.0
        RQRN2S(N,2,N5,N4)=0.0
        RQRHGS(N,2,N5,N4)=0.0
        RQRNH4(N,2,N5,N4)=0.0
        RQRNH3(N,2,N5,N4)=0.0
        RQRNO3(N,2,N5,N4)=0.0
        RQRNO2(N,2,N5,N4)=0.0
        RQRH1P(N,2,N5,N4)=0.0
        RQRH2P(N,2,N5,N4)=0.0
        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          DO  K=0,jcplx1
            RQROC(K,N,1,N5B,N4B)=0.0
            RQRON(K,N,1,N5B,N4B)=0.0
            RQROP(K,N,1,N5B,N4B)=0.0
            RQROA(K,N,1,N5B,N4B)=0.0
          ENDDO
      RQRCOS(N,1,N5B,N4B)=0.0
      RQRCHS(N,1,N5B,N4B)=0.0
      RQROXS(N,1,N5B,N4B)=0.0
      RQRNGS(N,1,N5B,N4B)=0.0
      RQRN2S(N,1,N5B,N4B)=0.0
      RQRHGS(N,1,N5B,N4B)=0.0
      RQRNH4(N,1,N5B,N4B)=0.0
      RQRNH3(N,1,N5B,N4B)=0.0
      RQRNO3(N,1,N5B,N4B)=0.0
      RQRNO2(N,1,N5B,N4B)=0.0
      RQRH1P(N,1,N5B,N4B)=0.0
      RQRH2P(N,1,N5B,N4B)=0.0
      ENDIF
      ENDIF
!
!     SNOW DRIFT
!
!     subroutine SnowdriftTransport(M)
!
!     QSM=snow transfer from watsub.f
!     RQS*=solute flux in snow transfer
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     VOLS=volume of snowpack from watsub.f
!     *W2=solute content of snowpack
!
      IF(NN.EQ.1)THEN
!
!     IF NO SNOW DRIFT THEN NO TRANSPORT
!
      IF(ABS(QSM(M,N,N5,N4)).LE.ZEROS2(N2,N1))THEN
      RQSCOS(N,N5,N4)=0.0
      RQSCHS(N,N5,N4)=0.0
      RQSOXS(N,N5,N4)=0.0
      RQSNGS(N,N5,N4)=0.0
      RQSN2S(N,N5,N4)=0.0
      RQSNH4(N,N5,N4)=0.0
      RQSNH3(N,N5,N4)=0.0
      RQSNO3(N,N5,N4)=0.0
      RQSH1P(N,N5,N4)=0.0
      RQSH2P(N,N5,N4)=0.0
!
!     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
!
      ELSEIF(QSM(M,N,N5,N4).GT.ZEROS2(N2,N1))THEN
      IF(VOLS(N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMAX1(0.0,AMIN1(VFLWX,QSM(M,N,N5,N4)/VOLS(N2,N1)))
      ELSE
      VFLW=VFLWX
      ENDIF
      RQSCOS(N,N5,N4)=VFLW*AMAX1(0.0,CO2W2(1,N2,N1))
      RQSCHS(N,N5,N4)=VFLW*AMAX1(0.0,CH4W2(1,N2,N1))
      RQSOXS(N,N5,N4)=VFLW*AMAX1(0.0,OXYW2(1,N2,N1))
      RQSNGS(N,N5,N4)=VFLW*AMAX1(0.0,ZNGW2(1,N2,N1))
      RQSN2S(N,N5,N4)=VFLW*AMAX1(0.0,ZN2W2(1,N2,N1))
      RQSNH4(N,N5,N4)=VFLW*AMAX1(0.0,ZN4W2(1,N2,N1))
      RQSNH3(N,N5,N4)=VFLW*AMAX1(0.0,ZN3W2(1,N2,N1))
      RQSNO3(N,N5,N4)=VFLW*AMAX1(0.0,ZNOW2(1,N2,N1))
      RQSH1P(N,N5,N4)=VFLW*AMAX1(0.0,Z1PW2(1,N2,N1))
      RQSH2P(N,N5,N4)=VFLW*AMAX1(0.0,ZHPW2(1,N2,N1))
!
!     IF DRIFT IS TO CURRENT FROM ADJACENT GRID CELL
!
      ELSEIF(QSM(M,N,N5,N4).LT.-ZEROS2(N2,N1))THEN
      IF(VOLS(N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AMIN1(0.0,AMAX1(-VFLWX,QSM(M,N,N5,N4)/VOLS(N5,N4)))
      ELSE
      VFLW=-VFLWX
      ENDIF
      RQSCOS(N,N5,N4)=VFLW*AMAX1(0.0,CO2W2(1,N5,N4))
      RQSCHS(N,N5,N4)=VFLW*AMAX1(0.0,CH4W2(1,N5,N4))
      RQSOXS(N,N5,N4)=VFLW*AMAX1(0.0,OXYW2(1,N5,N4))
      RQSNGS(N,N5,N4)=VFLW*AMAX1(0.0,ZNGW2(1,N5,N4))
      RQSN2S(N,N5,N4)=VFLW*AMAX1(0.0,ZN2W2(1,N5,N4))
      RQSNH4(N,N5,N4)=VFLW*AMAX1(0.0,ZN4W2(1,N5,N4))
      RQSNH3(N,N5,N4)=VFLW*AMAX1(0.0,ZN3W2(1,N5,N4))
      RQSNO3(N,N5,N4)=VFLW*AMAX1(0.0,ZNOW2(1,N5,N4))
      RQSH1P(N,N5,N4)=VFLW*AMAX1(0.0,Z1PW2(1,N5,N4))
      RQSH2P(N,N5,N4)=VFLW*AMAX1(0.0,ZHPW2(1,N5,N4))
      ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QSS=hourly solute in snow flux
!     RQS*=solute in snow flux
!
      XCOQSS(N,N5,N4)=XCOQSS(N,N5,N4)+RQSCOS(N,N5,N4)
      XCHQSS(N,N5,N4)=XCHQSS(N,N5,N4)+RQSCHS(N,N5,N4)
      XOXQSS(N,N5,N4)=XOXQSS(N,N5,N4)+RQSOXS(N,N5,N4)
      XNGQSS(N,N5,N4)=XNGQSS(N,N5,N4)+RQSNGS(N,N5,N4)
      XN2QSS(N,N5,N4)=XN2QSS(N,N5,N4)+RQSN2S(N,N5,N4)
      XN4QSS(N,N5,N4)=XN4QSS(N,N5,N4)+RQSNH4(N,N5,N4)
      XN3QSS(N,N5,N4)=XN3QSS(N,N5,N4)+RQSNH3(N,N5,N4)
      XNOQSS(N,N5,N4)=XNOQSS(N,N5,N4)+RQSNO3(N,N5,N4)
      XP1QSS(N,N5,N4)=XP1QSS(N,N5,N4)+RQSH1P(N,N5,N4)
      XP4QSS(N,N5,N4)=XP4QSS(N,N5,N4)+RQSH2P(N,N5,N4)
!     IF(I.EQ.118.AND.NX.EQ.3.AND.NY.EQ.4)THEN
!     WRITE(*,6969)'XOXQSS',I,J,N4,N5,N,M,MM,XOXQSS(N,N5,N4)
!    2,RQSOXS(N,N5,N4),VFLW,OXYW2(N2,N1),OXYW2(N5,N4)
!    3,QSM(M,N,N5,N4),VOLS(N2,N1),VOLS(N5,N4)
!6969  FORMAT(A8,7I4,20E12.4)
!     ENDIF
      ENDIF
4305  CONTINUE
4310  CONTINUE
      end subroutine OverlandFlowSnowdriftTransport
!------------------------------------------------------------------------------------------

  subroutine SurfaceGasVolatilDissol(M,NY,NX)
  implicit none

  integer, intent(in) :: M,NY, NX
  real(r8) :: VOLHGR(JY,JX),VOLN3R(JY,JX)
  real(r8) :: VOLN2R(JY,JX),VOLOXR(JY,JX)
  real(r8) :: VOLCOR(JY,JX),VOLCHR(JY,JX)
  real(r8) :: VOLNGR(JY,JX)
  real(r8) :: CNH3S0,CNH4S0,CO2G0,CH4G0
  real(r8) :: Z2GG0,Z2OG0,ZN3G0,H2GG0,OXYG0
!
!     VOLT=litter volume from hour1.f
!     VOLWM,VOLPM=micropore water volume, air volume from watsub.f
!     VOLW*=equivalent aqueous volume for gas
!     S*L=solubility of gas in water from hour1.f
!     C*G=gas concentration
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     R*DFG=surface gas volatilization
!     DFGS=rate constant for air-water gas exchange from watsub.f
!     *S2=litter solute content
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     R*DXR=gas exchange between atmosphere and surface litter water for gas flux calculations
!
  IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VOLPM(M,0,NY,NX).GT.ZEROS2(NY,NX) &
    .AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VOLWCO(0,NY,NX)=VOLWM(M,0,NY,NX)*SCO2L(0,NY,NX)
    VOLWCH(0,NY,NX)=VOLWM(M,0,NY,NX)*SCH4L(0,NY,NX)
    VOLWOX(0,NY,NX)=VOLWM(M,0,NY,NX)*SOXYL(0,NY,NX)
    VOLWNG(0,NY,NX)=VOLWM(M,0,NY,NX)*SN2GL(0,NY,NX)
    VOLWN2(0,NY,NX)=VOLWM(M,0,NY,NX)*SN2OL(0,NY,NX)
    VOLWN3(0,NY,NX)=VOLWM(M,0,NY,NX)*SNH3L(0,NY,NX)
    VOLWHG(0,NY,NX)=VOLWM(M,0,NY,NX)*SH2GL(0,NY,NX)
    VOLWXA(0,NY,NX)=14.0*VOLWM(M,0,NY,NX)
    CO2G0=CCO2G(0,NY,NX)*VOLPM(M,0,NY,NX)
    CH4G0=CCH4G(0,NY,NX)*VOLPM(M,0,NY,NX)
    OXYG0=COXYG(0,NY,NX)*VOLPM(M,0,NY,NX)
    Z2GG0=CZ2GG(0,NY,NX)*VOLPM(M,0,NY,NX)
    Z2OG0=CZ2OG(0,NY,NX)*VOLPM(M,0,NY,NX)
    ZN3G0=CNH3G(0,NY,NX)*VOLPM(M,0,NY,NX)
    H2GG0=CH2GG(0,NY,NX)*VOLPM(M,0,NY,NX)
    VOLCOR(NY,NX)=VOLWCO(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLCHR(NY,NX)=VOLWCH(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLOXR(NY,NX)=VOLWOX(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLNGR(NY,NX)=VOLWNG(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLN2R(NY,NX)=VOLWN2(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLN3R(NY,NX)=VOLWN3(0,NY,NX)+VOLPM(M,0,NY,NX)
    VOLHGR(NY,NX)=VOLWHG(0,NY,NX)+VOLPM(M,0,NY,NX)
    RCODFG(0,NY,NX)=DFGS(M,0,NY,NX) &
      *(AMAX1(ZEROS(NY,NX),CO2G0)*VOLWCO(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),CO2S2(0,NY,NX)+RCODXR) &
      *VOLPM(M,0,NY,NX))/VOLCOR(NY,NX)
    RCHDFG(0,NY,NX)=DFGS(M,0,NY,NX) &
      *(AMAX1(ZEROS(NY,NX),CH4G0)*VOLWCH(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),CH4S2(0,NY,NX)+RCHDXR) &
      *VOLPM(M,0,NY,NX))/VOLCHR(NY,NX)
    ROXDFG(0,NY,NX)=DFGS(M,0,NY,NX) &
      *(AMAX1(ZEROS(NY,NX),OXYG0)*VOLWOX(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),OXYS2(0,NY,NX)+ROXDXR) &
      *VOLPM(M,0,NY,NX))/VOLOXR(NY,NX)
    RNGDFG(0,NY,NX)=DFGS(M,0,NY,NX) &
      *(AMAX1(ZEROS(NY,NX),Z2GG0)*VOLWNG(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),Z2GS2(0,NY,NX)+RNGDXR) &
      *VOLPM(M,0,NY,NX))/VOLNGR(NY,NX)
    RN2DFG(0,NY,NX)=DFGS(M,0,NY,NX) &
      *(AMAX1(ZEROS(NY,NX),Z2OG0)*VOLWN2(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),Z2OS2(0,NY,NX)+RN2DXR) &
      *VOLPM(M,0,NY,NX))/VOLN2R(NY,NX)
    RN3DFG(0,NY,NX)=DFGS(M,0,NY,NX) &
      *(AMAX1(ZEROS(NY,NX),ZN3G0)*VOLWN3(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),ZNH3S2(0,NY,NX)+RN3DXR) &
      *VOLPM(M,0,NY,NX))/VOLN3R(NY,NX)
    CNH3S0=AMAX1(0.0,(ZNH3S2(0,NY,NX)+RN3DFG(0,NY,NX)))/VOLWXA(0,NY,NX)
    CNH4S0=AMAX1(0.0,ZNH4S2(0,NY,NX))/VOLWXA(0,NY,NX)
    RHGDFG(0,NY,NX)=DFGS(M,0,NY,NX) &
      *(AMAX1(ZEROS(NY,NX),H2GG0)*VOLWHG(0,NY,NX) &
      -AMAX1(ZEROS(NY,NX),H2GS2(0,NY,NX)+RHGDXR) &
      *VOLPM(M,0,NY,NX))/VOLHGR(NY,NX)

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly surface gas volatilization
!     R*DFG=surface gas volatilization
!
    XCODFG(0,NY,NX)=XCODFG(0,NY,NX)+RCODFG(0,NY,NX)
    XCHDFG(0,NY,NX)=XCHDFG(0,NY,NX)+RCHDFG(0,NY,NX)
    XOXDFG(0,NY,NX)=XOXDFG(0,NY,NX)+ROXDFG(0,NY,NX)
    XNGDFG(0,NY,NX)=XNGDFG(0,NY,NX)+RNGDFG(0,NY,NX)
    XN2DFG(0,NY,NX)=XN2DFG(0,NY,NX)+RN2DFG(0,NY,NX)
    XN3DFG(0,NY,NX)=XN3DFG(0,NY,NX)+RN3DFG(0,NY,NX)
    XHGDFG(0,NY,NX)=XHGDFG(0,NY,NX)+RHGDFG(0,NY,NX)
  ELSE
    RCODFG(0,NY,NX)=0.0
    RCHDFG(0,NY,NX)=0.0
    ROXDFG(0,NY,NX)=0.0
    RNGDFG(0,NY,NX)=0.0
    RN2DFG(0,NY,NX)=0.0
    RN3DFG(0,NY,NX)=0.0
    RHGDFG(0,NY,NX)=0.0
  ENDIF
  end subroutine SurfaceGasVolatilDissol
!------------------------------------------------------------------------------------------

      subroutine GasDiffusionConvection(M,NY,NX,FLQM)
      implicit none

      integer, intent(in) :: M, NY, NX
  real(r8),intent(in) :: FLQM(3,JD,JV,JH)
  real(r8) :: VFLW
  real(r8) :: VOLHGT(JY,JX),VOLNBT(JY,JX)
  real(r8) :: VOLN3T(JY,JX),VOLN2T(JY,JX)
  real(r8) :: VOLCOT(JY,JX),VOLCHT(JY,JX)
  real(r8) :: VOLOXT(JY,JX),VOLNGT(JY,JX)
  real(r8) :: CNH3S0,CCH4G2,CH2GG2
  real(r8) :: CNH3B0,CNH4B0,CNH4S0
  real(r8) :: CCO2G2,COXYG2,CZ2GG2,CZ2OG2,CNH3G2
  real(r8) :: DFVHGG,DFVN3G,DFVN2G,DFLG2
  real(r8) :: DNH3GQ,DH2GGQ
  real(r8) :: DCO2GQ,DCH4GQ,DOXYGQ,DZ2GGQ,DZ2OGQ

!
!     THETPM=air-filled porosity from watsub.f
!     BKDS=bulk density
!
      IF(THETPM(M,NU(NY,NX),NY,NX).GT.THETX &
      .AND.BKDS(NU(NY,NX),NY,NX).GT.ZERO)THEN
!
!     GASEOUS DIFFUSIVITIES
!
!     DFLG2=air-filled porosity effect on gaseous diffusivity
!     POROQ=Penman Water Linear Reduction tortuosity from starts.f
!     POROS=total porosity
!     DLYR=soil layer thickness
!     D*G=gaseous diffusivity in soil
!     *SGL2= gaseous diffusivity in air
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
      DFLG2=AMAX1(0.0,THETPM(M,NU(NY,NX),NY,NX))*POROQ &
      *THETPM(M,NU(NY,NX),NY,NX)/POROS(NU(NY,NX),NY,NX) &
      *AREA(3,NU(NY,NX),NY,NX)/AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
      DCO2G(3,NU(NY,NX),NY,NX)=DFLG2*CGSGL2(NU(NY,NX),NY,NX)
      DCH4G(3,NU(NY,NX),NY,NX)=DFLG2*CHSGL2(NU(NY,NX),NY,NX)
      DOXYG(3,NU(NY,NX),NY,NX)=DFLG2*OGSGL2(NU(NY,NX),NY,NX)
      DZ2GG(3,NU(NY,NX),NY,NX)=DFLG2*ZGSGL2(NU(NY,NX),NY,NX)
      DZ2OG(3,NU(NY,NX),NY,NX)=DFLG2*Z2SGL2(NU(NY,NX),NY,NX)
      DNH3G(3,NU(NY,NX),NY,NX)=DFLG2*ZHSGL2(NU(NY,NX),NY,NX)
      DH2GG(3,NU(NY,NX),NY,NX)=DFLG2*HGSGL2(NU(NY,NX),NY,NX)
!
!     SURFACE GAS CONCENTRATIONS
!
!     C*G2=gaseous concentration
!     *G2=gaseous content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     VOLPM=air-filled porosity
!
      CCO2G2=AMAX1(0.0,CO2G2(NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX))
      CCH4G2=AMAX1(0.0,CH4G2(NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX))
      COXYG2=AMAX1(0.0,OXYG2(NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX))
      CZ2GG2=AMAX1(0.0,Z2GG2(NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX))
      CZ2OG2=AMAX1(0.0,Z2OG2(NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX))
      CNH3G2=AMAX1(0.0,ZN3G2(NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX))
      CH2GG2=AMAX1(0.0,H2GG2(NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX))
!
!     EQUILIBRIUM CONCENTRATIONS AT SOIL SURFACE AT WHICH
!     GASEOUS DIFFUSION THROUGH SOIL SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM GASEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
!     D*GQ=gaseous diffusivity between soil surface and atmosphere
!     D*G=gaseous diffusivity in soil
!     PARG*=boundary layer conductance above soil surface
!     DFV*G=diffusive gas flux
!     C*E=atmospheric gas concentration from hour1.f
!     C*G2=gaseous concentration
!
      DCO2GQ=DCO2G(3,NU(NY,NX),NY,NX)*PARGCO(NY,NX) &
      /(DCO2G(3,NU(NY,NX),NY,NX)+PARGCO(NY,NX))
      DCH4GQ=DCH4G(3,NU(NY,NX),NY,NX)*PARGCH(NY,NX) &
      /(DCH4G(3,NU(NY,NX),NY,NX)+PARGCH(NY,NX))
      DOXYGQ=DOXYG(3,NU(NY,NX),NY,NX)*PARGOX(NY,NX) &
      /(DOXYG(3,NU(NY,NX),NY,NX)+PARGOX(NY,NX))
      DZ2GGQ=DZ2GG(3,NU(NY,NX),NY,NX)*PARGNG(NY,NX) &
      /(DZ2GG(3,NU(NY,NX),NY,NX)+PARGNG(NY,NX))
      DZ2OGQ=DZ2OG(3,NU(NY,NX),NY,NX)*PARGN2(NY,NX) &
      /(DZ2OG(3,NU(NY,NX),NY,NX)+PARGN2(NY,NX))
      DNH3GQ=DNH3G(3,NU(NY,NX),NY,NX)*PARGN3(NY,NX) &
      /(DNH3G(3,NU(NY,NX),NY,NX)+PARGN3(NY,NX))
      DH2GGQ=DH2GG(3,NU(NY,NX),NY,NX)*PARGH2(NY,NX) &
      /(DH2GG(3,NU(NY,NX),NY,NX)+PARGH2(NY,NX))
      DFVCOG=DCO2GQ*(CCO2E(NY,NX)-CCO2G2)
      DFVCHG=DCH4GQ*(CCH4E(NY,NX)-CCH4G2)
      DFVOXG=DOXYGQ*(COXYE(NY,NX)-COXYG2)
      DFVNGG=DZ2GGQ*(CZ2GE(NY,NX)-CZ2GG2)
      DFVN2G=DZ2OGQ*(CZ2OE(NY,NX)-CZ2OG2)
      DFVN3G=DNH3GQ*(CNH3E(NY,NX)-CNH3G2)
      DFVHGG=DH2GGQ*(CH2GE(NY,NX)-CH2GG2)
!
!     CONVECTIVE GAS TRANSFER DRIVEN BY SURFACE WATER FLUXES
!     FROM 'WATSUB' AND GAS CONCENTRATIONS IN THE SOIL SURFACE
!     OR THE ATMOSPHERE DEPENDING ON WATER FLUX DIRECTION
!
!     FLQM=total water flux into soil micropore+macropore from watsub.f
!     VOLPM=air-filled porosity
!     RFL*G=convective gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     *G2=gaseous content
!     C*E=atmospheric gas concentration from hour1.f
!
      IF(FLQM(3,NU(NY,NX),NY,NX).GT.0.0)THEN
      IF(VOLPM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=-AMAX1(0.0,AMIN1(VFLWX,FLQM(3,NU(NY,NX),NY,NX) &
      /VOLPM(M,NU(NY,NX),NY,NX)))
      ELSE
      VFLW=-VFLWX
      ENDIF
      RFLCOG=VFLW*AMAX1(0.0,CO2G2(NU(NY,NX),NY,NX))
      RFLCHG=VFLW*AMAX1(0.0,CH4G2(NU(NY,NX),NY,NX))
      RFLOXG=VFLW*AMAX1(0.0,OXYG2(NU(NY,NX),NY,NX))
      RFLNGG=VFLW*AMAX1(0.0,Z2GG2(NU(NY,NX),NY,NX))
      RFLN2G=VFLW*AMAX1(0.0,Z2OG2(NU(NY,NX),NY,NX))
      RFLN3G=VFLW*AMAX1(0.0,ZN3G2(NU(NY,NX),NY,NX))
      RFLH2G=VFLW*AMAX1(0.0,H2GG2(NU(NY,NX),NY,NX))
      ELSE
      RFLCOG=-FLQM(3,NU(NY,NX),NY,NX)*CCO2E(NY,NX)
      RFLCHG=-FLQM(3,NU(NY,NX),NY,NX)*CCH4E(NY,NX)
      RFLOXG=-FLQM(3,NU(NY,NX),NY,NX)*COXYE(NY,NX)
      RFLNGG=-FLQM(3,NU(NY,NX),NY,NX)*CZ2GE(NY,NX)
      RFLN2G=-FLQM(3,NU(NY,NX),NY,NX)*CZ2OE(NY,NX)
      RFLN3G=-FLQM(3,NU(NY,NX),NY,NX)*CNH3E(NY,NX)
      RFLH2G=-FLQM(3,NU(NY,NX),NY,NX)*CH2GE(NY,NX)
      ENDIF
!
!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE + CONVECTIVE FLUX
!
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux
!
      RCOFLG(3,NU(NY,NX),NY,NX)=DFVCOG+RFLCOG
      RCHFLG(3,NU(NY,NX),NY,NX)=DFVCHG+RFLCHG
      ROXFLG(3,NU(NY,NX),NY,NX)=DFVOXG+RFLOXG
      RNGFLG(3,NU(NY,NX),NY,NX)=DFVNGG+RFLNGG
      RN2FLG(3,NU(NY,NX),NY,NX)=DFVN2G+RFLN2G
      RN3FLG(3,NU(NY,NX),NY,NX)=DFVN3G+RFLN3G
      RHGFLG(3,NU(NY,NX),NY,NX)=DFVHGG+RFLH2G
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLG=hourly convective+diffusive gas flux
!
      XCOFLG(3,NU(NY,NX),NY,NX)=XCOFLG(3,NU(NY,NX),NY,NX) &
      +RCOFLG(3,NU(NY,NX),NY,NX)
      XCHFLG(3,NU(NY,NX),NY,NX)=XCHFLG(3,NU(NY,NX),NY,NX) &
      +RCHFLG(3,NU(NY,NX),NY,NX)
      XOXFLG(3,NU(NY,NX),NY,NX)=XOXFLG(3,NU(NY,NX),NY,NX) &
      +ROXFLG(3,NU(NY,NX),NY,NX)
      XNGFLG(3,NU(NY,NX),NY,NX)=XNGFLG(3,NU(NY,NX),NY,NX) &
      +RNGFLG(3,NU(NY,NX),NY,NX)
      XN2FLG(3,NU(NY,NX),NY,NX)=XN2FLG(3,NU(NY,NX),NY,NX) &
      +RN2FLG(3,NU(NY,NX),NY,NX)
      XN3FLG(3,NU(NY,NX),NY,NX)=XN3FLG(3,NU(NY,NX),NY,NX) &
      +RN3FLG(3,NU(NY,NX),NY,NX)
      XHGFLG(3,NU(NY,NX),NY,NX)=XHGFLG(3,NU(NY,NX),NY,NX) &
      +RHGFLG(3,NU(NY,NX),NY,NX)
!     IF(J.EQ.24)THEN
!     WRITE(*,3131)'ROXFLG',I,J,NX,NY,M,MM,XOXFLG(3,NU(NY,NX),NY,NX)
!    2,ROXFLG(3,NU(NY,NX),NY,NX),DFVOXG,RFLOXG,COXYE(NY,NX)
!    2,COXYG2,DOXYGQ,OXYG2(NU(NY,NX),NY,NX),FLQM(3,NU(NY,NX),NY,NX)
!    3,VFLW,DOXYG(3,NU(NY,NX),NY,NX),PARGOX(NY,NX)
!    4,THETPM(M,NU(NY,NX),NY,NX),VOLPM(M,NU(NY,NX),NY,NX)
!    5,DFGS(M,NU(NY,NX),NY,NX),FLWM(M,3,NU(NY,NX),NY,NX)
!    2+FLWHM(M,3,NU(NY,NX),NY,NX)
!     WRITE(*,3131)'RNGFLG',I,J,NX,NY,M,MM,XNGFLG(3,NU(NY,NX),NY,NX)
!    2,RNGFLG(3,NU(NY,NX),NY,NX),DFVNGG,RFLNGG,CZ2GE(NY,NX)
!    2,CZ2GG2,DZ2GGQ,Z2GG2(NU(NY,NX),NY,NX),FLQM(3,NU(NY,NX),NY,NX)
!    3,VFLW,DZ2GG(3,NU(NY,NX),NY,NX),PARGNG(NY,NX)
!    4,THETPM(M,NU(NY,NX),NY,NX),VOLPM(M,NU(NY,NX),NY,NX)
!3131  FORMAT(A8,6I4,30E12.4)
!     ENDIF
!     VOLWM=micropore water-filled porosity from watsub.f
!     VOLW*=equivalent aqueous volume for gas
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     S*L=solubility of gas in water from hour1.f
!     VOLPM=air-filled porosity
!     R*DFG=water-air gas flux
!     DFGS=rate constant for air-water gas exchange from watsub.f
!     *G2,*S2=gaseous,aqueous gas content
!     R*DXS=gas exchange between atmosphere and soil surface water
!
      IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLWCO(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX) &
      *SCO2L(NU(NY,NX),NY,NX)
      VOLWCH(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX) &
      *SCH4L(NU(NY,NX),NY,NX)
      VOLWOX(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX) &
      *SOXYL(NU(NY,NX),NY,NX)
      VOLWNG(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX) &
      *SN2GL(NU(NY,NX),NY,NX)
      VOLWN2(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX) &
      *SN2OL(NU(NY,NX),NY,NX)
      VOLWN3(NU(NY,NX),NY,NX)=VOLWMA(NU(NY,NX),NY,NX) &
      *SNH3L(NU(NY,NX),NY,NX)
      VOLWNB(NU(NY,NX),NY,NX)=VOLWMB(NU(NY,NX),NY,NX) &
      *SNH3L(NU(NY,NX),NY,NX)
      VOLWHG(NU(NY,NX),NY,NX)=VOLWM(M,NU(NY,NX),NY,NX) &
      *SH2GL(NU(NY,NX),NY,NX)
      VOLCOT(NY,NX)=VOLWCO(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
      VOLCHT(NY,NX)=VOLWCH(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
      VOLOXT(NY,NX)=VOLWOX(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
      VOLNGT(NY,NX)=VOLWNG(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
      VOLN2T(NY,NX)=VOLWN2(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
      VOLN3T(NY,NX)=VOLWN3(NU(NY,NX),NY,NX)+VOLPMA(NU(NY,NX),NY,NX)
      VOLNBT(NY,NX)=VOLWNB(NU(NY,NX),NY,NX)+VOLPMB(NU(NY,NX),NY,NX)
      VOLHGT(NY,NX)=VOLWHG(NU(NY,NX),NY,NX)+VOLPM(M,NU(NY,NX),NY,NX)
      RCODFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),CO2G2(NU(NY,NX),NY,NX)) &
      *VOLWCO(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,CO2S2(NU(NY,NX),NY,NX)+RCODXS) &
      *VOLPM(M,NU(NY,NX),NY,NX))/VOLCOT(NY,NX)
      RCHDFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),CH4G2(NU(NY,NX),NY,NX)) &
      *VOLWCH(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,CH4S2(NU(NY,NX),NY,NX)+RCHDXS) &
      *VOLPM(M,NU(NY,NX),NY,NX))/VOLCHT(NY,NX)
      ROXDFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),OXYG2(NU(NY,NX),NY,NX)) &
      *VOLWOX(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,OXYS2(NU(NY,NX),NY,NX)+ROXDXS) &
      *VOLPM(M,NU(NY,NX),NY,NX))/VOLOXT(NY,NX)
      RNGDFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),Z2GG2(NU(NY,NX),NY,NX)) &
      *VOLWNG(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,Z2GS2(NU(NY,NX),NY,NX)+RNGDXS) &
      *VOLPM(M,NU(NY,NX),NY,NX))/VOLNGT(NY,NX)
      RN2DFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),Z2OG2(NU(NY,NX),NY,NX)) &
      *VOLWN2(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,Z2OS2(NU(NY,NX),NY,NX)+RN2DXS) &
      *VOLPM(M,NU(NY,NX),NY,NX))/VOLN2T(NY,NX)
      IF(VOLN3T(NY,NX).GT.ZEROS2(NY,NX) &
      .AND.VOLWXA(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      RN3DFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),ZN3G2(NU(NY,NX),NY,NX)) &
      *VOLWN3(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,ZNH3S2(NU(NY,NX),NY,NX)+RN3DXS) &
      *VOLPMA(NU(NY,NX),NY,NX))/VOLN3T(NY,NX)
      CNH3S0=AMAX1(0.0,(ZNH3S2(NU(NY,NX),NY,NX) &
      +RN3DFG(NU(NY,NX),NY,NX))/VOLWXA(NU(NY,NX),NY,NX))
      CNH4S0=AMAX1(0.0,ZNH4S2(NU(NY,NX),NY,NX)) &
      /VOLWXA(NU(NY,NX),NY,NX)
      ELSE
      RN3DFG(NU(NY,NX),NY,NX)=0.0
      ENDIF
      IF(VOLNBT(NY,NX).GT.ZEROS2(NY,NX) &
      .AND.VOLWXB(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      RNBDFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),ZN3G2(NU(NY,NX),NY,NX)) &
      *VOLWNB(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,ZNH3B2(NU(NY,NX),NY,NX)+RNBDXS) &
      *VOLPMB(NU(NY,NX),NY,NX))/VOLNBT(NY,NX)
      CNH3B0=AMAX1(0.0,(ZNH3B2(NU(NY,NX),NY,NX) &
      +RNBDFG(NU(NY,NX),NY,NX))/VOLWXB(NU(NY,NX),NY,NX))
      CNH4B0=AMAX1(0.0,ZNH4B2(NU(NY,NX),NY,NX)) &
      /VOLWXB(NU(NY,NX),NY,NX)
      ELSE
      RNBDFG(NU(NY,NX),NY,NX)=0.0
      ENDIF
      RHGDFG(NU(NY,NX),NY,NX)=DFGS(M,NU(NY,NX),NY,NX) &
      *(AMAX1(ZEROS(NY,NX),H2GG2(NU(NY,NX),NY,NX)) &
      *VOLWHG(NU(NY,NX),NY,NX)-AMAX1(ZEROS(NY,NX) &
      ,H2GS2(NU(NY,NX),NY,NX)+RHGDXS) &
      *VOLPM(M,NU(NY,NX),NY,NX))/VOLHGT(NY,NX)
!     IF(I.EQ.121)THEN
!     WRITE(*,323)'RCODFG',I,J,M,MM,NX,NY,XCODFG(NU(NY,NX),NY,NX)
!    2,RCODFG(NU(NY,NX),NY,NX),DFGS(M,NU(NY,NX),NY,NX)
!    2,CO2G2(NU(NY,NX),NY,NX),VOLWCO(NU(NY,NX),NY,NX)
!    2,CO2S2(NU(NY,NX),NY,NX),RCODXS,VOLWM(M,NU(NY,NX),NY,NX)
!    4,VOLPM(M,NU(NY,NX),NY,NX),VOLCOT(NY,NX),SCO2L(NU(NY,NX),NY,NX)
!     5,RCOFLG(3,NU(NY,NX),NY,NX),DFVCOG,RFLCOG
!     WRITE(*,323)'ROXDFG',I,J,M,MM,NX,NY,XOXDFG(NU(NY,NX),NY,NX)
!    2,ROXDFG(NU(NY,NX),NY,NX),DFGS(M,NU(NY,NX),NY,NX)
!    2,OXYG2(NU(NY,NX),NY,NX),VOLWOX(NU(NY,NX),NY,NX)
!    2,OXYS2(NU(NY,NX),NY,NX),ROXDXS,VOLWM(M,NU(NY,NX),NY,NX)
!    4,VOLPM(M,NU(NY,NX),NY,NX),VOLOXT(NY,NX),SOXYL(NU(NY,NX),NY,NX)
!    5,ROXFLG(3,NU(NY,NX),NY,NX),DFVOXG,RFLOXG
!     WRITE(*,323)'RN3FLG',I,J,M,MM,NX,NY,RN3FLG(3,NU(NY,NX),NY,NX)
!    2,DNH3GQ,CNH3E(NY,NX),CNH3G2,FLQM(3,NU(NY,NX),NY,NX),CNH3GV
!    2,CNH3B2,ZNH3B2(NU(NY,NX),NY,NX),RNBDFG(NU(NY,NX),NY,NX)
!    3,DFGS(M,NU(NY,NX),NY,NX),ZN3G2B,VOLPMB(NU(NY,NX),NY,NX)
!    4,ZNH3B2(NU(NY,NX),NY,NX),VOLWNB(NU(NY,NX),NY,NX)
!    5,VOLWMB,SNH3L(NU(NY,NX),NY,NX)
!     WRITE(*,323)'RNGDFG',I,J,M,MM,NX,NY,RNGDFG(NU(NY,NX),NY,NX)
!    2,DFGS(M,NU(NY,NX),NY,NX),Z2GG2(NU(NY,NX),NY,NX)
!    3,VOLWNG(NU(NY,NX),NY,NX),Z2GS2(NU(NY,NX),NY,NX)
!    4,RNGDFS(NY,NX),VOLPM(M,NU(NY,NX),NY,NX),VOLNGT(NY,NX)
!     ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly water-air gas flux
!
      XCODFG(NU(NY,NX),NY,NX)=XCODFG(NU(NY,NX),NY,NX) &
      +RCODFG(NU(NY,NX),NY,NX)
      XCHDFG(NU(NY,NX),NY,NX)=XCHDFG(NU(NY,NX),NY,NX) &
      +RCHDFG(NU(NY,NX),NY,NX)
      XOXDFG(NU(NY,NX),NY,NX)=XOXDFG(NU(NY,NX),NY,NX) &
      +ROXDFG(NU(NY,NX),NY,NX)
      XNGDFG(NU(NY,NX),NY,NX)=XNGDFG(NU(NY,NX),NY,NX) &
      +RNGDFG(NU(NY,NX),NY,NX)
      XN2DFG(NU(NY,NX),NY,NX)=XN2DFG(NU(NY,NX),NY,NX) &
      +RN2DFG(NU(NY,NX),NY,NX)
      XN3DFG(NU(NY,NX),NY,NX)=XN3DFG(NU(NY,NX),NY,NX) &
      +RN3DFG(NU(NY,NX),NY,NX)
      XNBDFG(NU(NY,NX),NY,NX)=XNBDFG(NU(NY,NX),NY,NX) &
      +RNBDFG(NU(NY,NX),NY,NX)
      XHGDFG(NU(NY,NX),NY,NX)=XHGDFG(NU(NY,NX),NY,NX) &
      +RHGDFG(NU(NY,NX),NY,NX)
!     WRITE(*,3131)'ROXDFG',I,J,NX,NY,M,MM,XOXDFG(NU(NY,NX),NY,NX)
!    2,ROXDFG(NU(NY,NX),NY,NX),DFGS(M,NU(NY,NX),NY,NX)
!    2,AMAX1(ZEROS(NY,NX),OXYG2(NU(NY,NX),NY,NX))
!    3,VOLWOX(NU(NY,NX),NY,NX),AMAX1(ZEROS(NY,NX)
!    4,OXYS2(NU(NY,NX),NY,NX)),VOLPM(M,NU(NY,NX),NY,NX)
      ELSE
      RCODFG(NU(NY,NX),NY,NX)=0.0
      RCHDFG(NU(NY,NX),NY,NX)=0.0
      ROXDFG(NU(NY,NX),NY,NX)=0.0
      RNGDFG(NU(NY,NX),NY,NX)=0.0
      RN2DFG(NU(NY,NX),NY,NX)=0.0
      RN3DFG(NU(NY,NX),NY,NX)=0.0
      RNBDFG(NU(NY,NX),NY,NX)=0.0
      RHGDFG(NU(NY,NX),NY,NX)=0.0
      ENDIF
      ELSE
      RCOFLG(3,NU(NY,NX),NY,NX)=0.0
      RCHFLG(3,NU(NY,NX),NY,NX)=0.0
      ROXFLG(3,NU(NY,NX),NY,NX)=0.0
      RNGFLG(3,NU(NY,NX),NY,NX)=0.0
      RN2FLG(3,NU(NY,NX),NY,NX)=0.0
      RN3FLG(3,NU(NY,NX),NY,NX)=0.0
      RHGFLG(3,NU(NY,NX),NY,NX)=0.0
      RCODFG(NU(NY,NX),NY,NX)=0.0
      RCHDFG(NU(NY,NX),NY,NX)=0.0
      ROXDFG(NU(NY,NX),NY,NX)=0.0
      RN2DFG(NU(NY,NX),NY,NX)=0.0
      RNGDFG(NU(NY,NX),NY,NX)=0.0
      RN3DFG(NU(NY,NX),NY,NX)=0.0
      RNBDFG(NU(NY,NX),NY,NX)=0.0
      RHGDFG(NU(NY,NX),NY,NX)=0.0
      ENDIF
      end subroutine GasDiffusionConvection
!------------------------------------------------------------------------------------------

      subroutine SoluteFluxAdjacentCells(M,MX,NY,NX,NHE,NVS)
      implicit none

      integer, intent(in) :: M,MX, NY, NX, NHE, NVS
  real(r8) :: FLQM(3,JD,JV,JH)
  real(r8) :: VFLW,VOLWPA,VOLWOB
  real(r8) :: VOLWOA,VOLWPB
  real(r8) :: CNH3S0,CCH4G2,CCH4S1,CCO2S1
  real(r8) :: CCH4SH1,COXYSH1,CZ2OSH1,CZ2GSH1,CCH4SH2
  real(r8) :: CPO4BH1,CCO2SH2,COXYSH2,CZ2GSH2,CZ2OSH2
  real(r8) :: CH2GSH1,CNH4SH1,CNH3SH1,CNO3SH1,CNO2SH1,CP14SH1
  real(r8) :: CPO4B1,CCO2SH1,CCH4G1,CH2GG1,CH2GG2
  real(r8) :: CPO4SH1,CNH4BH1,CNH3BH1,CNO3BH1,CNO2BH1,CP14BH1
  real(r8) :: CH2GSH2,CNH4SH2,CNH3SH2,CNO3SH2,CNO2SH2,CP14SH2
  real(r8) :: CPO4SH2,CNH4BH2,CNH3BH2,CNO3BH2,CNO2BH2,CP14BH2
  real(r8) :: CND21,CND22,CND41,CNDO1,CNDG1,CNDH1,CNHG1,CNDC1
  real(r8) :: CNDC2,CND42,CNDO2,CNDG2,CNDH2,CNHG2
  real(r8) :: CNH3B0,CNH4B0,CNH4S0
  real(r8) :: CNH4B1,CNH3B1,CNO3B1,CNO2B1,CP14B1
  real(r8) :: COXYG1,CZ2GG1,CZ2OG1,CNH3G1
  real(r8) :: CCO2G2,COXYG2,CZ2GG2,CZ2OG2,CNH3G2
  real(r8) :: THETW1(JZ,JY,JX)
  real(r8) :: DLYR1,CCO2G1,DLYR2,DFHNGS,DFVHGG
  real(r8) :: CPO4BH2,DFHCHS,DFHCOS,DFHOXS,DFHHGS,DFVN3G
  real(r8) :: DFVN2G,DFHN4B,DFHN3B,DFHNOB,DFHN2B,DFHP1B,DFHPOB
  real(r8) :: DFHN2S,DFLG2,DFLGL,DVTGAS,DIFCQ,DIFCS,DIFHG,DIFN2
  real(r8) :: DFHNH4,DFHNH3,DFHNO3,DFHNO2,DFHP14,DFHPO4,DISPN
  real(r8) :: DIFNG,DIFOS,DIFOC,DIFON,DIFOP,DIFOA,DIFNH,DIFNO,DIFPO
  real(r8) :: SCH4X,SOXYX,SN2GX,SN2OX,SNH3X,SH2GX,SCO2X
  real(r8) :: TORTL,VH2GG2,VTATM,VTGAS
  real(r8) :: VOLH2B,VOLHOA,VOLHOB,VOLHPA,VOLHPB,VOLHMA,VOLHMB
  real(r8) :: VCO2G2,VCH4G2,VOXYG2,VZ2GG2,VZ2OG2,VNH3G2,VNHBG2
  real(r8) :: VOLW3A,VOLW3B,VOLH3A,VOLH3B,VOLW2A,VOLW2B,VOLH2A
  real(r8) :: VOLW4A,VOLW4B,VOLH4A,VOLH4B
  real(r8) :: VOLWHS,VOLWT

  integer :: IFLGB,N,L,K,LL
  integer :: N1,N2,N3,N4,N5,N6

! begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
      IFLGB=0
      DO 125 L=1,NL(NY,NX)
      N1=NX
      N2=NY
      N3=L
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
      DO 120 N=NCN(N2,N1),3
      IF(N.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      GO TO 120
      ELSE
      N4=NX+1
      N5=NY
      N6=L
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NY.EQ.NVS)THEN
      GO TO 120
      ELSE
      N4=NX
      N5=NY+1
      N6=L
      ENDIF
      ELSEIF(N.EQ.3)THEN
      IF(L.EQ.NL(NY,NX))THEN
      GO TO 120
      ELSE
      N4=NX
      N5=NY
      N6=L+1
      ENDIF
      ENDIF
      DO 1100 LL=N6,NL(NY,NX)
      IF(VOLX(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
      N6=LL
      GO TO 1101
      ENDIF
1100  CONTINUE
1101  CONTINUE
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
!     WATER CONTENTS AND WATER FLUXES 'FLQM' FROM 'WATSUB'
!
!     VOLX,VOLY=soil volume excluding rock, macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     VOLWM,VOLWHM=micropore,macropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!     FLPM=change in air volume
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4) &
      .AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
      IF(M.NE.MX)THEN
      VOLW4A=VOLWM(M,N3,N2,N1)*VLNH4(N3,N2,N1)
      VOLW4B=VOLWM(M,N3,N2,N1)*VLNHB(N3,N2,N1)
      VOLH4A=VOLWHM(M,N3,N2,N1)*VLNH4(N3,N2,N1)
      VOLH4B=VOLWHM(M,N3,N2,N1)*VLNHB(N3,N2,N1)
      VOLW3A=VOLWM(M,N3,N2,N1)*VLNO3(N3,N2,N1)
      VOLW3B=VOLWM(M,N3,N2,N1)*VLNOB(N3,N2,N1)
      VOLH3A=VOLWHM(M,N3,N2,N1)*VLNO3(N3,N2,N1)
      VOLH3B=VOLWHM(M,N3,N2,N1)*VLNOB(N3,N2,N1)
      VOLW2A=VOLWM(M,N3,N2,N1)*VLPO4(N3,N2,N1)
      VOLW2B=VOLWM(M,N3,N2,N1)*VLPOB(N3,N2,N1)
      VOLH2A=VOLWHM(M,N3,N2,N1)*VLPO4(N3,N2,N1)
      VOLH2B=VOLWHM(M,N3,N2,N1)*VLPOB(N3,N2,N1)
      VOLWMA(N6,N5,N4)=VOLWM(M,N6,N5,N4)*VLNH4(N6,N5,N4)
      VOLWMB(N6,N5,N4)=VOLWM(M,N6,N5,N4)*VLNHB(N6,N5,N4)
      VOLWXA(N6,N5,N4)=14.0*VOLWMA(N6,N5,N4)
      VOLWXB(N6,N5,N4)=14.0*VOLWMB(N6,N5,N4)
      VOLWOA=VOLWM(M,N6,N5,N4)*VLNO3(N6,N5,N4)
      VOLWOB=VOLWM(M,N6,N5,N4)*VLNOB(N6,N5,N4)
      VOLHOA=VOLWHM(M,N6,N5,N4)*VLNO3(N6,N5,N4)
      VOLHOB=VOLWHM(M,N6,N5,N4)*VLNOB(N6,N5,N4)
      VOLWPA=VOLWM(M,N6,N5,N4)*VLPO4(N6,N5,N4)
      VOLWPB=VOLWM(M,N6,N5,N4)*VLPOB(N6,N5,N4)
      VOLHPA=VOLWHM(M,N6,N5,N4)*VLPO4(N6,N5,N4)
      VOLHPB=VOLWHM(M,N6,N5,N4)*VLPOB(N6,N5,N4)
      VOLPMA(N6,N5,N4)=VOLPM(M,N6,N5,N4)*VLNH4(N6,N5,N4)
      VOLPMB(N6,N5,N4)=VOLPM(M,N6,N5,N4)*VLNHB(N6,N5,N4)
      THETW1(N3,N2,N1)=AMAX1(0.0,safe_adb(VOLWM(M,N3,N2,N1),VOLY(N3,N2,N1)))
      THETW1(N6,N5,N4)=AMAX1(0.0,safe_adb(VOLWM(M,N6,N5,N4),VOLY(N6,N5,N4)))
      FLVM(N6,N5,N4)=FLPM(M,N6,N5,N4)*XNPT
!
!     GASEOUS SOLUBILITIES
!
!     VOLWM=micropore water-filled porosity from watsub.f
!     VOLW*=equivalent aqueous volume for gas
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     S*L=solubility of gas in water from hour1.f
!     FLQM=total water flux into soil micropore+macropore from watsub.f
!
      IF(N.EQ.3)THEN
      VOLWCO(N6,N5,N4)=VOLWM(M,N6,N5,N4)*SCO2L(N6,N5,N4)
      VOLWCH(N6,N5,N4)=VOLWM(M,N6,N5,N4)*SCH4L(N6,N5,N4)
      VOLWOX(N6,N5,N4)=VOLWM(M,N6,N5,N4)*SOXYL(N6,N5,N4)
      VOLWNG(N6,N5,N4)=VOLWM(M,N6,N5,N4)*SN2GL(N6,N5,N4)
      VOLWN2(N6,N5,N4)=VOLWM(M,N6,N5,N4)*SN2OL(N6,N5,N4)
      VOLWN3(N6,N5,N4)=VOLWMA(N6,N5,N4)*SNH3L(N6,N5,N4)
      VOLWNB(N6,N5,N4)=VOLWMB(N6,N5,N4)*SNH3L(N6,N5,N4)
      VOLWHG(N6,N5,N4)=VOLWM(M,N6,N5,N4)*SH2GL(N6,N5,N4)
      ENDIF
      FLQM(N,N6,N5,N4)=(FLWM(M,N,N6,N5,N4)+FLWHM(M,N,N6,N5,N4))*XNPT
!
!     SOLUTE TRANSPORT IN MICROPORES
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     FLWM=water flux through soil micropore from watsub.f
!     VOLWM=micropore water-filled porosity from watsub.f
!     RFL*S=solute diffusive flux through micropore
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *S2,*B2=micropore solute content in non-band,band
!
      IF(FLWM(M,N,N6,N5,N4).GT.0.0)THEN
      IF(VOLWM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMAX1(0.0,AMIN1(VFLWX,FLWM(M,N,N6,N5,N4) &
      /VOLWM(M,N3,N2,N1)))
      ELSE
      VFLW=VFLWX
      ENDIF
      DO 9820 K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQC2(K,N3,N2,N1))
      RFLON(K)=VFLW*AMAX1(0.0,OQN2(K,N3,N2,N1))
      RFLOP(K)=VFLW*AMAX1(0.0,OQP2(K,N3,N2,N1))
      RFLOA(K)=VFLW*AMAX1(0.0,OQA2(K,N3,N2,N1))
9820  CONTINUE
      RFLCOS=VFLW*AMAX1(0.0,CO2S2(N3,N2,N1))
      RFLCHS=VFLW*AMAX1(0.0,CH4S2(N3,N2,N1))
      RFLOXS=VFLW*AMAX1(0.0,OXYS2(N3,N2,N1))
      RFLNGS=VFLW*AMAX1(0.0,Z2GS2(N3,N2,N1))
      RFLN2S=VFLW*AMAX1(0.0,Z2OS2(N3,N2,N1))
      RFLHGS=VFLW*AMAX1(0.0,H2GS2(N3,N2,N1))
      RFLNH4=VFLW*AMAX1(0.0,ZNH4S2(N3,N2,N1))
      RFLNH3=VFLW*AMAX1(0.0,ZNH3S2(N3,N2,N1))
      RFLNO3=VFLW*AMAX1(0.0,ZNO3S2(N3,N2,N1))
      RFLNO2=VFLW*AMAX1(0.0,ZNO2S2(N3,N2,N1))
      RFLP14=VFLW*AMAX1(0.0,H1PO42(N3,N2,N1))
      RFLPO4=VFLW*AMAX1(0.0,H2PO42(N3,N2,N1))
      RFLN4B=VFLW*AMAX1(0.0,ZNH4B2(N3,N2,N1))
      RFLN3B=VFLW*AMAX1(0.0,ZNH3B2(N3,N2,N1))
      RFLNOB=VFLW*AMAX1(0.0,ZNO3B2(N3,N2,N1))
      RFLN2B=VFLW*AMAX1(0.0,ZNO2B2(N3,N2,N1))
      RFLP1B=VFLW*AMAX1(0.0,H1POB2(N3,N2,N1))
      RFLPOB=VFLW*AMAX1(0.0,H2POB2(N3,N2,N1))
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS TO CURRENT FROM
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
      ELSE
      IF(VOLWM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AMIN1(0.0,AMAX1(-VFLWX,FLWM(M,N,N6,N5,N4) &
      /VOLWM(M,N6,N5,N4)))
      ELSE
      VFLW=-VFLWX
      ENDIF
      DO 9815 K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQC2(K,N6,N5,N4))
      RFLON(K)=VFLW*AMAX1(0.0,OQN2(K,N6,N5,N4))
      RFLOP(K)=VFLW*AMAX1(0.0,OQP2(K,N6,N5,N4))
      RFLOA(K)=VFLW*AMAX1(0.0,OQA2(K,N6,N5,N4))
9815  CONTINUE
      RFLCOS=VFLW*AMAX1(0.0,CO2S2(N6,N5,N4))
      RFLCHS=VFLW*AMAX1(0.0,CH4S2(N6,N5,N4))
      RFLOXS=VFLW*AMAX1(0.0,OXYS2(N6,N5,N4))
      RFLNGS=VFLW*AMAX1(0.0,Z2GS2(N6,N5,N4))
      RFLN2S=VFLW*AMAX1(0.0,Z2OS2(N6,N5,N4))
      RFLHGS=VFLW*AMAX1(0.0,H2GS2(N6,N5,N4))
      RFLNH4=VFLW*AMAX1(0.0,ZNH4S2(N6,N5,N4))
      RFLNH3=VFLW*AMAX1(0.0,ZNH3S2(N6,N5,N4))
      RFLNO3=VFLW*AMAX1(0.0,ZNO3S2(N6,N5,N4))
      RFLNO2=VFLW*AMAX1(0.0,ZNO2S2(N6,N5,N4))
      RFLP14=VFLW*AMAX1(0.0,H1PO42(N6,N5,N4))
      RFLPO4=VFLW*AMAX1(0.0,H2PO42(N6,N5,N4))
      RFLN4B=VFLW*AMAX1(0.0,ZNH4B2(N6,N5,N4))
      RFLN3B=VFLW*AMAX1(0.0,ZNH3B2(N6,N5,N4))
      RFLNOB=VFLW*AMAX1(0.0,ZNO3B2(N6,N5,N4))
      RFLN2B=VFLW*AMAX1(0.0,ZNO2B2(N6,N5,N4))
      RFLP1B=VFLW*AMAX1(0.0,H1POB2(N6,N5,N4))
      RFLPOB=VFLW*AMAX1(0.0,H2POB2(N6,N5,N4))
      ENDIF
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MICROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
      IF(THETW1(N3,N2,N1).GT.THETY(N3,N2,N1) &
      .AND.THETW1(N6,N5,N4).GT.THETY(N6,N5,N4) &
      .AND.VOLWM(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
      .AND.VOLWM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
!
!     VOLWM=micropore water-filled porosity from watsub.f
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
      DO 9810 K=0,jcplx1
      COQC1(K)=AMAX1(0.0,OQC2(K,N3,N2,N1)/VOLWM(M,N3,N2,N1))
      COQN1(K)=AMAX1(0.0,OQN2(K,N3,N2,N1)/VOLWM(M,N3,N2,N1))
      COQP1(K)=AMAX1(0.0,OQP2(K,N3,N2,N1)/VOLWM(M,N3,N2,N1))
      COQA1(K)=AMAX1(0.0,OQA2(K,N3,N2,N1)/VOLWM(M,N3,N2,N1))
      COQC2(K)=AMAX1(0.0,OQC2(K,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      COQN2(K)=AMAX1(0.0,OQN2(K,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      COQP2(K)=AMAX1(0.0,OQP2(K,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      COQA2(K)=AMAX1(0.0,OQA2(K,N6,N5,N4)/VOLWM(M,N6,N5,N4))
9810  CONTINUE
      CCO2S1=AMAX1(0.0,CO2S2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
      CCH4S1=AMAX1(0.0,CH4S2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
      COXYS1=AMAX1(0.0,OXYS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
      CZ2GS1=AMAX1(0.0,Z2GS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
      CZ2OS1=AMAX1(0.0,Z2OS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
      CH2GS1=AMAX1(0.0,H2GS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
      IF(VOLW4A.GT.ZEROS2(N2,N1))THEN
      CNH4S1=AMAX1(0.0,ZNH4S2(N3,N2,N1)/VOLW4A)
      CNH3S1=AMAX1(0.0,ZNH3S2(N3,N2,N1)/VOLW4A)
      ELSE
      CNH4S1=0.0
      CNH3S1=0.0
      ENDIF
      IF(VOLW3A.GT.ZEROS2(N2,N1))THEN
      CNO3S1=AMAX1(0.0,ZNO3S2(N3,N2,N1)/VOLW3A)
      CNO2S1=AMAX1(0.0,ZNO2S2(N3,N2,N1)/VOLW3A)
      ELSE
      CNO3S1=0.0
      CNO2S1=0.0
      ENDIF
      IF(VOLW2A.GT.ZEROS2(N2,N1))THEN
      CP14S1=AMAX1(0.0,H1PO42(N3,N2,N1)/VOLW2A)
      CPO4S1=AMAX1(0.0,H2PO42(N3,N2,N1)/VOLW2A)
      ELSE
      CP14S1=0.0
      CPO4S1=0.0
      ENDIF
      IF(VOLW4B.GT.ZEROS2(N2,N1))THEN
      CNH4B1=AMAX1(0.0,ZNH4B2(N3,N2,N1)/VOLW4B)
      CNH3B1=AMAX1(0.0,ZNH3B2(N3,N2,N1)/VOLW4B)
      ELSE
      CNH4B1=0.0
      CNH3B1=0.0
      ENDIF
      IF(VOLW3B.GT.ZEROS2(N2,N1))THEN
      CNO3B1=AMAX1(0.0,ZNO3B2(N3,N2,N1)/VOLW3B)
      CNO2B1=AMAX1(0.0,ZNO2B2(N3,N2,N1)/VOLW3B)
      ELSE
      CNO3B1=CNO3S1
      CNO2B1=CNO2S1
      ENDIF
      IF(VOLW2B.GT.ZEROS2(N2,N1))THEN
      CP14B1=AMAX1(0.0,H1POB2(N3,N2,N1)/VOLW2B)
      CPO4B1=AMAX1(0.0,H2POB2(N3,N2,N1)/VOLW2B)
      ELSE
      CP14B1=CP14S1
      CPO4B1=CPO4S1
      ENDIF
      CCO2S2=AMAX1(0.0,CO2S2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CCH4S2=AMAX1(0.0,CH4S2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      COXYS2=AMAX1(0.0,OXYS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CZ2GS2=AMAX1(0.0,Z2GS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CZ2OS2=AMAX1(0.0,Z2OS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CH2GS2=AMAX1(0.0,H2GS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      IF(VOLWMA(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      CNH3S2=AMAX1(0.0,ZNH3S2(N6,N5,N4)/VOLWMA(N6,N5,N4))
      CNH4S2=AMAX1(0.0,ZNH4S2(N6,N5,N4)/VOLWMA(N6,N5,N4))
      ELSE
      CNH3S2=0.0
      CNH4S2=0.0
      ENDIF
      IF(VOLWOA.GT.ZEROS2(N5,N4))THEN
      CNO3S2=AMAX1(0.0,ZNO3S2(N6,N5,N4)/VOLWOA)
      CNO2S2=AMAX1(0.0,ZNO2S2(N6,N5,N4)/VOLWOA)
      ELSE
      CNO3S2=0.0
      CNO2S2=0.0
      ENDIF
      IF(VOLWPA.GT.ZEROS2(N5,N4))THEN
      CP14S2=AMAX1(0.0,H1PO42(N6,N5,N4)/VOLWPA)
      CPO4S2=AMAX1(0.0,H2PO42(N6,N5,N4)/VOLWPA)
      ELSE
      CP14S2=0.0
      CPO4S2=0.0
      ENDIF
      IF(VOLWMB(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      CNH3B2=AMAX1(0.0,ZNH3B2(N6,N5,N4)/VOLWMB(N6,N5,N4))
      CNH4B2=AMAX1(0.0,ZNH4B2(N6,N5,N4)/VOLWMB(N6,N5,N4))
      ELSE
      CNH3B2=CNH3S2
      CNH4B2=CNH4S2
      ENDIF
      IF(VOLWOB.GT.ZEROS2(N5,N4))THEN
      CNO3B2=AMAX1(0.0,ZNO3B2(N6,N5,N4)/VOLWOB)
      CNO2B2=AMAX1(0.0,ZNO2B2(N6,N5,N4)/VOLWOB)
      ELSE
      CNO3B2=CNO3S2
      CNO2B2=CNO2S2
      ENDIF
      IF(VOLWPB.GT.ZEROS2(N5,N4))THEN
      CP14B2=AMAX1(0.0,H1POB2(N6,N5,N4)/VOLWPB)
      CPO4B2=AMAX1(0.0,H2POB2(N6,N5,N4)/VOLWPB)
      ELSE
      CP14B2=CP14S2
      CPO4B2=CPO4S2
      ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!
!     DLYR=soil layer thickness
!     TORT=micropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     FLWM=water flux through soil micropore from watsub.f
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
      TORTL=(TORT(M,N3,N2,N1)*DLYR1+TORT(M,N6,N5,N4)*DLYR2) &
      /(DLYR1+DLYR2)
      DISPN=DISP(N,N6,N5,N4) &
      *AMIN1(VFLWX,ABS(FLWM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
      DIFOC=(OCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFON=(ONSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFOP=(OPSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFOA=(OASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFNH=(ZNSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFNO=(ZOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFCS=(CLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFCQ=(CQSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFOS=(OLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFNG=(ZLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFN2=(ZVSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFHG=(HLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
      DO 9805 K=0,jcplx1
      DFVOC(K)=DIFOC*(COQC1(K)-COQC2(K))
      DFVON(K)=DIFON*(COQN1(K)-COQN2(K))
      DFVOP(K)=DIFOP*(COQP1(K)-COQP2(K))
      DFVOA(K)=DIFOA*(COQA1(K)-COQA2(K))
9805  CONTINUE
      DFVCOS=DIFCS*(CCO2S1-CCO2S2)
      DFVCHS=DIFCQ*(CCH4S1-CCH4S2)
      DFVOXS=DIFOS*(COXYS1-COXYS2)
      DFVNGS=DIFNG*(CZ2GS1-CZ2GS2)
      DFVN2S=DIFN2*(CZ2OS1-CZ2OS2)
      DFVHGS=DIFHG*(CH2GS1-CH2GS2)
      DFVNH4=DIFNH*(CNH4S1-CNH4S2)*AMIN1(VLNH4(N3,N2,N1) &
      ,VLNH4(N6,N5,N4))
      DFVNH3=DIFNH*(CNH3S1-CNH3S2)*AMIN1(VLNH4(N3,N2,N1) &
      ,VLNH4(N6,N5,N4))
      DFVNO3=DIFNO*(CNO3S1-CNO3S2)*AMIN1(VLNO3(N3,N2,N1) &
      ,VLNO3(N6,N5,N4))
      DFVNO2=DIFNO*(CNO2S1-CNO2S2)*AMIN1(VLNO3(N3,N2,N1) &
      ,VLNO3(N6,N5,N4))
      DFVP14=DIFPO*(CP14S1-CP14S2)*AMIN1(VLPO4(N3,N2,N1) &
      ,VLPO4(N6,N5,N4))
      DFVPO4=DIFPO*(CPO4S1-CPO4S2)*AMIN1(VLPO4(N3,N2,N1) &
      ,VLPO4(N6,N5,N4))
      DFVN4B=DIFNH*(CNH4B1-CNH4B2)*AMIN1(VLNHB(N3,N2,N1) &
      ,VLNHB(N6,N5,N4))
      DFVN3B=DIFNH*(CNH3B1-CNH3B2)*AMIN1(VLNHB(N3,N2,N1) &
      ,VLNHB(N6,N5,N4))
      DFVNOB=DIFNO*(CNO3B1-CNO3B2)*AMIN1(VLNOB(N3,N2,N1) &
      ,VLNOB(N6,N5,N4))
      DFVN2B=DIFNO*(CNO2B1-CNO2B2)*AMIN1(VLNOB(N3,N2,N1) &
      ,VLNOB(N6,N5,N4))
      DFVP1B=DIFPO*(CP14B1-CP14B2)*AMIN1(VLPOB(N3,N2,N1) &
      ,VLPOB(N6,N5,N4))
      DFVPOB=DIFPO*(CPO4B1-CPO4B2)*AMIN1(VLPOB(N3,N2,N1) &
      ,VLPOB(N6,N5,N4))
      ELSE
      DO 9905 K=0,jcplx1
      DFVOC(K)=0.0
      DFVON(K)=0.0
      DFVOP(K)=0.0
      DFVOA(K)=0.0
9905  CONTINUE
      DFVCOS=0.0
      DFVCHS=0.0
      DFVOXS=0.0
      DFVNGS=0.0
      DFVN2S=0.0
      DFVHGS=0.0
      DFVNH4=0.0
      DFVNH3=0.0
      DFVNO3=0.0
      DFVNO2=0.0
      DFVP14=0.0
      DFVPO4=0.0
      DFVN4B=0.0
      DFVN3B=0.0
      DFVNOB=0.0
      DFVN2B=0.0
      DFVP1B=0.0
      DFVPOB=0.0
      ENDIF
!
!     SOLUTE TRANSPORT IN MACROPORES
!
!     FLWHM=water flux through soil macropore from watsub.f
!
      IF(FLWHM(M,N,N6,N5,N4).GT.0.0)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN CURRENT
!     GRID CELL
!
!     VOLWHM=macropore water-filled porosity from watsub.f
!     VOLWAH=macropore porosity
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
      IF(VOLWHM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMAX1(0.0,AMIN1(VFLWX,FLWHM(M,N,N6,N5,N4) &
      /VOLWHM(M,N3,N2,N1)))
      ELSE
      VFLW=VFLWX
      ENDIF
!
!     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE
!
      IF(N.EQ.3.AND.VOLAH(N6,N5,N4).GT.VOLWHM(M,N6,N5,N4))THEN
      DO 9800 K=0,jcplx1
      RFHOC(K)=VFLW*AMAX1(0.0,(OQCH2(K,N3,N2,N1) &
      -AMIN1(0.0,ROCFXS(K,NU(N2,N1),N2,N1))))
      RFHON(K)=VFLW*AMAX1(0.0,(OQNH2(K,N3,N2,N1) &
      -AMIN1(0.0,RONFXS(K,NU(N2,N1),N2,N1))))
      RFHOP(K)=VFLW*AMAX1(0.0,(OQPH2(K,N3,N2,N1) &
      -AMIN1(0.0,ROPFXS(K,NU(N2,N1),N2,N1))))
      RFHOA(K)=VFLW*AMAX1(0.0,(OQAH2(K,N3,N2,N1) &
      -AMIN1(0.0,ROAFXS(K,NU(N2,N1),N2,N1))))
9800  CONTINUE
      RFHCOS=VFLW*AMAX1(0.0,(CO2SH2(N3,N2,N1) &
      -AMIN1(0.0,RCOFXS(NU(N2,N1),N2,N1))))
      RFHCHS=VFLW*AMAX1(0.0,(CH4SH2(N3,N2,N1) &
      -AMIN1(0.0,RCHFXS(NU(N2,N1),N2,N1))))
      RFHOXS=VFLW*AMAX1(0.0,(OXYSH2(N3,N2,N1) &
      -AMIN1(0.0,ROXFXS(NU(N2,N1),N2,N1))))
      RFHNGS=VFLW*AMAX1(0.0,(Z2GSH2(N3,N2,N1) &
      -AMIN1(0.0,RNGFXS(NU(N2,N1),N2,N1))))
      RFHN2S=VFLW*AMAX1(0.0,(Z2OSH2(N3,N2,N1) &
      -AMIN1(0.0,RN2FXS(NU(N2,N1),N2,N1))))
      RFHHGS=VFLW*AMAX1(0.0,(H2GSH2(N3,N2,N1) &
      -AMIN1(0.0,RHGFXS(NU(N2,N1),N2,N1))))
      RFHNH4=VFLW*AMAX1(0.0,(ZNH4H2(N3,N2,N1) &
      -AMIN1(0.0,RN4FXW(NU(N2,N1),N2,N1)*VLNH4(N3,N2,N1)))) &
      *VLNH4(N6,N5,N4)
      RFHNH3=VFLW*AMAX1(0.0,(ZNH3H2(N3,N2,N1) &
      -AMIN1(0.0,RN3FXW(NU(N2,N1),N2,N1)*VLNH4(N3,N2,N1)))) &
      *VLNH4(N6,N5,N4)
      RFHNO3=VFLW*AMAX1(0.0,(ZNO3H2(N3,N2,N1) &
      -AMIN1(0.0,RNOFXW(NU(N2,N1),N2,N1)*VLNO3(N3,N2,N1)))) &
      *VLNO3(N6,N5,N4)
      RFHNO2=VFLW*AMAX1(0.0,(ZNO2H2(N3,N2,N1) &
      -AMIN1(0.0,RNXFXS(NU(N2,N1),N2,N1)*VLNO3(N3,N2,N1)))) &
      *VLNO3(N6,N5,N4)
      RFHP14=VFLW*AMAX1(0.0,(H1P4H2(N3,N2,N1) &
      -AMIN1(0.0,RH1PXS(NU(N2,N1),N2,N1)*VLPO4(N3,N2,N1)))) &
      *VLPO4(N6,N5,N4)
      RFHPO4=VFLW*AMAX1(0.0,(H2P4H2(N3,N2,N1) &
      -AMIN1(0.0,RH2PXS(NU(N2,N1),N2,N1)*VLPO4(N3,N2,N1)))) &
      *VLPO4(N6,N5,N4)
      RFHN4B=VFLW*AMAX1(0.0,(ZN4BH2(N3,N2,N1) &
      -AMIN1(0.0,RN4FXB(NU(N2,N1),N2,N1)*VLNHB(N3,N2,N1)))) &
      *VLNHB(N6,N5,N4)
      RFHN3B=VFLW*AMAX1(0.0,(ZN3BH2(N3,N2,N1) &
      -AMIN1(0.0,RN3FXB(NU(N2,N1),N2,N1)*VLNHB(N3,N2,N1)))) &
      *VLNHB(N6,N5,N4)
      RFHNOB=VFLW*AMAX1(0.0,(ZNOBH2(N3,N2,N1) &
      -AMIN1(0.0,RNOFXB(NU(N2,N1),N2,N1)*VLNOB(N3,N2,N1)))) &
      *VLNOB(N6,N5,N4)
      RFHN2B=VFLW*AMAX1(0.0,(ZN2BH2(N3,N2,N1) &
      -AMIN1(0.0,RNXFXB(NU(N2,N1),N2,N1)*VLNOB(N3,N2,N1)))) &
      *VLNOB(N6,N5,N4)
      RFHP1B=VFLW*AMAX1(0.0,(H1PBH2(N3,N2,N1) &
      -AMIN1(0.0,RH1BXB(NU(N2,N1),N2,N1)*VLPOB(N3,N2,N1)))) &
      *VLPOB(N6,N5,N4)
      RFHPOB=VFLW*AMAX1(0.0,(H2PBH2(N3,N2,N1) &
      -AMIN1(0.0,RH2BXB(NU(N2,N1),N2,N1)*VLPOB(N3,N2,N1)))) &
      *VLPOB(N6,N5,N4)
!
!     OTHERWISE
!
      ELSE
      DO 9850 K=0,jcplx1
      RFHOC(K)=VFLW*AMAX1(0.0,OQCH2(K,N3,N2,N1))
      RFHON(K)=VFLW*AMAX1(0.0,OQNH2(K,N3,N2,N1))
      RFHOP(K)=VFLW*AMAX1(0.0,OQPH2(K,N3,N2,N1))
      RFHOA(K)=VFLW*AMAX1(0.0,OQAH2(K,N3,N2,N1))
9850  CONTINUE
      RFHCOS=VFLW*AMAX1(0.0,CO2SH2(N3,N2,N1))
      RFHCHS=VFLW*AMAX1(0.0,CH4SH2(N3,N2,N1))
      RFHOXS=VFLW*AMAX1(0.0,OXYSH2(N3,N2,N1))
      RFHNGS=VFLW*AMAX1(0.0,Z2GSH2(N3,N2,N1))
      RFHN2S=VFLW*AMAX1(0.0,Z2OSH2(N3,N2,N1))
      RFHHGS=VFLW*AMAX1(0.0,H2GSH2(N3,N2,N1))
      RFHNH4=VFLW*AMAX1(0.0,ZNH4H2(N3,N2,N1))*VLNH4(N6,N5,N4)
      RFHNH3=VFLW*AMAX1(0.0,ZNH3H2(N3,N2,N1))*VLNH4(N6,N5,N4)
      RFHNO3=VFLW*AMAX1(0.0,ZNO3H2(N3,N2,N1))*VLNO3(N6,N5,N4)
      RFHNO2=VFLW*AMAX1(0.0,ZNO2H2(N3,N2,N1))*VLNO3(N6,N5,N4)
      RFHP14=VFLW*AMAX1(0.0,H1P4H2(N3,N2,N1))*VLPO4(N6,N5,N4)
      RFHPO4=VFLW*AMAX1(0.0,H2P4H2(N3,N2,N1))*VLPO4(N6,N5,N4)
      RFHN4B=VFLW*AMAX1(0.0,ZN4BH2(N3,N2,N1))*VLNHB(N6,N5,N4)
      RFHN3B=VFLW*AMAX1(0.0,ZN3BH2(N3,N2,N1))*VLNHB(N6,N5,N4)
      RFHNOB=VFLW*AMAX1(0.0,ZNOBH2(N3,N2,N1))*VLNOB(N6,N5,N4)
      RFHN2B=VFLW*AMAX1(0.0,ZN2BH2(N3,N2,N1))*VLNOB(N6,N5,N4)
      RFHP1B=VFLW*AMAX1(0.0,H1PBH2(N3,N2,N1))*VLPOB(N6,N5,N4)
      RFHPOB=VFLW*AMAX1(0.0,H2PBH2(N3,N2,N1))*VLPOB(N6,N5,N4)
      ENDIF
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN ADJACENT
!     GRID CELL
!
      ELSEIF(FLWHM(M,N,N6,N5,N4).LT.0.0)THEN
      IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AMIN1(0.0,AMAX1(-VFLWX,FLWHM(M,N,N6,N5,N4) &
      /VOLWHM(M,N6,N5,N4)))
      ELSE
      VFLW=-VFLWX
      ENDIF
      DO 9665 K=0,jcplx1
      RFHOC(K)=VFLW*AMAX1(0.0,OQCH2(K,N6,N5,N4))
      RFHON(K)=VFLW*AMAX1(0.0,OQNH2(K,N6,N5,N4))
      RFHOP(K)=VFLW*AMAX1(0.0,OQPH2(K,N6,N5,N4))
      RFHOA(K)=VFLW*AMAX1(0.0,OQAH2(K,N6,N5,N4))
9665  CONTINUE
      RFHCOS=VFLW*AMAX1(0.0,CO2SH2(N6,N5,N4))
      RFHCHS=VFLW*AMAX1(0.0,CH4SH2(N6,N5,N4))
      RFHOXS=VFLW*AMAX1(0.0,OXYSH2(N6,N5,N4))
      RFHNGS=VFLW*AMAX1(0.0,Z2GSH2(N6,N5,N4))
      RFHN2S=VFLW*AMAX1(0.0,Z2OSH2(N6,N5,N4))
      RFHHGS=VFLW*AMAX1(0.0,H2GSH2(N6,N5,N4))
      RFHNH4=VFLW*AMAX1(0.0,ZNH4H2(N6,N5,N4))*VLNH4(N6,N5,N4)
      RFHNH3=VFLW*AMAX1(0.0,ZNH3H2(N6,N5,N4))*VLNH4(N6,N5,N4)
      RFHNO3=VFLW*AMAX1(0.0,ZNO3H2(N6,N5,N4))*VLNO3(N6,N5,N4)
      RFHNO2=VFLW*AMAX1(0.0,ZNO2H2(N6,N5,N4))*VLNO3(N6,N5,N4)
      RFHP14=VFLW*AMAX1(0.0,H1P4H2(N6,N5,N4))*VLPO4(N6,N5,N4)
      RFHPO4=VFLW*AMAX1(0.0,H2P4H2(N6,N5,N4))*VLPO4(N6,N5,N4)
      RFHN4B=VFLW*AMAX1(0.0,ZN4BH2(N6,N5,N4))*VLNHB(N6,N5,N4)
      RFHN3B=VFLW*AMAX1(0.0,ZN3BH2(N6,N5,N4))*VLNHB(N6,N5,N4)
      RFHNOB=VFLW*AMAX1(0.0,ZNOBH2(N6,N5,N4))*VLNOB(N6,N5,N4)
      RFHN2B=VFLW*AMAX1(0.0,ZN2BH2(N6,N5,N4))*VLNOB(N6,N5,N4)
      RFHP1B=VFLW*AMAX1(0.0,H1PBH2(N6,N5,N4))*VLPOB(N6,N5,N4)
      RFHPOB=VFLW*AMAX1(0.0,H2PBH2(N6,N5,N4))*VLPOB(N6,N5,N4)
      ELSE
!
!     NO MACROPORE FLUX
!
      DO 9795 K=0,jcplx1
      RFHOC(K)=0.0
      RFHON(K)=0.0
      RFHOP(K)=0.0
      RFHOA(K)=0.0
9795  CONTINUE
      RFHCOS=0.0
      RFHCHS=0.0
      RFHOXS=0.0
      RFHNGS=0.0
      RFHN2S=0.0
      RFHHGS=0.0
      RFHNH4=0.0
      RFHNH3=0.0
      RFHNO3=0.0
      RFHNO2=0.0
      RFHP14=0.0
      RFHPO4=0.0
      RFHN4B=0.0
      RFHN3B=0.0
      RFHNOB=0.0
      RFHN2B=0.0
      RFHP1B=0.0
      RFHPOB=0.0
      ENDIF
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MACROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VOLWHM=macropore water-filled porosity from watsub.f
!     THETY=hygroscopic water content
!     VOLAH=total macropore volume
!
      IF(VOLWHM(M,N3,N2,N1).GT.THETY(N3,N2,N1)*VOLAH(N3,N2,N1) &
      .AND.VOLWHM(M,N6,N5,N4).GT.THETY(N6,N5,N4)*VOLAH(N6,N5,N4))THEN
!
!     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
!
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     *H2=macropore solute content
!     VOLWHM=macropore water content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
      DO 9790 K=0,jcplx1
      COQCH1(K)=AMAX1(0.0,OQCH2(K,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      COQNH1(K)=AMAX1(0.0,OQNH2(K,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      COQPH1(K)=AMAX1(0.0,OQPH2(K,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      COQAH1(K)=AMAX1(0.0,OQAH2(K,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      COQCH2(K)=AMAX1(0.0,OQCH2(K,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      COQNH2(K)=AMAX1(0.0,OQNH2(K,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      COQPH2(K)=AMAX1(0.0,OQPH2(K,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      COQAH2(K)=AMAX1(0.0,OQAH2(K,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
9790  CONTINUE
      CCO2SH1=AMAX1(0.0,CO2SH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      CCH4SH1=AMAX1(0.0,CH4SH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      COXYSH1=AMAX1(0.0,OXYSH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      CZ2GSH1=AMAX1(0.0,Z2GSH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      CZ2OSH1=AMAX1(0.0,Z2OSH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      CH2GSH1=AMAX1(0.0,H2GSH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      IF(VOLH4A.GT.ZEROS2(N2,N1))THEN
      CNH4SH1=AMAX1(0.0,ZNH4H2(N3,N2,N1)/VOLH4A)
      CNH3SH1=AMAX1(0.0,ZNH3H2(N3,N2,N1)/VOLH4A)
      ELSE
      CNH4SH1=0.0
      CNH3SH1=0.0
      ENDIF
      IF(VOLH3A.GT.ZEROS2(N2,N1))THEN
      CNO3SH1=AMAX1(0.0,ZNO3H2(N3,N2,N1)/VOLH3A)
      CNO2SH1=AMAX1(0.0,ZNO2H2(N3,N2,N1)/VOLH3A)
      ELSE
      CNO3SH1=0.0
      CNO2SH1=0.0
      ENDIF
      IF(VOLH2A.GT.ZEROS2(N2,N1))THEN
      CP14SH1=AMAX1(0.0,H1P4H2(N3,N2,N1)/VOLH2A)
      CPO4SH1=AMAX1(0.0,H2P4H2(N3,N2,N1)/VOLH2A)
      ELSE
      CP14SH1=0.0
      CPO4SH1=0.0
      ENDIF
      IF(VOLH4B.GT.ZEROS2(N2,N1))THEN
      CNH4BH1=AMAX1(0.0,ZN4BH2(N3,N2,N1)/VOLH4B)
      CNH3BH1=AMAX1(0.0,ZN3BH2(N3,N2,N1)/VOLH4B)
      ELSE
      CNH4BH1=CNH4SH1
      CNH3BH1=CNH3SH1
      ENDIF
      IF(VOLH3B.GT.ZEROS2(N2,N1))THEN
      CNO3BH1=AMAX1(0.0,ZNOBH2(N3,N2,N1)/VOLH3B)
      CNO2BH1=AMAX1(0.0,ZN2BH2(N3,N2,N1)/VOLH3B)
      ELSE
      CNO3BH1=CNO3SH1
      CNO2BH1=CNO2SH1
      ENDIF
      IF(VOLH2B.GT.ZEROS2(N2,N1))THEN
      CP14BH1=AMAX1(0.0,H1PBH2(N3,N2,N1)/VOLH2B)
      CPO4BH1=AMAX1(0.0,H2PBH2(N3,N2,N1)/VOLH2B)
      ELSE
      CP14BH1=CP14SH1
      CPO4BH1=CPO4SH1
      ENDIF
      CCO2SH2=AMAX1(0.0,CO2SH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      CCH4SH2=AMAX1(0.0,CH4SH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      COXYSH2=AMAX1(0.0,OXYSH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      CZ2GSH2=AMAX1(0.0,Z2GSH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      CZ2OSH2=AMAX1(0.0,Z2OSH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      CH2GSH2=AMAX1(0.0,H2GSH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
      VOLHMA=VOLWHM(M,N6,N5,N4)*VLNH4(N6,N5,N4)
      IF(VOLHMA.GT.ZEROS2(N5,N4))THEN
      CNH4SH2=AMAX1(0.0,ZNH4H2(N6,N5,N4)/VOLHMA)
      CNH3SH2=AMAX1(0.0,ZNH3H2(N6,N5,N4)/VOLHMA)
      ELSE
      CNH4SH2=0.0
      CNH3SH2=0.0
      ENDIF
      VOLHOA=VOLWHM(M,N6,N5,N4)*VLNO3(N6,N5,N4)
      IF(VOLHOA.GT.ZEROS2(N5,N4))THEN
      CNO3SH2=AMAX1(0.0,ZNO3H2(N6,N5,N4)/VOLHOA)
      CNO2SH2=AMAX1(0.0,ZNO2H2(N6,N5,N4)/VOLHOA)
      ELSE
      CNO3SH2=0.0
      CNO2SH2=0.0
      ENDIF
      VOLHPA=VOLWHM(M,N6,N5,N4)*VLPO4(N6,N5,N4)
      IF(VOLHPA.GT.ZEROS2(N5,N4))THEN
      CP14SH2=AMAX1(0.0,H1P4H2(N6,N5,N4)/VOLHPA)
      CPO4SH2=AMAX1(0.0,H2P4H2(N6,N5,N4)/VOLHPA)
      ELSE
      CP14SH2=0.0
      CPO4SH2=0.0
      ENDIF
      VOLHMB=VOLWHM(M,N6,N5,N4)*VLNHB(N6,N5,N4)
      IF(VOLHMB.GT.ZEROS2(N5,N4))THEN
      CNH4BH2=AMAX1(0.0,ZN4BH2(N6,N5,N4)/VOLHMB)
      CNH3BH2=AMAX1(0.0,ZN3BH2(N6,N5,N4)/VOLHMB)
      ELSE
      CNH4BH2=CNH4SH2
      CNH3BH2=CNH3SH2
      ENDIF
      VOLHOB=VOLWHM(M,N6,N5,N4)*VLNOB(N6,N5,N4)
      IF(VOLHOB.GT.ZEROS2(N5,N4))THEN
      CNO3BH2=AMAX1(0.0,ZNOBH2(N6,N5,N4)/VOLHOB)
      CNO2BH2=AMAX1(0.0,ZN2BH2(N6,N5,N4)/VOLHOB)
      ELSE
      CNO3BH2=CNO3SH2
      CNO2BH2=CNO2SH2
      ENDIF
      VOLHPB=VOLWHM(M,N6,N5,N4)*VLPOB(N6,N5,N4)
      IF(VOLHPB.GT.ZEROS2(N5,N4))THEN
      CP14BH2=AMAX1(0.0,H1PBH2(N6,N5,N4)/VOLHPB)
      CPO4BH2=AMAX1(0.0,H2PBH2(N6,N5,N4)/VOLHPB)
      ELSE
      CP14BH2=CP14SH2
      CPO4BH2=CPO4SH2
      ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
!
!     DLYR=soil layer thickness
!     TORTH=macropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     FLWHM=water flux through soil macropore from watsub.f
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
      TORTL=(TORTH(M,N3,N2,N1)*DLYR1+TORTH(M,N6,N5,N4)*DLYR2) &
      /(DLYR1+DLYR2)
      DISPN=DISP(N,N6,N5,N4) &
      *AMIN1(VFLWX,ABS(FLWHM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
      DIFOC=(OCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFON=(ONSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFOP=(OPSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFOA=(OASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFNH=(ZNSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFNO=(ZOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFCS=(CLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFCQ=(CQSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFOS=(OLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFNG=(ZLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFN2=(ZVSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
      DIFHG=(HLSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MACROPORES
!
      DO 9785 K=0,jcplx1
      DFHOC(K)=DIFOC*(COQCH1(K)-COQCH2(K))
      DFHON(K)=DIFON*(COQNH1(K)-COQNH2(K))
      DFHOP(K)=DIFOP*(COQPH1(K)-COQPH2(K))
      DFHOA(K)=DIFOA*(COQAH1(K)-COQAH2(K))
!     WRITE(*,2121)'DFHOC',I,J,M,N4,N5,N6,K,DFHOC(K),OQCH2(K,N3,N2,N1)
!    2,OQCH2(K,N6,N5,N4),DIFOC,COQCH1(K),COQCH2(K)
!2121  FORMAT(A8,7I4,20E12.4)
9785  CONTINUE
      DFHCOS=DIFCS*(CCO2SH1-CCO2SH2)
      DFHCHS=DIFCQ*(CCH4SH1-CCH4SH2)
      DFHOXS=DIFOS*(COXYSH1-COXYSH2)
      DFHNGS=DIFNG*(CZ2GSH1-CZ2GSH2)
      DFHN2S=DIFN2*(CZ2OSH1-CZ2OSH2)
      DFHHGS=DIFNH*(CH2GSH1-CH2GSH2)
      DFHNH4=DIFNH*(CNH4SH1-CNH4SH2)*AMIN1(VLNH4(N3,N2,N1) &
      ,VLNH4(N6,N5,N4))
      DFHNH3=DIFNH*(CNH3SH1-CNH3SH2)*AMIN1(VLNH4(N3,N2,N1) &
      ,VLNH4(N6,N5,N4))
      DFHNO3=DIFNO*(CNO3SH1-CNO3SH2)*AMIN1(VLNO3(N3,N2,N1) &
      ,VLNO3(N6,N5,N4))
      DFHNO2=DIFNO*(CNO2SH1-CNO2SH2)*AMIN1(VLNO3(N3,N2,N1) &
      ,VLNO3(N6,N5,N4))
      DFHP14=DIFPO*(CP14SH1-CP14SH2)*AMIN1(VLPO4(N3,N2,N1) &
      ,VLPO4(N6,N5,N4))
      DFHPO4=DIFPO*(CPO4SH1-CPO4SH2)*AMIN1(VLPO4(N3,N2,N1) &
      ,VLPO4(N6,N5,N4))
      DFHN4B=DIFNH*(CNH4BH1-CNH4BH2)*AMIN1(VLNHB(N3,N2,N1) &
      ,VLNHB(N6,N5,N4))
      DFHN3B=DIFNH*(CNH3BH1-CNH3BH2)*AMIN1(VLNHB(N3,N2,N1) &
      ,VLNHB(N6,N5,N4))
      DFHNOB=DIFNO*(CNO3BH1-CNO3BH2)*AMIN1(VLNOB(N3,N2,N1) &
      ,VLNOB(N6,N5,N4))
      DFHN2B=DIFNO*(CNO2BH1-CNO2BH2)*AMIN1(VLNOB(N3,N2,N1) &
      ,VLNOB(N6,N5,N4))
      DFHP1B=DIFPO*(CP14BH1-CP14BH2)*AMIN1(VLPOB(N3,N2,N1) &
      ,VLPOB(N6,N5,N4))
      DFHPOB=DIFPO*(CPO4BH1-CPO4BH2)*AMIN1(VLPOB(N3,N2,N1) &
      ,VLPOB(N6,N5,N4))
      ELSE
      DO 9780 K=0,jcplx1
      DFHOC(K)=0.0
      DFHON(K)=0.0
      DFHOP(K)=0.0
      DFHOA(K)=0.0
9780  CONTINUE
      DFHCOS=0.0
      DFHCHS=0.0
      DFHOXS=0.0
      DFHNGS=0.0
      DFHN2S=0.0
      DFHHGS=0.0
      DFHNH4=0.0
      DFHNH3=0.0
      DFHNO3=0.0
      DFHNO2=0.0
      DFHP14=0.0
      DFHPO4=0.0
      DFHN4B=0.0
      DFHN3B=0.0
      DFHNOB=0.0
      DFHN2B=0.0
      DFHP1B=0.0
      DFHPOB=0.0
      ENDIF
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
      DO 9765 K=0,jcplx1
      ROCFLS(K,N,N6,N5,N4)=RFLOC(K)+DFVOC(K)
      RONFLS(K,N,N6,N5,N4)=RFLON(K)+DFVON(K)
      ROPFLS(K,N,N6,N5,N4)=RFLOP(K)+DFVOP(K)
      ROAFLS(K,N,N6,N5,N4)=RFLOA(K)+DFVOA(K)
      ROCFHS(K,N,N6,N5,N4)=RFHOC(K)+DFHOC(K)
      RONFHS(K,N,N6,N5,N4)=RFHON(K)+DFHON(K)
      ROPFHS(K,N,N6,N5,N4)=RFHOP(K)+DFHOP(K)
      ROAFHS(K,N,N6,N5,N4)=RFHOA(K)+DFHOA(K)
!     IF(N3.LE.3)THEN
!     WRITE(*,447)'ROCFLS',I,J,N4,N5,N6,M,MM,N,K
!    2,ROCFLS(K,N,N6,N5,N4),RFLOC(K),DFVOC(K),DIFOC
!    3,COQC1(K),COQC2(K),OQC2(K,N3,N2,N1),VOLWM(M,N3,N2,N1)
!    4,OQC2(K,N6,N5,N4),VOLWM(M,N6,N5,N4)
!    2,ROAFLS(K,N,N6,N5,N4),RFLOA(K),DFVOA(K),DIFOA
!    3,COQA1(K),COQA2(K),OQA2(K,N3,N2,N1),VOLWM(M,N3,N2,N1)
!    4,OQA2(K,N6,N5,N4),VOLWM(M,N6,N5,N4)
!447   FORMAT(A8,9I4,20E12.4)
!     ENDIF
9765  CONTINUE
      RCOFLS(N,N6,N5,N4)=RFLCOS+DFVCOS
      RCHFLS(N,N6,N5,N4)=RFLCHS+DFVCHS
      ROXFLS(N,N6,N5,N4)=RFLOXS+DFVOXS
      RNGFLS(N,N6,N5,N4)=RFLNGS+DFVNGS
      RN2FLS(N,N6,N5,N4)=RFLN2S+DFVN2S
      RHGFLS(N,N6,N5,N4)=RFLHGS+DFVHGS
      RN4FLW(N,N6,N5,N4)=RFLNH4+DFVNH4
      RN3FLW(N,N6,N5,N4)=RFLNH3+DFVNH3
      RNOFLW(N,N6,N5,N4)=RFLNO3+DFVNO3
      RNXFLS(N,N6,N5,N4)=RFLNO2+DFVNO2
      RH1PFS(N,N6,N5,N4)=RFLP14+DFVP14
      RH2PFS(N,N6,N5,N4)=RFLPO4+DFVPO4
      RN4FLB(N,N6,N5,N4)=RFLN4B+DFVN4B
      RN3FLB(N,N6,N5,N4)=RFLN3B+DFVN3B
      RNOFLB(N,N6,N5,N4)=RFLNOB+DFVNOB
      RNXFLB(N,N6,N5,N4)=RFLN2B+DFVN2B
      RH1BFB(N,N6,N5,N4)=RFLP1B+DFVP1B
      RH2BFB(N,N6,N5,N4)=RFLPOB+DFVPOB
      RCOFHS(N,N6,N5,N4)=RFHCOS+DFHCOS
      RCHFHS(N,N6,N5,N4)=RFHCHS+DFHCHS
      ROXFHS(N,N6,N5,N4)=RFHOXS+DFHOXS
      RNGFHS(N,N6,N5,N4)=RFHNGS+DFHNGS
      RN2FHS(N,N6,N5,N4)=RFHN2S+DFHN2S
      RHGFHS(N,N6,N5,N4)=RFHHGS+DFHHGS
      RN4FHW(N,N6,N5,N4)=RFHNH4+DFHNH4
      RN3FHW(N,N6,N5,N4)=RFHNH3+DFHNH3
      RNOFHW(N,N6,N5,N4)=RFHNO3+DFHNO3
      RNXFHS(N,N6,N5,N4)=RFHNO2+DFHNO2
      RH1PHS(N,N6,N5,N4)=RFHP14+DFHP14
      RH2PHS(N,N6,N5,N4)=RFHPO4+DFHPO4
      RN4FHB(N,N6,N5,N4)=RFHN4B+DFHN4B
      RN3FHB(N,N6,N5,N4)=RFHN3B+DFHN3B
      RNOFHB(N,N6,N5,N4)=RFHNOB+DFHNOB
      RNXFHB(N,N6,N5,N4)=RFHN2B+DFHN2B
      RH1BHB(N,N6,N5,N4)=RFHP1B+DFHP1B
      RH2BHB(N,N6,N5,N4)=RFHPOB+DFHPOB
!     IF(I.EQ.361.AND.NY.EQ.2)THEN
!     WRITE(*,443)'RN4FHB',I,J,M,MM,N1,N2,N3,N4,N5,N6,N
!    2,RN4FHB(N,N6,N5,N4),RFHN4B,DFHN4B,VFLW,ZN4BH2(N3,N2,N1)
!    2,RN4FXB(NU(N2,N1),N2,N1),VLNHB(N3,N2,N1),VLNHB(N6,N5,N4)
!    3,FLWHM(M,N,N6,N5,N4),VFLWX,VOLWHM(M,N3,N2,N1)
!    4,VOLWHM(M,N6,N5,N4)
!     WRITE(*,443)'RCOFLS',I,J,M,MM,N1,N2,N3,N4,N5,N6,N
!    2,RCOFLS(N,N6,N5,N4),RFLCOS,DFVCOS,DIFCS,CCO2S1,CCO2S2
!    3,CLSGL2(N6,N5,N4),TORTL,DISPN,XDPTH(N,N6,N5,N4)
!    4,CO2S2(N3,N2,N1),VOLWM(M,N3,N2,N1),VOLY(N3,N2,N1)
!    4,CO2S2(N6,N5,N4),VOLWM(M,N6,N5,N4),VOLY(N6,N5,N4)
!     WRITE(*,443)'ROXFLS',I,J,M,MM,N1,N2,N3,N4,N5,N6,N
!    2,ROXFLS(N,N6,N5,N4),RFLOXS,DFVOXS,DIFOS,COXYS1,COXYS2
!    3,OLSGL2(N6,N5,N4),TORTL,DISPN,XDPTH(N,N6,N5,N4)
!     WRITE(*,443)'RN3FLW',I,J,M,MM,N1,N2,N3,N4,N5,N6,N
!    2,RN3FLW(N,N6,N5,N4),RFLNH3,DFVNH3,DIFNH,CNH3S1,CNH3S2
!    3,ZHSGL2(N6,N5,N4),TORTL,DIFNH0,DIFNH1,DISPN,XDPTH(N,N6,N5,N4)
!    4,VFLW,ZNH3S2(N3,N2,N1),ZNH3S2(N6,N5,N4)
!     WRITE(*,443)'RH2PFS',I,J,M,MM,N1,N2,N3,N4,N5,N6,N
!    2,RH2PFS(N,N6,N5,N4),RFLPO4,DFVPO4,DIFPO,CPO4S1,CPO4S2
!    3,VLPO4(N3,N2,N1),VLPO4(N6,N5,N4),VOLW2A,VOLWPA
!    4,H2PO42(N3,N2,N1),H2PO42(N6,N5,N4)
!443   FORMAT(A8,11I4,20E12.4)
!     ENDIF
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
      DO 9755 K=0,jcplx1
      XOCFLS(K,N,N6,N5,N4)=XOCFLS(K,N,N6,N5,N4)+ROCFLS(K,N,N6,N5,N4)
      XONFLS(K,N,N6,N5,N4)=XONFLS(K,N,N6,N5,N4)+RONFLS(K,N,N6,N5,N4)
      XOPFLS(K,N,N6,N5,N4)=XOPFLS(K,N,N6,N5,N4)+ROPFLS(K,N,N6,N5,N4)
      XOAFLS(K,N,N6,N5,N4)=XOAFLS(K,N,N6,N5,N4)+ROAFLS(K,N,N6,N5,N4)
      XOCFHS(K,N,N6,N5,N4)=XOCFHS(K,N,N6,N5,N4)+ROCFHS(K,N,N6,N5,N4)
      XONFHS(K,N,N6,N5,N4)=XONFHS(K,N,N6,N5,N4)+RONFHS(K,N,N6,N5,N4)
      XOPFHS(K,N,N6,N5,N4)=XOPFHS(K,N,N6,N5,N4)+ROPFHS(K,N,N6,N5,N4)
      XOAFHS(K,N,N6,N5,N4)=XOAFHS(K,N,N6,N5,N4)+ROAFHS(K,N,N6,N5,N4)
9755  CONTINUE
      XCOFLS(N,N6,N5,N4)=XCOFLS(N,N6,N5,N4)+RCOFLS(N,N6,N5,N4)
      XCHFLS(N,N6,N5,N4)=XCHFLS(N,N6,N5,N4)+RCHFLS(N,N6,N5,N4)
      XOXFLS(N,N6,N5,N4)=XOXFLS(N,N6,N5,N4)+ROXFLS(N,N6,N5,N4)
      XNGFLS(N,N6,N5,N4)=XNGFLS(N,N6,N5,N4)+RNGFLS(N,N6,N5,N4)
      XN2FLS(N,N6,N5,N4)=XN2FLS(N,N6,N5,N4)+RN2FLS(N,N6,N5,N4)
      XHGFLS(N,N6,N5,N4)=XHGFLS(N,N6,N5,N4)+RHGFLS(N,N6,N5,N4)
      XN4FLW(N,N6,N5,N4)=XN4FLW(N,N6,N5,N4)+RN4FLW(N,N6,N5,N4)
      XN3FLW(N,N6,N5,N4)=XN3FLW(N,N6,N5,N4)+RN3FLW(N,N6,N5,N4)
      XNOFLW(N,N6,N5,N4)=XNOFLW(N,N6,N5,N4)+RNOFLW(N,N6,N5,N4)
      XNXFLS(N,N6,N5,N4)=XNXFLS(N,N6,N5,N4)+RNXFLS(N,N6,N5,N4)
      XH1PFS(N,N6,N5,N4)=XH1PFS(N,N6,N5,N4)+RH1PFS(N,N6,N5,N4)
      XH2PFS(N,N6,N5,N4)=XH2PFS(N,N6,N5,N4)+RH2PFS(N,N6,N5,N4)
      XN4FLB(N,N6,N5,N4)=XN4FLB(N,N6,N5,N4)+RN4FLB(N,N6,N5,N4)
      XN3FLB(N,N6,N5,N4)=XN3FLB(N,N6,N5,N4)+RN3FLB(N,N6,N5,N4)
      XNOFLB(N,N6,N5,N4)=XNOFLB(N,N6,N5,N4)+RNOFLB(N,N6,N5,N4)
      XNXFLB(N,N6,N5,N4)=XNXFLB(N,N6,N5,N4)+RNXFLB(N,N6,N5,N4)
      XH2BFB(N,N6,N5,N4)=XH2BFB(N,N6,N5,N4)+RH2BFB(N,N6,N5,N4)
      XCOFHS(N,N6,N5,N4)=XCOFHS(N,N6,N5,N4)+RCOFHS(N,N6,N5,N4)
      XCHFHS(N,N6,N5,N4)=XCHFHS(N,N6,N5,N4)+RCHFHS(N,N6,N5,N4)
      XOXFHS(N,N6,N5,N4)=XOXFHS(N,N6,N5,N4)+ROXFHS(N,N6,N5,N4)
      XNGFHS(N,N6,N5,N4)=XNGFHS(N,N6,N5,N4)+RNGFHS(N,N6,N5,N4)
      XN2FHS(N,N6,N5,N4)=XN2FHS(N,N6,N5,N4)+RN2FHS(N,N6,N5,N4)
      XHGFHS(N,N6,N5,N4)=XHGFHS(N,N6,N5,N4)+RHGFHS(N,N6,N5,N4)
      XN4FHW(N,N6,N5,N4)=XN4FHW(N,N6,N5,N4)+RN4FHW(N,N6,N5,N4)
      XN3FHW(N,N6,N5,N4)=XN3FHW(N,N6,N5,N4)+RN3FHW(N,N6,N5,N4)
      XNOFHW(N,N6,N5,N4)=XNOFHW(N,N6,N5,N4)+RNOFHW(N,N6,N5,N4)
      XNXFHS(N,N6,N5,N4)=XNXFHS(N,N6,N5,N4)+RNXFHS(N,N6,N5,N4)
      XH1PHS(N,N6,N5,N4)=XH1PHS(N,N6,N5,N4)+RH1PHS(N,N6,N5,N4)
      XH2PHS(N,N6,N5,N4)=XH2PHS(N,N6,N5,N4)+RH2PHS(N,N6,N5,N4)
      XN4FHB(N,N6,N5,N4)=XN4FHB(N,N6,N5,N4)+RN4FHB(N,N6,N5,N4)
      XN3FHB(N,N6,N5,N4)=XN3FHB(N,N6,N5,N4)+RN3FHB(N,N6,N5,N4)
      XNOFHB(N,N6,N5,N4)=XNOFHB(N,N6,N5,N4)+RNOFHB(N,N6,N5,N4)
      XNXFHB(N,N6,N5,N4)=XNXFHB(N,N6,N5,N4)+RNXFHB(N,N6,N5,N4)
      XH1BHB(N,N6,N5,N4)=XH1BHB(N,N6,N5,N4)+RH1BHB(N,N6,N5,N4)
      XH2BHB(N,N6,N5,N4)=XH2BHB(N,N6,N5,N4)+RH2BHB(N,N6,N5,N4)
!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
      IF(N.EQ.3)THEN
!
!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FINHM=macro-micropore water transfer from watsub.f
!     VOLWM,VOLWHM=micropore,macropore water volume
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
      IF(FINHM(M,N6,N5,N4).GT.0.0)THEN
      IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AMAX1(0.0,AMIN1(VFLWX,FINHM(M,N6,N5,N4) &
      /VOLWHM(M,N6,N5,N4)))
      ELSE
      VFLW=VFLWX
      ENDIF
      DO 9970 K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQCH2(K,N6,N5,N4))
      RFLON(K)=VFLW*AMAX1(0.0,OQNH2(K,N6,N5,N4))
      RFLOP(K)=VFLW*AMAX1(0.0,OQPH2(K,N6,N5,N4))
      RFLOA(K)=VFLW*AMAX1(0.0,OQAH2(K,N6,N5,N4))
9970  CONTINUE
      RFLCOS=VFLW*AMAX1(0.0,CO2SH2(N6,N5,N4))
      RFLCHS=VFLW*AMAX1(0.0,CH4SH2(N6,N5,N4))
      RFLOXS=VFLW*AMAX1(0.0,OXYSH2(N6,N5,N4))
      RFLNGS=VFLW*AMAX1(0.0,Z2GSH2(N6,N5,N4))
      RFLN2S=VFLW*AMAX1(0.0,Z2OSH2(N6,N5,N4))
      RFLHGS=VFLW*AMAX1(0.0,H2GSH2(N6,N5,N4))
      RFLNH4=VFLW*AMAX1(0.0,ZNH4H2(N6,N5,N4))*VLNH4(N6,N5,N4)
      RFLNH3=VFLW*AMAX1(0.0,ZNH3H2(N6,N5,N4))*VLNH4(N6,N5,N4)
      RFLNO3=VFLW*AMAX1(0.0,ZNO3H2(N6,N5,N4))*VLNO3(N6,N5,N4)
      RFLNO2=VFLW*AMAX1(0.0,ZNO2H2(N6,N5,N4))*VLNO3(N6,N5,N4)
      RFLP14=VFLW*AMAX1(0.0,H1P4H2(N6,N5,N4))*VLPO4(N6,N5,N4)
      RFLPO4=VFLW*AMAX1(0.0,H2P4H2(N6,N5,N4))*VLPO4(N6,N5,N4)
      RFLN4B=VFLW*AMAX1(0.0,ZN4BH2(N6,N5,N4))*VLNHB(N6,N5,N4)
      RFLN3B=VFLW*AMAX1(0.0,ZN3BH2(N6,N5,N4))*VLNHB(N6,N5,N4)
      RFLNOB=VFLW*AMAX1(0.0,ZNOBH2(N6,N5,N4))*VLNOB(N6,N5,N4)
      RFLN2B=VFLW*AMAX1(0.0,ZN2BH2(N6,N5,N4))*VLNOB(N6,N5,N4)
      RFLP1B=VFLW*AMAX1(0.0,H1PBH2(N6,N5,N4))*VLPOB(N6,N5,N4)
      RFLPOB=VFLW*AMAX1(0.0,H2PBH2(N6,N5,N4))*VLPOB(N6,N5,N4)
!
!     MICROPORE TO MACROPORE TRANSFER
!
      ELSEIF(FINHM(M,N6,N5,N4).LT.0.0)THEN
      IF(VOLWM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AMIN1(0.0,AMAX1(-VFLWX,FINHM(M,N6,N5,N4) &
      /VOLWM(M,N6,N5,N4)))
      ELSE
      VFLW=-VFLWX
      ENDIF
      DO 9965 K=0,jcplx1
      RFLOC(K)=VFLW*AMAX1(0.0,OQC2(K,N6,N5,N4))
      RFLON(K)=VFLW*AMAX1(0.0,OQN2(K,N6,N5,N4))
      RFLOP(K)=VFLW*AMAX1(0.0,OQP2(K,N6,N5,N4))
      RFLOA(K)=VFLW*AMAX1(0.0,OQA2(K,N6,N5,N4))
9965  CONTINUE
      RFLCOS=VFLW*AMAX1(0.0,CO2S2(N6,N5,N4))
      RFLCHS=VFLW*AMAX1(0.0,CH4S2(N6,N5,N4))
      RFLOXS=VFLW*AMAX1(0.0,OXYS2(N6,N5,N4))
      RFLNGS=VFLW*AMAX1(0.0,Z2GS2(N6,N5,N4))
      RFLN2S=VFLW*AMAX1(0.0,Z2OS2(N6,N5,N4))
      RFLHGS=VFLW*AMAX1(0.0,H2GS2(N6,N5,N4))
      RFLNH4=VFLW*AMAX1(0.0,ZNH4S2(N6,N5,N4))*VLNH4(N6,N5,N4)
      RFLNH3=VFLW*AMAX1(0.0,ZNH3S2(N6,N5,N4))*VLNH4(N6,N5,N4)
      RFLNO3=VFLW*AMAX1(0.0,ZNO3S2(N6,N5,N4))*VLNO3(N6,N5,N4)
      RFLNO2=VFLW*AMAX1(0.0,ZNO2S2(N6,N5,N4))*VLNO3(N6,N5,N4)
      RFLP14=VFLW*AMAX1(0.0,H1PO42(N6,N5,N4))*VLPO4(N6,N5,N4)
      RFLPO4=VFLW*AMAX1(0.0,H2PO42(N6,N5,N4))*VLPO4(N6,N5,N4)
      RFLN4B=VFLW*AMAX1(0.0,ZNH4B2(N6,N5,N4))*VLNHB(N6,N5,N4)
      RFLN3B=VFLW*AMAX1(0.0,ZNH3B2(N6,N5,N4))*VLNHB(N6,N5,N4)
      RFLNOB=VFLW*AMAX1(0.0,ZNO3B2(N6,N5,N4))*VLNOB(N6,N5,N4)
      RFLN2B=VFLW*AMAX1(0.0,ZNO2B2(N6,N5,N4))*VLNOB(N6,N5,N4)
      RFLP1B=VFLW*AMAX1(0.0,H1POB2(N6,N5,N4))*VLPOB(N6,N5,N4)
      RFLPOB=VFLW*AMAX1(0.0,H2POB2(N6,N5,N4))*VLPOB(N6,N5,N4)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
      ELSE
      DO 9960 K=0,jcplx1
      RFLOC(K)=0.0
      RFLON(K)=0.0
      RFLOP(K)=0.0
      RFLOA(K)=0.0
9960  CONTINUE
      RFLCOS=0.0
      RFLCHS=0.0
      RFLOXS=0.0
      RFLNGS=0.0
      RFLN2S=0.0
      RFLHGS=0.0
      RFLNH4=0.0
      RFLNH3=0.0
      RFLNO3=0.0
      RFLNO2=0.0
      RFLP14=0.0
      RFLPO4=0.0
      RFLN4B=0.0
      RFLN3B=0.0
      RFLNOB=0.0
      RFLN2B=0.0
      RFLP1B=0.0
      RFLPOB=0.0
      ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION
!     DIFFERENCES
!
!     VOLWM,VOLWHM=micropore,macropore water-filled porosity from watsub.f
!     DFV*S,DFV*B=diffusive solute flux between macro- and micropore in non-band,band
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     *2,*H2=solute content of micropores,macropores
!
      IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VOLWHS=AMIN1(XFRS*VOLT(N6,N5,N4),VOLWHM(M,N6,N5,N4))
      VOLWT=VOLWM(M,N6,N5,N4)+VOLWHS
      DO 9955 K=0,jcplx1
      DFVOC(K)=XNPH*(AMAX1(0.0,OQCH2(K,N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,OQC2(K,N6,N5,N4))*VOLWHS)/VOLWT
      DFVON(K)=XNPH*(AMAX1(0.0,OQNH2(K,N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,OQN2(K,N6,N5,N4))*VOLWHS)/VOLWT
      DFVOP(K)=XNPH*(AMAX1(0.0,OQPH2(K,N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,OQP2(K,N6,N5,N4))*VOLWHS)/VOLWT
      DFVOA(K)=XNPH*(AMAX1(0.0,OQAH2(K,N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,OQA2(K,N6,N5,N4))*VOLWHS)/VOLWT
9955  CONTINUE
      DFVCOS=XNPH*(AMAX1(0.0,CO2SH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,CO2S2(N6,N5,N4))*VOLWHS)/VOLWT
      DFVCHS=XNPH*(AMAX1(0.0,CH4SH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,CH4S2(N6,N5,N4))*VOLWHS)/VOLWT
      DFVOXS=XNPH*(AMAX1(0.0,OXYSH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,OXYS2(N6,N5,N4))*VOLWHS)/VOLWT
      DFVNGS=XNPH*(AMAX1(0.0,Z2GSH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,Z2GS2(N6,N5,N4))*VOLWHS)/VOLWT
      DFVN2S=XNPH*(AMAX1(0.0,Z2OSH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,Z2OS2(N6,N5,N4))*VOLWHS)/VOLWT
      DFVHGS=XNPH*(AMAX1(0.0,H2GSH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,H2GS2(N6,N5,N4))*VOLWHS)/VOLWT
      DFVNH4=XNPH*(AMAX1(0.0,ZNH4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNH4S2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNH4(N6,N5,N4)
      DFVNH3=XNPH*(AMAX1(0.0,ZNH3H2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNH3S2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNH4(N6,N5,N4)
      DFVNO3=XNPH*(AMAX1(0.0,ZNO3H2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNO3S2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNO3(N6,N5,N4)
      DFVNO2=XNPH*(AMAX1(0.0,ZNO2H2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNO2S2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNO3(N6,N5,N4)
      DFVP14=XNPH*(AMAX1(0.0,H1P4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,H1PO42(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLPO4(N6,N5,N4)
      DFVPO4=XNPH*(AMAX1(0.0,H2P4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,H2PO42(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLPO4(N6,N5,N4)
      DFVN4B=XNPH*(AMAX1(0.0,ZN4BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNH4B2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNHB(N6,N5,N4)
      DFVN3B=XNPH*(AMAX1(0.0,ZN3BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNH3B2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNHB(N6,N5,N4)
      DFVNOB=XNPH*(AMAX1(0.0,ZNOBH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNO3B2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNOB(N6,N5,N4)
      DFVN2B=XNPH*(AMAX1(0.0,ZN2BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,ZNO2B2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLNOB(N6,N5,N4)
      DFVP1B=XNPH*(AMAX1(0.0,H1PBH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,H1POB2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLPOB(N6,N5,N4)
      DFVPOB=XNPH*(AMAX1(0.0,H2PBH2(N6,N5,N4))*VOLWM(M,N6,N5,N4) &
      -AMAX1(0.0,H2POB2(N6,N5,N4))*VOLWHS)/VOLWT &
      *VLPOB(N6,N5,N4)
      ELSE
      DO 9975 K=0,jcplx1
      DFVOC(K)=0.0
      DFVON(K)=0.0
      DFVOP(K)=0.0
      DFVOA(K)=0.0
9975  CONTINUE
      DFVCOS=0.0
      DFVCHS=0.0
      DFVOXS=0.0
      DFVNGS=0.0
      DFVN2S=0.0
      DFVHGS=0.0
      DFVNH4=0.0
      DFVNH3=0.0
      DFVNO3=0.0
      DFVNO2=0.0
      DFVP14=0.0
      DFVPO4=0.0
      DFVN4B=0.0
      DFVN3B=0.0
      DFVNOB=0.0
      DFVN2B=0.0
      DFVP1B=0.0
      DFVPOB=0.0
      ENDIF
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
      DO 9950 K=0,jcplx1
      ROCFXS(K,N6,N5,N4)=RFLOC(K)+DFVOC(K)
      RONFXS(K,N6,N5,N4)=RFLON(K)+DFVON(K)
      ROPFXS(K,N6,N5,N4)=RFLOP(K)+DFVOP(K)
      ROAFXS(K,N6,N5,N4)=RFLOA(K)+DFVOA(K)
9950  CONTINUE
      RCOFXS(N6,N5,N4)=RFLCOS+DFVCOS
      RCHFXS(N6,N5,N4)=RFLCHS+DFVCHS
      ROXFXS(N6,N5,N4)=RFLOXS+DFVOXS
      RNGFXS(N6,N5,N4)=RFLNGS+DFVNGS
      RN2FXS(N6,N5,N4)=RFLN2S+DFVN2S
      RHGFXS(N6,N5,N4)=RFLHGS+DFVHGS
      RN4FXW(N6,N5,N4)=RFLNH4+DFVNH4
      RN3FXW(N6,N5,N4)=RFLNH3+DFVNH3
      RNOFXW(N6,N5,N4)=RFLNO3+DFVNO3
      RNXFXS(N6,N5,N4)=RFLNO2+DFVNO2
      RH1PXS(N6,N5,N4)=RFLP14+DFVP14
      RH2PXS(N6,N5,N4)=RFLPO4+DFVPO4
      RN4FXB(N6,N5,N4)=RFLN4B+DFVN4B
      RN3FXB(N6,N5,N4)=RFLN3B+DFVN3B
      RNOFXB(N6,N5,N4)=RFLNOB+DFVNOB
      RNXFXB(N6,N5,N4)=RFLN2B+DFVN2B
      RH1BXB(N6,N5,N4)=RFLP1B+DFVP1B
      RH2BXB(N6,N5,N4)=RFLPOB+DFVPOB
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FXS,X*FXB= hourly convective + diffusive solute flux between macro- and micropore in non-band,band
!     R*FXS,R*FXB=convective + diffusive solute flux between macro- and micropore in non-band,band
!
      DO 9945 K=0,jcplx1
      XOCFXS(K,N6,N5,N4)=XOCFXS(K,N6,N5,N4)+ROCFXS(K,N6,N5,N4)
      XONFXS(K,N6,N5,N4)=XONFXS(K,N6,N5,N4)+RONFXS(K,N6,N5,N4)
      XOPFXS(K,N6,N5,N4)=XOPFXS(K,N6,N5,N4)+ROPFXS(K,N6,N5,N4)
      XOAFXS(K,N6,N5,N4)=XOAFXS(K,N6,N5,N4)+ROAFXS(K,N6,N5,N4)
9945  CONTINUE
      XCOFXS(N6,N5,N4)=XCOFXS(N6,N5,N4)+RCOFXS(N6,N5,N4)
      XCHFXS(N6,N5,N4)=XCHFXS(N6,N5,N4)+RCHFXS(N6,N5,N4)
      XOXFXS(N6,N5,N4)=XOXFXS(N6,N5,N4)+ROXFXS(N6,N5,N4)
      XNGFXS(N6,N5,N4)=XNGFXS(N6,N5,N4)+RNGFXS(N6,N5,N4)
      XN2FXS(N6,N5,N4)=XN2FXS(N6,N5,N4)+RN2FXS(N6,N5,N4)
      XHGFXS(N6,N5,N4)=XHGFXS(N6,N5,N4)+RHGFXS(N6,N5,N4)
      XN4FXW(N6,N5,N4)=XN4FXW(N6,N5,N4)+RN4FXW(N6,N5,N4)
      XN3FXW(N6,N5,N4)=XN3FXW(N6,N5,N4)+RN3FXW(N6,N5,N4)
      XNOFXW(N6,N5,N4)=XNOFXW(N6,N5,N4)+RNOFXW(N6,N5,N4)
      XNXFXS(N6,N5,N4)=XNXFXS(N6,N5,N4)+RNXFXS(N6,N5,N4)
      XH1PXS(N6,N5,N4)=XH1PXS(N6,N5,N4)+RH1PXS(N6,N5,N4)
      XH2PXS(N6,N5,N4)=XH2PXS(N6,N5,N4)+RH2PXS(N6,N5,N4)
      XN4FXB(N6,N5,N4)=XN4FXB(N6,N5,N4)+RN4FXB(N6,N5,N4)
      XN3FXB(N6,N5,N4)=XN3FXB(N6,N5,N4)+RN3FXB(N6,N5,N4)
      XNOFXB(N6,N5,N4)=XNOFXB(N6,N5,N4)+RNOFXB(N6,N5,N4)
      XNXFXB(N6,N5,N4)=XNXFXB(N6,N5,N4)+RNXFXB(N6,N5,N4)
      XH1BXB(N6,N5,N4)=XH1BXB(N6,N5,N4)+RH1BXB(N6,N5,N4)
      XH2BXB(N6,N5,N4)=XH2BXB(N6,N5,N4)+RH2BXB(N6,N5,N4)
      ENDIF
      ENDIF
!
!     GASEOUS TRANSPORT FROM GASEOUS DIFFUSIVITY AND CONCENTRATION
!     DIFFERENCES BETWEEN ADJACENT GRID CELLS
!
!     THETPM,VOLPM=air-filled porosity,volume from watsub.f
!
      IF(THETPM(M,N3,N2,N1).GT.THETX &
      .AND.THETPM(M,N6,N5,N4).GT.THETX &
      .AND.VOLPM(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
      .AND.VOLPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
!
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
      DFLG2=2.0*AMAX1(0.0,THETPM(M,N3,N2,N1))*POROQ &
      *THETPM(M,N3,N2,N1)/POROS(N3,N2,N1) &
      *AREA(N,N3,N2,N1)/DLYR(N,N3,N2,N1)
      DFLGL=2.0*AMAX1(0.0,THETPM(M,N6,N5,N4))*POROQ &
      *THETPM(M,N6,N5,N4)/POROS(N6,N5,N4) &
      *AREA(N,N6,N5,N4)/DLYR(N,N6,N5,N4)
      CNDC1=DFLG2*CGSGL2(N3,N2,N1)
      CND41=DFLG2*CHSGL2(N3,N2,N1)
      CNDO1=DFLG2*OGSGL2(N3,N2,N1)
      CNDG1=DFLG2*ZGSGL2(N3,N2,N1)
      CND21=DFLG2*Z2SGL2(N3,N2,N1)
      CNDH1=DFLG2*ZHSGL2(N3,N2,N1)
      CNHG1=DFLG2*HGSGL2(N3,N2,N1)
      CNDC2=DFLGL*CGSGL2(N6,N5,N4)
      CND42=DFLGL*CHSGL2(N6,N5,N4)
      CNDO2=DFLGL*OGSGL2(N6,N5,N4)
      CNDG2=DFLGL*ZGSGL2(N6,N5,N4)
      CND22=DFLGL*Z2SGL2(N6,N5,N4)
      CNDH2=DFLGL*ZHSGL2(N6,N5,N4)
      CNHG2=DFLGL*HGSGL2(N6,N5,N4)
!
!     GASOUS CONDUCTANCES
!
!     D*G=gaseous diffusivity in soil
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
      DCO2G(N,N6,N5,N4)=(CNDC1*CNDC2)/(CNDC1+CNDC2)
      DCH4G(N,N6,N5,N4)=(CND41*CND42)/(CND41+CND42)
      DOXYG(N,N6,N5,N4)=(CNDO1*CNDO2)/(CNDO1+CNDO2)
      DZ2GG(N,N6,N5,N4)=(CNDG1*CNDG2)/(CNDG1+CNDG2)
      DZ2OG(N,N6,N5,N4)=(CND21*CND22)/(CND21+CND22)
      DNH3G(N,N6,N5,N4)=(CNDH1*CNDH2)/(CNDH1+CNDH2)
      DH2GG(N,N6,N5,N4)=(CNHG1*CNHG2)/(CNHG1+CNHG2)
!
!     GASEOUS CONCENTRATIONS FROM AIR-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
!
!     C*G1,C*G2=gaseous concentration in source,destination layer
!     *G2=gaseous content
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     VOLPM=air-filled porosity
!
      CCO2G1=AMAX1(0.0,CO2G2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      CCH4G1=AMAX1(0.0,CH4G2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      COXYG1=AMAX1(0.0,OXYG2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      CZ2GG1=AMAX1(0.0,Z2GG2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      CZ2OG1=AMAX1(0.0,Z2OG2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      CNH3G1=AMAX1(0.0,ZN3G2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      CH2GG1=AMAX1(0.0,H2GG2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      CCO2G2=AMAX1(0.0,CO2G2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
      CCH4G2=AMAX1(0.0,CH4G2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
      COXYG2=AMAX1(0.0,OXYG2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
      CZ2GG2=AMAX1(0.0,Z2GG2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
      CZ2OG2=AMAX1(0.0,Z2OG2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
      CNH3G2=AMAX1(0.0,ZN3G2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
      CH2GG2=AMAX1(0.0,H2GG2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
!
!     DIFFUSIVE GAS TRANSFER DRIVEN BY GAS CONCENTRATIONS IN
!     ADJACENT GRID CELLS
!
!     DFV*G=diffusive gas flux
!     C*G1,C*G2=gaseous concentration in source,destination layer
!
      DFVCOG=DCO2G(N,N6,N5,N4)*(CCO2G1-CCO2G2)
      DFVCHG=DCH4G(N,N6,N5,N4)*(CCH4G1-CCH4G2)
      DFVOXG=DOXYG(N,N6,N5,N4)*(COXYG1-COXYG2)
      DFVNGG=DZ2GG(N,N6,N5,N4)*(CZ2GG1-CZ2GG2)
      DFVN2G=DZ2OG(N,N6,N5,N4)*(CZ2OG1-CZ2OG2)
      DFVN3G=DNH3G(N,N6,N5,N4)*(CNH3G1-CNH3G2)
      DFVHGG=DH2GG(N,N6,N5,N4)*(CH2GG1-CH2GG2)
!
!     CONVECTIVE GAS TRANSFER DRIVEN BY SOIL WATER FLUXES
!     FROM 'WATSUB' AND GAS CONCENTRATIONS IN THE ADJACENT GRID CELLS
!     DEPENDING ON WATER FLUX DIRECTION
!
!     FLQM=total water flux into soil micropore+macropore from watsub.f
!     VOLPM=air-filled porosity
!     RFL*G=convective gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     *G2=gaseous content
!
      IF(FLQM(N,N6,N5,N4).GT.0.0)THEN
      IF(VOLPM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=-AMAX1(0.0,AMIN1(VFLWX,FLQM(N,N6,N5,N4) &
      /VOLPM(M,N6,N5,N4)))
      ELSE
      VFLW=-VFLWX
      ENDIF
      RFLCOG=VFLW*AMAX1(0.0,CO2G2(N6,N5,N4))
      RFLCHG=VFLW*AMAX1(0.0,CH4G2(N6,N5,N4))
      RFLOXG=VFLW*AMAX1(0.0,OXYG2(N6,N5,N4))
      RFLNGG=VFLW*AMAX1(0.0,Z2GG2(N6,N5,N4))
      RFLN2G=VFLW*AMAX1(0.0,Z2OG2(N6,N5,N4))
      RFLN3G=VFLW*AMAX1(0.0,ZN3G2(N6,N5,N4))
      RFLH2G=VFLW*AMAX1(0.0,H2GG2(N6,N5,N4))
      ELSE
      IF(VOLPM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=-AMIN1(0.0,AMAX1(-VFLWX,FLQM(N,N6,N5,N4) &
      /VOLPM(M,N3,N2,N1)))
      ELSE
      VFLW=VFLWX
      ENDIF
      RFLCOG=VFLW*AMAX1(0.0,CO2G2(N3,N2,N1))
      RFLCHG=VFLW*AMAX1(0.0,CH4G2(N3,N2,N1))
      RFLOXG=VFLW*AMAX1(0.0,OXYG2(N3,N2,N1))
      RFLNGG=VFLW*AMAX1(0.0,Z2GG2(N3,N2,N1))
      RFLN2G=VFLW*AMAX1(0.0,Z2OG2(N3,N2,N1))
      RFLN3G=VFLW*AMAX1(0.0,ZN3G2(N3,N2,N1))
      RFLH2G=VFLW*AMAX1(0.0,H2GG2(N3,N2,N1))
      ENDIF
!
!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE + CONVECTIVE FLUX
!
!     R*FLG=total convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!     DFV*G=diffusive gas flux
!     RFL*G=convective gas flux
!
      RCOFLG(N,N6,N5,N4)=DFVCOG+RFLCOG
      RCHFLG(N,N6,N5,N4)=DFVCHG+RFLCHG
      ROXFLG(N,N6,N5,N4)=DFVOXG+RFLOXG
      RNGFLG(N,N6,N5,N4)=DFVNGG+RFLNGG
      RN2FLG(N,N6,N5,N4)=DFVN2G+RFLN2G
      RN3FLG(N,N6,N5,N4)=DFVN3G+RFLN3G
      RHGFLG(N,N6,N5,N4)=DFVHGG+RFLH2G
!     IF(I.EQ.43)THEN
!     WRITE(*,3133)'ROXFL2',I,J,M,MM,N1,N2,N3,N,XOXFLG(N,N6,N5,N4)
!    2,ROXFLG(N,N6,N5,N4),DFVOXG,RFLOXG,COXYG1,COXYG2
!    3,OXYG2(N3,N2,N1),OXYG2(N6,N5,N4)
!    4,FLQM(N,N6,N5,N4),VFLW,DOXYG(N,N6,N5,N4)
!    5,THETPM(M,N3,N2,N1),THETPM(M,N6,N5,N4)
!    5,VOLPM(M,N3,N2,N1),VOLPM(M,N6,N5,N4)
!     WRITE(*,3133)'RNGFLG',I,J,M,MM,N4,N4,N6,N,RNGFLG(N,N6,N5,N4)
!    2,DFVNGG,RFLNGG,DZ2GG(N,N6,N5,N4),CZ2GG1,CZ2GG2
!3133  FORMAT(A8,8I4,20E12.4)
!     ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLG=hourly total convective+diffusive gas flux
!
      XCOFLG(N,N6,N5,N4)=XCOFLG(N,N6,N5,N4)+RCOFLG(N,N6,N5,N4)
      XCHFLG(N,N6,N5,N4)=XCHFLG(N,N6,N5,N4)+RCHFLG(N,N6,N5,N4)
      XOXFLG(N,N6,N5,N4)=XOXFLG(N,N6,N5,N4)+ROXFLG(N,N6,N5,N4)
      XNGFLG(N,N6,N5,N4)=XNGFLG(N,N6,N5,N4)+RNGFLG(N,N6,N5,N4)
      XN2FLG(N,N6,N5,N4)=XN2FLG(N,N6,N5,N4)+RN2FLG(N,N6,N5,N4)
      XN3FLG(N,N6,N5,N4)=XN3FLG(N,N6,N5,N4)+RN3FLG(N,N6,N5,N4)
      XHGFLG(N,N6,N5,N4)=XHGFLG(N,N6,N5,N4)+RHGFLG(N,N6,N5,N4)
      ELSE
      RCOFLG(N,N6,N5,N4)=0.0
      RCHFLG(N,N6,N5,N4)=0.0
      ROXFLG(N,N6,N5,N4)=0.0
      RNGFLG(N,N6,N5,N4)=0.0
      RN2FLG(N,N6,N5,N4)=0.0
      RN3FLG(N,N6,N5,N4)=0.0
      RHGFLG(N,N6,N5,N4)=0.0
      ENDIF
!
!     VOLATILIZATION-DISSOLUTION OF GASES IN SOIL
!     LAYER FROM GASEOUS CONCENTRATIONS VS. THEIR AQUEOUS
!     EQUIVALENTS DEPENDING ON SOLUBILITY FROM 'HOUR1'
!     AND TRANSFER COEFFICIENT 'DFGS' FROM 'WATSUB'
!
!     THETPM,VOLWPM=air-filled porosity,volume
!     R*DFG=water-air gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     DFGS=rate constant for air-water gas exchange from watsub.f
!     *G2,*S2=gaseous,aqueous gas content
!     VOLW*=equivalent aqueous volume for gas
!
      IF(N.EQ.3)THEN
      IF(THETPM(M,N6,N5,N4).GT.THETX)THEN
      RCODFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,CO2G2(N6,N5,N4))*VOLWCO(N6,N5,N4) &
      -CO2S2(N6,N5,N4)*VOLPM(M,N6,N5,N4)) &
      /(VOLWCO(N6,N5,N4)+VOLPM(M,N6,N5,N4))
      RCHDFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,CH4G2(N6,N5,N4))*VOLWCH(N6,N5,N4) &
      -CH4S2(N6,N5,N4)*VOLPM(M,N6,N5,N4)) &
      /(VOLWCH(N6,N5,N4)+VOLPM(M,N6,N5,N4))
      ROXDFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,OXYG2(N6,N5,N4))*VOLWOX(N6,N5,N4) &
      -OXYS2(N6,N5,N4)*VOLPM(M,N6,N5,N4)) &
      /(VOLWOX(N6,N5,N4)+VOLPM(M,N6,N5,N4))
      RNGDFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,Z2GG2(N6,N5,N4))*VOLWNG(N6,N5,N4) &
      -Z2GS2(N6,N5,N4)*VOLPM(M,N6,N5,N4)) &
      /(VOLWNG(N6,N5,N4)+VOLPM(M,N6,N5,N4))
      RN2DFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,Z2OG2(N6,N5,N4))*VOLWN2(N6,N5,N4) &
      -Z2OS2(N6,N5,N4)*VOLPM(M,N6,N5,N4)) &
      /(VOLWN2(N6,N5,N4)+VOLPM(M,N6,N5,N4))
      IF(VOLPMA(N6,N5,N4).GT.ZEROS2(N5,N4) &
      .AND.VOLWXA(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      RN3DFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,ZN3G2(N6,N5,N4))*VOLWN3(N6,N5,N4) &
      -ZNH3S2(N6,N5,N4)*VOLPMA(N6,N5,N4)) &
      /(VOLWN3(N6,N5,N4)+VOLPMA(N6,N5,N4))
      CNH3S0=AMAX1(0.0,(ZNH3S2(N6,N5,N4)+RN3DFG(N6,N5,N4)) &
      /VOLWXA(N6,N5,N4))
      CNH4S0=AMAX1(0.0,ZNH4S2(N6,N5,N4)) &
      /VOLWXA(N6,N5,N4)
      ELSE
      RN3DFG(N6,N5,N4)=0.0
      ENDIF
      IF(VOLPMB(N6,N5,N4).GT.ZEROS2(N5,N4) &
      .AND.VOLWXB(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      RNBDFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,ZN3G2(N6,N5,N4))*VOLWNB(N6,N5,N4) &
      -ZNH3B2(N6,N5,N4)*VOLPMB(N6,N5,N4)) &
      /(VOLWNB(N6,N5,N4)+VOLPMB(N6,N5,N4))
      CNH3B0=AMAX1(0.0,(ZNH3B2(N6,N5,N4)+RNBDFG(N6,N5,N4)) &
      /VOLWXB(N6,N5,N4))
      CNH4B0=AMAX1(0.0,ZNH4B2(N6,N5,N4))/VOLWXB(N6,N5,N4)
      ELSE
      RNBDFG(N6,N5,N4)=0.0
      ENDIF
      RHGDFG(N6,N5,N4)=DFGS(M,N6,N5,N4)*(AMAX1(ZEROS(N5,N4) &
      ,H2GG2(N6,N5,N4))*VOLWHG(N6,N5,N4) &
      -H2GS2(N6,N5,N4)*VOLPM(M,N6,N5,N4)) &
      /(VOLWHG(N6,N5,N4)+VOLPM(M,N6,N5,N4))
!     IF(I.EQ.121.AND.N6.EQ.2)THEN
!     WRITE(*,6666)'RCODFG',I,J,M,MM,N4,N5,N6,RCODFG(N6,N5,N4)
!    2,DFGS(M,N6,N5,N4),CO2G2(N6,N5,N4),VOLWCO(N6,N5,N4)
!    3,CO2S2(N6,N5,N4),VOLWM(M,N6,N5,N4),THETPM(M,N6,N5,N4)
!    4,SCO2L(N6,N5,N4),XCODFG(N6,N5,N4)
!     WRITE(*,6666)'RN3DFG',I,J,M,MM,N4,N5,N6,RN3DFG(N6,N5,N4)
!    2,DFGS(M,N6,N5,N4),ZN3S2A,VOLWN3(N6,N5,N4),ZNH3S2(N6,N5,N4)
!    3,VOLPMA(N6,N5,N4),RNBDFG(N6,N5,N4),ZN3S2B
!    4,VOLWNB(N6,N5,N4),ZNH3B2(N6,N5,N4),VOLPMB(N6,N5,N4)
!     WRITE(*,6666)'RCHDFG',I,J,M,MM,N4,N5,N6,RCHDFG(N6,N5,N4)
!    2,DFGS(M,N6,N5,N4),CH4G2(N6,N5,N4),VOLWCH(N6,N5,N4)
!    3,CH4S2(N6,N5,N4),VOLWM(M,N6,N5,N4),THETPM(M,N6,N5,N4)
!    4,SCH4L(N6,N5,N4),XCHDFG(N6,N5,N4)
!     WRITE(*,6666)'RNGDFG',I,J,M,MM,N4,N5,N6
!    2,RNGDFG(N6,N5,N4),DFGS(M,N6,N5,N4),Z2GG2(N6,N5,N4)
!    3,VOLWNG(N6,N5,N4),Z2GS2(N6,N5,N4),VOLPM(M,N6,N5,N4)
!    4,VOLWNG(N6,N5,N4),VOLPM(M,N6,N5,N4)
!6666  FORMAT(A8,7I4,20E12.4)
!     ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*DFG=hourly water-air gas flux
!
      XCODFG(N6,N5,N4)=XCODFG(N6,N5,N4)+RCODFG(N6,N5,N4)
      XCHDFG(N6,N5,N4)=XCHDFG(N6,N5,N4)+RCHDFG(N6,N5,N4)
      XOXDFG(N6,N5,N4)=XOXDFG(N6,N5,N4)+ROXDFG(N6,N5,N4)
      XNGDFG(N6,N5,N4)=XNGDFG(N6,N5,N4)+RNGDFG(N6,N5,N4)
      XN2DFG(N6,N5,N4)=XN2DFG(N6,N5,N4)+RN2DFG(N6,N5,N4)
      XN3DFG(N6,N5,N4)=XN3DFG(N6,N5,N4)+RN3DFG(N6,N5,N4)
      XNBDFG(N6,N5,N4)=XNBDFG(N6,N5,N4)+RNBDFG(N6,N5,N4)
      XHGDFG(N6,N5,N4)=XHGDFG(N6,N5,N4)+RHGDFG(N6,N5,N4)
      ELSE
      RCODFG(N6,N5,N4)=0.0
      RCHDFG(N6,N5,N4)=0.0
      ROXDFG(N6,N5,N4)=0.0
      RNGDFG(N6,N5,N4)=0.0
      RN2DFG(N6,N5,N4)=0.0
      RN3DFG(N6,N5,N4)=0.0
      RNBDFG(N6,N5,N4)=0.0
      RHGDFG(N6,N5,N4)=0.0
      ENDIF
      ENDIF
      ELSEIF(N.NE.3)THEN
      DCO2G(N,N6,N5,N4)=0.0
      DCH4G(N,N6,N5,N4)=0.0
      DOXYG(N,N6,N5,N4)=0.0
      DZ2GG(N,N6,N5,N4)=0.0
      DZ2OG(N,N6,N5,N4)=0.0
      DNH3G(N,N6,N5,N4)=0.0
      DH2GG(N,N6,N5,N4)=0.0
      DO 9750 K=0,jcplx1
      ROCFLS(K,N,N6,N5,N4)=0.0
      RONFLS(K,N,N6,N5,N4)=0.0
      ROPFLS(K,N,N6,N5,N4)=0.0
      ROAFLS(K,N,N6,N5,N4)=0.0
      ROCFHS(K,N,N6,N5,N4)=0.0
      RONFHS(K,N,N6,N5,N4)=0.0
      ROPFHS(K,N,N6,N5,N4)=0.0
      ROAFHS(K,N,N6,N5,N4)=0.0
9750  CONTINUE
      RCOFLS(N,N6,N5,N4)=0.0
      RCHFLS(N,N6,N5,N4)=0.0
      ROXFLS(N,N6,N5,N4)=0.0
      RNGFLS(N,N6,N5,N4)=0.0
      RN2FLS(N,N6,N5,N4)=0.0
      RHGFLS(N,N6,N5,N4)=0.0
      RN4FLW(N,N6,N5,N4)=0.0
      RN3FLW(N,N6,N5,N4)=0.0
      RNOFLW(N,N6,N5,N4)=0.0
      RNXFLS(N,N6,N5,N4)=0.0
      RH1PFS(N,N6,N5,N4)=0.0
      RH2PFS(N,N6,N5,N4)=0.0
      RN4FLB(N,N6,N5,N4)=0.0
      RN3FLB(N,N6,N5,N4)=0.0
      RNOFLB(N,N6,N5,N4)=0.0
      RNXFLB(N,N6,N5,N4)=0.0
      RH2BFB(N,N6,N5,N4)=0.0
      RCOFHS(N,N6,N5,N4)=0.0
      RCHFHS(N,N6,N5,N4)=0.0
      ROXFHS(N,N6,N5,N4)=0.0
      RNGFHS(N,N6,N5,N4)=0.0
      RN2FHS(N,N6,N5,N4)=0.0
      RHGFHS(N,N6,N5,N4)=0.0
      RN4FHW(N,N6,N5,N4)=0.0
      RN3FHW(N,N6,N5,N4)=0.0
      RNOFHW(N,N6,N5,N4)=0.0
      RNXFHS(N,N6,N5,N4)=0.0
      RH1PHS(N,N6,N5,N4)=0.0
      RH2PHS(N,N6,N5,N4)=0.0
      RN4FHB(N,N6,N5,N4)=0.0
      RN3FHB(N,N6,N5,N4)=0.0
      RNOFHB(N,N6,N5,N4)=0.0
      RNXFHB(N,N6,N5,N4)=0.0
      RH1BHB(N,N6,N5,N4)=0.0
      RH2BHB(N,N6,N5,N4)=0.0
      RCOFLG(N,N6,N5,N4)=0.0
      RCHFLG(N,N6,N5,N4)=0.0
      ROXFLG(N,N6,N5,N4)=0.0
      RNGFLG(N,N6,N5,N4)=0.0
      RN2FLG(N,N6,N5,N4)=0.0
      RN3FLG(N,N6,N5,N4)=0.0
      RHGFLG(N,N6,N5,N4)=0.0
      ENDIF
      ELSE
      DCO2G(N,N3,N2,N1)=0.0
      DCH4G(N,N3,N2,N1)=0.0
      DOXYG(N,N3,N2,N1)=0.0
      DZ2GG(N,N3,N2,N1)=0.0
      DZ2OG(N,N3,N2,N1)=0.0
      DNH3G(N,N3,N2,N1)=0.0
      DH2GG(N,N3,N2,N1)=0.0
      DO 9751 K=0,jcplx1
      ROCFLS(K,N,N3,N2,N1)=0.0
      RONFLS(K,N,N3,N2,N1)=0.0
      ROPFLS(K,N,N3,N2,N1)=0.0
      ROAFLS(K,N,N3,N2,N1)=0.0
      ROCFHS(K,N,N3,N2,N1)=0.0
      RONFHS(K,N,N3,N2,N1)=0.0
      ROPFHS(K,N,N3,N2,N1)=0.0
      ROAFHS(K,N,N3,N2,N1)=0.0
9751  CONTINUE
      RCOFLS(N,N3,N2,N1)=0.0
      RCHFLS(N,N3,N2,N1)=0.0
      ROXFLS(N,N3,N2,N1)=0.0
      RNGFLS(N,N3,N2,N1)=0.0
      RN2FLS(N,N3,N2,N1)=0.0
      RHGFLS(N,N3,N2,N1)=0.0
      RN4FLW(N,N3,N2,N1)=0.0
      RN3FLW(N,N3,N2,N1)=0.0
      RNOFLW(N,N3,N2,N1)=0.0
      RNXFLS(N,N3,N2,N1)=0.0
      RH1PFS(N,N3,N2,N1)=0.0
      RH2PFS(N,N3,N2,N1)=0.0
      RN4FLB(N,N3,N2,N1)=0.0
      RN3FLB(N,N3,N2,N1)=0.0
      RNOFLB(N,N3,N2,N1)=0.0
      RNXFLB(N,N3,N2,N1)=0.0
      RH2BFB(N,N3,N2,N1)=0.0
      RCOFHS(N,N3,N2,N1)=0.0
      RCHFHS(N,N3,N2,N1)=0.0
      ROXFHS(N,N3,N2,N1)=0.0
      RNGFHS(N,N3,N2,N1)=0.0
      RN2FHS(N,N3,N2,N1)=0.0
      RHGFHS(N,N3,N2,N1)=0.0
      RN4FHW(N,N3,N2,N1)=0.0
      RN3FHW(N,N3,N2,N1)=0.0
      RNOFHW(N,N3,N2,N1)=0.0
      RNXFHS(N,N3,N2,N1)=0.0
      RH1PHS(N,N3,N2,N1)=0.0
      RH2PHS(N,N3,N2,N1)=0.0
      RN4FHB(N,N3,N2,N1)=0.0
      RN3FHB(N,N3,N2,N1)=0.0
      RNOFHB(N,N3,N2,N1)=0.0
      RNXFHB(N,N3,N2,N1)=0.0
      RH1BHB(N,N3,N2,N1)=0.0
      RH2BHB(N,N3,N2,N1)=0.0
      RCOFLG(N,N3,N2,N1)=0.0
      RCHFLG(N,N3,N2,N1)=0.0
      ROXFLG(N,N3,N2,N1)=0.0
      RNGFLG(N,N3,N2,N1)=0.0
      RN2FLG(N,N3,N2,N1)=0.0
      RN3FLG(N,N3,N2,N1)=0.0
      RHGFLG(N,N3,N2,N1)=0.0
      ENDIF
120   CONTINUE
!
!     CHECK FOR BUBBLING IF THE SUM OF ALL GASEOUS EQUIVALENT
!     PARTIAL CONCENTRATIONS EXCEEDS ATMOSPHERIC PRESSURE
!
!     VOLWM=micropore water-filled porosity from watsub.f
!     VOLY=micropore volume
!     IFLGB=bubbling flag:0=enabled,1=disabled
!     S*L=solubility of gas in water from hour1.f
!
      IF(N3.GE.NUM(N2,N1).AND.M.NE.MX)THEN
      THETW1(N3,N2,N1)=AMAX1(0.0_r8,safe_adb(VOLWM(M,N3,N2,N1),VOLY(N3,N2,N1)))
      IF(THETW1(N3,N2,N1).GT.THETY(N3,N2,N1).AND.IFLGB.EQ.0)THEN
      SCO2X=12.0*SCO2L(N3,N2,N1)
      SCH4X=12.0*SCH4L(N3,N2,N1)
      SOXYX=32.0*SOXYL(N3,N2,N1)
      SN2GX=28.0*SN2GL(N3,N2,N1)
      SN2OX=28.0*SN2OL(N3,N2,N1)
      SNH3X=14.0*SNH3L(N3,N2,N1)
      SH2GX=2.0*SH2GL(N3,N2,N1)
!
!     GASEOUS EQUIVALENT PARTIAL CONCENTRATIONS
!
!     V*G2=molar gas concentration
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!     R*DFS=gas exchange between atmosphere and soil surface water
!
      IF(N3.EQ.NU(N2,N1))THEN
      VCO2G2=(CO2S2(N3,N2,N1)+RCODFS(NY,NX))/SCO2X
      VCH4G2=(CH4S2(N3,N2,N1)+RCHDFS(NY,NX))/SCH4X
      VOXYG2=(OXYS2(N3,N2,N1)+ROXDFS(NY,NX))/SOXYX
      VZ2GG2=(Z2GS2(N3,N2,N1)+RNGDFS(NY,NX))/SN2GX
      VZ2OG2=(Z2OS2(N3,N2,N1)+RN2DFS(NY,NX))/SN2OX
      VNH3G2=(ZNH3S2(N3,N2,N1)+RN3DFS(NY,NX))/SNH3X
      VNHBG2=(ZNH3B2(N3,N2,N1)+RNBDFS(NY,NX))/SNH3X
      VH2GG2=(H2GS2(N3,N2,N1)+RHGDFS(NY,NX))/SH2GX
      ELSE
      VCO2G2=CO2S2(N3,N2,N1)/SCO2X
      VCH4G2=CH4S2(N3,N2,N1)/SCH4X
      VOXYG2=OXYS2(N3,N2,N1)/SOXYX
      VZ2GG2=Z2GS2(N3,N2,N1)/SN2GX
      VZ2OG2=Z2OS2(N3,N2,N1)/SN2OX
      VNH3G2=ZNH3S2(N3,N2,N1)/SNH3X
      VNHBG2=ZNH3B2(N3,N2,N1)/SNH3X
      VH2GG2=H2GS2(N3,N2,N1)/SH2GX
      ENDIF
!
!     GASEOUS EQUIVALENT ATMOSPHERIC CONCENTRATION
!
!     VTATM=molar gas concentration at atmospheric pressure
!     VTGAS=total molar gas concentration
!
      VTATM=AMAX1(0.0,1.2194E+04*VOLWM(M,N3,N2,N1)/TKS(N3,N2,N1))
      VTGAS=VCO2G2+VCH4G2+VOXYG2+VZ2GG2+VZ2OG2+VNH3G2+VNHBG2+VH2GG2
!
!     PROPORTIONAL REMOVAL OF EXCESS AQUEOUS GASES
!
!     R*BBL=bubble flux
!     gas code:*CO*=CO2,*CH*=CH4,*OX*=O2,*NG*=N2,*N2*=N2O,*N3*=NH3,*HG*=H2
!     V*G2=molar gas concentration
!
      IF(VTGAS.GT.VTATM)THEN
      DVTGAS=0.5*(VTATM-VTGAS)
      RCOBBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VCO2G2/VTGAS)*SCO2X
      RCHBBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VCH4G2/VTGAS)*SCH4X
      ROXBBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VOXYG2/VTGAS)*SOXYX
      RNGBBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VZ2GG2/VTGAS)*SN2GX
      RN2BBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VZ2OG2/VTGAS)*SN2OX
      RN3BBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VNH3G2/VTGAS)*SNH3X
      RNBBBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VNHBG2/VTGAS)*SNH3X
      RHGBBL(N3,N2,N1)=AMIN1(0.0,DVTGAS*VH2GG2/VTGAS)*SH2GX
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*BBL=hourly bubble flux
!
      XCOBBL(N3,N2,N1)=XCOBBL(N3,N2,N1)+RCOBBL(N3,N2,N1)
      XCHBBL(N3,N2,N1)=XCHBBL(N3,N2,N1)+RCHBBL(N3,N2,N1)
      XOXBBL(N3,N2,N1)=XOXBBL(N3,N2,N1)+ROXBBL(N3,N2,N1)
      XNGBBL(N3,N2,N1)=XNGBBL(N3,N2,N1)+RNGBBL(N3,N2,N1)
      XN2BBL(N3,N2,N1)=XN2BBL(N3,N2,N1)+RN2BBL(N3,N2,N1)
      XN3BBL(N3,N2,N1)=XN3BBL(N3,N2,N1)+RN3BBL(N3,N2,N1)
      XNBBBL(N3,N2,N1)=XNBBBL(N3,N2,N1)+RNBBBL(N3,N2,N1)
      XHGBBL(N3,N2,N1)=XHGBBL(N3,N2,N1)+RHGBBL(N3,N2,N1)
      ELSE
      RCOBBL(N3,N2,N1)=0.0
      RCHBBL(N3,N2,N1)=0.0
      ROXBBL(N3,N2,N1)=0.0
      RNGBBL(N3,N2,N1)=0.0
      RN2BBL(N3,N2,N1)=0.0
      RN3BBL(N3,N2,N1)=0.0
      RNBBBL(N3,N2,N1)=0.0
      RHGBBL(N3,N2,N1)=0.0
      ENDIF
      ELSE
      IFLGB=1
      RCOBBL(N3,N2,N1)=0.0
      RCHBBL(N3,N2,N1)=0.0
      ROXBBL(N3,N2,N1)=0.0
      RNGBBL(N3,N2,N1)=0.0
      RN2BBL(N3,N2,N1)=0.0
      RN3BBL(N3,N2,N1)=0.0
      RNBBBL(N3,N2,N1)=0.0
      RHGBBL(N3,N2,N1)=0.0
      ENDIF

      ENDIF
125   CONTINUE
      end subroutine SoluteFluxAdjacentCells

end module InsideTranspMod
