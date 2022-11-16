module LateralTranspMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  USE abortutils, only : endrun
  use GridDataType
  use TFlxTypeMod
  use TracerIDMod
  use SoilBGCDataType
  use SoilWaterDataType
  use ChemTranspDataType
  use SnowDataType
  USE AqueChemDatatype
  USE EcoSIMCtrlDataType
  USE SedimentDataType
  USE MicrobialDataType
  USE SurfSoilDataType
  use SoilPropertyDataType
  use FlagDataType
  USE SoilHeatDataType
  use EcoSimConst
  use EcoSiMParDataMod, only : micpar
  use minimathmod , only : AZMAX1
  use EcoSIMConfig, only : jcplx => jcplxc,NFGs=>NFGsc
  use EcoSIMConfig, only : nlbiomcp=>nlbiomcpc
implicit none
  private
  character(len=*), parameter :: mod_filename = __FILE__

  public :: LateralTranspt
  contains

  subroutine LateralTranspt(I,J,NY,NX,LG)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(out) :: LG

  integer :: L,N,LX
  integer :: N1,N2,N3,N4,N5,N6
  integer :: N4B,N5B
  real(r8) :: VTATM,VTGAS,VH2GG2
  real(r8) :: VNH3G2,VZ2OG2
  real(r8) :: VCO2G2,VCH4G2,VOXYG2,VZ2GG2

! begin_execution
  call ZeroFluxArrays(NY,NX)
  call ZeroFluxAccumulators(NY,NX)

  LG=0
  LX=0

  D8575: DO L=NU(NY,NX),NL(NY,NX)
    !
    !     IDENTIFY LAYERS FOR BUBBLE FLUX TRANSFER
    !
    !     LG=lowest soil layer with gas phase
    !     V*G2=molar gas content
    !     *G=soil gas content
    !     VOLP=soil air-filled porosity
    !     VTATM=molar gas content at atmospheric pressure
    !     VTGAS=total molar gas contest
    !     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
    !             :*ZN3*=NH3,*H2G*=H2
!
    VCO2G2=trc_gasml(idg_CO2,L,NY,NX)/catomw
    VCH4G2=trc_gasml(idg_CH4,L,NY,NX)/catomw
    VOXYG2=trc_gasml(idg_O2,L,NY,NX)/32.0_r8
    VZ2GG2=trc_gasml(idg_N2,L,NY,NX)/28.0_r8
    VZ2OG2=trc_gasml(idg_N2O,L,NY,NX)/28.0_r8
    VNH3G2=trc_gasml(idg_NH3,L,NY,NX)/natomw
    VH2GG2=trc_gasml(idg_H2,L,NY,NX)/2.0_r8
    VTATM=AZMAX1(1.2194E+04*VOLP(L,NY,NX)/TKS(L,NY,NX))
    VTGAS=VCO2G2+VCH4G2+VOXYG2+VZ2GG2+VZ2OG2+VNH3G2+VH2GG2
    IF(THETP(L,NY,NX).LT.THETX.OR.VTGAS.GT.VTATM)LX=1
    IF(THETP(L,NY,NX).GE.THETX.AND.LX.EQ.0)LG=L

    VOLW1(L,NY,NX)=VOLW(L,NY,NX)
    VOLI1(L,NY,NX)=VOLI(L,NY,NX)
    VOLWH1(L,NY,NX)=VOLWH(L,NY,NX)
    VOLIH1(L,NY,NX)=VOLIH(L,NY,NX)

!
  !     NET WATER, HEAT, GAS, SOLUTE, SEDIMENT FLUX
  !
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !
    N1=NX
    N2=NY
    N3=L
    D8580: DO N=NCN(NY,NX),3
      IF(N.EQ.1)THEN
        N4=NX+1
        N5=NY
        N4B=NX-1
        N5B=NY
        N6=L
      ELSEIF(N.EQ.2)THEN
        N4=NX
        N5=NY+1
        N4B=NX
        N5B=NY-1
        N6=L
      ELSEIF(N.EQ.3)THEN
        N4=NX
        N5=NY
        N6=L+1
      ENDIF
!
!
      IF(L.EQ.NUM(N2,N1))THEN
        IF(N.NE.3)THEN

          call FluxesFromRunoff(N,N1,N2,N4,N5,N4B,N5B)
          !
          !     NET SALT FLUXES FROM RUNOFF AND SNOWPACK
          call SaltThruFluxRunoffAndSnowpack(N,N1,N2,N4,N5,N4B,N5B)
          !
    !
        ELSEIF(N.EQ.3)THEN
          call SaltFromRunoffSnowpack(N1,N2,NY,NX)
        ENDIF
!
        ! TOTAL FLUXES FROM SEDIMENT TRANSPORT
        call TotalFluxFromSedmentTransp(N,N1,N2,N4,N5,N4B,N5B,NY,NX)
      ENDIF
!
      call FluxBetweenGrids(N,N1,N2,N3,N4,N5,N6,NY,NX)
    ENDDO D8580
!
!     NET FREEZE-THAW
!
!     TTHAW,TTHAWH=net freeze-thaw flux in micropores,macropores
!     THTHAW=net freeze-thaw latent heat flux
!     THAW,THAWH=freeze-thaw flux in micropores,macropores from watsub.f
!     HTHAW=freeze-thaw latent heat flux from watsub.f
  !
    TTHAW(N3,N2,N1)=TTHAW(N3,N2,N1)+THAW(N3,N2,N1)
    TTHAWH(N3,N2,N1)=TTHAWH(N3,N2,N1)+THAWH(N3,N2,N1)
    THTHAW(N3,N2,N1)=THTHAW(N3,N2,N1)+HTHAW(N3,N2,N1)
  ENDDO D8575
  end subroutine LateralTranspt

!------------------------------------------------------------------------------------------

  subroutine ZeroFluxArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,K,M,NO,NGL
!     begin_execution
!
!     INITIALIZE NET WATER AND HEAT FLUXES FOR RUNOFF
!
  TQR(NY,NX)=0.0_r8
  THQR(NY,NX)=0.0_r8
  TQS(NY,NX)=0.0_r8
  TQW(NY,NX)=0.0_r8
  TQI(NY,NX)=0.0_r8
  THQS(NY,NX)=0.0_r8
!
!     INITIALIZE NET SOLUTE AND GAS FLUXES FOR RUNOFF
!
  D9960: DO K=1,micpar%n_litrsfk
    TOCQRS(K,NY,NX)=0.0_r8
    TONQRS(K,NY,NX)=0.0_r8
    TOPQRS(K,NY,NX)=0.0_r8
    TOAQRS(K,NY,NX)=0.0_r8
  ENDDO D9960
  TCOQRS(NY,NX)=0.0_r8
  TCHQRS(NY,NX)=0.0_r8
  TOXQRS(NY,NX)=0.0_r8
  TNGQRS(NY,NX)=0.0_r8
  TN2QRS(NY,NX)=0.0_r8
  THGQRS(NY,NX)=0.0_r8
  TN4QRS(NY,NX)=0.0_r8
  TN3QRS(NY,NX)=0.0_r8
  TNOQRS(NY,NX)=0.0_r8
  TNXQRS(NY,NX)=0.0_r8
  TP1QRS(NY,NX)=0.0_r8
  TPOQRS(NY,NX)=0.0_r8
  trcg_QSS(idg_beg:idg_end-1,NY,NX)=0.0_r8

  trcn_QSS(ids_nut_beg:ids_nuts_end,NY,NX)=0.0_r8
  DO  L=1,JS
    trcg_TBLS(idg_beg:idg_end-1,L,NY,NX)=0.0_r8
    trcn_TBLS(ids_nut_beg:ids_nuts_end,L,NY,NX)=0.0_r8
  ENDDO

  IF(ISALTG.NE.0)THEN
    TQRAL(NY,NX)=0.0_r8
    TQRFE(NY,NX)=0.0_r8
    TQRHY(NY,NX)=0.0_r8
    TQRCA(NY,NX)=0.0_r8
    TQRMG(NY,NX)=0.0_r8
    TQRNA(NY,NX)=0.0_r8
    TQRKA(NY,NX)=0.0_r8
    TQROH(NY,NX)=0.0_r8
    TQRSO(NY,NX)=0.0_r8
    TQRCL(NY,NX)=0.0_r8
    TQRC3(NY,NX)=0.0_r8
    TQRHC(NY,NX)=0.0_r8
    TQRAL1(NY,NX)=0.0_r8
    TQRAL2(NY,NX)=0.0_r8
    TQRAL3(NY,NX)=0.0_r8
    TQRAL4(NY,NX)=0.0_r8
    TQRALS(NY,NX)=0.0_r8
    TQRFE1(NY,NX)=0.0_r8
    TQRFE2(NY,NX)=0.0_r8
    TQRFE3(NY,NX)=0.0_r8
    TQRFE4(NY,NX)=0.0_r8
    TQRFES(NY,NX)=0.0_r8
    TQRCAO(NY,NX)=0.0_r8
    TQRCAC(NY,NX)=0.0_r8
    TQRCAH(NY,NX)=0.0_r8
    TQRCAS(NY,NX)=0.0_r8
    TQRMGO(NY,NX)=0.0_r8
    TQRMGC(NY,NX)=0.0_r8
    TQRMGH(NY,NX)=0.0_r8
    TQRMGS(NY,NX)=0.0_r8
    TQRNAC(NY,NX)=0.0_r8
    TQRNAS(NY,NX)=0.0_r8
    TQRKAS(NY,NX)=0.0_r8
    TQRH0P(NY,NX)=0.0_r8
    TQRH3P(NY,NX)=0.0_r8
    TQRF1P(NY,NX)=0.0_r8
    TQRF2P(NY,NX)=0.0_r8
    TQRC0P(NY,NX)=0.0_r8
    TQRC1P(NY,NX)=0.0_r8
    TQRC2P(NY,NX)=0.0_r8
    TQRM1P(NY,NX)=0.0_r8

    trcsa_TQS(idsa_beg:idsa_end,NY,NX)=0.0_r8
!
!     INITIALIZE NET SOLUTE AND GAS FLUXES FROM SNOWPACK DRIFT
!
    DO  L=1,JS
      trcsa_TBLS(idsa_beg:idsa_end,L,NY,NX)=0.0_r8
    ENDDO
  ENDIF
!
!     INITIALIZE NET SEDIMENT FLUXES FROM EROSION
!
  IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
    TSEDER(NY,NX)=0.0_r8
    TSANER(NY,NX)=0.0_r8
    TSILER(NY,NX)=0.0_r8
    TCLAER(NY,NX)=0.0_r8
    TCECER(NY,NX)=0.0_r8
    TAECER(NY,NX)=0.0_r8
    TNH4ER(NY,NX)=0.0_r8
    TNH3ER(NY,NX)=0.0_r8
    TNHUER(NY,NX)=0.0_r8
    TNO3ER(NY,NX)=0.0_r8
    TNH4EB(NY,NX)=0.0_r8
    TNH3EB(NY,NX)=0.0_r8
    TNHUEB(NY,NX)=0.0_r8
    TNO3EB(NY,NX)=0.0_r8
    TN4ER(NY,NX)=0.0_r8
    TNBER(NY,NX)=0.0_r8
    THYER(NY,NX)=0.0_r8
    TALER(NY,NX)=0.0_r8
    TFEER(NY,NX)=0.0_r8
    TCAER(NY,NX)=0.0_r8
    TMGER(NY,NX)=0.0_r8
    TNAER(NY,NX)=0.0_r8
    TKAER(NY,NX)=0.0_r8
    THCER(NY,NX)=0.0_r8
    TAL2ER(NY,NX)=0.0_r8
    TFE2ER(NY,NX)=0.0_r8
    TOH0ER(NY,NX)=0.0_r8
    TOH1ER(NY,NX)=0.0_r8
    TOH2ER(NY,NX)=0.0_r8
    TH1PER(NY,NX)=0.0_r8
    TH2PER(NY,NX)=0.0_r8
    TOH0EB(NY,NX)=0.0_r8
    TOH1EB(NY,NX)=0.0_r8
    TOH2EB(NY,NX)=0.0_r8
    TH1PEB(NY,NX)=0.0_r8
    TH2PEB(NY,NX)=0.0_r8
    TALOER(NY,NX)=0.0_r8
    TFEOER(NY,NX)=0.0_r8
    TCACER(NY,NX)=0.0_r8
    TCASER(NY,NX)=0.0_r8
    TALPER(NY,NX)=0.0_r8
    TFEPER(NY,NX)=0.0_r8
    TCPDER(NY,NX)=0.0_r8
    TCPHER(NY,NX)=0.0_r8
    TCPMER(NY,NX)=0.0_r8
    TALPEB(NY,NX)=0.0_r8
    TFEPEB(NY,NX)=0.0_r8
    TCPDEB(NY,NX)=0.0_r8
    TCPHEB(NY,NX)=0.0_r8
    TCPMEB(NY,NX)=0.0_r8

    TOMCER(1:nlbiomcp,1:NMICBSO,1:jcplx,NY,NX)=0.0_r8
    TOMNER(1:nlbiomcp,1:NMICBSO,1:jcplx,NY,NX)=0.0_r8
    TOMPER(1:nlbiomcp,1:NMICBSO,1:jcplx,NY,NX)=0.0_r8

    TOMCERff(1:nlbiomcp,1:NMICBSA,NY,NX)=0.0_r8
    TOMNERff(1:nlbiomcp,1:NMICBSA,NY,NX)=0.0_r8
    TOMPERff(1:nlbiomcp,1:NMICBSA,NY,NX)=0.0_r8

    TORCER(1:ndbiomcp,1:jcplx,NY,NX)=0.0_r8
    TORNER(1:ndbiomcp,1:jcplx,NY,NX)=0.0_r8
    TORPER(1:ndbiomcp,1:jcplx,NY,NX)=0.0_r8
    TOHCER(1:jcplx,NY,NX)=0.0_r8
    TOHNER(1:jcplx,NY,NX)=0.0_r8
    TOHPER(1:jcplx,NY,NX)=0.0_r8
    TOHAER(1:jcplx,NY,NX)=0.0_r8
    TOSCER(1:jsken,1:jcplx,NY,NX)=0.0_r8
    TOSAER(1:jsken,1:jcplx,NY,NX)=0.0_r8
    TOSNER(1:jsken,1:jcplx,NY,NX)=0.0_r8
    TOSPER(1:jsken,1:jcplx,NY,NX)=0.0_r8
  ENDIF
!
!     INITIALIZE NET SNOWPACK FLUXES WITHIN SNOWPACK
!
  DO  L=1,JS
    TFLWS(L,NY,NX)=0.0_r8
    TFLWW(L,NY,NX)=0.0_r8
    TFLWI(L,NY,NX)=0.0_r8
    THFLWW(L,NY,NX)=0.0_r8
  ENDDO
  end subroutine ZeroFluxArrays
!------------------------------------------------------------------------------------------

  subroutine ZeroFluxAccumulators(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,L
!     begin_execution
!
!     INITIALIZE WATER AND HEAT NET FLUX ACCUMULATORS WITHIN SOIL
!
  DO  L=NU(NY,NX),NL(NY,NX)
    TFLW(L,NY,NX)=0.0_r8
    TFLWX(L,NY,NX)=0.0_r8
    TFLWH(L,NY,NX)=0.0_r8
    THFLW(L,NY,NX)=0.0_r8
    TTHAW(L,NY,NX)=0.0_r8
    TTHAWH(L,NY,NX)=0.0_r8
    THTHAW(L,NY,NX)=0.0_r8
!
!     INITIALIZE GAS AND SOLUTE NET FLUX ACCUMULATORS WITHIN SOIL
!
    DO 8595 K=1,jcplx
      TOCFLS(K,L,NY,NX)=0.0_r8
      TONFLS(K,L,NY,NX)=0.0_r8
      TOPFLS(K,L,NY,NX)=0.0_r8
      TOAFLS(K,L,NY,NX)=0.0_r8
      TOCFHS(K,L,NY,NX)=0.0_r8
      TONFHS(K,L,NY,NX)=0.0_r8
      TOPFHS(K,L,NY,NX)=0.0_r8
      TOAFHS(K,L,NY,NX)=0.0_r8
8595  CONTINUE
    TCOFLS(L,NY,NX)=0.0_r8
      TCHFLS(L,NY,NX)=0.0_r8
      TOXFLS(L,NY,NX)=0.0_r8
      TNGFLS(L,NY,NX)=0.0_r8
      TN2FLS(L,NY,NX)=0.0_r8
      THGFLS(L,NY,NX)=0.0_r8
      TN4FLS(L,NY,NX)=0.0_r8
      TN3FLS(L,NY,NX)=0.0_r8
      TNOFLS(L,NY,NX)=0.0_r8
      TNXFLS(L,NY,NX)=0.0_r8
      TP1FLS(L,NY,NX)=0.0_r8
      TPOFLS(L,NY,NX)=0.0_r8
      TN4FLB(L,NY,NX)=0.0_r8
      TN3FLB(L,NY,NX)=0.0_r8
      TNOFLB(L,NY,NX)=0.0_r8
      TNXFLB(L,NY,NX)=0.0_r8
      TH1BFB(L,NY,NX)=0.0_r8
      TH2BFB(L,NY,NX)=0.0_r8
      TCOFHS(L,NY,NX)=0.0_r8
      TCHFHS(L,NY,NX)=0.0_r8
      TOXFHS(L,NY,NX)=0.0_r8
      TNGFHS(L,NY,NX)=0.0_r8
      TN2FHS(L,NY,NX)=0.0_r8
      THGFHS(L,NY,NX)=0.0_r8
      TN4FHS(L,NY,NX)=0.0_r8
      TN3FHS(L,NY,NX)=0.0_r8
      TNOFHS(L,NY,NX)=0.0_r8
      TNXFHS(L,NY,NX)=0.0_r8
      TP1FHS(L,NY,NX)=0.0_r8
      TPOFHS(L,NY,NX)=0.0_r8
      TN4FHB(L,NY,NX)=0.0_r8
      TN3FHB(L,NY,NX)=0.0_r8
      TNOFHB(L,NY,NX)=0.0_r8
      TNXFHB(L,NY,NX)=0.0_r8
      TH1BHB(L,NY,NX)=0.0_r8
      TH2BHB(L,NY,NX)=0.0_r8
      TCOFLG(L,NY,NX)=0.0_r8
      TCHFLG(L,NY,NX)=0.0_r8
      TOXFLG(L,NY,NX)=0.0_r8
      TNGFLG(L,NY,NX)=0.0_r8
      TN2FLG(L,NY,NX)=0.0_r8
      TNHFLG(L,NY,NX)=0.0_r8
      THGFLG(L,NY,NX)=0.0_r8
    IF(ISALTG.NE.0)THEN
      TALFLS(L,NY,NX)=0.0_r8
      TFEFLS(L,NY,NX)=0.0_r8
      THYFLS(L,NY,NX)=0.0_r8
      TCAFLS(L,NY,NX)=0.0_r8
      TMGFLS(L,NY,NX)=0.0_r8
      TNAFLS(L,NY,NX)=0.0_r8
      TKAFLS(L,NY,NX)=0.0_r8
      TOHFLS(L,NY,NX)=0.0_r8
      TSOFLS(L,NY,NX)=0.0_r8
      TCLFLS(L,NY,NX)=0.0_r8
      TC3FLS(L,NY,NX)=0.0_r8
      THCFLS(L,NY,NX)=0.0_r8
      TAL1FS(L,NY,NX)=0.0_r8
      TAL2FS(L,NY,NX)=0.0_r8
      TAL3FS(L,NY,NX)=0.0_r8
      TAL4FS(L,NY,NX)=0.0_r8
      TALSFS(L,NY,NX)=0.0_r8
      TFE1FS(L,NY,NX)=0.0_r8
      TFE2FS(L,NY,NX)=0.0_r8
      TFE3FS(L,NY,NX)=0.0_r8
      TFE4FS(L,NY,NX)=0.0_r8
      TFESFS(L,NY,NX)=0.0_r8
      TCAOFS(L,NY,NX)=0.0_r8
      TCACFS(L,NY,NX)=0.0_r8
      TCAHFS(L,NY,NX)=0.0_r8
      TCASFS(L,NY,NX)=0.0_r8
      TMGOFS(L,NY,NX)=0.0_r8
      TMGCFS(L,NY,NX)=0.0_r8
      TMGHFS(L,NY,NX)=0.0_r8
      TMGSFS(L,NY,NX)=0.0_r8
      TNACFS(L,NY,NX)=0.0_r8
      TNASFS(L,NY,NX)=0.0_r8
      TKASFS(L,NY,NX)=0.0_r8
      TH0PFS(L,NY,NX)=0.0_r8
      TH3PFS(L,NY,NX)=0.0_r8
      TF1PFS(L,NY,NX)=0.0_r8
      TF2PFS(L,NY,NX)=0.0_r8
      TC0PFS(L,NY,NX)=0.0_r8
      TC1PFS(L,NY,NX)=0.0_r8
      TC2PFS(L,NY,NX)=0.0_r8
      TM1PFS(L,NY,NX)=0.0_r8
      TH0BFB(L,NY,NX)=0.0_r8
      TH3BFB(L,NY,NX)=0.0_r8
      TF1BFB(L,NY,NX)=0.0_r8
      TF2BFB(L,NY,NX)=0.0_r8
      TC0BFB(L,NY,NX)=0.0_r8
      TC1BFB(L,NY,NX)=0.0_r8
      TC2BFB(L,NY,NX)=0.0_r8
      TM1BFB(L,NY,NX)=0.0_r8
      TALFHS(L,NY,NX)=0.0_r8
      TFEFHS(L,NY,NX)=0.0_r8
      THYFHS(L,NY,NX)=0.0_r8
      TCAFHS(L,NY,NX)=0.0_r8
      TMGFHS(L,NY,NX)=0.0_r8
      TNAFHS(L,NY,NX)=0.0_r8
      TKAFHS(L,NY,NX)=0.0_r8
      TOHFHS(L,NY,NX)=0.0_r8
      TSOFHS(L,NY,NX)=0.0_r8
      TCLFHS(L,NY,NX)=0.0_r8
      TC3FHS(L,NY,NX)=0.0_r8
      THCFHS(L,NY,NX)=0.0_r8
      TAL1HS(L,NY,NX)=0.0_r8
      TAL2HS(L,NY,NX)=0.0_r8
      TAL3HS(L,NY,NX)=0.0_r8
      TAL4HS(L,NY,NX)=0.0_r8
      TALSHS(L,NY,NX)=0.0_r8
      TFE1HS(L,NY,NX)=0.0_r8
      TFE2HS(L,NY,NX)=0.0_r8
      TFE3HS(L,NY,NX)=0.0_r8
      TFE4HS(L,NY,NX)=0.0_r8
      TFESHS(L,NY,NX)=0.0_r8
      TCAOHS(L,NY,NX)=0.0_r8
      TCACHS(L,NY,NX)=0.0_r8
      TCAHHS(L,NY,NX)=0.0_r8
      TCASHS(L,NY,NX)=0.0_r8
      TMGOHS(L,NY,NX)=0.0_r8
      TMGCHS(L,NY,NX)=0.0_r8
      TMGHHS(L,NY,NX)=0.0_r8
      TMGSHS(L,NY,NX)=0.0_r8
      TNACHS(L,NY,NX)=0.0_r8
      TNASHS(L,NY,NX)=0.0_r8
      TKASHS(L,NY,NX)=0.0_r8
      TH0PHS(L,NY,NX)=0.0_r8
      TH3PHS(L,NY,NX)=0.0_r8
      TF1PHS(L,NY,NX)=0.0_r8
      TF2PHS(L,NY,NX)=0.0_r8
      TC0PHS(L,NY,NX)=0.0_r8
      TC1PHS(L,NY,NX)=0.0_r8
      TC2PHS(L,NY,NX)=0.0_r8
      TM1PHS(L,NY,NX)=0.0_r8
      TH0BHB(L,NY,NX)=0.0_r8
      TH3BHB(L,NY,NX)=0.0_r8
      TF1BHB(L,NY,NX)=0.0_r8
      TF2BHB(L,NY,NX)=0.0_r8
      TC0BHB(L,NY,NX)=0.0_r8
      TC1BHB(L,NY,NX)=0.0_r8
      TC2BHB(L,NY,NX)=0.0_r8
      TM1BHB(L,NY,NX)=0.0_r8
    ENDIF
  ENDDO
  end subroutine ZeroFluxAccumulators

!------------------------------------------------------------------------------------------

  subroutine FluxesFromRunoff(N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B

  integer :: NN,K
  !     begin_execution
  !     NET WATER, SNOW AND HEAT FLUXES FROM RUNOFF
  !
  !     TQR,TQS,TQW,TQI=net water and snowpack snow,water,ice runoff
  !     THQR,THQS=net convective heat from surface water and snow runoff
  !     QR,HQR=runoff, convective heat from runoff from watsub.f
  !     QS,QW,QI=snow,water,ice transfer from watsub.f
  !     HQS=convective heat transfer from snow,water,ice transfer from watsub.f
  !
  D1202: DO NN=1,2
    TQR(N2,N1)=TQR(N2,N1)+QR(N,NN,N2,N1)
    THQR(N2,N1)=THQR(N2,N1)+HQR(N,NN,N2,N1)
    D8590: DO K=1,micpar%n_litrsfk
      TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)+XOCQRS(K,N,NN,N2,N1)
      TONQRS(K,N2,N1)=TONQRS(K,N2,N1)+XONQRS(K,N,NN,N2,N1)
      TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)+XOPQRS(K,N,NN,N2,N1)
      TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)+XOAQRS(K,N,NN,N2,N1)
    ENDDO D8590
    TCOQRS(N2,N1)=TCOQRS(N2,N1)+XCOQRS(N,NN,N2,N1)
    TCHQRS(N2,N1)=TCHQRS(N2,N1)+XCHQRS(N,NN,N2,N1)
    TOXQRS(N2,N1)=TOXQRS(N2,N1)+XOXQRS(N,NN,N2,N1)
    TNGQRS(N2,N1)=TNGQRS(N2,N1)+XNGQRS(N,NN,N2,N1)
    TN2QRS(N2,N1)=TN2QRS(N2,N1)+XN2QRS(N,NN,N2,N1)
    THGQRS(N2,N1)=THGQRS(N2,N1)+XHGQRS(N,NN,N2,N1)
    TN4QRS(N2,N1)=TN4QRS(N2,N1)+XN4QRW(N,NN,N2,N1)
    TN3QRS(N2,N1)=TN3QRS(N2,N1)+XN3QRW(N,NN,N2,N1)
    TNOQRS(N2,N1)=TNOQRS(N2,N1)+XNOQRW(N,NN,N2,N1)
    TNXQRS(N2,N1)=TNXQRS(N2,N1)+XNXQRS(N,NN,N2,N1)
    TP1QRS(N2,N1)=TP1QRS(N2,N1)+XP1QRW(N,NN,N2,N1)
    TPOQRS(N2,N1)=TPOQRS(N2,N1)+XP4QRW(N,NN,N2,N1)
    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
      TQR(N2,N1)=TQR(N2,N1)-QR(N,NN,N5,N4)
      THQR(N2,N1)=THQR(N2,N1)-HQR(N,NN,N5,N4)
      D8591: DO K=1,micpar%n_litrsfk
        TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)-XOCQRS(K,N,NN,N5,N4)
        TONQRS(K,N2,N1)=TONQRS(K,N2,N1)-XONQRS(K,N,NN,N5,N4)
        TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)-XOPQRS(K,N,NN,N5,N4)
        TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)-XOAQRS(K,N,NN,N5,N4)
      ENDDO D8591
      TCOQRS(N2,N1)=TCOQRS(N2,N1)-XCOQRS(N,NN,N5,N4)
      TCHQRS(N2,N1)=TCHQRS(N2,N1)-XCHQRS(N,NN,N5,N4)
      TOXQRS(N2,N1)=TOXQRS(N2,N1)-XOXQRS(N,NN,N5,N4)
      TNGQRS(N2,N1)=TNGQRS(N2,N1)-XNGQRS(N,NN,N5,N4)
      TN2QRS(N2,N1)=TN2QRS(N2,N1)-XN2QRS(N,NN,N5,N4)
      THGQRS(N2,N1)=THGQRS(N2,N1)-XHGQRS(N,NN,N5,N4)
      TN4QRS(N2,N1)=TN4QRS(N2,N1)-XN4QRW(N,NN,N5,N4)
      TN3QRS(N2,N1)=TN3QRS(N2,N1)-XN3QRW(N,NN,N5,N4)
      TNOQRS(N2,N1)=TNOQRS(N2,N1)-XNOQRW(N,NN,N5,N4)
      TNXQRS(N2,N1)=TNXQRS(N2,N1)-XNXQRS(N,NN,N5,N4)
      TP1QRS(N2,N1)=TP1QRS(N2,N1)-XP1QRW(N,NN,N5,N4)
      TPOQRS(N2,N1)=TPOQRS(N2,N1)-XP4QRW(N,NN,N5,N4)
    ENDIF
    !     WRITE(*,6631)'TQR',I,J,N1,N2,N4,N5,N,NN
    !    2,IFLBH(N,NN,N5,N4)
    !    2,TQR(N2,N1),QR(N,NN,N2,N1),QR(N,NN,N5,N4)
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQR(N2,N1)=TQR(N2,N1)-QR(N,NN,N5B,N4B)
      THQR(N2,N1)=THQR(N2,N1)-HQR(N,NN,N5B,N4B)
      D8592: DO K=1,micpar%n_litrsfk
        TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)-XOCQRS(K,N,NN,N5B,N4B)
        TONQRS(K,N2,N1)=TONQRS(K,N2,N1)-XONQRS(K,N,NN,N5B,N4B)
        TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)-XOPQRS(K,N,NN,N5B,N4B)
        TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)-XOAQRS(K,N,NN,N5B,N4B)
      ENDDO D8592
      TCOQRS(N2,N1)=TCOQRS(N2,N1)-XCOQRS(N,NN,N5B,N4B)
      TCHQRS(N2,N1)=TCHQRS(N2,N1)-XCHQRS(N,NN,N5B,N4B)
      TOXQRS(N2,N1)=TOXQRS(N2,N1)-XOXQRS(N,NN,N5B,N4B)
      TNGQRS(N2,N1)=TNGQRS(N2,N1)-XNGQRS(N,NN,N5B,N4B)
      TN2QRS(N2,N1)=TN2QRS(N2,N1)-XN2QRS(N,NN,N5B,N4B)
      THGQRS(N2,N1)=THGQRS(N2,N1)-XHGQRS(N,NN,N5B,N4B)
      TN4QRS(N2,N1)=TN4QRS(N2,N1)-XN4QRW(N,NN,N5B,N4B)
      TN3QRS(N2,N1)=TN3QRS(N2,N1)-XN3QRW(N,NN,N5B,N4B)
      TNOQRS(N2,N1)=TNOQRS(N2,N1)-XNOQRW(N,NN,N5B,N4B)
      TNXQRS(N2,N1)=TNXQRS(N2,N1)-XNXQRS(N,NN,N5B,N4B)
      TP1QRS(N2,N1)=TP1QRS(N2,N1)-XP1QRW(N,NN,N5B,N4B)
      TPOQRS(N2,N1)=TPOQRS(N2,N1)-XP4QRW(N,NN,N5B,N4B)
      !     WRITE(*,6631)'TQRB',I,J,N1,N2,N4B,N5B,N,NN
      !    2,IFLBH(N,NN,N5B,N4B)
      !    2,TQR(N2,N1),QR(N,NN,N5B,N4B)
!6631  FORMAT(A8,9I4,12E12.4)
    ENDIF
  ENDDO D1202
  TQS(N2,N1)=TQS(N2,N1)+QS(N,N2,N1)-QS(N,N5,N4)
  TQW(N2,N1)=TQW(N2,N1)+QW(N,N2,N1)-QW(N,N5,N4)
  TQI(N2,N1)=TQI(N2,N1)+QI(N,N2,N1)-QI(N,N5,N4)
  THQS(N2,N1)=THQS(N2,N1)+HQS(N,N2,N1)-HQS(N,N5,N4)
  !
  !     NET GAS AND SOLUTE FLUXES FROM RUNOFF AND SNOWPACK
  !
  !     T*QRS=net overland solute flux from runoff
  !     X*QRS=solute in runoff from trnsfr.f
  !     T*QSS=net overland solute flux from snowpack
  !     X*QSS=solute in snowpack flux from trnsfr.f
  !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
  !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
  !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
  !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
  !
  trcg_QSS(idg_CO2,N2,N1)=trcg_QSS(idg_CO2,N2,N1)+XCOQSS(N,N2,N1)-XCOQSS(N,N5,N4)
  trcg_QSS(idg_CH4,N2,N1)=trcg_QSS(idg_CH4,N2,N1)+XCHQSS(N,N2,N1)-XCHQSS(N,N5,N4)
  trcg_QSS(idg_O2,N2,N1)=trcg_QSS(idg_O2,N2,N1)+XOXQSS(N,N2,N1)-XOXQSS(N,N5,N4)
  trcg_QSS(idg_N2,N2,N1)=trcg_QSS(idg_N2,N2,N1)+XNGQSS(N,N2,N1)-XNGQSS(N,N5,N4)
  trcg_QSS(idg_N2O,N2,N1)=trcg_QSS(idg_N2O,N2,N1)+XN2QSS(N,N2,N1)-XN2QSS(N,N5,N4)
  trcn_QSS(ids_NH4,N2,N1)=trcn_QSS(ids_NH4,N2,N1)+XN4QSS(N,N2,N1)-XN4QSS(N,N5,N4)
  trcg_QSS(idg_NH3,N2,N1)=trcg_QSS(idg_NH3,N2,N1)+XN3QSS(N,N2,N1)-XN3QSS(N,N5,N4)
  trcn_QSS(ids_NO3,N2,N1)=trcn_QSS(ids_NO3,N2,N1)+XNOQSS(N,N2,N1)-XNOQSS(N,N5,N4)
  trcn_QSS(ids_H1PO4,N2,N1)=trcn_QSS(ids_H1PO4,N2,N1)+XP1QSS(N,N2,N1)-XP1QSS(N,N5,N4)
  trcn_QSS(ids_H2PO4,N2,N1)=trcn_QSS(ids_H2PO4,N2,N1)+XP4QSS(N,N2,N1)-XP4QSS(N,N5,N4)
  end subroutine FluxesFromRunoff
!------------------------------------------------------------------------------------------

  subroutine SaltThruFluxRunoffAndSnowpack(N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B

  integer :: NN
!     begin_execution
!
!     TQR*=net overland solute flux in runoff
!     XQR*=solute in runoff from trnsfrs.f
!     TQS*=net overland solute flux in snow drift
!     XQS*=solute in snow drift from trnsfrs.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  IF(ISALTG.NE.0)THEN
    D1203: DO NN=1,2
      TQRAL(N2,N1)=TQRAL(N2,N1)+XQRAL(N,NN,N2,N1)
      TQRFE(N2,N1)=TQRFE(N2,N1)+XQRFE(N,NN,N2,N1)
      TQRHY(N2,N1)=TQRHY(N2,N1)+XQRHY(N,NN,N2,N1)
      TQRCA(N2,N1)=TQRCA(N2,N1)+XQRCA(N,NN,N2,N1)
      TQRMG(N2,N1)=TQRMG(N2,N1)+XQRMG(N,NN,N2,N1)
      TQRNA(N2,N1)=TQRNA(N2,N1)+XQRNA(N,NN,N2,N1)
      TQRKA(N2,N1)=TQRKA(N2,N1)+XQRKA(N,NN,N2,N1)
      TQROH(N2,N1)=TQROH(N2,N1)+XQROH(N,NN,N2,N1)
      TQRSO(N2,N1)=TQRSO(N2,N1)+XQRSO(N,NN,N2,N1)
      TQRCL(N2,N1)=TQRCL(N2,N1)+XQRCL(N,NN,N2,N1)
      TQRC3(N2,N1)=TQRC3(N2,N1)+XQRC3(N,NN,N2,N1)
      TQRHC(N2,N1)=TQRHC(N2,N1)+XQRHC(N,NN,N2,N1)
      TQRAL1(N2,N1)=TQRAL1(N2,N1)+XQRAL1(N,NN,N2,N1)
      TQRAL2(N2,N1)=TQRAL2(N2,N1)+XQRAL2(N,NN,N2,N1)
      TQRAL3(N2,N1)=TQRAL3(N2,N1)+XQRAL3(N,NN,N2,N1)
      TQRAL4(N2,N1)=TQRAL4(N2,N1)+XQRAL4(N,NN,N2,N1)
      TQRALS(N2,N1)=TQRALS(N2,N1)+XQRALS(N,NN,N2,N1)
      TQRFE1(N2,N1)=TQRFE1(N2,N1)+XQRFE1(N,NN,N2,N1)
      TQRFE2(N2,N1)=TQRFE2(N2,N1)+XQRFE2(N,NN,N2,N1)
      TQRFE3(N2,N1)=TQRFE3(N2,N1)+XQRFE3(N,NN,N2,N1)
      TQRFE4(N2,N1)=TQRFE4(N2,N1)+XQRFE4(N,NN,N2,N1)
      TQRFES(N2,N1)=TQRFES(N2,N1)+XQRFES(N,NN,N2,N1)
      TQRCAO(N2,N1)=TQRCAO(N2,N1)+XQRCAO(N,NN,N2,N1)
      TQRCAC(N2,N1)=TQRCAC(N2,N1)+XQRCAC(N,NN,N2,N1)
      TQRCAH(N2,N1)=TQRCAH(N2,N1)+XQRCAH(N,NN,N2,N1)
      TQRCAS(N2,N1)=TQRCAS(N2,N1)+XQRCAS(N,NN,N2,N1)
      TQRMGO(N2,N1)=TQRMGO(N2,N1)+XQRMGO(N,NN,N2,N1)
      TQRMGC(N2,N1)=TQRMGC(N2,N1)+XQRMGC(N,NN,N2,N1)
      TQRMGH(N2,N1)=TQRMGH(N2,N1)+XQRMGH(N,NN,N2,N1)
      TQRMGS(N2,N1)=TQRMGS(N2,N1)+XQRMGS(N,NN,N2,N1)
      TQRNAC(N2,N1)=TQRNAC(N2,N1)+XQRNAC(N,NN,N2,N1)
      TQRNAS(N2,N1)=TQRNAS(N2,N1)+XQRNAS(N,NN,N2,N1)
      TQRKAS(N2,N1)=TQRKAS(N2,N1)+XQRKAS(N,NN,N2,N1)
      TQRH0P(N2,N1)=TQRH0P(N2,N1)+XQRH0P(N,NN,N2,N1)
      TQRH3P(N2,N1)=TQRH3P(N2,N1)+XQRH3P(N,NN,N2,N1)
      TQRF1P(N2,N1)=TQRF1P(N2,N1)+XQRF1P(N,NN,N2,N1)
      TQRF2P(N2,N1)=TQRF2P(N2,N1)+XQRF2P(N,NN,N2,N1)
      TQRC0P(N2,N1)=TQRC0P(N2,N1)+XQRC0P(N,NN,N2,N1)
      TQRC1P(N2,N1)=TQRC1P(N2,N1)+XQRC1P(N,NN,N2,N1)
      TQRC2P(N2,N1)=TQRC2P(N2,N1)+XQRC2P(N,NN,N2,N1)
      TQRM1P(N2,N1)=TQRM1P(N2,N1)+XQRM1P(N,NN,N2,N1)
      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
        TQRAL(N2,N1)=TQRAL(N2,N1)-XQRAL(N,NN,N5,N4)
        TQRFE(N2,N1)=TQRFE(N2,N1)-XQRFE(N,NN,N5,N4)
        TQRHY(N2,N1)=TQRHY(N2,N1)-XQRHY(N,NN,N5,N4)
        TQRCA(N2,N1)=TQRCA(N2,N1)-XQRCA(N,NN,N5,N4)
        TQRMG(N2,N1)=TQRMG(N2,N1)-XQRMG(N,NN,N5,N4)
        TQRNA(N2,N1)=TQRNA(N2,N1)-XQRNA(N,NN,N5,N4)
        TQRKA(N2,N1)=TQRKA(N2,N1)-XQRKA(N,NN,N5,N4)
        TQROH(N2,N1)=TQROH(N2,N1)-XQROH(N,NN,N5,N4)
        TQRSO(N2,N1)=TQRSO(N2,N1)-XQRSO(N,NN,N5,N4)
        TQRCL(N2,N1)=TQRCL(N2,N1)-XQRCL(N,NN,N5,N4)
        TQRC3(N2,N1)=TQRC3(N2,N1)-XQRC3(N,NN,N5,N4)
        TQRHC(N2,N1)=TQRHC(N2,N1)-XQRHC(N,NN,N5,N4)
        TQRAL1(N2,N1)=TQRAL1(N2,N1)+XQRAL1(N,NN,N2,N1)
        TQRAL2(N2,N1)=TQRAL2(N2,N1)+XQRAL2(N,NN,N2,N1)
        TQRAL3(N2,N1)=TQRAL3(N2,N1)+XQRAL3(N,NN,N2,N1)
        TQRAL4(N2,N1)=TQRAL4(N2,N1)+XQRAL4(N,NN,N2,N1)
        TQRALS(N2,N1)=TQRALS(N2,N1)+XQRALS(N,NN,N2,N1)
        TQRFE1(N2,N1)=TQRFE1(N2,N1)+XQRFE1(N,NN,N2,N1)
        TQRFE2(N2,N1)=TQRFE2(N2,N1)+XQRFE2(N,NN,N2,N1)
        TQRFE3(N2,N1)=TQRFE3(N2,N1)+XQRFE3(N,NN,N2,N1)
        TQRFE4(N2,N1)=TQRFE4(N2,N1)+XQRFE4(N,NN,N2,N1)
        TQRFES(N2,N1)=TQRFES(N2,N1)+XQRFES(N,NN,N2,N1)
        TQRCAO(N2,N1)=TQRCAO(N2,N1)+XQRCAO(N,NN,N2,N1)
        TQRCAC(N2,N1)=TQRCAC(N2,N1)+XQRCAC(N,NN,N2,N1)
        TQRCAH(N2,N1)=TQRCAH(N2,N1)+XQRCAH(N,NN,N2,N1)
        TQRCAS(N2,N1)=TQRCAS(N2,N1)+XQRCAS(N,NN,N2,N1)
        TQRMGO(N2,N1)=TQRMGO(N2,N1)+XQRMGO(N,NN,N2,N1)
        TQRMGC(N2,N1)=TQRMGC(N2,N1)+XQRMGC(N,NN,N2,N1)
        TQRMGH(N2,N1)=TQRMGH(N2,N1)+XQRMGH(N,NN,N2,N1)
        TQRMGS(N2,N1)=TQRMGS(N2,N1)+XQRMGS(N,NN,N2,N1)
        TQRNAC(N2,N1)=TQRNAC(N2,N1)+XQRNAC(N,NN,N2,N1)
        TQRNAS(N2,N1)=TQRNAS(N2,N1)+XQRNAS(N,NN,N2,N1)
        TQRKAS(N2,N1)=TQRKAS(N2,N1)+XQRKAS(N,NN,N2,N1)
        TQRH0P(N2,N1)=TQRH0P(N2,N1)+XQRH0P(N,NN,N2,N1)
        TQRH3P(N2,N1)=TQRH3P(N2,N1)+XQRH3P(N,NN,N2,N1)
        TQRF1P(N2,N1)=TQRF1P(N2,N1)+XQRF1P(N,NN,N2,N1)
        TQRF2P(N2,N1)=TQRF2P(N2,N1)+XQRF2P(N,NN,N2,N1)
        TQRC0P(N2,N1)=TQRC0P(N2,N1)+XQRC0P(N,NN,N2,N1)
        TQRC1P(N2,N1)=TQRC1P(N2,N1)+XQRC1P(N,NN,N2,N1)
        TQRC2P(N2,N1)=TQRC2P(N2,N1)+XQRC2P(N,NN,N2,N1)
        TQRM1P(N2,N1)=TQRM1P(N2,N1)+XQRM1P(N,NN,N2,N1)
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        TQRAL(N2,N1)=TQRAL(N2,N1)-XQRAL(N,NN,N5B,N4B)
        TQRFE(N2,N1)=TQRFE(N2,N1)-XQRFE(N,NN,N5B,N4B)
        TQRHY(N2,N1)=TQRHY(N2,N1)-XQRHY(N,NN,N5B,N4B)
        TQRCA(N2,N1)=TQRCA(N2,N1)-XQRCA(N,NN,N5B,N4B)
        TQRMG(N2,N1)=TQRMG(N2,N1)-XQRMG(N,NN,N5B,N4B)
        TQRNA(N2,N1)=TQRNA(N2,N1)-XQRNA(N,NN,N5B,N4B)
        TQRKA(N2,N1)=TQRKA(N2,N1)-XQRKA(N,NN,N5B,N4B)
        TQROH(N2,N1)=TQROH(N2,N1)-XQROH(N,NN,N5B,N4B)
        TQRSO(N2,N1)=TQRSO(N2,N1)-XQRSO(N,NN,N5B,N4B)
        TQRCL(N2,N1)=TQRCL(N2,N1)-XQRCL(N,NN,N5B,N4B)
        TQRC3(N2,N1)=TQRC3(N2,N1)-XQRC3(N,NN,N5B,N4B)
        TQRHC(N2,N1)=TQRHC(N2,N1)-XQRHC(N,NN,N5B,N4B)
        TQRAL1(N2,N1)=TQRAL1(N2,N1)-XQRAL1(N,NN,N5B,N4B)
        TQRAL2(N2,N1)=TQRAL2(N2,N1)-XQRAL2(N,NN,N5B,N4B)
        TQRAL3(N2,N1)=TQRAL3(N2,N1)-XQRAL3(N,NN,N5B,N4B)
        TQRAL4(N2,N1)=TQRAL4(N2,N1)-XQRAL4(N,NN,N5B,N4B)
        TQRALS(N2,N1)=TQRALS(N2,N1)-XQRALS(N,NN,N5B,N4B)
        TQRFE1(N2,N1)=TQRFE1(N2,N1)-XQRFE1(N,NN,N5B,N4B)
        TQRFE2(N2,N1)=TQRFE2(N2,N1)-XQRFE2(N,NN,N5B,N4B)
        TQRFE3(N2,N1)=TQRFE3(N2,N1)-XQRFE3(N,NN,N5B,N4B)
        TQRFE4(N2,N1)=TQRFE4(N2,N1)-XQRFE4(N,NN,N5B,N4B)
        TQRFES(N2,N1)=TQRFES(N2,N1)-XQRFES(N,NN,N5B,N4B)
        TQRCAO(N2,N1)=TQRCAO(N2,N1)-XQRCAO(N,NN,N5B,N4B)
        TQRCAC(N2,N1)=TQRCAC(N2,N1)-XQRCAC(N,NN,N5B,N4B)
        TQRCAH(N2,N1)=TQRCAH(N2,N1)-XQRCAH(N,NN,N5B,N4B)
        TQRCAS(N2,N1)=TQRCAS(N2,N1)-XQRCAS(N,NN,N5B,N4B)
        TQRMGO(N2,N1)=TQRMGO(N2,N1)-XQRMGO(N,NN,N5B,N4B)
        TQRMGC(N2,N1)=TQRMGC(N2,N1)-XQRMGC(N,NN,N5B,N4B)
        TQRMGH(N2,N1)=TQRMGH(N2,N1)-XQRMGH(N,NN,N5B,N4B)
        TQRMGS(N2,N1)=TQRMGS(N2,N1)-XQRMGS(N,NN,N5B,N4B)
        TQRNAC(N2,N1)=TQRNAC(N2,N1)-XQRNAC(N,NN,N5B,N4B)
        TQRNAS(N2,N1)=TQRNAS(N2,N1)-XQRNAS(N,NN,N5B,N4B)
        TQRKAS(N2,N1)=TQRKAS(N2,N1)-XQRKAS(N,NN,N5B,N4B)
        TQRH0P(N2,N1)=TQRH0P(N2,N1)-XQRH0P(N,NN,N5B,N4B)
        TQRH3P(N2,N1)=TQRH3P(N2,N1)-XQRH3P(N,NN,N5B,N4B)
        TQRF1P(N2,N1)=TQRF1P(N2,N1)-XQRF1P(N,NN,N5B,N4B)
        TQRF2P(N2,N1)=TQRF2P(N2,N1)-XQRF2P(N,NN,N5B,N4B)
        TQRC0P(N2,N1)=TQRC0P(N2,N1)-XQRC0P(N,NN,N5B,N4B)
        TQRC1P(N2,N1)=TQRC1P(N2,N1)-XQRC1P(N,NN,N5B,N4B)
        TQRC2P(N2,N1)=TQRC2P(N2,N1)-XQRC2P(N,NN,N5B,N4B)
        TQRM1P(N2,N1)=TQRM1P(N2,N1)-XQRM1P(N,NN,N5B,N4B)
      ENDIF
    ENDDO D1203
    trcsa_TQS(idsa_Al,N2,N1)=trcsa_TQS(idsa_Al,N2,N1)+XQSAL(N,N2,N1)-XQSAL(N,N5,N4)
    trcsa_TQS(idsa_Fe,N2,N1)=trcsa_TQS(idsa_Fe,N2,N1)+XQSFE(N,N2,N1)-XQSFE(N,N5,N4)
    trcsa_TQS(idsa_Hp,N2,N1)=trcsa_TQS(idsa_Hp,N2,N1)+XQSHY(N,N2,N1)-XQSHY(N,N5,N4)
    trcsa_TQS(idsa_Ca,N2,N1)=trcsa_TQS(idsa_Ca,N2,N1)+XQSCA(N,N2,N1)-XQSCA(N,N5,N4)
    trcsa_TQS(idsa_Mg,N2,N1)=trcsa_TQS(idsa_Mg,N2,N1)+XQSMG(N,N2,N1)-XQSMG(N,N5,N4)
    trcsa_TQS(idsa_Na,N2,N1)=trcsa_TQS(idsa_Na,N2,N1)+XQSNA(N,N2,N1)-XQSNA(N,N5,N4)
    trcsa_TQS(idsa_K,N2,N1)=trcsa_TQS(idsa_K,N2,N1)+XQSKA(N,N2,N1)-XQSKA(N,N5,N4)
    trcsa_TQS(idsa_OH,N2,N1)=trcsa_TQS(idsa_OH,N2,N1)+XQSOH(N,N2,N1)-XQSOH(N,N5,N4)
    trcsa_TQS(idsa_SO4,N2,N1)=trcsa_TQS(idsa_SO4,N2,N1)+XQSSO(N,N2,N1)-XQSSO(N,N5,N4)
    trcsa_TQS(idsa_Cl,N2,N1)=trcsa_TQS(idsa_Cl,N2,N1)+XQSCL(N,N2,N1)-XQSCL(N,N5,N4)
    trcsa_TQS(idsa_CO3,N2,N1)=trcsa_TQS(idsa_CO3,N2,N1)+XQSC3(N,N2,N1)-XQSC3(N,N5,N4)
    trcsa_TQS(idsa_HCO3,N2,N1)=trcsa_TQS(idsa_HCO3,N2,N1)+XQSHC(N,N2,N1)-XQSHC(N,N5,N4)
    trcsa_TQS(idsa_AlOH,N2,N1)=trcsa_TQS(idsa_AlOH,N2,N1)+XQSAL1(N,N2,N1)-XQSAL1(N,N5,N4)
    trcsa_TQS(idsa_AlOH2,N2,N1)=trcsa_TQS(idsa_AlOH2,N2,N1)+XQSAL2(N,N2,N1)-XQSAL2(N,N5,N4)
    trcsa_TQS(idsa_AlOH3,N2,N1)=trcsa_TQS(idsa_AlOH3,N2,N1)+XQSAL3(N,N2,N1)-XQSAL3(N,N5,N4)
    trcsa_TQS(idsa_AlOH4,N2,N1)=trcsa_TQS(idsa_AlOH4,N2,N1)+XQSAL4(N,N2,N1)-XQSAL4(N,N5,N4)
    trcsa_TQS(idsa_AlSO4,N2,N1)=trcsa_TQS(idsa_AlSO4,N2,N1)+XQSALS(N,N2,N1)-XQSALS(N,N5,N4)
    trcsa_TQS(idsa_FeOH,N2,N1)=trcsa_TQS(idsa_FeOH,N2,N1)+XQSFE1(N,N2,N1)-XQSFE1(N,N5,N4)
    trcsa_TQS(idsa_FeOH2,N2,N1)=trcsa_TQS(idsa_FeOH2,N2,N1)+XQSFE2(N,N2,N1)-XQSFE2(N,N5,N4)
    trcsa_TQS(idsa_FeOH3,N2,N1)=trcsa_TQS(idsa_FeOH3,N2,N1)+XQSFE3(N,N2,N1)-XQSFE3(N,N5,N4)
    trcsa_TQS(idsa_FeOH4,N2,N1)=trcsa_TQS(idsa_FeOH4,N2,N1)+XQSFE4(N,N2,N1)-XQSFE4(N,N5,N4)
    trcsa_TQS(idsa_FeSO4,N2,N1)=trcsa_TQS(idsa_FeSO4,N2,N1)+XQSFES(N,N2,N1)-XQSFES(N,N5,N4)
    trcsa_TQS(idsa_CaOH2,N2,N1)=trcsa_TQS(idsa_CaOH2,N2,N1)+XQSCAO(N,N2,N1)-XQSCAO(N,N5,N4)
    trcsa_TQS(idsa_CaCO3,N2,N1)=trcsa_TQS(idsa_CaCO3,N2,N1)+XQSCAC(N,N2,N1)-XQSCAC(N,N5,N4)
    trcsa_TQS(idsa_CaHCO3,N2,N1)=trcsa_TQS(idsa_CaHCO3,N2,N1)+XQSCAH(N,N2,N1)-XQSCAH(N,N5,N4)
    trcsa_TQS(idsa_CaSO4,N2,N1)=trcsa_TQS(idsa_CaSO4,N2,N1)+XQSCAS(N,N2,N1)-XQSCAS(N,N5,N4)
    trcsa_TQS(idsa_MgOH2,N2,N1)=trcsa_TQS(idsa_MgOH2,N2,N1)+XQSMGO(N,N2,N1)-XQSMGO(N,N5,N4)
    trcsa_TQS(idsa_MgCO3,N2,N1)=trcsa_TQS(idsa_MgCO3,N2,N1)+XQSMGC(N,N2,N1)-XQSMGC(N,N5,N4)
    trcsa_TQS(idsa_MgHCO3,N2,N1)=trcsa_TQS(idsa_MgHCO3,N2,N1)+XQSMGH(N,N2,N1)-XQSMGH(N,N5,N4)
    trcsa_TQS(idsa_MgSO4,N2,N1)=trcsa_TQS(idsa_MgSO4,N2,N1)+XQSMGS(N,N2,N1)-XQSMGS(N,N5,N4)
    trcsa_TQS(idsa_NaCO3,N2,N1)=trcsa_TQS(idsa_NaCO3,N2,N1)+XQSNAC(N,N2,N1)-XQSNAC(N,N5,N4)
    trcsa_TQS(idsa_NaSO4,N2,N1)=trcsa_TQS(idsa_NaSO4,N2,N1)+XQSNAS(N,N2,N1)-XQSNAS(N,N5,N4)
    trcsa_TQS(idsa_KSO4,N2,N1)=trcsa_TQS(idsa_KSO4,N2,N1)+XQSKAS(N,N2,N1)-XQSKAS(N,N5,N4)
    trcsa_TQS(idsa_H0PO4,N2,N1)=trcsa_TQS(idsa_H0PO4,N2,N1)+XQSH0P(N,N2,N1)-XQSH0P(N,N5,N4)
    trcsa_TQS(idsa_H3PO4,N2,N1)=trcsa_TQS(idsa_H3PO4,N2,N1)+XQSH3P(N,N2,N1)-XQSH3P(N,N5,N4)
    trcsa_TQS(idsa_FeHPO4,N2,N1)=trcsa_TQS(idsa_FeHPO4,N2,N1)+XQSF1P(N,N2,N1)-XQSF1P(N,N5,N4)
    trcsa_TQS(idsa_FeH2PO4,N2,N1)=trcsa_TQS(idsa_FeH2PO4,N2,N1)+XQSF2P(N,N2,N1)-XQSF2P(N,N5,N4)
    trcsa_TQS(idsa_CaPO4,N2,N1)=trcsa_TQS(idsa_CaPO4,N2,N1)+XQSC0P(N,N2,N1)-XQSC0P(N,N5,N4)
    trcsa_TQS(idsa_CaHPO4,N2,N1)=trcsa_TQS(idsa_CaHPO4,N2,N1)+XQSC1P(N,N2,N1)-XQSC1P(N,N5,N4)
    trcsa_TQS(idsa_CaH2PO4,N2,N1)=trcsa_TQS(idsa_CaH2PO4,N2,N1)+XQSC2P(N,N2,N1)-XQSC2P(N,N5,N4)
    trcsa_TQS(idsa_MgHPO4,N2,N1)=trcsa_TQS(idsa_MgHPO4,N2,N1)+XQSM1P(N,N2,N1)-XQSM1P(N,N5,N4)
  ENDIF
  end subroutine SaltThruFluxRunoffAndSnowpack
!------------------------------------------------------------------------------------------

  subroutine SaltFromRunoffSnowpack(N1,N2,NY,NX)
  implicit none
  integer, intent(in) :: N1,N2,NY,NX
  integer :: LS, LS2
!     begin_execution
!     NET WATER AND HEAT FLUXES THROUGH SNOWPACK
!
!     VHCPW,VHCPWX=current, minimum snowpack heat capacities
!     TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
!     THFLWW=convective heat fluxes of snow,water,ice in snowpack
!     XFLWS,XFLWW,XFLWI=snow,water,ice transfer from watsub.f
!     XHFLWW=convective heat flux from snow,water,ice transfer from watsub.f
!     FLSW,FLSWH,FLSWR=water flux from lowest snow layer to soil macropore,micropore,litter
!     HFLSW,HFLSWR=heat flux from lowest snow layer to soil,litter

  D1205: DO LS=1,JS
    IF(VHCPW(LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS.AND.VHCPW(LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
      TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1) &
      -XFLWS(LS2,N2,N1)
      TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1) &
      -XFLWW(LS2,N2,N1) &
      -FLSWR(LS,N2,N1)-FLSW(LS,N2,N1)-FLSWH(LS,N2,N1)
      TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1) &
      -XFLWI(LS2,N2,N1)
      THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1) &
      -XHFLWW(LS2,N2,N1) &
      -HFLSWR(LS,N2,N1)-HFLSW(LS,N2,N1)
!     IF(LS.EQ.5)THEN
!     WRITE(*,7754)'LS',I,J,N1,N2,LS,LS2,TFLWW(LS,N2,N1),XFLWW(LS,N2,N1)
!    2,XFLWW(LS2,N2,N1)
!7754  FORMAT(A8,6I4,100E14.6)
!     ENDIF
!
!     NET SOLUTE FLUXES THROUGH SNOWPACK
!
!     T*BLS=net solute flux in snowpack
!     X*BLS=solute flux in snowpack from trnsfr.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
      trcg_TBLS(idg_CO2,LS,N2,N1)=trcg_TBLS(idg_CO2,LS,N2,N1)+XCOBLS(LS,N2,N1) &
      -XCOBLS(LS2,N2,N1)
      trcg_TBLS(idg_CH4,LS,N2,N1)=trcg_TBLS(idg_CH4,LS,N2,N1)+XCHBLS(LS,N2,N1) &
      -XCHBLS(LS2,N2,N1)
      trcg_TBLS(idg_O2,LS,N2,N1)=trcg_TBLS(idg_O2,LS,N2,N1)+XOXBLS(LS,N2,N1) &
      -XOXBLS(LS2,N2,N1)
      trcg_TBLS(idg_N2,LS,N2,N1)=trcg_TBLS(idg_N2,LS,N2,N1)+XNGBLS(LS,N2,N1) &
      -XNGBLS(LS2,N2,N1)
      trcg_TBLS(idg_N2O,LS,N2,N1)=trcg_TBLS(idg_N2O,LS,N2,N1)+XN2BLS(LS,N2,N1) &
      -XN2BLS(LS2,N2,N1)
      trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1)+XN3BLW(LS,N2,N1) &
      -XN3BLW(LS2,N2,N1)

      trcn_TBLS(ids_NH4,LS,N2,N1)=trcn_TBLS(ids_NH4,LS,N2,N1)+XN4BLW(LS,N2,N1) &
      -XN4BLW(LS2,N2,N1)
      trcn_TBLS(ids_NO3,LS,N2,N1)=trcn_TBLS(ids_NO3,LS,N2,N1)+XNOBLW(LS,N2,N1) &
      -XNOBLW(LS2,N2,N1)
      trcn_TBLS(ids_H1PO4,LS,N2,N1)=trcn_TBLS(ids_H1PO4,LS,N2,N1)+XH1PBS(LS,N2,N1) &
      -XH1PBS(LS2,N2,N1)
      trcn_TBLS(ids_H2PO4,LS,N2,N1)=trcn_TBLS(ids_H2PO4,LS,N2,N1)+XH2PBS(LS,N2,N1) &
      -XH2PBS(LS2,N2,N1)
!
!     NET SALT FLUXES THROUGH SNOWPACK
!
!     T*BLS=net solute flux in snowpack
!     X*BLS=solute flux in snowpack from trnsfrs.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
      IF(ISALTG.NE.0)THEN
      trcsa_TBLS(idsa_Al,LS,N2,N1)=trcsa_TBLS(idsa_Al,LS,N2,N1)+XALBLS(LS,N2,N1) &
      -XALBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_Fe,LS,N2,N1)=trcsa_TBLS(idsa_Fe,LS,N2,N1)+XFEBLS(LS,N2,N1) &
      -XFEBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_Hp,LS,N2,N1)=trcsa_TBLS(idsa_Hp,LS,N2,N1)+XHYBLS(LS,N2,N1) &
      -XHYBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_Ca,LS,N2,N1)=trcsa_TBLS(idsa_Ca,LS,N2,N1)+XCABLS(LS,N2,N1) &
      -XCABLS(LS2,N2,N1)
      trcsa_TBLS(idsa_Mg,LS,N2,N1)=trcsa_TBLS(idsa_Mg,LS,N2,N1)+XMGBLS(LS,N2,N1) &
      -XMGBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_Na,LS,N2,N1)=trcsa_TBLS(idsa_Na,LS,N2,N1)+XNABLS(LS,N2,N1) &
      -XNABLS(LS2,N2,N1)
      trcsa_TBLS(idsa_K,LS,N2,N1)=trcsa_TBLS(idsa_K,LS,N2,N1)+XKABLS(LS,N2,N1) &
      -XKABLS(LS2,N2,N1)
      trcsa_TBLS(idsa_OH,LS,N2,N1)=trcsa_TBLS(idsa_OH,LS,N2,N1)+XOHBLS(LS,N2,N1) &
      -XOHBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_SO4,LS,N2,N1)=trcsa_TBLS(idsa_SO4,LS,N2,N1)+XSOBLS(LS,N2,N1) &
      -XSOBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_Cl,LS,N2,N1)=trcsa_TBLS(idsa_Cl,LS,N2,N1)+XCLBLS(LS,N2,N1) &
      -XCLBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_CO3,LS,N2,N1)=trcsa_TBLS(idsa_CO3,LS,N2,N1)+XC3BLS(LS,N2,N1) &
      -XC3BLS(LS2,N2,N1)
      trcsa_TBLS(idsa_HCO3,LS,N2,N1)=trcsa_TBLS(idsa_HCO3,LS,N2,N1)+XHCBLS(LS,N2,N1) &
      -XHCBLS(LS2,N2,N1)
      trcsa_TBLS(idsa_AlOH,LS,N2,N1)=trcsa_TBLS(idsa_AlOH,LS,N2,N1)+XAL1BS(LS,N2,N1) &
      -XAL1BS(LS2,N2,N1)
      trcsa_TBLS(idsa_AlOH2,LS,N2,N1)=trcsa_TBLS(idsa_AlOH2,LS,N2,N1)+XAL2BS(LS,N2,N1) &
      -XAL2BS(LS2,N2,N1)
      trcsa_TBLS(idsa_AlOH3,LS,N2,N1)=trcsa_TBLS(idsa_AlOH3,LS,N2,N1)+XAL3BS(LS,N2,N1) &
      -XAL3BS(LS2,N2,N1)
      trcsa_TBLS(idsa_AlOH4,LS,N2,N1)=trcsa_TBLS(idsa_AlOH4,LS,N2,N1)+XAL4BS(LS,N2,N1) &
      -XAL4BS(LS2,N2,N1)
      trcsa_TBLS(idsa_AlSO4,LS,N2,N1)=trcsa_TBLS(idsa_AlSO4,LS,N2,N1)+XALSBS(LS,N2,N1) &
      -XALSBS(LS2,N2,N1)
      trcsa_TBLS(idsa_FeOH,LS,N2,N1)=trcsa_TBLS(idsa_FeOH,LS,N2,N1)+XFE1BS(LS,N2,N1) &
      -XFE1BS(LS2,N2,N1)
      trcsa_TBLS(idsa_FeOH2,LS,N2,N1)=trcsa_TBLS(idsa_FeOH2,LS,N2,N1)+XFE2BS(LS,N2,N1) &
      -XFE2BS(LS2,N2,N1)
      trcsa_TBLS(idsa_FeOH3,LS,N2,N1)=trcsa_TBLS(idsa_FeOH3,LS,N2,N1)+XFE3BS(LS,N2,N1) &
      -XFE3BS(LS2,N2,N1)
      trcsa_TBLS(idsa_FeOH4,LS,N2,N1)=trcsa_TBLS(idsa_FeOH4,LS,N2,N1)+XFE4BS(LS,N2,N1) &
      -XFE4BS(LS2,N2,N1)
      trcsa_TBLS(idsa_FeSO4,LS,N2,N1)=trcsa_TBLS(idsa_FeSO4,LS,N2,N1)+XFESBS(LS,N2,N1) &
      -XFESBS(LS2,N2,N1)
      trcsa_TBLS(idsa_CaOH2,LS,N2,N1)=trcsa_TBLS(idsa_CaOH2,LS,N2,N1)+XCAOBS(LS,N2,N1) &
      -XCAOBS(LS2,N2,N1)
      trcsa_TBLS(idsa_CaCO3,LS,N2,N1)=trcsa_TBLS(idsa_CaCO3,LS,N2,N1)+XCACBS(LS,N2,N1) &
      -XCACBS(LS2,N2,N1)
      trcsa_TBLS(idsa_CaHCO3,LS,N2,N1)=trcsa_TBLS(idsa_CaHCO3,LS,N2,N1)+XCAHBS(LS,N2,N1) &
      -XCAHBS(LS2,N2,N1)
      trcsa_TBLS(idsa_CaSO4,LS,N2,N1)=trcsa_TBLS(idsa_CaSO4,LS,N2,N1)+XCASBS(LS,N2,N1) &
      -XCASBS(LS2,N2,N1)
      trcsa_TBLS(idsa_MgOH2,LS,N2,N1)=trcsa_TBLS(idsa_MgOH2,LS,N2,N1)+XMGOBS(LS,N2,N1) &
      -XMGOBS(LS2,N2,N1)
      trcsa_TBLS(idsa_MgCO3,LS,N2,N1)=trcsa_TBLS(idsa_MgCO3,LS,N2,N1)+XMGCBS(LS,N2,N1) &
      -XMGCBS(LS2,N2,N1)
      trcsa_TBLS(idsa_MgHCO3,LS,N2,N1)=trcsa_TBLS(idsa_MgHCO3,LS,N2,N1)+XMGHBS(LS,N2,N1) &
      -XMGHBS(LS2,N2,N1)
      trcsa_TBLS(idsa_MgSO4,LS,N2,N1)=trcsa_TBLS(idsa_MgSO4,LS,N2,N1)+XMGSBS(LS,N2,N1) &
      -XMGSBS(LS2,N2,N1)
      trcsa_TBLS(idsa_NaCO3,LS,N2,N1)=trcsa_TBLS(idsa_NaCO3,LS,N2,N1)+XNACBS(LS,N2,N1) &
      -XNACBS(LS2,N2,N1)
      trcsa_TBLS(idsa_NaSO4,LS,N2,N1)=trcsa_TBLS(idsa_NaSO4,LS,N2,N1)+XNASBS(LS,N2,N1) &
      -XNASBS(LS2,N2,N1)
      trcsa_TBLS(idsa_KSO4,LS,N2,N1)=trcsa_TBLS(idsa_KSO4,LS,N2,N1)+XKASBS(LS,N2,N1) &
      -XKASBS(LS2,N2,N1)
      trcsa_TBLS(idsa_H0PO4,LS,N2,N1)=trcsa_TBLS(idsa_H0PO4,LS,N2,N1)+XH0PBS(LS,N2,N1) &
      -XH0PBS(LS2,N2,N1)
      trcsa_TBLS(idsa_H3PO4,LS,N2,N1)=trcsa_TBLS(idsa_H3PO4,LS,N2,N1)+XH3PBS(LS,N2,N1) &
      -XH3PBS(LS2,N2,N1)
      trcsa_TBLS(idsa_FeHPO4,LS,N2,N1)=trcsa_TBLS(idsa_FeHPO4,LS,N2,N1)+XF1PBS(LS,N2,N1) &
      -XF1PBS(LS2,N2,N1)
      trcsa_TBLS(idsa_FeH2PO4,LS,N2,N1)=trcsa_TBLS(idsa_FeH2PO4,LS,N2,N1)+XF2PBS(LS,N2,N1) &
      -XF2PBS(LS2,N2,N1)
      trcsa_TBLS(idsa_CaPO4,LS,N2,N1)=trcsa_TBLS(idsa_CaPO4,LS,N2,N1)+XC0PBS(LS,N2,N1) &
      -XC0PBS(LS2,N2,N1)
      trcsa_TBLS(idsa_CaHPO4,LS,N2,N1)=trcsa_TBLS(idsa_CaHPO4,LS,N2,N1)+XC1PBS(LS,N2,N1) &
      -XC1PBS(LS2,N2,N1)
      trcsa_TBLS(idsa_CaH2PO4,LS,N2,N1)=trcsa_TBLS(idsa_CaH2PO4,LS,N2,N1)+XC2PBS(LS,N2,N1) &
      -XC2PBS(LS2,N2,N1)
      trcsa_TBLS(idsa_MgHPO4,LS,N2,N1)=trcsa_TBLS(idsa_MgHPO4,LS,N2,N1)+XM1PBS(LS,N2,N1) &
      -XM1PBS(LS2,N2,N1)
      ENDIF
!
!     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
!
      ELSE
      TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)
      TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1) &
      -FLSWR(LS,N2,N1)-FLSW(LS,N2,N1)-FLSWH(LS,N2,N1)
      TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)
      THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1) &
      -HFLSWR(LS,N2,N1)-HFLSW(LS,N2,N1)
!     IF(LS.EQ.5)THEN
!     WRITE(*,7755)'LS',I,J,N1,N2,LS,LS2,TFLWW(LS,N2,N1),XFLWW(LS,N2,N1)
!    2,FLSWR(LS,N2,N1),FLSW(LS,N2,N1),FLSWH(LS,N2,N1)
!7755  FORMAT(A8,6I4,100E14.6)
!     ENDIF
      trcg_TBLS(idg_CO2,LS,N2,N1)=trcg_TBLS(idg_CO2,LS,N2,N1)+XCOBLS(LS,N2,N1) &
      -trcs_XFLS(idg_CO2,3,0,N2,N1)-trcs_XFLS(idg_CO2,3,NUM(N2,N1),N2,N1) &
      -XCOFHS(3,NUM(N2,N1),N2,N1)
      trcg_TBLS(idg_CH4,LS,N2,N1)=trcg_TBLS(idg_CH4,LS,N2,N1)+XCHBLS(LS,N2,N1) &
      -trcs_XFLS(idg_CH4,3,0,N2,N1)-trcs_XFLS(idg_CH4,3,NUM(N2,N1),N2,N1) &
      -XCHFHS(3,NUM(N2,N1),N2,N1)
      trcg_TBLS(idg_O2,LS,N2,N1)=trcg_TBLS(idg_O2,LS,N2,N1)+XOXBLS(LS,N2,N1) &
      -trcs_XFLS(idg_O2,3,0,N2,N1)-trcs_XFLS(idg_O2,3,NUM(N2,N1),N2,N1) &
      -XOXFHS(3,NUM(N2,N1),N2,N1)
      trcg_TBLS(idg_N2,LS,N2,N1)=trcg_TBLS(idg_N2,LS,N2,N1)+XNGBLS(LS,N2,N1) &
      -trcs_XFLS(idg_N2,3,0,N2,N1)-trcs_XFLS(idg_N2,3,NUM(N2,N1),N2,N1) &
      -XNGFHS(3,NUM(N2,N1),N2,N1)
      trcg_TBLS(idg_N2O,LS,N2,N1)=trcg_TBLS(idg_N2O,LS,N2,N1)+XN2BLS(LS,N2,N1) &
      -trcs_XFLS(idg_N2O,3,0,N2,N1)-trcs_XFLS(idg_N2O,3,NUM(N2,N1),N2,N1) &
      -XN2FHS(3,NUM(N2,N1),N2,N1)
      trcn_TBLS(ids_NH4,LS,N2,N1)=trcn_TBLS(ids_NH4,LS,N2,N1)+XN4BLW(LS,N2,N1) &
      -trcs_XFLS(ids_NH4,3,0,N2,N1)-trcs_XFLS(ids_NH4,3,NUM(N2,N1),N2,N1) &
      -XN4FHW(3,NUM(N2,N1),N2,N1)-trcs_XFLS(ids_NH4B,3,NUM(N2,N1),N2,N1) &
      -XN4FHB(3,NUM(N2,N1),N2,N1)
      trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1)+XN3BLW(LS,N2,N1) &
      -trcs_XFLS(idg_NH3,3,0,N2,N1)-trcs_XFLS(idg_NH3,3,NUM(N2,N1),N2,N1) &
      -XN3FHW(3,NUM(N2,N1),N2,N1)-trcs_XFLS(idg_NH3B,3,NUM(N2,N1),N2,N1) &
      -XN3FHB(3,NUM(N2,N1),N2,N1)
      trcn_TBLS(ids_NO3,LS,N2,N1)=trcn_TBLS(ids_NO3,LS,N2,N1)+XNOBLW(LS,N2,N1) &
      -trcs_XFLS(ids_NO3,3,0,N2,N1)-trcs_XFLS(ids_NO3,3,NUM(N2,N1),N2,N1) &
      -XNOFHW(3,NUM(N2,N1),N2,N1)-trcs_XFLS(ids_NO3B,3,NUM(N2,N1),N2,N1) &
      -XNOFHB(3,NUM(N2,N1),N2,N1)
      trcn_TBLS(ids_H1PO4,LS,N2,N1)=trcn_TBLS(ids_H1PO4,LS,N2,N1)+XH1PBS(LS,N2,N1) &
      -trcs_XFLS(ids_H1PO4,3,0,N2,N1)-trcs_XFLS(ids_H1PO4,3,NUM(N2,N1),N2,N1) &
      -XH1PHS(3,NUM(N2,N1),N2,N1)-trcs_XFLS(ids_H1PO4B,3,NUM(N2,N1),N2,N1) &
      -XH1BHB(3,NUM(N2,N1),N2,N1)
      trcn_TBLS(ids_H2PO4,LS,N2,N1)=trcn_TBLS(ids_H2PO4,LS,N2,N1)+XH2PBS(LS,N2,N1) &
      -trcs_XFLS(ids_H2PO4,3,0,N2,N1)-trcs_XFLS(ids_H2PO4,3,NUM(N2,N1),N2,N1) &
      -XH2PHS(3,NUM(N2,N1),N2,N1)-trcs_XFLS(ids_H2PO4B,3,NUM(N2,N1),N2,N1) &
      -XH2BHB(3,NUM(N2,N1),N2,N1)
      IF(ISALTG.NE.0)THEN
      trcsa_TBLS(idsa_Al,LS,NY,NX)=trcsa_TBLS(idsa_Al,LS,NY,NX)+XALBLS(LS,NY,NX) &
      -XALFLS(3,0,N2,N1)-XALFLS(3,NUM(N2,N1),N2,N1) &
      -XALFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Fe,LS,NY,NX)=trcsa_TBLS(idsa_Fe,LS,NY,NX)+XFEBLS(LS,NY,NX) &
      -XFEFLS(3,0,N2,N1)-XFEFLS(3,NUM(N2,N1),N2,N1) &
      -XFEFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Hp,LS,NY,NX)=trcsa_TBLS(idsa_Hp,LS,NY,NX)+XHYBLS(LS,NY,NX) &
      -XHYFLS(3,0,N2,N1)-XHYFLS(3,NUM(N2,N1),N2,N1) &
      -XHYFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Ca,LS,NY,NX)=trcsa_TBLS(idsa_Ca,LS,NY,NX)+XCABLS(LS,NY,NX) &
      -XCAFLS(3,0,N2,N1)-XCAFLS(3,NUM(N2,N1),N2,N1) &
      -XCAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Mg,LS,NY,NX)=trcsa_TBLS(idsa_Mg,LS,NY,NX)+XMGBLS(LS,NY,NX) &
      -XMGFLS(3,0,N2,N1)-XMGFLS(3,NUM(N2,N1),N2,N1) &
      -XMGFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Na,LS,NY,NX)=trcsa_TBLS(idsa_Na,LS,NY,NX)+XNABLS(LS,NY,NX) &
      -XNAFLS(3,0,N2,N1)-XNAFLS(3,NUM(N2,N1),N2,N1) &
      -XNAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_K,LS,NY,NX)=trcsa_TBLS(idsa_K,LS,NY,NX)+XKABLS(LS,NY,NX) &
      -XKAFLS(3,0,N2,N1)-XKAFLS(3,NUM(N2,N1),N2,N1) &
      -XKAFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_OH,LS,NY,NX)=trcsa_TBLS(idsa_OH,LS,NY,NX)+XOHBLS(LS,NY,NX) &
      -XOHFLS(3,0,N2,N1)-XOHFLS(3,NUM(N2,N1),N2,N1) &
      -XOHFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_SO4,LS,NY,NX)=trcsa_TBLS(idsa_SO4,LS,NY,NX)+XSOBLS(LS,NY,NX) &
      -XSOFLS(3,0,N2,N1)-XSOFLS(3,NUM(N2,N1),N2,N1) &
      -XSOFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_Cl,LS,NY,NX)=trcsa_TBLS(idsa_Cl,LS,NY,NX)+XCLBLS(LS,NY,NX) &
      -XCLFLS(3,0,N2,N1)-XCLFLS(3,NUM(N2,N1),N2,N1) &
      -XCLFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CO3,LS,NY,NX)=trcsa_TBLS(idsa_CO3,LS,NY,NX)+XC3BLS(LS,NY,NX) &
      -XC3FLS(3,0,N2,N1)-XC3FLS(3,NUM(N2,N1),N2,N1) &
      -XC3FHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_HCO3,LS,NY,NX)=trcsa_TBLS(idsa_HCO3,LS,NY,NX)+XHCBLS(LS,NY,NX) &
      -XHCFLS(3,0,N2,N1)-XHCFLS(3,NUM(N2,N1),N2,N1) &
      -XHCFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH,LS,NY,NX)=trcsa_TBLS(idsa_AlOH,LS,NY,NX)+XAL1BS(LS,NY,NX) &
      -XAL1FS(3,0,N2,N1)-XAL1FS(3,NUM(N2,N1),N2,N1) &
      -XAL1HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH2,LS,NY,NX)=trcsa_TBLS(idsa_AlOH2,LS,NY,NX)+XAL2BS(LS,NY,NX) &
      -XAL2FS(3,0,N2,N1)-XAL2FS(3,NUM(N2,N1),N2,N1) &
      -XAL2HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH3,LS,NY,NX)=trcsa_TBLS(idsa_AlOH3,LS,NY,NX)+XAL3BS(LS,NY,NX) &
      -XAL3FS(3,0,N2,N1)-XAL3FS(3,NUM(N2,N1),N2,N1) &
      -XAL3HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlOH4,LS,NY,NX)=trcsa_TBLS(idsa_AlOH4,LS,NY,NX)+XAL4BS(LS,NY,NX) &
      -XAL4FS(3,0,N2,N1)-XAL4FS(3,NUM(N2,N1),N2,N1) &
      -XAL4HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_AlSO4,LS,NY,NX)=trcsa_TBLS(idsa_AlSO4,LS,NY,NX)+XALSBS(LS,NY,NX) &
      -XALSFS(3,0,N2,N1)-XALSFS(3,NUM(N2,N1),N2,N1) &
      -XALSHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH,LS,NY,NX)=trcsa_TBLS(idsa_FeOH,LS,NY,NX)+XFE1BS(LS,NY,NX) &
      -XFE1FS(3,0,N2,N1)-XFE1FS(3,NUM(N2,N1),N2,N1) &
      -XFE1HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH2,LS,NY,NX)=trcsa_TBLS(idsa_FeOH2,LS,NY,NX)+XFE2BS(LS,NY,NX) &
      -XFE2FS(3,0,N2,N1)-XFE2FS(3,NUM(N2,N1),N2,N1) &
      -XFE2HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH3,LS,NY,NX)=trcsa_TBLS(idsa_FeOH3,LS,NY,NX)+XFE3BS(LS,NY,NX) &
      -XFE3FS(3,0,N2,N1)-XFE3FS(3,NUM(N2,N1),N2,N1) &
      -XFE3HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeOH4,LS,NY,NX)=trcsa_TBLS(idsa_FeOH4,LS,NY,NX)+XFE4BS(LS,NY,NX) &
      -XFE4FS(3,0,N2,N1)-XFE4FS(3,NUM(N2,N1),N2,N1) &
      -XFE4HS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeSO4,LS,NY,NX)=trcsa_TBLS(idsa_FeSO4,LS,NY,NX)+XFESBS(LS,NY,NX) &
      -XFESFS(3,0,N2,N1)-XFESFS(3,NUM(N2,N1),N2,N1) &
      -XFESHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaOH2,LS,NY,NX)=trcsa_TBLS(idsa_CaOH2,LS,NY,NX)+XCAOBS(LS,NY,NX) &
      -XCAOFS(3,0,N2,N1)-XCAOFS(3,NUM(N2,N1),N2,N1) &
      -XCAOHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaCO3,LS,NY,NX)+XCACBS(LS,NY,NX) &
      -XCACFS(3,0,N2,N1)-XCACFS(3,NUM(N2,N1),N2,N1) &
      -XCACHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)=trcsa_TBLS(idsa_CaHCO3,LS,NY,NX)+XCAHBS(LS,NY,NX) &
      -XCAHFS(3,0,N2,N1)-XCAHFS(3,NUM(N2,N1),N2,N1) &
      -XALFHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaSO4,LS,NY,NX)=trcsa_TBLS(idsa_CaSO4,LS,NY,NX)+XCASBS(LS,NY,NX) &
      -XCASFS(3,0,N2,N1)-XCASFS(3,NUM(N2,N1),N2,N1) &
      -XCASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgOH2,LS,NY,NX)=trcsa_TBLS(idsa_MgOH2,LS,NY,NX)+XMGOBS(LS,NY,NX) &
      -XMGOFS(3,0,N2,N1)-XMGOFS(3,NUM(N2,N1),N2,N1) &
      -XMGOHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgCO3,LS,NY,NX)+XMGCBS(LS,NY,NX) &
      -XMGCFS(3,0,N2,N1)-XMGCFS(3,NUM(N2,N1),N2,N1) &
      -XMGCHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)=trcsa_TBLS(idsa_MgHCO3,LS,NY,NX)+XMGHBS(LS,NY,NX) &
      -XMGHFS(3,0,N2,N1)-XMGHFS(3,NUM(N2,N1),N2,N1) &
      -XMGHHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgSO4,LS,NY,NX)=trcsa_TBLS(idsa_MgSO4,LS,NY,NX)+XMGSBS(LS,NY,NX) &
      -XMGSFS(3,0,N2,N1)-XMGSFS(3,NUM(N2,N1),N2,N1) &
      -XMGSHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_NaCO3,LS,NY,NX)=trcsa_TBLS(idsa_NaCO3,LS,NY,NX)+XNACBS(LS,NY,NX) &
      -XNACFS(3,0,N2,N1)-XNACFS(3,NUM(N2,N1),N2,N1) &
      -XNACHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_NaSO4,LS,NY,NX)=trcsa_TBLS(idsa_NaSO4,LS,NY,NX)+XNASBS(LS,NY,NX) &
      -XNASFS(3,0,N2,N1)-XNASFS(3,NUM(N2,N1),N2,N1) &
      -XNASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_KSO4,LS,NY,NX)=trcsa_TBLS(idsa_KSO4,LS,NY,NX)+XKASBS(LS,NY,NX) &
      -XKASFS(3,0,N2,N1)-XKASFS(3,NUM(N2,N1),N2,N1) &
      -XKASHS(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_H0PO4,LS,NY,NX)=trcsa_TBLS(idsa_H0PO4,LS,NY,NX)+XH0PBS(LS,NY,NX) &
      -XH0PFS(3,0,N2,N1)-XH0PFS(3,NUM(N2,N1),N2,N1) &
      -XH0PHS(3,NUM(N2,N1),N2,N1)-XH0BFB(3,NUM(N2,N1),N2,N1) &
      -XH0BFB(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_H3PO4,LS,NY,NX)=trcsa_TBLS(idsa_H3PO4,LS,NY,NX)+XH3PBS(LS,NY,NX) &
      -XH3PFS(3,0,N2,N1)-XH3PHS(3,NUM(N2,N1),N2,N1) &
      -XH3PHS(3,NUM(N2,N1),N2,N1)-XH3BFB(3,NUM(N2,N1),N2,N1) &
      -XH3BFB(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)=trcsa_TBLS(idsa_FeHPO4,LS,NY,NX)+XF1PBS(LS,NY,NX) &
      -XF1PFS(3,0,N2,N1)-XF1PFS(3,NUM(N2,N1),N2,N1) &
      -XF1PHS(3,NUM(N2,N1),N2,N1)-XF1BFB(3,NUM(N2,N1),N2,N1) &
      -XF1BFB(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_FeH2PO4,LS,NY,NX)+XF2PBS(LS,NY,NX) &
      -XF2PFS(3,0,N2,N1)-XF2PFS(3,NUM(N2,N1),N2,N1) &
      -XF2PHS(3,NUM(N2,N1),N2,N1)-XF2BFB(3,NUM(N2,N1),N2,N1) &
      -XF2BFB(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaPO4,LS,NY,NX)+XC0PBS(LS,NY,NX) &
      -XC0PFS(3,0,N2,N1)-XC0PFS(3,NUM(N2,N1),N2,N1) &
      -XC0PHS(3,NUM(N2,N1),N2,N1)-XC0BFB(3,NUM(N2,N1),N2,N1) &
      -XC0BFB(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)=trcsa_TBLS(idsa_CaHPO4,LS,NY,NX)+XC1PBS(LS,NY,NX) &
      -XC1PFS(3,0,N2,N1)-XC1PFS(3,NUM(N2,N1),N2,N1) &
      -XC1PHS(3,NUM(N2,N1),N2,N1)-XC1BFB(3,NUM(N2,N1),N2,N1) &
      -XC1BFB(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)=trcsa_TBLS(idsa_CaH2PO4,LS,NY,NX)+XC2PBS(LS,NY,NX) &
      -XC2PFS(3,0,N2,N1)-XC2PFS(3,NUM(N2,N1),N2,N1) &
      -XC2PHS(3,NUM(N2,N1),N2,N1)-XC2BFB(3,NUM(N2,N1),N2,N1) &
      -XC2BFB(3,NUM(N2,N1),N2,N1)
      trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)=trcsa_TBLS(idsa_MgHPO4,LS,NY,NX)+XM1PBS(LS,NY,NX) &
      -XM1PFS(3,0,N2,N1)-XM1PFS(3,NUM(N2,N1),N2,N1) &
      -XM1PHS(3,NUM(N2,N1),N2,N1)-XM1BFB(3,NUM(N2,N1),N2,N1) &
      -XM1BFB(3,NUM(N2,N1),N2,N1)
      ENDIF
      ENDIF
!
!     WATER,GAS,SOLUTE,SALT FLUXES INTO SNOWPACK SURFACE
!
      ELSEIF(LS.EQ.1)THEN
      IF(abs(XFLWS(LS,N2,N1))>0._r8)THEN
      TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)
      TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1)
      TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)
      THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1)

      trcg_TBLS(idg_CO2,LS,N2,N1)=trcg_TBLS(idg_CO2,LS,N2,N1)+XCOBLS(LS,N2,N1)
      trcg_TBLS(idg_CH4,LS,N2,N1)=trcg_TBLS(idg_CH4,LS,N2,N1)+XCHBLS(LS,N2,N1)
      trcg_TBLS(idg_O2,LS,N2,N1)=trcg_TBLS(idg_O2,LS,N2,N1)+XOXBLS(LS,N2,N1)
      trcg_TBLS(idg_N2,LS,N2,N1)=trcg_TBLS(idg_N2,LS,N2,N1)+XNGBLS(LS,N2,N1)
      trcg_TBLS(idg_N2O,LS,N2,N1)=trcg_TBLS(idg_N2O,LS,N2,N1)+XN2BLS(LS,N2,N1)
      trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1)+XN3BLW(LS,N2,N1)

      trcn_TBLS(ids_NH4,LS,N2,N1)=trcn_TBLS(ids_NH4,LS,N2,N1)+XN4BLW(LS,N2,N1)
      trcn_TBLS(ids_NO3,LS,N2,N1)=trcn_TBLS(ids_NO3,LS,N2,N1)+XNOBLW(LS,N2,N1)
      trcn_TBLS(ids_H1PO4,LS,N2,N1)=trcn_TBLS(ids_H1PO4,LS,N2,N1)+XH1PBS(LS,N2,N1)
      trcn_TBLS(ids_H2PO4,LS,N2,N1)=trcn_TBLS(ids_H2PO4,LS,N2,N1)+XH2PBS(LS,N2,N1)

      IF(ISALTG.NE.0)THEN
      trcsa_TBLS(idsa_Al,LS,N2,N1)=trcsa_TBLS(idsa_Al,LS,N2,N1)+XALBLS(LS,N2,N1)
      trcsa_TBLS(idsa_Fe,LS,N2,N1)=trcsa_TBLS(idsa_Fe,LS,N2,N1)+XFEBLS(LS,N2,N1)
      trcsa_TBLS(idsa_Hp,LS,N2,N1)=trcsa_TBLS(idsa_Hp,LS,N2,N1)+XHYBLS(LS,N2,N1)
      trcsa_TBLS(idsa_Ca,LS,N2,N1)=trcsa_TBLS(idsa_Ca,LS,N2,N1)+XCABLS(LS,N2,N1)
      trcsa_TBLS(idsa_Mg,LS,N2,N1)=trcsa_TBLS(idsa_Mg,LS,N2,N1)+XMGBLS(LS,N2,N1)
      trcsa_TBLS(idsa_Na,LS,N2,N1)=trcsa_TBLS(idsa_Na,LS,N2,N1)+XNABLS(LS,N2,N1)
      trcsa_TBLS(idsa_K,LS,N2,N1)=trcsa_TBLS(idsa_K,LS,N2,N1)+XKABLS(LS,N2,N1)
      trcsa_TBLS(idsa_OH,LS,N2,N1)=trcsa_TBLS(idsa_OH,LS,N2,N1)+XOHBLS(LS,N2,N1)
      trcsa_TBLS(idsa_SO4,LS,N2,N1)=trcsa_TBLS(idsa_SO4,LS,N2,N1)+XSOBLS(LS,N2,N1)
      trcsa_TBLS(idsa_Cl,LS,N2,N1)=trcsa_TBLS(idsa_Cl,LS,N2,N1)+XCLBLS(LS,N2,N1)
      trcsa_TBLS(idsa_CO3,LS,N2,N1)=trcsa_TBLS(idsa_CO3,LS,N2,N1)+XC3BLS(LS,N2,N1)
      trcsa_TBLS(idsa_HCO3,LS,N2,N1)=trcsa_TBLS(idsa_HCO3,LS,N2,N1)+XHCBLS(LS,N2,N1)
      trcsa_TBLS(idsa_AlOH,LS,N2,N1)=trcsa_TBLS(idsa_AlOH,LS,N2,N1)+XAL1BS(LS,N2,N1)
      trcsa_TBLS(idsa_AlOH2,LS,N2,N1)=trcsa_TBLS(idsa_AlOH2,LS,N2,N1)+XAL2BS(LS,N2,N1)
      trcsa_TBLS(idsa_AlOH3,LS,N2,N1)=trcsa_TBLS(idsa_AlOH3,LS,N2,N1)+XAL3BS(LS,N2,N1)
      trcsa_TBLS(idsa_AlOH4,LS,N2,N1)=trcsa_TBLS(idsa_AlOH4,LS,N2,N1)+XAL4BS(LS,N2,N1)
      trcsa_TBLS(idsa_AlSO4,LS,N2,N1)=trcsa_TBLS(idsa_AlSO4,LS,N2,N1)+XALSBS(LS,N2,N1)
      trcsa_TBLS(idsa_FeOH,LS,N2,N1)=trcsa_TBLS(idsa_FeOH,LS,N2,N1)+XFE1BS(LS,N2,N1)
      trcsa_TBLS(idsa_FeOH2,LS,N2,N1)=trcsa_TBLS(idsa_FeOH2,LS,N2,N1)+XFE2BS(LS,N2,N1)
      trcsa_TBLS(idsa_FeOH3,LS,N2,N1)=trcsa_TBLS(idsa_FeOH3,LS,N2,N1)+XFE3BS(LS,N2,N1)
      trcsa_TBLS(idsa_FeOH4,LS,N2,N1)=trcsa_TBLS(idsa_FeOH4,LS,N2,N1)+XFE4BS(LS,N2,N1)
      trcsa_TBLS(idsa_FeSO4,LS,N2,N1)=trcsa_TBLS(idsa_FeSO4,LS,N2,N1)+XFESBS(LS,N2,N1)
      trcsa_TBLS(idsa_CaOH2,LS,N2,N1)=trcsa_TBLS(idsa_CaOH2,LS,N2,N1)+XCAOBS(LS,N2,N1)
      trcsa_TBLS(idsa_CaCO3,LS,N2,N1)=trcsa_TBLS(idsa_CaCO3,LS,N2,N1)+XCACBS(LS,N2,N1)
      trcsa_TBLS(idsa_CaHCO3,LS,N2,N1)=trcsa_TBLS(idsa_CaHCO3,LS,N2,N1)+XCAHBS(LS,N2,N1)
      trcsa_TBLS(idsa_CaSO4,LS,N2,N1)=trcsa_TBLS(idsa_CaSO4,LS,N2,N1)+XCASBS(LS,N2,N1)
      trcsa_TBLS(idsa_MgOH2,LS,N2,N1)=trcsa_TBLS(idsa_MgOH2,LS,N2,N1)+XMGOBS(LS,N2,N1)
      trcsa_TBLS(idsa_MgCO3,LS,N2,N1)=trcsa_TBLS(idsa_MgCO3,LS,N2,N1)+XMGCBS(LS,N2,N1)
      trcsa_TBLS(idsa_MgHCO3,LS,N2,N1)=trcsa_TBLS(idsa_MgHCO3,LS,N2,N1)+XMGHBS(LS,N2,N1)
      trcsa_TBLS(idsa_MgSO4,LS,N2,N1)=trcsa_TBLS(idsa_MgSO4,LS,N2,N1)+XMGSBS(LS,N2,N1)
      trcsa_TBLS(idsa_NaCO3,LS,N2,N1)=trcsa_TBLS(idsa_NaCO3,LS,N2,N1)+XNACBS(LS,N2,N1)
      trcsa_TBLS(idsa_NaSO4,LS,N2,N1)=trcsa_TBLS(idsa_NaSO4,LS,N2,N1)+XNASBS(LS,N2,N1)
      trcsa_TBLS(idsa_KSO4,LS,N2,N1)=trcsa_TBLS(idsa_KSO4,LS,N2,N1)+XKASBS(LS,N2,N1)
      trcsa_TBLS(idsa_H0PO4,LS,N2,N1)=trcsa_TBLS(idsa_H0PO4,LS,N2,N1)+XH0PBS(LS,N2,N1)
      trcsa_TBLS(idsa_H3PO4,LS,N2,N1)=trcsa_TBLS(idsa_H3PO4,LS,N2,N1)+XH3PBS(LS,N2,N1)
      trcsa_TBLS(idsa_FeHPO4,LS,N2,N1)=trcsa_TBLS(idsa_FeHPO4,LS,N2,N1)+XF1PBS(LS,N2,N1)
      trcsa_TBLS(idsa_FeH2PO4,LS,N2,N1)=trcsa_TBLS(idsa_FeH2PO4,LS,N2,N1)+XF2PBS(LS,N2,N1)
      trcsa_TBLS(idsa_CaPO4,LS,N2,N1)=trcsa_TBLS(idsa_CaPO4,LS,N2,N1)+XC0PBS(LS,N2,N1)
      trcsa_TBLS(idsa_CaHPO4,LS,N2,N1)=trcsa_TBLS(idsa_CaHPO4,LS,N2,N1)+XC1PBS(LS,N2,N1)
      trcsa_TBLS(idsa_CaH2PO4,LS,N2,N1)=trcsa_TBLS(idsa_CaH2PO4,LS,N2,N1)+XC2PBS(LS,N2,N1)
      trcsa_TBLS(idsa_MgHPO4,LS,N2,N1)=trcsa_TBLS(idsa_MgHPO4,LS,N2,N1)+XM1PBS(LS,N2,N1)
      ENDIF
      ENDIF
    ENDIF
  ENDDO D1205
  end subroutine SaltFromRunoffSnowpack
!------------------------------------------------------------------------------------------

  subroutine TotalFluxFromSedmentTransp(N,N1,N2,N4 &
      ,N5,N4B,N5B,NY,NX)
      implicit none
      integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B,NY,NX

      integer :: M,K,NO,NN,NGL
!     begin_execution
!
!     T*ER=net sediment flux
!     X*ER=sediment flux from erosion.f
!     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
!       :OMC,OMN,OMP=microbial C,N,P; ORC=microbial residue C,N,P
!       :OHC,OHN,OHP=adsorbed C,N,P; OSC,OSN,OSP=humus C,N,P
!       :NH4,NH3,NHU,NO3=fertilizer NH4,NH3,urea,NO3 in non-band
!       :NH4B,NH3B,NHUB,NO3B=fertilizer NH4,NH3,urea,NO3 in band
!       :XN4,XNB=adsorbed NH4 in non-band,band
!       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
!        =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
!       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
!       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
!       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
!       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
!       :PALO,PFEO=precip AlOH,FeOH
!       :PCAC,PCAS=precip CaCO3,CaSO4
!       :PALP,PFEP=precip AlPO4,FEPO4 in non-band
!       :PALPB,PFEPB=precip AlPO4,FEPO4 in band
!       :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
!       :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
!
  IF(N.NE.3.AND.(IERSNG.EQ.1.OR.IERSNG.EQ.3))THEN
    D9350: DO NN=1,2
      IF(ABS(XSEDER(N,NN,N2,N1)).GT.ZEROS(N2,N1) &
        .OR.ABS(XSEDER(N,NN,N5,N4)).GT.ZEROS(N5,N4))THEN
        TSEDER(N2,N1)=TSEDER(N2,N1)+XSEDER(N,NN,N2,N1)
        TSANER(N2,N1)=TSANER(N2,N1)+XSANER(N,NN,N2,N1)
        TSILER(N2,N1)=TSILER(N2,N1)+XSILER(N,NN,N2,N1)
        TCLAER(N2,N1)=TCLAER(N2,N1)+XCLAER(N,NN,N2,N1)
        TCECER(N2,N1)=TCECER(N2,N1)+XCECER(N,NN,N2,N1)
        TAECER(N2,N1)=TAECER(N2,N1)+XAECER(N,NN,N2,N1)
        TNH4ER(N2,N1)=TNH4ER(N2,N1)+XNH4ER(N,NN,N2,N1)
        TNH3ER(N2,N1)=TNH3ER(N2,N1)+XNH3ER(N,NN,N2,N1)
        TNHUER(N2,N1)=TNHUER(N2,N1)+XNHUER(N,NN,N2,N1)
        TNO3ER(N2,N1)=TNO3ER(N2,N1)+XNO3ER(N,NN,N2,N1)
        TNH4EB(N2,N1)=TNH4EB(N2,N1)+XNH4EB(N,NN,N2,N1)
        TNH3EB(N2,N1)=TNH3EB(N2,N1)+XNH3EB(N,NN,N2,N1)
        TNHUEB(N2,N1)=TNHUEB(N2,N1)+XNHUEB(N,NN,N2,N1)
        TNO3EB(N2,N1)=TNO3EB(N2,N1)+XNO3EB(N,NN,N2,N1)
        TN4ER(N2,N1)=TN4ER(N2,N1)+XN4ER(N,NN,N2,N1)
        TNBER(N2,N1)=TNBER(N2,N1)+XNBER(N,NN,N2,N1)
        THYER(N2,N1)=THYER(N2,N1)+XHYER(N,NN,N2,N1)
        TALER(N2,N1)=TALER(N2,N1)+XALER(N,NN,N2,N1)
        TFEER(N2,N1)=TFEER(N2,N1)+XFEER(N,NN,N2,N1)
        TCAER(N2,N1)=TCAER(N2,N1)+XCAER(N,NN,N2,N1)
        TMGER(N2,N1)=TMGER(N2,N1)+XMGER(N,NN,N2,N1)
        TNAER(N2,N1)=TNAER(N2,N1)+XNAER(N,NN,N2,N1)
        TKAER(N2,N1)=TKAER(N2,N1)+XKAER(N,NN,N2,N1)
        THCER(N2,N1)=THCER(N2,N1)+XHCER(N,NN,N2,N1)
        TAL2ER(N2,N1)=TAL2ER(N2,N1)+XAL2ER(N,NN,N2,N1)
        TFE2ER(N2,N1)=TFE2ER(N2,N1)+XFE2ER(N,NN,N2,N1)
        TOH0ER(N2,N1)=TOH0ER(N2,N1)+XOH0ER(N,NN,N2,N1)
        TOH1ER(N2,N1)=TOH1ER(N2,N1)+XOH1ER(N,NN,N2,N1)
        TOH2ER(N2,N1)=TOH2ER(N2,N1)+XOH2ER(N,NN,N2,N1)
        TH1PER(N2,N1)=TH1PER(N2,N1)+XH1PER(N,NN,N2,N1)
        TH2PER(N2,N1)=TH2PER(N2,N1)+XH2PER(N,NN,N2,N1)
        TOH0EB(N2,N1)=TOH0EB(N2,N1)+XOH0EB(N,NN,N2,N1)
        TOH1EB(N2,N1)=TOH1EB(N2,N1)+XOH1EB(N,NN,N2,N1)
        TOH2EB(N2,N1)=TOH2EB(N2,N1)+XOH2EB(N,NN,N2,N1)
        TH1PEB(N2,N1)=TH1PEB(N2,N1)+XH1PEB(N,NN,N2,N1)
        TH2PEB(N2,N1)=TH2PEB(N2,N1)+XH2PEB(N,NN,N2,N1)
        TALOER(N2,N1)=TALOER(N2,N1)+PALOER(N,NN,N2,N1)
        TFEOER(N2,N1)=TFEOER(N2,N1)+PFEOER(N,NN,N2,N1)
        TCACER(N2,N1)=TCACER(N2,N1)+PCACER(N,NN,N2,N1)
        TCASER(N2,N1)=TCASER(N2,N1)+PCASER(N,NN,N2,N1)
        TALPER(N2,N1)=TALPER(N2,N1)+PALPER(N,NN,N2,N1)
        TFEPER(N2,N1)=TFEPER(N2,N1)+PFEPER(N,NN,N2,N1)
        TCPDER(N2,N1)=TCPDER(N2,N1)+PCPDER(N,NN,N2,N1)
        TCPHER(N2,N1)=TCPHER(N2,N1)+PCPHER(N,NN,N2,N1)
        TCPMER(N2,N1)=TCPMER(N2,N1)+PCPMER(N,NN,N2,N1)
        TALPEB(N2,N1)=TALPEB(N2,N1)+PALPEB(N,NN,N2,N1)
        TFEPEB(N2,N1)=TFEPEB(N2,N1)+PFEPEB(N,NN,N2,N1)
        TCPDEB(N2,N1)=TCPDEB(N2,N1)+PCPDEB(N,NN,N2,N1)
        TCPHEB(N2,N1)=TCPHEB(N2,N1)+PCPHEB(N,NN,N2,N1)
        TCPMEB(N2,N1)=TCPMEB(N2,N1)+PCPMEB(N,NN,N2,N1)
        DO  K=1,jcplx
          DO NO=1,NFGs
            DO NGL=JGnio(NO),JGnfo(NO)
              DO M=1,nlbiomcp
                TOMCER(M,NGL,K,N2,N1)=TOMCER(M,NGL,K,N2,N1)+OMCER(M+(NGL-1)*nlbiomcp,K,N,NN,N2,N1)
                TOMNER(M,NGL,K,N2,N1)=TOMNER(M,NGL,K,N2,N1)+OMNER(M+(NGL-1)*nlbiomcp,K,N,NN,N2,N1)
                TOMPER(M,NGL,K,N2,N1)=TOMPER(M,NGL,K,N2,N1)+OMPER(M+(NGL-1)*nlbiomcp,K,N,NN,N2,N1)
              enddo
            enddo
          enddo
        ENDDO

        DO NO=1,NFGs
          DO NGL=JGniA(NO),JGnfA(NO)
            DO M=1,nlbiomcp
              TOMCERff(M,NGL,N2,N1)=TOMCERff(M,NGL,N2,N1)+OMCERff(M+(NGL-1)*nlbiomcp,N,NN,N2,N1)
              TOMNERff(M,NGL,N2,N1)=TOMNERff(M,NGL,N2,N1)+OMNERff(M+(NGL-1)*nlbiomcp,N,NN,N2,N1)
              TOMPERff(M,NGL,N2,N1)=TOMPERff(M,NGL,N2,N1)+OMPERff(M+(NGL-1)*nlbiomcp,N,NN,N2,N1)
            enddo
          enddo
        enddo


        D9375: DO K=1,jcplx
          D9370: DO M=1,ndbiomcp
            TORCER(M,K,N2,N1)=TORCER(M,K,N2,N1)+ORCER(M,K,N,NN,N2,N1)
            TORNER(M,K,N2,N1)=TORNER(M,K,N2,N1)+ORNER(M,K,N,NN,N2,N1)
            TORPER(M,K,N2,N1)=TORPER(M,K,N2,N1)+ORPER(M,K,N,NN,N2,N1)
          ENDDO D9370
          TOHCER(K,N2,N1)=TOHCER(K,N2,N1)+OHCER(K,N,NN,N2,N1)
          TOHNER(K,N2,N1)=TOHNER(K,N2,N1)+OHNER(K,N,NN,N2,N1)
          TOHPER(K,N2,N1)=TOHPER(K,N2,N1)+OHPER(K,N,NN,N2,N1)
          TOHAER(K,N2,N1)=TOHAER(K,N2,N1)+OHAER(K,N,NN,N2,N1)
          D9365: DO M=1,jsken
            TOSCER(M,K,N2,N1)=TOSCER(M,K,N2,N1)+OSCER(M,K,N,NN,N2,N1)
            TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)+OSAER(M,K,N,NN,N2,N1)
            TOSNER(M,K,N2,N1)=TOSNER(M,K,N2,N1)+OSNER(M,K,N,NN,N2,N1)
            TOSPER(M,K,N2,N1)=TOSPER(M,K,N2,N1)+OSPER(M,K,N,NN,N2,N1)
          ENDDO D9365
        ENDDO D9375
!     IF(NN.EQ.2)THEN
        TSEDER(N2,N1)=TSEDER(N2,N1)-XSEDER(N,NN,N5,N4)
        TSANER(N2,N1)=TSANER(N2,N1)-XSANER(N,NN,N5,N4)
        TSILER(N2,N1)=TSILER(N2,N1)-XSILER(N,NN,N5,N4)
        TCLAER(N2,N1)=TCLAER(N2,N1)-XCLAER(N,NN,N5,N4)
        TCECER(N2,N1)=TCECER(N2,N1)-XCECER(N,NN,N5,N4)
        TAECER(N2,N1)=TAECER(N2,N1)-XAECER(N,NN,N5,N4)
        TNH4ER(N2,N1)=TNH4ER(N2,N1)-XNH4ER(N,NN,N5,N4)
        TNH3ER(N2,N1)=TNH3ER(N2,N1)-XNH3ER(N,NN,N5,N4)
        TNHUER(N2,N1)=TNHUER(N2,N1)-XNHUER(N,NN,N5,N4)
        TNO3ER(N2,N1)=TNO3ER(N2,N1)-XNO3ER(N,NN,N5,N4)
        TNH4EB(N2,N1)=TNH4EB(N2,N1)-XNH4EB(N,NN,N5,N4)
        TNH3EB(N2,N1)=TNH3EB(N2,N1)-XNH3EB(N,NN,N5,N4)
        TNHUEB(N2,N1)=TNHUEB(N2,N1)-XNHUEB(N,NN,N5,N4)
        TNO3EB(N2,N1)=TNO3EB(N2,N1)-XNO3EB(N,NN,N5,N4)
        TN4ER(N2,N1)=TN4ER(N2,N1)-XN4ER(N,NN,N5,N4)
        TNBER(N2,N1)=TNBER(N2,N1)-XNBER(N,NN,N5,N4)
        THYER(N2,N1)=THYER(N2,N1)-XHYER(N,NN,N5,N4)
        TALER(N2,N1)=TALER(N2,N1)-XALER(N,NN,N5,N4)
        TFEER(N2,N1)=TFEER(N2,N1)-XFEER(N,NN,N5,N4)
        TCAER(N2,N1)=TCAER(N2,N1)-XCAER(N,NN,N5,N4)
        TMGER(N2,N1)=TMGER(N2,N1)-XMGER(N,NN,N5,N4)
        TNAER(N2,N1)=TNAER(N2,N1)-XNAER(N,NN,N5,N4)
        TKAER(N2,N1)=TKAER(N2,N1)-XKAER(N,NN,N5,N4)
        THCER(N2,N1)=THCER(N2,N1)-XHCER(N,NN,N5,N4)
        TAL2ER(N2,N1)=TAL2ER(N2,N1)-XAL2ER(N,NN,N5,N4)
        TFE2ER(N2,N1)=TFE2ER(N2,N1)-XFE2ER(N,NN,N5,N4)
        TOH0ER(N2,N1)=TOH0ER(N2,N1)-XOH0ER(N,NN,N5,N4)
        TOH1ER(N2,N1)=TOH1ER(N2,N1)-XOH1ER(N,NN,N5,N4)
        TOH2ER(N2,N1)=TOH2ER(N2,N1)-XOH2ER(N,NN,N5,N4)
        TH1PER(N2,N1)=TH1PER(N2,N1)-XH1PER(N,NN,N5,N4)
        TH2PER(N2,N1)=TH2PER(N2,N1)-XH2PER(N,NN,N5,N4)
        TOH0EB(N2,N1)=TOH0EB(N2,N1)-XOH0EB(N,NN,N5,N4)
        TOH1EB(N2,N1)=TOH1EB(N2,N1)-XOH1EB(N,NN,N5,N4)
        TOH2EB(N2,N1)=TOH2EB(N2,N1)-XOH2EB(N,NN,N5,N4)
        TH1PEB(N2,N1)=TH1PEB(N2,N1)-XH1PEB(N,NN,N5,N4)
        TH2PEB(N2,N1)=TH2PEB(N2,N1)-XH2PEB(N,NN,N5,N4)
        TALOER(N2,N1)=TALOER(N2,N1)-PALOER(N,NN,N5,N4)
        TFEOER(N2,N1)=TFEOER(N2,N1)-PFEOER(N,NN,N5,N4)
        TCACER(N2,N1)=TCACER(N2,N1)-PCACER(N,NN,N5,N4)
        TCASER(N2,N1)=TCASER(N2,N1)-PCASER(N,NN,N5,N4)
        TALPER(N2,N1)=TALPER(N2,N1)-PALPER(N,NN,N5,N4)
        TFEPER(N2,N1)=TFEPER(N2,N1)-PFEPER(N,NN,N5,N4)
        TCPDER(N2,N1)=TCPDER(N2,N1)-PCPDER(N,NN,N5,N4)
        TCPHER(N2,N1)=TCPHER(N2,N1)-PCPHER(N,NN,N5,N4)
        TCPMER(N2,N1)=TCPMER(N2,N1)-PCPMER(N,NN,N5,N4)
        TALPEB(N2,N1)=TALPEB(N2,N1)-PALPEB(N,NN,N5,N4)
        TFEPEB(N2,N1)=TFEPEB(N2,N1)-PFEPEB(N,NN,N5,N4)
        TCPDEB(N2,N1)=TCPDEB(N2,N1)-PCPDEB(N,NN,N5,N4)
        TCPHEB(N2,N1)=TCPHEB(N2,N1)-PCPHEB(N,NN,N5,N4)
        TCPMEB(N2,N1)=TCPMEB(N2,N1)-PCPMEB(N,NN,N5,N4)

        DO  K=1,jcplx
          DO  NO=1,NFGs
            DO NGL=JGnio(NO),JGnfo(NO)
              DO  M=1,nlbiomcp
                TOMCER(M,NGL,K,N2,N1)=TOMCER(M,NGL,K,N2,N1)-OMCER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
                TOMNER(M,NGL,K,N2,N1)=TOMNER(M,NGL,K,N2,N1)-OMNER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
                TOMPER(M,NGL,K,N2,N1)=TOMPER(M,NGL,K,N2,N1)-OMPER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
              enddo
            enddo
          enddo
        ENDDO


        DO  NO=1,NFGs
          DO  M=1,nlbiomcp
            DO NGL=JGniA(NO),JGnfA(NO)
              TOMCERff(M,NGL,N2,N1)=TOMCERff(M,NGL,N2,N1)-OMCERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
              TOMNERff(M,NGL,N2,N1)=TOMNERff(M,NGL,N2,N1)-OMNERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
              TOMPERff(M,NGL,N2,N1)=TOMPERff(M,NGL,N2,N1)-OMPERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
            enddo
          enddo
        enddo


        D7375: DO K=1,jcplx
          D7370: DO M=1,ndbiomcp
            TORCER(M,K,N2,N1)=TORCER(M,K,N2,N1)-ORCER(M,K,N,NN,N5,N4)
            TORNER(M,K,N2,N1)=TORNER(M,K,N2,N1)-ORNER(M,K,N,NN,N5,N4)
            TORPER(M,K,N2,N1)=TORPER(M,K,N2,N1)-ORPER(M,K,N,NN,N5,N4)
          ENDDO D7370
          TOHCER(K,N2,N1)=TOHCER(K,N2,N1)-OHCER(K,N,NN,N5,N4)
          TOHNER(K,N2,N1)=TOHNER(K,N2,N1)-OHNER(K,N,NN,N5,N4)
          TOHPER(K,N2,N1)=TOHPER(K,N2,N1)-OHPER(K,N,NN,N5,N4)
          TOHAER(K,N2,N1)=TOHAER(K,N2,N1)-OHAER(K,N,NN,N5,N4)
          D7365: DO M=1,jsken
            TOSCER(M,K,N2,N1)=TOSCER(M,K,N2,N1)-OSCER(M,K,N,NN,N5,N4)
            TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)-OSAER(M,K,N,NN,N5,N4)
            TOSNER(M,K,N2,N1)=TOSNER(M,K,N2,N1)-OSNER(M,K,N,NN,N5,N4)
            TOSPER(M,K,N2,N1)=TOSPER(M,K,N2,N1)-OSPER(M,K,N,NN,N5,N4)
          ENDDO D7365
        ENDDO D7375
!     ENDIF
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        IF(ABS(XSEDER(N,NN,N5B,N4B)).GT.ZEROS(N5,N4))THEN
          TSEDER(N2,N1)=TSEDER(N2,N1)-XSEDER(N,NN,N5B,N4B)
          TSANER(N2,N1)=TSANER(N2,N1)-XSANER(N,NN,N5B,N4B)
          TSILER(N2,N1)=TSILER(N2,N1)-XSILER(N,NN,N5B,N4B)
          TCLAER(N2,N1)=TCLAER(N2,N1)-XCLAER(N,NN,N5B,N4B)
          TCECER(N2,N1)=TCECER(N2,N1)-XCECER(N,NN,N5B,N4B)
          TAECER(N2,N1)=TAECER(N2,N1)-XAECER(N,NN,N5B,N4B)
          TNH4ER(N2,N1)=TNH4ER(N2,N1)-XNH4ER(N,NN,N5B,N4B)
          TNH3ER(N2,N1)=TNH3ER(N2,N1)-XNH3ER(N,NN,N5B,N4B)
          TNHUER(N2,N1)=TNHUER(N2,N1)-XNHUER(N,NN,N5B,N4B)
          TNO3ER(N2,N1)=TNO3ER(N2,N1)-XNO3ER(N,NN,N5B,N4B)
          TNH4EB(N2,N1)=TNH4EB(N2,N1)-XNH4EB(N,NN,N5B,N4B)
          TNH3EB(N2,N1)=TNH3EB(N2,N1)-XNH3EB(N,NN,N5B,N4B)
          TNHUEB(N2,N1)=TNHUEB(N2,N1)-XNHUEB(N,NN,N5B,N4B)
          TNO3EB(N2,N1)=TNO3EB(N2,N1)-XNO3EB(N,NN,N5B,N4B)
          TN4ER(N2,N1)=TN4ER(N2,N1)-XN4ER(N,NN,N5B,N4B)
          TNBER(N2,N1)=TNBER(N2,N1)-XNBER(N,NN,N5B,N4B)
          THYER(N2,N1)=THYER(N2,N1)-XHYER(N,NN,N5B,N4B)
          TALER(N2,N1)=TALER(N2,N1)-XALER(N,NN,N5B,N4B)
          TFEER(N2,N1)=TFEER(N2,N1)-XFEER(N,NN,N5B,N4B)
          TCAER(N2,N1)=TCAER(N2,N1)-XCAER(N,NN,N5B,N4B)
          TMGER(N2,N1)=TMGER(N2,N1)-XMGER(N,NN,N5B,N4B)
          TNAER(N2,N1)=TNAER(N2,N1)-XNAER(N,NN,N5B,N4B)
          TKAER(N2,N1)=TKAER(N2,N1)-XKAER(N,NN,N5B,N4B)
          THCER(N2,N1)=THCER(N2,N1)-XHCER(N,NN,N5B,N4B)
          TAL2ER(N2,N1)=TAL2ER(N2,N1)-XAL2ER(N,NN,N5B,N4B)
          TFE2ER(N2,N1)=TFE2ER(N2,N1)-XFE2ER(N,NN,N5B,N4B)
          TOH0ER(N2,N1)=TOH0ER(N2,N1)-XOH0ER(N,NN,N5B,N4B)
          TOH1ER(N2,N1)=TOH1ER(N2,N1)-XOH1ER(N,NN,N5B,N4B)
          TOH2ER(N2,N1)=TOH2ER(N2,N1)-XOH2ER(N,NN,N5B,N4B)
          TH1PER(N2,N1)=TH1PER(N2,N1)-XH1PER(N,NN,N5B,N4B)
          TH2PER(N2,N1)=TH2PER(N2,N1)-XH2PER(N,NN,N5B,N4B)
          TOH0EB(N2,N1)=TOH0EB(N2,N1)-XOH0EB(N,NN,N5B,N4B)
          TOH1EB(N2,N1)=TOH1EB(N2,N1)-XOH1EB(N,NN,N5B,N4B)
          TOH2EB(N2,N1)=TOH2EB(N2,N1)-XOH2EB(N,NN,N5B,N4B)
          TH1PEB(N2,N1)=TH1PEB(N2,N1)-XH1PEB(N,NN,N5B,N4B)
          TH2PEB(N2,N1)=TH2PEB(N2,N1)-XH2PEB(N,NN,N5B,N4B)
          TALOER(N2,N1)=TALOER(N2,N1)-PALOER(N,NN,N5B,N4B)
          TFEOER(N2,N1)=TFEOER(N2,N1)-PFEOER(N,NN,N5B,N4B)
          TCACER(N2,N1)=TCACER(N2,N1)-PCACER(N,NN,N5B,N4B)
          TCASER(N2,N1)=TCASER(N2,N1)-PCASER(N,NN,N5B,N4B)
          TALPER(N2,N1)=TALPER(N2,N1)-PALPER(N,NN,N5B,N4B)
          TFEPER(N2,N1)=TFEPER(N2,N1)-PFEPER(N,NN,N5B,N4B)
          TCPDER(N2,N1)=TCPDER(N2,N1)-PCPDER(N,NN,N5B,N4B)
          TCPHER(N2,N1)=TCPHER(N2,N1)-PCPHER(N,NN,N5B,N4B)
          TCPMER(N2,N1)=TCPMER(N2,N1)-PCPMER(N,NN,N5B,N4B)
          TALPEB(N2,N1)=TALPEB(N2,N1)-PALPEB(N,NN,N5B,N4B)
          TFEPEB(N2,N1)=TFEPEB(N2,N1)-PFEPEB(N,NN,N5B,N4B)
          TCPDEB(N2,N1)=TCPDEB(N2,N1)-PCPDEB(N,NN,N5B,N4B)
          TCPHEB(N2,N1)=TCPHEB(N2,N1)-PCPHEB(N,NN,N5B,N4B)
          TCPMEB(N2,N1)=TCPMEB(N2,N1)-PCPMEB(N,NN,N5B,N4B)

          D8380: DO K=1,jcplx
            DO  NO=1,NFGs
              DO NGL=JGnio(NO),JGnfo(NO)
                DO  M=1,nlbiomcp
                  TOMCER(M,NGL,K,N2,N1)=TOMCER(M,NGL,K,N2,N1)-OMCER(M+(NGL-1)*nlbiomcp,K,N,NN,N5B,N4B)
                  TOMNER(M,NGL,K,N2,N1)=TOMNER(M,NGL,K,N2,N1)-OMNER(M+(NGL-1)*nlbiomcp,K,N,NN,N5B,N4B)
                  TOMPER(M,NGL,K,N2,N1)=TOMPER(M,NGL,K,N2,N1)-OMPER(M+(NGL-1)*nlbiomcp,K,N,NN,N5B,N4B)
                enddo
              enddo
            enddo
          ENDDO D8380

          DO  NO=1,NFGs
            DO NGL=JGniA(NO),JGnfA(NO)
              DO  M=1,nlbiomcp
                TOMCERff(M,NGL,N2,N1)=TOMCERff(M,NGL,N2,N1)-OMCERff(M+(NGL-1)*nlbiomcp,N,NN,N5B,N4B)
                TOMNERff(M,NGL,N2,N1)=TOMNERff(M,NGL,N2,N1)-OMNERff(M+(NGL-1)*nlbiomcp,N,NN,N5B,N4B)
                TOMPERff(M,NGL,N2,N1)=TOMPERff(M,NGL,N2,N1)-OMPERff(M+(NGL-1)*nlbiomcp,N,NN,N5B,N4B)
              enddo
            enddo
          enddo


          D8375: DO K=1,jcplx
            D8370: DO M=1,ndbiomcp
              TORCER(M,K,N2,N1)=TORCER(M,K,N2,N1)-ORCER(M,K,N,NN,N5B,N4B)
              TORNER(M,K,N2,N1)=TORNER(M,K,N2,N1)-ORNER(M,K,N,NN,N5B,N4B)
              TORPER(M,K,N2,N1)=TORPER(M,K,N2,N1)-ORPER(M,K,N,NN,N5B,N4B)
            ENDDO D8370
            TOHCER(K,N2,N1)=TOHCER(K,N2,N1)-OHCER(K,N,NN,N5B,N4B)
            TOHNER(K,N2,N1)=TOHNER(K,N2,N1)-OHNER(K,N,NN,N5B,N4B)
            TOHPER(K,N2,N1)=TOHPER(K,N2,N1)-OHPER(K,N,NN,N5B,N4B)
            TOHAER(K,N2,N1)=TOHAER(K,N2,N1)-OHAER(K,N,NN,N5B,N4B)
            D8365: DO M=1,jsken
              TOSCER(M,K,N2,N1)=TOSCER(M,K,N2,N1)-OSCER(M,K,N,NN,N5B,N4B)
              TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)-OSAER(M,K,N,NN,N5B,N4B)
              TOSNER(M,K,N2,N1)=TOSNER(M,K,N2,N1)-OSNER(M,K,N,NN,N5B,N4B)
              TOSPER(M,K,N2,N1)=TOSPER(M,K,N2,N1)-OSPER(M,K,N,NN,N5B,N4B)
            ENDDO D8365
          ENDDO D8375
        ENDIF
      ENDIF
    ENDDO D9350
  ENDIF
  end subroutine TotalFluxFromSedmentTransp
!------------------------------------------------------------------------------------------

  subroutine FluxBetweenGrids(N,N1,N2,N3,N4,N5,N6,NY,NX)
  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,NY,NX
  integer, intent(inout) :: N6
  integer :: LL,K
  !     begin_execution
  !     NET HEAT, WATER FLUXES BETWEEN ADJACENT
  !     GRID CELLS
  !
  !     TFLW,TFLWH,TFLWH=net micropore,macropore water flux, heat flux
  !     FLW,FLWH,HFLW=micropore,macropore water flux, heat flux from watsub.f
  !     FLWNU,FLWHNU,HFLWNU=lake surface water flux, heat flux from watsub.f if lake surface disappears
  !
  IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
    DO 1200 LL=N6,NL(N5,N4)
      IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
        N6=LL
        GO TO 1201
      ENDIF
1200  CONTINUE
1201  CONTINUE
    IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(N3.EQ.NU(N2,N1).AND.N.EQ.3)THEN
        TFLW(N3,N2,N1)=TFLW(N3,N2,N1)+FLW(N,N3,N2,N1)-FLWNU(N5,N4)
        TFLWX(N3,N2,N1)=TFLWX(N3,N2,N1)+FLWX(N,N3,N2,N1)-FLWXNU(N5,N4)
        TFLWH(N3,N2,N1)=TFLWH(N3,N2,N1)+FLWH(N,N3,N2,N1)-FLWHNU(N5,N4)
        THFLW(N3,N2,N1)=THFLW(N3,N2,N1)+HFLW(N,N3,N2,N1)-HFLWNU(N5,N4)
        if(THFLW(N3,N2,N1)<-1.e10)then
          write(*,*)'THFLW(N3,N2,N1)+HFLW(N,N3,N2,N1)-HFLWNU(N5,N4)',&
            THFLW(N3,N2,N1),HFLW(N,N3,N2,N1),HFLWNU(N5,N4)
          write(*,*)'Ns=',N1,n2,n3,n4,n5
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ELSE
        TFLW(N3,N2,N1)=TFLW(N3,N2,N1)+FLW(N,N3,N2,N1)-FLW(N,N6,N5,N4)
        TFLWX(N3,N2,N1)=TFLWX(N3,N2,N1)+FLWX(N,N3,N2,N1)-FLWX(N,N6,N5,N4)
        TFLWH(N3,N2,N1)=TFLWH(N3,N2,N1)+FLWH(N,N3,N2,N1)-FLWH(N,N6,N5,N4)
        THFLW(N3,N2,N1)=THFLW(N3,N2,N1)+HFLW(N,N3,N2,N1)-HFLW(N,N6,N5,N4)
        if(THFLW(N3,N2,N1)<-1.e10)then
          write(*,*)'THFLW(N3,N2,N1)+HFLW(N,N3,N2,N1)-HFLW(N,N6,N5,N4)',&
            THFLW(N3,N2,N1),HFLW(N,N3,N2,N1),HFLW(N,N6,N5,N4)
          write(*,*)'Ns=',N,N1,n2,n3,n4,n5,n6
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ENDIF
      !     IF(N1.EQ.1.AND.N3.EQ.1)THEN
      !     WRITE(*,6632)'TFLW',I,J,N,N1,N2,N3,N4,N5,N6,NU(N2,N1)
      !    2,TFLW(N3,N2,N1),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4),FLWNU(N5,N4)
      !    3,THFLW(N3,N2,N1),HFLW(N,N3,N2,N1),HFLW(N,N6,N5,N4)
      !    2,HFLWNU(N5,N4),VOLW(N3,N2,N1)
!6632  FORMAT(A8,10I4,12E16.8)
      !     ENDIF
      !
      !     NET SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
      !
      !     T*FLS,T*FHS=net convective+diffusive solute flux through micropores,macropores
      !     X*FLS,X*FHS=convective+diffusive solute flux through micropores, macropores from trnsfr.f
      !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
      !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
      !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
      !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
      !
      DO 8585 K=1,jcplx
        TOCFLS(K,N3,N2,N1)=TOCFLS(K,N3,N2,N1)+XOCFLS(K,N,N3,N2,N1)-XOCFLS(K,N,N6,N5,N4)
        TONFLS(K,N3,N2,N1)=TONFLS(K,N3,N2,N1)+XONFLS(K,N,N3,N2,N1)-XONFLS(K,N,N6,N5,N4)
        TOPFLS(K,N3,N2,N1)=TOPFLS(K,N3,N2,N1)+XOPFLS(K,N,N3,N2,N1)-XOPFLS(K,N,N6,N5,N4)
        TOAFLS(K,N3,N2,N1)=TOAFLS(K,N3,N2,N1)+XOAFLS(K,N,N3,N2,N1)-XOAFLS(K,N,N6,N5,N4)
        TOCFHS(K,N3,N2,N1)=TOCFHS(K,N3,N2,N1)+XOCFHS(K,N,N3,N2,N1)-XOCFHS(K,N,N6,N5,N4)
        TONFHS(K,N3,N2,N1)=TONFHS(K,N3,N2,N1)+XONFHS(K,N,N3,N2,N1)-XONFHS(K,N,N6,N5,N4)
        TOPFHS(K,N3,N2,N1)=TOPFHS(K,N3,N2,N1)+XOPFHS(K,N,N3,N2,N1)-XOPFHS(K,N,N6,N5,N4)
        TOAFHS(K,N3,N2,N1)=TOAFHS(K,N3,N2,N1)+XOAFHS(K,N,N3,N2,N1)-XOAFHS(K,N,N6,N5,N4)
8585  CONTINUE
      TCOFLS(N3,N2,N1)=TCOFLS(N3,N2,N1)+trcs_XFLS(idg_CO2,N,N3,N2,N1)-trcs_XFLS(idg_CO2,N,N6,N5,N4)
      TCHFLS(N3,N2,N1)=TCHFLS(N3,N2,N1)+trcs_XFLS(idg_CH4,N,N3,N2,N1)-trcs_XFLS(idg_CH4,N,N6,N5,N4)
      TOXFLS(N3,N2,N1)=TOXFLS(N3,N2,N1)+trcs_XFLS(idg_O2,N,N3,N2,N1)-trcs_XFLS(idg_O2,N,N6,N5,N4)
      TNGFLS(N3,N2,N1)=TNGFLS(N3,N2,N1)+trcs_XFLS(idg_N2,N,N3,N2,N1)-trcs_XFLS(idg_N2,N,N6,N5,N4)
      TN2FLS(N3,N2,N1)=TN2FLS(N3,N2,N1)+trcs_XFLS(idg_N2O,N,N3,N2,N1)-trcs_XFLS(idg_N2O,N,N6,N5,N4)
      THGFLS(N3,N2,N1)=THGFLS(N3,N2,N1)+trcs_XFLS(idg_H2,N,N3,N2,N1)-trcs_XFLS(idg_H2,N,N6,N5,N4)
      TN4FLS(N3,N2,N1)=TN4FLS(N3,N2,N1)+trcs_XFLS(ids_NH4,N,N3,N2,N1)-trcs_XFLS(ids_NH4,N,N6,N5,N4)
      TN3FLS(N3,N2,N1)=TN3FLS(N3,N2,N1)+trcs_XFLS(idg_NH3,N,N3,N2,N1)-trcs_XFLS(idg_NH3,N,N6,N5,N4)
      TNOFLS(N3,N2,N1)=TNOFLS(N3,N2,N1)+trcs_XFLS(ids_NO3,N,N3,N2,N1)-trcs_XFLS(ids_NO3,N,N6,N5,N4)
      TNXFLS(N3,N2,N1)=TNXFLS(N3,N2,N1)+trcs_XFLS(ids_NO2,N,N3,N2,N1)-trcs_XFLS(ids_NO2,N,N6,N5,N4)
      TP1FLS(N3,N2,N1)=TP1FLS(N3,N2,N1)+trcs_XFLS(ids_H1PO4,N,N3,N2,N1)-trcs_XFLS(ids_H1PO4,N,N6,N5,N4)
      TPOFLS(N3,N2,N1)=TPOFLS(N3,N2,N1)+trcs_XFLS(ids_H2PO4,N,N3,N2,N1)-trcs_XFLS(ids_H2PO4,N,N6,N5,N4)
      TN4FLB(N3,N2,N1)=TN4FLB(N3,N2,N1)+trcs_XFLS(ids_NH4B,N,N3,N2,N1)-trcs_XFLS(ids_NH4B,N,N6,N5,N4)
      TN3FLB(N3,N2,N1)=TN3FLB(N3,N2,N1)+trcs_XFLS(idg_NH3B,N,N3,N2,N1)-trcs_XFLS(idg_NH3B,N,N6,N5,N4)
      TNOFLB(N3,N2,N1)=TNOFLB(N3,N2,N1)+trcs_XFLS(ids_NO3B,N,N3,N2,N1)-trcs_XFLS(ids_NO3B,N,N6,N5,N4)
      TNXFLB(N3,N2,N1)=TNXFLB(N3,N2,N1)+trcs_XFLS(ids_NO2B,N,N3,N2,N1)-trcs_XFLS(ids_NO2B,N,N6,N5,N4)
      TH1BFB(N3,N2,N1)=TH1BFB(N3,N2,N1)+trcs_XFLS(ids_H1PO4B,N,N3,N2,N1)-trcs_XFLS(ids_H1PO4B,N,N6,N5,N4)
      TH2BFB(N3,N2,N1)=TH2BFB(N3,N2,N1)+trcs_XFLS(ids_H2PO4B,N,N3,N2,N1)-trcs_XFLS(ids_H2PO4B,N,N6,N5,N4)
      TCOFHS(N3,N2,N1)=TCOFHS(N3,N2,N1)+XCOFHS(N,N3,N2,N1)-XCOFHS(N,N6,N5,N4)
      TCHFHS(N3,N2,N1)=TCHFHS(N3,N2,N1)+XCHFHS(N,N3,N2,N1)-XCHFHS(N,N6,N5,N4)
      TOXFHS(N3,N2,N1)=TOXFHS(N3,N2,N1)+XOXFHS(N,N3,N2,N1)-XOXFHS(N,N6,N5,N4)
      TNGFHS(N3,N2,N1)=TNGFHS(N3,N2,N1)+XNGFHS(N,N3,N2,N1)-XNGFHS(N,N6,N5,N4)
      TN2FHS(N3,N2,N1)=TN2FHS(N3,N2,N1)+XN2FHS(N,N3,N2,N1)-XN2FHS(N,N6,N5,N4)
      THGFHS(N3,N2,N1)=THGFHS(N3,N2,N1)+XHGFHS(N,N3,N2,N1)-XHGFHS(N,N6,N5,N4)
      TN4FHS(N3,N2,N1)=TN4FHS(N3,N2,N1)+XN4FHW(N,N3,N2,N1)-XN4FHW(N,N6,N5,N4)
      TN3FHS(N3,N2,N1)=TN3FHS(N3,N2,N1)+XN3FHW(N,N3,N2,N1)-XN3FHW(N,N6,N5,N4)
      TNOFHS(N3,N2,N1)=TNOFHS(N3,N2,N1)+XNOFHW(N,N3,N2,N1)-XNOFHW(N,N6,N5,N4)
      TNXFHS(N3,N2,N1)=TNXFHS(N3,N2,N1)+XNXFHS(N,N3,N2,N1)-XNXFHS(N,N6,N5,N4)
      TP1FHS(N3,N2,N1)=TP1FHS(N3,N2,N1)+XH1PHS(N,N3,N2,N1)-XH1PHS(N,N6,N5,N4)
      TPOFHS(N3,N2,N1)=TPOFHS(N3,N2,N1)+XH2PHS(N,N3,N2,N1)-XH2PHS(N,N6,N5,N4)
      TN4FHB(N3,N2,N1)=TN4FHB(N3,N2,N1)+XN4FHB(N,N3,N2,N1)-XN4FHB(N,N6,N5,N4)
      TN3FHB(N3,N2,N1)=TN3FHB(N3,N2,N1)+XN3FHB(N,N3,N2,N1)-XN3FHB(N,N6,N5,N4)
      TNOFHB(N3,N2,N1)=TNOFHB(N3,N2,N1)+XNOFHB(N,N3,N2,N1)-XNOFHB(N,N6,N5,N4)
      TNXFHB(N3,N2,N1)=TNXFHB(N3,N2,N1)+XNXFHB(N,N3,N2,N1)-XNXFHB(N,N6,N5,N4)
      TH1BHB(N3,N2,N1)=TH1BHB(N3,N2,N1)+XH1BHB(N,N3,N2,N1)-XH1BHB(N,N6,N5,N4)
      TH2BHB(N3,N2,N1)=TH2BHB(N3,N2,N1)+XH2BHB(N,N3,N2,N1)-XH2BHB(N,N6,N5,N4)
!
      !     NET GAS FLUXES BETWEEN ADJACENT GRID CELLS
      !
      !     T*FLG=net convective+diffusive gas flux
      !     X*FLG=convective+diffusive gas flux from trnsfr.f
      !     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!
      TCOFLG(N3,N2,N1)=TCOFLG(N3,N2,N1)+R3GasADTFlx(idg_CO2,N,N3,N2,N1)-R3GasADTFlx(idg_CO2,N,N6,N5,N4)
      TCHFLG(N3,N2,N1)=TCHFLG(N3,N2,N1)+R3GasADTFlx(idg_CH4,N,N3,N2,N1)-R3GasADTFlx(idg_CH4,N,N6,N5,N4)
      TOXFLG(N3,N2,N1)=TOXFLG(N3,N2,N1)+R3GasADTFlx(idg_O2,N,N3,N2,N1)-R3GasADTFlx(idg_O2,N,N6,N5,N4)
      TNGFLG(N3,N2,N1)=TNGFLG(N3,N2,N1)+R3GasADTFlx(idg_N2,N,N3,N2,N1)-R3GasADTFlx(idg_N2,N,N6,N5,N4)
      TN2FLG(N3,N2,N1)=TN2FLG(N3,N2,N1)+R3GasADTFlx(idg_N2O,N,N3,N2,N1)-R3GasADTFlx(idg_N2O,N,N6,N5,N4)
      TNHFLG(N3,N2,N1)=TNHFLG(N3,N2,N1)+R3GasADTFlx(idg_NH3,N,N3,N2,N1)-R3GasADTFlx(idg_NH3,N,N6,N5,N4)
      THGFLG(N3,N2,N1)=THGFLG(N3,N2,N1)+R3GasADTFlx(idg_H2,N,N3,N2,N1)-R3GasADTFlx(idg_H2,N,N6,N5,N4)
!
      !     NET SALT FLUXES BETWEEN ADJACENT GRID CELLS
      !
      !     T*FLS,T*FHS=net convective+diffusive solute flux through micropores,macropores
      !     X*FLS,X*FHS=convective+diffusive solute flux through micropores, macropores from trnsfrs.f
      !     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
      !          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
      !          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
      !          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
      !          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
      !          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
      !          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
      !     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
      !          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
      !          :*1=non-band,*B=band
!
      IF(ISALTG.NE.0)THEN
        TALFLS(N3,N2,N1)=TALFLS(N3,N2,N1)+XALFLS(N,N3,N2,N1)-XALFLS(N,N6,N5,N4)
        TFEFLS(N3,N2,N1)=TFEFLS(N3,N2,N1)+XFEFLS(N,N3,N2,N1)-XFEFLS(N,N6,N5,N4)
        THYFLS(N3,N2,N1)=THYFLS(N3,N2,N1)+XHYFLS(N,N3,N2,N1)-XHYFLS(N,N6,N5,N4)
        TCAFLS(N3,N2,N1)=TCAFLS(N3,N2,N1)+XCAFLS(N,N3,N2,N1)-XCAFLS(N,N6,N5,N4)
        TMGFLS(N3,N2,N1)=TMGFLS(N3,N2,N1)+XMGFLS(N,N3,N2,N1)-XMGFLS(N,N6,N5,N4)
        TNAFLS(N3,N2,N1)=TNAFLS(N3,N2,N1)+XNAFLS(N,N3,N2,N1)-XNAFLS(N,N6,N5,N4)
        TKAFLS(N3,N2,N1)=TKAFLS(N3,N2,N1)+XKAFLS(N,N3,N2,N1)-XKAFLS(N,N6,N5,N4)
        TOHFLS(N3,N2,N1)=TOHFLS(N3,N2,N1)+XOHFLS(N,N3,N2,N1)-XOHFLS(N,N6,N5,N4)
        TSOFLS(N3,N2,N1)=TSOFLS(N3,N2,N1)+XSOFLS(N,N3,N2,N1)-XSOFLS(N,N6,N5,N4)
        TCLFLS(N3,N2,N1)=TCLFLS(N3,N2,N1)+XCLFLS(N,N3,N2,N1)-XCLFLS(N,N6,N5,N4)
        TC3FLS(N3,N2,N1)=TC3FLS(N3,N2,N1)+XC3FLS(N,N3,N2,N1)-XC3FLS(N,N6,N5,N4)
        THCFLS(N3,N2,N1)=THCFLS(N3,N2,N1)+XHCFLS(N,N3,N2,N1)-XHCFLS(N,N6,N5,N4)
        TAL1FS(N3,N2,N1)=TAL1FS(N3,N2,N1)+XAL1FS(N,N3,N2,N1)-XAL1FS(N,N6,N5,N4)
        TAL2FS(N3,N2,N1)=TAL2FS(N3,N2,N1)+XAL2FS(N,N3,N2,N1)-XAL2FS(N,N6,N5,N4)
        TAL3FS(N3,N2,N1)=TAL3FS(N3,N2,N1)+XAL3FS(N,N3,N2,N1)-XAL3FS(N,N6,N5,N4)
        TAL4FS(N3,N2,N1)=TAL4FS(N3,N2,N1)+XAL4FS(N,N3,N2,N1)-XAL4FS(N,N6,N5,N4)
        TALSFS(N3,N2,N1)=TALSFS(N3,N2,N1)+XALSFS(N,N3,N2,N1)-XALSFS(N,N6,N5,N4)
        TFE1FS(N3,N2,N1)=TFE1FS(N3,N2,N1)+XFE1FS(N,N3,N2,N1)-XFE1FS(N,N6,N5,N4)
        TFE2FS(N3,N2,N1)=TFE2FS(N3,N2,N1)+XFE2FS(N,N3,N2,N1)-XFE2FS(N,N6,N5,N4)
        TFE3FS(N3,N2,N1)=TFE3FS(N3,N2,N1)+XFE3FS(N,N3,N2,N1)-XFE3FS(N,N6,N5,N4)
        TFE4FS(N3,N2,N1)=TFE4FS(N3,N2,N1)+XFE4FS(N,N3,N2,N1)-XFE4FS(N,N6,N5,N4)
        TFESFS(N3,N2,N1)=TFESFS(N3,N2,N1)+XFESFS(N,N3,N2,N1)-XFESFS(N,N6,N5,N4)
        TCAOFS(N3,N2,N1)=TCAOFS(N3,N2,N1)+XCAOFS(N,N3,N2,N1)-XCAOFS(N,N6,N5,N4)
        TCACFS(N3,N2,N1)=TCACFS(N3,N2,N1)+XCACFS(N,N3,N2,N1)-XCACFS(N,N6,N5,N4)
        TCAHFS(N3,N2,N1)=TCAHFS(N3,N2,N1)+XCAHFS(N,N3,N2,N1)-XCAHFS(N,N6,N5,N4)
        TCASFS(N3,N2,N1)=TCASFS(N3,N2,N1)+XCASFS(N,N3,N2,N1)-XCASFS(N,N6,N5,N4)
        TMGOFS(N3,N2,N1)=TMGOFS(N3,N2,N1)+XMGOFS(N,N3,N2,N1)-XMGOFS(N,N6,N5,N4)
        TMGCFS(N3,N2,N1)=TMGCFS(N3,N2,N1)+XMGCFS(N,N3,N2,N1)-XMGCFS(N,N6,N5,N4)
        TMGHFS(N3,N2,N1)=TMGHFS(N3,N2,N1)+XMGHFS(N,N3,N2,N1)-XMGHFS(N,N6,N5,N4)
        TMGSFS(N3,N2,N1)=TMGSFS(N3,N2,N1)+XMGSFS(N,N3,N2,N1)-XMGSFS(N,N6,N5,N4)
        TNACFS(N3,N2,N1)=TNACFS(N3,N2,N1)+XNACFS(N,N3,N2,N1)-XNACFS(N,N6,N5,N4)
        TNASFS(N3,N2,N1)=TNASFS(N3,N2,N1)+XNASFS(N,N3,N2,N1)-XNASFS(N,N6,N5,N4)
        TKASFS(N3,N2,N1)=TKASFS(N3,N2,N1)+XKASFS(N,N3,N2,N1)-XKASFS(N,N6,N5,N4)
        TH0PFS(N3,N2,N1)=TH0PFS(N3,N2,N1)+XH0PFS(N,N3,N2,N1)-XH0PFS(N,N6,N5,N4)
        TH3PFS(N3,N2,N1)=TH3PFS(N3,N2,N1)+XH3PFS(N,N3,N2,N1)-XH3PFS(N,N6,N5,N4)
        TF1PFS(N3,N2,N1)=TF1PFS(N3,N2,N1)+XF1PFS(N,N3,N2,N1)-XF1PFS(N,N6,N5,N4)
        TF2PFS(N3,N2,N1)=TF2PFS(N3,N2,N1)+XF2PFS(N,N3,N2,N1)-XF2PFS(N,N6,N5,N4)
        TC0PFS(N3,N2,N1)=TC0PFS(N3,N2,N1)+XC0PFS(N,N3,N2,N1)-XC0PFS(N,N6,N5,N4)
        TC1PFS(N3,N2,N1)=TC1PFS(N3,N2,N1)+XC1PFS(N,N3,N2,N1)-XC1PFS(N,N6,N5,N4)
        TC2PFS(N3,N2,N1)=TC2PFS(N3,N2,N1)+XC2PFS(N,N3,N2,N1)-XC2PFS(N,N6,N5,N4)
        TM1PFS(N3,N2,N1)=TM1PFS(N3,N2,N1)+XM1PFS(N,N3,N2,N1)-XM1PFS(N,N6,N5,N4)
        TH0BFB(N3,N2,N1)=TH0BFB(N3,N2,N1)+XH0BFB(N,N3,N2,N1)-XH0BFB(N,N6,N5,N4)
        TH3BFB(N3,N2,N1)=TH3BFB(N3,N2,N1)+XH3BFB(N,N3,N2,N1)-XH3BFB(N,N6,N5,N4)
        TF1BFB(N3,N2,N1)=TF1BFB(N3,N2,N1)+XF1BFB(N,N3,N2,N1)-XF1BFB(N,N6,N5,N4)
        TF2BFB(N3,N2,N1)=TF2BFB(N3,N2,N1)+XF2BFB(N,N3,N2,N1)-XF2BFB(N,N6,N5,N4)
        TC0BFB(N3,N2,N1)=TC0BFB(N3,N2,N1)+XC0BFB(N,N3,N2,N1)-XC0BFB(N,N6,N5,N4)
        TC1BFB(N3,N2,N1)=TC1BFB(N3,N2,N1)+XC1BFB(N,N3,N2,N1)-XC1BFB(N,N6,N5,N4)
        TC2BFB(N3,N2,N1)=TC2BFB(N3,N2,N1)+XC2BFB(N,N3,N2,N1)-XC2BFB(N,N6,N5,N4)
        TM1BFB(N3,N2,N1)=TM1BFB(N3,N2,N1)+XM1BFB(N,N3,N2,N1)-XM1BFB(N,N6,N5,N4)
        TALFHS(N3,N2,N1)=TALFHS(N3,N2,N1)+XALFHS(N,N3,N2,N1)-XALFHS(N,N6,N5,N4)
        TFEFHS(N3,N2,N1)=TFEFHS(N3,N2,N1)+XFEFHS(N,N3,N2,N1)-XFEFHS(N,N6,N5,N4)
        THYFHS(N3,N2,N1)=THYFHS(N3,N2,N1)+XHYFHS(N,N3,N2,N1)-XHYFHS(N,N6,N5,N4)
        TCAFHS(N3,N2,N1)=TCAFHS(N3,N2,N1)+XCAFHS(N,N3,N2,N1)-XCAFHS(N,N6,N5,N4)
        TMGFHS(N3,N2,N1)=TMGFHS(N3,N2,N1)+XMGFHS(N,N3,N2,N1)-XMGFHS(N,N6,N5,N4)
        TNAFHS(N3,N2,N1)=TNAFHS(N3,N2,N1)+XNAFHS(N,N3,N2,N1)-XNAFHS(N,N6,N5,N4)
        TKAFHS(N3,N2,N1)=TKAFHS(N3,N2,N1)+XKAFHS(N,N3,N2,N1)-XKAFHS(N,N6,N5,N4)
        TOHFHS(N3,N2,N1)=TOHFHS(N3,N2,N1)+XOHFHS(N,N3,N2,N1)-XOHFHS(N,N6,N5,N4)
        TSOFHS(N3,N2,N1)=TSOFHS(N3,N2,N1)+XSOFHS(N,N3,N2,N1)-XSOFHS(N,N6,N5,N4)
        TCLFHS(N3,N2,N1)=TCLFHS(N3,N2,N1)+XCLFHS(N,N3,N2,N1)-XCLFHS(N,N6,N5,N4)
        TC3FHS(N3,N2,N1)=TC3FHS(N3,N2,N1)+XC3FHS(N,N3,N2,N1)-XC3FHS(N,N6,N5,N4)
        THCFHS(N3,N2,N1)=THCFHS(N3,N2,N1)+XHCFHS(N,N3,N2,N1)-XHCFHS(N,N6,N5,N4)
        TAL1HS(N3,N2,N1)=TAL1HS(N3,N2,N1)+XAL1HS(N,N3,N2,N1)-XAL1HS(N,N6,N5,N4)
        TAL2HS(N3,N2,N1)=TAL2HS(N3,N2,N1)+XAL2HS(N,N3,N2,N1)-XAL2HS(N,N6,N5,N4)
        TAL3HS(N3,N2,N1)=TAL3HS(N3,N2,N1)+XAL3HS(N,N3,N2,N1)-XAL3HS(N,N6,N5,N4)
        TAL4HS(N3,N2,N1)=TAL4HS(N3,N2,N1)+XAL4HS(N,N3,N2,N1)-XAL4HS(N,N6,N5,N4)
        TALSHS(N3,N2,N1)=TALSHS(N3,N2,N1)+XALSHS(N,N3,N2,N1)-XALSHS(N,N6,N5,N4)
        TFE1HS(N3,N2,N1)=TFE1HS(N3,N2,N1)+XFE1HS(N,N3,N2,N1)-XFE1HS(N,N6,N5,N4)
        TFE2HS(N3,N2,N1)=TFE2HS(N3,N2,N1)+XFE2HS(N,N3,N2,N1)-XFE2HS(N,N6,N5,N4)
        TFE3HS(N3,N2,N1)=TFE3HS(N3,N2,N1)+XFE3HS(N,N3,N2,N1)-XFE3HS(N,N6,N5,N4)
        TFE4HS(N3,N2,N1)=TFE4HS(N3,N2,N1)+XFE4HS(N,N3,N2,N1)-XFE4HS(N,N6,N5,N4)
        TFESHS(N3,N2,N1)=TFESHS(N3,N2,N1)+XFESHS(N,N3,N2,N1)-XFESHS(N,N6,N5,N4)
        TCAOHS(N3,N2,N1)=TCAOHS(N3,N2,N1)+XCAOHS(N,N3,N2,N1)-XCAOHS(N,N6,N5,N4)
        TCACHS(N3,N2,N1)=TCACHS(N3,N2,N1)+XCACHS(N,N3,N2,N1)-XCACHS(N,N6,N5,N4)
        TCAHHS(N3,N2,N1)=TCAHHS(N3,N2,N1)+XCAHHS(N,N3,N2,N1)-XCAHHS(N,N6,N5,N4)
        TCASHS(N3,N2,N1)=TCASHS(N3,N2,N1)+XCASHS(N,N3,N2,N1)-XCASHS(N,N6,N5,N4)
        TMGOHS(N3,N2,N1)=TMGOHS(N3,N2,N1)+XMGOHS(N,N3,N2,N1)-XMGOHS(N,N6,N5,N4)
        TMGCHS(N3,N2,N1)=TMGCHS(N3,N2,N1)+XMGCHS(N,N3,N2,N1)-XMGCHS(N,N6,N5,N4)
        TMGHHS(N3,N2,N1)=TMGHHS(N3,N2,N1)+XMGHHS(N,N3,N2,N1)-XMGHHS(N,N6,N5,N4)
        TMGSHS(N3,N2,N1)=TMGSHS(N3,N2,N1)+XMGSHS(N,N3,N2,N1)-XMGSHS(N,N6,N5,N4)
        TNACHS(N3,N2,N1)=TNACHS(N3,N2,N1)+XNACHS(N,N3,N2,N1)-XNACHS(N,N6,N5,N4)
        TNASHS(N3,N2,N1)=TNASHS(N3,N2,N1)+XNASHS(N,N3,N2,N1)-XNASHS(N,N6,N5,N4)
        TKASHS(N3,N2,N1)=TKASHS(N3,N2,N1)+XKASHS(N,N3,N2,N1)-XKASHS(N,N6,N5,N4)
        TH0PHS(N3,N2,N1)=TH0PHS(N3,N2,N1)+XH0PHS(N,N3,N2,N1)-XH0PHS(N,N6,N5,N4)
        TH3PHS(N3,N2,N1)=TH3PHS(N3,N2,N1)+XH3PHS(N,N3,N2,N1)-XH3PHS(N,N6,N5,N4)
        TF1PHS(N3,N2,N1)=TF1PHS(N3,N2,N1)+XF1PHS(N,N3,N2,N1)-XF1PHS(N,N6,N5,N4)
        TF2PHS(N3,N2,N1)=TF2PHS(N3,N2,N1)+XF2PHS(N,N3,N2,N1)-XF2PHS(N,N6,N5,N4)
        TC0PHS(N3,N2,N1)=TC0PHS(N3,N2,N1)+XC0PHS(N,N3,N2,N1)-XC0PHS(N,N6,N5,N4)
        TC1PHS(N3,N2,N1)=TC1PHS(N3,N2,N1)+XC1PHS(N,N3,N2,N1)-XC1PHS(N,N6,N5,N4)
        TC2PHS(N3,N2,N1)=TC2PHS(N3,N2,N1)+XC2PHS(N,N3,N2,N1)-XC2PHS(N,N6,N5,N4)
        TM1PHS(N3,N2,N1)=TM1PHS(N3,N2,N1)+XM1PHS(N,N3,N2,N1)-XM1PHS(N,N6,N5,N4)
        TH0BHB(N3,N2,N1)=TH0BHB(N3,N2,N1)+XH0BHB(N,N3,N2,N1)-XH0BHB(N,N6,N5,N4)
        TH3BHB(N3,N2,N1)=TH3BHB(N3,N2,N1)+XH3BHB(N,N3,N2,N1)-XH3BHB(N,N6,N5,N4)
        TF1BHB(N3,N2,N1)=TF1BHB(N3,N2,N1)+XF1BHB(N,N3,N2,N1)-XF1BHB(N,N6,N5,N4)
        TF2BHB(N3,N2,N1)=TF2BHB(N3,N2,N1)+XF2BHB(N,N3,N2,N1)-XF2BHB(N,N6,N5,N4)
        TC0BHB(N3,N2,N1)=TC0BHB(N3,N2,N1)+XC0BHB(N,N3,N2,N1)-XC0BHB(N,N6,N5,N4)
        TC1BHB(N3,N2,N1)=TC1BHB(N3,N2,N1)+XC1BHB(N,N3,N2,N1)-XC1BHB(N,N6,N5,N4)
        TC2BHB(N3,N2,N1)=TC2BHB(N3,N2,N1)+XC2BHB(N,N3,N2,N1)-XC2BHB(N,N6,N5,N4)
        TM1BHB(N3,N2,N1)=TM1BHB(N3,N2,N1)+XM1BHB(N,N3,N2,N1)-XM1BHB(N,N6,N5,N4)
      ENDIF
    ELSE
      TFLW(N3,N2,N1)=0.0_r8
      TFLWX(N3,N2,N1)=0.0_r8
      TFLWH(N3,N2,N1)=0.0_r8
      THFLW(N3,N2,N1)=0.0_r8
      TTHAW(N3,N2,N1)=0.0_r8
      TTHAWH(N3,N2,N1)=0.0_r8
      THTHAW(N3,N2,N1)=0.0_r8
      DO 8596 K=1,jcplx
        TOCFLS(K,N3,N2,N1)=0.0_r8
        TONFLS(K,N3,N2,N1)=0.0_r8
        TOPFLS(K,N3,N2,N1)=0.0_r8
        TOAFLS(K,N3,N2,N1)=0.0_r8
        TOCFHS(K,N3,N2,N1)=0.0_r8
        TONFHS(K,N3,N2,N1)=0.0_r8
        TOPFHS(K,N3,N2,N1)=0.0_r8
        TOAFHS(K,N3,N2,N1)=0.0_r8
8596  CONTINUE
      TCOFLS(N3,N2,N1)=0.0_r8
      TCHFLS(N3,N2,N1)=0.0_r8
      TOXFLS(N3,N2,N1)=0.0_r8
      TNGFLS(N3,N2,N1)=0.0_r8
      TN2FLS(N3,N2,N1)=0.0_r8
      THGFLS(N3,N2,N1)=0.0_r8
      TN4FLS(N3,N2,N1)=0.0_r8
      TN3FLS(N3,N2,N1)=0.0_r8
      TNOFLS(N3,N2,N1)=0.0_r8
      TNXFLS(N3,N2,N1)=0.0_r8
      TP1FLS(N3,N2,N1)=0.0_r8
      TPOFLS(N3,N2,N1)=0.0_r8
      TN4FLB(N3,N2,N1)=0.0_r8
      TN3FLB(N3,N2,N1)=0.0_r8
      TNOFLB(N3,N2,N1)=0.0_r8
      TNXFLB(N3,N2,N1)=0.0_r8
      TH1BFB(N3,N2,N1)=0.0_r8
      TH2BFB(N3,N2,N1)=0.0_r8
      TCOFHS(N3,N2,N1)=0.0_r8
      TCHFHS(N3,N2,N1)=0.0_r8
      TOXFHS(N3,N2,N1)=0.0_r8
      TNGFHS(N3,N2,N1)=0.0_r8
      TN2FHS(N3,N2,N1)=0.0_r8
      THGFHS(N3,N2,N1)=0.0_r8
      TN4FHS(N3,N2,N1)=0.0_r8
      TN3FHS(N3,N2,N1)=0.0_r8
      TNOFHS(N3,N2,N1)=0.0_r8
      TNXFHS(N3,N2,N1)=0.0_r8
      TP1FHS(N3,N2,N1)=0.0_r8
      TPOFHS(N3,N2,N1)=0.0_r8
      TN4FHB(N3,N2,N1)=0.0_r8
      TN3FHB(N3,N2,N1)=0.0_r8
      TNOFHB(N3,N2,N1)=0.0_r8
      TNXFHB(N3,N2,N1)=0.0_r8
      TH1BHB(N3,N2,N1)=0.0_r8
      TH2BHB(N3,N2,N1)=0.0_r8
      TCOFLG(N3,N2,N1)=0.0_r8
      TCHFLG(N3,N2,N1)=0.0_r8
      TOXFLG(N3,N2,N1)=0.0_r8
      TNGFLG(N3,N2,N1)=0.0_r8
      TN2FLG(N3,N2,N1)=0.0_r8
      TNHFLG(N3,N2,N1)=0.0_r8
      THGFLG(N3,N2,N1)=0.0_r8
      IF(ISALTG.NE.0)THEN
        TALFLS(N3,N2,N1)=0.0_r8
        TFEFLS(N3,N2,N1)=0.0_r8
        THYFLS(N3,N2,N1)=0.0_r8
        TCAFLS(N3,N2,N1)=0.0_r8
        TMGFLS(N3,N2,N1)=0.0_r8
        TNAFLS(N3,N2,N1)=0.0_r8
        TKAFLS(N3,N2,N1)=0.0_r8
        TOHFLS(N3,N2,N1)=0.0_r8
        TSOFLS(N3,N2,N1)=0.0_r8
        TCLFLS(N3,N2,N1)=0.0_r8
        TC3FLS(N3,N2,N1)=0.0_r8
        THCFLS(N3,N2,N1)=0.0_r8
        TAL1FS(N3,N2,N1)=0.0_r8
        TAL2FS(N3,N2,N1)=0.0_r8
        TAL3FS(N3,N2,N1)=0.0_r8
        TAL4FS(N3,N2,N1)=0.0_r8
        TALSFS(N3,N2,N1)=0.0_r8
        TFE1FS(N3,N2,N1)=0.0_r8
        TFE2FS(N3,N2,N1)=0.0_r8
        TFE3FS(N3,N2,N1)=0.0_r8
        TFE4FS(N3,N2,N1)=0.0_r8
        TFESFS(N3,N2,N1)=0.0_r8
        TCAOFS(N3,N2,N1)=0.0_r8
        TCACFS(N3,N2,N1)=0.0_r8
        TCAHFS(N3,N2,N1)=0.0_r8
        TCASFS(N3,N2,N1)=0.0_r8
        TMGOFS(N3,N2,N1)=0.0_r8
        TMGCFS(N3,N2,N1)=0.0_r8
        TMGHFS(N3,N2,N1)=0.0_r8
        TMGSFS(N3,N2,N1)=0.0_r8
        TNACFS(N3,N2,N1)=0.0_r8
        TNASFS(N3,N2,N1)=0.0_r8
        TKASFS(N3,N2,N1)=0.0_r8
        TH0PFS(N3,N2,N1)=0.0_r8
        TH3PFS(N3,N2,N1)=0.0_r8
        TF1PFS(N3,N2,N1)=0.0_r8
        TF2PFS(N3,N2,N1)=0.0_r8
        TC0PFS(N3,N2,N1)=0.0_r8
        TC1PFS(N3,N2,N1)=0.0_r8
        TC2PFS(N3,N2,N1)=0.0_r8
        TM1PFS(N3,N2,N1)=0.0_r8
        TH0BFB(N3,N2,N1)=0.0_r8
        TH3BFB(N3,N2,N1)=0.0_r8
        TF1BFB(N3,N2,N1)=0.0_r8
        TF2BFB(N3,N2,N1)=0.0_r8
        TC0BFB(N3,N2,N1)=0.0_r8
        TC1BFB(N3,N2,N1)=0.0_r8
        TC2BFB(N3,N2,N1)=0.0_r8
        TM1BFB(N3,N2,N1)=0.0_r8
        TALFHS(N3,N2,N1)=0.0_r8
        TFEFHS(N3,N2,N1)=0.0_r8
        THYFHS(N3,N2,N1)=0.0_r8
        TCAFHS(N3,N2,N1)=0.0_r8
        TMGFHS(N3,N2,N1)=0.0_r8
        TNAFHS(N3,N2,N1)=0.0_r8
        TKAFHS(N3,N2,N1)=0.0_r8
        TOHFHS(N3,N2,N1)=0.0_r8
        TSOFHS(N3,N2,N1)=0.0_r8
        TCLFHS(N3,N2,N1)=0.0_r8
        TC3FHS(N3,N2,N1)=0.0_r8
        THCFHS(N3,N2,N1)=0.0_r8
        TAL1HS(N3,N2,N1)=0.0_r8
        TAL2HS(N3,N2,N1)=0.0_r8
        TAL3HS(N3,N2,N1)=0.0_r8
        TAL4HS(N3,N2,N1)=0.0_r8
        TALSHS(N3,N2,N1)=0.0_r8
        TFE1HS(N3,N2,N1)=0.0_r8
        TFE2HS(N3,N2,N1)=0.0_r8
        TFE3HS(N3,N2,N1)=0.0_r8
        TFE4HS(N3,N2,N1)=0.0_r8
        TFESHS(N3,N2,N1)=0.0_r8
        TCAOHS(N3,N2,N1)=0.0_r8
        TCACHS(N3,N2,N1)=0.0_r8
        TCAHHS(N3,N2,N1)=0.0_r8
        TCASHS(N3,N2,N1)=0.0_r8
        TMGOHS(N3,N2,N1)=0.0_r8
        TMGCHS(N3,N2,N1)=0.0_r8
        TMGHHS(N3,N2,N1)=0.0_r8
        TMGSHS(N3,N2,N1)=0.0_r8
        TNACHS(N3,N2,N1)=0.0_r8
        TNASHS(N3,N2,N1)=0.0_r8
        TKASHS(N3,N2,N1)=0.0_r8
        TH0PHS(N3,N2,N1)=0.0_r8
        TH3PHS(N3,N2,N1)=0.0_r8
        TF1PHS(N3,N2,N1)=0.0_r8
        TF2PHS(N3,N2,N1)=0.0_r8
        TC0PHS(N3,N2,N1)=0.0_r8
        TC1PHS(N3,N2,N1)=0.0_r8
        TC2PHS(N3,N2,N1)=0.0_r8
        TM1PHS(N3,N2,N1)=0.0_r8
        TH0BHB(N3,N2,N1)=0.0_r8
        TH3BHB(N3,N2,N1)=0.0_r8
        TF1BHB(N3,N2,N1)=0.0_r8
        TF2BHB(N3,N2,N1)=0.0_r8
        TC0BHB(N3,N2,N1)=0.0_r8
        TC1BHB(N3,N2,N1)=0.0_r8
        TC2BHB(N3,N2,N1)=0.0_r8
        TM1BHB(N3,N2,N1)=0.0_r8
      ENDIF
    ENDIF
  ENDIF
  end subroutine FluxBetweenGrids

end module LateralTranspMod
