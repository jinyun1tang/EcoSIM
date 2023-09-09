module LateralTranspMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
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
  use ErosionBalMod
  use SnowBalanceMod
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
  real(r8) :: VTATM,VTGAS
  real(r8) :: trcg_VOLG(idg_beg:idg_end)
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
    trcg_VOLG(idg_CO2)=trc_gasml(idg_CO2,L,NY,NX)/catomw
    trcg_VOLG(idg_CH4)=trc_gasml(idg_CH4,L,NY,NX)/catomw
    trcg_VOLG(idg_O2)=trc_gasml(idg_O2,L,NY,NX)/32.0_r8
    trcg_VOLG(idg_N2)=trc_gasml(idg_N2,L,NY,NX)/28.0_r8
    trcg_VOLG(idg_N2O)=trc_gasml(idg_N2O,L,NY,NX)/28.0_r8
    trcg_VOLG(idg_NH3)=trc_gasml(idg_NH3,L,NY,NX)/natomw
    trcg_VOLG(idg_H2)=trc_gasml(idg_H2,L,NY,NX)/2.0_r8

    VTATM=AZMAX1(1.2194E+04_r8*VLsoiAirP(L,NY,NX)/TKS(L,NY,NX))
!   NH3B does not have explicit gas species, so there is an inconsistency
!   with respect to the actual ebullition calculation, which involves
!   NH3B

    VTGAS=sum(trcg_VOLG(idg_beg:idg_end-1))

    IF(THETP(L,NY,NX).LT.THETX.OR.VTGAS.GT.VTATM)LX=1
    IF(THETP(L,NY,NX).GE.THETX.AND.LX.EQ.0)LG=L

    VLWatMicP1(L,NY,NX)=VLWatMicP(L,NY,NX)
    VLiceMicP1(L,NY,NX)=VLiceMicP(L,NY,NX)
    VLWatMacP1(L,NY,NX)=VLWatMacP(L,NY,NX)
    VLiceMacP1(L,NY,NX)=VLiceMacP(L,NY,NX)

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
      !exchange in the x direction
        N4=NX+1   !east
        N5=NY
        N4B=NX-1  !west
        N5B=NY
        N6=L
      ELSEIF(N.EQ.2)THEN
      !exchange in the y direction
        N4=NX
        N5=NY+1    !south
        N4B=NX
        N5B=NY-1   !north
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

          call FluxFromSnowRunoff(N,N1,N2,N4,N5,N4B,N5B)
          !
    !
        ELSEIF(N.EQ.3)THEN
          !vertical direction
          call SaltFromRunoffSnowpack(N1,N2,NY,NX)
        ENDIF
!
        ! TOTAL FLUXES FROM SEDIMENT TRANSPORT
        call TotalFluxFromSedmentTransp(N,N1,N2,N4,N5,N4B,N5B,NY,NX)
      ENDIF
!
      call FluxThruGrids(N,N1,N2,N3,N4,N5,N6,NY,NX)
    ENDDO D8580
!
!     NET FREEZE-THAW
!
!     WatFreezeThawMicP,WatFreezeThawMacP=net freeze-thaw flux in micropores,macropores
!     THeatFrezThaw=net freeze-thaw latent heat flux
!     THAW,THAWH=freeze-thaw flux in micropores,macropores from watsub.f
!     HTHAW=freeze-thaw latent heat flux from watsub.f
  !
    WatFreezeThawMicP(N3,N2,N1)=WatFreezeThawMicP(N3,N2,N1)+THAW(N3,N2,N1)
    WatFreezeThawMacP(N3,N2,N1)=WatFreezeThawMacP(N3,N2,N1)+THAWH(N3,N2,N1)
    THeatFrezThaw(N3,N2,N1)=THeatFrezThaw(N3,N2,N1)+HTHAW(N3,N2,N1)
  ENDDO D8575
  end subroutine LateralTranspt

!------------------------------------------------------------------------------------------

  subroutine ZeroRunoffArray(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,K
!
!     INITIALIZE NET WATER AND HEAT FLUXES FOR RUNOFF
!
  TQR(NY,NX)=0.0_r8
  THQR(NY,NX)=0.0_r8
!
!     INITIALIZE NET SOLUTE AND GAS FLUXES FOR RUNOFF
!
  D9960: DO K=1,micpar%n_litrsfk
    TOCQRS(K,NY,NX)=0.0_r8
    TONQRS(K,NY,NX)=0.0_r8
    TOPQRS(K,NY,NX)=0.0_r8
    TOAQRS(K,NY,NX)=0.0_r8
  ENDDO D9960
  end subroutine ZeroRunoffArray

!------------------------------------------------------------------------------------------

  subroutine ZeroFluxArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: M,NO,NGL
!     begin_execution

  call ZeroRunoffArray(NY,NX)

!     INITIALIZE NET SNOWPACK FLUXES WITHIN SNOWPACK

  call ZeroSnowArrays(NY,NX)

!
!     INITIALIZE NET SEDIMENT FLUXES FROM EROSION
!
  call ZeroErosionArray(NY,NX)

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
    TWatFlowCellMicP(L,NY,NX)=0.0_r8
    TFLWX(L,NY,NX)=0.0_r8
    TWaterFlowMacP(L,NY,NX)=0.0_r8
    THeatFlowSoiCell(L,NY,NX)=0.0_r8
    WatFreezeThawMicP(L,NY,NX)=0.0_r8
    WatFreezeThawMacP(L,NY,NX)=0.0_r8
    THeatFrezThaw(L,NY,NX)=0.0_r8
!
!     INITIALIZE GAS AND SOLUTE NET FLUX ACCUMULATORS WITHIN SOIL
!
    D8595: DO K=1,jcplx
      TOCFLS(K,L,NY,NX)=0.0_r8
      TONFLS(K,L,NY,NX)=0.0_r8
      TOPFLS(K,L,NY,NX)=0.0_r8
      TOAFLS(K,L,NY,NX)=0.0_r8
      TOCFHS(K,L,NY,NX)=0.0_r8
      TONFHS(K,L,NY,NX)=0.0_r8
      TOPFHS(K,L,NY,NX)=0.0_r8
      TOAFHS(K,L,NY,NX)=0.0_r8
    ENDDO D8595

    trcs_TFLS(ids_beg:ids_end,L,NY,NX)=0.0_r8

    trcs_TFHS(ids_beg:ids_end,L,NY,NX)=0.0_r8

    RTGasADFlx(idg_beg:idg_end-1,L,NY,NX)=0.0_r8

    IF(salt_model)THEN
      trcsa_TFLS(idsa_beg:idsab_end,L,NY,NX)=0.0_r8
      trcsa_TFHS(idsa_beg:idsab_end,L,NY,NX)=0.0_r8
    ENDIF
  ENDDO
  end subroutine ZeroFluxAccumulators

!------------------------------------------------------------------------------------------

  subroutine FluxesFromRunoff(N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B

  integer :: NN,K,NTG,NTN
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
    !water flux
    TQR(N2,N1)=TQR(N2,N1)+QR(N,NN,N2,N1)
    !heat flux
    THQR(N2,N1)=THQR(N2,N1)+HQR(N,NN,N2,N1)
    D8590: DO K=1,micpar%n_litrsfk
      TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)+XOCQRS(K,N,NN,N2,N1)
      TONQRS(K,N2,N1)=TONQRS(K,N2,N1)+XONQRS(K,N,NN,N2,N1)
      TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)+XOPQRS(K,N,NN,N2,N1)
      TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)+XOAQRS(K,N,NN,N2,N1)
    ENDDO D8590

    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
      !there is lateral runoff
      !water flux
      TQR(N2,N1)=TQR(N2,N1)-QR(N,NN,N5,N4)
      !heat flux
      THQR(N2,N1)=THQR(N2,N1)-HQR(N,NN,N5,N4)
      D8591: DO K=1,micpar%n_litrsfk
        TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)-XOCQRS(K,N,NN,N5,N4)
        TONQRS(K,N2,N1)=TONQRS(K,N2,N1)-XONQRS(K,N,NN,N5,N4)
        TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)-XOPQRS(K,N,NN,N5,N4)
        TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)-XOAQRS(K,N,NN,N5,N4)
      ENDDO D8591

    ENDIF

    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQR(N2,N1)=TQR(N2,N1)-QR(N,NN,N5B,N4B)
      THQR(N2,N1)=THQR(N2,N1)-HQR(N,NN,N5B,N4B)
      D8592: DO K=1,micpar%n_litrsfk
        TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)-XOCQRS(K,N,NN,N5B,N4B)
        TONQRS(K,N2,N1)=TONQRS(K,N2,N1)-XONQRS(K,N,NN,N5B,N4B)
        TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)-XOPQRS(K,N,NN,N5B,N4B)
        TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)-XOAQRS(K,N,NN,N5B,N4B)
      ENDDO D8592

      !     WRITE(*,6631)'TQRB',I,J,N1,N2,N4B,N5B,N,NN
      !    2,IFLBH(N,NN,N5B,N4B)
      !    2,TQR(N2,N1),QR(N,NN,N5B,N4B)
!6631  FORMAT(A8,9I4,12E12.4)
    ENDIF
  ENDDO D1202
  

  end subroutine FluxesFromRunoff
!------------------------------------------------------------------------------------------

  subroutine TotalFluxFromSedmentTransp(N,N1,N2,N4 &
    ,N5,N4B,N5B,NY,NX)
    implicit none
    integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B,NY,NX

    integer :: M,K,NO,NN,NGL,NTX,NTP
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
        TNH4ER(N2,N1)=TNH4ER(N2,N1)+XNH4ER(N,NN,N2,N1)
        TNH3ER(N2,N1)=TNH3ER(N2,N1)+XNH3ER(N,NN,N2,N1)
        TNHUER(N2,N1)=TNHUER(N2,N1)+XNHUER(N,NN,N2,N1)
        TNO3ER(N2,N1)=TNO3ER(N2,N1)+XNO3ER(N,NN,N2,N1)
        TNH4EB(N2,N1)=TNH4EB(N2,N1)+XNH4EB(N,NN,N2,N1)
        TNH3EB(N2,N1)=TNH3EB(N2,N1)+XNH3EB(N,NN,N2,N1)
        TNHUEB(N2,N1)=TNHUEB(N2,N1)+XNHUEB(N,NN,N2,N1)
        TNO3EB(N2,N1)=TNO3EB(N2,N1)+XNO3EB(N,NN,N2,N1)

        DO NTX=idx_beg,idx_end
          trcx_TER(NTX,N2,N1)=trcx_TER(NTX,N2,N1)+trcx_XER(NTX,N,NN,N2,N1)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          trcp_TER(NTP,N2,N1)=trcp_TER(NTP,N2,N1)+trcp_ER(NTP,N,NN,N2,N1)
        ENDDO

        DO  K=1,jcplx
          DO NO=1,NFGs
            DO NGL=JGnio(NO),JGnfo(NO)
              DO M=1,nlbiomcp
                TOMCER(M,NGL,K,N2,N1)=TOMCER(M,NGL,K,N2,N1) &
                  +OMCER(M+(NGL-1)*nlbiomcp,K,N,NN,N2,N1)
                TOMNER(M,NGL,K,N2,N1)=TOMNER(M,NGL,K,N2,N1) &
                  +OMNER(M+(NGL-1)*nlbiomcp,K,N,NN,N2,N1)
                TOMPER(M,NGL,K,N2,N1)=TOMPER(M,NGL,K,N2,N1) &
                  +OMPER(M+(NGL-1)*nlbiomcp,K,N,NN,N2,N1)
              enddo
            enddo
          enddo
        ENDDO

        DO NO=1,NFGs
          DO NGL=JGniA(NO),JGnfA(NO)
            DO M=1,nlbiomcp
              TOMCERff(M,NGL,N2,N1)=TOMCERff(M,NGL,N2,N1) &
                +OMCERff(M+(NGL-1)*nlbiomcp,N,NN,N2,N1)
              TOMNERff(M,NGL,N2,N1)=TOMNERff(M,NGL,N2,N1) &
                +OMNERff(M+(NGL-1)*nlbiomcp,N,NN,N2,N1)
              TOMPERff(M,NGL,N2,N1)=TOMPERff(M,NGL,N2,N1) &
                +OMPERff(M+(NGL-1)*nlbiomcp,N,NN,N2,N1)
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
        TNH4ER(N2,N1)=TNH4ER(N2,N1)-XNH4ER(N,NN,N5,N4)
        TNH3ER(N2,N1)=TNH3ER(N2,N1)-XNH3ER(N,NN,N5,N4)
        TNHUER(N2,N1)=TNHUER(N2,N1)-XNHUER(N,NN,N5,N4)
        TNO3ER(N2,N1)=TNO3ER(N2,N1)-XNO3ER(N,NN,N5,N4)
        TNH4EB(N2,N1)=TNH4EB(N2,N1)-XNH4EB(N,NN,N5,N4)
        TNH3EB(N2,N1)=TNH3EB(N2,N1)-XNH3EB(N,NN,N5,N4)
        TNHUEB(N2,N1)=TNHUEB(N2,N1)-XNHUEB(N,NN,N5,N4)
        TNO3EB(N2,N1)=TNO3EB(N2,N1)-XNO3EB(N,NN,N5,N4)

        DO NTX=idx_beg,idx_end
          trcx_TER(NTX,N2,N1)=trcx_TER(NTX,N2,N1)-trcx_XER(NTX,N,NN,N5,N4)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          trcp_TER(NTP,N2,N1)=trcp_TER(NTP,N2,N1)-trcp_ER(NTP,N,NN,N5,N4)
        ENDDO

        DO  K=1,jcplx
          DO  NO=1,NFGs
            DO NGL=JGnio(NO),JGnfo(NO)
              DO  M=1,nlbiomcp
                TOMCER(M,NGL,K,N2,N1)=TOMCER(M,NGL,K,N2,N1) &
                  -OMCER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
                TOMNER(M,NGL,K,N2,N1)=TOMNER(M,NGL,K,N2,N1) &
                  -OMNER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
                TOMPER(M,NGL,K,N2,N1)=TOMPER(M,NGL,K,N2,N1) &
                  -OMPER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
              enddo
            enddo
          enddo
        ENDDO


        DO  NO=1,NFGs
          DO  M=1,nlbiomcp
            DO NGL=JGniA(NO),JGnfA(NO)
              TOMCERff(M,NGL,N2,N1)=TOMCERff(M,NGL,N2,N1) &
                -OMCERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
              TOMNERff(M,NGL,N2,N1)=TOMNERff(M,NGL,N2,N1) &
                -OMNERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
              TOMPERff(M,NGL,N2,N1)=TOMPERff(M,NGL,N2,N1) &
                -OMPERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
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
          TNH4ER(N2,N1)=TNH4ER(N2,N1)-XNH4ER(N,NN,N5B,N4B)
          TNH3ER(N2,N1)=TNH3ER(N2,N1)-XNH3ER(N,NN,N5B,N4B)
          TNHUER(N2,N1)=TNHUER(N2,N1)-XNHUER(N,NN,N5B,N4B)
          TNO3ER(N2,N1)=TNO3ER(N2,N1)-XNO3ER(N,NN,N5B,N4B)
          TNH4EB(N2,N1)=TNH4EB(N2,N1)-XNH4EB(N,NN,N5B,N4B)
          TNH3EB(N2,N1)=TNH3EB(N2,N1)-XNH3EB(N,NN,N5B,N4B)
          TNHUEB(N2,N1)=TNHUEB(N2,N1)-XNHUEB(N,NN,N5B,N4B)
          TNO3EB(N2,N1)=TNO3EB(N2,N1)-XNO3EB(N,NN,N5B,N4B)

          DO NTX=idx_beg,idx_end
            trcx_TER(NTX,N2,N1)=trcx_TER(NTX,N2,N1)-trcx_XER(NTX,N,NN,N5B,N4B)
          ENDDO

          DO NTP=idsp_beg,idsp_end
            trcp_TER(NTP,N2,N1)=trcp_TER(NTP,N2,N1)-trcp_ER(NTP,N,NN,N5B,N4B)
          ENDDO

          D8380: DO K=1,jcplx
            DO  NO=1,NFGs
              DO NGL=JGnio(NO),JGnfo(NO)
                DO  M=1,nlbiomcp
                  TOMCER(M,NGL,K,N2,N1)=TOMCER(M,NGL,K,N2,N1) &
                    -OMCER(M+(NGL-1)*nlbiomcp,K,N,NN,N5B,N4B)
                  TOMNER(M,NGL,K,N2,N1)=TOMNER(M,NGL,K,N2,N1) &
                    -OMNER(M+(NGL-1)*nlbiomcp,K,N,NN,N5B,N4B)
                  TOMPER(M,NGL,K,N2,N1)=TOMPER(M,NGL,K,N2,N1) &
                    -OMPER(M+(NGL-1)*nlbiomcp,K,N,NN,N5B,N4B)
                enddo
              enddo
            enddo
          ENDDO D8380

          DO  NO=1,NFGs
            DO NGL=JGniA(NO),JGnfA(NO)
              DO  M=1,nlbiomcp
                TOMCERff(M,NGL,N2,N1)=TOMCERff(M,NGL,N2,N1) &
                  -OMCERff(M+(NGL-1)*nlbiomcp,N,NN,N5B,N4B)
                TOMNERff(M,NGL,N2,N1)=TOMNERff(M,NGL,N2,N1) &
                  -OMNERff(M+(NGL-1)*nlbiomcp,N,NN,N5B,N4B)
                TOMPERff(M,NGL,N2,N1)=TOMPERff(M,NGL,N2,N1) &
                  -OMPERff(M+(NGL-1)*nlbiomcp,N,NN,N5B,N4B)
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

  subroutine FluxThruGrids(N,N1,N2,N3,N4,N5,N6,NY,NX)
  implicit none
  integer, intent(in) :: N   !exchagne along direction
  integer, intent(in) :: NY,NX  !geophysical location
  integer, intent(in) :: N1,N2,N3  !source grid indices
  integer, intent(in) :: N4,N5 !dest grid indices  
  integer, intent(inout) :: N6
  integer :: LL,K,NTSA,NTS,NTG
  !     begin_execution
  !     NET HEAT, WATER FLUXES BETWEEN ADJACENT
  !     GRID CELLS
  !
  !     TFLW,TWaterFlowMacP,TWaterFlowMacP=net micropore,macropore water flux, heat flux
  !     FLW,FLWH,HFLW=micropore,macropore water flux, heat flux from watsub.f
  !     LakeSurfFlow,FLWHNU,LakeSurfHeatFlux=lake surface water flux, heat flux from watsub.f if lake surface disappears
  !when NCN /=3, it means lateral exchange is consdiered
  !N==3 means vertical direction
  IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
    D1200: DO LL=N6,NL(N5,N4)
      !modify the dest grid vertical location if needed
      IF(VLSoilPoreMicP(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
        N6=LL
        exit
      ENDIF
    ENDDO D1200

    IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(N3.EQ.NU(N2,N1).AND.N.EQ.3)THEN
        !vertical direction, source is at soil surface
        TWatFlowCellMicP(N3,N2,N1)=TWatFlowCellMicP(N3,N2,N1)+WaterFlowSoiMicP(N,N3,N2,N1)-LakeSurfFlow(N5,N4)
        TFLWX(N3,N2,N1)=TFLWX(N3,N2,N1)+FLWX(N,N3,N2,N1)-FLWXNU(N5,N4)
        TWaterFlowMacP(N3,N2,N1)=TWaterFlowMacP(N3,N2,N1)+WaterFlowMacP(N,N3,N2,N1)-FLWHNU(N5,N4)
        THeatFlowSoiCell(N3,N2,N1)=THeatFlowSoiCell(N3,N2,N1)+HeatFlow(N,N3,N2,N1)-LakeSurfHeatFlux(N5,N4)
        if(THeatFlowSoiCell(N3,N2,N1)<-1.e10)then
          write(*,*)'THeatFlowSoiCell(N3,N2,N1)+HeatFlow(N,N3,N2,N1)-LakeSurfHeatFlux(N5,N4)',&
            THeatFlowSoiCell(N3,N2,N1),HeatFlow(N,N3,N2,N1),LakeSurfHeatFlux(N5,N4)
          write(*,*)'Ns=',N1,n2,n3,n4,n5
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ELSE
        TWatFlowCellMicP(N3,N2,N1)=TWatFlowCellMicP(N3,N2,N1)+WaterFlowSoiMicP(N,N3,N2,N1)-WaterFlowSoiMicP(N,N6,N5,N4)
        TFLWX(N3,N2,N1)=TFLWX(N3,N2,N1)+FLWX(N,N3,N2,N1)-FLWX(N,N6,N5,N4)
        TWaterFlowMacP(N3,N2,N1)=TWaterFlowMacP(N3,N2,N1)+WaterFlowMacP(N,N3,N2,N1)-WaterFlowMacP(N,N6,N5,N4)
        THeatFlowSoiCell(N3,N2,N1)=THeatFlowSoiCell(N3,N2,N1)+HeatFlow(N,N3,N2,N1)-HeatFlow(N,N6,N5,N4)
        if(THeatFlowSoiCell(N3,N2,N1)<-1.e10)then
          write(*,*)'THeatFlowSoiCell(N3,N2,N1)+HeatFlow(N,N3,N2,N1)-HeatFlow(N,N6,N5,N4)',&
            THeatFlowSoiCell(N3,N2,N1),HeatFlow(N,N3,N2,N1),HeatFlow(N,N6,N5,N4)
          write(*,*)'Ns=',N,N1,n2,n3,n4,n5,n6
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ENDIF
      !     IF(N1.EQ.1.AND.N3.EQ.1)THEN
      !     WRITE(*,6632)'TFLW',I,J,N,N1,N2,N3,N4,N5,N6,NU(N2,N1)
      !    2,TWatFlowCellMicP(N3,N2,N1),WaterFlowSoiMicP(N,N3,N2,N1),WaterFlowSoiMicP(N,N6,N5,N4),LakeSurfFlow(N5,N4)
      !    3,THeatFlowSoiCell(N3,N2,N1),HeatFlow(N,N3,N2,N1),HeatFlow(N,N6,N5,N4)
      !    2,LakeSurfHeatFlux(N5,N4),VLWatMicP(N3,N2,N1)
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
      D8585: DO K=1,jcplx
        TOCFLS(K,N3,N2,N1)=TOCFLS(K,N3,N2,N1)+XOCFLS(K,N,N3,N2,N1)-XOCFLS(K,N,N6,N5,N4)
        TONFLS(K,N3,N2,N1)=TONFLS(K,N3,N2,N1)+XONFLS(K,N,N3,N2,N1)-XONFLS(K,N,N6,N5,N4)
        TOPFLS(K,N3,N2,N1)=TOPFLS(K,N3,N2,N1)+XOPFLS(K,N,N3,N2,N1)-XOPFLS(K,N,N6,N5,N4)
        TOAFLS(K,N3,N2,N1)=TOAFLS(K,N3,N2,N1)+XOAFLS(K,N,N3,N2,N1)-XOAFLS(K,N,N6,N5,N4)
        TOCFHS(K,N3,N2,N1)=TOCFHS(K,N3,N2,N1)+XOCFHS(K,N,N3,N2,N1)-XOCFHS(K,N,N6,N5,N4)
        TONFHS(K,N3,N2,N1)=TONFHS(K,N3,N2,N1)+XONFHS(K,N,N3,N2,N1)-XONFHS(K,N,N6,N5,N4)
        TOPFHS(K,N3,N2,N1)=TOPFHS(K,N3,N2,N1)+XOPFHS(K,N,N3,N2,N1)-XOPFHS(K,N,N6,N5,N4)
        TOAFHS(K,N3,N2,N1)=TOAFHS(K,N3,N2,N1)+XOAFHS(K,N,N3,N2,N1)-XOAFHS(K,N,N6,N5,N4)
      ENDDO D8585

      DO NTS=ids_beg,ids_end
        trcs_TFLS(NTS,N3,N2,N1)=trcs_TFLS(NTS,N3,N2,N1) &
          +trcs_XFLS(NTS,N,N3,N2,N1)-trcs_XFLS(NTS,N,N6,N5,N4)
        trcs_TFHS(NTS,N3,N2,N1)=trcs_TFHS(NTS,N3,N2,N1) &
          +trcs_XFHS(NTS,N,N3,N2,N1)-trcs_XFHS(NTS,N,N6,N5,N4)
      ENDDO
!
      !     NET GAS FLUXES BETWEEN ADJACENT GRID CELLS
      !
      !     T*FLG=net convective+diffusive gas flux
      !     X*FLG=convective+diffusive gas flux from trnsfr.f
      !     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!exclude NH3B
      DO NTG=idg_beg,idg_end-1
        RTGasADFlx(NTG,N3,N2,N1)=RTGasADFlx(NTG,N3,N2,N1)+R3GasADTFlx(NTG,N,N3,N2,N1)-R3GasADTFlx(NTG,N,N6,N5,N4)
      ENDDO
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
      IF(salt_model)THEN
        DO NTSA=idsa_beg,idsab_end
          trcsa_TFLS(NTSA,N3,N2,N1)=trcsa_TFLS(NTSA,N3,N2,N1) &
            +trcsa_XFLS(NTSA,N,N3,N2,N1)-trcsa_XFLS(NTSA,N,N6,N5,N4)
          trcsa_TFHS(NTSA,N3,N2,N1)=trcsa_TFHS(NTSA,N3,N2,N1) &
            +trcsa_XFHS(NTSA,N,N3,N2,N1)-trcsa_XFHS(NTSA,N,N6,N5,N4)
        ENDDO
      ENDIF
    ELSE
      TWatFlowCellMicP(N3,N2,N1)=0.0_r8
      TFLWX(N3,N2,N1)=0.0_r8
      TWaterFlowMacP(N3,N2,N1)=0.0_r8
      THeatFlowSoiCell(N3,N2,N1)=0.0_r8
      WatFreezeThawMicP(N3,N2,N1)=0.0_r8
      WatFreezeThawMacP(N3,N2,N1)=0.0_r8
      THeatFrezThaw(N3,N2,N1)=0.0_r8
      D8596: DO K=1,jcplx
        TOCFLS(K,N3,N2,N1)=0.0_r8
        TONFLS(K,N3,N2,N1)=0.0_r8
        TOPFLS(K,N3,N2,N1)=0.0_r8
        TOAFLS(K,N3,N2,N1)=0.0_r8
        TOCFHS(K,N3,N2,N1)=0.0_r8
        TONFHS(K,N3,N2,N1)=0.0_r8
        TOPFHS(K,N3,N2,N1)=0.0_r8
        TOAFHS(K,N3,N2,N1)=0.0_r8
      ENDDO D8596
      trcs_TFLS(ids_beg:ids_end,N3,N2,N1)=0.0_r8

      trcs_TFHS(ids_beg:ids_end,N3,N2,N1)=0.0_r8

      RTGasADFlx(idg_beg:idg_end-1,N3,N2,N1)=0.0_r8

      IF(salt_model)THEN
        trcsa_TFLS(idsa_beg:idsab_end,N3,N2,N1)=0.0_r8
        trcsa_TFHS(idsa_beg:idsab_end,N3,N2,N1)=0.0_r8
      ENDIF
    ENDIF
  ENDIF
  end subroutine FluxThruGrids

end module LateralTranspMod
