module BoundaryTranspMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb
  use GridConsts
  USE SoilPropertyDataType
  use TransfrDataMod
  use GridDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use IrrigationDataType
  use EcoSIMSolverPar
  use SnowDataType
  use SoilBGCDataType
  use ChemTranspDataType
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  public :: BoundaryFlux
  contains

!------------------------------------------------------------------------------------------

  subroutine BoundaryFlux(M,MX,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) :: M,MX,NHW, NHE, NVN, NVS

  integer :: NY,NX,L
  integer :: N1,N2,N3,N4,N5,N6,NN,N
  integer :: N4B,N5B,M1,M2,M3,M4,M5,M6
  real(r8) :: RCHQF,RCHGFU,RCHGFT,XN

!     N3,N2,N1=L,NY,NX of source grid cell
!     M6,M5,M4=L,NY,NX of destination grid cell
!
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU(NY,NX),NL(NY,NX)
        N1=NX
        N2=NY
        N3=L
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
        DO  N=1,3
          DO  NN=1,2
            IF(N.EQ.1)THEN
              N4=NX+1
              N5=NY
              N4B=NX-1
              N5B=NY
              N6=L
              IF(NN.EQ.1)THEN
                IF(NX.EQ.NHE)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX+1
                  M5=NY
                  M6=L
                  XN=-1.0
                  RCHQF=RCHQE(M2,M1)
                  RCHGFU=RCHGEU(M2,M1)
                  RCHGFT=RCHGET(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                IF(NX.EQ.NHW)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0
                  RCHQF=RCHQW(M5,M4)
                  RCHGFU=RCHGWU(M5,M4)
                  RCHGFT=RCHGWT(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.2)THEN
              N4=NX
              N5=NY+1
              N4B=NX
              N5B=NY-1
              N6=L
              IF(NN.EQ.1)THEN
                IF(NY.EQ.NVS)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY+1
                  M6=L
                  XN=-1.0
                  RCHQF=RCHQS(M2,M1)
                  RCHGFU=RCHGSU(M2,M1)
                  RCHGFT=RCHGST(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0
                  RCHQF=RCHQN(M5,M4)
                  RCHGFU=RCHGNU(M5,M4)
                  RCHGFT=RCHGNT(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.3)THEN
              N1=NX
              N2=NY
              N3=L
              N4=NX
              N5=NY
              N6=L+1
              IF(NN.EQ.1)THEN
                IF(L.EQ.NL(NY,NX))THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L+1
                  XN=-1.0
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                cycle
              ENDIF
            ENDIF

            IF(M.NE.MX)THEN
!
!     SURFACE SOLUTE TRANSPORT FROM BOUNDARY SURFACE
!     RUNOFF IN 'WATSUB' AND CONCENTRATIONS IN THE SURFACE SOIL LAYER
!
              call BoundaryRunoffandSnow(L,N,NN,M,MX,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
!
!     SOLUTE LOSS WITH SUBSURFACE MICROPORE WATER LOSS
!
            ENDIF
          ENDDO
!
!     NET GAS AND SOLUTE FLUXES IN EACH GRID CELL
!C
!     NET OVERLAND SOLUTE FLUX IN WATER
!
          call NetOverlandFlux(L,N,M,MX,NY,NX,N1,N2,N4B,N5B,N4,N5)
!
!     TOTAL SOLUTE FLUX IN MICROPORES AND MACROPORES
!
!     T*FLS=net convective + diffusive solute flux through micropores
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,R*FLB=convective + diffusive solute flux through micropores in non-band,band
!     T*FHS=net convective + diffusive solute flux through macropores
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,R*FHB=convective + diffusive solute flux through macropores in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
          IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
            call NetFluxMicroandMacropores(NY,NX,N,M,MX,N1,N2,N3,N4,N5,N6)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine BoundaryFlux
!------------------------------------------------------------------------------------------

  subroutine BoundaryRunoffandSnow(L,N,NN,M,MX,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
  implicit none

  integer, intent(in) :: L,N, NN, M,MX, N1, N2,N3, M1, M2,M3,M4, M5,M6
  real(r8), intent(in):: RCHQF
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K

! begin_execution
!
!     QRM =runoff from watsub.f
!     RQR*=solute in runoff
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  IF(M.NE.MX)THEN
    IF(L.EQ.NUM(M2,M1).AND.N.NE.3)THEN
      IF(IRCHG(NN,N,N2,N1).EQ.0.OR.test_aeqb(RCHQF,0.0_r8) &
        .OR.QRM(M,N2,N1).LE.ZEROS(N2,N1))THEN
        DO  K=0,jcplx1
          RQROC(K,N,NN,M5,M4)=0.0
          RQRON(K,N,NN,M5,M4)=0.0
          RQROP(K,N,NN,M5,M4)=0.0
          RQROA(K,N,NN,M5,M4)=0.0
        ENDDO
        RQRCOS(N,NN,M5,M4)=0.0
        RQRCHS(N,NN,M5,M4)=0.0
        RQROXS(N,NN,M5,M4)=0.0
        RQRNGS(N,NN,M5,M4)=0.0
        RQRN2S(N,NN,M5,M4)=0.0
        RQRHGS(N,NN,M5,M4)=0.0
        RQRNH4(N,NN,M5,M4)=0.0
        RQRNH3(N,NN,M5,M4)=0.0
        RQRNO3(N,NN,M5,M4)=0.0
        RQRNO2(N,NN,M5,M4)=0.0
        RQRH1P(N,NN,M5,M4)=0.0
        RQRH2P(N,NN,M5,M4)=0.0
      ELSE
!
!     SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
      IF((NN.EQ.1.AND.QRMN(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.2.AND.QRMN(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      FQRM=QRMN(M,N,NN,M5,M4)/QRM(M,N2,N1)
      DO 9540 K=0,jcplx1
      RQROC(K,N,NN,M5,M4)=RQROC0(K,N2,N1)*FQRM
      RQRON(K,N,NN,M5,M4)=RQRON0(K,N2,N1)*FQRM
      RQROP(K,N,NN,M5,M4)=RQROP0(K,N2,N1)*FQRM
      RQROA(K,N,NN,M5,M4)=RQROA0(K,N2,N1)*FQRM
9540  CONTINUE
      RQRCOS(N,NN,M5,M4)=RQRCOS0(N2,N1)*FQRM
      RQRCHS(N,NN,M5,M4)=RQRCHS0(N2,N1)*FQRM
      RQROXS(N,NN,M5,M4)=RQROXS0(N2,N1)*FQRM
      RQRNGS(N,NN,M5,M4)=RQRNGS0(N2,N1)*FQRM
      RQRN2S(N,NN,M5,M4)=RQRN2S0(N2,N1)*FQRM
      RQRHGS(N,NN,M5,M4)=RQRHGS0(N2,N1)*FQRM
      RQRNH4(N,NN,M5,M4)=RQRNH40(N2,N1)*FQRM
      RQRNH3(N,NN,M5,M4)=RQRNH30(N2,N1)*FQRM
      RQRNO3(N,NN,M5,M4)=RQRNO30(N2,N1)*FQRM
      RQRNO2(N,NN,M5,M4)=RQRNO20(N2,N1)*FQRM
      RQRH1P(N,NN,M5,M4)=RQRH1P0(N2,N1)*FQRM
      RQRH2P(N,NN,M5,M4)=RQRH2P0(N2,N1)*FQRM
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QRS=hourly solute in runoff
!     RQR*=solute in runoff
!
      DO 9565 K=0,jcplx1
      XOCQRS(K,N,NN,M5,M4)=XOCQRS(K,N,NN,M5,M4)+RQROC(K,N,NN,M5,M4)
      XONQRS(K,N,NN,M5,M4)=XONQRS(K,N,NN,M5,M4)+RQRON(K,N,NN,M5,M4)
      XOPQRS(K,N,NN,M5,M4)=XOPQRS(K,N,NN,M5,M4)+RQROP(K,N,NN,M5,M4)
      XOAQRS(K,N,NN,M5,M4)=XOAQRS(K,N,NN,M5,M4)+RQROA(K,N,NN,M5,M4)
9565  CONTINUE
      XCOQRS(N,NN,M5,M4)=XCOQRS(N,NN,M5,M4)+RQRCOS(N,NN,M5,M4)
      XCHQRS(N,NN,M5,M4)=XCHQRS(N,NN,M5,M4)+RQRCHS(N,NN,M5,M4)
      XOXQRS(N,NN,M5,M4)=XOXQRS(N,NN,M5,M4)+RQROXS(N,NN,M5,M4)
      XNGQRS(N,NN,M5,M4)=XNGQRS(N,NN,M5,M4)+RQRNGS(N,NN,M5,M4)
      XN2QRS(N,NN,M5,M4)=XN2QRS(N,NN,M5,M4)+RQRN2S(N,NN,M5,M4)
      XHGQRS(N,NN,M5,M4)=XHGQRS(N,NN,M5,M4)+RQRHGS(N,NN,M5,M4)
      XN4QRW(N,NN,M5,M4)=XN4QRW(N,NN,M5,M4)+RQRNH4(N,NN,M5,M4)
      XN3QRW(N,NN,M5,M4)=XN3QRW(N,NN,M5,M4)+RQRNH3(N,NN,M5,M4)
      XNOQRW(N,NN,M5,M4)=XNOQRW(N,NN,M5,M4)+RQRNO3(N,NN,M5,M4)
      XNXQRS(N,NN,M5,M4)=XNXQRS(N,NN,M5,M4)+RQRNO2(N,NN,M5,M4)
      XP1QRW(N,NN,M5,M4)=XP1QRW(N,NN,M5,M4)+RQRH1P(N,NN,M5,M4)
      XP4QRW(N,NN,M5,M4)=XP4QRW(N,NN,M5,M4)+RQRH2P(N,NN,M5,M4)
!
!     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
      ELSEIF((NN.EQ.2.AND.QRMN(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.1.AND.QRMN(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      DO 9629 K=0,jcplx1
      RQROC(K,N,NN,M5,M4)=0.0
      RQRON(K,N,NN,M5,M4)=0.0
      RQROP(K,N,NN,M5,M4)=0.0
      RQROA(K,N,NN,M5,M4)=0.0
9629  CONTINUE
      RQRCOS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CCOU
      RQRCHS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CCHU
      RQROXS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*COXU
      RQRNGS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CNNU
      RQRN2S(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CN2U
      RQRHGS(N,NN,M5,M4)=0.0
      RQRNH4(N,NN,M5,M4)=0.0
      RQRNH3(N,NN,M5,M4)=0.0
      RQRNO3(N,NN,M5,M4)=0.0
      RQRNO2(N,NN,M5,M4)=0.0
      RQRH1P(N,NN,M5,M4)=0.0
      RQRH2P(N,NN,M5,M4)=0.0
      XCOQRS(N,NN,M5,M4)=XCOQRS(N,NN,M5,M4)+RQRCOS(N,NN,M5,M4)
      XCHQRS(N,NN,M5,M4)=XCHQRS(N,NN,M5,M4)+RQRCHS(N,NN,M5,M4)
      XOXQRS(N,NN,M5,M4)=XOXQRS(N,NN,M5,M4)+RQROXS(N,NN,M5,M4)
      XNGQRS(N,NN,M5,M4)=XNGQRS(N,NN,M5,M4)+RQRNGS(N,NN,M5,M4)
      XN2QRS(N,NN,M5,M4)=XN2QRS(N,NN,M5,M4)+RQRN2S(N,NN,M5,M4)
      ELSE
      DO 9627 K=0,jcplx1
      RQROC(K,N,NN,M5,M4)=0.0
      RQRON(K,N,NN,M5,M4)=0.0
      RQROP(K,N,NN,M5,M4)=0.0
      RQROA(K,N,NN,M5,M4)=0.0
9627  CONTINUE
      RQRCOS(N,NN,M5,M4)=0.0
      RQRCHS(N,NN,M5,M4)=0.0
      RQROXS(N,NN,M5,M4)=0.0
      RQRNGS(N,NN,M5,M4)=0.0
      RQRN2S(N,NN,M5,M4)=0.0
      RQRHGS(N,NN,M5,M4)=0.0
      RQRNH4(N,NN,M5,M4)=0.0
      RQRNH3(N,NN,M5,M4)=0.0
      RQRNO3(N,NN,M5,M4)=0.0
      RQRNO2(N,NN,M5,M4)=0.0
      RQRH1P(N,NN,M5,M4)=0.0
      RQRH2P(N,NN,M5,M4)=0.0
      ENDIF
!     WRITE(*,1114)'RQR',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
!    2,QRM(M,N2,N1),QRMN(M,N,NN,M5,M4)
!    3,RQROC0(1,N2,N1),RQROC(1,N,NN,M5,M4),XOCQRS(1,N,NN,M5,M4)
!    3,RQROP0(1,N2,N1),RQROP(1,N,NN,M5,M4),XOPQRS(1,N,NN,M5,M4)
!1114  FORMAT(A8,11I4,30E12.4)
      ENDIF
!
!     BOUNDARY SNOW FLUX
!
      IF(NN.EQ.1)THEN
      RQSCOS(N,M5,M4)=0.0
      RQSCHS(N,M5,M4)=0.0
      RQSOXS(N,M5,M4)=0.0
      RQSNGS(N,M5,M4)=0.0
      RQSN2S(N,M5,M4)=0.0
      RQSNH4(N,M5,M4)=0.0
      RQSNH3(N,M5,M4)=0.0
      RQSNO3(N,M5,M4)=0.0
      RQSH1P(N,M5,M4)=0.0
      RQSH2P(N,M5,M4)=0.0
      ENDIF
      ENDIF
!     FLWM=water flux through soil micropore from watsub.f
!     VOLWM=micropore water-filled porosity from watsub.f
!     R*FLS=convective solute flux through micropores
!     R*FLW,R*FLB=convective solute flux through micropores in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(NCN(M2,M1).NE.3.OR.N.EQ.3)THEN
      IF(NN.EQ.1.AND.FLWM(M,N,M6,M5,M4).GT.0.0 &
      .OR.NN.EQ.2.AND.FLWM(M,N,M6,M5,M4).LT.0.0)THEN
      IF(VOLWM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,FLWM(M,N,M6,M5,M4) &
      /VOLWM(M,M3,M2,M1)))
      ELSE
      VFLW=0.0
      ENDIF
      DO 9520 K=0,jcplx1
      ROCFLS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQC2(K,M3,M2,M1))
      RONFLS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQN2(K,M3,M2,M1))
      ROPFLS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQP2(K,M3,M2,M1))
      ROAFLS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQA2(K,M3,M2,M1))
9520  CONTINUE
      RCOFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,CO2S2(M3,M2,M1))
      RCHFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,CH4S2(M3,M2,M1))
      ROXFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,OXYS2(M3,M2,M1))
      RNGFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,Z2GS2(M3,M2,M1))
      RN2FLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,Z2OS2(M3,M2,M1))
      RHGFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H2GS2(M3,M2,M1))
      RN4FLW(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNH4S2(M3,M2,M1)) &
      *VLNH4(M3,M2,M1)
      RN3FLW(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNH3S2(M3,M2,M1)) &
      *VLNH4(M3,M2,M1)
      RNOFLW(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNO3S2(M3,M2,M1)) &
      *VLNO3(M3,M2,M1)
      RNXFLS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNO2S2(M3,M2,M1)) &
      *VLNO3(M3,M2,M1)
      RH1PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H1PO42(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RH2PFS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H2PO42(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RN4FLB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNH4B2(M3,M2,M1)) &
      *VLNHB(M3,M2,M1)
      RN3FLB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNH3B2(M3,M2,M1)) &
      *VLNHB(M3,M2,M1)
      RNOFLB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNO3B2(M3,M2,M1)) &
      *VLNOB(M3,M2,M1)
      RNXFLB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNO2B2(M3,M2,M1)) &
      *VLNOB(M3,M2,M1)
      RH2BFB(N,M6,M5,M4)=VFLW*AMAX1(0.0,H2POB2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
!
!     SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
!
      ELSE
      DO 9515 K=0,jcplx1
      ROCFLS(K,N,M6,M5,M4)=0.0
      RONFLS(K,N,M6,M5,M4)=0.0
      ROPFLS(K,N,M6,M5,M4)=0.0
      ROAFLS(K,N,M6,M5,M4)=0.0
9515  CONTINUE
      RCOFLS(N,M6,M5,M4)=0.0
      RCHFLS(N,M6,M5,M4)=0.0
      ROXFLS(N,M6,M5,M4)=0.0
      RNGFLS(N,M6,M5,M4)=0.0
      RN2FLS(N,M6,M5,M4)=0.0
      RHGFLS(N,M6,M5,M4)=0.0
      RN4FLW(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN4U(M3,M2,M1) &
      *VLNH4(M3,M2,M1)
      RN3FLW(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN3U(M3,M2,M1) &
      *VLNH4(M3,M2,M1)
      RNOFLW(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNOU(M3,M2,M1) &
      *VLNO3(M3,M2,M1)
      RNXFLS(N,M6,M5,M4)=0.0
      RH1PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH1PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RH2PFS(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH2PU(M3,M2,M1) &
      *VLPO4(M3,M2,M1)
      RN4FLB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN4U(M3,M2,M1) &
      *VLNHB(M3,M2,M1)
      RN3FLB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN3U(M3,M2,M1) &
      *VLNHB(M3,M2,M1)
      RNOFLB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNOU(M3,M2,M1) &
      *VLNOB(M3,M2,M1)
      RNXFLB(N,M6,M5,M4)=0.0
      RH1BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH1PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      RH2BFB(N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH2PU(M3,M2,M1) &
      *VLPOB(M3,M2,M1)
      ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
!
!     FLWHM=water flux through soil macropore from watsub.f
!     VOLWHM=macropore water-filled porosity from watsub.f
!     RFH*S=solute diffusive flux through macropore
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
      IF(NN.EQ.1.AND.FLWHM(M,N,M6,M5,M4).GT.0.0 &
      .OR.NN.EQ.2.AND.FLWHM(M,N,M6,M5,M4).LT.0.0)THEN
      IF(VOLWHM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,FLWHM(M,N,M6,M5,M4) &
      /VOLWHM(M,M3,M2,M1)))
      ELSE
      VFLW=0.0
      ENDIF
      DO 9535 K=0,jcplx1
      ROCFHS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQCH2(K,M3,M2,M1))
      RONFHS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQNH2(K,M3,M2,M1))
      ROPFHS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQPH2(K,M3,M2,M1))
      ROAFHS(K,N,M6,M5,M4)=VFLW*AMAX1(0.0,OQAH2(K,M3,M2,M1))
9535  CONTINUE
      RCOFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,CO2SH2(M3,M2,M1))
      RCHFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,CH4SH2(M3,M2,M1))
      ROXFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,OXYSH2(M3,M2,M1))
      RNGFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,Z2GSH2(M3,M2,M1))
      RN2FHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,Z2OSH2(M3,M2,M1))
      RHGFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H2GSH2(M3,M2,M1))
      RN4FHW(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNH4H2(M3,M2,M1)) &
      *VLNH4(M3,M2,M1)
      RN3FHW(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNH3H2(M3,M2,M1)) &
      *VLNH4(M3,M2,M1)
      RNOFHW(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNO3H2(M3,M2,M1)) &
      *VLNO3(M3,M2,M1)
      RNXFHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNO2H2(M3,M2,M1)) &
      *VLNO3(M3,M2,M1)
      RH1PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H1P4H2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RH2PHS(N,M6,M5,M4)=VFLW*AMAX1(0.0,H2P4H2(M3,M2,M1)) &
      *VLPO4(M3,M2,M1)
      RN4FHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZN4BH2(M3,M2,M1)) &
      *VLNHB(M3,M2,M1)
      RN3FHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZN3BH2(M3,M2,M1)) &
      *VLNHB(M3,M2,M1)
      RNOFHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZNOBH2(M3,M2,M1)) &
      *VLNOB(M3,M2,M1)
      RNXFHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZN2BH2(M3,M2,M1)) &
      *VLNOB(M3,M2,M1)
      RH1BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,H1PBH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
      RH2BHB(N,M6,M5,M4)=VFLW*AMAX1(0.0,H2PBH2(M3,M2,M1)) &
      *VLPOB(M3,M2,M1)
!
!     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
!
      ELSE
      DO 9530 K=0,jcplx1
      ROCFHS(K,N,M6,M5,M4)=0.0
      RONFHS(K,N,M6,M5,M4)=0.0
      ROPFHS(K,N,M6,M5,M4)=0.0
      ROAFHS(K,N,M6,M5,M4)=0.0
9530  CONTINUE
      RCOFHS(N,M6,M5,M4)=0.0
      RCHFHS(N,M6,M5,M4)=0.0
      ROXFHS(N,M6,M5,M4)=0.0
      RNGFHS(N,M6,M5,M4)=0.0
      RN2FHS(N,M6,M5,M4)=0.0
      RN4FHW(N,M6,M5,M4)=0.0
      RHGFHS(N,M6,M5,M4)=0.0
      RN3FHW(N,M6,M5,M4)=0.0
      RNOFHW(N,M6,M5,M4)=0.0
      RNXFHS(N,M6,M5,M4)=0.0
      RH1PHS(N,M6,M5,M4)=0.0
      RH2PHS(N,M6,M5,M4)=0.0
      RN4FHB(N,M6,M5,M4)=0.0
      RN3FHB(N,M6,M5,M4)=0.0
      RNOFHB(N,M6,M5,M4)=0.0
      RNXFHB(N,M6,M5,M4)=0.0
      RH1BHB(N,M6,M5,M4)=0.0
      RH2BHB(N,M6,M5,M4)=0.0
      ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS,X*FLW,X*FLB=hourly solute flux in non-band,band micropores
!     X*FHS,X*FHW,X*FHB=hourly solute flux in non-band,band macropores
!     R*FLS,R*FLW,R*FLB=solute flux in non-band,band micropores
!     R*FHS,R*FHW,R*FHB=solute flux in non-band,band macropores
!
      DO 9555 K=0,jcplx1
      XOCFLS(K,N,M6,M5,M4)=XOCFLS(K,N,M6,M5,M4)+ROCFLS(K,N,M6,M5,M4)
      XONFLS(K,N,M6,M5,M4)=XONFLS(K,N,M6,M5,M4)+RONFLS(K,N,M6,M5,M4)
      XOPFLS(K,N,M6,M5,M4)=XOPFLS(K,N,M6,M5,M4)+ROPFLS(K,N,M6,M5,M4)
      XOAFLS(K,N,M6,M5,M4)=XOAFLS(K,N,M6,M5,M4)+ROAFLS(K,N,M6,M5,M4)
      XOCFHS(K,N,M6,M5,M4)=XOCFHS(K,N,M6,M5,M4)+ROCFHS(K,N,M6,M5,M4)
      XONFHS(K,N,M6,M5,M4)=XONFHS(K,N,M6,M5,M4)+RONFHS(K,N,M6,M5,M4)
      XOPFHS(K,N,M6,M5,M4)=XOPFHS(K,N,M6,M5,M4)+ROPFHS(K,N,M6,M5,M4)
      XOAFHS(K,N,M6,M5,M4)=XOAFHS(K,N,M6,M5,M4)+ROAFHS(K,N,M6,M5,M4)
9555  CONTINUE
      XCOFLS(N,M6,M5,M4)=XCOFLS(N,M6,M5,M4)+RCOFLS(N,M6,M5,M4)
      XCHFLS(N,M6,M5,M4)=XCHFLS(N,M6,M5,M4)+RCHFLS(N,M6,M5,M4)
      XOXFLS(N,M6,M5,M4)=XOXFLS(N,M6,M5,M4)+ROXFLS(N,M6,M5,M4)
      XNGFLS(N,M6,M5,M4)=XNGFLS(N,M6,M5,M4)+RNGFLS(N,M6,M5,M4)
      XN2FLS(N,M6,M5,M4)=XN2FLS(N,M6,M5,M4)+RN2FLS(N,M6,M5,M4)
      XHGFLS(N,M6,M5,M4)=XHGFLS(N,M6,M5,M4)+RHGFLS(N,M6,M5,M4)
      XN4FLW(N,M6,M5,M4)=XN4FLW(N,M6,M5,M4)+RN4FLW(N,M6,M5,M4)
      XN3FLW(N,M6,M5,M4)=XN3FLW(N,M6,M5,M4)+RN3FLW(N,M6,M5,M4)
      XNOFLW(N,M6,M5,M4)=XNOFLW(N,M6,M5,M4)+RNOFLW(N,M6,M5,M4)
      XNXFLS(N,M6,M5,M4)=XNXFLS(N,M6,M5,M4)+RNXFLS(N,M6,M5,M4)
      XH1PFS(N,M6,M5,M4)=XH1PFS(N,M6,M5,M4)+RH1PFS(N,M6,M5,M4)
      XH2PFS(N,M6,M5,M4)=XH2PFS(N,M6,M5,M4)+RH2PFS(N,M6,M5,M4)
      XN4FLB(N,M6,M5,M4)=XN4FLB(N,M6,M5,M4)+RN4FLB(N,M6,M5,M4)
      XN3FLB(N,M6,M5,M4)=XN3FLB(N,M6,M5,M4)+RN3FLB(N,M6,M5,M4)
      XNOFLB(N,M6,M5,M4)=XNOFLB(N,M6,M5,M4)+RNOFLB(N,M6,M5,M4)
      XNXFLB(N,M6,M5,M4)=XNXFLB(N,M6,M5,M4)+RNXFLB(N,M6,M5,M4)
      XH2BFB(N,M6,M5,M4)=XH2BFB(N,M6,M5,M4)+RH2BFB(N,M6,M5,M4)
      XCOFHS(N,M6,M5,M4)=XCOFHS(N,M6,M5,M4)+RCOFHS(N,M6,M5,M4)
      XCHFHS(N,M6,M5,M4)=XCHFHS(N,M6,M5,M4)+RCHFHS(N,M6,M5,M4)
      XOXFHS(N,M6,M5,M4)=XOXFHS(N,M6,M5,M4)+ROXFHS(N,M6,M5,M4)
      XNGFHS(N,M6,M5,M4)=XNGFHS(N,M6,M5,M4)+RNGFHS(N,M6,M5,M4)
      XN2FHS(N,M6,M5,M4)=XN2FHS(N,M6,M5,M4)+RN2FHS(N,M6,M5,M4)
      XHGFHS(N,M6,M5,M4)=XHGFHS(N,M6,M5,M4)+RHGFHS(N,M6,M5,M4)
      XN4FHW(N,M6,M5,M4)=XN4FHW(N,M6,M5,M4)+RN4FHW(N,M6,M5,M4)
      XN3FHW(N,M6,M5,M4)=XN3FHW(N,M6,M5,M4)+RN3FHW(N,M6,M5,M4)
      XNOFHW(N,M6,M5,M4)=XNOFHW(N,M6,M5,M4)+RNOFHW(N,M6,M5,M4)
      XNXFHS(N,M6,M5,M4)=XNXFHS(N,M6,M5,M4)+RNXFHS(N,M6,M5,M4)
      XH1PHS(N,M6,M5,M4)=XH1PHS(N,M6,M5,M4)+RH1PHS(N,M6,M5,M4)
      XH2PHS(N,M6,M5,M4)=XH2PHS(N,M6,M5,M4)+RH2PHS(N,M6,M5,M4)
      XN4FHB(N,M6,M5,M4)=XN4FHB(N,M6,M5,M4)+RN4FHB(N,M6,M5,M4)
      XN3FHB(N,M6,M5,M4)=XN3FHB(N,M6,M5,M4)+RN3FHB(N,M6,M5,M4)
      XNOFHB(N,M6,M5,M4)=XNOFHB(N,M6,M5,M4)+RNOFHB(N,M6,M5,M4)
      XNXFHB(N,M6,M5,M4)=XNXFHB(N,M6,M5,M4)+RNXFHB(N,M6,M5,M4)
      XH1BHB(N,M6,M5,M4)=XH1BHB(N,M6,M5,M4)+RH1BHB(N,M6,M5,M4)
      XH2BHB(N,M6,M5,M4)=XH2BHB(N,M6,M5,M4)+RH2BHB(N,M6,M5,M4)
      ENDIF
      ENDIF
!
!     GASOUS LOSS WITH SUBSURFACE MICROPORE WATER GAIN
!
!     FLWM,FLWHM=micropore,macropore water flux from watsub.f
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     VOLPM=air-filled porosity
!     R*FLG=convective gas flux
!     X*FLG=convective gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
      FLGM=(FLWM(M,N,M6,M5,M4)+FLWHM(M,N,M6,M5,M4))*XNPT
      IF(NN.EQ.1.AND.FLGM.LT.0.0 &
      .OR.NN.EQ.2.AND.FLGM.GT.0.0)THEN
      IF(VOLPM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=-AMAX1(-VFLWX,AMIN1(VFLWX,FLGM &
      /VOLPM(M,M3,M2,M1)))
      ELSE
      VFLW=0.0
      ENDIF
      RCOFLG(N,M6,M5,M4)=VFLW*AMAX1(0.0,CO2G2(M3,M2,M1))
      RCHFLG(N,M6,M5,M4)=VFLW*AMAX1(0.0,CH4G2(M3,M2,M1))
      ROXFLG(N,M6,M5,M4)=VFLW*AMAX1(0.0,OXYG2(M3,M2,M1))
      RNGFLG(N,M6,M5,M4)=VFLW*AMAX1(0.0,Z2GG2(M3,M2,M1))
      RN2FLG(N,M6,M5,M4)=VFLW*AMAX1(0.0,Z2OG2(M3,M2,M1))
      RN3FLG(N,M6,M5,M4)=VFLW*AMAX1(0.0,ZN3G2(M3,M2,M1))
      RHGFLG(N,M6,M5,M4)=VFLW*AMAX1(0.0,H2GG2(M3,M2,M1))
!
!     HOURLY GAS FLUX FOR USE IN REDIST.F
!
!     X*FLG=hourly convective gas flux
!
      XCOFLG(N,M6,M5,M4)=XCOFLG(N,M6,M5,M4)+RCOFLG(N,M6,M5,M4)
      XCHFLG(N,M6,M5,M4)=XCHFLG(N,M6,M5,M4)+RCHFLG(N,M6,M5,M4)
      XOXFLG(N,M6,M5,M4)=XOXFLG(N,M6,M5,M4)+ROXFLG(N,M6,M5,M4)
      XNGFLG(N,M6,M5,M4)=XNGFLG(N,M6,M5,M4)+RNGFLG(N,M6,M5,M4)
      XN2FLG(N,M6,M5,M4)=XN2FLG(N,M6,M5,M4)+RN2FLG(N,M6,M5,M4)
      XN3FLG(N,M6,M5,M4)=XN3FLG(N,M6,M5,M4)+RN3FLG(N,M6,M5,M4)
      XHGFLG(N,M6,M5,M4)=XHGFLG(N,M6,M5,M4)+RHGFLG(N,M6,M5,M4)
      ELSE
      RCOFLG(N,M6,M5,M4)=0.0
      RCHFLG(N,M6,M5,M4)=0.0
      ROXFLG(N,M6,M5,M4)=0.0
      RNGFLG(N,M6,M5,M4)=0.0
      RN2FLG(N,M6,M5,M4)=0.0
      RN3FLG(N,M6,M5,M4)=0.0
      RHGFLG(N,M6,M5,M4)=0.0
      ENDIF
      ELSE
      DO 9531 K=0,jcplx1
      ROCFHS(K,N,M6,M5,M4)=0.0
      RONFHS(K,N,M6,M5,M4)=0.0
      ROPFHS(K,N,M6,M5,M4)=0.0
      ROAFHS(K,N,M6,M5,M4)=0.0
9531  CONTINUE
      RCOFHS(N,M6,M5,M4)=0.0
      RCHFHS(N,M6,M5,M4)=0.0
      ROXFHS(N,M6,M5,M4)=0.0
      RNGFHS(N,M6,M5,M4)=0.0
      RN2FHS(N,M6,M5,M4)=0.0
      RN4FHW(N,M6,M5,M4)=0.0
      RHGFHS(N,M6,M5,M4)=0.0
      RN3FHW(N,M6,M5,M4)=0.0
      RNOFHW(N,M6,M5,M4)=0.0
      RNXFHS(N,M6,M5,M4)=0.0
      RH1PHS(N,M6,M5,M4)=0.0
      RH2PHS(N,M6,M5,M4)=0.0
      RN4FHB(N,M6,M5,M4)=0.0
      RN3FHB(N,M6,M5,M4)=0.0
      RNOFHB(N,M6,M5,M4)=0.0
      RNXFHB(N,M6,M5,M4)=0.0
      RH1BHB(N,M6,M5,M4)=0.0
      RH2BHB(N,M6,M5,M4)=0.0
      RCOFLG(N,M6,M5,M4)=0.0
      RCHFLG(N,M6,M5,M4)=0.0
      ROXFLG(N,M6,M5,M4)=0.0
      RNGFLG(N,M6,M5,M4)=0.0
      RN2FLG(N,M6,M5,M4)=0.0
      RN3FLG(N,M6,M5,M4)=0.0
      RHGFLG(N,M6,M5,M4)=0.0
      ENDIF
      end subroutine BoundaryRunoffandSnow
!------------------------------------------------------------------------------------------
!
      subroutine NetOverlandFlux(L,N,M,MX,NY,NX,N1,N2,N4B,N5B,N4,N5)
      implicit none

      integer, intent(in) :: L,N, M, MX,NY,NX,N1,N2,N4B,N5B,N4,N5
  integer :: NN,K,LS,LS2
!
!     TQR*=net overland solute flux
!     RQR*=overland solute flux
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
      IF(M.NE.MX)THEN
      IF(L.EQ.NUM(N2,N1))THEN
      IF(N.NE.3)THEN

      DO 1202 NN=1,2
      DO 9550 K=0,jcplx1
      TQROC(K,N2,N1)=TQROC(K,N2,N1)+RQROC(K,N,NN,N2,N1)
      TQRON(K,N2,N1)=TQRON(K,N2,N1)+RQRON(K,N,NN,N2,N1)
      TQROP(K,N2,N1)=TQROP(K,N2,N1)+RQROP(K,N,NN,N2,N1)
      TQROA(K,N2,N1)=TQROA(K,N2,N1)+RQROA(K,N,NN,N2,N1)
9550  CONTINUE
      TQRCOS(N2,N1)=TQRCOS(N2,N1)+RQRCOS(N,NN,N2,N1)
      TQRCHS(N2,N1)=TQRCHS(N2,N1)+RQRCHS(N,NN,N2,N1)
      TQROXS(N2,N1)=TQROXS(N2,N1)+RQROXS(N,NN,N2,N1)
      TQRNGS(N2,N1)=TQRNGS(N2,N1)+RQRNGS(N,NN,N2,N1)
      TQRN2S(N2,N1)=TQRN2S(N2,N1)+RQRN2S(N,NN,N2,N1)
      TQRHGS(N2,N1)=TQRHGS(N2,N1)+RQRHGS(N,NN,N2,N1)
      TQRNH4(N2,N1)=TQRNH4(N2,N1)+RQRNH4(N,NN,N2,N1)
      TQRNH3(N2,N1)=TQRNH3(N2,N1)+RQRNH3(N,NN,N2,N1)
      TQRNO3(N2,N1)=TQRNO3(N2,N1)+RQRNO3(N,NN,N2,N1)
      TQRNO2(N2,N1)=TQRNO2(N2,N1)+RQRNO2(N,NN,N2,N1)
      TQRH1P(N2,N1)=TQRH1P(N2,N1)+RQRH1P(N,NN,N2,N1)
      TQRH2P(N2,N1)=TQRH2P(N2,N1)+RQRH2P(N,NN,N2,N1)
      IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      DO 9551 K=0,jcplx1
      TQROC(K,N2,N1)=TQROC(K,N2,N1)-RQROC(K,N,NN,N5,N4)
      TQRON(K,N2,N1)=TQRON(K,N2,N1)-RQRON(K,N,NN,N5,N4)
      TQROP(K,N2,N1)=TQROP(K,N2,N1)-RQROP(K,N,NN,N5,N4)
      TQROA(K,N2,N1)=TQROA(K,N2,N1)-RQROA(K,N,NN,N5,N4)
9551  CONTINUE
      TQRCOS(N2,N1)=TQRCOS(N2,N1)-RQRCOS(N,NN,N5,N4)
      TQRCHS(N2,N1)=TQRCHS(N2,N1)-RQRCHS(N,NN,N5,N4)
      TQROXS(N2,N1)=TQROXS(N2,N1)-RQROXS(N,NN,N5,N4)
      TQRNGS(N2,N1)=TQRNGS(N2,N1)-RQRNGS(N,NN,N5,N4)
      TQRN2S(N2,N1)=TQRN2S(N2,N1)-RQRN2S(N,NN,N5,N4)
      TQRHGS(N2,N1)=TQRHGS(N2,N1)-RQRHGS(N,NN,N5,N4)
      TQRNH4(N2,N1)=TQRNH4(N2,N1)-RQRNH4(N,NN,N5,N4)
      TQRNH3(N2,N1)=TQRNH3(N2,N1)-RQRNH3(N,NN,N5,N4)
      TQRNO3(N2,N1)=TQRNO3(N2,N1)-RQRNO3(N,NN,N5,N4)
      TQRNO2(N2,N1)=TQRNO2(N2,N1)-RQRNO2(N,NN,N5,N4)
      TQRH1P(N2,N1)=TQRH1P(N2,N1)-RQRH1P(N,NN,N5,N4)
      TQRH2P(N2,N1)=TQRH2P(N2,N1)-RQRH2P(N,NN,N5,N4)
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO 9552 K=0,jcplx1
      TQROC(K,N2,N1)=TQROC(K,N2,N1)-RQROC(K,N,NN,N5B,N4B)
      TQRON(K,N2,N1)=TQRON(K,N2,N1)-RQRON(K,N,NN,N5B,N4B)
      TQROP(K,N2,N1)=TQROP(K,N2,N1)-RQROP(K,N,NN,N5B,N4B)
      TQROA(K,N2,N1)=TQROA(K,N2,N1)-RQROA(K,N,NN,N5B,N4B)
9552  CONTINUE
      TQRCOS(N2,N1)=TQRCOS(N2,N1)-RQRCOS(N,NN,N5B,N4B)
      TQRCHS(N2,N1)=TQRCHS(N2,N1)-RQRCHS(N,NN,N5B,N4B)
      TQROXS(N2,N1)=TQROXS(N2,N1)-RQROXS(N,NN,N5B,N4B)
      TQRNGS(N2,N1)=TQRNGS(N2,N1)-RQRNGS(N,NN,N5B,N4B)
      TQRN2S(N2,N1)=TQRN2S(N2,N1)-RQRN2S(N,NN,N5B,N4B)
      TQRHGS(N2,N1)=TQRHGS(N2,N1)-RQRHGS(N,NN,N5B,N4B)
      TQRNH4(N2,N1)=TQRNH4(N2,N1)-RQRNH4(N,NN,N5B,N4B)
      TQRNH3(N2,N1)=TQRNH3(N2,N1)-RQRNH3(N,NN,N5B,N4B)
      TQRNO3(N2,N1)=TQRNO3(N2,N1)-RQRNO3(N,NN,N5B,N4B)
      TQRNO2(N2,N1)=TQRNO2(N2,N1)-RQRNO2(N,NN,N5B,N4B)
      TQRH1P(N2,N1)=TQRH1P(N2,N1)-RQRH1P(N,NN,N5B,N4B)
      TQRH2P(N2,N1)=TQRH2P(N2,N1)-RQRH2P(N,NN,N5B,N4B)
      ENDIF
1202  CONTINUE
!
!     NET OVERLAND SOLUTE FLUX IN SNOW
!
!     TQS*=net solute flux in snow transfer
!     RQS*=solute flux in snow transfer
!
      TQSCOS(N2,N1)=TQSCOS(N2,N1)+RQSCOS(N,N2,N1)-RQSCOS(N,N5,N4)
      TQSCHS(N2,N1)=TQSCHS(N2,N1)+RQSCHS(N,N2,N1)-RQSCHS(N,N5,N4)
      TQSOXS(N2,N1)=TQSOXS(N2,N1)+RQSOXS(N,N2,N1)-RQSOXS(N,N5,N4)
      TQSNGS(N2,N1)=TQSNGS(N2,N1)+RQSNGS(N,N2,N1)-RQSNGS(N,N5,N4)
      TQSN2S(N2,N1)=TQSN2S(N2,N1)+RQSN2S(N,N2,N1)-RQSN2S(N,N5,N4)
      TQSNH4(N2,N1)=TQSNH4(N2,N1)+RQSNH4(N,N2,N1)-RQSNH4(N,N5,N4)
      TQSNH3(N2,N1)=TQSNH3(N2,N1)+RQSNH3(N,N2,N1)-RQSNH3(N,N5,N4)
      TQSNO3(N2,N1)=TQSNO3(N2,N1)+RQSNO3(N,N2,N1)-RQSNO3(N,N5,N4)
      TQSH1P(N2,N1)=TQSH1P(N2,N1)+RQSH1P(N,N2,N1)-RQSH1P(N,N5,N4)
      TQSH2P(N2,N1)=TQSH2P(N2,N1)+RQSH2P(N,N2,N1)-RQSH2P(N,N5,N4)
!
!     NET SOLUTE FLUX IN SNOWPACK
!
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     T*BLS=net solute flux in snowpack
!     R*BLS=solute flux in snowpack
!
      ELSEIF(N.EQ.3)THEN
      DO 1205 LS=1,JS
      IF(VHCPWM(M,LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS.AND.VHCPWM(M,LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
      TCOBLS(LS,N2,N1)=TCOBLS(LS,N2,N1)+RCOBLS(LS,N2,N1) &
      -RCOBLS(LS2,N2,N1)
      TCHBLS(LS,N2,N1)=TCHBLS(LS,N2,N1)+RCHBLS(LS,N2,N1) &
      -RCHBLS(LS2,N2,N1)
      TOXBLS(LS,N2,N1)=TOXBLS(LS,N2,N1)+ROXBLS(LS,N2,N1) &
      -ROXBLS(LS2,N2,N1)
      TNGBLS(LS,N2,N1)=TNGBLS(LS,N2,N1)+RNGBLS(LS,N2,N1) &
      -RNGBLS(LS2,N2,N1)
      TN2BLS(LS,N2,N1)=TN2BLS(LS,N2,N1)+RN2BLS(LS,N2,N1) &
      -RN2BLS(LS2,N2,N1)
      TN4BLW(LS,N2,N1)=TN4BLW(LS,N2,N1)+RN4BLW(LS,N2,N1) &
      -RN4BLW(LS2,N2,N1)
      TN3BLW(LS,N2,N1)=TN3BLW(LS,N2,N1)+RN3BLW(LS,N2,N1) &
      -RN3BLW(LS2,N2,N1)
      TNOBLW(LS,N2,N1)=TNOBLW(LS,N2,N1)+RNOBLW(LS,N2,N1) &
      -RNOBLW(LS2,N2,N1)
      TH1PBS(LS,N2,N1)=TH1PBS(LS,N2,N1)+RH1PBS(LS,N2,N1) &
      -RH1PBS(LS2,N2,N1)
      TH2PBS(LS,N2,N1)=TH2PBS(LS,N2,N1)+RH2PBS(LS,N2,N1) &
      -RH2PBS(LS2,N2,N1)
      ELSE
!
!     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
!
      TCOBLS(LS,N2,N1)=TCOBLS(LS,N2,N1)+RCOBLS(LS,N2,N1) &
      -RCOFLS(3,0,N2,N1)-RCOFLS(3,NUM(N2,N1),N2,N1)
!    3-RCOFHS(3,NUM(N2,N1),N2,N1)
      TCHBLS(LS,N2,N1)=TCHBLS(LS,N2,N1)+RCHBLS(LS,N2,N1) &
      -RCHFLS(3,0,N2,N1)-RCHFLS(3,NUM(N2,N1),N2,N1)
!    3-RCHFHS(3,NUM(N2,N1),N2,N1)
      TOXBLS(LS,N2,N1)=TOXBLS(LS,N2,N1)+ROXBLS(LS,N2,N1) &
      -ROXFLS(3,0,N2,N1)-ROXFLS(3,NUM(N2,N1),N2,N1)
!    3-ROXFHS(3,NUM(N2,N1),N2,N1)
      TNGBLS(LS,N2,N1)=TNGBLS(LS,N2,N1)+RNGBLS(LS,N2,N1) &
      -RNGFLS(3,0,N2,N1)-RNGFLS(3,NUM(N2,N1),N2,N1)
!    3-RNGFHS(3,NUM(N2,N1),N2,N1)
      TN2BLS(LS,N2,N1)=TN2BLS(LS,N2,N1)+RN2BLS(LS,N2,N1) &
      -RN2FLS(3,0,N2,N1)-RN2FLS(3,NUM(N2,N1),N2,N1)
!    3-RN2FHS(3,NUM(N2,N1),N2,N1)
      TN4BLW(LS,N2,N1)=TN4BLW(LS,N2,N1)+RN4BLW(LS,N2,N1) &
      -RN4FLW(3,0,N2,N1)-RN4FLW(3,NUM(N2,N1),N2,N1) &
      -RN4FLB(3,NUM(N2,N1),N2,N1)
!    3-RN4FHW(3,NUM(N2,N1),N2,N1)-RN4FHB(3,NUM(N2,N1),N2,N1)
      TN3BLW(LS,N2,N1)=TN3BLW(LS,N2,N1)+RN3BLW(LS,N2,N1) &
      -RN3FLW(3,0,N2,N1)-RN3FLW(3,NUM(N2,N1),N2,N1) &
      -RN3FLB(3,NUM(N2,N1),N2,N1)
!    3-RN3FHW(3,NUM(N2,N1),N2,N1)-RN3FHB(3,NUM(N2,N1),N2,N1)
      TNOBLW(LS,N2,N1)=TNOBLW(LS,N2,N1)+RNOBLW(LS,N2,N1) &
      -RNOFLW(3,0,N2,N1)-RNOFLW(3,NUM(N2,N1),N2,N1) &
      -RNOFLB(3,NUM(N2,N1),N2,N1)
!    3-RNOFHW(3,NUM(N2,N1),N2,N1)-RNOFHB(3,NUM(N2,N1),N2,N1)
      TH1PBS(LS,N2,N1)=TH1PBS(LS,N2,N1)+RH1PBS(LS,N2,N1) &
      -RH1PFS(3,0,N2,N1)-RH1PFS(3,NUM(N2,N1),N2,N1) &
      -RH1BFB(3,NUM(N2,N1),N2,N1)
!    3-RH1PHS(3,NUM(N2,N1),N2,N1)-RH1BHB(3,NUM(N2,N1),N2,N1)
      TH2PBS(LS,N2,N1)=TH2PBS(LS,N2,N1)+RH2PBS(LS,N2,N1) &
      -RH2PFS(3,0,N2,N1)-RH2PFS(3,NUM(N2,N1),N2,N1) &
      -RH2BFB(3,NUM(N2,N1),N2,N1)
!    3-RH2PHS(3,NUM(N2,N1),N2,N1)-RH2BHB(3,NUM(N2,N1),N2,N1)
      ENDIF
      ENDIF
1205  CONTINUE
      ENDIF
      ENDIF
      ENDIF
      end subroutine NetOverlandFlux
!------------------------------------------------------------------------------------------

  subroutine NetFluxMicroandMacropores(NY,NX,N,M,MX,N1,N2,N3,N4,N5,N6)
  implicit none

  integer, intent(in) :: NY,NX,N,M,N1,N2,N3,N4,N5,MX
  integer, intent(inout) :: N6
  integer :: K,LL

  DO LL=N6,NL(NY,NX)
    IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
      N6=LL
      exit
    ENDIF
  ENDDO

  IF(M.NE.MX)THEN
    IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      DO 9545 K=0,jcplx1
      TOCFLS(K,N3,N2,N1)=TOCFLS(K,N3,N2,N1)+ROCFLS(K,N,N3,N2,N1) &
      -ROCFLS(K,N,N6,N5,N4)
      TONFLS(K,N3,N2,N1)=TONFLS(K,N3,N2,N1)+RONFLS(K,N,N3,N2,N1) &
      -RONFLS(K,N,N6,N5,N4)
      TOPFLS(K,N3,N2,N1)=TOPFLS(K,N3,N2,N1)+ROPFLS(K,N,N3,N2,N1) &
      -ROPFLS(K,N,N6,N5,N4)
      TOAFLS(K,N3,N2,N1)=TOAFLS(K,N3,N2,N1)+ROAFLS(K,N,N3,N2,N1) &
      -ROAFLS(K,N,N6,N5,N4)
      TOCFHS(K,N3,N2,N1)=TOCFHS(K,N3,N2,N1)+ROCFHS(K,N,N3,N2,N1) &
      -ROCFHS(K,N,N6,N5,N4)
      TONFHS(K,N3,N2,N1)=TONFHS(K,N3,N2,N1)+RONFHS(K,N,N3,N2,N1) &
      -RONFHS(K,N,N6,N5,N4)
      TOPFHS(K,N3,N2,N1)=TOPFHS(K,N3,N2,N1)+ROPFHS(K,N,N3,N2,N1) &
      -ROPFHS(K,N,N6,N5,N4)
      TOAFHS(K,N3,N2,N1)=TOAFHS(K,N3,N2,N1)+ROAFHS(K,N,N3,N2,N1) &
      -ROAFHS(K,N,N6,N5,N4)
9545  CONTINUE
      TCOFLS(N3,N2,N1)=TCOFLS(N3,N2,N1)+RCOFLS(N,N3,N2,N1) &
      -RCOFLS(N,N6,N5,N4)
!     IF(I.EQ.361.AND.NY.EQ.2)THEN
!     WRITE(*,3378)'TCOFLS',I,J,M,MM,N1,N2,N3,N4,N5,N6,N
!    2,TCOFLS(N3,N2,N1),RCOFLS(N,N3,N2,N1)
!    2,RCOFLS(N,N6,N5,N4),XCOFLS(N,N3,N2,N1)
!    2,XCOFLS(N,N6,N5,N4)
!3378  FORMAT(A8,11I4,20E12.4)
!     ENDIF
      TCHFLS(N3,N2,N1)=TCHFLS(N3,N2,N1)+RCHFLS(N,N3,N2,N1) &
      -RCHFLS(N,N6,N5,N4)
      TOXFLS(N3,N2,N1)=TOXFLS(N3,N2,N1)+ROXFLS(N,N3,N2,N1) &
      -ROXFLS(N,N6,N5,N4)
      TNGFLS(N3,N2,N1)=TNGFLS(N3,N2,N1)+RNGFLS(N,N3,N2,N1) &
      -RNGFLS(N,N6,N5,N4)
      TN2FLS(N3,N2,N1)=TN2FLS(N3,N2,N1)+RN2FLS(N,N3,N2,N1) &
      -RN2FLS(N,N6,N5,N4)
      THGFLS(N3,N2,N1)=THGFLS(N3,N2,N1)+RHGFLS(N,N3,N2,N1) &
      -RHGFLS(N,N6,N5,N4)
      TN4FLW(N3,N2,N1)=TN4FLW(N3,N2,N1)+RN4FLW(N,N3,N2,N1) &
      -RN4FLW(N,N6,N5,N4)
      TN3FLW(N3,N2,N1)=TN3FLW(N3,N2,N1)+RN3FLW(N,N3,N2,N1) &
      -RN3FLW(N,N6,N5,N4)
      TNOFLW(N3,N2,N1)=TNOFLW(N3,N2,N1)+RNOFLW(N,N3,N2,N1) &
      -RNOFLW(N,N6,N5,N4)
      TNXFLS(N3,N2,N1)=TNXFLS(N3,N2,N1)+RNXFLS(N,N3,N2,N1) &
      -RNXFLS(N,N6,N5,N4)
      TH1PFS(N3,N2,N1)=TH1PFS(N3,N2,N1)+RH1PFS(N,N3,N2,N1) &
      -RH1PFS(N,N6,N5,N4)
      TH2PFS(N3,N2,N1)=TH2PFS(N3,N2,N1)+RH2PFS(N,N3,N2,N1) &
      -RH2PFS(N,N6,N5,N4)
      TN4FLB(N3,N2,N1)=TN4FLB(N3,N2,N1)+RN4FLB(N,N3,N2,N1) &
      -RN4FLB(N,N6,N5,N4)
      TN3FLB(N3,N2,N1)=TN3FLB(N3,N2,N1)+RN3FLB(N,N3,N2,N1) &
      -RN3FLB(N,N6,N5,N4)
      TNOFLB(N3,N2,N1)=TNOFLB(N3,N2,N1)+RNOFLB(N,N3,N2,N1) &
      -RNOFLB(N,N6,N5,N4)
      TNXFLB(N3,N2,N1)=TNXFLB(N3,N2,N1)+RNXFLB(N,N3,N2,N1) &
      -RNXFLB(N,N6,N5,N4)
      TH1BFB(N3,N2,N1)=TH1BFB(N3,N2,N1)+RH1BFB(N,N3,N2,N1) &
      -RH1BFB(N,N6,N5,N4)
      TH2BFB(N3,N2,N1)=TH2BFB(N3,N2,N1)+RH2BFB(N,N3,N2,N1) &
      -RH2BFB(N,N6,N5,N4)
      TCOFHS(N3,N2,N1)=TCOFHS(N3,N2,N1)+RCOFHS(N,N3,N2,N1) &
      -RCOFHS(N,N6,N5,N4)
      TCHFHS(N3,N2,N1)=TCHFHS(N3,N2,N1)+RCHFHS(N,N3,N2,N1) &
      -RCHFHS(N,N6,N5,N4)
      TOXFHS(N3,N2,N1)=TOXFHS(N3,N2,N1)+ROXFHS(N,N3,N2,N1) &
      -ROXFHS(N,N6,N5,N4)
      TNGFHS(N3,N2,N1)=TNGFHS(N3,N2,N1)+RNGFHS(N,N3,N2,N1) &
      -RNGFHS(N,N6,N5,N4)
      TN2FHS(N3,N2,N1)=TN2FHS(N3,N2,N1)+RN2FHS(N,N3,N2,N1) &
      -RN2FHS(N,N6,N5,N4)
      THGFHS(N3,N2,N1)=THGFHS(N3,N2,N1)+RHGFHS(N,N3,N2,N1) &
      -RHGFHS(N,N6,N5,N4)
      TN4FHW(N3,N2,N1)=TN4FHW(N3,N2,N1)+RN4FHW(N,N3,N2,N1) &
      -RN4FHW(N,N6,N5,N4)
      TN3FHW(N3,N2,N1)=TN3FHW(N3,N2,N1)+RN3FHW(N,N3,N2,N1) &
      -RN3FHW(N,N6,N5,N4)
      TNOFHW(N3,N2,N1)=TNOFHW(N3,N2,N1)+RNOFHW(N,N3,N2,N1) &
      -RNOFHW(N,N6,N5,N4)
      TNXFHS(N3,N2,N1)=TNXFHS(N3,N2,N1)+RNXFHS(N,N3,N2,N1) &
      -RNXFHS(N,N6,N5,N4)
      TH1PHS(N3,N2,N1)=TH1PHS(N3,N2,N1)+RH1PHS(N,N3,N2,N1) &
      -RH1PHS(N,N6,N5,N4)
      TH2PHS(N3,N2,N1)=TH2PHS(N3,N2,N1)+RH2PHS(N,N3,N2,N1) &
      -RH2PHS(N,N6,N5,N4)
      TN4FHB(N3,N2,N1)=TN4FHB(N3,N2,N1)+RN4FHB(N,N3,N2,N1) &
      -RN4FHB(N,N6,N5,N4)
      TN3FHB(N3,N2,N1)=TN3FHB(N3,N2,N1)+RN3FHB(N,N3,N2,N1) &
      -RN3FHB(N,N6,N5,N4)
      TNOFHB(N3,N2,N1)=TNOFHB(N3,N2,N1)+RNOFHB(N,N3,N2,N1) &
      -RNOFHB(N,N6,N5,N4)
      TNXFHB(N3,N2,N1)=TNXFHB(N3,N2,N1)+RNXFHB(N,N3,N2,N1) &
      -RNXFHB(N,N6,N5,N4)
      TH1BHB(N3,N2,N1)=TH1BHB(N3,N2,N1)+RH1BHB(N,N3,N2,N1) &
      -RH1BHB(N,N6,N5,N4)
      TH2BHB(N3,N2,N1)=TH2BHB(N3,N2,N1)+RH2BHB(N,N3,N2,N1) &
      -RH2BHB(N,N6,N5,N4)
      ELSE
      DO 9546 K=0,jcplx1
      TOCFLS(K,N3,N2,N1)=0.0
      TONFLS(K,N3,N2,N1)=0.0
      TOPFLS(K,N3,N2,N1)=0.0
      TOAFLS(K,N3,N2,N1)=0.0
      TOCFHS(K,N3,N2,N1)=0.0
      TONFHS(K,N3,N2,N1)=0.0
      TOPFHS(K,N3,N2,N1)=0.0
      TOAFHS(K,N3,N2,N1)=0.0
9546  CONTINUE
      TCOFLS(N3,N2,N1)=0.0
      TCHFLS(N3,N2,N1)=0.0
      TOXFLS(N3,N2,N1)=0.0
      TNGFLS(N3,N2,N1)=0.0
      TN2FLS(N3,N2,N1)=0.0
      THGFLS(N3,N2,N1)=0.0
      TN4FLW(N3,N2,N1)=0.0
      TN3FLW(N3,N2,N1)=0.0
      TNOFLW(N3,N2,N1)=0.0
      TNXFLS(N3,N2,N1)=0.0
      TH1PFS(N3,N2,N1)=0.0
      TH2PFS(N3,N2,N1)=0.0
      TN4FLB(N3,N2,N1)=0.0
      TN3FLB(N3,N2,N1)=0.0
      TNOFLB(N3,N2,N1)=0.0
      TNXFLB(N3,N2,N1)=0.0
      TH1BFB(N3,N2,N1)=0.0
      TH2BFB(N3,N2,N1)=0.0
      TCOFHS(N3,N2,N1)=0.0
      TCHFHS(N3,N2,N1)=0.0
      TOXFHS(N3,N2,N1)=0.0
      TNGFHS(N3,N2,N1)=0.0
      TN2FHS(N3,N2,N1)=0.0
      THGFHS(N3,N2,N1)=0.0
      TN4FHW(N3,N2,N1)=0.0
      TN3FHW(N3,N2,N1)=0.0
      TNOFHW(N3,N2,N1)=0.0
      TNXFHS(N3,N2,N1)=0.0
      TH1PHS(N3,N2,N1)=0.0
      TH2PHS(N3,N2,N1)=0.0
      TN4FHB(N3,N2,N1)=0.0
      TN3FHB(N3,N2,N1)=0.0
      TNOFHB(N3,N2,N1)=0.0
      TNXFHB(N3,N2,N1)=0.0
      TH1BHB(N3,N2,N1)=0.0
      TH2BHB(N3,N2,N1)=0.0
      ENDIF
      ENDIF
!
!     NET GAS FLUX
!
!     T*FLG=net convective+diffusive gas flux
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      TCOFLG(N3,N2,N1)=TCOFLG(N3,N2,N1)+RCOFLG(N,N3,N2,N1) &
      -RCOFLG(N,N6,N5,N4)
      TCHFLG(N3,N2,N1)=TCHFLG(N3,N2,N1)+RCHFLG(N,N3,N2,N1) &
      -RCHFLG(N,N6,N5,N4)
      TOXFLG(N3,N2,N1)=TOXFLG(N3,N2,N1)+ROXFLG(N,N3,N2,N1) &
      -ROXFLG(N,N6,N5,N4)
      TNGFLG(N3,N2,N1)=TNGFLG(N3,N2,N1)+RNGFLG(N,N3,N2,N1) &
      -RNGFLG(N,N6,N5,N4)
      TN2FLG(N3,N2,N1)=TN2FLG(N3,N2,N1)+RN2FLG(N,N3,N2,N1) &
      -RN2FLG(N,N6,N5,N4)
      TN3FLG(N3,N2,N1)=TN3FLG(N3,N2,N1)+RN3FLG(N,N3,N2,N1) &
      -RN3FLG(N,N6,N5,N4)
      THGFLG(N3,N2,N1)=THGFLG(N3,N2,N1)+RHGFLG(N,N3,N2,N1) &
      -RHGFLG(N,N6,N5,N4)
      ELSE
      TCOFLG(N3,N2,N1)=0.0
      TCHFLG(N3,N2,N1)=0.0
      TOXFLG(N3,N2,N1)=0.0
      TNGFLG(N3,N2,N1)=0.0
      TN2FLG(N3,N2,N1)=0.0
      TN3FLG(N3,N2,N1)=0.0
      THGFLG(N3,N2,N1)=0.0
      ENDIF
      end subroutine NetFluxMicroandMacropores

end module BoundaryTranspMod