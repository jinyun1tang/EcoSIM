module BoundaryTranspMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb,AZMAX1
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
      D9585: DO L=NU(NY,NX),NL(NY,NX)
        N1=NX
        N2=NY
        N3=L
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
        D9580: DO  N=NCN(NY,NX),3
          D9575: DO  NN=1,2
            IF(N.EQ.1)THEN
!direction x
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
!direction y
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
!vertical
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
          ENDDO D9575
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
        ENDDO D9580
      ENDDO D9585
    ENDDO
  ENDDO
  end subroutine BoundaryFlux
!------------------------------------------------------------------------------------------

  subroutine BoundaryRunoffandSnowXY(M,N,NN,N1,N2,M4,M5,RCHQF)
  implicit none
  integer, intent(in) :: M,N,NN,N1,N2,M4,M5
  real(r8), intent(in) :: RCHQF
  real(r8) :: FQRM
  integer :: K

  IF(IRCHG(NN,N,N2,N1).EQ.0.OR.test_aeqb(RCHQF,0.0_r8).OR.QRM(M,N2,N1).LE.ZEROS(N2,N1))THEN
    DO  K=1,jcplx
      RQROC(K,N,NN,M5,M4)=0.0_r8
      RQRON(K,N,NN,M5,M4)=0.0_r8
      RQROP(K,N,NN,M5,M4)=0.0_r8
      RQROA(K,N,NN,M5,M4)=0.0_r8
    ENDDO
    RQRCOS(N,NN,M5,M4)=0.0_r8
    RQRCHS(N,NN,M5,M4)=0.0_r8
    RQROXS(N,NN,M5,M4)=0.0_r8
    RQRNGS(N,NN,M5,M4)=0.0_r8
    RQRN2S(N,NN,M5,M4)=0.0_r8
    RQRHGS(N,NN,M5,M4)=0.0_r8
    RQRNH4(N,NN,M5,M4)=0.0_r8
    RQRNH3(N,NN,M5,M4)=0.0_r8
    RQRNO3(N,NN,M5,M4)=0.0_r8
    RQRNO2(N,NN,M5,M4)=0.0_r8
    RQRH1P(N,NN,M5,M4)=0.0_r8
    RQRH2P(N,NN,M5,M4)=0.0_r8
  ELSE
!
!     SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
    IF((NN.EQ.1.AND.QRMN(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.2.AND.QRMN(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      FQRM=QRMN(M,N,NN,M5,M4)/QRM(M,N2,N1)
      DO  K=1,jcplx
        RQROC(K,N,NN,M5,M4)=RQROC0(K,N2,N1)*FQRM
        RQRON(K,N,NN,M5,M4)=RQRON0(K,N2,N1)*FQRM
        RQROP(K,N,NN,M5,M4)=RQROP0(K,N2,N1)*FQRM
        RQROA(K,N,NN,M5,M4)=RQROA0(K,N2,N1)*FQRM
      enddo
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
      DO  K=1,jcplx
        XOCQRS(K,N,NN,M5,M4)=XOCQRS(K,N,NN,M5,M4)+RQROC(K,N,NN,M5,M4)
        XONQRS(K,N,NN,M5,M4)=XONQRS(K,N,NN,M5,M4)+RQRON(K,N,NN,M5,M4)
        XOPQRS(K,N,NN,M5,M4)=XOPQRS(K,N,NN,M5,M4)+RQROP(K,N,NN,M5,M4)
        XOAQRS(K,N,NN,M5,M4)=XOAQRS(K,N,NN,M5,M4)+RQROA(K,N,NN,M5,M4)
      enddo
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
      DO  K=1,jcplx
        RQROC(K,N,NN,M5,M4)=0.0_r8
        RQRON(K,N,NN,M5,M4)=0.0_r8
        RQROP(K,N,NN,M5,M4)=0.0_r8
        RQROA(K,N,NN,M5,M4)=0.0_r8
      enddo
      RQRCOS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CCOU
      RQRCHS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CCHU
      RQROXS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*COXU
      RQRNGS(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CNNU
      RQRN2S(N,NN,M5,M4)=QRMN(M,N,NN,M5,M4)*CN2U
      RQRHGS(N,NN,M5,M4)=0.0_r8
      RQRNH4(N,NN,M5,M4)=0.0_r8
      RQRNH3(N,NN,M5,M4)=0.0_r8
      RQRNO3(N,NN,M5,M4)=0.0_r8
      RQRNO2(N,NN,M5,M4)=0.0_r8
      RQRH1P(N,NN,M5,M4)=0.0_r8
      RQRH2P(N,NN,M5,M4)=0.0_r8
      XCOQRS(N,NN,M5,M4)=XCOQRS(N,NN,M5,M4)+RQRCOS(N,NN,M5,M4)
      XCHQRS(N,NN,M5,M4)=XCHQRS(N,NN,M5,M4)+RQRCHS(N,NN,M5,M4)
      XOXQRS(N,NN,M5,M4)=XOXQRS(N,NN,M5,M4)+RQROXS(N,NN,M5,M4)
      XNGQRS(N,NN,M5,M4)=XNGQRS(N,NN,M5,M4)+RQRNGS(N,NN,M5,M4)
      XN2QRS(N,NN,M5,M4)=XN2QRS(N,NN,M5,M4)+RQRN2S(N,NN,M5,M4)
    ELSE
      DO  K=1,jcplx
        RQROC(K,N,NN,M5,M4)=0.0_r8
        RQRON(K,N,NN,M5,M4)=0.0_r8
        RQROP(K,N,NN,M5,M4)=0.0_r8
        RQROA(K,N,NN,M5,M4)=0.0_r8
      enddo
      RQRCOS(N,NN,M5,M4)=0.0_r8
      RQRCHS(N,NN,M5,M4)=0.0_r8
      RQROXS(N,NN,M5,M4)=0.0_r8
      RQRNGS(N,NN,M5,M4)=0.0_r8
      RQRN2S(N,NN,M5,M4)=0.0_r8
      RQRHGS(N,NN,M5,M4)=0.0_r8
      RQRNH4(N,NN,M5,M4)=0.0_r8
      RQRNH3(N,NN,M5,M4)=0.0_r8
      RQRNO3(N,NN,M5,M4)=0.0_r8
      RQRNO2(N,NN,M5,M4)=0.0_r8
      RQRH1P(N,NN,M5,M4)=0.0_r8
      RQRH2P(N,NN,M5,M4)=0.0_r8
    ENDIF

  ENDIF
!
!     BOUNDARY SNOW FLUX
!
  IF(NN.EQ.1)THEN
    RQSCOS(N,M5,M4)=0.0_r8
    RQSCHS(N,M5,M4)=0.0_r8
    RQSOXS(N,M5,M4)=0.0_r8
    RQSNGS(N,M5,M4)=0.0_r8
    RQSN2S(N,M5,M4)=0.0_r8
    RQSNH4(N,M5,M4)=0.0_r8
    RQSNH3(N,M5,M4)=0.0_r8
    RQSNO3(N,M5,M4)=0.0_r8
    RQSH1P(N,M5,M4)=0.0_r8
    RQSH2P(N,M5,M4)=0.0_r8
  ENDIF
  end subroutine BoundaryRunoffandSnowXY
! ----------------------------------------------------------------------

  subroutine BoundaryRunoffandSnowZ(M,N,NN,M1,M2,M3,M4,M5,M6)
  implicit none
  integer, intent(in) :: M,N,NN,M1,M2,M3,M4,M5,M6

  integer :: K,NTG,NTN
  real(r8) :: VFLW

  IF(NN.EQ.1.AND.FLWM(M,N,M6,M5,M4).GT.0.0_r8 &
    .OR.NN.EQ.2.AND.FLWM(M,N,M6,M5,M4).LT.0.0_r8)THEN
    IF(VOLWM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,FLWM(M,N,M6,M5,M4)/VOLWM(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
    DO  K=1,jcplx
      ROCFLS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQC2(K,M3,M2,M1))
      RONFLS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQN2(K,M3,M2,M1))
      ROPFLS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQP2(K,M3,M2,M1))
      ROAFLS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQA2(K,M3,M2,M1))
    enddo
    !does not include NH3 and NH3B
    DO NTG=idg_beg,idg_end-2
      R3PoreSolFlx(NTG,N,M6,M5,M4)=VFLW*AZMAX1(trc_solml2(NTG,M3,M2,M1))
    ENDDo

    DO NTN=ids_nuts_beg,ids_nuts_end
      R3PoreSolFlx(NTN,N,M6,M5,M4)=VFLW*AZMAX1(trc_solml2(NTN,M3,M2,M1))*trcs_VLN(NTN,M3,M2,M1)
    ENDDO
!
!     SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
!
  ELSE
    DO  K=1,jcplx
      ROCFLS(K,N,M6,M5,M4)=0.0_r8
      RONFLS(K,N,M6,M5,M4)=0.0_r8
      ROPFLS(K,N,M6,M5,M4)=0.0_r8
      ROAFLS(K,N,M6,M5,M4)=0.0_r8
    enddo
    R3PoreSolFlx(idg_CO2,N,M6,M5,M4)=0.0_r8
    R3PoreSolFlx(idg_CH4,N,M6,M5,M4)=0.0_r8
    R3PoreSolFlx(idg_O2,N,M6,M5,M4)=0.0_r8
    R3PoreSolFlx(idg_N2,N,M6,M5,M4)=0.0_r8
    R3PoreSolFlx(idg_N2O,N,M6,M5,M4)=0.0_r8
    R3PoreSolFlx(idg_H2,N,M6,M5,M4)=0.0_r8

    R3PoreSolFlx(ids_NO2,N,M6,M5,M4)=0.0_r8
    R3PoreSolFlx(ids_NO2B,N,M6,M5,M4)=0.0_r8
!add irrigation flux
    R3PoreSolFlx(ids_NH4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN4U(M3,M2,M1)*trcs_VLN(ids_NH4,M3,M2,M1)
    R3PoreSolFlx(idg_NH3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN3U(M3,M2,M1)*trcs_VLN(idg_NH3,M3,M2,M1)
    R3PoreSolFlx(ids_NO3,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNOU(M3,M2,M1)*trcs_VLN(ids_NO3,M3,M2,M1)
    R3PoreSolFlx(ids_H1PO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4,M3,M2,M1)
    R3PoreSolFlx(ids_H2PO4,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH2PU(M3,M2,M1)*trcs_VLN(ids_H2PO4,M3,M2,M1)
    R3PoreSolFlx(ids_NH4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN4U(M3,M2,M1)*trcs_VLN(ids_NH4B,M3,M2,M1)
    R3PoreSolFlx(idg_NH3B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CN3U(M3,M2,M1)*trcs_VLN(idg_NH3B,M3,M2,M1)
    R3PoreSolFlx(ids_NO3B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CNOU(M3,M2,M1)*trcs_VLN(ids_NO3B,M3,M2,M1)
    R3PoreSolFlx(ids_H1PO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH1PU(M3,M2,M1)*trcs_VLN(ids_H1PO4B,M3,M2,M1)
    R3PoreSolFlx(ids_H2PO4B,N,M6,M5,M4)=FLWM(M,N,M6,M5,M4)*CH2PU(M3,M2,M1)*trcs_VLN(ids_H2PO4B,M3,M2,M1)
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
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,FLWHM(M,N,M6,M5,M4)/VOLWHM(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
    DO  K=1,jcplx
      ROCFHS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQCH2(K,M3,M2,M1))
      RONFHS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQNH2(K,M3,M2,M1))
      ROPFHS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQPH2(K,M3,M2,M1))
      ROAFHS(K,N,M6,M5,M4)=VFLW*AZMAX1(OQAH2(K,M3,M2,M1))
    enddo

    DO NTG=idg_beg,idg_end-2
      R3PoreSoHFlx(NTG,N,M6,M5,M4)=VFLW*AZMAX1(trc_soHml2(NTG,M3,M2,M1))
    ENDDO

    DO NTN=ids_nuts_beg,ids_nuts_end
      R3PoreSoHFlx(NTN,N,M6,M5,M4)=VFLW*AZMAX1(trc_soHml2(NTN,M3,M2,M1))*trcs_VLN(NTN,M3,M2,M1)
    ENDDO
!
!     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
!
  ELSE
    DO  K=1,jcplx
      ROCFHS(K,N,M6,M5,M4)=0.0_r8
      RONFHS(K,N,M6,M5,M4)=0.0_r8
      ROPFHS(K,N,M6,M5,M4)=0.0_r8
      ROAFHS(K,N,M6,M5,M4)=0.0_r8
    enddo
    R3PoreSoHFlx(ids_beg:ids_end,N,M6,M5,M4)=0.0_r8
  ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS,X*FLW,X*FLB=hourly solute flux in non-band,band micropores
!     X*FHS,X*FHW,X*FHB=hourly solute flux in non-band,band macropores
!     R*FLS,R*FLW,R*FLB=solute flux in non-band,band micropores
!     R*FHS,R*FHW,R*FHB=solute flux in non-band,band macropores
!
  DO  K=1,jcplx
    XOCFLS(K,N,M6,M5,M4)=XOCFLS(K,N,M6,M5,M4)+ROCFLS(K,N,M6,M5,M4)
    XONFLS(K,N,M6,M5,M4)=XONFLS(K,N,M6,M5,M4)+RONFLS(K,N,M6,M5,M4)
    XOPFLS(K,N,M6,M5,M4)=XOPFLS(K,N,M6,M5,M4)+ROPFLS(K,N,M6,M5,M4)
    XOAFLS(K,N,M6,M5,M4)=XOAFLS(K,N,M6,M5,M4)+ROAFLS(K,N,M6,M5,M4)
    XOCFHS(K,N,M6,M5,M4)=XOCFHS(K,N,M6,M5,M4)+ROCFHS(K,N,M6,M5,M4)
    XONFHS(K,N,M6,M5,M4)=XONFHS(K,N,M6,M5,M4)+RONFHS(K,N,M6,M5,M4)
    XOPFHS(K,N,M6,M5,M4)=XOPFHS(K,N,M6,M5,M4)+ROPFHS(K,N,M6,M5,M4)
    XOAFHS(K,N,M6,M5,M4)=XOAFHS(K,N,M6,M5,M4)+ROAFHS(K,N,M6,M5,M4)
  enddo

  trcs_XFLS(idg_CO2,N,M6,M5,M4)=trcs_XFLS(idg_CO2,N,M6,M5,M4)+R3PoreSolFlx(idg_CO2,N,M6,M5,M4)
  trcs_XFLS(idg_CH4,N,M6,M5,M4)=trcs_XFLS(idg_CH4,N,M6,M5,M4)+R3PoreSolFlx(idg_CH4,N,M6,M5,M4)
  trcs_XFLS(idg_O2,N,M6,M5,M4)=trcs_XFLS(idg_O2,N,M6,M5,M4)+R3PoreSolFlx(idg_O2,N,M6,M5,M4)
  trcs_XFLS(idg_N2,N,M6,M5,M4)=trcs_XFLS(idg_N2,N,M6,M5,M4)+R3PoreSolFlx(idg_N2,N,M6,M5,M4)
  trcs_XFLS(idg_N2O,N,M6,M5,M4)=trcs_XFLS(idg_N2O,N,M6,M5,M4)+R3PoreSolFlx(idg_N2O,N,M6,M5,M4)
  trcs_XFLS(idg_H2,N,M6,M5,M4)=trcs_XFLS(idg_H2,N,M6,M5,M4)+R3PoreSolFlx(idg_H2,N,M6,M5,M4)
  trcs_XFLS(ids_NH4,N,M6,M5,M4)=trcs_XFLS(ids_NH4,N,M6,M5,M4)+R3PoreSolFlx(ids_NH4,N,M6,M5,M4)
  trcs_XFLS(idg_NH3,N,M6,M5,M4)=trcs_XFLS(idg_NH3,N,M6,M5,M4)+R3PoreSolFlx(idg_NH3,N,M6,M5,M4)
  trcs_XFLS(ids_NO3,N,M6,M5,M4)=trcs_XFLS(ids_NO3,N,M6,M5,M4)+R3PoreSolFlx(ids_NO3,N,M6,M5,M4)
  trcs_XFLS(ids_NO2,N,M6,M5,M4)=trcs_XFLS(ids_NO2,N,M6,M5,M4)+R3PoreSolFlx(ids_NO2,N,M6,M5,M4)
  trcs_XFLS(ids_H1PO4,N,M6,M5,M4)=trcs_XFLS(ids_H1PO4,N,M6,M5,M4)+R3PoreSolFlx(ids_H1PO4,N,M6,M5,M4)
  trcs_XFLS(ids_H2PO4,N,M6,M5,M4)=trcs_XFLS(ids_H2PO4,N,M6,M5,M4)+R3PoreSolFlx(ids_H2PO4,N,M6,M5,M4)
  trcs_XFLS(ids_NH4B,N,M6,M5,M4)=trcs_XFLS(ids_NH4B,N,M6,M5,M4)+R3PoreSolFlx(ids_NH4B,N,M6,M5,M4)
  trcs_XFLS(idg_NH3B,N,M6,M5,M4)=trcs_XFLS(idg_NH3B,N,M6,M5,M4)+R3PoreSolFlx(idg_NH3B,N,M6,M5,M4)
  trcs_XFLS(ids_NO3B,N,M6,M5,M4)=trcs_XFLS(ids_NO3B,N,M6,M5,M4)+R3PoreSolFlx(ids_NO3B,N,M6,M5,M4)
  trcs_XFLS(ids_NO2B,N,M6,M5,M4)=trcs_XFLS(ids_NO2B,N,M6,M5,M4)+R3PoreSolFlx(ids_NO2B,N,M6,M5,M4)
  trcs_XFLS(ids_H1PO4B,N,M6,M5,M4)=trcs_XFLS(ids_H1PO4B,N,M6,M5,M4)+R3PoreSolFlx(ids_H1PO4B,N,M6,M5,M4)
  trcs_XFLS(ids_H2PO4B,N,M6,M5,M4)=trcs_XFLS(ids_H2PO4B,N,M6,M5,M4)+R3PoreSolFlx(ids_H2PO4B,N,M6,M5,M4)

  XCOFHS(N,M6,M5,M4)=XCOFHS(N,M6,M5,M4)+R3PoreSoHFlx(idg_CO2,N,M6,M5,M4)
  XCHFHS(N,M6,M5,M4)=XCHFHS(N,M6,M5,M4)+R3PoreSoHFlx(idg_CH4,N,M6,M5,M4)
  XOXFHS(N,M6,M5,M4)=XOXFHS(N,M6,M5,M4)+R3PoreSoHFlx(idg_O2,N,M6,M5,M4)
  XNGFHS(N,M6,M5,M4)=XNGFHS(N,M6,M5,M4)+R3PoreSoHFlx(idg_N2,N,M6,M5,M4)
  XN2FHS(N,M6,M5,M4)=XN2FHS(N,M6,M5,M4)+R3PoreSoHFlx(idg_N2O,N,M6,M5,M4)
  XHGFHS(N,M6,M5,M4)=XHGFHS(N,M6,M5,M4)+R3PoreSoHFlx(idg_H2,N,M6,M5,M4)
  XN4FHW(N,M6,M5,M4)=XN4FHW(N,M6,M5,M4)+R3PoreSoHFlx(ids_NH4,N,M6,M5,M4)
  XN3FHW(N,M6,M5,M4)=XN3FHW(N,M6,M5,M4)+R3PoreSoHFlx(idg_NH3,N,M6,M5,M4)
  XNOFHW(N,M6,M5,M4)=XNOFHW(N,M6,M5,M4)+R3PoreSoHFlx(ids_NO3,N,M6,M5,M4)
  XNXFHS(N,M6,M5,M4)=XNXFHS(N,M6,M5,M4)+R3PoreSoHFlx(ids_NO2,N,M6,M5,M4)
  XH1PHS(N,M6,M5,M4)=XH1PHS(N,M6,M5,M4)+R3PoreSoHFlx(ids_H1PO4,N,M6,M5,M4)
  XH2PHS(N,M6,M5,M4)=XH2PHS(N,M6,M5,M4)+R3PoreSoHFlx(ids_H2PO4,N,M6,M5,M4)
  XN4FHB(N,M6,M5,M4)=XN4FHB(N,M6,M5,M4)+R3PoreSoHFlx(ids_NH4B,N,M6,M5,M4)
  XN3FHB(N,M6,M5,M4)=XN3FHB(N,M6,M5,M4)+R3PoreSoHFlx(idg_NH3B,N,M6,M5,M4)
  XNOFHB(N,M6,M5,M4)=XNOFHB(N,M6,M5,M4)+R3PoreSoHFlx(ids_NO3B,N,M6,M5,M4)
  XNXFHB(N,M6,M5,M4)=XNXFHB(N,M6,M5,M4)+R3PoreSoHFlx(ids_NO2B,N,M6,M5,M4)
  XH1BHB(N,M6,M5,M4)=XH1BHB(N,M6,M5,M4)+R3PoreSoHFlx(ids_H1PO4B,N,M6,M5,M4)
  XH2BHB(N,M6,M5,M4)=XH2BHB(N,M6,M5,M4)+R3PoreSoHFlx(ids_H2PO4B,N,M6,M5,M4)
  end subroutine BoundaryRunoffandSnowZ

! ----------------------------------------------------------------------
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
      call BoundaryRunoffandSnowXY(M,N,NN,N1,N2,M4,M5,RCHQF)
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

        call BoundaryRunoffandSnowZ(M,N,NN,M1,M2,M3,M4,M5,M6)
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
    IF(NN.EQ.1.AND.FLGM.LT.0.0.OR.NN.EQ.2.AND.FLGM.GT.0.0)THEN
      IF(VOLPM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
        VFLW=-AMAX1(-VFLWX,AMIN1(VFLWX,FLGM/VOLPM(M,M3,M2,M1)))
      ELSE
        VFLW=0.0_r8
      ENDIF
      R3GasADFlx(idg_CO2,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(idg_CO2,M3,M2,M1))
      R3GasADFlx(idg_CH4,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(idg_CH4,M3,M2,M1))
      R3GasADFlx(idg_O2,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(idg_O2,M3,M2,M1))
      R3GasADFlx(idg_N2,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(idg_N2,M3,M2,M1))
      R3GasADFlx(idg_N2O,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(idg_N2O,M3,M2,M1))
      R3GasADFlx(idg_NH3,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(idg_NH3,M3,M2,M1))
      R3GasADFlx(idg_H2,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(idg_H2,M3,M2,M1))
!
!     HOURLY GAS FLUX FOR USE IN REDIST.F
!
!     X*FLG=hourly convective gas flux
!
      R3GasADTFlx(idg_CO2,N,M6,M5,M4)=R3GasADTFlx(idg_CO2,N,M6,M5,M4)+R3GasADFlx(idg_CO2,N,M6,M5,M4)
      R3GasADTFlx(idg_CH4,N,M6,M5,M4)=R3GasADTFlx(idg_CH4,N,M6,M5,M4)+R3GasADFlx(idg_CH4,N,M6,M5,M4)
      R3GasADTFlx(idg_O2,N,M6,M5,M4)=R3GasADTFlx(idg_O2,N,M6,M5,M4)+R3GasADFlx(idg_O2,N,M6,M5,M4)
      R3GasADTFlx(idg_N2,N,M6,M5,M4)=R3GasADTFlx(idg_N2,N,M6,M5,M4)+R3GasADFlx(idg_N2,N,M6,M5,M4)
      R3GasADTFlx(idg_N2O,N,M6,M5,M4)=R3GasADTFlx(idg_N2O,N,M6,M5,M4)+R3GasADFlx(idg_N2O,N,M6,M5,M4)
      R3GasADTFlx(idg_NH3,N,M6,M5,M4)=R3GasADTFlx(idg_NH3,N,M6,M5,M4)+R3GasADFlx(idg_NH3,N,M6,M5,M4)
      R3GasADTFlx(idg_H2,N,M6,M5,M4)=R3GasADTFlx(idg_H2,N,M6,M5,M4)+R3GasADFlx(idg_H2,N,M6,M5,M4)
    ELSE
      R3GasADFlx(idg_CO2,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_CH4,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_O2,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_N2,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_N2O,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_NH3,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_H2,N,M6,M5,M4)=0.0_r8
    ENDIF
  ELSE
      DO  K=1,jcplx
        ROCFHS(K,N,M6,M5,M4)=0.0_r8
        RONFHS(K,N,M6,M5,M4)=0.0_r8
        ROPFHS(K,N,M6,M5,M4)=0.0_r8
        ROAFHS(K,N,M6,M5,M4)=0.0_r8
      enddo
      R3PoreSoHFlx(idg_CO2,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(idg_CH4,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(idg_O2,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(idg_N2,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(idg_N2O,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_NH4,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(idg_H2,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(idg_NH3,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_NO3,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_NO2,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_H1PO4,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_H2PO4,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_NH4B,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(idg_NH3B,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_NO3B,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_NO2B,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_H1PO4B,N,M6,M5,M4)=0.0_r8
      R3PoreSoHFlx(ids_H2PO4B,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_CO2,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_CH4,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_O2,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_N2,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_N2O,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_NH3,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_H2,N,M6,M5,M4)=0.0_r8
  ENDIF
  end subroutine BoundaryRunoffandSnow
!------------------------------------------------------------------------------------------

  subroutine NetOverlandFluxXY(M,N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N4,N5,N4B,N5B

  integer :: NN,K

  DO NN=1,2
    DO  K=1,jcplx
      TQROC(K,N2,N1)=TQROC(K,N2,N1)+RQROC(K,N,NN,N2,N1)
      TQRON(K,N2,N1)=TQRON(K,N2,N1)+RQRON(K,N,NN,N2,N1)
      TQROP(K,N2,N1)=TQROP(K,N2,N1)+RQROP(K,N,NN,N2,N1)
      TQROA(K,N2,N1)=TQROA(K,N2,N1)+RQROA(K,N,NN,N2,N1)
    enddo
    trcg_TQR(idg_CO2,N2,N1)=trcg_TQR(idg_CO2,N2,N1)+RQRCOS(N,NN,N2,N1)
    trcg_TQR(idg_CH4,N2,N1)=trcg_TQR(idg_CH4,N2,N1)+RQRCHS(N,NN,N2,N1)
    trcg_TQR(idg_O2,N2,N1)=trcg_TQR(idg_O2,N2,N1)+RQROXS(N,NN,N2,N1)
    trcg_TQR(idg_N2,N2,N1)=trcg_TQR(idg_N2,N2,N1)+RQRNGS(N,NN,N2,N1)
    trcg_TQR(idg_N2O,N2,N1)=trcg_TQR(idg_N2O,N2,N1)+RQRN2S(N,NN,N2,N1)
    trcg_TQR(idg_H2,N2,N1)=trcg_TQR(idg_H2,N2,N1)+RQRHGS(N,NN,N2,N1)
    trcn_TQR(ids_NH4,N2,N1)=trcn_TQR(ids_NH4,N2,N1)+RQRNH4(N,NN,N2,N1)
    trcg_TQR(idg_NH3,N2,N1)=trcg_TQR(idg_NH3,N2,N1)+RQRNH3(N,NN,N2,N1)
    trcn_TQR(ids_NO3,N2,N1)=trcn_TQR(ids_NO3,N2,N1)+RQRNO3(N,NN,N2,N1)
    trcn_TQR(ids_NO2,N2,N1)=trcn_TQR(ids_NO2,N2,N1)+RQRNO2(N,NN,N2,N1)
    trcn_TQR(ids_H1PO4,N2,N1)=trcn_TQR(ids_H1PO4,N2,N1)+RQRH1P(N,NN,N2,N1)
    trcn_TQR(ids_H2PO4,N2,N1)=trcn_TQR(ids_H2PO4,N2,N1)+RQRH2P(N,NN,N2,N1)
    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      DO  K=1,jcplx
        TQROC(K,N2,N1)=TQROC(K,N2,N1)-RQROC(K,N,NN,N5,N4)
        TQRON(K,N2,N1)=TQRON(K,N2,N1)-RQRON(K,N,NN,N5,N4)
        TQROP(K,N2,N1)=TQROP(K,N2,N1)-RQROP(K,N,NN,N5,N4)
        TQROA(K,N2,N1)=TQROA(K,N2,N1)-RQROA(K,N,NN,N5,N4)
      enddo
      trcg_TQR(idg_CO2,N2,N1)=trcg_TQR(idg_CO2,N2,N1)-RQRCOS(N,NN,N5,N4)
      trcg_TQR(idg_CH4,N2,N1)=trcg_TQR(idg_CH4,N2,N1)-RQRCHS(N,NN,N5,N4)
      trcg_TQR(idg_O2,N2,N1)=trcg_TQR(idg_O2,N2,N1)-RQROXS(N,NN,N5,N4)
      trcg_TQR(idg_N2,N2,N1)=trcg_TQR(idg_N2,N2,N1)-RQRNGS(N,NN,N5,N4)
      trcg_TQR(idg_N2O,N2,N1)=trcg_TQR(idg_N2O,N2,N1)-RQRN2S(N,NN,N5,N4)
      trcg_TQR(idg_H2,N2,N1)=trcg_TQR(idg_H2,N2,N1)-RQRHGS(N,NN,N5,N4)
      trcn_TQR(ids_NH4,N2,N1)=trcn_TQR(ids_NH4,N2,N1)-RQRNH4(N,NN,N5,N4)
      trcg_TQR(idg_NH3,N2,N1)=trcg_TQR(idg_NH3,N2,N1)-RQRNH3(N,NN,N5,N4)
      trcn_TQR(ids_NO3,N2,N1)=trcn_TQR(ids_NO3,N2,N1)-RQRNO3(N,NN,N5,N4)
      trcn_TQR(ids_NO2,N2,N1)=trcn_TQR(ids_NO2,N2,N1)-RQRNO2(N,NN,N5,N4)
      trcn_TQR(ids_H1PO4,N2,N1)=trcn_TQR(ids_H1PO4,N2,N1)-RQRH1P(N,NN,N5,N4)
      trcn_TQR(ids_H2PO4,N2,N1)=trcn_TQR(ids_H2PO4,N2,N1)-RQRH2P(N,NN,N5,N4)
    ENDIF
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO  K=1,jcplx
        TQROC(K,N2,N1)=TQROC(K,N2,N1)-RQROC(K,N,NN,N5B,N4B)
        TQRON(K,N2,N1)=TQRON(K,N2,N1)-RQRON(K,N,NN,N5B,N4B)
        TQROP(K,N2,N1)=TQROP(K,N2,N1)-RQROP(K,N,NN,N5B,N4B)
        TQROA(K,N2,N1)=TQROA(K,N2,N1)-RQROA(K,N,NN,N5B,N4B)
      enddo
      trcg_TQR(idg_CO2,N2,N1)=trcg_TQR(idg_CO2,N2,N1)-RQRCOS(N,NN,N5B,N4B)
      trcg_TQR(idg_CH4,N2,N1)=trcg_TQR(idg_CH4,N2,N1)-RQRCHS(N,NN,N5B,N4B)
      trcg_TQR(idg_O2,N2,N1)=trcg_TQR(idg_O2,N2,N1)-RQROXS(N,NN,N5B,N4B)
      trcg_TQR(idg_N2,N2,N1)=trcg_TQR(idg_N2,N2,N1)-RQRNGS(N,NN,N5B,N4B)
      trcg_TQR(idg_N2O,N2,N1)=trcg_TQR(idg_N2O,N2,N1)-RQRN2S(N,NN,N5B,N4B)
      trcg_TQR(idg_H2,N2,N1)=trcg_TQR(idg_H2,N2,N1)-RQRHGS(N,NN,N5B,N4B)
      trcn_TQR(ids_NH4,N2,N1)=trcn_TQR(ids_NH4,N2,N1)-RQRNH4(N,NN,N5B,N4B)
      trcg_TQR(idg_NH3,N2,N1)=trcg_TQR(idg_NH3,N2,N1)-RQRNH3(N,NN,N5B,N4B)
      trcn_TQR(ids_NO3,N2,N1)=trcn_TQR(ids_NO3,N2,N1)-RQRNO3(N,NN,N5B,N4B)
      trcn_TQR(ids_NO2,N2,N1)=trcn_TQR(ids_NO2,N2,N1)-RQRNO2(N,NN,N5B,N4B)
      trcn_TQR(ids_H1PO4,N2,N1)=trcn_TQR(ids_H1PO4,N2,N1)-RQRH1P(N,NN,N5B,N4B)
      trcn_TQR(ids_H2PO4,N2,N1)=trcn_TQR(ids_H2PO4,N2,N1)-RQRH2P(N,NN,N5B,N4B)
    ENDIF
  enddo
!
!     NET OVERLAND SOLUTE FLUX IN SNOW
!
!     TQS*=net solute flux in snow transfer
!     RQS*=solute flux in snow transfer
!
  trcg_TQ(idg_O2,N2,N1)=trcg_TQ(idg_O2,N2,N1)+RQSCOS(N,N2,N1)-RQSCOS(N,N5,N4)
  trcg_TQ(idg_CH4,N2,N1)=trcg_TQ(idg_CH4,N2,N1)+RQSCHS(N,N2,N1)-RQSCHS(N,N5,N4)
  trcg_TQ(idg_O2,N2,N1)=trcg_TQ(idg_O2,N2,N1)+RQSOXS(N,N2,N1)-RQSOXS(N,N5,N4)
  trcg_TQ(idg_N2,N2,N1)=trcg_TQ(idg_N2,N2,N1)+RQSNGS(N,N2,N1)-RQSNGS(N,N5,N4)
  trcg_TQ(idg_N2O,N2,N1)=trcg_TQ(idg_N2O,N2,N1)+RQSN2S(N,N2,N1)-RQSN2S(N,N5,N4)
  trcn_TQ(ids_NH4,N2,N1)=trcn_TQ(ids_NH4,N2,N1)+RQSNH4(N,N2,N1)-RQSNH4(N,N5,N4)
  trcg_TQ(idg_NH3,N2,N1)=trcg_TQ(idg_NH3,N2,N1)+RQSNH3(N,N2,N1)-RQSNH3(N,N5,N4)
  trcn_TQ(ids_NO3,N2,N1)=trcn_TQ(ids_NO3,N2,N1)+RQSNO3(N,N2,N1)-RQSNO3(N,N5,N4)
  trcn_TQ(ids_H1PO4,N2,N1)=trcn_TQ(ids_H1PO4,N2,N1)+RQSH1P(N,N2,N1)-RQSH1P(N,N5,N4)
  trcn_TQ(ids_H2PO4,N2,N1)=trcn_TQ(ids_H2PO4,N2,N1)+RQSH2P(N,N2,N1)-RQSH2P(N,N5,N4)
  end subroutine NetOverlandFluxXY

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
!horizontal flow
        call NetOverlandFluxXY(M,N,N1,N2,N4,N5,N4B,N5B)
!
!     NET SOLUTE FLUX IN SNOWPACK
!
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     T*BLS=net solute flux in snowpack
!     R*BLS=solute flux in snowpack
!
      ELSEIF(N.EQ.3)THEN
! vertical direction
        call NetOverlandFluxZ(M,N1,N2,NY,NX)
      ENDIF
    ENDIF
  ENDIF
  end subroutine NetOverlandFlux
!------------------------------------------------------------------------------------------

  subroutine NetOverlandFluxZ(M,N1,N2,NY,NX)
  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX

  integer :: LS,LS2

  DO  LS=1,JS
    IF(VHCPWM(M,LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS.AND.VHCPWM(M,LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
        trcg_TBLS(idg_CO2,LS,N2,N1)=trcg_TBLS(idg_CO2,LS,N2,N1)+RCOBLS(LS,N2,N1)-RCOBLS(LS2,N2,N1)
        trcg_TBLS(idg_CH4,LS,N2,N1)=trcg_TBLS(idg_CH4,LS,N2,N1)+RCHBLS(LS,N2,N1)-RCHBLS(LS2,N2,N1)
        trcg_TBLS(idg_O2,LS,N2,N1)=trcg_TBLS(idg_O2,LS,N2,N1)+ROXBLS(LS,N2,N1)-ROXBLS(LS2,N2,N1)
        trcg_TBLS(idg_N2,LS,N2,N1)=trcg_TBLS(idg_N2,LS,N2,N1)+RNGBLS(LS,N2,N1)-RNGBLS(LS2,N2,N1)
        trcg_TBLS(idg_N2O,LS,N2,N1)=trcg_TBLS(idg_N2O,LS,N2,N1)+RN2BLS(LS,N2,N1)-RN2BLS(LS2,N2,N1)
        trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1)+RN3BLW(LS,N2,N1)-RN3BLW(LS2,N2,N1)

        trcn_TBLS(ids_NH4,LS,N2,N1)=trcn_TBLS(ids_NH4,LS,N2,N1)+RN4BLW(LS,N2,N1)-RN4BLW(LS2,N2,N1)
        trcn_TBLS(ids_NO3,LS,N2,N1)=trcn_TBLS(ids_NO3,LS,N2,N1)+RNOBLW(LS,N2,N1)-RNOBLW(LS2,N2,N1)
        trcn_TBLS(ids_H1PO4,LS,N2,N1)=trcn_TBLS(ids_H1PO4,LS,N2,N1)+RH1PBS(LS,N2,N1)-RH1PBS(LS2,N2,N1)
        trcn_TBLS(ids_H2PO4,LS,N2,N1)=trcn_TBLS(ids_H2PO4,LS,N2,N1)+RH2PBS(LS,N2,N1)-RH2PBS(LS2,N2,N1)
      ELSE
!
!     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
!
        trcg_TBLS(idg_CO2,LS,N2,N1)=trcg_TBLS(idg_CO2,LS,N2,N1)+RCOBLS(LS,N2,N1)-R3PoreSolFlx(idg_CO2,3,0,N2,N1)-R3PoreSolFlx(idg_CO2,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(idg_CO2,3,NUM(N2,N1),N2,N1)
        trcg_TBLS(idg_CH4,LS,N2,N1)=trcg_TBLS(idg_CH4,LS,N2,N1)+RCHBLS(LS,N2,N1)-R3PoreSolFlx(idg_CH4,3,0,N2,N1)-R3PoreSolFlx(idg_CH4,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(idg_CH4,3,NUM(N2,N1),N2,N1)
        trcg_TBLS(idg_O2,LS,N2,N1)=trcg_TBLS(idg_O2,LS,N2,N1)+ROXBLS(LS,N2,N1)-R3PoreSolFlx(idg_O2,3,0,N2,N1)-R3PoreSolFlx(idg_O2,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(idg_O2,3,NUM(N2,N1),N2,N1)
        trcg_TBLS(idg_N2,LS,N2,N1)=trcg_TBLS(idg_N2,LS,N2,N1)+RNGBLS(LS,N2,N1)-R3PoreSolFlx(idg_N2,3,0,N2,N1)-R3PoreSolFlx(idg_N2,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(idg_N2,3,NUM(N2,N1),N2,N1)
      trcg_TBLS(idg_N2O,LS,N2,N1)=trcg_TBLS(idg_N2O,LS,N2,N1)+RN2BLS(LS,N2,N1)-R3PoreSolFlx(idg_N2O,3,0,N2,N1)-R3PoreSolFlx(idg_N2O,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(idg_N2O,3,NUM(N2,N1),N2,N1)
        trcn_TBLS(ids_NH4,LS,N2,N1)=trcn_TBLS(ids_NH4,LS,N2,N1)+RN4BLW(LS,N2,N1) &
          -R3PoreSolFlx(ids_NH4,3,0,N2,N1)-R3PoreSolFlx(ids_NH4,3,NUM(N2,N1),N2,N1) &
          -R3PoreSolFlx(ids_NH4B,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(ids_NH4,3,NUM(N2,N1),N2,N1)-R3PoreSoHFlx(ids_NH4B,3,NUM(N2,N1),N2,N1)
        trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1)+RN3BLW(LS,N2,N1) &
          -R3PoreSolFlx(idg_NH3,3,0,N2,N1)-R3PoreSolFlx(idg_NH3,3,NUM(N2,N1),N2,N1) &
          -R3PoreSolFlx(idg_NH3B,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(idg_NH3,3,NUM(N2,N1),N2,N1)-R3PoreSoHFlx(idg_NH3B,3,NUM(N2,N1),N2,N1)
        trcn_TBLS(ids_NO3,LS,N2,N1)=trcn_TBLS(ids_NO3,LS,N2,N1)+RNOBLW(LS,N2,N1) &
          -R3PoreSolFlx(ids_NO3,3,0,N2,N1)-R3PoreSolFlx(ids_NO3,3,NUM(N2,N1),N2,N1) &
          -R3PoreSolFlx(ids_NO3B,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(ids_NO3,3,NUM(N2,N1),N2,N1)-R3PoreSoHFlx(ids_NO3B,3,NUM(N2,N1),N2,N1)
        trcn_TBLS(ids_H1PO4,LS,N2,N1)=trcn_TBLS(ids_H1PO4,LS,N2,N1)+RH1PBS(LS,N2,N1) &
          -R3PoreSolFlx(ids_H1PO4,3,0,N2,N1)-R3PoreSolFlx(ids_H1PO4,3,NUM(N2,N1),N2,N1) &
          -R3PoreSolFlx(ids_H1PO4B,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(ids_H1PO4,3,NUM(N2,N1),N2,N1)-R3PoreSoHFlx(ids_H1PO4B,3,NUM(N2,N1),N2,N1)
        trcn_TBLS(ids_H2PO4,LS,N2,N1)=trcn_TBLS(ids_H2PO4,LS,N2,N1)+RH2PBS(LS,N2,N1) &
          -R3PoreSolFlx(ids_H2PO4,3,0,N2,N1)-R3PoreSolFlx(ids_H2PO4,3,NUM(N2,N1),N2,N1) &
          -R3PoreSolFlx(ids_H2PO4B,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(ids_H2PO4,3,NUM(N2,N1),N2,N1)-R3PoreSoHFlx(ids_H2PO4B,3,NUM(N2,N1),N2,N1)
      ENDIF
    ENDIF
  enddo
  end subroutine NetOverlandFluxZ
!------------------------------------------------------------------------------------------

  subroutine NetFluxMicroandMacropores(NY,NX,N,M,MX,N1,N2,N3,N4,N5,N6)
  implicit none

  integer, intent(in) :: NY,NX,N,M,N1,N2,N3,N4,N5,MX
  integer, intent(inout) :: N6
  integer :: K,LL,NTS,NTG

  DO LL=N6,NL(NY,NX)
    IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
      N6=LL
      exit
    ENDIF
  ENDDO

  IF(M.NE.MX)THEN
    IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      DO  K=1,jcplx
        TOCFLS(K,N3,N2,N1)=TOCFLS(K,N3,N2,N1)+ROCFLS(K,N,N3,N2,N1)-ROCFLS(K,N,N6,N5,N4)
        TONFLS(K,N3,N2,N1)=TONFLS(K,N3,N2,N1)+RONFLS(K,N,N3,N2,N1)-RONFLS(K,N,N6,N5,N4)
        TOPFLS(K,N3,N2,N1)=TOPFLS(K,N3,N2,N1)+ROPFLS(K,N,N3,N2,N1)-ROPFLS(K,N,N6,N5,N4)
        TOAFLS(K,N3,N2,N1)=TOAFLS(K,N3,N2,N1)+ROAFLS(K,N,N3,N2,N1)-ROAFLS(K,N,N6,N5,N4)
        TOCFHS(K,N3,N2,N1)=TOCFHS(K,N3,N2,N1)+ROCFHS(K,N,N3,N2,N1)-ROCFHS(K,N,N6,N5,N4)
        TONFHS(K,N3,N2,N1)=TONFHS(K,N3,N2,N1)+RONFHS(K,N,N3,N2,N1)-RONFHS(K,N,N6,N5,N4)
        TOPFHS(K,N3,N2,N1)=TOPFHS(K,N3,N2,N1)+ROPFHS(K,N,N3,N2,N1)-ROPFHS(K,N,N6,N5,N4)
        TOAFHS(K,N3,N2,N1)=TOAFHS(K,N3,N2,N1)+ROAFHS(K,N,N3,N2,N1)-ROAFHS(K,N,N6,N5,N4)
      enddo
      DO NTS=ids_beg,ids_end
        R3PorTSolFlx(NTS,N3,N2,N1)=R3PorTSolFlx(NTS,N3,N2,N1)+R3PoreSolFlx(NTS,N,N3,N2,N1)-R3PoreSolFlx(NTS,N,N6,N5,N4)
        R3PorTSoHFlx(NTS,N3,N2,N1)=R3PorTSoHFlx(NTS,N3,N2,N1)+R3PoreSoHFlx(NTS,N,N3,N2,N1)-R3PoreSoHFlx(NTS,N,N6,N5,N4)

      ENDDO

    ELSE
      DO  K=1,jcplx
        TOCFLS(K,N3,N2,N1)=0.0_r8
        TONFLS(K,N3,N2,N1)=0.0_r8
        TOPFLS(K,N3,N2,N1)=0.0_r8
        TOAFLS(K,N3,N2,N1)=0.0_r8
        TOCFHS(K,N3,N2,N1)=0.0_r8
        TONFHS(K,N3,N2,N1)=0.0_r8
        TOPFHS(K,N3,N2,N1)=0.0_r8
        TOAFHS(K,N3,N2,N1)=0.0_r8
      enddo
      R3PorTSolFlx(ids_beg:ids_end,N3,N2,N1)=0.0_r8
      R3PorTSoHFlx(ids_beg:ids_end,N3,N2,N1)=0._r8

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
    DO NTG=idg_beg,idg_end
      RTGasADFlx(NTG,N3,N2,N1)=RTGasADFlx(NTG,N3,N2,N1)+R3GasADFlx(NTG,N,N3,N2,N1)-R3GasADFlx(NTG,N,N6,N5,N4)
    ENDDO
  ELSE
    RTGasADFlx(idg_beg:idg_end,N3,N2,N1)=0._r8
  ENDIF
  end subroutine NetFluxMicroandMacropores

end module BoundaryTranspMod
