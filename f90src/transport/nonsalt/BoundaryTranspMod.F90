module BoundaryTranspMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use GridConsts
  USE SoilPropertyDataType
  use TranspNoSaltDataMod
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
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

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
        N1=NX;N2=NY;N3=L
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
        D9580: DO  N=FlowDirIndicator(NY,NX),3
          D9575: DO  NN=1,2
            IF(N.EQ.1)THEN
              !WEST-EAST
              N4=NX+1
              N5=NY
              N4B=NX-1
              N5B=NY
              N6=L
              IF(NN.EQ.1)THEN
                !eastern boundary
                IF(NX.EQ.NHE)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX+1
                  M5=NY
                  M6=L
                  XN=-1.0
                  RCHQF=RechargEastSurf(M2,M1)
                  RCHGFU=RechargEastSubSurf(M2,M1)
                  RCHGFT=RechargRateEastWTBL(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                !western boundary
                IF(NX.EQ.NHW)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0
                  RCHQF=RechargWestSurf(M5,M4)
                  RCHGFU=RechargWestSubSurf(M5,M4)
                  RCHGFT=RechargRateWestWTBL(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.2)THEN
              !NORTH-SOUTH
              N4=NX
              N5=NY+1
              N4B=NX
              N5B=NY-1
              N6=L
              IF(NN.EQ.1)THEN
              ! southern boundary
                IF(NY.EQ.NVS)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY+1
                  M6=L
                  XN=-1.0_r8
                  RCHQF=RechargSouthSurf(M2,M1)
                  RCHGFU=RechargSouthSubSurf(M2,M1)
                  RCHGFT=RechargRateSouthWTBL(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
              ! northern boundary
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0_r8
                  RCHQF=RechargNorthSurf(M5,M4)
                  RCHGFU=RechargNorthSubSurf(M5,M4)
                  RCHGFT=RechargRateNorthWTBL(M5,M4)
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
              !lower boundary
                IF(L.EQ.NL(NY,NX))THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L+1
                  XN=-1.0_r8
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                !nothing for the upper boundary
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
          IF(FlowDirIndicator(N2,N1).NE.3.OR.N.EQ.3)THEN
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
  integer :: K,NTG,NTN,NTS

  IF(.not.XGridRunoffFlag(NN,N,N2,N1).OR.isclose(RCHQF,0.0_r8).OR.WatFlux4ErosionM(M,N2,N1).LE.ZEROS(N2,N1))THEN
    DO  K=1,jcplx
      RQROC(K,N,NN,M5,M4)=0.0_r8
      RQRON(K,N,NN,M5,M4)=0.0_r8
      RQROP(K,N,NN,M5,M4)=0.0_r8
      RQROA(K,N,NN,M5,M4)=0.0_r8
    ENDDO
    trcg_RQR(idg_beg:idg_end-1,N,NN,M5,M4)=0.0_r8
    trcn_RQR(ids_nut_beg:ids_nuts_end,N,NN,M5,M4)=0.0_r8
  ELSE
!
!     SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
    IF((NN.EQ.1.AND.QflxSurfRunoffM(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.2.AND.QflxSurfRunoffM(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      FQRM=QflxSurfRunoffM(M,N,NN,M5,M4)/WatFlux4ErosionM(M,N2,N1)
      DO  K=1,jcplx
        RQROC(K,N,NN,M5,M4)=RQROC0(K,N2,N1)*FQRM
        RQRON(K,N,NN,M5,M4)=RQRON0(K,N2,N1)*FQRM
        RQROP(K,N,NN,M5,M4)=RQROP0(K,N2,N1)*FQRM
        RQROA(K,N,NN,M5,M4)=RQROA0(K,N2,N1)*FQRM
      enddo

      DO NTG=idg_beg,idg_end-1
        trcg_RQR(NTG,N,NN,M5,M4)=trcg_RQR0(NTG,N2,N1)*FQRM
      ENDDO

      DO NTS=ids_nut_beg,ids_nuts_end
        trcn_RQR(NTS,N,NN,M5,M4)=trcn_RQR0(NTS,N2,N1)*FQRM
      ENDDO
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
      DO NTG=idg_beg,idg_end-1
        trcg_XRS(NTG,N,NN,M5,M4)=trcg_XRS(NTG,N,NN,M5,M4)+trcg_RQR(NTG,N,NN,M5,M4)
      ENDDO

      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_XRS(NTN,N,NN,M5,M4)=trcn_XRS(NTN,N,NN,M5,M4)+trcn_RQR(NTN,N,NN,M5,M4)
      ENDDO
!
!     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
    ELSEIF((NN.EQ.2.AND.QflxSurfRunoffM(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.1.AND.QflxSurfRunoffM(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      DO  K=1,jcplx
        RQROC(K,N,NN,M5,M4)=0.0_r8
        RQRON(K,N,NN,M5,M4)=0.0_r8
        RQROP(K,N,NN,M5,M4)=0.0_r8
        RQROA(K,N,NN,M5,M4)=0.0_r8
      enddo

      trcg_RQR(idg_CO2,N,NN,M5,M4)=QflxSurfRunoffM(M,N,NN,M5,M4)*CCOU
      trcg_RQR(idg_CH4,N,NN,M5,M4)=QflxSurfRunoffM(M,N,NN,M5,M4)*CCHU
      trcg_RQR(idg_O2,N,NN,M5,M4)=QflxSurfRunoffM(M,N,NN,M5,M4)*COXU
      trcg_RQR(idg_N2,N,NN,M5,M4)=QflxSurfRunoffM(M,N,NN,M5,M4)*CNNU
      trcg_RQR(idg_N2O,N,NN,M5,M4)=QflxSurfRunoffM(M,N,NN,M5,M4)*CN2U
      trcg_RQR(idg_H2,N,NN,M5,M4)=0.0_r8
      trcg_RQR(idg_NH3,N,NN,M5,M4)=0.0_r8

      trcn_RQR(ids_nut_beg:ids_nuts_end,N,NN,M5,M4)=0.0_r8

      DO NTG=idg_beg,idg_end-1
        trcg_XRS(NTG,N,NN,M5,M4)=trcg_XRS(NTG,N,NN,M5,M4)+trcg_RQR(NTG,N,NN,M5,M4)
      ENDDO
    ELSE
      DO  K=1,jcplx
        RQROC(K,N,NN,M5,M4)=0.0_r8
        RQRON(K,N,NN,M5,M4)=0.0_r8
        RQROP(K,N,NN,M5,M4)=0.0_r8
        RQROA(K,N,NN,M5,M4)=0.0_r8
      enddo
      trcg_RQR(idg_beg:idg_end-1,N,NN,M5,M4)=0.0_r8
      trcn_RQR(ids_nut_beg:ids_nuts_end,N,NN,M5,M4)=0.0_r8
    ENDIF

  ENDIF
!
!     BOUNDARY SNOW FLUX
!
  IF(NN.EQ.1)THEN
    trcg_RQS(idg_beg:idg_end-1,N,M5,M4)=0.0_r8
    RQSNH4(N,M5,M4)=0.0_r8
    RQSNO3(N,M5,M4)=0.0_r8
    RQSH1P(N,M5,M4)=0.0_r8
    RQSH2P(N,M5,M4)=0.0_r8
  ENDIF
  end subroutine BoundaryRunoffandSnowXY
! ----------------------------------------------------------------------

  subroutine BoundaryRunoffandSnowZ(M,N,NN,M1,M2,M3,M4,M5,M6)
  implicit none
  integer, intent(in) :: M,N,NN,M1,M2,M3,M4,M5,M6

  integer :: K,NTG,NTN,NTS
  real(r8) :: VFLW

  IF(NN.EQ.1.AND.WaterFlow2MicPM(M,N,M6,M5,M4).GT.0.0_r8 &
    .OR.NN.EQ.2.AND.WaterFlow2MicPM(M,N,M6,M5,M4).LT.0.0_r8)THEN
    IF(VLWatMicPM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MicPM(M,N,M6,M5,M4)/VLWatMicPM(M,M3,M2,M1)))
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
    R3PoreSolFlx(idg_beg:idg_end-2,N,M6,M5,M4)=0.0_r8

!add irrigation flux
    DO NTS=ids_nuts_beg,ids_nuts_end
      R3PoreSolFlx(NTS,N,M6,M5,M4)=WaterFlow2MicPM(M,N,M6,M5,M4) &
        *trcn_irrig(NTS,M3,M2,M1)*trcs_VLN(NTS,M3,M2,M1)
    ENDDO

  ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
!
!     WaterFlow2MacPM=water flux through soil macropore from watsub.f
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     RFH*S=solute diffusive flux through macropore
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  IF(NN.EQ.1.AND.WaterFlow2MacPM(M,N,M6,M5,M4).GT.0.0 &
    .OR.NN.EQ.2.AND.WaterFlow2MacPM(M,N,M6,M5,M4).LT.0.0)THEN
    IF(VLWatMacPM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MacPM(M,N,M6,M5,M4)/VLWatMacPM(M,M3,M2,M1)))
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

  DO NTS=ids_beg,ids_end
    trcs_3DTransp2MicP(NTS,N,M6,M5,M4)=trcs_3DTransp2MicP(NTS,N,M6,M5,M4)+R3PoreSolFlx(NTS,N,M6,M5,M4)
    trcs_3DTransp2MacP(NTS,N,M6,M5,M4)=trcs_3DTransp2MacP(NTS,N,M6,M5,M4)+R3PoreSoHFlx(NTS,N,M6,M5,M4)
  ENDDO


  end subroutine BoundaryRunoffandSnowZ

! ----------------------------------------------------------------------
  subroutine BoundaryRunoffandSnow(L,N,NN,M,MX,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
  implicit none

  integer, intent(in) :: L,N, NN, M,MX, N1, N2,N3, M1, M2,M3,M4, M5,M6
  real(r8), intent(in):: RCHQF
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K,NTG

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

!     WaterFlow2MicPM=water flux through soil micropore from watsub.f
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     R*FLS=convective solute flux through micropores
!     R*FLW,R*FLB=convective solute flux through micropores in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
    IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(FlowDirIndicator(M2,M1).NE.3.OR.N.EQ.3)THEN

        call BoundaryRunoffandSnowZ(M,N,NN,M1,M2,M3,M4,M5,M6)
      ENDIF
    ENDIF
!
!     GASOUS LOSS WITH SUBSURFACE MICROPORE WATER GAIN
!
!     WaterFlow2MicPM,WaterFlow2MacPM=micropore,macropore water flux from watsub.f
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!     VLsoiAirPM=air-filled porosity
!     R*FLG=convective gas flux
!     X*FLG=convective gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
    FLGM=(WaterFlow2MicPM(M,N,M6,M5,M4)+WaterFlow2MacPM(M,N,M6,M5,M4))*XNPT
    IF(NN.EQ.1.AND.FLGM.LT.0.0.OR.NN.EQ.2.AND.FLGM.GT.0.0)THEN
      IF(VLsoiAirPM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
        VFLW=-AMAX1(-VFLWX,AMIN1(VFLWX,FLGM/VLsoiAirPM(M,M3,M2,M1)))
      ELSE
        VFLW=0.0_r8
      ENDIF
!
!     HOURLY GAS FLUX FOR USE IN REDIST.F
!
!     X*FLG=hourly convective gas flux
!
      DO NTG=idg_beg,idg_end-1
        R3GasADFlx(NTG,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2(NTG,M3,M2,M1))
        R3GasADTFlx(NTG,N,M6,M5,M4)=R3GasADTFlx(NTG,N,M6,M5,M4)+R3GasADFlx(NTG,N,M6,M5,M4)
      ENDDO

    ELSE
      R3GasADFlx(idg_beg:idg_end-1,N,M6,M5,M4)=0.0_r8
    ENDIF
  ELSE
      DO  K=1,jcplx
        ROCFHS(K,N,M6,M5,M4)=0.0_r8
        RONFHS(K,N,M6,M5,M4)=0.0_r8
        ROPFHS(K,N,M6,M5,M4)=0.0_r8
        ROAFHS(K,N,M6,M5,M4)=0.0_r8
      enddo
      R3PoreSoHFlx(ids_beg:ids_end,N,M6,M5,M4)=0.0_r8
      R3GasADFlx(idg_beg:idg_end-1,N,M6,M5,M4)=0.0_r8
  ENDIF
  end subroutine BoundaryRunoffandSnow
!------------------------------------------------------------------------------------------

  subroutine NetOverlandFluxXY(M,N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N4,N5,N4B,N5B

  integer :: NN,K,NTG,NTS

  DO NN=1,2
    DO  K=1,jcplx
      TQROC(K,N2,N1)=TQROC(K,N2,N1)+RQROC(K,N,NN,N2,N1)
      TQRON(K,N2,N1)=TQRON(K,N2,N1)+RQRON(K,N,NN,N2,N1)
      TQROP(K,N2,N1)=TQROP(K,N2,N1)+RQROP(K,N,NN,N2,N1)
      TQROA(K,N2,N1)=TQROA(K,N2,N1)+RQROA(K,N,NN,N2,N1)
    enddo

    DO NTG=idg_beg,idg_end-1
      trcg_TQR(NTG,N2,N1)=trcg_TQR(NTG,N2,N1)+trcg_RQR(NTG,N,NN,N2,N1)
    ENDDO

    DO NTS=ids_nut_beg,ids_nuts_end
      trcn_TQR(NTS,N2,N1)=trcn_TQR(NTS,N2,N1)+trcn_RQR(NTS,N,NN,N2,N1)
    ENDDO

    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      DO  K=1,jcplx
        TQROC(K,N2,N1)=TQROC(K,N2,N1)-RQROC(K,N,NN,N5,N4)
        TQRON(K,N2,N1)=TQRON(K,N2,N1)-RQRON(K,N,NN,N5,N4)
        TQROP(K,N2,N1)=TQROP(K,N2,N1)-RQROP(K,N,NN,N5,N4)
        TQROA(K,N2,N1)=TQROA(K,N2,N1)-RQROA(K,N,NN,N5,N4)
      enddo

      DO NTG=idg_beg,idg_end-1
        trcg_TQR(NTG,N2,N1)=trcg_TQR(NTG,N2,N1)-trcg_RQR(NTG,N,NN,N5,N4)
      ENDDO

      DO NTS=ids_nut_beg,ids_nuts_end
        trcn_TQR(NTS,N2,N1)=trcn_TQR(NTS,N2,N1)-trcn_RQR(NTS,N,NN,N5,N4)
      ENDDO
    ENDIF
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO  K=1,jcplx
        TQROC(K,N2,N1)=TQROC(K,N2,N1)-RQROC(K,N,NN,N5B,N4B)
        TQRON(K,N2,N1)=TQRON(K,N2,N1)-RQRON(K,N,NN,N5B,N4B)
        TQROP(K,N2,N1)=TQROP(K,N2,N1)-RQROP(K,N,NN,N5B,N4B)
        TQROA(K,N2,N1)=TQROA(K,N2,N1)-RQROA(K,N,NN,N5B,N4B)
      enddo
      DO NTG=idg_beg,idg_end-1
        trcg_TQR(NTG,N2,N1)=trcg_TQR(NTG,N2,N1)-trcg_RQR(NTG,N,NN,N5B,N4B)
      ENDDO

      DO NTS=ids_nut_beg,ids_nuts_end
        trcn_TQR(NTS,N2,N1)=trcn_TQR(NTS,N2,N1)-trcn_RQR(NTS,N,NN,N5B,N4B)
      ENDDO
    ENDIF
  enddo
!
!     NET OVERLAND SOLUTE FLUX IN SNOW
!
!     TQS*=net solute flux in snow transfer
!     RQS*=solute flux in snow transfer
!
  DO NTG=idg_beg,idg_end-1
    trcg_TQ(NTG,N2,N1)=trcg_TQ(NTG,N2,N1)+trcg_RQS(NTG,N,N2,N1)-trcg_RQS(NTG,N,N5,N4)
  ENDDO

  trcn_TQ(ids_NH4,N2,N1)=trcn_TQ(ids_NH4,N2,N1)+RQSNH4(N,N2,N1)-RQSNH4(N,N5,N4)
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
!     VLSnowHeatCapM,VLHeatCapSnowMin=current,minimum volumetric heat capacity of snowpack
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
  integer :: NTG,NTN

  DO  LS=1,JS
    IF(VLSnowHeatCapM(M,LS,NY,NX).GT.VLHeatCapSnowMin(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS.AND.VLSnowHeatCapM(M,LS2,N2,N1).GT.VLHeatCapSnowMin(N2,N1))THEN
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1) &
            +trcg_RBLS(NTG,LS,N2,N1)-trcg_RBLS(NTG,LS2,N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1) &
            +trcn_RBLS(NTN,LS,N2,N1)-trcn_RBLS(NTN,LS2,N2,N1)
        ENDDO
      ELSE
!
!     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
!
    ! exclude NH3 and NH3B
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1) &
            +trcg_RBLS(NTG,LS,N2,N1)-R3PoreSolFlx(NTG,3,0,N2,N1) &
            -R3PoreSolFlx(NTG,3,NUM(N2,N1),N2,N1)
  !    3-R3PoreSoHFlx(NTG,3,NUM(N2,N1),N2,N1)
        ENDDO

        trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1) &
          -R3PoreSolFlx(idg_NH3B,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx(idg_NH3B,3,NUM(N2,N1),N2,N1)

        DO NTN=0,ids_nuts
          trcn_TBLS(ids_NH4+NTN,LS,N2,N1)=trcn_TBLS(ids_NH4+NTN,LS,N2,N1) &
            +trcn_RBLS(ids_NH4+NTN,LS,N2,N1)-R3PoreSolFlx(ids_NH4+NTN,3,0,N2,N1) &
            -R3PoreSolFlx(ids_NH4+NTN,3,NUM(N2,N1),N2,N1) &
            -R3PoreSolFlx(ids_NH4B+NTN,3,NUM(N2,N1),N2,N1)
!    -R3PoreSoHFlx(ids_NH4+NTN,3,NUM(N2,N1),N2,N1) &
!    -R3PoreSoHFlx(ids_NH4B+NTN,3,NUM(N2,N1),N2,N1)
        ENDDO
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
    IF(VLSoilPoreMicP(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
      N6=LL
      exit
    ENDIF
  ENDDO

  IF(M.NE.MX)THEN
    IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
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
        R3PorTSolFlx(NTS,N3,N2,N1)=R3PorTSolFlx(NTS,N3,N2,N1) &
          +R3PoreSolFlx(NTS,N,N3,N2,N1)-R3PoreSolFlx(NTS,N,N6,N5,N4)
        R3PorTSoHFlx(NTS,N3,N2,N1)=R3PorTSoHFlx(NTS,N3,N2,N1) &
          +R3PoreSoHFlx(NTS,N,N3,N2,N1)-R3PoreSoHFlx(NTS,N,N6,N5,N4)
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
  IF(VLSoilPoreMicP(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    DO NTG=idg_beg,idg_end-1
      RTGasADFlx(NTG,N3,N2,N1)=RTGasADFlx(NTG,N3,N2,N1) &
        +R3GasADFlx(NTG,N,N3,N2,N1)-R3GasADFlx(NTG,N,N6,N5,N4)
    ENDDO
  ELSE
    RTGasADFlx(idg_beg:idg_end-1,N3,N2,N1)=0._r8
  ENDIF
  end subroutine NetFluxMicroandMacropores

end module BoundaryTranspMod
