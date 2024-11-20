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
                  XN=-1.0    !going out
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
                  XN=1.0     !coming in
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
                  XN=-1.0_r8    !going out
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
                  XN=1.0_r8   !coming in
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
                  XN=-1.0_r8  !going out
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
  integer :: K,idg,NTN,ids,idom

  IF(.not.XGridRunoffFlag(NN,N,N2,N1).OR.isclose(RCHQF,0.0_r8).OR.WatFlux4ErosionM_2DH(M,N2,N1).LE.ZEROS(N2,N1))THEN
    DO  K=1,jcplx
      dom_2DFloXSurRunoffM(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
    ENDDO
    trcg_2DFloXSurRunoffM(idg_beg:idg_NH3,N,NN,M5,M4)          = 0.0_r8
    trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,N,NN,M5,M4) = 0.0_r8
  ELSE
!
!     SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
    IF((NN.EQ.1.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.2.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      FQRM=QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/WatFlux4ErosionM_2DH(M,N2,N1)
      DO  K=1,jcplx
        DO idom=idom_beg,idom_end
          dom_2DFloXSurRunoffM(idom,K,N,NN,M5,M4)=dom_FloXSurRunoff(idom,K,N2,N1)*FQRM
        ENDDO
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_2DFloXSurRunoffM(idg,N,NN,M5,M4)=trcg_FloXSurRunoff(idg,N2,N1)*FQRM
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_2DFloXSurRunoffM(ids,N,NN,M5,M4)=trcn_FloXSurRunoff(ids,N2,N1)*FQRM
      ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QRS=hourly solute in runoff
!     RQR*=solute in runoff
!
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          dom_2DFloXSurRunoff(idom,K,N,NN,M5,M4)=dom_2DFloXSurRunoff(idom,K,N,NN,M5,M4)+dom_2DFloXSurRunoffM(idom,K,N,NN,M5,M4)
        enddo
      enddo
      DO idg=idg_beg,idg_NH3
        trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)=trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)+trcg_2DFloXSurRunoffM(idg,N,NN,M5,M4)
      ENDDO

      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_FloXSurRunoff_2D(NTN,N,NN,M5,M4)=trcn_FloXSurRunoff_2D(NTN,N,NN,M5,M4)+trcn_2DFloXSurRunoffM(NTN,N,NN,M5,M4)
      ENDDO
!
!     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
    ELSEIF((NN.EQ.2.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.1.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      DO  K=1,jcplx
        dom_2DFloXSurRunoffM(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
      enddo

      trcg_2DFloXSurRunoffM(idg_CO2,N,NN,M5,M4) = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*CCOU
      trcg_2DFloXSurRunoffM(idg_CH4,N,NN,M5,M4) = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*CCHU
      trcg_2DFloXSurRunoffM(idg_O2,N,NN,M5,M4)  = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*COXU
      trcg_2DFloXSurRunoffM(idg_N2,N,NN,M5,M4)  = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*CNNU
      trcg_2DFloXSurRunoffM(idg_N2O,N,NN,M5,M4) = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*CN2U
      trcg_2DFloXSurRunoffM(idg_H2,N,NN,M5,M4)  = 0.0_r8
      trcg_2DFloXSurRunoffM(idg_NH3,N,NN,M5,M4) = 0.0_r8

      trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,N,NN,M5,M4)=0.0_r8

      DO idg=idg_beg,idg_NH3
        trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)=trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)+trcg_2DFloXSurRunoffM(idg,N,NN,M5,M4)
      ENDDO
    ELSE
      DO  K=1,jcplx
        dom_2DFloXSurRunoffM(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
      enddo
      trcg_2DFloXSurRunoffM(idg_beg:idg_NH3,N,NN,M5,M4)          = 0.0_r8
      trcn_2DFloXSurRunoffM(ids_nut_beg:ids_nuts_end,N,NN,M5,M4) = 0.0_r8
    ENDIF

  ENDIF
!
!     BOUNDARY SNOW FLUX
!
  IF(NN.EQ.1)THEN
    trcg_2DSnowDrift(idg_beg:idg_NH3,N,M5,M4) = 0.0_r8
    trcn_2DSnowDrift(ids_NH4,N,M5,M4)         = 0.0_r8
    trcn_2DSnowDrift(ids_NO3,N,M5,M4)         = 0.0_r8
    trcn_2DSnowDrift(ids_H1PO4,N,M5,M4)       = 0.0_r8
    trcn_2DSnowDrift(ids_H2PO4,N,M5,M4)       = 0.0_r8
  ENDIF
  end subroutine BoundaryRunoffandSnowXY
! ----------------------------------------------------------------------

  subroutine BoundaryRunoffandSnowZ(M,N,NN,M1,M2,M3,M4,M5,M6)
  implicit none
  integer, intent(in) :: M,N,NN,M1,M2,M3,M4,M5,M6
  integer :: K,idg,NTN,ids,idom
  real(r8) :: VFLW

  IF(NN.EQ.1.AND.WaterFlow2MicPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &
    .OR.NN.EQ.2.AND.WaterFlow2MicPM_3D(M,N,M6,M5,M4).LT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,M6,M5,M4)/VLWatMicPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_MicpTranspFlxM_3D(idom,K,N,M6,M5,M4)=VFLW*AZMAX1(DOM_MicP2(idom,K,M3,M2,M1))
      enddo
    enddo
    !does not include NH3 and NH3B
    DO idg=idg_beg,idg_end-2
      R3PoreSolFlx_3D(idg,N,M6,M5,M4)=VFLW*AZMAX1(trc_solml2_vr(idg,M3,M2,M1))
    ENDDo

    DO NTN=ids_nuts_beg,ids_nuts_end
      R3PoreSolFlx_3D(NTN,N,M6,M5,M4)=VFLW*AZMAX1(trc_solml2_vr(NTN,M3,M2,M1))*trcs_VLN_vr(NTN,M3,M2,M1)
    ENDDO
!
!     SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
!
  ELSE
    DO  K=1,jcplx
      DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    R3PoreSolFlx_3D(idg_beg:idg_end-2,N,M6,M5,M4)=0.0_r8

   !add irrigation flux
    DO ids=ids_nuts_beg,ids_nuts_end
      R3PoreSolFlx_3D(ids,N,M6,M5,M4)=WaterFlow2MicPM_3D(M,N,M6,M5,M4) &
        *trcn_irrig(ids,M3,M2,M1)*trcs_VLN_vr(ids,M3,M2,M1)
    ENDDO

  ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
!
!     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     RFH*S=solute diffusive flux through macropore
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  IF(NN.EQ.1.AND.WaterFlow2MacPM_3D(M,N,M6,M5,M4).GT.0.0 &
    .OR.NN.EQ.2.AND.WaterFlow2MacPM_3D(M,N,M6,M5,M4).LT.0.0)THEN
    IF(VLWatMacPM(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,M6,M5,M4)/VLWatMacPM(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_MacpTranspFlxM_3D(idom,K,N,M6,M5,M4)=VFLW*AZMAX1(DOM_MacP2(idom,K,M3,M2,M1))
      enddo
    enddo
   !exclude NH3
    DO idg=idg_beg,idg_end-2
      R3PoreSoHFlx_3D(idg,N,M6,M5,M4)=VFLW*AZMAX1(trc_soHml2_vr(idg,M3,M2,M1))
    ENDDO

    DO NTN=ids_nuts_beg,ids_nuts_end
      R3PoreSoHFlx_3D(NTN,N,M6,M5,M4)=VFLW*AZMAX1(trc_soHml2_vr(NTN,M3,M2,M1))*trcs_VLN_vr(NTN,M3,M2,M1)
    ENDDO
!
!     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
!
  ELSE
    DO  K=1,jcplx
      DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    R3PoreSoHFlx_3D(ids_beg:ids_end,N,M6,M5,M4)=0.0_r8
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
    do idom=idom_beg,idom_end
      DOM_MicpTransp_3D(idom,K,N,M6,M5,M4)=DOM_MicpTransp_3D(idom,K,N,M6,M5,M4) &
        +DOM_MicpTranspFlxM_3D(idom,K,N,M6,M5,M4)
      DOM_3DMacp_Transp_flx(idom,K,N,M6,M5,M4)=DOM_3DMacp_Transp_flx(idom,K,N,M6,M5,M4) &
        +DOM_MacpTranspFlxM_3D(idom,K,N,M6,M5,M4)
    enddo
  enddo

  DO ids=ids_beg,ids_end
    trcs_TransptMicP_3D(ids,N,M6,M5,M4)=trcs_TransptMicP_3D(ids,N,M6,M5,M4)+R3PoreSolFlx_3D(ids,N,M6,M5,M4)
    trcs_TransptMacP_3D(ids,N,M6,M5,M4)=trcs_TransptMacP_3D(ids,N,M6,M5,M4)+R3PoreSoHFlx_3D(ids,N,M6,M5,M4)
  ENDDO

  end subroutine BoundaryRunoffandSnowZ

! ----------------------------------------------------------------------
  subroutine BoundaryRunoffandSnow(L,N,NN,M,MX,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
  implicit none

  integer, intent(in) :: L,N, NN, M,MX, N1, N2,N3, M1, M2,M3,M4, M5,M6
  real(r8), intent(in):: RCHQF
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K,idg

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

!     WaterFlow2MicPM_3D=water flux through soil micropore from watsub.f
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     R*FLS=convective solute flux through micropores
!     R*FLW,R*FLB=convective solute flux through micropores in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
    IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(FlowDirIndicator(M2,M1).NE.3 .OR. N.EQ.3)THEN

        call BoundaryRunoffandSnowZ(M,N,NN,M1,M2,M3,M4,M5,M6)
      ENDIF
    ENDIF
!
!     GASOUS LOSS WITH SUBSURFACE MICROPORE WATER GAIN
!
!     WaterFlow2MicPM_3D,WaterFlow2MacPM_3D=micropore,macropore water flux from watsub.f
!     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
!     VLsoiAirPM=air-filled porosity
!     R*FLG=convective gas flux
!     X*FLG=convective gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
    FLGM=(WaterFlow2MicPM_3D(M,N,M6,M5,M4)+WaterFlow2MacPM_3D(M,N,M6,M5,M4))*dt_GasCyc
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
      DO idg=idg_beg,idg_NH3
        RGasADFlx_3D(idg,N,M6,M5,M4)=VFLW*AZMAX1(trc_gasml2_vr(idg,M3,M2,M1))
        Gas_3DAdvDif_Flx_vr(idg,N,M6,M5,M4)=Gas_3DAdvDif_Flx_vr(idg,N,M6,M5,M4)+RGasADFlx_3D(idg,N,M6,M5,M4)
      ENDDO

    ELSE
      RGasADFlx_3D(idg_beg:idg_NH3,N,M6,M5,M4)=0.0_r8
    ENDIF
  ELSE
      DO  K=1,jcplx
        DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
      enddo
      R3PoreSoHFlx_3D(ids_beg:ids_end,N,M6,M5,M4)=0.0_r8
      RGasADFlx_3D(idg_beg:idg_NH3,N,M6,M5,M4)=0.0_r8
  ENDIF
  end subroutine BoundaryRunoffandSnow
!------------------------------------------------------------------------------------------

  subroutine NetOverlandFluxXY(M,N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N4,N5,N4B,N5B

  integer :: NN,K,idg,ids,idom

  DO NN=1,2
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        dom_TFloXSurRunoff(idom,K,N2,N1)=dom_TFloXSurRunoff(idom,K,N2,N1) &
          +dom_2DFloXSurRunoffM(idom,K,N,NN,N2,N1)
      enddo
    enddo

    DO idg=idg_beg,idg_NH3
      trcg_TFloXSurRunoff(idg,N2,N1)=trcg_TFloXSurRunoff(idg,N2,N1)+trcg_2DFloXSurRunoffM(idg,N,NN,N2,N1)
    ENDDO

    DO ids=ids_nut_beg,ids_nuts_end
      trcn_TFloXSurRunoff_2D(ids,N2,N1)=trcn_TFloXSurRunoff_2D(ids,N2,N1)+trcn_2DFloXSurRunoffM(ids,N,NN,N2,N1)
    ENDDO

    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          dom_TFloXSurRunoff(idom,K,N2,N1)=dom_TFloXSurRunoff(idom,K,N2,N1) &
            -dom_2DFloXSurRunoffM(idom,K,N,NN,N5,N4)
        enddo
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_TFloXSurRunoff(idg,N2,N1)=trcg_TFloXSurRunoff(idg,N2,N1)-trcg_2DFloXSurRunoffM(idg,N,NN,N5,N4)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_TFloXSurRunoff_2D(ids,N2,N1)=trcn_TFloXSurRunoff_2D(ids,N2,N1)-trcn_2DFloXSurRunoffM(ids,N,NN,N5,N4)
      ENDDO
    ENDIF
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          dom_TFloXSurRunoff(idom,K,N2,N1)=dom_TFloXSurRunoff(idom,K,N2,N1)-dom_2DFloXSurRunoffM(idom,K,N,NN,N5B,N4B)
        enddo
      enddo
      DO idg=idg_beg,idg_NH3
        trcg_TFloXSurRunoff(idg,N2,N1)=trcg_TFloXSurRunoff(idg,N2,N1)-trcg_2DFloXSurRunoffM(idg,N,NN,N5B,N4B)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_TFloXSurRunoff_2D(ids,N2,N1)=trcn_TFloXSurRunoff_2D(ids,N2,N1)-trcn_2DFloXSurRunoffM(ids,N,NN,N5B,N4B)
      ENDDO
    ENDIF
  enddo
!
!     NET OVERLAND SOLUTE FLUX IN SNOW
!
!     TQS*=net solute flux in snow transfer
!     RQS*=solute flux in snow transfer
!
  DO idg=idg_beg,idg_NH3
    trcg_SnowDrift(idg,N2,N1)=trcg_SnowDrift(idg,N2,N1)+trcg_2DSnowDrift(idg,N,N2,N1)-trcg_2DSnowDrift(idg,N,N5,N4)
  ENDDO

  do ids=ids_nut_beg,ids_nuts_end
    trcn_SnowDrift(ids,N2,N1)=trcn_SnowDrift(ids,N2,N1)+trcn_2DSnowDrift(ids,N,N2,N1)-trcn_2DSnowDrift(ids,N,N5,N4)
  enddo
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
!     VLSnowHeatCapM,VLHeatCapSnowMin_col=current,minimum volumetric heat capacity of snowpack
!     T*BLS=net solute flux in snowpack
!     R*BLS=solute flux in snowpack
!
      ELSEIF(N.EQ.3)THEN
! vertical direction
        call NetOverlandVerticalFlux(M,N1,N2,NY,NX)
      ENDIF
    ENDIF
  ENDIF
  end subroutine NetOverlandFlux
!------------------------------------------------------------------------------------------

  subroutine NetOverlandVerticalFlux(M,N1,N2,NY,NX)
  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX

  integer :: LS,LS2
  integer :: idg,NTN

  DO  LS=1,JS
    IF(VLSnowHeatCapM_snvr(M,LS,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS .AND. VLSnowHeatCapM_snvr(M,LS2,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
        DO idg=idg_beg,idg_NH3
          trcg_TBLS_snvr(idg,LS,N2,N1)=trcg_TBLS_snvr(idg,LS,N2,N1) &
            +trcg_advW_snvr(idg,LS,N2,N1)-trcg_advW_snvr(idg,LS2,N2,N1)
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
        DO idg=idg_beg,idg_NH3
          trcg_TBLS_snvr(idg,LS,N2,N1)=trcg_TBLS_snvr(idg,LS,N2,N1) &
            +trcg_advW_snvr(idg,LS,N2,N1)-R3PoreSolFlx_3D(idg,3,0,N2,N1) &
            -R3PoreSolFlx_3D(idg,3,NUM(N2,N1),N2,N1)
  !    3-R3PoreSoHFlx_3D(idg,3,NUM(N2,N1),N2,N1)
        ENDDO

        trcg_TBLS_snvr(idg_NH3,LS,N2,N1)=trcg_TBLS_snvr(idg_NH3,LS,N2,N1) &
          -R3PoreSolFlx_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)
!    3-R3PoreSoHFlx_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)

! check trnsfr.f, loop 1205
        DO NTN=0,ids_nuts
          trcn_TBLS(ids_NH4+NTN,LS,N2,N1)=trcn_TBLS(ids_NH4+NTN,LS,N2,N1) &
            +trcn_RBLS(ids_NH4+NTN,LS,N2,N1)-R3PoreSolFlx_3D(ids_NH4+NTN,3,0,N2,N1) &
            -R3PoreSolFlx_3D(ids_NH4+NTN,3,NUM(N2,N1),N2,N1) &
            -R3PoreSolFlx_3D(ids_NH4B+NTN,3,NUM(N2,N1),N2,N1)
!    -R3PoreSoHFlx_3D(ids_NH4+NTN,3,NUM(N2,N1),N2,N1) &
!    -R3PoreSoHFlx_3D(ids_NH4B+NTN,3,NUM(N2,N1),N2,N1)
        ENDDO
      ENDIF
    ENDIF
  enddo
  end subroutine NetOverlandVerticalFlux
!------------------------------------------------------------------------------------------

  subroutine NetFluxMicroandMacropores(NY,NX,N,M,MX,N1,N2,N3,N4,N5,N6)
  implicit none

  integer, intent(in) :: NY,NX,N,M,N1,N2,N3,N4,N5,MX
  integer, intent(inout) :: N6
  integer :: K,LL,ids,idg,idom

  DO LL=N6,NL(NY,NX)
    IF(VLSoilPoreMicP_vr(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
      N6=LL
      exit
    ENDIF
  ENDDO

  IF(M.NE.MX)THEN
    IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Transp2Micp_vr(idom,K,N3,N2,N1)=DOM_Transp2Micp_vr(idom,K,N3,N2,N1) &
            +DOM_MicpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)
          DOM_Transp2Macp_flx(idom,K,N3,N2,N1)=DOM_Transp2Macp_flx(idom,K,N3,N2,N1) &
            +DOM_MacpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)
        enddo
      enddo
      DO ids=ids_beg,ids_end
        TR3MicPoreSolFlx_vr(ids,N3,N2,N1)=TR3MicPoreSolFlx_vr(ids,N3,N2,N1) &
          +R3PoreSolFlx_3D(ids,N,N3,N2,N1)-R3PoreSolFlx_3D(ids,N,N6,N5,N4)
        TR3MacPoreSolFlx_vr(ids,N3,N2,N1)=TR3MacPoreSolFlx_vr(ids,N3,N2,N1) &
          +R3PoreSoHFlx_3D(ids,N,N3,N2,N1)-R3PoreSoHFlx_3D(ids,N,N6,N5,N4)
      ENDDO

    ELSE
      DO  K=1,jcplx
        DOM_Transp2Micp_vr(idom_beg:idom_end,K,N3,N2,N1)=0.0_r8
        DOM_Transp2Macp_flx(idom_beg:idom_end,K,N3,N2,N1)=0.0_r8
      enddo
      TR3MicPoreSolFlx_vr(ids_beg:ids_end,N3,N2,N1)=0.0_r8
      TR3MacPoreSolFlx_vr(ids_beg:ids_end,N3,N2,N1)=0._r8

    ENDIF
  ENDIF
!
!     NET GAS FLUX
!
!     T*FLG=net convective+diffusive gas flux
!     R*FLG=convective+diffusive gas flux
!     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!
  IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    DO idg=idg_beg,idg_NH3
      Gas_AdvDif_Flx_vr(idg,N3,N2,N1)=Gas_AdvDif_Flx_vr(idg,N3,N2,N1) &
        +RGasADFlx_3D(idg,N,N3,N2,N1)-RGasADFlx_3D(idg,N,N6,N5,N4)
    ENDDO
  ELSE
    Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,N3,N2,N1)=0._r8
  ENDIF
  end subroutine NetFluxMicroandMacropores

end module BoundaryTranspMod
