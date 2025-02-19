module BoundaryTranspMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: isclose, AZMAX1
  use TracerPropMod, only: MolecularWeight
  use ClimForcDataType
  use DebugToolMod
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

  public :: XBoundaryFluxMM
  contains

!------------------------------------------------------------------------------------------

  subroutine XBoundaryFluxMM(M,MX,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) :: M,MX,NHW, NHE, NVN, NVS

  character(len=*), parameter :: subname='XBoundaryFluxMM'
  integer :: NY,NX,L
  integer :: N1,N2,N3,N4,N5,N6,NN,N
  integer :: N4B,N5B,M1,M2,M3,M4,M5,M6
  real(r8) :: RCHQF,RCHGFU,RCHGFT,XN


  call PrintInfo('beg '//subname)
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
            IF(N.EQ.iEastWestDirection)THEN
              !WEST-EAST
              N4  = NX+1  !right boundary
              N5  = NY
              N4B = NX-1  !left boundary
              N5B = NY
              N6  = L
              IF(NN.EQ.iOutflow)THEN   !eastward
                !eastern boundary
                IF(NX.EQ.NHE)THEN
                  M1 = NX
                  M2 = NY
                  M3 = L
                  !target grid
                  M4 = NX+1
                  M5 = NY
                  M6 = L
                  XN = -1.0    !going out
                  RCHQF  = RechargEastSurf(M2,M1)
                  RCHGFU = RechargEastSubSurf(M2,M1)
                  RCHGFT = RechargRateEastWTBL(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN  !west
                !western boundary
                IF(NX.EQ.NHW)THEN
                  M1 = NX
                  M2 = NY
                  M3 = L
                  !
                  M4 = NX
                  M5 = NY
                  M6 = L
                  XN = 1.0     !coming in
                  RCHQF  = RechargWestSurf(M5,M4)
                  RCHGFU = RechargWestSubSurf(M5,M4)
                  RCHGFT = RechargRateWestWTBL(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              !NORTH-SOUTH
              N4  = NX
              N5  = NY+1  !south boundary
              N4B = NX
              N5B = NY-1  !north boundary
              N6  = L
              IF(NN.EQ.iOutflow)THEN  !south
                ! southern boundary
                IF(NY.EQ.NVS)THEN
                  M1    = NX
                  M2    = NY
                  M3    = L
                  !target grid
                  M4    = NX
                  M5    = NY+1
                  M6    = L
                  XN    = -1.0_r8    !going out
                  RCHQF  = RechargSouthSurf(M2,M1)
                  RCHGFU = RechargSouthSubSurf(M2,M1)
                  RCHGFT = RechargRateSouthWTBL(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN  !north
                ! northern boundary
                IF(NY.EQ.NVN)THEN
                  M1 = NX
                  M2 = NY
                  M3 = L
                  !target
                  M4 = NX
                  M5 = NY
                  M6 = L
                  XN = 1.0_r8   !coming in
                  RCHQF  = RechargNorthSurf(M5,M4)
                  RCHGFU = RechargNorthSubSurf(M5,M4)
                  RCHGFT = RechargRateNorthWTBL(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN !vertical
              !source
              N1 = NX
              N2 = NY
              N3 = L
              !target
              N4 = NX
              N5 = NY
              N6 = L+1
              IF(NN.EQ.iOutflow)THEN
                !lower boundary
                IF(L.EQ.NL(NY,NX))THEN
                  M1=NX
                  M2=NY
                  M3=L
                  !target grid
                  M4=NX
                  M5=NY
                  M6=L+1
                  XN=-1.0_r8  !going out
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN
                !nothing for the upper boundary
                cycle
              ENDIF
            ENDIF
            
!     SURFACE SOLUTE TRANSPORT FROM BOUNDARY SURFACE
!     RUNOFF IN 'WATSUB' AND CONCENTRATIONS IN THE SURFACE SOIL LAYER
!
            call OverLandTracerFlowMM(L,N,NN,M,MX,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
!
          ENDDO D9575
          !
          !     NET OVERLAND SOLUTE FLUX IN WATER
          !
            IF(M.NE.MX)call NetTracerFlowOverLandM(L,N,M,MX,NY,NX,N1,N2,N4B,N5B,N4,N5)
!
!     TOTAL SOLUTE FLUX IN MICROPORES AND MACROPORES
!
          IF(FlowDirIndicator(N2,N1).NE.3 .OR. N.EQ.iVerticalDirection)THEN
            call NetTracerFlowXSoilPoresMM(NY,NX,N,M,MX,N1,N2,N3,N4,N5,N6)
          ENDIF
        ENDDO D9580
      ENDDO D9585
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine XBoundaryFluxMM
!------------------------------------------------------------------------------------------

  subroutine SurfTracerRunoffXYM(M,N,NN,N1,N2,M4,M5,RCHQF)

  implicit none
  integer, intent(in) :: M,N,NN,N1,N2,M4,M5
  real(r8), intent(in) :: RCHQF
  real(r8) :: FQRM
  integer :: K,idg,idn,ids,idom

  IF(.not.XGridRunoffFlag(NN,N,N2,N1) .OR. isclose(RCHQF,0.0_r8) .OR. SurfRunoffWatFluxM_2DH(M,N2,N1).LE.ZEROS(N2,N1))THEN
    DO  K=1,jcplx
      DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
    ENDDO
    trcg_FloXSurRof_flxM_2DH(idg_beg:idg_NH3,N,NN,M5,M4)          = 0.0_r8
    trcn_FloXSurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,NN,M5,M4) = 0.0_r8
  ELSE
!
!     SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
    IF((NN.EQ.iOutflow .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.iInflow .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      FQRM=QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/SurfRunoffWatFluxM_2DH(M,N2,N1)
      DO  K=1,jcplx
        DO idom=idom_beg,idom_end
          DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,M5,M4)=DOM_FloXSurRunoff_flxM(idom,K,N2,N1)*FQRM
        ENDDO
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_FloXSurRof_flxM_2DH(idg,N,NN,M5,M4)=trcg_FloXSurRunoff_flxM(idg,N2,N1)*FQRM
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_FloXSurRof_flxM_2DH(ids,N,NN,M5,M4)=trcn_FloXSurRunoff_flxM(ids,N2,N1)*FQRM
      ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_FloXSurRunoff_2D(idom,K,N,NN,M5,M4)=DOM_FloXSurRunoff_2D(idom,K,N,NN,M5,M4)+DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,M5,M4)
        enddo
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)=trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)+trcg_FloXSurRof_flxM_2DH(idg,N,NN,M5,M4)
      ENDDO

      DO idn=ids_nut_beg,ids_nuts_end
        trcn_FloXSurRunoff_2D(idn,N,NN,M5,M4)=trcn_FloXSurRunoff_2D(idn,N,NN,M5,M4)+trcn_FloXSurRof_flxM_2DH(idn,N,NN,M5,M4)
      ENDDO
!
!     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
    ELSEIF((NN.EQ.iInflow .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.iOutflow .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      DO  K=1,jcplx
        DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_FloXSurRof_flxM_2DH(idg,N,NN,M5,M4) = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*trcg_rain_mole_conc_col(idg,M5,M4)*MolecularWeight(idg)
      ENDDO

      trcn_FloXSurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,NN,M5,M4)=0.0_r8

      DO idg=idg_beg,idg_NH3
        trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)=trcg_FloXSurRunoff_2D(idg,N,NN,M5,M4)+trcg_FloXSurRof_flxM_2DH(idg,N,NN,M5,M4)
      ENDDO
    ELSE
      DO  K=1,jcplx
        DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
      enddo
      trcg_FloXSurRof_flxM_2DH(idg_beg:idg_NH3,N,NN,M5,M4)          = 0.0_r8
      trcn_FloXSurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,NN,M5,M4) = 0.0_r8
    ENDIF

  ENDIF
!
!     BOUNDARY SNOW FLUX
!
  IF(NN.EQ.iOutflow)THEN
    trcg_SnowDrift_flxM_2D(idg_beg:idg_NH3,N,M5,M4) = 0.0_r8
    trcn_SnowDrift_flxM_2D(ids_NH4,N,M5,M4)         = 0.0_r8
    trcn_SnowDrift_flxM_2D(ids_NO3,N,M5,M4)         = 0.0_r8
    trcn_SnowDrift_flxM_2D(ids_H1PO4,N,M5,M4)       = 0.0_r8
    trcn_SnowDrift_flxM_2D(ids_H2PO4,N,M5,M4)       = 0.0_r8
  ENDIF

  end subroutine SurfTracerRunoffXYM
! ----------------------------------------------------------------------

  subroutine BoundaryRunoffandSnowZM(M,N,NN,M1,M2,M3,M4,M5,M6)
  implicit none
  integer, intent(in) :: M,N,NN
  integer, intent(in) :: M1,M2,M3
  integer, intent(in) :: M4,M5,M6
  integer :: K,idg,idn,ids,idom
  real(r8) :: VFLW

  IF(NN.EQ.iOutflow .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &
    .OR. NN.EQ.iInflow .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).LT.0.0_r8)THEN

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

    DO idn=ids_beg,ids_end
      trcs_MicpTranspFlxM_3D(idn,N,M6,M5,M4)=VFLW*AZMAX1(trcs_solml2_vr(idn,M3,M2,M1))*trcs_VLN_vr(idn,M3,M2,M1)
    ENDDO
!
!     SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
!
  ELSE
    DO  K=1,jcplx
      DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    trcs_MicpTranspFlxM_3D(idg_beg:idg_end-2,N,M6,M5,M4)=0.0_r8

   !add irrigation flux
    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_MicpTranspFlxM_3D(ids,N,M6,M5,M4)=WaterFlow2MicPM_3D(M,N,M6,M5,M4) &
        *trcn_irrig_vr(ids,M3,M2,M1)*trcs_VLN_vr(ids,M3,M2,M1)
    ENDDO

  ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
!
  IF(NN.EQ.iOutflow .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &
    .OR. NN.EQ.iInflow .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).LT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,M6,M5,M4)/VLWatMacPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_MacpTranspFlxM_3D(idom,K,N,M6,M5,M4)=VFLW*AZMAX1(DOM_MacP2(idom,K,M3,M2,M1))
      enddo
    enddo

    DO idn=ids_beg,ids_end
      trcs_MacpTranspFlxM_3D(idn,N,M6,M5,M4)=VFLW*AZMAX1(trcs_soHml2_vr(idn,M3,M2,M1))*trcs_VLN_vr(idn,M3,M2,M1)
    ENDDO
!
!     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
!
  ELSE
    DO  K=1,jcplx
      DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    trcs_MacpTranspFlxM_3D(ids_beg:ids_end,N,M6,M5,M4)=0.0_r8
  ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTransp_3D(idom,K,N,M6,M5,M4)=DOM_MicpTransp_3D(idom,K,N,M6,M5,M4) &
        +DOM_MicpTranspFlxM_3D(idom,K,N,M6,M5,M4)
      DOM_Macp_Transp_flx_3D(idom,K,N,M6,M5,M4)=DOM_Macp_Transp_flx_3D(idom,K,N,M6,M5,M4) &
        +DOM_MacpTranspFlxM_3D(idom,K,N,M6,M5,M4)
    enddo
  enddo

  DO ids=ids_beg,ids_end
    trcs_TransptMicP_3D(ids,N,M6,M5,M4)=trcs_TransptMicP_3D(ids,N,M6,M5,M4)+trcs_MicpTranspFlxM_3D(ids,N,M6,M5,M4)
    trcs_TransptMacP_3D(ids,N,M6,M5,M4)=trcs_TransptMacP_3D(ids,N,M6,M5,M4)+trcs_MacpTranspFlxM_3D(ids,N,M6,M5,M4)
  ENDDO

  end subroutine BoundaryRunoffandSnowZM

! ----------------------------------------------------------------------
  subroutine OverLandTracerFlowMM(L,N,NN,M,MX,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
  implicit none

  integer, intent(in) :: L,N, NN, M,MX 
  integer, intent(in) :: N1, N2,N3
  integer, intent(in) :: M1, M2,M3,M4, M5,M6
  real(r8), intent(in):: RCHQF
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K,idg

! begin_execution
!
  IF(M.NE.MX)THEN
    IF(L.EQ.NUM(M2,M1) .AND. N.NE.iVerticalDirection)THEN
      call SurfTracerRunoffXYM(M,N,NN,N1,N2,M4,M5,RCHQF)
    ENDIF
!
    IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(FlowDirIndicator(M2,M1).NE.3 .OR. N.EQ.iVerticalDirection)THEN        
        call BoundaryRunoffandSnowZM(M,N,NN,M1,M2,M3,M4,M5,M6)
      ENDIF
    ENDIF
!
  ELSE
    DO  K=1,jcplx
      DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    trcs_MacpTranspFlxM_3D(ids_beg:ids_end,N,M6,M5,M4) = 0.0_r8
  ENDIF

!     GASOUS LOSS WITH SUBSURFACE MICROPORE WATER GAIN
!   dt_GasCyc=1/NPT, NPT is number of gas iterations per M
  FLGM=(WaterFlow2MicPM_3D(M,N,M6,M5,M4)+WaterFlow2MacPM_3D(M,N,M6,M5,M4))*dt_GasCyc

  IF(NN.EQ.iOutflow .AND. FLGM.LT.0.0_r8 .OR. NN.EQ.iInflow .AND. FLGM.GT.0.0_r8)THEN
    IF(VLsoiAirPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=-AMAX1(-VFLWX,AMIN1(VFLWX,FLGM/VLsoiAirPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
!
!     HOURLY GAS FLUX FOR USE IN REDIST.F
!
!     X*FLG=hourly convective gas flux
!
    DO idg=idg_beg,idg_NH3
      RGasADFlxMM_3D(idg,N,M6,M5,M4)        = VFLW*AZMAX1(trc_gasml2_vr(idg,M3,M2,M1))
      Gas_AdvDif_Flx_3D(idg,N,M6,M5,M4)   = Gas_AdvDif_Flx_3D(idg,N,M6,M5,M4)+RGasADFlxMM_3D(idg,N,M6,M5,M4)
    ENDDO
  ELSE
    RGasADFlxMM_3D(idg_beg:idg_NH3,N,M6,M5,M4)=0.0_r8
  ENDIF

  end subroutine OverLandTracerFlowMM
!------------------------------------------------------------------------------------------

  subroutine NetOverlandFluxXYM(M,N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N4,N5,N4B,N5B

  character(len=*), parameter :: subname='NetOverlandFluxXYM'
  integer :: NN,K,idg,ids,idom

  call PrintInfo('beg '//subname)
  DO NN=1,2
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1) &
          +DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N2,N1)
      enddo
    enddo

    DO idg=idg_beg,idg_NH3
      trcg_SurfRunoff_flx(idg,N2,N1)=trcg_SurfRunoff_flx(idg,N2,N1)+trcg_FloXSurRof_flxM_2DH(idg,N,NN,N2,N1)
    ENDDO

    DO ids=ids_nut_beg,ids_nuts_end
      trcn_SurfRunoff_flx(ids,N2,N1)=trcn_SurfRunoff_flx(ids,N2,N1)+trcn_FloXSurRof_flxM_2DH(ids,N,NN,N2,N1)
    ENDDO

    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1) &
            -DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N5,N4)
        enddo
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flx(idg,N2,N1)=trcg_SurfRunoff_flx(idg,N2,N1)-trcg_FloXSurRof_flxM_2DH(idg,N,NN,N5,N4)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flx(ids,N2,N1)=trcn_SurfRunoff_flx(ids,N2,N1)-trcn_FloXSurRof_flxM_2DH(ids,N,NN,N5,N4)
      ENDDO
    ENDIF

    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.iOutflow)THEN
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1)-DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N5B,N4B)
        enddo
      enddo
      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flx(idg,N2,N1)=trcg_SurfRunoff_flx(idg,N2,N1)-trcg_FloXSurRof_flxM_2DH(idg,N,NN,N5B,N4B)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flx(ids,N2,N1)=trcn_SurfRunoff_flx(ids,N2,N1)-trcn_FloXSurRof_flxM_2DH(ids,N,NN,N5B,N4B)
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
    trcg_SnowDrift_flxM(idg,N2,N1)=trcg_SnowDrift_flxM(idg,N2,N1)+trcg_SnowDrift_flxM_2D(idg,N,N2,N1)-trcg_SnowDrift_flxM_2D(idg,N,N5,N4)
  ENDDO

  do ids=ids_nut_beg,ids_nuts_end
    trcn_SnowDrift_flxM(ids,N2,N1)=trcn_SnowDrift_flxM(ids,N2,N1)+trcn_SnowDrift_flxM_2D(ids,N,N2,N1)-trcn_SnowDrift_flxM_2D(ids,N,N5,N4)
  enddo
  call PrintInfo('end '//subname)
  end subroutine NetOverlandFluxXYM

!------------------------------------------------------------------------------------------
!
  subroutine NetTracerFlowOverLandM(L,N,M,MX,NY,NX,N1,N2,N4B,N5B,N4,N5)
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

  IF(L.EQ.NUM(N2,N1))THEN
    IF(N.NE.iVerticalDirection)THEN
    !horizontal flow
      call NetOverlandFluxXYM(M,N,N1,N2,N4,N5,N4B,N5B)
    
    !     NET SOLUTE FLUX IN SNOWPACK
    ELSEIF(N.EQ.iVerticalDirection)THEN
      ! vertical direction
      call NetOverlandFluxZM(M,N1,N2,NY,NX)
    ENDIF
  ENDIF

  end subroutine NetTracerFlowOverLandM
!------------------------------------------------------------------------------------------

  subroutine NetOverlandFluxZM(M,N1,N2,NY,NX)
  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX

  integer :: LS,LS2
  integer :: idg,idn

  DO  LS=1,JS
    IF(VLSnowHeatCapM_snvr(M,LS,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS .AND. VLSnowHeatCapM_snvr(M,LS2,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
        DO idg=idg_beg,idg_NH3
          trcg_Aqua_flxM_snvr(idg,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg,LS,N2,N1) &
            +trcg_AquaAdv_flxM_snvr(idg,LS,N2,N1)-trcg_AquaAdv_flxM_snvr(idg,LS2,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_Aqua_flxM_snvr(idn,LS,N2,N1)=trcn_Aqua_flxM_snvr(idn,LS,N2,N1) &
            +trcn_AquaAdv_flxM_snvr(idn,LS,N2,N1)-trcn_AquaAdv_flxM_snvr(idn,LS2,N2,N1)
        ENDDO
      ELSE
!
!     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
!
        ! exclude NH3 and NH3B
        DO idg=idg_beg,idg_NH3
          trcg_Aqua_flxM_snvr(idg,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg,LS,N2,N1) &
            +trcg_AquaAdv_flxM_snvr(idg,LS,N2,N1)-trcs_MicpTranspFlxM_3D(idg,3,0,N2,N1) &
            -trcs_MicpTranspFlxM_3D(idg,3,NUM(N2,N1),N2,N1)
  !    3-trcs_MacpTranspFlxM_3D(idg,3,NUM(N2,N1),N2,N1)
        ENDDO

        trcg_Aqua_flxM_snvr(idg_NH3,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg_NH3,LS,N2,N1) &
          -trcs_MicpTranspFlxM_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)
!    3-trcs_MacpTranspFlxM_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)

! check trnsfr.f, loop 1205
        DO idn=0,ids_nuts
          trcn_Aqua_flxM_snvr(ids_NH4+idn,LS,N2,N1)=trcn_Aqua_flxM_snvr(ids_NH4+idn,LS,N2,N1) &
            +trcn_AquaAdv_flxM_snvr(ids_NH4+idn,LS,N2,N1)-trcs_MicpTranspFlxM_3D(ids_NH4+idn,3,0,N2,N1) &
            -trcs_MicpTranspFlxM_3D(ids_NH4+idn,3,NUM(N2,N1),N2,N1) &
            -trcs_MicpTranspFlxM_3D(ids_NH4B+idn,3,NUM(N2,N1),N2,N1)
!    -trcs_MacpTranspFlxM_3D(ids_NH4+idn,3,NUM(N2,N1),N2,N1) &
!    -trcs_MacpTranspFlxM_3D(ids_NH4B+idn,3,NUM(N2,N1),N2,N1)
        ENDDO
      ENDIF
    ENDIF
  enddo
  end subroutine NetOverlandFluxZM
!------------------------------------------------------------------------------------------

  subroutine NetTracerFlowXSoilPoresMM(NY,NX,N,M,MX,N1,N2,N3,N4,N5,N6)
  !
  !Description
  !Vertical or lateral flux between grid grids
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
          DOM_Transp2Micp_flxM_vr(idom,K,N3,N2,N1)=DOM_Transp2Micp_flxM_vr(idom,K,N3,N2,N1) &
            +DOM_MicpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)
          DOM_Transp2Macp_flxM_vr(idom,K,N3,N2,N1)=DOM_Transp2Macp_flxM_vr(idom,K,N3,N2,N1) &
            +DOM_MacpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)
        enddo
      enddo
      DO ids=ids_beg,ids_end
        trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1)=trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1) &
          +trcs_MicpTranspFlxM_3D(ids,N,N3,N2,N1)-trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)
        trcs_Transp2Macp_flxM_vr(ids,N3,N2,N1)=trcs_Transp2Macp_flxM_vr(ids,N3,N2,N1) &
          +trcs_MacpTranspFlxM_3D(ids,N,N3,N2,N1)-trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)
      ENDDO

    ELSE
      DO  K=1,jcplx
        DOM_Transp2Micp_flxM_vr(idom_beg:idom_end,K,N3,N2,N1)  = 0.0_r8
        DOM_Transp2Macp_flxM_vr(idom_beg:idom_end,K,N3,N2,N1) = 0.0_r8
      enddo
      trcs_Transp2Micp_flxM_vr(ids_beg:ids_end,N3,N2,N1) = 0.0_r8
      trcs_Transp2Macp_flxM_vr(ids_beg:ids_end,N3,N2,N1) = 0._r8
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
      Gas_AdvDif_FlxMM_vr(idg,N3,N2,N1)=Gas_AdvDif_FlxMM_vr(idg,N3,N2,N1)+RGasADFlxMM_3D(idg,N,N3,N2,N1)-RGasADFlxMM_3D(idg,N,N6,N5,N4)
    ENDDO
  ELSE
    Gas_AdvDif_FlxMM_vr(idg_beg:idg_NH3,N3,N2,N1)=0._r8
  ENDIF
  end subroutine NetTracerFlowXSoilPoresMM

end module BoundaryTranspMod
