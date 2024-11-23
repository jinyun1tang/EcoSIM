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
  use EcoSIMConfig, only : jcplx => jcplxc,NumMicbFunGrupsPerCmplx=>NumMicbFunGrupsPerCmplx
  use EcoSIMConfig, only : nlbiomcp=>NumLiveMicrbCompts
  use ErosionBalMod
  use SnowBalanceMod
  use SnowTransportMod
  use EcoSIMCtrlMod
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: XGridTranspt
  contains

  subroutine XGridTranspt(I,J,NY,NX,LG)
  !
  ! Energy and material transport across grids

  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(out) :: LG

  integer :: L,N,LX
  integer :: N1,N2,N3,N4,N5,N6
  integer :: N4B,N5B
  real(r8) :: VTATM,VTGAS
  real(r8) :: trcg_VOLG(idg_beg:idg_end)
  
! begin_execution
  if(lverb)write(*,*)'XGridTranspt'
  call ZeroFluxArrays(NY,NX)

  call ZeroFluxAccumulators(NY,NX)

  LG=0;LX=0

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
    trcg_VOLG(idg_CO2) = trc_gasml_vr(idg_CO2,L,NY,NX)/catomw
    trcg_VOLG(idg_CH4) = trc_gasml_vr(idg_CH4,L,NY,NX)/catomw
    trcg_VOLG(idg_O2)  = trc_gasml_vr(idg_O2,L,NY,NX)/32.0_r8
    trcg_VOLG(idg_N2)  = trc_gasml_vr(idg_N2,L,NY,NX)/28.0_r8
    trcg_VOLG(idg_N2O) = trc_gasml_vr(idg_N2O,L,NY,NX)/28.0_r8
    trcg_VOLG(idg_NH3) = trc_gasml_vr(idg_NH3,L,NY,NX)/natomw
    trcg_VOLG(idg_H2)  = trc_gasml_vr(idg_H2,L,NY,NX)/2.0_r8

    VTATM=AZMAX1(1.2194E+04_r8*VLsoiAirP_vr(L,NY,NX)/TKS_vr(L,NY,NX))
!   NH3B does not have explicit gas species, so there is an inconsistency
!   with respect to the actual ebullition calculation, which involves
!   NH3B

    VTGAS=sum(trcg_VOLG(idg_beg:idg_NH3))

    !air-concentration insignificant, or total gas volume > allwed air
    !LX==1, then too less air, or gas pressure > atmosphere
    !THETX minimum air-filled porosity
    IF(ThetaAir_vr(L,NY,NX).LT.THETX .OR. VTGAS.GT.VTATM)LX=1

    IF(ThetaAir_vr(L,NY,NX).GE.THETX .AND. LX.EQ.0)LG=L
    !make a copy of soil water/ice in micro- and macropores
    VLWatMicP1_vr(L,NY,NX) = VLWatMicP_vr(L,NY,NX)
    VLiceMicP1_vr(L,NY,NX) = VLiceMicP_vr(L,NY,NX)
    VLWatMacP1_vr(L,NY,NX) = VLWatMacP_vr(L,NY,NX)
    VLiceMacP1_vr(L,NY,NX) = VLiceMacP_vr(L,NY,NX)

!
  !     NET WATER, HEAT, GAS, SOLUTE, SEDIMENT FLUX
  !
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  ! source
    N1=NX;N2=NY;N3=L
    D8580: DO N=FlowDirIndicator(NY,NX),3
      IF(N.EQ.iEastWestDirection)THEN
        !exchange in the x direction, west-east
        N4=NX+1   !east
        N5=NY
        N4B=NX-1  !west
        N5B=NY
        N6=L
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        !exchange in the y direction, north-south
        N4=NX
        N5=NY+1    !south
        N4B=NX
        N5B=NY-1   !north
        N6=L
      ELSEIF(N.EQ.iVerticalDirection)THEN
        !vertical
        N4=NX
        N5=NY
        N6=L+1
      ENDIF
!
!
      IF(L.EQ.NUM(N2,N1))THEN
        !top layer runoff
        IF(N.NE.iVerticalDirection)THEN
          !horizontal exchange
          call OMH2OFluxesFromRunoff(N,N1,N2,N4,N5,N4B,N5B)

          call MassFluxThruSnowRunoff(N,N1,N2,N4,N5,N4B,N5B)
          !
        ELSEIF(N.EQ.iVerticalDirection)THEN
          !vertical direction
          call VerticalSaltFluxThruSnowpack(I,J,N1,N2,NY,NX)
        ENDIF
!
        ! TOTAL FLUXES FROM SEDIMENT TRANSPORT
        call TotalFluxFromSedmentTransp(N,N1,N2,N4,N5,N4B,N5B,NY,NX)
      ENDIF
!
      call FluxThruGrids(I,J,N,N1,N2,N3,N4,N5,N6,NY,NX)
    ENDDO D8580
!
!     NET FREEZE-THAW
!
!     WatIceThawMicP,WatIceThawMacP_vr=net freeze-thaw flux in micropores,macropores
!     THeatSoiThaw=net freeze-thaw latent heat flux
!     THAW,TLIceThawMacP=freeze-thaw flux in micropores,macropores from watsub.f
!     TLPhaseChangeHeat2Soi=freeze-thaw latent heat flux from watsub.f

    WatIceThawMicP_vr(N3,N2,N1) = WatIceThawMicP_vr(N3,N2,N1)+TLIceThawMicP(N3,N2,N1)
    WatIceThawMacP_vr(N3,N2,N1) = WatIceThawMacP_vr(N3,N2,N1)+TLIceThawMacP(N3,N2,N1)
    THeatSoiThaw_vr(N3,N2,N1)   = THeatSoiThaw_vr(N3,N2,N1)+TLPhaseChangeHeat2Soi(N3,N2,N1)
  ENDDO D8575
  end subroutine XGridTranspt

!------------------------------------------------------------------------------------------

  subroutine ZeroRunoffArray(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,K
!
!     INITIALIZE NET WATER AND HEAT FLUXES FOR RUNOFF
!
!
!     INITIALIZE NET SOLUTE AND GAS FLUXES FOR RUNOFF
!
  D9960: DO K=1,micpar%NumOfLitrCmplxs
    TOMQRS(idom_beg:idom_end,K,NY,NX)=0.0_r8
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
    TWatFlowCellMicP_vr(L,NY,NX)  = 0.0_r8
    TWatFlowCellMicPX_vr(L,NY,NX) = 0.0_r8
    TWatFlowCellMacP_vr(L,NY,NX)  = 0.0_r8
    THeatFlowCellSoil_vr(L,NY,NX) = 0.0_r8
    WatIceThawMicP_vr(L,NY,NX)    = 0.0_r8
    WatIceThawMacP_vr(L,NY,NX)    = 0.0_r8
    THeatSoiThaw_vr(L,NY,NX)      = 0.0_r8
!
!     INITIALIZE GAS AND SOLUTE NET FLUX ACCUMULATORS WITHIN SOIL
!
    D8595: DO K=1,jcplx
      DOM_Transp2Micp_vr(idom_beg:idom_end,K,L,NY,NX)  = 0.0_r8
      DOM_Transp2Macp_flx(idom_beg:idom_end,K,L,NY,NX) = 0.0_r8
    ENDDO D8595

    trcs_TransptMicP_vr(ids_beg:ids_end,L,NY,NX) = 0.0_r8
    trcs_TransptMacP_vr(ids_beg:ids_end,L,NY,NX) = 0.0_r8
    Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,L,NY,NX)   = 0.0_r8

    IF(salt_model)THEN
      trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,L,NY,NX) = 0.0_r8
      trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,L,NY,NX) = 0.0_r8
    ENDIF
  ENDDO
  end subroutine ZeroFluxAccumulators

!------------------------------------------------------------------------------------------

  subroutine OMH2OFluxesFromRunoff(N,N1,N2,N4,N5,N4B,N5B)
  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B

  integer :: NN,K,NTG,NTN,idom
  !     begin_execution
  !     NET WATER, SNOW AND HEAT FLUXES FROM RUNOFF
  !
  !     TXGridSurfRunoff_2DH,TQS,TQW,TQI=net water and snowpack snow,water,ice runoff
  !     THeatXGridBySurfRunoff_2DH,THQS=net convective heat from surface water and snow runoff
  !     QR,HQR=runoff, convective heat from runoff from watsub.f
  !     QS,QW,QI=snow,water,ice transfer from watsub.f
  !     HQS=convective heat transfer from snow,water,ice transfer from watsub.f
  !
  D1202: DO NN=1,2
    !water flux
    D8590: DO K=1,micpar%NumOfLitrCmplxs
      DO idom=idom_beg,idom_end
        TOMQRS(idom,K,N2,N1)=TOMQRS(idom,K,N2,N1)+dom_2DFloXSurRunoff(idom,K,N,NN,N2,N1)
      ENDDO
    ENDDO D8590

    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
      !there is lateral runoff
      !water flux
      D8591: DO K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end
          TOMQRS(idom,K,N2,N1)=TOMQRS(idom,K,N2,N1)-dom_2DFloXSurRunoff(idom,K,N,NN,N5,N4)
        ENDDO
      ENDDO D8591

    ENDIF

    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      IF(IFLBH(N,NN,N5B,N4B).EQ.1)then

        D8592: DO K=1,micpar%NumOfLitrCmplxs
          DO idom=idom_beg,idom_end
            TOMQRS(idom,K,N2,N1)=TOMQRS(idom,K,N2,N1)-dom_2DFloXSurRunoff(idom,K,N,NN,N5B,N4B)
          enddo
        ENDDO D8592
      ENDIF

    ENDIF
  ENDDO D1202

  end subroutine OMH2OFluxesFromRunoff
!------------------------------------------------------------------------------------------

  subroutine TotalFluxFromSedmentTransp(N,N1,N2,N4,N5,N4B,N5B,NY,NX)
    implicit none
    integer, intent(in) :: N   !direction of calculation
    integer, intent(in) :: N1,N2,N4,N5,N4B,N5B,NY,NX

    integer :: M,K,NO,NN,NGL,NTX,NTP,MID,NE,idom
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
!       :PCPM,PCPD,PCPH=precip CaH4P2O8,CaHPO4,apatite in non-band
!       :PCPMB,PCPDB,PCPHB=precip CaH4P2O8,CaHPO4,apatite in band
!
  IF(N.NE.3.AND.(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros))THEN
    !horizontal direction
    D9350: DO NN=1,2
      IF(ABS(cumSedErosion(N,NN,N2,N1)).GT.ZEROS(N2,N1) &
        .OR.ABS(cumSedErosion(N,NN,N5,N4)).GT.ZEROS(N5,N4))THEN
        !incoming from south or east grid 
        tErosionSedmLoss(N2,N1)=tErosionSedmLoss(N2,N1)+cumSedErosion(N,NN,N2,N1)
        TSANER(N2,N1)             = TSANER(N2,N1)+XSANER(N,NN,N2,N1)
        TSILER(N2,N1)             = TSILER(N2,N1)+XSILER(N,NN,N2,N1)
        TCLAER(N2,N1)             = TCLAER(N2,N1)+XCLAER(N,NN,N2,N1)
        TNH4Eros_col(N2,N1)       = TNH4Eros_col(N2,N1)+XNH4ER(N,NN,N2,N1)
        TNH3Eros_col(N2,N1)       = TNH3Eros_col(N2,N1)+XNH3ER(N,NN,N2,N1)
        TNUreaEros_col(N2,N1)     = TNUreaEros_col(N2,N1)+XNHUER(N,NN,N2,N1)
        TNO3Eros_col(N2,N1)       = TNO3Eros_col(N2,N1)+XNO3ER(N,NN,N2,N1)
        TNH4ErosBand_col(N2,N1)   = TNH4ErosBand_col(N2,N1)+XNH4EB(N,NN,N2,N1)
        TNH3ErosBand_col(N2,N1)   = TNH3ErosBand_col(N2,N1)+XNH3EB(N,NN,N2,N1)
        TNUreaErosBand_col(N2,N1) = TNUreaErosBand_col(N2,N1)+XNHUEB(N,NN,N2,N1)
        TNO3ErosBand_col(N2,N1)   = TNO3ErosBand_col(N2,N1)+XNO3EB(N,NN,N2,N1)

        DO NTX=idx_beg,idx_end
          trcx_TER(NTX,N2,N1)=trcx_TER(NTX,N2,N1)+trcx_XER(NTX,N,NN,N2,N1)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          trcp_TER(NTP,N2,N1)=trcp_TER(NTP,N2,N1)+trcp_ER(NTP,N,NN,N2,N1)
        ENDDO

        DO  K=1,jcplx
          DO NO=1,NumMicbFunGrupsPerCmplx
            DO NGL=JGnio(NO),JGnfo(NO)
              DO M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  TOMEERhetr(NE,MID,K,N2,N1)=TOMEERhetr(NE,MID,K,N2,N1)+OMEERhetr(NE,MID,K,N,NN,N2,N1)
                ENDDO  
              enddo
            enddo
          enddo
        ENDDO

        DO NO=1,NumMicbFunGrupsPerCmplx
          DO NGL=JGniA(NO),JGnfA(NO)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                TOMEERauto(NE,MID,N2,N1)=TOMEERauto(NE,MID,N2,N1)+OMEERauto(NE,MID,N,NN,N2,N1)
              ENDDO  
            enddo
          enddo
        enddo


        D9375: DO K=1,jcplx
          D9370: DO M=1,ndbiomcp
            DO NE=1,NumPlantChemElms
              TORMER(NE,M,K,N2,N1)=TORMER(NE,M,K,N2,N1)+ORMER(NE,M,K,N,NN,N2,N1)
            ENDDO
          ENDDO D9370
          DO idom=idom_beg,idom_end
            TOHMER(idom,K,N2,N1)=TOHMER(idom,K,N2,N1)+OHMER(idom,K,N,NN,N2,N1)
          ENDDO
          D9365: DO M=1,jsken
            TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)+OSAER(M,K,N,NN,N2,N1)
            DO NE=1,NumPlantChemElms
              TOSMER(NE,M,K,N2,N1)=TOSMER(NE,M,K,N2,N1)+OSMER(NE,M,K,N,NN,N2,N1) 
            ENDDO
          ENDDO D9365
        ENDDO D9375

!     IF(NN.EQ.2)THEN
!       outgoing
        tErosionSedmLoss(N2,N1)=tErosionSedmLoss(N2,N1)-cumSedErosion(N,NN,N5,N4)
        TSANER(N2,N1)             = TSANER(N2,N1)-XSANER(N,NN,N5,N4)
        TSILER(N2,N1)             = TSILER(N2,N1)-XSILER(N,NN,N5,N4)
        TCLAER(N2,N1)             = TCLAER(N2,N1)-XCLAER(N,NN,N5,N4)
        TNH4Eros_col(N2,N1)       = TNH4Eros_col(N2,N1)-XNH4ER(N,NN,N5,N4)
        TNH3Eros_col(N2,N1)       = TNH3Eros_col(N2,N1)-XNH3ER(N,NN,N5,N4)
        TNUreaEros_col(N2,N1)     = TNUreaEros_col(N2,N1)-XNHUER(N,NN,N5,N4)
        TNO3Eros_col(N2,N1)       = TNO3Eros_col(N2,N1)-XNO3ER(N,NN,N5,N4)
        TNH4ErosBand_col(N2,N1)   = TNH4ErosBand_col(N2,N1)-XNH4EB(N,NN,N5,N4)
        TNH3ErosBand_col(N2,N1)   = TNH3ErosBand_col(N2,N1)-XNH3EB(N,NN,N5,N4)
        TNUreaErosBand_col(N2,N1) = TNUreaErosBand_col(N2,N1)-XNHUEB(N,NN,N5,N4)
        TNO3ErosBand_col(N2,N1)   = TNO3ErosBand_col(N2,N1)-XNO3EB(N,NN,N5,N4)

        DO NTX=idx_beg,idx_end
          trcx_TER(NTX,N2,N1)=trcx_TER(NTX,N2,N1)-trcx_XER(NTX,N,NN,N5,N4)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          trcp_TER(NTP,N2,N1)=trcp_TER(NTP,N2,N1)-trcp_ER(NTP,N,NN,N5,N4)
        ENDDO

        DO  K=1,jcplx
          DO  NO=1,NumMicbFunGrupsPerCmplx
            DO NGL=JGnio(NO),JGnfo(NO)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  TOMEERhetr(NE,MID,K,N2,N1)=TOMEERhetr(NE,MID,K,N2,N1)-OMEERhetr(NE,MID,K,N,NN,N5,N4)
                ENDDO  
              enddo
            enddo
          enddo
        ENDDO


        DO  NO=1,NumMicbFunGrupsPerCmplx
          DO  M=1,nlbiomcp
            DO NGL=JGniA(NO),JGnfA(NO)
              MID=micpar%get_micb_id(M,NGL)   
              DO NE=1,NumPlantChemElms         
                TOMEERauto(NE,MID,N2,N1)=TOMEERauto(NE,MID,N2,N1)-OMEERauto(NE,MID,N,NN,N5,N4)
              ENDDO  
            enddo
          enddo
        enddo

        D7375: DO K=1,jcplx
          D7370: DO M=1,ndbiomcp
            DO NE=1,NumPlantChemElms
              TORMER(NE,M,K,N2,N1)=TORMER(NE,M,K,N2,N1)-ORMER(NE,M,K,N,NN,N5,N4)
            ENDDO
          ENDDO D7370
          DO idom=idom_beg,idom_end
            TOHMER(idom,K,N2,N1)=TOHMER(idom,K,N2,N1)-OHMER(idom,K,N,NN,N5,N4)
          ENDDO
          D7365: DO M=1,jsken
            TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)-OSAER(M,K,N,NN,N5,N4)   
            DO NE=1,NumPlantChemElms       
              TOSMER(NE,M,K,N2,N1)=TOSMER(NE,M,K,N2,N1)-OSMER(NE,M,K,N,NN,N5,N4)
            ENDDO
          ENDDO D7365
        ENDDO D7375
!     ENDIF
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        IF(ABS(cumSedErosion(N,NN,N5B,N4B)).GT.ZEROS(N5,N4))THEN
          tErosionSedmLoss(N2,N1)   = tErosionSedmLoss(N2,N1)-cumSedErosion(N,NN,N5B,N4B)
          TSANER(N2,N1)             = TSANER(N2,N1)-XSANER(N,NN,N5B,N4B)
          TSILER(N2,N1)             = TSILER(N2,N1)-XSILER(N,NN,N5B,N4B)
          TCLAER(N2,N1)             = TCLAER(N2,N1)-XCLAER(N,NN,N5B,N4B)
          TNH4Eros_col(N2,N1)       = TNH4Eros_col(N2,N1)-XNH4ER(N,NN,N5B,N4B)
          TNH3Eros_col(N2,N1)       = TNH3Eros_col(N2,N1)-XNH3ER(N,NN,N5B,N4B)
          TNUreaEros_col(N2,N1)     = TNUreaEros_col(N2,N1)-XNHUER(N,NN,N5B,N4B)
          TNO3Eros_col(N2,N1)       = TNO3Eros_col(N2,N1)-XNO3ER(N,NN,N5B,N4B)
          TNH4ErosBand_col(N2,N1)   = TNH4ErosBand_col(N2,N1)-XNH4EB(N,NN,N5B,N4B)
          TNH3ErosBand_col(N2,N1)   = TNH3ErosBand_col(N2,N1)-XNH3EB(N,NN,N5B,N4B)
          TNUreaErosBand_col(N2,N1) = TNUreaErosBand_col(N2,N1)-XNHUEB(N,NN,N5B,N4B)
          TNO3ErosBand_col(N2,N1)   = TNO3ErosBand_col(N2,N1)-XNO3EB(N,NN,N5B,N4B)

          DO NTX=idx_beg,idx_end
            trcx_TER(NTX,N2,N1)=trcx_TER(NTX,N2,N1)-trcx_XER(NTX,N,NN,N5B,N4B)
          ENDDO

          DO NTP=idsp_beg,idsp_end
            trcp_TER(NTP,N2,N1)=trcp_TER(NTP,N2,N1)-trcp_ER(NTP,N,NN,N5B,N4B)
          ENDDO

          D8380: DO K=1,jcplx
            DO  NO=1,NumMicbFunGrupsPerCmplx
              DO NGL=JGnio(NO),JGnfo(NO)
                DO  M=1,nlbiomcp
                  MID=micpar%get_micb_id(M,NGL)    
                  DO NE=1,NumPlantChemElms            
                    TOMEERhetr(NE,MID,K,N2,N1)=TOMEERhetr(NE,MID,K,N2,N1)-OMEERhetr(NE,MID,K,N,NN,N5B,N4B)
                  ENDDO  
                enddo
              enddo
            enddo
          ENDDO D8380

          DO  NO=1,NumMicbFunGrupsPerCmplx
            DO NGL=JGniA(NO),JGnfA(NO)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)      
                DO NE=1,NumPlantChemElms        
                  TOMEERauto(NE,MID,N2,N1)=TOMEERauto(NE,MID,N2,N1)-OMEERauto(NE,MID,N,NN,N5B,N4B)
                ENDDO  
              enddo
            enddo
          enddo

          D8375: DO K=1,jcplx
            D8370: DO M=1,ndbiomcp
              DO NE=1,NumPlantChemElms
                TORMER(NE,M,K,N2,N1)=TORMER(NE,M,K,N2,N1)-ORMER(NE,M,K,N,NN,N5B,N4B)
              ENDDO
            ENDDO D8370
            DO idom=idom_beg,idom_end
              TOHMER(idom,K,N2,N1)=TOHMER(idom,K,N2,N1)-OHMER(idom,K,N,NN,N5B,N4B)
            ENDDO
            D8365: DO M=1,jsken
              TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)-OSAER(M,K,N,NN,N5B,N4B)          
              DO NE=1,NumPlantChemElms  
                TOSMER(NE,M,K,N2,N1)=TOSMER(NE,M,K,N2,N1)-OSMER(NE,M,K,N,NN,N5B,N4B)
              ENDDO
            ENDDO D8365
          ENDDO D8375
        ENDIF
      ENDIF
    ENDDO D9350
  ENDIF
  end subroutine TotalFluxFromSedmentTransp
!------------------------------------------------------------------------------------------

  subroutine FluxThruGrids(I,J,N,N1,N2,N3,N4,N5,N6,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N          !exchagne along direction, 1 east-west, 2 north-south, 3 vertical
  integer, intent(in) :: NY,NX      !geophysical location
  integer, intent(in) :: N1,N2,N3   !source grid indices
  integer, intent(in) :: N4,N5      !dest grid indices  
  integer, intent(inout) :: N6
  integer :: LL,K,NTSA,NTS,NTG,idom

  !     begin_execution
  !     NET HEAT, WATER FLUXES BETWEEN ADJACENT
  !     GRID CELLS
  !
  !     TFLW,TWatFlowCellMacP_vr,TWatFlowCellMacP_vr=net micropore,macropore water flux, heat flux
  !     FLW,FLWH,HFLW=micropore,macropore water flux, heat flux from watsub.f
  !     LakeSurfFlowMicP,FLWHNU,LakeSurfHeatFlux=lake surface water flux, heat flux from watsub.f if lake surface disappears
  !when FlowDirIndicator /=3, it means lateral exchange is consdiered
  !N==3 means vertical direction

  IF(FlowDirIndicator(N2,N1).NE.3 .OR. N.EQ.iVerticalDirection)THEN
    !locate the vertical layer for the dest grid
    D1200: DO LL=N6,NL(N5,N4)
      !modify the dest grid vertical location if needed
      !by matching the vertical layer number between source and dest, if N/=3
      !if N==3, skip insignificant layers
      IF(VLSoilPoreMicP_vr(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
        N6=LL
        exit
      ENDIF
    ENDDO D1200

    IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(N3.EQ.NU(N2,N1) .AND. N.EQ.3)THEN      
        !vertical direction, source is at soil surface
        TWatFlowCellMicP_vr(N3,N2,N1)  = TWatFlowCellMicP_vr(N3,N2,N1)+WaterFlowSoiMicP_3D(N,N3,N2,N1)-LakeSurfFlowMicP_col(N5,N4)
        TWatFlowCellMicPX_vr(N3,N2,N1) = TWatFlowCellMicPX_vr(N3,N2,N1)+WaterFlowSoiMicPX(N,N3,N2,N1)-LakeSurfFlowMicPX(N5,N4)
        TWatFlowCellMacP_vr(N3,N2,N1)  = TWatFlowCellMacP_vr(N3,N2,N1)+WaterFlowMacP_3D(N,N3,N2,N1)-LakeSurfFlowMacP_col(N5,N4)
        THeatFlowCellSoil_vr(N3,N2,N1) = THeatFlowCellSoil_vr(N3,N2,N1)+HeatFlow2Soil_3D(N,N3,N2,N1)-LakeSurfHeatFlux_col(N5,N4)

        if(THeatFlowCellSoil_vr(N3,N2,N1)<-1.e10)then
          write(*,*)'THeatFlowCellSoil_vr(N3,N2,N1)+HeatFlow2Soil_3D(N,N3,N2,N1)-LakeSurfHeatFlux_col(N5,N4)',&
            THeatFlowCellSoil_vr(N3,N2,N1),HeatFlow2Soil_3D(N,N3,N2,N1),LakeSurfHeatFlux_col(N5,N4)
          write(*,*)'Ns=',N1,n2,n3,n4,n5
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ELSE
        TWatFlowCellMicP_vr(N3,N2,N1)  = TWatFlowCellMicP_vr(N3,N2,N1)+WaterFlowSoiMicP_3D(N,N3,N2,N1)-WaterFlowSoiMicP_3D(N,N6,N5,N4)
        TWatFlowCellMicPX_vr(N3,N2,N1) = TWatFlowCellMicPX_vr(N3,N2,N1)+WaterFlowSoiMicPX(N,N3,N2,N1)-WaterFlowSoiMicPX(N,N6,N5,N4)
        TWatFlowCellMacP_vr(N3,N2,N1)  = TWatFlowCellMacP_vr(N3,N2,N1)+WaterFlowMacP_3D(N,N3,N2,N1)-WaterFlowMacP_3D(N,N6,N5,N4)
        THeatFlowCellSoil_vr(N3,N2,N1) = THeatFlowCellSoil_vr(N3,N2,N1)+HeatFlow2Soil_3D(N,N3,N2,N1)-HeatFlow2Soil_3D(N,N6,N5,N4)

        if(THeatFlowCellSoil_vr(N3,N2,N1)<-1.e10)then
          write(*,*)'THeatFlowCellSoil_vr(N3,N2,N1)+HeatFlow2Soil_3D(N,N3,N2,N1)-HeatFlow2Soil_3D(N,N6,N5,N4)',&
            THeatFlowCellSoil_vr(N3,N2,N1),HeatFlow2Soil_3D(N,N3,N2,N1),HeatFlow2Soil_3D(N,N6,N5,N4)
          write(*,*)'Ns=',N,N1,n2,n3,n4,n5,n6
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ENDIF
      !     IF(N1.EQ.1.AND.N3.EQ.1)THEN
      !     WRITE(*,6632)'TFLW',I,J,N,N1,N2,N3,N4,N5,N6,NU(N2,N1)
      !    2,TWatFlowCellMicP_vr(N3,N2,N1),WaterFlowSoiMicP_3D(N,N3,N2,N1),WaterFlowSoiMicP_3D(N,N6,N5,N4),LakeSurfFlowMicP_col(N5,N4)
      !    3,THeatFlowCellSoil_vr(N3,N2,N1),HeatFlow2Soil_3D(N,N3,N2,N1),HeatFlow2Soil_3D(N,N6,N5,N4)
      !    2,LakeSurfHeatFlux_col(N5,N4),VLWatMicP_vr(N3,N2,N1)
!6632  FORMAT(A8,10I4,12E16.8)
      !     ENDIF
      !
      !     NET SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
      !
      !     T*FLS,T*FHS=net convective+diffusive solute flux through micropores,macropores
      !     X*FLS,X*FHS=convective+diffusive solute flux through micropores, macropores from TranspNoSalt.f
      !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
      !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
      !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
      !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
      !
      D8585: DO K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Transp2Micp_vr(idom,K,N3,N2,N1)=DOM_Transp2Micp_vr(idom,K,N3,N2,N1) &
            +DOM_MicpTransp_3D(idom,K,N,N3,N2,N1)-DOM_MicpTransp_3D(idom,K,N,N6,N5,N4)
          DOM_Transp2Macp_flx(idom,K,N3,N2,N1)=DOM_Transp2Macp_flx(idom,K,N3,N2,N1) &
            +DOM_3DMacp_Transp_flx(idom,K,N,N3,N2,N1)-DOM_3DMacp_Transp_flx(idom,K,N,N6,N5,N4)
        enddo
      ENDDO D8585

      DO NTS=ids_beg,ids_end
        trcs_TransptMicP_vr(NTS,N3,N2,N1)=trcs_TransptMicP_vr(NTS,N3,N2,N1) &
          +trcs_TransptMicP_3D(NTS,N,N3,N2,N1)-trcs_TransptMicP_3D(NTS,N,N6,N5,N4)
        trcs_TransptMacP_vr(NTS,N3,N2,N1)=trcs_TransptMacP_vr(NTS,N3,N2,N1) &
          +trcs_TransptMacP_3D(NTS,N,N3,N2,N1)-trcs_TransptMacP_3D(NTS,N,N6,N5,N4)
      ENDDO
!
      !     NET GAS FLUXES BETWEEN ADJACENT GRID CELLS
      !
      !     T*FLG=net convective+diffusive gas flux
      !     X*FLG=convective+diffusive gas flux from TranspNoSalt.f
      !     gas code:*CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*HG*=H2
!exclude NH3B
      DO NTG=idg_beg,idg_NH3
        Gas_AdvDif_Flx_vr(NTG,N3,N2,N1)=Gas_AdvDif_Flx_vr(NTG,N3,N2,N1) &
          +Gas_3DAdvDif_Flx_vr(NTG,N,N3,N2,N1)-Gas_3DAdvDif_Flx_vr(NTG,N,N6,N5,N4)
      ENDDO
!
      !     NET SALT FLUXES BETWEEN ADJACENT GRID CELLS
      !
      !     T*FLS,T*FHS=net convective+diffusive solute flux through micropores,macropores
      !     X*FLS,X*FHS=convective+diffusive solute flux through micropores, macropores from TranspSalt.f
      !     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
      !          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
      !          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
      !          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
      !          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
      !          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
      !          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
      !     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
      !          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
      !          :*1=non-band,*B=band
!
      IF(salt_model)THEN
        DO NTSA=idsalt_beg,idsaltb_end
          trcSalt_Flo2MicP_vr(NTSA,N3,N2,N1)=trcSalt_Flo2MicP_vr(NTSA,N3,N2,N1) &
            +trcSalt3DFlo2Cell(NTSA,N,N3,N2,N1)-trcSalt3DFlo2Cell(NTSA,N,N6,N5,N4)
          trcSalt_Flo2MacP_vr(NTSA,N3,N2,N1)=trcSalt_Flo2MacP_vr(NTSA,N3,N2,N1) &
            +trcSalt_XFHS(NTSA,N,N3,N2,N1)-trcSalt_XFHS(NTSA,N,N6,N5,N4)
        ENDDO
      ENDIF
    ELSE
      TWatFlowCellMicP_vr(N3,N2,N1)  = 0.0_r8
      TWatFlowCellMicPX_vr(N3,N2,N1) = 0.0_r8
      TWatFlowCellMacP_vr(N3,N2,N1)  = 0.0_r8
      THeatFlowCellSoil_vr(N3,N2,N1) = 0.0_r8
      WatIceThawMicP_vr(N3,N2,N1)    = 0.0_r8
      WatIceThawMacP_vr(N3,N2,N1)    = 0.0_r8
      THeatSoiThaw_vr(N3,N2,N1)      = 0.0_r8

      D8596: DO K=1,jcplx
        DOM_Transp2Micp_vr(idom_beg:idom_end,K,N3,N2,N1)  = 0.0_r8
        DOM_Transp2Macp_flx(idom_beg:idom_end,K,N3,N2,N1) = 0.0_r8
      ENDDO D8596

      trcs_TransptMicP_vr(ids_beg:ids_end,N3,N2,N1) = 0.0_r8
      trcs_TransptMacP_vr(ids_beg:ids_end,N3,N2,N1) = 0.0_r8
      Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,N3,N2,N1)   = 0.0_r8

      IF(salt_model)THEN
        trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,N3,N2,N1) = 0.0_r8
        trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,N3,N2,N1) = 0.0_r8
      ENDIF
    ENDIF
  ENDIF
  end subroutine FluxThruGrids

end module LateralTranspMod
