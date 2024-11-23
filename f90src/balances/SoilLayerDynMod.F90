module SoilLayerDynMod
! Description:
! subroutines to do soil relayering
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only : micpar
  use RootDataType
  use GridDataType
  USE EcoSIMCtrlDataType
  USE SoilWaterDataType
  use SoilPropertyDataType
  use FlagDataType
  use SOMDataType
  use MicrobialDataType
  USE SoilBGCDataType
  use AqueChemDatatype
  use FertilizerDataType
  use SoilHeatDataType
  use EcoSIMCtrlMod
  use SoilPhysDataType
  use SurfSoilDataType
  use SurfLitterDataType
  use EcosimConst
  use minimathmod, only : AZMAX1
  use UnitMod     , only : units
  use EcoSIMConfig, only : ndbiomcp => NumDeadMicrbCompts
  USE TFlxTypeMod , ONLY : tErosionSedmLoss,TDLYXF,TDayLenthPrevC,TDVOLI,TDORGC
implicit none


  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  real(r8) :: XVOLWP,WDPOBDL,WDNOBD1,WDPOBD0,WDPOBD1
  real(r8) :: WDNHBDL,WDNHBD0,WDNHBD1,WDNOBDL,WDNOBD0
  real(r8), PARAMETER :: ZEROC=0.1E-03_r8
  integer, parameter :: ist_water=0
  integer, parameter :: ist_soil=1
  integer, parameter :: ich_watlev =1
  integer, parameter :: ich_frzthaw=4
  integer, parameter :: ich_erosion=5
  integer, parameter :: ich_somloss=6

  public :: RelayerSoilProfile
  contains

  subroutine RelayerSoilProfile(NY,NX,DORGC,DVLiceMicP_vr)
  !
  !Description:
  !relayer the soil profiles
  implicit none
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: DORGC(JZ)  !change in organic matter, initial-final
  REAL(R8),INTENT(IN) :: DVLiceMicP_vr(JZ)  !change in ice volume, initial-final

  real(r8) :: CDPTHY(0:JZ),CDPTHX(JZ)
  integer :: IFLGL(0:JZ,6)  !flag for soil thickness change
  real(r8) :: DDLYRX(3)
  integer :: NN,K,M,N,NR,NZ,L
  integer :: L0,L1,NUX,ICHKL,NGL
  real(r8) :: FX,FY
  real(r8) :: FBO
  real(r8) :: FHO,FWO,FXVOLW
  real(r8) :: FO
  real(r8) :: ENGY0X,ENGY0,ENGY1X,ENGY1
  !     begin_execution
  !     SOIL SUBSIDENCE
  !
  if(lverb)write(*,*)'RelayerSoilProfile'
  IF(.not. erosion_model)return
  !soil relayering can occur due to freeze-thaw, soc change, and erosion

  call SoilRelayering(NY,NX,DORGC,DVLiceMicP_vr,CDPTHX,CDPTHY,IFLGL)

    !
    !     RECALCULATE SOIL LAYER THICKNESS
    !
  ICHKL=0
  D245: DO L=NU(NY,NX),NL(NY,NX)-1
    D230: DO NN=1,3

      call UpdateLayerThickness(NN,L,NY,NX,ICHKL,NUX,CDPTHX,CDPTHY,IFLGL,DDLYRX)

      !
      !     TRANSFER STATE VARIABLES BETWEEN LAYERS
      !
      !     IF(IFLGL(L,NN).EQ.1)THEN
      IF(ABS(DDLYRX(NN)).GT.ZERO)THEN
!L0,L1: target and source layers
        call getFLs(L,NN,NY,NX,NUX,DDLYRX,IFLGL,FX,FO,L1,L0)

        IF(FX.GT.ZERO)THEN
          iResetSoilProf_col(NY,NX)=itrue
          FY=1.0_r8-FX
          IF(FY.LE.ZERO2)FY=0.0_r8
          IF(SoiBulkDensity_vr(L0,NY,NX).LE.ZERO)THEN
!     TARGET POND LAYER
!
            call tgtPondLyr(L,L0,L1,NY,NX,NN,FX,FY,CDPTHY,IFLGL)
          ELSE
!
      !     MOVE ALL MATTER WITH CHANGES IN LAYER DEPTHS
      !
            IF(L0.NE.0)THEN
              call MoveFertMinerals(L,L0,L1,NY,NX,FX,FO)
            ENDIF
            call MoveHeatWat(L,L0,L1,NY,NX,FO,FX)

!
!     SOIL FERTILIZER
!
            call MoveFertSalt(L,L0,L1,NY,NX,FX,FWO)
!
!     SOIL SOLUTES IN BAND
!
            IF(L0.NE.0)THEN
              call MoveBandSolute(L,L0,L1,NY,NX,FX,FWO,FO)

              call MoveDisolvGas(L0,L1,NY,NX,FX,FWO)

              call MoveMacPoreSolute(L0,L1,NY,NX,FHO)
            ENDIF
!     SOIL ORGANIC MATTER
            call MoveSOM(L0,L1,L,NY,NX,FO,IFLGL)

            IF(NN.EQ.1)THEN
              IF(SoiBulkDensity_vr(L0,NY,NX).LE.ZERO.AND.SoiBulkDensity_vr(L1,NY,NX).LE.ZERO &
                .AND.VLWatMicP_vr(L0,NY,NX)+VLiceMicP_vr(L0,NY,NX).LE.ZEROS(NY,NX))THEN
                CumDepz2LayerBot_vr(L1,NY,NX) = CumDepz2LayerBot_vr(L0,NY,NX)
                CDPTHY(L1)                    = CDPTHY(L0)
              ENDIF
            ENDIF

          ENDIF
!
!     RESET LOWER LAYER NUMBER WITH EROSION
!
          IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
            IF(L.EQ.NL(NY,NX).AND.DLYR(3,L,NY,NX).GT.DLYRI_3D(3,L,NY,NX))THEN
              NL(NY,NX)=MIN(NLI(NY,NX),NL(NY,NX)+1)
            ENDIF
            IF(L.EQ.NL(NY,NX)-1.AND.CumDepz2LayerBot_vr(NL(NY,NX),NY,NX)-CumDepz2LayerBot_vr(L,NY,NX).LE.ZEROC)THEN
              CumDepz2LayerBot_vr(L,NY,NX)         = CumDepz2LayerBot_vr(L,NY,NX)+DLYR(3,NL(NY,NX),NY,NX)
              CDPTHY(L)                            = CDPTHY(L)+DLYR(3,NL(NY,NX),NY,NX)
              CumDepz2LayerBot_vr(NL(NY,NX),NY,NX) = CumDepz2LayerBot_vr(L,NY,NX)
              CDPTHY(NL(NY,NX))                    = CDPTHY(L)
              DLYR(3,NL(NY,NX),NY,NX)              = 0.0_r8
              NL(NY,NX)                            = L
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D230
  ENDDO D245
  end subroutine RelayerSoilProfile
!------------------------------------------------------------------------------------------
  subroutine SoilRelayering(NY,NX,DORGC,DVLiceMicP_vr,CDPTHX,CDPTHY,IFLGL)
!
! IFLGL: c1, ponding water, c2, pond disappear, c3, pond reappare, c4: freeze-thaw, c5: erosion, c6: som change
  implicit none

  integer, intent(in) :: NY,NX                 !column location
  REAL(R8),INTENT(IN) :: DVLiceMicP_vr(JZ)        !change in ice volume, final - initial
  real(r8),intent(in) :: DORGC(JZ)             !change in SOM, initial - final
  real(r8), intent(out) :: CDPTHX(JZ)          !copy of the old depths of each layer
  real(r8), intent(out) :: CDPTHY(0:JZ)        !
  integer , intent(out) :: IFLGL(0:JZ,6)       !flag of change type
  real(r8) :: DDLYX(0:JZ,6)                    !six process
  real(r8) :: DDLYR(0:JZ,6)                    !new layer thickness
  integer :: LX,LY,LL,NN,L
  real(r8) :: DDLEqv_OrgC,DDLEqv_Erosion
  real(r8) :: DDLWatEqv_IceMicP
  real(r8) :: DENSJ,DLEqv_MicP
  real(r8) :: DLYR_XMicP   !layer thickness excluding micropore
  integer  :: itoplyr_type          !surface layer type: 0 water, 1 soil  
! begin_execution

  !
  IF(SoiBulkDensity_vr(NU(NY,NX),NY,NX).LE.ZERO)THEN
    itoplyr_type=ist_water      !surface is water layer
  ELSE
    !it is a soil column
    itoplyr_type=ist_soil
  ENDIF

  DDLYX=0._r8
! starting from bottom up
  D225: DO LX=NL(NY,NX),NU(NY,NX),-1
    !make a copy of the depth, bottom of the layer
    CDPTHX(LX)  = CumDepz2LayerBot_vr(LX,NY,NX)
    CDPTHY(LX)  = CumDepz2LayerBot_vr(LX,NY,NX)
    IFLGL(LX,:) = 0
    !
    !     POND, from water to soil
    ! a layer is made up of soil micropores+macropores and solid materials (can be zero).
    ! compute the change
    IF(SoiBulkDensity_vr(LX,NY,NX).LE.ZERO)THEN  !current layer is water

      ! next layer is soil, or top layer is soil, i.e. water ponding on soil
      IF(SoiBulkDensity_vr(LX+1,NY,NX).GT.ZERO .OR. itoplyr_type.EQ.ist_soil)THEN        
        DLYR_XMicP           = DLYR(3,LX,NY,NX)-(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA(3,LX,NY,NX)
        DDLYX(LX,ich_watlev) = DLYR_XMicP+DDLYX(LX+1,ich_watlev)   !combined thickness current + next layer
        DDLYR(LX,ich_watlev) = DDLYX(LX+1,ich_watlev)              !make a copy of the next layer
        IFLGL(LX,ich_watlev) = 2                                        !
      ELSE        
        !next, current and top layers are all water, ponding water
        !DLYR_XMicP: non-micropore soil equivalent depth
        !DLYRI: initial water layer thickness, [m]

        DLYR_XMicP=DLYRI_3D(3,LX,NY,NX)-(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA(3,LX,NY,NX)

        !DLEqv_MicP: water+ice total thickness of next layer
        DLEqv_MicP=(VLWatMicP_vr(LX+1,NY,NX)+VLiceMicP_vr(LX+1,NY,NX))/AREA(3,LX,NY,NX)

        !there is expansion in layer LX, or next layer has water + ice
        IF(DLYR_XMicP.LT.-ZERO.OR.DLEqv_MicP.GT.ZERO)THEN
          DDLYX(LX,ich_watlev) = DLYR_XMicP+DDLYX(LX+1,ich_watlev)     !move next layer upward
          DDLYR(LX,ich_watlev) = AMIN1(DDLYX(LX+1,ich_watlev),DLEqv_MicP)   !

          IF(DLEqv_MicP.GT.ZERO)THEN
            ! already has significant micropore volume/thickness
            IFLGL(LX,ich_watlev)=1
          ELSE
            ! has no solid thickness
            IFLGL(LX,ich_watlev)=2
          ENDIF
        ELSE
          !shrink
          DLYR_XMicP           = DLYR(3,LX,NY,NX)-(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA(3,LX,NY,NX)
          DDLYX(LX,ich_watlev) = DLYR_XMicP+DDLYX(LX+1,ich_watlev)    !combine next layer to current
          DDLYR(LX,ich_watlev) = DDLYX(LX+1,ich_watlev)
          IFLGL(LX,ich_watlev) = 2
        ENDIF
      ENDIF

      !surface layer or layer above is still soil
      IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity_vr(LX-1,NY,NX).GT.ZERO)THEN
        DDLYX(LX-1,ich_watlev) = DDLYX(LX,ich_watlev)         !move upward
        DDLYR(LX-1,ich_watlev) = DDLYX(LX,ich_watlev)
        IFLGL(LX-1,ich_watlev) = 1
      ENDIF
      DDLYX(LX,ich_frzthaw) = 0.0_r8
      DDLYR(LX,ich_frzthaw) = 0.0_r8
      IFLGL(LX,ich_frzthaw) = 0
      DDLYX(LX,ich_erosion) = 0.0_r8
      DDLYR(LX,ich_erosion) = 0.0_r8
      IFLGL(LX,ich_erosion) = 0
      DDLYX(LX,ich_somloss) = 0.0_r8
      DDLYR(LX,ich_somloss) = 0.0_r8
      IFLGL(LX,ich_somloss) = 0
      !
      !     SOIL
      !
    ELSE
      !current layer is soil
      !     FREEZE-THAW
      !DVLiceMicP: change in ice volume, initial-final
      !DENSICE: mass density of ice, g/cm3
      !DENSJ: volume added per unit addition of ice (with respect to ice)
      !the change is put to layer thickness
      !DDLWatEqv_IceMicP: added thickness change
      ! 4: due to freeze-thaw
      IF(ABS(DVLiceMicP_vr(LX)).GT.ZEROS(NY,NX))THEN
        DENSJ=1._r8-DENSICE
        DDLWatEqv_IceMicP=DVLiceMicP_vr(LX)*DENSJ/AREA(3,LX,NY,NX)
        !bottom layer
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,ich_frzthaw) = DDLWatEqv_IceMicP
          DDLYR(LX,ich_frzthaw) = 0.0_r8
          IFLGL(LX,ich_frzthaw) = 0
         !not bottom layer
        ELSE
          !add depth of two layers
          DDLYX(LX,ich_frzthaw)=DDLWatEqv_IceMicP+DDLYX(LX+1,ich_frzthaw)
          !make a copy 
          DDLYR(LX,ich_frzthaw)=DDLYX(LX+1,ich_frzthaw)
          !    2+DLYRI_3D(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
          IFLGL(LX,ich_frzthaw)=0

          !top soil layer or layer above is water
          IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity_vr(LX-1,NY,NX).LE.ZERO)THEN
            DDLYX(LX-1,ich_frzthaw) = DDLYX(LX,ich_frzthaw)
            DDLYR(LX-1,ich_frzthaw) = DDLYX(LX,ich_frzthaw)
            !    2+DLYRI_3D(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
            IFLGL(LX-1,ich_frzthaw)=0
          ENDIF
        ENDIF
      ELSE
        ! no change in ice volume
        DDLWatEqv_IceMicP=0.0_r8
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,ich_frzthaw) = 0.0_r8
          DDLYR(LX,ich_frzthaw) = 0.0_r8
          IFLGL(LX,ich_frzthaw) = 0
        ELSE
          DDLYX(LX,ich_frzthaw) = DDLYX(LX+1,ich_frzthaw)
          DDLYR(LX,ich_frzthaw) = DDLYX(LX+1,ich_frzthaw)
          IFLGL(LX,ich_frzthaw) = 0
          IF(LX.EQ.NU(NY,NX))THEN
            DDLYX(LX-1,ich_frzthaw) = DDLYX(LX,ich_frzthaw)
            DDLYR(LX-1,ich_frzthaw) = DDLYX(LX,ich_frzthaw)
            IFLGL(LX-1,ich_frzthaw) = 0
          ENDIF
        ENDIF
      ENDIF

      !total change in ice volume
!      TDVOLI=TDVOLI+DVLiceMicP_vr(LX)
!      TDLYXF=TDLYXF+DDLWatEqv_IceMicP
      !
      !     EROSION model is on
      !
      IF((iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros) &
        .AND.ABS(tErosionSedmLoss(NY,NX)).GT.ZEROS(NY,NX))THEN
        IF(LX.EQ.NL(NY,NX))THEN
!  5: due to sediment erosion
!         total soil layer reduction due to erosion
          DDLEqv_Erosion        = -tErosionSedmLoss(NY,NX)/(SoilMicPMassLayerMX(NY,NX)/VLSoilPoreMicP_vr(NU(NY,NX),NY,NX))
          DDLYX(LX,ich_erosion) = DDLEqv_Erosion
          DDLYR(LX,ich_erosion) = DDLEqv_Erosion
          IFLGL(LX,ich_erosion) = 1
        ELSE
          DDLYX(LX,ich_erosion) = 0.0_r8
          DDLYR(LX,ich_erosion) = 0.0_r8
          IFLGL(LX,ich_erosion) = 0
        ENDIF
      ELSE
          DDLYX(LX,ich_erosion) = 0.0_r8
          DDLYR(LX,ich_erosion) = 0.0_r8
          IFLGL(LX,ich_erosion) = 0
      ENDIF

      !
      ! SOC GAIN OR LOSS, i.e. DORGC(LX) significantly non-zero
      ! SoilFracAsMacP: macropore fraction
      ! DDLEqv_OrgC: soil thickness added due to change in organic matter,
      ! keeping macropore fraction
      ! SoiBulkDensityt0_vr: initial bulk density,
      IF((iErosionMode.EQ.ieros_frzthawsom.OR.iErosionMode.EQ.ieros_frzthawsomeros)&
        .AND.ABS(DORGC(LX)).GT.ZEROS(NY,NX))THEN
        DDLEqv_OrgC = MWC2Soil*DORGC(LX)/((1.0_r8-SoilFracAsMacP_vr(LX,NY,NX))*SoiBulkDensityt0_vr(LX,NY,NX))/AREA(3,LX,NY,NX)

        ! obtain diagnostics only for NX==1
!        IF(NX.EQ.1)THEN
!          TDORGC=TDORGC+DORGC(LX)
!          TDayLenthPrevC=TDayLenthPrevC+DDLEqv_OrgC
!        ENDIF

        ! bottom layer, or next layer is water
        IF(LX.EQ.NL(NY,NX).OR.SoiBulkDensity_vr(LX+1,NY,NX).LE.ZERO)THEN
          DDLYX(LX,ich_somloss) = DDLEqv_OrgC
          DDLYR(LX,ich_somloss) = 0.0_r8
          IFLGL(LX,ich_somloss) = 1
        ELSE
          !add total change from current layer to next
          DDLYX(LX,ich_somloss) = DDLEqv_OrgC+DDLYX(LX+1,ich_somloss)
          DDLYR(LX,ich_somloss) = DDLYX(LX+1,ich_somloss)+DLYRI_3D(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
          IFLGL(LX,ich_somloss) = 1
          ! top layer, the layer right above is water
          IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity_vr(LX-1,NY,NX).LE.ZERO)THEN
            DDLYX(LX-1,ich_somloss) = DDLYX(LX,ich_somloss)
            DDLYR(LX-1,ich_somloss) = DDLYX(LX,ich_somloss)
            !    2+DLYRI_3D(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
            IFLGL(LX-1,ich_somloss)=1
          ENDIF
        ENDIF
      ELSE
        ! bottom layer
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,ich_somloss) = 0.0_r8
          DDLYR(LX,ich_somloss) = 0.0_r8
          IFLGL(LX,ich_somloss) = 0
        ELSE
          DDLYX(LX,ich_somloss) = DDLYX(LX+1,ich_somloss)
          DDLYR(LX,ich_somloss) = DDLYX(LX+1,ich_somloss)
          IFLGL(LX,ich_somloss) = 0
        ENDIF
      ENDIF
      DDLYX(LX,ich_watlev) = 0.0_r8
      DDLYR(LX,ich_watlev) = 0.0_r8
      IFLGL(LX,ich_watlev) = 0

    ENDIF

    ! apply the change
      !
      !     RESET SOIL LAYER DEPTHS
      !
    D200: DO NN=1,6
      !     IF(ABS(DDLYX(LX,NN)).GT.ZERO)iResetSoilProf_col(NY,NX)=1
!     c2, and c3 are not implemented
      IF(NN.NE.2.AND.NN.NE.3)THEN
        !
        !     POND
        !
        IF(SoiBulkDensity_vr(LX,NY,NX).LE.ZERO)THEN
          ! there are some changes
          IF(IFLGL(LX,NN).NE.0)THEN
            CumDepz2LayerBot_vr(LX,NY,NX)=CumDepz2LayerBot_vr(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX)=CDPTHY(LX)+DDLYR(LX,NN)
            !  not top layer
            IF(LX.NE.NU(NY,NX).AND.IFLGL(LX,ich_watlev).EQ.2)THEN
              DO  LL=LX-1,0,-1
                CumDepz2LayerBot_vr(LL,NY,NX)=CumDepz2LayerBot_vr(LL,NY,NX)+DDLYX(LX,NN)
                CDPTHY(LL)=CDPTHY(LL)+DDLYX(LX,NN)
              ENDDO
              DDLYX(LX,NN)=0.0_r8
            ENDIF
            !  top layer
            IF(LX.EQ.NU(NY,NX))THEN
              CumDepz2LayerBot_vr(LX-1,NY,NX)=CumDepz2LayerBot_vr(LX,NY,NX) &
                -(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA(3,LX,NY,NX)
              CDPTHY(LX-1)=CDPTHY(LX)-(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA(3,LX,NY,NX)
            ENDIF
          ENDIF
          !
          !  SOIL
          !
        ELSE
          !     IF(DDLYR(L,NN).NE.0.OR.DDLYR(LX,NN).NE.0)THEN
          !
          !     FREEZE-THAW
          !
          IF(NN.EQ.ich_frzthaw)THEN
            CumDepz2LayerBot_vr(LX,NY,NX)=CumDepz2LayerBot_vr(LX,NY,NX)+DDLYR(LX,NN)
            !     CDPTHY(LX)=CDPTHY(LX)+DDLYR(LX,NN)
! top layer
            IF(LX.EQ.NU(NY,NX))THEN
              CumDepz2LayerBot_vr(LX-1,NY,NX)=CumDepz2LayerBot_vr(LX-1,NY,NX)+DDLYR(LX-1,NN)
              !     CDPTHY(LX-1)=CDPTHY(LX-1)+DDLYR(LX-1,NN)

            ENDIF
          ENDIF
            !
            ! SET SURFACE ELEVATION FOR SOIL EROSION
            !
          IF(NN.EQ.ich_erosion.AND.IFLGL(LX,NN).EQ.1)THEN
            CumDepz2LayerBot_vr(LX,NY,NX) = CumDepz2LayerBot_vr(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX)                    = CDPTHY(LX)+DDLYR(LX,NN)

            IF(LX.EQ.NU(NY,NX))THEN
              CumDepz2LayerBot_vr(LX-1,NY,NX) = CumDepz2LayerBot_vr(LX-1,NY,NX)+DDLYR(LX,NN)
              CDPTHY(LX-1)                    = CDPTHY(LX-1)+DDLYR(LX,NN)
            ENDIF
          ENDIF
          !
          ! SET SOIL LAYER DEPTHS FOR CHANGES IN SOC
          !
          IF(NN.EQ.ich_somloss.AND.IFLGL(LX,NN).EQ.1)THEN
            CumDepz2LayerBot_vr(LX,NY,NX) = CumDepz2LayerBot_vr(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX)                    = CDPTHY(LX)+DDLYR(LX,NN)

            IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity_vr(LX-1,NY,NX).LE.ZERO)THEN
              CumDepz2LayerBot_vr(LX-1,NY,NX) = CumDepz2LayerBot_vr(LX-1,NY,NX)+DDLYR(LX-1,NN)
              CDPTHY(LX-1)                    = CDPTHY(LX-1)+DDLYR(LX-1,NN)

              IF(SoiBulkDensity_vr(LX-1,NY,NX).LE.ZERO)THEN
                DO  LY=LX-2,0,-1
                  IF(SoiBulkDensity_vr(LY+1,NY,NX).LE.ZERO)THEN
                    CumDepz2LayerBot_vr(LY,NY,NX) = CumDepz2LayerBot_vr(LY,NY,NX)+DDLYR(LX-1,NN)
                    CDPTHY(LY)                    = CDPTHY(LY)+DDLYR(LX-1,NN)
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D200
    VLSoilMicP_vr(LX,NY,NX)=VLSoilPoreMicP_vr(LX,NY,NX)
  ENDDO D225
  VLSoilMicP_vr(0,NY,NX)=VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)
  end subroutine SoilRelayering

!------------------------------------------------------------------------------------------


  subroutine getFLs(L,NN,NY,NX,NUX,DDLYRX,IFLGL,FX,FO,L1,L0)

  implicit none
  integer, intent(in) :: L,NN,NY,NX,NUX
  real(r8),intent(in) :: DDLYRX(3)
  integer, intent(in) :: IFLGL(0:JZ,6)
  real(r8), intent(out) :: FX,FO
  integer, intent(out) :: L1   !source
  integer, intent(out) :: L0   !target
  real(r8) :: DLEqv_MicP

! begin_execution
  IF(DDLYRX(NN).GT.ZERO)THEN
!
!     LAYERS MOVE DOWN (FX>0.), pond drys up
!
    IF(IFLGL(L,2).EQ.0)THEN
      L1 = L
      L0 = L+1
    ELSE
      L1 = NU(NY,NX)
      L0 = NUX
    ENDIF

    IF((SoiBulkDensity_vr(L,NY,NX).LE.ZERO.AND.IFLGL(L,ich_watlev).EQ.2) &
        .OR.(DLYR(3,L0,NY,NX).LE.ZEROC.AND.IFLGL(L,ich_somloss).EQ.1))THEN
      FX = 1.0_r8
      FO = 1.0_r8
    ELSE
      IF(SoiBulkDensity_vr(L0,NY,NX).LE.ZERO)THEN
        DLEqv_MicP=(VLWatMicP_vr(L0,NY,NX)+VLiceMicP_vr(L0,NY,NX))/AREA(3,L0,NY,NX)
        IF(DLEqv_MicP.GT.ZERO)THEN
          FX = AMIN1(1.0_r8,DDLYRX(NN)/DLEqv_MicP)
          FO = FX
        ELSE
          FX = 0.0_r8
          FO = 0.0_r8
        ENDIF
      ELSE
        IF(DLYR(3,L,NY,NX).GT.ZERO.AND.DLYR(3,L0,NY,NX).GT.ZERO)THEN
          FX = AMIN1(1.0,DDLYRX(NN)/DLYR(3,L,NY,NX))
          FO = AMIN1(1.0,DDLYRX(NN)/DLYR(3,L0,NY,NX))
        ELSE
          FX = 0.0_r8
          FO = 0.0_r8
        ENDIF
      ENDIF
    ENDIF
  ELSE
!
!     LAYERS MOVE UP (FX<=0.), pond builds up
!
    IF(IFLGL(L,3).EQ.0)THEN
      L1 = L+1
      L0 = L
    ELSE
      L1 = NU(NY,NX)
      L0 = 0
    ENDIF
    IF(SoiBulkDensity_vr(L0,NY,NX).LE.ZERO)THEN
      DLEqv_MicP=(VLWatMicP_vr(L0,NY,NX)+VLiceMicP_vr(L0,NY,NX))/AREA(3,L0,NY,NX)
      IF(DLEqv_MicP.GT.ZERO)THEN
        FX = AMIN1(1.0,-DDLYRX(NN)/DLEqv_MicP)
        FO = FX
      ELSE
        FX = 0.0_r8
        FO = 0.0_r8
      ENDIF
    ELSE
      IF(DLYR(3,L,NY,NX).GT.ZERO.AND.DLYR(3,L0,NY,NX).GT.ZERO)THEN
        FX = AMIN1(1.0,-DDLYRX(NN)/DLYR(3,L,NY,NX))
        FO = AMIN1(1.0,-DDLYRX(NN)/DLYR(3,L0,NY,NX))
      ELSE
        FX = 0.0_r8
        FO = 0.0_r8
      ENDIF
    ENDIF
  ENDIF
  end subroutine getFLs
!------------------------------------------------------------------------------------------

  subroutine UpdateLayerThickness(NN,L,NY,NX,ICHKL,NUX,CDPTHX,CDPTHY,IFLGL,DDLYRX)
  implicit none
  integer, intent(in) :: NN,L,NY,NX
  integer, intent(inout) :: ICHKL
  integer, intent(out) :: NUX    !old top layer index
  real(r8), intent(in) :: CDPTHX(JZ)
  real(r8), intent(inout) :: CDPTHY(0:JZ)
  integer, intent(inout) :: IFLGL(0:JZ,6)
  real(r8), intent(out) :: DDLYRX(3)
  real(r8) :: DDLYRY(JZ)
  real(r8) :: DLYR0
  real(r8) :: DLYRXX
  integer :: LL

! begin_execution
  NUX=0
  IF(NN.EQ.1)THEN
    DLYR(3,L,NY,NX) = CumDepz2LayerBot_vr(L,NY,NX)-CumDepz2LayerBot_vr(L-1,NY,NX) !current layer depth
    DLYRXX          = DLYR(3,L,NY,NX)
    IF(IFLGL(L,ich_watlev).EQ.0.AND.IFLGL(L+1,ich_watlev).NE.0)THEN
      DDLYRX(NN)=0.0_r8
      IF(SoiBulkDensity_vr(L,NY,NX).LE.ZERO)THEN
        DDLYRY(L)=DLYRI_3D(3,L,NY,NX)-DLYR(3,L,NY,NX)
      ELSE
        DDLYRY(L)=0.0_r8
      ENDIF
      ICHKL=1
    ELSEIF(IFLGL(L,ich_watlev).EQ.2.AND.(IFLGL(L+1,ich_watlev).EQ.0 .OR. DLYR(3,L,NY,NX).LE.DLYRI_3D(3,L,NY,NX)))THEN
      DDLYRX(NN)=0.0_r8
      IF(L.EQ.NU(NY,NX).OR.ICHKL.EQ.0)THEN
        DDLYRY(L)=0.0_r8
      ELSE
        DDLYRY(L)=DDLYRY(L-1)
      ENDIF
      IF(IFLGL(L,ich_watlev).EQ.2.AND.IFLGL(L+1,ich_watlev).EQ.0)ICHKL=0
    ELSE
      IF(ICHKL.EQ.0)THEN
        DDLYRX(NN) = DLYRI_3D(3,L,NY,NX)-DLYR(3,L,NY,NX)
        DDLYRY(L)  = DDLYRX(NN)
      ELSE
        DDLYRX(NN) = 0.0_r8
        DDLYRY(L)  = DDLYRY(L-1)
      ENDIF
    ENDIF
    CumDepz2LayerBot_vr(L,NY,NX)=CumDepz2LayerBot_vr(L,NY,NX)+DDLYRY(L)
  !     CDPTHY(L)=CDPTHY(L)+DDLYRY(L)
    DLYR(3,L,NY,NX)              = DLYR(3,L,NY,NX)+DDLYRY(L)
    SoiDepthMidLay_vr(L,NY,NX)   = 0.5_r8*(CumDepz2LayerBot_vr(L,NY,NX)+CumDepz2LayerBot_vr(L-1,NY,NX))
    CumSoilThickness_vr(L,NY,NX) = CumDepz2LayerBot_vr(L,NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
    IF(L.EQ.NL(NY,NX)-1)THEN
      DLYR(3,L+1,NY,NX)              = CumDepz2LayerBot_vr(L+1,NY,NX)-CumDepz2LayerBot_vr(L,NY,NX)
      SoiDepthMidLay_vr(L+1,NY,NX)   = 0.5_r8*(CumDepz2LayerBot_vr(L+1,NY,NX)+CumDepz2LayerBot_vr(L,NY,NX))
      CumSoilThickness_vr(L+1,NY,NX) = CumDepz2LayerBot_vr(L+1,NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)
    ENDIF
    IF(L.EQ.NU(NY,NX))THEN
      DPTHZ_vr(L,NY,NX)=0.5_r8*CumSoilThickness_vr(L,NY,NX)
    !     DDLYRX(NN)=DDLYRX(NN)+DDLYR(L,5)
    ELSE
      DPTHZ_vr(L,NY,NX)=0.5_r8*(CumSoilThickness_vr(L,NY,NX)+CumSoilThickness_vr(L-1,NY,NX))
    ENDIF
    IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
    !     DDLYRX(NN)=CumDepz2LayerBot_vr(L,NY,NX)-CDPTHX(L)
      DDLYRX(NN)=CDPTHY(L)-CDPTHX(L)
    ENDIF
!
!     RESET POND SURFACE LAYER NUMBER IF LOST TO EVAPORATION
!
  ELSEIF(NN.EQ.2)THEN
    IF((L.EQ.NU(NY,NX).AND.SoiBulkDensity_vr(NU(NY,NX),NY,NX).LE.ZERO) &
      .AND.(VHeatCapacity_vr(NU(NY,NX),NY,NX).LE.VHCPNX(NY,NX) &
      .OR.NUM(NY,NX).GT.NU(NY,NX)))THEN
      
      NUX=NU(NY,NX)
      DO LL=NUX+1,NL(NY,NX)
        IF(VLSoilPoreMicP_vr(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
          NU(NY,NX)         = LL
          DDLYRX(NN)        = DLYR(3,NUX,NY,NX)
          IFLGL(L,NN)       = 1
          DLYR(3,NUX,NY,NX) = 0.0_r8
          IF(SoiBulkDensity_vr(NUX,NY,NX).LE.ZERO)THEN
            VGeomLayer_vr(NUX,NY,NX)     = AREA(3,NUX,NY,NX)*DLYR(3,NUX,NY,NX)
            VLSoilPoreMicP_vr(NUX,NY,NX) = VGeomLayer_vr(NUX,NY,NX)*FracSoiAsMicP_vr(NUX,NY,NX)
          ENDIF
          exit
        ENDIF
      ENDDO
    ELSE
      DDLYRX(NN)=0.0_r8
      IFLGL(L,NN)=0
    ENDIF

!
!     RESET POND SURFACE LAYER NUMBER IF GAIN FROM PRECIPITATION
!
  ELSEIF(NN.EQ.3)THEN
    !CumSoilDeptht0 is the initial litter layer bottom
    !obtain the water exceeds litter layer water holding capacity
    XVOLWP=AZMAX1(VLWatMicP_vr(0,NY,NX)-VLWatheldCapSurf_col(NY,NX))
    IF(L.EQ.NU(NY,NX).AND.CumDepz2LayerBot_vr(0,NY,NX).GT.CumSoilDeptht0(NY,NX) &
      .AND.XVOLWP.GT.VLWatheldCapSurf_col(NY,NX)+VHCPNX(NY,NX)/cpw)THEN
          !     IF((SoiBulkDensity_vr(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))
          !    2.OR.(SoiBulkDensity_vr(L,NY,NX).LE.ZERO))THEN
      IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))THEN
        NU(NY,NX)               = NUI(NY,NX)
        NUM(NY,NX)              = NUI(NY,NX)
        DDLYRX(NN)              = (VLWatheldCapSurf_col(NY,NX)-XVOLWP)/AREA(3,0,NY,NX)
        IFLGL(L,NN)             = 1
        DLYR0                   = (AZMAX1(VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)-VWatLitRHoldCapcity_col(NY,NX))+VLitR_col(NY,NX))/AREA(3,0,NY,NX)
        DLYR(3,0,NY,NX)         = DLYR0+DDLYRX(NN)
        DLYR(3,NU(NY,NX),NY,NX) = DLYR(3,NU(NY,NX),NY,NX)-DDLYRX(NN)
        IF(L.GT.2)THEN
          DO LL=L-2,NU(NY,NX),-1
            CumDepz2LayerBot_vr(LL,NY,NX) = CumDepz2LayerBot_vr(L-1,NY,NX)
            CDPTHY(LL)                    = CDPTHY(L-1)
          ENDDO
        ENDIF
        CumDepz2LayerBot_vr(0,NY,NX)         = CumDepz2LayerBot_vr(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
        CDPTHY(0)                            = CDPTHY(NU(NY,NX))-DLYR(3,NU(NY,NX),NY,NX)
        SoiDepthMidLay_vr(NU(NY,NX),NY,NX)   = 0.5_r8*(CumDepz2LayerBot_vr(NU(NY,NX),NY,NX)+CumDepz2LayerBot_vr(0,NY,NX))
        CumSoilThickness_vr(NU(NY,NX),NY,NX) = DLYR(3,NU(NY,NX),NY,NX)
        DPTHZ_vr(NU(NY,NX),NY,NX)            = 0.5_r8*CumSoilThickness_vr(NU(NY,NX),NY,NX)
      ELSE
        DDLYRX(NN)  = 0.0_r8
        IFLGL(L,NN) = 0
      ENDIF
    ELSE
      DDLYRX(NN)  = 0.0_r8
      IFLGL(L,NN) = 0
    ENDIF
  ENDIF
  end subroutine UpdateLayerThickness
!------------------------------------------------------------------------------------------

  subroutine tgtPondLyr(L,L0,L1,NY,NX,NN,FX,FY,CDPTHY,IFLGL)
  implicit none
  integer, intent(in) :: L
  integer, intent(in) :: L0  !target
  integer, intent(in) :: L1  !source
  integer, intent(in) :: NY,NX,NN
  real(r8), intent(in) :: FX,FY
  real(r8), intent(inout) :: CDPTHY(0:JZ)
  integer, intent(in) ::  IFLGL(0:JZ,6)
  integer :: N,M,NZ,K,NGL,NR,NE,NTU,NTSA,NTSAB,NTG,NTP
  integer :: NTX,NTF,MID,idom
  real(r8) :: ENGY0,ENGY1
! begin_execution

  IF(L0.NE.0)THEN
    SAND(L1,NY,NX)                  = SAND(L1,NY,NX)+FX*SAND(L0,NY,NX)
    SILT(L1,NY,NX)                  = SILT(L1,NY,NX)+FX*SILT(L0,NY,NX)
    CLAY(L1,NY,NX)                  = CLAY(L1,NY,NX)+FX*CLAY(L0,NY,NX)
    trcx_solml_vr(idx_CEC,L1,NY,NX) = trcx_solml_vr(idx_CEC,L1,NY,NX)+FX*trcx_solml_vr(idx_CEC,L0,NY,NX)
    trcx_solml_vr(idx_AEC,L1,NY,NX) = trcx_solml_vr(idx_AEC,L1,NY,NX)+FX*trcx_solml_vr(idx_AEC,L0,NY,NX)
  ENDIF

  VLWatMicP_vr(L1,NY,NX)          = VLWatMicP_vr(L1,NY,NX)+FX*VLWatMicP_vr(L0,NY,NX)
  VLiceMicP_vr(L1,NY,NX)          = VLiceMicP_vr(L1,NY,NX)+FX*VLiceMicP_vr(L0,NY,NX)
  VLsoiAirP_vr(L1,NY,NX)          = VLsoiAirP_vr(L1,NY,NX)+FX*VLsoiAirP_vr(L0,NY,NX)
  VLMicP_vr(L1,NY,NX)             = VLMicP_vr(L1,NY,NX)+FX*VLMicP_vr(L0,NY,NX)
  VLSoilMicP_vr(L1,NY,NX)         = VLSoilMicP_vr(L1,NY,NX)+FX*VLSoilMicP_vr(L0,NY,NX)
  VLWatMicPX_vr(L1,NY,NX)         = VLWatMicP_vr(L1,NY,NX)
  ENGY1                           = VHeatCapacity_vr(L1,NY,NX)*TKS_vr(L1,NY,NX)
  ENGY0                           = VHeatCapacity_vr(L0,NY,NX)*TKS_vr(L0,NY,NX)
  ENGY1                           = ENGY1+FX*ENGY0
  VHeatCapacitySoilM_vr(L1,NY,NX) = VHeatCapacitySoilM_vr(L1,NY,NX)+FX*VHeatCapacitySoilM_vr(L0,NY,NX)
  VHeatCapacity_vr(L1,NY,NX)      = VHeatCapacitySoilM_vr(L1,NY,NX) &
    +cpw*(VLWatMicP_vr(L1,NY,NX)+VLWatMacP_vr(L1,NY,NX)) &
    +cpi*(VLiceMicP_vr(L1,NY,NX)+VLiceMacP_vr(L1,NY,NX))

  IF(VHeatCapacity_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L1,NY,NX)=ENGY1/VHeatCapacity_vr(L1,NY,NX)
  ELSE
    TKS_vr(L1,NY,NX)=TKS_vr(L0,NY,NX)
  ENDIF
  TCS(L1,NY,NX)=units%Kelvin2Celcius(TKS_vr(L1,NY,NX))

  DO NTF=ifertn_beg,ifertn_end
    FertN_soil_vr(NTF,L1,NY,NX)=FertN_soil_vr(NTF,L1,NY,NX)+FX*FertN_soil_vr(NTF,L0,NY,NX)
  ENDDO

  DO NTF=ifertnb_beg,ifertnb_end
    FertN_Band_vr(NTF,L1,NY,NX)=FertN_Band_vr(NTF,L1,NY,NX)+FX*FertN_Band_vr(NTF,L0,NY,NX)
  ENDDO

  DO NTU=ids_nuts_beg,ids_nuts_end
    if(NTU/=ids_H2PO4B .and. NTU/=ids_H1PO4B)THEN
      trc_solml_vr(NTU,L1,NY,NX)=trc_solml_vr(NTU,L1,NY,NX)+FX*trc_solml_vr(NTU,L0,NY,NX)
    ENDIF
  ENDDO

  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      trcSalt_solml_vr(NTSA,L1,NY,NX)=trcSalt_solml_vr(NTSA,L1,NY,NX)&
        +FX*trcSalt_solml_vr(NTSA,L0,NY,NX)
    ENDDO
  ENDIF

  IF(L0.NE.0)THEN
    trc_solml_vr(ids_H1PO4B,L1,NY,NX)=trc_solml_vr(ids_H1PO4B,L1,NY,NX) &
      +FX*trc_solml_vr(ids_H1PO4B,L0,NY,NX)
    trc_solml_vr(ids_H2PO4B,L1,NY,NX)=trc_solml_vr(ids_H2PO4B,L1,NY,NX) &
      +FX*trc_solml_vr(ids_H2PO4B,L0,NY,NX)

    IF(salt_model)THEN
      DO NTSAB=idsaltb_beg,idsaltb_end
        trcSalt_solml_vr(NTSAB,L1,NY,NX)=trcSalt_solml_vr(NTSAB,L1,NY,NX)+FX*trcSalt_solml_vr(NTSAB,L0,NY,NX)
      ENDDO
    ENDIF

    DO NTX=idx_beg+1,idx_cation_end
      trcx_solml_vr(NTX,L1,NY,NX)=trcx_solml_vr(NTX,L1,NY,NX)+FX*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO

    DO NTX=idx_AEC+1,idx_end
      trcx_solml_vr(NTX,L1,NY,NX)=trcx_solml_vr(NTX,L1,NY,NX)+FX*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO

    DO NTP=idsp_beg,idsp_end
      trcp_saltpml_vr(NTP,L1,NY,NX)=trcp_saltpml_vr(NTP,L1,NY,NX)+FX*trcp_saltpml_vr(NTP,L0,NY,NX)
    ENDDO

    DO NTG=idg_beg,idg_end-1
      trc_gasml_vr(NTG,L1,NY,NX)=trc_gasml_vr(NTG,L1,NY,NX)+FX*trc_gasml_vr(NTG,L0,NY,NX)
    ENDDO
  ENDIF
!exclude NH3 and NH3B, which are accounted in nutrients
  DO NTG=idg_beg,idg_end-2
    trc_solml_vr(NTG,L1,NY,NX)=trc_solml_vr(NTG,L1,NY,NX)+FX*trc_solml_vr(NTG,L0,NY,NX)
  ENDDO

  IF(IFLGL(L,3).EQ.0)THEN
    DO  K=1,jcplx
      DO  N=1,NumMicbFunGrupsPerCmplx
        DO  M=1,nlbiomcp
          DO NGL=JGnio(N),JGnfo(N)
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              mBiomeHeter_vr(NE,MID,K,L1,NY,NX)=mBiomeHeter_vr(NE,MID,K,L1,NY,NX)+FX*mBiomeHeter_vr(NE,MID,K,L0,NY,NX)
            ENDDO
          enddo
        enddo
      enddo
    ENDDO
    DO  N=1,NumMicbFunGrupsPerCmplx
      DO  M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor_vr(NE,MID,L1,NY,NX)=mBiomeAutor_vr(NE,MID,L1,NY,NX)+FX*mBiomeAutor_vr(NE,MID,L0,NY,NX)
          ENDDO
        enddo
      enddo
    enddo

    DO K=1,jcplx
      DO  M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          OMBioResdu_vr(NE,M,K,L1,NY,NX)=OMBioResdu_vr(NE,M,K,L1,NY,NX)+FX*OMBioResdu_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
      DO idom=idom_beg,idom_end
        DOM_vr(idom,K,L1,NY,NX)      = DOM_vr(idom,K,L1,NY,NX)+FX*DOM_vr(idom,K,L0,NY,NX)
        DOM_MacP_vr(idom,K,L1,NY,NX) = DOM_MacP_vr(idom,K,L1,NY,NX)+FX*DOM_MacP_vr(idom,K,L0,NY,NX)
        SorbedOM_vr(idom,K,L1,NY,NX) = SorbedOM_vr(idom,K,L1,NY,NX)+FX*SorbedOM_vr(idom,K,L0,NY,NX)
      enddo
      
      DO M=1,jsken
        SolidOMAct_vr(M,K,L1,NY,NX)=SolidOMAct_vr(M,K,L1,NY,NX)+FX*SolidOMAct_vr(M,K,L0,NY,NX)
        DO NE=1,NumPlantChemElms
          SolidOM_vr(NE,M,K,L1,NY,NX)=SolidOM_vr(NE,M,K,L1,NY,NX)+FX*SolidOM_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
          !
!     TARGET ROOT LAYER
!
  IF(L0.NE.0)THEN
    DO  NZ=1,NP(NY,NX)
      IF(RootMycoActiveBiomC_pvr(ipltroot,L0,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX) &
        .AND.RootMycoActiveBiomC_pvr(ipltroot,L1,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX))THEN
        DO N=1,MY(NZ,NY,NX)
          DO NTG=idg_beg,idg_end-1
            trcg_rootml_pvr(NTG,N,L1,NZ,NY,NX) = trcg_rootml_pvr(NTG,N,L1,NZ,NY,NX)+FX*trcg_rootml_pvr(NTG,N,L0,NZ,NY,NX)
            trcs_rootml_pvr(NTG,N,L1,NZ,NY,NX) = trcs_rootml_pvr(NTG,N,L1,NZ,NY,NX)+FX*trcs_rootml_pvr(NTG,N,L0,NZ,NY,NX)
          ENDDO
          DO  NR=1,NumRootAxes_pft(NZ,NY,NX)
            DO NE=1,NumPlantChemElms
              RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX)=RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX) &
                +FX*RootMyco1stStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
              RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX)=RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX) &
                +FX*RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
            ENDDO
            Root1stLen_rpvr(N,L1,NR,NZ,NY,NX)  = Root1stLen_rpvr(N,L1,NR,NZ,NY,NX)+FX*Root1stLen_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX)   = Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX)+FX*Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX) = Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX)+FX*Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)
          ENDDO
          DO NE=1,NumPlantChemElms
             RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX)=RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX)+FX* RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)
          ENDDO
          RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX) = RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX)+FX*RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)
          PopuRootMycoC_pvr(N,L1,NZ,NY,NX)       = PopuRootMycoC_pvr(N,L1,NZ,NY,NX)+FX* PopuRootMycoC_pvr(N,L0,NZ,NY,NX)
          RootProteinC_pvr(N,L1,NZ,NY,NX)        = RootProteinC_pvr(N,L1,NZ,NY,NX)+FX*RootProteinC_pvr(N,L0,NZ,NY,NX)
          Root1stXNumL_pvr(N,L1,NZ,NY,NX)        = Root1stXNumL_pvr(N,L1,NZ,NY,NX)+FX*Root1stXNumL_pvr(N,L0,NZ,NY,NX)
          Root2ndXNum_pvr(N,L1,NZ,NY,NX)         = Root2ndXNum_pvr(N,L1,NZ,NY,NX)+FX*Root2ndXNum_pvr(N,L0,NZ,NY,NX)
          RootLenPerPlant_pvr(N,L1,NZ,NY,NX)     = RootLenPerPlant_pvr(N,L1,NZ,NY,NX)+FX*RootLenPerPlant_pvr(N,L0,NZ,NY,NX)
          RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX) = RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX)+FX*RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)
          RootPoreVol_pvr(N,L1,NZ,NY,NX)         = RootPoreVol_pvr(N,L1,NZ,NY,NX)+FX*RootPoreVol_pvr(N,L0,NZ,NY,NX)
          RootVH2O_pvr(N,L1,NZ,NY,NX)            = RootVH2O_pvr(N,L1,NZ,NY,NX)+FX*RootVH2O_pvr(N,L0,NZ,NY,NX)
          Root1stRadius_pvr(N,L1,NZ,NY,NX)       = Root1stRadius_pvr(N,L1,NZ,NY,NX)+FX*Root1stRadius_pvr(N,L0,NZ,NY,NX)
          Root2ndRadius_pvr(N,L1,NZ,NY,NX)       = Root2ndRadius_pvr(N,L1,NZ,NY,NX)+FX*Root2ndRadius_pvr(N,L0,NZ,NY,NX)
          RootAreaPerPlant_pvr(N,L1,NZ,NY,NX)    = RootAreaPerPlant_pvr(N,L1,NZ,NY,NX)+FX*RootAreaPerPlant_pvr(N,L0,NZ,NY,NX)
          Root2ndAveLen_pvr(N,L1,NZ,NY,NX)       = Root2ndAveLen_pvr(N,L1,NZ,NY,NX)+FX*Root2ndAveLen_pvr(N,L0,NZ,NY,NX)
        ENDDO
        DO NE=1,NumPlantChemElms
          RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX)+FX*RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)
          RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX)+FX*RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!
          !     SOURCE POND LAYER
    !
  IF(L0.NE.0)THEN
    SAND(L0,NY,NX)                  = FY*SAND(L0,NY,NX)
    SILT(L0,NY,NX)                  = FY*SILT(L0,NY,NX)
    CLAY(L0,NY,NX)                  = FY*CLAY(L0,NY,NX)
    trcx_solml_vr(idx_CEC,L0,NY,NX) = FY*trcx_solml_vr(idx_CEC,L0,NY,NX)
    trcx_solml_vr(idx_AEC,L0,NY,NX) = FY*trcx_solml_vr(idx_AEC,L0,NY,NX)
  ENDIF
!     IF(SoiBulkDensity_vr(L0,NY,NX).LE.ZERO)THEN
!     VGeomLayer_vr(L0,NY,NX)=FY*VGeomLayer_vr(L0,NY,NX)
!     VLSoilPoreMicP_vr(L0,NY,NX)=FY*VLSoilPoreMicP_vr(L0,NY,NX)
!     ENDIF
  VLWatMicP_vr(L0,NY,NX)          = FY*VLWatMicP_vr(L0,NY,NX)
  VLiceMicP_vr(L0,NY,NX)          = FY*VLiceMicP_vr(L0,NY,NX)
  VLsoiAirP_vr(L0,NY,NX)          = FY*VLsoiAirP_vr(L0,NY,NX)
  VLMicP_vr(L0,NY,NX)             = FY*VLMicP_vr(L0,NY,NX)
  VLSoilMicP_vr(L0,NY,NX)         = FY*VLSoilMicP_vr(L0,NY,NX)
  VLWatMicPX_vr(L0,NY,NX)         = VLWatMicP_vr(L0,NY,NX)
  ENGY0                           = FY*ENGY0
  VHeatCapacitySoilM_vr(L0,NY,NX) = FY*VHeatCapacitySoilM_vr(L0,NY,NX)
  IF(L0.NE.0)THEN
    VHeatCapacity_vr(L0,NY,NX)=VHeatCapacitySoilM_vr(L0,NY,NX) &
      +cpw*(VLWatMicP_vr(L0,NY,NX)+VLWatMacP_vr(L0,NY,NX)) &
      +cpi*(VLiceMicP_vr(L0,NY,NX)+VLiceMacP_vr(L0,NY,NX))
  ELSE
    VHeatCapacity_vr(L0,NY,NX)=VHeatCapacitySoilM_vr(L0,NY,NX)+cpw*VLWatMicP_vr(L0,NY,NX)+cpi*VLiceMicP_vr(L0,NY,NX)
  ENDIF
  IF(VHeatCapacity_vr(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L0,NY,NX)=ENGY0/VHeatCapacity_vr(L0,NY,NX)
  ELSE
    TKS_vr(L0,NY,NX)=TKS_vr(L1,NY,NX)
  ENDIF
  TCS(L0,NY,NX)=units%Kelvin2Celcius(TKS_vr(L0,NY,NX))

  DO NTF=ifertn_beg,ifertn_end
    FertN_soil_vr(NTF,L0,NY,NX)=FY*FertN_soil_vr(NTF,L0,NY,NX)
  ENDDO

  DO NTF=ifertnb_beg,ifertnb_end
    FertN_Band_vr(NTF,L0,NY,NX)=FY*FertN_Band_vr(NTF,L0,NY,NX)
  ENDDO

  DO NTU=ids_nuts_beg,ids_nuts_end
    if(NTU/=ids_H1PO4B .and. NTU/=ids_H2PO4B)THEN
      trc_solml_vr(NTU,L0,NY,NX)=FY*trc_solml_vr(NTU,L0,NY,NX)
    ENDIF
  ENDDO
  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      trcSalt_solml_vr(NTSA,L0,NY,NX)=FY*trcSalt_solml_vr(NTSA,L0,NY,NX)
    ENDDO
  ENDIF
  IF(L0.NE.0)THEN
    trc_solml_vr(ids_H1PO4B,L0,NY,NX) = FY*trc_solml_vr(ids_H1PO4B,L0,NY,NX)
    trc_solml_vr(ids_H2PO4B,L0,NY,NX) = FY*trc_solml_vr(ids_H2PO4B,L0,NY,NX)
    IF(salt_model)THEN
      DO NTSAB=idsaltb_beg,idsaltb_end
        trcSalt_solml_vr(NTSAB,L0,NY,NX)=FY*trcSalt_solml_vr(NTSAB,L0,NY,NX)
      ENDDO
    ENDIF

    DO NTX=idx_beg+1,idx_cation_end
      trcx_solml_vr(NTX,L0,NY,NX)=FY*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO
    DO NTX=idx_AEC+1,idx_end
      trcx_solml_vr(NTX,L0,NY,NX)=FY*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO

    DO NTP=idsp_beg,idsp_end
      trcp_saltpml_vr(NTP,L0,NY,NX)=FY*trcp_saltpml_vr(NTP,L0,NY,NX)
    ENDDO

    DO NTG=idg_beg,idg_end-1
      trc_gasml_vr(NTG,L0,NY,NX)=FY*trc_gasml_vr(NTG,L0,NY,NX)
    ENDDO
  ENDIF
!exclude NH3 and NH3B, which are accounted in nutrients
  DO NTG=idg_beg,idg_end-2
    trc_solml_vr(NTG,L0,NY,NX)=FY*trc_solml_vr(NTG,L0,NY,NX)
  ENDDO
  IF(IFLGL(L,3).EQ.0)THEN
    DO  K=1,jcplx
       DO N=1,NumMicbFunGrupsPerCmplx
        DO M=1,nlbiomcp
          DO NGL=JGnio(N),JGnfo(N)
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              mBiomeHeter_vr(NE,MID,K,L0,NY,NX)=FY*mBiomeHeter_vr(NE,MID,K,L0,NY,NX)
            ENDDO
          ENDDO
        enddo
      enddo
    ENDDO

    DO N=1,NumMicbFunGrupsPerCmplx
      DO M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor_vr(NE,MID,L0,NY,NX)=FY*mBiomeAutor_vr(NE,MID,L0,NY,NX)
          ENDDO
        ENDDO
      enddo
    enddo

    DO K=1,jcplx
      DO  M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          OMBioResdu_vr(NE,M,K,L0,NY,NX)=FY*OMBioResdu_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
      do idom=idom_beg,idom_end
        DOM_vr(idom,K,L0,NY,NX)      = FY*DOM_vr(idom,K,L0,NY,NX)
        DOM_MacP_vr(idom,K,L0,NY,NX) = FY*DOM_MacP_vr(idom,K,L0,NY,NX)
        SorbedOM_vr(idom,K,L0,NY,NX) = FY*SorbedOM_vr(idom,K,L0,NY,NX)
      enddo
      DO  M=1,jsken
        SolidOMAct_vr(M,K,L0,NY,NX)=FY*SolidOMAct_vr(M,K,L0,NY,NX)
        do NE=1,NumPlantChemElms
          SolidOM_vr(NE,M,K,L0,NY,NX)=FY*SolidOM_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
          !
          !     SOURCE ROOT LAYER
          !
  IF(L0.NE.0)THEN
    DO  NZ=1,NP(NY,NX)
      IF(RootMycoActiveBiomC_pvr(ipltroot,L0,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX) &
        .AND.RootMycoActiveBiomC_pvr(ipltroot,L1,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX))THEN
        DO  N=1,MY(NZ,NY,NX)
          DO NTG=idg_beg,idg_end-1
            trcg_rootml_pvr(NTG,N,L0,NZ,NY,NX) = FY*trcg_rootml_pvr(NTG,N,L0,NZ,NY,NX)
            trcs_rootml_pvr(NTG,N,L0,NZ,NY,NX) = FY*trcs_rootml_pvr(NTG,N,L0,NZ,NY,NX)
          ENDDO
          DO NR=1,NumRootAxes_pft(NZ,NY,NX)
            DO NE=1,NumPlantChemElms
              RootMyco1stStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX) = FY*RootMyco1stStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
              RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX) = FY*RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
            ENDDO
            Root1stLen_rpvr(N,L0,NR,NZ,NY,NX)  = FY*Root1stLen_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)   = FY*Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX) = FY*Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)
          ENDDO
          DO NE=1,NumPlantChemElms
            RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)=FY* RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)
          ENDDO
          RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX) = FY*RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)
          PopuRootMycoC_pvr(N,L0,NZ,NY,NX)       = FY* PopuRootMycoC_pvr(N,L0,NZ,NY,NX)
          RootProteinC_pvr(N,L0,NZ,NY,NX)        = FY*RootProteinC_pvr(N,L0,NZ,NY,NX)
          Root1stXNumL_pvr(N,L0,NZ,NY,NX)        = FY*Root1stXNumL_pvr(N,L0,NZ,NY,NX)
          Root2ndXNum_pvr(N,L0,NZ,NY,NX)         = FY*Root2ndXNum_pvr(N,L0,NZ,NY,NX)
          RootLenPerPlant_pvr(N,L0,NZ,NY,NX)     = FY*RootLenPerPlant_pvr(N,L0,NZ,NY,NX)
          RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX) = FY*RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)
          RootPoreVol_pvr(N,L0,NZ,NY,NX)         = FY*RootPoreVol_pvr(N,L0,NZ,NY,NX)
          RootVH2O_pvr(N,L0,NZ,NY,NX)            = FY*RootVH2O_pvr(N,L0,NZ,NY,NX)
          Root1stRadius_pvr(N,L0,NZ,NY,NX)       = FY*Root1stRadius_pvr(N,L0,NZ,NY,NX)
          Root2ndRadius_pvr(N,L0,NZ,NY,NX)       = FY*Root2ndRadius_pvr(N,L0,NZ,NY,NX)
          RootAreaPerPlant_pvr(N,L0,NZ,NY,NX)    = FY*RootAreaPerPlant_pvr(N,L0,NZ,NY,NX)
          Root2ndAveLen_pvr(N,L0,NZ,NY,NX)       = FY*Root2ndAveLen_pvr(N,L0,NZ,NY,NX)
        ENDDO
        DO NE=1,NumPlantChemElms
          RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)=FY*RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)

          RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)=FY*RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  IF(NN.EQ.1)THEN
    IF(SoiBulkDensity_vr(L0,NY,NX).LE.ZERO.AND.SoiBulkDensity_vr(L1,NY,NX).LE.ZERO &
      .AND.VLWatMicP_vr(L0,NY,NX)+VLiceMicP_vr(L0,NY,NX).LE.ZEROS(NY,NX))THEN
      CumDepz2LayerBot_vr(L1,NY,NX)=CumDepz2LayerBot_vr(L0,NY,NX)
      CDPTHY(L1)=CDPTHY(L0)
    ENDIF
  ENDIF
  end subroutine tgtPondLyr

!------------------------------------------------------------------------------------------

  subroutine MoveSOM(L0,L1,L,NY,NX,FO,IFLGL)
!
  implicit none
  integer, intent(in) :: L0,L1,L
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: FO
  integer, intent(in) :: IFLGL(0:JZ,6)

  integer :: K,N,M,NGL,NR,NZ,NE,NTG,MID,idom
  real(r8) :: FXO,FRO
  real(r8) :: FXRTLG2,FXRTN2,FXEPOOLR,FXWTRTL
  real(r8) :: WTNDLE,FXEPOOLN
  real(r8) :: FXWTRT1E
  real(r8) :: FXWTRT2E,FXRTLG1
  real(r8) :: FXWTRTD,FXWSRTL,FRootPrimeAxsNum,FXRTNL,FXRTLGP,FXRTDNP
  real(r8) :: FXRTVLP,FXRTVLW,FXRRAD1,FXRRAD2,FXRootAreaPerPlant_pvr,FXRTLGA
  real(r8) :: FXOQN,FXOQP,FXOQA,FXOQCH,FXOQNH,FXOQPH,FXOQAH
  real(r8) :: FXOHC,FXOHN,FXOHP,FXOHA,FXOSC,FXOSA,FXOSN,FXOSP
  real(r8) :: FXOMC,FXOMN,FXOMP,FXORC,FXORN,FXORP,FXOQC
  real(r8) :: FXGA,FXGP
  
! begin_execution
  IF(IFLGL(L,3).EQ.0 .AND. L0.NE.0 &
    .AND. VLSoilPoreMicP_vr(L0,NY,NX).GT.ZEROS(NY,NX) &
    .AND. VLSoilPoreMicP_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    IF(L0.EQ.L.OR.CORGCI(L0,NY,NX).LE.ZERO)THEN
      FXO=FO
    ELSE
      FXO=AMIN1(0.5,FO*AMIN1(10.0,CORGCI(L1,NY,NX)/CORGCI(L0,NY,NX)))
    ENDIF

    DO  K=1,jcplx
      DO  N=1,NumMicbFunGrupsPerCmplx
        DO  M=1,nlbiomcp
          DO NGL=JGnio(N),JGnfo(N)
            MID=micpar%get_micb_id(M,NGL)          
            DO NE=1,NumPlantChemElms
              FXOMC                             = FXO*mBiomeHeter_vr(NE,MID,K,L0,NY,NX)
              mBiomeHeter_vr(NE,MID,K,L1,NY,NX) = mBiomeHeter_vr(NE,MID,K,L1,NY,NX)+FXOMC
              mBiomeHeter_vr(NE,MID,K,L0,NY,NX) = mBiomeHeter_vr(NE,MID,K,L0,NY,NX)-FXOMC
            ENDDO
          enddo
        enddo
      enddo
    ENDDO

    DO  N=1,NumMicbFunGrupsPerCmplx
      DO  M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            FXOMC                           = FXO*mBiomeAutor_vr(NE,MID,L0,NY,NX)
            mBiomeAutor_vr(NE,MID,L1,NY,NX) = mBiomeAutor_vr(NE,MID,L1,NY,NX)+FXOMC
            mBiomeAutor_vr(NE,MID,L0,NY,NX) = mBiomeAutor_vr(NE,MID,L0,NY,NX)-FXOMC
          ENDDO
        enddo
      enddo
    enddo

    DO  K=1,jcplx
      DO  M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          FXORC                          = FXO*OMBioResdu_vr(NE,M,K,L0,NY,NX)
          OMBioResdu_vr(NE,M,K,L1,NY,NX) = OMBioResdu_vr(NE,M,K,L1,NY,NX)+FXORC
          OMBioResdu_vr(NE,M,K,L0,NY,NX) = OMBioResdu_vr(NE,M,K,L0,NY,NX)-FXORC
        ENDDO
      ENDDO

      DO idom=idom_beg,idom_end
        FXOQC                   = FXO*DOM_vr(idom,K,L0,NY,NX)
        DOM_vr(idom,K,L1,NY,NX) = DOM_vr(idom,K,L1,NY,NX)+FXOQC
        DOM_vr(idom,K,L0,NY,NX) = DOM_vr(idom,K,L0,NY,NX)-FXOQC
      enddo

      IF(SoilFracAsMacP_vr(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP_vr(L0,NY,NX).GT.ZERO)THEN
        DO idom=idom_beg,idom_end
          FXOQCH                       = FXO*DOM_MacP_vr(idom,K,L0,NY,NX)
          DOM_MacP_vr(idom,K,L1,NY,NX) = DOM_MacP_vr(idom,K,L1,NY,NX)+FXOQCH
          DOM_MacP_vr(idom,K,L0,NY,NX) = DOM_MacP_vr(idom,K,L0,NY,NX)-FXOQCH
        ENDDO
      ENDIF
       DO idom=idom_beg,idom_end
        FXOHC                        = FXO*SorbedOM_vr(idom,K,L0,NY,NX)
        SorbedOM_vr(idom,K,L1,NY,NX) = SorbedOM_vr(idom,K,L1,NY,NX)+FXOHC
        SorbedOM_vr(idom,K,L0,NY,NX) = SorbedOM_vr(idom,K,L0,NY,NX)-FXOHC
      ENDDO

      DO M=1,jsken
        DO NE=1,NumPlantChemElms
          FXOSC                       = FXO*SolidOM_vr(NE,M,K,L0,NY,NX)
          SolidOM_vr(NE,M,K,L1,NY,NX) = SolidOM_vr(NE,M,K,L1,NY,NX)+FXOSC
          SolidOM_vr(NE,M,K,L0,NY,NX) = SolidOM_vr(NE,M,K,L0,NY,NX)-FXOSC
        ENDDO
      ENDDO
    ENDDO
!
!     ROOTS (why return?) 05/21/2022, jyt
!
    return
    DO NZ=1,NP(NY,NX)
      IF(RootMycoActiveBiomC_pvr(ipltroot,L0,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX) &
        .AND.RootMycoActiveBiomC_pvr(ipltroot,L1,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX))THEN
        IF(L0.EQ.L.OR.DPTHZ_vr(L1,NY,NX).LE.ZERO)THEN
          FRO=FO
        ELSE
          FRO=AMIN1(0.5,FO*DPTHZ_vr(L0,NY,NX)/DPTHZ_vr(L1,NY,NX))
        ENDIF

        DO  N=1,MY(NZ,NY,NX)
          DO NTG=idg_beg,idg_end-1
            FXGA                               = FRO*trcg_rootml_pvr(NTG,N,L0,NZ,NY,NX)
            trcg_rootml_pvr(NTG,N,L1,NZ,NY,NX) = trcg_rootml_pvr(NTG,N,L1,NZ,NY,NX)+FXGA
            trcg_rootml_pvr(NTG,N,L0,NZ,NY,NX) = trcg_rootml_pvr(NTG,N,L0,NZ,NY,NX)-FXGA

            FXGP                               = FRO*trcs_rootml_pvr(NTG,N,L0,NZ,NY,NX)
            trcs_rootml_pvr(NTG,N,L1,NZ,NY,NX) = trcs_rootml_pvr(NTG,N,L1,NZ,NY,NX)+FXGP
            trcs_rootml_pvr(NTG,N,L0,NZ,NY,NX) = trcs_rootml_pvr(NTG,N,L0,NZ,NY,NX)-FXGP
          ENDDO

          DO  NR=1,NumRootAxes_pft(NZ,NY,NX)
            DO NE=1,NumPlantChemElms
              FXWTRT1E                                       = FRO*RootMyco1stStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
              RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX) = RootMyco1stStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX)+FXWTRT1E
              RootMyco1stStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX) = RootMyco1stStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)-FXWTRT1E

              FXWTRT2E                                       = FRO*RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
              RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX) = RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX)+FXWTRT2E
              RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX) = RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)-FXWTRT2E
            ENDDO
            FXRTLG1                           = FRO*Root1stLen_rpvr(N,L0,NR,NZ,NY,NX)
            Root1stLen_rpvr(N,L1,NR,NZ,NY,NX) = Root1stLen_rpvr(N,L1,NR,NZ,NY,NX)+FXRTLG1
            Root1stLen_rpvr(N,L0,NR,NZ,NY,NX) = Root1stLen_rpvr(N,L0,NR,NZ,NY,NX)-FXRTLG1

            FXRTLG2                          = FRO*Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX) = Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX)+FXRTLG2
            Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX) = Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)-FXRTLG2

            FXRTN2                             = FRO*Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX) = Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX)+FXRTN2
            Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX) = Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)-FXRTN2
          ENDDO
          DO NE=1,NumPlantChemElms
             FXEPOOLR                                 = FRO*RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)
             RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX) = RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX)+FXEPOOLR
             RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX) = RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)-FXEPOOLR
          ENDDO

          FXWTRTL                                = FRO*RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)
          RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX) = RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX)+FXWTRTL
          RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX) = RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)-FXWTRTL

          FXWTRTD                          = FRO*PopuRootMycoC_pvr(N,L0,NZ,NY,NX)
          PopuRootMycoC_pvr(N,L1,NZ,NY,NX) = PopuRootMycoC_pvr(N,L1,NZ,NY,NX)+FXWTRTD
          PopuRootMycoC_pvr(N,L0,NZ,NY,NX) = PopuRootMycoC_pvr(N,L0,NZ,NY,NX)-FXWTRTD

          FXWSRTL                         = FRO*RootProteinC_pvr(N,L0,NZ,NY,NX)
          RootProteinC_pvr(N,L1,NZ,NY,NX) = RootProteinC_pvr(N,L1,NZ,NY,NX)+FXWSRTL
          RootProteinC_pvr(N,L0,NZ,NY,NX) = RootProteinC_pvr(N,L0,NZ,NY,NX)-FXWSRTL

          FRootPrimeAxsNum                = FRO*Root1stXNumL_pvr(N,L0,NZ,NY,NX)
          Root1stXNumL_pvr(N,L1,NZ,NY,NX) = Root1stXNumL_pvr(N,L1,NZ,NY,NX)+FRootPrimeAxsNum
          Root1stXNumL_pvr(N,L0,NZ,NY,NX) = Root1stXNumL_pvr(N,L0,NZ,NY,NX)-FRootPrimeAxsNum

          FXRTNL                         = FRO*Root2ndXNum_pvr(N,L0,NZ,NY,NX)
          Root2ndXNum_pvr(N,L1,NZ,NY,NX) = Root2ndXNum_pvr(N,L1,NZ,NY,NX)+FXRTNL
          Root2ndXNum_pvr(N,L0,NZ,NY,NX) = Root2ndXNum_pvr(N,L0,NZ,NY,NX)-FXRTNL

          FXRTLGP                            = FRO*RootLenPerPlant_pvr(N,L0,NZ,NY,NX)
          RootLenPerPlant_pvr(N,L1,NZ,NY,NX) = RootLenPerPlant_pvr(N,L1,NZ,NY,NX)+FXRTLGP
          RootLenPerPlant_pvr(N,L0,NZ,NY,NX) = RootLenPerPlant_pvr(N,L0,NZ,NY,NX)-FXRTLGP

          FXRTDNP                                = FRO*RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)
          RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX) = RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX)+FXRTDNP
          RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX) = RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)-FXRTDNP

          FXRTVLP                        = FRO*RootPoreVol_pvr(N,L0,NZ,NY,NX)
          RootPoreVol_pvr(N,L1,NZ,NY,NX) = RootPoreVol_pvr(N,L1,NZ,NY,NX)+FXRTVLP
          RootPoreVol_pvr(N,L0,NZ,NY,NX) = RootPoreVol_pvr(N,L0,NZ,NY,NX)-FXRTVLP

          FXRTVLW                     = FRO*RootVH2O_pvr(N,L0,NZ,NY,NX)
          RootVH2O_pvr(N,L1,NZ,NY,NX) = RootVH2O_pvr(N,L1,NZ,NY,NX)+FXRTVLW
          RootVH2O_pvr(N,L0,NZ,NY,NX) = RootVH2O_pvr(N,L0,NZ,NY,NX)-FXRTVLW

          FXRRAD1                          = FRO*Root1stRadius_pvr(N,L0,NZ,NY,NX)
          Root1stRadius_pvr(N,L1,NZ,NY,NX) = Root1stRadius_pvr(N,L1,NZ,NY,NX)+FXRRAD1
          Root1stRadius_pvr(N,L0,NZ,NY,NX) = Root1stRadius_pvr(N,L0,NZ,NY,NX)-FXRRAD1

          FXRRAD2                          = FRO*Root2ndRadius_pvr(N,L0,NZ,NY,NX)
          Root2ndRadius_pvr(N,L1,NZ,NY,NX) = Root2ndRadius_pvr(N,L1,NZ,NY,NX)+FXRRAD2
          Root2ndRadius_pvr(N,L0,NZ,NY,NX) = Root2ndRadius_pvr(N,L0,NZ,NY,NX)-FXRRAD2

          FXRootAreaPerPlant_pvr              = FRO*RootAreaPerPlant_pvr(N,L0,NZ,NY,NX)
          RootAreaPerPlant_pvr(N,L1,NZ,NY,NX) = RootAreaPerPlant_pvr(N,L1,NZ,NY,NX)+FXRootAreaPerPlant_pvr
          RootAreaPerPlant_pvr(N,L0,NZ,NY,NX) = RootAreaPerPlant_pvr(N,L0,NZ,NY,NX)-FXRootAreaPerPlant_pvr

          FXRTLGA                          = FRO*Root2ndAveLen_pvr(N,L0,NZ,NY,NX)
          Root2ndAveLen_pvr(N,L1,NZ,NY,NX) = Root2ndAveLen_pvr(N,L1,NZ,NY,NX)+FXRTLGA
          Root2ndAveLen_pvr(N,L0,NZ,NY,NX) = Root2ndAveLen_pvr(N,L0,NZ,NY,NX)-FXRTLGA
        ENDDO
!
!     ROOT NODULES
!
        DO NE=1,NumPlantChemElms
          WTNDLE                                 = FRO*RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)
          RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX)+WTNDLE
          RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX) = RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)-WTNDLE

          FXEPOOLN                               = FRO*RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)
          RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX)+FXEPOOLN
          RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX) = RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)-FXEPOOLN

        ENDDO
      ENDIF
    ENDDO
  ENDIF
  end subroutine MoveSOM
!------------------------------------------------------------------------------------------

  subroutine MoveMacPoreSolute(L0,L1,NY,NX,FHO)

  implicit none
  integer, intent(in) :: L0,L1,NY,NX
  real(r8), intent(in) :: FHO
  integer :: NTS,NTSA
  real(r8) :: FXSH,FXH

! begin_execution
!
!     SOIL MACROPORE N,P SOLUTES
!
  IF(SoilFracAsMacP_vr(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP_vr(L0,NY,NX).GT.ZERO)THEN

    DO NTS=ids_beg,ids_end
      FXSH                       = FHO*trc_soHml_vr(NTS,L0,NY,NX)
      trc_soHml_vr(NTS,L1,NY,NX) = trc_soHml_vr(NTS,L1,NY,NX)+FXSH
      trc_soHml_vr(NTS,L0,NY,NX) = trc_soHml_vr(NTS,L0,NY,NX)-FXSH
    ENDDO

!
!     SOIL MACROPORE SOLUBLE SALTS
!
    IF(salt_model)THEN
      DO NTSA=idsalt_beg,idsaltb_end
        FXH                             = FHO*trcSalt_soHml_vr(NTSA,L0,NY,NX)
        trcSalt_soHml_vr(NTSA,L1,NY,NX) = trcSalt_soHml_vr(NTSA,L1,NY,NX)+FXH
        trcSalt_soHml_vr(NTSA,L0,NY,NX) = trcSalt_soHml_vr(NTSA,L0,NY,NX)-FXH
      ENDDO
    ENDIF
!
  ENDIF
  end subroutine MoveMacPoreSolute

!------------------------------------------------------------------------------------------

  subroutine MoveBandSolute(L,L0,L1,NY,NX,FX,FWO,FO)
  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in) :: FX,FWO,FO
  real(r8) :: FAO,FCO
  integer  :: NTSAB,NTP,NTX
  real(r8) :: FXB,FXP,FXXC,FXXA
  real(r8) :: FXH1POB,FXH2POB

! only phosphrous is considered below because there is no gaseous phase.
  FXH1POB                           = FWO*trc_solml_vr(ids_H1PO4B,L0,NY,NX)
  trc_solml_vr(ids_H1PO4B,L1,NY,NX) = trc_solml_vr(ids_H1PO4B,L1,NY,NX)+FXH1POB
  trc_solml_vr(ids_H1PO4B,L0,NY,NX) = trc_solml_vr(ids_H1PO4B,L0,NY,NX)-FXH1POB

  FXH2POB                           = FWO*trc_solml_vr(ids_H2PO4B,L0,NY,NX)
  trc_solml_vr(ids_H2PO4B,L1,NY,NX) = trc_solml_vr(ids_H2PO4B,L1,NY,NX)+FXH2POB
  trc_solml_vr(ids_H2PO4B,L0,NY,NX) = trc_solml_vr(ids_H2PO4B,L0,NY,NX)-FXH2POB

  IF(salt_model)THEN
    DO NTSAB=idsaltb_beg,idsaltb_end
      FXB                              = FWO*trcSalt_solml_vr(NTSAB,L0,NY,NX)
      trcSalt_solml_vr(NTSAB,L1,NY,NX) = trcSalt_solml_vr(NTSAB,L1,NY,NX)+FXB
      trcSalt_solml_vr(NTSAB,L0,NY,NX) = trcSalt_solml_vr(NTSAB,L0,NY,NX)-FXB
    ENDDO
  ENDIF
!
!     SOIL ADSORBED CATIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L .OR. CEC_vr(L0,NY,NX).LE.ZERO)THEN
    FCO=FO
  ELSE
    FCO=AMIN1(0.5_r8,FO*CEC_vr(L1,NY,NX)/CEC_vr(L0,NY,NX))
  ENDIF
  DO NTX=idx_beg,idx_cation_end
    FXXC                        = FCO*trcx_solml_vr(NTX,L0,NY,NX)
    trcx_solml_vr(NTX,L1,NY,NX) = trcx_solml_vr(NTX,L1,NY,NX)+FXXC
    trcx_solml_vr(NTX,L0,NY,NX) = trcx_solml_vr(NTX,L0,NY,NX)-FXXC
  ENDDO
!
!     SOIL ADSORBED ANIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L.OR.AEC_vr(L0,NY,NX).LE.ZERO)THEN
    FAO=FO
  ELSE
    FAO=AMIN1(0.5_r8,FO*AEC_vr(L1,NY,NX)/AEC_vr(L0,NY,NX))
  ENDIF

  DO NTX=idx_cation_end+1,idx_end
    FXXA                        = FAO*trcx_solml_vr(NTX,L0,NY,NX)
    trcx_solml_vr(NTX,L1,NY,NX) = trcx_solml_vr(NTX,L1,NY,NX)+FXXA
    trcx_solml_vr(NTX,L0,NY,NX) = trcx_solml_vr(NTX,L0,NY,NX)-FXXA
  ENDDo
!
!     SOIL PRECIPITATES IN BAND, NON-BAND
!
  DO NTP=idsp_beg,idsp_end
    FXP                           = AMIN1(FX*trcp_saltpml_vr(NTP,L,NY,NX),trcp_saltpml_vr(NTP,L0,NY,NX))
    trcp_saltpml_vr(NTP,L1,NY,NX) = trcp_saltpml_vr(NTP,L1,NY,NX)+FXP
    trcp_saltpml_vr(NTP,L0,NY,NX) = trcp_saltpml_vr(NTP,L0,NY,NX)-FXP
  ENDDO

  end subroutine MoveBandSolute

!------------------------------------------------------------------------------------------


  Subroutine MoveDisolvGas(L0,L1,NY,NX,FX,FWO)
  implicit none
  integer, intent(in) :: L0,L1,NY,NX
  real(r8), intent(in) :: FX,FWO

  real(r8) :: FXG
  integer  :: NTG
!
!     SOIL GASEOUS GASES
! exclude NH3B, NH3
  DO NTG=idg_beg,idg_end-2
    FXG                        = FWO*trc_gasml_vr(NTG,L0,NY,NX)
    trc_gasml_vr(NTG,L1,NY,NX) = trc_gasml_vr(NTG,L1,NY,NX)+FXG
    trc_gasml_vr(NTG,L0,NY,NX) = trc_gasml_vr(NTG,L0,NY,NX)-FXG

    FXG                        = FWO*trc_solml_vr(NTG,L0,NY,NX)
    trc_solml_vr(NTG,L1,NY,NX) = trc_solml_vr(NTG,L1,NY,NX)+FXG
    trc_solml_vr(NTG,L0,NY,NX) = trc_solml_vr(NTG,L0,NY,NX)-FXG

  ENDDO
! add NH3
  NTG                        = idg_NH3
  FXG                        = FWO*trc_gasml_vr(NTG,L0,NY,NX)
  trc_gasml_vr(NTG,L1,NY,NX) = trc_gasml_vr(NTG,L1,NY,NX)+FXG
  trc_gasml_vr(NTG,L0,NY,NX) = trc_gasml_vr(NTG,L0,NY,NX)-FXG

  end Subroutine MoveDisolvGas

!------------------------------------------------------------------------------------------

  subroutine MoveFertSalt(L,L0,L1,NY,NX,FX,FWO)

  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in):: FX,FWO
  real(r8) :: FXZ,FXNUT,FXZN
  integer  :: NTS,NTSA,NTF

! begin_execution
  DO NTF=ifertn_beg,ifertn_end
    FXZN                        = AMIN1(FX*FertN_soil_vr(NTF,L,NY,NX),FertN_soil_vr(NTF,L0,NY,NX))
    FertN_soil_vr(NTF,L1,NY,NX) = FertN_soil_vr(NTF,L1,NY,NX)+FXZN
    FertN_soil_vr(NTF,L0,NY,NX) = FertN_soil_vr(NTF,L0,NY,NX)-FXZN
  ENDDO

  IF (L0>0) then
    DO NTF=ifertnb_beg,ifertnb_end
      FXZN                        = AMIN1(FX*FertN_Band_vr(NTF,L,NY,NX),FertN_Band_vr(NTF,L0,NY,NX))
      FertN_Band_vr(NTF,L1,NY,NX) = FertN_Band_vr(NTF,L1,NY,NX)+FXZN
      FertN_Band_vr(NTF,L0,NY,NX) = FertN_Band_vr(NTF,L0,NY,NX)-FXZN
    ENDDO
  endif
!
!     SOIL N,P SOLUTES IN BAND, NON-BAND
!
  DO NTS=ids_nuts_beg,ids_nuts_end
    if(NTS/=ids_H2PO4B .and. NTS/=ids_H1PO4B)THEN
      FXNUT                      = FWO*trc_solml_vr(NTS,L0,NY,NX)
      trc_solml_vr(NTS,L1,NY,NX) = trc_solml_vr(NTS,L1,NY,NX)+FXNUT
      trc_solml_vr(NTS,L0,NY,NX) = trc_solml_vr(NTS,L0,NY,NX)-FXNUT
    ENDIF
  ENDDO

!
!     SOIL SALT SOLUTES
!
  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      FXZ                             = FWO*trcSalt_solml_vr(NTSA,L0,NY,NX)
      trcSalt_solml_vr(NTSA,L1,NY,NX) = trcSalt_solml_vr(NTSA,L1,NY,NX)+FXZ
      trcSalt_solml_vr(NTSA,L0,NY,NX) = trcSalt_solml_vr(NTSA,L0,NY,NX)-FXZ
    ENDDO
  ENDIF
  end subroutine MoveFertSalt

!------------------------------------------------------------------------------------------

  subroutine MoveFertMinerals(L,L0,L1,NY,NX,FX,FO)
  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8),intent(in) :: FX,FO
  real(r8) :: FHO,FBO

  real(r8) :: FXSAND,FXSILT,FXCLAY,FXROCK
  real(r8) :: FXWDNHB,FXDPNHB,FXWDNOB
  real(r8) :: FXDPNOB,FXWDPOB,FXDPPOB
  real(r8) :: FXVOLWH,FXVOLIH,FXVOLAH
  IF(DLYR(3,L1,NY,NX).GT.ZERO.AND.DLYR(3,L0,NY,NX).GT.ZERO)THEN
!
!     SOIL FERTILIZER BANDS
!
    IF(IFNHB(NY,NX).EQ.1.AND.ROWN(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CumDepz2LayerBot_vr(L-1,NY,NX).LT.BandDepthNH4_col(NY,NX))THEN
        WDNHBDL                   = BandWidthNH4_vr(L,NY,NX)*DLYR(3,L,NY,NX)
        WDNHBD0                   = BandWidthNH4_vr(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDNHBD1                   = BandWidthNH4_vr(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDNHB                   = AMIN1(FX*WDNHBDL,WDNHBD0)
        WDNHBD1                   = WDNHBD1+FXWDNHB
        WDNHBD0                   = WDNHBD0-FXWDNHB
        BandWidthNH4_vr(L1,NY,NX) = WDNHBD1/DLYR(3,L1,NY,NX)
        BandWidthNH4_vr(L0,NY,NX) = WDNHBD0/DLYR(3,L0,NY,NX)
        IF(CumDepz2LayerBot_vr(L,NY,NX).GE.BandDepthNH4_col(NY,NX))THEN
          FXDPNHB                       = AMIN1(FX*BandThicknessNH4_vr(L,NY,NX),BandThicknessNH4_vr(L0,NY,NX))
          BandThicknessNH4_vr(L1,NY,NX) = BandThicknessNH4_vr(L1,NY,NX)+FXDPNHB
          BandThicknessNH4_vr(L0,NY,NX) = BandThicknessNH4_vr(L0,NY,NX)-FXDPNHB
        ENDIF
        trcs_VLN_vr(ids_NH4B,L1,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthNH4_vr(L1,NY,NX) &
          /ROWN(NY,NX)*BandThicknessNH4_vr(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        trcs_VLN_vr(ids_NH4B,L0,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthNH4_vr(L0,NY,NX) &
          /ROWN(NY,NX)*BandThicknessNH4_vr(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        trcs_VLN_vr(ids_NH4,L1,NY,NX) = 1.0_r8-trcs_VLN_vr(ids_NH4B,L1,NY,NX)
        trcs_VLN_vr(ids_NH4,L0,NY,NX) = 1.0_r8-trcs_VLN_vr(ids_NH4B,L0,NY,NX)

        trcs_VLN_vr(idg_NH3B,L1,NY,NX) = trcs_VLN_vr(ids_NH4B,L1,NY,NX)
        trcs_VLN_vr(idg_NH3B,L0,NY,NX) = trcs_VLN_vr(ids_NH4B,L0,NY,NX)
        trcs_VLN_vr(idg_NH3,L1,NY,NX)  = trcs_VLN_vr(ids_NH4,L1,NY,NX)
        trcs_VLN_vr(idg_NH3,L0,NY,NX)  = trcs_VLN_vr(ids_NH4,L0,NY,NX)
      ENDIF
    ENDIF
    IF(IFNOB(NY,NX).EQ.1.AND.ROWO(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CumDepz2LayerBot_vr(L-1,NY,NX).LT.BandDepthNO3_col(NY,NX))THEN
        WDNOBDL                   = BandWidthNO3_vr(L,NY,NX)*DLYR(3,L,NY,NX)
        WDNOBD0                   = BandWidthNO3_vr(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDNOBD1                   = BandWidthNO3_vr(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDNOB                   = AMIN1(FX*WDNOBDL,WDNOBD0)
        WDNOBD1                   = WDNOBD1+FXWDNOB
        WDNOBD0                   = WDNOBD0-FXWDNOB
        BandWidthNO3_vr(L1,NY,NX) = WDNOBD1/DLYR(3,L1,NY,NX)
        BandWidthNO3_vr(L0,NY,NX) = WDNOBD0/DLYR(3,L0,NY,NX)
        IF(CumDepz2LayerBot_vr(L,NY,NX).GE.BandDepthNO3_col(NY,NX))THEN
          FXDPNOB                       = AMIN1(FX*BandThicknessNO3_vr(L,NY,NX),BandThicknessNO3_vr(L0,NY,NX))
          BandThicknessNO3_vr(L1,NY,NX) = BandThicknessNO3_vr(L1,NY,NX)+FXDPNOB
          BandThicknessNO3_vr(L0,NY,NX) = BandThicknessNO3_vr(L0,NY,NX)-FXDPNOB
        ENDIF
        trcs_VLN_vr(ids_NO3B,L1,NY,NX)=AZMAX1(AMIN1(0.999_r8,BandWidthNO3_vr(L1,NY,NX) &
          /ROWO(NY,NX)*BandThicknessNO3_vr(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        trcs_VLN_vr(ids_NO3B,L0,NY,NX)=AZMAX1(AMIN1(0.999_r8,BandWidthNO3_vr(L0,NY,NX) &
          /ROWO(NY,NX)*BandThicknessNO3_vr(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        trcs_VLN_vr(ids_NO3,L1,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NO3B,L1,NY,NX)
        trcs_VLN_vr(ids_NO3,L0,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NO3B,L0,NY,NX)

        trcs_VLN_vr(ids_NO2,L1,NY,NX)  = trcs_VLN_vr(ids_NO3,L1,NY,NX)
        trcs_VLN_vr(ids_NO2,L0,NY,NX)  = trcs_VLN_vr(ids_NO3,L0,NY,NX)
        trcs_VLN_vr(ids_NO2B,L1,NY,NX) = trcs_VLN_vr(ids_NO3B,L1,NY,NX)
        trcs_VLN_vr(ids_NO2B,L0,NY,NX) = trcs_VLN_vr(ids_NO3B,L0,NY,NX)
      ENDIF
    ENDIF
    IF(IFPOB(NY,NX).EQ.1.AND.ROWP(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CumDepz2LayerBot_vr(L-1,NY,NX).LT.BandDepthPO4_col(NY,NX))THEN
        WDPOBDL                   = BandWidthPO4_vr(L,NY,NX)*DLYR(3,L,NY,NX)
        WDPOBD0                   = BandWidthPO4_vr(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDPOBD1                   = BandWidthPO4_vr(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDPOB                   = AMIN1(FX*WDPOBDL,WDPOBD0)
        WDPOBD1                   = WDPOBD1+FXWDPOB
        WDPOBD0                   = WDPOBD0-FXWDPOB
        BandWidthPO4_vr(L1,NY,NX) = WDPOBD1/DLYR(3,L1,NY,NX)
        BandWidthPO4_vr(L0,NY,NX) = WDPOBD0/DLYR(3,L0,NY,NX)
        IF(CumDepz2LayerBot_vr(L,NY,NX).GE.BandDepthPO4_col(NY,NX))THEN
          FXDPPOB                       = AMIN1(FX*BandThicknessPO4_vr(L,NY,NX),BandThicknessPO4_vr(L0,NY,NX))
          BandThicknessPO4_vr(L1,NY,NX) = BandThicknessPO4_vr(L1,NY,NX)+FXDPPOB
          BandThicknessPO4_vr(L0,NY,NX) = BandThicknessPO4_vr(L0,NY,NX)-FXDPPOB
        ENDIF
        trcs_VLN_vr(ids_H1PO4B,L1,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthPO4_vr(L1,NY,NX) &
          /ROWP(NY,NX)*BandThicknessPO4_vr(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        trcs_VLN_vr(ids_H1PO4B,L0,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthPO4_vr(L0,NY,NX) &
          /ROWP(NY,NX)*BandThicknessPO4_vr(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        trcs_VLN_vr(ids_H1PO4,L1,NY,NX)=1.0_r8-trcs_VLN_vr(ids_H1PO4B,L1,NY,NX)
        trcs_VLN_vr(ids_H1PO4,L0,NY,NX)=1.0_r8-trcs_VLN_vr(ids_H1PO4B,L0,NY,NX)

        trcs_VLN_vr(ids_H2PO4B,L1,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L1,NY,NX)
        trcs_VLN_vr(ids_H2PO4B,L0,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L0,NY,NX)
        trcs_VLN_vr(ids_H2PO4,L1,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L1,NY,NX)
        trcs_VLN_vr(ids_H2PO4,L0,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L0,NY,NX)
      ENDIF
    ENDIF
  ENDIF
!
!     SOIL MINERALS
!
  IF(L0.EQ.L.OR.SoiBulkDensityt0_vr(L0,NY,NX).LE.ZERO)THEN
    FBO=FX
  ELSE
    FBO=AMIN1(0.1,FX*SoiBulkDensityt0_vr(L1,NY,NX)/SoiBulkDensityt0_vr(L0,NY,NX))
  ENDIF
!     SoiBulkDensity_vr(L1,NY,NX)=(1.0-FO)*SoiBulkDensity_vr(L1,NY,NX)+FO*SoiBulkDensityt0_vr(L0,NY,NX)
  PH(L1,NY,NX)      = (1.0_r8-FO)*PH(L1,NY,NX)+FO*PH(L0,NY,NX)
  FXSAND            = FBO*SAND(L0,NY,NX)
  SAND(L1,NY,NX)    = SAND(L1,NY,NX)+FXSAND
  SAND(L0,NY,NX)    = SAND(L0,NY,NX)-FXSAND
  FXSILT            = FBO*SILT(L0,NY,NX)
  SILT(L1,NY,NX)    = SILT(L1,NY,NX)+FXSILT
  SILT(L0,NY,NX)    = SILT(L0,NY,NX)-FXSILT
  FXCLAY            = FBO*CLAY(L0,NY,NX)
  CLAY(L1,NY,NX)    = CLAY(L1,NY,NX)+FXCLAY
  CLAY(L0,NY,NX)    = CLAY(L0,NY,NX)-FXCLAY
  FXROCK            = FBO*ROCK_vr(L0,NY,NX)
  ROCK_vr(L1,NY,NX) = ROCK_vr(L1,NY,NX)+FXROCK
  ROCK_vr(L0,NY,NX) = ROCK_vr(L0,NY,NX)-FXROCK
!
!     SOIL WATER AND HEAT
!
  IF(SoilFracAsMacP_vr(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP_vr(L0,NY,NX).GT.ZERO)THEN
    IF(L0.EQ.L.OR.SoilFracAsMacPt0_vr(L0,NY,NX).LE.ZERO)THEN
      FHO=FO
    ELSE
      FHO=AMIN1(0.5_r8,FO*SoilFracAsMacPt0_vr(L1,NY,NX)/SoilFracAsMacPt0_vr(L0,NY,NX))
    ENDIF
    SoilFracAsMacP_vr(L1,NY,NX) = (1.0_r8-FO)*SoilFracAsMacP_vr(L1,NY,NX)+FO*SoilFracAsMacP_vr(L0,NY,NX)
    FXVOLWH                     = FHO*VLWatMacP_vr(L0,NY,NX)
    VLWatMacP_vr(L1,NY,NX)      = VLWatMacP_vr(L1,NY,NX)+FXVOLWH
    VLWatMacP_vr(L0,NY,NX)      = VLWatMacP_vr(L0,NY,NX)-FXVOLWH
    FXVOLIH                     = FHO*VLiceMacP_vr(L0,NY,NX)
    VLiceMacP_vr(L1,NY,NX)      = VLiceMacP_vr(L1,NY,NX)+FXVOLIH
    VLiceMacP_vr(L0,NY,NX)      = VLiceMacP_vr(L0,NY,NX)-FXVOLIH
    FXVOLAH                     = FHO*VLMacP_vr(L0,NY,NX)
    VLMacP_vr(L1,NY,NX)         = VLMacP_vr(L1,NY,NX)+FXVOLAH
    VLMacP_vr(L0,NY,NX)         = VLMacP_vr(L0,NY,NX)-FXVOLAH
  ENDIF
  end subroutine MoveFertMinerals

!------------------------------------------------------------------------------------------

  subroutine MoveHeatWat(L,L0,L1,NY,NX,FO,FX)
  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in) :: FO,FX
  real(r8) :: FWO,FXVOLW
  real(r8) :: FXENGY,FXVOLI,FXVLSoilMicP,FXVOLWX,FXVHCM
  real(r8) :: ENGY0,ENGY1

  IF(L0.EQ.L.OR.(L0>0 .and. POROSI_vr(L0,NY,NX).LE.ZERO))THEN
    FWO=FO
  ELSE
    FWO=AMIN1(0.5,FO*POROSI_vr(L1,NY,NX)/POROSI_vr(L0,NY,NX))
  ENDIF
!             FXSCNV=FWO*SatHydroCondVert_vr(L0,NY,NX)
!             SatHydroCondVert_vr(L1,NY,NX)=SatHydroCondVert_vr(L1,NY,NX)+FXSCNV
!             SatHydroCondVert_vr(L0,NY,NX)=SatHydroCondVert_vr(L0,NY,NX)-FXSCNV
!             FXSCNH=FWO*SatHydroCondHrzn_vr(L0,NY,NX)
!             SatHydroCondHrzn_vr(L1,NY,NX)=SatHydroCondHrzn_vr(L1,NY,NX)+FXSCNH
!             SatHydroCondHrzn_vr(L0,NY,NX)=SatHydroCondHrzn_vr(L0,NY,NX)-FXSCNH
  IF(L0.EQ.0)THEN
    FXVOLW=FX*AZMAX1(XVOLWP-VLWatheldCapSurf_col(NY,NX))
  ELSE
    FXVOLW=FWO*VLWatMicP_vr(L0,NY,NX)
  ENDIF
  VLWatMicP_vr(L1,NY,NX) = VLWatMicP_vr(L1,NY,NX)+FXVOLW
  VLWatMicP_vr(L0,NY,NX) = VLWatMicP_vr(L0,NY,NX)-FXVOLW
!     IF(VLiceMicP_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
  FXVOLI                 = FWO*VLiceMicP_vr(L0,NY,NX)
  VLiceMicP_vr(L1,NY,NX) = VLiceMicP_vr(L1,NY,NX)+FXVOLI
  VLiceMicP_vr(L0,NY,NX) = VLiceMicP_vr(L0,NY,NX)-FXVOLI
!     ENDIF
!     FXVOLA=FWO*VLMicP_vr(L0,NY,NX)
!     IF(L1.NE.NU(NY,NX))THEN
!     VLMicP_vr(L1,NY,NX)=VLMicP_vr(L1,NY,NX)+FXVOLA
!     ENDIF
!     IF(L0.NE.NU(NY,NX))THEN
!     VLMicP_vr(L0,NY,NX)=VLMicP_vr(L0,NY,NX)-FXVOLA
!     ENDIF
  FXVLSoilMicP                    = FWO*VLSoilMicP_vr(L0,NY,NX)
  VLSoilMicP_vr(L1,NY,NX)         = VLSoilMicP_vr(L1,NY,NX)+FXVLSoilMicP
  VLSoilMicP_vr(L0,NY,NX)         = VLSoilMicP_vr(L0,NY,NX)-FXVLSoilMicP
  FXVOLWX                         = FWO*VLWatMicPX_vr(L0,NY,NX)
  VLWatMicPX_vr(L1,NY,NX)         = VLWatMicPX_vr(L1,NY,NX)+FXVOLWX
  VLWatMicPX_vr(L0,NY,NX)         = VLWatMicPX_vr(L0,NY,NX)-FXVOLWX
  FXVHCM                          = FWO*VHeatCapacitySoilM_vr(L0,NY,NX)
  VHeatCapacitySoilM_vr(L1,NY,NX) = VHeatCapacitySoilM_vr(L1,NY,NX)+FXVHCM
  VHeatCapacitySoilM_vr(L0,NY,NX) = VHeatCapacitySoilM_vr(L0,NY,NX)-FXVHCM
  FXENGY                          = TKS_vr(L0,NY,NX)*(FXVHCM+cpw*FXVOLW+cpi*FXVOLI)
  ENGY1                           = VHeatCapacity_vr(L1,NY,NX)*TKS_vr(L1,NY,NX)+FXENGY
  ENGY0                           = VHeatCapacity_vr(L0,NY,NX)*TKS_vr(L0,NY,NX)-FXENGY
  VHeatCapacity_vr(L1,NY,NX)      = VHeatCapacity_vr(L1,NY,NX)+FXVHCM+cpw*FXVOLW+cpi*FXVOLI
  VHeatCapacity_vr(L0,NY,NX)      = VHeatCapacity_vr(L0,NY,NX)-FXVHCM-cpw*FXVOLW-cpi*FXVOLI
  IF(VHeatCapacity_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L1,NY,NX)=ENGY1/VHeatCapacity_vr(L1,NY,NX)
  ELSE
    TKS_vr(L1,NY,NX)=TKS_vr(L,NY,NX)
  ENDIF
  TCS(L1,NY,NX)=units%Kelvin2Celcius(TKS_vr(L1,NY,NX))
  IF(VHeatCapacity_vr(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L0,NY,NX)=ENGY0/VHeatCapacity_vr(L0,NY,NX)
  ELSE
    TKS_vr(L0,NY,NX)=TKS_vr(L,NY,NX)
  ENDIF
  TCS(L0,NY,NX)=units%Kelvin2Celcius(TKS_vr(L0,NY,NX))
  end subroutine MoveHeatWat

end module SoilLayerDynMod
