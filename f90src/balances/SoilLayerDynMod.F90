module SoilLayerDynMod
! Description:
! subroutines to do soil relayering
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use RootDataType
  use GridDataType
  use minimathmod, only : AZMAX1
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
  use SoilPhysDataType
  use SurfSoilDataType
  use SurfLitterDataType
  use EcosimConst
  use UnitMod     , only : units
  use EcoSIMConfig, only : ndbiomcp => ndbiomcpc
  USE TFlxTypeMod , ONLY : TSEDER,TDLYXF,TDYLXC,TDVOLI,TDORGC
implicit none


  private

  character(len=*), parameter :: mod_filename = __FILE__

  real(r8) :: DDLYRX(3)
  real(r8) :: XVOLWP,WDPOBDL,WDNOBD1,WDPOBD0,WDPOBD1
  real(r8) :: WDNHBDL,WDNHBD0,WDNHBD1,WDNOBDL,WDNOBD0
  real(r8), PARAMETER :: ZEROC=0.1E-03_r8
  integer, private, parameter :: ist_water=0
  integer, private, parameter :: ist_soil=1
  public :: RelayerSoilProfile
  contains


  subroutine RelayerSoilProfile(NY,NX,DORGC,DVLiceMicP,UDVLiceMicP,UDLYXF)
  !
  !Description:
  !relayer the soil profiles
  implicit none
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: DORGC(JZ,JY,JX)  !change in organic matter, initial-final
  REAL(R8),INTENT(IN) :: DVLiceMicP(JZ,JY,JX)  !change in ice volume, initial-final
  real(r8),intent(inout):: UDVLiceMicP,UDLYXF

  real(r8) :: CDPTHY(0:JZ,JY,JX),CDPTHX(JZ,JY,JX)
  integer :: IFLGL(0:JZ,6)  !flag for soil thickness change
  integer :: NN,K,M,N,NR,NZ,L
  integer :: L0,L1,NUX,ICHKLX,ICHKL,NGL
  real(r8) :: FX,FY
  real(r8) :: FBO
  real(r8) :: FHO,FWO,FXVOLW
  real(r8) :: FO
  real(r8) :: ENGY0X,ENGY0,ENGY1X,ENGY1
  !     begin_execution
  !     SOIL SUBSIDENCE
  !
  IF(IERSNG < 0)return
  !soil relayering can occur due to freeze-thaw, soc change, and erosion
  !  
    IF(SoiBulkDensity(NU(NY,NX),NY,NX).LE.ZERO)THEN
      ICHKLX=ist_water      !surface is water layer
    ELSE
      !it is a soil column
      ICHKLX=ist_soil
    ENDIF

    call SoilSubsidence(ICHKLX,NY,NX,DORGC,DVLiceMicP,UDLYXF,UDVLiceMicP,CDPTHX,CDPTHY,IFLGL)

    !
    !     RECALCULATE SOIL LAYER THICKNESS
    !
    ICHKL=0
    D245: DO L=NU(NY,NX),NL(NY,NX)-1
      D230: DO NN=1,3

        call getTSlyrthicks(NN,L,NY,NX,ICHKL,NUX,CDPTHX,CDPTHY,IFLGL)

        !
        !     TRANSFER STATE VARIABLES BETWEEN LAYERS
        !
        !     IF(IFLGL(L,NN).EQ.1)THEN
        IF(ABS(DDLYRX(NN)).GT.ZERO)THEN
!L0,L1: target and source layers
          call getFLs(L,NN,NY,NX,NUX,FX,FO,L1,L0,IFLGL)

          IF(FX.GT.ZERO)THEN
            IFLGS(NY,NX)=1
            FY=1.0_r8-FX
            IF(FY.LE.ZERO2)FY=0.0_r8
            IF(SoiBulkDensity(L0,NY,NX).LE.ZERO)THEN
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
                IF(SoiBulkDensity(L0,NY,NX).LE.ZERO.AND.SoiBulkDensity(L1,NY,NX).LE.ZERO &
                  .AND.VLWatMicP(L0,NY,NX)+VLiceMicP(L0,NY,NX).LE.ZEROS(NY,NX))THEN
                  CumDepth2LayerBottom(L1,NY,NX)=CumDepth2LayerBottom(L0,NY,NX)
                  CDPTHY(L1,NY,NX)=CDPTHY(L0,NY,NX)
                ENDIF
              ENDIF

            ENDIF
!
!     RESET LOWER LAYER NUMBER WITH EROSION
!
            IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
              IF(L.EQ.NL(NY,NX).AND.DLYR(3,L,NY,NX).GT.DLYRI(3,L,NY,NX))THEN
                NL(NY,NX)=MIN(NLI(NY,NX),NL(NY,NX)+1)
              ENDIF
              IF(L.EQ.NL(NY,NX)-1.AND.CumDepth2LayerBottom(NL(NY,NX),NY,NX)-CumDepth2LayerBottom(L,NY,NX).LE.ZEROC)THEN
                CumDepth2LayerBottom(L,NY,NX)=CumDepth2LayerBottom(L,NY,NX)+DLYR(3,NL(NY,NX),NY,NX)
                CDPTHY(L,NY,NX)=CDPTHY(L,NY,NX)+DLYR(3,NL(NY,NX),NY,NX)
                CumDepth2LayerBottom(NL(NY,NX),NY,NX)=CumDepth2LayerBottom(L,NY,NX)
                CDPTHY(NL(NY,NX),NY,NX)=CDPTHY(L,NY,NX)
                DLYR(3,NL(NY,NX),NY,NX)=0.0_r8
                NL(NY,NX)=L
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO D230
    ENDDO D245
  end subroutine RelayerSoilProfile
!------------------------------------------------------------------------------------------
  subroutine SoilSubsidence(ICHKLX,NY,NX,DORGC,DVLiceMicP,UDLYXF,UDVLiceMicP,&
    CDPTHX,CDPTHY,IFLGL)
!
! IFLGL: c1, ponding water, c2, pond disappear, c3, pond reappare, c4: freeze-thaw, c5: erosion, c6: som change
  implicit none
  integer, intent(in) :: ICHKLX   !surface layer type: 0 water, 1 soil
  integer, intent(in) :: NY,NX
  REAL(R8),INTENT(IN) :: DVLiceMicP(JZ,JY,JX)
  real(r8),intent(in) :: DORGC(JZ,JY,JX)
  real(r8), intent(inout) :: UDLYXF,UDVLiceMicP
  real(r8), intent(out) :: CDPTHX(JZ,JY,JX)
  real(r8), intent(inout) :: CDPTHY(0:JZ,JY,JX)
  integer, intent(inout) :: IFLGL(0:JZ,6)
  real(r8) :: DDLYX(0:JZ,6)
  real(r8) :: DDLYR(0:JZ,6)
  integer :: LX,LY,LL,NN,L
  real(r8) :: DDLYXC,DDLYXE
  real(r8) :: DDLYXF
  real(r8) :: DENSJ,DPTWI,DDLYXP
  integer, parameter :: ich_watlev =1
  integer, parameter :: ich_frzthaw=4  
  integer, parameter :: ich_erosion=5
  integer, parameter :: ich_socloss=6
! begin_execution

! starting from bottom up
  D225: DO LX=NL(NY,NX),NU(NY,NX),-1
    !make a copy of the depth, bottom of the layer
    CDPTHX(LX,NY,NX)=CumDepth2LayerBottom(LX,NY,NX)
    CDPTHY(LX,NY,NX)=CumDepth2LayerBottom(LX,NY,NX)
    !
    !     POND, from water to soil
    !
    ! compute the change
    IF(SoiBulkDensity(LX,NY,NX).LE.ZERO)THEN
      ! current layer is water, lower layer is soil, or surface is soil
      IF(SoiBulkDensity(LX+1,NY,NX).GT.ZERO.OR.ICHKLX.EQ.ist_soil)THEN
        !layer below is soil, layer above is also soil
        DDLYXP=DLYR(3,LX,NY,NX)-(VLWatMicP(LX,NY,NX)+VLiceMicP(LX,NY,NX))/AREA(3,LX,NY,NX)
        DDLYX(LX,ich_watlev)=DDLYXP+DDLYX(LX+1,ich_watlev)
        DDLYR(LX,ich_watlev)=DDLYX(LX+1,ich_watlev)
        IFLGL(LX,ich_watlev)=2
      ELSE
        !next layer, current and top are all water
        !DDLYXP: soil equivalent depth
        !DLYRI: initial soil layer thickness, [m]
        DDLYXP=DLYRI(3,LX,NY,NX)-(VLWatMicP(LX,NY,NX)+VLiceMicP(LX,NY,NX))/AREA(3,LX,NY,NX)

        !DPTWI: water+ice total thickness layer below
        DPTWI=(VLWatMicP(LX+1,NY,NX)+VLiceMicP(LX+1,NY,NX))/AREA(3,LX,NY,NX)
        !there is expansion in layer LX, or next layer has
        IF(DDLYXP.LT.-ZERO.OR.DPTWI.GT.ZERO)THEN
          DDLYX(LX,ich_watlev)=DDLYXP+DDLYX(LX+1,ich_watlev)
          DDLYR(LX,ich_watlev)=AMIN1(DDLYX(LX+1,ich_watlev),DPTWI)
          IF(DPTWI.GT.ZERO)THEN
! already has solid thickness
            IFLGL(LX,ich_watlev)=1
          ELSE
! has no solid thickness
            IFLGL(LX,ich_watlev)=2
          ENDIF
        ELSE
          !shrink
          DDLYXP=DLYR(3,LX,NY,NX)-(VLWatMicP(LX,NY,NX)+VLiceMicP(LX,NY,NX))/AREA(3,LX,NY,NX)
          DDLYX(LX,ich_watlev)=DDLYXP+DDLYX(LX+1,ich_watlev)
          DDLYR(LX,ich_watlev)=DDLYX(LX+1,ich_watlev)
          IFLGL(LX,ich_watlev)=2
        ENDIF
      ENDIF

      !surface layer or layer above is still soil
      IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity(LX-1,NY,NX).GT.ZERO)THEN
        DDLYX(LX-1,ich_watlev)=DDLYX(LX,ich_watlev)
        DDLYR(LX-1,ich_watlev)=DDLYX(LX,ich_watlev)
        IFLGL(LX-1,ich_watlev)=1
      ENDIF
      DDLYX(LX,ich_frzthaw)=0.0_r8
      DDLYR(LX,ich_frzthaw)=0.0_r8
      IFLGL(LX,ich_frzthaw)=0
      DDLYX(LX,ich_erosion)=0.0_r8
      DDLYR(LX,ich_erosion)=0.0_r8
      IFLGL(LX,ich_erosion)=0
      DDLYX(LX,ich_socloss)=0.0_r8
      DDLYR(LX,ich_socloss)=0.0_r8
      IFLGL(LX,ich_socloss)=0
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
      !DDLYXF: added thickness change
      ! 4: due to freeze-thaw
      IF(ABS(DVLiceMicP(LX,NY,NX)).GT.ZEROS(NY,NX))THEN
        DENSJ=1._r8-DENSICE
        DDLYXF=DVLiceMicP(LX,NY,NX)*DENSJ/AREA(3,LX,NY,NX)
        !bottom layer
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,ich_frzthaw)=DDLYXF
          DDLYR(LX,ich_frzthaw)=0.0_r8
          IFLGL(LX,ich_frzthaw)=0
         !layer not at the bottom
        ELSE
          DDLYX(LX,ich_frzthaw)=DDLYXF+DDLYX(LX+1,ich_frzthaw)
          DDLYR(LX,ich_frzthaw)=DDLYX(LX+1,ich_frzthaw)
          !    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
          IFLGL(LX,ich_frzthaw)=0
          !top soil layer or water layer
          IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity(LX-1,NY,NX).LE.ZERO)THEN
            DDLYX(LX-1,ich_frzthaw)=DDLYX(LX,ich_frzthaw)
            DDLYR(LX-1,ich_frzthaw)=DDLYX(LX,ich_frzthaw)
            !    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
            IFLGL(LX-1,ich_frzthaw)=0
          ENDIF
        ENDIF
      ELSE
        ! no change in ice volume
        DDLYXF=0.0_r8
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,ich_frzthaw)=0.0_r8
          DDLYR(LX,ich_frzthaw)=0.0_r8
          IFLGL(LX,ich_frzthaw)=0
        ELSE
          DDLYX(LX,ich_frzthaw)=DDLYX(LX+1,ich_frzthaw)
          DDLYR(LX,ich_frzthaw)=DDLYX(LX+1,ich_frzthaw)
          IFLGL(LX,ich_frzthaw)=0
          IF(LX.EQ.NU(NY,NX))THEN
            DDLYX(LX-1,ich_frzthaw)=DDLYX(LX,ich_frzthaw)
            DDLYR(LX-1,ich_frzthaw)=DDLYX(LX,ich_frzthaw)
            IFLGL(LX-1,ich_frzthaw)=0
          ENDIF
        ENDIF
      ENDIF

      !total change in ice volume
      TDVOLI=TDVOLI+DVLiceMicP(LX,NY,NX)
      TDLYXF=TDLYXF+DDLYXF
      UDVLiceMicP=UDVLiceMicP+DVLiceMicP(LX,NY,NX)
      UDLYXF=UDLYXF+DDLYXF
      !
      !     EROSION
      !
      IF((IERSNG.EQ.1.OR.IERSNG.EQ.3).AND.ABS(TSEDER(NY,NX)).GT.ZEROS(NY,NX))THEN
        IF(LX.EQ.NL(NY,NX))THEN
!  5: due to sediment erosion
!         total soil layer reduction due to erosion
          DDLYXE=-TSEDER(NY,NX)/(SoilMicPMassLayerMX(NY,NX)/VLSoilPoreMicP(NU(NY,NX),NY,NX))
          DDLYX(LX,ich_erosion)=DDLYXE
          DDLYR(LX,ich_erosion)=DDLYXE
          IFLGL(LX,ich_erosion)=1
        ELSE
          DDLYX(LX,ich_erosion)=0.0_r8
          DDLYR(LX,ich_erosion)=0.0_r8
          IFLGL(LX,ich_erosion)=0
        ENDIF
      ELSE
          DDLYX(LX,ich_erosion)=0.0_r8
          DDLYR(LX,ich_erosion)=0.0_r8
          IFLGL(LX,ich_erosion)=0
      ENDIF

      !
      ! SOC GAIN OR LOSS
      ! SoilFracAsMacP: macropore fraction
      ! DDLYXC: soil thickness added due to change in organic matter,
      ! keeping macropore fraction
      ! SoiBulkDensityt0: initial bulk density,
      IF((IERSNG.EQ.2.OR.IERSNG.EQ.3).AND.ABS(DORGC(LX,NY,NX)).GT.ZEROS(NY,NX))THEN
        DDLYXC=MWC2Soil*DORGC(LX,NY,NX)/((1.0_r8-SoilFracAsMacP(LX,NY,NX))*SoiBulkDensityt0(LX,NY,NX))/AREA(3,LX,NY,NX)
        ! obtain diagnostics only for NX==1
        IF(NX.EQ.1)THEN
          TDORGC=TDORGC+DORGC(LX,NY,NX)
          TDYLXC=TDYLXC+DDLYXC
        ENDIF
        ! bottom layer
        IF(LX.EQ.NL(NY,NX).OR.SoiBulkDensity(LX+1,NY,NX).LE.ZERO)THEN
          DDLYX(LX,ich_socloss)=DDLYXC
          DDLYR(LX,ich_socloss)=0.0_r8
          IFLGL(LX,ich_socloss)=1
        ELSE
          DDLYX(LX,ich_socloss)=DDLYXC+DDLYX(LX+1,ich_socloss)
          DDLYR(LX,ich_socloss)=DDLYX(LX+1,ich_socloss)+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
          IFLGL(LX,ich_socloss)=1
          ! top layer
          IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity(LX-1,NY,NX).LE.ZERO)THEN
            DDLYX(LX-1,ich_socloss)=DDLYX(LX,ich_socloss)
            DDLYR(LX-1,ich_socloss)=DDLYX(LX,ich_socloss)
            !    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
            IFLGL(LX-1,ich_socloss)=1
          ENDIF
        ENDIF
      ELSE
        ! bottom layer
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,ich_socloss)=0.0_r8
          DDLYR(LX,ich_socloss)=0.0_r8
          IFLGL(LX,ich_socloss)=0
        ELSE
          DDLYX(LX,ich_socloss)=DDLYX(LX+1,ich_socloss)
          DDLYR(LX,ich_socloss)=DDLYX(LX+1,ich_socloss)
          IFLGL(LX,ich_socloss)=0
        ENDIF
      ENDIF
      DDLYX(LX,ich_watlev)=0.0_r8
      DDLYR(LX,ich_watlev)=0.0_r8
      IFLGL(LX,ich_watlev)=0

    ENDIF

    ! apply the change
      !
      !     RESET SOIL LAYER DEPTHS
      !
    D200: DO NN=1,6
      !     IF(ABS(DDLYX(LX,NN)).GT.ZERO)IFLGS(NY,NX)=1
!     c2, and c3 are not implemented
      IF(NN.NE.2.AND.NN.NE.3)THEN
        !
        !     POND
        !
        IF(SoiBulkDensity(LX,NY,NX).LE.ZERO)THEN          
          ! there are some changes
          IF(IFLGL(LX,NN).NE.0)THEN
            CumDepth2LayerBottom(LX,NY,NX)=CumDepth2LayerBottom(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
            !  not top layer
            IF(LX.NE.NU(NY,NX).AND.IFLGL(LX,ich_watlev).EQ.2)THEN
              DO  LL=LX-1,0,-1
                CumDepth2LayerBottom(LL,NY,NX)=CumDepth2LayerBottom(LL,NY,NX)+DDLYX(LX,NN)
                CDPTHY(LL,NY,NX)=CDPTHY(LL,NY,NX)+DDLYX(LX,NN)
              ENDDO
              DDLYX(LX,NN)=0.0_r8
            ENDIF
            !  top layer
            IF(LX.EQ.NU(NY,NX))THEN
              CumDepth2LayerBottom(LX-1,NY,NX)=CumDepth2LayerBottom(LX,NY,NX)-(VLWatMicP(LX,NY,NX)+VLiceMicP(LX,NY,NX))/AREA(3,LX,NY,NX)
              CDPTHY(LX-1,NY,NX)=CDPTHY(LX,NY,NX)-(VLWatMicP(LX,NY,NX)+VLiceMicP(LX,NY,NX))/AREA(3,LX,NY,NX)
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
            CumDepth2LayerBottom(LX,NY,NX)=CumDepth2LayerBottom(LX,NY,NX)+DDLYR(LX,NN)
            !     CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
! top layer
            IF(LX.EQ.NU(NY,NX))THEN
              CumDepth2LayerBottom(LX-1,NY,NX)=CumDepth2LayerBottom(LX-1,NY,NX)+DDLYR(LX-1,NN)
              !     CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX-1,NN)

            ENDIF
          ENDIF
            !
            ! SET SURFACE ELEVATION FOR SOIL EROSION
            !
          IF(NN.EQ.ich_erosion.AND.IFLGL(LX,NN).EQ.1)THEN
            CumDepth2LayerBottom(LX,NY,NX)=CumDepth2LayerBottom(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)

            IF(LX.EQ.NU(NY,NX))THEN
              CumDepth2LayerBottom(LX-1,NY,NX)=CumDepth2LayerBottom(LX-1,NY,NX)+DDLYR(LX,NN)
              CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX,NN)
            ENDIF
          ENDIF
          !
          ! SET SOIL LAYER DEPTHS FOR CHANGES IN SOC
          !
          IF(NN.EQ.ich_socloss.AND.IFLGL(LX,NN).EQ.1)THEN
            CumDepth2LayerBottom(LX,NY,NX)=CumDepth2LayerBottom(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)

            IF(LX.EQ.NU(NY,NX).OR.SoiBulkDensity(LX-1,NY,NX).LE.ZERO)THEN
              CumDepth2LayerBottom(LX-1,NY,NX)=CumDepth2LayerBottom(LX-1,NY,NX)+DDLYR(LX-1,NN)
              CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX-1,NN)

              IF(SoiBulkDensity(LX-1,NY,NX).LE.ZERO)THEN
                DO  LY=LX-2,0,-1
                  IF(SoiBulkDensity(LY+1,NY,NX).LE.ZERO)THEN
                    CumDepth2LayerBottom(LY,NY,NX)=CumDepth2LayerBottom(LY,NY,NX)+DDLYR(LX-1,NN)
                    CDPTHY(LY,NY,NX)=CDPTHY(LY,NY,NX)+DDLYR(LX-1,NN)
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D200
    VLSoilMicP(LX,NY,NX)=VLSoilPoreMicP(LX,NY,NX)
  ENDDO D225
  VLSoilMicP(0,NY,NX)=VLWatMicP(0,NY,NX)+VLiceMicP(0,NY,NX)
  end subroutine SoilSubsidence

!------------------------------------------------------------------------------------------


  subroutine getFLs(L,NN,NY,NX,NUX,FX,FO,L1,L0,IFLGL)

  implicit none
  integer, intent(in) :: L,NN,NY,NX,NUX
  real(r8), intent(out) :: FX,FO
  integer, intent(out) :: L1   !source
  integer, intent(out) :: L0   !target
  integer, intent(in) :: IFLGL(0:JZ,6)
  real(r8) :: DPTWI

! begin_execution
  IF(DDLYRX(NN).GT.ZERO)THEN
!
!     LAYERS MOVE DOWN (FX>0.), pond drys up
!
    IF(IFLGL(L,2).EQ.0)THEN
      L1=L
      L0=L+1
    ELSE
      L1=NU(NY,NX)
      L0=NUX
    ENDIF

    IF((SoiBulkDensity(L,NY,NX).LE.ZERO.AND.IFLGL(L,1).EQ.2) &
        .OR.(DLYR(3,L0,NY,NX).LE.ZEROC.AND.IFLGL(L,6).EQ.1))THEN
      FX=1.0_r8
      FO=1.0_r8
    ELSE
      IF(SoiBulkDensity(L0,NY,NX).LE.ZERO)THEN
        DPTWI=(VLWatMicP(L0,NY,NX)+VLiceMicP(L0,NY,NX))/AREA(3,L0,NY,NX)
        IF(DPTWI.GT.ZERO)THEN
          FX=AMIN1(1.0_r8,DDLYRX(NN)/DPTWI)
          FO=FX
        ELSE
          FX=0.0_r8
          FO=0.0_r8
        ENDIF
      ELSE
        IF(DLYR(3,L,NY,NX).GT.ZERO.AND.DLYR(3,L0,NY,NX).GT.ZERO)THEN
          FX=AMIN1(1.0,DDLYRX(NN)/DLYR(3,L,NY,NX))
          FO=AMIN1(1.0,DDLYRX(NN)/DLYR(3,L0,NY,NX))
        ELSE
          FX=0.0_r8
          FO=0.0_r8
        ENDIF
      ENDIF
    ENDIF
  ELSE
!
!     LAYERS MOVE UP (FX<=0.), pond builds up
!
    IF(IFLGL(L,3).EQ.0)THEN
      L1=L+1
      L0=L
    ELSE
      L1=NU(NY,NX)
      L0=0
    ENDIF
    IF(SoiBulkDensity(L0,NY,NX).LE.ZERO)THEN
      DPTWI=(VLWatMicP(L0,NY,NX)+VLiceMicP(L0,NY,NX))/AREA(3,L0,NY,NX)
      IF(DPTWI.GT.ZERO)THEN
        FX=AMIN1(1.0,-DDLYRX(NN)/DPTWI)
        FO=FX
      ELSE
        FX=0.0_r8
        FO=0.0_r8
      ENDIF
    ELSE
      IF(DLYR(3,L,NY,NX).GT.ZERO.AND.DLYR(3,L0,NY,NX).GT.ZERO)THEN
        FX=AMIN1(1.0,-DDLYRX(NN)/DLYR(3,L,NY,NX))
        FO=AMIN1(1.0,-DDLYRX(NN)/DLYR(3,L0,NY,NX))
      ELSE
        FX=0.0_r8
        FO=0.0_r8
      ENDIF
    ENDIF
  ENDIF
  end subroutine getFLs
!------------------------------------------------------------------------------------------

  subroutine getTSlyrthicks(NN,L,NY,NX,ICHKL,NUX,CDPTHX,CDPTHY,IFLGL)
  implicit none
  integer, intent(in) :: NN,L,NY,NX
  integer, intent(inout) :: ICHKL
  integer, intent(out) :: NUX
  real(r8), intent(in) :: CDPTHX(JZ,JY,JX)
  real(r8), intent(inout) :: CDPTHY(0:JZ,JY,JX)
  integer, intent(inout) :: IFLGL(0:JZ,6)
  real(r8) :: DDLYRY(JZ)
  real(r8) :: DLYR0
  real(r8) :: DLYRXX
  integer :: LL

! begin_execution
  IF(NN.EQ.1)THEN
    DLYR(3,L,NY,NX)=CumDepth2LayerBottom(L,NY,NX)-CumDepth2LayerBottom(L-1,NY,NX)
    DLYRXX=DLYR(3,L,NY,NX)
    IF(IFLGL(L,1).EQ.0.AND.IFLGL(L+1,1).NE.0)THEN
      DDLYRX(NN)=0.0_r8
      IF(SoiBulkDensity(L,NY,NX).LE.ZERO)THEN
        DDLYRY(L)=DLYRI(3,L,NY,NX)-DLYR(3,L,NY,NX)
      ELSE
        DDLYRY(L)=0.0_r8
      ENDIF
      ICHKL=1
    ELSEIF(IFLGL(L,1).EQ.2.AND.(IFLGL(L+1,1).EQ.0 &
      .OR.DLYR(3,L,NY,NX).LE.DLYRI(3,L,NY,NX)))THEN
      DDLYRX(NN)=0.0_r8
      IF(L.EQ.NU(NY,NX).OR.ICHKL.EQ.0)THEN
        DDLYRY(L)=0.0_r8
      ELSE
        DDLYRY(L)=DDLYRY(L-1)
      ENDIF
      IF(IFLGL(L,1).EQ.2.AND.IFLGL(L+1,1).EQ.0)ICHKL=0
    ELSE
      IF(ICHKL.EQ.0)THEN
        DDLYRX(NN)=DLYRI(3,L,NY,NX)-DLYR(3,L,NY,NX)
        DDLYRY(L)=DDLYRX(NN)
      ELSE
        DDLYRX(NN)=0.0_r8
        DDLYRY(L)=DDLYRY(L-1)
      ENDIF
    ENDIF
    CumDepth2LayerBottom(L,NY,NX)=CumDepth2LayerBottom(L,NY,NX)+DDLYRY(L)
  !     CDPTHY(L,NY,NX)=CDPTHY(L,NY,NX)+DDLYRY(L)
    DLYR(3,L,NY,NX)=DLYR(3,L,NY,NX)+DDLYRY(L)
    DPTH(L,NY,NX)=0.5_r8*(CumDepth2LayerBottom(L,NY,NX)+CumDepth2LayerBottom(L-1,NY,NX))
    CDPTHZ(L,NY,NX)=CumDepth2LayerBottom(L,NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
    IF(L.EQ.NL(NY,NX)-1)THEN
      DLYR(3,L+1,NY,NX)=CumDepth2LayerBottom(L+1,NY,NX)-CumDepth2LayerBottom(L,NY,NX)
      DPTH(L+1,NY,NX)=0.5_r8*(CumDepth2LayerBottom(L+1,NY,NX)+CumDepth2LayerBottom(L,NY,NX))
      CDPTHZ(L+1,NY,NX)=CumDepth2LayerBottom(L+1,NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)
    ENDIF
    IF(L.EQ.NU(NY,NX))THEN
      DPTHZ(L,NY,NX)=0.5_r8*CDPTHZ(L,NY,NX)
    !     DDLYRX(NN)=DDLYRX(NN)+DDLYR(L,5)
    ELSE
      DPTHZ(L,NY,NX)=0.5_r8*(CDPTHZ(L,NY,NX)+CDPTHZ(L-1,NY,NX))
    ENDIF
    IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
    !     DDLYRX(NN)=CumDepth2LayerBottom(L,NY,NX)-CDPTHX(L,NY,NX)
      DDLYRX(NN)=CDPTHY(L,NY,NX)-CDPTHX(L,NY,NX)
    ENDIF
!
  !     RESET POND SURFACE LAYER NUMBER IF LOST TO EVAPORATION
      !
  ELSEIF(NN.EQ.2)THEN
    IF((L.EQ.NU(NY,NX).AND.SoiBulkDensity(NU(NY,NX),NY,NX).LE.ZERO) &
      .AND.(VHeatCapacity(NU(NY,NX),NY,NX).LE.VHCPNX(NY,NX) &
      .OR.NUM(NY,NX).GT.NU(NY,NX)))THEN
      NUX=NU(NY,NX)
      DO LL=NUX+1,NL(NY,NX)
        IF(VLSoilPoreMicP(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
          NU(NY,NX)=LL
          DDLYRX(NN)=DLYR(3,NUX,NY,NX)
          IFLGL(L,NN)=1
          DLYR(3,NUX,NY,NX)=0.0_r8
          IF(SoiBulkDensity(NUX,NY,NX).LE.ZERO)THEN
            VGeomLayer(NUX,NY,NX)=AREA(3,NUX,NY,NX)*DLYR(3,NUX,NY,NX)
            VLSoilPoreMicP(NUX,NY,NX)=VGeomLayer(NUX,NY,NX)*FMPR(NUX,NY,NX)
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
    XVOLWP=AZMAX1(VLWatMicP(0,NY,NX)-VOLWD(NY,NX))
    IF(L.EQ.NU(NY,NX).AND.CumDepth2LayerBottom(0,NY,NX).GT.CumSoilDeptht0(NY,NX) &
      .AND.XVOLWP.GT.VOLWD(NY,NX)+VHCPNX(NY,NX)/cpw)THEN
          !     IF((SoiBulkDensity(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))
          !    2.OR.(SoiBulkDensity(L,NY,NX).LE.ZERO))THEN
      IF(SoiBulkDensity(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))THEN
        NU(NY,NX)=NUI(NY,NX)
        NUM(NY,NX)=NUI(NY,NX)
        DDLYRX(NN)=(VOLWD(NY,NX)-XVOLWP)/AREA(3,0,NY,NX)
        IFLGL(L,NN)=1
        DLYR0=(AZMAX1(VLWatMicP(0,NY,NX)+VLiceMicP(0,NY,NX)-VWatLitrX(NY,NX)) &
          +VLitR(NY,NX))/AREA(3,0,NY,NX)
        DLYR(3,0,NY,NX)=DLYR0+DDLYRX(NN)
        DLYR(3,NU(NY,NX),NY,NX)=DLYR(3,NU(NY,NX),NY,NX)-DDLYRX(NN)
        IF(L.GT.2)THEN
          DO LL=L-2,NU(NY,NX),-1
            CumDepth2LayerBottom(LL,NY,NX)=CumDepth2LayerBottom(L-1,NY,NX)
            CDPTHY(LL,NY,NX)=CDPTHY(L-1,NY,NX)
          ENDDO
        ENDIF
        CumDepth2LayerBottom(0,NY,NX)=CumDepth2LayerBottom(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
        CDPTHY(0,NY,NX)=CDPTHY(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
        DPTH(NU(NY,NX),NY,NX)=0.5_r8*(CumDepth2LayerBottom(NU(NY,NX),NY,NX)+CumDepth2LayerBottom(0,NY,NX))
        CDPTHZ(NU(NY,NX),NY,NX)=DLYR(3,NU(NY,NX),NY,NX)
        DPTHZ(NU(NY,NX),NY,NX)=0.5_r8*CDPTHZ(NU(NY,NX),NY,NX)
      ELSE
        DDLYRX(NN)=0.0_r8
        IFLGL(L,NN)=0
      ENDIF
    ELSE
      DDLYRX(NN)=0.0_r8
      IFLGL(L,NN)=0
    ENDIF
  ENDIF
  end subroutine getTSlyrthicks
!------------------------------------------------------------------------------------------

  subroutine tgtPondLyr(L,L0,L1,NY,NX,NN,FX,FY,CDPTHY,IFLGL)
  implicit none
  integer, intent(in) :: L
  integer, intent(in) :: L0  !target
  integer, intent(in) :: L1  !source
  integer, intent(in) :: NY,NX,NN
  real(r8), intent(in) :: FX,FY
  real(r8), intent(inout) :: CDPTHY(0:JZ,JY,JX)
  integer, intent(in) ::  IFLGL(0:JZ,6)
  integer :: N,M,NZ,K,NGL,NR,NE,NTU,NTSA,NTSAB,NTG,NTP
  integer :: NTX,NTF
  real(r8) :: ENGY0,ENGY1
! begin_execution

  IF(L0.NE.0)THEN
    SAND(L1,NY,NX)=SAND(L1,NY,NX)+FX*SAND(L0,NY,NX)
    SILT(L1,NY,NX)=SILT(L1,NY,NX)+FX*SILT(L0,NY,NX)
    CLAY(L1,NY,NX)=CLAY(L1,NY,NX)+FX*CLAY(L0,NY,NX)
    trcx_solml(idx_CEC,L1,NY,NX)=trcx_solml(idx_CEC,L1,NY,NX) &
      +FX*trcx_solml(idx_CEC,L0,NY,NX)
    trcx_solml(idx_AEC,L1,NY,NX)=trcx_solml(idx_AEC,L1,NY,NX) &
      +FX*trcx_solml(idx_AEC,L0,NY,NX)
  ENDIF

  VLWatMicP(L1,NY,NX)=VLWatMicP(L1,NY,NX)+FX*VLWatMicP(L0,NY,NX)
  VLiceMicP(L1,NY,NX)=VLiceMicP(L1,NY,NX)+FX*VLiceMicP(L0,NY,NX)
  VLsoiAirP(L1,NY,NX)=VLsoiAirP(L1,NY,NX)+FX*VLsoiAirP(L0,NY,NX)
  VLMicP(L1,NY,NX)=VLMicP(L1,NY,NX)+FX*VLMicP(L0,NY,NX)
  VLSoilMicP(L1,NY,NX)=VLSoilMicP(L1,NY,NX)+FX*VLSoilMicP(L0,NY,NX)
  VLWatMicPX(L1,NY,NX)=VLWatMicP(L1,NY,NX)
  ENGY1=VHeatCapacity(L1,NY,NX)*TKS(L1,NY,NX)
  ENGY0=VHeatCapacity(L0,NY,NX)*TKS(L0,NY,NX)
  ENGY1=ENGY1+FX*ENGY0
  VHeatCapacitySoilM(L1,NY,NX)=VHeatCapacitySoilM(L1,NY,NX)+FX*VHeatCapacitySoilM(L0,NY,NX)
  VHeatCapacity(L1,NY,NX)=VHeatCapacitySoilM(L1,NY,NX) &
    +cpw*(VLWatMicP(L1,NY,NX)+VLWatMacP(L1,NY,NX)) &
    +cpi*(VLiceMicP(L1,NY,NX)+VLiceMacP(L1,NY,NX))

  IF(VHeatCapacity(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L1,NY,NX)=ENGY1/VHeatCapacity(L1,NY,NX)
  ELSE
    TKS(L1,NY,NX)=TKS(L0,NY,NX)
  ENDIF
  TCS(L1,NY,NX)=units%Kelvin2Celcius(TKS(L1,NY,NX))

  DO NTF=ifertn_beg,ifertn_end
    FertN_soil(NTF,L1,NY,NX)=FertN_soil(NTF,L1,NY,NX)+FX*FertN_soil(NTF,L0,NY,NX)
  ENDDO

  DO NTF=ifertnb_beg,ifertnb_end
    FertN_band(NTF,L1,NY,NX)=FertN_band(NTF,L1,NY,NX)+FX*FertN_band(NTF,L0,NY,NX)
  ENDDO

  DO NTU=ids_nuts_beg,ids_nuts_end
    if(NTU/=ids_H2PO4B .and. NTU/=ids_H1PO4B)THEN
      trc_solml(NTU,L1,NY,NX)=trc_solml(NTU,L1,NY,NX)+FX*trc_solml(NTU,L0,NY,NX)
    ENDIF
  ENDDO

  IF(salt_model)THEN
    DO NTSA=idsa_beg,idsa_end
      trcsa_solml(NTSA,L1,NY,NX)=trcsa_solml(NTSA,L1,NY,NX)&
        +FX*trcsa_solml(NTSA,L0,NY,NX)
    ENDDO
  ENDIF

  IF(L0.NE.0)THEN
    trc_solml(ids_H1PO4B,L1,NY,NX)=trc_solml(ids_H1PO4B,L1,NY,NX) &
      +FX*trc_solml(ids_H1PO4B,L0,NY,NX)
    trc_solml(ids_H2PO4B,L1,NY,NX)=trc_solml(ids_H2PO4B,L1,NY,NX) &
      +FX*trc_solml(ids_H2PO4B,L0,NY,NX)

    IF(salt_model)THEN
      DO NTSAB=idsab_beg,idsab_end
        trcsa_solml(NTSAB,L1,NY,NX)=trcsa_solml(NTSAB,L1,NY,NX)+FX*trcsa_solml(NTSAB,L0,NY,NX)
      ENDDO
    ENDIF

    DO NTX=idx_beg+1,idx_cation_end
      trcx_solml(NTX,L1,NY,NX)=trcx_solml(NTX,L1,NY,NX)+FX*trcx_solml(NTX,L0,NY,NX)
    ENDDO

    DO NTX=idx_AEC+1,idx_end
      trcx_solml(NTX,L1,NY,NX)=trcx_solml(NTX,L1,NY,NX)+FX*trcx_solml(NTX,L0,NY,NX)
    ENDDO

    DO NTP=idsp_beg,idsp_end
      trcp_salml(NTP,L1,NY,NX)=trcp_salml(NTP,L1,NY,NX)+FX*trcp_salml(NTP,L0,NY,NX)
    ENDDO

    DO NTG=idg_beg,idg_end-1
      trc_gasml(NTG,L1,NY,NX)=trc_gasml(NTG,L1,NY,NX)+FX*trc_gasml(NTG,L0,NY,NX)
    ENDDO
  ENDIF
!exclude NH3 and NH3B, which are accounted in nutrients
  DO NTG=idg_beg,idg_end-2
    trc_solml(NTG,L1,NY,NX)=trc_solml(NTG,L1,NY,NX)+FX*trc_solml(NTG,L0,NY,NX)
  ENDDO

  IF(IFLGL(L,3).EQ.0)THEN
    DO  K=1,jcplx
      DO  N=1,NFGs
        DO  M=1,nlbiomcp
          DO NGL=JGnio(N),JGnfo(N)
            OMC(M,NGL,K,L1,NY,NX)=OMC(M,NGL,K,L1,NY,NX)+FX*OMC(M,NGL,K,L0,NY,NX)
            OMN(M,NGL,K,L1,NY,NX)=OMN(M,NGL,K,L1,NY,NX)+FX*OMN(M,NGL,K,L0,NY,NX)
            OMP(M,NGL,K,L1,NY,NX)=OMP(M,NGL,K,L1,NY,NX)+FX*OMP(M,NGL,K,L0,NY,NX)
          enddo
        enddo
      enddo
    ENDDO
    DO  N=1,NFGs
      DO  M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          OMCff(M,NGL,L1,NY,NX)=OMCff(M,NGL,L1,NY,NX)+FX*OMCff(M,NGL,L0,NY,NX)
          OMNff(M,NGL,L1,NY,NX)=OMNff(M,NGL,L1,NY,NX)+FX*OMNff(M,NGL,L0,NY,NX)
          OMPff(M,NGL,L1,NY,NX)=OMPff(M,NGL,L1,NY,NX)+FX*OMPff(M,NGL,L0,NY,NX)
        enddo
      enddo
    enddo

    DO K=1,jcplx
      DO  M=1,ndbiomcp
        ORC(M,K,L1,NY,NX)=ORC(M,K,L1,NY,NX)+FX*ORC(M,K,L0,NY,NX)
        ORN(M,K,L1,NY,NX)=ORN(M,K,L1,NY,NX)+FX*ORN(M,K,L0,NY,NX)
        ORP(M,K,L1,NY,NX)=ORP(M,K,L1,NY,NX)+FX*ORP(M,K,L0,NY,NX)
      ENDDO
      OQC(K,L1,NY,NX)=OQC(K,L1,NY,NX)+FX*OQC(K,L0,NY,NX)
      OQN(K,L1,NY,NX)=OQN(K,L1,NY,NX)+FX*OQN(K,L0,NY,NX)
      OQP(K,L1,NY,NX)=OQP(K,L1,NY,NX)+FX*OQP(K,L0,NY,NX)
      OQA(K,L1,NY,NX)=OQA(K,L1,NY,NX)+FX*OQA(K,L0,NY,NX)
      OQCH(K,L1,NY,NX)=OQCH(K,L1,NY,NX)+FX*OQCH(K,L0,NY,NX)
      OQNH(K,L1,NY,NX)=OQNH(K,L1,NY,NX)+FX*OQNH(K,L0,NY,NX)
      OQPH(K,L1,NY,NX)=OQPH(K,L1,NY,NX)+FX*OQPH(K,L0,NY,NX)
      OQAH(K,L1,NY,NX)=OQAH(K,L1,NY,NX)+FX*OQAH(K,L0,NY,NX)
      OHC(K,L1,NY,NX)=OHC(K,L1,NY,NX)+FX*OHC(K,L0,NY,NX)
      OHN(K,L1,NY,NX)=OHN(K,L1,NY,NX)+FX*OHN(K,L0,NY,NX)
      OHP(K,L1,NY,NX)=OHP(K,L1,NY,NX)+FX*OHP(K,L0,NY,NX)
      OHA(K,L1,NY,NX)=OHA(K,L1,NY,NX)+FX*OHA(K,L0,NY,NX)
      DO M=1,jsken
        OSC(M,K,L1,NY,NX)=OSC(M,K,L1,NY,NX)+FX*OSC(M,K,L0,NY,NX)
        OSA(M,K,L1,NY,NX)=OSA(M,K,L1,NY,NX)+FX*OSA(M,K,L0,NY,NX)
        OSN(M,K,L1,NY,NX)=OSN(M,K,L1,NY,NX)+FX*OSN(M,K,L0,NY,NX)
        OSP(M,K,L1,NY,NX)=OSP(M,K,L1,NY,NX)+FX*OSP(M,K,L0,NY,NX)
      ENDDO
    ENDDO
  ENDIF
          !
!     TARGET ROOT LAYER
!
  IF(L0.NE.0)THEN
    DO  NZ=1,NP(NY,NX)
      IF(WTRTL(ipltroot,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.WTRTL(ipltroot,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        DO N=1,MY(NZ,NY,NX)
          DO NTG=idg_beg,idg_end-1
            trcg_rootml(NTG,N,L1,NZ,NY,NX)=trcg_rootml(NTG,N,L1,NZ,NY,NX)+FX*trcg_rootml(NTG,N,L0,NZ,NY,NX)
            trcs_rootml(NTG,N,L1,NZ,NY,NX)=trcs_rootml(NTG,N,L1,NZ,NY,NX)+FX*trcs_rootml(NTG,N,L0,NZ,NY,NX)
          ENDDO
          DO  NR=1,NRT(NZ,NY,NX)
            DO NE=1,npelms
              WTRT1E(NE,N,L1,NR,NZ,NY,NX)=WTRT1E(NE,N,L1,NR,NZ,NY,NX)+FX*WTRT1E(NE,N,L0,NR,NZ,NY,NX)
              WTRT2E(NE,N,L1,NR,NZ,NY,NX)=WTRT2E(NE,N,L1,NR,NZ,NY,NX)+FX*WTRT2E(NE,N,L0,NR,NZ,NY,NX)
            ENDDO
            RTLG1(N,L1,NR,NZ,NY,NX)=RTLG1(N,L1,NR,NZ,NY,NX)+FX*RTLG1(N,L0,NR,NZ,NY,NX)
            RTLG2(N,L1,NR,NZ,NY,NX)=RTLG2(N,L1,NR,NZ,NY,NX)+FX*RTLG2(N,L0,NR,NZ,NY,NX)
            RTN2(N,L1,NR,NZ,NY,NX)=RTN2(N,L1,NR,NZ,NY,NX)+FX*RTN2(N,L0,NR,NZ,NY,NX)
          ENDDO
          DO NE=1,npelms
            EPOOLR(NE,N,L1,NZ,NY,NX)=EPOOLR(NE,N,L1,NZ,NY,NX)+FX*EPOOLR(NE,N,L0,NZ,NY,NX)
          ENDDO
          WTRTL(N,L1,NZ,NY,NX)=WTRTL(N,L1,NZ,NY,NX)+FX*WTRTL(N,L0,NZ,NY,NX)
          WTRTD(N,L1,NZ,NY,NX)=WTRTD(N,L1,NZ,NY,NX)+FX*WTRTD(N,L0,NZ,NY,NX)
          WSRTL(N,L1,NZ,NY,NX)=WSRTL(N,L1,NZ,NY,NX)+FX*WSRTL(N,L0,NZ,NY,NX)
          RTN1(N,L1,NZ,NY,NX)=RTN1(N,L1,NZ,NY,NX)+FX*RTN1(N,L0,NZ,NY,NX)
          RTNL(N,L1,NZ,NY,NX)=RTNL(N,L1,NZ,NY,NX)+FX*RTNL(N,L0,NZ,NY,NX)
          RTLGP(N,L1,NZ,NY,NX)=RTLGP(N,L1,NZ,NY,NX)+FX*RTLGP(N,L0,NZ,NY,NX)
          RTDNP(N,L1,NZ,NY,NX)=RTDNP(N,L1,NZ,NY,NX)+FX*RTDNP(N,L0,NZ,NY,NX)
          RTVLP(N,L1,NZ,NY,NX)=RTVLP(N,L1,NZ,NY,NX)+FX*RTVLP(N,L0,NZ,NY,NX)
          RTVLW(N,L1,NZ,NY,NX)=RTVLW(N,L1,NZ,NY,NX)+FX*RTVLW(N,L0,NZ,NY,NX)
          RRAD1(N,L1,NZ,NY,NX)=RRAD1(N,L1,NZ,NY,NX)+FX*RRAD1(N,L0,NZ,NY,NX)
          RRAD2(N,L1,NZ,NY,NX)=RRAD2(N,L1,NZ,NY,NX)+FX*RRAD2(N,L0,NZ,NY,NX)
          RTARP(N,L1,NZ,NY,NX)=RTARP(N,L1,NZ,NY,NX)+FX*RTARP(N,L0,NZ,NY,NX)
          RTLGA(N,L1,NZ,NY,NX)=RTLGA(N,L1,NZ,NY,NX)+FX*RTLGA(N,L0,NZ,NY,NX)
        ENDDO
        DO NE=1,npelms
          WTNDLE(NE,L1,NZ,NY,NX)=WTNDLE(NE,L1,NZ,NY,NX)+FX*WTNDLE(NE,L0,NZ,NY,NX)
          EPOOLN(NE,L1,NZ,NY,NX)=EPOOLN(NE,L1,NZ,NY,NX)+FX*EPOOLN(NE,L0,NZ,NY,NX)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!
          !     SOURCE POND LAYER
    !
  IF(L0.NE.0)THEN
    SAND(L0,NY,NX)=FY*SAND(L0,NY,NX)
    SILT(L0,NY,NX)=FY*SILT(L0,NY,NX)
    CLAY(L0,NY,NX)=FY*CLAY(L0,NY,NX)
    trcx_solml(idx_CEC,L0,NY,NX)=FY*trcx_solml(idx_CEC,L0,NY,NX)
    trcx_solml(idx_AEC,L0,NY,NX)=FY*trcx_solml(idx_AEC,L0,NY,NX)
  ENDIF
!     IF(SoiBulkDensity(L0,NY,NX).LE.ZERO)THEN
!     VGeomLayer(L0,NY,NX)=FY*VGeomLayer(L0,NY,NX)
!     VLSoilPoreMicP(L0,NY,NX)=FY*VLSoilPoreMicP(L0,NY,NX)
!     ENDIF
  VLWatMicP(L0,NY,NX)=FY*VLWatMicP(L0,NY,NX)
  VLiceMicP(L0,NY,NX)=FY*VLiceMicP(L0,NY,NX)
  VLsoiAirP(L0,NY,NX)=FY*VLsoiAirP(L0,NY,NX)
  VLMicP(L0,NY,NX)=FY*VLMicP(L0,NY,NX)
  VLSoilMicP(L0,NY,NX)=FY*VLSoilMicP(L0,NY,NX)
  VLWatMicPX(L0,NY,NX)=VLWatMicP(L0,NY,NX)
  ENGY0=FY*ENGY0
  VHeatCapacitySoilM(L0,NY,NX)=FY*VHeatCapacitySoilM(L0,NY,NX)
  IF(L0.NE.0)THEN
    VHeatCapacity(L0,NY,NX)=VHeatCapacitySoilM(L0,NY,NX) &
      +cpw*(VLWatMicP(L0,NY,NX)+VLWatMacP(L0,NY,NX)) &
      +cpi*(VLiceMicP(L0,NY,NX)+VLiceMacP(L0,NY,NX))
  ELSE
    VHeatCapacity(L0,NY,NX)=VHeatCapacitySoilM(L0,NY,NX)+cpw*VLWatMicP(L0,NY,NX)+cpi*VLiceMicP(L0,NY,NX)
  ENDIF
  IF(VHeatCapacity(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L0,NY,NX)=ENGY0/VHeatCapacity(L0,NY,NX)
  ELSE
    TKS(L0,NY,NX)=TKS(L1,NY,NX)
  ENDIF
  TCS(L0,NY,NX)=units%Kelvin2Celcius(TKS(L0,NY,NX))

  DO NTF=ifertn_beg,ifertn_end
    FertN_soil(NTF,L0,NY,NX)=FY*FertN_soil(NTF,L0,NY,NX)
  ENDDO

  DO NTF=ifertnb_beg,ifertnb_end
    FertN_band(NTF,L0,NY,NX)=FY*FertN_band(NTF,L0,NY,NX)
  ENDDO

  DO NTU=ids_nuts_beg,ids_nuts_end
    if(NTU/=ids_H1PO4B .and. NTU/=ids_H2PO4B)THEN
      trc_solml(NTU,L0,NY,NX)=FY*trc_solml(NTU,L0,NY,NX)
    ENDIF
  ENDDO
  IF(salt_model)THEN
    DO NTSA=idsa_beg,idsa_end
      trcsa_solml(NTSA,L0,NY,NX)=FY*trcsa_solml(NTSA,L0,NY,NX)
    ENDDO
  ENDIF
  IF(L0.NE.0)THEN
    trc_solml(ids_H1PO4B,L0,NY,NX)=FY*trc_solml(ids_H1PO4B,L0,NY,NX)
    trc_solml(ids_H2PO4B,L0,NY,NX)=FY*trc_solml(ids_H2PO4B,L0,NY,NX)
    IF(salt_model)THEN
      DO NTSAB=idsab_beg,idsab_end
        trcsa_solml(NTSAB,L0,NY,NX)=FY*trcsa_solml(NTSAB,L0,NY,NX)
      ENDDO
    ENDIF

    DO NTX=idx_beg+1,idx_cation_end
      trcx_solml(NTX,L0,NY,NX)=FY*trcx_solml(NTX,L0,NY,NX)
    ENDDO
    DO NTX=idx_AEC+1,idx_end
      trcx_solml(NTX,L0,NY,NX)=FY*trcx_solml(NTX,L0,NY,NX)
    ENDDO

    DO NTP=idsp_beg,idsp_end
      trcp_salml(NTP,L0,NY,NX)=FY*trcp_salml(NTP,L0,NY,NX)
    ENDDO


    DO NTG=idg_beg,idg_end-1
      trc_gasml(NTG,L0,NY,NX)=FY*trc_gasml(NTG,L0,NY,NX)
    ENDDO
  ENDIF
!exclude NH3 and NH3B, which are accounted in nutrients
  DO NTG=idg_beg,idg_end-2
    trc_solml(NTG,L0,NY,NX)=FY*trc_solml(NTG,L0,NY,NX)
  ENDDO
  IF(IFLGL(L,3).EQ.0)THEN
    DO  K=1,jcplx
       DO N=1,NFGs
        DO M=1,nlbiomcp
          DO NGL=JGnio(N),JGnfo(N)
            OMC(M,NGL,K,L0,NY,NX)=FY*OMC(M,NGL,K,L0,NY,NX)
            OMN(M,NGL,K,L0,NY,NX)=FY*OMN(M,NGL,K,L0,NY,NX)
            OMP(M,NGL,K,L0,NY,NX)=FY*OMP(M,NGL,K,L0,NY,NX)
          ENDDO
        enddo
      enddo
    ENDDO

    DO N=1,NFGs
      DO M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          OMCff(M,NGL,L0,NY,NX)=FY*OMCff(M,NGL,L0,NY,NX)
          OMNff(M,NGL,L0,NY,NX)=FY*OMNff(M,NGL,L0,NY,NX)
          OMPff(M,NGL,L0,NY,NX)=FY*OMPff(M,NGL,L0,NY,NX)
        ENDDO
      enddo
    enddo

    DO K=1,jcplx
      DO  M=1,ndbiomcp
        ORC(M,K,L0,NY,NX)=FY*ORC(M,K,L0,NY,NX)
        ORN(M,K,L0,NY,NX)=FY*ORN(M,K,L0,NY,NX)
        ORP(M,K,L0,NY,NX)=FY*ORP(M,K,L0,NY,NX)
      ENDDO
      OQC(K,L0,NY,NX)=FY*OQC(K,L0,NY,NX)
      OQN(K,L0,NY,NX)=FY*OQN(K,L0,NY,NX)
      OQP(K,L0,NY,NX)=FY*OQP(K,L0,NY,NX)
      OQA(K,L0,NY,NX)=FY*OQA(K,L0,NY,NX)
      OQCH(K,L0,NY,NX)=FY*OQCH(K,L0,NY,NX)
      OQNH(K,L0,NY,NX)=FY*OQNH(K,L0,NY,NX)
      OQPH(K,L0,NY,NX)=FY*OQPH(K,L0,NY,NX)
      OQAH(K,L0,NY,NX)=FY*OQAH(K,L0,NY,NX)
      OHC(K,L0,NY,NX)=FY*OHC(K,L0,NY,NX)
      OHN(K,L0,NY,NX)=FY*OHN(K,L0,NY,NX)
      OHP(K,L0,NY,NX)=FY*OHP(K,L0,NY,NX)
      OHA(K,L0,NY,NX)=FY*OHA(K,L0,NY,NX)
      DO  M=1,jsken
        OSC(M,K,L0,NY,NX)=FY*OSC(M,K,L0,NY,NX)
        OSA(M,K,L0,NY,NX)=FY*OSA(M,K,L0,NY,NX)
        OSN(M,K,L0,NY,NX)=FY*OSN(M,K,L0,NY,NX)
        OSP(M,K,L0,NY,NX)=FY*OSP(M,K,L0,NY,NX)
      ENDDO
    ENDDO
  ENDIF
          !
          !     SOURCE ROOT LAYER
          !
  IF(L0.NE.0)THEN
    DO  NZ=1,NP(NY,NX)
      IF(WTRTL(ipltroot,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.WTRTL(ipltroot,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        DO  N=1,MY(NZ,NY,NX)
          DO NTG=idg_beg,idg_end-1
            trcg_rootml(NTG,N,L0,NZ,NY,NX)=FY*trcg_rootml(NTG,N,L0,NZ,NY,NX)
            trcs_rootml(NTG,N,L0,NZ,NY,NX)=FY*trcs_rootml(NTG,N,L0,NZ,NY,NX)
          ENDDO
          DO NR=1,NRT(NZ,NY,NX)
            DO NE=1,npelms
              WTRT1E(NE,N,L0,NR,NZ,NY,NX)=FY*WTRT1E(NE,N,L0,NR,NZ,NY,NX)
              WTRT2E(NE,N,L0,NR,NZ,NY,NX)=FY*WTRT2E(NE,N,L0,NR,NZ,NY,NX)
            ENDDO
            RTLG1(N,L0,NR,NZ,NY,NX)=FY*RTLG1(N,L0,NR,NZ,NY,NX)
            RTLG2(N,L0,NR,NZ,NY,NX)=FY*RTLG2(N,L0,NR,NZ,NY,NX)
            RTN2(N,L0,NR,NZ,NY,NX)=FY*RTN2(N,L0,NR,NZ,NY,NX)
          ENDDO
          DO NE=1,npelms
            EPOOLR(NE,N,L0,NZ,NY,NX)=FY*EPOOLR(NE,N,L0,NZ,NY,NX)
          ENDDO
          WTRTL(N,L0,NZ,NY,NX)=FY*WTRTL(N,L0,NZ,NY,NX)
          WTRTD(N,L0,NZ,NY,NX)=FY*WTRTD(N,L0,NZ,NY,NX)
          WSRTL(N,L0,NZ,NY,NX)=FY*WSRTL(N,L0,NZ,NY,NX)
          RTN1(N,L0,NZ,NY,NX)=FY*RTN1(N,L0,NZ,NY,NX)
          RTNL(N,L0,NZ,NY,NX)=FY*RTNL(N,L0,NZ,NY,NX)
          RTLGP(N,L0,NZ,NY,NX)=FY*RTLGP(N,L0,NZ,NY,NX)
          RTDNP(N,L0,NZ,NY,NX)=FY*RTDNP(N,L0,NZ,NY,NX)
          RTVLP(N,L0,NZ,NY,NX)=FY*RTVLP(N,L0,NZ,NY,NX)
          RTVLW(N,L0,NZ,NY,NX)=FY*RTVLW(N,L0,NZ,NY,NX)
          RRAD1(N,L0,NZ,NY,NX)=FY*RRAD1(N,L0,NZ,NY,NX)
          RRAD2(N,L0,NZ,NY,NX)=FY*RRAD2(N,L0,NZ,NY,NX)
          RTARP(N,L0,NZ,NY,NX)=FY*RTARP(N,L0,NZ,NY,NX)
          RTLGA(N,L0,NZ,NY,NX)=FY*RTLGA(N,L0,NZ,NY,NX)
        ENDDO
        DO NE=1,npelms
          WTNDLE(NE,L0,NZ,NY,NX)=FY*WTNDLE(NE,L0,NZ,NY,NX)

          EPOOLN(NE,L0,NZ,NY,NX)=FY*EPOOLN(NE,L0,NZ,NY,NX)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  IF(NN.EQ.1)THEN
    IF(SoiBulkDensity(L0,NY,NX).LE.ZERO.AND.SoiBulkDensity(L1,NY,NX).LE.ZERO &
      .AND.VLWatMicP(L0,NY,NX)+VLiceMicP(L0,NY,NX).LE.ZEROS(NY,NX))THEN
      CumDepth2LayerBottom(L1,NY,NX)=CumDepth2LayerBottom(L0,NY,NX)
      CDPTHY(L1,NY,NX)=CDPTHY(L0,NY,NX)
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

  integer :: K,N,M,NGL,NR,NZ,NE,NTG
  real(r8) :: FXO,FRO
  real(r8) :: FXRTLG2,FXRTN2,FXEPOOLR,FXWTRTL
  real(r8) :: FXWTNDLE,FXEPOOLN
  real(r8) :: FXWTRT1E
  real(r8) :: FXWTRT2E,FXRTLG1
  real(r8) :: FXWTRTD,FXWSRTL,FXRTN1,FXRTNL,FXRTLGP,FXRTDNP
  real(r8) :: FXRTVLP,FXRTVLW,FXRRAD1,FXRRAD2,FXRTARP,FXRTLGA
  real(r8) :: FXOQN,FXOQP,FXOQA,FXOQCH,FXOQNH,FXOQPH,FXOQAH
  real(r8) :: FXOHC,FXOHN,FXOHP,FXOHA,FXOSC,FXOSA,FXOSN,FXOSP
  real(r8) :: FXOMC,FXOMN,FXOMP,FXORC,FXORN,FXORP,FXOQC
  real(r8) :: FXGA,FXGP
! begin_execution
  IF(IFLGL(L,3).EQ.0.AND.L0.NE.0 &
    .AND.VLSoilPoreMicP(L0,NY,NX).GT.ZEROS(NY,NX) &
    .AND.VLSoilPoreMicP(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    IF(L0.EQ.L.OR.CORGCI(L0,NY,NX).LE.ZERO)THEN
      FXO=FO
    ELSE
      FXO=AMIN1(0.5,FO*AMIN1(10.0,CORGCI(L1,NY,NX) &
        /CORGCI(L0,NY,NX)))
    ENDIF

    DO  K=1,jcplx
      DO  N=1,NFGs
        DO  M=1,nlbiomcp
          DO NGL=JGnio(N),JGnfo(N)
            FXOMC=FXO*OMC(M,NGL,K,L0,NY,NX)
            OMC(M,NGL,K,L1,NY,NX)=OMC(M,NGL,K,L1,NY,NX)+FXOMC
            OMC(M,NGL,K,L0,NY,NX)=OMC(M,NGL,K,L0,NY,NX)-FXOMC
            FXOMN=FXO*OMN(M,NGL,K,L0,NY,NX)
            OMN(M,NGL,K,L1,NY,NX)=OMN(M,NGL,K,L1,NY,NX)+FXOMN
            OMN(M,NGL,K,L0,NY,NX)=OMN(M,NGL,K,L0,NY,NX)-FXOMN
            FXOMP=FXO*OMP(M,NGL,K,L0,NY,NX)
            OMP(M,NGL,K,L1,NY,NX)=OMP(M,NGL,K,L1,NY,NX)+FXOMP
            OMP(M,NGL,K,L0,NY,NX)=OMP(M,NGL,K,L0,NY,NX)-FXOMP
          enddo
        enddo
      enddo
    ENDDO

    DO  N=1,NFGs
      DO  M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          FXOMC=FXO*OMCff(M,NGL,L0,NY,NX)
          OMCff(M,NGL,L1,NY,NX)=OMCff(M,NGL,L1,NY,NX)+FXOMC
          OMCff(M,NGL,L0,NY,NX)=OMCff(M,NGL,L0,NY,NX)-FXOMC
          FXOMN=FXO*OMNff(M,NGL,L0,NY,NX)
          OMNff(M,NGL,L1,NY,NX)=OMNff(M,NGL,L1,NY,NX)+FXOMN
          OMNff(M,NGL,L0,NY,NX)=OMNff(M,NGL,L0,NY,NX)-FXOMN
          FXOMP=FXO*OMPff(M,NGL,L0,NY,NX)
          OMPff(M,NGL,L1,NY,NX)=OMPff(M,NGL,L1,NY,NX)+FXOMP
          OMPff(M,NGL,L0,NY,NX)=OMPff(M,NGL,L0,NY,NX)-FXOMP
        enddo
      enddo
    enddo

    DO  K=1,jcplx
      DO  M=1,ndbiomcp
        FXORC=FXO*ORC(M,K,L0,NY,NX)
        ORC(M,K,L1,NY,NX)=ORC(M,K,L1,NY,NX)+FXORC
        ORC(M,K,L0,NY,NX)=ORC(M,K,L0,NY,NX)-FXORC
        FXORN=FXO*ORN(M,K,L0,NY,NX)
        ORN(M,K,L1,NY,NX)=ORN(M,K,L1,NY,NX)+FXORN
        ORN(M,K,L0,NY,NX)=ORN(M,K,L0,NY,NX)-FXORN
        FXORP=FXO*ORP(M,K,L0,NY,NX)
        ORP(M,K,L1,NY,NX)=ORP(M,K,L1,NY,NX)+FXORP
        ORP(M,K,L0,NY,NX)=ORP(M,K,L0,NY,NX)-FXORP
      ENDDO
      FXOQC=FXO*OQC(K,L0,NY,NX)
      OQC(K,L1,NY,NX)=OQC(K,L1,NY,NX)+FXOQC
      OQC(K,L0,NY,NX)=OQC(K,L0,NY,NX)-FXOQC
      FXOQN=FXO*OQN(K,L0,NY,NX)
      OQN(K,L1,NY,NX)=OQN(K,L1,NY,NX)+FXOQN
      OQN(K,L0,NY,NX)=OQN(K,L0,NY,NX)-FXOQN
      FXOQP=FXO*OQP(K,L0,NY,NX)
      OQP(K,L1,NY,NX)=OQP(K,L1,NY,NX)+FXOQP
      OQP(K,L0,NY,NX)=OQP(K,L0,NY,NX)-FXOQP
      FXOQA=FXO*OQA(K,L0,NY,NX)
      OQA(K,L1,NY,NX)=OQA(K,L1,NY,NX)+FXOQA
      OQA(K,L0,NY,NX)=OQA(K,L0,NY,NX)-FXOQA
      IF(SoilFracAsMacP(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP(L0,NY,NX).GT.ZERO)THEN
        FXOQCH=FXO*OQCH(K,L0,NY,NX)
        OQCH(K,L1,NY,NX)=OQCH(K,L1,NY,NX)+FXOQCH
        OQCH(K,L0,NY,NX)=OQCH(K,L0,NY,NX)-FXOQCH
        FXOQNH=FXO*OQNH(K,L0,NY,NX)
        OQNH(K,L1,NY,NX)=OQNH(K,L1,NY,NX)+FXOQNH
        OQNH(K,L0,NY,NX)=OQNH(K,L0,NY,NX)-FXOQNH
        FXOQPH=FXO*OQPH(K,L0,NY,NX)
        OQPH(K,L1,NY,NX)=OQPH(K,L1,NY,NX)+FXOQPH
        OQPH(K,L0,NY,NX)=OQPH(K,L0,NY,NX)-FXOQPH
        FXOQAH=FXO*OQAH(K,L0,NY,NX)
        OQAH(K,L1,NY,NX)=OQAH(K,L1,NY,NX)+FXOQAH
        OQAH(K,L0,NY,NX)=OQAH(K,L0,NY,NX)-FXOQAH
      ENDIF
      FXOHC=FXO*OHC(K,L0,NY,NX)
      OHC(K,L1,NY,NX)=OHC(K,L1,NY,NX)+FXOHC
      OHC(K,L0,NY,NX)=OHC(K,L0,NY,NX)-FXOHC
      FXOHN=FXO*OHN(K,L0,NY,NX)
      OHN(K,L1,NY,NX)=OHN(K,L1,NY,NX)+FXOHN
      OHN(K,L0,NY,NX)=OHN(K,L0,NY,NX)-FXOHN
      FXOHP=FXO*OHP(K,L0,NY,NX)
      OHP(K,L1,NY,NX)=OHP(K,L1,NY,NX)+FXOHP
      OHP(K,L0,NY,NX)=OHP(K,L0,NY,NX)-FXOHP
      FXOHA=FXO*OHA(K,L0,NY,NX)
      OHA(K,L1,NY,NX)=OHA(K,L1,NY,NX)+FXOHA
      OHA(K,L0,NY,NX)=OHA(K,L0,NY,NX)-FXOHA
      DO M=1,jsken
        FXOSC=FXO*OSC(M,K,L0,NY,NX)
        OSC(M,K,L1,NY,NX)=OSC(M,K,L1,NY,NX)+FXOSC
        OSC(M,K,L0,NY,NX)=OSC(M,K,L0,NY,NX)-FXOSC
        FXOSA=FXO*OSA(M,K,L0,NY,NX)
        OSA(M,K,L1,NY,NX)=OSA(M,K,L1,NY,NX)+FXOSA
        OSA(M,K,L0,NY,NX)=OSA(M,K,L0,NY,NX)-FXOSA
        FXOSN=FXO*OSN(M,K,L0,NY,NX)
        OSN(M,K,L1,NY,NX)=OSN(M,K,L1,NY,NX)+FXOSN
        OSN(M,K,L0,NY,NX)=OSN(M,K,L0,NY,NX)-FXOSN
        FXOSP=FXO*OSP(M,K,L0,NY,NX)
        OSP(M,K,L1,NY,NX)=OSP(M,K,L1,NY,NX)+FXOSP
        OSP(M,K,L0,NY,NX)=OSP(M,K,L0,NY,NX)-FXOSP
      ENDDO
    ENDDO
!
!     ROOTS (why return?) 05/21/2022, jyt
!
    return
    DO NZ=1,NP(NY,NX)
      IF(WTRTL(ipltroot,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.WTRTL(ipltroot,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        IF(L0.EQ.L.OR.DPTHZ(L1,NY,NX).LE.ZERO)THEN
          FRO=FO
        ELSE
          FRO=AMIN1(0.5,FO*DPTHZ(L0,NY,NX)/DPTHZ(L1,NY,NX))
        ENDIF

        DO  N=1,MY(NZ,NY,NX)
          DO NTG=idg_beg,idg_end-1
            FXGA=FRO*trcg_rootml(NTG,N,L0,NZ,NY,NX)
            trcg_rootml(NTG,N,L1,NZ,NY,NX)=trcg_rootml(NTG,N,L1,NZ,NY,NX)+FXGA
            trcg_rootml(NTG,N,L0,NZ,NY,NX)=trcg_rootml(NTG,N,L0,NZ,NY,NX)-FXGA

            FXGP=FRO*trcs_rootml(NTG,N,L0,NZ,NY,NX)
            trcs_rootml(NTG,N,L1,NZ,NY,NX)=trcs_rootml(NTG,N,L1,NZ,NY,NX)+FXGP
            trcs_rootml(NTG,N,L0,NZ,NY,NX)=trcs_rootml(NTG,N,L0,NZ,NY,NX)-FXGP
          ENDDO

          DO  NR=1,NRT(NZ,NY,NX)
            DO NE=1,npelms
              FXWTRT1E=FRO*WTRT1E(NE,N,L0,NR,NZ,NY,NX)
              WTRT1E(NE,N,L1,NR,NZ,NY,NX)=WTRT1E(NE,N,L1,NR,NZ,NY,NX)+FXWTRT1E
              WTRT1E(NE,N,L0,NR,NZ,NY,NX)=WTRT1E(NE,N,L0,NR,NZ,NY,NX)-FXWTRT1E

              FXWTRT2E=FRO*WTRT2E(NE,N,L0,NR,NZ,NY,NX)
              WTRT2E(NE,N,L1,NR,NZ,NY,NX)=WTRT2E(NE,N,L1,NR,NZ,NY,NX)+FXWTRT2E
              WTRT2E(NE,N,L0,NR,NZ,NY,NX)=WTRT2E(NE,N,L0,NR,NZ,NY,NX)-FXWTRT2E
            ENDDO
            FXRTLG1=FRO*RTLG1(N,L0,NR,NZ,NY,NX)
            RTLG1(N,L1,NR,NZ,NY,NX)=RTLG1(N,L1,NR,NZ,NY,NX)+FXRTLG1
            RTLG1(N,L0,NR,NZ,NY,NX)=RTLG1(N,L0,NR,NZ,NY,NX)-FXRTLG1
            FXRTLG2=FRO*RTLG2(N,L0,NR,NZ,NY,NX)
            RTLG2(N,L1,NR,NZ,NY,NX)=RTLG2(N,L1,NR,NZ,NY,NX)+FXRTLG2
            RTLG2(N,L0,NR,NZ,NY,NX)=RTLG2(N,L0,NR,NZ,NY,NX)-FXRTLG2
            FXRTN2=FRO*RTN2(N,L0,NR,NZ,NY,NX)
            RTN2(N,L1,NR,NZ,NY,NX)=RTN2(N,L1,NR,NZ,NY,NX)+FXRTN2
            RTN2(N,L0,NR,NZ,NY,NX)=RTN2(N,L0,NR,NZ,NY,NX)-FXRTN2
          ENDDO
          DO NE=1,npelms
            FXEPOOLR=FRO*EPOOLR(NE,N,L0,NZ,NY,NX)
            EPOOLR(NE,N,L1,NZ,NY,NX)=EPOOLR(NE,N,L1,NZ,NY,NX)+FXEPOOLR
            EPOOLR(NE,N,L0,NZ,NY,NX)=EPOOLR(NE,N,L0,NZ,NY,NX)-FXEPOOLR
          ENDDO

          FXWTRTL=FRO*WTRTL(N,L0,NZ,NY,NX)
          WTRTL(N,L1,NZ,NY,NX)=WTRTL(N,L1,NZ,NY,NX)+FXWTRTL
          WTRTL(N,L0,NZ,NY,NX)=WTRTL(N,L0,NZ,NY,NX)-FXWTRTL
          FXWTRTD=FRO*WTRTD(N,L0,NZ,NY,NX)
          WTRTD(N,L1,NZ,NY,NX)=WTRTD(N,L1,NZ,NY,NX)+FXWTRTD
          WTRTD(N,L0,NZ,NY,NX)=WTRTD(N,L0,NZ,NY,NX)-FXWTRTD
          FXWSRTL=FRO*WSRTL(N,L0,NZ,NY,NX)
          WSRTL(N,L1,NZ,NY,NX)=WSRTL(N,L1,NZ,NY,NX)+FXWSRTL
          WSRTL(N,L0,NZ,NY,NX)=WSRTL(N,L0,NZ,NY,NX)-FXWSRTL
          FXRTN1=FRO*RTN1(N,L0,NZ,NY,NX)
          RTN1(N,L1,NZ,NY,NX)=RTN1(N,L1,NZ,NY,NX)+FXRTN1
          RTN1(N,L0,NZ,NY,NX)=RTN1(N,L0,NZ,NY,NX)-FXRTN1
          FXRTNL=FRO*RTNL(N,L0,NZ,NY,NX)
          RTNL(N,L1,NZ,NY,NX)=RTNL(N,L1,NZ,NY,NX)+FXRTNL
          RTNL(N,L0,NZ,NY,NX)=RTNL(N,L0,NZ,NY,NX)-FXRTNL
          FXRTLGP=FRO*RTLGP(N,L0,NZ,NY,NX)
          RTLGP(N,L1,NZ,NY,NX)=RTLGP(N,L1,NZ,NY,NX)+FXRTLGP
          RTLGP(N,L0,NZ,NY,NX)=RTLGP(N,L0,NZ,NY,NX)-FXRTLGP
          FXRTDNP=FRO*RTDNP(N,L0,NZ,NY,NX)
          RTDNP(N,L1,NZ,NY,NX)=RTDNP(N,L1,NZ,NY,NX)+FXRTDNP
          RTDNP(N,L0,NZ,NY,NX)=RTDNP(N,L0,NZ,NY,NX)-FXRTDNP
          FXRTVLP=FRO*RTVLP(N,L0,NZ,NY,NX)
          RTVLP(N,L1,NZ,NY,NX)=RTVLP(N,L1,NZ,NY,NX)+FXRTVLP
          RTVLP(N,L0,NZ,NY,NX)=RTVLP(N,L0,NZ,NY,NX)-FXRTVLP
          FXRTVLW=FRO*RTVLW(N,L0,NZ,NY,NX)
          RTVLW(N,L1,NZ,NY,NX)=RTVLW(N,L1,NZ,NY,NX)+FXRTVLW
          RTVLW(N,L0,NZ,NY,NX)=RTVLW(N,L0,NZ,NY,NX)-FXRTVLW
          FXRRAD1=FRO*RRAD1(N,L0,NZ,NY,NX)
          RRAD1(N,L1,NZ,NY,NX)=RRAD1(N,L1,NZ,NY,NX)+FXRRAD1
          RRAD1(N,L0,NZ,NY,NX)=RRAD1(N,L0,NZ,NY,NX)-FXRRAD1
          FXRRAD2=FRO*RRAD2(N,L0,NZ,NY,NX)
          RRAD2(N,L1,NZ,NY,NX)=RRAD2(N,L1,NZ,NY,NX)+FXRRAD2
          RRAD2(N,L0,NZ,NY,NX)=RRAD2(N,L0,NZ,NY,NX)-FXRRAD2
          FXRTARP=FRO*RTARP(N,L0,NZ,NY,NX)
          RTARP(N,L1,NZ,NY,NX)=RTARP(N,L1,NZ,NY,NX)+FXRTARP
          RTARP(N,L0,NZ,NY,NX)=RTARP(N,L0,NZ,NY,NX)-FXRTARP
          FXRTLGA=FRO*RTLGA(N,L0,NZ,NY,NX)
          RTLGA(N,L1,NZ,NY,NX)=RTLGA(N,L1,NZ,NY,NX)+FXRTLGA
          RTLGA(N,L0,NZ,NY,NX)=RTLGA(N,L0,NZ,NY,NX)-FXRTLGA
        ENDDO
!
!     ROOT NODULES
!
        DO NE=1,npelms
          FXWTNDLE=FRO*WTNDLE(NE,L0,NZ,NY,NX)
          WTNDLE(NE,L1,NZ,NY,NX)=WTNDLE(NE,L1,NZ,NY,NX)+FXWTNDLE
          WTNDLE(NE,L0,NZ,NY,NX)=WTNDLE(NE,L0,NZ,NY,NX)-FXWTNDLE

          FXEPOOLN=FRO*EPOOLN(NE,L0,NZ,NY,NX)
          EPOOLN(NE,L1,NZ,NY,NX)=EPOOLN(NE,L1,NZ,NY,NX)+FXEPOOLN
          EPOOLN(NE,L0,NZ,NY,NX)=EPOOLN(NE,L0,NZ,NY,NX)-FXEPOOLN

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
  IF(SoilFracAsMacP(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP(L0,NY,NX).GT.ZERO)THEN

    DO NTS=ids_beg,ids_end
      FXSH=FHO*trc_soHml(NTS,L0,NY,NX)
      trc_soHml(NTS,L1,NY,NX)=trc_soHml(NTS,L1,NY,NX)+FXSH
      trc_soHml(NTS,L0,NY,NX)=trc_soHml(NTS,L0,NY,NX)-FXSH
    ENDDO

!
!     SOIL MACROPORE SOLUBLE SALTS
!
    IF(salt_model)THEN
      DO NTSA=idsa_beg,idsab_end
        FXH=FHO*trcsa_soHml(NTSA,L0,NY,NX)
        trcsa_soHml(NTSA,L1,NY,NX)=trcsa_soHml(NTSA,L1,NY,NX)+FXH
        trcsa_soHml(NTSA,L0,NY,NX)=trcsa_soHml(NTSA,L0,NY,NX)-FXH
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
  FXH1POB=FWO*trc_solml(ids_H1PO4B,L0,NY,NX)
  trc_solml(ids_H1PO4B,L1,NY,NX)=trc_solml(ids_H1PO4B,L1,NY,NX)+FXH1POB
  trc_solml(ids_H1PO4B,L0,NY,NX)=trc_solml(ids_H1PO4B,L0,NY,NX)-FXH1POB

  FXH2POB=FWO*trc_solml(ids_H2PO4B,L0,NY,NX)
  trc_solml(ids_H2PO4B,L1,NY,NX)=trc_solml(ids_H2PO4B,L1,NY,NX)+FXH2POB
  trc_solml(ids_H2PO4B,L0,NY,NX)=trc_solml(ids_H2PO4B,L0,NY,NX)-FXH2POB

  IF(salt_model)THEN
    DO NTSAB=idsab_beg,idsab_end
      FXB=FWO*trcsa_solml(NTSAB,L0,NY,NX)
      trcsa_solml(NTSAB,L1,NY,NX)=trcsa_solml(NTSAB,L1,NY,NX)+FXB
      trcsa_solml(NTSAB,L0,NY,NX)=trcsa_solml(NTSAB,L0,NY,NX)-FXB
    ENDDO
  ENDIF
!
!     SOIL ADSORBED CATIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L.OR.CEC(L0,NY,NX).LE.ZERO)THEN
    FCO=FO
  ELSE
    FCO=AMIN1(0.5_r8,FO*CEC(L1,NY,NX)/CEC(L0,NY,NX))
  ENDIF
  DO NTX=idx_beg,idx_cation_end
    FXXC=FCO*trcx_solml(NTX,L0,NY,NX)
    trcx_solml(NTX,L1,NY,NX)=trcx_solml(NTX,L1,NY,NX)+FXXC
    trcx_solml(NTX,L0,NY,NX)=trcx_solml(NTX,L0,NY,NX)-FXXC
  ENDDO
!
!     SOIL ADSORBED ANIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L.OR.AEC(L0,NY,NX).LE.ZERO)THEN
    FAO=FO
  ELSE
    FAO=AMIN1(0.5_r8,FO*AEC(L1,NY,NX)/AEC(L0,NY,NX))
  ENDIF

  DO NTX=idx_cation_end+1,idx_end
    FXXA=FAO*trcx_solml(NTX,L0,NY,NX)
    trcx_solml(NTX,L1,NY,NX)=trcx_solml(NTX,L1,NY,NX)+FXXA
    trcx_solml(NTX,L0,NY,NX)=trcx_solml(NTX,L0,NY,NX)-FXXA
  ENDDo
!
!     SOIL PRECIPITATES IN BAND, NON-BAND
!
  DO NTP=idsp_beg,idsp_end
    FXP=AMIN1(FX*trcp_salml(NTP,L,NY,NX),trcp_salml(NTP,L0,NY,NX))
    trcp_salml(NTP,L1,NY,NX)=trcp_salml(NTP,L1,NY,NX)+FXP
    trcp_salml(NTP,L0,NY,NX)=trcp_salml(NTP,L0,NY,NX)-FXP
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
    FXG=FWO*trc_gasml(NTG,L0,NY,NX)
    trc_gasml(NTG,L1,NY,NX)=trc_gasml(NTG,L1,NY,NX)+FXG
    trc_gasml(NTG,L0,NY,NX)=trc_gasml(NTG,L0,NY,NX)-FXG

    FXG=FWO*trc_solml(NTG,L0,NY,NX)
    trc_solml(NTG,L1,NY,NX)=trc_solml(NTG,L1,NY,NX)+FXG
    trc_solml(NTG,L0,NY,NX)=trc_solml(NTG,L0,NY,NX)-FXG

  ENDDO
! add NH3
  NTG=idg_NH3
  FXG=FWO*trc_gasml(NTG,L0,NY,NX)
  trc_gasml(NTG,L1,NY,NX)=trc_gasml(NTG,L1,NY,NX)+FXG
  trc_gasml(NTG,L0,NY,NX)=trc_gasml(NTG,L0,NY,NX)-FXG

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
    FXZN=AMIN1(FX*FertN_soil(NTF,L,NY,NX),FertN_soil(NTF,L0,NY,NX))
    FertN_soil(NTF,L1,NY,NX)=FertN_soil(NTF,L1,NY,NX)+FXZN
    FertN_soil(NTF,L0,NY,NX)=FertN_soil(NTF,L0,NY,NX)-FXZN
  ENDDO

  IF (L0>0) then
    DO NTF=ifertnb_beg,ifertnb_end
      FXZN=AMIN1(FX*FertN_band(NTF,L,NY,NX),FertN_band(NTF,L0,NY,NX))
      FertN_band(NTF,L1,NY,NX)=FertN_band(NTF,L1,NY,NX)+FXZN
      FertN_band(NTF,L0,NY,NX)=FertN_band(NTF,L0,NY,NX)-FXZN
    ENDDO
  endif
!
!     SOIL N,P SOLUTES IN BAND, NON-BAND
!
  DO NTS=ids_nuts_beg,ids_nuts_end
    if(NTS/=ids_H2PO4B .and. NTS/=ids_H1PO4B)THEN
      FXNUT=FWO*trc_solml(NTS,L0,NY,NX)
      trc_solml(NTS,L1,NY,NX)=trc_solml(NTS,L1,NY,NX)+FXNUT
      trc_solml(NTS,L0,NY,NX)=trc_solml(NTS,L0,NY,NX)-FXNUT
    ENDIF
  ENDDO

!
!     SOIL SALT SOLUTES
!
  IF(salt_model)THEN
    DO NTSA=idsa_beg,idsa_end
      FXZ=FWO*trcsa_solml(NTSA,L0,NY,NX)
      trcsa_solml(NTSA,L1,NY,NX)=trcsa_solml(NTSA,L1,NY,NX)+FXZ
      trcsa_solml(NTSA,L0,NY,NX)=trcsa_solml(NTSA,L0,NY,NX)-FXZ
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
      IF(L.EQ.NU(NY,NX).OR.CumDepth2LayerBottom(L-1,NY,NX).LT.DPNH4(NY,NX))THEN
        WDNHBDL=WDNHB(L,NY,NX)*DLYR(3,L,NY,NX)
        WDNHBD0=WDNHB(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDNHBD1=WDNHB(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDNHB=AMIN1(FX*WDNHBDL,WDNHBD0)
        WDNHBD1=WDNHBD1+FXWDNHB
        WDNHBD0=WDNHBD0-FXWDNHB
        WDNHB(L1,NY,NX)=WDNHBD1/DLYR(3,L1,NY,NX)
        WDNHB(L0,NY,NX)=WDNHBD0/DLYR(3,L0,NY,NX)
        IF(CumDepth2LayerBottom(L,NY,NX).GE.DPNH4(NY,NX))THEN
          FXDPNHB=AMIN1(FX*DPNHB(L,NY,NX),DPNHB(L0,NY,NX))
          DPNHB(L1,NY,NX)=DPNHB(L1,NY,NX)+FXDPNHB
          DPNHB(L0,NY,NX)=DPNHB(L0,NY,NX)-FXDPNHB
        ENDIF
        trcs_VLN(ids_NH4B,L1,NY,NX)=AZMAX1(AMIN1(0.999,WDNHB(L1,NY,NX) &
          /ROWN(NY,NX)*DPNHB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        trcs_VLN(ids_NH4B,L0,NY,NX)=AZMAX1(AMIN1(0.999,WDNHB(L0,NY,NX) &
          /ROWN(NY,NX)*DPNHB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        trcs_VLN(ids_NH4,L1,NY,NX)=1.0_r8-trcs_VLN(ids_NH4B,L1,NY,NX)
        trcs_VLN(ids_NH4,L0,NY,NX)=1.0_r8-trcs_VLN(ids_NH4B,L0,NY,NX)

        trcs_VLN(idg_NH3B,L1,NY,NX)=trcs_VLN(ids_NH4B,L1,NY,NX)
        trcs_VLN(idg_NH3B,L0,NY,NX)=trcs_VLN(ids_NH4B,L0,NY,NX)
        trcs_VLN(idg_NH3,L1,NY,NX)=trcs_VLN(ids_NH4,L1,NY,NX)
        trcs_VLN(idg_NH3,L0,NY,NX)=trcs_VLN(ids_NH4,L0,NY,NX)
      ENDIF
    ENDIF
    IF(IFNOB(NY,NX).EQ.1.AND.ROWO(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CumDepth2LayerBottom(L-1,NY,NX).LT.DPNO3(NY,NX))THEN
        WDNOBDL=WDNOB(L,NY,NX)*DLYR(3,L,NY,NX)
        WDNOBD0=WDNOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDNOBD1=WDNOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDNOB=AMIN1(FX*WDNOBDL,WDNOBD0)
        WDNOBD1=WDNOBD1+FXWDNOB
        WDNOBD0=WDNOBD0-FXWDNOB
        WDNOB(L1,NY,NX)=WDNOBD1/DLYR(3,L1,NY,NX)
        WDNOB(L0,NY,NX)=WDNOBD0/DLYR(3,L0,NY,NX)
        IF(CumDepth2LayerBottom(L,NY,NX).GE.DPNO3(NY,NX))THEN
          FXDPNOB=AMIN1(FX*DPNOB(L,NY,NX),DPNOB(L0,NY,NX))
          DPNOB(L1,NY,NX)=DPNOB(L1,NY,NX)+FXDPNOB
          DPNOB(L0,NY,NX)=DPNOB(L0,NY,NX)-FXDPNOB
        ENDIF
        trcs_VLN(ids_NO3B,L1,NY,NX)=AZMAX1(AMIN1(0.999_r8,WDNOB(L1,NY,NX) &
          /ROWO(NY,NX)*DPNOB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        trcs_VLN(ids_NO3B,L0,NY,NX)=AZMAX1(AMIN1(0.999_r8,WDNOB(L0,NY,NX) &
          /ROWO(NY,NX)*DPNOB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        trcs_VLN(ids_NO3,L1,NY,NX)=1.0_r8-trcs_VLN(ids_NO3B,L1,NY,NX)
        trcs_VLN(ids_NO3,L0,NY,NX)=1.0_r8-trcs_VLN(ids_NO3B,L0,NY,NX)

        trcs_VLN(ids_NO2,L1,NY,NX)=trcs_VLN(ids_NO3,L1,NY,NX)
        trcs_VLN(ids_NO2,L0,NY,NX)=trcs_VLN(ids_NO3,L0,NY,NX)
        trcs_VLN(ids_NO2B,L1,NY,NX)=trcs_VLN(ids_NO3B,L1,NY,NX)
        trcs_VLN(ids_NO2B,L0,NY,NX)=trcs_VLN(ids_NO3B,L0,NY,NX)
      ENDIF
    ENDIF
    IF(IFPOB(NY,NX).EQ.1.AND.ROWP(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CumDepth2LayerBottom(L-1,NY,NX).LT.DPPO4(NY,NX))THEN
        WDPOBDL=WDPOB(L,NY,NX)*DLYR(3,L,NY,NX)
        WDPOBD0=WDPOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDPOBD1=WDPOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDPOB=AMIN1(FX*WDPOBDL,WDPOBD0)
        WDPOBD1=WDPOBD1+FXWDPOB
        WDPOBD0=WDPOBD0-FXWDPOB
        WDPOB(L1,NY,NX)=WDPOBD1/DLYR(3,L1,NY,NX)
        WDPOB(L0,NY,NX)=WDPOBD0/DLYR(3,L0,NY,NX)
        IF(CumDepth2LayerBottom(L,NY,NX).GE.DPPO4(NY,NX))THEN
          FXDPPOB=AMIN1(FX*DPPOB(L,NY,NX),DPPOB(L0,NY,NX))
          DPPOB(L1,NY,NX)=DPPOB(L1,NY,NX)+FXDPPOB
          DPPOB(L0,NY,NX)=DPPOB(L0,NY,NX)-FXDPPOB
        ENDIF
        trcs_VLN(ids_H1PO4B,L1,NY,NX)=AZMAX1(AMIN1(0.999,WDPOB(L1,NY,NX) &
          /ROWP(NY,NX)*DPPOB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        trcs_VLN(ids_H1PO4B,L0,NY,NX)=AZMAX1(AMIN1(0.999,WDPOB(L0,NY,NX) &
          /ROWP(NY,NX)*DPPOB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        trcs_VLN(ids_H1PO4,L1,NY,NX)=1.0_r8-trcs_VLN(ids_H1PO4B,L1,NY,NX)
        trcs_VLN(ids_H1PO4,L0,NY,NX)=1.0_r8-trcs_VLN(ids_H1PO4B,L0,NY,NX)

        trcs_VLN(ids_H2PO4B,L1,NY,NX)=trcs_VLN(ids_H1PO4B,L1,NY,NX)
        trcs_VLN(ids_H2PO4B,L0,NY,NX)=trcs_VLN(ids_H1PO4B,L0,NY,NX)
        trcs_VLN(ids_H2PO4,L1,NY,NX)=trcs_VLN(ids_H1PO4,L1,NY,NX)
        trcs_VLN(ids_H2PO4,L0,NY,NX)=trcs_VLN(ids_H1PO4,L0,NY,NX)
      ENDIF
    ENDIF
  ENDIF
!
!     SOIL MINERALS
!
  IF(L0.EQ.L.OR.SoiBulkDensityt0(L0,NY,NX).LE.ZERO)THEN
    FBO=FX
  ELSE
    FBO=AMIN1(0.1,FX*SoiBulkDensityt0(L1,NY,NX)/SoiBulkDensityt0(L0,NY,NX))
  ENDIF
!     SoiBulkDensity(L1,NY,NX)=(1.0-FO)*SoiBulkDensity(L1,NY,NX)+FO*SoiBulkDensityt0(L0,NY,NX)
  PH(L1,NY,NX)=(1.0_r8-FO)*PH(L1,NY,NX)+FO*PH(L0,NY,NX)
  FXSAND=FBO*SAND(L0,NY,NX)
  SAND(L1,NY,NX)=SAND(L1,NY,NX)+FXSAND
  SAND(L0,NY,NX)=SAND(L0,NY,NX)-FXSAND
  FXSILT=FBO*SILT(L0,NY,NX)
  SILT(L1,NY,NX)=SILT(L1,NY,NX)+FXSILT
  SILT(L0,NY,NX)=SILT(L0,NY,NX)-FXSILT
  FXCLAY=FBO*CLAY(L0,NY,NX)
  CLAY(L1,NY,NX)=CLAY(L1,NY,NX)+FXCLAY
  CLAY(L0,NY,NX)=CLAY(L0,NY,NX)-FXCLAY
  FXROCK=FBO*ROCK(L0,NY,NX)
  ROCK(L1,NY,NX)=ROCK(L1,NY,NX)+FXROCK
  ROCK(L0,NY,NX)=ROCK(L0,NY,NX)-FXROCK
!
!     SOIL WATER AND HEAT
!
  IF(SoilFracAsMacP(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP(L0,NY,NX).GT.ZERO)THEN
    IF(L0.EQ.L.OR.SoilFracAsMacPt0(L0,NY,NX).LE.ZERO)THEN
      FHO=FO
    ELSE
      FHO=AMIN1(0.5_r8,FO*SoilFracAsMacPt0(L1,NY,NX)/SoilFracAsMacPt0(L0,NY,NX))
    ENDIF
    SoilFracAsMacP(L1,NY,NX)=(1.0_r8-FO)*SoilFracAsMacP(L1,NY,NX)+FO*SoilFracAsMacP(L0,NY,NX)
    FXVOLWH=FHO*VLWatMacP(L0,NY,NX)
    VLWatMacP(L1,NY,NX)=VLWatMacP(L1,NY,NX)+FXVOLWH
    VLWatMacP(L0,NY,NX)=VLWatMacP(L0,NY,NX)-FXVOLWH
    FXVOLIH=FHO*VLiceMacP(L0,NY,NX)
    VLiceMacP(L1,NY,NX)=VLiceMacP(L1,NY,NX)+FXVOLIH
    VLiceMacP(L0,NY,NX)=VLiceMacP(L0,NY,NX)-FXVOLIH
    FXVOLAH=FHO*VLMacP(L0,NY,NX)
    VLMacP(L1,NY,NX)=VLMacP(L1,NY,NX)+FXVOLAH
    VLMacP(L0,NY,NX)=VLMacP(L0,NY,NX)-FXVOLAH
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

  IF(L0.EQ.L.OR.(L0>0 .and. POROSI(L0,NY,NX).LE.ZERO))THEN
    FWO=FO
  ELSE
    FWO=AMIN1(0.5,FO*POROSI(L1,NY,NX)/POROSI(L0,NY,NX))
  ENDIF
!             FXSCNV=FWO*SatHydroCondVert(L0,NY,NX)
!             SatHydroCondVert(L1,NY,NX)=SatHydroCondVert(L1,NY,NX)+FXSCNV
!             SatHydroCondVert(L0,NY,NX)=SatHydroCondVert(L0,NY,NX)-FXSCNV
!             FXSCNH=FWO*SatHydroCondHrzn(L0,NY,NX)
!             SatHydroCondHrzn(L1,NY,NX)=SatHydroCondHrzn(L1,NY,NX)+FXSCNH
!             SatHydroCondHrzn(L0,NY,NX)=SatHydroCondHrzn(L0,NY,NX)-FXSCNH
  IF(L0.EQ.0)THEN
    FXVOLW=FX*AZMAX1(XVOLWP-VOLWD(NY,NX))
  ELSE
    FXVOLW=FWO*VLWatMicP(L0,NY,NX)
  ENDIF
  VLWatMicP(L1,NY,NX)=VLWatMicP(L1,NY,NX)+FXVOLW
  VLWatMicP(L0,NY,NX)=VLWatMicP(L0,NY,NX)-FXVOLW
!     IF(VLiceMicP(L1,NY,NX).GT.ZEROS(NY,NX))THEN
  FXVOLI=FWO*VLiceMicP(L0,NY,NX)
  VLiceMicP(L1,NY,NX)=VLiceMicP(L1,NY,NX)+FXVOLI
  VLiceMicP(L0,NY,NX)=VLiceMicP(L0,NY,NX)-FXVOLI
!     ENDIF
!     FXVOLA=FWO*VLMicP(L0,NY,NX)
!     IF(L1.NE.NU(NY,NX))THEN
!     VLMicP(L1,NY,NX)=VLMicP(L1,NY,NX)+FXVOLA
!     ENDIF
!     IF(L0.NE.NU(NY,NX))THEN
!     VLMicP(L0,NY,NX)=VLMicP(L0,NY,NX)-FXVOLA
!     ENDIF
  FXVLSoilMicP=FWO*VLSoilMicP(L0,NY,NX)
  VLSoilMicP(L1,NY,NX)=VLSoilMicP(L1,NY,NX)+FXVLSoilMicP
  VLSoilMicP(L0,NY,NX)=VLSoilMicP(L0,NY,NX)-FXVLSoilMicP
  FXVOLWX=FWO*VLWatMicPX(L0,NY,NX)
  VLWatMicPX(L1,NY,NX)=VLWatMicPX(L1,NY,NX)+FXVOLWX
  VLWatMicPX(L0,NY,NX)=VLWatMicPX(L0,NY,NX)-FXVOLWX
  FXVHCM=FWO*VHeatCapacitySoilM(L0,NY,NX)
  VHeatCapacitySoilM(L1,NY,NX)=VHeatCapacitySoilM(L1,NY,NX)+FXVHCM
  VHeatCapacitySoilM(L0,NY,NX)=VHeatCapacitySoilM(L0,NY,NX)-FXVHCM
  FXENGY=TKS(L0,NY,NX)*(FXVHCM+cpw*FXVOLW+cpi*FXVOLI)
  ENGY1=VHeatCapacity(L1,NY,NX)*TKS(L1,NY,NX)+FXENGY
  ENGY0=VHeatCapacity(L0,NY,NX)*TKS(L0,NY,NX)-FXENGY
  VHeatCapacity(L1,NY,NX)=VHeatCapacity(L1,NY,NX)+FXVHCM+cpw*FXVOLW+cpi*FXVOLI
  VHeatCapacity(L0,NY,NX)=VHeatCapacity(L0,NY,NX)-FXVHCM-cpw*FXVOLW-cpi*FXVOLI
  IF(VHeatCapacity(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L1,NY,NX)=ENGY1/VHeatCapacity(L1,NY,NX)
  ELSE
    TKS(L1,NY,NX)=TKS(L,NY,NX)
  ENDIF
  TCS(L1,NY,NX)=units%Kelvin2Celcius(TKS(L1,NY,NX))
  IF(VHeatCapacity(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L0,NY,NX)=ENGY0/VHeatCapacity(L0,NY,NX)
  ELSE
    TKS(L0,NY,NX)=TKS(L,NY,NX)
  ENDIF
  TCS(L0,NY,NX)=units%Kelvin2Celcius(TKS(L0,NY,NX))
  end subroutine MoveHeatWat

end module SoilLayerDynMod
