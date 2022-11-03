module SoilLayerDynMod
! Description:
! subroutines to do soil relayering
  use data_kind_mod, only : r8 => SHR_KIND_R8
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
  use EcoSIMConfig, only : ndbiomcp => ndbiomcpc
  USE TFlxTypeMod , ONLY : TSEDER,TDLYXF,TDYLXC,TDVOLI,TDORGC
implicit none


  private

  character(len=*), parameter :: mod_filename = __FILE__

  real(r8) :: DDLYRX(3)
  real(r8) :: XVOLWP,WDPOBDL,WDNOBD1,WDPOBD0,WDPOBD1
  real(r8) :: WDNHBDL,WDNHBD0,WDNHBD1,WDNOBDL,WDNOBD0
  real(r8), PARAMETER :: ZEROC=0.1E-03_r8

  public :: RelayerSoilProfile
  contains


  subroutine RelayerSoilProfile(NY,NX,DORGC,DVOLI,UDVOLI,UDLYXF)
  !
  !Description:
  !relayer the soil profiles
  implicit none
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: DORGC(JZ,JY,JX)  !change in organic matter, initial-final
  REAL(R8),INTENT(IN) :: DVOLI(JZ,JY,JX)  !change in ice volume, initial-final
  real(r8),intent(inout):: UDVOLI,UDLYXF

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
  IF(IERSNG.GE.0)THEN
!  soil erosion flag is active
    IF(BKDS(NU(NY,NX),NY,NX).LE.ZERO)THEN
      ICHKLX=0
    ELSE
      !it is a soil column
      ICHKLX=1
    ENDIF

    call SoilSubsidence(ICHKLX,NY,NX,DORGC,DVOLI,UDLYXF,UDVOLI,CDPTHX,CDPTHY,IFLGL)

    !
    !     RECALCULATE SOIL LAYER THICKNESS
    !
    ICHKL=0
    DO 245 L=NU(NY,NX),NL(NY,NX)-1
      DO 230 NN=1,3

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
            IF(BKDS(L0,NY,NX).LE.ZERO)THEN
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

                call MoveMMPoreSolute(L0,L1,NY,NX,FHO)
              ENDIF
!     SOIL ORGANIC MATTER
              call MoveSOM(L0,L1,L,NY,NX,FO,IFLGL)

              IF(NN.EQ.1)THEN
                IF(BKDS(L0,NY,NX).LE.ZERO.AND.BKDS(L1,NY,NX).LE.ZERO &
                  .AND.VOLW(L0,NY,NX)+VOLI(L0,NY,NX).LE.ZEROS(NY,NX))THEN
                  CDPTH(L1,NY,NX)=CDPTH(L0,NY,NX)
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
              IF(L.EQ.NL(NY,NX)-1.AND.CDPTH(NL(NY,NX),NY,NX)-CDPTH(L,NY,NX).LE.ZEROC)THEN
                CDPTH(L,NY,NX)=CDPTH(L,NY,NX)+DLYR(3,NL(NY,NX),NY,NX)
                CDPTHY(L,NY,NX)=CDPTHY(L,NY,NX)+DLYR(3,NL(NY,NX),NY,NX)
                CDPTH(NL(NY,NX),NY,NX)=CDPTH(L,NY,NX)
                CDPTHY(NL(NY,NX),NY,NX)=CDPTHY(L,NY,NX)
                DLYR(3,NL(NY,NX),NY,NX)=0.0_r8
                NL(NY,NX)=L
              ENDIF
            ENDIF
          ENDIF
        ENDIF
230   CONTINUE
245   CONTINUE
  ENDIF

  end subroutine RelayerSoilProfile
!------------------------------------------------------------------------------------------
  subroutine SoilSubsidence(ICHKLX,NY,NX,DORGC,DVOLI,UDLYXF,UDVOLI,&
    CDPTHX,CDPTHY,IFLGL)
!
! IFLGL: c1, ponding water, c2, pond disappear, c3, pond reappare, c4: freeze-thaw, c5: erosion, c6: som change
  implicit none
  integer, intent(in) :: ICHKLX,NY,NX
  REAL(R8),INTENT(IN) :: DVOLI(JZ,JY,JX)
  real(r8),intent(in) :: DORGC(JZ,JY,JX)
  real(r8), intent(inout) :: UDLYXF,UDVOLI
  real(r8), intent(out) :: CDPTHX(JZ,JY,JX)
  real(r8), intent(inout) :: CDPTHY(0:JZ,JY,JX)
  integer, intent(inout) :: IFLGL(0:JZ,6)
  real(r8) :: DDLYX(0:JZ,6)
  real(r8) :: DDLYR(0:JZ,6)
  integer :: LX,LY,LL,NN,L
  real(r8) :: DDLYXC,DDLYXE
  real(r8) :: DDLYXF
  real(r8) :: DENSJ,DPTWI,DDLYXP

! begin_execution
! starting from bottom up
  DO  LX=NL(NY,NX),NU(NY,NX),-1
    !make a copy of the depth, bottom of the layer
    CDPTHX(LX,NY,NX)=CDPTH(LX,NY,NX)
    CDPTHY(LX,NY,NX)=CDPTH(LX,NY,NX)
    !
    !     POND, from water to soil
    !
! compute the change
    IF(BKDS(LX,NY,NX).LE.ZERO)THEN
      ! current layer is water
      IF(BKDS(LX+1,NY,NX).GT.ZERO.OR.ICHKLX.EQ.1)THEN
        !next layer is soil,
        DDLYXP=DLYR(3,LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
        DDLYX(LX,1)=DDLYXP+DDLYX(LX+1,1)
        DDLYR(LX,1)=DDLYX(LX+1,1)
        IFLGL(LX,1)=2
      ELSE
        !next layer is not soil
        !DDLYXP: soil equivalent depth
        !DLYRI: initial soil layer thickness, [m]
        DDLYXP=DLYRI(3,LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
        !DPTWI: water+ice total thickness
        DPTWI=(VOLW(LX+1,NY,NX)+VOLI(LX+1,NY,NX))/AREA(3,LX,NY,NX)
        !there is expansion in layer LX, or next layer has
        IF(DDLYXP.LT.-ZERO.OR.DPTWI.GT.ZERO)THEN
          DDLYX(LX,1)=DDLYXP+DDLYX(LX+1,1)
          DDLYR(LX,1)=AMIN1(DDLYX(LX+1,1),DPTWI)
          IF(DPTWI.GT.ZERO)THEN
! already has solid thickness
            IFLGL(LX,1)=1
          ELSE
! has no solid thickness
            IFLGL(LX,1)=2
          ENDIF
        ELSE
          DDLYXP=DLYR(3,LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
          DDLYX(LX,1)=DDLYXP+DDLYX(LX+1,1)
          DDLYR(LX,1)=DDLYX(LX+1,1)
          IFLGL(LX,1)=2
        ENDIF
      ENDIF
      !surface layer or layer above is still soil
      IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).GT.ZERO)THEN
        DDLYX(LX-1,1)=DDLYX(LX,1)
        DDLYR(LX-1,1)=DDLYX(LX,1)
        IFLGL(LX-1,1)=1
      ENDIF
      DDLYX(LX,4)=0.0_r8
      DDLYR(LX,4)=0.0_r8
      IFLGL(LX,4)=0
      DDLYX(LX,5)=0.0_r8
      DDLYR(LX,5)=0.0_r8
      IFLGL(LX,5)=0
      DDLYX(LX,6)=0.0_r8
      DDLYR(LX,6)=0.0_r8
      IFLGL(LX,6)=0
      !
      !     SOIL
      !
    ELSE
      !
      !     FREEZE-THAW
      !DVOLI: change in ice volume, initial-final
      !DENSI: mass density of ice, g/cm3
      !DENSJ: volume added per unit addition of ice (with respect to ice)
      !the change is put to layer thickness
      !DDLYXF: added thickness change
! 4: due to freeze-thaw
      IF(ABS(DVOLI(LX,NY,NX)).GT.ZEROS(NY,NX))THEN
        DENSJ=1._r8-DENSI
        DDLYXF=DVOLI(LX,NY,NX)*DENSJ/AREA(3,LX,NY,NX)
        !bottom layer
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,4)=DDLYXF
          DDLYR(LX,4)=0.0_r8
          IFLGL(LX,4)=0
        !layer not at the bottom
        ELSE
          DDLYX(LX,4)=DDLYXF+DDLYX(LX+1,4)
          DDLYR(LX,4)=DDLYX(LX+1,4)
          !    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
          IFLGL(LX,4)=0
          !top soil layer or water layer
          IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).LE.ZERO)THEN
            DDLYX(LX-1,4)=DDLYX(LX,4)
            DDLYR(LX-1,4)=DDLYX(LX,4)
            !    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
            IFLGL(LX-1,4)=0
          ENDIF
        ENDIF
! no change in ice volume
      ELSE

        DDLYXF=0.0_r8
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,4)=0.0_r8
          DDLYR(LX,4)=0.0_r8
          IFLGL(LX,4)=0
        ELSE
          DDLYX(LX,4)=DDLYX(LX+1,4)
          DDLYR(LX,4)=DDLYX(LX+1,4)
          IFLGL(LX,4)=0
          IF(LX.EQ.NU(NY,NX))THEN
            DDLYX(LX-1,4)=DDLYX(LX,4)
            DDLYR(LX-1,4)=DDLYX(LX,4)
            IFLGL(LX-1,4)=0
          ENDIF
        ENDIF
      ENDIF
      !total change in ice volume
      TDVOLI=TDVOLI+DVOLI(LX,NY,NX)
      TDLYXF=TDLYXF+DDLYXF
      UDVOLI=UDVOLI+DVOLI(LX,NY,NX)
      UDLYXF=UDLYXF+DDLYXF
      !
      !     EROSION
      !
      IF((IERSNG.EQ.1.OR.IERSNG.EQ.3).AND.ABS(TSEDER(NY,NX)).GT.ZEROS(NY,NX))THEN
        IF(LX.EQ.NL(NY,NX))THEN
!  5: due to sediment erosion
!         total soil layer reduction due to erosion
          DDLYXE=-TSEDER(NY,NX)/(BKVLNU(NY,NX)/VOLX(NU(NY,NX),NY,NX))
          DDLYX(LX,5)=DDLYXE
          DDLYR(LX,5)=DDLYXE
          IFLGL(LX,5)=1
        ELSE
          DDLYX(LX,5)=0.0_r8
          DDLYR(LX,5)=0.0_r8
          IFLGL(LX,5)=0
        ENDIF
      ELSE
          DDLYX(LX,5)=0.0_r8
          DDLYR(LX,5)=0.0_r8
          IFLGL(LX,5)=0
      ENDIF

      !
      !     SOC GAIN OR LOSS
      ! FHOL: macropore fraction
      ! DDLYXC: soil thickness added due to change in organic matter,
      ! keeping macropore fraction
      ! BKDSI: initial bulk density,
      IF((IERSNG.EQ.2.OR.IERSNG.EQ.3).AND.ABS(DORGC(LX,NY,NX)).GT.ZEROS(NY,NX))THEN
        DDLYXC=MWC2Soil*DORGC(LX,NY,NX) &
          /((1.0_r8-FHOL(LX,NY,NX))*BKDSI(LX,NY,NX))/AREA(3,LX,NY,NX)
! obtain diagnostics only for NX==1
        IF(NX.EQ.1)THEN
          TDORGC=TDORGC+DORGC(LX,NY,NX)
          TDYLXC=TDYLXC+DDLYXC
        ENDIF
! bottom layer
        IF(LX.EQ.NL(NY,NX).OR.BKDS(LX+1,NY,NX).LE.ZERO)THEN
          DDLYX(LX,6)=DDLYXC
          DDLYR(LX,6)=0.0_r8
          IFLGL(LX,6)=1
        ELSE
          DDLYX(LX,6)=DDLYXC+DDLYX(LX+1,6)
          DDLYR(LX,6)=DDLYX(LX+1,6)+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
          IFLGL(LX,6)=1
! top layer
          IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).LE.ZERO)THEN
            DDLYX(LX-1,6)=DDLYX(LX,6)
            DDLYR(LX-1,6)=DDLYX(LX,6)
            !    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
            IFLGL(LX-1,6)=1
          ENDIF
        ENDIF
      ELSE
! bottom layer
        IF(LX.EQ.NL(NY,NX))THEN
          DDLYX(LX,6)=0.0_r8
          DDLYR(LX,6)=0.0_r8
          IFLGL(LX,6)=0
        ELSE
          DDLYX(LX,6)=DDLYX(LX+1,6)
          DDLYR(LX,6)=DDLYX(LX+1,6)
          IFLGL(LX,6)=0
        ENDIF
      ENDIF
      DDLYX(LX,1)=0.0_r8
      DDLYR(LX,1)=0.0_r8
      IFLGL(LX,1)=0

    ENDIF
! apply the change
      !
      !     RESET SOIL LAYER DEPTHS
      !
    DO  NN=1,6
      !     IF(ABS(DDLYX(LX,NN)).GT.ZERO)IFLGS(NY,NX)=1
!     c2, and c3 are not implemented
      IF(NN.NE.2.AND.NN.NE.3)THEN
        !
        !     POND
        !
        IF(BKDS(LX,NY,NX).LE.ZERO)THEN
! there are some changes
          IF(IFLGL(LX,NN).NE.0)THEN
            CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
!           not top layer
            IF(LX.NE.NU(NY,NX).AND.IFLGL(LX,1).EQ.2)THEN
              DO  LL=LX-1,0,-1
                CDPTH(LL,NY,NX)=CDPTH(LL,NY,NX)+DDLYX(LX,NN)
                CDPTHY(LL,NY,NX)=CDPTHY(LL,NY,NX)+DDLYX(LX,NN)
              ENDDO
              DDLYX(LX,NN)=0.0_r8
            ENDIF
!           top layer
            IF(LX.EQ.NU(NY,NX))THEN
              CDPTH(LX-1,NY,NX)=CDPTH(LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
              CDPTHY(LX-1,NY,NX)=CDPTHY(LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
            ENDIF
          ENDIF
            !
            !     SOIL
            !
        ELSE
            !     IF(DDLYR(L,NN).NE.0.OR.DDLYR(LX,NN).NE.0)THEN
            !
            !     FREEZE-THAW
            !
          IF(NN.EQ.4)THEN
            CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
            !     CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
! top layer
            IF(LX.EQ.NU(NY,NX))THEN
              CDPTH(LX-1,NY,NX)=CDPTH(LX-1,NY,NX)+DDLYR(LX-1,NN)
              !     CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX-1,NN)

            ENDIF
          ENDIF
            !
            !     SET SURFACE ELEVATION FOR SOIL EROSION
            !
          IF(NN.EQ.5.AND.IFLGL(LX,NN).EQ.1)THEN
            CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)

            IF(LX.EQ.NU(NY,NX))THEN
              CDPTH(LX-1,NY,NX)=CDPTH(LX-1,NY,NX)+DDLYR(LX,NN)
              CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX,NN)
            ENDIF
          ENDIF
          !
          !     SET SOIL LAYER DEPTHS FOR CHANGES IN SOC
          !
          IF(NN.EQ.6.AND.IFLGL(LX,NN).EQ.1)THEN
            CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
            CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)

            IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).LE.ZERO)THEN
              CDPTH(LX-1,NY,NX)=CDPTH(LX-1,NY,NX)+DDLYR(LX-1,NN)
              CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX-1,NN)

              IF(BKDS(LX-1,NY,NX).LE.ZERO)THEN
                DO  LY=LX-2,0,-1
                  IF(BKDS(LY+1,NY,NX).LE.ZERO)THEN
                    CDPTH(LY,NY,NX)=CDPTH(LY,NY,NX)+DDLYR(LX-1,NN)
                    CDPTHY(LY,NY,NX)=CDPTHY(LY,NY,NX)+DDLYR(LX-1,NN)

                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    VOLY(LX,NY,NX)=VOLX(LX,NY,NX)
  ENDDO
  VOLY(0,NY,NX)=VOLW(0,NY,NX)+VOLI(0,NY,NX)
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

    IF((BKDS(L,NY,NX).LE.ZERO.AND.IFLGL(L,1).EQ.2) &
        .OR.(DLYR(3,L0,NY,NX).LE.ZEROC.AND.IFLGL(L,6).EQ.1))THEN
      FX=1.0_r8
      FO=1.0_r8
    ELSE
      IF(BKDS(L0,NY,NX).LE.ZERO)THEN
        DPTWI=(VOLW(L0,NY,NX)+VOLI(L0,NY,NX))/AREA(3,L0,NY,NX)
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
    IF(BKDS(L0,NY,NX).LE.ZERO)THEN
      DPTWI=(VOLW(L0,NY,NX)+VOLI(L0,NY,NX))/AREA(3,L0,NY,NX)
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
    DLYR(3,L,NY,NX)=CDPTH(L,NY,NX)-CDPTH(L-1,NY,NX)
    DLYRXX=DLYR(3,L,NY,NX)
    IF(IFLGL(L,1).EQ.0.AND.IFLGL(L+1,1).NE.0)THEN
      DDLYRX(NN)=0.0_r8
      IF(BKDS(L,NY,NX).LE.ZERO)THEN
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
    CDPTH(L,NY,NX)=CDPTH(L,NY,NX)+DDLYRY(L)
  !     CDPTHY(L,NY,NX)=CDPTHY(L,NY,NX)+DDLYRY(L)
    DLYR(3,L,NY,NX)=DLYR(3,L,NY,NX)+DDLYRY(L)
    DPTH(L,NY,NX)=0.5*(CDPTH(L,NY,NX)+CDPTH(L-1,NY,NX))
    CDPTHZ(L,NY,NX)=CDPTH(L,NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX)
    IF(L.EQ.NL(NY,NX)-1)THEN
      DLYR(3,L+1,NY,NX)=CDPTH(L+1,NY,NX)-CDPTH(L,NY,NX)
      DPTH(L+1,NY,NX)=0.5*(CDPTH(L+1,NY,NX)+CDPTH(L,NY,NX))
      CDPTHZ(L+1,NY,NX)=CDPTH(L+1,NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX)
    ENDIF
    IF(L.EQ.NU(NY,NX))THEN
      DPTHZ(L,NY,NX)=0.5*CDPTHZ(L,NY,NX)
    !     DDLYRX(NN)=DDLYRX(NN)+DDLYR(L,5)
    ELSE
      DPTHZ(L,NY,NX)=0.5*(CDPTHZ(L,NY,NX)+CDPTHZ(L-1,NY,NX))
    ENDIF
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
    !     DDLYRX(NN)=CDPTH(L,NY,NX)-CDPTHX(L,NY,NX)
      DDLYRX(NN)=CDPTHY(L,NY,NX)-CDPTHX(L,NY,NX)
    ENDIF
!
  !     RESET POND SURFACE LAYER NUMBER IF LOST TO EVAPORATION
      !
  ELSEIF(NN.EQ.2)THEN
    IF((L.EQ.NU(NY,NX).AND.BKDS(NU(NY,NX),NY,NX).LE.ZERO) &
      .AND.(VHCP(NU(NY,NX),NY,NX).LE.VHCPNX(NY,NX) &
      .OR.NUM(NY,NX).GT.NU(NY,NX)))THEN
      NUX=NU(NY,NX)
      DO LL=NUX+1,NL(NY,NX)
        IF(VOLX(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
          NU(NY,NX)=LL
          DDLYRX(NN)=DLYR(3,NUX,NY,NX)
          IFLGL(L,NN)=1
          DLYR(3,NUX,NY,NX)=0.0_r8
          IF(BKDS(NUX,NY,NX).LE.ZERO)THEN
            VOLT(NUX,NY,NX)=AREA(3,NUX,NY,NX)*DLYR(3,NUX,NY,NX)
            VOLX(NUX,NY,NX)=VOLT(NUX,NY,NX)*FMPR(NUX,NY,NX)
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
    XVOLWP=AZMAX1(VOLW(0,NY,NX)-VOLWD(NY,NX))
    IF(L.EQ.NU(NY,NX).AND.CDPTH(0,NY,NX).GT.CDPTHI(NY,NX) &
      .AND.XVOLWP.GT.VOLWD(NY,NX)+VHCPNX(NY,NX)/cpw)THEN
          !     IF((BKDS(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))
          !    2.OR.(BKDS(L,NY,NX).LE.ZERO))THEN
      IF(BKDS(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))THEN
        NU(NY,NX)=NUI(NY,NX)
        NUM(NY,NX)=NUI(NY,NX)
        DDLYRX(NN)=(VOLWD(NY,NX)-XVOLWP)/AREA(3,0,NY,NX)
        IFLGL(L,NN)=1
        DLYR0=(AZMAX1(VOLW(0,NY,NX)+VOLI(0,NY,NX)-VOLWRX(NY,NX)) &
          +VOLR(NY,NX))/AREA(3,0,NY,NX)
        DLYR(3,0,NY,NX)=DLYR0+DDLYRX(NN)
        DLYR(3,NU(NY,NX),NY,NX)=DLYR(3,NU(NY,NX),NY,NX)-DDLYRX(NN)
        IF(L.GT.2)THEN
          DO LL=L-2,NU(NY,NX),-1
            CDPTH(LL,NY,NX)=CDPTH(L-1,NY,NX)
            CDPTHY(LL,NY,NX)=CDPTHY(L-1,NY,NX)
          ENDDO
        ENDIF
        CDPTH(0,NY,NX)=CDPTH(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
        CDPTHY(0,NY,NX)=CDPTHY(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
        DPTH(NU(NY,NX),NY,NX)=0.5*(CDPTH(NU(NY,NX),NY,NX)+CDPTH(0,NY,NX))
        CDPTHZ(NU(NY,NX),NY,NX)=DLYR(3,NU(NY,NX),NY,NX)
        DPTHZ(NU(NY,NX),NY,NX)=0.5*CDPTHZ(NU(NY,NX),NY,NX)
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
  integer :: N,M,NZ,K,NGL,NR
  real(r8) :: ENGY0,ENGY1
! begin_execution

  IF(L0.NE.0)THEN
    SAND(L1,NY,NX)=SAND(L1,NY,NX)+FX*SAND(L0,NY,NX)
    SILT(L1,NY,NX)=SILT(L1,NY,NX)+FX*SILT(L0,NY,NX)
    CLAY(L1,NY,NX)=CLAY(L1,NY,NX)+FX*CLAY(L0,NY,NX)
    XCEC(L1,NY,NX)=XCEC(L1,NY,NX)+FX*XCEC(L0,NY,NX)
    XAEC(L1,NY,NX)=XAEC(L1,NY,NX)+FX*XAEC(L0,NY,NX)
  ENDIF

  VOLW(L1,NY,NX)=VOLW(L1,NY,NX)+FX*VOLW(L0,NY,NX)
  VOLI(L1,NY,NX)=VOLI(L1,NY,NX)+FX*VOLI(L0,NY,NX)
  VOLP(L1,NY,NX)=VOLP(L1,NY,NX)+FX*VOLP(L0,NY,NX)
  VOLA(L1,NY,NX)=VOLA(L1,NY,NX)+FX*VOLA(L0,NY,NX)
  VOLY(L1,NY,NX)=VOLY(L1,NY,NX)+FX*VOLY(L0,NY,NX)
  VOLWX(L1,NY,NX)=VOLW(L1,NY,NX)
  ENGY1=VHCP(L1,NY,NX)*TKS(L1,NY,NX)
  ENGY0=VHCP(L0,NY,NX)*TKS(L0,NY,NX)
  ENGY1=ENGY1+FX*ENGY0
  VHCM(L1,NY,NX)=VHCM(L1,NY,NX)+FX*VHCM(L0,NY,NX)
  VHCP(L1,NY,NX)=VHCM(L1,NY,NX) &
    +cpw*(VOLW(L1,NY,NX)+VOLWH(L1,NY,NX)) &
    +cpi*(VOLI(L1,NY,NX)+VOLIH(L1,NY,NX))
  IF(VHCP(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L1,NY,NX)=ENGY1/VHCP(L1,NY,NX)
  ELSE
    TKS(L1,NY,NX)=TKS(L0,NY,NX)
  ENDIF
  TCS(L1,NY,NX)=TKS(L1,NY,NX)-TC2K
  ZNH4FA(L1,NY,NX)=ZNH4FA(L1,NY,NX)+FX*ZNH4FA(L0,NY,NX)
  ZNH3FA(L1,NY,NX)=ZNH3FA(L1,NY,NX)+FX*ZNH3FA(L0,NY,NX)
  ZNHUFA(L1,NY,NX)=ZNHUFA(L1,NY,NX)+FX*ZNHUFA(L0,NY,NX)
  ZNO3FA(L1,NY,NX)=ZNO3FA(L1,NY,NX)+FX*ZNO3FA(L0,NY,NX)
  ZNH4FB(L1,NY,NX)=ZNH4FB(L1,NY,NX)+FX*ZNH4FB(L0,NY,NX)
  ZNH3FB(L1,NY,NX)=ZNH3FB(L1,NY,NX)+FX*ZNH3FB(L0,NY,NX)
  ZNHUFB(L1,NY,NX)=ZNHUFB(L1,NY,NX)+FX*ZNHUFB(L0,NY,NX)
  ZNO3FB(L1,NY,NX)=ZNO3FB(L1,NY,NX)+FX*ZNO3FB(L0,NY,NX)
  ZNH4S(L1,NY,NX)=ZNH4S(L1,NY,NX)+FX*ZNH4S(L0,NY,NX)
  ZNH4B(L1,NY,NX)=ZNH4B(L1,NY,NX)+FX*ZNH4B(L0,NY,NX)
  ZNH3S(L1,NY,NX)=ZNH3S(L1,NY,NX)+FX*ZNH3S(L0,NY,NX)
  ZNH3B(L1,NY,NX)=ZNH3B(L1,NY,NX)+FX*ZNH3B(L0,NY,NX)
  ZNO3S(L1,NY,NX)=ZNO3S(L1,NY,NX)+FX*ZNO3S(L0,NY,NX)
  ZNO3B(L1,NY,NX)=ZNO3B(L1,NY,NX)+FX*ZNO3B(L0,NY,NX)
  ZNO2S(L1,NY,NX)=ZNO2S(L1,NY,NX)+FX*ZNO2S(L0,NY,NX)
  ZNO2B(L1,NY,NX)=ZNO2B(L1,NY,NX)+FX*ZNO2B(L0,NY,NX)
  H1PO4(L1,NY,NX)=H1PO4(L1,NY,NX)+FX*H1PO4(L0,NY,NX)
  H2PO4(L1,NY,NX)=H2PO4(L1,NY,NX)+FX*H2PO4(L0,NY,NX)
  IF(ISALTG.NE.0)THEN
    ZAL(L1,NY,NX)=ZAL(L1,NY,NX)+FX*ZAL(L0,NY,NX)
    ZFE(L1,NY,NX)=ZFE(L1,NY,NX)+FX*ZFE(L0,NY,NX)
    ZHY(L1,NY,NX)=ZHY(L1,NY,NX)+FX*ZHY(L0,NY,NX)
    ZCA(L1,NY,NX)=ZCA(L1,NY,NX)+FX*ZCA(L0,NY,NX)
    ZMG(L1,NY,NX)=ZMG(L1,NY,NX)+FX*ZMG(L0,NY,NX)
    ZNA(L1,NY,NX)=ZNA(L1,NY,NX)+FX*ZNA(L0,NY,NX)
    ZKA(L1,NY,NX)=ZKA(L1,NY,NX)+FX*ZKA(L0,NY,NX)
    ZOH(L1,NY,NX)=ZOH(L1,NY,NX)+FX*ZOH(L0,NY,NX)
    ZSO4(L1,NY,NX)=ZSO4(L1,NY,NX)+FX*ZSO4(L0,NY,NX)
    ZCL(L1,NY,NX)=ZCL(L1,NY,NX)+FX*ZCL(L0,NY,NX)
    ZCO3(L1,NY,NX)=ZCO3(L1,NY,NX)+FX*ZCO3(L0,NY,NX)
    ZHCO3(L1,NY,NX)=ZHCO3(L1,NY,NX)+FX*ZHCO3(L0,NY,NX)
    ZALOH1(L1,NY,NX)=ZALOH1(L1,NY,NX)+FX*ZALOH1(L0,NY,NX)
    ZALOH2(L1,NY,NX)=ZALOH2(L1,NY,NX)+FX*ZALOH2(L0,NY,NX)
    ZALOH3(L1,NY,NX)=ZALOH3(L1,NY,NX)+FX*ZALOH3(L0,NY,NX)
    ZALOH4(L1,NY,NX)=ZALOH4(L1,NY,NX)+FX*ZALOH4(L0,NY,NX)
    ZALS(L1,NY,NX)=ZALS(L1,NY,NX)+FX*ZALS(L0,NY,NX)
    ZFEOH1(L1,NY,NX)=ZFEOH1(L1,NY,NX)+FX*ZFEOH1(L0,NY,NX)
    ZFEOH2(L1,NY,NX)=ZFEOH2(L1,NY,NX)+FX*ZFEOH2(L0,NY,NX)
    ZFEOH3(L1,NY,NX)=ZFEOH3(L1,NY,NX)+FX*ZFEOH3(L0,NY,NX)
    ZFEOH4(L1,NY,NX)=ZFEOH4(L1,NY,NX)+FX*ZFEOH4(L0,NY,NX)
    ZFES(L1,NY,NX)=ZFES(L1,NY,NX)+FX*ZFES(L0,NY,NX)
    ZCAO(L1,NY,NX)=ZCAO(L1,NY,NX)+FX*ZCAO(L0,NY,NX)
    ZCAC(L1,NY,NX)=ZCAC(L1,NY,NX)+FX*ZCAC(L0,NY,NX)
    ZCAH(L1,NY,NX)=ZCAH(L1,NY,NX)+FX*ZCAH(L0,NY,NX)
    ZCAS(L1,NY,NX)=ZCAS(L1,NY,NX)+FX*ZCAS(L0,NY,NX)
    ZMGO(L1,NY,NX)=ZMGO(L1,NY,NX)+FX*ZMGO(L0,NY,NX)
    ZMGC(L1,NY,NX)=ZMGC(L1,NY,NX)+FX*ZMGC(L0,NY,NX)
    ZMGH(L1,NY,NX)=ZMGH(L1,NY,NX)+FX*ZMGH(L0,NY,NX)
    ZMGS(L1,NY,NX)=ZMGS(L1,NY,NX)+FX*ZMGS(L0,NY,NX)
    ZNAC(L1,NY,NX)=ZNAC(L1,NY,NX)+FX*ZNAC(L0,NY,NX)
    ZNAS(L1,NY,NX)=ZNAS(L1,NY,NX)+FX*ZNAS(L0,NY,NX)
    ZKAS(L1,NY,NX)=ZKAS(L1,NY,NX)+FX*ZKAS(L0,NY,NX)
    H0PO4(L1,NY,NX)=H0PO4(L1,NY,NX)+FX*H0PO4(L0,NY,NX)
    H3PO4(L1,NY,NX)=H3PO4(L1,NY,NX)+FX*H3PO4(L0,NY,NX)
    ZFE1P(L1,NY,NX)=ZFE1P(L1,NY,NX)+FX*ZFE1P(L0,NY,NX)
    ZFE2P(L1,NY,NX)=ZFE2P(L1,NY,NX)+FX*ZFE2P(L0,NY,NX)
    ZCA0P(L1,NY,NX)=ZCA0P(L1,NY,NX)+FX*ZCA0P(L0,NY,NX)
    ZCA1P(L1,NY,NX)=ZCA1P(L1,NY,NX)+FX*ZCA1P(L0,NY,NX)
    ZCA2P(L1,NY,NX)=ZCA2P(L1,NY,NX)+FX*ZCA2P(L0,NY,NX)
    ZMG1P(L1,NY,NX)=ZMG1P(L1,NY,NX)+FX*ZMG1P(L0,NY,NX)
  ENDIF
  IF(L0.NE.0)THEN
    H1POB(L1,NY,NX)=H1POB(L1,NY,NX)+FX*H1POB(L0,NY,NX)
    H2POB(L1,NY,NX)=H2POB(L1,NY,NX)+FX*H2POB(L0,NY,NX)
    IF(ISALTG.NE.0)THEN
      H0POB(L1,NY,NX)=H0POB(L1,NY,NX)+FX*H0POB(L0,NY,NX)
      H3POB(L1,NY,NX)=H3POB(L1,NY,NX)+FX*H3POB(L0,NY,NX)
      ZFE1PB(L1,NY,NX)=ZFE1PB(L1,NY,NX)+FX*ZFE1PB(L0,NY,NX)
      ZFE2PB(L1,NY,NX)=ZFE2PB(L1,NY,NX)+FX*ZFE2PB(L0,NY,NX)
      ZCA0PB(L1,NY,NX)=ZCA0PB(L1,NY,NX)+FX*ZCA0PB(L0,NY,NX)
      ZCA1PB(L1,NY,NX)=ZCA1PB(L1,NY,NX)+FX*ZCA1PB(L0,NY,NX)
      ZCA2PB(L1,NY,NX)=ZCA2PB(L1,NY,NX)+FX*ZCA2PB(L0,NY,NX)
      ZMG1PB(L1,NY,NX)=ZMG1PB(L1,NY,NX)+FX*ZMG1PB(L0,NY,NX)
    ENDIF
    XN4(L1,NY,NX)=XN4(L1,NY,NX)+FX*XN4(L0,NY,NX)
    XNB(L1,NY,NX)=XNB(L1,NY,NX)+FX*XNB(L0,NY,NX)
    XHY(L1,NY,NX)=XHY(L1,NY,NX)+FX*XHY(L0,NY,NX)
    XAL(L1,NY,NX)=XAL(L1,NY,NX)+FX*XAL(L0,NY,NX)
    XFE(L1,NY,NX)=XFE(L1,NY,NX)+FX*XFE(L0,NY,NX)
    XCA(L1,NY,NX)=XCA(L1,NY,NX)+FX*XCA(L0,NY,NX)
    XMG(L1,NY,NX)=XMG(L1,NY,NX)+FX*XMG(L0,NY,NX)
    XNA(L1,NY,NX)=XNA(L1,NY,NX)+FX*XNA(L0,NY,NX)
    XKA(L1,NY,NX)=XKA(L1,NY,NX)+FX*XKA(L0,NY,NX)
    XHC(L1,NY,NX)=XHC(L1,NY,NX)+FX*XHC(L0,NY,NX)
    XALO2(L1,NY,NX)=XALO2(L1,NY,NX)+FX*XALO2(L0,NY,NX)
    XFEO2(L1,NY,NX)=XFEO2(L1,NY,NX)+FX*XFEO2(L0,NY,NX)
    XOH0(L1,NY,NX)=XOH0(L1,NY,NX)+FX*XOH0(L0,NY,NX)
    XOH1(L1,NY,NX)=XOH1(L1,NY,NX)+FX*XOH1(L0,NY,NX)
    XOH2(L1,NY,NX)=XOH2(L1,NY,NX)+FX*XOH2(L0,NY,NX)
    XH1P(L1,NY,NX)=XH1P(L1,NY,NX)+FX*XH1P(L0,NY,NX)
    XH2P(L1,NY,NX)=XH2P(L1,NY,NX)+FX*XH2P(L0,NY,NX)
    XOH0B(L1,NY,NX)=XOH0B(L1,NY,NX)+FX*XOH0B(L0,NY,NX)
    XOH1B(L1,NY,NX)=XOH1B(L1,NY,NX)+FX*XOH1B(L0,NY,NX)
    XOH2B(L1,NY,NX)=XOH2B(L1,NY,NX)+FX*XOH2B(L0,NY,NX)
    XH1PB(L1,NY,NX)=XH1PB(L1,NY,NX)+FX*XH1PB(L0,NY,NX)
    XH2PB(L1,NY,NX)=XH2PB(L1,NY,NX)+FX*XH2PB(L0,NY,NX)
    PALOH(L1,NY,NX)=PALOH(L1,NY,NX)+FX*PALOH(L0,NY,NX)
    PFEOH(L1,NY,NX)=PFEOH(L1,NY,NX)+FX*PFEOH(L0,NY,NX)
    PCACO(L1,NY,NX)=PCACO(L1,NY,NX)+FX*PCACO(L0,NY,NX)
    PCASO(L1,NY,NX)=PCASO(L1,NY,NX)+FX*PCASO(L0,NY,NX)
    PALPO(L1,NY,NX)=PALPO(L1,NY,NX)+FX*PALPO(L0,NY,NX)
    PFEPO(L1,NY,NX)=PFEPO(L1,NY,NX)+FX*PFEPO(L0,NY,NX)
    PCAPD(L1,NY,NX)=PCAPD(L1,NY,NX)+FX*PCAPD(L0,NY,NX)
    PCAPH(L1,NY,NX)=PCAPH(L1,NY,NX)+FX*PCAPH(L0,NY,NX)
    PCAPM(L1,NY,NX)=PCAPM(L1,NY,NX)+FX*PCAPM(L0,NY,NX)
    PALPB(L1,NY,NX)=PALPB(L1,NY,NX)+FX*PALPB(L0,NY,NX)
    PFEPB(L1,NY,NX)=PFEPB(L1,NY,NX)+FX*PFEPB(L0,NY,NX)
    PCPDB(L1,NY,NX)=PCPDB(L1,NY,NX)+FX*PCPDB(L0,NY,NX)
    PCPHB(L1,NY,NX)=PCPHB(L1,NY,NX)+FX*PCPHB(L0,NY,NX)
    PCPMB(L1,NY,NX)=PCPMB(L1,NY,NX)+FX*PCPMB(L0,NY,NX)
    CO2G(L1,NY,NX)=CO2G(L1,NY,NX)+FX*CO2G(L0,NY,NX)
    CH4G(L1,NY,NX)=CH4G(L1,NY,NX)+FX*CH4G(L0,NY,NX)
    OXYG(L1,NY,NX)=OXYG(L1,NY,NX)+FX*OXYG(L0,NY,NX)
    Z2GG(L1,NY,NX)=Z2GG(L1,NY,NX)+FX*Z2GG(L0,NY,NX)
    Z2OG(L1,NY,NX)=Z2OG(L1,NY,NX)+FX*Z2OG(L0,NY,NX)
    ZNH3G(L1,NY,NX)=ZNH3G(L1,NY,NX)+FX*ZNH3G(L0,NY,NX)
    H2GG(L1,NY,NX)=H2GG(L1,NY,NX)+FX*H2GG(L0,NY,NX)
  ENDIF
  CO2S(L1,NY,NX)=CO2S(L1,NY,NX)+FX*CO2S(L0,NY,NX)
  CH4S(L1,NY,NX)=CH4S(L1,NY,NX)+FX*CH4S(L0,NY,NX)
  OXYS(L1,NY,NX)=OXYS(L1,NY,NX)+FX*OXYS(L0,NY,NX)
  Z2GS(L1,NY,NX)=Z2GS(L1,NY,NX)+FX*Z2GS(L0,NY,NX)
  Z2OS(L1,NY,NX)=Z2OS(L1,NY,NX)+FX*Z2OS(L0,NY,NX)
  H2GS(L1,NY,NX)=H2GS(L1,NY,NX)+FX*H2GS(L0,NY,NX)
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
      IF(WTRTL(1,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.WTRTL(1,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        DO N=1,MY(NZ,NY,NX)
          CO2A(N,L1,NZ,NY,NX)=CO2A(N,L1,NZ,NY,NX)+FX*CO2A(N,L0,NZ,NY,NX)
          OXYA(N,L1,NZ,NY,NX)=OXYA(N,L1,NZ,NY,NX)+FX*OXYA(N,L0,NZ,NY,NX)
          CH4A(N,L1,NZ,NY,NX)=CH4A(N,L1,NZ,NY,NX)+FX*CH4A(N,L0,NZ,NY,NX)
          Z2OA(N,L1,NZ,NY,NX)=Z2OA(N,L1,NZ,NY,NX)+FX*Z2OA(N,L0,NZ,NY,NX)
          ZH3A(N,L1,NZ,NY,NX)=ZH3A(N,L1,NZ,NY,NX)+FX*ZH3A(N,L0,NZ,NY,NX)
          H2GA(N,L1,NZ,NY,NX)=H2GA(N,L1,NZ,NY,NX)+FX*H2GA(N,L0,NZ,NY,NX)
          CO2P(N,L1,NZ,NY,NX)=CO2P(N,L1,NZ,NY,NX)+FX*CO2P(N,L0,NZ,NY,NX)
          OXYP(N,L1,NZ,NY,NX)=OXYP(N,L1,NZ,NY,NX)+FX*OXYP(N,L0,NZ,NY,NX)
          CH4P(N,L1,NZ,NY,NX)=CH4P(N,L1,NZ,NY,NX)+FX*CH4P(N,L0,NZ,NY,NX)
          Z2OP(N,L1,NZ,NY,NX)=Z2OP(N,L1,NZ,NY,NX)+FX*Z2OP(N,L0,NZ,NY,NX)
          ZH3P(N,L1,NZ,NY,NX)=ZH3P(N,L1,NZ,NY,NX)+FX*ZH3P(N,L0,NZ,NY,NX)
          H2GP(N,L1,NZ,NY,NX)=H2GP(N,L1,NZ,NY,NX)+FX*H2GP(N,L0,NZ,NY,NX)
          DO  NR=1,NRT(NZ,NY,NX)
            WTRT1E(ielmc,N,L1,NR,NZ,NY,NX)=WTRT1E(ielmc,N,L1,NR,NZ,NY,NX)+FX*WTRT1E(ielmc,N,L0,NR,NZ,NY,NX)
            WTRT1E(ielmn,N,L1,NR,NZ,NY,NX)=WTRT1E(ielmn,N,L1,NR,NZ,NY,NX)+FX*WTRT1E(ielmn,N,L0,NR,NZ,NY,NX)
            WTRT1E(ielmp,N,L1,NR,NZ,NY,NX)=WTRT1E(ielmp,N,L1,NR,NZ,NY,NX)+FX*WTRT1E(ielmp,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmc,N,L1,NR,NZ,NY,NX)=WTRT2E(ielmc,N,L1,NR,NZ,NY,NX)+FX*WTRT2E(ielmc,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmn,N,L1,NR,NZ,NY,NX)=WTRT2E(ielmn,N,L1,NR,NZ,NY,NX)+FX*WTRT2E(ielmn,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmp,N,L1,NR,NZ,NY,NX)=WTRT2E(ielmp,N,L1,NR,NZ,NY,NX)+FX*WTRT2E(ielmp,N,L0,NR,NZ,NY,NX)
            RTLG1(N,L1,NR,NZ,NY,NX)=RTLG1(N,L1,NR,NZ,NY,NX)+FX*RTLG1(N,L0,NR,NZ,NY,NX)
            RTLG2(N,L1,NR,NZ,NY,NX)=RTLG2(N,L1,NR,NZ,NY,NX)+FX*RTLG2(N,L0,NR,NZ,NY,NX)
            RTN2(N,L1,NR,NZ,NY,NX)=RTN2(N,L1,NR,NZ,NY,NX)+FX*RTN2(N,L0,NR,NZ,NY,NX)
          ENDDO
          EPOOLR(ielmc,N,L1,NZ,NY,NX)=EPOOLR(ielmc,N,L1,NZ,NY,NX)+FX*EPOOLR(ielmc,N,L0,NZ,NY,NX)
          EPOOLR(ielmn,N,L1,NZ,NY,NX)=EPOOLR(ielmn,N,L1,NZ,NY,NX)+FX*EPOOLR(ielmn,N,L0,NZ,NY,NX)
          EPOOLR(ielmp,N,L1,NZ,NY,NX)=EPOOLR(ielmp,N,L1,NZ,NY,NX)+FX*EPOOLR(ielmp,N,L0,NZ,NY,NX)
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
        WTNDLE(L1,ielmc,NZ,NY,NX)=WTNDLE(L1,ielmc,NZ,NY,NX)+FX*WTNDLE(L0,ielmc,NZ,NY,NX)
        WTNDLE(L1,ielmn,NZ,NY,NX)=WTNDLE(L1,ielmn,NZ,NY,NX)+FX*WTNDLE(L0,ielmn,NZ,NY,NX)
        WTNDLE(L1,ielmp,NZ,NY,NX)=WTNDLE(L1,ielmp,NZ,NY,NX)+FX*WTNDLE(L0,ielmp,NZ,NY,NX)
        EPOOLN(L1,ielmc,NZ,NY,NX)=EPOOLN(L1,ielmc,NZ,NY,NX)+FX*EPOOLN(L0,ielmc,NZ,NY,NX)
        EPOOLN(L1,ielmn,NZ,NY,NX)=EPOOLN(L1,ielmn,NZ,NY,NX)+FX*EPOOLN(L0,ielmn,NZ,NY,NX)
        EPOOLN(L1,ielmp,NZ,NY,NX)=EPOOLN(L1,ielmp,NZ,NY,NX)+FX*EPOOLN(L0,ielmp,NZ,NY,NX)
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
    XCEC(L0,NY,NX)=FY*XCEC(L0,NY,NX)
    XAEC(L0,NY,NX)=FY*XAEC(L0,NY,NX)
  ENDIF
!     IF(BKDS(L0,NY,NX).LE.ZERO)THEN
!     VOLT(L0,NY,NX)=FY*VOLT(L0,NY,NX)
!     VOLX(L0,NY,NX)=FY*VOLX(L0,NY,NX)
!     ENDIF
  VOLW(L0,NY,NX)=FY*VOLW(L0,NY,NX)
  VOLI(L0,NY,NX)=FY*VOLI(L0,NY,NX)
  VOLP(L0,NY,NX)=FY*VOLP(L0,NY,NX)
  VOLA(L0,NY,NX)=FY*VOLA(L0,NY,NX)
  VOLY(L0,NY,NX)=FY*VOLY(L0,NY,NX)
  VOLWX(L0,NY,NX)=VOLW(L0,NY,NX)
  ENGY0=FY*ENGY0
  VHCM(L0,NY,NX)=FY*VHCM(L0,NY,NX)
  IF(L0.NE.0)THEN
    VHCP(L0,NY,NX)=VHCM(L0,NY,NX) &
      +cpw*(VOLW(L0,NY,NX)+VOLWH(L0,NY,NX)) &
      +cpi*(VOLI(L0,NY,NX)+VOLIH(L0,NY,NX))
  ELSE
    VHCP(L0,NY,NX)=VHCM(L0,NY,NX) &
      +cpw*VOLW(L0,NY,NX)+cpi*VOLI(L0,NY,NX)
  ENDIF
  IF(VHCP(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L0,NY,NX)=ENGY0/VHCP(L0,NY,NX)
  ELSE
    TKS(L0,NY,NX)=TKS(L1,NY,NX)
  ENDIF
  TCS(L0,NY,NX)=TKS(L0,NY,NX)-TC2K
  ZNH4FA(L0,NY,NX)=FY*ZNH4FA(L0,NY,NX)
  ZNH3FA(L0,NY,NX)=FY*ZNH3FA(L0,NY,NX)
  ZNHUFA(L0,NY,NX)=FY*ZNHUFA(L0,NY,NX)
  ZNO3FA(L0,NY,NX)=FY*ZNO3FA(L0,NY,NX)
  ZNH4FB(L0,NY,NX)=FY*ZNH4FB(L0,NY,NX)
  ZNH3FB(L0,NY,NX)=FY*ZNH3FB(L0,NY,NX)
  ZNHUFB(L0,NY,NX)=FY*ZNHUFB(L0,NY,NX)
  ZNO3FB(L0,NY,NX)=FY*ZNO3FB(L0,NY,NX)
  ZNH4S(L0,NY,NX)=FY*ZNH4S(L0,NY,NX)
  ZNH4B(L0,NY,NX)=FY*ZNH4B(L0,NY,NX)
  ZNH3S(L0,NY,NX)=FY*ZNH3S(L0,NY,NX)
  ZNH3B(L0,NY,NX)=FY*ZNH3B(L0,NY,NX)
  ZNO3S(L0,NY,NX)=FY*ZNO3S(L0,NY,NX)
  ZNO3B(L0,NY,NX)=FY*ZNO3B(L0,NY,NX)
  ZNO2S(L0,NY,NX)=FY*ZNO2S(L0,NY,NX)
  ZNO2B(L0,NY,NX)=FY*ZNO2B(L0,NY,NX)
  H1PO4(L0,NY,NX)=FY*H1PO4(L0,NY,NX)
  H2PO4(L0,NY,NX)=FY*H2PO4(L0,NY,NX)
  IF(ISALTG.NE.0)THEN
    ZAL(L0,NY,NX)=FY*ZAL(L0,NY,NX)
    ZFE(L0,NY,NX)=FY*ZFE(L0,NY,NX)
    ZHY(L0,NY,NX)=FY*ZHY(L0,NY,NX)
    ZCA(L0,NY,NX)=FY*ZCA(L0,NY,NX)
    ZMG(L0,NY,NX)=FY*ZMG(L0,NY,NX)
    ZNA(L0,NY,NX)=FY*ZNA(L0,NY,NX)
    ZKA(L0,NY,NX)=FY*ZKA(L0,NY,NX)
    ZOH(L0,NY,NX)=FY*ZOH(L0,NY,NX)
    ZSO4(L0,NY,NX)=FY*ZSO4(L0,NY,NX)
    ZCL(L0,NY,NX)=FY*ZCL(L0,NY,NX)
    ZCO3(L0,NY,NX)=FY*ZCO3(L0,NY,NX)
    ZHCO3(L0,NY,NX)=FY*ZHCO3(L0,NY,NX)
    ZALOH1(L0,NY,NX)=FY*ZALOH1(L0,NY,NX)
    ZALOH2(L0,NY,NX)=FY*ZALOH2(L0,NY,NX)
    ZALOH3(L0,NY,NX)=FY*ZALOH3(L0,NY,NX)
    ZALOH4(L0,NY,NX)=FY*ZALOH4(L0,NY,NX)
    ZALS(L0,NY,NX)=FY*ZALS(L0,NY,NX)
    ZFEOH1(L0,NY,NX)=FY*ZFEOH1(L0,NY,NX)
    ZFEOH2(L0,NY,NX)=FY*ZFEOH2(L0,NY,NX)
    ZFEOH3(L0,NY,NX)=FY*ZFEOH3(L0,NY,NX)
    ZFEOH4(L0,NY,NX)=FY*ZFEOH4(L0,NY,NX)
    ZFES(L0,NY,NX)=FY*ZFES(L0,NY,NX)
    ZCAO(L0,NY,NX)=FY*ZCAO(L0,NY,NX)
    ZCAC(L0,NY,NX)=FY*ZCAC(L0,NY,NX)
    ZCAH(L0,NY,NX)=FY*ZCAH(L0,NY,NX)
    ZCAS(L0,NY,NX)=FY*ZCAS(L0,NY,NX)
    ZMGO(L0,NY,NX)=FY*ZMGO(L0,NY,NX)
    ZMGC(L0,NY,NX)=FY*ZMGC(L0,NY,NX)
    ZMGH(L0,NY,NX)=FY*ZMGH(L0,NY,NX)
    ZMGS(L0,NY,NX)=FY*ZMGS(L0,NY,NX)
    ZNAC(L0,NY,NX)=FY*ZNAC(L0,NY,NX)
    ZNAS(L0,NY,NX)=FY*ZNAS(L0,NY,NX)
    ZKAS(L0,NY,NX)=FY*ZKAS(L0,NY,NX)
    H0PO4(L0,NY,NX)=FY*H0PO4(L0,NY,NX)
    H3PO4(L0,NY,NX)=FY*H3PO4(L0,NY,NX)
    ZFE1P(L0,NY,NX)=FY*ZFE1P(L0,NY,NX)
    ZFE2P(L0,NY,NX)=FY*ZFE2P(L0,NY,NX)
    ZCA0P(L0,NY,NX)=FY*ZCA0P(L0,NY,NX)
    ZCA1P(L0,NY,NX)=FY*ZCA1P(L0,NY,NX)
    ZCA2P(L0,NY,NX)=FY*ZCA2P(L0,NY,NX)
    ZMG1P(L0,NY,NX)=FY*ZMG1P(L0,NY,NX)
  ENDIF
  IF(L0.NE.0)THEN
    H1POB(L0,NY,NX)=FY*H1POB(L0,NY,NX)
    H2POB(L0,NY,NX)=FY*H2POB(L0,NY,NX)
    IF(ISALTG.NE.0)THEN
      H0POB(L0,NY,NX)=FY*H0POB(L0,NY,NX)
      H3POB(L0,NY,NX)=FY*H3POB(L0,NY,NX)
      ZFE1PB(L0,NY,NX)=FY*ZFE1PB(L0,NY,NX)
      ZFE2PB(L0,NY,NX)=FY*ZFE2PB(L0,NY,NX)
      ZCA0PB(L0,NY,NX)=FY*ZCA0PB(L0,NY,NX)
      ZCA1PB(L0,NY,NX)=FY*ZCA1PB(L0,NY,NX)
      ZCA2PB(L0,NY,NX)=FY*ZCA2PB(L0,NY,NX)
      ZMG1PB(L0,NY,NX)=FY*ZMG1PB(L0,NY,NX)
    ENDIF
    XN4(L0,NY,NX)=FY*XN4(L0,NY,NX)
    XNB(L0,NY,NX)=FY*XNB(L0,NY,NX)
    XHY(L0,NY,NX)=FY*XHY(L0,NY,NX)
    XAL(L0,NY,NX)=FY*XAL(L0,NY,NX)
    XFE(L0,NY,NX)=FY*XFE(L0,NY,NX)
    XCA(L0,NY,NX)=FY*XCA(L0,NY,NX)
    XMG(L0,NY,NX)=FY*XMG(L0,NY,NX)
    XNA(L0,NY,NX)=FY*XNA(L0,NY,NX)
    XKA(L0,NY,NX)=FY*XKA(L0,NY,NX)
    XHC(L0,NY,NX)=FY*XHC(L0,NY,NX)
    XALO2(L0,NY,NX)=FY*XALO2(L0,NY,NX)
    XFEO2(L0,NY,NX)=FY*XFEO2(L0,NY,NX)
    XOH0(L0,NY,NX)=FY*XOH0(L0,NY,NX)
    XOH1(L0,NY,NX)=FY*XOH1(L0,NY,NX)
    XOH2(L0,NY,NX)=FY*XOH2(L0,NY,NX)
    XH1P(L0,NY,NX)=FY*XH1P(L0,NY,NX)
    XH2P(L0,NY,NX)=FY*XH2P(L0,NY,NX)
    XOH0B(L0,NY,NX)=FY*XOH0B(L0,NY,NX)
    XOH1B(L0,NY,NX)=FY*XOH1B(L0,NY,NX)
    XOH2B(L0,NY,NX)=FY*XOH2B(L0,NY,NX)
    XH1PB(L0,NY,NX)=FY*XH1PB(L0,NY,NX)
    XH2PB(L0,NY,NX)=FY*XH2PB(L0,NY,NX)
    PALOH(L0,NY,NX)=FY*PALOH(L0,NY,NX)
    PFEOH(L0,NY,NX)=FY*PFEOH(L0,NY,NX)
    PCACO(L0,NY,NX)=FY*PCACO(L0,NY,NX)
    PCASO(L0,NY,NX)=FY*PCASO(L0,NY,NX)
    PALPO(L0,NY,NX)=FY*PALPO(L0,NY,NX)
    PFEPO(L0,NY,NX)=FY*PFEPO(L0,NY,NX)
    PCAPD(L0,NY,NX)=FY*PCAPD(L0,NY,NX)
    PCAPH(L0,NY,NX)=FY*PCAPH(L0,NY,NX)
    PCAPM(L0,NY,NX)=FY*PCAPM(L0,NY,NX)
    PALPB(L0,NY,NX)=FY*PALPB(L0,NY,NX)
    PFEPB(L0,NY,NX)=FY*PFEPB(L0,NY,NX)
    PCPDB(L0,NY,NX)=FY*PCPDB(L0,NY,NX)
    PCPHB(L0,NY,NX)=FY*PCPHB(L0,NY,NX)
    PCPMB(L0,NY,NX)=FY*PCPMB(L0,NY,NX)
    CO2G(L0,NY,NX)=FY*CO2G(L0,NY,NX)
    CH4G(L0,NY,NX)=FY*CH4G(L0,NY,NX)
    OXYG(L0,NY,NX)=FY*OXYG(L0,NY,NX)
    Z2GG(L0,NY,NX)=FY*Z2GG(L0,NY,NX)
    Z2OG(L0,NY,NX)=FY*Z2OG(L0,NY,NX)
    ZNH3G(L0,NY,NX)=FY*ZNH3G(L0,NY,NX)
    H2GG(L0,NY,NX)=FY*H2GG(L0,NY,NX)
  ENDIF
  CO2S(L0,NY,NX)=FY*CO2S(L0,NY,NX)
  CH4S(L0,NY,NX)=FY*CH4S(L0,NY,NX)
  OXYS(L0,NY,NX)=FY*OXYS(L0,NY,NX)
  Z2GS(L0,NY,NX)=FY*Z2GS(L0,NY,NX)
  Z2OS(L0,NY,NX)=FY*Z2OS(L0,NY,NX)
  H2GS(L0,NY,NX)=FY*H2GS(L0,NY,NX)
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
      IF(WTRTL(1,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.WTRTL(1,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        DO  N=1,MY(NZ,NY,NX)
          CO2A(N,L0,NZ,NY,NX)=FY*CO2A(N,L0,NZ,NY,NX)
          OXYA(N,L0,NZ,NY,NX)=FY*OXYA(N,L0,NZ,NY,NX)
          CH4A(N,L0,NZ,NY,NX)=FY*CH4A(N,L0,NZ,NY,NX)
          Z2OA(N,L0,NZ,NY,NX)=FY*Z2OA(N,L0,NZ,NY,NX)
          ZH3A(N,L0,NZ,NY,NX)=FY*ZH3A(N,L0,NZ,NY,NX)
          H2GA(N,L0,NZ,NY,NX)=FY*H2GA(N,L0,NZ,NY,NX)
          CO2P(N,L0,NZ,NY,NX)=FY*CO2P(N,L0,NZ,NY,NX)
          OXYP(N,L0,NZ,NY,NX)=FY*OXYP(N,L0,NZ,NY,NX)
          CH4P(N,L0,NZ,NY,NX)=FY*CH4P(N,L0,NZ,NY,NX)
          Z2OP(N,L0,NZ,NY,NX)=FY*Z2OP(N,L0,NZ,NY,NX)
          ZH3P(N,L0,NZ,NY,NX)=FY*ZH3P(N,L0,NZ,NY,NX)
          H2GP(N,L0,NZ,NY,NX)=FY*H2GP(N,L0,NZ,NY,NX)
          DO NR=1,NRT(NZ,NY,NX)
            WTRT1E(ielmc,N,L0,NR,NZ,NY,NX)=FY*WTRT1E(ielmc,N,L0,NR,NZ,NY,NX)
            WTRT1E(ielmn,N,L0,NR,NZ,NY,NX)=FY*WTRT1E(ielmn,N,L0,NR,NZ,NY,NX)
            WTRT1E(ielmp,N,L0,NR,NZ,NY,NX)=FY*WTRT1E(ielmp,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmc,N,L0,NR,NZ,NY,NX)=FY*WTRT2E(ielmc,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmn,N,L0,NR,NZ,NY,NX)=FY*WTRT2E(ielmn,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmp,N,L0,NR,NZ,NY,NX)=FY*WTRT2E(ielmp,N,L0,NR,NZ,NY,NX)
            RTLG1(N,L0,NR,NZ,NY,NX)=FY*RTLG1(N,L0,NR,NZ,NY,NX)
            RTLG2(N,L0,NR,NZ,NY,NX)=FY*RTLG2(N,L0,NR,NZ,NY,NX)
            RTN2(N,L0,NR,NZ,NY,NX)=FY*RTN2(N,L0,NR,NZ,NY,NX)
          ENDDO
          EPOOLR(ielmc,N,L0,NZ,NY,NX)=FY*EPOOLR(ielmc,N,L0,NZ,NY,NX)
          EPOOLR(ielmn,N,L0,NZ,NY,NX)=FY*EPOOLR(ielmn,N,L0,NZ,NY,NX)
          EPOOLR(ielmp,N,L0,NZ,NY,NX)=FY*EPOOLR(ielmp,N,L0,NZ,NY,NX)
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
        WTNDLE(L0,ielmc,NZ,NY,NX)=FY*WTNDLE(L0,ielmc,NZ,NY,NX)
        WTNDLE(L0,ielmn,NZ,NY,NX)=FY*WTNDLE(L0,ielmn,NZ,NY,NX)
        WTNDLE(L0,ielmp,NZ,NY,NX)=FY*WTNDLE(L0,ielmp,NZ,NY,NX)
        EPOOLN(L0,ielmc,NZ,NY,NX)=FY*EPOOLN(L0,ielmc,NZ,NY,NX)
        EPOOLN(L0,ielmn,NZ,NY,NX)=FY*EPOOLN(L0,ielmn,NZ,NY,NX)
        EPOOLN(L0,ielmp,NZ,NY,NX)=FY*EPOOLN(L0,ielmp,NZ,NY,NX)
      ENDIF
    ENDDO
  ENDIF
  IF(NN.EQ.1)THEN
    IF(BKDS(L0,NY,NX).LE.ZERO.AND.BKDS(L1,NY,NX).LE.ZERO &
      .AND.VOLW(L0,NY,NX)+VOLI(L0,NY,NX).LE.ZEROS(NY,NX))THEN
      CDPTH(L1,NY,NX)=CDPTH(L0,NY,NX)
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

  integer :: K,N,M,NGL,NR,NZ
  real(r8) :: FXO,FRO
  real(r8) :: FXRTLG2,FXRTN2,FXCPOOLR,FXZPOOLR,FXPPOOLR,FXWTRTL
  real(r8) :: FXWTNDL,FXWTNDLN,FXWTNDLP,FXCPOOLN,FXZPOOLN,FXPPOOLN
  real(r8) :: FXCO2P,FXOXYP,FXCH4P,FXZ2OP,FXZH3P,FXH2GP,FXWTRT1
  real(r8) :: FXWTR1N,FXWTR1P,FXWTRT2,FXWTR2N,FXWTR2P,FXRTLG1
  real(r8) :: FXWTRTD,FXWSRTL,FXRTN1,FXRTNL,FXRTLGP,FXRTDNP
  real(r8) :: FXRTVLP,FXRTVLW,FXRRAD1,FXRRAD2,FXRTARP,FXRTLGA
  real(r8) :: FXCO2A,FXOXYA,FXCH4A,FXZ2OA,FXZH3A,FXH2GA
  real(r8) :: FXOQN,FXOQP,FXOQA,FXOQCH,FXOQNH,FXOQPH,FXOQAH
  real(r8) :: FXOHC,FXOHN,FXOHP,FXOHA,FXOSC,FXOSA,FXOSN,FXOSP
  real(r8) :: FXOMC,FXOMN,FXOMP,FXORC,FXORN,FXORP,FXOQC

! begin_execution
  IF(IFLGL(L,3).EQ.0.AND.L0.NE.0 &
    .AND.VOLX(L0,NY,NX).GT.ZEROS(NY,NX) &
    .AND.VOLX(L1,NY,NX).GT.ZEROS(NY,NX))THEN
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
      IF(FHOL(L1,NY,NX).GT.ZERO.AND.FHOL(L0,NY,NX).GT.ZERO)THEN
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
      IF(WTRTL(1,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX) &
        .AND.WTRTL(1,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
        IF(L0.EQ.L.OR.DPTHZ(L1,NY,NX).LE.ZERO)THEN
          FRO=FO
        ELSE
          FRO=AMIN1(0.5,FO*DPTHZ(L0,NY,NX)/DPTHZ(L1,NY,NX))
        ENDIF
        DO  N=1,MY(NZ,NY,NX)
          FXCO2A=FRO*CO2A(N,L0,NZ,NY,NX)
          CO2A(N,L1,NZ,NY,NX)=CO2A(N,L1,NZ,NY,NX)+FXCO2A
          CO2A(N,L0,NZ,NY,NX)=CO2A(N,L0,NZ,NY,NX)-FXCO2A
          FXOXYA=FRO*OXYA(N,L0,NZ,NY,NX)
          OXYA(N,L1,NZ,NY,NX)=OXYA(N,L1,NZ,NY,NX)+FXOXYA
          OXYA(N,L0,NZ,NY,NX)=OXYA(N,L0,NZ,NY,NX)-FXOXYA
          FXCH4A=FRO*CH4A(N,L0,NZ,NY,NX)
          CH4A(N,L1,NZ,NY,NX)=CH4A(N,L1,NZ,NY,NX)+FXCH4A
          CH4A(N,L0,NZ,NY,NX)=CH4A(N,L0,NZ,NY,NX)-FXCH4A
          FXZ2OA=FRO*Z2OA(N,L0,NZ,NY,NX)
          Z2OA(N,L1,NZ,NY,NX)=Z2OA(N,L1,NZ,NY,NX)+FXZ2OA
          Z2OA(N,L0,NZ,NY,NX)=Z2OA(N,L0,NZ,NY,NX)-FXZ2OA
          FXZH3A=FRO*ZH3A(N,L0,NZ,NY,NX)
          ZH3A(N,L1,NZ,NY,NX)=ZH3A(N,L1,NZ,NY,NX)+FXZH3A
          ZH3A(N,L0,NZ,NY,NX)=ZH3A(N,L0,NZ,NY,NX)-FXZH3A
          FXH2GA=FRO*H2GA(N,L0,NZ,NY,NX)
          H2GA(N,L1,NZ,NY,NX)=H2GA(N,L1,NZ,NY,NX)+FXH2GA
          H2GA(N,L0,NZ,NY,NX)=H2GA(N,L0,NZ,NY,NX)-FXH2GA
          FXCO2P=FRO*CO2P(N,L0,NZ,NY,NX)
          CO2P(N,L1,NZ,NY,NX)=CO2P(N,L1,NZ,NY,NX)+FXCO2P
          CO2P(N,L0,NZ,NY,NX)=CO2P(N,L0,NZ,NY,NX)-FXCO2P
          FXOXYP=FRO*OXYP(N,L0,NZ,NY,NX)
          OXYP(N,L1,NZ,NY,NX)=OXYP(N,L1,NZ,NY,NX)+FXOXYP
          OXYP(N,L0,NZ,NY,NX)=OXYP(N,L0,NZ,NY,NX)-FXOXYP
          FXCH4P=FRO*CH4P(N,L0,NZ,NY,NX)
          CH4P(N,L1,NZ,NY,NX)=CH4P(N,L1,NZ,NY,NX)+FXCH4P
          CH4P(N,L0,NZ,NY,NX)=CH4P(N,L0,NZ,NY,NX)-FXCH4P
          FXZ2OP=FRO*Z2OP(N,L0,NZ,NY,NX)
          Z2OP(N,L1,NZ,NY,NX)=Z2OP(N,L1,NZ,NY,NX)+FXZ2OP
          Z2OP(N,L0,NZ,NY,NX)=Z2OP(N,L0,NZ,NY,NX)-FXZ2OP
          FXZH3P=FRO*ZH3P(N,L0,NZ,NY,NX)
          ZH3P(N,L1,NZ,NY,NX)=ZH3P(N,L1,NZ,NY,NX)+FXZH3P
          ZH3P(N,L0,NZ,NY,NX)=ZH3P(N,L0,NZ,NY,NX)-FXZH3P
          FXH2GP=FRO*H2GP(N,L0,NZ,NY,NX)
          H2GP(N,L1,NZ,NY,NX)=H2GP(N,L1,NZ,NY,NX)+FXH2GP
          H2GP(N,L0,NZ,NY,NX)=H2GP(N,L0,NZ,NY,NX)-FXH2GP
          DO  NR=1,NRT(NZ,NY,NX)
            FXWTRT1=FRO*WTRT1E(ielmc,N,L0,NR,NZ,NY,NX)
            WTRT1E(ielmc,N,L1,NR,NZ,NY,NX)=WTRT1E(ielmc,N,L1,NR,NZ,NY,NX)+FXWTRT1
            WTRT1E(ielmc,N,L0,NR,NZ,NY,NX)=WTRT1E(ielmc,N,L0,NR,NZ,NY,NX)-FXWTRT1
            FXWTR1N=FRO*WTRT1E(ielmn,N,L0,NR,NZ,NY,NX)
            WTRT1E(ielmn,N,L1,NR,NZ,NY,NX)=WTRT1E(ielmn,N,L1,NR,NZ,NY,NX)+FXWTR1N
            WTRT1E(ielmn,N,L0,NR,NZ,NY,NX)=WTRT1E(ielmn,N,L0,NR,NZ,NY,NX)-FXWTR1N
            FXWTR1P=FRO*WTRT1E(ielmp,N,L0,NR,NZ,NY,NX)
            WTRT1E(ielmp,N,L1,NR,NZ,NY,NX)=WTRT1E(ielmp,N,L1,NR,NZ,NY,NX)+FXWTR1P
            WTRT1E(ielmp,N,L0,NR,NZ,NY,NX)=WTRT1E(ielmp,N,L0,NR,NZ,NY,NX)-FXWTR1P
            FXWTRT2=FRO*WTRT2E(ielmc,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmc,N,L1,NR,NZ,NY,NX)=WTRT2E(ielmc,N,L1,NR,NZ,NY,NX)+FXWTRT2
            WTRT2E(ielmc,N,L0,NR,NZ,NY,NX)=WTRT2E(ielmc,N,L0,NR,NZ,NY,NX)-FXWTRT2
            FXWTR2N=FRO*WTRT2E(ielmn,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmn,N,L1,NR,NZ,NY,NX)=WTRT2E(ielmn,N,L1,NR,NZ,NY,NX)+FXWTR2N
            WTRT2E(ielmn,N,L0,NR,NZ,NY,NX)=WTRT2E(ielmn,N,L0,NR,NZ,NY,NX)-FXWTR2N
            FXWTR2P=FRO*WTRT2E(ielmp,N,L0,NR,NZ,NY,NX)
            WTRT2E(ielmp,N,L1,NR,NZ,NY,NX)=WTRT2E(ielmp,N,L1,NR,NZ,NY,NX)+FXWTR2P
            WTRT2E(ielmp,N,L0,NR,NZ,NY,NX)=WTRT2E(ielmp,N,L0,NR,NZ,NY,NX)-FXWTR2P
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
          FXCPOOLR=FRO*EPOOLR(ielmc,N,L0,NZ,NY,NX)
          EPOOLR(ielmc,N,L1,NZ,NY,NX)=EPOOLR(ielmc,N,L1,NZ,NY,NX)+FXCPOOLR
          EPOOLR(ielmc,N,L0,NZ,NY,NX)=EPOOLR(ielmc,N,L0,NZ,NY,NX)-FXCPOOLR
          FXZPOOLR=FRO*EPOOLR(ielmn,N,L0,NZ,NY,NX)
          EPOOLR(ielmn,N,L1,NZ,NY,NX)=EPOOLR(ielmn,N,L1,NZ,NY,NX)+FXZPOOLR
          EPOOLR(ielmn,N,L0,NZ,NY,NX)=EPOOLR(ielmn,N,L0,NZ,NY,NX)-FXZPOOLR
          FXPPOOLR=FRO*EPOOLR(ielmp,N,L0,NZ,NY,NX)
          EPOOLR(ielmp,N,L1,NZ,NY,NX)=EPOOLR(ielmp,N,L1,NZ,NY,NX)+FXPPOOLR
          EPOOLR(ielmp,N,L0,NZ,NY,NX)=EPOOLR(ielmp,N,L0,NZ,NY,NX)-FXPPOOLR
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
        FXWTNDL=FRO*WTNDLE(L0,ielmc,NZ,NY,NX)
        WTNDLE(L1,ielmc,NZ,NY,NX)=WTNDLE(L1,ielmc,NZ,NY,NX)+FXWTNDL
        WTNDLE(L0,ielmc,NZ,NY,NX)=WTNDLE(L0,ielmc,NZ,NY,NX)-FXWTNDL
        FXWTNDLN=FRO*WTNDLE(L0,ielmn,NZ,NY,NX)
        WTNDLE(L1,ielmn,NZ,NY,NX)=WTNDLE(L1,ielmn,NZ,NY,NX)+FXWTNDLN
        WTNDLE(L0,ielmn,NZ,NY,NX)=WTNDLE(L0,ielmn,NZ,NY,NX)-FXWTNDLN
        FXWTNDLP=FRO*WTNDLE(L0,ielmp,NZ,NY,NX)
        WTNDLE(L1,ielmp,NZ,NY,NX)=WTNDLE(L1,ielmp,NZ,NY,NX)+FXWTNDLP
        WTNDLE(L0,ielmp,NZ,NY,NX)=WTNDLE(L0,ielmp,NZ,NY,NX)-FXWTNDLP
        FXCPOOLN=FRO*EPOOLN(L0,ielmc,NZ,NY,NX)
        EPOOLN(L1,ielmc,NZ,NY,NX)=EPOOLN(L1,ielmc,NZ,NY,NX)+FXCPOOLN
        EPOOLN(L0,ielmc,NZ,NY,NX)=EPOOLN(L0,ielmc,NZ,NY,NX)-FXCPOOLN
        FXZPOOLN=FRO*EPOOLN(L0,ielmn,NZ,NY,NX)
        EPOOLN(L1,ielmn,NZ,NY,NX)=EPOOLN(L1,ielmn,NZ,NY,NX)+FXZPOOLN
        EPOOLN(L0,ielmn,NZ,NY,NX)=EPOOLN(L0,ielmn,NZ,NY,NX)-FXZPOOLN
        FXPPOOLN=FRO*EPOOLN(L0,ielmp,NZ,NY,NX)
        EPOOLN(L1,ielmp,NZ,NY,NX)=EPOOLN(L1,ielmp,NZ,NY,NX)+FXPPOOLN
        EPOOLN(L0,ielmp,NZ,NY,NX)=EPOOLN(L0,ielmp,NZ,NY,NX)-FXPPOOLN
      ENDIF
    ENDDO
  ENDIF
  end subroutine MoveSOM
!------------------------------------------------------------------------------------------

  subroutine MoveMMPoreSolute(L0,L1,NY,NX,FHO)

  implicit none
  integer, intent(in) :: L0,L1,NY,NX
  real(r8), intent(in) :: FHO
  real(r8) :: FXZSO4H
  real(r8) :: FXZNH4SH,FXZNH3SH,FXZNO3SH,FXZNO2SH,FXZNH4BH
  real(r8) :: FXZNH3BH,FXZNO3BH,FXZNO2BH,FXH1PO4H,FXH2PO4H,FXH1POBH
  real(r8) :: FXH2POBH,FXZALH,FXZFEH,FXZHYH,FXZCCH,FXZMAH,FXZNAH
  real(r8) :: FXZKAH,FXZOHH,FXZCLH,FXZCO3H,FXZHCO3H
  real(r8) :: FXZALO1H,FXZALO2H,FXZALO3H,FXZALO4H,FXZALSH,FXZFEO1H
  real(r8) :: FXZNACH,FXZNASH,FXZKASH,FXH0PO4H,FXH3PO4H,FXZFE1PH
  real(r8) :: FXZMG1BH,FXCO2SH,FXCH4SH,FXOXYSH,FXZ2GSH,FXZ2OSH
  real(r8) :: FXZCAHH,FXZCASH,FXZMGOH,FXZMGCH,FXZMGHH,FXZMGSH
  real(r8) :: FXZFE2PH,FXZCA0PH,FXZCA1PH,FXZCA2PH,FXZMG1PH,FXH0POBH
  real(r8) :: FXZFEO2H,FXZFEO3H,FXZFEO4H,FXZFESH,FXZCAOH,FXZCACH
  real(r8) :: FXH3POBH,FXZFE1BH,FXZFE2BH,FXZCA0BH,FXZCA1BH,FXZCA2BH

! begin_execution
!
!     SOIL MACROPORE N,P SOLUTES
!
  IF(FHOL(L1,NY,NX).GT.ZERO.AND.FHOL(L0,NY,NX).GT.ZERO)THEN
    FXZNH4SH=FHO*ZNH4SH(L0,NY,NX)
    ZNH4SH(L1,NY,NX)=ZNH4SH(L1,NY,NX)+FXZNH4SH
    ZNH4SH(L0,NY,NX)=ZNH4SH(L0,NY,NX)-FXZNH4SH
    FXZNH3SH=FHO*ZNH3SH(L0,NY,NX)
    ZNH3SH(L1,NY,NX)=ZNH3SH(L1,NY,NX)+FXZNH3SH
    ZNH3SH(L0,NY,NX)=ZNH3SH(L0,NY,NX)-FXZNH3SH
    FXZNO3SH=FHO*ZNO3SH(L0,NY,NX)
    ZNO3SH(L1,NY,NX)=ZNO3SH(L1,NY,NX)+FXZNO3SH
    ZNO3SH(L0,NY,NX)=ZNO3SH(L0,NY,NX)-FXZNO3SH
    FXZNO2SH=FHO*ZNO2SH(L0,NY,NX)
    ZNO2SH(L1,NY,NX)=ZNO2SH(L1,NY,NX)+FXZNO2SH
    ZNO2SH(L0,NY,NX)=ZNO2SH(L0,NY,NX)-FXZNO2SH
    FXZNH4BH=FHO*ZNH4BH(L0,NY,NX)
    ZNH4BH(L1,NY,NX)=ZNH4BH(L1,NY,NX)+FXZNH4BH
    ZNH4BH(L0,NY,NX)=ZNH4BH(L0,NY,NX)-FXZNH4BH
    FXZNH3BH=FHO*ZNH3BH(L0,NY,NX)
    ZNH3BH(L1,NY,NX)=ZNH3BH(L1,NY,NX)+FXZNH3BH
    ZNH3BH(L0,NY,NX)=ZNH3BH(L0,NY,NX)-FXZNH3BH
    FXZNO3BH=FHO*ZNO3BH(L0,NY,NX)
    ZNO3BH(L1,NY,NX)=ZNO3BH(L1,NY,NX)+FXZNO3BH
    ZNO3BH(L0,NY,NX)=ZNO3BH(L0,NY,NX)-FXZNO3BH
    FXZNO2BH=FHO*ZNO2BH(L0,NY,NX)
    ZNO2BH(L1,NY,NX)=ZNO2BH(L1,NY,NX)+FXZNO2BH
    ZNO2BH(L0,NY,NX)=ZNO2BH(L0,NY,NX)-FXZNO2BH
    FXH1PO4H=FHO*H1PO4H(L0,NY,NX)
    H1PO4H(L1,NY,NX)=H1PO4H(L1,NY,NX)+FXH1PO4H
    H1PO4H(L0,NY,NX)=H1PO4H(L0,NY,NX)-FXH1PO4H
    FXH2PO4H=FHO*H2PO4H(L0,NY,NX)
    H2PO4H(L1,NY,NX)=H2PO4H(L1,NY,NX)+FXH2PO4H
    H2PO4H(L0,NY,NX)=H2PO4H(L0,NY,NX)-FXH2PO4H
    FXH1POBH=FHO*H1POBH(L0,NY,NX)
    H1POBH(L1,NY,NX)=H1POBH(L1,NY,NX)+FXH1POBH
    H1POBH(L0,NY,NX)=H1POBH(L0,NY,NX)-FXH1POBH
    FXH2POBH=FHO*H2POBH(L0,NY,NX)
    H2POBH(L1,NY,NX)=H2POBH(L1,NY,NX)+FXH2POBH
    H2POBH(L0,NY,NX)=H2POBH(L0,NY,NX)-FXH2POBH
!
!     SOIL MACROPORE SOLUBLE SALTS
!
    IF(ISALTG.NE.0)THEN
      FXZALH=FHO*ZALH(L0,NY,NX)
      ZALH(L1,NY,NX)=ZALH(L1,NY,NX)+FXZALH
      ZALH(L0,NY,NX)=ZALH(L0,NY,NX)-FXZALH
      FXZFEH=FHO*ZFEH(L0,NY,NX)
      ZFEH(L1,NY,NX)=ZFEH(L1,NY,NX)+FXZFEH
      ZFEH(L0,NY,NX)=ZFEH(L0,NY,NX)-FXZFEH
      FXZHYH=FHO*ZHYH(L0,NY,NX)
      ZHYH(L1,NY,NX)=ZHYH(L1,NY,NX)+FXZHYH
      ZHYH(L0,NY,NX)=ZHYH(L0,NY,NX)-FXZHYH
      FXZCCH=FHO*ZCCH(L0,NY,NX)
      ZCCH(L1,NY,NX)=ZCCH(L1,NY,NX)+FXZCCH
      ZCCH(L0,NY,NX)=ZCCH(L0,NY,NX)-FXZCCH
      FXZMAH=FHO*ZMAH(L0,NY,NX)
      ZMAH(L1,NY,NX)=ZMAH(L1,NY,NX)+FXZMAH
      ZMAH(L0,NY,NX)=ZMAH(L0,NY,NX)-FXZMAH
      FXZNAH=FHO*ZNAH(L0,NY,NX)
      ZNAH(L1,NY,NX)=ZNAH(L1,NY,NX)+FXZNAH
      ZNAH(L0,NY,NX)=ZNAH(L0,NY,NX)-FXZNAH
      FXZKAH=FHO*ZKAH(L0,NY,NX)
      ZKAH(L1,NY,NX)=ZKAH(L1,NY,NX)+FXZKAH
      ZKAH(L0,NY,NX)=ZKAH(L0,NY,NX)-FXZKAH
      FXZOHH=FHO*ZOHH(L0,NY,NX)
      ZOHH(L1,NY,NX)=ZOHH(L1,NY,NX)+FXZOHH
      ZOHH(L0,NY,NX)=ZOHH(L0,NY,NX)-FXZOHH
      FXZSO4H=FHO*ZSO4H(L0,NY,NX)
      ZSO4H(L1,NY,NX)=ZSO4H(L1,NY,NX)+FXZSO4H
      ZSO4H(L0,NY,NX)=ZSO4H(L0,NY,NX)-FXZSO4H
      FXZCLH=FHO*ZCLH(L0,NY,NX)
      ZCLH(L1,NY,NX)=ZCLH(L1,NY,NX)+FXZCLH
      ZCLH(L0,NY,NX)=ZCLH(L0,NY,NX)-FXZCLH
      FXZCO3H=FHO*ZCO3H(L0,NY,NX)
      ZCO3H(L1,NY,NX)=ZCO3H(L1,NY,NX)+FXZCO3H
      ZCO3H(L0,NY,NX)=ZCO3H(L0,NY,NX)-FXZCO3H
      FXZHCO3H=FHO*ZHCO3H(L0,NY,NX)
      ZHCO3H(L1,NY,NX)=ZHCO3H(L1,NY,NX)+FXZHCO3H
      ZHCO3H(L0,NY,NX)=ZHCO3H(L0,NY,NX)-FXZHCO3H
      FXZALO1H=FHO*ZALO1H(L0,NY,NX)
      ZALO1H(L1,NY,NX)=ZALO1H(L1,NY,NX)+FXZALO1H
      ZALO1H(L0,NY,NX)=ZALO1H(L0,NY,NX)-FXZALO1H
      FXZALO2H=FHO*ZALO2H(L0,NY,NX)
      ZALO2H(L1,NY,NX)=ZALO2H(L1,NY,NX)+FXZALO2H
      ZALO2H(L0,NY,NX)=ZALO2H(L0,NY,NX)-FXZALO2H
      FXZALO3H=FHO*ZALO3H(L0,NY,NX)
      ZALO3H(L1,NY,NX)=ZALO3H(L1,NY,NX)+FXZALO3H
      ZALO3H(L0,NY,NX)=ZALO3H(L0,NY,NX)-FXZALO3H
      FXZALO4H=FHO*ZALO4H(L0,NY,NX)
      ZALO4H(L1,NY,NX)=ZALO4H(L1,NY,NX)+FXZALO4H
      ZALO4H(L0,NY,NX)=ZALO4H(L0,NY,NX)-FXZALO4H
      FXZALSH=FHO*ZALSH(L0,NY,NX)
      ZALSH(L1,NY,NX)=ZALSH(L1,NY,NX)+FXZALSH
      ZALSH(L0,NY,NX)=ZALSH(L0,NY,NX)-FXZALSH
      FXZFEO1H=FHO*ZFEO1H(L0,NY,NX)
      ZFEO1H(L1,NY,NX)=ZFEO1H(L1,NY,NX)+FXZFEO1H
      ZFEO1H(L0,NY,NX)=ZFEO1H(L0,NY,NX)-FXZFEO1H
      FXZFEO2H=FHO*ZFEO2H(L0,NY,NX)
      ZFEO2H(L1,NY,NX)=ZFEO2H(L1,NY,NX)+FXZFEO2H
      ZFEO2H(L0,NY,NX)=ZFEO2H(L0,NY,NX)-FXZFEO2H
      FXZFEO3H=FHO*ZFEO3H(L0,NY,NX)
      ZFEO3H(L1,NY,NX)=ZFEO3H(L1,NY,NX)+FXZFEO3H
      ZFEO3H(L0,NY,NX)=ZFEO3H(L0,NY,NX)-FXZFEO3H
      FXZFEO4H=FHO*ZFEO4H(L0,NY,NX)
      ZFEO4H(L1,NY,NX)=ZFEO4H(L1,NY,NX)+FXZFEO4H
      ZFEO4H(L0,NY,NX)=ZFEO4H(L0,NY,NX)-FXZFEO4H
      FXZFESH=FHO*ZFESH(L0,NY,NX)
      ZFESH(L1,NY,NX)=ZFESH(L1,NY,NX)+FXZFESH
      ZFESH(L0,NY,NX)=ZFESH(L0,NY,NX)-FXZFESH
      FXZCAOH=FHO*ZCAOH(L0,NY,NX)
      ZCAOH(L1,NY,NX)=ZCAOH(L1,NY,NX)+FXZCAOH
      ZCAOH(L0,NY,NX)=ZCAOH(L0,NY,NX)-FXZCAOH
      FXZCACH=FHO*ZCACH(L0,NY,NX)
      ZCACH(L1,NY,NX)=ZCACH(L1,NY,NX)+FXZCACH
      ZCACH(L0,NY,NX)=ZCACH(L0,NY,NX)-FXZCACH
      FXZCAHH=FHO*ZCAHH(L0,NY,NX)
      ZCAHH(L1,NY,NX)=ZCAHH(L1,NY,NX)+FXZCAHH
      ZCAHH(L0,NY,NX)=ZCAHH(L0,NY,NX)-FXZCAHH
      FXZCASH=FHO*ZCASH(L0,NY,NX)
      ZCASH(L1,NY,NX)=ZCASH(L1,NY,NX)+FXZCASH
      ZCASH(L0,NY,NX)=ZCASH(L0,NY,NX)-FXZCASH
      FXZMGOH=FHO*ZMGOH(L0,NY,NX)
      ZMGOH(L1,NY,NX)=ZMGOH(L1,NY,NX)+FXZMGOH
      ZMGOH(L0,NY,NX)=ZMGOH(L0,NY,NX)-FXZMGOH
      FXZMGCH=FHO*ZMGCH(L0,NY,NX)
      ZMGCH(L1,NY,NX)=ZMGCH(L1,NY,NX)+FXZMGCH
      ZMGCH(L0,NY,NX)=ZMGCH(L0,NY,NX)-FXZMGCH
      FXZMGHH=FHO*ZMGHH(L0,NY,NX)
      ZMGHH(L1,NY,NX)=ZMGHH(L1,NY,NX)+FXZMGHH
      ZMGHH(L0,NY,NX)=ZMGHH(L0,NY,NX)-FXZMGHH
      FXZMGSH=FHO*ZMGSH(L0,NY,NX)
      ZMGSH(L1,NY,NX)=ZMGSH(L1,NY,NX)+FXZMGSH
      ZMGSH(L0,NY,NX)=ZMGSH(L0,NY,NX)-FXZMGSH
      FXZNACH=FHO*ZNACH(L0,NY,NX)
      ZNACH(L1,NY,NX)=ZNACH(L1,NY,NX)+FXZNACH
      ZNACH(L0,NY,NX)=ZNACH(L0,NY,NX)-FXZNACH
      FXZNASH=FHO*ZNASH(L0,NY,NX)
      ZNASH(L1,NY,NX)=ZNASH(L1,NY,NX)+FXZNASH
      ZNASH(L0,NY,NX)=ZNASH(L0,NY,NX)-FXZNASH
      FXZKASH=FHO*ZKASH(L0,NY,NX)
      ZKASH(L1,NY,NX)=ZKASH(L1,NY,NX)+FXZKASH
      ZKASH(L0,NY,NX)=ZKASH(L0,NY,NX)-FXZKASH
      FXH0PO4H=FHO*H0PO4H(L0,NY,NX)
      H0PO4H(L1,NY,NX)=H0PO4H(L1,NY,NX)+FXH0PO4H
      H0PO4H(L0,NY,NX)=H0PO4H(L0,NY,NX)-FXH0PO4H
      FXH3PO4H=FHO*H3PO4H(L0,NY,NX)
      H3PO4H(L1,NY,NX)=H3PO4H(L1,NY,NX)+FXH3PO4H
      H3PO4H(L0,NY,NX)=H3PO4H(L0,NY,NX)-FXH3PO4H
      FXZFE1PH=FHO*ZFE1PH(L0,NY,NX)
      ZFE1PH(L1,NY,NX)=ZFE1PH(L1,NY,NX)+FXZFE1PH
      ZFE1PH(L0,NY,NX)=ZFE1PH(L0,NY,NX)-FXZFE1PH
      FXZFE2PH=FHO*ZFE2PH(L0,NY,NX)
      ZFE2PH(L1,NY,NX)=ZFE2PH(L1,NY,NX)+FXZFE2PH
      ZFE2PH(L0,NY,NX)=ZFE2PH(L0,NY,NX)-FXZFE2PH
      FXZCA0PH=FHO*ZCA0PH(L0,NY,NX)
      ZCA0PH(L1,NY,NX)=ZCA0PH(L1,NY,NX)+FXZCA0PH
      ZCA0PH(L0,NY,NX)=ZCA0PH(L0,NY,NX)-FXZCA0PH
      FXZCA1PH=FHO*ZCA1PH(L0,NY,NX)
      ZCA1PH(L1,NY,NX)=ZCA1PH(L1,NY,NX)+FXZCA1PH
      ZCA1PH(L0,NY,NX)=ZCA1PH(L0,NY,NX)-FXZCA1PH
      FXZCA2PH=FHO*ZCA2PH(L0,NY,NX)
      ZCA2PH(L1,NY,NX)=ZCA2PH(L1,NY,NX)+FXZCA2PH
      ZCA2PH(L0,NY,NX)=ZCA2PH(L0,NY,NX)-FXZCA2PH
      FXZMG1PH=FHO*ZMG1PH(L0,NY,NX)
      ZMG1PH(L1,NY,NX)=ZMG1PH(L1,NY,NX)+FXZMG1PH
      ZMG1PH(L0,NY,NX)=ZMG1PH(L0,NY,NX)-FXZMG1PH
      FXH0POBH=FHO*H0POBH(L0,NY,NX)
      H0POBH(L1,NY,NX)=H0POBH(L1,NY,NX)+FXH0POBH
      H0POBH(L0,NY,NX)=H0POBH(L0,NY,NX)-FXH0POBH
      FXH3POBH=FHO*H3POBH(L0,NY,NX)
      H3POBH(L1,NY,NX)=H3POBH(L1,NY,NX)+FXH3POBH
      H3POBH(L0,NY,NX)=H3POBH(L0,NY,NX)-FXH3POBH
      FXZFE1BH=FHO*ZFE1BH(L0,NY,NX)
      ZFE1BH(L1,NY,NX)=ZFE1BH(L1,NY,NX)+FXZFE1BH
      ZFE1BH(L0,NY,NX)=ZFE1BH(L0,NY,NX)-FXZFE1BH
      FXZFE2BH=FHO*ZFE2BH(L0,NY,NX)
      ZFE2BH(L1,NY,NX)=ZFE2BH(L1,NY,NX)+FXZFE2BH
      ZFE2BH(L0,NY,NX)=ZFE2BH(L0,NY,NX)-FXZFE2BH
      FXZCA0BH=FHO*ZCA0BH(L0,NY,NX)
      ZCA0BH(L1,NY,NX)=ZCA0BH(L1,NY,NX)+FXZCA0BH
      ZCA0BH(L0,NY,NX)=ZCA0BH(L0,NY,NX)-FXZCA0BH
      FXZCA1BH=FHO*ZCA1BH(L0,NY,NX)
      ZCA1BH(L1,NY,NX)=ZCA1BH(L1,NY,NX)+FXZCA1BH
      ZCA1BH(L0,NY,NX)=ZCA1BH(L0,NY,NX)-FXZCA1BH
      FXZCA2BH=FHO*ZCA2BH(L0,NY,NX)
      ZCA2BH(L1,NY,NX)=ZCA2BH(L1,NY,NX)+FXZCA2BH
      ZCA2BH(L0,NY,NX)=ZCA2BH(L0,NY,NX)-FXZCA2BH
      FXZMG1BH=FHO*ZMG1BH(L0,NY,NX)
      ZMG1BH(L1,NY,NX)=ZMG1BH(L1,NY,NX)+FXZMG1BH
      ZMG1BH(L0,NY,NX)=ZMG1BH(L0,NY,NX)-FXZMG1BH
    ENDIF
!
!     SOIL MACROPORE AQUEOUS GASES
!
    FXCO2SH=FHO*CO2SH(L0,NY,NX)
    CO2SH(L1,NY,NX)=CO2SH(L1,NY,NX)+FXCO2SH
    CO2SH(L0,NY,NX)=CO2SH(L0,NY,NX)-FXCO2SH
    FXCH4SH=FHO*CH4SH(L0,NY,NX)
    CH4SH(L1,NY,NX)=CH4SH(L1,NY,NX)+FXCH4SH
    CH4SH(L0,NY,NX)=CH4SH(L0,NY,NX)-FXCH4SH
    FXOXYSH=FHO*OXYSH(L0,NY,NX)
    OXYSH(L1,NY,NX)=OXYSH(L1,NY,NX)+FXOXYSH
    OXYSH(L0,NY,NX)=OXYSH(L0,NY,NX)-FXOXYSH
    FXZ2GSH=FHO*Z2GSH(L0,NY,NX)
    Z2GSH(L1,NY,NX)=Z2GSH(L1,NY,NX)+FXZ2GSH
    Z2GSH(L0,NY,NX)=Z2GSH(L0,NY,NX)-FXZ2GSH
    FXZ2OSH=FHO*Z2OSH(L0,NY,NX)
    Z2OSH(L1,NY,NX)=Z2OSH(L1,NY,NX)+FXZ2OSH
    Z2OSH(L0,NY,NX)=Z2OSH(L0,NY,NX)-FXZ2OSH
  ENDIF
  end subroutine MoveMMPoreSolute

!------------------------------------------------------------------------------------------

  subroutine MoveBandSolute(L,L0,L1,NY,NX,FX,FWO,FO)
  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in) :: FX,FWO,FO
  real(r8) :: FAO,FCO
  real(r8) :: FXZCA2PB,FXZMG1PB,FXXCEC,FXXN4,FXXNB,FXXHY
  real(r8) :: FXH0POB,FXH3POB,FXZFE1PB,FXZFE2PB,FXZCA0PB,FXZCA1PB
  real(r8) :: FXXH2P,FXXOH0B,FXXOH1B,FXXOH2B,FXXH1PB,FXXH2PB
  real(r8) :: FXXFEO2,FXXAEC,FXXOH0,FXXOH1,FXXOH2,FXXH1P
  real(r8) :: FXXAL,FXXFE,FXXCA,FXXMG,FXXNA,FXXKA,FXXHC,FXXALO2
  real(r8) :: FXPALOH,FXPFEOH,FXPCACO,FXPCASO,FXPALPO,FXPFEPO
  real(r8) :: FXPCAPD,FXPCAPH,FXPCAPM,FXPALPB,FXPFEPB,FXPCPDB
  real(r8) :: FXPCPHB,FXPCPMB
  real(r8) :: FXH1POB,FXH2POB

  FXH1POB=FWO*H1POB(L0,NY,NX)
  H1POB(L1,NY,NX)=H1POB(L1,NY,NX)+FXH1POB
  H1POB(L0,NY,NX)=H1POB(L0,NY,NX)-FXH1POB
  FXH2POB=FWO*H2POB(L0,NY,NX)
  H2POB(L1,NY,NX)=H2POB(L1,NY,NX)+FXH2POB
  H2POB(L0,NY,NX)=H2POB(L0,NY,NX)-FXH2POB
  IF(ISALTG.NE.0)THEN
    FXH0POB=FWO*H0POB(L0,NY,NX)
    H0POB(L1,NY,NX)=H0POB(L1,NY,NX)+FXH0POB
    H0POB(L0,NY,NX)=H0POB(L0,NY,NX)-FXH0POB
    FXH3POB=FWO*H3POB(L0,NY,NX)
    H3POB(L1,NY,NX)=H3POB(L1,NY,NX)+FXH3POB
    H3POB(L0,NY,NX)=H3POB(L0,NY,NX)-FXH3POB
    FXZFE1PB=FWO*ZFE1PB(L0,NY,NX)
    ZFE1PB(L1,NY,NX)=ZFE1PB(L1,NY,NX)+FXZFE1PB
    ZFE1PB(L0,NY,NX)=ZFE1PB(L0,NY,NX)-FXZFE1PB
    FXZFE2PB=FWO*ZFE2PB(L0,NY,NX)
    ZFE2PB(L1,NY,NX)=ZFE2PB(L1,NY,NX)+FXZFE2PB
    ZFE2PB(L0,NY,NX)=ZFE2PB(L0,NY,NX)-FXZFE2PB
    FXZCA0PB=FWO*ZCA0PB(L0,NY,NX)
    ZCA0PB(L1,NY,NX)=ZCA0PB(L1,NY,NX)+FXZCA0PB
    ZCA0PB(L0,NY,NX)=ZCA0PB(L0,NY,NX)-FXZCA0PB
    FXZCA1PB=FWO*ZCA1PB(L0,NY,NX)
    ZCA1PB(L1,NY,NX)=ZCA1PB(L1,NY,NX)+FXZCA1PB
    ZCA1PB(L0,NY,NX)=ZCA1PB(L0,NY,NX)-FXZCA1PB
    FXZCA2PB=FWO*ZCA2PB(L0,NY,NX)
    ZCA2PB(L1,NY,NX)=ZCA2PB(L1,NY,NX)+FXZCA2PB
    ZCA2PB(L0,NY,NX)=ZCA2PB(L0,NY,NX)-FXZCA2PB
    FXZMG1PB=FWO*ZMG1PB(L0,NY,NX)
    ZMG1PB(L1,NY,NX)=ZMG1PB(L1,NY,NX)+FXZMG1PB
    ZMG1PB(L0,NY,NX)=ZMG1PB(L0,NY,NX)-FXZMG1PB
  ENDIF
!
!     SOIL ADSORBED CATIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L.OR.CEC(L0,NY,NX).LE.ZERO)THEN
    FCO=FO
  ELSE
    FCO=AMIN1(0.5,FO*CEC(L1,NY,NX)/CEC(L0,NY,NX))
  ENDIF
  FXXCEC=FCO*XCEC(L0,NY,NX)
  XCEC(L1,NY,NX)=XCEC(L1,NY,NX)+FXXCEC
  XCEC(L0,NY,NX)=XCEC(L0,NY,NX)-FXXCEC
  FXXN4=FCO*XN4(L0,NY,NX)
  XN4(L1,NY,NX)=XN4(L1,NY,NX)+FXXN4
  XN4(L0,NY,NX)=XN4(L0,NY,NX)-FXXN4
  FXXNB=FCO*XNB(L0,NY,NX)
  XNB(L1,NY,NX)=XNB(L1,NY,NX)+FXXNB
  XNB(L0,NY,NX)=XNB(L0,NY,NX)-FXXNB
  FXXHY=FCO*XHY(L0,NY,NX)
  XHY(L1,NY,NX)=XHY(L1,NY,NX)+FXXHY
  XHY(L0,NY,NX)=XHY(L0,NY,NX)-FXXHY
  FXXAL=FCO*XAL(L0,NY,NX)
  XAL(L1,NY,NX)=XAL(L1,NY,NX)+FXXAL
  XAL(L0,NY,NX)=XAL(L0,NY,NX)-FXXAL
  FXXFE=FCO*XFE(L0,NY,NX)
  XFE(L1,NY,NX)=XFE(L1,NY,NX)+FXXFE
  XFE(L0,NY,NX)=XFE(L0,NY,NX)-FXXFE
  FXXCA=FCO*XCA(L0,NY,NX)
  XCA(L1,NY,NX)=XCA(L1,NY,NX)+FXXCA
  XCA(L0,NY,NX)=XCA(L0,NY,NX)-FXXCA
  FXXMG=FCO*XMG(L0,NY,NX)
  XMG(L1,NY,NX)=XMG(L1,NY,NX)+FXXMG
  XMG(L0,NY,NX)=XMG(L0,NY,NX)-FXXMG
  FXXNA=FCO*XNA(L0,NY,NX)
  XNA(L1,NY,NX)=XNA(L1,NY,NX)+FXXNA
  XNA(L0,NY,NX)=XNA(L0,NY,NX)-FXXNA
  FXXKA=FCO*XKA(L0,NY,NX)
  XKA(L1,NY,NX)=XKA(L1,NY,NX)+FXXKA
  XKA(L0,NY,NX)=XKA(L0,NY,NX)-FXXKA
  FXXHC=FCO*XHC(L0,NY,NX)
  XHC(L1,NY,NX)=XHC(L1,NY,NX)+FXXHC
  XHC(L0,NY,NX)=XHC(L0,NY,NX)-FXXHC
  FXXALO2=FCO*XALO2(L0,NY,NX)
  XALO2(L1,NY,NX)=XALO2(L1,NY,NX)+FXXALO2
  XALO2(L0,NY,NX)=XALO2(L0,NY,NX)-FXXALO2
  FXXFEO2=FCO*XFEO2(L0,NY,NX)
  XFEO2(L1,NY,NX)=XFEO2(L1,NY,NX)+FXXFEO2
  XFEO2(L0,NY,NX)=XFEO2(L0,NY,NX)-FXXFEO2
!
!     SOIL ADSORBED ANIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L.OR.AEC(L0,NY,NX).LE.ZERO)THEN
    FAO=FO
  ELSE
    FAO=AMIN1(0.5,FO*AEC(L1,NY,NX)/AEC(L0,NY,NX))
  ENDIF
  FXXAEC=FAO*XAEC(L0,NY,NX)
  XAEC(L1,NY,NX)=XAEC(L1,NY,NX)+FXXAEC
  XAEC(L0,NY,NX)=XAEC(L0,NY,NX)-FXXAEC
  FXXOH0=FAO*XOH0(L0,NY,NX)
  XOH0(L1,NY,NX)=XOH0(L1,NY,NX)+FXXOH0
  XOH0(L0,NY,NX)=XOH0(L0,NY,NX)-FXXOH0
  FXXOH1=FAO*XOH1(L0,NY,NX)
  XOH1(L1,NY,NX)=XOH1(L1,NY,NX)+FXXOH1
  XOH1(L0,NY,NX)=XOH1(L0,NY,NX)-FXXOH1
  FXXOH2=FAO*XOH2(L0,NY,NX)
  XOH2(L1,NY,NX)=XOH2(L1,NY,NX)+FXXOH2
  XOH2(L0,NY,NX)=XOH2(L0,NY,NX)-FXXOH2
  FXXH1P=FAO*XH1P(L0,NY,NX)
  XH1P(L1,NY,NX)=XH1P(L1,NY,NX)+FXXH1P
  XH1P(L0,NY,NX)=XH1P(L0,NY,NX)-FXXH1P
  FXXH2P=FAO*XH2P(L0,NY,NX)
  XH2P(L1,NY,NX)=XH2P(L1,NY,NX)+FXXH2P
  XH2P(L0,NY,NX)=XH2P(L0,NY,NX)-FXXH2P
  FXXOH0B=FAO*XOH0B(L0,NY,NX)
  XOH0B(L1,NY,NX)=XOH0B(L1,NY,NX)+FXXOH0B
  XOH0B(L0,NY,NX)=XOH0B(L0,NY,NX)-FXXOH0B
  FXXOH1B=FAO*XOH1B(L0,NY,NX)
  XOH1B(L1,NY,NX)=XOH1B(L1,NY,NX)+FXXOH1B
  XOH1B(L0,NY,NX)=XOH1B(L0,NY,NX)-FXXOH1B
  FXXOH2B=FAO*XOH2B(L0,NY,NX)
  XOH2B(L1,NY,NX)=XOH2B(L1,NY,NX)+FXXOH2B
  XOH2B(L0,NY,NX)=XOH2B(L0,NY,NX)-FXXOH2B
  FXXH1PB=FAO*XH1PB(L0,NY,NX)
  XH1PB(L1,NY,NX)=XH1PB(L1,NY,NX)+FXXH1PB
  XH1PB(L0,NY,NX)=XH1PB(L0,NY,NX)-FXXH1PB
  FXXH2PB=FAO*XH2PB(L0,NY,NX)
  XH2PB(L1,NY,NX)=XH2PB(L1,NY,NX)+FXXH2PB
  XH2PB(L0,NY,NX)=XH2PB(L0,NY,NX)-FXXH2PB
!
!     SOIL PRECIPITATES IN BAND, NON-BAND
!
  FXPALOH=AMIN1(FX*PALOH(L,NY,NX),PALOH(L0,NY,NX))
  PALOH(L1,NY,NX)=PALOH(L1,NY,NX)+FXPALOH
  PALOH(L0,NY,NX)=PALOH(L0,NY,NX)-FXPALOH
  FXPFEOH=AMIN1(FX*PFEOH(L,NY,NX),PFEOH(L0,NY,NX))
  PFEOH(L1,NY,NX)=PFEOH(L1,NY,NX)+FXPFEOH
  PFEOH(L0,NY,NX)=PFEOH(L0,NY,NX)-FXPFEOH
  FXPCACO=AMIN1(FX*PCACO(L,NY,NX),PCACO(L0,NY,NX))
  PCACO(L1,NY,NX)=PCACO(L1,NY,NX)+FXPCACO
  PCACO(L0,NY,NX)=PCACO(L0,NY,NX)-FXPCACO
  FXPCASO=AMIN1(FX*PCASO(L,NY,NX),PCASO(L0,NY,NX))
  PCASO(L1,NY,NX)=PCASO(L1,NY,NX)+FXPCASO
  PCASO(L0,NY,NX)=PCASO(L0,NY,NX)-FXPCASO
  FXPALPO=AMIN1(FX*PALPO(L,NY,NX),PALPO(L0,NY,NX))
  PALPO(L1,NY,NX)=PALPO(L1,NY,NX)+FXPALPO
  PALPO(L0,NY,NX)=PALPO(L0,NY,NX)-FXPALPO
  FXPFEPO=AMIN1(FX*PFEPO(L,NY,NX),PFEPO(L0,NY,NX))
  PFEPO(L1,NY,NX)=PFEPO(L1,NY,NX)+FXPFEPO
  PFEPO(L0,NY,NX)=PFEPO(L0,NY,NX)-FXPFEPO
  FXPCAPD=AMIN1(FX*PCAPD(L,NY,NX),PCAPD(L0,NY,NX))
  PCAPD(L1,NY,NX)=PCAPD(L1,NY,NX)+FXPCAPD
  PCAPD(L0,NY,NX)=PCAPD(L0,NY,NX)-FXPCAPD
  FXPCAPH=AMIN1(FX*PCAPH(L,NY,NX),PCAPH(L0,NY,NX))
  PCAPH(L1,NY,NX)=PCAPH(L1,NY,NX)+FXPCAPH
  PCAPH(L0,NY,NX)=PCAPH(L0,NY,NX)-FXPCAPH
  FXPCAPM=AMIN1(FX*PCAPM(L,NY,NX),PCAPM(L0,NY,NX))
  PCAPM(L1,NY,NX)=PCAPM(L1,NY,NX)+FXPCAPM
  PCAPM(L0,NY,NX)=PCAPM(L0,NY,NX)-FXPCAPM
  FXPALPB=AMIN1(FX*PALPB(L,NY,NX),PALPB(L0,NY,NX))
  PALPB(L1,NY,NX)=PALPB(L1,NY,NX)+FXPALPB
  PALPB(L0,NY,NX)=PALPB(L0,NY,NX)-FXPALPB
  FXPFEPB=AMIN1(FX*PFEPB(L,NY,NX),PFEPB(L0,NY,NX))
  PFEPB(L1,NY,NX)=PFEPB(L1,NY,NX)+FXPFEPB
  PFEPB(L0,NY,NX)=PFEPB(L0,NY,NX)-FXPFEPB
  FXPCPDB=AMIN1(FX*PCPDB(L,NY,NX),PCPDB(L0,NY,NX))
  PCPDB(L1,NY,NX)=PCPDB(L1,NY,NX)+FXPCPDB
  PCPDB(L0,NY,NX)=PCPDB(L0,NY,NX)-FXPCPDB
  FXPCPHB=AMIN1(FX*PCPHB(L,NY,NX),PCPHB(L0,NY,NX))
  PCPHB(L1,NY,NX)=PCPHB(L1,NY,NX)+FXPCPHB
  PCPHB(L0,NY,NX)=PCPHB(L0,NY,NX)-FXPCPHB
  FXPCPMB=AMIN1(FX*PCPMB(L,NY,NX),PCPMB(L0,NY,NX))
  PCPMB(L1,NY,NX)=PCPMB(L1,NY,NX)+FXPCPMB
  PCPMB(L0,NY,NX)=PCPMB(L0,NY,NX)-FXPCPMB

  end subroutine MoveBandSolute

!------------------------------------------------------------------------------------------


  Subroutine MoveDisolvGas(L0,L1,NY,NX,FX,FWO)
  implicit none
  integer, intent(in) :: L0,L1,NY,NX
  real(r8), intent(in) :: FX,FWO

  real(r8) :: FXH2GS,FXCO2G,FXCH4G,FXOXYG,FXZ2GG,FXZ2OG
  real(r8) :: FXZNH3G,FXH2GG,FXCO2S,FXCH4S,FXOXYS,FXZ2GS,FXZ2OS

!
!     SOIL GASEOUS GASES
!
  FXCO2G=FWO*CO2G(L0,NY,NX)
  CO2G(L1,NY,NX)=CO2G(L1,NY,NX)+FXCO2G
  CO2G(L0,NY,NX)=CO2G(L0,NY,NX)-FXCO2G
  FXCH4G=FWO*CH4G(L0,NY,NX)
  CH4G(L1,NY,NX)=CH4G(L1,NY,NX)+FXCH4G
  CH4G(L0,NY,NX)=CH4G(L0,NY,NX)-FXCH4G
  FXOXYG=FWO*OXYG(L0,NY,NX)
  OXYG(L1,NY,NX)=OXYG(L1,NY,NX)+FXOXYG
  OXYG(L0,NY,NX)=OXYG(L0,NY,NX)-FXOXYG
  FXZ2GG=FWO*Z2GG(L0,NY,NX)
  Z2GG(L1,NY,NX)=Z2GG(L1,NY,NX)+FXZ2GG
  Z2GG(L0,NY,NX)=Z2GG(L0,NY,NX)-FXZ2GG
  FXZ2OG=FWO*Z2OG(L0,NY,NX)
  Z2OG(L1,NY,NX)=Z2OG(L1,NY,NX)+FXZ2OG
  Z2OG(L0,NY,NX)=Z2OG(L0,NY,NX)-FXZ2OG
  FXZNH3G=FWO*ZNH3G(L0,NY,NX)
  ZNH3G(L1,NY,NX)=ZNH3G(L1,NY,NX)+FXZNH3G
  ZNH3G(L0,NY,NX)=ZNH3G(L0,NY,NX)-FXZNH3G
  FXH2GG=FWO*H2GG(L0,NY,NX)
  H2GG(L1,NY,NX)=H2GG(L1,NY,NX)+FXH2GG
  H2GG(L0,NY,NX)=H2GG(L0,NY,NX)-FXH2GG

  FXCO2S=FWO*CO2S(L0,NY,NX)
  CO2S(L1,NY,NX)=CO2S(L1,NY,NX)+FXCO2S
  CO2S(L0,NY,NX)=CO2S(L0,NY,NX)-FXCO2S
  FXCH4S=FWO*CH4S(L0,NY,NX)
  CH4S(L1,NY,NX)=CH4S(L1,NY,NX)+FXCH4S
  CH4S(L0,NY,NX)=CH4S(L0,NY,NX)-FXCH4S
  FXOXYS=FWO*OXYS(L0,NY,NX)
  OXYS(L1,NY,NX)=OXYS(L1,NY,NX)+FXOXYS
  OXYS(L0,NY,NX)=OXYS(L0,NY,NX)-FXOXYS
  FXZ2GS=FWO*Z2GS(L0,NY,NX)
  Z2GS(L1,NY,NX)=Z2GS(L1,NY,NX)+FXZ2GS
  Z2GS(L0,NY,NX)=Z2GS(L0,NY,NX)-FXZ2GS
  FXZ2OS=FWO*Z2OS(L0,NY,NX)
  Z2OS(L1,NY,NX)=Z2OS(L1,NY,NX)+FXZ2OS
  Z2OS(L0,NY,NX)=Z2OS(L0,NY,NX)-FXZ2OS
  FXH2GS=FWO*H2GS(L0,NY,NX)
  H2GS(L1,NY,NX)=H2GS(L1,NY,NX)+FXH2GS
  H2GS(L0,NY,NX)=H2GS(L0,NY,NX)-FXH2GS
  end Subroutine MoveDisolvGas

!------------------------------------------------------------------------------------------

  subroutine MoveFertSalt(L,L0,L1,NY,NX,FX,FWO)

  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in):: FX,FWO
  real(r8) :: FXZCA0P,FXZMG1P,FXZCA2P,FXZCA1P
  real(r8) :: FXZSO4,FXZOH,FXZCO3
  real(r8) :: FXZNHUFA,FXZNO3FA,FXZNH4FB,FXZNH3FB,FXZNHUFB,FXZNO3FB
  real(r8) :: FXZNH4S,FXZNH4B,FXZNH3S,FXZNH3B,FXZNO3S,FXZNO3B
  real(r8) :: FXZNO2S,FXZNO2B,FXH1PO4,FXH2PO4,FXZAL,FXZFE,FXZHY
  real(r8) :: FXZNH4FA,FXZNH3FA
  real(r8) :: FXZMGO,FXZMGC,FXZMGH,FXZMGS
  real(r8) :: FXZNAC,FXZNAS,FXZKAS,FXH0PO4,FXH3PO4,FXZFE1P,FXZFE2P
  real(r8) :: FXZHCO3,FXZALOH1,FXZALOH2,FXZALOH3,FXZALOH4,FXZALS
  real(r8) :: FXZMG,FXZNA,FXZKA,FXZCL
  real(r8) :: FXZFEOH1,FXZFEOH2,FXZFEOH3,FXZFEOH4,FXZFES,FXZCAO
  real(r8) :: FXZCAC,FXZCAH,FXZCAS,FXZCA

! begin_execution
  FXZNH4FA=AMIN1(FX*ZNH4FA(L,NY,NX),ZNH4FA(L0,NY,NX))
  ZNH4FA(L1,NY,NX)=ZNH4FA(L1,NY,NX)+FXZNH4FA
  ZNH4FA(L0,NY,NX)=ZNH4FA(L0,NY,NX)-FXZNH4FA
  FXZNH3FA=AMIN1(FX*ZNH3FA(L,NY,NX),ZNH3FA(L0,NY,NX))
  ZNH3FA(L1,NY,NX)=ZNH3FA(L1,NY,NX)+FXZNH3FA
  ZNH3FA(L0,NY,NX)=ZNH3FA(L0,NY,NX)-FXZNH3FA
  FXZNHUFA=AMIN1(FX*ZNHUFA(L,NY,NX),ZNHUFA(L0,NY,NX))
  ZNHUFA(L1,NY,NX)=ZNHUFA(L1,NY,NX)+FXZNHUFA
  ZNHUFA(L0,NY,NX)=ZNHUFA(L0,NY,NX)-FXZNHUFA
  FXZNO3FA=AMIN1(FX*ZNO3FA(L,NY,NX),ZNO3FA(L0,NY,NX))
  ZNO3FA(L1,NY,NX)=ZNO3FA(L1,NY,NX)+FXZNO3FA
  ZNO3FA(L0,NY,NX)=ZNO3FA(L0,NY,NX)-FXZNO3FA
  FXZNH4FB=AMIN1(FX*ZNH4FB(L,NY,NX),ZNH4FB(L0,NY,NX))
  ZNH4FB(L1,NY,NX)=ZNH4FB(L1,NY,NX)+FXZNH4FB
  ZNH4FB(L0,NY,NX)=ZNH4FB(L0,NY,NX)-FXZNH4FB
  FXZNH3FB=AMIN1(FX*ZNH3FB(L,NY,NX),ZNH3FB(L0,NY,NX))
  ZNH3FB(L1,NY,NX)=ZNH3FB(L1,NY,NX)+FXZNH3FB
  ZNH3FB(L0,NY,NX)=ZNH3FB(L0,NY,NX)-FXZNH3FB
  FXZNHUFB=AMIN1(FX*ZNHUFB(L,NY,NX),ZNHUFB(L0,NY,NX))
  ZNHUFB(L1,NY,NX)=ZNHUFB(L1,NY,NX)+FXZNHUFB
  ZNHUFB(L0,NY,NX)=ZNHUFB(L0,NY,NX)-FXZNHUFB
  FXZNO3FB=AMIN1(FX*ZNO3FB(L,NY,NX),ZNO3FB(L0,NY,NX))
  ZNO3FB(L1,NY,NX)=ZNO3FB(L1,NY,NX)+FXZNO3FB
  ZNO3FB(L0,NY,NX)=ZNO3FB(L0,NY,NX)-FXZNO3FB
!
!     SOIL N,P SOLUTES IN BAND, NON-BAND
!
  FXZNH4S=FWO*ZNH4S(L0,NY,NX)
  ZNH4S(L1,NY,NX)=ZNH4S(L1,NY,NX)+FXZNH4S
  ZNH4S(L0,NY,NX)=ZNH4S(L0,NY,NX)-FXZNH4S
  FXZNH4B=FWO*ZNH4B(L0,NY,NX)
  ZNH4B(L1,NY,NX)=ZNH4B(L1,NY,NX)+FXZNH4B
  ZNH4B(L0,NY,NX)=ZNH4B(L0,NY,NX)-FXZNH4B
  FXZNH3S=FWO*ZNH3S(L0,NY,NX)
  ZNH3S(L1,NY,NX)=ZNH3S(L1,NY,NX)+FXZNH3S
  ZNH3S(L0,NY,NX)=ZNH3S(L0,NY,NX)-FXZNH3S
  FXZNH3B=FWO*ZNH3B(L0,NY,NX)
  ZNH3B(L1,NY,NX)=ZNH3B(L1,NY,NX)+FXZNH3B
  ZNH3B(L0,NY,NX)=ZNH3B(L0,NY,NX)-FXZNH3B
  FXZNO3S=FWO*ZNO3S(L0,NY,NX)
  ZNO3S(L1,NY,NX)=ZNO3S(L1,NY,NX)+FXZNO3S
  ZNO3S(L0,NY,NX)=ZNO3S(L0,NY,NX)-FXZNO3S
  FXZNO3B=FWO*ZNO3B(L0,NY,NX)
  ZNO3B(L1,NY,NX)=ZNO3B(L1,NY,NX)+FXZNO3B
  ZNO3B(L0,NY,NX)=ZNO3B(L0,NY,NX)-FXZNO3B
  FXZNO2S=FWO*ZNO2S(L0,NY,NX)
  ZNO2S(L1,NY,NX)=ZNO2S(L1,NY,NX)+FXZNO2S
  ZNO2S(L0,NY,NX)=ZNO2S(L0,NY,NX)-FXZNO2S
  FXZNO2B=FWO*ZNO2B(L0,NY,NX)
  ZNO2B(L1,NY,NX)=ZNO2B(L1,NY,NX)+FXZNO2B
  ZNO2B(L0,NY,NX)=ZNO2B(L0,NY,NX)-FXZNO2B
  FXH1PO4=FWO*H1PO4(L0,NY,NX)
  H1PO4(L1,NY,NX)=H1PO4(L1,NY,NX)+FXH1PO4
  H1PO4(L0,NY,NX)=H1PO4(L0,NY,NX)-FXH1PO4
  FXH2PO4=FWO*H2PO4(L0,NY,NX)
  H2PO4(L1,NY,NX)=H2PO4(L1,NY,NX)+FXH2PO4
  H2PO4(L0,NY,NX)=H2PO4(L0,NY,NX)-FXH2PO4
!
!     SOIL SALT SOLUTES
!
  IF(ISALTG.NE.0)THEN
    FXZAL=FWO*ZAL(L0,NY,NX)
    ZAL(L1,NY,NX)=ZAL(L1,NY,NX)+FXZAL
    ZAL(L0,NY,NX)=ZAL(L0,NY,NX)-FXZAL
    FXZFE=FWO*ZFE(L0,NY,NX)
    ZFE(L1,NY,NX)=ZFE(L1,NY,NX)+FXZFE
    ZFE(L0,NY,NX)=ZFE(L0,NY,NX)-FXZFE
    FXZHY=FWO*ZHY(L0,NY,NX)
    ZHY(L1,NY,NX)=ZHY(L1,NY,NX)+FXZHY
    ZHY(L0,NY,NX)=ZHY(L0,NY,NX)-FXZHY
    FXZCA=FWO*ZCA(L0,NY,NX)
    ZCA(L1,NY,NX)=ZCA(L1,NY,NX)+FXZCA
    ZCA(L0,NY,NX)=ZCA(L0,NY,NX)-FXZCA
    FXZMG=FWO*ZMG(L0,NY,NX)
    ZMG(L1,NY,NX)=ZMG(L1,NY,NX)+FXZMG
    ZMG(L0,NY,NX)=ZMG(L0,NY,NX)-FXZMG
    FXZNA=FWO*ZNA(L0,NY,NX)
    ZNA(L1,NY,NX)=ZNA(L1,NY,NX)+FXZNA
    ZNA(L0,NY,NX)=ZNA(L0,NY,NX)-FXZNA
    FXZKA=FWO*ZKA(L0,NY,NX)
    ZKA(L1,NY,NX)=ZKA(L1,NY,NX)+FXZKA
    ZKA(L0,NY,NX)=ZKA(L0,NY,NX)-FXZKA
    FXZOH=FWO*ZOH(L0,NY,NX)
    ZOH(L1,NY,NX)=ZOH(L1,NY,NX)+FXZOH
    ZOH(L0,NY,NX)=ZOH(L0,NY,NX)-FXZOH
    FXZSO4=FWO*ZSO4(L0,NY,NX)
    ZSO4(L1,NY,NX)=ZSO4(L1,NY,NX)+FXZSO4
    ZSO4(L0,NY,NX)=ZSO4(L0,NY,NX)-FXZSO4
    FXZCL=FWO*ZCL(L0,NY,NX)
    ZCL(L1,NY,NX)=ZCL(L1,NY,NX)+FXZCL
    ZCL(L0,NY,NX)=ZCL(L0,NY,NX)-FXZCL
    FXZCO3=FWO*ZCO3(L0,NY,NX)
    ZCO3(L1,NY,NX)=ZCO3(L1,NY,NX)+FXZCO3
    ZCO3(L0,NY,NX)=ZCO3(L0,NY,NX)-FXZCO3
    FXZHCO3=FWO*ZHCO3(L0,NY,NX)
    ZHCO3(L1,NY,NX)=ZHCO3(L1,NY,NX)+FXZHCO3
    ZHCO3(L0,NY,NX)=ZHCO3(L0,NY,NX)-FXZHCO3
    FXZALOH1=FWO*ZALOH1(L0,NY,NX)
    ZALOH1(L1,NY,NX)=ZALOH1(L1,NY,NX)+FXZALOH1
    ZALOH1(L0,NY,NX)=ZALOH1(L0,NY,NX)-FXZALOH1
    FXZALOH2=FWO*ZALOH2(L0,NY,NX)
    ZALOH2(L1,NY,NX)=ZALOH2(L1,NY,NX)+FXZALOH2
    ZALOH2(L0,NY,NX)=ZALOH2(L0,NY,NX)-FXZALOH2
    FXZALOH3=FWO*ZALOH3(L0,NY,NX)
    ZALOH3(L1,NY,NX)=ZALOH3(L1,NY,NX)+FXZALOH3
    ZALOH3(L0,NY,NX)=ZALOH3(L0,NY,NX)-FXZALOH3
    FXZALOH4=FWO*ZALOH4(L0,NY,NX)
    ZALOH4(L1,NY,NX)=ZALOH4(L1,NY,NX)+FXZALOH4
    ZALOH4(L0,NY,NX)=ZALOH4(L0,NY,NX)-FXZALOH4
    FXZALS=FWO*ZALS(L0,NY,NX)
    ZALS(L1,NY,NX)=ZALS(L1,NY,NX)+FXZALS
    ZALS(L0,NY,NX)=ZALS(L0,NY,NX)-FXZALS
    FXZFEOH1=FWO*ZFEOH1(L0,NY,NX)
    ZFEOH1(L1,NY,NX)=ZFEOH1(L1,NY,NX)+FXZFEOH1
    ZFEOH1(L0,NY,NX)=ZFEOH1(L0,NY,NX)-FXZFEOH1
    FXZFEOH2=FWO*ZFEOH2(L0,NY,NX)
    ZFEOH2(L1,NY,NX)=ZFEOH2(L1,NY,NX)+FXZFEOH2
    ZFEOH2(L0,NY,NX)=ZFEOH2(L0,NY,NX)-FXZFEOH2
    FXZFEOH3=FWO*ZFEOH3(L0,NY,NX)
    ZFEOH3(L1,NY,NX)=ZFEOH3(L1,NY,NX)+FXZFEOH3
    ZFEOH3(L0,NY,NX)=ZFEOH3(L0,NY,NX)-FXZFEOH3
    FXZFEOH4=FWO*ZFEOH4(L0,NY,NX)
    ZFEOH4(L1,NY,NX)=ZFEOH4(L1,NY,NX)+FXZFEOH4
    ZFEOH4(L0,NY,NX)=ZFEOH4(L0,NY,NX)-FXZFEOH4
    FXZFES=FWO*ZFES(L0,NY,NX)
    ZFES(L1,NY,NX)=ZFES(L1,NY,NX)+FXZFES
    ZFES(L0,NY,NX)=ZFES(L0,NY,NX)-FXZFES
    FXZCAO=FWO*ZCAO(L0,NY,NX)
    ZCAO(L1,NY,NX)=ZCAO(L1,NY,NX)+FXZCAO
    ZCAO(L0,NY,NX)=ZCAO(L0,NY,NX)-FXZCAO
    FXZCAC=FWO*ZCAC(L0,NY,NX)
    ZCAC(L1,NY,NX)=ZCAC(L1,NY,NX)+FXZCAC
    ZCAC(L0,NY,NX)=ZCAC(L0,NY,NX)-FXZCAC
    FXZCAH=FWO*ZCAH(L0,NY,NX)
    ZCAH(L1,NY,NX)=ZCAH(L1,NY,NX)+FXZCAH
    ZCAH(L0,NY,NX)=ZCAH(L0,NY,NX)-FXZCAH
    FXZCAS=FWO*ZCAS(L0,NY,NX)
    ZCAS(L1,NY,NX)=ZCAS(L1,NY,NX)+FXZCAS
    ZCAS(L0,NY,NX)=ZCAS(L0,NY,NX)-FXZCAS
    FXZMGO=FWO*ZMGO(L0,NY,NX)
    ZMGO(L1,NY,NX)=ZMGO(L1,NY,NX)+FXZMGO
    ZMGO(L0,NY,NX)=ZMGO(L0,NY,NX)-FXZMGO
    FXZMGC=FWO*ZMGC(L0,NY,NX)
    ZMGC(L1,NY,NX)=ZMGC(L1,NY,NX)+FXZMGC
    ZMGC(L0,NY,NX)=ZMGC(L0,NY,NX)-FXZMGC
    FXZMGH=FWO*ZMGH(L0,NY,NX)
    ZMGH(L1,NY,NX)=ZMGH(L1,NY,NX)+FXZMGH
    ZMGH(L0,NY,NX)=ZMGH(L0,NY,NX)-FXZMGH
    FXZMGS=FWO*ZMGS(L0,NY,NX)
    ZMGS(L1,NY,NX)=ZMGS(L1,NY,NX)+FXZMGS
    ZMGS(L0,NY,NX)=ZMGS(L0,NY,NX)-FXZMGS
    FXZNAC=FWO*ZNAC(L0,NY,NX)
    ZNAC(L1,NY,NX)=ZNAC(L1,NY,NX)+FXZNAC
    ZNAC(L0,NY,NX)=ZNAC(L0,NY,NX)-FXZNAC
    FXZNAS=FWO*ZNAS(L0,NY,NX)
    ZNAS(L1,NY,NX)=ZNAS(L1,NY,NX)+FXZNAS
    ZNAS(L0,NY,NX)=ZNAS(L0,NY,NX)-FXZNAS
    FXZKAS=FWO*ZKAS(L0,NY,NX)
    ZKAS(L1,NY,NX)=ZKAS(L1,NY,NX)+FXZKAS
    ZKAS(L0,NY,NX)=ZKAS(L0,NY,NX)-FXZKAS
    FXH0PO4=FWO*H0PO4(L0,NY,NX)
    H0PO4(L1,NY,NX)=H0PO4(L1,NY,NX)+FXH0PO4
    H0PO4(L0,NY,NX)=H0PO4(L0,NY,NX)-FXH0PO4
    FXH3PO4=FWO*H3PO4(L0,NY,NX)
    H3PO4(L1,NY,NX)=H3PO4(L1,NY,NX)+FXH3PO4
    H3PO4(L0,NY,NX)=H3PO4(L0,NY,NX)-FXH3PO4
    FXZFE1P=FWO*ZFE1P(L0,NY,NX)
    ZFE1P(L1,NY,NX)=ZFE1P(L1,NY,NX)+FXZFE1P
    ZFE1P(L0,NY,NX)=ZFE1P(L0,NY,NX)-FXZFE1P
    FXZFE2P=FWO*ZFE2P(L0,NY,NX)
    ZFE2P(L1,NY,NX)=ZFE2P(L1,NY,NX)+FXZFE2P
    ZFE2P(L0,NY,NX)=ZFE2P(L0,NY,NX)-FXZFE2P
    FXZCA0P=FWO*ZCA0P(L0,NY,NX)
    ZCA0P(L1,NY,NX)=ZCA0P(L1,NY,NX)+FXZCA0P
    ZCA0P(L0,NY,NX)=ZCA0P(L0,NY,NX)-FXZCA0P
    FXZCA1P=FWO*ZCA1P(L0,NY,NX)
    ZCA1P(L1,NY,NX)=ZCA1P(L1,NY,NX)+FXZCA1P
    ZCA1P(L0,NY,NX)=ZCA1P(L0,NY,NX)-FXZCA1P
    FXZCA2P=FWO*ZCA2P(L0,NY,NX)
    ZCA2P(L1,NY,NX)=ZCA2P(L1,NY,NX)+FXZCA2P
    ZCA2P(L0,NY,NX)=ZCA2P(L0,NY,NX)-FXZCA2P
    FXZMG1P=FWO*ZMG1P(L0,NY,NX)
    ZMG1P(L1,NY,NX)=ZMG1P(L1,NY,NX)+FXZMG1P
    ZMG1P(L0,NY,NX)=ZMG1P(L0,NY,NX)-FXZMG1P
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
      IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPNH4(NY,NX))THEN
        WDNHBDL=WDNHB(L,NY,NX)*DLYR(3,L,NY,NX)
        WDNHBD0=WDNHB(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDNHBD1=WDNHB(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDNHB=AMIN1(FX*WDNHBDL,WDNHBD0)
        WDNHBD1=WDNHBD1+FXWDNHB
        WDNHBD0=WDNHBD0-FXWDNHB
        WDNHB(L1,NY,NX)=WDNHBD1/DLYR(3,L1,NY,NX)
        WDNHB(L0,NY,NX)=WDNHBD0/DLYR(3,L0,NY,NX)
        IF(CDPTH(L,NY,NX).GE.DPNH4(NY,NX))THEN
          FXDPNHB=AMIN1(FX*DPNHB(L,NY,NX),DPNHB(L0,NY,NX))
          DPNHB(L1,NY,NX)=DPNHB(L1,NY,NX)+FXDPNHB
          DPNHB(L0,NY,NX)=DPNHB(L0,NY,NX)-FXDPNHB
        ENDIF
        VLNHB(L1,NY,NX)=AZMAX1(AMIN1(0.999,WDNHB(L1,NY,NX) &
          /ROWN(NY,NX)*DPNHB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        VLNHB(L0,NY,NX)=AZMAX1(AMIN1(0.999,WDNHB(L0,NY,NX) &
          /ROWN(NY,NX)*DPNHB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        VLNH4(L1,NY,NX)=1.0-VLNHB(L1,NY,NX)
        VLNH4(L0,NY,NX)=1.0-VLNHB(L0,NY,NX)
      ENDIF
    ENDIF
    IF(IFNOB(NY,NX).EQ.1.AND.ROWO(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPNO3(NY,NX))THEN
        WDNOBDL=WDNOB(L,NY,NX)*DLYR(3,L,NY,NX)
        WDNOBD0=WDNOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDNOBD1=WDNOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDNOB=AMIN1(FX*WDNOBDL,WDNOBD0)
        WDNOBD1=WDNOBD1+FXWDNOB
        WDNOBD0=WDNOBD0-FXWDNOB
        WDNOB(L1,NY,NX)=WDNOBD1/DLYR(3,L1,NY,NX)
        WDNOB(L0,NY,NX)=WDNOBD0/DLYR(3,L0,NY,NX)
        IF(CDPTH(L,NY,NX).GE.DPNO3(NY,NX))THEN
          FXDPNOB=AMIN1(FX*DPNOB(L,NY,NX),DPNOB(L0,NY,NX))
          DPNOB(L1,NY,NX)=DPNOB(L1,NY,NX)+FXDPNOB
          DPNOB(L0,NY,NX)=DPNOB(L0,NY,NX)-FXDPNOB
        ENDIF
        VLNOB(L1,NY,NX)=AZMAX1(AMIN1(0.999,WDNOB(L1,NY,NX) &
          /ROWO(NY,NX)*DPNOB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        VLNOB(L0,NY,NX)=AZMAX1(AMIN1(0.999,WDNOB(L0,NY,NX) &
          /ROWO(NY,NX)*DPNOB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        VLNO3(L1,NY,NX)=1.0-VLNOB(L1,NY,NX)
          VLNO3(L0,NY,NX)=1.0-VLNOB(L0,NY,NX)
      ENDIF
    ENDIF
    IF(IFPOB(NY,NX).EQ.1.AND.ROWP(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPPO4(NY,NX))THEN
        WDPOBDL=WDPOB(L,NY,NX)*DLYR(3,L,NY,NX)
        WDPOBD0=WDPOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
        WDPOBD1=WDPOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
        FXWDPOB=AMIN1(FX*WDPOBDL,WDPOBD0)
        WDPOBD1=WDPOBD1+FXWDPOB
        WDPOBD0=WDPOBD0-FXWDPOB
        WDPOB(L1,NY,NX)=WDPOBD1/DLYR(3,L1,NY,NX)
        WDPOB(L0,NY,NX)=WDPOBD0/DLYR(3,L0,NY,NX)
        IF(CDPTH(L,NY,NX).GE.DPPO4(NY,NX))THEN
          FXDPPOB=AMIN1(FX*DPPOB(L,NY,NX),DPPOB(L0,NY,NX))
          DPPOB(L1,NY,NX)=DPPOB(L1,NY,NX)+FXDPPOB
          DPPOB(L0,NY,NX)=DPPOB(L0,NY,NX)-FXDPPOB
        ENDIF
        VLPOB(L1,NY,NX)=AZMAX1(AMIN1(0.999,WDPOB(L1,NY,NX) &
          /ROWP(NY,NX)*DPPOB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
        VLPOB(L0,NY,NX)=AZMAX1(AMIN1(0.999,WDPOB(L0,NY,NX) &
          /ROWP(NY,NX)*DPPOB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
        VLPO4(L1,NY,NX)=1.0-VLPOB(L1,NY,NX)
        VLPO4(L0,NY,NX)=1.0-VLPOB(L0,NY,NX)
      ENDIF
    ENDIF
  ENDIF
!
!     SOIL MINERALS
!
  IF(L0.EQ.L.OR.BKDSI(L0,NY,NX).LE.ZERO)THEN
    FBO=FX
  ELSE
    FBO=AMIN1(0.1,FX*BKDSI(L1,NY,NX)/BKDSI(L0,NY,NX))
  ENDIF
!     BKDS(L1,NY,NX)=(1.0-FO)*BKDS(L1,NY,NX)+FO*BKDSI(L0,NY,NX)
  PH(L1,NY,NX)=(1.0-FO)*PH(L1,NY,NX)+FO*PH(L0,NY,NX)
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
  IF(FHOL(L1,NY,NX).GT.ZERO.AND.FHOL(L0,NY,NX).GT.ZERO)THEN
    IF(L0.EQ.L.OR.FHOLI(L0,NY,NX).LE.ZERO)THEN
      FHO=FO
    ELSE
      FHO=AMIN1(0.5,FO*FHOLI(L1,NY,NX)/FHOLI(L0,NY,NX))
    ENDIF
    FHOL(L1,NY,NX)=(1.0-FO)*FHOL(L1,NY,NX)+FO*FHOL(L0,NY,NX)
    FXVOLWH=FHO*VOLWH(L0,NY,NX)
    VOLWH(L1,NY,NX)=VOLWH(L1,NY,NX)+FXVOLWH
    VOLWH(L0,NY,NX)=VOLWH(L0,NY,NX)-FXVOLWH
    FXVOLIH=FHO*VOLIH(L0,NY,NX)
    VOLIH(L1,NY,NX)=VOLIH(L1,NY,NX)+FXVOLIH
    VOLIH(L0,NY,NX)=VOLIH(L0,NY,NX)-FXVOLIH
    FXVOLAH=FHO*VOLAH(L0,NY,NX)
    VOLAH(L1,NY,NX)=VOLAH(L1,NY,NX)+FXVOLAH
    VOLAH(L0,NY,NX)=VOLAH(L0,NY,NX)-FXVOLAH
  ENDIF
  end subroutine MoveFertMinerals

!------------------------------------------------------------------------------------------

  subroutine MoveHeatWat(L,L0,L1,NY,NX,FO,FX)
  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in) :: FO,FX
  real(r8) :: FWO,FXVOLW
  real(r8) :: FXENGY,FXVOLI,FXVOLY,FXVOLWX,FXVHCM
  real(r8) :: ENGY0,ENGY1

  IF(L0.EQ.L.OR.(L0>0 .and. POROSI(L0,NY,NX).LE.ZERO))THEN
    FWO=FO
  ELSE
    FWO=AMIN1(0.5,FO*POROSI(L1,NY,NX)/POROSI(L0,NY,NX))
  ENDIF
!             FXSCNV=FWO*SCNV(L0,NY,NX)
!             SCNV(L1,NY,NX)=SCNV(L1,NY,NX)+FXSCNV
!             SCNV(L0,NY,NX)=SCNV(L0,NY,NX)-FXSCNV
!             FXSCNH=FWO*SCNH(L0,NY,NX)
!             SCNH(L1,NY,NX)=SCNH(L1,NY,NX)+FXSCNH
!             SCNH(L0,NY,NX)=SCNH(L0,NY,NX)-FXSCNH
  IF(L0.EQ.0)THEN
    FXVOLW=FX*AZMAX1(XVOLWP-VOLWD(NY,NX))
  ELSE
    FXVOLW=FWO*VOLW(L0,NY,NX)
  ENDIF
  VOLW(L1,NY,NX)=VOLW(L1,NY,NX)+FXVOLW
  VOLW(L0,NY,NX)=VOLW(L0,NY,NX)-FXVOLW
!     IF(VOLI(L1,NY,NX).GT.ZEROS(NY,NX))THEN
  FXVOLI=FWO*VOLI(L0,NY,NX)
  VOLI(L1,NY,NX)=VOLI(L1,NY,NX)+FXVOLI
  VOLI(L0,NY,NX)=VOLI(L0,NY,NX)-FXVOLI
!     ENDIF
!     FXVOLA=FWO*VOLA(L0,NY,NX)
!     IF(L1.NE.NU(NY,NX))THEN
!     VOLA(L1,NY,NX)=VOLA(L1,NY,NX)+FXVOLA
!     ENDIF
!     IF(L0.NE.NU(NY,NX))THEN
!     VOLA(L0,NY,NX)=VOLA(L0,NY,NX)-FXVOLA
!     ENDIF
  FXVOLY=FWO*VOLY(L0,NY,NX)
  VOLY(L1,NY,NX)=VOLY(L1,NY,NX)+FXVOLY
  VOLY(L0,NY,NX)=VOLY(L0,NY,NX)-FXVOLY
  FXVOLWX=FWO*VOLWX(L0,NY,NX)
  VOLWX(L1,NY,NX)=VOLWX(L1,NY,NX)+FXVOLWX
  VOLWX(L0,NY,NX)=VOLWX(L0,NY,NX)-FXVOLWX
  FXVHCM=FWO*VHCM(L0,NY,NX)
  VHCM(L1,NY,NX)=VHCM(L1,NY,NX)+FXVHCM
  VHCM(L0,NY,NX)=VHCM(L0,NY,NX)-FXVHCM
  FXENGY=TKS(L0,NY,NX)*(FXVHCM+cpw*FXVOLW+cpi*FXVOLI)
  ENGY1=VHCP(L1,NY,NX)*TKS(L1,NY,NX)+FXENGY
  ENGY0=VHCP(L0,NY,NX)*TKS(L0,NY,NX)-FXENGY
  VHCP(L1,NY,NX)=VHCP(L1,NY,NX)+FXVHCM+cpw*FXVOLW+cpi*FXVOLI
  VHCP(L0,NY,NX)=VHCP(L0,NY,NX)-FXVHCM-cpw*FXVOLW-cpi*FXVOLI
  IF(VHCP(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L1,NY,NX)=ENGY1/VHCP(L1,NY,NX)
  ELSE
    TKS(L1,NY,NX)=TKS(L,NY,NX)
  ENDIF
  TCS(L1,NY,NX)=TKS(L1,NY,NX)-TC2K
  IF(VHCP(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS(L0,NY,NX)=ENGY0/VHCP(L0,NY,NX)
  ELSE
    TKS(L0,NY,NX)=TKS(L,NY,NX)
  ENDIF
  TCS(L0,NY,NX)=TKS(L0,NY,NX)-TC2K
  end subroutine MoveHeatWat

end module SoilLayerDynMod
