module HydrologyDiagMod  
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: isclose, AZMAX1, AZMIN1
  use EcosimConst,   only: mGravAccelerat
  use SoilWaterDataType
  use GridDataType  
  use EcoSIMCtrlDataType
  use FlagDataType
  use SoilPropertyDataType
  use SoilPhysDataType
  implicit none
  private

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: DiagWaterTBLDepz

  contains
!------------------------------------------------------------------------------------------

  subroutine DiagWaterTBLDepz(I,J,NY,NX)
  !
  !Description:
  !Diagnose water table depth in the soil column
  implicit none
  integer, intent(in) :: I,J,NY,NX

  real(r8) :: THETPZ_vr(JZ)  !air-filled porosity in layer
  real(r8) :: PSIEquil       !equilibrium matric potential
  real(r8) :: THETW1
  real(r8) :: THETWM
  real(r8) :: THETPX
  real(r8) :: THETPW !     THETPW=minimum air-filled porosity for saturation (m3 m-3)  
  real(r8) :: THETWP  
  integer :: LL,L
  logical :: FoundWaterTable

!     begin_execution

  THETPW          = 0.01_r8
  THETWP          = 1.0_r8-THETPW
  FoundWaterTable = .false.

  DO L=NUI_col(NY,NX),NLI_col(NY,NX)
!     IDWaterTable=water table flag from site file
!     THETPZ,THETPW=current,minimum air-filled, porosity for water table
!     DPTH,ExtWaterTable=depth of soil layer midpoint, water table
!     PSIEquil=water potential in hydraulic equilibrium with layer below
!     THETW1,THETWP=water content at PSIEquil,minimum SWC for water table
!     DepzIntWTBL_col=water table depth
!
    IF(VLSoilPoreMicP_vr(L,NY,NX).LE.ZEROS(NY,NX))THEN
      THETW_vr(L,NY,NX)    = POROS_vr(L,NY,NX)
      THETI_vr(L,NY,NX)    = 0._r8
      ThetaAir_vr(L,NY,NX) = 0._r8
    ELSE
      THETW_vr(L,NY,NX)    = AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLWatMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)))
      THETI_vr(L,NY,NX)    = AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLiceMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)))
      ThetaAir_vr(L,NY,NX) = AZMAX1(VLsoiAirP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX))
    ENDIF

    THETPZ_vr(L)=AZMAX1(POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX)-THETI_vr(L,NY,NX))

    IF(IDWaterTable_col(NY,NX).NE.0)THEN
      IF(FoundWaterTable)exit
      
      IF(THETPZ_vr(L).LT.THETPW &         !layer is saturated
        .OR. L.EQ.NL_col(NY,NX))THEN      !bottom layer
        FoundWaterTable=.true.
        IF(SoilDepthMidLay_vr(L,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN   !current middle depth is above external water table
          D5705: DO LL=MIN(L+1,NL_col(NY,NX)),NL_col(NY,NX)
            IF(THETPZ_vr(LL).GE.THETPW .AND. LL.NE.NL_col(NY,NX))THEN !it is an unsaturated inner layer
              !air-filled pore greater minimum, i.e. not saturated
              FoundWaterTable=.false.
              exit
            ELSE IF(SoilDepthMidLay_vr(LL,NY,NX).GE.ExtWaterTable_col(NY,NX))THEN  !
              !current layer is lower than external water table
              exit
            ENDIF
          END DO D5705
        ENDIF

        !THETPW=saturation criterion for water table identification
        IF(FoundWaterTable)THEN
          IF(THETPZ_vr(L).GE.THETPW .AND. L.NE.NL_col(NY,NX))THEN !saturated and inside the hydrologically active zone
            !not bottom layer, saturated
            !PSIeqv in saturated layer
            PSIEquil = PSISoilMatricP_vr(L+1,NY,NX)-mGravAccelerat*(SoilDepthMidLay_vr(L+1,NY,NX)-SoilDepthMidLay_vr(L,NY,NX))
            THETWM   = THETWP*POROS_vr(L,NY,NX)          !threshold air-filled volume
            THETW1   = AMIN1(THETWM,EXP((LOGPSIAtSat(NY,NX)-LOG(-PSIEquil)) &
              *PSD_vr(L,NY,NX)/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))

            IF(THETWM.GT.THETW1)THEN
              THETPX                 = AMIN1(1.0_r8,AZMAX1((THETWM-THETW_vr(L,NY,NX))/(THETWM-THETW1)))
              DepzIntWTBL_col(NY,NX) = CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)*(1.0_r8-THETPX)     !internal water depth depth relative to current layer upper depth
            ELSE
              DepzIntWTBL_col(NY,NX) = CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)                     !set to upper edge of current layer
            ENDIF
          ELSEIF(L.GT.NU_col(NY,NX))THEN
            !not bottom layer, and not topsoil layer, partially saturated
            PSIEquil = PSISoilMatricP_vr(L,NY,NX)-mGravAccelerat*(SoilDepthMidLay_vr(L,NY,NX)-SoilDepthMidLay_vr(L-1,NY,NX))
            THETWM   = THETWP*POROS_vr(L-1,NY,NX)
            THETW1   = AMIN1(THETWM,EXP((LOGPSIAtSat(NY,NX)-LOG(-PSIEquil)) &
              *PSD_vr(L-1,NY,NX)/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L-1,NY,NX)))
            IF(THETWM.GT.THETW1)THEN
              THETPX                 = AMIN1(1.0_r8,AZMAX1((THETWM-THETW_vr(L-1,NY,NX))/(THETWM-THETW1)))
              DepzIntWTBL_col(NY,NX) = CumDepz2LayBottom_vr(L-1,NY,NX)-DLYR_3D(3,L-1,NY,NX)*(1.0_r8-THETPX)
            ELSE
              DepzIntWTBL_col(NY,NX)=CumDepz2LayBottom_vr(L-1,NY,NX)-DLYR_3D(3,L-1,NY,NX)
            ENDIF
          ELSE
            DepzIntWTBL_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)
          ENDIF
        ENDIF
      ENDIF      
    ENDIF
  END DO
  end subroutine DiagWaterTBLDepz
end module HydrologyDiagMod  