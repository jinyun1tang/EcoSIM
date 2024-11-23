module InitVegBGC

  use EcosimConst
  use GridConsts
  use GrosubPars, only : ibackward,iforward
  implicit none
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: InitIrradianceGeometry
  contains


!------------------------------------------------------------------------------------------
  subroutine InitIrradianceGeometry(YSIN,YCOS,SkyAzimuthAngle)
  use CanopyRadDataType
  implicit none

  real(r8), intent(out) :: YSIN(NumOfSkyAzimuthSects)
  real(r8), intent(out) :: YCOS(NumOfSkyAzimuthSects)
  real(r8), intent(out) :: SkyAzimuthAngle(NumOfSkyAzimuthSects)
  real(r8) :: ZAZI(NumOfLeafAzimuthSectors)
  real(r8) :: YAGL    !sector angle of sky azimuth
  real(r8) :: OMEGZ,OMEGY
  real(r8) :: ZAGL,DAZI

  integer :: L,M,N

  !     begin_execution
  !     SineLeafAngle,CosineLeafAngle=sine,cosine of leaf inclination class
  !     ZAZI=leaf azimuth class, it is pi because only on one side of the leaf
  !     SkyAzimuthAngle,YSIN,YCOS=sky azimuth,sine,cosine of sky azimuth
  !     OMEGA,OMEGX=incident aNGLe of diffuse radn at leaf,horizontal surface
  !     iScatteringDiffus:1=backscattering,2=forward scattering of sky radiation
  !

  D205: DO L=1,NumOfLeafAzimuthSectors
    ZAZI(L)=(L-0.5)*PICON/real(NumOfLeafAzimuthSectors,r8)
  ENDDO D205
  !NumOfSkyAzimuthSects: number of sky azimuth sectors
  !NumOfLeafAzimuthSectors: number of leaf azimuth sectors
  D230: DO N=1,NumOfSkyAzimuthSects
    SkyAzimuthAngle(N)   = PICON*(2*N-1)/real(NumOfSkyAzimuthSects,r8)
    YAGL                 = PICON/real(NumOfSkyAzimuthSects,r8)
    YSIN(N)              = SIN(YAGL)
    YCOS(N)              = COS(YAGL)
    TotSineSkyAngles_grd = TotSineSkyAngles_grd+YSIN(N)
    D225: DO L=1,NumOfLeafAzimuthSectors
      DAZI=COS(ZAZI(L)-SkyAzimuthAngle(N))
      DO  M=1,NumOfLeafZenithSectors
        OMEGY        = CosineLeafAngle(M)*YSIN(N)+SineLeafAngle(M)*YCOS(N)*DAZI
        OMEGA(N,M,L) = ABS(OMEGY)
        OMEGX(N,M,L) = OMEGA(N,M,L)/YSIN(N)
        IF(CosineLeafAngle(M).GT.YSIN(N))THEN
          OMEGZ=ACOS(OMEGY)
        ELSE
          OMEGZ=-ACOS(OMEGY)
        ENDIF
        IF(OMEGZ.GT.-PICON2h)THEN
          ZAGL=YAGL+2.0_r8*OMEGZ
        ELSE
          ZAGL=YAGL-2.0_r8*(PICON+OMEGZ)
        ENDIF
        IF(ZAGL.GT.0.0_r8 .AND. ZAGL.LT.PICON)THEN
          iScatteringDiffus(N,M,L)=ibackward
        ELSE
          iScatteringDiffus(N,M,L)=iforward
        ENDIF
      ENDDO
    ENDDO D225
  ENDDO D230
  end subroutine InitIrradianceGeometry

end module InitVegBGC
