module InitVegBGC

  use EcosimConst
  use GridConsts
  implicit none
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: InitIrradianceGeometry
  contains


!------------------------------------------------------------------------------------------
  subroutine InitIrradianceGeometry(YSIN,YCOS,SkyAzimuthAngle)
  use CanopyRadDataType
  implicit none


  real(r8), intent(out) :: YSIN(NumOfSkyAzimuthSectors)
  real(r8), intent(out) :: YCOS(NumOfSkyAzimuthSectors)
  real(r8), intent(out) :: SkyAzimuthAngle(NumOfSkyAzimuthSectors)
  real(r8) :: ZAZI(NumOfLeafAzimuthSectors)
  real(r8) :: YAGL
  real(r8) :: OMEGZ,OMEGY
  real(r8) :: ZAGL,DAZI

  integer :: L,M,N

  write(*,*) "In InitIrradianceGeometry"

  !     begin_execution
  !     SineLeafAngle,CosineLeafAngle=sine,cosine of leaf inclination class
  !     ZAZI=leaf azimuth class, it is pi because only on one side of the leaf
  !     SkyAzimuthAngle,YSIN,YCOS=sky azimuth,sine,cosine of sky azimuth
  !     OMEGA,OMEGX=incident aNGLe of diffuse radn at leaf,horizontal surface
  !     IALBY:1=backscattering,2=forward scattering of sky radiation
  !

  write(*,*) "First do loop over leaf azimuth sectors"
  D205: DO L=1,NumOfLeafAzimuthSectors
    write(*,*) "Setting ZAZI"
    write(*,*) "L = ", L, " of ", NumOfLeafAzimuthSectors
    ZAZI(L)=(L-0.5)*PICON/real(NumOfLeafAzimuthSectors,r8)
  ENDDO D205
  write(*,*) "Loop over sky azimuth sectors"
  !NumOfSkyAzimuthSectors: number of sky azimuth sectors
  !NumOfLeafAzimuthSectors: number of leaf azimuth sectors
  D230: DO N=1,NumOfSkyAzimuthSectors
    write(*,*) "Setting Yvars"
    write(*,*) "N = ", N, " of ", NumOfSkyAzimuthSectors
    SkyAzimuthAngle(N)=PICON*(2*N-1)/real(NumOfSkyAzimuthSectors,r8)
    YAGL=PICON/real(NumOfSkyAzimuthSectors,r8)
    YSIN(N)=SIN(YAGL)
    YCOS(N)=COS(YAGL)
    TYSIN=TYSIN+YSIN(N)
    D225: DO L=1,NumOfLeafAzimuthSectors
      write(*,*) "Setting DAZI"
      write(*,*) "L = ", L, " of ", NumOfLeafAzimuthSectors
      DAZI=COS(ZAZI(L)-SkyAzimuthAngle(N))
      DO  M=1,NumOfLeafZenithSectors
        write(*,*) "Setting Omega vars"
        write(*,*) "M = ", M, " of ", NumOfLeafZenithSectors
        OMEGY=CosineLeafAngle(M)*YSIN(N)+SineLeafAngle(M)*YCOS(N)*DAZI
        OMEGA(N,M,L)=ABS(OMEGY)
        OMEGX(N,M,L)=OMEGA(N,M,L)/YSIN(N)
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
        IF(ZAGL.GT.0.0_r8.AND.ZAGL.LT.PICON)THEN
          IALBY(N,M,L)=1
        ELSE
          IALBY(N,M,L)=2
        ENDIF
      ENDDO
    ENDDO D225
  ENDDO D230
  end subroutine InitIrradianceGeometry

end module InitVegBGC
