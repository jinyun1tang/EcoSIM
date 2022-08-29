module InitVegBGC

  use EcosimConst
  use GridConsts
  implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: InitIrradianceGeometry
  contains


!------------------------------------------------------------------------------------------
  subroutine InitIrradianceGeometry(YSIN,YCOS,YAZI)
  use CanopyRadDataType
  implicit none


  real(r8), intent(out) :: YSIN(JSA)
  real(r8), intent(out) :: YCOS(JSA)
  real(r8), intent(out) :: YAZI(JSA)
  real(r8) :: ZAZI(JLA)
  real(r8) :: YAGL
  real(r8) :: OMEGZ,OMEGY
  real(r8) :: ZAGL,DAZI

  integer :: L,M,N
  !     begin_execution
  !     ZSIN,ZCOS=sine,cosine of leaf inclination class
  !     ZAZI=leaf azimuth class, it is pi because only on one side of the leaf
  !     YAZI,YSIN,YCOS=sky azimuth,sine,cosine of sky azimuth
  !     OMEGA,OMEGX=incident aNGLe of diffuse radn at leaf,horizontal surface
  !     IALBY:1=backscattering,2=forward scattering of sky radiation
  !

  DO 205 L=1,JLA
    ZAZI(L)=(L-0.5)*PICON/real(JLA,r8)
205   CONTINUE
  !JSA: number of sky azimuth sectors
  !JLA: number of leaf azimuth sectors
  DO 230 N=1,JSA
    YAZI(N)=PICON*(2*N-1)/real(JSA,r8)
    YAGL=PICON/real(JSA,r8)
    YSIN(N)=SIN(YAGL)
    YCOS(N)=COS(YAGL)
    TYSIN=TYSIN+YSIN(N)
    DO 225 L=1,JLA
      DAZI=COS(ZAZI(L)-YAZI(N))
      DO  M=1,JLI
        OMEGY=ZCOS(M)*YSIN(N)+ZSIN(M)*YCOS(N)*DAZI
        OMEGA(N,M,L)=ABS(OMEGY)
        OMEGX(N,M,L)=OMEGA(N,M,L)/YSIN(N)
        IF(ZCOS(M).GT.YSIN(N))THEN
          OMEGZ=ACOS(OMEGY)
        ELSE
          OMEGZ=-ACOS(OMEGY)
        ENDIF
        IF(OMEGZ.GT.-PICON2)THEN
          ZAGL=YAGL+2.0*OMEGZ
        ELSE
          ZAGL=YAGL-2.0*(PICON+OMEGZ)
        ENDIF
        IF(ZAGL.GT.0.0.AND.ZAGL.LT.PICON)THEN
          IALBY(N,M,L)=1
        ELSE
          IALBY(N,M,L)=2
        ENDIF
      ENDDO
225 CONTINUE
230 CONTINUE
  end subroutine InitIrradianceGeometry

end module InitVegBGC
