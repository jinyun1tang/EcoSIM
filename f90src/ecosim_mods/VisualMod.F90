module VisualMod
!!
! Description:
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use fileUtil, only : open_safe
  implicit none

  private
  include "parameters.h"
  include "filec.h"
  include "files.h"
  include "blkc.h"
  include "blk1cp.h"
  include "blk1cr.h"
  include "blk1g.h"
  include "blk1n.h"
  include "blk1p.h"
  include "blk2a.h"
  include "blk2b.h"
  include "blk2c.h"
  include "blk3.h"
  include "blk5.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk9a.h"
  include "blk9b.h"
  include "blk9c.h"
  include "blk11a.h"
  include "blk11b.h"
  include "blk12a.h"
  include "blk12b.h"
  include "blk13a.h"
  include "blk13b.h"
  include "blk13c.h"
  include "blk14.h"
  include "blk15a.h"
  include "blk15b.h"
  include "blk16.h"
  include "blk18a.h"
  include "blk18b.h"

  character(len=*), parameter :: mod_filename = __FILE__
  real(r8) :: TCSNX,TTRN,TTLE,TTSH,TTGH,TTCO,TTCH

  integer :: L,NX,NY,N

  real(r8) :: OUT(100),SWC(JY,JX)
  real(r8), PARAMETER :: DEFAULT=-9999

  DATA ICHECK/0/

  integer, SAVE :: IYR1,IYR2,ICHECK

  public :: visual
  contains

  SUBROUTINE visual(I,J,NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

! begin_execution
  IF(ICHECK.EQ.0)THEN
    call OPEN_safe(16,PREFIX,'years','OLD',mod_filename,__LINE__)
    READ(16,*)IYR1,IYR2
!   WRITE(*,*)IYR1,IYR2
    CLOSE(16)
!   DO 8995 NX=NHW,NHE
!   DO 8990 NY=NVN,NVS
    DO 8995 NX=1,1
      DO 8990 NY=1,1
!     WRITE(19,90)'Model_name     ','Year           '
!    2,'Day            ','Hour           ','cum GPP  Spruce'
!    3,'cum GPP shrub  ','cum GPP moss   ','cum Ra Spruce  '
!    4,'cum Ra shrub   ','cum Ra moss    ','cum GPP Total  '
!    5,'cum NPP Total  ','cum Ra Total   ','cum Rh         '
!    6,'cum NEP        ','CO2 Flux       ','CH4 Flux       '
!    7,'Rn             ','LE             ','H              '
!    8,'SWC 0          ','SWC 1          ','SWC 2          '
!    8,'SWC 3          ','SWC 4          ','SWC 5          '
!    8,'SWC 6          ','SWC 7          ','SWC 8          '
!    8,'SWC 9          ','SWC 10         ','Temp 0         '
!    8,'Temp 1         ','Temp 2         ','Temp 3         '
!    8,'Temp 4         ','Temp 5         ','Temp 6         '
!    8,'Temp 7         ','Temp 8         ','Temp 9         '
!    8,'Temp 10        '
!     WRITE(19,90)'-              ','-              '
!    2,'-              ','-              ','kg C m-2 s-1   '
!    3,'kg C m-2 s-1   ','kg C m-2 s-1   ','kg C m-2 s-1   '
!    4,'kg C m-2 s-1   ','kg C m-2 s-1   ','kg C m-2 s-1   '
!    5,'kg C m-2 s-1   ','kg C m-2 s-1   ','kg C m-2 s-1   '
!    6,'kg C m-2 s-1   ','kg C m-2 s-1   ','kg C m-2 s-1   '
!    7,'W m-2          ','W m-2          ','W m-2          '
!    8,'m3 m-3         ','m3 m-3         ','m3 m-3         '
!    8,'m3 m-3         ','m3 m-3         ','m3 m-3         '
!    8,'m3 m-3         ','m3 m-3         ','m3 m-3         '
!    8,'m3 m-3         ','m3 m-3         ','m3 m-3         '
!    8,'m3 m-3         ','m3 m-3         ','K              '
!    8,'K              ','K              ','K              '
!    8,'K              ','K              ','K              '
!    8,'K              ','K              ','K              '
!    8,'K              '
!     WRITE(20,95)'Model_name     ','Year           '
!    2,'Day            ','cum Trspn Spruc','cum Trspn Shrub'
!    3,'cum Trspn Moss ','cum Evapotrnspn','cum Soil Evapn '
!    4,'cum Runoff     ','cum Disch/Rech ','Water Table    '
!    5,'Snowpack       ','LAI Total      ','Leaf Spruce    '
!    6,'Leaf Shrub     ','Leaf Moss      ','Wood Spruce    '
!    7,'Wood Shrub     ','Wood Moss      ','Root Spruce    '
!    8,'Root Shrub     ','Root Moss      ','SOC 0          '
!    9,'SOC 1          ','SOC 2          ','SOC 3          '
!    9,'SOC 4          ','SOC 5          ','SOC 6          '
!    9,'SOC 7          ','SOC 8          ','SOC 9          '
!    9,'SOC 10         '
!     WRITE(20,95)'-              ','-              '
!    2,'               ','mm             ','mm             '
!    3,'mm             ','mm             ','mm             '
!    4,'mm             ','mm             ','m              '
!    5,'m              ','m2 m-2         ','kg C m-2       '
!    6,'kg C m-2       ','kg C m-2       ','kg C m-2       '
!    6,'kg C m-2       ','kg C m-2       ','kg C m-2       '
!    6,'kg C m-2       ','kg C m-2       ','kg C m-2       '
!    6,'kg C m-2       ','kg C m-2       ','kg C m-2       '
!    6,'kg C m-2       ','kg C m-2       ','kg C m-2       '
!    6,'kg C m-2       ','kg C m-2       ','kg C m-2       '
!    6,'kg C m-2       '
8990  CONTINUE
8995  CONTINUE
    TCSNX=0.0
    ICHECK=1
  ENDIF
!
! SELECT YEARS
!
  IF(IYRC.GE.IYR1.AND.IYRC.LE.IYR2)THEN
    TTRN=0.0
    TTLE=0.0
    TTSH=0.0
    TTGH=0.0
    TTCO=0.0
    TTCH=0.0
    DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
        TTRN=TTRN+TRN(NY,NX)
        TTLE=TTLE+TLE(NY,NX)
        TTSH=TTSH+TSH(NY,NX)
        TTGH=TTGH+TGH(NY,NX)
        TTCO=TTCO+TCNET(NY,NX)
        TTCH=TTCH+HCH4G(NY,NX)
        IF(J.EQ.24)THEN
          IF(NU(NY,NX).EQ.7)THEN
            SWC(NY,NX)=(VOLW(8,NY,NX)+AMIN1(VOLAH(8,NY,NX) &
              ,VOLWH(8,NY,NX)))/VOLT(8,NY,NX)
          ELSEIF(NU(NY,NX).EQ.4)THEN
            SWC(NY,NX)=(VOLW(5,NY,NX)+AMIN1(VOLAH(5,NY,NX) &
              ,VOLWH(5,NY,NX)))/VOLT(5,NY,NX)
          ELSE
            SWC(NY,NX)=(VOLW(2,NY,NX)+AMIN1(VOLAH(2,NY,NX) &
              ,VOLWH(2,NY,NX)))/VOLT(2,NY,NX)
          ENDIF
        ENDIF
!     GO TO 9990
!     NZ=1
!
!     HOURLY AGGREGATION
!
!     IF(IFLGC(NZ,NY,NX).EQ.1)THEN
!     IF(I.EQ.IDAY0(NZ,NY,NX))THEN
!     ND=0
!     TCCTX=0.0
!     TCCMX=0.0
!     ENDIF
!     IF(J.EQ.1)THEN
!     TCCTD=0.0
!     TCCMD=-50.0
!     ENDIF
!     TCCTD=TCCTD+TCC(NZ,NY,NX)
!     TCCMD=AMAX1(TCCMD,TCC(NZ,NY,NX))
!     IF(J.EQ.24)THEN
!     TCCAD=TCCTD/24
!     TCCTX=TCCTX+TCCAD
!     TCCMX=TCCMX+TCCMD
!     ND=ND+1
!     ENDIF
!     IF(I.EQ.IDAYH(NZ,NY,NX)-1)THEN
!     TCCAS=TCCTX/ND
!     TCCMS=TCCMX/ND
!     ENDIF
!     ELSE
!     TCCAD=DEFAULT
!     TCCMD=DEFAULT
!     TCCAS=DEFAULT
!     TCCMS=DEFAULT
!     ENDIF
!
!     DAILY AND SEASONAL OUTPUT
!
!     IF(J.EQ.24)THEN
!     IYRZ=IYRC
!     XI=REAL(I)
!     XP=REAL(IDAY0(NZ,NY,NX))
!     XH=REAL(IDAYH(NZ,NY,NX))
!     IF(I.EQ.IDAY0(NZ,NY,NX))THEN
         DO 9980 N=1,100
           OUT(N)=0.0
9980     CONTINUE
!     ENDIF
!     IF(I.EQ.1)THEN
!     TCSNY=0.0
!     ENDIF
!     DTCSN=TCSN0(NZ,NY,NX)-TCSNY
!     TCSNX=TCSNX+DTCSN
!     TCSNY=TCSN0(NZ,NY,NX)
!     IF(I.EQ.IDAY0(NZ,NY,NX).OR.I.EQ.IDAYH(NZ,NY,NX))THEN
!     ICHKA=0
!     ICHKM=0
!     TCSNX=0.0
!     ENDIF
          OUT(1)=0.001*CARBN(1,NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(2)=0.001*CARBN(3,NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(3)=0.001*CARBN(2,NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(4)=-0.001*TCO2T(1,NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(5)=-0.001*TCO2T(3,NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(6)=-0.001*TCO2T(2,NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(7)=0.001*TGPP(NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(8)=0.001*TNPP(NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(9)=-0.001*TRAU(NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(10)=-0.001*THRE(NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(11)=0.001*TNBP(NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(12)=-0.001*HCO2G(NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(13)=-0.001*HCH4G(NY,NX)/(AREA(3,NU(NY,NX),NY,NX)*3600)
          OUT(14)=TRN(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
          OUT(15)=-TLE(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
          OUT(16)=-TSH(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
          L=1
          DO 60 N=17,27
            OUT(N)=(VOLW(L,NY,NX)+AMIN1(VOLAH(L,NY,NX) &
              ,VOLWH(L,NY,NX)))/VOLT(L,NY,NX)
            L=L+1
60        CONTINUE
          L=1
          DO 61 N=28,38
            OUT(N)=TKS(L,NY,NX)
            L=L+1
61        CONTINUE
          OUT(39)=-1000.0*CTRAN(1,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(40)=-1000.0*CTRAN(3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(41)=-1000.0*CTRAN(2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(42)=1000.0*UEVAP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(43)=OUT(42)-OUT(39)-OUT(40)-OUT(41)
          OUT(44)=1000.0*URUN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(45)=1000.0*UVOLO(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(46)=-(DPTHT(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
          OUT(47)=DPTHS(NY,NX)
          OUT(48)=ARLFC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(49)=0.001*WTLF(1,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(50)=0.001*WTLF(3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(51)=0.001*WTLF(2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(52)=0.001*WTSTK(1,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(53)=0.001*WTSTK(3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(54)=0.001*WTSTK(2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(55)=0.001*WTRT(1,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(56)=0.001*WTRT(3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          OUT(57)=0.001*WTRT(2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
          L=0
          DO 62 N=58,68
            OUT(N)=0.001*ORGC(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
            L=L+1
62        CONTINUE
!
!     WRITE OUTPUT
!
!     WRITE(19,91)'ECOSYS_HUMMOCK',IYRC,I,J,(OUT(N),N=1,38)
!     IF(J.EQ.24)THEN
!     WRITE(20,96)'ECOSYS_HUMMOCK',IYRC,I,(OUT(N),N=39,68)
!     ENDIF
!     ENDIF
90        FORMAT(50A16)
95        FORMAT(50A16)
91        FORMAT(A16,3I16,38E16.6)
96        FORMAT(A16,2I16,30E16.6)
9990    CONTINUE
9995  CONTINUE
!
!     WRITE LANDSCAPE OUTPUT
!
    WRITE(19,2025)'FLUXES',IYRC,I,J,TTRN*277.8/TAREA &
      ,TTLE*277.8/TAREA,TTSH*277.8/TAREA,TTGH*277.8/TAREA &
      ,TTCO*23.14815/TAREA,TTCH*23.14815/TAREA &
      ,((TCNET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815 &
      ,NX=NHW,NHE),NY=NVN,NVS),DEFAULT &
      ,((HCH4G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815 &
      ,NX=NHW,NHE),NY=NVN,NVS)
2025  FORMAT(A16,3I6,100E12.4)
    IF(J.EQ.24)THEN
      WRITE(20,2026)'SWC',IYRC,I,J,((DPTHS(NY,NX) &
      ,NX=NHW,NHE),NY=NVN,NVS),DEFAULT &
      ,((SWC(NY,NX),NX=NHW,NHE),NY=NVN,NVS),DEFAULT &
      ,((-(DPTHA(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX)) &
      ,NX=NHW,NHE),NY=NVN,NVS)
2026  FORMAT(A8,3I6,100E12.4)
    ENDIF
  ENDIF
  RETURN

  END subroutine visual
end module VisualMod