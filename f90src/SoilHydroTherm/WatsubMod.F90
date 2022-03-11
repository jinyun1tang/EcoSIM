module WatsubMod
!!
! Description:
! code to do water and enerby balance calculation

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils   , only : endrun, print_info
  use minimathmod  , only : test_aeqb, test_aneb,safe_adb,vapsat
  use EcosimConst
  use SOMDataType
  use WatsubPars
  use GridDataType
  implicit none

  private

  include "blkc.h"
  include "blk2a.h"
  include "blk2b.h"
  include "blk2c.h"
  include "blk5.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk10.h"
  include "blk11a.h"
  include "blk11b.h"
  include "blk13a.h"
  include "blk13b.h"
  include "blk13c.h"
  include "blk15a.h"
  include "blk15b.h"
  include "blk22a.h"
  include "blk22b.h"
  include "blk22c.h"
  include "blktest.h"

  character(len=*), parameter :: mod_filename = __FILE__

  real(r8) :: VOLWX1(JZ,JY,JX),VPQ(JY,JX),TKQ(JY,JX) &
    ,XVOLT(JY,JX),XVOLW(JY,JX),XVOLI(JY,JX),FMAC(JZ,JY,JX) &
    ,FGRD(JZ,JY,JX),VOLW1(0:JZ,JY,JX),VOLI1(0:JZ,JY,JX) &
    ,VHCP1(0:JZ,JY,JX),VHCP1A(JZ,JY,JX),VHCP1B(JZ,JY,JX) &
    ,TK1(0:JZ,JY,JX),TWFLXL(JZ,JY,JX),VOLW2(JZ,JY,JX) &
    ,VOLP1(0:JZ,JY,JX),TWFLXH(JZ,JY,JX),PRECM(JY,JX) &
    ,VOLS0(JS,JY,JX),VOLI0(JS,JY,JX),VOLW0(JS,JY,JX) &
    ,VOLS1(JS,JY,JX),DLYRS0(JS,JY,JX),VOLP1Z(JZ,JY,JX) &
    ,TK0(JS,JY,JX),AREAU(JZ,JY,JX),AREAUD(JZ,JY,JX),FLQ0S(JY,JX) &
    ,FLQ0I(JY,JX),FLQ0W(JY,JX),FLQ1(JY,JX),FLH1(JY,JX) &
    ,FLY1(JY,JX),HWFLQ0(JY,JX),HWFLQ1(JY,JX),HWFLY1(JY,JX) &
    ,RAR(JY,JX),RAGS(JY,JX),BAREW(JY,JX),CVRDW(JY,JX) &
    ,RARG(JY,JX),RAGR(JY,JX),RAGW(JY,JX),PAREW(JY,JX),PAREG(JY,JX) &
    ,RAG(JY,JX),PARSW(JY,JX),PARSG(JY,JX),PARER(JY,JX),PARSR(JY,JX) &
    ,QR1(2,2,JV,JH),HQR1(2,2,JV,JH),VOLPH1Z(JZ,JY,JX) &
    ,QS1(2,JV,JH),QW1(2,JV,JH),QI1(2,JV,JH),HQS1(2,JV,JH) &
    ,TQR1(JY,JX),THQR1(JY,JX),EVAPG(JY,JX),TTFLXL(JZ,JY,JX) &
    ,EVAPW(JY,JX),EVAPS(JY,JX),EVAPR(JY,JX),FLWRL(JY,JX) &
    ,HFLWRL(JY,JX),FINHL(JZ,JY,JX),FLWL(3,JD,JV,JH)
  real(r8) :: FLWHL(3,JD,JV,JH),HFLWL(3,JD,JV,JH) &
    ,TFLWL(JZ,JY,JX),TFLWHL(JZ,JY,JX),THFLWL(JZ,JY,JX) &
    ,WFLXL(JZ,JY,JX),TFLXL(JZ,JY,JX),AVCNHL(3,JD,JV,JH) &
    ,THRYW(JY,JX),THRMW(JY,JX),THRMS(JY,JX),THRMR(JY,JX) &
    ,THRYG(JY,JX),THRYR(JY,JX),RADXW(JY,JX),RADXG(JY,JX) &
    ,RADXR(JY,JX),FLWLX(3,JD,JV,JH),TFLWLX(JZ,JY,JX) &
    ,FLU1(JZ,JY,JX),HWFLU1(JZ,JY,JX),PSISM1(0:JZ,JY,JX) &
    ,ALTG(JY,JX),WFLXLH(JZ,JY,JX),DLYRR(JY,JX),WFLXR(JY,JX) &
    ,TFLXR(JY,JX),HCNDR(JY,JX),CNDH1(JZ,JY,JX),VOLA1(0:JZ,JY,JX) &
    ,THETWX(0:JZ,JY,JX),THETIX(0:JZ,JY,JX),THETPX(0:JZ,JY,JX) &
    ,VOLAH1(JZ,JY,JX),VOLWH1(JZ,JY,JX),VOLPH1(JZ,JY,JX) &
    ,VOLIH1(JZ,JY,JX),THETPY(0:JZ,JY,JX),FLWNX(JY,JX) &
    ,FLWXNX(JY,JX),FLWHNX(JY,JX),HFLWNX(JY,JX),PSISA1(JZ,JY,JX)
  real(r8) :: TQS1(JY,JX),TQW1(JY,JX),TQI1(JY,JX),THQS1(JY,JX) &
    ,TFLWS(JS,JY,JX),TFLWW(JS,JY,JX),TFLWI(JS,JY,JX) &
    ,THFLWW(JS,JY,JX),TFLX0(JS,JY,JX),WFLXS(JS,JY,JX) &
    ,WFLXI(JS,JY,JX),FLW0W(JS,JY,JX),FLW0S(JS,JY,JX) &
    ,FLW0I(JS,JY,JX),HFLW0W(JS,JY,JX) &
    ,VOLS0M(JS,JY,JX),VOLW0M(JS,JY,JX) &
    ,VOLI0M(JS,JY,JX),VHCPWMM(JS,JY,JX),TK0M(JS,JY,JX)

  real(r8) :: ALFZ,ALBW,ATCNVW,ATCNDW,ATCNVS,ATCNDS,ATCNVR,ATCNDR
  real(r8) :: ALBG,ALBR,AVCNDR,ALT1,ALT2,ALTB,ALTS1,ALTS2,AVCNDL
  real(r8) :: ATCNVL,ATCNDL,CNV1,CNV2,CNVR,CNDR,CND1,CNDL,CNVL
  real(r8) :: DFVR,DENSW1,DENSW2,DTKX,DTHW0,DTHA0,DTHW1,DTHA1
  real(r8) :: D,DTHW2,DTHA2,DTBLXX,DPTHH,DTBLYX,DPTHW1,DPTHW2
  real(r8) :: ENGYD,ENGYB,EFLXW,EVAPT2,EVAPW2,EVAPX2,EVAPS2
  real(r8) :: EFLXW2,ENGY0,EFLXG,EFLXR,EVAPR2,EFLXR2,ENGYR,ENGY1
  real(r8) :: FLWQW,FLWSW,FLWQB,FLWQBX,FLWQAX,FLWQAS,FLWQAH
  real(r8) :: FLQRS,FLQRH,FLYM,FLQM,FLHM,FLYM2,FKSAT,FLWLW,FLWLXW
  real(r8) :: FLWHLW,FLWRLW,FLWVLW,FLQ0S2,FLQ0W2,FLQ0I2,FLW0S2
  real(r8) :: FLW0W2,FLW0I2,FLWQX,FLWQM,FLVC,FLVX,FLVSS,FLW0T
  real(r8) :: FLWQGX,FLWQGS,FLWQGH,FLWQG,FLWQR,FCX,FCLX,FCDX
  real(r8) :: FLVS1,FLVSR,FLVR1,FLVSRX,FLVR1X,FLWLT,FLWRT,FVOLS0
  real(r8) :: FVOLI0,FLV1,FLV2,FLWLG,FLWLXG,FLWHLG,FLWRLG,FLWVLS
  real(r8) :: FLQX,FLQZ,FLQR,FLQ2,FLQHR,FLQL,FLWHX,FLVL,FLWT
  real(r8) :: FLWTH,FLWTHL,FLWU,FLWUL,FLWUX,FLWUH,FLWUHL,FINHX
  real(r8) :: FLWS,FLWW,FLWI,HFLWSW,HFLWQB,HFLWQA,HFLQR1,HWFLYM
  real(r8) :: HWFLQM,HWFLM2,HV,HFLXW,HFLWLW,HFLWRLW,HFLX02,HFLXW2
  real(r8) :: HWFLQ02,HFLW0W2,HFLWQM,HFLVSS,HFLWX,HFLWC,HFLWSS
  real(r8) :: HFLW0T,HFLWQG,HFLWQR,HFLVS1,HFLWS1,HFLVSR,HFLWSR
  real(r8) :: HFLVR1,HFLWR1,HFLVSRX,HFLWSRX,HFLVR1X,HFLWR1X
  real(r8) :: HFLWLT,HFLWRT,HFLX0,HFLXG,HFLXR,HWFLV1,HFLCR1
  real(r8) :: HWFLV2,HFLCR2,HFLXR2,HFLWLG,HFLWRLG,HFLQR,HFLQHR
  real(r8) :: HWFLQL,HWFLHL,HWFLVL,HWFLWL,HFLWSX,HFLWS,PAREX
  real(r8) :: PARSX,PARE,PARS,PSDX,PSISV1,PSIST0,PSIST1,PSISTL
  real(r8) :: PSISVL,PSISH1,PSISHL,PSISWD,PSISWT,PSISWTH,PSISUT
  real(r8) :: PSISUTH,PSISMX,Q,QRQ1,QSX,RADGX,RAR1,RAS,RASX
  real(r8) :: RASL,RFLXW,RFLX0,RFLXW2,RAGX,RA,RFLXG,RFLXR,RYLXW0
  real(r8) :: RYLXA0,RYLXW1,RYLXA1,RYLNW0,RYLNA0
  real(r8) :: RYLNW1,RYLNA1,RFLXR2,R,RYLXW2,RYLXA2,RYLNW2,RYLNA2
  real(r8) :: RCHQF,RCHGFU,RCHGFT,SFLXW,SFLXW2,SFLXG,SFLXR,SFLXR2
  real(r8) :: SS,TVOLWI,THRYX,THETPX0,THETPL,THETWA,THETWT,THETWH
  real(r8) :: THRMX,THETP1,TCND1W,THETP2,TCND2W,TCNDS,TKWX1
  real(r8) :: THETRR,TCNDR,TFLWSX,TFLWWX,TFLWIX,THFLWWX,TFLX1
  real(r8) :: TVOLWS,TFLX0X,THRMA,THRMZ,TCNDW0,TCNDA0,TCNDW1
  real(r8) :: TCNDA1,TCND1,THRMZ2,TFREEZ,TFLX,TCNDW2,TCNDA2
  real(r8) :: TCND2,TKLX,THETAX,TFLXH1,TFLXH,TK0XX,TKXX,UAG
  real(r8) :: VOLTX,VOLWRZ,VOLIRZ,VOLAT0,VOLWT,VOLAT,VFLXW,VP0
  real(r8) :: VFLXW2,VOLP01,VP1,VOLP02,VP2,VPY,VPR,VOLS0X,VOLW0X
  real(r8) :: VOLI0X,VHCPWMX,VOLWXG,VOLIXG,VFLXG,VFLXR,VOLWR2
  real(r8) :: VHCPR2,VOLW12,VHCP12,VFLXR2,VOLW1X,VHCP1X,VX,V
  real(r8) :: VPL,VOLP2,VOLPX2,VOLPH2,VOLP1X,VOLWH1X,VOLPH1X
  real(r8) :: VHCP1AX,VHCP1BX,VHCPXX,VHXX,WPX,WPLX,WTHET2,WFLXSX
  real(r8) :: WFLXIX,WTHET0,WTHET1,XNUSW0,XNUSA0,XNUSW1,XNUSA1
  real(r8) :: XNUSW2,XNUSA2,XN,Z3S

  integer :: ICHKL,IFLGH,IFLGU,IFLGUH,IFLGD,IFLGDH
  integer :: N6X(JY,JX)

  REAL(r8) :: RI,THETWR,THETW1,THETA1,THETAL,THETWL &
    ,TKR1,TKS1,TKY,TKW1,TK11,TK12,TK0X,TKXR,TK1X,TKX1,TFND1
  integer :: curday,curhour
  public :: watsub,InitWatsub
  contains

  subroutine InitWatsub

  implicit none

  call InitWatsubPars

  end subroutine InitWatsub
!------------------------------------------------------------------------------------------

  SUBROUTINE watsub(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CACULATES ENERGY BALANCES OF SNOW, RESIDUE
!     AND SOIL SURFACES, FREEZING, THAWING, AND HEAT AND WATER
!     TRANSFER THROUGH SOIL PROFILES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: K0,K1
  integer :: KL,L,L2,LL,M,MM,M1,M2,M3,M4,M5,M6,NX,NY
  integer :: N,N1,N2,N3,N4,N5,N6,NN,N4B,N5B,NUX

! begin_execution
!
  curday=i
  curhour=j
!  if(curday>=41)call print_info('PrepWaterEnergyBalance',__LINE__)
!  if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NVN,NHW)
  call PrepWaterEnergyBalance(I,J,NHW,NHE,NVN,NVS)
!
! DYNAMIC LOOP FOR FLUX CALCULATIONS
!
  DO 3320 M=1,NPH
    DO 9895 NX=NHW,NHE
      DO 9890 NY=NVN,NVS
!
!        if(curday>=41)call print_info('PrepForIterationM',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NY,NX)
        call PrepForIterationM(M,NY,NX)
!
!        if(curday>=41)call print_info('AtmosLandSurfaceExchange',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NY,NX)
        call AtmosLandSurfaceExchange(M,NY,NX)
!
!     CAPILLARY EXCHANGE OF WATER BETWEEN SOIL SURFACE AND RESIDUE
!        if(curday>=41)call print_info('SurfSoilResidueWaterCapillExch',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NY,NX)
        call SurfSoilResidueWaterCapillExch(M,NY,NX)

!        if(curday>=41)call print_info('InfiltrationRunoffPartition',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NY,NX)
        call InfiltrationRunoffPartition(M,NY,NX,N1,N2)
!
!        if(curday>=41)call print_info('LateralHydroExchange',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NY,NX)
        call LateralHydroExchange(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
!
!        if(curday>=41)call print_info('AccumWaterVaporHeatFluxes',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NY,NX)
        call AccumWaterVaporHeatFluxes(M,NY,NX)

!        if(curday>=41)call print_info('Subsurface3DFlow',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NY,NX)
        call Subsurface3DFlow(M,NY,NX,NHE,NVS)
9890  CONTINUE
9895  CONTINUE
!
!      if(curday>=41)call print_info('BOUNDARY WATER AND HEAT FLUXES',__LINE__)
!
!     XVOLT,XVOLW=excess water+ice,water in source grid cell
!     VOLP2,VOLPH2=air-filled porosity in micropores,macropores
!      if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NVN,NHW)
      call WaterHeatExchThruBoundaryFlow(M,NHW,NHE,NVN,NVS)
!
!      if(curday==41)call print_info('UPDATE STATE VARIABLES FROM FLUXES CALCULATED ABOVE',__LINE__)
!      if(curday>=175)then
!        write(*,*)'at line',__LINE__,'TK1(8,NVN,NHW)',TK1(8,NVN,NHW),M
!        if(TK1(8,NVN,NHW)>1.e3)then
!          call endrun(trim(mod_filename)//' at line',__LINE__)
!        endif
!      endif
      IF(M.NE.NPH)THEN
!
!        if(curday>=175)call print_info('UpdateStateSolutionNotNPH',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NVN,NHW)
        call UpdateStateSolutionNotNPH(M,NHW,NHE,NVN,NVS)
      ELSE
!
!        if(curday>=175)call print_info('UpdateFLuxAtNPH',__LINE__)
!        if(curday>=41)write(*,*)'HFLWRL(NY,NX)=',HFLWRL(NVN,NHW)
        call UpdateFLuxAtNPH(NHW,NHE,NVN,NVS)
      ENDIF
3320  CONTINUE
      RETURN

      END subroutine watsub
!------------------------------------------------------------------------------------------

  subroutine LocalCopyStateVars(I,NY,NX)
  implicit none
  integer, intent(in) :: I,NY,NX

  integer :: L,LWDPTH
! begin_execution
!
! ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
! FOR USE AT INTERNAL TIME STEP
!
! SET INITIAL SNOWPACK VALUES
!
! VOLS0,VOLSSL=snowpack snow content (water equivalent)
! VOLI0,VOLISSL=snowpack ice content
! VOLW0,VOLWSL=snowpack water content
! VOLS1,VOLSL=snowpack volume
! DLYRS0,DLYRS=snowpack depth
! VHCPWM,VHCPW=snowpack heat capacity
! TK0,TKW=snowpack temperature
!
  DO 60 L=1,JS
    VOLS0(L,NY,NX)=VOLSSL(L,NY,NX)
    VOLI0(L,NY,NX)=VOLISL(L,NY,NX)
    VOLW0(L,NY,NX)=VOLWSL(L,NY,NX)
    VOLS1(L,NY,NX)=VOLSL(L,NY,NX)
    DLYRS0(L,NY,NX)=DLYRS(L,NY,NX)
    VHCPWM(1,L,NY,NX)=VHCPW(L,NY,NX)
    TK0(L,NY,NX)=TKW(L,NY,NX)
60  CONTINUE
!
! SET INITIAL SOIL VALUES
!
! WFLXR,TFLXR=initialize surface litter freeze,thaw,latent heat
! CDPTH=depth to bottom of soil layer
! WDPTH,LWDPTH=depth,layer of subsurface irrigation
!
  WFLXR(NY,NX)=0.0_r8
  TFLXR(NY,NX)=0.0_r8
  DO 65 L=NUM(NY,NX),NL(NY,NX)
    IF(CDPTH(L,NY,NX).GE.WDPTH(I,NY,NX))THEN
      LWDPTH=L
      GO TO 55
    ENDIF
65  CONTINUE
55  CONTINUE
  DO 30 L=NUM(NY,NX),NL(NY,NX)
!
!   ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
!   FOR USE AT INTERNAL TIME STEP IN30    con SOIL LAYERS
!
!   PSISM1,PSISM=matric water potential
!   VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of micropores
!   VOLWX1=VOLW1 accounting for wetting front
!   VOLAH*,VOLWH*,VOLIH*,VOLPH*=pore,water,ice,air macropores
!   BKDS=bulk density
!   CCLAY=clay concentration
!   FVOLAH=parameter for clay effect on macropore volume
!   VOLX,VOLT=soil,total volumes
!   WP=wilting point
!   THETW*,THETI*,THETP*=water,ice,air-filled porosity
!   VHCP1,VHCM=volumetric heat capacities of total volume, solid
!   VHCP1A,VHCP1B=volumetric heat capacities of micropore,macropore
!
    PSISM1(L,NY,NX)=PSISM(L,NY,NX)
    VOLA1(L,NY,NX)=VOLA(L,NY,NX)
    VOLW1(L,NY,NX)=VOLW(L,NY,NX)
    VOLWX1(L,NY,NX)=VOLWX(L,NY,NX)
    VOLI1(L,NY,NX)=VOLI(L,NY,NX)
    VOLWH1(L,NY,NX)=VOLWH(L,NY,NX)
    VOLIH1(L,NY,NX)=VOLIH(L,NY,NX)
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLP1Z(L,NY,NX)=VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
      VOLP1(L,NY,NX)=AMAX1(0.0,VOLP1Z(L,NY,NX))
      if(abs(VOLP1Z(L,NY,NX))>1.e20_r8)then
        write(*,*)'aVOLA1=',VOLA1(L,NY,NX),L
        write(*,*)'VOLW1=',VOLW1(L,NY,NX)
        write(*,*)'VOLI1=',VOLI1(L,NY,NX)
      endif
    ELSE
      VOLP1Z(L,NY,NX)=0.0_r8
      VOLP1(L,NY,NX)=0.0_r8
    ENDIF
    VOLAH1(L,NY,NX)=AMAX1(0.0_r8,VOLAH(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX) &
      *(safe_adb(VOLW1(L,NY,NX),VOLY(L,NY,NX))-WP(L,NY,NX))*VOLT(L,NY,NX))
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLPH1Z(L,NY,NX)=VOLAH1(L,NY,NX)-VOLWH1(L,NY,NX)-VOLIH1(L,NY,NX)
      VOLPH1(L,NY,NX)=AMAX1(0.0,VOLPH1Z(L,NY,NX))
    ELSE
      VOLPH1Z(L,NY,NX)=0.0_r8
      VOLPH1(L,NY,NX)=0.0_r8
    ENDIF
    VOLWM(1,L,NY,NX)=VOLW1(L,NY,NX)
    VOLWHM(1,L,NY,NX)=VOLWH1(L,NY,NX)
    VOLPM(1,L,NY,NX)=VOLP1(L,NY,NX)+VOLPH1(L,NY,NX) &
      +THETPI*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
    VOLTX=VOLY(L,NY,NX)+VOLAH1(L,NY,NX)
    IF(VOLTX.GT.ZEROS2(NY,NX))THEN
      THETWX(L,NY,NX)=AMAX1(0.0,(VOLW1(L,NY,NX)+VOLWH1(L,NY,NX))/VOLTX)
      THETIX(L,NY,NX)=AMAX1(0.0,(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))/VOLTX)
      THETPX(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))/VOLTX)
    ELSE
      THETWX(L,NY,NX)=POROS(L,NY,NX)
      THETIX(L,NY,NX)=0.0_r8
      THETPX(L,NY,NX)=0.0_r8
    ENDIF
    THETPM(1,L,NY,NX)=THETPX(L,NY,NX)
    IF(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETPY(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX)) &
        /(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)))
!      if(curday>=41)then
!        write(*,*)'iTHETPY=',THETPY(L,NY,NX),L,VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)
!      endif
    ELSE
      THETPY(L,NY,NX)=0.0_r8
    ENDIF
    VHCP1(L,NY,NX)=VHCM(L,NY,NX)+cpw*(VOLW1(L,NY,NX) &
      +VOLWH1(L,NY,NX))+cpi*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
    VHCP1A(L,NY,NX)=VHCM(L,NY,NX)+cpw*VOLW1(L,NY,NX) &
      +cpi*VOLI1(L,NY,NX)
    VHCP1B(L,NY,NX)=cpw*VOLWH1(L,NY,NX)+cpi*VOLIH1(L,NY,NX)
!     IF(I.GT.331)THEN
!     WRITE(*,3376)'VOLWI',I,J,NX,NY,L,VOLW1(L,NY,NX)
!    2,VOLI1(L,NY,NX),VOLP1(L,NY,NX),VOLA1(L,NY,NX)
!    3,VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),VOLPH1(L,NY,NX)
!    3,VOLAH1(L,NY,NX),VOLT(L,NY,NX),VOLY(L,NY,NX)
!    4,THETWX(L,NY,NX),THETIX(L,NY,NX),THETPX(L,NY,NX)
!3376  FORMAT(A8,5I4,40E14.6)
!     ENDIF
!
!   MACROPOROSITY
!
!   FMAC,FGRD=macropore,micropore volume fractions
!   CNDH*=macropore hydraulic conductivity
!   TKS,TK1=soil temperature
!   FLU,HWFLU=subsurface water,convective heat fluxes
!   AREAU,AREAD=fractions of layer below natural,artifl water table
!
    IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      FMAC(L,NY,NX)=FHOL(L,NY,NX)*VOLAH1(L,NY,NX)/VOLAH(L,NY,NX)
      CNDH1(L,NY,NX)=CNDH(L,NY,NX)*(VOLAH1(L,NY,NX)/VOLAH(L,NY,NX))**2
    ELSE
      FMAC(L,NY,NX)=0.0_r8
      CNDH1(L,NY,NX)=0.0_r8
    ENDIF
    FGRD(L,NY,NX)=1.0-FMAC(L,NY,NX)
    TK1(L,NY,NX)=TKS(L,NY,NX)
    if(TK1(L,NY,NX)>1.e3_r8)then
      write(*,*)'TKS(L,NY,NX)',L,TKS(L,NY,NX)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
    IF(L.EQ.LWDPTH)THEN
      FLU(L,NY,NX)=PRECU(NY,NX)
      HWFLU(L,NY,NX)=cpw*TKA(NY,NX)*PRECU(NY,NX)
      FLU1(L,NY,NX)=FLU(L,NY,NX)*XNPH
      HWFLU1(L,NY,NX)=HWFLU(L,NY,NX)*XNPH
    ELSE
      FLU(L,NY,NX)=0.0_r8
      HWFLU(L,NY,NX)=0.0_r8
      FLU1(L,NY,NX)=0.0_r8
      HWFLU1(L,NY,NX)=0.0_r8
    ENDIF
    IF(CDPTH(L,NY,NX).GE.DTBLX(NY,NX))THEN
      AREAU(L,NY,NX)=AMIN1(1.0,AMAX1(0.0_r8,safe_adb(CDPTH(L,NY,NX)-DTBLX(NY,NX) &
        ,DLYR(3,L,NY,NX))))
    ELSE
      AREAU(L,NY,NX)=0.0_r8
    ENDIF
    IF(CDPTH(L,NY,NX).GE.DTBLY(NY,NX))THEN
      AREAUD(L,NY,NX)=AMIN1(1.0,AMAX1(0.0_r8,safe_adb(CDPTH(L,NY,NX)-DTBLY(NY,NX) &
        ,DLYR(3,L,NY,NX))))
    ELSE
      AREAUD(L,NY,NX)=0.0_r8
    ENDIF
30  CONTINUE
!
! ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
!     FOR USE AT INTERNAL TIME STEP IN SURFACE LITTER
!
!     THRMG=longwave emission from litter surface
!     VHCP1=volumetric heat capacity of litter
!     VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of litter
!     VOLWRX=maximum water retention by litter
!     XVOLT,XVOLW=free surface water+ice,water
!     VHCPRX=min heat capacity for litter water,heat fluxes
!     VOLR=litter volume
!     THETW*,THETI*,THETP*=water,ice,air concentrations
!     PSISM*=litter matric water potential
!
  THRMG(NY,NX)=0.0_r8
  VHCP1(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW(0,NY,NX) &
    +cpi*VOLI(0,NY,NX)
  VOLA1(0,NY,NX)=VOLA(0,NY,NX)
  VOLW1(0,NY,NX)=AMAX1(0.0,VOLW(0,NY,NX))
  VOLI1(0,NY,NX)=AMAX1(0.0,VOLI(0,NY,NX))
  VOLP1(0,NY,NX)=AMAX1(0.0,VOLA1(0,NY,NX)-VOLW1(0,NY,NX)-VOLI1(0,NY,NX))
  VOLWM(1,0,NY,NX)=VOLW1(0,NY,NX)
  VOLPM(1,0,NY,NX)=VOLP1(0,NY,NX)
  TVOLWI=VOLW1(0,NY,NX)+VOLI1(0,NY,NX)
  XVOLT(NY,NX)=AMAX1(0.0,TVOLWI-VOLWRX(NY,NX))
  IF(TVOLWI.GT.ZEROS(NY,NX))THEN
    VOLWRZ=VOLW1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
    VOLIRZ=VOLI1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
    XVOLW(NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)-VOLWRZ)
    XVOLI(NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)-VOLIRZ)
  ELSE
    XVOLW(NY,NX)=0.0_r8
    XVOLI(NY,NX)=0.0_r8
  ENDIF
  XVOLTM(1,NY,NX)=XVOLT(NY,NX)
  XVOLWM(1,NY,NX)=XVOLW(NY,NX)
  XVOLIM(1,NY,NX)=XVOLI(NY,NX)
  IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
    THETWX(0,NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)/VOLR(NY,NX))
    THETIX(0,NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)/VOLR(NY,NX))
    THETPX(0,NY,NX)=AMAX1(0.0,VOLP1(0,NY,NX)/VOLR(NY,NX)) &
      *AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
  ELSE
    THETWX(0,NY,NX)=0.0_r8
    THETIX(0,NY,NX)=0.0_r8
    THETPX(0,NY,NX)=1.0
  ENDIF
  THETPM(1,0,NY,NX)=THETPX(0,NY,NX)
  PSISM1(0,NY,NX)=PSISM(0,NY,NX)
  TK1(0,NY,NX)=TKS(0,NY,NX)
!  if(curday>=377)then
!    write(*,*)'LocalCopyStateVars'
!    write(*,*)'TK1(0,NY,NX)=',TK1(0,NY,NX)
!  endif
!     WRITE(*,7751)'THETPX',I,J,NX,NY
!    2,VOLW1(0,NY,NX),VOLI1(0,NY,NX),VOLP1(0,NY,NX),VOLA1(0,NY,NX)
!    3,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPX(0,NY,NX)
!    3,XVOLW(NY,NX),XVOLT(NY,NX),VOLWD(NY,NX),VOLWG(NY,NX),VOLR(NY,NX)
!    4,VOLW1(1,NY,NX),VOLI1(1,NY,NX),VOLP1(1,NY,NX)
!    5,VOLWRX(NY,NX)
!7751  FORMAT(A8,4I4,20E12.4)
  end subroutine LocalCopyStateVars
!------------------------------------------------------------------------------------------

  subroutine PrepWaterEnergyBalance(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L
!     begin_execution

  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS
      NUM(NY,NX)=NU(NY,NX)
!
!     ADJUST SURFACE ELEVATION USED IN RUNOFF FOR FREEZE-THAW, EROSION
!     AND SOC
!
!     ALTG,ALT=current,initial elevation of ground surface
!     CDPTH(NUM(NY,NX)-1,=depth of ground surface
!     ENGYP=cumulative rainfall energy impact on soil surface
!
      ALTG(NY,NX)=ALT(NY,NX)-CDPTH(NUM(NY,NX)-1,NY,NX)
      ENGYP(NY,NX)=ENGYP(NY,NX)*(1.0_r8-FENGYP)

      call LocalCopyStateVars(I,NY,NX)
!
!     SNOW AND RESIDUE COVERAGE OF SOIL SURFACE
!
!     FSNW,FSNX=fractions of snow,snow-free cover
!     DPTHS=snowpack depth
!     DPTHSX=minimum snowpack depth for full cover
!     BARE,CVRD=fractions of soil,litter cover
!     PRECA=precipitation+irrigation
!     PRECD,PRECB=direct,indirect precipn+irrign at soil surface
!
      FSNW(NY,NX)=AMIN1(1.0_r8,SQRT((DPTHS(NY,NX)/DPTHSX)))
      FSNX(NY,NX)=1.0_r8-FSNW(NY,NX)
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
        BARE(NY,NX)=AMIN1(1.0_r8,AMAX1(0.0_r8, &
          EXP(-0.8E-02_r8*(ORGC(0,NY,NX)/AREA(3,0,NY,NX)))))
      ELSE
        BARE(NY,NX)=1.0_r8
      ENDIF
      CVRD(NY,NX)=1.0_r8-BARE(NY,NX)
      PRECM(NY,NX)=1.0E+03_r8*PRECA(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      PRECD(NY,NX)=1.0E+03_r8*(PRECA(NY,NX)-TFLWCI(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
      PRECB(NY,NX)=1.0E+03_r8*(TFLWCI(NY,NX)-TFLWC(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
!     IF(PRECA(NY,NX).GT.0.0)THEN
!     WRITE(*,3112)'BARE',I,J,NX,NY,BARE(NY,NX)
!    2,FSNX(NY,NX),ORGC(0,NY,NX)/AREA(3,0,NY,NX),VOLWRX(NY,NX)
!    3,XVOLW(NY,NX),VOLWD(NY,NX)
!    4,PRECA(NY,NX),TFLWCI(NY,NX),TFLWC(NY,NX)
!    5,PRECA(NY,NX)*XNPH*1000*BARE(NY,NX),PRECD(NY,NX),PRECB(NY,NX)
!3112  FORMAT(A8,4I4,20E12.4)
!     ENDIF
!
!     RESIDUE WATER ABSORPTION CAPACITY
!
!     HCNDR=litter saturated hydraulic conductivity
!     DLYRR=litter depth
!
      HCNDR(NY,NX)=HCNDRR
      DLYRR(NY,NX)=AMAX1(2.5E-03_r8,DLYR(3,0,NY,NX))
!
!     DISTRIBUTION OF PRECIPITATION AND ITS HEAT AMONG SURFACE
!     RESIDUE, SOIL SURFACE, AND MACROPORES
!
!     PRECA,PRECW=rainfall+irrigation,snowfall (water equiv)
!     FLWQW=rainfall to snowpack
!     FLWSW=snowfall to snowpack
!     HFLWSW=convective heat flux to snowpack
!     FLWQB=precip to litter+soil surfaces
!     FLWQAX,FLWQBX=precip to soil,litter surfaces
!     HFLWQA,HFLWQB=convective heat flux to soil,litter surfaces
!     FLWQAS,FLWQAH=precip to soil micropores,macropores
!     TFLWC=canopy intercepted precipitation
      IF(PRECA(NY,NX).GT.0.0_r8.OR.PRECW(NY,NX).GT.0.0_r8)THEN
        FLWQW=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNW(NY,NX)
        FLWSW=PRECW(NY,NX)                                !snowfall
        HFLWSW=cps*TKA(NY,NX)*FLWSW+cpw*TKA(NY,NX)*FLWQW  !incoming heat flux from precipitations to snow-covered surface
        FLWQB=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNX(NY,NX)     !incoming precipitation to snow-free surface
        FLWQBX=FLWQB*CVRD(NY,NX)                          !water flux to snow-free coverd by litter
        HFLWQB=cpw*TKA(NY,NX)*FLWQBX                      !heat flux to snow-free surface covered by litter
        FLWQAX=FLWQB*BARE(NY,NX)                          !heat flux to snow-free surface not covered by litter
        HFLWQA=cpw*TKA(NY,NX)*FLWQAX
        FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)              !water flux to micropore
        FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)              !water flux to macropore
      ELSE
        FLWQW=-TFLWC(NY,NX)*FSNW(NY,NX)                   !
        FLWSW=0.0_r8
        HFLWSW=cpw*TKA(NY,NX)*FLWQW
        FLWQB=-TFLWC(NY,NX)*FSNX(NY,NX)
        FLWQBX=FLWQB*CVRD(NY,NX)
        HFLWQB=cpw*TKA(NY,NX)*FLWQBX
        FLWQAX=FLWQB*BARE(NY,NX)
        HFLWQA=cpw*TKA(NY,NX)*FLWQAX
        FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)
        FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)
      ENDIF
!
!     PRECIP ON SNOW ARRAYS EXPORTED TO TRNSFR.F, TRNSFRS.F
!     FOR SOLUTE FLUX CALCULATIONS
!
!     PRECW,PRECR,PRECQ,PRECI=snow,rain,snow+rain,irrigation
!     VHCPW,VHCPWX=current, minimum snowpack heat capacities
!     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!
      IF(PRECW(NY,NX).GT.0.0.OR.(PRECR(NY,NX).GT.0.0 &
        .AND.VHCPW(1,NY,NX).GT.VHCPWX(NY,NX)))THEN
        FLQRQ(NY,NX)=0.0_r8
        FLQRI(NY,NX)=0.0_r8
        FLQGQ(NY,NX)=PRECQ(NY,NX)
        FLQGI(NY,NX)=PRECI(NY,NX)
      ELSEIF((PRECQ(NY,NX).GT.0.0.OR.PRECI(NY,NX).GT.0.0) &
        .AND.VHCPW(1,NY,NX).LE.VHCPWX(NY,NX))THEN
        FLQRQ(NY,NX)=FLWQBX*PRECQ(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
        FLQRI(NY,NX)=FLWQBX*PRECI(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
        FLQGQ(NY,NX)=PRECQ(NY,NX)-FLQRQ(NY,NX)
        FLQGI(NY,NX)=PRECI(NY,NX)-FLQRI(NY,NX)
      ELSE
        FLQRQ(NY,NX)=0.0_r8
        FLQRI(NY,NX)=0.0_r8
        FLQGQ(NY,NX)=0.0_r8
        FLQGI(NY,NX)=0.0_r8
      ENDIF
!
!     GATHER PRECIPITATION AND MELTWATER FLUXES AND THEIR HEATS
!     AMONG ATMOSPHERE, SNOWPACK, RESIDUE AND SOIL SURFACES
!     INTO LOCAL ARRAYS FOR USE IN MASS AND ENERGY EXCHANGE
!     ALGORITHMS
!
!     XNPH=internal time step for fluxes through soil profile
!
!     FLW0S,FLQ0I,FLQ0W=snow,ice,water input to snowpack
!     HWFLQ0=convective heat flux to snowpack
!     FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
!     HWFLQ1,HWFLY1=convective heat flux to soil,litter surfaces
!
      FLQ0S(NY,NX)=FLWSW*XNPH
      FLQ0I(NY,NX)=0.0_r8
      FLQ0W(NY,NX)=FLWQW*XNPH
      HWFLQ0(NY,NX)=HFLWSW*XNPH
      FLQ1(NY,NX)=FLWQAS*XNPH
      FLH1(NY,NX)=FLWQAH*XNPH
      FLY1(NY,NX)=FLWQBX*XNPH
      HWFLQ1(NY,NX)=HFLWQA*XNPH
      HWFLY1(NY,NX)=HFLWQB*XNPH
!     IF(I.EQ.118.AND.NX.EQ.3.AND.NY.EQ.4)THEN
!     WRITE(*,4422)'FLQ0W',I,J,FLQ0W(NY,NX),FLWQW,XNPH
!     WRITE(*,4422)'FLY',I,J,PRECA(NY,NX),TFLWC(NY,NX),FLY1(NY,NX)
!    2,PSISM1(0,NY,NX),PSISM(0,NY,NX)
!    2,FLQ1(NY,NX),FLH1(NY,NX),FLWQBX
!    2,FLWQAS,FLWQAH
!    3,FGRD(NUM(NY,NX),NY,NX),FMAC(NUM(NY,NX),NY,NX)
!    4,FHOL(L,NY,NX),VOLAH1(L,NY,NX),VOLAH1(L,NY,NX)
!    5,FLWQAX,PRECA(NY,NX),TFLWC(NY,NX),FLWQBX
!    6,BARE(NY,NX),ORGC(0,NY,NX),XVOLW(NY,NX),VOLWG(NY,NX)
!    7,VOLW1(0,NY,NX),VOLWRX(NY,NX)
!4422  FORMAT(A8,2I4,60F18.6)
!     ENDIF
!
!     INITIALIZE PARAMETERS, FLUXES FOR ENERGY EXCHANGE
!     AT SNOW, RESIDUE AND SOIL SURFACES
!
!     RADGX=shortwave radiation at ground surface
!     RADXW,RADXG,RADXR= shortwave radn at snowpack,soil,litter
!     FRADG=fraction of shortwave radiation at ground surface
!     FSNW,FSNX=fractions of snow,snow-free cover
!     BARE,CVRD=fractions of soil,litter cover
!     XNPS=internal time step for fluxes through snowpack
!     THRYX=longwave radiation at ground surface
!     THRYW,THRYG,THRYR=longwave radn incident at snowpack,soil,litter
!     THRMW,THRMS,THRMR=longwave radn emitted by snowpack,soil,litter
!     EMMW,EMMS,EMMR=emissivity of snowpack,soil,litter surfaces
!
      RADGX=RADG(NY,NX)*XNPH
      RADXW(NY,NX)=RADGX*FSNW(NY,NX)*XNPS
      RADXG(NY,NX)=RADGX*FSNX(NY,NX)*BARE(NY,NX)
      RADXR(NY,NX)=RADGX*FSNX(NY,NX)*CVRD(NY,NX)*XNPR
      THRYX=(THS(NY,NX)*FRADG(NY,NX)+THRMCX(NY,NX))*XNPH
      THRYW(NY,NX)=THRYX*FSNW(NY,NX)*XNPS
      THRYG(NY,NX)=THRYX*FSNX(NY,NX)*BARE(NY,NX)
      THRYR(NY,NX)=THRYX*FSNX(NY,NX)*CVRD(NY,NX)*XNPR
      THRMW(NY,NX)=EMMW*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPY
      THRMS(NY,NX)=EMMS*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*BARE(NY,NX)*XNPH
      THRMR(NY,NX)=EMMR*2.04E-10_r8*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*CVRD(NY,NX)*XNPZ
!
!     AERODYNAMIC RESISTANCE OF CANOPY TO SNOW/RESIDUE/SOIL
!     SURFACE ENERGY EXCHANGE WITH ATMOSPHERE
!
!     ALFZ=parameter for canopy effect on windspeed
!     FRADG=fraction of shortwave radiation at ground surface
!     RAB,RAC=isothermal blr above canopy, canopy blr
!     ZT,ZS=canopy, surface roughness heights
!     UA,UAG=windspeeds above,below canopy
!     VPQ,VPA=vapor pressure within,above canopy
!     TKQ,TKA=temperature within,above canopy
!     TLEX,TSHX=net latent,sensible heat fluxes x blrs from prev hour
!     VAP=latent heat of evaporation
!     1.25E-03=heat capacity of air
!     AREA=surface area of grid cell
!
      ALFZ=2.0_r8*(1.0_r8-FRADG(NY,NX))
      IF(RAB(NY,NX).GT.ZERO.AND.ZT(NY,NX).GT.ZS(NY,NX) &
        .AND.ALFZ.GT.ZERO)THEN
        RAC(NY,NX)=AMIN1(RACX,AMAX1(0.0_r8,ZT(NY,NX)*EXP(ALFZ) &
          /(ALFZ/RAB(NY,NX))*AMAX1(0.0_r8,EXP(-ALFZ*ZS(NY,NX)/ZT(NY,NX)) &
          -EXP(-ALFZ*(ZD(NY,NX)+ZR(NY,NX))/ZT(NY,NX)))))
        UAG=UA(NY,NX)*EXP(-ALFZ)
      ELSE
        RAC(NY,NX)=0.0_r8
        UAG=UA(NY,NX)
      ENDIF
      VPQ(NY,NX)=VPA(NY,NX)-TLEX(NY,NX)/(VAP*AREA(3,NUM(NY,NX),NY,NX))
      TKQ(NY,NX)=TKA(NY,NX)-TSHX(NY,NX)/(1.25E-03_r8*AREA(3,NUM(NY,NX),NY,NX))
!     WRITE(*,3114)'RAC',I,J,NX,NY,RAC(NY,NX),FRADG(NY,NX)
!    2,RAB(NY,NX),ZT(NY,NX),ZS(NY,NX),ALFZ
!    3,VPQ(NY,NX),TKQ(NY,NX)
!3114  FORMAT(A8,4I4,20E12.4)
!
!     AERODYNAMIC RESISTANCE OF SNOWPACK, RESIDUE AND SOIL
!     SURFACES TO ENERGY EXCHANGE WITH ATMOSPHERE
!     Soil Sci. Soc. Am. J. 48:25-32
!
!     RAR=porosity-unlimited litter blr
!     DLYRR=litter depth
!     WGSGR=vapor diffusivity in litter
!     RAG,RAGW,RAGR=isothermal blrs at ground,snowpack,litter surfaces
!     RARX=boundary layer resistance (blr) of litter surface
!     THETPX*=air-filled porosity of litter
!     DFVR=porosity limitation to diffusion through litter
!     POROQ=litter tortuosity
!     RAR1=porosity-limited litter blr
!     PAREX,PARSX=conductances for latent,sensible heat fluxes
!     PAREW,PARSW=conductances for snowpack latent,sensible heatfluxes
!     PAREG,PARSG=conductances for soil latent,sensible heat fluxes
!     PARER,PARSR=conductances for litter latent,sensible heat fluxes
!     XNPR=internal time step for fluxes through litter
!
      RAR(NY,NX)=DLYRR(NY,NX)/WGSGR(NY,NX)
      RAG(NY,NX)=RAC(NY,NX)+RAB(NY,NX)
      RAGW(NY,NX)=RAG(NY,NX)
      RAGR(NY,NX)=RAG(NY,NX)+RARX
      RARG(NY,NX)=RAGR(NY,NX)
      THETPX0=AMAX1(ZERO2,THETPX(0,NY,NX))
      DFVR=THETPX0*POROQ*THETPX0/POROS(0,NY,NX)
      RAR1=RAG(NY,NX)+RAR(NY,NX)/DFVR
      PAREX=AREA(3,NUM(NY,NX),NY,NX)*XNPH               !conductance for latent heat flux
      PARSX=1.25E-03_r8*AREA(3,NUM(NY,NX),NY,NX)*XNPH   !conductance for sensible heat flux
      PAREW(NY,NX)=PAREX*FSNW(NY,NX)*XNPS
      PARSW(NY,NX)=PARSX*FSNW(NY,NX)*XNPS
      PAREG(NY,NX)=PAREX*FSNX(NY,NX)
      PARER(NY,NX)=PAREX*FSNX(NY,NX)*XNPR*CVRD(NY,NX)
      PARSG(NY,NX)=PARSX*FSNX(NY,NX)
      PARSR(NY,NX)=PARSX*FSNX(NY,NX)*XNPR*CVRD(NY,NX)
!     WRITE(*,3115)'RAR',I,J,NX,NY,RAR1,RAR(NY,NX),DFVR
!    2,THETPX(0,NY,NX),THETPX0,POROS(0,NY,NX),DLYRR(NY,NX)
!    3,RAG(NY,NX),PARER(NY,NX),PAREX,FSNX(NY,NX),XNPR,CVRD(NY,NX)
!    4,BARE(NY,NX)
!3115  FORMAT(A8,4I4,30E12.4)
!
!     BOUNDARY LAYER CONDUCTANCES FOR EXPORT TO TRNSFR.F
!
!     RAS,RASL=blrs of snowpack,snowpack layer
!     VOLS,VOLS1=volume of snowpack,snowpack layer
!     DLYRS=snowpack later depth
!     WGSGW=vapor diffusivity in snowpack
!     THETPL=snowpack air-filled porosity
!     THETPI=air content of ice
!     VOLS0,VOLI0,VOLW0,VOLS1=snow,ice,water,total volumes of snowpack
!     PARR=boundary layer conductance above litter,soil surfaces
!
      RAS=0.0_r8
      IF(VOLS(NY,NX).GT.ZEROS2(NY,NX))THEN
        DO 9775 L=1,JS
          IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            RASX=DLYRS(L,NY,NX)/WGSGW(L,NY,NX)
            THETPL=AMAX1(THETPI,1.0_r8-(VOLS0(L,NY,NX)+VOLI0(L,NY,NX) &
              +VOLW0(L,NY,NX))/VOLS1(L,NY,NX))
            RASL=RASX/AMAX1(ZERO,THETPL)**2.0_r8
            RAS=RAS+RASL
!     WRITE(*,3113)'RAS',I,J,NX,NY,L,ALFZ,RAS,RASL,RASX
!    2,DLYRS(L,NY,NX),WGSGW(L,NY,NX),THETPL,THETPI,VOLS0(L,NY,NX)
!    3,VOLI0(L,NY,NX),VOLW0(L,NY,NX),VOLS1(L,NY,NX),TKW(L,NY,NX)
!3113  FORMAT(A8,5I4,40E12.4)
          ENDIF
9775    CONTINUE
      ENDIF
      PARR(NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*XNPH/(RAGR(NY,NX)+RAS)   !this includes snow layer resistance
!     IF(NX.EQ.1)THEN
!     WRITE(*,3111)'RAC',I,J,NX,NY,ALFZ,RAC(NY,NX)
!    2,PARR(NY,NX),ZT(NY,NX),RAB(NY,NX)
!    2,RAGR(NY,NX),RAS,VOLS(NY,NX)
!    3,DLYRR(NY,NX),RAG(NY,NX),RAGR(NY,NX)
!    4,THETX,THETPX(0,NY,NX),VHCP1(0,NY,NX)
!    5,WGSGR(NY,NX),VOLW1(0,NY,NX)
!    5,VOLI1(0,NY,NX),VOLP1(0,NY,NX),VOLR(NY,NX),VOLA1(0,NY,NX)
!    4,TLEX(NY,NX),TSHX(NY,NX),RADG(NY,NX),THS(NY,NX)
!    5,FRADG(NY,NX),THRMCX(NY,NX),ZS(NY,NX)
!    6,XVOLW(NY,NX),VHCPRX(NY,NX)/cpw,VOLWD(NY,NX),ORGC(0,NY,NX)
!3111  FORMAT(A8,4I4,40E12.4)
!     ENDIF
9990  CONTINUE
9995  CONTINUE

  call InitSoilHydrauics(NHW,NHE,NVN,NVS)

  end subroutine PrepWaterEnergyBalance
!------------------------------------------------------------------------------------------

  subroutine InitSoilHydrauics(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: N,N1,N2,N3,N4,N5,N6
  integer :: NY,NX,L
!     begin_execution
!
!     INITIALIZE SOIL HYDRAULIC PARAMETERS IN LOCAL ARRAYS
!     FOR LATER USE IN WATER TRANSFER ALGORITHMS
!
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  DO 9985 NX=NHW,NHE
    DO 9980 NY=NVN,NVS
      DO 35 L=NUM(NY,NX),NL(NY,NX)
        DO 40 N=NCN(NY,NX),3
          N1=NX
          N2=NY
          N3=L
          IF(N.EQ.1)THEN
            IF(NX.EQ.NHE)THEN
              GO TO 50
            ELSE
              N4=NX+1
              N5=NY
              N6=L
            ENDIF
          ELSEIF(N.EQ.2)THEN
            IF(NY.EQ.NVS)THEN
              GO TO 50
            ELSE
              N4=NX
              N5=NY+1
              N6=L
            ENDIF
          ELSEIF(N.EQ.3)THEN
            IF(L.EQ.NL(NY,NX))THEN
              GO TO 50
            ELSE
              N4=NX
              N5=NY
              N6=L+1
            ENDIF
          ENDIF
!
    !     MACROPORE CONDUCTIVITY FROM 'HOUR1' AND GRAVITATIONAL
    !     GRADIENT USED TO CALCULATE MACROPORE FLOW FOR USE BELOW
    !
    !     CNDH1=macropore hydraulic conductivity
    !     AVCNHL=macropore hydraulic conductance
    !     DLYR=layer depth
    !
          IF(CNDH1(N3,N2,N1).GT.ZERO.AND.CNDH1(N6,N5,N4) &
            .GT.ZERO)THEN
            AVCNHL(N,N6,N5,N4)=2.0*CNDH1(N3,N2,N1)*CNDH1(N6,N5,N4) &
              /(CNDH1(N3,N2,N1)*DLYR(N,N6,N5,N4)+CNDH1(N6,N5,N4) &
              *DLYR(N,N3,N2,N1))
          ELSE
            AVCNHL(N,N6,N5,N4)=0.0_r8
          ENDIF
50      CONTINUE
40    CONTINUE
35  CONTINUE
9980  CONTINUE
9985  CONTINUE
  end subroutine InitSoilHydrauics
!------------------------------------------------------------------------------------------

  subroutine SnowSurfaceResidualIteration(L,M,NY,NX)
  implicit none
  integer, intent(in) :: L,M,NY,NX

  integer :: NN
! begin_execution
  DO 4000 NN=1,NPR
    !
    ! VAPOR FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! VP0,VPR,VPY=snowpack,litter, equilibrium vapor concentration
    ! TK0X,TKXR=snowpack,litter temperature
    ! PSISM1=litter matric water potential
    ! FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
    ! AREA=area of grid cell
    ! FSNW,CVRD=snow,litter cover fraction
    ! XNPQ=time step for flux calculation
    ! FLVRSX=snow-litter vapor flux
    ! HFLVRSX=convective heat flux from snow-litter vapor flux
    !
    !VPR=2.173E-03_r8/TKXR*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TKXR)) &
    !  *EXP(18.0_r8*PSISM1(0,NY,NX)/(8.3143_r8*TKXR))    !in residue vapor pressure
    VPR=vapsat(TKXR)*EXP(18.0_r8*PSISM1(0,NY,NX)/(8.3143_r8*TKXR))
    if(abs(VPR)>1.e20_r8)then
      write(*,*)'TKXR=',TKXR,TK1(0,NY,NX),TK1(NUM(NY,NX),NY,NX)
      write(*,*)'PSISM1(0,NY,NX)=',PSISM1(0,NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    IF(VOLP01.GT.ZEROS2(NY,NX).AND.THETPM(M,0,NY,NX).GT.THETX)THEN
      !VP0=2.173E-03_r8/TK0X*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TK0X))   !in snow vapor pressure
      VP0=vapsat(TK0X)
      FLVC=ATCNVR*(VP0-VPR)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ !snow <-> residue vapor flux
      VPY=(VP0*VOLP01+VPR*VOLPM(M,0,NY,NX))/(VOLP01+VOLPM(M,0,NY,NX))           !volume weighted vapor pressure
      FLVX=(VP0-VPY)*VOLP01*XNPC
      IF(FLVC.GE.0.0_r8)THEN
        FLVSRX=AMAX1(0.0_r8,AMIN1(FLVC,FLVX))
        HFLVSRX=(cpw*TK0X+VAP)*FLVSRX
      ELSE
        FLVSRX=AMIN1(0.0_r8,AMAX1(FLVC,FLVX))
        HFLVSRX=(cpw*TKXR+VAP)*FLVSRX
      ENDIF
    ELSE
      FLVSRX=0.0_r8
      HFLVSRX=0.0_r8
    ENDIF
    !
    ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! TKY=snow-litter equilibrium temperature
    ! HFLWC,HFLWX=snow-litter heat flux unltd,ltd by heat
    ! HFLWSRX=snow-litter heat flux
    ! VHCPWMM= volumetric heat capacity in snow layer
    TKY=(TK0X*VHCPWMM(L,NY,NX)+TKXR*VHCP1(0,NY,NX))/(VHCPWMM(L,NY,NX)+VHCP1(0,NY,NX))
    HFLWX=(TK0X-TKY)*VHCPWMM(L,NY,NX)*XNPC
    HFLWC=ATCNDR*(TK0X-TKXR)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
    IF(HFLWC.GE.0.0)THEN
      HFLWSRX=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
    ELSE
      HFLWSRX=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
    ENDIF
!     WRITE(*,7752)'TKXR',I,J,M,MM,NX,NY,L
!    2,FLVC,FLVX,VP1,VPR,VPY,ATCNVS,FSNW(NY,NX),BARE(NY,NX)
!    3,VOLP01,VOLPM(M,NUM(NY,NX),NY,NX),TK1X,TKXR
!    4,HFLVR1X,VHCP1(NUM(NY,NX),NY,NX),VHCP1(0,NY,NX)
!    3,VOLPM(M,0,NY,NX),VOLPM(M,NUM(NY,NX),NY,NX),VOLP01
!    2,THETPM(M,0,NY,NX),TKXR
!
!     VAPOR FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
!
!     THETPM,VOLPM=air-filled porosity,volume
!     VP1,VPY=soil,litter-soil equilibrium vapor concentration
!     TK1X=soil temperature
!     PSISV1=soil matric+osmotic water potentials
!     FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
!     FLVR1X=litter-soil vapor flux
!     HFLVR1X=convective heat of litter-soil vapor flux
!     TKXR,TK1X=interim calculation of litter,soil temperatures
!
    IF(VOLPM(M,0,NY,NX).GT.ZEROS(NY,NX) &
      .AND.VOLPM(M,NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      !VP1=2.173E-03/TK1X*0.61*EXP(5360.0*(3.661E-03-1.0/TK1X)) &
      !  *EXP(18.0*PSISV1/(8.3143*TK1X))
      VP1=vapsat(TK1X)*EXP(18.0*PSISV1/(8.3143*TK1X))
      FLVC=ATCNVS*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
      if(abs(FLVC)>1.e20_r8)then
        write(*,*)'ATCNVS=',ATCNVS,VPR,VP1
        write(*,*)'FSNW(NY,NX)*CVRD(NY,NX)=',FSNW(NY,NX),CVRD(NY,NX)
        write(*,*)'at line',__LINE__
        call endrun(trim(mod_filename)//'at line',__LINE__)
      endif
      VPY=(VPR*VOLPM(M,0,NY,NX)+VP1*VOLPM(M,NUM(NY,NX),NY,NX)) &
        /(VOLPM(M,0,NY,NX)+VOLPM(M,NUM(NY,NX),NY,NX))
      FLVX=(VPR-VPY)*VOLPM(M,0,NY,NX)*XNPC
      IF(FLVC.GE.0.0)THEN
        FLVR1X=AMAX1(0.0,AMIN1(FLVC,FLVX,VOLW0M(L,NY,NX)*XNPB))
        if(abs(FLVR1X)>1.0e20_r8)then
          write(*,*)'FLVC,FLVX,VOLW0M(L,NY,NX)*XNPB=',FLVC,FLVX,VOLW0M(L,NY,NX)*XNPB
          write(*,*)'at line',__LINE__
          call endrun(trim(mod_filename)//'at line',__LINE__)
        endif
        HFLVR1X=(cpw*TKXR+VAP)*FLVR1X
      ELSE
        FLVR1X=AMIN1(0.0,AMAX1(FLVC,FLVX))
        if(abs(FLVR1X)>1.0e20_r8)then
          write(*,*)'FLVC,FLVX=',FLVC,FLVX
          write(*,*)'at line',__LINE__
          call endrun(trim(mod_filename)//'at line',__LINE__)
        endif
        HFLVR1X=(cpw*TK1X+VAP)*FLVR1X
      ENDIF
    ELSE
      FLVR1X=0.0_r8
      HFLVR1X=0.0_r8
    ENDIF
    TKXR=TKXR-HFLVR1X/VHCP1(0,NY,NX)
    TK1X=TK1X+HFLVR1X/VHCP1(NUM(NY,NX),NY,NX)
!     IF(NY.EQ.6)THEN
!     WRITE(*,7752)'TK1X',I,J,M,MM,NX,NY,NN
!    2,FLVC,FLVX,VP1,VPR,VPY,ATCNVS,FSNW(NY,NX),BARE(NY,NX)
!    3,VOLP01,VOLPM(M,NUM(NY,NX),NY,NX),TK1X,TKXR,TK1(0,NY,NX)
!    4,HFLVR1X,VHCP1(NUM(NY,NX),NY,NX),VHCP1(0,NY,NX)
!    3,VOLPM(M,0,NY,NX),VOLPM(M,NUM(NY,NX),NY,NX),VOLP01
!    2,THETPM(M,0,NY,NX)
!     ENDIF
!
!     HEAT FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
!
!     TKY=litter-soil equilibrium temperature
!     HFLWC,HFLWX=litter-soil heat flux unltd,ltg by heat
!     HFLWR1X=litter-soil heat flux
!
    TKY=(TKXR*VHCP1(0,NY,NX)+TK1X*VHCP1(NUM(NY,NX),NY,NX)) &
      /(VHCP1(0,NY,NX)+VHCP1(NUM(NY,NX),NY,NX))
    HFLWX=(TKXR-TKY)*VHCP1(0,NY,NX)*XNPC
    HFLWC=ATCNDS*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
    IF(HFLWC.GE.0.0)THEN
      HFLWR1X=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
    ELSE
      HFLWR1X=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
    ENDIF
!
!     ACCUMULATE SNOW-LITTER, LITTER-SOIL HEAT FLUXES
!     WITHIN LONGER TIME STEP FOR SNOWPACK FLUX CALCULATIONS
!
    FLVSR=FLVSR+FLVSRX
    HFLVSR=HFLVSR+HFLVSRX
    HFLWSR=HFLWSR+HFLWSRX
    FLVR1=FLVR1+FLVR1X
    if(abs(FLVR1)>1.0e20_r8)then
      write(*,*)'FLVR1X=',FLVR1X
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HFLVR1=HFLVR1+HFLVR1X
    HFLWR1=HFLWR1+HFLWR1X
    TK0X=TK0X-HFLVSRX/VHCPWMM(L,NY,NX)
    TKXR=TKXR+(HFLVSRX-HFLWR1X)/VHCP1(0,NY,NX)
    TK1X=TK1X+HFLWR1X/VHCP1(NUM(NY,NX),NY,NX)
!     IF(I.EQ.53)THEN
!     WRITE(*,1114)'FLVR0',I,J,M,N,NX,NY
!    2,FLVSR,FLVSRX,FLVC,FLVX,ATCNVR,VP0,VPR,VPY,FSNW(NY,NX),CVRD(NY,NX)
!    3,XNPQ,THETPM(M,0,NY,NX),VOLP01,XNPX,XNPV
!    2,TK0M(1,NY,NX),TK1(0,NY,NX)
!    2,TK1(NUM(NY,NX),NY,NX),TK0X,TKXR,TK1X,FLVR1,HWFLVR1,FLVS1
!    4,HWFLVS1,HFLC0R1,HFLCR11,FLVR,HWFLVR,FLVS,HWFLVS
!    3,HFLC0R,HFLCR1,VPQ(NY,NX),VP0,VPR,VP1,PSISM1(0,NY,NX),PSISV1
!    5,AVCNVR,ATCNDR,AVCNVS,ATCNDS,VHCPWMM(L,NY,NX),VHCP1(0,NY,NX)
!    6,VHCP1(NUM(NY,NX),NY,NX),DLYRR(NY,NX),DLYRS0(NY,NX),CNV01,CNVR1
!    7,CNV11,CNV1,THETPX(NUM(NY,NX),NY,NX),POROQ
!    2,WGSGL(NUM(NY,NX),NY,NX),CVRD(NY,NX),HFLXR,HFLVS1,HFLCR1
!1114  FORMAT(A8,6I4,60E12.4)
!     ENDIF
4000  CONTINUE
  end subroutine SnowSurfaceResidualIteration
!------------------------------------------------------------------------------------------

  subroutine SnowPackIteration(M,NY,NX,ICHKL)
  implicit none
  integer, intent(in) :: M,NY,NX
  integer, intent(out) :: ICHKL
  integer :: L,L2
  ! begin_execution
  ! PHYSICAL AND HYDRAULIC PROPERTIES OF SNOWPACK INCLUDING
  ! AIR AND WATER-FILLED POROSITY, WATER POTENTIAL OF UNDERLYING
  ! SOIL SURFACE USED IN FLUX CALCULATIONS
  !
  ! VHCPW,VHCPWX=current, minimum snowpack heat capacities
  ! VOLS0M,VOLI0M,VOLW0M,VOLS1=snow,ice,water,total snowpack volume
  ! DENSS,DENSI,DENS0=snow,ice,minimum snow density
  ! AREA=area of grid cell
  ! VOLP01=snowpack air volume
  ! THETP1=snowpack air concentration
  ! CNV1=snowpack vapor conductivity
  ! VP1=snowpack vapor concentration
  ! TK0M=snowpack temperature
  ! WGSGW=vapor diffusivity
  ! DENSW1=snowpack density
  ICHKL=0
  DO 9880 L=1,JS
    IF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
      VOLS1(L,NY,NX)=VOLS0M(L,NY,NX)/DENSS(L,NY,NX) &
        +VOLW0M(L,NY,NX)+VOLI0M(L,NY,NX)
      DLYRS0(L,NY,NX)=VOLS1(L,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
      VOLP01=AMAX1(0.0_r8,VOLS1(L,NY,NX)-VOLS0M(L,NY,NX)-VOLI0M(L,NY,NX)-VOLW0M(L,NY,NX))
      THETP1=AMAX1(THETPI,VOLP01/VOLS1(L,NY,NX))
      CNV1=THETP1**2.0*WGSGW(L,NY,NX)
      VP1=vapsat(TK0M(L,NY,NX))
      !VP1=2.173E-03_r8/TK0M(L,NY,NX)*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TK0M(L,NY,NX)))
      IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        DENSW1=AMIN1(0.6_r8,(VOLS0M(L,NY,NX)+VOLW0M(L,NY,NX) &
          +VOLI0M(L,NY,NX)*DENSI)/VOLS1(L,NY,NX))
      ELSE
        DENSW1=DENS0(NY,NX)
      ENDIF
      !
      ! SNOW THERMAL CONDUCTIVITY FROM J GLACIOL 43:26-41
      !
      ! TCND1W=snow thermal conductivity
      ! FLWQX=porosity-unconstrained snow water flux
      !
      TCND1W=0.0036_r8*10**(2.650*DENSW1-1.652)
      !
      ! DISCHARGE OF MELTWATER AND ITS HEAT FROM SNOWPACK LAYER
      ! TO LOWER SNOWPACK LAYER
      !
      FLWQX=AMAX1(0.0_r8,AMAX1(0.0_r8,VOLW0M(L,NY,NX)) &
        -0.05_r8*AMAX1(0.0_r8,VOLS0M(L,NY,NX)))*XNPA
      !
      ! WATER AND HEAT FLUXES IN SNOWPACK
      !
      ! DLYRS0=snow layer thickness
      ! FLWQM=porosity-constrained snow water flux
      ! HFLWQM=convective heat flux from water flux
      !
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWMM(L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        !if L==JS-1, L2==JS, so top layer is treated here.
        VOLS1(L2,NY,NX)=VOLS0M(L2,NY,NX)/DENSS(L2,NY,NX) &
          +VOLW0M(L2,NY,NX)+VOLI0M(L2,NY,NX)
        DLYRS0(L2,NY,NX)=VOLS1(L2,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
        VOLP02=VOLS1(L2,NY,NX)-VOLS0M(L2,NY,NX)-VOLI0M(L2,NY,NX)-VOLW0M(L2,NY,NX)
        THETP2=AMAX1(THETPI,VOLP02/VOLS1(L2,NY,NX))
        FLWQM=AMIN1(THETP2,FLWQX)
        HFLWQM=cpw*TK0M(L,NY,NX)*FLWQM
        !
        ! VAPOR FLUX IN SNOWPACK
        !
        ! VOLP01,VOLP02=air-filled volumes of source, destination layers
        ! L2=destination layer
        ! CNV1,CNV2=vapor conductivities of source, destination layers
        ! VP1,VP2=vapor concentrations of source, destination layers
        ! TK0M=soil temperature
        ! ATCNVW=snow vapor conductance
        ! DLYRS0=snow layer thickness
        ! FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
        ! FLVSS,HFLVSS=vapor flux and its convective heat flux
        !
        IF(VOLP01.GT.ZEROS2(NY,NX).AND.VOLP02.GT.ZEROS2(NY,NX))THEN
          CNV2=THETP2**2.0*WGSGW(L2,NY,NX)
          !VP2=2.173E-03/TK0M(L2,NY,NX) &
          !  *0.61*EXP(5360.0*(3.661E-03-1.0/TK0M(L2,NY,NX)))
          VP2=vapsat(TK0M(L2,NY,NX))
          ATCNVW=2.0*CNV1*CNV2/(CNV1*DLYRS0(L2,NY,NX) &
            +CNV2*DLYRS0(L,NY,NX))
          FLVC=ATCNVW*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPY
          VPY=(VP1*VOLP01+VP2*VOLP02)/(VOLP01+VOLP02)
          FLVX=(VP1-VPY)*VOLP01*XNPA
          IF(FLVC.GE.0.0)THEN
            FLVSS=AMAX1(0.0,AMIN1(FLVC,FLVX,VOLW0M(L,NY,NX)*XNPX))
            HFLVSS=(cpw*TK0M(L,NY,NX)+VAP)*FLVSS
          ELSE
            FLVSS=AMIN1(0.0,AMAX1(FLVC,FLVX,-VOLW0M(L2,NY,NX)*XNPX))
            HFLVSS=(cpw*TK0M(L2,NY,NX)+VAP)*FLVSS
          ENDIF
        ELSE
          FLVSS=0.0_r8
          HFLVSS=0.0_r8
        ENDIF
        !
        ! HEAT FLUX IN SNOWPACK
        !
        ! DENSW2,TCNDW2=density,thermal conductivity in destination layer
        ! ATCNDW=thermal conductance
        ! DLYRS0=layer thickness
        ! TKY=equilibrium temperature
        ! HFLWX,HFLWC=heat-constrained,heat-unconstrained heat fluxes
        ! VHCPWMM,TK0M=volumetric heat capacity,temperature
        ! XNPX=time step for flux calculations
        ! FSNW=snow cover fraction
        ! XNPY=time step for snowpack flux calculations
        ! HFLWSS=snowpack heat flux
        ! FLW0S,FLQ0I,FLQ0W=snow,ice,water fluxes through snowpack
        ! HFLW0W=convective heat flux snow,water,ice fluxes
        !
        IF(VOLS1(L2,NY,NX).GT.ZEROS2(NY,NX))THEN
          DENSW2=AMIN1(0.6,(VOLS0M(L2,NY,NX)+VOLW0M(L2,NY,NX) &
            +VOLI0M(L2,NY,NX)*DENSI)/VOLS1(L2,NY,NX))
        ELSE
          DENSW2=DENS0(NY,NX)
        ENDIF
        TCND2W=0.0036_r8*10**(2.650_r8*DENSW2-1.652_r8)
        ATCNDW=2.0_r8*TCND1W*TCND2W/(TCND1W*DLYRS0(L2,NY,NX) &
          +TCND2W*DLYRS0(L,NY,NX))
        TKY=(TK0M(L,NY,NX)*VHCPWMM(L,NY,NX)+TK0M(L2,NY,NX) &
          *VHCPWMM(L2,NY,NX))/(VHCPWMM(L,NY,NX)+VHCPWMM(L2,NY,NX))
        HFLWX=(TK0M(L,NY,NX)-TKY)*VHCPWMM(L,NY,NX)*XNPA
        HFLWC=ATCNDW*(TK0M(L,NY,NX)-TK0M(L2,NY,NX)) &
          *AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPY
        IF(HFLWC.GE.0.0)THEN
          HFLWSS=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
        ELSE
          HFLWSS=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
        ENDIF
        FLW0T=FLWQM+FLVSS
        HFLW0T=HFLWQM+HFLVSS+HFLWSS
        FLW0S(L2,NY,NX)=0.0_r8
        FLW0W(L2,NY,NX)=FLW0T
        FLW0I(L2,NY,NX)=0.0_r8
        HFLW0W(L2,NY,NX)=HFLW0T
        FLQWM(M,L2,NY,NX)=FLQWM(M,L2,NY,NX)+FLWQM
        !     IF(NX.EQ.3.AND.NY.EQ.3.AND.L.EQ.1)THEN
        !     WRITE(*,7757)'FLW0',I,J,M,MM,L2,FLW0W(L2,NY,NX),FLW0T,FLWQM,FLVSS
        !    2,HFLW0W(L2,NY,NX),HFLW0T,HFLWQM,HFLVSS,HFLWSS,VP1,VP2,FLVX
        !    3,TCND1W,TCND2W,DENSW1,DENSW2,HFLXW2,VHCPWM2,THETP2,FLWQX
        !    2,VOLS0M(L,NY,NX),VOLW0M(L,NY,NX),TK0M(L2,NY,NX),TK0M(L,NY,NX)
        !7757  FORMAT(A8,5I4,30E14.6)
        !     ENDIF
        !
        ! DISCHARGE OF MELTWATER AND ITS HEAT FROM LOWEST SNOWPACK LAYER
        ! TO RESIDUE, SURFACE SOIL MICROPORES AND MACROPORES
        !
        ! FLWQX,FLWQR=porosity-unconstrained water flux to soil,litter
        ! FLWQGX,FLWQGS,FLWQGH=water flux to soil surface,
        ! micropores,macropores
        ! VOLP1,VOLPH1=air volumes of soil micropores,macropores
        ! FMAC,FGRD=macropore,micropore volume fractions
        ! HFLWQG,HFLWQR=convective heat fluxes to soil,litter
        ! THETWR,THETW1=litter, soil water concentration
        ! VOLWRX=litter water retention capacity
        ! PSISM1(0,PSISM1(NUM=litter,soil water potentials
        ! THETY=hygroscopic water concentration
        ! POROS=soil porosity
        ! FC,WP,FCL,WPL=field capacity,wilting point, log(FC),log(WP)
        ! FCI,WPI=FC,WP of ice
        ! THETIX=ice concentration
        ! BKVL=bulk density x volume of soil layer
        !
      ELSE
        !L==JS, top layer
        IF(ICHKL.EQ.0)THEN
          FLWQGX=FLWQX*BARE(NY,NX)
          FLWQGS=AMIN1(VOLP1(NUM(NY,NX),NY,NX)*XNPX &
            ,FLWQGX*FGRD(NUM(NY,NX),NY,NX))
          FLWQGH=AMIN1(VOLPH1(NUM(NY,NX),NY,NX)*XNPX &
            ,FLWQGX*FMAC(NUM(NY,NX),NY,NX))
          FLWQG=FLWQGS+FLWQGH
          HFLWQG=cpw*TK0M(L,NY,NX)*FLWQG
          FLWQR=FLWQX-FLWQG
          HFLWQR=cpw*TK0M(L,NY,NX)*FLWQR
          IF(VOLR(NY,NX).GT.ZEROS(NY,NX) &
            .AND.VOLW1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
            THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
            IF(THETWR.LT.FC(0,NY,NX))THEN
              PSISM1(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
                +((FCL(0,NY,NX)-LOG(THETWR)) &
                /FCD(0,NY,NX)*PSIMD(NY,NX))))
            ELSEIF(THETWR.LT.POROS0(NY,NX))THEN
              PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX) &
                +(((PSL(0,NY,NX)-LOG(THETWR)) &
                /PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
            ELSE
              THETWR=POROS0(NY,NX)
              PSISM1(0,NY,NX)=PSISE(0,NY,NX)
            ENDIF
          ELSE
            THETWR=POROS0(NY,NX)
            PSISM1(0,NY,NX)=PSISE(0,NY,NX)
          ENDIF
          THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX) &
            ,AMIN1(POROS(NUM(NY,NX),NY,NX) &
            ,safe_adb(VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX))))
          IF(BKVL(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
            IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
              PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
                +((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
                /FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
            ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
              PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX) &
                +(((PSL(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
                /PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*PSISD(NY,NX)))
            ELSE
              THETW1=POROS(NUM(NY,NX),NY,NX)
              PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
            ENDIF
          ELSEIF(VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
            FCX=FCI*THETIX(NUM(NY,NX),NY,NX)
            WPX=WPI*THETIX(NUM(NY,NX),NY,NX)
            FCLX=LOG(FCX)
            WPLX=LOG(WPX)
            PSDX=PSL(L,NY,NX)-FCLX
            FCDX=FCLX-WPLX
            IF(THETWX(NUM(NY,NX),NY,NX).LT.FCX)THEN
              PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
                +((FCLX-LOG(THETWX(NUM(NY,NX),NY,NX))) &
                /FCDX*PSIMD(NY,NX))))
            ELSEIF(THETWX(NUM(NY,NX),NY,NX) &
              .LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
              PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX) &
                +(((PSL(NUM(NY,NX),NY,NX)-LOG(THETWX(NUM(NY,NX),NY,NX))) &
                /PSDX)*PSISD(NY,NX)))
            ELSE
              THETW1=POROS(NUM(NY,NX),NY,NX)
              PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
            ENDIF
            !     WRITE(*,1119)'PSISMS',I,J,M,N,NX,NY,NUM(NY,NX)
            !    2,PSISM(NUM(NY,NX),NY,NX),THETW1,VOLW1(NUM(NY,NX),NY,NX)
            !    3,VOLX(NUM(NY,NX),NY,NX)
            !    2,THETW(NUM(NY,NX),NY,NX),THETI(NUM(NY,NX),NY,NX)
            !    3,FCX,WPX,POROS(NUM(NY,NX),NY,NX)
            !1119  FORMAT(A8,7I4,20E12.4)
          ELSE
            THETW1=POROS(NUM(NY,NX),NY,NX)
            PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
          ENDIF
          PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
!
          ! VAPOR FLUX BETWEEN SNOWPACK AND SOIL SURFACE
          !
          ! VOLP01,THETPM=air volume,concentration
          ! CNV1,CNV2=vapor conductances of source, destination layers
          ! VP1,VP2=vapor concentrations of source, destination layers
          ! POROS,POROQ=porosity, tortuosity
          ! WGSGL=vapor diffusivity
          ! TK0M,TK1=snow,soil surface temperature
          ! PSISV1=soil matric+osmotic potential
          ! ATCNVS=snow-soil vapor conductance
          ! DLYR=soil surface layer depth
          ! FLVC,FLVX=vapor flux unlimited,limited by vapor
          ! VPY=equilibrium vapor concentration
          ! XNPX=time step for flux calculations
          ! FLVS1,HFLVS1=vapor flux and its convective heat flux
          !
          IF(VOLP01.GT.ZEROS2(NY,NX).AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
            CNV2=WGSGL(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ &
              *THETPM(M,NUM(NY,NX),NY,NX)/POROS(NUM(NY,NX),NY,NX)
            !VP2=2.173E-03_r8/TK1(NUM(NY,NX),NY,NX) &
            !  *0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TK1(NUM(NY,NX),NY,NX))) &
            VP2=vapsat(TK1(NUM(NY,NX),NY,NX))*EXP(18.0_r8*PSISV1/(8.3143_r8*TK1(NUM(NY,NX),NY,NX)))
            ATCNVS=2.0_r8*CNV1*CNV2 &
              /(CNV1*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRS0(L,NY,NX))
            FLVC=ATCNVS*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX) &
              *FSNW(NY,NX)*BARE(NY,NX)*XNPY
            VPY=(VP1*VOLP01+VP2*VOLPM(M,NUM(NY,NX),NY,NX)) &
              /(VOLP01+VOLPM(M,NUM(NY,NX),NY,NX))
            FLVX=(VP1-VPY)*VOLP01*XNPA
            IF(FLVC.GE.0.0_r8)THEN
              FLVS1=AMAX1(0.0_r8,AMIN1(FLVC,FLVX))
              HFLVS1=(cpw*TK0M(L,NY,NX)+VAP)*FLVS1
            ELSE
              FLVS1=AMIN1(0.0_r8,AMAX1(FLVC,FLVX))
              HFLVS1=(cpw*TK1(NUM(NY,NX),NY,NX)+VAP)*FLVS1
            ENDIF
          ELSE
            CNV2=0.0_r8
            FLVS1=0.0_r8
            HFLVS1=0.0_r8
          ENDIF
          !
          ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE SOIL
          !
          ! WTHET2=multiplier for air concentration in thermal conductivity
          ! TCND1W,TCNDS=thermal conductivity of snowpack, soil surface
          ! STC,DTC=mineral component of thermal conductivity
          ! THETWX,THETIX,THETPX=soil surface water,ice,air concentrations
          ! BAREW=soil surface fraction
          ! ATCNDS=snowpack-soil thermal conductance
          ! TKWX1=interim snowpack temperature
          ! TKY=equilibrium temperature
          ! HFLWX,HFLWC=heat-constrained,heat-unconstrained heat fluxes
          ! XNPY=time step for snowpack flux calculations
          ! HFLWS1=snowpack-soil heat flux
          !
          WTHET2=1.467_r8-0.467_r8*THETPY(NUM(NY,NX),NY,NX)
          TCNDS=(STC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX) &
            *2.067E-03_r8+0.611_r8*THETIX(NUM(NY,NX),NY,NX)*7.844E-03_r8 &
            +WTHET2*THETPX(NUM(NY,NX),NY,NX)*9.050E-05_r8) &
            /(DTC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX) &
            +0.611_r8*THETIX(NUM(NY,NX),NY,NX) &
            +WTHET2*THETPX(NUM(NY,NX),NY,NX))
          IF(BARE(NY,NX).GT.ZERO)THEN
            ATCNDS=2.0_r8*TCND1W*TCNDS/(TCND1W*DLYR(3,NUM(NY,NX),NY,NX) &
              +TCNDS*DLYRS0(L,NY,NX))
          ELSE
            ATCNDS=0.0_r8
          ENDIF
          TKWX1=TK1(NUM(NY,NX),NY,NX)+HFLVS1/VHCP1(NUM(NY,NX),NY,NX)
          TKY=(TK0M(L,NY,NX)*VHCPWMM(L,NY,NX)+TKWX1*VHCP1(NUM(NY,NX),NY,NX)) &
            /(VHCPWMM(L,NY,NX)+VHCP1(NUM(NY,NX),NY,NX))
          HFLWX=(TK0M(L,NY,NX)-TKY)*VHCPWMM(L,NY,NX)*XNPA
          HFLWC=ATCNDS*(TK0M(L,NY,NX)-TKWX1)*AREA(3,NUM(NY,NX),NY,NX) &
            *FSNW(NY,NX)*BARE(NY,NX)*XNPY
          IF(HFLWC.GE.0.0_r8)THEN
            HFLWS1=AMAX1(0.0_r8,AMIN1(HFLWX,HFLWC))
          ELSE
            HFLWS1=AMIN1(0.0_r8,AMAX1(HFLWX,HFLWC))
          ENDIF
          ! IF(J.EQ.15.AND.M.EQ.NPH)THEN
          ! WRITE(*,1113)'HFLWS1',I,J,M,MM,L,FLVS1,FLVX,HFLVS1
          ! 2,HFLWS1,ATCNVS,VP1,VP2,CNV1,CNV2,PSISV1
          ! 3,HFLWX,HFLWC,ATCNDS,TKW(L,NY,NX),TK1(NUM(NY,NX),NY,NX)
          ! 4,THETPX(NUM(NY,NX),NY,NX),WGSGL(NUM(NY,NX),NY,NX)
          ! 5,VHCPWMM(L,NY,NX),TCND1W,TCNDS,PSISV1,TKY,TK0M(L,NY,NX),TKWX1
          ! 6,VOLP1(NUM(NY,NX),NY,NX),VOLPH1(NUM(NY,NX),NY,NX)
          ! 6,VOLT(NUM(NY,NX),NY,NX),VOLA1(NUM(NY,NX),NY,NX)
          ! 7,VOLW1(NUM(NY,NX),NY,NX),VOLI1(NUM(NY,NX),NY,NX)
          ! 8,POROS(NUM(NY,NX),NY,NX)
          !1113  FORMAT(A8,5I4,60E14.6)
          ! ENDIF
          !
          ! HEAT FLUX AMONG SNOWPACK, SURFACE RESIDUE AND SURFACE SOIL
          !
          ! FLVSR=snowpack-litter vapor flux
          ! HFLVSR,HFLWSR=snowpack-litter convective,conductive heat fluxes
          ! FLVS1=snowpack-soil vapor flux
          ! HFLVS1,HFLWS1=snowpack-soil convective,conductive heat fluxes
          ! VHCP1,VHCPRX=current,minimum litter heat capacities
          ! TK0X,TKXR,TK1X=snowpack,litter,soil temperatures
          ! CNVR,CNV1,CNV2=litter,snowpack,soil vapor conductivity
          ! THETP*,THETWX,THETIX=litter air,water,ice concentration
          ! POROS,POROQ=litter porosity, tortuosity
          ! CVRD=litter cover fraction
          ! WGSGR=litter vapor diffusivity
          ! ATCNVR,ATCNVS=snowpack-litter,litter-soil vapor conductance
          ! DLYRR,DLYRS0,DLYR=litter,snowpack,soil depths
          ! THETRR=dry litter concentration
          ! TCNDR,TCND1W,TCNDS=litter,snowpack,soil thermal conductivity
          ! ATCNDR,ATCNDS=snow-litter,litter-soil thermal conductance
          !
          FLVSR=0.0_r8
          HFLVSR=0.0_r8
          HFLWSR=0.0_r8
          FLVR1=0.0_r8
          HFLVR1=0.0_r8
          HFLWR1=0.0_r8
          IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
            TK0X=TK0M(L,NY,NX)
            TKXR=TK1(0,NY,NX)
            TK1X=TK1(NUM(NY,NX),NY,NX)
            CNVR=WGSGR(NY,NX)*THETPM(M,0,NY,NX)*POROQ &
              *THETPM(M,0,NY,NX)/POROS(0,NY,NX)
            IF(CVRD(NY,NX).GT.ZERO)THEN
              IF(CNV1.GT.ZERO.AND.CNVR.GT.ZERO)THEN
                ATCNVR=2.0*CNVR*CNV1/(CNV1*DLYRR(NY,NX)+CNVR*DLYRS0(L,NY,NX))
              ELSE
                ATCNVR=2.0*CNV1/(DLYRR(NY,NX)+DLYRS0(L,NY,NX))
              ENDIF
              IF(CNVR.GT.ZERO.AND.CNV2.GT.ZERO)THEN
                ATCNVS=2.0*CNVR*CNV2 &
                  /(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRR(NY,NX))
              ELSE
                ATCNVS=2.0*CNV2/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))
              ENDIF
              THETRR=AMAX1(0.0,1.0-THETPX(0,NY,NX)-THETWX(0,NY,NX)-THETIX(0,NY,NX))
              TCNDR=(0.779*THETRR*9.050E-04+0.622*THETWX(0,NY,NX) &
                *2.067E-03+0.380*THETIX(0,NY,NX)*7.844E-03+THETPX(0,NY,NX) &
                *9.050E-05)/(0.779*THETRR+0.622*THETWX(0,NY,NX) &
                +0.380*THETIX(0,NY,NX)+THETPX(0,NY,NX))
              IF(TCND1W.GT.ZERO.AND.TCNDR.GT.ZERO)THEN
                ATCNDR=2.0*TCND1W*TCNDR/(TCND1W*DLYRR(NY,NX)+TCNDR*DLYRS0(L,NY,NX))
              ELSE
                ATCNDR=0.0_r8
              ENDIF
              IF(TCNDR.GT.ZERO.AND.TCNDS.GT.ZERO)THEN
                ATCNDS=2.0*TCNDR*TCNDS/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX)+TCNDS*DLYRR(NY,NX))
              ELSE
                ATCNDS=0.0_r8
              ENDIF
            ELSE
              ATCNVR=0.0_r8
              ATCNVS=0.0_r8
              ATCNDR=0.0_r8
              ATCNDS=0.0_r8
            ENDIF
            !
            ! SHORTER TIME STEP FOR SURFACE RESIDUE FLUX CALCULATIONS
            call SnowSurfaceResidualIteration(L,M,NY,NX)

          ENDIF
          !
          ! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
          ! FOR LATER UPDATES TO STATE VARIABLES
          !
          ! FLWLT,FLWLW=total,accumulated water flux to soil micropores
          ! FLWLXW,FLWHLW=total,accumd snow-soil micropore,macropore water
          ! HFLWLT,HFLWLW=total,accumulated snow+litter heat flux to soil
          ! FLWRT,FLWRLW=total,accumulated snow+soil water flux to litter
          ! HFLWRT,HFLWRLW=total,accumulated snow+soil heat flux to litter
          ! FLQRM,FLQSM,FLQHM=total water flux to litter,soil micropore,macropore
          ! FLSW,FLSWH,FLSWR=water flux from lowest snow layer to soil macropore,micropore,litter
          ! HFLSW,HFLSWR=heat flux from lowest snow layer to soil,litter
!
          FLWLT=FLWQGS+FLVS1+FLVR1
          FLWLW=FLWLW+FLWLT
          if(abs(FLWLW)>1.e20_r8)then
            write(*,*)'FLWLW=',FLWQGS,FLVS1,FLVR1
            write(*,*)'at line',__LINE__
            call endrun(trim(mod_filename)//'at line',__LINE__)
          endif
          FLWLXW=FLWLXW+FLWQGS
          FLWHLW=FLWHLW+FLWQGH
          HFLWLT=HFLWQG+HFLVS1+HFLWS1+HFLVR1+HFLWR1
          HFLWLW=HFLWLW+HFLWLT
          FLWRT=FLWQR+FLVSR-FLVR1
          FLWRLW=FLWRLW+FLWRT
          HFLWRT=HFLWQR+HFLVSR+HFLWSR-HFLVR1-HFLWR1
          HFLWRLW=HFLWRLW+HFLWRT
          FLWVLW=0.0_r8
          FLQRM(M,NY,NX)=FLQRM(M,NY,NX)+FLWQR
          FLQSM(M,NY,NX)=FLQSM(M,NY,NX)+FLWQGS
          FLQHM(M,NY,NX)=FLQHM(M,NY,NX)+FLWQGH
          FLSW(L,NY,NX)=FLSW(L,NY,NX)+FLWLT
          FLSWH(L,NY,NX)=FLSWH(L,NY,NX)+FLWQGH
          HFLSW(L,NY,NX)=HFLSW(L,NY,NX)+HFLWLT
          FLSWR(L,NY,NX)=FLSWR(L,NY,NX)+FLWRT
          HFLSWR(L,NY,NX)=HFLSWR(L,NY,NX)+HFLWRT
          !     IF(I.EQ.53)THEN
          !     WRITE(*,7752)'FLWLW',I,J,M,MM,NX,NY,L
          !    2,FLWLW,FLWLT,FLWQGS,FLVS1,FLVR1,FLSW(L,NY,NX)
          !    2,FLWRLW,FLWRT,FLWQR,FLVSR,FLVR1,FLWQX,FLWQG,FLSWR(L,NY,NX)
          !    2,FLWQX,FLWQGX,BARE(NY,NX),VOLW0M(L,NY,NX),VOLS0M(L,NY,NX)
          !    2,HFLWLW,HFLWLT,HFLWQG,HFLVS1,HFLWS1,HFLVR1,HFLWR1
          !    3,HFLSW(L,NY,NX)
          !    3,HFLWRLW,HFLWRT,HFLWQR,HFLVSR,HFLWSR,HFLVR1,HFLWR1
          !    3,HFLSWR(L,NY,NX)
          !    2,VOLP01,THETPM(M,NUM(NY,NX),NY,NX),THETX
          !    3,CNV2,VP2,TK1(NUM(NY,NX),NY,NX),ATCNVS,FLVC,VPY,FLVX
          !    4,VP1,VP2,TK1X,PSISV1,HFLVR1X,HFLWR1X
          !    5,VHCP1(NUM(NY,NX),NY,NX)
          !    3,HFLWRLW,HFLWLW,VP0,VPR,VPY
          !    3,THETPX(NUM(NY,NX),NY,NX),FLWQX,BARE(NY,NX)
          !    4,VOLW0M(L,NY,NX),VOLS0M(L,NY,NX)
          !    2,HFLWX,HFLWC,ATCNDS,TK0M(L,NY,NX),TKWX1
          !    2,TCND1W,TCNDS,DLYR(3,NUM(NY,NX),NY,NX),DLYRS0(L,NY,NX)
          !    2,THETWX(NUM(NY,NX),NY,NX),THETIX(NUM(NY,NX),NY,NX)
          !    3,WTHET2,THETPX(NUM(NY,NX),NY,NX),VOLP1(NUM(NY,NX),NY,NX)
          !    2,VOLPH1(NUM(NY,NX),NY,NX),VOLA1(NUM(NY,NX),NY,NX)
          !    2,VOLW1(NUM(NY,NX),NY,NX),VOLI1(NUM(NY,NX),NY,NX),BARE(NY,NX)
          !    3,FLQRM(M,NY,NX),FLWQR,FLQSM(M,NY,NX),FLWQGS
          !    4,FLQHM(M,NY,NX),FLWQG,VOLS0M(L,NY,NX),VOLW0M(L,NY,NX)
          !    5,VOLI0M(L,NY,NX),DLYRS0(L,NY,NX),FLWQX,TK0M(L,NY,NX)
          !7752  FORMAT(A8,7I4,40E12.4)
          !     ENDIF
          ICHKL=1
        ENDIF
      ENDIF
    ENDIF
9880  CONTINUE
  end subroutine SnowPackIteration
!------------------------------------------------------------------------------------------

  subroutine SolveSnowpack(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX

  integer :: MM,L,L2
  !     begin_execution
  !     HEAT AND VAPOR FLUXES BETWEEN SNOWPACK AND ATMOSPHERE
  !
  !     VHCPWM=volumetric heat capacity of snowpack
  !     NPS=number of cycles for solving snowpack heat and water fluxes
  !     ALBW=snowpack albedo
  !     VOLS0M,VOLI0M,VOLW0M=snow,ice,water volumes
  !     RFLX0=net radiation input
  !     RADXW=shortwave radiation at snowpack surface
  !     THRYW=longwave radn incident at snowpack surface
  !     THRMX=longwave radn emitted by snowpack surface
  !     TK0M=snowpack surface temperature
  !     RFLXW2=net radiation

  DO 3000 MM=1,NPS

    ALBW=(0.85_r8*VOLS0M(1,NY,NX)+0.30_r8*VOLI0M(1,NY,NX)+0.06_r8*VOLW0M(1,NY,NX)) &
      /(VOLS0M(1,NY,NX)+VOLI0M(1,NY,NX)+VOLW0M(1,NY,NX))
    RFLX0=(1.0-ALBW)*RADXW(NY,NX)+THRYW(NY,NX)    !incoming radiation, short + longwave
    THRMX=THRMW(NY,NX)*TK0M(1,NY,NX)**4           !emitting longwave radiation,
    RFLXW2=RFLX0-THRMX                            !net radiation
    !
    !     AERODYNAMIC RESISTANCE ABOVE SNOWPACK INCLUDING
    !     RESISTANCE IMPOSED BY PLANT CANOPY
    !
    !     RI=Richardsons number
    !     RIB=isothermal RI
    !     TKQ=canopy air temperature
    !     RAGX,RA=snowpack blr
    !     RAG,RAGW=isothermal blrs at ground,snowpack surfaces
    !
    RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB(NY,NX)*(TKQ(NY,NX)-TK0M(1,NY,NX))))
    RAGX=AMAX1(RAM,0.8_r8*RAGW(NY,NX),AMIN1(1.2_r8*RAGW(NY,NX),&
      RAG(NY,NX)/(1.0_r8-10.0_r8*RI)))
    RAGW(NY,NX)=RAGX
    RA=RAGX
    !
    ! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
    !
    !     PARE,PARS=blcs for snowpack latent,sensible heat fluxes
    !     PAREW,PARSW=conductances for latent,sensible heat fluxes
    !     RZ=surface resistance
    !     VP0,VPQ=vapor pressure at snowpack surface, canopy air
    !     EVAPT2,EVAPW2,EVAPS2=evaporation total, water,snow
    !     XNPS=1/NPS
    !     EFLXW2=latent heat flux
    !     VAP,VAPS=latent heat of evaporation,sublimation
    !     VFLXW2=convective heat of evaporation flux
    !
    PARE=PAREW(NY,NX)/(RA+RZ)
    PARS=PARSW(NY,NX)/RA
    !VP0=2.173E-03_r8/TK0M(1,NY,NX)*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TK0M(1,NY,NX)))
    VP0=vapsat(TK0M(1,NY,NX))
    EVAPT2=PARE*(VPQ(NY,NX)-VP0)
    EVAPW2=AMAX1(EVAPT2,-AMAX1(0.0_r8,VOLW0M(1,NY,NX)*XNPA))
    EVAPX2=AMIN1(0.0_r8,EVAPT2-EVAPW2)
    EVAPS2=AMAX1(EVAPX2,-AMAX1(0.0_r8,VOLS0M(1,NY,NX)*XNPA))
    EFLXW2=EVAPW2*VAP+EVAPS2*VAPS
    IF(EVAPT2.LT.0.0_r8)THEN
      VFLXW2=(EVAPW2*cpw+EVAPS2*cps)*TK0M(1,NY,NX)
    ELSE
      VFLXW2=(EVAPW2*cpw+EVAPS2*cps)*TKQ(NY,NX)
    ENDIF
!
!     SOLVE FOR SNOWPACK SURFACE TEMPERATURE AT WHICH ENERGY
!     BALANCE OCCURS, SOLVE AND ACCUMULATE LATENT, SENSIBLE
!     STORAGE HEAT FLUXES AND EVAPORATION
!
!     SFLXW2,EFLXW2,RFLXW2=sensible,latent heat fluxes, net radiation
!     VFLXW2=convective heat flux from EFLXW2
!     HFLX02=storage heat flux
!     FLQ0S2,FLQ0W2,FLQ0I2=snow,water,ice input to snowpack
!     HWFLQ02=convective heat from snow,water,ice input to snowpack
!
    SFLXW2=PARS*(TKQ(NY,NX)-TK0M(1,NY,NX))
    HFLX02=RFLXW2+EFLXW2+SFLXW2
    HFLXW2=HFLX02+VFLXW2
    RFLXW=RFLXW+RFLXW2
    EFLXW=EFLXW+EFLXW2
    VFLXW=VFLXW+VFLXW2
    SFLXW=SFLXW+SFLXW2
    HFLXW=HFLXW+HFLXW2
    EVAPS(NY,NX)=EVAPS(NY,NX)+EVAPS2
    EVAPW(NY,NX)=EVAPW(NY,NX)+EVAPW2
    FLQ0S2=FLQ0S(NY,NX)*XNPS
    FLQ0W2=FLQ0W(NY,NX)*XNPS
    FLQ0I2=FLQ0I(NY,NX)*XNPS
    HWFLQ02=HWFLQ0(NY,NX)*XNPS
    FLW0S2=FLQ0S2+EVAPS2
    FLW0W2=FLQ0W2+EVAPW2
    FLW0I2=FLQ0I2
    HFLW0W2=HWFLQ02+HFLXW2
    FLW0S(1,NY,NX)=FLW0S2
    FLW0W(1,NY,NX)=FLW0W2
    FLW0I(1,NY,NX)=FLW0I2
    HFLW0W(1,NY,NX)=HFLW0W2
    FLQWM(M,1,NY,NX)=FLQWM(M,1,NY,NX)+FLQ0S2+FLQ0I2+FLQ0W2
!     IF(NX.EQ.3.AND.NY.EQ.3)THEN
!     WRITE(*,7759)'EVAP',I,J,M,MM,FLW0S2
!    2,FLQ0S2,EVAPS2,FLW0W2,FLQ0W2
!    3,FSNW(NY,NX),FLW0I2,FLQ0I2,RFLXW2,EFLXW2
!    4,SFLXW2,VFLXW2,RA,EVAPT2,EVAPX2,VPQ(NY,NX),VP0
!    5,VOLW0M(1,NY,NX),VOLS0M(1,NY,NX),VOLI0M(1,NY,NX)
!    6,HFLW0W(1,NY,NX),HFLW0W2,HWFLQ02,HFLXW2,RFLXW2,EFLXW2
!    7,SFLXW2,VFLXW2,TK0M(1,NY,NX),TKQ(NY,NX)
!    8,PARE,RA,RZ,EVAPS2,EVAPW2,EVAPT2
!7759  FORMAT(A8,4I4,40E14.6)
!     ENDIF
!

    call SnowPackIteration(M,NY,NX,ICHKL)
!
!     ACCUMULATE SNOWPACK FLUXES TO LONGER TIME STEP FOR
!     LITTER, SOIL FLUX CALCULATIONS
!
!     THRMG=total longwave emission
!     XFLWS,XFLWW,XFLWI=hourly accumulated snow,water,ice transfer
!     XHFLWW=hourly convective heat flux from snow,water,ice transfer
!     TFLWSX,TFLWWX,TFLWIX=net snow,water,ice transfer
!     THFLWWX=convective heat flux from net snow,water,ice transfer
!     TFLWS,TFLWW,TFLWI=accumulated net snow,water,ice transfer
!     THFLWW=convective heat flux from accumd snow,water,ice transfer
!
    THRMG(NY,NX)=THRMG(NY,NX)+THRMX
    DO 9860 L=1,JS
      XFLWS(L,NY,NX)=XFLWS(L,NY,NX)+FLW0S(L,NY,NX)
      XFLWW(L,NY,NX)=XFLWW(L,NY,NX)+FLW0W(L,NY,NX)
      XFLWI(L,NY,NX)=XFLWI(L,NY,NX)+FLW0I(L,NY,NX)
      XHFLWW(L,NY,NX)=XHFLWW(L,NY,NX)+HFLW0W(L,NY,NX)
      L2=MIN(JS,L+1)
!
!     IF WITHIN SNOWPACK
!
      IF(L.LT.JS.AND.VHCPWMM(L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        TFLWSX=FLW0S(L,NY,NX)-FLW0S(L2,NY,NX)
        TFLWWX=FLW0W(L,NY,NX)-FLW0W(L2,NY,NX)
        TFLWIX=FLW0I(L,NY,NX)-FLW0I(L2,NY,NX)
        THFLWWX=HFLW0W(L,NY,NX)-HFLW0W(L2,NY,NX)
        TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX
        TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX
        TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
        THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
!
!     IF AT BOTTOM OF SNOWPACK
!
      ELSEIF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
        TFLWSX=FLW0S(L,NY,NX)
        TFLWWX=FLW0W(L,NY,NX)-FLWRT-FLWLT-FLWQGH
        TFLWIX=FLW0I(L,NY,NX)
        THFLWWX=HFLW0W(L,NY,NX)-HFLWRT-HFLWLT
        TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX
        TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX
        TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
        THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
      ELSE
        TFLWSX=0.0_r8
        TFLWWX=0.0_r8
        TFLWIX=0.0_r8
        THFLWWX=0.0_r8
      ENDIF
!     IF(L.EQ.1)THEN
!     WRITE(*,7763)'TFLWW',I,J,M,MM,NX,NY,L,L2,TFLWW(L,NY,NX)
!    2,TFLWWX,FLW0W(L,NY,NX),FLW0W(L2,NY,NX),FLWRLW,FLWLW,FLWHLW
!    3,VOLW0(L,NY,NX),FLWRT,FLWLT,FLWQGH
!    2,THFLWW(L,NY,NX),THFLWWX
!    3,HFLW0W(L,NY,NX),HFLW0W(L2,NY,NX),HFLWRLW,HFLWLW
!    4,VHCPWMM(L,NY,NX)
!7763  FORMAT(A8,8I4,30E14.6)
!     ENDIF
!
!     FREEZE-THAW IN SNOWPACK FROM NET CHANGE IN SNOWPACK
!     HEAT STORAGE
!
!     VOLS0M,VOLW0M,VOLI0M=snow,water,ice volume
!     VHCPWMM,VHCPWMX,VHCPWX=previous,current,minimum heat capacity
!     TK0M=snowpack temperature
!     THFLWWX=net conductive+convective heat flux
!     TFLX1=unconstrained latent heat flux from freeze-thaw
!     FVOLS0,FVOLI0=fractions of total water in water,ice
!     TFLX0X=source-limited latent heat flux from freeze-thaw
!     WFLXSX,WFLXIX=freeze-thaw changes in water,ice
!     WFLXS,WFLXI=accumulated freeze-thaw
!     TFLX0=accumulated latent heat flux from freeze-thaw
!     XWFLXS,XWFLXI=hourly accumulated freeze-thaw
!     XTHAWW=hourly accumulated latent heat flux from freeze-thaw
!
      VOLS0X=AMAX1(0.0_r8,VOLS0M(L,NY,NX))
      VOLW0X=AMAX1(0.0_r8,VOLW0M(L,NY,NX))
      VOLI0X=AMAX1(0.0_r8,VOLI0M(L,NY,NX))
      ENGY0=VHCPWMM(L,NY,NX)*TK0M(L,NY,NX)
      VHCPWMX=cps*VOLS0X+cpw*VOLW0X+cpi*VOLI0X
      IF(VHCPWMX.GT.VHCPWX(NY,NX))THEN
        TK0X=(ENGY0+THFLWWX)/VHCPWMX
        IF((TK0X.LT.TFice.AND.VOLW0X.GT.ZERO*VOLS(NY,NX)) &
          .OR.(TK0X.GT.TFice.AND.VOLI0X+VOLS0X.GT.ZERO*VOLS(NY,NX)))THEN
          TFLX1=VHCPWMX*(TFice-TK0X)/2.7185*XNPX
          IF(TFLX1.LT.0.0_r8)THEN
            TVOLWS=VOLS0X+VOLI0X*DENSI
            IF(TVOLWS.GT.ZEROS2(NY,NX))THEN
              FVOLS0=VOLS0X/TVOLWS
              FVOLI0=VOLI0X*DENSI/TVOLWS
            ELSE
              FVOLS0=0.0_r8
              FVOLI0=0.0_r8
            ENDIF
            TFLX0X=AMAX1(-333.0*TVOLWS*XNPX,TFLX1)
            WFLXSX=-TFLX0X*FVOLS0/333.0
            WFLXIX=-TFLX0X*FVOLI0/333.0
          ELSE
            FVOLS0=0.0_r8
            FVOLI0=0.0_r8
            TFLX0X=AMIN1(333.0*VOLW0X*XNPX,TFLX1)
            WFLXSX=0.0_r8
            WFLXIX=-TFLX0X/333.0
          ENDIF
        ELSE
          TFLX1=0.0_r8
          FVOLS0=0.0_r8
          FVOLI0=0.0_r8
          TFLX0X=0.0_r8
          WFLXSX=0.0_r8
          WFLXIX=0.0_r8
        ENDIF
        WFLXS(L,NY,NX)=WFLXS(L,NY,NX)+WFLXSX
        WFLXI(L,NY,NX)=WFLXI(L,NY,NX)+WFLXIX
        TFLX0(L,NY,NX)=TFLX0(L,NY,NX)+TFLX0X
        XWFLXS(L,NY,NX)=XWFLXS(L,NY,NX)+WFLXSX
        XWFLXI(L,NY,NX)=XWFLXI(L,NY,NX)+WFLXIX
        XTHAWW(L,NY,NX)=XTHAWW(L,NY,NX)+TFLX0X
      ELSE
        TFLX0X=0.0_r8
        WFLXSX=0.0_r8
        WFLXIX=0.0_r8
      ENDIF
!     IF(L.EQ.5)THEN
!     WRITE(*,7758)'TFLX0',I,J,M,MM,NX,NY,L,TK0M(L,NY,NX)
!    4,TFLX0(L,NY,NX),WFLXS(L,NY,NX),WFLXI(L,NY,NX)
!    4,XTHAWW(L,NY,NX),XWFLXS(L,NY,NX),XWFLXI(L,NY,NX)
!    2,TK0X,TKW(L,NY,NX),VHCPWMX,TFLX1,VOLS0X
!    3,VOLW0X,VOLI0X,TFLWS(L,NY,NX),TFLWW(L,NY,NX)
!    4,TFLWI(L,NY,NX),FVOLS0,FVOLI0
!    5,TFLWW(L,NY,NX),THFLWW(L,NY,NX),FLW0W(L,NY,NX)
!7758  FORMAT(A8,7I4,30E14.6)
!     ENDIF
!
!     INTERNAL SNOWPACK SNOW, WATER, ICE, TEMPERATURE
!
!     VOLS0M,VOLW0M,VOLI0M=snow,water,ice volume
!     TFLWSX,TFLWWX,TFLWIX=net snow,water,ice transfer
!     THFLWWX=conductive+convective heat from snow,water,ice transfer
!     WFLXSX,WFLXIX=freeze-thaw changes in water,ice
!     TFLX0X=source-limited latent heat flux from freeze-thaw
!     DENSI=ice density
!     TK0M,TKA=snowpack,air temperature
!     VHCPWMM,VHCPWX=snowpack, minimum heat capacity
!
      VOLS0M(L,NY,NX)=VOLS0M(L,NY,NX)+TFLWSX-WFLXSX
      VOLW0M(L,NY,NX)=VOLW0M(L,NY,NX)+TFLWWX+WFLXSX+WFLXIX
      VOLI0M(L,NY,NX)=VOLI0M(L,NY,NX)-WFLXIX/DENSI
      ENGY0=VHCPWMM(L,NY,NX)*TK0M(L,NY,NX)
      VHCPWMM(L,NY,NX)=cps*VOLS0M(L,NY,NX)+cpw*VOLW0M(L,NY,NX) &
        +cpi*VOLI0M(L,NY,NX)
      IF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
        TK0M(L,NY,NX)=(ENGY0+THFLWWX+TFLX0X)/VHCPWMM(L,NY,NX)
      ELSEIF(L.EQ.1)THEN
        TK0M(L,NY,NX)=TKA(NY,NX)
      ELSE
        TK0M(L,NY,NX)=TK0M(L-1,NY,NX)
      ENDIF
!     IF(L.EQ.5)THEN
!     WRITE(*,7758)'TK0M',I,J,M,MM,NX,NY,L,TK0M(L,NY,NX)
!    2,THFLWWX,TFLX0X
!    3,VOLS0M(L,NY,NX),VOLW0M(L,NY,NX),VOLI0M(L,NY,NX),VOLS1(L,NY,NX)
!    2,WFLXSX,WFLXIX
!    3,TFLWSX,TFLWWX,TFLWIX
!    3,XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
!    4,XHFLWW(L,NY,NX),VHCPWMM(L,NY,NX)
!     ENDIF
9860  CONTINUE
3000  CONTINUE
  end subroutine SolveSnowpack
!------------------------------------------------------------------------------------------

  subroutine SoilSurfaceResidualIteration(M,NY,NX)

  implicit none
  integer, intent(in) :: M,NY,NX

  integer  :: NN
  real(r8) :: tk1pre
! begin_execution

  DO 5000 NN=1,NPR
    IF(VHCPR2.GT.VHCPRX(NY,NX))THEN
      !
      ! AERODYNAMIC RESISTANCE ABOVE RESIDUE INCLUDING
      ! RESISTANCE IMPOSED BY PLANT CANOPY
      !
      ! RI=Richardsons number
      ! RIB=isothermal RI
      ! TKQ,TKR1=canopy air,litter temperature
      ! RZ=surface resistance to evaporation, given as a prescribed parameter
      ! RAGX,RA=litter blr
      ! RAG,RAGR=isothermal blr at ground surface
      ! PARE,PARS=blcs for litter latent,sensible heat fluxes
      !
      RI=AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB(NY,NX)*(TKQ(NY,NX)-TKR1)))
      RAGX=AMAX1(RAM,0.8_r8*RAGR(NY,NX),AMIN1(1.2_r8*RAGR(NY,NX) &
        ,RARG(NY,NX)/(1.0_r8-10.0_r8*RI)))
      RAGR(NY,NX)=RAGX
      RA=RAGX
      PARE=PARER(NY,NX)/(RA+RZ)
      PARS=PARSR(NY,NX)/RA
!
!     NET RADIATION AT RESIDUE SURFACE
!
!     THRMZ2=longwave radiation emitted by litter
!     RFLXR2=litter net radiation
!     THETWR=litter water content
!     VOLWRX=maximum water retention by litter
!     PSISM1=litter matric water potential
!
      THRMZ2=THRMR(NY,NX)*TKR1**4
      RFLXR2=RFLX0-THRMZ2
      IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
        THETWR=AMIN1(VOLWRX(NY,NX),VOLWR2)/VOLR(NY,NX)
      ELSE
        THETWR=POROS0(NY,NX)
      ENDIF
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX) &
        .AND.VOLW1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
        THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
        IF(THETWR.LT.FC(0,NY,NX))THEN
          PSISM1(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
            +((FCL(0,NY,NX)-LOG(THETWR))/FCD(0,NY,NX)*PSIMD(NY,NX))))
        ELSEIF(THETWR.LT.POROS0(NY,NX))THEN
          PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX)+(((PSL(0,NY,NX)-LOG(THETWR)) &
            /PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
        ELSE
          THETWR=POROS0(NY,NX)
          PSISM1(0,NY,NX)=PSISE(0,NY,NX)
        ENDIF
      ELSE
        THETWR=POROS0(NY,NX)
        PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
!
!     VAPOR FLUX AT RESIDUE SURFACE
!
!     VPR,VP1,VPQ=vapor pressure in litter,soil,canopy air, m3/m-3
!     TKS1=soil temperature
!     EVAPR2=litter evaporation
!     EFLXR2=litter latent heat flux
!     VAP=latent heat of evaporation
!     VFLXR2=convective heat of evaporation flux
!
      !VPR=2.173E-03_r8/TKR1*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TKR1)) &
      VPR=vapsat(TKR1)*EXP(18.0_r8*PSISM1(0,NY,NX)/(8.3143_r8*TKR1))   !in litter
      if(abs(VPR)>1.e20_r8)then
        write(*,*)'TKR1=',TKR1,TK1(0,NY,NX),TK1(NUM(NY,NX),NY,NX)
        write(*,*)'PSISM1(0,NY,NX)=',PSISM1(0,NY,NX)
        call endrun(trim(mod_filename)//'at line',__LINE__)
      endif

      !VP1=2.173E-03_r8/TKS1*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TKS1)) &
      VP1=vapsat(TKS1)*EXP(18.0_r8*PSISV1/(8.3143_r8*TKS1))    !in soil
      EVAPR2=AMAX1(-AMAX1(0.0_r8,VOLWR2*XNPX),PARE*(VPQ(NY,NX)-VPR))
!      if(EVAPR2/=EVAPR2.OR.abs(EVAPR2)>1.e20_r8)then
!        write(*,*)'curday, curhour=',curday,curhour
!        write(*,*)'NN=',NN,M
!        write(*,*)'VOLWR2=',VOLWR2
!        write(*,*)'XNPX=',XNPX
!        write(*,*)'PARE=',PARE
!        write(*,*)'VPQ(NY,NX)=',VPQ(NY,NX)
!        write(*,*)'VPR=',VPR
!        write(*,*)'NUM=',NUM(NY,NX)
!        write(*,*)'TKA=',TKA(NY,NX)
!        write(*,*)'TKR1=',TKR1,TKS1,TK1(NUM(NY,NX),NY,NX)
!        write(*,*)'PSISM1=',PSISM1(0,NY,NX)
!        write(*,*)'ORGC=',ORGC(0,NY,NX)
!        call endrun(msg='NaN encounterd in '//mod_filename,line=__LINE__)
!      endif
      EFLXR2=EVAPR2*VAP             !energy flux
      VFLXR2=EVAPR2*cpw*TKR1        !energy flux
      !
      ! SOLVE FOR RESIDUE TO SOIL SURFACE HEAT FLUXES
      !
      ! FLVC,FLVX=vapor unconstrained,vapor constrained vapor flux
      ! XNPZ=time step for litter flux calculations
      ! VPY=equilibrium vapor concentration
      ! VOLPM=litter,soil air filled porosity
      ! FLV2=litter soil vapor flux
      ! HWFLV2=convective heat of litter soil vapor flux
      ! TKXR,TK1X=interim calculation of litter,soil temperatures
      ! TKY=equilibrium litter-soil temperature
      ! HFLWC,HFLWX=litter-soil heat flux unlimited,limited by heat
      ! HFLCR2=litter-soil heat flux
      ! THETPM: air-filled porosity
      IF(THETPM(M,0,NY,NX).GT.THETX.AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
        FLVC=ATCNVR*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*CVRD(NY,NX)*XNPZ
        if(abs(FLVC)>1.e20_r8)then
          write(*,*)'ATCNVR*(VPR-VP1)=',ATCNVR,VPR,VP1
          write(*,*)'FSNX(NY,NX),CVRD(NY,NX)=',FSNX(NY,NX),CVRD(NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
        VPY=(VPR*VOLPM(M,0,NY,NX)+VP1*VOLPM(M,NUM(NY,NX),NY,NX)) &
          /(VOLPM(M,0,NY,NX)+VOLPM(M,NUM(NY,NX),NY,NX))
        FLVX=(VPR-VPY)*VOLPM(M,0,NY,NX)*XNPB
        if(abs(FLVX)>1.e20_r8)then
          write(*,*)'(VPR-VPY)*VOLPM(M,0,NY,NX)',VPR,VPY,VOLPM(M,0,NY,NX)
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
        IF(FLVC.GE.0.0)THEN
          !from litter to soil
          FLV2=AMAX1(0.0,AMIN1(FLVC,FLVX))
          if(abs(FLV2)>1.e20_r8)then
            write(*,*)'FLVC,FLVX',FLVC,FLVX
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          HWFLV2=(cpw*TKR1+VAP)*FLV2
        ELSE
          !from soil to litter
          FLV2=AMIN1(0.0,AMAX1(FLVC,FLVX))
          if(abs(FLV2)>1.e20_r8)then
            write(*,*)'FLVC,FLVX',FLVC,FLVX
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          HWFLV2=(cpw*TKS1+VAP)*FLV2
        ENDIF
      ELSE
        FLV2=0.0_r8
        HWFLV2=0.0_r8
      ENDIF
      TKXR=TKR1-HWFLV2/VHCPR2         !update litter layer temperature
      TK1X=TKS1+HWFLV2/VHCP12         !update soil layer temperature
      TKY=(TKXR*VHCPR2+TK1X*VHCP12)/(VHCPR2+VHCP12)   !equilibrium temperature
      HFLWX=(TKXR-TKY)*VHCPR2*XNPB    !sensible heat flux
      HFLWC=ATCNDR*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)*CVRD(NY,NX)*XNPZ
      ! if(curday>=377)then
      !   write(*,*)'------'
      !   write(*,*)'curhour=',curhour
      !   write(*,*)'HFLWC =',HFLWC
      !   write(*,*)'TKXR  =',TKXR
      !   write(*,*)'TK1X  =',TK1X,NUM(NY,NX)
      !   write(*,*)'TKR1  =',TKR1,TK1(0,NY,NX)
      !   write(*,*)'HWFLV2=',HWFLV2
      !   write(*,*)'VHCPR2=',VHCPR2
      !   write(*,*)'TKS1  =',TKS1
      !   write(*,*)'VHCP12=',VHCP12
      ! endif
      IF(HFLWC.GE.0.0)THEN
        HFLCR2=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
        HFLCR2=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
!
!     SOLVE FOR RESIDUE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
!     SFLXR2,RFLXR2,EFLXR2=litter sensible,net radn,latent heat fluxes
!     HFLX02,HFLXR2=storage,total litter heat flux
!
      SFLXR2=PARS*(TKQ(NY,NX)-TKR1)    !sensible heat flux between canopy air and litter surface
      HFLX02=RFLXR2+EFLXR2+SFLXR2      !
      HFLXR2=HFLX02+VFLXR2
!      if(curday>=377)then
!        write(*,*)'HFLXR2=',HFLXR2
!        write(*,*)'HFLX02=',HFLX02
!        write(*,*)'VFLXR2=',VFLXR2
!        write(*,*)'RFLXR2=',RFLXR2
!        write(*,*)'EFLXR2=',EFLXR2
!        write(*,*)'SFLXR2=',SFLXR2
!      endif
!
!     AGGREGATE WATER AND ENERGY FLUXES FROM TIME STEP FOR LITTER
!     CALCULATIONS TO THAT FOR SOIL PROFILE
!
      EVAPR(NY,NX)=EVAPR(NY,NX)+EVAPR2
      RFLXR=RFLXR+RFLXR2
      EFLXR=EFLXR+EFLXR2
      VFLXR=VFLXR+VFLXR2
      SFLXR=SFLXR+SFLXR2
      HFLXR=HFLXR+HFLXR2
      ! write(*,*)'.........HFLXR',HFLXR,HFLXR2,TKQ(NY,NX),TKA(NY,NX),TKR1
      ! write(*,*)'HFLX02+VFLXR2,NN',HFLX02,VFLXR2,NN
      ! write(*,*)'RFLXR2+EFLXR2+SFLXR2',RFLXR2,EFLXR2,SFLXR2
      ! write(*,*)'PARS*(TKQ(NY,NX)-TKR1)',PARS,TKQ(NY,NX),TKR1
      ! call print_info(trim(mod_filename)//' at line',__LINE__)
      FLV1=FLV1+FLV2
      HWFLV1=HWFLV1+HWFLV2
      HFLCR1=HFLCR1+HFLCR2
      THRMZ=THRMZ+THRMZ2
    ELSE
      EVAPR2=0.0_r8
      RFLXR2=0.0_r8
      EFLXR2=0.0_r8
      VFLXR2=0.0_r8
      SFLXR2=0.0_r8
      HFLXR2=0.0_r8
      FLV2=0.0_r8
      HWFLV2=0.0_r8
      HFLCR2=0.0_r8
      THRMZ2=0.0_r8
    ENDIF
    VOLWR2=VOLWR2+FLYM2+EVAPR2-FLV2
    VOLW12=VOLW12+FLV2
    ENGYR=VHCPR2*TKR1
    ! VHCPR2: heat capacity, kJ/kg/Kelvin
    VHCPR2=cpo*ORGC(0,NY,NX)+cpw*VOLWR2+cpi*VOLI1(0,NY,NX)
    VHCP12=VHCP12+cpw*FLV2
    TKR1=(ENGYR+HFLXR2+HWFLM2-HWFLV2-HFLCR2)/VHCPR2
    !if((TKR1-TK1(NUM(NY,NX),NY,NX))*(TKR1-TKQ(NY,NX))>0._r8)then
     ! VHCP1(1,NY,NX)=VHCM(1,NY,NX)+cpw*(VOLW1(1,NY,NX) &
     !   +VOLWH1(1,NY,NX))+cpi*(VOLI1(1,NY,NX)+VOLIH1(1,NY,NX))
     ! tk1pre=TK1(1,NY,NX)
     ! TKR1=(TKR1*VHCPR2+VHCP1(1,NY,NX)*TK1(1,NY,NX))/(VHCPR2+VHCP1(1,NY,NX))
     ! TK1(1,NY,NX)=TKR1
     ! TK1(0,NY,NX)=TKR1
     ! write(*,*)'curday, NN=',curday,NN
     ! write(*,*)'curhr =',curhour,NPR
     ! write(*,*)'ENGYR =',ENGYR
     ! write(*,*)'HFLXR2=',HFLXR2
     ! write(*,*)'HWFLM2=',HWFLM2
     ! write(*,*)'HWFLV2=',HWFLV2
     ! write(*,*)'HFLCR2=',HFLCR2
     ! write(*,*)'VHCPR2=',VHCPR2
     ! write(*,*)'ORGC(0,NY,NX)=',ORGC(0,NY,NX)
     ! write(*,*)'VOLWR2=',VOLWR2
     ! write(*,*)'VOLI1(0,NY,NX)=',VOLI1(0,NY,NX)
     ! write(*,*)'TKR1  =',TKR1,TK1(0,NY,NX),TK1pre
!      call endrun(trim(mod_filename)//'at line',__LINE__)
    !endif
    TKS1=TKS1+(HWFLV2+HFLCR2)/VHCP12
!   IF(I.GT.350.AND.NX.EQ.1)THEN
!     WRITE(*,1111)'EFLXR2',I,J,M,NX,NY,NUM(NY,NX),NN
!    2,TKR1,TKS1,TKQ(NY,NX),RFLXR2,EFLXR2,SFLXR2,VFLXR2
!    3,HFLX02,HFLXR2,HWFLM2,HWFLV2,HFLCR2,VHCPR2,HFLXR
!    3,RA,RAGX,RAG(NY,NX),RAB(NY,NX),RAC(NY,NX)
!    4,RAR1,PARE,HFLWX,HFLWC,TKXR,TK1X,TKY
!    2,VPQ(NY,NX),VPR,EVAPR(NY,NX),EVAPR2,VOLWR2,XNPX
!    3,HFLCR2,HWFLV2,VHCPR2,VHCP12,RFLX0,THRMZ2,VHCPRX(NY,NX)
!    4,ALBR,RADXR(NY,NX),THRYR(NY,NX),THRMR(NY,NX),XNPZ
!    3,FLV1,FLV2,VPR,VP1,CNVR,CNV1,FLVC,FLVX,XNPZ,XNPR
!    3,PSISM1(0,NY,NX),PSISV1,THETWR,VOLWRX(NY,NX),ORGC(0,NY,NX)
!    4,VHCPRX(NY,NX),PARS,PARE,RA,RZ,RI,TKQ(NY,NX),VOLW1(0,NY,NX)
!    5,VOLW1(NUM(NY,NX),NY,NX),VOLT(NUM(NY,NX),NY,NX),FLV1
!    5,CNVR,CNV1,VOLX(0,NY,NX),POROQ,WGSGR(NY,NX)
!    5,ATCNDR,TCNDR
!    6,TCND1,STC(NUM(NY,NX),NY,NX),THETWX(NUM(NY,NX),NY,NX)
!    2,THETIX(NUM(NY,NX),NY,NX),WTHET1,THETPX(NUM(NY,NX),NY,NX),TCNDA1
!    4,DTC(NUM(NY,NX),NY,NX),VOLP1(0,NY,NX),VOLR(NY,NX)
!    7,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPY(0,NY,NX),ORGC(0,NY,NX)
!    6,DLYR(3,0,NY,NX),DLYR(3,NUM(NY,NX),NY,NX)
!    8,CVRD(NY,NX),XVOLW(NY,NX),VOLWG(NY,NX)
!1111  FORMAT(A8,7I4,100E12.4)
!     ENDIF
5000  CONTINUE
  end subroutine SoilSurfaceResidualIteration
!------------------------------------------------------------------------------------------

  subroutine SoilSurfaceEnergyBalance(M,NY,NX)
!!
! Description:
  implicit none
  integer, intent(in) :: M,NY,NX
! begin_execution
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!
  EVAPR(NY,NX)=0.0_r8
  RFLXR=0.0_r8
  EFLXR=0.0_r8
  VFLXR=0.0_r8
  SFLXR=0.0_r8
  HFLXR=0.0_r8
  FLV1=0.0_r8
  HWFLV1=0.0_r8
  HFLCR1=0.0_r8
  THRMZ=0.0_r8
!
! NET RADIATION AT RESIDUE SURFACE
!
! ALBR=litter albedo
! BKVL=litter mass
! VOLW1,VOLI1=water,ice volume in litter
! RADXR,THRYR=incoming shortwave,longwave radiation
! TKR1,TKS1=litter,soil temperature
! VOLWR2,VOLW12=litter,soil water volume
! VHCPR2,VHCP12=litter,soil heat capacity
!
!
  ALBR=(0.20_r8*BKVL(0,NY,NX)+0.06_r8*VOLW1(0,NY,NX)+0.30_r8 &
    *VOLI1(0,NY,NX))/(BKVL(0,NY,NX)+VOLW1(0,NY,NX)+VOLI1(0,NY,NX))
  RFLX0=(1.0_r8-ALBR)*RADXR(NY,NX)+THRYR(NY,NX)  !radiation incident on litter layer
  TKR1=TK1(0,NY,NX)                              !kelvin, initial litter layer temperature
  VOLWR2=VOLW1(0,NY,NX)                          !volumetric water content
  VHCPR2=VHCP1(0,NY,NX)                          !heat capacity, initialized with residual layer
  TKS1=TK1(NUM(NY,NX),NY,NX)
  VOLW12=VOLW1(NUM(NY,NX),NY,NX)
  VHCP12=VHCP1(NUM(NY,NX),NY,NX)
  !
  ! THERMAL CONDUCTIVITY BETWEEN SURFACE RESIDUE AND SOIL SURFACE
  !
  ! CNVR,CNV1=litter,soil vapor conductivity
  ! THETPM=litter air concentration
  ! POROS,POROQ=litter porosity, tortuosity
  ! WGSGR,WGSGL=litter,soil vapor diffusivity
  ! CVRD=litter cover fraction
  ! ATCNVR=litter-soil vapor conductance
  ! DLYRR,DLYR=litter,soil depths
  ! THETRR=dry litter concentration
  ! DTH*,RYL*,DNU*,TRB*=turbulence effects on thermal conductivity
  ! WTHET0,WTHET1=multiplier for air concn in thermal conductivity
  ! TCNDW*,TCNDA*=thermal conductivity of water,air
  ! TCNDR,TCND1=litter,soil thermal conductivity
  ! ATCNDR=litter-soil thermal conductance
  !
  CNVR=WGSGR(NY,NX)*THETPM(M,0,NY,NX)*POROQ*THETPM(M,0,NY,NX)/POROS(0,NY,NX)
  CNV1=WGSGL(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ &
    *THETPM(M,NUM(NY,NX),NY,NX)/POROS(NUM(NY,NX),NY,NX)
  IF(CVRD(NY,NX).GT.ZERO)THEN
    IF(CNVR.GT.ZERO.AND.CNV1.GT.ZERO)THEN
      ATCNVR=2.0_r8*CNVR*CNV1/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV1*DLYRR(NY,NX))
    ELSE
      ATCNVR=2.0_r8*CNVR/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))*CVRD(NY,NX)
    ENDIF
  ELSE
    ATCNVR=0.0_r8
  ENDIF
  THETRR=AMAX1(0.0_r8,1.0_r8-THETPX(0,NY,NX)-THETWX(0,NY,NX)-THETIX(0,NY,NX))
  DTKX=ABS(TK1(0,NY,NX)-TK1(NUM(NY,NX),NY,NX))*1.0E-06_r8
  DTHW0=AMAX1(0.0_r8,THETWX(0,NY,NX)-TRBW)**3
  DTHA0=AMAX1(0.0_r8,THETPX(0,NY,NX)-TRBA)**3
  DTHW1=AMAX1(0.0_r8,THETWX(NUM(NY,NX),NY,NX)-TRBW)**3
  DTHA1=AMAX1(0.0_r8,THETPX(NUM(NY,NX),NY,NX)-TRBA)**3
  RYLXW0=DTKX*DTHW0
  RYLXA0=DTKX*DTHA0
  RYLXW1=DTKX*DTHW1
  RYLXA1=DTKX*DTHA1
  RYLNW0=AMIN1(1.0E+04_r8,RYLXW*RYLXW0)
  RYLNA0=AMIN1(1.0E+04_r8,RYLXA*RYLXA0)
  RYLNW1=AMIN1(1.0E+04_r8,RYLXW*RYLXW1)
  RYLNA1=AMIN1(1.0E+04_r8,RYLXA*RYLXA1)
  XNUSW0=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW0**0.25_r8/DNUSW)
  XNUSA0=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA0**0.25_r8/DNUSA)
  XNUSW1=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW1**0.25_r8/DNUSW)
  XNUSA1=AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA1**0.25_r8/DNUSA)
  TCNDW0=2.067E-03_r8*XNUSW0
  TCNDA0=9.050E-05_r8*XNUSA0
  TCNDW1=2.067E-03_r8*XNUSW1
  TCNDA1=9.050E-05_r8*XNUSA1
  WTHET0=1.467_r8-0.467_r8*THETPY(0,NY,NX)
  WTHET1=1.467_r8-0.467_r8*THETPY(NUM(NY,NX),NY,NX)
  TCNDR=(0.779_r8*THETRR*9.050E-04_r8+0.622_r8*THETWX(0,NY,NX)*TCNDW0 &
    +0.380_r8*THETIX(0,NY,NX)*7.844E-03_r8+WTHET0*THETPX(0,NY,NX)*TCNDA0) &
    /(0.779_r8*THETRR+0.622_r8*THETWX(0,NY,NX) &
    +0.380_r8*THETIX(0,NY,NX)+WTHET0*THETPX(0,NY,NX))

!  if(curday>=40)then
!    write(*,*)'curhour=',curhour
!    write(*,*)'THETIX =',THETIX(NUM(NY,NX),NY,NX)
!    write(*,*)'WTHET1 =',WTHET1
!    write(*,*)'THETPX =',THETPX(NUM(NY,NX),NY,NX)
!  endif
  TCND1=(STC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX)*TCNDW1 &
    +0.611_r8*THETIX(NUM(NY,NX),NY,NX)*7.844E-03_r8 &
    +WTHET1*THETPX(NUM(NY,NX),NY,NX)*TCNDA1) &
    /(DTC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX) &
    +0.611_r8*THETIX(NUM(NY,NX),NY,NX)+WTHET1*THETPX(NUM(NY,NX),NY,NX))
  ATCNDR=2.0_r8*TCNDR*TCND1/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX) &
    +TCND1*DLYRR(NY,NX))
!
! SMALLER TIME STEP FOR SOLVING SURFACE RESIDUE ENERGY EXCHANGE
!
  call SoilSurfaceResidualIteration(M,NY,NX)
  end subroutine SoilSurfaceEnergyBalance
!------------------------------------------------------------------------------------------

  subroutine PrepSoilSurfaceEnerbyBalance(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
! begin_execution
!
! PHYSICAL AND HYDRAULIC PROPERTIES OF SOIL SURFACE INCLUDING
! AIR AND WATER-FILLED POROSITY, AND WATER POTENTIAL USED IN
! FLUX CALCULATIONS
!
! THETW1,THETY=current,hygroscopic water concentration
! POROS=porosity
! VOLW1,VOLXI=volume of micropore water
! VOLXI=soil volume less macropore,rock
! FC,WP=water contents at field capacity,wilting point
! FCL,WPL=log FC,WP
! FCD,PSD=FCL-WPL,log(POROS)-FCL
! PSISM1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
! PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
! PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
! SRP=parameter for deviation from linear log-log water retention
! function from hour1.f
! PSISO=osmotic potential
!
  IF(BKVL(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX) &
      ,AMIN1(POROS(NUM(NY,NX),NY,NX) &
      ,safe_adb(VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX))))
    IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
        +((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
        /FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
    ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX) &
        +(((PSL(NUM(NY,NX),NY,NX)-LOG(THETW1)) &
        /PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*PSISD(NY,NX)))
    ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
    ENDIF
!
!   SURFACE WATER LAYER
!
!   THETIX,THETWX=ice,water concentration
!   FCI,WPI=ice field capacity,wilting point
!   PSISM1=matric water potential
!
  ELSEIF(VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX).and.&
    THETIX(NUM(NY,NX),NY,NX)>ZEROS2(NY,NX))THEN
    FCX=FCI*THETIX(NUM(NY,NX),NY,NX)
    WPX=WPI*THETIX(NUM(NY,NX),NY,NX)
    FCLX=LOG(FCX)
    WPLX=LOG(WPX)
    PSDX=PSL(NU(NY,NX),NY,NX)-FCLX
    FCDX=FCLX-WPLX
    IF(THETWX(NUM(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX) &
        +((FCLX-LOG(THETWX(NUM(NY,NX),NY,NX)))/FCDX*PSIMD(NY,NX))))
    ELSEIF(THETWX(NUM(NY,NX),NY,NX).LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX) &
        +(((PSL(NUM(NY,NX),NY,NX)-LOG(THETWX(NUM(NY,NX),NY,NX))) &
        /PSDX)*PSISD(NY,NX)))
    ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
    ENDIF
!
!   WRITE(*,1119)'PSISMG',I,J,M,N,NX,NY,NUM(NY,NX)
!    2,PSISM(NUM(NY,NX),NY,NX)
!    2,THETW(NUM(NY,NX),NY,NX),THETI(NUM(NY,NX),NY,NX)
!    3,FCX,WPX,POROS(NUM(NY,NX),NY,NX)
  ELSE
    THETW1=POROS(NUM(NY,NX),NY,NX)
    PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
  ENDIF
  PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
!
! IF(NX.EQ.4.AND.NY.EQ.5)THEN
!     WRITE(*,3232)'PSISV1',I,J,M,NX,NY,NUM(NY,NX),PSISV1
!    2,PSISM1(NUM(NY,NX),NY,NX),PSISO(NUM(NY,NX),NY,NX)
!    3,THETWX(NUM(NY,NX),NY,NX),THETW1,POROS(NUM(NY,NX),NY,NX)
!    4,PSL(NUM(NY,NX),NY,NX),LOG(THETW1),PSD(NUM(NY,NX),NY,NX)
!    5,VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX)
!    5,VOLX(NUM(NY,NX),NY,NX)
!    5,SRP(NUM(NY,NX),NY,NX)
!3232  FORMAT(A8,6I4,20E14.6)
! ENDIF
!
! SOIL SURFACE ALBEDO, NET RADIATION
!
! VOLW1,VOLI1=water,ice volume in micopores
! VOLWH1,VOLIH1=water,ice volume in macopores
! ALBG,ALBS=albedo of ground surface,soil
! BKVL=soil mass
! RADXG,THRYG,RFLXG=incoming shortwave,longwave,net radiation
! THRMA,THRMS=emitted longwave radiation, emissivity
! TK1=soil temperature
!
  VOLWXG=VOLW1(NUM(NY,NX),NY,NX)+VOLWH1(NUM(NY,NX),NY,NX)
  VOLIXG=VOLI1(NUM(NY,NX),NY,NX)+VOLIH1(NUM(NY,NX),NY,NX)
  IF(VOLWXG+VOLIXG.GT.ZEROS2(NY,NX))THEN
    ALBG=(ALBS(NY,NX)*BKVL(NUM(NY,NX),NY,NX)+0.06*VOLWXG &
      +0.30*VOLIXG)/(BKVL(NUM(NY,NX),NY,NX)+VOLWXG+VOLIXG)
  ELSE
    ALBG=ALBS(NY,NX)
  ENDIF
  RFLX0=(1.0-ALBG)*RADXG(NY,NX)+THRYG(NY,NX)
  THRMA=THRMS(NY,NX)*TK1(NUM(NY,NX),NY,NX)**4
  RFLXG=RFLX0-THRMA
!
! AERODYNAMIC RESISTANCE ABOVE SOIL SURFACE INCLUDING
! RESISTANCE IMPOSED BY PLANT CANOPY
!
! THETPX*=air-filled porosity of soil
! DFVR=porosity limitation to diffusion through soil
! POROQ=soil tortuosity
! RAR1=porosity-limited litter blr
! RAGZ=combined soil+litter blr
! RI=Richardsons number
! RIB=isothermal RI
! TKQ,TK1=canopy air,soil temperature
! RAGZ,RA=soil+litter blr
! RAGS=isothermal blr at ground surface
!
  THETPX0=AMAX1(ZERO,THETPX(0,NY,NX))
  DFVR=THETPX0*POROQ*THETPX0/POROS(0,NY,NX)
  RAR1=RAG(NY,NX)+RAR(NY,NX)/DFVR
  RI=AMAX1(-0.3,AMIN1(0.075,RIB(NY,NX)*(TKQ(NY,NX)-TK1(NUM(NY,NX),NY,NX))))
  RAGX=AMAX1(RAM,0.8*RAGS(NY,NX),AMIN1(1.2*RAGS(NY,NX),RAR1/(1.0-10.0*RI)))
  RAGS(NY,NX)=RAGX
  RA=RAGR(NY,NX)+RAGS(NY,NX)
! IF(I.EQ.63.AND.NX.EQ.1)THEN
!     WRITE(*,7776)'RAGX',I,J,M,NX,NY,RAGZ,BARE(NY,NX),RAG(NY,NX)
!    2,CVRDW(NY,NX),RAR1,RI,RIB(NY,NX),TKQ(NY,NX),TK1(NUM(NY,NX),NY,NX)
!    3,TK1(0,NY,NX),RAGX,RAM,RAGS(NY,NX),RA
!    4,RAR(NY,NX),DFVR,THETPX0,POROQ,THETPX(0,NY,NX)
!    5,DLYRR(NY,NX),WGSGR(NY,NX)
!7776  FORMAT(A8,5I4,30E12.4)
! ENDIF
!
! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
!
! PARE,PARS=blcs for latent,sensible heat fluxes over soil
! PAREG,PARSG=conductances for latent,sensible heat fluxes
! RZ=minimum surface resistance
! VP1,VPQ=vapor pressure at soil surface, canopy air
! EVAPG=evaporation
! EFLXG=latent heat flux
! XH=rate constant
! VOLW2=soil water volume
! VAP=latent heat of evaporation
! VFLXG=convective heat of evaporation flux
!
  PARE=PAREG(NY,NX)/(RA+RZ)
  PARS=PARSG(NY,NX)/RA
  TKX1=TK1(NUM(NY,NX),NY,NX)
  !VP1=2.173E-03/TKX1 &
  !    *0.61*EXP(5360.0*(3.661E-03-1.0/TKX1)) &
  VP1=vapsat(TKX1)*EXP(18.0*PSISV1/(8.3143*TKX1))
  EVAPG(NY,NX)=AMAX1(PARE*(VPQ(NY,NX)-VP1) &
      ,-AMAX1(0.0,VOLW2(NUM(NY,NX),NY,NX)*XNPX))
  EFLXG=EVAPG(NY,NX)*VAP
  IF(EVAPG(NY,NX).LT.0.0)THEN
    VFLXG=EVAPG(NY,NX)*cpw*TK1(NUM(NY,NX),NY,NX)
  ELSE
    VFLXG=EVAPG(NY,NX)*cpw*TKQ(NY,NX)
  ENDIF
  VOLW2(NUM(NY,NX),NY,NX)=VOLW2(NUM(NY,NX),NY,NX)+EVAPG(NY,NX)
!
! SOLVE FOR SOIL SURFACE TEMPERATURE AT WHICH ENERGY
! BALANCE OCCURS, SOLVE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
!
! SFLXG,EFLXG,RFLXG=sensible,latent heat fluxes, net radiation
! VFLXG=convective heat flux from EFLXG
! HFLXG=storage heat flux
!
  SFLXG=PARS*(TKQ(NY,NX)-TK1(NUM(NY,NX),NY,NX))
  HFLX0=RFLXG+EFLXG+SFLXG
  HFLXG=HFLX0+VFLXG
! IF(I.GT.331)THEN
!     WRITE(*,1112)'EFLXG',I,J,M,NX,NY,NUM(NY,NX),EVAPG(NY,NX)
!    2,TK1(NUM(NY,NX),NY,NX),TK1(0,NY,NX),DLYRR(NY,NX)
!    2,RFLXG,EFLXG,SFLXG,VFLXG,HFLX0,HFLXG,RA,RAC(NY,NX),RAG(NY,NX)
!    3,RAGZ,RAR1,RAGX,RI,RAGS(NY,NX),BKVL(NUM(NY,NX),NY,NX)
!    3,VOLW2(NUM(NY,NX),NY,NX),VOLI1(NUM(NY,NX),NY,NX),RFLX0,ALBG
!    4,RADXG(NY,NX),THRYG(NY,NX),THRMA,THRYW(NY,NX),THS(NY,NX)
!    5,VPQ(NY,NX),VP1,FRADG(NY,NX),THRMCX(NY,NX)
!    6,PSISO(NUM(NY,NX),NY,NX)
!    6,FLQM,PARE,PARS,PARSG(NY,NX),HWFLQM
!    7,THETPY(NUM(NY,NX),NY,NX),RAR(NY,NX),THETPY(0,NY,NX)
!    8,VHCP1(0,NY,NX),VHCPRX(NY,NX),VHCP1(NUM(NY,NX),NY,NX)
!    3,TKQ(NY,NX),BARE(NY,NX),THETWX(NUM(NY,NX),NY,NX)
!    5,PSISM1(NUM(NY,NX),NY,NX),THETW1,VOLW1(NUM(NY,NX),NY,NX)
!1112  FORMAT(A8,6I4,60E12.4)
! ENDIF
!
  end subroutine PrepSoilSurfaceEnerbyBalance
!------------------------------------------------------------------------------------------

  subroutine PrepForIterationM(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
  integer :: L
! begin_execution
! INITIALIZE NET SURFACE FLUX ACCUMULATORS
!
! TQR1,TQS1,TQW1,TQI1=net water and snowpack snow,water,ice runoff
! THQR1,THQS1=net convective heat from surface water and snow runoff
! BAREW,CVRDW=fractions of soil,litter cover including free water+ice
! RAGS= boundary layer resistance at soil surface
! PARG=boundary layer conductance above soil surface
!
  TQR1(NY,NX)=0.0_r8
  THQR1(NY,NX)=0.0_r8
  TQS1(NY,NX)=0.0_r8
  TQW1(NY,NX)=0.0_r8
  TQI1(NY,NX)=0.0_r8
  THQS1(NY,NX)=0.0_r8
  IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
    BAREW(NY,NX)=AMAX1(0.0,BARE(NY,NX)-AMIN1(1.0,AMAX1(0.0,XVOLT(NY,NX)/VOLWD(NY,NX))))
  ELSE
    BAREW(NY,NX)=1.0
  ENDIF
  CVRDW(NY,NX)=1.0-BAREW(NY,NX)
  RAGS(NY,NX)=1.0/(BAREW(NY,NX)/RAGR(NY,NX)+CVRDW(NY,NX)/RAR1)
  PARG(M,NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*XNPH/(RAGS(NY,NX)+RAS)
!
! REDISTRIBUTE INCOMING PRECIPITATION
! BETWEEN RESIDUE AND SOIL SURFACE
!
! BKDS=bulk density
! FLQRS,FLQRH=water flux from soil micropores,macropores to litter
! FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
! VOLP1,VOLPH1=air-filled microporosity,macroporosity
! HFLQR1=convective heat flux from soil to litter
! FLYM,HWFLYM=total water flux, convective heat flux to litter
! FLQM,FLHM=total water flux to soil micropores, macropores
! HWFLQM=total convective heat flux to soil micropores, macropores
! XNPR=time step for litter water,heat flux calculations
!
  IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    FLQRS=AMAX1(0.0,FLQ1(NY,NX)-VOLP1(NUM(NY,NX),NY,NX))
    FLQRH=AMAX1(0.0,FLH1(NY,NX)-VOLPH1(NUM(NY,NX),NY,NX))
    HFLQR1=cpw*TKA(NY,NX)*(FLQRS+FLQRH)
    FLYM=FLY1(NY,NX)+FLQRS+FLQRH
    HWFLYM=HWFLY1(NY,NX)+HFLQR1
    FLQM=FLQ1(NY,NX)-FLQRS
    FLHM=FLH1(NY,NX)-FLQRH
    HWFLQM=HWFLQ1(NY,NX)-HFLQR1
  ELSE
    FLYM=FLY1(NY,NX)
    HWFLYM=HWFLY1(NY,NX)
    FLQM=FLQ1(NY,NX)
    FLHM=FLH1(NY,NX)
    HWFLQM=HWFLQ1(NY,NX)
  ENDIF
  FLYM2=FLYM*XNPR
  HWFLM2=HWFLYM*XNPR
!
! WATER GAS EXCHANGE COEFFICIENTS IN SURFACE LITTER
!
! VOLA1,VOLI1,VOLW1,VOLPM=total,ice-,water-,air-filled porosity
! TFND1=temperature effect on gas diffusivity
! DFGS=rate constant for air-water gas exchange
! Z1R,Z2RW,Z2RD,Z3RX=parameters for litter air-water gas transfers
! XNPD=time step for gas transfer calculations
! TORT=tortuosity for aqueous diffusivity
!
  VOLAT0=VOLA1(0,NY,NX)-VOLI1(0,NY,NX)
  IF(VOLAT0.GT.ZEROS2(NY,NX).AND.VOLPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    THETWA=AMAX1(0.0,AMIN1(1.0_r8,VOLW1(0,NY,NX)/VOLAT0))
    TFND1=(TK1(0,NY,NX)/298.15)**6
    IF(THETWA.GT.Z3R)THEN
      DFGS(M,0,NY,NX)=AMAX1(0.0_r8,TFND1*XNPD/((Z1R**-1)*EXP(Z2RW*(THETWA-Z3R))))
    ELSE
      DFGS(M,0,NY,NX)=AMIN1(1.0_r8,TFND1*XNPD/((Z1R**-1)*EXP(Z2RD*(THETWA-Z3R))))
    ENDIF
  ELSE
    DFGS(M,0,NY,NX)=0.0_r8
  ENDIF
  IF(VOLWRX(NY,NX).GT.ZEROS(NY,NX))THEN
    THETWT=AMIN1(1.0,VOLW(0,NY,NX)/VOLWRX(NY,NX))
  ELSE
    THETWT=1.0
  ENDIF
  TORT(M,0,NY,NX)=0.7*THETWT**2
!
! KINETIC ENERGY OF DIRECT RAINFALL AND THROUGHFALL
!
! PRECD,PRECB=direct,indirect precipn+irrign at soil surface
! ENGYD,ENGYB=energy impact of direct,indirect precipn+irrign at soil surface
! VOLWG=ground surface water retention capacity
! XVOLW=free surface water
! ZT=canopy height
! ENGYPM=total energy impact for use in erosion.f
! ENGYP=cumulative rainfall energy impact on soil surface
! FKSAT=reduction in soil surface Ksat from rainfall energy impact
!
  IF(PRECD(NY,NX).GT.ZERO)THEN
    ENGYD=AMAX1(0.0,8.95+8.44*LOG(PRECM(NY,NX)))
  ELSE
    ENGYD=0.0
  ENDIF
  IF(PRECB(NY,NX).GT.ZERO)THEN
    ENGYB=AMAX1(0.0,15.8*SQRT(AMIN1(2.5,ZT(NY,NX)))-5.87)
  ELSE
    ENGYB=0.0
  ENDIF
  IF(ENGYD+ENGYB.GT.ZERO)THEN
    HV=1.0E+03*AMAX1(0.0,XVOLT(NY,NX)-VOLWG(NY,NX))/AREA(3,NU(NY,NX),NY,NX)
    ENGYPM(M,NY,NX)=(ENGYD*PRECD(NY,NX)+ENGYB*PRECB(NY,NX)) &
      *EXP(-2.0*HV)*BARE(NY,NX)*XNPH
    ENGYP(NY,NX)=ENGYP(NY,NX)+ENGYPM(M,NY,NX)
  ELSE
    ENGYPM(M,NY,NX)=0.0
  ENDIF
  FKSAT=EXP(-2.0E-03*(CSILT(NU(NY,NX),NY,NX) &
      +CCLAY(NU(NY,NX),NY,NX))*ENGYP(NY,NX))
! IF(ENGYD+ENGYB.GT.ZERO)THEN
!     WRITE(*,1117)'FKSAT',I,J,M,NX,NY,FKSAT,ENGYP(NY,NX)
!    2,ENGYPM(M,NY,NX),ENGYD,PRECD(NY,NX),ENGYB,PRECB(NY,NX)
!    3,PRECM(NY,NX),HV,XVOLWM(M,NY,NX),XVOLIM(M,NY,NX)
!    4,XVOLT(NY,NX),VOLWG(NY,NX),ORGC(0,NY,NX),BARE(NY,NX)
!    5,CCLAY(NU(NY,NX),NY,NX),CSILT(NU(NY,NX),NY,NX),ZT(NY,NX)
!1117  FORMAT(A8,5I4,20E12.4)
!  ENDIF
!
!  SNOWPACK FLUX ACCUMULATORS
!
!  TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
!  THFLWW=convective heat fluxes of snow,water,ice in snowpack
!
  DO 9875 L=1,JS
    TFLWS(L,NY,NX)=0.0
    TFLWW(L,NY,NX)=0.0
    TFLWI(L,NY,NX)=0.0
    THFLWW(L,NY,NX)=0.0
9875  CONTINUE
!
! SURFACE FLUX ACCUMULATORS
!
! TWFLXL,TWFLXH=total freeze-thaw in micropores,macropores
! TTFLXL=total latent heat from freeze-thaw
! TFLWL,TFLWHL=net water flux in micropores,macropores
! THFLWL=net heat flux
!
  DO 9885 L=NUM(NY,NX),NL(NY,NX)
    TWFLXL(L,NY,NX)=0.0
    TWFLXH(L,NY,NX)=0.0
    TTFLXL(L,NY,NX)=0.0
    TFLWL(L,NY,NX)=0.0
    TFLWLX(L,NY,NX)=0.0
    TFLWHL(L,NY,NX)=0.0
    THFLWL(L,NY,NX)=0.0
    VOLW2(L,NY,NX)=VOLW1(L,NY,NX)
!
!   GAS EXCHANGE COEFFICIENTS SOIL LAYERS
!
!   VOLA1,VOLI1,VOLW1=total,ice-,water-filled microporosity
!   VOLAH1,VOLIH1,VOLWH1=total,ice-,water-filled macroporosity
!   VOLPM=air-filled porosity
!   TFND1=temperature effect on gas diffusivity
!   DFGS=rate constant for air-water gas exchange
!   Z1S,Z2SW,Z2SD,Z3SX=parameters for soil air-water gas transfers
!   XNPD=time step for gas transfer calculations
!   TORT,TORTH=tortuosity for aqueous diffn in micropores,macropres
!
    VOLWT=VOLW1(L,NY,NX)+VOLWH1(L,NY,NX)
    VOLAT=VOLA1(L,NY,NX)+VOLAH1(L,NY,NX) &
      -VOLI1(L,NY,NX)-VOLIH1(L,NY,NX)
    IF(VOLAT.GT.ZEROS2(NY,NX) &
      .AND.VOLPM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWA=AMAX1(0.0,AMIN1(1.0,VOLWT/VOLAT))
      TFND1=(TK1(L,NY,NX)/298.15)**6
      Z3S=AMAX1(Z3SX,FC(L,NY,NX)/POROS(L,NY,NX))
      IF(THETWA.GT.Z3S)THEN
        DFGS(M,L,NY,NX)=AMAX1(0.0_r8,TFND1*XNPD/((Z1S**(-1))*EXP(Z2SW*(THETWA-Z3S))))
      ELSE
        DFGS(M,L,NY,NX)=AMIN1(1.0_r8,TFND1*XNPD/((Z1S**(-1))*EXP(Z2SD*(THETWA-Z3S))))
      ENDIF
    ELSE
      DFGS(M,L,NY,NX)=0.0
    ENDIF
!   IF(I.EQ.121.AND.L.EQ.2)THEN
!     WRITE(*,3371)'DFGS',I,J,M,NX,NY,L,DFGS(M,L,NY,NX)
!    2,THETWA,VOLWT,VOLAT,VOLW1(L,NY,NX),VOLA1(L,NY,NX)
!    3,VOLWH1(L,NY,NX),VOLAH1(L,NY,NX),BKDS(L,NY,NX),FHOL(L,NY,NX)
!    3,Z3S,Z2S*(THETWA-Z3S),EXP(Z2S*(THETWA-Z3S)),Z1S**-1
!    4,(Z1S**-1)*EXP(Z2S*(THETWA-Z3S))
!3371  FORMAT(A8,6I4,20E14.6)
!   ENDIF
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      THETWT=safe_adb(VOLWM(M,L,NY,NX),VOLY(L,NY,NX))
      TORT(M,L,NY,NX)=0.7*THETWT**2*(1.0-FHOL(L,NY,NX))
    ELSE
      TORT(M,L,NY,NX)=0.7
    ENDIF
    IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWH=VOLWHM(M,L,NY,NX)/VOLAH1(L,NY,NX)
      TORTH(M,L,NY,NX)=AMIN1(1.0,2.8*THETWH**3)*FHOL(L,NY,NX)
    ELSE
      TORTH(M,L,NY,NX)=0.0
    ENDIF
9885  CONTINUE
! IF(NX.EQ.4.AND.NY.EQ.5)THEN
!     WRITE(*,3132)'FLQR1',I,J,M,NX,NY,FLY1(NY,NX),FLQ1(NY,NX)
!    2,VHCPWM(M,1,NY,NX),VHCPWX(NY,NX)
!    2,FLH1(NY,NX),FLYM,FLQM,FLHM,FLQR1
!    3,FMAC(NUM(NY,NX),NY,NX),FGRD(NUM(NY,NX),NY,NX)
!    5,VOLAH1(NUM(NY,NX),NY,NX),FVOLAH,CCLAY(NUM(NY,NX),NY,NX)
!    4,VOLW1(NUM(NY,NX),NY,NX),VOLX(NUM(NY,NX),NY,NX),WP(L,NY,NX)
!    2,VOLT(NUM(NY,NX),NY,NX),VOLAH1(NUM(NY,NX),NY,NX)
!    5,VOLWRX(NY,NX),VOLW1(0,NY,NX),VOLI1(0,NY,NX)
!    6,PSISM1(0,NY,NX),PSISM1(NUM(NY,NX),NY,NX)
!3132  FORMAT(A8,5I4,40E12.4)
! ENDIF
!
! ENERGY EXCHANGE VARIABLES AT SNOW SURFACE IF PRESENT
!
! RFLXW,EFLXW,VFLXW,SFLXW,HFLXW=netradn,latent,convective,sensible
! and storage heat fluxes
! FLWLW,FLWHLW=water from snowpack to soil micropores,macropores
! HFLWLW=conv heat from snowpack to soil micropores,macropores
! FLWRLW=water flux from snowpack to litter
! HFLWRLW=convective heat flux from snowpack to litter
! FLWVLW=snowpack-litter water flux accounting for wetting front
!
  RFLXW=0.0_r8
  EFLXW=0.0_r8
  VFLXW=0.0_r8
  SFLXW=0.0_r8
  HFLXW=0.0_r8
  FLWLW=0.0_r8
  FLWLXW=0.0_r8
  FLWHLW=0.0_r8
  HFLWLW=0.0_r8
  FLWRLW=0.0_r8
  HFLWRLW=0.0_r8
  FLWVLW=0.0_r8
!
! EVAPS,EVAPW=evaporation from soil,snowpack surfaces
! FLQRM,FLQSM,FLQHM=water into litter,soil micropores,micropores for use in trnsfr.f
!
  EVAPS(NY,NX)=0.0_r8
  EVAPW(NY,NX)=0.0_r8
  FLQRM(M,NY,NX)=0.0_r8
  FLQSM(M,NY,NX)=0.0_r8
  FLQHM(M,NY,NX)=0.0_r8
!
! FLUX VARIABLES IN SNOWPACK
!
! TFLX0=latent heat from freeze-thaw
! WFLXS,WFLXI=freeze-thaw between snow,ice and water
! FLW0S,FLW0W,FLW0I=snow,water,ice fluxes
! HFLW0W=convective heat flux from snow,water,ice fluxes
! FLQWM=snowpack water flux
! VOLS0M,VOLW0M,VOLI0M=snow,water,ice contents
! VHCPWMM=volumetric heat capacity
! TK0M=temperature
!
  DO 9765 L=1,JS
    TFLX0(L,NY,NX)=0.0_r8
    WFLXS(L,NY,NX)=0.0_r8
    WFLXI(L,NY,NX)=0.0_r8
    FLW0S(L,NY,NX)=0.0_r8
    FLW0W(L,NY,NX)=0.0_r8
    FLW0I(L,NY,NX)=0.0_r8
    HFLW0W(L,NY,NX)=0.0_r8
    FLQWM(M,L,NY,NX)=0.0_r8
    VOLS0M(L,NY,NX)=VOLS0(L,NY,NX)
    VOLW0M(L,NY,NX)=VOLW0(L,NY,NX)
    VOLI0M(L,NY,NX)=VOLI0(L,NY,NX)
    VHCPWMM(L,NY,NX)=VHCPWM(M,L,NY,NX)
    TK0M(L,NY,NX)=TK0(L,NY,NX)
9765  CONTINUE
!
  end subroutine PrepForIterationM
!------------------------------------------------------------------------------------------

  subroutine SummaryAfterEnergyBalance(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
! begin_execution
!
! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
! FOR LATER UPDATES TO STATE VARIABLES
!
! FLWLG,FLWHLG=water flux from atm to soil micropores,macropores
! HFLWLG=convective heat flux from atm to soil
! FLWRLG=water flux from atm to litter
! HFLWRLG=convective heat flux from atm to litter
! FLWVLS=water flux within soil accounting for wetting front
!
  FLWLG=FLQM+EVAPG(NY,NX)+FLV1
  if(abs(FLWLG)>1.e20_r8)then
    write(*,*)'FLQM+EVAPG(NY,NX)+FLV1',FLQM,EVAPG(NY,NX),FLV1
    call endrun(trim(mod_filename)//' at line',__LINE__)
  endif
  FLWLXG=FLQM+EVAPG(NY,NX)+FLV1
  FLWHLG=FLHM
  HFLWLG=HWFLQM+HFLXG+HWFLV1+HFLCR1
  FLWRLG=FLYM+EVAPR(NY,NX)-FLV1
  HFLWRLG=HWFLYM+HFLXR-HWFLV1-HFLCR1
!  write(*,*)'HFLWRLGHFLWRLGHFLWRLGHFLWRLG',HWFLYM,HFLXR,HWFLV1,HFLCR1
  FLWVLS=(VOLW1(NUM(NY,NX),NY,NX)-VOLWX1(NUM(NY,NX),NY,NX))*XNPH
! IF(IYRC.EQ.1813.AND.I.EQ.231)THEN
!     WRITE(*,7754)'FLWLG',I,J,M,NX,NY
!    2,FLWLG,FLQM,EVAPG(NY,NX),FLV1
!    3,FLWRLG,FLYM,EVAPR(NY,NX),FLV1
!    4,HFLWLG,HWFLQM,HFLXG,HWFLV1,HFLCR1
!    5,HFLWRLG,HWFLYM,HFLXR
!     ENDIF
!
! GENERATE NEW SNOWPACK
!
! XFLWS,XFLWW,XFLWI=hourly snow,water,ice transfer
! FLQ0S,FLQ0W,FLQ0I=snow,water,ice input to snowpack
! XHFLWW=hourly convective heat flux from snow,water,ice transfer
! HWFLQ0=convective heat flux from snow,water,ice to snowpack
!
  IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX).AND.FLQ0S(NY,NX).GT.ZEROS(NY,NX))THEN
    XFLWS(1,NY,NX)=XFLWS(1,NY,NX)+FLQ0S(NY,NX)
    XFLWW(1,NY,NX)=XFLWW(1,NY,NX)+FLQ0W(NY,NX)
    XFLWI(1,NY,NX)=XFLWI(1,NY,NX)+FLQ0I(NY,NX)
    XHFLWW(1,NY,NX)=XHFLWW(1,NY,NX)+HWFLQ0(NY,NX)
!     WRITE(*,4422)'INIT',I,J,FLQ0S(NY,NX),FLQ0W(NY,NX)
!    3,FLQ0I(NY,NX),HWFLQ0(NY,NX),XFLWS(1,NY,NX),XFLWW(1,NY,NX)
!    2,XFLWI(1,NY,NX),XHFLWW(1,NY,NX),HFLWL(3,NUM(NY,NX),NY,NX)
!    3,HFLW(3,NUM(NY,NX),NY,NX),FSNX(NY,NX),VHCP1(NUM(NY,NX),NY,NX)
!    4*TK1(NUM(NY,NX),NY,NX),HFLWRL(NY,NX),HFLWR(NY,NX)
!    5,VHCP1(0,NY,NX)*TK1(0,NY,NX),HEATH(NY,NX),RFLXG,RFLXR,RFLXW
!    2,SFLXG,SFLXR,SFLXW,EFLXG,EFLXR,EFLXW,VFLXG,VFLXR,VFLXW
  ENDIF
  THRMG(NY,NX)=THRMG(NY,NX)+THRMA+THRMZ
  end subroutine SummaryAfterEnergyBalance
!------------------------------------------------------------------------------------------

  subroutine SurfSoilResidueWaterCapillExch(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX

  integer :: K0,K1
! begin_execution
! CNDR,HCNDR=current,saturated litter hydraulic conductivity
! PSISE,PSISM1(0,=air entry,current litter water potential
! VOLW1(0,VOLWRX=current,maximum litter water volume
! CND1,HCND=soil hydraulic conductivity
! FKSAT=reduction in soil surface Ksat from rainfall energy impact
! K1=soil relative water-filled porosity
! THETWX,POROS=soil water content,porosity
! AVCNDR=litter-soil hydraulic conductance
! DLYRR,DLYR=litter,soil thicknesses
! FLQR,FLQL=litter-soil water flux unltd,ltd by water
! XNPH=time step of flux calculations
! CVRD=fraction of litter cover
! PSISM1(NUM=soil matric water potential
! VOLW1(NUM=soil water volume
! HFLQL=convective heat from litter-soil water flux
! FLWL,HFLWL=micropore water,heat flux
! FLWRL,HFLWRL=total litter water,heat flux
! FLWRM=litter-soil water flux for solute transfer in trnsfr.f
! CND1,CNDL=hydraulic conductivity of source,destination layer
! HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
!
  IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
    IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
    ELSE
      THETWR=POROS0(NY,NX)
    ENDIF
    THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX) &
      ,AMIN1(POROS(NUM(NY,NX),NY,NX) &
      ,safe_adb(VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX))))
    K0=MAX(1,MIN(100,INT(100.0*(AMAX1(0.0,POROS0(NY,NX) &
      -THETWR))/POROS0(NY,NX))+1))
    K1=MAX(1,MIN(100,INT(100.0*(AMAX1(0.0,POROS(NUM(NY,NX),NY,NX) &
      -THETW1))/POROS(NUM(NY,NX),NY,NX))+1))
    CNDR=HCND(3,K0,0,NY,NX)
    CND1=HCND(3,K1,NUM(NY,NX),NY,NX)*FKSAT
    AVCNDR=2.0*CNDR*CND1/(CNDR*DLYR(3,NUM(NY,NX),NY,NX) &
      +CND1*DLYRR(NY,NX))
    PSIST0=PSISM1(0,NY,NX)+PSISH(0,NY,NX) &
      +PSISO(0,NY,NX)
    PSIST1=PSISM1(NUM(NY,NX),NY,NX)+PSISH(NUM(NY,NX),NY,NX) &
      +PSISO(NUM(NY,NX),NY,NX)
    FLQX=AVCNDR*(PSIST0-PSIST1) &
      *AREA(3,NUM(NY,NX),NY,NX)*CVRDW(NY,NX)*XNPH
    IF(FLQX.GE.0.0)THEN
      IF(THETWR.GT.THETS(0,NY,NX))THEN
        FLQZ=FLQX+AMIN1((THETWR-THETS(0,NY,NX)) &
          *VOLR(NY,NX),AMAX1(0.0,(THETS(NUM(NY,NX),NY,NX)-THETW1) &
          *VOLY(NUM(NY,NX),NY,NX)))*XNPX
      ELSE
        FLQZ=FLQX
      ENDIF
      FLQR=AMAX1(0.0,AMIN1(FLQZ,VOLW1(0,NY,NX)*XNPX,VOLP1(NUM(NY,NX),NY,NX)))
      FLQ2=AMAX1(0.0,AMIN1(FLQX,VOLW1(0,NY,NX)*XNPX,VOLP1(NUM(NY,NX),NY,NX)))
    ELSE
      IF(THETW1.GT.THETS(NUM(NY,NX),NY,NX))THEN
        FLQZ=FLQX+AMAX1((THETS(NUM(NY,NX),NY,NX)-THETW1) &
          *VOLY(NUM(NY,NX),NY,NX),AMIN1(0.0,(THETWR-THETS(0,NY,NX)) &
          *VOLR(NY,NX)))*XNPX
      ELSE
        FLQZ=FLQX
      ENDIF
      FLQR=AMIN1(0.0,AMAX1(FLQZ,-VOLW1(NUM(NY,NX),NY,NX)*XNPX,-VOLP1(0,NY,NX)))
      FLQ2=AMIN1(0.0,AMAX1(FLQX,-VOLW1(NUM(NY,NX),NY,NX)*XNPX,-VOLP1(0,NY,NX)))
    ENDIF
    IF(VOLP1Z(NUM(NY,NX),NY,NX).LT.0.0)THEN
      FLQR=FLQR+AMIN1(0.0,AMAX1(-VOLW1(NUM(NY,NX),NY,NX)*XNPX &
        ,VOLP1Z(NUM(NY,NX),NY,NX)))
      FLQ2=FLQ2+AMIN1(0.0,AMAX1(-VOLW1(NUM(NY,NX),NY,NX)*XNPX &
        ,VOLP1Z(NUM(NY,NX),NY,NX)))
    ENDIF
    IF(FLQR.GT.0.0)THEN
      HFLQR=cpw*TK1(0,NY,NX)*FLQR
    ELSE
      HFLQR=cpw*TK1(NUM(NY,NX),NY,NX)*FLQR
    ENDIF
    FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
    if(abs(FLWL(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qFLWL(3,NUM(NY,NX),NY,NX)=',FLQR
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
    FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
    HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR
!    if(curday>=41)then
    !  write(*,*)'curhour=======',curday,curhour
    !  write(*,*)'HFLWRL(NY,NX)-HFLQR=',HFLWRL(NY,NX),HFLQR
    !  write(*,*)'tk1[0,1],FLQR',TK1(0,NY,NX),TK1(NUM(NY,NX),NY,NX),FLQR
    !  write(*,*)'at line',__LINE__
!    endif
    FLWRM(M,NY,NX)=FLQR
!   IF(I.GT.350.AND.NX.EQ.1)THEN
!     WRITE(*,4322)'FLQR',I,J,M,NX,NY,NUM(NY,NX),K0,K1
!    2,FLQR,FLWRL(NY,NX),FLWL(3,NUM(NY,NX),NY,NX),FLQX,FLQZ,FLQ2
!    2,PSISM1(0,NY,NX),PSISM1(NUM(NY,NX),NY,NX),PSIST0,PSIST1
!    3,CVRDW(NY,NX),CNDR,CND1,AVCNDR,FKSAT
!    3,POROS0(NY,NX),VOLW1(0,NY,NX)
!    2,VOLP1(NUM(NY,NX),NY,NX),FLWL(3,NUM(NY,NX),NY,NX)
!    3,VOLP1Z(NUM(NY,NX),NY,NX),THETWX(0,NY,NX)
!    4,THETWX(NUM(NY,NX),NY,NX)
!    4,THETS(0,NY,NX),THETS(NUM(NY,NX),NY,NX)
!    4,THETWR,THETW1,XVOLT(NY,NX)
!    3,HFLQR,HFLWRL(NY,NX),HFLWL(3,NUM(NY,NX),NY,NX)
!    6,THETWR,VHCP1(0,NY,NX),VHCPRX(NY,NX)
!    2,FLWLYH,VOLX(NUM(NY,NX),NY,NX),VOLA1(NUM(NY,NX),NY,NX)
!    3,VOLP1(NUM(NY,NX),NY,NX),VOLW1(NUM(NY,NX),NY,NX)
!    3,VOLI1(NUM(NY,NX),NY,NX),VOLP1(0,NY,NX),VOLW1(0,NY,NX)
!    3,VOLI1(0,NY,NX),PSISM1(0,NY,NX),VOLP1Z(NUM(NY,NX),NY,NX)
!    4,VOLX(NUM(NY,NX),NY,NX),VOLAH(NUM(NY,NX),NY,NX)
!    4,PSISM1(NUM(NY,NX),NY,NX),AVCNDR
!    2,VOLAH1(NUM(NY,NX),NY,NX),VOLPH1(NUM(NY,NX),NY,NX)
!    2,VOLWH1(NUM(NY,NX),NY,NX),VOLIH1(NUM(NY,NX),NY,NX)
!4322  FORMAT(A8,8I4,40E12.4)
!     ENDIF
  ELSE
    FLQR=XVOLW(NY,NX)*XNPX
    HFLQR=cpw*TK1(0,NY,NX)*FLQR
    FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
    if(abs(FLWL(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
      write(*,*)'qrFLWL(3,NUM(NY,NX),NY,NX)=',XVOLW(NY,NX),XNPX
      write(*,*)'at line',__LINE__
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif
    HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
    FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
    HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR
!    if(curday>=41)then
    !  write(*,*)'=======HFLWRL(NY,NX)-HFLQR1=',HFLWRL(NY,NX),HFLQR,TK1(0,NY,NX)
    !  write(*,*)'at line',__LINE__
!    endif
    FLWRM(M,NY,NX)=FLQR
!   IF(I.GT.350.AND.NX.EQ.1)THEN
!     WRITE(*,4323)'FLQ0',I,J,M,NX,NY
!    2,FLQR,XVOLW(NY,NX),XNPX,FLWRL(NY,NX)
!4323  FORMAT(A8,5I4,12E12.4)
!   ENDIF
  ENDIF
! IF((I/10)*10.EQ.I)THEN
!     WRITE(*,4321)'HCNDR',I,J,M,NX,NY,K1,AVCND1,CNDR,CND1,DLYRR(NY,NX)
!    2,PSISM1(0,NY,NX),PSISM1(NUM(NY,NX),NY,NX),FLQL,HFLQL
!    3,VOLWR2,EVAPR(NY,NX),VOLWRX(NY,NX)-VOLW1(0,NY,NX)
!    2-VOLI1(0,NY,NX),VOLW1(NUM(NY,NX),NY,NX),VOLW1(0,NY,NX)
!    4,VOLP1(NUM(NY,NX),NY,NX),POROS(NUM(NY,NX),NY,NX)
!    5,VOLWG(NY,NX),FLYM,HCNDR(NY,NX),PSISE(0,NY,NX),PSISM1(0,NY,NX)
!    6,THETWR,VHCP1(0,NY,NX),VHCPRX(NY,NX)
!    7,VOLPH1(NUM(NY,NX),NY,NX),VOLW1(0,NY,NX),VOLWRX(NY,NX)
!4321  FORMAT(A8,6I4,30E12.4)
!     ENDIF
  end subroutine SurfSoilResidueWaterCapillExch
!------------------------------------------------------------------------------------------

  subroutine InfiltrationRunoffPartition(M,NY,NX,N1,N2)
  implicit none
  integer, intent(in) :: M,NY,NX
  integer, intent(out):: N1,N2
!     begin_execution
!     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
!     OF THE LITTER IS EXCEEDED
!
!     VOLPH1=air-filled macroporosity
!     FINHR,HFINHR=water,convective heat from litter to macropores
!     VOLW1(0,VOLWRX=current,maximum litter water volume
!     FLWL,HFLWL=micropore water,heat flux
!     FLWRL,HFLWRL=total litter water,heat flux
!
  IF(VOLPH1(NUM(NY,NX),NY,NX).GT.0.0 &
    .AND.XVOLW(NY,NX).GT.0.0)THEN
    FLQHR=AMIN1(XVOLW(NY,NX)*XNPX,VOLPH1(NUM(NY,NX),NY,NX))
    HFLQHR=FLQHR*cpw*TK1(0,NY,NX)
    FLWHL(3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)+FLQHR
    HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQHR
    FLWRL(NY,NX)=FLWRL(NY,NX)-FLQHR
    HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQHR
!    if(curday>=41)then
!      write(*,*)'HFLWRL(NY,NX)xx=',HFLWRL(NY,NX),HFLQHR
!      write(*,*)'at line',__LINE__
!    endif
!     IF(I.GT.350.AND.NX.EQ.1)THEN
!     WRITE(*,4357)'FLQHR',I,J,M,NX,NY,NUM(NY,NX),FLQHR,FLWRL(NY,NX)
!    2,FLWHL(3,NUM(NY,NX),NY,NX),VOLPH1(NUM(NY,NX),NY,NX)
!    3,XVOLW(NY,NX),VOLW1(0,NY,NX),VOLWRX(NY,NX)
!    4,HFLQHR,HFLWRL(NY,NX),HFLWL(3,NUM(NY,NX),NY,NX),TK1(0,NY,NX)
!4357  FORMAT(A8,6I4,40E12.4)
!     ENDIF
  ENDIF
!
!     FREEZE-THAW IN RESIDUE SURFACE FROM NET CHANGE IN RESIDUE
!     SURFACE HEAT STORAGE
!
!     TFREEZ=litter freezing temperature
!     PSISM1=litter water potential
!     VOLW1*,VOLI1=litter water,ice volume
!     VHCP1*=litter volumetric heat capacity
!     TK1*=litter temperature
!     ORGC=litter organic C
!     HFLWRL=total litter conductive, convective heat flux
!     TFLX1,TFLX=unltd,ltd latent heat from freeze-thaw
!     TFLXR,WFLXR=litter water,latent heat flux from freeze-thaw
!
  TFREEZ=-9.0959E+04/(PSISM1(0,NY,NX)-333.0)
  VOLW1X=AMAX1(0.0,VOLW1(0,NY,NX)+FLWRL(NY,NX))
  ENGYR=VHCP1(0,NY,NX)*TK1(0,NY,NX)
  VHCP1X=cpo*ORGC(0,NY,NX)+cpw*VOLW1X &
  +cpi*VOLI1(0,NY,NX)
  IF(VHCP1X.GT.ZEROS(NY,NX))THEN
    TK1X=(ENGYR+HFLWRL(NY,NX))/VHCP1X
  ELSE
    TK1X=TK1(0,NY,NX)
  ENDIF
  IF((TK1X.LT.TFREEZ &
    .AND.VOLW1(0,NY,NX).GT.ZERO*VOLT(0,NY,NX)) &
    .OR.(TK1X.GT.TFREEZ &
    .AND.VOLI1(0,NY,NX).GT.ZERO*VOLT(0,NY,NX)))THEN
    TFLX1=VHCP1(0,NY,NX)*(TFREEZ-TK1X) &
      /((1.0+TFREEZ*6.2913E-03)*(1.0-0.10*PSISM1(0,NY,NX)))*XNPX
    IF(TFLX1.LT.0.0)THEN
      TFLX=AMAX1(-333.0*DENSI*VOLI1(0,NY,NX)*XNPX,TFLX1)
    ELSE
      TFLX=AMIN1(333.0*VOLW1X*XNPX,TFLX1)
    ENDIF
    TFLXR(NY,NX)=TFLX
    WFLXR(NY,NX)=-TFLX/333.0
!     IF(NX.EQ.1.AND.NY.EQ.1.AND.test_aneb(TFLX1,0.0_r8))THEN
!     WRITE(*,5352)'TFLXR',I,J,M,NX,NY,WFLXR(NY,NX),TFLXR(NY,NX)
!    2,TFLX1,TFREEZ,TK1(0,NY,NX),TK1X
!    2,TFREEZ-TK1X,FLWRL(NY,NX),HFLWRL(NY,NX)
!    2,THETWR,TFLX,VOLI1(0,NY,NX),VOLW1(0,NY,NX),VHCP1(0,NY,NX)
!    3,PSISM1(0,NY,NX),THAWR(NY,NX),HTHAWR(NY,NX)
!5352  FORMAT(A8,5I4,30E12.4)
!     ENDIF
  ELSE
    WFLXR(NY,NX)=0.0
    TFLXR(NY,NX)=0.0
  ENDIF
!
!     THICKNESS OF WATER FILMS IN LITTER AND SOIL SURFACE
!     FROM WATER POTENTIALS FOR GAS EXCHANGE IN TRNSFR.F
!
  IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
    FILM(M,0,NY,NX)=AMAX1(1.0E-06,EXP(-13.650-0.857*LOG(-PSISM1(0,NY,NX))))
  ELSE
    FILM(M,0,NY,NX)=1.0E-03
  ENDIF
  FILM(M,NUM(NY,NX),NY,NX)=AMAX1(1.0E-06 &
    ,EXP(-13.650-0.857*LOG(-PSISM1(NUM(NY,NX),NY,NX))))
!
!     OVERLAND FLOW WHEN WATER STORAGE CAPACITY
!     OF THE SOIL SURFACE PLUS MACROPORES IS EXCEEDED
!
  N1=NX
  N2=NY
!
!     SURFACE WATER FLUX
!
!     N2,N1=NY,NX of source grid cell
!     XVOLT,XVOLW=excess water+ice,water in source grid cell
!     VOLWG=ground surface water retention capacity
!     VX=ponded water volume above surface retention capacity
!     D,R,S,V=depth,perimeter,slope,velocity of runoff
!     DIST=distance between source,destination
!     ZM=surface roughness height for runoff
!     Q=runoff from Mannings equation
!     QRM,QRV=runoff,velocity for erosion, solute transfer
!
  IF(XVOLT(N2,N1).GT.VOLWG(N2,N1))THEN
    VX=XVOLT(N2,N1)-VOLWG(N2,N1)
    D=VX/AREA(3,0,N2,N1)
    R=D/2.828
    V=R**0.67_r8*safe_adb(SQRT(SLOPE(0,N2,N1)),ZM(N2,N1))
    Q=V*D*AREA(3,NUM(N2,N1),N2,N1)*3.6E+03*XNPH
    VOLW1X=AMAX1(0.0,VOLW1(0,N2,N1)+WFLXR(N2,N1))
    QRM(M,N2,N1)=AMIN1(Q,VX*XNPX,VOLW1X*XNPX)*XVOLW(N2,N1)/XVOLT(N2,N1)
    QRV(M,N2,N1)=V
!     IF(I.EQ.232)THEN
!     WRITE(*,5554)'QRINT',I,J,M,N1,N2,QRM(M,N2,N1),QRV(M,N2,N1)
!    2,Q,V,D,VX*XNPX,SLOPE(0,N2,N1),ZM(N2,N1)
!    3,XVOLW(N2,N1),XVOLT(N2,N1),VOLWG(N2,N1),VOLW1(0,N2,N1)
!    4,WFLXR(N2,N1)
!5554  FORMAT(A8,5I4,20E12.4)
!     ENDIF
  ELSE
    QRM(M,N2,N1)=0.0
    QRV(M,N2,N1)=0.0
  ENDIF
  end subroutine InfiltrationRunoffPartition
!------------------------------------------------------------------------------------------

  subroutine LateralHydroExchange(M,NY,NX,NHE,NHW,NVS,NVN,N1,N2)
  implicit none
  integer, intent(in) :: M,NY,NX,NHE,NHW,NVS,NVN
  integer, intent(in) :: N1,N2
  integer :: N,NN,N4,N5,N4B,N5B
!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO 4310 N=1,2
    DO 4305 NN=1,2
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE.AND.NN.EQ.1 &
          .OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
          GO TO 4305
        ELSE
          N4=NX+1
          N5=NY
          N4B=NX-1
          N5B=NY
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS.AND.NN.EQ.1.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
          GO TO 4305
        ELSE
          N4=NX
          N5=NY+1
          N4B=NX
          N5B=NY-1
        ENDIF
      ENDIF
!
!     ELEVATION OF EACH PAIR OF ADJACENT GRID CELLS
!
!     XVOLT,XVOLW=excess water+ice,water in destination grid cell
!     ALT1,ALT2=elevation of source,destination
!     QRQ1=equilibrium runoff
!     QR1,HQR1=runoff, convective heat from runoff
!     QR,HQR=hourly-accumulated runoff, convective heat from runoff
!
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        ALT1=ALTG(N2,N1)+XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
!
!     EAST OR SOUTH RUNOFF
!
        IF(NN.EQ.1)THEN
          ALT2=ALTG(N5,N4)+XVOLT(N5,N4)/AREA(3,NU(N5,N4),N5,N4)
          IF(ALT1.GT.ALT2)THEN
            QRQ1=AMAX1(0.0,((ALT1-ALT2)*AREA(3,NUM(N2,N1),N2,N1) &
              *AREA(3,NU(N5,N4),N5,N4)-XVOLT(N5,N4)*AREA(3,NUM(N2,N1),N2,N1) &
              +XVOLT(N2,N1)*AREA(3,NU(N5,N4),N5,N4)) &
              /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5,N4),N5,N4)))
            QR1(N,2,N5,N4)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
            HQR1(N,2,N5,N4)=cpw*TK1(0,N2,N1)*QR1(N,2,N5,N4)
            QR(N,2,N5,N4)=QR(N,2,N5,N4)+QR1(N,2,N5,N4)
            HQR(N,2,N5,N4)=HQR(N,2,N5,N4)+HQR1(N,2,N5,N4)
            QRMN(M,N,2,N5,N4)=QR1(N,2,N5,N4)
            IFLBM(M,N,2,N5,N4)=0
!     IF(NX.EQ.1.AND.NY.EQ.1)THEN
!     WRITE(*,5555)'QRFOR',I,J,M,N1,N2,N4,N5,N,NN
!    2,QRM(M,N2,N1),QR1(N,2,N5,N4),QR(N,2,N5,N4)
!    2,ALT1,ALT2,ALTG(N2,N1),ALTG(N5,N4),ALT(N2,N1),ALT(N5,N4)
!    3,QRQ1,FSLOPE(N,N2,N1),QR1(2,2,4,1)
!5555  FORMAT(A8,9I4,30E12.4)
!     ENDIF
          ELSE
            QR1(N,2,N5,N4)=0.0
            HQR1(N,2,N5,N4)=0.0
            QRMN(M,N,2,N5,N4)=0.0
            IFLBM(M,N,2,N5,N4)=1
          ENDIF
        ENDIF
!
!     WEST OR NORTH RUNOFF
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            ALTB=ALTG(N5B,N4B)+XVOLT(N5B,N4B)/AREA(3,NU(N5,N4B),N5B,N4B)
            IF(ALT1.GT.ALTB)THEN
              QRQ1=AMAX1(0.0,((ALT1-ALTB)*AREA(3,NUM(N2,N1),N2,N1) &
                *AREA(3,NU(N5B,N4B),N5B,N4B)-XVOLT(N5B,N4B) &
                *AREA(3,NUM(N2,N1),N2,N1) &
                +XVOLT(N2,N1)*AREA(3,NU(N5B,N4B),N5B,N4B)) &
                /(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5B,N4B),N5B,N4B)))
              QR1(N,1,N5B,N4B)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
              HQR1(N,1,N5B,N4B)=cpw*TK1(0,N2,N1)*QR1(N,1,N5B,N4B)
              QR(N,1,N5B,N4B)=QR(N,1,N5B,N4B)+QR1(N,1,N5B,N4B)
              HQR(N,1,N5B,N4B)=HQR(N,1,N5B,N4B)+HQR1(N,1,N5B,N4B)
              QRMN(M,N,1,N5B,N4B)=QR1(N,1,N5B,N4B)
              IFLBM(M,N,1,N5B,N4B)=1
!     WRITE(*,5555)'QRBAK',I,J,M,N1,N2,N4B,N5B,N,NN
!    2,QRM(M,N2,N1),QR1(N,1,N5B,N4B),QR(N,1,N5B,N4B)
!    2,ALT1,ALTB,ALTG(N2,N1),ALTG(N5B,N4B),QRQ1,FSLOPE(N,N2,N1)
            ELSE
              QR1(N,1,N5B,N4B)=0.0
              HQR1(N,1,N5B,N4B)=0.0
              QRMN(M,N,1,N5B,N4B)=0.0
              IFLBM(M,N,1,N5B,N4B)=0
            ENDIF
          ENDIF
        ENDIF
      ELSE
        QR1(N,2,N5,N4)=0.0
        HQR1(N,2,N5,N4)=0.0
        QRMN(M,N,2,N5,N4)=0.0
        IFLBM(M,N,2,N5,N4)=0
        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          QR1(N,1,N5B,N4B)=0.0
          HQR1(N,1,N5B,N4B)=0.0
          QRMN(M,N,1,N5B,N4B)=0.0
          IFLBM(M,N,1,N5B,N4B)=0
        ENDIF
      ENDIF
!     WRITE(*,5557)'QRFORA',I,J,M,N1,N2,N4,N5,N,NN,IFLBM(M,N,2,N5,N4)
!    2,QRM(M,N2,N1),QRMN(M,N,2,N5,N4),QR1(N,2,N5,N4),QR(N,2,N5,N4)
!5557  FORMAT(A8,10I4,30E12.4)
!     IF(N4B.GT.0.AND.N5B.GT.0)THEN
!     WRITE(*,5557)'QRBAKA',I,J,M,N1,N2,N4B,N5B,N,NN
!    2,IFLBM(M,N,1,N5B,N4B)
!    2,QRM(M,N2,N1),QRMN(M,N,1,N5B,N4B),QR1(N,1,N5B,N4B)
!    3,QR(N,1,N5B,N4B)
!     ENDIF
!
!     SNOW REDISTRIBUTION FROM SNOWPACK
!
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     ALTS1,ALTS2=elevation of source,destination snowpack surfaces
!     SS,DIST=slope,distance between source,destination
!     QSX=transfer fraction
!     QS1,QW1,QI1=snow,water,ice transfer
!     HQS1=convective heat transfer from snow,water,ice transfer
!     VOLS0,VOLW0,VOLI0=snow,water,ice volume
!     DPTHSX=minimum snowpack depth for full cover
!     QS,QW,QI=hourly-accumulated snow,water,ice transfer
!     HQS=hourly-accumd convective heat from snow,water,ice transfer
!     QSM=snow transfer for solute flux calculation
!
      IF(NN.EQ.1)THEN
        ALTS1=ALTG(N2,N1)+DPTHS(N2,N1)
        ALTS2=ALTG(N5,N4)+DPTHS(N5,N4)
        SS=(ALTS1-ALTS2)/DIST(N,NU(N5,N4),N5,N4)
        QSX=SS/AMAX1(1.0,DIST(N,NU(N5,N4),N5,N4))*XNPH
        IF(SS.GT.0.0.AND.DPTHS(N2,N1).GT.DPTHSX)THEN
          QS1(N,N5,N4)=QSX*VOLS0(1,N2,N1)
          QW1(N,N5,N4)=QSX*VOLW0(1,N2,N1)
          QI1(N,N5,N4)=QSX*VOLI0(1,N2,N1)
          HQS1(N,N5,N4)=TK0(1,N2,N1)*(cps*QS1(N,N5,N4) &
            +cpw*QW1(N,N5,N4)+cpi*QI1(N,N5,N4))
        ELSEIF(SS.LT.0.0.AND.DPTHS(N5,N4).GT.DPTHSX)THEN
          QS1(N,N5,N4)=QSX*VOLS0(1,N5,N4)
          QW1(N,N5,N4)=QSX*VOLW0(1,N5,N4)
          QI1(N,N5,N4)=QSX*VOLI0(1,N5,N4)
          HQS1(N,N5,N4)=TK0(1,N5,N4)*(cps*QS1(N,N5,N4) &
            +cpw*QW1(N,N5,N4)+cpi*QI1(N,N5,N4))
        ELSE
          QS1(N,N5,N4)=0.0_r8
          QW1(N,N5,N4)=0.0_r8
          QI1(N,N5,N4)=0.0_r8
          HQS1(N,N5,N4)=0.0_r8
        ENDIF
        QS(N,N5,N4)=QS(N,N5,N4)+QS1(N,N5,N4)
        QW(N,N5,N4)=QW(N,N5,N4)+QW1(N,N5,N4)
        QI(N,N5,N4)=QI(N,N5,N4)+QI1(N,N5,N4)
        HQS(N,N5,N4)=HQS(N,N5,N4)+HQS1(N,N5,N4)
        QSM(M,N,N5,N4)=QS1(N,N5,N4)
!     IF(NX.EQ.2.AND.NY.EQ.5)THEN
!     WRITE(*,5556)'QS1',I,J,M,N1,N2,N4,N5,N,QSX,QS1(N,N5,N4)
!    2,QW1(N,N5,N4),QI1(N,N5,N4),VOLS0(N3,N2,N1),VOLW0(N3,N2,N1)
!    3,VOLI0(N3,N2,N1),ALTS1,ALTS2,ALTG(N2,N1),ALTG(N5,N4)
!    4,DIST(N,NU(N5,N4),N5,N4),SS,DLYRS0(N2,N1),DLYRS0(N5,N4)
!    5,VOLS1(N2,N1),VOLS1(N5,N4),VOLWG(N2,N1),VOLWG(N5,N4)
!5556  FORMAT(A8,8I4,30E12.4)
!     ENDIF
      ENDIF
4305  CONTINUE
4310  CONTINUE
  end subroutine LateralHydroExchange
!------------------------------------------------------------------------------------------

  subroutine Subsurface3DFlow(M,NY,NX,NHE,NVS)
  implicit none
  integer, intent(in) :: M,NY,NX,NHE,NVS

  integer :: N,N1,N2,N3,N4,N5,N6,L,LL,K1,KL
  !     begin_execution
  !
  !     WATER AND ENERGY TRANSFER THROUGH SOIL PROFILE
  !
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !
  IFLGH=0
  DO 4400 L=1,NL(NY,NX)
    N1=NX
    N2=NY
    N3=L
    !
    !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
    !
    DO 4320 N=NCN(N2,N1),3
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE)THEN
          GO TO 4320
        ELSE
          N4=NX+1
          N5=NY
          N6=L
    !
          !     ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
          !
          !     IF(N2.EQ.2.AND.(N1.EQ.2.OR.N1.EQ.3).AND.L.LE.15)THEN
          !     GO TO 4320
          !     ENDIF
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS)THEN
          GO TO 4320
        ELSE
          N4=NX
          N5=NY+1
          N6=L
    !
          !     ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
          !
          !     IF(N1.EQ.3.AND.(N2.EQ.1.OR.N2.EQ.2).AND.L.LE.15)THEN
          !     GO TO 4320
          !     ENDIF
          !
          !     END ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
          !
        ENDIF
      ELSEIF(N.EQ.3)THEN
        IF(L.EQ.NL(NY,NX))THEN
          GO TO 4320
        ELSE
          N4=NX
          N5=NY
          N6=L+1
        ENDIF
      ENDIF
!
!     SKIP NON-EXISTENT DESTINATION SOIL LAYERS
!
      DO 1100 LL=N6,NL(NY,NX)
        IF(VOLX(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
          N6=LL
          GO TO 1101
        ENDIF
1100  CONTINUE
1101  CONTINUE
      IF(N3.EQ.NU(N2,N1))N6X(N2,N1)=N6
      !
      !     POROSITIES 'THETP*', WATER CONTENTS 'THETA*', AND POTENTIALS
      !     'PSIS*' FOR EACH GRID CELL
      !
      !     THETA1,THETAL=micropore water concn in source,destination cells
      !     THETY=hygroscopic water concentration
      !     POROS=soil porosity
      !     VOLXI=soil volume excluding rock, macropore
      !
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
        IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4) &
          .AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
          THETA1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1) &
            ,safe_adb(VOLW1(N3,N2,N1),VOLY(N3,N2,N1))))
          THETAL=AMAX1(THETY(N6,N5,N4),AMIN1(POROS(N6,N5,N4) &
            ,safe_adb(VOLW1(N6,N5,N4),VOLY(N6,N5,N4))))
          !
          !     WATER POTENTIAL OF UPPER LAYER
          !
          !     BKVL=soil mass
          !     FC,WP=water contents at field capacity,wilting point
          !     FCL,WPL=log FC,WP
          !     FCD,PSD=FCL-WPL,log(POROS)-FCL
          !     PSISA1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
          !     PSIMX,PSIMD,PSIMS=log water potential at FC,WP,saturation
          !     PSISD=PSIMX-PSIMS
          !     SRP=parameter for deviation from linear log-log water retention
          !     PSISO=osmotic potential
          !
          IF(BKVL(N3,N2,N1).GT.ZEROS(NY,NX))THEN
            IF(THETA1.LT.FC(N3,N2,N1))THEN
              PSISA1(N3,N2,N1)=AMAX1(PSIHY,-EXP(PSIMX(N2,N1) &
                +((FCL(N3,N2,N1)-LOG(THETA1)) &
                /FCD(N3,N2,N1)*PSIMD(N2,N1))))
            ELSEIF(THETA1.LT.POROS(N3,N2,N1)-DTHETW)THEN
              PSISA1(N3,N2,N1)=-EXP(PSIMS(N2,N1) &
                +(((PSL(N3,N2,N1)-LOG(THETA1)) &
                /PSD(N3,N2,N1))**SRP(N3,N2,N1)*PSISD(N2,N1)))
            ELSE
              PSISA1(N3,N2,N1)=PSISE(N3,N2,N1)
            ENDIF
            !
            !     SUBSURFCE UPPER WATER LAYER
            !
            !     THETIX,THETWX=ice,water concentration
            !     FCI,WPI=ice field capacity,wilting point
            !     PSISA1=matric water potential
            !
          ELSEIF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1).and.&
            THETIX(N3,N2,N1)>ZEROS2(N2,N1))THEN
            FCX=FCI*THETIX(N3,N2,N1)
            WPX=WPI*THETIX(N3,N2,N1)
            FCLX=LOG(FCX)
            WPLX=LOG(WPX)
            PSDX=PSL(N3,N2,N1)-FCLX
            FCDX=FCLX-WPLX
            IF(THETWX(N3,N2,N1).LT.FCX)THEN
              PSISA1(N3,N2,N1)=AMAX1(PSIHY,-EXP(PSIMX(N2,N1) &
                +((FCLX-LOG(THETWX(N3,N2,N1))) &
                /FCDX*PSIMD(NY,NX))))
            ELSEIF(THETWX(N3,N2,N1).LT.POROS(N3,N2,N1)-DTHETW)THEN
              PSISA1(N3,N2,N1)=-EXP(PSIMS(N2,N1) &
                +(((PSL(N3,N2,N1)-LOG(THETWX(N3,N2,N1))) &
                /PSDX)*PSISD(N2,N1)))
            ELSE
              PSISA1(N3,N2,N1)=PSISE(N3,N2,N1)
            ENDIF
          ELSE
            PSISA1(N3,N2,N1)=PSISE(N3,N2,N1)
          ENDIF
          !     IF(N1.EQ.4.AND.N2.EQ.1.AND.N3.EQ.11)THEN
          !     WRITE(*,1119)'PSISA1',I,J,M,N,N1,N2,N3,PSISA1(N3,N2,N1)
          !    2,THETWX(N3,N2,N1),THETIX(N3,N2,N1),FCX,WPX
          !    3,BKVL(N3,N2,N1),VOLX(N3,N2,N1)
          !     ENDIF
          !
          !     WATER POTENTIAL OF LOWER LAYER
!
          IF(BKVL(N6,N5,N4).GT.ZEROS(NY,NX))THEN
            IF(THETAL.LT.FC(N6,N5,N4))THEN
              PSISA1(N6,N5,N4)=AMAX1(PSIHY,-EXP(PSIMX(N5,N4) &
                +((FCL(N6,N5,N4)-LOG(THETAL)) &
                /FCD(N6,N5,N4)*PSIMD(N5,N4))))
            ELSEIF(THETAL.LT.POROS(N6,N5,N4)-DTHETW)THEN
              PSISA1(N6,N5,N4)=-EXP(PSIMS(N5,N4) &
                +(((PSL(N6,N5,N4)-LOG(THETAL)) &
                /PSD(N6,N5,N4))**SRP(N6,N5,N4)*PSISD(N5,N4)))
            ELSE
              PSISA1(N6,N5,N4)=PSISE(N6,N5,N4)
            ENDIF
            !
            !     SUBSURFCE LOWER WATER LAYER
            !
          ELSEIF(VOLX(N6,N5,N4).GT.ZEROS2(N5,N4).and.&
            THETIX(N6,N5,N4)>ZEROS2(N5,N4))THEN
            FCX=FCI*THETIX(N6,N5,N4)
            WPX=WPI*THETIX(N6,N5,N4)
            FCLX=LOG(FCX)
            WPLX=LOG(WPX)
            PSDX=PSL(N6,N5,N4)-FCLX
            FCDX=FCLX-WPLX
            IF(THETWX(N6,N5,N4).LT.FCX)THEN
              PSISA1(N6,N5,N4)=AMAX1(PSIHY,-EXP(PSIMX(N5,N4) &
                +((FCLX-LOG(THETWX(N6,N5,N4))) &
                /FCDX*PSIMD(NY,NX))))
            ELSEIF(THETWX(N6,N5,N4).LT.POROS(N6,N5,N4)-DTHETW)THEN
              PSISA1(N6,N5,N4)=-EXP(PSIMS(NY,NX) &
                +(((PSL(N6,N5,N4)-LOG(THETWX(N6,N5,N4))) &
                /PSDX)*PSISD(NY,NX)))
            ELSE
              PSISA1(N6,N5,N4)=PSISE(N6,N5,N4)
            ENDIF
          ELSE
            PSISA1(N6,N5,N4)=PSISE(N6,N5,N4)
          ENDIF
          !     IF(N1.EQ.4.AND.N2.EQ.1.AND.N3.EQ.11)THEN
          !     WRITE(*,1119)'PSISAL',I,J,M,N,N4,N5,N6,PSISA1(N6,N5,N4)
          !    2,THETWX(N6,N5,N4),THETIX(N6,N5,N4),FCX,WPX
          !    3,BKVL(N6,N5,N4),VOLX(N6,N5,N4)
          !     ENDIF
          !
          !     ACCOUNT FOR WETTING FRONTS WHEN CALCULATING WATER CONTENTS,
          !     MATRIC WATER POTENTIALS AND HYDRAULIC CONDUCTIVITIES USED
          !     IN WATER FLUX CALCULATIONS
          !
          !     THETW1,THETWL=water concentrations in source,destination cells
          !     CND1,CNDL=hydraulic conductivities in source,destination cells
          !     FKSAT=reduction in soil surface Ksat from rainfall energy impact
          !     PSISM1=soil matric potential
          !     VOLWX1=VOLW1 accounting for wetting front
          !
          !     DARCY FLOW IF BOTH CELLS ARE SATURATED
          !     (CURRENT WATER POTENTIAL > AIR ENTRY WATER POTENTIAL)
          !
          IF(PSISA1(N3,N2,N1).GT.PSISA(N3,N2,N1) &
            .AND.PSISA1(N6,N5,N4).GT.PSISA(N6,N5,N4))THEN
            THETW1=THETA1
            THETWL=THETAL
            K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1) &
              -THETW1)/POROS(N3,N2,N1))+1))
            KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4) &
              -THETWL)/POROS(N6,N5,N4))+1))
            PSISM1(N3,N2,N1)=PSISA1(N3,N2,N1)
            PSISM1(N6,N5,N4)=PSISA1(N6,N5,N4)
            !
            !     GREEN-AMPT FLOW IF ONE LAYER IS SATURATED
            !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POENTIAL)
            !
            !     GREEN-AMPT FLOW IF SOURCE CELL SATURATED
            !
          ELSEIF(PSISA1(N3,N2,N1).GT.PSISA(N3,N2,N1))THEN
            THETW1=THETA1
            THETWL=AMAX1(THETY(N6,N5,N4),AMIN1(POROS(N6,N5,N4) &
              ,safe_adb(VOLWX1(N6,N5,N4),VOLY(N6,N5,N4))))
            K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1) &
              -THETW1)/POROS(N3,N2,N1))+1))
            KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4) &
              -AMIN1(THETS(N6,N5,N4),THETWL))/POROS(N6,N5,N4))+1))
            PSISM1(N3,N2,N1)=PSISA1(N3,N2,N1)
            IF(BKVL(N6,N5,N4).GT.ZEROS(NY,NX))THEN
              IF(THETWL.LT.FC(N6,N5,N4))THEN
                PSISM1(N6,N5,N4)=AMAX1(PSIHY,-EXP(PSIMX(N5,N4) &
                  +((FCL(N6,N5,N4)-LOG(THETWL)) &
                  /FCD(N6,N5,N4)*PSIMD(N5,N4))))
              ELSEIF(THETWL.LT.POROS(N6,N5,N4)-DTHETW)THEN
                PSISM1(N6,N5,N4)=-EXP(PSIMS(N5,N4) &
                  +(((PSL(N6,N5,N4)-LOG(THETWL)) &
                  /PSD(N6,N5,N4))**SRP(N6,N5,N4)*PSISD(N5,N4)))
              ELSE
                THETWL=POROS(N6,N5,N4)
                PSISM1(N6,N5,N4)=PSISE(N6,N5,N4)
              ENDIF
            ELSE
              THETWL=POROS(N6,N5,N4)
              PSISM1(N6,N5,N4)=PSISE(N6,N5,N4)
            ENDIF
            !     IF(N3.EQ.NUM(NY,NX))THEN
            !     WRITE(*,1116)'GA',I,J,M,N1,N2,N3,N4,N5,N6,N
            !    3,PSISM1(N3,N2,N1),PSISM1(N6,N5,N4)
            !    3,VOLW1(N3,N2,N1),VOLW1(N6,N5,N4)
            !    3,VOLWX1(N3,N2,N1),VOLWX1(N6,N5,N4)
            !    6,THETW1,THETWL
!1116  FORMAT(A8,10I4,100E12.4)
            !     ENDIF
            !
            !     GREEN-AMPT FLOW IF ADJACENT CELL SATURATED
!
          ELSEIF(PSISA1(N6,N5,N4).GT.PSISA(N6,N5,N4))THEN
            THETW1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1) &
              ,safe_adb(VOLWX1(N3,N2,N1),VOLY(N3,N2,N1))))
            THETWL=THETAL
            K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1) &
              -AMIN1(THETS(N3,N2,N1),THETW1))/POROS(N3,N2,N1))+1))
            KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4) &
              -THETWL)/POROS(N6,N5,N4))+1))
            IF(BKVL(N3,N2,N1).GT.ZEROS(NY,NX))THEN
              IF(THETW1.LT.FC(N3,N2,N1))THEN
                PSISM1(N3,N2,N1)=AMAX1(PSIHY,-EXP(PSIMX(N2,N1) &
                  +((FCL(N3,N2,N1)-LOG(THETW1)) &
                  /FCD(N3,N2,N1)*PSIMD(N2,N1))))
              ELSEIF(THETW1.LT.POROS(N3,N2,N1)-DTHETW)THEN
                PSISM1(N3,N2,N1)=-EXP(PSIMS(N2,N1) &
                  +(((PSL(N3,N2,N1)-LOG(THETW1)) &
                  /PSD(N3,N2,N1))**SRP(N3,N2,N1)*PSISD(N2,N1)))
              ELSE
                THETW1=POROS(N3,N2,N1)
                PSISM1(N3,N2,N1)=PSISE(N3,N2,N1)
              ENDIF
            ELSE
              THETW1=POROS(N3,N2,N1)
              PSISM1(N3,N2,N1)=PSISE(N3,N2,N1)
            ENDIF
            !
            !     RICHARDS FLOW IF NEITHER CELL IS SATURATED
            !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POTENTIAL)
            !
          ELSE
            THETW1=THETA1
            THETWL=THETAL
            K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1) &
              -THETW1)/POROS(N3,N2,N1))+1))
            KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4) &
              -THETWL)/POROS(N6,N5,N4))+1))
            PSISM1(N3,N2,N1)=PSISA1(N3,N2,N1)
            PSISM1(N6,N5,N4)=PSISA1(N6,N5,N4)
          ENDIF
          !
          !     HYDRAULIC CONUCTIVITY
          !
          !     CND1,CNDL=hydraulic conductivity of source,destination layer
          !     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
          !
          IF(N3.EQ.NUM(NY,NX))THEN
            CND1=HCND(N,K1,N3,N2,N1)*FKSAT
          ELSE
            CND1=HCND(N,K1,N3,N2,N1)
          ENDIF
          CNDL=HCND(N,KL,N6,N5,N4)
          !
          !     TOTAL SOIL WATER POTENTIAL = MATRIC, GRAVIMETRIC + OSMOTIC
          !
          !     PSISM1,PSISH,PSISO=soil matric,gravitational,osmotic potentials
    !
          PSIST1=PSISM1(N3,N2,N1)+PSISH(N3,N2,N1)+PSISO(N3,N2,N1)
          PSISTL=PSISM1(N6,N5,N4)+PSISH(N6,N5,N4)+PSISO(N6,N5,N4)
          PSISV1=PSISM1(N3,N2,N1)+PSISO(N3,N2,N1)
          PSISVL=PSISM1(N6,N5,N4)+PSISO(N6,N5,N4)
          !     IF(N6.EQ.12)THEN
          !     WRITE(*,7272)'PSIM',I,J,M,N1,N2,N3,N4,N5,N6
          !    2,PSISM1(N3,N2,N1),PSISM1(N6,N5,N4),PSISA1(N3,N2,N1)
          !    3,PSISA1(N6,N5,N4),THETWL,THETAL
          !    2,PSIMX(N5,N4),FCL(N6,N5,N4),FCD(N6,N5,N4),PSIMD(N5,N4)
          !    3,POROS(N6,N5,N4),PSIMS(N5,N4),PSL(N6,N5,N4),PSD(N6,N5,N4)
          !    4,SRP(N6,N5,N4),PSISD(N5,N4),PSISE(N6,N5,N4)
          !    5,THETY(N6,N5,N4),POROS(N6,N5,N4),VOLW1(N6,N5,N4),VOLX(N6,N5,N4)
!7272  FORMAT(A8,9I4,30E12.4)
          !     ENDIF
          !
          !     HYDRAULIC CONDUCTIVITY FROM CURRENT WATER CONTENT
          !     AND LOOKUP ARRAY GENERATED IN 'HOUR1'
          !
          !     CND1,CNDL=hydraulic conductivities in source,destination cells
          !     FKSAT=reduction in soil surface Ksat from rainfall energy impact
          !     AVCNDL=source-destination hydraulic conductance
          !     DLYR=layer thickness
          !
          IF(CND1.GT.ZERO.AND.CNDL.GT.ZERO)THEN
            AVCNDL=2.0*CND1*CNDL/(CND1*DLYR(N,N6,N5,N4) &
              +CNDL*DLYR(N,N3,N2,N1))
          ELSE
            AVCNDL=0.0
          ENDIF
          !
          !     WATER FLUX FROM WATER POTENTIALS, HYDRAULIC CONDUCTIVITY
          !     CONSTRAINED BY WATER POTENTIAL GRADIENT, COUPLED WITH
          !     CONVECTIVE HEAT FLUX FROM WATER FLUX
          !
          !     FLQX,FLQL=micropore water flux unlimited,limited by source water
          !     XNPH=time step of flux calculations
          !     VOLW2,VOLP1=water,air contents of source,destination micropores
          !     HWFLWL=convective heat flux from micropore water flux
          !     VOLP1Z=excess water+ice relative to porosity
          !
          FLQX=AVCNDL*(PSIST1-PSISTL)*AREA(N,N3,N2,N1)*XNPH
          IF(FLQX.GE.0.0)THEN
            IF(THETW1.GT.THETS(N3,N2,N1))THEN
              FLQZ=FLQX+AMIN1((THETW1-THETS(N3,N2,N1)) &
                *VOLY(N3,N2,N1),AMAX1(0.0,(THETS(N6,N5,N4)-THETWL) &
                *VOLY(N6,N5,N4)))*XNPX
            ELSE
              FLQZ=FLQX
            ENDIF
            FLQL=AMAX1(0.0,AMIN1(FLQZ,VOLW2(N3,N2,N1)*XNPX &
              ,VOLP1(N6,N5,N4)*XNPX))
            FLQ2=AMAX1(0.0,AMIN1(FLQX,VOLW2(N3,N2,N1)*XNPX &
              ,VOLP1(N6,N5,N4)*XNPX))
            !     FLQL1=(THETW1-THETS(N3,N2,N1))*VOLY(N3,N2,N1)
            !     FLQL2=(THETS(N6,N5,N4)-THETWL)*VOLY(N6,N5,N4)
            !     FLQL3=FLQX+AMIN1(FLQL1,AMAX1(0.0,FLQL2))*XNPX
            !     FLQL4=AMAX1(0.0,AMIN1(FLQL3,VOLP1(N6,N5,N4)*XNPX))
          ELSE
            IF(THETWL.GT.THETS(N6,N5,N4))THEN
              FLQZ=FLQX+AMAX1((THETS(N6,N5,N4)-THETWL) &
                *VOLY(N6,N5,N4),AMIN1(0.0,(THETW1-THETS(N3,N2,N1)) &
                *VOLY(N3,N2,N1)))*XNPX
            ELSE
              FLQZ=FLQX
            ENDIF
            FLQL=AMIN1(0.0,AMAX1(FLQZ,-VOLW2(N6,N5,N4)*XNPX,-VOLP1(N3,N2,N1)*XNPX))
            FLQ2=AMIN1(0.0,AMAX1(FLQX,-VOLW2(N6,N5,N4)*XNPX,-VOLP1(N3,N2,N1)*XNPX))
            !     FLQL1=(THETS(N6,N5,N4)-THETWL)*VOLY(N6,N5,N4)
            !     FLQL2=(THETW1-THETS(N3,N2,N1))*VOLY(N3,N2,N1)
            !     FLQL3=FLQX+AMAX1(FLQL1,AMIN1(0.0,FLQL2))*XNPX
            !     FLQL4=AMIN1(0.0,AMAX1(FLQL3,-VOLP1(N3,N2,N1)*XNPX))
          ENDIF
          IF(N.EQ.3.AND.VOLP1Z(N6,N5,N4).LT.0.0)THEN
            FLQL=FLQL+AMIN1(0.0,AMAX1(-VOLW2(N6,N5,N4)*XNPX,VOLP1Z(N6,N5,N4)))
            FLQ2=FLQ2+AMIN1(0.0,AMAX1(-VOLW2(N6,N5,N4)*XNPX,VOLP1Z(N6,N5,N4)))
          ENDIF
          IF(FLQL.GT.0.0)THEN
            HWFLQL=cpw*TK1(N3,N2,N1)*FLQL
          ELSE
            HWFLQL=cpw*TK1(N6,N5,N4)*FLQL
          ENDIF
          VOLW2(N3,N2,N1)=VOLW2(N3,N2,N1)-FLQL
          VOLW2(N6,N5,N4)=VOLW2(N6,N5,N4)+FLQL
!
          !     MACROPORE FLOW FROM POISEUILLE FLOW IF MACROPORES PRESENT
          !
          !     PSISH1,PSISHL=macropore total water potl in source,destination
          !     DLYR=layer thickness
          !     VOLWH1,VOLPH1=macropore water,air content
    !
          IF(VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1) &
            .AND.VOLAH1(N6,N5,N4).GT.ZEROS2(N5,N4).AND.IFLGH.EQ.0)THEN
            PSISH1=PSISH(N3,N2,N1)+0.0098*DLYR(3,N3,N2,N1) &
              *(AMIN1(1.0,AMAX1(0.0,VOLWH1(N3,N2,N1)/VOLAH1(N3,N2,N1)))-0.5)
            PSISHL=PSISH(N6,N5,N4)+0.0098*DLYR(3,N6,N5,N4) &
              *(AMIN1(1.0,AMAX1(0.0,VOLWH1(N6,N5,N4)/VOLAH1(N6,N5,N4)))-0.5)
            !
            !     MACROPORE FLOW IF GRAVITATIONAL GRADIENT IS POSITIVE
            !     AND MACROPORE POROSITY EXISTS IN ADJACENT CELL
            !
            !     FLWHX,FLWHL=macropore water flux unltd,ltd by source water
            !     XNPH=time step of flux calculations
            !     VOLW2,VOLP1=water,air contents of source,destination micropores
            !     HWFLHL=convective heat flux from micropore water flux
            !
            FLWHX=AVCNHL(N,N6,N5,N4)*(PSISH1-PSISHL)*AREA(N,N3,N2,N1)*XNPH
            IF(N.NE.3)THEN
              IF(PSISH1.GT.PSISHL)THEN
                FLWHL(N,N6,N5,N4)=AMAX1(0.0,AMIN1(AMIN1(VOLWH1(N3,N2,N1) &
                  ,VOLPH1(N6,N5,N4))*XNPX,FLWHX))
              ELSEIF(PSISH1.LT.PSISHL)THEN
                FLWHL(N,N6,N5,N4)=AMIN1(0.0,AMAX1(AMAX1(-VOLWH1(N6,N5,N4) &
                  ,-VOLPH1(N3,N2,N1))*XNPX,FLWHX))
              ELSE
                FLWHL(N,N6,N5,N4)=0.0
              ENDIF
            ELSE
              FLWHL(N,N6,N5,N4)=AMAX1(0.0,AMIN1(AMIN1(VOLWH1(N3,N2,N1)*XNPX &
                +FLWHL(N,N3,N2,N1),VOLPH1(N6,N5,N4)*XNPX),FLWHX))
            ENDIF
            IF(N.EQ.3)THEN
              FLWHL(N,N6,N5,N4)=FLWHL(N,N6,N5,N4)+AMIN1(0.0,VOLPH1Z(N6,N5,N4))
            ENDIF
            FLWHM(M,N,N6,N5,N4)=FLWHL(N,N6,N5,N4)
            !     IF(N4.EQ.1)THEN
            !     WRITE(*,5478)'FLWH',I,J,M,N1,N2,N3,IFLGH
            !    2,FLHM,FLWHX,FLWHL(N,N3,N2,N1),FLWHL(N,N6,N5,N4)
            !    2,AVCNHL(N,N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
            !    3,VOLPH1(N3,N2,N1),VOLPH1(N6,N5,N4),VOLWH1(N3,N2,N1)
            !    4,VOLWH1(N6,N5,N4),VOLAH1(N3,N2,N1),VOLAH1(N6,N5,N4)
            !    5,DLYR(N,N6,N5,N4),DLYR(N,N3,N2,N1),AREA(N,N3,N2,N1)
            !    7,CNDH1(N3,N2,N1),CNDH1(N6,N5,N4),XNPH,HWFLHL
!5478  FORMAT(A8,7I4,30E12.4)
            !     ENDIF
          ELSE
            FLWHL(N,N6,N5,N4)=0.0
            FLWHM(M,N,N6,N5,N4)=0.0
            IF(VOLPH1(N6,N5,N4).LE.0.0)IFLGH=1
          ENDIF
          IF(FLWHL(N,N6,N5,N4).GT.0.0)THEN
            HWFLHL=cpw*TK1(N3,N2,N1)*FLWHL(N,N6,N5,N4)
          ELSE
            HWFLHL=cpw*TK1(N6,N5,N4)*FLWHL(N,N6,N5,N4)
          ENDIF
!
          !     VAPOR PRESSURE AND DIFFUSIVITY IN EACH GRID CELL
          !
          !     THETPM,THETX=current, minimum air-filled porosity
          !     TK11,TK12=interim soil temperature in source,destination
          !     VP1,VPL=vapor concentration in source,destination
          !     PSISV1,PSISVL=matric+osmotic water potl in source,destination
          !     CNV1,CNV2=vapor conductivities of source, destination
          !     POROS,POROQ=porosity, tortuosity
          !     WGSGL=vapor diffusivity
          !     ATCNVL=source,destination vapor conductance
          !     DLYR=soil layer depth
          !     FLVC,FLVX=vapor flux unlimited,limited by vapor
          !     VPY=equilibrium vapor concentration
          !     XNPX=time step for flux calculations
          !     FLVL,HWFLVL=vapor flux and its convective heat flux
!
          IF(THETPM(M,N3,N2,N1).GT.THETX &
            .AND.THETPM(M,N6,N5,N4).GT.THETX)THEN
            TK11=TK1(N3,N2,N1)
            TK12=TK1(N6,N5,N4)
            !VP1=2.173E-03/TK11 &
            !*0.61*EXP(5360.0*(3.661E-03-1.0/TK11)) &
            VP1=vapsat(TK11)*EXP(18.0*PSISV1/(8.3143*TK11))
            !VPL=2.173E-03/TK12 &
            !*0.61*EXP(5360.0*(3.661E-03-1.0/TK12)) &
            VPL=vapsat(TK12)*EXP(18.0*PSISVL/(8.3143*TK12))
            CNV1=WGSGL(N3,N2,N1)*THETPM(M,N3,N2,N1)*POROQ &
              *THETPM(M,N3,N2,N1)/POROS(N3,N2,N1)
            CNVL=WGSGL(N6,N5,N4)*THETPM(M,N6,N5,N4)*POROQ &
              *THETPM(M,N6,N5,N4)/POROS(N6,N5,N4)
            ATCNVL=2.0*CNV1*CNVL &
              /(CNV1*DLYR(N,N6,N5,N4)+CNVL*DLYR(N,N3,N2,N1))
            !
            !     VAPOR FLUX FROM VAPOR PRESSURE AND DIFFUSIVITY,
            !     AND CONVECTIVE HEAT FLUX FROM VAPOR FLUX
        !
            FLVC=ATCNVL*(VP1-VPL)*AREA(N,N3,N2,N1)*XNPH
            VPY=(VP1*VOLPM(M,N3,N2,N1)+VPL*VOLPM(M,N6,N5,N4)) &
              /(VOLPM(M,N3,N2,N1)+VOLPM(M,N6,N5,N4))
            FLVX=(VP1-VPY)*VOLPM(M,N3,N2,N1)*XNPX
            IF(FLVC.GE.0.0)THEN
              FLVL=AMAX1(0.0,AMIN1(FLVC,FLVX))
              HWFLVL=(cpw*TK1(N3,N2,N1)+VAP)*FLVL
            ELSE
              FLVL=AMIN1(0.0,AMAX1(FLVC,FLVX))
              HWFLVL=(cpw*TK1(N6,N5,N4)+VAP)*FLVL
            ENDIF
          ELSE
            FLVL=0.0
            HWFLVL=0.0
          ENDIF
          !
          !     FLWL=total water+vapor flux to destination
          !     FLWLX=total unsaturated water+vapor flux to destination
          !     HWFLWL=total convective heat flux from water+vapor flux
          !
          FLWL(N,N6,N5,N4)=FLQL+FLVL
          if(abs(FLWL(N,N6,N5,N4))>1.e20)then
            write(*,*)'FLQL+FLVL=',FLQL,FLVL
            write(*,*)'at line',__LINE__
            call endrun(trim(mod_filename)//'at line',__LINE__)
          endif
          FLWLX(N,N6,N5,N4)=FLQ2+FLVL
          HWFLWL=HWFLQL+HWFLVL
          !     IF(I.EQ.232.AND.N3.LE.NUM(N2,N1))THEN
          !     WRITE(*,1115)'FLWL',I,J,M,N1,N2,N3,N4,N5,N6,N,K1,KL
          !    2,FLWL(N,N3,N2,N1),FLWL(N,N6,N5,N4),FLQX,FLQZ,FLQL,FLQ2,FLVL
          !    3,PSIST1,PSISTL
          !    4,AVCNDL,CND1,CNDL,FKSAT
          !    2,VOLP1(N3,N2,N1),VOLP1(N6,N5,N4),VOLW1(N3,N2,N1)
          !    3,VOLWX1(N3,N2,N1),VOLW1(N6,N5,N4),VOLWX1(N6,N5,N4)
          !    6,THETW1,THETS(N3,N2,N1),THETWL,THETS(N6,N5,N4)
          !    7,PSISA1(N3,N2,N1),PSISA1(N6,N5,N4),PSISM1(N3,N2,N1)
          !    7,PSISM1(N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
          !    3,VOLY(N3,N2,N1),VOLY(N6,N5,N4)
          !    3,FLVL,FLVX,VP1,VPL,VPY,CNV1,CNVL,ATCNVL
          !    4,VOLA1(N6,N5,N4),VOLI1(N6,N5,N4),SCNV(N3,N2,N1)
          !    5,SCNV(N6,N5,N4),VOLP1(N3,N2,N1),VOLP1(N6,N5,N4)
          !    7,VOLP1Z(N6,N5,N4),VOLT(N3,N2,N1),VOLT(N6,N5,N4)
          !    8,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4),AREA(N,N3,N2,N1)
          !    3,HWFLWL,HWFLQL,HWFLVL
          !    9,TK1(N3,N2,N1),TK1(N6,N5,N4),TKY
          !    9,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),POROS(N6,N5,N4)
          !    9,VOLP1(N3,N2,N1),VOLX(N3,N2,N1),VOLA1(N6,N5,N4),VOLA1(N3,N2,N1)
          !    8,VOLP1(N6,N5,N4),VOLX(N6,N5,N4),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
          !    1,THETPM(M,N3,N2,N1),THETPM(M,N6,N5,N4),THETX
          !    4,FLQL1,FLQL2,FLQL3,FLQL4
!1115  FORMAT(A8,12I4,100E12.4)
          !     ENDIF
          !
          !     THERMAL CONDUCTIVITY IN EACH GRID CELL
          !
          !     DTH*,RYL*,DNU*,TRB*,XNUS*=turbulence effects on thermal conductivity
          !     THETWX,THETPX=water,air concentration
          !     TCNDW*,TCNDA*=thermal conductivity of water,air
          !     TCND1,TCNDL=soil thermal conductivity in source,destination
          !     WTHET*=multiplier for air concn in thermal conductivity
          !     ATCNDL=source-destination thermal conductance
      !
          DTKX=ABS(TK1(N3,N2,N1)-TK1(N6,N5,N4))*1.0E-06
          IF(BKDS(N3,N2,N1).GT.ZERO.OR.THETWX(N3,N2,N1) &
            +THETIX(N3,N2,N1).GT.ZERO)THEN
            DTHW1=AMAX1(0.0,THETWX(N3,N2,N1)-TRBW)**3
            DTHA1=AMAX1(0.0,THETPX(N3,N2,N1)-TRBA)**3
            RYLXW1=DTKX*DTHW1
            RYLXA1=DTKX*DTHA1
            RYLNW1=AMIN1(1.0E+04,RYLXW*RYLXW1)
            RYLNA1=AMIN1(1.0E+04,RYLXA*RYLXA1)
            XNUSW1=AMAX1(1.0,0.68+0.67*RYLNW1**0.25/DNUSW)
            XNUSA1=AMAX1(1.0,0.68+0.67*RYLNA1**0.25/DNUSA)
            TCNDW1=2.067E-03*XNUSW1
            TCNDA1=9.050E-05*XNUSA1
            WTHET1=1.467-0.467*THETPY(N3,N2,N1)
            TCND1=(STC(N3,N2,N1)+THETWX(N3,N2,N1)*TCNDW1 &
              +0.611*THETIX(N3,N2,N1)*7.844E-03 &
              +WTHET1*THETPX(N3,N2,N1)*TCNDA1) &
              /(DTC(N3,N2,N1)+THETWX(N3,N2,N1)+0.611*THETIX(N3,N2,N1) &
              +WTHET1*THETPX(N3,N2,N1))
          ELSE
            TCND1=0.0
          ENDIF
          IF(BKDS(N6,N5,N4).GT.ZERO.OR.THETWX(N6,N5,N4) &
            +THETIX(N6,N5,N4).GT.ZERO)THEN
            DTHW2=AMAX1(0.0,THETWX(N6,N5,N4)-TRBW)**3
            DTHA2=AMAX1(0.0,THETPX(N6,N5,N4)-TRBA)**3
            RYLXW2=DTKX*DTHW2
            RYLXA2=DTKX*DTHA2
            RYLNW2=AMIN1(1.0E+04,RYLXW*RYLXW2)
            RYLNA2=AMIN1(1.0E+04,RYLXA*RYLXA2)
            XNUSW2=AMAX1(1.0,0.68+0.67*RYLNW2**0.25/DNUSW)
            XNUSA2=AMAX1(1.0,0.68+0.67*RYLNA2**0.25/DNUSA)
            TCNDW2=2.067E-03*XNUSW2
            TCNDA2=9.050E-05*XNUSA2
            WTHET2=1.467-0.467*THETPY(N6,N5,N4)
            TCND2=(STC(N6,N5,N4)+THETWX(N6,N5,N4)*TCNDW2 &
              +0.611*THETIX(N6,N5,N4)*7.844E-03 &
              +WTHET2*THETPX(N6,N5,N4)*TCNDA2) &
              /(DTC(N6,N5,N4)+THETWX(N6,N5,N4)+0.611*THETIX(N6,N5,N4) &
              +WTHET2*THETPX(N6,N5,N4))
          ELSE
            TCND2=0.0
          ENDIF
          ATCNDL=(2.0*TCND1*TCND2)/(TCND1*DLYR(N,N6,N5,N4) &
            +TCND2*DLYR(N,N3,N2,N1))
          !
          !     HEAT FLOW FROM THERMAL CONDUCTIVITY AND TEMPERATURE GRADIENT
          !
          !     VHCP1,VHCPW=volumetric heat capacity of soil,snowpack
          !     TK1X,TKLX=interim temperatures of source,destination
          !     HWFLVL,HFLXG=convective heat from soil vapor flux
          !     HFLXG=storage heat flux from snowpack
          !     TKY=equilibrium source-destination temperature
          !     HFLWC,HFLWX=source-destination heat flux unltd,ltd by heat
          !     ATCNDL=source-destination thermal conductance
          !     HFLWSX=source-destination conductive heat flux
          !     HFLWL=total conductive+convective source-destination heat flux
          !
          IF(VHCP1(N3,N2,N1).GT.VHCPNX(NY,NX))THEN
            IF(N3.EQ.NUM(NY,NX).AND.VHCPW(1,N2,N1).LE.VHCPWX(N2,N1))THEN
              TK1X=TK1(N3,N2,N1)-(HWFLVL-HFLXG)/VHCP1(N3,N2,N1)
              if(abs(TK1X)>1.e5_r8)then
                write(*,*)'TK1(N3,N2,N1)-HWFLVL/VHCP1(N3,N2,N1)',&
                  TK1(N3,N2,N1),HWFLVL,HFLXG,VHCP1(N3,N2,N1)
                write(*,*)'N1,n2,n3',N1,N2,N3
                call endrun(trim(mod_filename)//' at line',__LINE__)
              endif
            ELSE
              TK1X=TK1(N3,N2,N1)-HWFLVL/VHCP1(N3,N2,N1)
              if(abs(TK1X)>1.e5_r8)then
                write(*,*)'TK1(N3,N2,N1)-HWFLVL/VHCP1(N3,N2,N1)',&
                  TK1(N3,N2,N1),HWFLVL,VHCP1(N3,N2,N1)
                write(*,*)'N1,n2,n3',N1,N2,N3
                call endrun(trim(mod_filename)//' at line',__LINE__)
              endif
            ENDIF
          ELSE
            TK1X=TK1(N3,N2,N1)
          ENDIF

          IF(VHCP1(N6,N5,N4).GT.ZEROS(NY,NX))THEN
            TKLX=TK1(N6,N5,N4)+HWFLVL/VHCP1(N6,N5,N4)
          ELSE
            TKLX=TK1(N6,N5,N4)
          ENDIF
          TKY=(VHCP1(N3,N2,N1)*TK1X+VHCP1(N6,N5,N4)*TKLX) &
            /(VHCP1(N3,N2,N1)+VHCP1(N6,N5,N4))
          HFLWX=(TK1X-TKY)*VHCP1(N3,N2,N1)*XNPX
          HFLWC=ATCNDL*(TK1X-TKLX)*AREA(N,N3,N2,N1)*XNPH
          if(HFLWX<-1.e10_r8)then
            write(*,*)'(TK1X-TKY)*VHCP1(N3,N2,N1)',TK1X,TKY,TKLX,VHCP1(N3,N2,N1)
            write(*,*)'N1,N2,N3',n1,n2,n3
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          IF(HFLWC.GE.0.0)THEN
            HFLWSX=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
          ELSE
            HFLWSX=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
            if(HFLWSX<-1.e10_r8)then
              write(*,*)'HFLWX,HFLWC',HFLWX,HFLWC
              call endrun(trim(mod_filename)//' at line',__LINE__)
            endif
          ENDIF
          HFLWL(N,N6,N5,N4)=HWFLWL+HWFLHL+HFLWSX
          if(HFLWL(N,N6,N5,N4)<-1.e10_r8)then
            write(*,*)'HWFLWL+HWFLHL+HFLWSX',HWFLWL,HWFLHL,HFLWSX
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          !     IF(N3.EQ.1)THEN
          !     WRITE(*,8765)'HFLWL',I,J,M,N1,N2,N3,N4,N5,N6,N
          !    2,HFLWL(N,N3,N2,N1),HFLWL(N,N6,N5,N4)
          !    2,HWFLWL,HWFLQL,HWFLVL,HWFLHL,HFLWC,HFLWX,HFLWSX
          !    3,ATCNDL,TK1X,TKLX,TKY,HWFLVL,TK1(N3,N2,N1),TK1(N6,N5,N4)
          !    2,TCND1,TCND2,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4)
          !    4,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),VOLY(N3,N2,N1)
          !    5,VOLY(N6,N5,N4),VOLW1(N3,N2,N1),VOLW1(N6,N5,N4)
          !    3,THETPX(N3,N2,N1),THETIX(N3,N2,N1),THETWX(N3,N2,N1)
          !    3,THETPX(N6,N5,N4),THETIX(N6,N5,N4),THETWX(N6,N5,N4)
          !    3,THETPY(N3,N2,N1),THETPY(N6,N5,N4)
          !    4,RYLNA2,XNUSA2,XNUSW2
          !    4,STC(N6,N5,N4),DTC(N6,N5,N4)
          !    2,WTHET2,TCNDA2,TCNDW2
!8765  FORMAT(A8,10I4,60E14.6)
          !     ENDIF
          !
          !     TOTAL WATER, VAPOR AND HEAT FLUXES
          !
          !     FLW,FLWX,FLWH=total water flux through micropores,macropores
          !     HFLW=total heat flux
          !     FLWM=water flux used for solute flux calculations in trnsfr.f
          !
          FLW(N,N6,N5,N4)=FLW(N,N6,N5,N4)+FLWL(N,N6,N5,N4)
          FLWX(N,N6,N5,N4)=FLWX(N,N6,N5,N4)+FLWLX(N,N6,N5,N4)
          FLWH(N,N6,N5,N4)=FLWH(N,N6,N5,N4)+FLWHL(N,N6,N5,N4)
          HFLW(N,N6,N5,N4)=HFLW(N,N6,N5,N4)+HFLWL(N,N6,N5,N4)
          if(HFLW(N,N6,N5,N4)<-1.e10_r8)then
            write(*,*)'HFLWL(N,N6,N5,N4)',HFLWL(N,N6,N5,N4)
            write(*,*)'Ns=',N,N4,n5,n6
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          FLWM(M,N,N6,N5,N4)=FLWL(N,N6,N5,N4)
          IF(N.EQ.3)THEN
            !     IF(I.EQ.55)THEN
            !     WRITE(*,1115)'FLWL2',I,J,M,N1,N2,N3,N4,N5,N6,N,FLWL(N,N3,N2,N1)
            !    2,FLWL(N,N6,N5,N4),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
            !    3,FLQL,FLVL,FLQX,FLVX,HFLWX
            !    3,CND1,CNDL,AVCNDL,AVCNVL,VP1,VPL,PSIST1,PSISTL
            !    4,UAG,VOLA1(N6,N5,N4),VOLI1(N6,N5,N4),SCNV(N6,N5,N4)
            !    5,VOLP1(N3,N2,N1),VOLP1(N6,N5,N4),TKY
            !    7,TK1(N3,N2,N1),TK1(N6,N5,N4),VOLT(N3,N2,N1),VOLT(N6,N5,N4)
            !    8,VOLW1(N6,N5,N4),VOLP1(N6,N5,N4),VOLX(N6,N5,N4),VOLW1(N3,N2,N1)
            !    9,VOLP1(N3,N2,N1),VOLX(N3,N2,N1),VOLA1(N6,N5,N4),VOLA1(N3,N2,N1)
            !    6,THETW1,THETWL,PSISA1(N3,N2,N1)
            !    7,PSISA1(N6,N5,N4),PSISM1(N3,N2,N1)
            !    7,PSISM1(N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
            !    8,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4),AREA(N,N3,N2,N1)
            !    9,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),POROS(N6,N5,N4)
            !     ENDIF
            !
            !     WATER FILM THICKNESS FOR CALCULATING GAS EXCHANGE IN TRNSFR.F
            !
            FILM(M,N6,N5,N4)=AMAX1(1.0E-06,EXP(-13.833-0.857*LOG(-PSISA1(N6,N5,N4))))
          ENDIF
        ELSEIF(N.NE.3)THEN
          FLWL(N,N6,N5,N4)=0.0
          FLWLX(N,N6,N5,N4)=0.0
          FLWHL(N,N6,N5,N4)=0.0
          HFLWL(N,N6,N5,N4)=0.0
          FLWM(M,N,N6,N5,N4)=0.0
          FLWHM(M,N,N6,N5,N4)=0.0
        ENDIF
      ELSE
        IF(N.EQ.3)THEN
          FLWL(N,N3,N2,N1)=0.0
          FLWLX(N,N3,N2,N1)=0.0
          FLWHL(N,N3,N2,N1)=0.0
          HFLWL(N,N3,N2,N1)=0.0
          FLWHM(M,N,N3,N2,N1)=0.0
          FLWHM(M,N,N3,N2,N1)=0.0
        ELSE
          FLWL(N,N6,N5,N4)=0.0
          FLWLX(N,N6,N5,N4)=0.0
          FLWHL(N,N6,N5,N4)=0.0
          HFLWL(N,N6,N5,N4)=0.0
          FLWM(M,N,N6,N5,N4)=0.0
          FLWHM(M,N,N6,N5,N4)=0.0
        ENDIF
        !     IF(I.EQ.336)THEN
        !     WRITE(*,1115)'FLWLX',I,J,M,N1,N2,N3,N4,N5,N6,N
        !    2,FLWL(N,N3,N2,N1),FLW(N,N3,N2,N1)
        !    2,FLWL(N,N6,N5,N4),FLW(N,N6,N5,N4)
        !    3,VOLX(N3,N2,N1),VOLX(N6,N5,N4)
        !     ENDIF
      ENDIF
4320  CONTINUE
4400  CONTINUE
  end subroutine Subsurface3DFlow
!------------------------------------------------------------------------------------------

  subroutine WaterHeatExchThruBoundaryFlow(M,NHW,NHE,NVN,NVS)
!     boundary flow involes exchange with external water table, and through
!     tile drainage
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS

  integer :: NY,NX
  integer :: L,LL
  integer :: N,NN,N1,N2,N3,N4,N5,N4B,N5B,N6
  integer :: M1,M2,M3,M4,M5,M6,K1,KL
!     begin_execution

  DO 9595 NX=NHW,NHE
    DO 9590 NY=NVN,NVS

      DO 9585 L=NUM(NY,NX),NL(NY,NX)
        VOLP2=VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
        VOLPX2=VOLP2
        VOLPH2=VOLAH1(L,NY,NX)-VOLWH1(L,NY,NX)-VOLIH1(L,NY,NX)
!
!     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO WATER TABLE
!
!     IDTBL=water table flag
!     DPTH,DTBLX=depth to layer midpoint,natural water table
!     PSISM1,PSISE=matric,air entry water potential
!     DTBLXX=equilibrium water potential with natural water table
!     DPTHA=active layer depth
!     IFLGU=micropore discharge flag to natural water table
!
        IF(IDTBL(NY,NX).NE.0.AND.DPTH(L,NY,NX).LT.DTBLX(NY,NX))THEN
          IF(PSISM1(L,NY,NX).GT.0.0098*(DPTH(L,NY,NX)-DTBLX(NY,NX)))THEN
            IFLGU=0
            DO 9565 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
              DTBLXX=DTBLX(NY,NX)+PSISE(LL,NY,NX)/0.0098
              IF(DPTH(LL,NY,NX).LT.DTBLXX)THEN
                IF((PSISM1(LL,NY,NX).LE.0.0098*(DPTH(LL,NY,NX)-DTBLXX) &
                  .AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
                  IFLGU=1
                ENDIF
              ENDIF
9565        CONTINUE
          ELSE
            IFLGU=1
          ENDIF
        ELSE
          IFLGU=1
        ENDIF
!
!     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO WATER TABLE
!
!     VOLAH1,VOLWH1,VOLIH1=macropore volume,water,ice content
!     DPTHH depth to layer macropore water
!     CDPTH=depth to layer bottom
!     DLYR=layer thickness
!     IFLGUH=macropore discharge flag to natural water table
!
        IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          DPTHH=CDPTH(L,NY,NX)-(VOLWH1(L,NY,NX)+VOLIH1(L,NY,NX)) &
            /VOLAH1(L,NY,NX)*DLYR(3,L,NY,NX)
        ELSE
          DPTHH=CDPTH(L,NY,NX)
        ENDIF
        IF(IDTBL(NY,NX).NE.0.AND.DPTHH.LT.DTBLX(NY,NX) &
          .AND.VOLWH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          IFLGUH=0
    !     DO 9566 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
    !     IF(DPTH(LL,NY,NX).LT.DTBLX(NY,NX))THEN
    !     IF(VOLAH1(LL,NY,NX).LE.ZEROS(NY,NX))THEN
    !     IFLGUH=1
    !     ENDIF
    !     ENDIF
!9566  CONTINUE
        ELSE
          IFLGUH=1
        ENDIF
!     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
!     WRITE(*,9567)'IFLGU',I,J,M,NX,NY,L,IFLGU,IFLGUH,PSISM1(L,NY,NX)
!    2,0.0098*(DPTH(L,NY,NX)-DTBLX(NY,NX)),DTBLXX
!    2,DPTH(L,NY,NX),DTBLX(NY,NX),DTBLZ(NY,NX),DTBLI(NY,NX)
!    3,VOLAH1(L,NY,NX),VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),CDPTH(L,NY,NX)
!    4,DLYR(3,L,NY,NX),DTBLZ(NY,NX),DPTHH,THETX,DPTHA(NY,NX)
!9567  FORMAT(A8,8I4,30E12.4)
!     ENDIF
!
!     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO TILE DRAIN
!
!     IDTBL=water table flag
!     DPTH,DTBLY=depth to layer midpoint, artificial water table
!     PSISM1,PSISE=soil,air entry matric potential
!     DTBLYX=equilibrium water potential with artificial water table
!     IFLGD=micropore discharge flag to artificial water table
!
        IF(IDTBL(NY,NX).GE.3.AND.DPTH(L,NY,NX).LT.DTBLY(NY,NX))THEN
          IF(PSISM1(L,NY,NX).GT.0.0098*(DPTH(L,NY,NX)-DTBLY(NY,NX)))THEN
            IFLGD=0
            IF(L.LT.NL(NY,NX))THEN
              DO 9568 LL=L+1,NL(NY,NX)
                DTBLYX=DTBLY(NY,NX)+PSISE(LL,NY,NX)/0.0098
                IF(DPTH(LL,NY,NX).LT.DTBLYX)THEN
                  IF((PSISM1(LL,NY,NX).LE.0.0098*(DPTH(LL,NY,NX)-DTBLYX) &
                    .AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
                    IFLGD=1
                  ENDIF
                ENDIF
9568          CONTINUE
            ENDIF
          ELSE
            IFLGD=1
          ENDIF
        ELSE
          IFLGD=1
        ENDIF
!
!     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO TILE DRAIN
!
!     VOLAH1,VOLWH1,VOLIH1=macropore volume,water,ice content
!     CDPTH=depth to layer bottom
!     DLYR=layer thickness
!     IFLGDH=macropore discharge flag to artificial water table
!
        IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          DPTHH=CDPTH(L,NY,NX)-(VOLWH1(L,NY,NX)+VOLIH1(L,NY,NX)) &
            /VOLAH1(L,NY,NX)*DLYR(3,L,NY,NX)
        ELSE
          DPTHH=CDPTH(L,NY,NX)
        ENDIF
        IF(IDTBL(NY,NX).GE.3.AND.DPTHH.LT.DTBLY(NY,NX) &
          .AND.VOLWH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          IFLGDH=0
          IF(L.LT.NL(NY,NX))THEN
            DO 9569 LL=L+1,NL(NY,NX)
              IF(DPTH(LL,NY,NX).LT.DTBLY(NY,NX))THEN
                IF(VOLAH1(LL,NY,NX).LE.ZEROS(NY,NX))THEN
                  IFLGDH=1
                ENDIF
              ENDIF
9569        CONTINUE
          ENDIF
        ELSE
          IFLGDH=1
        ENDIF
    !     IF(L.EQ.12)THEN
!     WRITE(*,9567)'IFLGD',I,J,M,NX,NY,L,IFLGD,IFLGDH,PSISM1(L,NY,NX)
!    2,PSISE(L,NY,NX)+0.0098*(DPTH(L,NY,NX)-DTBLY(NY,NX))
!    2,DPTH(L,NY,NX),DTBLY(NY,NX),DTBLD(NY,NX),DTBLDI(NY,NX)
!    3,VOLAH1(L,NY,NX),VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),CDPTH(L,NY,NX)
!    4,DLYR(3,L,NY,NX),DPTHH,THETX,DPTHA(NY,NX),RCHGFT
!     ENDIF
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
!     N3,N2,N1=L,NY,NX of source grid cell
!     M6,M5,M4=L,NY,NX of destination grid cell
!
        N1=NX
        N2=NY
        N3=L
!
!     LOCATE EXTERNAL BOUNDARIES
!
        DO 9580 N=1,3
          DO 9575 NN=1,2
            IF(N.EQ.1)THEN
              N4=NX+1
              N5=NY
              N4B=NX-1
              N5B=NY
              N6=L
              IF(NN.EQ.1)THEN
                IF(NX.EQ.NHE)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX+1
                  M5=NY
                  M6=L
                  XN=-1.0
                  RCHQF=RCHQE(M2,M1)
                  RCHGFU=RCHGEU(M2,M1)
                  RCHGFT=RCHGET(M2,M1)
                ELSE
                  GO TO 9575
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                IF(NX.EQ.NHW)THEN
                  M1=NX+1
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0
                  RCHQF=RCHQW(M5,M4)
                  RCHGFU=RCHGWU(M5,M4)
                  RCHGFT=RCHGWT(M5,M4)
                ELSE
                  GO TO 9575
                ENDIF
              ENDIF
            ELSEIF(N.EQ.2)THEN
              N4=NX
              N5=NY+1
              N4B=NX
              N5B=NY-1
              N6=L
              IF(NN.EQ.1)THEN
                IF(NY.EQ.NVS)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY+1
                  M6=L
                  XN=-1.0
                  RCHQF=RCHQS(M2,M1)
                  RCHGFU=RCHGSU(M2,M1)
                  RCHGFT=RCHGST(M2,M1)
                ELSE
                  GO TO 9575
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY+1
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0
                  RCHQF=RCHQN(M5,M4)
                  RCHGFU=RCHGNU(M5,M4)
                  RCHGFT=RCHGNT(M5,M4)
                ELSE
                  GO TO 9575
                ENDIF
              ENDIF
            ELSEIF(N.EQ.3)THEN
              N4=NX
              N5=NY
              N6=L+1
              IF(NN.EQ.1)THEN
                IF(L.EQ.NL(NY,NX))THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L+1
                  XN=-1.0
                  RCHGFU=RCHGD(M2,M1)
                  RCHGFT=1.0
                ELSE
                  GO TO 9575
                ENDIF
              ELSEIF(NN.EQ.2)THEN
                GO TO 9575
              ENDIF
            ENDIF
!
!     BOUNDARY SURFACE RUNOFF DEPENDING ON ASPECT, SLOPE
!     VELOCITY, HYDRAULIC RADIUS AND SURFACE WATER STORAGE
!
!     CDPTH,CDPTHI=current,initial surface elevation
!     BKDS=bulk density
!     IRCHG,RCHQ*=runoff boundary flags
!
            IF(L.EQ.NUM(N2,N1).AND.N.NE.3 &
              .AND.(CDPTH(NU(N2,N1)-1,N2,N1).LE.CDPTHI(N2,N1) &
              .OR.BKDS(NUI(N2,N1),N2,N1).GT.ZERO))THEN
!     WRITE(*,7744)'QRCHK',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
!    1,IRCHG(NN,N,N2,N1),RCHQF,QRM(M,N2,N1)
              IF(IRCHG(NN,N,N2,N1).EQ.0.OR.test_aeqb(RCHQF,0._r8) &
                .OR.ABS(QRM(M,N2,N1)).LT.ZEROS(N2,N1))THEN
                QR1(N,NN,M5,M4)=0.0
                HQR1(N,NN,M5,M4)=0.0
              ELSE
                !
                ! SURFACE BOUNDARY WATER FLUX
                !
                ! DPTHW1,DPTHW2=surface water depth of source,destination
                ! ALT1,ALT2=elevation of source,destination
                ! XVOLT=excess surface water+ice
                ! VOLWG=ground surface water retention capacity
                ! DTBLX=natural water table depth
                ! QR1,HQR1=runoff, convective heat from runoff
                ! QR,HQR=hourly-accumulated runoff, convective heat from runoff
                ! QRM,QRV=runoff,velocity for erosion, solute transfer
                ! XN=direction
                !
                ! RUNOFF
                !
                DPTHW1=XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
                DPTHW2=VOLWG(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
                ALT1=ALTG(N2,N1)+DPTHW1
                ALT2=ALTG(N2,N1)+DPTHW2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)
                IF(ALT1.GT.ALT2 &
                  .AND.CDPTH(NU(N2,N1)-1,N2,N1)-DPTHW1.LT.DTBLX(N2,N1))THEN
                  QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
                  HQR1(N,NN,M5,M4)=cpw*TK1(0,N2,N1)*QR1(N,NN,M5,M4)
                  QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
                  HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
                ! F(N1.EQ.1.AND.N2.EQ.1)THEN
                !     WRITE(*,7744)'QRBND',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
                !    1,IRCHG(NN,N,N2,N1),QRM(M,N2,N1)
                !    2,QR1(N,NN,M5,M4),QR(N,NN,M5,M4)
                !    2,ALT1,ALT2,ALTG(N2,N1),FSLOPE(N,N2,N1)
                !    3,VOLWG(N2,N1),VOLWRX(N2,N1),ZM(N2,N1),ZS(N2,N1)
                !    4,VOLW1(0,N2,N1),VOLI1(0,N2,N1),DLYR(N,NUM(N2,N1),N2,N1)
                !    5,XVOLT(N2,N1)-VOLWG(N2,N1),DTBLX(N2,N1)
                !    6,DTBLX(N2,N1)-CDPTH(NU(N2,N1)-1,N2,N1)+DPTHW1
                ! 7,XVOLTM(M,N2,N1),XVOLWM(M,N2,N1)
!7744  FORMAT(A8,12I4,40E12.4)
                ! ENDIF
                !
                ! RUNON
!
              ELSEIF(CDPTH(NU(N2,N1)-1,N2,N1)-DPTHW1.GT.DTBLX(N2,N1))THEN
                VX=AMIN1(0.0,(DTBLX(N2,N1)-CDPTH(NU(N2,N1)-1,N2,N1)+DPTHW1) &
                  *AREA(3,NUM(N2,N1),N2,N1))
                QRM(M,N2,N1)=VX*XNPX
                QRV(M,N2,N1)=0.0
                QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
                HQR1(N,NN,M5,M4)=cpw*TK1(0,N2,N1)*QR1(N,NN,M5,M4)
                QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
                HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
                !     WRITE(*,7744)'QRBNB',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
                !    1,IRCHG(NN,N,N2,N1),QRM(M,N2,N1)
                !    2,QR1(N,NN,M5,M4),QR(N,NN,M5,M4),QRMN(M,N,NN,M5,M4)
                !    2,ALTG(N2,N1),FSLOPE(N,N2,N1)
                !    3,VOLWG(N2,N1),VOLWRX(N2,N1),ZM(N2,N1),ZS(N2,N1)
                !    4,VOLW1(0,N2,N1),VOLI1(0,N2,N1),DLYR(N,NUM(N2,N1),N2,N1)
                !    5,XVOLT(N2,N1)-VOLWG(N2,N1),DTBLX(N2,N1)
                !    6,CDPTH(NU(N2,N1)-1,N2,N1),DPTHW1
                !    7,XVOLTM(M,N2,N1),XVOLWM(M,N2,N1)
              ELSE
                QR1(N,NN,M5,M4)=0.0
                HQR1(N,NN,M5,M4)=0.0
              ENDIF
              QRMN(M,N,NN,M5,M4)=QR1(N,NN,M5,M4)
              IFLBM(M,N,NN,M5,M4)=0
!     WRITE(*,7745)'QRB',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
!    2,IFLBM(M,N,NN,M5,M4)
!    2,QRM(M,N2,N1),QRMN(M,N,NN,M5,M4),QR(N,NN,M5,M4)
!7745  FORMAT(A8,12I4,40E12.4)
!
        !     BOUNDARY SNOW FLUX
        !
        !     QS1,QW1,QI1=snow,water,ice transfer
        !     HQS1=convective heat transfer from snow,water,ice transfer
        !     QS,QW,QI=cumulative hourly snow,water,ice transfer
        !     HQS=cumulative hourly convective heat transfer from snow,water,ice transfer
!
              IF(NN.EQ.1)THEN
                QS1(N,M5,M4)=0.0
                QW1(N,M5,M4)=0.0
                QI1(N,M5,M4)=0.0
                HQS1(N,M5,M4)=0.0
!     QS(N,M5,M4)=QS(N,M5,M4)+QS1(N,M5,M4)
!     QW(N,M5,M4)=QW(N,M5,M4)+QW1(N,M5,M4)
!     QI(N,M5,M4)=QI(N,M5,M4)+QI1(N,M5,M4)
!     HQS(N,M5,M4)=HQS(N,M5,M4)+HQS1(N,M5,M4)
                QSM(M,N,M5,M4)=QS1(N,M5,M4)
              ENDIF
            ENDIF
          ELSE
            IF(N.NE.3)THEN
              QR1(N,NN,M5,M4)=0.0
              HQR1(N,NN,M5,M4)=0.0
            ENDIF
          ENDIF
!
          ! BOUNDARY SUBSURFACE WATER AND HEAT TRANSFER DEPENDING
          ! ON LEVEL OF WATER TABLE
!
          IF(VOLX(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
            IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
!
              ! IF NO WATER TABLE
              !
              ! IDTBL=water table flag
              ! THETA1,THETAX=water content ahead,behind wetting front
              ! K1,KL=pore water class ahead,behind wetting front
              ! CND1,CNDL=hydraulic conductivity ahead,behind wetting front
              ! FKSAT=reduction in soil surface Ksat from rainfall energy impact
              ! FLWL,FLWLX=lower boundary micropore water flux
              ! FLWHL=lower boundary macropore water flux
              ! HFLWL=convective heat from lower boundary water flux
              ! XH,XN,XNPH=rate constant,direction indicator,time step
              ! SLOPE=sin(vertical slope)=1
              ! RCHG*=boundary flags
!
              IF(IDTBL(N2,N1).EQ.0.OR.N.EQ.3)THEN
                THETA1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1) &
                  ,safe_adb(VOLW1(N3,N2,N1),VOLY(N3,N2,N1))))
                THETAX=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1) &
                  ,safe_adb(VOLWX1(N3,N2,N1),VOLY(N3,N2,N1))))
                K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1) &
                  -THETA1)/POROS(N3,N2,N1))+1))
                KL=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1) &
                  -THETAX)/POROS(N3,N2,N1))+1))
                IF(N3.EQ.NUM(NY,NX))THEN
                  CND1=HCND(N,K1,N3,N2,N1)*FKSAT
                ELSE
                  CND1=HCND(N,K1,N3,N2,N1)
                ENDIF
                CNDL=HCND(N,KL,N3,N2,N1)
                FLWL(N,M6,M5,M4)=AMIN1(VOLW1(N3,N2,N1)*XNPX &
                  ,XN*0.0098*(-ABS(SLOPE(N,N2,N1)))*CND1*AREA(3,N3,N2,N1)) &
                  *RCHGFU*RCHGFT*XNPH
                if(abs(FLWL(N,M6,M5,M4))>1.e20)then
                  write(*,*)'VOLW1(N3,N2,N1)*XNPX=',VOLW1(N3,N2,N1),XNPX
                  write(*,*)'XN=',XN,CND1
                  write(*,*)'RCHGFU*RCHGFT*XNPH=',RCHGFU,RCHGFT,XNPH
                  write(*,*)'at line',__LINE__
                  call endrun(trim(mod_filename)//'at line',__LINE__)
                endif
                FLWLX(N,M6,M5,M4)=FLWL(N,M6,M5,M4)
                FLWHL(N,M6,M5,M4)=AMIN1(VOLWH1(L,NY,NX)*XNPX &
                 ,XN*0.0098*(-ABS(SLOPE(N,N2,N1)))*CNDH1(L,NY,NX)*AREA(3,N3,N2,N1)) &
                 *RCHGFU*RCHGFT*XNPH
                HFLWL(N,M6,M5,M4)=cpw*TK1(N3,N2,N1)*(FLWL(N,M6,M5,M4)+FLWHL(N,M6,M5,M4))
                if(HFLWL(N,M6,M5,M4)<-1.e10_r8)then
                  write(*,*)'cpw*TK1(N3,N2,N1)*(FLWL(N,M6,M5,M4)+FLWHL(N,M6,M5,M4))',&
                    cpw,TK1(N3,N2,N1),FLWL(N,M6,M5,M4),FLWHL(N,M6,M5,M4)
                  call endrun(trim(mod_filename)//' at line',__LINE__)
                endif
                !     IF(I.EQ.336)THEN
                !     WRITE(*,4443)'ABV',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,XN
                !    2,FLWL(N,M6,M5,M4)
                !    2,VOLP2,RCHGFU,VOLX(N3,N2,N1),VOLW1(N3,N2,N1)
                !    3,VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1),VOLPH2,VOLI1(N3,N2,N1)
                !    4,VOLIH1(N3,N2,N1),VOLP1(N3,N2,N1),HFLWL(N,M6,M5,M4)
                !    5,PSISM1(N3,N2,N1),PSISE(N3,N2,N1),FLWHL(N,M6,M5,M4),DTBLD(N2,N1)
                !    6,SLOPE(N,N2,N1)
                !4443  FORMAT(A8,11I4,30E12.4)
                !     ENDIF
              ELSE
!
          !     MICROPORE DISCHARGE ABOVE WATER TABLE
          !
          !     IFLGU=micropore discharge flag to natural water table
          !     PSISWD=water potential from water table slope
          !     XN,RCHG*=direction indicator,boundary flag
          !     SLOPE=sin(lateral slope)
          !     DLYR=layer width
          !     DTBLG=water table slope
          !     PSISWT=water potential driving micropore discharge
          !     PSISA1,PSISO=matric,osmotic water potential
          !     DPTH,DTBLX=depth to layer midpoint,natural water table
          !     DPTHT=depth to internal water table
          !     FLWL=micropore discharge to natural water table
          !     HFLWL=convective heat from discharge to natural water table
          !     HCND=saturated hydraulic conductivity
          !     AREAU=fraction of layer below natural water table
!
                IF(IFLGU.EQ.0.AND.test_aneb(RCHGFT,0._r8))THEN
                  PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1) &
                    *(1.0-DTBLG(N2,N1))
                  PSISWT=AMIN1(0.0,-PSISA1(N3,N2,N1)-0.03*PSISO(N3,N2,N1) &
                    +0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1)) &
                    -0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1)))
                  IF(PSISWT.LT.0.0)PSISWT=PSISWT-PSISWD
                  FLWT=PSISWT*HCND(N,1,N3,N2,N1)*AREA(N,N3,N2,N1) &
                    *(1.0-AREAU(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
                  FLWL(N,M6,M5,M4)=XN*FLWT
                  if(abs(FLWL(N,M6,M5,M4))>1.e20)then
                    write(*,*)'XN*FLWT=',XN,FLWT
                    write(*,*)'at line',__LINE__
                    call endrun(trim(mod_filename)//'at line',__LINE__)
                  endif
                  FLWLX(N,M6,M5,M4)=XN*FLWT
                  HFLWL(N,M6,M5,M4)=cpw*TK1(N3,N2,N1)*XN*FLWT
                  if(HFLWL(N,M6,M5,M4)<-1.e10_r8)then
                    write(*,*)'cpw*TK1(N3,N2,N1)*XN*FLWT',cpw,TK1(N3,N2,N1),FLWT
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif
                !     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
                !     WRITE(*,4445)'DISCHMI',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU
                !    2,XN,FLWL(N,M6,M5,M4),FLWT,PSISWT,PSISWD,HCND(N,1,N3,N2,N1)
                !    3,AREA(N,N3,N2,N1),AREAU(N3,N2,N1),RCHGFU,RCHGFT
                !    4,PSISE(N3,N2,N1),PSISA1(N3,N2,N1),DPTH(N3,N2,N1),DTBLX(N2,N1)
                !    2,PSISO(N3,N2,N1),0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1))
                !    3,0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1))
                !4445  FORMAT(A8,12I4,30E14.6)
                  !     ENDIF
                ELSE
                  FLWL(N,M6,M5,M4)=0.0
                  FLWLX(N,M6,M5,M4)=0.0
                  HFLWL(N,M6,M5,M4)=0.0
                ENDIF
!
          !     MACROPORE DISCHARGE ABOVE WATER TABLE
          !
        !     IFLGUH=macropore discharge flag to natural water table
        !     PSISWD=water potential from water table slope
        !     XN,RCHG*=direction indicator,boundary flag
        !     SLOPE=sin(lateral slope)
        !     DLYR=layer width
        !     DTBLG=water table slope
        !     PSISWTH=water potential driving macropore discharge
        !     PSISO=osmotic water potential
        !     DPTHH,DTBLX=depth to layer macropore water,natural water table
        !     DPTHT=depth to internal water table
        !     FLWTH,FLWTHL=macropore discharge unltd,ltd by macropore water
        !     CNDH1=macropore hydraulic conductivity
        !     FLWHL=macropore discharge to natural water table
        !     HFLWL=convective heat from discharge to natural water table
        !     HCND=saturated hydraulic conductivity
          !     AREAU=fraction of layer below natural water table
!
                IF(IFLGUH.EQ.0.AND.test_aneb(RCHGFT,0._r8) &
                  .AND.VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
                  PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1) &
                    *(1.0-DTBLG(N2,N1))
                  PSISWTH=-0.03*PSISO(N3,N2,N1) &
                    +0.0098*(DPTHH-DTBLX(N2,N1)) &
                    -0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
                  IF(PSISWTH.LT.0.0)PSISWTH=PSISWTH-PSISWD
                  FLWTH=PSISWTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1) &
                    *(1.0-AREAU(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
                  FLWTHL=AMAX1(FLWTH,AMIN1(0.0,-(VOLWH1(N3,N2,N1)*XNPX &
                    +FLWHL(3,N3,N2,N1)-FLWHL(3,N3+1,N2,N1))))
                  FLWHL(N,M6,M5,M4)=XN*FLWTHL
                  HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+cpw*TK1(N3,N2,N1)*XN*FLWTHL
                  if(HFLWL(N,M6,M5,M4)<-1.e10_r8)then
                    write(*,*)'cpw*TK1(N3,N2,N1)*XN*FLWTHL',TK1(N3,N2,N1),FLWTHL
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif
                    !     WRITE(*,4446)'DISCHMA',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGUH
                    !    2,XN,FLWHL(N,M6,M5,M4),FLWTHL,FLWTH,PSISWTH,CNDH1(N3,N2,N1)
                    !    3,DPTH(N3,N2,N1),DLYR(3,N3,N2,N1),DPTHH,VOLWH1(N3,N2,N1)
                    !    4,VOLIH1(L,NY,NX),VOLAH1(N3,N2,N1),DTBLX(N2,N1),PSISWD
                    !4446  FORMAT(A8,12I4,30E14.6)
                ELSE
                  FLWHL(N,M6,M5,M4)=0.0
                ENDIF
!
                !     MICROPORE DISCHARGE ABOVE TILE DRAIN
                !
                !     IFLGD=micropore discharge flag to artificial water table
                !     PSISWD=water potential from water table slope
                !     XN,RCHG*=direction indicator,boundary flag
                !     SLOPE=sin(lateral slope)
                !     DLYR=layer width
                !     DTBLG=water table slope
                !     PSISWT=water potential driving micropore discharge
                !     PSISA1,PSISO=matric,osmotic water potential
                !     DPTH,DTBLY=depth to layer midpoint,artificial water table
                !     DPTHT=depth to internal water table
                !     FLWL=micropore discharge to natural+artificial water table
                !     HFLWL=convective heat from dischg to natural+artifl water table
                !     HCND=saturated hydraulic conductivity
                !     AREAUD=fraction of layer below artificial water table
!
                IF(IFLGD.EQ.0.AND.test_aneb(RCHGFT,0._r8))THEN
                  PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1) &
                    *(1.0-DTBLG(N2,N1))
                  PSISWT=AMIN1(0.0,-PSISA1(N3,N2,N1)-0.03*PSISO(N3,N2,N1) &
                    +0.0098*(DPTH(N3,N2,N1)-DTBLY(N2,N1)) &
                    -0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1)))
                  IF(PSISWT.LT.0.0)PSISWT=PSISWT-PSISWD
                  FLWT=PSISWT*HCND(N,1,N3,N2,N1)*AREA(N,N3,N2,N1) &
                    *(1.0-AREAUD(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
                  FLWL(N,M6,M5,M4)=FLWL(N,M6,M5,M4)+XN*FLWT
                  if(abs(FLWL(N,M6,M5,M4))>1.e20)then
                    write(*,*)'aXN*FLWT=',XN,FLWT
                    write(*,*)'at line',__LINE__
                    call endrun(trim(mod_filename)//'at line',__LINE__)
                  endif
                  FLWLX(N,M6,M5,M4)=FLWLX(N,M6,M5,M4)+XN*FLWT
                  HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+cpw*TK1(N3,N2,N1)*XN*FLWT
                  if(HFLWL(N,M6,M5,M4)<-1.e10)then
                    write(*,*)'cpw*TK1(N3,N2,N1)*XN*FLWT',TK1(N3,N2,N1),FLWT
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif
            !     IF(L.EQ.12)THEN
            !     WRITE(*,4445)'DISCHMD',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGD
            !    2,XN,FLWL(N,M6,M5,M4),FLWT,PSISWT,PSISWD
            !    3,0.0098*(DPTH(N3,N2,N1)-DTBLY(N2,N1)),RCHGFU,RCHGFT
            !    3,HCND(N,1,N3,N2,N1),AREA(N,N3,N2,N1),AREAUD(N3,N2,N1)
            !    4,PSISE(N3,N2,N1),PSISA1(N3,N2,N1),PSISO(N3,N2,N1)
            !    5,DPTH(N3,N2,N1),DTBLY(N2,N1),DPTHT(N2,N1)
            !    3,0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1))
            !     ENDIF
                ENDIF
            !
            !     MACROPORE DISCHARGE ABOVE TILE DRAIN
            !
            !     IFLGDH=macropore discharge flag to artificial water table
            !     PSISWD=water potential from water table slope
            !     XN,RCHG*=direction indicator,boundary flag
            !     SLOPE=sin(lateral slope)
            !     DLYR=layer width
            !     DTBLG=water table slope
            !     PSISWTH=water potential driving macropore discharge
            !     PSISO=osmotic water potential
            !     DPTHH,DTBLY=depth to layer macropore water,artificl water table
            !     DPTHT=depth to internal water table
            !     FLWTH,FLWTHL=macropore discharge unltd,ltd by macropore water
            !     CNDH1=macropore hydraulic conductivity
            !     FLWHL=macropore discharge to artificial water table
            !     HFLWL=convective heat from discharge to artificial water table
            !     HCND=saturated hydraulic conductivity
            !     AREAUD=fraction of layer below artificial water table
!
                IF(IFLGDH.EQ.0.AND.test_aneb(RCHGFT,0._r8) &
                  .AND.VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
                  PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1) &
                    *(1.0-DTBLG(N2,N1))
                  PSISWTH=-0.03*PSISO(N3,N2,N1) &
                    +0.0098*(DPTHH-DTBLY(N2,N1)) &
                    -0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
                  IF(PSISWTH.LT.0.0)PSISWTH=PSISWTH-PSISWD
                  FLWTH=PSISWTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1) &
                    *(1.0-AREAUD(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
                  FLWTHL=AMAX1(FLWTH,AMIN1(0.0 &
                    ,-(VOLWH1(N3,N2,N1)*XNPX+FLWHL(3,N3,N2,N1)-FLWHL(3,N3+1,N2,N1))))
                  FLWHL(N,M6,M5,M4)=FLWHL(N,M6,M5,M4)+XN*FLWTHL
                  HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+cpw*TK1(N3,N2,N1)*XN*FLWTHL
                  if(HFLWL(N,M6,M5,M4)<-1.e10_r8)then
                    write(*,*)'cpw*TK1(N3,N2,N1)*XN*FLWTHL',TK1(N3,N2,N1),FLWTHL
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif
                  !     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
            !     WRITE(*,4446)'DISCHDH',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGDH
            !    2,XN,FLWHL(N,M6,M5,M4),FLWTHL,FLWTH,PSISWTH,CNDH1(N3,N2,N1)
            !    3,DPTH(N3,N2,N1),DLYR(3,N3,N2,N1),DPTHH,VOLWH1(N3,N2,N1)
            !    4,VOLIH1(L,NY,NX),VOLAH1(N3,N2,N1),DTBLY(N2,N1),PSISWD
            !    2,0.0098*(DPTHH-DTBLY(N2,N1))
            !    3,0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
            !     ENDIF
                ENDIF
!
            !     MICROPORE RECHARGE BELOW WATER TABLE
            !
            !     DPTH,DTBLX=depth to layer midpoint,natural water table
            !     DPTHA=active layer depth
            !     VOLP2=air volume
            !     PSISWD=water potential from water table slope
            !     XN,RCHG*=direction indicator,boundary flag
            !     SLOPE=sin(lateral slope)
            !     DLYR=layer width
            !     DTBLG=water table slope
            !     PSISUT=water potential driving micropore recharge
            !     PSISA1,PSISO=matric,osmotic water potential
            !     DPTH,DTBLX=depth to layer midpoint,natural water table
            !     FLWU,FLWUL=micropore recharge unltd,ltd by micropore air volume
            !     FLWL=micropore recharge from natural water table
            !     HFLWL=convective heat from recharge from natural water table
            !     HCND=saturated hydraulic conductivity
            !     AREAU=fraction of layer below natural water table
            !
                IF(DPTH(N3,N2,N1).GE.DTBLX(N2,N1) &
                  .AND.DPTHA(N2,N1).GT.DTBLX(N2,N1) &
                  .AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1) &
                  .AND.(VOLP2.GT.ZEROS2(N2,N1).OR.BKDS(N3,N2,N1).LE.ZERO) &
                  .AND.VOLP1Z(N3,N2,N1).GT.0.0 &
                  .AND.test_aneb(RCHGFT,0._r8))THEN
                  PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1) &
                    *(1.0-DTBLG(N2,N1))
                  PSISUT=AMAX1(0.0,-PSISA1(N3,N2,N1)-0.03*PSISO(N3,N2,N1) &
                    +0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1)))
                  IF(PSISUT.GT.0.0)PSISUT=PSISUT+PSISWD
                  FLWU=PSISUT*HCND(N,1,N3,N2,N1)*AREA(N,N3,N2,N1) &
                    *AREAU(N3,N2,N1)/(RCHGFU+1.0)*RCHGFT*XNPH
                  IF(BKDS(N3,N2,N1).GT.ZERO)THEN
                    FLWUL=AMIN1(FLWU,VOLP2)
                    FLWUX=AMIN1(FLWU,VOLPX2)
                  ELSE
                    FLWUL=FLWU
                    FLWUX=FLWU
                  ENDIF
                  FLWL(N,M6,M5,M4)=FLWL(N,M6,M5,M4)+XN*FLWUL
                  if(abs(FLWL(N,M6,M5,M4))>1.e20_r8)then
                    write(*,*)'bXN*FLWUL=',XN,FLWUL
                    write(*,*)'at line',__LINE__
                    call endrun(trim(mod_filename)//'at line',__LINE__)
                  endif
                  FLWLX(N,M6,M5,M4)=FLWLX(N,M6,M5,M4)+XN*FLWUX
                  HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+cpw*TK1(N3,N2,N1)*XN*FLWUL
                  if(HFLWL(N,M6,M5,M4)<-1.e10)then
                    write(*,*)'cpw*TK1(N3,N2,N1)*XN*FLWUL',TK1(N3,N2,N1),FLWUL
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif
                !     IF(J.EQ.15)THEN
                !     WRITE(*,4444)'RECHGMI',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU
                !    2,XN,FLWL(N,M6,M5,M4),FLWU,FLWUL
                !    2,PSISUT,PSISA1(N3,N2,N1),PSISO(N3,N2,N1)
                !    3,0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1)),RCHGFU,RCHGFT
                !    4,HCND(N,1,N3,N2,N1),AREA(N,N3,N2,N1),AREAU(N3,N2,N1)
                !    5,VOLP2,VOLA1(L,NY,NX),VOLW1(L,NY,NX),VOLI1(L,NY,NX),DTBLX(N2,N1)
                !    6,DPTH(N3,N2,N1),DPTHA(N2,N1),DPTHT(N2,N1)
                !    6,DTBLD(N2,N1),VOLW1(N3,N2,N1),VOLI1(N3,N2,N1)
                !    7,VOLX(N3,N2,N1),VOLP1(N3,N2,N1),VOLP1Z(N3,N2,N1)
                !    8,FLWL(3,N3,N2,N1),FLWL(3,N6,N5,N4),QR1(N,NN,M5,M4)
                !    9,DLYR(N,N3,N2,N1),DLYR(3,N3,N2,N1),PSISWD
                !    1,CDPTH(N3,N2,N1),AREA(N,N3,N2,N1),SLOPE(N,N2,N1)
                !4444  FORMAT(A8,12I4,40E14.6)
                !     ENDIF
                  VOLP2=VOLP2-XN*FLWL(N,M6,M5,M4)
                  VOLPX2=VOLPX2-XN*FLWLX(N,M6,M5,M4)
                ENDIF
!
                !     MACROPORE RECHARGE BELOW WATER TABLE
                !
                !     PSISWD=water potential from water table slope
                !     XN,RCHG*=direction indicator,boundary flag
                !     SLOPE=sin(lateral slope)
                !     DLYR=layer width
                !     DTBLG=water table slope
                !     PSISUTH=water potential driving macropore recharge
                !     PSISO=osmotic water potential
                !     DPTHH,DTBLX=depth to layer macropore water,natural water table
                !     DPTHT=depth to internal water table
                !     CNDH1=macropore hydraulic conductivity
                !     FLWUH,FLWUHL=macropore recharge unltd,ltd by macropore air volume
                !     FLWHL=macropore discharge to natural water table
                !     HFLWL=convective heat from discharge to natural water table
                !     HCND=saturated hydraulic conductivity
                !     AREAU=fraction of layer below natural water table
!
                IF(DPTHH.GT.DTBLX(N2,N1) &
                  .AND.DPTHA(N2,N1).GT.DTBLX(N2,N1) &
                  .AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1) &
                  .AND.VOLPH2.GT.ZEROS2(NY,NX) &
                  .AND.test_aneb(RCHGFT,0.0_r8))THEN
                  PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1) &
                    *(1.0-DTBLG(N2,N1))
                  PSISUTH=-0.03*PSISO(N3,N2,N1) &
                    +0.0098*(DPTHH-DTBLX(N2,N1))
                  IF(PSISUTH.GT.0.0)PSISUTH=PSISUTH+PSISWD
                  FLWUH=PSISUTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1) &
                    *AREAU(N3,N2,N1)/(RCHGFU+1.0)*RCHGFT*XNPH
                  FLWUHL=AMIN1(FLWUH,VOLPH2*XNPX)
                  FLWHL(N,M6,M5,M4)=FLWHL(N,M6,M5,M4)+XN*FLWUHL
                  HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+cpw*TK1(N3,N2,N1)*XN*FLWUHL
                  if(HFLWL(N,M6,M5,M4)<-1.e10_r8)then
                    write(*,*)'cpw*TK1(N3,N2,N1)*XN*FLWUHL',TK1(N3,N2,N1),FLWUHL
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif
                !     IF(M6.EQ.7)THEN
                !     WRITE(*,4447)'RECHGMA',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU
                !    2,XN,FLWHL(N,M6,M5,M4),FLWUH,FLWUHL,DPTHH,PSISUTH,VOLPH2
                !    3,CNDH1(N3,N2,N1),DTBLX(N2,N1),CDPTH(N3,N2,N1),DPTHT(N2,N1)
                !    6,DTBLD(N2,N1),DPTH(N3,N2,N1),VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1)
                !    8,FLWHL(3,N3,N2,N1),FLWHL(3,N3+1,N2,N1),RCHGFU,AREA(N,N3,N2,N1)
                !    9,DLYR(N,N3,N2,N1),DLYR(3,N3,N2,N1),PSISWD
                !    1,SLOPE(N,N2,N1),AREAU(N3,N2,N1)
                !4447  FORMAT(A8,12I4,40E14.6)
                !     ENDIF
                  VOLPH2=VOLPH2-XN*FLWHL(N,M6,M5,M4)
                ENDIF
              ENDIF
!
        !     SUBSURFACE HEAT SOURCE/SINK
        !
        !     HFLWL=heat flux across lower boundary
        !     TK1=lower boundary soil temperature
        !     TKSD=deep source/sink temperature from geothermal flux
        !     TCNDG=thermal conductivity below lower boundary
        !     DPTHSK,CDPTH=depth of thermal sink/source, lower boundary
        !
              IF(N.EQ.3.AND.IETYP(N2,N1).NE.-2)THEN
                HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+(TK1(N3,N2,N1) &
                  -TKSD(N2,N1))*TCNDG/(DPTHSK(N2,N1)-CDPTH(N3,N2,N1)) &
                  *AREA(N,N3,N2,N1)*XNPH
              ENDIF
              FLW(N,M6,M5,M4)=FLW(N,M6,M5,M4)+FLWL(N,M6,M5,M4)
              FLWX(N,M6,M5,M4)=FLWX(N,M6,M5,M4)+FLWLX(N,M6,M5,M4)
              FLWH(N,M6,M5,M4)=FLWH(N,M6,M5,M4)+FLWHL(N,M6,M5,M4)
              HFLW(N,M6,M5,M4)=HFLW(N,M6,M5,M4)+HFLWL(N,M6,M5,M4)
              FLWM(M,N,M6,M5,M4)=FLWL(N,M6,M5,M4)
              FLWHM(M,N,M6,M5,M4)=FLWHL(N,M6,M5,M4)
            ENDIF
          ELSE
            FLWL(N,M6,M5,M4)=0.0
            FLWLX(N,M6,M5,M4)=0.0
            FLWHL(N,M6,M5,M4)=0.0
            HFLWL(N,M6,M5,M4)=0.0
            FLWM(M,N,M6,M5,M4)=0.0
            FLWHM(M,N,M6,M5,M4)=0.0
          ENDIF
9575    CONTINUE
    !
    !     NET WATER AND HEAT FLUXES IN RUNOFF AND SNOW DRIFT
    !
    !     TQR1,THQR1=net runoff,convective heat from runoff
    !     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
    !     QR1,HQR1=runoff, convective heat from runoff
    !     QS1,QW1,QI1=snow,water,ice transfer
    !     HQS1=convective heat transfer from snow,water,ice transfer
!
        IF(L.EQ.NUM(N2,N1).AND.N.NE.3)THEN
          DO 1202 NN=1,2
            TQR1(N2,N1)=TQR1(N2,N1)+QR1(N,NN,N2,N1)
            THQR1(N2,N1)=THQR1(N2,N1)+HQR1(N,NN,N2,N1)
            IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
              TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5,N4)
              THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5,N4)
            ENDIF
            !     IF(I.GT.350.AND.NX.EQ.1)THEN
        !     WRITE(*,6631)'TQR1',I,J,M,N1,N2,N4,N5,N,NN
        !    2,IFLBM(M,N,NN,N5,N4),TQR1(N2,N1),THQR1(N2,N1)
        !    2,QR1(N,NN,N2,N1),QR1(N,NN,N5,N4)
        !    3,QR(N,NN,N2,N1),QR(N,NN,N5,N4)
        !    2,HQR1(N,NN,N2,N1),HQR1(N,NN,N5,N4)
        !    3,HQR(N,NN,N2,N1),HQR(N,NN,N5,N4)
            !6631  FORMAT(A8,10I4,12E12.4)
            !     ENDIF
            IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
              TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5B,N4B)
              THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5B,N4B)
        !     IF(I.GT.350.AND.NX.EQ.1)THEN
        !     WRITE(*,6631)'TQRB1',I,J,M,N1,N2,N4B,N5B,N,NN
        !    2,IFLBM(M,N,NN,N5B,N4B),TQR1(N2,N1),THQR1(N2,N1)
        !    2,QR1(N,NN,N5B,N4B),HQR1(N,NN,N5B,N4B)
        !    2,QR(N,NN,N5B,N4B),HQR(N,NN,N5B,N4B)
        !     ENDIF
            ENDIF
            IF(M.EQ.NPH)THEN
              IFLBH(N,NN,N5,N4)=IFLBM(M,N,NN,N5,N4)
              IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
                IFLBH(N,NN,N5B,N4B)=IFLBM(M,N,NN,N5B,N4B)
              ENDIF
            ENDIF
1202      CONTINUE
          TQS1(N2,N1)=TQS1(N2,N1)+QS1(N,N2,N1)-QS1(N,N5,N4)
          TQW1(N2,N1)=TQW1(N2,N1)+QW1(N,N2,N1)-QW1(N,N5,N4)
          TQI1(N2,N1)=TQI1(N2,N1)+QI1(N,N2,N1)-QI1(N,N5,N4)
          THQS1(N2,N1)=THQS1(N2,N1)+HQS1(N,N2,N1)-HQS1(N,N5,N4)
        ENDIF
!
        !     NET WATER AND HEAT FLUXES THROUGH SOIL AND SNOWPACK
        !
        !     TFLWL,THFLWL=net water micropore,macropore flux
        !     THFLWL=net convective+conductive heat flux
        !     FLWL =micropore water,heat flux
        !     FLWHL=macropore water,heat flux
        !     HFLWL=soil heat flux
!
        IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
          DO 1200 LL=N6,NL(N5,N4)
            IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
              N6=LL
              GO TO 1201
            ENDIF
1200      CONTINUE
1201      CONTINUE
          IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            TFLWL(N3,N2,N1)=TFLWL(N3,N2,N1)+FLWL(N,N3,N2,N1)-FLWL(N,N6,N5,N4)
            if(abs(TFLWL(N3,N2,N1))>1.e20_r8)then
              write(*,*)'TFLWL(N3,N2,N1)=',TFLWL(N3,N2,N1)
              write(*,*)'FLWL(N,N3,N2,N1)=',FLWL(N,N3,N2,N1),FLWL(N,N6,N5,N4)
              write(*,*)'N3,n2,n1',n,n3,n2,n1
              write(*,*)'at line',__LINE__
              call endrun(trim(mod_filename)//'at line',__LINE__)
            endif
            TFLWLX(N3,N2,N1)=TFLWLX(N3,N2,N1)+FLWLX(N,N3,N2,N1)-FLWLX(N,N6,N5,N4)
            TFLWHL(N3,N2,N1)=TFLWHL(N3,N2,N1)+FLWHL(N,N3,N2,N1)-FLWHL(N,N6,N5,N4)
            THFLWL(N3,N2,N1)=THFLWL(N3,N2,N1)+HFLWL(N,N3,N2,N1)-HFLWL(N,N6,N5,N4)
    !     IF(N3.EQ.1)THEN
    !     WRITE(*,3378)'THFLW',I,J,M,N1,N2,N3,N4,N5,N6,N
    !    2,TFLWL(N3,N2,N1),FLWL(N,N3,N2,N1),FLWL(N,N6,N5,N4)
    !    2,THFLWL(N3,N2,N1),HFLWL(N,N3,N2,N1),HFLWL(N,N6,N5,N4)
    !    2,FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
    !    2,HFLW(N,N3,N2,N1),HFLW(N,N6,N5,N4)
    !3378  FORMAT(A8,10I4,20E14.6)
    !     ENDIF
          ELSE
            TFLWL(N3,N2,N1)=0.0
            TFLWLX(N3,N2,N1)=0.0
            TFLWHL(N3,N2,N1)=0.0
            THFLWL(N3,N2,N1)=0.0
          ENDIF
        ENDIF
9580  CONTINUE
!
!     INFILTRATION OF WATER FROM MACROPORES INTO MICROPORES
!
!     VOLWH1=macropore volume
!     FINHX,FINHL=macro-micropore transfer unltd,ltd by water,air volume
!     FINHM=macro-micropore transfer for use in trnsfr.f
!     HCND=saturated hydraulic conductivity
!     PSISE,PSISA1=air entry,matric water potentials
!     PHOL,HRAD=path length between,radius of macropores from hour1.f
!     XNPH=time step
!     VOLW1X,VOLP1X=current micropore water,air volume
!     VOLWH1X,VOLPH1X=current macropore water,air volume
!
      IF(VOLWH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
        FINHX=6.283*HCND(2,1,N3,N2,N1)*AREA(3,N3,N2,N1) &
          *(PSISE(N3,N2,N1)-PSISA1(N3,N2,N1)) &
          /LOG(PHOL(N3,N2,N1)/HRAD(N3,N2,N1))*XNPH
        VOLW1X=VOLW1(N3,N2,N1)+TFLWL(N3,N2,N1)+FLU1(N3,N2,N1)
        VOLP1X=AMAX1(0.0,VOLA1(N3,N2,N1)-VOLW1X-VOLI1(N3,N2,N1))
        VOLWH1X=VOLWH1(N3,N2,N1)+TFLWHL(N3,N2,N1)
        VOLPH1X=AMAX1(0.0,VOLAH1(N3,N2,N1)-VOLWH1X-VOLIH1(N3,N2,N1))
        IF(FINHX.GT.0.0)THEN
          FINHL(N3,N2,N1)=AMAX1(0.0,AMIN1(FINHX,VOLWH1X,VOLP1X))
        ELSE
          FINHL(N3,N2,N1)=AMIN1(0.0,AMAX1(FINHX,-VOLPH1X,-VOLW1X))
        ENDIF
        FINHM(M,N3,N2,N1)=FINHL(N3,N2,N1)
        FINH(N3,N2,N1)=FINH(N3,N2,N1)+FINHL(N3,N2,N1)
!     IF(NX.EQ.1.AND.NY.EQ.1)THEN
!     WRITE(*,3366)'FINHL',I,J,M,N4,N5,N6,IFLGH,FINHL(N3,N2,N1)
!    3,FINHX,VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1),VOLP1(N3,N2,N1)
!    4,PSISA1(N3,N2,N1),HCND(2,1,N3,N2,N1),PHOL(N3,N2,N1)
!    5,HRAD(N3,N2,N1)
!3366  FORMAT(A8,7I4,20E12.4)
!     ENDIF
      ELSE
        FINHL(N3,N2,N1)=0.0
        FINHM(M,N3,N2,N1)=0.0
      ENDIF
!
!     FREEZE-THAW IN SOIL LAYER MICROPORE FROM NET CHANGE IN SOIL
!     LAYER HEAT STORAGE
!
!     TFREEZ=micropore freezing temperature
!     PSISA1,PSISO=micropore matric,osmotic potential
!     VOLW1*,VOLI1=micropore water,ice volume
!     VOLWH1*,VOLIH1=macropore water,ice volume
!     VHCP1X,VHCP1AX,VHCP1BX=total soil,micropore,macropore heat capacity
!     VHCM=soil solid volumetric heat capacity
!     TK1*=soil temperature
!     THFLWL=total soil conductive, convective heat flux
!     HWFLU1=subsurface convective heat flux
!     TFLX1,TFLX=latent heat from micro freeze-thaw unltd,ltd by water,ice
!     WFLXL=soil water flux from micropore freeze-thaw
!     TFLXH1,TFLXL=latent heat from macro freeze-thaw unltd,ltd by water,ice
!     WFLXLH=soil water flux from macropore freeze-thaw
!
      PSISMX=PSISA1(N3,N2,N1)+PSISO(N3,N2,N1)
      TFREEZ=-9.0959E+04/(PSISMX-333.0)
      VOLW1X=VOLW1(N3,N2,N1)+TFLWL(N3,N2,N1)+FINHL(N3,N2,N1)+FLU1(N3,N2,N1)
      VOLWH1X=VOLWH1(N3,N2,N1)+TFLWHL(N3,N2,N1)-FINHL(N3,N2,N1)
      ENGY1=VHCP1(N3,N2,N1)*TK1(N3,N2,N1)
      VHCP1X=VHCM(N3,N2,N1)+cpw*(VOLW1X+VOLWH1X) &
        +cpi*(VOLI1(N3,N2,N1)+VOLIH1(N3,N2,N1))
      IF(VHCP1X.GT.ZEROS(NY,NX))THEN
        TK1X=(ENGY1+THFLWL(N3,N2,N1)+HWFLU1(N3,N2,N1))/VHCP1X
        IF((TK1X.LT.TFREEZ &
          .AND.VOLW1(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1)) &
          .OR.(TK1X.GT.TFREEZ &
          .AND.VOLI1(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1)))THEN
          VHCP1AX=VHCM(N3,N2,N1)+cpw*VOLW1X+cpi*VOLI1(N3,N2,N1)
          TFLX1=VHCP1AX*(TFREEZ-TK1X)/((1.0+6.2913E-03*TFREEZ) &
            *(1.0-0.10*PSISMX))*XNPX
          IF(TFLX1.LT.0.0_r8)THEN
            TFLX=AMAX1(-333.0*DENSI*VOLI1(N3,N2,N1)*XNPX,TFLX1)
          ELSE
            TFLX=AMIN1(333.0*VOLW1X*XNPX,TFLX1)
          ENDIF
            WFLXL(N3,N2,N1)=-TFLX/333.0
          ELSE
            TFLX=0.0
            WFLXL(N3,N2,N1)=0.0
          ENDIF
        ELSE
          TFLX=0.0
          WFLXL(N3,N2,N1)=0.0
        ENDIF
!
!     FREEZE-THAW IN SOIL LAYER MACROPORE FROM NET CHANGE IN SOIL
!     LAYER HEAT STORAGE
!
        IF((TK1X.LT.TFice.AND.VOLWH1(N3,N2,N1) &
          .GT.ZERO*VOLT(N3,N2,N1)).OR.(TK1X.GT.TFice &
          .AND.VOLIH1(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1)))THEN
          VHCP1BX=cpw*VOLWH1X+cpi*VOLIH1(L,NY,NX)
          TFLXH1=VHCP1BX*(TFREEZ-TK1X)/((1.0_r8+6.2913E-03_r8*TFREEZ) &
            *(1.0_r8-0.10_r8*PSISMX))*XNPX
          IF(TFLXH1.LT.0.0_r8)THEN
            TFLXH=AMAX1(-333.0_r8*DENSI*VOLIH1(N3,N2,N1)*XNPX,TFLXH1)
          ELSE
            TFLXH=AMIN1(333.0_r8*VOLWH1X*XNPX,TFLXH1)
          ENDIF
          WFLXLH(N3,N2,N1)=-TFLXH/333.0_r8
        ELSE
          TFLXH=0.0_r8
          WFLXLH(N3,N2,N1)=0.0_r8
        ENDIF
        TFLXL(N3,N2,N1)=TFLX+TFLXH
          !
          !     TOTAL AND HOURLY ACCUMULATED FREEZE-THAW FLUXES
        !
        !     THAW,THAWH=hourly accumulated freeze-thaw flux in micropores,macropores
        !     HTHAW=hourly accumulated freeze-thaw latent heat flux
        !     TWFLXL,TWFLXH=total accumulated freeze-thaw in micropores,macropores
          !     TTFLXL=total latent heat flux
!
        THAW(N3,N2,N1)=THAW(N3,N2,N1)+WFLXL(N3,N2,N1)
        THAWH(N3,N2,N1)=THAWH(N3,N2,N1)+WFLXLH(N3,N2,N1)
        HTHAW(N3,N2,N1)=HTHAW(N3,N2,N1)+TFLXL(N3,N2,N1)
        TWFLXL(N3,N2,N1)=TWFLXL(N3,N2,N1)+WFLXL(N3,N2,N1)
        TWFLXH(N3,N2,N1)=TWFLXH(N3,N2,N1)+WFLXLH(N3,N2,N1)
        TTFLXL(N3,N2,N1)=TTFLXL(N3,N2,N1)+TFLXL(N3,N2,N1)
!     IF(N3.EQ.1)THEN
!     WRITE(*,4359)'TFLX',I,J,M,N1,N2,N3,TWFLXL(N3,N2,N1)
!    2,WFLXL(N3,N2,N1),THAW(N3,N2,N1),VHCP1X,VHCP1AX
!    2,TFLX1,TFLX,TFREEZ,TK1(N3,N2,N1),TK1X,VHCP1(N3,N2,N1)
!    3,TFLXL(N3,N2,N1),TTFLXL(N3,N2,N1)
!    2,TFLWL(N3,N2,N1),FINHL(N3,N2,N1),FLU1(N3,N2,N1),THFLWL(N3,N2,N1)
!    4,HWFLU1(N3,N2,N1),VOLW1X,VOLW1(N3,N2,N1),VOLI1(N3,N2,N1)
!    5,PSISA1(N3,N2,N1),PSISM(N3,N2,N1),PSISO(N3,N2,N1),PSISMX
!    5,VOLWH1(N3,N2,N1),VOLIH1(N3,N2,N1)
!    5,FGRD(N3,N2,N1),FMAC(N3,N2,N1)
!4359  FORMAT(A8,6I4,40E14.6)
!     ENDIF
!
!     DISSIPATE WETTING FRONT
!
!     VOLW1=soil micropore water content
!     VOLWX1=soil micropore water content behind wetting front
!     FLWVL=water flux from wetted to drier soil
!
9585  CONTINUE
9590  CONTINUE
9595  CONTINUE

  end subroutine WaterHeatExchThruBoundaryFlow
!------------------------------------------------------------------------------------------

  subroutine UpdateStateSolutionNotNPH(M,NHW,NHE,NVN,NVS)
  !
  !Description
  ! Early exit of the watsub solver
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX
  integer :: L,NUX,LL,Ls

  real(r8) :: tk1pres

! begin_execution
!  if(curday>=176)then
!    write(*,*)'UpdateStateSolutionNotNPH at line',__LINE__,M,tk1(8,NVN,NHW)
!  endif
  DO 9795 NX=NHW,NHE
    DO 9790 NY=NVN,NVS
      !     SNOWPACK WATER, ICE, SNOW AND TEMPERATURE
      !
      !     VOLS0,VOLW0,VOLI0=snow,water,ice volumes in snowpack
      !     TFLWS,TFLWW=net snow,water flux
      !     WFLXS,WFLXI=snow-water,ice-water freeze-thaw flux
      !     DENSI=ice density
      !     VHCPWM=snowpack volumetric heat capacity
      !     TK0=snowpack temperature
      !     THFLWW=total snowpack conductive+convective heat flux
      !     TFLX0=snowpack latent heat flux from freeze-thaw
      !
      DO 9780 L=1,JS
        VOLS0(L,NY,NX)=VOLS0(L,NY,NX)+TFLWS(L,NY,NX)-WFLXS(L,NY,NX)
        VOLW0(L,NY,NX)=VOLW0(L,NY,NX)+TFLWW(L,NY,NX)+WFLXS(L,NY,NX)+WFLXI(L,NY,NX)
        VOLI0(L,NY,NX)=VOLI0(L,NY,NX)-WFLXI(L,NY,NX)/DENSI
        ENGY0=VHCPWM(M,L,NY,NX)*TK0(L,NY,NX)
        VHCPWM(M+1,L,NY,NX)=cps*VOLS0(L,NY,NX)+cpw*VOLW0(L,NY,NX) &
          +cpi*VOLI0(L,NY,NX)
        IF(VHCPWM(M+1,L,NY,NX).GT.VHCPWX(NY,NX))THEN
          TK0(L,NY,NX)=(ENGY0+THFLWW(L,NY,NX)+TFLX0(L,NY,NX))/VHCPWM(M+1,L,NY,NX)
        ELSEIF(L.EQ.1)THEN
          TK0(L,NY,NX)=TKA(NY,NX)
        ELSE
          TK0(L,NY,NX)=TK0(L-1,NY,NX)
        ENDIF
        !       IF(L.EQ.1)THEN
        !         WRITE(*,7753)'TKM',I,J,M,NX,NY,L,TK0(L,NY,NX)
        !          2,THFLWW(L,NY,NX),TFLX0(L,NY,NX),HFLW0W(L,NY,NX)
        !          3,HFLWRLW,HFLWLW
        !          3,VOLS0(L,NY,NX),VOLW0(L,NY,NX),VOLI0(L,NY,NX),VOLS1(L,NY,NX)
        !          2,WFLXS(L,NY,NX),WFLXI(L,NY,NX)
        !          3,TFLWS(L,NY,NX),TFLWW(L,NY,NX),TFLWI(L,NY,NX)
        !          3,XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
        !          4,XHFLWW(L,NY,NX),HFLSW(L,NY,NX)
        !          5,VHCPWM(M+1,L,NY,NX)
!7753  FORMAT(A8,6I4,80E14.6)
        !       ENDIF
9780  CONTINUE
!      if(curday>=176)then
!        write(*,*)'line',__LINE__,'tk1',TK1(8,ny,nx),TK1(9,ny,nx),M
!      endif
      !
      !     SNOW RUNOFF
      !
      !     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
      !
      VOLS0(1,NY,NX)=VOLS0(1,NY,NX)+TQS1(NY,NX)
      VOLW0(1,NY,NX)=VOLW0(1,NY,NX)+TQW1(NY,NX)
      VOLI0(1,NY,NX)=VOLI0(1,NY,NX)+TQI1(NY,NX)
      ENGY0=VHCPWM(M+1,1,NY,NX)*TK0(1,NY,NX)
      VHCPWM(M+1,1,NY,NX)=cps*VOLS0(1,NY,NX)+cpw*VOLW0(1,NY,NX)+cpi*VOLI0(1,NY,NX)
      IF(VHCPWM(M+1,1,NY,NX).GT.VHCPWX(NY,NX))THEN
        TK0(1,NY,NX)=(ENGY0+THQS1(NY,NX))/VHCPWM(M+1,1,NY,NX)
      ELSE
        TK0(1,NY,NX)=TKA(NY,NX)
      ENDIF
!
      !     IF SNOWPACK DISAPPEARS ALL MATERIAL,HEAT TRANSFERRED TO SOIL SURFACE
      !
      !     VHCP1,VHCP1A,VHCP1P=total soil,soil+micropore,macropore heat capacity
      !     TK1=soil surface temperature, why not to litter layer
      !
      IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX).AND.TKA(NY,NX).GT.TFice)THEN
        FLWS=VOLS0(1,NY,NX)
        FLWW=VOLW0(1,NY,NX)
        FLWI=VOLI0(1,NY,NX)
        HFLWS=(cpw*FLWW+cps*FLWS+cpi*FLWI)*TK0(1,NY,NX)
        VOLS0(1,NY,NX)=VOLS0(1,NY,NX)-FLWS
        VOLW0(1,NY,NX)=VOLW0(1,NY,NX)-FLWW
        VOLI0(1,NY,NX)=VOLI0(1,NY,NX)-FLWI
        VOLW1(NUM(NY,NX),NY,NX)=VOLW1(NUM(NY,NX),NY,NX)+FLWW
        VOLI1(NUM(NY,NX),NY,NX)=VOLI1(NUM(NY,NX),NY,NX)+FLWI+FLWS/DENSI
        ENGY1=VHCP1(NUM(NY,NX),NY,NX)*TK1(NUM(NY,NX),NY,NX)
        VHCP1(NUM(NY,NX),NY,NX)=VHCM(NUM(NY,NX),NY,NX) &
          +cpw*(VOLW1(NUM(NY,NX),NY,NX)+VOLWH1(NUM(NY,NX),NY,NX)) &
          +cpi*(VOLI1(NUM(NY,NX),NY,NX)+VOLIH1(NUM(NY,NX),NY,NX))
        VHCP1A(NUM(NY,NX),NY,NX)=VHCM(NUM(NY,NX),NY,NX) &
          +cpw*VOLW1(NUM(NY,NX),NY,NX)+cpi*VOLI1(NUM(NY,NX),NY,NX)
        VHCP1B(NUM(NY,NX),NY,NX)=cpw*VOLWH1(NUM(NY,NX),NY,NX) &
          +cpi*VOLIH1(NUM(NY,NX),NY,NX)
        IF(VHCP1(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
          TK1(NUM(NY,NX),NY,NX)=(ENGY1+HFLWS)/VHCP1(NUM(NY,NX),NY,NX)
        ELSE
          TK1(NUM(NY,NX),NY,NX)=TKA(NY,NX)
        ENDIF
      ENDIF
!      if(curday>=176)then
!        write(*,*)'line',__LINE__,'tk1',TK1(8,ny,nx),TK1(9,ny,nx),NUM(NY,NX),NL(NY,NX)
!      endif
      !
      ! SURFACE RESIDUE WATER AND TEMPERATURE
      !
      ! XVOLT,XVOLW=free water+ice,water in litter layer
      ! VOLWM,VOLPM=surface water,air content for use in trnsfr.f
      ! VOLWRX=maximum water retention by litter
      ! VHCP1=volumetric heat capacity of litter
      ! VOLA1,VOLW1,VOLI1,VOLP1=pore,water,ice,air volumes of litter
      ! VOLWRX=maximum water retention by litter
      ! TFLXR,WFLXR=litter water,latent heat flux from freeze-thaw
      ! VOLR=dry litter volume
      ! THETWX,THETIX,THETPX=water,ice,air concentrations
      ! VHCP1=volumetric heat capacity of litter
      ! TK1=litter temperature
      ! HFLWRL,TFLXR,THQR1=litter total cond+conv,latent,runoff heat flux
      !
      ! if(curday>=41)then
      !   write(*,*)'VOLW1(0,NY,NX)=',VOLW1(0,NY,NX),FLWRL(NY,NX),WFLXR(NY,NX),TQR1(NY,NX)
      !   write(*,*)'VOLI1(0,NY,NX)=',VOLI1(0,NY,NX),WFLXR(NY,NX)
      !   write(*,*)'at line',__LINE__
      ! endif
      VOLW1(0,NY,NX)=AMAX1(0.0_r8,VOLW1(0,NY,NX)+FLWRL(NY,NX)+WFLXR(NY,NX)+TQR1(NY,NX))
      VOLI1(0,NY,NX)=AMAX1(0.0_r8,VOLI1(0,NY,NX)-WFLXR(NY,NX)/DENSI)
      VOLP1(0,NY,NX)=AMAX1(0.0_r8,VOLA1(0,NY,NX)-VOLW1(0,NY,NX)-VOLI1(0,NY,NX))
      VOLWM(M+1,0,NY,NX)=VOLW1(0,NY,NX)
      VOLPM(M+1,0,NY,NX)=VOLP1(0,NY,NX)
      TVOLWI=VOLW1(0,NY,NX)+VOLI1(0,NY,NX)
      XVOLT(NY,NX)=AMAX1(0.0_r8,TVOLWI-VOLWRX(NY,NX))
      IF(TVOLWI.GT.ZEROS(NY,NX))THEN
        VOLWRZ=VOLW1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
        VOLIRZ=VOLI1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
        XVOLW(NY,NX)=AMAX1(0.0_r8,VOLW1(0,NY,NX)-VOLWRZ)
        XVOLI(NY,NX)=AMAX1(0.0_r8,VOLI1(0,NY,NX)-VOLIRZ)
      ELSE
        XVOLW(NY,NX)=0.0_r8
        XVOLI(NY,NX)=0.0_r8
      ENDIF
      XVOLTM(M+1,NY,NX)=XVOLT(NY,NX)
      XVOLWM(M+1,NY,NX)=XVOLW(NY,NX)
      XVOLIM(M+1,NY,NX)=XVOLI(NY,NX)
      IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
        THETWX(0,NY,NX)=AMAX1(0.0_r8,VOLW1(0,NY,NX)/VOLR(NY,NX))
        THETIX(0,NY,NX)=AMAX1(0.0_r8,VOLI1(0,NY,NX)/VOLR(NY,NX))
        THETPX(0,NY,NX)=AMAX1(0.0_r8,VOLP1(0,NY,NX)/VOLR(NY,NX)) &
          *AMAX1(0.0_r8,(1.0_r8-XVOLT(NY,NX)/VOLWD(NY,NX)))
      ELSE
        THETWX(0,NY,NX)=0.0_r8
        THETIX(0,NY,NX)=0.0_r8
        THETPX(0,NY,NX)=1.0_r8
      ENDIF
      THETPM(M+1,0,NY,NX)=THETPX(0,NY,NX)
      VHCPXX=VHCP1(0,NY,NX)              !heat capacity
      TK0XX=TK1(0,NY,NX)                 !residual temperature
      ENGYR=VHCP1(0,NY,NX)*TK1(0,NY,NX)  !initial energy content
      VHCP1(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW1(0,NY,NX)+cpi*VOLI1(0,NY,NX)  !update heat capacity
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
!        if(curday>=176)write(*,*)trim(mod_filename)//'at line',__LINE__,'tk',TK1(NUM(NY,NX),NY,NX)
        tk1pres=TK1(0,NY,NX)
        TK1(0,NY,NX)=(ENGYR+HFLWRL(NY,NX)+TFLXR(NY,NX) &
          +THQR1(NY,NX))/VHCP1(0,NY,NX)
        Ls=NUM(NY,NX)
!        if(curday>=176)then
!          write(*,*)trim(mod_filename)//'at line',__LINE__,'tk',TK1(0,NY,NX)
!          write(*,*)'ENGYR,HFLWRL(NY,NX),TFLXR(NY,NX),THQR1(NY,NX),VHCP1(0,NY,NX)'
!          write(*,*)ENGYR,HFLWRL(NY,NX),TFLXR(NY,NX),THQR1(NY,NX),VHCP1(0,NY,NX)
!          write(*,*)VHCM(Ls,NY,NX),VOLW1(Ls,NY,NX)+VOLWH1(Ls,NY,NX),VOLI1(Ls,NY,NX)+VOLIH1(Ls,NY,NX)
!          write(*,*)BKDS(Ls,NY,NX),VHCPRX(NY,NX),VHCP1(0,NY,NX)
!        endif
        if(abs(TK1(0,NY,NX)-TK1(Ls,NY,NX))>15._r8)TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
          !too great difference between temepratures in residual layer and top soil layer
!          VHCP1(Ls,NY,NX)=VHCM(Ls,NY,NX)+cpw*(VOLW1(Ls,NY,NX) &
!            +VOLWH1(Ls,NY,NX))+cpi*(VOLI1(Ls,NY,NX)+VOLIH1(Ls,NY,NX))
!          TK1(0,NY,NX)=(TK1(0,NY,NX)*VHCP1(0,NY,NX)+TK1(Ls,NY,NX)*VHCP1(Ls,NY,NX)) &
!            /(VHCP1(0,NY,NX)+VHCP1(1,NY,NX))
!          TK1(Ls,NY,NX)=TK1(0,NY,NX)
!          if(TK1(0,NY,NX)<0._r8)then
!            write(*,*)'at line',__LINE__,'TK1',TK1(0,NY,NX)
!          endif
!          write(*,*)
!          write(*,*)'ENGYR+HFLWRL(NY,NX)+TFLXR(NY,NX)+THQR1(NY,NX),VHCP1(0,NY,NX),ENGYR=',&
!            (ENGYR+HFLWRL(NY,NX)+TFLXR(NY,NX)+THQR1(NY,NX)),VHCP1(0,NY,NX)
!          write(*,*)'HFLWRL(NY,NX),TFLXR(NY,NX),THQR1(NY,NX)=',HFLWRL(NY,NX),TFLXR(NY,NX),THQR1(NY,NX)
!          write(*,*)'tk1pres, TK1(0,NY,NX)=',tk1pres,TK1(0,NY,NX)
!          write(*,*)'ORGC(0,NY,NX),VOLW1(0,NY,NX),VOLI1(0,NY,NX)=',cpo*ORGC(0,NY,NX),cpw*VOLW1(0,NY,NX) &
!            ,cpi*VOLI1(0,NY,NX)
!          write(*,*)'at line',__LINE__
!          if(TK1(0,NY,NX)<200._r8)call endrun(trim(mod_filename)//'at line',__LINE__)
!        endif
      ELSE
        TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
      ENDIF
!      if(curday>=176)write(*,*)trim(mod_filename)//'at line',__LINE__,'tk',TK1(NUM(NY,NX),NY,NX)

      ! if(curday>=377)then
      !   write(*,*)'TK1(0,NY,NX)=',TK1(0,NY,NX)
      !   write(*,*)'VHCP1(0,NY,NX).GT.VHCPRX(NY,NX)',&
      !    VHCP1(0,NY,NX).GT.VHCPRX(NY,NX)
      ! endif
      ! IF(I.GT.350.AND.NX.EQ.1)THEN
      !   WRITE(*,7754)'VOLW0',I,J,M,NX,NY,NUM(NY,NX)
      !  2,VOLW1(0,NY,NX),VOLI1(0,NY,NX),VOLP1(0,NY,NX)
      !  3,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPX(0,NY,NX)
      !  2,FLWRL(NY,NX),WFLXR(NY,NX),TQR1(NY,NX),EVAPR(NY,NX)
      !  4,XVOLW(NY,NX),XVOLI(NY,NX),XVOLT(NY,NX)
      !  5,TVOLWI,VOLWRX(NY,NX)
      !  5,FLWRLG,FLWRLW,FLYM,FLQR
      !  3,TK1(0,NY,NX),HFLWRL(NY,NX),TFLXR(NY,NX)
      !  3,THQR1(NY,NX),ENGYR,VHCP1(0,NY,NX),VHCPRX(NY,NX)
      !  2,HFLWRLG,HFLWRLW,HWFLYM,HFLXR
      !  3,VOLPM(M,0,NY,NX),VOLPM(M,NUM(NY,NX),NY,NX)
      !  4,VOLI1(NUM(NY,NX),NY,NX),TK1(NUM(NY,NX),NY,NX)
      !7754  FORMAT(A8,6I4,60E12.4)
      !       ENDIF
      !
      ! SOIL LAYER WATER, ICE AND TEMPERATURE
      !
      ! VOLW1,VOLI1=micropore water,ice volume
      ! VOLWX1=micropore water volume behind wetting front
      ! VOLWH1,VOLIH1=macropore water,ice volume
      ! TFLWL=net water flux
      ! FINHL=micropore-macropore flux
      ! TWFLXL,TWFLXH=total accumulated freeze-thaw in micropores,macropores
      ! FLU1=subsurface water input
      ! DENSI=ice density
      ! VOLA1,VOLAH1=micropore,macropore volume
      ! VOLP1,VOLPH1=micropore,macropore air volume
      ! VOLWM,VOLWHM,VOLPM,FLPM=micropore,macropore water volume, air volume
      ! and change in air volume for use in trnsfr.f
      ! THETWX,THETIX,THETPX,THETPY=bulk water,ice,air concn,air-filled porosity
      ! THETPM=air concentration for use in trnsfr.f
      ! FMAC,FGRD=macropore,micropore fraction
      ! CNDH1=maropore hydraulic conductivity
      ! VHCP1,VHCM=volumetric heat capacities of total volume, solid
      ! VHCP1A,VHCP1B=volumetric heat capacities of soil+micropore,macropore
      ! TK1=soil temperature
      !
      !if(curday>=176)then
      !  write(*,*)'line',__LINE__,'tk1',TK1(8,ny,nx),TK1(9,ny,nx),M
      !endif
      DO 9785 L=NUM(NY,NX),NL(NY,NX)
!        if(curday>=176)then
!          write(*,*)'at line',__LINE__,'tk',TK1(L,NY,NX),L
!        endif
        IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VOLW1(L,NY,NX)=VOLW1(L,NY,NX)+TFLWL(L,NY,NX)+FINHL(L,NY,NX)+TWFLXL(L,NY,NX)+FLU1(L,NY,NX)
          if(abs(VOLW1(L,NY,NX))>1.e20_r8)then
            write(*,*)'VOLW1(L,NY,NX)=',VOLW1(L,NY,NX),L
            write(*,*)'TFLWL=',TFLWL(L,NY,NX)
            write(*,*)'FINHL=',FINHL(L,NY,NX)
            write(*,*)'TWFLXL=',TWFLXL(L,NY,NX)+FLU1(L,NY,NX)
            write(*,*)'at line',__LINE__
            call endrun(trim(mod_filename)//'at line',__LINE__)
          endif
          VOLWX1(L,NY,NX)=VOLWX1(L,NY,NX)+TFLWLX(L,NY,NX)+FINHL(L,NY,NX)+TWFLXL(L,NY,NX)+FLU1(L,NY,NX)
          VOLWX1(L,NY,NX)=AMIN1(VOLW1(L,NY,NX),VOLWX1(L,NY,NX))
          VOLI1(L,NY,NX)=VOLI1(L,NY,NX)-TWFLXL(L,NY,NX)/DENSI
          VOLWH1(L,NY,NX)=VOLWH1(L,NY,NX)+TFLWHL(L,NY,NX)-FINHL(L,NY,NX)+TWFLXH(L,NY,NX)
          VOLIH1(L,NY,NX)=VOLIH1(L,NY,NX)-TWFLXH(L,NY,NX)/DENSI
          IF(BKDS(L,NY,NX).GT.ZERO)THEN
            VOLP1Z(L,NY,NX)=VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
            if(abs(VOLP1Z(L,NY,NX))>1.e20_r8)then
              write(*,*)'VOLA1=',VOLA1(L,NY,NX),L
              write(*,*)'VOLW1=',VOLW1(L,NY,NX)
              write(*,*)'VOLI1=',VOLI1(L,NY,NX)
            endif
            VOLP1(L,NY,NX)=AMAX1(0.0,VOLP1Z(L,NY,NX))
            VOLPH1Z(L,NY,NX)=VOLAH1(L,NY,NX)-VOLWH1(L,NY,NX)-VOLIH1(L,NY,NX)
            VOLPH1(L,NY,NX)=AMAX1(0.0,VOLPH1Z(L,NY,NX))
            VOLAH1(L,NY,NX)=AMAX1(0.0,VOLAH(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX) &
              *(safe_adb(VOLW1(L,NY,NX),VOLY(L,NY,NX))-WP(L,NY,NX))*VOLT(L,NY,NX))
          ELSE
            VOLP1Z(L,NY,NX)=0.0
            VOLP1(L,NY,NX)=0.0
            VOLPH1Z(L,NY,NX)=0.0
            VOLPH1(L,NY,NX)=0.0
            VOLA1(L,NY,NX)=VOLW1(L,NY,NX)+VOLI1(L,NY,NX)
            VOLAH1(L,NY,NX)=0.0
          ENDIF
          VOLWM(M+1,L,NY,NX)=VOLW1(L,NY,NX)
          VOLWHM(M+1,L,NY,NX)=VOLWH1(L,NY,NX)
          VOLPM(M+1,L,NY,NX)=VOLP1(L,NY,NX)+VOLPH1(L,NY,NX)+THETPI*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
          FLPM(M,L,NY,NX)=VOLPM(M,L,NY,NX)-VOLPM(M+1,L,NY,NX)
          VOLTX=VOLY(L,NY,NX)+VOLAH1(L,NY,NX)
          THETWX(L,NY,NX)=AMAX1(0.0,(VOLW1(L,NY,NX)+VOLWH1(L,NY,NX))/VOLTX)
          THETIX(L,NY,NX)=AMAX1(0.0,(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))/VOLTX)
          THETPX(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))/VOLTX)
          THETPM(M+1,L,NY,NX)=THETPX(L,NY,NX)
          IF(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            THETPY(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX)) &
              /(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)))
!              if(curday>=41)then
!                write(*,*)'aTHETPY=',THETPY(L,NY,NX),L,VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)
!                write(*,*)'VOLP1 =',VOLP1(L,NY,NX),VOLPH1(L,NY,NX)
!              endif
          ELSE
            THETPY(L,NY,NX)=0.0_r8
          ENDIF
          IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            FMAC(L,NY,NX)=FHOL(L,NY,NX)*VOLAH1(L,NY,NX)/VOLAH(L,NY,NX)
            CNDH1(L,NY,NX)=CNDH(L,NY,NX)*(VOLAH1(L,NY,NX)/VOLAH(L,NY,NX))**2
          ELSE
            FMAC(L,NY,NX)=0.0_r8
            CNDH1(L,NY,NX)=0.0_r8
          ENDIF
          FGRD(L,NY,NX)=1.0_r8-FMAC(L,NY,NX)
          TKXX=TK1(L,NY,NX)
          VHXX=VHCP1(L,NY,NX)
          ENGY1=VHCP1(L,NY,NX)*TK1(L,NY,NX)
          if(TK1(L,NY,NX)>1.e3_r8.or.TK1(L,NY,NX)<0._r8)then
            write(*,*)'L=',L,NY,NX,NUM(NY,NX)
            write(*,*)'BKDS(L,NY,NX)=',BKDS(L,NY,NX)
            write(*,*)'VHCP1(L,NY,NX),TK1(L,NY,NX)',L,VHCP1(L,NY,NX),TK1(L,NY,NX)
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          VHCP1(L,NY,NX)=VHCM(L,NY,NX)+cpw*(VOLW1(L,NY,NX)&
            +VOLWH1(L,NY,NX))+cpi*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
          VHCP1A(L,NY,NX)=VHCM(L,NY,NX)+cpw*VOLW1(L,NY,NX)+cpi*VOLI1(L,NY,NX)
          VHCP1B(L,NY,NX)=cpw*VOLWH1(L,NY,NX)+cpi*VOLIH1(L,NY,NX)
          !
          !         BEGIN ARTIFICIAL SOIL WARMING
          !
          !         THFLWL=THFLWL incremented for soil warming
          !         TKSZ=temperature used to calculate additional heat flux for warming
          !
          !         IF(NX.EQ.3.AND.NY.EQ.2.AND.L.GT.NUM(NY,NX)
          !           3.AND.L.LE.17.AND.I.GE.152.AND.I.LE.304)THEN
          !           THFLWL(L,NY,NX)=THFLWL(L,NY,NX)
          !             2+(TKSZ(I,J,L)-TK1(L,NY,NX))*VHCP1(L,NY,NX)*XNPH
          !             WRITE(*,3379)'TKSZ',I,J,M,NX,NY,L,TKSZ(I,J,L)
          !               2,TK1(L,NY,NX),VHCP1(L,NY,NX),THFLWL(L,NY,NX)
          !3379  FORMAT(A8,6I4,12E12.4)
          !           ENDIF
          !
          !           END ARTIFICIAL SOIL WARMING
!
          IF(VHCP1(L,NY,NX).GT.ZEROS(NY,NX))THEN
            TK1(L,NY,NX)=(ENGY1+THFLWL(L,NY,NX)+TTFLXL(L,NY,NX) &
              +HWFLU1(L,NY,NX))/VHCP1(L,NY,NX)
            if(abs(TK1(L,NY,NX))>1.e3_r8)then
              write(*,*)'ENGY1+THFLWL(L,NY,NX)+TTFLXL(L,NY,NX)',&
                ENGY1,THFLWL(L,NY,NX),TTFLXL(L,NY,NX)
              write(*,*)'HWFLU1(L,NY,NX)),VHCP1(L,NY,NX)',&
                HWFLU1(L,NY,NX),VHCP1(L,NY,NX)
              call endrun(trim(mod_filename)//' at line',__LINE__)
            endif
          ELSEIF(L.EQ.1)THEN
            TK1(L,NY,NX)=TKA(NY,NX)
          ELSE
            TK1(L,NY,NX)=TK1(L-1,NY,NX)
          ENDIF
        ELSE
          VOLWM(M+1,L,NY,NX)=0.0
          VOLWHM(M+1,L,NY,NX)=0.0
          VOLPM(M+1,L,NY,NX)=0.0
          FLPM(M,L,NY,NX)=VOLPM(M,L,NY,NX)
          THETPM(M+1,L,NY,NX)=0.0
        ENDIF
        !         IF(L.LE.NUM(NY,NX)+1)THEN
        !           WRITE(*,3377)'VOLW1',I,J,M,NX,NY,L,N6X(NY,NX)
        !             2,VOLW1(L,NY,NX),VOLWX1(L,NY,NX)
        !             3,VOLI1(L,NY,NX),VOLA1(L,NY,NX),VOLP1(L,NY,NX)
        !             2,TFLWL(L,NY,NX),TFLWLX(L,NY,NX)
        !             3,TWFLXL(L,NY,NX),FINHL(L,NY,NX),FLU1(L,NY,NX)
        !             3,VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),VOLAH1(L,NY,NX),VOLPH1(L,NY,NX)
        !             4,VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
        !             5,VOLPM(M,L,NY,NX),VOLPM(M+1,L,NY,NX)
        !             5,PSISM1(L,NY,NX),THETPX(L,NY,NX)
        !             6,FLWL(3,L,NY,NX),FLWL(3,L+1,NY,NX)
        !             7,FLWL(2,L,NY,NX),FLWL(2,L,NY+1,NX)
        !             8,FLWL(1,L,NY,NX),FLWL(1,L,NY,NX+1)
        !             6,FLWLX(3,L,NY,NX),FLWLX(3,L+1,NY,NX)
        !             7,FLWLX(2,L,NY,NX),FLWLX(2,L,NY+1,NX)
        !             8,FLWLX(1,L,NY,NX),FLWLX(1,L,NY,NX+1)
        !             6,FLW(3,L,NY,NX),FLW(3,L+1,NY,NX)
        !             7,FLW(2,L,NY,NX),FLW(2,L,NY+1,NX)
        !             8,FLW(1,L,NY,NX),FLW(1,L,NY,NX+1)
        !             9,WFLXL(L,NY,NX),THAW(L,NY,NX)
        !             9,FLPM(M,L,NY,NX),FLSW(L,NY,NX)
        !           WRITE(*,3377)'VOLWH1',I,J,M,NX,NY,L,N6X(NY,NX),VOLWH1(L,NY,NX)
        !             2,TFLWHL(L,NY,NX),FINHL(L,NY,NX),VOLIH1(L,NY,NX)
        !             4,TWFLXH(L,NY,NX),TQR1(NY,NX),VOLPH1(L,NY,NX)
        !             6,FLWHL(3,L,NY,NX),FLWHL(3,L+1,NY,NX)
        !             7,FLWHL(2,L,NY,NX),FLWHL(2,L,NY+1,NX)
        !             8,FLWHL(1,L,NY,NX),FLWHL(1,L,NY,NX+1)
        !             3,VHCM(L,NY,NX),VOLW1(L,NY,NX),VOLWH1(L,NY,NX),VOLI1(L,NY,NX)
        !             4,THETW(L,NY,NX),THETI(L,NY,NX),FINHL(L,NY,NX),THQR1(NY,NX)
        !             5,HFLWL(3,L,NY,NX),HFLWL(3,N6X(NY,NX),NY,NX)
        !             6,HFLWL(1,L,NY,NX),HFLWL(1,L,NY,NX+1)
        !           WRITE(*,3377)'TK1',I,J,M,NX,NY,L,N6X(NY,NX),TK1(L,NY,NX)
        !             2,THFLWL(L,NY,NX),TTFLXL(L,NY,NX),HWFLU1(L,NY,NX)
        !             3,VHCP1(L,NY,NX),HFLWL(3,L,NY,NX),HFLWL(3,L+1,NY,NX)
        !             4,HFLWLW,HFLWLG,TKXX,VHXX,ENGY1,TK1(0,NY,NX)
        !3377  FORMAT(A8,7I4,40E12.4)
        !         ENDIF

9785  CONTINUE
      !
      !       RESET SURFACE LAYER NUMBER AND TRANSFER ALL WATER TO SOIL SURFACE LAYER
      !       IF LAKE SURFACE LAYER IS LOST TO EVAPORATION
      !
      !       NUM=new surface layer number after complete lake evaporation
      !       FLWNU,FLWHNU,HFLWNU=lake surface water flux, heat flux if lake surface disappears
!
      IF(BKDS(NUM(NY,NX),NY,NX).LE.ZERO &
        .AND.VHCP1(NUM(NY,NX),NY,NX).LE.VHCPNX(NY,NX))THEN
        NUX=NUM(NY,NX)
        DO 9970 LL=NUX+1,NL(NY,NX)
          IF(VOLX(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
            NUM(NY,NX)=LL
            FLWNX(NY,NX)=FLW(3,NUM(NY,NX),NY,NX)
            FLWXNX(NY,NX)=FLWX(3,NUM(NY,NX),NY,NX)
            FLWHNX(NY,NX)=FLWH(3,NUM(NY,NX),NY,NX)
            HFLWNX(NY,NX)=HFLW(3,NUM(NY,NX),NY,NX)
            ! WRITE(*,5598)'SURFM',I,J,M,NX,NY,LL,NUX,NUM(NY,NX)
            ! 2,FLWNX(NY,NX)
            ! 2,VOLW1(NUX,NY,NX),VOLW1(NUM(NY,NX),NY,NX)
            ! 2,VHCP1(NUX,NY,NX),VHCP1(NUM(NY,NX),NY,NX)
            ! 2,TK1(NUX,NY,NX),TK1(NUM(NY,NX),NY,NX)
            ! 3,FLW(3,NUX,NY,NX),FLW(3,NUM(NY,NX),NY,NX)
            ! 3,HFLW(3,NUX,NY,NX),HFLW(3,NUM(NY,NX),NY,NX)
            ! 4,VHCPNX(NY,NX),VHCP1(0,NY,NX)
            !5598          FORMAT(A8,8I4,20E12.4)
            GO TO 9971
          ENDIF
9970    CONTINUE
      ENDIF
9971  CONTINUE
9790  CONTINUE
9795  CONTINUE

  end subroutine UpdateStateSolutionNotNPH
!------------------------------------------------------------------------------------------

  subroutine NoSurfResidueForEnergyBalance(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
! begin_execution

  TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
!  if(curday>=377)then
!    write(*,*)'TK1(0,NY,NX)=',TK1(0,NY,NX)
!    write(*,*)'NoSurfResidueForEnergyBalanc'
!  endif
  EVAPR(NY,NX)=0.0_r8
  RFLXR=0.0_r8
  EFLXR=0.0_r8
  VFLXR=0.0_r8
  SFLXR=0.0_r8
  HFLXR=0.0_r8
  FLV1=0.0_r8
  HWFLV1=0.0_r8
  HFLCR1=0.0_r8
  THRMZ=0.0_r8
  end subroutine NoSurfResidueForEnergyBalance
!------------------------------------------------------------------------------------------

  subroutine SnowCoveredTopSoilFlux(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

! begin_execution
  RFLXG=0.0_r8
  EFLXG=0.0_r8
  VFLXG=0.0_r8
  SFLXG=0.0_r8
  HFLXG=0.0_r8
  RFLXR=0.0_r8
  EFLXR=0.0_r8
  VFLXR=0.0_r8
  SFLXR=0.0_r8
  HFLXR=0.0_r8
  FLWLG=0.0_r8
  FLWLXG=0.0_r8
  FLWHLG=0.0_r8
  HFLWLG=0.0_r8
  FLWRLG=0.0_r8
  HFLWRLG=0.0_r8
  FLWVLS=0.0_r8
  EVAPG(NY,NX)=0.0_r8
  EVAPR(NY,NX)=0.0_r8
  end subroutine SnowCoveredTopSoilFlux
!------------------------------------------------------------------------------------------

  subroutine AccumWaterVaporHeatFluxes(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX
! begin_execution
! HOURLY-ACCUMULATED WATER, VAPOR AND HEAT FLUXES THROUGH
! SURFACE RESIDUE AND SOIL SURFACE
!
! THAWR,HTHAWR=litter water,heat fluxes from freeze-thaw
! FLW,FLWH,HFLW=soil micropore,macropore,heat fluxes
! FLWR,HFLWR=litter water,heat fluxes
! FLSW,FLSWH=water from snowpack to soil micropores,macropores
! HFLSW=convective heat from snowpack to soil
! FLSWR=water flux from snowpack to litter
! HFLSWR=convective heat flux from snowpack to litter
! HEATI,HEATE,HEATS,HEATG=net radiation,latent,sensible,storage heat
! TEVAPG=total evaporation
! FLWM,FLWHM=water flux into soil micropore,macropore for use in trnsfr.f
! VOLWX1=VOLW1 accounting for wetting front
!
  THAWR(NY,NX)=THAWR(NY,NX)+WFLXR(NY,NX)
  HTHAWR(NY,NX)=HTHAWR(NY,NX)+TFLXR(NY,NX)
  FLW(3,NUM(NY,NX),NY,NX)=FLW(3,NUM(NY,NX),NY,NX) &
    +FLWL(3,NUM(NY,NX),NY,NX)
  FLWX(3,NUM(NY,NX),NY,NX)=FLWX(3,NUM(NY,NX),NY,NX) &
    +FLWLX(3,NUM(NY,NX),NY,NX)
  FLWH(3,NUM(NY,NX),NY,NX)=FLWH(3,NUM(NY,NX),NY,NX) &
    +FLWHL(3,NUM(NY,NX),NY,NX)
  HFLW(3,NUM(NY,NX),NY,NX)=HFLW(3,NUM(NY,NX),NY,NX) &
    +HFLWL(3,NUM(NY,NX),NY,NX)
  FLWR(NY,NX)=FLWR(NY,NX)+FLWRL(NY,NX)
  HFLWR(NY,NX)=HFLWR(NY,NX)+HFLWRL(NY,NX)
!  if(HFLWR(NY,NX)<-1.4_r8)then
!    write(*,*)
!    write(*,*)'******HFLWR(NY,NX)+HFLWRL(NY,NX)',HFLWR(NY,NX),HFLWRL(NY,NX)
!    call print_info(trim(mod_filename)//' at line',__LINE__)
!  endif
  HEATI(NY,NX)=HEATI(NY,NX)+RFLXG+RFLXR+RFLXW
  HEATS(NY,NX)=HEATS(NY,NX)+SFLXG+SFLXR+SFLXW
  HEATE(NY,NX)=HEATE(NY,NX)+EFLXG+EFLXR+EFLXW
  HEATV(NY,NX)=HEATV(NY,NX)+VFLXG+VFLXR+VFLXW
  HEATH(NY,NX)=HEATH(NY,NX)+RFLXG+RFLXR+RFLXW &
    +SFLXG+SFLXR+SFLXW+EFLXG+EFLXR+EFLXW+VFLXG+VFLXR+VFLXW
  TEVAPG(NY,NX)=TEVAPG(NY,NX)+EVAPG(NY,NX)+EVAPR(NY,NX) &
    +EVAPS(NY,NX)+EVAPW(NY,NX)
!  if(TEVAPG(NY,NX)/=TEVAPG(NY,NX))then
!    write(*,*)'EVAPG=',EVAPG(NY,NX)
!    write(*,*)'EVAPR=',EVAPR(NY,NX)
!    write(*,*)'EVAPS=',EVAPS(NY,NX)
!    write(*,*)'EVAPW=',EVAPW(NY,NX)
!    call endrun(msg='NaN encounterd in '//mod_filename,line=__LINE__)
!  endif
  FLWM(M,3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)
  FLWHM(M,3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)
  end subroutine AccumWaterVaporHeatFluxes
!------------------------------------------------------------------------------------------

  subroutine UpdateFLuxAtNPH(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
! begin_execution
  DO 9695 NX=NHW,NHE
    DO 9690 NY=NVN,NVS
      IF(NUM(NY,NX).EQ.NU(NY,NX))THEN
        FLWNU(NY,NX)=FLW(3,N6X(NY,NX),NY,NX)
        FLWXNU(NY,NX)=FLWX(3,N6X(NY,NX),NY,NX)
        FLWHNU(NY,NX)=FLWH(3,N6X(NY,NX),NY,NX)
        HFLWNU(NY,NX)=HFLW(3,N6X(NY,NX),NY,NX)
      ELSE
        FLWNU(NY,NX)=FLWNX(NY,NX)
        FLWXNU(NY,NX)=FLWXNX(NY,NX)
        FLWHNU(NY,NX)=FLWHNX(NY,NX)
        HFLWNU(NY,NX)=HFLWNX(NY,NX)
      ENDIF
9690  CONTINUE
9695  CONTINUE
  end subroutine UpdateFLuxAtNPH
!------------------------------------------------------------------------------------------

  subroutine ExposedSoilFlux(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX

! begin_execution
! Watch out for L, is its value defined?
  call PrepSoilSurfaceEnerbyBalance(M,NY,NX)

! ENERGY BALANCE AT RESIDUE SURFACE
!
  IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
!
    call SoilSurfaceEnergyBalance(M,NY,NX)
!
!   IF NO SURFACE RESIDUE
!
  ELSE
    call NoSurfResidueForEnergyBalance(NY,NX)
  ENDIF
!  write(*,*)'SummaryAfterEnergyBalance'
  call SummaryAfterEnergyBalance(NY,NX)
  end subroutine ExposedSoilFlux
!------------------------------------------------------------------------------------------

  subroutine AtmosLandSurfaceExchange(M,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX

! begin_execution

  IF(VHCPWM(M,1,NY,NX).GT.VHCPWX(NY,NX))THEN
!   VHCPW,VHCPWX=current, minimum snowpack heat capacities
    call SolveSnowpack(M,NY,NX)
  ENDIF
!
! ENERGY EXCHANGE AT SOIL SURFACE IF EXPOSED UNDER SNOWPACK
! FSNW,FSNX=fractions of snow,snow-free cover
  IF(FSNX(NY,NX).GT.0.0_r8.AND.(BKDS(NUM(NY,NX),NY,NX).GT.ZERO &
    .OR.VHCP1(NUM(NY,NX),NY,NX).GT.VHCPNX(NY,NX)))THEN
!
    call ExposedSoilFlux(M,NY,NX)
!
  ELSE
!
    call SnowCoveredTopSoilFlux(NY,NX)
  ENDIF
!
! AGGREGATE RESIDUE AND SOIL SURFACE FLUXES BENEATH SNOW
! AND ATMOSPHERE
!
! FLWL,FLWLX=total water flux into soil micropores
! FLWHL=total water flux into soil macropores
! HFLWL=total heat flux into soil
! FLWRL,FLWLX=total water flux into litter
! HFLWRL=total heat flux into litter
! FLWV*=total internal vapor flux in soil
!
  FLWL(3,NUM(NY,NX),NY,NX)=FLWLW+FLWLG
  if(abs(FLWL(3,NUM(NY,NX),NY,NX))>1.e20_r8)then
    write(*,*)'FLWL(3,NUM(NY,NX),NY,NX)=',FLWLW,FLWLG
    write(*,*)'at line',__LINE__
    call endrun(trim(mod_filename)//'at line',__LINE__)
  endif
  FLWLX(3,NUM(NY,NX),NY,NX)=FLWLXW+FLWLXG
  FLWHL(3,NUM(NY,NX),NY,NX)=FLWHLW+FLWHLG
  HFLWL(3,NUM(NY,NX),NY,NX)=HFLWLW+HFLWLG
  FLWRL(NY,NX)=FLWRLW+FLWRLG
  HFLWRL(NY,NX)=HFLWRLW+HFLWRLG
!  if(curday>=41)then
!    write(*,*)'wwwwwwww5445HFLWRL(NY,NX)=',HFLWRLW,HFLWRLG
!    write(*,*)'at line',__LINE__
!  endif
! IF(I.GT.350.AND.NX.EQ.1)THEN
!   WRITE(*,7756)'FLWT',I,J,M,NX,NY,NUM(NY,NX)
!    2,FLWL(3,NUM(NY,NX),NY,NX),FLWLW,FLWLG
!    4,FLWHL(3,NUM(NY,NX),NY,NX),FLWHLW,FLWHLG
!    5,HFLWL(3,NUM(NY,NX),NY,NX),HFLWLW,HFLWLG
!    6,FLWRL(NY,NX),FLWRLW,FLWRLG
!    7,HFLWRL(NY,NX),HFLWRLW,HFLWRLG
!    9,HFLWQG,HFLVS1,HFLWS1,HFLVR1,HFLWR1
!    1,HFLWQR,HFLVSR,HFLWSR,HFLVR1,HFLWR1
!    8,FLQ0S(NY,NX),FLQ0W(NY,NX),FLQ0I(NY,NX),FLQM,FLHM,FLYM
!    9,PRECA(NY,NX)*XNPH,PRECW(NY,NX)*XNPH
!    7,RFLXW,EFLXW,SFLXW,-HFLXW+VFLXW,RFLXG,EFLXG,SFLXG,-HFLXG+VFLXG
!    7,RFLXR,EFLXR,SFLXR,-HFLXR+VFLXR
!    8,FLW0S(1,NY,NX),FLQ0S(NY,NX),EVAPS(NY,NX)
!    9,FLW0W(1,NY,NX),FLQ0W(NY,NX),EVAPW(NY,NX)
!    1,FLW0I(1,NY,NX),FLQ0I(NY,NX)
!    2,HFLW0W(1,NY,NX),HWFLQ0(NY,NX),HFLXW
!    3,XFLWS(1,NY,NX),XFLWW(1,NY,NX),XFLWI(1,NY,NX),XHFLWW(1,NY,NX)
!    8,FSNW(NY,NX),FSNX(NY,NX),CVRD(NY,NX),BARE(NY,NX)
!    9,THETPM(M,NUM(NY,NX),NY,NX)
!7756  FORMAT(A8,6I4,60E12.4)
! ENDIF
  end subroutine AtmosLandSurfaceExchange


end module WatsubMod
