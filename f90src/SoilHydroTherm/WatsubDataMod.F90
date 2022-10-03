module WatsubDataMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  VOLWX1(:,:,:)                      !
  real(r8),allocatable ::  VPQ(:,:)                           !
  real(r8),allocatable ::  TKQ(:,:)                           !
  real(r8),allocatable ::  XVOLT(:,:)                         !
  real(r8),allocatable ::  XVOLW(:,:)                         !
  real(r8),allocatable ::  XVOLI(:,:)                         !
  real(r8),allocatable ::  FMAC(:,:,:)                        !
  real(r8),allocatable ::  FGRD(:,:,:)                        !
  real(r8),allocatable ::  VOLW1(:,:,:)                       !
  real(r8),allocatable ::  VOLI1(:,:,:)                       !
  real(r8),allocatable ::  VHCP1(:,:,:)                       !
  real(r8),allocatable ::  VHCP1A(:,:,:)                      !
  real(r8),allocatable ::  VHCP1B(:,:,:)                      !
  real(r8),allocatable ::  TK1(:,:,:)                         !
  real(r8),allocatable ::  TWFLXL(:,:,:)                      !
  real(r8),allocatable ::  VOLW2(:,:,:)                       !
  real(r8),allocatable ::  VOLP1(:,:,:)                       !
  real(r8),allocatable ::  TWFLXH(:,:,:)                      !
  real(r8),allocatable ::  PRECM(:,:)                         !
  real(r8),allocatable ::  VOLS0(:,:,:)                       !
  real(r8),allocatable ::  VOLI0(:,:,:)                       !
  real(r8),allocatable ::  VOLW0(:,:,:)                       !
  real(r8),allocatable ::  VOLS1(:,:,:)                       !
  real(r8),allocatable ::  DLYRS0(:,:,:)                      !
  real(r8),allocatable ::  VOLP1Z(:,:,:)                      !
  real(r8),allocatable ::  TK0(:,:,:)                         !
  real(r8),allocatable ::  AREAU(:,:,:)                       !
  real(r8),allocatable ::  AREAUD(:,:,:)                      !
  real(r8),allocatable ::  FLQ0S(:,:)                         !
  real(r8),allocatable ::  FLQ0I(:,:)                         !
  real(r8),allocatable ::  FLQ0W(:,:)                         !
  real(r8),allocatable ::  FLQ1(:,:)                          !
  real(r8),allocatable ::  FLH1(:,:)                          !
  real(r8),allocatable ::  FLY1(:,:)                          !
  real(r8),allocatable ::  HWFLQ0(:,:)                        !
  real(r8),allocatable ::  HWFLQ1(:,:)                        !
  real(r8),allocatable ::  HWFLY1(:,:)                        !
  real(r8),allocatable ::  RAR(:,:)                           !
  real(r8),allocatable ::  RAGS(:,:)                          !
  real(r8),allocatable ::  BAREW(:,:)                         !
  real(r8),allocatable ::  CVRDW(:,:)                         !
  real(r8),allocatable ::  RARG(:,:)                          !
  real(r8),allocatable ::  RAGR(:,:)                          !
  real(r8),allocatable ::  RAGW(:,:)                          !
  real(r8),allocatable ::  PAREW(:,:)                         !
  real(r8),allocatable ::  PAREG(:,:)                         !
  real(r8),allocatable ::  RAG(:,:)                           !
  real(r8),allocatable ::  PARSW(:,:)                         !
  real(r8),allocatable ::  PARSG(:,:)                         !
  real(r8),allocatable ::  PARER(:,:)                         !
  real(r8),allocatable ::  PARSR(:,:)                         !
  real(r8),allocatable ::  QR1(:,:,:,:)                       !
  real(r8),allocatable ::  HQR1(:,:,:,:)                      !
  real(r8),allocatable ::  VOLPH1Z(:,:,:)                     !
  real(r8),allocatable ::  QS1(:,:,:)                         !
  real(r8),allocatable ::  QW1(:,:,:)                         !
  real(r8),allocatable ::  QI1(:,:,:)                         !
  real(r8),allocatable ::  HQS1(:,:,:)                        !
  real(r8),allocatable ::  TQR1(:,:)                          !
  real(r8),allocatable ::  THQR1(:,:)                         !
  real(r8),allocatable ::  EVAPG(:,:)                         !
  real(r8),allocatable ::  TTFLXL(:,:,:)                      !
  real(r8),allocatable ::  EVAPW(:,:)                         !
  real(r8),allocatable ::  EVAPS(:,:)                         !
  real(r8),allocatable ::  EVAPR(:,:)                         !
  real(r8),allocatable ::  FLWRL(:,:)                         !
  real(r8),allocatable ::  HFLWRL(:,:)                        !
  real(r8),allocatable ::  FINHL(:,:,:)                       !
  real(r8),allocatable ::  FLWL(:,:,:,:)                      !
  real(r8),allocatable ::  FLWHL(:,:,:,:)                     !
  real(r8),allocatable ::  HFLWL(:,:,:,:)                     !
  real(r8),allocatable ::  TFLWL(:,:,:)                       !
  real(r8),allocatable ::  TFLWHL(:,:,:)                      !
  real(r8),allocatable ::  THFLWL(:,:,:)                      !
  real(r8),allocatable ::  WFLXL(:,:,:)                       !
  real(r8),allocatable ::  TFLXL(:,:,:)                       !
  real(r8),allocatable ::  AVCNHL(:,:,:,:)                    !
  real(r8),allocatable ::  THRYW(:,:)                         !
  real(r8),allocatable ::  THRMW(:,:)                         !
  real(r8),allocatable ::  THRMS(:,:)                         !
  real(r8),allocatable ::  THRMR(:,:)                         !
  real(r8),allocatable ::  THRYG(:,:)                         !
  real(r8),allocatable ::  THRYR(:,:)                         !
  real(r8),allocatable ::  RADXW(:,:)                         !
  real(r8),allocatable ::  RADXG(:,:)                         !
  real(r8),allocatable ::  RADXR(:,:)                         !
  real(r8),allocatable ::  FLWLX(:,:,:,:)                     !
  real(r8),allocatable ::  TFLWLX(:,:,:)                      !
  real(r8),allocatable ::  FLU1(:,:,:)                        !
  real(r8),allocatable ::  HWFLU1(:,:,:)                      !
  real(r8),allocatable ::  PSISM1(:,:,:)                      !
  real(r8),allocatable ::  ALTG(:,:)                          !
  real(r8),allocatable ::  WFLXLH(:,:,:)                      !
  real(r8),allocatable ::  DLYRR(:,:)                         !
  real(r8),allocatable ::  WFLXR(:,:)                         !
  real(r8),allocatable ::  TFLXR(:,:)                         !
  real(r8),allocatable ::  HCNDR(:,:)                         !
  real(r8),allocatable ::  CNDH1(:,:,:)                       !
  real(r8),allocatable ::  VOLA1(:,:,:)                       !
  real(r8),allocatable ::  THETWX(:,:,:)                      !
  real(r8),allocatable ::  THETIX(:,:,:)                      !
  real(r8),allocatable ::  THETPX(:,:,:)                      !
  real(r8),allocatable ::  VOLAH1(:,:,:)                      !
  real(r8),allocatable ::  VOLWH1(:,:,:)                      !
  real(r8),allocatable ::  VOLPH1(:,:,:)                      !
  real(r8),allocatable ::  VOLIH1(:,:,:)                      !
  real(r8),allocatable ::  THETPY(:,:,:)                      !
  real(r8),allocatable ::  FLWNX(:,:)                         !
  real(r8),allocatable ::  FLWXNX(:,:)                        !
  real(r8),allocatable ::  FLWHNX(:,:)                        !
  real(r8),allocatable ::  HFLWNX(:,:)                        !
  real(r8),allocatable ::  PSISA1(:,:,:)                      !
  real(r8),allocatable ::  TQS1(:,:)                          !
  real(r8),allocatable ::  TQW1(:,:)                          !
  real(r8),allocatable ::  TQI1(:,:)                          !
  real(r8),allocatable ::  THQS1(:,:)                         !
  real(r8),allocatable ::  TFLWS(:,:,:)                       !
  real(r8),allocatable ::  TFLWW(:,:,:)                       !
  real(r8),allocatable ::  TFLWI(:,:,:)                       !
  real(r8),allocatable ::  THFLWW(:,:,:)                      !
  real(r8),allocatable ::  TFLX0(:,:,:)                       !
  real(r8),allocatable ::  WFLXS(:,:,:)                       !
  real(r8),allocatable ::  WFLXI(:,:,:)                       !
  real(r8),allocatable ::  FLW0W(:,:,:)                       !
  real(r8),allocatable ::  FLW0S(:,:,:)                       !
  real(r8),allocatable ::  FLW0I(:,:,:)                       !
  real(r8),allocatable ::  HFLW0W(:,:,:)                      !
  real(r8),allocatable ::  VOLS0M(:,:,:)                      !
  real(r8),allocatable ::  VOLW0M(:,:,:)                      !
  real(r8),allocatable ::  VOLI0M(:,:,:)                      !
  real(r8),allocatable ::  VHCPWMM(:,:,:)                     !
  real(r8),allocatable ::  TK0M(:,:,:)                        !
  integer, allocatable ::  N6X(:,:)

!----------------------------------------------------------------------
  public :: InitWatsubData
contains
  subroutine InitWatSubData

  implicit none

  allocate(N6X(JY,JX));         N6X=0
  allocate(VOLWX1(JZ,JY,JX));   VOLWX1=0._r8
  allocate(VPQ(JY,JX));         VPQ=0._r8
  allocate(TKQ(JY,JX));         TKQ=0._r8
  allocate(XVOLT(JY,JX));       XVOLT=0._r8
  allocate(XVOLW(JY,JX));       XVOLW=0._r8
  allocate(XVOLI(JY,JX));       XVOLI=0._r8
  allocate(FMAC(JZ,JY,JX));     FMAC=0._r8
  allocate(FGRD(JZ,JY,JX));     FGRD=0._r8
  allocate(VOLW1(0:JZ,JY,JX));  VOLW1=0._r8
  allocate(VOLI1(0:JZ,JY,JX));  VOLI1=0._r8
  allocate(VHCP1(0:JZ,JY,JX));  VHCP1=0._r8
  allocate(VHCP1A(JZ,JY,JX));   VHCP1A=0._r8
  allocate(VHCP1B(JZ,JY,JX));   VHCP1B=0._r8
  allocate(TK1(0:JZ,JY,JX));    TK1=0._r8
  allocate(TWFLXL(JZ,JY,JX));   TWFLXL=0._r8
  allocate(VOLW2(JZ,JY,JX));    VOLW2=0._r8
  allocate(VOLP1(0:JZ,JY,JX));  VOLP1=0._r8
  allocate(TWFLXH(JZ,JY,JX));   TWFLXH=0._r8
  allocate(PRECM(JY,JX));       PRECM=0._r8
  allocate(VOLS0(JS,JY,JX));    VOLS0=0._r8
  allocate(VOLI0(JS,JY,JX));    VOLI0=0._r8
  allocate(VOLW0(JS,JY,JX));    VOLW0=0._r8
  allocate(VOLS1(JS,JY,JX));    VOLS1=0._r8
  allocate(DLYRS0(JS,JY,JX));   DLYRS0=0._r8
  allocate(VOLP1Z(JZ,JY,JX));   VOLP1Z=0._r8
  allocate(TK0(JS,JY,JX));      TK0=0._r8
  allocate(AREAU(JZ,JY,JX));    AREAU=0._r8
  allocate(AREAUD(JZ,JY,JX));   AREAUD=0._r8
  allocate(FLQ0S(JY,JX));       FLQ0S=0._r8
  allocate(FLQ0I(JY,JX));       FLQ0I=0._r8
  allocate(FLQ0W(JY,JX));       FLQ0W=0._r8
  allocate(FLQ1(JY,JX));        FLQ1=0._r8
  allocate(FLH1(JY,JX));        FLH1=0._r8
  allocate(FLY1(JY,JX));        FLY1=0._r8
  allocate(HWFLQ0(JY,JX));      HWFLQ0=0._r8
  allocate(HWFLQ1(JY,JX));      HWFLQ1=0._r8
  allocate(HWFLY1(JY,JX));      HWFLY1=0._r8
  allocate(RAR(JY,JX));         RAR=0._r8
  allocate(RAGS(JY,JX));        RAGS=0._r8
  allocate(BAREW(JY,JX));       BAREW=0._r8
  allocate(CVRDW(JY,JX));       CVRDW=0._r8
  allocate(RARG(JY,JX));        RARG=0._r8
  allocate(RAGR(JY,JX));        RAGR=0._r8
  allocate(RAGW(JY,JX));        RAGW=0._r8
  allocate(PAREW(JY,JX));       PAREW=0._r8
  allocate(PAREG(JY,JX));       PAREG=0._r8
  allocate(RAG(JY,JX));         RAG=0._r8
  allocate(PARSW(JY,JX));       PARSW=0._r8
  allocate(PARSG(JY,JX));       PARSG=0._r8
  allocate(PARER(JY,JX));       PARER=0._r8
  allocate(PARSR(JY,JX));       PARSR=0._r8
  allocate(QR1(2,2,JV,JH));     QR1=0._r8
  allocate(HQR1(2,2,JV,JH));    HQR1=0._r8
  allocate(VOLPH1Z(JZ,JY,JX));  VOLPH1Z=0._r8
  allocate(QS1(2,JV,JH));       QS1=0._r8
  allocate(QW1(2,JV,JH));       QW1=0._r8
  allocate(QI1(2,JV,JH));       QI1=0._r8
  allocate(HQS1(2,JV,JH));      HQS1=0._r8
  allocate(TQR1(JY,JX));        TQR1=0._r8
  allocate(THQR1(JY,JX));       THQR1=0._r8
  allocate(EVAPG(JY,JX));       EVAPG=0._r8
  allocate(TTFLXL(JZ,JY,JX));   TTFLXL=0._r8
  allocate(EVAPW(JY,JX));       EVAPW=0._r8
  allocate(EVAPS(JY,JX));       EVAPS=0._r8
  allocate(EVAPR(JY,JX));       EVAPR=0._r8
  allocate(FLWRL(JY,JX));       FLWRL=0._r8
  allocate(HFLWRL(JY,JX));      HFLWRL=0._r8
  allocate(FINHL(JZ,JY,JX));    FINHL=0._r8
  allocate(FLWL(3,JD,JV,JH));   FLWL=0._r8
  allocate(FLWHL(3,JD,JV,JH));  FLWHL=0._r8
  allocate(HFLWL(3,JD,JV,JH));  HFLWL=0._r8
  allocate(TFLWL(JZ,JY,JX));    TFLWL=0._r8
  allocate(TFLWHL(JZ,JY,JX));   TFLWHL=0._r8
  allocate(THFLWL(JZ,JY,JX));   THFLWL=0._r8
  allocate(WFLXL(JZ,JY,JX));    WFLXL=0._r8
  allocate(TFLXL(JZ,JY,JX));    TFLXL=0._r8
  allocate(AVCNHL(3,JD,JV,JH)); AVCNHL=0._r8
  allocate(THRYW(JY,JX));       THRYW=0._r8
  allocate(THRMW(JY,JX));       THRMW=0._r8
  allocate(THRMS(JY,JX));       THRMS=0._r8
  allocate(THRMR(JY,JX));       THRMR=0._r8
  allocate(THRYG(JY,JX));       THRYG=0._r8
  allocate(THRYR(JY,JX));       THRYR=0._r8
  allocate(RADXW(JY,JX));       RADXW=0._r8
  allocate(RADXG(JY,JX));       RADXG=0._r8
  allocate(RADXR(JY,JX));       RADXR=0._r8
  allocate(FLWLX(3,JD,JV,JH));  FLWLX=0._r8
  allocate(TFLWLX(JZ,JY,JX));   TFLWLX=0._r8
  allocate(FLU1(JZ,JY,JX));     FLU1=0._r8
  allocate(HWFLU1(JZ,JY,JX));   HWFLU1=0._r8
  allocate(PSISM1(0:JZ,JY,JX)); PSISM1=0._r8
  allocate(ALTG(JY,JX));        ALTG=0._r8
  allocate(WFLXLH(JZ,JY,JX));   WFLXLH=0._r8
  allocate(DLYRR(JY,JX));       DLYRR=0._r8
  allocate(WFLXR(JY,JX));       WFLXR=0._r8
  allocate(TFLXR(JY,JX));       TFLXR=0._r8
  allocate(HCNDR(JY,JX));       HCNDR=0._r8
  allocate(CNDH1(JZ,JY,JX));    CNDH1=0._r8
  allocate(VOLA1(0:JZ,JY,JX));  VOLA1=0._r8
  allocate(THETWX(0:JZ,JY,JX)); THETWX=0._r8
  allocate(THETIX(0:JZ,JY,JX)); THETIX=0._r8
  allocate(THETPX(0:JZ,JY,JX)); THETPX=0._r8
  allocate(VOLAH1(JZ,JY,JX));   VOLAH1=0._r8
  allocate(VOLWH1(JZ,JY,JX));   VOLWH1=0._r8
  allocate(VOLPH1(JZ,JY,JX));   VOLPH1=0._r8
  allocate(VOLIH1(JZ,JY,JX));   VOLIH1=0._r8
  allocate(THETPY(0:JZ,JY,JX)); THETPY=0._r8
  allocate(FLWNX(JY,JX));       FLWNX=0._r8
  allocate(FLWXNX(JY,JX));      FLWXNX=0._r8
  allocate(FLWHNX(JY,JX));      FLWHNX=0._r8
  allocate(HFLWNX(JY,JX));      HFLWNX=0._r8
  allocate(PSISA1(JZ,JY,JX));   PSISA1=0._r8
  allocate(TQS1(JY,JX));        TQS1=0._r8
  allocate(TQW1(JY,JX));        TQW1=0._r8
  allocate(TQI1(JY,JX));        TQI1=0._r8
  allocate(THQS1(JY,JX));       THQS1=0._r8
  allocate(TFLWS(JS,JY,JX));    TFLWS=0._r8
  allocate(TFLWW(JS,JY,JX));    TFLWW=0._r8
  allocate(TFLWI(JS,JY,JX));    TFLWI=0._r8
  allocate(THFLWW(JS,JY,JX));   THFLWW=0._r8
  allocate(TFLX0(JS,JY,JX));    TFLX0=0._r8
  allocate(WFLXS(JS,JY,JX));    WFLXS=0._r8
  allocate(WFLXI(JS,JY,JX));    WFLXI=0._r8
  allocate(FLW0W(JS,JY,JX));    FLW0W=0._r8
  allocate(FLW0S(JS,JY,JX));    FLW0S=0._r8
  allocate(FLW0I(JS,JY,JX));    FLW0I=0._r8
  allocate(HFLW0W(JS,JY,JX));   HFLW0W=0._r8
  allocate(VOLS0M(JS,JY,JX));   VOLS0M=0._r8
  allocate(VOLW0M(JS,JY,JX));   VOLW0M=0._r8
  allocate(VOLI0M(JS,JY,JX));   VOLI0M=0._r8
  allocate(VHCPWMM(JS,JY,JX));  VHCPWMM=0._r8
  allocate(TK0M(JS,JY,JX));     TK0M=0._r8
  end subroutine InitWatSubData

!----------------------------------------------------------------------
  subroutine DestructWatSubData
  use abortutils, only : destroy
  implicit none

  call destroy(N6X)
  call destroy(VOLWX1)
  call destroy(VPQ)
  call destroy(TKQ)
  call destroy(XVOLT)
  call destroy(XVOLW)
  call destroy(XVOLI)
  call destroy(FMAC)
  call destroy(FGRD)
  call destroy(VOLW1)
  call destroy(VOLI1)
  call destroy(VHCP1)
  call destroy(VHCP1A)
  call destroy(VHCP1B)
  call destroy(TK1)
  call destroy(TWFLXL)
  call destroy(VOLW2)
  call destroy(VOLP1)
  call destroy(TWFLXH)
  call destroy(PRECM)
  call destroy(VOLS0)
  call destroy(VOLI0)
  call destroy(VOLW0)
  call destroy(VOLS1)
  call destroy(DLYRS0)
  call destroy(VOLP1Z)
  call destroy(TK0)
  call destroy(AREAU)
  call destroy(AREAUD)
  call destroy(FLQ0S)
  call destroy(FLQ0I)
  call destroy(FLQ0W)
  call destroy(FLQ1)
  call destroy(FLH1)
  call destroy(FLY1)
  call destroy(HWFLQ0)
  call destroy(HWFLQ1)
  call destroy(HWFLY1)
  call destroy(RAR)
  call destroy(RAGS)
  call destroy(BAREW)
  call destroy(CVRDW)
  call destroy(RARG)
  call destroy(RAGR)
  call destroy(RAGW)
  call destroy(PAREW)
  call destroy(PAREG)
  call destroy(RAG)
  call destroy(PARSW)
  call destroy(PARSG)
  call destroy(PARER)
  call destroy(PARSR)
  call destroy(QR1)
  call destroy(HQR1)
  call destroy(VOLPH1Z)
  call destroy(QS1)
  call destroy(QW1)
  call destroy(QI1)
  call destroy(HQS1)
  call destroy(TQR1)
  call destroy(THQR1)
  call destroy(EVAPG)
  call destroy(TTFLXL)
  call destroy(EVAPW)
  call destroy(EVAPS)
  call destroy(EVAPR)
  call destroy(FLWRL)
  call destroy(HFLWRL)
  call destroy(FINHL)
  call destroy(FLWL)
  call destroy(FLWHL)
  call destroy(HFLWL)
  call destroy(TFLWL)
  call destroy(TFLWHL)
  call destroy(THFLWL)
  call destroy(WFLXL)
  call destroy(TFLXL)
  call destroy(AVCNHL)
  call destroy(THRYW)
  call destroy(THRMW)
  call destroy(THRMS)
  call destroy(THRMR)
  call destroy(THRYG)
  call destroy(THRYR)
  call destroy(RADXW)
  call destroy(RADXG)
  call destroy(RADXR)
  call destroy(FLWLX)
  call destroy(TFLWLX)
  call destroy(FLU1)
  call destroy(HWFLU1)
  call destroy(PSISM1)
  call destroy(ALTG)
  call destroy(WFLXLH)
  call destroy(DLYRR)
  call destroy(WFLXR)
  call destroy(TFLXR)
  call destroy(HCNDR)
  call destroy(CNDH1)
  call destroy(VOLA1)
  call destroy(THETWX)
  call destroy(THETIX)
  call destroy(THETPX)
  call destroy(VOLAH1)
  call destroy(VOLWH1)
  call destroy(VOLPH1)
  call destroy(VOLIH1)
  call destroy(THETPY)
  call destroy(FLWNX)
  call destroy(FLWXNX)
  call destroy(FLWHNX)
  call destroy(HFLWNX)
  call destroy(PSISA1)
  call destroy(TQS1)
  call destroy(TQW1)
  call destroy(TQI1)
  call destroy(THQS1)
  call destroy(TFLWS)
  call destroy(TFLWW)
  call destroy(TFLWI)
  call destroy(THFLWW)
  call destroy(TFLX0)
  call destroy(WFLXS)
  call destroy(WFLXI)
  call destroy(FLW0W)
  call destroy(FLW0S)
  call destroy(FLW0I)
  call destroy(HFLW0W)
  call destroy(VOLS0M)
  call destroy(VOLW0M)
  call destroy(VOLI0M)
  call destroy(VHCPWMM)
  call destroy(TK0M)
  end subroutine DestructWatSubData

end module WatsubDataMod
