module WatsubDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  VOLWX1(:,:,:)                      !

  real(r8),allocatable ::  XVOLT(:,:)                         !
  real(r8),allocatable ::  XVOLW(:,:)                         !
  real(r8),allocatable ::  XVOLI(:,:)                         !

  real(r8),allocatable ::  TWFLXL(:,:,:)                      !
  real(r8),allocatable ::  VOLW2(:,:,:)                       !

  real(r8),allocatable ::  TWFLXH(:,:,:)                      !
  real(r8),allocatable ::  PRECM(:,:)                         !
  real(r8),allocatable ::  VOLP1Z(:,:,:)                      !

  real(r8),allocatable ::  AREAU(:,:,:)                       !
  real(r8),allocatable ::  AREAUD(:,:,:)                      !


  real(r8),allocatable ::  FLQ1(:,:)                          !
  real(r8),allocatable ::  FLH1(:,:)                          !
  real(r8),allocatable ::  FLY1(:,:)                          !

  real(r8),allocatable ::  HWFLQ1(:,:)                        !
  real(r8),allocatable ::  HWFLY1(:,:)                        !
  real(r8),allocatable ::  RAR(:,:)                           !
  real(r8),allocatable ::  RAGS(:,:)                          !
  real(r8),allocatable ::  BAREW(:,:)                         !
  real(r8),allocatable ::  CVRDW(:,:)                         !
  real(r8),allocatable ::  RARG(:,:)                          !
  real(r8),allocatable ::  RAGR(:,:)                          !

  real(r8),allocatable ::  PAREG(:,:)                         !

  real(r8),allocatable ::  PARSG(:,:)                         !
  real(r8),allocatable ::  PARER(:,:)                         !
  real(r8),allocatable ::  PARSR(:,:)                         !

  real(r8),allocatable ::  VOLPH1Z(:,:,:)                     !
  real(r8),allocatable ::  EVAPG(:,:)                         !
  real(r8),allocatable ::  TTFLXL(:,:,:)                      !


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

  real(r8),allocatable ::  THRMS(:,:)                         !
  real(r8),allocatable ::  THRMR(:,:)                         !
  real(r8),allocatable ::  THRYG(:,:)                         !
  real(r8),allocatable ::  THRYR(:,:)                         !

  real(r8),allocatable ::  RADXG(:,:)                         !
  real(r8),allocatable ::  RADXR(:,:)                         !
  real(r8),allocatable ::  FLWLX(:,:,:,:)                     !
  real(r8),allocatable ::  TFLWLX(:,:,:)                      !
  real(r8),allocatable ::  FLU1(:,:,:)                        !
  real(r8),allocatable ::  HWFLU1(:,:,:)                      !

  real(r8),allocatable ::  WFLXLH(:,:,:)                      !

  real(r8),allocatable ::  WFLXR(:,:)                         !
  real(r8),allocatable ::  TFLXR(:,:)                         !
  real(r8),allocatable ::  HCNDR(:,:)                         !
  real(r8),allocatable ::  CNDH1(:,:,:)                       !
  real(r8),allocatable ::  VOLA1(:,:,:)                       !
  real(r8),allocatable ::  VOLAH1(:,:,:)                      !

  real(r8),allocatable ::  FLWNX(:,:)                         !
  real(r8),allocatable ::  FLWXNX(:,:)                        !
  real(r8),allocatable ::  FLWHNX(:,:)                        !
  real(r8),allocatable ::  HFLWNX(:,:)                        !
  real(r8),allocatable ::  PSISA1(:,:,:)                      !


  integer, allocatable ::  N6X(:,:)

!----------------------------------------------------------------------
  public :: InitWatsubData
contains
  subroutine InitWatSubData

  implicit none

  allocate(N6X(JY,JX));         N6X=0
  allocate(VOLWX1(JZ,JY,JX));   VOLWX1=0._r8

  allocate(XVOLT(JY,JX));       XVOLT=0._r8
  allocate(XVOLW(JY,JX));       XVOLW=0._r8
  allocate(XVOLI(JY,JX));       XVOLI=0._r8


  allocate(TWFLXL(JZ,JY,JX));   TWFLXL=0._r8
  allocate(VOLW2(JZ,JY,JX));    VOLW2=0._r8

  allocate(TWFLXH(JZ,JY,JX));   TWFLXH=0._r8
  allocate(PRECM(JY,JX));       PRECM=0._r8
  allocate(VOLP1Z(JZ,JY,JX));   VOLP1Z=0._r8

  allocate(AREAU(JZ,JY,JX));    AREAU=0._r8
  allocate(AREAUD(JZ,JY,JX));   AREAUD=0._r8

  allocate(FLQ1(JY,JX));        FLQ1=0._r8
  allocate(FLH1(JY,JX));        FLH1=0._r8
  allocate(FLY1(JY,JX));        FLY1=0._r8

  allocate(HWFLQ1(JY,JX));      HWFLQ1=0._r8
  allocate(HWFLY1(JY,JX));      HWFLY1=0._r8
  allocate(RAR(JY,JX));         RAR=0._r8
  allocate(RAGS(JY,JX));        RAGS=0._r8
  allocate(BAREW(JY,JX));       BAREW=0._r8
  allocate(CVRDW(JY,JX));       CVRDW=0._r8
  allocate(RARG(JY,JX));        RARG=0._r8
  allocate(RAGR(JY,JX));        RAGR=0._r8

  allocate(PAREG(JY,JX));       PAREG=0._r8

  allocate(PARSG(JY,JX));       PARSG=0._r8
  allocate(PARER(JY,JX));       PARER=0._r8
  allocate(PARSR(JY,JX));       PARSR=0._r8
  allocate(VOLPH1Z(JZ,JY,JX));  VOLPH1Z=0._r8
  allocate(EVAPG(JY,JX));       EVAPG=0._r8
  allocate(TTFLXL(JZ,JY,JX));   TTFLXL=0._r8

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

  allocate(THRMS(JY,JX));       THRMS=0._r8
  allocate(THRMR(JY,JX));       THRMR=0._r8
  allocate(THRYG(JY,JX));       THRYG=0._r8
  allocate(THRYR(JY,JX));       THRYR=0._r8

  allocate(RADXG(JY,JX));       RADXG=0._r8
  allocate(RADXR(JY,JX));       RADXR=0._r8
  allocate(FLWLX(3,JD,JV,JH));  FLWLX=0._r8
  allocate(TFLWLX(JZ,JY,JX));   TFLWLX=0._r8
  allocate(FLU1(JZ,JY,JX));     FLU1=0._r8
  allocate(HWFLU1(JZ,JY,JX));   HWFLU1=0._r8

  allocate(WFLXLH(JZ,JY,JX));   WFLXLH=0._r8

  allocate(WFLXR(JY,JX));       WFLXR=0._r8
  allocate(TFLXR(JY,JX));       TFLXR=0._r8
  allocate(HCNDR(JY,JX));       HCNDR=0._r8
  allocate(CNDH1(JZ,JY,JX));    CNDH1=0._r8
  allocate(VOLA1(0:JZ,JY,JX));  VOLA1=0._r8

  allocate(VOLAH1(JZ,JY,JX));   VOLAH1=0._r8


  allocate(FLWNX(JY,JX));       FLWNX=0._r8
  allocate(FLWXNX(JY,JX));      FLWXNX=0._r8
  allocate(FLWHNX(JY,JX));      FLWHNX=0._r8
  allocate(HFLWNX(JY,JX));      HFLWNX=0._r8
  allocate(PSISA1(JZ,JY,JX));   PSISA1=0._r8

  end subroutine InitWatSubData

!----------------------------------------------------------------------
  subroutine DestructWatSubData
  use abortutils, only : destroy
  implicit none

  call destroy(N6X)
  call destroy(VOLWX1)
  call destroy(XVOLT)
  call destroy(XVOLW)
  call destroy(XVOLI)

  call destroy(TWFLXL)
  call destroy(VOLW2)

  call destroy(TWFLXH)
  call destroy(PRECM)
  call destroy(VOLP1Z)

  call destroy(AREAU)
  call destroy(AREAUD)

  call destroy(FLQ1)
  call destroy(FLH1)
  call destroy(FLY1)

  call destroy(HWFLQ1)
  call destroy(HWFLY1)
  call destroy(RAR)
  call destroy(RAGS)
  call destroy(BAREW)
  call destroy(CVRDW)
  call destroy(RARG)
  call destroy(RAGR)
  call destroy(PAREG)

  call destroy(PARSG)
  call destroy(PARER)
  call destroy(PARSR)
  call destroy(VOLPH1Z)
  call destroy(EVAPG)
  call destroy(TTFLXL)

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

  call destroy(THRMS)
  call destroy(THRMR)
  call destroy(THRYG)
  call destroy(THRYR)

  call destroy(RADXG)
  call destroy(RADXR)
  call destroy(FLWLX)
  call destroy(TFLWLX)
  call destroy(FLU1)
  call destroy(HWFLU1)

  call destroy(WFLXLH)

  call destroy(WFLXR)
  call destroy(TFLXR)
  call destroy(HCNDR)
  call destroy(CNDH1)
  call destroy(VOLA1)

  call destroy(VOLAH1)

  call destroy(FLWNX)
  call destroy(FLWXNX)
  call destroy(FLWHNX)
  call destroy(HFLWNX)
  call destroy(PSISA1)

  end subroutine DestructWatSubData

end module WatsubDataMod
