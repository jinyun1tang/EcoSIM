module SnowPhysData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod

  implicit none

  character(len=*), private, parameter :: mod_filename=__FILE__

  real(r8),allocatable ::  trcg_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcn_TBLS(:,:,:,:)
  real(r8),allocatable ::  trcsa_TBLS(:,:,:,:)                      !
  real(r8),allocatable ::  TQS(:,:)                           !
  real(r8),allocatable :: trcsa_TQS(:,:,:)  
  real(r8),allocatable ::  trcn_TQR(:,:,:)                        !
  real(r8),allocatable ::  trcsa_TQR(:,:,:)                         !
  real(r8),allocatable ::  trcg_QSS(:,:,:)
  real(r8),allocatable ::  trcn_QSS(:,:,:)
  real(r8),allocatable ::  trcg_TQR(:,:,:)                        !

  real(r8),allocatable ::  THQS(:,:)                          !
  real(r8),allocatable ::  TQI(:,:)                           !
  real(r8),allocatable ::  TQW(:,:)                           !
  real(r8),allocatable ::  TQR1(:,:)                          !
  real(r8),allocatable ::  THQR1(:,:)                         !
  real(r8),allocatable ::  TQS1(:,:)                          !
  real(r8),allocatable ::  TQW1(:,:)                          !
  real(r8),allocatable ::  TQI1(:,:)                          !
  real(r8),allocatable ::  THQS1(:,:)                         !
  real(r8),allocatable ::  TK0M(:,:,:)                        !
  real(r8),allocatable ::  VHCPWMM(:,:,:)                     !  
  real(r8),allocatable ::  VOLI0M(:,:,:)                      !
  real(r8),allocatable ::  VOLW0M(:,:,:)                      !
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
  real(r8),allocatable ::  VOLS0(:,:,:)                       !
  real(r8),allocatable ::  VOLI0(:,:,:)                       !
  real(r8),allocatable ::  VOLW0(:,:,:)                       !
  real(r8),allocatable ::  VOLS1(:,:,:)                       !
  real(r8),allocatable ::  DLYRS0(:,:,:)                      !
   
  public :: InitSnowPhysData
  public :: DestructSnowPhysData
  contains
!----------------------------------------------------------------------  
  subroutine  InitSnowPhysData
  implicit none

  allocate(trcg_TBLS(idg_beg:idg_end-1,JS,JY,JX));        trcg_TBLS=0._r8
  allocate(trcn_TBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_TBLS=0._r8
  allocate(trcsa_TBLS(idsa_beg:idsa_end,JS,JY,JX));       trcsa_TBLS=0._r8  
  allocate(trcsa_TQS(idsa_beg:idsa_end,JY,JX));           trcsa_TQS=0._r8  
  allocate(trcn_TQR(ids_nut_beg:ids_nuts_end,JY,JX));     trcn_TQR=0._r8  
  allocate(trcsa_TQR(idsa_beg:idsa_end,JY,JX));           trcsa_TQR=0._r8  
  allocate(trcg_QSS(idg_beg:idg_end-1,JY,JX));            trcg_QSS=0._r8
  allocate(trcn_QSS(ids_nut_beg:ids_nuts_end,JY,JX));trcn_QSS=0._r8
  allocate(trcg_TQR(idg_beg:idg_end-1,JY,JX));      trcg_TQR=0._r8

  allocate(TQI(JY,JX));         TQI=0._r8
  allocate(TQW(JY,JX));         TQW=0._r8
  allocate(TQS(JY,JX));         TQS=0._r8
  allocate(THQS(JY,JX));        THQS=0._r8
  allocate(TK0M(JS,JY,JX));     TK0M=0._r8
  allocate(THQR1(JY,JX));       THQR1=0._r8
  allocate(THQS1(JY,JX));       THQS1=0._r8
  allocate(TQR1(JY,JX));        TQR1=0._r8
  allocate(TQS1(JY,JX));        TQS1=0._r8
  allocate(TQW1(JY,JX));        TQW1=0._r8
  allocate(TQI1(JY,JX));        TQI1=0._r8
  allocate(VHCPWMM(JS,JY,JX));  VHCPWMM=0._r8
  allocate(VOLI0M(JS,JY,JX));   VOLI0M=0._r8
  allocate(VOLW0M(JS,JY,JX));   VOLW0M=0._r8
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
  allocate(VOLS0(JS,JY,JX));    VOLS0=0._r8
  allocate(VOLI0(JS,JY,JX));    VOLI0=0._r8
  allocate(VOLW0(JS,JY,JX));    VOLW0=0._r8
  allocate(VOLS1(JS,JY,JX));    VOLS1=0._r8
  allocate(DLYRS0(JS,JY,JX));   DLYRS0=0._r8

  end subroutine  InitSnowPhysData
!----------------------------------------------------------------------  
  subroutine DestructSnowPhysData
  use abortutils, only : destroy
  implicit none

  call destroy(trcg_TBLS)
  call destroy(trcn_TBLS)
  call destroy(trcsa_TBLS)
  call destroy(trcn_TQR)  
  call destroy(trcsa_TQR)
  call destroy(trcg_QSS)
  call destroy(trcn_QSS)  
  call destroy(trcg_TQR)
  call destroy(TQW)  
  call destroy(TQR1)
  call destroy(TQS1)
  call destroy(TQW1)
  call destroy(TQI1)
  call destroy(TQI)  
  call destroy(THQS)  
  call destroy(THQS1)
  call destroy(TK0M)  
  call destroy(VHCPWMM)
  call destroy(VOLI0M)
  call destroy(VOLW0M)
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
  call destroy(VOLS0)
  call destroy(VOLI0)
  call destroy(VOLW0)
  call destroy(VOLS1)
  call destroy(DLYRS0)
  call destroy(trcsa_TQS)
  call destroy(TQS)
  end subroutine DestructSnowPhysData
end module SnowPhysData