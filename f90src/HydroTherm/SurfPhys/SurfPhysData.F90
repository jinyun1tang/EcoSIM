module SurfPhysData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  save
  character(len=*), private, parameter :: mod_filename=__FILE__

  real(r8),allocatable ::  XVOLT(:,:)                         !
  real(r8),allocatable ::  XVWatMicP(:,:)                         !
  real(r8),allocatable ::  XViceMicP(:,:)                         !
  real(r8),allocatable ::  THRMR(:,:)                         !
  real(r8),allocatable ::  VOLA10(:,:)                        !
  real(r8),allocatable ::  THRYR(:,:)                         !
  real(r8),allocatable ::  THRMS(:,:)                         !
  real(r8),allocatable ::  RADXR(:,:)                         !
  real(r8),allocatable ::  THRYG(:,:)                         !
  real(r8),allocatable ::  RADXG(:,:)                         !
  real(r8),allocatable ::  RAGR(:,:)                          !
  real(r8),allocatable ::  PAREG(:,:)                         !
  real(r8),allocatable ::  PARER(:,:)                         !
  real(r8),allocatable ::  RARG(:,:)                          !
  real(r8),allocatable ::  PARSR(:,:)                         !
  real(r8),allocatable ::  RAR(:,:)                           !
  real(r8),allocatable ::  FLWRL(:,:)                         !
  real(r8),allocatable ::  HFLWRL(:,:)                        !  
  real(r8),allocatable ::  PARSG(:,:)                         !
  real(r8),allocatable ::  EVAPR(:,:)                         !
  real(r8),allocatable ::  EVAPG(:,:)                         !  
  real(r8),allocatable ::  WFLXR(:,:)                         !
  real(r8),allocatable ::  RAGS(:,:)                          !    
  real(r8),allocatable ::  TFLXR(:,:)                         !  
  real(r8),allocatable ::  CVRDW(:,:)                         !
  real(r8),allocatable ::  FLH1(:,:)                          !
  real(r8),allocatable ::  PRECM(:,:)                         !
  real(r8),allocatable ::  FLQ1(:,:)                          !
  real(r8),allocatable ::  HWFLQ1(:,:)                        !
  real(r8),allocatable ::  HWFLY1(:,:)                        !
  real(r8),allocatable ::  FLY1(:,:)                          !
  real(r8),allocatable ::  BAREW(:,:)                         !
  real(r8),allocatable ::  HCNDR(:,:)                         !

  public :: InitSurfPhysData,DestructSurfPhysData

  contains
!------------------------------------------------------------------------------------------  

  subroutine InitSurfPhysData  
  implicit none

  allocate(XVOLT(JY,JX));       XVOLT=0._r8
  allocate(XVWatMicP(JY,JX));       XVWatMicP=0._r8
  allocate(XViceMicP(JY,JX));       XViceMicP=0._r8
  allocate(THRMR(JY,JX));       THRMR=0._r8
  allocate(VOLA10(JY,JX));      VOLA10=0._r8
  allocate(THRYR(JY,JX));       THRYR=0._r8
  allocate(THRMS(JY,JX));       THRMS=0._r8
  allocate(RADXR(JY,JX));       RADXR=0._r8    
  allocate(THRYG(JY,JX));       THRYG=0._r8
  allocate(RADXG(JY,JX));       RADXG=0._r8
  allocate(RAGR(JY,JX));        RAGR=0._r8
  allocate(PAREG(JY,JX));       PAREG=0._r8
  allocate(PARER(JY,JX));       PARER=0._r8
  allocate(RARG(JY,JX));        RARG=0._r8
  allocate(PARSR(JY,JX));       PARSR=0._r8
  allocate(RAR(JY,JX));         RAR=0._r8  
  allocate(FLWRL(JY,JX));       FLWRL=0._r8
  allocate(HFLWRL(JY,JX));      HFLWRL=0._r8  
  allocate(PARSG(JY,JX));       PARSG=0._r8  
  allocate(EVAPR(JY,JX));       EVAPR=0._r8  
  allocate(EVAPG(JY,JX));       EVAPG=0._r8  
  allocate(WFLXR(JY,JX));       WFLXR=0._r8  
  allocate(RAGS(JY,JX));        RAGS=0._r8  
  allocate(TFLXR(JY,JX));       TFLXR=0._r8  
  allocate(CVRDW(JY,JX));       CVRDW=0._r8
  allocate(FLH1(JY,JX));        FLH1=0._r8
  allocate(PRECM(JY,JX));       PRECM=0._r8
  allocate(FLQ1(JY,JX));        FLQ1=0._r8
  allocate(HWFLQ1(JY,JX));      HWFLQ1=0._r8
  allocate(HWFLY1(JY,JX));      HWFLY1=0._r8  
  allocate(FLY1(JY,JX));        FLY1=0._r8  
  allocate(BAREW(JY,JX));       BAREW=0._r8  
  allocate(HCNDR(JY,JX));       HCNDR=0._r8  
  end subroutine InitSurfPhysData  
!------------------------------------------------------------------------------------------  

  subroutine DestructSurfPhysData
  use abortutils, only : destroy
  implicit none

  call destroy(XVOLT)
  call destroy(XVWatMicP)
  call destroy(XViceMicP)
  call destroy(THRMR)
  call destroy(VOLA10)
  call destroy(THRYR)
  call destroy(THRMS)
  call destroy(RADXR)  
  call destroy(THRYG)
  call destroy(RADXG)
  call destroy(RAGR)
  call destroy(PAREG)
  call destroy(PARER)
  call destroy(RARG)
  call destroy(PARSR)
  call destroy(RAR)  
  call destroy(FLWRL)
  call destroy(HFLWRL)  
  call destroy(PARSG)  
  call destroy(EVAPR)
  call destroy(EVAPG)
  call destroy(WFLXR)  
  call destroy(RAGS) 
  call destroy(TFLXR)
  call destroy(CVRDW)
  call destroy(FLH1)
  call destroy(PRECM)
  call destroy(FLQ1)
  call destroy(HWFLQ1)
  call destroy(HWFLY1)
  call destroy(FLY1)
  call destroy(BAREW)
  call destroy(HCNDR)  
  end subroutine DestructSurfPhysData

end module SurfPhysData
