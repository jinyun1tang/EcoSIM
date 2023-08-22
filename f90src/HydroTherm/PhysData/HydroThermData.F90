module HydroThermData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__
  real(r8),allocatable ::  THETWX(:,:,:)                      !
  real(r8),allocatable ::  PSISM1(:,:,:)                      !  
  real(r8),allocatable ::  TK1(:,:,:)                         !  
  real(r8),allocatable ::  DLYRR(:,:)                         !  
  real(r8),allocatable ::  THETIX(:,:,:)                      !  
  real(r8),allocatable ::  THETPX(:,:,:)                      !  
  real(r8),allocatable ::  VolHeatCapacity(:,:,:)                       !
  real(r8),allocatable ::  THETPY(:,:,:)                      !
  real(r8),allocatable ::  VOLW1(:,:,:)                       !    
  real(r8),allocatable ::  VOLI1(:,:,:)                       !  
  real(r8),allocatable ::  FLQ0S(:,:)                         !  
  real(r8),allocatable ::  FLQ0W(:,:)                         !  
  real(r8),allocatable ::  VOLPH1(:,:,:)                      !  
  real(r8),allocatable ::  EVAPW(:,:)                         !  
  real(r8),allocatable ::  EVAPS(:,:)                         !  
  real(r8),allocatable ::  EVAPSN(:,:)
  real(r8),allocatable ::  FMAC(:,:,:)                        !  
  real(r8),allocatable ::  HWFLQ0(:,:)                        !  
  real(r8),allocatable ::  VPQ(:,:)                           !  
  real(r8),allocatable ::  PARSW(:,:)                         !  
  real(r8),allocatable ::  FLQ0I(:,:)                         !  
  real(r8),allocatable ::  TKQ(:,:)                           !  
  real(r8),allocatable ::  THRMW(:,:)                         !  
  real(r8),allocatable ::  RAGW(:,:)                          !  
  real(r8),allocatable ::  THRYW(:,:)                         !  
  real(r8),allocatable ::  VOLP1(:,:,:)                       !  
  real(r8),allocatable ::  TK0(:,:,:)                         !  
  real(r8),allocatable ::  VOLWH1(:,:,:)                      !
  real(r8),allocatable ::  FGRD(:,:,:)                        !
  real(r8),allocatable ::  VolHeatCapacityA(:,:,:)                      !
  real(r8),allocatable ::  VolHeatCapacityB(:,:,:)                      !
  real(r8),allocatable ::  RAG(:,:)                           !
  real(r8),allocatable ::  PAREW(:,:)                         ! 
  real(r8),allocatable ::  ALTG(:,:)                          ! 
  real(r8),allocatable ::  VOLIH1(:,:,:)                      !
  real(r8),allocatable ::  RADXW(:,:)                         !  
  real(r8),allocatable ::  QR1(:,:,:,:)                       !  
  real(r8),allocatable ::  HQR1(:,:,:,:)                      !  
  real(r8),allocatable ::  QS1(:,:,:)                         !  
  real(r8),allocatable ::  QW1(:,:,:)                         !
  real(r8),allocatable ::  QI1(:,:,:)                         !  
  real(r8),allocatable ::  HQS1(:,:,:)                        !  
  real(r8),allocatable ::  HFLWL(:,:,:,:)                     !  
  real(r8),allocatable ::  FLWL(:,:,:,:)                      !  
  real(r8),allocatable ::  FLWHL(:,:,:,:)                     !  
  real(r8),allocatable ::  FLWLX(:,:,:,:)                     !  
  real(r8),allocatable ::  VOLW2(:,:,:)                       !
  real(r8),allocatable ::  VOLP1Z(:,:,:)                      !
  real(r8),allocatable ::  VOLWX1(:,:,:)                      !

  public :: InitHydroThermData
  public :: DestructHydroThermData
  contains

!------------------------------------------------------------------------------------------
  subroutine InitHydroThermData
  implicit none


  allocate(THETWX(0:JZ,JY,JX)); THETWX=0._r8
  allocate(PSISM1(0:JZ,JY,JX)); PSISM1=0._r8
  allocate(TK1(0:JZ,JY,JX));    TK1=0._r8    
  allocate(DLYRR(JY,JX));       DLYRR=0._r8  
  allocate(THETIX(0:JZ,JY,JX)); THETIX=0._r8 
  allocate(THETPX(0:JZ,JY,JX)); THETPX=0._r8   
  allocate(VolHeatCapacity(0:JZ,JY,JX));  VolHeatCapacity=0._r8  
  allocate(THETPY(0:JZ,JY,JX)); THETPY=0._r8  
  allocate(VOLW1(0:JZ,JY,JX));  VOLW1=0._r8  
  allocate(VOLI1(0:JZ,JY,JX));  VOLI1=0._r8
  allocate(FLQ0S(JY,JX));       FLQ0S=0._r8
  allocate(FLQ0W(JY,JX));       FLQ0W=0._r8
  allocate(VOLPH1(JZ,JY,JX));   VOLPH1=0._r8 
  allocate(EVAPSN(JY,JX));      EVAPSN=0._r8  
  allocate(EVAPW(JY,JX));       EVAPW=0._r8    
  allocate(EVAPS(JY,JX));       EVAPS=0._r8  
  allocate(FMAC(JZ,JY,JX));     FMAC=0._r8  
  allocate(HWFLQ0(JY,JX));      HWFLQ0=0._r8  
  allocate(VPQ(JY,JX));         VPQ=0._r8  
  allocate(PARSW(JY,JX));       PARSW=0._r8
  allocate(FLQ0I(JY,JX));       FLQ0I=0._r8    
  allocate(TKQ(JY,JX));         TKQ=0._r8  
  allocate(THRMW(JY,JX));       THRMW=0._r8  
  allocate(RAGW(JY,JX));        RAGW=0._r8  
  allocate(THRYW(JY,JX));       THRYW=0._r8  
  allocate(VOLP1(0:JZ,JY,JX));  VOLP1=0._r8  
  allocate(TK0(JS,JY,JX));      TK0=0._r8  
  allocate(VOLWH1(JZ,JY,JX));   VOLWH1=0._r8
  allocate(FGRD(JZ,JY,JX));     FGRD=0._r8
  allocate(VolHeatCapacityA(JZ,JY,JX));   VolHeatCapacityA=0._r8
  allocate(VolHeatCapacityB(JZ,JY,JX));   VolHeatCapacityB=0._r8  
  allocate(RAG(JY,JX));         RAG=0._r8
  allocate(PAREW(JY,JX));       PAREW=0._r8
  allocate(ALTG(JY,JX));        ALTG=0._r8  
  allocate(VOLIH1(JZ,JY,JX));   VOLIH1=0._r8  
  allocate(RADXW(JY,JX));       RADXW=0._r8  
  allocate(QR1(2,2,JV,JH));     QR1=0._r8  
  allocate(HQR1(2,2,JV,JH));    HQR1=0._r8
  allocate(QS1(2,JV,JH));       QS1=0._r8
  allocate(QW1(2,JV,JH));       QW1=0._r8
  allocate(QI1(2,JV,JH));       QI1=0._r8  
  allocate(HQS1(2,JV,JH));      HQS1=0._r8
  allocate(HFLWL(3,JD,JV,JH));  HFLWL=0._r8
  allocate(FLWL(3,JD,JV,JH));   FLWL=0._r8  
  allocate(FLWHL(3,JD,JV,JH));  FLWHL=0._r8
  allocate(FLWLX(3,JD,JV,JH));  FLWLX=0._r8
  allocate(VOLW2(JZ,JY,JX));    VOLW2=0._r8
  allocate(VOLP1Z(JZ,JY,JX));   VOLP1Z=0._r8
  allocate(VOLWX1(JZ,JY,JX));   VOLWX1=0._r8

  end subroutine InitHydroThermData

!------------------------------------------------------------------------------------------
  subroutine DestructHydroThermData
  use abortutils, only : destroy  
  implicit none

  call destroy(QR1)
  call destroy(THETWX)
  call destroy(PSISM1)  
  call destroy(TK1)  
  call destroy(DLYRR)  
  call destroy(THETIX)  
  call destroy(THETPX)  
  call destroy(VolHeatCapacity)  
  call destroy(THETPY)  
  call destroy(VOLW1)  
  call destroy(VOLI1)
  call destroy(FLQ0S)
  call destroy(FLQ0W)
  call destroy(VOLPH1)    
  call destroy(EVAPW)  
  call destroy(EVAPS)  
  call destroy(EVAPSN)  
  call destroy(FMAC)  
  call destroy(HWFLQ0)
  call destroy(VPQ)    
  call destroy(PARSW)  
  call destroy(FLQ0I)
  call destroy(TKQ)
  call destroy(THRMW)
  call destroy(RAGW) 
  call destroy(THRYW)
  call destroy(VOLP1)
  call destroy(TK0)  
  call destroy(VOLWH1)
  call destroy(FGRD)
  call destroy(VolHeatCapacityA)
  call destroy(VolHeatCapacityB)  
  call destroy(RAG)
  call destroy(PAREW)
  call destroy(ALTG)
  call destroy(VOLIH1)
  call destroy(RADXW)
  call destroy(HQR1)    
  call destroy(QS1)  
  call destroy(QW1)
  call destroy(QI1)
  call destroy(HQS1)    
  call destroy(HFLWL)  
  call destroy(FLWL)  
  call destroy(FLWHL)
  call destroy(FLWLX)
  call destroy(VOLW2)
  call destroy(VOLP1Z)
  call destroy(VOLWX1)

  end subroutine DestructHydroThermData
end module HydroThermData